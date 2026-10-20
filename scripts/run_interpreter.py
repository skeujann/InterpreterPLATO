"""CLI for one audited PLATO-only mechanism hypothesis or a JSON batch."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
import argparse
import json
import math
import sys
from typing import Any, Dict, List, Sequence, Tuple

PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE.parent))

from InterpreterPLATO import ExplanationPipeline
from InterpreterPLATO.connectivity import validate_path_search
from InterpreterPLATO.targets import CACHE_SCHEMA_VERSION, PlatoPredictionCache


REQUIRED_BATCH_FIELDS = (
    "drug1_name",
    "drug2_name",
    "drug1_smiles",
    "drug2_smiles",
)
ALLOWED_BATCH_FIELDS = {
    *REQUIRED_BATCH_FIELDS,
    "cell_line",
    "predicted_synergy",
    "advanced",
    "disease_label",
    "model_fidelity_evidence",
}
FORBIDDEN_TARGET_FIELDS = {
    "database",
    "local_targets",
    "drug1_targets",
    "drug2_targets",
}


def add_path_arguments(command: argparse.ArgumentParser) -> None:
    command.add_argument("--path-max-length", type=int, default=5, help="Maximum PPI hops per evidence path (0-20)")
    command.add_argument("--max-evidence-paths", type=int, default=10, help="Maximum displayed target-to-disease paths (0-100)")
    command.add_argument("--max-cross-drug-paths", type=int, default=5, help="Maximum displayed cross-drug PPI paths (0-100)")


def path_options_from_args(args: argparse.Namespace) -> Dict[str, int]:
    return validate_path_search({
        "max_length": args.path_max_length,
        "max_paths": args.max_evidence_paths,
        "max_cross_drug_paths": args.max_cross_drug_paths,
    })


def add_context_arguments(command: argparse.ArgumentParser) -> None:
    add_path_arguments(command)
    context = command.add_mutually_exclusive_group(required=True)
    context.add_argument(
        "--context-config",
        help="Hash-pinned production context JSON/YAML (paths are relative to this file)",
    )
    context.add_argument(
        "--demo",
        action="store_true",
        help="Use bundled demo context; never implied and not suitable for production",
    )
    command.add_argument(
        "--plato-cache",
        default=str(PACKAGE / "data" / "processed" / "plato_cache"),
        help=f"Directory containing audited PLATO cache schema v{CACHE_SCHEMA_VERSION} entries",
    )
    command.add_argument(
        "--plato-release",
        default="r37",
        help="Expected PLATO release label (default: r37)",
    )
    command.add_argument(
        "--trusted-plato-sha256",
        action="append",
        default=[],
        help=(
            "Independent SHA-256 of a PLATO raw export. Repeat for both drugs; "
            "cache self-hashes and user attestation are not enough."
        ),
    )
    command.add_argument(
        "--quick",
        action="store_true",
        help=(
            "Use reduced network randomization budgets for iterative debugging; "
            "the result records an inadequate-budget warning"
        ),
    )


def _load_batch(parser: argparse.ArgumentParser, source_arg: str) -> List[Dict[str, Any]]:
    source = Path(source_arg).expanduser()
    if not source.is_file():
        parser.error(f"batch input is not a readable file: {source}")
    try:
        value = json.loads(
            source.read_text(encoding="utf-8-sig"),
            parse_constant=lambda item: (_ for _ in ()).throw(
                ValueError(f"non-finite JSON constant: {item}")
            ),
        )
    except (OSError, UnicodeError, json.JSONDecodeError, ValueError) as exc:
        parser.error(f"batch input is not strict UTF-8 JSON: {exc}")
    if not isinstance(value, list):
        parser.error("batch input root must be a JSON list")
    if not value:
        parser.error("batch input list must contain at least one drug pair")

    rows: List[Dict[str, Any]] = []
    for index, raw in enumerate(value):
        if not isinstance(raw, Mapping):
            parser.error(f"batch row {index} must be a JSON object")
        row = dict(raw)
        forbidden = sorted(FORBIDDEN_TARGET_FIELDS & set(row))
        if forbidden:
            parser.error(
                f"batch row {index} contains disabled target fields: {', '.join(forbidden)}"
            )
        unknown = sorted(set(row) - ALLOWED_BATCH_FIELDS)
        if unknown:
            parser.error(f"batch row {index} contains unknown fields: {', '.join(unknown)}")
        missing = [
            field
            for field in REQUIRED_BATCH_FIELDS
            if not isinstance(row.get(field), str) or not row[field].strip()
        ]
        if missing:
            parser.error(
                f"batch row {index} requires non-empty strings for: {', '.join(missing)}"
            )
        for field in ("cell_line", "disease_label"):
            if field in row and not isinstance(row[field], str):
                parser.error(f"batch row {index} field {field!r} must be a string")
        if "advanced" in row and not isinstance(row["advanced"], bool):
            parser.error(f"batch row {index} field 'advanced' must be true or false")
        if "predicted_synergy" in row:
            try:
                finite = not isinstance(row["predicted_synergy"], bool) and math.isfinite(
                    float(row["predicted_synergy"])
                )
            except (TypeError, ValueError):
                finite = False
            if not finite:
                parser.error(
                    f"batch row {index} field 'predicted_synergy' must be a finite number"
                )
        if "model_fidelity_evidence" in row and row["model_fidelity_evidence"] is not None:
            if not isinstance(row["model_fidelity_evidence"], Mapping):
                parser.error(
                    f"batch row {index} field 'model_fidelity_evidence' must be an object or null"
                )
        rows.append(row)
    return rows


def _prepare_output(parser: argparse.ArgumentParser, output_arg: str) -> Path:
    destination = Path(output_arg).expanduser()
    if destination.exists() and destination.is_dir():
        parser.error(f"output path is a directory: {destination}")
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        parser.error(f"unable to create output directory {destination.parent}: {exc}")
    return destination


def _preflight_cache(
    parser: argparse.ArgumentParser,
    cache_arg: str,
    queries: Sequence[Tuple[str, str]],
    trusted_export_sha256: Sequence[str],
) -> Path:
    cache_dir = Path(cache_arg).expanduser()
    if not cache_dir.is_dir():
        parser.error(
            f"PLATO cache directory does not exist: {cache_dir}. "
            "Archive each manually downloaded response before running the interpreter."
        )
    cache = PlatoPredictionCache(cache_dir)
    failures: List[str] = []
    trusted = {
        str(value).strip().lower()
        for value in trusted_export_sha256
        if len(str(value).strip()) == 64
        and all(char in "0123456789abcdefABCDEF" for char in str(value).strip())
    }
    seen = set()
    for label, smiles in queries:
        key = (label, smiles)
        if key in seen:
            continue
        seen.add(key)
        try:
            _, metadata = cache.load(smiles)
        except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
            failures.append(f"{label}: unreadable cache entry ({exc})")
            continue
        status = str(metadata.get("status", ""))
        if not status.startswith("cache_hit"):
            failures.append(f"{label}: cache miss")
            continue
        if str(metadata.get("cache_schema_version", "")) != CACHE_SCHEMA_VERSION:
            failures.append(
                f"{label}: cache schema {metadata.get('cache_schema_version')!r}; "
                f"current schema is {CACHE_SCHEMA_VERSION}"
            )
        if metadata.get("schema_fixture_verified") is not True:
            failures.append(
                f"{label}: the response schema is not verified against the controlled r37 fixture"
            )
        source = metadata.get("source", {})
        source_sha256 = ""
        if isinstance(source, Mapping):
            source_sha256 = str(source.get("sha256", "")).strip().lower()
        if not source_sha256 or source_sha256 not in trusted:
            failures.append(
                f"{label}: PLATO raw export SHA-256 is not independently allowlisted"
            )
    if failures:
        parser.error(
            "PLATO cache preflight failed; no score was produced:\n  - "
            + "\n  - ".join(failures)
        )
    return cache_dir


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Generate research mechanism hypotheses from verified PLATO cache entries",
    )
    sub = parser.add_subparsers(dest="command", required=True)
    one = sub.add_parser(
        "summary",
        help="Generate a hypothesis for one named drug pair (primary workflow)",
        description="Analyze one drug pair; --drug1 and --drug2 are the primary inputs.",
    )
    add_context_arguments(one)
    one.add_argument("--drug1", required=True, help="First drug name")
    one.add_argument("--drug2", required=True, help="Second drug name")
    one.add_argument("--drug1-smiles", required=True, help="Exact SMILES submitted to PLATO")
    one.add_argument("--drug2-smiles", required=True, help="Exact SMILES submitted to PLATO")
    one.add_argument(
        "--cell-line",
        default="",
        help="Explicit loaded cell line (required with --context-config)",
    )
    one.add_argument("--disease-label", default="")
    one.add_argument("--synergy", type=float, default=0.0)
    one.add_argument("--output")
    one.add_argument("--no-advanced", action="store_true")

    batch = sub.add_parser(
        "batch", help="Generate hypotheses for a validated JSON list of combinations"
    )
    add_context_arguments(batch)
    batch.add_argument("--input", required=True)
    batch.add_argument("--output", required=True)

    args = parser.parse_args()
    try:
        path_options = path_options_from_args(args)
    except ValueError as exc:
        parser.error(str(exc))
    rows: List[Dict[str, Any]] = []
    if args.command == "batch":
        rows = _load_batch(parser, args.input)
        if not args.demo:
            missing_cell_lines = [
                index for index, row in enumerate(rows)
                if not isinstance(row.get("cell_line"), str) or not row["cell_line"].strip()
            ]
            if missing_cell_lines:
                parser.error(
                    "production batch rows require an explicit non-empty cell_line; "
                    f"missing in rows {missing_cell_lines[:20]}"
                )
        queries = [
            (f"row {index} drug1", row["drug1_smiles"])
            for index, row in enumerate(rows)
        ] + [
            (f"row {index} drug2", row["drug2_smiles"])
            for index, row in enumerate(rows)
        ]
    else:
        if not args.demo and not str(args.cell_line).strip():
            parser.error("--cell-line is required with --context-config")
        queries = [(args.drug1, args.drug1_smiles), (args.drug2, args.drug2_smiles)]

    cache_dir = _preflight_cache(
        parser, args.plato_cache, queries, args.trusted_plato_sha256
    )
    pipeline_config: Dict[str, Any] = {
        "path_search": path_options,
        "run_mode": "research_hypothesis",
        "generated_from_test_fixture": False,
        "target_mode": "plato_only",
        "plato": {
            "cache_dir": str(cache_dir),
            "provider_release": args.plato_release,
            "trusted_export_sha256": args.trusted_plato_sha256,
        },
    }
    if args.quick:
        pipeline_config["network"] = {
            "randomizations": 100,
            "sensitivity_randomizations": 20,
        }
        pipeline_config["credibility"] = {"stability_randomizations": 25}
    try:
        pipe = ExplanationPipeline(config=pipeline_config)
        if args.demo:
            pipe.load_demo_resources()
        else:
            pipe.load_production_resources(args.context_config)
    except (FileNotFoundError, ImportError, OSError, RuntimeError, TypeError, ValueError, json.JSONDecodeError) as exc:
        parser.error(f"unable to initialize interpreter context: {exc}")

    if args.command == "batch":
        destination = _prepare_output(parser, args.output)
        try:
            results = pipe.explain_batch(rows, str(destination))
        except (KeyError, OSError, RuntimeError, TypeError, ValueError) as exc:
            parser.error(f"unable to explain batch: {exc}")
        for result in results:
            result.print_summary()
        return 0

    try:
        result = pipe.explain_combination(
            drug1_name=args.drug1,
            drug2_name=args.drug2,
            drug1_smiles=args.drug1_smiles,
            drug2_smiles=args.drug2_smiles,
            cell_line=args.cell_line,
            predicted_synergy=args.synergy,
            advanced=not args.no_advanced,
            disease_label=args.disease_label,
        )
    except (KeyError, OSError, RuntimeError, TypeError, ValueError) as exc:
        parser.error(f"unable to explain drug pair: {exc}")
    result.print_summary()
    if args.output:
        destination = _prepare_output(parser, args.output)
        try:
            pipe.save_result(result, destination)
        except OSError as exc:
            parser.error(f"unable to write output {destination}: {exc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
