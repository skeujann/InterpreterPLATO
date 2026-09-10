"""Sequentially capture two PLATO target-fishing responses and explain a pair."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import sys
import time
from typing import Any, Dict

PACKAGE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE.parent))

from InterpreterPLATO import ExplanationPipeline

try:  # Running as a script keeps this directory on sys.path.
    from .query_plato_once import (
        DEFAULT_ENDPOINT,
        DEFAULT_MAX_RESPONSE_BYTES,
        capture_plato_single_query,
    )
except ImportError:  # pragma: no cover - direct CLI execution path.
    from query_plato_once import (
        DEFAULT_ENDPOINT,
        DEFAULT_MAX_RESPONSE_BYTES,
        capture_plato_single_query,
    )


def _slug(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9]+", "_", value).strip("_").lower()
    return cleaned or "compound"


def _path_text(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value.resolve())
    if isinstance(value, dict):
        return {key: _path_text(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_path_text(item) for item in value]
    return value


def _normalize_drug_name(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", value.strip().lower()).strip()


def _normalize_trusted_hashes(values: list[str]) -> list[str]:
    normalized: set[str] = set()
    for value in values:
        digest = str(value or "").strip().lower()
        if not re.fullmatch(r"[0-9a-f]{64}", digest):
            raise ValueError(
                "--trusted-plato-sha256 values must be complete 64-character SHA-256 digests"
            )
        normalized.add(digest)
    return sorted(normalized)


def _read_smiles_rows(path: Path) -> list[tuple[str, str]]:
    if not path.is_file():
        return []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return []
        fields = {field.strip().lower(): field for field in reader.fieldnames}
        name_field = next(
            (fields[key] for key in ("drug_name", "name", "drug", "compound") if key in fields),
            None,
        )
        smiles_field = next(
            (fields[key] for key in ("smiles", "canonical_smiles", "canonical_smile") if key in fields),
            None,
        )
        if not name_field or not smiles_field:
            return []
        rows: list[tuple[str, str]] = []
        for row in reader:
            name = str(row.get(name_field, "")).strip()
            smiles = str(row.get(smiles_field, "")).strip()
            if name and smiles:
                rows.append((name, smiles))
        return rows


def resolve_smiles_from_local_tables(
    drug_name: str,
    lookup_files: list[Path] | None = None,
) -> dict[str, Any]:
    """Resolve a drug name to one local SMILES, failing closed on ambiguity."""
    normalized = _normalize_drug_name(drug_name)
    if not normalized:
        raise ValueError("drug name is empty")
    hits: list[dict[str, str]] = []
    # Rebuilds and formal runs must not depend on an ambient workspace layout.
    # Only caller-declared lookup tables are consulted.
    for path in lookup_files or []:
        for observed_name, smiles in _read_smiles_rows(path):
            if _normalize_drug_name(observed_name) == normalized:
                hits.append({
                    "source": str(path),
                    "drug_name": observed_name,
                    "smiles": smiles,
                })
    if not hits:
        raise ValueError(
            f"no SMILES found for {drug_name!r}; provide the matching --drug*-smiles "
            "value or an explicit --smiles-lookup file"
        )
    distinct: dict[str, dict[str, str]] = {}
    for hit in hits:
        distinct.setdefault(hit["smiles"], hit)
    if len(distinct) > 1:
        sources = "; ".join(
            f"{hit['drug_name']} -> {hit['smiles']} ({hit['source']})"
            for hit in distinct.values()
        )
        raise ValueError(f"ambiguous local SMILES for {drug_name!r}: {sources}")
    selected = next(iter(distinct.values()))
    return {
        "drug_name": drug_name,
        "smiles": selected["smiles"],
        "source": selected["source"],
        "matched_name": selected["drug_name"],
        "hit_count": len(hits),
    }


def _write_json(path: Path, value: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(_path_text(value), ensure_ascii=False, indent=2, allow_nan=False),
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Run the pair-level workflow: PLATO single query for drug 1, "
            "PLATO single query for drug 2, then an audited research mechanism hypothesis."
        )
    )
    parser.add_argument("--drug1", required=True, help="First drug name")
    parser.add_argument("--drug2", required=True, help="Second drug name")
    parser.add_argument("--drug1-smiles", help="Exact SMILES submitted for drug 1")
    parser.add_argument("--drug2-smiles", help="Exact SMILES submitted for drug 2")
    parser.add_argument(
        "--smiles-lookup",
        action="append",
        default=[],
        help=(
            "CSV with drug_name/name and smiles columns. Repeatable. "
            "Required for name lookup when an exact --drug*-smiles value is "
            "omitted; no SMILES table is bundled or discovered implicitly."
        ),
    )
    parser.add_argument(
        "--resolve-only",
        action="store_true",
        help="Resolve drug names to SMILES and stop before contacting PLATO",
    )
    parser.add_argument("--endpoint", default=DEFAULT_ENDPOINT)
    parser.add_argument("--plato-release", default="r37")
    parser.add_argument("--timeout", type=int, default=60)
    parser.add_argument(
        "--max-response-bytes",
        type=int,
        default=DEFAULT_MAX_RESPONSE_BYTES,
        help="Hard streaming limit applied independently to each PLATO response",
    )
    parser.add_argument(
        "--trusted-plato-sha256",
        action="append",
        default=[],
        help=(
            "SHA-256 from an independently authenticated PLATO export allowlist, "
            "obtained before this capture. Repeat for all accepted pair exports. "
            "Same-run digests are capture checksums, not trust anchors."
        ),
    )
    parser.add_argument(
        "--delay-seconds",
        type=float,
        default=3.0,
        help="Polite delay between the two sequential PLATO form submissions",
    )
    parser.add_argument(
        "--context-config",
        default=str(PACKAGE / "data" / "production_context.json"),
        help="Hash-pinned production context JSON",
    )
    parser.add_argument("--demo", action="store_true", help="Use bundled demo context")
    parser.add_argument(
        "--cell-line",
        default="",
        help="Explicit production expression cell line (required unless --demo)",
    )
    parser.add_argument("--disease-label", default="")
    parser.add_argument("--synergy", type=float, default=0.0)
    parser.add_argument("--no-advanced", action="store_true")
    parser.add_argument("--quick", action="store_true")
    parser.add_argument(
        "--run-dir",
        help=(
            "Run-specific output directory. Defaults to "
            "InterpreterPLATO/outputs/research/plato_pair_runs/<timestamp>_<drug1>__<drug2>"
        ),
    )
    parser.add_argument("--output", help="Explanation JSON output path")
    args = parser.parse_args()

    if args.delay_seconds < 0 or args.delay_seconds > 60:
        parser.error("--delay-seconds must be between 0 and 60")
    try:
        trusted_export_sha256 = _normalize_trusted_hashes(args.trusted_plato_sha256)
    except ValueError as exc:
        parser.error(str(exc))
    lookup_files = [Path(item).expanduser() for item in args.smiles_lookup]
    try:
        drug1_resolution = (
            {
                "drug_name": args.drug1,
                "smiles": args.drug1_smiles,
                "source": "cli_argument",
                "matched_name": args.drug1,
                "hit_count": 1,
            }
            if args.drug1_smiles
            else resolve_smiles_from_local_tables(args.drug1, lookup_files or None)
        )
        drug2_resolution = (
            {
                "drug_name": args.drug2,
                "smiles": args.drug2_smiles,
                "source": "cli_argument",
                "matched_name": args.drug2,
                "hit_count": 1,
            }
            if args.drug2_smiles
            else resolve_smiles_from_local_tables(args.drug2, lookup_files or None)
        )
    except ValueError as exc:
        parser.error(str(exc))
    drug1_smiles = str(drug1_resolution["smiles"])
    drug2_smiles = str(drug2_resolution["smiles"])
    if args.resolve_only:
        print(json.dumps(
            {
                "drug1": drug1_resolution,
                "drug2": drug2_resolution,
            },
            ensure_ascii=False,
            indent=2,
            allow_nan=False,
        ))
        return 0
    if not args.demo and not str(args.cell_line).strip():
        parser.error("--cell-line is required for production context; use --demo only for demo data")

    captured_at = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    run_dir = Path(
        args.run_dir
        or PACKAGE
        / "outputs"
        / "research"
        / "plato_pair_runs"
        / f"{captured_at.replace(':', '').replace('+', 'Z')}_{_slug(args.drug1)}__{_slug(args.drug2)}"
    ).expanduser()
    raw_dir = run_dir / "raw_plato"
    cache_dir = run_dir / "plato_cache"
    output_path = Path(args.output).expanduser() if args.output else run_dir / "explanation.json"
    manifest_path = run_dir / "run_manifest.json"

    try:
        print(f"[1/3] Querying PLATO for {args.drug1} ...")
        first = capture_plato_single_query(
            drug_name=args.drug1,
            smiles=drug1_smiles,
            endpoint=args.endpoint,
            release=args.plato_release,
            timeout=args.timeout,
            max_response_bytes=args.max_response_bytes,
            raw_dir=raw_dir / "drug1",
            cache_dir=cache_dir,
        )
        if args.delay_seconds:
            print(f"[1/3] Waiting {args.delay_seconds:g}s before the second query ...")
            time.sleep(args.delay_seconds)
        print(f"[2/3] Querying PLATO for {args.drug2} ...")
        second = capture_plato_single_query(
            drug_name=args.drug2,
            smiles=drug2_smiles,
            endpoint=args.endpoint,
            release=args.plato_release,
            timeout=args.timeout,
            max_response_bytes=args.max_response_bytes,
            raw_dir=raw_dir / "drug2",
            cache_dir=cache_dir,
        )
    except Exception as exc:  # noqa: BLE001 - CLI should preserve raw evidence and stop.
        parser.error(str(exc))

    capture_sha256 = sorted({first["raw_sha256"], second["raw_sha256"]})
    missing_trust_anchors = sorted(set(capture_sha256) - set(trusted_export_sha256))
    trust_anchor_match = not missing_trust_anchors and bool(capture_sha256)
    manifest: Dict[str, Any] = {
        "schema_version": "1.0",
        "workflow": "sequential_plato_single_query_pair_then_explain",
        "captured_at_utc": captured_at,
        "run_dir": run_dir,
        "raw_dir": raw_dir,
        "cache_dir": cache_dir,
        "query_order": [args.drug1, args.drug2],
        "delay_seconds": args.delay_seconds,
        "plato_release": args.plato_release,
        "capture_sha256": capture_sha256,
        "capture_sha256_semantics": "same_run_integrity_checksum_not_a_trust_anchor",
        "trusted_export_sha256": trusted_export_sha256,
        "trusted_export_sha256_source": (
            "caller_supplied_external_allowlist" if trusted_export_sha256 else "not_configured"
        ),
        "trust_anchor_independence_status": (
            "caller_asserted_external_not_verifiable_by_helper"
            if trusted_export_sha256 else "missing"
        ),
        "trusted_export_sha256_match": trust_anchor_match,
        "missing_trusted_export_sha256": missing_trust_anchors,
        "content_authentication_status": (
            "caller_supplied_external_hash_allowlist_match"
            if trust_anchor_match else "same_run_capture_not_independently_authenticated"
        ),
        "smiles_resolution": {
            "drug1": drug1_resolution,
            "drug2": drug2_resolution,
        },
        "plato_results": [first, second],
        "context_config": args.context_config,
        "demo": bool(args.demo),
        "explanation_output": output_path,
        "run_mode": "demo" if args.demo else "research_hypothesis",
        "formal_result_eligible": False,
        "formal_result_eligibility_reasons": [
            "thesis_research_hypothesis_only",
            *(["demo_context_not_for_formal_results"] if args.demo else []),
        ],
    }
    _write_json(manifest_path, manifest)

    if not trust_anchor_match:
        manifest["explanation_status"] = "not_run_without_external_trust_anchor"
        manifest["formal_result_eligible"] = False
        manifest["formal_result_eligibility_reasons"] = [
            "thesis_research_hypothesis_only",
            "missing_independent_plato_trust_anchor",
        ]
        manifest["explanation_error"] = (
            "Research-hypothesis generation requires every captured PLATO response SHA-256 "
            "to match a caller-supplied, independently authenticated allowlist."
        )
        _write_json(manifest_path, manifest)
        print(
            "PLATO captures were archived, but research-hypothesis generation was not run: "
            "the same-run response digests are not independent trust anchors.",
            file=sys.stderr,
        )
        print(manifest_path)
        return 2

    pipeline_config: Dict[str, Any] = {
        "run_mode": "research_hypothesis",
        "generated_from_test_fixture": False,
        "target_mode": "plato_only",
        "plato": {
            "cache_dir": str(cache_dir),
            "provider_release": args.plato_release,
            "trusted_export_sha256": trusted_export_sha256,
        },
    }
    if args.quick:
        pipeline_config["network"] = {
            "randomizations": 100,
            "sensitivity_randomizations": 20,
        }
        pipeline_config["credibility"] = {"stability_randomizations": 25}

    try:
        print("[3/3] Generating the drug-pair mechanism hypothesis ...")
        pipe = ExplanationPipeline(config=pipeline_config)
        if args.demo:
            pipe.load_demo_resources()
        else:
            pipe.load_production_resources(args.context_config)
        result = pipe.explain_combination(
            drug1_name=args.drug1,
            drug2_name=args.drug2,
            drug1_smiles=drug1_smiles,
            drug2_smiles=drug2_smiles,
            cell_line=args.cell_line,
            predicted_synergy=args.synergy,
            advanced=not args.no_advanced,
            disease_label=args.disease_label,
        )
        result.print_summary()
        output_record = pipe.save_result(result, output_path)
    except Exception as exc:  # noqa: BLE001 - raw/caches already preserved.
        manifest["explanation_status"] = "failed"
        manifest["explanation_error"] = str(exc)
        _write_json(manifest_path, manifest)
        parser.error(
            f"PLATO captures succeeded and manifest was saved at {manifest_path}, "
            f"but explanation failed: {exc}"
        )

    manifest["explanation_status"] = "research_hypothesis_generated"
    manifest["run_mode"] = "demo" if args.demo else "research_hypothesis"
    manifest["formal_result_eligible"] = False
    manifest["formal_result_eligibility_reasons"] = [
        "thesis_research_hypothesis_only",
        *(["demo_context_not_for_formal_results"] if args.demo else []),
    ]
    manifest["explanation_output_sha256"] = output_record["sha256"]
    manifest["explanation_output_sha256_sidecar"] = output_record["sha256_sidecar"]
    _write_json(manifest_path, manifest)
    print(output_path)
    print(manifest_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
