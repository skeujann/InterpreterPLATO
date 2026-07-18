"""Plan, verify, and rebuild InterpreterPLATO data without network access."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from typing import Any, Dict, Mapping, Optional, Sequence


PACKAGE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE_ROOT = PACKAGE_ROOT.parent
DEFAULT_REGISTRY = PACKAGE_ROOT / "data" / "source_registry.json"
CHUNK_SIZE = 8 * 1024 * 1024


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_registry(path: Path) -> Dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    inputs = payload.get("inputs")
    if not isinstance(inputs, list) or not inputs:
        raise ValueError("registry inputs must be a non-empty list")
    identifiers = [str(item.get("id", "")) for item in inputs if isinstance(item, dict)]
    if len(identifiers) != len(inputs) or any(not item for item in identifiers):
        raise ValueError("every registry input requires a non-empty id")
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("registry input ids must be unique")
    return payload


def _parse_overrides(values: Sequence[str]) -> Dict[str, Path]:
    overrides: Dict[str, Path] = {}
    for value in values:
        identifier, separator, raw_path = value.partition("=")
        if not separator or not identifier.strip() or not raw_path.strip():
            raise ValueError(f"--input must use ID=PATH: {value!r}")
        if identifier in overrides:
            raise ValueError(f"duplicate --input override: {identifier}")
        overrides[identifier] = Path(raw_path).expanduser().resolve()
    return overrides


def resolve_inputs(
    registry: Mapping[str, Any], overrides: Mapping[str, Path]
) -> Dict[str, tuple[Path, Mapping[str, Any]]]:
    resolved: Dict[str, tuple[Path, Mapping[str, Any]]] = {}
    known = {str(item["id"]) for item in registry["inputs"]}
    unknown = sorted(set(overrides) - known)
    if unknown:
        raise ValueError(f"unknown input override(s): {', '.join(unknown)}")
    for item in registry["inputs"]:
        identifier = str(item["id"])
        if identifier in overrides:
            path = overrides[identifier]
        else:
            base = WORKSPACE_ROOT if item.get("path_base") == "workspace" else PACKAGE_ROOT
            path = (base / str(item["default_path"])).resolve()
        resolved[identifier] = (path, item)
    return resolved


def verify_inputs(resolved: Mapping[str, tuple[Path, Mapping[str, Any]]]) -> bool:
    valid = True
    for identifier, (path, record) in resolved.items():
        if not path.is_file():
            print(f"MISSING {identifier}: {path}")
            valid = False
            continue
        expected_size = int(record["size_bytes"])
        actual_size = path.stat().st_size
        if actual_size != expected_size:
            print(f"FAIL {identifier}: size expected {expected_size}, got {actual_size}: {path}")
            valid = False
            continue
        actual_hash = _sha256(path)
        expected_hash = str(record["sha256"]).lower()
        if actual_hash != expected_hash:
            print(f"FAIL {identifier}: SHA-256 expected {expected_hash}, got {actual_hash}: {path}")
            valid = False
            continue
        print(f"OK   {identifier}: {path}")
    return valid


def _run(command: Sequence[str]) -> None:
    print("RUN  " + subprocess.list2cmdline(list(command)))
    subprocess.run(list(command), cwd=str(PACKAGE_ROOT), check=True)


def rebuild(resolved: Mapping[str, tuple[Path, Mapping[str, Any]]], component: str) -> None:
    paths = {identifier: str(value[0]) for identifier, value in resolved.items()}
    if component in {"context", "all"}:
        _run([
            sys.executable,
            str(PACKAGE_ROOT / "scripts" / "build_context_artifacts.py"),
            "all",
            "--hgnc", paths["hgnc"],
            "--string-links", paths["string_links"],
            "--string-info", paths["string_info"],
            "--depmap-expression", paths["depmap_expression"],
            "--depmap-model", paths["depmap_model"],
            "--depmap-profiles", paths["depmap_profiles"],
            "--score-combos", paths["scoreri_cell_lines"],
            "--reactome", paths["reactome"],
            "--cosmic", paths["cosmic"],
        ])
    if component in {"target-mapping", "all"}:
        _run([
            sys.executable,
            str(PACKAGE_ROOT / "scripts" / "build_chembl_target_mapping.py"),
            "--sqlite", paths["chembl_sqlite"],
            "--hgnc", paths["hgnc"],
            "--archive", paths["chembl_archive"],
            "--checksums", paths["chembl_checksums"],
        ])
    _run([
        sys.executable,
        str(PACKAGE_ROOT / "scripts" / "verify_manifest.py"),
        "--scope", "runtime",
    ])


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Explicit, offline source verification and runtime rebuild workflow"
    )
    parser.add_argument("command", choices=("plan", "verify", "rebuild"))
    parser.add_argument("--registry", type=Path, default=DEFAULT_REGISTRY)
    parser.add_argument(
        "--input", action="append", default=[], metavar="ID=PATH",
        help="Override one registry path; repeat for multiple inputs.",
    )
    parser.add_argument(
        "--component", choices=("context", "target-mapping", "all"), default="all",
        help="Component rebuilt by the rebuild command (default: all).",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(argv)
    try:
        registry = load_registry(args.registry.expanduser().resolve())
        overrides = _parse_overrides(args.input)
        resolved = resolve_inputs(registry, overrides)
        if args.command == "plan":
            for identifier, (path, record) in resolved.items():
                acquisition = record.get("acquisition_url")
                if not acquisition and record.get("path_base") == "package" and path.is_file():
                    acquisition = "included package-local registry"
                print(f"{identifier}\n  expected: {path}")
                print(f"  acquire:  {acquisition or 'provide explicitly'}")
                print(f"  policy:   {record.get('redistribution', 'review upstream terms')}")
            return 0
        if not verify_inputs(resolved):
            return 1
        if args.command == "rebuild":
            rebuild(resolved, args.component)
        return 0
    except (OSError, ValueError, json.JSONDecodeError, subprocess.CalledProcessError) as exc:
        print(f"ERROR: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
