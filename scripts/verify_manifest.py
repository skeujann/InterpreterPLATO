"""Verify local artifacts recorded in a scoped data lock.

The verifier is intentionally read-only: it never downloads data or rewrites the
manifest.  Every manifest object that contains ``relative_path`` together with
``sha256`` and/or ``local_sha256`` is checked against the selected project root.
Runtime verification is the default so ordinary startup checks do not hash the
roughly 35 GiB source bundle.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Sequence, Tuple


DEFAULT_ROOT = Path(__file__).resolve().parents[1]
SCOPE_MANIFESTS = {
    "runtime": Path("data/runtime.lock.json"),
    "source": Path("data/source.lock.json"),
}
SHA256_RE = re.compile(r"^[0-9a-fA-F]{64}$")
HASH_CHUNK_SIZE = 8 * 1024 * 1024


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Read-only verification of files pinned in a provenance manifest."
    )
    parser.add_argument(
        "--root",
        type=Path,
        default=DEFAULT_ROOT,
        help="Project root used to resolve relative_path entries (default: InterpreterPLATO).",
    )
    parser.add_argument(
        "--scope",
        choices=("runtime", "source", "all"),
        default="runtime",
        help=(
            "Lock scope to verify (default: runtime). The source scope hashes "
            "the full multi-gigabyte source bundle."
        ),
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        help=(
            "Custom manifest path; relative paths are resolved under --root. "
            "When supplied, it overrides --scope."
        ),
    )
    return parser.parse_args(argv)


def iter_artifacts(node: Any, location: str = "$") -> Iterator[Tuple[str, Dict[str, Any]]]:
    """Yield every object carrying a local path and at least one local digest."""
    if isinstance(node, dict):
        has_digest = any(key in node for key in ("sha256", "local_sha256"))
        if "relative_path" in node and has_digest:
            yield location, node
        for key, value in node.items():
            yield from iter_artifacts(value, f"{location}.{key}")
    elif isinstance(node, list):
        for index, value in enumerate(node):
            yield from iter_artifacts(value, f"{location}[{index}]")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(HASH_CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve_artifact_path(root: Path, relative_path: str) -> Path:
    relative = Path(relative_path)
    if relative.is_absolute():
        raise ValueError("relative_path must not be absolute")

    resolved = (root / relative).resolve()
    try:
        resolved.relative_to(root)
    except ValueError as exc:
        raise ValueError("relative_path escapes --root") from exc
    return resolved


def verify_artifact(root: Path, location: str, artifact: Dict[str, Any]) -> List[str]:
    errors: List[str] = []
    raw_relative_path = artifact.get("relative_path")
    if not isinstance(raw_relative_path, str) or not raw_relative_path.strip():
        return [f"{location}: relative_path must be a non-empty string"]

    display_path = raw_relative_path.replace("\\", "/")
    expected_hashes: List[Tuple[str, str]] = []
    for key in ("sha256", "local_sha256"):
        if key not in artifact:
            continue
        value = artifact[key]
        if not isinstance(value, str) or not SHA256_RE.fullmatch(value):
            errors.append(f"{display_path}: {key} must be exactly 64 hexadecimal characters")
        else:
            expected_hashes.append((key, value.lower()))

    if not expected_hashes:
        if not errors:
            errors.append(f"{display_path}: no usable sha256/local_sha256 value")
        return errors

    try:
        path = resolve_artifact_path(root, raw_relative_path)
    except ValueError as exc:
        errors.append(f"{display_path}: {exc}")
        return errors

    if not path.is_file():
        errors.append(f"{display_path}: file is missing or is not a regular file")
        return errors

    try:
        before = path.stat()
        expected_size = artifact.get("size_bytes")
        if expected_size is not None:
            if isinstance(expected_size, bool) or not isinstance(expected_size, int) or expected_size < 0:
                errors.append(f"{display_path}: size_bytes must be a non-negative integer")
            elif before.st_size != expected_size:
                errors.append(
                    f"{display_path}: size mismatch (expected {expected_size}, got {before.st_size})"
                )

        actual_hash = sha256_file(path)
        after = path.stat()
    except OSError as exc:
        errors.append(f"{display_path}: could not read file ({exc})")
        return errors

    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        errors.append(f"{display_path}: file changed during verification")
        return errors

    for key, expected_hash in expected_hashes:
        if actual_hash != expected_hash:
            errors.append(
                f"{display_path}: {key} mismatch (expected {expected_hash}, got {actual_hash})"
            )
    return errors


def load_manifest(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _manifest_paths(root: Path, scope: str, manifest: Optional[Path]) -> List[Path]:
    if manifest is not None:
        selected = manifest.expanduser()
        return [(root / selected if not selected.is_absolute() else selected).resolve()]
    scopes = ("runtime", "source") if scope == "all" else (scope,)
    return [(root / SCOPE_MANIFESTS[item]).resolve() for item in scopes]


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    root = args.root.expanduser().resolve()

    if not root.is_dir():
        print(f"ERROR: --root is not a directory: {root}", file=sys.stderr)
        return 2

    manifest_paths = _manifest_paths(root, args.scope, args.manifest)
    total = 0
    failed = 0
    for manifest_path in manifest_paths:
        try:
            manifest = load_manifest(manifest_path)
        except (OSError, json.JSONDecodeError) as exc:
            print(f"ERROR: unable to load manifest {manifest_path}: {exc}", file=sys.stderr)
            return 2

        artifacts = list(iter_artifacts(manifest))
        if not artifacts:
            print(
                f"ERROR: manifest contains no verifiable relative_path/hash entries: {manifest_path}",
                file=sys.stderr,
            )
            return 2

        print(f"LOCK {manifest_path}")
        total += len(artifacts)
        for location, artifact in artifacts:
            errors = verify_artifact(root, location, artifact)
            display_path = str(artifact.get("relative_path", location)).replace("\\", "/")
            if errors:
                failed += 1
                for error in errors:
                    print(f"FAIL {error}")
            else:
                print(f"OK   {display_path}")

    passed = total - failed
    print(f"Verified {total} artifact record(s): {passed} passed, {failed} failed.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
