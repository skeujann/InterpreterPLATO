"""Offline tests for scoped data-lock verification."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from InterpreterPLATO.scripts import bootstrap_data, verify_manifest


def _write_lock(path: Path, artifact_path: str, payload: bytes) -> None:
    path.write_text(json.dumps({
        "schema_version": "1.0",
        "artifacts": [{
            "relative_path": artifact_path,
            "size_bytes": len(payload),
            "sha256": hashlib.sha256(payload).hexdigest(),
        }],
    }), encoding="utf-8")


def test_default_scope_verifies_runtime_only(tmp_path: Path) -> None:
    data = tmp_path / "data"
    data.mkdir()
    runtime_payload = b"runtime"
    source_payload = b"source"
    (tmp_path / "runtime.bin").write_bytes(runtime_payload)
    (tmp_path / "source.bin").write_bytes(source_payload)
    _write_lock(data / "runtime.lock.json", "runtime.bin", runtime_payload)
    _write_lock(data / "source.lock.json", "source.bin", source_payload)

    assert verify_manifest.main(["--root", str(tmp_path)]) == 0
    (tmp_path / "source.bin").write_bytes(b"tampered-source")
    assert verify_manifest.main(["--root", str(tmp_path)]) == 0
    assert verify_manifest.main(["--root", str(tmp_path), "--scope", "source"]) == 1


def test_repository_locks_have_verifiable_artifacts() -> None:
    root = Path(__file__).resolve().parents[1]
    for name in ("runtime.lock.json", "source.lock.json"):
        manifest = verify_manifest.load_manifest(root / "data" / name)
        artifacts = list(verify_manifest.iter_artifacts(manifest))
        assert artifacts, name
        assert all("relative_path" in artifact for _, artifact in artifacts)


def test_lock_sizes_and_standalone_registry_are_consistent() -> None:
    root = Path(__file__).resolve().parents[1]
    runtime = verify_manifest.load_manifest(root / "data" / "runtime.lock.json")
    runtime_records = [item for _, item in verify_manifest.iter_artifacts(runtime)]
    assert runtime["manifest_kind"] == "runtime_integrity_lock"
    assert sum(int(item["size_bytes"]) for item in runtime_records) == runtime["total_size_bytes"]
    assert runtime["runtime_artifacts"][0]["relative_path"].endswith(".tsv.gz")

    source = verify_manifest.load_manifest(root / "data" / "source.lock.json")
    assert sum(int(item["size_bytes"]) for item in source["artifacts"]) == source["source_size_bytes"]
    assert sum(int(item["size_bytes"]) for item in source["build_registries"]) == source["build_registry_size_bytes"]

    registry = bootstrap_data.load_registry(root / "data" / "source_registry.json")
    resolved = bootstrap_data.resolve_inputs(registry, {})
    score_path, _ = resolved["scoreri_cell_lines"]
    assert score_path == (root / "data" / "scoreri_cell_lines.tsv").resolve()
