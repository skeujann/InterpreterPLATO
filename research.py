"""Thesis-oriented provenance and output-contract helpers.

The research contract deliberately labels every current result as a biological
mechanism hypothesis.  It does not weaken any PLATO integrity gate and it does
not create a production approval path.
"""

from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
from importlib import metadata as importlib_metadata
from pathlib import Path
from typing import Any, Dict, Mapping
import json
import os
import platform
import subprocess
import sys
import tempfile

from .mechanism import NETWORK_METHOD_VERSION, STATISTICAL_CONTRACT_VERSION


RESEARCH_CONTRACT_VERSION = "1.0"
DEFAULT_RUN_MODE = "research_hypothesis"
ALLOWED_RUN_MODES = {DEFAULT_RUN_MODE, "demo", "test"}
TARGET_RANKING_POLICY_VERSION = "context_target_ranking_v1"
MECHANISM_POLICY_VERSION = "rule_based_mechanism_v2_frozen_statistical_targets"


def canonical_json_bytes(value: Any) -> bytes:
    """Serialize JSON-like data deterministically for audit hashing."""
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
        default=str,
    ).encode("utf-8")


def stable_sha256(value: Any) -> str:
    return sha256(canonical_json_bytes(value)).hexdigest()


def _code_snapshot(package_root: Path) -> Dict[str, Any]:
    """Hash runtime Python sources, including dirty working-tree contents."""
    candidates = list(package_root.glob("*.py"))
    scripts = package_root / "scripts"
    if scripts.is_dir():
        candidates.extend(scripts.rglob("*.py"))
    files = sorted({path.resolve() for path in candidates if path.is_file()})
    digest = sha256()
    records = []
    for path in files:
        relative = path.relative_to(package_root.resolve()).as_posix()
        content = path.read_bytes()
        content_sha256 = sha256(content).hexdigest()
        digest.update(relative.encode("utf-8"))
        digest.update(b"\0")
        digest.update(content)
        digest.update(b"\0")
        records.append({"path": relative, "sha256": content_sha256})
    return {
        "algorithm": "sha256_ordered_relative_path_and_file_bytes_v1",
        "sha256": digest.hexdigest(),
        "file_count": len(records),
        "files": records,
    }


def _git_identity(package_root: Path) -> Dict[str, Any]:
    """Best-effort Git identity; source hashing remains authoritative locally."""
    configured = str(os.environ.get("INTERPRETERPLATO_CODE_COMMIT", "")).strip()
    result: Dict[str, Any] = {
        "commit": configured,
        "commit_source": "environment" if configured else "unavailable",
        "working_tree_dirty": None,
    }
    try:
        commit_process = subprocess.run(
            ["git", "-C", str(package_root), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=3,
            check=False,
        )
        if commit_process.returncode == 0 and commit_process.stdout.strip():
            result["commit"] = commit_process.stdout.strip()
            result["commit_source"] = "git"
        dirty_process = subprocess.run(
            [
                "git", "-C", str(package_root), "status", "--porcelain",
                "--untracked-files=normal", "--", ".",
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
            timeout=3,
            check=False,
        )
        if dirty_process.returncode == 0:
            result["working_tree_dirty"] = bool(dirty_process.stdout.strip())
    except (OSError, subprocess.SubprocessError):
        pass
    return result


def _dependency_versions() -> Dict[str, str]:
    distributions = {
        "interpreterplato": "interpreterplato",
        "networkx": "networkx",
        "numpy": "numpy",
        "pandas": "pandas",
        "PyYAML": "PyYAML",
        "rdkit": "rdkit",
    }
    versions: Dict[str, str] = {}
    for label, distribution in distributions.items():
        try:
            versions[label] = importlib_metadata.version(distribution)
        except importlib_metadata.PackageNotFoundError:
            versions[label] = "not_installed_as_distribution"
    return versions


def plato_resolution_audit(resolution: Mapping[str, Any]) -> Dict[str, Any]:
    cache = resolution.get("cache", {})
    cache = cache if isinstance(cache, Mapping) else {}
    source = cache.get("source", {})
    source = source if isinstance(source, Mapping) else {}
    query = cache.get("query", {})
    query = query if isinstance(query, Mapping) else {}
    return {
        "provider": str(resolution.get("provider", "")),
        "release": str(resolution.get("provider_release", "")),
        "cache_schema_version": str(cache.get("cache_schema_version", "")),
        "parser_version": str(cache.get("parser_version", "")),
        "adapter": str(cache.get("adapter", "")),
        "schema_fingerprint": str(cache.get("schema_fingerprint", "")),
        "schema_fixture_verified": cache.get("schema_fixture_verified") is True,
        "query_binding_status": str(query.get("query_binding_status", "")),
        "canonical_isomeric_smiles": str(query.get("canonical_isomeric_smiles", "")),
        "inchikey": str(query.get("inchikey", "")),
        "raw_response_sha256": str(source.get("sha256", "")),
        "payload_sha256": str(source.get("payload_sha256", "")),
        "normalized_targets_sha256": str(source.get("normalized_targets_sha256", "")),
        "audit_metadata_sha256": str(source.get("audit_metadata_sha256", "")),
    }


def build_research_audit(
    *,
    package_root: Path,
    config: Mapping[str, Any],
    run_mode: str,
    generated_from_test_fixture: bool,
    generated_at_utc: str,
    runtime_manifest: Mapping[str, Any],
    artifact_inventory: Mapping[str, Any],
    plato_inputs: Mapping[str, Any],
    credibility_formula_version: str,
    plato_parser_version: str,
) -> Dict[str, Any]:
    """Build a self-contained audit record for one research result."""
    return {
        "research_contract_version": RESEARCH_CONTRACT_VERSION,
        "generated_at_utc": generated_at_utc,
        "run_mode": run_mode,
        "generated_from_test_fixture": bool(generated_from_test_fixture),
        "code": {
            "git": _git_identity(package_root),
            "snapshot": _code_snapshot(package_root),
        },
        "configuration": {
            "sha256": stable_sha256(config),
            "hash_semantics": "canonical_json_of_effective_runtime_configuration",
        },
        "policies": {
            "statistical_contract": STATISTICAL_CONTRACT_VERSION,
            "network_statistics": NETWORK_METHOD_VERSION,
            "target_ranking": TARGET_RANKING_POLICY_VERSION,
            "mechanism": MECHANISM_POLICY_VERSION,
            "credibility": credibility_formula_version,
            "plato_parser": plato_parser_version,
        },
        "runtime": {
            "python_version": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "platform": platform.platform(),
            "executable": str(Path(sys.executable).resolve()),
            "dependencies": _dependency_versions(),
            "runtime_manifest": deepcopy(dict(runtime_manifest)),
            "artifact_inventory": deepcopy(dict(artifact_inventory)),
        },
        "plato_inputs": deepcopy(dict(plato_inputs)),
        "randomization": {
            "random_seed": config.get("network", {}).get("random_seed"),
            "network_randomizations": config.get("network", {}).get("randomizations"),
            "sensitivity_randomizations": config.get("network", {}).get("sensitivity_randomizations"),
            "stability_randomizations": config.get("credibility", {}).get("stability_randomizations"),
        },
    }


def attach_result_payload_sha256(payload: Mapping[str, Any]) -> Dict[str, Any]:
    """Hash the complete canonical result except for the hash value itself."""
    result = deepcopy(dict(payload))
    raw_audit = result.get("audit", {})
    audit = deepcopy(dict(raw_audit)) if isinstance(raw_audit, Mapping) else {}
    audit["result_payload_sha256_scope"] = (
        "canonical_json_of_full_result_excluding_audit.result_payload_sha256"
    )
    audit.pop("result_payload_sha256", None)
    result["audit"] = audit
    audit["result_payload_sha256"] = stable_sha256(result)
    return result


def write_json_with_sha256(path: str | Path, payload: Any) -> Dict[str, Any]:
    """Atomically write pretty JSON and an exact-byte SHA-256 sidecar."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    data = json.dumps(
        payload,
        ensure_ascii=False,
        indent=2,
        allow_nan=False,
    ).encode("utf-8")
    digest = sha256(data).hexdigest()
    sidecar = Path(f"{destination}.sha256")

    def atomic_write(target: Path, content: bytes) -> None:
        handle = tempfile.NamedTemporaryFile(
            mode="wb",
            dir=target.parent,
            prefix=f".{target.name}.",
            suffix=".tmp",
            delete=False,
        )
        temporary = Path(handle.name)
        try:
            with handle:
                handle.write(content)
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(temporary, target)
        except Exception:
            try:
                temporary.unlink(missing_ok=True)
            except OSError:
                pass
            raise

    atomic_write(destination, data)
    atomic_write(sidecar, f"{digest}  {destination.name}\n".encode("utf-8"))
    return {
        "path": str(destination.resolve()),
        "sha256": digest,
        "size_bytes": len(data),
        "sha256_sidecar": str(sidecar.resolve()),
    }
