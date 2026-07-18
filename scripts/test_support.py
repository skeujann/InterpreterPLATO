"""Shared helpers for strict PLATO-only regression tests."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import InterpreterPLATO.targets as target_module
from InterpreterPLATO import ExplanationPipeline
from InterpreterPLATO.targets import PlatoPredictionCache, inspect_plato_payload, parse_plato_json


DEFAULT_TARGETS = (
    {
        "target_symbol": "EGFR",
        "target_name": "Epidermal growth factor receptor",
        "organism": "Homo sapiens",
        "taxid": "9606",
        "target_type": "single protein",
        "score": 12.4,
        "reliable": "yes",
        "max_similarity": 0.82,
        "action": "inhibitor",
    },
)


def plato_payload(
    smiles: str,
    target_records: Sequence[Mapping[str, Any]] = DEFAULT_TARGETS,
    release: str = "r37",
) -> dict[str, Any]:
    return {
        "provider": "PLATO",
        "release": release,
        "method": "target_fishing",
        "query_smiles": smiles,
        "targets": [dict(item) for item in target_records],
    }


def allow_test_schema(payload: Mapping[str, Any]) -> str:
    """Allowlist one schema fingerprint for this test process only."""
    report = inspect_plato_payload(payload, parse_plato_json(payload))
    fingerprint = str(report["schema_fingerprint"])
    fingerprints = set(target_module.VERIFIED_PLATO_SCHEMA_FINGERPRINTS)
    fingerprints.add(fingerprint)
    target_module.VERIFIED_PLATO_SCHEMA_FINGERPRINTS = tuple(sorted(fingerprints))
    return fingerprint


def write_payload(path: Path, payload: Mapping[str, Any]) -> str:
    data = json.dumps(
        payload, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False,
    ).encode("utf-8")
    path.write_bytes(data)
    return hashlib.sha256(data).hexdigest()


def archive_verified_plato_cache(
    cache_dir: Path,
    smiles: str,
    drug_name: str,
    target_records: Sequence[Mapping[str, Any]] = DEFAULT_TARGETS,
    release: str = "r37",
) -> str:
    payload = plato_payload(smiles, target_records, release=release)
    allow_test_schema(payload)
    source = cache_dir / f"{drug_name}.json"
    source_sha256 = write_payload(source, payload)
    cache = PlatoPredictionCache(cache_dir)
    cache.import_file(
        smiles=smiles,
        result_file=source,
        drug_name=drug_name,
        provider_release=release,
        confirmed_real_export=True,
    )
    return source_sha256


def strict_pipeline(cache_dir: Path, trusted_hashes: Iterable[str], **config: Any) -> ExplanationPipeline:
    return ExplanationPipeline(config={
        "plato": {
            "cache_dir": str(cache_dir),
            "trusted_export_sha256": list(trusted_hashes),
        },
        "advanced": False,
        "network": {
            "randomizations": 5,
            "sensitivity_randomizations": 3,
            "random_seed": 7,
        },
        "credibility": {
            "stability_randomizations": 3,
        },
        **config,
    })
