"""Test-only source exceptions must not become integrity or scientific approvals."""
from copy import deepcopy
import json

import pytest

from InterpreterPLATO import development_inputs as dev
from InterpreterPLATO.engineering import file_sha
from InterpreterPLATO.targets import PlatoPredictionCache, inspect_plato_payload


@pytest.fixture
def test_input(tmp_path, monkeypatch):
    # Minimal synthetic POST response; the fixture schema is enabled only here.
    payload = [{"chemblId": "CHEMBL203", "nome": "Epidermal growth factor receptor:Homo sapiens",
                "score": 10., "simili": [], "tid": 1, "validato": True}]
    raw = tmp_path / "fixture.json"
    raw.write_text(json.dumps(payload), encoding="utf-8")
    fingerprint = inspect_plato_payload(payload)["schema_fingerprint"]
    monkeypatch.setattr(dev, "OBSERVED_SCHEMAS", frozenset({fingerprint}))
    cache = tmp_path / "cache"
    evidence = {"endpoint": "https://prometheus.farmacia.uniba.it/plato/pred.html",
                "http_method": "POST", "response_sha256": file_sha(raw), "submitted_smiles": "CCO",
                "form": {"tipoout": "json", "tipopred": "mussel"}}
    PlatoPredictionCache(cache).import_file("CCO", raw, "fixture", "r37", query_binding_evidence=evidence)
    inventory = {str(path.resolve()): file_sha(path) for path in cache.rglob("*") if path.is_file()}
    config = {"run_mode": "test", "generated_from_test_fixture": True,
              "plato": {"cache_dir": str(cache), "trusted_export_sha256": []},
              "test_source_review_exception": dev.test_policy(inventory)}
    pipe = dev.DevelopmentTestPipeline(config)
    fused, resolution = pipe.targets.resolve("fixture", "CCO")
    assert fused
    return pipe, fused, resolution


def test_approval_exceptions_are_local_and_do_not_change_source_status(test_input):
    pipe, fused, resolution = test_input
    before = deepcopy(resolution)
    strict = pipe._plato_cache_integrity_failures(resolution, "r37")
    assert strict and set(strict) <= dev.APPROVAL_FAILURES
    pipe._validate_plato_only_resolution("fixture", "CCO", fused, resolution)
    assert resolution == before
    assert pipe._plato_cache_integrity_failures(resolution, "r37") == strict
    assert pipe.review_waivers
    assert pipe.config["plato"]["trusted_export_sha256"] == []
    assert not pipe.config["test_source_review_exception"]["production_schema_approved"]


@pytest.mark.parametrize("change", ["identity", "raw_hash", "raw_bytes", "binding", "schema", "mode", "zero_targets"])
def test_data_failures_are_not_waived(test_input, change):
    pipe, fused, resolution = test_input
    resolution = deepcopy(resolution)
    cache = resolution["cache"]
    if change == "identity":
        cache["query"]["inchikey"] = "WRONG"
    elif change == "raw_hash":
        cache["source"]["sha256"] = "0" * 64
    elif change == "raw_bytes":
        from pathlib import Path
        path = Path(pipe.config["plato"]["cache_dir"]) / cache["source"]["raw_blob_relative_path"]
        path.write_bytes(b"[]")
    elif change == "binding":
        cache["query"]["single_query_evidence"]["submitted_smiles"] = "CCN"
    elif change == "schema":
        cache["schema_fingerprint"] = "unknown"
    elif change == "mode":
        pipe.config["run_mode"] = "research_hypothesis"
    elif change == "zero_targets":
        fused = []
    with pytest.raises(ValueError):
        pipe._validate_plato_only_resolution("fixture", "CCO", fused, resolution)
    assert pipe._test_validation_active is False


def test_non_test_construction_refused():
    with pytest.raises(ValueError, match="requires test mode"):
        dev.DevelopmentTestPipeline({"run_mode": "research_hypothesis"})
