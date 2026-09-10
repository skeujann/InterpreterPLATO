"""Regression tests for the thesis research-hypothesis output contract."""

from __future__ import annotations

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path

import pytest

from InterpreterPLATO import ExplanationPipeline
from InterpreterPLATO.research import stable_sha256
from InterpreterPLATO.scripts.test_support import (
    archive_verified_plato_cache,
    strict_pipeline,
)


def _test_result(tmp_path: Path):
    cache_dir = tmp_path / "plato_cache"
    cache_dir.mkdir()
    first_hash = archive_verified_plato_cache(cache_dir, "CCO", "drug-a")
    second_hash = archive_verified_plato_cache(
        cache_dir,
        "CCN",
        "drug-b",
        [{
            "target_symbol": "RRM1",
            "target_name": "Ribonucleotide reductase catalytic subunit M1",
            "organism": "Homo sapiens",
            "taxid": "9606",
            "target_type": "single protein",
            "score": 11.5,
            "reliable": "yes",
            "max_similarity": 0.77,
        }],
    )
    pipeline = strict_pipeline(cache_dir, [first_hash, second_hash])
    pipeline.load_demo_resources()
    result = pipeline.explain_combination(
        "drug-a",
        "drug-b",
        "CCO",
        "CCN",
        predicted_synergy=4.2,
        advanced=False,
    )
    return pipeline, result, first_hash, second_hash


def test_result_is_explicitly_scoped_and_auditable(tmp_path: Path) -> None:
    _, result, first_hash, second_hash = _test_result(tmp_path)
    payload = result.to_dict()

    assert payload["research_contract_version"] == "1.0"
    assert payload["run_mode"] == "test"
    assert payload["formal_result_eligible"] is False
    assert payload["generated_from_test_fixture"] is True
    assert "synthetic_test_fixture_not_for_research_conclusions" in (
        payload["formal_result_eligibility_reasons"]
    )
    assert payload["interpretation_scope"] == {
        "claim": "auditable_post_hoc_biological_mechanism_hypothesis",
        "not_a_model_causal_explanation": True,
        "not_experimental_validation": True,
        "not_clinical_efficacy_probability": True,
    }

    audit = payload["audit"]
    assert audit["run_mode"] == "test"
    assert audit["generated_from_test_fixture"] is True
    assert len(audit["code"]["snapshot"]["sha256"]) == 64
    assert audit["code"]["snapshot"]["file_count"] > 0
    assert len(audit["configuration"]["sha256"]) == 64
    assert audit["policies"]["target_ranking"]
    assert audit["policies"]["mechanism"]
    assert audit["policies"]["credibility"]
    assert audit["policies"]["plato_parser"]
    assert audit["plato_inputs"]["drug1"]["raw_response_sha256"] == first_hash
    assert audit["plato_inputs"]["drug2"]["raw_response_sha256"] == second_hash
    assert audit["randomization"]["random_seed"] == 7

    claimed_hash = audit["result_payload_sha256"]
    unhashed = deepcopy(payload)
    unhashed["audit"].pop("result_payload_sha256")
    assert claimed_hash == stable_sha256(unhashed)
    assert result.to_dict() == payload


def test_saved_result_has_exact_byte_hash_sidecar(tmp_path: Path) -> None:
    pipeline, result, _, _ = _test_result(tmp_path)
    destination = tmp_path / "result.json"

    record = pipeline.save_result(result, destination)
    output_bytes = destination.read_bytes()
    expected = sha256(output_bytes).hexdigest()

    assert record["sha256"] == expected
    assert json.loads(output_bytes)["formal_result_eligible"] is False
    assert Path(record["sha256_sidecar"]).read_text(encoding="utf-8") == (
        f"{expected}  {destination.name}\n"
    )


def test_invalid_research_mode_and_fixture_claims_fail_closed() -> None:
    with pytest.raises(ValueError, match="run_mode"):
        ExplanationPipeline(config={"run_mode": "production"})
    with pytest.raises(ValueError, match="generated_from_test_fixture"):
        ExplanationPipeline(config={
            "run_mode": "research_hypothesis",
            "generated_from_test_fixture": True,
        })
