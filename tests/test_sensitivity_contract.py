"""Fail-closed sensitivity checks, without real PLATO input or model claims."""

from __future__ import annotations

from copy import deepcopy
from typing import Any, Dict

import pytest

import InterpreterPLATO.interpreter as interpreter_module
from InterpreterPLATO import ExplanationPipeline
from InterpreterPLATO.advanced import AdvancedResult
from InterpreterPLATO.mechanism import MechanismResult


def _fidelity(supplied: Dict[str, Any]) -> Dict[str, Any]:
    pipeline = ExplanationPipeline()
    return pipeline._evaluate_model_fidelity(
        "drug-a", "drug-b", "CCO", "CCN", "A549", 10.0,
        AdvancedResult(), supplied,
    )


def _identified_perturbations() -> list:
    return [
        {"kind": "remove_target", "gene": "EGFR", "prediction_after": 8.0},
        {"kind": "remove_target", "gene": "MTOR", "prediction_after": 7.0},
        {"kind": "remove_module", "module": "signaling", "prediction_after": 6.0},
    ]


def test_caller_cannot_rebind_original_prediction() -> None:
    result = _fidelity({
        "baseline_prediction": 100.0,
        "provided_prediction": 100.0,
        "perturbations": _identified_perturbations(),
    })
    assert result["provided_prediction"] == 10.0
    assert result["caller_provided_prediction"] == 100.0
    assert result["baseline_difference"] == 90.0
    assert result["caller_prediction_consistent"] is False
    assert result["baseline_consistent"] is False
    assert result["meets_internal_sensitivity_thresholds"] is False
    assert result["interpretation_type"] == "post_hoc_biological_hypothesis"


def test_mismatching_claim_fails_even_when_baseline_matches() -> None:
    result = _fidelity({
        "baseline_prediction": 10.0,
        "provided_prediction": 100.0,
        "perturbations": _identified_perturbations(),
    })
    assert result["baseline_difference"] == 0.0
    assert result["baseline_consistent"] is False
    assert result["meets_internal_sensitivity_thresholds"] is False


def test_semantically_duplicate_perturbations_are_not_three_experiments() -> None:
    repeated = {"kind": "remove_target", "gene": "EGFR", "prediction_after": 8.0}
    result = _fidelity({
        "baseline_prediction": 10.0,
        "perturbations": [
            repeated,
            {"prediction_after": 8.0, "perturbation": {"gene": "EGFR", "kind": "remove_target"}},
            {**repeated, "perturbation_id": "different-label-same-gene"},
        ],
    })
    assert result["unique_perturbation_count"] == 1
    assert result["duplicate_perturbation_count"] == 2
    assert result["rejected_perturbations"] == 2
    assert result["meets_internal_sensitivity_thresholds"] is False


def test_kind_only_records_remain_descriptive_not_unique() -> None:
    result = _fidelity({
        "baseline_prediction": 10.0,
        "perturbations": [
            {"kind": kind, "prediction_after": 12.0}
            for kind in ("remove_target", "remove_module", "remove_path")
        ],
    })
    assert result["status"] == "evaluated"
    assert result["support_rate"] == 0.0
    assert len(result["perturbations"]) == 3
    assert result["unique_perturbation_count"] == 0
    assert result["unidentified_perturbation_count"] == 3
    assert result["meets_internal_sensitivity_thresholds"] is False
    assert result["evidence_source"] == "caller_supplied_unverified"


def test_identified_distinct_perturbations_can_pass_internal_checks_only() -> None:
    result = _fidelity({
        "baseline_prediction": 10.0,
        "perturbations": _identified_perturbations(),
    })
    assert result["unique_perturbation_count"] == 3
    assert result["perturbation_type_count"] == 2
    assert result["meets_internal_sensitivity_thresholds"] is True
    assert result["evidence_source"] == "caller_supplied_unverified"
    assert result["independently_validated"] is False


def _network_result(budget: int = 1000, **updates: Any) -> MechanismResult:
    cheng = {
        "valid": True,
        "quality_status": "adequate",
        "n_randomizations": budget,
        "n_valid_randomizations": budget,
        "random_seed": 42,
        "distance_method": "confidence_weighted_target_to_nearest_fixed_disease",
        "null_method": "degree_matched_fixed_disease_module_disjoint_complete_draws",
        "method_version": "weighted_fixed_disease_closest_v2",
        "null_A": {"std": 1.0},
        "null_B": {"std": 1.2},
        **updates,
    }
    return MechanismResult(
        mechanism_type="complementary_exposure", confidence=0.7,
        description="synthetic network mechanism", metrics={"cheng": cheng},
    )


def _stability(monkeypatch: pytest.MonkeyPatch, baseline: MechanismResult, reruns: list) -> Dict[str, Any]:
    pipeline = ExplanationPipeline(config={
        "credibility": {"stability_randomizations": 1000, "top_k_stability": [5, 20]},
    })
    outputs = iter(reruns)
    monkeypatch.setattr(interpreter_module, "classify_mechanism", lambda *args, **kwargs: next(outputs))
    return pipeline._top_k_stability(
        [object()] * 20, [object()] * 20, baseline, top_k=10, disease_label="test",
    )


def test_low_budget_baseline_cannot_borrow_rerun_budget(monkeypatch: pytest.MonkeyPatch) -> None:
    result = _stability(monkeypatch, _network_result(5), [_network_result(), _network_result()])
    assert result["descriptive_stable"] is True
    assert result["primary_baseline"]["network_randomizations"] == 5
    assert result["randomization_budget_sufficient"] is False
    assert result["protocol_budget_matches"] is False
    assert result["meets_full_robustness_criteria"] is False


@pytest.mark.parametrize("invalid", [
    {"valid": False},
    {"quality_status": "limited"},
    {"null_A": {"std": 0.0}},
    {"n_valid_randomizations": 999},
])
def test_one_invalid_null_prevents_full_stability(monkeypatch: pytest.MonkeyPatch, invalid: Dict[str, Any]) -> None:
    result = _stability(monkeypatch, _network_result(), [_network_result(**invalid), _network_result()])
    assert result["descriptive_stable"] is True
    assert result["network_null_quality_adequate"] is False
    assert result["meets_full_robustness_criteria"] is False


def test_adequate_but_different_protocols_are_not_full_stability(monkeypatch: pytest.MonkeyPatch) -> None:
    result = _stability(monkeypatch, _network_result(2000), [_network_result(), _network_result()])
    assert result["randomization_budget_sufficient"] is True
    assert result["network_null_quality_adequate"] is True
    assert result["common_network_protocol"] is False
    assert result["meets_full_robustness_criteria"] is False


def test_nonnetwork_mechanism_has_only_descriptive_stability(monkeypatch: pytest.MonkeyPatch) -> None:
    baseline = _network_result()
    baseline.mechanism_type = "shared_target_co_modulation"
    result = _stability(monkeypatch, baseline, [deepcopy(baseline), deepcopy(baseline)])
    assert result["descriptive_stable"] is True
    assert result["statistical_network_scope"] is False
    assert result["meets_full_robustness_criteria"] is False


def test_adequate_common_protocol_alone_is_not_statistical_perturbation(monkeypatch: pytest.MonkeyPatch) -> None:
    result = _stability(monkeypatch, _network_result(), [_network_result(), _network_result()])
    assert result["descriptive_stable"] is True
    assert result["randomization_budget_sufficient"] is True
    assert result["network_null_quality_adequate"] is True
    assert result["common_network_protocol"] is True
    assert result["meets_full_robustness_criteria"] is False
    assert result["statistical_target_set_variation"] is False
    assert "not_external_validation" in result["full_robustness_scope"]


def test_frozen_statistical_inputs_cannot_earn_network_robustness_points(monkeypatch: pytest.MonkeyPatch) -> None:
    baseline = _network_result(statistical_target_sets={"drug1": ["EGFR", "MTOR"], "drug2": ["RRM1"]})
    result = _stability(monkeypatch, baseline, [deepcopy(baseline), deepcopy(baseline)])
    assert result["descriptive_stable"] is True
    assert result["common_network_protocol"] is True
    assert result["statistical_target_set_variation"] is False
    assert result["meets_full_robustness_criteria"] is False
    assert "top_k_does_not_perturb_statistical_target_sets" in result["limitations"]


def test_missing_method_version_prevents_common_protocol(monkeypatch: pytest.MonkeyPatch) -> None:
    baseline = _network_result(method_version=None)
    result = _stability(monkeypatch, baseline, [deepcopy(baseline), deepcopy(baseline)])
    assert result["common_network_protocol"] is False
    assert result["meets_full_robustness_criteria"] is False
