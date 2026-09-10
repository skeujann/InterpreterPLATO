"""Offline protocol and isolation checks for interpreter-only sensitivity."""

from __future__ import annotations

from copy import deepcopy
from typing import Any

import pytest

import InterpreterPLATO.advanced as advanced
from InterpreterPLATO.context import CancerContext, ContextScore
from InterpreterPLATO.mechanism import MechanismResult


def _result(mechanism: str = "complementary_exposure", **updates: Any) -> MechanismResult:
    return MechanismResult(
        mechanism_type=mechanism, confidence=0.7 if mechanism == "complementary_exposure" else 0.5,
        description="synthetic interpreter result",
        metrics={"cheng": {
            "valid": True,
            "quality_status": "adequate",
            "n_randomizations": 1000,
            "n_valid_randomizations": 1000,
            "random_seed": 42,
            "method_version": "weighted_fixed_disease_closest_v2",
            "distance_method": "confidence_weighted_target_to_nearest_fixed_disease",
            "null_method": "degree_matched_fixed_disease_module_disjoint_complete_draws",
            "null_A": {"std": 1.0},
            "null_B": {"std": 1.2},
            **updates,
        }},
    )


def _compare(baseline: MechanismResult, perturbed: MechanismResult, **options: Any) -> dict:
    return advanced._counterfactual_comparison(
        baseline, perturbed, **{"network_randomizations": 1000, "random_seed": 42, **options},
    )


def _target(symbol: str) -> ContextScore:
    return ContextScore(
        target_id=symbol, target_name=symbol, final_score=1.0, evidence_score=1.0,
        expression_score=None, disease_proximity_score=0.0, pathway_score=1.0,
        centrality_score=0.0, evidence_status="predicted_only", source="PLATO",
        source_evidence={},
    )


@pytest.mark.parametrize("side", ["baseline", "perturbed"])
@pytest.mark.parametrize("updates,reason", [
    ({"valid": False}, "invalid_network_statistics"),
    ({"quality_status": "limited"}, "inadequate_null_quality"),
    ({"n_randomizations": 10}, "insufficient_requested_randomizations"),
    ({"n_valid_randomizations": 999}, "insufficient_valid_randomizations"),
    ({"n_valid_randomizations": 1001}, "invalid_acceptance_fraction"),
    ({"n_randomizations": 2000, "n_valid_randomizations": 1000}, "invalid_acceptance_fraction"),
    ({"n_valid_randomizations": float("nan")}, "insufficient_valid_randomizations"),
    ({"null_A": {"std": 0.0}}, "degenerate_or_missing_null_variance"),
    ({"null_B": {"std": float("inf")}}, "degenerate_or_missing_null_variance"),
    ({"null_A": None}, "degenerate_or_missing_null_variance"),
    ({"random_seed": 43}, "configured_seed_mismatch"),
    ({"method_version": ""}, "missing_method_version"),
    ({"distance_method": None}, "missing_distance_method"),
    ({"null_method": ""}, "missing_null_method"),
])
def test_actual_null_failure_prevents_evaluable_comparison(side: str, updates: dict, reason: str) -> None:
    results = {"baseline": _result(), "perturbed": _result("overlapping_exposure")}
    results[side].metrics["cheng"].update(updates)
    comparison = _compare(**results)
    assert comparison["comparison_evaluable"] is False
    assert comparison["comparison_status"] == "not_evaluable"
    assert comparison["mechanism_changed"] is None
    assert comparison["descriptive_mechanism_changed"] is True
    assert f"{side}:{reason}" in comparison["comparison_audit"]["reason_codes"]


@pytest.mark.parametrize("field,value", [
    ("method_version", "different_version"),
    ("distance_method", "different_distance"),
    ("null_method", "different_null"),
    ("random_seed", 43),
    ("n_randomizations", 2000),
])
def test_individually_adequate_but_different_protocols_fail(field: str, value: Any) -> None:
    comparison = _compare(_result(), _result(**{field: value}))
    assert comparison["comparison_evaluable"] is False
    assert f"protocol_mismatch:{field}" in comparison["comparison_audit"]["reason_codes"]


def test_identical_actual_protocol_must_also_match_sensitivity_configuration() -> None:
    result = _result(n_randomizations=2000, n_valid_randomizations=2000)
    comparison = _compare(result, deepcopy(result))
    assert comparison["comparison_evaluable"] is False
    assert "baseline:configured_budget_mismatch" in comparison["comparison_audit"]["reason_codes"]
    assert "perturbed:configured_budget_mismatch" in comparison["comparison_audit"]["reason_codes"]


def test_no_network_statistics_is_not_an_evaluable_comparison() -> None:
    result = _result()
    result.metrics = {}
    comparison = _compare(result, deepcopy(result))
    assert comparison["comparison_evaluable"] is False
    assert comparison["descriptive_mechanism_changed"] is False


@pytest.mark.parametrize("changed", [True, False])
def test_adequate_matched_realizations_allow_interpreter_comparison_only(changed: bool) -> None:
    comparison = _compare(
        _result(n_randomizations=2000, n_valid_randomizations=1800),
        _result("overlapping_exposure" if changed else "complementary_exposure",
                n_randomizations=2000, n_valid_randomizations=1900),
        network_randomizations=2000,
    )
    assert comparison["comparison_evaluable"] is True
    assert comparison["mechanism_changed"] is changed
    assert comparison["comparison_audit"]["reason_codes"] == []
    assert comparison["comparison_audit"]["baseline"]["valid_fraction"] == 0.9
    assert comparison["comparison_audit"]["scope"] == "interpreter_sensitivity_not_model_fidelity_or_biological_causality"


def test_target_and_module_reports_keep_descriptive_changes_without_borrowing_budget(monkeypatch: pytest.MonkeyPatch) -> None:
    context = CancerContext()
    context.pathways = {"module": {"A", "C"}}
    monkeypatch.setattr(advanced, "classify_mechanism", lambda *args, **kwargs: _result("overlapping_exposure"))
    reports = advanced.counterfactual_analysis(
        [_target("A"), _target("B")], [_target("C"), _target("D")], context,
        _result(n_randomizations=10, n_valid_randomizations=10),
    )
    assert {item["kind"] for item in reports} == {"target", "module"}
    for report in reports:
        assert report["comparison_evaluable"] is False
        assert report["mechanism_changed"] is None
        assert report["descriptive_mechanism_changed"] is True
        assert report["score_delta"] == pytest.approx(-0.2)
        assert "baseline:insufficient_valid_randomizations" in report["comparison_audit"]["reason_codes"]
        assert "not_model_fidelity_or_biological_causality" in report["analysis_semantics"]


def _context() -> CancerContext:
    context = CancerContext()
    context.ppi.add_edge("A", "B1", details={"source": ["synthetic"]}, distance=1.0)
    context.ppi.add_edge("B1", "D", distance=1.0)
    context.ppi.add_edge("A", "B2", distance=1.0)
    context.ppi.add_edge("B2", "D", distance=1.0)
    context.ppi.graph["audit"] = {"items": ["baseline"]}
    context.expression = {"cell": {"A": 1.0}}
    context.pathways = {"module": {"A", "D"}}
    context.disease_sets = {"disease": {"D"}}
    context.default_disease_genes = {"D"}
    context.resource_metadata = {"ppi": {"nested": ["baseline"]}}
    context._disease_distance_cache[(0, "disease")] = {"A": 2.0}
    return context


def test_masked_snapshot_isolates_nested_graph_and_resource_objects() -> None:
    context = _context()
    masked = advanced._masked_path_context(context, ["B1"])
    assert "B1" not in masked.ppi and "B1" in context.ppi
    assert masked._disease_distance_cache == {}
    masked.ppi.graph["audit"]["items"].append("mutation")
    masked.expression["cell"]["A"] = 99.0
    masked.pathways["module"].add("MUTATION")
    masked.disease_sets["disease"].add("MUTATION")
    masked.default_disease_genes.add("MUTATION")
    masked.resource_metadata["ppi"]["nested"].append("mutation")
    assert context.ppi.graph["audit"]["items"] == ["baseline"]
    assert context.expression == {"cell": {"A": 1.0}}
    assert context.pathways == {"module": {"A", "D"}}
    assert context.disease_sets == {"disease": {"D"}}
    assert context.default_disease_genes == {"D"}
    assert context.resource_metadata == {"ppi": {"nested": ["baseline"]}}


def test_path_reports_share_actual_protocol_gate_and_each_mask_is_isolated(monkeypatch: pytest.MonkeyPatch) -> None:
    context = _context()
    paths = [{"path": ["A", bridge, "D"], "bridge_genes": [bridge]} for bridge in ("B1", "B2")]
    monkeypatch.setattr(advanced, "dynamic_hypotheses", lambda *args, **kwargs: [])
    monkeypatch.setattr(advanced, "local_evidence_paths", lambda *args, **kwargs: paths)
    monkeypatch.setattr(advanced, "build_evidence_graph", lambda *args, **kwargs: {"nodes": [], "edges": []})
    monkeypatch.setattr(advanced, "counterfactual_analysis", lambda *args, **kwargs: [])
    snapshots = []

    def classify(left: Any, right: Any, masked: CancerContext, **options: Any) -> MechanismResult:
        assert masked.expression == {"cell": {"A": 1.0}}
        assert masked.resource_metadata["ppi"]["nested"] == ["baseline"]
        snapshots.append(set(masked.ppi.nodes))
        masked.expression["cell"]["A"] = 99.0
        masked.resource_metadata["ppi"]["nested"].append("mutation")
        return _result("overlapping_exposure", null_B={"std": 0.0})

    monkeypatch.setattr(advanced, "classify_mechanism", classify)
    result = advanced.explain_advanced("drug-a", "drug-b", [_target("A")], [_target("D")], context, _result())
    assert snapshots == [{"A", "B2", "D"}, {"A", "B1", "D"}]
    assert context.expression == {"cell": {"A": 1.0}}
    assert context.resource_metadata["ppi"]["nested"] == ["baseline"]
    assert len(result.counterfactuals) == 2
    for report in result.counterfactuals:
        assert report["kind"] == "path"
        assert report["comparison_evaluable"] is False
        assert report["mechanism_changed"] is None
        assert report["descriptive_mechanism_changed"] is True
        assert report["score_delta"] == pytest.approx(-0.2)
        assert "perturbed:degenerate_or_missing_null_variance" in report["comparison_audit"]["reason_codes"]
    assert all(edge["type"] != "changes_mechanism" for edge in result.evidence_graph["edges"])
