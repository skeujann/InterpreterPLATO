"""Offline regression tests for the frozen-input statistical contract."""

from __future__ import annotations

from pathlib import Path

import networkx as nx
import pytest

import InterpreterPLATO.mechanism as mechanism_module
import InterpreterPLATO.interpreter as interpreter_module
from InterpreterPLATO.context import CancerContext, ContextScore
from InterpreterPLATO.mechanism import (
    NETWORK_METHOD_VERSION,
    MechanismResult,
    _cheng_metrics,
    classify_mechanism,
)
from InterpreterPLATO.scripts.test_support import archive_verified_plato_cache, strict_pipeline
from InterpreterPLATO.batch_statistics import apply_batch_fdr


def _score(symbol: str) -> ContextScore:
    return ContextScore(
        target_id=symbol, target_name=symbol, final_score=1.0,
        evidence_score=1.0, expression_score=None, disease_proximity_score=0.0,
        pathway_score=1.0, centrality_score=0.0, evidence_status="predicted_only",
        source="PLATO", source_evidence={},
    )


def _cycle(size: int = 80) -> nx.Graph:
    graph = nx.relabel_nodes(nx.cycle_graph(size), lambda value: f"G{value}")
    nx.set_edge_attributes(graph, 1.0, "distance")
    return graph


def test_custom_distance_explicitly_differs_from_disease_to_target_reference() -> None:
    graph = nx.path_graph([f"G{i}" for i in range(8)])
    nx.set_edge_attributes(graph, 1.0, "distance")
    result = _cheng_metrics(graph, {"G0", "G7"}, {"G1"}, {"G0"}, samples=0)
    assert result["d_A_obs"] == 3.5
    assert min(nx.shortest_path_length(graph, "G0", target) for target in ("G0", "G7")) == 0
    assert result["reference_method_reproduction"] is False
    assert result["null_exchangeability_validated"] is False
    assert result["method_version"] == NETWORK_METHOD_VERSION
    assert result["biological_interaction_tested"] is False


@pytest.mark.parametrize("invalid_draw,reason", [
    ({"G3", "ISOLATED"}, "incomplete_reachability"),
    ({"G3"}, "incomplete_size"),
])
def test_partial_null_draw_is_rejected_as_a_whole(monkeypatch, invalid_draw, reason) -> None:
    graph = _cycle()
    graph.add_node("ISOLATED")
    monkeypatch.setattr(mechanism_module, "_degree_matched_sample", lambda *args: (invalid_draw, 0.0))
    result = _cheng_metrics(graph, {"G0", "G1"}, {"G2", "G3"}, {"G0"}, samples=20)
    assert result["n_valid_randomizations"] == 0
    assert result["quality_status"] == "limited"
    assert result["statistically_supported"] is False
    for audit in result["randomization_audit"].values():
        assert audit["accepted"] == 0
        assert audit["rejected"] == 20
        assert audit["sample_size"] == 2
        assert audit["rejection_reasons"][reason] == 20


def test_disconnected_cycles_cannot_turn_truncated_samples_into_adequate_null() -> None:
    graph = nx.disjoint_union(nx.cycle_graph(100), nx.cycle_graph(100))
    graph = nx.relabel_nodes(graph, lambda value: f"G{value}")
    nx.set_edge_attributes(graph, 1.0, "distance")
    result = _cheng_metrics(
        graph, {f"G{i}" for i in range(10)}, {f"G{i}" for i in range(10, 20)},
        {"G0"}, samples=1000,
    )
    assert result["valid"] is True
    assert result["quality_status"] == "limited"
    assert result["n_valid_randomizations"] < 1000
    assert result["near_A"] is False and result["near_B"] is False
    for audit in result["randomization_audit"].values():
        assert audit["accepted"] + audit["rejected"] == 1000
        assert audit["rejection_reasons"]["incomplete_reachability"] > 900


@pytest.mark.parametrize("left,right,disease", [
    ({"G0", "ABSENT", "ABSENT2"}, {"G2"}, {"G0"}),
    ({"G0"}, {"G2", "ABSENT", "ABSENT2"}, {"G0"}),
    ({"G0"}, {"G2"}, {"G0", "ABSENT", "ABSENT2"}),
    ({"G0", "ISOLATED"}, {"G2"}, {"G0"}),
])
def test_insufficient_mapping_or_incomplete_reachability_fails_closed(left, right, disease) -> None:
    graph = _cycle()
    graph.add_node("ISOLATED")
    result = _cheng_metrics(graph, left, right, disease, samples=1000)
    assert result["valid"] is False
    assert result["statistically_supported"] is False
    assert result["both_drugs_disease_proximity_supported"] is False
    assert result["n_valid_randomizations"] == 0


def test_nonfinite_observed_distance_fails_closed() -> None:
    result = _cheng_metrics(_cycle(), {"G0"}, {"G1"}, {"G2"}, disease_distances={"G0": float("inf"), "G1": 1.0})
    assert result["valid"] is False


@pytest.mark.parametrize("budget", [0, 1, 999])
def test_insufficient_actual_null_cannot_claim_proximity_support(budget) -> None:
    result = _cheng_metrics(_cycle(), {"G0"}, {"G1"}, {"G0"}, samples=budget)
    assert result["quality_status"] == "limited"
    assert result["statistically_supported"] is False
    assert result["near_A"] is False and result["near_B"] is False


def test_degenerate_null_cannot_claim_proximity_support() -> None:
    graph = nx.complete_graph([f"G{i}" for i in range(12)])
    nx.set_edge_attributes(graph, 1.0, "distance")
    result = _cheng_metrics(graph, {"G0"}, {"G0"}, {"G0"}, samples=1000)
    assert result["null_A"]["std"] == 0.0
    assert result["quality_status"] == "limited"
    assert result["statistically_supported"] is False


def test_randomization_is_exactly_drug_swap_and_graph_insertion_invariant() -> None:
    graph = _cycle()
    left, right, disease = {"G0", "G1"}, {"G4", "G5", "G6"}, {"G0", "G10"}
    original = _cheng_metrics(graph, left, right, disease, samples=1000)
    swapped = _cheng_metrics(graph, right, left, disease, samples=1000)
    for first, second in (("z_A", "z_B"), ("p_A", "p_B"), ("null_A", "null_B"), ("d_A_obs", "d_B_obs")):
        assert original[first] == swapped[second]
        assert original[second] == swapped[first]
    assert original["s_AB"] == swapped["s_AB"]
    assert original["pattern"] == swapped["pattern"]
    assert original["randomization_audit"]["drug1"] == swapped["randomization_audit"]["drug2"]
    reordered = nx.Graph()
    reordered.add_nodes_from(reversed(list(graph.nodes)))
    reordered.add_edges_from(reversed(list(graph.edges(data=True))))
    assert original == _cheng_metrics(reordered, left, right, disease, samples=1000)


def test_identical_target_sets_use_identical_null_draws() -> None:
    result = _cheng_metrics(_cycle(), {"G0", "G1"}, {"G0", "G1"}, {"G0"}, samples=1000)
    assert result["null_A"] == result["null_B"]
    assert result["p_A"] == result["p_B"]
    assert result["randomization_audit"]["drug1"] == result["randomization_audit"]["drug2"]


def test_network_and_ora_inputs_do_not_depend_on_context_top_k_or_order() -> None:
    context = CancerContext()
    context.load_ppi(_cycle())
    context.default_disease_genes = {"G0"}
    context.load_pathways({
        "SHARED": {"G0", "G1", "G2", "G3", "G4", "G5"},
        **{f"NULL{i}": {f"N{i}_{j}" for j in range(5)} for i in range(99)},
    })
    left, right = [_score(f"G{i}") for i in range(3)], [_score(f"G{i}") for i in range(3, 6)]
    first = classify_mechanism(left, right, context, top_k=1, network_randomizations=20)
    second = classify_mechanism(list(reversed(left)), list(reversed(right)), context, top_k=3, network_randomizations=20)
    assert first.metrics["cheng"] == second.metrics["cheng"]
    assert first.metrics["cheng"]["statistical_target_sets"] == {"drug1": ["G0", "G1", "G2"], "drug2": ["G3", "G4", "G5"]}
    for key in ("significant_drug1", "significant_drug2", "significant_shared"):
        assert first.metrics["pathway_enrichment"][key] == second.metrics["pathway_enrichment"][key]


@pytest.mark.parametrize("stable_flag", [True, False])
def test_pipeline_preserves_full_statistical_targets_and_versions(tmp_path: Path, monkeypatch, stable_flag) -> None:
    cache = tmp_path / "cache"
    cache.mkdir()
    def records(symbols):
        return [{"target_symbol": symbol, "target_name": symbol, "organism": "Homo sapiens", "taxid": "9606", "target_type": "single protein", "score": 12.0, "reliable": "yes", "max_similarity": 0.8} for symbol in symbols]
    left, right = ["EGFR", "MTOR", "AKT1"], ["RRM1", "RRM2", "TYMS"]
    hashes = [archive_verified_plato_cache(cache, "CCO", "a", records(left)), archive_verified_plato_cache(cache, "CCN", "b", records(right))]
    pipeline = strict_pipeline(cache, hashes, top_k_targets=1)
    pipeline.load_demo_resources()
    monkeypatch.setattr(pipeline, "_top_k_stability", lambda *args, **kwargs: {"evaluable": True, "descriptive_stable": True, "stable": stable_flag, "agreement_fraction": 1.0, "meets_full_robustness_criteria": False})
    result = pipeline.explain_combination("a", "b", "CCO", "CCN", advanced=False)
    payload = result.to_dict()
    assert len(payload["top_targets_drug1"]) == 1
    assert payload["network_statistics"]["statistical_target_sets"] == {"drug1": sorted(left), "drug2": sorted(right)}
    assert payload["result_schema_version"] == "5.4"
    assert payload["statistical_contract_version"] == "1.1"
    assert payload["audit"]["policies"]["network_statistics"] == NETWORK_METHOD_VERSION
    assert payload["audit"]["policies"]["statistical_contract"] == "1.1"
    assert payload["credibility"]["components"]["robustness"]["points"] == 0.0
    assert "mechanism_sensitive_to_top_k" not in payload["credibility"]["warnings"]
    assert payload["generated_from_test_fixture"] is True
    assert payload["formal_result_eligible"] is False


@pytest.mark.parametrize("pattern,kind,near_a,near_b,expected", [
    ("P4", "single_exposure", True, False, True),
    ("P3", "overlapping_single_exposure", False, True, True),
    ("P1", "high_overlap", True, True, True),
    ("P2", "complementary_exposure", True, True, True),
    ("P4", "single_exposure", True, True, False),
    ("P2", "complementary_exposure", True, False, False),
    ("P4", "single_exposure", False, False, False),
])
def test_primary_network_support_matches_pattern_without_duplicate_convergence_bonus(
    tmp_path, monkeypatch, pattern, kind, near_a, near_b, expected,
) -> None:
    cache = tmp_path / "cache"
    cache.mkdir()
    hashes = [archive_verified_plato_cache(cache, smiles, name) for smiles, name in (("CCO", "a"), ("CCN", "b"))]
    pipeline = strict_pipeline(cache, hashes)
    pipeline.load_demo_resources()
    network = {
        "valid": True, "quality_status": "adequate", "pattern": pattern,
        "near_A": near_a, "near_B": near_b,
        "p_A": 0.001 if near_a else 0.8, "p_B": 0.001 if near_b else 0.8,
        "statistically_supported": near_a or near_b,
        "pair_interaction_supported": near_a and near_b,
        "both_drugs_disease_proximity_supported": near_a and near_b,
        "n_randomizations": 1000, "n_valid_randomizations": 1000,
        "null_A": {"std": 1.0}, "null_B": {"std": 1.0},
    }
    baseline = MechanismResult(kind, 0.6, "Synthetic pattern fixture", metrics={"cheng": network})
    monkeypatch.setattr(interpreter_module, "classify_mechanism", lambda *args, **kwargs: baseline)
    result = pipeline.explain_combination("a", "b", "CCO", "CCN", advanced=False)
    assert result.credibility["mechanism_supported"] is expected
    component = result.credibility["components"]["mechanism_support"]
    assert component["inputs"]["convergent_mechanism_evidence"] is False
    assert component["value"] == pytest.approx(0.45 * 0.6 + 0.35 * expected)
    apply_batch_fdr([result])
    assert result.credibility["single_pair_mechanism_supported"] is expected
    assert result.credibility["batch_mechanism_supported"] is expected
    assert result.to_dict()["network_both_drugs_disease_proximity_supported"] is (near_a and near_b)
