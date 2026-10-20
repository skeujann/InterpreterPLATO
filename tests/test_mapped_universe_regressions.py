"""Software controls for mapped-universe inference and annotation-count invariance.

These deliberately synthetic controls are not drug mechanism validation cases.
"""
from __future__ import annotations

import networkx as nx
import pytest

import InterpreterPLATO.mechanism as mechanism_module
from InterpreterPLATO.context import CancerContext, ContextScore
from InterpreterPLATO.mechanism import _cheng_metrics, classify_mechanism


def cycle(size=400):
    graph = nx.relabel_nodes(nx.cycle_graph(size), lambda n: f"G{n}")
    nx.set_edge_attributes(graph, 1.0, "distance")
    return graph


def score(symbol):
    return ContextScore(
        target_id=symbol, target_name=symbol, final_score=1.0,
        evidence_score=1.0, expression_score=None, disease_proximity_score=0.0,
        pathway_score=1.0, centrality_score=0.0, evidence_status="predicted_only",
        source="PLATO", source_evidence={},
    )


@pytest.mark.parametrize("side", ["drug1", "drug2", "disease_module"])
def test_partial_mapping_matches_explicit_intersection_and_reports_exclusions(side):
    graph = cycle()
    sets = {"drug1": {"G1", "G2"}, "drug2": {"G3", "G4"}, "disease_module": {"G0"}}
    reference = _cheng_metrics(graph, *sets.values(), samples=1000)
    sets[side].add("UNMAPPED")
    result = _cheng_metrics(graph, *sets.values(), samples=1000)
    assert result["valid"] is True
    assert result["quality_status"] == "adequate"
    for key in ("d_A_obs", "d_B_obs", "d_AB", "s_AB", "p_A", "p_B", "z_A", "z_B", "null_A", "null_B", "randomization_audit"):
        assert result[key] == reference[key]
    assert result["excluded_unmapped_genes"][side] == ["UNMAPPED"]
    assert result["network_coverage"][side]["input"] == len(sets[side])
    assert result["network_coverage"][side]["mapping_fraction"] == (len(sets[side])-1)/len(sets[side])
    assert result["analysis_universe"] == reference["analysis_universe"]
    assert "network_analysis_conditional_on_mapped_genes" in result["warnings"]
    assert "UNMAPPED" in sets[side]  # No in-place deletion of caller input.


def test_strict_coverage_is_still_available_and_empty_intersection_cannot_pass():
    graph = cycle()
    strict = _cheng_metrics(graph, {"G1", "UNMAPPED"}, {"G2"}, {"G0"}, min_mapping_coverage=1.0)
    empty = _cheng_metrics(graph, {"UNMAPPED"}, {"G2"}, {"G0"}, min_mapping_coverage=0.0)
    for result in (strict, empty):
        assert result["valid"] is False
        assert result["n_valid_randomizations"] == 0
        assert result["statistically_supported"] is False
    assert "low_network_mapping_coverage" in strict["warnings"]


@pytest.mark.parametrize("left,right,near_a,near_b,pattern", [
    ({"G1", "G2"}, {"G3", "G4"}, True, True, "P2"),
    ({"G1", "G2"}, {"G200", "G201"}, True, False, "P4"),
    ({"G199", "G200"}, {"G202", "G203"}, False, False, "P6"),
])
def test_support_discordant_and_negative_network_controls(left, right, near_a, near_b, pattern):
    result = _cheng_metrics(cycle(), left, right, {"G0"}, samples=1000, random_seed=42)
    assert result["quality_status"] == "adequate"
    assert result["n_valid_randomizations"] == 1000
    assert result["near_A"] is near_a and result["near_B"] is near_b
    assert result["pattern"] == pattern
    assert result["both_drugs_disease_proximity_supported"] is (near_a and near_b)
    assert result["biological_interaction_tested"] is False


def annotation_context(monkeypatch, count, shared=False):
    context = CancerContext()
    context.load_ppi(cycle())
    context.default_disease_genes = {"G0"}
    context.load_pathways({f"S{i}": {"G1", "G2", "G3", "G4", "G5"} for i in range(count)})
    def enrichment(targets):
        prefix = "S" if shared else ("A" if "G1" in targets else "B")
        return [{"pathway": f"{prefix}{i}", "q": 0.01} for i in range(count)]
    # Freeze q-values: this tests the score rule, separately from ORA testing.
    monkeypatch.setattr(context, "pathway_enrichment", enrichment)
    return context


@pytest.mark.parametrize("shared,qualified,far", [
    (False, False, False), (False, False, True), (True, False, False), (True, True, False),
])
def test_more_pathway_annotations_do_not_inflate_rule_score(monkeypatch, shared, qualified, far):
    if qualified:
        monkeypatch.setattr(mechanism_module, "qualified_target_action_status", lambda item: "inhibitory")
    results = []
    for count in (1, 5):
        context = annotation_context(monkeypatch, count, shared)
        result = classify_mechanism(
            [score("G1"), score("G2")],
            [score("G200"), score("G201")] if far else [score("G3"), score("G4")],
            context, network_randomizations=0,
        )
        results.append(result)
    assert results[0].mechanism_type == results[1].mechanism_type
    assert results[0].confidence == results[1].confidence
    assert results[0].score_semantics == "uncalibrated_evidence_strength"
    if not shared:
        assert len(results[0].complementary_pathways) == 1
        assert len(results[1].complementary_pathways) == 5
        expected = 0.35 if far else 0.525  # Bidirectional nearest mean is 1.5.
        assert results[0].confidence == pytest.approx(expected)
    else:
        assert results[0].confidence == pytest.approx(0.60 if qualified else 0.52)


def test_shared_target_count_does_not_add_bonus_beyond_jaccard(monkeypatch):
    context = CancerContext()
    context.load_ppi(cycle())
    context.default_disease_genes = {"G0"}
    monkeypatch.setattr(mechanism_module, "qualified_target_action_status", lambda item: "inhibitory")
    results = []
    for symbols in (("G1",), ("G1", "G2", "G3")):
        profile = [score(symbol) for symbol in symbols]
        results.append(classify_mechanism(profile, profile, context, network_randomizations=0))
    assert all(result.mechanism_type == "target_overlap_or_redundancy" for result in results)
    assert results[0].confidence == results[1].confidence == pytest.approx(0.90)
