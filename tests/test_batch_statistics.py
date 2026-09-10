"""Batch FDR must respect mechanism-specific requirements and planned tests."""

from copy import deepcopy
from types import SimpleNamespace

import pytest

from InterpreterPLATO.batch_statistics import apply_batch_fdr
from InterpreterPLATO.mechanism import MechanismResult


def _result(mechanism_type="single_exposure", p_a=0.001, p_b=0.8, *, name="a", quality="adequate", local=True):
    both = mechanism_type in {"high_overlap", "complementary_exposure"}
    return SimpleNamespace(
        drug1_name=name,
        drug2_name="b",
        cell_line="A549",
        mechanism=MechanismResult(
            mechanism_type=mechanism_type,
            confidence=0.6,
            description="Synthetic statistical fixture",
            credibility_score=75.0,
            metrics={"cheng": {
                "valid": True,
                "quality_status": quality,
                "n_randomizations": 1000,
                "n_valid_randomizations": 1000,
                "null_A": {"std": 1.0},
                "null_B": {"std": 1.2},
                "pattern": "P2" if both else "P4",
                "near_A": True,
                "near_B": both,
                "z_A": -3.0,
                "z_B": -2.0 if both else 0.1,
                "p_A": p_a,
                "p_B": p_b,
                "p_A_adjusted": min(1.0, 2.0 * p_a) if isinstance(p_a, (int, float)) else p_a,
                "p_B_adjusted": min(1.0, 2.0 * p_b) if isinstance(p_b, (int, float)) else p_b,
                "statistically_supported": True,
                "pair_interaction_supported": both,
            }},
        ),
        credibility={
            "final_score": 75.0,
            "evidence_credibility_score": 75.0,
            "mechanism_supported": local,
            "score_profile": {"name": "uncalibrated_local_rules", "formula_version": "fixture"},
        },
    )


@pytest.mark.parametrize("mechanism_type", ["single_exposure", "overlapping_single_exposure", "indirect_exposure"])
def test_single_exposure_only_requires_originally_exposed_side(mechanism_type):
    result = _result(mechanism_type)
    local = deepcopy(result.mechanism.metrics["cheng"])
    apply_batch_fdr([result])
    cheng = result.mechanism.metrics["cheng"]
    assert cheng["batch_q_A"] == pytest.approx(0.002)
    assert cheng["batch_q_B"] == pytest.approx(0.8)
    assert cheng["batch_required_sides"] == ["A"]
    assert cheng["batch_pattern_supported"] is True
    assert result.credibility["batch_mechanism_supported"] is True
    assert result.credibility["mechanism_supported"] is True
    assert {key: cheng[key] for key in local} == local
    assert result.credibility["final_score"] == result.credibility["batch_adjusted_score"] == 75.0
    assert result.mechanism.credibility_score == 75.0
    assert result.credibility["score_profile"]["batch_FDR_adjustment_applied"] is False
    assert result.credibility["score_profile"]["name"] == "uncalibrated_local_rules"


def test_two_drug_pattern_failure_clears_canonical_but_preserves_local_support():
    result = _result("complementary_exposure")
    apply_batch_fdr([result])
    assert result.mechanism.metrics["cheng"]["batch_required_sides"] == ["A", "B"]
    assert result.mechanism.metrics["cheng"]["pair_interaction_supported"] is True
    assert result.mechanism.metrics["cheng"]["batch_pattern_supported"] is False
    assert result.credibility["single_pair_mechanism_supported"] is True
    assert result.credibility["mechanism_supported"] is False
    assert result.credibility["batch_mechanism_supported"] is False
    assert result.credibility["final_score"] == 75.0
    assert "network_support_not_significant_after_batch_FDR" in result.credibility["warnings"]
    assert result.mechanism.mechanism_type == "complementary_exposure"


def test_inadequate_null_is_retained_in_family_but_can_never_support():
    adequate = _result(name="adequate")
    inadequate = _result(name="inadequate", quality="limited", p_a=0.0, p_b=0.0)
    apply_batch_fdr([adequate, inadequate])
    cheng = inadequate.mechanism.metrics["cheng"]
    assert cheng["batch_hypothesis_count"] == 4
    assert cheng["batch_family"]["evaluable_hypothesis_count"] == 2
    assert cheng["batch_family"]["status"] == "partially_evaluable"
    assert cheng["batch_q_A"] == cheng["batch_p_A_used"] == 1.0
    assert cheng["batch_test_A_evaluable"] is False
    assert cheng["batch_pattern_status"] == "unevaluable"
    assert inadequate.credibility["mechanism_supported"] is False
    assert adequate.mechanism.metrics["cheng"]["batch_q_A"] == pytest.approx(0.004)


@pytest.mark.parametrize("field", ["n_randomizations", "n_valid_randomizations"])
@pytest.mark.parametrize("invalid_count", [None, True, "1000", 1000.0, 1000.5, float("nan"), float("inf"), -1, 999])
def test_adequate_label_cannot_override_invalid_actual_draw_counts(field, invalid_count):
    result = _result(p_a=0.0, p_b=0.0)
    cheng = result.mechanism.metrics["cheng"]
    cheng[field] = invalid_count
    apply_batch_fdr([result])
    assert cheng["batch_family"]["status"] == "unevaluable"
    assert cheng["batch_p_A_used"] == cheng["batch_p_B_used"] == 1.0
    assert cheng["batch_test_A_status"] == "invalid_or_inadequate_network_null"
    assert result.credibility["mechanism_supported"] is False


@pytest.mark.parametrize("requested,accepted", [(1000, 1001), (2000, 1799), (10000, 1000)])
def test_adequate_label_requires_consistent_counts_and_acceptance_fraction(requested, accepted):
    result = _result()
    cheng = result.mechanism.metrics["cheng"]
    cheng.update({"n_randomizations": requested, "n_valid_randomizations": accepted})
    apply_batch_fdr([result])
    assert cheng["batch_test_A_evaluable"] is False
    assert cheng["batch_test_B_evaluable"] is False
    assert cheng["batch_pattern_supported"] is False


def test_exact_ninety_percent_acceptance_with_enough_draws_is_evaluable():
    result = _result()
    cheng = result.mechanism.metrics["cheng"]
    cheng.update({"n_randomizations": 2000, "n_valid_randomizations": 1800})
    apply_batch_fdr([result])
    assert cheng["batch_test_A_evaluable"] is True
    assert cheng["batch_pattern_supported"] is True


@pytest.mark.parametrize("side", ["A", "B"])
@pytest.mark.parametrize("invalid_std", [None, True, "bad", float("nan"), float("inf"), float("-inf"), 0.0, -1.0])
def test_adequate_label_cannot_override_invalid_null_variation(side, invalid_std):
    result = _result(p_a=0.0, p_b=0.0)
    cheng = result.mechanism.metrics["cheng"]
    cheng[f"null_{side}"]["std"] = invalid_std
    apply_batch_fdr([result])
    assert cheng["batch_test_A_evaluable"] is False
    assert cheng["batch_test_B_evaluable"] is False
    assert cheng["batch_q_A"] == cheng["batch_q_B"] == 1.0


@pytest.mark.parametrize("missing_field", ["n_randomizations", "n_valid_randomizations", "null_A", "null_B"])
def test_legacy_adequate_label_without_null_audit_is_unevaluable(missing_field):
    result = _result()
    cheng = result.mechanism.metrics["cheng"]
    del cheng[missing_field]
    apply_batch_fdr([result])
    assert cheng["batch_family"]["status"] == "unevaluable"
    assert result.credibility["mechanism_supported"] is False


@pytest.mark.parametrize("malformed_null", [None, [], "adequate", {"mean": 1.0}])
def test_malformed_null_report_fails_closed(malformed_null):
    result = _result()
    cheng = result.mechanism.metrics["cheng"]
    cheng["null_A"] = malformed_null
    apply_batch_fdr([result])
    assert cheng["batch_family"]["status"] == "unevaluable"


@pytest.mark.parametrize("bad_p", [None, float("nan"), float("inf"), -0.1, 1.1, "bad", True])
def test_invalid_probabilities_are_planned_unevaluable_tests(bad_p):
    result = _result(p_a=bad_p, p_b=bad_p)
    apply_batch_fdr([result])
    cheng = result.mechanism.metrics["cheng"]
    assert cheng["batch_family"]["status"] == "unevaluable"
    assert cheng["batch_hypothesis_count"] == 2
    assert cheng["batch_q_A"] == cheng["batch_q_B"] == 1.0
    assert cheng["batch_pattern_supported"] is False


def test_family_identity_and_q_values_are_order_invariant():
    first = [_result(name="alpha"), _result(name="beta", p_a=0.01)]
    second = deepcopy(list(reversed(first)))
    apply_batch_fdr(first)
    apply_batch_fdr(second)
    by_name = {item.drug1_name: item.mechanism.metrics["cheng"] for item in first}
    for result in second:
        actual = result.mechanism.metrics["cheng"]
        expected = by_name[result.drug1_name]
        assert actual["batch_family"] == expected["batch_family"]
        assert actual["batch_q_A"] == expected["batch_q_A"]
    duplicate = deepcopy(first)
    duplicate.append(deepcopy(first[0]))
    apply_batch_fdr(duplicate)
    assert duplicate[0].mechanism.metrics["cheng"]["batch_family_id"] != by_name["alpha"]["batch_family_id"]


def test_batch_never_upgrades_previously_unsupported_mechanism():
    result = _result(local=False)
    apply_batch_fdr([result])
    assert result.mechanism.metrics["cheng"]["batch_pattern_supported"] is True
    assert result.credibility["single_pair_mechanism_supported"] is False
    assert result.credibility["batch_mechanism_supported"] is False
    assert result.credibility["mechanism_supported"] is False


def test_non_network_primary_support_remains_local():
    result = _result("vertical_pathway_coblockade")
    apply_batch_fdr([result])
    assert result.credibility["batch_mechanism_supported"] is None
    assert result.credibility["mechanism_supported"] is True
    assert "warnings" not in result.credibility


def test_legacy_local_near_fallback_identifies_side_b():
    result = _result(p_a=0.8, p_b=0.001)
    cheng = result.mechanism.metrics["cheng"]
    del cheng["near_A"]
    del cheng["near_B"]
    cheng["z_A"], cheng["z_B"] = 0.1, -3.0
    apply_batch_fdr([result])
    assert cheng["batch_required_sides"] == ["B"]
    assert cheng["batch_pattern_supported"] is True


def test_explicit_local_near_failure_cannot_be_overridden_by_small_p():
    result = _result()
    result.mechanism.metrics["cheng"]["near_A"] = False
    apply_batch_fdr([result])
    assert result.mechanism.metrics["cheng"]["batch_pattern_supported"] is False


def test_missing_network_evidence_is_included_as_unevaluable():
    result = _result()
    result.mechanism.metrics.clear()
    apply_batch_fdr([result])
    assert result.mechanism.metrics["cheng"]["batch_family"]["status"] == "unevaluable"
    assert result.credibility["mechanism_supported"] is False


def test_repeated_application_preserves_single_pair_snapshot_and_deduplicates_warnings():
    result = _result("complementary_exposure")
    apply_batch_fdr([result])
    first = deepcopy(result.credibility)
    apply_batch_fdr([result])
    assert result.credibility == first


def test_empty_batch_is_safe():
    results = []
    assert apply_batch_fdr(results) is results


def _identified_result():
    result = _result()
    result.metadata = {
        "disease_context_labels": ["NSCLC"],
        "audit": {
            "plato_inputs": {
                "drug1": {"canonical_isomeric_smiles": "CCO", "raw_response_sha256": "a" * 64},
                "drug2": {"canonical_isomeric_smiles": "CCN", "raw_response_sha256": "b" * 64},
            },
            "runtime": {"artifact_inventory": {
                "ppi": {"sha256": "1" * 64},
                "expression": {"sha256": "2" * 64},
                "disease": {"sha256": "3" * 64},
                "pathways": {"sha256": ["4" * 64]},
            }},
        },
    }
    result.mechanism.metrics["cheng"].update({
        "statistical_target_sets": {"drug1": ["EGFR", "ERBB2"], "drug2": ["MTOR"]},
        "method_version": "frozen_method_v2",
    })
    return result


def _family(result):
    apply_batch_fdr([result])
    return result.mechanism.metrics["cheng"]["batch_family"]


def test_names_only_identity_is_explicitly_incomplete():
    family = _family(_result())
    assert family["identity_complete"] is False
    member = family["members"][0]
    assert member["method_version"] is None
    assert member["disease_context_labels"] is None
    assert member["drugs"][0]["statistical_targets"] is None
    assert "method_version" in member["missing_identity_fields"]


def test_different_disease_same_drug_names_has_distinct_family_identity():
    first = _identified_result()
    second = deepcopy(first)
    second.metadata["disease_context_labels"] = ["breast_cancer"]
    assert _family(first)["identity_complete"] is True
    assert _family(first)["id"] != _family(second)["id"]


def test_different_frozen_targets_same_drug_names_has_distinct_family_identity():
    first = _identified_result()
    second = deepcopy(first)
    second.mechanism.metrics["cheng"]["statistical_target_sets"]["drug1"] = ["KRAS"]
    assert _family(first)["id"] != _family(second)["id"]


@pytest.mark.parametrize("changed_field", ["structure", "source", "artifact", "method"])
def test_family_identity_tracks_input_and_method_changes(changed_field):
    first = _identified_result()
    second = deepcopy(first)
    if changed_field == "structure":
        second.metadata["audit"]["plato_inputs"]["drug1"]["canonical_isomeric_smiles"] = "CCCCO"
    elif changed_field == "source":
        second.metadata["audit"]["plato_inputs"]["drug1"]["raw_response_sha256"] = "c" * 64
    elif changed_field == "artifact":
        second.metadata["audit"]["runtime"]["artifact_inventory"]["ppi"]["sha256"] = "5" * 64
    else:
        second.mechanism.metrics["cheng"]["method_version"] = "new_method"
    assert _family(first)["id"] != _family(second)["id"]


def test_family_identity_excludes_outcomes_and_machine_local_inventory_fields():
    first = _identified_result()
    second = deepcopy(first)
    second.mechanism.metrics["cheng"].update({"p_A": 0.9, "z_A": 0.1, "near_A": False})
    second.metadata["audit"]["generated_at_utc"] = "different time"
    second.metadata["audit"]["runtime"]["artifact_inventory"]["PLATO_cache"] = {"root": "different/path"}
    second.metadata["audit"]["runtime"]["artifact_inventory"]["ppi"]["path"] = "different/path"
    second.mechanism.metrics["cheng"]["statistical_target_sets"]["drug1"].reverse()
    assert _family(first)["id"] == _family(second)["id"]


@pytest.mark.parametrize("missing_field", ["method", "disease", "targets"])
def test_missing_inference_identity_never_claims_complete(missing_field):
    result = _identified_result()
    if missing_field == "method":
        result.mechanism.metrics["cheng"].pop("method_version")
    elif missing_field == "disease":
        result.metadata.pop("disease_context_labels")
    else:
        result.mechanism.metrics["cheng"].pop("statistical_target_sets")
    assert _family(result)["identity_complete"] is False
