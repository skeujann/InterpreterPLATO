"""Batch-level multiplicity control without rewriting local evidence or scores."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from hashlib import sha256
import json
import math
from numbers import Integral
from typing import Any, Optional


BATCH_METHOD = "Benjamini-Hochberg_across_all_planned_pair_by_drug_tests"
BATCH_SCOPE = "all_results_in_this_declared_batch_not_automatically_the_entire_study"
NETWORK_TYPES = {
    "high_overlap", "complementary_exposure", "single_exposure",
    "overlapping_single_exposure", "indirect_exposure",
}
PATTERN_TYPES = {
    "P1": "high_overlap",
    "P2": "complementary_exposure",
    "P3": "overlapping_single_exposure",
    "P4": "single_exposure",
    "P4_edge": "overlapping_single_exposure",
}


def _finite_number(value: Any) -> Optional[float]:
    if value is None or isinstance(value, bool):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError, OverflowError):
        return None
    return number if math.isfinite(number) else None


def _probability(value: Any) -> Optional[float]:
    number = _finite_number(value)
    return number if number is not None and 0.0 <= number <= 1.0 else None


def _adequate_null(cheng: Mapping[str, Any]) -> bool:
    if cheng.get("valid") is not True or cheng.get("quality_status") != "adequate":
        return False
    requested = cheng.get("n_randomizations")
    accepted = cheng.get("n_valid_randomizations")
    if any(isinstance(value, bool) or not isinstance(value, Integral) for value in (requested, accepted)):
        return False
    requested, accepted = int(requested), int(accepted)
    if requested < 1000 or accepted < 1000 or accepted > requested or accepted * 10 < requested * 9:
        return False
    for side in ("A", "B"):
        null = cheng.get(f"null_{side}")
        if not isinstance(null, Mapping):
            return False
        std = _finite_number(null.get("std"))
        if std is None or std <= 0.0:
            return False
    return True


def _local_near(cheng: Mapping[str, Any], side: str, alpha: float) -> bool:
    if f"near_{side}" in cheng:
        return cheng[f"near_{side}"] is True
    z_score = _finite_number(cheng.get(f"z_{side}"))
    adjusted_p = _probability(cheng.get(f"p_{side}_adjusted"))
    return bool(z_score is not None and z_score < 0.0 and adjusted_p is not None and adjusted_p <= alpha)


def _required_sides(pattern_type: str, cheng: Mapping[str, Any], alpha: float) -> list[str]:
    if pattern_type in {"high_overlap", "complementary_exposure"}:
        return ["A", "B"]
    if pattern_type in {"single_exposure", "overlapping_single_exposure", "indirect_exposure"}:
        near = [side for side in ("A", "B") if _local_near(cheng, side, alpha)]
        return near if len(near) == 1 else []
    return []


def _append_unique(target: dict[str, Any], key: str, value: str) -> None:
    target[key] = list(dict.fromkeys([*target.get(key, []), value]))


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _identity_text(value: Any) -> Optional[str]:
    if value is None:
        return None
    return str(value).strip() or None


def _identity_list(value: Any) -> Optional[list[str]]:
    if not isinstance(value, (list, tuple, set, frozenset)):
        return None
    return sorted({text for item in value if (text := _identity_text(item)) is not None})


def _artifact_identity(inventory: Mapping[str, Any]) -> dict[str, Any]:
    identity = {}
    for resource_name in ("ppi", "expression", "disease", "pathways"):
        resource = _mapping(inventory.get(resource_name))
        raw_hash = resource.get("sha256")
        identity[resource_name] = {
            "sha256": _identity_list(raw_hash) if isinstance(raw_hash, (list, tuple, set, frozenset)) else _identity_text(raw_hash),
            "artifact_kind": _identity_text(resource.get("artifact_kind")),
            "artifact_kinds": _identity_list(resource.get("artifact_kinds")),
            "version": _identity_text(resource.get("version")),
            "mapping_version": _identity_text(resource.get("mapping_version")),
            "cell_line": _identity_text(resource.get("cell_line")),
            "label": _identity_text(resource.get("label")),
        }
    return identity


def _family_member(result: Any) -> dict[str, Any]:
    metadata = _mapping(getattr(result, "metadata", None))
    audit = _mapping(metadata.get("audit"))
    inputs = _mapping(audit.get("plato_inputs"))
    cheng = _mapping(result.mechanism.metrics.get("cheng"))
    target_sets = _mapping(cheng.get("statistical_target_sets"))
    runtime = _mapping(audit.get("runtime"))
    artifacts = _artifact_identity(_mapping(runtime.get("artifact_inventory")))
    drugs = []
    missing = []
    for side in ("drug1", "drug2"):
        source = _mapping(inputs.get(side))
        drug = {
            "name": _identity_text(getattr(result, f"{side}_name", None)),
            "canonical_isomeric_smiles": _identity_text(source.get("canonical_isomeric_smiles")),
            "inchikey": _identity_text(source.get("inchikey")),
            "raw_response_sha256": _identity_text(source.get("raw_response_sha256")),
            "payload_sha256": _identity_text(source.get("payload_sha256")),
            "normalized_targets_sha256": _identity_text(source.get("normalized_targets_sha256")),
            "provider": _identity_text(source.get("provider")),
            "release": _identity_text(source.get("release")),
            "parser_version": _identity_text(source.get("parser_version")),
            "schema_fingerprint": _identity_text(source.get("schema_fingerprint")),
            "statistical_targets": _identity_list(target_sets.get(side)),
        }
        drugs.append(drug)
    drugs.sort(key=lambda item: json.dumps(item, sort_keys=True, ensure_ascii=False))
    for index, drug in enumerate(drugs):
        for field in ("name", "canonical_isomeric_smiles", "raw_response_sha256", "statistical_targets"):
            if not drug[field]:
                missing.append(f"drugs[{index}].{field}")
    cell_line = _identity_text(getattr(result, "cell_line", None))
    disease_labels = _identity_list(metadata.get("disease_context_labels"))
    method_version = _identity_text(cheng.get("method_version"))
    for field, value in (("cell_line", cell_line), ("disease_context_labels", disease_labels), ("method_version", method_version)):
        if not value:
            missing.append(field)
    for resource_name, identity in artifacts.items():
        if not identity["sha256"]:
            missing.append(f"runtime_artifacts.{resource_name}.sha256")
    return {
        "drugs": drugs,
        "cell_line": cell_line,
        "disease_context_labels": disease_labels,
        "method_version": method_version,
        "runtime_artifacts": artifacts,
        "missing_identity_fields": missing,
    }


def apply_batch_fdr(results: Sequence[Any], *, alpha: float = 0.05) -> Sequence[Any]:
    """Apply BH to every planned pair-by-drug test in an explicitly scoped batch.

    Unavailable tests remain in the family at p = 1 but cannot establish support.
    Numeric credibility scores remain unchanged: batch support is a separate
    inferential state, not an arbitrary numerical confidence penalty.
    """
    threshold = _probability(alpha)
    if threshold is None or not 0.0 < threshold < 1.0:
        raise ValueError("alpha must be a finite probability strictly between zero and one")
    if not results:
        return results

    members = [_family_member(result) for result in results]
    members.sort(key=lambda item: json.dumps(item, sort_keys=True, ensure_ascii=False))
    family_identity = {
        "version": "planned_pair_by_drug_family_v2",
        "method": BATCH_METHOD,
        "scope": BATCH_SCOPE,
        "alpha": threshold,
        "members": members,
        "tests_per_pair": 2,
    }
    family_id = sha256(json.dumps(
        family_identity, sort_keys=True, ensure_ascii=False, separators=(",", ":"),
    ).encode("utf-8")).hexdigest()

    hypotheses = []
    for result in results:
        cheng = result.mechanism.metrics.get("cheng", {})
        if not isinstance(cheng, dict):
            cheng = {}
            result.mechanism.metrics["cheng"] = cheng
        elif "cheng" not in result.mechanism.metrics:
            result.mechanism.metrics["cheng"] = cheng
        quality_ok = _adequate_null(cheng)
        for side in ("A", "B"):
            p_value = _probability(cheng.get(f"p_{side}"))
            evaluable = quality_ok and p_value is not None
            reason = (
                "evaluable" if evaluable else
                "invalid_or_inadequate_network_null" if not quality_ok else
                "missing_or_invalid_p_value"
            )
            hypotheses.append({
                "cheng": cheng,
                "side": side,
                "p": p_value if evaluable else 1.0,
                "evaluable": evaluable,
                "reason": reason,
            })

    total = len(hypotheses)
    order = sorted(range(total), key=lambda index: hypotheses[index]["p"])
    q_values = [1.0] * total
    running = 1.0
    for rank_index in range(total - 1, -1, -1):
        original_index = order[rank_index]
        running = min(running, hypotheses[original_index]["p"] * total / (rank_index + 1))
        q_values[original_index] = min(1.0, running)

    evaluable_count = sum(item["evaluable"] for item in hypotheses)
    family_status = (
        "evaluable" if evaluable_count == total else
        "partially_evaluable" if evaluable_count else "unevaluable"
    )
    family = {
        **family_identity,
        "id": family_id,
        "pair_count": len(results),
        "hypothesis_count": total,
        "evaluable_hypothesis_count": evaluable_count,
        "unevaluable_hypothesis_count": total - evaluable_count,
        "status": family_status,
        "missing_test_policy": "retain_planned_hypothesis_with_p_1_and_never_treat_as_support",
        "identity_complete": all(not member["missing_identity_fields"] for member in members),
        "identity_semantics": "drug_structures_source_exports_frozen_targets_context_artifacts_and_method_not_p_values_or_observed_support",
        "study_protocol_requirement": "freeze_the_complete_study_family_before_inspecting_results_do_not_split_batches_to_select_significance",
        "dependence_assumption": "BH_requires_independence_or_suitable_positive_dependence_not_verified_by_this_function",
    }
    for index, hypothesis in enumerate(hypotheses):
        cheng, side = hypothesis["cheng"], hypothesis["side"]
        cheng[f"batch_q_{side}"] = q_values[index]
        cheng[f"batch_p_{side}_used"] = hypothesis["p"]
        cheng[f"batch_test_{side}_evaluable"] = hypothesis["evaluable"]
        cheng[f"batch_test_{side}_status"] = hypothesis["reason"]
        cheng[f"batch_near_{side}"] = bool(
            hypothesis["evaluable"] and q_values[index] <= threshold
            and _local_near(cheng, side, threshold)
        )

    for result in results:
        cheng = result.mechanism.metrics["cheng"]
        primary_type = result.mechanism.mechanism_type
        network_primary = primary_type in NETWORK_TYPES
        pattern_type = primary_type if network_primary else PATTERN_TYPES.get(cheng.get("pattern"), "")
        required = _required_sides(pattern_type, cheng, threshold)
        supported = bool(required) and all(cheng[f"batch_near_{side}"] for side in required)
        pattern_status = (
            "not_applicable" if pattern_type not in NETWORK_TYPES else
            "unevaluable" if not required or not all(cheng[f"batch_test_{side}_evaluable"] for side in required) else
            "supported" if supported else "not_supported"
        )
        cheng.update({
            "batch_multiple_testing_method": BATCH_METHOD,
            "batch_hypothesis_count": total,
            "batch_family_id": family_id,
            "batch_family": dict(family),
            "batch_required_sides": required,
            "batch_pattern_supported": supported,
            "batch_pattern_status": pattern_status,
            "batch_support_semantics": "local_pattern_requirements_and_evaluable_BH_q_values_no_new_biological_claim_from_nonsignificance",
        })

        credibility = result.credibility
        credibility.setdefault("single_pair_mechanism_supported", bool(credibility.get("mechanism_supported", False)))
        local_supported = credibility["single_pair_mechanism_supported"]
        credibility["batch_mechanism_supported"] = bool(local_supported and supported) if network_primary else None
        if network_primary:
            credibility["mechanism_supported"] = bool(local_supported and supported)
        credibility.setdefault("base_score", credibility.get("final_score", credibility.get("evidence_credibility_score")))
        credibility["batch_adjusted_score"] = credibility.get("final_score", credibility.get("evidence_credibility_score"))
        profile = dict(credibility.get("score_profile", {}))
        profile.update({
            "batch_FDR_applied": True,
            "batch_FDR_adjustment_applied": False,
            "batch_multiple_testing_method": BATCH_METHOD,
            "batch_hypothesis_count": total,
            "batch_family_id": family_id,
            "batch_support_semantics": "separate_batch_inference_local_rule_score_unchanged_not_a_probability",
        })
        credibility["score_profile"] = profile
        if network_primary and not supported:
            limitation = (
                "network_support_unevaluable_after_batch_FDR" if pattern_status == "unevaluable" else
                "network_support_not_significant_after_batch_FDR"
            )
            _append_unique(credibility, "limitations", limitation)
            _append_unique(credibility, "warnings", limitation)
            credibility.setdefault("limitation_severity", {})[limitation] = "high"
            result.mechanism.uncertainty = list(dict.fromkeys([*result.mechanism.uncertainty, limitation]))
    return results
