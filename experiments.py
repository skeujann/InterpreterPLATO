"""Explicit, non-destructive engineering controls, applied AFTER input approval.

These are interventions on this rule-based interpreter, not on ScoreRI or cells.
Resources, PLATO trust, and eligibility criteria are never weakened here.
"""
from __future__ import annotations

from dataclasses import replace
from typing import Any, Mapping

POLICY_VERSION = "engineering_controls_v1"
VARIANTS = ("full", "no_expression_ranking", "raw_plato_order", "drop_top_target")


def validate_experiment(value: Any) -> dict:
    if not isinstance(value, Mapping) or set(value) - {"variant"}:
        raise ValueError("engineering_experiment accepts only a variant")
    variant = value.get("variant", "full")
    if variant not in VARIANTS:
        raise ValueError(f"Unknown engineering variant: {variant!r}")
    return {"variant": variant}


def ranking_weights(config: Mapping[str, Any]) -> dict:
    variant = validate_experiment(config.get("engineering_experiment", {}))["variant"]
    weights = dict(config["context_weights"])
    if variant == "no_expression_ranking":
        weights["expression"] = 0.0
    return weights


def apply_control(ranked1, ranked2, config):
    """Return fresh lists and a complete before/after target-set audit.

    raw_plato_order changes ORDER ONLY. final_score remains the descriptive
    context score; provider scores are never converted into probabilities.
    drop_top_target removes the first context-ranked target on BOTH sides.
    All downstream mechanism/network/pathway analyses consume retained lists.
    """
    variant = validate_experiment(config.get("engineering_experiment", {}))["variant"]
    audit = {"policy_version": POLICY_VERSION, "variant": variant,
             "not_a_model_or_biological_intervention": True,
             "resource_files_modified": False, "targets": {}}
    outputs = []
    for side, original in (("drug1", ranked1), ("drug2", ranked2)):
        ranked = list(original)
        if variant == "raw_plato_order":
            def provider_score(item):
                return float(item.source_evidence["plato"]["score"])
            ranked.sort(key=lambda item: (-provider_score(item), item.target_id))
            ranked = [replace(item, ordering_policy="raw_plato_score_descending",
                              ordering_score=provider_score(item),
                              scoring_notes=item.scoring_notes + [
                                  "engineering_order_control_context_scores_retained"])
                      for item in ranked]
        elif variant == "drop_top_target":
            if len(ranked) < 2:
                raise ValueError(f"drop_top_target requires >=2 eligible targets for {side}")
            ranked = ranked[1:]
        audit["targets"][side] = {
            "before": [item.target_id for item in original],
            "after": [item.target_id for item in ranked],
            "removed": [item.target_id for item in original
                        if item.target_id not in {target.target_id for target in ranked}],
        }
        outputs.append(ranked)
    audit["scope"] = {
        "full": "unmodified_context_ranking",
        "no_expression_ranking": "expression_weight_zero_only_annotations_and_quality_checks_retained",
        "raw_plato_order": "provider_order_only_context_scores_and_resources_retained",
        "drop_top_target": "one_context_top_target_removed_per_drug_downstream_recomputed",
    }[variant]
    return outputs[0], outputs[1], audit
