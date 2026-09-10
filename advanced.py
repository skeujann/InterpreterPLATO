"""Optional advanced explanations: dynamics, paths, graphs, counterfactuals."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from copy import deepcopy
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple
from hashlib import sha256
import json
import math

import networkx as nx

from .context import CancerContext, ContextScore, edge_distance, gene
from .mechanism import (
    MIN_ADEQUATE_RANDOMIZATIONS,
    MechanismResult,
    classify_mechanism,
    qualified_target_action_status,
)


PROCESS_GENES: Dict[str, Set[str]] = {
    "replication_stress": {"RRM1", "RRM2", "TYMS", "DHFR", "POLA1", "ATR", "CHEK1"},
    "dna_damage": {"TOP1", "TOP2A", "PARP1", "ATM", "ATR", "BRCA1", "BRCA2"},
    "dna_repair": {"PARP1", "BRCA1", "BRCA2", "MGMT", "RAD51", "WEE1"},
    "mitotic_stress": {"TUBB", "TUBA1B", "AURKA", "AURKB", "PLK1", "KIF11"},
    "survival_signaling": {"EGFR", "ERBB2", "SRC", "AKT1", "PIK3CA", "MTOR", "MAPK1"},
    "apoptosis_regulation": {"BCL2", "BCL2L1", "MCL1", "MDM2", "TP53", "BAX"},
    "epigenetic_reprogramming": {"DNMT1", "DNMT3A", "DNMT3B", "HDAC1", "HDAC2", "EZH2"},
}


@dataclass
class AdvancedResult:
    dynamic_hypotheses: List[Dict[str, Any]] = field(default_factory=list)
    local_paths: List[Dict[str, Any]] = field(default_factory=list)
    evidence_graph: Dict[str, Any] = field(default_factory=dict)
    counterfactuals: List[Dict[str, Any]] = field(default_factory=list)
    semantics: Dict[str, str] = field(default_factory=lambda: {
        "dynamic_hypotheses": "rule-generated hypotheses requiring temporal validation",
        "local_paths": "weighted network evidence paths, not causal chains",
        "evidence_graph": "provenance organization, not independent validation",
        "counterfactuals": "interpreter sensitivity analyses, not biological causal effects",
    })

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def _processes(profile: Sequence[ContextScore], top_k: int = 10) -> Dict[str, List[str]]:
    selected = {gene(item.target_id) for item in profile[:top_k]}
    return {
        process: sorted(selected & members)
        for process, members in PROCESS_GENES.items()
        if selected & members
    }


def dynamic_hypotheses(
    drug1_name: str,
    drug2_name: str,
    drug1_targets: Sequence[ContextScore],
    drug2_targets: Sequence[ContextScore],
    *,
    top_k: int = 10,
) -> List[Dict[str, Any]]:
    left = _processes(drug1_targets, top_k=top_k)
    right = _processes(drug2_targets, top_k=top_k)
    hypotheses: List[Dict[str, Any]] = []
    profiles = {drug1_name: drug1_targets, drug2_name: drug2_targets}
    rules = [
        (
            "replication_stress", "mitotic_stress", "cell_cycle_program_co_modulation",
            "Replication-associated and mitotic programs are both represented; direction and schedule require validation.",
        ),
        (
            "dna_damage", "dna_repair", "damage_repair_program_co_modulation",
            "DNA-damage-associated and repair programs are both represented; reduced repair capacity is not inferred without directional evidence.",
        ),
        (
            "survival_signaling", "apoptosis_regulation", "survival_apoptosis_program_co_modulation",
            "Survival-signaling and apoptosis-regulatory programs are both represented; an apoptotic-threshold change is not inferred.",
        ),
        (
            "epigenetic_reprogramming", "dna_damage", "epigenetic_damage_program_co_modulation",
            "Epigenetic and DNA-damage-associated programs are both represented; sensitization requires experimental validation.",
        ),
    ]
    for first, second, neutral_state_type, neutral_explanation in rules:
        orientations = []
        if first in left and second in right:
            orientations.append((drug1_name, drug2_name, left[first], right[second]))
        if first in right and second in left:
            orientations.append((drug2_name, drug1_name, right[first], left[second]))
        for priming_drug, responding_drug, priming_targets, responding_targets in orientations:
            priming_index = {gene(item.target_id): item for item in profiles[priming_drug]}
            responding_index = {gene(item.target_id): item for item in profiles[responding_drug]}
            priming_directions = {
                target: qualified_target_action_status(priming_index[target])
                for target in priming_targets if target in priming_index
            }
            responding_directions = {
                target: qualified_target_action_status(responding_index[target])
                for target in responding_targets if target in responding_index
            }
            priming_inhibitory = bool(priming_directions) and all(
                value == "inhibitory" for value in priming_directions.values()
            )
            responding_inhibitory = bool(responding_directions) and all(
                value == "inhibitory" for value in responding_directions.values()
            )
            state_type = neutral_state_type
            explanation = neutral_explanation
            directional_claim_supported = False
            if first == "dna_damage" and second == "dna_repair" and priming_inhibitory and responding_inhibitory:
                state_type = "directional_damage_repair_hypothesis"
                explanation = (
                    "Traceable inhibitory local actions occur in both DNA-damage-associated and repair programs; "
                    "reduced repair capacity remains a temporal hypothesis."
                )
                directional_claim_supported = True
            elif first == "survival_signaling" and priming_inhibitory:
                explanation = (
                    "Traceable inhibitory local action supports survival-signaling suppression, paired with "
                    "non-directional apoptosis-regulatory evidence; an apoptotic-threshold change is not inferred."
                )
                directional_claim_supported = True
            elif first == "replication_stress" and priming_inhibitory and responding_inhibitory:
                state_type = "directional_cell_cycle_program_hypothesis"
                explanation = (
                    "Traceable inhibitory local actions span replication-associated and mitotic programs; "
                    "sequential vulnerability remains a schedule-dependent hypothesis."
                )
                directional_claim_supported = True
            hypotheses.append({
                "state_type": state_type,
                "priming_drug": priming_drug,
                "responding_drug": responding_drug,
                "priming_process": first,
                "responding_process": second,
                "priming_targets": priming_targets,
                "responding_targets": responding_targets,
                "confidence": 0.65 if directional_claim_supported else 0.5,
                "hypothesis_score": 0.65 if directional_claim_supported else 0.5,
                "score_semantics": "uncalibrated_rule_strength",
                "calibrated_probability": None,
                "evidence_type": "dynamic_hypothesis",
                "description": explanation,
                "direction_evidence": {
                    "priming_targets": priming_directions,
                    "responding_targets": responding_directions,
                    "local_action_required_for_directional_claim": True,
                },
                "directional_claim_supported": directional_claim_supported,
                "warning": "schedule_dependency_possible",
                "warnings": ["schedule_dependency_possible", "temporal_evidence_not_supplied"],
            })
    return hypotheses


def local_evidence_paths(
    drug1_name: str,
    drug2_name: str,
    drug1_targets: Sequence[ContextScore],
    drug2_targets: Sequence[ContextScore],
    context: CancerContext,
    *,
    top_k: int = 10,
    max_length: int = 5,
    disease_label: str = "",
) -> List[Dict[str, Any]]:
    graph = context.ppi
    if graph.number_of_nodes() == 0:
        return []
    disease = context.disease_genes(disease_label)
    destinations = [item for item in disease if item in graph]
    if not destinations:
        return []
    max_degree = max((degree for _, degree in graph.degree()), default=1)
    candidates: List[Dict[str, Any]] = []
    try:
        _, nearest_disease_paths = nx.multi_source_dijkstra(
            graph, destinations, weight=edge_distance
        )
    except (nx.NetworkXError, nx.NodeNotFound):
        return []
    for drug_name, profile in ((drug1_name, drug1_targets), (drug2_name, drug2_targets)):
        for target in profile[:top_k]:
            source = gene(target.target_id)
            disease_to_source = nearest_disease_paths.get(source)
            if disease_to_source is None:
                continue
            path = list(reversed(disease_to_source))
            destination = path[-1]
            if len(path) - 1 > max_length:
                continue
            edge_pairs = list(zip(path, path[1:]))
            edge_data = [graph.get_edge_data(left, right, default={}) for left, right in edge_pairs]
            edge_costs = [edge_distance(left, right, item) for (left, right), item in zip(edge_pairs, edge_data)]
            edge_confidences = [
                item.get("confidence") if item.get("confidence_available") is True else None
                for item in edge_data
            ]
            edge_support_factors = [
                min(1.0, max(0.0, float(value))) if value is not None and math.isfinite(float(value)) else 0.5
                for value in edge_confidences
            ]
            weighted_cost = sum(edge_costs)
            confidence_geomean = (
                math.exp(sum(math.log(max(value, 1e-6)) for value in edge_support_factors) / len(edge_support_factors))
                if edge_support_factors else 1.0
            )
            hub_values = [
                math.log1p(graph.degree(node)) / math.log1p(max_degree)
                for node in path[1:-1]
            ] if max_degree > 0 else []
            hub_penalty = float(sum(hub_values) / len(hub_values)) if hub_values else 0.0
            hop_factor = 1.0 / (1.0 + len(path) - 1)
            score = target.final_score * confidence_geomean * hop_factor * math.exp(-0.25 * hub_penalty)
            candidates.append({
                "drug": drug_name,
                "source_target": source,
                "destination": destination,
                "path": path,
                "length": len(path) - 1,
                "bridge_genes": path[1:-1],
                "score": float(score),
                "weighted_cost": float(weighted_cost),
                "edge_confidences": edge_confidences,
                "edge_support_factors": edge_support_factors,
                "edge_confidence_geomean": float(confidence_geomean),
                "edge_confidence_coverage": sum(value is not None for value in edge_confidences) / max(1, len(edge_confidences)),
                "hub_penalty": hub_penalty,
                "score_components": {
                    "target_score": target.final_score,
                    "edge_confidence_geomean": float(confidence_geomean),
                    "hop_factor": hop_factor,
                    "hub_factor": math.exp(-0.25 * hub_penalty),
                },
                "weighting_method": "nearest_disease_confidence_path_with_hub_penalty",
                "weight_semantics": "confidence_weighted_when_available; topology-only edges use a conservative support factor",
                "source_evidence_status": target.evidence_status,
            })
    candidates.sort(key=lambda item: item["score"], reverse=True)
    unique: List[Dict[str, Any]] = []
    seen = set()
    for item in candidates:
        key = (item["drug"], tuple(item["path"]))
        if key not in seen:
            seen.add(key)
            unique.append(item)
    if top_k <= 1:
        return unique[:max(0, top_k)]
    selected: List[Dict[str, Any]] = []
    # With room for both drugs, reserve one best path per drug and then fill
    # the remaining slots globally by score.  A one-slot result is purely the
    # highest scoring path and is not biased toward drug1.
    for drug_name in (drug1_name, drug2_name):
        first = next((item for item in unique if item["drug"] == drug_name), None)
        if first is not None and first not in selected:
            selected.append(first)
    for item in unique:
        if item not in selected:
            selected.append(item)
        if len(selected) >= top_k:
            break
    return selected[:top_k]


def build_evidence_graph(
    drug1_name: str,
    drug2_name: str,
    drug1_targets: Sequence[ContextScore],
    drug2_targets: Sequence[ContextScore],
    context: CancerContext,
    paths: Sequence[Mapping[str, Any]],
    dynamics: Sequence[Mapping[str, Any]],
) -> Dict[str, Any]:
    nodes: Dict[str, Dict[str, Any]] = {}
    edges: List[Dict[str, Any]] = []

    def add_node(node_id: str, node_type: str, **attrs: Any) -> None:
        nodes.setdefault(node_id, {"id": node_id, "type": node_type, **attrs})

    for drug_name, profile in ((drug1_name, drug1_targets), (drug2_name, drug2_targets)):
        drug_id = f"drug:{drug_name}"
        add_node(drug_id, "drug", name=drug_name)
        for target in profile:
            target_id = f"gene:{gene(target.target_id)}"
            add_node(target_id, "gene", name=target.target_name)
            relation = "known_and_predicted_target" if target.evidence_status == "known_and_predicted" else (
                "known_target" if target.evidence_status == "known_only" else "PLATO_predicted_target"
            )
            edges.append({
                "source": drug_id,
                "target": target_id,
                "type": relation,
                "confidence": target.evidence_score,
                "context_score": target.final_score,
                "provenance": target.source_evidence,
            })
            if target.source_evidence.get("plato"):
                plato_details = target.source_evidence.get("plato")
                evidence_identity = json.dumps(
                    {"drug": drug_name, "target": gene(target.target_id), "details": plato_details},
                    ensure_ascii=False, sort_keys=True, separators=(",", ":"), default=str,
                )
                evidence_hash = sha256(evidence_identity.encode("utf-8")).hexdigest()[:16]
                evidence_id = f"evidence:PLATO:{evidence_hash}:{gene(target.target_id)}"
                add_node(
                    evidence_id,
                    "prediction_evidence",
                    provider="PLATO",
                    drug=drug_name,
                    details=plato_details,
                )
                edges.append({"source": evidence_id, "target": target_id, "type": "supports_prediction"})

    for name, members in context.pathways_for(
        [item.target_id for item in drug1_targets] + [item.target_id for item in drug2_targets]
    ).items():
        pathway_id = f"pathway:{name}"
        add_node(
            pathway_id,
            "pathway",
            name=name,
            source=context.pathway_sources.get(name, ""),
            version=context.pathway_versions.get(name, ""),
            taxid=context.pathway_taxids.get(name, ""),
            identifier_namespace=context.pathway_identifier_namespaces.get(name, ""),
            mapping_version=context.pathway_mapping_versions.get(name, ""),
            content_sha256=context.pathway_content_hashes.get(name, ""),
            artifact_kind=context.pathway_artifact_kinds.get(name, ""),
            artifact_sha256=context.pathway_artifact_hashes.get(name, ""),
        )
        for member in members:
            add_node(f"gene:{member}", "gene", name=member)
            edges.append({"source": f"gene:{member}", "target": pathway_id, "type": "belongs_to_pathway"})

    for disease_name, members in context.disease_sets.items():
        disease_id = f"disease_module:{disease_name}"
        add_node(disease_id, "disease_module", name=disease_name)
        for member in sorted(members):
            if f"gene:{member}" in nodes:
                edges.append({"source": f"gene:{member}", "target": disease_id, "type": "belongs_to_disease_module"})

    for item in paths:
        sequence = list(item.get("path", []))
        for left, right in zip(sequence, sequence[1:]):
            add_node(f"gene:{left}", "gene", name=left)
            add_node(f"gene:{right}", "gene", name=right)
            attributes = context.ppi.get_edge_data(left, right, default={})
            distance = edge_distance(left, right, attributes)
            edge = {
                "source": f"gene:{left}",
                "target": f"gene:{right}",
                "type": "PPI",
                "path_score": item.get("score", 0.0),
                "confidence": attributes.get("confidence") if attributes.get("confidence_available") is True else None,
                "confidence_available": attributes.get("confidence_available", False),
                "weight_semantics": "confidence" if attributes.get("confidence_available") is True else "topology_or_explicit_distance_only",
                "distance": distance if math.isfinite(distance) else None,
                "source_dataset": context.resource_metadata.get("ppi", {}).get("source", ""),
            }
            if edge not in edges:
                edges.append(edge)

    for index, item in enumerate(dynamics):
        state_id = f"state:{item.get('state_type', index)}"
        add_node(state_id, "cell_state", name=item.get("state_type", "dynamic_state"))
        edges.append({
            "source": f"drug:{item.get('priming_drug')}",
            "target": state_id,
            "type": "hypothesized_priming_relation",
            "evidence_type": "rule_generated_hypothesis",
            "confidence": item.get("confidence", 0.0),
            "directional_claim_supported": bool(item.get("directional_claim_supported")),
        })
        edges.append({
            "source": state_id,
            "target": f"drug:{item.get('responding_drug')}",
            "type": "hypothesized_response_relation",
            "evidence_type": "rule_generated_hypothesis",
            "confidence": item.get("confidence", 0.0),
            "directional_claim_supported": bool(item.get("directional_claim_supported")),
        })
    return {
        "nodes": list(nodes.values()),
        "edges": edges,
        "semantics": "provenance graph; graph presence alone is not independent validation",
    }


def _counterfactual_comparison(
    baseline: MechanismResult,
    perturbed: MechanismResult,
    *,
    network_randomizations: int,
    random_seed: int,
) -> Dict[str, Any]:
    """Audit actual null realizations before interpreting a sensitivity comparison."""
    reasons: List[str] = []
    audits: Dict[str, Dict[str, Any]] = {}
    protocol_fields = ("method_version", "distance_method", "null_method")

    def finite_number(value: Any) -> Optional[float]:
        if isinstance(value, bool):
            return None
        try:
            number = float(value)
        except (TypeError, ValueError, OverflowError):
            return None
        return number if math.isfinite(number) else None

    for label, result in (("baseline", baseline), ("perturbed", perturbed)):
        network = result.metrics.get("cheng", {})
        network = network if isinstance(network, Mapping) else {}
        requested = finite_number(network.get("n_randomizations"))
        accepted = finite_number(network.get("n_valid_randomizations"))
        seed = finite_number(network.get("random_seed"))
        null_stds = {}
        for arm in ("A", "B"):
            null = network.get(f"null_{arm}", {})
            null_stds[arm] = finite_number(null.get("std")) if isinstance(null, Mapping) else None
        audit = {
            "valid": network.get("valid") is True,
            "quality_status": network.get("quality_status"),
            "n_randomizations": requested,
            "n_valid_randomizations": accepted,
            "valid_fraction": accepted / requested if requested and accepted is not None else None,
            "random_seed": seed,
            "null_stds": null_stds,
            **{key: network.get(key) for key in protocol_fields},
        }
        audits[label] = audit
        if not audit["valid"]:
            reasons.append(f"{label}:invalid_network_statistics")
        if audit["quality_status"] != "adequate":
            reasons.append(f"{label}:inadequate_null_quality")
        if requested is None or not requested.is_integer() or requested < MIN_ADEQUATE_RANDOMIZATIONS:
            reasons.append(f"{label}:insufficient_requested_randomizations")
        if accepted is None or not accepted.is_integer() or accepted < MIN_ADEQUATE_RANDOMIZATIONS:
            reasons.append(f"{label}:insufficient_valid_randomizations")
        if requested is None or accepted is None or requested <= 0 or not 0.9 * requested <= accepted <= requested:
            reasons.append(f"{label}:invalid_acceptance_fraction")
        if any(value is None or value <= 0 for value in null_stds.values()):
            reasons.append(f"{label}:degenerate_or_missing_null_variance")
        if requested != network_randomizations:
            reasons.append(f"{label}:configured_budget_mismatch")
        if seed is None or not seed.is_integer() or seed != random_seed:
            reasons.append(f"{label}:configured_seed_mismatch")
        for key in protocol_fields:
            if not isinstance(audit[key], str) or not audit[key].strip():
                reasons.append(f"{label}:missing_{key}")
    for key in (*protocol_fields, "random_seed", "n_randomizations"):
        if audits["baseline"][key] != audits["perturbed"][key]:
            reasons.append(f"protocol_mismatch:{key}")
    evaluable = not reasons
    descriptive_changed = perturbed.mechanism_type != baseline.mechanism_type
    return {
        "descriptive_mechanism_changed": descriptive_changed,
        "mechanism_changed": descriptive_changed if evaluable else None,
        "comparison_evaluable": evaluable,
        "comparison_status": "evaluated" if evaluable else "not_evaluable",
        "comparison_audit": {
            "reason_codes": reasons,
            "baseline": audits["baseline"],
            "perturbed": audits["perturbed"],
            "configured_network_randomizations": network_randomizations,
            "configured_random_seed": random_seed,
            "minimum_accepted_randomizations": MIN_ADEQUATE_RANDOMIZATIONS,
            "minimum_accepted_fraction": 0.9,
            "scope": "interpreter_sensitivity_not_model_fidelity_or_biological_causality",
        },
    }


def _masked_path_context(context: CancerContext, bridges: Sequence[str]) -> CancerContext:
    """Build an isolated resource snapshot with fresh caches for graph masking."""
    masked = CancerContext()
    for attribute in (
        "ppi", "expression", "pathways", "pathway_sources", "pathway_versions",
        "pathway_taxids", "pathway_identifier_namespaces", "pathway_mapping_versions",
        "pathway_content_hashes", "pathway_artifact_kinds", "pathway_artifact_hashes",
        "disease_sets", "default_disease_genes", "default_disease_label",
        "resource_metadata", "is_demo_context",
    ):
        setattr(masked, attribute, deepcopy(getattr(context, attribute)))
    masked.ppi.remove_nodes_from([node for node in bridges if node in masked.ppi])
    return masked


def counterfactual_analysis(
    drug1_targets: Sequence[ContextScore],
    drug2_targets: Sequence[ContextScore],
    context: CancerContext,
    baseline: MechanismResult,
    *,
    max_targets_per_drug: int = 5,
    top_k: int = 10,
    disease_label: str = "",
    network_randomizations: int = 1000,
    random_seed: int = 42,
    min_network_coverage: float = 0.5,
) -> List[Dict[str, Any]]:
    reports: List[Dict[str, Any]] = []
    for drug_key, left, right in (
        ("drug1", drug1_targets, drug2_targets),
        ("drug2", drug2_targets, drug1_targets),
    ):
        for removed in left[:max_targets_per_drug]:
            retained = [item for item in left if item.target_id != removed.target_id]
            if not retained:
                continue
            options = {
                "top_k": top_k,
                "disease_label": disease_label,
                "network_randomizations": network_randomizations,
                "random_seed": random_seed,
                "min_network_coverage": min_network_coverage,
            }
            perturbed = classify_mechanism(retained, right, context, **options) if drug_key == "drug1" else classify_mechanism(right, retained, context, **options)
            reports.append({
                "perturbation": "remove_target",
                "kind": "target",
                "drug": drug_key,
                "target_id": removed.target_id,
                "target_name": removed.target_name,
                "target_evidence_status": removed.evidence_status,
                "before_mechanism": baseline.mechanism_type,
                "after_mechanism": perturbed.mechanism_type,
                "before_confidence": baseline.confidence,
                "after_confidence": perturbed.confidence,
                "confidence_delta": perturbed.confidence - baseline.confidence,
                "score_delta": perturbed.confidence - baseline.confidence,
                **_counterfactual_comparison(
                    baseline, perturbed,
                    network_randomizations=network_randomizations, random_seed=random_seed,
                ),
                "analysis_semantics": "interpreter_rule_sensitivity_not_model_fidelity_or_biological_causality",
                "network_randomizations": network_randomizations,
            })
    # Module suppression: remove target members of an implicated pathway from
    # the frozen fused target profiles, without re-querying PLATO.
    implicated = context.pathways_for(
        [item.target_id for item in drug1_targets] + [item.target_id for item in drug2_targets]
    )
    for pathway, members in list(implicated.items())[:5]:
        member_set = set(members)
        retained1 = [item for item in drug1_targets if gene(item.target_id) not in member_set]
        retained2 = [item for item in drug2_targets if gene(item.target_id) not in member_set]
        if not retained1 or not retained2:
            continue
        perturbed = classify_mechanism(
            retained1, retained2, context,
            top_k=top_k,
            disease_label=disease_label,
            network_randomizations=network_randomizations,
            random_seed=random_seed,
            min_network_coverage=min_network_coverage,
        )
        reports.append({
            "perturbation": "suppress_module",
            "kind": "module",
            "module": pathway,
            "removed_targets": sorted(member_set),
            "before_mechanism": baseline.mechanism_type,
            "after_mechanism": perturbed.mechanism_type,
            "before_confidence": baseline.confidence,
            "after_confidence": perturbed.confidence,
            "confidence_delta": perturbed.confidence - baseline.confidence,
            "score_delta": perturbed.confidence - baseline.confidence,
            **_counterfactual_comparison(
                baseline, perturbed,
                network_randomizations=network_randomizations, random_seed=random_seed,
            ),
            "analysis_semantics": "interpreter_rule_sensitivity_not_model_fidelity_or_biological_causality",
            "network_randomizations": network_randomizations,
        })
    reports.sort(
        key=lambda item: (bool(item.get("mechanism_changed")), abs(item["confidence_delta"])),
        reverse=True,
    )
    return reports


def explain_advanced(
    drug1_name: str,
    drug2_name: str,
    drug1_targets: Sequence[ContextScore],
    drug2_targets: Sequence[ContextScore],
    context: CancerContext,
    mechanism: MechanismResult,
    *,
    top_k: int = 10,
    disease_label: str = "",
    sensitivity_randomizations: int = 1000,
    random_seed: int = 42,
    min_network_coverage: float = 0.5,
) -> AdvancedResult:
    dynamics = dynamic_hypotheses(
        drug1_name, drug2_name, drug1_targets, drug2_targets, top_k=top_k
    )
    paths = local_evidence_paths(
        drug1_name, drug2_name, drug1_targets, drug2_targets, context,
        top_k=top_k,
        disease_label=disease_label,
    )
    graph = build_evidence_graph(drug1_name, drug2_name, drug1_targets, drug2_targets, context, paths, dynamics)
    counterfactuals = counterfactual_analysis(
        drug1_targets,
        drug2_targets,
        context,
        mechanism,
        top_k=top_k,
        disease_label=disease_label,
        network_randomizations=sensitivity_randomizations,
        random_seed=random_seed,
        min_network_coverage=min_network_coverage,
    )
    # Path masking uses a graph snapshot so predicted targets cannot silently
    # re-enter through a fresh PLATO resolution.
    for item in paths[:5]:
        bridges = list(item.get("bridge_genes", []))
        if not bridges:
            continue
        masked_context = _masked_path_context(context, bridges)
        perturbed = classify_mechanism(
            drug1_targets,
            drug2_targets,
            masked_context,
            top_k=top_k,
            disease_label=disease_label,
            network_randomizations=sensitivity_randomizations,
            random_seed=random_seed,
            min_network_coverage=min_network_coverage,
        )
        counterfactuals.append({
            "perturbation": "mask_path",
            "kind": "path",
            "path": item.get("path", []),
            "masked_bridge_genes": bridges,
            "before_mechanism": mechanism.mechanism_type,
            "after_mechanism": perturbed.mechanism_type,
            "before_confidence": mechanism.confidence,
            "after_confidence": perturbed.confidence,
            "confidence_delta": perturbed.confidence - mechanism.confidence,
            "score_delta": perturbed.confidence - mechanism.confidence,
            **_counterfactual_comparison(
                mechanism, perturbed,
                network_randomizations=sensitivity_randomizations, random_seed=random_seed,
            ),
            "analysis_semantics": "interpreter_network_path_sensitivity_not_model_fidelity_or_biological_causality",
            "network_randomizations": sensitivity_randomizations,
        })
    for index, item in enumerate(counterfactuals):
        node_id = f"counterfactual:{index}"
        graph["nodes"].append({"id": node_id, "type": "counterfactual", "details": item})
        graph["edges"].append({
            "source": node_id,
            "target": "mechanism:baseline",
            "type": "changes_mechanism" if item.get("mechanism_changed") else "tests_sensitivity_of",
            "confidence_delta": item.get("confidence_delta", 0.0),
        })
    graph["nodes"].append({"id": "mechanism:baseline", "type": "mechanism", "name": mechanism.mechanism_type})
    return AdvancedResult(dynamics, paths, graph, counterfactuals)
