"""Compact drug-pair mechanism classification with explicit evidence."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple
from bisect import bisect_left
import math
import re

import networkx as nx
import numpy as np

from .context import CancerContext, ContextScore, edge_distance, gene


REPLICATION = {"RRM1", "RRM2", "TYMS", "DHFR", "POLA1", "ATR", "CHEK1"}
MITOSIS = {"TUBB", "TUBA1B", "AURKA", "AURKB", "PLK1", "KIF11"}
DNA_DAMAGE = {"TOP1", "TOP2A", "PARP1", "ATM", "ATR", "BRCA1", "BRCA2"}
DNA_REPAIR = {"PARP1", "BRCA1", "BRCA2", "MGMT", "RAD51", "WEE1"}
MIN_ADEQUATE_RANDOMIZATIONS = 1000
INHIBITORY_ACTION_TERMS = (
    "inhibitor", "inhibition", "antagonist", "blocker", "negative modulator",
    "negative allosteric", "degrader", "degradation", "suppressor", "inactivator", "inhibitory",
)
ACTIVATING_ACTION_TERMS = (
    "agonist", "activator", "activation", "positive modulator", "positive allosteric", "inducer", "activating",
)
QUANTITATIVE_ACTIVITY_TYPES = {"ic50", "ec50", "ac50", "xc50", "ki", "kd", "potency"}
MIN_PCHEMBL_FOR_QUALIFIED_DIRECTION = 5.0


def _action_direction(values: Sequence[str]) -> str:
    normalized = " | ".join(value.strip().lower() for value in values if value.strip())
    normalized = re.sub(r"[-_/]+", " ", normalized)
    controlled_terms = INHIBITORY_ACTION_TERMS + ACTIVATING_ACTION_TERMS
    if any(
        re.search(rf"\b(?:not(?:\s+an?)?|no|non)\s+{re.escape(term)}\b", normalized)
        for term in controlled_terms
    ):
        return "unknown"
    inverse_agonist = bool(re.search(r"\binverse\s+agonist\b", normalized))
    remainder = re.sub(r"\binverse\s+agonist\b", "", normalized)
    inhibitory = inverse_agonist or any(
        re.search(rf"\b{re.escape(term)}\b", remainder) for term in INHIBITORY_ACTION_TERMS
    )
    activating = any(
        re.search(rf"\b{re.escape(term)}\b", remainder) for term in ACTIVATING_ACTION_TERMS
    )
    if inhibitory and activating:
        return "conflict"
    if inhibitory:
        return "inhibitory"
    if activating:
        return "activating"
    return "unknown"


def _wilson_interval(successes: int, trials: int, z: float = 1.96) -> Optional[List[float]]:
    if trials <= 0:
        return None
    proportion = successes / trials
    denominator = 1.0 + (z * z) / trials
    centre = (proportion + (z * z) / (2.0 * trials)) / denominator
    margin = z * math.sqrt(
        proportion * (1.0 - proportion) / trials + (z * z) / (4.0 * trials * trials)
    ) / denominator
    return [max(0.0, centre - margin), min(1.0, centre + margin)]


def target_action_status(item: ContextScore) -> str:
    """Reduce local/PLATO action assertions to a conservative direction."""
    evidence = item.source_evidence or {}
    local_values: List[str] = []
    local = evidence.get("local")
    if isinstance(local, Mapping):
        local_values.extend(str(local.get(key) or "") for key in ("action", "action_type", "mechanism_of_action"))
        assertions = local.get("assertions", [])
        if isinstance(assertions, list):
            for assertion in assertions:
                if isinstance(assertion, Mapping):
                    local_values.extend(str(assertion.get(key) or "") for key in ("action", "action_type", "mechanism_of_action"))
    plato = evidence.get("plato")
    plato_values = [str(plato.get("action") or "")] if isinstance(plato, Mapping) else []

    local_direction = _action_direction(local_values)
    plato_direction = _action_direction(plato_values)
    if local_direction == "conflict":
        return "conflict"
    if (
        local_direction in {"inhibitory", "activating"}
        and plato_direction in {"inhibitory", "activating"}
        and local_direction != plato_direction
    ):
        return "conflict"
    # PLATO action is retained as provenance but cannot establish direction
    # until a release-specific action field is validated.
    return local_direction


def qualified_target_action_status(item: ContextScore) -> str:
    """Direction from the same artifact-backed quantitative assertion only."""
    evidence = item.source_evidence or {}
    local = evidence.get("local")
    if not isinstance(local, Mapping) or local.get("input_provenance_status") != "sqlite_artifact_audited":
        return "unknown"
    provenance = local.get("provenance") if isinstance(local.get("provenance"), Mapping) else {}
    provenance_hash = str(provenance.get("database_sha256") or "").lower()
    try:
        schema_version = int(provenance.get("schema_user_version") or 0)
    except (TypeError, ValueError):
        return "unknown"
    if not (
        len(provenance_hash) == 64
        and all(character in "0123456789abcdef" for character in provenance_hash)
        and schema_version > 0
        and str(provenance.get("dataset_version") or "").strip()
        and str(provenance.get("mapping_version") or "").strip()
    ):
        return "unknown"
    assertions = local.get("assertions", [])
    directions: List[str] = []
    if isinstance(assertions, list):
        for assertion in assertions:
            if not isinstance(assertion, Mapping):
                continue
            known_action = str(assertion.get("known_action") or "").strip().lower()
            if known_action not in {"yes", "true", "1", "known"}:
                continue
            activity_type = re.sub(r"[^a-z0-9]", "", str(assertion.get("activity_type") or "").lower())
            try:
                potency = float(assertion.get("pchembl_value"))
            except (TypeError, ValueError, OverflowError):
                continue
            if not (
                str(assertion.get("source") or "").strip()
                and activity_type in QUANTITATIVE_ACTIVITY_TYPES
                and math.isfinite(potency)
                and MIN_PCHEMBL_FOR_QUALIFIED_DIRECTION <= potency <= 14.0
                and str(assertion.get("stable_record_id") or "").strip()
                and str(assertion.get("provenance_hash") or "").lower() == provenance_hash
            ):
                continue
            direction = _action_direction([
                str(assertion.get(key) or "")
                for key in ("action", "action_type", "mechanism_of_action")
            ])
            if direction != "unknown":
                directions.append(direction)
    unique = set(directions)
    if "conflict" in unique or {"inhibitory", "activating"} <= unique:
        return "conflict"
    return next(iter(unique)) if len(unique) == 1 else "unknown"


def _action_summary(items: Sequence[ContextScore]) -> Dict[str, Any]:
    statuses = {gene(item.target_id): target_action_status(item) for item in items}
    qualified_statuses = {gene(item.target_id): qualified_target_action_status(item) for item in items}
    total = max(1, len(items))
    return {
        "by_target": statuses,
        "direction_known_fraction": sum(value in {"inhibitory", "activating"} for value in statuses.values()) / total,
        "inhibitory_fraction": sum(value == "inhibitory" for value in statuses.values()) / total,
        "activating_fraction": sum(value == "activating" for value in statuses.values()) / total,
        "conflict_count": sum(value == "conflict" for value in statuses.values()),
        "qualified_by_target": qualified_statuses,
        "qualified_direction_known_fraction": sum(
            value in {"inhibitory", "activating"} for value in qualified_statuses.values()
        ) / total,
        "qualified_conflict_count": sum(value == "conflict" for value in qualified_statuses.values()),
    }


@dataclass
class MechanismResult:
    mechanism_type: str
    confidence: float
    description: str
    shared_targets: List[str] = field(default_factory=list)
    shared_pathways: List[str] = field(default_factory=list)
    complementary_pathways: List[str] = field(default_factory=list)
    metrics: Dict[str, Any] = field(default_factory=dict)
    evidence: List[Dict[str, Any]] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    score_semantics: str = "uncalibrated_evidence_strength"
    calibrated_probability: Optional[float] = None
    credibility_score: Optional[float] = None
    coverage: Dict[str, Any] = field(default_factory=dict)
    uncertainty: List[str] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


def _target_set(profile: Sequence[ContextScore], top_k: int) -> Set[str]:
    return {gene(item.target_id) for item in profile[:top_k] if gene(item.target_id)}


def _directed_closest_distance(
    graph: nx.Graph,
    left: Set[str],
    right: Set[str],
    distance_lookup: Optional[Callable[[str, str], Optional[float]]] = None,
) -> Optional[float]:
    mapped_left = {item for item in left if item in graph}
    mapped_right = {item for item in right if item in graph}
    if not mapped_left or not mapped_right:
        return None
    if distance_lookup is not None:
        values = []
        for source in mapped_left:
            nearest = min(
                (
                    value for value in (
                        distance_lookup(source, target) for target in mapped_right
                    )
                    if value is not None
                ),
                default=None,
            )
            if nearest is not None:
                values.append(float(nearest))
    else:
        try:
            distances = nx.multi_source_dijkstra_path_length(graph, mapped_right, weight=edge_distance)
        except (nx.NetworkXError, ValueError):
            return None
        values = [float(distances[item]) for item in mapped_left if item in distances]
    return float(np.mean(values)) if values else None


def _mean_cross_distance(
    graph: nx.Graph,
    left: Set[str],
    right: Set[str],
    distance_lookup: Optional[Callable[[str, str], Optional[float]]] = None,
) -> Optional[float]:
    forward = _directed_closest_distance(graph, left, right, distance_lookup)
    reverse = _directed_closest_distance(graph, right, left, distance_lookup)
    values = [value for value in (forward, reverse) if value is not None]
    return float(np.mean(values)) if values else None


def _mean_internal_distance(
    graph: nx.Graph,
    nodes: Set[str],
    distance_lookup: Optional[Callable[[str, str], Optional[float]]] = None,
) -> Optional[float]:
    if len(nodes) < 2:
        return 0.0
    distances: List[float] = []
    mapped = {item for item in nodes if item in graph}
    for source in mapped:
        others = mapped - {source}
        if not others:
            continue
        if distance_lookup is not None:
            nearest = min(
                (
                    value for value in (distance_lookup(source, item) for item in others)
                    if value is not None
                ),
                default=None,
            )
        else:
            try:
                lengths = nx.single_source_dijkstra_path_length(graph, source, weight=edge_distance)
            except nx.NodeNotFound:
                continue
            nearest = min((lengths[item] for item in others if item in lengths), default=None)
        if nearest is not None:
            distances.append(float(nearest))
    return sum(distances) / len(distances) if distances else None


def _degree_matched_sample(
    graph: nx.Graph,
    template: Set[str],
    rng: np.random.Generator,
    degree_groups: Mapping[int, Sequence[str]],
    sorted_degrees: Sequence[int],
    nodes: Sequence[str],
) -> Tuple[Set[str], float]:
    selected: List[str] = []
    errors: List[float] = []

    def draw(group: Sequence[str], banned: Set[str]) -> Optional[str]:
        if not group:
            return None
        for _ in range(min(12, max(2, len(group)))):
            candidate = group[int(rng.integers(0, len(group)))]
            if candidate not in banned:
                return candidate
        return next((candidate for candidate in group if candidate not in banned), None)

    for original in sorted(template):
        original_degree = graph.degree(original)
        banned = set(template) | set(selected)
        candidate = draw(degree_groups.get(original_degree, []), banned)
        if candidate is None:
            position = bisect_left(sorted_degrees, original_degree)
            left_index, right_index = position - 1, position
            while candidate is None and (left_index >= 0 or right_index < len(sorted_degrees)):
                left_degree = sorted_degrees[left_index] if left_index >= 0 else None
                right_degree = sorted_degrees[right_index] if right_index < len(sorted_degrees) else None
                left_error = abs(math.log1p(original_degree) - math.log1p(left_degree)) if left_degree is not None else float("inf")
                right_error = abs(math.log1p(original_degree) - math.log1p(right_degree)) if right_degree is not None else float("inf")
                if left_error <= right_error:
                    candidate = draw(degree_groups.get(left_degree, []), banned)
                    left_index -= 1
                else:
                    candidate = draw(degree_groups.get(right_degree, []), banned)
                    right_index += 1
        if candidate is None:
            # Never recycle an observed/template target into its own null
            # sample.  If the graph cannot provide a disjoint degree-matched
            # draw, return an unusable sample so the caller downgrades the
            # randomization quality instead of contaminating the null.
            candidate = draw(nodes, set(template) | set(selected))
        if candidate is None:
            return set(), float("inf")
        selected.append(candidate)
        errors.append(abs(math.log1p(original_degree) - math.log1p(graph.degree(candidate))))
    return set(selected), float(np.mean(errors)) if errors else 0.0


def _network_set_coverage(nodes: Set[str], graph: nx.Graph, distances: Optional[Mapping[str, float]] = None) -> Dict[str, Any]:
    mapped = {item for item in nodes if item in graph}
    reachable = mapped if distances is None else {item for item in mapped if item in distances}
    return {
        "input": len(nodes),
        "mapped": len(mapped),
        "mapping_fraction": len(mapped) / max(1, len(nodes)),
        "reachable": len(reachable),
        "reachability_fraction": len(reachable) / max(1, len(mapped)),
        "effective_fraction": len(reachable) / max(1, len(nodes)),
    }


def _cheng_metrics(
    graph: nx.Graph,
    left: Set[str],
    right: Set[str],
    disease: Set[str],
    samples: int = 1000,
    random_seed: int = 42,
    min_mapping_coverage: float = 0.5,
    max_degree_match_log_error: float = 0.35,
    disease_distances: Optional[Mapping[str, float]] = None,
    distance_lookup: Optional[Callable[[str, str], Optional[float]]] = None,
    degree_sampling_index: Optional[Mapping[str, Any]] = None,
) -> Dict[str, Any]:
    mapped_left = {item for item in left if item in graph}
    mapped_right = {item for item in right if item in graph}
    mapped_disease = {item for item in disease if item in graph}
    coverage = {
        "drug1": _network_set_coverage(left, graph),
        "drug2": _network_set_coverage(right, graph),
        "disease_module": _network_set_coverage(disease, graph),
    }
    if not mapped_left or not mapped_right or not mapped_disease:
        return {
            "valid": False,
            "pattern": "insufficient_network_support",
            "quality_status": "insufficient",
            "network_coverage": coverage,
            "warnings": ["network_mapping_incomplete"],
        }
    disease_distances = (
        disease_distances
        if disease_distances is not None
        else nx.multi_source_dijkstra_path_length(graph, mapped_disease, weight=edge_distance)
    )
    coverage["drug1"] = _network_set_coverage(left, graph, disease_distances)
    coverage["drug2"] = _network_set_coverage(right, graph, disease_distances)
    observed_values_a = [disease_distances[item] for item in mapped_left if item in disease_distances]
    observed_values_b = [disease_distances[item] for item in mapped_right if item in disease_distances]
    observed_a = float(np.mean(observed_values_a)) if observed_values_a else None
    observed_b = float(np.mean(observed_values_b)) if observed_values_b else None
    d_aa = _mean_internal_distance(graph, mapped_left, distance_lookup)
    d_bb = _mean_internal_distance(graph, mapped_right, distance_lookup)
    d_ab = _mean_cross_distance(graph, mapped_left, mapped_right, distance_lookup)
    if None in (observed_a, observed_b, d_aa, d_bb, d_ab):
        return {
            "valid": False,
            "pattern": "insufficient_network_support",
            "quality_status": "insufficient",
            "network_coverage": coverage,
            "warnings": ["network_components_disconnected"],
        }

    requested = max(0, int(samples))
    rng = np.random.default_rng(random_seed)
    if degree_sampling_index is not None:
        graph_nodes = list(degree_sampling_index.get("nodes", ()))
        degree_groups = {
            int(degree): list(group)
            for degree, group in dict(degree_sampling_index.get("degree_groups", {})).items()
        }
        sorted_degrees = list(degree_sampling_index.get("sorted_degrees", ()))
    else:
        graph_nodes = list(graph.nodes())
        degree_groups = {}
        for node in graph_nodes:
            degree_groups.setdefault(graph.degree(node), []).append(node)
        sorted_degrees = sorted(degree_groups)
    random_a: List[float] = []
    random_b: List[float] = []
    match_errors_a: List[float] = []
    match_errors_b: List[float] = []
    for _ in range(requested):
        sample_a, error_a = _degree_matched_sample(
            graph, mapped_left, rng, degree_groups, sorted_degrees, graph_nodes
        )
        sample_b, error_b = _degree_matched_sample(
            graph, mapped_right, rng, degree_groups, sorted_degrees, graph_nodes
        )
        values_a = [disease_distances[item] for item in sample_a if item in disease_distances]
        values_b = [disease_distances[item] for item in sample_b if item in disease_distances]
        if values_a and values_b:
            random_a.append(float(np.mean(values_a)))
            random_b.append(float(np.mean(values_b)))
            match_errors_a.append(error_a)
            match_errors_b.append(error_b)
    n_valid = min(len(random_a), len(random_b))
    std_a = float(np.std(random_a, ddof=1)) if len(random_a) > 1 else 0.0
    std_b = float(np.std(random_b, ddof=1)) if len(random_b) > 1 else 0.0
    mean_a = float(np.mean(random_a)) if random_a else None
    mean_b = float(np.mean(random_b)) if random_b else None
    z_a = (float(observed_a) - mean_a) / std_a if mean_a is not None and std_a > 0 else None
    z_b = (float(observed_b) - mean_b) / std_b if mean_b is not None and std_b > 0 else None
    hits_a = 1 + sum(value <= float(observed_a) for value in random_a)
    hits_b = 1 + sum(value <= float(observed_b) for value in random_b)
    p_a = hits_a / (len(random_a) + 1) if random_a else None
    p_b = hits_b / (len(random_b) + 1) if random_b else None
    p_a_interval = _wilson_interval(hits_a, len(random_a) + 1) if random_a else None
    p_b_interval = _wilson_interval(hits_b, len(random_b) + 1) if random_b else None
    p_a_adjusted = min(1.0, 2.0 * p_a) if p_a is not None else None
    p_b_adjusted = min(1.0, 2.0 * p_b) if p_b is not None else None
    p_a_adjusted_interval = [min(1.0, 2.0 * value) for value in p_a_interval] if p_a_interval else None
    p_b_adjusted_interval = [min(1.0, 2.0 * value) for value in p_b_interval] if p_b_interval else None
    s_ab = float(d_ab) - (float(d_aa) + float(d_bb)) / 2.0
    near_a_raw = bool(
        z_a is not None and p_a_adjusted is not None and p_a_adjusted_interval
        and z_a < 0.0 and p_a_adjusted <= 0.05 and p_a_adjusted_interval[1] <= 0.05
    )
    near_b_raw = bool(
        z_b is not None and p_b_adjusted is not None and p_b_adjusted_interval
        and z_b < 0.0 and p_b_adjusted <= 0.05 and p_b_adjusted_interval[1] <= 0.05
    )
    mapping_ok = all(
        coverage[key]["mapping_fraction"] >= min_mapping_coverage
        for key in ("drug1", "drug2", "disease_module")
    ) and all(
        coverage[key]["effective_fraction"] >= min_mapping_coverage
        for key in ("drug1", "drug2")
    )
    randomization_ok = (
        requested >= MIN_ADEQUATE_RANDOMIZATIONS
        and n_valid >= math.ceil(0.9 * requested)
    )
    null_ok = std_a > 0 and std_b > 0
    mean_match_error_a = float(np.mean(match_errors_a)) if match_errors_a else None
    mean_match_error_b = float(np.mean(match_errors_b)) if match_errors_b else None
    degree_match_ok_a = mean_match_error_a is not None and mean_match_error_a <= max_degree_match_log_error
    degree_match_ok_b = mean_match_error_b is not None and mean_match_error_b <= max_degree_match_log_error
    degree_match_ok = degree_match_ok_a and degree_match_ok_b
    near_a = near_a_raw and degree_match_ok_a
    near_b = near_b_raw and degree_match_ok_b
    overlap = bool(mapped_left & mapped_right) or s_ab < 0.0
    if near_a and near_b:
        pattern = "P1" if overlap else "P2"
    elif near_a != near_b:
        pattern = "P4_edge" if overlap else "P4"
    elif overlap:
        pattern = "P5"
    else:
        pattern = "P6"
    warnings: List[str] = []
    if not mapping_ok:
        warnings.append("low_network_mapping_coverage")
    if not randomization_ok:
        warnings.append("insufficient_valid_randomizations")
    if requested < MIN_ADEQUATE_RANDOMIZATIONS:
        warnings.append("randomizations_below_adequate_minimum")
    if not null_ok:
        warnings.append("degenerate_degree_matched_null")
    if not degree_match_ok:
        warnings.append("poor_degree_match_quality")
    quality = "adequate" if mapping_ok and randomization_ok and null_ok and degree_match_ok else "limited"
    return {
        "valid": True,
        "pattern": pattern,
        "quality_status": quality,
        "statistically_supported": bool(near_a or near_b),
        "pair_interaction_supported": bool(near_a and near_b),
        "z_A": z_a,
        "z_B": z_b,
        "p_A": p_a,
        "p_B": p_b,
        "p_A_adjusted": p_a_adjusted,
        "p_B_adjusted": p_b_adjusted,
        "p_A_monte_carlo_95ci": p_a_interval,
        "p_B_monte_carlo_95ci": p_b_interval,
        "p_A_adjusted_monte_carlo_95ci": p_a_adjusted_interval,
        "p_B_adjusted_monte_carlo_95ci": p_b_adjusted_interval,
        "multiple_testing_method": "Bonferroni_for_two_drug_disease_tests",
        "s_AB": s_ab,
        "d_A_obs": observed_a,
        "d_B_obs": observed_b,
        "d_AA": d_aa,
        "d_BB": d_bb,
        "d_AB": d_ab,
        "n_randomizations": requested,
        "n_valid_randomizations": n_valid,
        "minimum_randomizations_for_adequate_quality": MIN_ADEQUATE_RANDOMIZATIONS,
        "random_seed": random_seed,
        "degree_match_mean_log_error": (
            (mean_match_error_a + mean_match_error_b) / 2.0
            if mean_match_error_a is not None and mean_match_error_b is not None else None
        ),
        "degree_match_mean_log_error_A": mean_match_error_a,
        "degree_match_mean_log_error_B": mean_match_error_b,
        "degree_match_ok_A": degree_match_ok_a,
        "degree_match_ok_B": degree_match_ok_b,
        "max_degree_match_log_error": max_degree_match_log_error,
        "null_A": {
            "mean": mean_a,
            "std": std_a,
            "q025": float(np.quantile(random_a, 0.025)) if random_a else None,
            "q975": float(np.quantile(random_a, 0.975)) if random_a else None,
        },
        "null_B": {
            "mean": mean_b,
            "std": std_b,
            "q025": float(np.quantile(random_b, 0.025)) if random_b else None,
            "q975": float(np.quantile(random_b, 0.975)) if random_b else None,
        },
        "network_coverage": coverage,
        "distance_method": "confidence_weighted_symmetric_closest_set",
        "null_method": "degree_matched_fixed_disease_module",
        "warnings": warnings,
    }


def _functional_pattern(left: Set[str], right: Set[str], shared_pathways: Sequence[str]) -> Optional[Tuple[str, float, str]]:
    if ((left & REPLICATION) and (right & MITOSIS)) or ((right & REPLICATION) and (left & MITOSIS)):
        return (
            "cell_cycle_phase_co_modulation",
            0.62,
            "The target profiles map to replication/S-phase and mitotic/M-phase programs; effect direction is not yet established.",
        )
    if ((left & DNA_DAMAGE) and (right & DNA_REPAIR)) or ((right & DNA_DAMAGE) and (left & DNA_REPAIR)):
        return (
            "damage_repair_co_modulation",
            0.60,
            "The target profiles map to DNA-damage and repair programs; reinforcement requires action-direction evidence.",
        )
    if shared_pathways:
        return (
            "shared_pathway_co_modulation",
            min(0.70, 0.52 + 0.03 * len(shared_pathways)),
            "Distinct targets from both drugs are independently enriched in one or more common pathways after BH-FDR correction; action direction is not established.",
        )
    return None


def classify_mechanism(
    drug1_targets: Sequence[ContextScore],
    drug2_targets: Sequence[ContextScore],
    context: CancerContext,
    *,
    top_k: int = 10,
    disease_label: str = "",
    network_randomizations: int = 1000,
    random_seed: int = 42,
    min_network_coverage: float = 0.5,
) -> MechanismResult:
    if int(top_k) <= 0:
        raise ValueError("top_k must be positive")
    requested_randomizations = int(network_randomizations)
    if requested_randomizations < 0:
        raise ValueError("network_randomizations must be non-negative")
    try:
        requested_mapping_coverage = float(min_network_coverage)
    except (TypeError, ValueError):
        raise ValueError("min_network_coverage must be a finite number") from None
    if not math.isfinite(requested_mapping_coverage) or not 0.0 <= requested_mapping_coverage <= 1.0:
        raise ValueError("min_network_coverage must be between 0 and 1")
    effective_mapping_coverage = max(0.5, requested_mapping_coverage)
    selected1 = list(drug1_targets[:top_k])
    selected2 = list(drug2_targets[:top_k])
    left = _target_set(selected1, top_k)
    right = _target_set(selected2, top_k)
    if not left or not right:
        return MechanismResult(
            mechanism_type="unknown",
            confidence=0.0,
            description="At least one drug has no resolved human targets.",
            warnings=["missing_target_profile"],
        )

    shared = sorted(left & right)
    union = left | right
    jaccard = len(shared) / len(union) if union else 0.0
    pathways_left = context.pathways_for(left)
    pathways_right = context.pathways_for(right)
    descriptive_shared_pathways = sorted(set(pathways_left) & set(pathways_right))
    enrichment_left_records = context.pathway_enrichment(left)
    enrichment_right_records = context.pathway_enrichment(right)
    enrichment_left = {
        str(item["pathway"]): item
        for item in enrichment_left_records
        if float(item.get("q", 1.0)) <= 0.05
    }
    enrichment_right = {
        str(item["pathway"]): item
        for item in enrichment_right_records
        if float(item.get("q", 1.0)) <= 0.05
    }
    shared_pathways = sorted(set(enrichment_left) & set(enrichment_right))
    unique_left = sorted(set(enrichment_left) - set(enrichment_right))
    unique_right = sorted(set(enrichment_right) - set(enrichment_left))
    complementary = [f"{a} <> {b}" for a in unique_left[:3] for b in unique_right[:3]][:5]
    functional = _functional_pattern(left, right, shared_pathways)
    action_left = _action_summary(selected1)
    action_right = _action_summary(selected2)

    def inhibitory_drivers(summary: Mapping[str, Any], drivers: Set[str]) -> bool:
        statuses = summary.get("qualified_by_target", {})
        selected_statuses = [statuses.get(item, "unknown") for item in drivers]
        return bool(selected_statuses) and all(value == "inhibitory" for value in selected_statuses)

    if functional and functional[0] == "cell_cycle_phase_co_modulation":
        if left & REPLICATION and right & MITOSIS:
            driver_left, driver_right = left & REPLICATION, right & MITOSIS
        else:
            driver_left, driver_right = left & MITOSIS, right & REPLICATION
        directional = bool(
            not (driver_left & driver_right)
            and inhibitory_drivers(action_left, driver_left)
            and inhibitory_drivers(action_right, driver_right)
        )
        if directional:
            functional = (
                "cell_cycle_phase_complementarity", 0.84,
                "Directionally supported inhibition of replication/S-phase and mitotic/M-phase programs provides complementary cell-cycle evidence.",
            )
    elif functional and functional[0] == "damage_repair_co_modulation":
        if left & DNA_DAMAGE and right & DNA_REPAIR:
            driver_left, driver_right = left & DNA_DAMAGE, right & DNA_REPAIR
        else:
            driver_left, driver_right = left & DNA_REPAIR, right & DNA_DAMAGE
        directional = bool(
            not (driver_left & driver_right)
            and inhibitory_drivers(action_left, driver_left)
            and inhibitory_drivers(action_right, driver_right)
        )
        if directional:
            functional = (
                "replication_damage_reinforcement", 0.80,
                "Directionally supported inhibition across DNA-damage and repair programs provides a reinforcement hypothesis.",
            )
    elif functional and functional[0] == "shared_pathway_co_modulation":
        shared_pathway_genes = set().union(*(context.pathways.get(name, set()) for name in shared_pathways))
        if (
            inhibitory_drivers(action_left, left & shared_pathway_genes)
            and inhibitory_drivers(action_right, right & shared_pathway_genes)
        ):
            functional = (
                "vertical_pathway_coblockade",
                min(0.90, 0.60 + 0.05 * len(shared_pathways)),
                "Directionally supported inhibitory targets from both drugs converge on a shared pathway.",
            )
    cross_distance = _mean_cross_distance(
        context.ppi, left, right, context.target_distance
    )
    disease = context.disease_genes(disease_label)
    cheng = _cheng_metrics(
        context.ppi,
        left,
        right,
        disease,
        samples=requested_randomizations,
        random_seed=random_seed,
        min_mapping_coverage=effective_mapping_coverage,
        disease_distances=context.disease_distances(disease_label),
        distance_lookup=context.target_distance,
        degree_sampling_index=context.degree_sampling_index(),
    )
    disease_coverage_left = len(left & disease) / max(1, len(left))
    disease_coverage_right = len(right & disease) / max(1, len(right))

    metrics = {
        "target_jaccard": jaccard,
        "n_targets_drug1": len(left),
        "n_targets_drug2": len(right),
        "mean_cross_target_distance": cross_distance,
        "disease_coverage_drug1": disease_coverage_left,
        "disease_coverage_drug2": disease_coverage_right,
        "direct_disease_target_fraction_drug1": disease_coverage_left,
        "direct_disease_target_fraction_drug2": disease_coverage_right,
        "network_coverage": cheng.get("network_coverage", {}),
        "cheng": cheng,
        "functional_pattern": functional[0] if functional else "unknown",
        "functional_pattern_score": functional[1] if functional else None,
        "functional_pattern_description": functional[2] if functional else "",
        "pathway_enrichment": {
            "method": "one_sided_hypergeometric_BH_FDR",
            "fdr_threshold": 0.05,
            "pathway_size_bounds": {"minimum": 5, "maximum": 500},
            "descriptive_shared_pathway_count": len(descriptive_shared_pathways),
            "descriptive_shared_pathways": descriptive_shared_pathways,
            "significant_drug1": list(enrichment_left.values()),
            "significant_drug2": list(enrichment_right.values()),
            "significant_shared": [
                {
                    "pathway": name,
                    "drug1": enrichment_left[name],
                    "drug2": enrichment_right[name],
                    "max_q": max(
                        float(enrichment_left[name]["q"]),
                        float(enrichment_right[name]["q"]),
                    ),
                }
                for name in shared_pathways
            ],
        },
        "target_action": {"drug1": action_left, "drug2": action_right},
    }
    evidence: List[Dict[str, Any]] = []
    if shared:
        evidence.append({"type": "shared_targets", "targets": shared})
    if shared_pathways:
        evidence.append({
            "type": "shared_pathways",
            "pathways": shared_pathways,
            "statistical_test": "one_sided_hypergeometric_BH_FDR",
            "fdr_threshold": 0.05,
            "details": metrics["pathway_enrichment"]["significant_shared"],
        })
    elif descriptive_shared_pathways:
        evidence.append({
            "type": "descriptive_shared_pathways_not_enriched",
            "pathways": descriptive_shared_pathways,
            "semantics": "raw overlap only; not used to select the primary mechanism",
        })
    if complementary:
        evidence.append({
            "type": "distinct_pathway_annotation_pairs",
            "pairs": complementary,
            "semantics": "descriptive Cartesian pairing; not tested pathway complementarity",
        })

    network_map = {
        "P1": ("high_overlap", 0.70, "Both drugs occupy overlapping network neighborhoods near the disease module."),
        "P2": ("complementary_exposure", 0.90, "Both drugs are close to the disease module in distinct network neighborhoods."),
        "P4": ("single_exposure", 0.60, "Only one drug is directly close to the disease module."),
        "P4_edge": ("overlapping_single_exposure", 0.50, "Only one drug is close to the disease module, with overlap between the target neighborhoods."),
    }
    if jaccard >= 0.3 or len(shared) >= 2:
        shared_directions = [
            (
                action_left["qualified_by_target"].get(target, "unknown"),
                action_right["qualified_by_target"].get(target, "unknown"),
            )
            for target in shared
        ]
        if shared_directions and all(
            left_direction == right_direction and left_direction in {"inhibitory", "activating"}
            for left_direction, right_direction in shared_directions
        ):
            mechanism = "target_overlap_or_redundancy"
            confidence = min(0.95, 0.55 + 0.35 * jaccard + 0.03 * len(shared))
            description = "The drugs share directionally concordant contextualized targets, supporting an overlap/redundancy interpretation."
        elif any(
            {left_direction, right_direction} == {"inhibitory", "activating"}
            for left_direction, right_direction in shared_directions
        ):
            mechanism = "opposing_shared_target_modulation"
            confidence = min(0.70, 0.45 + 0.20 * jaccard)
            description = "The drugs share targets but have opposing documented action directions on at least one shared target."
        else:
            mechanism = "shared_target_co_modulation"
            confidence = min(0.68, 0.42 + 0.20 * jaccard)
            description = "The drugs share contextualized targets, but action direction is incomplete or unresolved."
    elif functional and functional[0] in {
        "cell_cycle_phase_complementarity", "replication_damage_reinforcement",
        "cell_cycle_phase_co_modulation", "damage_repair_co_modulation",
    }:
        mechanism, confidence, description = functional
        evidence.append({"type": "functional_mechanism", "pattern": mechanism})
    elif (
        cheng.get("valid")
        and cheng.get("quality_status") == "adequate"
        and cheng.get("statistically_supported")
        and cheng.get("pattern") in network_map
    ):
        mechanism, confidence, description = network_map[cheng["pattern"]]
        evidence.append({"type": "cheng_network_pattern", **cheng})
    elif shared_pathways:
        mechanism, confidence, description = functional or (
            "shared_pathway_co_modulation", 0.52,
            "Distinct targets from both drugs map to common pathways; action direction is not established.",
        )
    elif complementary and cross_distance is not None and cross_distance <= 3.0:
        mechanism = "network_close_distinct_processes"
        confidence = min(0.70, 0.45 + 0.05 * (3.0 - cross_distance) + 0.02 * len(complementary))
        description = "The target profiles have distinct pathway annotations and are close in the PPI network; pathway complementarity has not been tested."
    elif complementary:
        mechanism = "distinct_annotated_processes"
        confidence = min(0.55, 0.35 + 0.02 * len(complementary))
        description = "The target profiles have distinct pathway annotations; no pathway interaction or enrichment test was performed."
    elif cross_distance is not None and cross_distance <= 2.0:
        mechanism = "network_proximal_targets"
        confidence = min(0.72, 0.48 + 0.10 * (2.0 - cross_distance))
        description = "The target sets are distinct but close in the PPI network."
    else:
        mechanism = "unknown"
        confidence = 0.25 if (shared or shared_pathways or cross_distance is not None) else 0.0
        description = "The available target and cancer-context evidence does not support one dominant mechanism."

    warnings: List[str] = []
    if context.ppi.number_of_nodes() == 0:
        warnings.append("missing_ppi_context")
    if not disease:
        warnings.append("missing_disease_context")
    if not cheng.get("valid"):
        warnings.append("insufficient_network_support")
    elif not cheng.get("statistically_supported"):
        warnings.append("network_result_not_statistically_supported")
    warnings.extend(cheng.get("warnings", []))
    if not context.pathways:
        warnings.append("missing_pathway_context")
    if all(item.evidence_status == "predicted_only" for item in list(drug1_targets[:top_k]) + list(drug2_targets[:top_k])):
        warnings.append("mechanism_depends_only_on_external_predictions")
    warnings.append("mechanism_score_not_calibrated_probability")

    if any(
        sum(item.evidence_status == "predicted_only" for item in items) / max(1, len(items)) > 0.5
        for items in (selected1, selected2)
    ):
        warnings.append("mechanism_mostly_depends_on_external_predictions")

    def fraction(items: Sequence[ContextScore], attribute: str) -> float:
        return sum(bool(getattr(item, attribute, False)) for item in items) / max(1, len(items))

    coverage = {
        "drug1": {
            "targets": len(selected1),
            "ppi_mapping_fraction": fraction(selected1, "ppi_mapped"),
            "expression_fraction": fraction(selected1, "expression_available"),
            "expression_measured_fraction": fraction(selected1, "expression_measured"),
            "expression_detected_fraction": fraction(selected1, "expression_detected"),
            "pathway_mapping_fraction": fraction(selected1, "pathway_mapped"),
            "action_direction_fraction": action_left["direction_known_fraction"],
            "action_conflict_count": action_left["conflict_count"],
            "qualified_action_direction_fraction": action_left["qualified_direction_known_fraction"],
            "qualified_action_conflict_count": action_left["qualified_conflict_count"],
        },
        "drug2": {
            "targets": len(selected2),
            "ppi_mapping_fraction": fraction(selected2, "ppi_mapped"),
            "expression_fraction": fraction(selected2, "expression_available"),
            "expression_measured_fraction": fraction(selected2, "expression_measured"),
            "expression_detected_fraction": fraction(selected2, "expression_detected"),
            "pathway_mapping_fraction": fraction(selected2, "pathway_mapped"),
            "action_direction_fraction": action_right["direction_known_fraction"],
            "action_conflict_count": action_right["conflict_count"],
            "qualified_action_direction_fraction": action_right["qualified_direction_known_fraction"],
            "qualified_action_conflict_count": action_right["qualified_conflict_count"],
        },
        "network": cheng.get("network_coverage", {}),
        "disease_context_genes": len(disease),
    }
    if min(coverage["drug1"]["ppi_mapping_fraction"], coverage["drug2"]["ppi_mapping_fraction"]) < effective_mapping_coverage:
        warnings.append("low_ppi_target_coverage")
    if min(coverage["drug1"]["expression_fraction"], coverage["drug2"]["expression_fraction"]) < 0.5:
        warnings.append("low_target_expression_coverage")
    if min(coverage["drug1"]["action_direction_fraction"], coverage["drug2"]["action_direction_fraction"]) < 0.5:
        warnings.append("low_target_action_direction_coverage")
    if coverage["drug1"]["action_conflict_count"] or coverage["drug2"]["action_conflict_count"]:
        warnings.append("conflicting_target_action_directions")
    warnings = list(dict.fromkeys(warnings))

    return MechanismResult(
        mechanism_type=mechanism,
        confidence=float(max(0.0, min(confidence, 1.0))),
        description=description,
        shared_targets=shared,
        shared_pathways=shared_pathways,
        complementary_pathways=complementary,
        metrics=metrics,
        evidence=evidence,
        warnings=warnings,
        credibility_score=None,
        coverage=coverage,
        uncertainty=list(warnings),
    )
