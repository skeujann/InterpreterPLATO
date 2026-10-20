"""Bounded, deterministic evidence paths and explicit connectivity diagnostics."""

from __future__ import annotations

from collections import Counter, deque
import heapq
import math
from typing import Any, Dict, List, Mapping, Sequence

import networkx as nx

from .context import CancerContext, ContextScore, edge_distance, gene


PATH_METHOD = "hop_bounded_weighted_paths_v1"
DEFAULT_PATH_SEARCH = {"max_length": 5, "max_paths": 10, "max_cross_drug_paths": 5}


def validate_path_search(options: Mapping[str, Any]) -> Dict[str, int]:
    if not isinstance(options, Mapping):
        raise ValueError("path_search must be a mapping")
    unknown = set(options) - set(DEFAULT_PATH_SEARCH)
    if unknown:
        raise ValueError(f"Unknown path_search options: {sorted(unknown)}")
    values = {**DEFAULT_PATH_SEARCH, **options}
    for key, value in values.items():
        upper = 20 if key == "max_length" else 100
        if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value <= upper:
            raise ValueError(f"path_search.{key} must be an integer in [0, {upper}]")
    return values


def _hop_distances(graph: nx.Graph, destinations: Sequence[str]) -> Dict[str, int]:
    distances = {node: 0 for node in sorted(set(destinations)) if node in graph}
    queue = deque(distances)
    while queue:
        node = queue.popleft()
        for neighbor in graph[node]:
            if neighbor not in distances:
                distances[neighbor] = distances[node] + 1
                queue.append(neighbor)
    return distances


def bounded_weighted_paths(
    graph: nx.Graph, sources: Sequence[str], destinations: Sequence[str], max_length: int,
    *, hop_distances: Mapping[str, int] | None = None,
) -> Dict[str, List[str]]:
    """Minimize cost subject to a hop bound, using (node, hops) states.

    A single label per node is insufficient: a cheaper, longer prefix can
    exhaust the hop budget. Positive costs ensure final optimal walks are
    simple. Equal-cost results prefer fewer hops, then lexical path order.
    """
    if isinstance(max_length, bool) or not isinstance(max_length, int) or max_length < 0:
        raise ValueError("max_length must be a non-negative integer")
    hops = _hop_distances(graph, destinations) if hop_distances is None else hop_distances
    wanted = {source for source in sources if hops.get(source, max_length + 1) <= max_length}
    if not wanted:
        return {}
    labels = {(node, 0): (0.0, (node,)) for node in sorted(set(destinations)) if node in graph}
    heap = [(0.0, 0, path) for _, path in labels.values()]
    heapq.heapify(heap)
    result: Dict[str, List[str]] = {}
    while heap and wanted:
        cost, length, path = heapq.heappop(heap)
        node = path[0]
        if labels.get((node, length)) != (cost, path):
            continue
        if node in wanted:
            result[node] = list(path)
            wanted.remove(node)
        if length >= max_length or not wanted:
            continue
        for neighbor, attributes in graph[node].items():
            distance = edge_distance(node, neighbor, attributes)
            if not math.isfinite(distance) or distance <= 0:
                continue
            candidate = (cost + distance, (neighbor,) + path)
            state = (neighbor, length + 1)
            if state not in labels or candidate < labels[state]:
                labels[state] = candidate
                heapq.heappush(heap, (candidate[0], length + 1, candidate[1]))
    return result


def _profile(profile: Sequence[ContextScore]) -> Dict[str, ContextScore]:
    result = {}
    for target in sorted(profile, key=lambda item: (-item.final_score, gene(item.target_id))):
        result.setdefault(gene(target.target_id), target)
    return result


def path_record(path: List[str], graph: nx.Graph, target_score: float, max_degree: int) -> Dict[str, Any]:
    pairs = list(zip(path, path[1:]))
    attributes = [graph[left][right] for left, right in pairs]
    confidences = [item.get("confidence") if item.get("confidence_available") is True else None for item in attributes]
    factors = [min(1.0, max(0.0, float(value))) if value is not None and math.isfinite(float(value)) else 0.5 for value in confidences]
    confidence = math.exp(sum(math.log(max(value, 1e-6)) for value in factors) / len(factors)) if factors else 1.0
    hubs = [math.log1p(graph.degree(node)) / math.log1p(max_degree) for node in path[1:-1]] if max_degree > 0 else []
    penalty = sum(hubs) / len(hubs) if hubs else 0.0
    hop_factor = 1.0 / len(path)
    return {
        "source_target": path[0], "destination": path[-1], "path": path,
        "length": len(path) - 1, "bridge_genes": path[1:-1],
        "score": float(target_score * confidence * hop_factor * math.exp(-0.25 * penalty)),
        "weighted_cost": float(sum(edge_distance(a, b, data) for (a, b), data in zip(pairs, attributes))),
        "edge_confidences": confidences, "edge_support_factors": factors,
        "edge_confidence_geomean": confidence,
        "edge_confidence_coverage": sum(value is not None for value in confidences) / max(1, len(confidences)),
        "ppi_edge_count": len(pairs), "edge_confidence_applicable": bool(pairs),
        "hub_penalty": penalty,
        "score_components": {"target_score": target_score, "edge_confidence_geomean": confidence, "hop_factor": hop_factor, "hub_factor": math.exp(-0.25 * penalty)},
        "weighting_method": PATH_METHOD,
        "weight_semantics": "confidence_weighted_when_available; topology-only edges use a conservative support factor",
        "score_semantics": "uncalibrated_display_priority_not_mechanism_support",
        "path_semantics": "candidate_context_path_not_direct_pharmacological_causation",
    }


def _order(item: Mapping[str, Any]) -> tuple:
    return (-item["score"], item["length"], tuple(item["path"]), item.get("drug_index", 0))


def disease_paths(
    drug_names: Sequence[str], profiles: Sequence[Sequence[ContextScore]], context: CancerContext,
    *, max_paths: int, max_length: int, disease_label: str,
) -> tuple[List[Dict[str, Any]], Dict[str, Any]]:
    graph = context.ppi
    disease = context.disease_genes(disease_label)
    destinations = sorted(disease & set(graph))
    indexed = [_profile(profile) for profile in profiles]
    hops = _hop_distances(graph, destinations)
    paths = bounded_weighted_paths(graph, [node for profile in indexed for node in profile], destinations, max_length, hop_distances=hops)
    max_degree = max((degree for _, degree in graph.degree()), default=0)
    candidates = []
    for index, (name, profile) in enumerate(zip(drug_names, indexed), start=1):
        for node, target in profile.items():
            if node in paths:
                candidates.append({
                    **path_record(paths[node], graph, target.final_score, max_degree),
                    "drug": name, "drug_index": index, "path_kind": "target_to_disease",
                    "source_evidence_status": target.evidence_status,
                    "source_pharmacology": target.pharmacology,
                    "direct_disease_membership": len(paths[node]) == 1,
                })
    candidates.sort(key=_order)
    selected = []
    if max_paths >= 2:
        for index in (1, 2):
            first = next((item for item in candidates if item["drug_index"] == index), None)
            if first is not None:
                selected.append(first)
    for item in candidates:
        if len(selected) >= max_paths:
            break
        if item not in selected:
            selected.append(item)
    selected_keys = {(item["drug_index"], item["source_target"]) for item in selected}
    reports = []
    for index, (name, profile) in enumerate(zip(drug_names, indexed), start=1):
        records = []
        for node in profile:
            if not graph:
                status = "missing_ppi"
            elif node not in graph:
                status = "target_not_in_ppi"
            elif not disease:
                status = "missing_disease_context"
            elif not destinations:
                status = "disease_not_in_ppi"
            elif node not in hops:
                status = "disconnected_from_disease"
            elif hops[node] > max_length:
                status = "path_length_limit"
            elif node not in paths:
                status = "no_finite_weighted_path"
            elif (index, node) not in selected_keys:
                status = "display_limit"
            else:
                status = "shown"
            records.append({"target": node, "status": status, "minimum_hops": hops.get(node), "ppi_mapped": node in graph})
        reports.append({
            "drug": name, "drug_index": index, "eligible_target_count": len(profile),
            "mapped_target_count": sum(node in graph for node in profile),
            "reachable_target_count": sum(node in hops for node in profile),
            "within_length_target_count": sum(node in paths for node in profile),
            "displayed_target_count": sum(item["status"] == "shown" for item in records),
            "status_counts": dict(Counter(item["status"] for item in records)), "targets": records,
        })
    report = {
        "method": PATH_METHOD, "status": "available" if candidates else "no_disease_paths",
        "source_scope": "all_eligible_targets_before_display_top_k",
        "max_length": max_length, "max_paths": max_paths,
        "ppi_nodes": graph.number_of_nodes(), "ppi_edges": graph.number_of_edges(),
        "component_count": nx.number_connected_components(graph),
        "disease_gene_count": len(disease), "mapped_disease_gene_count": len(destinations),
        "unmapped_disease_genes": sorted(disease - set(graph)),
        "candidate_path_count": len(candidates), "displayed_path_count": len(selected),
        "drugs": reports,
        "semantics": "descriptive_reachability_not_statistical_or_causal_support",
    }
    return selected, report


def cross_drug_paths(
    drug_names: Sequence[str], profiles: Sequence[Sequence[ContextScore]], context: CancerContext,
    *, max_paths: int, max_length: int,
) -> tuple[List[Dict[str, Any]], Dict[str, Any]]:
    graph = context.ppi
    left, right = [_profile(profile) for profile in profiles]
    left_hops = _hop_distances(graph, list(right))
    right_hops = _hop_distances(graph, list(left))
    forward = bounded_weighted_paths(graph, list(left), list(right), max_length, hop_distances=left_hops)
    backward = bounded_weighted_paths(graph, list(right), list(left), max_length, hop_distances=right_hops)
    routes = {tuple(path) for path in forward.values()} | {tuple(reversed(path)) for path in backward.values()}
    max_degree = max((degree for _, degree in graph.degree()), default=0)
    candidates = []
    for route in sorted(routes):
        source, destination = left[route[0]], right[route[-1]]
        candidates.append({
            **path_record(list(route), graph, math.sqrt(max(0.0, source.final_score * destination.final_score)), max_degree),
            "drug": drug_names[0], "destination_drug": drug_names[1],
            "path_kind": "cross_drug_target", "shared_target": len(route) == 1,
            "source_pharmacology": source.pharmacology, "destination_pharmacology": destination.pharmacology,
            "source_evidence_status": source.evidence_status, "destination_evidence_status": destination.evidence_status,
        })
    candidates.sort(key=_order)
    if not graph:
        status = "missing_ppi"
    elif not set(left) & set(graph) or not set(right) & set(graph):
        status = "target_mapping_incomplete"
    elif not any(node in left_hops for node in left):
        status = "disconnected_target_sets"
    elif not candidates:
        status = "path_length_limit"
    else:
        status = "connected"
    return candidates[:max_paths], {
        "status": status, "max_length": max_length, "max_paths": max_paths,
        "search_scope": "nearest_opposite_target_from_each_eligible_target_in_both_directions",
        "shared_targets": sorted(set(left) & set(right) & set(graph)),
        "candidate_path_count": len(candidates), "displayed_path_count": min(max_paths, len(candidates)),
        "reachable_drug1_targets": sorted(set(left) & set(left_hops)),
        "reachable_drug2_targets": sorted(set(right) & set(right_hops)),
        "unmapped_drug1_targets": sorted(set(left) - set(graph)),
        "unmapped_drug2_targets": sorted(set(right) - set(graph)),
        "semantics": "PPI_reachability_not_synergy_or_causal_evidence",
    }
