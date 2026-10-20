"""Cancer-context resources and target ranking.

The module intentionally owns all expression, PPI, pathway, and disease-state
logic so the orchestrator does not need to know individual file formats.
"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from dataclasses import asdict, dataclass, field
from functools import wraps
from pathlib import Path
from threading import RLock
from types import MappingProxyType
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple, TypeVar
from copy import deepcopy
from hashlib import sha256
import json
import math
import pickle

import networkx as nx
import numpy as np
import pandas as pd

from .targets import FusedTarget


def _missing_label(value: Any) -> bool:
    """Return whether a scalar identifier is a true missing value.

    ``str(np.nan)`` is ``"nan"`` and used to silently create a fake gene or
    cell-line label.  Keep literal strings intact, but reject scalar NA/NaN
    values at resource boundaries.
    """
    if value is None:
        return True
    if isinstance(value, str):
        return False
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        return False
    return bool(missing) if isinstance(missing, (bool, np.bool_)) else False


def gene(value: Any) -> str:
    return "" if _missing_label(value) else str(value).strip().upper()


def _required_gene(value: Any, location: str) -> str:
    if _missing_label(value):
        raise ValueError(f"missing/NaN gene label in {location}")
    normalized = gene(value)
    if not normalized:
        raise ValueError(f"blank gene label in {location}")
    return normalized


def _required_label(value: Any, location: str) -> str:
    if _missing_label(value):
        raise ValueError(f"missing/NaN label in {location}")
    normalized = str(value).strip()
    if not normalized:
        raise ValueError(f"blank label in {location}")
    return normalized


_MethodReturn = TypeVar("_MethodReturn")


def _synchronized(method: Callable[..., _MethodReturn]) -> Callable[..., _MethodReturn]:
    """Serialize access to mutable resources and their derived caches."""
    @wraps(method)
    def wrapper(self: "CancerContext", *args: Any, **kwargs: Any) -> _MethodReturn:
        with self._state_lock:
            return method(self, *args, **kwargs)
    return wrapper


def _finite_float(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _hypergeometric_tail(
    hits: int,
    population_size: int,
    success_population: int,
    draws: int,
) -> float:
    """Exact one-sided P(X >= hits), evaluated in log space."""
    if hits <= 0:
        return 1.0
    upper = min(success_population, draws)
    lower = max(0, draws - (population_size - success_population))
    if hits <= lower:
        return 1.0
    if hits > upper or population_size <= 0:
        return 0.0

    def log_choose(n: int, k: int) -> float:
        if k < 0 or k > n:
            return float("-inf")
        return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

    denominator = log_choose(population_size, draws)
    log_probabilities = [
        log_choose(success_population, count)
        + log_choose(population_size - success_population, draws - count)
        - denominator
        for count in range(max(hits, lower), upper + 1)
    ]
    finite = [value for value in log_probabilities if math.isfinite(value)]
    if not finite:
        return 0.0
    maximum = max(finite)
    tail = math.exp(maximum) * sum(math.exp(value - maximum) for value in finite)
    return min(1.0, max(0.0, tail))


def _bh_fdr(p_values: Sequence[float]) -> List[float]:
    """Benjamini-Hochberg adjusted p-values in original order."""
    count = len(p_values)
    adjusted = [1.0] * count
    running = 1.0
    for rank, index in reversed(list(enumerate(
        sorted(range(count), key=lambda item: (p_values[item], item)),
        start=1,
    ))):
        candidate = min(1.0, float(p_values[index]) * count / rank)
        running = min(running, candidate)
        adjusted[index] = running
    return adjusted


def _stable_sha256(value: Any) -> str:
    payload = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), default=str)
    return sha256(payload.encode("utf-8")).hexdigest()


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _normalized_confidence(value: Any, scale: str = "auto") -> Optional[float]:
    """Normalize common PPI confidence scales to 0..1 when possible."""
    try:
        score = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(score) or score < 0:
        return None
    if scale == "string_1000" and score > 1.0:
        score /= 1000.0
    elif scale == "percent_100" and score > 1.0:
        score /= 100.0
    elif scale == "auto" and score > 100.0:
        score /= 1000.0
    elif scale == "auto" and score > 1.0:
        score /= 100.0
    return min(max(score, 0.0), 1.0)


def _edge_confidence(
    attributes: Mapping[str, Any],
    scale: str = "auto",
    weight_semantics: str = "unknown",
) -> tuple[Optional[float], bool]:
    if "confidence" in attributes:
        return _normalized_confidence(attributes.get("confidence"), scale), True
    if "combined_score" in attributes:
        return _normalized_confidence(attributes.get("combined_score"), "string_1000"), True
    if "score" in attributes:
        return _normalized_confidence(attributes.get("score"), scale), True
    if "weight" in attributes and weight_semantics == "confidence":
        return _normalized_confidence(attributes.get("weight"), scale), True
    return None, False


def edge_distance(_: Any, __: Any, attributes: Mapping[str, Any]) -> float:
    """Return a path cost from an edge while preserving explicit distances."""
    try:
        explicit = float(attributes.get("distance"))
        if math.isfinite(explicit) and explicit > 0:
            return explicit
    except (TypeError, ValueError):
        pass
    confidence, supplied = _edge_confidence(attributes)
    if confidence is None:
        confidence = 0.05
    return 1.0 / max(confidence, 0.05)


@dataclass
class ContextScore:
    target_id: str
    target_name: str
    final_score: float
    evidence_score: float
    expression_score: Optional[float]
    disease_proximity_score: float
    pathway_score: float
    centrality_score: float
    evidence_status: str
    source: str
    source_evidence: Dict[str, Any]
    expression_value: Optional[float] = None
    expression_available: bool = False
    expression_measured: bool = False
    expression_detected: bool = False
    expression_missing: bool = True
    expression_background_complete: Optional[bool] = None
    available_scoring_dimensions: List[str] = field(default_factory=list)
    unavailable_scoring_dimensions: Dict[str, str] = field(default_factory=dict)
    effective_weights: Dict[str, float] = field(default_factory=dict)
    local_evidence_tier: str = "not_applicable"
    ppi_mapped: bool = False
    pathway_mapped: bool = False
    disease_context_available: bool = False
    scoring_notes: List[str] = field(default_factory=list)
    pharmacology: Dict[str, Any] = field(default_factory=dict)
    ordering_policy: str = "context_score_descending"
    ordering_score: Optional[float] = None

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    @property
    def contribution_score(self) -> float:
        """Compatibility name used by the previous interpreter."""
        return self.final_score

    @property
    def confidence(self) -> str:
        """Compatibility alias; this is a ranking tier, not a probability."""
        return "High" if self.final_score >= 0.65 else "Medium" if self.final_score >= 0.4 else "Low"

    @property
    def ranking_tier(self) -> str:
        return self.confidence


class CancerContext:
    """In-memory cancer context with compact, explicit loaders."""

    def __init__(self):
        self._state_lock = RLock()
        self.ppi: nx.Graph = nx.Graph()
        self.expression: Dict[str, Dict[str, float]] = {}
        self.pathways: Dict[str, Set[str]] = {}
        self.pathway_sources: Dict[str, str] = {}
        self.pathway_versions: Dict[str, str] = {}
        self.pathway_taxids: Dict[str, str] = {}
        self.pathway_identifier_namespaces: Dict[str, str] = {}
        self.pathway_mapping_versions: Dict[str, str] = {}
        self.pathway_content_hashes: Dict[str, str] = {}
        self.pathway_artifact_kinds: Dict[str, str] = {}
        self.pathway_artifact_hashes: Dict[str, str] = {}
        self.disease_sets: Dict[str, Set[str]] = {}
        self.default_disease_genes: Set[str] = set()
        self.default_disease_label: str = ""
        self.resource_metadata: Dict[str, Dict[str, Any]] = {
            "ppi": {}, "expression": {}, "pathways": {}, "disease": {},
        }
        self._runtime_mutations: Set[str] = set()
        self._runtime_findings: List[str] = []
        self._last_runtime_audit: Dict[str, Any] = {
            "status": "not_run", "mutated_resources": [], "findings": [],
        }
        self._resource_provenance_baselines: Dict[str, Dict[str, Any]] = {}
        self._resource_generations: Dict[str, int] = {
            "ppi": 0, "expression": 0, "pathways": 0, "disease": 0,
        }
        self._disease_distance_cache: Dict[Tuple[int, str], Dict[str, float]] = {}
        self._target_distance_cache: Dict[Tuple[int, str, str], Optional[float]] = {}
        self._target_source_distance_cache: Dict[Tuple[int, str], Dict[str, float]] = {}
        self._gene_to_pathways_cache: Optional[Dict[str, Tuple[str, ...]]] = None
        self._pathway_gene_universe_cache: Optional[Set[str]] = None
        self._ppi_rank_metadata_cache: Optional[Dict[str, Any]] = None
        self._ppi_public_topology_cache: Optional[Dict[str, Any]] = None
        self.is_demo_context = False

    def _invalidate_caches(self, resource: str) -> None:
        self._resource_generations[resource] = self._resource_generations.get(resource, 0) + 1
        if resource == "ppi":
            self._disease_distance_cache.clear()
            self._target_distance_cache.clear()
            self._target_source_distance_cache.clear()
            self._ppi_rank_metadata_cache = None
            self._ppi_public_topology_cache = None
        elif resource == "disease":
            self._disease_distance_cache.clear()
        elif resource == "pathways":
            self._gene_to_pathways_cache = None
            self._pathway_gene_universe_cache = None

    def _resource_reloaded(self, resource: str) -> None:
        self._invalidate_caches(resource)
        # A mutation is sticky for this CancerContext instance. Loading an
        # unrelated cell line/pathway/label must not clear it; create a fresh
        # pipeline after an integrity failure.
        if resource in self._runtime_mutations:
            return
        self._resource_provenance_baselines[resource] = deepcopy(
            self.resource_metadata.get(resource, {})
        )

    def _reject_changed_artifact(self, resource: str, path: Path) -> None:
        finding = f"{resource}:artifact_changed_while_loading"
        self._runtime_mutations.add(resource)
        self._runtime_findings = list(dict.fromkeys(self._runtime_findings + [finding]))
        raise ValueError(f"{resource.capitalize()} artifact changed while it was being loaded: {path}")

    def _audit_before_resource_load(self, resource: str) -> None:
        if resource in self._resource_provenance_baselines:
            self.audit_runtime_resources()

    def _record_expression_metadata(
        self,
        lines: Iterable[str],
        source: str,
        version: str,
        background_complete: Optional[bool],
        detection_threshold: float,
        taxid: str,
        identifier_namespace: str,
        mapping_version: str,
        artifact_kind: str,
        artifact_sha256: str,
    ) -> None:
        if background_complete is not None and not isinstance(background_complete, bool):
            raise TypeError("background_complete must be True, False, or None")
        previous = self.resource_metadata.get("expression", {})
        declared = dict(previous.get("background_complete_by_cell_line", {}))
        sources = dict(previous.get("source_by_cell_line", {}))
        versions = dict(previous.get("version_by_cell_line", {}))
        thresholds = dict(previous.get("detection_threshold_by_cell_line", {}))
        taxids = dict(previous.get("taxid_by_cell_line", {}))
        namespaces = dict(previous.get("identifier_namespace_by_cell_line", {}))
        mapping_versions = dict(previous.get("mapping_version_by_cell_line", {}))
        artifact_kinds = dict(previous.get("artifact_kind_by_cell_line", {}))
        artifact_hashes = dict(previous.get("artifact_sha256_by_cell_line", {}))
        for line in lines:
            declared[str(line)] = background_complete
            sources[str(line)] = source
            versions[str(line)] = version
            thresholds[str(line)] = detection_threshold
            taxids[str(line)] = str(taxid)
            namespaces[str(line)] = str(identifier_namespace)
            mapping_versions[str(line)] = str(mapping_version)
            artifact_kinds[str(line)] = artifact_kind
            artifact_hashes[str(line)] = artifact_sha256
        self.resource_metadata["expression"] = {
            "source": source,
            "version": version,
            "cell_lines": sorted(self.expression),
            "gene_counts": {line: len(values) for line, values in self.expression.items()},
            "background_complete_by_cell_line": declared,
            "source_by_cell_line": sources,
            "version_by_cell_line": versions,
            "detection_threshold_by_cell_line": thresholds,
            "taxid_by_cell_line": taxids,
            "identifier_namespace_by_cell_line": namespaces,
            "mapping_version_by_cell_line": mapping_versions,
            "artifact_kind_by_cell_line": artifact_kinds,
            "artifact_sha256_by_cell_line": artifact_hashes,
            "content_sha256_by_cell_line": {
                line: _stable_sha256(sorted(values.items()))
                for line, values in self.expression.items()
            },
        }

    @staticmethod
    def _ppi_topology_metadata(graph: nx.Graph) -> Dict[str, Any]:
        nodes = graph.number_of_nodes()
        largest = max((len(component) for component in nx.connected_components(graph)), default=0)
        return {
            "nodes": nodes,
            "edges": graph.number_of_edges(),
            "largest_connected_component_nodes": largest,
            "largest_connected_component_fraction": largest / max(1, nodes),
            "content_sha256": _stable_sha256({
                "nodes": sorted(str(node) for node in graph.nodes()),
                "edges": sorted(
                    (
                        min(str(left), str(right)),
                        max(str(left), str(right)),
                        data.get("confidence"),
                        data.get("distance"),
                        data.get("confidence_available"),
                        data.get("confidence_warning"),
                    )
                    for left, right, data in graph.edges(data=True)
                ),
            }),
        }

    @staticmethod
    def _prepare_ppi(
        graph: nx.Graph,
        confidence_scale: str = "auto",
        weight_semantics: str = "auto",
        taxid: str = "",
        identifier_namespace: str = "gene_symbol",
        mapping_version: str = "",
    ) -> nx.Graph:
        if graph.is_multigraph():
            raise TypeError("PPI MultiGraph is not supported; collapse parallel edges explicitly")
        if graph.is_directed():
            raise TypeError("PPI graph must be undirected")
        node_mapping = {
            node: _required_gene(node, "PPI graph node")
            for node in graph.nodes()
        }
        prepared = nx.relabel_nodes(graph, node_mapping, copy=True)
        prepared.remove_edges_from(list(nx.selfloop_edges(prepared)))
        invalid_confidence_edges = 0
        assumed_unweighted_edges = 0
        distance_only_edges = 0
        ambiguous_weight_edges = 0
        for _, _, attributes in prepared.edges(data=True):
            if (
                weight_semantics == "distance"
                and "distance" not in attributes
                and _finite_float(attributes.get("weight")) is not None
                and float(attributes["weight"]) > 0
            ):
                attributes["distance"] = float(attributes["weight"])
            confidence, supplied = _edge_confidence(
                attributes,
                confidence_scale,
                "confidence" if weight_semantics == "confidence" else "unknown",
            )
            if "weight" in attributes and weight_semantics == "auto" and not any(
                key in attributes for key in ("confidence", "combined_score", "score", "distance")
            ):
                ambiguous_weight_edges += 1
            explicit_distance = _finite_float(attributes.get("distance"))
            explicit_distance = explicit_distance if explicit_distance is not None and explicit_distance > 0 else None
            if confidence is None and supplied:
                confidence = 0.05
                invalid_confidence_edges += 1
                attributes["confidence"] = confidence
                attributes["confidence_warning"] = "invalid_value_conservative_floor"
                attributes["confidence_available"] = False
                attributes["distance"] = 1.0 / confidence
            elif confidence is None and explicit_distance is not None:
                distance_only_edges += 1
                attributes["confidence"] = None
                attributes["confidence_available"] = False
                attributes["confidence_warning"] = "explicit_distance_without_confidence"
                attributes["distance"] = explicit_distance
            elif confidence is None:
                assumed_unweighted_edges += 1
                attributes["confidence"] = None
                attributes["confidence_available"] = False
                attributes["confidence_warning"] = "unweighted_topology_only"
                attributes["distance"] = 20.0
            else:
                attributes["confidence"] = confidence
                attributes["confidence_available"] = True
                if explicit_distance is None:
                    attributes["distance"] = 1.0 / max(confidence, 0.05)
        prepared.graph["confidence_audit"] = {
            "scale": confidence_scale,
            "invalid_confidence_edges": invalid_confidence_edges,
            "assumed_unweighted_edges": assumed_unweighted_edges,
            "distance_only_edges": distance_only_edges,
            "total_edges": prepared.number_of_edges(),
            "invalid_confidence_fraction": invalid_confidence_edges / max(1, prepared.number_of_edges()),
            "unweighted_edge_fraction": assumed_unweighted_edges / max(1, prepared.number_of_edges()),
            "ambiguous_weight_edges": ambiguous_weight_edges,
            "weight_semantics": weight_semantics,
        }
        return prepared

    @_synchronized
    def load_ppi(
        self,
        source: Any,
        source_name: str = "",
        version: str = "",
        confidence_scale: str = "auto",
        weight_semantics: str = "auto",
        taxid: str = "",
        identifier_namespace: str = "gene_symbol",
        mapping_version: str = "",
        allow_unsafe_pickle: bool = False,
    ) -> "CancerContext":
        self._audit_before_resource_load("ppi")
        if confidence_scale not in {"auto", "string_1000", "percent_100"}:
            raise ValueError("confidence_scale must be auto, string_1000, or percent_100")
        if weight_semantics not in {"auto", "confidence", "distance", "topology"}:
            raise ValueError("weight_semantics must be auto, confidence, distance, or topology")
        if isinstance(source, (nx.Graph, nx.DiGraph, nx.MultiGraph, nx.MultiDiGraph)):
            self.ppi = self._prepare_ppi(source, confidence_scale, weight_semantics)
            self.resource_metadata["ppi"] = {
                "source": source_name or "in_memory",
                "version": version,
                "taxid": str(taxid),
                "identifier_namespace": identifier_namespace,
                "mapping_version": mapping_version,
                "artifact_kind": "in_memory_user_declared",
                "artifact_sha256": "",
                "deserialization_safety": "in_memory_no_deserialization",
                "unsafe_deserialization": False,
                **self._ppi_topology_metadata(self.ppi),
                "confidence_audit": self.ppi.graph.get("confidence_audit", {}),
            }
            self._resource_reloaded("ppi")
            return self
        path = Path(source)
        if path.suffix.lower() in {".pkl", ".pickle"}:
            if not allow_unsafe_pickle:
                raise ValueError(
                    "PPI pickle loading is disabled because pickle can execute code. "
                    "Use a .tsv/.tsv.gz edge table, or explicitly opt in with "
                    "allow_unsafe_pickle=True only for a separately authenticated trusted artifact."
                )
            artifact_hash_before = _file_sha256(path)
            with path.open("rb") as handle:
                graph = pickle.load(handle)
            if isinstance(graph, dict) and isinstance(graph.get("graph"), nx.Graph):
                graph = graph["graph"]
            if not isinstance(graph, (nx.Graph, nx.DiGraph, nx.MultiGraph, nx.MultiDiGraph)):
                raise TypeError(f"PPI pickle does not contain a networkx graph: {path}")
            artifact_hash_after = _file_sha256(path)
            if artifact_hash_before != artifact_hash_after:
                self._reject_changed_artifact("ppi", path)
            self.ppi = self._prepare_ppi(graph, confidence_scale, weight_semantics)
            embedded = graph.graph.get("metadata", {}) if isinstance(graph.graph, dict) else {}
            self.resource_metadata["ppi"] = {
                "source": source_name or str(path.resolve()),
                "version": version or str(embedded.get("version", "")),
                "taxid": str(taxid or embedded.get("taxid", "")),
                "identifier_namespace": identifier_namespace,
                "mapping_version": mapping_version,
                "artifact_kind": "unsafe_pickle_file",
                "artifact_sha256": artifact_hash_before,
                "deserialization_safety": "unsafe_pickle_explicit_opt_in",
                "unsafe_deserialization": True,
                **self._ppi_topology_metadata(self.ppi),
                "confidence_audit": self.ppi.graph.get("confidence_audit", {}),
            }
            self._resource_reloaded("ppi")
            return self
        artifact_hash_before = _file_sha256(path)
        lower_name = path.name.lower()
        separator = "\t" if lower_name.endswith((".tsv", ".tsv.gz")) else (
            "," if lower_name.endswith((".csv", ".csv.gz")) else None
        )
        frame = pd.read_csv(
            path,
            sep=separator,
            engine="python" if separator is None else "c",
            compression="infer",
        )
        lowered = {str(column).strip().lower(): column for column in frame.columns}
        left = next((lowered[name] for name in ("gene_a", "protein_a", "source", "gene1") if name in lowered), None)
        right = next((lowered[name] for name in ("gene_b", "protein_b", "target", "gene2") if name in lowered), None)
        weight = next((lowered[name] for name in ("confidence", "weight", "score", "combined_score") if name in lowered), None)
        if not left or not right:
            raise ValueError("PPI table requires gene_a/gene_b (or compatible) columns")
        graph = nx.Graph()
        for row_number, row in enumerate(frame.to_dict(orient="records"), start=2):
            a = _required_gene(row.get(left), f"PPI table row {row_number} column {left}")
            b = _required_gene(row.get(right), f"PPI table row {row_number} column {right}")
            if a == b:
                continue
            raw_weight = row.get(weight) if weight else None
            parsed_weight = _finite_float(raw_weight)
            parse_failed = bool(weight) and parsed_weight is None
            edge_weight = parsed_weight if parsed_weight is not None else 1.0
            column_name = str(weight).lower() if weight else ""
            known_string_score = column_name == "combined_score"
            known_confidence_column = column_name in {"confidence", "combined_score"}
            use_as_confidence = known_confidence_column or (bool(weight) and weight_semantics == "confidence")
            use_as_distance = bool(weight) and weight_semantics == "distance"
            confidence = None if parse_failed or not use_as_confidence else _normalized_confidence(
                edge_weight,
                "string_1000" if known_string_score else confidence_scale,
            )
            warning = ""
            if confidence is None and weight and (parse_failed or use_as_confidence):
                confidence = 0.05
                warning = "invalid_value_conservative_floor"
                confidence_available = False
                edge_distance_value = 1.0 / confidence
            elif confidence is None:
                confidence_available = False
                warning = "explicit_distance_without_confidence" if use_as_distance and not parse_failed else (
                    "ambiguous_weight_topology_only" if weight and weight_semantics == "auto" else "unweighted_topology_only"
                )
                edge_distance_value = edge_weight if use_as_distance and not parse_failed and edge_weight > 0 else 20.0
            else:
                confidence_available = True
                edge_distance_value = 1.0 / max(confidence, 0.05)
            graph.add_edge(
                a, b,
                weight=edge_weight,
                confidence=confidence,
                confidence_available=confidence_available,
                distance=edge_distance_value,
                confidence_warning=warning,
            )
        self.ppi = graph
        artifact_hash_after = _file_sha256(path)
        if artifact_hash_before != artifact_hash_after:
            self._reject_changed_artifact("ppi", path)
        self.resource_metadata["ppi"] = {
            "source": source_name or str(path.resolve()),
            "version": version,
            "taxid": str(taxid),
            "identifier_namespace": identifier_namespace,
            "mapping_version": mapping_version,
            "artifact_kind": "safe_edge_table_file",
            "artifact_sha256": artifact_hash_before,
            "deserialization_safety": "non_executable_tabular_parser",
            "unsafe_deserialization": False,
            **self._ppi_topology_metadata(graph),
            "confidence_audit": {
                "scale": "string_1000" if weight and str(weight).lower() == "combined_score" else confidence_scale,
                "invalid_confidence_edges": sum(
                    data.get("confidence_warning") == "invalid_value_conservative_floor"
                    for _, _, data in graph.edges(data=True)
                ),
                "assumed_unweighted_edges": sum(
                    data.get("confidence_warning") in {"unweighted_topology_only", "ambiguous_weight_topology_only"}
                    for _, _, data in graph.edges(data=True)
                ),
                "distance_only_edges": 0,
                "total_edges": graph.number_of_edges(),
                "invalid_confidence_fraction": sum(
                    data.get("confidence_warning") == "invalid_value_conservative_floor"
                    for _, _, data in graph.edges(data=True)
                ) / max(1, graph.number_of_edges()),
                "unweighted_edge_fraction": sum(
                    data.get("confidence_warning") in {"unweighted_topology_only", "ambiguous_weight_topology_only"}
                    for _, _, data in graph.edges(data=True)
                ) / max(1, graph.number_of_edges()),
                "ambiguous_weight_edges": sum(
                    data.get("confidence_warning") == "ambiguous_weight_topology_only"
                    for _, _, data in graph.edges(data=True)
                ),
                "weight_semantics": weight_semantics,
            },
        }
        self._resource_reloaded("ppi")
        return self

    @_synchronized
    def load_expression(
        self,
        source: Any,
        cell_line: str = "default",
        source_name: str = "",
        version: str = "",
        background_complete: Optional[bool] = None,
        detection_threshold: float = 0.0,
        taxid: str = "",
        identifier_namespace: str = "gene_symbol",
        mapping_version: str = "",
    ) -> "CancerContext":
        self._audit_before_resource_load("expression")
        threshold = _finite_float(detection_threshold)
        if threshold is None:
            raise ValueError("detection_threshold must be finite")
        if isinstance(source, Mapping):
            if source and all(isinstance(value, Mapping) for value in source.values()):
                self.expression = {}
                for line, values in source.items():
                    line_label = _required_label(line, "in-memory expression cell line")
                    cleaned: Dict[str, float] = {}
                    for key, value in values.items():
                        gene_label = _required_gene(
                            key, f"in-memory expression profile {line_label!r}"
                        )
                        number = _finite_float(value)
                        if number is not None:
                            cleaned[gene_label] = number
                    self.expression[line_label] = cleaned
                updated_lines = list(self.expression)
            else:
                line_label = _required_label(cell_line, "in-memory expression cell line")
                cleaned = {}
                for key, value in source.items():
                    gene_label = _required_gene(
                        key, f"in-memory expression profile {line_label!r}"
                    )
                    number = _finite_float(value)
                    if number is not None:
                        cleaned[gene_label] = number
                self.expression[line_label] = cleaned
                updated_lines = [line_label]
            self._record_expression_metadata(
                updated_lines, source_name or "in_memory", version, background_complete, threshold,
                str(taxid), identifier_namespace, mapping_version,
                "in_memory_user_declared", "",
            )
            self._resource_reloaded("expression")
            return self
        expression_path = Path(source)
        expression_artifact_hash = _file_sha256(expression_path)
        frame = pd.read_csv(expression_path)
        lowered = {column.lower(): column for column in frame.columns}
        gene_col = next((lowered[name] for name in ("gene", "gene_symbol", "target", "symbol") if name in lowered), None)
        value_col = next((lowered[name] for name in ("expression", "expression_value", "value", "tpm") if name in lowered), None)
        line_col = next((lowered[name] for name in ("cell_line", "cell", "model") if name in lowered), None)
        if gene_col and value_col:
            updated_lines: Set[str] = set()
            normalized_genes = np.asarray([
                _required_gene(value, f"expression row {row_number} column {gene_col}")
                for row_number, value in enumerate(frame[gene_col].tolist(), start=2)
            ], dtype=object)
            if line_col:
                line_values = pd.Series([
                    _required_label(value, f"expression row {row_number} column {line_col}")
                    for row_number, value in enumerate(frame[line_col].tolist(), start=2)
                ], index=frame.index, dtype=object)
            else:
                requested_line = _required_label(cell_line, "expression cell_line argument")
                line_values = pd.Series(requested_line, index=frame.index, dtype=object)
            incoming_lines = set(line_values.tolist())
            for line in incoming_lines:
                self.expression[line] = {}
            numeric_values = pd.to_numeric(frame[value_col], errors="coerce").to_numpy(
                dtype=np.float64, na_value=np.nan,
            )
            normalized_lines = line_values.to_numpy(dtype=object)
            valid = (normalized_genes != "") & np.isfinite(numeric_values)
            valid_rows = pd.DataFrame({
                "line": normalized_lines[valid],
                "gene": normalized_genes[valid],
                "value": numeric_values[valid],
            })
            for line, rows in valid_rows.groupby("line", sort=False):
                line = str(line)
                # dict preserves the prior loader's last-row-wins behavior for
                # duplicate gene/cell-line observations.
                self.expression.setdefault(line, {}).update(dict(zip(
                    rows["gene"].tolist(), rows["value"].tolist(),
                )))
                updated_lines.add(line)
            self._record_expression_metadata(
                updated_lines, source_name or str(Path(source).resolve()), version, background_complete, threshold,
                str(taxid), identifier_namespace, mapping_version,
                "file", expression_artifact_hash,
            )
            if expression_artifact_hash != _file_sha256(expression_path):
                self._reject_changed_artifact("expression", expression_path)
            self._resource_reloaded("expression")
            return self
        # Wide format: first column is a gene, remaining columns are cell lines.
        gene_col = frame.columns[0]
        normalized_genes = np.asarray([
            _required_gene(value, f"expression row {row_number} column {gene_col}")
            for row_number, value in enumerate(frame[gene_col].tolist(), start=2)
        ], dtype=object)
        valid_genes = normalized_genes != ""
        numeric_values = frame.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").to_numpy(
            dtype=np.float64, na_value=np.nan,
        )
        loaded_lines: List[str] = []
        for column_index, raw_line in enumerate(frame.columns[1:]):
            line = _required_label(raw_line, "wide expression column header")
            loaded_lines.append(line)
            line_values = numeric_values[:, column_index]
            valid = valid_genes & np.isfinite(line_values)
            # Converting through lists keeps plain Python strings/floats in the
            # public expression mapping, matching the row-wise implementation.
            self.expression[line] = dict(zip(
                normalized_genes[valid].tolist(), line_values[valid].tolist(),
            ))
        self._record_expression_metadata(
            loaded_lines,
            source_name or str(Path(source).resolve()),
            version,
            background_complete,
            threshold,
            str(taxid),
            identifier_namespace,
            mapping_version,
            "file",
            expression_artifact_hash,
        )
        if expression_artifact_hash != _file_sha256(expression_path):
            self._reject_changed_artifact("expression", expression_path)
        self._resource_reloaded("expression")
        return self

    @_synchronized
    def load_pathways(
        self,
        source: Any,
        source_name: str = "",
        version: str = "",
        taxid: str = "",
        identifier_namespace: str = "gene_symbol",
        mapping_version: str = "",
    ) -> "CancerContext":
        self._audit_before_resource_load("pathways")
        if isinstance(source, Mapping):
            for raw_name, genes in source.items():
                name = _required_label(raw_name, "in-memory pathway name")
                if isinstance(genes, (str, bytes)):
                    raise TypeError(f"pathway {name!r} genes must be an iterable of gene identifiers, not a string")
                selected = {
                    _required_gene(item, f"in-memory pathway {name!r}")
                    for item in genes
                }
                self.pathways[name] = selected
                self.pathway_sources[name] = source_name or "in_memory"
                self.pathway_versions[name] = version
                self.pathway_taxids[name] = str(taxid)
                self.pathway_identifier_namespaces[name] = identifier_namespace
                self.pathway_mapping_versions[name] = mapping_version
                self.pathway_content_hashes[name] = _stable_sha256(sorted(selected))
                self.pathway_artifact_kinds[name] = "in_memory_user_declared"
                self.pathway_artifact_hashes[name] = ""
            self.resource_metadata["pathways"] = {
                "source": source_name or "in_memory",
                "version": version,
                "taxid": str(taxid),
                "identifier_namespace": identifier_namespace,
                "mapping_version": mapping_version,
                "pathway_count": len(self.pathways),
                "gene_union_count": len(set().union(*self.pathways.values())) if self.pathways else 0,
                "sources": sorted(set(self.pathway_sources.values())),
                "versions": sorted(set(self.pathway_versions.values())),
                "provenance_by_pathway": {
                    name: {
                        "source": self.pathway_sources.get(name, ""),
                        "version": self.pathway_versions.get(name, ""),
                        "taxid": self.pathway_taxids.get(name, ""),
                        "identifier_namespace": self.pathway_identifier_namespaces.get(name, ""),
                        "mapping_version": self.pathway_mapping_versions.get(name, ""),
                        "content_sha256": self.pathway_content_hashes.get(name, ""),
                        "artifact_kind": self.pathway_artifact_kinds.get(name, ""),
                        "artifact_sha256": self.pathway_artifact_hashes.get(name, ""),
                    }
                    for name in sorted(self.pathways)
                },
            }
            self._resource_reloaded("pathways")
            return self
        path = Path(source)
        if path.suffix.lower() == ".json":
            artifact_hash = _file_sha256(path)
            payload = json.loads(path.read_text(encoding="utf-8"))
            self.load_pathways(
                payload, source_name or str(path.resolve()), version,
                taxid, identifier_namespace, mapping_version,
            )
            if artifact_hash != _file_sha256(path):
                self._reject_changed_artifact("pathways", path)
            for name in payload:
                self.pathway_artifact_kinds[str(name)] = "file"
                self.pathway_artifact_hashes[str(name)] = artifact_hash
                provenance = self.resource_metadata["pathways"]["provenance_by_pathway"].get(str(name), {})
                provenance.update({"artifact_kind": "file", "artifact_sha256": artifact_hash})
            self.resource_metadata["pathways"].update({
                "artifact_kind": "file",
                "artifact_sha256": artifact_hash,
            })
            self._resource_reloaded("pathways")
            return self
        if path.suffix.lower() == ".gmt":
            artifact_hash = _file_sha256(path)
            for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
                parts = line.split("\t")
                if len(parts) >= 3:
                    name = _required_label(parts[0], f"GMT pathway name in {path}")
                    selected = {
                        _required_gene(item, f"GMT pathway {name!r} in {path}")
                        for item in parts[2:]
                    }
                    self.pathways[name] = selected
                    self.pathway_sources[name] = source_name or path.stem
                    self.pathway_versions[name] = version
                    self.pathway_taxids[name] = str(taxid)
                    self.pathway_identifier_namespaces[name] = identifier_namespace
                    self.pathway_mapping_versions[name] = mapping_version
                    self.pathway_content_hashes[name] = _stable_sha256(sorted(selected))
                    self.pathway_artifact_kinds[name] = "file"
                    self.pathway_artifact_hashes[name] = artifact_hash
            if artifact_hash != _file_sha256(path):
                self._reject_changed_artifact("pathways", path)
            self.resource_metadata["pathways"] = {
                "source": source_name or str(path.resolve()),
                "version": version,
                "taxid": str(taxid),
                "identifier_namespace": identifier_namespace,
                "mapping_version": mapping_version,
                "artifact_kind": "file",
                "artifact_sha256": artifact_hash,
                "pathway_count": len(self.pathways),
                "gene_union_count": len(set().union(*self.pathways.values())) if self.pathways else 0,
                "sources": sorted(set(self.pathway_sources.values())),
                "versions": sorted(set(self.pathway_versions.values())),
                "provenance_by_pathway": {
                    name: {
                        "source": self.pathway_sources.get(name, ""),
                        "version": self.pathway_versions.get(name, ""),
                        "taxid": self.pathway_taxids.get(name, ""),
                        "identifier_namespace": self.pathway_identifier_namespaces.get(name, ""),
                        "mapping_version": self.pathway_mapping_versions.get(name, ""),
                        "content_sha256": self.pathway_content_hashes.get(name, ""),
                        "artifact_kind": self.pathway_artifact_kinds.get(name, ""),
                        "artifact_sha256": self.pathway_artifact_hashes.get(name, ""),
                    }
                    for name in sorted(self.pathways)
                },
            }
            self._resource_reloaded("pathways")
            return self
        raise ValueError(f"Unsupported pathway format: {path}")

    @_synchronized
    def load_disease_genes(
        self,
        genes: Iterable[str],
        name: str = "default",
        source_name: str = "",
        version: str = "",
        taxid: str = "",
        identifier_namespace: str = "gene_symbol",
        mapping_version: str = "",
    ) -> "CancerContext":
        self._audit_before_resource_load("disease")
        if isinstance(genes, (str, bytes)):
            raise TypeError("disease_genes must be an iterable of gene identifiers, not a string")
        name = _required_label(name, "disease label")
        selected = {
            _required_gene(item, f"disease gene set {name!r}")
            for item in genes
        }
        self.disease_sets[name] = selected
        if (
            name == "default"
            or not self.default_disease_label
            or name == self.default_disease_label
        ):
            self.default_disease_genes = selected
            self.default_disease_label = name
        previous = self.resource_metadata.get("disease", {})
        source_by_label = dict(previous.get("source_by_label", {}))
        version_by_label = dict(previous.get("version_by_label", {}))
        taxid_by_label = dict(previous.get("taxid_by_label", {}))
        namespace_by_label = dict(previous.get("identifier_namespace_by_label", {}))
        mapping_by_label = dict(previous.get("mapping_version_by_label", {}))
        hash_by_label = dict(previous.get("content_sha256_by_label", {}))
        artifact_kind_by_label = dict(previous.get("artifact_kind_by_label", {}))
        artifact_hash_by_label = dict(previous.get("artifact_sha256_by_label", {}))
        source_by_label[name] = source_name or "user_supplied"
        version_by_label[name] = version
        taxid_by_label[name] = str(taxid)
        namespace_by_label[name] = identifier_namespace
        mapping_by_label[name] = mapping_version
        hash_by_label[name] = _stable_sha256(sorted(selected))
        artifact_kind_by_label[name] = "in_memory_user_declared"
        artifact_hash_by_label[name] = ""
        self.resource_metadata["disease"] = {
            "source": source_name or "user_supplied",
            "version": version,
            "source_by_label": source_by_label,
            "version_by_label": version_by_label,
            "taxid_by_label": taxid_by_label,
            "identifier_namespace_by_label": namespace_by_label,
            "mapping_version_by_label": mapping_by_label,
            "content_sha256_by_label": hash_by_label,
            "artifact_kind_by_label": artifact_kind_by_label,
            "artifact_sha256_by_label": artifact_hash_by_label,
            "labels": sorted(self.disease_sets),
            "gene_count": len(selected),
            "gene_counts": {label: len(values) for label, values in self.disease_sets.items()},
            "gene_union_count": len(set().union(*self.disease_sets.values())) if self.disease_sets else 0,
        }
        self._resource_reloaded("disease")
        return self

    @_synchronized
    def load_disease_genes_file(
        self,
        source: str | Path,
        name: str = "default",
        source_name: str = "",
        version: str = "",
        taxid: str = "",
        identifier_namespace: str = "gene_symbol",
        mapping_version: str = "",
        gene_column: str = "",
    ) -> "CancerContext":
        path = Path(source)
        artifact_hash = _file_sha256(path)
        if path.suffix.lower() in {".csv", ".tsv"}:
            frame = pd.read_csv(path, sep="\t" if path.suffix.lower() == ".tsv" else ",")
            lowered = {str(column).lower(): column for column in frame.columns}
            selected_column = (
                lowered.get(gene_column.lower()) if gene_column else
                next((lowered[key] for key in ("gene", "gene_symbol", "symbol", "target") if key in lowered), frame.columns[0])
            )
            if selected_column is None:
                raise ValueError(f"Disease gene column not found: {gene_column}")
            genes = frame[selected_column].tolist()
        else:
            genes = []
            for line in path.read_text(encoding="utf-8", errors="ignore").splitlines():
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                value = stripped.split("\t", 1)[0].split(",", 1)[0].strip()
                if gene(value) not in {"GENE", "GENE_SYMBOL", "SYMBOL"}:
                    genes.append(value)
        self.load_disease_genes(
            genes, name, source_name or str(path.resolve()), version, taxid,
            identifier_namespace, mapping_version,
        )
        if artifact_hash != _file_sha256(path):
            self._reject_changed_artifact("disease", path)
        metadata = self.resource_metadata["disease"]
        metadata["artifact_kind_by_label"][name] = "file"
        metadata["artifact_sha256_by_label"][name] = artifact_hash
        metadata.update({"artifact_kind": "file", "artifact_sha256": artifact_hash})
        self._resource_reloaded("disease")
        return self

    @_synchronized
    def audit_runtime_resources(self) -> Dict[str, Any]:
        """Recompute loaded-resource identities and detect direct post-load mutation."""
        current_findings: List[str] = []
        changed_resources: Set[str] = set()

        for resource, baseline in self._resource_provenance_baselines.items():
            if baseline != self.resource_metadata.get(resource, {}):
                current_findings.append(f"{resource}:metadata_changed_after_load")
                self._runtime_mutations.add(resource)
                changed_resources.add(resource)

        ppi_metadata = self.resource_metadata.get("ppi", {})
        if ppi_metadata.get("content_sha256"):
            observed_ppi = self._ppi_topology_metadata(self.ppi)
            if any(
                ppi_metadata.get(key) != observed_ppi.get(key)
                for key in ("nodes", "edges", "content_sha256")
            ):
                current_findings.append("ppi:content_changed_after_load")
                self._runtime_mutations.add("ppi")
                changed_resources.add("ppi")
            confidence_audit = dict(ppi_metadata.get("confidence_audit", {}))
            total_edges = self.ppi.number_of_edges()
            invalid_edges = sum(
                data.get("confidence_warning") == "invalid_value_conservative_floor"
                for _, _, data in self.ppi.edges(data=True)
            )
            unweighted_edges = sum(
                data.get("confidence_warning") in {"unweighted_topology_only", "ambiguous_weight_topology_only"}
                or (
                    data.get("confidence_available") is False
                    and data.get("confidence") is None
                    and data.get("confidence_warning") != "explicit_distance_without_confidence"
                )
                for _, _, data in self.ppi.edges(data=True)
            )
            distance_only_edges = sum(
                data.get("confidence_warning") == "explicit_distance_without_confidence"
                for _, _, data in self.ppi.edges(data=True)
            )
            confidence_audit.update({
                "invalid_confidence_edges": invalid_edges,
                "assumed_unweighted_edges": unweighted_edges,
                "distance_only_edges": distance_only_edges,
                "total_edges": total_edges,
                "invalid_confidence_fraction": invalid_edges / max(1, total_edges),
                "unweighted_edge_fraction": unweighted_edges / max(1, total_edges),
            })
            ppi_metadata.update(observed_ppi)
            ppi_metadata["confidence_audit"] = confidence_audit

        expression_metadata = self.resource_metadata.get("expression", {})
        expected_expression = expression_metadata.get("content_sha256_by_cell_line", {})
        observed_expression = {
            line: _stable_sha256(sorted(values.items()))
            for line, values in self.expression.items()
        }
        if expected_expression and expected_expression != observed_expression:
            current_findings.append("expression:content_changed_after_load")
            self._runtime_mutations.add("expression")
            changed_resources.add("expression")
        if expression_metadata:
            expression_metadata["content_sha256_by_cell_line"] = observed_expression
            expression_metadata["cell_lines"] = sorted(self.expression)
            expression_metadata["gene_counts"] = {
                line: len(values) for line, values in self.expression.items()
            }

        expected_pathways = dict(self.pathway_content_hashes)
        observed_pathways = {
            name: _stable_sha256(sorted(values)) for name, values in self.pathways.items()
        }
        observed_pathway_provenance = {
            name: {
                "source": self.pathway_sources.get(name, ""),
                "version": self.pathway_versions.get(name, ""),
                "taxid": self.pathway_taxids.get(name, ""),
                "identifier_namespace": self.pathway_identifier_namespaces.get(name, ""),
                "mapping_version": self.pathway_mapping_versions.get(name, ""),
                "content_sha256": observed_pathways.get(name, ""),
                "artifact_kind": self.pathway_artifact_kinds.get(name, ""),
                "artifact_sha256": self.pathway_artifact_hashes.get(name, ""),
            }
            for name in sorted(self.pathways)
        }
        baseline_pathway_provenance = self._resource_provenance_baselines.get(
            "pathways", {}
        ).get("provenance_by_pathway", {})
        if baseline_pathway_provenance and baseline_pathway_provenance != observed_pathway_provenance:
            current_findings.append("pathways:provenance_changed_after_load")
            self._runtime_mutations.add("pathways")
            changed_resources.add("pathways")
        if expected_pathways and expected_pathways != observed_pathways:
            current_findings.append("pathways:content_changed_after_load")
            self._runtime_mutations.add("pathways")
            changed_resources.add("pathways")
        self.pathway_content_hashes = observed_pathways
        pathway_metadata = self.resource_metadata.get("pathways", {})
        if pathway_metadata:
            pathway_metadata["pathway_count"] = len(self.pathways)
            pathway_metadata["gene_union_count"] = len(set().union(*self.pathways.values())) if self.pathways else 0
            pathway_metadata["provenance_by_pathway"] = observed_pathway_provenance

        disease_metadata = self.resource_metadata.get("disease", {})
        expected_disease = disease_metadata.get("content_sha256_by_label", {})
        observed_disease = {
            label: _stable_sha256(sorted(values)) for label, values in self.disease_sets.items()
        }
        default_diverged = bool(
            self.default_disease_label
            and self.default_disease_label in self.disease_sets
            and self.default_disease_genes != self.disease_sets[self.default_disease_label]
        )
        if expected_disease and (expected_disease != observed_disease or default_diverged):
            current_findings.append("disease:content_or_default_selection_changed_after_load")
            self._runtime_mutations.add("disease")
            changed_resources.add("disease")
        if disease_metadata:
            disease_metadata["content_sha256_by_label"] = observed_disease
            disease_metadata["labels"] = sorted(self.disease_sets)
            disease_metadata["gene_counts"] = {
                label: len(values) for label, values in self.disease_sets.items()
            }
            disease_metadata["gene_union_count"] = len(set().union(*self.disease_sets.values())) if self.disease_sets else 0

        for resource in changed_resources:
            self._invalidate_caches(resource)
        for resource in self._resource_provenance_baselines:
            if resource not in self._runtime_mutations:
                self._resource_provenance_baselines[resource] = deepcopy(
                    self.resource_metadata.get(resource, {})
                )
        self._runtime_findings = list(dict.fromkeys(self._runtime_findings + current_findings))
        report = {
            "status": "mutated_after_load" if self._runtime_mutations else "verified_current_snapshot",
            "mutated_resources": sorted(self._runtime_mutations),
            "findings": list(self._runtime_findings),
            "audit_semantics": "content identities recomputed before explanation; detected mutations remain sticky for this context and require a new pipeline",
        }
        self._last_runtime_audit = deepcopy(report)
        return report

    def _pathway_index(self) -> Dict[str, Tuple[str, ...]]:
        if self._gene_to_pathways_cache is None:
            mutable: Dict[str, List[str]] = {}
            for pathway_name in sorted(self.pathways):
                for symbol in self.pathways[pathway_name]:
                    mutable.setdefault(symbol, []).append(pathway_name)
            self._gene_to_pathways_cache = {
                symbol: tuple(names) for symbol, names in mutable.items()
            }
        return self._gene_to_pathways_cache

    def _pathway_gene_universe(self) -> Set[str]:
        if self._pathway_gene_universe_cache is None:
            self._pathway_gene_universe_cache = set(self._pathway_index())
        return self._pathway_gene_universe_cache

    @property
    @_synchronized
    def pathway_gene_universe_size(self) -> int:
        """Number of distinct genes in the currently loaded pathway universe."""
        return len(self._pathway_gene_universe())

    @property
    @_synchronized
    def gene_to_pathways(self) -> Mapping[str, Tuple[str, ...]]:
        """Read-only cached inverted pathway index for high-volume consumers."""
        return MappingProxyType(self._pathway_index())

    @_synchronized
    def pathway_names_for_gene(self, target: str) -> Tuple[str, ...]:
        return self._pathway_index().get(gene(target), ())

    @_synchronized
    def pathway_enrichment(
        self,
        targets: Iterable[str],
        *,
        min_size: int = 5,
        max_size: int = 500,
    ) -> List[Dict[str, Any]]:
        """One-sided pathway over-representation with BH-FDR correction.

        The query is restricted to the loaded pathway gene universe.  FDR is
        calculated across every pathway passing the size bounds; the returned
        compact list contains pathways with at least one query hit.
        """
        if isinstance(targets, (str, bytes)):
            raise TypeError("targets must be an iterable of gene identifiers, not a string")
        if isinstance(min_size, bool) or isinstance(max_size, bool):
            raise TypeError("min_size and max_size must be integers")
        min_size, max_size = int(min_size), int(max_size)
        if min_size < 1 or max_size < min_size:
            raise ValueError("pathway size bounds require 1 <= min_size <= max_size")
        selected = {
            _required_gene(item, "pathway enrichment query")
            for item in targets
        }
        universe = self._pathway_gene_universe()
        query = selected & universe
        population_size = len(universe)
        if not query or population_size == 0:
            return []

        tested: List[Tuple[str, Set[str], int, float]] = []
        for name in sorted(self.pathways):
            members = self.pathways[name] & universe
            background_size = len(members)
            if not min_size <= background_size <= max_size:
                continue
            hits = members & query
            p_value = _hypergeometric_tail(
                len(hits), population_size, background_size, len(query)
            )
            tested.append((name, hits, background_size, p_value))
        q_values = _bh_fdr([item[3] for item in tested])
        tested_count = len(tested)
        results = [
            {
                "pathway": name,
                "hits": sorted(hits),
                "hit_count": len(hits),
                "background_size": background_size,
                "query_size": len(query),
                "universe_size": population_size,
                "tested_pathway_count": tested_count,
                "p": float(p_value),
                "q": float(q_value),
            }
            for (name, hits, background_size, p_value), q_value in zip(tested, q_values)
            if hits
        ]
        results.sort(key=lambda item: (item["q"], item["p"], item["pathway"]))
        return results

    def _disease_distance_map(self, disease: Set[str]) -> Dict[str, float]:
        mapped = sorted(item for item in disease if item in self.ppi)
        if not mapped:
            return {}
        cache_key = (
            self._resource_generations["ppi"],
            _stable_sha256(mapped),
        )
        cached = self._disease_distance_cache.get(cache_key)
        if cached is None:
            cached = dict(nx.multi_source_dijkstra_path_length(
                self.ppi, mapped, weight=edge_distance,
            ))
            self._disease_distance_cache[cache_key] = cached
        return cached

    @_synchronized
    def disease_distances(self, label: str = "") -> Mapping[str, float]:
        """Read-only multi-source PPI distances from a disease module."""
        return MappingProxyType(self._disease_distance_map(self.disease_genes(label)))

    @_synchronized
    def target_distances(self, target: str) -> Mapping[str, float]:
        """Read-only weighted distances from one target, cached per source."""
        target_gene = gene(target)
        if not target_gene or target_gene not in self.ppi:
            return MappingProxyType({})
        cache_key = (self._resource_generations["ppi"], target_gene)
        cached = self._target_source_distance_cache.get(cache_key)
        if cached is None:
            cached = dict(nx.single_source_dijkstra_path_length(
                self.ppi, target_gene, weight=edge_distance,
            ))
            self._target_source_distance_cache[cache_key] = cached
        return MappingProxyType(cached)

    @_synchronized
    def target_distance(self, left: str, right: str) -> Optional[float]:
        """Cached weighted shortest-path distance for a canonical target pair."""
        left_gene = gene(left)
        right_gene = gene(right)
        if not left_gene or not right_gene or left_gene not in self.ppi or right_gene not in self.ppi:
            return None
        if left_gene == right_gene:
            return 0.0
        first, second = sorted((left_gene, right_gene))
        cache_key = (self._resource_generations["ppi"], first, second)
        if cache_key not in self._target_distance_cache:
            distance = self.target_distances(first).get(second)
            self._target_distance_cache[cache_key] = (
                float(distance) if distance is not None else None
            )
        return self._target_distance_cache[cache_key]

    def _ppi_rank_metadata(self) -> Dict[str, Any]:
        if self._ppi_rank_metadata_cache is None:
            degrees = dict(self.ppi.degree())
            degree_groups: Dict[int, List[str]] = {}
            for node, degree_value in degrees.items():
                degree_groups.setdefault(int(degree_value), []).append(str(node))
            frozen_groups = {
                degree_value: tuple(sorted(nodes))
                for degree_value, nodes in degree_groups.items()
            }
            self._ppi_rank_metadata_cache = {
                "max_degree": max(degrees.values(), default=1),
                "degrees": degrees,
                "nodes": tuple(sorted(str(node) for node in self.ppi.nodes())),
                "degree_groups": frozen_groups,
                "sorted_degrees": tuple(sorted(frozen_groups)),
            }
        return self._ppi_rank_metadata_cache

    @_synchronized
    def degree_sampling_index(self) -> Mapping[str, Any]:
        """Read-only cached node/degree buckets for degree-matched sampling."""
        metadata = self._ppi_rank_metadata()
        return MappingProxyType({
            "nodes": metadata["nodes"],
            "degree_groups": MappingProxyType(metadata["degree_groups"]),
            "sorted_degrees": metadata["sorted_degrees"],
        })

    @_synchronized
    def ppi_topology_metadata(self) -> Dict[str, Any]:
        """Cached serializable PPI topology identity for non-audit consumers."""
        if self._ppi_public_topology_cache is None:
            self._ppi_public_topology_cache = self._ppi_topology_metadata(self.ppi)
        return deepcopy(self._ppi_public_topology_cache)

    @_synchronized
    def disease_genes(self, label: str = "") -> Set[str]:
        if _missing_label(label):
            raise ValueError("disease label must not be NaN")
        if label:
            return set(self.disease_sets.get(label, set()))
        return set(self.default_disease_genes)

    @_synchronized
    def expression_for(self, cell_line: str) -> Dict[str, float]:
        _, values, _ = self.expression_profile(cell_line)
        return values

    @_synchronized
    def expression_profile(self, cell_line: str) -> tuple[str, Dict[str, float], bool]:
        """Return used label, profile, and whether a fallback was required."""
        if _missing_label(cell_line):
            raise ValueError("expression cell-line label must not be NaN")
        if cell_line in self.expression:
            return cell_line, self.expression[cell_line], False
        lowered = str(cell_line).strip().lower()
        for name, values in self.expression.items():
            if name.lower() == lowered:
                return name, values, False
        if "default" in self.expression:
            return "default", self.expression["default"], True
        return "", {}, bool(cell_line)

    @_synchronized
    def context_status(
        self,
        cell_line: str = "",
        disease_label: str = "",
        min_expression_background_genes: int = 1000,
        audit_resources: bool = True,
    ) -> Dict[str, Any]:
        if min_expression_background_genes < 1:
            raise ValueError("min_expression_background_genes must be positive")
        used_line, expression, fallback = self.expression_profile(cell_line)
        disease = self.disease_genes(disease_label)
        expression_metadata = self.resource_metadata.get("expression", {})
        declared_complete = expression_metadata.get("background_complete_by_cell_line", {}).get(used_line)
        detection_threshold = expression_metadata.get("detection_threshold_by_cell_line", {}).get(used_line, 0.0)
        expression_gene_count = len(expression)
        expression_values = np.asarray(list(expression.values()), dtype=float) if expression else np.asarray([], dtype=float)
        expression_unique_count = len(set(float(value) for value in expression_values))
        expression_std = float(np.std(expression_values)) if expression_values.size else None
        expression_iqr = (
            float(np.quantile(expression_values, 0.75) - np.quantile(expression_values, 0.25))
            if expression_values.size else None
        )
        expression_background_informative = bool(
            expression_gene_count >= min_expression_background_genes
            and expression_unique_count >= 10
            and expression_std is not None
            and expression_std > 1e-8
            and expression_iqr is not None
            and expression_iqr > 1e-8
        )
        expression_background_score_available = bool(
            expression
            and declared_complete is True
            and expression_background_informative
        )
        expression_detected_gene_count = sum(
            float(value) > float(detection_threshold)
            for value in expression.values()
            if math.isfinite(value)
        )
        background_scope = (
            "missing" if not expression else
            "declared_complete" if declared_complete is True else
            "declared_partial" if declared_complete is False else
            "unspecified"
        )
        warnings: List[str] = []
        if self.ppi.number_of_nodes() == 0:
            warnings.append("missing_ppi_context")
        if not expression:
            warnings.append("missing_expression_context")
        elif fallback:
            warnings.append("cell_line_expression_fallback")
        if expression:
            if expression_gene_count < min_expression_background_genes:
                warnings.append("expression_background_too_small_for_percentile")
            if not expression_background_informative:
                warnings.append("expression_background_low_dynamic_range")
            if declared_complete is False:
                warnings.append("expression_background_declared_partial")
            elif declared_complete is None:
                warnings.append("expression_background_completeness_unspecified")
            if not expression_background_score_available:
                warnings.append("expression_scoring_dimension_unavailable")
        if not self.pathways:
            warnings.append("missing_pathway_context")
        if not disease:
            warnings.append("missing_disease_context")
            if disease_label and disease_label not in self.disease_sets:
                warnings.append("disease_label_not_found")
        if self.is_demo_context:
            warnings.append("demo_context_not_for_production")
        runtime_integrity = (
            self.audit_runtime_resources() if audit_resources else deepcopy(self._last_runtime_audit)
        )
        warnings.extend(f"runtime_resource_mutation:{item}" for item in runtime_integrity["mutated_resources"])
        available_dimensions = ["evidence"]
        if expression_background_score_available:
            available_dimensions.append("expression")
        if disease and self.ppi.number_of_nodes() > 0:
            available_dimensions.append("disease")
        if self.pathways:
            available_dimensions.append("pathway")
        if self.ppi.number_of_nodes() > 0:
            available_dimensions.append("centrality")
        all_dimensions = {"evidence", "expression", "disease", "pathway", "centrality"}
        return {
            "requested_cell_line": cell_line,
            "used_expression_cell_line": used_line,
            "expression_fallback": fallback,
            "expression_gene_count": expression_gene_count,
            "expression_measured_gene_count": expression_gene_count,
            "expression_detected_gene_count": expression_detected_gene_count,
            "expression_background_scope": background_scope,
            "expression_background_complete": declared_complete,
            "expression_background_informative": expression_background_informative,
            "expression_background_score_available": expression_background_score_available,
            "expression_unique_value_count": expression_unique_count,
            "expression_standard_deviation": expression_std,
            "expression_iqr": expression_iqr,
            "expression_detection_threshold": detection_threshold,
            "min_expression_background_genes": min_expression_background_genes,
            "disease_label": disease_label or self.default_disease_label or "default",
            "used_disease_label": disease_label if disease_label in self.disease_sets else (
                self.default_disease_label if not disease_label and disease else ""
            ),
            "disease_gene_count": len(disease),
            "available_scoring_dimensions": available_dimensions,
            "unavailable_scoring_dimensions": sorted(all_dimensions - set(available_dimensions)),
            "is_demo": self.is_demo_context,
            "warnings": warnings,
            "runtime_integrity": runtime_integrity,
            "resources": deepcopy(self.resource_metadata),
        }

    @staticmethod
    def _local_evidence_tier(target: FusedTarget) -> str:
        local = target.local if isinstance(target.local, Mapping) else {}
        assertions = local.get("assertions", []) if isinstance(local, Mapping) else []
        if not isinstance(assertions, list):
            assertions = []
        def valid_quantitative(item: Any) -> bool:
            if not isinstance(item, Mapping):
                return False
            pchembl = _finite_float(item.get("pchembl_value"))
            return bool(
                str(item.get("source") or "").strip()
                and str(item.get("activity_type") or "").strip()
                and pchembl is not None
                and 0.0 <= pchembl <= 14.0
            )

        quantitative = any(valid_quantitative(item) for item in assertions)
        strong_quantitative = any(
            valid_quantitative(item) and float(item.get("pchembl_value")) >= 5.0
            for item in assertions if isinstance(item, Mapping)
        )
        provenance = local.get("provenance", {}) if isinstance(local.get("provenance"), Mapping) else {}
        provenance_hash = str(provenance.get("database_sha256") or "").lower()
        artifact_audited = bool(
            local.get("input_provenance_status") == "sqlite_artifact_audited"
            and len(provenance_hash) == 64
            and all(character in "0123456789abcdef" for character in provenance_hash)
            and int(provenance.get("schema_user_version") or 0) > 0
            and bool(str(provenance.get("dataset_version") or "").strip())
            and bool(str(provenance.get("mapping_version") or "").strip())
            and any(
                valid_quantitative(item)
                and bool(str(item.get("stable_record_id") or "").strip())
                and str(item.get("provenance_hash") or "").lower() == provenance_hash
                for item in assertions
            )
        )
        traceable = bool(str(local.get("source") or "").strip()) or any(
            isinstance(item, Mapping)
            and bool(str(item.get("source") or item.get("evidence") or "").strip())
            for item in assertions
        )
        if quantitative and artifact_audited:
            return "artifact_backed_quantitative" if strong_quantitative else "artifact_backed_weak_quantitative"
        if quantitative:
            return "declared_quantitative"
        if traceable:
            return "declared_nonquantitative"
        return "unverified_local_assertion"

    @classmethod
    def _evidence_prior(
        cls,
        target: FusedTarget,
        plato_applicability_threshold: Optional[float] = None,
    ) -> float:
        if target.evidence_status in {"known_and_predicted", "known_only"}:
            # An external alignment remains useful provenance, but does not
            # numerically boost a known target without a validated PLATO
            # reliability/applicability calibration for the active release.
            return {
                "artifact_backed_quantitative": 0.80,
                "artifact_backed_weak_quantitative": 0.58,
                "declared_quantitative": 0.62,
                "declared_nonquantitative": 0.58,
                "unverified_local_assertion": 0.50,
            }[cls._local_evidence_tier(target)]
        if target.plato:
            # PLATO's provider score is not assumed to be a calibrated
            # probability. Reliability contributes only an ordinal prior.
            applicability_supported = bool(
                plato_applicability_threshold is not None
                and target.plato.applicability_score is not None
                and target.plato.applicability_score >= plato_applicability_threshold
            )
            if target.plato.reliable is True:
                return 0.50 if applicability_supported else 0.35
            if target.plato.reliable is False:
                return 0.20 if applicability_supported else 0.15
            return 0.35 if applicability_supported else 0.25
        return 0.25

    def _disease_proximity(
        self,
        target: str,
        disease: Set[str],
        distances: Optional[Mapping[str, float]] = None,
    ) -> float:
        target = gene(target)
        if not target or target not in self.ppi or not disease:
            return 0.0
        distance_map = distances if distances is not None else self._disease_distance_map(disease)
        nearest = distance_map.get(target)
        return 1.0 / (1.0 + float(nearest)) if nearest is not None else 0.0

    def _pathway_relevance(self, target: str, disease: Set[str]) -> float:
        target = gene(target)
        membership_names = self._pathway_index().get(target, ())
        if not membership_names:
            return 0.0
        if not disease:
            return 0.5
        overlaps = [
            len(self.pathways[name] & disease) / max(1, len(disease))
            for name in membership_names
        ]
        return min(1.0, 0.5 + 2.5 * max(overlaps, default=0.0))

    @_synchronized
    def rank_targets(
        self,
        targets: Sequence[FusedTarget],
        cell_line: str = "",
        disease_label: str = "",
        weights: Optional[Mapping[str, float]] = None,
        plato_applicability_threshold: Optional[float] = None,
    ) -> List[ContextScore]:
        if plato_applicability_threshold is not None:
            threshold = _finite_float(plato_applicability_threshold)
            if threshold is None:
                raise ValueError("plato_applicability_threshold must be finite or None")
            plato_applicability_threshold = threshold
        allowed_weights = {"evidence", "expression", "disease", "pathway", "centrality"}
        weights = {
            key: value for key, value in dict(weights or {
                "evidence": 0.43,
                "expression": 0.24,
                "disease": 0.22,
                "pathway": 0.11,
                "centrality": 0.0,
            }).items()
            if key in allowed_weights
        }
        cleaned_weights: Dict[str, float] = {}
        for key, value in weights.items():
            number = _finite_float(value)
            if number is None:
                raise ValueError(f"context weight {key!r} must be finite")
            cleaned_weights[key] = max(number, 0.0)
        total = sum(cleaned_weights.values()) or 1.0
        weights = {key: value / total for key, value in cleaned_weights.items()}
        used_line, expression, _ = self.expression_profile(cell_line)
        expression_background = sorted(value for value in expression.values() if math.isfinite(value))
        expression_background_informative = bool(
            len(expression_background) >= 10
            and len(set(expression_background)) >= 10
            and float(np.std(expression_background)) > 1e-8
            and float(np.quantile(expression_background, 0.75) - np.quantile(expression_background, 0.25)) > 1e-8
        )
        expression_metadata = self.resource_metadata.get("expression", {})
        background_complete = expression_metadata.get(
            "background_complete_by_cell_line", {}
        ).get(used_line)
        expression_background_usable = bool(
            background_complete is True and expression_background_informative
        )
        detection_threshold = expression_metadata.get(
            "detection_threshold_by_cell_line", {}
        ).get(used_line, 0.0)
        detection_threshold = _finite_float(detection_threshold)
        detection_threshold = detection_threshold if detection_threshold is not None else 0.0
        disease = self.disease_genes(disease_label)
        disease_distances = self._disease_distance_map(disease)
        mapped_disease_available = bool(disease_distances)
        ppi_rank_metadata = self._ppi_rank_metadata()
        max_degree = int(ppi_rank_metadata["max_degree"])
        degrees: Mapping[str, int] = ppi_rank_metadata["degrees"]
        pathway_index = self._pathway_index()
        dimension_order = ("evidence", "expression", "disease", "pathway", "centrality")

        ranked: List[ContextScore] = []
        for index, target in enumerate(targets):
            target_gene = gene(target.target_id)
            evidence = self._evidence_prior(target, plato_applicability_threshold)
            local_evidence_tier = (
                self._local_evidence_tier(target)
                if target.evidence_status in {"known_and_predicted", "known_only"}
                else "external_prediction"
            )
            expression_value = expression.get(target_gene)
            expression_measured = expression_value is not None and math.isfinite(expression_value)
            expression_detected = bool(expression_measured and float(expression_value) > float(detection_threshold))
            expression_available = bool(expression_measured and expression_background_usable)
            expr = (
                (
                    bisect_left(expression_background, float(expression_value))
                    + bisect_right(expression_background, float(expression_value))
                ) / (2.0 * len(expression_background))
                if expression_available and expression_detected and expression_background else
                0.0 if expression_available else None
            )
            proximity = self._disease_proximity(target_gene, disease, disease_distances)
            pathway = self._pathway_relevance(target_gene, disease)
            centrality = (
                math.log1p(degrees[target_gene]) / math.log1p(max_degree)
                if target_gene in self.ppi and max_degree > 0 else 0.0
            )
            dimension_available = {
                "evidence": True,
                "expression": expression_available,
                "disease": bool(mapped_disease_available and target_gene in self.ppi),
                "pathway": bool(self.pathways),
                "centrality": target_gene in self.ppi,
            }
            unavailable_reasons: Dict[str, str] = {}
            if not expression_available:
                unavailable_reasons["expression"] = (
                    "expression_missing_for_target" if not expression_measured else
                    "expression_background_declared_partial" if background_complete is False else
                    "expression_background_completeness_unspecified" if background_complete is None else
                    "expression_background_uninformative"
                )
            if not dimension_available["disease"]:
                unavailable_reasons["disease"] = (
                    "disease_gene_set_missing_or_unmapped" if not mapped_disease_available
                    else "target_unmapped_in_ppi"
                )
            if not dimension_available["pathway"]:
                unavailable_reasons["pathway"] = "pathway_resource_missing"
            if not dimension_available["centrality"]:
                unavailable_reasons["centrality"] = "target_unmapped_in_ppi"

            available_weight = sum(
                weights.get(name, 0.0)
                for name in dimension_order
                if dimension_available[name]
            )
            effective_weights = {
                name: weights.get(name, 0.0) / available_weight
                for name in dimension_order
                if dimension_available[name] and weights.get(name, 0.0) > 0.0 and available_weight > 0.0
            }
            dimension_values: Dict[str, Optional[float]] = {
                "evidence": evidence,
                "expression": expr,
                "disease": proximity,
                "pathway": pathway,
                "centrality": centrality,
            }
            final = sum(
                effective_weights[name] * float(dimension_values[name])
                for name in effective_weights
                if dimension_values[name] is not None
            )
            ranked.append(ContextScore(
                target_id=target.target_id,
                target_name=target.target_name,
                final_score=float(final),
                evidence_score=float(evidence),
                expression_score=float(expr) if expr is not None else None,
                disease_proximity_score=float(proximity),
                pathway_score=float(pathway),
                centrality_score=float(centrality),
                evidence_status=target.evidence_status,
                source=target.source,
                source_evidence=target.to_dict(),
                expression_value=float(expression_value) if expression_measured else None,
                expression_available=expression_available,
                expression_measured=expression_measured,
                expression_detected=expression_detected,
                expression_missing=not expression_measured,
                expression_background_complete=background_complete,
                available_scoring_dimensions=[
                    name for name in dimension_order if dimension_available[name]
                ],
                unavailable_scoring_dimensions=unavailable_reasons,
                effective_weights=effective_weights,
                local_evidence_tier=local_evidence_tier,
                ppi_mapped=target_gene in self.ppi,
                pathway_mapped=target_gene in pathway_index,
                disease_context_available=bool(disease),
                scoring_notes=[
                    "plato_score_not_used_as_probability",
                    (
                        "PLATO_applicability_threshold_met"
                        if target.plato and plato_applicability_threshold is not None
                        and target.plato.applicability_score is not None
                        and target.plato.applicability_score >= plato_applicability_threshold
                        else "PLATO_applicability_unqualified_or_unconfigured"
                    ) if target.plato else "no_PLATO_prediction",
                    "external_alignment_not_used_as_probability_boost" if target.evidence_status == "known_and_predicted" else "single_source_evidence_prior",
                    "expression_ranked_against_complete_cell_line_profile" if expression_available else (
                        unavailable_reasons.get("expression", "expression_dimension_unavailable")
                    ),
                    "available_dimensions_renormalized" if len(effective_weights) < sum(
                        weights.get(name, 0.0) > 0.0 for name in dimension_order
                    ) else "configured_dimension_weights_retained",
                    "centrality_descriptive_only" if weights.get("centrality", 0.0) == 0.0 else "centrality_weighted",
                ],
            ))
        ranked.sort(key=lambda item: item.final_score, reverse=True)
        return ranked

    @_synchronized
    def pathways_for(self, targets: Iterable[str]) -> Dict[str, List[str]]:
        if isinstance(targets, (str, bytes)):
            raise TypeError("targets must be an iterable of gene identifiers, not a string")
        selected = {
            _required_gene(item, "pathway membership query")
            for item in targets
        }
        hits: Dict[str, List[str]] = {}
        index = self._pathway_index()
        for symbol in sorted(selected):
            for name in index.get(symbol, ()):
                hits.setdefault(name, []).append(symbol)
        return {name: hits[name] for name in sorted(hits)}
