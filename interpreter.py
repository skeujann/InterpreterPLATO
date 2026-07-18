"""Public orchestration API for the compact PLATO-aware interpreter."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from copy import deepcopy
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple
import hashlib
import json
import math

import pandas as pd

from .advanced import AdvancedResult, explain_advanced
from .context import CancerContext, ContextScore
from .mechanism import MechanismResult, classify_mechanism
from .targets import (
    CACHE_SCHEMA_VERSION,
    FusedTarget,
    PlatoPredictionCache,
    ChemblTargetMapper,
    TargetResolver,
    VERIFIED_QUERY_BINDING_STATUSES,
    standardize_smiles,
)


@dataclass
class DrugTargetExplanation:
    drug_name: str
    drug_smiles: str
    all_targets: List[ContextScore]
    top_targets: List[ContextScore]
    resolution: Dict[str, Any]

    @property
    def high_confidence_targets(self) -> List[ContextScore]:
        return [target for target in self.all_targets if target.confidence == "High"]

    def to_dict(self) -> Dict[str, Any]:
        return {
            "drug_name": self.drug_name,
            "drug_smiles": self.drug_smiles,
            "targets": [target.to_dict() for target in self.all_targets],
            "top_targets": [target.to_dict() for target in self.top_targets],
            "resolution": self.resolution,
        }


@dataclass
class ModuleHit:
    module_id: str
    pathway_name: str
    module_support_score: float
    pathway_source: str = ""
    pathway_version: str = ""
    pathway_taxid: str = ""
    pathway_identifier_namespace: str = ""
    pathway_mapping_version: str = ""
    pathway_content_sha256: str = ""
    pathway_artifact_kind: str = ""
    pathway_artifact_sha256: str = ""
    module_type: str = "pathway"
    num_genes: int = 0
    pathway_gene_count: int = 0
    matched_genes_drug1: List[str] = field(default_factory=list)
    matched_genes_drug2: List[str] = field(default_factory=list)
    drug1_support_score: Optional[float] = None
    drug2_support_score: Optional[float] = None
    drug1_enrichment_p: Optional[float] = None
    drug1_enrichment_q: Optional[float] = None
    drug2_enrichment_p: Optional[float] = None
    drug2_enrichment_q: Optional[float] = None
    enrichment_supported: bool = False
    shared_enrichment_supported: bool = False

    def to_dict(self) -> Dict[str, Any]:
        return {
            **asdict(self),
            "score_semantics": "descriptive_target_support_mean_ranked_after_pathway_enrichment",
            "statistical_enrichment_test": "one_sided_hypergeometric",
            "multiple_testing_method": "Benjamini-Hochberg",
            "fdr_threshold": 0.05,
        }


@dataclass
class ModuleExplanation:
    top_modules: List[ModuleHit] = field(default_factory=list)
    top_disease_modules: List[ModuleHit] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        return {
            "top_modules": [item.to_dict() for item in self.top_modules],
            "top_disease_modules": [item.to_dict() for item in self.top_disease_modules],
        }


@dataclass
class CombinationExplanation:
    drug1_name: str
    drug2_name: str
    cell_line: str
    predicted_synergy: float
    drug1_targets: DrugTargetExplanation
    drug2_targets: DrugTargetExplanation
    modules: ModuleExplanation
    mechanism: MechanismResult
    advanced: AdvancedResult
    metadata: Dict[str, Any] = field(default_factory=dict)
    credibility: Dict[str, Any] = field(default_factory=dict)
    model_fidelity: Dict[str, Any] = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())

    def warning_flags(self) -> List[str]:
        flags = list(self.mechanism.warnings)
        for profile in (self.drug1_targets, self.drug2_targets):
            cache_status = profile.resolution.get("cache", {}).get("status")
            if cache_status == "cache_miss":
                flags.append(f"PLATO_cache_miss:{profile.drug_name}")
            if cache_status == "cache_error":
                flags.append(f"PLATO_cache_error:{profile.drug_name}")
            if profile.resolution.get("release_mismatch"):
                flags.append(f"PLATO_release_mismatch:{profile.drug_name}")
            local_status = profile.resolution.get("local", {}).get("status")
            if local_status in {"database_error", "database_missing", "drug_identity_conflict"}:
                flags.append(f"local_target_{local_status}:{profile.drug_name}")
            for warning in profile.resolution.get("cache", {}).get("parser_warnings", []):
                flags.append(f"PLATO:{profile.drug_name}:{warning}")
            if all(target.evidence_status == "predicted_only" for target in profile.all_targets):
                flags.append(f"external_predictions_only:{profile.drug_name}")
        flags.extend(self.metadata.get("context_status", {}).get("warnings", []))
        flags.extend(self.credibility.get("warnings", []))
        if self.model_fidelity.get("status") != "evaluated":
            flags.append("model_fidelity_not_evaluated")
        return list(dict.fromkeys(flags))

    def _summary(self) -> str:
        left = ", ".join(target.target_id for target in self.drug1_targets.top_targets[:3]) or "none"
        right = ", ".join(target.target_id for target in self.drug2_targets.top_targets[:3]) or "none"
        primary_reason = ", ".join(self.credibility.get("limitations", [])[-2:])
        credibility_score = self.credibility.get("evidence_credibility_score")
        credibility_text = "invalid" if credibility_score is None else f"{float(credibility_score):.1f}/100"
        return (
            f"{self.drug1_name} ({left}) + {self.drug2_name} ({right}): "
            f"{self.mechanism.mechanism_type}, mechanism_score={self.mechanism.confidence:.3f} "
            f"(uncalibrated), evidence_credibility_score={credibility_text}. "
            f"{self.mechanism.description} Main limitations: {primary_reason or 'none reported'}."
        )

    @staticmethod
    def _selected_metadata_value(
        metadata: Mapping[str, Any],
        mapping_key: str,
        selected_key: str,
        fallback_key: str = "",
    ) -> Any:
        values = metadata.get(mapping_key, {})
        if isinstance(values, Mapping) and selected_key in values:
            return deepcopy(values[selected_key])
        return deepcopy(metadata.get(fallback_key)) if fallback_key else None

    def _compact_context_status(self) -> Dict[str, Any]:
        """Project the runtime audit onto resources used by this explanation.

        ``CancerContext.context_status`` intentionally contains a complete
        in-memory audit snapshot.  Keeping that snapshot on the result object
        is useful for credibility decisions, but serializing all cell lines
        and all pathway records into every drug-pair JSON is redundant.  This
        projection leaves the internal snapshot untouched and emits only the
        selected expression/disease records, the PPI summary, and provenance
        for pathways that actually occur in the reported module hits.
        """
        raw = self.metadata.get("context_status", {})
        if not isinstance(raw, Mapping):
            return {}
        compact = {
            str(key): deepcopy(value)
            for key, value in raw.items()
            if key != "resources"
        }
        resources = raw.get("resources", {})
        if not isinstance(resources, Mapping):
            compact["resources"] = {}
            compact["artifact_sha256"] = {
                "ppi": "", "expression": "", "pathways": "", "disease": "",
            }
            return compact

        ppi_raw = resources.get("ppi", {})
        ppi = deepcopy(dict(ppi_raw)) if isinstance(ppi_raw, Mapping) else {}

        used_line = str(raw.get("used_expression_cell_line", ""))
        expression_raw = resources.get("expression", {})
        expression_raw = expression_raw if isinstance(expression_raw, Mapping) else {}
        expression = {
            "used_cell_line": used_line,
            "source": self._selected_metadata_value(
                expression_raw, "source_by_cell_line", used_line, "source"
            ),
            "version": self._selected_metadata_value(
                expression_raw, "version_by_cell_line", used_line, "version"
            ),
            "taxid": self._selected_metadata_value(
                expression_raw, "taxid_by_cell_line", used_line, "taxid"
            ),
            "identifier_namespace": self._selected_metadata_value(
                expression_raw, "identifier_namespace_by_cell_line", used_line,
                "identifier_namespace",
            ),
            "mapping_version": self._selected_metadata_value(
                expression_raw, "mapping_version_by_cell_line", used_line,
                "mapping_version",
            ),
            "gene_count": self._selected_metadata_value(
                expression_raw, "gene_counts", used_line
            ),
            "background_complete": self._selected_metadata_value(
                expression_raw, "background_complete_by_cell_line", used_line
            ),
            "detection_threshold": self._selected_metadata_value(
                expression_raw, "detection_threshold_by_cell_line", used_line
            ),
            "content_sha256": self._selected_metadata_value(
                expression_raw, "content_sha256_by_cell_line", used_line
            ),
            "artifact_kind": self._selected_metadata_value(
                expression_raw, "artifact_kind_by_cell_line", used_line,
                "artifact_kind",
            ),
            "artifact_sha256": self._selected_metadata_value(
                expression_raw, "artifact_sha256_by_cell_line", used_line,
                "artifact_sha256",
            ),
        }

        used_disease = str(raw.get("used_disease_label", ""))
        disease_raw = resources.get("disease", {})
        disease_raw = disease_raw if isinstance(disease_raw, Mapping) else {}
        disease = {
            "used_label": used_disease,
            "source": self._selected_metadata_value(
                disease_raw, "source_by_label", used_disease, "source"
            ),
            "version": self._selected_metadata_value(
                disease_raw, "version_by_label", used_disease, "version"
            ),
            "taxid": self._selected_metadata_value(
                disease_raw, "taxid_by_label", used_disease, "taxid"
            ),
            "identifier_namespace": self._selected_metadata_value(
                disease_raw, "identifier_namespace_by_label", used_disease,
                "identifier_namespace",
            ),
            "mapping_version": self._selected_metadata_value(
                disease_raw, "mapping_version_by_label", used_disease,
                "mapping_version",
            ),
            "gene_count": self._selected_metadata_value(
                disease_raw, "gene_counts", used_disease
            ),
            "content_sha256": self._selected_metadata_value(
                disease_raw, "content_sha256_by_label", used_disease
            ),
            "artifact_kind": self._selected_metadata_value(
                disease_raw, "artifact_kind_by_label", used_disease,
                "artifact_kind",
            ),
            "artifact_sha256": self._selected_metadata_value(
                disease_raw, "artifact_sha256_by_label", used_disease,
                "artifact_sha256",
            ),
        }

        pathways_raw = resources.get("pathways", {})
        pathways_raw = pathways_raw if isinstance(pathways_raw, Mapping) else {}
        pathway_provenance = pathways_raw.get("provenance_by_pathway", {})
        pathway_provenance = (
            pathway_provenance if isinstance(pathway_provenance, Mapping) else {}
        )
        hit_items = list(self.modules.top_modules) + list(self.modules.top_disease_modules)
        hit_names = list(dict.fromkeys(
            str(item.pathway_name) for item in hit_items if str(item.pathway_name)
        ))
        hit_by_name = {str(item.pathway_name): item for item in hit_items}
        provenance_by_hit: Dict[str, Dict[str, Any]] = {}
        for name in hit_names:
            stored = pathway_provenance.get(name)
            if isinstance(stored, Mapping):
                provenance_by_hit[name] = deepcopy(dict(stored))
                continue
            item = hit_by_name[name]
            provenance_by_hit[name] = {
                "source": item.pathway_source,
                "version": item.pathway_version,
                "taxid": item.pathway_taxid,
                "identifier_namespace": item.pathway_identifier_namespace,
                "mapping_version": item.pathway_mapping_version,
                "content_sha256": item.pathway_content_sha256,
                "artifact_kind": item.pathway_artifact_kind,
                "artifact_sha256": item.pathway_artifact_sha256,
            }
        pathway_summary_keys = (
            "source", "version", "taxid", "identifier_namespace",
            "mapping_version", "artifact_kind", "artifact_sha256",
            "pathway_count", "gene_union_count",
        )
        pathways = {
            key: deepcopy(pathways_raw[key])
            for key in pathway_summary_keys if key in pathways_raw
        }
        pathways.update({
            "hit_pathway_count": len(hit_names),
            "hit_pathways": hit_names,
            "provenance_by_hit_pathway": provenance_by_hit,
        })

        compact["resources"] = {
            "ppi": ppi,
            "expression": expression,
            "pathways": pathways,
            "disease": disease,
        }
        compact["artifact_sha256"] = {
            "ppi": str(ppi.get("artifact_sha256") or ""),
            "expression": str(expression.get("artifact_sha256") or ""),
            "pathways": str(pathways.get("artifact_sha256") or ""),
            "disease": str(disease.get("artifact_sha256") or ""),
        }
        return compact

    def to_dict(self) -> Dict[str, Any]:
        compact_context_status = self._compact_context_status()
        serialized_metadata = deepcopy(self.metadata)
        if "context_status" in serialized_metadata:
            serialized_metadata["context_status"] = compact_context_status
        serialized_credibility = deepcopy(self.credibility)
        coverage = serialized_credibility.get("coverage")
        if isinstance(coverage, dict) and "context" in coverage:
            coverage["context"] = deepcopy(compact_context_status)
        dynamics = self.advanced.dynamic_hypotheses
        dynamic_type = dynamics[0].get("state_type", "unknown_dynamic_state") if dynamics else "unknown_dynamic_state"
        dynamic_confidence = dynamics[0].get("confidence", 0.0) if dynamics else 0.0
        cheng = self.mechanism.metrics.get("cheng", {})
        functional_type = self.mechanism.metrics.get("functional_pattern", "unknown")
        network_type_map = {
            "P1": "high_overlap",
            "P2": "complementary_exposure",
            "P4": "single_exposure",
            "P4_edge": "overlapping_single_exposure",
        }
        network_supported = bool(
            cheng.get("valid")
            and cheng.get("quality_status") == "adequate"
            and cheng.get("statistically_supported")
        )
        network_type = network_type_map.get(cheng.get("pattern"), "unknown") if network_supported else "insufficient_network_support"
        network_score = self.mechanism.confidence if self.mechanism.mechanism_type == network_type else 0.0
        functional_score = self.mechanism.metrics.get("functional_pattern_score")
        mechanism_score_tier = "high" if self.mechanism.confidence >= 0.75 else "medium" if self.mechanism.confidence >= 0.45 else "low"
        return {
            "result_schema_version": "5.1",
            "drug1_name": self.drug1_name,
            "drug2_name": self.drug2_name,
            "cell_line": self.cell_line,
            "predicted_synergy": self.predicted_synergy,
            "timestamp": self.timestamp,
            "drug1_targets": self.drug1_targets.to_dict(),
            "drug2_targets": self.drug2_targets.to_dict(),
            "top_targets_drug1": [(target.target_name, target.final_score) for target in self.drug1_targets.top_targets],
            "top_targets_drug2": [(target.target_name, target.final_score) for target in self.drug2_targets.top_targets],
            "top_targets_drug1_detail": [target.to_dict() for target in self.drug1_targets.top_targets],
            "top_targets_drug2_detail": [target.to_dict() for target in self.drug2_targets.top_targets],
            "target_source": self.metadata.get("target_source", "mixed"),
            "target_source_detail": {
                "drug1": self.metadata.get("target_source_drug1", ""),
                "drug2": self.metadata.get("target_source_drug2", ""),
            },
            "target_resolution_detail": {
                "drug1": self.drug1_targets.resolution,
                "drug2": self.drug2_targets.resolution,
            },
            "network_pattern": cheng.get("pattern", "insufficient_network_support"),
            "network_pair_interaction_supported": bool(cheng.get("pair_interaction_supported", False)),
            "mechanism_type": self.mechanism.mechanism_type,
            "mechanism_score": self.mechanism.confidence,
            "mechanism_score_semantics": self.mechanism.score_semantics,
            "calibrated_probability": self.mechanism.calibrated_probability,
            "evidence_credibility_score": serialized_credibility.get("evidence_credibility_score"),
            "credibility_status": serialized_credibility.get("status", "insufficient"),
            "validation": {
                "integrity_valid": bool(serialized_credibility.get("integrity_valid", False)),
                "mechanism_supported": bool(serialized_credibility.get("mechanism_supported", False)),
                "externally_validated": bool(serialized_credibility.get("externally_validated", False)),
            },
            "mechanism_confidence": self.mechanism.confidence,
            "network_mechanism_type": network_type,
            "network_mechanism_confidence": network_score,
            "network_mechanism_score_semantics": "uncalibrated_evidence_strength",
            "network_calibrated_probability": None,
            "mechanism": self.mechanism.to_dict(),
            "functional_mechanism_type": functional_type,
            "functional_mechanism_confidence": float(functional_score) if functional_score is not None else 0.0,
            "functional_mechanism": {
                "mechanism_type": functional_type,
                "confidence": float(functional_score) if functional_score is not None else 0.0,
                "score_semantics": "uncalibrated_rule_strength",
                "calibrated_probability": None,
            },
            "dynamic_state_type": dynamic_type,
            "dynamic_state_confidence": dynamic_confidence,
            "dynamic_state_evidence": self.advanced.dynamic_hypotheses,
            "local_evidence": {"top_paths": self.advanced.local_paths},
            "top_evidence_paths": self.advanced.local_paths,
            "top_evidence_edges": [edge for edge in self.advanced.evidence_graph.get("edges", []) if edge.get("type") == "PPI"],
            "evidence_graphs": {"explanation_evidence_graph": self.advanced.evidence_graph},
            "counterfactual_evidence": self.advanced.counterfactuals,
            "sensitivity_evidence": self.advanced.counterfactuals,
            "top_modules": [(item.pathway_name, item.module_support_score) for item in self.modules.top_modules],
            "top_modules_detail": [item.to_dict() for item in self.modules.top_modules],
            "top_disease_modules": [(item.pathway_name, item.module_support_score) for item in self.modules.top_disease_modules],
            "top_disease_modules_detail": [item.to_dict() for item in self.modules.top_disease_modules],
            "module_proximities": {item.module_id: item.module_support_score for item in self.modules.top_modules},
            "disease_module_proximities": {item.module_id: item.module_support_score for item in self.modules.top_disease_modules},
            "disease_context_kind": "disease_gene_set" if self.metadata.get("disease_context_size", 0) else "none",
            "disease_context_source": self.metadata.get("disease_context_source", "loaded"),
            "disease_context_labels": self.metadata.get("disease_context_labels", []),
            "disease_context_size": self.metadata.get("disease_context_size", 0),
            "mechanism_score_tier": mechanism_score_tier,
            "credibility": serialized_credibility,
            "model_fidelity": self.model_fidelity,
            "interpretation_type": self.model_fidelity.get("interpretation_type", "post_hoc_biological_hypothesis"),
            "integrated_explanation": {
                "primary_layer": "target_and_context_evidence",
                "primary_mechanism_type": self.mechanism.mechanism_type,
                "primary_confidence": self.mechanism.confidence,
                "primary_mechanism_score": self.mechanism.confidence,
                "score_semantics": self.mechanism.score_semantics,
                "calibrated_probability": None,
                "evidence_credibility_score": serialized_credibility.get("evidence_credibility_score"),
                "credibility_status": serialized_credibility.get("status", "insufficient"),
                "summary": self._summary(),
            },
            "warning_flags": self.warning_flags(),
            "metadata": serialized_metadata,
        }

    def print_summary(self) -> None:
        print(self._summary())
        if self.warning_flags():
            print("Warnings:", ", ".join(self.warning_flags()))


class SynergyInterpreter:
    """Small public facade around target, context, mechanism, and advanced layers."""

    _PRODUCT_TARGET_MODE = "plato_only"

    def __init__(self, config: Optional[Mapping[str, Any]] = None, config_path: Optional[str] = None):
        package = Path(__file__).resolve().parent
        loaded: Dict[str, Any] = {}
        if isinstance(config, (str, Path)) and not config_path:
            config_path = str(config)
            config = None
        if config_path:
            path = Path(config_path)
            text = path.read_text(encoding="utf-8")
            if path.suffix.lower() == ".json":
                loaded = json.loads(text)
            else:
                try:
                    import yaml
                except ImportError as exc:  # pragma: no cover
                    raise ImportError("pyyaml is required for YAML configuration") from exc
                loaded = yaml.safe_load(text) or {}
        explicit_config = self._merge_config(loaded, dict(config or {}))
        self.config = self._merge_config(self.default_config(package), explicit_config)
        target_mode = str(self.config.get("target_mode", "plato_only")).strip().lower()
        if target_mode != self._PRODUCT_TARGET_MODE:
            raise ValueError(
                "SynergyInterpreter is a strict PLATO-only product API; "
                "hybrid/local target modes are not supported"
            )
        if self.config.get("database") not in (None, ""):
            raise ValueError("A local target database is forbidden in strict PLATO-only mode")
        # The operational mode is a class invariant.  ``config`` remains a
        # mutable mapping because production context loading records hashes in
        # it, so every public scoring call also re-checks that it was not
        # mutated into a less strict state.
        self.config["target_mode"] = self._PRODUCT_TARGET_MODE
        self.config["database"] = None
        plato = self.config["plato"]
        if not plato.get("enabled", True):
            raise ValueError("PLATO cannot be disabled when target_mode='plato_only'")
        if plato.get("human_only") is not True:
            raise ValueError("plato.human_only must be true in strict PLATO-only mode")
        if plato.get("reliable_only") is not True:
            raise ValueError("plato.reliable_only must be true in strict PLATO-only mode")
        max_targets = int(plato.get("max_targets", 10000))
        if max_targets < 100:
            raise ValueError("plato.max_targets must be at least 100 in strict PLATO-only mode")
        self.targets = TargetResolver(
            database=None,
            plato_cache=plato.get("cache_dir"),
            human_only=True,
            reliable_only=True,
            max_plato_targets=max_targets,
            provider_release=plato.get("provider_release", "r37"),
            target_mapping=plato.get("target_mapping_path"),
            target_mapping_sha256=str(plato.get("target_mapping_sha256", "")),
        )
        self.context = CancerContext()
        self.combination_response_rows: List[Dict[str, Any]] = []
        self.prediction_evaluator: Optional[Callable[[Mapping[str, Any]], Any]] = None
        self._batch_context_audit_active = False
        self._production_context_loaded = False

    @staticmethod
    def default_config(package: Path) -> Dict[str, Any]:
        return {
            "target_mode": "plato_only",
            "database": None,
            "default_cell_line": "default",
            "top_k_targets": 10,
            "advanced": True,
            "context_weights": {
                "evidence": 0.43,
                "expression": 0.24,
                "disease": 0.22,
                "pathway": 0.11,
                "centrality": 0.0,
            },
            "network": {
                "randomizations": 1000,
                "sensitivity_randomizations": 1000,
                "random_seed": 42,
                "min_mapping_coverage": 0.5,
            },
            "credibility": {
                "top_k_stability": [5, 10, 20],
                # Stability claims must use the same minimum Monte-Carlo
                # quality budget as the primary network analysis.  Debug
                # callers may lower this value, but the result is then marked
                # descriptive rather than stable.
                "stability_randomizations": 1000,
                "max_model_perturbations": 10,
                "higher_synergy_is_stronger": True,
                "min_model_perturbations": 3,
                "min_relative_prediction_change": 0.05,
                "baseline_tolerance_fraction": 0.01,
                "max_top_k_mechanism_score_delta": 0.15,
                "min_expression_background_genes": 1000,
                "full_expression_background_genes": 5000,
                "trusted_artifact_sha256": {},
            },
            "plato": {
                "enabled": True,
                "cache_dir": str(package / "data" / "processed" / "plato_cache"),
                "provider_release": "r37",
                "human_only": True,
                "reliable_only": True,
                # r37 currently advertises fewer than this many targets.  The
                # high ceiling avoids silently taking only a top-20 subset;
                # the runtime gate fails if an export ever exceeds it.
                "max_targets": 10000,
                # Independent trust anchor.  Import-time/self-recorded hashes
                # are not sufficient: production exports must be explicitly
                # allowlisted here after independent acquisition/verification.
                "trusted_export_sha256": [],
                "target_mapping_path": str(package / "data" / "processed" / "context" / "chembl_r37_target_hgnc_mapping.tsv.gz"),
                "target_mapping_sha256": "",
                "score_min": 0.0,
                "score_max": 13.0,
                "allowed_reliable_values": ["yes", "true", "1", "reliable", "high"],
                "allowed_target_types": [
                    "single protein", "protein complex", "protein family",
                    "protein protein interaction",
                ],
                "min_applicability_for_corroboration": None,
                "min_applicability_for_prediction": None,
            },
        }

    @staticmethod
    def _merge_config(*items: Mapping[str, Any]) -> Dict[str, Any]:
        merged: Dict[str, Any] = {}
        for item in items:
            for key, value in item.items():
                if isinstance(value, Mapping) and isinstance(merged.get(key), Mapping):
                    merged[key] = {**merged[key], **value}
                else:
                    merged[key] = value
        return merged

    def load_demo_resources(self) -> Dict[str, Any]:
        processed = Path(__file__).resolve().parent / "data" / "processed"
        ppi = processed / "demo_ppi_network.pkl"
        expression = processed / "demo_expression.csv"
        pathways = processed / "demo_pathways.json"
        if ppi.exists():
            # The legacy bundled demo is never accepted as production input.
            # Keep the unsafe pickle opt-in explicit until the tiny demo asset
            # is migrated; production loading rejects pickle unconditionally.
            self.context.load_ppi(
                ppi,
                source_name="built_in_demo",
                allow_unsafe_pickle=True,
            )
        if expression.exists():
            self.context.load_expression(expression, source_name="built_in_demo")
        if pathways.exists():
            self.context.load_pathways(pathways, "demo")
        default_disease = {"TP53", "EGFR", "KRAS", "PIK3CA", "AKT1", "MAPK1", "BRAF", "MTOR"}
        self.context.load_disease_genes(default_disease, source_name="built_in_demo")
        self.context.is_demo_context = True
        self._production_context_loaded = False
        return {
            "ppi_nodes": self.context.ppi.number_of_nodes(),
            "ppi_edges": self.context.ppi.number_of_edges(),
            "cell_lines": len(self.context.expression),
            "pathways": len(self.context.pathways),
            "disease_genes": len(self.context.default_disease_genes),
        }

    @staticmethod
    def _production_file_sha256(path: Path) -> str:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    def load_production_resources(
        self,
        config_or_path: Mapping[str, Any] | str | Path,
    ) -> Dict[str, Any]:
        """Load a complete, hash-pinned production cancer context.

        File paths in a JSON/YAML configuration are resolved relative to that
        configuration file.  All four artifacts are validated before any is
        loaded, and a failed load restores the previous context.
        """
        config_path: Optional[Path] = None
        if isinstance(config_or_path, (str, Path)):
            config_path = Path(config_or_path).expanduser().resolve()
            if not config_path.is_file():
                raise FileNotFoundError(f"Production context config not found: {config_path}")
            text = config_path.read_text(encoding="utf-8")
            if config_path.suffix.lower() == ".json":
                payload = json.loads(text)
            else:
                try:
                    import yaml
                except ImportError as exc:  # pragma: no cover
                    raise ImportError("pyyaml is required for YAML production context configuration") from exc
                payload = yaml.safe_load(text)
            base_dir = config_path.parent
        elif isinstance(config_or_path, Mapping):
            payload = dict(config_or_path)
            declared_base = payload.get("base_dir")
            base_dir = (
                Path(str(declared_base)).expanduser().resolve()
                if declared_base
                else Path.cwd().resolve()
            )
        else:
            raise TypeError("config_or_path must be a mapping or a JSON/YAML path")

        if not isinstance(payload, Mapping):
            raise TypeError("Production context configuration root must be a mapping")
        resources = payload.get("resources", payload)
        if not isinstance(resources, Mapping):
            raise TypeError("Production context 'resources' must be a mapping")

        required_resources = ("ppi", "expression", "pathways", "disease")
        missing_resources = [name for name in required_resources if name not in resources]
        if missing_resources:
            raise ValueError(
                "Production context requires ppi, expression, pathways, and disease resources; "
                f"missing: {', '.join(missing_resources)}"
            )
        required_fields = (
            "path", "source", "version", "taxid",
            "identifier_namespace", "mapping_version", "sha256",
        )
        prepared: Dict[str, Dict[str, Any]] = {}
        for resource_name in required_resources:
            raw = resources[resource_name]
            if not isinstance(raw, Mapping):
                raise TypeError(f"Production resource '{resource_name}' must be a mapping")
            missing_fields = [
                field_name for field_name in required_fields
                if field_name not in raw or not str(raw[field_name]).strip()
            ]
            if missing_fields:
                raise ValueError(
                    f"Production resource '{resource_name}' is missing required fields: "
                    f"{', '.join(missing_fields)}"
                )
            expected_hash = str(raw["sha256"]).strip().lower()
            if len(expected_hash) != 64 or any(char not in "0123456789abcdef" for char in expected_hash):
                raise ValueError(
                    f"Production resource '{resource_name}' sha256 must be exactly 64 hexadecimal characters"
                )
            taxid = str(raw["taxid"]).strip()
            if taxid != "9606":
                raise ValueError(
                    f"Production resource '{resource_name}' must be human taxid 9606, got {taxid!r}"
                )
            namespace = str(raw["identifier_namespace"]).strip().lower()
            if namespace not in {"hgnc_symbol", "gene_symbol"}:
                raise ValueError(
                    f"Production resource '{resource_name}' must use an HGNC gene-symbol namespace; "
                    f"got {raw['identifier_namespace']!r}"
                )
            artifact_path = Path(str(raw["path"])).expanduser()
            if not artifact_path.is_absolute():
                artifact_path = base_dir / artifact_path
            artifact_path = artifact_path.resolve()
            if not artifact_path.is_file():
                raise FileNotFoundError(
                    f"Production resource '{resource_name}' file not found: {artifact_path}"
                )
            if resource_name == "ppi" and artifact_path.suffix.lower() in {".pkl", ".pickle"}:
                raise ValueError(
                    "Production PPI resources must use a safe tabular format; "
                    "pickle deserialization is forbidden"
                )
            actual_hash = self._production_file_sha256(artifact_path)
            if actual_hash != expected_hash:
                raise ValueError(
                    f"Production resource '{resource_name}' SHA-256 mismatch for {artifact_path}: "
                    f"expected {expected_hash}, observed {actual_hash}"
                )
            prepared[resource_name] = {
                **dict(raw),
                "path": artifact_path,
                "taxid": taxid,
                "sha256": expected_hash,
            }
        mapping_versions = {
            str(prepared[name]["mapping_version"]).strip()
            for name in required_resources
        }
        if len(mapping_versions) != 1:
            raise ValueError(
                "Production context resources must use the same HGNC mapping release; "
                f"got {sorted(mapping_versions)}"
            )
        context_mapping_version = next(iter(mapping_versions))
        if "target_mapping" in resources:
            resource_name = "target_mapping"
            raw = resources[resource_name]
            if not isinstance(raw, Mapping):
                raise TypeError("Production resource 'target_mapping' must be a mapping")
            missing_fields = [
                field_name for field_name in required_fields
                if field_name not in raw or not str(raw[field_name]).strip()
            ]
            if missing_fields:
                raise ValueError(
                    "Production resource 'target_mapping' is missing required fields: "
                    f"{', '.join(missing_fields)}"
                )
            expected_hash = str(raw["sha256"]).strip().lower()
            if len(expected_hash) != 64 or any(char not in "0123456789abcdef" for char in expected_hash):
                raise ValueError(
                    "Production resource 'target_mapping' sha256 must be exactly 64 hexadecimal characters"
                )
            taxid = str(raw["taxid"]).strip()
            if taxid != "9606":
                raise ValueError(
                    f"Production resource 'target_mapping' must be human taxid 9606, got {taxid!r}"
                )
            namespace = str(raw["identifier_namespace"]).strip().lower()
            if namespace != "chembl_target_id_to_hgnc_symbol":
                raise ValueError(
                    "Production target_mapping must translate ChEMBL target IDs to HGNC symbols"
                )
            if context_mapping_version.lower() not in str(raw["mapping_version"]).strip().lower():
                raise ValueError(
                    "Production target_mapping and context resources must declare the same HGNC release"
                )
            artifact_path = Path(str(raw["path"])).expanduser()
            if not artifact_path.is_absolute():
                artifact_path = base_dir / artifact_path
            artifact_path = artifact_path.resolve()
            if not artifact_path.is_file():
                raise FileNotFoundError(
                    f"Production resource 'target_mapping' file not found: {artifact_path}"
                )
            actual_hash = self._production_file_sha256(artifact_path)
            if actual_hash != expected_hash:
                raise ValueError(
                    f"Production resource 'target_mapping' SHA-256 mismatch for {artifact_path}: "
                    f"expected {expected_hash}, observed {actual_hash}"
                )
            prepared[resource_name] = {
                **dict(raw),
                "path": artifact_path,
                "taxid": taxid,
                "sha256": expected_hash,
            }

        runtime_manifest_verified = False
        runtime_manifest_path: Optional[Path] = None
        runtime_manifest_sha256 = ""
        if config_path is not None:
            configured_manifest = payload.get("runtime_manifest", "runtime.lock.json")
            runtime_manifest_path = Path(str(configured_manifest)).expanduser()
            if not runtime_manifest_path.is_absolute():
                runtime_manifest_path = config_path.parent / runtime_manifest_path
            runtime_manifest_path = runtime_manifest_path.resolve()
            if not runtime_manifest_path.is_file():
                raise FileNotFoundError(
                    "A version-controlled runtime integrity manifest is required for production: "
                    f"{runtime_manifest_path}"
                )
            runtime_manifest_sha256 = self._production_file_sha256(runtime_manifest_path)
            manifest_payload = json.loads(runtime_manifest_path.read_text(encoding="utf-8"))
            if not isinstance(manifest_payload, Mapping) or str(
                manifest_payload.get("manifest_kind", manifest_payload.get("lock_kind", ""))
            ) != "runtime_integrity_lock":
                raise ValueError("runtime manifest must declare manifest_kind=runtime_integrity_lock")

            manifest_root_value = str(manifest_payload.get("root", "..")).strip() or ".."
            manifest_root = (runtime_manifest_path.parent / manifest_root_value).resolve()
            manifest_records: Dict[Path, Dict[str, Any]] = {}

            def collect_manifest_records(value: Any) -> None:
                if isinstance(value, Mapping):
                    relative = value.get("relative_path")
                    digest = str(value.get("sha256", "")).strip().lower()
                    if relative and len(digest) == 64:
                        resolved = (manifest_root / str(relative)).resolve()
                        previous = manifest_records.get(resolved)
                        if previous and str(previous.get("sha256", "")).lower() != digest:
                            raise ValueError(f"runtime manifest has conflicting hashes for {resolved}")
                        manifest_records[resolved] = dict(value)
                    for nested in value.values():
                        collect_manifest_records(nested)
                elif isinstance(value, list):
                    for nested in value:
                        collect_manifest_records(nested)

            collect_manifest_records(manifest_payload)
            required_runtime_paths = {config_path.resolve()} | {
                item["path"].resolve() for item in prepared.values()
            }
            missing_from_manifest = sorted(
                str(path) for path in required_runtime_paths if path not in manifest_records
            )
            if missing_from_manifest:
                raise ValueError(
                    "runtime manifest does not cover all production inputs: "
                    + ", ".join(missing_from_manifest)
                )
            for path in required_runtime_paths:
                record = manifest_records[path]
                expected = str(record.get("sha256", "")).strip().lower()
                actual = self._production_file_sha256(path)
                if actual != expected:
                    raise ValueError(
                        f"runtime manifest SHA-256 mismatch for {path}: expected {expected}, observed {actual}"
                    )
                expected_size = record.get("size_bytes")
                if expected_size is not None and int(expected_size) != path.stat().st_size:
                    raise ValueError(f"runtime manifest size mismatch for {path}")
            runtime_manifest_verified = True

        previous_context = self.context
        previous_target_mapper = self.targets.target_mapper
        previous_production_state = self._production_context_loaded
        self.context = CancerContext()
        try:
            item = prepared["ppi"]
            self.load_ppi_network(
                ppi_file=item["path"],
                source_name=str(item["source"]),
                version=str(item["version"]),
                confidence_scale=str(item.get("confidence_scale", "auto")),
                weight_semantics=str(item.get("weight_semantics", "auto")),
                taxid=item["taxid"],
                identifier_namespace=str(item["identifier_namespace"]),
                mapping_version=str(item["mapping_version"]),
            )

            item = prepared["expression"]
            self.load_gene_expression(
                expression_file=item["path"],
                cell_line=str(item.get("cell_line", "default")),
                source_name=str(item["source"]),
                version=str(item["version"]),
                background_complete=item.get("background_complete"),
                detection_threshold=float(item.get("detection_threshold", 0.0)),
                taxid=item["taxid"],
                identifier_namespace=str(item["identifier_namespace"]),
                mapping_version=str(item["mapping_version"]),
            )

            item = prepared["pathways"]
            self.load_pathway_data(
                pathway_file=item["path"],
                source_name=str(item["source"]),
                version=str(item["version"]),
                taxid=item["taxid"],
                identifier_namespace=str(item["identifier_namespace"]),
                mapping_version=str(item["mapping_version"]),
            )

            item = prepared["disease"]
            disease_label = str(item.get("name", "default"))
            self.load_disease_genes_file(
                item["path"],
                name=disease_label,
                source_name=str(item["source"]),
                version=str(item["version"]),
                taxid=item["taxid"],
                identifier_namespace=str(item["identifier_namespace"]),
                mapping_version=str(item["mapping_version"]),
                gene_column=str(item.get("gene_column", "")),
            )
            if "target_mapping" in prepared:
                item = prepared["target_mapping"]
                self.targets.target_mapper = ChemblTargetMapper(
                    item["path"],
                    expected_sha256=item["sha256"],
                    source=str(item["source"]),
                    release=str(item["version"]),
                )
                if self.targets.target_mapper.status != "loaded":
                    raise ValueError(
                        "Production target_mapping resource did not load: "
                        f"{self.targets.target_mapper.status}; "
                        f"{self.targets.target_mapper.load_warnings}"
                    )

            observed_hashes = {
                "ppi": {str(self.context.resource_metadata.get("ppi", {}).get("artifact_sha256", ""))},
                "expression": {
                    str(value) for value in self.context.resource_metadata.get("expression", {})
                    .get("artifact_sha256_by_cell_line", {}).values()
                },
                "pathways": {
                    str(value) for value in self.context.pathway_artifact_hashes.values()
                },
                "disease": {
                    str(self.context.resource_metadata.get("disease", {})
                        .get("artifact_sha256_by_label", {}).get(disease_label, ""))
                },
            }
            for resource_name, hashes in observed_hashes.items():
                hashes.discard("")
                expected = prepared[resource_name]["sha256"]
                if not hashes or hashes != {expected}:
                    raise ValueError(
                        f"Production resource '{resource_name}' did not retain the verified artifact identity; "
                        f"expected {expected}, loaded {sorted(hashes)}"
                    )
            if self.context.ppi.number_of_nodes() == 0 or self.context.ppi.number_of_edges() == 0:
                raise ValueError("Production PPI resource loaded an empty network")
            if not self.context.expression or not any(self.context.expression.values()):
                raise ValueError("Production expression resource loaded no gene-expression values")
            if not self.context.pathways:
                raise ValueError("Production pathway resource loaded no pathways")
            if not self.context.disease_genes(disease_label):
                raise ValueError(f"Production disease resource loaded no genes for label {disease_label!r}")
            self.context.is_demo_context = False
        except Exception:
            self.context = previous_context
            self.targets.target_mapper = previous_target_mapper
            self._production_context_loaded = previous_production_state
            raise

        credibility = self.config.setdefault("credibility", {})
        if not isinstance(credibility, dict):
            self.context = previous_context
            self.targets.target_mapper = previous_target_mapper
            self._production_context_loaded = previous_production_state
            raise TypeError("credibility configuration must be a mutable mapping")
        trusted = credibility.get("trusted_artifact_sha256", {})
        if not isinstance(trusted, Mapping):
            self.context = previous_context
            self.targets.target_mapper = previous_target_mapper
            self._production_context_loaded = previous_production_state
            raise TypeError("credibility.trusted_artifact_sha256 must be a mapping")
        observed_artifacts = {
            **{name: prepared[name]["sha256"] for name in required_resources},
            **({"target_mapping": prepared["target_mapping"]["sha256"]} if "target_mapping" in prepared else {}),
        }
        credibility["verified_artifact_sha256"] = observed_artifacts
        credibility["runtime_manifest"] = {
            "verified": runtime_manifest_verified,
            "path": str(runtime_manifest_path) if runtime_manifest_path else "",
            "sha256": runtime_manifest_sha256,
            "trust_semantics": (
                "version_controlled_runtime_lock"
                if runtime_manifest_verified else "integrity_only_not_independently_trusted"
            ),
        }
        credibility["trusted_artifact_sha256"] = {
            **dict(trusted),
            **(observed_artifacts if runtime_manifest_verified else {}),
        }
        self._production_context_loaded = True

        return {
            "mode": "production",
            "config_path": str(config_path) if config_path else "",
            "runtime_manifest": credibility["runtime_manifest"],
            "resources": {
                name: {
                    "path": str(prepared[name]["path"]),
                    "source": str(prepared[name]["source"]),
                    "version": str(prepared[name]["version"]),
                    "sha256": prepared[name]["sha256"],
                }
                for name in required_resources
            } | ({
                "target_mapping": {
                    "path": str(prepared["target_mapping"]["path"]),
                    "source": str(prepared["target_mapping"]["source"]),
                    "version": str(prepared["target_mapping"]["version"]),
                    "sha256": prepared["target_mapping"]["sha256"],
                }
            } if "target_mapping" in prepared else {}),
            "ppi_nodes": self.context.ppi.number_of_nodes(),
            "ppi_edges": self.context.ppi.number_of_edges(),
            "cell_lines": len(self.context.expression),
            "pathways": len(self.context.pathways),
            "disease_label": disease_label,
            "disease_genes": len(self.context.disease_genes(disease_label)),
        }

    def load_ppi_network(
        self,
        ppi_network: Any = None,
        ppi_file: Any = None,
        source_name: str = "",
        version: str = "",
        confidence_scale: str = "auto",
        weight_semantics: str = "auto",
        taxid: str = "",
        identifier_namespace: str = "gene_symbol",
        mapping_version: str = "",
    ) -> None:
        self.context.load_ppi(
            ppi_network if ppi_network is not None else ppi_file,
            source_name=source_name,
            version=version,
            confidence_scale=confidence_scale,
            weight_semantics=weight_semantics,
            taxid=taxid,
            identifier_namespace=identifier_namespace,
            mapping_version=mapping_version,
        )

    def load_gene_expression(
        self,
        gene_expression: Any = None,
        expression_file: Any = None,
        cell_line: str = "default",
        expression_data: Any = None,
        source_name: str = "",
        version: str = "",
        background_complete: Optional[bool] = None,
        detection_threshold: float = 0.0,
        taxid: str = "",
        identifier_namespace: str = "gene_symbol",
        mapping_version: str = "",
    ) -> None:
        source = gene_expression if gene_expression is not None else expression_data
        self.context.load_expression(
            source if source is not None else expression_file,
            cell_line,
            source_name=source_name,
            version=version,
            background_complete=background_complete,
            detection_threshold=detection_threshold,
            taxid=taxid,
            identifier_namespace=identifier_namespace,
            mapping_version=mapping_version,
        )

    def load_pathway_data(
        self,
        pathway_data: Any = None,
        pathway_file: Any = None,
        source_name: str = "",
        version: str = "",
        taxid: str = "",
        identifier_namespace: str = "gene_symbol",
        mapping_version: str = "",
    ) -> None:
        self.context.load_pathways(
            pathway_data if pathway_data is not None else pathway_file,
            source_name=source_name,
            version=version,
            taxid=taxid,
            identifier_namespace=identifier_namespace,
            mapping_version=mapping_version,
        )

    def load_disease_genes(
        self,
        disease_genes: Iterable[str],
        name: str = "default",
        source_name: str = "",
        version: str = "",
        taxid: str = "",
        identifier_namespace: str = "gene_symbol",
        mapping_version: str = "",
    ) -> None:
        self.context.load_disease_genes(
            disease_genes,
            name,
            source_name=source_name,
            version=version,
            taxid=taxid,
            identifier_namespace=identifier_namespace,
            mapping_version=mapping_version,
        )

    def load_disease_genes_file(
        self,
        disease_file: str | Path,
        name: str = "default",
        source_name: str = "",
        version: str = "",
        taxid: str = "",
        identifier_namespace: str = "gene_symbol",
        mapping_version: str = "",
        gene_column: str = "",
    ) -> None:
        self.context.load_disease_genes_file(
            disease_file, name, source_name, version, taxid,
            identifier_namespace, mapping_version, gene_column,
        )

    def register_local_targets(self, drug_name: str, targets: Sequence[Any]) -> None:
        del drug_name, targets
        raise ValueError("Local target registration is not available in strict PLATO-only mode")

    def _assert_strict_plato_only_runtime(self) -> None:
        """Detect post-construction attempts to weaken the product invariant."""
        if str(self.config.get("target_mode", "")).strip().lower() != self._PRODUCT_TARGET_MODE:
            raise RuntimeError("target_mode was mutated; strict PLATO-only scoring is locked")
        if self.config.get("database") not in (None, ""):
            raise RuntimeError("A local target database was attached after initialization")
        plato = self.config.get("plato")
        if not isinstance(plato, Mapping):
            raise RuntimeError("PLATO configuration is missing or invalid")
        if (
            plato.get("enabled") is not True
            or plato.get("human_only") is not True
            or plato.get("reliable_only") is not True
        ):
            raise RuntimeError(
                "PLATO enabled/human_only/reliable_only safeguards were weakened after initialization"
            )
        configured_max = int(plato.get("max_targets", 0) or 0)
        if configured_max < 100:
            raise RuntimeError("plato.max_targets was reduced below the strict minimum of 100")
        if self.targets.database is not None:
            raise RuntimeError("TargetResolver has a forbidden local database")
        if self.targets.local_map:
            raise RuntimeError("TargetResolver contains forbidden registered local targets")
        if (
            self.targets.cache is None
            or self.targets.human_only is not True
            or self.targets.reliable_only is not True
            or int(self.targets.max_plato_targets) < 100
        ):
            raise RuntimeError("TargetResolver strict PLATO-only safeguards were weakened")
        if int(self.targets.max_plato_targets) != configured_max:
            raise RuntimeError("Configured and active PLATO target limits differ")

    @staticmethod
    def _configured_trusted_export_hashes(value: Any) -> set[str]:
        """Normalize independent export trust anchors without trusting cache metadata."""
        if isinstance(value, str):
            candidates: Iterable[Any] = [value]
        elif isinstance(value, Mapping):
            candidates = value.values()
        else:
            candidates = value or []
        normalized: set[str] = set()
        for candidate in candidates:
            nested = candidate if isinstance(candidate, (list, tuple, set)) else [candidate]
            for item in nested:
                digest = str(item or "").strip().lower()
                if len(digest) == 64 and all(char in "0123456789abcdef" for char in digest):
                    normalized.add(digest)
        return normalized

    @staticmethod
    def _plato_cache_integrity_failures(
        resolution: Mapping[str, Any],
        expected_release: str,
        trusted_export_sha256: Any = (),
        expected_structure: Optional[Mapping[str, Any]] = None,
    ) -> List[str]:
        """Return fail-closed PLATO cache issues that make pair scoring unsafe."""
        cache = resolution.get("cache", {})
        if not isinstance(cache, Mapping):
            return ["cache_metadata_missing"]
        failures: List[str] = []
        status = str(cache.get("status", ""))
        if not status.startswith("cache_hit"):
            failures.append(status or "cache_not_loaded")
        if str(cache.get("cache_schema_version", "")) != CACHE_SCHEMA_VERSION:
            failures.append("cache_schema_not_current")
        if cache.get("schema_fixture_verified") is not True:
            failures.append("PLATO_r37_schema_fixture_not_verified")
        adapter = str(cache.get("adapter", "")).strip().lower()
        if not adapter or "generic" in adapter:
            failures.append("generic_or_missing_PLATO_adapter")
        if (
            str(cache.get("provider_name", "")).strip().upper() != "PLATO"
            or cache.get("provider_validation_status") != "valid_plato"
        ):
            failures.append("PLATO_provider_not_schema_verified")
        if (
            str(cache.get("provider_method", "")).strip().lower().replace("-", "_")
            != "target_fishing"
            or cache.get("method_validation_status") != "valid_target_fishing"
        ):
            failures.append("PLATO_target_fishing_method_not_schema_verified")
        if cache.get("raw_blob_integrity_status") != "verified_immutable_raw_blob":
            failures.append("PLATO_raw_export_blob_not_immutable")
        query = cache.get("query", {})
        if (
            not isinstance(query, Mapping)
            or query.get("valid") is not True
            or query.get("standardization_status") != "rdkit_canonical_identity"
            or query.get("query_binding_status") not in VERIFIED_QUERY_BINDING_STATUSES
            or not str(query.get("canonical_isomeric_smiles", "")).strip()
            or not str(query.get("inchikey", "")).strip()
        ):
            failures.append("payload_query_structure_not_verified")
        elif expected_structure is not None and (
            str(query.get("canonical_isomeric_smiles"))
            != str(expected_structure.get("canonical_isomeric_smiles"))
            or str(query.get("inchikey")).upper() != str(expected_structure.get("inchikey")).upper()
        ):
            failures.append("payload_query_structure_identity_mismatch")
        source = cache.get("source", {})
        required_hashes = ("sha256", "payload_sha256", "normalized_targets_sha256", "audit_metadata_sha256")
        if not isinstance(source, Mapping) or any(not str(source.get(key, "")).strip() for key in required_hashes):
            failures.append("cache_audit_hashes_incomplete")
        source_sha256 = str(source.get("sha256", "")).strip().lower() if isinstance(source, Mapping) else ""
        trusted = SynergyInterpreter._configured_trusted_export_hashes(trusted_export_sha256)
        if not source_sha256 or source_sha256 not in trusted:
            failures.append("PLATO_export_sha256_not_independently_trusted")
        actual_release = str(resolution.get("provider_release", "")).strip()
        release_source = str(cache.get("release_source", "")).strip()
        release_source_verified = release_source == "payload" or (
            release_source == "verified_schema_fixture"
            and cache.get("schema_fixture_verified") is True
        )
        if (
            not actual_release
            or actual_release.lower() == "unknown"
            or actual_release.lower() != str(expected_release).strip().lower()
            or not release_source_verified
            or resolution.get("release_mismatch") is True
        ):
            failures.append("PLATO_release_not_payload_verified")
        raw_warnings = cache.get("parser_warnings", [])
        if raw_warnings is None:
            raw_warnings = []
        if not isinstance(raw_warnings, (list, tuple, set)):
            failures.append("cache_parser_warnings_invalid")
            raw_warnings = [raw_warnings]
        warnings = [str(item) for item in raw_warnings]
        integrity_names = {
            "cached_payload_hash_mismatch",
            "cached_payload_query_structure_mismatch",
            "cached_query_binding_status_mismatch",
            "normalized_targets_hash_mismatch",
            "normalized_targets_hash_missing",
            "parser_output_changed_since_import",
            "cache_query_structure_mismatch",
            "audit_metadata_hash_missing",
            "audit_metadata_hash_mismatch",
            "cache_parser_version_differs_from_runtime",
            "plato_query_binding_user_asserted_unverified",
        }
        if any(
            warning in integrity_names
            or warning.startswith("cached_parser_report_")
            or warning.startswith("cached_provider_")
            for warning in warnings
        ):
            failures.append("cache_or_parser_integrity_failure")
        try:
            parsed_targets = int(cache.get("parsed_targets", 0) or 0)
        except (TypeError, ValueError, OverflowError):
            parsed_targets = 0
            failures.append("cache_parsed_target_count_invalid")
        if parsed_targets < 1:
            failures.append("PLATO_export_contains_no_targets")
        return list(dict.fromkeys(failures))

    def _validate_plato_only_resolution(
        self,
        drug_name: str,
        smiles: str,
        fused: Sequence[FusedTarget],
        resolution: Mapping[str, Any],
    ) -> None:
        if not str(smiles).strip():
            raise ValueError(
                f"PLATO-only requires the exact submitted SMILES for {drug_name}"
            )
        structure = standardize_smiles(smiles)
        if (
            not str(structure.get("standardization_status", "")).startswith("rdkit_canonical")
            or not str(structure.get("canonical_isomeric_smiles", "")).strip()
            or not str(structure.get("inchikey", "")).strip()
        ):
            raise ValueError(
                f"Strict PLATO-only requires RDKit-validated canonical SMILES and InChIKey for {drug_name}"
            )
        plato_config = self.config["plato"]
        failures = self._plato_cache_integrity_failures(
            resolution,
            plato_config.get("provider_release", "r37"),
            plato_config.get("trusted_export_sha256", []),
            structure,
        )
        parsed_targets = int(resolution.get("cache", {}).get("parsed_targets", 0) or 0)
        if parsed_targets > int(self.targets.max_plato_targets):
            failures.append("PLATO_export_exceeds_untruncated_target_limit")
        if failures:
            raise ValueError(
                f"Unable to explain {drug_name} in strict PLATO-only mode: "
                + ", ".join(failures)
            )
        if not fused:
            raise ValueError(
                f"Unable to explain {drug_name}: PLATO returned no eligible explicit-human protein targets"
            )
        if any(
            target.local is not None
            or target.plato is None
            or target.source != "PLATO"
            or target.evidence_status != "predicted_only"
            for target in fused
            ):
            raise ValueError(f"Non-PLATO target evidence reached strict mode for {drug_name}")
        score_min = self._numeric_score(plato_config.get("score_min", 0.0))
        score_max = self._numeric_score(plato_config.get("score_max", 13.0))
        if score_min is None or score_max is None or score_min > score_max:
            raise ValueError("plato.score_min/score_max must define a finite ascending interval")
        allowed_types = {
            " ".join(str(value).strip().lower().replace("_", " ").replace("-", " ").split())
            for value in plato_config.get("allowed_target_types", [])
            if str(value).strip()
        }
        allowed_reliable = {
            str(value).strip().lower() for value in plato_config.get("allowed_reliable_values", [])
            if str(value).strip()
        }
        target_failures: List[str] = []
        for target in fused:
            evidence = target.plato
            if evidence is None:
                target_failures.append(f"{target.target_id}:missing_PLATO_evidence")
                continue
            organism = " ".join(str(evidence.organism).strip().lower().split())
            taxid = str(evidence.taxid).strip()
            inferred = "organism_inferred_from_target_name" in evidence.parse_warnings
            if not (taxid == "9606" or (organism in {"human", "homo sapiens", "h sapiens"} and not inferred)):
                target_failures.append(f"{target.target_id}:human_origin_not_explicit")
            normalized_type = " ".join(
                str(evidence.target_type).strip().lower().replace("_", " ").replace("-", " ").split()
            )
            if not allowed_types or normalized_type not in allowed_types or not evidence.is_protein_target():
                target_failures.append(f"{target.target_id}:target_type_not_allowed_protein")
            score = self._numeric_score(evidence.score)
            if score is None or not float(score_min) <= score <= float(score_max):
                target_failures.append(f"{target.target_id}:PLATO_score_out_of_range_or_missing")
            reliable_raw = str(evidence.reliable_raw).strip().lower()
            if evidence.reliable is not True or reliable_raw not in allowed_reliable:
                target_failures.append(f"{target.target_id}:PLATO_reliability_not_explicitly_reliable")
        if target_failures:
            raise ValueError(
                f"Unable to explain {drug_name}: invalid PLATO target fields: "
                + ", ".join(target_failures[:20])
                + (f" (+{len(target_failures) - 20} more)" if len(target_failures) > 20 else "")
            )

    def _filter_plato_applicability(
        self,
        drug_name: str,
        fused: Sequence[FusedTarget],
        resolution: Dict[str, Any],
    ) -> List[FusedTarget]:
        """Apply a configured applicability-domain threshold before mechanism ranking."""
        threshold_raw = self.config["plato"].get("min_applicability_for_prediction")
        if threshold_raw is None:
            resolution["mechanism_applicability_filter"] = {
                "status": "not_configured",
                "input_targets": len(fused),
                "eligible_targets": len(fused),
                "excluded_targets": 0,
            }
            return list(fused)
        threshold = self._numeric_score(threshold_raw)
        if threshold is None or not 0.0 <= threshold <= 1.0:
            raise ValueError("min_applicability_for_prediction must be null or within [0, 1]")
        eligible = [
            target for target in fused
            if target.plato is not None
            and self._numeric_score(target.plato.applicability_score) is not None
            and float(target.plato.applicability_score) >= threshold
        ]
        resolution["mechanism_applicability_filter"] = {
            "status": "applied",
            "threshold": threshold,
            "input_targets": len(fused),
            "eligible_targets": len(eligible),
            "excluded_targets": len(fused) - len(eligible),
            "missing_applicability_excluded": sum(
                target.plato is None or self._numeric_score(target.plato.applicability_score) is None
                for target in fused
            ),
        }
        if not eligible:
            raise ValueError(
                f"Unable to explain {drug_name}: no reliable PLATO targets meet the configured "
                f"applicability threshold {threshold}"
            )
        return eligible

    def artifact_hash_inventory(
        self,
        cell_line: str = "",
        disease_label: str = "",
    ) -> Dict[str, Any]:
        """Return allowlist candidates; callers must verify them independently."""
        used_line, _, _ = self.context.expression_profile(
            str(cell_line).strip() or self.config.get("default_cell_line", "default")
        )
        used_disease = (
            disease_label if disease_label in self.context.disease_sets
            else self.context.default_disease_label if not disease_label else ""
        )
        expression_metadata = self.context.resource_metadata.get("expression", {})
        disease_metadata = self.context.resource_metadata.get("disease", {})
        target_inventory = {
            "PLATO_cache": {
                "root": str(self.targets.cache.root) if self.targets.cache else "",
                "hash_scope": "pair_specific_raw_exports",
                "source_authenticity": (
                    "requires_independent_trusted_export_sha256_match; "
                    "cache self-consistency and user attestation are insufficient"
                ),
                "trusted_export_sha256": sorted(self._configured_trusted_export_hashes(
                    self.config.get("plato", {}).get("trusted_export_sha256", [])
                )),
            }
        }
        return {
            **target_inventory,
            "ppi": {
                "sha256": self.context.resource_metadata.get("ppi", {}).get("artifact_sha256", ""),
                "artifact_kind": self.context.resource_metadata.get("ppi", {}).get("artifact_kind", ""),
            },
            "expression": {
                "sha256": expression_metadata.get("artifact_sha256_by_cell_line", {}).get(used_line, ""),
                "artifact_kind": expression_metadata.get("artifact_kind_by_cell_line", {}).get(used_line, ""),
                "cell_line": used_line,
            },
            "disease": {
                "sha256": disease_metadata.get("artifact_sha256_by_label", {}).get(used_disease, ""),
                "artifact_kind": disease_metadata.get("artifact_kind_by_label", {}).get(used_disease, ""),
                "label": used_disease,
            },
            "pathways": {
                "sha256": sorted({value for value in self.context.pathway_artifact_hashes.values() if value}),
                "artifact_kinds": sorted(set(self.context.pathway_artifact_kinds.values())),
            },
            "warning": "inventory only; source authenticity is not verified until compared with independent trusted records",
        }

    def set_prediction_evaluator(self, evaluator: Optional[Callable[[Mapping[str, Any]], Any]]) -> None:
        """Attach a model-specific baseline/perturbation evaluator.

        The callback receives a JSON-like request. It should return a numeric
        synergy score or ``{"score": value, ...}``. Without this callback the
        result is explicitly labelled as a post-hoc biological hypothesis.
        """
        if evaluator is not None and not callable(evaluator):
            raise TypeError("prediction evaluator must be callable or None")
        self.prediction_evaluator = evaluator

    @staticmethod
    def _source_label(targets: Sequence[FusedTarget]) -> str:
        statuses = {target.evidence_status for target in targets}
        if "known_and_predicted" in statuses or ({"known_only", "predicted_only"} <= statuses):
            return "local+PLATO"
        if statuses == {"predicted_only"}:
            return "PLATO"
        return "local"

    def _rank_modules(
        self,
        ranked1: Sequence[ContextScore],
        ranked2: Sequence[ContextScore],
        disease_label: str = "",
    ) -> ModuleExplanation:
        score1 = {item.target_id.upper(): item.final_score for item in ranked1}
        score2 = {item.target_id.upper(): item.final_score for item in ranked2}
        enrichment1 = {
            str(item["pathway"]): item
            for item in self.context.pathway_enrichment(score1)
        }
        enrichment2 = {
            str(item["pathway"]): item
            for item in self.context.pathway_enrichment(score2)
        }
        hits: List[ModuleHit] = []
        for name, genes in self.context.pathways.items():
            matched1 = sorted(member for member in genes if member in score1)
            matched2 = sorted(member for member in genes if member in score2)
            if not matched1 and not matched2:
                continue
            support1 = sum(score1[member] for member in matched1) / len(matched1) if matched1 else None
            support2 = sum(score2[member] for member in matched2) / len(matched2) if matched2 else None
            available = [value for value in (support1, support2) if value is not None]
            score = sum(available) / len(available)
            stat1 = enrichment1.get(name)
            stat2 = enrichment2.get(name)
            supported1 = bool(stat1 and float(stat1.get("q", 1.0)) <= 0.05)
            supported2 = bool(stat2 and float(stat2.get("q", 1.0)) <= 0.05)
            if not (supported1 or supported2):
                continue
            hits.append(ModuleHit(
                module_id=name,
                pathway_name=name,
                module_support_score=float(score),
                pathway_source=self.context.pathway_sources.get(name, ""),
                pathway_version=self.context.pathway_versions.get(name, ""),
                pathway_taxid=self.context.pathway_taxids.get(name, ""),
                pathway_identifier_namespace=self.context.pathway_identifier_namespaces.get(name, ""),
                pathway_mapping_version=self.context.pathway_mapping_versions.get(name, ""),
                pathway_content_sha256=self.context.pathway_content_hashes.get(name, ""),
                pathway_artifact_kind=self.context.pathway_artifact_kinds.get(name, ""),
                pathway_artifact_sha256=self.context.pathway_artifact_hashes.get(name, ""),
                num_genes=len(set(matched1) | set(matched2)),
                pathway_gene_count=len(genes),
                matched_genes_drug1=matched1,
                matched_genes_drug2=matched2,
                drug1_support_score=support1,
                drug2_support_score=support2,
                drug1_enrichment_p=float(stat1["p"]) if stat1 else None,
                drug1_enrichment_q=float(stat1["q"]) if stat1 else None,
                drug2_enrichment_p=float(stat2["p"]) if stat2 else None,
                drug2_enrichment_q=float(stat2["q"]) if stat2 else None,
                enrichment_supported=True,
                shared_enrichment_supported=supported1 and supported2,
            ))
        hits.sort(key=lambda item: (
            item.shared_enrichment_supported,
            -max(
                item.drug1_enrichment_q if item.drug1_enrichment_q is not None else 1.0,
                item.drug2_enrichment_q if item.drug2_enrichment_q is not None else 1.0,
            ),
            item.module_support_score,
        ), reverse=True)
        disease = self.context.disease_genes(disease_label)
        disease_hits = [item for item in hits if self.context.pathways.get(item.module_id, set()) & disease]
        return ModuleExplanation(hits[:10], disease_hits[:10])

    @staticmethod
    def _numeric_score(value: Any) -> Optional[float]:
        if isinstance(value, Mapping):
            value = value.get("score", value.get("predicted_synergy"))
        try:
            score = float(value)
        except (TypeError, ValueError):
            return None
        return score if math.isfinite(score) else None

    @staticmethod
    def _json_safe(value: Any) -> Any:
        if isinstance(value, Mapping):
            return {str(key): SynergyInterpreter._json_safe(item) for key, item in value.items()}
        if isinstance(value, (list, tuple, set)):
            return [SynergyInterpreter._json_safe(item) for item in value]
        if isinstance(value, float):
            return value if math.isfinite(value) else None
        if value is None or isinstance(value, (str, int, bool)):
            return value
        return str(value)

    def _evaluate_model_fidelity(
        self,
        drug1_name: str,
        drug2_name: str,
        drug1_smiles: str,
        drug2_smiles: str,
        cell_line: str,
        predicted_synergy: float,
        advanced: AdvancedResult,
        supplied: Optional[Mapping[str, Any]],
    ) -> Dict[str, Any]:
        if supplied is not None:
            baseline = self._numeric_score(supplied.get("baseline_prediction", predicted_synergy))
            stronger_is_higher = bool(self.config.get("credibility", {}).get("higher_synergy_is_stronger", True))
            credibility_config = self.config.get("credibility", {})
            min_perturbations = max(3, int(credibility_config.get("min_model_perturbations", 3)))
            configured_relative_change = float(
                credibility_config.get("min_relative_prediction_change", 0.05)
            )
            configured_tolerance = float(
                credibility_config.get("baseline_tolerance_fraction", 0.01)
            )
            if not math.isfinite(configured_relative_change) or configured_relative_change < 0:
                raise ValueError("min_relative_prediction_change must be finite and non-negative")
            if not math.isfinite(configured_tolerance) or configured_tolerance < 0:
                raise ValueError("baseline_tolerance_fraction must be finite and non-negative")
            min_relative_change = max(0.05, configured_relative_change)
            tolerance_fraction = min(0.01, configured_tolerance)
            evaluations: List[Dict[str, Any]] = []
            rejected = 0
            raw_perturbations = supplied.get("perturbations", [])
            if not isinstance(raw_perturbations, list):
                raw_perturbations = []
                rejected += 1
            for item in raw_perturbations:
                if not isinstance(item, Mapping) or baseline is None:
                    rejected += 1
                    continue
                score_after = self._numeric_score(
                    item.get("prediction_after", item.get("score_after", item.get("score")))
                )
                if score_after is None:
                    rejected += 1
                    continue
                delta = score_after - baseline
                perturbation = self._json_safe(
                    item.get("perturbation", item.get("kind", "precomputed"))
                )
                kind = str(item.get("kind") or "").strip()
                if not kind and isinstance(perturbation, Mapping):
                    kind = str(perturbation.get("kind") or "").strip()
                kind = kind or "precomputed"
                relative_change = abs(delta) / max(1.0, abs(baseline))
                supports = delta < 0 if stronger_is_higher else delta > 0
                evaluations.append({
                    "perturbation": perturbation,
                    "kind": kind,
                    "prediction_after": score_after,
                    "prediction_delta": delta,
                    "relative_prediction_change": relative_change,
                    "supports_explanation": supports,
                    "material_support": bool(supports and relative_change >= min_relative_change),
                })
            status = "evaluated" if baseline is not None and evaluations else "partial"
            provided_prediction = self._numeric_score(
                supplied.get("provided_prediction", predicted_synergy)
            )
            provided_prediction = predicted_synergy if provided_prediction is None else provided_prediction
            baseline_difference = (
                baseline - provided_prediction if baseline is not None else None
            )
            baseline_tolerance = tolerance_fraction * max(1.0, abs(provided_prediction))
            baseline_consistent = bool(
                baseline_difference is not None
                and abs(baseline_difference) <= baseline_tolerance
            )
            support_rate = (
                sum(item["supports_explanation"] for item in evaluations) / len(evaluations)
                if evaluations else None
            )
            material_support_rate = (
                sum(item["material_support"] for item in evaluations) / len(evaluations)
                if evaluations else None
            )
            perturbation_type_count = len({item["kind"] for item in evaluations})
            meets_thresholds = bool(
                status == "evaluated"
                and baseline_consistent
                and len(evaluations) >= min_perturbations
                and perturbation_type_count >= 2
                and float(material_support_rate or 0.0) >= 2.0 / 3.0
            )
            return {
                "status": status,
                "method": str(supplied.get("method") or "").strip() or "precomputed_model_attribution",
                "interpretation_type": "model_linked_sensitivity" if status == "evaluated" else "post_hoc_biological_hypothesis",
                "baseline_prediction": baseline if baseline is not None else predicted_synergy,
                "provided_prediction": provided_prediction,
                "baseline_difference": baseline_difference,
                "baseline_tolerance": baseline_tolerance,
                "baseline_consistent": baseline_consistent,
                "higher_synergy_is_stronger": stronger_is_higher,
                "support_rate": support_rate,
                "material_support_rate": material_support_rate,
                "minimum_perturbations": min_perturbations,
                "minimum_relative_prediction_change": min_relative_change,
                "unique_perturbation_count": len(evaluations),
                "perturbation_type_count": perturbation_type_count,
                "perturbations": evaluations,
                "rejected_perturbations": rejected,
                "evidence_source": "caller_supplied_unverified",
                "meets_internal_sensitivity_thresholds": meets_thresholds,
                "independently_validated": False,
                "independent_validation_status": "not_performed",
                "warning": "caller_supplied_fidelity_evidence_not_independently_verified",
            }
        if self.prediction_evaluator is None:
            return {
                "status": "not_evaluated",
                "method": "none",
                "interpretation_type": "post_hoc_biological_hypothesis",
                "baseline_prediction": predicted_synergy,
                "support_rate": None,
                "perturbations": [],
                "evidence_source": "none",
                "independently_validated": False,
                "warning": "predicted_synergy_stored_but_no_model_ablation_evaluator_attached",
            }

        common = {
            "drug1_name": drug1_name,
            "drug2_name": drug2_name,
            "drug1_smiles": drug1_smiles,
            "drug2_smiles": drug2_smiles,
            "cell_line": cell_line,
        }
        try:
            baseline_raw = self.prediction_evaluator({**common, "kind": "baseline", "predicted_synergy": predicted_synergy})
            baseline = self._numeric_score(baseline_raw)
        except Exception as exc:
            return {
                "status": "error",
                "method": "prediction_evaluator",
                "interpretation_type": "post_hoc_biological_hypothesis",
                "baseline_prediction": predicted_synergy,
                "perturbations": [],
                "evidence_source": "runtime_evaluator",
                "independently_validated": False,
                "error": f"{type(exc).__name__}: {exc}",
            }
        if baseline is None:
            return {
                "status": "error",
                "method": "prediction_evaluator",
                "interpretation_type": "post_hoc_biological_hypothesis",
                "baseline_prediction": predicted_synergy,
                "perturbations": [],
                "evidence_source": "runtime_evaluator",
                "independently_validated": False,
                "error": "evaluator_did_not_return_a_finite_baseline_score",
            }

        limit = max(0, int(self.config.get("credibility", {}).get("max_model_perturbations", 10)))
        stronger_is_higher = bool(self.config.get("credibility", {}).get("higher_synergy_is_stronger", True))
        evaluations: List[Dict[str, Any]] = []
        errors: List[str] = []
        seen_perturbations = set()
        for perturbation in advanced.counterfactuals[:limit]:
            identity = (
                perturbation.get("kind"),
                perturbation.get("drug"),
                perturbation.get("target_id"),
                perturbation.get("module"),
                tuple(perturbation.get("path", [])),
            )
            if identity in seen_perturbations:
                continue
            seen_perturbations.add(identity)
            try:
                raw = self.prediction_evaluator({**common, "kind": "perturbation", "perturbation": perturbation})
                score_after = self._numeric_score(raw)
            except Exception as exc:
                errors.append(f"{perturbation.get('kind', 'unknown')}:{type(exc).__name__}:{exc}")
                continue
            if score_after is None:
                errors.append(f"{perturbation.get('kind', 'unknown')}:non_finite_score")
                continue
            delta = score_after - baseline
            supports = delta < 0 if stronger_is_higher else delta > 0
            evaluations.append({
                "perturbation": perturbation,
                "prediction_after": score_after,
                "prediction_delta": delta,
                "supports_explanation": supports,
            })
        status = "evaluated" if evaluations else "partial"
        credibility_config = self.config.get("credibility", {})
        min_perturbations = max(3, int(credibility_config.get("min_model_perturbations", 3)))
        configured_relative_change = float(credibility_config.get("min_relative_prediction_change", 0.05))
        configured_tolerance = float(credibility_config.get("baseline_tolerance_fraction", 0.01))
        if not math.isfinite(configured_relative_change) or configured_relative_change < 0:
            raise ValueError("min_relative_prediction_change must be finite and non-negative")
        if not math.isfinite(configured_tolerance) or configured_tolerance < 0:
            raise ValueError("baseline_tolerance_fraction must be finite and non-negative")
        min_relative_change = max(0.05, configured_relative_change)
        tolerance_fraction = min(0.01, configured_tolerance)
        baseline_tolerance = tolerance_fraction * max(1.0, abs(predicted_synergy))
        baseline_difference = baseline - predicted_synergy
        baseline_consistent = abs(baseline_difference) <= baseline_tolerance
        for item in evaluations:
            item["relative_prediction_change"] = abs(item["prediction_delta"]) / max(1.0, abs(baseline))
            item["material_support"] = bool(
                item["supports_explanation"] and item["relative_prediction_change"] >= min_relative_change
            )
        material_support_rate = (
            sum(item["material_support"] for item in evaluations) / len(evaluations) if evaluations else None
        )
        perturbation_type_count = len({
            item.get("perturbation", {}).get("kind")
            for item in evaluations
            if isinstance(item.get("perturbation"), Mapping)
        })
        meets_internal_sensitivity_thresholds = bool(
            status == "evaluated"
            and baseline_consistent
            and len(evaluations) >= min_perturbations
            and perturbation_type_count >= 2
            and float(material_support_rate or 0.0) >= 2.0 / 3.0
        )
        return {
            "status": status,
            "method": "user_attached_prediction_evaluator",
            "interpretation_type": "model_linked_sensitivity" if status == "evaluated" else "post_hoc_biological_hypothesis",
            "baseline_prediction": baseline,
            "provided_prediction": predicted_synergy,
            "baseline_difference": baseline_difference,
            "baseline_tolerance": baseline_tolerance,
            "baseline_consistent": baseline_consistent,
            "higher_synergy_is_stronger": stronger_is_higher,
            "support_rate": sum(item["supports_explanation"] for item in evaluations) / len(evaluations) if evaluations else None,
            "material_support_rate": material_support_rate,
            "minimum_perturbations": min_perturbations,
            "minimum_relative_prediction_change": min_relative_change,
            "unique_perturbation_count": len(evaluations),
            "perturbation_type_count": perturbation_type_count,
            "perturbations": evaluations,
            "errors": errors,
            "evidence_source": "runtime_evaluator",
            "meets_internal_sensitivity_thresholds": meets_internal_sensitivity_thresholds,
            "independently_validated": False,
            "independent_validation_status": "not_performed",
        }

    def _top_k_stability(
        self,
        ranked1: Sequence[ContextScore],
        ranked2: Sequence[ContextScore],
        baseline: MechanismResult,
        top_k: int,
        disease_label: str,
    ) -> Dict[str, Any]:
        config = self.config.get("credibility", {})
        alternatives = sorted({
            int(value) for value in config.get("top_k_stability", [5, 10, 20])
            if int(value) > 0 and int(value) != top_k
        })
        # Put the actual primary run first so an equivalent effective target
        # count cannot cause it to be discarded as a duplicate.
        requested = [top_k, *alternatives]
        network = self.config.get("network", {})
        stability_randomizations = int(config.get("stability_randomizations", 1000))
        if stability_randomizations < 0:
            raise ValueError("stability_randomizations must be non-negative")
        configured_score_tolerance = float(config.get("max_top_k_mechanism_score_delta", 0.15))
        if not math.isfinite(configured_score_tolerance) or configured_score_tolerance < 0:
            raise ValueError("max_top_k_mechanism_score_delta must be a finite non-negative number")
        score_tolerance = min(0.15, configured_score_tolerance)
        results: List[Dict[str, Any]] = []
        seen_effective = set()
        for value in requested:
            effective = (min(value, len(ranked1)), min(value, len(ranked2)))
            if effective in seen_effective:
                continue
            seen_effective.add(effective)
            is_baseline = value == top_k
            mechanism = baseline if is_baseline else classify_mechanism(
                ranked1,
                ranked2,
                self.context,
                top_k=value,
                disease_label=disease_label,
                network_randomizations=stability_randomizations,
                random_seed=int(network.get("random_seed", 42)),
                min_network_coverage=float(network.get("min_mapping_coverage", 0.5)),
            )
            results.append({
                "top_k": value,
                "effective_targets": effective,
                "mechanism_type": mechanism.mechanism_type,
                "mechanism_score": mechanism.confidence,
                "network_randomizations": mechanism.metrics.get("cheng", {}).get("n_randomizations", 0),
                "is_primary_baseline": is_baseline,
            })
        labels = {item["mechanism_type"] for item in results}
        scores = [float(item["mechanism_score"]) for item in results]
        score_delta = max(scores) - min(scores) if scores else 0.0
        evaluable = len(results) >= 2
        label_stable = len(labels) <= 1
        score_stable = score_delta <= score_tolerance
        baseline_agreement = all(
            item["mechanism_type"] == baseline.mechanism_type
            and abs(float(item["mechanism_score"]) - float(baseline.confidence)) <= score_tolerance
            for item in results
        )
        randomization_budget_sufficient = stability_randomizations >= 1000
        descriptive_stable = bool(
            evaluable and label_stable and score_stable and baseline_agreement
        )
        stable = bool(descriptive_stable and randomization_budget_sufficient)
        return {
            "evaluable": evaluable,
            "stable": stable,
            "descriptive_stable": descriptive_stable,
            "label_stable": label_stable,
            "score_stable": score_stable,
            "agrees_with_primary_baseline": baseline_agreement,
            "primary_baseline": {
                "top_k": top_k,
                "mechanism_type": baseline.mechanism_type,
                "mechanism_score": baseline.confidence,
                "network_randomizations": baseline.metrics.get("cheng", {}).get("n_randomizations", 0),
            },
            "mechanism_score_delta": score_delta,
            "maximum_allowed_score_delta": score_tolerance,
            "tested": results,
            "agreement_fraction": max((sum(item["mechanism_type"] == label for item in results) for label in labels), default=0) / max(1, len(results)),
            "network_randomizations_per_test": stability_randomizations,
            "randomization_budget_sufficient": randomization_budget_sufficient,
            "meets_full_robustness_criteria": stable,
        }

    def _build_credibility(
        self,
        profile1: DrugTargetExplanation,
        profile2: DrugTargetExplanation,
        mechanism: MechanismResult,
        context_status: Mapping[str, Any],
        model_fidelity: Mapping[str, Any],
        stability: Mapping[str, Any],
    ) -> Dict[str, Any]:
        profiles = (profile1, profile2)
        plato_only_mode = True
        selected = [profile1.top_targets, profile2.top_targets]
        all_selected = selected[0] + selected[1]
        # Target-source quality is assessed on each drug's complete eligible
        # PLATO collection, before top-k mechanism selection.  Top-k/context
        # ranking must not manufacture better target-source quality.
        eligible_by_drug = [profile1.all_targets, profile2.all_targets]
        all_eligible = eligible_by_drug[0] + eligible_by_drug[1]
        local_supported_fractions = [
            sum(item.evidence_status != "predicted_only" for item in items) / max(1, len(items))
            for items in selected
        ]
        artifact_local_fractions = [
            sum(item.local_evidence_tier == "artifact_backed_quantitative" for item in items) / max(1, len(items))
            for items in selected
        ]
        both_have_local_majority = min(local_supported_fractions, default=0.0) >= 0.5
        external_only = [
            profile.drug_name for profile, items in zip(profiles, selected)
            if items and all(item.evidence_status == "predicted_only" for item in items)
        ]
        plato_items_by_drug = [
            [item for item in items if item.source_evidence.get("plato")]
            for items in eligible_by_drug
        ]
        plato_items = plato_items_by_drug[0] + plato_items_by_drug[1]
        aligned_items = [item for item in all_selected if item.evidence_status == "known_and_predicted"]
        aligned = len(aligned_items)
        plato_config = self.config.get("plato", {})
        applicability_threshold_raw = plato_config.get("min_applicability_for_corroboration")
        applicability_threshold = (
            None if applicability_threshold_raw is None
            else self._numeric_score(applicability_threshold_raw)
        )
        if applicability_threshold_raw is not None and applicability_threshold is None:
            raise ValueError("min_applicability_for_corroboration must be finite or null")

        cache_hits = [
            profile.resolution.get("cache", {}) for profile in profiles
            if str(profile.resolution.get("cache", {}).get("status", "")).startswith("cache_hit")
        ]
        schema_fixtures_verified = bool(cache_hits) and all(
            cache.get("schema_fixture_verified") is True for cache in cache_hits
        )
        strict_plato_cache_ready = bool(
            plato_only_mode
            and len(cache_hits) == len(profiles)
            and all(
                not self._plato_cache_integrity_failures(
                    profile.resolution,
                    plato_config.get("provider_release", "r37"),
                    plato_config.get("trusted_export_sha256", []),
                )
                for profile in profiles
            )
        )

        def qualified_corroboration(item: ContextScore) -> bool:
            evidence = item.source_evidence.get("plato") or {}
            applicability = self._numeric_score(evidence.get("applicability_score"))
            return bool(
                schema_fixtures_verified
                and applicability_threshold is not None
                and evidence.get("reliable") is True
                and applicability is not None
                and applicability >= applicability_threshold
            )

        qualified_aligned = sum(qualified_corroboration(item) for item in aligned_items)
        applicability_known = sum(
            item.source_evidence.get("plato", {}).get("applicability_score") is not None for item in plato_items
        )
        ppi_coverage = min(
            mechanism.coverage.get("drug1", {}).get("ppi_mapping_fraction", 0.0),
            mechanism.coverage.get("drug2", {}).get("ppi_mapping_fraction", 0.0),
        )
        expression_coverage = min(
            mechanism.coverage.get("drug1", {}).get("expression_fraction", 0.0),
            mechanism.coverage.get("drug2", {}).get("expression_fraction", 0.0),
        )
        cheng = mechanism.metrics.get("cheng", {})
        parser_warnings = list(dict.fromkeys(
            warning for cache in cache_hits for warning in cache.get("parser_warnings", [])
        ))
        releases_consistent = not any(profile.resolution.get("release_mismatch") for profile in profiles)
        fusion_stats = [profile.resolution.get("fusion", {}) for profile in profiles]
        local_resolution = [profile.resolution.get("local", {}) for profile in profiles]
        local_species_resolved = not any(
            stats.get("local_origin_unspecified", 0) for stats in fusion_stats
        )
        local_target_types_resolved = not any(
            stats.get("local_target_type_unspecified", 0) or stats.get("local_target_type_unknown", 0)
            for stats in fusion_stats
        )
        def explicitly_human(item: ContextScore) -> bool:
            evidence = item.source_evidence.get("plato")
            if not isinstance(evidence, Mapping):
                return False
            organism = " ".join(str(evidence.get("organism", "")).strip().lower().split())
            taxid = str(evidence.get("taxid", "")).strip()
            warnings_for_target = set(evidence.get("parse_warnings", []) or [])
            return bool(
                taxid == "9606"
                or (
                    organism in {"human", "homo sapiens", "h sapiens"}
                    and "organism_inferred_from_target_name" not in warnings_for_target
                )
            )

        all_selected_human = bool(all_eligible) and all(explicitly_human(item) for item in all_eligible)
        resources = context_status.get("resources", {})
        ppi_metadata = resources.get("ppi", {})
        ppi_audit = ppi_metadata.get("confidence_audit", {})
        credibility_config = self.config.get("credibility", {})
        min_expression = max(
            1000, int(credibility_config.get("min_expression_background_genes", 1000))
        )
        full_expression = max(
            5000, int(credibility_config.get("full_expression_background_genes", 5000))
        )
        if min_expression < 1 or full_expression < min_expression:
            raise ValueError("expression background thresholds must be positive and full >= minimum")
        expression_gene_count = int(context_status.get("expression_gene_count", 0))
        expression_declared_complete = context_status.get("expression_background_complete")
        expression_informative = context_status.get("expression_background_informative") is True
        expression_background_ready = bool(
            expression_gene_count >= min_expression
            and expression_declared_complete is True
            and expression_informative
        )
        expression_background_full = bool(
            expression_gene_count >= full_expression
            and expression_declared_complete is True
            and expression_informative
        )
        action_coverage = min(
            mechanism.coverage.get("drug1", {}).get("qualified_action_direction_fraction", 0.0),
            mechanism.coverage.get("drug2", {}).get("qualified_action_direction_fraction", 0.0),
        )
        action_conflicts = sum(
            int(mechanism.coverage.get(key, {}).get("qualified_action_conflict_count", 0))
            for key in ("drug1", "drug2")
        )
        ppi_nodes = int(ppi_metadata.get("nodes", 0))
        ppi_edges = int(ppi_metadata.get("edges", 0))
        ppi_component_fraction = float(ppi_metadata.get("largest_connected_component_fraction", 0.0) or 0.0)
        ppi_source = str(ppi_metadata.get("source", "")).strip()
        ppi_version = str(ppi_metadata.get("version", "")).strip()
        ppi_taxid = str(ppi_metadata.get("taxid", "")).strip()
        ppi_namespace = str(ppi_metadata.get("identifier_namespace", "")).strip().lower()
        ppi_mapping_version = str(ppi_metadata.get("mapping_version", "")).strip()
        ppi_content_hash = str(ppi_metadata.get("content_sha256", "")).strip()
        ppi_quality_ready = bool(
            ppi_nodes >= 1000
            and ppi_edges >= 3000
            and ppi_component_fraction >= 0.5
            and ppi_source not in {"", "in_memory", "built_in_demo"}
            and ppi_version
            and ppi_taxid == "9606"
            and ppi_namespace in {"gene_symbol", "hgnc_symbol"}
            and ppi_mapping_version
            and ppi_content_hash
            and float(ppi_audit.get("invalid_confidence_fraction", 1.0)) == 0.0
            and float(ppi_audit.get("unweighted_edge_fraction", 1.0)) <= 0.2
        )

        strengths: List[str] = []
        limitations: List[str] = ["mechanism_score_is_uncalibrated"]
        warnings: List[str] = []
        if qualified_aligned:
            strengths.append("qualified_local_and_PLATO_target_corroboration")
        elif aligned:
            limitations.append("PLATO_alignment_not_counted_as_qualified_corroboration")
            if applicability_threshold is None:
                limitations.append("PLATO_applicability_threshold_not_configured")
        if stability.get("stable"):
            strengths.append("mechanism_stable_across_tested_top_k")
        elif not stability.get("evaluable"):
            limitations.append("top_k_stability_not_evaluable_with_available_target_counts")
        elif stability.get("descriptive_stable") and not stability.get("randomization_budget_sufficient"):
            limitations.append("top_k_stability_randomization_budget_insufficient")
            warnings.append("top_k_stability_not_claimed_below_randomization_budget")
        else:
            limitations.append("mechanism_sensitive_to_top_k")
            warnings.append("mechanism_sensitive_to_top_k")
        if cheng.get("quality_status") == "adequate" and cheng.get("statistically_supported"):
            strengths.append("degree_matched_network_support")
        else:
            limitations.append("network_support_not_statistically_adequate")
        if model_fidelity.get("status") == "evaluated":
            strengths.append("prediction_model_perturbations_evaluated")
            if not model_fidelity.get("meets_internal_sensitivity_thresholds"):
                limitations.append("model_fidelity_does_not_meet_internal_sensitivity_thresholds")
        else:
            limitations.append("prediction_model_fidelity_not_evaluated")
        if context_status.get("is_demo"):
            limitations.append("demo_context_not_for_production")
        if context_status.get("expression_fallback"):
            limitations.append("requested_cell_line_used_fallback_expression")
        if not expression_background_ready:
            limitations.append("expression_background_below_minimum_context_threshold")
        if expression_declared_complete is False:
            limitations.append("expression_background_declared_partial")
        elif expression_declared_complete is None:
            limitations.append("expression_background_completeness_not_declared")
        if not expression_informative:
            limitations.append("expression_background_low_dynamic_range")
        local_artifacts_audited = all(
            item.get("status") == "ok"
            and item.get("schema") == "rich"
            and bool(item.get("database_sha256"))
            and int(item.get("schema_user_version") or 0) > 0
            and bool(str(item.get("dataset_version") or "").strip())
            and bool(str(item.get("mapping_version") or "").strip())
            and item.get("database_hash_status") in {"sqlite_static_snapshot", "sqlite_snapshot_with_wal"}
            for item in local_resolution
        )
        if plato_only_mode:
            limitations.append("single_source_PLATO_prediction_without_curated_target_evidence")
            limitations.append("PLATO_action_direction_not_release_validated")
            if strict_plato_cache_ready:
                strengths.append("strict_PLATO_cache_integrity_checks_passed")
        else:
            if min(artifact_local_fractions, default=0.0) < 0.5:
                limitations.append("insufficient_artifact_backed_quantitative_local_target_evidence")
            if not local_artifacts_audited:
                limitations.append("local_target_evidence_not_resolved_from_hashed_rich_database")
            if action_coverage < 0.5:
                limitations.append("insufficient_local_action_direction_evidence")
        if action_conflicts:
            limitations.append("conflicting_local_and_external_action_directions")
            warnings.append("target_action_direction_conflict")
        if external_only and not plato_only_mode:
            limitations.append("one_or_more_drugs_depend_only_on_external_predictions")
        elif not plato_only_mode and min(local_supported_fractions, default=0.0) < 0.5:
            limitations.append("mechanism_majority_depends_on_external_predictions")
        if plato_items and applicability_known < len(plato_items):
            limitations.append("PLATO_applicability_domain_incomplete")
        if cache_hits and not schema_fixtures_verified:
            limitations.append("PLATO_schema_not_verified_by_controlled_fixture")
        if cache_hits and not strict_plato_cache_ready:
            limitations.append("PLATO_export_not_independently_authenticated_and_schema_verified")
        if not releases_consistent:
            limitations.append("PLATO_release_mismatch")
        if not plato_only_mode and any(stats.get("local_origin_unspecified", 0) for stats in fusion_stats):
            limitations.append("local_target_species_assumed_for_legacy_records")
        if not plato_only_mode and any(stats.get("local_target_type_unspecified", 0) for stats in fusion_stats):
            limitations.append("local_target_type_unspecified")
        if not plato_only_mode and any(stats.get("local_target_type_unknown", 0) for stats in fusion_stats):
            limitations.append("local_target_type_unrecognized")
        if any(stats.get("ambiguous_alignments", 0) for stats in fusion_stats):
            limitations.append("ambiguous_target_alignments_retained_unmerged")
            warnings.append("ambiguous_target_alignment")
        if any(stats.get("local_identifier_conflicts", 0) for stats in fusion_stats):
            limitations.append("conflicting_local_target_identifiers_retained_unmerged")
            warnings.append("local_target_identifier_conflict")
        if any(stats.get("target_symbol_disagreements", 0) for stats in fusion_stats):
            limitations.append("target_symbol_disagreement_on_strong_identifier_match")
            warnings.append("target_symbol_disagreement")
        if ppi_audit.get("invalid_confidence_edges", 0):
            limitations.append("PPI_contains_invalid_confidence_values")
            warnings.append("invalid_PPI_confidence_values")
        if ppi_audit.get("assumed_unweighted_edges", 0):
            limitations.append("PPI_contains_topology_only_edges")
        if not ppi_quality_ready:
            limitations.append("PPI_resource_below_production_credibility_threshold")
        if not plato_only_mode and any(item.get("status") == "database_error" for item in local_resolution):
            limitations.append("local_target_database_error")
            warnings.append("local_target_database_error")
        if not plato_only_mode and any(item.get("status") == "database_missing" for item in local_resolution):
            limitations.append("local_target_database_missing")
        if not plato_only_mode and any(item.get("status") == "drug_identity_conflict" for item in local_resolution):
            limitations.append("local_drug_identity_conflict")
            warnings.append("local_drug_identity_conflict")
        limitations.extend(parser_warnings)

        network_coverage = cheng.get("network_coverage", {})
        disease_mapping_fraction = network_coverage.get("disease_module", {}).get("mapping_fraction", 0.0)
        disease_mapped_count = int(network_coverage.get("disease_module", {}).get("mapped", 0))
        network_pattern_to_mechanism = {
            "P1": "high_overlap",
            "P2": "complementary_exposure",
            "P4": "single_exposure",
            "P4_edge": "overlapping_single_exposure",
        }
        network_context_hit = bool(
            cheng.get("quality_status") == "adequate"
            and cheng.get("pair_interaction_supported")
            and cheng.get("pattern") in {"P1", "P2"}
            and network_pattern_to_mechanism.get(cheng.get("pattern")) == mechanism.mechanism_type
        )
        directional_functional_types = {
            "vertical_pathway_coblockade",
        }
        target_overlap_supported = bool(
            mechanism.mechanism_type == "target_overlap_or_redundancy" and mechanism.shared_targets
        )
        functional_supported = mechanism.mechanism_type in directional_functional_types
        pathway_evidence_traceable = bool(
            not mechanism.shared_pathways
            or all(
                self.context.pathway_sources.get(name, "") not in {"", "in_memory", "built_in_demo"}
                and bool(self.context.pathway_versions.get(name, ""))
                and self.context.pathway_taxids.get(name, "") == "9606"
                and self.context.pathway_identifier_namespaces.get(name, "").lower() in {"gene_symbol", "hgnc_symbol"}
                and bool(self.context.pathway_mapping_versions.get(name, ""))
                and bool(self.context.pathway_content_hashes.get(name, ""))
                and self.context.pathway_artifact_kinds.get(name, "") == "file"
                and bool(self.context.pathway_artifact_hashes.get(name, ""))
                for name in mechanism.shared_pathways
            )
        )
        mechanism_supported = bool(network_context_hit or target_overlap_supported or functional_supported)
        local_genes1 = {
            item.target_id.upper() for item in selected[0] if item.evidence_status != "predicted_only"
        }
        local_genes2 = {
            item.target_id.upper() for item in selected[1] if item.evidence_status != "predicted_only"
        }
        if target_overlap_supported:
            shared_set = set(mechanism.shared_targets)
            action_by_drug = mechanism.metrics.get("target_action", {})
            driver_directions_supported = all(
                action_by_drug.get(key, {}).get("qualified_by_target", {}).get(target) == "inhibitory"
                for key in ("drug1", "drug2")
                for target in shared_set
            )
            mechanism_local_driver_supported = bool(
                shared_set <= local_genes1
                and shared_set <= local_genes2
                and all(
                    item.local_evidence_tier == "artifact_backed_quantitative"
                    for items in selected
                    for item in items
                    if item.target_id.upper() in shared_set
                )
                and driver_directions_supported
            )
        elif network_context_hit:
            action_by_drug = mechanism.metrics.get("target_action", {})
            driver_directions_supported = all(
                action_by_drug.get(key, {}).get("qualified_by_target", {}).get(item.target_id.upper()) == "inhibitory"
                for key, items in zip(("drug1", "drug2"), selected)
                for item in items
            )
            mechanism_local_driver_supported = driver_directions_supported and all(
                item.evidence_status != "predicted_only"
                and item.local_evidence_tier == "artifact_backed_quantitative"
                for item in all_selected
            )
        else:
            shared_pathway_genes = set().union(*(
                self.context.pathways.get(name, set()) for name in mechanism.shared_pathways
            )) if mechanism.shared_pathways else set()
            action_by_drug = mechanism.metrics.get("target_action", {})
            qualified_pathway_drivers = []
            for key, items in zip(("drug1", "drug2"), selected):
                statuses = action_by_drug.get(key, {}).get("qualified_by_target", {})
                qualified_pathway_drivers.append(any(
                    item.target_id.upper() in shared_pathway_genes
                    and item.local_evidence_tier == "artifact_backed_quantitative"
                    and statuses.get(item.target_id.upper()) == "inhibitory"
                    for item in items
                ))
            mechanism_local_driver_supported = bool(
                functional_supported
                and pathway_evidence_traceable
                and all(qualified_pathway_drivers)
            )
            driver_directions_supported = bool(functional_supported and all(qualified_pathway_drivers))
        pathway_disease_hit = any(
            self.context.pathways.get(name, set()) & self.context.disease_genes(
                str(context_status.get("used_disease_label", ""))
            )
            for name in mechanism.shared_pathways
        )
        direct_disease_hit = bool(
            mechanism.metrics.get("direct_disease_target_fraction_drug1", 0.0) > 0
            and mechanism.metrics.get("direct_disease_target_fraction_drug2", 0.0) > 0
        )
        cancer_relevance_hit = bool(network_context_hit or pathway_disease_hit or direct_disease_hit)
        expression_metadata = resources.get("expression", {})
        used_line = str(context_status.get("used_expression_cell_line", ""))
        expression_source = str(
            expression_metadata.get("source_by_cell_line", {}).get(
                used_line, expression_metadata.get("source", "")
            )
        ).strip()
        expression_version = str(
            expression_metadata.get("version_by_cell_line", {}).get(
                used_line, expression_metadata.get("version", "")
            )
        ).strip()
        expression_taxid = str(
            expression_metadata.get("taxid_by_cell_line", {}).get(
                used_line, expression_metadata.get("taxid", "")
            )
        ).strip()
        expression_namespace = str(
            expression_metadata.get("identifier_namespace_by_cell_line", {}).get(
                used_line, expression_metadata.get("identifier_namespace", "")
            )
        ).strip().lower()
        expression_mapping_version = str(
            expression_metadata.get("mapping_version_by_cell_line", {}).get(
                used_line, expression_metadata.get("mapping_version", "")
            )
        ).strip()
        expression_content_hash = str(
            expression_metadata.get("content_sha256_by_cell_line", {}).get(
                used_line, expression_metadata.get("content_sha256", "")
            )
        ).strip()
        expression_artifact_kind = str(
            expression_metadata.get("artifact_kind_by_cell_line", {}).get(
                used_line, expression_metadata.get("artifact_kind", "")
            )
        ).strip()
        expression_artifact_hash = str(
            expression_metadata.get("artifact_sha256_by_cell_line", {}).get(
                used_line, expression_metadata.get("artifact_sha256", "")
            )
        ).strip().lower()
        disease_metadata = resources.get("disease", {})
        used_disease_label = str(context_status.get("used_disease_label", ""))
        disease_source = str(
            disease_metadata.get("source_by_label", {}).get(
                used_disease_label, disease_metadata.get("source", "")
            )
        ).strip()
        disease_version = str(
            disease_metadata.get("version_by_label", {}).get(
                used_disease_label, disease_metadata.get("version", "")
            )
        ).strip()
        disease_taxid = str(
            disease_metadata.get("taxid_by_label", {}).get(
                used_disease_label, disease_metadata.get("taxid", "")
            )
        ).strip()
        disease_namespace = str(
            disease_metadata.get("identifier_namespace_by_label", {}).get(
                used_disease_label, disease_metadata.get("identifier_namespace", "")
            )
        ).strip().lower()
        disease_mapping_version = str(
            disease_metadata.get("mapping_version_by_label", {}).get(
                used_disease_label, disease_metadata.get("mapping_version", "")
            )
        ).strip()
        disease_content_hash = str(
            disease_metadata.get("content_sha256_by_label", {}).get(
                used_disease_label, disease_metadata.get("content_sha256", "")
            )
        ).strip()
        disease_artifact_kind = str(
            disease_metadata.get("artifact_kind_by_label", {}).get(
                used_disease_label, disease_metadata.get("artifact_kind", "")
            )
        ).strip()
        disease_artifact_hash = str(
            disease_metadata.get("artifact_sha256_by_label", {}).get(
                used_disease_label, disease_metadata.get("artifact_sha256", "")
            )
        ).strip().lower()
        ppi_artifact_kind = str(ppi_metadata.get("artifact_kind", "")).strip()
        ppi_artifact_hash = str(ppi_metadata.get("artifact_sha256", "")).strip().lower()
        contextual_resources_traceable = bool(
            expression_source not in {"", "in_memory", "built_in_demo"}
            and expression_version
            and expression_taxid == "9606"
            and expression_namespace in {"gene_symbol", "hgnc_symbol"}
            and expression_mapping_version
            and expression_content_hash
            and expression_artifact_kind == "file"
            and expression_artifact_hash
            and disease_source not in {"", "user_supplied", "built_in_demo"}
            and disease_version
            and disease_taxid == "9606"
            and disease_namespace in {"gene_symbol", "hgnc_symbol"}
            and disease_mapping_version
            and disease_content_hash
            and disease_artifact_kind == "file"
            and disease_artifact_hash
            and ppi_artifact_kind in {"file", "safe_edge_table_file"}
            and ppi_artifact_hash
        )
        trusted_hashes = credibility_config.get("trusted_artifact_sha256", {})
        if not isinstance(trusted_hashes, Mapping):
            raise TypeError("credibility.trusted_artifact_sha256 must be a mapping")

        def hash_allowlisted(kind: str, actual_hashes: Sequence[str]) -> bool:
            configured = trusted_hashes.get(kind, [])
            values = [configured] if isinstance(configured, str) else list(configured or [])
            allowed = {str(value).strip().lower() for value in values if str(value).strip()}
            actual = [str(value).strip().lower() for value in actual_hashes if str(value).strip()]
            return bool(actual) and all(value in allowed for value in actual)

        local_artifact_hashes = sorted({
            str(item.get("database_sha256") or "").lower() for item in local_resolution
            if item.get("database_sha256")
        })
        pathway_artifact_hashes = sorted({
            self.context.pathway_artifact_hashes.get(name, "").lower()
            for name in mechanism.shared_pathways
            if self.context.pathway_artifact_hashes.get(name, "")
        })
        contextual_allowlist_status = {
            "ppi": hash_allowlisted("ppi", [ppi_artifact_hash]),
            "expression": hash_allowlisted("expression", [expression_artifact_hash]),
            "disease": hash_allowlisted("disease", [disease_artifact_hash]),
            "pathways": (
                hash_allowlisted("pathways", pathway_artifact_hashes)
                if functional_supported else None
            ),
        }
        trusted_plato_hashes = self._configured_trusted_export_hashes(
            plato_config.get("trusted_export_sha256", [])
        )
        plato_source_hashes = [
            str(cache.get("source", {}).get("sha256", "")).strip().lower()
            for cache in cache_hits
            if isinstance(cache.get("source"), Mapping)
        ]
        plato_export_hashes_allowlisted = bool(
            len(plato_source_hashes) == len(profiles)
            and all(value and value in trusted_plato_hashes for value in plato_source_hashes)
        )
        artifact_allowlist_status = {
            "PLATO_raw_export": plato_export_hashes_allowlisted,
            **contextual_allowlist_status,
        }
        applicable_allowlist_results = [
            value for value in artifact_allowlist_status.values() if value is not None
        ]
        artifact_hashes_allowlisted = bool(applicable_allowlist_results) and all(
            value is True for value in applicable_allowlist_results
        )
        if artifact_hashes_allowlisted:
            strengths.append("required_artifact_hashes_match_explicit_allowlist")
        runtime_integrity = context_status.get("runtime_integrity", {})
        mutated_resources = set(runtime_integrity.get("mutated_resources", []))
        required_mutations = mutated_resources & {"ppi", "expression", "disease"}
        if functional_supported and "pathways" in mutated_resources:
            required_mutations.add("pathways")
        runtime_resources_clean = not required_mutations
        context_ready = (
            not context_status.get("is_demo")
            and not context_status.get("expression_fallback")
            and context_status.get("disease_gene_count", 0) >= 5
            and disease_mapped_count >= 5
            and disease_mapping_fraction >= 0.5
            and ppi_coverage >= 0.5
            and expression_coverage >= 0.5
            and expression_background_ready
            and ppi_quality_ready
            and contextual_resources_traceable
            and cancer_relevance_hit
            and runtime_resources_clean
        )
        if not context_ready:
            limitations.append("cancer_context_not_effectively_mapped_to_both_drugs")
        if not contextual_resources_traceable:
            limitations.append("context_resource_source_or_version_not_traceable")
        if not artifact_hashes_allowlisted:
            limitations.append("artifact_hash_not_allowlisted")
        if not runtime_resources_clean:
            limitations.append("context_resource_mutated_after_provenance_was_recorded")
            warnings.extend(
                f"runtime_resource_mutation:{resource}" for resource in sorted(required_mutations)
            )
        if disease_mapped_count < 5:
            limitations.append("too_few_mapped_disease_genes")
        if not mechanism_supported:
            limitations.append("primary_mechanism_not_directly_supported")
        if not mechanism_local_driver_supported:
            limitations.append(
                "primary_mechanism_lacks_release_validated_PLATO_driver_direction"
                if plato_only_mode
                else "primary_mechanism_not_supported_by_qualified_local_drivers"
            )
        if mechanism_supported and not driver_directions_supported:
            limitations.append("mechanism_driver_direction_requires_signed_disease_role_evidence")
        if mechanism.shared_pathways and not pathway_evidence_traceable:
            limitations.append("pathway_evidence_source_or_version_not_traceable")
        if not all_selected_human:
            limitations.append("selected_targets_not_all_explicitly_human")
        limitations.append("independent_validation_not_performed")

        integrity_warning_names = {
            "cached_payload_hash_mismatch",
            "cached_payload_query_structure_mismatch",
            "cached_query_binding_status_mismatch",
            "normalized_targets_hash_mismatch",
            "normalized_targets_hash_missing",
            "parser_output_changed_since_import",
            "cache_query_structure_mismatch",
            "audit_metadata_hash_missing",
            "audit_metadata_hash_mismatch",
        }
        integrity_issue = any(
            warning in integrity_warning_names
            or warning.startswith("cached_parser_report_")
            or warning.startswith("cached_provider_")
            or warning == "cache_parser_version_differs_from_runtime"
            for warning in parser_warnings
        )
        identity_conflict = not plato_only_mode and any(
            item.get("status") == "drug_identity_conflict" for item in local_resolution
        )
        invalid_reasons: List[str] = []
        if required_mutations:
            invalid_reasons.append("runtime_context_resource_mutation")
        if integrity_issue:
            invalid_reasons.append("PLATO_cache_or_parser_integrity_failure")
        if identity_conflict:
            invalid_reasons.append("local_drug_identity_conflict")

        def clip01(value: Any) -> float:
            number = self._numeric_score(value)
            if number is None:
                return 0.0
            return max(0.0, min(1.0, float(number)))

        explicit_human_fraction = (
            sum(explicitly_human(item) for item in all_eligible) / len(all_eligible)
            if all_eligible else 0.0
        )
        qualified_corroboration_fraction = qualified_aligned / max(1, len(all_selected))
        min_local_fraction = min(local_supported_fractions, default=0.0)
        min_artifact_local_fraction = min(artifact_local_fractions, default=0.0)
        plato_count = max(1, len(plato_items))
        target_evidence_per_drug: Dict[str, Dict[str, Any]] = {}
        for index, (profile, items) in enumerate(zip(profiles, plato_items_by_drug), start=1):
            count = len(items)
            denominator = max(1, count)
            identifier_fraction = sum(
                bool(
                    str(item.source_evidence.get("plato", {}).get("target_symbol", "")).strip()
                    or str(item.source_evidence.get("plato", {}).get("uniprot_id", "")).strip()
                )
                for item in items
            ) / denominator
            reliability_reported_fraction = sum(
                item.source_evidence.get("plato", {}).get("reliable") is not None
                for item in items
            ) / denominator
            reliable_fraction = sum(
                item.source_evidence.get("plato", {}).get("reliable") is True
                for item in items
            ) / denominator
            applicability_reported_fraction = sum(
                item.source_evidence.get("plato", {}).get("applicability_score") is not None
                for item in items
            ) / denominator
            score_reported_fraction = sum(
                self._numeric_score(item.source_evidence.get("plato", {}).get("score")) is not None
                for item in items
            ) / denominator
            human_fraction = sum(explicitly_human(item) for item in items) / denominator
            cache_ready_for_drug = not self._plato_cache_integrity_failures(
                profile.resolution,
                plato_config.get("provider_release", "r37"),
                plato_config.get("trusted_export_sha256", []),
            )
            value = clip01(
                0.25 * human_fraction
                + 0.20 * identifier_fraction
                + 0.15 * reliability_reported_fraction
                + 0.15 * reliable_fraction
                + 0.10 * applicability_reported_fraction
                + 0.05 * score_reported_fraction
                + 0.10 * float(cache_ready_for_drug)
            ) if count else 0.0
            target_evidence_per_drug[f"drug{index}"] = {
                "drug_name": profile.drug_name,
                "value": round(value, 6),
                "parsed_target_count": int(
                    profile.resolution.get("cache", {}).get("parsed_targets", 0) or 0
                ),
                "eligible_target_count": count,
                "mechanism_applicability_filter": dict(
                    profile.resolution.get("mechanism_applicability_filter", {})
                ),
                "explicit_human_fraction": human_fraction,
                "PLATO_identifier_fraction": identifier_fraction,
                "PLATO_reliability_reported_fraction": reliability_reported_fraction,
                "PLATO_reliable_fraction": reliable_fraction,
                "PLATO_applicability_reported_fraction": applicability_reported_fraction,
                "PLATO_score_reported_fraction": score_reported_fraction,
                "strict_PLATO_cache_ready": cache_ready_for_drug,
            }
        target_evidence_value = min(
            (float(item["value"]) for item in target_evidence_per_drug.values()),
            default=0.0,
        )
        # Pooled values are presentation-only.  They never determine the
        # target component, which uses the weaker drug via minimum aggregation.
        plato_identifier_fraction = sum(
            bool(
                str(item.source_evidence.get("plato", {}).get("target_symbol", "")).strip()
                or str(item.source_evidence.get("plato", {}).get("uniprot_id", "")).strip()
            ) for item in plato_items
        ) / plato_count
        plato_reliability_reported_fraction = sum(
            item.source_evidence.get("plato", {}).get("reliable") is not None for item in plato_items
        ) / plato_count
        plato_reliable_fraction = sum(
            item.source_evidence.get("plato", {}).get("reliable") is True for item in plato_items
        ) / plato_count
        plato_applicability_fraction = applicability_known / plato_count
        plato_score_reported_fraction = sum(
            self._numeric_score(item.source_evidence.get("plato", {}).get("score")) is not None
            for item in plato_items
        ) / plato_count

        conflict_free = 1.0 if action_conflicts == 0 else 0.0
        direction_support_value = 0.0 if plato_only_mode else clip01(
            0.55 * action_coverage
            + 0.35 * float(bool(mechanism_local_driver_supported))
            + 0.10 * conflict_free
        )

        disease_effective_mapping = clip01(disease_mapping_fraction) * clip01(
            disease_mapped_count / 5.0
        )
        expression_background_quality = (
            clip01(expression_gene_count / max(1, full_expression))
            if expression_declared_complete is True and expression_informative else 0.0
        )
        cancer_context_value = clip01(
            0.25 * clip01(ppi_coverage)
            + 0.25 * clip01(expression_coverage)
            + 0.20 * disease_effective_mapping
            + 0.15 * expression_background_quality
            + 0.15 * float(bool(cancer_relevance_hit))
        )

        convergent_mechanism_evidence = bool(
            network_context_hit
            or (target_overlap_supported and driver_directions_supported)
            or (functional_supported and pathway_evidence_traceable and pathway_disease_hit)
        )
        mechanism_support_value = clip01(
            0.45 * clip01(mechanism.confidence)
            + 0.35 * float(bool(mechanism_supported))
            + 0.20 * float(convergent_mechanism_evidence)
        )

        artifact_allowlist_fraction = sum(
            value is True for value in applicable_allowlist_results
        ) / max(1, len(applicable_allowlist_results))
        plato_integrity_or_not_applicable = float(strict_plato_cache_ready)
        if plato_only_mode:
            provenance_integrity_value = clip01(
                0.35 * float(strict_plato_cache_ready)
                + 0.25 * float(bool(contextual_resources_traceable))
                + 0.20 * float(bool(ppi_quality_ready))
                + 0.20 * artifact_allowlist_fraction
            )
        else:
            provenance_integrity_value = clip01(
                0.30 * artifact_allowlist_fraction
                + 0.25 * float(bool(contextual_resources_traceable))
                + 0.20 * float(bool(local_artifacts_audited))
                + 0.15 * float(bool(ppi_quality_ready))
                + 0.10 * plato_integrity_or_not_applicable
            )

        if stability.get("evaluable"):
            agreement = clip01(stability.get("agreement_fraction"))
            tolerance = max(0.0, float(stability.get("maximum_allowed_score_delta", 0.0) or 0.0))
            delta = max(0.0, float(stability.get("mechanism_score_delta", 0.0) or 0.0))
            score_stability = 1.0 if tolerance == 0.0 and delta == 0.0 else (
                max(0.0, 1.0 - delta / tolerance) if tolerance > 0 else 0.0
            )
            randomization_factor = clip01(
                float(stability.get("network_randomizations_per_test", 0) or 0) / 1000.0
            )
            robustness_value = clip01(agreement * score_stability * randomization_factor)
        else:
            robustness_value = 0.0

        model_status = str(model_fidelity.get("status", ""))
        model_source = str(model_fidelity.get("evidence_source", ""))
        model_fidelity_value = 0.0
        model_contradiction = False
        if model_status == "evaluated":
            minimum_perturbations = max(
                1, int(model_fidelity.get("minimum_perturbations", credibility_config.get("min_model_perturbations", 3)) or 1)
            )
            perturbation_count = int(
                model_fidelity.get("unique_perturbation_count", len(model_fidelity.get("perturbations", []))) or 0
            )
            perturbation_types = int(model_fidelity.get("perturbation_type_count", 0) or 0)
            support_rate = clip01(
                model_fidelity.get("material_support_rate", model_fidelity.get("support_rate", 0.0))
            )
            baseline_consistent = float(model_fidelity.get("baseline_consistent", True) is True)
            model_fidelity_value = clip01(
                baseline_consistent
                * clip01(perturbation_count / minimum_perturbations)
                * clip01(perturbation_types / 2.0)
                * support_rate
            )
            if model_source == "caller_supplied_unverified":
                # Caller-supplied attribution is useful as a falsification
                # check, but it is not an independent trust anchor.  It may
                # trigger a contradiction penalty below; it must never add
                # positive credibility points.
                model_fidelity_value = 0.0
            model_contradiction = bool(
                model_fidelity.get("support_rate") is not None
                and float(model_fidelity.get("support_rate") or 0.0) < 0.5
            )
        independent_validation_value = 0.0

        target_component_inputs = {
            "aggregation": "minimum_across_drugs",
            "selection_stage": "complete_eligible_PLATO_targets_before_top_k_context_ranking",
            "per_drug": target_evidence_per_drug,
        }
        direction_component_inputs = {
            "status": "not_scored_without_release_validated_PLATO_action_semantics",
            "PLATO_action_is_provenance_only": True,
            "action_conflict_count": action_conflicts,
        }
        provenance_component_inputs = {
            "artifact_allowlist_fraction": artifact_allowlist_fraction,
            "context_resources_traceable": contextual_resources_traceable,
            "PPI_quality_ready": ppi_quality_ready,
            "PLATO_cache_integrity": strict_plato_cache_ready,
            "PLATO_raw_export_sha256_allowlisted": plato_export_hashes_allowlisted,
            "platform_export_authenticity": (
                "independent_trusted_export_sha256_match"
                if plato_export_hashes_allowlisted else "not_authenticated"
            ),
        }

        component_specs = {
            "target_evidence": (0.22, target_evidence_value, target_component_inputs),
            "direction_support": (0.14, direction_support_value, direction_component_inputs),
            "cancer_context": (0.18, cancer_context_value, {
                "PPI_mapping_fraction": ppi_coverage,
                "expression_mapping_fraction": expression_coverage,
                "disease_effective_mapping_fraction": disease_effective_mapping,
                "expression_background_quality": expression_background_quality,
                "cancer_relevance_hit": cancer_relevance_hit,
            }),
            "mechanism_support": (0.14, mechanism_support_value, {
                "mechanism_rule_strength": mechanism.confidence,
                "primary_mechanism_supported": mechanism_supported,
                "convergent_mechanism_evidence": convergent_mechanism_evidence,
            }),
            "provenance_integrity": (0.12, provenance_integrity_value, provenance_component_inputs),
            "robustness": (0.06, robustness_value, {
                "top_k_stability_evaluable": bool(stability.get("evaluable")),
                "agreement_fraction": stability.get("agreement_fraction"),
                "mechanism_score_delta": stability.get("mechanism_score_delta"),
                "randomizations": stability.get("network_randomizations_per_test"),
            }),
            "model_fidelity": (0.04, model_fidelity_value, {
                "status": model_status,
                "evidence_source": model_source,
                "contradicts_explanation": model_contradiction,
            }),
            "independent_validation": (0.10, independent_validation_value, {
                "status": "not_performed",
                "required_evidence": "independent mechanism-specific validation with negative controls",
            }),
        }
        components: Dict[str, Dict[str, Any]] = {}
        component_status_overrides = {
            "direction_support": "not_applicable",
            "robustness": None if stability.get("evaluable") else "unmeasured",
            "model_fidelity": (
                "unverified_caller_evidence_noncontributing"
                if model_source == "caller_supplied_unverified"
                else None if model_status == "evaluated" else "unmeasured"
            ),
            "independent_validation": "unmeasured",
        }
        for name, (weight, value, inputs) in component_specs.items():
            max_points = 100.0 * weight
            component_value = clip01(value)
            status_override = component_status_overrides.get(name)
            components[name] = {
                "value": round(component_value, 6),
                "weight": weight,
                "max_points": max_points,
                "points": round(max_points * component_value, 6),
                "status": status_override or (
                    "contributing" if component_value > 0 else "not_demonstrated"
                ),
                "inputs": inputs,
            }
        raw_score = round(sum(item["points"] for item in components.values()), 6)
        weight_gap_points = round(max(0.0, 100.0 - raw_score), 6)
        weight_gap_breakdown: Dict[str, float] = {}
        for component in components.values():
            status = str(component["status"])
            bucket = "partially_demonstrated" if status == "contributing" else status
            gap = max(0.0, float(component["max_points"]) - float(component["points"]))
            weight_gap_breakdown[bucket] = weight_gap_breakdown.get(bucket, 0.0) + gap
        weight_gap_breakdown = {
            key: round(value, 6) for key, value in weight_gap_breakdown.items()
        }

        target_factor_candidates: List[Tuple[str, float]] = []
        if not plato_only_mode and not local_species_resolved:
            target_factor_candidates.append(("local_target_species_unresolved", 0.60))
        if not plato_only_mode and len(external_only) >= 2:
            target_factor_candidates.append(("both_drugs_external_predictions_only", 0.55))
        elif not plato_only_mode and len(external_only) == 1:
            target_factor_candidates.append(("one_drug_external_predictions_only", 0.75))
        if not plato_only_mode and not local_target_types_resolved:
            target_factor_candidates.append(("local_target_type_unresolved", 0.80))
        if any(
            stats.get("ambiguous_alignments", 0)
            or stats.get("local_identifier_conflicts", 0)
            or stats.get("target_symbol_disagreements", 0)
            for stats in fusion_stats
        ):
            target_factor_candidates.append(("target_alignment_or_identifier_conflict", 0.80))
        target_quality_factor = min(
            (value for _, value in target_factor_candidates), default=1.0
        )
        context_factor_candidates: List[Tuple[str, float]] = []
        if context_status.get("is_demo"):
            context_factor_candidates.append(("demo_context", 0.50))
        if context_status.get("expression_fallback"):
            context_factor_candidates.append(("expression_cell_line_fallback", 0.80))
        context_quality_factor = min(
            (value for _, value in context_factor_candidates), default=1.0
        )
        contradiction_factor = 0.70 if action_conflicts else 1.0
        quality_factor_product = (
            target_quality_factor * context_quality_factor * contradiction_factor
        )

        score_status = "valid"
        final_score: Optional[float] = round(raw_score * quality_factor_product, 1)
        if not selected[0] or not selected[1]:
            score_status = "insufficient"
            final_score = 0.0
        if invalid_reasons:
            score_status = "invalid"
            final_score = None
            warnings.extend(f"credibility_invalid:{reason}" for reason in invalid_reasons)
        if model_contradiction and final_score is not None:
            pre_adjustment = final_score
            final_score = round(max(0.0, final_score - 4.0), 1)
            model_adjustment = {
                "code": "model_fidelity_contradicts_explanation",
                "delta": round(final_score - pre_adjustment, 1),
            }
        else:
            model_adjustment = None

        limitations = list(dict.fromkeys(limitations + invalid_reasons))
        informational_limitations = {
            "mechanism_score_is_uncalibrated",
            "prediction_model_fidelity_not_evaluated",
            "independent_validation_not_performed",
            "PLATO_applicability_threshold_not_configured",
            "top_k_stability_not_evaluable_with_available_target_counts",
            "top_k_stability_randomization_budget_insufficient",
            "model_fidelity_does_not_meet_internal_sensitivity_thresholds",
        }
        if not network_context_hit:
            informational_limitations.add("network_support_not_statistically_adequate")
        if ppi_quality_ready and float(ppi_audit.get("unweighted_edge_fraction", 0.0)) <= 0.2:
            informational_limitations.add("PPI_contains_topology_only_edges")
        if not functional_supported:
            informational_limitations.add("pathway_evidence_source_or_version_not_traceable")
        if both_have_local_majority:
            informational_limitations.update({
                "PLATO_alignment_not_counted_as_qualified_corroboration",
                "PLATO_applicability_domain_incomplete",
                "PLATO_schema_not_verified_by_controlled_fixture",
                "PLATO_real_export_status_is_user_declared_only",
                "plato_query_binding_user_asserted_unverified",
            })
        high_severity_limitations = [
            item for item in limitations if item not in informational_limitations
        ]
        warnings.extend(
            f"credibility_limitation:{item}" for item in high_severity_limitations
        )
        warnings = list(dict.fromkeys(warnings))
        adjustments = [model_adjustment] if model_adjustment else []
        coverage_report: Dict[str, Any] = {
            "target_mode": "plato_only" if plato_only_mode else "hybrid",
            "mechanism": mechanism.coverage,
            "action_direction_fraction": action_coverage,
            "action_conflict_count": action_conflicts,
            "PLATO_applicability_fraction": (
                applicability_known / max(1, len(plato_items)) if plato_items else None
            ),
            "PLATO_applicability_reported_fraction": (
                applicability_known / max(1, len(plato_items)) if plato_items else None
            ),
            "top_k_stability": stability,
            "context": dict(context_status),
            "artifact_allowlist": artifact_allowlist_status,
            "source_authenticity": "user_attested_not_cryptographically_authenticated",
        }
        if plato_only_mode:
            coverage_report.update({
                "strict_PLATO_cache_ready": strict_plato_cache_ready,
                "PLATO_identifier_fraction": plato_identifier_fraction,
                "PLATO_reliability_reported_fraction": plato_reliability_reported_fraction,
                "PLATO_reliable_fraction": plato_reliable_fraction,
                "PLATO_score_reported_fraction": plato_score_reported_fraction,
                "PLATO_schema_fixture_verified": schema_fixtures_verified,
            })
        else:
            coverage_report.update({
                "source_alignment_fraction": aligned / max(1, len(all_selected)),
                "source_agreement_fraction": qualified_aligned / max(1, len(all_selected)),
                "qualified_source_corroboration_fraction": qualified_corroboration_fraction,
                "PLATO_corroboration_applicability_threshold": applicability_threshold,
                "PLATO_schema_fixture_verified": schema_fixtures_verified,
                "local_supported_fraction": {
                    profile1.drug_name: local_supported_fractions[0],
                    profile2.drug_name: local_supported_fractions[1],
                },
                "artifact_backed_quantitative_local_fraction": {
                    profile1.drug_name: artifact_local_fractions[0],
                    profile2.drug_name: artifact_local_fractions[1],
                },
            })
        report = {
            "evidence_credibility_score": final_score,
            "status": score_status,
            "integrity_valid": score_status != "invalid",
            "mechanism_supported": bool(mechanism_supported),
            "externally_validated": bool(
                model_fidelity.get("independently_validated") is True
            ),
            "validation_semantics": {
                "integrity_valid": "cache, parser, and context integrity gates passed",
                "mechanism_supported": "configured mechanism evidence rules were met",
                "externally_validated": "independent experimental or external validation was supplied",
            },
            "formula_version": "evidence_credibility_v4_verified_plato_only",
            "score_profile": {
                "name": "single_pair",
                "formula_version": "evidence_credibility_v4_verified_plato_only",
                "target_mode": "plato_only",
                "batch_FDR_applied": False,
            },
            "scale": {"minimum": 0.0, "maximum": 100.0, "higher_is_stronger": True},
            "score_semantics": (
                "0-100 continuous evidence support index; not a probability of "
                "biological truth, treatment response, or clinical success"
            ),
            "not_probability": True,
            "not_accuracy_guarantee": True,
            "calibrated_accuracy": None,
            "calibrated_probability": None,
            "raw_score": raw_score,
            "quality_factors": {
                "target_quality": {
                    "value": target_quality_factor,
                    "reasons": [code for code, _ in target_factor_candidates],
                },
                "context_quality": {
                    "value": context_quality_factor,
                    "reasons": [code for code, _ in context_factor_candidates],
                },
                "contradiction": {
                    "value": contradiction_factor,
                    "reasons": ["action_direction_conflict"] if action_conflicts else [],
                },
                "combined": round(quality_factor_product, 6),
            },
            "final_score": final_score,
            "base_score": final_score,
            "batch_adjusted_score": None,
            "missing_weight_points": weight_gap_points,
            "missing_weight_points_breakdown": weight_gap_breakdown,
            "components": components,
            "adjustments": adjustments,
            "invalid_reasons": invalid_reasons,
            "independent_validation_status": "not_performed",
            "coverage": coverage_report,
            "strengths": list(dict.fromkeys(strengths)),
            "limitations": limitations,
            "limitation_severity": {
                item: ("high" if item in high_severity_limitations else "informational")
                for item in limitations
            },
            "warnings": warnings,
        }
        mechanism.credibility_score = final_score
        mechanism.calibrated_probability = None
        mechanism.uncertainty = list(dict.fromkeys(mechanism.uncertainty + limitations))
        return report

    def explain_combination(
        self,
        drug1_name: str,
        drug2_name: str,
        drug1_smiles: str = "",
        drug2_smiles: str = "",
        cell_line: str = "",
        drug1_targets: Optional[Sequence[Any]] = None,
        drug2_targets: Optional[Sequence[Any]] = None,
        predicted_synergy: float = 0.0,
        advanced: Optional[bool] = None,
        disease_label: str = "",
        model_fidelity_evidence: Optional[Mapping[str, Any]] = None,
    ) -> CombinationExplanation:
        self._assert_strict_plato_only_runtime()
        synergy_value = self._numeric_score(predicted_synergy)
        if synergy_value is None:
            raise ValueError("predicted_synergy must be finite")
        if drug1_targets is not None or drug2_targets is not None:
            raise ValueError(
                "Caller-supplied targets are disabled when target_mode='plato_only'"
            )
        fused1, resolution1 = self.targets.resolve(drug1_name, drug1_smiles, drug1_targets)
        fused2, resolution2 = self.targets.resolve(drug2_name, drug2_smiles, drug2_targets)
        self._validate_plato_only_resolution(
            drug1_name, drug1_smiles, fused1, resolution1
        )
        self._validate_plato_only_resolution(
            drug2_name, drug2_smiles, fused2, resolution2
        )
        fused1 = self._filter_plato_applicability(drug1_name, fused1, resolution1)
        fused2 = self._filter_plato_applicability(drug2_name, fused2, resolution2)
        if not fused1:
            raise ValueError(f"Unable to explain {drug1_name}: no eligible cached PLATO human protein targets")
        if not fused2:
            raise ValueError(f"Unable to explain {drug2_name}: no eligible cached PLATO human protein targets")
        effective_cell = str(cell_line).strip() or self.config.get("default_cell_line", "default")
        if self._production_context_loaded:
            used_line, expression_profile, fallback = self.context.expression_profile(effective_cell)
            if not expression_profile or fallback or not used_line:
                available = sorted(self.context.expression)
                preview = ", ".join(available[:12])
                suffix = f" (+{len(available) - 12} more)" if len(available) > 12 else ""
                raise ValueError(
                    f"Production explanation requires an explicit loaded cell line; "
                    f"{effective_cell!r} is unavailable. Available: {preview}{suffix}"
                )
        if not self._batch_context_audit_active:
            self.context.audit_runtime_resources()
        ranked1 = self.context.rank_targets(
            fused1,
            effective_cell,
            disease_label=disease_label,
            weights=self.config.get("context_weights"),
            plato_applicability_threshold=self.config.get("plato", {}).get("min_applicability_for_prediction"),
        )
        ranked2 = self.context.rank_targets(
            fused2,
            effective_cell,
            disease_label=disease_label,
            weights=self.config.get("context_weights"),
            plato_applicability_threshold=self.config.get("plato", {}).get("min_applicability_for_prediction"),
        )
        top_k = int(self.config.get("top_k_targets", 10))
        if top_k <= 0:
            raise ValueError("top_k_targets must be positive")
        profile1 = DrugTargetExplanation(drug1_name, drug1_smiles, ranked1, ranked1[:top_k], resolution1)
        profile2 = DrugTargetExplanation(drug2_name, drug2_smiles, ranked2, ranked2[:top_k], resolution2)
        network = self.config.get("network", {})
        mechanism = classify_mechanism(
            profile1.top_targets,
            profile2.top_targets,
            self.context,
            top_k=top_k,
            disease_label=disease_label,
            network_randomizations=int(network.get("randomizations", 1000)),
            random_seed=int(network.get("random_seed", 42)),
            min_network_coverage=float(network.get("min_mapping_coverage", 0.5)),
        )
        modules = self._rank_modules(profile1.top_targets, profile2.top_targets, disease_label)
        enable_advanced = self.config.get("advanced", True) if advanced is None else advanced
        advanced_result = explain_advanced(
            drug1_name,
            drug2_name,
            profile1.top_targets,
            profile2.top_targets,
            self.context,
            mechanism,
            top_k=top_k,
            disease_label=disease_label,
            sensitivity_randomizations=int(network.get("sensitivity_randomizations", 1000)),
            random_seed=int(network.get("random_seed", 42)),
            min_network_coverage=float(network.get("min_mapping_coverage", 0.5)),
        ) if enable_advanced else AdvancedResult()
        model_fidelity = self._evaluate_model_fidelity(
            drug1_name,
            drug2_name,
            drug1_smiles,
            drug2_smiles,
            effective_cell,
            synergy_value,
            advanced_result,
            model_fidelity_evidence,
        )
        min_expression_background = max(
            1000,
            int(self.config.get("credibility", {}).get("min_expression_background_genes", 1000)),
        )
        context_status = self.context.context_status(
            effective_cell,
            disease_label,
            min_expression_background_genes=min_expression_background,
            audit_resources=False,
        )
        stability = self._top_k_stability(ranked1, ranked2, mechanism, top_k, disease_label)
        credibility = self._build_credibility(
            profile1, profile2, mechanism, context_status, model_fidelity, stability
        )
        source1, source2 = self._source_label(fused1), self._source_label(fused2)
        release1, release2 = resolution1.get("provider_release", "unknown"), resolution2.get("provider_release", "unknown")
        effective_release = release1 if release1 == release2 else "mixed"
        used_disease_label = str(context_status.get("used_disease_label", ""))
        disease_metadata = context_status.get("resources", {}).get("disease", {})
        used_disease_source = disease_metadata.get("source_by_label", {}).get(
            used_disease_label, "unknown"
        )
        return CombinationExplanation(
            drug1_name=drug1_name,
            drug2_name=drug2_name,
            cell_line=effective_cell,
            predicted_synergy=synergy_value,
            drug1_targets=profile1,
            drug2_targets=profile2,
            modules=modules,
            mechanism=mechanism,
            advanced=advanced_result,
            metadata={
                "target_mode": self.config.get("target_mode", "plato_only"),
                "target_source": source1 if source1 == source2 else "mixed",
                "target_source_drug1": source1,
                "target_source_drug2": source2,
                "PLATO_release": effective_release,
                "PLATO_releases": {"drug1": release1, "drug2": release2},
                "PLATO_expected_release": self.config["plato"].get("provider_release", "r37"),
                "advanced_enabled": bool(enable_advanced),
                "disease_context_size": len(self.context.disease_genes(used_disease_label)) if used_disease_label else 0,
                "disease_context_labels": [used_disease_label] if used_disease_label else [],
                "disease_context_source": used_disease_source,
                "context_status": context_status,
            },
            credibility=credibility,
            model_fidelity=model_fidelity,
        )

    def explain_batch(self, combinations: Sequence[Mapping[str, Any]], save_path: Optional[str] = None) -> List[CombinationExplanation]:
        self._assert_strict_plato_only_runtime()
        forbidden = {"database", "local_targets", "drug1_targets", "drug2_targets"}
        for index, item in enumerate(combinations):
            present = sorted(forbidden & set(item))
            if present:
                raise ValueError(
                    f"PLATO-only batch row {index} contains forbidden local target fields: "
                    + ", ".join(present)
                )
        self.context.audit_runtime_resources()
        self._batch_context_audit_active = True
        try:
            results = [self.explain_combination(**dict(item)) for item in combinations]
        finally:
            self._batch_context_audit_active = False
        hypotheses: List[Tuple[float, CombinationExplanation, str]] = []
        for result in results:
            cheng = result.mechanism.metrics.get("cheng", {})
            for key in ("A", "B"):
                value = self._numeric_score(cheng.get(f"p_{key}"))
                if value is not None:
                    hypotheses.append((value, result, key))
        ordered = sorted(enumerate(hypotheses), key=lambda item: item[1][0])
        q_values: Dict[int, float] = {}
        running = 1.0
        total = len(ordered)
        for rank_index in range(total - 1, -1, -1):
            original_index, (p_value, _, _) = ordered[rank_index]
            rank = rank_index + 1
            running = min(running, p_value * total / rank)
            q_values[original_index] = min(1.0, running)
        for original_index, (_, result, key) in enumerate(hypotheses):
            cheng = result.mechanism.metrics.get("cheng", {})
            cheng[f"batch_q_{key}"] = q_values[original_index]
            cheng["batch_multiple_testing_method"] = "Benjamini-Hochberg_across_all_pair_by_drug_tests"
            cheng["batch_hypothesis_count"] = total
        for result in results:
            credibility = result.credibility
            base_score = credibility.get("base_score", credibility.get("final_score"))
            credibility["base_score"] = base_score
            credibility["batch_adjusted_score"] = credibility.get("final_score")
            profile = dict(credibility.get("score_profile", {}))
            profile.update({
                "name": "batch_BH_FDR",
                "formula_version": "evidence_credibility_v4_verified_plato_only",
                "target_mode": "plato_only",
                "batch_FDR_applied": True,
                "batch_FDR_adjustment_applied": False,
                "batch_multiple_testing_method": "Benjamini-Hochberg_across_all_pair_by_drug_tests",
                "batch_hypothesis_count": total,
            })
            credibility["score_profile"] = profile
        network_types = {
            "high_overlap", "complementary_exposure", "single_exposure",
            "overlapping_single_exposure",
        }
        for result in results:
            cheng = result.mechanism.metrics.get("cheng", {})
            if result.mechanism.mechanism_type not in network_types:
                continue
            q_supported = all(
                self._numeric_score(cheng.get(key)) is not None and float(cheng[key]) <= 0.05
                for key in ("batch_q_A", "batch_q_B")
            )
            score = self._numeric_score(
                result.credibility.get("evidence_credibility_score")
            )
            if score is not None and not q_supported:
                limitation = "network_support_not_significant_after_batch_FDR"
                delta = -min(6.0, score)
                adjusted = round(max(0.0, score + delta), 1)
                result.credibility["evidence_credibility_score"] = adjusted
                result.credibility["final_score"] = adjusted
                result.credibility["batch_adjusted_score"] = adjusted
                result.credibility.setdefault("score_profile", {})[
                    "batch_FDR_adjustment_applied"
                ] = True
                result.credibility.setdefault("adjustments", []).append({
                    "code": "batch_FDR_not_significant",
                    "delta": round(delta, 1),
                    "method": "Benjamini-Hochberg_across_all_pair_by_drug_tests",
                })
                result.credibility.setdefault("limitations", []).append(limitation)
                result.credibility.setdefault("warnings", []).append(
                    f"credibility_adjustment:{limitation}"
                )
                result.credibility.setdefault("limitation_severity", {})[limitation] = "high"
                result.mechanism.credibility_score = adjusted
                result.mechanism.uncertainty = list(dict.fromkeys(
                    result.mechanism.uncertainty + [limitation]
                ))
                for key in ("limitations", "warnings"):
                    values = result.credibility.get(key, [])
                    result.credibility[key] = list(dict.fromkeys(values))
        if save_path:
            self.save_results(results, save_path)
        return results

    @staticmethod
    def save_results(results: Sequence[CombinationExplanation], path: str | Path) -> None:
        destination = Path(path)
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(
            json.dumps([result.to_dict() for result in results], ensure_ascii=False, indent=2, allow_nan=False),
            encoding="utf-8",
        )


ExplanationPipeline = SynergyInterpreter


class CounterfactualHypothesisEngine:
    """Compatibility wrapper exposing counterfactuals already computed from a frozen target snapshot."""

    def __init__(self, pipeline: SynergyInterpreter, **_: Any):
        self.pipeline = pipeline

    def generate_hypotheses(self, **kwargs: Any) -> Dict[str, Any]:
        explanation = self.pipeline.explain_combination(**kwargs)
        return {
            "baseline": explanation.to_dict(),
            "perturbations": explanation.advanced.counterfactuals,
            "hypotheses": [item for item in explanation.advanced.counterfactuals if item.get("mechanism_changed")],
            "metadata": {"method": "frozen_fused_target_snapshot"},
        }
