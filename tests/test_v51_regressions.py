"""Focused offline regressions for the InterpreterPLATO 5.1 trust boundary."""

from __future__ import annotations

import hashlib
import json
from itertools import permutations
from pathlib import Path
from typing import Any, Dict

import networkx as nx
import pytest

import InterpreterPLATO.interpreter as interpreter_module
from InterpreterPLATO import ExplanationPipeline
from InterpreterPLATO.context import CancerContext, ContextScore
from InterpreterPLATO.mechanism import MechanismResult, classify_mechanism
from InterpreterPLATO.scripts.test_support import (
    archive_verified_plato_cache,
    strict_pipeline,
)
from InterpreterPLATO.targets import (
    ChemblTargetMapper,
    fuse_local_and_plato_targets,
    parse_plato_json,
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _plato_record(symbol: str, **updates: Any) -> Dict[str, Any]:
    record: Dict[str, Any] = {
        "target_symbol": symbol,
        "target_name": symbol,
        "organism": "Homo sapiens",
        "taxid": "9606",
        "target_type": "single protein",
        "score": 12.0,
        "reliable": "yes",
        "max_similarity": 0.8,
    }
    record.update(updates)
    return record


def _context_score(symbol: str) -> ContextScore:
    return ContextScore(
        target_id=symbol,
        target_name=symbol,
        final_score=1.0,
        evidence_score=1.0,
        expression_score=None,
        disease_proximity_score=0.0,
        pathway_score=1.0,
        centrality_score=0.0,
        evidence_status="predicted_only",
        source="PLATO",
        source_evidence={},
    )


def test_hgnc_hyphen_survives_mapping_and_matches_all_contexts(
    tmp_path: Path,
) -> None:
    mapping_path = tmp_path / "target_mapping.tsv"
    mapping_path.write_text(
        "\t".join((
            "target_chembl_id", "target_name", "chembl_target_type",
            "component_uniprot_id", "canonical_hgnc_symbol", "hgnc_id",
            "mapping_status", "target_component_count",
        ))
        + "\n"
        + "\t".join((
            "CHEMBL2632", "HLA class I A", "SINGLE PROTEIN", "P04439",
            "HLA-A", "HGNC:4931", "mapped", "1",
        ))
        + "\n",
        encoding="utf-8",
    )
    mapper = ChemblTargetMapper(
        mapping_path,
        expected_sha256=_sha256(mapping_path),
        source="synthetic ChEMBL",
        release="r37",
    )
    parsed = parse_plato_json({"targets": [{
        **_plato_record(""),
        "target_name": "HLA class I A",
        "target_chembl_id": "CHEMBL2632",
    }]})
    enriched, metadata = mapper.enrich(parsed)

    assert metadata["mapped_targets"] == 1
    assert enriched[0].target_symbol == "HLA-A"
    assert enriched[0].target_id == "HLA-A"
    assert enriched[0].raw["chembl_target_mapping"]["target_symbol_lookup_key"] == "HLAA"

    context = CancerContext()
    graph = nx.Graph()
    graph.add_edge("HLA-A", "B2M", confidence=0.9, distance=1.0 / 0.9)
    context.load_ppi(graph)
    context.load_expression(
        {"HLA-A": 5.0, "B2M": 4.0},
        cell_line="A549",
        background_complete=True,
    )
    context.load_pathways({
        "MHC-I": {"HLA-A", "B2M", "TAP1", "TAP2", "ERAP1"},
    })

    assert enriched[0].target_id in context.ppi
    assert enriched[0].target_id in context.expression_for("A549")
    assert context.pathway_names_for_gene(enriched[0].target_id) == ("MHC-I",)


def test_duplicate_plato_permutations_are_deterministic_and_conflicts_fail_closed() -> None:
    conservative_records = [
        _plato_record(
            "HLA-A", target_chembl_id="CHEMBL2632", uniprot_id="P04439",
            score=12.0, reliable="yes", rank=1, max_similarity=0.8,
        ),
        _plato_record(
            "HLA-A", target_chembl_id="CHEMBL2632", uniprot_id="P04439",
            score=9.0, reliable="no", rank=4, max_similarity=0.4,
        ),
    ]
    merged = [
        parse_plato_json({"targets": list(order)})[0].to_dict()
        for order in permutations(conservative_records)
    ]
    assert merged[0] == merged[1]
    assert merged[0]["target_symbol"] == "HLA-A"
    assert merged[0]["score"] == 9.0
    assert merged[0]["reliable"] is None
    assert merged[0]["rank"] == 4

    conflicting_records = [
        _plato_record(
            "HLA-A", target_chembl_id="CHEMBL-CONFLICT", uniprot_id="P04439",
        ),
        _plato_record(
            "HLA-B", target_chembl_id="CHEMBL-CONFLICT", uniprot_id="P01889",
        ),
    ]
    conflict_results = [
        parse_plato_json({"targets": list(order)})
        for order in permutations(conflicting_records)
    ]
    assert conflict_results[0][0].to_dict() == conflict_results[1][0].to_dict()
    conflict = conflict_results[0][0]
    assert conflict.target_symbol == ""
    assert conflict.target_type == ""
    assert "duplicate_plato_identity_conflict_failed_closed" in conflict.parse_warnings

    fused, stats = fuse_local_and_plato_targets(
        [], conflict_results[0], human_only=True, reliable_only=False,
    )
    assert fused == []
    assert stats["excluded_non_protein"] == 1


def test_caller_supplied_model_fidelity_can_penalize_but_never_add(
    tmp_path: Path,
) -> None:
    cache_dir = tmp_path / "plato_cache"
    cache_dir.mkdir()
    first_hash = archive_verified_plato_cache(
        cache_dir, "CCO", "drug-a", [_plato_record("EGFR")],
    )
    second_hash = archive_verified_plato_cache(
        cache_dir, "CCN", "drug-b", [_plato_record("RRM1")],
    )

    def explain(fidelity: Dict[str, Any] | None = None) -> Dict[str, Any]:
        pipeline = strict_pipeline(cache_dir, [first_hash, second_hash])
        pipeline.load_demo_resources()
        return pipeline.explain_combination(
            "drug-a",
            "drug-b",
            "CCO",
            "CCN",
            predicted_synergy=4.2,
            advanced=False,
            model_fidelity_evidence=fidelity,
        ).to_dict()["credibility"]

    supportive = {
        "baseline_prediction": 4.2,
        "provided_prediction": 4.2,
        "perturbations": [
            {"kind": "remove_target", "prediction_after": 3.0},
            {"kind": "remove_pathway", "prediction_after": 2.5},
            {"kind": "mask_context", "prediction_after": 3.5},
        ],
    }
    contradictory = {
        "baseline_prediction": 4.2,
        "provided_prediction": 4.2,
        "perturbations": [
            {"kind": "remove_target", "prediction_after": 5.0},
            {"kind": "remove_pathway", "prediction_after": 6.0},
            {"kind": "mask_context", "prediction_after": 5.5},
        ],
    }

    baseline = explain()
    supported = explain(supportive)
    contradicted = explain(contradictory)

    assert supported["final_score"] <= baseline["final_score"]
    assert supported["components"]["model_fidelity"]["points"] == 0.0
    assert (
        supported["components"]["model_fidelity"]["status"]
        == "unverified_caller_evidence_noncontributing"
    )
    assert contradicted["final_score"] < baseline["final_score"]
    assert any(
        item["code"] == "model_fidelity_contradicts_explanation"
        and item["delta"] < 0
        for item in contradicted["adjustments"]
    )


def test_failed_production_target_mapper_load_rolls_back_context_and_mapper(
    tmp_path: Path,
) -> None:
    ppi = tmp_path / "ppi.tsv"
    ppi.write_text(
        "gene_a\tgene_b\tcombined_score\nEGFR\tKRAS\t900\n",
        encoding="utf-8",
    )
    expression = tmp_path / "expression.csv"
    expression.write_text(
        "gene,A549\nEGFR,5.1\nKRAS,4.2\nTP53,6.0\n",
        encoding="utf-8",
    )
    pathways = tmp_path / "pathways.gmt"
    pathways.write_text(
        "signaling\ttest\tEGFR\tKRAS\tTP53\tBRAF\tAKT1\n",
        encoding="utf-8",
    )
    disease = tmp_path / "disease.txt"
    disease.write_text("EGFR\nKRAS\n", encoding="utf-8")
    invalid_mapping = tmp_path / "invalid_target_mapping.tsv"
    invalid_mapping.write_text("wrong\tcolumns\n1\t2\n", encoding="utf-8")

    def resource(path: Path, namespace: str = "hgnc_symbol", **extra: Any) -> Dict[str, Any]:
        return {
            "path": path.name,
            "source": "synthetic",
            "version": "test",
            "taxid": "9606",
            "identifier_namespace": namespace,
            "mapping_version": "HGNC-2026-04-01",
            "sha256": _sha256(path),
            **extra,
        }

    configuration = {
        "runtime_manifest": "runtime.lock.json",
        "resources": {
            "ppi": resource(
                ppi, confidence_scale="auto", weight_semantics="confidence",
            ),
            "expression": resource(expression, background_complete=True),
            "pathways": resource(pathways),
            "disease": resource(disease, name="pan_cancer"),
            "target_mapping": resource(
                invalid_mapping, "chembl_target_id_to_hgnc_symbol",
            ),
        },
    }
    config_path = tmp_path / "production.json"
    config_path.write_text(
        json.dumps(configuration, sort_keys=True), encoding="utf-8",
    )
    artifacts = []
    for path in (
        config_path, ppi, expression, pathways, disease, invalid_mapping,
    ):
        artifacts.append({
            "relative_path": path.name,
            "size_bytes": path.stat().st_size,
            "sha256": _sha256(path),
        })
    (tmp_path / "runtime.lock.json").write_text(
        json.dumps({
            "manifest_kind": "runtime_integrity_lock",
            "root": ".",
            "artifacts": artifacts,
        }, sort_keys=True),
        encoding="utf-8",
    )

    pipeline = ExplanationPipeline()
    original_context = pipeline.context
    original_mapper = pipeline.targets.target_mapper
    with pytest.raises(ValueError, match="target_mapping resource did not load"):
        pipeline.load_production_resources(config_path)

    assert pipeline.context is original_context
    assert pipeline.targets.target_mapper is original_mapper
    assert pipeline._production_context_loaded is False


def test_sub_1000_randomizations_cannot_claim_stability(monkeypatch: pytest.MonkeyPatch) -> None:
    pipeline = ExplanationPipeline(config={
        "credibility": {
            "stability_randomizations": 999,
            "top_k_stability": [5, 20],
        },
    })
    baseline = MechanismResult(
        mechanism_type="shared_pathway_co_modulation",
        confidence=0.6,
        description="synthetic stable mechanism",
        metrics={"cheng": {"n_randomizations": 999}},
    )
    monkeypatch.setattr(
        interpreter_module,
        "classify_mechanism",
        lambda *args, **kwargs: baseline,
    )

    result = pipeline._top_k_stability(
        [object()] * 20,
        [object()] * 20,
        baseline,
        top_k=10,
        disease_label="",
    )

    assert result["evaluable"] is True
    assert result["descriptive_stable"] is True
    assert result["randomization_budget_sufficient"] is False
    assert result["stable"] is False
    assert result["meets_full_robustness_criteria"] is False


def test_shared_pathway_requires_bh_fdr_support_to_drive_mechanism() -> None:
    left = [_context_score("A1"), _context_score("A2")]
    right = [_context_score("B1"), _context_score("B2")]
    shared_members = {"A1", "A2", "B1", "B2", "X"}

    raw_only_context = CancerContext()
    raw_only_context.load_pathways({"SHARED": shared_members})
    raw_only = classify_mechanism(
        left, right, raw_only_context, network_randomizations=0,
    )

    assert raw_only.mechanism_type == "unknown"
    assert raw_only.shared_pathways == []
    assert (
        raw_only.metrics["pathway_enrichment"]["descriptive_shared_pathways"]
        == ["SHARED"]
    )
    assert raw_only.metrics["pathway_enrichment"]["significant_shared"] == []
    assert any(
        item["type"] == "descriptive_shared_pathways_not_enriched"
        for item in raw_only.evidence
    )

    enriched_context = CancerContext()
    pathways = {"SHARED": shared_members}
    for index in range(99):
        pathways[f"NULL_{index:02d}"] = {
            f"N{index}_{member}" for member in range(5)
        }
    enriched_context.load_pathways(pathways)
    enriched = classify_mechanism(
        left, right, enriched_context, network_randomizations=0,
    )

    significant = enriched.metrics["pathway_enrichment"]["significant_shared"]
    assert significant and significant[0]["pathway"] == "SHARED"
    assert significant[0]["max_q"] <= 0.05
    assert enriched.shared_pathways == ["SHARED"]
    assert enriched.mechanism_type == "shared_pathway_co_modulation"
