"""Regression test for pair-local production result serialization."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path
import sys
import tempfile

import networkx as nx


sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from InterpreterPLATO import ExplanationPipeline
from test_support import archive_verified_plato_cache, strict_pipeline


def _large_context_status(base: dict) -> tuple[dict, dict[str, str]]:
    hashes = {
        "ppi": "1" * 64,
        "expression": "2" * 64,
        "pathways": "3" * 64,
        "disease": "4" * 64,
    }
    lines = ["A549"] + [f"CELL_{index:02d}" for index in range(1, 70)]
    expression = {
        "source": "DepMap Public",
        "version": "26Q1",
        "cell_lines": lines,
        "gene_counts": {line: 19_214 for line in lines},
        "background_complete_by_cell_line": {line: True for line in lines},
        "source_by_cell_line": {line: "DepMap Public" for line in lines},
        "version_by_cell_line": {line: "26Q1" for line in lines},
        "detection_threshold_by_cell_line": {line: 1.0 for line in lines},
        "taxid_by_cell_line": {line: "9606" for line in lines},
        "identifier_namespace_by_cell_line": {line: "hgnc_symbol" for line in lines},
        "mapping_version_by_cell_line": {line: "HGNC-2026-04-01" for line in lines},
        "artifact_kind_by_cell_line": {line: "file" for line in lines},
        "artifact_sha256_by_cell_line": {line: hashes["expression"] for line in lines},
        "content_sha256_by_cell_line": {
            line: f"{index + 10:064x}" for index, line in enumerate(lines)
        },
    }
    provenance = {}
    for index in range(2_840):
        name = "PATHWAY_HIT" if index == 0 else f"PATHWAY_{index:04d}"
        provenance[name] = {
            "source": "Reactome",
            "version": "96",
            "taxid": "9606",
            "identifier_namespace": "hgnc_symbol",
            "mapping_version": "HGNC-2026-04-01",
            "content_sha256": f"{index + 1:064x}",
            "artifact_kind": "file",
            "artifact_sha256": hashes["pathways"],
        }
    pathways = {
        "source": "Reactome",
        "version": "96",
        "taxid": "9606",
        "identifier_namespace": "hgnc_symbol",
        "mapping_version": "HGNC-2026-04-01",
        "artifact_kind": "file",
        "artifact_sha256": hashes["pathways"],
        "pathway_count": 2_840,
        "gene_union_count": 11_000,
        "sources": ["Reactome"],
        "versions": ["96"],
        "provenance_by_pathway": provenance,
    }
    disease = {
        "source": "COSMIC Cancer Gene Census",
        "version": "v101-GRCh38-Tier1",
        "source_by_label": {"pan_cancer": "COSMIC Cancer Gene Census"},
        "version_by_label": {"pan_cancer": "v101-GRCh38-Tier1"},
        "taxid_by_label": {"pan_cancer": "9606"},
        "identifier_namespace_by_label": {"pan_cancer": "hgnc_symbol"},
        "mapping_version_by_label": {"pan_cancer": "HGNC-2026-04-01"},
        "content_sha256_by_label": {"pan_cancer": "5" * 64},
        "artifact_kind_by_label": {"pan_cancer": "file"},
        "artifact_sha256_by_label": {"pan_cancer": hashes["disease"]},
        "labels": ["pan_cancer"],
        "gene_counts": {"pan_cancer": 584},
        "artifact_kind": "file",
        "artifact_sha256": hashes["disease"],
    }
    context = deepcopy(base)
    context.update({
        "requested_cell_line": "A549",
        "used_expression_cell_line": "A549",
        "used_disease_label": "pan_cancer",
        "disease_label": "pan_cancer",
        "disease_gene_count": 584,
        "resources": {
            "ppi": {
                "source": "STRING",
                "version": "12.0",
                "taxid": "9606",
                "identifier_namespace": "hgnc_symbol",
                "mapping_version": "HGNC-2026-04-01",
                "artifact_kind": "file",
                "artifact_sha256": hashes["ppi"],
                "nodes": 16_036,
                "edges": 234_299,
                "largest_connected_component_nodes": 15_722,
                "largest_connected_component_fraction": 0.98,
                "content_sha256": "6" * 64,
                "confidence_audit": {"invalid_confidence_edges": 0},
            },
            "expression": expression,
            "pathways": pathways,
            "disease": disease,
        },
    })
    return context, hashes


def main() -> int:
    with tempfile.TemporaryDirectory() as directory:
        root = Path(directory)
        cache_dir = root / "plato_cache"
        cache_dir.mkdir()
        sha1 = archive_verified_plato_cache(cache_dir, "CCO", "drug-a", [{
            "target_symbol": "EGFR",
            "target_name": "EGFR",
            "organism": "Homo sapiens",
            "taxid": "9606",
            "target_type": "single protein",
            "score": 12.0,
            "reliable": "yes",
            "max_similarity": 0.80,
        }])
        sha2 = archive_verified_plato_cache(cache_dir, "CCN", "drug-b", [{
            "target_symbol": "RRM1",
            "target_name": "RRM1",
            "organism": "Homo sapiens",
            "taxid": "9606",
            "target_type": "single protein",
            "score": 11.0,
            "reliable": "yes",
            "max_similarity": 0.75,
        }])
        pipeline = strict_pipeline(
            cache_dir,
            [sha1, sha2],
            network={"randomizations": 5, "random_seed": 7},
            credibility={"stability_randomizations": 3},
        )
        graph = nx.Graph()
        graph.add_edges_from([("EGFR", "GRB2"), ("GRB2", "RRM1")])
        pipeline.load_ppi_network(graph)
        pipeline.load_gene_expression({"EGFR": 8.0, "RRM1": 7.0}, cell_line="A549")
        pipeline.load_pathway_data({
            "PATHWAY_HIT": ["EGFR", "RRM1", "GRB2", "AKT1", "MAPK1"],
            "BACKGROUND_PATHWAY": [f"BACKGROUND_{index}" for index in range(300)],
        })
        pipeline.load_disease_genes(["EGFR", "RRM1"], name="pan_cancer")
        result = pipeline.explain_combination(
            "drug-a", "drug-b",
            drug1_smiles="CCO",
            drug2_smiles="CCN",
            cell_line="A549",
            disease_label="pan_cancer",
        )

        full_context, hashes = _large_context_status(
            result.metadata["context_status"]
        )
        result.metadata["context_status"] = full_context
        result.credibility["coverage"]["context"] = full_context

        compact_payload = result.to_dict()
        legacy_payload = deepcopy(compact_payload)
        legacy_payload["metadata"]["context_status"] = deepcopy(full_context)
        legacy_payload["credibility"]["coverage"]["context"] = deepcopy(full_context)
        legacy_size = len(json.dumps(
            legacy_payload, ensure_ascii=False, indent=2, allow_nan=False
        ).encode("utf-8"))
        compact_size = len(json.dumps(
            compact_payload, ensure_ascii=False, indent=2, allow_nan=False
        ).encode("utf-8"))

        serialized_context = compact_payload["metadata"]["context_status"]
        assert serialized_context["artifact_sha256"] == hashes
        assert serialized_context["resources"]["expression"]["used_cell_line"] == "A549"
        assert "cell_lines" not in serialized_context["resources"]["expression"]
        assert serialized_context["resources"]["disease"]["used_label"] == "pan_cancer"
        pathway_output = serialized_context["resources"]["pathways"]
        assert pathway_output["hit_pathways"] == ["PATHWAY_HIT"]
        assert set(pathway_output["provenance_by_hit_pathway"]) == {"PATHWAY_HIT"}
        assert "provenance_by_pathway" not in pathway_output
        assert compact_payload["credibility"]["coverage"]["context"] == serialized_context

        # Serialization must not weaken or mutate the complete in-memory audit.
        assert len(result.metadata["context_status"]["resources"]["expression"]["cell_lines"]) == 70
        assert len(result.metadata["context_status"]["resources"]["pathways"]["provenance_by_pathway"]) == 2_840
        assert len(result.credibility["coverage"]["context"]["resources"]["expression"]["cell_lines"]) == 70
        assert len(result.credibility["coverage"]["context"]["resources"]["pathways"]["provenance_by_pathway"]) == 2_840
        assert compact_size < legacy_size * 0.15

    print(
        "compact result serialization passed: "
        f"legacy={legacy_size} bytes, compact={compact_size} bytes, "
        f"reduction={(1 - compact_size / legacy_size) * 100:.1f}%"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
