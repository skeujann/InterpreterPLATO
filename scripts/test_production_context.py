"""Offline regression test for the hash-pinned single-pair production path."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import tempfile


sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from InterpreterPLATO import ExplanationPipeline
from test_support import archive_verified_plato_cache, strict_pipeline


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resource(path: Path, source: str, version: str, **extra: object) -> dict[str, object]:
    return {
        "path": path.name,
        "source": source,
        "version": version,
        "taxid": "9606",
        "identifier_namespace": "hgnc_symbol",
        "mapping_version": "HGNC-2026-04-01",
        "sha256": sha256(path),
        **extra,
    }


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="interpreterplato-production-") as temporary:
        root = Path(temporary)

        ppi = root / "ppi.tsv"
        ppi.write_text(
            "gene_a\tgene_b\tcombined_score\nEGFR\tKRAS\t900\n",
            encoding="utf-8",
        )

        expression = root / "expression.csv"
        expression.write_text(
            "gene,A549\nEGFR,5.1\nKRAS,4.2\nTP53,6.0\n",
            encoding="utf-8",
        )
        pathways = root / "pathways.gmt"
        pathways.write_text("EGFR-KRAS signaling\tR-HSA-TEST\tEGFR\tKRAS\n", encoding="utf-8")
        disease = root / "disease.txt"
        disease.write_text("EGFR\nKRAS\n", encoding="utf-8")

        config = {
            "schema_version": 1,
            "runtime_manifest": "runtime.lock.json",
            "resources": {
                "ppi": resource(
                    ppi, "synthetic STRING", "test",
                    confidence_scale="auto", weight_semantics="confidence",
                ),
                "expression": resource(
                    expression, "synthetic DepMap", "test",
                    background_complete=True, detection_threshold=0.0,
                ),
                "pathways": resource(pathways, "synthetic Reactome", "test"),
                "disease": resource(
                    disease, "synthetic COSMIC", "test", name="pan_cancer_tier1",
                ),
            },
        }
        config_path = root / "production.json"
        config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
        runtime_manifest = {
            "schema_version": "1.0",
            "manifest_kind": "runtime_integrity_lock",
            "root": ".",
            "artifacts": [
                {
                    "relative_path": path.name,
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
                for path in (config_path, ppi, expression, pathways, disease)
            ],
        }
        (root / "runtime.lock.json").write_text(
            json.dumps(runtime_manifest, indent=2), encoding="utf-8",
        )

        cache_dir = root / "plato_cache"
        cache_dir.mkdir()
        sha1 = archive_verified_plato_cache(cache_dir, "CCO", "EGFR inhibitor", [{
            "target_symbol": "EGFR",
            "target_name": "EGFR",
            "organism": "Homo sapiens",
            "taxid": "9606",
            "target_type": "single protein",
            "score": 12.0,
            "reliable": "yes",
            "max_similarity": 0.80,
        }])
        sha2 = archive_verified_plato_cache(cache_dir, "CCN", "KRAS inhibitor", [{
            "target_symbol": "KRAS",
            "target_name": "KRAS",
            "organism": "Homo sapiens",
            "taxid": "9606",
            "target_type": "single protein",
            "score": 11.0,
            "reliable": "yes",
            "max_similarity": 0.75,
        }])
        pipe = strict_pipeline(
            cache_dir,
            [sha1, sha2],
            network={"randomizations": 10, "sensitivity_randomizations": 5},
            credibility={"stability_randomizations": 5},
        )
        inventory = pipe.load_production_resources(config_path)
        assert inventory["mode"] == "production"
        assert inventory["cell_lines"] == 1
        assert pipe.context.is_demo_context is False
        assert set(pipe.config["credibility"]["trusted_artifact_sha256"]) >= {
            "ppi", "expression", "pathways", "disease",
        }

        result = pipe.explain_combination(
            drug1_name="EGFR inhibitor",
            drug2_name="KRAS inhibitor",
            drug1_smiles="CCO",
            drug2_smiles="CCN",
            cell_line="A549",
            disease_label="pan_cancer_tier1",
            predicted_synergy=10.0,
            advanced=False,
        )
        payload = result.to_dict()
        assert payload["integrated_explanation"]["summary"]
        assert "demo_context_not_for_production" not in payload["warning_flags"]
        assert payload["cell_line"] == "A549"
        assert payload["credibility"]["components"]["provenance_integrity"]["inputs"][
            "context_resources_traceable"
        ] is True

        expression.write_text(expression.read_text(encoding="utf-8") + "BRAF,3.0\n", encoding="utf-8")
        rejected = False
        try:
            ExplanationPipeline().load_production_resources(config_path)
        except ValueError as exc:
            rejected = "SHA-256 mismatch" in str(exc)
        assert rejected, "tampered production artifact was not rejected"

    print("production context tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
