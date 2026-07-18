"""Strict PLATO-only contract tests for the compact interpreter."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile
from typing import Callable
from collections import Counter


sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from InterpreterPLATO import ExplanationPipeline
from InterpreterPLATO.targets import (
    PlatoPredictionCache,
    inspect_plato_payload,
    fuse_local_and_plato_targets,
    parse_plato_json,
)

from test_support import (
    archive_verified_plato_cache,
    allow_test_schema,
    plato_payload,
    strict_pipeline,
    write_payload,
)


def expect_error(fragment: str, func: Callable[[], object]) -> str:
    try:
        func()
    except Exception as exc:  # noqa: BLE001 - regression helper checks public errors.
        message = str(exc)
        assert fragment.lower() in message.lower(), message
        return message
    raise AssertionError(f"expected error containing {fragment!r}")


def test_target_parsing_filters() -> None:
    payload = plato_payload("CCO", [
        {
            "target_symbol": "EGFR",
            "target_name": "EGFR",
            "organism": "Homo sapiens",
            "taxid": "9606",
            "target_type": "single protein",
            "score": 12.1,
            "reliable": "yes",
        },
        {
            "target_symbol": "MOUSE1",
            "target_name": "Mouse target",
            "organism": "Mus musculus",
            "taxid": "10090",
            "target_type": "single protein",
            "score": 11.0,
            "reliable": "yes",
        },
        {
            "target_symbol": "MIR21",
            "target_name": "miR-21",
            "organism": "Homo sapiens",
            "taxid": "9606",
            "target_type": "nucleic acid",
            "score": 10.0,
            "reliable": "yes",
        },
        {
            "target_symbol": "KRAS",
            "target_name": "KRAS",
            "organism": "Homo sapiens",
            "taxid": "9606",
            "target_type": "single protein",
            "score": 9.0,
            "reliable": "no",
        },
    ])
    parsed = parse_plato_json(payload)
    fused, stats = fuse_local_and_plato_targets(
        [], parsed, human_only=True, reliable_only=True, max_plato_targets=100,
    )
    assert [item.target_id for item in fused] == ["EGFR"]
    assert stats["excluded_non_human"] == 1
    assert stats["excluded_non_protein"] == 1
    assert stats["excluded_unreliable"] == 1

    bad_score = parse_plato_json(plato_payload("CCO", [{
        "target_symbol": "BRAF",
        "target_name": "BRAF",
        "organism": "Homo sapiens",
        "taxid": "9606",
        "target_type": "single protein",
        "score": 99,
        "reliable": "yes",
    }]))[0]
    assert bad_score.score is None
    assert "plato_score_out_of_range_0_13" in bad_score.parse_warnings


def test_real_plato_r37_fixture_contract() -> None:
    fixture = Path(__file__).resolve().parents[1] / "data" / "raw" / "plato" / "r37" / "golden_fixture" / "gemcitabine_target_fishing_response_2026-07-16.json"
    payload = json.loads(fixture.read_text(encoding="utf-8-sig"))
    parsed = parse_plato_json(payload)
    report = inspect_plato_payload(payload, parsed)
    assert len(parsed) == 60
    assert Counter(target.reliable for target in parsed) == Counter({False: 44, True: 16})
    assert parsed[0].target_chembl_id == "CHEMBL1997"
    assert parsed[0].organism == "Homo sapiens"
    assert parsed[0].target_type == "single protein"
    assert parsed[0].score == 13.0
    assert report["schema_fixture_verified"] is True
    assert report["provider_name"] == "PLATO"
    assert report["provider_method"] == "target_fishing"
    assert report["reported_release"] == "r37"
    assert report["reported_release_source"] == "verified_schema_fixture"
    with tempfile.TemporaryDirectory(prefix="interpreterplato-fixture-cache-") as cache_root:
        raw_sha256 = "d2fc1cecf4784261e6674ec21559dd8cbfb60c282ed06c76817884102524484d"
        smiles = "c1cn(c(=O)nc1N)[C@H]2C([C@@H]([C@H](O2)CO)O)(F)F"
        cache = PlatoPredictionCache(cache_root)
        cache.import_file(
            smiles=smiles,
            result_file=fixture,
            drug_name="Gemcitabine",
            provider_release="r37",
            confirmed_real_export=True,
            query_binding_evidence={
                "endpoint": "https://prometheus.farmacia.uniba.it/plato/pred.html",
                "http_method": "POST",
                "submitted_smiles": smiles,
                "response_sha256": raw_sha256,
                "form": {
                    "smiles": smiles,
                    "tipoout": "json",
                    "tipopred": "mussel",
                },
            },
        )
        _, metadata = cache.load(smiles)
        assert metadata["query"]["query_binding_status"] == "plato_single_query_post_verified"
        assert "plato_query_binding_user_asserted_unverified" not in metadata["parser_warnings"]


def test_strict_api_rejects_local_and_mutable_modes(root: Path) -> None:
    expect_error(
        "strict PLATO-only",
        lambda: ExplanationPipeline(config={"target_mode": "hybrid"}),
    )

    pipe = strict_pipeline(root / "empty-cache-a", [])
    expect_error(
        "caller-supplied targets",
        lambda: pipe.explain_combination(
            "drug-a", "drug-b", "CCO", "CCN", drug1_targets=[("EGFR", "EGFR")],
        ),
    )

    mutated = strict_pipeline(root / "empty-cache-b", [])
    mutated.config["target_mode"] = "hybrid"
    expect_error(
        "mutated",
        lambda: mutated.explain_combination("drug-a", "drug-b", "CCO", "CCN"),
    )

    registered = strict_pipeline(root / "empty-cache-c", [])
    registered.targets.register_local_targets("drug-a", [("EGFR", "EGFR")])
    expect_error(
        "registered local targets",
        lambda: registered.explain_combination("drug-a", "drug-b", "CCO", "CCN"),
    )


def test_unverified_and_untrusted_plato_exports_fail(root: Path) -> None:
    unverified = root / "unverified"
    unverified.mkdir()
    source = unverified / "handwritten.json"
    source_hash = write_payload(source, plato_payload("CCO"))
    PlatoPredictionCache(unverified).import_file(
        smiles="CCO",
        result_file=source,
        drug_name="ethanol",
        provider_release="r37",
        confirmed_real_export=True,
    )
    pipe = strict_pipeline(unverified, [source_hash])
    expect_error(
        "schema_fixture_not_verified",
        lambda: pipe.explain_combination("ethanol", "ethylamine", "CCO", "CCN"),
    )

    trusted = root / "trusted"
    trusted.mkdir()
    sha1 = archive_verified_plato_cache(trusted, "CCO", "ethanol")
    sha2 = archive_verified_plato_cache(trusted, "CCN", "ethylamine", [{
        "target_symbol": "RRM1",
        "target_name": "Ribonucleotide reductase catalytic subunit M1",
        "organism": "Homo sapiens",
        "taxid": "9606",
        "target_type": "single protein",
        "score": 11.5,
        "reliable": "yes",
        "max_similarity": 0.77,
    }])
    assert sha1 != sha2
    expect_error(
        "independently_trusted",
        lambda: strict_pipeline(trusted, []).explain_combination(
            "ethanol", "ethylamine", "CCO", "CCN",
        ),
    )


def test_strict_positive_fixture_and_integrity(root: Path) -> None:
    cache_dir = root / "positive"
    cache_dir.mkdir()
    sha1 = archive_verified_plato_cache(cache_dir, "CCO", "ethanol")
    sha2 = archive_verified_plato_cache(cache_dir, "CCN", "ethylamine", [{
        "target_symbol": "RRM1",
        "target_name": "Ribonucleotide reductase catalytic subunit M1",
        "organism": "Homo sapiens",
        "taxid": "9606",
        "target_type": "single protein",
        "score": 11.5,
        "reliable": "yes",
        "max_similarity": 0.77,
    }])

    pipe = strict_pipeline(cache_dir, [sha1, sha2])
    pipe.load_demo_resources()
    result = pipe.explain_combination(
        "ethanol",
        "ethylamine",
        "CCO",
        "CCN",
        predicted_synergy=4.2,
        advanced=False,
    )
    payload = result.to_dict()
    assert payload["result_schema_version"] == "5.1"
    assert payload["metadata"]["target_mode"] == "plato_only"
    assert payload["credibility"]["formula_version"] == "evidence_credibility_v4_verified_plato_only"
    assert payload["credibility"]["score_profile"]["name"] == "single_pair"
    assert payload["credibility"]["base_score"] == payload["credibility"]["final_score"]
    assert payload["credibility"]["batch_adjusted_score"] is None
    target_inputs = payload["credibility"]["components"]["target_evidence"]["inputs"]
    assert target_inputs["aggregation"] == "minimum_across_drugs"
    assert set(target_inputs["per_drug"]) == {"drug1", "drug2"}
    assert payload["drug1_targets"]["resolution"]["cache"]["schema_fixture_verified"] is True

    batch_pipe = strict_pipeline(cache_dir, [sha1, sha2])
    batch_pipe.load_demo_resources()
    batch_result = batch_pipe.explain_batch([{
        "drug1_name": "ethanol",
        "drug2_name": "ethylamine",
        "drug1_smiles": "CCO",
        "drug2_smiles": "CCN",
        "predicted_synergy": 4.2,
        "advanced": False,
    }])[0]
    batch_credibility = batch_result.to_dict()["credibility"]
    assert batch_credibility["score_profile"]["name"] == "batch_BH_FDR"
    assert batch_credibility["score_profile"]["batch_FDR_applied"] is True
    assert "base_score" in batch_credibility
    assert "batch_adjusted_score" in batch_credibility

    cache = PlatoPredictionCache(cache_dir)
    payload_variant = plato_payload("CCO", [{
        "target_symbol": "BRAF",
        "target_name": "BRAF",
        "organism": "Homo sapiens",
        "taxid": "9606",
        "target_type": "single protein",
        "score": 10.0,
        "reliable": "yes",
    }])
    allow_test_schema(payload_variant)
    changed_source = cache_dir / "changed.json"
    write_payload(changed_source, payload_variant)
    expect_error(
        "refusing silent overwrite",
        lambda: cache.import_file(
            smiles="CCO",
            result_file=changed_source,
            drug_name="ethanol",
            provider_release="r37",
            confirmed_real_export=True,
        ),
    )

    _, metadata = cache.load("CCO")
    raw_blob = cache_dir / metadata["source"]["raw_blob_relative_path"]
    raw_blob.write_text(json.dumps({"provider": "PLATO"}), encoding="utf-8")
    expect_error(
        "raw_export_blob",
        lambda: strict_pipeline(cache_dir, [sha1, sha2]).explain_combination(
            "ethanol", "ethylamine", "CCO", "CCN",
        ),
    )


def test_invalid_smiles_fail_closed(root: Path) -> None:
    pipe = strict_pipeline(root / "empty-cache-d", [])
    expect_error(
        "RDKit-validated",
        lambda: pipe.explain_combination("bad", "drug-b", "not-a-smiles", "CCN"),
    )


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="interpreterplato-contract-") as temporary:
        root = Path(temporary)
        test_target_parsing_filters()
        fixture = (
            Path(__file__).resolve().parents[1]
            / "data" / "raw" / "plato" / "r37" / "golden_fixture"
            / "gemcitabine_target_fishing_response_2026-07-16.json"
        )
        if fixture.is_file():
            test_real_plato_r37_fixture_contract()
        else:
            print("external PLATO fixture absent: fixture contract skipped")
        test_strict_api_rejects_local_and_mutable_modes(root)
        test_unverified_and_untrusted_plato_exports_fail(root)
        test_strict_positive_fixture_and_integrity(root)
        test_invalid_smiles_fail_closed(root)
    print("strict PLATO-only contract tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
