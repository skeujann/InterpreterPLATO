"""Structure grouping changes lookup, never historical request provenance."""
from copy import deepcopy
import json

import pytest

from InterpreterPLATO import ExplanationPipeline, development_inputs as dev
from InterpreterPLATO.engineering import file_sha
from InterpreterPLATO.targets import PlatoPredictionCache, stereo_agnostic_smiles, inspect_plato_payload
from InterpreterPLATO.scripts.test_support import archive_verified_plato_cache, strict_pipeline, DEFAULT_TARGETS

RIGHT = "F[C@H](Cl)Br"
LEFT = "F[C@@H](Cl)Br"
UNSPECIFIED = "FC(Cl)Br"


@pytest.mark.parametrize("a,b", [(RIGHT, LEFT), (RIGHT, UNSPECIFIED), ("F/C=C/F", "F/C=C\\F")])
def test_stereo_is_removed(a, b):
    assert stereo_agnostic_smiles(a) == stereo_agnostic_smiles(b)


@pytest.mark.parametrize("a,b", [("CCO", "COC"), ("[13CH3]CO", "CCO"),
                                 ("C[NH3+]", "CN"), ("CCO.Cl", "CCO")])
def test_other_chemical_differences_remain(a, b):
    assert stereo_agnostic_smiles(a) != stereo_agnostic_smiles(b)


def test_reuse_validates_requested_group_and_retains_original_query(tmp_path):
    cache = tmp_path / "cache"
    cache.mkdir()
    raw_hash = archive_verified_plato_cache(cache, RIGHT, "right")
    before = {path: path.read_bytes() for path in cache.rglob("*") if path.is_file()}
    pipe = strict_pipeline(cache, [raw_hash])
    assert pipe.config["plato"]["ignore_stereochemistry"] is True
    for smiles in (LEFT, RIGHT, UNSPECIFIED):
        fused, resolution = pipe.targets.resolve("compound", smiles)
        pipe._validate_plato_only_resolution("compound", smiles, fused, resolution)
        metadata = resolution["cache"]
        assert metadata["query"]["submitted_smiles"] == RIGHT
        assert metadata["structure_matching"]["requested_structure"]["submitted_smiles"] == smiles
        assert len(fused) == 1
    assert before == {path: path.read_bytes() for path in cache.rglob("*") if path.is_file()}
    # The underlying archive reader keeps exact semantics.
    assert PlatoPredictionCache(cache).load(LEFT)[1]["status"] == "cache_miss"


def test_group_selection_is_input_independent_and_target_order_independent(tmp_path):
    a = dict(DEFAULT_TARGETS[0])
    b = {**a, "target_symbol": "TP53", "target_name": "TP53"}
    archive_verified_plato_cache(tmp_path, RIGHT, "right", (a, b))
    archive_verified_plato_cache(tmp_path, LEFT, "left", (b, a))
    cache = PlatoPredictionCache(tmp_path)
    matches = [cache.load_ignoring_stereochemistry(smiles)[1] for smiles in (RIGHT, LEFT, UNSPECIFIED)]
    assert len({item["cache_path"] for item in matches}) == 1
    assert all(item["structure_matching"]["n_cached_structures"] == 2 for item in matches)


def test_conflicting_prediction_scores_are_not_arbitrarily_merged(tmp_path):
    archive_verified_plato_cache(tmp_path, RIGHT, "right")
    archive_verified_plato_cache(tmp_path, LEFT, "left", [{**DEFAULT_TARGETS[0], "score": 1.}])
    cache = PlatoPredictionCache(tmp_path)
    for smiles in (RIGHT, LEFT, UNSPECIFIED):
        targets, metadata = cache.load_ignoring_stereochemistry(smiles)
        assert not targets
        assert metadata["status"] == "cache_error_stereoisomer_prediction_conflict"
        assert len(metadata["structure_matching"]["members"]) == 2


def test_other_group_cannot_be_reused_and_exact_mode_can_be_restored(tmp_path):
    digest = archive_verified_plato_cache(tmp_path, RIGHT, "right")
    assert PlatoPredictionCache(tmp_path).load_ignoring_stereochemistry("CCO")[1]["status"] == "cache_miss"
    pipe = strict_pipeline(tmp_path, [digest], plato={"cache_dir": str(tmp_path),
                           "trusted_export_sha256": [digest], "ignore_stereochemistry": False})
    fused, resolution = pipe.targets.resolve("compound", LEFT)
    with pytest.raises(ValueError, match="cache_miss"):
        pipe._validate_plato_only_resolution("compound", LEFT, fused, resolution)


def test_bad_group_member_is_not_hidden_by_an_intact_exact_entry(tmp_path):
    archive_verified_plato_cache(tmp_path, RIGHT, "right")
    archive_verified_plato_cache(tmp_path, LEFT, "left")
    cache = PlatoPredictionCache(tmp_path)
    _, meta = cache.load(LEFT)
    (tmp_path / meta["source"]["raw_blob_relative_path"]).write_bytes(b"[]")
    targets, metadata = cache.load_ignoring_stereochemistry(RIGHT)
    assert not targets and metadata["status"] == "cache_error_structure_group_integrity"


def test_real_input_policy_keeps_original_post_binding_on_cross_stereo_reuse(tmp_path, monkeypatch):
    payload = [{"chemblId": "CHEMBL203", "nome": "Epidermal growth factor receptor:Homo sapiens",
                "score": 10., "simili": [], "tid": 1, "validato": True}]
    raw = tmp_path / "raw.json"
    raw.write_text(json.dumps(payload), encoding="utf-8")
    monkeypatch.setattr(dev, "OBSERVED_SCHEMAS", frozenset({inspect_plato_payload(payload)["schema_fingerprint"]}))
    cache = tmp_path / "cache"
    evidence = {"endpoint": "https://prometheus.farmacia.uniba.it/plato/pred.html", "http_method": "POST",
                "response_sha256": file_sha(raw), "submitted_smiles": RIGHT,
                "form": {"tipoout": "json", "tipopred": "mussel"}}
    PlatoPredictionCache(cache).import_file(RIGHT, raw, "right", "r37", query_binding_evidence=evidence)
    pins = {str(path.resolve()): file_sha(path) for path in cache.rglob("*") if path.is_file()}
    pipe = dev.DevelopmentTestPipeline({"run_mode": "test", "generated_from_test_fixture": True,
                                       "plato": {"cache_dir": str(cache)},
                                       "test_source_review_exception": dev.test_policy(pins)})
    fused, resolution = pipe.targets.resolve("compound", LEFT)
    pipe._validate_plato_only_resolution("compound", LEFT, fused, resolution)
    original = resolution["cache"]["query"]["single_query_evidence"]
    assert original["submitted_smiles"] == RIGHT
    bad = deepcopy(resolution)
    bad["cache"]["query"]["single_query_evidence"]["submitted_smiles"] = LEFT
    with pytest.raises(ValueError, match="test_official_request_binding_missing"):
        pipe._validate_plato_only_resolution("compound", LEFT, fused, bad)


def test_stereo_policy_requires_boolean_and_cannot_drift_after_construction():
    with pytest.raises(TypeError, match="ignore_stereochemistry"):
        ExplanationPipeline({"plato": {"ignore_stereochemistry": "false"}})
    pipe = ExplanationPipeline()
    pipe.config["plato"]["ignore_stereochemistry"] = False
    with pytest.raises(RuntimeError, match="stereochemistry policies differ"):
        pipe._assert_strict_plato_only_runtime()
