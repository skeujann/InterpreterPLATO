"""Offline software acceptance only; fixtures never validate real mechanisms."""
from copy import deepcopy
from pathlib import Path

import networkx as nx
import pytest

from InterpreterPLATO import ExplanationPipeline
from InterpreterPLATO import engineering as eng
from InterpreterPLATO import targets
from InterpreterPLATO.experiments import apply_control, validate_experiment
from InterpreterPLATO.research import _code_snapshot, _dependency_versions, write_json_with_sha256
from InterpreterPLATO.scripts.run_engineering import main
from InterpreterPLATO.scripts.test_support import archive_verified_plato_cache, strict_pipeline


@pytest.fixture
def fixture_suite(tmp_path):
    cache = tmp_path / "cache"
    cache.mkdir()
    hashes = []
    for name, smiles, genes in (("a", "CCO", ("EGFR", "TP53")), ("b", "CCN", ("RRM1", "KRAS"))):
        records = [{"target_symbol": gene, "target_name": gene, "organism": "Homo sapiens",
                    "taxid": "9606", "target_type": "single protein", "score": 8.0 + index,
                    "reliable": "yes"} for index, gene in enumerate(genes)]
        hashes.append(archive_verified_plato_cache(cache, smiles, name, records))
    pipe = strict_pipeline(cache, hashes, pharmacology_roles={"enabled": False})
    graph = nx.Graph([("EGFR", "TP53"), ("TP53", "KRAS"), ("KRAS", "RRM1")])
    pipe.context.load_ppi(graph)
    pipe.context.load_expression({"TEST": {**{f"G{i}": float(i) for i in range(30)},
                                          "EGFR": 29., "TP53": 0., "RRM1": 25., "KRAS": 1.}},
                                 background_complete=True)
    pipe.context.load_pathways({"TEST_PATH": {"EGFR", "TP53", "KRAS", "RRM1"}})
    pipe.context.load_disease_genes({"EGFR"})
    selection = {"selection_id": "synthetic-dev", "drug_structures": {
        "raw_smiles_by_dataset_id": {"a": "CCO", "b": "CCN"}}, "selected": [{
            "candidate_id": "DEV-01", "drug_a": "a", "drug_b": "b", "cell_line_name": "TEST",
            "pred_syn_mean": 2., "role": "development_only"}]}
    cases = eng.selection_cases(selection)
    report = {"effective_config": deepcopy(pipe.config), "inventory": {},
              "source_snapshot": _code_snapshot(eng.PACKAGE), "dependencies": _dependency_versions(),
              "global_issues": [], "upstream_provenance": {"selection_id": "synthetic-dev"},
              "cases": [{"case_id": "DEV-01", "ready": True, "issues": []}]}
    return pipe, selection, cases, report


def test_template_never_approves_and_bad_case_ids_fail(fixture_suite):
    _, selection, _, _ = fixture_suite
    review = eng.review_template(selection)
    assert {item["status"] for item in review["drugs"].values()} == {"pending"}
    for invalid in ("../escape", "a/b", "CON/child", "CON", "LPT1"):
        bad = deepcopy(selection)
        bad["selected"][0]["candidate_id"] = invalid
        with pytest.raises(ValueError, match="candidate_id"):
            eng.selection_cases(bad)
    duplicate = deepcopy(selection)
    duplicate["selected"] *= 2
    with pytest.raises(ValueError, match="duplicate"):
        eng.selection_cases(duplicate)


def test_audit_pending_and_mismatch_are_blocked(fixture_suite):
    pipe, selection, cases, _ = fixture_suite
    review = eng.review_template(selection)
    result = eng.audit_cases(cases, review, pipe, True, input_policy="strict")
    assert not result[0]["ready"]
    assert {issue["code"] for issue in result[0]["issues"]} == {"chemical_identity_review_pending"}
    for item in review["drugs"].values():
        item.update(status="reviewed", reviewed_by="synthetic tester", reviewed_on="2026-10-13",
                    identity_reference="synthetic fixture only")
    assert eng.audit_cases(cases, review, pipe, True, input_policy="strict")[0]["ready"]
    review["drugs"]["a"]["reviewed_smiles"] = "CCC"
    result = eng.audit_cases(cases, review, pipe, True, input_policy="strict")
    assert not result[0]["ready"]
    assert result[0]["issues"][0]["code"] == "reviewed_structure_mismatch"


def test_review_cannot_approve_production_schema(fixture_suite, monkeypatch):
    pipe, selection, cases, _ = fixture_suite
    review = eng.review_template(selection)
    for item in review["drugs"].values():
        item.update(status="reviewed", reviewed_by="tester", reviewed_on="2026-10-13", identity_reference="fixture")
    monkeypatch.setattr(targets, "VERIFIED_PLATO_SCHEMA_FINGERPRINTS", ())
    result = eng.audit_cases(cases, review, pipe, True, input_policy="strict")
    assert not result[0]["ready"]
    assert any(issue["code"] == "plato_input_rejected" for issue in result[0]["issues"])


def test_test_mode_accepts_pending_and_keeps_identity_notes(fixture_suite):
    pipe, selection, cases, _ = fixture_suite
    reviews = eng.review_template(selection)
    reviews["drugs"]["a"]["notes"] = "Reference has different stereochemistry"
    row = eng.audit_cases(cases, reviews, pipe, True)[0]
    assert row["ready"] and not row["issues"]
    assert row["inputs"]["drug1"]["observations"][0]["code"] == "chemical_identity_review_pending"
    assert row["inputs"]["drug1"]["identity_reference_record"]["notes"] == reviews["drugs"]["a"]["notes"]
    assert not row["inputs"]["drug1"]["manual_acceptance_required"]
    assert eng.audit_cases(cases, {"drugs": {}}, pipe, True)[0]["ready"]


def test_default_cli_needs_no_reviews_or_production_approval(fixture_suite, tmp_path, monkeypatch):
    pipe, selection, _, _ = fixture_suite
    selection_path = tmp_path / "selection.json"
    write_json_with_sha256(selection_path, selection)
    monkeypatch.setattr(eng, "resource_check", lambda _: ([], {}))
    monkeypatch.setattr(pipe, "load_production_resources", lambda _: {})
    monkeypatch.setattr(eng, "DevelopmentTestPipeline", lambda config: pipe)
    # The fixture has its synthetic allowlist; hide it until the actual parser runs.
    allowed = targets.VERIFIED_PLATO_SCHEMA_FINGERPRINTS
    monkeypatch.setattr(targets, "VERIFIED_PLATO_SCHEMA_FINGERPRINTS", ())
    original = eng.audit_cases
    def audit(*args):
        monkeypatch.setattr(targets, "VERIFIED_PLATO_SCHEMA_FINGERPRINTS", allowed)
        return original(*args)
    monkeypatch.setattr(eng, "audit_cases", audit)
    out = tmp_path / "cli"
    assert main(["preflight", "--selection", str(selection_path), "--plato-cache", str(tmp_path),
                 "--output-dir", str(out)]) == 0
    report = eng.read_json(out / "preflight.json")
    assert report["input_policy"] == "test" and report["all_ready"]
    assert not report["manual_acceptance_required"] and not report["global_issues"]


def test_interventions_recompute_and_preserve_raw_resources(fixture_suite):
    pipe, _, cases, _ = fixture_suite
    args = cases[0]["input"]
    full = pipe.explain_combination(**args)
    context_before = deepcopy(pipe.context.expression)
    pipe.config["engineering_experiment"] = {"variant": "no_expression_ranking"}
    noexpr = pipe.explain_combination(**args)
    assert any(item.effective_weights.get("expression", 0) > 0 for item in full.drug1_targets.all_targets)
    assert all(item.effective_weights.get("expression", 0) == 0 for item in noexpr.drug1_targets.all_targets)
    assert all(item.expression_measured for item in noexpr.drug1_targets.all_targets)
    pipe.config["engineering_experiment"] = {"variant": "raw_plato_order"}
    raw = pipe.explain_combination(**args)
    assert raw.drug1_targets.all_targets[0].target_id == "TP53"
    assert raw.drug1_targets.to_dict()["top_targets_semantics"].startswith("provider_score")
    assert {item.target_id: item.final_score for item in raw.drug1_targets.all_targets} == {
        item.target_id: item.final_score for item in full.drug1_targets.all_targets}
    pipe.config["engineering_experiment"] = {"variant": "drop_top_target"}
    dropped = pipe.explain_combination(**args)
    assert len(dropped.drug1_targets.all_targets) == 1
    assert dropped.drug1_targets.all_targets[0].target_id == full.drug1_targets.all_targets[1].target_id
    audit = dropped.metadata["audit"]["engineering_experiment"]
    assert audit["targets"]["drug1"]["removed"] == [full.drug1_targets.all_targets[0].target_id]
    assert full.mechanism.metrics != dropped.mechanism.metrics
    assert pipe.context.expression == context_before


def test_single_target_deletion_is_refused(fixture_suite):
    pipe, _, cases, _ = fixture_suite
    full = pipe.explain_combination(**cases[0]["input"])
    with pytest.raises(ValueError, match=">=2"):
        apply_control(full.drug1_targets.all_targets[:1], full.drug2_targets.all_targets,
                      {"engineering_experiment": {"variant": "drop_top_target"}})


def test_recorded_score_does_not_change_biology_and_repeat_is_deterministic(fixture_suite):
    pipe, _, cases, _ = fixture_suite
    args = cases[0]["input"]
    first = pipe.explain_combination(**args)
    repeat = pipe.explain_combination(**args)
    changed = pipe.explain_combination(**{**args, "predicted_synergy": 999.})
    for result in (repeat, changed):
        assert first.drug1_targets.to_dict() == result.drug1_targets.to_dict()
        assert first.mechanism.metrics == result.mechanism.metrics
        assert first.modules.to_dict() == result.modules.to_dict()
    assert changed.predicted_synergy == 999.


def test_top_k_preserves_statistics_and_swap_preserves_side_identity(fixture_suite):
    pipe, _, cases, _ = fixture_suite
    args = cases[0]["input"]
    first = pipe.explain_combination(**args)
    pipe.config["top_k_targets"] = 1
    display = pipe.explain_combination(**args)
    assert len(display.drug1_targets.top_targets) == 1
    assert display.mechanism.metrics["statistical_target_sets"] == first.mechanism.metrics["statistical_target_sets"]
    for key in ("d_A_obs", "d_B_obs", "d_AB", "s_AB", "p_A", "p_B"):
        assert display.mechanism.metrics["cheng"].get(key) == first.mechanism.metrics["cheng"].get(key)
    swapped = pipe.explain_combination(**{**args,
        "drug1_name": args["drug2_name"], "drug1_smiles": args["drug2_smiles"],
        "drug2_name": args["drug1_name"], "drug2_smiles": args["drug1_smiles"]})
    assert swapped.drug1_targets.to_dict() == display.drug2_targets.to_dict()
    assert swapped.drug2_targets.to_dict() == display.drug1_targets.to_dict()
    assert swapped.mechanism.metrics["cheng"].get("d_AB") == display.mechanism.metrics["cheng"].get("d_AB")


def test_runner_resume_and_hash_tampering(fixture_suite, tmp_path):
    pipe, _, cases, report = fixture_suite
    root = tmp_path / "run"
    result = eng.run_suite(root, pipe, cases, report, ("full", "no_expression_ranking"))
    assert result["complete"] and result["n_planned"] == 2
    assert result["biological_validation_performed"] is False
    assert result["comparisons_to_full"][0]["statistical_target_sets_unchanged"]
    assert result["comparisons_to_full"][0]["status"] == "descriptive_comparison"
    eng.run_suite(root, pipe, cases, report, ("full", "no_expression_ranking"), resume=True)
    assert len(list(root.rglob("attempt_*"))) == 2
    payload = eng.read_json(next(root.rglob("explanation.json")))
    assert payload["generated_from_test_fixture"] is True
    assert payload["audit"]["engineering_case"]["case_id"] == "DEV-01"
    artifact = next(root.rglob("explanation.json.report.md"))
    artifact.write_text("tampered", encoding="utf-8")
    assert eng.summarize(root)["counts"] == {"corrupt": 1, "complete": 1}
    with pytest.raises(ValueError, match="corrupt"):
        eng.run_suite(root, pipe, cases, report, ("full", "no_expression_ranking"), resume=True)


def test_all_blocked_denominator_and_changed_identity(fixture_suite, tmp_path):
    pipe, _, cases, report = fixture_suite
    cases = [dict(cases[0], case_id=f"DEV-{i}") for i in range(6)]
    report["cases"] = [{"case_id": case["case_id"], "ready": False,
                        "issues": [eng.issue("pending", "test", "retain denominator")]} for case in cases]
    result = eng.run_suite(tmp_path / "run", pipe, cases, report)
    assert result["counts"] == {"blocked": 6}
    assert result["n_planned"] == 6 and not result["complete"]
    report["effective_config"]["top_k_targets"] = 17
    with pytest.raises(ValueError, match="changed"):
        eng.run_suite(tmp_path / "run", pipe, cases, report, resume=True)


def test_failure_does_not_discard_later_cases(fixture_suite, tmp_path, monkeypatch):
    pipe, _, cases, report = fixture_suite
    original = pipe.explain_combination
    calls = []
    def sometimes(**kwargs):
        calls.append(kwargs)
        if len(calls) == 1:
            raise RuntimeError("simulated failure")
        return original(**kwargs)
    monkeypatch.setattr(pipe, "explain_combination", sometimes)
    result = eng.run_suite(tmp_path / "run", pipe, cases, report, ("full", "no_expression_ranking"))
    assert result["counts"] == {"failed": 1, "complete": 1}
    assert not result["complete"]


def test_changed_input_after_preflight_and_output_lock(fixture_suite, tmp_path):
    pipe, _, cases, report = fixture_suite
    path = tmp_path / "input.json"
    path.write_text("before", encoding="utf-8")
    report["inventory"][str(path)] = eng.file_sha(path)
    path.write_text("after", encoding="utf-8")
    with pytest.raises(ValueError, match="Input changed"):
        eng.run_suite(tmp_path / "run", pipe, cases, report)
    report["inventory"] = {}
    root = tmp_path / "run"
    root.mkdir()
    (root / ".runner.lock").write_text("other", encoding="utf-8")
    with pytest.raises(FileExistsError):
        eng.run_suite(root, pipe, cases, report)
    assert (root / ".runner.lock").read_text() == "other"


def test_bad_experiment_and_strict_json(tmp_path):
    for invalid in ({"variant": "anything"}, {"variant": "full", "bypass": True}, "full"):
        with pytest.raises(ValueError):
            validate_experiment(invalid)
    path = tmp_path / "bad.json"
    for text in ('{"x": NaN}', '{"x": 1, "x": 2}'):
        path.write_text(text, encoding="utf-8")
        with pytest.raises(ValueError):
            eng.read_json(path)
    with pytest.raises(ValueError):
        ExplanationPipeline(config={"engineering_experiment": {"variant": "no_gate"}})


def test_cli_prepare_preserves_existing_reviews(fixture_suite, tmp_path):
    _, selection, _, _ = fixture_suite
    path = tmp_path / "selection.json"
    review = tmp_path / "reviews.json"
    write_json_with_sha256(path, selection)
    argv = ["prepare", "--selection", str(path), "--output", str(review)]
    assert main(argv) == 0
    before = review.read_bytes()
    assert main(argv) == 2
    assert review.read_bytes() == before
