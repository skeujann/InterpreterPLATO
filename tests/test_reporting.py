"""Protect scientific distinctions at the human-readable output boundary."""

from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

import pytest

from InterpreterPLATO.reporting import build_explanation_layers, render_payload, write_readable_report
from InterpreterPLATO.scripts.render_explanation import main
from InterpreterPLATO.research import stable_sha256


def payload():
    return {
        "drug1_name": "a", "drug2_name": "b", "cell_line": "cell", "run_mode": "test",
        "drug1_targets": {"targets": [{"target_id": "X", "evidence_status": "predicted_only"}]},
        "drug2_targets": {"targets": [{"target_id": "X", "evidence_status": "predicted_only"}]},
        "mechanism_type": "network_proximal_targets", "mechanism_score": .5,
        "mechanism": {"shared_targets": [], "shared_pathways": [], "metrics": {}},
        "local_evidence": {
            "cross_drug_paths": [{"drug": "a", "destination_drug": "b", "path": ["X"], "length": 0}],
            "top_paths": [{"drug": "a", "path": ["X"], "length": 0}],
        },
        "network_connectivity": {"cross_drug": {"candidate_path_count": 3}},
        "network_statistics": {
            "valid": True, "quality_status": "adequate", "near_A": True, "near_B": False,
            "p_A": .01, "p_B": .4, "p_A_adjusted": .02, "p_B_adjusted": .8,
        },
    }


def test_full_shared_targets_and_zero_hop_are_not_lost_to_top_k_classification():
    data = payload()
    before = deepcopy(data)
    view = build_explanation_layers(data)
    findings = {item["id"]: item for item in view["findings"]}
    assert findings["shared_targets"]["targets"] == ["X"]
    assert findings["shared_targets"]["classification_subset"] == []
    assert findings["cross_drug"]["counts"] == {"displayed": 1, "nonzero": 0, "zero_hop": 1}
    assert findings["disease_paths"]["counts"]["zero_hop"] == 1
    text = render_payload(data)
    assert "共享候选靶点零跳" in text and "疾病集合成员零跳" in text
    assert "0 条非零跳双药 PPI" in text
    assert data == before


def test_structure_policy_is_displayed_and_single_pair_does_not_claim_batch_correction():
    data = payload()
    data["audit"] = {"plato_inputs": {"drug1": {"structure_matching": {
        "policy": "ignore_stereochemistry", "cross_stereochemistry_reuse": True}}}}
    text = render_payload(data)
    assert "忽略立体标记" in text and "复用同一连接结构" in text
    assert "批次校正后的疾病邻近结论" not in text
    assert "实际校正方法" in text


def test_nonzero_routes_and_display_cap_do_not_claim_all_routes_found():
    data = payload()
    data["local_evidence"]["cross_drug_paths"].append({"drug": "a", "destination_drug": "b", "path": ["A", "B", "C"]})
    text = render_payload(data)
    assert "1 条非零跳双药 PPI" in text and "2 条 PPI 边" in text
    assert "3 条候选记录" in text
    data["local_evidence"]["cross_drug_paths"] = []
    text = render_payload(data)
    assert "当前展示名额为 0" in text
    assert "不连通" not in text


def test_disabled_search_is_not_reported_as_disconnected():
    data = payload()
    data["local_evidence"] = {}
    data["network_connectivity"] = {"status": "disabled"}
    text = render_payload(data)
    assert "路径搜索未启用" in text
    assert "不连通" not in text


def test_batch_q_overrules_local_near_in_presented_statistics():
    data = payload()
    data["mechanism_type"] = "single_exposure"
    data["network_statistics"].update(batch_q_A=.12, batch_q_B=.8,
        batch_near_A=False, batch_near_B=False, batch_test_A_evaluable=True, batch_test_B_evaluable=True)
    sides = build_explanation_layers(data)["statistics"]["sides"]
    assert sides[0]["near"] is False
    assert "未达到" in sides[0]["status"]
    assert "批次 q=0.12" in render_payload(data)


def test_unevaluable_placeholder_q_is_not_an_observed_negative_result():
    data = payload()
    data["network_statistics"].update(quality_status="limited", batch_q_A=1., batch_q_B=1.,
        batch_near_A=False, batch_near_B=False, batch_test_A_evaluable=False, batch_test_B_evaluable=False,
        warnings=["insufficient_valid_randomizations"],
        randomization_audit={"drug1": {"requested": 2000, "accepted": 1403,
            "rejection_reasons": {"incomplete_reachability": 597}}})
    text = render_payload(data)
    assert "70.15%" in text and "597 次" in text
    assert "占位值" in text and "批次 q=1" not in text
    assert "未达到当前疾病邻近显著性标准" not in text
    assert "本次未作疾病邻近显著性判定" in text
    assert "机制支持不足" not in text


def test_shared_pathway_annotations_are_not_promoted_to_enrichment():
    data = payload()
    data["mechanism"]["evidence"] = [{"type": "descriptive_shared_pathways_not_enriched", "pathways": ["P"]}]
    view = build_explanation_layers(data)
    item = next(f for f in view["findings"] if f["id"] == "pathways")
    assert item["significant_shared"] == [] and item["descriptive_shared"] == ["P"]
    assert "另有 1 项共同通路注释" in render_payload(data)


def test_explicit_failure_is_not_fabricated_into_an_explanation():
    data = {"drug1": "a", "drug2": "b", "error": "No eligible human protein targets",
            "eligible_targets": {"drug1": 8, "drug2": 0}}
    view = build_explanation_layers(data)
    assert view["status"] == "not_generated" and view["findings"] == []
    text = render_payload(data)
    assert "No eligible human protein targets" in text and "靶点 0 个" in text
    assert "疾病邻近统计" not in text


def test_saved_report_hash_and_cli_leave_inputs_untouched(tmp_path):
    data = payload()
    source = tmp_path / "result.json"
    source.write_text(json.dumps(data), encoding="utf-8")
    before = source.read_bytes()
    assert main(["--input", str(source)]) == 0
    report = tmp_path / "result.json.report.md"
    assert report.exists()
    assert source.read_bytes() == before
    assert report.with_name(report.name + ".sha256").read_text().split()[0] == sha256(report.read_bytes()).hexdigest()
    with pytest.raises(SystemExit):
        main(["--input", str(source), "--output", str(source)])
    assert source.read_bytes() == before


def test_live_api_and_saved_report_share_the_view_and_preserve_statistics(tmp_path):
    from InterpreterPLATO.scripts.test_support import archive_verified_plato_cache, strict_pipeline
    cache = tmp_path / "cache"
    cache.mkdir()
    hashes = [archive_verified_plato_cache(cache, "CCO", "a"), archive_verified_plato_cache(cache, "CCN", "b")]
    pipe = strict_pipeline(cache, hashes)
    pipe.load_demo_resources()
    result = pipe.explain_combination("a", "b", "CCO", "CCN", advanced=False)
    before = deepcopy(result.mechanism.to_dict())
    data = result.to_dict()
    assert result.mechanism.to_dict() == before
    assert "explanation_layers" in data
    assert data["network_statistics"] == before["metrics"]["cheng"]
    digest = data["audit"].pop("result_payload_sha256")
    assert digest == stable_sha256(data)
    destination = tmp_path / "new.json"
    record = pipe.save_result(result, destination)
    saved = json.loads(destination.read_text(encoding="utf-8"))
    assert saved["integrated_explanation"]["summary"] == destination.with_name("new.json.report.md").read_text(encoding="utf-8")
    assert record["report_path"].endswith("new.json.report.md")
    assert result.mechanism.to_dict() == before


def test_batch_report_preserves_each_pair_and_single_sided_identity(tmp_path):
    data = payload()
    data["network_statistics"].update(batch_q_A=.042, batch_q_B=.989,
        batch_near_A=True, batch_near_B=False, batch_test_A_evaluable=True, batch_test_B_evaluable=True)
    output = tmp_path / "batch.md"
    write_readable_report(output, [data, {"drug1": "c", "drug2": "d", "error": "missing targets"}])
    text = output.read_text(encoding="utf-8")
    assert "a：显著邻近" in text and "b：未达到" in text
    assert "批次 q=0.042" in text and "# c + d" in text


def test_standalone_cli_runs_from_an_unrelated_directory(tmp_path):
    source = tmp_path / "saved.json"
    source.write_text(json.dumps(payload()), encoding="utf-8")
    script = Path(__file__).resolve().parents[1] / "scripts" / "render_explanation.py"
    process = subprocess.run([sys.executable, "-B", str(script), "--input", str(source)],
                             cwd=tmp_path, capture_output=True, text=True)
    assert process.returncode == 0, process.stderr
    assert (tmp_path / "saved.json.report.md").is_file()
