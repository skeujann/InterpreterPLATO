"""Offline development-case audit and sequential, identity-checked execution.

No network queries, trust grants, selection changes, or biological validation.
The denominator is the frozen case x variant plan, including blocked cases.
"""
from __future__ import annotations

from collections import Counter
from copy import deepcopy
from datetime import date, datetime, timezone
from hashlib import sha256
from pathlib import Path
import json
import math
import os
import platform
import re
import sys
import time

from . import targets
from .experiments import VARIANTS
from .interpreter import ExplanationPipeline
from .development_inputs import DevelopmentTestPipeline, test_policy
from .research import (_code_snapshot, _dependency_versions, plato_resolution_audit,
                       stable_sha256, write_json_with_sha256)

CONTRACT = "engineering_validation_v1"
PACKAGE = Path(__file__).resolve().parent
SAFE_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]{0,79}\Z")
RESERVED_IDS = {"con", "prn", "aux", "nul", *(f"com{i}" for i in range(1, 10)),
                *(f"lpt{i}" for i in range(1, 10))}
CONTROLS = (*VARIANTS, "seed_43", "seed_44", "top_k_5", "top_k_20")


def utc_now():
    return datetime.now(timezone.utc).isoformat()


def file_sha(path):
    digest = sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    def invalid(value):
        raise ValueError(f"Non-finite JSON constant: {value}")
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"Duplicate JSON key: {key}")
            result[key] = value
        return result
    return json.loads(Path(path).read_text(encoding="utf-8-sig"),
                      parse_constant=invalid, object_pairs_hook=unique)


def selection_cases(selection):
    """Preserve original candidate IDs/order and upstream lineage, no reselecting."""
    if not isinstance(selection, dict) or not selection.get("selection_id"):
        raise ValueError("A named ScoreRI selection is required")
    structures = selection["drug_structures"]["raw_smiles_by_dataset_id"]
    rows = selection["selected"]
    if not isinstance(rows, list) or not rows:
        raise ValueError("selection.selected must be a nonempty list")
    cases, seen = [], set()
    for row in rows:
        case_id = row["candidate_id"]
        if (not isinstance(case_id, str) or not SAFE_ID.fullmatch(case_id)
                or case_id.lower() in RESERVED_IDS or case_id.casefold() in seen):
            raise ValueError(f"Unsafe or duplicate candidate_id: {case_id!r}")
        seen.add(case_id.casefold())
        score = row["pred_syn_mean"]
        if isinstance(score, bool) or not isinstance(score, (int, float)) or not math.isfinite(score):
            raise ValueError(f"{case_id}: predicted synergy must be finite")
        if row.get("role") != "development_only":
            raise ValueError("This runner accepts development_only cases, not independent holdout claims")
        cases.append({"case_id": case_id, "role": "development_only", "upstream": row,
                      "input": {"drug1_name": row["drug_a"], "drug2_name": row["drug_b"],
                                "drug1_smiles": structures.get(row["drug_a"], ""),
                                "drug2_smiles": structures.get(row["drug_b"], ""),
                                "cell_line": row["cell_line_name"],
                                "predicted_synergy": score}})
    return cases


def review_template(selection):
    cases = selection_cases(selection)
    drugs = {}
    for case in cases:
        for side in ("drug1", "drug2"):
            name = case["input"][side + "_name"]
            drugs[name] = {"status": "pending", "reviewed_smiles": case["input"][side + "_smiles"],
                           "reviewed_by": "", "reviewed_on": "", "identity_reference": "",
                           "notes": "核实名称、盐型及立体化学；此记录不批准 PLATO schema 或原始响应。"}
    return {"contract": CONTRACT, "selection_id": selection["selection_id"],
            "selection_payload_sha256": stable_sha256(selection), "drugs": drugs}


def issue(code, detail, action):
    return {"code": code, "detail": str(detail), "action": action}


def resource_check(context_path):
    """Cheap all-file diagnostics before the authoritative production loader."""
    findings, inventory = [], {}
    try:
        config = read_json(context_path)
        inventory[str(context_path.resolve())] = file_sha(context_path)
        resources = config["resources"]
        for name in ("ppi", "expression", "pathways", "disease", "target_mapping"):
            record = resources.get(name, {})
            path = (context_path.parent / record.get("path", "")).resolve()
            if not record or not path.is_file():
                findings.append(issue("resource_missing", f"{name}: {path}",
                                      "恢复锁定版本资源；不要删除资源来模拟消融。"))
                continue
            inventory[str(path)] = file_sha(path)
            if inventory[str(path)] != record.get("sha256"):
                findings.append(issue("resource_hash_mismatch", name, "恢复锁定文件并核对来源，不要改哈希凑通过。"))
        manifest = config.get("runtime_manifest")
        if not manifest:
            findings.append(issue("runtime_manifest_missing", context_path, "配置受版本管理的 runtime manifest。"))
        else:
            path = (context_path.parent / manifest).resolve()
            if path.is_file():
                inventory[str(path)] = file_sha(path)
            else:
                findings.append(issue("runtime_manifest_missing", path, "恢复 runtime.lock.json。"))
    except (OSError, ValueError, TypeError, KeyError) as exc:
        findings.append(issue("context_config_invalid", exc, "提供有效的 production_context.json。"))
    return findings, inventory


def audit_cases(cases, reviews, pipeline, context_ready, input_policy="test"):
    if input_policy not in ("test", "strict"):
        raise ValueError("input_policy must be test or strict")
    results, drug_audits = [], {}
    for case in cases:
        findings, resolved = [], {}
        for side in ("drug1", "drug2"):
            name, smiles = (case["input"][side + suffix] for suffix in ("_name", "_smiles"))
            key = (name, smiles)
            if key not in drug_audits:
                errors, observations = [], []
                structure = targets.standardize_smiles(smiles)
                if not structure.get("inchikey") or not str(structure.get("standardization_status", "")).startswith("rdkit_canonical"):
                    errors.append(issue("invalid_structure_or_rdkit_missing", name, "核对 SMILES，并使用含 RDKit 的项目环境。"))
                review = reviews.get("drugs", {}).get(name, {})
                approved = (review.get("status") == "reviewed" and all(
                    isinstance(review.get(field), str) and review[field].strip()
                    for field in ("reviewed_by", "reviewed_on", "identity_reference", "reviewed_smiles")))
                if approved:
                    try:
                        date.fromisoformat(review["reviewed_on"])
                    except ValueError:
                        approved = False
                review_findings = errors if input_policy == "strict" else observations
                if not approved:
                    review_findings.append(issue("chemical_identity_review_pending", name,
                                        "测试方式不要求填写签字或审批记录；按原提交结构运行。" if input_policy == "test"
                                        else "核实候选名与结构、盐型/立体化学，填写来源和审核记录。"))
                elif not targets.structure_identity_matches(targets.standardize_smiles(review["reviewed_smiles"]),
                        {**structure, "ignore_stereochemistry": bool(pipeline and pipeline.config["plato"].get("ignore_stereochemistry"))}):
                    review_findings.append(issue("reviewed_structure_mismatch", name, "保留原提交结构与参考结构的差异记录。"))
                audit = {"structure": structure, "observations": observations,
                         "identity_reference_record": deepcopy(review),
                         "manual_acceptance_required": input_policy == "strict"}
                if pipeline is not None:
                    try:
                        fused, resolution = pipeline.targets.resolve(name, smiles)
                        audit["plato"] = plato_resolution_audit(resolution)
                        audit["cache_status"] = resolution.get("cache", {}).get("status")
                        audit["eligible_target_count"] = len(fused)
                        pipeline._validate_plato_only_resolution(name, smiles, fused, resolution)
                        if not pipeline._filter_plato_applicability(name, fused, resolution):
                            raise ValueError("No targets remain after applicability filtering")
                    except (OSError, ValueError, TypeError, RuntimeError, KeyError) as exc:
                        errors.append(issue("plato_input_rejected", exc,
                                            "检查具体的数据缺失、结构绑定、文件完整性或可用靶点问题。"))
                drug_audits[key] = (errors, audit)
            errors, audit = drug_audits[key]
            findings.extend({**error, "drug": name} for error in errors)
            resolved[side] = audit
        if context_ready:
            used, profile, fallback = pipeline.context.expression_profile(case["input"]["cell_line"])
            if not used or not profile or fallback:
                findings.append(issue("cell_context_unavailable", case["input"]["cell_line"],
                                      "检查细胞系名称映射，不允许 default 表达替代。"))
        results.append({"case_id": case["case_id"], "issues": findings, "inputs": resolved,
                        "ready": context_ready and not findings})
    return results


def preflight(selection_path, reviews_path, context_path, cache_path, trusted_hashes=(), input_policy="test",
              ignore_stereochemistry=True):
    if input_policy not in ("test", "strict"):
        raise ValueError("input_policy must be test or strict")
    selection_path = Path(selection_path).resolve()
    reviews_path = Path(reviews_path).resolve() if reviews_path else None
    context_path, cache_path = Path(context_path).resolve(), Path(cache_path).resolve()
    selection = read_json(selection_path)
    reviews = read_json(reviews_path) if reviews_path else {"drugs": {}}
    if (not isinstance(reviews, dict) or not isinstance(reviews.get("drugs"), dict)
            or any(not isinstance(value, dict) for value in reviews["drugs"].values())):
        raise ValueError("Source reviews must contain a drugs mapping of review records")
    cases = selection_cases(selection)
    findings, inventory = resource_check(context_path)
    inventory[str(selection_path)] = file_sha(selection_path)
    if reviews_path:
        inventory[str(reviews_path)] = file_sha(reviews_path)
    review_matches = (reviews.get("selection_id") == selection["selection_id"] and
                      reviews.get("selection_payload_sha256") == stable_sha256(selection))
    if input_policy == "strict" and (reviews.get("selection_id") != selection["selection_id"] or
            reviews.get("selection_payload_sha256") != stable_sha256(selection)):
        findings.append(issue("review_selection_mismatch", reviews_path, "基于当前 selection 重新登记，勿沿用其他清单审核。"))
    if input_policy == "strict" and not targets.VERIFIED_PLATO_SCHEMA_FINGERPRINTS:
        findings.append(issue("production_schema_allowlist_empty", "没有已批准的生产 PLATO schema",
                              "取得真实单查询 fixture、结构与 release 证明并审核解析；不自动放开白名单。"))
    if not cache_path.is_dir():
        findings.append(issue("cache_directory_missing", cache_path, "恢复原始单药 PLATO 缓存。"))
    cache_inventory = {str(path.resolve()): file_sha(path) for path in sorted(cache_path.rglob("*")) if path.is_file()}
    inventory.update(cache_inventory)
    trusted_hashes = sorted(set(trusted_hashes))
    if any(not re.fullmatch(r"[0-9a-f]{64}", value) for value in trusted_hashes):
        raise ValueError("Trusted PLATO SHA-256 values must be 64 lowercase hex characters")
    config = {"plato": {"cache_dir": str(cache_path), "trusted_export_sha256": trusted_hashes,
                        "ignore_stereochemistry": ignore_stereochemistry}}
    if input_policy == "test":
        config.update(run_mode="test", generated_from_test_fixture=False,
                      test_source_review_exception=test_policy(cache_inventory),
                      network={"randomizations": 2000, "sensitivity_randomizations": 2000, "random_seed": 42},
                      credibility={"stability_randomizations": 2000})
    pipeline, context_ready = None, False
    try:
        pipeline_type = DevelopmentTestPipeline if input_policy == "test" else ExplanationPipeline
        pipeline = pipeline_type(config=config)
        if not any(item["code"].startswith(("resource_", "context_config", "runtime_manifest")) for item in findings):
            pipeline.load_production_resources(context_path)
            context_ready = True
    except (OSError, ValueError, TypeError, ImportError, RuntimeError) as exc:
        findings.append(issue("production_context_rejected", exc, "按生产加载器错误修复资源/环境，保留全部锁定校验。"))
    case_audits = audit_cases(cases, reviews if review_matches else {"drugs": {}}, pipeline,
                             context_ready, input_policy)
    for record in case_audits:
        record["ready"] = record["ready"] and not findings
    config = deepcopy(pipeline.config) if pipeline else config
    role_path = Path(config.get("pharmacology_roles", {}).get("path", ""))
    if role_path.is_file():
        inventory[str(role_path.resolve())] = file_sha(role_path)
    report = {"contract": CONTRACT, "checked_at": utc_now(), "selection_id": selection["selection_id"],
              "scope": "development_engineering_only", "formal_result_eligible": False,
              "input_policy": input_policy, "manual_acceptance_required": input_policy == "strict",
              "reference_records_match_selection": review_matches if reviews_path else None,
              "all_ready": bool(case_audits) and all(item["ready"] for item in case_audits),
              "n_cases": len(cases), "global_issues": findings, "cases": case_audits,
              "dependencies": _dependency_versions(), "inventory": inventory,
              "effective_config": config,
              "upstream_provenance": {key: selection.get(key) for key in (
                  "selection_id", "protocol_id", "suite_id", "checkpoint_kind", "training_seeds", "sources", "policy")},
              "upstream_checkpoint_reverified_this_run": False,
              "source_snapshot": _code_snapshot(PACKAGE)}
    return pipeline, cases, report


def verified_json(path):
    path = Path(path)
    sidecar = Path(str(path) + ".sha256")
    expected = sidecar.read_text(encoding="utf-8").split()[0]
    if file_sha(path) != expected:
        raise ValueError(f"Hash mismatch: {path}")
    return read_json(path)


def verify_receipt(root, receipt, suite_id):
    if receipt.get("suite_id") != suite_id or receipt.get("state") != "complete":
        raise ValueError("Completed receipt identity mismatch")
    expected_names = {"explanation.json", "explanation.json.report.md"}
    if set(receipt.get("artifacts", {})) != expected_names:
        raise ValueError("Completed receipt must cover JSON and readable report")
    for name, digest in receipt["artifacts"].items():
        path = root / name
        if file_sha(path) != digest or Path(str(path) + ".sha256").read_text(encoding="utf-8").split()[0] != digest:
            raise ValueError(f"Result artifact changed: {path}")
    payload = read_json(root / "explanation.json")
    audit = payload["audit"]
    case_audit = audit.get("engineering_case", {})
    if (case_audit.get("suite_id") != suite_id or
            any(case_audit.get(key) != value for key, value in receipt["task"].items())):
        raise ValueError("Explanation case identity mismatch")
    claimed = audit.pop("result_payload_sha256")
    if stable_sha256(payload) != claimed:
        raise ValueError("Internal explanation payload hash mismatch")


def summarize(root):
    """Verify the frozen plan and every successful receipt; never ignore failures."""
    root = Path(root)
    manifest = verified_json(root / "manifest.json")
    if stable_sha256(manifest["identity"]) != manifest["suite_id"]:
        raise ValueError("Suite identity hash mismatch")
    rows, problems = [], []
    for task in manifest["identity"]["tasks"]:
        directory = root / "tasks" / task["case_id"] / task["variant"]
        attempts = sorted(directory.glob("attempt_*"))
        row = {**task, "state": "pending"}
        if attempts:
            latest = attempts[-1]
            row["attempt"] = latest.name
            try:
                receipt = verified_json(latest / "status.json")
                if receipt.get("suite_id") != manifest["suite_id"] or receipt.get("task") != task:
                    raise ValueError("Task receipt does not match frozen plan")
                row.update({key: receipt[key] for key in ("state", "elapsed_seconds", "error", "issues", "metrics") if key in receipt})
                if row["state"] == "complete":
                    verify_receipt(latest, receipt, manifest["suite_id"])
            except (OSError, ValueError, TypeError, KeyError, IndexError) as exc:
                row["state"] = "corrupt"
                problems.append({**task, "error": str(exc)})
        rows.append(row)
    counts = dict(Counter(row["state"] for row in rows))
    comparisons = []
    baselines = {row["case_id"]: row for row in rows if row["variant"] == "full" and row["state"] == "complete"}
    for row in rows:
        if row["variant"] == "full":
            continue
        base = baselines.get(row["case_id"])
        comparison = {"case_id": row["case_id"], "variant": row["variant"], "status": "not_evaluable"}
        if base and row["state"] == "complete":
            left, right = base["metrics"], row["metrics"]
            comparison.update(status="descriptive_comparison", mechanism_unchanged=(
                left["mechanism_type"] == right["mechanism_type"]))
            comparison["top_target_jaccard"] = {}
            for side in ("drug1", "drug2"):
                first, second = set(left["top_targets"][side]), set(right["top_targets"][side])
                comparison["top_target_jaccard"][side] = len(first & second) / len(first | second) if first | second else None
            comparison["statistical_target_sets_unchanged"] = left["statistical_target_sets"] == right["statistical_target_sets"]
            comparison["network_delta"] = {
                key: (right["network"].get(key) - left["network"].get(key)
                      if all(isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)
                             for value in (left["network"].get(key), right["network"].get(key))) else None)
                for key in ("d_A_obs", "d_B_obs", "d_AB", "s_AB", "p_A", "p_B")}
            comparison["both_nulls_adequate"] = all(
                item["network"].get("quality_status") == "adequate" for item in (left, right))
        comparisons.append(comparison)
    report = {"contract": CONTRACT, "checked_at": utc_now(), "suite_id": manifest["suite_id"],
              "scope": "development_engineering_only", "n_cases": len(manifest["identity"]["cases"]),
              "n_planned": len(rows), "counts": counts,
              "complete": bool(rows) and counts.get("complete", 0) == len(rows) and not problems,
              "formal_result_eligible": False, "biological_validation_performed": False,
              "batch_fdr_applied": False,
              "statistical_scope": "descriptive_case_comparisons_no_cross_case_significance_claim",
              "problems": problems, "tasks": rows, "comparisons_to_full": comparisons}
    write_json_with_sha256(root / "summary.json", report)
    return report


def run_suite(root, pipeline, cases, preflight_report, variants=("full",), resume=False):
    """Never overwrite attempts. Success is skipped only after artifact verification."""
    root = Path(root).resolve()
    if not variants or len(set(variants)) != len(variants) or set(variants) - set(CONTROLS):
        raise ValueError("Variants must be unique supported controls")
    tasks = [{"case_id": case["case_id"], "variant": variant} for case in cases for variant in variants]
    identity = {"contract": CONTRACT, "cases": cases, "tasks": tasks,
                "config": preflight_report["effective_config"], "inventory": preflight_report["inventory"],
                "source_snapshot": preflight_report["source_snapshot"],
                "dependencies": preflight_report["dependencies"],
                "runtime": {"python": sys.version, "system": platform.system(), "machine": platform.machine()},
                "preflight": {key: preflight_report[key] for key in ("cases", "global_issues", "upstream_provenance")}}
    suite_id = stable_sha256(identity)
    def verify_live_inputs():
        for path, expected in identity["inventory"].items():
            if file_sha(path) != expected:
                raise ValueError(f"Input changed since preflight: {path}")
        if _code_snapshot(PACKAGE)["sha256"] != identity["source_snapshot"]["sha256"]:
            raise ValueError("Runtime source changed since preflight")
    verify_live_inputs()
    root.mkdir(parents=True, exist_ok=True)
    lock = root / ".runner.lock"
    with lock.open("x", encoding="utf-8") as handle:
        handle.write(str(os.getpid()))
    try:
        manifest_path = root / "manifest.json"
        if manifest_path.exists():
            old = verified_json(manifest_path)
            if not resume or old.get("suite_id") != suite_id:
                raise ValueError("Output exists or source/config/input changed; use a NEW output directory")
            summary = summarize(root)
            if summary["problems"]:
                raise ValueError("Resume refused: saved artifacts corrupt; retain them and use a new directory")
        else:
            if any(path != lock for path in root.iterdir()):
                raise ValueError("Output directory is not empty and has no frozen manifest")
            write_json_with_sha256(manifest_path, {"suite_id": suite_id, "created_at": utc_now(), "identity": identity})
        write_json_with_sha256(root / "preflight.json", preflight_report)
        case_map = {case["case_id"]: case for case in cases}
        audits = {item["case_id"]: item for item in preflight_report["cases"]}
        for task in tasks:
            directory = root / "tasks" / task["case_id"] / task["variant"]
            previous = sorted(directory.glob("attempt_*"))
            if previous and verified_json(previous[-1] / "status.json")["state"] == "complete":
                continue  # all completed receipts already verified by summarize()
            number = max([int(path.name.split("_")[-1]) for path in previous] or [0]) + 1
            attempt = directory / f"attempt_{number:04d}"
            attempt.mkdir(parents=True, exist_ok=False)
            started = time.monotonic()
            receipt = {"suite_id": suite_id, "task": task, "state": "running", "started_at": utc_now(),
                       "case_provenance": case_map[task["case_id"]]}
            write_json_with_sha256(attempt / "status.json", receipt)
            try:
                audit = audits[task["case_id"]]
                if not audit["ready"]:
                    receipt.update(state="blocked", issues=preflight_report["global_issues"] + audit["issues"])
                else:
                    verify_live_inputs()
                    pipeline.config = deepcopy(preflight_report["effective_config"])
                    variant = task["variant"]
                    pipeline.config["engineering_experiment"] = {"variant": variant if variant in VARIANTS else "full"}
                    if variant.startswith("seed_"):
                        pipeline.config["network"]["random_seed"] = int(variant[5:])
                    if variant.startswith("top_k_"):
                        pipeline.config["top_k_targets"] = int(variant[6:])
                    result = pipeline.explain_combination(**case_map[task["case_id"]]["input"])
                    result.metadata["audit"]["engineering_case"] = {"suite_id": suite_id, **task,
                        "selection": preflight_report["upstream_provenance"],
                        "candidate": case_map[task["case_id"]]["upstream"],
                        "identity_notes": {side: {key: value.get(key) for key in (
                            "observations", "identity_reference_record", "manual_acceptance_required")}
                            for side, value in audit.get("inputs", {}).items()}}
                    pipeline.save_result(result, attempt / "explanation.json")
                    network = result.mechanism.metrics.get("cheng", {})
                    receipt.update(state="complete", artifacts={name: file_sha(attempt / name) for name in (
                        "explanation.json", "explanation.json.report.md")}, metrics={
                        "mechanism_type": result.mechanism.mechanism_type,
                        "statistical_target_sets": result.mechanism.metrics.get("statistical_target_sets", {}),
                        "network": {key: network.get(key) for key in (
                            "d_A_obs", "d_B_obs", "d_AB", "s_AB", "p_A", "p_B", "quality_status", "n_valid_randomizations")},
                        "top_targets": {"drug1": [x.target_id for x in result.drug1_targets.top_targets],
                                        "drug2": [x.target_id for x in result.drug2_targets.top_targets]},
                        "target_counts": [len(result.drug1_targets.all_targets), len(result.drug2_targets.all_targets)]})
            except Exception as exc:
                receipt.update(state="failed", error=f"{type(exc).__name__}: {exc}")
            finally:
                receipt.update(finished_at=utc_now(), elapsed_seconds=round(time.monotonic() - started, 3))
                write_json_with_sha256(attempt / "status.json", receipt)
            summarize(root)
        return summarize(root)
    finally:
        lock.unlink()
