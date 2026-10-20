"""Read-only, finding-first presentation of live or archived explanations.

This module does not query targets, search a graph, or recompute statistics.
Legacy inference fields remain authoritative; presentation is an additive view.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Mapping, Sequence
import hashlib
import math


PRESENTATION_VERSION = "1.0"
CLASSIFICATION_LABELS = {
    "network_close_distinct_processes": "两组候选靶点网络接近，具有不同的通路注释",
    "distinct_annotated_processes": "两组候选靶点具有不同的富集通路注释",
    "network_proximal_targets": "两组候选靶点在 PPI 网络中接近",
    "shared_pathway_co_modulation": "两侧候选靶点共同富集于部分通路",
    "shared_target_co_modulation": "排序子集共享候选靶点，作用方向尚未确定",
    "target_overlap_or_redundancy": "排序子集共享靶点，已记录的作用方向一致",
    "opposing_shared_target_modulation": "共享靶点中存在相反的已记录作用方向",
    "cell_cycle_phase_co_modulation": "候选靶点分别涉及复制期与有丝分裂期过程",
    "cell_cycle_phase_complementarity": "复制期与有丝分裂期过程具有抑制方向注释",
    "damage_repair_co_modulation": "候选靶点分别涉及 DNA 损伤与修复过程",
    "replication_damage_reinforcement": "DNA 损伤与修复过程具有抑制方向注释",
    "vertical_pathway_coblockade": "共同富集通路中的候选靶点具有抑制方向注释",
    "high_overlap": "单对规则分类：两侧疾病邻近，并存在靶点或网络邻域重叠",
    "complementary_exposure": "单对规则分类：两侧疾病邻近，网络邻域分离",
    "single_exposure": "单对规则分类：一侧疾病邻近，网络邻域分离",
    "overlapping_single_exposure": "单对规则分类：一侧疾病邻近，并存在靶点或网络邻域重叠",
    "unknown": "当前规则未选择单一主分类",
}
ROLE_LABELS = {
    "direct_pharmacodynamic": "直接药效作用文献注释",
    "family_pharmacodynamic": "蛋白家族层面的药效注释（family evidence；未确认具体亚型）",
    "drug_handling": "代谢、转运或细胞内保留注释",
    "context_association": "上下文关联注释",
    "unresolved": "尚未匹配药理角色注释",
}
CONNECTIVITY_LABELS = {
    "disabled": "路径搜索未启用",
    "not_evaluated": "未执行路径搜索",
    "missing_ppi": "未加载 PPI 网络",
    "missing_disease_context": "未加载疾病基因集合",
    "disease_not_in_ppi": "疾病基因未映射到 PPI 网络",
    "target_mapping_incomplete": "靶点映射后无法完成双侧搜索",
    "disconnected_target_sets": "两组已映射靶点位于不连通的网络部分",
    "path_length_limit": "连接超出当前路径长度上限",
    "target_not_in_ppi": "靶点未映射到 PPI",
    "disconnected_from_disease": "靶点所在网络部分与疾病集合不连通",
    "display_limit": "路径已找到，但未进入当前展示名额",
}


def _number(value: Any) -> str:
    if value is None:
        return "未计算"
    if isinstance(value, (int, float)) and math.isfinite(value):
        return f"{value:.6g}"
    return str(value)


def _safe(value: Any) -> str:
    """Keep provider names and error text from injecting Markdown structure."""
    text = " ".join(str(value).split())
    for char in ("\\", "`", "*", "_", "[", "]", "<", ">", "|", "#"):
        text = text.replace(char, "\\" + char)
    return text


def _path_records(paths: Sequence[Mapping[str, Any]]) -> list:
    return [{key: path.get(key) for key in (
        "drug", "destination_drug", "path", "length", "ppi_edge_count",
        "shared_target", "direct_disease_membership", "source_evidence_status",
        "destination_evidence_status",
    )} for path in paths]


def _path_counts(records: Sequence[Mapping[str, Any]]) -> Dict[str, int]:
    zero = sum(len(record.get("path") or []) == 1 for record in records)
    nonzero = sum(len(record.get("path") or []) > 1 for record in records)
    return {"displayed": len(records), "nonzero": nonzero, "zero_hop": zero}


def build_explanation_layers(payload: Mapping[str, Any]) -> Dict[str, Any]:
    """Build a JSON-serializable view without mutating the source payload."""
    names = [str(payload.get(f"drug{i}_name", payload.get(f"drug{i}", f"药物{i}"))) for i in (1, 2)]
    view: Dict[str, Any] = {
        "presentation_version": PRESENTATION_VERSION,
        "drug_names": names,
        "cell_line": payload.get("cell_line", ""),
        "status": "completed",
        "findings": [],
        "statistics": {},
        "run_mode": payload.get("run_mode", "未记录"),
        "source_result_sha256": payload.get("audit", {}).get("result_payload_sha256"),
        "structure_matching": {side: entry.get("structure_matching", {})
                               for side, entry in payload.get("audit", {}).get("plato_inputs", {}).items()},
    }
    if "mechanism_type" not in payload:
        if not payload.get("error"):
            raise ValueError("Expected an explanation JSON or an explicit failure record")
        reason = str(payload["error"])
        if "PLATO returned no eligible explicit-human protein targets" in reason:
            missing = [name for i, name in enumerate(names, 1)
                       if payload.get("eligible_targets", {}).get(f"drug{i}") == 0]
            if missing:
                reason = "、".join(missing) + " 的 PLATO 结果中没有通过当前筛选的人类蛋白候选靶点"
        view.update(status="not_generated", reason=reason, original_error=str(payload["error"]),
                    eligible_targets=dict(payload.get("eligible_targets", {})))
        return view

    def finding(key: str, title: str, text: str, **details: Any) -> None:
        view["findings"].append({"id": key, "title": title, "text": text, **details})

    profiles = [payload.get(f"drug{i}_targets", {}) for i in (1, 2)]
    pharmacology = payload.get("pharmacology", {})
    target_rows = []
    for i, (name, profile) in enumerate(zip(names, profiles), 1):
        targets = profile.get("targets", [])
        roles = pharmacology.get(f"drug{i}", profile.get("pharmacology_roles", {}))
        groups = roles.get("groups", {})
        target_rows.append({
            "drug": name, "count": len(targets),
            "predicted_only_count": sum(t.get("evidence_status") == "predicted_only" for t in targets),
            "roles": {key: [str(t["target_id"]) for t in groups.get(key, [])] for key in ROLE_LABELS},
            "literature_actions": [
                {"target": t["target_id"], **{key: record.get(key) for key in (
                    "action", "process", "evidence_scope", "experimental_context", "reference_ids",
                )}}
                for t in targets for record in t.get("pharmacology", {}).get("records", [])
            ],
        })
    finding("targets", "候选靶点与药理角色",
            "；".join(f"{r['drug']}：{r['count']} 个候选靶点" for r in target_rows), drugs=target_rows)

    local = payload.get("local_evidence", {})
    connectivity = payload.get("network_connectivity", local.get("connectivity", {}))
    cross = connectivity.get("cross_drug", {})
    cross_paths = _path_records(local.get("cross_drug_paths", []))
    counts = _path_counts(cross_paths)
    if counts["displayed"]:
        text = (f"展示 {counts['nonzero']} 条非零跳双药 PPI 连接，"
                f"以及 {counts['zero_hop']} 条共享候选靶点零跳记录")
    elif cross.get("candidate_path_count", 0):
        text = f"搜索得到 {cross['candidate_path_count']} 条候选记录，当前展示名额为 0"
    else:
        status = cross.get("status", connectivity.get("status", "not_evaluated"))
        text = CONNECTIVITY_LABELS.get(status, "当前输出没有双药路径记录")
    finding("cross_drug", "双药连接", text, counts=counts, paths=cross_paths,
            candidate_count=cross.get("candidate_path_count"),
            search_scope=cross.get("search_scope"), max_length=cross.get("max_length"))

    target_sets = [{str(t["target_id"]) for t in profile.get("targets", [])} for profile in profiles]
    shared = sorted(target_sets[0] & target_sets[1])
    finding("shared_targets", "共享候选靶点",
            f"完整候选列表共有 {len(shared)} 个共享靶点" if shared else "完整候选列表未发现共享靶点",
            targets=shared, scope="all_eligible_target_identifiers",
            classification_subset=payload.get("mechanism", {}).get("shared_targets", []))

    disease_paths = _path_records(local.get("top_paths", payload.get("top_evidence_paths", [])))
    disease_counts = _path_counts(disease_paths)
    disease_text = (f"展示 {disease_counts['nonzero']} 条非零跳靶点—疾病 PPI 路径，"
                    f"以及 {disease_counts['zero_hop']} 条靶点本身属于疾病集合的零跳记录")
    if not disease_paths:
        status = connectivity.get("status", "not_evaluated")
        disease_text = CONNECTIVITY_LABELS.get(status, "当前输出没有靶点—疾病路径记录")
        if connectivity.get("candidate_path_count", 0):
            disease_text = f"搜索得到 {connectivity['candidate_path_count']} 条靶点—疾病候选记录，当前未展示"
    diagnostics = []
    for drug in connectivity.get("drugs", []):
        diagnostics.append({
            "drug": drug.get("drug"),
            "eligible": drug.get("eligible_target_count"), "mapped": drug.get("mapped_target_count"),
            "reachable": drug.get("reachable_target_count"),
            "within_length": drug.get("within_length_target_count"),
            "exceptions": [{"target": t.get("target"), "status": t.get("status")}
                           for t in drug.get("targets", []) if t.get("status") not in {"shown", "display_limit"}],
        })
    finding("disease_paths", "靶点与疾病连接", disease_text,
            paths=disease_paths, counts=disease_counts, diagnostics=diagnostics)

    mechanism = payload.get("mechanism", {})
    metrics = mechanism.get("metrics", {})
    common = mechanism.get("shared_pathways", [])
    descriptive = next((e.get("pathways", []) for e in mechanism.get("evidence", [])
                        if e.get("type") == "descriptive_shared_pathways_not_enriched"), [])
    finding("pathways", "通路与功能",
            (f"发现 {len(common)} 项两侧分别通过 BH 校正的共同富集通路注释" if common else
             f"未检出两侧共同显著富集的通路；另有 {len(descriptive)} 项共同通路注释") ,
            significant_shared=list(common), descriptive_shared=list(descriptive),
            functional_type=payload.get("functional_mechanism_type", "unknown"),
            distinct_annotation_pairs=list(mechanism.get("complementary_pathways", [])),
            enrichment_details=metrics.get("pathway_enrichment", {}).get("significant_shared", []))
    primary = str(payload["mechanism_type"])
    finding("classification", "规则分类", CLASSIFICATION_LABELS.get(primary, primary),
            code=primary, score=payload.get("mechanism_score", payload.get("mechanism_confidence")),
            score_semantics="uncalibrated_rule_strength_not_probability")
    dynamics = payload.get("dynamic_state_evidence", [])
    finding("dynamic", "动态状态",
            "未输出动态或时序状态假说" if not dynamics else
            f"输出 {len(dynamics)} 项动态状态假说：" + "、".join(str(d.get("state_type", "未命名")) for d in dynamics),
            hypotheses=dynamics)

    network = payload.get("network_statistics", metrics.get("cheng", {}))
    has_batch = "batch_multiple_testing_method" in network or "batch_q_A" in network
    sides = []
    for i, side in enumerate(("A", "B"), 1):
        evaluable = (network.get(f"batch_test_{side}_evaluable", False) if has_batch else
                     bool(network.get("valid") and network.get("quality_status") == "adequate"))
        near = network.get(f"batch_near_{side}" if has_batch else f"near_{side}") if evaluable else None
        audit = network.get("randomization_audit", {}).get(f"drug{i}", {})
        requested, accepted = audit.get("requested"), audit.get("accepted")
        sides.append({
            "drug": names[i - 1], "evaluable": bool(evaluable), "near": near,
            "status": ("显著邻近当前疾病模块" if near else "未达到当前疾病邻近显著性标准")
                      if evaluable else "本次未作疾病邻近显著性判定",
            "p": network.get(f"p_{side}"), "z": network.get(f"z_{side}"),
            "adjusted_p": network.get(f"p_{side}_adjusted"),
            "batch_q": network.get(f"batch_q_{side}") if has_batch else None,
            "batch_q_is_placeholder": bool(has_batch and not evaluable),
            "requested": requested, "accepted": accepted,
            "acceptance_rate": accepted / requested if requested and accepted is not None else None,
            "rejection_reasons": dict(audit.get("rejection_reasons", {})),
            "evaluation_reason": network.get(f"batch_test_{side}_status", ""),
        })
    view["statistics"] = {
        "scope": "靶点集合与疾病模块的邻近性检验",
        "disease_labels": list(payload.get("disease_context_labels", [])),
        "disease_source": payload.get("disease_context_source", ""),
        "disease_gene_count": payload.get("disease_context_size"),
        "multiple_testing_method": network.get("batch_multiple_testing_method" if has_batch else "multiple_testing_method"),
        "batch_applied": has_batch, "sides": sides,
        "quality_status": network.get("quality_status", "not_evaluated"),
        "warnings": list(network.get("warnings", [])),
        "valid": network.get("valid"),
    }
    return view


def render_explanation(view: Mapping[str, Any]) -> str:
    """Render the additive view as Chinese Markdown, also readable in a terminal."""
    names = view["drug_names"]
    lines = [f"# {_safe(names[0])} + {_safe(names[1])}", "",
             f"细胞系：{_safe(view.get('cell_line') or '未指定')}；运行模式：{_safe(view.get('run_mode'))}。", ""]
    if any(item.get("policy") == "ignore_stereochemistry" for item in view.get("structure_matching", {}).values()):
        lines += ["结构匹配：忽略立体标记，按分子连接关系识别；保留同位素、电荷及盐组分。", ""]
        for index, side in enumerate(("drug1", "drug2")):
            match = view.get("structure_matching", {}).get(side, {})
            if match.get("cross_stereochemistry_reuse"):
                lines += [f"{_safe(names[index])}：复用同一连接结构的已归档 PLATO 结果；原提交构型记录见审计字段。", ""]
    if view["status"] != "completed":
        lines += ["本次未生成双药解释。", "", f"原因：{_safe(view['reason'])}。"]
        for i, name in enumerate(names, 1):
            count = view.get("eligible_targets", {}).get(f"drug{i}")
            if count is not None:
                lines += ["", f"{_safe(name)}：通过筛选的候选靶点 {count} 个。"]
        return "\n".join(lines) + "\n"

    for finding in view["findings"]:
        key = finding["id"]
        lines += [f"## {finding['title']}", "", _safe(finding["text"]) + "。", ""]
        if key == "targets":
            for drug in finding["drugs"]:
                roles = [f"{ROLE_LABELS[k]}：{'、'.join(ids[:12])}" + (f"等，共 {len(ids)} 个" if len(ids) > 12 else "")
                         for k, ids in drug["roles"].items() if ids]
                lines.append(f"- {_safe(drug['drug'])}：仅预测证据的候选靶点 {drug['predicted_only_count']} 个。")
                lines.extend(f"  - {_safe(role)}。" for role in roles)
                for action in drug["literature_actions"]:
                    direction = {"inhibition": "抑制", "activation": "激活",
                                 "microtubule_stabilization": "微管稳定化",
                                 "polyglutamylation_substrate": "多聚谷氨酸化底物",
                                 "uptake_substrate": "摄取转运底物",
                                 "efflux_substrate": "外排转运底物",
                                 "metabolism_substrate": "代谢底物"}.get(action["action"], action["action"])
                    scope = "家族层面" if action["evidence_scope"] == "target_family" else "蛋白层面"
                    lines.append(f"  - 文献注释 {_safe(action['target'])}：{_safe(direction)}；{_safe(action['process'])}（{scope}）。")
        elif key in {"cross_drug", "disease_paths"}:
            for path in finding["paths"]:
                genes = path.get("path") or []
                if not genes:
                    continue
                kind = ("共享候选靶点零跳" if key == "cross_drug" else "疾病集合成员零跳") if len(genes) == 1 else f"{len(genes)-1} 条 PPI 边"
                pair = str(path.get("drug") or "")
                if key == "cross_drug":
                    pair += " / " + str(path.get("destination_drug") or "")
                lines.append(f"- {_safe(pair)}：{' — '.join(_safe(g) for g in genes)}（{kind}）。")
            if key == "cross_drug" and finding.get("candidate_count") is not None:
                lines += ["", f"本次搜索得到 {finding['candidate_count']} 条候选记录；以上为选出的展示记录。"]
            for drug in finding.get("diagnostics", []):
                lines.append(f"- {_safe(drug['drug'])}：{drug['mapped']}/{drug['eligible']} 个候选靶点映射到 PPI，{drug['reachable']} 个可达疾病集合，{drug['within_length']} 个在长度上限内。")
                for item in drug["exceptions"]:
                    lines.append(f"  - {_safe(item['target'])}：{_safe(CONNECTIVITY_LABELS.get(item['status'], item['status']))}。")
        elif key == "shared_targets":
            if finding["targets"]:
                lines.append("、".join(_safe(t) for t in finding["targets"]) + "。")
            lines += ["", "此处统计完整候选列表；主分类中的共享靶点字段使用排序子集，范围不同。"]
        elif key == "pathways":
            lines.extend(f"- {_safe(path)}。" for path in finding["significant_shared"])
            if finding["distinct_annotation_pairs"]:
                lines += ["", "不同富集通路的注释配对：", ""]
                lines.extend(f"- {_safe(pair)}。" for pair in finding["distinct_annotation_pairs"])
            if not finding["significant_shared"] and finding["descriptive_shared"]:
                examples = finding["descriptive_shared"][:5]
                lines += ["", "共同注释示例（最多列出 5 项）：" + "、".join(_safe(p) for p in examples) + "。"]
            lines += ["", "以上记录通路归属与富集情况；激活、抑制及通路间的作用方向按独立证据记录。"]
        elif key == "classification":
            lines += [f"规则标签：{_safe(finding['code'])}；规则分值：{_number(finding['score'])}（非概率）。",
                      "", "该标签保留单对分析的分类，疾病邻近结论及实际校正方法见下方统计结果。"]
        lines.append("")

    stats = view["statistics"]
    labels = "、".join(str(x) for x in stats["disease_labels"]) or str(stats["disease_source"] or "未记录名称")
    lines += ["## 疾病邻近统计", "", f"疾病集合：{_safe(labels)}；基因数：{_number(stats['disease_gene_count'])}。",
              "", "检验对象是各侧靶点集合与疾病模块的距离。路径记录单独展示。", "",
              f"多重检验方法：{_safe(stats['multiple_testing_method'] or '未记录')}。", ""]
    for side in stats["sides"]:
        values = f"p={_number(side['p'])}"
        if stats["batch_applied"]:
            values += ("；批次 q 为占位值，不用于显著性判定" if side["batch_q_is_placeholder"] else f"；批次 q={_number(side['batch_q'])}")
        else:
            values += f"；单对校正 p={_number(side['adjusted_p'])}"
        lines.append(f"- {_safe(side['drug'])}：{side['status']}（{values}）。")
        if side["acceptance_rate"] is not None:
            lines.append(f"  - 随机抽样接受 {side['accepted']}/{side['requested']} 次，接受率 {side['acceptance_rate']:.2%}。")
        rejected = side["rejection_reasons"]
        for code, count in rejected.items():
            if count:
                reason = {"incomplete_reachability": "随机集合中存在无法到达疾病集合的节点",
                          "incomplete_size": "随机集合未达到要求大小"}.get(code, code)
                lines.append(f"  - {_safe(reason)}：拒绝 {count} 次。")
    reasons = {
        "network_components_disconnected": "观测靶点集合中存在与疾病集合不连通的节点，本次未计算整组距离检验",
        "insufficient_valid_randomizations": "有效随机抽样数量或接受率未达到当前统计执行门槛",
        "randomizations_below_adequate_minimum": "请求的随机抽样次数未达到当前统计执行门槛",
        "degenerate_degree_matched_null": "随机背景距离分布没有足够的变异",
        "poor_degree_match_quality": "随机集合的节点度数匹配未达到当前门槛",
        "low_network_mapping_coverage": "靶点或疾病基因的网络映射覆盖未达到当前门槛",
    }
    notes = [reasons[code] for code in stats["warnings"] if code in reasons]
    if notes:
        lines += ["", "本次统计执行情况：", ""] + [f"- {note}。" for note in notes]
    lines += ["", "候选靶点、药理文献注释、PPI 连接和疾病邻近检验按各自证据层展示。", ""]
    return "\n".join(lines)


def render_payload(payload: Mapping[str, Any]) -> str:
    report = render_explanation(build_explanation_layers(payload))
    experiment = payload.get("audit", {}).get("engineering_experiment", {})
    variant = experiment.get("variant", "full")
    if variant != "full":
        note = {
            "raw_plato_order": "仅按 PLATO 原始分数重排；表内上下文分数保留，不代表按该分数排序。不是无上下文基线。",
            "no_expression_ranking": "仅移除表达排序权重；表达注释和质量检查保留。不是全资源消融。",
            "drop_top_target": "每药删除一个上下文首位靶点后重算下游；这是软件输入扰动，不是药理敲除实验。",
        }.get(variant, "工程控制，不代表独立机制验证。")
        report = f"工程对照：`{variant}`。{note}\n\n" + report
    return report


def write_readable_report(path: str | Path, payload: Any) -> Dict[str, str]:
    """Save a companion report without overwriting a JSON destination."""
    destination = Path(path)
    reports = payload if isinstance(payload, list) else [payload]
    content = "\n---\n\n".join(render_payload(record) for record in reports)
    destination.parent.mkdir(parents=True, exist_ok=True)
    raw = content.encode("utf-8")
    destination.write_bytes(raw)
    digest = hashlib.sha256(raw).hexdigest()
    sidecar = Path(str(destination) + ".sha256")
    sidecar.write_text(f"{digest}  {destination.name}\n", encoding="utf-8")
    return {"report_path": str(destination), "report_sha256": digest,
            "report_sha256_sidecar": str(sidecar)}


def companion_report_path(path: str | Path) -> Path:
    source = Path(path)
    return source.with_name(source.name + ".report.md")
