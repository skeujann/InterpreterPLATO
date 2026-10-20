# 路径连通性与输出合同（5.4.1）

本功能组织现有生物证据，适用于检查药物—靶点—PPI—疾病链路和双药靶点的可达性。
它不调用协同预测模型，不生成表达扰动实验，不改变 PLATO 候选、PPI 原始边或网络统计判据。

## 搜索和排序

方法标识为 `hop_bounded_weighted_paths_v1`。在 `(基因, 已走跳数)` 状态空间上搜索，
以现有 `edge_distance` 为代价，在指定最大跳数内选最小总代价路线。相同代价依次按
较少跳数、基因符号字典序破同分。不能先求无限制加权最短路再按长度丢弃，否则会漏掉
稍贵但合法的短路。

搜索输入是两个药物的完整合格靶点集，不再因显示 top-k 截断搜索输入。每个靶点至多
给出一条到最近疾病基因的受限加权路径；双药连接从两侧分别找最近的对侧靶点并合并去重。
这些是代表性最优路径，不是所有简单路径枚举。展示排序仍考虑靶点分数、边置信度、跳数
与枢纽惩罚；双药路径的靶点项取两端分数的几何均值，仅用于展示排序。

同一基因同时属于靶点和疾病集时，保留零跳记录并标注 `direct_disease_membership=true`；
共享药物靶点标注 `shared_target=true`。零跳路径的 `ppi_edge_count=0`、
`edge_confidence_applicable=false`；为兼容旧消费者，`edge_confidence_coverage` 仍是数值 0。
没有置信度的真实 PPI 边继续使用既有保守因子，并明确记录置信度缺失。

## API 和 CLI

```python
pipe = ExplanationPipeline(config={
    "path_search": {"max_length": 5, "max_paths": 10, "max_cross_drug_paths": 5}
})
```

对应 CLI 选项对 `summary`、`batch` 和 `plato-query-pair` 都有效；查询入口在请求 PLATO 前
检查参数，并把实际选项写入运行清单：

```text
--path-max-length 5 --max-evidence-paths 10 --max-cross-drug-paths 5
```

最大跳数允许 0–20，两个显示预算允许 0–100；错误类型、负数、未知键均报错。
预算为零表示不展示该类路径，诊断仍保留已找到的候选计数。`advanced=false` 时报告
`network_connectivity.status=disabled`，不冒充已经评估。

## 输出

- `local_evidence.top_paths` / `top_evidence_paths`：疾病路径，兼容既有字段。
- `local_evidence.cross_drug_paths`：双药靶点之间的 PPI 路径，单独报告。
- `local_evidence.connectivity` / `network_connectivity`：完整靶点计数、PPI 映射、
  可达性、长度内候选与展示计数，以及逐靶点原因。两处为同一报告。
- `evidence_graphs.explanation_evidence_graph`：疾病路径和双药路径共用去重后的无向 PPI 边；
  边的 `path_kinds` 保留用途。路径节点加入后再连接选定疾病基因集，避免端点遗漏。
  机制摘要与两种药物之间的 `summarizes_combination_member` 边只表示报告归属，
  标记为 `report_structure_only`，不属于 PPI 或药理证据。

疾病路径原因包括 `missing_ppi`、`target_not_in_ppi`、`missing_disease_context`、
`disease_not_in_ppi`、`disconnected_from_disease`、`path_length_limit`、
`no_finite_weighted_path`、`display_limit`、`shown`。`display_limit` 说明已找到路线，
不能写成网络不连通。双药诊断独立于疾病集，缺少疾病上下文时仍能显示双药现有连接。

## 解释边界

连通不等于药物协同、直接药效作用、因果链或模型贡献。通路归属、疾病归属和规则假说
属于不同边类型，不能当作额外 PPI 边计算连接。PPI 边只来自已加载网络；外部调用传入
不存在的路径节点或边时拒绝构造证据图。原有全靶点网络统计、随机背景、支持门槛保持独立。
