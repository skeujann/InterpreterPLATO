# 输出展示合同 1.0

5.4.2 新增 `explanation_layers`，中文摘要、终端和 Markdown 报告共用该展示层。
`result_schema_version` 仍为 5.4；机器使用的分类、分值、validation、统计和审计字段保持原义。

## 展示顺序

1. 候选靶点与药理角色：完整候选数、仅预测证据数、文献注释角色和已有作用方向。
2. 双药连接：实际展示的 PPI 路径和共享靶点零跳分开计数。搜索候选数不称作所有可能路径数。
3. 共享候选靶点：完整候选列表的标识符交集；另存主分类排序子集的交集，避免混用范围。
4. 靶点—疾病连接：非零跳路径、疾病集合成员零跳、映射及断连诊断。
5. 通路与功能：两侧分别显著富集的共同通路、原始共同注释和不同通路配对分别记录。
6. 规则分类及动态状态：保留分类来源和规则分语义。动态假说为空时明确没有时序输出。
7. 疾病邻近统计：疾病集合名称与大小、每侧 p、校正值、抽样接受率和拒绝原因。

“未达到疾病邻近显著性标准”只用于可以评估的检验；未运行、断连、随机背景门槛未达到
等状态报告具体原因，不当成阴性结果。批量结果读取 `batch_near_*` 与
`batch_test_*_evaluable`，不把未校正结果覆盖到批次结果上。不可评估假设的占位 q 值
不当成有效检验结果。原始数值继续保存在 JSON。

单对分类标签保留计算时的含义，可能在批次校正后不再达到显著性标准，所以与最终统计
分别展示。路径中使用无方向连接符，PPI 连接不写成抑制/激活链。

## 文件与接口

- `CombinationExplanation.print_summary()`：中文分层摘要。
- `CombinationExplanation.to_dict()`：新增 `explanation_layers`，中文摘要进入
  `integrated_explanation.summary`；其余原有结果字段不变。
- `save_result` / `save_results`：写原有 JSON 与 SHA-256，并新增
  `原文件名.report.md` 与 `.sha256`，返回 `report_path` 等字段。
- 在项目目录执行 `python scripts/render_explanation.py --input FILE.json --output FILE.md`：
  离线重排旧结果。接受单个解释、解释列表、或带 `error` 的失败记录。
  此命令不运行靶点查询、解释器推断或上游预测模型，也不改写输入。
  已安装为包时也可通过 `python -m InterpreterPLATO.scripts.render_explanation` 调用。

统计方法、90% 接受率门槛、模型权重、靶点筛选和评分均未由展示层修改。
