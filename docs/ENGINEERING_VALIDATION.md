# 毕业论文工程验证：当前执行方案

版本 `engineering_validation_v1`；更新日期 **2026-10-13**。

当前目标是验证模块**可用、可追溯、对输入变化作出合理响应**。
按用户要求，5.5.1 默认恢复原测试方式，不要求人工验收、来源签字或生产 schema 批准。
5.5.2 起默认忽略立体标记进行结构匹配；同位素、电荷和盐组分仍区分。
原始 PLATO 提交结构完整保留，详见 [结构匹配说明](STRUCTURE_MATCHING.md)。
不把本方案写成真实机制准确率验证。
旧版 12 开发 / 24 留出 / 双人盲评方案保留为可选后续研究，不是当前交付条件。

2026-10-13 切换测试方式之前的真实材料采集与阻断历史，见
[本地开发案例说明](../research/README.md)（原六例来源准备记录由数据保管者单独提供）。
其中的人工审批要求不再是默认测试入口的前置条件；结构差异与数据缺失记录继续保留。

## 1. 先冻结分母和解释边界

使用 `research/score_selected_dev6_20261009/selection.json` 中既有六例，保持候选编号、
顺序、ScoreRI 三种子预测、样本行号和预测/检查点来源哈希。不按解释漂亮程度换案例。
它们是按预测分数选出的开发案例，不是独立留出或随机抽样；重复药物不能视为独立证据。
本次解释不重新执行 ScoreRI 推理，也不重新校验上游全部检查点；保留原选择登记的证明链。

`complete=true` 只表示冻结计划中的软件任务全部完成且产物哈希通过。
所有解释仍然 `formal_result_eligible=false`。失败、受阻、靶点不足、未知药理作用也进入分母。
不能把软件测试通过数或成功运行比例写成机制真实性、泛化性能或协同因果解释。

## 2. 本地命令（从 InterpreterPLATO 目录执行）

Windows 使用项目解释器 `.venv\Scripts\python.exe` 替换下列 `python`。
Linux 使用已按项目锁文件配置的环境。所有命令均离线、顺序执行，不向 PLATO 发请求。

直接运行，无需 `prepare`、审核表或可信来源哈希参数：

```text
python scripts/run_engineering.py run --selection research/score_selected_dev6_20261009/selection.json --output-dir outputs/dev6_test_v1
```

默认 `--input-policy test` 只豁免原测试方式的六项来源批准门禁。
缓存原始字节哈希、单查询 POST 与准确结构绑定、已见响应格式、可靠人源蛋白筛选、
上下文资源完整性继续自动检查。哈希冻结表示文件未变，不等于独立来源批准。
结果记录 `run_mode=test`，保持原始来源状态；不会改写 production schema 白名单。
随机化、敏感性和稳定性预算均沿用原测试的 2000 次，随机种子为 42；90% 门槛不变。

`--reviews` 是可选参考资料：在测试模式只保留身份检查结果及立体化学差异，
不要求填写签字，也不因 `pending` 拦截；原输入 SMILES 不会自动替换为参考结构。
参考资料只有绑定同一 selection 时才附加到案例；否则标记匹配失败，不引用为该案例的身份记录。

仅查看自动输入诊断而不计算解释：

```text
python scripts/run_engineering.py preflight --selection research/score_selected_dev6_20261009/selection.json --output-dir outputs/dev6_preflight_v1
```

`preflight.json` 的每个实际数据问题都有 `code/detail/action`。
返回码 2 表示至少一例没有完成，并不要求人工验收。缺少缓存、文件损坏或零可靠靶点的案例
仍保留在六例分母中，其他能运行的案例继续执行。

如需旧的严格来源流程，可显式使用 `--input-policy strict --reviews <审核文件>`，
并提供已审核的 `--trusted-plato-sha256`。该模式仍检查生产 schema；默认测试无需此步骤。

可选的解释器对照计划（普通重跑不需要启用）：

```text
python scripts/run_engineering.py run --selection research/score_selected_dev6_20261009/selection.json --output-dir outputs/dev6_controls_v1 --variants full no_expression_ranking raw_plato_order drop_top_target seed_43 seed_44
```

这是一组六例 × 六条件的 **36 个解释任务**，不是 36 次模型训练。
不加载 ScoreRI 模型，不执行模型表达特征替换实验。

## 3. 每种控制到底改变什么

| 条件 | 实际变化 | 不得声称 |
| --- | --- | --- |
| `full` | 原完整流程 | 已证明真实机制 |
| `no_expression_ranking` | 表达排序权重设为零，其他可用维度重归一化；表达注释与质量检查保留 | 去掉全部表达资源的消融 |
| `raw_plato_order` | 按 PLATO 原始 score 降序（同分按 target_id）；上下文数值与资源保留 | 完全无上下文 / PLATO-only 机制基线 |
| `drop_top_target` | 每药删除一个完整上下文排序的首位靶点；网络、通路、规则机制、敏感性重新计算 | 真实靶点敲除或 ScoreRI 梯度/模型消融 |
| `seed_43` / `seed_44` | 只改变网络/敏感性随机种子，基准 full 为 42 | 生物重复实验 |
| `top_k_5` / `top_k_20`（可选） | 改显示和探索性规则所用 top-k；完整统计靶点集合不变 | 完整靶点集合扰动 |

删除控制不足两个合格靶点时明确失败，不能返回伪造的空机制。
原始 PLATO 缓存、PPI、表达文件不作修改。before/after/removed 集合写入结果审计。
“合理响应”不是要求所有扰动都改变机制类别：不受此条件影响的完整集合统计应保持一致；
真正移除靶点须重算，类别可以仍相同。随机种子变化不应更改输入靶点身份。

## 4. 输出和续跑

输出包含 `manifest.json`、`preflight.json`、`summary.json` 和各自 `.sha256`。
每个 `tasks/<case_id>/<variant>/attempt_0001/` 保存 `status.json`；成功时另存
`explanation.json`、中文报告及哈希。状态含 `blocked/failed/complete`，中断可能保留 `running`。
保存完整源选择信息、代码快照、配置、依赖版本、资源/缓存/审核登记哈希、耗时和异常。
软件版本为 5.5.2，解释结果仍兼容 schema 5.4，新增结构匹配策略及审计字段。

```text
python scripts/run_engineering.py summarize --output-dir outputs/dev6_controls_v1
```

相同命令加 `--resume` 可继续：只有计划身份一致、成功回执和 JSON/中文报告哈希均匹配才跳过。
失败/受阻的再次执行保存新 attempt，不删除旧记录。更换源码、资源、缓存、审核记录、
条件或配置须用**新目录**，避免混算。发现产物损坏时拒绝续跑，不覆盖旧证据。
`.runner.lock` 阻止同目录并行；硬终止留下锁时，先确认没有活跃进程，再人工移走该锁。

工程汇总仅描述案例和控制，不汇总跨案例显著性或独立证据支持率；
`batch_fdr_applied=false`。若以后做跨案例推断，须另外预定义完整检验家族并校正，
不能把此处逐案例 q 值当成全实验 q 值。
`comparisons_to_full` 记录 top 靶点集合 Jaccard、完整统计集合是否改变、机制类别是否改变、
网络数值差异和零模型充分性；任何一侧未完成均为 `not_evaluable`，不以零变化填补缺失。

可选的软件与研究检查，不是默认运行前置条件：

1. 离线回归：输入拒绝、顺序执行、完整分母、失败保留、续跑/篡改、真实重算。
2. 六例真实输入诊断：保留结构差异，检查缓存、上下文及文件完整性。
3. 六例完整输出：每例从 ScoreRI 编号追溯到 PLATO、上下文和解释产物。
4. 对照/敏感性：保留不变、变化和无法评估的结果，不只摘漂亮案例。
5. 论文展示：流程图、六例状态表、典型案例证据链、变化表、失败与限制。

药理作用缺乏可追溯资料时保持 `unknown`。扩充角色库是后续人工文献整理，
不得把靶点预测、网络连接或模型高分改写成已知抑制/激活方向。
