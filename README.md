# InterpreterPLATO 5.5.2

本项目是独立仓库，与 [ScoreRI](https://github.com/skeujann/ScoreRI) 分开维护。所有命令均从本仓库根目录执行；输入药物对和预测分数由调用方提供，运行时不需要相邻的 ScoreRI 源码。

## 安装与测试

```bash
python -m pip install -r requirements-dev.lock
python -m pip install --no-deps .
interpreter-plato --help
python -m pytest tests
```

仓库保留源码、测试、轻量配置和数据哈希清单。原始数据、处理后上下文、PLATO 缓存、六例开发选择及来源记录、输出和虚拟环境需在本地单独准备；详见 [数据准备](docs/DATA_BOOTSTRAP.md) 和 [本地开发案例说明](research/README.md)。

InterpreterPLATO 是一个药物组合生物学解释器。它的作用不是重新训练一个靶点预测模型，而是把一对药物组合的靶点证据、癌症细胞系上下文、PPI 网络、通路和疾病基因集合整合起来，生成可审计的机制解释。

## 默认忽略立体信息（5.5.2）

解释器默认使用 `plato.ignore_stereochemistry=true`：去除原子手性和双键立体标记，
按相同分子连接关系匹配 PLATO 缓存。同位素、电荷、键级及盐组分保持区分。
缓存仍保存真实查询的完整结构；同组存在多份预测时，只有靶点解析结果一致才复用，
结果冲突会明确列出，不把不同靶点集合拼成一份预测。药理文献注释仍保留原文所对应的准确结构。
工程入口可用 `--respect-stereochemistry` 恢复精确构型匹配；其余入口可在配置中将上述开关设为 false。
说明与验证边界见 [结构匹配说明](docs/STRUCTURE_MATCHING.md)。

## 六例开发案例运行（5.5.1）

`scripts/run_engineering.py` 默认恢复原测试方式：**不要求人工验收、来源签字或生产 schema 批准**。
无需 `prepare` 或 `--reviews` 即可执行；原始文件哈希、查询结构绑定及可靠靶点筛选继续生效。
保留逐例输入诊断、冻结计划、顺序运行、
失败/拒绝保留、哈希验证续跑及完整分母汇总。既有六个 ScoreRI 候选保持开发案例身份。
显式控制含表达排序权重关闭、PLATO 分数排序、真实靶点删除、随机种子与 top-k；
不同控制的严格范围及可执行命令见 [工程验证方案](docs/ENGINEERING_VALIDATION.md)。
默认完整解释算法不变，测试方式沿用 2000 次随机化；结果标记 `run_mode=test`，不伪造来源批准。
只有显式选择 `--input-policy strict` 时才启用人工来源门禁。主解释器的严格模式保持原行为。

## 发现优先的中文输出（5.4.2）

终端与 `integrated_explanation.summary` 先列出候选靶点、双药 PPI 路径、共享靶点、
靶点—疾病路径及通路注释，再单独呈现各药的疾病邻近检验。共享靶点零跳与非零跳 PPI
路径分别计数；显示数量与搜索候选数量分别记录。未生成解释时直接给出具体原因。
`explanation_layers` 提供对应的结构化展示字段，原有推断、统计、审计字段继续保留。

`save_result` / `save_results` 在 JSON 旁自动写出 `原文件名.report.md` 及其 `.sha256`。
单对和批量 CLI 使用相同的中文输出，批次结果显示校正后的 q 值；不可评估假设的占位 q
不显示为有效显著性结果。此版本不修改分类规则、评分、随机抽样或 90% 接受率门槛。

已有 JSON 可直接生成新报告，无需重新请求 PLATO 或运行模型：

```powershell
python scripts/render_explanation.py --input explanation.json --output explanation_report.md
```

详见 [输出字段与语义](docs/OUTPUT_PRESENTATION.md)。

## 论文研究定位（当前唯一输出语义）

本项目当前定位为“论文中的可审计机制假说模块”。它生成的是基于给定药物对、
PLATO 靶点和癌症上下文的 **post-hoc 生物学机制假说**，不是模型因果解释、实验验证、
临床疗效概率或真实作用机制证明。`predicted_synergy` 默认只是随结果记录的上游输入；
默认未接入上游模型重算，不得写成“解释了模型为何给出该分数”；调用方提供的扰动表即使
通过内部一致性检查，仍属于未独立核验的证据，不能替代模型行为与归因验证。

每个新结果都必须声明：

- `research_contract_version=1.0`；
- `run_mode=research_hypothesis | demo | test`；
- `formal_result_eligible=false` 及不可升级的原因；
- `generated_from_test_fixture`，防止合成测试进入论文结果；
- 代码快照、Git 状态、有效配置、依赖、数据清单、PLATO 输入和策略版本的审计记录；
- 结果内部 canonical JSON 哈希，以及落盘文件的精确字节 `.sha256` sidecar。

即使未来 PLATO schema 和来源完整性门禁通过，当前研究合同也不会自动把输出升级为
“正式机制结论”。门禁回答的是“输入证据是否允许进入分析”，研究合同回答的是
“分析结果最多能声称什么”。详细边界和论文报告要求见
[论文机制假说模块范围](docs/THESIS_INTERPRETATION_SCOPE.md)。

当前统计合同为 `statistical_contract_version=1.1`。网络分析采用明确版本化的自定义
`weighted_fixed_disease_closest_v3_mapped_universe`，不是 Cheng 等人的严格方法复现，也未完成外部校准。
完整合格 PLATO 靶点集分别与冻结网络/通路背景取交集用于统计；上下文 top-k 仅用于展示和探索性规则选择。
网络输出列出未映射基因、实际分析集合和原始覆盖率；覆盖率不足或已映射成员不可达时仍不建立支持。
机制分数不再因共享通路条数或描述性通路配对数量增加而加分，分数仍是未校准的规则值。
精确定义、零模型、支持判据和批量检验边界见[统计方法与推断边界](docs/STATISTICAL_METHODS.md)。
后续独立评测遵循[拟议验证协议](docs/VALIDATION_PROTOCOL.md)，其状态是
`proposed_not_executed`，不能将软件回归测试当作真实机制验证。

当前版本采用严格的 `PLATO-only` 靶点策略：

- 两个药物的候选靶点都必须来自 PLATO 单分子 target-fishing 查询结果。
- 每个药物都需要单独在 PLATO 平台查询、单独下载结果、单独本地归档。
- 本地 PPI、DepMap 表达、Reactome 通路、COSMIC 疾病基因集合只用于癌症上下文重排和机制解释。
- 本地上下文数据不能补充、替代或伪造药物靶点。
- 程序不做批量访问、不爬取 PLATO 网页；如需程序查询，只允许用户触发的单次 PLATO 表单请求，并保存完整审计证据。

也就是说，这里的 `PLATO-only` 指“正式药物靶点来源只认 PLATO”，不是指本程序拥有独立靶点预测模型、公开 API 授权或批量查询能力；PLATO 仍然是外部平台。

## 药物特异的文献角色（5.4）

文献角色层只注释已进入分析的 PLATO 候选，不增加候选或改变上下文排序、网络/富集统计。
摘要优先呈现有文献依据的药效靶点，并将代谢、转运、胞内保留及角色未知的候选分开报告。
角色按完整 InChIKey 和基因匹配，不能根据基因名称前缀、蛋白类型或 PPI 连接推断。
首批资源含两种药物、七条关系；TUBB3 只保留微管家族级证据，不宣称已验证该亚型结合。
未匹配关系保持未知。文献依据、适用范围、配置和扩展格式见
[文献角色方法](docs/PHARMACOLOGY_ROLES.md)。

## 路径连通与模块衔接（5.4.1）

在限定跳数内搜索最小加权代价路径，保留满足长度限制的替代路线；从完整合格靶点集搜索，
展示条数单独控制。补充双药靶点之间的 PPI 路径，并在路径节点齐全后连接选定疾病模块。
`local_evidence.connectivity` / `network_connectivity` 报告未映射、真正断连、长度受限及仅因展示
预算而未显示的靶点。直接命中疾病基因的零跳路径单独标记，不能解释成观测到 PPI 边。

API 的 `path_search` 配置与命令行 `--path-max-length`、`--max-evidence-paths`、
`--max-cross-drug-paths` 同步，默认分别为 5、10、5。JSON schema 保持兼容的 5.4，新增字段
采用 `hop_bounded_weighted_paths_v1` 方法标识。连通报告不增加靶点、不补造网络边，
也不为机制统计加分。具体字段、限制与示例见[连通性说明](docs/CONNECTIVITY.md)。

## 当前生产状态

真实 PLATO 输入进入研究机制假说分析，必须先同时满足以下输入门禁；通过这些门禁仍不
产生已验证的“正式机制结论”：

- PLATO 缓存 schema 为 `3`。
- 两个药物的原始 PLATO 响应字节流都被不可变保存，并记录 SHA-256。
- PLATO 响应中的查询结构必须能和提交给 PLATO 的 SMILES 对上。
- 结构绑定依赖 RDKit 解析、规范 SMILES 和 InChIKey。
- release 必须由 payload 明示，或由未来经独立审计、显式批准并进入 production
  allowlist 的 r37 fixture 绑定；当前没有这样的 fixture。
- 响应 schema 必须通过上述 production-approved fixture；当前 gate 关闭。
  受限原件与可发布的脱敏派生副本必须分开管理。
- 两个药物都必须保留至少一个可靠、明确人源、蛋白类靶点。
- PPI、表达、通路、疾病等上下文文件必须通过配置中的 SHA-256 校验。
- 两个 PLATO 原始导出文件的 SHA-256 都必须来自独立可信 allowlist。

当前观察到的 schema fingerprint 来自 Gemcitabine PLATO r37 target-fishing fixture
candidate，可用于回归检查 `chemblId / nome / score / validato / simili` parser contract。
它的 manifest 仍标记为 `raw_json_and_pdf_fixture_candidate` 和
`not_yet_allowlisted_for_scoring`；production schema allowlist 当前为空，因此该 candidate
会得到 `schema_fixture_verified=false`，不能绑定 PLATO/r37 来源或打开正式计分。仓库当前把
fixture manifest 和 target-fishing JSON 作为唯一的小型受控数据例外显式跟踪；PDF、普通
查询响应和其他受限原件不随代码提交。fixture 位于：

```text
InterpreterPLATO/data/raw/plato/r37/golden_fixture/
```

因此当前版本不能把任何 PLATO 导入用于正式计分。未来 fixture 完成独立来源审核后，必须
在同一个 reviewed change 中同步更新审批状态、provenance、production fingerprint
allowlist 和测试。`pred.html` 这类 JSON 响应本身没有
显式 query SMILES 字段，所以正式解释器仍要求额外满足查询结构绑定、每个药物
PLATO 原始导出 SHA-256 独立 allowlist、可靠人源蛋白靶点等门禁。手写 JSON、
手动填写 release、或者声明“这个文件来自 PLATO”，都不能打开正式计分门禁。

当前实现生成的结果（包括未来通过输入门禁的研究假说）会声明：

- `result_schema_version=5.4`
- `statistical_contract_version=1.1`
- `credibility.formula_version=evidence_credibility_v5_frozen_statistical_targets`
- `network_statistics.method_version=weighted_fixed_disease_closest_v3_mapped_universe`

旧 schema 3/4/5.0/5.1/5.2 以及旧公式 v1/v2/v3/v4 的结果只保留历史审计意义，必须从原始输入
按当前协议重新运行，不能通过改版本号升级，也不能与当前版本分数直接比较。
`evidence_credibility_score` 是未校准的证据评价规则分，不是靶点预测正确概率、临床疗效
概率或实验成功率；输入可审计、规则支持、敏感性和独立验证应分开报告。

## 一对药物：依次查询 PLATO 并进入机制解释

解释器的主任务是一对药物组合。低频调试时可以直接用 pair-level 脚本：

```powershell
python scripts/query_pair_and_explain.py `
  --drug1 "Drug-A" `
  --drug2 "Drug-B" `
  --smiles-lookup C:\datasets\drug_smiles.csv `
  --cell-line A549 `
  --disease-label pan_cancer_tier1 `
  --quick
```

如果不传 `--drug1-smiles / --drug2-smiles`，必须用 `--smiles-lookup PATH` 显式指定查询表；InterpreterPLATO 不分发 SMILES 表，也不会扫描工作区相邻项目。也可以直接传入两个 SMILES。解析测试可以先跑：

```powershell
python scripts/query_pair_and_explain.py `
  --drug1 "Gemcitabine" `
  --drug2 "Paclitaxel" `
  --smiles-lookup C:\datasets\drug_smiles.csv `
  --resolve-only
```

如果同一药名对应多个不同 SMILES，脚本会失败并要求你手动指定 SMILES。

流程是严格串行的：

1. 解析或使用 drug1/drug2 的精确 SMILES；
2. 对 drug1 向 PLATO prediction 表单发起一次 target-fishing JSON 请求；
3. 等待 `--delay-seconds`，默认 3 秒；
4. 对 drug2 发起一次请求；
5. 分别保存两个药物的 raw JSON、metadata sidecar、SHA-256 和 cache；
6. 校验调用方预先提供的独立可信 SHA-256 allowlist；匹配后才进入解释。

默认输出目录：

```text
InterpreterPLATO/outputs/research/plato_pair_runs/<timestamp>_<drug1>__<drug2>/
```

目录里会保存 `raw_plato/`、`plato_cache/` 和 `run_manifest.json`；只有外部
allowlist 等完整性门禁全部通过时才会生成研究机制假说 `explanation.json` 及
`explanation.json.sha256`。`pred.html` 是 PLATO
表单页面路径；脚本请求的是 `tipoout=json`，所以本地保存和解析的是 JSON 响应，
不是普通 HTML 页面。

pair 脚本不会把本次捕获的 SHA-256 自动当成可信 allowlist。同次计算的哈希只证明
保存后的字节未变化，不能证明响应来源。未传入或未匹配
`--trusted-plato-sha256` 时，脚本会保留两次查询的 raw、sidecar、cache 和 manifest，
随后以退出码 2 停止，不生成机制假说。首次采集建议先归档，再由独立审核记录批准
两个原始响应哈希。只有 fixture 审核完成并启用 production schema gate 后，才能使用现有
cache 和 `run_interpreter.py` 生成研究假说及未校准规则分，不能升级为正式机制结论。

## 单个药物：查询并归档一个 PLATO 结果

如果你希望程序替你完成网页表单的“单次查询”，使用：

```powershell
python scripts/query_plato_once.py `
  --drug-name "Compound-A" `
  --smiles "提交给 PLATO 的精确 SMILES" `
  --release r37 `
  --print-report
```

这个脚本的约束很硬：一次运行只发起一次 `POST`，只允许 PLATO 官方 target-fishing 表单端点，不循环、不并发、不做批量。它会同时保存：

- PLATO 原始 JSON 响应；
- 单次查询 metadata sidecar；
- 响应字节流 SHA-256；
- 提交表单、endpoint、HTTP 方法和 SMILES 结构绑定证据；
- 本地 PLATO cache v3 归档。

这种来源会被标记为 `query_binding_status=plato_single_query_post_verified`。它解决的是“这个本地 JSON 是否来自这次提交的 SMILES”这个审计问题；正式解释仍然需要两种药物各自的 PLATO cache、可信 SHA-256 allowlist、人源蛋白靶点和上下文 manifest 全部过门槛。

## 手动归档一个 PLATO 结果

先在 PLATO 网页对单个分子完成查询，并下载结果文件。然后用本地导入脚本归档：

```powershell
python scripts/import_plato_result.py `
  --drug-name "Compound-A" `
  --smiles "提交给 PLATO 的精确 SMILES" `
  --result "C:/path/plato-result.json" `
  --release r37 `
  --attest-downloaded-from-plato `
  --print-report
```

`--attest-downloaded-from-plato` 只记录用户声明：这个文件是用户从 PLATO 下载的。它不是 PLATO 数字签名，不是来源认证，也不是 schema 验证。

因此，导入成功只代表“本地归档成功”，不代表“可以正式计分”。只有 production-approved、
fingerprint 已进入 production allowlist 的 fixture（当前不存在）、结构绑定、release、原始
导出 SHA-256、可靠人源蛋白靶点等门禁全部通过后，才允许进入正式解释。

## 运行单对药物解释

在来源已审计的 fixture 获得生产批准、production allowlist 已启用，并且两个合格 PLATO
cache 都准备好之后，才可以运行：

```powershell
python scripts/run_interpreter.py summary `
  --context-config data/production_context.json `
  --plato-cache data/processed/plato_cache `
  --plato-release r37 `
  --trusted-plato-sha256 "<drug-a-raw-export-sha256>" `
  --trusted-plato-sha256 "<drug-b-raw-export-sha256>" `
  --drug1 "Drug-A" --drug2 "Drug-B" `
  --drug1-smiles "SMILES-A" --drug2-smiles "SMILES-B" `
  --cell-line A549 --disease-label pan_cancer_tier1 `
  --output outputs/research/drug-a__drug-b.json
```

CLI 会在正式运行前预检输入和缓存状态。以下情况会直接失败，并且不会写出正式分数：

- PLATO cache 缺失。
- cache schema 不是当前版本 `3`。
- r37 fixture 未验证。
- 查询 SMILES 和 payload 结构不一致。
- 缺少独立可信的 PLATO 原始导出 SHA-256。
- batch JSON 根类型或字段非法。
- 输入中出现 `database`、`local_targets`、`drug1_targets`、`drug2_targets` 等本地靶点字段。

Python 公共入口是：

```python
from InterpreterPLATO import ExplanationPipeline
```

生产调用不接受 `drug1_targets`、`drug2_targets` 或 SQLite 靶点库。包顶层也不再导出本地靶点 resolver 或 local/PLATO 融合工具。

## 批量解释

批量输入必须是 JSON list，每一项至少包含：

```text
drug1_name, drug2_name, drug1_smiles, drug2_smiles
```

允许的可选字段包括：

```text
cell_line, predicted_synergy, advanced, disease_label, model_fidelity_evidence
```

批量模式把单对规则分与批量统计支持分开：

- `base_score` 保留单对基础分；兼容字段 `batch_adjusted_score` 与本地分数相等，不再施加固定扣分。
- `score_profile.name=single_pair`；`batch_FDR_applied=true` 表示做了批量检验，
  `batch_FDR_adjustment_applied=false` 表示没有把它改写成分数调整。
- `single_pair_mechanism_supported` 与 `batch_mechanism_supported` 分开保存；非网络主机制的
  后者为 `null`，不能用药物－疾病邻近性检验替它背书。
- `network_statistics.batch_q_A/B`、`batch_required_sides`、`batch_pattern_status` 和
  `batch_family` 保存校正结果、模式所需药物侧及检验家族身份。

每个传入批次的 N 对结果保留全部 2N 个计划药物－疾病检验；无效或不可评估的检验按
`p=1` 留在家族中且不能建立支持。双药邻近模式要求两侧均满足支持条件，单药暴露模式仅
要求原本邻近的那一侧。检验家族记录药物结构、来源哈希、固定靶点、上下文 artifact 和方法，
不以观察到的显著性定义身份。缺失身份字段会被显式标记，不能当作完整审计。

本函数只校正传入批次，不自动覆盖整篇论文。研究者必须预先冻结完整研究家族，不能拆分
批次挑选显著结果。BH 的独立性或适当正相关条件及当前自定义零模型的可交换性未经验证，
因此这里报告的是探索性的条件零模型支持，不宣称已经实现校准的总体错误率控制；
不显著也不等于证明没有机制或没有协同。详见[统计合同](docs/STATISTICAL_METHODS.md)。

## 数据与上下文

当前本地已经能加载并校验以下癌症上下文层：

- PPI：STRING human PPI。
- 表达：DepMap 表达矩阵。
- 通路：Reactome。
- 疾病集合：COSMIC Cancer Gene Census。
- 基因映射：固定版本 HGNC。
- 靶点 ID 翻译：ChEMBL r37 target ID → HGNC symbol。

这些数据只参与癌症上下文重排、ID 翻译和机制解释，不参与靶点补全。

ChEMBL r37 在这里不是第二个靶点预测来源。它只做：

```text
PLATO target_chembl_id → ChEMBL r37 target/component → HGNC 2026-04-01 symbol
```

也就是把 PLATO 已经给出的 `CHEMBL target id` 翻译成 PPI、DepMap、Reactome、COSMIC 能识别的人类基因符号。构建命令：

```powershell
python scripts/build_chembl_target_mapping.py
```

输入原始文件位置：

```text
InterpreterPLATO/data/raw/target_mapping/chembl/r37/
```

输出：

```text
InterpreterPLATO/data/processed/context/chembl_r37_target_hgnc_mapping.tsv.gz
InterpreterPLATO/data/processed/context/chembl_r37_target_hgnc_mapping_qc.json
```

当前策略只把唯一单组件 target 直接补成一个 HGNC symbol；多组件复合体/家族先保守跳过直接映射，避免把一个复合体误压成单基因机制。

主要配置文件：

- `InterpreterPLATO/data/production_context.json`
- `InterpreterPLATO/data/runtime.lock.json`（默认快速校验范围）
- `InterpreterPLATO/data/source.lock.json`（约 35.45 GiB 原始快照）
- `InterpreterPLATO/data/manifest.lock.json`（详细来源与许可证元数据的兼容清单）
- `InterpreterPLATO/data/processed/context/`

历史本地 demo 靶点文件已经归档到：

```text
InterpreterPLATO/data/processed/archive/pre_plato_only/
```

旧输出也已经归档到：

```text
InterpreterPLATO/outputs/archive/pre_plato_only/
InterpreterPLATO/outputs/archive/pre_research_contract_20260718/
```

这些文件只保留审计历史，不能作为当前生产结果或回归基线。

## 可复现性

`pyproject.toml` 约束 Python 3.9--3.12；`requirements.lock` 和 `requirements-dev.lock` 固定本次审核环境的依赖版本。RDKit 是必需依赖，因为严格结构绑定依赖 RDKit 完成 SMILES 解析、规范化和 InChIKey 生成。每项数据的官方来源、当前 release、获取限制和 checksum 证据见[数据来源与获取说明书](docs/DATA_SOURCES.md)；首次准备当前锁定数据见 [Data bootstrap](docs/DATA_BOOTSTRAP.md)；升级数据版本、重建、审核、更新 locks 和回滚见[数据更新操作手册](docs/DATA_UPDATE_GUIDE.md)。

研究运行的 `audit` 会自动记录关键版本与哈希；论文归档时还建议同时保存：

- Python 版本。
- 解析后的依赖版本，例如 `pip freeze`。
- 代码版本。
- 上下文 manifest。
- PLATO cache schema/parser 版本。
- 所有输入和输出文件的 SHA-256。

常用检查命令：

```powershell
python scripts/build_context_artifacts.py all
python scripts/verify_manifest.py --scope runtime
python scripts/verify_manifest.py --scope source
python -m pytest tests
```

## 当前最重要的限制

程序记录了本地 Gemcitabine r37 fixture candidate 的 schema fingerprint，但 production
allowlist 当前为空，所以 `schema_fixture_verified=false`，正式评分会 fail closed。fixture
独立审核并获批后，还必须在同一 reviewed change 中更新 manifest/provenance、代码 allowlist
和测试；两个药物也要分别通过 query binding 与独立 raw-SHA allowlist。runtime lock 只
校验 fixture 字节身份，不能替代这些信任门禁。
