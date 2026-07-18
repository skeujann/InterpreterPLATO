# InterpreterPLATO 5.1

InterpreterPLATO 是一个药物组合生物学解释器。它的作用不是重新训练一个靶点预测模型，而是把一对药物组合的靶点证据、癌症细胞系上下文、PPI 网络、通路和疾病基因集合整合起来，生成可审计的机制解释。

当前版本采用严格的 `PLATO-only` 靶点策略：

- 两个药物的候选靶点都必须来自 PLATO 单分子 target-fishing 查询结果。
- 每个药物都需要单独在 PLATO 平台查询、单独下载结果、单独本地归档。
- 本地 PPI、DepMap 表达、Reactome 通路、COSMIC 疾病基因集合只用于癌症上下文重排和机制解释。
- 本地上下文数据不能补充、替代或伪造药物靶点。
- 程序不做批量访问、不爬取 PLATO 网页；如需程序查询，只允许用户触发的单次 PLATO 表单请求，并保存完整审计证据。

也就是说，这里的 `PLATO-only` 指“正式药物靶点来源只认 PLATO”，不是指本程序拥有独立靶点预测模型、公开 API 授权或批量查询能力；PLATO 仍然是外部平台。

## 当前生产状态

正式结果必须同时满足以下条件：

- PLATO 缓存 schema 为 `3`。
- 两个药物的原始 PLATO 响应字节流都被不可变保存，并记录 SHA-256。
- PLATO 响应中的查询结构必须能和提交给 PLATO 的 SMILES 对上。
- 结构绑定依赖 RDKit 解析、规范 SMILES 和 InChIKey。
- 响应 payload 中必须能验证 release 和 target-fishing 方法。
- 响应 schema 必须通过受控、真实、脱敏的 PLATO r37 golden fixture 验证。
- 两个药物都必须保留至少一个可靠、明确人源、蛋白类靶点。
- PPI、表达、通路、疾病等上下文文件必须通过配置中的 SHA-256 校验。
- 两个 PLATO 原始导出文件的 SHA-256 都必须来自独立可信 allowlist。

本地审计数据包可以放置 Gemcitabine 的 PLATO r37 target-fishing JSON 响应，作为
解析器的 golden fixture，用来确认 `chemblId / nome / score / validato / simili`
这套真实字段结构。代码仓库不提交原始响应或 PDF；本地放置位置为：

```text
InterpreterPLATO/data/raw/plato/r37/golden_fixture/
```

但这不等于所有 PLATO 导入都可以正式计分。`pred.html` 这类 JSON 响应本身没有
显式 query SMILES 字段，所以正式解释器仍要求额外满足查询结构绑定、每个药物
PLATO 原始导出 SHA-256 独立 allowlist、可靠人源蛋白靶点等门禁。手写 JSON、
手动填写 release、或者声明“这个文件来自 PLATO”，都不能打开正式计分门禁。

未来通过全部门禁的正式结果会声明：

- `result_schema_version=5.1`
- `credibility.formula_version=evidence_credibility_v4_verified_plato_only`

旧 schema 3/4/5.0、旧公式 v1/v2/v3 的结果只是历史审计产物，不能和 v5.1 分数直接比较。`evidence_credibility_score` 表示证据完整性和可审计性，不是靶点预测正确概率、临床疗效概率或实验成功率。

## 一对药物：依次查询 PLATO 并进入机制解释

解释器的主任务是一对药物组合。低频调试时可以直接用 pair-level 脚本：

```powershell
python InterpreterPLATO/scripts/query_pair_and_explain.py `
  --drug1 "Drug-A" `
  --drug2 "Drug-B" `
  --cell-line A549 `
  --disease-label pan_cancer_tier1 `
  --quick
```

如果不传 `--drug1-smiles / --drug2-smiles`，必须用 `--smiles-lookup PATH` 显式指定查询表；InterpreterPLATO 不分发 SMILES 表，也不会扫描工作区相邻项目。也可以直接传入两个 SMILES。解析测试可以先跑：

```powershell
python InterpreterPLATO/scripts/query_pair_and_explain.py `
  --drug1 "Gemcitabine" `
  --drug2 "Paclitaxel" `
  --resolve-only
```

如果同一药名对应多个不同 SMILES，脚本会失败并要求你手动指定 SMILES。

流程是严格串行的：

1. 解析或使用 drug1/drug2 的精确 SMILES；
2. 对 drug1 向 PLATO prediction 表单发起一次 target-fishing JSON 请求；
3. 等待 `--delay-seconds`，默认 3 秒；
4. 对 drug2 发起一次请求；
5. 分别保存两个药物的 raw JSON、metadata sidecar、SHA-256 和 cache；
6. 用这两个 cache 进入 PLATO-only 协同机制解释。

默认输出目录：

```text
InterpreterPLATO/outputs/plato_pair_runs/<timestamp>_<drug1>__<drug2>/
```

目录里会有 `raw_plato/`、`plato_cache/`、`explanation.json` 和 `run_manifest.json`。`pred.html` 是 PLATO 表单页面路径；脚本请求的是 `tipoout=json`，所以本地保存和解析的是 JSON 响应，不是普通 HTML 页面。

pair 脚本会把本次捕获的两个原始响应 SHA-256 作为本次运行 allowlist，并在 `run_manifest.json` 中记录 `trusted_export_sha256_source=same_run_single_query_capture`。这适合低频调试和可追踪运行；如果要做正式论文级审计，建议把两个 SHA-256 固化到独立清单后再复跑解释。

## 单个药物：查询并归档一个 PLATO 结果

如果你希望程序替你完成网页表单的“单次查询”，使用：

```powershell
python InterpreterPLATO/scripts/query_plato_once.py `
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
python InterpreterPLATO/scripts/import_plato_result.py `
  --drug-name "Compound-A" `
  --smiles "提交给 PLATO 的精确 SMILES" `
  --result "C:/path/plato-result.json" `
  --release r37 `
  --attest-downloaded-from-plato `
  --print-report
```

`--attest-downloaded-from-plato` 只记录用户声明：这个文件是用户从 PLATO 下载的。它不是 PLATO 数字签名，不是来源认证，也不是 schema 验证。

因此，导入成功只代表“本地归档成功”，不代表“可以正式计分”。只有 controlled fixture、结构绑定、release、原始导出 SHA-256、可靠人源蛋白靶点等门禁全部通过后，才允许进入正式解释。

## 运行单对药物解释

在 verified fixture 和两个合格 PLATO cache 都准备好之后，可以运行：

```powershell
python InterpreterPLATO/scripts/run_interpreter.py summary `
  --context-config InterpreterPLATO/data/production_context.json `
  --plato-cache InterpreterPLATO/data/processed/plato_cache `
  --plato-release r37 `
  --trusted-plato-sha256 "<drug-a-raw-export-sha256>" `
  --trusted-plato-sha256 "<drug-b-raw-export-sha256>" `
  --drug1 "Drug-A" --drug2 "Drug-B" `
  --drug1-smiles "SMILES-A" --drug2-smiles "SMILES-B" `
  --cell-line A549 --disease-label pan_cancer_tier1 `
  --output InterpreterPLATO/outputs/drug-a__drug-b.json
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

批量模式会记录单对基础分和批量 FDR 调整后的分数：

- `base_score`：单对药物组合在独立解释时的基础分。
- `batch_adjusted_score`：批量 Benjamini-Hochberg FDR 调整后的分数。
- `score_profile.name`：`single_pair` 或 `batch_BH_FDR`。

这样同一对药物在单独运行和批量运行时，即使分数不同，也能追踪差异来源。

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
python InterpreterPLATO/scripts/build_chembl_target_mapping.py
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
```

这些文件只保留审计历史，不能作为当前生产结果或回归基线。

## 可复现性

`pyproject.toml` 约束 Python 3.9--3.12；`requirements.lock` 和 `requirements-dev.lock` 固定本次审核环境的依赖版本。RDKit 是必需依赖，因为严格结构绑定依赖 RDKit 完成 SMILES 解析、规范化和 InChIKey 生成。完整的数据获取、外部 ScoreRI 输入和离线重建流程见 `docs/DATA_BOOTSTRAP.md`。

正式运行时建议同时记录：

- Python 版本。
- 解析后的依赖版本，例如 `pip freeze`。
- 代码版本。
- 上下文 manifest。
- PLATO cache schema/parser 版本。
- 所有输入和输出文件的 SHA-256。

常用检查命令：

```powershell
python InterpreterPLATO/scripts/build_context_artifacts.py all
python InterpreterPLATO/scripts/verify_manifest.py --scope runtime
python InterpreterPLATO/scripts/verify_manifest.py --scope source
python -m pytest InterpreterPLATO/tests
```

## 当前最重要的限制

现在还不能把任意 PLATO-like JSON 当成真实 PLATO r37 输出。下一步要真正进入生产计分，需要获得一份真实、可审计、脱敏后的 PLATO r37 导出样本，用它固定 schema fingerprint，并把每个药物原始导出的 SHA-256 写入可信 allowlist。
