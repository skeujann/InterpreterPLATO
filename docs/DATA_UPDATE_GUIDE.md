# InterpreterPLATO 数据更新操作手册

本文说明如何在不把普通原始数据提交到 Git 的前提下，更新 InterpreterPLATO 的癌症
上下文、ChEMBL→HGNC 靶点映射、ScoreRI 细胞系清单和 PLATO 审计数据；当前已跟踪的
小型受控 fixture 例外会在下文单独说明。

每项数据的官方入口、当前锁定 release、容器/解压关系、访问限制和 checksum 证据等级
统一记录在[数据来源与获取说明书](DATA_SOURCES.md)。本手册只描述升级操作；获取新文件
前应先核对来源说明书，不能仅凭同名文件或本地 SHA 判断上游身份。

所有命令默认从 `InterpreterPLATO` 项目根目录运行：

```powershell
Set-Location C:\path\to\ScoredRI\InterpreterPLATO
```

仓库以代码、构建脚本、配置和完整性清单为主。`data/raw/`、`data/processed/`、
`outputs/`、普通药物的 PLATO 响应和模型文件都应留在本地数据包，不得因为更新数据而
提交到 Git。当前仓库有且仅有一个经过显式放行的受控数据例外：
`data/raw/plato/r37/golden_fixture/` 下的 fixture manifest 和 100280-byte target-fishing
JSON；两者是当前 parser schema anchor candidate，并由 source/runtime locks 固定，但这
不等于已通过生产来源审计。production schema allowlist 当前为空，该 candidate 不能打开
正式计分；每药物 raw-SHA allowlist 是审核通过后仍需满足的另一道门禁。PDF、其他查询响应和受限
原件不随代码提交。受限原件与可发布的脱敏派生副本必须分开管理；新增或替换例外
必须重新经过许可、脱敏和独立发布审查，不能把普通查询数据仿照该目录提交。

## 1. 先理解数据边界

```text
官方或授权来源
    ↓ 手工下载、记录来源；有则校验上游 checksum，否则明确证据等级
data/raw/<来源>/<版本>/                 不可变源快照
    ↓ 离线 builder，先写 candidate 目录
data/processed/candidates/<更新批次>/   候选产物和 QC
    ↓ 人工比较覆盖率、拒绝项和生物学差异
data/processed/releases/<更新批次>/     已审核运行时包
    ↓ 更新 production_context.json
data/runtime.lock.json                  生产加载完整性门禁

每个药物的 PLATO 单次查询
    ↓ 保存原始响应、sidecar、SHA-256
data/raw/... + data/processed/plato_cache/
    ↓ 独立 raw-SHA allowlist + query binding + 已审计 fixture 的 production schema allowlist
正式机制计分（当前 schema allowlist 为空，因此 fail closed）
```

各文件的职责如下：

| 文件或目录 | 职责 | 是否自动生成 |
|---|---|---|
| `data/source_registry.json` | 当前重建所需的 12 个输入、默认路径、大小和 SHA-256 | 否，人工审核更新 |
| `data/source.lock.json` | 完整 `data/raw` 快照和小型构建注册表的字节身份 | 否，人工审核更新 |
| `data/manifest.lock.json` | 来源、可用的上游 checksum、许可证、证据等级和构建元数据 | 否，人工审核更新 |
| `data/production_context.json` | 生产实际加载的 5 个运行时资源 | 否，人工审核更新 |
| `data/runtime.lock.json` | 生产配置、5 个运行时资源和受控 PLATO fixture 的字节身份 | 否，人工审核更新 |
| `scripts/verify_manifest.py` | 只读校验现有 lock | 从不改文件 |
| `scripts/bootstrap_data.py plan/verify` | 只读列出和校验 registry 输入 | 从不下载、从不改 lock |
| `scripts/bootstrap_data.py rebuild` | 复现当前 pin，直接写 live `processed/context` | 会改产物，不改 lock |

项目目前没有 lock 生成器。不要把一次构建观察到的哈希直接抄进 lock 来“消除报错”；
必须先确认来源、版本、构建代码和 QC 差异都符合预期。

## 2. 区分四种更新

| 更新类型 | 正确入口 |
|---|---|
| 恢复或复现当前已锁版本 | `bootstrap_data.py verify/rebuild` |
| HGNC、STRING、DepMap、Reactome、COSMIC 或 ChEMBL 升级 | 新版本目录 → 更新代码 pin/registry → candidate 构建 → QC → config/locks |
| ScoreRI 细胞系集合变化 | 重建 `scoreri_cell_lines.tsv` → 审核别名 → 重建 DepMap 表达层 |
| 新增某个药物的 PLATO 响应 | 单次查询并归档；通常不修改 context locks |

`bootstrap_data.py rebuild` 只适合复现当前版本，不适合直接升级：它先用旧
`source_registry.json` 校验全部 12 个输入，构建后又用旧 `runtime.lock.json`
校验输出。升级数据时，旧 pin 会阻止新源，旧 runtime lock 也必然拒绝新产物。
此外，即使指定 `--component context` 或 `target-mapping`，wrapper 仍会先校验
全部 12 个 registry 输入。

## 3. 建立更新前基线

建议使用 Python 3.9–3.12，并安装锁定依赖：

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install --requirement requirements.lock
.venv\Scripts\python -m pip install --no-deps --editable .
. .\.venv\Scripts\Activate.ps1
python -c "import sys; print(sys.executable)"
```

激活命令会在当前 PowerShell 会话中启用虚拟环境，最后一行打印实际解释器供核对；后续
文中的 `python` 都应指向
`.venv\Scripts\python.exe`。如果执行策略禁止激活，则不要跳过环境绑定：把后续每个
`python` 明确替换为 `.venv\Scripts\python`。

更新前先保存当前代码版本、配置和验证结果：

```powershell
git status --short
git rev-parse HEAD
python scripts\verify_manifest.py --scope runtime
python scripts\bootstrap_data.py plan
python scripts\bootstrap_data.py verify
```

完整源快照约 35 GiB，下面的校验会读取整个数据包，应预留时间和磁盘 I/O：

```powershell
python scripts\verify_manifest.py --scope source
python scripts\verify_manifest.py --manifest data\manifest.lock.json
```

`bootstrap_data.py verify` 只检查当前构建所需的 12 个 registry 输入；
`verify_manifest.py --scope source` 检查 `source.lock.json` 记录的完整 raw 快照，
两者不是等价命令。`--scope all` 只包含 runtime 和 source，不会自动校验
`manifest.lock.json`，后者必须通过 `--manifest` 显式选择。

verifier 只检查“已经列入 lock 的文件”，不会发现 raw 目录中新出现但尚未登记的文件。
用下面的只读比较检查 source lock 与实际 raw 路径集合；无输出才表示一致：

```powershell
$root = (Resolve-Path .).Path
$lock = Get-Content data\source.lock.json -Raw | ConvertFrom-Json
$locked = @($lock.artifacts.relative_path | Sort-Object)
$actual = @(Get-ChildItem data\raw -Recurse -File |
  ForEach-Object { $_.FullName.Substring($root.Length + 1).Replace('\','/') } |
  Sort-Object)
Compare-Object $locked $actual
```

升级前还应保留旧 `data/processed` 运行时包。不要在唯一一份 live 目录上试构建。

## 4. 当前五层生产上下文是怎么做出来的

| 生产资源 | 原始输入 | 构建器 | 运行时输出 |
|---|---|---|---|
| PPI | HGNC + STRING human info/full links | `build_context_artifacts.py all` | `string_v<STRING>_hgnc_<HGNC>.tsv.gz` |
| 表达 | HGNC + DepMap expression/Model/Profiles + ScoreRI cell-line registry | 同上 | `depmap_<RELEASE>_scoreri_log2tpm.csv.gz` |
| 通路 | HGNC + Reactome GMT | 同上 | `reactome_v<RELEASE>_hgnc_<HGNC>.gmt` |
| 疾病集 | HGNC + COSMIC Cancer Gene Census | 同上 | `cosmic_cgc_<RELEASE>_pan_cancer_tier1.txt` |
| 靶点 ID 翻译 | ChEMBL SQLite + HGNC | `build_chembl_target_mapping.py` | `chembl_<RELEASE>_target_hgnc_mapping.tsv.gz` |

ChEMBL 只把 PLATO 已给出的 ChEMBL target ID 翻译为 HGNC symbol，不能作为第二个
药物靶点来源。

当前锁定基线为 HGNC `2026-04-01`、STRING `v12.0`、DepMap expression
`26Q1`、Reactome `v96`、COSMIC CGC `v101`、ChEMBL `r37`，以及 81 个
ScoreRI 细胞系名称。DepMap 的 expression、`Model.csv` 和 `OmicsProfiles.csv`
已统一落在 `data/raw/cell_line_context/depmap/26Q1/`；其中 metadata 两个文件于
2026-07-23 从 26Q1 release 页面重新取得，size/SHA 与先前锁定快照一致。

主要转换规则：

| 层 | 规范化与过滤规则 |
|---|---|
| HGNC | 唯一 Entrez → approved symbol 精确匹配 → 唯一 alias/previous symbol；歧义项拒绝，不猜测 |
| STRING | 只处理 human 文件，`combined_score >= 700`；HGNC 映射后拒绝未映射端点和 self-loop，重复边取最大分 |
| DepMap | 基因统一映射 HGNC；细胞系只接受唯一匹配或显式审核 alias；只取可验证 default profile；保留源 `log2(TPM+1)` 值 |
| Reactome | gene 统一映射 HGNC，拒绝未映射/歧义成员，去重 membership，映射后空 pathway 不进入运行时 GMT |
| COSMIC | gene 统一映射 HGNC，拒绝歧义项；生产疾病集只取审核后的 pan-cancer Tier 1 symbols |
| ChEMBL target map | 按 `GENE_SYMBOL` → `GENE_SYMBOL_OTHER` → 唯一 UniProt 映射 HGNC；正式直接映射只接受唯一单组件 target |

当前输入的最低格式要求：

- STRING info 为 TSV，包含 `#string_protein_id`、`preferred_name`；links 包含
  `protein1 protein2 combined_score`。默认阈值为 700/1000，映射后的重复边保留最大分。
- ScoreRI registry 为 TSV/TXT 或 CSV，必须含 `cell_line_name`。
- DepMap `Model.csv` 必须包含 `ModelID`、`CellLineName`、
  `StrippedCellLineName`、`CCLEName`；expression 必须包含
  `SequencingID`、`ModelID`、`IsDefaultEntryForModel`；Profiles 必须包含
  `ModelID`、`SequencingID`、`DataType`、`IsDefaultEntryForModel`。三份文件必须
  来自同一 release。
- Reactome GMT 每行至少为 pathway name、stable ID、gene symbols 三列。
- COSMIC 必须包含 `GENE_SYMBOL`、`TIER`、`ROLE_IN_CANCER`、
  `TUMOUR_TYPES_SOMATIC`；Tier 1 结果不得为空。COSMIC 是许可受限数据。
- ChEMBL `checksums.txt` 必须是包含 `file`、`SHA256` 两列的 tab 分隔表，并能按
  archive 文件名找到 SHA-256；SQLite schema 也必须与 builder 查询兼容。

## 5. 获取新版本源数据

先按[数据来源与获取说明书](DATA_SOURCES.md)确认准确上游文件、账号/许可边界、容器
留存方式和当前 provenance 缺口。特别注意：STRING 当前 registry 锁定 gunzip 后文本；
Reactome、COSMIC 和 ChEMBL 当前构建也需要预先准备各自的提取产物。builder 和
`bootstrap_data.py` 都不会下载数据，也不会替新 release 生成可信来源记录。

先查看当前预期文件、官方入口和再分发提示：

```powershell
python scripts\bootstrap_data.py plan
```

这个命令不会联网下载。应由有权限的用户从官方版本页手工获取数据，并执行：

1. 把新版本放进新的版本目录，不覆盖旧 raw 文件；不要使用无版本的 `latest/`。
2. 保存官方 README、LICENSE、required attribution 和 checksum 文件。
3. 对容器文件和最终输入分别留档，例如同时保存 ChEMBL archive 和解压后的 SQLite。
4. 不把账号、Cookie、临时签名 URL 或访问令牌写入 registry、manifest 或 README。
5. COSMIC 等许可受限数据只在授权环境中下载和使用。

单文件的本地身份可用 PowerShell 只读计算：

```powershell
$path = "data\raw\<source>\<version>\<file>"
(Get-Item -LiteralPath $path).Length
(Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant()
```

如果上游发布了 checksum，必须先独立验证上游 checksum，再把本地 SHA-256 记录进
`source_registry.json` 和 provenance manifest。只有“本地文件自己算出的哈希”并不能
证明它来自正确上游。

`--input ID=PATH` 只覆盖路径，不改变 registry 中的 expected size/SHA。因此它只适合
验证“当前已 pin 的同一字节副本被移到其他磁盘”，不能用旧 registry 验证新 release：

```powershell
python scripts\bootstrap_data.py verify `
  --input hgnc=D:\data\hgnc_complete_set.txt `
  --input string_links=D:\data\9606.protein.links.full.txt
```

新 release 必须先建立并人工审核包含新 size/SHA 的 candidate registry，再验证全部
12 项：

```powershell
python scripts\bootstrap_data.py verify `
  --registry data\source_registry.candidate.json `
  --input hgnc=D:\data\new-release\hgnc_complete_set.txt
```

candidate registry 通过并完成后续审核，才可提升为正式
`data/source_registry.json`。

## 6. 升级 release 时同步更新代码 pin

不能只给 builder 传一个新路径，却继续把输出标成旧 release。至少审核：

- `scripts/build_context_artifacts.py` 中的 `HGNC_RELEASE`、
  `HGNC_EXPECTED_SHA256`、`STRING_RELEASE`、`DEPMAP_RELEASE`、
  `REACTOME_RELEASE`、`COSMIC_RELEASE`、默认输入路径和派生输出名。
- `scripts/build_chembl_target_mapping.py` 中的 `HGNC_RELEASE`、
  `CHEMBL_RELEASE`、默认 archive/SQLite/checksum 路径和输出名。
- `data/source_registry.json` 中相应输入的 path、source（当前 schema 把版本写在
  source 文本和部分路径中，并无独立 version 字段）、size、SHA-256、acquisition URL
  和许可说明。
- 测试里固定的 release、文件名、哈希和生物学基线。

CLI 始终启用 HGNC 硬哈希检查；更新 HGNC 时必须更新
`HGNC_EXPECTED_SHA256`。HGNC 一旦升级，PPI、表达、Reactome、COSMIC 和
ChEMBL target map 必须使用同一 HGNC 快照整批重建，不能混用 mapping version。

ChEMBL builder 自身只用 `checksums.txt` 校验 archive；它不会认证
`checksums.txt`，也不会把 HGNC 或解压后的 SQLite 与 registry 自动比对。因此构建前
的 candidate registry 校验与官方 archive checksum 校验不能省略，正式构建也不应使用
`--skip-archive-check`；完整 source lock 在 QC 通过后按第 11 节提升。

## 7. 必要时更新 ScoreRI 细胞系清单

如果 ScoreRI 组合表或细胞系集合变化，必须在 candidate 构建前更新
`data/scoreri_cell_lines.tsv`。下面的 PowerShell 以 ordinal 规则 trim、去空、去重、
排序，并写出无 BOM 的 UTF-8/LF 文件：

```powershell
$combo = "D:\reviewed-data\ScoreRI\combo_data.tsv"
$output = "data\scoreri_cell_lines.tsv"
$rows = @(Import-Csv -LiteralPath $combo -Delimiter "`t" -Encoding UTF8)
if ($rows.Count -eq 0 -or
    $rows[0].PSObject.Properties.Name -notcontains "cell_line_name") {
  throw "ScoreRI table is empty or missing cell_line_name: $combo"
}
$names = @($rows |
  ForEach-Object { ([string]$_.cell_line_name).Trim() } |
  Where-Object { $_ })
if ($names.Count -eq 0) { throw "no non-empty cell_line_name values: $combo" }
$unique = [System.Collections.Generic.HashSet[string]]::new(
  [System.StringComparer]::Ordinal
)
foreach ($name in $names) { [void]$unique.Add($name) }
$ordered = [string[]]$unique
[System.Array]::Sort($ordered, [System.StringComparer]::Ordinal)
$lines = [string[]](@("cell_line_name") + $ordered)
$outputPath = [System.IO.Path]::GetFullPath((Join-Path (Get-Location) $output))
[System.IO.File]::WriteAllText(
  $outputPath,
  (($lines -join "`n") + "`n"),
  [System.Text.UTF8Encoding]::new($false)
)

$ordered.Count
(Get-FileHash -LiteralPath $combo -Algorithm SHA256).Hash.ToLowerInvariant()
(Get-Item -LiteralPath $outputPath).Length
(Get-FileHash -LiteralPath $outputPath -Algorithm SHA256).Hash.ToLowerInvariant()
```

随后：

1. 记录新版组合表 SHA-256 和唯一名称数。
2. 更新 candidate source registry 中 cell-line registry 的 size/SHA。
3. 记录候选 `source.lock.json` 需要更新的 build registry、`derived_from` 和汇总字段；
   此时不要改正式 lock，等 QC 通过后再按第 11 节一次性提升。
4. 重新审核 `EXPLICIT_SCORE_CELL_LINE_ALIASES` 和
   `EXPLICIT_SCORE_CELL_LINE_REJECTIONS`，不能用字符串相似度猜 DepMap 模型。
5. 用新的 `--expected-score-cell-lines <N>` 重建表达层并检查每个拒绝原因。

如果唯一名称数不再是 81，还要同步更新 wrapper、测试和文档中的预期值。DepMap
expression、Model、Profiles 必须来自同一 release；当前基线三份文件均位于 26Q1
版本目录，下一次升级也必须先补齐同 release 的三份输入和 provenance 后再构建。

代码 pin、candidate registry 和可选的 ScoreRI 清单全部定稿后，在构建前验证新
registry 的全部 12 项：

```powershell
python scripts\bootstrap_data.py verify `
  --registry data\source_registry.candidate.json
```

## 8. 在 candidate 目录离线构建

新版本升级不要先写 `data/processed/context`。创建一个不会被生产配置引用的目录：

```powershell
$releaseTag = "2026-07-19_review"
$candidate = "data\processed\candidates\$releaseTag"
if (Test-Path -LiteralPath $candidate) {
  throw "candidate already exists; choose a new releaseTag: $candidate"
}
New-Item -ItemType Directory -Path $candidate
```

先把下面所有占位符替换成 candidate registry 已审核的路径；不要照抄旧 release 路径：

```powershell
$hgnc = "data\raw\gene_mapping\hgnc\<HGNC_RELEASE>\<HGNC_FILE>"
$stringLinks = "data\raw\ppi\string\<STRING_LINKS_FILE>"
$stringInfo = "data\raw\ppi\string\<STRING_INFO_FILE>"
$depmapExpression = "data\raw\cell_line_context\depmap\<DEPMAP_RELEASE>\<EXPRESSION_FILE>"
$depmapModel = "data\raw\cell_line_context\depmap\<DEPMAP_RELEASE>\Model.csv"
$depmapProfiles = "data\raw\cell_line_context\depmap\<DEPMAP_RELEASE>\OmicsProfiles.csv"
$reactome = "data\raw\pathways\reactome\<REACTOME_GMT>"
$cosmic = "data\raw\disease\cosmic\<COSMIC_RELEASE>\<COSMIC_TSV>"
$chemblRoot = "data\raw\target_mapping\chembl\<CHEMBL_RELEASE>"
$chemblSqlite = "$chemblRoot\<CHEMBL_SQLITE>"
$chemblArchive = "$chemblRoot\<CHEMBL_ARCHIVE>"
$chemblChecksums = "$chemblRoot\checksums.txt"
$expectedCellLines = 81  # 改成已审核的 registry 唯一名称数
```

然后运行完整 context 构建：

```powershell
python scripts\build_context_artifacts.py all `
  --output-dir $candidate `
  --hgnc $hgnc `
  --string-links $stringLinks `
  --string-info $stringInfo `
  --string-score-threshold 700 `
  --depmap-expression $depmapExpression `
  --depmap-model $depmapModel `
  --depmap-profiles $depmapProfiles `
  --score-combos data\scoreri_cell_lines.tsv `
  --expected-score-cell-lines $expectedCellLines `
  --reactome $reactome `
  --cosmic $cosmic
```

随后把 ChEMBL target mapping 也写入同一个 candidate 目录：

```powershell
python scripts\build_chembl_target_mapping.py `
  --sqlite $chemblSqlite `
  --hgnc $hgnc `
  --archive $chemblArchive `
  --checksums $chemblChecksums `
  --output "$candidate\chembl_<CHEMBL_RELEASE>_target_hgnc_mapping.tsv.gz" `
  --qc "$candidate\chembl_<CHEMBL_RELEASE>_target_hgnc_mapping_qc.json"
```

建议一次运行 context 的 `all`。如果连续运行 `ppi`、`depmap` 等子命令，
`build_qc.json` 最终只描述最后一次命令，应把每次 QC 另行保留。部分 gzip 输出采用
`mtime=0` 和逐文件原子替换，但整套 bundle 不是事务；中途失败仍可能留下新旧混合文件，
这也是必须使用 candidate 目录的原因。

复现当前 pin 时才使用 wrapper：

```powershell
python scripts\bootstrap_data.py rebuild --component all
```

它直接写 live `data/processed/context`，且没有 candidate `--output-dir`、STRING
threshold 或 expected-cell-count 透传能力。

## 9. 审核 QC，而不只是审核哈希

至少检查并与上一版逐项比较：

| QC 文件 | 必看指标 |
|---|---|
| `hgnc_mapper_qc.json` | approved symbols、唯一/歧义 alias、HGNC 输入 SHA |
| `string_ppi_qc.json` | 节点、边、最大连通分量、阈值过滤、未映射和重复边 |
| `depmap_expression_qc.json` | 请求细胞系数、默认 profile 可用数、缺失/拒绝项、基因数 |
| `reactome_qc.json` | pathway 数、membership 数、未映射/歧义基因、空 pathway |
| `cosmic_qc.json` | 输入/保留行、Tier 1 gene 数、拒绝项 |
| `chembl_*_qc.json` | SQLite/HGNC SHA、mapped/unmapped targets、唯一单组件 targets |

同时审阅所有 mapping/rejection sidecar。当前基线中 ScoreRI 请求 81 个细胞系，
70 个有可验证默认表达 profile，`SR` 被明确拒绝为“不可猜测”；这类声明缺口不能在
升级时被静默改写。任何数量变化都要有来源和解释，不能仅因为 `status=complete` 就接受。

QC 通常包含构建时间和本机绝对路径，不能当成字节可复现的运行时 trust anchor。
原始 QC 不应原地“清理”；如需发布，应另建 sanitized derivative、计算新 SHA 并记录
派生关系。真正进入生产配置的是 5 个已审核运行时文件。

## 10. 提升 candidate 并更新生产配置

QC 通过后，停止并发生产加载，把 candidate 作为一个不可变、版本化 release 保留；
不要先删除旧 release。推荐路径：

```text
data/processed/releases/<releaseTag>/
```

先把 candidate 整目录复制为新的 release；目标已经存在时必须停止，不能覆盖或合并：

```powershell
$releaseTag = "<SAME_RELEASE_TAG_USED_FOR_CANDIDATE>"
$candidate = "data\processed\candidates\$releaseTag"
$releaseRoot = "data\processed\releases"
$release = Join-Path $releaseRoot $releaseTag
if (-not (Test-Path -LiteralPath $candidate -PathType Container)) {
  throw "candidate does not exist: $candidate"
}
if (Test-Path -LiteralPath $release) {
  throw "release already exists; never overwrite or merge it: $release"
}
New-Item -ItemType Directory -Path $releaseRoot -Force | Out-Null
Copy-Item -LiteralPath $candidate -Destination $release -Recurse -ErrorAction Stop
```

用审核记录中预先确认的 SHA-256 核对最终目录内恰好要进入生产配置的 5 个文件。下面的
文件名和 digest 都必须替换，不能用复制后临时算出的值反过来当“预期值”：

```powershell
$expectedRuntime = [ordered]@{
  "<PPI_FILENAME>"            = "<APPROVED_PPI_SHA256>"
  "<EXPRESSION_FILENAME>"     = "<APPROVED_EXPRESSION_SHA256>"
  "<PATHWAYS_FILENAME>"       = "<APPROVED_PATHWAYS_SHA256>"
  "<DISEASE_FILENAME>"        = "<APPROVED_DISEASE_SHA256>"
  "<TARGET_MAPPING_FILENAME>" = "<APPROVED_TARGET_MAPPING_SHA256>"
}
if ($expectedRuntime.Count -ne 5) { throw "expected exactly five runtime files" }
foreach ($entry in $expectedRuntime.GetEnumerator()) {
  $path = Join-Path $release $entry.Key
  if (-not (Test-Path -LiteralPath $path -PathType Leaf)) {
    throw "missing promoted runtime file: $path"
  }
  $actual = (Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant()
  if ($actual -ne $entry.Value.ToLowerInvariant()) {
    throw "promoted runtime hash mismatch: $path"
  }
}
```

复制失败、文件缺失或任一 digest 不匹配时，停止提升并隔离这个不完整的新目录；不要更新
正式 config/locks，也不要覆盖旧 release。五项全部通过后才执行下面的切换准备。

切换前把当前配置与 lock 成对保存到被 Git 忽略的审计目录：

```powershell
$backup = "data\processed\release_manifests\$(Get-Date -Format yyyyMMdd_HHmmss)"
New-Item -ItemType Directory -Path $backup
Copy-Item -LiteralPath data\production_context.json -Destination $backup
Copy-Item -LiteralPath data\runtime.lock.json -Destination $backup
```

然后更新 `data/production_context.json` 的 5 个 resource：

- `path`：相对于 `data/`，例如 `processed/releases/<tag>/<file>`；
- `source` 和 `version`：真实上游版本；
- `taxid` 和 `identifier_namespace`；
- `mapping_version`：所有层使用同一 HGNC release；
- `sha256`：最终提升文件的 SHA-256。

这里有两套不同的相对路径基准，不能互相照抄：

| 文件 | 路径基准 | 示例 |
|---|---|---|
| `data/production_context.json` 的 resource `path` | 相对于 config 所在的 `data/` | `processed/releases/<tag>/x.gz` |
| `data/runtime.lock.json` 的 `relative_path` | 当前 lock 的 `root: ".."` 指向项目根 | `data/processed/releases/<tag>/x.gz` |

当前 runtime lock 位于 `data/`，生产 loader 会读取它的 `root`，因此必须保留
`"root": ".."` 或进行一次明确的 lock-layout 迁移。standalone verifier 则忽略
manifest 内的 `root`，直接以命令行 `--root` 解析 `relative_path`。

先把所有新文件放到最终路径，再切换配置。保留旧配置、旧 runtime lock 和旧 release，
这样回滚只需恢复上一组已验证的配置/lock，而不必重新构建 35 GiB 数据。

## 11. 按固定顺序更新三类 lock

lock 是人工审查后的信任锚，建议按下面顺序更新：

1. `source_registry.json`：当前构建所需 12 项的 active path、source 中的版本标识、大小、SHA、URL、政策。
2. `source.lock.json`：当前保留的完整 `data/raw` inventory；更新
   `source_artifact_count`、`source_size_bytes/GiB`、build registry 和来源摘要。
3. `manifest.lock.json`：详细来源、可用的上游 checksum、许可证、verification status、
   构建版本和相关 sidecar 元数据。
4. 定稿 5 个运行时文件及其最终路径。
5. 更新 `production_context.json` 的路径、版本、mapping version 和 SHA。
6. 最后计算定稿后的 `production_context.json` 自身大小/SHA，并更新
   `runtime.lock.json`。

`runtime.lock.json` 必须同时覆盖：

- `production_configuration`；
- PPI、expression、pathways、disease、target mapping 五个 runtime artifacts；
- 仓库显式跟踪的受控 PLATO fixture manifest 和 JSON；fresh clone 已包含这两个小型例外；
- 各分组字节数和 `total_size_bytes`。

其中 `total_size_bytes` 是 lock 中所有带 `relative_path` 和 digest 的记录大小之和。
`production_context.json` 必须先定稿再计算哈希，否则会形成自相矛盾的 lock。
当前分组汇总公式为：

```text
runtime_artifact_bytes = production_configuration.size_bytes
                       + 五个 runtime_artifacts 的 size_bytes
fixture_bytes          = offline_test_fixtures 的 size_bytes 总和
total_size_bytes       = runtime_artifact_bytes + fixture_bytes
```

`verify_manifest.py` 不检查这些汇总字段；现有测试只检查总和中的一部分一致性，所以
更新人必须独立重算并审阅三个值。

提升 live 前，建议组装一个与项目相同相对路径结构的验证根目录，放入候选
`data/production_context.json`、`data/runtime.lock.json`、5 个 runtime 文件和锁中
声明的 fixture，然后执行：

```powershell
python scripts\verify_manifest.py --root <CANDIDATE_ROOT> --scope runtime
```

这一步只验证 lock 声明的字节。随后必须实际加载候选配置，交叉检查 config 内部路径、
hash、taxid、namespace、HGNC mapping version、target mapper 和 runtime lock 语义：

```powershell
$candidateConfig = "<CANDIDATE_ROOT>\data\production_context.json"
python -c "import sys; from pathlib import Path; sys.path.insert(0, str(Path.cwd().parent)); from InterpreterPLATO import ExplanationPipeline; print(ExplanationPipeline().load_production_resources(sys.argv[1]))" $candidateConfig
```

standalone verifier 使用命令行的 `--root` 解析相对路径，不读取 manifest 的 `root`；
production loader 会读取 runtime lock 的 `root`。两项都通过，才表示候选同时满足字节
身份和加载语义。

## 12. 最终验证和发布门槛

先做快速生产校验，再做完整源校验：

```powershell
python scripts\verify_manifest.py --scope runtime
python scripts\bootstrap_data.py verify
python scripts\verify_manifest.py --scope source
python scripts\verify_manifest.py --manifest data\manifest.lock.json
python -m pytest
```

还应实际加载 `production_context.json`，对代表性细胞系运行 smoke test。哈希全部通过只
代表文件身份一致，不代表 PPI 覆盖、表达细胞系覆盖、通路统计或疾病模块在生物学上合理。

提交前确认数据边界：

```powershell
git status --short --ignored
git check-ignore -v data\raw\<file>
git check-ignore -v data\processed\<file>
git check-ignore -v outputs\<file>
```

允许提交的是代码、文档、配置、registry 和 lock 元数据；不要提交普通 raw/processed
字节、普通查询的 PLATO 响应、PDF、缓存、模型权重或运行输出。当前唯一的数据例外是
Git 已显式跟踪且 locks 已固定的 fixture manifest 和 target-fishing JSON。新增或替换
fixture 仍须通过许可、脱敏和独立发布审查，并同步更新忽略规则、测试及 locks。

## 13. 更新普通药物的 PLATO 数据

一个新药物响应不是癌症 context 更新，通常不修改 `runtime.lock.json` 或
`source.lock.json`。推荐使用单次查询 helper，因为它能绑定本次提交的精确 SMILES：

```powershell
$run = "D:\secure\plato-runs\20260719-study001"

python scripts\query_plato_once.py `
  --drug-name "Study-Drug-1" `
  --smiles "<EXACT_SMILES_1>" `
  --release r37 `
  --raw-dir "$run\raw\drug1" `
  --cache-dir "$run\cache"

Start-Sleep -Seconds 3

python scripts\query_plato_once.py `
  --drug-name "Study-Drug-2" `
  --smiles "<EXACT_SMILES_2>" `
  --release r37 `
  --raw-dir "$run\raw\drug2" `
  --cache-dir "$run\cache"
```

它只提交一次 HTTPS target-fishing 请求，并保存：

- 原始 JSON 字节；
- metadata sidecar；
- response SHA-256 和大小；
- endpoint、HTTP method、提交表单和 query binding；
- schema v3 PLATO cache。

原始响应不可编辑。相同 cache key 已存在而内容变化时，代码会拒绝静默覆盖；新 release、
新 parser 或新实验应使用独立 cache/run 目录并保留旧档。精确 SMILES 会发送给外部
PLATO 服务；专有结构必须先获得授权。命令行参数可能进入进程列表和 PowerShell 历史，
`--print-report` 还可能把 SMILES 和绝对路径写入终端日志。

手工下载只能用于证据归档：

```powershell
python scripts\import_plato_result.py `
  --drug-name "Drug-A" `
  --smiles "<当时提交的精确 SMILES>" `
  --result "D:\plato\drug-a.json" `
  --release r37 `
  --attest-downloaded-from-plato `
  --print-report
```

`--attest-downloaded-from-plato` 只是用户声明，不是平台认证。当前观察到的 r37 JSON
本身没有 query SMILES，手工导入脚本也不接受可审计的 transport sidecar，所以这种导入
不能满足正式 query-binding gate。

同一次查询刚算出的 SHA-256 只能证明保存后的字节没变，不能自动成为可信 allowlist。
独立审核记录不应只保存一列裸哈希；至少绑定 study/compound ID、canonical InChIKey、
PLATO release、原始字节数、审核人、审核时间和代码 commit。正式计分再显式传入两个
已批准的原始响应哈希：

`--trusted-plato-sha256` 只接受原始响应字节的 SHA-256，即 cache
`source.sha256`/`raw_blobs/<sha>.json` 的身份。`payload_sha256`、normalized-target
hash、cache envelope hash、schema fingerprint 或 target-mapping hash 都不能替代它。
allowlist 应按实验最小化，由独立审核人重新计算并批准，不应无限扩张为全局清单。

```powershell
python scripts\run_interpreter.py summary `
  --context-config data\production_context.json `
  --plato-cache "$run\cache" `
  --trusted-plato-sha256 <DRUG_A_RAW_SHA256> `
  --trusted-plato-sha256 <DRUG_B_RAW_SHA256> `
  --drug1 "Drug-A" --drug2 "Drug-B" `
  --drug1-smiles "<SMILES-A>" --drug2-smiles "<SMILES-B>" `
  --cell-line A549 `
  --disease-label pan_cancer_tier1 `
  --output "$run\explanation.json"
```

pair helper 如果没有获得匹配的外部 allowlist，会保留两次 raw/cache 和
`run_manifest.json`，随后以 exit code 2 停止，不生成正式 explanation。若第二次查询在
manifest 写出前失败，第一次查询的 raw/cache 仍会保留，应按 run 目录检查残留证据。
`--quick` 只适合调试，不能替代正式稳定性和敏感性分析。

PLATO 原始内容至少可能存在于 raw 文件、cache 的 `raw_blobs/` 和 cache envelope 的
`raw_result` 中；只清理一个目录并不等于完成脱敏。原始审计文件不能原地修改，应另建
`sanitized_derivative`，计算新 SHA，并明确它不再具有原始 allowlist 身份。需要保护的
内容可能包括专有 SMILES/InChIKey、药物名称、supporting ligands、绝对路径、研究 ID、
细胞系和疾病标签。普通查询数据不提交 Git，也不自动加入 source/runtime lock。

## 14. 更新受控 PLATO schema fixture

只有在受控 r37 fixture 新增、替换，或 PLATO schema 确实改变时，才执行本节。普通药物
查询不应修改 fingerprint allowlist 或 context locks。

当前 manifest 仍明确记录 `status=raw_json_and_pdf_fixture_candidate` 和
`production_use=not_yet_allowlisted_for_scoring`，production schema allowlist 因此保持为空。
当前 candidate 仍可用于 parser 回归，但必须得到 `schema_fixture_verified=false`，严格解释器
会 fail closed。不能把“可解析”写成“来源已审计”或“生产已放行”。

manifest notes 中现有的 “raw-file trust anchors” 只能按字节完整性标识理解，不能当成独立
来源认证、production schema approval 或每药物 raw-SHA allowlist 条目。

只有独立审计完成并留下可核查的审核状态/证据后，才能在同一个 reviewed change 中同步
更新 fixture approval、provenance、`VERIFIED_PLATO_SCHEMA_FINGERPRINTS` 和 contract tests。
不得只改 manifest，代码也不会因 payload 自称 r37 而自动批准 candidate。

1. 保存真实原始响应和独立获取记录；原始字节永不脱敏覆盖。
2. 另做脱敏副本。只能改变敏感 scalar value，不能随意删 key、改嵌套或改数据类型，
   否则 schema fingerprint 会改变。
3. 在受限 provenance manifest 中记录 provider、release、method、query compound、精确
   SMILES、捕获时间、HTTP 证据、文件大小和 SHA-256；敏感凭据不得进入 manifest。公开的
   脱敏 manifest 只能使用化名，且不得暴露专有 SMILES；它与受限原件是两个独立产物。
4. 用解析器检查字段覆盖、记录数、warnings 和 fingerprint：

```powershell
$fixture = "<FULL_PATH_TO_CONTROLLED_FIXTURE_JSON>"
python -c "import json,sys; from pathlib import Path; from InterpreterPLATO.targets import inspect_plato_payload; p=json.loads(Path(sys.argv[1]).read_text(encoding='utf-8-sig')); r=inspect_plato_payload(p); print(r['schema_fingerprint']); print(r['warnings'])" $fixture
```

对当前 candidate，预期 `schema_fingerprint_allowlisted=false`、
`schema_fixture_verified=false`，并出现 schema 未验证 warning；如果它在独立审核前返回
`true`，必须停止发布。

5. 只有独立审计确认它确实代表目标 release 的 target-fishing schema 后，才能更新
   `VERIFIED_PLATO_SCHEMA_FINGERPRINTS`。不能因为某个文件被拒绝就把它的 fingerprint
   加进去。
6. 同步更新 fixture manifest、contract tests、`source.lock.json`、
   `runtime.lock.json` 和详细 provenance manifest。
7. 重新跑完整 pytest 和 production preflight。

schema fingerprint 只覆盖 JSON key path 和 value kind，不证明每份未来文件来自 PLATO；
每个药物仍必须通过 query binding、独立 raw SHA allowlist、人源蛋白靶点和 context
完整性门禁。

## 15. 失败与回滚

任何步骤失败时：

1. 停止提升，不要修改 hash 来匹配失败产物。
2. 保留 candidate、QC 和日志用于分析，不覆盖旧 release。
3. 恢复上一组 `production_context.json`、`runtime.lock.json` 和对应 runtime release。
4. 重新运行 runtime 校验，确认生产恢复。
5. 如果 source acquisition 错误，修复或重新获取 raw；不要让错误文件进入 source lock。
6. 如果 PLATO 第二次查询失败，保留第一次 raw/cache，标记本次 pair run 不完整。

生产 loader 会在内存中原子替换完整 context，但 builder 和整套磁盘 bundle 不是事务。
因此，candidate 构建、版本化 release 和保留上一版，是实际可用的回滚机制。

## 16. 更新完成检查表

- [ ] 新 raw 位于独立版本目录，旧 raw 未被覆盖。
- [ ] 可用的上游 checksum、许可证和 attribution 已保存并核验；无上游 checksum 的来源已明确证据等级和限制。
- [ ] release 常量、HGNC 硬哈希、默认路径和输出名已同步。
- [ ] source registry 的 active inputs 全部通过。
- [ ] candidate 中 context 和 target mapping 均构建成功。
- [ ] 全部 QC 和 rejection sidecar 已与上一版人工比较。
- [ ] ScoreRI 细胞系数量、别名和拒绝项已审核。
- [ ] 五个 production resources 使用同一 HGNC mapping version。
- [ ] production config 先定稿，runtime lock 后计算。
- [ ] runtime/source/provenance 三类 lock 与汇总字段一致。
- [ ] runtime、source、显式 provenance manifest 和 pytest 全部通过。
- [ ] 代表性 production load/smoke test 通过。
- [ ] 未完成 fixture 独立审核时，production schema allowlist 为空且严格评分 fail closed。
- [ ] Git 中没有 raw、processed、普通 PLATO 响应、模型或输出；受控 fixture 例外已单独审批。
- [ ] 回滚所需的上一版配置、lock 和 runtime release 仍可用。
