# InterpreterPLATO 数据来源与获取说明书

本文说明 InterpreterPLATO 的数据来自哪里、当前锁定了哪些版本、下载后需要保留
哪些原件、不同 checksum 能证明什么，以及哪些文件只是本地归档、并未参与生产构建。
版本升级、candidate 构建、QC、提升和回滚的具体操作见
[数据更新操作手册](DATA_UPDATE_GUIDE.md)；只复现当前 pin 时见
[Data bootstrap](DATA_BOOTSTRAP.md)。

本文记录的是 2026-07-22 的当前基线。机器可读的文件路径、大小和 SHA-256 以
`data/source_registry.json`、`data/source.lock.json`、
`data/manifest.lock.json` 和 `data/runtime.lock.json` 为准；手册与 lock 不一致时必须
停止构建并审核差异，不能直接把观察到的新哈希抄回 lock。

> 许可提示：本文只记录项目的获取和审计边界，不替代任何上游许可。下载、使用或
> 再分发前，应由数据使用者检查对应 release 的官方条款、引用要求和机构授权。

## 1. 先区分四类数据

| 类别 | 当前范围 | 用途 | 是否随代码仓库发布 |
|---|---:|---|---|
| 最小重建输入 | `source_registry.json` 的 12 项，共约 35.04 GiB | 重建五个生产上下文文件 | 否；仅发布路径、大小、哈希和来源元数据 |
| 完整源审计包 | `source.lock.json` 的 40 个 raw 文件，共约 35.45 GiB | 原件、容器、许可文件及历史/备用数据审计 | 否 |
| 快速运行时包 | 生产配置、五个上下文文件和两个 PLATO candidate fixture 文件，共 7.64 MiB | 加载上下文和离线 parser 回归 | code-only 发布不含这些字节 |
| 每药物 PLATO 原始响应 | 每次单分子查询的 raw JSON、sidecar 和 SHA-256 | 唯一正式候选靶点来源 | 否；按查询单独归档和审批 |

“在 `data/raw` 中”不等于“已接入生产”。当前真正被 builder 消费的只有 registry
中的 12 项；OncoKB、GO、GDSC、PRISM、MSigDB 等虽然在完整 source lock 中有字节
记录，但当前不参与 PPI、表达、通路、疾病集合或 target-ID mapping 的构建。

当前也没有在仓库中公布 runtime/source bundle 的下载 URL、制品签名或数据管理员
申请入口。因此，fresh code-only clone 可以安装代码、查看获取计划和运行不依赖外部
数据的测试，但不能仅凭仓库恢复当前完整数据包，也不能通过
`verify_manifest.py --scope source`。这不是安装错误。若项目方以后分发数据包，至少
必须同时公布：制品名称和版本、适用代码 commit、包自身 SHA-256 或签名、访问方式、
解包根目录、上游许可清单和撤回办法。

## 2. 数据如何进入生产上下文

| 生产层 | 上游输入 | 当前运行时产物 | 作用边界 |
|---|---|---|---|
| 基因标识 | HGNC 2026-04-01 | `gene_map.tsv` 等审核产物 | 给其余五路提供统一 approved-symbol 命名空间 |
| PPI | HGNC + STRING v12.0 human info/full links | `string_v12.0_hgnc_2026-04-01.tsv.gz` | 网络邻近和机制上下文 |
| 表达 | HGNC + DepMap expression/Model/Profiles + ScoreRI cell-line registry | `depmap_26Q1_scoreri_log2tpm.csv.gz` | ScoreRI 细胞系表达上下文 |
| 通路 | HGNC + Reactome v96 GMT | `reactome_v96_hgnc_2026-04-01.gmt` | 通路富集和共同机制 |
| 疾病集合 | HGNC + COSMIC CGC v101 | `cosmic_cgc_v101_pan_cancer_tier1.txt` | pan-cancer Tier 1 疾病模块 |
| 靶点 ID 翻译 | HGNC + ChEMBL r37 target/component 表 | `chembl_r37_target_hgnc_mapping.tsv.gz` | 只把 PLATO 的 ChEMBL target ID 翻译为 HGNC symbol |
| 药物候选靶点 | 每个药物各自的 PLATO 单次查询响应 | `data/processed/plato_cache/` 中的审核缓存 | 唯一正式候选靶点来源；不由上述本地数据补全 |

所有基因上下文最终统一到同一个 HGNC release。approved symbol 中本来存在的连字符
是标识的一部分，不能用“删掉所有连字符”的方式对齐；builder 应通过 approved
symbol、唯一 Entrez、唯一 alias/previous symbol 等受审计映射处理来源差异，歧义项
必须拒绝。

ChEMBL 在这里不是第二套药物靶点数据库。代码只执行：

```text
PLATO target_chembl_id
  -> ChEMBL r37 target/component
  -> HGNC 2026-04-01 approved symbol
```

不得把 ChEMBL drug-target association 导入为额外靶点证据。

## 3. 哪个文件回答哪个审计问题

| 文件 | 回答的问题 | 不证明什么 |
|---|---|---|
| `data/source_registry.json` | 当前重建所需 12 项输入在哪里、应有多大、SHA-256 是什么 | 不代表完整 raw 目录，也不自动认证来源 |
| `data/source.lock.json` | 完整 40 文件 source snapshot 的字节身份是什么 | 不代表 40 项都参与生产，也不是完整许可判定 |
| `data/manifest.lock.json` | 主要生产来源的 release、官方入口、许可提示、容器和派生关系是什么 | 当前并未逐项覆盖所有 40 个 raw 文件的完整 provenance |
| `data/production_context.json` | 生产实际加载哪五个上下文文件 | 不记录完整重建源 |
| `data/runtime.lock.json` | 配置、五个运行时文件和离线 candidate fixture 的字节是否一致 | 不批准 PLATO fixture，也不证明生物学机制成立 |

快速查看当前 12 项来源：

```powershell
python scripts\bootstrap_data.py plan
```

这个命令只读 registry，不联网、不下载、不解压，也不生成 lock。

## 4. 当前 12 项最小重建输入

下表中的“证据等级”定义见第 6 节。

| Registry ID | 锁定来源 | 当前注册文件 | 证据等级 |
|---|---|---|---|
| `hgnc` | HGNC Complete Set，2026-04-01 quarterly archive | `hgnc_complete_set_2026-04-01.txt` | B：官方日期化 URL + 本地 SHA |
| `string_info` | STRING v12.0 human protein info | 解压后的 `9606.protein.info.v12.0.txt` | B：官方版本列表 + 本地 SHA |
| `string_links` | STRING v12.0 human full links | 解压后的 `9606.protein.links.full.v12.0.txt` | B：官方版本列表 + 本地 SHA |
| `depmap_expression` | DepMap Public 26Q1 expression | `OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv` | B：release 文件页 + 本地 SHA |
| `depmap_model` | DepMap Public 26Q1 Model metadata | `data/raw/cell_line_context/depmap/26Q1/Model.csv` | B：release 文件页 + 本地 SHA |
| `depmap_profiles` | DepMap Public 26Q1 OmicsProfiles metadata | `data/raw/cell_line_context/depmap/26Q1/OmicsProfiles.csv` | B：release 文件页 + 本地 SHA |
| `scoreri_cell_lines` | ScoreRI 历史组合表的唯一细胞系名 | `data/scoreri_cell_lines.tsv` | D：本地派生 registry |
| `reactome` | Reactome v96 pathways | ZIP 唯一成员 `ReactomePathways.gmt` | A2：官方版本 URL 容器/成员复核 |
| `cosmic` | COSMIC CGC v101 GRCh38 | tar -> `.tsv.gz` -> `.tsv` | A1：官方 archive SHA 已匹配 |
| `chembl_archive` | ChEMBL r37 SQLite archive | `chembl_37_sqlite.tar.gz` | A1：官方 checksum 已匹配 |
| `chembl_checksums` | ChEMBL r37 checksum 表 | `checksums.txt` | B：官方 release 元数据 + 本地 SHA |
| `chembl_sqlite` | ChEMBL r37 解压 SQLite | `chembl_37.db` | C：本地字节已锁；当前记录未独立证明它是该 archive 的成员 |

当前 registry 的精确大小和 SHA-256 如下。升级时不要手改本段代替 registry 审核：

```text
hgnc               17030334  ad19966c08ccebd7b0fad25d553a5db8770a188fef9cf93fc24dd74796c45d34
string_info         6257901  a65e52748022d6a8e953843645a5032b094017225d95de9c6cb6fd4577d0a7ad
string_links     1043993447  d9da3020817fff2d7be6f960848be49e780b22ce04e73bcf63c4874ca6ddfde3
depmap_expression 305007605  0377be80c525fde98cbd2c6e8b06bdf2a4014a9683eb70182c1f8649d711021a
depmap_model          697455  ea4e0b2a3bc806f81df62689a5ae75f1a100135727a3d7b8a4c7ccc8815183f8
depmap_profiles       671558  978f9ff9a11214e9d923619c20d2bc36bf4ca3e30fd46a6d9a424ea0e7c49bd2
scoreri_cell_lines       569  139da95d7150be6af73baa217a2e627d21767240cd2f59d9db9776a1b7742db4
reactome             1017374  cf17607f842c1742d3ecb96e6a4b03f11889c39f44005a37c31d12e4c006b65e
cosmic                152821  c7cc8514db46411bc88e01efde97070cb229f19175ac76207a8cca69b1f108b0
chembl_archive    5764252857  33c203740555f96067710cdfc1c3c55d890660e5908ec5cbf5817492c290d281
chembl_checksums          790  e949808982cc9f0d08fa16b42ebe75b4e04b511d8b96693c973890e6e024e62a
chembl_sqlite     30480314368  4be13df3b68e25dcd0bff44bf094033b5aebe98f415acdc8c1cdf380e0c15142
```

## 5. 逐源获取和落盘要求

### 5.1 HGNC：共同基因命名基线

- 官方归档入口：[HGNC archived data](https://www.genenames.org/download/archive/)
- 当前日期化原件：[hgnc_complete_set_2026-04-01.txt](https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/hgnc_complete_set_2026-04-01.txt)
- 当前落盘：
  `data/raw/gene_mapping/hgnc/2026-04-01/hgnc_complete_set_2026-04-01.txt`
- 访问与许可：公开入口；当前官方
  [HGNC licence](https://www.genenames.org/about/license/) 将数据置于 CC0，仍建议保留
  HGNC 来源和版本归属；再分发时以所用 release 的条款为准。
- 当前验证：官方日期化 URL 和文件名固定 release；没有匹配到独立发布的上游
  checksum，本地 SHA 只能固定取得后的字节。

HGNC 被 PPI、表达、Reactome、COSMIC 和 ChEMBL mapping 五路共同复用。升级 HGNC
必须整批重建这五路，不能让一个运行时包混用两个 HGNC release。

### 5.2 STRING：human PPI

- 官方下载页：[STRING downloads](https://string-db.org/cgi/download.pl)
- 上游压缩文件：
  `9606.protein.info.v12.0.txt.gz`、
  `9606.protein.links.full.v12.0.txt.gz`
- 当前 registry 锁定的是 gunzip 后的两个 `.txt`，分别放在
  `data/raw/ppi/string/`。
- 当前验证：官方 human/taxon 版本列表确认文件名；没有匹配到独立上游 checksum，
  对解压后字节记录本地 SHA。
- 当前官方 [STRING licensing](https://string-db.org/cgi/access?footer_active_subpage=licensing)
  标注 CC BY 4.0；使用和再分发时同时保留 release、引用和归属。

当前锁定布局要求先得到与 registry 完全相同的解压文本。builder 本身能够读取
`.gz`，但仅换成 `.gz` 路径不会自动改变 registry 的 expected size/SHA，也不等于
复现当前 pin。容器、解压工具/版本、解压后大小和 SHA 都应写入获取记录。

### 5.3 DepMap：表达矩阵和模型元数据

- 当前 expression 文件页：[DepMap Public 26Q1 expression](https://depmap.org/portal/data_page/?file=OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv&release=DepMap+Public+26Q1)
- 当前 Model 文件页：[DepMap Public 26Q1 Model](https://depmap.org/portal/data_page/?file=Model.csv&release=DepMap+Public+26Q1)
- 当前 OmicsProfiles 文件页：[DepMap Public 26Q1 OmicsProfiles](https://depmap.org/portal/data_page/?file=OmicsProfiles.csv&release=DepMap+Public+26Q1)
- 当前 release 入口：[DepMap Public 26Q1 all data](https://depmap.org/portal/data_page/?release=DepMap+Public+26Q1&tab=allData)
- 当前落盘：
  `data/raw/cell_line_context/depmap/26Q1/OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv`、
  `data/raw/cell_line_context/depmap/26Q1/Model.csv`、
  `data/raw/cell_line_context/depmap/26Q1/OmicsProfiles.csv`
- 数据语义：源值为 `log2(TPM + 1)` 的 human protein-coding expression；builder
  不应再次取对数。
- 访问、条款和引用以所选 release 页面当时的提示为准。

`Model.csv` 和 `OmicsProfiles.csv` 已于 2026-07-23 从 DepMap Public 26Q1 release
页面重新取得；本地 size/SHA 与先前锁定快照完全一致，因此当前表达层可以把三份
DepMap 输入作为同一 26Q1 批次处理。没有匹配到独立上游 checksum，证据等级仍为
B：release 页面确认文件身份，本地 SHA 固定取得后的字节。

新升级仍必须从同一个 DepMap release 同批获取 expression、`Model.csv` 和
`OmicsProfiles.csv`，保存 release 页面、下载时间、文件元数据和条款快照。新文件统一
放到 `data/raw/cell_line_context/depmap/<release>/`，三份文件交叉验证 default RNA
profile 后再构建；不能仅凭同名文件或本地 SHA 宣称 release 一致。

当前表达基线覆盖 ScoreRI registry 的 70 个 verified-default 细胞系；10 个已映射模型
没有 verified default expression，`SR` 仍明确未解析。这个缺口必须保留在 QC 中，不能
用近似字符串猜测模型。

### 5.4 ScoreRI：细胞系构建 registry

- 代码来源：[skeujann/ScoredRI](https://github.com/skeujann/ScoredRI)
- 当前 `data/scoreri_cell_lines.tsv` 不是表达数据，也不是完整组合表；它只保留历史
  `ScoreRI/data/combo_data.tsv` 的 `cell_line_name` 列排序去重后的 81 个标识。
- 原组合表记录的 SHA-256：
  `52367a1a8fcdd12db9ca9bfc1285820d4295c6d4212d1a7317dfa737c7418c4d`。
- 原组合表不随 InterpreterPLATO 代码发布；只有经审核的 cell-line registry 参与
  expression 构建。
- ScoreRI 仓库的代码许可不能自动证明 `combo_data.tsv` 的数据许可；发布 registry 或
  重新取得原表前需要单独确认数据所有者和使用条款。

当前仓库没有记录该原表可公开下载的 artifact URL、tag/commit 或数据管理员申请入口。
如果需要重建 81 个名字，应由原数据保管者提供原表和上述哈希，并把实际 ScoreRI code
commit、文件相对路径、取得时间和许可一起记入新 provenance。

当前 `depmap_expression_qc.json`/`build_qc.json` 已由 2026-07-23 的完整 `all` 重建
刷新：builder 读取 569-byte `data/scoreri_cell_lines.tsv`，不再依赖外部
`ScoreRI/data/combo_data.tsv`。表达矩阵主产物 hash 与旧基线一致；QC 中的输入路径和
输入哈希反映当前可独立复现的 registry。

### 5.5 Reactome：通路集合

- 官方数据入口：[Reactome download data](https://reactome.org/download-data/)
- 当前版本容器：[Reactome v96 GMT ZIP](https://download.reactome.org/96/ReactomePathways.gmt.zip)
- ZIP 当前应只有本项目使用的成员 `ReactomePathways.gmt`；registry 锁定的是提取后
  GMT，落在 `data/raw/pathways/reactome/ReactomePathways.gmt`。
- 当前验证：从官方 `/96/` URL 取得容器，记录容器本地 SHA，并确认提取成员的 SHA
  与本地 GMT 一致；这不是“上游另行公布 checksum”。
- 当前官方 [Reactome licence](https://reactome.org/license) 将 annotation/data files
  置于 CC0，并鼓励归属；仍应记录 release 和引用。

当前无版本目录是 legacy 布局。下一次升级应放到
`data/raw/pathways/reactome/<release>/`，不要覆盖 v96 字节。
当前 source lock 只保留提取后的 GMT；ZIP 的大小、SHA 和 member-match 观察记录在
`manifest.lock.json`，ZIP 字节本身未进入 40 文件 raw snapshot。

### 5.6 COSMIC：Cancer Gene Census

- 官方 release 入口：[COSMIC CGC v101](https://cancer.sanger.ac.uk/cosmic/download/cosmic/v101/cancergenecensus)
- 官方条款：[COSMIC terms](https://www.cosmickb.org/terms/)
- 当前 release：v101、GRCh38。
- 原始容器：`Cosmic_CancerGeneCensus_Tsv_v101_GRCh38.tar`。
- 容器内保留 README 和 `Cosmic_CancerGeneCensus_v101_GRCh38.tsv.gz`；registry 消费
  再解压后的 `.tsv`。
- 当前官方 archive SHA-256 与本地 tar 一致：
  `2e9565c45e89985fb1c096f29af7cc67166ffaf4eb8964f004986f653e22c755`。

COSMIC 是账号/许可受限来源。必须由获授权用户在授权环境中下载；raw tar、gz、TSV、
账号、Cookie、访问令牌和临时签名 URL 都不能随代码提交或写入公开 manifest。可以
记录稳定 landing page、release、下载时间、官方 checksum、操作者角色和本地哈希，
但不要记录凭据。

### 5.7 ChEMBL：PLATO target ID 翻译

- 官方 release 目录：[ChEMBL 37 files](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_37/)
- 同目录应取得：`chembl_37_sqlite.tar.gz`、`checksums.txt`、`LICENSE`、
  `REQUIRED.ATTRIBUTION` 和 `schema_documentation.txt`。
- 保留原 archive，同时解压出
  `chembl_37/chembl_37_sqlite/chembl_37.db`；builder 不负责解压。
- 正式构建不得使用 `--skip-archive-check`。
- r37 发布记录日期为 2026-05-29；精确许可和 attribution 以本 release 随包的
  `LICENSE` 和 `REQUIRED.ATTRIBUTION` 为准。

ChEMBL builder 会用 `checksums.txt` 验证 archive，但不能认证 `checksums.txt` 本身，
也不会仅凭查询 SQLite 证明 DB 确实来自该 archive。因此需要三层记录：官方 release
目录取得的 checksum 文件、本地 archive 匹配结果、解压 SQLite 的独立 registry/lock
哈希。当前 manifest 把 DB 标为 extracted SQLite，但没有独立 member-binding 记录，
所以这仍是已知 provenance 缺口。应从已验证 archive 向空的版本目录执行一次干净解压
并记录工具/命令和成员关系；不能把另行取得的同名 `chembl_37.db` 与已验证 archive
拼在一起。

只读核对 archive 的示例：

```powershell
$checksums = Import-Csv -LiteralPath `
  data\raw\target_mapping\chembl\r37\checksums.txt -Delimiter "`t"
$rows = @($checksums | Where-Object { $_.file -eq "chembl_37_sqlite.tar.gz" })
if ($rows.Count -ne 1) { throw "expected exactly one ChEMBL archive checksum" }
$expected = ([string]$rows[0].SHA256).Trim().ToLowerInvariant()
$actual = (Get-FileHash -LiteralPath `
  data\raw\target_mapping\chembl\r37\chembl_37_sqlite.tar.gz `
  -Algorithm SHA256).Hash.ToLowerInvariant()
if ($actual -ne $expected) { throw "ChEMBL archive checksum mismatch" }
```

### 5.8 PLATO：每个药物的原始靶点响应

- helper 当前允许的默认表单 endpoint：
  `https://prometheus.farmacia.uniba.it/plato/pred.html`。
- 一次运行只允许一个用户触发的单分子 target-fishing 请求；不循环、不并发、不批量
  爬取。
- 每个药物分别保存 raw response bytes、query sidecar、提交的精确 SMILES、HTTP/表单
  元数据、captured time 和 SHA-256。
- 用户声明“这是从 PLATO 下载的”只是一条 attestation，不是平台签名或来源认证。

仓库当前没有保存可独立审核的 PLATO 官方 release 说明、服务条款快照或 r37 schema
规范。现有 JSON 也没有 release 或 query-SMILES 字段，脚本的 `--release` 是调用方
提供值；因此目录名、manifest 中的 `reported_release=r37` 和调用方参数都不能证明
provider release。Gemcitabine JSON/PDF 只是 candidate，production schema allowlist
为空，不能打开正式计分。未来必须先独立审核一个真实、脱敏、可追溯的 r37 导出，固定
schema fingerprint；随后每个药物的 raw export SHA-256 还要分别进入独立可信
allowlist。runtime lock 只能证明 candidate 字节未变，不能代替这些信任门槛。

## 6. Checksum 和来源证据等级

| 等级 | 含义 | 当前例子 | 能证明什么 |
|---|---|---|---|
| A1 | 上游明确公布 checksum，且本地容器匹配 | COSMIC archive、ChEMBL archive | 所得容器与上游公布身份一致 |
| A2 | 从官方版本 URL 观察并锁定容器、成员和解压产物关系 | Reactome ZIP/member | 本次官方版本路径和转换链可追踪，但不应说成另有官方 checksum |
| B | 官方日期化/版本化 URL 或 release 文件页 + 本地 SHA | HGNC、STRING、DepMap expression | 固定本次取得字节和版本页面；本地 SHA 本身不是来源认证 |
| C | 只有本地 size/SHA，或容器到提取产物的绑定尚未独立记录 | DepMap legacy Model/Profiles、当前 ChEMBL SQLite | 只能识别受控旧副本，不能单独证明 release 或容器成员关系 |
| D | 从另一受控输入确定性派生 | ScoreRI cell-line registry | 能证明派生字节；还需保留上游输入哈希和转换规则 |

无论等级如何，都应同时记录最终 builder 输入的大小和 SHA-256。压缩或归档来源还应
分别记录：容器哈希、成员名称和成员哈希、解压后输入哈希、转换是否仅为无损解压。

本地计算：

```powershell
$path = "data\raw\<source>\<release>\<file>"
(Get-Item -LiteralPath $path).Length
(Get-FileHash -LiteralPath $path -Algorithm SHA256).Hash.ToLowerInvariant()
```

没有上游 checksum 时，获取记录至少补充稳定 URL、release、取得时间、原始文件名、
文件大小、可用的 HTTP `Last-Modified`/ETag、解压步骤和本地 SHA，并明确标记为 B 或 C，
不要写成“official checksum verified”。

## 7. 完整 source lock 中其余文件

完整 source lock 的 40 项比 registry 多 29 个 raw 文件，分成三类。
`manifest.lock.json` 当前只覆盖其中 18 个 raw 路径；另有 22 个路径没有逐文件
provenance 条目，扣除项目内部 README 和 ChEMBL 安装说明后仍有 18 个上游或伴随文件
缺记录。因此 `provenance_companion=manifest.lock.json` 只能表示关联文件，不能理解为
“40 项来源和许可证已经完整审计”。

### 7.1 生产来源的审计配套文件

- COSMIC 的 tar、README 和 `.tsv.gz`；
- ChEMBL 的 LICENSE、required attribution、schema documentation 和 SQLite 安装说明；
- DepMap `ModelCondition.csv`；
- raw 目录内的说明文件。

这些文件可能对 provenance 或人工审核重要，但不是当前 12 项 registry 中的直接 builder
输入。不要为了让 registry 数量一致而删除它们。

### 7.2 尚未接入生产的本地数据

| 本地目录/文件 | 官方入口 | 当前状态 |
|---|---|---|
| `cancer_driver_context/oncokb/cancerGeneList.tsv` | [OncoKB Cancer Genes](https://www.oncokb.org/cancer-genes)、[terms](https://www.oncokb.org/terms) | release/date 未知；文件内的 “COSMIC CGC (v99)” 只是字段，不是 OncoKB release；未参与 disease layer |
| `functional_annotations/go/goa_human.gaf.gz` | [GO annotation downloads](https://geneontology.org/docs/download-go-annotations/)、[citation/licence](https://geneontology.org/docs/go-citation-policy/) | GAF header 记录 GO 2026-03-01、UniProt GOA 2026-03-04；legacy 名称已变化，未来应登记不可变月度 release，不能用 current URL 覆盖 |
| `monotherapy_response/gdsc/*27Oct23.xlsx` | [GDSC bulk downloads](https://www.cancerrxgene.org/downloads/bulk_download) | 文件名和官方路径对应 release 8.5 / 27Oct23，但当前 manifest 未留下载/checksum 记录；未被 builder 消费 |
| `monotherapy_response/prism/*` | [PRISM Repurposing 19Q4](https://depmap.org/portal/data_page/?release=PRISM+Repurposing+19Q4&tab=allData)、[Figshare archive](https://figshare.com/articles/dataset/PRISM_Repurposing_19Q3_Primary_Screen/9393293) | 内容对应 19Q4 / Figshare v4；应补 DOI、version、file IDs 和上游 checksum，未被 builder 消费 |
| `pathways/hallmark/h.all.v2026.1.Hs.symbols.gmt` | [MSigDB](https://www.gsea-msigdb.org/gsea/msigdb/index.jsp) | v2026.1.Hs，source-lock 归档，未替代 Reactome；下载需按官方注册和条款执行 |
| `pathways/msigdb/c5.go.bp.v2026.1.Hs.symbols.gmt` | [MSigDB license terms](https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp) | v2026.1.Hs，source-lock 归档，未参与当前 pathway layer；子集合条款和引用需逐项保留 |

对于这些文件，`source.lock` 只证明“本地当时存在这些字节”。在补齐 release、稳定官方
URL、许可、取得记录和独立 checksum 证据之前，不得把它们静默接入生产，也不要用
官网的 `latest` 覆盖旧文件。若要启用，必须新增 registry/manifest 条目、明确 builder
语义、HGNC 统一规则、QC、测试和 runtime 配置，并经过一次完整 review。

当前 registry schema 没有独立 `release`/`version` 字段，版本只编码在 `source` 文本
和部分路径中。升级时必须让两者与 builder pin 保持一致；后续若扩展 schema，应新增并
校验结构化 release 字段，不能在文档里假定该字段已经存在。

### 7.3 PLATO candidate fixture

完整本地审计包还包含 Gemcitabine target-fishing JSON、两份 PDF 和 fixture manifest。
code-only 发布不包含这些字节；普通药物查询响应也不得仿照该目录直接提交。candidate
只用于 parser contract 和人工审核，不是 production-approved 数据，也不是任意
PLATO-like JSON 的通行证。

## 8. 新数据的目录和 provenance 记录

新 release 一律使用不可变版本目录：

```text
data/raw/<用途>/<来源>/<release>/
  README 或 acquisition record
  上游 LICENSE / terms snapshot / required attribution
  上游 checksum 文件（如有）
  原始 archive 或压缩文件
  builder 最终消费的无损解压产物
```

Reactome 的无版本目录是 legacy，不应成为新版本模板；DepMap metadata 已改为随
expression 放入 release 目录，后续也应保持这个布局。不要使用 `latest/`，不要覆盖旧
release，也不要把账号、Cookie、token 或签名 URL 写进记录。

每次获取至少记录：

```text
source / dataset / release / organism / genome build（如适用）
official landing page / exact stable source URL
retrieved_at_utc / operator role / authorized environment
original filename / archive members / decompression-only transformation
size_bytes / local_sha256
upstream checksum + algorithm + checksum source（如有）
license URL / local license file / required attribution
verification_status / known provenance gaps
```

## 9. 从官方来源准备一次新版本

1. 先运行 `python scripts\bootstrap_data.py plan`，确认当前 12 项和访问限制。
2. 从官方版本页或获授权数据保管者取得新文件；保存 release、许可、checksum 和获取
   记录，不覆盖旧 raw。
3. 完成必要的无损解压。当前锁定布局需要 STRING `.txt`、Reactome GMT、COSMIC TSV
   和 ChEMBL SQLite；builder 不会替你完成 Reactome/COSMIC/ChEMBL 容器展开。
4. 先核对上游 checksum，再计算容器、成员和最终输入的本地 SHA。
5. 建立 `data/source_registry.candidate.json`，填入新路径、size、SHA、来源和分发提示。
6. 用 candidate registry 验证全部 12 项；`--input` 只换路径，不会把旧 hash 变成新
   release 的 hash。
7. 在 `data/processed/candidates/<批次>/` 构建，比较覆盖率、拒绝项、生物学统计和旧
   基线；不要先覆盖 live context。
8. QC 通过后，在同一个 reviewed change 中同步 builder release pins、registry、
   provenance manifest、production config、source/runtime locks、测试和本手册。
9. 按[数据更新操作手册](DATA_UPDATE_GUIDE.md)完成提升、验证和回滚演练。

当前 `bootstrap_data.py rebuild` 是“重建已 pin 版本”的入口：它会先按旧 registry
校验全部 12 项，再写 live context，最后对旧 runtime lock 验证。因此它不是新 release
的 candidate 升级工具。

## 10. 数据包和 Git 发布边界

- 不提交 `data/raw/`、`data/processed/`、普通 PLATO response、COSMIC 原件或 ChEMBL
  数据库。
- code-only 仓库只提交代码、测试、说明、配置样例、registry/lock 元数据和经审核的
  小型标识 registry。
- lock 中出现路径和哈希不表示相应字节可以公开再分发。
- 受限原件与可公开的脱敏派生物必须使用不同制品、不同哈希和不同访问控制。
- 发布前执行 `git status --short`，逐项检查 staged files；不能使用宽泛的
  `git add data`。

数据包的推荐验证顺序：

```powershell
# 快速运行时包；约 7.64 MiB
python scripts\verify_manifest.py --scope runtime

# 12 个当前 builder 输入
python scripts\bootstrap_data.py verify

# 完整 40 文件源审计包；约 35.45 GiB，耗时和磁盘 I/O 较大
python scripts\verify_manifest.py --scope source

# 详细 provenance companion；需要显式指定
python scripts\verify_manifest.py --manifest data\manifest.lock.json
```

如果 code-only clone 在 `--scope source` 下报告文件缺失，先确认是否实际获得了完整
source bundle；不要用空文件、同名新 release 或重新计算的新哈希“修复”lock。
