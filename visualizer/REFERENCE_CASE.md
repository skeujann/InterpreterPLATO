# BRAF / MEK 文献实例

首页默认展示 **Dabrafenib + Trametinib**。这个组合有可查的原始研究，靶点层级清楚，
适合用一个小网络演示关系筛选、路径高亮、原始属性、来源查看和离线导出。
组合与 BRAF V600 突变黑色素瘤的研究背景见
[Flaherty et al., NEJM 2012, PMID 23020132](https://pubmed.ncbi.nlm.nih.gov/23020132/)。

这是一份**人工整理的公共注释实例，不是 InterpreterPLATO 或 PLATO 的运行结果**。
没有生成协同分数、疾病邻近 p 值、通路富集、细胞系预测或患者层面的结论。
图中的路径和位置不能继承临床试验的有效性结论。

## 在界面中查看

1. 双击 `index.html`，或打包后的单文件 HTML。
2. 默认“疾病主路径”显示双药、MEK1 / MEK2、BRAF 和疾病模块共 6 个节点，
   以及 3 条靶点关联、2 条已有 PPI、1 条疾病背景归属。
   桌面从左向右经 BRAF 读到最右侧疾病模块；窄画布使用紧凑布局，疾病模块在下方。
   BRAF 是关键中间节点，不代表共同直接靶点。疾病模块仅为文献背景，阅读顺序不代表作用方向或因果终点。
3. 展开“三步读懂这张图”。“1 靶点”高亮三条文献药理靶点关联：Dabrafenib—BRAF，Trametinib—MEK1/MEK2。
4. “2 连接”高亮 `MAP2K1 — BRAF` 与 `MAP2K2 — BRAF` 两条既有无向互作，不表示激活、抑制或协同。
5. “完整网络”显示全部 11 个节点、20 条边、6 条选定路径。“3 背景”展开完整网络并高亮通路成员与文献背景注释，继续使用无箭头点线。
6. 点击“查看来源与限制”打开来源卡片；“实例 JSON”下载完整的可审阅文件。

实例布局只用于随代码打包的这份对象，普通上传文件仍使用通用布局。
可选择另外四种布局。切换到“合成测试 · 功能边界”仍能演示未映射、断连和共享零跳。

## 关系来源

冻结整理日期：**2026-10-21**。文件是
[`examples/dabrafenib-trametinib.json`](examples/dabrafenib-trametinib.json)，附精确字节 SHA-256。

| 内容 | 来源与范围 |
| --- | --- |
| Dabrafenib—BRAF | [NCI Drug Dictionary](https://www.cancer.gov/publications/dictionaries/cancer-drug/def/dabrafenib)，选定药理靶点，不是完整靶点谱 |
| Trametinib—MEK1/MEK2 | [NCI Dictionary](https://www.cancer.gov/publications/dictionaries/cancer-terms/def/mekinist)，选定药理靶点，不是 PLATO 预测 |
| 九条无向蛋白互作 | UniProtKB 七个人类条目的 `INTERACTION` 记录；保留其 IntAct 标识和实验计数 |
| 七个通路成员 | [Reactome R-HSA-5673001](https://reactome.org/content/detail/R-HSA-5673001) 的参考实体参与者子集；不是富集分析 |
| BRAF 与黑色素瘤背景 | 上述原始研究的 BRAF V600 背景，人工整理的**单基因文献注释**；不是 COSMIC 或计算得到的疾病基因集 |

预先选定七个 canonical UniProt accession：BRAF/P15056、MAP2K1/Q02750、MAP2K2/P36507、
MAPK1/P28482、MAPK3/P27361、RAF1/P04049、KSR1/Q8IVT5。
只保留两端 accession 均精确匹配、`organismDiffer=false`、有实验记录的非自环互作，
按无向蛋白对去重，不折叠异构体，不从 Reactome 通路邻接或生物学常识补出互作。
选定路径只复用这些边，不声称枚举了全部路径。

| 无向蛋白对 | UniProt `numberOfExperiments` |
| --- | ---: |
| BRAF — MAP2K1 | 65 |
| BRAF — MAP2K2 | 13 |
| BRAF — RAF1 | 73 |
| KSR1 — MAP2K2 | 10 |
| KSR1 — RAF1 | 5 |
| MAP2K1 — MAPK1 | 2 |
| MAP2K1 — MAPK3 | 2 |
| MAP2K1 — RAF1 | 39 |
| MAP2K2 — RAF1 | 7 |

实验计数按单个 UniProt 注释记录原值呈现，不把两端重复记录相加，不作置信度或边权。
互作注释不一定代表直接结合，也未在该突变、组织或临床场景逐边验证。
记录的通路参与不等于药物组合的通路富集，网络连通不等于协同、因果或临床机制。

`reference_provenance` 保留所选 UniProt 原始字段、条目版本、官方 API URL 和原始响应
SHA-256，以及 Reactome 参与者查询 URL、所选成员和响应哈希。完整上游响应保存在本地
整理快照中，不随轻量查看器分发；在线条目可能更新，不能把未来 API 响应与旧哈希不一致
解释成此实例已被验证或证伪。公开 JSON 中的所选事实与版本可独立审阅。

## 署名与许可

- **UniProt Consortium**：选取和整理的 UniProtKB 标识符、互作事实及 IntAct 交叉引用。
  原数据按 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) 提供，
  见 [UniProt 许可说明](https://www.uniprot.org/help/license)。修改包括字段选择、无向去重、
  中文界面说明和显示别名；未转载原始功能描述正文，也不代表 UniProt / IntAct 背书。
- **Reactome**：只选取通路标识和七个参与蛋白事实，按
  [CC0](https://creativecommons.org/publicdomain/zero/1.0/) 使用，
  见 [Reactome 许可说明](https://reactome.org/license)。不复制图示、品牌素材或数据库 dump。
- **NCI 与原始研究**：引用药理靶点事实、组合研究身份和链接，使用自行撰写的说明，
  不转载论文正文、摘要、表格、图片或说明书。

这些第三方注释的许可不改变项目代码在根目录 `NOTICE` 中的授权状态。
界面、可下载 JSON 和单文件分发均保留来源与许可链接。

## 维护

更新公开事实 JSON 与其 SHA-256 后执行：

```bash
node visualizer/prepare-reference.cjs
npm test --prefix visualizer
node visualizer/build.cjs
```

`prepare-reference.cjs` 验证冻结 JSON 的哈希并生成同内容的 `reference.js`，不发起网络请求。
查看器与 CI 均使用冻结数据；测试比较两份内容，并检查每条 PPI 的记录、路径与标识边界。
