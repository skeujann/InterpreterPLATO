# InterpreterPLATO Visualizer

轻量、纯前端、可离线运行的解释 JSON 查看器。无需部署模型、Python 服务或数据库。

## 直接使用

下载仓库源码后，双击本目录的 `index.html`，并保留同目录的 JS、CSS 和 `vendor/` 文件。
浏览器会先显示一个明确标注的合成概念示例；点击 **打开解释 JSON**，或拖入本地文件，
即可查看已有解释。支持单条解释、解释数组，以及 `{"results": [...]}`。

文件只在当前浏览器内通过 File API 读取。应用不向服务器上传数据、不查询 PLATO、
不下载原始生物医学资源，也不重新运行解释器。Cytoscape.js 已随目录打包，没有 CDN、
远程字体、分析脚本或运行时安装步骤。已有 Python 环境也不是查看器的运行前提。

## 功能

- 根据 `evidence_graphs.explanation_evidence_graph.nodes/edges` 展示既有证据图。
- 节点搜索、节点类别与边类别筛选、四种布局、缩放和平移；搜索后按 Enter 查看节点详情。
- 高亮 `local_evidence.cross_drug_paths` 与 `local_evidence.top_paths` 的已有 PPI 边。
- 区分共享节点零跳、疾病成员零跳和非零跳路径。
- 悬停／点击查看原始属性与来源；独立诊断面板保留上游未映射、断连和展示限制信息。
- 选择批次中的单条记录；失败记录仍可查看诊断，不伪造空机制。
- 导出带运行模式和解释边界标识的 PNG；可下载内置合成示例 JSON。

上限为 20 MiB JSON、100 条解释、每条 5,000 个节点和 20,000 条边。
这是面向结果文件的轻量查看器，不是原始数据库或大规模全网络浏览平台。
超过 400 个可见节点时，默认力导向布局会改用网格布局，避免不必要的前端计算。

## 证据与视觉语义

| 关系 | 样式 | 解释 |
| --- | --- | --- |
| `PPI` | 青绿色实线、无箭头 | 输入文件记录的无向 PPI 连接 |
| 候选靶点关联 | 金色虚线、无箭头 | 文件中的药物—候选靶点关系 |
| 疾病／通路集合归属 | 紫色点线、无箭头 | 集合成员关系，不是 PPI 或作用方向 |
| 证据／报告结构 | 灰色虚线、默认隐藏 | 证据支持、报告归属或敏感性记录 |
| 动态假说 | 粉色虚线、默认隐藏 | 规则生成假说，不解释为已验证的因果作用 |
| 未识别类型 | 灰色关系、默认隐藏 | 保留原类型并提示，不归类成 PPI |

查看器不根据路径列表、节点位置、相近名称或示意图推断、补全任何关系。
一条路径缺少节点或对应 PPI 边时，路径仍列出诊断，但禁用高亮；不自动补边。
无证据图的旧结果会明确提示缺失，不从靶点列表拼出一个网络。
重复节点 ID、悬空边、有向／自环／非基因 PPI 作为输入错误拒绝展示。

上游 `build_evidence_graph` 会检查路径是否存在于其加载的 PPI 图中；前端只能呈现
JSON 中的记录，**无法独立核验上传文件的真实性、上游原始数据或再分发授权**。
界面和 PNG 均保留这一边界；图中连通不构成真实机制、协同、激活／抑制或因果证明。

内置示例的药物、`SYN_*` 节点、PPI、通路与诊断均为手工构造的合成占位：
`run_mode=test`、`generated_from_test_fixture=true`、`formal_result_eligible=false`、
`metadata.conceptual_example=true`、`metadata.calculation_performed=false`。
它只演示交互功能，明确标注“概念示意（非真实计算结果）”，不冒充真实 PLATO 预测。
真实输入的测试模式会显示为测试模式，不会因 `run_mode=test` 就被改称合成数据。

## 可选：构建单文件 HTML

开发者安装 Node.js 22 或更高版本后，从仓库根目录执行：

```bash
node visualizer/build.cjs
```

输出到 `visualizer/dist/InterpreterPLATO-Visualizer.html`。将这个单文件交给其他人后，
对方双击即可打开，无需 Node.js、Python、其他资源文件或网络连接。
第三方库的 MIT 许可证同时保留在单文件中。

## 开发检查

```bash
npm ci --prefix visualizer --ignore-scripts
npm test --prefix visualizer
```

运行时只有已打包的 Cytoscape.js；`jsdom` 仅用于开发期 DOM 交互测试。
测试涵盖图的完整性、路径不补边、零跳、无箭头视觉规则、输入上限、批次失败记录、
本地文件读取、文本安全、布局交互、模式标识和单文件构建。
DOM 与 Cytoscape headless 测试不能替代真实浏览器的 Canvas 和视觉验收。

## 第三方依赖

Cytoscape.js **3.34.3**，MIT 许可证；官方来源、npm 完整性摘要和文件 SHA-256
保存在 `vendor/dependency.json`，许可证在 `vendor/CYTOSCAPE-LICENSE.txt`。
可选构建在打包前验证第三方文件的哈希，不静默下载或更新依赖。
保留其许可证不改变 InterpreterPLATO 项目本身在 `NOTICE` 中声明的授权状态。
