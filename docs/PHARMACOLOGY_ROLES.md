# 文献角色注释：方法、证据与边界

5.4 增加 `drug_specific_literature_roles_v1`。目的：让药物特异的药效关系与药物代谢、转运、胞内保留分开进入叙述，并让读者看见角色未知及家族级证据。

## 文献依据

| 来源 | 用于本实现的依据 | 不支持的推论 |
| --- | --- | --- |
| [Gregori-Puigjané et al., PNAS, 2012](https://pmc.ncbi.nlm.nih.gov/articles/PMC3396511/)；DOI 10.1073/pnas.1204524109 | 生物活性关系和作用机制靶点需区分，预测后的文献与实验确认另行处理 | 任意结合或相似性预测都等于治疗作用机制 |
| [Drug and Clinical Candidate Drug Data in ChEMBL, J Med Chem, 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12516679/)；DOI 10.1021/acs.jmedchem.5c00920 | 药物作用机制有独立的策展记录及引用；一般 bioactivity 不能直接替代 MoA | 所有数据库相互作用都是主要药效靶点 |
| [Izbicka et al., Cancer Chemother Pharmacol, 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2728224/)；DOI 10.1007/s00280-009-0954-4 | Pralatrexate 的 DHFR 抑制与摄取、FPGS 多聚谷氨酸化承担不同环节 | 这些单药实验已经验证 HCT116 中的组合机制 |
| [Marre et al., Cancer Res, 1996](https://pubmed.ncbi.nlm.nih.gov/8640817/) | 肝微粒体与肝细胞实验支持 CYP3A 参与 docetaxel 生物转化 | CYP3A 参与代谢意味着 docetaxel 的主要抗癌作用是抑制 CYP3A |
| [Shirakawa et al., Jpn J Cancer Res, 1999](https://pmc.ncbi.nlm.nih.gov/articles/PMC5926029/)；DOI 10.1111/j.1349-7006.1999.tb00723.x | 人 P-gp 转染体系支持 docetaxel 是 ABCB1 的转运底物 | docetaxel 本身必然通过抑制 ABCB1 逆转耐药 |
| [OATP1B3/MRP2 与 docetaxel 毒性动力学研究, 2014](https://pmc.ncbi.nlm.nih.gov/articles/PMC4120017/) | 转运体表达体系显示 OATP1B3 介导摄取 | 所有肿瘤和细胞系中的摄取程度相同 |
| [Docetaxel 药品说明书，12.1/12.3 节](https://www.dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid=654319ae-e5f9-4a6a-b566-201923153b6e) | 微管稳定作用与 CYP3A4 代谢；补充精确药物的药理来源 | 将微管家族的说明升级为人 TUBB3 亚型结合验证 |
| [Azzoli et al., Clin Cancer Res, 2007](https://pubmed.ncbi.nlm.nih.gov/17473201/)；DOI 10.1158/1078-0432.CCR-06-1754 | 提供 Pralatrexate/紫杉类组合研究背景 | 作为本次 HCT116 机制评分或外部验证输入 |

这是有针对性的检索，不是系统综述。角色划分是依据文献制定的工程策略；论文没有验证本实现的规则、分数或性能。

## 角色是药物—蛋白关系的属性

不能按 `CYP`、`SLC`、`ABC` 前缀或蛋白类型自动判定。酶或转运体也可能是另一种药物的直接药效靶点。同一药物—蛋白关系可以有多种角色，输出允许重叠分组。

- `pharmacodynamic_target`：文献指定的药效作用关系；`exact_protein` 与 `target_family` 分组展示。
- `drug_metabolism`：药物的代谢关系。
- `drug_transport`：药物摄取或外排。
- `drug_activation_retention`：药物活化或细胞内保留环节。
- `context_association`：仅有上下文关联，不提升为直接药效靶点。
- 未匹配记录：`unresolved_role` / `unknown`。未知不表示没有作用。

首批资源只有两种精确结构、七条关系：Pralatrexate 的 DHFR、SLC19A1、FPGS，以及 docetaxel 的 TUBB3、CYP3A4、ABCB1、SLCO1B3。DHFR 是药效靶点；FPGS 在本条文献记录中描述胞内多聚谷氨酸化/保留环节，不由此断言它没有其他作用。TUBB3 仅映射到已知的微管家族作用，不把家族陈述当成该亚型的直接结合证据。其他 PLATO 候选（例如 EGFR、ERBB2）没有匹配角色依据时保持未知，不被硬编码为无关。

## 结构、来源及范围

`data/pharmacology_roles_v1.json` 保存参考文献和记录。每条记录要求：稳定 ID、药物名、完整 InChIKey、精确基因符号、角色、证据粒度、作用类型、过程、简短证据说明、引用 ID 和实验上下文。

运行时先对提交结构执行既有标准化，再按**完整 InChIKey + 精确候选基因**匹配。不按药名、结构前缀、蛋白名称或通配符匹配；盐型、立体异构体或别的结构未匹配时不自动转移注释。只有现存合格 PLATO 候选才能获得注释，目录中有而 PLATO 中没有的靶点不会被加入。

资源的预期 SHA-256 写在默认配置中；自定义资源也必须提供哈希。缺失引用、重复关系、非法角色/粒度、部分结构标识、通配符和字节改变均拒绝加载。哈希用于确认运行所用字节，不证明科学结论正确。首次资源属于开发期文献策展，不宣称经过独立盲评。

药物特异的角色、引用和实验上下文写入 `ContextScore.pharmacology`，与原始 `source_evidence` 分开。PLATO 源证据和 `predicted_only` 状态不被篡改，也不借文献角色绕过 PLATO 输入检查或原有方向/可信度门槛。

## 对解释和统计的影响

摘要中的代表项改为有文献支持的药效关系；家族项显式标记 `[family evidence]`。没有已分配药效角色时会说明缺失，不把上下文排名第一名自动当成主要药理靶点。完整 `top_targets` 继续保留上下文排序，并声明其语义。

`pharmacology` 输出药效、家族药效、药物处置、上下文关联和未知分组，并为显著通路附上角色组成。若通路的全部命中候选都只有药物处置角色，则标记 `drug_handling_only=true`；含未知、直接药效或缺失成员时不能这样标记。每个模块、证据路径和药物—候选边也保留角色信息。

不根据“ADME”字符串删除通路；通路标签本身不是药物机制，成员角色才是解释依据。这里不计算新的“药效靶点富集 p 值”：原有 q 值仍对应完整 PLATO 候选集合，`pharmacodynamic_enrichment_tested=false`。

角色层不改变上下文排名及其数值、不添加候选、不重选统计集合、不改变原有机制标签、机制分或可信度分。原始网络模式和候选上下文描述作为单独字段保留，避免将文献锚定与统计支持混为一谈。网络方法仍为 v3，统计合同仍为 1.1；新增的是结果 schema 5.4 和角色策略 v1。

证据分组可以提高叙述可读性，但不是药效或组合因果链的验证。剂量、时间、细胞暴露、靶点占有率、组织表达和组合干预证据仍需单独获得。使用这些文献设计注释后，不能再将相同文献当作独立验证集。

## 配置与扩展

默认启用随包附带的资源。`pharmacology_roles.enabled=false` 可关闭；对照运行应使用相同网络预算、种子、靶点及上下文，比较统计和数值是否一致。替换资源时使用 `pharmacology_roles.path` 与 `pharmacology_roles.sha256`，同时更新版本、引用和记录范围。

运行不自动联网补注释。新增药物需要检索直接支持该药物—蛋白关系的来源，记录结构身份与适用范围，并保留未知及互相矛盾的证据，不能只收录支持预期结果的关系。当前框架保留多条注释，不提供自动的冲突裁决或作用强度合成。
