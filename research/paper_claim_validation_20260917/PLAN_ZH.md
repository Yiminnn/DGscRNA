# 下一步：检验原始 R workflow 的发表主张

这是针对用户“结合 PI 和 reviewer，如何证明现有 paper 效果最优”的具体补实验方案。本轮只核对原始要求、既有结果并制定方案，未启动新的模型拟合。近似复现已满足用户当前要求；历史权重逐细胞恢复不再作为后续实验的前置条件。

**后续用户修订：Darmanis 不再作为新增主验证的 golden dataset。** 用户指出其原标签与目标注释不匹配。主验证改为 **PTC + 已在使用的 GSE274546 GBM**；Darmanis 的既有运行和所有分数保留为补充与局限案例，不新增六档 HVG、大规模调参或所谓细粒度最优性实验。现有证据支持标签分辨率/适用性不足，不支持直接判定全部原始测序数据质量差。详细依据见 [DARMANIS_SCOPE_UPDATE_ZH.md](DARMANIS_SCOPE_UPDATE_ZH.md)。

## 1. 明确要检验的三个主张

1. **流程选择：** 在 PTC 和 GBM 预先列出的候选中，HVG → PCA → UMAP → R HDBSCAN 的终端注释表现是否最高，或与最高者接近。
2. **方法贡献：** 合适的 marker、context 选择和终端 DL/refinement 分别带来多少改善。
3. **整套方法比较：** 在测试标签不参与调优、竞争方法获得适当输入/reference 的条件下，DG-scRNA 是否优于强竞争方法。

这三个主张需要不同对照。内部聚类消融的最高分不能代替与 annotation tools 的比较；按同一批标签挑到的最高分不能代替可部署的选择规则。

原始方法锚点为 LogNormalize、VST/CCA2000（多批次时）、ScaleData、PCA30、uwot UMAP2、R dbscan::hdbscan(minPts=50)、原 density score 与 10-epoch MLP，源码终端阈值为 0.90。稿件所写 0.70 另列敏感性。所有注释结果仍使用终端 DL；合法 no-op 与实际训练分别记录。

**历史事实需保留：** 原始 PTC 八样本共同整合，NMT 最终选 Thyroid/PCA-SNN/none，TTU 选 Pubmed/UMAP-HDBSCAN/mean。不能把“原论文流程有效”改写成“原论文全部组都选 UMAP-HDBSCAN”。新 NMT/TTU 分组整合是另外的实验条件。

## 2. PI 与 reviewer 的真实要求

| 来源 | 原要求 | 本方案对应 |
|---|---|---|
| 08-13 PI 转录，行 63–90、95–104 | 一个 golden dataset；固定 marker 检验降维/聚类；原标签与不同 marker 的注释同图比较 | A、B；按用户后续修订使用 GSE274546，通过与性能排名无关的标签/QC 核验后选展示样本 |
| 08-04 PI B1；08-13 行 119–123 | GBM/浸润细胞相关 marker；从其他已发表数据提取 reference markers | B；保留组织、疾病和测量技术来源 |
| 07-30 PI A3、A4、A6 | 消融和竞争方法分别成文；response letter 引用补充文件 | A/B 与 C 分开交付 |
| Reviewer 2-2、2-5；3-C3 | 模块消融、超参数敏感性、DL 是否传播错误 | A、B、D |
| Reviewer 3-B2 | marker 调优不能读取测试标签；对手应有公平 context 匹配机会 | C，优先级高于扩展更多参数 |
| Reviewer 3-B1、B3、B4 | 多类指标/Unknown；修复 SignacX；强基线与配对统计 | C |
| Reviewer 2-6 | 批次校正与生物背景的关系 | D，复用现有 PTC 结果并补具体缺口 |
| Reviewer 2-9、3-C1 | 实测时间、内存、>100k 单次运行的规模曲线 | E |
| Reviewer 3-A | 实际 MLP 不接收图信息，需改 framing 或实现真实 GNN | 方法表述按现有实现修正；保留原算法不需为名称重造模型 |

来源为本地已保存原文：[reviewer](../../paper/comments.md)、[08-13 转录](../meeting_0813/transcript_full.txt)、[07-30 任务](../meeting_0730/ACTION_ITEMS_0730.md)、[08-04 任务](../meeting_0804/ACTION_ITEMS_0804.md)。旧任务文档的完成状态是当时快照，本表只使用其要求，不将旧缺口自动认定为今天仍未完成。

## 3. A：先做一个严格 R 的流程比较

### A0：冻结评价和执行定义

- 先冻结细胞、基因 ID、marker 文件、reference 来源、注释词表与语义映射。GBM 保持数据集真实标签层级；同义词、父子类型的映射规则对所有方法一致，粗细层级分表。不能看排名后选择某个层级作为主指标。
- PTC 保留原论文 non-T-positive F1/binary AUC 用于原结果比较；同步报告 T-positive F1、TCR-positive recall、coverage 及 Unknown-as-error 指标。原稿 F1 不得继续误称为 T-positive F1，TCR 未检出不得当作已经验证的非 T。
- 原 PTC checkpoint、fresh all8 CCA 和 fresh NMT/TTU CCA 独立成表。原始 all8 结果用于复现锚点；用户指定的两组独立整合用于分组消融。
- GBM 主指标为与固定注释层级匹配的多类 macro-F1，逐类 F1、混淆矩阵、Unknown/coverage 同时报告。ARI/NMI 是聚类结果，不能替代终端注释 F1。

### A1：GBM 的关键缺口

直接在现有 **GSE274546** 队列补原始 R 实现的六档 HVG：500、1000、2000、3000、5000、all。之前完整的大网格属于另一 Python benchmark，不能直接填入原 R 表。沿用事先冻结的 97 样本/55 患者主队列和 121 样本/59 患者完整队列报告范围，不按新 R 分数重新筛选样本。

拟合前先核验 counts/细胞 ID 与原始作者标签的对应、真实患者 ID、标签层级和主要类别支持数。原作者派生细胞状态与通用细胞类型分开评价；无法由 marker 词表支持的状态不冒充通用分类 gold label。marker 缺类属于方法/reference 覆盖限制，不以此宣布整个数据集无效。该检查是准入核验，不预先假设 GSE274546 的每个标签都适合评分。

每档比较 PCA30-SNN、PCA30-HDBSCAN、PCA30→UMAP2-SNN、PCA30→UMAP2-HDBSCAN，固定 marker/cutoff、评分及 DL 条件。另加 HVG 表达空间直接聚类的无降维对照，单独说明其距离与高维计算成本。先核实每臂实际执行的 normalization、feature identities 和输入宽度，再扩展。

PI 要求逐样本讲清楚，因此保留每位 donor/sample 的输出，避免只展示跨患者 pooled 值。已完成 pooled CCA2000 作为明确标注的另一个作用范围保留，不把新 per-sample 结果与之冒充单因素比较。

从已冻结队列中按细胞量、患者/样本结构和标签支持选择资源 pilot，不依据 DG-scRNA 或 UMAP/HDBSCAN 的优劣选择。pilot 核验通过后扩展并行 SLURM。主结果覆盖固定队列全部患者；主图展示样本在新性能评分前确定并说明选择规则。Darmanis 已有结果独立保留，不把两个 GBM 队列合成一个匿名“GBM”总分。

### A2：特征选择的两个不同实验

1. **原流程 HVG 预算：** 同时改变该 R 分支实际使用的 integration/geometry/scoring/DL 基因；回答哪条完整配置表现最好。
2. **单独几何 HVG：** 固定同一表达矩阵、scoring/DL 特征清单及 marker 分母，只改变 PCA/UMAP 输入；回答 HVG 的几何作用。原流程默认得分较高并不能替代这个机制对照。

PTC 已有广泛的两类结果，应先复用。NMT 已有 geometry-only 对照的固定 2000 scoring genes 不含 Thyroid 唯一 T marker CD3D，导致 T 种子缺失；这组全零不能用于判定降维无效。新增机制对照采用统一、预先记录的 marker-retention 规则（例如固定 scoring 基因与完整 marker 基因的并集），在两边使用同一规则，DL 输入固定。该分支标为 marker 保留消融，不改写原基线。

### A3：同等机会的有限敏感性

保留原参数作为不可替换的历史锚点。主网格使用相同固定设置；另给各聚类路线相近的调优预算。补 UMAP 2/10/30 维、R HDBSCAN minPts 与 SNN resolution 的局部曲线，不用调过的 UMAP 对比完全没有调参机会的对手。

重复固定候选的 seeds 0、1、2、3、42，分开表示 embedding 和 MLP 随机性；优先核验原配置、HVG5000 候选和各竞争路线。原 checkpoint 的 DL-seed 检验不能冒称重新做了整合/UMAP 的端到端种子检验。种子只是算法重复，不增加患者数。

UMAP 聚类不必限定为可视化的二维，这是增加 10/30D 敏感性的理由，不能据其官方示例认定本数据上一定更好：[官方说明](https://umap-learn.readthedocs.io/en/latest/clustering.html)。R 的 HDBSCAN 参数按实际 minPts 接口设置：[R 文档](https://search.r-project.org/CRAN/refmans/dbscan/html/hdbscan.html)。

## 4. B：固定流程后检验 marker 与 DL 的贡献

- 对 PTC 复用全部 17 个原库和 3 个 cutoff。GBM 保留正常脑、GBM、浸润免疫/血管及 AllHuman 对照，并核验可用的独立文献/reference-derived markers；源研究含目标病例时明确记录，不称独立验证。
- CellMarker 的不同组织子集不完全等同于 PI 要求的不同发表研究 reference；补齐来源核对后再比较，不能靠换库名称填格。
- 使用 marker 选择方式（固定 vs 在训练患者中选择）× refinement（无 vs 有）的 2×2 比较。保持比较所需的词表、cutoff 规则和候选预算透明；不要同时偷偷改变 cutoff、marker universe 和 DL，再将总变化归给一个模块。
- 固定完整原算法的 density 公式。原代码没有的 alpha/Jaccard 混合或 epsilon 迭代先修正稿件描述，不为解释稿件而虚构一次“参数消融”。
- 保留全库结果分布、panel/细胞类覆盖与 marker 丢失。Unknown 独立分析其表达、QC/doublet 线索与 TCR 支持，不能直接给 Unknown 编造新细胞类型。
- 一张固定显示坐标的图依次展示原标签、各 reference-marker 的初始注释、最终 DL 注释及 Unknown；参与实际聚类的高维 embedding 与二维展示坐标明确区分。

## 5. C：发表优越性最关键的一步——公平选择与竞争方法

在同一患者所有组织/样本保持同一 fold 的前提下，做患者留出。真实患者 ID 先由原始 metadata 核验，不能把 sample 名当作独立患者。PTC 和 GBM 可在已有队列内完成此步骤，不要求额外独立 Pu 队列。

训练患者用于选择 HVG/marker/cutoff/聚类参数；留下患者的真值仅用于评分。MLP 使用 marker 生成的标签，和用于配置选择的 curated/TCR 标签用途分开记录。优先每个目标样本独立执行无监督和 marker/DL 注释，满足 PI 的使用情境；若某项保留全队列无监督整合，必须明确为 transductive、仅隔离评价标签的实验，不能写成未见患者独立重训。

目前所有病例的分数都已用于历史探索，因此这属于固定方案后的回顾性患者交叉验证；不能把它包装成从未被研究者看过的新盲测集。现有同标签最大值保留为探索性上界，单独列出实际可用的选择规则及其与上界的差距。

优先对齐 scType、scCATCH、SCINA 等可用 marker-based 对手，复用 SingleR/CellTypist/CHETAH/scmap 等既有运行中满足同一评价协议的预测；核实 SignacX 输入/label mapping 后决定是否保留。Reviewer 要求的 GNN 对手作为独立方法比较核验可运行性，不用给 DG-scRNA 添加 GNN 来改变待验证原算法。

所有工具使用同一保留细胞和固定评价词表，但按各工具要求使用合适的 counts/normalized 输入。marker 方法获得同等合适的 marker 候选和选择机会；reference-based 方法明确其获得的标签/参考信息，按信息条件分组并完整展示，不隐去强对手。

**优先纠正旧比较的 coverage 分母问题：** `results/methods/GBM_DOC2.md` 中 CHETAH/scmap/DG-scRNA 的部分分数只评价被调用细胞，coverage 不同。主指标应覆盖全部测试细胞，Unknown 作为未正确分类纳入分母；已调用细胞上的条件性指标只能作为附表。任何工具的失败、no-op、无支持类别均透明记录，不能把安装或输入错误当科学失败。

统计按患者配对并报告效应量和区间，多次比较进行预先规定的校正；报告与最强有效对手的绝对差和相对差，不只挑最弱基线。PTC 患者数量较少，展示全部患者，不把细胞或随机种子当独立生物重复。

## 6. D/E：完成 reviewer 的必要补充

- **Batch：** 复用 PTC 无校正/CCA/Harmony 与生物保留诊断。CCA 与 Harmony 改变算法及作用空间，应解释为实际工作流比较，不能称纯“校正先后位置”效应。共享细胞类型内考察混合，防止将组织状态消除当改善。
- **参数：** DNN 宽度、epoch/训练曲线、0.70/0.90 置信阈值等按实际源码逐项敏感性；已有训练历史和概率能回答的问题直接复用，不重复拟合。
- **规模：** 完整 HCL 分成 59 个 tissue 的运行不等于单次 >100k。若保留该主张，在可用数据上按固定抽样规则做单次规模曲线（例如 10k/30k/50k/100k/>100k），记录完整 pipeline wall time、CPU/GPU peak memory、冷启动/后续重复和硬件；与有效对手同时报告。优先 pilot，再决定 SLURM 内存与并发。
- **方法文字：** 按实际 MLP calibration pipeline 解释图、聚类和 DL 的连接；0.70/0.90、F1 正类、binary AUC、Accuracy 来源和细胞数需与可计算定义一致。最优性实验不能修复错误的方法描述。

## 7. 预先规定如何下结论

- 原 HVG2000/UMAP-HDBSCAN 在约定主终点及条件中最高：可写 **best among the evaluated configurations**，注明数据集、条件和终点。
- 只有 HVG5000 或不同 marker 最高：称优化后的配置；原配置作为原方法锚点保留，不声称其具体参数也第一。
- 若差异小且区间宽：报告接近最高及区间。若要正式声称 non-inferior，需在验证运行前依据科学意义规定容许差值；“不显著”不等于已证明非劣效。
- 若 NMT 仍由 SNN 更好：保留原始按组分支，限定 UMAP-HDBSCAN 的优势条件；不剔除 NMT，也不靠合并总分隐藏例外。
- 若竞争工具更好：报告其优势与所需 reference 信息，准确限定 DG-scRNA 的贡献。内部消融最优不自动等于现有方法中最优。

## 8. 执行次序与交付

**A0 → GSE274546 的标签准入/原 R pilot → A1/A2 → B/C → PTC 缺口与 D/E。** Darmanis 不再是前置任务。优先修正评价和选择协议，复用所有可验证缓存。科学计算、数值核验、绘图及统计全部经 SLURM，并行数组限制并发；OOM/timeout 先调整资源，不改变科学条件。

交付沿用原 notebook 和 OneDrive 目录：Supplement A 为一个 golden GBM 的流程/marker/DL 消融；Supplement B 为公平多类方法比较；资源/敏感性及完整网格作补充表。原 English decision tree 标明历史默认、调优候选、已验证与未验证节点；response letter 每项指向对应表图。旧数据和原 notebook 内容保留，不另做网页。

本方案不预先承诺获胜结果。它将用户希望保留的论文核心拆成可直接接受或否定的具体实验，并让结果对应 PI 与 reviewer 的实际要求。
