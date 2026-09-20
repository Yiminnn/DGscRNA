# Reviewer 补齐计划 · 待确认

状态：2026-09-20，仅制定计划；本轮未提交新实验、未改已有结果。

目标：在现有版本上关闭 Reviewer 1/2/3 的全部 35 项要求。关闭可以是补足实验、整理已完成证据、修正原文，或对缺乏独立验证的主张作有依据的限定；不把“程序跑完”当作 reviewer 要求已满足。来源均为 2026-07-09 决定信，逐项原文位置、现状、工作包和验收条件见 `REVIEWER_COVERAGE.tsv`。当前 PI 指令优先采用 08-13，再保留不冲突的 08-04、07-30、07-21 要求。

## 固定边界

- 以原始 Rmd、已核验 R helper 和恢复的原始数据为准。原锚点为 LogNormalize、HVG/CCA 2000、PCA30、uwot UMAP2、R HDBSCAN minPts50、实际 density score、256/128 MLP、10 epochs、终端阈值 0.90。HVG3000 和阈值 0.70 是已命名的对照；不把 `min.cells.to.block=3000` 写成 HVG 数。
- GBM 单样本 RNA 分支的 scoring 使用全部 RNA genes；HVG 影响 geometry 和 DL 输入。PTC CCA 分支的 integration/geometry/scoring/DL 特征分别记录。不能把不同阶段的“全基因/2000基因”混为一个开关。
- 先补齐 GBM，再进入 PTC 和公共数据的剩余计算。沿用 `DGscRNA` 仓库的 `align-r-reference`，保留现有未跟踪文件与其他分支。
- 原 `notebooks/dgscrna_results.ipynb` 继续追加，所有旧输出保留；Paseo HTML 只是阅读导出。继续整理原 OneDrive 目录，不另建替代结果本或新版网页。
- DG-scRNA 的方法结果采用最终 DL/refinement；初始 marker-only、clustering 和合法 no-op/invalid 单列。
- PTC 新分析固定 NMT＝MT-1/2＋N-1/2；TTU＝TU-1/2＋T-1/2。历史 all8 基线单列；NMT 的 Thyroid/PCA-SNN/none 与 TTU 的 Pubmed/UMAP-HDBSCAN/mean 原始选择保留。S2、S3 是不同的原始终点，不能合并成一套人为标签。
- 不新增独立 Pu/GSE184362 验证队列；已在原库里的 GSE184362 markers 保留。Darmanis/brain_GBM 已有结果保留，不作为新增调参或主验证队列。
- 主张限定为预先列出的数据、候选和终点，不预设 UMAP/HDBSCAN 或原始 2000 HVG 必须排名第一。全部有效比较结果纳入完整补充资料。

## SignacX：已按用户确认锁定

**用户已确认：SignacX 按原稿保留 0 T cell；scType 保留实际结果。**

1. 历史 paper 表继续保留 SignacX 的原始 `None` 条目，并明确标注 0 predicted T cells。原 S3 native labels 和 binary flags 不改。0 T calls 不等于 0 accuracy/AUC，也不等于程序没有输出。
2. scType 原稿总体 F1＝0.9392，不改为零。历史 F1/AUC、严格 T-positive 评价、普通 accuracy 分别写明定义；无法定位来源的历史 Accuracy 保留为 paper-reported，不能冒充已复算值。
3. R3-B3 补输入、原函数、版本、native labels 与映射审计。旧结果确实无 T/NK calls；目前不能把它归因于已经证实的“列选错”或“graph 缺失”。
4. 后续缓存核验/公平新比较独立成表，不覆盖历史 0T 表，也不把新预测强制改为零。若某个 SignacX 运行不能满足可比输入和词表要求，则保留其历史行和原因，但不作为支持新优越性 headline 的有效对手。

详细依据：`SCTYPE_AUDIT.md`（文件名沿用最初问题，结论已按用户澄清更新）。

## 已完成，直接复用

| 已有证据 | 当前边界 |
|---|---|
| GBM 121 样本、59 患者；97 样本/55 患者主队列；六档预算×四路线，共 726 单元、2,904 partitions 和终端结果/图 | 原 R 主比较已经完成；不把历史 Python 7×3 当成原 R 终端比较 |
| GBM 363 几何控制、180 表示控制、54 MLP 控制 | 表示/参数控制仅覆盖三个预先按规模选择的样本，不能改称全队列鲁棒性 |
| GBM 患者留出选库、marker×DL 四方、五个对手及等预算比较、Unknown/错误传播 | 保留信息条件和标签支持限制；不是 11 个公共集全部对手均已公平运行 |
| PTC 原基线对账、17 原库×3 cutoff、分组/all8、CCA/NONE/Harmony、seed/marker-retention、四患者统计 | 后续对手缓存重评不等于同信息重拟合；不能用 TCR 验证非 T 类或 T 亚型 |
| 11 个公共数据集原 R 四路线、多组织 marker、终端结果 | 100 个分析单元/12,048 终端条件的旧 campaign 包含 PTC；HCL 分 59 个 tissue，不等于单次 60 万细胞 |
| 3 方法×5 规模×3 独立进程＝45 次 CPU 资源测试，最大单次 120k | CPU wall-time/MaxRSS 已有；GPU 显存为 N/A，完整 grid 搜索成本尚需另外汇总 |
| PTC S2/S3 逐样本 module 表达、患者组成、T-core/TCR 定义敏感性 | 是支持性表达分析，不能替代独立亚型真值、真实 TCR 假阴性率或功能验证 |

## A. GBM 流程缺口

对应 R2-2、R2-5；衔接 PI-06/07/15 和最近 meeting 的 DGCyTOF 比较。

**A1：完成 R 下游统一的 embedding 对照。** 以冻结的 121 样本为范围，主报告保留 97 样本/55 患者。

- 主预算 HVG2000；预先指定 HVG5000 为敏感性，不重新按全测试结果挑预算。六档预算的原 R 四路线旧结果全部保留。
- 七种表示：noDR、PCA2、FA2、ICA2、Isomap2、UMAP2、t-SNE2；三种聚类：K-means、GMM、R HDBSCAN。另保留已有 PCA30/SNN、PCA30/HDBSCAN、PCA30→UMAP2/SNN、PCA30→UMAP2/HDBSCAN 原锚点。
- 统一对照的表示输入固定为同一 scaled-HVG 矩阵；noDR 就是该矩阵。原 R 的 PCA30→UMAP2 作为原锚点单列，不能将 direct-HVG UMAP 偷换成历史流程。这个比较同时呈现配对表示和原流程分支，而不声称每项都是仅改变一个数值参数。
- 复用相同细胞、scoring genes、DL genes、CM2_glioma_other、mean cutoff 规则和终端定义。改变聚类后 DEG、seed 与数值 cutoff 可以随原算法改变，记录为整条分支的总效应。
- K-means/GMM 的 K 候选固定为 5/10/15/20/30/40，由当前 fold 的训练患者终端 macro-F1 选择，同分取较小 K；不读取测试患者标签。小样本不满足 K 的分支按预定合法性规则记录。噪声/未知细胞保留在分母。
- 两个预算的主比较为 121×2×7×3＝5,082 个 selected partition 条件。包含上述 K 候选的计算清单上限为 121×2×7×(6＋6＋1)＝22,022 个候选条件（缓存复用前）；额外候选成本计入搜索账本。只使用固定 context，不对每格重做 16 库×3 cutoff 的全搜索。完全一致的缓存通过输入/参数 hash 后复用。
- 每种 partition 都给图和终端 DL 注释图；同一展示坐标不能被误写为实际聚类输入。队列统计按患者配对，单样本图保留 NL022 等事前展示规则。

**A2：原 R no-clustering 对照。** 在上述冻结队列和 PTC 两组增加去掉图/簇信息的 cell-wise marker seed→同一 DL 支线。原 density scorer 依赖 cluster DEG，去聚类后不能直接调用同一接口：先冻结每细胞 marker 评分、种子/Unknown 规则、训练内阈值选择与 marker universe，再执行。该支线检验“cluster-DEG 种子构造”的贡献；明确替换了种子生成机制，不能宣称纯粹移除了 clustering 而其他算法完全不变。

验收：原 R 锚点回归一致；固定清单每项有有效结果或有依据的结构性状态；OOM/代码错误必须修复重跑；终端标签、图、患者配对表、资源记录齐全。

## B. 邻居数与训练验证

对应 R2-5、R3-S2/S4。复用已完成的 minPts25/50/100、SNN resolution、UMAP2/10/30、MLP 宽度、5/10/20 epoch 和五 seed 结果。

- 分别改变 SNN `k.param=10/20/40` 和 UMAP `n.neighbors=15/30/60`；20、30 为现原实现锚点。一次只改一个旋钮。
- 首先覆盖既有三个 GBM 规模 pilot（TKU4163、NL022、SN040）×HVG2000/5000，以及 PTC 两组各自的原 marker/route。若要将敏感性写成全 GBM 结论，再按冻结规则扩展，不用三个样本冒充全队列。
- 以相同初始化/训练划分增加逐 epoch 的训练与 held-out pseudo-label validation loss/accuracy，最长 30 epoch；保存 5/10/20/30 checkpoint。原 10-epoch endpoint 不变；现有宽度实验不重跑。
- seeds 固定 0、1、42。GBM 使用 UMAP2_HDBSCAN_R＋CM2_glioma_other/mean；PTC 使用 NMT Thyroid/PCA-SNN/none 和 TTU Pubmed/UMAP-HDBSCAN/mean。GBM 三 pilot×两预算×三 seed＋PTC 两组×三 seed，共 24 个学习曲线条件。合法无可训练分支明确记 no-op，不伪造曲线；GBM 备用 CM2_primary_all_context/mean 仅在原条件无训练时独立增加，另记条件数和原因。
- pseudo-label validation 只检验对 marker 种子的泛化，不能叫独立生物真值验证。epoch/阈值选择只看训练患者和内部验证，不用测试真值早停。原代码未实现的 α、ε通过修正文稿解决，不新造算法来补表。

验收：真实日志、逐 epoch 验证曲线、原锚点/新增 checkpoint 指标、seed 变异、适用范围齐全；不凭“不显著”宣称参数全面最优。

## C. 公平方法比较

对应 R2-1/3、R3-B1/B2/B3/B4/S7；先 GBM，随后 PTC 和公共集按独立队列并行。

- 逐格盘点既有预测，只重跑输入、词表、细胞集合、参考信息或选择协议不满足新版对比的缺格。已有 GBM 五对手和等预算比较保留，不整套重跑。
- 公共集清单沿用：HCL、baron_human、blood_DLBCL、brain_GBM、breast_TNBC、colorectal、immune_ALL_human、kidney_ccRCC、muraro、segerstolpe、xin。brain_GBM/Darmanis仅整理已有输出；其他数据补必要对手缺格。PBMC/BM 子集沿用 immune_ALL_human，四套 pancreas 已存在，不因 reviewer 举例再另下载同类队列。
- 方法清单包括 DG-scRNA、scType、scCATCH、SCINA，以及适用的现有 SingleR、CellTypist、CHETAH、scmap、scDeepSort、SignacX 缓存/实现。按工具适用性冻结 dataset×method 矩阵；不将无匹配 reference 或无输出类别的模型伪装成同信息有效对手。scDeepSort 已满足一个实际 GNN 对照，不为 DG 添加 GNN 来改变原算法；需要其他组织模型时仅用来源清楚且无目标数据泄漏的已有 reference。
- 每个数据保留多套解剖/疾病/邻近组织/免疫血管/AllHuman marker 候选。PTC 使用原 17 库，包括外甲状腺组织与 HPA 等；GBM 主选择保留已审计 13 库，来源与作者标签重叠的库单列。
- 同一患者所有样本进入同一 fold。训练患者选择 marker/参数；测试患者只评价。无 donor ID 的数据按可验证 study/batch 划分，并明确不能做患者层级推断；若连独立 study/batch 都没有，则使用预固定配置或另一数据集选参，不能虚构独立留出。只能描述性分析的结果明确标注。队列内已看过标签的分析标为回顾性，不能改称全新盲测。
- marker 方法获得相同合适候选与可比选择预算；reference/pretrained 方法分组报告其信息条件。固定-marker 表与训练内选择表并列，不把训练标签选库包装成无标签自动最优选择。
- PTC 新公平比较固定同一 92,404 历史细胞交集/清单与预声明 TCR 终点，NMT、TTU 分组及患者边界不变；按工具要求使用 counts/normalized 等正确输入。原 paper 表（含 SignacX 0T）独立保留。
- 先用既有标签留出的 transductive 结果作复用表；若保留“未见患者直接泛化”表述，再补留出患者不参与整合拟合的对应训练/投影或逐样本运行，单独命名。不能以未看测试标签冒充完全 inductive。

验收：每个方法有输入/信息条件说明、最终逐细胞预测、各类及 macro-F1/weighted-F1、coverage/Unknown、混淆矩阵和正确分母。GBM 按患者 bootstrap/配对检验/Holm；PTC 展示全部四患者及 exact 配对结果。历史 F1 的正类、hard-call AUC 与现代概率 AUC 不混用。

## D. Marker、Unknown 与 PTC 生物证据

对应 R1-2、R2-4/7/8、R3-C2/S5/S6；衔接 PI-02/03/04/12/13。

- 将选中库、来源研究、组织/疾病/技术类型和具体基因表与当前 terminal endpoint 对齐。补真实表达 dotplot/violin；同坐标列原作者标签、不同 reference/marker 初始标签、最终 DL、Unknown。旧分数分布图复用并补投稿版。
- 公共数据与 GBM 的选库不只看名称；保留正常/疾病/免疫/血管/不匹配库作为预声明对照。RNA-supported 与混合 RNA/protein 来源先列覆盖和来源，不把未确认的 assay 记录称 RNA-only。
- PTC Unknown 做逐样本、同一可支持细胞大类内 DE/QC、marker module、TCR 支持、doublet 分数关联；样本/患者是推断单位。GBM Unknown 已完成部分复用。
- PTC 按 S2/S3 分开给逐样本/逐患者组成、全分母与 N↔T、MT↔TU 配对展示；每类配对只有两个患者，不把八样本当八独立患者。
- 对保留的 CD4/Treg/CD8 等亚群描述补 marker 证据和已存在的 TCR/T-core 敏感性结果。计算关联不能证明独立亚型真值、因果进展、新群体或治疗靶点；相应措辞改为描述性/假设性。
- 原 TCR capture 的独立假阴性率没有 gold truth。复用“转录上清楚 T 但无 productive TCR”的代理比例和阈值敏感性，明确其依赖定义；不制造真实漏检率。
- 免疫治疗讨论补原始文献依据，与现有细胞状态相连；不承诺新湿实验/临床验证，不恢复已排除的独立 Pu 队列。

验收：图表均可追溯到具体细胞/label endpoint；统计单位、效应/区间及局限清楚；每条保留的生物结论都有对应证据或被适当限定。

## E. Batch 的方法间对照

对应 R2-6；衔接最新 PI 的多样本 PTC 要求。

- 复用两组 CCA/NONE/Harmony 和 all8/分组结果；在共享类别与每患者/样本内比较 batch mixing、生物保留及误注释。
- 补标准工具在可比输入条件下的 batch-associated error/Unknown 变化，与 DG 并列；主要复用 C 的新公平预测，不重复拟合。
- 技术批次与组织来源混杂时不能声称识别了独立技术因果效应。CCA/NONE/Harmony 是实际工作流比较，不写成纯“先后顺序”实验。
- 当前 35 项 reviewer 的关闭不要求虚构 CCA-after-UMAP。若另保留 PI 的纯校正位置主张，需另定同算法、合法空间的配对设计；本计划先用已明确的生物保留/方法间误注释回答 R2-6，并在回复中说明范围。

验收：两组各自有同端点方法比较、逐样本误差、共享类生物保留和混杂说明。

## F. 跨工具一致性与资源账本

对应 R1-1、R2-9、R3-C1/S1。

- 复用已保存逐细胞预测，补真正跨工具 pairwise agreement、全工具一致比例与一致/不一致细胞的真值表现；Unknown 与可支持类别分层。旧单工具 top-K reference 投票不当作跨工具共识。
- 汇总 ARI/NMI/FMI，严格区分 partition metric 与最终注释 metric；一致不等于正确。
- 全 grid 成本按 prepare/integration、embedding/cluster、DEG、marker scoring、DL、汇总列清：唯一计算、缓存命中、实际训练、no-op、失败重试、CPU-hours、并行 elapsed、peak memory。不能把 45 次单 workflow scaling 当完整选库成本。
- CPU 120k 曲线不重跑。CPU 实现的 GPU memory 填 N/A；方法和正文取消未经实测的 GPU 优势承诺。这满足原版本 CPU 资源报告的真实范围。
- 若确认后希望保留具体 GPU 加速主张，另加同数据/同参数的 DG CPU+GPU-MLP 模式五规模×三重复（15 次）实测；R 前处理仍是 CPU，单独列 GPU peak memory 与整流程耗时。此项是可选扩展，不默认改变原执行版本。

验收：完整搜索账本与单 workflow 曲线分开；硬件、实际资源、重复、缓存/冷启动限制可核验；无未实测 GPU 数据。

## G. 投稿稿件、图表与逐条回复

对应 R1-3、R2-7/8、R3-A、全部 M 项及已有实验待整理条目。

- 在已有 manuscript/response 的副本上修订，保留原稿。按原实现描述 context-aware marker calibration＋expression MLP；不称真实 GNN/message passing/generative model。
- 统一 HVG2000、PCA30、MLP256/128、LogNormalize、10 epochs、0.90/0.70、实际 scorer、训练划分、指标正类和 accuracy 定义；错误描述不通过改原算法来迁就。
- 主图包含完整 preprocessing→注释/终端流程；消融条件单列补充，不挤入主图。统一图号、表号、工具/引用/资源映射。
- 110,497→92,404 形成逐阶段 barcode 台账和 flowchart，分别记录来源、过滤条件、排除数量和对应输入文件，不用差值推断过滤原因。
- 输出 Supplement A（GBM/流程与消融）、Supplement B（公平方法比较）、完整附表/资源与来源、35 条 point-by-point response。每条回复绑定图表和源文件，未完成/受限项明确写原因与修订后的主张。原 Figure S1 与 Tables S1–S4 必须逐一验收：关键图可直接查看，XLSX/CSV 可打开，说明和新编号明确；不能仅以新的 Supplement A/B 替代旧请求。
- 最终文稿逐页检查：图号、清晰度、图例、指标、细胞/患者数、引用、主张范围和历史/新结果版本。

## H. 来源、可复现包与发布

对应 R3-M3、R3-T1/T2/T3；这些是单独的完成关卡。

- 湿实验真实阈值（>80/>85）、strainer（70/40 μm）、建库 5′/3′、DoubletFinder 版本不能靠新代码猜。先找现有原始实验记录；仍有冲突则列给用户确认来源。
- GEO 本地准备包已有，先核验 FASTQ/BAM/CRAM 的实际位置/可用性、GEX/TCR manifest、元数据与 cell-count 台账。只有取得 accession 和 reviewer access 并核验能用，才关闭数据公开项。
- 一键 Table 2/GBM 示例、冻结环境/依赖、输入 hashes、预期输出、故障说明和 figure/table 重建入口在指定 branch 完成；用干净环境实际执行验证后再标 reproducible。
- 准备可发布归档包及 DOI 元数据；正式 GEO/Zenodo 提交、外部发布前只在最后一步核对用户指定账号及具体发布内容。不得访问任何未确认身份的连接。仅有本地 ZIP/branch 不等于已有 DOI。
- 原 Accuracy 定义/脚本仍找不到时，历史值保留明确来源缺口，新文稿主表使用可复算且定义清楚的指标，不伪称已精确复现。

## 执行与交付

1. **确认后先冻结配置与映射（P0）。** 核验原锚点、数据/方法缺格、排除标准、SignacX 历史表、marker 候选、统计规则和缓存 hash。
2. **先完成 GBM（P1）。** A＋GBM 的 B/C/D/F；完成每个新增 condition 的图与终端结果，追加原 notebook 后再进入后续大任务。
3. **PTC 与公共基准并行补缺（P2）。** C/D/E，按独立数据/样本/条件 SLURM arrays 运行。PTC 大内存阶段严格限制并发，保存可恢复中间件。
4. **交付收尾（P3）。** G/H；代码、表图、原 notebook、Paseo 阅读导出与原 OneDrive 目录一致。来源检索和修稿草案可与计算并行，不把未完成实验写成已得结论。

计算起始资源沿用已有成功脚本：GBM/评分/终端常规 4 CPU、16–32 GB，数组先上限 8；PTC prepare 8 CPU、192 GB、并发至多 2，原 R DEG/评分先用 128 GB、4 workers。旧 batch MaxRSS 不是所有子进程内存之和，不能据此直接压低配额。正式提交前按队列配额和 pilot 的完整任务内存观测下调并发或调整内存；不因内存不足删条件、删细胞或更换 marker。所有科学计算、指标、图和 notebook 执行均走 SLURM。监测 OOM/timeout，按原配置恢复；只重复失败或输入不一致的部分。

新增结果建议独立子目录 `results/hvg_ptc_20260916_v1/reviewer_completion_20260920/`，旧结果只引用。实际运行量和预计时间在第一批 pilot 的 wall-time/MaxRSS 之后报告；不在无测量依据时承诺一晚完成。

**最终验收：35/35 条都有来源→工作→证据→回复定位。** 科学支持、反例、未能独立验证与工作完成分别记录。湿实验记录、GEO、永久 DOI 未拿到时，清楚保留这些未闭环项，不能用“计算全部跑完”替代。

本计划等待用户确认后执行。默认不新增独立 Pu、不重新调 Darmanis、不改原算法成 GNN、不额外启动 GPU 优势项目；保留用户已确认的 SignacX 历史 0T。
