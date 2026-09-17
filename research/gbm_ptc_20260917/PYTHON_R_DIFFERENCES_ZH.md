# 当前 Python 与原始 R workflow：差异、原因和解决方式

检查基点：Python 包 `align-r-reference` 的 `4bf17c4`；仓库 `examples/R/`；归档 PTC 的 R/source.R、R/source.py；本轮实际执行的 GBM 和 PTC 重建脚本。代码审查日期：2026-09-17。

**结论：简易 Python 包与论文原始流程的差别较大，足以改变初始训练标签和最终注释。主要差别是算法、数据层、默认参数和实现缺口，不能笼统称为 R/Python 平台固有差异。可以解决其中的大部分；严格复现 Seurat CCA 时可以由 Python 调用固定版本的 R 后端。** 本审查没有把这些差别都重新跑成患者级消融，因此不能给出“整体差了多少 F1”的因果数字。

## 先区分三份代码

| 实现 | 实际含义 |
|---|---|
| 仓库 `dgscrna/` 简易 Python 包 | Scanpy 预处理和 DEG，Python 聚类，带 dropout 的 PyTorch DL，内部伪标签分数选择输出 |
| 当前 GBM 实验执行器 | 显式保存 lognorm；VST/HVG、PCA/UMAP等网格；单独控制 scoring/DL 特征、噪声和失败状态；DL沿用包实现并设15epochs等参数 |
| 当前 PTC 原始基线重建 | 读取原始八样本 checkpoint/分区；R原始风格 DEG 和 density score；匹配原 source.py 的 PyTorch DL；10epochs、CCA2000输入等 |

“PTC 重建与原函数的权重/概率完全一致”指第三行在同输入、同seed、同当前环境下的核验，不代表第一行的整条流程与 R 已一致。原始所谓“R流程”的 DL 本来就由 Python/PyTorch 实现。

## 逐环节对照

| 环节 | 原始 R/PTC 路线 | 当前简易 Python 包 | 性质与解决方式 |
|---|---|---|---|
| QC / 双细胞 | 显式计算 percent.mt，过滤，回归mt，DoubletFinder；历史 checkpoint 固定保留细胞 | mt列存在才过滤；函数没有建立mt基因标记；无对应DoubletFinder或mt回归 | 方法/实现差别；先固定相同保留细胞，显式匹配QC规则。换Scrublet也属于换算法 |
| normalization / 表达层 | RNA、integrated data、scale.data各有用途 | normalize_total1e4、log1p后直接scale `.X`，没有保存lognorm层 | 可修实现问题；保存 counts/lognorm/scaled，逐阶段显式指定输入 |
| HVG | 每样本VST2k；SelectIntegrationFeatures；最终CCA checkpoint 2k | dispersion型 `flavor='seurat'` 默认阈值，没有固定2k；DL又按HVG标记选列 | 方法/参数差别；用counts上的VST风格、固定数量和批次选择规则，再核对基因ID列表。不是Python不能取同样HVG |
| batch correction | CCA anchors + IntegrateData，产生校正表达assay | 辅助函数提供Harmony/BBKNN；主pipeline不自动调用；文档列出的scVI没有实现分支 | 算法/能力差别；CCA可调用同一R后端。Harmony/BBKNN应作为明确替代方案，不能宣称等价 |
| PCA/UMAP | archived integrated2k，PCA30；历史uwot UMAP为neighbors30、cosine、min.dist0.3 | Scanpy/Python UMAP；预处理neighbors15，UMAP未显式匹配原参数 | 先解决可调参数/输入差异，再处理不同UMAP/近邻后端。相同seed数值不保证不同实现同一坐标 |
| graph clustering | Seurat SNN图 + FindClusters；也有PCA与UMAP两条路线 | Scanpy邻接图 + Leiden/Louvain | 图构建和社区算法选择不同；相同resolution不意味着相同分区。可导入同一图或调用同后端 |
| HDBSCAN | R dbscan，minPts50；0为noise且原评分路径仍处理该编号 | Python hdbscan默认50/50；noise转Unknown并排除评分。GBM网格常用15/15，另有敏感性臂 | 参数、noise政策和实现细节都需匹配；先冻结同一embedding比较partition和noise，不能直接比cluster编号 |
| DEG | 原始默认assay为integrated；FindAllMarkers的v4规则、表达比例预筛、logFC及p值定义 | Scanpy Wilcoxon默认仅top100；优先lognorm，否则退回X；输出/校正及logFC定义不同 | 统计定义和输入差别，可解决；需要显式匹配所有筛选及校正，或复用同一R DEG表 |
| density score | 交集DEG的log2FC总和 / 完整panel长度；单基因×0.8；并列Undecided；none/mean/0.5 | 核心公式、并列和cutoff已基本对齐；额外默认padj<0.05，输入DEG不同 | 公式差距较小，但上游输入会放大差异。R density_score自身仅过滤logFC>1；FindAllMarkers另有返回筛选，不能误说R完全没有p值筛选 |
| DL architecture/loss | 256/128、LeakyReLU、无dropout；Softmax输出后再传CrossEntropy | 同宽度，但每隐藏层Dropout0.3；raw logits传CrossEntropy | 实质性模型/目标函数差异，两者都是PyTorch，完全不是R/Python固有差异。legacy复现与现代改进应显式分模式 |
| DL训练 | torch random_split90/10、固定顺序不shuffle；输出类别含Undecided；模型初始化未固定 | sklearn stratified split；shuffle=True；类别只取known；固定模型seed；可自动CUDA | 可解决的协议差别。还需匹配输入表达值、基因顺序和训练池，不能只匹配网络宽度 |
| confidence | 概率先四位小数round，再与0.90比较 | 原始float概率直接与0.90比较 | 可精确对齐；只解释边界个案，不能解释全部TTU差异 |
| 候选选择/评价 | 原notebook对NMT/TTU选择不同marker/partition/cutoff；论文TCR指标有独立定义 | `optimal_annotation`按known伪标签上的train/test weighted-F1均值选候选 | 评价目标不同。训练一致性高不等于论文TCR表现或生物学准确率高；需要明确选择目标和评价数据 |

核心源码：[预处理/整合](../../dgscrna/core/preprocessing.py)、[聚类/DEG](../../dgscrna/core/clustering.py)、[评分](../../dgscrna/core/marker_scoring.py)、[DL训练](../../dgscrna/core/deep_learning.py)、[模型](../../dgscrna/models/deep_model.py)、[pipeline](../../dgscrna/core/utils.py)、[归档R](reference/ptc_archive/source.R)、[归档DL](reference/ptc_archive/source.py)。

仓库 examples/R 与归档源码并非同哈希。R差异包括DoubletFinder接口名、density_score参数/默认clusterings；归档DL增加h5ad输入支持等。实际方法应结合原notebook选择、保存对象command log与显式调用参数判断，不能只凭某个草稿函数的默认值。归档R封装还引用全局 `panc8` 等状态，所以参考模式也应固定实际输入而不是不加检查地执行整份历史脚本。

## batch effect：哪些不同可以消除

**CCA、Harmony和BBKNN校正的是不同对象。** Seurat IntegrateData生成校正表达assay；Harmony调整PCA坐标；BBKNN构建跨批次平衡的邻接图。这决定了后续聚类、marker评分和DL会不会接收到校正后的值。[Seurat IntegrateData](https://satijalab.org/seurat/reference/integratedata)、[Scanpy Harmony](https://scanpy.readthedocs.io/en/stable/generated/scanpy.external.pp.harmony_integrate.html)、[BBKNN官方实现](https://github.com/Teichlab/bbknn)。

- 要复现原始CCA流程：Python通过Rscript执行固定Seurat版本的VST、anchors、IntegrateData，导出明确的cell/gene索引、integrated表达、PCA/UMAP/分区，再回到Python DL。跨语言传递本身可以严格核对；本项目已验证原loom与重训矩阵全部值/顺序一致。
- 要比较batch算法：固定细胞、特征及评价定义。若只比较几何校正的作用，固定marker scoring与DL的表达输入，只替换产生图/embedding的步骤；若CCA同时改变评分和DL表达层，则应报告整条组合流程的差别。
- 使用Harmony后，需要显式让neighbors/UMAP/HDBSCAN消费 `X_pca_harmony` 并重建依赖的图和embedding。当前辅助函数仅写入该坐标，包后续HDBSCAN仍读取 `X_pca`/现有 `X_umap`；这属于可修代码连接问题。
- 使用BBKNN时，Leiden应读取生成的图；HDBSCAN本身使用坐标，不能仅因BBKNN图存在就认定HDBSCAN已经得到批次校正。

HVG也有可用的对应实现：Scanpy的`seurat_v3`/`seurat_v3_paper`使用counts，后者配合batch_key模拟SelectIntegrationFeatures的跨批次排序；仍应实际核对选中特征和版本，而不是只匹配“2000”这个数字。[官方源码说明](https://github.com/scverse/scanpy/blob/main/src/scanpy/preprocessing/_highly_variable_genes.py)。

## 已通过SLURM确认的实现问题

小型合成数据和受控调用路径检查保存在 [probe结果](validation/package_probe_results.json)。其中两个probe以mock隔离数据流，不衡量Harmony或聚类性能，也不使用患者数据。

1. **默认预处理 → DEG层丢失：** 在合成数据上，lognorm没有保存，mt指标没有建立，返回的200个logFC全部非有限值。已有单元测试的toy对象事先手动保存lognorm，不能覆盖该完整入口问题。
2. **多个聚类的DEG覆盖：** `utils.py`先给每种聚类调用find_markers，但都写`rank_genes_groups`；之后才逐种评分。两种聚类的probe中，Leiden评分读取了HDBSCAN的DEG。应使用每个partition独立的key，并把对应key传给评分。
3. **Harmony结果未传递给HDBSCAN：** probe确认`X_pca_harmony`存在，但PCA模式HDBSCAN消费的是未校正`X_pca`。
4. **Unknown可成为训练类别：** 包的training mask只排除Undecided；当噪声已被评分为Unknown时，Unknown会作为known类训练。当前GBM执行器显式排除noise/abstentions，因此不能把包入口问题直接归到已完成的GBM结果上。

这些是实现问题，可以修复。本分支先保存现状和证据，未悄悄改变包行为。

## 哪些才是平台/数值实现差异

Seurat/AnnData容器、矩阵转置、稀疏格式、基因ID和顺序属于接口适配，通常能消除。浮点精度、BLAS/线程、CPU/GPU、近邻搜索、UMAP优化与版本差异可能使独立实现难以逐位一致；可以通过固定版本、精度、线程、随机状态和共用后端缩小差异，再用阶段级指标验收。不能事先保证所有不同后端逐位相等，也不能把上面的算法替换归入这个类别。

当前最有力的实例是：原始R流程的DL实际也是Python；在同一输入和seed下，重建DL与原完整函数的权重和概率完全一致。这说明DL协议可以对齐。它仍不恢复未记录的历史训练池/模型状态，所以历史PTC逐细胞差异与平台差异不能混为一谈。

## 建议的修复顺序

1. 修数据流问题：counts/lognorm/scaled分层、partition独立DEG、显式corrected representation、noise/Unknown训练政策和终端状态。
2. 建立显式R-reference模式：固定QC细胞、四类基因集合（anchor、geometry、scoring、DL）、CCA/SNN/R-DEG和原DL协议；需要严格一致的步骤共用R后端。保留Python-native模式作为独立方法选择。
3. 逐阶段验收：细胞和基因索引 → HVG → integrated表达 → embedding/图 → partition及noise → DEG及density矩阵 → 初始训练池 → DL概率/最终标签。分区比较使用允许簇编号置换的指标。
4. 基线验收后再做CCA/Harmony/BBKNN及其他消融。使用预设的成对种子和相同下游输入；不通过挑seed、改reference标签或切换指标制造“最优”。

这里列出的是可执行的对齐路线。本次只归档当前版本并完成审查，没有将这些改动混入已经跑完的实验版本。
