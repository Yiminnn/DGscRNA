# 下一轮执行顺序：GSE274546 + PTC

范围：Darmanis 不再新增拟合；不新增独立 Pu 队列。保留既有结果，主实验使用原始 R 工作流和终端 DL/refinement。本文是计划，不代表新增实验已经提交或完成。

## 1. 首批：GBM 标签核验与原 R pilot

- GSE274546 沿用固定的 97 样本/55 患者主队列，121 样本/59 患者完整队列另报。
- 核对原始 counts、细胞 ID、患者 ID、原作者标签层级；固定可公平比较的细胞大类，肿瘤状态单列。marker vocabulary 缺类记录为覆盖限制。
- 固定 marker 来源、基因/蛋白证据类型、评分词表和最终指标。分开正常脑、GBM、浸润免疫、血管、AllHuman 及来源可核验的外部文献/reference markers。
- 按细胞量小/中/大选择 3 个资源 pilot 样本，选择依据不使用 DG-scRNA 分数。先跑原 HVG2000/PCA30/UMAP2/HDBSCAN50 与 PCA30/SNN 对照，均到最终 DL。
- 验收输入/细胞顺序、marker 覆盖、各阶段基因清单、完整终端状态和内存耗时。达到这些条件再扩展；不要求与旧 Python 分数一致。

## 2. 主实验：原 R 的六档特征 × 四条路线

HVG500/1000/2000/3000/5000/all × PCA30/UMAP2 × SNN/HDBSCAN，即每样本 24 个核心配置。UMAP 均在 PCA30 之后。基线参数沿用原 R；固定 marker/cutoff，保留每个患者的结果。

同时区分整条流程的特征预算和只改变几何特征的对照。后者固定评分/DL 表达及基因清单，补 2000/5000/all 的必要配对；不将 marker 丢失误归因于几何。无降维、UMAP 2/10/30D 和聚类参数作为事先列出的有限敏感性，各路线有可比的选择机会。

GBM 主报告为固定标签层级的终端 macro-F1、逐类 F1、Unknown/coverage；ARI/NMI 单列为聚类评价。

## 3. PTC 复用与模块贡献

- 复用已完成的 17 个原 marker 库 × 3 cutoff 与 HVG/聚类结果；原始八样本基线和新 NMT/TTU 整合分别保留。
- 原 NMT 的 SNN 与原 TTU 的 UMAP-HDBSCAN 均保留。近似复现足够，不继续以历史权重缺失阻塞。
- 对原配置、HVG5000/UMAP-HDBSCAN 和预先锁定的有竞争力 SNN 配置补相同的 5 个种子；明确哪些只变 MLP，哪些重跑表示/聚类。种子不作为患者重复。
- 补固定/训练患者选择 marker × 无/有 DL 的四方分析；能用既有预测、概率和训练历史完成的分析不重复拟合。
- 修复 NMT 几何对照中评分基因缺 CD3D 的机制解释：保留原臂，新增统一 marker-retention 规则的受控对照，不能据原臂全零判断几何优劣。

## 4. 发表验证：患者留出与竞争方法

- 在训练患者中选择 HVG/marker/cutoff 等配置，留出患者标签只评分；同一患者不同组织保持同一 fold。区分逐样本独立拟合与全队列无监督整合的 transductive 标签留出。
- 先对齐强 marker-based 对手及原论文对手，再复用或补齐 reference/GNN 方法。保证合适的输入、相当的 marker/reference 机会，并记录监督信息差异。
- 统一全测试细胞分母，Unknown 纳入错误并另报 coverage。旧表仅在已调用细胞上评分的结果需要重新评价。
- 输出每患者配对差和不确定区间；PTC 保留原稿 F1/AUC，并同步报告 T-positive、TCR recall 和 Unknown 相关指标。
- 已有 reviewer 公共队列结果作为多类支持；规模测试、batch、Unknown 和真实源码参数敏感性按 reviewer 缺口补齐，不先重做整个 100-unit campaign。

## 5. 交付与结论

在原 notebook 追加：固定原 R 配置与各候选的比较、HVG/marker/DL 的独立贡献、患者留出下与强对手的比较。每种聚类保留图，主图以同一显示坐标展示原标签/marker/最终 DL；英文决策树标记各消融节点。代码仅进入 `align-r-reference`，新结果整理至原授权 OneDrive 目录，不做新版网页。

所有数值分析、拟合、检验与绘图经 SLURM；按 pilot 实测内存设置并行数组、checkpoint 和 OOM/timeout 恢复。不在 pilot 前承诺完成时长。

最终将原配置第一、优化后配置第一、近似并列、组别例外分别报告；保留真实 NMT 分支。方法优越性的结论来自上述同条件和公平评价，而不是只保留最高的同标签网格分数。
