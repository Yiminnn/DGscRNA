# Darmanis 的用途修订

用户最新意见：Darmanis 原始标签与目标任务不匹配，不适合作为证明既有 paper 流程最优的主样本。接受这个任务范围修订：**不再围绕 Darmanis 新增大规模 HVG/聚类调优；主验证转到 PTC + GSE274546。** 先前将 PI 早期点名直接转化为今天的主实验优先级，并不合适。

## 已核实的问题

1. 本地策展文件 `data_bench/brain_GBM/metadata.json` 只有 7 个细胞大类：myeloid、neoplastic、OPC、oligodendrocyte、astrocyte、vascular lymphangioblast、neuron。不能把它当作具有完整免疫亚型、肿瘤细胞状态标签的数据集。
2. `handoff/mc/gbm_labels.py` 已记录旧 12 簇映射的限制：部分肿瘤簇高度患者特异，C04 的 125 个细胞全部来自 BT_S6；通用 marker 词表并不表达这个患者特异状态。因此细簇 ID 与 marker cell-type 名称无法自然一一对应。
3. `handoff/b1/gbm_panel_meta.json` 和 `handoff/b1/B1_B3_REPORT.md` 已记录当时 GBM-only marker 库缺少 OPC，而数据中有 406 个 OPC。这个结构性缺类是 **marker/reference 的覆盖缺口**，不能归因于样本测序差，也不能靠改变聚类算法解决。该记录适用于当时的库，不自动推广到其他 marker 库。
4. 旧映射源码还记录了自身问题：`non-neuron` 被正则归为 neuron，Neutrophil 被归为 Other。当前 campaign 使用另一套冻结映射；旧脚本问题不能自动归因于当前所有结果，更不能归咎于原论文全部标签错误。
5. 原研究来自 4 位患者的肿瘤核心和周围组织；Periphery 不等于健康对照。原研究关注侵袭肿瘤细胞与空间相关异质性，原论文对粗类注释给出了表达与参考比较证据。现有读到的资料没有证实“所有原始表达数据不可用”。来源：[Darmanis et al., 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5810554/)、[GSE84465](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84465)。

这次只读取原始元数据、已完成的 QC/标签审计和原论文，没有在登录节点加载矩阵或做数值计算。没有新增 QC 实验，因此不以未经检验的 dropout、线粒体比例或测序深度数值作质量判决。

## 后续如何处理

- 移出新增主实验和主图案例，不以其细粒度排名判断原始 HVG/UMAP/HDBSCAN 是否最优。
- 既有完整结果、原始数据及失败/限制记录全部保留。若报告其 7 大类结果，统一语义规则与 coverage 分母，并明确其有限任务范围。
- 不因 DG-scRNA 在某数据集得分低而删除该数据集的既有 benchmark，也不将 marker 库缺类一般化为数据不可用。
- GSE274546 作为当前已在使用的 GBM 主队列，先核验原作者标签层级、cell ID 对应、患者与主要类别支持，再补原 R 工作流。沿用既定主/完整队列，不按新模型排名选样本。
- PTC 保持既有分组、原 marker roster、原始基线与终端 DL 口径。Reviewer 多类证据继续使用已经纳入的其他公开队列，不依赖 Darmanis 一套数据。
