"""Keep one clear entry point in the existing result directory, retaining its old index."""
from pathlib import Path
import shutil
from common import OUT,ROOT,require_slurm


def update(stage,phase):
    require_slurm();stage=Path(stage)
    assert phase in ['GBM_core','GBM_full','PTC_full']
    rel=OUT.relative_to(ROOT).as_posix()
    index=stage/'README.md';backup=stage/'README_before_paper_claim_validation_20260917.md'
    if index.exists() and not backup.exists():shutil.copy2(index,backup)
    def link(path,label,required=True):
        if required:assert (stage/path).exists(),f'Delivery index target missing: {path}'
        return f'[{label}]({path})'
    phases={'GBM_core':'GBM 原 R 主网格已完成；比较方法、机制控制和规模测试另行验收。',
            'GBM_full':'GBM 主网格、机制控制、比较方法和重复规模测试已完成；PTC 后续补实验待完成。',
            'PTC_full':'本轮 GBM 与 PTC 补充实验已完成，原 notebook 与历史结果保留。'}
    rows=[link('notebooks/dgscrna_results.ipynb','原结果 notebook：全部新内容继续追加于原文件'),
          link(rel+'/summary/GBM_CORE_REPORT_ZH.md','GBM 主网格中文结果')+' · '+link(rel+'/summary/GBM_CORE_REPORT.md','English core report'),
          link(rel+'/summary/workflow_decision_tree.pdf','英文完整流程与 ablation 决策树：PDF')+' · '+link(rel+'/summary/workflow_decision_tree.svg','可编辑 SVG'),
          link(rel+'/summary/all_clustering_figure_index.csv','GBM 全部 2,904 个聚类结果的图索引'),
          link(rel+'/summary/primary_fixed_marker_24.csv','固定 marker 的六档特征 × 四路线')+' · '+link(rel+'/summary/fixed_marker_patient_paired.csv','配对患者差与区间'),
          link(rel+'/summary/marker_context_distributions.pdf','不同 marker 来源的患者分数分布'),
          link(rel+'/summary/patient_heldout_summary.csv','训练患者选择配置后的留出评价')+' · '+link(rel+'/summary/marker_DL_factorial_summary.csv','marker × DL 四方比较')]
    if phase in ['GBM_full','PTC_full']:
        rows += [link(rel+'/GBM_full_summary/GBM_FULL_REPORT_ZH.md','GBM 完整补实验中文结果')+' · '+link(rel+'/GBM_full_summary/GBM_FULL_REPORT.md','English full GBM report'),
                 link(rel+'/workflow_choice_summary/workflow_node_evidence.csv','流程各节点能支持什么结论')+' · '+link(rel+'/workflow_choice_summary/marker_DL_patient_paired.csv','marker × DL 患者配对贡献'),
                 link(rel+'/comparison_summary/paired_patient_comparisons.csv','GBM 与全部比较方法的绝对/相对增益及配对检验'),
                 link(rel+'/controls_summary/geometry_and_DL.png','几何与 DL 特征作用')+' · '+link(rel+'/controls_summary/MLP_training_learning_curves.pdf','实际 DL 学习曲线'),
                 link(rel+'/unknown_expression/INTERPRETATION.md','Unknown 表达分析与解释范围'),
                 link(rel+'/scalability_summary/resource_repeated_mean_SD.csv','单次输入 10k–120k：三方法各重复三次的耗时与内存')]
    if phase=='PTC_full':
        rows += [link(rel+'/PTC_summary/PTC_FOLLOWUP_REPORT_ZH.md','PTC 补实验中文结果')+' · '+link(rel+'/PTC_summary/PTC_FOLLOWUP_REPORT.md','English PTC report'),
                 link(rel+'/PTC_summary/claim_summary_by_group.csv','NMT 与 TTU 的原配置和实际结果'),
                 link(rel+'/PTC_summary/all_new_clustering_figures.csv','PTC 新增 44 个聚类条件的图索引'),
                 link(rel+'/PTC_comparator_replay/ENDPOINT_CORRECTIONS.md','SignacX 与历史比较表的评价端点修正')]
    if phase!='GBM_core':
        folder='PTC_summary' if phase=='PTC_full' else 'GBM_full_summary'
        rows += [link(rel+'/'+folder+'/REVIEWER_PI_EVIDENCE.md','PI / reviewer 要求对应的证据与剩余限制')]
    receipt_folder={'GBM_core':'summary','GBM_full':'GBM_full_summary','PTC_full':'PTC_summary'}[phase]
    rows += [link(rel+'/'+receipt_folder+'/DELIVERY_RECEIPT.json','本阶段上传与下载校验记录',required=False),
             link('README_before_paper_claim_validation_20260917.md','此前 100 单元原 R campaign、reviewer 数据集及历史 PTC 结果入口'),
             link('handoff/paper_claim_validation_20260917/README.md','执行代码范围与复现说明')]
    text='# DG-scRNA 实验结果目录\n\n'+phases[phase]+'\n\n'
    text+='\n'.join('- '+row for row in rows)+'\n\n'
    text+='最终 DG-scRNA 使用 terminal DL/refinement 结果；marker-only 是消融。无需训练、无法训练及实际训练分别记录。聚类 ARI/NMI 与注释指标分开；Unknown 保留在分母。PTC 保留 NMT 的 SNN 与 TTU 的 UMAP-HDBSCAN 原分支，历史 non-T F1、严格 T/TCR 检出一致性和保存标签一致率分别解释。\n\n'
    text+='所有新图与表均位于原目录树，未新建网页或替代 notebook。逐细胞明细压缩归档的恢复位置见 [目录与归档说明](DELIVERY_LAYOUT.md)。旧结果未删除；历史模型权重及原 Accuracy 计算来源缺失的限制继续保留。\n'
    index.write_text(text)
    layout=f'''# 目录与明细归档

原 notebook 为 `notebooks/dgscrna_results.ipynb`，保留已执行输出。HPC 源码含原运行路径；在其他机器重新执行前需要配置路径和环境。

本轮目录为 `{rel}/`。报告、汇总 CSV、PNG/PDF 和英文决策树直接可读；逐细胞标签、训练记录及来源校验值另行打包。模型权重、完整概率 NPZ 和大矩阵保留在 HPC，不因交付压缩归档而删除。

| 归档 | 解压目标（相对本轮目录） | 内容 |
|---|---|---|
| `artifacts/<sample>_all_core_labels_metrics.tar.gz` | `GBM/<sample>/` | 六档特征、四路线、全部 marker/cutoff 的初始与终端标签、评价及训练记录 |
| `artifacts/<sample>_controls_and_comparators.tar.gz` | 本轮目录本身 | 比较方法、geometry-only、MLP 和表示控制明细 |
| `artifacts/single_run_resource_artifacts.tar.gz` | `scalability/` | 各方法/规模的第一轮资源测试记录 |
| `artifacts/resource_repetition_artifacts.tar.gz` | `scalability_repeats/` | 另外两次独立进程运行 |
| `PTC_followups/artifacts/<unit>.tar.gz` | `PTC_followups/` | PTC 新控制和旧网格重新评价的逐细胞明细 |
| `source/*source*.tar.gz` | 单独的源码目录 | 实际执行的固定源码快照；保留各阶段差异 |

尚未达到该阶段交付门槛的归档可能不存在；完成状态以 README 所链接的阶段报告和校验 receipt 为准。每个阶段先复制，再完整下载比对，并读回 receipt 验证。早期 receipt 记录的是当时 notebook 的散列；后续追加内容会改变 notebook，当前状态应看最新阶段 receipt。

早期 `r_reference_campaign_20260917` 的归档结构保持不变，见保留的历史 README。不要把不同真值口径、不同信息条件或不同分母的分数直接合并。
'''
    (stage/'DELIVERY_LAYOUT.md').write_text(layout)
    return ['README.md',backup.name,'DELIVERY_LAYOUT.md']
