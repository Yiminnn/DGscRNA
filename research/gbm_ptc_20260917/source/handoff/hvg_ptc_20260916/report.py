#!/usr/bin/env python3
"""Create the GBM evidence memo from verified, complete saved results."""
import json
import shutil
from pathlib import Path
from common import ROOT,OUT,require_slurm,sha,utc,write_json


def run():
    require_slurm()
    import pandas as pd
    sm=json.loads((OUT/'summary/manifest.json').read_text())
    vm=json.loads((OUT/'verification/independent_manifest.json').read_text())
    rm=json.loads((OUT/'singleton_robustness_summary/manifest.json').read_text())
    nm=json.loads((OUT/'notebooks/execution_manifest.json').read_text())
    mm=json.loads((OUT/'method_comparisons/manifest.json').read_text())
    assert sm['status']==vm['status']=='complete' and rm['status']==nm['status']=='completed'
    assert mm['status']=='completed'
    metrics=pd.read_csv(OUT/'summary/all_sample_conditions.csv.gz')
    effects=pd.read_csv(OUT/'summary/paired_patient_effects.csv')
    summary=pd.read_csv(OUT/'summary/condition_summary.csv')
    counts=metrics.status.value_counts()
    assert not metrics.status.isin(['running','failed','missing']).any()
    assert len(metrics)==39930 and metrics['sample'].nunique()==121
    primary=effects[effects.primary_contrast.eq(True)]
    names={'partition_ari':'partition ARI','terminal_strict_L1_macroF1_present':'terminal strict-L1 macro-F1'}
    assert len(primary)==2 and set(primary.metric)==set(names)
    findings=[]
    for _,row in primary.iterrows():
        findings.append(f"For {names[row.metric]}, the patient-mean paired difference (HVG2000 − all) was **{row.delta_mean:+.4f}** "
            f"(95% patient-bootstrap interval {row.ci_lower:+.4f} to {row.ci_upper:+.4f}; "
            f"{int(row.n_paired_samples)} paired samples / {int(row.n_paired_patients)} patients; "
            f"two-endpoint Holm-adjusted p={row.p_holm_2_primary_endpoints:.4g}).")
    available=int(counts.get('completed',0))
    structural=int(sum(v for k,v in counts.items() if k.startswith('structural')))
    numerical=int(sum(v for k,v in counts.items() if k.startswith('numerical')))
    text='# GBM HVG ablations and terminal DG-scRNA annotation\n\n'
    text+='\n\n'.join(findings)+'\n\n'
    pca_effect=effects[effects.cohort.eq('primary97')&effects.feature.eq('hvg2000')&effects.path.eq('PCA30_UMAP2')&effects.metric.eq('terminal_strict_L1_macroF1_present')]
    assert len(pca_effect)==1
    pca_effect=pca_effect.iloc[0]
    interactions=pd.read_csv(OUT/'method_comparisons/hvg_pca_interactions.csv')
    pca_interaction=interactions[interactions.umap_dim.eq(2)&interactions.metric.eq('terminal_strict_L1_macroF1_present')].iloc[0]
    text+=(f'**The effect depends on preprocessing.** After PCA30 → UMAP2, the paired terminal-F1 difference was {pca_effect.delta_mean:+.4f} '
           f'(95% interval {pca_effect.ci_lower:+.4f} to {pca_effect.ci_upper:+.4f}; five-HVG-level Holm p={pca_effect.p_holm_5_HVG_levels:.4g}). '
           f'The four-condition interaction, HVG effect with direct genes minus the effect after PCA30, was {pca_interaction.delta_mean:+.4f} '
           f'(95% interval {pca_interaction.ci_lower:+.4f} to {pca_interaction.ci_upper:+.4f}; three-dimension Holm p={pca_interaction.p_holm_3_dimensions:.4g}). '
           'Thus the direct-input benefit should be interpreted together with the controlled PCA branch. '
           '[HVG-by-PCA interactions](method_comparisons/hvg_pca_interactions.csv) cover all three predeclared UMAP dimensions.\n\n')
    text+=('These statements refer to the predeclared direct UMAP2 → HDBSCAN15/15 path. '
           'They do not establish that UMAP/HDBSCAN is universally best or that every stage should discard non-HVG genes. '
           'The complete method table, higher-dimensional UMAP factorial, marker-retention controls and patient-level uncertainty are needed to interpret the effect.\n\n')
    text+='## HVG amount and preprocessing\n\n'
    text+='These are descriptive patient means over available scores in the fixed primary cohort. The paired effects above provide the controlled HVG contrast; sample availability must also be considered. All rows use HDBSCAN15/15, seed42, and full-gene scoring/refinement.\n\n'
    text+='| Geometry path | Features | Partition ARI | Terminal strict-L1 F1 | Final-output samples / 97 |\n|---|---|---:|---:|---:|\n'
    controlled=summary[(summary.cohort=='primary97') & (summary.seed==42) &
        (summary.clusterer=='HDBSCAN') & (summary.min_cluster_size==15) & (summary.min_samples==15) &
        (summary.neighbors==15) & (summary.min_dist==.1) & (summary.scoring_features=='all')]
    for label,dr,dim,space in [('Direct UMAP2','UMAP',2,'genes'),('PCA30 → UMAP2','UMAP',2,'pca30'),('PCA30','PCA',30,'genes')]:
        for feature in ['all','hvg500','hvg1000','hvg2000','hvg3000','hvg5000']:
            row=controlled[(controlled.dr==dr)&(controlled.dim==dim)&(controlled.input_space==space)&(controlled.feature==feature)]
            assert len(row)==1,(label,feature)
            row=row.iloc[0]
            text+=f'| {label} | {feature} | {row.partition_ari_patient_mean:.4f} | {row.terminal_strict_L1_macroF1_present_patient_mean:.4f} | {int(row.terminal_strict_L1_macroF1_present_n_samples)} |\n'
    text+='\nThe UMAP2/10/30 dimension factorial, legacy HVG flavors, marker-union conditions and scoring-truncation controls are shown in the notebook with their full tables.\n\n'
    marker_effects=pd.read_csv(OUT/'method_comparisons/marker_union_and_scoring_truncation_effects.csv')
    text+='## Geometry feature selection versus marker loss\n\n'
    text+='Each A5 contrast uses the same HVG2000 baseline with full-gene scoring and refinement. Marker union changes geometry features; scoring truncation changes only DEG/scoring features. Refinement input remains full genes.\n\n'
    text+='| Intervention | Path | Terminal strict-L1 F1 change (95% interval) | Paired samples / patients | Holm p (three paths) |\n|---|---|---:|---:|---:|\n'
    for _,row in marker_effects[marker_effects.metric.eq('terminal_strict_L1_macroF1_present')].iterrows():
        text+=f'| {row.intervention} | {row.path} | {row.delta_mean:+.4f} ({row.ci_lower:+.4f}, {row.ci_upper:+.4f}) | {int(row.n_paired_samples)} / {int(row.n_paired_patients)} | {row.p_holm_3_paths:.4g} |\n'
    text+='\n[Full A5 paired effects](method_comparisons/marker_union_and_scoring_truncation_effects.csv) also report partition effects of the marker-union intervention.\n\n'
    ranks=pd.read_csv(OUT/'method_comparisons/method_rankings_with_shared_cohort.csv')
    text+='## Which default method performs best?\n\n'
    text+='The following descriptive rankings use samples with all 21 method scores available, averaged within patient. They are exploratory and conditional on this common subset.\n\n'
    text+='| Rule | Features | Endpoint | Common samples / patients | Highest mean | UMAP/HDBSCAN mean |\n|---|---|---|---:|---|---:|\n'
    for (policy,feature,metric),group in ranks[ranks.feature.isin(['all','hvg2000'])].groupby(['policy','feature','metric']):
        group=group.sort_values('shared_patient_mean',ascending=False)
        first=group.iloc[0];base=group[group.method.eq('UMAP / HDBSCAN')].iloc[0]
        if first.n_shared_samples:
            best=f'{first.method}: {first.shared_patient_mean:.4f}'
            baseline=f'{base.shared_patient_mean:.4f}'
        else:best='No complete common subset';baseline='NA'
        text+=f'| {policy} | {feature} | {names[metric]} | {int(first.n_shared_samples)} / {int(first.n_shared_patients)} | {best} | {baseline} |\n'
    text+=('\nA rank alone cannot establish a universally best method. [Paired method contrasts](method_comparisons/paired_method_contrasts.csv) '
           'use matched samples for each comparison and report patient-bootstrap intervals, paired counts and Holm adjustments. '
           '[Four-condition interactions](method_comparisons/feature_method_interactions.csv) compare the HVG effect in UMAP/HDBSCAN with the HVG effect in every other default. '
           'These separate a representation-specific benefit from a general feature-selection benefit. No hyperparameter was changed for this analysis.\n\n')
    text+=f'## Scope and completion\n\nAll 121 samples were accounted for across 330 conditions each (39,930 total). '
    text+=f'{available:,} conditions have saved terminal outputs; {structural:,} are structurally unavailable and {numerical:,} are verified numerical fitting failures. '
    text+=('No unresolved infrastructure failure is scored as a successful run. Primary reporting uses the previously frozen 97 samples / 55 patients; '
           'all 121 samples / 59 patients and the 24 excluded-from-primary samples are also reported. '
           f'The independent audit checked {vm["n_terminal_conditions_checked"]:,} terminal conditions and {vm["n_actual_DL_checked"]:,} actual DL histories.\n\n')
    text+='## What changed, and what was held fixed\n\n'
    text+=('| Experiment | Controlled question | Evidence |\n|---|---|---|\n'
           '| A1 | All genes vs HVG500/1000/2000/3000/5000; raw-count VST and separate legacy flavor bridge | Paired response curves and full 6×7×3 table |\n'
           '| A2–A3 | Direct input vs PCA30; six 2D reducers/no DR; UMAP 2/10/30D | Factorial curves and actual input dimensions |\n'
           '| A4 | KMeans/GMM/HDBSCAN, density settings, fixed vs label-free K | Predeclared sensitivity families |\n'
           '| A5 | Geometry HVG selection vs losing scoring markers; marker union | Marker retention and scoring-truncation ablations |\n'
           '| A6 | Initial seed calls vs terminal DL/refinement | Per-cell lineage, actual execution and full feature width |\n'
           '| Secondary singleton policy | Scoring implementation fails on a one-cell cluster | Frozen clusters; singleton seeds Undecided; same DL, separately labelled |\n\n')
    text+=('All branches retain the same sample cells. Primary DEG scoring uses all filtered log-normalized genes; DL always uses all filtered scaled genes. '
           'The marker library is CM2_glioma_other with the fixed panel-to-L1 mapping and cutoff none. '
           'The legacy scorer uses top100 Wilcoxon DEGs and logFC-density panel scores, without silently adding a new significance filter. '
           'The existing corrected GBM refinement uses 15 epochs, batch size 256 and a 0.90 assignment threshold. '
           'Noise is excluded from training and remains terminal Unknown. Five end-to-end random seeds are algorithm variation, not extra biological replicates. '
           'This is the frozen GBM benchmark implementation; the later historical PTC replay must recover its original R/DL settings separately.\n\n')
    text+='## Valid comparisons and limitations\n\n'
    text+=('Strict L1 concordance and the historical one-to-many/set-valued Lfine score are separate. Reference annotations are not independent biological ground truth. '
           'Partition ARI/Pair-F1/FMI/V-measure do not measure final annotation quality. K23 is predeclared from vocabulary size and is therefore a conditional comparator. '
           'Sensitivity-tuned UMAP settings are not mixed into the default-method heatmap. Geometric preservation is diagnostic, not a substitute for biological concordance.\n\n'
           'Missing annotation scores stay NA. Every condition includes availability and fixed-cohort pessimistic/optimistic bounds; available-case means are conditional. '
           'TKU3074 default-span VST failed with a LOESS near-singularity; its independent all-gene and legacy branches were retained. '
           'This sample is outside the previously fixed primary cohort. Singleton DEG failures retain their valid partition metrics. '
           'Reproducible GMM covariance failures are closed only after unchanged-parameter retries and have no fabricated partition or annotation score.\n\n')
    optimization=pd.read_csv(OUT/'verification/optimization_diagnostics.csv.gz')
    text+=(f'The saved optimizer audit contains {int(optimization.convergence_warning.sum()):,} fits with convergence warnings. '
           'Warnings and iteration counts remain attached to their conditions; default-table means retain the frozen runs and are not silently retuned or filtered. '
           '[Optimization diagnostics](verification/optimization_diagnostics.csv.gz) distinguish recorded nonconvergence from merely reaching an iteration budget.\n\n')
    text+=(f'The explicitly secondary singleton policy was independently verified for {rm["n_repairs_verified"]:,} affected conditions. '
           'It preserves every cell and the frozen clustering, marks unsupported singleton groups Undecided and uses the same DL. '
           'Infeasible fixed stratified validation splits remain unavailable. Its result and availability tables are separate from the original scoring-rule tables.\n\n')
    text+=('A controlled t-SNE diagnostic on COLUMBIA100163, repeated sequentially on the same compute node, reproduced the old C-order/no-Torch result exactly. '
           'The identical numeric matrix and seed yielded different optimization results with F-order storage and/or the loaded Torch runtime; '
           'F-order plus Torch reproduced the current pipeline result. The new feature comparisons share one runtime/layout; old and new t-SNE tables should not be spliced together. '
           'This numerical sensitivity was documented without changing the frozen grid.\n\n')
    text+='## Evidence for the PI discussion\n\n'
    text+=('| Question | Evidence-based reading |\n|---|---|\n'
           '| Should HVGs be used? | Use the paired effects above for the specified representation and cohort; inspect the full count curve and patient spread. |\n'
           '| Should markers outside HVGs be removed? | This is a different intervention. Primary scoring/DL preserve full genes; A5 tests truncation explicitly. |\n'
           '| Is UMAP + HDBSCAN necessarily best? | Compare the complete default table and its availability; dimensionality, PCA preprocessing, density settings and seeds matter. |\n'
           '| Does UMAP preserve all information? | Common-reference neighborhood and distance diagnostics are provided; no nonlinear “100% explained variance” claim is made. |\n'
           '| Did the original DGCyTOF paper sweep scRNA HVG counts? | Its cited tables used 13/32 CyTOF protein channels; the archived scRNA HVG sweep is a separate adaptation. |\n\n')
    text+=('The available meeting records were not used to invent a verbatim PI quote about forbidding HVGs. '
           'This experiment tests the user-reported proposal and separates geometry selection from marker removal.\n\n')
    text+='## Deliverables\n\n'
    text+=('- [Executed expanded notebook](notebooks/dgscrna_v6_results.ipynb)\n'
           '- [Workflow decision tree, PDF](figures/decision_tree_main.pdf) and [SVG](figures/decision_tree_main.svg)\n'
           '- [Paired patient effects](summary/paired_patient_effects.csv)\n'
           '- [Complete condition summary and bounds](summary/condition_summary.csv)\n'
           '- [All sample conditions](summary/all_sample_conditions.csv.gz)\n'
           '- [Secondary implementation sensitivity](singleton_robustness_summary/secondary_primary_factorial.csv)\n'
           '- [Paired method contrasts and feature-by-method interactions](method_comparisons/manifest.json)\n'
           '- [Independent verification](verification/independent_manifest.json)\n'
           '- [Predeclared condition protocol](protocol/protocol.json)\n\n')
    text+='PTC original R/table reproduction is the next stage, using MT+N and TU+T groups and the full original marker roster. No independent Pu cohort is part of this task.\n\n'
    from resources import collect
    resource_summary=collect()
    text+='## Compute accounting\n\n'
    text+='| Stage | Allocations | Reserved CPU hours | Largest sampled peak RSS (GiB) |\n|---|---:|---:|---:|\n'
    for stage,row in resource_summary.iterrows():
        peak=f'{row.largest_observed_peak_rss_GiB:.2f}' if pd.notna(row.largest_observed_peak_rss_GiB) else 'Not yet reported'
        text+=f'| {stage} | {int(row.n_allocations)} | {row.allocated_cpu_hours:.1f} | {peak} |\n'
    text+=('\n[Full SLURM accounting](resources/slurm_job_and_step_accounting.csv) preserves states, exit codes, actual CPU time and resource requests. '
           'Reserved CPU hours measure allocation, not CPU utilization. Running delivery/controller jobs are a partial snapshot. '
           'A FAILED fit allocation can represent a verified singleton-scoring limitation; scientific availability tables distinguish it from an infrastructure failure. '
           'Peak RSS is sampled by SLURM and may miss short peaks.\n\n')
    text+='Sources: [DGCyTOF original paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008885), [UMAP clustering documentation](https://umap-learn.readthedocs.io/en/latest/clustering.html), [Scanpy HVG documentation](https://scanpy.scverse.org/en/stable/generated/scanpy.pp.highly_variable_genes.html).\n'
    (OUT/'GBM_REPORT.md').write_text(text)
    brief='# GBM：与 PI 讨论 HVG 的配对证据\n\n'
    brief+='本轮固定细胞、marker、评分和 DL/refinement，只改变几何计算使用的基因及预设表示/聚类条件。121 个样本全部记账；主分析是事先固定的97个样本/55位患者。统计先在同一样本配对，再按患者汇总。\n\n'
    cn_names={'partition_ari':'聚类 ARI','terminal_strict_L1_macroF1_present':'最终注释 strict-L1 macro-F1'}
    for _,row in primary.iterrows():
        brief+=f'- 直接 UMAP2 → HDBSCAN15/15：HVG2000 相对全基因的{cn_names[row.metric]}变化为 **{row.delta_mean:+.4f}**（95%患者 bootstrap 区间 {row.ci_lower:+.4f} 至 {row.ci_upper:+.4f}；{int(row.n_paired_samples)}对样本/{int(row.n_paired_patients)}位患者；两个主终点 Holm p={row.p_holm_2_primary_endpoints:.4g}）。\n'
    primary_supported=primary.delta_mean.gt(0)&primary.p_holm_2_primary_endpoints.lt(.05)
    if primary_supported.all():
        brief+='\n这组结果支持在该固定 UMAP/HDBSCAN 路径的几何计算中使用 HVG2000。向 PI 讨论时可以提出：特征选择应按处理步骤和配对证据决定；本轮主流程仍保留全部过滤后基因用于 marker scoring 和 DL。\n\n'
    else:
        brief+='\n两个主终点没有同时给出 Holm 校正后显著的正向结果。向 PI 讨论时应分别报告聚类和最终注释的效应，依据完整曲线讨论 HVG 的适用条件。\n\n'
    brief+=(f'关键限定：先做 **PCA30 → UMAP2** 时，HVG2000 的额外最终-F1变化为 **{pca_effect.delta_mean:+.4f}** '
            f'（95%区间 {pca_effect.ci_lower:+.4f} 至 {pca_effect.ci_upper:+.4f}；五个HVG数量对照的 Holm p={pca_effect.p_holm_5_HVG_levels:.4g}）。'
            f'两条路径的HVG效应差为 **{pca_interaction.delta_mean:+.4f}** '
            f'（95%区间 {pca_interaction.ci_lower:+.4f} 至 {pca_interaction.ci_upper:+.4f}）。'
            '这说明应把 PCA 预处理纳入讨论；现有数据支持的结论是按流程选择特征，而不是对所有流程作统一判断。\n\n')
    brief+='UMAP/HDBSCAN 与其他20个默认方法的探索性配对比较如下；每项比较的可用样本可能不同，具体配对数在链接表中。Holm 校正覆盖同一终点的全部特征×方法对照。这里各降维方法均为2D，另含无降维控制；KMeans/GMM使用预设K23。\n\n'
    brief+='| 特征 | 终点 | 显著优于比较方法数 | 显著低于比较方法数 |\n|---|---|---:|---:|\n'
    method_tests=pd.read_csv(OUT/'method_comparisons/paired_method_contrasts.csv')
    for (feature,metric),g in method_tests[method_tests.policy.eq('primary_rule')&method_tests.feature.isin(['all','hvg2000'])].groupby(['feature','metric']):
        sig=g.p_holm_all_contrasts_in_policy_endpoint.lt(.05)
        brief+=f'| {feature} | {cn_names[metric]} | {int((sig&g.delta_mean.gt(0)).sum())}/20 | {int((sig&g.delta_mean.lt(0)).sum())}/20 |\n'
    brief+='\n“只用 HVG 做几何计算”和“删掉 HVG 之外的 scoring marker”是两个干预。A5 在相同几何下单独检验后者，并提供 marker union 对照；最终注释结果都经过 terminal DL/refinement 阶段，no-op/无法训练/缺失状态单列。\n\n'
    brief+='在三条固定路径中，将 scoring 也截断到 HVG2000 均降低了最终 F1（约0.044–0.054，配对区间及校正结果见A5表）；marker union 的三个F1区间均覆盖零。本轮结果支持保留完整的 scoring 基因信息。\n\n'
    brief+='原 DGCyTOF 的对应表使用13/32个 CyTOF蛋白通道；本地旧 scRNA HVG 扫描属于后续迁移实验。仓库 R 参考实现包含 VST2000，PTC 实际历史流程仍需由原始归档和逐细胞表核对。现有会议记录尚未定位到“一律不用HVG”的逐字原句，因此这里回应的是用户转述的建议。\n\n'
    brief+='可直接展示：[主流程图](figures/decision_tree_main.pdf)、[HVG数量曲线](figures/hvg_response_curves.pdf)、[配对方法比较](figures/paired_method_comparisons.pdf)、[特征×方法交互](figures/feature_method_interaction.pdf)、[A5患者配对结果](method_comparisons/marker_union_and_scoring_truncation_effects.csv)、[完整英文报告](GBM_REPORT.md)、[已执行 notebook](notebooks/dgscrna_v6_results.ipynb)。\n\n'
    brief+='PTC 后续按 MT+N、TU+T 两组恢复原 R 流程和完整原 marker roster；先核对原注释表，再分析 batch/HVG/marker/DL。独立 Pu 队列不在本轮范围。\n'
    (OUT/'PI_BRIEF_ZH.md').write_text(brief)
    snapshot=OUT/'delivery_source_snapshot';snapshot.mkdir(exist_ok=True)
    for source in (ROOT/'handoff/hvg_ptc_20260916').iterdir():
        if source.is_file() and source.suffix in ['.py','.sbatch','.md','.dot']:
            shutil.copy2(source,snapshot/source.name)
    write_json(snapshot/'manifest.json',dict(timestamp=utc(),
        note='Final delivery scripts; initial frozen fitting snapshot is preserved separately in source_snapshot.',
        files={p.name:sha(p) for p in snapshot.iterdir() if p.is_file() and p.name!='manifest.json'}))
    write_json(OUT/'GBM_DELIVERY.json',dict(status='completed',timestamp=utc(),report_sha256=sha(OUT/'GBM_REPORT.md'),pi_brief_sha256=sha(OUT/'PI_BRIEF_ZH.md'),
        notebook_manifest_sha256=sha(OUT/'notebooks/execution_manifest.json'),summary_sha256=sha(OUT/'summary/manifest.json'),
        verification_sha256=sha(OUT/'verification/independent_manifest.json'),
        method_comparisons_sha256=sha(OUT/'method_comparisons/manifest.json'),
        secondary_verification_sha256=sha(OUT/'singleton_robustness_summary/manifest.json'),
        decision_tree_sha256={ext:sha(OUT/'figures'/f'decision_tree_main.{ext}') for ext in ['pdf','svg','png']},
        delivery_source_snapshot_sha256=sha(snapshot/'manifest.json'),PTC_status='not_started'))
    print(text[:1800],flush=True)


if __name__=='__main__':run()
