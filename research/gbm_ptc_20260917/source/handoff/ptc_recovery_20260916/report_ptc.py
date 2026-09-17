"""Write evidence-based PTC and combined PI reports from complete results."""
import json
import os
from pathlib import Path
from ptc_common import ROOT,BASE,RECOVERY,GROUPS,require_slurm,sha,utc,write_json

def run():
    require_slurm()
    import pandas as pd
    summary=json.loads((BASE/'summary/manifest.json').read_text())
    assert summary['status']=='complete'
    h=pd.read_csv(BASE/'summary/HVG_paired_patient_effects.csv')
    interaction=pd.read_csv(BASE/'summary/HVG_by_PCA_interaction.csv')
    batch=pd.read_csv(BASE/'summary/matched_batch_effects.csv')
    census=pd.read_csv(BASE/'summary/terminal_execution_census.csv')
    caches=pd.read_csv(BASE/'summary/unique_model_census.csv')
    qc=pd.read_csv(BASE/'summary/raw_counts_QC_parity.csv')
    panels=pd.read_csv(BASE/'summary/patient_heldout_panel_summary.csv')
    bio=pd.read_csv(BASE/'batch_biology/shared_lineage_metrics.csv')
    condition_means=pd.read_csv(BASE/'summary/condition_summary.csv.gz',low_memory=False)
    batch_means=condition_means[condition_means['mode'].isin(['matched','pooled'])&condition_means.analysis_scope.isin(GROUPS)&
      condition_means.scoring_assay.eq('RNA')&condition_means.library.eq('CellMarker_AllTissues')&condition_means.cutoff.eq('mean')&
      condition_means.space.eq('UMAP2')&condition_means.clusterer.eq('HDBSCAN_R')&condition_means.stage.eq('terminal_DL090')]
    assert len(batch_means)==8
    def table(frame,columns=None):
        shown=frame[columns].copy() if columns else frame.copy()
        for col in shown.columns:
            if pd.api.types.is_float_dtype(shown[col]):
                shown[col]=shown[col].map(lambda x:'NA' if pd.isna(x) else f'{x:.4f}')
        shown=shown.astype(object).where(pd.notna(shown),'NA')
        return shown.to_markdown(index=False,disable_numparse=True)
    main=h[h.path.eq('direct_UMAP2')&h.feature.eq('hvg2000')]
    name={'cluster_S2_native_ARI':'clustering ARI vs archived S2','S2_macro_F1_reference_present':'terminal broad macro-F1 vs S2',
          'productive_strict_recall':'productive TCR-positive recall','productive_strict_apparent_F1':'apparent productive-TCR binary F1',
          'productive_strict_detection_yield':'TCR detection yield among predicted T cells','coverage':'annotation coverage'}
    cols=['path','feature','metric','n_paired_samples','n_paired_patients','delta_mean','ci_lower','ci_upper','p_exact','p_holm']
    findings=[]
    for r in main.itertuples():
        if r.n_paired_samples==0:
            findings.append(f'- {name[r.metric]}: **NA**, no sample has defined values in both arms. No T calls in an arm makes detection yield undefined; this is not imputed as zero.')
            continue
        findings.append(f'- {name[r.metric]}: HVG2000 − all **{r.delta_mean:+.4f}**, descriptive patient-bootstrap interval '
          f'[{r.ci_lower:+.4f}, {r.ci_upper:+.4f}], exact paired p={r.p_exact:.3f}; {r.n_paired_samples} samples / {r.n_paired_patients} patients.')
    text='# PTC reconstruction and complete terminal DG-scRNA experiments\n\n'
    text+=('All **11,696 terminal annotation conditions** have independently verified outputs, covering '
           '**49,776 sample-stage rows**, **1,288 physical partitions**, eight samples and four patients. '
           f'There are {summary["n_unique_caches"]:,} unique terminal caches and {summary["n_unique_trained_models"]:,} distinct trained model caches; '
           'every trained model has an independent saved-weight forward check. Training reuse and no-op/no-known-label states are counted separately. Marker-only outputs are explicitly ablations.\n\n')
    text+='## What the HVG experiment shows\n\n'
    text+=('The fixed primary transfer path is direct UMAP2 → HDBSCAN15/15, AllTissues/mean, full RNA marker scoring and the same group-selected RNA2000 DL input. '
           'HVG selection changes geometry only. The effect therefore does not conflate geometry feature selection with loss of scoring genes.\n\n')
    text+='\n'.join(findings)+'\n\n'
    pca=h[h.path.eq('PCA30_UMAP2')&h.feature.eq('hvg2000')&h.metric.eq('S2_macro_F1_reference_present')].iloc[0]
    text+=(f'After **PCA30 → UMAP2**, the paired terminal macro-F1 effect was **{pca.delta_mean:+.4f}** '
           f'(descriptive interval {pca.ci_lower:+.4f} to {pca.ci_upper:+.4f}). '
           'The feature-selection response therefore depends on the preprocessing path. HVG2000 was the fixed primary contrast; '
           'the full count curve is shown rather than replacing it with whichever count scores highest afterward.\n\n')
    text+=('These are an exploratory reconstruction on four patients. There are only 16 sign-flip permutations, so the smallest two-sided p-value is 0.125. '
           'An interval excluding zero here must not be represented as a conventional significant confirmatory result. '
           'S2/S3 are prior model outputs, so ARI/F1 against them are concordance, not independent accuracy. '
           'TCR detection gives orthogonal positive evidence; non-detection is not an established negative cell identity. '
           'The p_holm column corrects available HVG levels within each path/endpoint; the two predefined ARI/F1 endpoints also have a separate correction in the CSV. '
           'NA denotes no applicable cells, for example detection yield when no T cells are called.\n\n')
    text+=table(h[h.feature.eq('hvg2000')],cols)+'\n\n'
    bygroup=pd.read_csv(BASE/'summary/HVG_by_group_patient_effects.csv')
    text+='Group-specific effects retain possible tissue/context heterogeneity instead of hiding it in the combined mean:\n\n'
    text+=table(bygroup[bygroup.feature.eq('hvg2000')&bygroup.metric.isin(['S2_macro_F1_reference_present','productive_strict_recall'])],['group']+cols)+'\n\n'
    text+='### Does PCA change the HVG effect?\n\n'
    text+=table(interaction,['metric','n_paired_samples','n_paired_patients','delta_mean','ci_lower','ci_upper','p_exact','p_holm'])+'\n\n'
    text+=('The interaction is [HVG2000 − all] with direct genes minus the same effect after PCA30. '
           'The complete six-level feature curve, seven representations, three clusterers and five-seed controls are retained; '
           'no favorable subset defines a universal winning method. Seeds are algorithmic variation, not additional patients. '
           'KMeans uses n_init10; GMM retains n_init1/max_iter100, diagonal covariance and reg1e-4. '
           '[All convergence warnings and actual parameters](verification/optimization_diagnostics/all_diagnostics.csv.gz) are preserved rather than tuned away.\n\n'
           '[Full factorial](summary/primary_factorial_patient_means.csv) · [Paired method comparisons](summary/paired_method_comparisons.csv) · '
           '[Seed variability](summary/seed_variability.csv) · [Numerical-recovery exclusion sensitivity](summary/numerical_recovery_sensitivity.csv).\n\n')
    text+='## Original data and annotation replay\n\n'
    text+=('Raw 10x contains 110,497 cells. Applying the recovered source thresholds nFeature>200 and mitochondrial fraction<15% to full raw-gene QC metadata leaves 99,896; '
           'the source 7.5% expected doublet counts total 7,492, exactly yielding 92,404 saved cells across all eight samples. '
           'Every retained cell has exactly matching **all 33,694 raw-gene counts, nCount, nFeature and mitochondrial fraction** after '
           'the 15 standard Seurat underscore-to-dash gene-name substitutions. This establishes raw-data and fixed-cell provenance; '
           'it does not recover the original DoubletFinder random state, pK or fitted doublet identities independently.\n\n')
    text+=table(qc,['sample','n_raw_cells','n_after_raw_reference_QC','raw_expected_doublets','n_archived_S2','n_retained_cells_raw_count_exact','gene3_first_QC_shortfall_vs_raw_QC'])+'\n\n'
    text+=('A separate candidate that first removes genes detected in fewer than 3 cells changes QC for one cell each in MT-2, N-1 and N-2. '
           'Its 99,893 pre-doublet cells are not substituted for the full-raw-gene 99,896 definition that matches the archived QC fields and retained totals. '
           'Both candidate definitions remain in the audit; this is evidence of order-sensitive preprocessing, not independent proof of the original doublet model.\n\n')
    text+=('Independent R and Python checks reproduce S2 initial, terminal native and terminal general annotations for all 92,404 cells. '
           'S3 native terminal and detailed-T annotations also exactly match the recovered CSV. '
           'The archived S2 refinement changed none of the initially known cells and resolved 23,036 initially Undecided cells. '
           '**This is exact archived-output replay.** No historical trained DG-scRNA weights or model-initialization seed were recovered, '
           'so bitwise historical retraining equivalence is unavailable. The only model.pth belongs to a FashionMNIST tutorial.\n\n'
           'Original native-label simplification is reapplied: NCOMMREFF drops one prefix field, cancer labels drop three, other labels drop two, '
           'joining remaining fields with `+` so CD8+ remains intact. New strict T/NK/NKT/ambiguous-lymphoid definitions were frozen before new outcome evaluation.\n\n')
    text+='## Two-group and batch reconstruction\n\n'
    text+=('**MTN** contains MT-1, MT-2, N-1, N-2 (48,255 cells); **TUT** contains TU-1, TU-2, T-1, T-2 (44,149). '
           'Patient1=N-1/T-1, Patient2=N-2/T-2, Patient3=TU-1/MT-1, Patient4=TU-2/MT-2. '
           'Each group therefore has one sample from each patient. The original saved checkpoint instead integrated all eight samples together; '
           'the requested two-group workflow is a controlled new intervention. No independent Pu/GSE184362 expression cohort is used.\n\n'
           'No correction, original-style CCA and Harmony use matched cells and the same group-selected 2,000 features. '
           'Matched R single-sample controls use those same selected genes, RNA normalization, DL inputs, UMAP and clustering settings. '
           'This separates pooling from correction without confusing R/Python implementations or different HDBSCAN settings. '
           'Both PCA30 and UMAP2 spaces are tested with SNN resolution 0.5 and R HDBSCAN minPts50. '
           'Sample is the correction covariate, with tissue/biology confounding acknowledged. '
           'CCA corrects expression before PCA; Harmony corrects the matched RNA PCA30 representation before UMAP. '
           'Their comparison changes both correction method and stage, so it does not identify a pure correction-position causal effect. '
           'No undefined UMAP-then-CCA operation is used to fill a comparison cell.\n\n')
    text+=('The fixed representative path gives different orderings across concordance and TCR endpoints. '
           'These four-patient means use AllTissues/mean, UMAP2/R HDBSCAN50, full-RNA scoring and RNA2000 DL throughout; '
           'they do not select the best marker context for each correction:\n\n')
    text+=table(batch_means,['analysis_scope','correction','S2_macro_F1_reference_present','productive_strict_recall','productive_strict_apparent_F1','coverage'])+'\n\n'
    show=batch[batch.space.eq('UMAP2')&batch.clusterer.eq('HDBSCAN_R')&batch.metric.isin(['S2_macro_F1_reference_present','productive_strict_apparent_F1','coverage'])]
    text+=table(show,['group','contrast','metric','delta_mean','ci_lower','ci_upper','p_exact','p_holm'])+'\n\n'
    text+=('[All matched batch contrasts](summary/matched_batch_effects.csv) include PCA30/SNN as well. '
           '[Legacy integrated-assay contrasts](summary/legacy_integrated_effects.csv) separately change CCA scoring and DL features; '
           'they do not enter the full-RNA main correction comparison.\n\n')
    text+='### Mixing and RNA-state preservation\n\n'
    text+=('A fixed balanced selection of up to 500 cells per sample within each shared lineage uses 30-neighbor sample entropy '
           'and 15-neighbor within-sample neighborhood retention versus uncorrected data. Productive-TCR-positive cells and shared S2 lineages '
           'are evaluated separately. RNA proliferation, interferon, stress, CD4, CD8 and Treg modules provide descriptive state diagnostics; '
           'they are not independent truth and may share annotation markers.\n\n')
    text+=table(bio[bio.scope.eq('productive_TCR_positive')],['group','correction','space','n_cells','mean_batch_entropy','within_sample_neighbor_retention_vs_NONE'])+'\n\n'
    text+=('The prespecified [RNA marker evidence](figures/ptc_RNA_marker_evidence.pdf) shows CD3D/E/TRAC, NK/B/myeloid/epithelial and other markers '
           'for called T, Unknown and other called cells stratified by productive TCR detection. '
           '[Patient/tissue composition](figures/ptc_patient_tissue_composition.pdf) compares initial and terminal calls; '
           '[full-cell embeddings](figures/ptc_embedding_TCR_lineage.pdf) display terminal calls and detection separately. '
           '[Unknown/QC strata](biology_detail/TCR_stratified_Unknown_QC_by_sample.csv) preserve exact source QC; historical doublet scores were unavailable.\n\n')
    text+='## Original marker roster, refinement and context selection\n\n'
    text+=('All **17 symbol-based libraries × three cutoffs** are restored: CellMarker Lymph, Lymph node, Lymphoid tissue, Thyroid, Epithelium, '
           'Blood, Thymus and AllTissues; Pubmed_34663816/GSE184362-derived markers; HPA_allThyroid and seven HPA evidence strata. '
           'Symbols, duplicates and full panel denominators remain original. The separate Ensembl marker RDS is not interchangeable. '
           'The manuscript refers to 16 libraries; the recovered native roster contains 17, and the complete actual roster is reported.\n\n'
           'The score uses explicit Seurat-v4-compatible Wilcoxon and log2(mean(expm1(x))+1) fold changes; '
           'logfc.threshold=0.25, min.pct=0.1, raw p<0.01 precede scoring log2FC>1. '
           'Density sums hit logFC over the full marker-panel denominator, with the original 0.8 singleton-marker penalty; '
           'ties remain Undecided. Cutoffs are none, mean and 0.5. A marker-union DE testing optimization is exactly checked '
           'against full-gene tests, without changing normalization or denominators.\n\n'
           'The archived MLP is 256/128 LeakyReLU, Softmax before CrossEntropy, Adamax1e-3, ten epochs, batch size 256, '
           '90/10 known-cell split seed 42. New model initialization is explicitly seed 42. Original known calls remain unchanged; '
           'rounded four-decimal confidence uses 0.90, with 0.70 from the same probabilities as manuscript sensitivity. '
           'The independent verifier reconstructs every trained cache forward pass and both terminal thresholds, cell orders, input hashes and splits.\n\n')
    vocabulary=pd.read_csv(BASE/'vocabulary_audit/original_library_vocabulary.csv')
    text+=('Libraries differ in what cell types they can represent. The [sample-specific vocabulary-capacity audit](vocabulary_audit/reference_vocabulary_capacity_by_sample.csv) '
           'reports supported S2/S3 classes and the theoretical reference-present macro-F1 ceiling under perfect supported calls plus Unknown abstention. '
           'This is representational capacity, not achieved accuracy; a favorable binary-T score need not imply broad multi-class coverage.\n\n')
    text+=table(vocabulary,['library','n_native_panels','n_distinct_broad_lineages','contains_strict_T'])+'\n\n'
    text+=table(census)+'\n\nUnique model/cache census:\n\n'+table(caches)+'\n\n'
    text+=('For context selection, three patients choose among the 17×3 panels by patient-mean apparent productive-TCR binary F1; '
           'the fourth patient provides the held-out endpoints. Deterministic ties use library name then mean/none/0.5. '
           '**This is patient-held-out panel selection on transductively fitted predictions**, not fully inductive validation: '
           'embeddings, clusterings and DL models were not fitted afresh with the held-out patient excluded.\n\n')
    chosen=panels[panels.space.eq('UMAP2')&panels.clusterer.eq('HDBSCAN_R')]
    text+=table(chosen,['group','correction','policy','productive_strict_apparent_F1','productive_strict_recall','coverage'])+'\n\n'
    text+='## Source and implementation discrepancies retained in the evidence\n\n'
    text+=('| Item | Recovered evidence / treatment |\n|---|---|\n'
      '| Normalization | Manuscript TMM versus archived Seurat LogNormalize10k; reconstruction follows the executed source/checkpoint. |\n'
      '| Feature selection | Manuscript paragraph600 and reference code use HVGs. A blanket statement that the original workflow never used HVGs is not supported. |\n'
      '| DL cutoff | Code0.90 and manuscript0.70 both reported from the same saved probabilities. |\n'
      '| Model | Executed dense MLP, not graph/RL reinforcement learning; user “RL” means terminal DL/refinement. |\n'
      '| Historical integration | All eight samples in one CCA versus the explicitly requested MTN/TUT groups. |\n'
      '| TCR definitions | S2 any-contig35,727; productive/high-confidence34,807; pairedTRA/TRB30,888; S3 supplied36,184. |\n'
      '| S3 source versions | 5,435 per-cell TCR differences and3,901 DG binary-flag differences remain unresolved; keep native calls and literal flags separate. |\n'
      '| SCINA names | All623 apparent label differences are UTF8/Latin1 γδ encoding only; compact spellings are normalized through the pre-frozen ontology. |\n'
      '| Library count | Native RDS17 versus manuscript16; all native libraries retained. |\n'
      '| Noise/index adapter | PTC R noise0 is scored as an observed group. Actual cluster IDs are iterated rather than assuming contiguous0..L−1; exact original density parity was verified. |\n'
      '| GMM numerical issue | N-2/HVG2000/t-SNE2 float32 failure reproduced, precision-only float64 recovery at unchanged parameters; full and exclusion-sensitivity results retained. |\n'
      '| R HDBSCAN crash | Official1.2.6 fixes a root write, but MTN48,255 still crosses a 32-bit index-product boundary. Private1.2.6.9001 casts indices toR_xlen_t;100,004 boundary checks, small-data exactness and all24 TUT old/new partitions pass. Not OOM. |\n'
      '| Refinement runtime recovery | Empty validation-index dtype and cross-node cache publication fixed without model changes; duplicate results compared and original artifacts retained. |\n\n')
    text+=('Primary software evidence: [Seurat v4-compatible marker testing](https://satijalab.org/seurat/reference/findallmarkers), '
           '[dbscan1.2.6 release](https://github.com/mhahsler/dbscan/releases/tag/dbscan_1.2.6), '
           '[MST implementation](https://github.com/mhahsler/dbscan/blob/dbscan_1.2.6/src/mst.cpp), '
           '[lower-triangle index macro](https://github.com/mhahsler/dbscan/blob/dbscan_1.2.6/src/lt.h). '
           'Original manuscript/source discrepancies are documented against `paper/submission_v16/DG_scRNA_04232026_V16_cell_report.docx` '
           'and the staged native scripts, not inferred from a modern vignette.\n\n')
    text+='## Files and reproducibility\n\n'
    text+=('- [Complete decision tree](figures/decision_tree_complete.pdf) and [figure inventory](figures/figure_manifest.json).\n'
           '- [Complete condition census](summary/condition_verification.csv.gz), [all sample-stage endpoints](summary/all_sample_stage_metrics.csv.gz), '
           '[actual features](verification/geometry_census/actual_features_and_geometry.csv), [physical partitions](verification/geometry_census/physical_partitions.csv).\n'
           '- [Archived competitors with corrected names](archived_comparators/historical_TCR_endpoint_replay_corrected_names.csv); '
           'original inputs/reference conditions are retained and are not presented as a matched ranking against the new pipeline.\n'
           '- [English protocol](../../../handoff/ptc_recovery_20260916/CONTROLLED_PROTOCOL.md), '
           '[analysis plan](../../../handoff/ptc_recovery_20260916/FINAL_ANALYSIS_PLAN.md), '
           '[source/table reconciliation](../ptc_recovery/R_baseline_replay/replay.json).\n'
           '- New and preserved GBM chapters appear in the executed central notebook. All numerical work, extraction, R, fitting, verification and figures ran via SLURM.\n')
    (BASE/'PTC_REPORT.md').write_text(text)
    # Chinese briefing answers the user's PI question without preselecting its direction.
    gbm=BASE.parent
    zh='# GBM + PTC：给 PI 的实验结论与讨论材料\n\n'
    zh+=('所有 DG-scRNA 主结果都使用 **DL/refinement 后的最终输出**；marker-only 是明确命名的消融。'
         'GBM 已先完整交付，PTC 原始注释核对通过后才开展新实验。\n\n'
         '## 可以怎样回应“不用 HVG”\n\n'
         '现有证据支持把问题拆到具体阶段。GBM 在固定 direct UMAP2→HDBSCAN15/15 路径中，HVG2000 相对全基因：'
         '患者配对 ARI +0.1468、最终 strict-L1 macro-F1 +0.1570；97 个样本、55 位患者，两个主终点 Holm 校正后均有明确证据。'
         '但 PCA30→UMAP2 路径的最终 F1 差异仅−0.0020，区间跨0。'
         '因此可以反驳“所有路径都不需要 HVG”这一笼统建议，不能把结果说成“HVG 或 UMAP+HDBSCAN 永远最好”。\n\n'
         '还必须区分 geometry 和 marker scoring：GBM 把 scoring 限为 HVG2000 后，三个固定路径的最终 F1 下降约0.044–0.054。'
         '更有证据的建议是按流程验证 geometry 特征选择，同时保留完整 marker scoring 基因。'
         'DGCyTOF 原文 Table4/5 比较的是13/32个蛋白通道；repo 中另有 scRNA HVG 扫描，不能把二者混为原文 HVG 实验。\n\n'
         '没有找到 PI“一律不取 HVG”的逐字转录；以下检验的是用户转述的建议，不制造 PI 引语。\n\n'
         '## PTC 完整补充\n\n'
         '92,404 个细胞的原始 counts、QC 和 S2/S3 注释列逐细胞复核一致。历史模型权重缺失，所以准确表述是“原始输出复核一致”，'
         '不是“从头训练逐位复现”。新实验按 MT+N 和 TU+T 两组，包含原始17套 marker×3 cutoff、NONE/CCA/Harmony、'
         '匹配参数的 R 单样本对照、完整 HVG×降维×聚类、refinement/阈值消融和 TCR/RNA/batch 生物学诊断。\n\n'
         f'共 **11,696 个 terminal 条件、49,776 行样本/阶段指标、{summary["n_unique_trained_models"]:,} 个不同训练模型 cache** 完成独立核验。'
         '以下为 PTC 的 HVG2000−all 患者配对效应：\n\n')
    direct=h[h.path.eq('direct_UMAP2')&h.feature.eq('hvg2000')&h.metric.eq('S2_macro_F1_reference_present')].iloc[0]
    zh+=(f'PTC 在 direct UMAP2/HDBSCAN 路径的最终 S2 broad macro-F1 配对差异为 **{direct.delta_mean:+.4f}**；'
         f'PCA30→UMAP2 后为 **{pca.delta_mean:+.4f}**。这里同样应讲清楚具体路径，不能说所有情况下都要或都不要 HVG。\n\n')
    zh+=table(h[h.feature.eq('hvg2000')],cols)+'\n\n'
    zh+=('PTC 只有4位患者，精确双侧检验最小 p=0.125；不能用92,404个细胞或5个seed放大显著性。'
         'S2/S3 是旧模型注释，一致性不等于独立真值准确率。TCR 阴性也不等于非T细胞。'
         'batch 结果要同时看 TCR recall、coverage、RNA 状态和邻域保留；不能仅凭混合更好判定胜出。\n\n'
         '固定 AllTissues/mean、R UMAP2/HDBSCAN50 的 batch 对照如下。不同终点可以给出不同排序，所以既要看旧注释一致性，也要看 TCR 和 RNA 支持：\n\n')
    zh+=table(batch_means,['analysis_scope','correction','S2_macro_F1_reference_present','productive_strict_recall','productive_strict_apparent_F1','coverage'])+'\n\n'
    zh+=('### 与原稿必须统一的地方\n\n'
         '- 原稿 TMM 与实际 LogNormalize10k、0.70 与代码0.90、16套与实际17套 marker、原来8样本一起 CCA 与新两组设计，都有独立标注。\n'
         '- S3 的 TCR 二值列有5,435个细胞与 S2 不同，DG二值T标记有3,901个差异；保留版本差异，不改标签凑一致。SCINA的623条差异已证实只是字符编码。\n'
         '- R HDBSCAN 两次段错误是软件索引问题，已隔离修复并做边界与完整 TUT 一致性核验；没有把失败算成方法劣势。\n'
         '- held-out panel 选择只对选择环节留出患者；预测是 transductive 拟合，不冒充完全独立外部验证。\n\n'
         '建议 meeting 顺序：完整 decision tree → GBM HVG×PCA 四方比较 → geometry/scoring 分离 → PTC 原始核对 → '
         '匹配 R 的单样本/pooled/batch 对照 → marker 与 refinement 消融 → 小患者数与来源差异限制。\n\n'
         '[PTC 全报告](PTC_REPORT.md) · [完整流程主图](figures/decision_tree_complete.pdf) · '
         '[GBM 报告](../GBM_REPORT.md) · [执行后的 notebook](../notebooks/dgscrna_v7_results.ipynb)\n')
    (BASE/'PI_BRIEF_ZH.md').write_text(zh)
    write_json(BASE/'report_manifest.json',dict(status='complete',source_sha256=sha(Path(__file__)),
      summary_manifest_sha256=sha(BASE/'summary/manifest.json'),outputs={n:sha(BASE/n) for n in ['PTC_REPORT.md','PI_BRIEF_ZH.md']},
      completed_at=utc(),job=os.environ['SLURM_JOB_ID']))
    print('PTC report and combined Chinese PI briefing written',flush=True)

if __name__=='__main__':run()
