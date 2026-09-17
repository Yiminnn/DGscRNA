"""Generate the scientific report and execute only appended cells of the existing notebook."""
import base64,hashlib,json,os,sys,shutil,fcntl
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'
TAG='r_reference_campaign_20260917'

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def markdown(frame):
    def fmt(v):
        if isinstance(v,float):return '' if v!=v else f'{v:.4f}'
        return str(v).replace('|',' / ').replace('\n',' ')
    return '\n'.join(['| '+' | '.join(map(str,frame.columns))+' |','| '+' | '.join(['---']*len(frame.columns))+' |',
        *['| '+' | '.join(fmt(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None)]])

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import pandas as pd
    import nbformat
    from nbclient import NotebookClient
    from evaluation_rules import canonical,broad
    dest=OUT/'summary';s=json.loads((dest/'campaign_summary.json').read_text())
    partial='--allow-partial' in sys.argv
    if not partial:assert s['status']=='complete',s['status']
    inventory=pd.read_csv(dest/'completion_inventory.csv')
    roster=json.loads((OUT/'markers/marker_roster.json').read_text())
    marker_rows=[]
    for item in roster['units']:
        marker_rows.append(dict(dataset=item['dataset'],unit=item['unit'],primary='; '.join(item['primary']),
          related='; '.join(item['support']),n_libraries=len(item['libraries']),
          libraries='; '.join(item['libraries']),selection_basis=item['selection_basis']))
    marker_table=pd.DataFrame(marker_rows)
    marker_table.to_csv(dest/'marker_context_roster.csv',index=False)
    reviewer_contexts=marker_table[marker_table.dataset.ne('HCL')][['dataset','primary','related','n_libraries']]
    inputs=[json.loads(p.read_text()) for p in OUT.glob('inputs/*/*/input_manifest.json')]
    input_table=pd.DataFrame([dict(dataset=m['dataset'],unit=m['unit'],n_cells=m['n_cells'],n_genes=m['n_genes'],
      n_batches=len(m['batch_sizes']),input_semantics=m['input_semantics'],source=m['source']) for m in inputs])
    input_table.to_csv(dest/'input_roster.csv',index=False)
    cohort=input_table.groupby('dataset',as_index=False).agg(cells=('n_cells','sum'),analysis_units=('unit','nunique'))
    assert len(cohort)==11 and int(cohort.set_index('dataset').loc['HCL','cells'])==599926
    cohort.to_csv(dest/'cohort_sizes.csv',index=False)
    allm=pd.read_csv(dest/'all_annotation_metrics.csv.gz',low_memory=False)
    best=pd.read_csv(dest/'benchmark_descriptive_maxima.csv')
    b=allm[allm.dataset.ne('PTC')&allm.dataset.ne('HCL')&allm.stage.eq('final090')&allm.endpoint.eq('common_lineage')]
    fixed=b[b.library.eq('CM2_primary_normal')&b.cutoff.eq('mean')]
    comparison=fixed.pivot(index='dataset',columns='route',values='macro_F1').reset_index() if len(fixed) else pd.DataFrame()
    descriptive=best[best.dataset.ne('HCL')&best.endpoint.eq('common_lineage')][['dataset','n_tested','descriptive_maximum_macro_F1','maximum_route','maximum_library','maximum_cutoff']]
    descriptive.to_csv(dest/'reviewer_descriptive_maxima_common_lineage.csv',index=False)
    p=pd.read_csv(dest/'PTC_fixed_historical_marker_contexts.csv')
    select=((p.scope.eq('NMT')&p.route.isin(['PCA30_SNN','seurat_clusters'])) |
            (p.scope.eq('TTU')&p.route.isin(['UMAP2_HDBSCAN_R','hdbscan.UMAP_clusters'])))
    columns=['unit','scope','F1_T','AUC_binary','accuracy','unknown_fraction','saved_native_concordance','saved_broad_lineage_concordance']
    columns += [c for c in ['accuracy_unknown_as_error','macro_F1_T_nonT_unknown_as_error','coverage'] if c in p.columns]
    selected=p[select&p.stage.eq('final090')&p.endpoint.eq('strict_T_name_rule')][columns]
    selected=selected.rename(columns={'accuracy':'historical_binary_mapping_accuracy'})
    selected.to_csv(dest/'PTC_fixed_paper_route_terminal.csv',index=False)
    statuses=pd.DataFrame(list(s['terminal_status_counts'].items()),columns=['terminal_state','conditions'])
    ptc_maxima=pd.read_csv(dest/'PTC_descriptive_maxima_by_family.csv')
    ptc_maxima_table=ptc_maxima[ptc_maxima.metric.eq('macro_F1_T_nonT_unknown_as_error')][
        ['group','family','value','unit','route','library','cutoff','unknown_fraction']]
    # A complete mapping audit exposes resolution mismatches; it does not change evaluation.
    pm=pd.read_csv(OUT/'markers/native_panel_metadata.csv',keep_default_na=False)
    mapping=[];truthrows=[]
    for dataset in sorted(input_table.dataset.unique()):
        for row in pm.itertuples(index=False):
            mapping.append(dict(dataset=dataset,native=row.native,cell_name=row.cell_name,
                curated_semantic=canonical(row.cell_name),common_lineage=broad(row.cell_name,dataset)))
        for m in inputs:
            if m['dataset']!=dataset:continue
            ref=pd.read_csv(OUT/'inputs'/dataset/m['unit']/'evaluation_only.csv.gz',keep_default_na=False)
            for label,n in ref.truth.value_counts().items():truthrows.append(dict(dataset=dataset,unit=m['unit'],truth=label,n_cells=int(n),curated_semantic=canonical(label),common_lineage=broad(label,dataset)))
    pd.DataFrame(mapping).to_csv(dest/'native_marker_label_mapping.csv.gz',index=False)
    pd.DataFrame(truthrows).to_csv(dest/'curated_truth_label_mapping.csv',index=False)
    report=f'''# R-reference PTC ablation and reviewer datasets — 2026-09-17

Status: **{s['status']}**. Evaluated {s['evaluated_units']}/{s['expected_analysis_units']} analysis units; plotted {s['plotted_units']}; audited {s.get('audited_units',0)}. Completed terminal annotation conditions: **{s['terminal_annotation_conditions_evaluated']}**. Clustering evaluations: **{s['clustering_conditions_evaluated']}**.

This report extends the existing notebook and preserves its GBM results. The full grid, including poor results, no-op models, untrainable conditions and Unknown calls, remains available. No independent Pu validation cohort was added.

## Cohorts and marker selection

{markdown(cohort)}

The reviewer cohort roster follows the previously approved human datasets in `handoff/deck_datasets_provenance.md`. HCL is analyzed in all 59 original tissue groups, retaining 599,926 cells. Its aggregate is a **tissue-conditional** analysis; it is not a pooled 600k-cell CCA/HDBSCAN scalability benchmark. No claim of that benchmark is made here.

For each reviewer unit the CellMarker candidates were fixed from sampled anatomy and normal/disease context before scoring: primary tissue; disease-specific primary panels where appropriate; relevant blood/lymphoid, vascular, stromal or sampled extranodal tissues; related unions; AllHuman. Distinct native normal/cancer/tissue/type panels and full gene denominators are preserved. PTC retains all 17 original archived libraries, including extra-thyroid CellMarker tissues, HPA and the GSE184362-derived Pubmed library. See `marker_context_roster.csv`, `markers/marker_roster.json` and `markers/native_panel_metadata.csv` for exact choices and evidence PMIDs.

{markdown(reviewer_contexts)}

The library count includes non-empty primary, related-tissue, union and AllHuman contexts. Empty database contexts are recorded and excluded before fitting; they are not silently replaced. HCL's 59 tissue-specific rosters are listed individually in the complete roster.

CellMarker source: [CellMarker 2.0](https://bio-bigdata.hrbmu.edu.cn/CellMarker2.0/index.html), [database publication](https://doi.org/10.1093/nar/gkac947). Cached workbook SHA256: `{roster['source_sha256']}`. Source studies may occur in marker evidence; these are not independent-reference generalization experiments.

## What is actually compared

- PTC archived all-eight-sample reference: four clustering branches × 17 marker libraries × three density cutoffs = 204 terminal arms. The selected NMT Thyroid/PCA-SNN/none and TTU Pubmed/UMAP-HDBSCAN/mean routes reproduce the previously validated refit, including probabilities. They do **not** recover the historical model weights or erase the previously documented differences from Sup labels.
- Seventeen newly prepared PTC conditions compare joint CCA gene budgets (500/1k/2k/3k/5k/all), no correction/Harmony, and all-eight versus separate NMT/TTU integration.
- Twelve additional conditions isolate geometry gene budget on a fixed all-gene CCA expression matrix. Scoring and DL use the same fixed 2,000 genes and byte-identical input. This contrast is conditional on that all-gene CCA fit.
- Reviewer units use reference 2,000-gene preparation, four clustering branches, all prespecified marker contexts and all three cutoffs. Single-batch units retain RNA scoring; adaptive dimension/neighbor caps for small batches are logged. This campaign does not replace the earlier GBM HVG grid with a new R HVG sweep.
- Marker-only, terminal confidence 0.90, and 0.70 sensitivity are all reported. “DG-scRNA final” means terminal DL/refinement. A valid no-op is distinct from executed training; an absent scientific result is never scored as a failure of the method.

QC cells, normalization, R-compatible DEG statistics, the original density formula and MLP parameters are fixed. “All genes” means genes detected in at least three cells in every included PTC sample, without variance ranking. Changing the full CCA budget changes anchor, geometry, scoring and DL genes together; it cannot isolate a geometry mechanism. RNA versus Harmony holds the scoring/DL expression fixed. CCA versus RNA changes expression and geometry together; `verification/design_audit.json` records actual gene identities and cell order.

Muraro uses inspected raw.X integer counts. Xin supplies RPKM. Immune_ALL uses its prepared matrix from the published mixed UMI/full-length count layer; fractional values are retained. None is silently substituted with scaled X or rounded into UMI counts. See `input_roster.csv` and input manifests.

## Reviewer results: fixed context and complete grid

Fixed illustrative context: primary-normal CellMarker, density cutoff mean, terminal confidence 0.90. All four routes are retained below; the context was not chosen from a performance maximum.

{markdown(comparison)}

The following maxima are **descriptive best-on-these-labels values** over the entire measured marker/route/cutoff grid. They are not independently selected operating points and cannot establish a global optimum.

{markdown(descriptive)}

The primary curated-semantic and secondary common-lineage metrics are different resolutions, with different effective class counts. They must not be mixed. Unknowns and predictions outside the truth vocabulary remain in the denominators; unmatched labels are not merged into a mutually correct “Other” class. `native_marker_label_mapping.csv.gz` and `curated_truth_label_mapping.csv` expose every transformation. Per-class results and confusion matrices are retained for each unit. HCL summaries explicitly report the number of tissue groups and distinguish equal-tissue mean F1, cell-weighted mean tissue F1 and pooled accuracy.

Where at least three donor labels are available, `donor_heldout_label_selection.csv` chooses a route/library/cutoff using only the other donors' labels and scores the held-out labels. This is a transductive fixed-cohort sensitivity analysis: all cells entered unsupervised fitting, so it is not an unseen-donor training/refit experiment. Do not treat random seeds or tissue partitions as independent patients.

## PTC: fixed paper-selected marker and clustering routes

The table keeps the paper-selected NMT and TTU routes fixed across interventions. TCR detection is an imperfect independent positive proxy; absence of detected TCR is not established non-T ground truth. Strict T-name and historical broad T/NK compatibility endpoints are both saved. Agreement with archived native/broad labels is **concordance**, not accuracy against independent truth.

The historical binary mapping treats Unknown as non-T. Its ordinary accuracy and non-T F1 can therefore reward abstention. `PTC_abstention_aware_metrics.csv.gz` additionally treats Unknown as an error in accuracy and as a false negative for its true class, reporting coverage and two-class macro F1. T-positive F1 is unchanged by this treatment. A completely unresolved arm has zero abstention-aware accuracy and macro F1, even if its historical non-T F1 appears favorable.

{markdown(selected)}

`PTC_paired_ablation_changes.csv.gz` matches route, marker, cutoff, scope and endpoint between interventions. `DL_vs_initial_paired_changes.csv.gz` isolates the annotation changes after DL on the same fitted inputs. Full condition-level metrics are in `all_annotation_metrics.csv.gz`; all clustering results are in the separate `all_clustering_metrics.csv`.

`PTC_descriptive_maxima_by_family.csv` reports the best observed T-positive F1 and abstention-aware metrics separately for archived geometry, fresh all-eight integration, joint CCA budgets, RNA/Harmony and isolated geometry. Keeping those families separate prevents a conditional geometry control from being presented as the unchanged original pipeline.

The earlier paper restoration diagnosis remains valid: original saved Sup endpoints were identified; original F1/AUC definitions were reconciled; historical DL weights and the manuscript Accuracy computation were not recovered. New ablations use the validated refit and explicitly named fresh integrations. They are not a claim of exact historical retraining.

## Terminal states and verification

{markdown(statuses)}

Every complete analysis unit has model/NPZ/history checksum checks, prediction-export equality, known-label retention, exact threshold reconstruction, and cell-order checks. The original selected PTC refits have zero probability difference from the prior validated run. Early reviewer scoring jobs started before immutable per-job source guards; their saved outputs are separately replayed with the unmodified original density function. That replay verifies outputs and does not retroactively reconstruct transient script bytes. Later R and terminal jobs retain executed source copies. Failed attempts and resource-only retries are preserved in SLURM logs and submission records.

Two implementation safeguards are explicit: R dbscan uses the previously validated 64-bit MST-index patch for cohorts above the 32-bit index limit, and density scores index the observed cluster IDs rather than assuming that noise label 0 must exist. The scoring formula is preserved; arbitrary ID remapping is used when replaying the legacy function in validation. Fixed curated cells, logged small-batch adaptations and current pinned software versions remain part of this reference execution, not a claim of historical binary identity.

## How to use the results in the paper

Use “best among the evaluated configurations on [dataset], under [endpoint and marker context]” only where the table supports it. A ranking on the same labels is descriptive optimization. Do not call every node optimal: untested QC, resolution, minPts, UMAP and MLP settings remain fixed; the grid is not fully factorial. The isolated geometry controls support a narrower feature-budget interpretation than the full-workflow budget sweep. A favorable clustering ARI/NMI does not establish superior terminal annotation.

Notebook: `notebooks/dgscrna_results.ipynb` in the existing delivery root. The English decision tree is `summary/workflow_decision_tree.svg` (editable), with PDF and PNG exports. Every analysis unit has an all-four-clustering atlas, fixed-context marker-to-DL plots, and complete marker-grid heatmaps in its `figures` directory.
'''
    (dest/'RESULTS_AND_INTERPRETATION.md').write_text(report)
    (dest/'RESULTS_ZH.md').write_text(f'''# PTC 与 reviewer 数据集：R 原始流程重跑

状态：**{s['status']}**；完成评价 {s['evaluated_units']}/99 个分析单元，{s['terminal_annotation_conditions_evaluated']} 个最终注释组合；图已完成 {s['plotted_units']}，完整性审计 {s.get('audited_units',0)}。

PTC 保留原始 17 个 marker 库，NMT 与 TTU 的论文选定 marker 路线分别固定报告。除整体 CCA/HVG 数量比较外，另有固定校正表达、marker 打分与 DL 输入、只改变聚类特征的对照。RNA/Harmony 的 DL 输入完全一致；CCA/RNA 同时改变表达和聚类空间。

Reviewer 沿用已确认的 11 套人类数据，先按真实采样组织固定 CellMarker 候选，再报告全部组合。HCL 保留全部 599,926 个细胞，按 59 个原始组织组分别运行；不能把它写成一次 60 万细胞整体聚类的性能证明。

所有最终 DG-scRNA 数值都来自 terminal DL，保留 Unknown、无需训练及无法训练的状态。原始 saved labels 的对齐属于 concordance；TCR 的阴性不等于确定非 T 细胞。原论文历史权重/Accuracy 来源尚未恢复，这次的新实验不改变此前诊断。

PTC 同时给出原始二分类指标和 Unknown 记为错误的指标：原规则把 Unknown 算作非 T，不能仅凭这种规则下较高的非 T F1/accuracy 判断方法更优。

## 实际采用的 reviewer marker 组织

下表列出主组织、相关组织及实际非空 marker 库数量；同时保留相关组织联合库与 AllHuman 对照。癌症数据另有主组织 normal、disease 与 all-context 库。HCL 的 59 个组织单独列在完整 roster 中。

{markdown(reviewer_contexts)}

## 固定 marker 后的四种聚类路线

以下均固定 primary-normal CellMarker、mean cutoff 和最终 DL 0.90 阈值，数值为 common-lineage macro F1。因此可以直接比较路线；不可与细粒度 curated-semantic F1 混用。

{markdown(comparison)}

## PTC 各对照家族内的观测最好结果

下表按 Unknown 记为错误的 T/non-T macro F1 排序，分别保留 NMT、TTU 及各类干预。它描述已经测到的组合，使用了同一套评价标签；各家族的输入不同，不能混称为原始流程的单一最优配置。TCR 检出仍只是评价代理。

{markdown(ptc_maxima_table)}

目前只能根据完整表格说明「在这个数据集、marker 和指标下，所测组合中的最好结果」，不能把整条流程所有步骤都写成最优。固定参数、未测交互及语义标签分辨率必须交代。全部结果与解释见 [RESULTS_AND_INTERPRETATION.md](RESULTS_AND_INTERPRETATION.md)，实际基因与输入一致性见 `../verification/design_audit.json`。
''')
    if '--report-only' in sys.argv:
        print('REPORT_WRITTEN',s['status'],flush=True);return
    path=ROOT/'notebooks/dgscrna_results.ipynb'
    with (OUT/'notebook_update.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        original_hash=sha(path)
        nb=nbformat.read(path,as_version=4)
        backup=OUT/'notebook_before_campaign.ipynb'
        if not backup.exists():shutil.copy2(path,backup)
        baseline=[c for c in nb.cells if TAG not in c.metadata.get('tags',[])]
        old=nbformat.read(backup,as_version=4)
        assert len(baseline)>=len(old.cells)
        for i,cell in enumerate(old.cells):assert cell==baseline[i],f'Existing cell {i} changed outside campaign'
        new=[]
        def md(text):new.append(nbformat.v4.new_markdown_cell(text,metadata={'tags':[TAG]}))
        def code(text):new.append(nbformat.v4.new_code_cell(text,metadata={'tags':[TAG]}))
        md(f'''<a id="r-reference-campaign-20260917"></a>

# R-reference PTC ablations and reviewer cohorts — 2026-09-17

**{s['status']}**: {s['evaluated_units']}/99 analysis units evaluated; {s['terminal_annotation_conditions_evaluated']} terminal annotation conditions. This section extends the original notebook; every preceding GBM/PTC cell is preserved.

PTC: 17 original libraries, 30 preparation/reference units including isolated geometry controls. Reviewer: the 11 approved human cohorts, with prespecified CellMarker contexts; HCL is conditional on its 59 original tissue groups. Original saved labels, validated refit, and fresh integrations remain distinct. Historical model weights and the manuscript Accuracy provenance remain unresolved.

All metrics retain abstentions. Descriptive maxima are not independent operating-point selection. Full protocol, label mappings, condition tables and audit records are available under `results/hvg_ptc_20260916_v1/r_reference_campaign_20260917/`.''')
        code('''import os
from pathlib import Path
import pandas as pd
from IPython.display import display, Image, Markdown
assert os.environ.get("SLURM_JOB_ID"), "Run scientific notebook cells inside a SLURM allocation"
ROOT = Path.cwd().parent if Path.cwd().name == "notebooks" else Path.cwd()
CAMPAIGN = ROOT / "results/hvg_ptc_20260916_v1/r_reference_campaign_20260917"
def show_figure(relative, width=1000):
    display(Image(filename=str(CAMPAIGN / relative), width=width))
display(pd.read_csv(CAMPAIGN / "summary/cohort_sizes.csv"))
show_figure("summary/workflow_decision_tree.png")''')
        md('## Complete marker roster and endpoints\n\nCandidates were fixed from anatomy and disease context before scoring. Native context prefixes and full gene denominators remain unchanged. Curated-semantic F1 and common-lineage F1 have different class resolutions; PTC TCR agreement and archived-label concordance have different interpretations.')
        code('display(pd.read_csv(CAMPAIGN / "summary/marker_context_roster.csv"))\ndisplay(pd.read_csv(CAMPAIGN / "summary/completion_inventory.csv"))')
        for image_name,title in [('reviewer_fixed_context_comparison','Reviewer cohorts: fixed primary-normal / mean context'),('PTC_joint_CCA_budget','PTC: joint CCA feature-budget intervention'),('PTC_geometry_only','PTC: isolated geometry features on fixed CCA expression')]:
            if (dest/(image_name+'.png')).exists():
                md('## '+title)
                code(f'show_figure("summary/{image_name}.png")')
        md('## Full-grid descriptive optima and fixed PTC routes\n\nThese observed maxima use the same evaluation labels. They support a within-dataset comparison, not a global-optimum claim. The PTC table keeps the historical group-specific route fixed.')
        code('display(pd.read_csv(CAMPAIGN / "summary/reviewer_descriptive_maxima_common_lineage.csv"))\ndisplay(pd.read_csv(CAMPAIGN / "summary/PTC_fixed_paper_route_terminal.csv"))\ndisplay(pd.read_csv(CAMPAIGN / "summary/PTC_descriptive_maxima_by_family.csv"))')
        md('## Every clustering branch and terminal marker grid\n\nEach atlas uses one saved UMAP for visual comparison of all four clustering partitions. Numbers identify clusters; legends identify biological labels and donors. Annotation illustration uses a fixed context, while the heatmaps retain the entire marker/cutoff grid. Individual PDF exports and marker-to-DL panels are saved alongside each atlas.')
        ready=inventory[inventory.evaluated&inventory.plotted].copy()
        ready['sort_group']=ready.unit.map(lambda u:0 if u=='brain_GBM' else (1 if u.startswith('PTC_') else (3 if u.startswith('HCL__') else 2)))
        for row in ready.sort_values(['sort_group','unit']).itertuples(index=False):
            prep=Path(row.directory);rel=prep.relative_to(OUT)
            md(f'### {row.unit}\n\n{int(row.n_cells):,} cells; correction `{row.correction}`. Four clustering branches and the complete terminal grid follow.')
            figures=[str(rel/'figures/all_cluster_branches.png')]
            figures += [str(p.relative_to(OUT)) for p in sorted((prep/'figures').glob('terminal_marker_grid_*.png'))]
            if row.unit in ['brain_GBM','PTC_archived_CCA2000']:figures.append(str(rel/'figures/marker_to_terminal_fixed_context.png'))
            code('\n'.join('show_figure('+repr(f)+')' for f in figures))
        md('## Interpretation and audit\n\n'+(dest/'RESULTS_ZH.md').read_text().split('\n',1)[1])
        code('display(pd.read_csv(CAMPAIGN / "summary/all_terminal_statuses.csv").groupby(["dataset", "dl_status"]).size().rename("conditions").to_frame())\ndisplay(pd.read_csv(CAMPAIGN / "verification/actual_feature_and_input_contrasts.csv"))')
        campaign=nbformat.v4.new_notebook(cells=new,metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'}})
        NotebookClient(campaign,timeout=1200,kernel_name='python3',resources={'metadata':{'path':str(ROOT/'notebooks')}},allow_errors=False).execute()
        assert sha(path)==original_hash,'Notebook changed during execution; refusing to overwrite external edits'
        banner=nbformat.v4.new_markdown_cell(
          f'**2026-09-17 最新实验状态：{s["status"]}。** 已完成评价 {s["evaluated_units"]}/99 个分析单元。'
          ' [本轮 R 参考实验、完整流程图与所有聚类图](#r-reference-campaign-20260917)位于文末；'
          '本轮状态以该节为准。此前所有 GBM/PTC 单元和结果完整保留。',metadata={'tags':[TAG]})
        nb.cells=[banner]+baseline+campaign.cells
        temp=path.with_suffix('.ipynb.campaign_part');nbformat.write(nb,temp);temp.replace(path)
        current=nbformat.read(path,as_version=4)
        assert current.cells[1:1+len(baseline)]==baseline
        record=dict(status=s['status'],job=os.environ['SLURM_JOB_ID'],notebook=str(path),
          old_cells_preserved=len(baseline),campaign_cells=len(campaign.cells)+1,all_new_code_cells_executed=True,
          errors=0,backup=str(backup),previous_sha256=original_hash,current_sha256=sha(path),
          no_existing_GBM_or_PTC_cells_removed=True,report_sha256=sha(dest/'RESULTS_AND_INTERPRETATION.md'))
        (dest/'notebook_update_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
        print(json.dumps(record,indent=2),flush=True)

if __name__=='__main__':run()
