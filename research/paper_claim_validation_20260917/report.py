"""Extend the existing notebook only; retain all earlier GBM/PTC cells and outputs."""
import fcntl
import json
import os
from pathlib import Path
import shutil
from common import ROOT, OUT, FEATURES, ROUTES, require_slurm, sha, write_json, checked, utc
TAG='paper_claim_validation_20260917'

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import nbformat
    from nbclient import NotebookClient
    import workflow
    import diagnostics
    import plot_unit
    d=OUT/'summary'
    assert checked(d,'aggregate_manifest.json','AGGREGATE_COMPLETE')
    audit=json.loads((OUT/'protocol/input_audit.json').read_text())
    agg=json.loads((d/'aggregate_manifest.json').read_text())
    assert agg['completed_units']==726 and agg['n_clustering_results']==2904
    for budget in FEATURES:
        plot_unit.run(OUT/'GBM'/audit['display_sample']/budget,refresh_marker_legend=True)
    diagnostics.run()
    means=pd.read_csv(d/'patient_aggregates.csv',dtype={'cutoff':str})
    fixed=means[(means.cohort=='primary97')&(means.library=='CM2_glioma_other')&(means.cutoff=='mean')&(means.stage=='terminal090')&(means.family=='native_R_budget')]
    assert len(fixed)==24
    fixed.to_csv(d/'primary_fixed_marker_24.csv',index=False)
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained')
    colors=['#0072B2','#D55E00','#009E73','#CC79A7']
    for route,color in zip(ROUTES,colors):
        g=fixed[fixed.route==route].set_index('budget').loc[FEATURES]
        ax.plot(range(6),g.macroF1_present_mean,marker='o',color=color,label=route)
    ax.set_xticks(range(6),['500','1,000','2,000','3,000','5,000','All'])
    ax.set(xlabel='Native R feature budget',ylabel='Patient-weighted terminal L1 macro-F1',
        title='GSE274546: fixed CM2_glioma_other / mean cutoff')
    ax.legend(frameon=False,fontsize=8);ax.spines[['top','right']].set_visible(False)
    for ext in ['png','pdf']:fig.savefig(d/f'fixed_marker_HVG_routes.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)
    paired=pd.read_csv(d/'fixed_marker_patient_paired.csv')
    a=paired[paired.cohort=='primary97'].sort_values(['budget','route']).copy()
    a=a[~((a.budget=='hvg2000')&(a.route=='UMAP2_HDBSCAN_R'))]
    fig,ax=plt.subplots(figsize=(9,8),layout='constrained')
    ax.errorbar(a.mean_delta,range(len(a)),xerr=np.vstack([a.mean_delta-a.CI95_low,a.CI95_high-a.mean_delta]),fmt='o',ms=4,color='#0072B2')
    ax.axvline(0,color='#555',ls='--',lw=1)
    ax.set_yticks(range(len(a)),a.budget+' / '+a.route,fontsize=8)
    ax.set(xlabel='Terminal macro-F1 difference versus HVG2000 / UMAP2 / HDBSCAN',
        title='Paired patients; 95% bootstrap intervals; fixed marker and cutoff')
    ax.spines[['top','right']].set_visible(False)
    for ext in ['png','pdf']:fig.savefig(d/f'paired_HVG_route_effects.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)
    cv=pd.read_csv(d/'patient_heldout_test_results.csv',dtype={'cutoff':str})
    cv_summary=cv.groupby(['cohort','evidence','selection'])[['macroF1_present','coverage','unknown_rate']].agg(['mean','std','count'])
    cv_summary.columns=['_'.join(c) for c in cv_summary.columns]
    cv_summary.reset_index().to_csv(d/'patient_heldout_summary.csv',index=False)
    factorial=pd.read_csv(d/'marker_DL_factorial_heldout.csv')
    fac=factorial.groupby(['cohort','marker_selection','stage'])[['macroF1_present','coverage']].agg(['mean','std','count'])
    fac.columns=['_'.join(c) for c in fac.columns];fac.reset_index().to_csv(d/'marker_DL_factorial_summary.csv',index=False)
    inventory=[]
    for sample in (OUT/'protocol/cohort.csv').read_text().splitlines()[1:]:
        sid=sample.split(',')[0]
        for budget in FEATURES:
            prep=OUT/'GBM'/sid/budget
            assert checked(prep/'figures','manifest.json','FIGURES_COMPLETE')
            inventory.append(dict(sample=sid,budget=budget,routes=';'.join(ROUTES),
                atlas_png=str((prep/'figures/all_cluster_routes.png').relative_to(ROOT)),
                atlas_pdf=str((prep/'figures/all_cluster_routes.pdf').relative_to(ROOT)),
                metrics=str((prep/'evaluation/metrics.csv').relative_to(ROOT))))
    pd.DataFrame(inventory).to_csv(d/'all_clustering_figure_index.csv',index=False)
    workflow.run()
    rank=fixed.sort_values('macroF1_present_mean',ascending=False)
    best=rank.iloc[0];original=fixed[(fixed.budget=='hvg2000')&(fixed.route=='UMAP2_HDBSCAN_R')].iloc[0]
    status=pd.read_csv(d/'terminal_status_counts.csv')
    report=f'''# Native-R GSE274546 core experiment

Completed: 121 samples / 59 patients, 429,305 author-filtered cells; the historical primary cohort is 97 samples / 55 patients. All 726 sample-feature units and 2,904 clustering outputs are included. Each route has 16 frozen marker libraries and three density cutoffs, giving 139,392 terminal conditions. These include valid no-op and untrainable states; they are not 139,392 successful neural-network training runs.

The original RNA branch uses native VST genes for PCA and normalized DL input, while DEG scoring keeps all eligible RNA genes. There is no batch correction within a single GBM sample. This differs from the PTC integrated-assay branch where changing the integration feature budget can also change scoring genes. The original10-epoch MLP preserves known native-panel labels and fills Undecided only. A 0.90 threshold is primary; 0.70 is a separately reported sensitivity.

At fixed CM2_glioma_other / mean cutoff, the highest primary-cohort point estimate is {best.budget}/{best.route}: patient-weighted terminal L1 macro-F1 {best.macroF1_present_mean:.6f}. The original HVG2000/UMAP2/HDBSCAN configuration is {original.macroF1_present_mean:.6f}. Paired patient effects and intervals, rather than the highest point alone, determine the strength of this comparison. No equivalence or noninferiority margin was prespecified.

Configuration selection uses five frozen patient folds. Training-patient labels choose HVG/route/marker/cutoff; heldout labels only score. The cohort was previously explored, so this is retrospective label-heldout validation, not a never-seen cohort. Every target sample is fitted independently. Author L1 macro-F1 uses all cells; Unknown and unmapped predictions count as errors. Fixed11 and collapsed-neuron10 are secondary, with no outcome-driven choice of label granularity. Malignant state ARI is a clustering evaluation, not final subtype annotation.

CARE_TME and BrainAtlas112 participated in reference-label construction; UNION_all includes them. The primary marker-selection table excludes these three, while the complete concordance table retains them. The fixed glioma marker vocabulary lacks several author classes; missing classes remain in the denominator.

All clustering configurations have PNG/PDF atlases in all_clustering_figure_index.csv. Each sample uses fixed HVG2000 UMAP display coordinates for every atlas. Display coordinates never substitute for the actual clustering input. The median-size primary sample {audit['display_sample']} was chosen before model scoring; it also has all-marker before/after-DL figures.

This report completes the GBM core stage. Geometry-only and parameter controls, PTC repeat/marker-retention analyses, fair comparator benchmarks and additional reviewer resource evidence have separate completion gates. It does not assert those follow-ups are complete. Internal workflow comparisons do not establish superiority to other annotation methods. Historical NMT selected PCA-SNN and remains an explicit exception.
'''
    (d/'GBM_CORE_REPORT.md').write_text(report)
    (d/'GBM_CORE_REPORT_ZH.md').write_text(f'''# GSE274546 原 R 核心实验\n\n已完成121样本/59患者、429305细胞；历史主队列97样本/55患者保留。六档特征×四聚类路线×16库×3cutoff均到终端DL，实际训练与no-op分开记录。\n\n固定旧marker/mean时，最高点估计为{best.budget}/{best.route}，患者平均终端macro-F1={best.macroF1_present_mean:.6f}；原HVG2000/UMAP2/HDBSCAN={original.macroF1_present_mean:.6f}。完整配对差与区间见fixed_marker_patient_paired.csv，不根据单个最大值宣称全局最优。\n\n新结果追加至原notebook；所有2904个聚类结果均有图。主图样本{audit['display_sample']}按细胞数中位数事先选定。标签、Unknown及marker缺类口径固定；患者留出为回顾性验证。PTC历史NMT-SNN分支不删除。后续geometry-only、竞争方法、PTC补实验和规模验证单独验收，本报告不将它们写成已完成。\n''')
    path=ROOT/'notebooks/dgscrna_results.ipynb'
    with (OUT/'notebook.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        before=sha(path);nb=nbformat.read(path,as_version=4)
        baseline=[c for c in nb.cells if TAG not in c.get('metadata',{}).get('tags',[])]
        backup=OUT/'notebook_before_claim_validation.ipynb'
        if not backup.exists():shutil.copy2(path,backup)
        cells=[]
        def md(s):cells.append(nbformat.v4.new_markdown_cell(s,metadata={'tags':[TAG]}))
        def code(s):cells.append(nbformat.v4.new_code_cell(s,metadata={'tags':[TAG]}))
        md('<a id="paper-claim-validation-20260917"></a>\n# GSE274546: native-R workflow validation\n\n'+report.split('\n',1)[1])
        code('''from pathlib import Path
import pandas as pd
from IPython.display import display, Image, Markdown
C=Path('/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/paper_claim_validation_20260917')
def table(name): return pd.read_csv(C/'summary'/name)
def figure(rel): display(Image(filename=str(C/rel)))
display(pd.read_csv(C/'protocol/cohort.csv'))
display(pd.read_csv(C/'verification/pilot_checks.csv'))
display(pd.read_csv(C/'markers/coverage_vocabulary.csv'))''')
        md('## Full decision tree and fixed-marker HVG comparison\n\nThe native-R single-sample scoring universe is all eligible RNA genes. Actual geometry and DL feature lists are saved per unit.')
        code('figure("summary/workflow_decision_tree.png")\nfigure("summary/fixed_marker_HVG_routes.png")\nfigure("summary/paired_HVG_route_effects.png")\ndisplay(table("primary_fixed_marker_24.csv"))\ndisplay(table("fixed_marker_patient_paired.csv"))')
        md('## Patient-heldout selection and marker × DL contribution\n\nThe same selected marker library is used on both sides of each no-DL/with-DL pair. Fixed-cutoff factorial and broader configuration selection are separate tables.')
        code('display(table("patient_heldout_summary.csv"))\ndisplay(table("patient_heldout_selection.csv"))\ndisplay(table("marker_DL_factorial_summary.csv"))\ndisplay(table("terminal_status_counts.csv"))')
        md('## Marker-source distributions at the fixed original geometry\n\nEvery dot is one primary-cohort patient after averaging that patient\'s samples. All 16 libraries retain the same HVG2000/UMAP2/HDBSCAN geometry and mean cutoff. Orange libraries overlap author-label construction and remain concordance results; they are excluded from primary marker selection.')
        code('figure("summary/marker_context_distributions.png")\ndisplay(table("marker_context_summary.csv"))')
        md('## Per-class errors, Unknown and actual refinement changes\n\nAll cells remain in each confusion-matrix denominator. Wrong known marker seeds are preserved by the historical algorithm; the audit separately records correct and incorrect newly filled cells.')
        code('figure("summary/per_class_original2000.png")\nfigure("summary/display_sample_confusions.png")\ndisplay(table("original2000_fixed_marker_per_class_mean.csv"))\ndisplay(table("fixed_marker_DL_change_audit.csv"))')
        md('## Prespecified display sample: all six budgets and all four clusterers')
        for budget in FEATURES:
            md(f'### {audit["display_sample"]}: {budget}')
            code(f'figure("GBM/{audit["display_sample"]}/{budget}/figures/all_cluster_routes.png")')
            links='\n'.join(f'- [{route}: all marker-only and terminal labels](../results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/GBM/{audit["display_sample"]}/{budget}/figures/{route}_all_markers.pdf)' for route in ROUTES)
            md(links)
        md('## Every sample: fixed-marker table and clustering atlas\n\nAll six feature budgets have four-route PNG/PDF figures. The original-2000 atlas is embedded below for every sample; adjacent links expose every other budget. No sample is removed because its score is poor.')
        code('ALL=pd.read_csv(C/"summary/all_annotation_metrics.csv.gz",dtype={"cutoff":str})\nFIXED=ALL[(ALL.library=="CM2_glioma_other")&(ALL.cutoff=="mean")&(ALL.stage=="terminal090")&(ALL.family=="native_R_budget")]')
        cohort=pd.read_csv(OUT/'protocol/cohort.csv')
        for row in cohort.itertuples():
            links=' | '.join(f'[{b} PNG](../results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/GBM/{row.sample}/{b}/figures/all_cluster_routes.png)' for b in FEATURES)
            md(f'### {row.sample} — {row.patient}\n\n{row.n_cells:,} cells; primary cohort: {row.primary}.\n\n{links}')
            code(f'display(FIXED[FIXED["sample"]=={row.sample!r}][["budget","route","macroF1_present","coverage","dl_status"]])\nfigure("GBM/{row.sample}/hvg2000/figures/all_cluster_routes.png")')
        campaign=nbformat.v4.new_notebook(cells=cells,metadata={'kernelspec':{'display_name':'Python 3','language':'python','name':'python3'}})
        NotebookClient(campaign,timeout=1200,kernel_name='python3',resources={'metadata':{'path':str(ROOT/'notebooks')}},allow_errors=False).execute()
        assert sha(path)==before,'Notebook changed externally during update'
        banner=nbformat.v4.new_markdown_cell('**原 R 的 GSE274546 核心实验已完成。** [六档HVG、完整聚类图、患者留出与marker/DL对照](#paper-claim-validation-20260917)追加在文末；原有GBM/PTC内容全部保留。',metadata={'tags':[TAG]})
        nb.cells=[banner]+baseline+campaign.cells
        tmp=path.with_suffix('.ipynb.claim_part');nbformat.write(nb,tmp);tmp.replace(path)
        current=nbformat.read(path,as_version=4)
        assert current.cells[1:1+len(baseline)]==baseline
        errors=sum(o.get('output_type')=='error' for c in campaign.cells for o in c.get('outputs',[]))
        assert errors==0
        write_json(d/'notebook_update_manifest.json',dict(status='completed_GBM_core',original_cells_preserved=len(baseline),
            added_cells=len(campaign.cells)+1,errors=errors,previous_sha256=before,current_sha256=sha(path),
            job=os.environ['SLURM_JOB_ID'],completed_at=utc(),source_sha256=sha(__file__)))
    print('GBM_CORE_NOTEBOOK_UPDATED',sha(path),flush=True)

if __name__=='__main__':run()
