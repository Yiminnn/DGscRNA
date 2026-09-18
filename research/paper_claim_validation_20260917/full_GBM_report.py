"""Extend the original notebook with completed GBM controls and fair comparisons."""
import fcntl
import json
import os
import shutil
from common import ROOT,OUT,require_slurm,checked,complete,sha,write_json,utc
TAG='paper_claim_validation_20260917_followups'

def run():
    require_slurm()
    import pandas as pd
    import nbformat
    from nbclient import NotebookClient
    import workflow
    assert (OUT/'summary/GBM_CORE_DELIVERED.json').exists()
    folders=['controls_summary','comparison_summary','comparison_summary/diagnostics',
             'unknown_summary','marker_evidence_summary','scalability_summary']
    for folder in folders:assert checked(OUT/folder),folder
    dest=OUT/'GBM_full_summary';dest.mkdir(exist_ok=True)
    workflow.run()
    core=pd.read_csv(OUT/'summary/primary_fixed_marker_24.csv')
    best=core.sort_values('macroF1_present_mean',ascending=False).iloc[0]
    original=core[(core.budget=='hvg2000')&(core.route=='UMAP2_HDBSCAN_R')].iloc[0]
    held=pd.read_csv(OUT/'comparison_summary/patient_heldout_summary.csv')
    primary=held[held.cohort=='primary97'].set_index('method')
    marker=primary.loc[['DG-scRNA','scType','scCATCH','SCINA']].sort_values('macroF1_present_mean',ascending=False)
    text=f'''# Completed GBM evidence package

The native-R core covers 121 samples / 59 patients and 429,305 cells, with the historical primary 97 samples / 55 patients reported separately. All 726 feature-budget units and 2,904 clustering outputs have terminal predictions, evaluation tables and figures. The 139,392 terminal conditions include trained, no-op and untrainable states; these are not 139,392 successful neural-network fits.

At fixed historical glioma marker/mean cutoff, HVG2000/UMAP2/HDBSCAN has primary patient-weighted terminal macro-F1 {original.macroF1_present_mean:.6f}. The highest point among the 24 fixed-marker configurations is {best.budget}/{best.route}, {best.macroF1_present_mean:.6f}. Paired-patient intervals and Holm-adjusted comparisons determine the strength of evidence. Selecting a maximum is not proof of a global optimum, equivalence or noninferiority.

All 363 geometry-only controls hold RNA scoring and normalized HVG2000 DL expression fixed while changing geometry. Their 2000-gene anchor exactly equals the original arm. Native versus fixed DL genes at the same geometry isolates the DL-input effect. Another 54 MLP controls and 180 representation controls cover the three count-selected size pilots. Model seeds and embedding seeds are varied separately; the original default is retained. These repeats describe algorithm sensitivity, not additional patients or all-cohort optimal dimensions.

The matched-partition patient-heldout comparison uses terminal DG-scRNA, scType and scCATCH on the same HVG2000/UMAP2/HDBSCAN partition, plus per-cell SCINA with the same marker roster. Each receives its recorded threshold opportunities. CARE_TME, BrainAtlas112 and UNION_all are excluded from primary marker selection because they overlap author-label construction. The highest heldout point among marker-information methods is {marker.index[0]} ({marker.iloc[0].macroF1_present_mean:.6f}); DG-scRNA is {primary.loc['DG-scRNA','macroF1_present_mean']:.6f}. All competitors remain in the tables. Tuned 24-configuration DG/scType results are not substituted into this fixed-partition comparison.

SingleR references contain only labelled training-patient cells, excluding all samples of each heldout patient; its primary mean is {primary.loc['SingleR','macroF1_present_mean']:.6f}. Published pretrained Brain scDeepSort is an atlas-GNN condition, with mean {primary.loc['scDeepSort','macroF1_present_mean']:.6f}. These reference-information conditions differ from marker methods. The pretrained atlas lacks a malignant class; missing vocabulary remains in the denominator. Every primary score includes all cells, including Unknown/unmapped calls. These are retrospective patient-label-heldout folds in an explored cohort, not a never-seen external test cohort.

Unknown is analyzed using counts, detected genes and visible mitochondrial percentage within author cell class. These are descriptive QC associations; doublets and novel cell types are not inferred. The original DL preserves known marker seeds; retained seed errors and newly filled wrong predictions are reported separately. Marker tables expose recovered source/species/assay records and unresolved links without manufacturing independence.

The resource experiment contains 15 cold-input runs: three methods on the same nested 10k / 30k / 50k / 100k / 120k pooled-cell inputs. Each point is one workflow run. Actual wall time, hardware and scheduler peak RSS are recorded. This is one timing measurement per point, not repeated runtime inference. Pooled input is not biological batch-correction evidence; see RESOURCE_INTERPRETATION.md for preprocessing and memory-measurement limits.

The original notebook is extended, with all clustering figures indexed and the English decision tree updated. This completes the GBM evidence stage only. Historical PTC NMT-SNN and TTU-UMAP/HDBSCAN anchors remain separate; new PTC repeats and marker-retention controls follow delivery. No new Darmanis or independent Pu cohort is included.
'''
    (dest/'GBM_FULL_REPORT.md').write_text(text)
    (dest/'GBM_FULL_REPORT_ZH.md').write_text(f'''# GBM 完整补实验

原R主网格覆盖121样本/59患者/429305细胞，97样本/55患者主队列另报。固定旧marker/mean时，原HVG2000/UMAP2/HDBSCAN患者平均终端macro-F1={original.macroF1_present_mean:.6f}；24配置最高点为{best.budget}/{best.route}，{best.macroF1_present_mean:.6f}。需结合配对区间与多重比较，不能把最大值写成全局最优。

363个geometry-only、54个MLP和180个表示/聚类控制完成。参数与种子控制仅限事先按细胞数选定的3个样本。固定原聚类条件下，marker方法留出患者最高点为{marker.index[0]}（{marker.iloc[0].macroF1_present_mean:.6f}），DG-scRNA={primary.loc['DG-scRNA','macroF1_present_mean']:.6f}。SingleR和scDeepSort使用不同reference信息，结果单列保留。

单次规模验证为三方法×五档细胞量，最大120000细胞；不是多个小样本合计。Unknown、marker来源和DL新增错误表已更新。沿用原notebook和OneDrive目录；PTC后续控制单独推进，本文件不宣称整项研究完成。
''')
    methods='''# Method description corrections supported by executed source

- Describe DG-scRNA as context-aware marker annotation with expression-MLP refinement. The terminal model receives normalized expression, not graph edges or learned graph embeddings. SNN/HDBSCAN supplies upstream partitions.
- Single-sample native R uses LogNormalize 10,000, VST, ScaleData(max 10), PCA30 and uwot UMAP(cosine, 30 neighbors,min_dist0.3), followed by R HDBSCAN(minPts 50) or SNN/Louvain(resolution 0.5). RNA DEG scoring retains all eligible genes; DL uses normalized selected genes.
- Density sums original DEG log2FC above 1 divided by full panel length, with singleton factor 0.8. Ties are Undecided. The implementation has no alpha/Jaccard mixture or epsilon iteration.
- The preserved legacy MLP is 256/128 LeakyReLU, Softmax passed to CrossEntropyLoss, Adamax 0.001, batch 256, 10 epochs and the archived 90/10 split. The unusual Softmax/loss combination is retained for reproduction. Known labels are preserved; only Undecided cells can be filled.
- Confidence 0.90 is the implemented primary threshold. 0.70 is a sensitivity using the same saved probabilities.
- Annotation macro-F1 and partition ARI/NMI are distinct. State the all-cell denominator and Unknown handling. PTC keeps historical non-T-positive F1 and hard-call binary AUC definitions; unresolved historical Accuracy remains a provenance issue.
- Patient folds are retrospective. Bootstrap intervals condition on frozen heldout predictions; no formal noninferiority margin was specified.

These are proposed wording corrections stored with the experiment; the manuscript itself is preserved.
'''
    (dest/'METHODS_CORRECTIONS.md').write_text(methods)
    path=ROOT/'notebooks/dgscrna_results.ipynb'
    with (OUT/'notebook.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        before=sha(path);nb=nbformat.read(path,as_version=4)
        baseline=[c for c in nb.cells if TAG not in c.get('metadata',{}).get('tags',[])]
        backup=OUT/'notebook_before_GBM_followups.ipynb'
        if not backup.exists():shutil.copy2(path,backup)
        cells=[]
        def md(s):cells.append(nbformat.v4.new_markdown_cell(s,metadata={'tags':[TAG]}))
        def code(s):cells.append(nbformat.v4.new_code_cell(s,metadata={'tags':[TAG]}))
        md('<a id="gbm-native-r-followups"></a>\n# GBM: completed controls and fair comparisons\n\n'+text.split('\n',1)[1])
        code('''from pathlib import Path
import pandas as pd
from IPython.display import display, Image, Markdown
F=Path('/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/paper_claim_validation_20260917')
def extra_table(rel): display(pd.read_csv(F/rel))
def extra_figure(rel): display(Image(filename=str(F/rel)))
extra_figure('summary/workflow_decision_tree.png')''')
        md('## Geometry and DL input have different effects\n\nFull RNA scoring is fixed. The 2000 anchor is checked exactly; patient contrasts separate geometry from DL-input changes.')
        code("extra_figure('controls_summary/geometry_and_DL.png')\nextra_table('controls_summary/geometry_vs_DL_paired_effects.csv')")
        md('## Parameter and seed sensitivity\n\nThree prespecified size pilots. MLP and embedding seeds are separate; neither is a patient replicate. Every representation condition has a clustering/terminal PNG and PDF.')
        code("extra_figure('controls_summary/MLP_controls.png')\nextra_figure('controls_summary/representation_controls.png')\nextra_table('controls_summary/representation_figure_index.csv')")
        index=pd.read_csv(OUT/'controls_summary/representation_figure_index.csv')
        for (sample,budget),group in index.groupby(['sample','budget'],sort=False):
            links='\n'.join(f'- [{r.name}](../results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/{r.pdf})' for r in group.itertuples())
            md(f'### {sample} / {budget}: every representation-clustering result\n\n'+links)
        md('## Matched marker opportunities and distinct reference conditions\n\nAll cluster methods use the original partition. Training patients select marker and threshold. Full 24-configuration tuning is a separate experiment.')
        code("extra_figure('comparison_summary/diagnostics/heldout_method_comparison.png')\nextra_figure('comparison_summary/diagnostics/heldout_class_F1.png')\nextra_figure('comparison_summary/diagnostics/display_sample_methods.png')\nextra_table('comparison_summary/patient_heldout_summary.csv')\nextra_table('comparison_summary/paired_patient_comparisons.csv')\nextra_table('comparison_summary/training_patient_choices.csv')\nextra_table('comparison_summary/SCINA_solver_counts.csv')")
        md('### Equal tuning budget: DG-scRNA and scType\n\nBoth select among 24 geometries and 13 primary marker libraries with 3 thresholds on training patients. This secondary table does not include the fixed-partition scCATCH result. SCINA native numerical failures use independently audited boundary guards; every solver state and dropped signature count is recorded.')
        code("extra_table('comparison_summary/equal24_DG_scType_comparison.csv')\nextra_table('comparison_summary/equal24_DG_scType_choices.csv')")
        md('## Unknown, seed errors and refinement outcomes\n\nQC associations condition on original class; they do not diagnose doublets or novel types.')
        code("extra_figure('unknown_summary/Unknown_QC_association.png')\nextra_table('unknown_summary/retained_seed_and_DL_new_errors.csv')\nextra_table('marker_evidence_summary/library_source_assay_coverage_audit.csv')\ndisplay(Markdown((F/'marker_evidence_summary/MARKER_EVIDENCE.md').read_text()))")
        md('## Single-run resource measurements, 10k through 120k cells\n\nFrozen nested pooled counts; resource evidence, not new batch-correction or biological validation.')
        code("extra_figure('scalability_summary/resource_curves.png')\nextra_table('scalability_summary/cold_pipeline_resource_curve.csv')\ndisplay(Markdown((F/'scalability_summary/RESOURCE_INTERPRETATION.md').read_text()))")
        for n in [10000,30000,50000,100000,120000]:code(f"extra_figure('scalability/{n}/clustering.png')")
        md('## Wording supported by the implementation\n\n'+methods.split('\n',1)[1])
        addition=nbformat.v4.new_notebook(cells=cells,metadata={'kernelspec':{'display_name':'Python3','language':'python','name':'python3'}})
        NotebookClient(addition,timeout=1200,kernel_name='python3',resources={'metadata':{'path':str(ROOT/'notebooks')}},allow_errors=False).execute()
        assert sha(path)==before,'Notebook changed externally'
        banner=nbformat.v4.new_markdown_cell('**GBM原R主网格、机制控制、公平比较与单次规模实验已完成。** [完整补实验](#gbm-native-r-followups)追加在原notebook；PTC后续单独推进。',metadata={'tags':[TAG]})
        nb.cells=[banner]+baseline+addition.cells
        tmp=path.with_suffix('.ipynb.full_GBM_part');nbformat.write(nb,tmp);tmp.replace(path)
        current=nbformat.read(path,as_version=4);assert current.cells[1:1+len(baseline)]==baseline
        errors=sum(o.get('output_type')=='error' for c in addition.cells for o in c.get('outputs',[]));assert errors==0
        write_json(dest/'notebook_manifest.json',dict(status='completed',old_cells_preserved=len(baseline),added_cells=len(addition.cells)+1,
            errors=errors,previous_sha256=before,current_sha256=sha(path)))
    write_json(dest/'manifest.json',dict(status='completed_GBM_PTC_followups_pending',stage_manifests={v:sha(OUT/v/'manifest.json') for v in folders},
        files={p.name:sha(p) for p in dest.iterdir() if p.suffix in ['.md','.json'] and p.name!='manifest.json'},
        notebook_sha256=sha(path),job=os.environ['SLURM_JOB_ID'],source_sha256=sha(__file__),completed_at=utc()))
    complete(dest)

if __name__=='__main__':run()
