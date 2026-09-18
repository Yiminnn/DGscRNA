"""Append completed PTC controls to the canonical notebook, preserving old cells."""
import fcntl
import json
import os
from pathlib import Path
import shutil
from common import ROOT,OUT,sha,checked,write_json,complete,utc
from ptc_followup_common import PTC,ANCHORS,require_ptc

TAG='paper_claim_validation_20260917_PTC_followups'

def run():
    require_ptc()
    import pandas as pd
    import nbformat
    from nbclient import NotebookClient
    import workflow
    dest=OUT/'PTC_summary';assert checked(dest)
    workflow.run()
    anchors=pd.read_csv(dest/'original_anchor_patient_metrics.csv')
    tuned=pd.read_csv(dest/'heldout_workflow_patient_metrics.csv')
    paper=pd.read_csv(dest/'original_paper_endpoint_reconciliation.csv')
    core=pd.read_csv(dest/'fresh24_patient_summary.csv')
    rows=[]
    for group,a in ANCHORS.items():
        p=paper[(paper.group==group)&(paper.patient=='ALL')&(paper.stage=='final090')].iloc[0]
        b=anchors[anchors.group==group].F1_T.mean();t=tuned[tuned.group==group].F1_T.mean()
        candidates=core[(core.group==group)&(core.marker_choice=='training_patient_selected_context')]
        best=candidates.sort_values('F1_T',ascending=False,kind='stable').iloc[0]
        rows.append(dict(group=group,historical_route=a['route'],historical_library=a['library'],historical_cutoff=a['cutoff'],
            reconstructed_paper_nonT_F1=float(p.F1_nonT),reconstructed_paper_T_F1=float(p.F1_T),
            reconstructed_paper_hardcall_AUC=float(p.AUC_binary),
            original_anchor_patient_mean_productive_TCR_F1=float(b),heldout_workflow_patient_mean_productive_TCR_F1=float(t),
            best_observed_fresh_budget=str(best.budget),best_observed_fresh_route=best.route,
            best_observed_fresh_TCR_F1=float(best.F1_T)))
    result=pd.DataFrame(rows);result.to_csv(dest/'claim_summary_by_group.csv',index=False)
    details='\n'.join(f"- {r.group}: archived {r.historical_route}, {r.historical_library}/{r.historical_cutoff}; "
        f"reconstructed paper non-T F1 {r.reconstructed_paper_nonT_F1:.6f}, T F1 {r.reconstructed_paper_T_F1:.6f}, "
        f"hard-call AUC {r.reconstructed_paper_hardcall_AUC:.6f}. Four-patient mean strict-T/productive-TCR F1 "
        f"is {r.original_anchor_patient_mean_productive_TCR_F1:.6f}; training-patient-selected fresh workflow is "
        f"{r.heldout_workflow_patient_mean_productive_TCR_F1:.6f}." for r in result.itertuples())
    text='''# PTC: original anchors, patient-label holdout and controlled follow-ups

The completed PTC native-R grid is reused, and missing seed and marker-retention
controls are added. The archived all-eight-sample CCA2000 fit and new within-group
NMT/TTU fits remain separate. All original 17 marker libraries and three cutoffs
are retained. The selected original NMT SNN and TTU UMAP-HDBSCAN anchors are not
replaced by a universal UMAP-HDBSCAN story.

'''+details+'''

These endpoints have different meanings. Paper non-T-positive F1 and hard-call
binary AUC preserve the original evaluation. Strict-T F1 against productive,
high-confidence TCR detection is an assay-concordance endpoint, with TCR-positive
recall and detection yield reported separately. Undetected TCR does not establish
non-T identity. S2/S3 native-label agreement is concordance. The old reported
Accuracy row still lacks its original calculation source; it is not silently
replaced by ordinary accuracy and does not block the accepted approximate replay.

There are four patients, each with one NMT and one TTU sample. Both samples are
excluded from label-based marker/configuration selection when that patient is
held out. The primary selection objective was frozen as mean training-patient
strict-T/productive-TCR F1. Expression integration and pseudo-label refinement
remain transductive within the full group. This is retrospective patient-label
holdout in a previously explored dataset, not an unseen-patient expression test.

The marker x DL analysis uses one selected context for initial and terminal calls,
so DL contribution is not confounded by independently choosing the best marker for
each stage. A stage-specific selection sensitivity is named separately. All four
paired patient differences, 256 bootstrap resamples and 16 sign flips are retained.
The smallest attainable two-sided sign-flip p-value is 0.125. Intervals condition
on the saved predictions and do not treat cells or seeds as independent patients.

Model-seed controls vary initialization only (0,1,2,3,42), retaining split42.
Representation controls vary PCA/UMAP seeds and rerun clustering, with fixed CCA
expression, feature identities, SNN seed0 and MLP seed42. They do not measure
integration-seed variation. Default controls reproduce saved terminal labels,
splits and rounded confidence; numerical probability tolerances are recorded.
Historical trained weights have not been recovered.

The marker-retention experiment preserves the exact original geometry and cluster
partitions while using a single scoring universe per group: fixed2000 union all
17-library marker genes, intersected with eligible CCAall genes. Geometry2000,
geometry5000 and geometryall share that score matrix and the same DL2000 matrix.
All 17x3 contexts reach terminal DL. This tests the missing-marker mechanism and
geometry conditional on a fixed CCAall fit; it is not a substitute for original
CCA2000. No isolated CD3D insertion or outcome-based label editing is performed.

All 22 new control units have both clustering routes plotted (44 clustering
conditions); all 50 MLP task groups have saved terminal outputs and evaluation.
The pre-existing native-R grid and batch-biology comparisons are retained.
Historical competitor outputs keep their original input/reference conditions;
the matched GBM comparison does not retroactively make these PTC comparisons fair.
No new Darmanis or independent Pu validation cohort was run.

A publication claim must name the evaluated cohorts, metric, marker-selection
policy and search space. An observed maximum is not proof of a global optimum;
retaining NMT exceptions, endpoint tradeoffs and Unknown coverage is part of the
evidence. Use the actual tables below to distinguish a fixed original configuration
from a tuned workflow and from a merely higher point estimate.
'''
    (dest/'PTC_FOLLOWUP_REPORT.md').write_text(text)
    (dest/'PTC_FOLLOWUP_REPORT_ZH.md').write_text('''# PTC 补充实验完成

保留原八样本 CCA2000 与新 NMT/TTU 分组整合，保留 NMT 的 SNN 和 TTU 的 UMAP-HDBSCAN 原始分支。新增内容为四患者留出 marker/流程选择、marker×DL 四方、分别改变 MLP 与 PCA/UMAP 的五种子对照，以及统一 marker-retention 机制对照。复用已完成的30个原R单元；新增22个控制单元（44个聚类条件）和50组MLP任务。

论文原 non-T F1/AUC、严格T/TCR检出一致性、TCR recall 和 Unknown coverage 分开报告。TCR未检出不等于确定非T；S2/S3一致率不等于独立准确率。4个患者不能由细胞数或种子数扩充，精确双侧符号翻转检验最小p值为0.125。缺失历史权重及Accuracy原始计算来源继续保留记录，不妨碍用户接受的近似复现。

实际分组数值见 claim_summary_by_group.csv。最高点、训练患者选出的流程和固定原配置是不同结论；不能在存在NMT例外时写成所有步骤在两组始终最优。
''')
    notebook=ROOT/'notebooks/dgscrna_results.ipynb'
    with (OUT/'notebook.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        before=sha(notebook);nb=nbformat.read(notebook,as_version=4)
        baseline=[c for c in nb.cells if TAG not in c.get('metadata',{}).get('tags',[])]
        backup=PTC/'notebook_before_PTC_followups.ipynb'
        if not backup.exists():shutil.copy2(notebook,backup)
        cells=[]
        def md(s):cells.append(nbformat.v4.new_markdown_cell(s,metadata={'tags':[TAG]}))
        def code(s):cells.append(nbformat.v4.new_code_cell(s,metadata={'tags':[TAG]}))
        md('<a id="ptc-r-followups"></a>\n'+text)
        code('''from pathlib import Path
import pandas as pd
from IPython.display import display, Image
P=Path('/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/paper_claim_validation_20260917')
def ptc_table(name): display(pd.read_csv(P/'PTC_summary'/name))
def ptc_figure(name): display(Image(filename=str(P/name)))
ptc_figure('summary/workflow_decision_tree.png')
ptc_table('claim_summary_by_group.csv')''')
        md('## Fresh group workflows and patient-label holdout\n\nFixed original markers and training-patient marker choices receive separate curves. All tests retain Unknown cells.')
        code("ptc_figure('PTC_summary/PTC_workflow_comparison.png')\nptc_table('fresh24_patient_summary.csv')\nptc_table('paired_patient_contrasts.csv')\nptc_table('paired_patient_differences.csv')")
        md('## Marker and DL contributions\n\nThe main four-way analysis uses the same terminal-selected marker in both stages; patient differences and the interaction are explicit.')
        code("ptc_table('marker_DL_paired_effects.csv')\nptc_table('marker_DL_patient_differences.csv')")
        md('## Initialization versus representation stability\n\nSame five seeds, separate perturbations. Each point aggregates the same four patients.')
        code("ptc_figure('PTC_summary/PTC_MLP_seed_stability.png')\nptc_figure('PTC_summary/PTC_representation_seed_stability.png')\nptc_table('seed_patient_mean_summary.csv')")
        md('## Uniform marker-retention mechanism\n\nSame geometry, partitions and DL input as the earlier geometry-only controls. Only the scoring gene universe changes under one rule shared by all budgets.')
        code("ptc_figure('PTC_summary/PTC_marker_retention_mechanism.png')\nptc_table('marker_retention_matched_summary.csv')")
        md('## Original paper endpoint reconciliation\n\nOriginal non-T-positive F1 and hard-call AUC remain separate from strict-T productive-TCR selection. No historical Accuracy value is fabricated.')
        code("ptc_table('original_paper_endpoint_reconciliation.csv')")
        md('## Every new clustering result\n\nFixed group display coordinates; saved lineage and TCR are evaluation overlays only. Each panel shows the fixed historical group marker before and after terminal DL.')
        figures=pd.read_csv(dest/'all_new_clustering_figures.csv')
        for row in figures.itertuples():
            md(f"### {row.name} / {row.route}\n\n[PDF](../results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/{row.pdf})")
            code(f'ptc_figure({row.png!r})')
        addition=nbformat.v4.new_notebook(cells=cells,metadata={'kernelspec':{'display_name':'Python3','language':'python','name':'python3'}})
        NotebookClient(addition,timeout=1200,kernel_name='python3',resources={'metadata':{'path':str(ROOT/'notebooks')}},allow_errors=False).execute()
        assert sha(notebook)==before,'Original notebook changed externally'
        banner=nbformat.v4.new_markdown_cell('**GBM及PTC本轮补充实验已完成。** [GBM](#gbm-native-r-followups) / [PTC](#ptc-r-followups)；均追加于原notebook，旧细胞和历史结果保留。',metadata={'tags':[TAG]})
        nb.cells=[banner]+baseline+addition.cells
        tmp=notebook.with_suffix('.ipynb.PTC_part');nbformat.write(nb,tmp);tmp.replace(notebook)
        reread=nbformat.read(notebook,as_version=4);assert reread.cells[1:1+len(baseline)]==baseline
        errors=sum(o.get('output_type')=='error' for c in addition.cells for o in c.get('outputs',[]));assert errors==0
        write_json(dest/'notebook_manifest.json',dict(status='completed',old_cells_preserved=len(baseline),
            added_cells=len(addition.cells)+1,errors=errors,before_sha256=before,after_sha256=sha(notebook),job=os.environ['SLURM_JOB_ID']))

if __name__=='__main__':run()
