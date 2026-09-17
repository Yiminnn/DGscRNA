"""Publish passed computational checks while retaining the unresolved paper gate."""
from pathlib import Path
import os,json,hashlib
from datetime import datetime,timezone
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/ptc_paper_baseline'

def run():
    assert os.environ.get('SLURM_JOB_ID')
    m=json.loads((OUT/'evaluation_full_parallel/manifest.json').read_text())
    v=json.loads((OUT/'input_and_terminal_verification.json').read_text())
    assert v['status']=='passed' and v['loom_float32_differing_values']==0
    assert len(m['full_gene_marker_union_parity'])==2
    for p in m['full_gene_marker_union_parity']:
        assert p['initial_exact']==p['terminal_exact']==p['n_cells']==92404
        assert p['max_absolute_statistic_delta']<=1e-12
    home=json.loads((OUT/'home_original_sup_manifest.json').read_text())
    assert home['archive_identity']['status']=='byte_identical_SHA256'
    assert home['table_identity']['status']=='all_selected_original_final_columns_exact'
    literal=[]
    for route in ['NMT_Thyroid_Seurat_none','TTU_Pubmed_UMAPHDBSCAN_mean']:
        value=json.loads((OUT/'literal_original_DL'/route/'manifest.json').read_text())
        assert value['original_function_body_modified'] is False
        assert value['n_terminal_equal_wrapper']==value['n_cells']==92404
        assert value['model_max_absolute_weight_delta']==value['pool_max_absolute_probability_delta']==0
        literal.append(value)
    path=OUT/'STATUS.md';s=path.read_text()
    s=s.replace('## Active full-gene verification','## Completed full-gene verification')
    s=s.replace('The same full-gene statistical calls now use four forked workers inside each four-CPU/64-GB SLURM allocation. Both routes continue in parallel.',
                'The same full-gene statistical calls used four forked workers inside each four-CPU/64-GB SLURM allocation. Both routes completed.')
    s=s.replace('Baseline agreement and the Table 2 metric definition are still pending; no further ablation claim is released.',
                'The source-defined F1/AUC calculation is recovered. Exact historical terminal retraining and the source of the Accuracy row remain unreconciled; no further ablation claim is released.')
    s+='\n\nFull-gene array **7339596** and independent comparison **7339600** completed successfully. All marker-union DEG statistics match full-gene statistics to 1e-12; initial and terminal calls are identical in all 92,404 cells for both routes. Thus the computational gene-test optimization does not explain the historical DL mismatches. Full-gene models, predictions and lineage remain in `replay_selected_routes_full_parallel/`; verification is in `evaluation_full_parallel/marker_union_full_gene_parity.csv`. The paper-workflow gate remains **not fully passed**, because the recorded terminal differences and Accuracy-source issue persist.\n'
    s+='\n\n**Unmodified full DL function check (SLURM7339744):** both routes were also executed directly through the archived `run_dgscrna`, with the original eight-worker loaders and seed42. All 92,404 terminal labels, all trained weights and all pool probabilities exactly match the reconstruction in each route (maximum numeric delta0). Historical selected-group differences remain NMT16 and TTU1641. This excludes the refinement-wrapper implementation as their cause for this fixed input/seed/environment. The original function body was not edited; a constructor wrapper retained its model for saving. The initial attempts failed on a missing optional import and on modern AnnData converting the test fixture to categorical labels. Isolated dependencies and string-preserving fixture serialization resolved those compatibility issues; failed logs/fixtures remain preserved.\n'
    path.write_text(s)
    zh=OUT/'PAPER_RECONCILIATION_ZH.md';s=zh.read_text()
    s=s.replace('另有全基因 R 核验正在运行，见 `STATUS.md` 和 SLURM 7339596/7339600；它核实计算加速没有改变原始统计规则。',
        '全基因 R 核验已通过（SLURM 7339596/7339600）：两条路线在全部 92,404 个细胞的初始标签和最终 DL 标签上完全一致，marker 对应的五项 DEG 统计与全基因计算在 1e-12 内一致。因此，上述历史 DL 差异不是 marker-union 计算加速造成的。')
    s+='\n\n补充核验（SLURM7339744）：直接运行归档的完整 `run_dgscrna`，保留原始八个数据加载 worker、10轮训练和同一 seed42；两条路线的全部权重、概率、92,404 个最终标签与当前重跑实现完全一致。所以当前剩余历史差异不能归因于 DL 改写。原始函数未修改；测试环境补充缺失依赖，并保持输入标签为原始字符串类型。\n'
    zh.write_text(s)
    home_report=OUT/'HOME_ORIGINAL_FILES_ZH.md'
    hs=home_report.read_text()
    hs+='\n\nSLURM7339700 已完成：14 组列对应全部一致；home 与 scratch 两份完整归档均为 70,383,705,339 字节，SHA256 均为 `'+home['archive_identity']['files'][0]['sha256']+'`。因此当前恢复的原始结果确实来自用户指定的 home/work 归档。\n'
    home_report.write_text(hs)
    report=dict(status='all_requested_checkpoint_reconstruction_checks_finished; original_paper_gate_not_fully_passed',
        completed_utc=datetime.now(timezone.utc).isoformat(),job=os.environ['SLURM_JOB_ID'],
        full_gene_fit_array='7339596',full_gene_parity_job='7339600',
        original_DG_Table2_F1_AUC_entries_recovered=8,
        full_gene_marker_union_parity=True,exact_historical_terminal_reproduction=False,
        remaining_native_mismatches={r['group']:r['n_terminal_mismatch'] for r in m['routes']},
        further_ablation_started=False,original_Accuracy_source_unresolved=True,
        home_archive_SHA256_identical=True,original_Sup_final_files_located=True,
        original_Sup_column_pairs_exact=home['table_identity']['n_compared_columns'],
        unmodified_full_DL_function_parity=True,literal_DL_verifications=literal,
        missing_historical_provenance=['original model weights/initialization state','Table 2 Accuracy evaluation source/denominator'],
        evaluation_manifest_sha256=hashlib.sha256((OUT/'evaluation_full_parallel/manifest.json').read_bytes()).hexdigest())
    (OUT/'RECONSTRUCTION_STATUS.json').write_text(json.dumps(report,indent=2)+'\n')
    current=ROOT/'handoff/ptc_paper_baseline_20260916/CURRENT_TASK.md'
    current.write_text('''# PTC baseline reconciliation — verified scientific state

User constraints: original paper R workflow and terminal DL outputs first; no new PTC ablations or optimality claims until baseline reconciliation. Sup files are the original final used labels. Preserve both S2 and S3 final endpoints. No independent Pu cohort, new webpage or replacement notebook. All scientific computation remains through SLURM. Only yimin identity; OneDrive ownership and the existing delivery folder were explicitly confirmed. No subagents or Remote Control.

## Completed evidence

- Original notebook choices: NMT=Thyroid/Seurat/none/final DL; TTU=NCOMMREFF Pubmed34663816/UMAP-HDBSCAN/mean/final DL. Original CCA checkpoint integrates all eight samples; the two groups select different terminal outputs. New independent NMT integration had lost the only Thyroid T-cell marker CD3D and cannot replace this baseline.
- Home/work source7339700: `/users/PCON0080/yimin/work/_archives/tcr.tar.gz` is SHA256-identical to the scratch archive (`3c89e0b5e8bc67f3c61690c9f80de180c582cba65156dff4aa941a5266b21e4b`). S2 original file=`tcr/scripts/annotations.xlsx`; S3 native/detailed final file=`tcr/rawdata/data_with_validation.csv`. All14 checked column pairs agree for all92,404 cells. See `HOME_ORIGINAL_FILES_ZH.md` and `home_original_sup_manifest.json` under the baseline results.
- Historical metric recovery7339541: all8 DG Table2 F1/AUC entries match4decimals, using original any-filtered-contig validation (35,727 positives), original broad T/NK list, and original MLmetrics default F1 positive class0. Overall F1=.9518588, binary AUC=.939177; T-positiveF1=.9244205; ordinarybinaryaccuracy=.941182.
- Actual terminal fits7339544 and7339596 (marker-union and full2000gene tests) both ran. Full-gene verification7339600 passed: all projected DEG statistics agree to1e-12; all92,404 initial and terminal labels are identical for both routes. Full NMT job completed26:20, TTU40:29, noOOM.
- Input/saved-model verification7339582: all92,404x2,000 float32 values and ordering match the original loom exactly; saved model probabilities reproduce exactly.
- Literal archived `run_dgscrna`7339744: original function body/eight-worker loaders, same seed42 and fixed input reproduce all weights/probabilities/terminal labels of the reconstruction exactly. Earlier missing-import7339707 and categorical-fixture7339723 failures were resolved and preserved; isolated dependency install7339719 completed. No model/label tuning.
- Per-cell discrepancy audit7339591: NMT16 and TTU1641 native final differences remain against original S3 selections, all in the DL pool. TTU1315 historical-called->Unknown,165 reverse,161 called-type changes. New F1=.951692, AUC=.938690. Historical model state unavailable; cannot attribute its differences uniquely to random seed.
- Accuracy consistency7339608: the4 Table2 Accuracy entries cannot be ordinary binary accuracy of the same prediction vector as the recovered F1. DG F1=.9519 requires accuracy>=.908124 even allowing rounding, above manuscript.8922. No target/label changes permitted to force agreement.
- Original manuscript recovery7339775: `original_writing_note/paper.md` has the same8 DG F1/AUC values, but Table2 has no Accuracy row and its metric definitions omit ordinaryaccuracy. This draft cannot supply the V16 calculation.

## Remaining scientific limits

The original workflow gate is NOT fully passed: exact historical terminal retraining differs (NMT16/TTU1641), and the V16 Accuracy source/denominator remains unresolved. Existing source code itself saves final CSV, not model weights. Original final labels ARE recovered; do not ask the user for those labels again. Further PTC ablation/optimality interpretation stays paused.

## Delivery verification

Canonical notebook=`notebooks/dgscrna_results.ipynb`. Only markdown cells0,126,130 may change; all other150 cells/GBM figures are preserved against `notebook_before_baseline_correction.ipynb`. Check `notebook_correction_manifest.json` and actual notebook/backup.

Scientific report=`results/hvg_ptc_20260916_v1/ptc_paper_baseline/STATUS.md`; Chinese report=`PAPER_RECONCILIATION_ZH.md`; machine record=`RECONSTRUCTION_STATUS.json`. Final publication job7339754 runs notebook correction, staging, rclone copy and full-download check. Its authoritative receipt is `results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916/PTC_baseline_correction_upload_receipt_final.json`; inspect actual job exit and receipt before claiming delivery. It includes original final-label files and literal/full-gene results, excluding duplicated input matrices. Remote=`onedrive:work_od/share/dgscrna_GSE274546_TKU3186/01_report/GBM_PTC_results_20260916`.

The corrective /goal remains subject to its own completion audit. Verified computation and publication do not automatically mean exact original-paper reproduction. Preserve remaining uncertainty and do not declare a passed baseline solely from old-label parity or recovered F1/AUC.
''')
    claude=ROOT/'CLAUDE.md'
    c=claude.read_text()
    c=c.replace('The full-gene replacement is7339596 (NMT complete, TTU running), comparison7339600, final notebook/OneDrive update7339754.',
        'Full-gene7339596 and comparison7339600 completed with exact parity; final notebook/OneDrive update7339754 must be checked against its final receipt before claiming delivery.')
    claude.write_text(c)
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':run()
