"""Correct the existing notebook's baseline claims; preserve all other cells."""
from pathlib import Path
import json,os,hashlib,shutil
from datetime import datetime,timezone

ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/ptc_paper_baseline'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import pandas as pd
    path=ROOT/'notebooks/dgscrna_results.ipynb'
    backup=OUT/'notebook_before_baseline_correction.ipynb'
    if not backup.exists():shutil.copy2(path,backup)
    before=json.loads(path.read_text());after=json.loads(path.read_text())
    original=json.loads(backup.read_text())
    assert len(before['cells'])==153
    assert [before['cells'][i]['id'] for i in [0,126,130]]==['93716c97','85888a74','a71eb454']
    link='../results/hvg_ptc_20260916_v1/ptc_paper_baseline/STATUS.md'
    intro=('# DG-scRNA: complete GBM results and PTC baseline reconciliation\n\n'
      'This is the existing full results notebook. The GBM chapters and their complete clustering figure atlases are preserved. '
      'The existing two-group PTC experiments remain available below, but they are not yet a validated reproduction of the paper’s selected PTC workflow. '
      'Saved-label parity established table provenance only. The PTC chapter now records the original group-specific routes, recovered paper metrics, '
      'and the remaining rerun/Accuracy discrepancies. DG-scRNA annotation results mean terminal DL/refinement. '
      f'[Current paper-baseline evidence]({link}).')
    metrics=pd.read_csv(OUT/'paper_Table2_historical_source_comparison.csv')
    dg=metrics[(metrics.method=='DG') & (metrics.metric!='Accuracy')]
    assert dg.agrees_at_reported_4_decimals.all() and len(dg)==8
    md=('# PTC: original paper baseline first; two-group experiments retained\n\n'
      '**Correction following user review (2026-09-16): the original-workflow gate was not passed by replaying saved S2 labels. '
      'Do not use the later AllTissues/mean reconstruction to judge optimality of the paper’s selected workflow.**\n\n'
      '| Selected samples | Original marker | Original partition | Density cutoff | Endpoint |\n'
      '|---|---|---|---|---|\n'
      '| N-1/N-2/MT-1/MT-2 | CellMarker Thyroid | Seurat SNN | none | terminal DL |\n'
      '| T-1/T-2/TU-1/TU-2 | NCOMMREFF / Pubmed_34663816 | UMAP-HDBSCAN | mean | terminal DL |\n\n'
      'These choices are explicit in the original notebook’s `origin` and `annotation_map_orig`. '
      'The archived checkpoint integrates all eight samples together; the two groups select different annotation outputs. '
      'Rebuilding CCA independently in NMT and TTU changes that baseline. In particular, the rebuilt NMT feature set omits CD3D, '
      'the sole Thyroid T-cell marker, while the original checkpoint retains it.\n\n'
      '**Recovered Table 2 calculation:** original CSV validation flags (35,727 filtered-contig-positive cells), '
      'the original broad T/NK mapping, and the unmodified R `MLmetrics::F1_Score` default reproduce every DG F1/AUC entry at four decimals. '
      'The default treats class 0 as positive: these published F1 values are non-T-class F1, not weighted multiclass F1.\n\n'
      '| Tissue | Paper F1 | Recovered F1 | Paper AUC | Recovered AUC |\n|---|---:|---:|---:|---:|\n')
    for scope in ['Normal adjacent','Primary tumor','Metastasis','Overall']:
        f=dg[(dg.scope==scope)&(dg.metric=='F1 score')].iloc[0]
        a=dg[(dg.scope==scope)&(dg.metric=='AUC-ROC')].iloc[0]
        md+=f'| {scope} | {f.paper:.4f} | {f.reconstructed:.6f} | {a.paper:.4f} | {a.reconstructed:.6f} |\n'
    md+=('\nThe same original predictions give overall T-positive F1 0.9244205 and ordinary binary accuracy 0.941182. '
      'The paper’s Accuracy 0.8922 remains unexplained; it is not silently replaced. The S3 `T_cell` flag is a different vector '
      'and coincides exactly with the broad final T/NK calls, so it is not used as independently verified TCR truth. '
      'T+TU native labels match the saved Pubmed/UMAP-HDBSCAN/mean terminal branch for all 44,149 cells.\n\n')
    full_evaluation=OUT/'evaluation_full_parallel/manifest.json'
    evaluation=full_evaluation if full_evaluation.exists() else OUT/'evaluation_marker_union/manifest.json'
    if evaluation.exists():
        e=json.loads(evaluation.read_text())
        md+=('**New terminal rerun from the original checkpoint:**\n\n'
             '| Route | Cells | Exact native terminal calls | Terminal mismatches | DL status |\n'
             '|---|---:|---:|---:|---|\n')
        for r in e['routes']:
            md+=f"| {r['group']} | {r['n_cells']} | {r['n_terminal_exact']} | {r['n_terminal_mismatch']} | {r['dl_status']} |\n"
        md+='\nThis rerun uses explicit seed 42 and the original scoring/DL rules; historical initialization is unavailable. '
        if full_evaluation.exists():
            assert all(r['initial_exact']==r['terminal_exact']==r['n_cells'] and r['max_absolute_statistic_delta']<=1e-12 for r in e['full_gene_marker_union_parity'])
            md+='Full-gene R verification passed: all projected DEG statistics agree to 1e-12, and initial and terminal labels are identical to the marker-union computation across all 92,404 cells in each route. '
        else:
            md+='Full-gene verification is tracked separately. '
    if (OUT/'paper_accuracy_consistency_manifest.json').exists():
        md+='\n\nThe Accuracy discrepancy is a reporting inconsistency under the recovered binary endpoint: even allowing four-decimal rounding, F1 0.9519 requires ordinary accuracy at least 0.908124, above the paper’s 0.8922. Its different source/denominator or reporting error is still unresolved. '
        if (OUT/'original_writing_note/manifest.json').exists():
            md+='The recovered archive manuscript `tcr/scripts/paper.md` contains the same eight DG F1/AUC values in Table2, but no Accuracy row or ordinary-accuracy formula; it therefore does not supply the missing V16 calculation. '
    literal_paths=[OUT/'literal_original_DL'/r/'manifest.json' for r in ['NMT_Thyroid_Seurat_none','TTU_Pubmed_UMAPHDBSCAN_mean']]
    if all(p.exists() for p in literal_paths):
        literal=[json.loads(p.read_text()) for p in literal_paths]
        assert all(r['n_terminal_equal_wrapper']==r['n_cells']==92404 and r['model_max_absolute_weight_delta']==r['pool_max_absolute_probability_delta']==0 for r in literal)
        md+='\n\nDirect execution of the unmodified archived `run_dgscrna` (SLURM7339744, same input/seed42, original eight-worker loaders) exactly reproduces every model weight, pool probability and terminal label from the reconstruction for both routes. The remaining historical differences are therefore not caused by the refinement wrapper on this fixed input and environment. '
    diagnostic_path=OUT/'recovery_diagnosis/initialization_replicates/SUMMARY.json'
    if diagnostic_path.exists():
        diagnostic=json.loads(diagnostic_path.read_text())
        assert diagnostic['route_fits_completed']==12 and diagnostic['no_seed_selected']
        ranges=diagnostic['ranges']
        md+=('\n\n**Initialization diagnosis (2026-09-17):** five prespecified model seeds (0–4), with fixed original-checkpoint inputs, '
             'reconstructed initial calls, split seed42 and training settings, were run for both selected routes. '
             'A separate seed42 parity control reproduces the existing weights, probabilities and labels exactly. '
             f"The five combined runs span source-defined F1 {ranges['F1_class0']['min']:.6f}–{ranges['F1_class0']['max']:.6f} "
             f"and AUC {ranges['AUC_binary']['min']:.6f}–{ranges['AUC_binary']['max']:.6f}. "
             'These are observed five-seed ranges, not confidence intervals or proof of the historical cause. '
             'No seed is selected as the replacement historical baseline. The old pre-DL training-label snapshot remains unavailable in the checked cohort files. '
             '[Diagnosis and restoration path](../results/hvg_ptc_20260916_v1/ptc_paper_baseline/RESTORATION_DIAGNOSIS_ZH.md). ')
    md+=(f'\n[Baseline status, source links and discrepancy ledger]({link}). '
         'Later PTC ablation interpretation stays paused until the baseline discrepancies are resolved. '
         'The preserved analyses below remain separate two-group experiments. No independent Pu cohort is run.')
    provenance=('## Archived-output identity and original-workflow reconstruction\n\n'
      'R and Python verified exact S2 initial/native/general and S3 native/detailed-label identity against saved outputs. '
      'That check did not rerun the original scoring or train the DL model. The original checkpoint and S3 Seurat partition also agree in all 92,404 cells. '
      'The original selected-route reconstruction and Table 2 reconciliation are now tracked above; historical model initialization/weights remain unavailable. '
      'All raw-gene counts and QC previously matched after 15 Seurat underscore-to-dash symbol substitutions. '
      'The raw QC plus expected 7.5% doublet count reproduces saved cell totals, but does not independently recover historical stochastic doublet identities.')
    home_path=OUT/'home_original_sup_manifest.json'
    if home_path.exists():
        home=json.loads(home_path.read_text())
        assert home['archive_identity']['status']=='byte_identical_SHA256'
        assert home['table_identity']['status']=='all_selected_original_final_columns_exact'
        provenance+=('\n\n**Original final-label files located in home/work:** the user confirmed that Sup contains the final used labels. '
          '`/users/PCON0080/yimin/work/_archives/tcr.tar.gz` is byte-identical (full SHA256) to the recovered scratch archive. '
          '`tcr/scripts/annotations.xlsx` preserves the S2 final native/general labels; '
          '`tcr/rawdata/data_with_validation.csv` preserves the S3 final native/detailed labels. '
          f"All {home['table_identity']['n_compared_columns']} checked column pairs agree for all 92,404 cells. "
          'Both Sup endpoints are retained verbatim in `original_final_sup_references.csv.gz`. '
          '[Original files and column map](../results/hvg_ptc_20260916_v1/ptc_paper_baseline/HOME_ORIGINAL_FILES_ZH.md).')
    for i,s in [(0,intro),(126,md),(130,provenance)]:after['cells'][i]['source']=s.splitlines(keepends=True)
    for i in range(153):
        if i not in [0,126,130]:assert before['cells'][i]==after['cells'][i]==original['cells'][i]
    tmp=path.with_suffix('.baseline_update.tmp')
    tmp.write_text(json.dumps(after,ensure_ascii=False,indent=1)+'\n');tmp.replace(path)
    report=dict(job=os.environ['SLURM_JOB_ID'],updated_utc=datetime.now(timezone.utc).isoformat(),
        notebook=str(path),backup_sha256=sha(backup),current_sha256=sha(path),
        n_cells=153,changed_markdown_cells=[0,126,130],unchanged_other_cells=150,
        all_GBM_scientific_cells_and_outputs_preserved=True,
        scope='Correct original-baseline claim in the same notebook; no new notebook version and no webpage')
    (OUT/'notebook_correction_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':run()
