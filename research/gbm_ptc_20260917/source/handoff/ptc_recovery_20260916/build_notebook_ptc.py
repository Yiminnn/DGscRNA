"""Execute new PTC chapters and preserve all completed GBM scientific cells."""
import json
import os
import shutil
import sys
from pathlib import Path
from ptc_common import ROOT,BASE,RECOVERY,require_slurm,sha,utc,write_json

def run():
    require_slurm()
    import nbformat as nbf
    from nbclient import NotebookClient
    assert (BASE/'summary/COMPLETE').exists() and (BASE/'figures/COMPLETE').exists()
    summary=json.loads((BASE/'summary/manifest.json').read_text())
    figs=json.loads((BASE/'figures/figure_manifest.json').read_text())
    assert summary['status']==figs['status']=='complete'
    assert json.loads((BASE/'report_manifest.json').read_text())['status']=='complete'
    out=BASE.parent/'notebooks';out.mkdir(exist_ok=True)
    original=out/'dgscrna_v6_results.ipynb';old=nbf.read(original,as_version=4)
    assert len(old.cells)==98
    assert not any(o.output_type=='error' for c in old.cells if c.cell_type=='code' for o in c.get('outputs',[]))
    cells=[]
    def md(s):cells.append(nbf.v4.new_markdown_cell(s))
    def code(s):cells.append(nbf.v4.new_code_cell(s))
    md('# PTC: original-source recovery and complete two-group reconstruction\n\n'
       'This completed section supersedes the PTC-pending statements in the dated GBM v6 and archived v5 narrative. '
       'Every headline DG-scRNA annotation is the terminal DL/refinement output. “RL” in the user request means this final stage, not reinforcement learning. '
       'Marker-only and confidence0.70 appear as explicitly named ablations. '
       'MTN=MT-1/MT-2/N-1/N-2 and TUT=TU-1/TU-2/T-1/T-2; no independent Pu expression cohort is included.')
    code('''from pathlib import Path
import os, json, hashlib
import pandas as pd
from IPython.display import display, Markdown, Image
assert os.environ.get("SLURM_JOB_ID"), "Scientific notebook execution requires SLURM. Saved results can be viewed anywhere."
PTC = Path("/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/ptc_experiments")
PTC_REC = PTC.parent / "ptc_recovery"
ptc_manifest = json.loads((PTC / "summary/manifest.json").read_text())
assert ptc_manifest["status"] == "complete"
assert ptc_manifest["n_terminal_conditions"] == 11696
assert ptc_manifest["n_sample_stage_rows"] == 49776
for name, digest in ptc_manifest["outputs"].items():
    assert hashlib.sha256((PTC / "summary" / name).read_bytes()).hexdigest() == digest, name
ptc_conditions = pd.read_csv(PTC / "summary/condition_verification.csv.gz")
ptc_metrics = pd.read_csv(PTC / "summary/all_sample_stage_metrics.csv.gz",low_memory=False)
ptc_hvg = pd.read_csv(PTC / "summary/HVG_paired_patient_effects.csv")
display(pd.Series({k:ptc_manifest[k] for k in ["n_cells","n_samples","n_patients","n_evaluation_units","n_terminal_conditions","n_sample_stage_rows","n_physical_partitions","n_unique_caches","n_unique_trained_models"]},name="Independently verified PTC scope"))
def ptc_figure(name):
    display(Image(filename=str(PTC / "figures" / (name + ".png")), width=1400))
''')
    md('## Complete workflow and tested decisions\n\n'
       'The left lane preserves the completed GBM design; the right lane shows original-source recovery, fixed groups, matched R controls, '
       'correction, the full original marker roster and the terminal classifier. '
       'Archived reference annotations and TCR outcomes never enter geometry or DL fitting. Marker-derived seed calls train the classifier. Panel selection is a later evaluation procedure.')
    code('ptc_figure("decision_tree_complete")')
    md('## Exact original-output replay and raw-count provenance\n\n'
       'R and Python reproduce all S2 initial/native/general annotations and all S3 native/detailed-T annotations cell by cell. '
       'The original trained weights and initialization seed are unavailable: this is exact archived-output replay, not a claim of bitwise retraining. '
       'All raw-gene counts and QC match after 15 Seurat underscore-to-dash symbol substitutions. '
       'The reference QC plus expected 7.5% doublet counts yields the exact saved cell totals, but does not independently recover the historical stochastic doublet identities.')
    code('''display(pd.read_json(PTC_REC / "R_baseline_replay/replay.json",typ="series"))
display(pd.read_csv(PTC / "summary/raw_counts_QC_parity.csv")[["sample","n_raw_cells","n_after_raw_reference_QC","raw_expected_doublets","n_archived_S2","n_retained_cells_raw_count_exact","n_retained_cells_raw_nCount_exact","n_retained_cells_raw_nFeature_exact","n_retained_cells_raw_mt_exact","gene3_first_QC_shortfall_vs_raw_QC"]])
ptc_figure("ptc_source_QC_audit")''')
    md('## HVG amount, PCA preprocessing and the complete method factorial\n\n'
       'The primary transfer contrast is HVG2000 minus all under direct UMAP2/HDBSCAN15/15, fixed AllTissues/mean, full RNA scoring and constant RNA2000 DL features. '
       'Patient-paired sample differences are averaged within patient. With four patients there are 16 exact sign flips and 256 exhaustive size-four bootstrap resamples. '
       'The smallest two-sided p-value is 0.125; bootstrap intervals are descriptive. '
       'S2/S3 concordance measures similarity to historical predictions, not independent accuracy; TCR-positive recall supplies complementary positive evidence.')
    code('''display(ptc_hvg[ptc_hvg.feature.eq("hvg2000")][["path","metric","n_paired_samples","n_paired_patients","delta_mean","ci_lower","ci_upper","p_exact","p_holm"]])
display(pd.read_csv(PTC / "summary/HVG_by_PCA_interaction.csv"))
display(pd.read_csv(PTC / "summary/HVG_by_group_patient_effects.csv").query("feature == 'hvg2000'"))
ptc_figure("ptc_HVG_PCA_response")
ptc_figure("ptc_complete_factorial")''')
    md('All six feature levels × seven representations × three clusterers use the same eight samples. '
       'The complete paired method comparisons are exploratory, with Holm correction over 20 competing defaults per feature/endpoint. '
       'Five seeds measure computational variation and are not additional patients. '
       'A reproduced GMM float32 covariance failure in N-2/HVG2000/t-SNE2 was recovered using float64 at unchanged K/regularization/seed; '
       'the affected condition and its exclusion sensitivity remain explicit.')
    code('''display(pd.read_csv(PTC / "verification/geometry_census/actual_features_and_geometry.csv").groupby("feature")[["actual_geometry_genes","actual_DL_genes"]].agg(["min","max","count"]))
ptc_methods = pd.read_csv(PTC / "summary/paired_method_comparisons.csv")
display(ptc_methods[ptc_methods.metric.eq("S2_macro_F1_reference_present")])
display(pd.read_csv(PTC / "summary/seed_variability.csv",header=[0,1],index_col=[0,1,2,3]))
display(pd.read_csv(PTC / "summary/numerical_recovery_sensitivity.csv"))''')
    code('''ptc_optimization = pd.read_csv(PTC / "verification/optimization_diagnostics/all_diagnostics.csv.gz")
display(ptc_optimization.groupby(["stage","method","feature"])[["convergence_warning","iteration_budget_reached"]].agg(["sum","count"]))''')
    md('## Matched single-sample, pooling and batch correction\n\n'
       'Each group has one sample from each of four patients. MTN and TUT are the requested new groups; the archived S2 checkpoint had all eight samples in one CCA. '
       'New matched R single-sample controls use exactly the pooled group-selected 2000 genes, full-RNA normalization/scoring, fixed DL inputs, '
       'UMAP30-neighbor/cosine/min.dist0.3 settings, SNN0.5 and HDBSCANminPts50. '
       'Pooling and correction can therefore be compared without changing implementation and density settings simultaneously. '
       'CCA acts on expression before PCA; Harmony acts on PCA30 before UMAP. Method and correction stage change together, so this is not a pure causal position test. '
       'Sample, patient and tissue are confounded; correction quality cannot be determined from mixing alone.')
    code('''ptc_batch = pd.read_csv(PTC / "summary/matched_batch_effects.csv")
display(ptc_batch[ptc_batch.metric.isin(["S2_macro_F1_reference_present","productive_strict_apparent_F1","coverage"])][["group","space","clusterer","contrast","metric","delta_mean","ci_lower","ci_upper","p_exact","p_holm"]])
ptc_figure("ptc_matched_pooling_batch")
display(pd.read_csv(PTC / "summary/legacy_integrated_effects.csv"))''')
    md('## Original marker libraries and cutoff conditions\n\n'
       'The actual recovered roster contains 17 symbol libraries, although the manuscript lists 16. '
       'Exact symbols, panel sizes and duplicate-marker denominators are retained. GSE184362/Pubmed34663816 is used as an original marker library only. '
       'The Ensembl-marker RDS is a different file and is not applied to symbol counts. '
       'The primary score retains all RNA genes regardless of geometry HVGs. The legacy CCA integrated-assay scorer and DL input are a separately labelled intervention.')
    code('''display(ptc_conditions.groupby(["mode","library","cutoff"]).size().rename("terminal conditions").unstack(fill_value=0))
ptc_retention = pd.read_csv(PTC / "evaluation_reference/geometry_marker_retention.csv.gz")
display(ptc_retention.groupby(["feature","library"])[["full_panel_denominator","scoring_retained","geometry_retained","scoring_retained_but_geometry_lost"]].mean())
display(pd.read_csv(PTC / "vocabulary_audit/original_library_vocabulary.csv"))
display(pd.read_csv(PTC / "vocabulary_audit/reference_vocabulary_capacity_by_sample.csv"))
ptc_figure("ptc_marker_roster_TCR")
ptc_figure("ptc_marker_roster_coverage")''')
    md('## Real terminal DL, no-op states and threshold sensitivity\n\n'
       'PTC preserves the archived MLP256/128, LeakyReLU, Softmax before CrossEntropy, Adamax1e-3, ten epochs and batch size 256. '
       'The known-cell split is90/10 with seed42; new model initialization is also fixed at42. '
       'Known calls are never overwritten. Four-decimal rounded confidence 0.90 is the primary terminal rule; 0.70 is evaluated from the same saved probabilities. '
       'Empty pools are no-op; no-known-label cases retain Undecided; single-known-class training is explicit. '
       'An independent model implementation verifies saved weights, all pool probabilities, exact terminal labels and both thresholds for every used model cache.')
    code('''display(pd.read_csv(PTC / "summary/terminal_execution_census.csv"))
display(pd.read_csv(PTC / "summary/unique_model_census.csv"))
ptc_figure("ptc_refinement_threshold")
ptc_refinement = pd.read_csv(PTC / "summary/refinement_and_threshold_effects.csv.gz")
display(ptc_refinement[ptc_refinement["mode"].eq("pooled") & ptc_refinement.library.eq("CellMarker_AllTissues") & ptc_refinement.cutoff.eq("mean") & ptc_refinement.space.eq("UMAP2")])''')
    md('## Orthogonal evidence, Unknown cells and preserved biology\n\n'
       'Productive high-confidence TCR is positive evidence, not a complete binary gold standard. Strict T excludes NK, NKT and ambiguous lymphoid names; '
       'the historical broad mapping is reported as a separate sensitivity. TCR-negative predicted T cells are accompanied by RNA-core support, not automatically called false positives. '
       'Unknown coverage, RNA support, nCount, nFeature and mitochondrial fraction remain visible. '
       'Balanced shared-lineage mixing and within-sample neighbor retention are paired with RNA state diagnostics; RNA modules can overlap annotation markers and are descriptive.')
    code('''ptc_main = ptc_metrics[ptc_metrics.stage.eq("terminal_DL090") & ptc_metrics["mode"].eq("pooled") & ptc_metrics.scoring_assay.eq("RNA") & ptc_metrics.library.eq("CellMarker_AllTissues") & ptc_metrics.cutoff.eq("mean")]
display(ptc_main.groupby(["group","correction","space","clusterer"])[["productive_strict_recall","productive_strict_detection_yield","productive_strict_apparent_F1","productive_permissive_apparent_F1","coverage","undetected_predicted_T_RNA_T_support","Unknown_RNA_T_support","Unknown_median_nCount","Unknown_median_nFeature","Unknown_median_percent_mt"]].mean())
display(pd.read_csv(PTC / "batch_biology/shared_lineage_availability.csv"))
ptc_figure("ptc_batch_biology")''')
    md('The following prespecified marker panels and per-patient/tissue composition plots expand the scalar metrics. '
       'All cells are shown in the pooled embeddings; each correction has its own coordinate system. '
       'RNA evidence is stratified into called T, Unknown and other called cells, each with productive TCR detected/not detected. '
       'The historical doublet model/score was not recovered; the QC table contains exact source metadata and does not invent new doublet labels.')
    code('''ptc_figure("ptc_RNA_marker_evidence")
ptc_figure("ptc_patient_tissue_composition")
ptc_figure("ptc_embedding_TCR_lineage")
ptc_unknown_QC = pd.read_csv(PTC / "biology_detail/TCR_stratified_Unknown_QC_by_sample.csv")
display(ptc_unknown_QC[ptc_unknown_QC.stage.eq("terminal_DL090") & ptc_unknown_QC.stratum.str.startswith("Unknown")])''')
    md('## Patient-held-out context choice\n\n'
       'Within a fixed group/correction/space/clusterer, three patients select among 17×3 library/cutoff settings by mean apparent productive-TCR binary F1. '
       'The fourth patient supplies held-out endpoints. Ties are deterministic. '
       '**The predictions were fitted transductively**: expression features, embedding, clustering and DL were not refitted excluding the held-out patient. '
       'This validates the context-selection step only, not a fully inductive model or independent external cohort.')
    code('''display(pd.read_csv(PTC / "summary/patient_heldout_panel_selection.csv").query("policy == 'heldout_selected'"))
display(pd.read_csv(PTC / "summary/patient_heldout_panel_summary.csv"))
ptc_figure("ptc_heldout_context")''')
    md('## Archived competitors and unresolved source versions\n\n'
       'Archived competitors retain their original input/reference conditions. They are not ranked as a fair matched contest against the new two-group reconstruction. '
       'The 623 apparent SCINA label differences are all UTF8/Latin1 gamma-delta encoding artifacts; compact cell names are mapped through the frozen name ontology. '
       'S3 TCR still differs in 5,435 cells from the S2 any-contig definition, and its DG binary T flag differs in 3,901 cells from the original source CSV. '
       'Literal source flags and native-name-based calls are separate sensitivities. The mismatches are preserved, not repaired by choosing favorable labels.')
    code('''display(pd.read_csv(PTC / "archived_comparators/SCINA_623_encoding_parity.csv"))
ptc_archived = pd.read_csv(PTC / "archived_comparators/historical_TCR_endpoint_replay_corrected_names.csv")
display(ptc_archived[ptc_archived["sample"].eq("ALL") & ptc_archived.TCR_definition.eq("TCR_cell_high_confidence_productive_TCR")])
ptc_figure("ptc_archived_comparator_audit")''')
    md('## Numerical recovery and complete provenance\n\n'
       'The R dbscan1.2.5/1.2.6 crashes are distinguished from OOM. A private1.2.6.9001 build fixes a 32-bit lower-triangle index product for the 48,255-cell MTN group. '
       '100,004 boundary-index checks and small-data exactness tests pass; all 24 full TUT partitions match the previous builds exactly. '
       'The empty-validation-index and cross-node cache-publication fixes preserve the original model and outputs. '
       'All raw/source hashes, scripts, checkpoints, retries and SLURM jobs remain available; shared environments were not upgraded.')
    code('''display(pd.read_csv(PTC / "verification/recovery_parity/TUT_all_partitions_exact.csv"))
display(pd.read_json(PTC / "verification/geometry_census/manifest.json",typ="series"))
display(Markdown((PTC / "PI_BRIEF_ZH.md").read_text()))
print("Complete English PTC report:", PTC / "PTC_REPORT.md")
print("Complete source and analysis protocol:", PTC.parent.parent.parent / "handoff/ptc_recovery_20260916")''')
    kernel_root=BASE.parent/'jupyter';kernel=kernel_root/'kernels/dgscrna_ptc';kernel.mkdir(parents=True,exist_ok=True)
    write_json(kernel/'kernel.json',dict(argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'],display_name='DG-scRNA PTC pinned environment',language='python'))
    os.environ['JUPYTER_PATH']=str(kernel_root)+(os.pathsep+os.environ['JUPYTER_PATH'] if os.environ.get('JUPYTER_PATH') else '')
    section=nbf.v4.new_notebook(cells=cells,metadata={'kernelspec':{'name':'dgscrna_ptc','display_name':'DG-scRNA PTC pinned environment','language':'python'}})
    NotebookClient(section,timeout=1200,kernel_name='dgscrna_ptc',resources={'metadata':{'path':str(out)}}).execute()
    assert not any(o.output_type=='error' for c in section.cells if c.cell_type=='code' for o in c.get('outputs',[]))
    section_path=out/'ptc_complete_executed.ipynb';nbf.write(section,section_path)
    intro=nbf.v4.new_markdown_cell('# DG-scRNA v7: complete GBM and PTC evidence\n\n'
        'GBM was completed and independently verified first. The newly executed PTC section now follows the current GBM chapters and precedes the dated v5 archive. '
        'The 98 preserved v6 cells retain their original outputs and wording; earlier statements that PTC is pending are superseded by the completed PTC section. '
        'All DG-scRNA headline results use terminal DL/refinement. See the complete workflow, source/QC reconciliation and the Chinese PI briefing in the PTC chapter.')
    position=next(i for i,c in enumerate(old.cells) if c.cell_type=='markdown' and '# Archived v5 results' in c.source)
    original_cells=list(old.cells)
    old.cells=[intro]+old.cells[:position]+section.cells+old.cells[position:]
    old.metadata['ptc_completion']=dict(source_GBM_notebook_sha256=sha(original),preserved_GBM_archive_cells=98,
        newly_executed_PTC_cells=len(section.cells),summary_sha256=sha(BASE/'summary/manifest.json'),completed_at=utc())
    target=out/'dgscrna_v7_results.ipynb';nbf.write(old,target)
    preserved=old.cells[1:position+1]+old.cells[position+1+len(section.cells):]
    assert preserved==original_cells
    central=ROOT/'notebooks/dgscrna_results.ipynb';archive=central.parent/'archive';archive.mkdir(exist_ok=True)
    if central.exists() and sha(central)!=sha(target):
        backup=archive/f'dgscrna_results_{sha(central)[:12]}.ipynb'
        if not backup.exists():shutil.copy2(central,backup)
    shutil.copy2(target,central)
    write_json(out/'ptc_execution_manifest.json',dict(status='complete',notebook=str(target),notebook_sha256=sha(target),
        workspace_notebook=str(central),workspace_notebook_sha256=sha(central),n_cells=len(old.cells),
        newly_executed_PTC_cells=len(section.cells),preserved_GBM_archive_cells=len(preserved),
        standalone_PTC_notebook_sha256=sha(section_path),source_GBM_notebook_sha256=sha(original),
        source_sha256=sha(Path(__file__)),summary_sha256=sha(BASE/'summary/manifest.json'),
        figure_manifest_sha256=sha(BASE/'figures/figure_manifest.json'),completed_at=utc(),job=os.environ['SLURM_JOB_ID']))
    print('Executed PTC section; preserved 98 original GBM/archive cells; published',target,flush=True)

if __name__=='__main__':run()
