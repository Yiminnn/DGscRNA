#!/usr/bin/env python3
"""Execute the expanded GBM notebook, preserving the sealed v5 package as an archive."""
import json
import os
from pathlib import Path
import shutil
import sys
from common import ROOT,OUT,require_slurm,sha,utc,write_json


def build():
    require_slurm()
    import nbformat as nbf
    from nbclient import NotebookClient
    import pandas as pd
    summary=json.loads((OUT/'summary/manifest.json').read_text())
    verify=json.loads((OUT/'verification/independent_manifest.json').read_text())
    figures=json.loads((OUT/'figures/figure_manifest.json').read_text())
    assert summary['status']=='complete' and verify['status']=='complete' and figures['status']=='complete'
    assert (OUT/'singleton_robustness_summary/COMPLETE').exists()
    assert (OUT/'method_comparisons/COMPLETE').exists()
    assert summary['n_conditions']==39930 and summary['n_samples']==121
    old=ROOT/'results/g274_v5_delivery_v1/package_v3'
    archive=OUT/'archive_v5'
    if not archive.exists():shutil.copytree(old,archive)
    dest=OUT/'notebooks';dest.mkdir(parents=True,exist_ok=True)
    cells=[]
    def md(s):cells.append(nbf.v4.new_markdown_cell(s))
    def code(s):cells.append(nbf.v4.new_code_cell(s))
    md('# DG-scRNA v6: GBM HVG ablations through terminal annotation\n\n'
       'This update reports the complete predeclared GBM experiment. All DG-scRNA annotation endpoints use the terminal DL/refinement output. '
       'Clustering concordance is a separate endpoint. Marker-only outputs appear only in explicitly named ablations. '
       'The original v5 notebook follows as a dated archive; its earlier PTC table excerpt is not the requested R reproduction. '
       'PTC work starts after this GBM delivery and will be appended once the original-table reproduction gate is resolved.')
    code('''from pathlib import Path
import os, json, hashlib
import numpy as np
import pandas as pd
from IPython.display import display, Markdown, Image
assert os.environ.get("SLURM_JOB_ID"), "Notebook execution requires SLURM; saved outputs can be viewed anywhere."
candidates = [Path.cwd(), Path.cwd().parent, Path("/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1")]
NEW_RUN = next(p.resolve() for p in candidates if (p / "summary/manifest.json").exists())
sm = json.loads((NEW_RUN / "summary/manifest.json").read_text())
vm = json.loads((NEW_RUN / "verification/independent_manifest.json").read_text())
assert sm["status"] == vm["status"] == "complete"
for name, digest in sm["outputs"].items():
    assert hashlib.sha256((NEW_RUN / "summary" / name).read_bytes()).hexdigest() == digest, name
metrics = pd.read_csv(NEW_RUN / "summary/all_sample_conditions.csv.gz")
condition_summary = pd.read_csv(NEW_RUN / "summary/condition_summary.csv")
paired_effects = pd.read_csv(NEW_RUN / "summary/paired_patient_effects.csv")
assert len(metrics) == 39930 and metrics["sample"].nunique() == 121
assert metrics.loc[metrics.evaluable.eq(True), "sample"].nunique() == 97
assert metrics.loc[metrics.evaluable.eq(True), "patient"].nunique() == 55
display(pd.Series({"samples":121, "primary samples":97, "primary patients":55,
                   "planned annotation conditions":39930,
                   "saved terminal outputs":vm["n_terminal_conditions_checked"],
                   "actual DL histories independently checked":vm["n_actual_DL_checked"]}, name="Verified scope"))
display(metrics.groupby("status").size().rename("conditions").to_frame())
def new_figure(name, caption=None):
    if caption: display(Markdown(caption))
    display(Image(filename=str(NEW_RUN / "figures" / (name + ".png")), width=1300))
''')
    md('## Main workflow and ablation map\n\n'
       'A1 changes geometry features; A2 changes PCA preprocessing; A3 changes representation and dimensionality; '
       'A4 changes clustering and its separately named sensitivities; A5 tests marker retention and scoring truncation; '
       'A6 compares seed calls with the terminal refinement output. Reference labels never enter fitting. '
       'K23 is predeclared from the fixed vocabulary size and is a conditional comparison; label-free K selection is separate.')
    code('new_figure("decision_tree_main")')
    md('## Primary hypothesis: HVG2000 versus all genes\n\n'
       'The fixed primary path is direct UMAP2 → HDBSCAN15/15. Paired sample differences are averaged within patient; patients are the independent units. '
       'The two primary endpoints are partition ARI and terminal strict-L1 macro-F1 over present reference classes. '
       'The legacy set-valued Lfine score is retained as a separate historical concordance measure. '
       'Intervals use 10,000 patient-bootstrap replicates. Both the two-endpoint primary Holm correction and the five-HVG-level secondary corrections are saved. '
       'This is a controlled reanalysis of previously inspected data, not a new independent validation cohort.')
    code('''primary_test = paired_effects[paired_effects.primary_contrast.eq(True)]
display(primary_test[["metric","n_paired_samples","n_paired_patients","delta_mean","ci_lower","ci_upper","p_value","p_holm_2_primary_endpoints"]])
display(paired_effects[(paired_effects.cohort == "primary97") & (paired_effects.path == "direct_UMAP2")][
    ["feature","metric","n_paired_samples","n_paired_patients","delta_mean","ci_lower","ci_upper","p_holm_5_HVG_levels"]])''')
    caption_by_name={r['name']:r['caption'] for r in figures['figures']}
    for name,title in [('hvg_response_curves','HVG amount and terminal annotation'),
                       ('primary_factorial_heatmaps','Complete 6 × 7 × 3 primary comparison'),
                       ('paired_method_comparisons','Does UMAP/HDBSCAN outperform the other defaults?'),
                       ('feature_method_interaction','Does the HVG effect depend on the method?'),
                       ('umap_dimension_factorial','PCA preprocessing and UMAP dimensions'),
                       ('marker_retention_and_scoring_ablation','Geometry selection versus scoring-marker loss'),
                       ('common_reference_geometry_diagnostics','Geometry preservation in a common reference')]:
        if name in caption_by_name:
            md('## '+title+'\n\n'+caption_by_name[name]);code(f'new_figure({name!r})')
        if name=='primary_factorial_heatmaps':
            code('''primary_table = condition_summary[(condition_summary.cohort == "primary97") & condition_summary.families.str.contains("E1_primary_table")].copy()
primary_table["method"] = primary_table.dr + " / " + primary_table.clusterer
display(primary_table.pivot(index="feature", columns="method", values="terminal_strict_L1_macroF1_present_n_samples"))
display(primary_table[["feature","method","n_expected_samples","n_terminal_outputs","n_structural_unavailable",
    "terminal_strict_L1_macroF1_present_patient_mean","terminal_strict_L1_macroF1_present_fixed_cohort_lower",
    "terminal_strict_L1_macroF1_present_fixed_cohort_upper"]])''')
        if name=='umap_dimension_factorial':
            md('The HVG effect is conditional on PCA preprocessing. The following exploratory four-condition interaction compares [HVG2000 − all] with direct genes against the same feature contrast after PCA30, at each predeclared UMAP dimension. Positive values indicate a larger HVG response in the direct-input workflow; Holm correction spans the three dimensions per endpoint.')
            code('display(pd.read_csv(NEW_RUN / "method_comparisons/hvg_pca_interactions.csv"))')
    md('## HVG definition, marker union and the historical geometry bridge\n\n'
       'The primary HVG sweep uses count-based VST. The legacy dispersion-based selections at 2,000 and 5,000 genes are separate controls. '
       'The marker-union condition adds detected panel markers to HVG2000 for geometry while keeping scoring and refinement unchanged. '
       'SCANPY_UMAP is a separate legacy workflow with PCA30 and its recorded neighbor/UMAP defaults; its result cannot be attributed solely to the HVG flavor.')
    code('''features = pd.read_csv(NEW_RUN / "summary/actual_feature_inventory.csv")
display(features.groupby("feature")[["actual_geometry_genes","final_DL_genes"]].agg(["min","median","max","count"]))
bridge = condition_summary[(condition_summary.cohort == "primary97") &
    condition_summary.feature.isin(["all","hvg2000","hvg5000","hvg2000_markers","seurat2000","seurat5000"]) &
    (condition_summary.clusterer == "HDBSCAN") & (condition_summary.seed == 42) &
    (condition_summary.min_cluster_size == 15) & (condition_summary.min_samples == 15) &
    (condition_summary.scoring_features == "all") &
    (((condition_summary.dr == "PCA") & (condition_summary.dim == 30)) |
     ((condition_summary.dr == "UMAP") & (condition_summary.dim == 2) & (condition_summary.neighbors == 15) & (condition_summary.min_dist == .1)) |
     (condition_summary.dr == "SCANPY_UMAP"))]
display(bridge[["feature","dr","dim","input_space","min_dist",
    "partition_ari_patient_mean","terminal_strict_L1_macroF1_present_patient_mean",
    "terminal_strict_L1_macroF1_present_n_samples"]].sort_values(["dr","input_space","feature"]))
display(pd.read_csv(NEW_RUN / "method_comparisons/marker_union_and_scoring_truncation_effects.csv"))''')
    md('## DL execution and structural limits\n\n'
       'This GBM experiment keeps the existing corrected refinement fixed at 15 epochs, batch size 256 and probability threshold 0.90, using all filtered scaled genes. '
       'Its settings are not silently substituted for the historical PTC R/DL implementation. '
       'An empty pool is a valid terminal no-op. All-noise or insufficient-contrast outputs are recorded as structural abstention. '
       'Insufficient training classes retain known calls and leave the pool Unknown, without claiming successful DL. '
       'The frozen Wilcoxon scorer cannot score singleton clusters: these conditions retain their partition metrics but have no fabricated final annotation score. '
       'TKU3074 has a default-span VST LOESS singularity and is outside the predefined primary cohort; its all-gene and legacy-flavor controls remain available. '
       'Missing scores are accompanied by fixed-denominator bounds. Comparing only available conditions is explicitly conditional.')
    code('''display(metrics.groupby(["feature","dl_status"], dropna=False).size().unstack(fill_value=0))
display(metrics.loc[metrics.status.str.startswith("structural"), ["sample","feature","dr","clusterer","structural_reason"]].groupby(["feature","dr","clusterer","structural_reason"], dropna=False).size().rename("n"))
dl_executed = metrics[metrics.training_executed.eq(True)]
display(dl_executed[["n_training","n_pool","n_dl_assigned"]].describe())
display(pd.read_csv(NEW_RUN / "verification/sample_gates.csv"))''')
    md('### Paired method evidence and refinement effects\n\n'
       'The following method and feature-by-method contrasts are exploratory analyses of the frozen factorial. '
       'The four-condition interaction distinguishes an HVG effect specific to UMAP/HDBSCAN from an effect shared by other methods. '
       'Shared-cohort rankings require availability for all 21 defaults and can be based on few samples; paired contrasts and fixed-cohort bounds accompany them. '
       'Refinement effects compare terminal outputs with explicitly named initial-call ablations; no-op outputs are counted without claiming DL was trained.')
    code('''rankings = pd.read_csv(NEW_RUN / "method_comparisons/method_rankings_with_shared_cohort.csv")
display(rankings[(rankings.feature.isin(["all","hvg2000"])) & (rankings.metric == "terminal_strict_L1_macroF1_present")].sort_values(["policy","feature","shared_patient_mean"], ascending=[True,True,False]))
method_tests = pd.read_csv(NEW_RUN / "method_comparisons/paired_method_contrasts.csv")
display(method_tests[(method_tests.policy == "primary_rule") & (method_tests.feature.isin(["all","hvg2000"]))])
display(pd.read_csv(NEW_RUN / "method_comparisons/feature_method_interactions.csv").query("feature == 'hvg2000'"))
display(pd.read_csv(NEW_RUN / "method_comparisons/terminal_refinement_effect.csv"))''')
    md('## Seed variation, label-free K and patient-held-out HVG choice\n\n'
       'The five seeds quantify algorithm variation; they are not five biological replicates. K is selected by silhouette, without reference labels, in a separately named sensitivity. '
       'For the held-out-patient HVG sensitivity, other patients choose the feature count, with a minimum 90% availability in training patients. '
       'These analyses do not mix tuned configurations into the default-method primary heatmap.')
    code('''display(pd.read_csv(NEW_RUN / "summary/seed_variability.csv", header=[0,1], index_col=list(range(7))).head(20))
selected_k = pd.read_csv(NEW_RUN / "summary/label_free_selected_K.csv")
display(selected_k.groupby(["feature","clusterer","k"]).size().rename("samples").to_frame())
heldout = pd.read_csv(NEW_RUN / "summary/leave_one_patient_out_HVG_selection.csv")
display(heldout)
display(heldout.groupby("metric")[["heldout_score","heldout_all","heldout_hvg2000"]].agg(["mean","count"]))''')
    md('## Secondary implementation sensitivity: singleton abstention\n\n'
       'This sensitivity was added after encountering unsupported singleton DEG contrasts. It does not replace the primary scoring rule. '
       'Every verified singleton-scoring failure is treated uniformly: preserve all cells and clusters, score supported groups, seed singleton groups Undecided, '
       'then use the same DL/refinement. No K, embedding, feature, marker or classifier parameter is changed. '
       'An infeasible fixed 10% stratified DL validation split remains unavailable. Both mean concordance and availability are shown.')
    code('''new_figure("secondary_singleton_robustness")
display(pd.read_csv(NEW_RUN / "singleton_robustness_summary/verification_counts.csv"))
display(pd.read_csv(NEW_RUN / "singleton_robustness_summary/secondary_primary_factorial.csv"))''')
    md('## Reproducibility and relation to the original comparison\n\n'
       'Each fit preserves cells, actual geometry/scoring/DL widths, code hashes, environment versions, seeds and SLURM provenance. '
       'The independent verifier reconstructs strict annotation metrics from exact confusion counts and checks real DL history widths and epochs. '
       'The all-gene branch was rerun in the current controlled pipeline; old and new partition parity is reported, including any nonlinear numerical sensitivity. '
       'The original DGCyTOF Table 4/5 used 13/32 protein channels, not an scRNA HVG-count sweep '
       '([original paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008885)). '
       'The historical scRNA HVG sweep is retained below as prior work, separate from this paired GSE274546 experiment.')
    code('''parity = pd.read_csv(NEW_RUN / "verification/legacy_allgene_partition_parity.csv")
display(parity.groupby(["dr","clusterer"]).agg(n=("sample","count"), exact=("exact_labels","sum"), minimum_partition_ari=("partition_ari","min"), mean_partition_ari=("partition_ari","mean")))
optimization = pd.read_csv(NEW_RUN / "verification/optimization_diagnostics.csv.gz")
display(optimization.groupby(["stage","method","feature"])[["convergence_warning","iteration_budget_reached"]].agg(["sum","count"]))
display(pd.read_json(NEW_RUN / "protocol/protocol.json", typ="series"))
print("Scientific scripts:", NEW_RUN.parent.parent / "handoff/hvg_ptc_20260916")
print("Full cell-level outputs remain under", NEW_RUN / "fits")''')
    md('---\n\n# Archived v5 results\n\nThe following cells reproduce the earlier sealed report package. They keep their original scopes and limitations; the new complete GBM comparison above supersedes cross-table inferences about HVG. The archived PTC section remains a prior table audit until the requested original R reproduction is added.')
    oldbook=nbf.read(old/'notebooks/v5_results.ipynb',as_version=4)
    for c in oldbook.cells:
        if c.cell_type=='code' and 'candidates = [Path.cwd(), Path.cwd().parent]' in c.source:
            c.source=c.source.replace('candidates = [Path.cwd(), Path.cwd().parent]',
                                      'candidates = [NEW_RUN / "archive_v5"]')
        cells.append(c)
    kernel_root=OUT/'jupyter'
    kernel_path=kernel_root/'kernels/dgscrna_hvg'
    kernel_path.mkdir(parents=True,exist_ok=True)
    write_json(kernel_path/'kernel.json',dict(argv=[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'],
                                            display_name='DG-scRNA pinned environment',language='python'))
    os.environ['JUPYTER_PATH']=str(kernel_root)+(os.pathsep+os.environ['JUPYTER_PATH'] if os.environ.get('JUPYTER_PATH') else '')
    book=nbf.v4.new_notebook(cells=cells,metadata={'kernelspec':{'name':'dgscrna_hvg','display_name':'DG-scRNA pinned environment','language':'python'}})
    target=dest/'dgscrna_v6_results.ipynb'
    nbf.write(book,target)
    client=NotebookClient(book,timeout=1200,kernel_name='dgscrna_hvg',resources={'metadata':{'path':str(dest)}})
    client.execute()
    nbf.write(book,target)
    assert not any(o.output_type=='error' for c in book.cells if c.cell_type=='code' for o in c.get('outputs',[]))
    # A stable workspace entry point is a copy, preserving all older notebooks and sealed packages.
    central=ROOT/'notebooks';central.mkdir(exist_ok=True)
    central_target=central/'dgscrna_results.ipynb'
    if central_target.exists() and sha(central_target)!=sha(target):
        prior=central/'archive';prior.mkdir(exist_ok=True)
        backup=prior/f'dgscrna_results_{sha(central_target)[:12]}.ipynb'
        if not backup.exists():shutil.copy2(central_target,backup)
    shutil.copy2(target,central_target)
    write_json(dest/'execution_manifest.json',dict(timestamp=utc(),status='completed',
        notebook_sha256=sha(target),n_cells=len(book.cells),summary_sha256=sha(OUT/'summary/manifest.json'),
        verification_sha256=sha(OUT/'verification/independent_manifest.json'),
        source_script_sha256=sha(Path(__file__)),workspace_notebook=str(central/'dgscrna_results.ipynb')))
    print(f'Executed {target} ({len(book.cells)} cells)',flush=True)


if __name__=='__main__':build()
