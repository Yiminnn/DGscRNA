"""Validate all figure coverage and extend the existing full notebook on SLURM."""
from pathlib import Path
import copy
import hashlib
import json
import os
import shutil
import subprocess
import re
from datetime import datetime, timezone

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
BASE = ROOT / 'results/hvg_ptc_20260916_v1'
DEST = BASE / 'gbm_cluster_figures'
assert os.environ.get('SLURM_JOB_ID')


def sha(p):
    with Path(p).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def save(p, obj):
    p.parent.mkdir(exist_ok=True, parents=True)
    tmp = p.with_suffix(p.suffix + '.part')
    tmp.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    tmp.replace(p)


def cellhash(c):
    return hashlib.sha256(json.dumps(c, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


def main():
    import pandas as pd
    import nbformat as nbf
    from nbclient import NotebookClient
    samples = (BASE / 'protocol/samples.txt').read_text().split()
    gs = json.loads((BASE / 'protocol/geometries.json').read_text())
    conditions = {g['geometry_id'] + '/' + a['arm_id'] for g in gs for a in g['arms']}
    assert len(samples) == 121 and len(conditions) == 330
    tables, roster, figure_hashes = [], [], {}
    for sample in samples:
        folder = DEST / 'samples' / sample
        assert (folder / 'COMPLETE').read_text().strip() == sha(folder / 'manifest.json')
        m = json.loads((folder / 'manifest.json').read_text())
        assert m['status'] == 'complete' and m['n_conditions'] == 330
        assert m['plotted_all_cells'] and m['fit_and_annotation_sources_read_only']
        for name, expected in m['outputs'].items():
            assert sha(folder / name) == expected, (sample, name)
            figure_hashes[str((folder / name).relative_to(DEST))] = expected
        table = pd.read_csv(folder / 'condition_figure_index.csv')
        metrics = pd.read_csv(BASE / 'evaluation' / sample / 'metrics.csv').set_index('condition')
        assert set(table.condition) == conditions and not table.condition.duplicated().any()
        check = table.set_index('condition').loc[metrics.index]
        assert (check.cluster_plotted == metrics.partition_ari.notna()).all()
        assert (check.terminal_plotted == metrics.status.eq('completed')).all()
        assert (check.status == metrics.status).all()
        assert (table.n_cells == m['n_cells']).all()
        assert m['n_cluster_plots'] == int(table.cluster_plotted.sum())
        assert m['n_terminal_plots'] == int(table.terminal_plotted.sum())
        assert table.pdf_page.min() == 1 and table.pdf_page.max() == m['atlas_pages'] == 28
        for endpoint in ['clusters', 'terminal']:
            name = endpoint + '_all_conditions.pdf'
            info = subprocess.check_output(['pdfinfo', str(folder / name)], text=True)
            assert int(re.search(r'^Pages:\s+(\d+)', info, re.M).group(1)) == 28
            table[endpoint + '_pdf'] = 'samples/' + sample + '/' + name
        tables.append(table)
        roster.append(dict(sample=sample, patient=str(metrics.patient.iloc[0]),
            primary_cohort=bool(metrics.evaluable.iloc[0]), n_cells=m['n_cells'],
            conditions=330, cluster_plots=m['n_cluster_plots'], terminal_plots=m['n_terminal_plots'],
            clusters_pdf='samples/' + sample + '/clusters_all_conditions.pdf',
            terminal_pdf='samples/' + sample + '/terminal_all_conditions.pdf',
            pages_per_atlas=28, plot_job=m['slurm_job']))
    combined = pd.concat(tables, ignore_index=True)
    sample_table = pd.DataFrame(roster)
    assert len(combined) == 39930 and combined.terminal_plotted.sum() == 35211
    assert sample_table.primary_cohort.sum() == 97
    combined.to_csv(DEST / 'all_condition_figure_index.csv.gz', index=False)
    sample_table.to_csv(DEST / 'sample_atlas_index.csv', index=False)
    tk = DEST / 'samples/TKU3186'
    primary = sorted(tk.glob('primary_*.png'))
    secondary = sorted(tk.glob('secondary_clusters_page_*.png'))
    assert len(primary) == 12 and len(secondary) == 18
    summary = dict(status='complete', samples=121, conditions=39930,
        cluster_plots=int(combined.cluster_plotted.sum()), terminal_plots=35211,
        unavailable_cluster_panels=int((~combined.cluster_plotted).sum()),
        unavailable_terminal_panels=int((~combined.terminal_plotted).sum()),
        per_sample_atlases=242, total_atlas_pages=121*28*2,
        all_cells_plotted=True, refitting_performed=False,
        example='TKU3186: fixed previous notebook example, not selected by results',
        notebook_primary_conditions=126, notebook_additional_conditions=204,
        source_protocol_sha256=sha(BASE / 'protocol/geometries.json'),
        original_scientific_delivery_sha256=sha(BASE / 'EXPERIMENT_DELIVERY.json'),
        validation='Exact 121 x 330 condition census, source checksums, all-cell alignment, saved cluster/final vectors, independent availability parity, PDF page counts and full figure hashes.',
        slurm_job=os.environ['SLURM_JOB_ID'], completed_at=datetime.now(timezone.utc).isoformat(),
        outputs={name:sha(DEST / name) for name in ['all_condition_figure_index.csv.gz', 'sample_atlas_index.csv']})
    save(DEST / 'figure_hashes.json', figure_hashes)
    save(DEST / 'coverage_manifest.json', summary)

    baseline = DEST / 'baseline'
    baseline.mkdir(exist_ok=True)
    source_nb = BASE / 'notebooks/dgscrna_v7_results.ipynb'
    v6_path = BASE / 'notebooks/dgscrna_v6_results.ipynb'
    current = ROOT / 'notebooks/dgscrna_results.ipynb'
    expected_before = json.loads((BASE / 'notebooks/ptc_execution_manifest.json').read_text())['notebook_sha256']
    assert sha(source_nb) == expected_before
    backup = baseline / 'dgscrna_results_before_GBM_figures.ipynb'
    if not backup.exists():
        assert sha(current) == expected_before
        shutil.copy2(current, backup)
    assert sha(backup) == expected_before
    previous = nbf.read(backup, as_version=4)
    v6 = nbf.read(v6_path, as_version=4)
    assert len(previous.cells) == 126 and len(v6.cells) == 98
    assert previous.cells[1:36] == v6.cells[:35]
    assert previous.cells[63:] == v6.cells[35:]
    ptc = copy.deepcopy(previous.cells[36:63])
    assert len(ptc) == 27 and ptc[0].source.startswith('# PTC:')

    cells = []
    def md(text): cells.append(nbf.v4.new_markdown_cell(text))
    def code(text): cells.append(nbf.v4.new_code_cell(text))
    md('# GBM: every saved clustering condition and terminal annotation\n\n'
       'This section extends the full previous GBM notebook. Its original tables, figures and archived analyses remain below; '
       'the completed PTC chapter is preserved without changing any cell content or output. '
       'All **121 samples × 330 conditions** have a figure location, including explicitly unavailable panels. '
       'The figures read frozen embeddings, cluster IDs and final DL/refinement calls; no model was retrained and no condition was selected for visual appeal. '
       'TKU3186 is the previously fixed example. Its full 6 × 7 × 3 primary comparison and all remaining cluster conditions are displayed here; '
       'both complete atlases are provided for every sample.')
    code('''from pathlib import Path
import os, json
import pandas as pd
from IPython.display import display, Markdown, Image
assert os.environ.get("SLURM_JOB_ID"), "Execute through SLURM; saved notebook outputs can be viewed anywhere."
_candidates = [Path.cwd() / "gbm_cluster_figures", Path.cwd().parent / "results/hvg_ptc_20260916_v1/gbm_cluster_figures", Path("/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/gbm_cluster_figures")]
GBM_ATLAS = next(p for p in _candidates if (p / "coverage_manifest.json").exists())
atlas_manifest = json.loads((GBM_ATLAS / "coverage_manifest.json").read_text())
atlas_index = pd.read_csv(GBM_ATLAS / "all_condition_figure_index.csv.gz")
sample_atlases = pd.read_csv(GBM_ATLAS / "sample_atlas_index.csv")
assert len(atlas_index) == 39930 and len(sample_atlases) == 121
display(pd.Series({k:atlas_manifest[k] for k in ["samples", "conditions", "cluster_plots", "terminal_plots", "unavailable_cluster_panels", "unavailable_terminal_panels", "per_sample_atlases", "total_atlas_pages"]}, name="Complete figure coverage"))
ATLAS_LINK = "../results/hvg_ptc_20260916_v1/gbm_cluster_figures/"
def show_gbm_image(name):
    display(Image(filename=str(GBM_ATLAS / "samples/TKU3186" / (name + ".png")), width=1500))
''')
    md('## Reading the cluster and final-annotation figures\n\n'
       'Each primary matrix has seven representation rows and four columns: archived broad L1 reference, KMeans K23, GMM K23/diagonal and HDBSCAN15/15. '
       'The first matrix colors the saved cluster IDs; the second colors **terminal DL/refinement** outputs. Gray denotes noise or final Unknown/abstention. '
       'Numeric cluster labels are fit-specific and distinguish repeated colors. Every cell is plotted, with equal coordinate aspect. '
       'Partition ARI uses the archived fine reference; terminal macro-F1 uses strict L1. These reference annotations are not independent ground truth.\n\n'
       'For a native 2D fit, the saved fitting coordinates are shown directly. For a fit above 2D, only axes 1–2 are displayed, while the original clustering used every fitted dimension. '
       'No-DR clustering used the full selected gene space; its panels use one frozen all-gene PCA2 display for comparison. This display never changes the saved clusters. '
       'A valid cluster result can have an unavailable final annotation, for example when the frozen scoring rule cannot handle a singleton cluster. '
       'Those cases retain the cluster plot and a clearly marked unavailable annotation panel. No marker-only result substitutes for a missing terminal output.')
    for feature in ['all', 'hvg500', 'hvg1000', 'hvg2000', 'hvg3000', 'hvg5000']:
        title = 'All filtered genes' if feature == 'all' else 'HVG ' + feature[3:]
        md('## ' + title + ': every primary representation and clusterer\n\n'
           'The same cells and frozen marker library/cutoff are used across this matrix. '
           'HVG selection changes geometry features; primary scoring and DL retain all filtered genes. '
           'All rows use seed42. Compare the visual partitions with the final annotations and the patient-level comparisons above.')
        code(f'show_gbm_image("primary_{feature}_clusters")')
        code(f'show_gbm_image("primary_{feature}_terminal")')
    md('## All additional GBM clustering conditions for TKU3186\n\n'
       'The following atlas pages include every remaining condition: PCA preprocessing, UMAP dimensionality, seeds, neighbors/min_dist, '
       'HDBSCAN density settings, GMM covariance, label-free K candidates, legacy HVG definitions, marker union and scoring truncation. '
       'The first displayed page overlaps the final six primary conditions to preserve the original page numbering. '
       'All 330 terminal-annotation panels are also available in the linked final-output atlas. '
       'Separate scoring conditions may share identical cluster vectors; they remain indexed separately.\n\n'
       '[TKU3186 cluster atlas](../results/hvg_ptc_20260916_v1/gbm_cluster_figures/samples/TKU3186/clusters_all_conditions.pdf) · '
       '[TKU3186 terminal annotation atlas](../results/hvg_ptc_20260916_v1/gbm_cluster_figures/samples/TKU3186/terminal_all_conditions.pdf) · '
       '[Condition and page index](../results/hvg_ptc_20260916_v1/gbm_cluster_figures/samples/TKU3186/condition_figure_index.csv)')
    code('''for path in sorted((GBM_ATLAS / "samples/TKU3186").glob("secondary_clusters_page_*.png")):
    display(Image(filename=str(path), width=1500))
''')
    md('## Complete figure directory for all 121 GBM samples\n\n'
       'Each PDF contains 28 pages covering all 330 planned conditions. The cohort column preserves the fixed primary97 / secondary24 distinction. '
       'The two plotted-count columns count actual saved results; unavailable panels remain in their original positions. '
       'Use the full condition index to locate any feature, representation, clustering parameter, random seed or scoring condition by PDF page and panel.\n\n'
       '[All-condition index, CSV.gz](../results/hvg_ptc_20260916_v1/gbm_cluster_figures/all_condition_figure_index.csv.gz) · '
       '[Sample-level index, CSV](../results/hvg_ptc_20260916_v1/gbm_cluster_figures/sample_atlas_index.csv)')
    code('''lines = ["| Sample | Patient | Cohort | Cells | Cluster plots | Final plots | Figures |", "|---|---|---|---:|---:|---:|---|"]
for row in sample_atlases.itertuples():
    cohort = "primary97" if row.primary_cohort else "secondary24"
    links = f"[clusters]({ATLAS_LINK}{row.clusters_pdf}) / [terminal DL]({ATLAS_LINK}{row.terminal_pdf})"
    lines.append(f"| {row.sample} | {row.patient} | {cohort} | {row.n_cells} | {row.cluster_plots} | {row.terminal_plots} | {links} |")
display(Markdown("\\n".join(lines)))
''')
    md('## Figure provenance and preservation\n\n'
       'The figure ledger covers all 39,930 conditions and matches the original per-condition evaluation availability. '
       'Saved source checksums and cell alignment were checked before plotting. PDF page counts and image hashes were checked after export. '
       'Each sample folder includes the source hashes, figure manifest and condition-to-page mapping. '
       'Original model results and PTC results were opened read-only. The previously executed GBM/archive and PTC notebook cells are retained with their existing outputs.')
    code('display(pd.Series({k:atlas_manifest[k] for k in ["all_cells_plotted", "refitting_performed", "example", "validation", "source_protocol_sha256", "original_scientific_delivery_sha256"]}, name="Figure provenance"))')
    extension = nbf.v4.new_notebook(cells=cells, metadata=copy.deepcopy(v6.metadata))
    NotebookClient(extension, timeout=300, kernel_name='dgscrna_hvg',
                   resources={'metadata': {'path': str(BASE)}}).execute()
    assert not any(o.output_type == 'error' for c in extension.cells if c.cell_type == 'code' for o in c.get('outputs', []))
    nbf.write(extension, DEST / 'GBM_added_figure_cells.ipynb')

    intro = nbf.v4.new_markdown_cell('# DG-scRNA: complete GBM results with preserved PTC analysis\n\n'
        'This is the existing full results notebook with GBM clustering figures added. '
        'All 98 cells of the previous complete GBM/archive notebook are preserved with their outputs. '
        'The GBM extension shows every primary method/HVG combination for TKU3186, the remaining clustering conditions and links to complete atlases for all 121 samples. '
        'The 27-cell completed PTC chapter follows the GBM material unchanged. '
        'Earlier dated statements that PTC is pending are superseded by that preserved completed chapter. '
        'DG-scRNA annotation results always mean terminal DL/refinement output.')
    result = nbf.v4.new_notebook(cells=[intro] + copy.deepcopy(v6.cells[:35]) + extension.cells +
                               copy.deepcopy(v6.cells[35:]) + ptc,
                               metadata=copy.deepcopy(previous.metadata))
    result.metadata['gbm_figure_extension'] = dict(original_gbm_cells=98, preserved_ptc_cells=27,
        added_cells=len(extension.cells), source_gbm_sha256=sha(v6_path),
        previous_notebook_sha256=expected_before, coverage_manifest_sha256=sha(DEST / 'coverage_manifest.json'))
    nbf.validate(result)
    assert [cellhash(c) for c in result.cells[1:36] + result.cells[36 + len(extension.cells):-27]] == [cellhash(c) for c in v6.cells]
    assert [cellhash(c) for c in result.cells[-27:]] == [cellhash(c) for c in previous.cells[36:63]]
    assert not any(o.output_type == 'error' for c in result.cells if c.cell_type == 'code' for o in c.get('outputs', []))
    candidate = DEST / 'dgscrna_results.ipynb'
    nbf.write(result, candidate)
    reread = nbf.read(candidate, as_version=4)
    assert reread.cells == result.cells
    assert sha(current) in [expected_before, sha(candidate)], 'Canonical notebook changed externally; preserve it'
    temporary = current.with_suffix('.ipynb.part')
    shutil.copy2(candidate, temporary)
    temporary.replace(current)
    receipt = dict(status='NOTEBOOK_AND_FIGURES_VERIFIED', canonical_notebook=str(current),
        notebook_sha256=sha(current), previous_notebook_sha256=expected_before,
        preserved_gbm_archive_cells=98, preserved_ptc_cells=27, added_gbm_cells=len(extension.cells),
        notebook_cells=len(result.cells), preserved_ptc_cell_hashes=[cellhash(c) for c in ptc],
        source_gbm_cell_hashes=[cellhash(c) for c in v6.cells],
        coverage=summary, canonical_and_delivery_copy_equal=sha(current)==sha(candidate),
        job=os.environ['SLURM_JOB_ID'], completed_at=datetime.now(timezone.utc).isoformat())
    save(DEST / 'notebook_verification.json', receipt)
    (DEST / 'COMPLETE').write_text(sha(DEST / 'notebook_verification.json') + '\n')
    print(json.dumps({k:v for k,v in receipt.items() if k not in ['preserved_ptc_cell_hashes', 'source_gbm_cell_hashes']}), flush=True)


if __name__ == '__main__':
    main()
