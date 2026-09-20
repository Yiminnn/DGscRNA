"""Insert/refresh saved GBM evidence within the existing v5 notebook chapters.

All pre-campaign cells remain exactly unchanged. Only this script's tagged
chapter additions are replaced on later milestones, after saving the previous version.
"""
from pathlib import Path
from datetime import datetime, timezone
import copy
import csv
import hashlib
import json
import os
import re
import shutil
import subprocess

assert os.environ.get('SLURM_JOB_ID'), 'Notebook execution requires SLURM'
import nbformat as nbf
from nbclient import NotebookClient

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
CAMP = ROOT / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
OUT = CAMP / 'notebook'
TARGET = ROOT / 'notebooks/dgscrna_results.ipynb'
TAG = 'reviewer_completion_20260920'
OUT.mkdir(parents=True, exist_ok=True)


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def checked(path, manifest='manifest.json', flag='COMPLETE'):
    try:
        return (path / flag).read_text().strip() == sha(path / manifest)
    except OSError:
        return False


def link(path):
    return '../' + str(path.relative_to(ROOT))


def cellhash(cells):
    return hashlib.sha256(json.dumps(cells, sort_keys=True, ensure_ascii=False).encode()).hexdigest()


lock = OUT / 'update.lock'
lock.mkdir()
try:
    before = sha(TARGET)
    book = nbf.read(TARGET, as_version=4)
    original = [c for c in book.cells if TAG not in c.metadata.get('tags', [])]
    assert len(original) == 1067, 'Original v5 notebook has changed; audit before editing'
    if len(book.cells) == len(original):
        assert before == '9c38f579348b1052523cfdae839847389daf5e7b5f7b0eb32c12e66183214afb', 'Reaudit changed baseline'
    original_hash = cellhash(original)
    backup = OUT / f'before_{before[:16]}.ipynb'
    if not backup.exists():
        shutil.copy2(TARGET, backup)
    assert sha(backup) == before
    cells = []
    chapter = '9'

    def md(text):
        cells.append(nbf.v4.new_markdown_cell(text, metadata={'tags': [TAG], 'v5_option_a_chapter': chapter, 'v5_option_a_role': 'reviewer_addition'}))

    def code(text):
        cells.append(nbf.v4.new_code_cell(text, metadata={'tags': [TAG], 'v5_option_a_chapter': chapter, 'v5_option_a_role': 'reviewer_addition'}))

    at = datetime.now(timezone.utc).isoformat()
    md(f'''<a id="reviewer-completion-20260920"></a>
### GBM reviewer completion: controlled representations

Snapshot: {at}. New evidence is integrated into the existing v5 chapters; all 1,067 prior cells and outputs are preserved.

The original R/helper anchor uses HVG2000, PCA30→UMAP, R HDBSCAN, the original density scorer and terminal DL/refinement. The seven-representation comparison uses identical scaled-HVG input; direct-HVG UMAP is separate from the original anchor. HVG5000 is a predefined sensitivity analysis.

New A1 scope: 121 samples × 2 budgets × 7 representations × 13 clustering candidates; K is selected using training patients and terminal0.90 is the primary endpoint. A2 replaces cluster-DEG seeds with cell-wise marker seeds and retains the same DL. It is a seed-mechanism comparison. Cohort rankings await complete results and verification.

Current execution details: {link(ROOT / 'results/hvg_ptc_20260916_v1/paper_review_text_20260920/index.html')}
''')
    chapter = '1'
    code('''from pathlib import Path
import os, json
import pandas as pd
from IPython.display import display, Image, Markdown
assert os.environ.get('SLURM_JOB_ID'), 'Execute within SLURM'
REVIEW_ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
REVIEW_CAMP=REVIEW_ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
def review_image(path):
    display(Image(filename=str(path), width=1250))
''')
    chapter = '9'
    code('''review_status=[]
for folder,title in [('embedding','A1 embedding'),('no_clustering','A2 cell-wise seeds'),('controls','B neighbors / learning'),('comparison/scDeepSort_LogNormalize','C corrected scDeepSort')]:
    item=json.loads((REVIEW_CAMP/folder/'status.json').read_text())
    review_status.append({'Experiment':title,'Status':item.get('status'),'Processed':item.get('completed'),'Remaining':item.get('remaining'),'Recorded at':item.get('updated_at')})
display(pd.DataFrame(review_status))
''')

    evidence = CAMP / 'evidence'
    verification = json.loads((evidence / 'verification.json').read_text())
    assert verification['status'] == 'passed'
    manifest = json.loads((evidence / 'manifest.json').read_text())
    images = ['NL022_selected_marker_dotplot.png', 'NL022_selected_marker_violin.png',
              'NL022_markers_initial_terminal_same_coordinates.png']
    for name in images:
        assert sha(evidence / name) == manifest['files'][name]
    chapter = '4'
    md('''### Selected-marker RNA expression in NL022

NL022 is the predefined median-size display sample with 3,628 cells. RNA LogNormalize matches the original R export numerically. Complete tables retain all 1,364 available marker genes; 27 display genes were selected by overall detection frequency without using author-label contrasts.

Dot size shows the fraction of expressing cells; color shows mean log1p expression over every cell in the group. Violin plots and initial/terminal labels on shared coordinates follow. Because these marker genes also contribute to annotation, the figures describe expression consistency rather than independent biological validation or optimality.
''')
    for name in images:
        code(f"review_image(REVIEW_CAMP/'evidence/{name}')")
    md(f'''Complete marker tables and provenance:

{link(evidence / 'NL022_selected_library_marker_inventory.csv')}

{link(evidence / 'NL022_all_selected_marker_expression_by_endpoint.csv.gz')}

{link(evidence / 'NL022_selected_library_gene_source_records.csv.gz')}
''')

    chapter = '9'
    md('''### Controlled embedding examples and complete figure index

The following examples use the predefined small sample TKU4163/HVG2000/PCA2 with K-means K10 and HDBSCAN. They were not chosen by observed performance. Every candidate receives partition and terminal-label figures, including legitimate states in which DL was not executed. Original UMAP coordinates are shared for display; model fitting uses the representation stated in each title.
''')
    demo = CAMP / 'embedding/TKU4163/hvg2000/PCA2'
    proof = CAMP / 'no_clustering/a1_first_review/TKU4163/hvg2000/PCA2/validation.json'
    pdict = json.loads(proof.read_text())
    assert pdict['status'] == 'passed' and pdict['figures_manifest_sha256'] == sha(demo / 'figures/manifest.json')
    fm = json.loads((demo / 'figures/manifest.json').read_text())
    for name in ['KMeans_K10.png', 'HDBSCAN_R.png']:
        assert sha(demo / 'figures' / name) == fm['files'][name]
        code(f"review_image(REVIEW_CAMP/'embedding/TKU4163/hvg2000/PCA2/figures/{name}')")
    spec = json.loads((CAMP / 'protocol/embedding.json').read_text())
    records = []
    for sample in spec['samples']:
        for budget in spec['budgets']:
            for space in spec['spaces']:
                unit = CAMP / 'embedding' / sample / budget / space
                if not checked(unit, 'fit_manifest.json', 'FIT_COMPLETE') or not checked(unit / 'evaluation') or not checked(unit / 'figures'):
                    continue
                fm = json.loads((unit / 'figures/manifest.json').read_text())
                for condition in fm['conditions']:
                    for name in condition['files']:
                        path = unit / 'figures' / name
                        assert sha(path) == fm['files'][name]
                        records.append(dict(sample=sample, budget=budget, representation=space,
                            candidate=Path(condition['condition']['dest']).name,
                            dl_status=condition['dl_status'], training_executed=condition['training_executed'],
                            figure_relative_to_notebook=link(path), sha256=fm['files'][name]))
    index = OUT / 'all_completed_A1_figures.csv'
    with index.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]));writer.writeheader();writer.writerows(records)
    md(f'''This snapshot indexes {len(records)} PNG/PDF files. Units enter the index after terminal fitting, evaluation and every figure are complete.

{link(index)}

The original six-budget/four-route results and historical Python figures remain in their chapters. Complete new cohort statistics will be integrated after verification.
''')
    chapter = '5'
    md('''### scDeepSort input correction

The previous scDeepSort comparison used raw counts; the official model requires LogNormalize. Input values have been independently checked at small, medium and large sample sizes. The corrected121-sample run is in progress. Old predictions and agreement tables remain historical versions; corrected method comparisons, patient statistics and cross-tool agreement require separate verification.

The pending comparison does not establish universal optimality of the original workflow. Further PTC computation awaits the complete GBM deliverable. The historical SignacX None / zero predicted T-cell entry is retained.
''')
    addition = nbf.v4.new_notebook(cells=cells, metadata={'kernelspec':{
        'name':'dgscrna_hvg','display_name':'DG-scRNA pinned environment','language':'python'}})
    os.environ['JUPYTER_PATH'] = str(ROOT / 'results/hvg_ptc_20260916_v1/jupyter')
    NotebookClient(addition, timeout=300, kernel_name='dgscrna_hvg',
        resources={'metadata':{'path':str(ROOT / 'notebooks')}}, allow_errors=False).execute()
    assert all(o.get('output_type') != 'error' for c in addition.cells for o in c.get('outputs', []))
    assert sha(TARGET) == before, 'Notebook modified during preparation; refuse overwrite'
    inserted = []
    groups = {ch:[c for c in addition.cells if c.metadata['v5_option_a_chapter']==ch] for ch in ('1','4','5','9')}
    ends = {ch:max(i for i,c in enumerate(original) if c.metadata.get('v5_option_a_chapter')==ch) for ch in groups}
    for i, cell in enumerate(original):
        inserted.append(cell)
        for ch in groups:
            if i == ends[ch]:
                inserted.extend(groups[ch])
    book.cells = inserted
    nbf.validate(book)
    temp = TARGET.with_name(TARGET.name + '.reviewer.tmp')
    nbf.write(book, temp);temp.replace(TARGET)
    after = nbf.read(TARGET, as_version=4)
    retained = [c for c in after.cells if TAG not in c.metadata.get('tags', [])]
    assert cellhash(retained) == original_hash
    subprocess.run(['/users/PCON0080/yimin/miniforge3/bin/jupyter','trust',str(TARGET)],check=True)
    result = dict(status='completed', updated_at=at, job=os.environ['SLURM_JOB_ID'],
        step=os.environ.get('SLURM_STEP_ID'), original_cells=len(original), original_cells_sha256=original_hash,
        all_original_cells_and_outputs_unchanged=True, inserted_cells=len(addition.cells),
        total_cells=len(after.cells), input_notebook_sha256=before, notebook_sha256=sha(TARGET),
        backup=str(backup), all_A1_figure_index_rows=len(records),
        models_rerun=False, scope='Verified marker expression and first A1 examples; full GBM analyses remain in progress')
    (OUT / 'manifest.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(result,ensure_ascii=False),flush=True)
finally:
    lock.rmdir()
