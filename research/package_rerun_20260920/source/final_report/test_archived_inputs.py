"""Retrospective function tests only; never impersonates a completed new campaign."""
from pathlib import Path
import json
import os
import shutil
import sys
import hashlib

assert os.environ.get('SLURM_JOB_ID')
CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
import run as report
import revalidate
from render_B import render
import numpy as np
import pandas as pd
import nbformat

ROOT = report.ROOT
OUT = report.CAMP / 'final_report_validation' / (sys.argv[1] if len(sys.argv) == 2 else 'archived_functions_v1')
OUT.mkdir(parents=True)
before = report.sha(ROOT / 'notebooks/dgscrna_results.ipynb')
bound_inputs = {}
report.old_inputs(bound_inputs)
read = lambda p: report.read_table(p, True)
a2 = report.OLD / 'no_clustering_lfine_v1/summary'
old_candidates = read(a2 / 'all_candidate_metrics.csv')
old_candidates = old_candidates[old_candidates.stage.isin(report.STAGES)].copy()
old_anchors = read(a2 / 'all_original_route_metrics.csv')
old_anchors = old_anchors[old_anchors.stage.isin(report.STAGES)].copy()
choices, pairs, stats = [read(a2 / name) for name in ['training_patient_lambda_choices.csv', 'patient_paired_results.csv', 'patient_paired_summary.csv']]
folds = read(report.OLD / 'no_clustering/patient_folds.csv')
a2out = OUT / 'A2'; a2out.mkdir()
a2proof, newchoices = revalidate.revalidate_a2(old_candidates.copy(), old_anchors.copy(), old_candidates, old_anchors,
    choices, pairs, stats, folds, a2out)
badpairs = pairs.copy()
badpairs.loc[0, 'lfine_macroF1_cellwise'] = np.nextafter(badpairs.loc[0, 'lfine_macroF1_cellwise'], np.inf)
try:
    revalidate.revalidate_a2(old_candidates.copy(), old_anchors.copy(), old_candidates, old_anchors,
        choices, badpairs, stats, folds, OUT / 'must_not_exist_A2')
except RuntimeError as error:
    assert 'Exact inference vector differs' in str(error), str(error)
else:
    raise RuntimeError('A2 one-ULP patient tamper passed')

c = report.OLD / 'comparison_lfine_v1/summary'
candidates = read(c / 'all_eligible_candidate_metrics.csv.gz')
fresh = candidates[candidates.method.eq('DG-scRNA') & candidates.primary_candidate].copy()
cchoices, heldout, summary, contrasts = [read(c / name) for name in ['training_patient_choices.csv',
    'patient_heldout_results.csv', 'patient_heldout_summary.csv', 'paired_patient_comparisons.csv']]
cfolds = read(ROOT / 'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/protocol/patient_folds.csv')
cout = OUT / 'C'; cout.mkdir()
cproof, comparison = revalidate.revalidate_c(fresh, candidates, cchoices, heldout, summary, contrasts, cfolds, cout)
badheldout = heldout.copy()
ix = badheldout[badheldout.method.eq('DG-scRNA')].index[0]
badheldout.loc[ix, 'lfine_macroF1'] = np.nextafter(badheldout.loc[ix, 'lfine_macroF1'], np.inf)
try:
    revalidate.revalidate_c(fresh, candidates, cchoices, badheldout, summary, contrasts, cfolds, OUT / 'must_not_exist_C')
except RuntimeError as error:
    assert 'Exact inference vector differs' in str(error), str(error)
else:
    raise RuntimeError('C one-ULP patient tamper passed')

original = read(report.OLD / 'controls_lfine_v1/full/terminal_metrics.csv.gz')
tasks = json.loads((CODE.parent / 'extension_manager/CAMPAIGN_CANDIDATE_FULL_v3.tasks.json').read_text())
data = original.copy()
for index, row in data.iterrows():
    candidates_tasks = [t for t in tasks if t['family'] == row.task and t['sample'] == row['sample'] and t['budget'] == row.budget]
    if row.task == 'neighbors':
        candidates_tasks = [t for t in candidates_tasks if t['configuration']['snn_k'] == row.snn_k
            and t['configuration']['umap_neighbors'] == row.umap_neighbors
            and row.route == t['configuration']['space'] + '_' + t['configuration']['method']]
        assert len(candidates_tasks) == 1, (row.to_dict(), candidates_tasks)
        task = candidates_tasks[0]
        data.loc[index, 'configuration'] = task['configuration']['name']
        data.loc[index, 'family'] = 'neighbor_control'
    else:
        candidates_tasks = [t for t in candidates_tasks if t['configuration']['model_seed'] == row.model_seed]
        assert len(candidates_tasks) == 1
        task = candidates_tasks[0]
        data.loc[index, 'configuration'] = f'seed{int(row.model_seed)}_epoch{int(row.epochs):02d}'
        data.loc[index, 'family'] = 'learning_control'
    data.loc[index, 'campaign_task'] = task['index']
mapped = report.b_rows(data.drop(columns=['task', 'control', 'snn_k', 'umap_neighbors', 'model_seed', 'epochs', 'row_key']), tasks)
for name in ['task', 'snn_k', 'umap_neighbors', 'model_seed', 'epochs', 'lfine_macroF1', 'coverage']:
    assert np.array_equal(original[name].to_numpy(), mapped[name].to_numpy(), equal_nan=True) if name != 'task' else original[name].tolist() == mapped[name].tolist()
descriptions, checkpoints = render(mapped, OUT / 'B')
assert len(descriptions) == 48 and len(checkpoints) == 48
assert descriptions.default_relation.eq('lower_than_an_alternative').sum() == 26

# Presentation-only fixture: the canonical notebook and real completion gates stay untouched.
texts = report.snippets(descriptions, stats, newchoices, comparison, OUT)
fixture = OUT / 'notebook_fixture'
(fixture / 'notebooks').mkdir(parents=True)
source = nbformat.read(ROOT / 'notebooks/dgscrna_results.ipynb', as_version=4)
source.metadata['packaged_core_rerun'] = {'fixture': True}
source.metadata['packaged_A1_selection'] = {'fixture': True}
nbformat.write(source, fixture / 'notebooks/dgscrna_results.ipynb')
archive = fixture / 'old'; (archive / 'notebook').mkdir(parents=True)
report.ROOT, report.OLD = fixture, archive
report.write(OUT / 'manifest.json', {'status': 'archived_function_fixture_only'})
try:
    report.update_notebook(OUT, texts, False, expected_before='mismatched_digest')
except RuntimeError as error:
    assert 'Notebook changed during report revalidation' in str(error)
else:
    raise RuntimeError('Concurrent notebook mismatch accepted')
report.update_notebook(OUT, texts, False, expected_before=report.sha(fixture / 'notebooks/dgscrna_results.ipynb'))
candidate = nbformat.read(OUT / 'dgscrna_results_candidate.ipynb', as_version=4)
assert len(candidate.cells) == len(source.cells) == 20
assert [i for i, (old, new) in enumerate(zip(source.cells, candidate.cells)) if old != new] == [10, 15, 16, 17]
assert report.sha(ROOT / 'notebooks/dgscrna_results.ipynb') == before
result = dict(status='passed_archived_function_regression', new_campaign_results_claimed=False,
    A2=a2proof, C=cproof, B_rows=len(mapped), B_descriptive_sweeps=len(descriptions), B_checkpoint_rows=len(checkpoints),
    one_ULP_inferential_vector_tamper_rejected=['A2', 'C'], notebook_cells=20, changed_fixture_cells=[10,15,16,17],
    bound_original_input_files=len(bound_inputs), concurrent_notebook_change_rejected=True,
    canonical_notebook_unchanged=True, canonical_notebook_sha256=before,
    source_hashes={str(p): report.sha(p) for p in CODE.glob('*.py')},
    job=os.environ['SLURM_JOB_ID'], step=os.environ.get('SLURM_STEP_ID'))
report.write(OUT / 'verification.json', result)
print(json.dumps(result, indent=2))
