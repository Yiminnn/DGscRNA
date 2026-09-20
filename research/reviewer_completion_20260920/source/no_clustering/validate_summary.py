"""Read-only acceptance checks for frozen GBM A2 results; run after aggregation."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import sys

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT = ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/no_clustering'
CODE = OUT/'source_v2'
sys.path.insert(0, str(CODE))
from a2_common import verify_source, verify_terminal, checked, sha, write_json, require_slurm


def run():
    require_slurm()
    import numpy as np
    import pandas as pd

    protocol = verify_source()
    summary = OUT/'summary'
    assert checked(summary)
    manifest = json.loads((summary/'manifest.json').read_text())
    assert manifest['source_bundle_sha256'] == protocol['source_bundle_sha256']
    for name, digest in manifest['outputs'].items():
        assert sha(summary/name) == digest, name
    assert manifest['n_samples'] == 121 and manifest['n_candidate_conditions'] == 1210
    assert manifest['n_sample_result_rows_by_cohort'] == {'primary': 388, 'all': 484}
    assert manifest['n_contrasts'] == 32

    metrics = pd.read_csv(summary/'all_candidate_metrics.csv')
    assert len(metrics) == 3630
    identity = ['sample', 'budget', 'arm_id', 'stage']
    assert not metrics.duplicated(identity).any()
    assert set(metrics['sample']) == set(protocol['samples'])
    assert set(metrics['stage']) == {'marker_only', 'terminal090', 'terminal070'}
    assert metrics.primary.isin([True, False]).all()
    for field in ['macroF1_present', 'macroF1_fixed11', 'weightedF1', 'accuracy', 'coverage', 'unknown_rate', 'mapped_coverage']:
        assert np.isfinite(metrics[field]).all() and metrics[field].between(0, 1).all(), field
    np.testing.assert_allclose(metrics.coverage + metrics.unknown_rate, 1, atol=1e-12, rtol=0)

    input_hashes = {}
    cell_counts = {}
    for sample in protocol['samples']:
        seeds = []
        for budget in protocol['budgets']:
            unit = OUT/'GBM'/sample/budget
            cfg = json.loads((unit/'config.json').read_text())
            assert checked(unit, 'fit_manifest.json', 'FIT_COMPLETE')
            fit = json.loads((unit/'fit_manifest.json').read_text())
            assert fit['input_signature'] == cfg['input_signature']
            assert fit['reference_labels_used_for_fit'] is False
            evaluation = unit/'evaluation'
            assert checked(evaluation)
            em = json.loads((evaluation/'manifest.json').read_text())
            for name, digest in em['outputs'].items():
                assert sha(evaluation/name) == digest
            assert em['every_candidate_plotted'] and em['all_cell_denominator']
            assert em['n_metric_rows'] == 15 and em['n_candidates'] == 5
            n_cells = protocol['inputs'][sample+'/'+budget]['n_cells']
            assert em['n_cells'] == n_cells
            frame = metrics[(metrics['sample'] == sample) & (metrics.budget == budget)]
            assert len(frame) == 15 and frame.n_cells.eq(n_cells).all()
            seed = pd.read_csv(unit/'cellwise_seed/initial_calls.csv.gz', dtype=str, keep_default_na=False)
            assert len(seed) == n_cells and not seed.cell_id.duplicated().any()
            seeds.append(seed)
            for arm in protocol['arms']:
                aid = arm['id']
                tm = verify_terminal(unit/'cellwise_seed', aid)
                assert tm['DL_features'] == protocol['inputs'][sample+'/'+budget]['DL_features']
                with np.load(unit/'cellwise_seed/terminal'/aid/'terminal.npz') as z:
                    assert np.array_equal(z['initial'], seed[aid])
                    known = z['initial'] != 'Undecided'
                    for final in ['final090', 'final070']:
                        assert len(z[final]) == n_cells
                        assert np.array_equal(z[final][known], z['initial'][known])
                assert fit['terminal_manifests'][aid] == sha(unit/'cellwise_seed/terminal'/aid/'terminal_manifest.json')
            input_hashes[sample+'/'+budget] = sha(evaluation/'manifest.json')
            cell_counts[sample] = n_cells
        assert seeds[0].equals(seeds[1]), f'All-RNA seeds differ across HVG budgets: {sample}'

    folds = pd.read_csv(OUT/'patient_folds.csv')
    assert sha(OUT/'patient_folds.csv') == protocol['patient_folds_sha256']
    terminal = metrics[metrics.stage == 'terminal090'].merge(folds[['sample', 'patient', 'fold']], on=['sample', 'patient'], validate='many_to_one')
    choices = pd.read_csv(summary/'training_patient_lambda_choices.csv')
    assert len(choices) == 20
    assert choices.test_labels_used_for_selection.eq(False).all()
    selected = pd.read_csv(summary/'selected_and_fixed_sample_results.csv')
    assert selected.groupby('cohort').size().to_dict() == {'all': 484, 'primary': 388}
    for row in choices.itertuples():
        eligible = terminal[terminal.primary.eq(True)] if row.cohort == 'primary' else terminal
        eligible = eligible[eligible.budget == row.budget]
        training = eligible[eligible.fold != row.fold]
        patient = training.groupby(['patient', 'lambda_value']).macroF1_present.mean()
        scores = patient.groupby('lambda_value').mean().reset_index(name='score')
        scores['distance'] = (scores.lambda_value - 1).abs()
        expected = scores.sort_values(['score', 'distance', 'lambda_value'], ascending=[False, True, True], kind='stable').iloc[0]
        assert row.lambda_value == expected.lambda_value
        np.testing.assert_allclose(row.training_patient_macroF1, expected.score, rtol=0, atol=1e-12)
        held = selected[(selected.cohort == row.cohort) & (selected.budget == row.budget) & (selected.selection == 'training_patient_selected') & (selected.fold == row.fold)]
        expected_samples = set(eligible.loc[eligible.fold == row.fold, 'sample'])
        assert set(held['sample']) == expected_samples
        assert held.lambda_value.eq(row.lambda_value).all()
    assert selected.loc[selected.selection == 'fixed_lambda1', 'lambda_value'].eq(1).all()

    paired = pd.read_csv(summary/'patient_paired_summary.csv')
    assert len(paired) == 32 and paired.groupby('cohort').size().eq(16).all()
    assert paired.loc[paired.cohort == 'primary', 'n_patients'].eq(55).all()
    assert paired.loc[paired.cohort == 'all', 'n_patients'].eq(59).all()
    assert paired.p_value.between(0, 1).all() and paired.p_holm.between(0, 1).all()
    assert (paired.p_holm >= paired.p_value - 1e-12).all()
    assert (paired.CI025 <= paired.CI975).all()
    np.testing.assert_allclose(paired.mean_cellwise_macroF1 - paired.mean_original_macroF1,
                               paired.mean_delta_cellwise_minus_original, rtol=0, atol=1e-12)
    largest = sorted(cell_counts.items(), key=lambda item: (-item[1], item[0]))[0]
    report = dict(status='passed', scope='GBM A2 output and selection acceptance checks; not whole work-package A',
                  n_samples=121, n_evaluation_units=242, n_candidates=1210, n_metric_rows=3630,
                  every_candidate_preserved_and_plotted=True, all_cells_preserved=True,
                  all_RNA_seed_identity_across_HVG_budgets=True, terminal_known_seeds_preserved=True,
                  patient_selector_recomputed_from_training_patients=True, n_contrasts=32,
                  largest_sample={'sample': largest[0], 'n_cells': largest[1]},
                  source_bundle_sha256=protocol['source_bundle_sha256'],
                  validation_source_sha256=sha(Path(__file__)), summary_manifest_sha256=sha(summary/'manifest.json'),
                  evaluation_manifest_hashes=input_hashes, whole_work_package_A_complete=False,
                  job=os.environ['SLURM_JOB_ID'], completed_at=datetime.now(timezone.utc).isoformat())
    write_json(OUT/'validation.json', report)
    print(json.dumps({key: value for key, value in report.items() if key != 'evaluation_manifest_hashes'}, indent=2), flush=True)


if __name__ == '__main__':
    run()
