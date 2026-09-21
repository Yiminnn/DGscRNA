#!/usr/bin/env python3
"""Independent A1 adapter pilot comparison; old outputs are verification-only."""
from pathlib import Path
import argparse
import importlib.util
import json
import os
import subprocess
import sys

CODE = Path(__file__).resolve().parent
ROOT = CODE.parents[2]
SPACES = ['noDR', 'PCA2', 'FA2', 'ICA2', 'Isomap2', 'UMAP2', 'TSNE2']
PARTITIONS = [f'{method}_K{k:02d}' for method in ['KMeans', 'GMM'] for k in [5, 10, 15, 20, 30, 40]] + ['HDBSCAN_R']


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def main():
    assert os.environ.get('SLURM_JOB_ID')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--actual', type=Path, required=True, help='New sample/budget directory')
    parser.add_argument('--expected', type=Path, required=True, help='Historical A1 sample/budget; verification only')
    parser.add_argument('--rscript', required=True)
    parser.add_argument('--evaluation-python', required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--spaces', nargs='+', choices=SPACES, default=SPACES)
    args = parser.parse_args()
    spaces = args.spaces
    assert len(spaces) == len(set(spaces))
    import numpy as np
    import pandas as pd
    helper = module('independent_packaged_terminal_check', CODE.parent / 'verify_packaged_parity.py')
    adapter = module('adapter_contract_only', CODE / 'run.py')
    source_manifest = adapter.verify_source()
    sha = helper.sha
    actual, expected, out = args.actual.resolve(), args.expected.resolve(), args.out.resolve()
    assert not out.exists(), 'Preserve any previous independent verification'
    out.mkdir(parents=True)
    pairs, geometry_checks, cases, pending, sources = [], [], [], [], {}
    for name in ['scaled_HVG.float64.bin', 'features.txt', 'cells.csv']:
        assert sha(actual / 'geometry' / name) == sha(expected / 'geometry' / name), name
        geometry_checks.append(dict(artifact=name, exact_bytes=True, sha256=sha(actual / 'geometry' / name)))
    evaluation_proofs = {}
    for space in spaces:
        new, old = actual / space, expected / space
        assert adapter.checked(new, 'fit_manifest.json', 'FIT_COMPLETE')
        assert adapter.checked(old, 'fit_manifest.json', 'FIT_COMPLETE')
        cfg = json.loads((new / 'config.json').read_text())
        fit = json.loads((new / 'fit_manifest.json').read_text())
        assert fit['config'] == cfg and fit['installed_package_terminal'] and not fit['old_models_or_caches_used']
        assert cfg['space'] == space and len(cfg['conditions']) == 13
        assert fit['source_sha256'] == sha(CODE / 'run.py')
        for folder in [new, old]:
            assert adapter.checked(folder, 'representation.json', 'REPRESENTATION_COMPLETE')
        representation = json.loads((new / 'representation.json').read_text())
        if space != 'noDR':
            left = pd.read_csv(new / 'embedding.csv', dtype={'cell_id': str}, keep_default_na=False)
            right = pd.read_csv(old / 'embedding.csv', dtype={'cell_id': str}, keep_default_na=False)
            pd.testing.assert_frame_equal(left, right, check_exact=True)
            assert sha(new / 'embedding.csv') == representation['embedding_sha256']
            geometry_checks.append(dict(space=space, coordinates_exact=True))
        if space == 'ICA2':
            assert representation['comparator_policy_id'] == 'ICA2_parallel_caps_then_deflation_v2'
            policy = json.loads((CODE / 'protocol/embedding_convergence_repair_20260920_v2.json').read_text())
            attempts = representation['attempts']
            assert 1 <= len(attempts) <= 4 and attempts[-1]['accepted']
            assert not any(attempt['accepted'] for attempt in attempts[:-1])
            for attempt, planned in zip(attempts, policy['attempts']):
                assert attempt['algorithm'] == planned['algorithm'] and attempt['cap'] == planned['max_iter']
                assert sha(attempt['path']) == attempt['sha256']
            assert representation['actual_solver'] == attempts[-1]['algorithm']
            assert representation['actual_iteration_cap'] == attempts[-1]['cap']
            assert representation['fallback_used'] == (attempts[-1]['algorithm'] == 'deflation')
            if representation['fallback_used']:
                residuals = representation['fixed_point_residuals']
                assert len(residuals) == 2 and np.isfinite(residuals).all() and max(residuals) < 1e-4
            assert representation['policy_sha256'] == sha(CODE / 'protocol/embedding_convergence_repair_20260920_v2.json')
        for name in PARTITIONS:
            newp, oldp = new / name, old / name
            for folder in [newp, oldp]:
                assert adapter.checked(folder, 'partition_manifest.json', 'PARTITION_COMPLETE')
                pm = json.loads((folder / 'partition_manifest.json').read_text())
                assert sha(folder / 'clusters.csv') == pm['clusters_sha256']
                assert adapter.checked(folder, 'score_manifest.json', 'SCORE_COMPLETE')
                sm = json.loads((folder / 'score_manifest.json').read_text())
                assert sm['reference_labels_used_for_fit'] is False
                assert sha(folder / 'initial_calls.csv.gz') == sm['initial_sha256']
            for filename in ['clusters.csv', 'cells.csv', 'initial_calls.csv.gz', 'cluster_calls.csv.gz', 'marker_retention.csv.gz']:
                left = pd.read_csv(newp / filename, dtype=str, keep_default_na=False)
                right = pd.read_csv(oldp / filename, dtype=str, keep_default_na=False)
                pd.testing.assert_frame_equal(left, right, check_exact=True)
            nt = newp / 'terminal/L00_mean'
            terminal = json.loads((nt / 'terminal_manifest.json').read_text())
            assert Path(terminal['cache_directory']).resolve().is_relative_to(new / 'DL_cache')
            training = json.loads((nt / 'training_manifest.json').read_text())
            assert Path(training['provenance']['first_condition']).resolve().is_relative_to(new)
            assert training['provenance']['reference_labels_used_for_fit'] is False
            assert training['provenance']['input_signature']['torch_threads'] == 4
            pairs.append(dict(space=space, partition=name, actual=str(newp), expected=str(oldp)))
            pending.append((space, name, nt, oldp / 'terminal/L00_mean', cfg))
    # Verify all R DEG tables/density matrices in one isolated R process.
    adapter.write_json(out / 'R_pairs.json', pairs)
    env = os.environ.copy()
    for key in ['R_LIBS', 'DGSCRNA_REFERENCE_R_LIB']:
        env.pop(key, None)
    env.update(R_LIBS_USER='', R_LIBS_SITE='', R_ENVIRON_USER='/dev/null', R_PROFILE_USER='/dev/null')
    subprocess.run([args.rscript, '--vanilla', str(CODE / 'verify_R.R'), str(out / 'R_pairs.json'),
        str(out / 'R_comparison.json')], check=True, env=env)
    r_results = json.loads((out / 'R_comparison.json').read_text())
    assert len(r_results) == 13 * len(spaces) and all(row['DEG_exact'] and row['density_exact'] for row in r_results)
    terminal_checks = []
    by_space = {space: [] for space in spaces}
    configs = {}
    for space, name, actual_terminal, expected_terminal, cfg in pending:
        terminal = json.loads((actual_terminal / 'terminal_manifest.json').read_text())
        record = dict(route=f'{space}_{name}', library='CM2_glioma_other', cutoff='mean',
            dl_status=terminal['dl_status'], training_executed=terminal['training_executed'])
        case = helper.compare_terminal(actual_terminal, expected_terminal, record,
            'L00_mean', 'L00_mean', terminal_checks, fresh=False)
        by_space[space].append(case)
        cases.append(case)
        configs[space] = cfg
    source_hashes = {str(CODE / name): sha(CODE / name) for name in
        ['run.py', 'adaptive_ica.py', 'export_geometry.R', 'native_geometry.R', 'score_candidates.R',
         'SOURCE_PROTOCOL_MANIFEST.json', 'verify.py', 'verify_R.R']}
    source_hashes |= {str(CODE / 'protocol' / name): digest for name, digest in source_manifest['frozen_protocols'].items()}
    for name in ['evaluate_extension_lfine.py', 'verify_packaged_parity.py']:
        source_hashes[str(CODE.parent / name)] = sha(CODE.parent / name)
    assert source_hashes[str(CODE.parent / 'evaluate_extension_lfine.py')] == '424fb152cee5d95990b0eae724761601199230cb55512e22bd8de5cfcdfd25a9'
    for space in spaces:
        dest = out / space
        dest.mkdir()
        cfg = configs[space]
        proof = dict(status='passed_exact', sample=cfg['sample'], budget=cfg['budget'], space=space,
            conditions=by_space[space], terminal_conditions=13, coordinates_exact=True,
            partition_labels_exact=True, DEG_density_exact=True, installed_package_DL_exact=True,
            old_models_or_caches_used=False, reference_labels_used_for_fit=False,
            source_hashes=source_hashes, actual=str(actual / space), expected=str(expected / space),
            actual_fit_manifest_sha256=sha(actual / space / 'fit_manifest.json'), job=os.environ['SLURM_JOB_ID'])
        adapter.write_json(dest / 'parity.json', proof)
        specification = dict(sample=cfg['sample'], budget=cfg['budget'], family='seven_space_control',
            configuration=space, expected_conditions=13, source_hashes=source_hashes,
            parity_receipt=dict(path=str(dest / 'parity.json'), sha256=sha(dest / 'parity.json')),
            conditions=[dict(route=f'{space}_{name}', library='CM2_glioma_other', cutoff='mean',
                terminal_directory=str(actual / space / name / 'terminal/L00_mean')) for name in PARTITIONS])
        adapter.write_json(dest / 'evaluation_spec.json', specification)
        subprocess.run([args.evaluation_python, '-s', str(CODE.parent / 'evaluate_extension_lfine.py'),
            '--spec', str(dest / 'evaluation_spec.json'), '--out', str(dest / 'evaluation')], check=True)
        assert adapter.checked(dest / 'evaluation')
        evaluation = json.loads((dest / 'evaluation/manifest.json').read_text())
        assert evaluation['n_conditions'] == 13 and evaluation['n_threshold_rows'] == evaluation['n_valid'] == 26
        assert sha(dest / 'evaluation/metrics.csv.gz') == evaluation['outputs']['metrics.csv.gz']
        evaluation_proofs[space] = dict(parity_path=str(dest / 'parity.json'),
            parity_sha256=sha(dest / 'parity.json'), specification_path=str(dest / 'evaluation_spec.json'),
            specification_sha256=sha(dest / 'evaluation_spec.json'),
            manifest_path=str(dest / 'evaluation/manifest.json'), manifest_sha256=sha(dest / 'evaluation/manifest.json'),
            metrics_path=str(dest / 'evaluation/metrics.csv.gz'), metrics_sha256=sha(dest / 'evaluation/metrics.csv.gz'))
    result = dict(status='passed_exact', sample=configs[spaces[0]]['sample'], budget=configs[spaces[0]]['budget'],
        spaces=spaces, n_representations=len(spaces), n_partitions=13 * len(spaces),
        n_terminal_conditions=13 * len(spaces), n_Lfine_threshold_rows=26 * len(spaces),
        geometry_checks=geometry_checks, terminal_conditions=cases, R_comparison_sha256=sha(out / 'R_comparison.json'),
        terminal_checks=terminal_checks, source_hashes=source_hashes,
        evaluation_proofs=evaluation_proofs,
        scope='Only the listed sample/budget/representation units have been independently verified; this is not a full1694-unit completion claim',
        ICA_scope='The TKU4163/HVG2000 pilot exercises the successful canonical parallel5000 branch; other branches require their own observed-unit verification',
        native_UMAP='Direct scaled-HVG uwot comparator, not the default PCA30-to-UMAP workflow',
        old_models_or_caches_used=False, job=os.environ['SLURM_JOB_ID'], completed_at=adapter.utc())
    adapter.write_json(out / 'manifest.json', result)
    adapter.complete(out)
    print(json.dumps({k: result[k] for k in ['status', 'n_representations', 'n_partitions', 'n_terminal_conditions', 'n_Lfine_threshold_rows']}, indent=2))


if __name__ == '__main__':
    main()
