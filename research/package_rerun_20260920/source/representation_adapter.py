"""Fresh native-R representation controls using the released DL backend.

Runs the archived R numerical bodies with only explicit input/output paths,
isolated libraries and identifier-safe CSV reads. No truth enters fitting.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import ast
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import traceback

HERE = Path(__file__).resolve().parent
ROUTES = ['PCA30_SNN', 'PCA30_HDBSCAN_R', 'UMAP2_SNN', 'UMAP2_HDBSCAN_R']


def need(value, message):
    if not value:
        raise RuntimeError(message)


def load(path):
    return json.loads(Path(path).read_text())


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def write_new(path, value):
    with Path(path).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def check_sources():
    lock = load(HERE / 'representation_sources.json')
    need(lock['status'] == 'frozen_for_bounded_validation', 'Representation source lock is not frozen')
    for name, record in lock['files'].items():
        need(sha(record['path']) == record['sha256'], 'Frozen source changed: ' + name)
    derivation = load(HERE / 'representation_derivation.json')
    for name in ['prepare_representation_R.R', 'score_representation_R.R']:
        recipe = derivation[name]
        text = Path(recipe['original_path']).read_text()
        need(sha(recipe['original_path']) == recipe['original_sha256'], 'Archived R source changed')
        for edit in recipe['replacements']:
            need(text.count(edit['before']) == 1, 'Unrecognized R source transformation')
            text = text.replace(edit['before'], edit['after'])
        need(text == (HERE / name).read_text() and sha(HERE / name) == recipe['portable_sha256'],
             'R numerical body changed beyond its recorded IO/library/identifier edits')
    recipe = derivation['configuration_function']
    original = Path(recipe['original_path']).read_text()
    portable = (HERE / 'representation_configurations.py').read_text()
    def function(text):
        node = next(n for n in ast.parse(text).body if isinstance(n, ast.FunctionDef) and n.name == 'configurations')
        return ast.get_source_segment(text, node)
    need(function(original) == function(portable), 'Frozen scientific configuration function changed')
    audit = load(lock['files']['roundtrip_audit']['path'])
    need(audit['status'] == 'passed_exact_annotation_with_explained_membership_roundtrip' and
         audit['archived_representation_all_diagnostics_exact'] is True and
         audit['all48_seeds_DEG_density_exact'] is True,
         'Default core diagnostic exception lacks its independent causal proof')
    return lock


def execute(command, environment, log, expected_text=None):
    with log.open('x') as handle:
        result = subprocess.run(command, env=environment, cwd='/tmp', stdout=handle, stderr=subprocess.STDOUT)
    need(result.returncode == 0, 'Scientific subprocess failed; inspect ' + str(log))
    if expected_text:
        need(expected_text in log.read_text(), 'Missing scientific verification marker: ' + str(log))


def score_parity(verifier, actual, expected, label, runtime, environment, logs, check, diagnostic_roundtrip=None):
    verifier.receipt(actual, 'score_manifest.json', 'SCORE_COMPLETE')
    verifier.receipt(expected, 'score_manifest.json', 'SCORE_COMPLETE')
    am, em = load(actual / 'score_manifest.json'), load(expected / 'score_manifest.json')
    aa, ea = verifier.arm_map(am), verifier.arm_map(em)
    need(set(aa) == set(ea) and len(aa) == 48, 'Expected all48 marker/cutoff seed conditions')
    diagnostic = None
    for name in ['cells.csv', 'clusters.csv', 'density_diagnostics.csv']:
        need((actual / name).exists() == (expected / name).exists(), 'Clustering artifact presence differs')
        if (actual / name).exists():
            left, right = verifier.frame(actual / name), verifier.frame(expected / name)
            if name == 'density_diagnostics.csv' and diagnostic_roundtrip is not None:
                import numpy as np
                audit = Path(diagnostic_roundtrip)
                verifier.frame_equal(left, verifier.frame(audit / 'rds_membership.csv'), label + '/audited_RDS_membership', check)
                verifier.frame_equal(right, verifier.frame(audit / 'csv_membership.csv'), label + '/audited_core_CSV_membership', check)
                verifier.frame_equal(left.drop(columns='membership'), right.drop(columns='membership'),
                                     label + '/cell_and_noise_identity', check)
                delta = np.abs(left.membership.astype(float).to_numpy() - right.membership.astype(float).to_numpy())
                diagnostic = dict(comparison='deterministic_CSV_roundtrip_expected_difference', membership_exact=False,
                                  changed_cells=int((delta > 0).sum()), max_abs=float(delta.max()),
                                  RDS_and_CSV_memberships_each_exact_to_independent_recomputation=True,
                                  causal_audit_sha256=sha(audit / 'verification.json'))
            else:
                verifier.frame_equal(left, right, label + '/' + name, check)
    ai, ei = verifier.frame(actual / 'initial_calls.csv.gz'), verifier.frame(expected / 'initial_calls.csv.gz')
    verifier.frame_equal(ai[['cell_id']], ei[['cell_id']], label + '/initial_cell_order', check)
    ac, ec = verifier.frame(actual / 'cluster_calls.csv.gz'), verifier.frame(expected / 'cluster_calls.csv.gz')
    ar, er = verifier.frame(actual / 'marker_retention.csv.gz'), verifier.frame(expected / 'marker_retention.csv.gz')
    for key, (_, arm) in aa.items():
        _, reference_arm = ea[key]
        verifier.frame_equal(ai[[arm['seed_column']]].set_axis(['initial'], axis=1),
                             ei[[reference_arm['seed_column']]].set_axis(['initial'], axis=1),
                             label + '/' + str(key) + '/initial', check)
        left = ac[(ac.library == key[0]) & (ac.cutoff == key[1])].drop(columns='arm_id').reset_index(drop=True)
        right = ec[(ec.library == key[0]) & (ec.cutoff == key[1])].drop(columns='arm_id').reset_index(drop=True)
        verifier.frame_equal(left, right, label + '/' + str(key) + '/cluster_calls', check)
    libraries = sorted({key[0] for key in aa})
    for library in libraries:
        verifier.frame_equal(ar[ar.library == library].reset_index(drop=True), er[er.library == library].reset_index(drop=True),
                             label + '/' + library + '/marker_retention', check)
    path = logs / (label + '_libraries.json')
    write_new(path, libraries)
    log = logs / (label + '_DEG_density.log')
    execute([runtime['rscript'], '--vanilla', str(HERE / 'verify_packaged_r_artifacts.R'),
             str(actual), str(expected), str(path)], environment, log, 'R_DEG_AND_NAMED_DENSITY_OBJECTS_EXACT')
    return dict(reference=str(expected), initial_arms_exact=48, density_libraries_exact=16,
                DEG_exact=True, diagnostic_roundtrip=diagnostic, R_verification_log_sha256=sha(log))


def run(args):
    need(os.environ.get('SLURM_JOB_ID'), 'Representation fitting/evaluation requires SLURM')
    need(int(os.environ.get('SLURM_CPUS_PER_TASK', 0)) >= 4, 'Four allocated CPUs required')
    need(sys.flags.no_user_site and not sys.flags.optimize and not os.environ.get('PYTHONPATH'),
         'Use installed Python -s without PYTHONPATH or -O')
    os.chdir('/tmp')
    lock = check_sources()
    worker = module('_representation_core_gate', HERE / 'core_task.py')
    gate, _, package = worker.validate_gate(args.gate)
    gate_hash = sha(args.gate)
    protocol = module('_archived_representation_config', HERE / 'representation_configurations.py')
    if args.mode == 'default':
        configurations = protocol.configurations('TKU4163', 'hvg2000', include_defaults=True)
    elif args.mode == 'changed':
        configurations = [c for c in protocol.configurations('TKU4163', 'hvg2000')
                          if c['space'] == 'UMAP10' and c['method'] == 'HDBSCAN_R']
    else:
        need(args.index is not None, 'Task index required')
        manifest_path = HERE / 'representation_tasks.json'
        need(sha(manifest_path) == lock['task_manifest_expected_sha256'], 'Representation task manifest changed')
        tasks = load(manifest_path)
        need(len(tasks) == 180 and 0 <= args.index < len(tasks), 'Invalid180-group task roster/index')
        configurations = [tasks[args.index]['configuration']]
    if args.mode != 'default':
        need(args.default_proof is not None, 'Default parity must pass before varied representation fitting')
        default = load(args.default_proof)
        need(default['status'] in {'passed', 'passed_exact'} and default['terminal_conditions'] == 8 and
             default['mode'] == 'default' and default['source_lock_sha256'] == sha(HERE / 'representation_sources.json') and
             default['gate_sha256'] == gate_hash, 'Default parity belongs to another source/protocol')
    sample, budget = configurations[0]['sample'], configurations[0]['budget']
    source = (args.source or Path(gate['output_root']) / 'core' / sample / budget).resolve()
    source_helper = module('_representation_fresh_source_validator', HERE / 'mlp_controls_adapter.py')
    prep, snapshots = source_helper.verify_source(source, gate, gate_hash, dict(sample=sample, budget=budget, routes=ROUTES))
    markers = source / 'markers/libraries.json'
    need(load(markers) == load(gate['markers']['path']) and len(load(markers)) == 16, 'Marker contents differ from frozen core')
    pm = load(prep / 'prepare_manifest.json')
    need(pm['expression_sha256'] == sha(prep / 'expression_PCA30.rds'), 'Fresh source expression changed')
    for path in [prep / 'expression_PCA30.rds', markers]:
        snapshots[str(path)] = sha(path)
    from dgscrna.reference.runner import doctor
    need(doctor(gate['runtime']['rscript'], gate['runtime'].get('reference_r_lib')) == gate['runtime']['fingerprint'],
         'Installed runtime differs from frozen release')
    out = args.out.resolve()
    need(out.is_relative_to(Path(gate['output_root']).resolve() / 'extensions') and not out.exists(),
         'Use a fresh extension directory; preserve existing outputs')
    out.mkdir(parents=True)
    logs = out / 'logs'; logs.mkdir()
    environment = os.environ.copy()
    for name in ['PYTHONPATH', 'R_LIBS', 'R_LIBS_USER', 'R_LIBS_SITE', 'DGSCRNA_REFERENCE_R_LIB',
                 'DGSCRNA_EXAMPLE_OUT', 'DGSCRNA_ONLY_ROUTE']:
        environment.pop(name, None)
    environment.update(DGSCRNA_REFERENCE_OUT=str(out), DGSCRNA_REQUIRE_SLURM='1',
                       DGSCRNA_DL_CACHE_ROOT=str(out / 'DL_cache'), DGSCRNA_DEG_WORKERS='4',
                       R_ENVIRON_USER=os.devnull, R_PROFILE_USER=os.devnull,
                       OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', PYTHONNOUSERSITE='1')
    if gate['runtime'].get('reference_r_lib'):
        environment['DGSCRNA_REFERENCE_R_LIB'] = gate['runtime']['reference_r_lib']
    config = dict(mode=args.mode, sample=sample, budget=budget, configurations=configurations,
                  source=str(source), gate_sha256=gate_hash, source_lock_sha256=sha(HERE / 'representation_sources.json'),
                  package_root=str(package), immutable_source_hashes=snapshots, job=os.environ['SLURM_JOB_ID'],
                  step=os.environ.get('SLURM_STEP_ID'), started_at=datetime.now(timezone.utc).isoformat())
    write_new(out / 'run_config.json', config)
    try:
        verifier = module('_representation_independent_verifier', HERE / 'verify_packaged_parity.py')
        selected, prepared = [], []
        for scientific in configurations:
            target = out / scientific['name']
            target.mkdir()
            cfg = dict(scientific, source_prepare=str(prep), output=str(target), marker_path=str(markers))
            config_path = target / 'config.json'; write_new(config_path, cfg)
            execute([gate['runtime']['rscript'], '--vanilla', str(HERE / 'prepare_representation_R.R'), str(config_path)],
                    environment, logs / (scientific['name'] + '_prepare.log'))
            execute([gate['runtime']['rscript'], '--vanilla', str(HERE / 'score_representation_R.R'), sample,
                     str(target), str(config_path)], environment, logs / (scientific['name'] + '_score.log'))
            route = scientific['space'] + '_' + scientific['method']
            score = target / route
            sm = load(score / 'score_manifest.json')
            arms = verifier.arm_map(sm)
            need(len(arms) == 48 and sm['reference_labels_used_for_fit'] is False, 'Scoring incomplete or truth used')
            for library in protocol.LIBRARIES:
                aid, _ = arms[(library, 'mean')]
                execute([gate['runtime']['python'], '-s', '-m', 'dgscrna.reference.backend.terminal', str(score), aid],
                        environment, logs / (scientific['name'] + '_' + aid + '_DL.log'))
                td = score / 'terminal' / aid
                tm = load(td / 'terminal_manifest.json')
                need(Path(tm['cache_directory']).resolve().is_relative_to(out / 'DL_cache'), 'External DL cache reused')
                training = load(td / 'training_manifest.json')
                need(Path(training['provenance']['first_condition']).resolve().is_relative_to(out), 'Old cache first condition reused')
                selected.append(dict(scientific=scientific, score=str(score), arm=aid, route=route, library=library,
                                     cutoff='mean', terminal_directory=str(td), dl_status=tm['dl_status'],
                                     training_executed=tm['training_executed']))
            prepared.append((scientific, target, score))
        # All new fits finish before archived numerical results or truth are opened.
        stage_checks, checks, conditions, new_core_conditions, old_core_conditions = [], [], [], [], []
        for scientific, target, score in prepared:
            route = scientific['space'] + '_' + scientific['method']
            old_prep = (Path(gate['reference_root']) / 'GBM' / sample / budget if args.mode == 'default' else
                        Path(gate['reference_root']) / 'GBM_representation_controls' / sample / budget / scientific['name'])
            for name in ['cells.csv', 'selected_features.txt', 'geometry_features.txt', 'scoring_features.txt', 'DL_features.txt']:
                need(sha(target / name) == sha(prep / name) == sha(old_prep / name), 'Native feature/cell input changed: ' + name)
            actual_pm, expected_pm = load(target / 'prepare_manifest.json'), load(old_prep / 'prepare_manifest.json')
            need(actual_pm['DL_binary_sha256'] == pm['DL_binary_sha256'] == expected_pm['DL_binary_sha256'] and
                 sha(actual_pm['DL_binary']) == pm['DL_binary_sha256'], 'Native normalized DL matrix changed')
            input_log = logs / (scientific['name'] + '_RNA_input_parity.log')
            execute([gate['runtime']['rscript'], '--vanilla', str(HERE / 'verify_representation_inputs.R'),
                     str(target), str(prep), str(old_prep), scientific['space'], args.mode], environment, input_log,
                    'RNA_ASSAY_AND_GENE_ORDER_EXACT_CONTROL_COORDINATES_EXACT')
            representation_prep = (Path(gate['reference_root']) / 'GBM_representation_controls' / sample / budget / scientific['name'])
            comparisons = [('archive_representation', representation_prep / route)]
            if args.mode == 'default':
                execute([gate['runtime']['rscript'], '--vanilla', str(HERE / 'verify_representation_inputs.R'),
                         str(target), str(prep), str(representation_prep), scientific['space'], 'representation'],
                        environment, logs / (scientific['name'] + '_archived_representation_RNA_parity.log'),
                        'RNA_ASSAY_AND_GENE_ORDER_EXACT_CONTROL_COORDINATES_EXACT')
                comparisons.extend([('old_core', old_prep / route), ('new_core', prep / route)])
            for reference_type, expected_score in comparisons:
                label = scientific['name'] + '_' + reference_type
                diagnostic_roundtrip = (Path(lock['files']['roundtrip_audit']['path']).parent
                                        if args.mode == 'default' and route == 'UMAP2_HDBSCAN_R'
                                        and reference_type in {'old_core', 'new_core'} else None)
                stage = score_parity(verifier, score, expected_score, label, gate['runtime'], environment, logs, checks,
                                     diagnostic_roundtrip=diagnostic_roundtrip)
                stage_checks.append(dict(configuration=scientific['name'], reference_type=reference_type,
                                         RNA_and_gene_order_exact=True, normalized_DL_exact=True,
                                         control_coordinates_exact=True, input_R_log_sha256=sha(input_log), **stage))
                reference_arms = verifier.arm_map(load(expected_score / 'score_manifest.json'))
                for record in [r for r in selected if r['scientific'] == scientific]:
                    eid, _ = reference_arms[(record['library'], record['cutoff'])]
                    result = verifier.compare_terminal(Path(record['terminal_directory']), expected_score / 'terminal' / eid,
                                                       record, record['arm'], eid, checks, fresh=False)
                    destination = (conditions if reference_type == 'archive_representation' else
                                   new_core_conditions if reference_type == 'new_core' else old_core_conditions)
                    destination.append(result)
        expected_conditions = len(configurations) * 2
        need(len(conditions) == expected_conditions, 'Incomplete terminal comparison roster')
        need(all(sha(path) == digest for path, digest in snapshots.items()), 'Immutable fresh source changed')
        check_sources(); worker.validate_gate(args.gate)
        parity = dict(status='passed' if args.mode == 'default' else 'passed_exact', mode=args.mode, sample=sample, budget=budget,
                      terminal_conditions=expected_conditions, configurations=configurations, stage_checks=stage_checks,
                      conditions=conditions, new_core_conditions=new_core_conditions, old_core_conditions=old_core_conditions,
                      artifact_checks=checks, archived_representation_all_artifacts_exact=True,
                      annotation_and_marker_density_exact_to_core=args.mode == 'default',
                      core_membership_diagnostic_classification='independently_reproduced_CSV_roundtrip_difference' if args.mode == 'default' else 'not_applicable',
                      source_lock_sha256=sha(HERE / 'representation_sources.json'), gate_sha256=gate_hash,
                      normalized_RNA_genes_and_DL_unchanged=True, no_old_cache_reuse=True,
                      reference_labels_used_for_fit=False, truth_opened_before_parity=False,
                      job=os.environ['SLURM_JOB_ID'], step=os.environ.get('SLURM_STEP_ID'),
                      completed_at=datetime.now(timezone.utc).isoformat())
        write_new(out / 'parity.json', parity)
        specification = dict(sample=sample, budget=budget, family='representation_control',
                             configuration='default_parity_four_routes' if args.mode == 'default' else configurations[0]['name'],
                             expected_conditions=expected_conditions,
                             conditions=[{k:r[k] for k in ['route', 'library', 'cutoff', 'terminal_directory']} for r in selected],
                             parity_receipt=dict(path=str(out / 'parity.json'), sha256=sha(out / 'parity.json')),
                             source_hashes={r['path']:r['sha256'] for r in lock['files'].values()})
        write_new(out / 'evaluation_spec.json', specification)
        execute([gate['runtime']['python'], '-s', str(HERE / 'evaluate_extension_lfine.py'), '--spec',
                 str(out / 'evaluation_spec.json'), '--out', str(out / 'evaluation')], environment, out / 'evaluation.log')
        verifier.receipt(out / 'evaluation', 'manifest.json', 'COMPLETE')
        evaluation = load(out / 'evaluation/manifest.json')
        need(evaluation['status'] == 'completed' and evaluation['n_conditions'] == expected_conditions and
             evaluation['n_valid'] == evaluation['n_threshold_rows'] == expected_conditions * 2,
             'Incomplete Lfine evaluation')
        need(evaluation['script_sha256'] == lock['files']['extension_evaluator']['sha256'] and
             evaluation['specification_sha256'] == sha(out / 'evaluation_spec.json') and
             evaluation['parity_receipt'] == specification['parity_receipt'] and
             evaluation['outputs']['metrics.csv.gz'] == sha(out / 'evaluation/metrics.csv.gz'),
             'Lfine output receipt differs from verified input')
        check_sources()
        receipt = dict(parity, lfine_valid_threshold_rows=expected_conditions * 2,
                       lfine_manifest_sha256=sha(out / 'evaluation/manifest.json'),
                       lfine_metrics_sha256=sha(out / 'evaluation/metrics.csv.gz'),
                       parity_sha256=sha(out / 'parity.json'), config_sha256=sha(out / 'run_config.json'))
        write_new(out / 'manifest.json', receipt)
        (out / 'COMPLETE').write_text(sha(out / 'manifest.json') + '\n')
        print(json.dumps(dict(status=receipt['status'], mode=args.mode, terminal_conditions=expected_conditions,
                              lfine_threshold_rows=expected_conditions * 2, output=str(out))), flush=True)
    except BaseException:
        write_new(out / 'FAILURE.json', dict(status='failed_preserved', traceback=traceback.format_exc(),
                                            job=os.environ['SLURM_JOB_ID']))
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gate', required=True, type=Path)
    parser.add_argument('--mode', choices=['default', 'changed', 'task'], required=True)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--index', type=int)
    parser.add_argument('--default-proof', type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args(); args.gate = args.gate.resolve()
    run(args)
