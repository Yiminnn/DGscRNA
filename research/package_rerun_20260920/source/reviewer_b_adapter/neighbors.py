"""Fresh Reviewer B neighbor conditions using frozen R bodies and package DL."""
from pathlib import Path
import argparse
import os
import shutil
import sys

import support as shared


def run_neighbor(gate, cfg, pilot_root=None):
    root = Path(gate['output_root']).resolve()
    base = root / ('pilots/reviewer_b_default_TKU4163_hvg2000' if pilot_root else 'extensions/reviewer_b')
    dest = (base / 'neighbors' / cfg['sample'] / cfg['budget'] / cfg['name']).resolve()
    shared.need(dest.is_relative_to(root) and not dest.is_relative_to(root / 'core') and not dest.exists(),
                f'Preserve existing/unsafe neighbor destination: {dest}')
    source_root, source, provenance = shared.source(gate, cfg, pilot_root)
    before = shared.geometry.file_inventory(source_root)
    route = cfg['space'] + '_' + cfg['method']
    default = cfg['snn_k'] == 20 and cfg['umap_neighbors'] == 30
    reference = (Path(gate['reference_root']) / 'GBM' / cfg['sample'] / cfg['budget'] / route if default else
                 shared.OLD / 'neighbors' / cfg['sample'] / cfg['budget'] / cfg['name'] / route)
    old_score = shared.checked(reference, 'score_manifest.json', 'SCORE_COMPLETE')
    dest.mkdir(parents=True)
    env = shared.environment(gate, dest)
    env['DGSCRNA_ONLY_ROUTE'] = route
    try:
        (dest / 'markers').mkdir()
        shutil.copyfile(gate['markers']['path'], dest / 'markers/libraries.json')
        prepared = dest / 'prepare'
        config = dict(cfg, dest=str(prepared), source_prepare=str(source))
        shared.write_new(dest / 'config.json', config)
        shared.write_new(dest / 'input_provenance.json', dict(gate_sha256=gate['_gate_sha256'], source=provenance,
                    source_snapshot_sha256=shared.sha(shared.HERE / 'SOURCE_DERIVATION.json'),
                    runner_sha256=shared.sha(__file__), support_sha256=shared.sha(shared.__file__),
                    default_recomputed=default, all_conditions_recomputed=True,
                    no_old_training_cache=True, old_L1_evaluation_used=False))
        shared.run([gate['runtime']['rscript'], '--vanilla', str(shared.HERE / 'prepare_neighbors.R'), str(dest / 'config.json')],
                   dest / 'prepare.log', env)
        shared.run([gate['runtime']['rscript'], '--vanilla', str(shared.HERE / 'score_neighbors.R'), cfg['sample'],
                    str(prepared), str(dest / 'config.json')], dest / 'score.log', env)
        actual = prepared / route
        score = shared.checked(actual, 'score_manifest.json', 'SCORE_COMPLETE')
        aarms, earms = shared.geometry.arm_map(score), shared.geometry.arm_map(old_score)
        shared.need(set(aarms) == set(earms) and len(aarms) == 48, 'Full 48-arm marker scoring roster changed')
        verifier = shared.geometry.load_verifier(gate)
        checks = []
        for filename in ['cells.csv', 'clusters.csv', 'initial_calls.csv.gz']:
            verifier.frame_equal(verifier.frame(actual / filename), verifier.frame(reference / filename), filename, checks)
        libraries = list(shared.load(dest / 'markers/libraries.json'))
        shared.write_new(dest / 'verify_libraries.json', libraries)
        shared.run([gate['runtime']['rscript'], '--vanilla',
                    gate['launcher_files']['verify_packaged_r_artifacts.R']['path'],
                    str(actual), str(reference), str(dest / 'verify_libraries.json')], dest / 'R_object_parity.log', env)
        shared.need('R_DEG_AND_NAMED_DENSITY_OBJECTS_EXACT' in (dest / 'R_object_parity.log').read_text(),
                    'R DEG/density object verifier did not confirm parity')
        conditions, evaluation_conditions = [], []
        for library in shared.LIBRARIES:
            key = (library, 'mean')
            aid, eid = aarms[key], earms[key]
            shared.run([gate['runtime']['python'], '-s', '-m', 'dgscrna.reference.backend.terminal',
                        str(actual), aid], dest / (aid + '.log'), env)
            td = actual / 'terminal' / aid
            terminal = shared.checked(td, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
            shared.need(Path(terminal['cache_directory']).resolve().is_relative_to(dest / 'DL_cache'), 'External DL cache reused')
            training = shared.load(td / 'training_manifest.json')
            shared.need(Path(training['provenance']['first_condition']).resolve().is_relative_to(dest), 'Cached fit originated elsewhere')
            record = dict(route=route, library=library, cutoff='mean', dl_status=terminal['dl_status'],
                          training_executed=terminal['training_executed'])
            conditions.append(verifier.compare_terminal(td, reference / 'terminal' / eid, record, aid, eid, checks, fresh=False))
            evaluation_conditions.append(dict(route=route, library=library, cutoff='mean', terminal_directory=str(td)))
        shared.need(shared.geometry.file_inventory(source_root) == before, 'Accepted core was modified')
        proof = dict(status='passed_exact', sample=cfg['sample'], budget=cfg['budget'], family='neighbor_control',
                     configuration=cfg['name'], config=cfg, terminal_conditions=2, conditions=conditions, checks=checks,
                     all48_initial_arms_exact=True, DEG_and_density_objects_exact=True,
                     core_unchanged=True, fresh_unit_cache=True, gate_sha256=gate['_gate_sha256'],
                     input_provenance_sha256=shared.sha(dest / 'input_provenance.json'), job=os.environ['SLURM_JOB_ID'])
        shared.write_new(dest / 'parity.json', proof)
        evaluation = shared.evaluate(gate, dest, cfg, 'neighbor_control', evaluation_conditions, dest / 'parity.json', env)
        shared.need(shared.geometry.file_inventory(source_root) == before, 'Core changed during evaluation')
        shared.write_new(dest / 'acceptance.json', dict(status='passed_exact', sample=cfg['sample'], budget=cfg['budget'],
                    configuration=cfg['name'], family='neighbor_control', terminal_conditions=2,
                    parity_sha256=shared.sha(dest / 'parity.json'), Lfine=evaluation,
                    source_snapshot_sha256=shared.sha(shared.HERE / 'SOURCE_DERIVATION.json'),
                    core_unchanged=True, gate_sha256=gate['_gate_sha256'], job=os.environ['SLURM_JOB_ID']))
        (dest / 'NEIGHBOR_VERIFIED_COMPLETE').write_text(shared.sha(dest / 'acceptance.json') + '\n')
        print('NEIGHBOR_VERIFIED_COMPLETE', cfg['sample'], cfg['budget'], cfg['name'], flush=True)
        return dest
    except BaseException as error:
        shared.write_new(dest / 'failure_preserved.json', dict(status='failed_preserved', config=cfg,
                          error=repr(error), gate_sha256=gate['_gate_sha256'], job=os.environ['SLURM_JOB_ID']))
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gate', required=True, type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--pilot-full-run', type=Path)
    mode.add_argument('--index', type=int, help='Index in the 66-condition neighbor roster')
    args = parser.parse_args()
    gate, _ = shared.context(args.gate.resolve())
    configs = [t for t in shared.tasks() if t['task'] == 'neighbors']
    shared.need(len(configs) == 66, 'Neighbor roster changed')
    if args.pilot_full_run:
        selected = [c for c in configs if c['sample'] == 'TKU4163' and c['budget'] == 'hvg2000'
                    and c['snn_k'] == 20 and c['umap_neighbors'] == 30]
        shared.need(len(selected) == 3, 'Default parity requires three distinct default routes')
        outputs = [run_neighbor(gate, c, args.pilot_full_run.resolve()) for c in selected]
        pilot = Path(gate['output_root']) / 'pilots/reviewer_b_default_TKU4163_hvg2000'
        shared.write_new(pilot / 'neighbors_parity.json', dict(status='passed_exact', default_routes=3,
                         initial_arms_per_route=48, terminal_conditions=6, Lfine_threshold_rows=12,
                         receipts={str(p / 'acceptance.json'): shared.sha(p / 'acceptance.json') for p in outputs},
                         source_snapshot_sha256=shared.sha(shared.HERE / 'SOURCE_DERIVATION.json')))
        (pilot / 'NEIGHBOR_DEFAULT_PARITY_COMPLETE').write_text(shared.sha(pilot / 'neighbors_parity.json') + '\n')
    else:
        shared.need(0 <= args.index < 66, 'Invalid neighbor index')
        pilot = Path(gate['output_root']) / 'pilots/reviewer_b_default_TKU4163_hvg2000'
        shared.checked(pilot, 'neighbors_parity.json', 'NEIGHBOR_DEFAULT_PARITY_COMPLETE')
        run_neighbor(gate, configs[args.index])


if __name__ == '__main__':
    main()
