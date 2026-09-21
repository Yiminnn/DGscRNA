"""Fresh 30-epoch Reviewer B trajectories with historical validation instrumentation.

Training/validation targets are marker pseudo-labels. Lfine truth is opened only
by the separate evaluator after all checkpoint predictions are frozen.
"""
from pathlib import Path
import argparse
import importlib.util
import os
import sys

import support as shared


def arrays_equal(actual, expected, *, fields=None):
    import numpy as np
    with np.load(actual, allow_pickle=False) as a, np.load(expected, allow_pickle=False) as b:
        if fields is None:
            shared.need(set(a.files) == set(b.files), 'Checkpoint fields changed')
            fields = a.files
        for key in fields:
            shared.need(a[key].dtype == b[key].dtype, f'Checkpoint dtype changed: {key}')
            if np.issubdtype(a[key].dtype, np.inexact):
                np.testing.assert_allclose(a[key], b[key], rtol=0, atol=0, equal_nan=True, err_msg=key)
            else:
                np.testing.assert_array_equal(a[key], b[key], err_msg=key)


def weights_equal(actual, expected):
    import torch
    a = torch.load(actual, map_location='cpu', weights_only=True)
    b = torch.load(expected, map_location='cpu', weights_only=True)
    shared.need(set(a) == set(b), 'Weight names changed')
    for name in a:
        shared.need(a[name].dtype == b[name].dtype and torch.equal(a[name], b[name]), f'Weight differs: {name}')


def wrap_terminal(directory, cells, training, arm, score, epoch):
    """Expose unchanged learning checkpoint arrays through the shared evaluation IO."""
    import numpy as np
    with np.load(directory / 'terminal.npz', allow_pickle=False) as arrays:
        initial = arrays['initial']
        pool = arrays['pool_indices']
        known = np.flatnonzero(initial != 'Undecided')
        shared.need(np.array_equal(pool, np.flatnonzero(initial == 'Undecided')), 'Checkpoint unresolved pool changed')
        predictions = cells[['cell_id']].copy()
        for stage in ['initial', 'final090', 'final070']:
            predictions[stage] = arrays[stage]
        for stage in ['final090', 'final070']:
            shared.need(np.array_equal(arrays[stage][known], initial[known]), 'Known pseudo-labels changed')
        if training['training_executed']:
            probabilities = arrays['probabilities']
            shared.need(probabilities.shape == (len(pool), len(arrays['classes'])) and np.isfinite(probabilities).all(),
                        'Invalid saved checkpoint probabilities')
            rounded = np.asarray([round(value, 4) for value in probabilities.max(axis=1)], dtype=np.float32)
            calls = arrays['classes'][probabilities.argmax(axis=1)]
            for stage, threshold in [('final090', .9), ('final070', .7)]:
                shared.need(np.array_equal(arrays[stage][pool], np.where(rounded >= threshold, calls, 'Unknown')),
                            'Checkpoint threshold reconstruction failed')
        n_called = {stage: int((~np.isin(arrays[stage], ['Unknown', 'Undecided'])).sum())
                    for stage in ['final090', 'final070']}
    predictions.to_csv(directory / 'predictions.csv.gz', index=False)
    manifest = dict(status='completed', terminal_valid=True, family='learning_control', arm=arm,
                    epoch=epoch, dl_status=training['dl_status'], training_executed=training['training_executed'],
                    n_known=training['n_known'], n_pool=training['n_pool'], n_training_classes=training['n_training_classes'],
                    n_final_called090=n_called['final090'], n_final_called070=n_called['final070'],
                    DL_sha256=score['DL_binary_sha256'], DL_features=score['DL_features'],
                    terminal_sha256=shared.sha(directory / 'terminal.npz'),
                    predictions_sha256=shared.sha(directory / 'predictions.csv.gz'),
                    known_labels_unchanged=True, threshold_reconstructed=True, reference_labels_used_for_fit=False,
                    identical_result_reused=False, validation_target='held-out marker pseudo-labels, not independent biological truth',
                    job=os.environ['SLURM_JOB_ID'])
    shared.write_new(directory / 'terminal_manifest.json', manifest)
    (directory / 'TERMINAL_COMPLETE').write_text(shared.sha(directory / 'terminal_manifest.json') + '\n')
    return manifest


def run_learning(gate, cfg, library, pilot_root=None):
    import numpy as np
    import pandas as pd
    import torch
    shared.validate_sources()
    source_hashes = {str(p): shared.sha(p) for p in
                     [Path(__file__), Path(shared.__file__), shared.HERE / 'SOURCE_DERIVATION.json',
                      shared.HERE / 'refine_learning.py']}
    root = Path(gate['output_root']).resolve()
    base = root / ('pilots/reviewer_b_default_TKU4163_hvg2000' if pilot_root else 'extensions/reviewer_b')
    dest = (base / 'learning' / cfg['sample'] / cfg['budget'] / f"seed{cfg['model_seed']}" / library).resolve()
    shared.need(dest.is_relative_to(root) and not dest.is_relative_to(root / 'core') and not dest.exists(),
                f'Preserve existing/unsafe learning destination: {dest}')
    source_root, prepared, provenance = shared.source(gate, cfg, pilot_root)
    before = shared.geometry.file_inventory(source_root)
    source = prepared / 'UMAP2_HDBSCAN_R'
    score = shared.checked(source, 'score_manifest.json', 'SCORE_COMPLETE')
    aid = shared.geometry.arm_map(score)[(library, 'mean')]
    arm = score['arms'][aid]
    reference = (shared.OLD / ('learning_parity' if pilot_root else 'learning') / cfg['sample'] / cfg['budget'] /
                 f"seed{cfg['model_seed']}" / library)
    cells = pd.read_csv(source / 'cells.csv', dtype=str, keep_default_na=False)
    initial = pd.read_csv(source / 'initial_calls.csv.gz', dtype=str, keep_default_na=False)
    shared.need(np.array_equal(cells.cell_id, initial.cell_id), 'Seed/cell order mismatch')
    binary = Path(score['DL_binary'])
    shared.need(shared.sha(binary) == score['DL_binary_sha256'], 'Normalized DL expression changed')
    x = np.memmap(binary, mode='r', dtype='<f4', shape=(len(cells), int(score['DL_features'])))
    shared.need(np.isfinite(x).all(), 'Nonfinite normalized DL expression')
    dest.mkdir(parents=True)
    environment = shared.environment(gate, dest)
    # Only this direct training process needs the common output context. No old
    # module/cache is imported, and training always creates a new destination.
    os.environ['DGSCRNA_REFERENCE_OUT'] = str(dest)
    os.environ['DGSCRNA_REQUIRE_SLURM'] = '1'
    shared.write_new(dest / 'protocol.json', dict(config=cfg, library=library, primary=library == shared.LIBRARIES[0],
                     source=provenance, gate_sha256=gate['_gate_sha256'], seed_column=aid,
                     DL_sha256=score['DL_binary_sha256'], source_score_sha256=shared.sha(source / 'score_manifest.json'),
                     source_snapshot_sha256=shared.sha(shared.HERE / 'SOURCE_DERIVATION.json'),
                     runner_sha256=shared.sha(__file__), support_sha256=shared.sha(shared.__file__),
                     source_hashes=source_hashes,
                     no_old_training_cache=True, old_L1_evaluation_used=False,
                     validation_target='held-out marker pseudo-labels, not independent biological truth'))
    try:
        torch.set_num_threads(4)
        spec = importlib.util.spec_from_file_location('fresh_reviewer_b_learning', shared.HERE / 'refine_learning.py')
        refine = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(refine)
        refine.PARAMS.update(epochs=cfg['epochs'], model_seed=cfg['model_seed'], split_seed=42,
                             input='Native R normalized selected RNA; original fixed genes/order')
        refine.train_cache(x, initial[aid].to_numpy(dtype=str), dest,
                           dict(source_score=str(source), score_sha256=shared.sha(source / 'score_manifest.json'),
                                DL_sha256=shared.sha(binary), library=library, cutoff='mean', seed_column=aid,
                                reference_labels_used_for_fit=False, split_seed=42, source_sha256=shared.sha(refine.__file__)))
        training = shared.checked(dest, 'training_manifest.json', 'COMPLETE')
        shared.need(training['terminal_valid'] is True, 'Invalid training state')
        old_training = shared.checked(reference, 'training_manifest.json', 'COMPLETE')
        for key in ['dl_status', 'training_executed', 'prediction_executed', 'n_cells', 'n_known', 'n_pool',
                    'n_training_classes', 'classes', 'input_width', 'n_train', 'n_validation', 'known_seed_validation_accuracy']:
            shared.need(training.get(key) == old_training.get(key), f'Learning metadata differs: {key}')
        arrays_equal(dest / 'terminal.npz', reference / 'terminal.npz')
        if training['training_executed']:
            shared.need(shared.load(dest / 'training_history.json') == shared.load(reference / 'training_history.json'),
                        'Per-epoch training/held-out validation history differs')
            weights_equal(dest / 'model_state.pt', reference / 'model_state.pt')
            epochs = [e for e in (5, 10, 20, 30) if e <= cfg['epochs']]
            pd.DataFrame(shared.load(dest / 'training_history.json')).to_csv(dest / 'learning_history.csv', index=False)
        else:
            shared.need(training['dl_status'] in {'no_op_all_initially_known', 'no_known_labels_archived_Undecided_terminal',
                                                 'structural_insufficient_known_split'}, 'Unexpected no-training state')
            epochs = [0]
        evaluations, conditions = [], []
        for epoch in epochs:
            actual = dest / f'epoch{epoch:02d}' if epoch else dest
            expected = reference / f'epoch{epoch:02d}' if epoch else reference
            arrays_equal(actual / 'terminal.npz', expected / 'terminal.npz')
            if epoch:
                weights_equal(actual / 'model_state.pt', expected / 'model_state.pt')
            terminal = wrap_terminal(actual, cells, training, arm, score, epoch)
            # Epoch10 must also equal the independently archived same-seed
            # original-width 10-epoch model, not merely the extended trajectory.
            if epoch == 10:
                control = 'original' if cfg['model_seed'] == 42 else f"model_seed{cfg['model_seed']}"
                original_route = (Path(gate['reference_root']) / 'GBM_DL_controls' / cfg['sample'] / cfg['budget'] /
                                  control / 'UMAP2_HDBSCAN_R')
                original_score = shared.checked(original_route, 'score_manifest.json', 'SCORE_COMPLETE')
                original_aid = shared.geometry.arm_map(original_score)[(library, 'mean')]
                original10 = original_route / 'terminal' / original_aid
                shared.checked(original10, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
                arrays_equal(actual / 'terminal.npz', original10 / 'terminal.npz',
                             fields=['initial', 'final090', 'final070', 'classes', 'probabilities', 'train_indices', 'validation_indices'])
                weights_equal(actual / 'model_state.pt', original10 / 'model_state.pt')
            condition = dict(route='UMAP2_HDBSCAN_R', library=library, cutoff='mean', status='passed_exact',
                             epoch=epoch, actual_manifest_sha256=shared.sha(actual / 'terminal_manifest.json'),
                             reference_terminal_sha256=shared.sha(expected / 'terminal.npz'),
                             arrays_exact=True, weights_exact=bool(epoch), histories_exact=training['training_executed'],
                             same_seed_original10_exact=epoch == 10)
            proof = dict(status='passed_exact', sample=cfg['sample'], budget=cfg['budget'], family='learning_control',
                         model_seed=cfg['model_seed'], epoch=epoch, conditions=[condition], terminal_conditions=1,
                         truth_used_for_fitting=False, heldout_targets='marker pseudo-labels',
                         gate_sha256=gate['_gate_sha256'], job=os.environ['SLURM_JOB_ID'])
            shared.write_new(actual / 'checkpoint_parity.json', proof)
            evaluation_cfg = dict(cfg, name=f"seed{cfg['model_seed']}_epoch{epoch:02d}")
            evaluation_conditions = [dict(route='UMAP2_HDBSCAN_R', library=library, cutoff='mean', terminal_directory=str(actual))]
            evaluation = shared.evaluate(gate, actual, evaluation_cfg, 'learning_control', evaluation_conditions,
                                         actual / 'checkpoint_parity.json', environment)
            evaluations.append(dict(epoch=epoch, **evaluation))
            conditions.append(condition)
        shared.need(shared.geometry.file_inventory(source_root) == before, 'Accepted core changed')
        shared.need(shared.sha(binary) == score['DL_binary_sha256'], 'DL input changed during training')
        shared.validate_sources()
        shared.need(all(shared.sha(p) == digest for p, digest in source_hashes.items()),
                    'Learning adapter source changed during execution')
        accepted = dict(status='passed_exact', sample=cfg['sample'], budget=cfg['budget'], model_seed=cfg['model_seed'],
                        epochs=cfg['epochs'], library=library, primary=library == shared.LIBRARIES[0],
                        training_executed=training['training_executed'], dl_status=training['dl_status'],
                        conditions=conditions, Lfine=evaluations, core_unchanged=True,
                        source_hashes=source_hashes,
                        training_manifest_sha256=shared.sha(dest / 'training_manifest.json'),
                        validation_target='held-out marker pseudo-labels, not independent biological truth',
                        gate_sha256=gate['_gate_sha256'], job=os.environ['SLURM_JOB_ID'])
        shared.write_new(dest / 'acceptance.json', accepted)
        (dest / 'LEARNING_VERIFIED_COMPLETE').write_text(shared.sha(dest / 'acceptance.json') + '\n')
        print('LEARNING_VERIFIED_COMPLETE', cfg['sample'], cfg['budget'], cfg['model_seed'], library, flush=True)
        return dest, training
    except BaseException as error:
        shared.write_new(dest / 'failure_preserved.json', dict(status='failed_preserved', config=cfg, library=library,
                          error=repr(error), gate_sha256=gate['_gate_sha256'], job=os.environ['SLURM_JOB_ID']))
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gate', required=True, type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--pilot-full-run', type=Path)
    mode.add_argument('--index', type=int, help='Index in the 18 primary learning-trajectory roster')
    args = parser.parse_args()
    gate, _ = shared.context(args.gate.resolve())
    configs = [c for c in shared.tasks() if c['task'] == 'learning']
    shared.need(len(configs) == 18, 'Learning roster changed')
    if args.pilot_full_run:
        cfg = dict(task='learning', sample='TKU4163', budget='hvg2000', model_seed=42, epochs=10)
        dest, _ = run_learning(gate, cfg, shared.LIBRARIES[0], args.pilot_full_run.resolve())
        pilot = Path(gate['output_root']) / 'pilots/reviewer_b_default_TKU4163_hvg2000'
        shared.write_new(pilot / 'learning_parity.json', dict(status='passed_exact', epochs=10,
                         acceptance=str(dest / 'acceptance.json'), acceptance_sha256=shared.sha(dest / 'acceptance.json'),
                         source_snapshot_sha256=shared.sha(shared.HERE / 'SOURCE_DERIVATION.json')))
        (pilot / 'LEARNING_DEFAULT_PARITY_COMPLETE').write_text(shared.sha(pilot / 'learning_parity.json') + '\n')
    else:
        shared.need(0 <= args.index < 18, 'Invalid learning trajectory index')
        pilot = Path(gate['output_root']) / 'pilots/reviewer_b_default_TKU4163_hvg2000'
        shared.checked(pilot, 'learning_parity.json', 'LEARNING_DEFAULT_PARITY_COMPLETE')
        cfg = configs[args.index]
        _, training = run_learning(gate, cfg, shared.LIBRARIES[0])
        if not training['training_executed']:
            run_learning(gate, cfg, shared.LIBRARIES[1])


if __name__ == '__main__':
    main()
