"""Prespecified initialization diagnostic: seeds 0..4, plus seed42 parity control.

Keep all historical-route inputs, initial calls, splits, batches and thresholds
fixed. Never select a winning seed or replace archived annotation references.
"""
from pathlib import Path
import hashlib
import json
import os
import random
import sys
import time

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT = ROOT / 'results/hvg_ptc_20260916_v1/ptc_paper_baseline'
DEST = OUT / 'recovery_diagnosis/initialization_replicates'
ROUTES = ['NMT_Thyroid_Seurat_none', 'TTU_Pubmed_UMAPHDBSCAN_mean']
SEEDS = [0, 1, 2, 3, 4, 42]


def general(label):
    if label.startswith('NCOMMREFF+'):
        return label.split('+', 1)[1]
    if label.startswith('cancer+'):
        return '+'.join(label.split('+')[3:])
    return label


def fit():
    import numpy as np
    import pandas as pd
    import torch
    from torch.utils.data import TensorDataset, DataLoader, random_split
    sys.path.insert(0, str(ROOT / 'handoff/ptc_recovery_20260916'))
    from refine import build_model

    start = time.perf_counter()
    task = int(os.environ['SLURM_ARRAY_TASK_ID'])
    route, seed = ROUTES[task // len(SEEDS)], SEEDS[task % len(SEEDS)]
    group = 'NMT' if route.startswith('NMT') else 'TTU'
    dest = DEST / route / f'seed_{seed}'
    dest.mkdir(parents=True, exist_ok=False)
    ref = OUT / 'replay_selected_routes_full_parallel' / route
    cells = pd.read_csv(ref / 'cells.csv', dtype=str, keep_default_na=False)
    z = np.load(ref / 'terminal_DL/terminal.npz', allow_pickle=False)
    initial = z['initial']
    known = np.flatnonzero(initial != 'Undecided')
    pool = z['pool_indices']
    classes = z['classes']
    x = np.memmap(ref / 'DL_archived_CCA2000.float32.bin', mode='r', dtype='<f4', shape=(len(cells), 2000))
    torch.set_num_threads(4)
    torch.set_num_interop_threads(1)
    torch.use_deterministic_algorithms(True)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    lookup = {label:i for i,label in enumerate(classes)}
    dataset = TensorDataset(torch.from_numpy(np.array(x[known], dtype=np.float32, copy=True)),
                            torch.tensor([lookup[v] for v in initial[known]], dtype=torch.long))
    ntrain = round(len(known) * .90)
    train, val = random_split(dataset, [ntrain, len(known)-ntrain], generator=torch.Generator().manual_seed(42))
    assert np.array_equal(known[np.asarray(train.indices)], z['train_indices'])
    assert np.array_equal(known[np.asarray(val.indices)], z['validation_indices'])
    loader = DataLoader(train, batch_size=256, shuffle=False, num_workers=0)
    model = build_model(2000, len(classes))
    optimizer = torch.optim.Adamax(model.parameters(), lr=.001)
    criterion = torch.nn.CrossEntropyLoss()
    history = []
    for epoch in range(10):
        total = 0
        correct = 0
        loss_sum = 0.
        for samples, targets in loader:
            predictions = model(samples)
            loss = criterion(predictions, targets)
            assert torch.isfinite(loss)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            total += len(targets)
            correct += int(predictions.argmax(1).eq(targets).sum())
            loss_sum += float(loss.detach()) * len(targets)
        history.append(dict(epoch=epoch+1, train_accuracy=correct/total, train_loss=loss_sum/total))
    with torch.no_grad():
        val_correct = 0
        for samples, targets in DataLoader(val, batch_size=256, shuffle=False, num_workers=0):
            val_correct += int(model(samples).argmax(1).eq(targets).sum())
        probs = np.concatenate([model(torch.from_numpy(np.array(x[pool[lo:lo+256]], dtype=np.float32, copy=True))).numpy()
                                for lo in range(0, len(pool), 256)])
    rounded = np.asarray([round(value, 4) for value in probs.max(1)], dtype=np.float32)
    final = initial.copy()
    final[pool] = np.where(rounded >= .9, classes[probs.argmax(1)], 'Unknown')
    assert np.array_equal(final[known], initial[known])
    torch.save(model.state_dict(), dest / 'model_state.pt')
    np.savez_compressed(dest / 'terminal.npz', final090=final, pool_indices=pool, probabilities=probs, classes=classes)
    pd.DataFrame(history).to_csv(dest / 'training_history.csv', index=False)
    paper = pd.read_csv(OUT / 'paper_baseline_reference.csv.gz', dtype=str, keep_default_na=False).set_index('cell_id')
    selected = cells.group.eq(group).to_numpy()
    old = paper.loc[cells.cell_id, 'paper_final_native'].to_numpy()
    unknown = ['Unknown', 'Undecided', 'No_Annotation']
    old_unknown = np.isin(old, unknown)
    new_unknown = np.isin(final, unknown)
    mismatch = old != final
    record = dict(job=os.environ['SLURM_JOB_ID'], route=route, group=group, model_seed=seed,
                  diagnostic_seeds=[0,1,2,3,4], parity_control_seed=42, split_seed=42, threshold=.9,
                  n_selected_cells=int(selected.sum()), n_mismatches=int((mismatch & selected).sum()),
                  n_unknown=int((new_unknown & selected).sum()),
                  old_called_to_Unknown=int((~old_unknown & new_unknown & selected).sum()),
                  old_Unknown_to_called=int((old_unknown & ~new_unknown & selected).sum()),
                  changed_called_type=int((mismatch & ~old_unknown & ~new_unknown & selected).sum()),
                  known_seed_validation_accuracy=val_correct/len(val),
                  input_manifest_sha256=hashlib.sha256((ref / 'score_manifest.json').read_bytes()).hexdigest(),
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  no_reference_labels_changed=True, no_seed_selection=True)
    if seed == 42:
        assert np.array_equal(final, z['final090'])
        assert np.array_equal(probs, z['probabilities'])
        saved = torch.load(ref / 'terminal_DL/model_state.pt', weights_only=True)
        assert all(torch.equal(value, saved[key]) for key,value in model.state_dict().items())
        record['seed42_control_weights_probabilities_terminal_exact'] = True
    record['elapsed_seconds'] = time.perf_counter() - start
    (dest / 'manifest.json').write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2), flush=True)


def summarize():
    import numpy as np
    import pandas as pd
    from sklearn.metrics import accuracy_score, f1_score, roc_auc_score
    cells = pd.read_csv(OUT / 'replay_selected_routes_full_parallel' / ROUTES[0] / 'cells.csv', dtype=str)
    yframe = pd.read_csv(OUT / 'original_DG_binary_pairs_for_R.csv.gz').set_index('cell_id').loc[cells.cell_id]
    y = yframe.truth.to_numpy()
    old_binary = yframe.prediction.to_numpy()
    paper = pd.read_csv(OUT / 'paper_baseline_reference.csv.gz', dtype=str, keep_default_na=False).set_index('cell_id').loc[cells.cell_id]
    old = paper.paper_final_native.to_numpy()
    t_names = set(json.loads((OUT / 'historical_T_names_from_vignette.json').read_text()))
    manifests = []
    metrics = []
    combined = {}
    for seed in SEEDS:
        labels = np.full(len(cells), '', dtype=old.dtype)
        for route in ROUTES:
            path = DEST / route / f'seed_{seed}'
            record = json.loads((path / 'manifest.json').read_text())
            manifests.append(record)
            selected = cells.group.eq(record['group']).to_numpy()
            z = np.load(path / 'terminal.npz', allow_pickle=False)
            labels[selected] = z['final090'][selected]
        assert (labels != '').all()
        combined[seed] = labels
        binary = np.asarray([int(general(value) in t_names) for value in labels])
        scopes = [('Overall', np.ones(len(cells), dtype=bool))]
        scopes += [(scope, yframe.scope.eq(scope).to_numpy()) for scope in sorted(yframe.scope.unique())]
        for scope,mask in scopes:
            metrics.append(dict(model_seed=seed, scope=scope, F1_class0=f1_score(y[mask], binary[mask], pos_label=0),
                                AUC_binary=roc_auc_score(y[mask], binary[mask]), accuracy=accuracy_score(y[mask], binary[mask]),
                                historical_F1_class0=f1_score(y[mask], old_binary[mask], pos_label=0),
                                historical_AUC_binary=roc_auc_score(y[mask], old_binary[mask]),
                                n_terminal_mismatches=int((labels[mask] != old[mask]).sum())))
    pd.DataFrame(manifests).to_csv(DEST / 'route_summary.csv', index=False)
    pd.DataFrame(metrics).to_csv(DEST / 'metrics_all_seeds.csv', index=False)
    stack = np.stack([combined[seed] for seed in SEEDS if seed != 42])
    audit = cells.copy()
    audit['historical_native'] = old
    audit['n_of_5_seeds_matching_historical'] = (stack == old).sum(0)
    audit['all_5_seeds_agree'] = (stack == stack[0]).all(0)
    audit.to_csv(DEST / 'per_cell_stability.csv.gz', index=False)
    overall = pd.DataFrame(metrics).query("scope == 'Overall' and model_seed != 42")
    ranges = {key:dict(min=float(overall[key].min()), max=float(overall[key].max())) for key in ['F1_class0','AUC_binary','accuracy']}
    counts = {group:dict(n_cells=int(mask.sum()), n_all_5_seeds_agree=int(audit.loc[mask,'all_5_seeds_agree'].sum()),
                        n_no_seed_matches_historical=int(audit.loc[mask,'n_of_5_seeds_matching_historical'].eq(0).sum()))
              for group in ['NMT','TTU'] for mask in [cells.group.eq(group)]}
    report = dict(job=os.environ['SLURM_JOB_ID'], diagnostic_seeds=[0,1,2,3,4], parity_control_seed=42,
                  route_fits_completed=len(manifests), ranges=ranges, per_group_stability=counts,
                  no_seed_selected=True, no_original_labels_modified=True,
                  interpretation='Observed five-seed ranges are diagnostics, not confidence intervals or proof of historical cause.')
    (DEST / 'SUMMARY.json').write_text(json.dumps(report, indent=2) + '\n')
    print(pd.DataFrame(manifests)[['route','model_seed','n_mismatches','n_unknown']].to_string(index=False), flush=True)
    print(overall.to_string(index=False), flush=True)
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    assert os.environ.get('SLURM_JOB_ID'), 'Scientific diagnostics require SLURM'
    summarize() if len(sys.argv)>1 and sys.argv[1]=='summarize' else fit()
