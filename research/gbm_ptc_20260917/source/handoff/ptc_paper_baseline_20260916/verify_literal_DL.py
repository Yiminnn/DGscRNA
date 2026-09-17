"""Call the archived run_dgscrna end to end on the same fixed R-scored input.

This checks a remaining implementation concern in the historical baseline. The
archived function, training loop, eight-worker data loaders and scalar cutoff
comparison are unmodified. The constructor wrapper only retains its model object
for saving/comparison. No hyperparameter or label search is performed.
"""
from pathlib import Path
import hashlib
import importlib.util
import json
import os
import random

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT = ROOT / 'results/hvg_ptc_20260916_v1/ptc_paper_baseline'


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def run():
    assert os.environ.get('SLURM_JOB_ID')
    import anndata as ad
    import numpy as np
    import pandas as pd
    import torch
    route = ['NMT_Thyroid_Seurat_none', 'TTU_Pubmed_UMAPHDBSCAN_mean'][int(os.environ['SLURM_ARRAY_TASK_ID'])]
    ref = OUT / 'replay_selected_routes_marker_union' / route
    dest = OUT / 'literal_original_DL' / route
    assert not dest.exists() or not any(dest.iterdir()), 'Preserve existing attempt outputs'
    dest.mkdir(parents=True, exist_ok=True)
    source = ROOT / 'results/hvg_ptc_20260916_v1/ptc_recovery/archive/tcr/ptc_val/scripts/DGscRNA-Share/R/source.py'
    source_hash = sha(source)
    spec = importlib.util.spec_from_file_location('archived_dgscrna_source', source)
    original = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(original)
    torch.set_num_threads(4)
    torch.set_num_interop_threads(1)
    initial = pd.read_csv(ref / 'initial_calls.csv', dtype=str, keep_default_na=False)
    cells = pd.read_csv(ref / 'cells.csv', dtype=str, keep_default_na=False)
    assert np.array_equal(initial.cell_id, cells.cell_id)
    x = np.memmap(ref / 'DL_archived_CCA2000.float32.bin', mode='r', dtype='<f4', shape=(len(cells), 2000))
    column = 'seurat_clusters_CellMarker_Thyroid_none' if route.startswith('NMT') else 'hdbscan.UMAP_clusters_NCOMMREFF_mean'
    obs = pd.DataFrame({column: initial.initial.to_numpy()}, index=cells.cell_id.to_numpy())
    input_path = dest / 'fixed_original_R_scores.h5ad'
    # Loom supplied string annotations. Prevent modern AnnData from silently
    # converting the fixture to categorical values that cannot accept Unknown.
    ad.AnnData(X=np.asarray(x), obs=obs).write_h5ad(input_path, convert_strings_to_categoricals=False)
    input_hash = sha(input_path)
    held = []
    constructor = original.DeepModel

    def retain_model(*args, **kwargs):
        model = constructor(*args, **kwargs)
        held.append(model)
        return model

    original.DeepModel = retain_model
    random.seed(42)
    np.random.seed(42)
    torch.manual_seed(42)
    torch.use_deterministic_algorithms(True)
    os.chdir(dest)
    original.run_dgscrna(str(dest), input_path.name, study_sets=True, EPOCHS=10)
    assert sha(source) == source_hash and sha(input_path) == input_hash
    assert len(held) == 1
    model = held[0]
    torch.save(model.state_dict(), dest / 'literal_model_state.pt')
    saved = pd.read_csv(dest / ('annotated_' + input_path.name + '.csv'), keep_default_na=False)
    native = saved[column + '_DGscRNA'].to_numpy(dtype=str)
    known = initial.initial.ne('Undecided').to_numpy()
    assert np.array_equal(native[known], initial.initial.to_numpy()[known])
    assert len(native) == len(cells)
    z = np.load(ref / 'terminal_DL/terminal.npz', allow_pickle=False)
    with torch.no_grad():
        probs = np.concatenate([model(torch.from_numpy(np.array(x[part], dtype=np.float32, copy=True))).numpy()
                                for part in np.array_split(z['pool_indices'], np.arange(256, len(z['pool_indices']), 256))])
    np.savez_compressed(dest / 'literal_probabilities.npz', probabilities=probs,
                        pool_indices=z['pool_indices'], classes=z['classes'])
    baseline = torch.load(ref / 'terminal_DL/model_state.pt', weights_only=True)
    delta = max(float(torch.max(torch.abs(value - baseline[key]))) for key, value in model.state_dict().items())
    compare = cells.copy()
    compare['literal_original_terminal'] = native
    compare['refinement_wrapper_terminal'] = z['final090']
    compare['exact'] = native == z['final090']
    compare.to_csv(dest / 'literal_vs_wrapper.csv.gz', index=False)
    compare[~compare.exact].to_csv(dest / 'literal_vs_wrapper_differences.csv', index=False)
    paper = pd.read_csv(OUT / 'paper_baseline_reference.csv.gz', keep_default_na=False).set_index('cell_id')
    selected = cells.group.eq('NMT' if route.startswith('NMT') else 'TTU').to_numpy()
    historical = paper.loc[cells.cell_id, 'paper_final_native'].to_numpy()
    manifest = dict(job=os.environ['SLURM_JOB_ID'], route=route,
        source=str(source), source_sha256=source_hash, input_sha256=input_hash,
        original_function_body_modified=False, constructor_wrapper='retain unmodified model for saving only',
        model_seed=42, historical_seed_unknown=True, data_loader_workers=8,
        n_cells=len(cells), n_terminal_equal_wrapper=int(compare.exact.sum()),
        n_terminal_different_wrapper=int((~compare.exact).sum()),
        model_max_absolute_weight_delta=delta,
        pool_max_absolute_probability_delta=float(np.max(np.abs(probs - z['probabilities']))),
        n_selected_group=int(selected.sum()),
        n_selected_historical_exact=int((native[selected] == historical[selected]).sum()),
        n_selected_historical_mismatch=int((native[selected] != historical[selected]).sum()),
        script_sha256=sha(Path(__file__)))
    (dest / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == '__main__':
    run()
