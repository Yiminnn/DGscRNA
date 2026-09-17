"""Audit historical seed columns and existing DL confidence; never refit or relabel."""
from pathlib import Path
import hashlib
import json
import os
import re

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
BASE = ROOT / 'results/hvg_ptc_20260916_v1'
OUT = BASE / 'ptc_paper_baseline'
ARCHIVE = BASE / 'ptc_recovery/archive/tcr'
ROUTES = ['NMT_Thyroid_Seurat_none', 'TTU_Pubmed_UMAPHDBSCAN_mean']


def run():
    assert os.environ.get('SLURM_JOB_ID'), 'Scientific audits require SLURM'
    import h5py
    import numpy as np
    import pandas as pd

    dest = OUT / 'recovery_diagnosis'
    dest.mkdir(exist_ok=True)
    parent = OUT / 'replay_selected_routes_full_parallel'
    cells = pd.read_csv(parent / ROUTES[0] / 'cells.csv', dtype=str, keep_default_na=False)
    assert cells.cell_id.is_unique and cells.original_cell_id.is_unique
    meta = pd.read_csv(ARCHIVE / 'rawdata/metadata.txt', sep='\t', dtype=str)
    samplemap = dict(zip(meta.sc_ID, meta.Sample))
    initial = {}
    for route in ROUTES:
        d = pd.read_csv(parent / route / 'initial_calls.csv', dtype=str, keep_default_na=False)
        assert d.cell_id.equals(cells.cell_id)
        initial[route] = d.initial.to_numpy()

    def route_for(name):
        if re.search(r'^seurat_clusters_(?:CellMarker_)?Thyroid_none(?:$|_)', name):
            return ROUTES[0]
        if re.search(r'^hdbscan\.UMAP_clusters_NCOMMREFF_mean(?:$|_)', name):
            return ROUTES[1]
        return None

    def canonical_ids(frame):
        for key in ['cell_id', 'CellID', 'X', 'Unnamed: 0']:
            if key in frame and frame[key].is_unique:
                if set(frame[key]) == set(cells.cell_id):
                    return frame[key], 'exact canonical cell_id'
                if set(frame[key]) == set(cells.original_cell_id):
                    mapping = dict(zip(cells.original_cell_id, cells.cell_id))
                    return frame[key].map(mapping), 'exact original cell_id'
        for samplecol in ['sample.name', 'sample', 'orig.ident']:
            if samplecol not in frame:
                continue
            sample = frame[samplecol].map(samplemap).fillna(frame[samplecol])
            for key in ['CellID', 'X', 'Unnamed: 0', 'cell_id']:
                if key not in frame:
                    continue
                barcode = frame[key].str.extract(r'([ACGT]{12,})', expand=False)
                ids = sample + '_' + barcode
                if ids.is_unique and set(ids) == set(cells.cell_id):
                    return ids, 'exact sample + barcode'
        return None, 'no verified exact join'

    inventory = []
    comparisons = []

    def audit_frame(path, frame, all_columns):
        candidates = [c for c in all_columns if route_for(c)]
        ids, method = canonical_ids(frame)
        inventory.append(dict(source=str(path.relative_to(ARCHIVE)), n_rows=len(frame),
                              n_columns=len(all_columns), candidate_columns=candidates,
                              identity_columns=[c for c in frame if c not in candidates], join=method))
        if ids is None:
            return
        frame.index = ids
        frame = frame.loc[cells.cell_id]
        for column in candidates:
            # Final/general columns are listed in the inventory only. Comparing
            # their simplified labels with native pre-DL calls is not a seed audit.
            if 'DGCyTOF' in column or 'DGscRNA' in column:
                continue
            route = route_for(column)
            old = frame[column].to_numpy()
            new = initial[route]
            known_old = old != 'Undecided'
            known_new = new != 'Undecided'
            terminal_column = ('DGCyTOF' in column or 'DGscRNA' in column)
            record = dict(source=str(path.relative_to(ARCHIVE)), column=column, route=route,
                          column_role='terminal (cannot establish initial equality)' if terminal_column else 'pre-DL named column',
                          n_cells=len(frame), n_same_as_current_initial=int((old == new).sum()),
                          n_historical_known=int(known_old.sum()), n_current_known=int(known_new.sum()),
                          n_historical_known_current_Undecided=int((known_old & ~known_new).sum()),
                          n_historical_Undecided_current_known=int((~known_old & known_new).sum()),
                          n_conflicting_known_types=int((known_old & known_new & (old != new)).sum()))
            comparisons.append(record)
            if not terminal_column:
                mismatch = old != new
                rows = cells.loc[mismatch].copy()
                rows['archived_column_value'] = old[mismatch]
                rows['current_initial'] = new[mismatch]
                token = hashlib.sha256((str(path) + column).encode()).hexdigest()[:12]
                rows.to_csv(dest / f'initial_differences_{token}.csv.gz', index=False)

    for relative in ['rawdata/data_with_validation.csv', 'rawdata/data_with_validation+3_cell_types 2.csv',
                     'ptc_val/scripts/DGscRNA-Share/annotated_integrated_data.loom.csv',
                     'ptc_val/ptc_batch/annotated_integrated_data.loom.csv',
                     'ptc_val/scripts/DGscRNA-Share/data/ptc1/annotated_integrated_data.loom.csv']:
        path = ARCHIVE / relative
        columns = pd.read_csv(path, nrows=0).columns.tolist()
        selected = [c for c in columns if route_for(c) or c in ['cell_id','CellID','X','Unnamed: 0','sample.name','sample','orig.ident']]
        frame = pd.read_csv(path, usecols=selected, dtype=str, keep_default_na=False)
        audit_frame(path, frame, columns)

    path = ARCHIVE / 'rawdata/integrated_data.loom'
    with h5py.File(path, 'r') as h:
        columns = list(h['col_attrs'])
        selected = [c for c in columns if route_for(c) or c in ['cell_id','CellID','X','sample.name','sample','orig.ident']]
        values = {}
        for c in selected:
            dataset = h['col_attrs'][c]
            values[c] = dataset.asstr()[:] if h5py.check_string_dtype(dataset.dtype) else dataset[:].astype(str)
        audit_frame(path, pd.DataFrame(values), columns)

    pd.DataFrame(comparisons, columns=['source','column','route','column_role','n_cells',
        'n_same_as_current_initial','n_historical_known','n_current_known',
        'n_historical_known_current_Undecided','n_historical_Undecided_current_known',
        'n_conflicting_known_types']).to_csv(dest / 'historical_initial_column_comparison.csv', index=False)
    (dest / 'historical_column_inventory.json').write_text(json.dumps(inventory, indent=2) + '\n')

    reference = pd.read_csv(OUT / 'paper_baseline_reference.csv.gz', dtype=str, keep_default_na=False).set_index('cell_id')
    historical = reference.loc[cells.cell_id, 'paper_final_native'].to_numpy()
    rows = []
    details = []
    for route in ROUTES:
        group = 'NMT' if route.startswith('NMT') else 'TTU'
        z = np.load(parent / route / 'terminal_DL/terminal.npz', allow_pickle=False)
        pool = z['pool_indices']
        selected = cells.group.to_numpy()[pool] == group
        index = pool[selected]
        probs = z['probabilities'][selected]
        sorted_probs = np.sort(probs, axis=1)
        frame = cells.iloc[index].copy()
        frame['historical'] = historical[index]
        frame['new_terminal'] = z['final090'][index]
        frame['new_argmax'] = z['classes'][probs.argmax(axis=1)]
        frame['confidence'] = probs.max(axis=1)
        frame['margin'] = sorted_probs[:, -1] - sorted_probs[:, -2]
        frame['mismatch'] = frame.historical.ne(frame.new_terminal)
        frame['group'] = group
        details.append(frame)
        known_old = ~frame.historical.isin(['Unknown', 'Undecided', 'No_Annotation'])
        new_unknown = frame.new_terminal.isin(['Unknown', 'Undecided', 'No_Annotation'])
        masks = {'all_DL_pool':np.ones(len(frame), dtype=bool),
                 'all_mismatches':frame.mismatch.to_numpy(),
                 'historical_called_to_Unknown':(known_old & new_unknown).to_numpy()}
        for label in sorted(frame.historical.unique()):
            masks['historical_called_to_Unknown:' + label] = (known_old & new_unknown & frame.historical.eq(label)).to_numpy()
        for category, mask in masks.items():
            g = frame.loc[mask]
            if not len(g):
                continue
            rows.append(dict(group=group, category=category, n_cells=len(g),
                             confidence_min=float(g.confidence.min()), confidence_q25=float(g.confidence.quantile(.25)),
                             confidence_median=float(g.confidence.median()), confidence_q75=float(g.confidence.quantile(.75)),
                             confidence_max=float(g.confidence.max()), n_confidence_below_070=int(g.confidence.lt(.7).sum()),
                             n_confidence_085_to_090=int((g.confidence.ge(.85)&g.confidence.lt(.9)).sum()),
                             n_argmax_same_as_historical=int(g.new_argmax.eq(g.historical).sum())))
    pd.DataFrame(rows).to_csv(dest / 'existing_model_confidence_summary.csv', index=False)
    pd.concat(details).to_csv(dest / 'existing_model_pool_confidence.csv.gz', index=False)
    manifest = dict(job=os.environ['SLURM_JOB_ID'], no_refitting=True, no_labels_modified=True,
                    source_comparisons=len(comparisons),
                    exact_pre_DL_columns=[r for r in comparisons if r['column_role']=='pre-DL named column' and r['n_same_as_current_initial']==r['n_cells']],
                    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (dest / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(pd.DataFrame(comparisons).to_string(index=False), flush=True)
    print(pd.DataFrame(rows).query("category in ['all_DL_pool','all_mismatches','historical_called_to_Unknown']").to_string(index=False), flush=True)
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == '__main__':
    run()
