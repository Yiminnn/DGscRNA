#!/usr/bin/env python3
"""Inventory the staged historical PTC objects without fitting or changing labels."""
import json
import os
from pathlib import Path
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'hvg_ptc_20260916'))
from common import ROOT, OUT, require_slurm, sha, utc, write_json, runtime_record


def plain(value):
    import numpy as np
    if isinstance(value, bytes):
        return value.decode(errors='replace')
    if isinstance(value, np.generic):
        return plain(value.item())
    if isinstance(value, (tuple, list, np.ndarray)):
        return [plain(v) for v in value]
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    return value


def run():
    require_slurm()
    import h5py
    import numpy as np
    import pandas as pd
    from anndata.io import read_elem

    base = OUT / 'ptc_recovery'
    assert (base / 'ARCHIVE_STAGED').exists()
    manifest = json.loads((base / 'archive_manifest.json').read_text())
    assert manifest['status'] == 'completed'
    dest = base / 'inventory'
    dest.mkdir(exist_ok=True)
    archive = base / 'archive'
    started = time.monotonic()
    artifacts = []
    h5paths = sorted(archive.rglob('*.h5ad'))
    h5paths += [ROOT / 'data_ptc/tcr/data_processed/adata_all_raw.h5ad']
    misnamed_h5 = []
    for candidate in sorted(archive.rglob('*.ipynb')):
        with candidate.open('rb') as source:
            if source.read(8) == b'\x89HDF\r\n\x1a\n':
                misnamed_h5.append(candidate)
    h5paths += misnamed_h5
    for i, path in enumerate(h5paths):
        key = f'object_{i:02d}_{path.stem}'
        row = dict(source=str(path), key=key, bytes=path.stat().st_size)
        prior = dest / f'{key}.json'
        if prior.exists():
            saved = json.loads(prior.read_text())
            if saved.get('source') == str(path) and saved.get('bytes') == path.stat().st_size:
                artifacts.append(saved)
                print('Reusing completed schema inventory', key, flush=True)
                continue
        print('Inspecting HDF5 schema', path, flush=True)
        with h5py.File(path, 'r') as handle:
            row['root_keys'] = list(handle.keys())
            row['attrs'] = plain(dict(handle.attrs))
            if 'obs' not in handle or 'var' not in handle:
                row['status'] = 'non_standard_h5ad_requires_schema_review'
                row['root_structure'] = {name: dict(kind=type(obj).__name__,
                    shape=list(obj.shape) if isinstance(obj,h5py.Dataset) else None,
                    keys=list(obj.keys()) if isinstance(obj,h5py.Group) else None,
                    attrs=plain(dict(obj.attrs))) for name,obj in handle.items()}
                artifacts.append(row)
                write_json(dest / f'{key}.json', row)
                print('Nonstandard HDF5 schema recorded for review:', key, flush=True)
                continue
            obs = read_elem(handle['obs'])
            var = read_elem(handle['var'])
            obs.to_csv(dest / f'{key}.obs.csv.gz')
            var.to_csv(dest / f'{key}.var.csv.gz')
            row.update(n_obs=len(obs), n_vars=len(var), obs_unique=bool(obs.index.is_unique),
                       var_unique=bool(var.index.is_unique), obs_columns=list(obs.columns),
                       obs_first_ids=obs.index[:5].astype(str).tolist(),
                       var_first_ids=var.index[:5].astype(str).tolist())
            row['obs_summary'] = {}
            for col in obs.columns:
                vals = obs[col]
                entry = dict(dtype=str(vals.dtype), non_null=int(vals.notna().sum()),
                             n_unique=int(vals.nunique()))
                if vals.nunique() <= 200:
                    entry['counts'] = {str(k): int(v) for k, v in vals.value_counts(dropna=False).items()}
                else:
                    entry['examples'] = vals.dropna().astype(str).head(5).tolist()
                row['obs_summary'][str(col)] = entry
            row['matrices'] = {}
            matrix_paths = ['X']
            if 'raw/X' in handle:
                matrix_paths.append('raw/X')
                raw_var = read_elem(handle['raw/var'])
                if isinstance(raw_var, dict):
                    row['raw_var_legacy_dict_keys'] = list(raw_var)
                    # Legacy SeuratDisk exports can omit AnnData dataframe encoding.
                    if all(isinstance(v,np.ndarray) and v.ndim==1 for v in raw_var.values()):
                        raw_var = pd.DataFrame(raw_var)
                        index_name = handle['raw/var'].attrs.get('_index','_index')
                        if isinstance(index_name,bytes):index_name=index_name.decode()
                        if index_name in raw_var:raw_var=raw_var.set_index(index_name)
                    else:
                        row['raw_var_status'] = 'legacy_nested_dict_requires_schema_review'
                if isinstance(raw_var,pd.DataFrame):
                    raw_var.to_csv(dest / f'{key}.raw_var.csv.gz')
                    row['raw_n_vars'] = len(raw_var)
            if 'layers' in handle:
                matrix_paths += [f'layers/{name}' for name in handle['layers']]
            for name in matrix_paths:
                if name not in handle:
                    continue
                obj = handle[name]
                entry = dict(attrs=plain(dict(obj.attrs)))
                if isinstance(obj, h5py.Group):
                    entry['keys'] = list(obj.keys())
                    entry['shape'] = plain(obj.attrs.get('shape'))
                    if 'data' not in obj:
                        row['matrices'][name] = entry
                        continue
                    data = obj['data']
                    entry.update(dtype=str(data.dtype), stored_values=int(data.size))
                    # Diagnostic samples do not establish full-matrix count semantics.
                    values = np.concatenate([data[:min(50000, data.size)], data[max(0, data.size-50000):]])
                else:
                    entry.update(shape=list(obj.shape), dtype=str(obj.dtype))
                    values = np.asarray(obj[:min(256, obj.shape[0]), :min(256, obj.shape[1])]).ravel()
                entry['diagnostic_sample_only'] = True
                entry['sampled_values'] = int(values.size)
                if values.size:
                    entry.update(sample_min=float(np.nanmin(values)), sample_max=float(np.nanmax(values)),
                                 sample_finite_fraction=float(np.isfinite(values).mean()),
                                 sample_integer_fraction=float(np.isclose(values, np.round(values), atol=1e-6).mean()))
                row['matrices'][name] = entry
        artifacts.append(row)
        write_json(dest / f'{key}.json', row)
        print('Inventoried', key, row['n_obs'], row['n_vars'], flush=True)

    loom_rows = []
    for i, path in enumerate(sorted(archive.rglob('*.loom'))):
        key = f'loom_{i:02d}_{path.parent.name}_{path.stem}'
        row = dict(source=str(path), key=key, bytes=path.stat().st_size)
        with h5py.File(path, 'r') as handle:
            row['root_keys'] = list(handle.keys())
            row['attrs'] = plain(dict(handle.attrs))
            matrix = handle['matrix']
            row['matrix_shape_genes_by_cells'] = list(matrix.shape)
            values = np.asarray(matrix[:min(256, matrix.shape[0]), :min(256, matrix.shape[1])])
            row['matrix_diagnostic_sample_only'] = dict(n_values=int(values.size),
                minimum=float(np.nanmin(values)), maximum=float(np.nanmax(values)),
                integer_fraction=float(np.isclose(values, np.round(values), atol=1e-6).mean()))
            for axis, size in [('col_attrs', matrix.shape[1]), ('row_attrs', matrix.shape[0])]:
                row[axis] = {name: dict(shape=list(obj.shape), dtype=str(obj.dtype))
                             for name, obj in handle[axis].items()}
                columns = {}
                for name, obj in handle[axis].items():
                    if obj.ndim == 1 and len(obj) == size:
                        columns[name] = obj.asstr()[:] if h5py.check_string_dtype(obj.dtype) else obj[:]
                pd.DataFrame(columns).to_csv(dest / f'{key}.{axis}.csv.gz', index=False)
            row['layers'] = {name: dict(shape=list(obj.shape), dtype=str(obj.dtype))
                             for name, obj in handle.get('layers', {}).items()}
        loom_rows.append(row)
        write_json(dest / f'{key}.json', row)
        print('Inventoried', key, row['matrix_shape_genes_by_cells'], flush=True)
    write_json(dest / 'loom_inventory.json', loom_rows)

    table_rows = []
    for i, path in enumerate(sorted(archive.rglob('*.csv'))):
        try:
            table = pd.read_csv(path, low_memory=False)
            info = dict(source=str(path), n_rows=len(table), columns=list(table.columns),
                        first_rows=json.loads(table.head(3).to_json(orient='records')))
            info['low_cardinality'] = {
                str(c): {str(k): int(v) for k, v in table[c].value_counts(dropna=False).items()}
                for c in table.columns if table[c].nunique() <= 30
            }
        except Exception as exc:
            info = dict(source=str(path), error=repr(exc))
        table_rows.append(info)
    write_json(dest / 'csv_inventory.json', table_rows)

    books = []
    workbook_paths = sorted((ROOT / 'paper/submission_v16').glob('Supplementary table*.xlsx'))
    workbook_paths += [archive / 'tcr/scripts/annotations.xlsx']
    for path in workbook_paths:
        excel = pd.ExcelFile(path)
        books.append(dict(source=str(path), sha256=sha(path), sheets={
            name: json.loads(pd.read_excel(excel, sheet_name=name, header=None, nrows=10).to_json(orient='values'))
            for name in excel.sheet_names
        }))
    write_json(dest / 'supplementary_headers.json', books)

    sources = dest / 'notebook_sources'
    sources.mkdir(exist_ok=True)
    notebooks = []
    for i, path in enumerate(sorted(archive.rglob('*.ipynb'))):
        if path in misnamed_h5:
            notebooks.append(dict(source=str(path),status='misnamed_HDF5_not_a_notebook',
                                  original_sha256=sha(path),HDF5_metadata_in_objects=True))
            continue
        try:
            notebook = json.loads(path.read_text())
        except (UnicodeDecodeError,json.JSONDecodeError) as exc:
            notebooks.append(dict(source=str(path),status='invalid_notebook_requires_review',
                original_sha256=sha(path),error=repr(exc)))
            continue
        cells = notebook.get('cells', [])
        target = sources / f'{i:02d}_{path.stem}.txt'
        target.write_text('\n\n'.join(
            f'# SOURCE CELL {j} ({cell.get("cell_type")})\n' + ''.join(cell.get('source', []))
            for j, cell in enumerate(cells)
        ))
        output_path = sources / f'{i:02d}_{path.stem}.outputs.txt'
        outputs = []
        for j, cell in enumerate(cells):
            for output in cell.get('outputs', []):
                value = output.get('text', output.get('data', {}).get('text/plain', []))
                if output.get('output_type') == 'error':
                    value = output.get('traceback', [])
                value = ''.join(value) if isinstance(value, list) else str(value)
                if value:
                    outputs.append(f'# OUTPUT CELL {j}; execution_count={cell.get("execution_count")}\n{value}')
        output_path.write_text('\n\n'.join(outputs))
        notebooks.append(dict(source=str(path), extracted_source=str(target), n_cells=len(cells),
                              extracted_text_outputs=str(output_path),
                              original_sha256=sha(path), source_text_sha256=sha(target),
                              output_text_sha256=sha(output_path)))
    write_json(dest / 'notebooks.json', notebooks)
    write_json(dest / 'manifest.json', dict(status='completed', timestamp=utc(),
        archive_manifest_sha256=sha(base / 'archive_manifest.json'), objects=artifacts,
        source_sha256=sha(Path(__file__)), seconds=time.monotonic()-started,
        **runtime_record()))
    (dest / 'COMPLETE').write_text(sha(dest / 'manifest.json')+'\n')


if __name__ == '__main__':
    run()
