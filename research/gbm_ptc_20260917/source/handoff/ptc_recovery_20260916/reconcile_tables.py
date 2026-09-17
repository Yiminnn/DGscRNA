#!/usr/bin/env python3
"""Identify per-cell provenance of supplied PTC tables without changing their labels.

This is a diagnostic source audit, not a fitted annotation model or a passed
reproduction gate. Exact agreement candidates must still be traced to source code.
"""
import gc
import json
from pathlib import Path
import re
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'hvg_ptc_20260916'))
from common import ROOT, OUT, require_slurm, sha, utc, write_json, runtime_record


def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    base = OUT / 'ptc_recovery'
    assert (base / 'inventory/COMPLETE').exists()
    dest = base / 'table_reconciliation'
    dest.mkdir(exist_ok=True)
    archive = base / 'archive'
    started = time.monotonic()
    meta_path = archive / 'tcr/rawdata/metadata.txt'
    meta = pd.read_csv(meta_path, sep='\t').set_index('Sample')
    sample_map = dict(zip(meta.sc_ID.astype(str), meta.index.astype(str)))
    sample_map.update({s: s for s in meta.index})
    meta.to_csv(dest / 'original_sample_patient_map.csv')
    s2path = ROOT / 'paper/submission_v16/Supplementary table S2_T cell type identification_DGscRNA.xlsx'
    s2 = pd.read_excel(s2path, header=3)
    s2.columns = s2.columns.astype(str).str.strip()
    assert 'Sample_ID' in s2
    s2 = s2.loc[s2['Sample_ID'].notna()].copy()
    s2 = s2.set_index('Sample_ID')
    assert s2.index.is_unique
    s2.to_csv(dest / 'S2_as_supplied.csv.gz')
    targets = [c for c in ['Cell_type_annotation_CellMarker2.0', 'DG_scRNA_Finalized_Cell_Types', 'Cell_types_general',
                           'is_T_cell_Real', 'Prediction True and False'] if c in s2]
    assert len(targets) == 5, s2.columns.tolist()
    id_pattern = re.compile(r'^(?:MT|TU|N|T)-[12]_[ACGT]+$')

    def canonical_ids(frame):
        # Source notebooks explicitly join sample metadata with the barcode prefix.
        idcols = [c for c in ['Sample_ID', 'CellID', 'Cell ID', 'cell_id', '_index',
                              'Unnamed: 0', 'index', 'barcode'] if c in frame]
        idcols.extend([c for c in frame.columns[:2] if c not in idcols])
        samplecols = [c for c in ['sample.name', 'Sample.name', 'sampleid', 'sample',
                                  'orig.ident', 'Sample', 'sample_id'] if c in frame]
        candidates = []
        for c in idcols:
            vals = frame[c].astype(str)
            direct = vals.str.match(id_pattern)
            if direct.any():
                ids = vals.where(direct)
                candidates.append((int(ids.isin(s2.index).sum()), ids, f'direct:{c}'))
            barcode = vals.str.extract(r'([ACGT]{12,})', expand=False)
            for sc in samplecols:
                samples = frame[sc].astype(str).map(sample_map)
                ids = samples + '_' + barcode
                candidates.append((int(ids.isin(s2.index).sum()), ids, f'{sc}+barcode({c})'))
        if not candidates:
            return None, 'no_supported_id_fields'
        # Selection uses cell IDs alone, never annotation agreement.
        count, ids, rule = max(candidates, key=lambda x: x[0])
        return (ids if count else None), rule

    paths = [archive / 'tcr/rawdata/data_with_validation.csv',
             archive / 'tcr/rawdata/data_with_validation+3_cell_types 2.csv',
             archive / 'tcr/rawdata/subtypes_032523.csv']
    inventory_manifest = json.loads((base / 'inventory/manifest.json').read_text())
    for obj in inventory_manifest['objects']:
        # Several public validation cohorts reuse P1/T1-style names and 10x barcodes.
        # Their coincidental barcode matches are not evidence of this cohort.
        if '/ptc_val/' not in obj['source']:
            path = base / 'inventory' / f'{obj["key"]}.obs.csv.gz'
            if path.exists():paths.append(path)
    for obj in json.loads((base / 'inventory/loom_inventory.json').read_text()):
        if '/rawdata/' in obj['source']:
            paths.append(base / 'inventory' / f'{obj["key"]}.col_attrs.csv.gz')
    paths += [archive / 'tcr/scripts/annotations.xlsx']
    paths += sorted((base / 'inventory_workspace').glob('*.meta.csv'))
    joins, comparisons = [], []
    for i, path in enumerate(dict.fromkeys(paths)):
        print('Auditing table', path, flush=True)
        frame = pd.read_excel(path) if path.suffix == '.xlsx' else pd.read_csv(path, low_memory=False)
        frame.columns = frame.columns.astype(str).str.strip()
        ids, rule = canonical_ids(frame)
        item = dict(source=str(path), key=f'table_{i:02d}', rows=len(frame),
                    columns=len(frame.columns), id_rule=rule, source_sha256=sha(path))
        if ids is None:
            item['status'] = 'no_target_cohort_overlap'
            joins.append(item)
            continue
        item.update(n_valid_ids=int(ids.notna().sum()), n_duplicate_valid_ids=int(ids[ids.notna()].duplicated().sum()))
        if item['n_duplicate_valid_ids']:
            item['status'] = 'duplicate_ids_requires_resolution'
            joins.append(item)
            continue
        frame.index = ids
        frame = frame.loc[frame.index.notna()]
        common = s2.index.intersection(frame.index)
        item.update(status='joined', n_joined=len(common), n_missing_S2=len(s2)-len(common),
                    n_extra=int((~frame.index.isin(s2.index)).sum()))
        joins.append(item)
        if not len(common):
            continue
        candidate_columns = [c for c in frame.columns
            if re.search(r'annot|cell.?type|DGscRNA|DGCyTOF|is_T_cell|label|validation_t_cell|^T_cells$', c, re.I)]
        table_matches = []
        for target in targets:
            expected = s2.loc[common, target]
            numeric = target in ['is_T_cell_Real', 'Prediction True and False']
            if numeric:
                expected = pd.to_numeric(expected, errors='coerce')
            else:
                expected = expected.astype('string')
            for col in candidate_columns:
                observed = frame.loc[common, col]
                observed = pd.to_numeric(observed, errors='coerce') if numeric else observed.astype('string')
                both = expected.notna() & observed.notna()
                equal = (expected == observed).fillna(False)
                count = int(equal.sum())
                if count == 0:
                    continue
                row = dict(source=str(path), key=item['key'], target=target, column=col,
                    n_joined=len(common), n_both_nonmissing=int(both.sum()),
                    n_exact=count, n_disagree=int((both & ~equal).sum()),
                    exact_fraction=count / len(common),
                    all_S2_cells_exact=bool(count==len(s2)),
                    terminal_name_hint=bool(re.search(r'DGscRNA|DGCyTOF|Final', col)))
                comparisons.append(row)
                table_matches.append(row)
                if row['all_S2_cells_exact']:
                    frame.loc[s2.index, [col]].to_csv(dest / f'{item["key"]}.{len(table_matches):04d}.exact_column.csv.gz')
        pd.DataFrame(comparisons).to_csv(dest / 'candidate_column_agreement.csv', index=False)
        write_json(dest / 'table_joins.json', joins)
        print('Joined', len(common), 'S2 cells; exact full-cohort columns',
              sum(r['all_S2_cells_exact'] for r in table_matches), flush=True)
        del frame
        gc.collect()
    write_json(dest / 'table_joins.json', joins)
    result = pd.DataFrame(comparisons)
    result.to_csv(dest / 'candidate_column_agreement.csv', index=False)
    best = result.sort_values(['target','exact_fraction','n_joined'], ascending=[True,False,False]).groupby('target').head(15)
    best.to_csv(dest / 'highest_agreement_candidates.csv', index=False)
    write_json(dest / 'manifest.json', dict(status='diagnostic_completed_not_reproduction_gate',
        timestamp=utc(), n_S2=len(s2), targets=targets, n_source_tables=len(joins),
        S2_sha256=sha(s2path), metadata_sha256=sha(meta_path),
        source_sha256=sha(Path(__file__)), seconds=time.monotonic()-started,
        **runtime_record()))
    (dest / 'COMPLETE').write_text(sha(dest / 'manifest.json')+'\n')


if __name__ == '__main__':
    run()
