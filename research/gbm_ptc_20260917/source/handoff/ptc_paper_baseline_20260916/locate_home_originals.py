"""Verify the home/work archive and map original final outputs to supplied Sup files.

This is file provenance and exact-label reconciliation, never model fitting or
label remapping. Execute in SLURM because spreadsheet comparisons are numerical.
"""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
import re

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
BASE = ROOT / 'results/hvg_ptc_20260916_v1'
REC = BASE / 'ptc_recovery'
OUT = BASE / 'ptc_paper_baseline'
ARCHIVE = REC / 'archive/tcr'
HOME_ARCHIVE = Path('/users/PCON0080/yimin/work/_archives/tcr.tar.gz')
SCRATCH_ARCHIVE = Path('/fs/scratch/PCON0080/yimin/_tcr_stage/tcr.tar.gz')


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def verify_archives():
    paths = [HOME_ARCHIVE, SCRATCH_ARCHIVE]
    records = []
    with ThreadPoolExecutor(max_workers=2) as pool:
        hashes = list(pool.map(sha, paths))
    for path, digest in zip(paths, hashes):
        records.append(dict(path=str(path), bytes=path.stat().st_size,
                            sha256=digest))
    assert hashes[0] == hashes[1], 'Home archive differs; inspect before assuming identity'
    return dict(status='byte_identical_SHA256', files=records)


def verify_tables():
    import pandas as pd
    metadata = pd.read_csv(ARCHIVE / 'rawdata/metadata.txt', sep='\t')
    samplemap = dict(zip(metadata.sc_ID, metadata.Sample))
    samplemap.update({s: s for s in metadata.Sample})
    s2path = ROOT / 'paper/submission_v16/Supplementary table S2_T cell type identification_DGscRNA.xlsx'
    s3path = ROOT / 'paper/submission_v16/Supplementary table S3_Method comparison in cell type annotation.xlsx'
    s2 = pd.read_excel(s2path, header=3, keep_default_na=False)
    s2.columns = s2.columns.str.strip()
    s2 = s2[s2.Sample_ID.str.match(r'^(?:MT|TU|N|T)-[12]_[ACGT]+$')].set_index('Sample_ID')
    assert s2.index.is_unique and len(s2) == 92404

    def canonical(frame):
        frame = frame.copy()
        frame.columns = frame.columns.astype(str).str.strip()
        fields = [c for c in ['Sample_ID', 'CellID', 'Cell ID', 'cell_id', '_index',
                             'Unnamed: 0', 'index', 'barcode', 'X'] if c in frame]
        fields += [c for c in frame.columns[:2] if c not in fields]
        samples = [c for c in ['sample.name', 'Sample.name', 'sampleid', 'sample',
                              'orig.ident', 'Sample', 'sample_id'] if c in frame]
        candidates = []
        for col in fields:
            values = frame[col].astype(str)
            ids = values.where(values.str.match(r'^(?:MT|TU|N|T)-[12]_[ACGT]+$'))
            candidates.append((int(ids.isin(s2.index).sum()), ids, f'direct:{col}'))
            barcode = values.str.extract(r'([ACGT]{12,})', expand=False)
            for sc in samples:
                ids = frame[sc].map(samplemap) + '_' + barcode
                candidates.append((int(ids.isin(s2.index).sum()), ids, f'{sc}+barcode({col})'))
        count, ids, rule = max(candidates, key=lambda x: x[0])
        assert count == len(s2), (count, rule)
        frame.index = ids
        frame = frame[frame.index.notna()]
        assert frame.index.is_unique
        return frame.loc[s2.index], rule

    s3 = pd.read_excel(s3path, header=5, keep_default_na=False)
    s3.columns = s3.columns.str.strip()
    s3 = s3[s3['sample.name'].isin(samplemap)]
    s3, s3rule = canonical(s3)
    records, inputs, joins = [], [s2path, s3path], []

    def compare(supname, target, path, pairs, member):
        inputs.append(path)
        source = pd.read_excel(path, keep_default_na=False) if path.suffix == '.xlsx' else pd.read_csv(path, keep_default_na=False)
        source, rule = canonical(source)
        joins.append(dict(local_file=str(path), archive_member=member,
                          id_rule=rule, n_cells=len(source)))
        for target_col, source_col in pairs:
            expected, observed = target[target_col], source[source_col]
            if target_col in ['is_T_cell_Real', 'Prediction True and False'] or target_col.endswith('_ct_T_cells'):
                expected = pd.to_numeric(expected)
                observed = pd.to_numeric(observed)
            else:
                expected = expected.astype(str)
                observed = observed.astype(str)
            equal = expected.eq(observed)
            records.append(dict(supplement=supname, supplement_column=target_col,
                original_file=str(path), original_archive_member=member, original_column=source_col,
                n_cells=len(equal), n_exact=int(equal.sum()), n_mismatch=int((~equal).sum()),
                transformation='cell barcode/sample alignment only; no label mapping'))
        return source

    s2pairs = [('Cell_type_annotation_CellMarker2.0', 'cell_type_annotation_CellMarker2.0'),
               ('DG_scRNA_Finalized_Cell_Types', 'DG_scRNA_Finalized_Cell_Types'),
               ('Cell_types_general', 'cell_types_general'),
               ('is_T_cell_Real', 'is_T_cell_Real')]
    annotation = compare('S2', s2, ARCHIVE / 'scripts/annotations.xlsx', s2pairs, 'tcr/scripts/annotations.xlsx')
    native_pairs = [(a, 'DGCyTOF_Finalized_Cell_Types' if b == 'DG_scRNA_Finalized_Cell_Types' else b)
                    for a, b in s2pairs if a != 'is_T_cell_Real']
    compare('S2', s2, REC / 'inventory_workspace/object_01.meta.csv', native_pairs,
            'tcr/rawdata/integrated_data_final_annotation.Rdata (metadata exported in SLURM)')
    compare('S2', s2, REC / 'inventory/loom_02_rawdata_integrated_data.col_attrs.csv.gz', native_pairs,
            'tcr/rawdata/integrated_data.loom (column attributes exported in SLURM)')
    s3pairs = [('DGscRNA_prediction_cell_type', 'Final_DGCyTOF'),
               ('DGscRNA_prediction_detailed_T_cell', 'Final_DGCyTOF_General_w_T_cell')]
    old = compare('S3', s3, ARCHIVE / 'rawdata/data_with_validation.csv', s3pairs,
                  'tcr/rawdata/data_with_validation.csv')
    compare('S3', s3, ARCHIVE / 'rawdata/data_with_validation+3_cell_types 2.csv', s3pairs,
            'tcr/rawdata/data_with_validation+3_cell_types 2.csv')
    result = pd.DataFrame(records)
    result.to_csv(OUT / 'home_original_sup_column_crosswalk.csv', index=False)
    assert result.n_mismatch.eq(0).all(), result.to_string(index=False)
    # Preserve both final endpoints verbatim: S2 and S3 contain distinct selections.
    references = pd.DataFrame(index=s2.index)
    references['S2_final_native'] = s2.DG_scRNA_Finalized_Cell_Types
    references['S2_final_general'] = s2.Cell_types_general
    references['S3_final_native'] = s3.DGscRNA_prediction_cell_type
    references['S3_final_detailed'] = s3.DGscRNA_prediction_detailed_T_cell
    references['original_validation_t_cell'] = old.validation_t_cell
    references['S2_is_T_cell_Real'] = s2.is_T_cell_Real
    references['S3_T_cell_as_supplied'] = s3.T_cell
    references.to_csv(OUT / 'original_final_sup_references.csv.gz', index_label='cell_id')
    return dict(status='all_selected_original_final_columns_exact', n_cells=len(s2),
                n_compared_columns=len(result), joins=joins, S3_id_rule=s3rule,
                user_confirmation='Sup files are original final used labels; retain both S2 and S3 endpoints',
                model_weights_recovered=False,
                inputs=[dict(path=str(p), sha256=sha(p)) for p in dict.fromkeys(inputs)],
                outputs=['home_original_sup_column_crosswalk.csv', 'original_final_sup_references.csv.gz'])


def run():
    assert os.environ.get('SLURM_JOB_ID'), 'Run through SLURM'
    # Archive reads and table checks run concurrently inside the allocation.
    with ThreadPoolExecutor(max_workers=2) as pool:
        archive_future = pool.submit(verify_archives)
        tables = verify_tables()
        print(json.dumps(tables, indent=2), flush=True)
        archives = archive_future.result()
    report = dict(job=os.environ['SLURM_JOB_ID'], completed_utc=datetime.now(timezone.utc).isoformat(),
                  archive_identity=archives, table_identity=tables, script_sha256=sha(Path(__file__)))
    (OUT / 'home_original_sup_manifest.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(archives, indent=2), flush=True)


if __name__ == '__main__':
    run()
