"""Preserve the full cohort and add a CCA-feasible, donor-size-only sensitivity."""
import copy,fcntl,hashlib,json,os,shutil
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'
UNIT='colorectal_CCA28_ge31'

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for chunk in iter(lambda:f.read(8*1024*1024),b''):h.update(chunk)
    return h.hexdigest()

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    source=OUT/'inputs/colorectal/colorectal';dest=source.parent/UNIT
    dest.mkdir(exist_ok=True)
    parent=json.loads((source/'input_manifest.json').read_text())
    fit=pd.read_csv(source/'cells_fit.csv',dtype=str)
    sizes=fit.batch.value_counts()
    keep=fit.batch.map(sizes).ge(31).to_numpy()
    excluded=fit.loc[~keep].copy()
    assert len(excluded)==3 and set(excluded.batch)=={'HTA8_6004'}
    assert int(keep.sum())==47104 and fit.loc[keep,'batch'].nunique()==28
    if not (dest/'EXPORT_COMPLETE').exists():
        x=sp.csr_matrix((np.memmap(source/'x.bin',dtype='<f8',mode='r'),
                        np.memmap(source/'i.bin',dtype='<i4',mode='r'),
                        np.memmap(source/'p.bin',dtype='<i4',mode='r')),
                       shape=(parent['n_cells'],parent['n_genes']),copy=False)[keep].tocsr()
        x.data.astype('<f8',copy=False).tofile(dest/'x.bin')
        x.indices.astype('<i4',copy=False).tofile(dest/'i.bin')
        x.indptr.astype('<i4',copy=False).tofile(dest/'p.bin')
        shutil.copy2(source/'genes.csv',dest/'genes.csv')
        fit.loc[keep].to_csv(dest/'cells_fit.csv',index=False)
        truth=pd.read_csv(source/'evaluation_only.csv.gz',keep_default_na=False,dtype=str)
        assert truth.cell_id.tolist()==fit.cell_id.tolist()
        truth.loc[keep].to_csv(dest/'evaluation_only.csv.gz',index=False)
        excluded['reason']='donor has fewer than 31 cells, insufficient for the reference 30-dimensional CCA'
        excluded.to_csv(dest/'excluded_small_donor_cells.csv',index=False)
        record=copy.deepcopy(parent)
        record.update(unit=UNIT,path=str(dest),scope='CCA_feasibility_donors_ge31',n_cells=x.shape[0],nnz=x.nnz,
            batch_sizes=fit.loc[keep,'batch'].value_counts().to_dict(),derived_control_of='colorectal',
            curator_cell_set_preserved=False,original_full_cohort_preserved_at=str(source),
            exclusion_rule='Keep donor groups with at least 31 curated cells, before inspecting annotation performance.',
            excluded_donor='HTA8_6004',excluded_cells=3,parent_manifest_sha256=sha(source/'input_manifest.json'),
            script_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],
            files={name:sha(dest/name) for name in ['x.bin','i.bin','p.bin','genes.csv','cells_fit.csv']})
        (dest/'input_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
        (dest/'EXPORT_COMPLETE').write_text(sha(dest/'input_manifest.json')+'\n')
    marker_source=OUT/'markers/colorectal.json';marker_target=OUT/'markers'/f'{UNIT}.json'
    shutil.copy2(marker_source,marker_target);assert sha(marker_source)==sha(marker_target)
    roster_path=OUT/'markers/marker_roster.json'
    with (OUT/'dispatcher.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        roster=json.loads(roster_path.read_text())
        if not any(r['unit']==UNIT for r in roster['units']):
            backup=roster_path.with_name('marker_roster_before_colorectal_CCA_control.json')
            if not backup.exists():shutil.copy2(roster_path,backup)
            row=copy.deepcopy(next(r for r in roster['units'] if r['unit']=='colorectal'))
            row.update(unit=UNIT,scope='CCA_feasibility_donors_ge31',derived_control_of='colorectal',
                selection_basis='Byte-identical pre-existing colorectal marker contexts; donor-size-only CCA feasibility control, no marker reselection.')
            roster['units'].append(row)
            temp=roster_path.with_suffix('.part');temp.write_text(json.dumps(roster,indent=2)+'\n');temp.replace(roster_path)
    print('EXPORTED',UNIT,47104,'cells; 28 donors; unchanged marker libraries',flush=True)

if __name__=='__main__':run()
