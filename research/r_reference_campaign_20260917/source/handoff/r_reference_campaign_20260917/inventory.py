"""Read dataset schemas and marker metadata only inside a SLURM allocation."""
import hashlib
import json
import os
from pathlib import Path

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT = ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'

def main():
    assert os.environ.get('SLURM_JOB_ID')
    import anndata as ad
    import numpy as np
    import pandas as pd
    OUT.mkdir(parents=True, exist_ok=True)
    paths = {k: ROOT/f'data_bench/{k}/{k}.h5ad' for k in
             ['brain_GBM','breast_TNBC','colorectal','kidney_ccRCC','blood_DLBCL']}
    paths.update({k: ROOT/f'data_deck/{k}/{k}.h5ad' for k in
                  ['baron_human','muraro','segerstolpe','xin','immune_ALL_human','HCL']})
    records = {}
    for key, path in paths.items():
        a = ad.read_h5ad(path, backed='r')
        meta = {}
        for c in a.obs:
            if any(t in c.lower() for t in ['batch','sample','donor','patient','tissue','gt','cell_type','annotation','stage']):
                counts = a.obs[c].astype(str).value_counts()
                meta[c] = {'n_unique':len(counts), 'values':counts.to_dict() if len(counts)<=120 else counts.head(8).to_dict()}
        record = dict(path=str(path), shape=list(a.shape), obs_columns=list(a.obs),
                      metadata=meta, var_columns=list(a.var), genes=a.var_names[:5].tolist(),
                      layers=list(a.layers), raw_shape=list(a.raw.shape) if a.raw is not None else None,
                      uns_keys=list(a.uns), obsm_keys=list(a.obsm), sample_X_dtype=str(a.X.dtype))
        if a.raw is not None:
            record['raw_var_columns']=list(a.raw.var)
            record['raw_genes']=a.raw.var_names[:5].tolist()
        records[key] = record
        print(key, a.shape, flush=True)
        a.file.close()
    (OUT/'dataset_inventory.json').write_text(json.dumps(records,indent=2,default=str)+'\n')
    marker_file=ROOT/'handoff/refdb/Cell_marker_Human.xlsx'
    markers=pd.read_excel(marker_file)
    summary={'source':str(marker_file),'sha256':hashlib.sha256(marker_file.read_bytes()).hexdigest(),
             'columns':list(markers),'n_rows':len(markers),'job':os.environ['SLURM_JOB_ID']}
    for col in ['species','tissue_class','tissue_type','cell_type','cancer_type']:
        if col in markers:
            summary[col]=markers[col].fillna('<NA>').astype(str).value_counts().to_dict()
    (OUT/'cellmarker_inventory.json').write_text(json.dumps(summary,indent=2,default=str)+'\n')
    print('CellMarker',summary['n_rows'],list(markers),flush=True)

if __name__=='__main__': main()
