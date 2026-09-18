"""Nested label-free pooled GBM counts for single-run resource measurements only."""
import hashlib
import json
import os
from common import OUT, require_slurm, write_json, complete, checked, sha, utc
SIZES=[10000,30000,50000,100000,120000]

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    dest=OUT/'scalability';dest.mkdir(exist_ok=True)
    if checked(dest,'input_manifest.json','INPUTS_COMPLETE'):return
    cohort=pd.read_csv(OUT/'protocol/cohort.csv');rows=[];genes=set()
    for sample in cohort['sample']:
        src=OUT/'inputs'/sample
        cells=pd.read_csv(src/'cells_fit.csv',dtype=str,keep_default_na=False)[['cell_id']]
        cells['sample']=sample;cells['source_row']=np.arange(len(cells));cells['pooled_id']=sample+'|'+cells.cell_id
        cells['selection_hash']=cells.pooled_id.map(lambda s:hashlib.sha256(('resource_curve_20260917|'+s).encode()).hexdigest())
        rows.append(cells);genes.update(pd.read_csv(src/'genes.csv',dtype=str).gene)
    allcells=pd.concat(rows,ignore_index=True).sort_values('selection_hash',kind='stable').head(max(SIZES))
    allcells['size_order']=np.arange(len(allcells));assert len(allcells)==max(SIZES)
    genes=sorted(genes);lookup={g:i for i,g in enumerate(genes)};blocks=[];metadata=[];hashes={}
    for sample,g in allcells.groupby('sample',sort=True):
        src=OUT/'inputs'/sample;m=json.loads((src/'input_manifest.json').read_text())
        assert checked(src,'input_manifest.json','INPUT_COMPLETE')
        X=sp.csr_matrix((np.fromfile(src/'x.bin',dtype='<f8'),np.fromfile(src/'i.bin',dtype='<i4'),np.fromfile(src/'p.bin',dtype='<i4')),shape=(m['n_cells'],m['n_genes']))
        X=X[g.source_row.to_numpy()].copy()
        conversion=np.asarray([lookup[v] for v in pd.read_csv(src/'genes.csv',dtype=str).gene],dtype=np.int32)
        X.indices=conversion[X.indices];X._shape=(X.shape[0],len(genes));X.sort_indices()
        blocks.append(X);metadata.append(g);hashes[sample]=sha(src/'input_manifest.json')
    full=sp.vstack(blocks,format='csr');obs=pd.concat(metadata,ignore_index=True);del blocks
    order=np.argsort(obs.size_order.to_numpy());full=full[order];obs=obs.iloc[order].reset_index(drop=True)
    manifests=[]
    for n in SIZES:
        sample=f'SCALE_{n}';p=OUT/'inputs'/sample;p.mkdir(exist_ok=True)
        if checked(p,'input_manifest.json','INPUT_COMPLETE'):
            manifests.append(str(p/'input_manifest.json'));continue
        X=full[:n].copy();seen=np.asarray((X>0).sum(axis=0)).ravel()>=3
        X=X[:,seen].tocsr();gg=np.asarray(genes)[seen]
        assert X.nnz<2**31 and (np.asarray(X.sum(axis=1)).ravel()>0).all()
        X.data.astype('<f8').tofile(p/'x.bin');X.indices.astype('<i4').tofile(p/'i.bin');X.indptr.astype('<i4').tofile(p/'p.bin')
        pd.DataFrame({'gene':gg}).to_csv(p/'genes.csv',index=False)
        pd.DataFrame({'cell_id':obs.pooled_id.iloc[:n],'batch':'pooled_resource_test'}).to_csv(p/'cells_fit.csv',index=False)
        obs.iloc[:n].drop(columns='selection_hash').to_csv(p/'source_cell_identity.csv',index=False)
        files={name:sha(p/name) for name in ['x.bin','i.bin','p.bin','genes.csv','cells_fit.csv']}
        write_json(p/'input_manifest.json',dict(status='completed',sample=sample,patient='pooled_runtime_only',primary=False,
            n_cells=n,n_genes=X.shape[1],nnz=X.nnz,fitting_files=files,evaluation_files={},
            input_semantics='Nested barcode-hash sample from all121frozen GBM samples. Original nonnegative integer counts; genes seen in>=3pooled cells. Source identities retained. One uncorrected pooled batch for runtime only.',
            no_accuracy_or_batch_biology_claim=True,no_reference_labels_loaded=True,source_manifests=hashes,
            source_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
        complete(p,'input_manifest.json','INPUT_COMPLETE');manifests.append(str(p/'input_manifest.json'))
        del X
    write_json(dest/'input_manifest.json',dict(status='completed',sizes=SIZES,nested_selection=True,
        label_free_selection=True,purpose='Single-run resource curve, not additional biological validation or a claim that pooled samples are one batch.',
        inputs=manifests,job=os.environ['SLURM_JOB_ID'],source_sha256=sha(__file__),completed_at=utc()))
    complete(dest,'input_manifest.json','INPUTS_COMPLETE')

if __name__=='__main__':run()
