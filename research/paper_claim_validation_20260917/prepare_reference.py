"""Fold-specific labelled references; heldout-patient labels are never opened here."""
import hashlib
import json
import os
import sys
from common import OUT, L1, require_slurm, checked, complete, sha, write_json, utc

def run(fold):
    require_slurm()
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    fold=int(fold);dest=OUT/'reference_inputs'/f'fold{fold}';dest.mkdir(parents=True,exist_ok=True)
    if checked(dest):return
    cohort=pd.read_csv(OUT/'protocol/cohort.csv')
    folds=pd.read_csv(OUT/'protocol/patient_folds.csv')[['patient','fold']].drop_duplicates()
    cohort=cohort.merge(folds,on='patient',validate='many_to_one')
    training=cohort[cohort.primary & cohort.fold.ne(fold)]
    tests=cohort[cohort.fold.eq(fold)]
    assert not set(training.patient)&set(tests.patient)
    meta=[];hashes={}
    for row in training.itertuples():
        p=OUT/'evaluation_inputs'/row.sample/'truth.csv.gz'
        im=json.loads((OUT/'inputs'/row.sample/'input_manifest.json').read_text())
        assert sha(p)==im['evaluation_files']['truth.csv.gz']
        t=pd.read_csv(p,dtype=str,keep_default_na=False)[['cell_id','L1']]
        t['sample']=row.sample;t['patient']=row.patient;t['source_row']=np.arange(len(t))
        t['reference_id']=row.sample+'|'+t.cell_id
        t['selection_hash']=t.reference_id.map(lambda v:hashlib.sha256(('SingleR-reference-v1|'+v).encode()).hexdigest())
        meta.append(t);hashes[row.sample]=sha(p)
    all_meta=pd.concat(meta,ignore_index=True)
    # Reference information is deliberately balanced and bounded independently
    # of tool performance. Rare classes retain every available training cell.
    selected=all_meta.sort_values('selection_hash',kind='stable').groupby(['patient','L1'],sort=True).head(100)
    selected=selected.sort_values(['sample','source_row'],kind='stable')
    genes=set()
    for sample in training['sample']:
        genes.update(pd.read_csv(OUT/'inputs'/sample/'genes.csv',dtype=str).gene)
    genes=sorted(genes);index={g:i for i,g in enumerate(genes)}
    blocks=[];records=[];sources={}
    for sample,g in selected.groupby('sample',sort=True):
        src=OUT/'inputs'/sample;im=json.loads((src/'input_manifest.json').read_text())
        assert checked(src,'input_manifest.json','INPUT_COMPLETE')
        data=np.fromfile(src/'x.bin',dtype='<f8');col=np.fromfile(src/'i.bin',dtype='<i4');ptr=np.fromfile(src/'p.bin',dtype='<i4')
        X=sp.csr_matrix((data,col,ptr),shape=(im['n_cells'],im['n_genes']))
        X=X[g.source_row.to_numpy()].copy()
        totals=np.asarray(X.sum(axis=1)).ravel();assert (totals>0).all()
        X.data=np.log1p(X.data*np.repeat(10000/totals,np.diff(X.indptr)))
        local=pd.read_csv(src/'genes.csv',dtype=str).gene
        conversion=np.asarray([index[v] for v in local],dtype=np.int32)
        X.indices=conversion[X.indices];X._shape=(X.shape[0],len(genes));X.sort_indices()
        blocks.append(X);records.append(g);sources[sample]=sha(src/'input_manifest.json')
    X=sp.vstack(blocks,format='csr');obs=pd.concat(records,ignore_index=True)
    assert X.shape==(len(obs),len(genes)) and not obs.reference_id.duplicated().any()
    assert set(obs.patient).isdisjoint(tests.patient) and set(obs.L1)<=set(L1)
    X.data.astype('<f8').tofile(dest/'x.bin');X.indices.astype('<i4').tofile(dest/'i.bin');X.indptr.astype('<i4').tofile(dest/'p.bin')
    obs.drop(columns='selection_hash').to_csv(dest/'reference_cells.csv',index=False)
    pd.DataFrame({'gene':genes}).to_csv(dest/'genes.csv',index=False)
    tests[['sample','patient','fold','primary']].to_csv(dest/'test_samples.csv',index=False)
    write_json(dest/'manifest.json',dict(status='completed',fold=fold,n_cells=X.shape[0],n_genes=X.shape[1],nnz=X.nnz,
        training_samples=training['sample'].tolist(),training_patients=sorted(set(training.patient)),
        heldout_patients=sorted(set(tests.patient)),reference_truth_sha256=hashes,input_manifests=sources,
        support=obs.L1.value_counts().to_dict(),selection='At most100cells per training-patient x L1 class by frozen cell-ID hash. All rare-class cells kept. Training labels from primary97 only; all samples of heldout patients excluded.',
        normalization='log1p counts / eligible-gene library total x10000, before selecting reference cells',
        supervision='Curated training-patient labels; distinct information condition from marker-only DG-scRNA.',
        files={p.name:sha(p) for p in dest.iterdir() if p.is_file()},job=os.environ['SLURM_JOB_ID'],source_sha256=sha(__file__),completed_at=utc()))
    complete(dest)

if __name__=='__main__':run(sys.argv[1] if len(sys.argv)>1 else os.environ['SLURM_ARRAY_TASK_ID'])
