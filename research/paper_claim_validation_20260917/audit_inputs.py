"""Audit the frozen GSE274546 cohort; export label-free sparse input for native R."""
import argparse
import gc
import json
import os
from pathlib import Path
from common import ROOT, OUT, INPUTS, L1, require_slurm, samples, sha, write_json, complete, checked, utc

def audit(sample):
    require_slurm()
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    src = INPUTS/sample
    dest = OUT/'inputs'/sample
    dest.mkdir(parents=True, exist_ok=True)
    if checked(dest, 'input_manifest.json', 'INPUT_COMPLETE'):
        print('INPUT_CACHED', sample, flush=True)
        return
    m = json.loads((src/'manifest.json').read_text())
    # Historical input completion is a text receipt, not a manifest digest.
    # Its artifacts are verified against the frozen manifest below.
    assert (src/'COMPLETE').is_file()
    assert m['status']=='completed' and m['HVG_selection'] is False and m['gold_used_for_fitting'] is False
    for name in ['counts_gene_filtered.h5ad', 'labels.csv', 'lfine22.json', 'qc_cells.csv']:
        assert sha(src/name)==m['artifacts'][name]['sha256'], (sample, name)
    rawobs = ROOT/'data_bench/GSE274546/mtx'/sample/'obs.csv'
    assert sha(rawobs)==m['sha256'][str(rawobs)]
    obs = pd.read_csv(rawobs, keep_default_na=False, dtype=str)
    labels = pd.read_csv(src/'labels.csv', keep_default_na=False, dtype=str)
    a = ad.read_h5ad(src/'counts_gene_filtered.h5ad')
    assert a.shape==(m['n_cells'],m['n_genes_retained'])
    assert not len(a.obs.columns) and not len(a.var.columns) and not a.uns
    assert not len(a.layers) and not len(a.obsm) and not len(a.obsp) and a.raw is None
    assert a.obs_names.is_unique and a.var_names.is_unique
    assert list(a.obs_names)==list(obs.CellID)==list(labels.CellID)
    assert set(obs.Sample)=={sample} and obs.Patient.nunique()==1 and not obs.Patient.isin(['','NA']).any()
    assert set(obs.L1)<=set(L1), set(obs.L1)-set(L1)
    x = sp.csr_matrix(a.X, dtype=np.float64)
    x.sum_duplicates(); x.eliminate_zeros(); x.sort_indices()
    assert np.isfinite(x.data).all() and (x.data>=0).all() and np.array_equal(x.data, np.round(x.data))
    assert (x.getnnz(axis=0)>=3).all() and (np.asarray(x.sum(axis=1)).ravel()>0).all()
    assert x.nnz < np.iinfo(np.int32).max
    x.data.astype('<f8').tofile(dest/'x.bin')
    x.indices.astype('<i4').tofile(dest/'i.bin')
    x.indptr.astype('<i4').tofile(dest/'p.bin')
    pd.DataFrame({'cell_id':a.obs_names.astype(str), 'batch':sample}).to_csv(dest/'cells_fit.csv',index=False)
    pd.DataFrame({'gene':a.var_names.astype(str)}).to_csv(dest/'genes.csv',index=False)
    ed = OUT/'evaluation_inputs'/sample
    ed.mkdir(parents=True, exist_ok=True)
    truth = obs.copy().rename(columns={'CellID':'cell_id'})
    truth['lfine_original'] = labels.lfine_original.to_numpy()
    truth['lfine_eval23'] = labels.lfine_eval23.to_numpy()
    truth.to_csv(ed/'truth.csv.gz',index=False,compression='gzip')
    qc = pd.read_csv(src/'qc_cells.csv')
    assert list(qc.CellID)==list(a.obs_names)
    qc.to_csv(ed/'qc.csv.gz',index=False,compression='gzip')
    record = dict(status='completed', sample=sample, patient=obs.Patient.iloc[0],
        n_cells=a.n_obs,n_genes=a.n_vars,nnz=int(x.nnz), primary=bool(m['evaluable']),
        source_manifest_sha256=sha(src/'manifest.json'), rawobs_sha256=sha(rawobs),
        counts_sha256=sha(src/'counts_gene_filtered.h5ad'),
        input_semantics='author-filtered nonnegative integer RNA counts; genes detected in >=3 cells',
        truth_labels_excluded_from_fit=True, cell_order_identical=True,
        primary_rule=m['evaluability_rule'], primary_rule_reused_not_reselected=True,
        L1_counts={k:int(v) for k,v in obs.L1.value_counts().items()},
        malignant_state_counts={k:int(v) for k,v in obs.loc[obs.L1.eq('Malignant'),'MalState'].value_counts().items()},
        fitting_files={p.name:sha(p) for p in [dest/n for n in ['x.bin','i.bin','p.bin','cells_fit.csv','genes.csv']]},
        evaluation_files={p.name:sha(p) for p in ed.iterdir() if p.is_file()},
        source_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],completed_at=utc())
    write_json(dest/'input_manifest.json',record)
    complete(dest,'input_manifest.json','INPUT_COMPLETE')
    print('INPUT_COMPLETE',sample,a.shape,record['patient'],record['primary'],flush=True)
    del a,x,obs,labels,qc
    gc.collect()

def aggregate():
    require_slurm()
    import pandas as pd
    records=[]
    for sample in samples():
        d=OUT/'inputs'/sample
        assert checked(d,'input_manifest.json','INPUT_COMPLETE'),sample
        records.append(json.loads((d/'input_manifest.json').read_text()))
    t=pd.DataFrame([{k:r[k] for k in ['sample','patient','n_cells','n_genes','nnz','primary']} for r in records])
    assert len(t)==121 and t.patient.nunique()==59
    primary=t[t.primary].sort_values(['n_cells','sample'])
    assert len(primary)==97 and primary.patient.nunique()==55
    pilot=primary.iloc[[0,len(primary)//2,len(primary)-1]].copy()
    pilot['selection']=['minimum_cell_count','median_cell_count','maximum_cell_count']
    d=OUT/'protocol';d.mkdir(parents=True,exist_ok=True)
    t.to_csv(d/'cohort.csv',index=False)
    pilot.to_csv(d/'pilot.csv',index=False)
    (d/'pilot_samples.txt').write_text('\n'.join(pilot['sample'])+'\n')
    t.sort_values(['primary','n_cells','sample'],ascending=[False,True,True]).to_csv(d/'sample_order.csv',index=False)
    # One patient always stays in one fold. Folds use ID only, never outcomes.
    import hashlib
    patients=sorted(t.patient.unique(),key=lambda p:hashlib.sha256(('dgscrna-R-reference-20260917:'+p).encode()).hexdigest())
    fold={p:i%5 for i,p in enumerate(patients)}
    t.assign(fold=t.patient.map(fold)).to_csv(d/'patient_folds.csv',index=False)
    manifest=dict(status='completed',n_samples=len(t),n_patients=t.patient.nunique(),
        primary_samples=len(primary),primary_patients=primary.patient.nunique(),n_cells=int(t.n_cells.sum()),
        pilots=pilot.to_dict('records'),display_sample=pilot.iloc[1]['sample'],
        pilot_selection_used_prediction_scores=False,
        cohort_selection='Reuse historical evaluability mask; report all121 separately; no ranking-based exclusions',
        evaluation='Author L1 strict11 macro-F1 over present truth classes; fixed11 and collapsed-neuron10 secondary. Unknown/unmapped are errors on all cells. Malignant states separate.',
        fit_scope='Each sample independently normalized and fitted. Single batch; no CCA or patient pooling.',
        retrospective_selection='Frozen 5-fold patient-label-heldout retrospective evaluation; cohort previously explored',
        per_sample_input_sha256={s:sha(OUT/'inputs'/s/'input_manifest.json') for s in samples()},
        source_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],completed_at=utc())
    write_json(d/'input_audit.json',manifest)
    complete(d,'input_audit.json','INPUT_AUDIT_COMPLETE')
    print(json.dumps(manifest,ensure_ascii=False),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['audit','aggregate']);p.add_argument('--shards',type=int,default=4)
    a=p.parse_args()
    if a.mode=='aggregate':aggregate()
    else:
        i=int(os.environ.get('SLURM_ARRAY_TASK_ID','0'))
        for s in samples()[i::a.shards]:audit(s)
