"""GBM-only corrected-input replay; prepare -> unchanged predictor -> evaluate."""
from pathlib import Path
import os,sys,json,subprocess,time
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
sys.path.insert(0,str(ROOT/'handoff/paper_claim_validation_20260917'))
from common import OUT,L1,require_slurm,sha,write_json,utc,checked,complete
from prediction_helpers import map_labels
CODE=ROOT/'handoff/reviewer_completion_20260920/comparison'
DEST=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/comparison/scDeepSort_LogNormalize'
MODEL=Path('/fs/scratch/PCON0080/yimin/tools/deepsort-1.0/deepsort-pretrained')

def run_sample(sample):
    require_slurm()
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    from sklearn.metrics import f1_score,precision_recall_fscore_support,confusion_matrix
    src=OUT/'inputs'/sample;prep=OUT/'GBM'/sample/'hvg2000';dest=DEST/sample;dest.mkdir(parents=True,exist_ok=True)
    if checked(dest):return json.loads((dest/'manifest.json').read_text())
    started=time.monotonic();im=json.loads((src/'input_manifest.json').read_text())
    assert checked(src,'input_manifest.json','INPUT_COMPLETE')
    for name in ['x.bin','i.bin','p.bin','genes.csv','cells_fit.csv']:assert sha(src/name)==im['fitting_files'][name]
    X=sp.csr_matrix((np.fromfile(src/'x.bin',dtype='<f8'),np.fromfile(src/'i.bin',dtype='<i4'),np.fromfile(src/'p.bin',dtype='<i4')),
        shape=(im['n_cells'],im['n_genes']))
    assert np.isfinite(X.data).all() and (X.data>=0).all() and np.allclose(X.data,np.rint(X.data))
    X.eliminate_zeros();totals=np.asarray(X.sum(1)).ravel();assert (totals>0).all()
    # Normalize the original eligible RNA universe, before model gene intersection.
    X.data=np.log1p(X.data*np.repeat(10000/totals,np.diff(X.indptr)))
    symbols=pd.read_csv(prep/'Seurat_gene_names.csv',dtype=str,keep_default_na=False)
    input_genes=pd.read_csv(src/'genes.csv',dtype=str,keep_default_na=False).gene
    assert list(input_genes)==list(symbols.source) and symbols.Seurat.is_unique
    cells=pd.read_csv(src/'cells_fit.csv',dtype=str,keep_default_na=False).cell_id
    assert cells.is_unique and len(cells)==im['n_cells']
    pos={g:i for i,g in enumerate(symbols.Seurat)}
    features=(prep/'DL_features.txt').read_text().splitlines()
    ref=np.memmap(prep/'DL.float32.bin',dtype='<f4',mode='r',shape=(len(cells),len(features)))
    got=X[:,[pos[g] for g in features]].toarray().astype(np.float32)
    difference=float(np.max(np.abs(got-ref)))
    np.testing.assert_allclose(got,ref,rtol=1e-6,atol=1e-6)
    allowed=set((MODEL/'human/statistics/Brain_genes.txt').read_text().splitlines())
    keep=symbols.Seurat.isin(allowed).to_numpy();assert int(keep.sum())>1000
    # Match names exactly to the frozen model vocabulary; do not invent aliases or truth mappings.
    shared=symbols.Seurat[keep]
    pd.DataFrame(X[:,keep].T.toarray(),index=shared,columns=cells).to_csv(dest/'lognorm_model_input.csv')
    pd.DataFrame({'cell_id':cells}).to_csv(dest/'cells.csv',index=False)
    model_sources=[MODEL/'human/models/human-Brain.pt',MODEL/'human/statistics/Brain_genes.txt',
        MODEL/'human/statistics/Brain_cell_type.txt',MODEL/'human/graphs/human_Brain_data.npz',MODEL/'celltype2subtype.xlsx']
    write_json(dest/'input_manifest.json',dict(sample=sample,n_cells=len(cells),n_original_RNA_genes=X.shape[1],n_shared_genes=int(keep.sum()),
        normalization='log1p(count / original eligible RNA library total * 10000) before model gene intersection',
        normalized_exactly_once=True,model_input_dtype='float64',model_input_integers_cast=False,
        R_parity_n_cells=len(cells),R_parity_n_genes=len(features),R_parity_max_abs_difference=difference,
        input_manifest_sha256=sha(src/'input_manifest.json'),R_expression_sha256=sha(prep/'DL.float32.bin'),
        input_csv_sha256=sha(dest/'lognorm_model_input.csv'),cells_sha256=sha(dest/'cells.csv'),
        pretrained_sources={str(p):sha(p) for p in model_sources},
        source_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],prepared_at=utc(),no_truth_labels_read=True))
    del X,got,ref
    subprocess.run([str(OUT/'vendor_envs/scDeepSort/bin/python'),str(CODE/'scdeepsort_native_predict.py'),str(dest)],check=True)
    # Prediction is complete before loading any author truth.
    pm=json.loads((dest/'prediction_manifest.json').read_text());assert sha(dest/'predictions.csv.gz')==pm['predictions_sha256']
    pred=pd.read_csv(dest/'predictions.csv.gz',dtype=str,keep_default_na=False)
    truth=pd.read_csv(OUT/'evaluation_inputs'/sample/'truth.csv.gz',dtype=str,keep_default_na=False)
    old=pd.read_csv(OUT/'comparators/scDeepSort'/sample/'predictions.csv.gz',dtype=str,keep_default_na=False)
    assert list(pred.cell_id)==list(truth.cell_id)==list(old.cell_id)==list(cells)
    metrics=[];perclass=[];mapped=truth[['cell_id','L1']].copy()
    for condition,raw in [('historical_raw_counts',old),('corrected_LogNormalize',pred)]:
        p=map_labels('scDeepSort','published_human_Brain_atlas',raw.default)
        y=truth.L1.to_numpy();supported=np.isin(p,L1);unknown=raw.default.eq('Unknown').to_numpy()
        labels=sorted(set(y));precision,recall,f1,support=precision_recall_fscore_support(y,p,labels=L1,zero_division=0)
        metrics.append(dict(sample=sample,condition=condition,n_cells=len(y),accuracy=float((p==y).mean()),
            macroF1_present=float(f1[support>0].mean()),weightedF1=float(f1_score(y,p,labels=labels,average='weighted',zero_division=0)),
            coverage=float((~unknown).mean()),mapped_coverage=float(supported.mean()),unknown_rate=float(unknown.mean()),
            off_vocabulary_rate=float((~unknown&~supported).mean())))
        perclass.extend(dict(sample=sample,condition=condition,label=lab,precision=float(pr),recall=float(re),F1=float(ff),support=int(su))
            for lab,pr,re,ff,su in zip(L1,precision,recall,f1,support))
        mapped[condition]=p
        all_labels=L1+sorted(set(p)-set(L1))
        cf=pd.DataFrame(confusion_matrix(y,p,labels=all_labels),index=all_labels,columns=all_labels)
        cf.rename_axis('author_L1').to_csv(dest/f'{condition}_confusion.csv')
    pd.DataFrame(metrics).to_csv(dest/'metrics.csv',index=False);pd.DataFrame(perclass).to_csv(dest/'per_class.csv',index=False)
    mapped.to_csv(dest/'mapped_labels.csv.gz',index=False)
    changed=int((old.default!=pred.default).sum())
    record=dict(status='completed',sample=sample,n_cells=len(cells),input_R_parity_max_abs_difference=difference,
        changed_native_predictions_vs_raw_counts=changed,new_training=False,prediction_only=True,
        annotation_endpoint='official scDeepSort default; unsure_rate2 / 21 labels; not DG0.90',
        vocabulary_limit='Normal Brain model lacks Malignant; generic/fetal neurons remain ambiguous under frozen L1 mapping',
        wall_seconds=time.monotonic()-started,job=os.environ['SLURM_JOB_ID'],completed_at=utc(),
        old_predictions_untouched_sha256=sha(OUT/'comparators/scDeepSort'/sample/'predictions.csv.gz'),
        files={p.name:sha(p) for p in dest.iterdir() if p.is_file() and p.name not in ['manifest.json','COMPLETE']})
    write_json(dest/'manifest.json',record);complete(dest)
    print('CORRECTED_GBM_COMPLETE',json.dumps(record),flush=True)
    return record

def run():
    require_slurm();DEST.mkdir(parents=True,exist_ok=True)
    samples=sys.argv[1:] or ['TKU4163','NL022']
    result=[run_sample(s) for s in samples]
    write_json(DEST/'pilot_manifest.json',dict(status='pilot_completed',samples=result,full_121_complete=False,
        no_PTC_or_public_fits=True,job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    write_json(DEST/'status.json',dict(stage='C_GBM_SCDEEPSORT',status='pilot_completed',updated_at=utc(),jobs=[os.environ['SLURM_JOB_ID']],
        completed=[s+' official-input prediction and full-cell evaluation' for s in samples],
        remaining=['Root decision on measured full121 replay schedule','Refresh full comparisons only after full corrected predictions'],
        evidence=[str(DEST/'pilot_manifest.json')]))

if __name__=='__main__':run()
