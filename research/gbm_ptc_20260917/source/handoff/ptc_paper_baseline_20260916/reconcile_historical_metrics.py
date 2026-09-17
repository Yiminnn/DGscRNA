"""Reproduce source-defined historical metrics without changing native labels."""
from pathlib import Path
import os,json,hashlib
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
BASE=ROOT/'results/hvg_ptc_20260916_v1'
OUT=BASE/'ptc_paper_baseline'

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import pandas as pd
    import numpy as np
    from sklearn.metrics import f1_score,accuracy_score,roc_auc_score,confusion_matrix
    raw=BASE/'ptc_recovery/archive/tcr/rawdata'
    d=pd.read_csv(raw/'data_with_validation.csv',keep_default_na=False)
    names=pd.read_csv(raw/'metadata.txt',sep='\t')
    samplemap=dict(zip(names.sc_ID,names.Sample))
    d['sample']=d['sample.name'].map(samplemap)
    d['cell_id']=d['sample']+'_'+d.iloc[:,0].str.extract(r'([ACGT]{12,})',expand=False)
    d=d.set_index('cell_id')
    d['scope']=np.select([d['sample'].str.startswith('N-'),d['sample'].str.startswith('MT-')],
                         ['Normal adjacent','Metastasis'],default='Primary tumor')
    T_names=json.loads((OUT/'historical_T_names_from_vignette.json').read_text())
    y=pd.to_numeric(d.validation_t_cell).astype(int)
    # Fixed original-vignette collapse; NK belongs to this historical broad list.
    predictions={'DG':d.Final_DGCyTOF_General.isin(T_names).astype(int)}
    for method in ['SCINA','scCATCH','scType']:
        predictions[method]=pd.to_numeric(d[method+'_ct_T_cells']).astype(int)
    groups=[('Overall',d.index)]+list(d.groupby('scope').groups.items())+list(d.groupby('sample').groups.items())
    rows=[]
    for method,pred in predictions.items():
        for scope,ids in groups:
            yy=y.loc[ids];pp=pred.loc[ids]
            tn,fp,fn,tp=confusion_matrix(yy,pp,labels=[0,1]).ravel()
            rows.append(dict(method=method,scope=scope,n=len(ids),TP=int(tp),TN=int(tn),FP=int(fp),FN=int(fn),
                F1_source_default_positive0=f1_score(yy,pp,pos_label=0),
                F1_T_positive1=f1_score(yy,pp,pos_label=1),
                AUC_from_binary_calls=roc_auc_score(yy,pp),accuracy=accuracy_score(yy,pp)))
    metrics=pd.DataFrame(rows)
    metrics.to_csv(OUT/'historical_original_validation_metrics.csv',index=False)
    targets=pd.read_csv(OUT/'paper_Table2_target_comparison.csv')
    columns={'F1 score':'F1_source_default_positive0','AUC-ROC':'AUC_from_binary_calls','Accuracy':'accuracy'}
    comparison=[]
    for r in targets.itertuples(index=False):
        value=metrics.loc[(metrics.method==r.method)&(metrics.scope==r.scope),columns[r.metric]].item()
        comparison.append(dict(method=r.method,scope=r.scope,metric=r.metric,paper=r.paper,
            reconstructed=value,agrees_at_reported_4_decimals=round(value,4)==r.paper,
            truth='Original CSV validation_t_cell (independently matched to any filtered TCR contig)',
            definition=columns[r.metric]))
    comparison=pd.DataFrame(comparison)
    comparison.to_csv(OUT/'paper_Table2_historical_source_comparison.csv',index=False)
    print(comparison.to_string(index=False),flush=True)
    # Small explicit export permits independent verification in the original R metric library.
    pd.DataFrame({'cell_id':d.index,'scope':d.scope,'truth':y,'prediction':predictions['DG']}).to_csv(
        OUT/'original_DG_binary_pairs_for_R.csv.gz',index=False)
    # Branch provenance: preserve both raw general labels and the fixed dictionary
    # used by the archived notebook, instead of evaluating literal capitalization.
    extended=pd.read_csv(raw/'data_with_validation+3_cell_types 2.csv',keep_default_na=False)
    extended['sample']=extended['sample.name'].map(samplemap)
    extended['cell_id']=extended['sample']+'_'+extended['X'].str.extract(r'([ACGT]{12,})',expand=False)
    extended=extended.set_index('cell_id').loc[d.index]
    native=d.Final_DGCyTOF
    pub_general=extended['hdbscan.UMAP_clusters_NCOMMREFF_mean_DGCyTOF_General']
    # For NCOMMREFF the original Python parser discards exactly this prefix.
    native_general=native.str.replace(r'^NCOMMREFF\+','',regex=True)
    mask=d.scope.eq('Primary tumor')
    equal=native_general[mask].eq(pub_general[mask])
    branch=pd.DataFrame({'cell_id':d.index[mask],'historical_native':native[mask],
        'source_parsed_general':native_general[mask],'archived_Pubmed_branch_general':pub_general[mask],
        'agrees':equal})
    branch.to_csv(OUT/'TTU_original_branch_provenance.csv.gz',index=False)
    manifest=dict(job=os.environ['SLURM_JOB_ID'],paper_rows=len(comparison),
        matching_rows=int(comparison.agrees_at_reported_4_decimals.sum()),
        matching_DG_rows=int(comparison.loc[comparison.method.eq('DG'),'agrees_at_reported_4_decimals'].sum()),
        TTU_branch_cells=len(branch),TTU_branch_exact=int(equal.sum()),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        status='historical_metric_reconciliation; full workflow gate remains pending')
    (OUT/'historical_metric_reconciliation_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2),flush=True)

if __name__=='__main__':run()
