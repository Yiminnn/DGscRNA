"""Independent confusion-count parity for the report's prediction reader/mapping."""
import os
from common import OUT,L1,require_slurm,write_json,sha,utc

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    from prediction_helpers import read_native,map_labels
    sample='TKU4163';truth=pd.read_csv(OUT/'evaluation_inputs'/sample/'truth.csv.gz',dtype=str)
    proofs=[]
    for method in ['DG-scRNA','scType','scCATCH','SCINA','SingleR','scDeepSort']:
        path=OUT/'GBM'/sample/'hvg2000/evaluation' if method=='DG-scRNA' else OUT/'comparators'/method/sample/'evaluation'
        metrics=pd.read_csv(path/'metrics.csv',dtype={'cutoff':str})
        saved=pd.read_csv(path/'confusions.csv.gz',dtype={'cutoff':str})
        if method=='DG-scRNA':
            metrics=metrics[metrics.stage=='terminal090'];saved=saved[saved.stage=='terminal090']
        if method in ['DG-scRNA','scType','scCATCH']:
            metrics=metrics[(metrics.budget=='hvg2000')&(metrics.route=='UMAP2_HDBSCAN_R')]
        if method in ['DG-scRNA','scType','scCATCH','SCINA']:
            metrics=metrics[metrics.library.isin(['CM2_glioma_other','CM2_primary_all_context'])]
        for c in metrics.itertuples():
            pred,source=read_native(method,sample,c.budget,c.route,c.library,c.cutoff)
            assert np.array_equal(pred.cell_id,truth.cell_id)
            p=map_labels(method,c.library,pred.prediction)
            got=pd.DataFrame({'truth':truth.L1,'prediction':p}).value_counts().sort_index()
            ref=saved.copy()
            for key in ['budget','route','library','cutoff']:ref=ref[ref[key]==getattr(c,key)]
            ref=ref.set_index(['truth','prediction']).n.sort_index()
            assert list(ref.index)==list(got.index) and np.array_equal(ref.to_numpy(),got.to_numpy()),(method,c.library,c.cutoff)
            proofs.append(dict(method=method,library=c.library,cutoff=c.cutoff,n_cells=len(p),source_sha256=sha(source)))
    write_json(OUT/'verification/selected_prediction_mapping_parity.json',dict(status='passed',sample=sample,
        exact_all_cell_confusion_counts=True,conditions=proofs,source_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],completed_at=utc()))

if __name__=='__main__':run()
