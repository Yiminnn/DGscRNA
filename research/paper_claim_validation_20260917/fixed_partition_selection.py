"""Freeze DG marker/cutoff choices independently of unfinished comparator fits."""
import os
from common import OUT,require_slurm,checked,sha,write_json,complete,utc

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    assert checked(OUT/'summary','aggregate_manifest.json','AGGREGATE_COMPLETE')
    dest=OUT/'DG_fixed_partition_selection';dest.mkdir(exist_ok=True)
    if checked(dest):return
    path=OUT/'summary/all_annotation_metrics.csv.gz'
    dg=pd.read_csv(path,dtype={'cutoff':str})
    dg=dg[(dg.stage=='terminal090')&(dg.family=='native_R_budget')&
          (dg.budget=='hvg2000')&(dg.route=='UMAP2_HDBSCAN_R')&
          ~dg.library.isin(['CARE_TME','BrainAtlas112','UNION_all'])].copy()
    dg['method']='DG-scRNA'
    folds=pd.read_csv(OUT/'protocol/patient_folds.csv')[['patient','fold']].drop_duplicates()
    spec=['method','budget','route','library','cutoff']
    measures=['macroF1_present','macroF1_fixed11','accuracy','coverage','unknown_rate','mapped_coverage','off_vocabulary_rate']
    assert not dg.duplicated(['sample']+spec).any()
    choices=[];predictions=[]
    for name,c in [('primary97',dg[dg.primary]),('all121',dg)]:
        patient=c.groupby(['patient']+spec)[measures].mean().reset_index().merge(folds,on='patient',validate='many_to_one')
        assert patient.patient.nunique()==(55 if name=='primary97' else 59)
        for fold in sorted(patient.fold.unique()):
            train=patient[patient.fold!=fold];test=patient[patient.fold==fold]
            rank=train.groupby(spec)[measures].mean().reset_index().sort_values(
                ['macroF1_present']+spec,ascending=[False]+[True]*len(spec),kind='stable')
            assert len(rank)==13*3
            best=rank.iloc[0];mask=np.ones(len(test),dtype=bool)
            for key in spec:mask &= test[key].eq(best[key]).to_numpy()
            held=test[mask];assert held.patient.nunique()==test.patient.nunique()
            assert set(train.patient).isdisjoint(held.patient)
            choices.append(dict(cohort=name,fold=int(fold),n_training_patients=train.patient.nunique(),
                n_test_patients=held.patient.nunique(),training_mean_F1=float(best.macroF1_present),
                n_candidate_configs=len(rank),**{key:best[key] for key in spec}))
            for row in held.to_dict('records'):predictions.append(dict(cohort=name,**row))
    pd.DataFrame(predictions).to_csv(dest/'patient_heldout_results.csv',index=False)
    pd.DataFrame(choices).to_csv(dest/'training_patient_choices.csv',index=False)
    write_json(dest/'manifest.json',dict(status='completed',n_heldout_patient_rows=len(predictions),
        source_scope='Same original fixed-partition DG selection as aggregate_comparators; only the scheduling dependency is separated.',
        unchanged_candidates_and_ties=True,no_new_fits=True,
        aggregate_manifest_sha256=sha(OUT/'summary/aggregate_manifest.json'),
        annotation_metrics_sha256=sha(path),folds_sha256=sha(OUT/'protocol/patient_folds.csv'),
        files={p.name:sha(p) for p in dest.glob('*.csv')},
        source_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest)
    print('DG_FIXED_PARTITION_SELECTION_COMPLETE',len(predictions),flush=True)

if __name__=='__main__':run()
