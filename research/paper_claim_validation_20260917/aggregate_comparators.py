"""Retrospective patient-heldout matched-marker and separate-reference comparisons."""
import json
import os
from common import OUT, require_slurm, checked, complete, sha, write_json, utc
METHODS=['scType','scCATCH','SCINA','SingleR','scDeepSort']

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    from scipy.stats import wilcoxon
    assert checked(OUT/'summary','aggregate_manifest.json','AGGREGATE_COMPLETE')
    cohort=pd.read_csv(OUT/'protocol/cohort.csv');frames=[];manifests={}
    for method in METHODS:
        for sample in cohort['sample']:
            src=OUT/'comparators'/method/sample/'evaluation'
            assert checked(src),str(src)
            frames.append(pd.read_csv(src/'metrics.csv',dtype={'cutoff':str}))
            manifests[str(src.relative_to(OUT))]=sha(src/'manifest.json')
    tools=pd.concat(frames,ignore_index=True)
    dg=pd.read_csv(OUT/'summary/all_annotation_metrics.csv.gz',dtype={'cutoff':str})
    dg=dg[(dg.stage=='terminal090')&(dg.family=='native_R_budget')].copy();dg['method']='DG-scRNA'
    shared=['method','sample','patient','primary','budget','route','library','cutoff','macroF1_present','macroF1_fixed11','accuracy','coverage','unknown_rate','mapped_coverage','off_vocabulary_rate']
    allrows=pd.concat([dg[shared],tools[shared]],ignore_index=True)
    dest=OUT/'comparison_summary';dest.mkdir(exist_ok=True)
    allrows.to_csv(dest/'all_comparator_metrics.csv.gz',index=False,compression='gzip')
    fixed_geometry=(allrows.budget=='hvg2000')&(allrows.route=='UMAP2_HDBSCAN_R')
    keep=fixed_geometry | allrows.method.isin(['SCINA','SingleR','scDeepSort'])
    keep &= ~((allrows.method=='SingleR')&(allrows.cutoff!='default'))
    matched=allrows[keep & ~allrows.library.isin(['CARE_TME','BrainAtlas112','UNION_all'])].copy()
    folds=pd.read_csv(OUT/'protocol/patient_folds.csv')[['patient','fold']].drop_duplicates()
    spec=['method','budget','route','library','cutoff']
    measures=['macroF1_present','macroF1_fixed11','accuracy','coverage','unknown_rate','mapped_coverage','off_vocabulary_rate']
    selected=[];predictions=[];fixed=[]
    for name,c in [('primary97',matched[matched.primary]),('all121',matched)]:
        patient=c.groupby(['patient']+spec)[measures].mean().reset_index().merge(folds,on='patient',validate='many_to_one')
        for method,g in patient.groupby('method'):
            assert g.patient.nunique()==(55 if name=='primary97' else 59)
            for fold in sorted(g.fold.unique()):
                train=g[g.fold!=fold];test=g[g.fold==fold]
                rank=train.groupby(spec)[measures].mean().reset_index().sort_values(['macroF1_present']+spec,ascending=[False]+[True]*len(spec),kind='stable')
                best=rank.iloc[0];mask=np.ones(len(test),dtype=bool)
                for key in spec:mask &= test[key].eq(best[key]).to_numpy()
                held=test[mask];assert held.patient.nunique()==test.patient.nunique()
                selected.append(dict(cohort=name,fold=int(fold),n_training_patients=train.patient.nunique(),
                    n_test_patients=held.patient.nunique(),training_mean_F1=float(best.macroF1_present),
                    n_candidate_configs=len(rank),**{key:best[key] for key in spec}))
                for row in held.to_dict('records'):predictions.append(dict(cohort=name,**row))
            anchor=g[(g.library=='CM2_glioma_other')&(g.cutoff==('mean' if method=='DG-scRNA' else 'default' if method in ['scType','scCATCH'] else 'overlap1'))]
            for row in anchor.to_dict('records'):fixed.append(dict(cohort=name,**row))
    cv=pd.DataFrame(predictions);cv.to_csv(dest/'patient_heldout_results.csv',index=False)
    pd.DataFrame(selected).to_csv(dest/'training_patient_choices.csv',index=False)
    pd.DataFrame(fixed).to_csv(dest/'fixed_glioma_marker_anchors.csv',index=False)
    summary=cv.groupby(['cohort','method'])[measures].agg(['mean','std','count']);summary.columns=['_'.join(v) for v in summary.columns]
    summary.reset_index().to_csv(dest/'patient_heldout_summary.csv',index=False)
    rng=np.random.default_rng(20260917);contrasts=[];differences=[]
    for name,c in cv.groupby('cohort'):
        baseline=c[c.method=='DG-scRNA'][['patient','macroF1_present','coverage']].rename(columns={'macroF1_present':'DG_F1','coverage':'DG_coverage'})
        for method,g in c[c.method!='DG-scRNA'].groupby('method'):
            paired=g.merge(baseline,on='patient',validate='one_to_one')
            delta=(paired.macroF1_present-paired.DG_F1).to_numpy()
            draws=delta[rng.integers(0,len(delta),size=(10000,len(delta)))].mean(axis=1)
            contrasts.append(dict(cohort=name,method=method,comparison='competitor minus DG-scRNA; matched fixed partition for cluster-based marker methods',
                n_patients=len(delta),mean_delta=float(delta.mean()),CI95_low=float(np.quantile(draws,.025)),CI95_high=float(np.quantile(draws,.975)),
                p=float(wilcoxon(delta).pvalue) if np.any(np.abs(delta)>1e-14) else 1.,
                mean_coverage_delta=float((paired.coverage-paired.DG_coverage).mean()),
                information='labelled training-patient reference' if method=='SingleR' else 'published atlas pretrained GNN' if method=='scDeepSort' else 'matched marker libraries'))
            for row,diff in zip(paired.itertuples(),delta):differences.append(dict(cohort=name,method=method,patient=row.patient,delta=float(diff)))
    stats=pd.DataFrame(contrasts);stats['p_Holm']=1.
    for name,g in stats.groupby('cohort'):
        ix=g.sort_values('p').index
        stats.loc[ix,'p_Holm']=np.minimum(1,np.maximum.accumulate(stats.loc[ix,'p'].to_numpy()*np.arange(len(ix),0,-1)))
    stats.to_csv(dest/'paired_patient_comparisons.csv',index=False)
    pd.DataFrame(differences).to_csv(dest/'paired_patient_differences.csv',index=False)
    write_json(dest/'manifest.json',dict(status='completed',methods=['DG-scRNA']+METHODS,
        primary_scope='Original native-R HVG2000/UMAP2/HDBSCAN partition shared by cluster-based marker methods. Marker/cutoff selected on training patients; SingleR default pruned and scDeepSort published default use distinct reference information.',
        not_compared_as_equal='DG/scType24-configuration tuning is not contrasted with fixed-partition scCATCH as an equal tuning comparison.',
        uncertainty='Paired patient bootstrap conditional on frozen cross-validated predictions; model selection not re-fitted inside bootstrap.',
        retrospective=True,all_cells_in_denominator=True,source_manifests=manifests,
        files={p.name:sha(p) for p in dest.iterdir() if p.suffix in ['.csv','.gz']},
        job=os.environ['SLURM_JOB_ID'],source_sha256=sha(__file__),completed_at=utc()))
    complete(dest)

if __name__=='__main__':run()
