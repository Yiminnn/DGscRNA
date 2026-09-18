"""Patient-based native-R comparisons and frozen retrospective heldout selection."""
import argparse
import json
import os
from common import OUT, FEATURES, ROUTES, require_slurm, samples, checked, sha, write_json, complete, utc

def run(partial=False):
    require_slurm()
    import numpy as np
    import pandas as pd
    from scipy.stats import wilcoxon
    frames=[];partitions=[];missing=[]
    for sample in samples():
        for budget in FEATURES:
            d=OUT/'GBM'/sample/budget/'evaluation'
            if not checked(d):missing.append(sample+'/'+budget);continue
            frames.append(pd.read_csv(d/'metrics.csv',dtype={'cutoff':str}))
            partitions.append(pd.read_csv(d/'clustering.csv'))
    if missing and not partial:raise RuntimeError(f'{len(missing)} core units remain incomplete')
    assert frames,'No complete units to aggregate'
    df=pd.concat(frames,ignore_index=True);cp=pd.concat(partitions,ignore_index=True)
    dest=OUT/'summary';dest.mkdir(exist_ok=True)
    df.to_csv(dest/'all_annotation_metrics.csv.gz',index=False,compression='gzip')
    cp.to_csv(dest/'all_clustering_metrics.csv.gz',index=False,compression='gzip')
    spec=['budget','route','library','cutoff','stage','family','seed']
    measures=['macroF1_present','macroF1_fixed11','accuracy','coverage','unknown_rate','mapped_coverage','off_vocabulary_rate','coarse10_macroF1_present']
    assert not df.duplicated(['sample']+spec).any()
    native=df[(df.family=='native_R_budget')&(df.seed==42)].copy()
    folds=pd.read_csv(OUT/'protocol/patient_folds.csv')
    patient_folds=folds[['patient','fold']].drop_duplicates()
    assert not patient_folds.patient.duplicated().any()
    all_summary=[];paired=[];paired_diffs=[];cv=[];choices=[];factorial=[]
    rng=np.random.default_rng(20260917)
    config=['budget','route','library','cutoff']
    base=dict(budget='hvg2000',route='UMAP2_HDBSCAN_R',library='CM2_glioma_other',cutoff='mean')
    for cohort,t in [('primary97',native[native.primary]),('all121',native)]:
        patient=t.groupby(['patient']+spec,dropna=False)[measures].mean().reset_index()
        means=patient.groupby(spec,dropna=False)[measures].agg(['mean','std','count']).reset_index()
        means.columns=['_'.join(v for v in key if v) for key in means.columns]
        means.insert(0,'cohort',cohort);all_summary.append(means)
        final=patient[(patient.stage=='terminal090')].copy()
        fixed=final[(final.library==base['library'])&(final.cutoff==base['cutoff'])]
        baseline=fixed[(fixed.budget==base['budget'])&(fixed.route==base['route'])][['patient','macroF1_present','coverage']]
        baseline=baseline.rename(columns={'macroF1_present':'baseline_F1','coverage':'baseline_coverage'})
        for (budget,route),g in fixed.groupby(['budget','route']):
            a=g.merge(baseline,on='patient',validate='one_to_one')
            delta=(a.macroF1_present-a.baseline_F1).to_numpy()
            if not len(delta):continue
            boot=delta[rng.integers(0,len(delta),size=(10000,len(delta)))].mean(axis=1)
            p=float(wilcoxon(delta).pvalue) if np.any(np.abs(delta)>1e-14) else 1.
            paired.append(dict(cohort=cohort,budget=budget,route=route,n_patients=len(a),
                mean_delta=float(delta.mean()),median_delta=float(np.median(delta)),
                CI95_low=float(np.quantile(boot,.025)),CI95_high=float(np.quantile(boot,.975)),p_wilcoxon=p,
                mean_coverage_delta=float((a.coverage-a.baseline_coverage).mean()),
                baseline='hvg2000/UMAP2_HDBSCAN_R/CM2_glioma_other/mean',
                inference_unit='patient; samples averaged within patient; seeds never treated as patients'))
            for row,difference in zip(a.itertuples(),delta):paired_diffs.append(dict(cohort=cohort,budget=budget,route=route,patient=row.patient,delta=float(difference)))
        if missing:continue
        # The database-only analysis excludes libraries explicitly used to construct author labels.
        for evidence,allowed in [('database_or_external',~final.library.isin(['CARE_TME','BrainAtlas112','UNION_all'])),
                                 ('all_libraries_concordance',np.ones(len(final),dtype=bool))]:
            cand=final[allowed].merge(patient_folds,on='patient',validate='many_to_one')
            for fold in sorted(cand.fold.unique()):
                train=cand[cand.fold!=fold];test=cand[cand.fold==fold]
                objectives=[('all_routes',train)] + [(r,train[train.route==r]) for r in ROUTES]
                for selection,training in objectives:
                    ranking=training.groupby(config).macroF1_present.mean().reset_index()
                    ranking=ranking.sort_values(['macroF1_present']+config,ascending=[False]+[True]*len(config),kind='stable')
                    best=ranking.iloc[0]
                    mask=np.ones(len(test),dtype=bool)
                    for k in config:mask &= test[k].eq(best[k]).to_numpy()
                    held=test[mask]
                    assert held.patient.nunique()==test.patient.nunique()
                    choice=dict(cohort=cohort,evidence=evidence,fold=int(fold),selection=selection,
                        **{k:best[k] for k in config},training_mean_F1=float(best.macroF1_present),
                        n_training_patients=train.patient.nunique(),n_test_patients=test.patient.nunique())
                    choices.append(choice)
                    for row in held.to_dict('records'):
                        cv.append(dict(cohort=cohort,evidence=evidence,selection=selection,**row))
        # Matched 2x2: marker selection uses training terminal scores; the same selected
        # library is then used for both marker-only and DL cells of the factorial.
        f=patient[(patient.budget=='hvg2000')&(patient.route=='UMAP2_HDBSCAN_R')&(patient.cutoff=='mean')]
        f=f[~f.library.isin(['CARE_TME','BrainAtlas112','UNION_all'])].merge(patient_folds,on='patient',validate='many_to_one')
        for fold in sorted(f.fold.unique()):
            train=f[(f.fold!=fold)&(f.stage=='terminal090')]
            ranking=train.groupby('library').macroF1_present.mean().sort_values(ascending=False,kind='stable')
            chosen=ranking.index[0]
            for selection,lib in [('fixed_marker','CM2_glioma_other'),('training_patient_selected_marker',chosen)]:
                test=f[(f.fold==fold)&(f.library==lib)&f.stage.isin(['marker_only','terminal090'])]
                for row in test.to_dict('records'):factorial.append(dict(cohort=cohort,marker_selection=selection,**row))
    pd.concat(all_summary,ignore_index=True).to_csv(dest/'patient_aggregates.csv',index=False)
    stat=pd.DataFrame(paired)
    # Holm adjustment within each cohort's prespecified 23 non-baseline comparisons.
    stat['p_Holm']=1.
    for cohort in stat.cohort.unique():
        ix=stat.index[(stat.cohort==cohort)&~((stat.budget=='hvg2000')&(stat.route=='UMAP2_HDBSCAN_R'))]
        ordered=stat.loc[ix].sort_values('p_wilcoxon').index
        adj=np.maximum.accumulate(stat.loc[ordered,'p_wilcoxon'].to_numpy()*np.arange(len(ordered),0,-1))
        stat.loc[ordered,'p_Holm']=np.minimum(adj,1.)
    stat.to_csv(dest/'fixed_marker_patient_paired.csv',index=False)
    pd.DataFrame(paired_diffs).to_csv(dest/'fixed_marker_patient_differences.csv',index=False)
    pd.DataFrame(choices).to_csv(dest/'patient_heldout_selection.csv',index=False)
    pd.DataFrame(cv).to_csv(dest/'patient_heldout_test_results.csv',index=False)
    pd.DataFrame(factorial).to_csv(dest/'marker_DL_factorial_heldout.csv',index=False)
    df[df.stage=='terminal090'].groupby(['dl_status','training_executed']).size().rename('n').reset_index().to_csv(dest/'terminal_status_counts.csv',index=False)
    m=dict(status='partial' if missing else 'completed',expected_units=726,completed_units=len(frames),missing_units=missing,
        n_annotation_metric_rows=len(df),n_clustering_results=len(cp),
        primary='patient-weighted mean of per-sample L1 macro-F1 over present author classes; full-cell denominator',
        multiplicity='Holm across23 fixed-marker configurations versus original default separately for primary/full cohort',
        uncertainty='10000 paired patient bootstrap resamples; Wilcoxon paired sensitivity; no noninferiority claim',
        selection='5 folds frozen by patient ID; training labels choose config; heldout labels score only; retrospective cohort',
        scope='Within GSE274546 candidates. Does not itself establish superiority to other annotation methods.',
        job=os.environ['SLURM_JOB_ID'],source_sha256=sha(__file__),completed_at=utc())
    write_json(dest/'aggregate_manifest.json',m)
    if not missing:complete(dest,'aggregate_manifest.json','AGGREGATE_COMPLETE')
    print('AGGREGATED',len(frames),'missing',len(missing),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--partial',action='store_true');a=p.parse_args();run(a.partial)
