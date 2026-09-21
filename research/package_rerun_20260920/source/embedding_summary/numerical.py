"""Unchanged A1 patient selection/statistics; two terminal endpoints only."""
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
MEASURES=['lfine_macroF1', 'coverage', 'legacy_called_coverage', 'abstain_rate', 'offvocab_rate', 'acc_on_called', 'lfine_scored_class_cell_fraction', 'marker_vocab_oracle_upper_bound']

def select(allmetrics, anchors, eligibility, folds, spec, protocol):
    foldmap=folds.drop_duplicates("patient").set_index("patient").fold.to_dict()
    selected_frames=[];choices=[];patients_all=[];paired=[];rng=np.random.default_rng(protocol['bootstrap_seed'])
    for cohort in ['primary97','all121']:
        data=allmetrics[allmetrics.primary] if cohort=='primary97' else allmetrics
        n_samples=97 if cohort=='primary97' else 121;n_patients=55 if cohort=='primary97' else 59
        assert data['sample'].nunique()==n_samples and data.patient.nunique()==n_patients
        terminal=data[data.stage.eq('terminal090')]
        candidates=terminal.groupby(['patient','budget','space','method','k'])[MEASURES].mean().reset_index()
        candidates['fold']=candidates.patient.map(foldmap)
        selected=[]
        for (budget,space,method),g in candidates.groupby(['budget','space','method']):
            for fold in sorted(g.fold.unique()):
                train=g[g.fold.ne(fold)];test=g[g.fold.eq(fold)]
                rankings=train.groupby('k').lfine_macroF1.mean().reset_index().sort_values(['lfine_macroF1','k'],ascending=[False,True],kind='stable')
                expected={0} if method=='HDBSCAN_R' else set(spec['K'])
                assert set(rankings.k)==expected and rankings.lfine_macroF1.notna().all(),'No candidate may be dropped from Lfine K selection'
                eligible_sets=[set(q.loc[q.lfine_macroF1.notna(),'patient']) for _,q in train.groupby('k')]
                assert all(ids==eligible_sets[0] for ids in eligible_sets),'Candidate-dependent patient missingness blocks K selection'
                best=int(rankings.iloc[0].k)
                held=test[test.k.eq(best)];assert held.patient.nunique()==test.patient.nunique()
                train_ids=set(train.patient);held_ids=set(held.patient);assert not train_ids&held_ids
                choices.append(dict(cohort=cohort,budget=budget,space=space,space_display='ICA2 adaptive' if space=='ICA2' else space,method=method,fold=int(fold),selected_k=best,
                    selection_metric='lfine_macroF1',selection_endpoint='terminal090',training_patient_mean_lfine_macroF1=float(rankings.iloc[0].lfine_macroF1),
                    n_training_patients=len(train_ids),n_training_patients_Lfine_eligible=len(eligible_sets[0]),n_test_patients=len(held_ids),
                    n_test_patients_Lfine_eligible=int(held.lfine_macroF1.notna().sum()),exact_tie_rule='smaller_K',reused_L1_K=False))
                take=data[data.budget.eq(budget)&data.space.eq(space)&data.method.eq(method)&data.k.eq(best)&data.patient.isin(held_ids)].copy()
                take['fold']=int(fold);take['cohort']=cohort;selected.append(take)
        chosen=pd.concat(selected,ignore_index=True)
        assert len(chosen)==n_samples*2*7*3*2 and not chosen.duplicated(['sample','budget','space','method','stage']).any()
        selected_frames.append(chosen)
        patientrows=[]
        for keys,g in chosen.groupby(['cohort','patient','budget','space','method','stage']):
            row=dict(zip(['cohort','patient','budget','space','method','stage'],keys));row['n_samples_total']=len(g)
            row['space_display']='ICA2 adaptive' if row['space']=='ICA2' else row['space']
            for metric in MEASURES:
                row[metric]=g[metric].mean();row[metric+'_n_samples']=int(g[metric].notna().sum())
            patientrows.append(row)
        patients=pd.DataFrame(patientrows);patients_all.append(patients)
        for budget in spec['budgets']:
            ref=anchors[anchors.budget.eq(budget)&anchors.stage.eq('terminal090')&(anchors.primary if cohort=='primary97' else True)]
            ref=ref.groupby('patient')[MEASURES].mean().reset_index()
            current=patients[patients.budget.eq(budget)&patients.stage.eq('terminal090')]
            for (space,method),g in current.groupby(['space','method']):
                join=g.merge(ref[['patient','lfine_macroF1','coverage']],on='patient',validate='one_to_one',suffixes=('','_anchor')).sort_values('patient')
                assert len(join)==n_patients
                assert join.lfine_macroF1.notna().equals(join.lfine_macroF1_anchor.notna()),'Candidate and anchor must have the same truth-only eligibility'
                valid=join[join.lfine_macroF1.notna()];assert len(valid)>=2
                delta=(valid.lfine_macroF1-valid.lfine_macroF1_anchor).to_numpy()
                boot=delta[rng.integers(0,len(delta),size=(protocol['bootstrap_replicates'],len(delta)))].mean(1)
                paired.append(dict(cohort=cohort,budget=budget,space=space,space_display='ICA2 adaptive' if space=='ICA2' else space,method=method,n_patients_total=n_patients,n_patients_Lfine=len(valid),
                    n_patients_Lfine_ineligible=n_patients-len(valid),n_samples_total=n_samples,
                    n_samples_Lfine=sum(e['lfine_metric_eligible'] for e in eligibility.values() if cohort=='all121' or e['primary']),
                    candidate_mean_F1=valid.lfine_macroF1.mean(),anchor_mean_F1=valid.lfine_macroF1_anchor.mean(),mean_delta=delta.mean(),
                    CI95_low=np.quantile(boot,.025),CI95_high=np.quantile(boot,.975),
                    p_wilcoxon=wilcoxon(delta).pvalue if (np.abs(delta)>1e-14).any() else 1.,
                    candidate_coverage_all_patients=join.coverage.mean(),anchor_coverage_all_patients=join.coverage_anchor.mean(),
                    n_patients_coverage=int(join.coverage.notna().sum()),selection_metric='lfine_macroF1'))
    stats=pd.DataFrame(paired);stats['p_Holm']=1.
    for _,g in stats.groupby(['cohort','budget']):
        assert len(g)==21
        order=g.sort_values('p_wilcoxon',kind='stable').index
        stats.loc[order,'p_Holm']=np.minimum(1.,np.maximum.accumulate(stats.loc[order,'p_wilcoxon'].to_numpy()*np.arange(21,0,-1)))
    assert len(choices)==420 and len(stats)==84
    return pd.DataFrame(choices), pd.concat(selected_frames), pd.concat(patients_all), stats

def independently_select(data, anchors, eligibility, folds, spec, protocol):
    foldmap=folds.drop_duplicates("patient").set_index("patient").fold.to_dict()
    choices=[];selected_frames=[];patient_frames=[];comparisons=[];rng=np.random.default_rng(protocol['bootstrap_seed'])
    for cohort in ['primary97','all121']:
        part=data[data.primary] if cohort=='primary97' else data
        n_samples=97 if cohort=='primary97' else 121;n_patients=55 if cohort=='primary97' else 59
        assert part['sample'].nunique()==n_samples and part.patient.nunique()==n_patients
        terminal=part[part.stage.eq('terminal090')]
        means=terminal.groupby(['patient','budget','space','method','k']).lfine_macroF1.mean().reset_index()
        means['fold']=means.patient.map(foldmap)
        these=[]
        for budget,space,method in sorted(set(zip(means.budget,means.space,means.method))):
            unit=means[means.budget.eq(budget)&means.space.eq(space)&means.method.eq(method)]
            for fold in sorted(unit.fold.unique()):
                training=unit[unit.fold.ne(fold)];held=unit[unit.fold.eq(fold)]
                scores={int(k):float(score) for k,score in training.groupby('k').lfine_macroF1.mean().items()}
                assert set(scores)==({0} if method=='HDBSCAN_R' else set(spec['K'])) and all(np.isfinite(v) for v in scores.values())
                k=min(scores,key=lambda x:(-scores[x],x))
                valid_sets=[set(g.loc[g.lfine_macroF1.notna(),'patient']) for _,g in training.groupby('k')]
                assert all(v==valid_sets[0] for v in valid_sets)
                selected_held=held[held.k.eq(k)];train_ids=set(training.patient);test_ids=set(selected_held.patient)
                assert not train_ids&test_ids and len(train_ids|test_ids)==n_patients
                row=dict(cohort=cohort,budget=budget,space=space,space_display='ICA2 adaptive' if space=='ICA2' else space,method=method,fold=int(fold),selected_k=k,
                    selection_metric='lfine_macroF1',selection_endpoint='terminal090',training_patient_mean_lfine_macroF1=scores[k],
                    n_training_patients=len(train_ids),n_training_patients_Lfine_eligible=len(valid_sets[0]),n_test_patients=len(test_ids),
                    n_test_patients_Lfine_eligible=int(selected_held.lfine_macroF1.notna().sum()),exact_tie_rule='smaller_K',reused_L1_K=False)
                choices.append(row);these.append(row)
        membership=part.copy();membership['fold']=membership.patient.map(foldmap);membership['cohort']=cohort
        chosen=pd.DataFrame(these)[['budget','space','method','fold','selected_k']].rename(columns={'selected_k':'k'})
        selected=membership.merge(chosen,on=['budget','space','method','fold','k'],how='inner',validate='many_to_one')
        assert len(selected)==n_samples*84 and not selected.duplicated(['sample','budget','space','method','stage']).any()
        selected_frames.append(selected)
        grouping=['cohort','patient','budget','space','method','stage']
        patients=selected.groupby(grouping)[MEASURES].mean()
        patients['n_samples_total']=selected.groupby(grouping).size()
        for metric in MEASURES:patients[metric+'_n_samples']=selected.groupby(grouping)[metric].count()
        patients=patients.reset_index();patients['space_display']=patients.space.map(lambda s:'ICA2 adaptive' if s=='ICA2' else s);patient_frames.append(patients)
        for budget in spec['budgets']:
            ref=anchors[anchors.budget.eq(budget)&anchors.stage.eq('terminal090')&(anchors.primary if cohort=='primary97' else True)]
            ref=ref.groupby('patient')[['lfine_macroF1','coverage']].mean()
            current=patients[patients.budget.eq(budget)&patients.stage.eq('terminal090')]
            for space,method in sorted(set(zip(current.space,current.method))):
                g=current[current.space.eq(space)&current.method.eq(method)].set_index('patient').sort_index()
                assert set(g.index)==set(ref.index) and len(g)==n_patients
                ref2=ref.loc[g.index];valid=g.lfine_macroF1.notna();assert valid.equals(ref2.lfine_macroF1.notna())
                delta=(g.loc[valid,'lfine_macroF1']-ref2.loc[valid,'lfine_macroF1']).to_numpy();assert len(delta)>=2
                samples=rng.integers(len(delta),size=(protocol['bootstrap_replicates'],len(delta)))
                distribution=np.mean(delta[samples],axis=1)
                comparisons.append(dict(cohort=cohort,budget=budget,space=space,space_display='ICA2 adaptive' if space=='ICA2' else space,method=method,n_patients_total=n_patients,n_patients_Lfine=len(delta),
                    n_patients_Lfine_ineligible=n_patients-len(delta),n_samples_total=n_samples,
                    n_samples_Lfine=sum(e['lfine_metric_eligible'] for e in eligibility.values() if cohort=='all121' or e['primary']),
                    candidate_mean_F1=g.loc[valid,'lfine_macroF1'].mean(),anchor_mean_F1=ref2.loc[valid,'lfine_macroF1'].mean(),mean_delta=delta.mean(),
                    CI95_low=np.quantile(distribution,.025),CI95_high=np.quantile(distribution,.975),
                    p_wilcoxon=wilcoxon(delta).pvalue if (np.abs(delta)>1e-14).any() else 1.,
                    candidate_coverage_all_patients=g.coverage.mean(),anchor_coverage_all_patients=ref2.coverage.mean(),
                    n_patients_coverage=int(g.coverage.notna().sum()),selection_metric='lfine_macroF1'))
    stats=pd.DataFrame(comparisons);stats['p_Holm']=1.
    for _,group in stats.groupby(['cohort','budget']):
        assert len(group)==21;running=0.
        for rank,index in enumerate(sorted(group.index,key=lambda i:stats.at[i,'p_wilcoxon'])):
            running=max(running,(21-rank)*stats.at[index,'p_wilcoxon']);stats.at[index,'p_Holm']=min(1.,running)
    return pd.DataFrame(choices), pd.concat(selected_frames), pd.concat(patient_frames), stats

def equal(a,b,keys,columns=None):
    import numpy as np
    import pandas as pd
    assert not a.duplicated(keys).any() and not b.duplicated(keys).any()
    a=a.set_index(keys).sort_index();b=b.set_index(keys).sort_index();assert a.index.equals(b.index)
    for key in list(a.columns) if columns is None else columns:
        if pd.api.types.is_numeric_dtype(a[key]) and not pd.api.types.is_bool_dtype(a[key]):
            np.testing.assert_allclose(a[key].to_numpy(float),b[key].to_numpy(float),rtol=0,atol=1e-12,equal_nan=True,err_msg=key)
        else:assert a[key].fillna('<NA>').astype(str).equals(b[key].fillna('<NA>').astype(str)),key
