"""Independent full Lfine K selection, denominator and paired-statistic acceptance."""
import os
import common as c
def equal(a,b,keys,columns=None):
    import numpy as np
    import pandas as pd
    assert not a.duplicated(keys).any() and not b.duplicated(keys).any()
    a=a.set_index(keys).sort_index();b=b.set_index(keys).sort_index();assert a.index.equals(b.index)
    for key in list(a.columns) if columns is None else columns:
        if pd.api.types.is_numeric_dtype(a[key]) and not pd.api.types.is_bool_dtype(a[key]):
            np.testing.assert_allclose(a[key].to_numpy(float),b[key].to_numpy(float),rtol=0,atol=1e-12,equal_nan=True,err_msg=key)
        else:assert a[key].fillna('<NA>').astype(str).equals(b[key].fillna('<NA>').astype(str)),key

def main():
    c.require_slurm();contract=c.contract();links=c.archive_gate(contract,full=True)
    import numpy as np
    import pandas as pd
    from scipy.stats import wilcoxon
    summary=c.OUT/'summary';sm=c.checked_outputs(summary)
    assert sm['status']=='completed_pending_independent_selection_verification'
    frames=[];anchorframes=[];proofs={};eligibility={}
    for shard in c.js(c.OUT/'shards.json'):
        sample,budget=shard['sample'],shard['budget']
        for space in shard['spaces']:
            c.unit_acceptance(sample,budget,space,links)
            ep=c.OUT/'evaluation'/sample/budget/space;vp=c.OUT/'verification'/sample/budget/space
            em=c.checked_outputs(ep);vm=c.checked_outputs(vp);assert em['n_valid']==39 and vm['status']=='passed' and vm['n_unavailable']==0
            for p in [ep/'manifest.json',vp/'manifest.json']:proofs[str(p)]=c.sha(p)
            frames.append(pd.read_csv(ep/'metrics.csv.gz',dtype={'cutoff':str}));anchorframes.append(pd.read_csv(ep/'anchor_metrics.csv',dtype={'cutoff':str}))
            e=c.js(ep/'eligibility.json')
            if sample in eligibility:assert eligibility[sample]==e
            else:eligibility[sample]=e
    data=pd.concat(frames,ignore_index=True);anchors=pd.concat(anchorframes,ignore_index=True).drop_duplicates()
    assert len(data)==66066 and len(anchors)==726 and len(proofs)==3388
    keys=['sample','budget','space','method','k','stage']
    equal(data,pd.read_csv(summary/'all_candidate_metrics.csv.gz',dtype={'cutoff':str}),keys)
    equal(anchors,pd.read_csv(summary/'anchor_metrics.csv.gz',dtype={'cutoff':str}),keys)
    assert eligibility==c.js(summary/'sample_eligibility.json')
    folds=pd.read_csv(c.CAMP/'protocol/embedding_patient_folds.csv');foldmap=folds.drop_duplicates('patient').set_index('patient').fold.to_dict()
    spec=c.js(c.CAMP/'protocol/embedding.json');protocol=c.js(c.OUT/'protocol.json')
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
        assert len(selected)==n_samples*126 and not selected.duplicated(['sample','budget','space','method','stage']).any()
        selected_frames.append(selected)
        grouping=['cohort','patient','budget','space','method','stage']
        patients=selected.groupby(grouping)[c.MEASURES].mean()
        patients['n_samples_total']=selected.groupby(grouping).size()
        for metric in c.MEASURES:patients[metric+'_n_samples']=selected.groupby(grouping)[metric].count()
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
    equal(pd.DataFrame(choices),pd.read_csv(summary/'patient_fold_K_choices.csv'),['cohort','budget','space','method','fold'])
    equal(pd.concat(selected_frames),pd.read_csv(summary/'selected_sample_metrics.csv.gz',dtype={'cutoff':str}),['cohort','sample','budget','space','method','stage'])
    equal(pd.concat(patient_frames),pd.read_csv(summary/'patient_metrics.csv.gz'),['cohort','patient','budget','space','method','stage'])
    equal(stats,pd.read_csv(summary/'paired_vs_original_anchor.csv'),['cohort','budget','space','method'])
    assert len(choices)==420 and len(stats)==84
    out=c.OUT/'summary_verification';out.mkdir(parents=True,exist_ok=True)
    stats.to_csv(out/'independent_paired_statistics.csv',index=False);pd.DataFrame(choices).to_csv(out/'independent_K_choices.csv',index=False)
    c.write(out/'independent_unit_proofs.json',proofs)
    manifest=dict(status='passed',contract_sha256=c.sha(c.OUT/'contract.json'),endpoint=contract['endpoint'],
        inputs={str(summary/'manifest.json'):c.sha(summary/'manifest.json'),str(c.ARCHIVE/'summary/manifest.json'):c.sha(c.ARCHIVE/'summary/manifest.json')},
        outputs={name:c.sha(out/name) for name in ['independent_paired_statistics.csv','independent_K_choices.csv','independent_unit_proofs.json']},
        n_representations=1694,n_candidate_rows=66066,n_unique_partitions=22022,n_selection_rows=420,n_paired_contrasts=84,
        all_cohorts_and_denominators_verified=True,Lfine_training_patient_selection_verified=True,L1_K_not_reused=True,
        metric_and_per_class_verification_required=True,paired_statistics_and_Holm_verified=True,
        inference='Conditional on selected predictions; retrospective patient-label-held-out selection',
        no_fitting=True,job=os.environ['SLURM_JOB_ID'],completed_at=c.utc())
    c.write(out/'manifest.json',manifest);c.complete(out);print(manifest,flush=True)
if __name__=='__main__':main()
