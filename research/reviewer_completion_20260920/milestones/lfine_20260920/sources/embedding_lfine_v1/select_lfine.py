"""Select K using Lfine training-patient scores; require the complete frozen scope."""
import os
import common as c

def main():
    c.require_slurm();contract=c.contract();links=c.archive_gate(contract,full=True)
    import numpy as np
    import pandas as pd
    from scipy.stats import wilcoxon
    out=c.OUT/'summary';assert not c.checked(out),'Preserve the completed summary'
    spec=c.js(c.CAMP/'protocol/embedding.json');protocol=c.js(c.OUT/'protocol.json')
    folds=pd.read_csv(c.CAMP/'protocol/embedding_patient_folds.csv')
    assert len(folds)==121 and folds['sample'].nunique()==121 and folds.primary.dtype==bool and folds.primary.sum()==97
    assert folds.groupby('patient').fold.nunique().eq(1).all()
    foldmap=folds.drop_duplicates('patient').set_index('patient').fold.to_dict()
    frames=[];anchorframes=[];inputs={str(c.ARCHIVE/'summary/manifest.json'):c.sha(c.ARCHIVE/'summary/manifest.json')};eligibility={}
    for shard in c.js(c.OUT/'shards.json'):
        sample,budget=shard['sample'],shard['budget']
        for space in shard['spaces']:
            c.unit_acceptance(sample,budget,space,links)
            directory=c.OUT/'evaluation'/sample/budget/space;verification=c.OUT/'verification'/sample/budget/space
            em=c.checked_outputs(directory);vm=c.checked_outputs(verification)
            assert em['status']=='completed' and vm['status']=='passed' and em['n_valid']==39 and vm['n_unavailable']==0
            for p in [directory/'manifest.json',verification/'manifest.json']:inputs[str(p)]=c.sha(p)
            data=pd.read_csv(directory/'metrics.csv.gz',dtype={'cutoff':str});anchor=pd.read_csv(directory/'anchor_metrics.csv',dtype={'cutoff':str})
            assert len(data)==39 and len(anchor)==3 and data.terminal_valid.all() and anchor.terminal_valid.all()
            assert data.primary.dtype==bool and anchor.primary.dtype==bool
            e=c.js(directory/'eligibility.json')
            if sample in eligibility:assert eligibility[sample]==e
            else:eligibility[sample]=e
            frames.append(data);anchorframes.append(anchor)
    allmetrics=pd.concat(frames,ignore_index=True);anchors=pd.concat(anchorframes,ignore_index=True)
    assert len(allmetrics)==66066 and len(eligibility)==121
    assert not allmetrics.duplicated(['sample','budget','space','method','k','stage']).any()
    assert allmetrics.primary.sum()==97*2*7*13*3 and allmetrics.patient.nunique()==59
    # Anchors are deliberately checked alongside every unit. All seven copies
    # must be identical before reducing to the fixed 242 anchors × three stages.
    for _,g in anchors.groupby(['sample','budget','stage']):
        assert len(g)==7 and len(g.drop_duplicates())==1
    anchors=anchors.drop_duplicates().reset_index(drop=True);assert len(anchors)==726
    assert allmetrics.lfine_macroF1.notna().equals(allmetrics.lfine_metric_eligible)
    assert anchors.lfine_macroF1.notna().equals(anchors.lfine_metric_eligible)
    for s,e in eligibility.items():
        assert e['patient']==folds.set_index('sample').loc[s,'patient'] and e['primary']==bool(folds.set_index('sample').loc[s,'primary'])
    selected_frames=[];choices=[];patients_all=[];paired=[];rng=np.random.default_rng(protocol['bootstrap_seed'])
    for cohort in ['primary97','all121']:
        data=allmetrics[allmetrics.primary] if cohort=='primary97' else allmetrics
        n_samples=97 if cohort=='primary97' else 121;n_patients=55 if cohort=='primary97' else 59
        assert data['sample'].nunique()==n_samples and data.patient.nunique()==n_patients
        terminal=data[data.stage.eq('terminal090')]
        candidates=terminal.groupby(['patient','budget','space','method','k'])[c.MEASURES].mean().reset_index()
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
        assert len(chosen)==n_samples*2*7*3*3 and not chosen.duplicated(['sample','budget','space','method','stage']).any()
        selected_frames.append(chosen)
        patientrows=[]
        for keys,g in chosen.groupby(['cohort','patient','budget','space','method','stage']):
            row=dict(zip(['cohort','patient','budget','space','method','stage'],keys));row['n_samples_total']=len(g)
            row['space_display']='ICA2 adaptive' if row['space']=='ICA2' else row['space']
            for metric in c.MEASURES:
                row[metric]=g[metric].mean();row[metric+'_n_samples']=int(g[metric].notna().sum())
            patientrows.append(row)
        patients=pd.DataFrame(patientrows);patients_all.append(patients)
        for budget in spec['budgets']:
            ref=anchors[anchors.budget.eq(budget)&anchors.stage.eq('terminal090')&(anchors.primary if cohort=='primary97' else True)]
            ref=ref.groupby('patient')[c.MEASURES].mean().reset_index()
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
    out.mkdir(parents=True,exist_ok=True)
    allmetrics.to_csv(out/'all_candidate_metrics.csv.gz',index=False);anchors.to_csv(out/'anchor_metrics.csv.gz',index=False)
    pd.concat(selected_frames).to_csv(out/'selected_sample_metrics.csv.gz',index=False)
    pd.DataFrame(choices).to_csv(out/'patient_fold_K_choices.csv',index=False)
    pd.concat(patients_all).to_csv(out/'patient_metrics.csv.gz',index=False);stats.to_csv(out/'paired_vs_original_anchor.csv',index=False)
    c.write(out/'sample_eligibility.json',eligibility)
    for path,digest in inputs.items():assert c.sha(path)==digest,path
    manifest=dict(status='completed_pending_independent_selection_verification',endpoint=contract['endpoint'],
        n_representations=1694,n_candidate_rows=66066,n_unique_partitions=22022,n_anchor_rows=726,n_selection_rows=420,n_paired_contrasts=84,
        selection_recomputed_with_Lfine_training_patients=True,reused_L1_selected_K=False,
        all_cells_retained=True,eligibility='Truth-derived support20 class eligibility only; candidate-dependent missingness rejected',
        scope='Two HVG budgets, fixed marker context; not an all-gene fit or strict fine-label annotation',
        inference='Retrospective held-out patient K selection, conditional-on-selected-predictions paired inference',
        contract_sha256=c.sha(c.OUT/'contract.json'),inputs=inputs,
        outputs={p.name:c.sha(p) for p in out.iterdir() if p.suffix in ['.csv','.gz','.json'] and p.name!='manifest.json'},
        job=os.environ['SLURM_JOB_ID'],completed_at=c.utc(),no_fitting=True)
    c.write(out/'manifest.json',manifest);c.complete(out);print(manifest,flush=True)
if __name__=='__main__':main()
