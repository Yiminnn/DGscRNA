"""Independent patient selection and paired inference; reads only validated fits."""
from pathlib import Path
import os,json
import common as c
from unit_extra import manual_classes,close

def equal_table(expected,actual,keys,columns=None,tol=1e-12):
    import numpy as np
    import pandas as pd
    assert not expected.duplicated(keys).any() and not actual.duplicated(keys).any(),keys
    a=expected.set_index(keys).sort_index();b=actual.set_index(keys).sort_index()
    assert a.index.equals(b.index),(keys,len(a),len(b))
    cols=list(a.columns) if columns is None else columns
    assert set(cols)<=set(b.columns)
    for name in cols:
        if pd.api.types.is_numeric_dtype(a[name]) and not pd.api.types.is_bool_dtype(a[name]):
            np.testing.assert_allclose(a[name].to_numpy(float),b[name].to_numpy(float),rtol=0,atol=tol,equal_nan=True,err_msg=name)
        else:assert a[name].fillna('<NA>').astype(str).equals(b[name].fillna('<NA>').astype(str)),name

def main():
    c.require_slurm();contract=c.contract()
    import numpy as np
    import pandas as pd
    from scipy.stats import wilcoxon
    tasks=c.js(c.OUT/'tasks.json');assert len(tasks)==1694
    summary=c.EMBED/'summary_adaptive_v1';assert c.checked(summary),'No partial summary is accepted'
    sm=c.js(summary/'manifest.json')
    assert sm['source_sha256']==c.sha(c.CAMP/'source_snapshots/embedding_presentation_v1/select_and_summarize.py')
    assert sm['adaptive_policy_sha256']==c.sha(c.POLICY) and sm['original_summary_preserved'] is True
    assert sm['status']=='completed' and sm['all_ICA_rows_use_adaptive_policy_label'] is True
    assert sm['protocol_sha256']==c.sha(c.CAMP/'protocol/embedding_selection.json')
    assert sm['presentation_source_manifest_sha256']==c.sha(c.CAMP/'source_snapshots/embedding_presentation_v1/SOURCE_MANIFEST.json')
    assert sm['adaptive_policy_name']==c.js(c.POLICY)['comparator_name']
    assert sm['expected_candidate_partitions']==22022 and sm['selected_all121_conditions']==5082
    expected_files={'all_candidate_metrics.csv.gz','all_candidate_clustering.csv.gz','selected_sample_metrics.csv.gz',
        'patient_fold_K_choices.csv','patient_metrics.csv.gz','paired_vs_original_anchor.csv','ICA_solver_execution.csv',
        'paired_embedding_comparison.png','paired_embedding_comparison.pdf'}
    assert set(sm['files'])==expected_files
    for name,digest in sm['files'].items():assert c.sha(summary/name)==digest,name
    for path,digest in sm['input_manifests'].items():assert c.sha(path)==digest,path
    spec=c.js(c.CAMP/'protocol/embedding.json');selection=c.js(c.CAMP/'protocol/embedding_selection.json')
    foldfile=c.CAMP/'protocol/embedding_patient_folds.csv';assert c.sha(foldfile)==selection['patient_folds_sha256']
    foldrows=pd.read_csv(foldfile);assert len(foldrows)==121 and foldrows['sample'].nunique()==121
    assert foldrows.primary.dtype==bool and foldrows.primary.sum()==97
    assert foldrows.patient.nunique()==59 and foldrows.loc[foldrows.primary,'patient'].nunique()==55
    assert foldrows.groupby('patient').fold.nunique().eq(1).all()
    foldmap=foldrows.drop_duplicates('patient').set_index('patient').fold.to_dict()
    metadata=foldrows.set_index('sample');frames=[];clusterframes=[];solverrows=[];proofs={};state_counts={}
    for task in tasks:
        assert c.valid_report(task,deep=True),'Every one of 1694 unit proofs and raw artifact hashes is required'
        out=c.reportpath(task);proof=c.js(out/'manifest.json');proofs[str(out/'manifest.json')]=c.sha(out/'manifest.json')
        assert proof['n_conditions']==13 and proof['n_metric_rows']==39
        assert proof['patient']==metadata.loc[task['sample'],'patient'] and proof['primary']==bool(metadata.loc[task['sample'],'primary'])
        assert proof['n_cells_each']==int(metadata.loc[task['sample'],'n_cells'])
        context=[]
        for route in c.ROUTES:
            method='HDBSCAN_R' if route=='HDBSCAN_R' else route.split('_K')[0]
            context.append(dict(route=route,method=method,k=0 if method=='HDBSCAN_R' else int(route.split('_K')[1])))
        context=pd.DataFrame(context)
        frame=pd.read_csv(out/'validated_metrics.csv.gz',dtype={'cutoff':str}).merge(context,on='route',validate='many_to_one')
        frame['space']=task['space'];frame['space_display']='ICA2 adaptive' if task['space']=='ICA2' else task['space'];frames.append(frame)
        cp=pd.read_csv(out/'validated_clustering.csv.gz').merge(context,on='route',validate='one_to_one');cp['space']=task['space'];cp['space_display']=frame.space_display.iloc[0];clusterframes.append(cp)
        if task['space']=='ICA2':solverrows.append(c.js(out/'solver.json'))
        for row in proof['states']:state_counts[row['dl_status']]=state_counts.get(row['dl_status'],0)+1
    allmetrics=pd.concat(frames,ignore_index=True);allclusters=pd.concat(clusterframes,ignore_index=True)
    assert len(allmetrics)==66066 and len(allclusters)==22022 and len(solverrows)==242
    actual=pd.read_csv(summary/'all_candidate_metrics.csv.gz',dtype={'cutoff':str})
    equal_table(allmetrics,actual,['sample','budget','space','method','k','stage'])
    equal_table(allclusters,pd.read_csv(summary/'all_candidate_clustering.csv.gz'),['sample','budget','space','method','k'])
    sr=pd.DataFrame(solverrows)
    equal_table(sr,pd.read_csv(summary/'ICA_solver_execution.csv'),['sample','budget','space'],columns=[x for x in sr.columns if x not in ['sample','budget','space','old_failure_path','old_failure_sha256']])
    # Independently reconstruct the original R anchor's two reported measures
    # from its frozen native final090 labels, before any patient aggregation.
    mapping=pd.read_csv(c.REFERENCE/'markers/panel_L1_mapping.csv',dtype=str,keep_default_na=False)
    m=mapping[mapping.library.eq('CM2_glioma_other')];lookup=dict(zip(m.panel,m.L1));anchorrows=[]
    for sample in spec['samples']:
        truth=pd.read_csv(c.REFERENCE/'evaluation_inputs'/sample/'truth.csv.gz',dtype=str,keep_default_na=False);y=truth.L1.to_numpy()
        for budget in spec['budgets']:
            prep=c.REFERENCE/'GBM'/sample/budget;ev=c.js(prep/'evaluation/manifest.json');assert c.checked(prep/'evaluation')
            assert c.sha(prep/'evaluation/metrics.csv')==ev['outputs']['metrics.csv']
            rows=pd.read_csv(prep/'evaluation/metrics.csv',dtype={'cutoff':str})
            chosen=rows[rows.route.eq('UMAP2_HDBSCAN_R')&rows.library.eq('CM2_glioma_other')&rows.cutoff.eq('mean')&rows.stage.eq('terminal090')&rows.family.eq('native_R_budget')]
            assert len(chosen)==1;row=chosen.iloc[0]
            route=prep/'UMAP2_HDBSCAN_R';score=c.js(route/'score_manifest.json');assert c.checked(route,'score_manifest.json','SCORE_COMPLETE')
            arms=[a for a,v in score['arms'].items() if v['library']=='CM2_glioma_other' and v['cutoff']=='mean'];assert len(arms)==1
            td=route/'terminal'/arms[0];tm=c.js(td/'terminal_manifest.json');assert c.checked(td,'terminal_manifest.json','TERMINAL_COMPLETE')
            assert tm['score_manifest_sha256']==c.sha(route/'score_manifest.json') and tm['terminal_sha256']==c.sha(td/'terminal.npz')
            cells=pd.read_csv(route/'cells.csv',dtype=str,keep_default_na=False).cell_id.to_numpy();assert np.array_equal(cells,truth.cell_id)
            with np.load(td/'terminal.npz',allow_pickle=False) as z:native=z['final090'].copy()
            pred=np.asarray([lookup.get(v,'Unknown' if v in ['Unknown','Undecided','Noise',''] else 'UNMAPPABLE') for v in native])
            classes=manual_classes(y,pred,c.L1);f1=np.asarray([v['F1'] for v in classes]);support=np.asarray([v['support'] for v in classes])
            ff=f1[support>0].mean();coverage=(~np.isin(native,['Unknown','Undecided','Noise',''])).mean()
            close(ff,row.macroF1_present,(sample,budget,'anchorF1'));close(coverage,row.coverage,(sample,budget,'anchorCoverage'))
            assert row.patient==metadata.loc[sample,'patient'] and row.primary==metadata.loc[sample,'primary'] and row.n_cells==len(cells)
            anchorrows.append(dict(sample=sample,budget=budget,patient=row.patient,primary=bool(row.primary),macroF1_present=row.macroF1_present,coverage=row.coverage,
                terminal_sha256=tm['terminal_sha256'],evaluation_manifest_sha256=c.sha(prep/'evaluation/manifest.json')))
    anchors=pd.DataFrame(anchorrows)
    choices=[];selected_frames=[];patient_frames=[];paired=[];rng=np.random.default_rng(20260920)
    for cohort in ['primary97','all121']:
        data=allmetrics[allmetrics.primary] if cohort=='primary97' else allmetrics
        expected_samples=97 if cohort=='primary97' else 121;expected_patients=55 if cohort=='primary97' else 59
        assert data['sample'].nunique()==expected_samples and data.patient.nunique()==expected_patients
        terminal=data[data.stage.eq('terminal090')]
        patient_candidates=terminal.groupby(['patient','budget','space','method','k'])[c.MEASURES].mean().reset_index()
        patient_candidates['fold']=patient_candidates.patient.map(foldmap)
        rows=[]
        for budget,space,method in sorted(set(zip(patient_candidates.budget,patient_candidates.space,patient_candidates.method))):
            group=patient_candidates[(patient_candidates.budget==budget)&(patient_candidates.space==space)&(patient_candidates.method==method)]
            for fold in sorted(group.fold.unique()):
                train=group[group.fold.ne(fold)];held=group[group.fold.eq(fold)]
                scores={int(k):float(v) for k,v in train.groupby('k').macroF1_present.mean().items()}
                assert set(scores)==({0} if method=='HDBSCAN_R' else set(c.KS))
                selected_k=min(scores,key=lambda k:(-scores[k],k))
                train_patients=set(train.patient);test_patients=set(held.patient)
                assert train_patients.isdisjoint(test_patients) and len(train_patients|test_patients)==expected_patients
                choices.append(dict(cohort=cohort,budget=budget,space=space,space_display='ICA2 adaptive' if space=='ICA2' else space,method=method,
                    fold=int(fold),selected_k=selected_k,training_patient_mean_macroF1=scores[selected_k],n_training_patients=len(train_patients),n_test_patients=len(test_patients)))
                take=data[(data.budget==budget)&(data.space==space)&(data.method==method)&(data.k==selected_k)&data.patient.isin(test_patients)].copy()
                take['fold']=int(fold);take['cohort']=cohort;rows.append(take)
        selected=pd.concat(rows,ignore_index=True)
        assert len(selected)==expected_samples*2*7*3*3 and not selected.duplicated(['sample','budget','space','method','stage']).any()
        selected_frames.append(selected)
        patients=selected.groupby(['cohort','patient','budget','space','method','stage'])[c.MEASURES].mean().reset_index()
        patients['space_display']=patients.space.map(lambda x:'ICA2 adaptive' if x=='ICA2' else x);patient_frames.append(patients)
        for budget in spec['budgets']:
            ref=anchors[anchors.budget.eq(budget)&(anchors.primary if cohort=='primary97' else True)].groupby('patient')[['macroF1_present','coverage']].mean()
            current=patients[patients.budget.eq(budget)&patients.stage.eq('terminal090')]
            for space,method in sorted(set(zip(current.space,current.method))):
                group=current[current.space.eq(space)&current.method.eq(method)].set_index('patient').sort_index()
                assert set(group.index)==set(ref.index) and len(group)==expected_patients
                ref_here=ref.loc[group.index];delta=(group.macroF1_present-ref_here.macroF1_present).to_numpy()
                boot=np.mean(delta[rng.integers(len(delta),size=(10000,len(delta)))],axis=1)
                paired.append(dict(cohort=cohort,budget=budget,space=space,space_display='ICA2 adaptive' if space=='ICA2' else space,method=method,n_patients=len(group),
                    candidate_mean_F1=group.macroF1_present.mean(),anchor_mean_F1=ref_here.macroF1_present.mean(),mean_delta=delta.mean(),
                    CI95_low=np.quantile(boot,.025),CI95_high=np.quantile(boot,.975),p_wilcoxon=wilcoxon(delta).pvalue if (np.abs(delta)>1e-14).any() else 1.,
                    candidate_coverage=group.coverage.mean(),anchor_coverage=ref_here.coverage.mean()))
    stats=pd.DataFrame(paired);stats['p_Holm']=1.
    for _,group in stats.groupby(['cohort','budget']):
        assert len(group)==21
        running=0.
        for rank,index in enumerate(sorted(group.index,key=lambda i:stats.at[i,'p_wilcoxon'])):
            running=max(running,(21-rank)*stats.at[index,'p_wilcoxon']);stats.at[index,'p_Holm']=min(1.,running)
    equal_table(pd.DataFrame(choices),pd.read_csv(summary/'patient_fold_K_choices.csv'),['cohort','budget','space','method','fold'])
    equal_table(pd.concat(selected_frames),pd.read_csv(summary/'selected_sample_metrics.csv.gz',dtype={'cutoff':str}),['cohort','sample','budget','space','method','stage'])
    equal_table(pd.concat(patient_frames),pd.read_csv(summary/'patient_metrics.csv.gz'),['cohort','patient','budget','space','method','stage'])
    equal_table(stats,pd.read_csv(summary/'paired_vs_original_anchor.csv'),['cohort','budget','space','method'])
    assert len(stats)==84 and len(choices)==420 and sm['n_patient_contrasts']==84 and sm['n_selection_rows']==420
    out=c.OUT/'summary';out.mkdir(exist_ok=True)
    stats.to_csv(out/'independent_paired_statistics.csv',index=False);pd.DataFrame(choices).to_csv(out/'independent_K_choices.csv',index=False)
    anchors.to_csv(out/'independent_anchor_checks.csv',index=False)
    c.write(out/'unit_proof_hashes.json',proofs)
    report=dict(status='passed',scope='Archived original L1 scientific protocol; independent full A1 acceptance',n_representations=1694,n_candidates=22022,n_metric_rows=66066,
        cohorts=dict(primary97=dict(samples=97,patients=55),all121=dict(samples=121,patients=59)),
        n_patient_fold_choices=420,n_paired_contrasts=84,Holm_family_size=21,bootstrap_replicates=10000,bootstrap_seed=20260920,
        old_figure_files_verified=44044,adaptive_figure_files_verified=6292,n_ICA_solver_rows=242,DL_state_counts=state_counts,
        independent_all_metrics_verified=True,independent_patient_K_selection=True,independent_paired_statistics_and_Holm=True,
        all_supported_unknown_and_off_vocabulary_cells_retained=True,conditional_on_selected_predictions_inference=True,
        source_summary_manifest_sha256=c.sha(summary/'manifest.json'),validation_contract_sha256=c.sha(c.OUT/'contract.json'),
        outputs={name:c.sha(out/name) for name in ['independent_paired_statistics.csv','independent_K_choices.csv','independent_anchor_checks.csv','unit_proof_hashes.json']},
        visual_scope='Artifact hashes and image formats verified globally; representative visual review remains separately documented',
        no_model_fits=True,job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),completed_at=c.utc())
    c.write(out/'manifest.json',report);c.complete(out);print(json.dumps(report,indent=2))

if __name__=='__main__':main()
