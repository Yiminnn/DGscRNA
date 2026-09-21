"""Full independent A1 unit acceptance, extending the proven three-pilot audit."""
from pathlib import Path
import contextlib,json,os
import common as c
from unit_extra import manual_classes,partition_scores,close,check_figures,solver_audit,terminal_state

def run(task,contract):
    c.require_slurm()
    import numpy as np
    import pandas as pd
    import pilot_baseline as baseline
    unit=c.unitpath(task);out=c.reportpath(task)
    assert c.ready(task),'Incomplete source unit; no missing condition may be silently dropped'
    if c.valid_report(task,deep=True):return c.js(out/'manifest.json')
    before=c.signature(task);out.mkdir(parents=True,exist_ok=True)
    with (out/'pilot_baseline.log').open('w') as log,contextlib.redirect_stdout(log):
        baseline.run(task['sample'],task['budget'],task['space'])
    cfg=c.js(unit/'config.json');prep=Path(cfg['prep']);rep=c.js(unit/'representation.json')
    assert cfg['dest']==str(unit) and cfg['prep']==str(c.REFERENCE/'GBM'/task['sample']/task['budget'])
    assert cfg['frozen_protocol_sha256']==c.sha(c.CAMP/'protocol/embedding.json') and cfg['reference_labels_used_for_fit'] is False
    pm=c.js(prep/'prepare_manifest.json');im=c.js(c.REFERENCE/'inputs'/task['sample']/'input_manifest.json')
    assert type(im['primary']) is bool
    cells=pd.read_csv(prep/'cells.csv',dtype=str,keep_default_na=False).cell_id.to_numpy()
    truth=pd.read_csv(c.REFERENCE/'evaluation_inputs'/task['sample']/'truth.csv.gz',dtype=str,keep_default_na=False)
    assert np.array_equal(cells,truth.cell_id) and len(set(cells))==len(cells)
    assert c.sha(c.REFERENCE/'evaluation_inputs'/task['sample']/'truth.csv.gz')==im['evaluation_files']['truth.csv.gz']
    gm=c.js(Path(cfg['geometry'])/'manifest.json')
    assert c.sha(Path(cfg['geometry'])/'cells.csv')==gm['cells_sha256']
    assert c.sha(Path(cfg['geometry'])/'features.txt')==gm['features_sha256']
    assert gm['n_cells']==len(cells) and rep['n_cells']==len(cells)
    mapping=pd.read_csv(c.REFERENCE/'markers/panel_L1_mapping.csv',dtype=str,keep_default_na=False)
    lookup=dict(zip(mapping.loc[mapping.library.eq('CM2_glioma_other'),'panel'],mapping.loc[mapping.library.eq('CM2_glioma_other'),'L1']))
    recorded=pd.read_csv(unit/'evaluation/metrics.csv',dtype={'cutoff':str})
    evaluation=c.js(unit/'evaluation/manifest.json')
    assert evaluation['source_sha256']==c.sha(c.CAMP/'source_snapshots/embedding_v8/legacy/evaluate.py')
    perclass=pd.read_csv(unit/'evaluation/per_class.csv.gz',dtype={'cutoff':str})
    confusions=pd.read_csv(unit/'evaluation/confusions.csv.gz',dtype={'cutoff':str})
    clustering=pd.read_csv(unit/'evaluation/clustering.csv')
    assert pd.api.types.is_bool_dtype(recorded.primary.dtype) and recorded.primary.eq(im['primary']).all()
    assert recorded.patient.eq(im['patient']).all() and recorded.status.eq('completed').all()
    assert not recorded.duplicated(['route','stage']).any() and not perclass.duplicated(['route','stage','label']).any()
    assert len(perclass)==13*3*11 and not clustering.route.duplicated().any()
    solver=solver_audit(unit,cfg,rep,contract)
    artifacts={prep/'prepare_manifest.json',prep/'cells.csv',prep/'expression_PCA30.rds',
        Path(cfg['geometry'])/'manifest.json',Path(cfg['geometry'])/'cells.csv',Path(cfg['geometry'])/'features.txt',Path(gm['binary']),
        c.REFERENCE/'inputs'/task['sample']/'input_manifest.json',c.REFERENCE/'evaluation_inputs'/task['sample']/'truth.csv.gz',
        c.REFERENCE/'markers/panel_L1_mapping.csv',c.REFERENCE/'markers/libraries.json',unit/'cells.csv'}
    artifacts.update(unit/name for name in before)
    artifacts.update(unit/'evaluation'/name for name in evaluation['outputs'])
    if task['space']!='noDR':artifacts.add(unit/'embedding.csv')
    if task['space']=='ICA2':
        artifacts.update(Path(a['path']) for a in rep.get('attempts',[]))
        if (unit/'CONVERGENCE_FAILURE.json').is_file():artifacts.add(unit/'CONVERGENCE_FAILURE.json')
    old_fig=check_figures(unit/'figures',cfg['conditions'])
    assert old_fig['source_sha256']==c.sha(c.CAMP/'source_snapshots/embedding_v8/evaluate_and_plot.py')
    assert old_fig['n_cells']==len(cells)
    artifacts.add(Path(old_fig['display']));artifacts.update(unit/'figures'/name for name in old_fig['files'])
    if task['space']=='ICA2':
        adaptive=check_figures(unit/'figures_adaptive_v1',cfg['conditions'])
        for key,path in [('representation_manifest_sha256',unit/'representation.json'),('fit_manifest_sha256',unit/'fit_manifest.json'),
            ('evaluation_manifest_sha256',unit/'evaluation/manifest.json'),('policy_sha256',c.POLICY),
            ('presentation_source_manifest_sha256',c.CAMP/'source_snapshots/embedding_presentation_v1/SOURCE_MANIFEST.json')]:assert adaptive[key]==c.sha(path),key
        assert adaptive['space_display']=='ICA2 adaptive' and adaptive['old_figures_preserved'] is True and adaptive['metrics_recomputed'] is False
        assert adaptive['source_sha256']==c.sha(c.CAMP/'source_snapshots/embedding_presentation_v1/render_adaptive_figures.py') and adaptive['n_cells']==len(cells)
        artifacts.add(Path(adaptive['display']));artifacts.update(unit/'figures_adaptive_v1'/name for name in adaptive['files'])
        for key in ['actual_solver','actual_iteration_cap','canonical_parallel5000_status','fallback_used']:assert adaptive[key]==solver[key]
    y=truth.L1.to_numpy();coarse=lambda a:np.asarray(['Neuron' if v in ['Excitatory neuron','Inhibitory neuron','AMBIGUOUS_NEURON'] else v for v in a])
    yc=coarse(y);coarse_labels=[v for v in c.L1 if v not in ['Excitatory neuron','Inhibitory neuron']]+['Neuron']
    recomputed=[];states=[];reclusters=[]
    for condition in cfg['conditions']:
        route=Path(condition['dest']);name=route.name;td=route/'terminal/L00_mean'
        assert route==unit/name and name in c.ROUTES
        assert condition['route']==task['space']+'_'+name
        assert condition['method']==('HDBSCAN_R' if name=='HDBSCAN_R' else name.split('_K')[0])
        assert condition['k']==(None if name=='HDBSCAN_R' else int(name.split('_K')[1]))
        cl=pd.read_csv(route/'clusters.csv',dtype=str,keep_default_na=False)
        initial=pd.read_csv(route/'initial_calls.csv.gz',dtype=str,keep_default_na=False)
        pred=pd.read_csv(td/'predictions.csv.gz',dtype=str,keep_default_na=False)
        tm=c.js(td/'terminal_manifest.json');sm=c.js(route/'score_manifest.json')
        artifacts.update(route/name for name in ['partition_manifest.json','clusters.csv','cells.csv','score_manifest.json','score_input_fingerprint.json','initial_calls.csv.gz'])
        artifacts.add(Path(sm['DL_binary']));artifacts.add(Path(sm['execution_source']))
        artifacts.update(td/name for name in ['terminal_manifest.json','training_manifest.json','terminal.npz','predictions.csv.gz'])
        artifacts.update(td/name for name in c.js(td/'training_manifest.json')['outputs'])
        assert sm['source_sha256']==c.sha(c.CAMP/'source_snapshots/embedding_v8/score_candidates.R')
        assert set(sm['arms'])=={'L00_mean'} and sm['arms']['L00_mean']['library']=='CM2_glioma_other' and sm['arms']['L00_mean']['cutoff']=='mean'
        part=c.js(route/'partition_manifest.json')
        assert part['status']=='completed' and part['route']==condition['route'] and part['method']==condition['method'] and part['k']==condition['k']
        if condition['method']=='HDBSCAN_R':
            assert part['implementation']=='dbscan::hdbscan' and part['minPts']==50 and part['noise_label']==0 and part['noise_scored_as_cluster'] is True
        else:
            expected=dict(n_clusters=condition['k'],n_init=10,max_iter=300,tol=1e-4,algorithm='lloyd',random_state=42) if condition['method']=='KMeans' else dict(n_components=condition['k'],covariance_type='diag',reg_covar=1e-4,max_iter=1000,tol=1e-3,n_init=1,random_state=42)
            for field,value in expected.items():assert part['params'][field]==value,(name,field)
            if condition['method']=='GMM':assert part['converged_'] is True
        assert tm['DL_features']==pm['features']['DL'] and sm['reference_labels_used_for_fit'] is False
        state=terminal_state(td,initial.L00_mean.to_numpy(),pred,tm,cells,contract)
        states.append(dict(route=name,**state))
        known=initial.L00_mean.to_numpy()!='Undecided'
        for stage,column in c.STAGES.items():
            row=recorded[recorded.route.eq(name)&recorded.stage.eq(stage)].iloc[0]
            native=pred[column].to_numpy();mapped=np.asarray([lookup.get(v,'Unknown' if v in ['Unknown','Undecided','Noise',''] else 'UNMAPPABLE') for v in native])
            classes=manual_classes(y,mapped,c.L1);f1=np.asarray([r['F1'] for r in classes]);support=np.asarray([r['support'] for r in classes])
            coarse_classes=manual_classes(yc,coarse(mapped),coarse_labels);cf=np.asarray([r['F1'] for r in coarse_classes]);cs=np.asarray([r['support'] for r in coarse_classes])
            abstain=np.isin(native,['Unknown','Undecided','Noise','']);supported=np.isin(mapped,c.L1)
            values=dict(macroF1_present=f1[support>0].mean(),macroF1_fixed11=f1.mean(),weightedF1=np.average(f1,weights=support),
                accuracy=(mapped==y).mean(),coverage=(~abstain).mean(),unknown_rate=abstain.mean(),mapped_coverage=supported.mean(),
                off_vocabulary_rate=(~supported&~abstain).mean(),coarse10_macroF1_present=cf[cs>0].mean(),coarse10_macroF1_fixed10=cf.mean(),
                n_known=int(known.sum()),n_pool=int((~known).sum()),n_training_classes=tm['n_training_classes'],
                n_new_correct=int(((~known)&(mapped==y)).sum()),n_new_incorrect=int(((~known)&(~abstain)&(mapped!=y)).sum()),
                n_initially_wrong_retained=int((known&(mapped!=y)).sum()))
            for field,value in values.items():close(value,row[field],(name,stage,field))
            context=dict(sample=task['sample'],patient=im['patient'],primary=im['primary'],budget=task['budget'],route=name,
                library='CM2_glioma_other',cutoff='mean',arm_id='L00_mean',stage=stage,family='native_R_budget',seed=42,n_cells=len(cells),DL_features=pm['features']['DL'],
                status='completed',dl_status=tm['dl_status'],training_executed=tm['training_executed'],reference_overlap=False)
            for field,value in context.items():assert row[field]==value,(name,stage,field,row[field],value)
            recomputed.append(dict(**context,**values))
            sub=perclass[perclass.route.eq(name)&perclass.stage.eq(stage)].set_index('label')
            assert set(sub.index)==set(c.L1)
            for item in classes:
                for field in ['precision','recall','F1','support']:close(item[field],sub.loc[item['label'],field],(name,stage,item['label'],field))
            expected=CounterPairs(y,mapped)
            sub=confusions[confusions.route.eq(name)&confusions.stage.eq(stage)]
            actual={(r.truth,r.prediction):int(r.n) for r in sub.itertuples()}
            assert len(actual)==len(sub) and actual==expected and sum(actual.values())==len(cells)
        scores=partition_scores(y,cl.cluster);fine=partition_scores(truth.lfine_original,cl.cluster)
        malignant=y=='Malignant';ms=partition_scores(truth.loc[malignant,'MalState'],cl.loc[malignant,'cluster'])['ARI'] if malignant.sum()>1 else None
        values=dict(n_cells=len(cells),n_clusters=cl.loc[cl.cluster.ne('0') if condition['method']=='HDBSCAN_R' else cl.cluster.notna(),'cluster'].nunique(),
            n_partition_ids=cl.cluster.nunique(),noise_rate=cl.cluster.eq('0').mean() if condition['method']=='HDBSCAN_R' else 0.,
            L1_ARI=scores['ARI'],L1_NMI=scores['NMI'],L1_FMI=scores['FMI'],fine_ARI=fine['ARI'],malignant_state_ARI=ms)
        row=clustering[clustering.route.eq(name)].iloc[0]
        for field,value in values.items():close(value,row[field],(name,field))
        reclusters.append(dict(sample=task['sample'],patient=im['patient'],primary=im['primary'],budget=task['budget'],route=name,evaluation_type='clustering',**values))
    assert c.signature(task)==before,'Source manifests changed during validation'
    # Preserve independent recomputation and the parity-verified recorded floats.
    # The latter avoids CSV-rounding artifacts in exact-tie K reconstruction.
    pd.DataFrame(recomputed).to_csv(out/'independent_metrics.csv.gz',index=False,compression='gzip')
    recorded.to_csv(out/'validated_metrics.csv.gz',index=False,compression='gzip')
    pd.DataFrame(reclusters).to_csv(out/'independent_clustering.csv.gz',index=False,compression='gzip')
    clustering.to_csv(out/'validated_clustering.csv.gz',index=False,compression='gzip')
    c.write(out/'solver.json',solver)
    c.write(out/'source_artifacts.json',{str(path):baseline.sha(path) for path in sorted(artifacts)})
    baseline_path=c.OUT/'pilot_baseline'/task['sample']/task['budget']/task['space']/'validation.json'
    report=dict(status='passed',**task,n_conditions=13,n_metric_rows=39,n_cells_each=len(cells),patient=im['patient'],primary=im['primary'],
        states=states,input_manifests=before,validation_contract_sha256=c.sha(c.OUT/'contract.json'),pilot_baseline_proof_sha256=c.sha(baseline_path),
        all_initial_terminal090_terminal070_metrics=True,all_per_class_confusions_coarse_and_partition_metrics=True,
        training_splits_and_thresholds_reconstructed=True,invalid_terminal_states_accepted=0,all_cells_in_denominator=True,
        old_figure_files_verified=26,adaptive_figure_files_verified=26 if task['space']=='ICA2' else 0,
        plot_artifact_checks='All hashes and PNG/PDF format checks; does not claim visual review of every figure',
        no_model_fits=True,job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),completed_at=c.utc(),
        outputs={name:c.sha(out/name) for name in ['independent_metrics.csv.gz','validated_metrics.csv.gz','independent_clustering.csv.gz','validated_clustering.csv.gz','solver.json','source_artifacts.json']})
    c.write(out/'manifest.json',report);c.complete(out);return report

def CounterPairs(a,b):
    from collections import Counter
    return dict(Counter(zip(a,b)))
