"""PTC paired-patient evidence, preserving four patients and all counterexamples."""
import itertools
import json
import os
from pathlib import Path
from common import OUT,OLD,sha,checked,write_json,complete,utc
from ptc_followup_common import PTC,ANCHORS,CONTROL_ROUTES,require_ptc,original_arm

def terminal_execution_tables(mlp_tasks,preparations):
    """One row per requested follow-up condition, without patient/stage duplication."""
    require_ptc()
    import pandas as pd
    from ptc_control_job import context_arms
    conditions=[]
    for cfg in mlp_tasks:
        assert checked(Path(cfg['dest']))
        conditions.extend((cfg,Path(cfg['source']),Path(cfg['dest'])/'terminal'/aid,aid,cfg['seed'])
                          for aid in cfg['arm_ids'])
    for cfg in preparations:
        assert checked(Path(cfg['dest'])/'evaluation')
        for route in CONTROL_ROUTES:
            source=Path(cfg['dest'])/route
            conditions.extend((cfg,source,source/'terminal'/aid,aid,42) for aid in context_arms(cfg,route))
    nontraining={'no_op_all_initially_known':'no_op_all_initially_known',
        'no_known_labels_archived_Undecided_terminal':'untrainable_no_known_labels',
        'structural_insufficient_known_split':'untrainable_insufficient_known_split'}
    rows=[]
    for cfg,source,d,aid,seed in conditions:
        assert checked(d,'terminal_manifest.json','TERMINAL_COMPLETE')
        tm=json.loads((d/'terminal_manifest.json').read_text())
        training=json.loads((d/'training_manifest.json').read_text())
        assert tm['source']==str(source) and tm['model_seed']==seed and tm['split_seed']==42
        assert tm['score_manifest_sha256']==sha(source/'score_manifest.json')
        assert tm['training_manifest_sha256']==sha(d/'training_manifest.json')
        for key in ['dl_status','training_executed','n_known','n_pool','n_training_classes']:
            assert tm[key]==training[key],(str(d),key)
        status=tm['dl_status'];trained=tm['training_executed'];reused=tm['identical_result_reused']
        assert isinstance(trained,bool) and isinstance(reused,bool)
        assert trained==(status in ['trained','trained_single_known_class'])
        assert trained or status in nontraining,status
        action='cached_terminal_reuse' if reused else ('fresh_training' if trained else nontraining[status])
        if not training['terminal_valid']:action='invalid_terminal'
        rows.append(dict(family=cfg['kind'],group=cfg['group'],control_name=cfg['name'],
            source_unit=source.parent.name,route=source.name,arm_id=aid,
            library=tm['arm']['library'],cutoff=tm['arm']['cutoff'],model_seed=seed,split_seed=42,
            dl_status=status,condition_action=action,cached_training_executed=trained,
            identical_result_reused=reused,training_executed_in_this_condition=trained and not reused,
            terminal_valid=training['terminal_valid'],n_cells=training['n_cells'],n_known=tm['n_known'],
            n_pool=tm['n_pool'],n_training_classes=tm['n_training_classes'],
            terminal_directory=str(d),terminal_manifest_sha256=sha(d/'terminal_manifest.json'),
            training_manifest_sha256=tm['training_manifest_sha256'],cache_key=tm['cache_key'],
            condition_job=tm['job'],original_training_or_noop_job=training['job']))
    ledger=pd.DataFrame(rows)
    assert len(ledger)==len(conditions) and ledger.terminal_directory.is_unique
    summary=ledger.groupby(['family','group','dl_status','condition_action','terminal_valid'],dropna=False).size().rename('n_conditions').reset_index()
    assert int(summary.n_conditions.sum())==len(ledger)
    return ledger,summary

def paired(delta):
    import numpy as np
    v=np.asarray(delta,dtype=float);assert v.shape==(4,) and np.isfinite(v).all()
    indices=np.asarray(list(itertools.product(range(4),repeat=4)))
    draws=v[indices].mean(1)
    signs=np.asarray(list(itertools.product([-1,1],repeat=4)))
    observed=abs(float(v.mean()));permuted=np.abs((signs*v).mean(1))
    return dict(n_patients=4,mean_delta=float(v.mean()),CI95_low=float(np.quantile(draws,.025)),
        CI95_high=float(np.quantile(draws,.975)),exact_sign_flip_p=float((permuted>=observed-1e-14).mean()),
        bootstrap_resamples=256,sign_patterns=16)

def run():
    require_ptc()
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    dest=OUT/'PTC_summary';dest.mkdir(exist_ok=True)
    assert checked(PTC/'selection')
    existing=pd.read_csv(PTC/'selection/all_existing_patient_metrics.csv.gz',keep_default_na=False)
    choices=pd.read_csv(PTC/'selection/marker_choices_by_patient.csv',keep_default_na=False)
    workflow=pd.read_csv(PTC/'selection/workflow_choices_by_patient.csv',keep_default_na=False)
    primary=lambda f:f[(f.truth_definition=='productive_TCR')&(f.endpoint=='strict_T_name_rule')&(f.patient!='ALL')]
    metric_columns=['F1_T','F1_nonT','F1_nonT_unknown_as_error','macro_F1_unknown_as_error',
                    'accuracy_unknown_as_error','TCR_positive_recall','TCR_detection_yield_in_predicted_T','coverage']
    finals=primary(existing);finals=finals[finals.stage=='final090']
    chosen_rows=[];baseline_rows=[];core_rows=[];contrasts=[];contrast_differences=[]
    for group,a in ANCHORS.items():
        original=finals[(finals.group==group)&(finals.unit=='PTC_archived_CCA2000')&
            (finals.route==a['route'])&(finals.library==a['library'])&(finals.cutoff==a['cutoff'])]
        assert len(original)==4
        baseline_rows.append(original)
        for row in workflow[workflow.group==group].itertuples():
            q=finals[(finals.group==group)&(finals.patient==row.held_patient)&(finals.unit==row.unit)&
                (finals.route==row.route)&(finals.arm_id==row.arm_id)];assert len(q)==1;chosen_rows.append(q)
        for budget in ['500','1000','2000','3000','5000','all']:
            unit=f'PTC_{group}_CCA{budget}'
            for route in ['PCA30_SNN','PCA30_HDBSCAN_R','UMAP2_SNN','UMAP2_HDBSCAN_R']:
                d=finals[(finals.group==group)&(finals.unit==unit)&(finals.route==route)]
                for mode in ['fixed_historical_context','training_patient_selected_context']:
                    if mode=='fixed_historical_context':selected=d[(d.library==a['library'])&(d.cutoff==a['cutoff'])].copy()
                    else:
                        c=choices[(choices.group==group)&(choices.unit==unit)&(choices.route==route)]
                        selected=d.merge(c[['held_patient','arm_id']],left_on=['patient','arm_id'],right_on=['held_patient','arm_id'],validate='one_to_one')
                    assert len(selected)==4
                    selected=selected.copy();selected['marker_choice']=mode;selected['budget']=budget
                    core_rows.append(selected)
    baseline=pd.concat(baseline_rows,ignore_index=True);tuned=pd.concat(chosen_rows,ignore_index=True);core=pd.concat(core_rows,ignore_index=True)
    baseline.to_csv(dest/'original_anchor_patient_metrics.csv',index=False)
    tuned.to_csv(dest/'heldout_workflow_patient_metrics.csv',index=False)
    core.to_csv(dest/'fresh24_patient_metrics.csv',index=False)
    core_summary=core.groupby(['group','marker_choice','budget','route'])[metric_columns].mean().reset_index()
    core_summary.to_csv(dest/'fresh24_patient_summary.csv',index=False)
    for group in ANCHORS:
        left=tuned[tuned.group==group];right=baseline[baseline.group==group]
        p=left.merge(right,on='patient',suffixes=('_candidate','_original'),validate='one_to_one').sort_values('patient')
        contrasts.append(dict(group=group,contrast='heldout_fresh24_workflow_minus_fixed_archived_anchor',metric='F1_T',
            **paired(p.F1_T_candidate-p.F1_T_original)))
        for row in p.itertuples():
            contrast_differences.append(dict(group=group,patient=row.patient,
                contrast='heldout_fresh24_workflow_minus_fixed_archived_anchor',metric='F1_T',
                delta=row.F1_T_candidate-row.F1_T_original))
    # Matched marker x DL factorial: both stages use the terminal-selected marker.
    factorial=pd.read_csv(PTC/'selection/marker_DL_factorial.csv.gz',keep_default_na=False)
    interactions=[];patient_effects=[]
    for (group,unit,route),d in factorial.groupby(['group','unit','route'],sort=True):
        wide=d.pivot(index='patient',columns=['marker_choice','stage'],values='F1_T').sort_index()
        assert wide.shape==(4,4) and not wide.isna().any().any()
        f0=wide[('fixed_historical_context','initial')];f1=wide[('fixed_historical_context','final090')]
        s0=wide[('training_patient_selected_context','initial')];s1=wide[('training_patient_selected_context','final090')]
        for effect,v in [('DL_given_fixed_marker',f1-f0),('DL_given_selected_marker',s1-s0),
                         ('marker_selection_without_DL',s0-f0),('marker_selection_with_DL',s1-f1),
                         ('marker_DL_interaction',(s1-s0)-(f1-f0))]:
            interactions.append(dict(group=group,unit=unit,route=route,effect=effect,**paired(v)))
            for patient,value in v.items():patient_effects.append(dict(group=group,unit=unit,route=route,effect=effect,patient=patient,delta=float(value)))
    pd.DataFrame(interactions).to_csv(dest/'marker_DL_paired_effects.csv',index=False)
    pd.DataFrame(patient_effects).to_csv(dest/'marker_DL_patient_differences.csv',index=False)
    # New seed and retention conditions are separate from the reused native grid.
    preparations=json.loads((PTC/'selection/preparation_tasks.json').read_text())
    mlp_tasks=json.loads((PTC/'selection/MLP_tasks.json').read_text())
    ledger,execution_summary=terminal_execution_tables(mlp_tasks,preparations)
    ledger.to_csv(dest/'terminal_execution_ledger.csv.gz',index=False)
    execution_summary.to_csv(dest/'terminal_execution_summary.csv',index=False)
    execution_overview=dict(n_requested_conditions=len(ledger),
        n_fresh_training=int(ledger.training_executed_in_this_condition.sum()),
        n_cached_terminal_reuse=int(ledger.identical_result_reused.sum()),
        n_invalid_terminal=int((~ledger.terminal_valid).sum()),
        dl_status_counts={k:int(v) for k,v in ledger.dl_status.value_counts().items()},
        condition_action_counts={k:int(v) for k,v in ledger.condition_action.value_counts().items()},
        counting_unit='requested follow-up terminal condition, not patient, metric stage, or unique model',
        scope='50 MLP task groups and 22 new control units; archived seed42 MLP parity pilots are included once within the 50 groups; separate scoring-only parity pilots and reused original grid are excluded',
        training_executed_note='The raw training flag describes the cached result; fresh training additionally requires identical_result_reused=False')
    write_json(dest/'terminal_execution_overview.json',execution_overview)
    sources={};seed_frames=[];retention_frames=[];retention_states=[];figure_rows=[]
    for cfg in mlp_tasks:
        d=Path(cfg['dest']);assert checked(d);sources[str(d)]=sha(d/'manifest.json')
        f=pd.read_csv(d/'metrics.csv.gz',keep_default_na=False);f=primary(f);f=f[f.stage=='final090'].copy()
        f['seed_family']='MLP';f['seed']=cfg['seed'];f['base_unit']=f.unit;seed_frames.append(f)
    for cfg in preparations:
        d=Path(cfg['dest']);e=d/'evaluation';assert checked(e);sources[str(e)]=sha(e/'manifest.json')
        f=pd.read_csv(e/'metrics.csv.gz',keep_default_na=False)
        for route in CONTROL_ROUTES:
            assert (d/route/'clustering_and_terminal.png').exists()
            figure_rows.append(dict(**cfg,route=route,png=str((d/route/'clustering_and_terminal.png').relative_to(OUT)),
                pdf=str((d/route/'clustering_and_terminal.pdf').relative_to(OUT))))
        if cfg['kind']=='retention':
            retention_frames.append(f);retention_states.append(pd.read_csv(e/'terminal_statuses.csv.gz',keep_default_na=False))
        else:
            f=primary(f);f=f[f.stage=='final090'].copy();f['seed_family']='representation';f['seed']=cfg['seed']
            f['base_unit']=f"PTC_{cfg['group']}_CCA{cfg['budget']}";seed_frames.append(f)
    # Existing default42 is the fifth representation seed, with no refit silently substituted.
    for group in ANCHORS:
        for budget in ['2000','5000']:
            f=finals[(finals.group==group)&(finals.unit==f'PTC_{group}_CCA{budget}')&finals.route.isin(CONTROL_ROUTES)].copy()
            f['seed_family']='representation';f['seed']=42;f['base_unit']=f.unit;seed_frames.append(f)
    seeds=pd.concat(seed_frames,ignore_index=True);selected_seeds=[]
    for (family,group,base,route,seed),d in seeds.groupby(['seed_family','group','base_unit','route','seed'],sort=True):
        a=ANCHORS[group]
        fixed=d[(d.library==a['library'])&(d.cutoff==a['cutoff'])].copy()
        assert len(fixed)==4,(family,group,base,route,seed,len(fixed))
        fixed['marker_choice']='fixed_historical_context';selected_seeds.append(fixed)
        if base=='PTC_archived_CCA2000':continue
        c=choices[(choices.group==group)&(choices.unit==base)&(choices.route==route)]
        selected=d.merge(c[['held_patient','arm_id']],left_on=['patient','arm_id'],right_on=['held_patient','arm_id'],validate='one_to_one')
        assert len(selected)==4
        selected['marker_choice']='training_patient_selected_context';selected_seeds.append(selected)
    seeds=pd.concat(selected_seeds,ignore_index=True)
    seeds.to_csv(dest/'seed_patient_metrics.csv',index=False)
    seed_summary=seeds.groupby(['seed_family','group','base_unit','route','marker_choice','seed'])[metric_columns].mean().reset_index()
    seed_summary.to_csv(dest/'seed_patient_mean_summary.csv',index=False)
    retained=pd.concat(retention_frames,ignore_index=True);rstates=pd.concat(retention_states,ignore_index=True)
    retained.to_csv(dest/'marker_retention_all_metrics.csv.gz',index=False);rstates.to_csv(dest/'marker_retention_terminal_states.csv.gz',index=False)
    mechanisms=[]
    for group,a in ANCHORS.items():
        for budget in ['2000','5000','all']:
            oldunit=f'PTC_{group}_GEOMETRY{budget}_FIXED_CCAall_DL2000'
            newunit=f'PTC_{group}_geometry{budget}_uniform_marker_retention'
            oldstate=pd.read_csv(PTC/'existing_grid_evaluation'/oldunit/'terminal_statuses.csv.gz',keep_default_na=False)
            for route in CONTROL_ROUTES:
                for name,df,st,unit in [('fixed2000_scoring',existing,oldstate,oldunit),('uniform_all17_marker_retention',retained,rstates,newunit)]:
                    take=primary(df);take=take[(take.unit==unit)&(take.route==route)&(take.library==a['library'])&
                        (take.cutoff==a['cutoff'])&(take.stage=='final090')]
                    stat=st[(st.unit==unit)&(st.route==route)&(st.library==a['library'])&(st.cutoff==a['cutoff'])&(st.patient!='ALL')]
                    joined=take.merge(stat[['patient','n_strict_T_seeds','n_known']],on='patient',validate='one_to_one')
                    assert len(joined)==4
                    for row in joined.itertuples():mechanisms.append(dict(group=group,budget=budget,route=route,scoring_rule=name,
                        patient=row.patient,F1_T=row.F1_T,coverage=row.coverage,TCR_positive_recall=row.TCR_positive_recall,
                        n_strict_T_seeds=row.n_strict_T_seeds,n_known=row.n_known,dl_status=row.dl_status))
    mechanism=pd.DataFrame(mechanisms);mechanism.to_csv(dest/'marker_retention_matched_patient.csv',index=False)
    mechanism.groupby(['group','budget','route','scoring_rule'])[['F1_T','coverage','TCR_positive_recall','n_strict_T_seeds']].mean().reset_index().to_csv(dest/'marker_retention_matched_summary.csv',index=False)
    for (group,budget,route),d in mechanism.groupby(['group','budget','route']):
        wide=d.pivot(index='patient',columns='scoring_rule',values='F1_T').sort_index()
        contrasts.append(dict(group=group,budget=budget,route=route,contrast='uniform_retention_minus_fixed2000_scoring',metric='F1_T',
            **paired(wide.uniform_all17_marker_retention-wide.fixed2000_scoring)))
        for patient,row in wide.iterrows():
            contrast_differences.append(dict(group=group,budget=budget,route=route,patient=patient,
                contrast='uniform_retention_minus_fixed2000_scoring',metric='F1_T',
                delta=row.uniform_all17_marker_retention-row.fixed2000_scoring))
    comparisons=pd.DataFrame(contrasts);comparisons['p_Holm_within_contrast_family']=1.
    for _,family in comparisons.groupby('contrast'):
        ix=family.sort_values('exact_sign_flip_p',kind='stable').index
        comparisons.loc[ix,'p_Holm_within_contrast_family']=np.minimum(1,np.maximum.accumulate(comparisons.loc[ix,'exact_sign_flip_p'].to_numpy()*np.arange(len(ix),0,-1)))
    comparisons.to_csv(dest/'paired_patient_contrasts.csv',index=False)
    pd.DataFrame(contrast_differences).to_csv(dest/'paired_patient_differences.csv',index=False)
    pd.DataFrame(figure_rows).to_csv(dest/'all_new_clustering_figures.csv',index=False)
    # Do not discard historical paper endpoint or silently replace its non-T-positive F1.
    paper=[]
    for group,a in ANCHORS.items():
        paper.append(existing[(existing.group==group)&(existing.unit=='PTC_archived_CCA2000')&
            (existing.route==a['route'])&(existing.library==a['library'])&(existing.cutoff==a['cutoff'])&
            (existing.truth_definition=='paper_original')&(existing.endpoint=='paper_broad_T_compatibility')])
    pd.concat(paper,ignore_index=True).to_csv(dest/'original_paper_endpoint_reconciliation.csv',index=False)
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    colors=['#0072B2','#D55E00','#009E73','#CC79A7'];order=['500','1000','2000','3000','5000','all']
    for i,group in enumerate(ANCHORS):
        for j,mode in enumerate(['fixed_historical_context','training_patient_selected_context']):
            ax=axes[i,j];d=core_summary[(core_summary.group==group)&(core_summary.marker_choice==mode)]
            for (route,view),color in zip(d.groupby('route',sort=False),colors):
                ax.plot(range(6),view.set_index('budget').F1_T.reindex(order),marker='o',label=route,color=color)
            ax.set_xticks(range(6),['500','1k','2k','3k','5k','All']);ax.set_ylim(0,1)
            ax.set(title=group+' / '+mode.replace('_',' '),xlabel='Fresh within-group CCA feature budget',ylabel='Patient-mean strict-T / productive-TCR F1')
    axes[0,0].legend(fontsize=7,frameon=False)
    for ext in ['png','pdf']:fig.savefig(dest/f'PTC_workflow_comparison.{ext}',dpi=210,bbox_inches='tight')
    plt.close(fig)
    for family in ['MLP','representation']:
        fig,axes=plt.subplots(1,2,figsize=(13,5),layout='constrained')
        for ax,group in zip(axes,ANCHORS):
            view=seed_summary[(seed_summary.seed_family==family)&(seed_summary.group==group)]
            for (base,route,mode),d in view.groupby(['base_unit','route','marker_choice'],sort=True):
                label=base.replace('PTC_'+group+'_','')+'/'+route+('/CV marker' if mode.startswith('training') else '/fixed marker')
                ax.plot(d.seed,d.F1_T,'o-',linewidth=.8,markersize=3,label=label)
            ax.set(title=group,xlabel=family+' seed',ylabel='Four-patient mean strict-T / TCR F1',ylim=(0,1))
            ax.legend(fontsize=5,frameon=False,loc='best')
        fig.suptitle('Seeds are algorithm repeats; four patients remain the biological units')
        for ext in ['png','pdf']:fig.savefig(dest/f'PTC_{family}_seed_stability.{ext}',dpi=210,bbox_inches='tight')
        plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for i,group in enumerate(ANCHORS):
        for j,route in enumerate(CONTROL_ROUTES):
            ax=axes[i,j];d=mechanism[(mechanism.group==group)&(mechanism.route==route)]
            for rule,view in d.groupby('scoring_rule',sort=False):
                avg=view.groupby('budget').F1_T.mean().reindex(['2000','5000','all'])
                ax.plot(range(3),avg,'o-',label=rule)
            ax.set_xticks(range(3),['2k','5k','All']);ax.set_ylim(0,1)
            ax.set(title=group+' / '+route,xlabel='Geometry genes; same CCAall fit and DL2000',ylabel='Patient-mean strict-T / TCR F1')
    axes[0,0].legend(fontsize=7,frameon=False)
    for ext in ['png','pdf']:fig.savefig(dest/f'PTC_marker_retention_mechanism.{ext}',dpi=210,bbox_inches='tight')
    plt.close(fig)
    write_json(dest/'manifest.json',dict(status='completed',n_patients=4,old_units_reused=30,
        new_control_units=22,new_clustering_conditions=44,MLP_task_groups=50,
        terminal_execution=execution_overview,
        reference_weights_recovered=False,default_parity=sha(PTC/'verification/default_parity.json'),
        seed_and_control_sources=sources,selection_manifest_sha256=sha(PTC/'selection/manifest.json'),
        selection_is_transductive_patient_label_holdout=True,bootstrap_is_conditional_on_frozen_predictions=True,
        TCR_absence_is_not_established_nonT=True,job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest)

if __name__=='__main__':run()
