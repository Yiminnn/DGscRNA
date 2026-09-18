"""Freeze PTC training-patient marker choices before new seed sensitivity fits."""
import json
import os
from pathlib import Path
from common import OLD,ROUTES,sha,checked,write_json,complete,utc
from ptc_followup_common import PTC,ANCHORS,SEEDS,CONTROL_ROUTES,require_ptc,old_units,original_arm

def run():
    require_ptc()
    import numpy as np
    import pandas as pd
    dest=PTC/'selection';dest.mkdir(parents=True,exist_ok=True)
    if checked(dest):return
    tables=[];sources={}
    for prep in old_units():
        d=PTC/'existing_grid_evaluation'/prep.name;assert checked(d)
        tables.append(pd.read_csv(d/'metrics.csv.gz',keep_default_na=False))
        sources[str(d/'manifest.json')]=sha(d/'manifest.json')
    frame=pd.concat(tables,ignore_index=True)
    primary=frame[(frame.truth_definition=='productive_TCR')&(frame.endpoint=='strict_T_name_rule')&(frame.patient!='ALL')]
    final=primary[primary.stage.eq('final090')]
    patients=['Patient1','Patient2','Patient3','Patient4']
    fold_rows=[];factorial=[];stage_choices=[];global_rows=[];seed_contexts={}
    for group,anchor in ANCHORS.items():
        units=[p.name for p in old_units() if p.name.startswith('PTC_'+group+'_')]
        units+=['PTC_archived_CCA2000','PTC_ALL8_CCA2000']
        for unit in units:
            for route in ROUTES:
                d=final[(final.group==group)&(final.unit==unit)&(final.route==route)]
                assert len(d)==51*4,(group,unit,route,len(d))
                assert not d.duplicated(['patient','arm_id']).any()
                fixed=d[(d.library==anchor['library'])&(d.cutoff==anchor['cutoff'])].arm_id.unique()
                assert len(fixed)==1;fixed=fixed[0]
                selected_ids={fixed}
                for held in patients:
                    train=d[d.patient!=held]
                    ranking=train.groupby('arm_id',sort=True).agg(objective=('F1_T','mean'),n=('patient','nunique')).reset_index()
                    assert ranking.n.eq(3).all()
                    winner=ranking.sort_values(['objective','arm_id'],ascending=[False,True],kind='stable').iloc[0]
                    chosen=str(winner.arm_id);selected_ids.add(chosen)
                    test=d[(d.patient==held)&(d.arm_id==chosen)].iloc[0]
                    fold_rows.append(dict(group=group,unit=unit,route=route,held_patient=held,arm_id=chosen,
                        library=test.library,cutoff=test.cutoff,training_mean_F1_T=float(winner.objective),
                        test_F1_T=float(test.F1_T),test_coverage=float(test.coverage),n_training_patients=3))
                    for mode,aid in [('fixed_historical_context',fixed),('training_patient_selected_context',chosen)]:
                        for stage in ['initial','final090']:
                            q=primary[(primary.group==group)&(primary.unit==unit)&(primary.route==route)&
                                (primary.patient==held)&(primary.arm_id==aid)&(primary.stage==stage)]
                            assert len(q)==1
                            row=q.iloc[0].to_dict();row.update(marker_choice=mode,selection_stage='final090',held_patient=held)
                            factorial.append(row)
                    # A separately identified sensitivity, never substituted into the matched 2x2.
                    for stage in ['initial','final090']:
                        z=primary[(primary.group==group)&(primary.unit==unit)&(primary.route==route)&
                            (primary.patient!=held)&(primary.stage==stage)]
                        ranked=z.groupby('arm_id',sort=True).F1_T.mean().sort_values(ascending=False,kind='stable')
                        aid=str(ranked.index[0])
                        q=primary[(primary.group==group)&(primary.unit==unit)&(primary.route==route)&
                            (primary.patient==held)&(primary.arm_id==aid)&(primary.stage==stage)].iloc[0]
                        stage_choices.append(dict(group=group,unit=unit,route=route,held_patient=held,stage=stage,
                            arm_id=aid,training_mean_F1_T=float(ranked.iloc[0]),test_F1_T=float(q.F1_T)))
                if unit in [f'PTC_{group}_CCA2000',f'PTC_{group}_CCA5000'] and route in CONTROL_ROUTES:
                    path=OLD/'PTC_ablation'/unit/route
                    seed_contexts[str(path)]=dict(group=group,unit=unit,route=route,arm_ids=sorted(selected_ids),
                        basis='Historical group context plus all four fold-specific seed42 training-patient choices')
        # Equal 24-route/budget opportunities: markers/cutoffs also selected on three patients.
        candidates=[]
        for budget in ['500','1000','2000','3000','5000','all']:
            for route in ROUTES:
                sub=final[(final.group==group)&(final.unit==f'PTC_{group}_CCA{budget}')&(final.route==route)].copy()
                sub['configuration_order']=len(candidates);candidates.append(sub)
        grid=pd.concat(candidates,ignore_index=True)
        for held in patients:
            scores=grid[grid.patient!=held].groupby(['configuration_order','unit','route','arm_id'],sort=True).F1_T.mean().reset_index()
            best=scores.sort_values(['F1_T','configuration_order','arm_id'],ascending=[False,True,True],kind='stable').iloc[0]
            q=grid[(grid.patient==held)&(grid.unit==best.unit)&(grid.route==best.route)&(grid.arm_id==best.arm_id)].iloc[0]
            global_rows.append(dict(group=group,held_patient=held,unit=best.unit,route=best.route,arm_id=best.arm_id,
                library=q.library,cutoff=q.cutoff,training_mean_F1_T=float(best.F1_T),test_F1_T=float(q.F1_T),
                test_coverage=float(q.coverage),test_TCR_recall=float(q.TCR_positive_recall)))
    pd.DataFrame(fold_rows).to_csv(dest/'marker_choices_by_patient.csv',index=False)
    pd.DataFrame(factorial).to_csv(dest/'marker_DL_factorial.csv.gz',index=False)
    pd.DataFrame(stage_choices).to_csv(dest/'stage_specific_selection_sensitivity.csv',index=False)
    pd.DataFrame(global_rows).to_csv(dest/'workflow_choices_by_patient.csv',index=False)
    frame.to_csv(dest/'all_existing_patient_metrics.csv.gz',index=False)
    write_json(dest/'frozen_seed_contexts.json',seed_contexts)
    # Frozen compute manifests. No held-out labels are available to any fit script.
    mlp=[]
    for path,context in seed_contexts.items():
        for seed in SEEDS:
            name=f"{context['unit']}_{context['route']}_MLPseed{seed}"
            mlp.append(dict(name=name,source=path,group=context['group'],arm_ids=context['arm_ids'],seed=seed,
                dest=str(PTC/'MLP_seeds'/name),kind='MLP_seed'))
    for group,anchor in ANCHORS.items():
        source=OLD/'PTC_archived_CCA2000'/anchor['archive']
        for seed in SEEDS:
            name=f'PTC_archived_{group}_{anchor["route"]}_MLPseed{seed}'
            mlp.append(dict(name=name,source=str(source),group=group,arm_ids=[original_arm(source,group)],seed=seed,
                dest=str(PTC/'MLP_seeds'/name),kind='MLP_seed'))
    preparations=[]
    for group in ANCHORS:
        for budget in ['2000','5000']:
            for seed in [0,1,2,3]:
                name=f'PTC_{group}_CCA{budget}_representation_seed{seed}'
                preparations.append(dict(name=name,group=group,budget=budget,seed=seed,kind='representation',
                    dest=str(PTC/'representation'/name)))
        for budget in ['2000','5000','all']:
            name=f'PTC_{group}_geometry{budget}_uniform_marker_retention'
            preparations.append(dict(name=name,group=group,budget=budget,seed=42,kind='retention',
                dest=str(PTC/'marker_retention'/name)))
    assert len(mlp)==50 and len(preparations)==22 and len(seed_contexts)==8
    write_json(dest/'MLP_tasks.json',mlp);write_json(dest/'preparation_tasks.json',preparations)
    write_json(dest/'manifest.json',dict(status='complete',objective='Patient-mean strict-T F1 versus productive high-confidence TCR detection',
        heldout='Both samples of each of four patients excluded from label-based selection; frozen full-group expression integration is transductive',
        evaluation_sources=sources,seed_contexts_sha256=sha(dest/'frozen_seed_contexts.json'),
        seed_results_used_for_selection=False,n_MLP_tasks=len(mlp),n_preparations=len(preparations),
        tie_breaking='Frozen HVG order then route order then arm_id; within a route, arm_id',
        job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest)

if __name__=='__main__':run()
