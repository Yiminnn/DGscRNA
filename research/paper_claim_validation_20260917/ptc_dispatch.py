"""Metadata-only PTC continuation, strictly after verified full GBM delivery."""
import json
import subprocess
import time
from pathlib import Path
from common import OUT,OLD,sha,checked,write_json,utc
from dispatch import freeze,submit,states
from ptc_followup_common import PTC,ANCHORS,CONTROL_ROUTES,old_units,original_arm

def tick():
    path=PTC/'dispatch_state.json';PTC.mkdir(parents=True,exist_ok=True)
    state=json.loads(path.read_text()) if path.exists() else dict(jobs={},phase='waiting_for_verified_full_GBM_delivery')
    gate=OUT/'GBM_full_summary/GBM_FULL_DELIVERED.json'
    if not gate.exists():
        state['last_check']=utc();write_json(path,state);return False
    g=json.loads(gate.read_text())
    assert g['PTC_compute_may_start'] and g['receipt_sha256']==sha(OUT/'GBM_full_summary/DELIVERY_RECEIPT.json')
    active=states()
    if 'science_source' not in state:state['science_source']=str(freeze('PTC_followups'))
    source=Path(state['science_source']);counts={}
    for job,(name,status) in active.items():counts[name]=counts.get(name,0)+1
    def launch(key,script,args,name,mem,wall,cpus,limit,done):
        r=state['jobs'].setdefault(key,{})
        if done:r['status']='complete';return True
        if r.get('job') in active:return False
        if r.get('job') and time.time()-r.get('submitted',0)<180:return False
        if r.get('job'):
            lines=subprocess.check_output(['sacct','-j',r['job'],'-X','-n','-P','--format=JobID,State,ExitCode'],text=True).splitlines()
            final=[v.split('|') for v in lines if v.startswith(r['job']+'|')]
            if final and final[0][1].split()[0] in ['FAILED','OUT_OF_MEMORY','TIMEOUT','CANCELLED','COMPLETED','NODE_FAIL']:
                r.update(status='needs_review',accounting=final[0])
            return False
        if counts.get(name,0)>=limit:return False
        r.update(job=submit(source,script,args,[f'--job-name={name}',f'--cpus-per-task={cpus}',f'--mem={mem}',f'--time={wall}']),
            submitted=time.time(),source=str(source),status='submitted')
        counts[name]=counts.get(name,0)+1;write_json(path,state);return False
    def config_path(cfg):
        f=PTC/'configurations'/(cfg['name']+'.json')
        if f.exists():assert json.loads(f.read_text())==cfg
        else:write_json(f,cfg)
        return f
    # Saved-output evaluation is parallel and does not rerun completed fits.
    all_old=True
    for prep in old_units():
        all_old=launch('old_eval/'+prep.name,'ptc_evaluate.py',[prep],'claim_PTC_saved_eval','16G','02:00:00',2,8,
            checked(PTC/'existing_grid_evaluation'/prep.name)) and all_old
    if all_old:
        launch('selection','ptc_select.py',[],'claim_PTC_select','24G','02:00:00',2,1,checked(PTC/'selection'))
    # Native scoring and terminal pilots before releasing the new scientific arms.
    pilot_scoring=[];pilot_mlp=[]
    for group,budget in [('NMT','2000'),('TTU','5000')]:
        name=f'PTC_{group}_CCA{budget}_default_parity'
        cfg=dict(name=name,group=group,budget=budget,seed=42,kind='default_parity',dest=str(PTC/'default_parity'/name))
        f=config_path(cfg);d=Path(cfg['dest']);route=ANCHORS[group]['route']
        ready=launch('prepare/'+name,'ptc_control_job.py',['prepare',f],'claim_PTC_prepare','128G','04:00:00',4,2,
            checked(d,'prepare_manifest.json','PREPARED'))
        if ready:
            launch('score/'+name+'/'+route,'ptc_control_job.py',['score',f,route],'claim_PTC_score','128G','06:00:00',4,4,
                (d/route/'default_parity.json').exists())
        pilot_scoring.append(d/route/'default_parity.json')
        a=ANCHORS[group];s=OLD/'PTC_archived_CCA2000'/a['archive']
        name=f'PTC_archived_{group}_{a["route"]}_MLPseed42'
        cfg=dict(name=name,source=str(s),group=group,arm_ids=[original_arm(s,group)],seed=42,
            dest=str(PTC/'MLP_seeds'/name),kind='MLP_seed')
        f=config_path(cfg);d=Path(cfg['dest']);pilot_mlp.append(d)
        launch('MLP/'+name,'ptc_control_job.py',['MLP',f],'claim_PTC_MLP','32G','04:00:00',4,4,checked(d))
    pilot_ready=all(p.exists() and json.loads(p.read_text())['status']=='passed' for p in pilot_scoring) and all(checked(p) for p in pilot_mlp)
    proof=PTC/'verification/default_parity.json'
    if pilot_ready and not proof.exists():
        write_json(proof,dict(status='passed',scoring_records=[str(p) for p in pilot_scoring],MLP_records=[str(p) for p in pilot_mlp],
            source=str(source),score_sha256=sha(source/'ptc_score_R.R'),terminal_sha256=sha(source/'ptc_terminal.py'),completed_at=utc()))
    n_prepared=n_scores=n_terminal=n_controls=n_mlp=0
    scoring_limit=4
    marker_pilot=False
    if pilot_ready and checked(PTC/'selection'):
        preparations=json.loads((PTC/'selection/preparation_tasks.json').read_text())
        mlp=json.loads((PTC/'selection/MLP_tasks.json').read_text())
        pilot=next(c for c in preparations if c['kind']=='retention' and c['group']=='NMT' and c['budget']=='2000')
        marker_pilot_source=Path(pilot['dest'])/'UMAP2_HDBSCAN_R'
        if checked(marker_pilot_source,'score_manifest.json','SCORE_COMPLETE'):
            aid=original_arm(marker_pilot_source,'NMT')
            marker_pilot=checked(marker_pilot_source/'terminal'/aid,'terminal_manifest.json','TERMINAL_COMPLETE')
        # After preparation and the MLP series finish, their allocations are free.
        # Raise only scoring concurrency; retain every task's memory and threads.
        if marker_pilot and all(checked(Path(c['dest']),'prepare_manifest.json','PREPARED') for c in preparations) and all(checked(Path(c['dest'])) for c in mlp):
            scoring_limit=8
        contexts=json.loads((PTC/'selection/frozen_seed_contexts.json').read_text())
        for cfg in [pilot]+[c for c in preparations if c!=pilot]:
            if cfg!=pilot and not marker_pilot:continue
            f=config_path(cfg);d=Path(cfg['dest']);name=cfg['name']
            prepared=launch('prepare/'+name,'ptc_control_job.py',['prepare',f],'claim_PTC_prepare','128G','04:00:00',4,2,
                checked(d,'prepare_manifest.json','PREPARED'))
            n_prepared+=int(prepared)
            if not prepared:continue
            all_terminal=True
            for route in CONTROL_ROUTES:
                s=d/route
                scored=launch('score/'+name+'/'+route,'ptc_control_job.py',['score',f,route],'claim_PTC_score','128G','06:00:00',4,scoring_limit,
                    checked(s,'score_manifest.json','SCORE_COMPLETE') and (cfg['kind']!='retention' or (s/'retention_invariants.json').exists()))
                n_scores+=int(scored)
                if not scored:all_terminal=False;continue
                if cfg['kind']=='retention':arms=list(json.loads((s/'score_manifest.json').read_text())['arms'])
                else:arms=contexts[str(OLD/'PTC_ablation'/f"PTC_{cfg['group']}_CCA{cfg['budget']}"/route)]['arm_ids']
                # Run the fixed context first to expose the marker-retention resource pilot promptly.
                fixed=original_arm(s,cfg['group']);arms=[fixed]+[a for a in arms if a!=fixed]
                for aid in arms:
                    if cfg==pilot and not marker_pilot and (route!='UMAP2_HDBSCAN_R' or aid!=fixed):
                        all_terminal=False;continue
                    done=launch('terminal/'+name+'/'+route+'/'+aid,'ptc_control_job.py',['terminal',f,route,aid],
                        'claim_PTC_terminal','32G','04:00:00',4,24,checked(s/'terminal'/aid,'terminal_manifest.json','TERMINAL_COMPLETE'))
                    n_terminal+=int(done);all_terminal=done and all_terminal
            if all_terminal:
                done=launch('finish/'+name,'ptc_control_job.py',['finish',f],'claim_PTC_post','16G','01:00:00',2,4,checked(d/'evaluation'))
                n_controls+=int(done)
        for cfg in mlp:
            f=config_path(cfg)
            done=launch('MLP/'+cfg['name'],'ptc_control_job.py',['MLP',f],'claim_PTC_MLP','32G','06:00:00',4,4,checked(Path(cfg['dest'])))
            n_mlp+=int(done)
    state.update(phase='PTC_followup_execution',last_check=utc(),old_grid_evaluated=all_old,default_parity_passed=pilot_ready,
        marker_retention_resource_pilot_passed=marker_pilot,n_control_preparations=n_prepared,n_control_scores=n_scores,
        scoring_concurrency_limit=scoring_limit,
        n_control_terminal_arms=n_terminal,n_control_units_complete=n_controls,expected_control_units=22,
        n_MLP_tasks_complete=n_mlp,expected_MLP_tasks=50,
        review_required={k:r['accounting'] for k,r in state['jobs'].items() if r.get('status')=='needs_review'})
    all_done=n_controls==22 and n_mlp==50 and not state['review_required']
    if all_done:
        state['phase']='PTC_followup_fits_complete_reporting_pending'
        if not state.get('finalizer_job'):
            final_source=freeze('PTC_finalize')
            state['finalizer_job']=submit(final_source,'ptc_finalize.py',[],['--job-name=claim_PTC_finalize','--cpus-per-task=4','--mem=32G','--time=06:00:00'])
            state['finalizer_source']=str(final_source)
    write_json(path,state)
    print('PTC_PROGRESS',state['phase'],n_controls,'/22 controls',n_mlp,'/50 MLP',len(state['review_required']),'review',flush=True)
    return all_done

if __name__=='__main__':
    while not tick():time.sleep(45)
