"""Resumable GBM control/comparator scheduling independent of core array execution."""
import argparse
import csv
import fcntl
import json
import time
from pathlib import Path
from common import OUT, checked, sha, write_json, utc
from dispatch import freeze, submit, states
from sccatch_resource_gate import covering_pilot, DIMENSIONS

def tick():
    with (OUT/'aux_dispatch.lock').open('a') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:
            print('AUX_WAITING_FOR_EXISTING_TICK',flush=True)
            return False
        path=OUT/'aux_dispatch_state.json'
        state=json.loads(path.read_text()) if path.exists() else dict(jobs={},phase='GBM_auxiliary_running')
        samples=list(csv.DictReader((OUT/'protocol/sample_order.csv').open()))
        active=states();source=freeze('auxiliary')
        n_active=sum(v[0] in ['claim_GBM_geometry','claim_GBM_scType','claim_GBM_SCINA','claim_GBM_scCATCH','claim_GBM_DL_control'] for v in active.values())
        core=json.loads((OUT/'dispatch_state.json').read_text())
        limit=48 if core.get('n_evaluated_plotted')==726 else 32
        type_ready=checked(OUT/'comparators/scType/TKU4163','audit_sensitivity_manifest.json','AUDIT_SENSITIVITY_COMPLETE')
        scina_ready=all(checked(OUT/f'comparators/SCINA/{s}/L{i:02d}') for s in ['TKU4163','NL022','SN040'] for i in range(16))
        catch_source=sha(source/'sccatch_R.R')
        catch_audit=OUT/'comparators/scCATCH/TKU4163/hvg2000/PCA30_SNN'
        catch_audited=checked(catch_audit,'audit_manifest.json','AUDIT_COMPLETE')
        if catch_audited:catch_audited=json.loads((catch_audit/'audit_manifest.json').read_text())['source_sha256']==catch_source
        catch_completed={}
        for sample in ['NL022','SN040']:
            pilot=OUT/f'comparators/scCATCH/{sample}/hvg2000/UMAP2_HDBSCAN_R'
            if checked(pilot,'cohort_manifest.json','COHORT_COMPLETE'):
                manifest=json.loads((pilot/'cohort_manifest.json').read_text())
                if manifest['source_sha256']==catch_source:catch_completed[sample]=manifest
        catch_ready=catch_audited and len(catch_completed)==2
        profile_dir=OUT/'verification/scCATCH_resource_profiles';catch_profiles={};catch_bounds={}
        if checked(profile_dir):
            profiles=json.loads((profile_dir/'manifest.json').read_text())
            if profiles['r_source_sha256']==catch_source and profiles['marker_manifest_sha256']==sha(OUT/'markers/scCATCH/manifest.json'):
                catch_profiles=profiles['profiles']
                catch_bounds={s:catch_profiles[s] for s,m in catch_completed.items()
                              if catch_profiles[s]['n_cells']==m['n_cells']}
        catch_task_gates={}
        dl_ready=(OUT/'verification/DL_control_default_parity.json').exists()
        pilots=(OUT/'protocol/pilot_samples.txt').read_text().split()
        tasks=[]
        for row in samples:
            sample=row['sample'];base=OUT/'GBM'/sample
            for budget in ['hvg2000','hvg5000','all']:
                prep=base/budget
                if (prep/'pipeline_manifest.json').exists() and checked(base/'hvg2000','prepare_manifest.json','PREPARED'):
                    tasks.append((f'geometry/{sample}/{budget}','geometry_only.py',[sample,budget],
                        'claim_GBM_geometry',prep/'evaluation_geometry_only_DL2000','00:50:00','24G',4))
            if type_ready and all((base/b/'pipeline_manifest.json').exists() for b in ['hvg500','hvg1000','hvg2000','hvg3000','hvg5000','all']):
                tasks.append((f'scType/{sample}','comparator_job.py',['scType',sample],
                    'claim_GBM_scType',OUT/'comparators/scType'/sample/'evaluation','01:00:00','24G',2))
            if all(checked(OUT/'comparators/SCINA'/sample/f'L{i:02d}') for i in range(16)):
                tasks.append((f'SCINA/{sample}','evaluate_comparator.py',['SCINA',sample],
                    'claim_GBM_SCINA',OUT/'comparators/SCINA'/sample/'evaluation','00:30:00','8G',2))
            catch_gate=None
            if catch_audited and (base/'hvg2000/UMAP2_HDBSCAN_R/SCORE_COMPLETE').exists():
                if catch_ready:catch_gate=dict(policy='all_three_native_pilots_complete')
                elif sample in catch_profiles:
                    profile=catch_profiles[sample];prep=base/'hvg2000'
                    if (profile['prepare_manifest_sha256']==sha(prep/'prepare_manifest.json') and
                        profile['score_manifest_sha256']==sha(prep/'UMAP2_HDBSCAN_R/score_manifest.json')):
                        covered=covering_pilot(profile,catch_bounds)
                        if covered:
                            catch_gate=dict(policy='bounded_by_completed_native_pilot',pilot=covered,
                                pilot_job=catch_completed[covered]['job'],
                                sample_dimensions={k:profile[k] for k in DIMENSIONS},
                                pilot_dimensions={k:catch_bounds[covered][k] for k in DIMENSIONS},
                                profile_manifest_sha256=sha(profile_dir/'manifest.json'))
            if catch_gate:
                # The largest native pairwise-DEG pilot already exceeded three
                # hours. Reserve its bounded recovery allowance for cohort jobs;
                # keep the audited implementation and memory request unchanged.
                tasks.append((f'scCATCH/{sample}','comparator_job.py',['scCATCH',sample],
                    'claim_GBM_scCATCH',OUT/'comparators/scCATCH'/sample/'evaluation','12:00:00','24G',2))
                catch_task_gates[sample]=catch_gate
            if dl_ready and sample in pilots:
                for budget in ['hvg2000','hvg5000']:
                    if not (base/budget/'pipeline_manifest.json').exists():continue
                    for config in ['original','model_seed0','model_seed1','model_seed2','model_seed3','width_small','width_large','epochs5','epochs20']:
                        tasks.append((f'DL_control/{sample}/{budget}/{config}','dl_controls.py',[sample,budget,config],
                            'claim_GBM_DL_control',OUT/'GBM_DL_controls'/sample/budget/config,'01:00:00','16G',4))
        for key,script,args,name,dest,wall,mem,cpus in tasks:
            record=state['jobs'].setdefault(key,{})
            if checked(dest):record['status']='complete';continue
            if record.get('job') in active:continue
            if record.get('job') and time.time()-record['submitted']<180:continue
            if record.get('job'):
                # No blind retries of scientific failures. The agent reviews the log,
                # preserves attempts, and explicitly clears the job for a corrected retry.
                import subprocess
                raw=subprocess.check_output(['sacct','-j',record['job'],'-X','-n','-P','--format=JobID,State,ExitCode'],text=True)
                final=[s.split('|') for s in raw.splitlines() if s.startswith(record['job']+'|')]
                if final and final[0][1].split()[0] in ['FAILED','TIMEOUT','OUT_OF_MEMORY','CANCELLED','COMPLETED','NODE_FAIL']:
                    record['status']='needs_review';record['accounting']=final[0]
                continue
            if n_active>=limit:continue
            jid=submit(source,script,args,[f'--job-name={name}',f'--cpus-per-task={cpus}',f'--mem={mem}',f'--time={wall}'])
            record.update(job=jid,source=str(source),submitted=time.time(),status='submitted');n_active+=1
            if name=='claim_GBM_scCATCH':record['resource_gate']=catch_task_gates[args[1]]
            write_json(path,state)
        complete_count=sum(r.get('status')=='complete' for r in state['jobs'].values())
        state.update(last_check=utc(),expected=780,ready_tasks=len(tasks),completed=complete_count,concurrency_limit=limit,
            scType_parity_gate=type_ready,SCINA_three_size_pilot_gate=scina_ready,
            scCATCH_three_size_gate=catch_ready,DL_default_parity_gate=dl_ready,
            scCATCH_resource_eligible_samples=len(catch_task_gates),
            review_required={k:r['accounting'] for k,r in state['jobs'].items() if r.get('status')=='needs_review'})
        if complete_count==780:state['phase']='GBM_geometry_three_marker_methods_MLP_controls_complete'
        write_json(path,state)
        print('AUX_PROGRESS',complete_count,'/780',len(state['review_required']),'review',flush=True)
        return complete_count==780

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--loop',action='store_true');a=p.parse_args()
    while True:
        if tick() or not a.loop:break
        time.sleep(45)
