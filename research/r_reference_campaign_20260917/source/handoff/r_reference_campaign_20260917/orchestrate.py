"""Metadata-only SLURM dispatcher. Never loads scientific matrices or changes algorithms."""
import fcntl,json,os,subprocess,time
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
CODE=ROOT/'handoff/r_reference_campaign_20260917'
OUT=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'
ROUTES=['PCA30_SNN','PCA30_HDBSCAN_R','UMAP2_SNN','UMAP2_HDBSCAN_R']

def write(path,value):
    temp=path.with_name(path.name+'.part');temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)

def submit(options,script,args=()):
    cmd=['sbatch','--parsable',*options,str(CODE/script),*map(str,args)]
    result=subprocess.run(cmd,cwd=ROOT,text=True,capture_output=True)
    if result.returncode:raise RuntimeError(result.stderr)
    job=result.stdout.strip().split(';')[0]
    with (OUT/'submission_events.jsonl').open('a') as f:f.write(json.dumps({'time':time.time(),'job':job,'command':cmd})+'\n')
    print('SUBMITTED',job,script,list(args),flush=True);return job

def tick():
    OUT.mkdir(exist_ok=True)
    with (OUT/'dispatcher.lock').open('w') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        path=OUT/'dispatch_state.json'
        state=json.loads(path.read_text()) if path.exists() else {'units':{},'errors':[]}
        rows=subprocess.check_output(['squeue','-u','yimin','-h','-r','-o','%i|%j|%T'],text=True).splitlines()
        activeids={r.split('|')[0].split('_')[0] for r in rows}
        current_terminal_arrays={r.split('|')[0].split('_')[0] for r in rows if '|Rref_terminal|' in r}
        active_scores=sum('|Rref_score|' in r for r in rows)
        active_plots=sum('|Rref_plot|' in r for r in rows)
        active_audits=sum('|Rref_artifact_audit|' in r for r in rows)
        active_abstention=sum('|Rref_PTC_abstention|' in r for r in rows)
        queued=len(rows)
        prepared=list((OUT/'benchmark').glob('*/reference_CCA2000/PREPARED'))
        prepared+=list((OUT/'PTC_ablation').glob('*/PREPARED'))
        def priority(flag):
            p=str(flag)
            return (0 if '/PTC_ablation/' in p else (2 if '/HCL__' in p else 1),p)
        for flag in sorted(prepared,key=priority):
            prep=flag.parent
            m=json.loads((prep/'prepare_manifest.json').read_text())
            unit=m['unit'];u=state['units'].setdefault(unit,{})
            u['prepared']=str(prep)
            if not (prep/'GRID_SCORED').exists():
                # The initial GBM score pilot is registered externally.
                jid=u.get('score_job')
                if jid and jid not in activeids:
                    u['needs_score_job_review']=True
                    continue
                if not jid and active_scores<16 and queued<850:
                    if unit.endswith('CCAall'):mem,limit,partition='128G','03:00:00','nextgen'
                    elif '_GEOMETRY' in unit:mem,limit,partition='64G','03:00:00','nextgen'
                    elif unit.startswith('PTC_'):
                        mem,limit,partition=('96G','03:00:00','nextgen') if m['n_cells']>60000 else ('64G','03:00:00','nextgen')
                    else:
                        mem='12G' if m['n_cells']<5000 else ('24G' if m['n_cells']<15000 else ('32G' if m['n_cells']<30000 else '48G'))
                        limit,partition='02:00:00','nextgen'
                    u['score_job']=submit(['--mem='+mem,'--time='+limit,'--partition='+partition],'reference_score.sbatch',[unit,prep])
                    active_scores+=1;queued+=1
                continue
            u['grid_scored']=True
            expected=[]
            for route in ROUTES:
                source=prep/route
                sm=json.loads((source/'score_manifest.json').read_text())
                for aid in sm['arms']:
                    expected.append({'source':str(source),'arm_id':aid})
            missing=[t for t in expected if not (Path(t['source'])/'terminal'/t['arm_id']/'TERMINAL_COMPLETE').exists()]
            u['expected_terminal_conditions']=len(expected)
            u['completed_terminal_conditions']=len(expected)-len(missing)
            if not missing:
                u['terminal_complete']=True
                if (prep/'evaluation/COMPLETE').exists():
                    u['evaluated']=True
                    if m['dataset']=='PTC':
                        if (prep/'evaluation/ABSTENTION_COMPLETE').exists():u['abstention_audited']=True
                        elif not u.get('abstention_job') and active_abstention<6 and queued<850:
                            u['abstention_job']=submit([], 'ptc_abstention.sbatch',[prep]);queued+=1;active_abstention+=1
                        elif u.get('abstention_job') and u['abstention_job'] not in activeids:
                            u['needs_abstention_job_review']=True
                    if (prep/'verification/AUDIT_COMPLETE').exists():u['audited']=True
                    elif not u.get('audit_job') and active_audits<6 and queued<850:
                        u['audit_job']=submit([], 'audit_unit.sbatch',[prep,m['dataset']]);queued+=1;active_audits+=1
                    elif u.get('audit_job') and u['audit_job'] not in activeids:
                        u['needs_audit_job_review']=True
                    if (prep/'figures/FIGURES_COMPLETE').exists():u['plotted']=True
                    elif not u.get('plot_job') and active_plots<8 and queued<850:
                        u['plot_job']=submit([], 'plot_unit.sbatch',[prep,m['dataset']]);queued+=1;active_plots+=1
                    elif u.get('plot_job') and u['plot_job'] not in activeids:
                        u['needs_plot_job_review']=True
                elif not u.get('evaluation_job') and queued<850:
                    u['evaluation_job']=submit([], 'evaluate_unit.sbatch',[prep,m['dataset']]);queued+=1
                elif u.get('evaluation_job') and u['evaluation_job'] not in activeids:
                    u['needs_evaluation_job_review']=True
                continue
            jid=u.get('terminal_job')
            if jid and jid not in activeids:
                u['needs_terminal_job_review']=True
                continue
            if not jid and len(current_terminal_arrays)<5 and queued+len(missing)<850:
                taskfile=OUT/'task_lists'/(unit+'.terminal.json');taskfile.parent.mkdir(exist_ok=True)
                write(taskfile,missing)
                mem='32G' if unit.endswith('CCAall') else '12G'
                u['terminal_job']=submit([f'--array=0-{len(missing)-1}%12','--mem='+mem],
                                         'terminal.sbatch',[taskfile])
                current_terminal_arrays.add(u['terminal_job']);queued+=len(missing)
        state['updated_at_epoch']=time.time()
        expected=set((OUT/'benchmark_units.txt').read_text().splitlines())
        expected.update(s['unit'] for s in json.loads((CODE/'ptc_ablation_conditions.json').read_text()))
        expected.update(f'PTC_{g}_GEOMETRY{h}_FIXED_CCAall_DL2000' for g in ['NMT','TTU'] for h in ['500','1000','2000','3000','5000','all'])
        archived=OUT/'PTC_archived_CCA2000'
        complete=all(all(state['units'].get(unit,{}).get(key,False) for key in ['evaluated','audited','plotted']) for unit in expected)
        complete=complete and all(state['units'].get(unit,{}).get('abstention_audited',False) for unit in expected if unit.startswith('PTC_'))
        complete=complete and all((archived/p).exists() for p in ['evaluation/COMPLETE','verification/AUDIT_COMPLETE','figures/FIGURES_COMPLETE','evaluation/ABSTENTION_COMPLETE'])
        complete=complete and len(list((OUT/'verification/pre_guard_density').glob('*.json')))==16
        if complete and not state.get('finalization_job'):
            state['finalization_job']=submit([], 'finalize_campaign.sbatch')
            activeids.add(state['finalization_job'])
        if state.get('finalization_job') and state['finalization_job'] not in activeids and not (OUT/'summary/DELIVERY_RECEIPT.json').exists():
            state['finalization_needs_review']=True
        write(path,state)
        print('STATUS',len(prepared),'prepared',sum(u.get('grid_scored',False) for u in state['units'].values()),
              'scored',sum(u.get('terminal_complete',False) for u in state['units'].values()),'terminal units',flush=True)
        return state

if __name__=='__main__':
    import sys
    if '--loop' in sys.argv:
        assert os.environ.get('SLURM_JOB_ID')
        for _ in range(1440):
            try:tick()
            except Exception as exc:print('DISPATCH_ERROR',repr(exc),flush=True)
            time.sleep(45)
    else:tick()
