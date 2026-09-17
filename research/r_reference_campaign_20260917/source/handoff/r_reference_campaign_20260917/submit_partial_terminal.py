"""Submit metadata-only DL arrays for immutable, completed scoring routes."""
import argparse
import fcntl
import hashlib
import json
import subprocess
from pathlib import Path
from orchestrate import OUT, ROUTES, submit, write


def run():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('units',nargs='+')
    parser.add_argument('--concurrency',type=int,default=12)
    args=parser.parse_args()
    assert 1<=args.concurrency<=24
    assert subprocess.check_output(['id','-un'],text=True).strip()=='yimin'
    with (OUT/'dispatcher.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        path=OUT/'dispatch_state.json'
        state=json.loads(path.read_text())
        rows=subprocess.check_output(['squeue','-u','yimin','-h','-r','-o','%i|%j|%T'],text=True).splitlines()
        activeids={r.split('|')[0].split('_')[0] for r in rows}
        arrays={r.split('|')[0].split('_')[0] for r in rows if '|Rref_terminal|' in r}
        queued=len(rows)
        for unit in args.units:
            u=state['units'][unit]
            prep=Path(u['prepared'])
            if (prep/'GRID_SCORED').exists() or u.get('terminal_job'):
                print('Full-grid dispatcher already owns',unit,flush=True);continue
            if any(job in activeids for job in u.get('partial_terminal_jobs',[])):
                print('Existing partial array still active',unit,flush=True);continue
            tasks=[]
            for route in ROUTES:
                source=prep/route
                if not (source/'SCORE_COMPLETE').exists():continue
                manifest=source/'score_manifest.json'
                assert hashlib.sha256(manifest.read_bytes()).hexdigest()==(source/'SCORE_COMPLETE').read_text().strip()
                for aid in json.loads(manifest.read_text())['arms']:
                    if not (source/'terminal'/aid/'TERMINAL_COMPLETE').exists():
                        tasks.append({'source':str(source),'arm_id':aid})
            if not tasks:
                print('No unassigned completed-route arms',unit,flush=True);continue
            if len(arrays)>=5 or queued+len(tasks)>=850:
                print('Wait for campaign queue capacity',unit,flush=True);continue
            assignment_hash=hashlib.sha256(json.dumps(tasks,sort_keys=True).encode()).hexdigest()
            taskfile=OUT/'task_lists'/(unit+'.partial_terminal.'+assignment_hash[:16]+'.json')
            if taskfile.exists():assert json.loads(taskfile.read_text())==tasks
            else:write(taskfile,tasks)
            mem='32G' if unit.endswith('CCAall') else '12G'
            job=submit([f'--array=0-{len(tasks)-1}%{args.concurrency}','--mem='+mem],
                       'terminal.sbatch',[taskfile])
            u.setdefault('partial_terminal_jobs',[]).append(job)
            u.setdefault('partial_terminal_taskfiles',{})[job]=str(taskfile)
            write(path,state)
            arrays.add(job);activeids.add(job);queued+=len(tasks)
            print('PARTIAL_TERMINAL',unit,job,len(tasks),'arms',flush=True)


if __name__=='__main__':run()
