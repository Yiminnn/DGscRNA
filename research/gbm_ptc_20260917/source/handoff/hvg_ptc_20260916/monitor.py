#!/usr/bin/env python3
"""Operational metadata only; safe on the login node (no scientific arrays or metrics)."""
import json
from collections import Counter
import subprocess
from common import OUT, samples, utc, write_json

prep=Counter()
for s in samples():
    p=OUT/'prepared'/s/'manifest.json'
    prep[json.loads(p.read_text()).get('status','unknown') if p.exists() else 'missing']+=1
tasks=Counter();geometries=Counter();failures=[];elapsed=[]
for p in (OUT/'task_status').glob('*/*.json'):
    try:m=json.loads(p.read_text())
    except (OSError,json.JSONDecodeError):continue
    tasks[m.get('status','running')]+=1
    geometries.update(g['status'] for g in m.get('geometries',[]))
    failures.extend(dict(sample=m['sample'],feature=m['feature'],**f) for f in m.get('failures',[]))
    if m.get('status')=='completed':elapsed.append(m.get('elapsed_seconds',0))
queue=subprocess.run(['squeue','-u','yimin','-h','-o','%j|%T|%i|%R'],capture_output=True,text=True,check=True).stdout
queue_rows=[r.split('|',3) for r in queue.splitlines() if r.startswith('gbm_hvg_')]
state=dict(timestamp=utc(),preparation=dict(prep),sample_feature_tasks=dict(tasks),
           completed_geometry_metadata=dict(geometries),recorded_failures=failures,
           queue=queue_rows,complete_task_elapsed_seconds=elapsed,
           submitted={'7328008':{'offset':0,'range':'0-725','max_parallel':208},
                      '7328106':{'offset':726,'range':'0-199','max_parallel':48},
                      '7328647':{'offset':926,'range':'0-72,74-162','max_parallel':32}},
           not_yet_submitted_task_ids=[],
           structurally_unavailable_feature_tasks=[547,548,549,550,551,999])
write_json(OUT/'operational_status.json',state)
print(json.dumps({k:v for k,v in state.items() if k not in ['queue','complete_task_elapsed_seconds','recorded_failures']},indent=2))
print('Recorded geometry failures:',len(failures),'(full details saved in operational_status.json)')
print('Queue entries:',dict(Counter((r[0],r[1]) for r in queue_rows)))
