#!/usr/bin/env python3
"""Recheck reproducible numerical fit errors with unchanged parameters; never rescue/tune."""
import json
import sys
import traceback
from pathlib import Path
from common import OUT,ROOT,require_slurm,sha,utc,write_json,runtime_record
require_slurm()
import numpy as np
import torch
from threadpoolctl import threadpool_limits
torch.set_num_threads(4);torch.set_num_interop_threads(1);threadpool_limits(8)
sys.path.insert(0,str(ROOT/'handoff/g274_table4'))
from fit import clusters

geometries={g['geometry_id']:g for g in json.loads((OUT/'protocol/geometries.json').read_text())}
checked=[]
for taskfile in (OUT/'task_status').glob('*/*.json'):
    task=json.loads(taskfile.read_text())
    if task.get('status')!='failed':continue
    sample=task['sample']
    for failed in task.get('failures',[]):
        g=geometries[failed['geometry_id']]
        for arm in g['arms']:
            ap=OUT/'fits'/sample/g['geometry_id']/arm['arm_id']
            mp=ap/'manifest.json'
            if not mp.exists() or (ap/'COMPLETE').exists():continue
            am=json.loads(mp.read_text())
            if am.get('status')!='failed' or 'ill-defined empirical covariance' not in am.get('error',''):continue
            assert arm['clusterer']=='GMM'
            dest=OUT/'numerical_failure_checks'/sample/g['geometry_id']/arm['arm_id']
            dest.mkdir(parents=True,exist_ok=True)
            if (dest/'verification.json').exists():
                prior=json.loads((dest/'verification.json').read_text())
                if prior['original_fit_manifest_sha256']==sha(mp):continue
            ep=ap.parent/'embedding.npy'
            em=json.loads((ap.parent/'embedding_manifest.json').read_text())
            assert sha(ep)==em['sha256']
            z=np.load(ep,allow_pickle=False)
            attempts=[]
            for i in range(3):
                try:
                    labels,info=clusters(z,arm,g['seed'])
                    np.save(dest/f'retry_{i}_labels.npy',labels,allow_pickle=False)
                    attempts.append(dict(status='succeeded',clustering=info))
                except ValueError as exc:
                    attempts.append(dict(status='failed',error=str(exc),traceback=traceback.format_exc()))
            result=dict(sample=sample,geometry_id=g['geometry_id'],arm_id=arm['arm_id'],timestamp=utc(),
                original_fit_manifest_sha256=sha(mp),embedding_sha256=sha(ep),seed=g['seed'],arm=arm,
                attempts=attempts,parameters_changed=False,precision=str(z.dtype),
                reproducible_covariance_failure=all(x['status']=='failed' and 'ill-defined empirical covariance' in x['error'] for x in attempts),
                **runtime_record())
            write_json(dest/'verification.json',result)
            checked.append({k:result[k] for k in ['sample','geometry_id','arm_id','reproducible_covariance_failure']})
            print(checked[-1],flush=True)
write_json(OUT/'numerical_failure_checks/latest_check.json',dict(timestamp=utc(),checked=checked))
