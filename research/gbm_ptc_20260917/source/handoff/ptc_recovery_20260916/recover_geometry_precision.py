"""Verify float32 GMM covariance failure, then recover only precision at fixed parameters."""
import json
import os
from pathlib import Path
import sys
import traceback
from ptc_common import BASE,ROOT,selected_task,geometry_dir,sha,utc,write_json,require_slurm

require_slurm()
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
sys.path.insert(0,str(ROOT/'handoff/g274_table4'))
sys.path.insert(0,str(ROOT/'handoff/hvg_ptc_20260916'))
from fit import clusters
t=selected_task();dest=geometry_dir(t)
m=json.loads((dest/'geometry_manifest.json').read_text())
assert m['task']==t and not (dest/'GEOMETRY_COMPLETE').exists()
z=np.load(dest/'embedding.npy',allow_pickle=False)
assert z.dtype==np.float32 and np.isfinite(z).all()
cells=(dest/'cells.txt').read_text().splitlines()
recovery=dict(task=t,job=os.environ['SLURM_JOB_ID'],started_at=utc(),
  source_sha256=sha(Path(__file__)),embedding_sha256=sha(dest/'embedding.npy'),
  fixed_parameters='Same K23, diagonal covariance,reg_covar1e-4,seed42, same fitted embedding',
  changed_parameter='Only input numerical precision for GMM:float32 to float64',
  reporting_policy='Flag all recovered results; also report primary comparisons excluding precision recoveries')
with threadpool_limits(limits=4):
  for method in t['clusterers']:
    part=dest/method;part.mkdir(exist_ok=True)
    if (part/'CLUSTER_COMPLETE').exists():
      m['partitions'][method]=json.loads((part/'cluster_manifest.json').read_text());continue
    arm=dict(clusterer=method,k=23,covariance='diag',min_cluster_size=15,min_samples=15,families=[])
    if method=='GMM':
      try:
        clusters(z,arm,t['seed'])
      except ValueError as error:
        assert 'ill-defined empirical covariance' in str(error)
        recovery['float32_failure_independently_reproduced']=str(error)
      else:
        raise RuntimeError('Original numerical failure did not reproduce; investigate before recovery')
      labels,info=clusters(z.astype(np.float64),arm,t['seed'])
      precision='float64'
    else:
      labels,info=clusters(z,arm,t['seed']);precision='float32'
    pd.DataFrame({'cell_id':cells,'cluster':labels}).to_csv(part/'clusters.csv',index=False)
    cm=dict(task=t,clusterer=method,arm=arm,fitting=info,numerical_precision=precision,
       precision_recovery=method=='GMM',noise_label=-1 if method=='HDBSCAN' else None,
       noise_scored_as_observed_cluster=True,clusters_sha256=sha(part/'clusters.csv'))
    write_json(part/'cluster_manifest.json',cm)
    (part/'CLUSTER_COMPLETE').write_text(sha(part/'cluster_manifest.json')+'\n')
    m['partitions'][method]=cm
recovery['completed_at']=utc();write_json(dest/'NUMERICAL_RECOVERY.json',recovery)
m.update(status='completed',completed_at=utc(),representation=json.loads((dest/'embedding_info.json').read_text()),
         numerical_precision_recovery=recovery)
write_json(dest/'geometry_manifest.json',m)
(dest/'GEOMETRY_COMPLETE').write_text(sha(dest/'geometry_manifest.json')+'\n')
print('Verified numerical precision recovery complete',t['task_id'])
