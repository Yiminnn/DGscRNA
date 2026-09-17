"""Original Harmony generation, isolated from installed Harmony2."""
import json
import os
import sys
import time
from ptc_common import BASE, RECOVERY, require_slurm, sha, utc, write_json

require_slurm()
import numpy as np
import pandas as pd
sys.path.insert(0,str(RECOVERY/'vendor_harmony'))
import harmonypy
assert harmonypy.__version__=='0.0.10'
group = sys.argv[1] if len(sys.argv)>1 else ['MTN','TUT'][int(os.environ['SLURM_ARRAY_TASK_ID'])]
dest=BASE/'prepared'/group
pca=pd.read_csv(dest/'NONE_PCA30.csv',index_col=0)
meta=pd.read_csv(dest/'cells.csv',index_col=0).loc[pca.index,['sample_id']]
assert pca.shape[1]==30 and meta.sample_id.nunique()==4 and not pca.isna().any().any()
start=time.perf_counter()
params=dict(theta=2,lamb=1,sigma=0.1,tau=0,block_size=0.05,max_iter_harmony=10,
            max_iter_kmeans=20,epsilon_cluster=1e-5,epsilon_harmony=1e-4,random_state=42)
model=harmonypy.run_harmony(pca.to_numpy(),meta,['sample_id'],**params)
z=np.asarray(model.Z_corr.T)
assert z.shape==pca.shape and np.isfinite(z).all()
pd.DataFrame(z,index=pca.index,columns=pca.columns).to_csv(dest/'HARMONY_PCA30.csv')
write_json(dest/'harmony_manifest.json',dict(status='completed',group=group,version=harmonypy.__version__,
    implementation=str(harmonypy.__file__),params=params,job=os.environ['SLURM_JOB_ID'],
    elapsed_seconds=time.perf_counter()-start,completed_at=utc(),
    input_sha256=sha(dest/'NONE_PCA30.csv'),output_sha256=sha(dest/'HARMONY_PCA30.csv'),
    objectives=list(map(float,model.objective_harmony)),kmeans_rounds=list(map(int,model.kmeans_rounds)),
    annotation_fields_in_fit=False,batch_covariate='sample_id'))
(dest/'HARMONY_COMPLETE').write_text(sha(dest/'harmony_manifest.json')+'\n')
