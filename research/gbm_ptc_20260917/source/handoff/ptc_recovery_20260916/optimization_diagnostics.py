"""Collect frozen optimization warnings and actual method parameters."""
import json
import os
from pathlib import Path
from ptc_common import BASE,task_list,geometry_dir,require_slurm,sha,utc,write_json
require_slurm()
import pandas as pd
assert (BASE/'verification/geometry_census/COMPLETE').exists()
rows=[]
for t in task_list():
    d=geometry_dir(t);m=json.loads((d/'geometry_manifest.json').read_text())
    stages=[('representation',t['dr'],m['representation'])]
    stages.extend(('clustering',method,c['fitting']) for method,c in m['partitions'].items())
    for stage,method,info in stages:
        params=info.get('params',{});warnings=info.get('warnings',[])
        it=info.get('n_iter');budget=params.get('max_iter')
        # sklearn TSNE stores the final zero-based loop index in n_iter_.
        # Other fitted methods here record an iteration count.
        zero_based=stage=='representation' and method=='TSNE'
        iterations=None if it is None else it+int(zero_based)
        rows.append(dict(task_index=t['task_id'],sample=t['sample'],group=t['group'],feature=t['feature'],
            dr=t['dr'],seed=t['seed'],input_space=t['input_space'],stage=stage,method=method,
            converged=info.get('converged'),n_iter=it,max_iter=budget,
            n_iter_is_zero_based=zero_based,iterations_completed=iterations,
            iteration_budget_reached=bool(iterations is not None and budget is not None and iterations>=budget),
            convergence_warning=any('converg' in str(w).lower() for w in warnings),
            warnings=json.dumps(warnings),seconds=info.get('seconds'),params=json.dumps(params,sort_keys=True)))
out=BASE/'verification/optimization_diagnostics';out.mkdir(exist_ok=True)
df=pd.DataFrame(rows);df.to_csv(out/'all_diagnostics.csv.gz',index=False)
df.groupby(['stage','method','feature'])[['iteration_budget_reached','convergence_warning']].agg(['sum','count']).to_csv(out/'summary.csv')
write_json(out/'manifest.json',dict(status='complete',n_rows=len(df),n_convergence_warnings=int(df.convergence_warning.sum()),
    n_iteration_budget_reached=int(df.iteration_budget_reached.sum()),source_sha256=sha(Path(__file__)),
    note='KMeans uses n_init10; GMM uses sklearn default n_init1 and max_iter100, diagonal covariance/reg1e-4. sklearn TSNE n_iter_ is the final zero-based loop index, so completed iterations are n_iter_+1. Reaching a budget is distinct from a recorded convergence warning. Parameters unchanged throughout.',
    job=os.environ['SLURM_JOB_ID'],completed_at=utc(),outputs={p.name:sha(p) for p in out.glob('*.csv*')}))
(out/'COMPLETE').write_text(sha(out/'manifest.json')+'\n')
print('Frozen numerical optimization metadata inventoried',flush=True)
