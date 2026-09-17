import os
import json
from pathlib import Path
from ptc_common import BASE,sha,utc,require_slurm
import refine
require_slurm()
import torch
torch.set_num_threads(4);torch.set_num_interop_threads(1)
tasks=json.loads((BASE/'protocol/pooled_partition_tasks.json').read_text())
t=tasks[int(os.environ['SLURM_ARRAY_TASK_ID'])]
gd=BASE/'pooled'/t['group']/t['correction']
score=gd/f'{t["space"]}_{t["clusterer"]}'/f'score_{t["assay"]}'
refine.refine_score(score,{**t,'dispatch_sha256':sha(Path(__file__))},sha(Path(refine.__file__)))
matching=[d for d in tasks if d['group']==t['group'] and d['correction']==t['correction']]
if all((gd/f'{d["space"]}_{d["clusterer"]}'/f'score_{d["assay"]}'/'REFINEMENT_COMPLETE').exists() for d in matching):
    (gd/'ANNOTATION_COMPLETE').write_text(utc()+'\n')
print('Pooled partition terminal refinement complete',t['task_id'])
