import os
from pathlib import Path
from ptc_common import BASE,sha,utc,require_slurm
import refine
require_slurm()
import torch
torch.set_num_threads(4);torch.set_num_interop_threads(1)
i=int(os.environ['SLURM_ARRAY_TASK_ID'])
sample=['MT-1','MT-2','N-1','N-2','TU-1','TU-2','T-1','T-2'][i]
group='MTN' if i<4 else 'TUT'
gd=BASE/'matched_samples'/sample
for space in ['PCA30','UMAP2']:
  for method in ['SNN','HDBSCAN_R']:
    context=dict(sample=sample,group=group,correction='per_sample_NONE',space=space,
      clusterer=method,family='matched_R_single',seed=42,dispatch_sha256=sha(Path(__file__)))
    refine.refine_score(gd/f'{space}_{method}'/'score_RNA',context,sha(Path(refine.__file__)))
(gd/'ANNOTATION_COMPLETE').write_text(utc()+'\n')
