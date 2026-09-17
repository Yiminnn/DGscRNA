#!/usr/bin/env python3
"""Isolate array layout / loaded runtimes in legacy t-SNE numerical parity."""
import argparse
import json
import sys
from common import OUT,ROOT,require_slurm,write_json,sha,utc,runtime_record
require_slurm()
parser=argparse.ArgumentParser();parser.add_argument('--torch-runtime',type=int,required=True)
parser.add_argument('--run-label',default='initial')
args=parser.parse_args()
import numpy as np
if args.torch_runtime:
    import torch
    torch.set_num_threads(4);torch.set_num_interop_threads(1)
from threadpoolctl import threadpool_limits,threadpool_info
threadpool_limits(8)
from sklearn.metrics import adjusted_rand_score
sys.path.insert(0,str(ROOT/'handoff/g274_table4'))
from fit_core import FitConfig,reduce_matrix,cluster_matrix
sample='COLUMBIA100163'
old=ROOT/'results/g274_cohort/fit_v1'/sample
x=np.load(old/'prepare/X.npy',allow_pickle=False)
oldz=np.load(old/'dr/TSNE/Z.npy',allow_pickle=False)
dest=OUT/'verification/tsne_numerical_diagnostic'/args.run_label/f'torch_{args.torch_runtime}'
dest.mkdir(parents=True,exist_ok=True)
rows=[];config=FitConfig(gmm_covariance='diag')
for order in ['C','F']:
    matrix=np.array(x,order=order,copy=True);matrix.setflags(write=False)
    z,meta=reduce_matrix(matrix,'TSNE',config)
    np.save(dest/f'coordinates_{order}.npy',z,allow_pickle=False)
    row=dict(sample=sample,array_order=order,torch_runtime_loaded=bool(args.torch_runtime),
        same_numeric_input=bool(np.array_equal(x,matrix)),coordinates_exactly_equal=bool(np.array_equal(oldz,z)),
        coordinate_max_abs_difference=float(np.abs(oldz-z).max()),reducer=meta,cluster_parity={})
    for cl in ['KMeans','GMM','HDBSCAN']:
        labels,cm=cluster_matrix(z,cl,'TSNE',config)
        prior=np.load(old/'arms'/f'TSNE__{cl}'/'labels.npy',allow_pickle=False)
        row['cluster_parity'][cl]=dict(ari=float(adjusted_rand_score(prior,labels)),
                                        exact=bool(np.array_equal(prior,labels)))
    rows.append(row)
    print(json.dumps({k:v for k,v in row.items() if k!='reducer'}),flush=True)
write_json(dest/'diagnostic.json',dict(timestamp=utc(),rows=rows,source_scaled_hash=sha(old/'prepare/X.npy'),
    legacy_code_sha256=sha(ROOT/'handoff/g274_table4/fit_core.py'),diagnostic_code_sha256=sha(__file__),
    threadpools=threadpool_info(),**runtime_record()))
