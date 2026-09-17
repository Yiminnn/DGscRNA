#!/usr/bin/env python3
"""Label-free geometric diagnostics in a common all-gene reference, after fit."""
import argparse
import json
import sys
from common import ROOT, OUT, require_slurm, sha, utc, write_json, samples


def run(sample, partial=False):
    require_slurm()
    import numpy as np
    import pandas as pd
    from sklearn.metrics import pairwise_distances
    from threadpoolctl import threadpool_limits
    sys.path.insert(0,str(ROOT/'handoff/g274_table4'))
    from metrics_core import sample_index_pairs, pair_distances, distance_residual_variance
    from resolve import structural_resolution
    threadpool_limits(8)
    prep=OUT/'prepared'/sample
    x=np.load(prep/'scaled_all.npy',mmap_mode='r',allow_pickle=False)
    out=OUT/'quality'/sample
    out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(20260916)
    cells=np.sort(rng.choice(len(x),min(2000,len(x)),replace=False))
    np.save(out/'sampled_cell_indices.npy',cells,allow_pickle=False)
    n=len(cells); k=20
    assert n>2*k
    d=pairwise_distances(x[cells],metric='euclidean',n_jobs=1)
    np.fill_diagonal(d,np.inf)
    ranks=np.argsort(d,axis=1,kind='stable')
    original_neighbors=ranks[:,:k]
    inverse=np.empty(ranks.shape,dtype=np.int32)
    inverse[np.arange(n)[:,None],ranks]=np.arange(1,n+1,dtype=np.int32)
    del d,ranks
    pairs=sample_index_pairs(len(x),max_pairs=20000,seed=20260916)
    np.save(out/'sampled_distance_pairs.npy',np.asarray(pairs,dtype=np.int32),allow_pickle=False)
    original_dist=pair_distances(x,pairs,batch_size=128)
    rows=[]; missing=[]
    gs=json.loads((OUT/'protocol/geometries.json').read_text())
    for g in gs:
        dest=OUT/'fits'/sample/g['geometry_id']
        if not (dest/'COMPLETE').exists():
            resolutions=[structural_resolution(sample,g,arm) if not (dest/arm['arm_id']/'COMPLETE').exists()
                         else {'status':'fit_complete'} for arm in g['arms']]
            if any(r is None for r in resolutions):
                missing.append(g['geometry_id']);continue
            if all(r['status']=='structural_preprocessing_failure' for r in resolutions):
                rows.append(dict(sample=sample,geometry_id=g['geometry_id'],feature=g['feature'],dr=g['dr'],
                    dim=g['dim'],input_space=g['input_space'],seed=g['seed'],status='structural_preprocessing_failure'))
                continue
        else:
            assert (dest/'COMPLETE').read_text().strip()==sha(dest/'manifest.json')
        em=json.loads((dest/'embedding_manifest.json').read_text())
        assert em['prepared_manifest_sha256']==sha(prep/'manifest.json')
        if g['dr']=='none':
            idx=np.load(prep/f'indices_{g["feature"]}.npy',allow_pickle=False)
            z=np.asarray(x[:,idx])
        else:
            assert sha(dest/'embedding.npy')==em['sha256']
            z=np.load(dest/'embedding.npy',allow_pickle=False)
        dz=pairwise_distances(z[cells],metric='euclidean',n_jobs=1)
        np.fill_diagonal(dz,np.inf)
        embedded_neighbors=np.argsort(dz,axis=1,kind='stable')[:,:k]
        penalty=np.maximum(inverse[np.arange(n)[:,None],embedded_neighbors]-k,0).sum(dtype=np.float64)
        trust=1-2*penalty/(n*k*(2*n-3*k-1))
        overlap=np.mean([len(set(a)&set(b))/k for a,b in zip(original_neighbors,embedded_neighbors)])
        rv=distance_residual_variance(original_dist,pair_distances(z,pairs,batch_size=128))
        rows.append(dict(sample=sample,geometry_id=g['geometry_id'],feature=g['feature'],dr=g['dr'],
            dim=g['dim'],input_space=g['input_space'],seed=g['seed'],neighbors=g['neighbors'],min_dist=g['min_dist'],status='completed',
            n_neighbor_cells=n,neighbor_k=k,trustworthiness_common_allgenes=float(trust),
            neighbor_overlap_common_allgenes=float(overlap),distance_rv_common_allgenes=rv['rv_e'],
            n_distance_pairs=len(pairs),reference='all filtered scaled genes',sampling_seed=20260916))
        pd.DataFrame(rows).to_csv(out/'geometry_quality.csv',index=False)
    meta=dict(sample=sample,status='completed' if not missing else 'incomplete',timestamp=utc(),
        n_expected=len(gs),n_complete=len(rows),missing=missing,reference='same all-gene scaled expression for all features',
        neighborhood_scope='same label-free random subset of at most 2000 cells; k=20',
        distance_scope='same 20000 uniformly sampled unordered pairs across full sample',
        outputs={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name not in ['manifest.json','COMPLETE']})
    write_json(out/'manifest.json',meta)
    if not missing:
        (out/'COMPLETE').write_text(sha(out/'manifest.json')+'\n')
    if missing and not partial:
        raise RuntimeError(f'{len(missing)} representations not complete')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=int);p.add_argument('--sample');p.add_argument('--partial',action='store_true')
    a=p.parse_args();run(a.sample or samples()[a.index],a.partial)
