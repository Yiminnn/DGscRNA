"""Annotation-free PTC geometry grid, using the frozen GBM transfer implementations."""
import gc
import json
import os
from pathlib import Path
import sys
import time
import traceback
from ptc_common import ROOT, BASE, require_slurm, selected_task, geometry_dir, sha, utc, write_json

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    from threadpoolctl import threadpool_limits
    sys.path.insert(0, str(ROOT/'handoff/g274_table4'))
    sys.path.insert(0, str(ROOT/'handoff/hvg_ptc_20260916'))
    from fit import representation, clusters
    t = selected_task()
    dest = geometry_dir(t)
    dest.mkdir(parents=True, exist_ok=True)
    prep = BASE/'samples'/t['sample']
    assert (prep/'PREPARED').exists()
    inputs = {n: sha(prep/n) for n in ['PREPARED', 'geometry_all.json', 'geometry_genes.txt',
                                      'HVG_rank_top5000.txt', 'cells.csv']}
    sources = {str(p): sha(p) for p in [Path(__file__), Path(__file__).with_name('ptc_common.py'),
               ROOT/'handoff/hvg_ptc_20260916/fit.py', ROOT/'handoff/g274_table4/fit_core.py',
               BASE/'protocol/single_sample_geometries.json']}
    if (dest/'GEOMETRY_COMPLETE').exists():
        m = json.loads((dest/'geometry_manifest.json').read_text())
        assert (dest/'GEOMETRY_COMPLETE').read_text().strip()==sha(dest/'geometry_manifest.json')
        assert m['sources']==sources and m['inputs']==inputs
        return
    start = time.perf_counter()
    info = json.loads((prep/'geometry_all.json').read_text())
    genes = (prep/'geometry_genes.txt').read_text().splitlines()
    cells = pd.read_csv(prep/'cells.csv', index_col=0).index
    assert len(cells)==info['n_cells'] and len(genes)==info['n_genes']
    mm = np.memmap(prep/'geometry_all.float32.bin', dtype='<f4', mode='r',
                   shape=(len(cells), len(genes)))
    if t['feature']=='all':
        idx = np.arange(len(genes))
    else:
        ranked = (prep/'HVG_rank_top5000.txt').read_text().splitlines()
        selected = ranked[:int(t['feature'][3:])]
        positions = {g:i for i,g in enumerate(genes)}
        idx = np.array([positions[g] for g in selected])
        assert len(idx)==int(t['feature'][3:])
    # A physical subset is passed to every reducer, including no-reduction controls.
    x = np.ascontiguousarray(mm[:,idx], dtype=np.float32)
    assert np.isfinite(x).all()
    np.save(dest/'feature_indices.npy', idx, allow_pickle=False)
    (dest/'cells.txt').write_text('\n'.join(cells)+'\n')
    manifest = dict(task=t, started_at=utc(), job=os.environ['SLURM_JOB_ID'],
        inputs=inputs, sources=sources, annotation_fields_in_fit=False,
        geometry_feature_width=x.shape[1], n_cells=len(x), partitions={}, status='running')
    write_json(dest/'geometry_manifest.json', manifest)
    with threadpool_limits(limits=min(4, int(os.environ.get('SLURM_CPUS_PER_TASK',4)))):
        if (dest/'embedding.npy').exists() and (dest/'embedding_info.json').exists():
            z = np.load(dest/'embedding.npy', allow_pickle=False)
            ri = json.loads((dest/'embedding_info.json').read_text())
        else:
            z, ri = representation(x, t)
            if t['dr']!='none':
                np.save(dest/'embedding.npy', z, allow_pickle=False)
                write_json(dest/'embedding_info.json', ri)
        for clusterer in t['clusterers']:
            part = dest/clusterer
            part.mkdir(exist_ok=True)
            if (part/'CLUSTER_COMPLETE').exists():
                cm = json.loads((part/'cluster_manifest.json').read_text())
                assert (part/'CLUSTER_COMPLETE').read_text().strip()==sha(part/'cluster_manifest.json')
                assert cm['task']==t
                manifest['partitions'][clusterer] = cm
                continue
            arm = dict(clusterer=clusterer, k=23, covariance='diag',
                       min_cluster_size=15, min_samples=15, families=[])
            cl, ci = clusters(z, arm, t['seed'])
            pd.DataFrame({'cell_id':cells, 'cluster':cl}).to_csv(part/'clusters.csv', index=False)
            cm = dict(task=t, clusterer=clusterer, arm=arm, fitting=ci,
                noise_label=-1 if clusterer=='HDBSCAN' else None,
                noise_scored_as_observed_cluster=True,
                note='PTC scoring retains noise as an observed group, as in historical R; Python HDBSCAN 15/15 is the GBM geometry transfer, separate from R minPts50',
                clusters_sha256=sha(part/'clusters.csv'))
            write_json(part/'cluster_manifest.json', cm)
            (part/'CLUSTER_COMPLETE').write_text(sha(part/'cluster_manifest.json')+'\n')
            manifest['partitions'][clusterer] = cm
            write_json(dest/'geometry_manifest.json', manifest)
        manifest.update(representation=ri, status='completed', elapsed_seconds=time.perf_counter()-start,
                        completed_at=utc())
    write_json(dest/'geometry_manifest.json', manifest)
    (dest/'GEOMETRY_COMPLETE').write_text(sha(dest/'geometry_manifest.json')+'\n')
    print(t['sample'],t['geometry_id'],'complete',flush=True)

if __name__=='__main__':
    try:
        run()
    except Exception:
        t=selected_task(); dest=geometry_dir(t);dest.mkdir(parents=True,exist_ok=True)
        write_json(dest/'geometry_error.json',dict(task=t,error=traceback.format_exc(),time=utc()))
        raise
