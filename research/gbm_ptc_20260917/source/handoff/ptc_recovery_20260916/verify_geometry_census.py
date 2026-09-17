"""Independent feature masks, cell order, source and physical partition census."""
import json
import os
import time
from pathlib import Path
from ptc_common import BASE,GROUPS,require_slurm,sha,utc,write_json,task_list,geometry_dir

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    tasks=task_list()
    while any(not (geometry_dir(t)/'GEOMETRY_COMPLETE').exists() for t in tasks):
        time.sleep(30)
    dest=BASE/'verification/geometry_census';dest.mkdir(parents=True,exist_ok=True)
    rows=[];partitions=[]
    source_cache={}
    for t in tasks:
        d=geometry_dir(t);m=json.loads((d/'geometry_manifest.json').read_text())
        assert (d/'GEOMETRY_COMPLETE').read_text().strip()==sha(d/'geometry_manifest.json')
        assert m['task']==t and m['status']=='completed' and not m['annotation_fields_in_fit']
        p=BASE/'samples'/t['sample']
        for name,h in m['inputs'].items():assert sha(p/name)==h
        for path,h in m['sources'].items():
            if path not in source_cache:source_cache[path]=sha(path)
            assert source_cache[path]==h
        genes=(p/'geometry_genes.txt').read_text().splitlines()
        rank=(p/'HVG_rank_top5000.txt').read_text().splitlines()
        index=np.load(d/'feature_indices.npy',allow_pickle=False)
        desired=genes if t['feature']=='all' else rank[:int(t['feature'][3:])]
        assert list(np.asarray(genes)[index])==desired and len(set(index))==len(index)
        assert m['geometry_feature_width']==len(index)
        cells=pd.read_csv(p/'cells.csv',index_col=0).index.astype(str).to_list()
        assert (d/'cells.txt').read_text().splitlines()==cells and m['n_cells']==len(cells)
        if t['dr']!='none':
            z=np.load(d/'embedding.npy',mmap_mode='r',allow_pickle=False)
            assert z.shape==(len(cells),t['dim']) and np.isfinite(z).all()
        for method in t['clusterers']:
            part=d/method;cm=json.loads((part/'cluster_manifest.json').read_text())
            assert (part/'CLUSTER_COMPLETE').read_text().strip()==sha(part/'cluster_manifest.json')
            assert cm['clusters_sha256']==sha(part/'clusters.csv')
            cl=pd.read_csv(part/'clusters.csv',dtype=str)
            assert cl.cell_id.to_list()==cells and cl.cluster.notna().all()
            partitions.append(dict(partition_path=str(part.relative_to(BASE)),n_cells=len(cells),
              n_clusters=cl.cluster.nunique(),precision_recovery=bool(cm.get('precision_recovery',False)),
              cluster_sha256=sha(part/'clusters.csv')))
        rows.append(dict(**t,actual_geometry_genes=len(index),actual_DL_genes=2000,
            n_cells=len(cells),feature_indices_sha256=sha(d/'feature_indices.npy'),
            geometry_manifest_sha256=sha(d/'geometry_manifest.json')))
    for top in ['pooled','matched_samples']:
        for manifest in sorted((BASE/top).glob('**/cluster_manifest.json')):
            part=manifest.parent;cm=json.loads(manifest.read_text())
            assert (part/'CLUSTER_COMPLETE').read_text().strip()==sha(manifest)
            assert cm['clusters_sha256']==sha(part/'clusters.csv')
            cl=pd.read_csv(part/'clusters.csv',dtype=str)
            assert cl.cell_id.is_unique and cl.cluster.notna().all()
            partitions.append(dict(partition_path=str(part.relative_to(BASE)),n_cells=len(cl),
              n_clusters=cl.cluster.nunique(),precision_recovery=False,cluster_sha256=sha(part/'clusters.csv')))
    assert len(rows)==496 and len(partitions)==1288
    assert len({x['partition_path'] for x in partitions})==1288
    pd.DataFrame(rows).to_csv(dest/'actual_features_and_geometry.csv',index=False)
    pd.DataFrame(partitions).to_csv(dest/'physical_partitions.csv',index=False)
    write_json(dest/'manifest.json',dict(status='complete',n_geometries=496,n_physical_partitions=1288,
      n_precision_recoveries=sum(x['precision_recovery'] for x in partitions),sources_checked=source_cache,
      outputs={p.name:sha(p) for p in dest.glob('*.csv')},job=os.environ['SLURM_JOB_ID'],
      source_sha256=sha(Path(__file__)),completed_at=utc()))
    (dest/'COMPLETE').write_text(sha(dest/'manifest.json')+'\n')
    print('All 496 feature masks and 1288 physical partitions independently checked',flush=True)

if __name__=='__main__':run()
