"""Balanced shared-lineage batch mixing with RNA-state and local-neighbor preservation."""
from pathlib import Path
import json
import os
from ptc_common import BASE,GROUPS,require_slurm,sha,utc,write_json
require_slurm()
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits

def neighbors(x,k):
    k=min(k,len(x)-1)
    raw=NearestNeighbors(n_neighbors=k+1,metric='euclidean',n_jobs=4).fit(x).kneighbors(x,return_distance=False)
    return np.stack([a[a!=i][:k] for i,a in enumerate(raw)])

ref=pd.read_csv(BASE/'evaluation_reference/reference_cells.csv.gz',index_col=0)
dest=BASE/'batch_biology';dest.mkdir(exist_ok=True)
rows=[];availability=[];selected=[];per_sample=[]
rng=np.random.default_rng(42)
states=['Proliferation_mean_log1p','Interferon_mean_log1p','Stress_mean_log1p',
        'CD4_mean_log1p','CD8_mean_log1p','Treg_mean_log1p']
with threadpool_limits(limits=4):
 for group,samples in GROUPS.items():
  r=ref[ref.group.eq(group)]
  scopes={'productive_TCR_positive':r.TCR_cell_high_confidence_productive_TCR.astype(bool)}
  scopes.update({f'S2_broad_{name}':r.S2_terminal_broad.eq(name) for name in sorted(r.S2_terminal_broad.unique())})
  embeddings={(correction,space):pd.read_csv(BASE/'pooled'/group/correction/f'{space}.csv',index_col=0)
              for correction in ['NONE','CCA','HARMONY'] for space in ['PCA30','UMAP2']}
  for scope,mask in scopes.items():
   counts=r[mask].groupby('sample').size().reindex(samples,fill_value=0)
   n=min(500,int(counts.min()))
   availability.append(dict(group=group,scope=scope,n_per_sample=n,eligible=n>=20,
      **{f'n_{sample}':int(counts[sample]) for sample in samples}))
   if n<20:continue
   cells=[]
   for sample in samples:
    possible=r.index[mask&r['sample'].eq(sample)].to_numpy()
    chosen=np.sort(rng.choice(possible,size=n,replace=False));cells.extend(chosen)
    selected.extend(dict(group=group,scope=scope,sample=sample,cell_id=c) for c in chosen)
   rr=r.loc[cells];batch=pd.Categorical(rr['sample'],categories=samples).codes
   baseline={}
   for space in ['PCA30','UMAP2']:
    bx=embeddings['NONE',space].loc[cells].to_numpy()
    baseline[space]={s:neighbors(bx[batch==j],15) for j,s in enumerate(samples)}
   for correction in ['NONE','CCA','HARMONY']:
    for space in ['PCA30','UMAP2']:
     x=embeddings[correction,space].loc[cells].to_numpy();assert np.isfinite(x).all()
     nn=neighbors(x,30)
     probabilities=np.stack([(batch[nn]==j).mean(1) for j in range(4)],axis=1)
     entropy=-(probabilities*np.log(np.maximum(probabilities,1e-15))).sum(1)/np.log(4)
     inverse_simpson=1/(probabilities**2).sum(1)
     other=(batch[nn]!=batch[:,None]).mean(1)
     record=dict(group=group,scope=scope,correction=correction,space=space,n_cells=len(cells),
       n_per_sample=n,n_neighbors=nn.shape[1],mean_batch_entropy=float(entropy.mean()),
       mean_batch_inverse_simpson=float(inverse_simpson.mean()),mean_other_sample_neighbor_fraction=float(other.mean()))
     overlaps=[]
     for j,sample in enumerate(samples):
      ix=batch==j;within=neighbors(x[ix],15);orig=baseline[space][sample]
      overlap=np.array([len(set(a)&set(b))/len(a) for a,b in zip(within,orig)])
      overlaps.extend(overlap)
      per_sample.append(dict(group=group,scope=scope,correction=correction,space=space,sample=sample,
          patient=str(rr.loc[rr['sample'].eq(sample),'patient'].iloc[0]),n_cells=int(ix.sum()),
          mean_batch_entropy=float(entropy[ix].mean()),mean_batch_inverse_simpson=float(inverse_simpson[ix].mean()),
          mean_other_sample_neighbor_fraction=float(other[ix].mean()),
          within_sample_neighbor_retention_vs_NONE=float(overlap.mean())))
     record['within_sample_neighbor_retention_vs_NONE']=float(np.mean(overlaps))
     for state in states:
      values=rr[state].to_numpy();smoothed=values[nn].mean(1);variance=float(values.var())
      record[state+'_neighbor_spearman']=float(spearmanr(values,smoothed).statistic) if variance>0 and smoothed.var()>0 else np.nan
      record[state+'_neighbor_normalized_MSE']=float(np.mean((values-smoothed)**2)/variance) if variance>0 else np.nan
     rows.append(record)
     print(group,scope,correction,space,'complete',flush=True)
pd.DataFrame(rows).to_csv(dest/'shared_lineage_metrics.csv',index=False)
pd.DataFrame(per_sample).to_csv(dest/'shared_lineage_metrics_by_sample.csv',index=False)
pd.DataFrame(availability).to_csv(dest/'shared_lineage_availability.csv',index=False)
pd.DataFrame(selected).to_csv(dest/'fixed_balanced_cells.csv.gz',index=False)
write_json(dest/'manifest.json',dict(status='completed',source_sha256=sha(Path(__file__)),
  sample_cap_per_lineage=500,min_required_each_sample=20,seed=42,neighbor_metric='Euclidean',
  scopes='Productive TCR-positive cells and archived S2 broad lineages, separately',
  interpretation='Descriptive batch mixing and biological state preservation; sample,patient and tissue are confounded. RNA modules are not independent validation labels.',
  n_metric_rows=len(rows),completed_at=utc(),job=os.environ['SLURM_JOB_ID']))
(dest/'COMPLETE').write_text(sha(dest/'manifest.json')+'\n')
