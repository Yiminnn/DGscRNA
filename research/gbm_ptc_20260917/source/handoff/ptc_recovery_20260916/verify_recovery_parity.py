"""Full-cohort parity of the dbscan index fix and unresolved table-field audit."""
from pathlib import Path
import json
import os
from ptc_common import BASE,RECOVERY,require_slurm,sha,write_json,utc
require_slurm()
import pandas as pd
import numpy as np
from sklearn.metrics import adjusted_rand_score
dest=BASE/'verification/recovery_parity';dest.mkdir(parents=True,exist_ok=True)
rows=[]
for old in ['prior_dbscan125','prior_dbscan126']:
 for correction in ['NONE','CCA','HARMONY']:
  for space in ['PCA30','UMAP2']:
   for method in ['SNN','HDBSCAN_R']:
    rel=Path('pooled/TUT')/correction/f'{space}_{method}'/'clusters.csv'
    a=pd.read_csv(BASE/old/rel,dtype=str).set_index('cell_id').cluster
    b=pd.read_csv(BASE/rel,dtype=str).set_index('cell_id').cluster.loc[a.index]
    equal=int(a.eq(b).sum());ari=float(adjusted_rand_score(a,b))
    rows.append(dict(prior=old,correction=correction,space=space,method=method,n_cells=len(a),
       exact_labels=equal,n_differences=len(a)-equal,ARI=ari))
    assert equal==len(a) and ari==1.0
pd.DataFrame(rows).to_csv(dest/'TUT_all_partitions_exact.csv',index=False)
differences=[]
for path in (RECOVERY/'s3_reconciliation').glob('*.mismatches.csv.gz'):
 d=pd.read_csv(path,index_col=0,keep_default_na=False)
 missing=lambda s:s.astype(str).str.strip().str.lower().isin(['','na','nan','none'])
 both=missing(d.S3)&missing(d.archive)
 differences.append(dict(file=path.name,n_flagged=len(d),both_represent_missing=int(both.sum()),
      literal_equal=int(d.S3.eq(d.archive).sum()),non_missing_disagreements=int((~both&~d.S3.eq(d.archive)).sum())))
pd.DataFrame(differences).to_csv(dest/'S3_disagreement_missingness_diagnosis.csv',index=False)
write_json(dest/'manifest.json',dict(status='passed',n_full_TUT_comparisons=len(rows),
  all_TUT_partitions_exact=True,source_sha256=sha(Path(__file__)),completed_at=utc(),job=os.environ['SLURM_JOB_ID'],
  conclusion='Both official pre-fix successful TUT outputs match the index64 reconstruction cell-for-cell; no clustering parameters changed. Missing-value representation is separated from true source-column disagreement.'))
(dest/'COMPLETE').write_text(sha(dest/'manifest.json')+'\n')
print('Full TUT partition parity passed; S3 missingness audit complete')
