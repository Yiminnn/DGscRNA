"""SLURM-only numerical regression against successful frozen v6 ICA geometry."""
from pathlib import Path
from datetime import datetime,timezone
import json
import os
import shutil
import sys
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parent))
import run as core
import cache_compatibility
core.require_slurm()
cache_compatibility.own_source(core)
sample,budget='TKU4163','hvg2000'
reference=core.CAMPAIGN/'embedding_geometry_smoke'/sample/budget/'ICA2'
assert core.checked(reference,'representation.json','REPRESENTATION_COMPLETE')
old=json.loads((reference/'representation.json').read_text())
assert old['params']['algorithm']=='parallel' and old['params']['max_iter']==5000
assert not any(w.startswith('ConvergenceWarning:') for w in old['warnings'])
base=core.CAMPAIGN/'embedding_v8_regression'/sample/budget
dest=base/'fresh_ICA2'
assert not dest.exists(), 'Preserve prior regression attempts'
dest.mkdir(parents=True)
geometry=core.OUTPUT/sample/budget/'geometry'
cfg=dict(sample=sample,budget=budget,space='ICA2',prep=str(core.REFERENCE/'GBM'/sample/budget),
    dest=str(dest),geometry=str(geometry),embedding=str(dest/'embedding.csv'),
    hdbscan_dest=str(dest/'HDBSCAN_R'),frozen_protocol_sha256=core.sha(core.CAMPAIGN/'protocol/embedding.json'),
    source_sha256=core.sha(core.CODE/'run.py'),reference_labels_used_for_fit=False,
    purpose='Successful-parallel numerical regression only; no DEG, DL or truth reads')
core.write_json(dest/'config.json',cfg)
x,cells,gm=core.get_geometry(cfg)
z=core.make_embedding(cfg,x,cells,gm)
rep=json.loads((dest/'representation.json').read_text())
assert len(rep['attempts'])==1 and rep['actual_solver']=='parallel' and rep['actual_iteration_cap']==5000
assert rep['canonical_parallel5000_status']=='converged' and rep['fallback_used'] is False
assert rep['n_iter_']==old['n_iter_']
assert core.sha(dest/'embedding.csv')==core.sha(reference/'embedding.csv')
conditions=core.make_partitions(cfg,z,cells)
assert len(conditions)==13
records={}
for condition in conditions:
    directory=Path(condition['dest']);name=directory.name
    assert core.checked(reference/name,'partition_manifest.json','PARTITION_COMPLETE')
    a=pd.read_csv(directory/'clusters.csv',dtype=str,keep_default_na=False)
    b=pd.read_csv(reference/name/'clusters.csv',dtype=str,keep_default_na=False)
    assert a.equals(b),name
    records[name]=dict(clusters_sha256=core.sha(directory/'clusters.csv'),
        reference_clusters_sha256=core.sha(reference/name/'clusters.csv'),
        partition_manifest_sha256=core.sha(directory/'partition_manifest.json'),
        cell_order_and_labels_exact=True)
# Consume an isolated byte-exact copy of the old representation. Only the new
# test tree receives consumer receipts; original gate/production caches do not.
reuse=base/'reuse_v6_ICA2';reuse.mkdir()
for name in ('representation.json','REPRESENTATION_COMPLETE','embedding.csv'):
    shutil.copy2(reference/name,reuse/name)
reuse_cfg={**cfg,'dest':str(reuse),'embedding':str(reuse/'embedding.csv')}
before=core.sha(reuse/'representation.json')
reused=core.make_embedding(reuse_cfg,x,cells,gm)
expected=pd.read_csv(reference/'embedding.csv')[['x','y']].to_numpy()
np.testing.assert_array_equal(reused,expected)
assert core.sha(reuse/'representation.json')==before==core.sha(reference/'representation.json')
receipts=list((reuse/'cache_consumers_v8').glob('*.json'))
assert len(receipts)==1
consumer=json.loads(receipts[0].read_text())
assert consumer['producer_source_sha256']==old['source_sha256']
assert consumer['consumer_source_sha256']==core.sha(core.CODE/'run.py')
assert core.sha(consumer['producer_manifest'])==before
try:
    cache_compatibility.accepted_producer(core,'0'*64)
    raise AssertionError('Unapproved producer accepted')
except AssertionError as exc:
    assert 'Unapproved cache producer' in str(exc)
report=dict(status='passed',scope='Successful ICA branch and source-pinned cache reuse; no production activation',
    sample=sample,budget=budget,n_cells=len(cells),n_features=x.shape[1],
    old_source_sha256=old['source_sha256'],new_source_sha256=core.sha(core.CODE/'run.py'),
    source_manifest_sha256=core.sha(core.CODE/'SOURCE_MANIFEST.json'),
    protocol_sha256=cfg['frozen_protocol_sha256'],
    adaptive_policy_sha256=core.sha(core.CAMPAIGN/'protocol/embedding_convergence_repair_20260920_v2.json'),
    original_parallel5000_constructor_used=True,n_iter=rep['n_iter_'],
    embedding_csv_bytes_exact=True,all_13_partition_labels_and_cells_exact=True,
    partitions=records,old_cache_representation_bytes_unchanged=True,
    producer_and_consumer_sources_separate=True,unapproved_producer_rejected=True,
    consumer_receipt_sha256=core.sha(receipts[0]),author_labels_read=False,models_DL_fitted=False,
    job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),
    completed_at=datetime.now(timezone.utc).isoformat())
core.write_json(base/'validation.json',report)
print(json.dumps(report,indent=2),flush=True)
