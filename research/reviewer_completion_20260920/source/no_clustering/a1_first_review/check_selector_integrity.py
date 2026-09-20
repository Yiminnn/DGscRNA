"""SLURM metadata checks for final-selector integrity guards; no model fitting."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os

assert os.environ.get('SLURM_JOB_ID')
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
CAMPAIGN=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
OUT=CAMPAIGN/'no_clustering/a1_first_review/selector_integrity'
SOURCE=ROOT/'handoff/reviewer_completion_20260920/embedding/select_and_summarize.py'
spec=importlib.util.spec_from_file_location('a1_selector_checked',SOURCE)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
import pandas as pd
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
frozen=CAMPAIGN/'source_snapshots/embedding_v6'
frozen_manifest=json.loads((frozen/'SOURCE_MANIFEST.json').read_text())
assert sha(frozen/'select_and_summarize.py')==frozen_manifest['select_and_summarize.py']
units=[CAMPAIGN/'embedding/TKU4163/hvg2000/PCA2/evaluation',
       ROOT/'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/GBM/TKU4163/hvg2000/evaluation']
positive=[]
for unit in units:
    module.verify_evaluation(unit)
    frame=pd.read_csv(unit/'metrics.csv')
    module.assert_primary_boolean(frame)
    positive.append(dict(directory=str(unit),manifest_sha256=sha(unit/'manifest.json'),rows=len(frame),primary_dtype=str(frame.primary.dtype)))
rejected=[]
for label,frame in [('strings',pd.DataFrame({'primary':['True','False']})),
                    ('numeric',pd.DataFrame({'primary':[1,0]})),
                    ('nullable_missing',pd.DataFrame({'primary':pd.Series([True,pd.NA],dtype='boolean')}))]:
    try:module.assert_primary_boolean(frame)
    except AssertionError:rejected.append(label)
    else:raise AssertionError('Did not reject '+label)
fixture=OUT/'corruption_fixture';fixture.mkdir(parents=True,exist_ok=True)
(fixture/'metrics.csv').write_text('primary,value\nTrue,1\n')
manifest={'outputs':{'metrics.csv':sha(fixture/'metrics.csv')}}
(fixture/'manifest.json').write_text(json.dumps(manifest))
(fixture/'COMPLETE').write_text(sha(fixture/'manifest.json')+'\n')
module.verify_evaluation(fixture,['metrics.csv'])
(fixture/'metrics.csv').write_text('primary,value\nTrue,2\n')
try:module.verify_evaluation(fixture,['metrics.csv'])
except AssertionError as exc:
    assert 'Changed evaluation output' in str(exc)
    corruption_rejected=True
else:raise AssertionError('Corrupted output was accepted')
result=dict(status='passed',source=str(SOURCE),source_sha256=sha(SOURCE),frozen_v6_unchanged=True,
            real_evaluations=positive,invalid_primary_cases_rejected=rejected,corrupted_metric_rejected=corruption_rejected,
            scope='Integrity and boolean guards only; no final statistics or fitting executed',
            job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),completed_at=datetime.now(timezone.utc).isoformat())
(OUT/'validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2),flush=True)
