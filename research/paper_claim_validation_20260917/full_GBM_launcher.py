"""Observe independent completion gates, then freeze and submit the GBM finalizer."""
import csv
import json
import time
from common import OUT,checked,write_json,utc
from dispatch import freeze,submit

def ready():
    checks={'core_delivered':(OUT/'summary/GBM_CORE_DELIVERED.json').exists()}
    for label,file,expected in [('auxiliary','aux_dispatch_state.json',780),('representation','representation_dispatch_state.json',180),('scalability','scalability_dispatch_state.json',15)]:
        path=OUT/file
        state=json.loads(path.read_text()) if path.exists() else {}
        checks[label]=state.get('completed')==expected and not state.get('review_required')
    samples=list(csv.DictReader((OUT/'protocol/cohort.csv').open()))
    for method in ['scType','scCATCH','SCINA','SingleR','scDeepSort']:
        checks[method]=all(checked(OUT/'comparators'/method/r['sample']/'evaluation') for r in samples)
    write_json(OUT/'full_GBM_completion_gates.json',dict(checks=checks,ready=all(checks.values()),checked_at=utc()))
    return all(checks.values())

if __name__=='__main__':
    while not ready():time.sleep(45)
    path=OUT/'full_GBM_finalizer_submission.json'
    if not path.exists():
        source=freeze('full_GBM_finalize')
        job=submit(source,'full_GBM_finalize.py',[],['--job-name=claim_GBM_full_finalize','--cpus-per-task=4','--mem=32G','--time=06:00:00'])
        write_json(path,dict(job=job,source=str(source),submitted_at=utc()))
