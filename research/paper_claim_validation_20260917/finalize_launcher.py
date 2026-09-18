"""Run after the core controller exits; freeze the latest reviewed finalizer sources."""
import json
from common import OUT, write_json, utc
from dispatch import freeze, submit
if __name__=='__main__':
    state=json.loads((OUT/'dispatch_state.json').read_text())
    assert state['n_evaluated_plotted']==726 and not state['review_required']
    source=freeze('core_finalizer')
    jid=submit(source,'core_finalize.py',options=['--job-name=claim_GBM_delivery','--cpus-per-task=4','--mem=24G','--time=04:00:00'])
    write_json(OUT/'core_finalizer_submission.json',dict(job=jid,source=str(source),submitted_at=utc()))
