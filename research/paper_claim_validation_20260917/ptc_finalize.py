"""Complete PTC analysis, canonical notebook and verified existing-directory delivery."""
import json
import os
from common import OUT,sha,checked,write_json,utc
from ptc_followup_common import require_ptc

def run():
    require_ptc()
    import ptc_aggregate,ptc_comparator_replay,ptc_report,ptc_delivery
    if not checked(OUT/'PTC_summary'):ptc_aggregate.run()
    if not checked(OUT/'PTC_comparator_replay'):ptc_comparator_replay.run()
    if not (OUT/'PTC_summary/notebook_manifest.json').exists():ptc_report.run()
    ptc_delivery.run()
    receipt=json.loads((OUT/'PTC_summary/DELIVERY_RECEIPT.json').read_text())
    write_json(OUT/'PTC_summary/PTC_FULL_DELIVERED.json',dict(status='accepted_GBM_PTC_followups_completed',
        receipt_sha256=sha(OUT/'PTC_summary/DELIVERY_RECEIPT.json'),
        notebook_sha256=receipt['notebook_sha256'],
        job=os.environ['SLURM_JOB_ID'],completed_at=utc()))

if __name__=='__main__':run()
