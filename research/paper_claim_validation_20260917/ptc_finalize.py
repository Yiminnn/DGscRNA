"""Complete PTC analysis, canonical notebook and verified existing-directory delivery."""
import json
import os
from pathlib import Path
from common import OUT,sha,checked,write_json,utc
from ptc_followup_common import PTC,require_ptc

def run():
    require_ptc()
    import ptc_aggregate,ptc_comparator_replay,ptc_report,ptc_delivery
    if not checked(OUT/'PTC_summary'):ptc_aggregate.run()
    if not checked(OUT/'PTC_comparator_replay'):ptc_comparator_replay.run()
    if not (OUT/'PTC_summary/notebook_manifest.json').exists():ptc_report.run()
    ptc_delivery.run()
    # Audit the actual first package in this SLURM allocation before publishing
    # the overall completion gate; no scientific fits or metrics are rerun.
    import verify_PTC_followup_delivery,deliver_PTC_followup_audit
    audit=PTC/'verification/PTC_followup_archive_content_audit'
    if not checked(audit):
        write_json(PTC/'verification/PTC_followup_content_audit_submission.json',dict(
            job=os.environ['SLURM_JOB_ID'],source=str(Path(__file__).resolve().parent),
            mode='Audit phase inside the existing finalizer allocation; no separate job submitted',
            main_delivery_returned_successfully=True,started_at=utc()))
        verify_PTC_followup_delivery.run()
    if not (audit/'AUDIT_REMOTE_VERIFIED.json').exists():deliver_PTC_followup_audit.run()
    audit_remote=json.loads((audit/'AUDIT_REMOTE_VERIFIED.json').read_text())
    assert audit_remote['status']=='verified'
    assert audit_remote['receipt_sha256']==sha(audit/'AUDIT_DELIVERY_RECEIPT.json')
    receipt=json.loads((OUT/'PTC_summary/DELIVERY_RECEIPT.json').read_text())
    write_json(OUT/'PTC_summary/PTC_FULL_DELIVERED.json',dict(status='accepted_GBM_PTC_followups_completed',
        receipt_sha256=sha(OUT/'PTC_summary/DELIVERY_RECEIPT.json'),
        content_audit_manifest_sha256=sha(audit/'manifest.json'),
        content_audit_delivery_receipt_sha256=sha(audit/'AUDIT_DELIVERY_RECEIPT.json'),
        notebook_sha256=receipt['notebook_sha256'],
        job=os.environ['SLURM_JOB_ID'],completed_at=utc()))

if __name__=='__main__':run()
