"""Complete all GBM follow-ups, original notebook and verified delivery in order."""
import json
import os
from common import OUT,require_slurm,checked,sha,write_json,utc

def run():
    require_slurm()
    import aggregate_controls,aggregate_comparators,comparison_diagnostics,unknown_diagnostics
    import marker_evidence_audit,aggregate_scalability,full_GBM_report,deliver_GBM_followups,legacy_coverage_audit,unknown_expression,workflow_choice_summary,verify_core_delivery
    assert (OUT/'summary/GBM_CORE_DELIVERED.json').exists()
    audit=OUT/'verification/core_archive_content_audit'
    if not checked(audit):verify_core_delivery.run()
    audited=json.loads((audit/'manifest.json').read_text())
    assert audited['core_delivery_manifest_sha256']==sha(OUT/'summary/core_delivery_manifest.json')
    assert audited['receipt_sha256']==sha(OUT/'summary/DELIVERY_RECEIPT.json')
    for folder,module in [('controls_summary',aggregate_controls),('workflow_choice_summary',workflow_choice_summary),('comparison_summary',aggregate_comparators),
                          ('comparison_summary/diagnostics',comparison_diagnostics),('unknown_summary',unknown_diagnostics),
                          ('unknown_expression',unknown_expression),
                          ('marker_evidence_summary',marker_evidence_audit),('scalability_summary',aggregate_scalability),
                          ('legacy_coverage_audit',legacy_coverage_audit)]:
        if not checked(OUT/folder):module.run()
    # Refresh source-evidence wording without changing the frozen marker roster.
    marker_evidence_audit.run()
    if not checked(OUT/'GBM_full_summary'):full_GBM_report.run()
    deliver_GBM_followups.run()
    write_json(OUT/'GBM_full_summary/GBM_FULL_DELIVERED.json',dict(status='GBM_complete_PTC_followups_pending',
        receipt_sha256=sha(OUT/'GBM_full_summary/DELIVERY_RECEIPT.json'),PTC_compute_may_start=True,
        job=os.environ['SLURM_JOB_ID'],completed_at=utc()))

if __name__=='__main__':run()
