"""Publish the PTC archive-content audit as a separately verified supplement."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

from common import ROOT, OUT, checked, sha, write_json, utc
from delivery import STAGE, REMOTE
from ptc_followup_common import PTC, require_ptc


def run():
    require_ptc()
    dest = PTC / 'verification/PTC_followup_archive_content_audit'
    assert checked(dest)
    audit = json.loads((dest / 'manifest.json').read_text())
    submission = PTC / 'verification/PTC_followup_content_audit_submission.json'
    source = Path(json.loads(submission.read_text())['source']) / 'verify_PTC_followup_delivery.py'
    assert sha(source) == audit['source_sha256']
    assert audit['delivery_manifest_sha256'] == sha(OUT / 'PTC_summary/delivery_manifest.json')
    shutil.copy2(source, dest / 'audit_execution_source.py')
    shutil.copy2(__file__, dest / 'delivery_execution_source.py')
    files = [dest / name for name in ['manifest.json', 'COMPLETE',
             'audit_execution_source.py', 'delivery_execution_source.py']] + [submission]
    records = []
    for path in files:
        relative = path.relative_to(ROOT)
        target = STAGE / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        records.append(dict(path=str(relative), sha256=sha(target)))
    listing = dest / 'upload_files.txt'
    listing.write_text('\n'.join(r['path'] for r in records) + '\n')
    subprocess.run(['rclone', 'copy', str(STAGE), REMOTE, '--files-from', str(listing),
                    '--transfers', '2', '--checkers', '4'], check=True)
    subprocess.run(['rclone', 'check', str(STAGE), REMOTE, '--files-from', str(listing),
                    '--download', '--one-way', '--checkers', '2'], check=True)
    receipt = dest / 'AUDIT_DELIVERY_RECEIPT.json'
    write_json(receipt, dict(status='delivered_and_verified', files=records,
        audit_manifest_sha256=sha(dest / 'manifest.json'), remote=REMOTE,
        scope='PTC archive-content audit supplement; the full PTC delivery has its own receipt.',
        verification='Copy and full-download check exited zero; receipt readback verified separately.',
        job=os.environ['SLURM_JOB_ID'], source_sha256=sha(__file__), completed_at=utc()))
    remote_receipt = REMOTE + '/' + str(receipt.relative_to(ROOT))
    subprocess.run(['rclone', 'copyto', str(receipt), remote_receipt], check=True)
    with tempfile.TemporaryDirectory(dir=dest) as folder:
        downloaded = Path(folder) / 'receipt.json'
        subprocess.run(['rclone', 'copyto', remote_receipt, str(downloaded)], check=True)
        assert sha(downloaded) == sha(receipt)
    shutil.copy2(receipt, STAGE / receipt.relative_to(ROOT))
    write_json(dest / 'AUDIT_REMOTE_VERIFIED.json', dict(status='verified',
        receipt_sha256=sha(receipt), job=os.environ['SLURM_JOB_ID'], completed_at=utc()))
    print('PTC_FOLLOWUP_AUDIT_SUPPLEMENT_DELIVERED', flush=True)


if __name__ == '__main__':
    run()
