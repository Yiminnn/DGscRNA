"""Record verified payload delivery; separately verify the uploaded receipt itself."""
import datetime
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from stage_delivery import ROOT, OUT, PARENT, STAGE, REMOTE, sha


def atomic_json(path, value):
    temporary = path.with_name(path.name + '.part')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def record_payload_delivery():
    plan = PARENT / 'R_reference_campaign_20260917_delivery_manifest.json'
    manifest = json.loads(plan.read_text())
    assert manifest['job'] == os.environ['SLURM_JOB_ID']
    record = dict(status='delivered_and_verified',
        time=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        job=os.environ['SLURM_JOB_ID'], remote=REMOTE, n_files=manifest.get('n_upload_files',manifest['n_files']),
        payload_files=manifest['n_files'],receipt_itself_excluded_from_file_count=True,
        total_bytes=manifest['total_bytes']+(plan.stat().st_size if manifest.get('delivered_manifest') else 0),
        manifest_sha256=sha(plan),delivered_manifest=manifest.get('delivered_manifest'),
        notebook_sha256=sha(ROOT / 'notebooks/dgscrna_results.ipynb'),
        verification='rclone copy and one-way full-download check both exited 0; no differences',
        receipt_upload_verification='Separate REMOTE_RECEIPT_UPLOADED.json marker is required.',
        old_results_preserved=True, website_modified=False, new_version_notebook_created=False)
    path = OUT / 'summary/DELIVERY_RECEIPT.json'
    atomic_json(path, record)
    target = STAGE / path.relative_to(ROOT)
    shutil.copy2(path, target)
    print(json.dumps(record, indent=2), flush=True)


def upload_and_verify_receipt(receipt, marker, remote_path, runner=subprocess.run):
    """Never publish the success marker until an exact remote readback matches."""
    receipt = Path(receipt)
    marker = Path(marker)
    record = json.loads(receipt.read_text())
    assert record['job'] == os.environ['SLURM_JOB_ID']
    receipt_hash = sha(receipt)
    # An earlier successful marker must not describe a newly attempted finalizer.
    if marker.exists():
        previous_hash = sha(marker)
        preserved = marker.with_name(marker.stem + '.previous_' + previous_hash[:12] + '.json')
        if preserved.exists():
            assert sha(preserved) == previous_hash
            marker.unlink()
        else:
            marker.rename(preserved)
    runner(['rclone', 'copyto', str(receipt), remote_path], check=True)
    with tempfile.TemporaryDirectory(prefix='receipt_readback_', dir=marker.parent) as directory:
        downloaded = Path(directory) / 'DELIVERY_RECEIPT.json'
        runner(['rclone', 'copyto', remote_path, str(downloaded)], check=True)
        assert sha(downloaded) == receipt_hash, 'Remote receipt hash differs after upload'
    assert sha(receipt) == receipt_hash, 'Local receipt changed during verification'
    success = dict(status='remote_receipt_uploaded_and_verified',
        time=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        job=os.environ['SLURM_JOB_ID'], remote_receipt=remote_path,
        receipt_sha256=receipt_hash, manifest_sha256=record['manifest_sha256'],
        verification='Receipt upload and independent full-download SHA256 comparison succeeded.',
        finalizer_exit_requires_independent_SLURM_check=True)
    atomic_json(marker, success)
    return success


def self_test():
    """Exercise copy failure, corrupted readback and successful remote receipt gates."""
    with tempfile.TemporaryDirectory(prefix='receipt_verification_test_') as directory:
        root = Path(directory)
        receipt = root / 'DELIVERY_RECEIPT.json'
        marker = root / 'REMOTE_RECEIPT_UPLOADED.json'
        remote = root / 'remote_receipt.json'
        atomic_json(receipt, {'job': os.environ['SLURM_JOB_ID'], 'manifest_sha256': 'fixture'})
        def fail_copy(command, **kwargs):
            raise subprocess.CalledProcessError(1, command)
        try:
            upload_and_verify_receipt(receipt, marker, str(remote), fail_copy)
        except subprocess.CalledProcessError:
            pass
        else:
            raise AssertionError('Upload failure was not propagated')
        assert not marker.exists()
        def corrupt_download(command, **kwargs):
            if command[2] == str(receipt):
                shutil.copy2(command[2], command[3])
            else:
                Path(command[3]).write_text('corrupted remote content')
        try:
            upload_and_verify_receipt(receipt, marker, str(remote), corrupt_download)
        except AssertionError as exc:
            assert 'Remote receipt hash differs' in str(exc)
        else:
            raise AssertionError('Corrupted remote receipt was accepted')
        assert not marker.exists()
        def local_copy(command, **kwargs):
            shutil.copy2(command[2], command[3])
        upload_and_verify_receipt(receipt, marker, str(remote), local_copy)
        assert json.loads(marker.read_text())['receipt_sha256'] == sha(receipt)
        # Retrying a previously successful upload must invalidate its active marker.
        try:
            upload_and_verify_receipt(receipt, marker, str(remote), fail_copy)
        except subprocess.CalledProcessError:
            pass
        assert not marker.exists() and len(list(root.glob('REMOTE_RECEIPT_UPLOADED.previous_*.json'))) == 1
    print('Receipt verification self-test passed: upload failure, mismatched readback, success, stale-marker retry.', flush=True)


if __name__ == '__main__':
    assert os.environ.get('SLURM_JOB_ID')
    if '--self-test' in sys.argv:
        self_test()
    elif '--upload-receipt' in sys.argv:
        receipt = OUT / 'summary/DELIVERY_RECEIPT.json'
        plan = PARENT / 'R_reference_campaign_20260917_delivery_manifest.json'
        assert json.loads(receipt.read_text())['manifest_sha256'] == sha(plan)
        success = upload_and_verify_receipt(receipt, OUT / 'summary/REMOTE_RECEIPT_UPLOADED.json',
            REMOTE + '/' + str(receipt.relative_to(ROOT)))
        print(json.dumps(success, indent=2), flush=True)
    else:
        record_payload_delivery()
