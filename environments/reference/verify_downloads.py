"""Verify the pristine package archives used for a separate copy installation."""
from pathlib import Path
import hashlib
import json
import os
import sys
from urllib.parse import urlparse

here = Path(__file__).resolve().parent
prefix_root = Path(sys.argv[1]).resolve()
cache = prefix_root/'package_cache'
output = Path(sys.argv[2])
assert not output.exists()
records = []
for runtime, name in [('r', 'r-linux-64.explicit.txt'), ('python', 'python-linux-64.explicit.txt')]:
    provenance = json.loads((here/(name+'.provenance.json')).read_text())
    metadata = {v['name']: v for p in (prefix_root/runtime/'conda-meta').glob('*.json')
                for v in [json.loads(p.read_text())]}
    assert set(metadata) == {v['name'] for v in provenance['packages']}
    for package in provenance['packages']:
        installed = metadata[package['name']]
        assert installed['version'] == package['version'] and installed['build'] == package['build']
        archive = cache/Path(urlparse(package['url']).path).name
        assert archive.is_file(), archive
        with archive.open('rb') as stream:
            md5 = hashlib.file_digest(stream, 'md5').hexdigest()
        with archive.open('rb') as stream:
            sha256 = hashlib.file_digest(stream, 'sha256').hexdigest()
        assert md5 == package['md5'], archive
        if package.get('sha256'):
            assert sha256 == package['sha256'], archive
        records.append(dict(runtime=runtime, name=package['name'], version=package['version'],
                            archive=str(archive), md5=md5, sha256=sha256,
                            installed_prefix=str(prefix_root/runtime), source_url=package['url']))
proof = dict(status='pristine_download_checksums_and_installed_builds_verified',
             package_cache=str(cache), records=records, package_records=len(records),
             separate_new_prefixes=True, copy_install=True,
             job=os.environ.get('SLURM_JOB_ID', 'local'), step=os.environ.get('SLURM_STEP_ID'))
output.write_text(json.dumps(proof, indent=2)+'\n')
print('PRISTINE_PACKAGE_PROOF_COMPLETE', len(records), flush=True)
