"""Recover the previously unselected original writing note for metric provenance."""
from pathlib import Path
import hashlib
import json
import os
import tarfile

assert os.environ.get('SLURM_JOB_ID')
ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT = ROOT / 'results/hvg_ptc_20260916_v1/ptc_paper_baseline'
source = Path('/fs/scratch/PCON0080/yimin/_tcr_stage/tcr.tar.gz')
member_name = 'tcr/scripts/paper.md'
dest = OUT / 'original_writing_note'
dest.mkdir(exist_ok=True)
with source.open('rb') as handle, tarfile.open(fileobj=handle, mode='r|gz') as archive:
    for member in archive:
        if member.name != member_name:
            continue
        assert member.isfile() and member.size < 20 * 1024 * 1024
        content = archive.extractfile(member).read()
        path = dest / 'paper.md'
        if path.exists():
            assert path.read_bytes() == content
        else:
            path.write_bytes(content)
        record = dict(job=os.environ['SLURM_JOB_ID'], source_archive=str(source),
            member=member_name, archive_mtime=member.mtime, bytes=len(content),
            sha256=hashlib.sha256(content).hexdigest(), compressed_bytes_read=handle.tell())
        (dest / 'manifest.json').write_text(json.dumps(record, indent=2) + '\n')
        print(json.dumps(record, indent=2), flush=True)
        break
    else:
        raise RuntimeError('Expected source note absent')
