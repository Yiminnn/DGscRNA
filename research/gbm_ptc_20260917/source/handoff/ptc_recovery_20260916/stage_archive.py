#!/usr/bin/env python3
"""Stage exact selected PTC archive members after verified GBM delivery only."""
import hashlib
import json
import os
from pathlib import Path,PurePosixPath
import sys
import tarfile
import time
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'hvg_ptc_20260916'))
from common import ROOT,OUT,require_slurm,sha,utc,write_json,runtime_record


def run():
    require_slurm()
    gate=OUT/'GBM_DELIVERY.json'
    assert gate.exists(),'GBM delivery must finish before PTC archive extraction'
    delivery=json.loads(gate.read_text())
    assert delivery['status']=='completed'
    assert delivery['report_sha256']==sha(OUT/'GBM_REPORT.md')
    assert delivery['notebook_manifest_sha256']==sha(OUT/'notebooks/execution_manifest.json')
    archive=Path('/fs/scratch/PCON0080/yimin/_tcr_stage/tcr.tar.gz')
    roster=Path(__file__).with_name('archive_members.txt')
    wanted=set(roster.read_text().splitlines());assert wanted
    for name in wanted:
        p=PurePosixPath(name)
        assert not p.is_absolute() and '..' not in p.parts and p.parts[0]=='tcr'
    target=OUT/'ptc_recovery';target.mkdir(exist_ok=True)
    dest=target/'archive';dest.mkdir(exist_ok=True)
    start=time.monotonic();before=archive.stat();found={}
    prior_log=target/'archive_extraction.jsonl'
    if prior_log.exists():
        for line in prior_log.read_text().splitlines():
            row=json.loads(line)
            path=dest/row['member']
            if path.exists() and path.stat().st_size==row['size'] and sha(path)==row['sha256']:
                found[row['member']]=row
    remaining=wanted-set(found)
    print('Requested archive members',len(wanted),'already verified',len(found),flush=True)
    if remaining:
        last_report=time.monotonic()
        with archive.open('rb') as compressed, tarfile.open(fileobj=compressed,mode='r|gz',bufsize=1024*1024) as tar:
            for member in tar:
                if time.monotonic()-last_report>=60:
                    print('Archive scan:',compressed.tell(),'/',before.st_size,'compressed bytes;',
                          len(found),'/',len(wanted),'selected members staged; current',member.name,flush=True)
                    last_report=time.monotonic()
                if member.name not in remaining:continue
                assert member.isfile(),f'Unexpected non-regular archive member: {member.name}'
                path=dest/member.name;path.parent.mkdir(parents=True,exist_ok=True)
                tmp=path.with_name(path.name+'.part')
                h=hashlib.sha256();n=0
                with tar.extractfile(member) as inp,tmp.open('wb') as output:
                    while True:
                        block=inp.read(8*1024*1024)
                        if not block:break
                        output.write(block);h.update(block);n+=len(block)
                assert n==member.size
                tmp.replace(path)
                row=dict(member=member.name,size=n,sha256=h.hexdigest(),archive_mtime=member.mtime,
                         extracted_at=utc(),slurm_job_id=os.environ['SLURM_JOB_ID'])
                with prior_log.open('a') as log:log.write(json.dumps(row)+'\n')
                found[member.name]=row;remaining.remove(member.name)
                print(len(found),'/',len(wanted),member.name,n,'bytes',flush=True)
                if not remaining:break
    after=archive.stat()
    assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns),'Source archive changed'
    assert not remaining,sorted(remaining)
    write_json(target/'archive_manifest.json',dict(status='completed',timestamp=utc(),
        gbm_delivery_sha256=sha(gate),source_archive=str(archive),source_archive_size=before.st_size,
        source_archive_mtime_ns=before.st_mtime_ns,selection_sha256=sha(roster),
        extractor_sha256=sha(Path(__file__)),n_members=len(found),members=found,
        seconds=time.monotonic()-start,**runtime_record()))
    (target/'ARCHIVE_STAGED').write_text(sha(target/'archive_manifest.json')+'\n')
    print('PTC archive staged with per-file provenance',flush=True)


if __name__=='__main__':run()
