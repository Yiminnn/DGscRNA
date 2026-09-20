"""One-CPU SLURM watcher for independently validating the two large A1 pilots."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import subprocess
import sys
import time

ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
CAMP=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
OUT=CAMP/'no_clustering/a1_first_review'
WATCH=OUT/'watcher'
SOURCE=OUT/'verifier_v2'
TASKS=['NL022','SN040']


def sha(path):
    with path.open('rb') as handle:return hashlib.file_digest(handle,'sha256').hexdigest()


def write(path,value):
    tmp=path.with_name(path.name+f'.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)


def checked(directory,manifest='manifest.json',flag='COMPLETE'):
    try:return (directory/flag).read_text().strip()==sha(directory/manifest)
    except OSError:return False


def matches(sample):
    base=CAMP/'embedding'/sample/'hvg2000/PCA2'
    path=OUT/sample/'hvg2000/PCA2/validation.json'
    if not path.exists():return False
    proof=json.loads(path.read_text())
    return (proof['status']=='passed' and proof['validation_source_sha256']==sha(SOURCE/'validate_unit.py')
            and proof['fit_manifest_sha256']==sha(base/'fit_manifest.json')
            and proof['evaluation_manifest_sha256']==sha(base/'evaluation/manifest.json')
            and proof['figures_manifest_sha256']==sha(base/'figures/manifest.json'))


def main():
    assert os.environ.get('SLURM_JOB_ID')
    WATCH.mkdir(exist_ok=True)
    lock=WATCH/'ACTIVE'
    lock.mkdir()  # Fail closed if a second watcher is started.
    write(lock/'owner.json',dict(job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),pid=os.getpid()))
    deadline=time.monotonic()+6*24*3600
    try:
        manifest=json.loads((SOURCE/'SOURCE_MANIFEST.json').read_text())
        while time.monotonic()<deadline and not (WATCH/'STOP').exists():
            for name,digest in manifest.items():assert sha(SOURCE/name)==digest
            for sample in TASKS:
                if matches(sample):continue
                unit=CAMP/'embedding'/sample/'hvg2000/PCA2'
                ready=(checked(unit,'fit_manifest.json','FIT_COMPLETE') and checked(unit/'evaluation') and checked(unit/'figures'))
                if not ready:continue
                print('INDEPENDENT_A1_VALIDATION_START',sample,flush=True)
                with (WATCH/(sample+'.log')).open('a') as log:
                    subprocess.run([sys.executable,str(SOURCE/'validate_unit.py'),'--sample',sample,'--budget','hvg2000','--space','PCA2'],
                                   stdout=log,stderr=subprocess.STDOUT,check=True,timeout=3600)
                assert matches(sample)
                print('INDEPENDENT_A1_VALIDATION_PASSED',sample,flush=True)
            complete=[sample for sample in TASKS if matches(sample)]
            state=dict(stage='A1_PILOT_INDEPENDENT',work_package='A',status='passed' if len(complete)==2 else 'waiting',
                       completed=complete,remaining=[sample for sample in TASKS if sample not in complete],
                       numerical_validation_only=True,whole_work_package_A_complete=False,
                       job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),
                       updated_at=datetime.now(timezone.utc).isoformat())
            write(WATCH/'status.json',state)
            if len(complete)==2:
                write(WATCH/'manifest.json',dict(**state,proofs={sample:sha(OUT/sample/'hvg2000/PCA2/validation.json') for sample in TASKS}))
                (WATCH/'COMPLETE').write_text(sha(WATCH/'manifest.json')+'\n')
                return
            time.sleep(45)
        raise RuntimeError('Independent pilot watcher stopped before both validation proofs completed')
    except Exception as exc:
        write(WATCH/'status.json',dict(stage='A1_PILOT_INDEPENDENT',work_package='A',status='failed',error=repr(exc),
              whole_work_package_A_complete=False,job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),
              updated_at=datetime.now(timezone.utc).isoformat()))
        raise
    finally:
        (lock/'owner.json').unlink();lock.rmdir()


if __name__=='__main__':main()
