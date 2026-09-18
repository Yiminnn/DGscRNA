"""Copy derived results into the existing authorized directory; verify remote bytes."""
import ast
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
from common import ROOT, CODE, OUT, FEATURES, ROUTES, require_slurm, checked, sha, write_json, utc

STAGE=ROOT/'results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916/GBM_PTC_results_20260916'
REMOTE='onedrive:work_od/share/dgscrna_GSE274546_TKU3186/01_report/GBM_PTC_results_20260916'

def run():
    require_slurm()
    assert checked(OUT/'summary','aggregate_manifest.json','AGGREGATE_COMPLETE')
    nb=json.loads((OUT/'summary/notebook_update_manifest.json').read_text())
    assert nb['current_sha256']==sha(ROOT/'notebooks/dgscrna_results.ipynb')
    files=[]
    def copy(p):
        target=STAGE/p.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(p,target);files.append(str(target.relative_to(STAGE)))
    def bundle(paths,target,relative):
        target.parent.mkdir(parents=True,exist_ok=True);tmp=target.with_name(target.name+'.part')
        with tarfile.open(tmp,'w:gz',compresslevel=1) as arc:
            for p in sorted(paths):arc.add(p,arcname=str(p.relative_to(relative)),recursive=False)
        tmp.replace(target);files.append(str(target.relative_to(STAGE)))
    copy(ROOT/'notebooks/dgscrna_results.ipynb')
    source=Path(__file__).resolve().parent
    codefiles=[p for p in source.iterdir() if p.suffix in ['.py','.R','.sbatch']]
    for p in codefiles:
        if p.suffix=='.py':ast.parse(p.read_text(),filename=str(p))
        if p.suffix=='.sbatch':subprocess.run(['bash','-n',str(p)],check=True)
    # Delivered executable sources are the actual immutable finalizer bundle.
    bundle(codefiles,STAGE/OUT.relative_to(ROOT)/'source/finalizer_source.tar.gz',source)
    for p in CODE.glob('*.md'):copy(p)
    sourcebundles=[p for p in (OUT/'source_snapshots').rglob('*') if p.is_file() and p.suffix in ['.py','.R','.sbatch','.json']]
    bundle(sourcebundles,STAGE/OUT.relative_to(ROOT)/'source/executed_source_snapshots.tar.gz',OUT/'source_snapshots')
    # Preserve third-party source/license attribution; do not package environments.
    bundle([p for p in (CODE/'vendor_sources').iterdir() if p.is_file()],
        STAGE/OUT.relative_to(ROOT)/'source/comparator_sources.tar.gz',CODE/'vendor_sources')
    for folder in ['summary','protocol','markers','verification']:
        for p in (OUT/folder).rglob('*'):
            if not p.is_file() or p.suffix in ['.rds','.bin','.npz','.lock']:continue
            if p.name in ['DELIVERY_RECEIPT.json','REMOTE_RECEIPT_UPLOADED.json','GBM_CORE_DELIVERED.json','core_delivery_manifest.json']:continue
            copy(p)
    cohort=list(csv.DictReader((OUT/'protocol/cohort.csv').open()))
    artifacts=[]
    for r in cohort:
        sample=r['sample'];base=OUT/'GBM'/sample;paths=[]
        for budget in FEATURES:
            unit=base/budget
            assert checked(unit/'evaluation') and checked(unit/'figures','manifest.json','FIGURES_COMPLETE')
            for p in (unit/'figures').glob('*'):
                if p.is_file():copy(p)
            for p in (unit/'evaluation').glob('*'):
                if p.is_file():paths.append(p)
            for pattern in ['*manifest.json','*.txt','cells.csv','PCA30.csv','UMAP2.csv','*features.csv','*genes.csv']:
                paths.extend(p for p in unit.glob(pattern) if p.is_file())
            for route in ROUTES:
                src=unit/route
                for name in ['score_manifest.json','clusters.csv','cells.csv','initial_calls.csv.gz','SCORE_COMPLETE']:
                    if (src/name).exists():paths.append(src/name)
                for p in (src/'terminal').glob('*/*'):
                    if p.is_file() and p.name in ['predictions.csv.gz','terminal_manifest.json','training_manifest.json','training_history.json','history.csv','TERMINAL_COMPLETE']:
                        paths.append(p)
        bundle(set(paths),STAGE/OUT.relative_to(ROOT)/'artifacts'/f'{sample}_all_core_labels_metrics.tar.gz',base)
        artifacts.append(dict(sample=sample,units=6,clustering_results=24,terminal_conditions=1152))
    # Scheduler accounting includes failed attempts, not just successful fits.
    state=json.loads((OUT/'dispatch_state.json').read_text());jobs={state['main_array']}
    jobs.update(str(j) for j in [7362290,7362298,7362321,7362326,7362345,7362359,7362361,7362447,7362470,7362487,7362488,7362561])
    for r in state['units'].values():
        for field in ['post_job','retry_job']:
            if r.get(field):jobs.add(r[field])
    resource=OUT/'resources';resource.mkdir(exist_ok=True)
    proc=subprocess.run(['sacct','-j',','.join(sorted(jobs)),'-P','--format=JobID,JobName,State,ExitCode,ElapsedRaw,AllocCPUS,ReqMem,MaxRSS,NodeList,Start,End'],text=True,capture_output=True,check=True)
    (resource/'core_slurm_accounting.psv').write_text(proc.stdout);copy(resource/'core_slurm_accounting.psv')
    write_json(resource/'manifest.json',dict(status='collected',job_ids=sorted(jobs),
        statement='Actual allocation/job/step records; core per-sample jobs do not prove single-run >100k scalability.',
        job=os.environ['SLURM_JOB_ID'],completed_at=utc()));copy(resource/'manifest.json')
    manifest=OUT/'summary/core_delivery_manifest.json'
    records=[dict(path=s,sha256=sha(STAGE/s),bytes=(STAGE/s).stat().st_size) for s in sorted(set(files))]
    write_json(manifest,dict(status='staged',scope='GBM core only',files=records,n_files=len(records),
        total_bytes=sum(v['bytes'] for v in records),sample_artifacts=artifacts,remote=REMOTE,
        old_results_preserved=True,job=os.environ['SLURM_JOB_ID'],created_at=utc()))
    copy(manifest)
    listing=OUT/'core_upload_files.txt';listing.write_text('\n'.join(sorted(set(files)))+'\n')
    subprocess.run(['rclone','copy',str(STAGE),REMOTE,'--files-from',str(listing),'--transfers','4','--checkers','8',
        '--log-file',str(OUT/'core_upload.log'),'--log-level','INFO'],check=True)
    subprocess.run(['rclone','check',str(STAGE),REMOTE,'--files-from',str(listing),'--download','--one-way','--checkers','4',
        '--log-file',str(OUT/'core_remote_check.log'),'--log-level','INFO'],check=True)
    receipt=OUT/'summary/DELIVERY_RECEIPT.json'
    write_json(receipt,dict(status='delivered_and_verified',scope='GBM core only',n_files=len(set(files)),
        manifest_sha256=sha(manifest),notebook_sha256=sha(ROOT/'notebooks/dgscrna_results.ipynb'),
        remote=REMOTE,verification='rclone copy then full-download one-way check both exited zero',
        job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    remote_receipt=REMOTE+'/'+str(receipt.relative_to(ROOT))
    subprocess.run(['rclone','copyto',str(receipt),remote_receipt],check=True)
    with tempfile.TemporaryDirectory(dir=OUT,prefix='receipt_readback_') as td:
        got=Path(td)/'receipt.json';subprocess.run(['rclone','copyto',remote_receipt,str(got)],check=True)
        assert sha(got)==sha(receipt)
    copy(receipt)
    write_json(OUT/'summary/REMOTE_RECEIPT_UPLOADED.json',dict(status='remote_receipt_uploaded_and_verified',
        receipt_sha256=sha(receipt),job=os.environ['SLURM_JOB_ID'],completed_at=utc()))

if __name__=='__main__':run()
