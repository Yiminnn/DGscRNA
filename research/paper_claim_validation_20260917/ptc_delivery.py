"""Deliver PTC additions to the same authorized OneDrive directory and verify them."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
from common import ROOT,CODE,OUT,sha,checked,write_json,utc
from ptc_followup_common import PTC,require_ptc
from delivery import STAGE,REMOTE

def run():
    require_ptc();assert checked(OUT/'PTC_summary')
    assert (OUT/'PTC_summary/notebook_manifest.json').exists()
    files=[]
    def copy(p):
        target=STAGE/p.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
        files.append(str(target.relative_to(STAGE)))
    def bundle(paths,name,relative):
        target=STAGE/OUT.relative_to(ROOT)/'PTC_followups/artifacts'/name;target.parent.mkdir(parents=True,exist_ok=True)
        tmp=target.with_name(target.name+'.part')
        with tarfile.open(tmp,'w:gz',compresslevel=1) as arc:
            for p in sorted(set(paths)):arc.add(p,arcname=str(p.relative_to(relative)),recursive=False)
        tmp.replace(target);files.append(str(target.relative_to(STAGE)))
    copy(ROOT/'notebooks/dgscrna_results.ipynb')
    for p in (OUT/'PTC_summary').iterdir():
        if p.is_file() and p.name not in ['delivery_manifest.json','DELIVERY_RECEIPT.json','PTC_FULL_DELIVERED.json']:copy(p)
    for p in (OUT/'summary').glob('workflow_*'):copy(p)
    for p in CODE.glob('*.md'):copy(p)
    for directory in ['selection','verification','configurations']:
        for p in (PTC/directory).rglob('*'):
            if p.is_file():copy(p)
    selected={'predictions.csv.gz','terminal_manifest.json','training_manifest.json','training_history.json','manifest.json',
        'metrics.csv.gz','terminal_statuses.csv.gz','clusters.csv','initial_calls.csv.gz','cluster_calls.csv.gz',
        'marker_retention.csv.gz','score_manifest.json','prepare_manifest.json','cells.csv','PCA30.csv','UMAP2.csv',
        'DL_features.txt','scoring_features.txt','geometry_features.txt','markers_outside_eligible_CCA.txt',
        'SCORE_COMPLETE','COMPLETE','TERMINAL_COMPLETE','PREPARED','default_parity.json','retention_invariants.json','sessionInfo.txt'}
    for family in ['representation','marker_retention','MLP_seeds','default_parity','existing_grid_evaluation']:
        for d in (PTC/family).iterdir():
            if not d.is_dir():continue
            paths=[]
            for p in d.rglob('*'):
                if not p.is_file():continue
                if p.suffix in ['.png','.pdf']:copy(p)
                elif p.name in selected:paths.append(p)
            bundle(paths,f'{d.name}.tar.gz',PTC)
    source=Path(__file__).resolve().parent
    source_directories={source,Path(json.loads((PTC/'dispatch_state.json').read_text())['science_source'])}
    bundle([p for s in source_directories for p in s.iterdir() if p.suffix in ['.py','.R','.sbatch','.json']],
        'PTC_immutable_sources.tar.gz',OUT/'source_snapshots')
    copy(PTC/'dispatch_state.json')
    manifest=OUT/'PTC_summary/delivery_manifest.json'
    records=[dict(path=f,sha256=sha(STAGE/f),bytes=(STAGE/f).stat().st_size) for f in sorted(set(files))]
    write_json(manifest,dict(status='staged',remote=REMOTE,scope='PTC additions; completed GBM and historical results preserved',
        files=records,n_files=len(records),total_bytes=sum(r['bytes'] for r in records),job=os.environ['SLURM_JOB_ID']))
    copy(manifest)
    listing=PTC/'upload_files.txt';listing.write_text('\n'.join(sorted(set(files)))+'\n')
    subprocess.run(['rclone','copy',str(STAGE),REMOTE,'--files-from',str(listing),'--transfers','4','--checkers','8',
        '--log-file',str(PTC/'upload.log'),'--log-level','INFO'],check=True)
    subprocess.run(['rclone','check',str(STAGE),REMOTE,'--files-from',str(listing),'--download','--one-way','--checkers','4',
        '--log-file',str(PTC/'remote_check.log'),'--log-level','INFO'],check=True)
    receipt=OUT/'PTC_summary/DELIVERY_RECEIPT.json'
    write_json(receipt,dict(status='delivered_and_verified',scope='GBM and PTC accepted follow-up plan',remote=REMOTE,
        n_files=len(set(files)),manifest_sha256=sha(manifest),notebook_sha256=sha(ROOT/'notebooks/dgscrna_results.ipynb'),
        verification='rclone full-download one-way check and receipt readback',job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    remote=REMOTE+'/'+str(receipt.relative_to(ROOT));subprocess.run(['rclone','copyto',str(receipt),remote],check=True)
    with tempfile.TemporaryDirectory(dir=PTC,prefix='receipt_') as d:
        got=Path(d)/'receipt.json';subprocess.run(['rclone','copyto',remote,str(got)],check=True);assert sha(got)==sha(receipt)
    copy(receipt)

if __name__=='__main__':run()
