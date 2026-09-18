"""Add GBM follow-ups to the same authorized OneDrive directory and verify bytes."""
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import tarfile
import tempfile
from common import ROOT,CODE,OUT,require_slurm,checked,sha,write_json,utc
from delivery import STAGE,REMOTE

def run():
    require_slurm();assert checked(OUT/'GBM_full_summary')
    files=[]
    def copy(source):
        target=STAGE/source.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source,target);files.append(str(target.relative_to(STAGE)))
    def bundle(paths,relative,name):
        target=STAGE/OUT.relative_to(ROOT)/name;target.parent.mkdir(parents=True,exist_ok=True)
        temporary=target.with_name(target.name+'.part')
        with tarfile.open(temporary,'w:gz',compresslevel=1) as arc:
            for path in sorted(set(paths)):arc.add(path,arcname=str(path.relative_to(relative)),recursive=False)
        temporary.replace(target);files.append(str(target.relative_to(STAGE)))
    copy(ROOT/'notebooks/dgscrna_results.ipynb')
    folders=['controls_summary','comparison_summary','unknown_summary','unknown_expression','marker_evidence_summary','scalability_summary','legacy_coverage_audit','GBM_full_summary']
    for name in folders:
        for p in (OUT/name).rglob('*'):
            if p.is_file() and p.name not in ['delivery_manifest.json','DELIVERY_RECEIPT.json','REMOTE_RECEIPT_UPLOADED.json','GBM_FULL_DELIVERED.json']:copy(p)
    for p in (OUT/'summary').glob('workflow_*'):copy(p)
    for p in CODE.glob('*.md'):copy(p)
    snapshots=[p for p in (OUT/'source_snapshots').rglob('*') if p.is_file() and p.suffix in ['.py','.R','.sbatch','.json']]
    bundle(snapshots,OUT/'source_snapshots','source/full_GBM_executed_source_snapshots.tar.gz')
    source=Path(__file__).resolve().parent
    bundle([p for p in source.iterdir() if p.suffix in ['.py','.R','.sbatch']],source,'source/full_GBM_finalizer_source.tar.gz')
    selected={'predictions.csv.gz','terminal_manifest.json','training_manifest.json','training_history.json','manifest.json',
        'fit_manifest.json','cohort_manifest.json','cohort_predictions.csv.gz','metrics.csv','per_class.csv.gz','confusions.csv.gz',
        'clusters.csv','initial_calls.csv.gz','config.json','design.json','score_manifest.json','prepare_manifest.json','SCORE_COMPLETE',
        'COMPLETE','TERMINAL_COMPLETE','COHORT_COMPLETE','control_feature_genes.txt','control_embedding.csv'}
    cohort=list(csv.DictReader((OUT/'protocol/cohort.csv').open()))
    for row in cohort:
        sample=row['sample'];paths=[]
        for method in ['scType','scCATCH','SCINA','SingleR','scDeepSort']:
            for p in (OUT/'comparators'/method/sample).rglob('*'):
                if p.is_file() and p.name in selected:paths.append(p)
        for family in ['GBM_DL_controls','GBM_representation_controls']:
            for p in (OUT/family/sample).rglob('*'):
                if not p.is_file():continue
                if p.suffix in ['.png','.pdf']:copy(p)
                elif p.name in selected:
                    if p.name=='control_embedding.csv' and p.parent.name.startswith('RNA_noDR'):continue
                    paths.append(p)
        for budget in ['hvg2000','hvg5000','all']:
            paths.extend(p for p in (OUT/'GBM'/sample/budget/'evaluation_geometry_only_DL2000').iterdir() if p.is_file())
            for p in (OUT/'GBM'/sample/budget).glob('*/terminal_geometry_only_DL2000/*/*'):
                if p.is_file() and p.name in selected:paths.append(p)
        bundle(paths,OUT,f'artifacts/{sample}_controls_and_comparators.tar.gz')
    resources=[]
    for p in (OUT/'scalability').rglob('*'):
        if p.is_file() and p.suffix in ['.png','.pdf']:copy(p)
        elif p.is_file() and (p.name in selected or p.suffix=='.json'):resources.append(p)
    bundle(resources,OUT/'scalability','artifacts/single_run_resource_artifacts.tar.gz')
    resources=[]
    for p in (OUT/'scalability_repeats').rglob('*'):
        if p.is_file() and p.suffix in ['.png','.pdf']:copy(p)
        elif p.is_file() and (p.name in selected or p.suffix=='.json'):resources.append(p)
    bundle(resources,OUT/'scalability_repeats','artifacts/resource_repetition_artifacts.tar.gz')
    for name in ['submissions.jsonl','dispatch_state.json','aux_dispatch_state.json','representation_dispatch_state.json','scalability_dispatch_state.json','SCINA_library_dispatch_state.json','analysis_dispatch_state.json']:
        copy(OUT/name)
    for p in (OUT/'verification').rglob('*'):
        if p.is_file() and p.suffix in ['.json','.csv','.md','.txt']:copy(p)
    records=[dict(path=n,sha256=sha(STAGE/n),bytes=(STAGE/n).stat().st_size) for n in sorted(set(files))]
    manifest=OUT/'GBM_full_summary/delivery_manifest.json'
    write_json(manifest,dict(status='staged',scope='GBM follow-ups; existing core and historical results retained',
        n_files=len(records),total_bytes=sum(r['bytes'] for r in records),files=records,remote=REMOTE,job=os.environ['SLURM_JOB_ID']))
    copy(manifest)
    listing=OUT/'full_GBM_upload_files.txt';listing.write_text('\n'.join(sorted(set(files)))+'\n')
    subprocess.run(['rclone','copy',str(STAGE),REMOTE,'--files-from',str(listing),'--transfers','4','--checkers','8',
                    '--log-file',str(OUT/'full_GBM_upload.log'),'--log-level','INFO'],check=True)
    subprocess.run(['rclone','check',str(STAGE),REMOTE,'--files-from',str(listing),'--download','--one-way','--checkers','4',
                    '--log-file',str(OUT/'full_GBM_remote_check.log'),'--log-level','INFO'],check=True)
    receipt=OUT/'GBM_full_summary/DELIVERY_RECEIPT.json'
    write_json(receipt,dict(status='delivered_and_verified',scope='GBM completed; PTC follow-ups pending',remote=REMOTE,
        n_files=len(set(files)),manifest_sha256=sha(manifest),notebook_sha256=sha(ROOT/'notebooks/dgscrna_results.ipynb'),
        verification='rclone copy and full-download one-way check exited zero',job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    remote=REMOTE+'/'+str(receipt.relative_to(ROOT))
    subprocess.run(['rclone','copyto',str(receipt),remote],check=True)
    with tempfile.TemporaryDirectory(dir=OUT,prefix='full_GBM_receipt_') as directory:
        got=Path(directory)/'receipt.json';subprocess.run(['rclone','copyto',remote,str(got)],check=True);assert sha(got)==sha(receipt)
    copy(receipt)
    write_json(OUT/'GBM_full_summary/REMOTE_RECEIPT_UPLOADED.json',dict(status='verified',receipt_sha256=sha(receipt),completed_at=utc()))

if __name__=='__main__':run()
