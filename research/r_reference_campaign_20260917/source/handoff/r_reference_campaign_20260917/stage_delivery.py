"""Package derived results in the already authorized OneDrive directory layout."""
import ast,fcntl,gzip,hashlib,json,os,shutil,subprocess,tarfile
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'
CODE=ROOT/'handoff/r_reference_campaign_20260917'
PARENT=ROOT/'results/hvg_ptc_20260916_v1/onedrive_existing_results_20260916'
STAGE=PARENT/'GBM_PTC_results_20260916'
REMOTE='onedrive:work_od/share/dgscrna_GSE274546_TKU3186/01_report/GBM_PTC_results_20260916'

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import pandas as pd
    summary=json.loads((OUT/'summary/campaign_summary.json').read_text());assert summary['status']=='complete'
    nb=json.loads((OUT/'summary/notebook_update_manifest.json').read_text());assert nb['status']=='complete'
    assert sha(ROOT/'notebooks/dgscrna_results.ipynb')==nb['current_sha256']
    inv=pd.read_csv(OUT/'summary/completion_inventory.csv');assert len(inv)==99 and inv.audited.all()
    design=json.loads((OUT/'verification/design_audit.json').read_text());assert len(design['comparisons'])==14
    early=list((OUT/'verification/pre_guard_density').glob('*.json'));assert len(early)==16
    assert all(json.loads(p.read_text())['status']=='passed' for p in early)
    code_files=sorted(p for p in CODE.rglob('*') if p.is_file() and '__pycache__' not in str(p) and p.name!='CURRENT_TASK.md' and not p.name.endswith(('.pyc','.next')))
    for p in code_files:
        if p.suffix=='.py':ast.parse(p.read_text(),filename=str(p))
        if p.suffix=='.sbatch':subprocess.run(['bash','-n',str(p)],check=True)
    rs=[str(p) for p in code_files if p.suffix=='.R']
    subprocess.run(['/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript','-e',
       'for(p in commandArgs(trailingOnly=TRUE))parse(p)',*rs],check=True)
    checked=dict(status='passed',python_files=sum(p.suffix=='.py' for p in code_files),
      R_files=len(rs),shell_files=sum(p.suffix=='.sbatch' for p in code_files),job=os.environ['SLURM_JOB_ID'])
    (OUT/'verification/source_syntax_checks.json').write_text(json.dumps(checked,indent=2)+'\n')
    source_index=[dict(path=str(p.relative_to(ROOT)),sha256=sha(p)) for p in code_files]
    (OUT/'summary/source_manifest.json').write_text(json.dumps(source_index,indent=2)+'\n')
    files=[];artifacts=[]
    def copy(p):
        rel=p.relative_to(ROOT);target=STAGE/rel;target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(p,target);files.append(str(rel))
    def bundle(paths,target,relative_root):
        target.parent.mkdir(parents=True,exist_ok=True)
        temp=target.with_name(target.name+'.part')
        with tarfile.open(temp,'w:gz',compresslevel=1) as archive:
            for p in sorted(paths):archive.add(p,arcname=str(p.relative_to(relative_root)),recursive=False)
        temp.replace(target);files.append(str(target.relative_to(STAGE)))
    copy(ROOT/'notebooks/dgscrna_results.ipynb')
    for p in code_files:copy(p)
    for rel in ['handoff/ptc_recovery_20260916/refine.py','handoff/ptc_recovery_20260916/ptc_common.py',
                'handoff/ptc_recovery_20260916/label_rules.py','handoff/harmonize.py','handoff/deck_datasets_provenance.md']:
        copy(ROOT/rel)
    for sub in ['summary','markers','verification']:
        for p in sorted((OUT/sub).rglob('*')):
            if p.is_file() and p.name not in ['DELIVERY_RECEIPT.json']:copy(p)
    for p in OUT.glob('*.json'):
        if p.name not in ['dispatch_state.json']:copy(p)
    for p in OUT.glob('*.txt'):copy(p)
    for p in OUT.glob('*events.jsonl'):copy(p)
    for dataset_dir in sorted((OUT/'inputs').iterdir()):
        if dataset_dir.is_dir():
            selected=[p for p in dataset_dir.rglob('*') if p.is_file() and p.suffix not in ['.bin']]
            target=STAGE/dataset_dir.relative_to(ROOT)/'evaluation_inputs_and_manifests.tar.gz'
            bundle(selected,target,dataset_dir)
    for row in inv.itertuples(index=False):
        prep=Path(row.directory)
        paths=[]
        for p in prep.rglob('*'):
            if not p.is_file():continue
            relative=p.relative_to(prep);parts=relative.parts
            # Dense expression/anchor inputs and model weights remain on HPC.
            if p.name in ['expression_PCA30.rds','RNA_normalized.rds','anchors.rds','DL.float32.bin','model_state.pt','terminal.npz']:
                continue
            if p.name.endswith(('.native_order_backup','.part')):continue
            if 'figures' in parts or 'verification' in parts:
                copy(p);continue
            paths.append(p)
            if len(parts)==1 or (parts[0]=='evaluation' and len(parts)==2):copy(p)
        target=STAGE/prep.relative_to(ROOT)/'terminal_results_and_audits.tar.gz'
        bundle(paths,target,prep)
        for p in prep.glob('*/terminal/*/training_manifest.json'):
            m=json.loads(p.read_text())
            for name in ['terminal.npz','model_state.pt']:
                if name in m['outputs']:artifacts.append(dict(path=str(p.parent/name),sha256=m['outputs'][name],
                    availability='retained on HPC; terminal label CSV and manifest included in per-unit delivery archive'))
        print('PACKAGED',row.unit,flush=True)
    # Executed per-job sources are compactly archived instead of thousands of loose files.
    bundle([p for p in (OUT/'execution_sources').rglob('*') if p.is_file()],
      STAGE/OUT.relative_to(ROOT)/'executed_sources.tar.gz',OUT/'execution_sources')
    local_record=OUT/'summary/local_retained_model_artifacts.jsonl.gz'
    with gzip.open(local_record,'wt') as f:
        for r in artifacts:f.write(json.dumps(r)+'\n')
    copy(local_record)
    layout=OUT/'summary/DELIVERY_LAYOUT.md'
    layout.write_text('''# Delivery layout

The existing `notebooks/dgscrna_results.ipynb` includes all prior cells plus the executed campaign section.
`summary/` contains the full grid, paired ablation differences, marker roster, label mappings and interpretation.
Each analysis unit has directly viewable figures, aggregate evaluation tables and an artifact-audit report.
Its `terminal_results_and_audits.tar.gz` contains all initial/final per-cell annotation CSVs, confidence values,
training histories/manifests, cluster assignments, DEG/density outputs and per-condition confusion matrices.
Extract that archive inside the unit directory to restore the detailed result tree.

Models and probability NPZ files remain on HPC, with original paths/checksums in
`local_retained_model_artifacts.jsonl.gz`; exported final labels for every arm are included in the archives.
Dense expression/anchor matrices and raw h5ad files are not duplicated into OneDrive.
Input sources and checksums remain in per-cohort input-manifest archives. Executed source copies are in
`executed_sources.tar.gz`; reviewable workflow code is under `handoff/r_reference_campaign_20260917/`.
No website or replacement-version notebook is generated, and no old result files are deleted.
''')
    copy(layout)
    files=sorted(set(files))
    listpath=PARENT/'R_reference_campaign_20260917_upload_files.txt';listpath.write_text('\n'.join(files)+'\n')
    manifest=[dict(path=rel,bytes=(STAGE/rel).stat().st_size,sha256=sha(STAGE/rel)) for rel in files]
    manifest_path=PARENT/'R_reference_campaign_20260917_delivery_manifest.json'
    manifest_path.write_text(json.dumps(dict(remote=REMOTE,status='staged',n_files=len(files),
       total_bytes=sum(r['bytes'] for r in manifest),files=manifest,job=os.environ['SLURM_JOB_ID']),indent=2)+'\n')
    print('STAGED',len(files),'files',sum(r['bytes'] for r in manifest),'bytes',flush=True)

if __name__=='__main__':run()
