"""Stage the correction in the already-authorized existing OneDrive directory."""
from pathlib import Path
import json,hashlib,shutil
from datetime import datetime,timezone
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
BASE=ROOT/'results/hvg_ptc_20260916_v1'
OUT=BASE/'ptc_paper_baseline'
DELIVERY=BASE/'onedrive_existing_results_20260916'
STAGE=DELIVERY/'GBM_PTC_results_20260916'

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def run():
    # File copies, metadata and checksums only; scientific outputs already ran in SLURM.
    notebook=json.loads((OUT/'notebook_correction_manifest.json').read_text())
    assert sha(ROOT/'notebooks/dgscrna_results.ipynb')==notebook['current_sha256']
    delivery_path=BASE/'EXPERIMENT_DELIVERY.json'
    m=json.loads(delivery_path.read_text())
    m['completed_scope']='Historical GBM and two-group PTC experiment census; not a passed original-paper PTC gate'
    m['PTC_paper_baseline_status']='under_reconciliation; all eight DG paper F1/AUC calculations recovered; terminal refit differences and Accuracy discrepancy documented'
    m['PTC_paper_baseline_report']=str(OUT/'STATUS.md')
    m['notebook_baseline_prior_sha256']=notebook['backup_sha256']
    m['notebook_sha256']=notebook['current_sha256']
    m['notebook_updated_at']=notebook['updated_utc']
    m['notebook_baseline_correction_manifest']=str(OUT/'notebook_correction_manifest.json')
    delivery_path.write_text(json.dumps(m,indent=2)+'\n')
    paths=[ROOT/'notebooks/dgscrna_results.ipynb',delivery_path,
        ROOT/'handoff/plan_20260916_hvg_ptc/PTC_BASELINE_GATE.md',
        BASE/'original_R_workflow_review/REVIEW_ZH.md',BASE/'ptc_experiments/PTC_REPORT.md',
        ROOT/'results/ptc/S2_RECONCILIATION.md',
        BASE/'ptc_recovery/inventory/notebook_sources/06_ptc_paper.txt',
        BASE/'ptc_recovery/inventory_workspace/object_01.commands.txt',
        BASE/'ptc_recovery/archive/tcr/scripts/annotations.xlsx',
        BASE/'ptc_recovery/archive/tcr/rawdata/data_with_validation.csv']
    paths.extend(p for p in OUT.iterdir() if p.is_file() and p.suffix!='.ipynb')
    for name in ['evaluation_marker_union','replay_selected_routes_marker_union']:
        paths.extend(p for p in (OUT/name).rglob('*') if p.is_file() and p.suffix!='.bin')
    if (OUT/'evaluation_full_parallel/manifest.json').exists():
        for name in ['evaluation_full_parallel','replay_selected_routes_full_parallel']:
            paths.extend(p for p in (OUT/name).rglob('*') if p.is_file() and p.suffix!='.bin')
    if (OUT/'literal_original_DL').exists():
        paths.extend(p for p in (OUT/'literal_original_DL').rglob('*')
                     if p.is_file() and p.suffix != '.h5ad')
    if (OUT/'original_writing_note/manifest.json').exists():
        paths.extend(p for p in (OUT/'original_writing_note').iterdir() if p.is_file())
    paths.extend(p for p in (ROOT/'handoff/ptc_paper_baseline_20260916').rglob('*')
                 if p.is_file() and '__pycache__' not in p.parts)
    records=[]
    for source in sorted(set(paths)):
        rel=source.relative_to(ROOT);target=STAGE/rel
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(source,target)
        h=sha(source);assert h==sha(target)
        records.append(dict(path=str(rel),bytes=source.stat().st_size,sha256=h))
    (DELIVERY/'PTC_baseline_correction_upload_files.txt').write_text('\n'.join(r['path'] for r in records)+'\n')
    report=dict(staged_utc=datetime.now(timezone.utc).isoformat(),files=records,
        n_files=len(records),total_bytes=sum(r['bytes'] for r in records),
        remote='onedrive:work_od/share/dgscrna_GSE274546_TKU3186/01_report/GBM_PTC_results_20260916',
        scope='Original PTC baseline correction in existing notebook and directory; no raw expression matrices duplicated')
    (DELIVERY/'PTC_baseline_correction_staging.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='files'},indent=2))

if __name__=='__main__':run()
