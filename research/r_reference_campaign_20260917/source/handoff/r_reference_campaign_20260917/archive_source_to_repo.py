"""Copy code/protocol only to align-r-reference; never copy patient/result matrices."""
import hashlib,json,shutil,subprocess
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
CODE=ROOT/'handoff/r_reference_campaign_20260917'
REPO=ROOT/'DGscRNA'
DEST=REPO/'research/r_reference_campaign_20260917'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def run():
    assert subprocess.check_output(['git','branch','--show-current'],cwd=REPO,text=True).strip()=='align-r-reference'
    suffixes={'.py','.R','.sbatch','.md','.json'}
    paths=[p for p in CODE.rglob('*') if p.is_file() and p.suffix in suffixes and '__pycache__' not in p.parts and p.name!='CURRENT_TASK.md']
    paths += [ROOT/'handoff/ptc_recovery_20260916'/n for n in ['refine.py','ptc_common.py','label_rules.py']]
    paths += [ROOT/'handoff/harmonize.py']
    records=[]
    for source in sorted(paths):
        relative=source.relative_to(ROOT);target=DEST/'source'/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
        assert sha(source)==sha(target)
        records.append(dict(source=str(relative),snapshot=str(target.relative_to(DEST)),sha256=sha(source),bytes=source.stat().st_size))
    for name in ['README.md','PLAN.md']:shutil.copy2(CODE/name,DEST/name)
    manifest=dict(branch='align-r-reference',scope='Added research/r_reference_campaign_20260917 only; installed Python package and earlier research snapshot unchanged.',
      source_files=len(records),files=records,
      dependencies=['Prior source archive: research/gbm_ptc_20260917/reference/ptc_archive/source.R and source.py',
                    'R dbscan 64-bit MST-index patch described in the prior reproducibility archive',
                    'Site-specific R/Python environments, raw data and cached CellMarker workbook (not included)'],
      excluded=['patient metadata','per-cell labels','expression/anchor matrices','model weights','rendered notebooks',
                'result tables','credentials','OneDrive files','transient current-task notes'],
      portability='Preserves executed absolute HPC paths; research execution archive, not a portable standalone package release.')
    (DEST/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print('ARCHIVED_CODE',len(records),DEST,flush=True)

if __name__=='__main__':run()
