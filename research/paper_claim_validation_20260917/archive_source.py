"""Archive reviewed execution code only on the user-selected repository branch."""
import json
from pathlib import Path
import shutil
import subprocess
from common import ROOT, CODE, sha, utc

def run():
    repo=ROOT/'DGscRNA'
    assert subprocess.check_output(['git','branch','--show-current'],cwd=repo,text=True).strip()=='align-r-reference'
    dest=repo/'research/paper_claim_validation_20260917';dest.mkdir(parents=True,exist_ok=True)
    files=[p for p in CODE.iterdir() if p.suffix in ['.py','.R','.sbatch','.md','.json'] and p.name!='CURRENT_TASK.md']
    files.extend(p for p in (CODE/'vendor_sources').iterdir() if p.is_file() and p.name in
        ['sctype_score_original.R','scType_LICENSE','scCATCH_LICENSE.md','SOURCE_PROVENANCE.json'])
    records=[]
    for p in files:
        target=dest/p.relative_to(CODE);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
        assert sha(target)==sha(p)
        records.append(dict(path=str(p.relative_to(CODE)),sha256=sha(p)))
    write=dict(scope='Only research/paper_claim_validation_20260917; existing package, other research archives and other branches unchanged.',
        branch='align-r-reference',status='execution_source; full experimental completion is not implied',
        files=records,archived_at=utc(),
        excluded=['data','per-cell labels','patient metadata','model weights','environments','credentials','rendered notebooks'],
        dependencies='Absolute site paths preserved. Original data and pinned comparator source archives are installed outside Git. See README and vendor source provenance.')
    (dest/'SOURCE_MANIFEST.json').write_text(json.dumps(write,indent=2)+'\n')
    print(dest,len(files))

if __name__=='__main__':run()
