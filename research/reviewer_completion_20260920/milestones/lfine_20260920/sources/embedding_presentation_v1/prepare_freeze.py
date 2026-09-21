"""SLURM metadata-only freeze; does not activate jobs or modify old artifacts."""
from pathlib import Path
import ast,hashlib,json,os,shutil
assert os.environ.get('SLURM_JOB_ID')
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
CODE=Path(__file__).resolve().parent
FROZEN=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/source_snapshots/embedding_presentation_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for p in CODE.glob('*.py'):ast.parse(p.read_text(),filename=str(p))
base=ROOT/'handoff/reviewer_completion_20260920/embedding/select_and_summarize.py'
def funcs(p):
    source=p.read_text();return {n.name:ast.get_source_segment(source,n) for n in ast.parse(source).body if isinstance(n,ast.FunctionDef)}
a,b=funcs(base),funcs(CODE/'select_and_summarize.py')
for name in ('verify_evaluation','assert_primary_boolean'):assert a[name]==b[name]
assert not FROZEN.exists();FROZEN.mkdir()
for p in CODE.iterdir():
    if p.is_file():shutil.copy2(p,FROZEN/p.name)
manifest={p.name:sha(p) for p in FROZEN.iterdir() if p.is_file()}
(FROZEN/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(dict(status='frozen_not_activated',path=str(FROZEN),manifest_sha256=sha(FROZEN/'SOURCE_MANIFEST.json'),
    preserved_evaluation_and_boolean_guards=True,job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'))),flush=True)
