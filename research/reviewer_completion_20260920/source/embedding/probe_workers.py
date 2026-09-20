"""Measure four-worker original-R scoring and require exact two-worker parity."""
from pathlib import Path
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
import shutil
import subprocess

assert os.environ.get('SLURM_JOB_ID')
ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
CAMP = ROOT / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
SOURCE = CAMP / 'source_snapshots/embedding_v6'
ORIGINAL = CAMP / 'embedding/SN040/hvg2000/PCA2'
OUT = CAMP / 'embedding_worker_probe/SN040/hvg2000/PCA2_workers4'
R = '/fs/scratch/PCON0080/yimin/mamba_envs/deconv_r2/bin/Rscript'


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


OUT.mkdir(parents=True, exist_ok=False)
config = json.loads((ORIGINAL / 'config.json').read_text())
condition = next(c.copy() for c in config['conditions'] if c['method']=='KMeans' and c['k']==5)
old = Path(condition['dest'])
original_manifest = json.loads((old / 'score_manifest.json').read_text())
assert original_manifest['DEG_workers'] == 2
assert (old / 'SCORE_COMPLETE').read_text().strip() == sha(old / 'score_manifest.json')
dest = OUT / old.name;dest.mkdir()
shutil.copy2(old / 'clusters.csv', dest / 'clusters.csv')
condition['dest'] = str(dest)
config.update(dest=str(OUT), conditions=[condition],
              purpose='Operational worker-count resource probe; unchanged scientific parameters and numerical source')
(OUT / 'config.json').write_text(json.dumps(config, indent=2)+'\n')
env = os.environ.copy();env['DGSCRNA_DEG_WORKERS']='4'
subprocess.run([R, str(SOURCE / 'score_candidates.R'), str(OUT / 'config.json')],env=env,check=True)
current_manifest = json.loads((dest / 'score_manifest.json').read_text())
assert current_manifest['DEG_workers'] == 4
assert (dest / 'SCORE_COMPLETE').read_text().strip() == sha(dest / 'score_manifest.json')
assert current_manifest['source_sha256'] == original_manifest['source_sha256']
for filename in ['initial_calls.csv.gz', 'cluster_calls.csv.gz', 'marker_retention.csv.gz']:
    with gzip.open(dest / filename, 'rb') as a, gzip.open(old / filename, 'rb') as b:
        assert a.read() == b.read(), f'Worker-count numerical mismatch: {filename}'
subprocess.run([R,'-e',
    "a<-commandArgs(TRUE);for(f in c('DEG.rds','density_scores.rds')){x<-readRDS(file.path(a[1],f));y<-readRDS(file.path(a[2],f));stopifnot(identical(x,y))};cat('DEG_AND_DENSITY_EXACT_PARITY\\n')",
    str(dest),str(old)],env=env,check=True)
cg = Path('/sys/fs/cgroup') / Path('/proc/self/cgroup').read_text().strip().split('::',1)[1].lstrip('/')
job_group = next(p for p in [cg,*cg.parents] if p.name == 'job_'+os.environ['SLURM_JOB_ID'])
peak = int((job_group / 'memory.peak').read_text())
record = dict(status='passed', completed_at=datetime.now(timezone.utc).isoformat(),
    job=os.environ['SLURM_JOB_ID'], source_sha256=sha(__file__), numeric_scorer_sha256=sha(SOURCE / 'score_candidates.R'),
    sample='SN040', n_cells=10135, n_features_scoring=26230, budget='hvg2000', candidate='PCA2_KMeans_K05',
    workers_before=2, workers_after=4, original_score_seconds=original_manifest['elapsed_seconds'],
    workers4_score_seconds=current_manifest['elapsed_seconds'],
    score_seed_and_retention_tables_byte_equal_after_decompression=True,
    DEG_and_density_R_objects_identical=True, author_truth_used=False,
    job_cgroup_peak_bytes=peak, job_cgroup_path=str(job_group),
    memory_scope='Whole isolated SLURM job cgroup, including all R child processes and charged page cache; not batch MaxRSS.',
    source_manifests={str(old / 'score_manifest.json'):sha(old / 'score_manifest.json'),str(dest / 'score_manifest.json'):sha(dest / 'score_manifest.json')},
    operational_acceptance=bool(peak < 28*1024**3 and current_manifest['elapsed_seconds'] < original_manifest['elapsed_seconds']),
    limitation='One largest-cell sample / one five-cluster candidate; larger K can alter peak memory and runtime.')
(OUT / 'validation.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record),flush=True)
