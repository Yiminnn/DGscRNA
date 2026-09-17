"""Shared provenance/IO for the user-approved HVG experiment; no compute on import."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import sys

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT = ROOT / 'results/hvg_ptc_20260916_v1'
VENDOR = OUT / 'vendor'
if VENDOR.is_dir():
    sys.path.insert(0, str(VENDOR))
INPUTS = ROOT / 'results/g274_cohort/inputs_v1'
FEATURES = ['all', 'hvg500', 'hvg1000', 'hvg2000', 'hvg3000', 'hvg5000']
EXTRA_FEATURES = ['seurat2000', 'seurat5000', 'hvg2000_markers']
SEEDS = [42, 7, 17, 29, 101]
MARKER = 'CM2_glioma_other'
MARKERS = ROOT / 'handoff/markers_v3'


def require_slurm():
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Scientific computation requires a SLURM allocation')
    if os.environ.get('PYTHONOPTIMIZE', '0') != '0':
        raise RuntimeError('Assertions must remain enabled')


def utc():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def clean(x):
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [clean(v) for v in x]
    if hasattr(x, 'tolist'):
        return clean(x.tolist())
    if isinstance(x, float):
        import math
        return x if math.isfinite(x) else None
    if isinstance(x, Path):
        return str(x)
    return x


def write_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f'.tmp.{os.getpid()}')
    tmp.write_text(json.dumps(clean(data), indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    tmp.replace(path)


def key(data):
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(',', ':')).encode()).hexdigest()[:16]


def version_record():
    from importlib.metadata import version, PackageNotFoundError
    import sys
    d = {'python': sys.version}
    for name in ['numpy', 'scipy', 'pandas', 'scanpy', 'anndata', 'scikit-learn',
                 'scikit-misc', 'umap-learn', 'hdbscan', 'torch', 'dgscrna']:
        try:
            d[name] = version(name)
        except PackageNotFoundError:
            d[name] = 'not_installed'
    return d


def runtime_record():
    return {'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
            'array_job': os.environ.get('SLURM_ARRAY_JOB_ID'),
            'array_task': os.environ.get('SLURM_ARRAY_TASK_ID'),
            'host': os.uname().nodename,
            'cpus': os.environ.get('SLURM_CPUS_PER_TASK')}


def samples():
    return (ROOT / 'handoff/g274/cohort_samples.txt').read_text().split()
