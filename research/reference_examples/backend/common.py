"""IO and frozen paths; imports perform no scientific computation."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os

ROOT = Path(os.environ.get('DGSCRNA_SITE_ROOT',str(Path.cwd())))
CODE = Path(__file__).resolve().parent
OUT = Path(os.environ['DGSCRNA_EXAMPLE_OUT'])
INPUTS = ROOT/'results/g274_cohort/inputs_v1'
OLD = ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'
PYTHON = os.environ.get('DGSCRNA_PYTHON','python')
RSCRIPT = os.environ.get('DGSCRNA_RSCRIPT','Rscript')
FEATURES = ['hvg500', 'hvg1000', 'hvg2000', 'hvg3000', 'hvg5000', 'all']
ROUTES = ['PCA30_SNN', 'PCA30_HDBSCAN_R', 'UMAP2_SNN', 'UMAP2_HDBSCAN_R']
L1 = ['Malignant', 'TAM', 'Lymphocyte', 'Oligodendrocyte', 'Astrocyte', 'OPC',
      'Excitatory neuron', 'Inhibitory neuron', 'Endothel', 'Pericyte', 'Other']

def require_slurm():
    assert os.environ.get('SLURM_JOB_ID'), 'All scientific computation requires SLURM'
    assert not os.environ.get('PYTHONOPTIMIZE', '0') != '0'

def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()

def utc():
    return datetime.now(timezone.utc).isoformat()

def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f'.part.{os.getpid()}')
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False)+'\n')
    tmp.replace(path)

def samples():
    return (ROOT/'handoff/g274/cohort_samples.txt').read_text().split()

def complete(dest, manifest='manifest.json', flag='COMPLETE'):
    (Path(dest)/flag).write_text(sha(Path(dest)/manifest)+'\n')

def checked(dest, manifest='manifest.json', flag='COMPLETE'):
    p = Path(dest)
    return (p/flag).exists() and (p/flag).read_text().strip() == sha(p/manifest)
