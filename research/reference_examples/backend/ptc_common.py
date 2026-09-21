"""PTC experiment IO and provenance; importing this module performs no science."""
from pathlib import Path
import hashlib
import json
import os
from datetime import datetime, timezone

ROOT = Path(os.environ.get('DGSCRNA_SITE_ROOT',str(Path.cwd())))
BASE = Path(os.environ['DGSCRNA_EXAMPLE_OUT'])
RECOVERY = ROOT/'results/hvg_ptc_20260916_v1/ptc_recovery'
PYTHON = os.environ.get('DGSCRNA_PYTHON','python')
RSCRIPT = os.environ.get('DGSCRNA_RSCRIPT','Rscript')
GROUPS = {'MTN': ['MT-1', 'MT-2', 'N-1', 'N-2'], 'TUT': ['TU-1', 'TU-2', 'T-1', 'T-2']}
MARKERS = RECOVERY/'archive/tcr/ptc_val/scripts/DGscRNA-Share/data/full_marker_symbol.RDS'

def utc():
    return datetime.now(timezone.utc).isoformat()

def require_slurm():
    assert os.environ.get('SLURM_JOB_ID'), 'Scientific computation requires SLURM'
    assert not os.environ.get('PYTHONOPTIMIZE', '0') != '0'

def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def write_json(path, obj):
    import math
    def safe(value):
        if isinstance(value, float) and not math.isfinite(value):
            return 'NaN' if math.isnan(value) else ('+Infinity' if value>0 else '-Infinity')
        if isinstance(value, dict):
            return {str(k): safe(v) for k,v in value.items()}
        if isinstance(value, (list,tuple)):
            return [safe(v) for v in value]
        return value
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name+'.part')
    temp.write_text(json.dumps(safe(obj), indent=2, ensure_ascii=False, allow_nan=False)+'\n')
    temp.replace(path)

def task_list():
    return json.loads((BASE/'protocol/single_sample_geometries.json').read_text())

def selected_task():
    import sys
    index = int(sys.argv[1] if len(sys.argv)>1 else os.environ['SLURM_ARRAY_TASK_ID'])
    return task_list()[index]

def geometry_dir(task):
    return BASE/'single_sample'/task['sample']/task['geometry_id']
