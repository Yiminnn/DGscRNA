#!/bin/bash
#SBATCH --job-name=claim_scDeepSort_install
#SBATCH --account=PCON0080
#SBATCH --qos=ascend-default
#SBATCH --partition=nextgen
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/fs/scratch/PCON0080/yimin/dgscrna/logs/claim_scDeepSort_install_%j.out
set -euo pipefail
claim_root=/fs/scratch/PCON0080/yimin/dgscrna
claim_env=$claim_root/results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/vendor_envs/scDeepSort
export UV_PYTHON_INSTALL_DIR=$claim_root/results/hvg_ptc_20260916_v1/paper_claim_validation_20260917/vendor_envs/managed_python
uv_bin=/users/PCON0080/yimin/.local/bin/uv
if [ ! -x "$claim_env/bin/python" ]; then
  "$uv_bin" venv --python 3.8 --python-preference only-managed "$claim_env"
fi
"$uv_bin" pip install --python "$claim_env/bin/python" numpy==1.23.5 scipy==1.10.1 pandas==1.5.3 scikit-learn==1.2.2 dgl==0.4.3.post2 xlrd==1.2.0 openpyxl==3.1.5 tqdm==4.67.1
"$uv_bin" pip install --python "$claim_env/bin/python" --index-url https://download.pytorch.org/whl/cpu torch==1.13.1+cpu
"$claim_env/bin/python" - <<'PY'
from pathlib import Path
import sys,sysconfig,shutil,json,os,hashlib
base=Path('/fs/scratch/PCON0080/yimin/dgscrna/results/hvg_ptc_20260916_v1/paper_claim_validation_20260917')
env=base/'vendor_envs/scDeepSort'
original=Path('/fs/scratch/PCON0080/yimin/tools/deepsort-1.0')
dest=Path(sysconfig.get_paths()['purelib'])/'deepsort'
shutil.copytree(original/'deepsort',dest,dirs_exist_ok=True)
link=env/'deepsort-pretrained'
if not link.exists():link.symlink_to(original/'deepsort-pretrained',target_is_directory=True)
assert link.resolve()==(original/'deepsort-pretrained').resolve()
import torch,dgl
from deepsort import DeepSortPredictor
model=DeepSortPredictor(species='human',tissue='Brain')
files={str(p.relative_to(original)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (original/'deepsort').glob('*.py')}
(base/'scDeepSort_install.json').write_text(json.dumps(dict(status='installed_and_checkpoint_loaded',job=os.environ['SLURM_JOB_ID'],python=sys.version,torch=torch.__version__,dgl=dgl.__version__,source_files=files,pretrained_source=str(link.resolve())),indent=2)+'\n')
PY
