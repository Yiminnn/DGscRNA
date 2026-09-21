"""Record actual numerical package paths, versions and CPU Torch availability."""
from pathlib import Path
import json
import os
import sys
assert os.environ.get('SLURM_JOB_ID')
import numpy
import pandas
import scipy
import torch
prefix = Path(sys.argv[1]).resolve()
assert Path(sys.prefix).resolve() == prefix
packages = {m.__name__: {'version': m.__version__, 'path': str(Path(m.__file__).resolve())}
            for m in [numpy, pandas, scipy, torch]}
assert all(Path(item['path']).is_relative_to(prefix) for item in packages.values())
torch.set_num_threads(4)
assert torch.isfinite(torch.zeros((2, 2)) @ torch.ones((2, 2))).all()
Path(sys.argv[2]).write_text(json.dumps(dict(status='all_numerical_packages_in_new_prefix',
    prefix=str(prefix), python=sys.version, packages=packages, torch_CPU_test=True,
    job=os.environ['SLURM_JOB_ID']), indent=2)+'\n')
