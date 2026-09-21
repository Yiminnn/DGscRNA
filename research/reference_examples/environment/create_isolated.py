"""Install exact package locks with Conda's explicit isolated configuration search path."""
from pathlib import Path
import json
import os
import sys

assert os.environ.get('SLURM_JOB_ID'), 'Environment installation requires SLURM'
root = Path(sys.argv[1]).resolve()
here = Path(__file__).resolve().parent
cache = root/'package_cache'
os.environ['CONDA_PKGS_DIRS'] = str(cache)
os.environ['CONDA_SAFETY_CHECKS'] = 'enabled'
os.environ['CONDA_ALWAYS_COPY'] = 'true'
from conda.base.context import context, reset_context
from conda.cli.python_api import run_command, Commands
import conda

# An env-variable list alone merges with ~/.condarc and can still read its caches.
# The documented API search_path=() excludes site/user configuration without editing it.
reset_context(search_path=())
assert tuple(context.pkgs_dirs) == (str(cache),), context.pkgs_dirs
assert context.always_copy
record = dict(conda_version=conda.__version__, conda_python=sys.executable,
              configuration_search_path=[], package_cache=str(cache),
              copy_install=True, safety_checks=str(context.safety_checks),
              job=os.environ['SLURM_JOB_ID'])
(root/'installer_configuration.json').write_text(json.dumps(record, indent=2)+'\n')
for runtime, lock in [('r', 'r-linux-64.explicit.txt'), ('python', 'python-linux-64.explicit.txt')]:
    prefix = root/runtime
    assert not prefix.exists(), prefix
    _, _, status = run_command(Commands.CREATE, '--copy', '--prefix', str(prefix), '--file', str(here/lock),
                              search_path=(), stdout=None, stderr=None, use_exception_handler=False)
    assert status == 0
    assert tuple(context.pkgs_dirs) == (str(cache),)
print('ISOLATED_PREFIXES_INSTALLED', root, flush=True)
