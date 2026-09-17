"""纯stdlib启动保护；必须先于numpy/sklearn/numba导入。"""
import os
from pathlib import Path
import sys

THREAD_KEYS = ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMBA_NUM_THREADS')


def worker_environment(environ):
    env = dict(environ)
    allocation = int(env.get('SLURM_CPUS_PER_TASK', '1'))
    requested = int(env.get('G274_THREADS', str(allocation)))
    if allocation < 1 or requested < 1:
        raise ValueError('线程数与SLURM分配CPU数必须为正整数')
    limit = min(allocation, requested)
    for key in THREAD_KEYS:
        env[key] = str(limit)
    env.update(PYTHONHASHSEED='42', OMP_DYNAMIC='FALSE', MKL_DYNAMIC='FALSE')
    return env


def bootstrap():
    """重新exec以确保PYTHONHASHSEED对解释器启动生效，而非只改运行时字典。"""
    env = worker_environment(os.environ)
    keys = (*THREAD_KEYS, 'PYTHONHASHSEED', 'OMP_DYNAMIC', 'MKL_DYNAMIC')
    if any(os.environ.get(key) != env[key] for key in keys):
        os.execve(sys.executable, [sys.executable, str(Path(sys.argv[0]).resolve()), *sys.argv[1:]], env)
