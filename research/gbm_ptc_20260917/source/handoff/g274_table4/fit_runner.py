"""隔离执行Table4拟合并留存证据；不导入旧脚本，不读取labels.csv。

主CLI只允许SLURM；API允许至多128个细胞的合成测试。每个stage独立进程，
记录的是该进程生命周期的ru_maxrss，不冒称算法净增量；SLURM MaxRSS另行记录。
metrics故意保持not_run，等待独立测量契约，不猜测Fscore/NPE/RV/ED。
"""
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
from importlib.metadata import PackageNotFoundError, version
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import time
import traceback

from runtime_guard import THREAD_KEYS, bootstrap, worker_environment

if __name__ == '__main__':
    bootstrap()

import numpy as np

from fit_core import (CLUSTER_ORDER, DR_ORDER, FitConfig, cluster_matrix,
                      prepare_matrix, reduce_matrix)


SCHEMA_VERSION = 2
SOURCE_FILES = ('fit_core.py', 'fit_runner.py', 'runtime_guard.py')


class IntegrityError(RuntimeError):
    """冻结输入或运行依赖与预期hash不一致。"""


def verify_dependencies(expected):
    if not expected:
        raise IntegrityError('缺少expected dependency hashes')
    for name, digest in expected.items():
        path = Path(name)
        if not path.is_file() or sha256(path) != digest:
            raise IntegrityError(f'artifact哈希不匹配：{name}')


def expected_dependencies(out, manifest, dr=None):
    expected = {str(out / 'source' / name): digest
                for name, digest in manifest['source_sha256'].items()}
    expected.update({str(Path(manifest['input_dir']) / name): info['sha256']
                     for name, info in manifest['input_files'].items()})
    for name, info in manifest.get('prepare', {}).get('artifacts', {}).items():
        expected[str(out / name)] = info['sha256']
    if dr is not None:
        for name, info in manifest['dr'][dr].get('artifacts', {}).items():
            expected[str(out / name)] = info['sha256']
    return expected


def _json_value(value):
    if isinstance(value, dict):
        return {str(k): _json_value(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(v) for v in value]
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return str(value)
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    return repr(value)


def write_json(path, obj):
    path = Path(path)
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(_json_value(obj), ensure_ascii=False, indent=2,
                              allow_nan=False) + '\n')
    tmp.replace(path)


def sha256(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def _utc():
    return datetime.now(timezone.utc).isoformat()


def _local_size_guard(n_cells):
    if not os.environ.get('SLURM_JOB_ID') and n_cells > 128:
        raise RuntimeError('实际研究矩阵只能在SLURM中拟合；本地仅允许≤128细胞微型测试')


def load_clean_counts(input_dir):
    """只读无annotation的counts对象，标签对齐由后续评估独立处理。"""
    import anndata as ad

    input_dir = Path(input_dir)
    if not (input_dir / 'COMPLETE').is_file():
        raise ValueError('输入probe尚未COMPLETE')
    manifest = json.loads((input_dir / 'manifest.json').read_text())
    if (manifest.get('status') != 'completed' or manifest.get('HVG_selection') is not False
            or manifest.get('gold_used_for_fitting') is not False):
        raise ValueError('输入manifest不符合无HVG/无gold契约')
    expected = (manifest['n_cells'], manifest['n_genes_retained'])
    _local_size_guard(expected[0])
    a = ad.read_h5ad(input_dir / 'counts_gene_filtered.h5ad')
    if a.shape != expected:
        raise ValueError('h5ad形状与冻结manifest不符')
    if (len(a.obs.columns) or len(a.var.columns) or a.uns or len(a.obsm) or len(a.varm)
            or len(a.layers) or len(a.obsp) or len(a.varp) or a.raw is not None):
        raise ValueError('拟合输入含annotation或额外数据；禁止历史labels/cluster泄漏')
    if not a.obs_names.is_unique or not a.var_names.is_unique:
        raise ValueError('细胞或基因ID重复')
    if any(not str(x) or '\n' in str(x) or '\t' in str(x) for x in [*a.obs_names, *a.var_names]):
        raise ValueError('细胞或基因ID不能含空白ID/换行/tab')
    return a


def _artifact_info(path, output_dir):
    path = Path(path)
    return str(path.relative_to(output_dir)), {'sha256': sha256(path),
                                              'size_bytes': path.stat().st_size}


def worker(request_path):
    """单stage工作进程；异常也产生result，进程被杀则由父进程补failed记录。"""
    request = json.loads(Path(request_path).read_text())
    out = Path(request['output_dir'])
    dest = Path(request_path).parent
    started = time.perf_counter()
    source = Path(__file__).resolve().parent
    result = {'status': 'running', 'started_at': _utc(), 'pid': os.getpid(),
              'kind': request['kind'], 'requested_config': request['config'],
              'worker_source_sha256': {name: sha256(source / name)
                                       for name in SOURCE_FILES},
              'peak_rss_scope': 'isolated_worker_lifetime_including_load_and_libraries'}
    write_json(dest / 'result.json', result)
    code = 0
    try:
        verify_dependencies(request.get('expected_dependencies'))
        import sklearn
        from threadpoolctl import threadpool_info, threadpool_limits
        limit = int(os.environ['OMP_NUM_THREADS'])
        controller = threadpool_limits(limits=limit)
        result.update(thread_limit=limit,
                      runtime_environment={k: os.environ.get(k) for k in (*THREAD_KEYS, 'PYTHONHASHSEED')})
        config = FitConfig(**request['config'])
        artifacts = []
        if request['kind'] == 'prepare':
            a = load_clean_counts(request['input_dir'])
            t_fit = time.perf_counter()
            x = prepare_matrix(a.X, config.target_sum, config.max_value)
            result['fit_seconds'] = time.perf_counter() - t_fit
            path = dest / 'X.npy'
            np.save(path, x, allow_pickle=False)
            (dest / 'cell_ids.tsv').write_text('\n'.join(a.obs_names) + '\n')
            (dest / 'genes.tsv').write_text('\n'.join(a.var_names) + '\n')
            result.update(n_cells=x.shape[0], n_genes=x.shape[1],
                          matrix_file=str(path.relative_to(out)),
                          preprocessing='normalize_total_log1p_scale',
                          gene_selection='all_input_genes; no_HVG; no_additional_filter',
                          params={'target_sum': config.target_sum, 'max_value': config.max_value})
            artifacts = [path, dest / 'cell_ids.tsv', dest / 'genes.tsv']
        elif request['kind'] == 'dr':
            source = out / 'prepare' / 'X.npy'
            x = np.load(source, mmap_mode='r', allow_pickle=False)
            _local_size_guard(x.shape[0])
            t_fit = time.perf_counter()
            z, details = reduce_matrix(x, request['dr'], config)
            result['fit_seconds'] = time.perf_counter() - t_fit
            result.update(details)
            path = dest / 'Z.npy'
            if request['dr'] == 'none':
                # 同一文件硬链接，保留真实全维Z且不重复复制约456MB矩阵。
                os.link(source, path)
            else:
                np.save(path, z, allow_pickle=False)
            result.update(n_cells=z.shape[0], actual_dimension=z.shape[1],
                          coordinates_file=str(path.relative_to(out)),
                          cell_ids_file='prepare/cell_ids.tsv',
                          cell_ids_sha256=sha256(out / 'prepare' / 'cell_ids.tsv'),
                          fit_input_sha256=sha256(source))
            artifacts = [path]
        elif request['kind'] == 'cluster':
            source = out / 'dr' / request['dr'] / 'Z.npy'
            z = np.load(source, mmap_mode='r', allow_pickle=False)
            _local_size_guard(z.shape[0])
            t_fit = time.perf_counter()
            labels, details = cluster_matrix(z, request['clusterer'], request['dr'], config)
            result['fit_seconds'] = time.perf_counter() - t_fit
            result.update(details)
            path = dest / 'labels.npy'
            np.save(path, labels, allow_pickle=False)
            result.update(labels_file=str(path.relative_to(out)),
                          actual_fit_dimension=z.shape[1], fit_input_sha256=sha256(source),
                          cell_ids_file='prepare/cell_ids.tsv',
                          cell_ids_sha256=sha256(out / 'prepare' / 'cell_ids.tsv'))
            artifacts = [path]
        else:
            raise ValueError(f"未知stage：{request['kind']}")
        verify_dependencies(request.get('expected_dependencies'))
        result['dependency_verification'] = 'before_and_after'
        result['consumed_dependencies'] = request['expected_dependencies']
        result['threadpools'] = threadpool_info()
        if any(pool['num_threads'] > limit for pool in result['threadpools']):
            raise RuntimeError('实际线程池超过分配上限')
        result['artifacts'] = dict(_artifact_info(path, out) for path in artifacts)
        result['status'] = 'completed'
    except BaseException as exc:
        code = 1
        result.update(status='failed', error_type=type(exc).__name__, error=str(exc),
                      traceback=traceback.format_exc())
    finally:
        result.update(ended_at=_utc(), elapsed_seconds=time.perf_counter() - started,
                      peak_rss_bytes=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) * 1024)
        write_json(dest / 'result.json', result)
    return code


def _seal_stage_metadata(result, dest, out):
    # 只把JSON自身hash加入父manifest，不再写回result.json，避免self-hash循环。
    artifacts = result.setdefault('artifacts', {})
    for name in ('result.json', 'request.json', 'worker.log'):
        path = dest / name
        if path.is_file():
            relative, info = _artifact_info(path, out)
            artifacts[relative] = info
    return result


def _dispatch(out, relative_dir, request, stage_timeout_s):
    dest = out / relative_dir
    dest.mkdir(parents=True, exist_ok=False)
    request.update(output_dir=str(out))
    write_json(dest / 'request.json', request)
    t0 = time.perf_counter()
    try:
        verify_dependencies(request.get('expected_dependencies'))
    except IntegrityError as exc:
        result = {'status': 'failed', 'error_type': 'IntegrityError', 'error': str(exc),
                  'fit_seconds': None, 'peak_rss_bytes': None,
                  'elapsed_seconds': time.perf_counter() - t0}
        write_json(dest / 'result.json', result)
        return _seal_stage_metadata(result, dest, out)
    with (dest / 'worker.log').open('w') as log:
        try:
            proc = subprocess.run([sys.executable, str(out / 'source' / 'fit_runner.py'),
                                   '--worker', str(dest / 'request.json')],
                                  stdout=log, stderr=subprocess.STDOUT, timeout=stage_timeout_s,
                                  env=worker_environment(os.environ), check=False)
            failure = {'error_type': 'WorkerExit', 'error': f'worker returncode={proc.returncode}',
                       'returncode': proc.returncode}
        except subprocess.TimeoutExpired:
            proc = None
            failure = {'error_type': 'StageTimeout', 'error': f'stage超过{stage_timeout_s}秒',
                       'returncode': None}
    path = dest / 'result.json'
    result = json.loads(path.read_text()) if path.exists() else {}
    if proc is None or proc.returncode != 0 or result.get('status') != 'completed':
        if result.get('status') != 'failed':
            result.update(status='failed', **failure)
        result.setdefault('peak_rss_bytes', None)
        result.setdefault('fit_seconds', None)
        result.setdefault('elapsed_seconds', time.perf_counter() - t0)
        write_json(path, result)
    result['worker_wall_seconds'] = time.perf_counter() - t0
    result['log_file'] = str((dest / 'worker.log').relative_to(out))
    return _seal_stage_metadata(result, dest, out)


def _versions():
    result = {'python': sys.version}
    for package in ('numpy', 'scipy', 'scikit-learn', 'anndata', 'scanpy', 'umap-learn', 'hdbscan'):
        try:
            result[package] = version(package)
        except PackageNotFoundError:
            result[package] = 'not_installed'
    return result


def run_fit(input_dir, output_dir, config, dr_methods=DR_ORDER,
            clusterers=CLUSTER_ORDER, stage_timeout_s=None):
    """产物只写新目录；每格独立状态，未选择的格明确not_run，不伪装零噪声。"""
    input_dir, out = Path(input_dir).resolve(), Path(output_dir).resolve()
    if not dr_methods or not clusterers or len(set(dr_methods)) != len(dr_methods) or len(set(clusterers)) != len(clusterers):
        raise ValueError('至少指定一个无重复的DR和clusterer')
    if set(dr_methods) - set(DR_ORDER) or set(clusterers) - set(CLUSTER_ORDER):
        raise ValueError('DR或clusterer不在冻结21格名单')
    contract = json.loads((input_dir / 'manifest.json').read_text())
    _local_size_guard(contract['n_cells'])
    out.mkdir(parents=True, exist_ok=False)
    original_source = Path(__file__).resolve().parent
    source = out / 'source'
    source.mkdir()
    for name in SOURCE_FILES:
        shutil.copy2(original_source / name, source / name)
    manifest = {
        'schema_version': SCHEMA_VERSION, 'status': 'running', 'started_at': _utc(),
        'input_dir': str(input_dir), 'config': asdict(config), 'gold_used_for_fitting': False,
        'selected_dr': list(dr_methods), 'selected_clusterers': list(clusterers),
        'expected_dr': list(DR_ORDER), 'expected_clusterers': list(CLUSTER_ORDER),
        'stage_timeout_seconds': stage_timeout_s, 'versions': _versions(),
        'source_sha256': {name: sha256(source / name) for name in SOURCE_FILES},
        'input_files': {name: {'sha256': sha256(input_dir / name),
                              'size_bytes': (input_dir / name).stat().st_size}
                        for name in ('counts_gene_filtered.h5ad', 'manifest.json', 'COMPLETE')},
        'slurm_job_id': os.environ.get('SLURM_JOB_ID'),
        'slurm_maxrss_bytes': None,
        'threads': {k: os.environ.get(k) for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                                                 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS', 'PYTHONHASHSEED')},
        'metrics': {'status': 'not_run', 'reason': 'pending_metric_contract'},
        'dr': {d: {'status': 'not_run', 'coordinates_file': None,
                   'actual_dimension': None} for d in DR_ORDER},
        'arms': {f'{d}__{c}': {'status': 'not_run', 'dr': d, 'clusterer': c,
                               'labels_file': None, 'n_noise': None, 'n_clusters': None,
                               'noise_frac': None} for d in DR_ORDER for c in CLUSTER_ORDER}}
    write_json(out / 'manifest.json', manifest)
    request = {'kind': 'prepare', 'input_dir': str(input_dir), 'config': asdict(config),
               'expected_dependencies': expected_dependencies(out, manifest)}
    manifest['prepare'] = _dispatch(out, 'prepare', request, stage_timeout_s)
    write_json(out / 'manifest.json', manifest)
    if manifest['prepare']['status'] == 'completed':
        for dr in dr_methods:
            manifest['dr'][dr]['status'] = 'running'
            write_json(out / 'manifest.json', manifest)
            request = {'kind': 'dr', 'dr': dr, 'config': asdict(config),
                       'expected_dependencies': expected_dependencies(out, manifest)}
            state = _dispatch(out, f'dr/{dr}', request, stage_timeout_s)
            manifest['dr'][dr].update(state)
            if state['status'] != 'completed':
                for cl in clusterers:
                    manifest['arms'][f'{dr}__{cl}'].update(status='blocked', reason='DR_failed')
                write_json(out / 'manifest.json', manifest)
                continue
            for cl in clusterers:
                key = f'{dr}__{cl}'
                manifest['arms'][key]['status'] = 'running'
                write_json(out / 'manifest.json', manifest)
                request = {'kind': 'cluster', 'dr': dr, 'clusterer': cl, 'config': asdict(config),
                           'expected_dependencies': expected_dependencies(out, manifest, dr)}
                state = _dispatch(out, f'arms/{key}', request, stage_timeout_s)
                manifest['arms'][key].update(state)
                write_json(out / 'manifest.json', manifest)
    else:
        for dr in dr_methods:
            manifest['dr'][dr].update(status='blocked', reason='prepare_failed')
            for cl in clusterers:
                manifest['arms'][f'{dr}__{cl}'].update(status='blocked', reason='prepare_failed')
    canvas = ('dr/PCA/Z.npy' if manifest['dr']['PCA']['status'] == 'completed' else None)
    manifest['common_display_coordinates'] = canvas
    manifest['dr']['none'].update(display_only=True, display_coordinates=canvas,
                                   display_note='仅共同PCA2画布；实际聚类在完整scaled表达空间')
    states = [s['status'] for s in manifest['arms'].values()]
    manifest['status'] = ('completed' if all(s == 'completed' for s in states)
                          else 'failed' if manifest['prepare']['status'] != 'completed' else 'partial')
    final_dependencies = expected_dependencies(out, manifest)
    for state in [*manifest['dr'].values(), *manifest['arms'].values()]:
        for name, info in state.get('artifacts', {}).items():
            final_dependencies[str(out / name)] = info['sha256']
    try:
        verify_dependencies(final_dependencies)
        manifest['final_integrity'] = {'status': 'verified', 'expected_dependencies': final_dependencies}
    except IntegrityError as exc:
        manifest['status'] = 'failed'
        manifest['final_integrity'] = {'status': 'failed', 'error': str(exc)}
    manifest['ended_at'] = _utc()
    write_json(out / 'manifest.json', manifest)
    if manifest['status'] == 'completed':
        (out / 'FIT_COMPLETE').write_text('21 arms fitted; evaluation remains independent\n')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', help=argparse.SUPPRESS)
    parser.add_argument('--input-dir', type=Path)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--gmm-covariance', choices=('full', 'diag', 'full_2d_diag_none'))
    parser.add_argument('--dr', nargs='+', choices=DR_ORDER, default=list(DR_ORDER))
    parser.add_argument('--clusterers', nargs='+', choices=CLUSTER_ORDER, default=list(CLUSTER_ORDER))
    parser.add_argument('--stage-timeout-seconds', type=float)
    args = parser.parse_args()
    if args.worker:
        return worker(args.worker)
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('正式拟合必须在SLURM作业内运行')
    if args.input_dir is None or args.output_dir is None or args.gmm_covariance is None:
        parser.error('--input-dir、--output-dir、--gmm-covariance必须显式提供')
    manifest = run_fit(args.input_dir, args.output_dir, FitConfig(gmm_covariance=args.gmm_covariance),
                       args.dr, args.clusterers, args.stage_timeout_seconds)
    print(json.dumps({'status': manifest['status'], 'output_dir': str(args.output_dir)}, ensure_ascii=False))
    return 0 if (manifest.get('final_integrity', {}).get('status') == 'verified'
                 and all(manifest['arms'][f'{d}__{c}']['status'] == 'completed'
                         for d in args.dr for c in args.clusterers)) else 1


if __name__ == '__main__':
    raise SystemExit(main())
