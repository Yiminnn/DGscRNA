#!/usr/bin/env python3
"""A1 research adapter: installed R-reference preparation + installed terminal DL.

The reducer/partition numerical functions and adaptive ICA are the frozen A1
sources. Only IO/environment wiring changes. This is not the package's default
PCA30-to-UMAP workflow. It never reads author truth or historical model results.
"""
from pathlib import Path
import argparse
import ast
from datetime import datetime, timezone
import fcntl
import hashlib
import importlib.metadata
import json
import os
import subprocess
import sys
import time
from types import SimpleNamespace
import warnings

CODE = Path(__file__).resolve().parent
ROOT = CODE.parents[2]
CAMPAIGN = CODE  # Frozen protocol root used by the unchanged adaptive helper.
SPACES = ['noDR', 'PCA2', 'FA2', 'ICA2', 'Isomap2', 'UMAP2', 'TSNE2']
KS = [5, 10, 15, 20, 30, 40]


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + f'.part.{os.getpid()}')
    tmp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')
    tmp.replace(path)


def checked(dest, manifest='manifest.json', flag='COMPLETE'):
    directory = Path(dest)
    return (directory / flag).is_file() and (directory / manifest).is_file() and (
        directory / flag).read_text().strip() == sha(directory / manifest)


def complete(dest, manifest='manifest.json', flag='COMPLETE'):
    (Path(dest) / flag).write_text(sha(Path(dest) / manifest) + '\n')


def log(message, **kwargs):
    print(json.dumps(dict(time=utc(), event=message, **kwargs)), flush=True)


def call_r(name, *args):
    env = os.environ.copy()
    for key in ['R_LIBS', 'PYTHONPATH', 'DGSCRNA_REFERENCE_R_LIB']:
        env.pop(key, None)
    env.update(R_LIBS_USER='', R_LIBS_SITE='', R_ENVIRON_USER='/dev/null', R_PROFILE_USER='/dev/null')
    subprocess.run([RSCRIPT, '--vanilla', str(CODE / name), *map(str, args)], check=True, env=env)


def verify_source():
    manifest = json.loads((CODE / 'SOURCE_PROTOCOL_MANIFEST.json').read_text())
    for relative, digest in manifest['original_sources'].items():
        assert sha(ROOT / relative) == digest, f'Original source changed: {relative}'
    for name, digest in manifest['adapted_files'].items():
        assert sha(CODE / name) == digest, f'Adapted source changed: {name}'
    for name, digest in manifest['frozen_protocols'].items():
        assert sha(CODE / 'protocol' / name) == digest, f'Protocol changed: {name}'
    for name, changes in manifest['R_substitutions'].items():
        restored = (CODE / name).read_text()
        for change in reversed(changes):
            assert restored.count(change['after']) == 1
            restored = restored.replace(change['after'], change['before'])
        assert restored == (CODE / 'original' / name).read_text()
    original_path = CODE / 'original' / 'run.py'
    assert sha(original_path) == manifest['original_sources'][str((ROOT /
        'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/source_snapshots/embedding_v8/run.py').relative_to(ROOT))]
    return manifest


def strict_representation_cache(core, cfg, saved):
    assert saved['source_sha256'] == sha(__file__)
    assert saved['sample'] == cfg['sample'] and saved['budget'] == cfg['budget'] and saved['space'] == cfg['space']
    if cfg['space'] == 'ICA2':
        assert saved['policy_sha256'] == sha(CODE / 'protocol/embedding_convergence_repair_20260920_v2.json')
        assert saved['adaptive_helper_sha256'] == sha(CODE / 'adaptive_ica.py')


cache_compatibility = SimpleNamespace(representation=strict_representation_cache)


def load_frozen_functions(manifest):
    source = CODE / 'original/run.py'
    tree = ast.parse(source.read_text())
    names = set(manifest['preserved_python_functions'])
    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    assert {node.name for node in functions} == names and len(functions) == len(names)
    # No imports or original root/output configuration are executed.
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(source), 'exec'), globals())


def input_contract(package_run):
    assert checked(package_run, 'run_manifest.json', 'COMPLETE'), 'Packaged input run must be complete'
    manifest = json.loads((package_run / 'run_manifest.json').read_text())
    config = json.loads((package_run / 'run_config.json').read_text())
    cfg = config['config']
    assert manifest['status'] == 'completed' and manifest['reference_labels_used_for_fit'] is False
    assert cfg['preset'] == 'gbm-reference' and cfg['features'] in {'2000', '5000'}
    sample, budget = cfg['sample'], 'hvg' + cfg['features']
    specification = json.loads((CODE / 'protocol/embedding.json').read_text())
    assert sample in specification['samples'] and budget in specification['budgets']
    assert config['input_sha256'][cfg['markers']] == sha(CODE / 'protocol/libraries.json')
    proof = package_run / 'packaged_parity.json'
    if not proof.exists():
        proof = package_run / 'parity.json'
    assert proof.exists(), 'A separately accepted packaged input is required before A1 fitting'
    accepted = json.loads(proof.read_text())
    assert accepted['status'] == 'passed_exact'
    assert accepted['run_manifest_sha256'] == sha(package_run / 'run_manifest.json')
    assert accepted['run_config_sha256'] == sha(package_run / 'run_config.json')
    prep = package_run / 'GBM' / sample / budget
    assert checked(prep, 'prepare_manifest.json', 'PREPARED')
    pm = json.loads((prep / 'prepare_manifest.json').read_text())
    assert sha(prep / 'expression_PCA30.rds') == pm['expression_sha256']
    assert Path(pm['DL_binary']).resolve().is_relative_to(package_run)
    assert sha(pm['DL_binary']) == pm['DL_binary_sha256']
    return dict(sample=sample, budget=budget, prep=str(prep),
        package_run=str(package_run), package_version=manifest['package_version'],
        package_run_manifest_sha256=sha(package_run / 'run_manifest.json'),
        package_config_sha256=sha(package_run / 'run_config.json'),
        package_parity_sha256=sha(proof), prepare_manifest_sha256=sha(prep / 'prepare_manifest.json'),
        DL_binary_sha256=pm['DL_binary_sha256'], installed_reference_sources=config['sources'])


def installed_contract(expected):
    import dgscrna.reference as package
    root = Path(package.__file__).parent
    assert root.is_relative_to(Path(sys.prefix)), 'Use the dedicated installed-wheel runtime'
    for relative, digest in expected.items():
        assert sha(root / relative) == digest, f'Installed package differs from accepted input: {relative}'
    versions = {name: importlib.metadata.version(name) for name in
                ['dgscrna', 'numpy', 'scipy', 'pandas', 'scikit-learn', 'torch', 'joblib', 'threadpoolctl']}
    assert versions['scikit-learn'] == '1.9.0'
    return dict(python=sys.executable, package_root=str(root), package_source_hashes=expected,
        package_versions=versions, BLAS_environment={name: os.environ.get(name) for name in
        ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS']},
        torch_threads=min(4, int(os.environ['SLURM_CPUS_PER_TASK'])), Rscript=RSCRIPT)


def configuration(inputs, runtime, space):
    dest = OUTPUT / inputs['sample'] / inputs['budget'] / space
    cfg = dict(sample=inputs['sample'], budget=inputs['budget'], space=space,
        prep=inputs['prep'], dest=str(dest), geometry=str(dest.parent / 'geometry'),
        embedding=str(dest / 'embedding.csv'), hdbscan_dest=str(dest / 'HDBSCAN_R'),
        frozen_protocol_sha256=sha(CODE / 'protocol/embedding.json'),
        adaptive_ICA_policy_sha256=sha(CODE / 'protocol/embedding_convergence_repair_20260920_v2.json'),
        markers=str(CODE / 'protocol/libraries.json'), source_sha256=sha(__file__),
        adapter_source_manifest_sha256=sha(CODE / 'SOURCE_PROTOCOL_MANIFEST.json'),
        reference_labels_used_for_fit=False, input_contract=inputs, runtime=runtime)
    dest.mkdir(parents=True, exist_ok=True)
    if (dest / 'config.json').exists():
        old = json.loads((dest / 'config.json').read_text())
        assert {k: v for k, v in old.items() if k != 'conditions'} == cfg, 'Resume contract changed'
        return old
    write_json(dest / 'config.json', cfg)
    return cfg


def fit_space(cfg):
    from dgscrna.reference.backend import terminal
    dest = Path(cfg['dest'])
    lock = dest / '.fit.lock'
    lock.mkdir()  # No automatic reclaim; a stale lock requires checking its scheduler job.
    write_json(lock / 'owner.json', dict(job=os.environ['SLURM_JOB_ID'], pid=os.getpid(), host=os.uname().nodename))
    try:
        if checked(dest, 'fit_manifest.json', 'FIT_COMPLETE'):
            manifest = json.loads((dest / 'fit_manifest.json').read_text())
            assert manifest['source_sha256'] == sha(__file__) and manifest['config'] == cfg
            for condition in cfg['conditions']:
                path = Path(condition['dest']) / 'terminal/L00_mean'
                assert checked(path, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
                tm = json.loads((path / 'terminal_manifest.json').read_text())
                assert sha(path / 'predictions.csv.gz') == tm['predictions_sha256']
                assert sha(path / 'terminal_manifest.json') == manifest['terminal_manifest_hashes'][condition['route']]
            return
        x, cells, geometry = get_geometry(cfg)
        z = make_embedding(cfg, x, cells, geometry)
        cfg['conditions'] = make_partitions(cfg, z, cells)
        assert len(cfg['conditions']) == 13
        write_json(dest / 'config.json', cfg)
        call_r('score_candidates.R', dest / 'config.json')
        os.environ['DGSCRNA_DL_CACHE_ROOT'] = str(dest / 'DL_cache')
        for condition in cfg['conditions']:
            assert sha(Path(condition['dest']) / 'cells.csv') == sha(Path(cfg['prep']) / 'cells.csv'), 'DL cell order changed'
            terminal.finish_route(condition['dest'], 'L00_mean')
        hashes = {}
        for condition in cfg['conditions']:
            path = Path(condition['dest']) / 'terminal/L00_mean'
            assert checked(path, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
            tm = json.loads((path / 'terminal_manifest.json').read_text())
            assert Path(tm['cache_directory']).resolve().is_relative_to(dest / 'DL_cache')
            assert tm['reference_labels_used_for_fit'] is False and tm['terminal_valid']
            hashes[condition['route']] = sha(path / 'terminal_manifest.json')
        write_json(dest / 'fit_manifest.json', dict(status='fit_complete_evaluation_pending', config=cfg,
            terminal_manifest_hashes=hashes, installed_package_terminal=True, old_models_or_caches_used=False,
            reference_labels_used_for_fit=False, source_sha256=sha(__file__),
            job=os.environ['SLURM_JOB_ID'], completed_at=utc()))
        complete(dest, 'fit_manifest.json', 'FIT_COMPLETE')
        log('fit_complete', sample=cfg['sample'], budget=cfg['budget'], space=cfg['space'], n_conditions=13)
    finally:
        (lock / 'owner.json').unlink()
        lock.rmdir()


def main():
    global OUTPUT, RSCRIPT, adaptive_ica
    assert os.environ.get('SLURM_JOB_ID'), 'Scientific fitting requires SLURM'
    assert not sys.flags.optimize
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package-run', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--rscript', required=True)
    parser.add_argument('--spaces', nargs='+', choices=SPACES, default=SPACES)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    OUTPUT, RSCRIPT = args.out.resolve(), str(Path(args.rscript).resolve())
    allowed = ROOT / 'results/hvg_ptc_20260916_v1/package_reference_rerun_20260920'
    assert OUTPUT.is_relative_to(allowed) and OUTPUT != allowed
    os.environ['DGSCRNA_REFERENCE_OUT'] = str(OUTPUT)
    os.environ['DGSCRNA_REQUIRE_SLURM'] = '1'
    os.environ['DGSCRNA_DEG_WORKERS'] = str(min(4, int(os.environ['SLURM_CPUS_PER_TASK'])))
    for name in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS']:
        assert os.environ.get(name) == '1', f'Set {name}=1 before starting Python'
    manifest = verify_source()
    import adaptive_ica
    assert sha(Path(adaptive_ica.__file__)) == manifest['adapted_files']['adaptive_ica.py']
    load_frozen_functions(manifest)
    inputs = input_contract(args.package_run.resolve())
    runtime = installed_contract(inputs['installed_reference_sources'])
    from dgscrna.reference.backend import terminal
    terminal.threads()
    for space in args.spaces:
        dest = OUTPUT / inputs['sample'] / inputs['budget'] / space
        if (dest / 'config.json').exists() and not args.resume:
            raise RuntimeError(f'Existing result requires explicit --resume: {dest}')
        fit_space(configuration(inputs, runtime, space))
    verify_source()


if __name__ == '__main__':
    main()
