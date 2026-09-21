"""Shared release/input checks for bounded Reviewer B research adapters."""
from pathlib import Path
import ast
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import uuid

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
sys.path.insert(0, str(PARENT))
import core_task
import geometry_control as geometry

need, sha, load, checked, write_new = geometry.need, geometry.sha, geometry.load, geometry.checked, geometry.write_new
LIBRARIES = ('CM2_glioma_other', 'CM2_primary_all_context')
SAMPLES = ('TKU4163', 'NL022', 'SN040')
BUDGETS = ('hvg2000', 'hvg5000')
OLD = HERE.parents[2] / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/controls'


def tasks():
    result = []
    for sample in SAMPLES:
        for budget in BUDGETS:
            configs = [(space, 'SNN', k, 30, 'SNN_k') for space in ('PCA30', 'UMAP2') for k in (10, 20, 40)]
            configs += [('UMAP2', method, 20, n, 'UMAP_neighbors')
                        for n in (15, 30, 60) for method in ('SNN', 'HDBSCAN_R')
                        if not (n == 30 and method == 'SNN')]
            for space, method, k, n, kind in configs:
                result.append(dict(task='neighbors', sample=sample, budget=budget, space=space,
                                   method=method, snn_k=k, umap_neighbors=n, minPts=50, resolution=.5,
                                   embedding_seed=42, kind=kind, name=f'{space}_{method}_k{k}_n{n}'))
            for seed in (0, 1, 42):
                result.append(dict(task='learning', sample=sample, budget=budget, model_seed=seed, epochs=30))
    need(len(result) == 84, 'Reviewer B roster changed')
    return result


def validate_sources():
    derivation = load(HERE / 'SOURCE_DERIVATION.json')
    for name, record in derivation['files'].items():
        need(sha(record['source']) == record['source_sha256'] and sha(HERE / name) == record['adapter_sha256'],
             f'Frozen Reviewer B source changed: {name}')
    need(sha(derivation['protocol']['path']) == derivation['protocol']['sha256'], 'Reviewer B protocol changed')
    original = Path(derivation['files']['refine_learning.py']['source']).read_text()
    actual = (HERE / 'refine_learning.py').read_text()
    def components(text):
        return [ast.dump(n, include_attributes=False) for n in ast.parse(text).body
                if (isinstance(n, ast.FunctionDef) and n.name in {'build_model', 'train_cache'})
                or (isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == 'PARAMS' for t in n.targets))]
    need(components(original) == components(actual), 'Learning arithmetic or initialization changed')
    for name, start in [('prepare_neighbors.R', "if(cfg$embedding_seed!=42L)"),
                        ('score_neighbors.R', 'deg_workers<-')]:
        old = Path(derivation['files'][name]['source']).read_text()
        new = (HERE / name).read_text()
        need(old[old.index(start):] == new[new.index(start):], f'Frozen R numerical body changed: {name}')
    return derivation


def context(gate_path):
    need(os.environ.get('SLURM_JOB_ID'), 'Reviewer B computation requires SLURM')
    need(sys.flags.no_user_site and not sys.flags.optimize and not os.environ.get('PYTHONPATH'),
         'Use frozen Python -s, without -O or PYTHONPATH')
    need(int(os.environ.get('SLURM_CPUS_PER_TASK', 0)) >= 4, 'Allocate four CPUs for exact Torch parity')
    validate_sources()
    gate, _, package = core_task.validate_gate(gate_path)
    from dgscrna.reference.runner import doctor
    need(doctor(gate['runtime']['rscript'], gate['runtime'].get('reference_r_lib')) == gate['runtime']['fingerprint'],
         'Runtime differs from released pilots')
    gate['_gate_sha256'] = sha(gate_path)
    return gate, package


def source(gate, cfg, pilot_root=None):
    if pilot_root:
        need(cfg['sample'] == 'TKU4163' and cfg['budget'] == 'hvg2000', 'Pilot is restricted to TKU4163/HVG2000')
        root = Path(pilot_root).resolve()
    else:
        root = Path(gate['output_root']) / 'core' / cfg['sample'] / cfg['budget']
    prep, provenance = geometry.source_run(root, gate, cfg['sample'], cfg['budget'], pilot=bool(pilot_root))
    return root, prep, provenance


def environment(gate, destination):
    env = os.environ.copy()
    for name in ['PYTHONPATH', 'DGSCRNA_DL_CACHE_ROOT', 'DGSCRNA_REFERENCE_OUT', 'DGSCRNA_EXAMPLE_OUT',
                 'R_LIBS', 'R_LIBS_USER', 'R_LIBS_SITE', 'DGSCRNA_REFERENCE_R_LIB']:
        env.pop(name, None)
    env.update(DGSCRNA_REFERENCE_OUT=str(destination), DGSCRNA_REQUIRE_SLURM='1',
               DGSCRNA_DL_CACHE_ROOT=str(destination / 'DL_cache' / uuid.uuid4().hex),
               DGSCRNA_DEG_WORKERS='2', PYTHONNOUSERSITE='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
               MKL_NUM_THREADS='1', R_LIBS_USER='', R_LIBS_SITE='', R_ENVIRON_USER=os.devnull, R_PROFILE_USER=os.devnull)
    if gate['runtime'].get('reference_r_lib'):
        env['DGSCRNA_REFERENCE_R_LIB'] = gate['runtime']['reference_r_lib']
    return env


def run(command, log, env):
    with log.open('x') as stream:
        result = subprocess.run(command, env=env, cwd='/tmp', stdout=stream, stderr=subprocess.STDOUT)
    need(result.returncode == 0, f'Worker failed ({result.returncode}); inspect {log}')


def evaluate(gate, dest, cfg, family, conditions, proof_path, env):
    helper = PARENT / 'evaluate_extension_lfine.py'
    need(sha(helper) == '424fb152cee5d95990b0eae724761601199230cb55512e22bd8de5cfcdfd25a9', 'Shared Lfine evaluator changed')
    spec = dict(sample=cfg['sample'], budget=cfg['budget'], family=family,
                configuration=cfg['name'], expected_conditions=len(conditions), conditions=conditions,
                parity_receipt=dict(path=str(proof_path), sha256=sha(proof_path)),
                source_hashes={str(p): sha(p) for p in [helper, HERE / 'support.py', HERE / 'SOURCE_DERIVATION.json']})
    write_new(dest / 'evaluation_spec.json', spec)
    run([gate['runtime']['python'], '-s', str(helper), '--spec', str(dest / 'evaluation_spec.json'),
         '--out', str(dest / 'evaluation')], dest / 'Lfine_evaluation.log', env)
    manifest = checked(dest / 'evaluation', 'manifest.json', 'COMPLETE')
    need(manifest['status'] == 'completed' and manifest['n_conditions'] == len(conditions)
         and manifest['n_threshold_rows'] == manifest['n_valid'] == 2 * len(conditions), 'Incomplete Lfine output')
    need(sha(dest / 'evaluation/metrics.csv.gz') == manifest['outputs']['metrics.csv.gz'], 'Lfine payload changed')
    return dict(manifest=str(dest / 'evaluation/manifest.json'), manifest_sha256=sha(dest / 'evaluation/manifest.json'),
                metrics_sha256=manifest['outputs']['metrics.csv.gz'], valid_threshold_rows=manifest['n_valid'])
