"""A1 fitting: identical R scaled-HVG, explicit reducers, native R scorer and DL.

No author labels are read anywhere in this fitting entry point. Evaluation and
patient-fold K selection occur in separate code after candidates are frozen.
"""
from pathlib import Path
import argparse
import fcntl
import json
import os
import shutil
import subprocess
import sys
import time
import warnings

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE))
sys.path.insert(0, str(CODE / 'legacy'))
import adaptive_ica
import cache_compatibility
from common import ROOT, OUT as REFERENCE, RSCRIPT, require_slurm, sha, utc, write_json, checked, complete
CAMPAIGN = ROOT / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
OUTPUT = CAMPAIGN / 'embedding'
SPACES = ['noDR', 'PCA2', 'FA2', 'ICA2', 'Isomap2', 'UMAP2', 'TSNE2']
KS = [5, 10, 15, 20, 30, 40]


def log(message, **kwargs):
    print(json.dumps(dict(time=utc(), event=message, **kwargs)), flush=True)


def call_r(name, *args):
    subprocess.run([RSCRIPT, str(CODE / name), *map(str, args)], check=True)


def config(sample, budget, space):
    assert sample in (ROOT / 'handoff/g274/cohort_samples.txt').read_text().split()
    assert budget in ['hvg2000', 'hvg5000'] and space in SPACES + ['anchor_parity']
    dest = OUTPUT / sample / budget / space
    dest.mkdir(parents=True, exist_ok=True)
    result = dict(sample=sample, budget=budget, space=space,
                  prep=str(REFERENCE / 'GBM' / sample / budget), dest=str(dest),
                  geometry=str(OUTPUT / sample / budget / 'geometry'),
                  embedding=str(dest / 'embedding.csv'),
                  hdbscan_dest=str(dest / 'HDBSCAN_R'),
                  frozen_protocol_sha256=sha(CAMPAIGN / 'protocol/embedding.json'),
                  source_sha256=sha(__file__), reference_labels_used_for_fit=False)
    if (dest / 'config.json').exists():
        previous = json.loads((dest / 'config.json').read_text())
        for name in ['sample', 'budget', 'space', 'prep', 'dest', 'geometry', 'frozen_protocol_sha256']:
            assert previous[name] == result[name], f'Changed immutable condition: {name}'
        cache_compatibility.record_consumer(sys.modules[__name__],previous,'configuration',
            dest/'config.json',previous['source_sha256'])
        # Preserve the configuration's producing source; consumer provenance is
        # separate. Completed v5/v6 caches are not rewritten as current-source products.
        return previous
    result.update(configuration_origin_source_sha256=sha(__file__),
        adaptive_ICA_policy_sha256=sha(CAMPAIGN/'protocol/embedding_convergence_repair_20260920_v2.json'))
    write_json(dest / 'config.json', result)
    return result


def get_geometry(cfg):
    import numpy as np
    import pandas as pd
    directory = Path(cfg['geometry'])
    directory.mkdir(exist_ok=True)
    # Exactly one exporter per sample/budget; science jobs are also submitted
    # per sample/budget for initial pilot. Cross-node duplicates must compare.
    with (directory / 'export.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not checked(directory):
            call_r('export_geometry.R', Path(cfg['dest']) / 'config.json')
    assert checked(directory)
    m = json.loads((directory / 'manifest.json').read_text())
    assert m['prepare_manifest_sha256'] == sha(Path(cfg['prep']) / 'prepare_manifest.json')
    assert sha(m['binary']) == m['binary_sha256']
    assert sha(directory / 'cells.csv') == m['cells_sha256']
    x = np.memmap(m['binary'], mode='r', dtype='<f8', shape=(m['n_cells'], m['n_features']))
    assert np.isfinite(x).all()
    cells = pd.read_csv(directory / 'cells.csv', dtype=str, keep_default_na=False).cell_id.to_numpy()
    return x, cells, m


def model_metadata(model, caught, seconds):
    import numpy as np
    output = dict(implementation=type(model).__module__ + '.' + type(model).__name__,
                  params=model.get_params(deep=False), elapsed_seconds=seconds,
                  warnings=[f'{w.category.__name__}: {w.message}' for w in caught])
    for name in ['n_iter_', 'converged_', 'kl_divergence_']:
        value = getattr(model, name, None)
        if value is not None:
            output[name] = value.item() if isinstance(value, np.generic) else value
    return output


def fit_model(model, x, predict=False):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        t0 = time.monotonic()
        z = model.fit_predict(x) if predict else model.fit_transform(x)
    return z, model_metadata(model, caught, time.monotonic() - t0)


def make_embedding(cfg, x, cells, geometry_manifest):
    import numpy as np
    import pandas as pd
    from sklearn.decomposition import PCA, FactorAnalysis, FastICA
    from sklearn.manifold import Isomap, TSNE
    import sklearn
    dest = Path(cfg['dest']);space = cfg['space']
    if checked(dest, 'representation.json', 'REPRESENTATION_COMPLETE'):
        m = json.loads((dest / 'representation.json').read_text())
        cache_compatibility.representation(sys.modules[__name__],cfg,m)
        assert m['geometry_binary_sha256'] == geometry_manifest['binary_sha256']
        assert m['protocol_sha256'] == cfg['frozen_protocol_sha256']
        if space == 'noDR':
            return x
        assert sha(cfg['embedding']) == m['embedding_sha256']
        frame = pd.read_csv(cfg['embedding'], keep_default_na=False)
        assert np.array_equal(frame.cell_id.astype(str), cells)
        return frame[['x', 'y']].to_numpy()
    if space == 'noDR':
        z = x;metadata = dict(implementation='identity', params={})
    elif space == 'UMAP2':
        t0 = time.monotonic()
        call_r('native_geometry.R', 'umap', dest / 'config.json')
        frame = pd.read_csv(cfg['embedding'], keep_default_na=False)
        assert np.array_equal(frame.cell_id.astype(str), cells)
        z = frame[['x', 'y']].to_numpy()
        metadata = json.loads(Path(cfg['embedding'] + '.params.json').read_text())
        metadata['elapsed_seconds'] = time.monotonic() - t0
    elif space == 'ICA2':
        z,metadata = adaptive_ica.fit(sys.modules[__name__],cfg,x)
    else:
        if space == 'PCA2':
            model = PCA(n_components=2, svd_solver='randomized', iterated_power=7, random_state=42)
        elif space == 'FA2':
            model = FactorAnalysis(n_components=2, max_iter=5000, tol=0.01, random_state=42)
        elif space == 'Isomap2':
            model = Isomap(n_components=2, n_neighbors=15, n_jobs=1)
        elif space == 'TSNE2':
            model = TSNE(n_components=2, perplexity=30, init='pca', learning_rate='auto',
                         max_iter=1000, random_state=42, n_jobs=1)
        else:
            raise ValueError(space)
        z, metadata = fit_model(model, x)
        metadata['sklearn_version'] = sklearn.__version__
        from sklearn.exceptions import ConvergenceWarning
        if any(w.startswith('ConvergenceWarning:') for w in metadata['warnings']):
            write_json(dest / 'CONVERGENCE_FAILURE.json', metadata)
            raise RuntimeError(f'{space}: unconverged representation; requires repair, not a valid result')
    assert z.shape == (len(cells), x.shape[1] if space == 'noDR' else 2) and np.isfinite(z).all()
    if space not in ['UMAP2', 'noDR']:
        pd.DataFrame(dict(cell_id=cells, x=z[:, 0], y=z[:, 1])).to_csv(cfg['embedding'], index=False)
    metadata.update(status='completed', space=space, sample=cfg['sample'], budget=cfg['budget'],
                    geometry_binary_sha256=geometry_manifest['binary_sha256'],
                    protocol_sha256=cfg['frozen_protocol_sha256'], n_cells=len(cells),
                    embedding_sha256=None if space == 'noDR' else sha(cfg['embedding']),
                    job=os.environ['SLURM_JOB_ID'], reference_labels_used_for_fit=False,
                    source_sha256=sha(__file__), completed_at=utc())
    write_json(dest / 'representation.json', metadata)
    complete(dest, 'representation.json', 'REPRESENTATION_COMPLETE')
    return z


def make_partitions(cfg, z, cells):
    import numpy as np
    import pandas as pd
    from sklearn.cluster import KMeans
    from sklearn.mixture import GaussianMixture
    dest = Path(cfg['dest']);conditions = []
    for method in ['KMeans', 'GMM', 'HDBSCAN_R']:
        for k in KS if method != 'HDBSCAN_R' else [None]:
            name = method if k is None else f'{method}_K{k:02d}'
            directory = dest / name;directory.mkdir(exist_ok=True)
            item = dict(route=cfg['space'] + '_' + name, dest=str(directory), method=method, k=k)
            conditions.append(item)
            if checked(directory, 'partition_manifest.json', 'PARTITION_COMPLETE'):
                prior = json.loads((directory / 'partition_manifest.json').read_text())
                assert prior['protocol_sha256'] == cfg['frozen_protocol_sha256']
                assert prior['representation_manifest_sha256'] == sha(dest / 'representation.json')
                assert prior['clusters_sha256'] == sha(directory / 'clusters.csv')
                continue
            t0 = time.monotonic()
            if k is not None and len(cells) < k:
                write_json(directory / 'STRUCTURALLY_INAPPLICABLE.json', dict(reason='n_cells < K', k=k, n_cells=len(cells)))
                raise RuntimeError('Unexpected sample smaller than frozen K; needs explicit manifest accounting')
            if method == 'HDBSCAN_R':
                call_r('native_geometry.R', 'hdbscan', dest / 'config.json')
                frame = pd.read_csv(directory / 'clusters.csv', dtype=str, keep_default_na=False)
                assert np.array_equal(frame.cell_id, cells)
                metadata = json.loads((directory / 'cluster_params.json').read_text())
            else:
                if method == 'KMeans':
                    model = KMeans(n_clusters=k, n_init=10, max_iter=300, tol=1e-4, algorithm='lloyd', random_state=42)
                else:
                    model = GaussianMixture(n_components=k, covariance_type='diag', reg_covar=1e-4,
                                            max_iter=1000, tol=1e-3, n_init=1, random_state=42)
                labels, metadata = fit_model(model, z, predict=True)
                if not bool(getattr(model, 'converged_', True)):
                    write_json(directory / 'CONVERGENCE_FAILURE.json', metadata)
                    raise RuntimeError(f'{name}: unconverged GMM requires repair')
                assert labels.shape == (len(cells),)
                pd.DataFrame(dict(cell_id=cells, cluster=labels)).to_csv(directory / 'clusters.csv', index=False)
            metadata.update(item, status='completed', elapsed_total_seconds=time.monotonic() - t0,
                            protocol_sha256=cfg['frozen_protocol_sha256'],
                            representation_manifest_sha256=sha(dest / 'representation.json'),
                            clusters_sha256=sha(directory / 'clusters.csv'),
                            n_cells=len(cells), job=os.environ['SLURM_JOB_ID'], completed_at=utc(),
                            reference_labels_used_for_fit=False)
            write_json(directory / 'partition_manifest.json', metadata)
            complete(directory, 'partition_manifest.json', 'PARTITION_COMPLETE')
            log('partition_complete', sample=cfg['sample'], budget=cfg['budget'], route=item['route'])
    return conditions


def finish(cfg):
    import terminal
    terminal.OUT = CAMPAIGN
    os.environ['DGSCRNA_DL_CACHE_ROOT'] = str(OUTPUT / 'DL_cache')
    for condition in cfg['conditions']:
        terminal.finish_route(condition['dest'], 'L00_mean')


def parity(sample, budget):
    import numpy as np
    import pandas as pd
    cfg = config(sample, budget, 'anchor_parity')
    dest = Path(cfg['dest']);reference = Path(cfg['prep']) / 'UMAP2_HDBSCAN_R'
    condition = dict(route='UMAP2_HDBSCAN_R', dest=str(dest / 'UMAP2_HDBSCAN_R'), method='HDBSCAN_R', k=None)
    route = Path(condition['dest']);route.mkdir(exist_ok=True)
    shutil.copy2(reference / 'clusters.csv', route / 'clusters.csv')
    cfg['conditions'] = [condition];write_json(dest / 'config.json', cfg)
    call_r('score_candidates.R', dest / 'config.json')
    new = pd.read_csv(route / 'initial_calls.csv.gz', dtype=str, keep_default_na=False)
    old = pd.read_csv(reference / 'initial_calls.csv.gz', dtype=str, keep_default_na=False)
    assert np.array_equal(new.cell_id, old.cell_id)
    sm = json.loads((reference / 'score_manifest.json').read_text())
    old_id = next(a for a, m in sm['arms'].items() if m['library'] == 'CM2_glioma_other' and m['cutoff'] == 'mean')
    assert np.array_equal(new.L00_mean, old[old_id]), 'Native R initial labels differ'
    finish(cfg)
    with np.load(route / 'terminal/L00_mean/terminal.npz') as a, np.load(reference / f'terminal/{old_id}/terminal.npz') as b:
        for key in a.files:
            if key == 'probabilities':
                np.testing.assert_allclose(a[key], b[key], rtol=1e-6, atol=1e-7)
            elif key == 'confidence_rounded':
                np.testing.assert_allclose(a[key], b[key], rtol=0, atol=0, equal_nan=True)
            else:
                assert np.array_equal(a[key], b[key]), key
    write_json(dest / 'PARITY_PASSED.json', dict(status='passed', sample=sample, budget=budget,
        source=str(reference), original_score_sha256=sha(reference / 'score_manifest.json'),
        initial_labels_exact=True, terminal_labels_splits_exact=True, probability_rtol=1e-6,
        source_sha256=sha(__file__), job=os.environ['SLURM_JOB_ID'], completed_at=utc()))
    log('anchor_parity_passed', sample=sample, budget=budget)


def run(sample, budget, space):
    cfg = config(sample, budget, space);dest = Path(cfg['dest'])
    if checked(dest, 'fit_manifest.json', 'FIT_COMPLETE'):
        previous = json.loads((dest / 'fit_manifest.json').read_text())
        cache_compatibility.validate_completed_fit(sys.modules[__name__],cfg,previous)
        return
    x, cells, gm = get_geometry(cfg)
    z = make_embedding(cfg, x, cells, gm)
    cfg['conditions'] = make_partitions(cfg, z, cells)
    write_json(dest / 'config.json', cfg)
    call_r('score_candidates.R', dest / 'config.json')
    finish(cfg)
    manifests = {c['route']: sha(Path(c['dest']) / 'terminal/L00_mean/terminal_manifest.json') for c in cfg['conditions']}
    write_json(dest / 'fit_manifest.json', dict(status='fit_complete_evaluation_pending', config=cfg,
        terminal_manifest_hashes=manifests, reference_labels_used_for_fit=False,
        source_sha256=sha(__file__), job=os.environ['SLURM_JOB_ID'], completed_at=utc()))
    complete(dest, 'fit_manifest.json', 'FIT_COMPLETE')
    log('fit_complete', sample=sample, budget=budget, space=space, n_conditions=len(manifests))


def main():
    require_slurm()
    cache_compatibility.own_source(sys.modules[__name__])
    import terminal
    terminal.threads()
    parser = argparse.ArgumentParser();parser.add_argument('tasks');parser.add_argument('mode', nargs='?', default='run')
    args = parser.parse_args()
    tasks = json.loads(Path(args.tasks).read_text())
    task = tasks[int(os.environ.get('SLURM_ARRAY_TASK_ID', '0'))]
    try:
        if args.mode == 'parity':
            parity(task['sample'], task['budget'])
        elif args.mode == 'geometry':
            get_geometry(config(task['sample'], task['budget'], 'PCA2'))
        else:
            for space in task.get('spaces', [task.get('space')]):
                run(task['sample'], task['budget'], space)
                subprocess.run([sys.executable, str(CODE / 'evaluate_and_plot.py'),
                    str(OUTPUT / task['sample'] / task['budget'] / space)], check=True)
    except Exception as exc:
        directory = OUTPUT / task['sample'] / task['budget'] / 'failures';directory.mkdir(parents=True, exist_ok=True)
        write_json(directory / (os.environ['SLURM_JOB_ID'] + '_' + os.environ.get('SLURM_ARRAY_TASK_ID', '0') + '.json'),
                   dict(task=task, mode=args.mode, error=repr(exc), job=os.environ['SLURM_JOB_ID'], time=utc()))
        raise


if __name__ == '__main__':
    main()
