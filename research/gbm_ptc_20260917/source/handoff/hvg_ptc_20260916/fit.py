#!/usr/bin/env python3
"""Label-free GBM fits and terminal annotation. Reference evaluation is a separate job."""
import argparse
from contextlib import redirect_stdout, redirect_stderr
import gc
import json
import os
from pathlib import Path
import resource
import sys
import time
import traceback
import warnings

from common import (ROOT, OUT, MARKERS, MARKER, require_slurm, sha, utc,
                    version_record, runtime_record, write_json, key)


def code_sources():
    paths = [Path(__file__), Path(__file__).with_name('common.py'),
             OUT/'protocol/geometries.json', ROOT/'handoff/g274_table4/cgo_grid.py',
             ROOT/'handoff/gbm/darmanis_region_final.py',
             ROOT/'DGscRNA/dgscrna/core/deep_learning.py',
             ROOT/'DGscRNA/dgscrna/models/deep_model.py',
             MARKERS/f'{MARKER}.csv', MARKERS/'mapping_L1_v3.csv']
    return {str(p): sha(p) for p in paths}


def fit_metadata(model, caught):
    from fit_core import _fit_metadata
    result = _fit_metadata(model, caught)
    # HDBSCAN exposes a joblib.Memory object even when caching is disabled.
    # Preserve its actual configuration instead of trying to JSON-serialize the object.
    memory = result['params'].get('memory')
    if memory is not None and hasattr(memory, 'location'):
        result['params']['memory'] = {'type': type(memory).__name__, 'location': memory.location}
    return result


def prepared(sample):
    import numpy as np
    import anndata as ad
    import scipy.sparse as sp
    dest = OUT/'prepared'/sample
    m = json.loads((dest/'manifest.json').read_text())
    assert (dest/'PREPARED').read_text().strip() == sha(dest/'manifest.json')
    assert m['status'] == 'completed' and not m['gold_used_for_preprocessing']
    for name, h in m['outputs'].items():
        assert sha(dest/name) == h, f'Preparation checksum mismatch: {name}'
    cells, genes = [(dest/n).read_text().splitlines() for n in ['cells.tsv', 'genes.tsv']]
    a = ad.AnnData(sp.load_npz(dest/'lognorm.npz'))
    a.obs_names, a.var_names = cells, genes
    a.layers['lognorm'] = a.X
    a.obsm['X_scaled'] = np.load(dest/'scaled_all.npy', mmap_mode='r', allow_pickle=False)
    assert a.shape == (len(cells), len(genes))
    assert not len(a.obs.columns) and not len(a.var.columns)
    assert a.obs_names.is_unique and a.var_names.is_unique
    return a, m, dest


def representation(x, g):
    import numpy as np
    from sklearn.decomposition import PCA, FactorAnalysis, FastICA
    from sklearn.manifold import Isomap, TSNE
    started = time.perf_counter()
    info = {'geometry_feature_width': x.shape[1], 'input_shape': list(x.shape)}
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        if g['input_space'] == 'pca30':
            pre = PCA(n_components=30, svd_solver='arpack', random_state=g['seed'])
            x = np.asarray(pre.fit_transform(x), dtype=np.float32)
            info['pre_pca'] = fit_metadata(pre, [])
            info['pre_pca']['explained_variance_ratio'] = pre.explained_variance_ratio_.tolist()
        info['reducer_input_shape'] = list(x.shape)
        kwargs = dict(n_components=g['dim'])
        if g['dr'] == 'none':
            z, model = x, None
        elif g['dr'] == 'SCANPY_UMAP':
            import anndata as ad
            import scipy.sparse as sp
            import scanpy as sc
            b = ad.AnnData(sp.csr_matrix((len(x), 1), dtype=np.float32))
            b.obsm['X_pca'] = x
            sc.pp.neighbors(b, n_neighbors=g['neighbors'], use_rep='X_pca', random_state=g['seed'])
            sc.tl.umap(b, n_components=g['dim'], min_dist=g['min_dist'], random_state=g['seed'])
            z, model = b.obsm['X_umap'], None
            info['scanpy_neighbors_params'] = b.uns['neighbors']['params']
            info['scanpy_umap_params'] = b.uns['umap']['params']
        else:
            if g['dr'] == 'PCA':
                model = PCA(**kwargs, random_state=g['seed'],
                            svd_solver='arpack' if g['dim'] == 30 else 'auto')
            elif g['dr'] == 'FA':
                model = FactorAnalysis(**kwargs, max_iter=5000, random_state=g['seed'])
            elif g['dr'] == 'ICA':
                model = FastICA(**kwargs, max_iter=5000, tol=1e-4, whiten='unit-variance', random_state=g['seed'])
            elif g['dr'] == 'Isomap':
                model = Isomap(**kwargs, n_neighbors=g['neighbors'], n_jobs=1)
            elif g['dr'] == 'UMAP':
                import umap
                model = umap.UMAP(**kwargs, n_neighbors=g['neighbors'], min_dist=g['min_dist'],
                                  random_state=g['seed'], n_jobs=1)
            elif g['dr'] == 'TSNE':
                model = TSNE(**kwargs, random_state=g['seed'], init='pca', n_jobs=1)
            else:
                raise ValueError(g['dr'])
            z = model.fit_transform(x)
        info.update(fit_metadata(model, caught) if model is not None else
                    {'warnings': [f'{w.category.__name__}: {w.message}' for w in caught]})
    z = np.asarray(z, dtype=np.float32)
    assert z.shape == (len(x), g['dim'] if g['dim'] is not None else x.shape[1])
    assert np.isfinite(z).all()
    info.update(seconds=time.perf_counter()-started, output_shape=list(z.shape))
    if model is not None and g['dr'] == 'PCA':
        info['explained_variance_ratio'] = model.explained_variance_ratio_.tolist()
    return z, info


def clusters(z, arm, seed):
    import numpy as np
    from sklearn.cluster import KMeans
    from sklearn.mixture import GaussianMixture
    from fit_core import summarize_labels
    if arm['clusterer'] == 'KMeans':
        model = KMeans(n_clusters=arm['k'], n_init=10, random_state=seed)
    elif arm['clusterer'] == 'GMM':
        model = GaussianMixture(n_components=arm['k'], covariance_type=arm['covariance'],
                                reg_covar=1e-4, random_state=seed)
    else:
        import hdbscan
        model = hdbscan.HDBSCAN(min_cluster_size=arm['min_cluster_size'],
                                min_samples=arm['min_samples'], core_dist_n_jobs=1)
    start = time.perf_counter()
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        cl = np.asarray(model.fit_predict(z), dtype=np.int64)
    info = {**fit_metadata(model, caught), **summarize_labels(cl),
            'seconds': time.perf_counter()-start}
    _, size = np.unique(cl[cl != -1], return_counts=True)
    info['cluster_sizes'] = size.tolist()
    if 'E2_K_selection' in arm['families']:
        from sklearn.metrics import silhouette_score, calinski_harabasz_score, davies_bouldin_score
        info['label_free_selection'] = dict(
            silhouette=float(silhouette_score(z, cl, sample_size=min(2000, len(z)), random_state=42)),
            calinski_harabasz=float(calinski_harabasz_score(z, cl)),
            davies_bouldin=float(davies_bouldin_score(z, cl)))
    return cl, info


def annotate(a, cluster_ids, idx, arm, seed, dest):
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import anndata as ad
    import scipy.sparse as sp
    from cgo_grid import load_panels, score_panels, calls_at_cutoff, seed_labels, final_prediction
    from dgscrna.core.deep_learning import train_deep_model
    cl = np.asarray(['Noise' if c == -1 else str(c) for c in cluster_ids])
    order = [c for c in pd.unique(cl) if c != 'Noise']
    start = time.perf_counter()
    if len(order) < 2:
        # Predeclared structural abstention: no supported marker contrasts; never a fitted classifier.
        final = np.full(a.n_obs, 'Unknown')
        line = np.where(cl == 'Noise', 'noise_terminal', 'no_cluster_contrast')
        status = dict(dl_status='terminal_abstention', cluster_status='all_noise' if not order else 'lt2_clusters',
                      final_valid=False, terminal_valid=True, training_attempted=False,
                      training_executed=False, prediction_executed=False,
                      n_training=0, n_training_classes=0, n_pool=int((cl != 'Noise').sum()),
                      n_noise=int((cl == 'Noise').sum()), n_dl_assigned=0,
                      n_final_called=0, lineage_counts=pd.Series(line).value_counts().to_dict())
        status.update(seconds=time.perf_counter()-start, actual_dl_feature_width=a.n_vars)
        np.savez_compressed(dest/'predictions.npz', cluster=cluster_ids, seed=final, final=final, lineage=line)
        return status
    scoring_idx = np.arange(a.n_vars) if arm['scoring_features'] == 'all' else idx
    b = ad.AnnData(a.layers['lognorm'][:, scoring_idx].copy())
    b.obs_names = a.obs_names
    b.var_names = a.var_names[scoring_idx]
    b.obs['_cl'] = pd.Categorical(cl)
    b.layers['lognorm'] = b.X
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        sc.tl.rank_genes_groups(b, '_cl', groups=order, method='wilcoxon',
                                n_genes=100, layer='lognorm', use_raw=False)
    deg = {c: dict(zip(b.uns['rank_genes_groups']['names'][c],
                      b.uns['rank_genes_groups']['logfoldchanges'][c])) for c in order}
    pd.DataFrame([dict(cluster=c, gene=g, logfoldchange=float(v))
                  for c, d in deg.items() for g, v in d.items()]).to_csv(dest/'top100_degs.csv.gz', index=False)
    S, names = score_panels(deg, load_panels(b), order)
    pd.DataFrame(S, index=names, columns=order).to_csv(dest/'panel_scores.csv.gz')
    mapping = pd.read_csv(MARKERS/'mapping_L1_v3.csv')
    mapping = mapping[mapping.marker_set.eq(MARKER)].set_index('panel')['gold_class'].to_dict()
    raw_calls = dict(zip(order, calls_at_cutoff(S, names, arm['cutoff'])))
    initial = seed_labels(cl, raw_calls, mapping)
    del b
    def traced_train(subset, annotation_key, **kwargs):
        kwargs['random_state'] = seed
        assert subset.obsm['X_scaled'].shape[1] == a.n_vars
        assert not subset.obs['cluster'].eq('Noise').any()
        result = train_deep_model(subset, annotation_key, **kwargs)
        write_json(dest/'training_history.json', {k: v for k, v in result.items() if k != 'model'})
        return result
    with (dest/'refinement.log').open('w') as log, redirect_stdout(log), redirect_stderr(log):
        final, lineage, info = final_prediction(a, cl, initial, train_fn=traced_train)
    if info['dl_status'] == 'dl_error':
        # Technical failure remains resumable; its current prediction is never scored as final.
        write_json(dest/'refinement_error.json', info)
        raise RuntimeError(f'Refinement failed: {info.get("error", "unknown")}')
    # Store fixed-width Unicode, not pandas' object arrays (no pickle-dependent artifacts).
    final = np.asarray(final.to_numpy(), dtype=str)
    lineage = np.asarray(lineage.to_numpy(), dtype=str)
    assert (final[cl == 'Noise'] == 'Unknown').all()
    retained = ~np.isin(initial, ['Unknown', 'Undecided', 'Noise']) & (cl != 'Noise')
    assert np.array_equal(final[retained], initial[retained])
    info.update(terminal_valid=True, cluster_status='ok',
                n_final_called=int((~np.isin(final, ['Unknown', 'Undecided', 'Noise'])).sum()),
                seconds=time.perf_counter()-start, actual_scoring_feature_width=len(scoring_idx),
                actual_dl_feature_width=a.n_vars, raw_cluster_calls=raw_calls,
                deg_warnings=[str(w.message) for w in caught],
                terminal_note='structurally untrainable; known calls retained, pool Unknown'
                    if info['dl_status'] == 'lt2_training_classes' else 'terminal output of corrected refinement')
    np.savez_compressed(dest/'predictions.npz', cluster=cluster_ids,
                        seed=np.asarray(initial, dtype=str), final=final, lineage=lineage)
    return info


def valid_cache(path, sources, inputs):
    if not (path/'COMPLETE').exists():
        return None
    m = json.loads((path/'manifest.json').read_text())
    assert (path/'COMPLETE').read_text().strip() == sha(path/'manifest.json'), 'Manifest changed after completion'
    assert m['sources'] == sources and m['prepared_manifest_sha256'] == inputs, 'Cached scientific code/input changed'
    for n, h in m['outputs'].items():
        assert sha(path/n) == h, f'Cached result changed: {path/n}'
    return m


def run(sample, feature, selected=None):
    require_slurm()
    import numpy as np
    import torch
    from threadpoolctl import threadpool_limits
    torch.set_num_threads(min(4, int(os.environ.get('SLURM_CPUS_PER_TASK', 1))))
    torch.set_num_interop_threads(1)
    sys.path.insert(0, str(ROOT/'handoff/g274_table4'))
    a, prep, src = prepared(sample)
    idx = np.load(src/f'indices_{feature}.npy', allow_pickle=False)
    # The actual geometry array is explicitly subset; AnnData.obsm cannot bypass the mask.
    x = np.asarray(a.obsm['X_scaled'][:, idx], dtype=np.float32)
    assert x.shape == (a.n_obs, prep['feature_sets'][feature]['n_features'])
    x.setflags(write=False)
    inputs = sha(src/'manifest.json')
    sources = code_sources()
    gs = [g for g in json.loads((OUT/'protocol/geometries.json').read_text())
          if g['feature'] == feature and (selected is None or g['geometry_id'] in selected)]
    assert gs
    status = dict(sample=sample, feature=feature, started_at=utc(), **runtime_record(),
                  sources=sources, prepared_manifest_sha256=inputs, failures=[], geometries=[])
    taskdest = OUT/'task_status'/sample
    taskdest.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter()
    limits = threadpool_limits(limits=min(8, int(os.environ.get('SLURM_CPUS_PER_TASK', 1))))
    for g in gs:
        dest = OUT/'fits'/sample/g['geometry_id']
        dest.mkdir(parents=True, exist_ok=True)
        cache = valid_cache(dest, sources, inputs)
        if cache is not None:
            status['geometries'].append(dict(geometry_id=g['geometry_id'], status='verified_cache'))
            continue
        base = dict(sample=sample, geometry=g, sources=sources, prepared_manifest_sha256=inputs,
                    started_at=utc(), status='running', gold_used_for_fitting=False,
                    cell_ids_sha256=sha(src/'cells.tsv'), feature_ids_sha256=sha(src/f'indices_{feature}.npy'),
                    actual_geometry_feature_width=x.shape[1], versions=version_record(), **runtime_record())
        if (dest/'manifest.json').exists():
            old = json.loads((dest/'manifest.json').read_text())
            archive = dest/f'attempt_{old.get("slurm_job_id", "unknown")}.json'
            if not archive.exists():
                write_json(archive, old)
        write_json(dest/'manifest.json', base)
        try:
            # Keep completed embeddings on retry only when their recorded code/input hashes match.
            embpath = dest/'embedding_manifest.json'
            em = json.loads(embpath.read_text()) if embpath.exists() else None
            if em and em['sources'] == sources and em['prepared_manifest_sha256'] == inputs:
                if g['dr'] == 'none':
                    z = x
                else:
                    assert sha(dest/'embedding.npy') == em['sha256']
                    z = np.load(dest/'embedding.npy', allow_pickle=False)
                ri = em['reducer_info']
            else:
                z, ri = representation(x, g)
                if g['dr'] != 'none':
                    np.save(dest/'embedding.npy', z, allow_pickle=False)
                write_json(embpath, dict(sources=sources, prepared_manifest_sha256=inputs,
                    sha256=sha(dest/'embedding.npy') if g['dr'] != 'none' else None, reducer_info=ri))
            base['reducer'] = ri
            partition_cache = {}
            failed_arms = []
            for arm in g['arms']:
                ap = dest/arm['arm_id']
                ap.mkdir(parents=True, exist_ok=True)
                am = valid_cache(ap, sources, inputs)
                if am is not None:
                    continue
                am = dict(sample=sample, geometry_id=g['geometry_id'], arm=arm, sources=sources,
                          prepared_manifest_sha256=inputs, started_at=utc(),
                          status='running', gold_used_for_fitting=False, **runtime_record())
                write_json(ap/'manifest.json', am)
                try:
                    pk = key({k: arm[k] for k in ['clusterer', 'k', 'covariance', 'min_cluster_size', 'min_samples']})
                    if pk not in partition_cache:
                        partition_cache[pk] = clusters(z, arm, g['seed'])
                    cl, ci = partition_cache[pk]
                    # Persist clustering before any annotation work, to retain a valid partition on DL failure.
                    np.save(ap/'clusters.npy', cl, allow_pickle=False)
                    am['clustering'] = ci
                    write_json(ap/'manifest.json', am)
                    ai = annotate(a, cl, idx, arm, g['seed'], ap)
                    am.update(annotation=ai, status='completed', finished_at=utc(),
                              outputs={p.name: sha(p) for p in ap.iterdir() if p.is_file()
                                       and p.name not in ['manifest.json', 'COMPLETE']})
                    write_json(ap/'manifest.json', am)
                    (ap/'COMPLETE').write_text(sha(ap/'manifest.json')+'\n')
                    print(f'{sample} {feature} {g["dr"]}{g["dim"]} {g["input_space"]} '
                          f'{arm["clusterer"]} {ai["dl_status"]} '
                          f'clusters={ci["n_clusters"]} noise={ci["noise_frac"]:.3f} '
                          f'{time.perf_counter()-start:.0f}s', flush=True)
                except Exception as exc:
                    am.update(status='failed', error=str(exc), traceback=traceback.format_exc(), finished_at=utc())
                    write_json(ap/'manifest.json', am)
                    failed_arms.append(arm['arm_id'])
                    print(am['traceback'], flush=True)
                gc.collect()
            if failed_arms:
                raise RuntimeError(f'{len(failed_arms)} arms failed: {failed_arms}')
            base.update(status='completed', finished_at=utc(),
                        peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024,
                        outputs={p.name: sha(p) for p in dest.iterdir() if p.is_file()
                                 and p.name not in ['manifest.json', 'COMPLETE'] and not p.name.startswith('attempt_')})
            write_json(dest/'manifest.json', base)
            (dest/'COMPLETE').write_text(sha(dest/'manifest.json')+'\n')
            status['geometries'].append(dict(geometry_id=g['geometry_id'], status='completed'))
            del z, partition_cache
            gc.collect()
        except Exception as exc:
            base.update(status='failed', finished_at=utc(), error=str(exc), traceback=traceback.format_exc())
            write_json(dest/'manifest.json', base)
            status['failures'].append(dict(geometry_id=g['geometry_id'], error=str(exc)))
            print(base['traceback'], flush=True)
        status.update(elapsed_seconds=time.perf_counter()-start, updated_at=utc())
        write_json(taskdest/f'{feature}.json', status)
    status.update(status='completed' if not status['failures'] else 'failed', finished_at=utc(),
                  elapsed_seconds=time.perf_counter()-start)
    write_json(taskdest/f'{feature}.json', status)
    if status['failures']:
        raise RuntimeError(f'{sample}/{feature}: {len(status["failures"])} geometries failed')


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--sample')
    p.add_argument('--feature')
    p.add_argument('--task', type=int)
    p.add_argument('--select-file', type=Path)
    args = p.parse_args()
    if args.task is not None:
        task = json.loads((OUT/'protocol/tasks.json').read_text())[args.task]
        args.sample, args.feature = task['sample'], task['feature']
    selected = args.select_file.read_text().split() if args.select_file else None
    run(args.sample, args.feature, selected)
