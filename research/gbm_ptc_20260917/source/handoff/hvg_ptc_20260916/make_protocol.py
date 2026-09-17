#!/usr/bin/env python3
"""Enumerate approved conditions using only stdlib; does not access expression/labels."""
import json
from collections import Counter
from pathlib import Path
from common import FEATURES, SEEDS, OUT, ROOT, key, samples, write_json


def build():
    geometries = {}

    def add(feature, dr, dim=2, input_space='genes', seed=42, neighbors=15, min_dist=.1,
            clusterers=('HDBSCAN',), mcs=15, ms=15, k=23, score_features='all', family='E2'):
        g = dict(feature=feature, dr=dr, dim=None if dr == 'none' else dim,
                 input_space=input_space, seed=seed, neighbors=neighbors, min_dist=min_dist)
        gid = key(g)
        if gid not in geometries:
            geometries[gid] = dict(geometry_id=gid, **g, arms={})
        for cl in clusterers:
            a = dict(clusterer=cl, k=k if cl != 'HDBSCAN' else None,
                     covariance='diag' if cl == 'GMM' else None,
                     min_cluster_size=mcs if cl == 'HDBSCAN' else None,
                     min_samples=ms if cl == 'HDBSCAN' else None,
                     scoring_features=score_features, marker='CM2_glioma_other', cutoff='none')
            aid = key(a)
            row = geometries[gid]['arms'].setdefault(aid, dict(arm_id=aid, **a, families=[]))
            if family not in row['families']:
                row['families'].append(family)

    # Original Table 4 shape, all three clusterers, six fixed feature levels.
    for f in FEATURES:
        for dr in ['PCA', 'FA', 'ICA', 'Isomap', 'UMAP', 'TSNE', 'none']:
            add(f, dr, clusterers=('KMeans', 'GMM', 'HDBSCAN'), family='E1_primary_table')
    for f in FEATURES:
        for space in ['genes', 'pca30']:
            for dim in [2, 10, 30]:
                add(f, 'UMAP', dim=dim, input_space=space, family='E2_factorial')
        add(f, 'PCA', dim=30, family='E2_factorial')
    # Fixed, predeclared seed comparison, not post-hoc winners.
    seed_paths = [('PCA', 2, 'genes', 'HDBSCAN'), ('ICA', 2, 'genes', 'HDBSCAN'),
                  ('UMAP', 2, 'genes', 'HDBSCAN'), ('UMAP', 2, 'pca30', 'HDBSCAN'),
                  ('PCA', 30, 'genes', 'HDBSCAN'), ('none', 2, 'genes', 'KMeans')]
    for f in ['all', 'hvg2000']:
        for seed in SEEDS:
            for dr, dim, space, cl in seed_paths:
                add(f, dr, dim=dim, input_space=space, seed=seed, clusterers=(cl,), family='E2_seeds')
        for space in ['genes', 'pca30']:
            for nn in [15, 30, 50]:
                add(f, 'UMAP', dim=10, input_space=space, neighbors=nn, family='E2_neighbors')
            for md in [0., .1, .5]:
                add(f, 'UMAP', dim=10, input_space=space, min_dist=md, family='E2_min_dist')
            for mcs in [5, 15, 30, 60]:
                add(f, 'UMAP', dim=10, input_space=space, mcs=mcs, family='E2_HDB_mcs')
            for ms in [5, 15, 30]:
                add(f, 'UMAP', dim=10, input_space=space, ms=ms, family='E2_HDB_ms')
        for mcs in [5, 15, 30, 60]:
            add(f, 'PCA', dim=30, mcs=mcs, family='E2_HDB_mcs')
        for ms in [5, 15, 30]:
            add(f, 'PCA', dim=30, ms=ms, family='E2_HDB_ms')
        for k in [10, 15, 23, 30, 40]:
            add(f, 'PCA', dim=30, k=k, clusterers=('KMeans', 'GMM'), family='E2_K_selection')
    for f in ['seurat2000', 'seurat5000']:
        for dr in ['PCA', 'FA', 'ICA', 'Isomap', 'UMAP', 'TSNE', 'none']:
            add(f, dr, clusterers=('KMeans', 'GMM', 'HDBSCAN'), family='E2_HVG_flavor_bridge')
        add(f, 'UMAP', input_space='pca30', family='E2_HVG_flavor_bridge')
        add(f, 'PCA', dim=30, family='E2_HVG_flavor_bridge')
        add(f, 'SCANPY_UMAP', input_space='pca30', min_dist=.5, family='E2_legacy_geometry_bridge')
    for dr, dim, space in [('UMAP', 2, 'genes'), ('UMAP', 2, 'pca30'), ('PCA', 30, 'genes')]:
        add('hvg2000_markers', dr, dim=dim, input_space=space, family='E3_marker_union')
        add('hvg2000', dr, dim=dim, input_space=space, score_features='geometry', family='E3_scoring_truncation')
    rows = []
    for g in geometries.values():
        g['arms'] = list(g['arms'].values())
        rows.append(g)
    rows.sort(key=lambda g: (FEATURES.index(g['feature']) if g['feature'] in FEATURES else 6,
                             g['feature'], g['seed'] != 42, g['geometry_id']))
    return rows


if __name__ == '__main__':
    dest = OUT / 'protocol'
    dest.mkdir(parents=True, exist_ok=True)
    rows = build()
    roster = samples()
    assert len(roster) == 121 and len(set(roster)) == 121
    families = Counter(f for g in rows for a in g['arms'] for f in a['families'])
    features = list(dict.fromkeys(g['feature'] for g in rows))
    write_json(dest / 'geometries.json', rows)
    primary_order = ['all', 'hvg2000', 'hvg500', 'hvg1000', 'hvg3000', 'hvg5000']
    task_pairs = [(s, f) for s in roster for f in primary_order]
    task_pairs += [(s, f) for s in roster for f in features if f not in primary_order]
    tasks = [dict(task_id=i, sample=s, feature=f) for i, (s, f) in enumerate(task_pairs)]
    write_json(dest / 'tasks.json', tasks)
    (dest / 'samples.txt').write_text('\n'.join(roster) + '\n')
    summary = dict(status='approved_with_user_amendment', independent_Pu_cohort=False,
        priority='GBM_delivery_then_PTC', n_samples=len(roster), n_geometries_per_sample=len(rows),
        n_partition_conditions_per_sample=sum(len(g['arms']) for g in rows),
        n_total_conditions=len(roster)*sum(len(g['arms']) for g in rows),
        family_conditions_per_sample=families, n_sample_feature_tasks=len(tasks),
        seeds=SEEDS, primary_features=FEATURES,
        stage='terminal_DL_refinement', stage_clarification='user explicitly confirmed RL means DL/refinement final output',
        dl_feature_policy='all gene-filtered scaled genes, fixed across geometry feature arms',
        annotation_protocol='legacy top100 Wilcoxon/logFC density scores and L1 seed mapping; corrected terminal Noise semantics',
        no_reference_labels_in_fit=True, final_scores_only=True,
        classifier_status_required=['executed','empty_pool','invalid','terminal_abstention'],
        initial_primary_hypothesis='HVG2000 versus all, direct UMAP2 HDBSCAN15/15',
        geometry_note='direct UMAP min_dist0.1; separate Scanpy PCA30/UMAP2 min_dist0.5 legacy bridge',
        source_plan=str(ROOT/'handoff/plan_20260916_hvg_ptc/USER_UPDATE.md'))
    write_json(dest/'protocol.json',summary)
    print(json.dumps(summary, indent=2))
