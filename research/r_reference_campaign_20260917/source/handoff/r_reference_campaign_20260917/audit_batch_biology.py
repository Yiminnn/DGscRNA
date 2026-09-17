"""Post-fit PTC neighborhood diagnostics; expression and annotation fits are untouched.

All data access and numerical operations require SLURM. Exact Euclidean neighbors
are queried for the same seeded, population/sample-stratified cells in all arms.
"""
import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT = ROOT / 'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'
DEST = OUT / 'verification/batch_biology'
UNITS = [f'PTC_{g}_{c}' for g in ['NMT', 'TTU']
         for c in ['CCA2000', 'NONE_RNA2000', 'HARMONY_RNA2000']]
EXCLUDED = {'Unknown', 'Other', 'Lymphoid_ambiguous', 'Tumor_unspecified'}
K = 30
MIN_BATCH_CELLS = 10
QUERY_CAP = 128
SEED = 42


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def sha_ids(values):
    return hashlib.sha256(('\n'.join(map(str, values)) + '\n').encode()).hexdigest()


def exact_neighbors(coordinates, reference_indices, query_indices, k=K):
    """Exact cKDTree queries with self removed by identity, including tied zeros.

    Reference input is sorted by canonical global index. Distances and returned
    global indices break returned ties deterministically; no coordinate jitter.
    A tie at the kth boundary is counted and reported, never interpreted as a
    unique graph. The result is deterministic in the recorded scipy version.
    """
    import numpy as np
    from scipy.spatial import cKDTree
    r = np.sort(np.asarray(reference_indices, dtype=np.int64))
    q = np.asarray(query_indices, dtype=np.int64)
    assert len(r) > k and np.isin(q, r).all()
    tree = cKDTree(coordinates[r])
    ask = min(k + 2, len(r))
    distance, index = tree.query(coordinates[q], k=ask, eps=0, p=2,
                                 workers=int(os.environ.get('SLURM_CPUS_PER_TASK', '1')))
    indices = r[index]
    selected = np.empty((len(q), k), dtype=np.int64)
    distances = np.empty((len(q), k), dtype=np.float64)
    tied_boundary = np.zeros(len(q), dtype=bool)
    for row, qi in enumerate(q):
        mask = indices[row] != qi
        ix, dd = indices[row, mask], distance[row, mask]
        order = np.lexsort((ix, dd))
        ix, dd = ix[order], dd[order]
        assert len(ix) >= k
        selected[row], distances[row] = ix[:k], dd[:k]
        if len(dd) > k:
            tied_boundary[row] = dd[k - 1] == dd[k]
    assert not (selected == q[:, None]).any()
    return selected, distances, tied_boundary


def self_test():
    """Check scientific edge cases against direct distance calculation."""
    import numpy as np
    rng = np.random.default_rng(SEED)
    x = rng.normal(size=(91, 7))
    q = np.array([0, 3, 15, 90])
    actual, _, _ = exact_neighbors(x, np.arange(len(x)), q, k=9)
    for row, qi in enumerate(q):
        d = ((x - x[qi]) ** 2).sum(axis=1)
        d[qi] = np.inf
        expected = np.lexsort((np.arange(len(x)), d))[:9]
        assert np.array_equal(actual[row], expected)
    duplicate = np.zeros((45, 2))
    actual, distance, tied = exact_neighbors(duplicate, np.arange(45), np.arange(45), k=30)
    assert not (actual == np.arange(45)[:, None]).any()
    assert not distance.any() and tied.all()
    # Enumerating every other cell must give normalized mixing exactly one,
    # including a strongly imbalanced batch composition.
    batch = np.array(['a'] * 7 + ['b'] * 3)
    cross = np.array([(batch[np.arange(10) != i] != batch[i]).mean() for i in range(10)])
    expectation = np.array([(10 - (batch == batch[i]).sum()) / 9 for i in range(10)])
    np.testing.assert_allclose(cross / expectation, 1)
    assert np.isfinite(cross / expectation).all()
    return {'status': 'passed', 'checks': ['exact neighbors versus brute-force distances',
            'self removal with duplicate coordinates and tied boundaries',
            'leave-self-out mixing expectation with imbalanced batches']}


def query_selection(meta, cap):
    import numpy as np
    rows = []
    for (prefix, population, sample), sub in meta.groupby(['prefix', 'population', 'sample'], sort=True):
        if population in EXCLUDED:
            continue
        indices = sub.index.to_numpy(dtype=np.int64)
        seed_bytes = hashlib.sha256(f'{SEED}|{prefix}|{population}|{sample}'.encode()).digest()
        rng = np.random.default_rng(int.from_bytes(seed_bytes[:8], 'little'))
        chosen = np.sort(rng.choice(indices, size=min(cap, len(indices)), replace=False))
        for idx in chosen:
            rows.append({'index': int(idx), 'cell_id': meta.loc[idx, 'cell_id'],
                         'prefix': prefix, 'population': population, 'sample': sample,
                         'stratum_cells': len(indices), 'stratum_queries': len(chosen),
                         'sampling_weight': len(indices) / len(chosen)})
    import pandas as pd
    return pd.DataFrame(rows).sort_values('index').reset_index(drop=True)


def audit_unit(unit, pilot=False):
    import numpy as np
    import pandas as pd
    import scipy
    started = time.monotonic()
    assert unit in UNITS
    prep = OUT / 'PTC_ablation' / unit
    preparation = json.loads((prep / 'prepare_manifest.json').read_text())
    baseline = ROOT / 'results/hvg_ptc_20260916_v1/ptc_paper_baseline'
    refpath = baseline / 'paper_baseline_reference.csv.gz'
    tcrpath = baseline / 'original_DG_binary_pairs_for_R.csv.gz'
    cells = pd.read_csv(prep / 'cells.csv', dtype=str, keep_default_na=False)
    ref = pd.read_csv(refpath, dtype=str, keep_default_na=False).set_index('cell_id')
    tcr = pd.read_csv(tcrpath, keep_default_na=False).set_index('cell_id')
    assert cells.cell_id.is_unique and ref.index.is_unique and tcr.index.is_unique
    aligned = ref.loc[cells.cell_id].reset_index()
    assert aligned.cell_id.tolist() == cells.cell_id.tolist()
    assert aligned['sample'].tolist() == cells.batch.tolist()
    sys.path.insert(0, str(ROOT / 'handoff/ptc_recovery_20260916'))
    from label_rules import broad_lineage
    meta = pd.DataFrame({'cell_id': cells.cell_id, 'sample': cells.batch,
                         'native': aligned.paper_final_native})
    meta['population'] = meta.native.map({v: broad_lineage(v) for v in meta.native.unique()})
    meta['prefix'] = meta['sample'].str.split('-').str[0]
    meta['paper_tissue'] = meta.prefix.map({'N': 'Normal adjacent', 'MT': 'Metastasis',
                                          'T': 'Primary tumor', 'TU': 'Primary tumor'})
    assert meta.paper_tissue.notna().all()
    meta['tcr_detected'] = tcr.loc[cells.cell_id, 'truth'].astype(int).to_numpy()
    assert set(meta.tcr_detected) <= {0, 1}
    cap = 16 if pilot else QUERY_CAP
    queries = query_selection(meta, cap)
    dest = DEST / ('pilot' if pilot else 'units') / unit
    dest.mkdir(parents=True, exist_ok=True)
    queries.to_csv(dest / 'query_cells.csv.gz', index=False)
    meta.groupby(['prefix', 'population', 'sample'], sort=True).size().rename('n_cells').reset_index().to_csv(
        dest / 'fixed_population_counts.csv', index=False)
    records, eligibility, graphs = [], [], {}
    inputs = {str(p): sha(p) for p in [prep / 'cells.csv', prep / 'prepare_manifest.json',
                                      prep / 'DL_features.txt', refpath, tcrpath]}
    sample = meta['sample'].to_numpy()
    population = meta.population.to_numpy()
    native = meta.native.to_numpy()
    detected = meta.tcr_detected.to_numpy()
    paper_tissue = meta.paper_tissue.to_numpy()
    graph_number = 0
    for space in ['PCA30', 'UMAP2']:
        path = prep / (space + '.csv')
        frame = pd.read_csv(path, index_col=0)
        assert frame.index.tolist() == cells.cell_id.tolist(), f'{unit}/{space} cell order changed'
        coordinates = frame.to_numpy(dtype=np.float64)
        assert coordinates.shape[1] == (30 if space == 'PCA30' else 2)
        assert np.isfinite(coordinates).all()
        inputs[str(path)] = sha(path)
        for prefix, sub in meta.groupby('prefix', sort=True):
            reference_indices = sub.index.to_numpy(dtype=np.int64)
            selected_queries = queries[queries.prefix.eq(prefix)]
            qi = selected_queries['index'].to_numpy(dtype=np.int64)
            neighbors, distance, tied = exact_neighbors(coordinates, reference_indices, qi)
            graph_number += 1
            graphs[f'g{graph_number}_query'] = qi
            graphs[f'g{graph_number}_neighbors'] = neighbors
            graphs[f'g{graph_number}_distance'] = distance
            graphs[f'g{graph_number}_name'] = np.array(f'{space}|label_retention|{prefix}')
            sizes = sub.population.value_counts()
            native_sizes = sub.native.value_counts()
            n = len(sub)
            for row, query in enumerate(selected_queries.itertuples(index=False)):
                i = query.index
                broad_same = float((population[neighbors[row]] == population[i]).mean())
                native_same = float((native[neighbors[row]] == native[i]).mean())
                broad_chance = (int(sizes[population[i]]) - 1) / (n - 1)
                native_chance = (int(native_sizes[native[i]]) - 1) / (n - 1)
                tcr_fraction = float(detected[neighbors[row]].mean())
                tcr_chance = (int(sub.tcr_detected.sum()) - int(detected[i])) / (n - 1)
                records.append(dict(unit=unit, space=space, diagnostic='label_retention',
                    cell_id=query.cell_id, prefix=prefix, population=query.population,
                    sample=query.sample, sampling_weight=query.sampling_weight,
                    stratum_cells=query.stratum_cells, stratum_queries=query.stratum_queries,
                    archived_broad_purity=broad_same, archived_native_purity=native_same,
                    broad_chance_purity=broad_chance, native_chance_purity=native_chance,
                    broad_chance_adjusted_purity=(broad_same - broad_chance) / (1 - broad_chance)
                        if broad_chance < 1 else None,
                    native_chance_adjusted_purity=(native_same - native_chance) / (1 - native_chance)
                        if native_chance < 1 else None,
                    tcr_detected=int(detected[i]), neighbor_TCR_detection_fraction=tcr_fraction,
                    TCR_detection_chance=tcr_chance,
                    TCR_detection_enrichment=tcr_fraction / tcr_chance if tcr_chance > 0 else None,
                    kth_distance_tied=bool(tied[row])))
        for (prefix, pop), sub in meta.groupby(['prefix', 'population'], sort=True):
            counts = sub['sample'].value_counts()
            reason = ('excluded_ambiguous_archived_population' if pop in EXCLUDED else
                      'not_shared_across_two_samples' if len(counts) < 2 else
                      'fewer_than_10_cells_in_a_sample' if int(counts.min()) < MIN_BATCH_CELLS else
                      'too_few_population_cells_for_k30' if len(sub) <= K else 'eligible')
            eligibility.append(dict(unit=unit, space=space, prefix=prefix, population=pop,
                n_cells=len(sub), n_samples=len(counts), min_cells_per_sample=int(counts.min()),
                sample_counts=json.dumps(counts.to_dict(), sort_keys=True), status=reason))
            if reason != 'eligible':
                continue
            selected_queries = queries[queries.prefix.eq(prefix) & queries.population.eq(pop)]
            qi = selected_queries['index'].to_numpy(dtype=np.int64)
            neighbors, distance, tied = exact_neighbors(coordinates, sub.index.to_numpy(), qi)
            graph_number += 1
            graphs[f'g{graph_number}_query'] = qi
            graphs[f'g{graph_number}_neighbors'] = neighbors
            graphs[f'g{graph_number}_distance'] = distance
            graphs[f'g{graph_number}_name'] = np.array(f'{space}|batch_mixing|{prefix}|{pop}')
            for row, query in enumerate(selected_queries.itertuples(index=False)):
                observed = float((sample[neighbors[row]] != query.sample).mean())
                expectation = (len(sub) - int(counts[query.sample])) / (len(sub) - 1)
                records.append(dict(unit=unit, space=space, diagnostic='batch_mixing',
                    cell_id=query.cell_id, prefix=prefix, population=pop, sample=query.sample,
                    sampling_weight=query.sampling_weight, stratum_cells=query.stratum_cells,
                    stratum_queries=query.stratum_queries, observed_cross_sample_fraction=observed,
                    expected_cross_sample_fraction=expectation,
                    normalized_cross_sample_mixing=observed / expectation,
                    kth_distance_tied=bool(tied[row])))
        # Biological tissue preservation is a distinct, confounded diagnostic.
        # NMT has normal-adjacent/metastatic tissue; TTU has only primary tumor.
        for pop, sub in meta.groupby('population', sort=True):
            if pop in EXCLUDED or sub.paper_tissue.nunique() < 2 or len(sub) <= K:
                continue
            selected_queries = queries[queries.population.eq(pop)]
            qi = selected_queries['index'].to_numpy(dtype=np.int64)
            neighbors, distance, tied = exact_neighbors(coordinates, sub.index.to_numpy(), qi)
            counts = sub.paper_tissue.value_counts()
            graph_number += 1
            graphs[f'g{graph_number}_query'] = qi
            graphs[f'g{graph_number}_neighbors'] = neighbors
            graphs[f'g{graph_number}_distance'] = distance
            graphs[f'g{graph_number}_name'] = np.array(f'{space}|tissue_retention|{pop}')
            for row, query in enumerate(selected_queries.itertuples(index=False)):
                same = float((paper_tissue[neighbors[row]] == paper_tissue[query.index]).mean())
                chance = (int(counts[paper_tissue[query.index]]) - 1) / (len(sub) - 1)
                records.append(dict(unit=unit, space=space, diagnostic='tissue_retention',
                    cell_id=query.cell_id, prefix=query.prefix, population=pop, sample=query.sample,
                    sampling_weight=query.sampling_weight, stratum_cells=query.stratum_cells,
                    stratum_queries=query.stratum_queries, paper_tissue=paper_tissue[query.index],
                    paper_tissue_purity=same, paper_tissue_chance_purity=chance,
                    paper_tissue_chance_adjusted_purity=(same - chance) / (1 - chance),
                    kth_distance_tied=bool(tied[row])))
        print(unit, space, 'diagnostics complete', round(time.monotonic() - started, 1), 'seconds', flush=True)
    pd.DataFrame(records).to_csv(dest / 'query_diagnostics.csv.gz', index=False)
    pd.DataFrame(eligibility).to_csv(dest / 'mixing_population_eligibility.csv', index=False)
    graphs['cell_ids'] = cells.cell_id.to_numpy(dtype=str)
    np.savez_compressed(dest / 'exact_neighbor_graphs.npz', **graphs)
    report = dict(status='pilot_complete' if pilot else 'complete', unit=unit,
        job=os.environ['SLURM_JOB_ID'], elapsed_seconds=time.monotonic() - started,
        k=K, seed=SEED, max_queries_per_population_sample=cap, n_cells=len(cells), n_queries=len(queries),
        query_cell_ids_sha256=sha_ids(queries.cell_id), cell_order_sha256=sha_ids(cells.cell_id),
        correction=preparation['correction'], features=preparation['features'],
        DL_input_sha256=preparation['DL_binary_sha256'],
        DL_feature_ids_sha256=sha(prep / 'DL_features.txt'),
        coordinate_metric='exact Euclidean; unchanged stored coordinates; no scaling or jitter',
        inputs=inputs, source_sha256=sha(__file__),
        label_rules_sha256=sha(ROOT / 'handoff/ptc_recovery_20260916/label_rules.py'),
        versions={'numpy': np.__version__, 'pandas': pd.__version__, 'scipy': scipy.__version__},
        evaluation_only=True, no_expression_or_annotation_refit=True,
        self_test=self_test(),
        outputs={p.name: sha(p) for p in dest.iterdir() if p.is_file() and p.name not in ['manifest.json', 'COMPLETE']})
    (dest / 'manifest.json').write_text(json.dumps(report, indent=2) + '\n')
    (dest / 'COMPLETE').write_text(sha(dest / 'manifest.json') + '\n')
    print(json.dumps(report, indent=2), flush=True)


def summarize():
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    summaries, strata, inputs, manifests = [], [], {}, []
    metrics = {'batch_mixing': ['observed_cross_sample_fraction', 'expected_cross_sample_fraction',
                              'normalized_cross_sample_mixing'],
               'label_retention': ['archived_broad_purity', 'archived_native_purity',
                                  'broad_chance_adjusted_purity', 'native_chance_adjusted_purity'],
               'tissue_retention': ['paper_tissue_purity', 'paper_tissue_chance_adjusted_purity']}
    for unit in UNITS:
        dest = DEST / 'units' / unit
        manifest = json.loads((dest / 'manifest.json').read_text())
        assert manifest['status'] == 'complete' and manifest['self_test']['status'] == 'passed'
        assert (dest / 'COMPLETE').read_text().strip() == sha(dest / 'manifest.json')
        for filename, checksum in manifest['outputs'].items():
            assert sha(dest / filename) == checksum
        manifests.append(manifest)
        inputs[unit] = sha(dest / 'manifest.json')
        df = pd.read_csv(dest / 'query_diagnostics.csv.gz')
        eligibility = pd.read_csv(dest / 'mixing_population_eligibility.csv')
        for space in ['PCA30', 'UMAP2']:
            row = dict(unit=unit, group=unit.split('_')[1], condition=unit.split('_', 2)[2],
                       space=space, n_cells=manifest['n_cells'], queried_cells=manifest['n_queries'],
                       query_fraction=manifest['n_queries'] / manifest['n_cells'])
            for diagnostic, names in metrics.items():
                frame = df[df.space.eq(space) & df.diagnostic.eq(diagnostic)]
                if not len(frame):
                    for name in names:
                        row[name + '_population_sample_balanced'] = None
                        row[name + '_cell_weighted'] = None
                    continue
                grouped = frame.groupby(['prefix', 'population', 'sample'], sort=True)
                table = grouped[names].mean().reset_index()
                table['unit'], table['space'], table['diagnostic'] = unit, space, diagnostic
                table['query_cells'] = grouped.size().to_numpy()
                table['represented_cells'] = grouped.stratum_cells.first().to_numpy()
                strata.append(table)
                # Average samples within a population/prefix, then equal weight
                # population/prefix strata. This is not a patient-level CI.
                means = table.groupby(['prefix', 'population'], sort=True)[names].mean().mean()
                for name in names:
                    row[name + '_population_sample_balanced'] = float(means[name])
                    valid = frame[name].notna()
                    row[name + '_cell_weighted'] = float(np.average(frame.loc[valid, name],
                        weights=frame.loc[valid, 'sampling_weight'])) if valid.any() else None
                row[diagnostic + '_query_cells'] = len(frame)
                row[diagnostic + '_tied_neighbor_boundary_fraction'] = float(frame.kth_distance_tied.mean())
            eligible = eligibility[eligibility.space.eq(space) & eligibility.status.eq('eligible')]
            row['mixing_eligible_population_prefixes'] = len(eligible)
            row['mixing_represented_cells'] = int(eligible.n_cells.sum())
            row['mixing_represented_cell_fraction'] = row['mixing_represented_cells'] / row['n_cells']
            known = df[df.space.eq(space) & df.diagnostic.eq('label_retention')]
            row['label_retention_represented_cells'] = int(known.groupby(['prefix', 'population', 'sample']).stratum_cells.first().sum())
            row['label_retention_represented_cell_fraction'] = row['label_retention_represented_cells'] / row['n_cells']
            positive = known[known.tcr_detected.eq(1)]
            row['TCR_positive_queries'] = len(positive)
            for name in ['neighbor_TCR_detection_fraction', 'TCR_detection_enrichment']:
                row['TCR_positive_' + name] = float(np.average(positive[name], weights=positive.sampling_weight)) if len(positive) else None
            summaries.append(row)
    for group in ['NMT', 'TTU']:
        subset = [m for m in manifests if m['unit'].startswith('PTC_' + group + '_')]
        assert len({m['query_cell_ids_sha256'] for m in subset}) == 1, 'Correction arms used different query cells'
        assert len({m['cell_order_sha256'] for m in subset}) == 1, 'Correction arms used different cell order'
        assert len({m['label_rules_sha256'] for m in subset}) == 1
        matched = [m for m in subset if m['unit'].endswith(('NONE_RNA2000', 'HARMONY_RNA2000'))]
        assert len(matched) == 2
        assert len({m['DL_input_sha256'] for m in matched}) == 1
        assert len({m['DL_feature_ids_sha256'] for m in matched}) == 1
    frame = pd.DataFrame(summaries)
    frame.to_csv(DEST / 'PTC_batch_biology_unit_summary.csv', index=False)
    pd.concat(strata, ignore_index=True).to_csv(DEST / 'PTC_batch_biology_population_sample_summary.csv', index=False)
    differences = []
    value_columns = [c for c in frame if c.endswith(('_population_sample_balanced', '_cell_weighted'))]
    for (group, space), subset in frame.groupby(['group', 'space']):
        subset = subset.set_index('condition')
        for baseline, condition, meaning in [('NONE_RNA2000', 'HARMONY_RNA2000', 'matched RNA inputs; geometry correction'),
                                             ('CCA2000', 'NONE_RNA2000', 'combined expression/gene-universe/geometry comparison')]:
            row = dict(group=group, space=space, baseline=baseline, condition=condition, interpretation=meaning)
            for column in value_columns:
                row['delta_' + column] = subset.loc[condition, column] - subset.loc[baseline, column]
            differences.append(row)
    pd.DataFrame(differences).to_csv(DEST / 'PTC_batch_biology_paired_changes.csv', index=False)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 8, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(2, 2, figsize=(8.3, 6.0), sharex=True, sharey='row')
    conditions = ['NONE_RNA2000', 'HARMONY_RNA2000', 'CCA2000']
    plotted = [('normalized_cross_sample_mixing_population_sample_balanced', 'Normalized cross-sample mixing\n(observed / composition expectation)'),
               ('archived_broad_purity_population_sample_balanced', 'Archived broad-label\nneighbor concordance')]
    for col, group in enumerate(['NMT', 'TTU']):
        for row, (metric, label) in enumerate(plotted):
            ax = axes[row, col]
            for space, color, symbol in [('PCA30', '#0072B2', 'o'), ('UMAP2', '#D55E00', 's')]:
                selected = frame[frame.group.eq(group) & frame.space.eq(space)].set_index('condition').loc[conditions]
                ax.plot(range(3), selected[metric], marker=symbol, color=color, label=space)
            if row == 0:
                ax.axhline(1, color='#555555', lw=.8, ls='--')
                ax.set_ylim(0, max(1.05, float(frame[metric].max()) * 1.05))
                ax.set_title(group + ' | fixed populations and biological strata')
            else:
                ax.set_ylim(0, 1)
            ax.set_xticks(range(3), ['RNA / none', 'RNA / Harmony', 'CCA2000'])
            ax.grid(axis='y', alpha=.2)
            ax.spines[['top', 'right']].set_visible(False)
            if col == 0:
                ax.set_ylabel(label)
    axes[0, 0].legend(loc='best', frameon=False)
    fig.suptitle('PTC correction diagnostics | same query cells, exact k = 30 neighbors', fontsize=11)
    fig.text(.5, .017, 'Mixing is evaluated within archived population × MT/N/TU/T strata; 1 is exchangeable-batch expectation.\n'
             'Archived labels are concordance references. Higher mixing alone is not a biological-quality ranking.',
             ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .065, 1, .965))
    for ext in ['png', 'pdf']:
        fig.savefig(DEST / ('PTC_batch_biology_comparison.' + ext), dpi=250, bbox_inches='tight')
    plt.close(fig)
    interpretation = '''# PTC batch mixing and biological-label retention

This is a post-fit diagnostic on six existing correction controls and two saved spaces (PCA30 and UMAP2). No expression, features, clusters, seed annotations or DL outputs are refitted or modified.

Queries are selected with seed 42, up to 128 cells per archived broad population × sample, identically across the three correction arms. Each query uses exact Euclidean k=30 neighbors among all cells in its stated reference subset, excluding itself. Query coverage, represented-cell coverage, excluded populations and k-boundary ties are reported. Saved neighbor indices/distances permit reconstruction. Sampling is for diagnostic cost, not a new analysis cohort or a biological replicate.

**Mixing:** search only inside the same archived broad population and the same MT/N/TU/T sample prefix. Require at least two represented samples, at least ten cells in every represented sample, and more than 30 cells overall. Unknown/Other/ambiguous lymphoid/generic tumor populations are excluded from the primary mixing estimate and their counts remain visible. For query i, expected cross-sample fraction = (population-stratum size − own-sample size)/(population-stratum size − 1). Divide the observed fraction by this leave-self-out expectation. Zero means no cross-sample neighbors; one is the exchangeable sample-label expectation; values above one can occur and are not automatically better. Average sample means within each population/prefix, then average those strata equally. A separately reported cell-weighted estimate uses stratum-size/query-count weights. There is no patient confidence interval or significance claim.

**Biological-label retention:** search across all cells within the same sample-prefix stratum. Report same archived broad and native label fractions, plus chance-adjusted purity (observed − random-composition expectation)/(1 − expectation). Unknown/ambiguous cells remain possible neighbors but are excluded as primary queries. Archived annotations are concordance references inherited from the historical workflow, not independent truth; they can favor historical geometry. Negative adjusted purity is allowed. Broad and native label resolutions are separate.

**Tissue-state preservation:** for NMT, search within each fixed population across samples and report preservation of Normal adjacent versus Metastasis, with the same chance adjustment. TTU contains only Primary tumor, so this diagnostic is uninformative and remains missing. Tissue state and sample are confounded: tissue purity cannot establish removal of technical effects without loss of true state. For mixing, the separate T/TU prefixes are retained conservatively even though both are primary tumors; MT and N are never treated as biological states that should be indiscriminately mixed.

**TCR support:** within the label-retention graph, report neighboring TCR-detection fraction/enrichment for TCR-positive query cells, weighted back to sampled population/sample strata. These summaries are conditional on the archived populations included as queries. TCR-undetected is not established non-T truth; this is detection homophily, not T/non-T annotation accuracy.

RNA→Harmony uses fixed scoring/DL inputs and tests geometry correction. CCA→RNA changes expression correction, feature identity and geometry together and is a combined-workflow comparison. Distances are used only within each saved space; PCA and UMAP are not assumed to preserve the same metric. These diagnostics support a conditional trade-off assessment, not a global optimum or a ranking solely by the mixing ratio.
'''
    (DEST / 'INTERPRETATION.md').write_text(interpretation)
    report = dict(status='complete', expected_units=UNITS, n_units=len(manifests), n_unit_spaces=len(frame),
        same_queries_and_cell_order_across_correction_arms=True, k=K, seed=SEED, query_cap=QUERY_CAP,
        population_sample_balanced_summary=True, no_independent_truth_claim=True,
        job=os.environ['SLURM_JOB_ID'], source_sha256=sha(__file__), unit_manifest_sha256=inputs,
        outputs={p.name: sha(p) for p in DEST.iterdir() if p.is_file() and p.name not in ['manifest.json', 'COMPLETE']})
    (DEST / 'manifest.json').write_text(json.dumps(report, indent=2) + '\n')
    (DEST / 'COMPLETE').write_text(sha(DEST / 'manifest.json') + '\n')
    print(json.dumps(report, indent=2), flush=True)


def main():
    assert os.environ.get('SLURM_JOB_ID'), 'All diagnostic computation requires SLURM'
    parser = argparse.ArgumentParser()
    parser.add_argument('--unit', choices=UNITS)
    parser.add_argument('--pilot', action='store_true')
    parser.add_argument('--summarize', action='store_true')
    args = parser.parse_args()
    if args.summarize:
        summarize()
    else:
        unit = args.unit or UNITS[int(os.environ['SLURM_ARRAY_TASK_ID'])]
        audit_unit(unit, pilot=args.pilot)


if __name__ == '__main__':
    main()
