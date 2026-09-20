"""Frozen patient-fold K selection after all A1 predictions/figures are complete."""
from pathlib import Path
import json
import os
import sys
CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE / 'legacy'))
from common import ROOT, OUT as REFERENCE, require_slurm, sha, checked, complete, write_json, utc
CAMPAIGN = ROOT / 'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
OUT = CAMPAIGN / 'embedding'


def verify_evaluation(directory, output_names=None):
    """Reject stale/corrupt consumed outputs even when COMPLETE itself matches."""
    directory = Path(directory)
    assert checked(directory), f'Incomplete evaluation manifest: {directory}'
    manifest = json.loads((directory / 'manifest.json').read_text())
    names = manifest['outputs'] if output_names is None else output_names
    for name in names:
        assert name in manifest['outputs'], f'Unregistered evaluation output: {directory / name}'
        assert sha(directory / name) == manifest['outputs'][name], f'Changed evaluation output: {directory / name}'
    return manifest


def assert_primary_boolean(frame):
    import pandas as pd
    assert pd.api.types.is_bool_dtype(frame['primary'].dtype), 'Primary cohort flags must be parsed booleans'
    assert frame['primary'].isin([True, False]).all(), 'Primary cohort flags must not contain missing values'


def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    from scipy.stats import wilcoxon
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    spec = json.loads((CAMPAIGN / 'protocol/embedding.json').read_text())
    choice_protocol = json.loads((CAMPAIGN / 'protocol/embedding_selection.json').read_text())
    foldfile = CAMPAIGN / 'protocol/embedding_patient_folds.csv'
    assert sha(foldfile) == choice_protocol['patient_folds_sha256']
    folds = pd.read_csv(foldfile)[['patient', 'fold']].drop_duplicates()
    assert not folds.patient.duplicated().any()
    frames, partitions, sources, missing = [], [], {}, []
    for sample in spec['samples']:
        for budget in spec['budgets']:
            for space in spec['spaces']:
                directory = OUT / sample / budget / space
                if not checked(directory, 'fit_manifest.json', 'FIT_COMPLETE') or not checked(directory / 'evaluation') or not checked(directory / 'figures'):
                    missing.append(f'{sample}/{budget}/{space}')
                    continue
                cfg = json.loads((directory / 'config.json').read_text())
                verify_evaluation(directory / 'evaluation')
                context = pd.DataFrame([dict(route=Path(c['dest']).name, method=c['method'], k=c['k'] if c['k'] is not None else 0) for c in cfg['conditions']])
                m = pd.read_csv(directory / 'evaluation/metrics.csv', dtype={'cutoff': str})
                assert_primary_boolean(m)
                assert len(m) == 39 and m.library.eq('CM2_glioma_other').all() and m.cutoff.eq('mean').all()
                m = m.merge(context, on='route', validate='many_to_one').assign(space=space)
                frames.append(m)
                cp = pd.read_csv(directory / 'evaluation/clustering.csv').merge(context, on='route', validate='one_to_one').assign(space=space)
                assert len(cp) == 13
                partitions.append(cp)
                sources[str(directory / 'evaluation/manifest.json')] = sha(directory / 'evaluation/manifest.json')
    summary = OUT / 'summary';summary.mkdir(exist_ok=True)
    if missing:
        write_json(summary / 'INCOMPLETE.json', dict(missing=missing, complete_representations=len(frames), expected=1694, time=utc()))
        raise RuntimeError(f'{len(missing)} representations remain incomplete; no partial ranking released')
    all_metrics = pd.concat(frames, ignore_index=True)
    assert_primary_boolean(all_metrics)
    assert len(all_metrics) == 3 * 22022
    assert not all_metrics.duplicated(['sample', 'budget', 'space', 'method', 'k', 'stage']).any()
    all_metrics.to_csv(summary / 'all_candidate_metrics.csv.gz', index=False, compression='gzip')
    pd.concat(partitions, ignore_index=True).to_csv(summary / 'all_candidate_clustering.csv.gz', index=False, compression='gzip')
    measures = ['macroF1_present', 'macroF1_fixed11', 'weightedF1', 'accuracy', 'coverage', 'unknown_rate', 'mapped_coverage', 'off_vocabulary_rate', 'coarse10_macroF1_present']
    selected_rows, choices, patient_rows, comparison_rows = [], [], [], []
    rng = np.random.default_rng(20260920)
    for cohort, input_rows in [('primary97', all_metrics[all_metrics.primary]), ('all121', all_metrics)]:
        terminal = input_rows[input_rows.stage.eq('terminal090')]
        patients = terminal.groupby(['patient', 'budget', 'space', 'method', 'k'])[measures].mean().reset_index()
        patients = patients.merge(folds, on='patient', validate='many_to_one')
        assert patients.fold.notna().all()
        chosen_samples = []
        for (budget, space, method), group in patients.groupby(['budget', 'space', 'method']):
            for fold in sorted(group.fold.unique()):
                train, test = group[group.fold != fold], group[group.fold == fold]
                ranking = train.groupby('k').macroF1_present.mean().reset_index().sort_values(['macroF1_present', 'k'], ascending=[False, True], kind='stable')
                k = int(ranking.iloc[0].k)
                if method == 'HDBSCAN_R':
                    assert k == 0 and len(ranking) == 1
                else:
                    assert set(ranking.k) == set(spec['K'])
                held = test[test.k == k]
                assert held.patient.nunique() == test.patient.nunique()
                choices.append(dict(cohort=cohort, budget=budget, space=space, method=method, fold=int(fold), selected_k=k,
                                    training_patient_mean_macroF1=float(ranking.iloc[0].macroF1_present),
                                    n_training_patients=train.patient.nunique(), n_test_patients=test.patient.nunique()))
                selected = input_rows[(input_rows.budget == budget) & (input_rows.space == space) & (input_rows.method == method) & (input_rows.k == k) & input_rows.patient.isin(held.patient)].copy()
                selected['fold'] = int(fold);selected['cohort'] = cohort
                chosen_samples.append(selected)
        selected = pd.concat(chosen_samples, ignore_index=True)
        assert not selected.duplicated(['sample', 'budget', 'space', 'method', 'stage']).any()
        expected = (97 if cohort == 'primary97' else 121) * 2 * 7 * 3
        assert selected.stage.eq('terminal090').sum() == expected
        selected_rows.append(selected)
        patient = selected.groupby(['cohort', 'patient', 'budget', 'space', 'method', 'stage'])[measures].mean().reset_index()
        patient_rows.append(patient)
        for budget in spec['budgets']:
            anchors = []
            for sample in spec['samples']:
                anchor_path = REFERENCE / 'GBM' / sample / budget / 'evaluation/metrics.csv'
                verify_evaluation(anchor_path.parent, ['metrics.csv'])
                sources[str(anchor_path)] = sha(anchor_path)
                sources[str(anchor_path.parent / 'manifest.json')] = sha(anchor_path.parent / 'manifest.json')
                old = pd.read_csv(anchor_path, dtype={'cutoff': str})
                assert_primary_boolean(old)
                old = old[(old.route == 'UMAP2_HDBSCAN_R') & (old.library == 'CM2_glioma_other') & (old.cutoff == 'mean') & (old.stage == 'terminal090') & (old.family == 'native_R_budget')]
                if cohort == 'primary97':
                    old = old[old.primary]
                anchors.append(old)
            anchor = pd.concat(anchors).groupby('patient')[['macroF1_present', 'coverage']].mean().rename(columns={'macroF1_present': 'anchor_F1', 'coverage': 'anchor_coverage'}).reset_index()
            current = patient[(patient.budget == budget) & (patient.stage == 'terminal090')]
            for (space, method), group in current.groupby(['space', 'method']):
                joined = group.merge(anchor, on='patient', validate='one_to_one')
                assert len(joined) == len(anchor)
                delta = (joined.macroF1_present - joined.anchor_F1).to_numpy()
                bootstrap = delta[rng.integers(0, len(delta), size=(10000, len(delta)))].mean(axis=1)
                comparison_rows.append(dict(cohort=cohort, budget=budget, space=space, method=method, n_patients=len(delta),
                    candidate_mean_F1=float(joined.macroF1_present.mean()), anchor_mean_F1=float(joined.anchor_F1.mean()),
                    mean_delta=float(delta.mean()), CI95_low=float(np.quantile(bootstrap, .025)), CI95_high=float(np.quantile(bootstrap, .975)),
                    p_wilcoxon=float(wilcoxon(delta).pvalue) if np.any(np.abs(delta) > 1e-14) else 1.,
                    candidate_coverage=float(joined.coverage.mean()), anchor_coverage=float(joined.anchor_coverage.mean()),
                    reference='Same-budget original native R PCA30->UMAP2/HDBSCAN; CM2_glioma_other/mean; terminal090',
                    interpretation='Patient-fold selected K; retrospective conditional-on-selected-predictions paired inference'))
    selected = pd.concat(selected_rows, ignore_index=True)
    selected.to_csv(summary / 'selected_sample_metrics.csv.gz', index=False, compression='gzip')
    pd.DataFrame(choices).to_csv(summary / 'patient_fold_K_choices.csv', index=False)
    pd.concat(patient_rows, ignore_index=True).to_csv(summary / 'patient_metrics.csv.gz', index=False, compression='gzip')
    stats = pd.DataFrame(comparison_rows)
    stats['p_Holm'] = 1.
    for _, group in stats.groupby(['cohort', 'budget']):
        assert len(group) == 21
        ordered = group.sort_values('p_wilcoxon', kind='stable').index
        stats.loc[ordered, 'p_Holm'] = np.minimum(1., np.maximum.accumulate(stats.loc[ordered, 'p_wilcoxon'].to_numpy() * np.arange(len(ordered), 0, -1)))
    stats.to_csv(summary / 'paired_vs_original_anchor.csv', index=False)
    fig, axes = plt.subplots(1, 2, figsize=(12, 9), layout='constrained', sharey=True)
    ordering = [(space, method) for space in spec['spaces'] for method in spec['clusterers']]
    for ax, budget in zip(axes, spec['budgets']):
        table = stats[(stats.cohort == 'primary97') & (stats.budget == budget)].set_index(['space', 'method']).loc[ordering]
        y = np.arange(len(table));values = table.mean_delta.to_numpy()
        ax.errorbar(values, y, xerr=np.vstack([values-table.CI95_low, table.CI95_high-values]), fmt='o', markersize=4, capsize=2, color='#27647b')
        ax.axvline(0, color='#888', linewidth=.8);ax.set_yticks(y, [f'{space} / {method}' for space, method in ordering], fontsize=8)
        ax.set_title(budget + ': 55 primary-cohort patients');ax.set_xlabel('Mean macro-F1 difference from original R anchor (95% CI)')
        ax.grid(axis='x', color='#eeeeee');ax.spines[['top', 'right']].set_visible(False)
    axes[0].invert_yaxis()
    fig.suptitle('Controlled representations and clustering: patient-fold K selection\nFixed marker context; all cells retained; all contrasts shown', fontsize=12)
    for suffix in ['png', 'pdf']:
        fig.savefig(summary / f'paired_embedding_comparison.{suffix}', dpi=180, bbox_inches='tight')
    plt.close(fig)
    write_json(summary / 'manifest.json', dict(status='completed', expected_candidate_partitions=22022,
        selected_all121_conditions=5082, n_selection_rows=len(choices), n_patient_contrasts=len(stats),
        protocol_sha256=sha(CAMPAIGN / 'protocol/embedding_selection.json'), source_sha256=sha(__file__),
        input_manifests=sources, job=os.environ['SLURM_JOB_ID'], completed_at=utc(),
        files={p.name: sha(p) for p in summary.iterdir() if p.suffix in ['.csv', '.gz', '.png', '.pdf']}))
    complete(summary)


if __name__ == '__main__':
    run()
