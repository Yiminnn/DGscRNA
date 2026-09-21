"""Verify tied-rank inference on losslessly recorded producer patient vectors.

The original independent groupby reduction remains checked at its existing
precision. Wilcoxon is independently recomputed on the exact producer vectors:
machine-level differences between groupby.mean and Series.mean alter ties.
No producer arithmetic, p value, rank, selection rule or tolerance is changed.
"""
from pathlib import Path
import hashlib
import json


def verify(numerical, actual, independent, anchors, spec, out):
    import numpy as np
    import pandas as pd
    from scipy.stats import wilcoxon
    choices, selected, patients, statistics = actual
    reference_choices, reference_selected, reference_patients, reference_statistics = independent
    keys = [['cohort', 'budget', 'space', 'method', 'fold'],
        ['cohort', 'sample', 'budget', 'space', 'method', 'stage'],
        ['cohort', 'patient', 'budget', 'space', 'method', 'stage']]
    for left, right, key in zip(actual[:3], independent[:3], keys):
        numerical.equal(left, right, key)
    comparison_keys = ['cohort', 'budget', 'space', 'method']
    non_p = [name for name in statistics.columns if name not in comparison_keys + ['p_wilcoxon', 'p_Holm']]
    numerical.equal(statistics, reference_statistics, comparison_keys, non_p)
    out = Path(out)
    assert not out.exists(), 'Preserve inference verification attempts'
    out.mkdir(parents=True)
    arrays, entries, expected = {}, [], []
    for cohort in ['primary97', 'all121']:
        for budget in spec['budgets']:
            raw = selected[selected.cohort.eq(cohort) & selected.budget.eq(budget) & selected.stage.eq('terminal090')]
            current = patients[patients.cohort.eq(cohort) & patients.budget.eq(budget) & patients.stage.eq('terminal090')]
            anchor_samples = anchors[anchors.budget.eq(budget) & anchors.stage.eq('terminal090') &
                (anchors.primary if cohort == 'primary97' else True)]
            # This is the original anchor reduction contract, rebuilt from rows.
            anchor_patient = anchor_samples.groupby('patient')[['lfine_macroF1', 'coverage']].mean()
            for space, method in sorted(set(zip(current.space, current.method))):
                canonical = current[current.space.eq(space) & current.method.eq(method)].set_index('patient').sort_index()
                raw_group = raw[raw.space.eq(space) & raw.method.eq(method)]
                assert set(canonical.index) == set(anchor_patient.index) == set(raw_group.patient)
                # Rebuild producer candidate means independently from selected
                # sample rows, using exactly its Series.mean reduction contract.
                rebuilt = np.asarray([raw_group.loc[raw_group.patient.eq(patient), 'lfine_macroF1'].mean()
                    for patient in canonical.index], dtype=np.float64)
                np.testing.assert_array_equal(rebuilt, canonical.lfine_macroF1.to_numpy(dtype=np.float64))
                reference = anchor_patient.loc[canonical.index]
                valid = canonical.lfine_macroF1.notna()
                assert valid.equals(reference.lfine_macroF1.notna())
                ident = f'contrast_{len(entries):03d}'
                arrays[ident + '_patients'] = canonical.index[valid].to_numpy(dtype=str)
                arrays[ident + '_candidate'] = rebuilt[valid.to_numpy()]
                arrays[ident + '_anchor'] = reference.loc[valid, 'lfine_macroF1'].to_numpy(dtype=np.float64)
                entries.append(dict(key=ident, cohort=cohort, budget=budget, space=space, method=method,
                    n_patients_total=len(canonical), n_patients_Lfine=int(valid.sum())))
    assert len(entries) == 84
    path = out / 'canonical_patient_vectors.npz'
    np.savez_compressed(path, **arrays)
    # Actual inference reads the saved float64 artifact, not an approximate CSV
    # or the producer's p values. NPZ round-trip must preserve every bit/value.
    with np.load(path, allow_pickle=False) as saved:
        for key, value in arrays.items():
            np.testing.assert_array_equal(saved[key], value)
        for entry in entries:
            ident = entry['key']
            delta = saved[ident + '_candidate'] - saved[ident + '_anchor']
            p = wilcoxon(delta).pvalue if (np.abs(delta) > 1e-14).any() else 1.
            expected.append({**{name: entry[name] for name in comparison_keys}, 'p_wilcoxon': float(p)})
    expected = pd.DataFrame(expected)
    expected['p_Holm'] = 1.
    for _, group in expected.groupby(['cohort', 'budget']):
        assert len(group) == 21
        running = 0.
        for rank, index in enumerate(sorted(group.index, key=lambda i: expected.at[i, 'p_wilcoxon'])):
            running = max(running, (21 - rank) * expected.at[index, 'p_wilcoxon'])
            expected.at[index, 'p_Holm'] = min(1., running)
    numerical.equal(statistics, expected, comparison_keys, ['p_wilcoxon', 'p_Holm'])
    left = statistics.set_index(comparison_keys).sort_index()
    right = reference_statistics.set_index(comparison_keys).sort_index()
    observed = {name: float(np.abs(left[name] - right[name]).max()) for name in ['p_wilcoxon', 'p_Holm']}
    expected.to_csv(out / 'independent_canonical_p_values.csv', index=False)
    (out / 'contrast_order.json').write_text(json.dumps(entries, indent=2) + '\n')
    report = dict(status='passed_canonical_inference', contrasts=84,
        independent_selection_patient_means_and_non_p_statistics_atol=1e-12,
        producer_patient_means_exactly_rebuilt_from_selected_sample_rows=True,
        canonical_float64_vectors_losslessly_saved=True,
        Wilcoxon_and_Holm_independently_recomputed=True,
        canonical_vs_producer_max_abs={name: float(np.abs(left[name] - expected.set_index(comparison_keys).sort_index()[name]).max())
            for name in ['p_wilcoxon', 'p_Holm']},
        original_independent_reduction_p_max_abs=observed,
        reason='Original producer uses Series.mean for patient groups; independent groupby.mean may round differently and change Wilcoxon ties. Independent means and all non-p statistics retain the original tolerance; inference uses lossless canonical producer vectors.',
        production_arithmetic_changed=False, ranks_rounded=False, tolerance_relaxed=False,
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        outputs={p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()})
    (out / 'manifest.json').write_text(json.dumps(report, indent=2) + '\n')
    return report
