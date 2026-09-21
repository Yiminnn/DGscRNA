"""Rebuild frozen A2/C inference inputs; never refit or recompute inference."""
import numpy as np
import pandas as pd
from canonical_C import SPEC, MEASURES, ordered_patient_means, ordered_training_rank


def need(value, message):
    if not value:
        raise RuntimeError(message)


def align(old, fresh, keys):
    """Restore the original input row order before its specified reductions."""
    need(not old.duplicated(keys).any() and not fresh.duplicated(keys).any(), 'Duplicate conditions')
    want = old[keys].copy()
    got = want.merge(fresh, on=keys, how='left', validate='one_to_one', sort=False, indicator=True)
    need(got._merge.eq('both').all(), 'Fresh condition is missing')
    need(len(got) == len(old), 'Condition roster changed')
    return got.drop(columns='_merge')


def compare_rows(old, fresh, keys, fields, exact=False):
    fresh = align(old, fresh, keys)
    for name in ['patient', 'primary', 'n_cells', 'lfine_n_classes', 'truth_sha256']:
        if name in old and name in fresh:
            need(old[name].reset_index(drop=True).equals(fresh[name]), 'Endpoint metadata differs: ' + name)
    deviations = {}
    for name in fields:
        a, b = old[name].to_numpy(dtype=float), fresh[name].to_numpy(dtype=float)
        need(np.array_equal(np.isnan(a), np.isnan(b)), 'Missing metric differs: ' + name)
        finite = np.isfinite(a) & np.isfinite(b)
        deviations[name] = float(np.max(np.abs(a[finite] - b[finite]), initial=0))
        if exact:
            need(np.array_equal(a, b, equal_nan=True), 'Exact inference vector differs: ' + name)
        else:
            np.testing.assert_allclose(a, b, rtol=0, atol=1e-12, equal_nan=True)
    return fresh, deviations


def a2_reductions(metrics, reference, folds, budgets):
    """Original A2 group means, lambda ranking and patient joins; no new CI/p."""
    terminal = metrics[metrics.stage.eq('terminal090')].merge(
        folds[['sample', 'patient', 'fold']], on=['sample', 'patient'], validate='many_to_one')
    old = reference[reference.stage.eq('terminal090')]
    choices, pairs, selected_rows = [], [], []
    for cohort in ['primary', 'all']:
        subset = terminal[terminal.primary.eq(True)] if cohort == 'primary' else terminal
        baseline = old[old.primary.eq(True)] if cohort == 'primary' else old
        expected_samples = 97 if cohort == 'primary' else 121
        need(subset['sample'].nunique() == baseline['sample'].nunique() == expected_samples,
             'A2 cohort mismatch')
        for budget in budgets:
            data = subset[subset.budget == budget]
            pat = data.groupby(['patient', 'fold', 'lambda_value'])[['lfine_macroF1', 'coverage']].mean().reset_index()
            selected = []
            for fold in sorted(pat.fold.unique()):
                train = pat[pat.fold != fold]
                test = pat[pat.fold == fold]
                scores = train.groupby('lambda_value').lfine_macroF1.mean().reset_index(name='training_patient_lfine_macroF1')
                need(len(scores) == 5, 'A2 lambda grid incomplete')
                scores['distance_from_primary'] = (scores.lambda_value - 1).abs()
                rank = scores.sort_values(['training_patient_lfine_macroF1', 'distance_from_primary', 'lambda_value'],
                    ascending=[False, True, True], kind='stable')
                chosen = rank.iloc[0]
                choices.append(dict(cohort=cohort, budget=budget, fold=int(fold), lambda_value=float(chosen.lambda_value),
                    training_patient_lfine_macroF1=float(chosen.training_patient_lfine_macroF1),
                    n_training_patients=train.patient.nunique(), n_test_patients=test.patient.nunique(),
                    test_labels_used_for_selection=False))
                selected.append(data[(data.fold == fold) & (data.lambda_value == chosen.lambda_value)].assign(
                    selection='training_patient_selected', cohort=cohort))
            selections = [data[data.lambda_value == 1].assign(selection='fixed_lambda1', cohort=cohort), pd.concat(selected)]
            for chosen in selections:
                need(chosen['sample'].nunique() == len(chosen) == expected_samples, 'A2 selected sample mismatch')
                selected_rows.append(chosen)
                patient = chosen.groupby('patient')[['lfine_macroF1', 'coverage']].mean()
                mode = chosen.selection.iloc[0]
                for route, group in baseline[baseline.budget == budget].groupby('route'):
                    orig = group.groupby('patient')[['lfine_macroF1', 'coverage']].mean()
                    match = patient.join(orig, how='outer', lsuffix='_cellwise', rsuffix='_original')
                    need(not match.isna().any().any(), 'Incomplete A2 patient pair')
                    for patient_id, row in match.iterrows():
                        pairs.append(dict(cohort=cohort, budget=budget, selection=mode,
                            reference_route=route, patient=patient_id, **row.to_dict()))
    return pd.DataFrame(choices), pd.DataFrame(pairs), pd.concat(selected_rows, ignore_index=True)


def revalidate_a2(fresh, anchors, old_candidates, old_anchors, old_choices, old_pairs, old_stats, folds, out):
    need(len(fresh) == 2420 and len(anchors) == 1936, 'Incomplete fresh A2 endpoints')
    keys = ['sample', 'budget', 'lambda_value', 'stage']
    fresh, candidate_drift = compare_rows(old_candidates, fresh, keys, ['lfine_macroF1', 'coverage'])
    anchors, anchor_drift = compare_rows(old_anchors, anchors, ['sample', 'budget', 'route', 'stage'],
        ['lfine_macroF1', 'coverage'])
    choices, pairs, selected = a2_reductions(fresh, anchors, folds, ['hvg2000', 'hvg5000'])
    need(len(choices) == 20 and len(old_stats) == 32, 'A2 selection/inference roster mismatch')
    choices, _ = compare_rows(old_choices, choices, ['cohort', 'budget', 'fold'],
        ['lambda_value', 'training_patient_lfine_macroF1', 'n_training_patients', 'n_test_patients'], exact=True)
    pairkeys = ['cohort', 'budget', 'selection', 'reference_route', 'patient']
    pairfields = ['lfine_macroF1_cellwise', 'coverage_cellwise', 'lfine_macroF1_original', 'coverage_original']
    pairs, pair_drift = compare_rows(old_pairs, pairs, pairkeys, pairfields, exact=True)
    # These exact float64 vectors are the original inferential inputs. No rounding
    # or relaxed p-value tolerance can authorize retention of old inference.
    np.savez(out / 'A2_exact_patient_vectors.npz', **{field: pairs[field].to_numpy() for field in pairfields})
    pairs.to_csv(out / 'A2_revalidated_patient_pairs.csv', index=False)
    choices.to_csv(out / 'A2_revalidated_lambda_choices.csv', index=False)
    selected.to_csv(out / 'A2_fresh_selected_samples.csv.gz', index=False)
    old_stats.to_csv(out / 'A2_retained_original_inference.csv', index=False)
    fresh.to_csv(out / 'A2_fresh_candidate_metrics.csv.gz', index=False)
    state = selected.groupby(['cohort', 'budget', 'selection', 'dl_status', 'training_executed']).size().rename('samples').reset_index()
    state.to_csv(out / 'A2_fresh_execution_states.csv', index=False)
    return dict(status='exact_patient_inputs_and_choices', candidate_rows=len(fresh), anchor_rows=len(anchors),
        choice_rows=len(choices), patient_pair_rows=len(pairs), retained_contrasts=len(old_stats),
        sample_metric_max_abs_delta={'candidate': candidate_drift, 'anchor': anchor_drift},
        patient_vector_max_abs_delta=pair_drift, inference_recomputed=False), choices


def revalidate_c(fresh, old_candidates, old_choices, old_heldout, old_summary, old_contrasts, folds, out):
    # The compact comparison is the frozen original partition (39 DG candidates).
    # The separate equal24 search and every non-DG prediction remain archived.
    old = old_candidates[old_candidates.method.eq('DG-scRNA') & old_candidates.primary_candidate].copy()
    need(len(old) == 121 * 39, 'Original C candidate roster changed')
    fresh = fresh.assign(method='DG-scRNA')
    fresh, sample_drift = compare_rows(old, fresh, ['sample'] + SPEC, MEASURES)
    choices, heldout = [], []
    for cohort, data in [('primary97', fresh[fresh.primary]), ('all121', fresh)]:
        patient = ordered_patient_means(data, folds[['patient', 'fold']].drop_duplicates())
        need(patient.patient.nunique() == (55 if cohort == 'primary97' else 59), 'C patient roster changed')
        for fold in sorted(patient.fold.unique()):
            rank = ordered_training_rank(patient, fold)
            need(len(rank) == 39 and rank.lfine_macroF1.notna().all(), 'C training roster changed')
            best = rank.iloc[0]
            chosen = patient[patient.fold.eq(fold)]
            for field in SPEC:
                chosen = chosen[chosen[field].eq(best[field])]
            need(chosen.patient.is_unique, 'Duplicate held-out C patient')
            choices.append(dict(cohort=cohort, fold=int(fold), training_patient_lfine_macroF1=float(best.lfine_macroF1),
                **{field: best[field] for field in SPEC}))
            for row in chosen.to_dict('records'):
                heldout.append(dict(cohort=cohort, selection='fixed_partition', **row))
    choices, heldout = pd.DataFrame(choices), pd.DataFrame(heldout)
    previous = old_choices[old_choices.method.eq('DG-scRNA')]
    choices = align(previous, choices, ['cohort', 'fold', 'method'])
    for field in SPEC + ['training_patient_lfine_macroF1']:
        need(np.array_equal(previous[field].to_numpy(), choices[field].to_numpy()), 'C selected configuration/score differs: ' + field)
    reference = old_heldout[old_heldout.method.eq('DG-scRNA')]
    heldout, vector_drift = compare_rows(reference, heldout, ['cohort', 'method', 'patient'], MEASURES, exact=True)
    np.savez(out / 'C_exact_DG_patient_vectors.npz', **{field: heldout[field].to_numpy() for field in MEASURES})
    choices.to_csv(out / 'C_revalidated_DG_fold_choices.csv', index=False)
    heldout.to_csv(out / 'C_fresh_DG_heldout_patients.csv', index=False)
    fresh.to_csv(out / 'C_fresh_DG_candidate_metrics.csv.gz', index=False)
    rows = []
    for method in ['DG-scRNA', 'scType', 'scCATCH', 'SCINA', 'SingleR', 'scDeepSort']:
        data = heldout if method == 'DG-scRNA' else old_heldout
        data = data[data.cohort.eq('primary97') & data.method.eq(method)]
        need(len(data) == 55 and data.patient.is_unique, 'C compact patient roster changed')
        previous = old_summary[old_summary.cohort.eq('primary97') & old_summary.method.eq(method)]
        need(len(previous) == 1, 'C original summary mismatch')
        row = dict(method=method, n_patients=55, input_origin='fresh packaged terminal' if method == 'DG-scRNA' else 'archived fixed predictions')
        for field in ['lfine_macroF1', 'coverage', 'offvocab_rate']:
            row[field] = float(data[field].mean())
            np.testing.assert_allclose(row[field], previous.iloc[0][field + '_mean'], rtol=0, atol=1e-12)
        if method != 'DG-scRNA':
            contrast = old_contrasts[old_contrasts.cohort.eq('primary97') & old_contrasts.method.eq(method)]
            need(len(contrast) == 1, 'C contrast missing')
            row.update({field: float(contrast.iloc[0][field]) for field in ['mean_delta', 'CI95_low', 'CI95_high', 'p_Holm']})
        rows.append(row)
    pd.DataFrame(rows).to_csv(out / 'C_current_DG_archived_comparators.csv', index=False)
    return dict(status='exact_DG_patient_inputs_and_choices', fresh_DG_candidate_rows=len(fresh),
        choice_rows=len(choices), exact_DG_patient_rows=len(heldout), non_DG_refitted=False,
        sample_metric_max_abs_delta=sample_drift, patient_vector_max_abs_delta=vector_drift,
        inference_recomputed=False, equal24_search='archived separate analysis; not refreshed here'), rows
