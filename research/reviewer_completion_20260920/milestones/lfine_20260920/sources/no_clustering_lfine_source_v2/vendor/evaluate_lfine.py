#!/usr/bin/env python3
"""Re-score frozen GBM terminal calls using v5's Lfine set-valued rule.

No expression matrices, model fits, or changes to archived results. SLURM only.
The existing native-R semantic map is frozen, including the Non-neuron correction.
"""
from pathlib import Path
import argparse
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
import os
import sys
import time

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
NATIVE = ROOT/'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917'
OUT = ROOT/'results/hvg_ptc_20260916_v1/lfine_compact_20260920'
BUDGETS = ['all', 'hvg500', 'hvg1000', 'hvg2000', 'hvg3000', 'hvg5000']
ROUTES = ['PCA30_SNN', 'PCA30_HDBSCAN_R', 'UMAP2_SNN', 'UMAP2_HDBSCAN_R']
FIXED = 'CM2_glioma_other'
ABSTAIN = {'Unknown', 'Undecided', 'Noise', 'nan', 'None', ''}
SCORE_FIELDS = ['lfine_macroF1', 'coverage', 'legacy_called_coverage', 'abstain_rate',
                'offvocab_rate', 'acc_on_called', 'n_distinct_calls', 'n_classes_hit',
                'lfine_n_classes', 'lfine_scored_class_cell_fraction',
                'marker_vocab_oracle_upper_bound', 'n_reference_lfine_disagreements']
SOURCE_PATHS = [Path(__file__), ROOT/'handoff/g274/grid_sets.py',
                ROOT/'handoff/g274_table4/v5_final_annotations.py',
                ROOT/'handoff/markers_v3/mapping_L1_v3.csv',
                NATIVE/'markers/panel_L1_mapping.csv', NATIVE/'markers/manifest.json',
                NATIVE/'protocol/cohort.csv']


def sha(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def initialize():
    global np, pd, v5, helpers, mappings, libraries
    import numpy as np
    import pandas as pd
    sys.path.insert(0, str(ROOT/'handoff/g274_table4'))
    import v5_final_annotations as v5
    helpers = v5.legacy_helpers()
    mapping = pd.read_csv(NATIVE/'markers/panel_L1_mapping.csv', dtype=str, keep_default_na=False)
    mappings = {lib: dict(zip(g.panel, g.L1)) for lib, g in mapping.groupby('library')}
    libraries = json.loads((NATIVE/'markers/manifest.json').read_text())['libraries']
    helpers['panel_labels'] = {lib: set(mappings[lib].values()) for lib in libraries}


def audit_mapping():
    current = pd.read_csv(NATIVE/'markers/panel_L1_mapping.csv', dtype=str, keep_default_na=False)
    old = pd.read_csv(ROOT/'handoff/markers_v3/mapping_L1_v3.csv', dtype=str, keep_default_na=False)
    audit = current.merge(old[['marker_set', 'panel', 'gold_class']],
                          left_on=['library', 'panel'], right_on=['marker_set', 'panel'],
                          how='left', validate='one_to_one')
    audit['old_v5_panel_available'] = audit.gold_class.notna()
    audit['exact_old_v5_parent'] = audit.old_v5_panel_available & audit.L1.eq(audit.gold_class)
    audit['has_lfine_targets'] = audit.L1.isin(helpers['LFINE_PREFIX'])
    audit['parent_used'] = audit.L1
    audit['map_policy'] = np.where(audit.old_v5_panel_available,
        np.where(audit.exact_old_v5_parent, 'same_as_v5', 'frozen_native_R_negation_correction'),
        'existing_native_R_prespecified_semantic_mapping')
    changed = audit[audit.old_v5_panel_available & ~audit.exact_old_v5_parent]
    assert len(changed) == 3
    assert changed.panel.str.contains('Non-neuron', regex=False).all()
    assert changed.L1.eq('UNMAPPABLE').all() and changed.gold_class.eq('AMBIGUOUS_NEURON').all()
    assert audit[audit.library.eq(FIXED)].exact_old_v5_parent.all()
    audit.to_csv(OUT/'mapping_audit.csv.gz', index=False)
    summary = audit.groupby('library', sort=False).agg(
        n_panels=('panel', 'size'), n_old_v5_panels=('old_v5_panel_available', 'sum'),
        n_exact_old_v5_parent=('exact_old_v5_parent', 'sum'),
        n_panels_with_lfine_targets=('has_lfine_targets', 'sum')).reset_index()
    summary.to_csv(OUT/'mapping_summary.csv', index=False)
    (OUT/'lfine_target_prefixes.json').write_text(json.dumps(helpers['LFINE_PREFIX'], indent=2)+'\n')
    return dict(original_library_count=8, added_library_count=8,
                original_panel_changes=changed[['library', 'panel', 'gold_class', 'L1']].to_dict('records'),
                fixed_marker_identical_to_v5=True, reference_outcomes_used_for_mapping=False)


def sample_run(cohort_row):
    sample, patient, primary = cohort_row['sample'], cohort_row['patient'], bool(cohort_row['primary'])
    truth_path = NATIVE/'evaluation_inputs'/sample/'truth.csv.gz'
    im = json.loads((NATIVE/'inputs'/sample/'input_manifest.json').read_text())
    truth_hash = sha(truth_path)
    assert truth_hash == im['evaluation_files']['truth.csv.gz']
    truth = pd.read_csv(truth_path, dtype=str, keep_default_na=False).rename(
        columns={'cell_id': 'CellID', 'lfine_original': 'Lfine'})
    assert truth.CellID.is_unique and truth.Patient.eq(patient).all()
    scope = v5.label_scope(truth, helpers)
    assert np.array_equal(scope[0], truth.Lfine.to_numpy())
    tasks = [(b, r, FIXED) for b in BUDGETS for r in ROUTES]
    tasks += [('all', 'UMAP2_HDBSCAN_R', lib) for lib in libraries if lib != FIXED]
    assert len(tasks) == 39 and len(set(tasks)) == 39
    rows, per_class = [], []
    score_manifests = {}
    for budget, route, library in tasks:
        src = NATIVE/'GBM'/sample/budget/route
        if (budget, route) not in score_manifests:
            sp = src/'score_manifest.json'
            sh = sha(sp)
            assert (src/'SCORE_COMPLETE').read_text().strip() == sh
            score_manifests[(budget, route)] = (json.loads(sp.read_text()), sh)
        score, score_hash = score_manifests[(budget, route)]
        arm_ids = [aid for aid, a in score['arms'].items()
                   if a['library'] == library and str(a['cutoff']) == 'mean']
        assert len(arm_ids) == 1
        aid = arm_ids[0]
        terminal = src/'terminal'/aid
        mp = terminal/'terminal_manifest.json'
        pp = terminal/'predictions.csv.gz'
        row = dict(sample=sample, patient=patient, primary=primary, budget=budget, route=route,
                   library=library, cutoff='mean', stage='terminal090', arm_id=aid,
                   family='native_R_budget', n_cells=len(truth),
                   lfine_n_classes=len(scope[1]), status='unavailable', terminal_valid=False,
                   score_manifest_sha256=score_hash, truth_sha256=truth_hash,
                   predictions_path=str(pp.relative_to(ROOT)),
                   terminal_manifest_path=str(mp.relative_to(ROOT)),
                   core_hvg24=(library == FIXED),
                   fullgene_marker16=(budget == 'all' and route == 'UMAP2_HDBSCAN_R'))
        if not mp.exists() or not (terminal/'TERMINAL_COMPLETE').exists():
            row['unavailable_reason'] = 'missing_terminal_manifest_or_completion'
            rows.append(row)
            continue
        mh = sha(mp)
        tm = json.loads(mp.read_text())
        assert (terminal/'TERMINAL_COMPLETE').read_text().strip() == mh
        assert tm['score_manifest_sha256'] == score_hash
        assert tm['arm']['library'] == library and tm['arm']['cutoff'] == 'mean'
        row.update(terminal_manifest_sha256=mh, dl_status=tm['dl_status'],
                   training_executed=tm['training_executed'], DL_features=tm['DL_features'],
                   n_training_classes=tm['n_training_classes'], n_pool=tm['n_pool'])
        if not tm['terminal_valid'] or tm['status'] != 'completed':
            row['unavailable_reason'] = 'invalid_terminal_state'
            rows.append(row)
            continue
        ph = sha(pp)
        assert ph == tm['predictions_sha256']
        pred = pd.read_csv(pp, dtype=str, keep_default_na=False,
                           usecols=['cell_id', 'initial', 'final090'])
        assert np.array_equal(pred.cell_id, truth.CellID)
        known = pred.initial.ne('Undecided')
        assert np.array_equal(pred.loc[known, 'initial'], pred.loc[known, 'final090'])
        mapped = np.asarray([mappings[library].get(v, 'Unknown' if v in ABSTAIN else 'UNMAPPABLE')
                             for v in pred.final090], dtype=object)
        values = v5.metrics(mapped, truth, scope, helpers, library)
        row.update({k: values[k] for k in SCORE_FIELDS})
        assert row['n_reference_lfine_disagreements'] == 0
        row.update(status='completed', terminal_valid=True, predictions_sha256=ph,
                   source_native_labels_preserved=True)
        rows.append(row)
        if sample == 'TKU3186':
            lf, classes, targets = scope
            ok = np.asarray([g in targets.get(q, ()) for g, q in zip(lf, mapped)])
            for c in classes:
                g = lf == c
                has = np.asarray([c in targets.get(q, ()) for q in mapped])
                tp, fn, fp = int((g & ok).sum()), int((g & ~ok).sum()), int((~ok & has & ~g).sum())
                per_class.append(dict(sample=sample, budget=budget, route=route, library=library,
                    cutoff='mean', stage='terminal090', lfine_class=c, support=int(g.sum()),
                    TP=tp, FN=fn, FP=fp, F1=2*tp/(2*tp+fn+fp) if 2*tp+fn+fp else 0.0))
    pd.DataFrame(rows).to_csv(OUT/'sample_metrics'/f'{sample}.csv', index=False)
    if per_class:
        pd.DataFrame(per_class).to_csv(OUT/'TKU3186_perclass.csv', index=False)
    return dict(sample=sample, n_conditions=len(rows), n_valid=sum(r['terminal_valid'] for r in rows),
                n_cells=len(truth), n_scored_classes=len(scope[1]), compose_mismatches=0)


def summarize(df, keys):
    rows = []
    scopes = {'primary97': df[df.primary], 'all121': df, 'sensitivity24': df[~df.primary],
              'TKU3186': df[df['sample'].eq('TKU3186')]}
    metrics = ['lfine_macroF1', 'coverage', 'legacy_called_coverage', 'offvocab_rate',
               'marker_vocab_oracle_upper_bound', 'lfine_scored_class_cell_fraction']
    for scope, data in scopes.items():
        for labels, g in data.groupby(keys, sort=False):
            if not isinstance(labels, tuple):
                labels = (labels,)
            row = dict(scope=scope, **dict(zip(keys, labels)), n_samples_expected=g['sample'].nunique(),
                       n_patients_expected=g.patient.nunique(), n_samples_valid=int(g.terminal_valid.sum()),
                       n_samples_unavailable=int((~g.terminal_valid).sum()),
                       dl_status_counts=json.dumps(g.dl_status.fillna('not_recorded').value_counts().to_dict(), sort_keys=True),
                       n_training_executed=int(g.training_executed.fillna(False).sum()))
            for metric in metrics:
                valid = g.loc[g.terminal_valid & g[metric].notna()]
                patient = valid.groupby('patient')[metric].mean()
                row.update({metric+'_patient_mean': float(patient.mean()) if len(patient) else np.nan,
                            metric+'_sample_mean': float(valid[metric].mean()) if len(valid) else np.nan,
                            metric+'_n_samples': len(valid), metric+'_n_patients': len(patient)})
            rows.append(row)
    return pd.DataFrame(rows)


def verify_historical():
    """Exact historical saved endpoint parity, not a new model run."""
    pc = pd.read_csv(ROOT/'results/cohort_v4/TKU3186/percell.csv', dtype=str, keep_default_na=False)
    saved = pd.read_csv(ROOT/'results/g274_v5/final_annotations/TKU3186_final_24.csv', dtype={'cutoff': str})
    scope = v5.label_scope(pc, helpers)
    checks = []
    old_map = pd.read_csv(ROOT/'handoff/markers_v3/mapping_L1_v3.csv', keep_default_na=False)
    old_helpers = dict(helpers)
    old_helpers['panel_labels'] = {m: set(g.gold_class) for m, g in old_map.groupby('marker_set')}
    for row in saved.itertuples():
        got = v5.metrics(pc[row.prediction_column].to_numpy(), pc, scope, old_helpers, row.marker_set)
        assert np.isclose(got['lfine_macroF1'], row.lfine_macroF1, rtol=0, atol=1e-12)
        checks.append(dict(marker=row.marker_set, cutoff=row.cutoff, difference=got['lfine_macroF1']-row.lfine_macroF1))
    pd.DataFrame(checks).to_csv(OUT/'historical_v5_endpoint_parity.csv', index=False)
    return len(checks)


def main():
    assert os.environ.get('SLURM_JOB_ID'), 'Scientific evaluation requires SLURM.'
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=4)
    args = parser.parse_args()
    start = time.monotonic()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'sample_metrics').mkdir(exist_ok=True)
    initialize()
    hashes_before = {str(p.relative_to(ROOT)): sha(p) for p in SOURCE_PATHS}
    mapping_audit = audit_mapping()
    parity_n = verify_historical()
    cohort = pd.read_csv(NATIVE/'protocol/cohort.csv')
    assert len(cohort) == 121 and cohort.primary.sum() == 97
    with ProcessPoolExecutor(max_workers=args.workers, initializer=initialize) as pool:
        audits = []
        for i, audit in enumerate(pool.map(sample_run, cohort.to_dict('records'))):
            audits.append(audit)
            if (i+1) % 10 == 0 or i+1 == len(cohort):
                print(f'EVALUATED {i+1}/{len(cohort)} samples, {time.monotonic()-start:.1f}s', flush=True)
    data = pd.concat([pd.read_csv(OUT/'sample_metrics'/f'{s}.csv') for s in cohort['sample']], ignore_index=True)
    assert len(data) == 4719 and not data.duplicated(['sample', 'budget', 'route', 'library']).any()
    assert data.core_hvg24.sum() == 2904 and data.fullgene_marker16.sum() == 1936
    core, marker = data[data.core_hvg24].copy(), data[data.fullgene_marker16].copy()
    data.to_csv(OUT/'metrics_all_conditions.csv.gz', index=False)
    core.to_csv(OUT/'metrics_hvg24.csv.gz', index=False)
    marker.to_csv(OUT/'metrics_fullgene_markers.csv.gz', index=False)
    summarize(core, ['budget', 'route']).to_csv(OUT/'summary_hvg24.csv', index=False)
    summarize(marker, ['library']).to_csv(OUT/'summary_fullgene_markers.csv', index=False)
    core[core['sample'].eq('TKU3186')].to_csv(OUT/'TKU3186_hvg24.csv', index=False)
    marker[marker['sample'].eq('TKU3186')].to_csv(OUT/'TKU3186_fullgene_markers.csv', index=False)
    pd.DataFrame(audits).to_csv(OUT/'sample_validation.csv', index=False)
    per_class = pd.read_csv(OUT/'TKU3186_perclass.csv')
    class_mean = per_class.groupby(['budget', 'route', 'library']).F1.mean()
    tk = data[data['sample'].eq('TKU3186')].set_index(['budget', 'route', 'library'])
    np.testing.assert_allclose(class_mean.loc[tk.index], tk.lfine_macroF1, rtol=0, atol=1e-12)
    hashes_after = {str(p.relative_to(ROOT)): sha(p) for p in SOURCE_PATHS}
    assert hashes_before == hashes_after
    outputs = {str(p.relative_to(OUT)): sha(p) for p in OUT.glob('*') if p.is_file() and p.suffix in ['.csv', '.gz', '.json'] and p.name != 'manifest.json'}
    manifest = dict(status='completed' if data.terminal_valid.all() else 'completed_with_unavailable_conditions',
        timestamp=datetime.now(timezone.utc).isoformat(), job=os.environ['SLURM_JOB_ID'],
        host=os.uname().nodename, workers=args.workers, seconds=time.monotonic()-start,
        n_conditions=len(data), n_valid=int(data.terminal_valid.sum()), core_hvg24_rows=len(core),
        fullgene_marker16_rows=len(marker), n_samples=121, n_primary_samples=97,
        n_total_patients=cohort.patient.nunique(), n_primary_patients=cohort.loc[cohort.primary, 'patient'].nunique(),
        source_hashes=hashes_after, outputs=outputs, mapping_audit=mapping_audit,
        historical_v5_endpoint_checks=parity_n, historical_v5_endpoint_checks_passed=parity_n,
        all_sample_compose_lfine_matches_saved_truth=True, TKU3186_perclass_reproduces_macro=True,
        no_fitting=True, no_expression_matrix_loading=True, no_reference_outcome_mapping=True,
        endpoint='v5 set-valued Lfine macro-F1, sample-present classes with support>=20 excluding Other/nan; all cells retained in TP/FP/FN; target sets include every observed reference class',
        cohort_aggregation='Patient mean = equal weight per patient after within-patient available-sample averaging; sample mean also supplied. Fixed 97/24 historical split preserved.',
        compatibility_limit='Broad native-panel meanings are mapped through the frozen semantic parents into compatible Lfine target sets. This is not strict one-to-one fine-type prediction.',
        invalid_terminal_policy='Unavailable/invalid terminal results remain NA and do not become marker-only results',
        display_policy='No author-L1 performance metric emitted; internal source columns are only ontology provenance')
    (OUT/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(json.dumps({k:manifest[k] for k in ['status', 'seconds', 'n_conditions', 'n_valid', 'historical_v5_endpoint_checks_passed']}, indent=2), flush=True)


if __name__ == '__main__':
    main()
