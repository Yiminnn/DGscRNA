#!/usr/bin/env python3
"""Read-only final-prediction audit. All processing runs on SLURM; never fit/train.

Reuses the exact legacy n>=20 one-to-many scoring helpers by AST extraction,
not by importing analysis modules with result-writing or fitting side effects.
"""
import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import subprocess

import numpy as np
import pandas as pd

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT = ROOT / 'results/g274_v5/final_annotations'
COHORT = ROOT / 'results/cohort_v4'
SAMPLE = 'TKU3186'
CUTS = ('none', 'mean', '0.5')
HASHES = {}


def record(path):
    path = Path(path)
    if str(path) not in HASHES:
        with path.open('rb') as f:
            HASHES[str(path)] = hashlib.file_digest(f, 'sha256').hexdigest()
    return path


def dump(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2, ensure_ascii=False, default=str) + '\n')


def csv(path, **kwargs):
    return pd.read_csv(record(path), **kwargs)


def legacy_helpers():
    source = record(ROOT / 'handoff/g274/grid_sets.py')
    tree = ast.parse(source.read_text())
    names = {'compose_lfine', 'lfine_targets', 'macro_f1_lfine', 'LFINE_PREFIX', 'SETS'}
    nodes = [n for n in tree.body if
             (isinstance(n, ast.FunctionDef) and n.name in names) or
             (isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id in names for t in n.targets))]
    assert len(nodes) == 5
    namespace = {'np': np, 'pd': pd}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(source), 'exec'), namespace)
    return namespace


def execution(row):
    if row is None:
        return 'missing_log'
    note = str(row['dl_note'])
    if note == 'ok':
        return 'executed' if str(row['dl_ok']).lower() == 'true' else 'inconsistent_log'
    if note == 'empty_pool':
        return 'no-op' if int(row['n_refine_pool']) == 0 and str(row['dl_ok']).lower() == 'false' else 'inconsistent_log'
    if note == 'lt2_classes':
        return 'fallback_lt2_classes'
    if note.startswith('error'):
        return 'fallback_error'
    return 'unrecognized_log'


def label_scope(pc, helpers):
    lf = helpers['compose_lfine'](pc)
    counts = pd.Series(lf).value_counts()
    classes = [c for c in counts.index if counts[c] >= 20 and c not in ('Other', 'nan')]
    targets = helpers['lfine_targets'](list(counts.index))
    return lf, classes, targets


def metrics(pred, pc, scope, helpers, marker):
    lf, classes, targets = scope
    # Legacy final-grid coverage/correctness: broad L1 predictions receive Lfine
    # one-to-many credit; all cells remain in TP/FP/FN, only class averaging is n>=20.
    abstain = np.isin(pred, ['Unknown', 'Undecided', 'Noise', 'nan', 'None', ''])
    offvocab = (~abstain) & ~np.isin(pred, list(helpers['LFINE_PREFIX']))
    called = ~(abstain | offvocab)
    ok = np.array([g in targets.get(q, ()) for g, q in zip(lf, pred)])
    malignant_gold = pc.L1.astype(str).to_numpy() == 'Malignant'
    malignant_pred = pred == 'Malignant'
    tp = int((malignant_gold & malignant_pred).sum())
    fp = int((~malignant_gold & malignant_pred).sum())
    fn = int((malignant_gold & ~malignant_pred).sum())
    # Prediction-vocabulary ceiling under the existing generous credit rule, not
    # a newly optimized cluster assignment or a strict-label accuracy ceiling.
    reachable = set().union(*(targets.get(label, set()) for label in helpers['panel_labels'][marker]))
    return dict(lfine_macroF1=helpers['macro_f1_lfine'](pred, lf, classes, targets),
                coverage=float((~abstain).mean()), legacy_called_coverage=float(called.mean()),
                abstain_rate=float(abstain.mean()), offvocab_rate=float(offvocab.mean()),
                acc_on_called=float(ok[called].mean()) if called.any() else np.nan,
                malignant_F1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else np.nan,
                n_distinct_calls=int(len(set(pred[called]))),
                n_classes_hit=sum(bool(((lf == c) & ok).any()) for c in classes),
                n_cells=len(pc), lfine_n_classes=len(classes),
                lfine_scored_class_cell_fraction=float(np.isin(lf, classes).mean()),
                marker_vocab_oracle_upper_bound=sum(c in reachable for c in classes)/len(classes) if classes else np.nan,
                n_reference_lfine_disagreements=int((lf != pc.Lfine.astype(str).to_numpy()).sum()))


def self_check(h):
    lf = np.array(['Malignant_OPC', 'TAM_MD', 'Other'], dtype=object)
    targets = h['lfine_targets'](lf)
    assert targets['OPC'] == {'Malignant_OPC'}
    assert h['macro_f1_lfine'](np.array(['OPC', 'TAM', 'Other']), lf, list(lf[:2]), targets) == 1
    assert h['macro_f1_lfine'](np.array(['Unknown']*3), lf, list(lf[:2]), targets) == 0
    for note, flag, pool, status in [('ok', True, 2, 'executed'), ('empty_pool', False, 0, 'no-op'),
                                    ('lt2_classes', False, 2, 'fallback_lt2_classes'),
                                    ('error:test', False, 2, 'fallback_error')]:
        assert execution(dict(dl_note=note, dl_ok=flag, n_refine_pool=pool)) == status
    assert execution(dict(dl_note='ok', dl_ok=False, n_refine_pool=2)) == 'inconsistent_log'
    assert execution(None) == 'missing_log'


def main():
    assert os.environ.get('SLURM_JOB_ID'), 'Data processing and plotting require SLURM'
    args = argparse.ArgumentParser()
    args.add_argument('--inspect', action='store_true')
    inspect = args.parse_args().inspect
    OUT.mkdir(parents=True, exist_ok=True)
    helpers = legacy_helpers()
    self_check(helpers)
    sets = helpers['SETS']
    mapping = csv(ROOT / 'handoff/markers_v3/mapping_L1_v3.csv')
    helpers['panel_labels'] = {m: set(mapping.loc[mapping.marker_set.eq(m), 'gold_class']) for m in sets}
    record(ROOT / 'handoff/g274/annotate_v3.py')
    record(ROOT / 'handoff/g274/cohort_v4.sbatch')
    pc = csv(COHORT / SAMPLE / 'percell.csv', usecols=lambda c: c in ['CellID', 'cluster_final', 'L1', 'L2', 'L3', 'MalState', 'Lfine'] or c.startswith('dl_'))
    log = csv(COHORT / SAMPLE / 'run_log.csv', dtype={'cutoff': str})
    umap = np.load(record(COHORT / SAMPLE / 'umap.npy'), allow_pickle=False)
    assert pc.CellID.is_unique and umap.shape == (len(pc), 2) and np.isfinite(umap).all()
    assert set(log.mcs) == {15}
    assert not log.duplicated(['marker_set', 'cutoff']).any()
    assert set(zip(log.marker_set, log.cutoff)) == {(s, c) for s in sets for c in CUTS}
    scope = label_scope(pc, helpers)
    final = []
    for marker in sets:
        for cutoff in CUTS:
            column = f'dl_{marker}__{cutoff}'
            row = log[(log.marker_set == marker) & (log.cutoff == cutoff)].iloc[0]
            final.append(dict(sample=SAMPLE, marker_set=marker, cutoff=cutoff,
                              stage='after DL', prediction_column=column,
                              execution_status=execution(row), dl_note=row.dl_note,
                              n_refine_pool=int(row.n_refine_pool), n_panels=int(row.n_panels), mcs=int(row.mcs),
                              source_percell=str(COHORT / SAMPLE / 'percell.csv'),
                              **metrics(pc[column].astype(str).to_numpy(), pc, scope, helpers, marker)))
    final = pd.DataFrame(final)
    assert len(final) == 24
    final.to_csv(OUT / 'TKU3186_final_24.csv', index=False)
    ref = csv(ROOT / 'results/refsel/TKU3186_full_grid.csv', dtype={'cutoff': str})
    ref = ref.loc[ref.stage.eq('after DL')].copy()
    assert len(ref) == 24 and not ref.duplicated(['marker_set', 'cutoff']).any()
    comparison = []
    for row in final.to_dict('records'):
        old = ref[(ref.marker_set == row['marker_set']) & (ref.cutoff == row['cutoff'])].iloc[0]
        for key in ['lfine_macroF1', 'coverage', 'acc_on_called', 'malignant_F1', 'n_distinct_calls', 'n_classes_hit']:
            match = bool(np.isclose(row[key], old[key], rtol=0, atol=0.000050001, equal_nan=True))
            comparison.append(dict(marker_set=row['marker_set'], cutoff=row['cutoff'], metric=key,
                                   percell_final_value=row[key], refsel_after_DL_value=old[key],
                                   difference=row[key]-old[key], matches_4dp=match))
    pd.DataFrame(comparison).to_csv(OUT / 'TKU3186_refsel_after_DL_check.csv', index=False)
    # Compare only like-named final columns, IDs, reference and coordinates.
    frozen = ROOT / 'handoff/FROZEN_single_sample_TKU3186/results_mcs15'
    fpc = csv(frozen / 'percell.csv', usecols=list(pc.columns))
    fumap = np.load(record(frozen / 'umap.npy'), allow_pickle=False)
    obs = csv(ROOT / 'data_bench/GSE274546/mtx/TKU3186/obs.csv', usecols=['CellID', 'L1', 'L3', 'MalState'])
    ordered = list(pc.CellID) == list(fpc.CellID)
    obs_ordered = list(pc.CellID) == list(obs.CellID)
    lineage = dict(n_cells=len(pc), umap_shape=list(umap.shape),
                   cohort_v4_frozen_cell_order_equal=ordered,
                   cohort_v4_frozen_umap_array_equal=bool(np.array_equal(umap, fumap)),
                   cohort_v4_frozen_column_mismatches={c: int((pc[c].fillna('<NA>').astype(str) != fpc[c].fillna('<NA>').astype(str)).sum()) for c in pc.columns} if ordered else None,
                   input_obs_cell_order_equal=obs_ordered,
                   input_obs_composed_lfine_mismatches=int((scope[0] != helpers['compose_lfine'](obs)).sum()) if obs_ordered else None,
                   stored_lfine_vs_legacy_compose_mismatches=int((scope[0] != pc.Lfine.astype(str).to_numpy()).sum()),
                   refsel_after_DL_metric_checks=len(comparison),
                   refsel_after_DL_metric_matches=sum(r['matches_4dp'] for r in comparison),
                   reference_labels='compose_lfine(percell L1/L3/MalState); not a different run',
                   umap_lineage_limit='npy contains no CellIDs; row alignment supported by same-directory paired writer and frozen-array/ordered-CellID check, not embedded identity metadata',
                   projection='Existing HVG2000 / scaled PCA30 / neighbors15 / UMAP / HDBSCAN min_cluster_size=min_samples=15; no recomputation',
                   refsel_table_policy='Only stage=after DL retained; mismatches are not substituted into plots or current metrics')
    dump('TKU3186_lineage.json', lineage)
    labels = sorted(set(pc[[f'dl_{m}__none' for m in sets]].astype(str).to_numpy().ravel()))
    dump('inspection.json', dict(labels=labels, lineage=lineage, final=final.to_dict('records')))
    print(json.dumps(dict(labels=labels, lineage=lineage, final_none=final[final.cutoff.eq('none')].to_dict('records')), default=str), flush=True)
    if inspect:
        return
    cohort_audit(sets, helpers)
    plot_final(pc, umap, sets, final)
    unchanged = {}
    for path, digest in HASHES.items():
        with open(path, 'rb') as f:
            unchanged[path] = hashlib.file_digest(f, 'sha256').hexdigest() == digest
    assert all(unchanged.values()), 'A read-only input changed during the audit'
    dump('source_hashes.json', HASHES)
    dump('self_check.json', dict(status='passed', slurm_job_id=os.environ['SLURM_JOB_ID'],
                                host=os.uname().nodename, final_rows=24,
                                read_only_sources_checked=len(unchanged), all_sources_unchanged=True,
                                refsel_matches=sum(r['matches_4dp'] for r in comparison),
                                refsel_checks=len(comparison), no_fit_no_train_no_projection=True))
    print('FINAL_ANNOTATIONS_COMPLETE', flush=True)


def cohort_audit(sets, helpers):
    old_path = ROOT / 'results/ablation/SUPP_per_sample_markerset.csv'
    old = csv(old_path).set_index('sample', verify_integrity=True)
    expected = record(ROOT / 'handoff/g274/cohort_samples.txt').read_text().splitlines()
    samples = sorted(set(expected) | set(old.index) | {p.name for p in COHORT.iterdir() if p.is_dir() and not p.name.startswith('.')})
    rows, inventory, logs = [], [], []
    for sample in samples:
        directory = COHORT / sample
        logpath, pcpath = directory / 'run_log.csv', directory / 'percell.csv'
        run = csv(logpath, dtype={'cutoff': str}) if logpath.exists() else pd.DataFrame()
        if not run.empty:
            assert not run.duplicated(['marker_set', 'cutoff']).any(), sample
            for item in run.to_dict('records'):
                logs.append(dict(**item, execution_status=execution(item), source_run_log=str(logpath)))
        skipped = directory / 'UNEVALUABLE.csv'
        status = 'available'
        reason = ''
        if skipped.exists():
            status, reason = 'skip_UNEVALUABLE', csv(skipped).to_json(orient='records')
        elif not pcpath.exists():
            status, reason = 'missing_percell', 'No persisted final per-cell file'
        item = dict(sample=sample, in_old_wide=sample in old.index, sample_status=status,
                    reason=reason, has_run_log=logpath.exists(), has_percell=pcpath.exists(),
                    source_percell=str(pcpath))
        pc = None
        if status == 'available':
            wanted = {'CellID', 'L1', 'L3', 'MalState', 'Lfine', 'cluster_final'} | {f'dl_{m}__none' for m in sets}
            pc = csv(pcpath, usecols=lambda c: c in wanted)
            assert pc.CellID.is_unique, sample
            scope = label_scope(pc, helpers)
            item.update(n_cells=len(pc), lfine_n_classes=len(scope[1]),
                        stored_lfine_mismatches=int((scope[0] != pc.Lfine.astype(str).to_numpy()).sum()),
                        n_clusters=int(pc.cluster_final[pc.cluster_final.ne('Noise')].nunique()),
                        n_noise=int(pc.cluster_final.eq('Noise').sum()))
        inventory.append(item)
        for marker in sets:
            column = f'dl_{marker}__none'
            lr = None
            if not run.empty:
                hit = run[run.marker_set.eq(marker) & run.cutoff.eq('none')]
                lr = hit.iloc[0] if len(hit) == 1 else None
            row = dict(sample=sample, marker_set=marker, cutoff='none', stage='after DL',
                       prediction_column=column, source_percell=str(pcpath),
                       execution_status=execution(lr), dl_note=str(lr.dl_note) if lr is not None else '',
                       mcs=int(lr.mcs) if lr is not None else np.nan,
                       in_old_wide=sample in old.index, final_available=False,
                       old_wide_value=float(old.loc[sample, marker]) if sample in old.index else np.nan,
                       comparison_status='unverified_' + status)
            if pc is not None and column in pc:
                row.update(metrics(pc[column].astype(str).to_numpy(), pc, scope, helpers, marker))
                row['final_available'] = True
                if sample not in old.index:
                    row['comparison_status'] = 'not_in_old_wide'
                elif not np.isfinite(row['old_wide_value']) or not np.isfinite(row['lfine_macroF1']):
                    row['comparison_status'] = 'unverified_missing_score'
                else:
                    row['difference'] = row['lfine_macroF1'] - row['old_wide_value']
                    row['comparison_status'] = 'matched_4dp' if abs(row['difference']) <= 0.000050001 else 'mismatch'
            elif pc is not None:
                row['comparison_status'] = 'unverified_missing_final_column'
            rows.append(row)
        print(f'COHORT {sample} {status}', flush=True)
    result, inv, runlogs = pd.DataFrame(rows), pd.DataFrame(inventory), pd.DataFrame(logs)
    inv.to_csv(OUT / 'cohort_sample_inventory.csv', index=False)
    runlogs.to_csv(OUT / 'cohort_run_log_statuses.csv', index=False)
    result.to_csv(OUT / 'cohort_none_final_vs_old_wide.csv', index=False)
    result[result.in_old_wide & result.comparison_status.ne('matched_4dp')].to_csv(OUT / 'cohort_old_wide_differences.csv', index=False)
    finalcols = [c for c in result if c not in ('old_wide_value', 'difference', 'comparison_status')]
    result.loc[result.final_available, finalcols].to_csv(OUT / 'cohort_none_final_metrics.csv', index=False)
    statuses = runlogs.groupby(['cutoff', 'execution_status']).size().rename('n_arms').reset_index()
    statuses.to_csv(OUT / 'cohort_execution_counts_all_cutoffs.csv', index=False)
    aggregate = []
    for marker in sets:
        sub = result[result.marker_set.eq(marker)]
        available = sub[sub.final_available]
        oldscope = sub[sub.in_old_wide]
        aggregate.append(dict(marker_set=marker, cutoff='none', n_final_samples=len(available),
                              n_old_wide_samples=len(oldscope),
                              n_old_wide_matched=int(oldscope.comparison_status.eq('matched_4dp').sum()),
                              n_old_wide_mismatch=int(oldscope.comparison_status.eq('mismatch').sum()),
                              n_old_wide_unverified=int(oldscope.comparison_status.str.startswith('unverified_').sum()),
                              n_executed=int(available.execution_status.eq('executed').sum()),
                              n_no_op=int(available.execution_status.eq('no-op').sum()),
                              n_fallback_lt2=int(available.execution_status.eq('fallback_lt2_classes').sum()),
                              n_fallback_error=int(available.execution_status.eq('fallback_error').sum()),
                              n_missing_or_inconsistent_log=int((~available.execution_status.isin(['executed','no-op','fallback_lt2_classes','fallback_error'])).sum()),
                              final_lfine_macroF1_median=float(available.lfine_macroF1.median()),
                              final_coverage_median=float(available.coverage.median())))
    pd.DataFrame(aggregate).to_csv(OUT / 'cohort_none_final_summary.csv', index=False)
    oldscope = result[result.in_old_wide]
    matched = oldscope.comparison_status.eq('matched_4dp')
    verified_samples = oldscope.groupby('sample').comparison_status.apply(lambda s: len(s) == 8 and s.eq('matched_4dp').all())
    summary = dict(expected_samples=len(set(expected)), sample_directories_audited=len(samples),
                   old_wide_samples=len(old), old_wide_marker_values=len(oldscope),
                   final_samples=int(inv.sample_status.eq('available').sum()),
                   final_arms=int(result.final_available.sum()),
                   skipped_samples=inv.loc[inv.sample_status.ne('available'), ['sample','sample_status','reason']].to_dict('records'),
                   old_wide_matches=int(matched.sum()),
                   old_wide_mismatches=int(oldscope.comparison_status.eq('mismatch').sum()),
                   old_wide_unverified=int(oldscope.comparison_status.str.startswith('unverified_').sum()),
                   old_samples_all_8_match=int(verified_samples.sum()),
                   none_final_execution_counts=result.loc[result.final_available, 'execution_status'].value_counts().to_dict(),
                   all_cutoff_runlog_counts=statuses.to_dict('records'),
                   stored_lfine_mismatch_samples=inv.loc[inv.stored_lfine_mismatches.fillna(0).gt(0), ['sample','stored_lfine_mismatches']].to_dict('records'),
                   read_scope='All cohort_v4 run logs; only references, IDs, clusters and eight dl_*__none columns in percell. Other cutoffs NOT rescored across cohort.',
                   scoring_source=str(ROOT / 'handoff/g274/grid_sets.py'),
                   prediction_writer=str(ROOT / 'handoff/g274/annotate_v3.py'),
                   prediction_job=str(ROOT / 'handoff/g274/cohort_v4.sbatch'),
                   old_wide_source=str(old_path),
                   old_wide_producer='Not found among workspace .py/.sh/.sbatch sources; fig_selection.py and sample_inventory.py only consume it',
                   scoring_rule='Legacy n>=20 class average; Other/nan excluded from averaged classes only; all cells in score; all present Lfine classes build target sets. No new 97-rule.',
                   interpretation='Final means persisted dl_ output, not proof of neural update. ok=executed; empty_pool=no-op; error/lt2_classes=fallback, never executed.',
                   old_wide_retention='Missing stage field is not evidence of before-DL. Keep original as historical report record. Only matched sample-marker values are empirically verified against this final source; mismatch and missing-source cells cannot support an unqualified final-only cohort claim.',
                   whole_old_wide_verified=bool(matched.all()),
                   claim_limit='Reference concordance under generous one-to-many credit, not independent accuracy. BrainAtlas112 and CARE_TME overlap reference-label construction. Matching does not validate automatic marker selection superiority.')
    dump('cohort_final_audit.json', summary)
    print('COHORT_SUMMARY ' + json.dumps(summary, default=str), flush=True)


def plot_final(pc, umap, sets, final):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from PIL import Image

    # Reuse the current v5 paint helper and its fixed three-color/shape encoding;
    # unlike importing a fitting script, AST extraction cannot execute its main.
    source = record(ROOT / 'handoff/g274_table4/plot_results.py')
    tree = ast.parse(source.read_text())
    selected = [n for n in tree.body if
                (isinstance(n, ast.FunctionDef) and n.name == 'paint') or
                (isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id in {'COLORS','SHAPES'} for t in n.targets))]
    assert len(selected) == 3
    ns = {'np': np, 'Line2D': Line2D}
    exec(compile(ast.Module(body=selected, type_ignores=[]), str(source), 'exec'), ns)
    # The skill's /tmp path is login-node local; stage its unchanged source in
    # the allowed output directory and evaluate it as ESM on the SLURM node.
    validator = (OUT / 'palette_validator_source.txt').read_text()
    validation = subprocess.run(['node', '--input-type=module', '-', ','.join(ns['COLORS']),
                                 '--mode', 'light', '--pairs', 'all'],
                                input="process.argv[1] = 'validate_palette.js';\n" + validator,
                                text=True, capture_output=True)
    (OUT / 'palette_validation.txt').write_text(validation.stdout + validation.stderr)
    assert validation.returncode == 0 and 'ALL CHECKS PASS' in validation.stdout, validation.stderr + validation.stdout
    plt.rcParams.update({'font.size': 10, 'figure.facecolor': '#fcfcfb', 'axes.facecolor': '#fcfcfb',
                         'axes.spines.top': False, 'axes.spines.right': False, 'text.color': '#0b0b0b'})
    columns = [f'dl_{m}__none' for m in sets]
    order = sorted(set(pc[columns].astype(str).to_numpy().ravel()))
    assert len(order) <= len(ns['COLORS']) * len(ns['SHAPES'])
    palette = {label: dict(color=ns['COLORS'][i % 3], shape=ns['SHAPES'][i // 3]) for i, label in enumerate(order)}
    dump('prediction_palette.json', dict(source=str(source), encoding=palette,
                                        note='Fixed full-vocabulary hue x shape encoding across all eight panels; never recolored by per-panel rank. Table provides exact cell labels. Static light-theme PNG.'))
    cells = pc[['CellID'] + columns].copy()
    cells['UMAP1'], cells['UMAP2'] = umap[:, 0], umap[:, 1]
    cells.to_csv(OUT / 'TKU3186_final_none_cells.csv', index=False)
    counts, index, bounds = [], [], []
    panel_rows = final[final.cutoff.eq('none')].set_index('marker_set')

    def check_text(fig, name):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        width, height = fig.canvas.get_width_height()
        for text in fig.findobj(matplotlib.text.Text):
            if not text.get_visible() or not text.get_text():
                continue
            box = text.get_window_extent(renderer)
            assert box.x0 >= -1 and box.y0 >= -1 and box.x1 <= width+1 and box.y1 <= height+1, (name, text.get_text(), box)
            bounds.append(dict(figure=name, text=text.get_text(), bounds=list(box.extents)))

    caption = ('Input HVG2000 / scaled PCA30 / neighbors15 / UMAP / HDB15 (min_cluster_size=min_samples=15).\n'
               'All 4,384 cells; identical existing UMAP; no re-projection. Only dl_<marker>__none final predictions are displayed.\n'
               'Final does not imply a neural update: executed = run_log ok; no-op = empty_pool. No before-stage scores.\n'
               'Scores: legacy Lfine one-to-many concordance, n>=20 (22 classes); not independent accuracy.\n'
               '* BrainAtlas112 / CARE_TME contributed to the reference; oracle upper bound is marker-vocabulary reachability, not a fitted optimum.')
    for marker in sets:
        row = panel_rows.loc[marker]
        pred = pc[f'dl_{marker}__none'].astype(str).to_numpy()
        assert len(pred) == len(umap)
        counts.extend(dict(marker_set=marker, cutoff='none', label=label, n_cells=int((pred == label).sum()), **palette[label]) for label in order)
        fig, ax = plt.subplots(figsize=(12, 8))
        handles = ns['paint'](ax, umap, pred, order)
        ax.legend(handles=handles, loc='upper left', bbox_to_anchor=(1.02, 1), fontsize=10, frameon=False)
        title = marker + (' *' if marker in ('BrainAtlas112','CARE_TME') else '')
        ax.set_title(f'{SAMPLE} | {title}\nFinal prediction, cutoff=none | {row.execution_status} ({row.dl_note})', loc='left', fontsize=13, pad=36)
        ax.text(0, 1.012, f'Lfine macro-F1 {row.lfine_macroF1:.4f} | coverage {row.coverage:.2%} | oracle bound {row.marker_vocab_oracle_upper_bound:.3f}', transform=ax.transAxes, fontsize=10)
        ax.set_xlabel('Existing UMAP 1'); ax.set_ylabel('Existing UMAP 2')
        fig.text(.025, .025, caption, fontsize=8.5, va='bottom', linespacing=1.5)
        fig.subplots_adjust(left=.075, right=.67, bottom=.23, top=.83)
        name = f'TKU3186_final_{marker}__none.png'
        check_text(fig, name)
        fig.savefig(OUT / name, dpi=160)
        plt.close(fig)
        index.append(dict(marker_set=marker, cutoff='none', prediction_column=f'dl_{marker}__none',
                          execution_status=row.execution_status, n_cells=len(pred),
                          source_percell=str(COHORT / SAMPLE / 'percell.csv'), source_umap=str(COHORT / SAMPLE / 'umap.npy'),
                          figure=str(OUT / name), legend_labels=len(handles)))
    fig, axes = plt.subplots(2, 4, figsize=(23, 13))
    for ax, marker in zip(axes.ravel(), sets):
        row = panel_rows.loc[marker]
        ns['paint'](ax, umap, pc[f'dl_{marker}__none'].astype(str).to_numpy(), order)
        ax.set_title(marker + (' *' if marker in ('BrainAtlas112','CARE_TME') else '') +
                     f'\n{row.execution_status} | F1 {row.lfine_macroF1:.4f} | coverage {row.coverage:.2%}\nOracle bound {row.marker_vocab_oracle_upper_bound:.3f}', fontsize=11, loc='left')
    legend = [Line2D([], [], marker=palette[label]['shape'], markerfacecolor=palette[label]['color'],
                     markeredgecolor=palette[label]['color'], linestyle='none', label=label, markersize=8) for label in order]
    fig.legend(handles=legend, loc='lower center', bbox_to_anchor=(.5, .17), ncol=6, fontsize=11, frameon=False)
    fig.suptitle('TKU3186 | Eight marker sets | Final dl_ predictions | cutoff=none', fontsize=20, y=.97)
    fig.text(.035, .02, caption, fontsize=11, va='bottom', linespacing=1.5)
    fig.subplots_adjust(left=.035, right=.985, top=.87, bottom=.28, hspace=.36, wspace=.2)
    name = 'TKU3186_eight_markers_final_none.png'
    check_text(fig, name)
    fig.savefig(OUT / name, dpi=150)
    plt.close(fig)
    pd.DataFrame(counts).to_csv(OUT / 'TKU3186_final_label_counts.csv', index=False)
    pd.DataFrame(index).to_csv(OUT / 'figure_index.csv', index=False)
    assert pd.DataFrame(counts).groupby('marker_set').n_cells.sum().eq(len(pc)).all()
    dump('figure_text_bounds.json', bounds)
    with Image.open(OUT / name) as im:
        im.verify()
    # One audit figure for the cohort validation result, without before-stage metrics.
    summary = pd.read_csv(OUT / 'cohort_none_final_summary.csv')
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5), sharey=True)
    for ax, key, title in zip(axes, ['n_old_wide_matched','n_old_wide_mismatch','n_old_wide_unverified'],
                             ['Matches final dl_ output', 'Differs from final dl_ output', 'No final output to verify']):
        ax.barh(summary.marker_set, summary[key], height=.45, color=ns['COLORS'][0])
        for i, value in enumerate(summary[key]):
            ax.text(value+1, i, str(value), va='center')
        ax.set_xlim(0, int(summary.n_old_wide_samples.max())+12)
        ax.set_title(title, fontsize=11)
        ax.set_xlabel('Old-table sample-marker values')
    axes[0].invert_yaxis()
    fig.suptitle('Old cohort table verification | Eight dl_*__none final columns | Legacy n>=20 rule', fontsize=14)
    fig.text(.02, .02, 'Numerical match tolerance: 0.000050001 (old table rounded to 4 decimals). Missing stage is not evidence of a before-DL result.\n'
             'Unverified entries remain historical report records; mismatches are source differences, not DL gains or before-stage metrics.', fontsize=9)
    fig.subplots_adjust(left=.18, right=.97, top=.85, bottom=.18, wspace=.2)
    check_text(fig, 'cohort_final_verification.png')
    fig.savefig(OUT / 'cohort_final_verification.png', dpi=150)
    plt.close(fig)
    print('PLOTS_COMPLETE 8 individual + overview + cohort audit', flush=True)


if __name__ == '__main__':
    main()
