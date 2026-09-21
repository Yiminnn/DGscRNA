#!/usr/bin/env python3
"""Frozen B saved-prediction Lfine evaluation. All execution requires SLURM."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import importlib.util
import json
import os
import time

ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
HERE = Path(__file__).resolve().parent
BASE = ROOT/'results/hvg_ptc_20260916_v1'
OLD = BASE/'reviewer_completion_20260920/controls'
NATIVE = BASE/'paper_claim_validation_20260917'
COMPACT = BASE/'lfine_compact_20260920'
OUT = BASE/'reviewer_completion_20260920/controls_lfine_v1'
LIBRARIES = ['CM2_glioma_other', 'CM2_primary_all_context']
PILOTS = ['TKU4163', 'NL022', 'SN040']


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def write_json(path, value):
    path = Path(path)
    tmp = path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')
    tmp.replace(path)


def utc():
    return datetime.now(timezone.utc).isoformat()


def verify_freeze():
    frozen = json.loads((HERE/'freeze.json').read_text())
    for name, digest in frozen['files'].items():
        assert sha(ROOT/name) == digest, ('frozen source changed', name)
    return frozen


def complete(path, manifest, flag):
    digest = sha(path/manifest)
    assert (path/flag).read_text().strip() == digest, (path, flag)
    return json.loads((path/manifest).read_text()), digest


def main():
    assert os.environ.get('SLURM_JOB_ID'), 'Use SLURM for scientific evaluation.'
    parser = argparse.ArgumentParser()
    parser.add_argument('--mode', choices=['pilot', 'full'], required=True)
    args = parser.parse_args()
    frozen = verify_freeze()
    if args.mode == 'full':
        authorization = json.loads((OUT/'FULL_AUTHORIZED.json').read_text())
        assert authorization['protocol_sha256'] == sha(HERE/'PROTOCOL.md')
        pilot = json.loads((OUT/'pilot/validation.json').read_text())
        assert pilot['status'] == 'passed' and pilot['n_rows'] == 24
    start = time.monotonic()
    dest = OUT/args.mode
    dest.mkdir(exist_ok=False)
    job = {'job': os.environ['SLURM_JOB_ID'], 'step': os.environ.get('SLURM_STEP_ID')}
    status = dict(work_package='B', stage='B_LFINE', status='running',
        updated_at=utc(), jobs=[job], completed=0, remaining=24 if args.mode=='pilot' else 408,
        summary='Saved-terminal Lfine compatible-target-set evaluation; no fitting.',
        evidence=[str(HERE/'PROTOCOL.md'), str(HERE/'freeze.json')])
    write_json(OUT/'status.json', status)
    import numpy as np
    import pandas as pd
    spec = importlib.util.spec_from_file_location('frozen_compact',
        HERE/'reference_sources/compact_evaluate_lfine.py')
    provider = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(provider)
    provider.initialize()
    v5, helpers, mappings = provider.v5, provider.helpers, provider.mappings
    fields = provider.SCORE_FIELDS
    tasks = json.loads((OLD/'tasks.json').read_text())
    assert len(tasks) == 84
    old_summary, old_summary_hash = complete(OLD/'summary', 'manifest.json', 'COMPLETE')
    cohort = pd.read_csv(NATIVE/'protocol/cohort.csv').set_index('sample')
    truth, scopes, truth_hashes = {}, {}, {}
    for sample in PILOTS:
        tp = NATIVE/'evaluation_inputs'/sample/'truth.csv.gz'
        im = json.loads((NATIVE/'inputs'/sample/'input_manifest.json').read_text())
        assert sha(tp) == im['evaluation_files']['truth.csv.gz']
        t = pd.read_csv(tp, dtype=str, keep_default_na=False).rename(
            columns={'cell_id':'CellID', 'lfine_original':'Lfine'})
        assert t.CellID.is_unique and t.Patient.eq(cohort.loc[sample, 'patient']).all()
        scope = v5.label_scope(t, helpers)
        np.testing.assert_array_equal(scope[0], t.Lfine.to_numpy())
        truth[sample], scopes[sample], truth_hashes[sample] = t, scope, sha(tp)
    mapping = pd.read_csv(NATIVE/'markers/panel_L1_mapping.csv', dtype=str, keep_default_na=False)
    old_map = pd.read_csv(ROOT/'handoff/markers_v3/mapping_L1_v3.csv', dtype=str, keep_default_na=False)
    audit = mapping[mapping.library.isin(LIBRARIES)].merge(
        old_map[['marker_set','panel','gold_class']], left_on=['library','panel'],
        right_on=['marker_set','panel'], how='left', validate='one_to_one')
    audit['old_v5_panel_available'] = audit.gold_class.notna()
    audit['exact_old_v5_parent'] = audit.old_v5_panel_available & audit.L1.eq(audit.gold_class)
    audit['has_lfine_targets'] = audit.L1.isin(helpers['LFINE_PREFIX'])
    assert audit[audit.library.eq(LIBRARIES[0])].exact_old_v5_parent.all()
    audit.rename(columns={'L1':'semantic_parent'}).to_csv(dest/'mapping_audit.csv.gz', index=False)
    write_json(dest/'lfine_target_prefixes.json', helpers['LFINE_PREFIX'])
    rows, class_rows, parity, artifacts = [], [], [], {}
    compact = pd.read_csv(COMPACT/'metrics_hvg24.csv.gz') if args.mode=='pilot' else None

    def record(path, expected=None):
        digest = sha(path)
        assert expected is None or digest == expected, ('artifact hash', str(path))
        artifacts[str(path.relative_to(ROOT))] = digest
        return digest

    def score_prediction(cfg, pp, terminal_path, provenance, library, epoch, thresholds, route):
        sample, budget = cfg['sample'], cfg['budget']
        pred = pd.read_csv(pp, dtype=str, keep_default_na=False)
        t, scope = truth[sample], scopes[sample]
        np.testing.assert_array_equal(pred.cell_id, t.CellID)
        assert pred.cell_id.is_unique
        known = pred.initial.ne('Undecided')
        with np.load(terminal_path, allow_pickle=False) as terminal:
            for field in ['initial', 'final090', 'final070']:
                np.testing.assert_array_equal(pred[field].to_numpy(), terminal[field])
                if field != 'initial':
                    np.testing.assert_array_equal(pred.loc[known,'initial'], pred.loc[known,field])
        for digits in thresholds:
            column, stage = 'final'+digits, 'terminal'+digits
            mapped = np.asarray([mappings[library].get(v,
                'Unknown' if v in provider.ABSTAIN else 'UNMAPPABLE') for v in pred[column]], dtype=object)
            values = v5.metrics(mapped, t, scope, helpers, library)
            row = dict(task=cfg['task'], sample=sample, patient=cohort.loc[sample,'patient'],
                primary=bool(cohort.loc[sample,'primary']), budget=budget, route=route,
                control=cfg.get('name', 'learning'), library=library, cutoff='mean', stage=stage,
                model_seed=cfg.get('model_seed',42), epochs=epoch,
                snn_k=cfg.get('snn_k',20), umap_neighbors=cfg.get('umap_neighbors',30),
                n_cells=len(t), status='completed', terminal_valid=True,
                dl_status=provenance['dl_status'], training_executed=provenance['training_executed'],
                n_training_classes=provenance['n_training_classes'], n_pool=provenance['n_pool'],
                backup=False, reference_labels_used_for_fit=False,
                truth_sha256=truth_hashes[sample], predictions_path=str(pp.relative_to(ROOT)),
                predictions_sha256=record(pp), terminal_sha256=record(terminal_path),
                source_native_labels_preserved=True,
                **{key:values[key] for key in fields})
            assert row['n_reference_lfine_disagreements'] == 0
            row_key = '|'.join(str(row[k]) for k in ['task','sample','budget','control','library','model_seed','epochs','stage'])
            row['row_key'] = row_key
            lf, classes, targets = scope
            ok = np.asarray([g in targets.get(q, ()) for g,q in zip(lf,mapped)])
            f1s = []
            for label in classes:
                g = lf == label
                has = np.asarray([label in targets.get(q, ()) for q in mapped])
                tp, fn, fp = int((g&ok).sum()), int((g&~ok).sum()), int((~ok&has&~g).sum())
                f1 = 2*tp/(2*tp+fn+fp) if 2*tp+fn+fp else 0.0
                f1s.append(f1)
                class_rows.append(dict(row_key=row_key, lfine_class=label, support=int(g.sum()), TP=tp,FN=fn,FP=fp,F1=f1))
            independent = float(np.mean(f1s)) if f1s else float('nan')
            np.testing.assert_allclose(independent, row['lfine_macroF1'], rtol=0,atol=1e-12,equal_nan=True)
            if args.mode == 'pilot':
                reference = compact[compact['sample'].eq(sample)&compact.budget.eq(budget)&compact.route.eq(route)&compact.library.eq(library)]
                assert len(reference)==1 and reference.iloc[0].terminal_valid
                ref = reference.iloc[0]
                for field in fields:
                    np.testing.assert_allclose(row[field], ref[field], rtol=0,atol=1e-12,equal_nan=True)
                refpp = ROOT/ref.predictions_path
                record(refpp, ref.predictions_sha256)
                native = pd.read_csv(refpp,dtype=str,keep_default_na=False)
                for field in ['cell_id','initial',column]:
                    np.testing.assert_array_equal(pred[field],native[field])
                parity.append(dict(row_key=row_key, n_metrics=len(fields), all_metrics_equal_1e12=True,
                    native_calls_identical=True, reference_predictions_path=ref.predictions_path,
                    lfine_macroF1_difference=row['lfine_macroF1']-ref.lfine_macroF1))
            rows.append(row)

    for cfg in tasks:
        if cfg['task'] == 'neighbors':
            if args.mode=='pilot' and (cfg['snn_k']!=20 or cfg['umap_neighbors']!=30):
                continue
            p = OLD/'neighbors'/cfg['sample']/cfg['budget']/cfg['name']
            manifest, mh = complete(p,'manifest.json','COMPLETE')
            record(p/'manifest.json', old_summary['verified_artifacts'][str(p/'manifest.json')])
            for key,value in cfg.items(): assert manifest['config'][key]==value,(key,p)
            route = cfg['space']+'_'+cfg['method']
            score = p/route
            sm, sh = complete(score,'score_manifest.json','SCORE_COMPLETE')
            record(score/'score_manifest.json',sh)
            assert not sm['reference_labels_used_for_fit']
            for library in LIBRARIES[:1] if args.mode=='pilot' else LIBRARIES:
                arms=[aid for aid,arm in sm['arms'].items() if arm['library']==library and str(arm['cutoff'])=='mean']
                assert len(arms)==1
                td=score/'terminal'/arms[0]
                tm, th = complete(td,'terminal_manifest.json','TERMINAL_COMPLETE')
                assert tm['terminal_valid'] and tm['status']=='completed'
                assert tm['score_manifest_sha256']==sh and not tm['reference_labels_used_for_fit']
                assert tm['arm']['library']==library and str(tm['arm']['cutoff'])=='mean'
                record(td/'terminal_manifest.json',th)
                record(td/'training_manifest.json',tm['training_manifest_sha256'])
                record(td/'terminal.npz',tm['terminal_sha256'])
                record(td/'predictions.csv.gz',tm['predictions_sha256'])
                score_prediction(cfg,td/'predictions.csv.gz',td/'terminal.npz',tm,library,10,
                    ['090'] if args.mode=='pilot' else ['090','070'],route)
        else:
            assert cfg['task']=='learning'
            if args.mode=='pilot' and cfg['model_seed']!=42: continue
            p=OLD/'learning'/cfg['sample']/cfg['budget']/f"seed{cfg['model_seed']}"/LIBRARIES[0]
            manifest,mh=complete(p,'manifest.json','LEARNING_COMPLETE')
            assert manifest['config']==cfg and not manifest['backup']
            record(p/'manifest.json',old_summary['verified_artifacts'][str(p/'manifest.json')])
            tr=json.loads((p/'training_manifest.json').read_text())
            record(p/'training_manifest.json',manifest['training_manifest_sha256'])
            assert tr['terminal_valid'] and tr['training_executed']
            assert not tr['provenance']['reference_labels_used_for_fit']
            assert tr['provenance']['library']==LIBRARIES[0] and tr['provenance']['cutoff']=='mean'
            assert tr['params']['model_seed']==cfg['model_seed'] and tr['params']['epochs']==30
            assert tr['params']['split_seed']==42
            record(Path(tr['provenance']['source_score'])/'score_manifest.json',tr['provenance']['score_sha256'])
            for cp in tr['checkpoints']:
                if args.mode=='pilot' and cp['epoch']!=10: continue
                td=p/f"epoch{cp['epoch']:02d}"
                record(td/'terminal.npz',cp['terminal_sha256'])
                record(td/'predictions.csv.gz',old_summary['verified_artifacts'][str(td/'predictions.csv.gz')])
                score_prediction(cfg,td/'predictions.csv.gz',td/'terminal.npz',tr,LIBRARIES[0],cp['epoch'],
                    ['090'] if args.mode=='pilot' else ['090','070'],'UMAP2_HDBSCAN_R')
    result=pd.DataFrame(rows)
    expected=24 if args.mode=='pilot' else 408
    assert len(result)==expected and result.row_key.is_unique
    assert not any(c in result for c in ['accuracy','macroF1_present','malignant_F1'])
    assert result.terminal_valid.all() and result.n_reference_lfine_disagreements.eq(0).all()
    result.to_csv(dest/'terminal_metrics.csv.gz',index=False)
    pd.DataFrame(class_rows).to_csv(dest/'per_class_counts.csv.gz',index=False)
    if args.mode=='pilot': pd.DataFrame(parity).to_csv(dest/'compact_parity.csv',index=False)
    write_json(dest/'artifact_hashes.json',artifacts)
    verify_freeze()
    for name,digest in artifacts.items(): assert sha(ROOT/name)==digest,('artifact changed during evaluation',name)
    validation=dict(status='passed',mode=args.mode,endpoint='Lfine compatible-target-set metrics',
        completed_at=utc(),**job,elapsed_seconds=time.monotonic()-start,n_rows=len(result),
        n_neighbor_rows=int(result.task.eq('neighbors').sum()),n_checkpoint_rows=int(result.task.eq('learning').sum()),
        n_independent_perclass_rows=len(class_rows),n_compact_parity_rows=len(parity),
        n_compact_fields=len(fields),no_fitting=True,no_expression_matrix_loading=True,
        no_best_selection=True,no_L1_performance_values=True,all_native_labels_preserved=True,
        all_known_seed_labels_retained=True,all_terminal_hashes_verified=True,
        sources_unchanged=True,freeze_sha256=sha(HERE/'freeze.json'),protocol_sha256=sha(HERE/'PROTOCOL.md'),
        old_summary_manifest_sha256=old_summary_hash,
        limitation='Three frozen GBM pilots, HVG2000/5000 only; compatible target sets are not strict fine-type predictions.',
        files={p.name:sha(p) for p in dest.iterdir() if p.is_file()})
    write_json(dest/'validation.json',validation)
    status.update(status='awaiting_validation',updated_at=utc(),completed=len(result),remaining=0 if args.mode=='full' else 408,
        summary='24 default/checkpoint endpoints match frozen compact; awaiting parent review for full408.' if args.mode=='pilot' else '408 terminal results reevaluated; awaiting independent validation.',
        evidence=[str(dest/'validation.json'),str(dest/'terminal_metrics.csv.gz'),str(HERE/'freeze.json')])
    write_json(OUT/'status.json',status)
    print(json.dumps(validation,indent=2),flush=True)


if __name__=='__main__':
    main()
