#!/usr/bin/env python3
"""Independent reference evaluation AFTER predictions have been saved; SLURM only."""
import argparse
import json
from pathlib import Path
import sys
from common import ROOT, OUT, INPUTS, MARKER, require_slurm, sha, utc, write_json, samples

L1 = ['Malignant','TAM','Lymphocyte','Oligodendrocyte','Astrocyte','OPC',
      'Excitatory neuron','Inhibitory neuron','Endothel','Pericyte','Other']
ABSTAIN = ['Unknown','Undecided','Noise','nan','None','']


def evaluate(sample, partial=False):
    require_slurm()
    import numpy as np
    import pandas as pd
    from sklearn.metrics import precision_recall_fscore_support, cohen_kappa_score
    sys.path.insert(0, str(ROOT/'handoff/g274_table4'))
    from cgo_grid import reference_frame, load_helpers, label_scope, score_prediction
    from metrics_core import partition_metrics, macro_f1_22
    from resolve import structural_resolution
    out = OUT/'evaluation'/sample
    out.mkdir(parents=True, exist_ok=True)
    pc = reference_frame(sample)
    goldpath = INPUTS/sample/'labels.csv'
    gold = pd.read_csv(goldpath)
    manifest = json.loads((INPUTS/sample/'manifest.json').read_text())
    assert sha(goldpath) == manifest['artifacts']['labels.csv']['sha256']
    rawobs = ROOT/'data_bench/GSE274546/mtx'/sample/'obs.csv'
    assert sha(rawobs) == manifest['sha256'][str(rawobs)]
    cells = (OUT/'prepared'/sample/'cells.tsv').read_text().splitlines()
    assert list(pc.CellID) == list(gold.CellID) == cells
    assert len(set(pc.Patient.astype(str))) == 1
    helper = load_helpers()
    scope = label_scope(pc, helper)
    lfine22 = json.loads((INPUTS/sample/'lfine22.json').read_text())['classes']
    labels = gold.lfine_eval23.to_numpy(dtype=str)
    truth = pc.L1.to_numpy(dtype=str)
    assert set(truth) <= set(L1), f'Unspecified L1 labels: {set(truth)-set(L1)}'
    gs = json.loads((OUT/'protocol/geometries.json').read_text())
    rows, confusions, filehashes, partition_cache, annotation_cache = [], [], {}, {}, {}
    missing, failed, structural = [], [], []
    for g in gs:
        for arm in g['arms']:
            cid = f'{g["geometry_id"]}/{arm["arm_id"]}'
            dest = OUT/'fits'/sample/g['geometry_id']/arm['arm_id']
            row = dict(sample=sample, patient=str(pc.Patient.iloc[0]),
                       evaluable=bool(manifest['evaluable']), condition=cid,
                       geometry_id=g['geometry_id'], arm_id=arm['arm_id'],
                       feature=g['feature'], dr=g['dr'], dim=g['dim'],
                       input_space=g['input_space'], seed=g['seed'], neighbors=g['neighbors'],
                       min_dist=g['min_dist'], **{k:v for k,v in arm.items() if k not in ['arm_id','families']},
                       families=';'.join(arm['families']), n_cells=len(cells))
            if not (dest/'COMPLETE').exists():
                resolution = structural_resolution(sample,g,arm)
                if resolution is not None:
                    structural.append(cid)
                    row.update(status=resolution['status'],error=resolution.get('error',resolution['reason']),
                               final_valid=False,terminal_valid=False,training_executed=False,
                               structural_reason=resolution['reason'])
                    if resolution['partition_available']:
                        cl=np.load(dest/'clusters.npy',allow_pickle=False)
                        pm=partition_metrics(labels,cl)
                        row.update({f'partition_{k}':pm[k] for k in ['pair_f1','ari','fmi','v_measure']})
                        row.update(noise_rate=pm['noise_rate'],n_clusters=pm['n_clusters'])
                        am=json.loads((dest/'manifest.json').read_text())
                        ci=am.get('clustering',{})
                        row.update({f'internal_{k}':v for k,v in ci.get('label_free_selection',{}).items()})
                        row['clustering_seconds']=ci.get('seconds')
                    rows.append(row)
                    continue
                note = json.loads((dest/'manifest.json').read_text()) if (dest/'manifest.json').exists() else {}
                row.update(status=note.get('status','missing'), error=note.get('error'))
                (failed if note.get('status')=='failed' else missing).append(cid)
                rows.append(row)
                continue
            am = json.loads((dest/'manifest.json').read_text())
            assert sha(dest/'manifest.json') == (dest/'COMPLETE').read_text().strip()
            assert am['status'] == 'completed' and not am['gold_used_for_fitting']
            for name,h in am['outputs'].items():
                assert sha(dest/name) == h, f'Result checksum mismatch: {dest/name}'
            filehashes[cid] = sha(dest/'manifest.json')
            p = np.load(dest/'predictions.npz', allow_pickle=False)
            cl, final, initial, lineage = [p[k] for k in ['cluster','final','seed','lineage']]
            assert cl.shape == final.shape == initial.shape == lineage.shape == (len(cells),)
            np.testing.assert_array_equal(cl, np.load(dest/'clusters.npy',allow_pickle=False))
            assert (final[cl==-1] == 'Unknown').all()
            retained = (~np.isin(initial, ABSTAIN)) & (cl!=-1)
            np.testing.assert_array_equal(initial[retained], final[retained])
            ai, ci = am['annotation'], am['clustering']
            assert ai['terminal_valid']
            row.update(status='completed', dl_status=ai['dl_status'], final_valid=ai['final_valid'],
                       terminal_valid=ai['terminal_valid'], training_executed=ai['training_executed'],
                       n_pool=ai['n_pool'], n_training=ai['n_training'],
                       n_training_classes=ai['n_training_classes'], n_dl_assigned=ai['n_dl_assigned'],
                       noise_rate=ci['noise_frac'], n_clusters=ci['n_clusters'],
                       clustering_seconds=ci['seconds'], annotation_seconds=ai['seconds'])
            gm = json.loads((dest.parent/'manifest.json').read_text())
            row['reduction_seconds'] = gm.get('reducer',{}).get('seconds')
            row['geometry_features'] = gm.get('actual_geometry_feature_width')
            row['dl_features'] = ai['actual_dl_feature_width']
            if 'label_free_selection' in ci:
                row.update({f'internal_{k}':v for k,v in ci['label_free_selection'].items()})
            cluster_key = cl.tobytes()
            if cluster_key not in partition_cache:
                pm = partition_metrics(labels, cl)
                aux = macro_f1_22(labels, cl, lfine22)
                partition_cache[cluster_key] = {**{f'partition_{k}':pm[k] for k in ['pair_f1','ari','fmi','v_measure']},
                    'hungarian_auxiliary':aux}
            pm = partition_cache[cluster_key]
            row.update({k:v for k,v in pm.items() if k!='hungarian_auxiliary'})
            aux = pm['hungarian_auxiliary']
            for k in ['macro_f1_22','macro_f1','coverage']:
                if k in aux:
                    row[f'partition_hungarian_{k}'] = aux[k]
            for name,pred in [('terminal',final),('marker_only_ablation',initial)]:
                pk = tuple(pred)
                if pk not in annotation_cache:
                    legacy = score_prediction(pred, pc, scope, helper)
                    strict = precision_recall_fscore_support(truth, pred, labels=L1, zero_division=0)
                    present = strict[3] > 0
                    vals = dict(legacy,
                        lfine_fixed22_setF1=helper['macro_f1_lfine'](pred, scope[0], lfine22, scope[2]),
                        strict_L1_macroF1_fixed11=float(strict[2].mean()),
                        strict_L1_macroF1_present=float(strict[2][present].mean()),
                        strict_L1_accuracy=float((pred==truth).mean()),
                        strict_L1_kappa=float(cohen_kappa_score(truth,pred)),
                        strict_L1_called_accuracy=float((pred[~np.isin(pred, ABSTAIN)]==truth[~np.isin(pred, ABSTAIN)]).mean())
                            if (~np.isin(pred, ABSTAIN)).any() else None)
                    annotation_cache[pk] = vals
                row.update({f'{name}_{k}':v for k,v in annotation_cache[pk].items()})
            # Compact, exact count table supports a second implementation of strict metrics.
            counts = pd.DataFrame({'gold':truth,'prediction':final}).value_counts().reset_index(name='n')
            for record in counts.to_dict('records'):
                confusions.append(dict(condition=cid, **record))
            rows.append(row)
    table = pd.DataFrame(rows)
    assert len(table)==330 and not table.condition.duplicated().any()
    table.to_csv(out/'metrics.csv',index=False)
    pd.DataFrame(confusions).to_csv(out/'terminal_L1_confusions.csv.gz',index=False)
    state = dict(sample=sample, timestamp=utc(), expected=len(table), complete=int(table.status.eq('completed').sum()),
        structural_unavailable=structural,n_structural_unavailable=len(structural),
        missing=missing, failed=failed, status='completed' if not missing and not failed else 'incomplete',
        evaluation_source_sha256=sha(Path(__file__)), gold_sha256=sha(goldpath), rawobs_sha256=sha(rawobs),
        result_manifest_hashes=filehashes,
        outputs={n:sha(out/n) for n in ['metrics.csv','terminal_L1_confusions.csv.gz']})
    write_json(out/'manifest.json',state)
    if state['status']=='completed':
        (out/'COMPLETE').write_text(sha(out/'manifest.json')+'\n')
    print(json.dumps({k:v for k,v in state.items() if k not in ['result_manifest_hashes','missing','failed','structural_unavailable']}),flush=True)
    if state['status']!='completed' and not partial:
        raise RuntimeError(f'{sample}: {len(missing)} missing and {len(failed)} failed conditions')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--sample')
    p.add_argument('--index',type=int)
    p.add_argument('--partial',action='store_true')
    args=p.parse_args()
    evaluate(args.sample or samples()[args.index],args.partial)
