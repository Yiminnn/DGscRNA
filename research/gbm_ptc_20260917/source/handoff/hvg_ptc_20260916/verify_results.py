#!/usr/bin/env python3
"""Independent count-based audit and legacy partition parity; no model fitting."""
import argparse
import json
import sys
from pathlib import Path
from common import OUT,ROOT,INPUTS,samples,require_slurm,sha,utc,write_json


def verify(partial=False):
    require_slurm()
    import numpy as np
    import pandas as pd
    from sklearn.metrics import adjusted_rand_score
    from evaluate import L1
    from resolve import structural_resolution
    out=OUT/'verification';out.mkdir(parents=True,exist_ok=True)
    protocol=json.loads((OUT/'protocol/geometries.json').read_text())
    expected={(g['geometry_id'],a['arm_id']) for g in protocol for a in g['arms']}
    assert len(expected)==330
    checked=[];baseline=[];gates=[];missing=[];histories=[];optimization=[]
    def optimization_record(meta,**identity):
        warnings=meta.get('warnings',[])
        n_iter=meta.get('n_iter');maximum=meta.get('params',{}).get('max_iter')
        optimization.append(dict(**identity,n_iter=n_iter,max_iter=maximum,
            converged=meta.get('converged'),
            convergence_warning=any('converg' in text.lower() for text in warnings),
            iteration_budget_reached=n_iter is not None and maximum is not None and n_iter>=maximum,
            warning_count=len(warnings),warnings=' | '.join(warnings)))
    for sample in samples():
        ep=OUT/'evaluation'/sample
        if not (ep/'COMPLETE').exists():
            missing.append(sample);continue
        em=json.loads((ep/'manifest.json').read_text())
        assert sha(ep/'manifest.json')==(ep/'COMPLETE').read_text().strip()
        assert sha(ep/'metrics.csv')==em['outputs']['metrics.csv']
        assert sha(ep/'terminal_L1_confusions.csv.gz')==em['outputs']['terminal_L1_confusions.csv.gz']
        table=pd.read_csv(ep/'metrics.csv').set_index('condition')
        counts=pd.read_csv(ep/'terminal_L1_confusions.csv.gz')
        assert set(tuple(c.split('/')) for c in table.index)==expected
        pm=json.loads((OUT/'prepared'/sample/'manifest.json').read_text())
        n=pm['n_cells']
        for cid,c in counts.groupby('condition'):
            assert c.n.sum()==n and (c.n>0).all()
            tp=c[c.gold.eq(c.prediction)].groupby('gold').n.sum().reindex(L1,fill_value=0)
            gold=c.groupby('gold').n.sum().reindex(L1,fill_value=0)
            pred=c.groupby('prediction').n.sum().reindex(L1,fill_value=0)
            f1=(2*tp/(gold+pred).replace(0,np.nan)).fillna(0)
            present=gold.gt(0)
            known=~c.prediction.isin(['Unknown','Undecided','Noise','nan','None',''])
            actual=dict(terminal_strict_L1_macroF1_fixed11=float(f1.mean()),
                        terminal_strict_L1_macroF1_present=float(f1[present].mean()),
                        terminal_strict_L1_accuracy=float(tp.sum()/n),terminal_coverage=float(c.loc[known,'n'].sum()/n))
            for metric,value in actual.items():
                assert np.isclose(table.loc[cid,metric],value,atol=1e-12,rtol=1e-12),(sample,cid,metric)
            checked.append(dict(sample=sample,condition=cid,n_cells=n,count_metrics_verified=True))
        for g in protocol:
            embedding_manifest=OUT/'fits'/sample/g['geometry_id']/'embedding_manifest.json'
            if embedding_manifest.exists():
                gm=json.loads(embedding_manifest.read_text())
                optimization_record(gm['reducer_info'],sample=sample,feature=g['feature'],
                    geometry_id=g['geometry_id'],arm_id=None,stage='representation',method=g['dr'],seed=g['seed'])
            for a in g['arms']:
                cid=g['geometry_id']+'/'+a['arm_id'];ap=OUT/'fits'/sample/g['geometry_id']/a['arm_id']
                row=table.loc[cid]
                am=json.loads((ap/'manifest.json').read_text()) if (ap/'manifest.json').exists() else None
                if am and 'clustering' in am:
                    optimization_record(am['clustering'],sample=sample,feature=g['feature'],
                        geometry_id=g['geometry_id'],arm_id=a['arm_id'],stage='clustering',method=a['clusterer'],seed=g['seed'])
                # A valid partition is independently auditable even when the
                # downstream legacy DEG scorer cannot annotate singleton groups.
                cp=ap/'clusters.npy'
                if g['feature']=='all' and 'E1_primary_table' in a['families'] and cp.exists():
                    old=ROOT/'results/g274_cohort/fit_v1'/sample/'arms'/f'{g["dr"]}__{a["clusterer"]}'/'labels.npy'
                    if old.exists():
                        prior=np.load(old,allow_pickle=False)
                        current=np.load(cp,allow_pickle=False)
                        assert prior.shape==current.shape==(n,)
                        baseline.append(dict(sample=sample,dr=g['dr'],clusterer=a['clusterer'],
                            annotation_status=row.status,
                            exact_labels=bool(np.array_equal(prior,current)),
                            partition_ari=float(adjusted_rand_score(prior,current)),
                            same_scaled_matrix=pm['old_allgene_scaled_match']))
                if row.status!='completed':
                    r=structural_resolution(sample,g,a)
                    assert r is not None and r['status']==row.status
                    assert pd.isna(row.get('terminal_strict_L1_macroF1_present'))
                    continue
                assert am is not None
                assert sha(ap/'manifest.json')==(ap/'COMPLETE').read_text().strip()
                p=np.load(ap/'predictions.npz',allow_pickle=False)
                ai=am['annotation'];cl=p['cluster'];noise=cl==-1
                assert p['final'].dtype.kind=='U' and p['lineage'].dtype.kind=='U'
                assert (p['final'][noise]=='Unknown').all()
                assert ai['actual_dl_feature_width']==pm['n_genes']
                if ai['training_executed']:
                    hist=json.loads((ap/'training_history.json').read_text())
                    assert hist['input_dim']==pm['n_genes'] and hist['num_classes']==ai['n_training_classes']
                    assert len(hist['train_losses'])==15
                    histories.append(dict(sample=sample,condition=cid,verified_full_width=True,
                                          input_dim=hist['input_dim'],n_epochs=15))
        gates.append(dict(sample=sample,expected=330,n_terminal=int(table.status.eq('completed').sum()),
                          n_structural=int(table.status.str.startswith('structural').sum()),
                          n_numerical=int(table.status.str.startswith('numerical').sum())))
    pd.DataFrame(checked).to_csv(out/'independent_count_metrics.csv.gz',index=False)
    pd.DataFrame(histories).to_csv(out/'actual_DL_width_and_epochs.csv.gz',index=False)
    pd.DataFrame(baseline).to_csv(out/'legacy_allgene_partition_parity.csv',index=False)
    pd.DataFrame(gates).to_csv(out/'sample_gates.csv',index=False)
    pd.DataFrame(optimization).to_csv(out/'optimization_diagnostics.csv.gz',index=False)
    if missing and not partial:raise RuntimeError(f'Missing sample evaluation: {len(missing)}')
    status=dict(timestamp=utc(),status='complete' if not missing else 'partial',n_samples=len(gates),
                n_terminal_conditions_checked=len(checked),n_actual_DL_checked=len(histories),missing_samples=missing,
                verification_source_sha256=sha(Path(__file__)),
                outputs={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name not in ['independent_manifest.json','COMPLETE']})
    write_json(out/'independent_manifest.json',status)
    if not missing:(out/'COMPLETE').write_text(sha(out/'independent_manifest.json')+'\n')
    print(json.dumps({k:v for k,v in status.items() if k not in ['outputs','missing_samples']}),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--partial',action='store_true');verify(p.parse_args().partial)
