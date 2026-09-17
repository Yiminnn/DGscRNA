#!/usr/bin/env python3
"""Verify the secondary singleton policy independently; never overwrite primary tables."""
import argparse
import json
from pathlib import Path
from common import ROOT,OUT,samples,sha,require_slurm,utc,write_json


def run(partial=False):
    require_slurm()
    import numpy as np
    import pandas as pd
    from evaluate import L1
    done=[];missing=[];merged=[]
    target=OUT/'singleton_robustness_summary';target.mkdir(parents=True,exist_ok=True)
    for sample in samples():
        p=OUT/'singleton_robustness'/sample
        if not (p/'COMPLETE').exists():missing.append(sample);continue
        m=json.loads((p/'manifest.json').read_text())
        assert m['status']=='completed' and sha(p/'manifest.json')==(p/'COMPLETE').read_text().strip()
        for name,h in m['outputs'].items():assert sha(p/name)==h
        primary=pd.read_csv(OUT/'evaluation'/sample/'metrics.csv').set_index('condition')
        secondary=pd.read_csv(p/'metrics.csv').set_index('condition')
        counts=pd.read_csv(p/'confusions.csv.gz')
        expected=primary.index[primary.status.eq('structural_annotation_failure')]
        assert set(secondary.index)==set(expected) and len(secondary)==m['expected_repairs']
        primary['primary_status']=primary.status
        primary['secondary_policy']='unchanged_primary'
        for col in ['training_executed','final_valid','terminal_valid']:
            if col in primary:primary[col]=primary[col].astype('boolean')
        for cid,row in secondary.iterrows():
            ap=p/Path(cid);state=json.loads((ap/'manifest.json').read_text())
            old=OUT/'fits'/sample/Path(cid)
            assert state['provenance']['primary_manifest_sha256']==sha(old/'manifest.json')
            assert state['provenance']['primary_clusters_sha256']==sha(old/'clusters.npy')
            assert state['prediction_sha256']==sha(ap/'predictions.npz')
            pred=np.load(ap/'predictions.npz',allow_pickle=False)
            cl=pred['cluster'];seed=pred['seed'];final=pred['final']
            np.testing.assert_array_equal(cl,np.load(old/'clusters.npy',allow_pickle=False))
            ids,sizes=np.unique(cl[cl!=-1],return_counts=True);singletons=ids[sizes==1]
            assert (seed[np.isin(cl,singletons)]=='Undecided').all()
            assert (final[cl==-1]=='Unknown').all()
            known=(cl!=-1)&~np.isin(seed,['Unknown','Undecided','Noise','UNMAPPABLE','nan','None',''])
            np.testing.assert_array_equal(final[known],seed[known])
            if row.status=='completed':
                c=counts[counts.condition.eq(cid)]
                assert c.n.sum()==len(cl)
                actual=pd.DataFrame({'gold':[], 'prediction':[], 'n':[]})
                tp=c[c.gold.eq(c.prediction)].groupby('gold').n.sum().reindex(L1,fill_value=0)
                gold=c.groupby('gold').n.sum().reindex(L1,fill_value=0)
                calls=c.groupby('prediction').n.sum().reindex(L1,fill_value=0)
                f1=(2*tp/(gold+calls).replace(0,np.nan)).fillna(0)
                assert np.isclose(row.terminal_strict_L1_macroF1_present,f1[gold>0].mean(),atol=1e-12)
                assert np.isclose(row.terminal_strict_L1_macroF1_fixed11,f1.mean(),atol=1e-12)
                assert np.isclose(row.terminal_strict_L1_accuracy,tp.sum()/len(cl),atol=1e-12)
                if state['annotation']['training_executed']:
                    hist=json.loads((ap/'training_history.json').read_text())
                    assert hist['input_dim']==state['actual_DL_width'] and len(hist['train_losses'])==15
            elif row.status=='structural_DL_unavailable':
                evidence=state['training_split_evidence']
                labels,nlabels=np.unique(seed[known],return_counts=True)
                assert evidence['n_training']==known.sum() and evidence['n_classes']==len(labels)
                assert evidence['n_test']==int(np.ceil(.1*known.sum()))
                assert evidence['n_test']<len(labels) or known.sum()-evidence['n_test']<len(labels) or nlabels.min()<2
                assert pd.isna(row.get('terminal_strict_L1_macroF1_present'))
            else:raise AssertionError((sample,cid,row.status))
            primary.loc[cid,'secondary_policy']='singleton_seed_Undecided'
            primary.loc[cid,'status']=row.status
            for col in secondary.columns:
                if col.startswith('terminal_') or col in ['dl_status','training_executed','n_dl_assigned','final_valid']:
                    primary.loc[cid,col]=row.get(col)
            if row.status!='completed':
                primary.loc[cid,['training_executed','final_valid','terminal_valid']]=False
        merged.append(primary.reset_index())
        done.append(dict(sample=sample,n_repairs=len(secondary),n_success=int(secondary.status.eq('completed').sum()),
                         n_structural_DL=int(secondary.status.eq('structural_DL_unavailable').sum())))
    pd.DataFrame(done).to_csv(target/'verification_counts.csv',index=False)
    if merged:
        data=pd.concat(merged,ignore_index=True)
        data.to_csv(target/'secondary_all_conditions.csv.gz',index=False)
        primary=data[data.evaluable.eq(True)&data.families.str.contains('E1_primary_table')]
        metric='terminal_strict_L1_macroF1_present'
        rows=[]
        for cid,d in primary.groupby('condition'):
            patient=d.groupby('patient')[metric].mean()
            first=d.iloc[0]
            rows.append(dict(condition=cid,feature=first.feature,dr=first.dr,clusterer=first.clusterer,
                patient_mean=patient.mean(),n_available_samples=int(d[metric].notna().sum()),
                n_expected_samples=len(d),n_available_patients=int(patient.notna().sum()),
                n_secondary_repairs=int(d.secondary_policy.eq('singleton_seed_Undecided').sum()),
                fixed_cohort_lower=d.assign(value=d[metric].fillna(0)).groupby('patient').value.mean().mean(),
                fixed_cohort_upper=d.assign(value=d[metric].fillna(1)).groupby('patient').value.mean().mean()))
        pd.DataFrame(rows).to_csv(target/'secondary_primary_factorial.csv',index=False)
    state=dict(status='completed' if not missing else 'partial',timestamp=utc(),n_samples=len(done),
        n_repairs_verified=sum(r['n_repairs'] for r in done),missing_samples=missing,source_sha256=sha(Path(__file__)),
        outputs={p.name:sha(p) for p in target.iterdir() if p.is_file() and p.name not in ['manifest.json','COMPLETE']})
    write_json(target/'manifest.json',state)
    if not missing:(target/'COMPLETE').write_text(sha(target/'manifest.json')+'\n')
    print({k:v for k,v in state.items() if k not in ['outputs','missing_samples']},flush=True)
    if missing and not partial:raise RuntimeError(f'Secondary policy missing for {len(missing)} samples')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--partial',action='store_true');run(p.parse_args().partial)
