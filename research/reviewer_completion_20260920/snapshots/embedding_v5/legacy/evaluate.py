"""Author labels enter here, after all native-panel predictions have been frozen."""
import argparse
import json
import os
from pathlib import Path
from common import OUT, L1, ROUTES, require_slurm, sha, write_json, complete, checked, utc

def run(prep,partial=False,family='native_R_budget',only_arms=None):
    require_slurm()
    import numpy as np
    import pandas as pd
    from sklearn.metrics import precision_recall_fscore_support, adjusted_rand_score, normalized_mutual_info_score, fowlkes_mallows_score
    prep=Path(prep);pm=json.loads((prep/'prepare_manifest.json').read_text())
    sample=pm['sample'];budget=pm['budget']
    assert family in ['native_R_budget','geometry_only_fixed_DL2000']
    d=prep/('evaluation' if family=='native_R_budget' else 'evaluation_geometry_only_DL2000');d.mkdir(exist_ok=True)
    if checked(d) and not partial:return
    im=json.loads((OUT/'inputs'/sample/'input_manifest.json').read_text())
    ep=OUT/'evaluation_inputs'/sample/'truth.csv.gz'
    assert sha(ep)==im['evaluation_files']['truth.csv.gz']
    truth=pd.read_csv(ep,keep_default_na=False,dtype=str)
    cells=pd.read_csv(prep/'cells.csv',keep_default_na=False,dtype=str)
    assert np.array_equal(cells.cell_id,truth.cell_id)
    mp=pd.read_csv(OUT/'markers/panel_L1_mapping.csv',keep_default_na=False,dtype=str)
    lookups={lib:dict(zip(g.panel,g.L1)) for lib,g in mp.groupby('library')}
    y=truth.L1.to_numpy(dtype=str)
    coarse=lambda x:np.asarray(['Neuron' if v in ['Excitatory neuron','Inhibitory neuron','AMBIGUOUS_NEURON'] else v for v in x])
    coarse_labels=[v for v in L1 if v not in ['Excitatory neuron','Inhibitory neuron']]+['Neuron']
    yc=coarse(y)
    rows=[];classes=[];confusions=[];partitions=[];files={};missing=[]
    for route in ROUTES:
        src=prep/route
        if not checked(src,'score_manifest.json','SCORE_COMPLETE'):
            missing.append(route+'/SCORE');continue
        sm=json.loads((src/'score_manifest.json').read_text())
        cl=pd.read_csv(src/'clusters.csv',keep_default_na=False,dtype=str)
        assert np.array_equal(cl.cell_id,truth.cell_id)
        partitions.append(dict(sample=sample,patient=im['patient'],primary=im['primary'],budget=budget,route=route,
            evaluation_type='clustering',n_cells=len(y),n_clusters=cl.loc[~cl.cluster.eq('0') if 'HDBSCAN' in route else cl.cluster.notna(),'cluster'].nunique(),n_partition_ids=cl.cluster.nunique(),noise_rate=float(cl.cluster.eq('0').mean()) if 'HDBSCAN' in route else 0.,
            L1_ARI=float(adjusted_rand_score(y,cl.cluster)),L1_NMI=float(normalized_mutual_info_score(y,cl.cluster)),L1_FMI=float(fowlkes_mallows_score(y,cl.cluster)),
            fine_ARI=float(adjusted_rand_score(truth.lfine_original,cl.cluster)),
            malignant_state_ARI=float(adjusted_rand_score(truth.loc[truth.L1.eq('Malignant'),'MalState'],cl.loc[truth.L1.eq('Malignant'),'cluster'])) if truth.L1.eq('Malignant').sum()>1 else None))
        roots=[(family,src/('terminal' if family=='native_R_budget' else 'terminal_geometry_only_DL2000'))]
        for family,terminal in roots:
            for aid,arm in sm['arms'].items():
                if only_arms is not None and aid not in only_arms:continue
                td=terminal/aid
                if not checked(td,'terminal_manifest.json','TERMINAL_COMPLETE'):
                    missing.append(f'{route}/{family}/{aid}');continue
                tm=json.loads((td/'terminal_manifest.json').read_text())
                assert tm['score_manifest_sha256']==sha(src/'score_manifest.json')
                assert sha(td/'terminal.npz')==tm['terminal_sha256']
                files[str(td.relative_to(prep))]=sha(td/'terminal_manifest.json')
                z=np.load(td/'terminal.npz',allow_pickle=False)
                known=z['initial']!='Undecided'
                assert np.array_equal(z['initial'][known],z['final090'][known])
                mapping=lookups[arm['library']]
                for stage,key in [('marker_only','initial'),('terminal090','final090'),('terminal070','final070')]:
                    native=z[key]
                    pred=np.asarray([mapping.get(v,'Unknown' if v in ['Unknown','Undecided','Noise',''] else 'UNMAPPABLE') for v in native])
                    precision,recall,f1,support=precision_recall_fscore_support(y,pred,labels=L1,zero_division=0)
                    present=support>0
                    pc=coarse(pred)
                    cp,cr,cf,cs=precision_recall_fscore_support(yc,pc,labels=coarse_labels,zero_division=0)
                    context=dict(sample=sample,patient=im['patient'],primary=im['primary'],budget=budget,route=route,
                        library=arm['library'],cutoff=arm['cutoff'],arm_id=aid,stage=stage,family=family,
                        seed=pm['seed'],n_cells=len(y),DL_features=tm['DL_features'])
                    abstain=np.isin(native,['Unknown','Undecided','Noise',''])
                    mapped=np.isin(pred,L1)
                    rows.append(dict(**context,status='completed',dl_status=tm['dl_status'],training_executed=tm['training_executed'],
                        macroF1_present=float(f1[present].mean()),macroF1_fixed11=float(f1.mean()),
                        weightedF1=float((f1*support).sum()/support.sum()),accuracy=float((pred==y).mean()),coverage=float((~abstain).mean()),
                        unknown_rate=float(abstain.mean()),mapped_coverage=float(mapped.mean()),
                        off_vocabulary_rate=float((~mapped & ~abstain).mean()),
                        coarse10_macroF1_present=float(cf[cs>0].mean()),coarse10_macroF1_fixed10=float(cf.mean()),
                        n_known=tm['n_known'],n_pool=tm['n_pool'],n_training_classes=tm['n_training_classes'],
                        n_new_correct=int(((~known)&(pred==y)).sum()),
                        n_new_incorrect=int(((~known)&(~abstain)&(pred!=y)).sum()),
                        n_initially_wrong_retained=int((known&(pred!=y)).sum()),
                        reference_overlap=arm['library'] in ['CARE_TME','BrainAtlas112','UNION_all']))
                    for i,label in enumerate(L1):
                        classes.append(dict(**context,label=label,precision=float(precision[i]),recall=float(recall[i]),F1=float(f1[i]),support=int(support[i])))
                    ct=pd.DataFrame({'truth':y,'prediction':pred}).value_counts().reset_index(name='n')
                    for item in ct.to_dict('records'):confusions.append(dict(**context,**item))
    pd.DataFrame(rows).to_csv(d/'metrics.csv',index=False)
    pd.DataFrame(classes).to_csv(d/'per_class.csv.gz',index=False,compression='gzip')
    pd.DataFrame(confusions).to_csv(d/'confusions.csv.gz',index=False,compression='gzip')
    pd.DataFrame(partitions).to_csv(d/'clustering.csv',index=False)
    m=dict(status='partial' if missing else 'completed',sample=sample,budget=budget,family=family,only_arms=only_arms,n_metric_rows=len(rows),
        missing=missing,terminal_manifests=files,truth_sha256=sha(ep),mapping_sha256=sha(OUT/'markers/panel_L1_mapping.csv'),
        primary='macro-F1 over present original-author L1 classes; every cell retained and Unknown/unmapped errors',
        secondary='fixed11 and collapsed-neuron10, reported separately without choosing granularity from ranking',
        malignant_states='Partition ARI only; no claim of state-specific terminal annotation',
        source_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],completed_at=utc(),
        outputs={n:sha(d/n) for n in ['metrics.csv','per_class.csv.gz','confusions.csv.gz','clustering.csv']})
    write_json(d/'manifest.json',m)
    if not missing:complete(d)
    elif not partial:raise RuntimeError(f'{sample}/{budget}: {len(missing)} missing terminal conditions')
    print('EVALUATED',sample,budget,len(rows),'missing',len(missing),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('prep');p.add_argument('--partial',action='store_true');a=p.parse_args()
    run(a.prep,a.partial)
