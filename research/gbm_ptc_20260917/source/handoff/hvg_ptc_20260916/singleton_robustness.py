#!/usr/bin/env python3
"""Explicit secondary scoring policy: unsupported singleton clusters seed Undecided.

Frozen primary fits are never modified. Only verified singleton-DEG failures are eligible.
No cells, markers, embedding, cluster IDs, DL features or hyperparameters are changed.
"""
import argparse
from contextlib import redirect_stdout,redirect_stderr
import json
from pathlib import Path
import sys
import time
import traceback
from common import ROOT,OUT,MARKERS,MARKER,require_slurm,sha,utc,write_json,samples,runtime_record


def run(sample):
    require_slurm()
    import numpy as np
    import pandas as pd
    import anndata as ad
    import scanpy as sc
    import torch
    from threadpoolctl import threadpool_limits
    torch.set_num_threads(4);torch.set_num_interop_threads(1);threadpool_limits(8)
    sys.path.insert(0,str(ROOT/'handoff/g274_table4'))
    from fit import prepared
    from resolve import structural_resolution
    from cgo_grid import load_panels,score_panels,calls_at_cutoff,seed_labels,final_prediction
    from dgscrna.core.deep_learning import train_deep_model
    from darmanis_region_final import ABSTAIN
    source_hash=sha(Path(__file__))
    source_eval=OUT/'evaluation'/sample
    assert (source_eval/'COMPLETE').exists(), 'Wait for the complete frozen-primary sample evaluation'
    out=OUT/'singleton_robustness'/sample;out.mkdir(parents=True,exist_ok=True)
    protocol=json.loads((OUT/'protocol/geometries.json').read_text())
    todo=[]
    for g in protocol:
        for arm in g['arms']:
            primary=OUT/'fits'/sample/g['geometry_id']/arm['arm_id']
            if (primary/'COMPLETE').exists():continue
            r=structural_resolution(sample,g,arm)
            if r and r['reason']=='legacy_DEG_cannot_score_singleton_clusters':todo.append((g,arm,r))
    a,pm,prep=prepared(sample)
    mapping=pd.read_csv(MARKERS/'mapping_L1_v3.csv')
    mapping=mapping[mapping.marker_set.eq(MARKER)].set_index('panel')['gold_class'].to_dict()
    outputs=[]
    for g,arm,resolution in todo:
        dest=out/g['geometry_id']/arm['arm_id'];dest.mkdir(parents=True,exist_ok=True)
        primary=OUT/'fits'/sample/g['geometry_id']/arm['arm_id']
        provenance=dict(source_script_sha256=source_hash,primary_manifest_sha256=sha(primary/'manifest.json'),
                        primary_clusters_sha256=sha(primary/'clusters.npy'),prepared_manifest_sha256=sha(prep/'manifest.json'))
        if (dest/'COMPLETE').exists():
            previous=json.loads((dest/'manifest.json').read_text())
            assert previous['provenance']==provenance
            assert sha(dest/'predictions.npz')==previous['prediction_sha256']
            outputs.append(previous);continue
        started=time.perf_counter()
        state=dict(sample=sample,geometry_id=g['geometry_id'],arm_id=arm['arm_id'],provenance=provenance,
                   status='running',policy='secondary_singleton_seed_Undecided',gold_used_for_fitting=False,
                   started_at=utc(),**runtime_record())
        write_json(dest/'manifest.json',state)
        try:
            cl_ids=np.load(primary/'clusters.npy',allow_pickle=False)
            cl=np.asarray(['Noise' if x==-1 else str(x) for x in cl_ids])
            counts=pd.Series(cl).value_counts()
            all_groups=[c for c in pd.unique(cl) if c!='Noise']
            supported=[c for c in all_groups if counts[c]>=2]
            unsupported=[c for c in all_groups if counts[c]==1]
            assert len(unsupported)>0 and supported
            idx=np.load(prep/f'indices_{g["feature"]}.npy',allow_pickle=False)
            scoring_idx=np.arange(a.n_vars) if arm['scoring_features']=='all' else idx
            b=ad.AnnData(a.layers['lognorm'][:,scoring_idx].copy())
            b.obs_names=a.obs_names;b.var_names=a.var_names[scoring_idx]
            b.obs['_cl']=pd.Categorical(cl);b.layers['lognorm']=b.X
            # Singleton cells remain in the reference rest and in the final prediction roster.
            sc.tl.rank_genes_groups(b,'_cl',groups=supported,method='wilcoxon',n_genes=100,layer='lognorm',use_raw=False)
            deg={c:dict(zip(b.uns['rank_genes_groups']['names'][c],b.uns['rank_genes_groups']['logfoldchanges'][c])) for c in supported}
            S,names=score_panels(deg,load_panels(b),supported)
            pd.DataFrame(S,index=names,columns=supported).to_csv(dest/'supported_panel_scores.csv.gz')
            calls=dict(zip(supported,calls_at_cutoff(S,names,arm['cutoff'])))
            calls.update({c:'Undecided' for c in unsupported})
            seed=seed_labels(cl,calls,mapping)
            assert all(seed[i]=='Undecided' for i,c in enumerate(cl) if c in unsupported)
            def traced_train(sub,key,**kwargs):
                kwargs['random_state']=g['seed']
                model=train_deep_model(sub,key,**kwargs)
                write_json(dest/'training_history.json',{k:v for k,v in model.items() if k!='model'})
                return model
            with (dest/'refinement.log').open('w') as log,redirect_stdout(log),redirect_stderr(log):
                final,lineage,info=final_prediction(a,cl,seed,train_fn=traced_train)
            outcome='completed' if info['dl_status']!='dl_error' else 'failed'
            if outcome=='failed':
                training_seed=seed[(cl_ids!=-1)&~np.isin(seed,list(ABSTAIN))]
                labels,label_counts=np.unique(training_seed,return_counts=True)
                n_test=int(np.ceil(.1*len(training_seed)))
                split_impossible=(n_test<len(labels) or len(training_seed)-n_test<len(labels)
                                  or (len(label_counts)>0 and label_counts.min()<2))
                message=info.get('error','')
                if split_impossible and ('test_size' in message or 'least populated class' in message):
                    outcome='structural_DL_unavailable'
                    state['unsupported_reason']='fixed_stratified_DL_split_not_feasible'
                    state['training_split_evidence']=dict(n_training=len(training_seed),n_test=n_test,
                        n_classes=len(labels),label_counts=dict(zip(labels,label_counts.tolist())),error=message)
            final=np.asarray(final,dtype=str);lineage=np.asarray(lineage,dtype=str)
            assert len(final)==pm['n_cells'] and (final[cl_ids==-1]=='Unknown').all()
            np.savez_compressed(dest/'predictions.npz',cluster=cl_ids,seed=np.asarray(seed,dtype=str),final=final,lineage=lineage)
            state.update(status=outcome,
                annotation=info,n_singleton_cells=len(unsupported),supported_groups=supported,unsupported_groups=unsupported,
                actual_scoring_width=len(scoring_idx),actual_DL_width=a.n_vars,
                seconds=time.perf_counter()-started,prediction_sha256=sha(dest/'predictions.npz'),finished_at=utc())
            write_json(dest/'manifest.json',state)
            if state['status']!='failed':(dest/'COMPLETE').write_text(sha(dest/'manifest.json')+'\n')
            outputs.append(state)
            print(sample,g['feature'],g['dr'],arm['clusterer'],len(unsupported),info['dl_status'],flush=True)
        except Exception as exc:
            state.update(status='failed',error=str(exc),traceback=traceback.format_exc(),finished_at=utc())
            write_json(dest/'manifest.json',state);outputs.append(state)
            print(state['traceback'],flush=True)
    # Reference labels are loaded only after every repaired terminal prediction is saved.
    from cgo_grid import reference_frame,load_helpers,label_scope,score_prediction
    from sklearn.metrics import precision_recall_fscore_support
    from evaluate import L1
    pc=reference_frame(sample);helper=load_helpers();scope=label_scope(pc,helper)
    assert list(pc.CellID)==list(a.obs_names)
    rows=[];confusions=[]
    for state in outputs:
        dest=out/state['geometry_id']/state['arm_id']
        row=dict(sample=sample,condition=state['geometry_id']+'/'+state['arm_id'],status=state['status'],
                 policy='secondary_singleton_seed_Undecided',n_cells=a.n_obs)
        if state['status']=='completed':
            p=np.load(dest/'predictions.npz',allow_pickle=False);pred=p['final'];truth=pc.L1.to_numpy(dtype=str)
            score=score_prediction(pred,pc,scope,helper)
            pr,re,f1,support=precision_recall_fscore_support(truth,pred,labels=L1,zero_division=0)
            row.update({f'terminal_{k}':v for k,v in score.items()})
            row.update(terminal_strict_L1_macroF1_present=float(f1[support>0].mean()),
                       terminal_strict_L1_macroF1_fixed11=float(f1.mean()),terminal_strict_L1_accuracy=float((pred==truth).mean()),
                       dl_status=state['annotation']['dl_status'],training_executed=state['annotation']['training_executed'],
                       n_dl_assigned=state['annotation']['n_dl_assigned'],n_singleton_cells=state['n_singleton_cells'],
                       final_valid=state['annotation']['final_valid'],terminal_valid=True)
            count=pd.DataFrame({'gold':truth,'prediction':pred}).value_counts().reset_index(name='n')
            confusions.extend(dict(condition=row['condition'],**r) for r in count.to_dict('records'))
        rows.append(row)
    cols=['sample','condition','status','policy','n_cells']
    pd.DataFrame(rows,columns=None if rows else cols).to_csv(out/'metrics.csv',index=False)
    pd.DataFrame(confusions,columns=None if confusions else ['condition','gold','prediction','n']).to_csv(out/'confusions.csv.gz',index=False)
    failures=[r['condition'] for r in rows if r['status']=='failed']
    manifest=dict(sample=sample,status='completed' if not failures else 'failed',expected_repairs=len(todo),
        completed_repairs=sum(r['status']=='completed' for r in rows),
        structural_DL_unavailable=sum(r['status']=='structural_DL_unavailable' for r in rows),
        failed_conditions=failures,source_script_sha256=source_hash,timestamp=utc(),
        primary_evaluation_sha256=sha(source_eval/'manifest.json'),
        outputs={n:sha(out/n) for n in ['metrics.csv','confusions.csv.gz']})
    write_json(out/'manifest.json',manifest)
    if not failures:(out/'COMPLETE').write_text(sha(out/'manifest.json')+'\n')
    if failures:raise RuntimeError(f'{sample}: {len(failures)} repair conditions failed')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--sample');p.add_argument('--index',type=int)
    x=p.parse_args();run(x.sample or samples()[x.index])
