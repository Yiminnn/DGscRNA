"""Original PTC refinement with isolated model-seed controls and audited caching."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import uuid
from common import checked, sha, write_json, complete, utc
from ptc_followup_common import PTC, require_ptc
import refine_controls as refine
from terminal import preserve_equivalent_duplicate

def fit(source,arm_id,dest,seed=42,parity_reference=None):
    require_ptc()
    import numpy as np
    import pandas as pd
    import torch
    source=Path(source);dest=Path(dest);dest.mkdir(parents=True,exist_ok=True)
    assert checked(source,'score_manifest.json','SCORE_COMPLETE')
    m=json.loads((source/'score_manifest.json').read_text());arm=m['arms'][arm_id]
    if checked(dest,'terminal_manifest.json','TERMINAL_COMPLETE'):
        tm=json.loads((dest/'terminal_manifest.json').read_text())
        assert tm['score_manifest_sha256']==sha(source/'score_manifest.json') and tm['model_seed']==seed
        return tm
    assert m['initial_sha256']==sha(source/'initial_calls.csv.gz')
    cells=pd.read_csv(source/'cells.csv',dtype=str,keep_default_na=False)
    initial=pd.read_csv(source/'initial_calls.csv.gz',usecols=['cell_id',arm['seed_column']],dtype=str,keep_default_na=False)
    assert np.array_equal(initial.cell_id,cells.cell_id)
    binary=Path(m['DL_binary']);width=int(m['DL_features']);dlhash=m['DL_binary_sha256']
    assert binary.stat().st_size==len(cells)*width*4 and sha(binary)==dlhash
    x=np.memmap(binary,mode='r',dtype='<f4',shape=(len(cells),width))
    assert np.isfinite(x).all()
    refine.PARAMS.update(input=m.get('DL_input_description',m.get('geometry','original PTC CCA expression')),
                         architecture=[256,128],epochs=10,model_seed=int(seed),split_seed=42)
    labels=initial[arm['seed_column']].to_numpy(dtype=str)
    signature=dict(DL_sha256=dlhash,n_cells=len(cells),n_features=width,
        cell_order_sha256=hashlib.sha256('\n'.join(cells.cell_id).encode()).hexdigest(),
        params=dict(refine.PARAMS),training_source_sha256=sha(Path(refine.__file__)),
        torch_version=torch.__version__,torch_threads=torch.get_num_threads())
    key=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()+b'\0'+json.dumps(labels.tolist(),ensure_ascii=False).encode()).hexdigest()
    caches=PTC/'DL_cache';caches.mkdir(parents=True,exist_ok=True);cache=caches/key
    with (caches/(key+'.lock')).open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        reused=checked(cache,'training_manifest.json','COMPLETE')
        if not reused:
            tmp=caches/(key+'.partial.'+uuid.uuid4().hex);tmp.mkdir()
            refine.train_cache(x,labels,tmp,dict(cache_key=key,input_signature=signature,
                first_condition=str(dest),reference_labels_used_for_fit=False,
                historical_model_weights_recovered=False,driver_source_sha256=sha(__file__)))
            try:tmp.rename(cache)
            except FileExistsError:
                preserve_equivalent_duplicate(tmp,cache);reused=True
        tm=json.loads((cache/'training_manifest.json').read_text())
        assert tm['provenance']['input_signature']==signature
        with np.load(cache/'terminal.npz',allow_pickle=False) as z:
            assert np.array_equal(z['initial'],labels)
            known=labels!='Undecided';pool=~known
            for name in ['final090','final070']:assert np.array_equal(z[name][known],labels[known])
            if tm['training_executed']:
                p=z['probabilities'];np.testing.assert_allclose(p.sum(1),1,rtol=1e-5,atol=1e-6)
                assert np.isfinite(p).all() and len(p)==pool.sum()
                confidence=np.asarray([round(v,4) for v in p.max(1)],dtype=np.float32)
                raw=z['classes'][p.argmax(1)]
                for stage,cut in [('final090',.9),('final070',.7)]:
                    assert np.array_equal(z[stage][pool],np.where(confidence>=cut,raw,'Unknown'))
            parity=None
            if parity_reference is not None:
                original=Path(parity_reference)
                assert seed==42 and (original/'TERMINAL_COMPLETE').exists()
                with np.load(original/'terminal.npz',allow_pickle=False) as ref:
                    for field in ['initial','final090','final070','classes','pool_indices','train_indices','validation_indices']:
                        assert np.array_equal(z[field],ref[field]), (str(original),field)
                    np.testing.assert_allclose(z['confidence_rounded'],ref['confidence_rounded'],rtol=0,atol=0,equal_nan=True)
                    np.testing.assert_allclose(z['probabilities'],ref['probabilities'],rtol=1e-6,atol=1e-7)
                    parity=dict(reference=str(original),terminal_labels_splits_rounded_confidence_exact=True,
                        probabilities_bitwise_equal=bool(np.array_equal(z['probabilities'],ref['probabilities'])),
                        probability_rtol=1e-6,probability_atol=1e-7)
            result=cells[['cell_id']].copy()
            for name in ['initial','final090','final070','lineage','confidence_rounded']:result[name]=z[name]
        for name in ['terminal.npz','training_manifest.json','training_history.json','model_state.pt']:
            if (cache/name).exists() and not (dest/name).exists():os.link(cache/name,dest/name)
        path=dest/('predictions.csv.gz.part.'+os.environ['SLURM_JOB_ID']);result.to_csv(path,index=False,compression='gzip');path.replace(dest/'predictions.csv.gz')
        manifest=dict(status='completed',arm=arm,model_seed=int(seed),split_seed=42,
            score_manifest_sha256=sha(source/'score_manifest.json'),source=str(source),
            cache_key=key,cache_directory=str(cache),identical_result_reused=reused,
            training_manifest_sha256=sha(cache/'training_manifest.json'),terminal_sha256=sha(cache/'terminal.npz'),
            predictions_sha256=sha(dest/'predictions.csv.gz'),DL_sha256=dlhash,DL_features=width,
            dl_status=tm['dl_status'],training_executed=tm['training_executed'],
            n_known=tm['n_known'],n_pool=tm['n_pool'],n_training_classes=tm['n_training_classes'],
            default_parity=parity,reference_labels_used_for_fit=False,job=os.environ['SLURM_JOB_ID'],completed_at=utc())
        write_json(dest/'terminal_manifest.json',manifest);complete(dest,'terminal_manifest.json','TERMINAL_COMPLETE')
        print('PTC_TERMINAL',str(dest),tm['dl_status'],'seed',seed,flush=True)
        return manifest
