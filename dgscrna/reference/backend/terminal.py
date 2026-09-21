"""Native-panel terminal DL with exact-input caching and independent output checks."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import sys
import uuid
try:
    from .common import OUT, require_slurm, sha, write_json, complete, checked, utc, execution_id
    from . import legacy_refine as refine
except ImportError:
    from common import OUT, require_slurm, sha, write_json, complete, checked, utc, execution_id
    import legacy_refine as refine

def preserve_equivalent_duplicate(temp,cache):
    """GPFS flock can be node-local: verify before accepting atomic publication."""
    require_slurm()
    import numpy as np
    temp=Path(temp);cache=Path(cache)
    assert checked(temp,'training_manifest.json','COMPLETE')
    assert checked(cache,'training_manifest.json','COMPLETE')
    left=json.loads((temp/'training_manifest.json').read_text())
    right=json.loads((cache/'training_manifest.json').read_text())
    assert left['provenance']['cache_key']==right['provenance']['cache_key']
    assert left['provenance']['input_signature']==right['provenance']['input_signature']
    for path,manifest in [(temp,left),(cache,right)]:
        assert sha(path/'terminal.npz')==manifest['outputs']['terminal.npz']
    with np.load(temp/'terminal.npz',allow_pickle=False) as a,np.load(cache/'terminal.npz',allow_pickle=False) as b:
        assert set(a.files)==set(b.files)
        for name in a.files:
            if name=='probabilities':np.testing.assert_allclose(a[name],b[name],rtol=1e-6,atol=1e-7)
            elif name=='confidence_rounded':np.testing.assert_allclose(a[name],b[name],rtol=0,atol=0,equal_nan=True)
            else:assert np.array_equal(a[name],b[name]),f'Concurrent deterministic result differs: {name}'
    write_json(temp/'DUPLICATE_TERMINAL_EQUIVALENCE.json',dict(status='passed',published_cache=str(cache),
        terminal_labels_splits_rounded_confidence_exact=True,probability_rtol=1e-6,probability_atol=1e-7,
        published_manifest_sha256=sha(cache/'training_manifest.json'),duplicate_manifest_sha256=sha(temp/'training_manifest.json'),
        note='Verified duplicate retained; atomic canonical publication reused',job=execution_id()))
    directory=cache.parent/'duplicate_publications';directory.mkdir(exist_ok=True)
    temp.rename(directory/temp.name)

def finish_route(source, only_arm=None, dl_prep=None):
    require_slurm()
    import numpy as np
    import pandas as pd
    import torch
    source=Path(source)
    assert checked(source,'score_manifest.json','SCORE_COMPLETE')
    m=json.loads((source/'score_manifest.json').read_text())
    assert sha(source/'initial_calls.csv.gz')==m['initial_sha256']
    initial=pd.read_csv(source/'initial_calls.csv.gz',keep_default_na=False,dtype=str)
    cells=pd.read_csv(source/'cells.csv',keep_default_na=False,dtype=str)
    assert np.array_equal(initial.cell_id,cells.cell_id)
    binary=Path(m['DL_binary']);width=int(m['DL_features'])
    dlhash=m['DL_binary_sha256']
    family='native_R_budget'
    if dl_prep is not None:
        dl_prep=Path(dl_prep)
        dm=json.loads((dl_prep/'prepare_manifest.json').read_text())
        assert np.array_equal(pd.read_csv(dl_prep/'cells.csv',keep_default_na=False,dtype=str).cell_id,cells.cell_id)
        binary=Path(dm['DL_binary']);width=int(dm['features']['DL']);dlhash=dm['DL_binary_sha256']
        family='geometry_only_fixed_DL2000'
    assert binary.stat().st_size==len(cells)*width*4 and sha(binary)==dlhash
    x=np.memmap(binary,mode='r',dtype='<f4',shape=(len(cells),width))
    assert np.isfinite(x).all()
    refine.PARAMS['input']=('Native R integrated expression, selected genes in recorded order; not scaled expression'
                            if m.get('assay')=='integrated' else
                            'Native R normalized RNA, selected genes in recorded order; not scaled expression')
    source_hash=sha(Path(refine.__file__))
    signature=dict(DL_sha256=dlhash,n_cells=len(cells),n_features=width,
        cell_order_sha256=hashlib.sha256('\n'.join(cells.cell_id).encode()).hexdigest(),
        training_source_sha256=source_hash,params=refine.PARAMS,
        torch_version=torch.__version__,torch_threads=torch.get_num_threads())
    signature_bytes=json.dumps(signature,sort_keys=True).encode()
    cache_root=Path(os.environ.get('DGSCRNA_DL_CACHE_ROOT',str(OUT/'DL_cache')))
    assert cache_root.is_relative_to(OUT),'DL cache override must stay inside this result campaign'
    cache_root.mkdir(parents=True,exist_ok=True)
    terminal_root=source/('terminal' if dl_prep is None else 'terminal_geometry_only_DL2000')
    terminal_root.mkdir(exist_ok=True)
    arms=[only_arm] if only_arm else list(m['arms'])
    for aid in arms:
        arm=m['arms'][aid];dest=terminal_root/aid;dest.mkdir(exist_ok=True)
        if checked(dest,'terminal_manifest.json','TERMINAL_COMPLETE'):
            old=json.loads((dest/'terminal_manifest.json').read_text())
            assert old['score_manifest_sha256']==sha(source/'score_manifest.json')
            assert old['DL_sha256']==dlhash
            continue
        labels=initial[arm['seed_column']].to_numpy(dtype=str)
        key=hashlib.sha256(signature_bytes+b'\0'+json.dumps(labels.tolist(),ensure_ascii=False).encode()).hexdigest()
        cache=cache_root/key
        with (cache_root/(key+'.lock')).open('a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            reused=checked(cache,'training_manifest.json','COMPLETE')
            if not reused:
                tmp=cache_root/(key+'.partial.'+uuid.uuid4().hex)
                tmp.mkdir()
                refine.train_cache(x,labels,tmp,dict(cache_key=key,input_signature=signature,
                    first_condition=str(dest),reference_labels_used_for_fit=False,
                    driver_source_sha256=sha(__file__),job=execution_id()))
                assert checked(tmp,'training_manifest.json','COMPLETE')
                try:
                    tmp.rename(cache)
                except FileExistsError:
                    preserve_equivalent_duplicate(tmp,cache)
                    reused=True
            cm=json.loads((cache/'training_manifest.json').read_text())
            assert cm['provenance']['cache_key']==key
            z=np.load(cache/'terminal.npz',allow_pickle=False)
            assert np.array_equal(z['initial'],labels)
            known=labels!='Undecided';pool=~known
            for name in ['final090','final070']:
                assert z[name].shape==labels.shape and np.array_equal(z[name][known],labels[known])
            if cm['training_executed']:
                prob=z['probabilities']
                assert prob.shape==(int(pool.sum()),len(z['classes'])) and np.isfinite(prob).all()
                np.testing.assert_allclose(prob.sum(axis=1),1,rtol=1e-5,atol=1e-6)
                rounded=np.asarray([round(v,4) for v in prob.max(axis=1)],dtype=np.float32)
                pred=z['classes'][prob.argmax(axis=1)]
                for name,threshold in [('final090',.9),('final070',.7)]:
                    expected=np.where(rounded>=threshold,pred,'Unknown')
                    assert np.array_equal(z[name][pool],expected)
            # Every condition retains its terminal predictions; large model/probability
            # artifacts are shared only for byte-identical inputs, labels and parameters.
            for name in ['terminal.npz','training_manifest.json','training_history.json','model_state.pt']:
                if (cache/name).exists() and not (dest/name).exists():os.link(cache/name,dest/name)
            result=cells[['cell_id']].copy()
            for name in ['initial','final090','final070','lineage','confidence_rounded']:result[name]=z[name]
            tmp=dest/('predictions.csv.gz.part.'+execution_id())
            result.to_csv(tmp,index=False,compression='gzip');tmp.replace(dest/'predictions.csv.gz')
            tm=dict(status='completed',terminal_valid=True,arm=arm,family=family,
                cache_key=key,cache_directory=str(cache),identical_result_reused=reused,
                dl_status=cm['dl_status'],training_executed=cm['training_executed'],
                n_known=cm['n_known'],n_pool=cm['n_pool'],n_training_classes=cm['n_training_classes'],
                n_final_called090=cm['n_final_called090'],n_final_called070=cm['n_final_called070'],
                score_manifest_sha256=sha(source/'score_manifest.json'),DL_sha256=dlhash,DL_features=width,
                training_manifest_sha256=sha(cache/'training_manifest.json'),
                terminal_sha256=sha(dest/'terminal.npz'),predictions_sha256=sha(dest/'predictions.csv.gz'),
                known_labels_unchanged=True,threshold_reconstructed=True,
                reference_labels_used_for_fit=False,job=execution_id(),completed_at=utc())
            write_json(dest/'terminal_manifest.json',tm);complete(dest,'terminal_manifest.json','TERMINAL_COMPLETE')
            print('TERMINAL_COMPLETE',source.parent.name,source.name,aid,cm['dl_status'],'cached',reused,flush=True)
    if all(checked(terminal_root/aid,'terminal_manifest.json','TERMINAL_COMPLETE') for aid in m['arms']):
        (terminal_root/'ALL_TERMINAL_COMPLETE').write_text(sha(source/'score_manifest.json')+'\n')

def threads():
    require_slurm()
    import torch
    torch.set_num_threads(min(4,int(os.environ.get('SLURM_CPUS_PER_TASK',4))))
    torch.set_num_interop_threads(1)

if __name__=='__main__':
    threads()
    finish_route(sys.argv[1],None if len(sys.argv)<3 or sys.argv[2]=='all' else sys.argv[2],
                 sys.argv[3] if len(sys.argv)>3 else None)
