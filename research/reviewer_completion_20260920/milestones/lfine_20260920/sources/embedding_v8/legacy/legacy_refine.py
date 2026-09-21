"""Faithful PTC MLP refinement with explicit seeds, histories, weights and terminal states."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import random
import shutil
import sys
import time
from ptc_common import BASE, RECOVERY, require_slurm, sha, utc, write_json, task_list, geometry_dir

PARAMS=dict(input='fixed group-selected RNA2000 or separately named CCA2000',
    architecture=[256,128],activation='LeakyReLU(0.01)',output='Softmax(dim=1)',
    loss='CrossEntropyLoss on softmax probabilities (historical implementation)',
    optimizer='Adamax',lr=0.001,epochs=10,batch_size=256,shuffle=False,
    split_seed=42,model_seed=42,known_train_fraction=0.90,confidence_round_decimals=4,
    primary_threshold=0.90,sensitivity_threshold=0.70,num_workers=0,
    num_workers_note='Deterministic in-memory tensor loading; preserves historical order and batch contents')

def build_model(n_features,n_classes):
    import torch.nn as nn
    class DeepModel(nn.Module):
        def __init__(self):
            super().__init__()
            self.activation=nn.LeakyReLU()
            self.output_activation=nn.Softmax(dim=1)
            self.fc1=nn.Linear(n_features,256)
            self.fc2=nn.Linear(256,128)
            self.fc3=nn.Linear(128,n_classes)
        def forward(self,x):
            return self.output_activation(self.fc3(self.activation(self.fc2(self.activation(self.fc1(x))))))
    return DeepModel()

def train_cache(x,initial,dest,provenance):
    import numpy as np
    import torch
    from torch.utils.data import TensorDataset,DataLoader,random_split
    start=time.perf_counter()
    known=np.flatnonzero(initial!='Undecided')
    pool=np.flatnonzero(initial=='Undecided')
    classes=np.asarray(sorted(set(initial)),dtype=str)
    final90=initial.copy();final70=initial.copy()
    predicted=np.full(len(initial),'',dtype=initial.dtype)
    confidence=np.full(len(initial),np.nan,dtype=np.float32)
    lineage=np.asarray(np.where(initial=='Undecided','unresolved','retained_initial'),dtype='U32')
    history=[]
    info=dict(params=PARAMS,provenance=provenance,started_at=utc(),job=os.environ['SLURM_JOB_ID'],
      n_cells=len(initial),n_known=len(known),n_pool=len(pool),n_training_classes=len(set(initial[known])),
      classes=classes.tolist(),input_width=x.shape[1],training_executed=False,prediction_executed=False,
      terminal_valid=True)
    probabilities=np.empty((0,len(classes)),dtype=np.float32)
    train_index=np.array([],dtype=int);val_index=np.array([],dtype=int)
    if len(pool)==0:
        info['dl_status']='no_op_all_initially_known'
    elif len(known)==0:
        info['dl_status']='no_known_labels_archived_Undecided_terminal'
    elif len(known)<2:
        info['dl_status']='structural_insufficient_known_split'
        final90[pool]='Unknown';final70[pool]='Unknown'
    else:
        random.seed(42);np.random.seed(42);torch.manual_seed(42)
        torch.use_deterministic_algorithms(True)
        # Includes Undecided in output vocabulary exactly as the archived implementation.
        lookup={name:i for i,name in enumerate(classes)}
        labels=np.array([lookup[v] for v in initial[known]],dtype=np.int64)
        X=torch.from_numpy(np.array(x[known],dtype=np.float32,copy=True))
        y=torch.from_numpy(labels)
        dataset=TensorDataset(X,y)
        ntrain=round(len(known)*0.90)
        train,val=random_split(dataset,[ntrain,len(known)-ntrain],generator=torch.Generator().manual_seed(42))
        train_index=known[np.asarray(train.indices,dtype=np.int64)];val_index=known[np.asarray(val.indices,dtype=np.int64)]
        loader=DataLoader(train,batch_size=256,shuffle=False,num_workers=0)
        model=build_model(x.shape[1],len(classes))
        criterion=torch.nn.CrossEntropyLoss()
        optimizer=torch.optim.Adamax(model.parameters(),lr=1e-3)
        for epoch in range(10):
            total=0;correct=0;weighted_loss=0.0;legacy_loss=0.0
            for samples,targets in loader:
                probabilities_batch=model(samples)
                loss=criterion(probabilities_batch,targets)
                assert torch.isfinite(loss)
                optimizer.zero_grad();loss.backward();optimizer.step()
                n=len(targets);total+=n
                correct+=int((probabilities_batch.argmax(1)==targets).sum())
                weighted_loss+=float(loss.detach())*n
                legacy_loss+=float(loss.detach())*256
            history.append(dict(epoch=epoch+1,n=total,mean_loss=weighted_loss/total,
                                legacy_batch_size_weighted_sum=legacy_loss,accuracy=correct/total))
        with torch.no_grad():
            val_outputs=[];val_targets=[]
            for samples,targets in DataLoader(val,batch_size=256,shuffle=False,num_workers=0):
                val_outputs.append(model(samples));val_targets.append(targets)
            if val_outputs:
                vo=torch.cat(val_outputs);vy=torch.cat(val_targets)
                info['known_seed_validation_accuracy']=float((vo.argmax(1)==vy).float().mean())
            else:
                info['known_seed_validation_accuracy']=None
            chunks=[]
            for lo in range(0,len(pool),256):
                batch=torch.from_numpy(np.array(x[pool[lo:lo+256]],dtype=np.float32,copy=True))
                chunks.append(model(batch).numpy())
            probabilities=np.concatenate(chunks)
        assert probabilities.shape==(len(pool),len(classes)) and np.isfinite(probabilities).all()
        raw=probabilities.max(axis=1)
        rounded=np.asarray([round(v,4) for v in raw],dtype=np.float32)
        pool_calls=classes[probabilities.argmax(axis=1)]
        predicted[pool]=pool_calls;confidence[pool]=rounded
        for threshold,result in [(0.9,final90),(0.7,final70)]:
            accepted=rounded>=threshold
            result[pool[accepted]]=pool_calls[accepted]
            result[pool[~accepted]]='Unknown'
        lineage[pool]=np.where(rounded>=0.9,'DL_accepted','DL_low_confidence')
        info.update(dl_status='trained_single_known_class' if len(set(labels))==1 else 'trained',
                    training_executed=True,prediction_executed=True,n_train=len(train_index),n_validation=len(val_index))
        torch.save(model.state_dict(),dest/'model_state.pt')
        write_json(dest/'training_history.json',history)
    assert np.array_equal(final90[known],initial[known]) and np.array_equal(final70[known],initial[known])
    np.savez_compressed(dest/'terminal.npz',initial=initial,final090=final90,final070=final70,
       classes=classes,pool_indices=pool,probabilities=probabilities,predicted=predicted,
       confidence_rounded=confidence,lineage=lineage,train_indices=train_index,validation_indices=val_index)
    info.update(n_final_called090=int((~np.isin(final90,['Unknown','Undecided'])).sum()),
                n_final_called070=int((~np.isin(final70,['Unknown','Undecided'])).sum()),
                elapsed_seconds=time.perf_counter()-start,completed_at=utc())
    info['outputs']={p.name:sha(p) for p in dest.iterdir() if p.is_file()}
    write_json(dest/'training_manifest.json',info)
    (dest/'COMPLETE').write_text(sha(dest/'training_manifest.json')+'\n')

def refine_score(score,context,source_hash):
    import numpy as np
    import pandas as pd
    assert (score/'SCORE_COMPLETE').read_text().strip()==sha(score/'score_manifest.json')
    m=json.loads((score/'score_manifest.json').read_text())
    assert m['initial_sha256']==sha(score/'initial_calls.csv.gz')
    seeds=pd.read_csv(score/'initial_calls.csv.gz',keep_default_na=False,dtype=str).set_index('cell_id')
    group=context['group']
    prep=BASE/'prepared'/group
    gcells=pd.read_csv(prep/'cells.csv',index_col=0).index
    positions=gcells.get_indexer(seeds.index)
    assert (positions>=0).all() and not seeds.index.has_duplicates
    space='CCA2000' if m['assay']=='integrated' else 'RNA2000'
    binary=prep/f'DL_{space}.float32.bin'
    assert binary.stat().st_size==len(gcells)*2000*4
    mm=np.memmap(binary,mode='r',dtype='<f4',shape=(len(gcells),2000))
    x=np.array(mm[positions],dtype=np.float32,copy=True)
    assert np.isfinite(x).all()
    input_record=dict(group=group,sample=context.get('sample'),space=space,n_cells=len(x),
        binary_sha256=sha(binary),cell_order_sha256=hashlib.sha256('\n'.join(seeds.index).encode()).hexdigest(),
        feature_sha256=sha(prep/'integration_features.txt'))
    source_record=dict(refinement_source_sha256=source_hash,
       original_source_sha256=sha(RECOVERY/'archive/tcr/ptc_val/scripts/DGscRNA-Share/R/source.py'))
    input_signature=json.dumps(dict(input=input_record,source=source_record,params=PARAMS),sort_keys=True).encode()
    cache_root=BASE/'refinement_cache';cache_root.mkdir(exist_ok=True)
    for aid,arm in m['arms'].items():
        dest=score/'refinement'/aid;dest.mkdir(parents=True,exist_ok=True)
        if (dest/'TERMINAL_COMPLETE').exists():
            old=json.loads((dest/'terminal_manifest.json').read_text())
            assert (dest/'TERMINAL_COMPLETE').read_text().strip()==sha(dest/'terminal_manifest.json')
            # V2 only fixes empty validation-index dtype and cross-node publication.
            # Existing V1 terminal results keep their original provenance and hashes.
            assert old['source_record']['original_source_sha256']==source_record['original_source_sha256']
            assert old['source_record']['refinement_source_sha256'] in {
                source_hash,'cfb83584e2d86c0c3ebcd4101108fd74d1ad520dcfe3dce477cb55ebaf497179'}
            continue
        initial=seeds[aid].to_numpy(dtype=str)
        key=hashlib.sha256(input_signature+b'\0'+json.dumps(initial.tolist(),ensure_ascii=False).encode()).hexdigest()
        cache=cache_root/key
        with open(cache_root/f'{key}.lock','a') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            reused=(cache/'COMPLETE').exists()
            if not reused:
                temp=cache_root/f'{key}.tmp.{os.environ["SLURM_JOB_ID"]}.{os.getpid()}'
                temp.mkdir(exist_ok=False)
                provenance=dict(cache_key=key,input=input_record,source=source_record,first_condition=str(dest),
                    initial_labels_sha256=hashlib.sha256(json.dumps(initial.tolist(),ensure_ascii=False).encode()).hexdigest())
                train_cache(x,initial,temp,provenance)
                try:
                    temp.rename(cache)
                except FileExistsError:
                    # Some shared filesystems implement flock locally to a node.
                    # Atomic publication chooses one complete result; verify the
                    # concurrent duplicate before referencing that published result.
                    assert (cache/'COMPLETE').exists()
                    with np.load(temp/'terminal.npz',allow_pickle=False) as a, np.load(cache/'terminal.npz',allow_pickle=False) as b:
                        for name in ['initial','classes','pool_indices','final090','final070','train_indices','validation_indices']:
                            assert np.array_equal(a[name],b[name]),f'Concurrent deterministic result differs: {name}'
                        np.testing.assert_allclose(a['probabilities'],b['probabilities'],rtol=1e-6,atol=1e-7)
                    write_json(temp/'DUPLICATE_TERMINAL_EQUIVALENCE.json',dict(
                        published_cache=str(cache),all_terminal_labels_and_splits_exact=True,
                        probability_rtol=1e-6,probability_atol=1e-7,
                        note='Cross-node duplicate computation retained; canonical complete publication reused'))
                    duplicate_root=cache_root/'duplicate_publications';duplicate_root.mkdir(exist_ok=True)
                    temp.rename(duplicate_root/temp.name)
                    reused=True
            cm=json.loads((cache/'training_manifest.json').read_text())
            assert (cache/'COMPLETE').read_text().strip()==sha(cache/'training_manifest.json')
            assert cm['provenance']['cache_key']==key
            if not (dest/'terminal.npz').exists():
                os.link(cache/'terminal.npz',dest/'terminal.npz')
        terminal=dict(status='completed',terminal_valid=True,context=context,arm=arm,
            dl_status=cm['dl_status'],model_training_executed_for_cache=cm['training_executed'],
            identical_result_reused=reused,cache_key=key,cache_directory=str(cache),
            cache_manifest_sha256=sha(cache/'training_manifest.json'),
            source_record=source_record,input=input_record,score_manifest_sha256=sha(score/'score_manifest.json'),
            terminal_sha256=sha(dest/'terminal.npz'),n_known=cm['n_known'],n_pool=cm['n_pool'],
            n_training_classes=cm['n_training_classes'],n_final_called090=cm['n_final_called090'],
            n_final_called070=cm['n_final_called070'],completed_at=utc(),job=os.environ['SLURM_JOB_ID'],
            reuse_rule='Only byte-identical input identity/cell order, initial labels, model seed and training source share a cache')
        write_json(dest/'terminal_manifest.json',terminal)
        (dest/'TERMINAL_COMPLETE').write_text(sha(dest/'terminal_manifest.json')+'\n')
        print(group,context.get('sample','pooled'),aid,cm['dl_status'],'reused' if reused else 'new',flush=True)
    (score/'REFINEMENT_COMPLETE').write_text(sha(score/'score_manifest.json')+'\n')

def run():
    require_slurm()
    import torch
    torch.set_num_threads(min(4,int(os.environ.get('SLURM_CPUS_PER_TASK',4))))
    torch.set_num_interop_threads(1)
    mode=sys.argv[1]
    task_id=int(sys.argv[2] if len(sys.argv)>2 else os.environ['SLURM_ARRAY_TASK_ID'])
    source_hash=sha(Path(__file__))
    if mode in ['single','pilot']:
        t=task_list()[task_id];gd=geometry_dir(t)
        for method in t['clusterers']:
            if mode=='pilot' and method!='HDBSCAN':continue
            refine_score(gd/method/'score_RNA',{**t,'clusterer':method},source_hash)
        if mode=='single':(gd/'ANNOTATION_COMPLETE').write_text(utc()+'\n')
    elif mode=='pooled':
        group=['MTN','TUT'][task_id//3];correction=['NONE','CCA','HARMONY'][task_id%3]
        gd=BASE/'pooled'/group/correction
        for space in ['PCA30','UMAP2']:
            for method in ['SNN','HDBSCAN_R']:
                for assay in (['RNA','integrated'] if correction=='CCA' else ['RNA']):
                    context=dict(group=group,correction=correction,space=space,clusterer=method,seed=42,
                        family='legacy_integrated_scoring_and_DL' if assay=='integrated' else 'controlled_batch')
                    refine_score(gd/f'{space}_{method}'/f'score_{assay}',context,source_hash)
        (gd/'ANNOTATION_COMPLETE').write_text(utc()+'\n')
    else:raise ValueError(mode)

if __name__=='__main__':run()
