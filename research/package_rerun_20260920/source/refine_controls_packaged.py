# Parameterized sensitivity copy of unchanged legacy_refine.py; original source SHA256: e32daabbf2e883666deace944020ebdab406f705f24200ae1cb3a8df9860d3c3
# Only model width, epoch count and explicit seeds are configurable. Native defaults require numerical parity.
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
from dgscrna.reference.backend.common import require_slurm, sha, utc, write_json

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
            self.fc1=nn.Linear(n_features,PARAMS["architecture"][0])
            self.fc2=nn.Linear(PARAMS["architecture"][0],PARAMS["architecture"][1])
            self.fc3=nn.Linear(PARAMS["architecture"][1],n_classes)
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
        random.seed(PARAMS["model_seed"]);np.random.seed(PARAMS["model_seed"]);torch.manual_seed(PARAMS["model_seed"])
        torch.use_deterministic_algorithms(True)
        # Includes Undecided in output vocabulary exactly as the archived implementation.
        lookup={name:i for i,name in enumerate(classes)}
        labels=np.array([lookup[v] for v in initial[known]],dtype=np.int64)
        X=torch.from_numpy(np.array(x[known],dtype=np.float32,copy=True))
        y=torch.from_numpy(labels)
        dataset=TensorDataset(X,y)
        ntrain=round(len(known)*0.90)
        train,val=random_split(dataset,[ntrain,len(known)-ntrain],generator=torch.Generator().manual_seed(PARAMS["split_seed"]))
        train_index=known[np.asarray(train.indices,dtype=np.int64)];val_index=known[np.asarray(val.indices,dtype=np.int64)]
        loader=DataLoader(train,batch_size=256,shuffle=False,num_workers=0)
        model=build_model(x.shape[1],len(classes))
        criterion=torch.nn.CrossEntropyLoss()
        optimizer=torch.optim.Adamax(model.parameters(),lr=1e-3)
        for epoch in range(PARAMS["epochs"]):
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
