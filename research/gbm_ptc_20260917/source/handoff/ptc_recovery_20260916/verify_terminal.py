"""Independent per-model/per-cell verification. Does not import fitting/refinement code."""
from pathlib import Path
import hashlib
import json
import os
import uuid
from functools import lru_cache
from ptc_common import BASE,sha,utc

@lru_cache(maxsize=4)
def binary_info(group,space):
    import pandas as pd
    prep=BASE/'prepared'/group
    path=prep/f'DL_{space}.float32.bin'
    cells=pd.read_csv(prep/'cells.csv',index_col=0).index
    assert cells.is_unique and path.stat().st_size==len(cells)*2000*4
    return path,cells,sha(path),sha(prep/'integration_features.txt')

def independent_model(state,n_features,n_classes):
    import torch
    import torch.nn as nn
    class Model(nn.Module):
      def __init__(self):
        super().__init__()
        self.fc1=nn.Linear(n_features,256)
        self.fc2=nn.Linear(256,128)
        self.fc3=nn.Linear(128,n_classes)
      def forward(self,x):
        x=torch.nn.functional.leaky_relu(self.fc1(x),negative_slope=0.01)
        x=torch.nn.functional.leaky_relu(self.fc2(x),negative_slope=0.01)
        return torch.softmax(self.fc3(x),dim=1)
    model=Model();model.load_state_dict(state);model.eval();return model

def verify_terminal(dest,initial,cells):
    import numpy as np
    import torch
    assert os.environ.get('SLURM_JOB_ID')
    tm=json.loads((dest/'terminal_manifest.json').read_text())
    assert (dest/'TERMINAL_COMPLETE').read_text().strip()==sha(dest/'terminal_manifest.json')
    assert tm['terminal_valid'] and tm['status']=='completed'
    assert sha(dest/'terminal.npz')==tm['terminal_sha256']
    cache=Path(tm['cache_directory'])
    assert cache.parent==BASE/'refinement_cache'
    assert sha(cache/'training_manifest.json')==tm['cache_manifest_sha256']
    cm=json.loads((cache/'training_manifest.json').read_text())
    assert (cache/'COMPLETE').read_text().strip()==tm['cache_manifest_sha256']
    assert cm['provenance']['cache_key']==tm['cache_key']
    assert cm['provenance']['input']==tm['input']
    assert cm['provenance']['source']==tm['source_record']
    assert cm['params']['epochs']==10 and cm['params']['batch_size']==256
    assert cm['params']['optimizer']=='Adamax' and cm['params']['lr']==0.001
    assert cm['params']['model_seed']==42 and cm['params']['split_seed']==42
    assert cm['params']['output']=='Softmax(dim=1)'
    assert tm['input']['cell_order_sha256']==hashlib.sha256('\n'.join(cells).encode()).hexdigest()
    z=np.load(dest/'terminal.npz',allow_pickle=False)
    assert np.array_equal(z['initial'],initial)
    assert len(initial)==cm['n_cells']==len(cells)
    known=np.flatnonzero(initial!='Undecided');pool=np.flatnonzero(initial=='Undecided')
    assert np.array_equal(z['pool_indices'],pool)
    assert np.array_equal(z['classes'],np.array(sorted(set(initial)),dtype=str))
    assert np.array_equal(z['final090'][known],initial[known])
    assert np.array_equal(z['final070'][known],initial[known])
    vdir=BASE/'verification/cache';vdir.mkdir(parents=True,exist_ok=True)
    verification=vdir/f'{tm["cache_key"]}.json'
    if verification.exists():
      vm=json.loads(verification.read_text())
      assert vm['cache_manifest_sha256']==tm['cache_manifest_sha256'] and vm['status']=='passed'
    else:
      for name,h in cm['outputs'].items():assert sha(cache/name)==h,(cache,name)
      path,gcells,binary_sha,feature_sha=binary_info(tm['input']['group'],tm['input']['space'])
      assert binary_sha==tm['input']['binary_sha256'] and feature_sha==tm['input']['feature_sha256']
      positions=gcells.get_indexer(cells);assert (positions>=0).all()
      if cm['training_executed']:
        assert len(known)>=2 and len(pool)>0
        history=json.loads((cache/'training_history.json').read_text())
        assert len(history)==10 and [h['epoch'] for h in history]==list(range(1,11))
        assert all(np.isfinite(h['mean_loss']) for h in history)
        train=z['train_indices'];val=z['validation_indices']
        assert len(train)==round(0.9*len(known))
        assert set(train).isdisjoint(val) and set(train)|set(val)==set(known)
        split=torch.randperm(len(known),generator=torch.Generator().manual_seed(42)).numpy()
        assert np.array_equal(train,known[split[:len(train)]])
        assert np.array_equal(val,known[split[len(train):]])
        state=torch.load(cache/'model_state.pt',map_location='cpu',weights_only=True)
        model=independent_model(state,2000,len(z['classes']))
        mm=np.memmap(path,dtype='<f4',mode='r',shape=(len(gcells),2000))
        predicted=[]
        with torch.no_grad():
          for lo in range(0,len(pool),256):
            x=np.asarray(mm[positions[pool[lo:lo+256]]],dtype=np.float32)
            predicted.append(model(torch.from_numpy(x)).numpy())
        p=np.concatenate(predicted)
        assert np.isfinite(p).all()
        np.testing.assert_allclose(p,z['probabilities'],rtol=1e-5,atol=1e-6)
        np.testing.assert_allclose(p.sum(1),1,rtol=0,atol=1e-5)
        # Exact terminal replay uses the recorded probabilities, not a rounded display.
        recorded=z['probabilities'];confidence=np.array([round(v,4) for v in recorded.max(1)],dtype=np.float32)
        names=z['classes'][recorded.argmax(1)]
        assert np.array_equal(z['predicted'][pool],names)
        assert np.array_equal(z['confidence_rounded'][pool],confidence)
        for cutoff,key in [(0.9,'final090'),(0.7,'final070')]:
          expected=np.where(confidence>=cutoff,names,'Unknown')
          assert np.array_equal(z[key][pool],expected)
        assert np.array_equal(z['lineage'][pool],np.where(confidence>=0.9,'DL_accepted','DL_low_confidence'))
        max_error=float(np.max(np.abs(p-recorded)))
      else:
        if len(pool)==0:assert cm['dl_status']=='no_op_all_initially_known'
        elif len(known)==0:
          assert cm['dl_status']=='no_known_labels_archived_Undecided_terminal'
          assert (z['final090']=='Undecided').all() and (z['final070']=='Undecided').all()
        else:
          assert len(known)==1 and cm['dl_status']=='structural_insufficient_known_split'
          assert (z['final090'][pool]=='Unknown').all()
        max_error=None
      vm=dict(status='passed',cache_key=tm['cache_key'],cache_manifest_sha256=tm['cache_manifest_sha256'],
         model_forward_checked=cm['training_executed'],max_probability_difference=max_error,
         n_cells=len(cells),n_known=len(known),n_pool=len(pool),dl_status=cm['dl_status'],
         verifier_source_sha256=sha(Path(__file__)),job=os.environ['SLURM_JOB_ID'],completed_at=utc())
      # Immutable write-once cache: atomic hard-link publication works across nodes.
      tmp=vdir/f'.{tm["cache_key"]}.{uuid.uuid4().hex}.tmp';tmp.write_text(json.dumps(vm,indent=2)+'\n')
      try:os.link(tmp,verification)
      except FileExistsError:
        prior=json.loads(verification.read_text())
        assert prior['status']=='passed' and prior['cache_manifest_sha256']==tm['cache_manifest_sha256']
      finally:tmp.unlink()
    return tm,z
