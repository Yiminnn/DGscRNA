"""Independent numerical checks for the archived model class and terminal bookkeeping."""
import ast
import json
import os
from pathlib import Path
import tempfile
from ptc_common import BASE, RECOVERY, require_slurm, sha, write_json
from refine import build_model, train_cache

require_slurm()
import numpy as np
import torch
import torch.nn as nn
torch.set_num_threads(1)
source=RECOVERY/'archive/tcr/ptc_val/scripts/DGscRNA-Share/R/source.py'
tree=ast.parse(source.read_text())
node=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='DeepModel')
namespace={'nn':nn}
exec(compile(ast.Module(body=[node],type_ignores=[]),str(source),'exec'),namespace)
torch.manual_seed(42)
reference=namespace['DeepModel'](5,3)
actual=build_model(5,3);actual.load_state_dict(reference.state_dict())
x=torch.tensor([[0,1,2,3,4],[1,0,-1,2,0],[0.1,0.2,0.3,0.4,0.5]],dtype=torch.float32)
y=torch.tensor([0,1,2])
assert torch.equal(reference(x),actual(x))
criterion=nn.CrossEntropyLoss()
for model in [reference,actual]:
    opt=torch.optim.Adamax(model.parameters(),lr=1e-3)
    loss=criterion(model(x),y);opt.zero_grad();loss.backward();opt.step()
for name,v in reference.state_dict().items():assert torch.equal(v,actual.state_dict()[name]),name
dest=BASE/'verification/implementation_v2';dest.mkdir(parents=True,exist_ok=True)
cases={
  'all_known':np.array(['T cell']*6+['B cell']*6),
  'none_known':np.array(['Undecided']*12),
  'single_known':np.array(['T cell']*10+['Undecided']*2),
  'two_known':np.array(['T cell']*5+['B cell']*5+['Undecided']*2),
  'three_known_empty_validation':np.array(['T cell']*3+['Undecided']*9)}
rng=np.random.default_rng(99);X=rng.normal(size=(12,5)).astype(np.float32)
states={}
for name,initial in cases.items():
  folder=dest/name;folder.mkdir(exist_ok=True)
  train_cache(X,initial,folder,{'verification_case':name})
  z=np.load(folder/'terminal.npz',allow_pickle=False)
  known=initial!='Undecided'
  assert np.array_equal(z['initial'],initial)
  assert np.array_equal(z['final090'][known],initial[known])
  assert np.array_equal(z['final070'][known],initial[known])
  info=json.loads((folder/'training_manifest.json').read_text());states[name]=info['dl_status']
  if info['training_executed']:
    assert len(json.loads((folder/'training_history.json').read_text()))==10
    assert len(set(z['train_indices'])&set(z['validation_indices']))==0
    assert set(z['train_indices'])|set(z['validation_indices'])==set(np.flatnonzero(known))
    model=build_model(5,len(z['classes']))
    model.load_state_dict(torch.load(folder/'model_state.pt',weights_only=True))
    with torch.no_grad():p=model(torch.from_numpy(X[z['pool_indices']])).numpy()
    np.testing.assert_allclose(p,z['probabilities'],rtol=0,atol=0)
    expected=z['classes'][p.argmax(1)]
    confidence=np.asarray([round(v,4) for v in p.max(1)],dtype=np.float32)
    expected90=np.where(confidence>=0.9,expected,'Unknown')
    expected70=np.where(confidence>=0.7,expected,'Unknown')
    assert np.array_equal(z['final090'][z['pool_indices']],expected90)
    assert np.array_equal(z['final070'][z['pool_indices']],expected70)
assert states=={'all_known':'no_op_all_initially_known',
 'none_known':'no_known_labels_archived_Undecided_terminal',
 'single_known':'trained_single_known_class','two_known':'trained',
 'three_known_empty_validation':'trained_single_known_class'}
write_json(dest/'verification.json',dict(status='passed',
  checks=['Archived class exact forward pass','Archived class exact first Adamax update',
          'All-known/no-known/single-known/two-known terminal cases','Known labels immutable',
          '10 training epochs','Disjoint exhaustive train/validation split',
          'Saved model reproduces pool probabilities','Both confidence thresholds replay'],
  states=states,source_sha256=sha(source),refinement_sha256=sha(Path(__file__).with_name('refine.py')),
  job=os.environ['SLURM_JOB_ID']))
(dest/'COMPLETE').write_text(sha(dest/'verification.json')+'\n')
print('Archived MLP and terminal-state checks passed')
