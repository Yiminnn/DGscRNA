"""Official unchanged scDeepSort predictor; isolated Python3.8 environment."""
from pathlib import Path
import os,sys,json,time,hashlib
assert os.environ.get('SLURM_JOB_ID')
os.environ['DGLBACKEND']='pytorch'
import numpy as np
import pandas as pd
import torch,dgl,deepsort
from deepsort import DeepSortPredictor

dest=Path(sys.argv[1]);m=json.loads((dest/'input_manifest.json').read_text())
def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()
inp=dest/'lognorm_model_input.csv'
assert sha(inp)==m['input_csv_sha256']
started=time.monotonic()
model=DeepSortPredictor(species='human',tissue='Brain',unsure_rate=2.)
raw=model.predict(str(inp))
cells=pd.read_csv(dest/'cells.csv',dtype=str).cell_id.tolist()
assert raw['index'].astype(str).tolist()==cells
raw=raw.rename(columns={'index':'cell_id'})
raw['default']=raw.cell_type.replace({'unsure':'Unknown'})
raw.to_csv(dest/'predictions.csv.gz',index=False,compression='gzip')
proof=dict(status='predicted',sample=m['sample'],n_cells=len(cells),n_shared_genes=m['n_shared_genes'],
    input='Seurat-compatible LogNormalize before model gene intersection; float64 CSV; no integer cast',
    native_unsure_rate=2.,threshold='2 / 21 classes',training_executed=False,
    pretrained_model='published human Brain',versions=dict(python=sys.version,numpy=np.__version__,pandas=pd.__version__,
        torch=torch.__version__,dgl=dgl.__version__),elapsed_inference_seconds=time.monotonic()-started,
    source_sha256=sha(__file__),official_predictor_source=str(Path(deepsort.__file__).with_name('predict.py')),
    official_predictor_sha256=sha(Path(deepsort.__file__).with_name('predict.py')),
    predictions_sha256=sha(dest/'predictions.csv.gz'),input_manifest_sha256=sha(dest/'input_manifest.json'),
    job=os.environ['SLURM_JOB_ID'])
(dest/'prediction_manifest.json').write_text(json.dumps(proof,indent=2)+'\n')
print('CORRECTED_NATIVE_PREDICTION',json.dumps(proof),flush=True)
