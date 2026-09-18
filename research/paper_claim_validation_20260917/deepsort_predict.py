"""Published scDeepSort1.0 CPU predictor; standalone Python3.8 compatibility."""
import hashlib
import json
import os
from pathlib import Path
import sys
import time
assert os.environ.get('SLURM_JOB_ID')
os.environ['DGLBACKEND']='pytorch'
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch
import dgl
from deepsort import DeepSortPredictor

ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
BASE=ROOT/'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917'
MODEL=Path('/fs/scratch/PCON0080/yimin/tools/deepsort-1.0/deepsort-pretrained')
def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def run(sample):
    dest=BASE/'comparators/scDeepSort'/sample;dest.mkdir(parents=True,exist_ok=True)
    if (dest/'COMPLETE').exists():return
    assert json.loads((BASE/'scDeepSort_install.json').read_text())['status']=='installed_and_checkpoint_loaded'
    src=BASE/'inputs'/sample;m=json.loads((src/'input_manifest.json').read_text())
    assert (src/'INPUT_COMPLETE').read_text().strip()==sha(src/'input_manifest.json')
    X=sp.csr_matrix((np.fromfile(src/'x.bin',dtype='<f8'),np.fromfile(src/'i.bin',dtype='<i4'),
        np.fromfile(src/'p.bin',dtype='<i4')),shape=(m['n_cells'],m['n_genes']))
    genes=pd.read_csv(src/'genes.csv',dtype=str).gene
    cells=pd.read_csv(src/'cells_fit.csv',dtype=str).cell_id
    allowed=set((MODEL/'human/statistics/Brain_genes.txt').read_text().splitlines())
    keep=genes.isin(allowed).to_numpy();shared=genes[keep]
    assert not shared.duplicated().any() and keep.sum()>1000
    X=X[:,keep].astype(np.int32)
    inp=dest/'temporary_published_counts_input.csv'
    pd.DataFrame(X.T.toarray(),index=shared,columns=cells).to_csv(inp)
    input_hash=sha(inp);start=time.monotonic()
    model=DeepSortPredictor(species='human',tissue='Brain')
    raw=model.predict(str(inp))
    assert list(raw['index'].astype(str))==list(cells)
    raw=raw.rename(columns={'index':'cell_id'})
    raw['default']=raw.cell_type.replace({'unsure':'Unknown'})
    raw.to_csv(dest/'predictions.csv.gz',index=False,compression='gzip')
    record=dict(status='completed',method='scDeepSort',sample=sample,n_cells=len(cells),n_overlap_genes=int(keep.sum()),
        input='Original eligible integer counts. Only genes supported by the published checkpoint are exported; the tool itself performs this same filter.',
        input_csv_sha256=input_hash,input_manifest_sha256=sha(src/'input_manifest.json'),
        model_sha256=sha(MODEL/'human/models/human-Brain.pt'),
        pretrained_label_sha256=sha(MODEL/'human/statistics/Brain_cell_type.txt'),
        n_pretrained_classes=21,unsure_rate=2.,threshold=2./21,
        supervision='Published human cell atlas Brain checkpoint; no local author labels or marker tuning',
        vocabulary_limitation='Checkpoint lacks a malignant class. Generic/fetal neurons cannot be relabelled as excitatory/inhibitory using evaluation labels.',
        versions=dict(python=sys.version,torch=torch.__version__,dgl=dgl.__version__,numpy=np.__version__),
        elapsed_seconds=time.monotonic()-start,job=os.environ['SLURM_JOB_ID'],source_sha256=sha(__file__),
        predictions_sha256=sha(dest/'predictions.csv.gz'))
    (dest/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    (dest/'COMPLETE').write_text(sha(dest/'manifest.json')+'\n')
    # This is a deterministic temporary export, not an original/raw data file.
    inp.unlink()

if __name__=='__main__':run(sys.argv[1])
