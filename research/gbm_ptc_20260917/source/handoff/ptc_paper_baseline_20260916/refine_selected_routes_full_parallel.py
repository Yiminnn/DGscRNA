"""Terminal DL for the notebook's original routes; reference labels are never fit."""
from pathlib import Path
import os
import sys
import json
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
sys.path.insert(0,str(ROOT/'handoff/ptc_recovery_20260916'))
import refine

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import numpy as np
    import pandas as pd
    import torch
    torch.set_num_threads(4)
    torch.set_num_interop_threads(1)
    task=int(os.environ['SLURM_ARRAY_TASK_ID'])
    route=['NMT_Thyroid_Seurat_none','TTU_Pubmed_UMAPHDBSCAN_mean'][task]
    dest=ROOT/'results/hvg_ptc_20260916_v1/ptc_paper_baseline/replay_selected_routes_full_parallel'/route
    assert (dest/'SCORE_COMPLETE').read_text().strip()==refine.sha(dest/'score_manifest.json')
    m=json.loads((dest/'score_manifest.json').read_text())
    calls=pd.read_csv(dest/'initial_calls.csv',keep_default_na=False,dtype=str)
    cells=pd.read_csv(dest/'cells.csv',keep_default_na=False,dtype=str)
    assert np.array_equal(cells.cell_id,calls.cell_id)
    binary=dest/'DL_archived_CCA2000.float32.bin'
    assert binary.stat().st_size==len(calls)*2000*4
    x=np.memmap(binary,mode='r',dtype='<f4',shape=(len(calls),2000))
    assert np.isfinite(x).all()
    out=dest/'terminal_DL';out.mkdir(exist_ok=True)
    assert not (out/'COMPLETE').exists()
    refine.PARAMS['input']='archived all-eight-sample CCA2000 integrated assay; original cell and feature order'
    refine.PARAMS['sensitivity_threshold']=0.70
    refine.PARAMS['sensitivity_note']='Stored probability-derived diagnostic only; no extra fit and not used to declare baseline parity'
    refine.train_cache(x,calls.initial.to_numpy(dtype=str),out,dict(route=route,
        score_manifest_sha256=refine.sha(dest/'score_manifest.json'),
        original_source_sha256=refine.sha(ROOT/'results/hvg_ptc_20260916_v1/ptc_recovery/archive/tcr/ptc_val/scripts/DGscRNA-Share/R/source.py'),
        initialization='Historical model seed/weights unavailable; explicitly seeded 42 reconstruction',
        training_scope=m['annotation_fit_scope'],reference_annotations_used_for_training=False))
    z=np.load(out/'terminal.npz',allow_pickle=False)
    result=cells.copy()
    result['initial_native']=z['initial']
    result['terminal_native_090']=z['final090']
    result['terminal_native_070_diagnostic']=z['final070']
    result['lineage']=z['lineage']
    result['confidence_rounded']=z['confidence_rounded']
    result['selected_for_paper_group']=result.group.eq(m['selected_group'])
    result.to_csv(dest/'terminal_predictions.csv.gz',index=False)
    (dest/'TERMINAL_COMPLETE').write_text(refine.sha(out/'training_manifest.json')+'\n')
    print(route,json.loads((out/'training_manifest.json').read_text())['dl_status'],flush=True)

if __name__=='__main__':run()
