"""Run one marker/cutoff condition through the validated original-style MLP."""
import json
import fcntl
import os
import sys
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'
# Freeze the driver before an array task loads or trains anything.
if __name__=='__main__' and not os.environ.get('DGSCRNA_TERMINAL_EXECUTION_SOURCE'):
    import shutil
    snapshot_dir=OUT/'execution_sources'/os.environ['SLURM_JOB_ID']
    snapshot_dir.mkdir(parents=True,exist_ok=True)
    snapshot_script=snapshot_dir/'terminal.py'
    if not snapshot_script.exists():shutil.copy2(__file__,snapshot_script)
    os.environ['DGSCRNA_TERMINAL_EXECUTION_SOURCE']=str(snapshot_script)
    os.execv(sys.executable,[sys.executable,str(snapshot_script),*sys.argv[1:]])
sys.path.insert(0,str(ROOT/'handoff/ptc_recovery_20260916'))
import refine

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import numpy as np
    import pandas as pd
    import torch
    torch.set_num_threads(min(4,int(os.environ.get('SLURM_CPUS_PER_TASK',4))))
    torch.set_num_interop_threads(1)
    task=int(os.environ['SLURM_ARRAY_TASK_ID'])
    if sys.argv[1]=='PTC':
        cols=['seurat_clusters','seurat.UMAP_clusters','hdbscan_clusters','hdbscan.UMAP_clusters']
        source=OUT/'PTC_archived_CCA2000'/cols[task//51]
        armid=f'L{(task%51)//3:02d}_{["none","mean","p050"][task%3]}'
    else:
        spec=json.loads(Path(sys.argv[1]).read_text())[task]
        source=Path(spec['source']);armid=spec['arm_id']
    dest=source/'terminal'/armid;dest.mkdir(parents=True,exist_ok=True)
    with (dest/'terminal.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX)
        finish_arm(source,armid,dest)

def finish_arm(source,armid,dest):
    import numpy as np
    import pandas as pd
    import torch
    assert (source/'SCORE_COMPLETE').read_text().strip()==refine.sha(source/'score_manifest.json')
    m=json.loads((source/'score_manifest.json').read_text())
    arm=m['arms'][armid]
    trained=(dest/'COMPLETE').exists()
    if trained:
        previous=json.loads((dest/'training_manifest.json').read_text())
        assert (dest/'COMPLETE').read_text().strip()==refine.sha(dest/'training_manifest.json')
        assert previous['provenance']['score_manifest_sha256']==refine.sha(source/'score_manifest.json')
        if (dest/'TERMINAL_COMPLETE').exists():
            assert (dest/'TERMINAL_COMPLETE').read_text().strip()==refine.sha(dest/'training_manifest.json')
            print('Already complete',dest,flush=True);return
    assert m['initial_sha256']==refine.sha(source/'initial_calls.csv.gz')
    initial=pd.read_csv(source/'initial_calls.csv.gz',usecols=['cell_id',arm['seed_column']],keep_default_na=False,dtype=str)
    cells=pd.read_csv(source/'cells.csv',keep_default_na=False,dtype=str)
    assert np.array_equal(cells.cell_id,initial.cell_id)
    binary=Path(m['DL_binary']);nfeatures=int(m['DL_features'])
    assert binary.stat().st_size==len(cells)*nfeatures*4
    assert refine.sha(binary)==m['DL_binary_sha256']
    x=np.memmap(binary,mode='r',dtype='<f4',shape=(len(cells),nfeatures))
    assert np.isfinite(x).all()
    refine.PARAMS['input']=m.get('DL_input_description',m.get('geometry','R reference expression'))
    labels=initial[arm['seed_column']].to_numpy(dtype=str)
    if not trained:
        refine.train_cache(x,labels,dest,dict(arm=arm,source=str(source),
            score_manifest_sha256=refine.sha(source/'score_manifest.json'),
            training_source_sha256=refine.sha(Path(refine.__file__)),
            driver_source_sha256=refine.sha(Path(__file__)),execution_source=str(Path(__file__)),
            task_assignment_source=sys.argv[1],
            task_assignment_sha256=refine.sha(Path(sys.argv[1])) if sys.argv[1]!='PTC' else None,
            torch_num_threads=torch.get_num_threads(),
            historical_model_weights_recovered=False,evaluation_labels_used_for_fit=False))
    z=np.load(dest/'terminal.npz',allow_pickle=False)
    result=cells[['cell_id']].copy()
    for name in ['initial','final090','final070','lineage','confidence_rounded']:result[name]=z[name]
    temporary=dest/('predictions.csv.gz.part.'+os.environ['SLURM_JOB_ID'])
    result.to_csv(temporary,index=False,compression='gzip');temporary.replace(dest/'predictions.csv.gz')
    temporary=dest/('TERMINAL_COMPLETE.part.'+os.environ['SLURM_JOB_ID'])
    temporary.write_text(refine.sha(dest/'training_manifest.json')+'\n');temporary.replace(dest/'TERMINAL_COMPLETE')
    print(source.name,armid,json.loads((dest/'training_manifest.json').read_text())['dl_status'],flush=True)

if __name__=='__main__':run()
