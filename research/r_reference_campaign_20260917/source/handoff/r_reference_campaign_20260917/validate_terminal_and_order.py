import hashlib,json,os
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
BASE=ROOT/'results/hvg_ptc_20260916_v1'
OUT=BASE/'r_reference_campaign_20260917'

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import numpy as np
    import pandas as pd
    checks=[]
    spec=json.loads((OUT/'ptc_selected_initial_parity.json').read_text())['checks']
    for name,cl in [('NMT_Thyroid_Seurat_none','seurat_clusters'),('TTU_Pubmed_UMAPHDBSCAN_mean','hdbscan.UMAP_clusters')]:
        old=BASE/'ptc_paper_baseline/replay_selected_routes_full_parallel'/name/'terminal_DL/terminal.npz'
        new=OUT/'PTC_archived_CCA2000'/cl/'terminal'/spec[name]['arm']/'terminal.npz'
        a=np.load(old,allow_pickle=False);b=np.load(new,allow_pickle=False)
        exact={k:bool(np.array_equal(a[k],b[k])) for k in ['initial','final090','final070','classes','train_indices','validation_indices']}
        assert all(exact.values()),(name,exact)
        delta=float(np.max(np.abs(a['probabilities']-b['probabilities']))) if a['probabilities'].size else 0.0
        assert delta<1e-6
        checks.append(dict(route=name,exact=exact,max_probability_delta=delta))
    orders=[]
    for mpath in sorted((OUT/'PTC_ablation').glob('*/prepare_manifest.json')):
        m=json.loads(mpath.read_text());prep=mpath.parent
        if m['group']=='ALL8':
            expected=pd.read_csv(BASE/'ptc_paper_baseline/replay_selected_routes_full_parallel/NMT_Thyroid_Seurat_none/cells.csv').cell_id.to_numpy()
        else:
            oldgroup={'NMT':'MTN','TTU':'TUT'}[m['group']]
            expected=pd.read_csv(BASE/'ptc_experiments/prepared'/oldgroup/'cells.csv',index_col=0).index.to_numpy()
        got=pd.read_csv(prep/'cells.csv').cell_id.to_numpy()
        orders.append(dict(unit=m['unit'],same_cell_set=set(expected)==set(got),same_cell_order=bool(np.array_equal(expected,got))))
    report=dict(job=os.environ['SLURM_JOB_ID'],selected_terminal_parity=checks,PTC_cell_order_audit=orders)
    (OUT/'selected_terminal_and_cell_order_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':run()
