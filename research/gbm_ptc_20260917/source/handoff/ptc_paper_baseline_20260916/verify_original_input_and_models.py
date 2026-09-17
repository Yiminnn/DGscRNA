"""Independently check archived loom input and saved terminal model predictions."""
from pathlib import Path
import os,json,sys,hashlib
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
BASE=ROOT/'results/hvg_ptc_20260916_v1'
OUT=BASE/'ptc_paper_baseline'

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import h5py,numpy as np,pandas as pd,torch
    torch.set_num_threads(2)
    sys.path.insert(0,str(ROOT/'handoff/ptc_recovery_20260916'))
    from verify_terminal import independent_model
    parent=OUT/'replay_selected_routes_marker_union'
    first=parent/'NMT_Thyroid_Seurat_none'
    cells=pd.read_csv(first/'cells.csv',dtype=str)
    features=(first/'integration_features.txt').read_text().splitlines()
    binary=first/'DL_archived_CCA2000.float32.bin'
    x=np.memmap(binary,mode='r',dtype='<f4',shape=(len(cells),len(features)))
    source=BASE/'ptc_recovery/archive/tcr/rawdata/integrated_data.loom'
    maxerr=0.;mismatches=0
    with h5py.File(source,'r') as h:
        assert np.array_equal(h['col_attrs/CellID'].asstr()[:],cells.original_cell_id)
        assert np.array_equal(h['row_attrs/Gene'].asstr()[:],features)
        assert h['matrix'].shape==(2000,92404)
        for lo in range(0,len(cells),512):
            stop=min(lo+512,len(cells))
            old=h['matrix'][:,lo:stop].T.astype(np.float32)
            new=x[lo:stop]
            maxerr=max(maxerr,float(np.max(np.abs(old-new))))
            mismatches+=int(np.count_nonzero(old!=new))
    rows=[]
    for route in ['NMT_Thyroid_Seurat_none','TTU_Pubmed_UMAPHDBSCAN_mean']:
        d=parent/route
        assert (d/'TERMINAL_COMPLETE').exists()
        assert (d/'integration_features.txt').read_bytes()==(first/'integration_features.txt').read_bytes()
        assert (d/'cells.csv').read_bytes().splitlines()[0]==(first/'cells.csv').read_bytes().splitlines()[0]
        m=json.loads((d/'terminal_DL/training_manifest.json').read_text())
        score=json.loads((d/'score_manifest.json').read_text())
        assert score['DL_input_sha256']==json.loads((first/'score_manifest.json').read_text())['DL_input_sha256']
        z=np.load(d/'terminal_DL/terminal.npz',allow_pickle=False)
        initial=pd.read_csv(d/'initial_calls.csv',keep_default_na=False).initial.to_numpy(dtype=str)
        assert np.array_equal(z['initial'],initial)
        known=np.flatnonzero(initial!='Undecided');pool=np.flatnonzero(initial=='Undecided')
        assert np.array_equal(pool,z['pool_indices'])
        assert np.array_equal(z['final090'][known],initial[known])
        assert m['training_executed'] and m['params']['epochs']==10
        state=torch.load(d/'terminal_DL/model_state.pt',weights_only=True,map_location='cpu')
        model=independent_model(state,2000,len(z['classes']))
        prob=[]
        with torch.no_grad():
            for lo in range(0,len(pool),256):
                xx=np.array(x[pool[lo:lo+256]],dtype=np.float32,copy=True)
                prob.append(model(torch.from_numpy(xx)).numpy())
        prob=np.concatenate(prob)
        np.testing.assert_allclose(prob,z['probabilities'],rtol=1e-5,atol=1e-6)
        rounded=np.asarray([round(v,4) for v in z['probabilities'].max(1)],dtype=np.float32)
        names=z['classes'][z['probabilities'].argmax(1)]
        expected=np.where(rounded>=0.9,names,'Unknown')
        assert np.array_equal(expected,z['final090'][pool])
        rows.append(dict(route=route,n_cells=len(initial),n_known=len(known),n_pool=len(pool),
            max_saved_weight_forward_delta=float(np.max(np.abs(prob-z['probabilities']))),
            terminal_replay_exact=True))
    report=dict(job=os.environ['SLURM_JOB_ID'],matrix_shape=[92404,2000],
        loom_float32_input_max_difference=maxerr,loom_float32_differing_values=mismatches,
        source_cell_and_feature_order_exact=True,model_verification=rows,
        status='passed' if mismatches==0 else 'input_discrepancy',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (OUT/'input_and_terminal_verification.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
    assert mismatches==0

if __name__=='__main__':run()
