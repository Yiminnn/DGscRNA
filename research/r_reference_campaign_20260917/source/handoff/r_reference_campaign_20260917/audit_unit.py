"""Verify inputs, every saved terminal artifact, and threshold/retained-label invariants."""
import hashlib,json,os,sys
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import numpy as np
    import pandas as pd
    prep=Path(sys.argv[1]);dataset=sys.argv[2]
    dest=prep/'verification';dest.mkdir(exist_ok=True)
    assert (prep/'evaluation/COMPLETE').exists()
    paths=sorted(prep.glob('*/score_manifest.json'));assert len(paths)==4
    checked_binary={};rows=[];order_checks=[]
    ids=None
    for path in paths:
        source=path.parent;m=json.loads(path.read_text())
        assert (source/'SCORE_COMPLETE').read_text().strip()==sha(path)
        if m.get('execution_source'):
            assert sha(m['execution_source'])==m['source_sha256']
        assert sha(source/'initial_calls.csv.gz')==m['initial_sha256']
        cells=pd.read_csv(source/'cells.csv',dtype=str,keep_default_na=False)
        initial=pd.read_csv(source/'initial_calls.csv.gz',dtype=str,keep_default_na=False)
        if ids is None:ids=cells.cell_id.to_numpy()
        assert np.array_equal(ids,cells.cell_id) and np.array_equal(ids,initial.cell_id)
        assert len(set(ids))==len(ids)==m['n_cells']
        binary=Path(m['DL_binary'])
        if str(binary) not in checked_binary:
            assert binary.stat().st_size==len(ids)*m['DL_features']*4
            checked_binary[str(binary)]=sha(binary)
        assert checked_binary[str(binary)]==m['DL_binary_sha256']
        for aid,arm in m['arms'].items():
            td=source/'terminal'/aid
            tm=json.loads((td/'training_manifest.json').read_text())
            for flag in ['COMPLETE','TERMINAL_COMPLETE']:assert (td/flag).read_text().strip()==sha(td/'training_manifest.json')
            assert tm['provenance']['score_manifest_sha256']==sha(path)
            assert not tm['provenance']['evaluation_labels_used_for_fit']
            for filename,value in tm['outputs'].items():assert sha(td/filename)==value,(td,filename)
            with np.load(td/'terminal.npz',allow_pickle=False) as z:
                before=initial[arm['seed_column']].to_numpy(dtype=str)
                assert np.array_equal(z['initial'],before)
                known=before!='Undecided';pool=np.flatnonzero(~known)
                assert known.sum()==tm['n_known'] and len(pool)==tm['n_pool']
                assert np.array_equal(z['pool_indices'],pool)
                assert np.array_equal(z['final090'][known],before[known])
                assert np.array_equal(z['final070'][known],before[known])
                if tm['training_executed']:
                    assert len(set(z['train_indices']) & set(z['validation_indices']))==0
                    assert set(z['train_indices'])|set(z['validation_indices'])==set(np.flatnonzero(known))
                    p=z['probabilities'];assert p.shape==(len(pool),len(z['classes'])) and np.isfinite(p).all()
                    np.testing.assert_allclose(p.sum(1),1,atol=1e-5)
                    rounded=np.asarray([round(v,4) for v in p.max(1)],dtype=np.float32)
                    assert np.array_equal(rounded,z['confidence_rounded'][pool])
                    called=z['classes'][p.argmax(1)]
                    for threshold,key in [(.9,'final090'),(.7,'final070')]:
                        expected=before.copy();expected[pool]=np.where(rounded>=threshold,called,'Unknown')
                        assert np.array_equal(expected,z[key])
                else:assert len(z['probabilities'])==0
                for key,field in [('final090','n_final_called090'),('final070','n_final_called070')]:
                    assert int((~np.isin(z[key],['Unknown','Undecided'])).sum())==tm[field]
                predictions=pd.read_csv(td/'predictions.csv.gz',dtype=str,keep_default_na=False)
                assert np.array_equal(predictions.cell_id,ids)
                for key in ['initial','final090','final070']:assert np.array_equal(predictions[key],z[key])
            rows.append(dict(route=source.name,arm_id=aid,dl_status=tm['dl_status'],artifact_hashes_match=True,initial_known_retained=True,threshold_reconstruction_exact=True))
    if dataset=='PTC':
        if prep.name=='PTC_archived_CCA2000' or prep.name.startswith('PTC_ALL8'):
            expected=pd.read_csv(ROOT/'results/hvg_ptc_20260916_v1/ptc_paper_baseline/replay_selected_routes_full_parallel/NMT_Thyroid_Seurat_none/cells.csv').cell_id.to_numpy()
        else:
            group='MTN' if '_NMT_' in prep.name else 'TUT'
            expected=pd.read_csv(ROOT/'results/hvg_ptc_20260916_v1/ptc_experiments/prepared'/group/'cells.csv',index_col=0).index.to_numpy()
        assert np.array_equal(ids,expected),'PTC fixed cell-order contract violated'
        order_checks.append('PTC canonical cell order exact')
    else:
        expected=pd.read_csv(OUT/'inputs'/dataset/prep.parent.name/'cells_fit.csv',dtype=str).cell_id.to_numpy()
        assert set(ids)==set(expected),'Reviewer cell set changed'
        order_checks.append('All prepared reviewer cells retained exactly once')
    pd.DataFrame(rows).to_csv(dest/'terminal_artifact_audit.csv',index=False)
    report=dict(status='passed',job=os.environ['SLURM_JOB_ID'],dataset=dataset,unit=str(prep),
      terminal_arms_audited=len(rows),clustering_routes=4,cell_checks=order_checks,DL_binary_hashes=checked_binary,
      all_model_npz_history_hashes_match=True,all_prediction_csvs_match=True,evaluation_labels_excluded_from_fit=True,
      script_sha256=sha(__file__))
    (dest/'audit_manifest.json').write_text(json.dumps(report,indent=2)+'\n')
    (dest/'AUDIT_COMPLETE').write_text(sha(dest/'audit_manifest.json')+'\n')
    print('AUDITED',prep,len(rows),flush=True)

if __name__=='__main__':run()
