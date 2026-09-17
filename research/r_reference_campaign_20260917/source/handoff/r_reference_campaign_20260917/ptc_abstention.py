"""Complement historical binary PTC metrics with explicit unknown-as-error metrics."""
import hashlib,json,os,sys
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'
sys.path.insert(0,str(ROOT/'handoff/ptc_recovery_20260916'))
from label_rules import strict_T
from evaluation_rules import ptc_general

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import numpy as np
    import pandas as pd
    from functools import lru_cache
    strict=lru_cache(maxsize=65536)(strict_T)
    prep=Path(sys.argv[1]);dest=prep/'evaluation';assert (dest/'COMPLETE').exists()
    ref=pd.read_csv(ROOT/'results/hvg_ptc_20260916_v1/ptc_paper_baseline/paper_baseline_reference.csv.gz',keep_default_na=False).set_index('cell_id')
    truth=pd.read_csv(ROOT/'results/hvg_ptc_20260916_v1/ptc_paper_baseline/original_DG_binary_pairs_for_R.csv.gz').set_index('cell_id')
    names=set(json.loads((ROOT/'results/hvg_ptc_20260916_v1/ptc_paper_baseline/historical_T_names_from_vignette.json').read_text()))
    rows=[];sources=[]
    for path in sorted(prep.glob('*/score_manifest.json')):
        source=path.parent;m=json.loads(path.read_text());sources.append(hashlib.sha256(path.read_bytes()).hexdigest())
        ids=pd.read_csv(source/'cells.csv',dtype=str).cell_id
        metadata=ref.loc[ids];y=truth.loc[ids,'truth'].to_numpy(dtype=int)
        scopes=[('all',np.ones(len(y),dtype=bool))]
        for col in ['group','sample']:
            scopes += [(value,metadata[col].eq(value).to_numpy()) for value in metadata[col].unique()]
        for aid,arm in m['arms'].items():
            with np.load(source/'terminal'/aid/'terminal.npz',allow_pickle=False) as z:
                for stage in ['initial','final090','final070']:
                    labels=z[stage];unique,index=np.unique(labels,return_inverse=True)
                    missing=np.isin(labels,['Undecided','Unknown','No_Annotation'])
                    for endpoint,called in [('strict_T_name_rule',np.asarray([strict(str(v)) for v in unique],dtype=bool)[index]),
                        ('paper_broad_T_compatibility',np.asarray([ptc_general(str(v)) in names for v in unique],dtype=bool)[index])]:
                        for scope,mask in scopes:
                            yy=y[mask];pp=called[mask];unknown=missing[mask];resolved=~unknown
                            tp=int(((yy==1)&pp&resolved).sum());fp=int(((yy==0)&pp&resolved).sum())
                            tn=int(((yy==0)&~pp&resolved).sum());fn_called=int(((yy==1)&~pp&resolved).sum())
                            fn_T=int((yy==1).sum())-tp
                            fn_nonT=int((yy==0).sum())-tn
                            f1T=2*tp/(2*tp+fp+fn_T) if 2*tp+fp+fn_T else 0.0
                            f1N=2*tn/(2*tn+fn_called+fn_nonT) if 2*tn+fn_called+fn_nonT else 0.0
                            n=len(yy);ncalled=int(resolved.sum())
                            rows.append(dict(dataset='PTC',unit=prep.name,route=source.name,arm_id=aid,
                              library=arm['library'],cutoff=arm['cutoff'],stage=stage,scope=scope,endpoint=endpoint,
                              n_cells=n,n_called=ncalled,coverage=ncalled/n,
                              accuracy_unknown_as_error=(tp+tn)/n,
                              accuracy_among_called=(tp+tn)/ncalled if ncalled else None,
                              F1_T_unknown_as_error=f1T,F1_nonT_unknown_as_error=f1N,
                              macro_F1_T_nonT_unknown_as_error=(f1T+f1N)/2,
                              unknown_true_T=int((unknown&(yy==1)).sum()),unknown_true_nonT=int((unknown&(yy==0)).sum()),
                              TP_called=tp,TN_called=tn,FP_called=fp,FN_called=fn_called))
    frame=pd.DataFrame(rows);frame.to_csv(dest/'abstention_aware_metrics.csv.gz',index=False)
    original=pd.read_csv(dest/'metrics.csv.gz')
    keys=['unit','route','arm_id','stage','scope','endpoint']
    comparison=frame.merge(original,on=keys,validate='one_to_one',suffixes=('_new','_old'))
    np.testing.assert_allclose(comparison.F1_T_unknown_as_error,comparison.F1_T,atol=1e-12)
    assert (comparison.accuracy_unknown_as_error<=comparison.accuracy+1e-12).all()
    report=dict(status='passed',unit=prep.name,job=os.environ['SLURM_JOB_ID'],rows=len(frame),
      unknown_policy='Unknown is neither a called T nor a called non-T; it is an error in accuracy and a false negative for its true class.',
      T_positive_F1_unchanged_from_binary_mapping=True,source_manifests_sha256=sources,
      script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (dest/'abstention_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    (dest/'ABSTENTION_COMPLETE').write_text(hashlib.sha256((dest/'abstention_audit.json').read_bytes()).hexdigest()+'\n')
    print('ABSTENTION_AUDITED',prep,len(frame),flush=True)

if __name__=='__main__':run()
