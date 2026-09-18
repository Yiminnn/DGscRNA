"""Patient-level PTC assay endpoints, keeping historical and strict mappings apart."""
import json
import os
from pathlib import Path
import sys
from common import sha,write_json,complete,checked,utc
from ptc_followup_common import PTC,PAPER,ROUTE_MAP,require_ptc,read_reference,old_units
from ptc_label_rules import strict_T,broad_lineage,ptc_general

def binary_metrics(y,p,unknown):
    import numpy as np
    y=np.asarray(y,dtype=bool);p=np.asarray(p,dtype=bool);u=np.asarray(unknown,dtype=bool)
    assert not (p&u).any()
    tp=int((y&p).sum());fp=int((~y&p).sum());fn=int((y&~p).sum());tn=int((~y&~p).sum())
    tc=int((~y&~p&~u).sum());fc=int((y&~p&~u).sum());n=len(y);pos=int(y.sum());neg=n-pos
    safe=lambda a,b:a/b if b else 0.
    f1t=safe(2*tp,2*tp+fp+fn);f1n=safe(2*tc,2*tc+fc+neg-tc)
    return dict(n_cells=n,n_detected=pos,n_predicted_T=int(p.sum()),TP=tp,FP=fp,FN=fn,TN=tn,
        F1_T=f1t,F1_nonT=safe(2*tn,2*tn+fp+fn),F1_nonT_unknown_as_error=f1n,
        macro_F1_unknown_as_error=(f1t+f1n)/2,accuracy=(tp+tn)/n,accuracy_unknown_as_error=(tp+tc)/n,
        AUC_binary=.5*(tp/pos+tn/neg) if pos and neg else None,
        TCR_positive_recall=safe(tp,pos),TCR_detection_yield_in_predicted_T=safe(tp,int(p.sum())),
        coverage=float((~u).mean()),unknown_fraction=float(u.mean()),
        unknown_detected=int((y&u).sum()),unknown_undetected=int((~y&u).sum()))

def evaluate_source(source,ref,terminal_roots=None,context=None):
    require_ptc()
    import numpy as np
    import pandas as pd
    source=Path(source);m=json.loads((source/'score_manifest.json').read_text())
    assert checked(source,'score_manifest.json','SCORE_COMPLETE')
    cells=pd.read_csv(source/'cells.csv',dtype=str,keep_default_na=False)
    assert cells.cell_id.is_unique and set(cells.cell_id)<=set(ref.index)
    metadata=ref.loc[cells.cell_id];rows=[];states=[]
    names=set(json.loads((PAPER/'historical_T_names_from_vignette.json').read_text()))
    scopes=[]
    for group,gm in metadata.groupby('group',sort=True):
        scopes.append((group,'ALL',group,metadata.index.isin(gm.index)))
        for patient,pm in gm.groupby('patient',sort=True):
            assert pm['sample'].nunique()==1
            scopes.append((group,patient,pm['sample'].iloc[0],metadata.index.isin(pm.index)))
    if terminal_roots is None:terminal_roots={aid:source/'terminal'/aid for aid in m['arms']}
    common=dict(unit=source.parent.name,route=ROUTE_MAP.get(source.name,source.name),**(context or {}))
    for aid,td in terminal_roots.items():
        td=Path(td);arm=m['arms'][aid]
        tm=json.loads((td/'training_manifest.json').read_text())
        if (td/'terminal_manifest.json').exists():
            assert checked(td,'terminal_manifest.json','TERMINAL_COMPLETE')
            terminal_manifest=json.loads((td/'terminal_manifest.json').read_text())
            assert terminal_manifest['score_manifest_sha256']==sha(source/'score_manifest.json')
        else:
            assert (td/'TERMINAL_COMPLETE').read_text().strip()==sha(td/'training_manifest.json')
            assert tm['provenance']['score_manifest_sha256']==sha(source/'score_manifest.json')
        with np.load(td/'terminal.npz',allow_pickle=False) as z:
            initial=z['initial'];known=initial!='Undecided'
            assert len(initial)==len(metadata)
            assert np.array_equal(z['final090'][known],initial[known])
            for group,patient,scope,mask in scopes:
                lab=initial[mask];unique,count=np.unique(lab,return_counts=True)
                t_seeds=sum(int(n) for v,n in zip(unique,count) if strict_T(str(v)))
                states.append(dict(**common,arm_id=aid,library=arm['library'],cutoff=arm['cutoff'],
                    group=group,patient=patient,scope=scope,n_cells=int(mask.sum()),n_known=int((lab!='Undecided').sum()),
                    n_strict_T_seeds=t_seeds,dl_status=tm['dl_status'],training_executed=tm['training_executed'],
                    terminal_directory=str(td),training_manifest_sha256=sha(td/'training_manifest.json')))
            for stage in ['initial','final090','final070']:
                lab=z[stage];unique,index=np.unique(lab,return_inverse=True)
                strict=np.asarray([strict_T(str(v)) for v in unique],dtype=bool)[index]
                broad=np.asarray([ptc_general(str(v)) in names for v in unique],dtype=bool)[index]
                unknown=np.isin(lab,['Unknown','Undecided','No_Annotation',''])
                same=lab==metadata.paper_native.to_numpy(dtype=str)
                endpoints=[('productive_TCR','strict_T_name_rule',metadata.TCR_cell_high_confidence_productive_TCR.to_numpy(),strict),
                    ('paper_original','paper_broad_T_compatibility',metadata.paper_truth.to_numpy(),broad),
                    ('paper_original','strict_T_name_rule',metadata.paper_truth.to_numpy(),strict),
                    ('any_filtered_contig','strict_T_name_rule',metadata.TCR_any_filtered_contig.to_numpy(),strict),
                    ('S3_supplied','strict_T_name_rule',metadata.TCR_S3_supplied.to_numpy(),strict)]
                for truth,endpoint,y,p in endpoints:
                    for group,patient,scope,mask in scopes:
                        rows.append(dict(**common,arm_id=aid,library=arm['library'],cutoff=arm['cutoff'],stage=stage,
                            group=group,patient=patient,scope=scope,truth_definition=truth,endpoint=endpoint,
                            saved_native_concordance=float(same[mask].mean()),dl_status=tm['dl_status'],
                            **binary_metrics(y[mask],p[mask],unknown[mask])))
    return rows,states

def run_old(prep):
    require_ptc()
    import numpy as np
    import pandas as pd
    prep=Path(prep);dest=PTC/'existing_grid_evaluation'/prep.name;dest.mkdir(parents=True,exist_ok=True)
    if checked(dest):return
    ref=read_reference();rows=[];states=[];sources={}
    paths=sorted(prep.glob('*/score_manifest.json'));assert len(paths)==4
    for path in paths:
        r,s=evaluate_source(path.parent,ref);rows.extend(r);states.extend(s);sources[str(path)]=sha(path)
    frame=pd.DataFrame(rows);original=pd.read_csv(prep/'evaluation/metrics.csv.gz')
    original['route']=original.route.map(lambda v:ROUTE_MAP.get(v,v))
    keys=['route','arm_id','stage','scope','endpoint']
    got=frame[frame.truth_definition.eq('paper_original')]
    merged=got.merge(original,on=keys,validate='one_to_one',suffixes=('_new','_old'))
    assert len(merged)==len(got),(len(merged),len(got))
    for metric in ['F1_T','F1_nonT','AUC_binary','unknown_fraction','accuracy']:
        np.testing.assert_allclose(merged[metric+'_new'],merged[metric+'_old'],atol=1e-12,rtol=0,equal_nan=True)
    frame.to_csv(dest/'metrics.csv.gz',index=False);pd.DataFrame(states).to_csv(dest/'terminal_statuses.csv.gz',index=False)
    write_json(dest/'manifest.json',dict(status='complete',unit=prep.name,n_metric_rows=len(frame),
        historical_metrics_parity_rows=len(merged),historical_metric_max_tolerance=1e-12,
        sources=sources,reference_labels_used_for_fit=False,refit_performed=False,
        TCR_absence_is_not_established_nonT=True,rule_sha256=sha(Path(__file__).with_name('ptc_label_rules.py')),
        job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest);print('PTC_OLD_EVALUATION',prep.name,len(frame),flush=True)

if __name__=='__main__':
    if sys.argv[1]=='array':run_old(old_units()[int(os.environ['SLURM_ARRAY_TASK_ID'])])
    else:run_old(sys.argv[1])
