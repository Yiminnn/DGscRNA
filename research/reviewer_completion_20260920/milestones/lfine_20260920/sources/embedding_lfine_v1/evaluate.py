"""Evaluate one sample/budget's seven frozen A1 spaces and original-R anchor."""
from pathlib import Path
import argparse,os
import common as c

def run(index,space):
    c.require_slurm();contract=c.contract();links=c.archive_gate(contract)
    import numpy as np
    import pandas as pd
    from endpoint import initialize,score
    task=c.js(c.OUT/'shards.json')[index];sample,budget=task['sample'],task['budget']
    assert space in task['spaces']
    out=c.OUT/'evaluation'/sample/budget/space
    if c.checked(out):return c.checked_outputs(out)
    out.mkdir(parents=True,exist_ok=True);lock=out/'ACTIVE';lock.mkdir()
    c.write(lock/'owner.json',dict(job=os.environ['SLURM_JOB_ID'],pid=os.getpid()))
    try:
        helpers,exact,lookup=initialize()
        tp=c.NATIVE/'evaluation_inputs'/sample/'truth.csv.gz';impath=c.NATIVE/'inputs'/sample/'input_manifest.json';im=c.js(impath)
        assert c.sha(tp)==im['evaluation_files']['truth.csv.gz']
        truth=pd.read_csv(tp,dtype=str,keep_default_na=False).rename(columns={'cell_id':'CellID','lfine_original':'Lfine'})
        assert truth.CellID.is_unique and truth.Patient.eq(im['patient']).all() and type(im['primary']) is bool
        scope=exact['label_scope'](truth,helpers);assert np.array_equal(scope[0],truth.Lfine.to_numpy())
        inputs={str(tp):c.sha(tp),str(impath):c.sha(impath),str(c.ARCHIVE/'contract.json'):c.sha(c.ARCHIVE/'contract.json')}
        counts=pd.Series(scope[0]).value_counts()
        eligibility=dict(sample=sample,patient=im['patient'],primary=im['primary'],n_cells=len(truth),
            n_scored_classes=len(scope[1]),lfine_metric_eligible=bool(scope[1]),truth_sha256=c.sha(tp),
            class_support={str(k):int(v) for k,v in counts.items()},macro_classes=list(scope[1]),
            compatible_targets={k:sorted(v) for k,v in scope[2].items()},all_cells_retained=True)
        rows=[];classes=[];anchors=[]
        def consume(route,context,collection,accepted_artifacts=None):
            smp=route/'score_manifest.json';assert c.checked(route,'score_manifest.json','SCORE_COMPLETE')
            sm=c.js(smp);inputs[str(smp)]=c.sha(smp)
            if accepted_artifacts is not None:assert accepted_artifacts[str(smp)]==c.sha(smp)
            arms=[k for k,v in sm['arms'].items() if v['library']==c.FIXED and str(v['cutoff'])=='mean'];assert len(arms)==1
            td=route/'terminal'/arms[0];tmp=td/'terminal_manifest.json'
            base=dict(sample=sample,patient=im['patient'],primary=im['primary'],budget=budget,library=c.FIXED,cutoff='mean',
                arm_id=arms[0],family='native_R_budget',n_cells=len(truth),lfine_n_classes=len(scope[1]),
                lfine_metric_eligible=bool(scope[1]),truth_sha256=c.sha(tp),**context)
            if not c.checked(td,'terminal_manifest.json','TERMINAL_COMPLETE'):
                for stage in c.STAGES:collection.append(dict(**base,stage=stage,status='unavailable',terminal_valid=False,unavailable_reason='missing_terminal_completion'))
                return
            tm=c.js(tmp);assert tm['score_manifest_sha256']==c.sha(smp)
            if accepted_artifacts is not None:assert accepted_artifacts[str(tmp)]==c.sha(tmp)
            inputs[str(tmp)]=c.sha(tmp)
            base.update(dl_status=tm['dl_status'],training_executed=tm['training_executed'],DL_features=tm['DL_features'],
                n_known=tm['n_known'],n_pool=tm['n_pool'],n_training_classes=tm['n_training_classes'],terminal_manifest_sha256=c.sha(tmp))
            if tm['status']!='completed' or tm['terminal_valid'] is not True:
                for stage in c.STAGES:collection.append(dict(**base,stage=stage,status='invalid',terminal_valid=False,unavailable_reason='invalid_terminal_state'))
                return
            pp=td/'predictions.csv.gz';zp=td/'terminal.npz'
            assert c.sha(pp)==tm['predictions_sha256'] and c.sha(zp)==tm['terminal_sha256']
            if accepted_artifacts is not None:
                assert accepted_artifacts[str(pp)]==c.sha(pp) and accepted_artifacts[str(zp)]==c.sha(zp)
            inputs[str(pp)]=c.sha(pp);inputs[str(zp)]=c.sha(zp)
            pred=pd.read_csv(pp,dtype=str,keep_default_na=False)
            assert np.array_equal(pred.cell_id,truth.CellID)
            with np.load(zp,allow_pickle=False) as z:
                for column in c.STAGES.values():assert np.array_equal(pred[column],z[column])
            known=pred.initial.ne('Undecided')
            for column in ['final090','final070']:assert np.array_equal(pred.loc[known,'initial'],pred.loc[known,column])
            for stage,column in c.STAGES.items():
                values,perclass=score(pred[column].to_numpy(),truth,scope,helpers,exact,lookup)
                context2=dict(**base,stage=stage,status='completed',terminal_valid=True,predictions_sha256=c.sha(pp),
                    predictions_path=str(pp),source_native_labels_preserved=True)
                scored=dict(context2);scored.update(values);collection.append(scored)
                for row in perclass:classes.append(dict(sample=sample,budget=budget,space=context['space'],route=context['route'],stage=stage,**row))
        for space in [space]:
            proof=unit_accept=c.unit_acceptance(sample,budget,space,links);inputs[str(proof)]=c.sha(proof)
            accepted_artifacts=c.js(proof.parent/'source_artifacts.json')
            inputs[str(proof.parent/'source_artifacts.json')]=c.sha(proof.parent/'source_artifacts.json')
            for path in [tp,impath,c.NATIVE/'markers/panel_L1_mapping.csv']:
                assert accepted_artifacts[str(path)]==c.sha(path),path
            unit=c.EMBED/sample/budget/space;cfg=c.js(unit/'config.json');rep=c.js(unit/'representation.json')
            for name in ['config.json','representation.json','fit_manifest.json']:inputs[str(unit/name)]=c.sha(unit/name)
            for condition in cfg['conditions']:
                consume(Path(condition['dest']),dict(space=space,space_display='ICA2 adaptive' if space=='ICA2' else space,
                    route=Path(condition['dest']).name,method=condition['method'],k=condition['k'] if condition['k'] is not None else 0,
                    actual_solver=rep.get('actual_solver',rep.get('params',{}).get('algorithm')) if space=='ICA2' else None),rows,accepted_artifacts)
        consume(c.NATIVE/'GBM'/sample/budget/'UMAP2_HDBSCAN_R',dict(space='original_R_anchor',space_display='Original R PCA30 to UMAP2',
            route='UMAP2_HDBSCAN_R',method='HDBSCAN_R',k=0,actual_solver=None),anchors)
        data=pd.DataFrame(rows);anchor=pd.DataFrame(anchors);perclass=pd.DataFrame(classes)
        assert len(data)==39 and len(anchor)==3 and not data.duplicated(['space','route','stage']).any()
        for frame in [data,anchor]:
            for field in c.FIELDS:
                if field not in frame:frame[field]=np.nan
        data.to_csv(out/'metrics.csv.gz',index=False);anchor.to_csv(out/'anchor_metrics.csv',index=False)
        perclass.to_csv(out/'per_class.csv.gz',index=False);c.write(out/'eligibility.json',eligibility)
        for path,digest in inputs.items():assert c.sha(path)==digest,path
        manifest=dict(status='completed' if data.terminal_valid.all() and anchor.terminal_valid.all() else 'completed_with_unavailable',
            sample=sample,budget=budget,space=space,n_rows=len(data),n_anchor_rows=3,n_valid=int(data.terminal_valid.sum()),
            contract_sha256=c.sha(c.OUT/'contract.json'),inputs=inputs,outputs={name:c.sha(out/name) for name in ['metrics.csv.gz','anchor_metrics.csv','per_class.csv.gz','eligibility.json']},
            all_cells_retained=True,endpoint=contract['endpoint'],no_fitting=True,job=os.environ['SLURM_JOB_ID'],completed_at=c.utc())
        c.write(out/'manifest.json',manifest);c.complete(out);return manifest
    finally:
        (lock/'owner.json').unlink();lock.rmdir()

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--index',type=int,default=int(os.environ.get('SLURM_ARRAY_TASK_ID','-1')));p.add_argument('--space');a=p.parse_args()
    assert 0<=a.index<242
    for space in ([a.space] if a.space else c.js(c.OUT/'shards.json')[a.index]['spaces']):print(run(a.index,space),flush=True)
