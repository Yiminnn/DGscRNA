"""A2 Lfine re-evaluation from immutable seed/terminal calls; never fit or remap."""
import argparse,json,os
from pathlib import Path
from common import *


def sample_run(sample,protocol,snapshot,provider):
    import numpy as np,pandas as pd
    assert sample in protocol['samples']
    record=snapshot['samples'][sample];verify_files(record['files'])
    cohort=pd.read_csv(REFERENCE/'protocol/cohort.csv');bool_column(cohort)
    row=cohort[cohort['sample']==sample];assert len(row)==1
    patient=str(row.patient.iloc[0]);primary=bool(row.primary.iloc[0])
    truth=pd.read_csv(REFERENCE/'evaluation_inputs'/sample/'truth.csv.gz',dtype=str,keep_default_na=False).rename(columns={'cell_id':'CellID','lfine_original':'Lfine'})
    assert truth.CellID.is_unique and truth.Patient.eq(patient).all()
    scope=provider.v5.label_scope(truth,provider.helpers)
    assert np.array_equal(scope[0],truth.Lfine.to_numpy())
    mapping=provider.mappings[FIXED]
    dest=OUT/'samples'/sample;dest.mkdir(parents=True,exist_ok=True)
    if checked(dest):
        m=read(dest/'manifest.json')
        assert m['protocol_sha256']==sha(OUT/'protocol_v2.json') and m['input_record']==record
        for name,digest in m['outputs'].items():assert sha(dest/name)==digest
        return
    metrics=[];parity=[];counts=[]
    compact=pd.read_csv(COMPACT/'metrics_hvg24.csv.gz');compact=compact[compact['sample']==sample]
    def evaluate_source(unit,budget,family,arm,lambda_value=None):
        source=Path(unit['directory']);score=read(source/'score_manifest.json')
        assert checked(source,'score_manifest.json','SCORE_COMPLETE')
        initial=pd.read_csv(source/'initial_calls.csv.gz',dtype=str,keep_default_na=False)
        assert np.array_equal(initial.cell_id,truth.CellID)
        assert sha(source/'initial_calls.csv.gz')==score['initial_sha256']
        aid=arm['id'];terminal=source/'terminal'/aid
        tm=read(terminal/'terminal_manifest.json') if (terminal/'terminal_manifest.json').is_file() else None
        base=dict(sample=sample,patient=patient,primary=primary,budget=budget,family=family,
            route='cellwise_seed' if family=='cellwise_seed' else source.name,library=FIXED,
            lambda_value=lambda_value,arm_id=aid,n_cells=len(truth),lfine_n_classes=len(scope[1]),
            score_manifest_sha256=sha(source/'score_manifest.json'),
            truth_sha256=sha(REFERENCE/'evaluation_inputs'/sample/'truth.csv.gz'),
            training_executed=tm['training_executed'] if tm else False,
            dl_status=tm['dl_status'] if tm else 'missing_terminal',
            terminal_valid=bool(tm and tm['terminal_valid'] and tm['status']=='completed'),
            terminal_manifest_sha256=sha(terminal/'terminal_manifest.json') if tm else None,
            original_prediction_source=str(terminal.relative_to(ROOT)))
        pred=None
        if tm:
            assert checked(terminal,'terminal_manifest.json','TERMINAL_COMPLETE')
            assert tm['score_manifest_sha256']==sha(source/'score_manifest.json')
            assert tm['DL_sha256']==score['DL_binary_sha256']
            assert tm['arm']['library']==FIXED
            if base['terminal_valid']:
                assert sha(terminal/'predictions.csv.gz')==tm['predictions_sha256']
                assert sha(terminal/'terminal.npz')==tm['terminal_sha256']
                pred=pd.read_csv(terminal/'predictions.csv.gz',dtype=str,keep_default_na=False)
                assert np.array_equal(pred.cell_id,truth.CellID) and np.array_equal(pred.initial,initial[aid])
                known=pred.initial.ne('Undecided')
                for column in ('final090','final070'):assert np.array_equal(pred.loc[known,column],pred.loc[known,'initial'])
                base['predictions_sha256']=tm['predictions_sha256']
        for stage,column in STAGES.items():
            item=dict(base,stage=stage,endpoint='v5 set-valued Lfine compatibility macro-F1')
            if stage!='initial' and pred is None:
                item.update(status='unavailable',unavailable_reason='missing_or_invalid_terminal_no_seed_substitution')
                metrics.append(item);continue
            values=initial[aid].to_numpy() if stage=='initial' else pred[column].to_numpy()
            semantic=np.asarray([mapping.get(v,'Unknown' if v in provider.ABSTAIN else 'UNMAPPABLE') for v in values],dtype=object)
            scored=provider.v5.metrics(semantic,truth,scope,provider.helpers,FIXED)
            item.update({key:scored[key] for key in FIELDS})
            assert item['n_reference_lfine_disagreements']==0
            item.update(status='completed' if len(scope[1]) else 'no_eligible_lfine_classes',
                initial_diagnostic_only=stage=='initial',source_native_labels_preserved=True)
            metrics.append(item)
            # Explicit confusion counts allow an independent metric reconstruction.
            lf,classes,targets=scope
            ok=np.asarray([gold in targets.get(call,()) for gold,call in zip(lf,semantic)])
            for label in classes:
                gold=lf==label;has=np.asarray([label in targets.get(call,()) for call in semantic])
                tp=int((gold&ok).sum());fn=int((gold&~ok).sum());fp=int((~ok&has&~gold).sum())
                counts.append(dict(sample=sample,budget=budget,family=family,route=item['route'],lambda_value=lambda_value,
                    stage=stage,lfine_class=label,support=int(gold.sum()),TP=tp,FP=fp,FN=fn))
            if family=='original_cluster_DEG' and stage=='terminal090':
                previous=compact[(compact.budget==budget)&(compact.route==source.name)&(compact.library==FIXED)]
                assert len(previous)==1
                for field in FIELDS:
                    assert np.isclose(item[field],previous.iloc[0][field],rtol=0,atol=1e-12,equal_nan=True),(sample,budget,source.name,field)
                assert previous.iloc[0].predictions_sha256==item['predictions_sha256']
                parity.append(dict(sample=sample,budget=budget,route=source.name,all_fields_match=True,
                    predictions_sha256=item['predictions_sha256'],lfine_macroF1=item['lfine_macroF1']))
    for budget in protocol['budgets']:
        unit=record['units'][budget];fitpath=OLD/'GBM'/sample/budget
        assert checked(fitpath,'fit_manifest.json','FIT_COMPLETE')
        fit=read(fitpath/'fit_manifest.json');cfg=read(fitpath/'config.json')
        assert fit['source_bundle_sha256']==protocol['original_A2_source_bundle_sha256']
        assert cfg['arms']==protocol['arms']
        assert fit['input_signature']==cfg['input_signature']
        for arm in protocol['arms']:
            terminal=Path(unit['cellwise']['directory'])/'terminal'/arm['id']
            if (terminal/'terminal_manifest.json').is_file():
                assert fit['terminal_manifests'][arm['id']]==sha(terminal/'terminal_manifest.json')
            evaluate_source(unit['cellwise'],budget,'cellwise_seed',arm,arm['lambda'])
        for route in ROUTES:
            spec=unit['reference'][route]
            evaluate_source(spec,budget,'original_cluster_DEG',dict(id=spec['arm_id']))
    frame=pd.DataFrame(metrics);bool_column(frame)
    assert len(frame)==54 and frame[frame.family.eq('cellwise_seed')].shape[0]==30
    assert not frame.duplicated(['sample','budget','family','route','lambda_value','stage']).any()
    frame.to_csv(dest/'metrics.csv',index=False)
    pd.DataFrame(counts).to_csv(dest/'per_class_counts.csv.gz',index=False,compression='gzip')
    pd.DataFrame(parity).to_csv(dest/'compact_anchor_parity.csv',index=False)
    assert len(parity)==8
    verify_files(record['files'])
    write(dest/'manifest.json',dict(status='completed' if frame.status.eq('completed').all() else 'completed_with_unavailable',
        sample=sample,n_rows=54,n_candidate_rows=30,n_reference_rows=24,n_cells=len(truth),
        all_cells_retained=True,no_fitting=True,no_expression_loading=True,
        protocol_sha256=sha(OUT/'protocol_v2.json'),source_manifest_sha256=sha(CODE/'SOURCE_MANIFEST.json'),
        input_record=record,outputs={p.name:sha(p) for p in dest.iterdir() if p.name in ('metrics.csv','per_class_counts.csv.gz','compact_anchor_parity.csv')},
        job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),completed_at=utc()))
    complete(dest);print('A2_LFINE_SAMPLE_COMPLETE',sample,len(frame),flush=True)


def main():
    p=verify();snapshot=read(OUT/'inputs.json');provider=load_provider()
    parser=argparse.ArgumentParser();parser.add_argument('--sample');parser.add_argument('--all',action='store_true');args=parser.parse_args()
    assert bool(args.sample)!=args.all
    requested=p['samples'] if args.all else [args.sample]
    approval=read(OUT/'APPROVED.json')
    assert set(requested)<=set(approval['allowed_samples']), 'Samples exceed reviewed execution phase'
    for sample in requested:sample_run(sample,p,snapshot,provider)
    verify()
if __name__=='__main__':main()
