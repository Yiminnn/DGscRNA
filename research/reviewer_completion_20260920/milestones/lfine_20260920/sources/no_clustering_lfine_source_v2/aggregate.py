"""Frozen-fold selection and original-route comparisons after every A2 unit exists."""
import json
import os
from common import *
write_json=write


def run():
    protocol=verify()
    assert read(OUT/'APPROVED.json').get('allow_aggregate') is True, 'Aggregation awaits reviewed full-evaluation phase'
    import numpy as np
    import pandas as pd
    from scipy.stats import wilcoxon
    from statsmodels.stats.multitest import multipletests
    dest=OUT/'summary';dest.mkdir(exist_ok=True)
    folds_path=OLD/'patient_folds.csv'
    assert sha(folds_path)==protocol['patient_folds_sha256']
    folds=pd.read_csv(folds_path)
    patientfold=folds[['patient','fold']].drop_duplicates()
    assert not patientfold.patient.duplicated().any()
    measurements=[];reference=[];hashes={}
    for sample in protocol['samples']:
        unit=OUT/'samples'/sample
        assert checked(unit),f'Incomplete Lfine sample: {sample}'
        em=read(unit/'manifest.json')
        assert em['protocol_sha256']==sha(OUT/'protocol_v2.json')
        assert em['source_manifest_sha256']==sha(CODE/'SOURCE_MANIFEST.json')
        for name,digest in em['outputs'].items():assert sha(unit/name)==digest
        frame=pd.read_csv(unit/'metrics.csv');bool_column(frame)
        assert len(frame)==54
        measurements.append(frame[frame.family.eq('cellwise_seed')])
        reference.append(frame[frame.family.eq('original_cluster_DEG')])
        hashes[str(unit/'manifest.json')]=sha(unit/'manifest.json')
    metrics=pd.concat(measurements,ignore_index=True)
    reference_all=pd.concat(reference,ignore_index=True)
    old=reference_all[reference_all.stage.eq('terminal090')]
    bool_column(metrics);bool_column(old)
    assert len(metrics)==3630 and len(reference_all)==2904 and len(old)==968
    assert set(metrics.lambda_value)=={0,.5,1,1.5,2}
    # Never drop an invalid/no-class/missing sample and quietly rank a subset.
    unavailable=pd.concat([metrics[~metrics.status.eq('completed')],reference_all[~reference_all.status.eq('completed')]])
    unavailable.to_csv(dest/'unavailable_conditions.csv',index=False)
    if len(unavailable) or not metrics.lfine_macroF1.notna().all() or not reference_all.lfine_macroF1.notna().all():
        write(dest/'INCOMPLETE.json',dict(status='not_ranked',reason='Unavailable endpoint values retained; no silent exclusion',
            n_unavailable=len(unavailable),time=utc()))
        raise RuntimeError('Resolve unavailable endpoints before complete paired inference')
    reference_all.to_csv(dest/'all_original_route_metrics.csv',index=False)
    metrics.to_csv(dest/'all_candidate_metrics.csv',index=False)
    metrics.groupby(['budget','lambda_value','stage','dl_status','training_executed']).size().rename('conditions').reset_index().to_csv(dest/'execution_states.csv',index=False)
    folds=folds[['sample','patient','fold']]
    terminal=metrics[metrics.stage=='terminal090'].merge(folds,on=['sample','patient'],validate='many_to_one')
    choices=[];heldout=[];pairs=[];summary=[]
    rng=np.random.default_rng(42)
    for cohort in ['primary','all']:
        subset=terminal[terminal.primary.eq(True)] if cohort=='primary' else terminal
        baseline=old[old.primary.eq(True)] if cohort=='primary' else old
        expected_samples=protocol['primary_n_samples'] if cohort=='primary' else len(protocol['samples'])
        assert subset['sample'].nunique()==baseline['sample'].nunique()==expected_samples
        for budget in protocol['budgets']:
            data=subset[subset.budget==budget]
            pat=data.groupby(['patient','fold','lambda_value'])[['lfine_macroF1','coverage']].mean().reset_index()
            selected=[]
            for fold in sorted(pat.fold.unique()):
                train=pat[pat.fold!=fold];test=pat[pat.fold==fold]
                scores=train.groupby('lambda_value').lfine_macroF1.mean().reset_index(name='training_patient_lfine_macroF1')
                assert len(scores)==5
                scores['distance_from_primary']=(scores.lambda_value-1).abs()
                rank=scores.sort_values(['training_patient_lfine_macroF1','distance_from_primary','lambda_value'],ascending=[False,True,True],kind='stable')
                chosen=rank.iloc[0]
                choices.append(dict(cohort=cohort,budget=budget,fold=int(fold),lambda_value=float(chosen.lambda_value),
                    training_patient_lfine_macroF1=float(chosen.training_patient_lfine_macroF1),n_training_patients=train.patient.nunique(),
                    n_test_patients=test.patient.nunique(),test_labels_used_for_selection=False,
                    all_training_candidate_scores=scores.to_json(orient='records')))
                selected.append(data[(data.fold==fold)&(data.lambda_value==chosen.lambda_value)].assign(selection='training_patient_selected',cohort=cohort))
            selections=[data[data.lambda_value==1].assign(selection='fixed_lambda1',cohort=cohort),pd.concat(selected)]
            for selected in selections:
                assert selected['sample'].nunique()==expected_samples and len(selected)==expected_samples
                heldout.append(selected)
                patient=selected.groupby('patient')[['lfine_macroF1','coverage']].mean()
                mode=selected.selection.iloc[0]
                for route,group in baseline[baseline.budget==budget].groupby('route'):
                    orig=group.groupby('patient')[['lfine_macroF1','coverage']].mean()
                    match=patient.join(orig,how='outer',lsuffix='_cellwise',rsuffix='_original')
                    assert not match.isna().any().any()
                    delta=(match.lfine_macroF1_cellwise-match.lfine_macroF1_original).to_numpy()
                    boot=delta[rng.integers(0,len(delta),size=(2000,len(delta)))].mean(axis=1)
                    p=1.0 if np.all(delta==0) else float(wilcoxon(delta,zero_method='pratt',alternative='two-sided',method='auto').pvalue)
                    summary.append(dict(cohort=cohort,budget=budget,selection=mode,reference_route=route,n_patients=len(delta),
                        mean_cellwise_lfine_macroF1=float(match.lfine_macroF1_cellwise.mean()),
                        mean_original_lfine_macroF1=float(match.lfine_macroF1_original.mean()),
                        mean_delta_cellwise_minus_original=float(delta.mean()),CI025=float(np.quantile(boot,.025)),CI975=float(np.quantile(boot,.975)),
                        p_value=p,mean_coverage_delta=float((match.coverage_cellwise-match.coverage_original).mean()),
                        interval_scope='conditional on frozen predictions and selected configurations; no refitting'))
                    for patient_id,row in match.iterrows():pairs.append(dict(cohort=cohort,budget=budget,selection=mode,reference_route=route,patient=patient_id,**row.to_dict()))
    summary=pd.DataFrame(summary)
    for cohort,indices in summary.groupby('cohort').groups.items():
        assert len(indices)==16
        summary.loc[indices,'p_holm']=multipletests(summary.loc[indices,'p_value'],method='holm')[1]
    pd.DataFrame(choices).to_csv(dest/'training_patient_lambda_choices.csv',index=False)
    pd.concat(heldout,ignore_index=True).to_csv(dest/'selected_and_fixed_sample_results.csv',index=False)
    pd.DataFrame(pairs).to_csv(dest/'patient_paired_results.csv',index=False)
    summary.to_csv(dest/'patient_paired_summary.csv',index=False)
    report='''# A2 Lfine endpoint revision — HVG2000/5000 only

The user requested Lfine-only presentation after the original L1 evaluation.
This is a separately frozen, retrospective endpoint revision. Existing fits,
native predictions, marker mappings, and original L1 results remain unchanged.
No result is relabelled as a newly fitted all-gene model.

The exact v5/compact set-valued Lfine compatibility macro-F1 averages reference
classes with at least20cells excluding Other/nan, with every cell retained in
TP/FP/FN and target sets built from every observed class. Unknown/Undecided and
unmappable predictions receive no compatible-class credit. Broad native panels
are not presented as one-to-one fine-state predictions.

Fixedλ1 and training-patient-selectedλ are separate. Selection now maximizes
training-patient mean terminal090 Lfine score on the same patientfolds; exact
ties choose nearest1 then smallerλ. Held-out patients never select their λ.
Initial calls are diagnostic only; terminal070 is retained as sensitivity.
All no-op/invalid execution states remain explicit. An unavailable endpoint
blocks a complete ranking instead of silently dropping a condition or patient.

Each contrast uses the same marker, budget, cells and endpoint in four original
cluster-DEG routes. Patient bootstrap2000, seed42, WilcoxonPratt/two-sided and
Holm16 contrasts/cohort are unchanged from original A2. Intervals condition on
saved predictions/selected configurations and exclude refitting/selection
uncertainty. All32 contrasts are reported; no outcome-selected best result.

This is an explicit seed-mechanism replacement, not pure clustering deletion,
not a proof of global optimality, and not part of the all-gene main baseline.
'''
    (dest/'REPORT.md').write_text(report)
    files=[p for p in dest.iterdir() if p.is_file() and p.name not in ['manifest.json','COMPLETE']]
    write_json(dest/'manifest.json',dict(status='completed',source_manifest_sha256=protocol['source_manifest_sha256'],
        protocol_sha256=sha(OUT/'protocol_v2.json'),endpoint=protocol['endpoint'],
        original_L1_results_preserved=True,no_model_fitting=True,HVG2000_5000_only=True,
        n_samples=len(protocol['samples']),n_budgets=2,n_candidate_conditions=1210,selection_modes=['fixed_lambda1','training_patient_selected'],
        n_sample_result_rows_by_cohort=pd.concat(heldout).groupby('cohort').size().astype(int).to_dict(),
        every_candidate_retained=True,selection_uses_training_patients_only=True,patient_folds_sha256=sha(folds_path),
        n_pairs=len(pairs),n_contrasts=len(summary),input_manifest_hashes=hashes,
        outputs={p.name:sha(p) for p in files},job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    verify();complete(dest);print('A2_LFINE_AGGREGATION_COMPLETE',len(summary),flush=True)


if __name__=='__main__':run()
