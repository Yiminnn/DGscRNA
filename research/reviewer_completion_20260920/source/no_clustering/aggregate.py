"""Frozen-fold selection and original-route comparisons after every A2 unit exists."""
import json
import os
from a2_common import CODE, OUT, REFERENCE, verify_source, require_slurm, checked, sha, write_json, complete, utc


def run():
    require_slurm();protocol=verify_source()
    import numpy as np
    import pandas as pd
    from scipy.stats import wilcoxon
    from statsmodels.stats.multitest import multipletests
    dest=OUT/'summary';dest.mkdir(exist_ok=True)
    folds_path=OUT/'patient_folds.csv'
    assert sha(folds_path)==protocol['patient_folds_sha256']
    folds=pd.read_csv(folds_path)
    patientfold=folds[['patient','fold']].drop_duplicates()
    assert not patientfold.patient.duplicated().any()
    measurements=[];reference=[];hashes={}
    for sample in protocol['samples']:
        for budget in protocol['budgets']:
            unit=OUT/'GBM'/sample/budget/'evaluation'
            assert checked(unit),f'Incomplete A2 unit: {sample}/{budget}'
            em=json.loads((unit/'manifest.json').read_text())
            assert em['source_bundle_sha256']==protocol['source_bundle_sha256']
            for name,expected in em['outputs'].items():assert sha(unit/name)==expected
            frame=pd.read_csv(unit/'metrics.csv')
            assert len(frame)==15 and set(frame.lambda_value)=={0,.5,1,1.5,2}
            measurements.append(frame);hashes[str(unit/'manifest.json')]=sha(unit/'manifest.json')
            old=REFERENCE/'GBM'/sample/budget/'evaluation'
            assert checked(old)
            om=json.loads((old/'manifest.json').read_text())
            assert sha(old/'metrics.csv')==om['outputs']['metrics.csv']
            oldframe=pd.read_csv(old/'metrics.csv',dtype={'cutoff':str})
            oldframe=oldframe[(oldframe.library=='CM2_glioma_other')&(oldframe.cutoff=='mean')&(oldframe.stage=='terminal090')]
            assert len(oldframe)==4
            reference.append(oldframe);hashes[str(old/'manifest.json')]=sha(old/'manifest.json')
    metrics=pd.concat(measurements,ignore_index=True);old=pd.concat(reference,ignore_index=True)
    assert metrics.primary.isin([True,False]).all() and old.primary.isin([True,False]).all(), 'Primary cohort flags must be actual parsed booleans'
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
            pat=data.groupby(['patient','fold','lambda_value'])[['macroF1_present','coverage']].mean().reset_index()
            selected=[]
            for fold in sorted(pat.fold.unique()):
                train=pat[pat.fold!=fold];test=pat[pat.fold==fold]
                scores=train.groupby('lambda_value').macroF1_present.mean().reset_index(name='training_patient_macroF1')
                assert len(scores)==5
                scores['distance_from_primary']=(scores.lambda_value-1).abs()
                rank=scores.sort_values(['training_patient_macroF1','distance_from_primary','lambda_value'],ascending=[False,True,True],kind='stable')
                chosen=rank.iloc[0]
                choices.append(dict(cohort=cohort,budget=budget,fold=int(fold),lambda_value=float(chosen.lambda_value),
                    training_patient_macroF1=float(chosen.training_patient_macroF1),n_training_patients=train.patient.nunique(),
                    n_test_patients=test.patient.nunique(),test_labels_used_for_selection=False,
                    all_training_candidate_scores=scores.to_json(orient='records')))
                selected.append(data[(data.fold==fold)&(data.lambda_value==chosen.lambda_value)].assign(selection='training_patient_selected',cohort=cohort))
            selections=[data[data.lambda_value==1].assign(selection='fixed_lambda1',cohort=cohort),pd.concat(selected)]
            for selected in selections:
                assert selected['sample'].nunique()==expected_samples and len(selected)==expected_samples
                heldout.append(selected)
                patient=selected.groupby('patient')[['macroF1_present','coverage']].mean()
                mode=selected.selection.iloc[0]
                for route,group in baseline[baseline.budget==budget].groupby('route'):
                    orig=group.groupby('patient')[['macroF1_present','coverage']].mean()
                    match=patient.join(orig,how='outer',lsuffix='_cellwise',rsuffix='_original')
                    assert not match.isna().any().any()
                    delta=(match.macroF1_present_cellwise-match.macroF1_present_original).to_numpy()
                    boot=delta[rng.integers(0,len(delta),size=(2000,len(delta)))].mean(axis=1)
                    p=1.0 if np.all(delta==0) else float(wilcoxon(delta,zero_method='pratt',alternative='two-sided',method='auto').pvalue)
                    summary.append(dict(cohort=cohort,budget=budget,selection=mode,reference_route=route,n_patients=len(delta),
                        mean_cellwise_macroF1=float(match.macroF1_present_cellwise.mean()),
                        mean_original_macroF1=float(match.macroF1_present_original.mean()),
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
    report='''# A2 GBM cell-wise seed replacement

All five predefined λ candidates are retained. Primaryλ1 and training-patient-selected
λ are separate. Terminal0.90 is the primary endpoint; marker-only and terminal0.70
are available in all_candidate_metrics.csv. No-op and unavailable learned-refinement
states are listed explicitly in execution_states.csv; they were not discarded.

The selector uses the existing patient folds, means samples within patients before
averaging training patients, and never uses the held-out patients' labels to select
λ. Exact ties choose nearest1, then smallerλ. The analysis remains retrospective.

Each original-route comparison fixes the marker library and original mean cutoff,
budget, cells and terminal endpoint. Positive delta means the cell-wise replacement
has higher macro-F1. This is a comparison of seed-construction mechanisms, not pure
clustering deletion or a proof that all no-clustering methods are inferior.

Intervals use patient bootstrap with 2,000 replicates and condition on the frozen
predictions and selected configurations. Wilcoxon tests use patient differences;
Holm correction is applied to16 predefined contrasts within each cohort. These
intervals do not include refitting or selection uncertainty. Every sample/budget
has a saved five-candidate seed/terminal figure grid on shared original UMAP axes;
that UMAP is only a display and is not an A2 model input.
'''
    (dest/'REPORT.md').write_text(report)
    files=[p for p in dest.iterdir() if p.is_file() and p.name not in ['manifest.json','COMPLETE']]
    write_json(dest/'manifest.json',dict(status='completed',source_bundle_sha256=protocol['source_bundle_sha256'],
        n_samples=len(protocol['samples']),n_budgets=2,n_candidate_conditions=1210,selection_modes=['fixed_lambda1','training_patient_selected'],
        n_sample_result_rows_by_cohort=pd.concat(heldout).groupby('cohort').size().astype(int).to_dict(),
        every_candidate_retained=True,selection_uses_training_patients_only=True,patient_folds_sha256=sha(folds_path),
        n_pairs=len(pairs),n_contrasts=len(summary),input_manifest_hashes=hashes,
        outputs={p.name:sha(p) for p in files},job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest);print('A2_AGGREGATION_COMPLETE',len(summary),flush=True)


if __name__=='__main__':run()
