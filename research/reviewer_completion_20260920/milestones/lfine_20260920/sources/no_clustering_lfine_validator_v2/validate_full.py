"""Independent saved-label Lfine counts, patient selection, and paired inference.

The evaluator/provider and aggregator are never called. Counts use contingency
tables; selection uses explicit patient dictionaries; Holm uses its definition.
No fits, expression matrices, notebook writes, or endpoint changes.
"""
from pathlib import Path
from collections import Counter, defaultdict
import json, os, sys, math
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/no_clustering_lfine_v1'
REF=ROOT/'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917'
COMPACT=ROOT/'results/hvg_ptc_20260916_v1/lfine_compact_20260920'
sys.path.insert(0,str(OUT/'source_v2'))
import common as c

def close(a,b):
    assert np.isclose(a,b,rtol=0,atol=1e-12,equal_nan=True),(a,b)

def key(row):
    return (row.sample,row.budget,row.family,row.route,None if pd.isna(row.lambda_value) else float(row.lambda_value),row.stage)

def verify_sample(sample, protocol, inputs, compact, prefixes, mapper):
    unit=OUT/'samples'/sample
    assert c.checked(unit),sample
    manifest=c.read(unit/'manifest.json')
    assert manifest['status']=='completed'
    assert manifest['protocol_sha256']==c.sha(OUT/'protocol_v2.json')
    assert manifest['source_manifest_sha256']==c.sha(OUT/'source_v2/SOURCE_MANIFEST.json')
    assert manifest['input_record']==inputs['samples'][sample]
    c.verify_files(manifest['input_record']['files'])
    for name,digest in manifest['outputs'].items():assert c.sha(unit/name)==digest
    frame=pd.read_csv(unit/'metrics.csv');counts=pd.read_csv(unit/'per_class_counts.csv.gz')
    anchors=pd.read_csv(unit/'compact_anchor_parity.csv');c.bool_column(frame)
    for col in ['training_executed','terminal_valid','source_native_labels_preserved']:
        assert pd.api.types.is_bool_dtype(frame[col]) and frame[col].notna().all()
    assert len(frame)==54 and len(anchors)==8 and anchors.all_fields_match.eq(True).all()
    expected={(sample,b,'cellwise_seed','cellwise_seed',float(a['lambda']),s)
              for b in protocol['budgets'] for a in protocol['arms'] for s in c.STAGES}
    expected|={(sample,b,'original_cluster_DEG',r,None,s) for b in protocol['budgets'] for r in c.ROUTES for s in c.STAGES}
    assert set(map(key,frame.itertuples()))==expected
    assert frame.status.eq('completed').all() and frame.terminal_valid.all()
    truth=pd.read_csv(REF/'evaluation_inputs'/sample/'truth.csv.gz',dtype=str,keep_default_na=False)
    assert truth.cell_id.is_unique and truth.Patient.nunique()==1
    assert frame.patient.eq(truth.Patient.iloc[0]).all()
    labels=truth.lfine_original.tolist();support=Counter(labels)
    # Independently reconstruct the existing reference endpoint from raw columns.
    bad={'nan','','None','NA'}
    recomposed=[f'Malignant_{m}' if m not in bad else l3 if l3 not in bad else l1
                for l1,l3,m in zip(truth.L1,truth.L3,truth.MalState)]
    assert recomposed==labels
    classes={label for label,n in support.items() if n>=20 and label not in {'Other','nan'}}
    assert classes,'An empty endpoint scope must not be silently averaged'
    targets={parent:{label for label in support if any(label.startswith(p) for p in prefs)} for parent,prefs in prefixes.items()}
    reachable=set().union(*(targets.get(parent,set()) for parent in set(mapper.values())))
    abstain={'Unknown','Undecided','Noise','nan','None',''}
    count_index={}
    for row in counts.itertuples():
        identifier=key(row)+(row.lfine_class,)
        assert identifier not in count_index
        count_index[identifier]=(row.TP,row.FP,row.FN,row.support)
    cached={};terminals={};checked_count=0;recomputed=[]
    for row in frame.itertuples():
        td=ROOT/row.original_prediction_source
        if td not in terminals:
            tm=c.read(td/'terminal_manifest.json');terminals[td]=tm
            assert c.checked(td,'terminal_manifest.json','TERMINAL_COMPLETE')
            assert c.sha(td/'predictions.csv.gz')==tm['predictions_sha256']
            assert c.sha(td/'terminal.npz')==tm['terminal_sha256']
            prediction=pd.read_csv(td/'predictions.csv.gz',dtype=str,keep_default_na=False)
            assert prediction.cell_id.tolist()==truth.cell_id.tolist()
            with np.load(td/'terminal.npz',allow_pickle=False) as arrays:
                for column in ['initial','final090','final070']:
                    assert arrays[column].tolist()==prediction[column].tolist()
            known=prediction.initial.ne('Undecided')
            for column in ['final090','final070']:
                assert prediction.loc[known,column].tolist()==prediction.loc[known,'initial'].tolist()
            cached[td/'predictions.csv.gz']=prediction
        tm=terminals[td]
        assert row.training_executed==tm['training_executed'] and row.dl_status==tm['dl_status']
        assert row.terminal_manifest_sha256==c.sha(td/'terminal_manifest.json')
        path=td.parent.parent/'initial_calls.csv.gz' if row.stage=='initial' else td/'predictions.csv.gz'
        column=row.arm_id if row.stage=='initial' else c.STAGES[row.stage]
        if path not in cached:cached[path]=pd.read_csv(path,dtype=str,keep_default_na=False)
        calls=cached[path];assert calls.cell_id.tolist()==truth.cell_id.tolist()
        semantic=[mapper.get(label,'Unknown' if label in abstain else 'UNMAPPABLE') for label in calls[column]]
        table=Counter(zip(labels,semantic));calls_count=Counter(semantic)
        correct=sum(n for (gold,pred),n in table.items() if gold in targets.get(pred,set()))
        fs=[];hits=0
        for label in sorted(classes):
            tp=sum(n for (gold,pred),n in table.items() if gold==label and gold in targets.get(pred,set()))
            fn=support[label]-tp
            fp=sum(n for (gold,pred),n in table.items() if gold!=label and gold not in targets.get(pred,set()) and label in targets.get(pred,set()))
            assert count_index[key(row)+(label,)]==(tp,fp,fn,support[label])
            checked_count+=1;hits+=tp>0
            fs.append(2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else 0.)
        n=len(labels);n_abstain=sum(v for label,v in calls_count.items() if label in abstain)
        n_off=sum(v for label,v in calls_count.items() if label not in abstain and label not in prefixes)
        n_called=n-n_abstain-n_off
        metrics=dict(lfine_macroF1=sum(fs)/len(fs),coverage=1-n_abstain/n,
            legacy_called_coverage=n_called/n,abstain_rate=n_abstain/n,offvocab_rate=n_off/n,
            acc_on_called=correct/n_called if n_called else np.nan,
            n_distinct_calls=sum(label in prefixes and label not in abstain for label in calls_count),
            n_classes_hit=hits,lfine_n_classes=len(classes),
            lfine_scored_class_cell_fraction=sum(support[label] for label in classes)/n,
            marker_vocab_oracle_upper_bound=sum(label in reachable for label in classes)/len(classes),
            n_reference_lfine_disagreements=0)
        assert row.n_cells==n
        for field,value in metrics.items():close(value,getattr(row,field))
        if row.family=='original_cluster_DEG' and row.stage=='terminal090':
            anchor=compact[(compact['sample']==sample)&(compact.budget==row.budget)&(compact.route==row.route)&(compact.library==c.FIXED)]
            assert len(anchor)==1
            for field in c.FIELDS:close(metrics[field],anchor.iloc[0][field])
            assert row.predictions_sha256==anchor.iloc[0].predictions_sha256
        recomputed.append(dict(sample=sample,budget=row.budget,family=row.family,route=row.route,
            lambda_value=None if pd.isna(row.lambda_value) else row.lambda_value,stage=row.stage,**metrics))
    assert checked_count==len(counts)
    return frame,dict(sample=sample,n_cells=len(labels),n_metric_rows=54,n_class_count_rows=checked_count,
        n_compact_anchors=8,manifest_sha256=c.sha(unit/'manifest.json')),recomputed

def verify_selection_and_statistics(protocol,candidates,reference):
    summary=OUT/'summary';manifest=c.read(summary/'manifest.json');assert c.checked(summary)
    assert manifest['source_manifest_sha256']==c.sha(OUT/'source_v2/SOURCE_MANIFEST.json')
    assert manifest['protocol_sha256']==c.sha(OUT/'protocol_v2.json')
    assert manifest['n_samples']==121 and manifest['n_candidate_conditions']==1210 and manifest['n_contrasts']==32
    assert manifest['n_sample_result_rows_by_cohort']=={'all':484,'primary':388}
    for name,digest in manifest['outputs'].items():assert c.sha(summary/name)==digest
    for path,digest in manifest['input_manifest_hashes'].items():assert c.sha(path)==digest
    assert len(manifest['input_manifest_hashes'])==121
    assert pd.read_csv(summary/'unavailable_conditions.csv').empty
    for name,expected in [('all_candidate_metrics.csv',candidates),('all_original_route_metrics.csv',reference)]:
        actual=pd.read_csv(summary/name);c.bool_column(actual)
        pd.testing.assert_frame_equal(actual,expected.reset_index(drop=True),check_exact=False,rtol=0,atol=1e-12)
    folds=pd.read_csv(c.OLD/'patient_folds.csv')
    assert c.sha(c.OLD/'patient_folds.csv')==protocol['patient_folds_sha256']
    assert set(folds['sample'])==set(protocol['samples']) and not folds['sample'].duplicated().any()
    pf=folds[['patient','fold']].drop_duplicates();assert pf.patient.is_unique
    fold_by_sample=folds.set_index('sample').fold.to_dict()
    choices=pd.read_csv(summary/'training_patient_lambda_choices.csv');assert len(choices)==20
    assert pd.api.types.is_bool_dtype(choices.test_labels_used_for_selection)
    assert choices.test_labels_used_for_selection.eq(False).all()
    saved_selected=pd.read_csv(summary/'selected_and_fixed_sample_results.csv');c.bool_column(saved_selected)
    paired=pd.read_csv(summary/'patient_paired_results.csv')
    comparisons=pd.read_csv(summary/'patient_paired_summary.csv')
    assert len(comparisons)==32 and not comparisons.duplicated(['cohort','budget','selection','reference_route']).any()
    rng=np.random.default_rng(42);results=[];accounting=[];verified_selected=0;verified_pairs=0
    for cohort in ['primary','all']:
        terminal=candidates[candidates.stage.eq('terminal090')]
        orig=reference[reference.stage.eq('terminal090')]
        if cohort=='primary':terminal=terminal[terminal.primary];orig=orig[orig.primary]
        expected_n=97 if cohort=='primary' else 121
        expected_patients=55 if cohort=='primary' else 59
        assert terminal['sample'].nunique()==orig['sample'].nunique()==expected_n
        for budget in protocol['budgets']:
            data=terminal[terminal.budget.eq(budget)]
            patient_values=defaultdict(list);rows_by_sample={}
            for row in data.itertuples():
                patient_values[(row.patient,fold_by_sample[row.sample],row.lambda_value)].append(row.lfine_macroF1)
                rows_by_sample[(row.sample,row.lambda_value)]=row
            patient_means={k:sum(values)/len(values) for k,values in patient_values.items()}
            picked={}
            for fold in sorted(set(fold_by_sample.values())):
                scores={lam:sum(v for (pid,f,l),v in patient_means.items() if f!=fold and l==lam)/sum(f!=fold and l==lam for pid,f,l in patient_means)
                        for lam in [0.,.5,1.,1.5,2.]}
                chosen=min(scores,key=lambda lam:(-scores[lam],abs(lam-1),lam));picked[fold]=chosen
                saved=choices[(choices.cohort==cohort)&(choices.budget==budget)&(choices.fold==fold)]
                assert len(saved)==1 and saved.iloc[0].lambda_value==chosen
                close(saved.iloc[0].training_patient_lfine_macroF1,scores[chosen])
                training={pid for pid,f,lam in patient_means if f!=fold}
                test={pid for pid,f,lam in patient_means if f==fold}
                assert training.isdisjoint(test)
                assert saved.iloc[0].n_training_patients==len(training) and saved.iloc[0].n_test_patients==len(test)
                fullscores=json.loads(saved.iloc[0].all_training_candidate_scores)
                assert {v['lambda_value'] for v in fullscores}==set(scores)
                # pandas DataFrame.to_json defaults to ten decimal places. This
                # diagnostic field is rounded; the selected mean and all CSV
                # values above/below retain the strict 1e-12 acceptance bound.
                for value in fullscores:
                    assert abs(value['training_patient_lfine_macroF1']-scores[value['lambda_value']])<=5.1e-11
            for mode in ['fixed_lambda1','training_patient_selected']:
                chosen_rows={sample:rows_by_sample[(sample,1 if mode=='fixed_lambda1' else picked[fold_by_sample[sample]])]
                             for sample in sorted(set(data['sample']))}
                saved=saved_selected[(saved_selected.cohort==cohort)&(saved_selected.budget==budget)&(saved_selected.selection==mode)]
                assert len(saved)==expected_n and set(saved['sample'])==set(chosen_rows)
                for row in saved.itertuples():
                    expected=chosen_rows[row.sample]
                    assert row.lambda_value==expected.lambda_value and row.fold==fold_by_sample[row.sample]
                    assert row.training_executed==expected.training_executed and row.dl_status==expected.dl_status
                    for field in c.FIELDS:close(getattr(row,field),getattr(expected,field))
                verified_selected+=len(saved)
                states=Counter((row.lambda_value,row.dl_status,bool(row.training_executed)) for row in chosen_rows.values())
                for (lam,state,trained),n in sorted(states.items()):
                    accounting.append(dict(cohort=cohort,budget=budget,selection=mode,lambda_value=lam,dl_status=state,training_executed=trained,n_samples=n))
                pat=defaultdict(list)
                for row in chosen_rows.values():pat[row.patient].append((row.lfine_macroF1,row.coverage))
                avg={pid:np.mean(values,axis=0) for pid,values in pat.items()}
                for route in sorted(c.ROUTES):
                    base=orig[(orig.budget==budget)&(orig.route==route)];by_pat=defaultdict(list)
                    for row in base.itertuples():by_pat[row.patient].append((row.lfine_macroF1,row.coverage))
                    old={pid:np.mean(values,axis=0) for pid,values in by_pat.items()}
                    ids=sorted(avg);assert set(ids)==set(old) and len(ids)==expected_patients
                    delta=np.array([avg[pid][0]-old[pid][0] for pid in ids])
                    boots=np.array([delta[draw].mean() for draw in rng.integers(len(ids),size=(2000,len(ids)))])
                    p=1. if np.all(delta==0) else wilcoxon(delta,zero_method='pratt',alternative='two-sided',method='auto').pvalue
                    want=dict(n_patients=len(ids),mean_cellwise_lfine_macroF1=np.mean([avg[pid][0] for pid in ids]),
                        mean_original_lfine_macroF1=np.mean([old[pid][0] for pid in ids]),mean_delta_cellwise_minus_original=delta.mean(),
                        CI025=np.quantile(boots,.025),CI975=np.quantile(boots,.975),p_value=p,
                        mean_coverage_delta=np.mean([avg[pid][1]-old[pid][1] for pid in ids]))
                    subset=comparisons[(comparisons.cohort==cohort)&(comparisons.budget==budget)&(comparisons.selection==mode)&(comparisons.reference_route==route)]
                    assert len(subset)==1
                    for field,value in want.items():close(value,subset.iloc[0][field])
                    paired_subset=paired[(paired.cohort==cohort)&(paired.budget==budget)&(paired.selection==mode)&(paired.reference_route==route)]
                    assert len(paired_subset)==len(ids) and set(paired_subset.patient)==set(ids)
                    for row in paired_subset.itertuples():
                        for field,value in [('lfine_macroF1_cellwise',avg[row.patient][0]),('coverage_cellwise',avg[row.patient][1]),
                                            ('lfine_macroF1_original',old[row.patient][0]),('coverage_original',old[row.patient][1])]:close(getattr(row,field),value)
                    verified_pairs+=len(paired_subset)
                    results.append(dict(cohort=cohort,budget=budget,selection=mode,reference_route=route,**want))
        cohort_results=[r for r in results if r['cohort']==cohort]
        assert len(cohort_results)==16
        ordered=sorted(enumerate(cohort_results),key=lambda t:t[1]['p_value']);current=0.
        for rank,(index,row) in enumerate(ordered):
            current=max(current,min(1.,(16-rank)*row['p_value']));row['p_holm']=current
            saved=comparisons[(comparisons.cohort==cohort)&(comparisons.budget==row['budget'])&(comparisons.selection==row['selection'])&(comparisons.reference_route==row['reference_route'])]
            close(current,saved.iloc[0].p_holm)
    assert verified_selected==len(saved_selected)==872 and verified_pairs==len(paired)==1824
    return results,accounting

def main():
    assert os.environ.get('SLURM_JOB_ID')
    source_manifest=Path(__file__).parent/'SOURCE_MANIFEST.json'
    sm=c.read(source_manifest)
    for name,digest in sm['files'].items():assert c.sha(Path(__file__).parent/name)==digest
    for path,digest in sm['pinned_evidence'].items():assert c.sha(ROOT/path)==digest
    protocol=c.verify();approval=c.read(OUT/'APPROVED.json')
    assert approval['allow_aggregate'] is True and set(approval['allowed_samples'])==set(protocol['samples'])
    inputs=c.read(OUT/'inputs.json');compact=pd.read_csv(COMPACT/'metrics_hvg24.csv.gz')
    prefixes=c.read(COMPACT/'lfine_target_prefixes.json')
    mapping=pd.read_csv(REF/'markers/panel_L1_mapping.csv',dtype=str,keep_default_na=False)
    mapping=mapping[mapping.library.eq(c.FIXED)];assert mapping.panel.is_unique
    mapper=dict(zip(mapping.panel,mapping.L1));frames=[];receipts=[];metrics=[]
    for sample in protocol['samples']:
        frame,receipt,rows=verify_sample(sample,protocol,inputs,compact,prefixes,mapper)
        frames.append(frame);receipts.append(receipt);metrics.extend(rows)
        print('A2_LFINE_INDEPENDENT_SAMPLE',sample,len(receipts),flush=True)
    allrows=pd.concat(frames,ignore_index=True)
    candidates=allrows[allrows.family.eq('cellwise_seed')].copy()
    reference=allrows[allrows.family.eq('original_cluster_DEG')].copy()
    assert len(candidates)==3630 and len(reference)==2904
    results,accounting=verify_selection_and_statistics(protocol,candidates,reference)
    dest=OUT/'full_validation';dest.mkdir(exist_ok=True)
    assert not (dest/'validation.json').exists(),'Preserve completed independent proof'
    pd.DataFrame(metrics).to_csv(dest/'independent_metric_rows.csv.gz',index=False)
    pd.DataFrame(receipts).to_csv(dest/'sample_roster.csv',index=False)
    pd.DataFrame(results).to_csv(dest/'independent_paired_summary.csv',index=False)
    pd.DataFrame(accounting).to_csv(dest/'selected_training_accounting.csv',index=False)
    c.verify()
    report=dict(status='passed',scope='Independent A2 Lfine saved-prediction replay and conditional patient inference; HVG2000/5000 only',
        n_samples=121,n_cells=sum(r['n_cells'] for r in receipts),n_metric_rows=6534,n_candidate_conditions=1210,
        n_perclass_count_rows=sum(r['n_class_count_rows'] for r in receipts),n_compact_anchors=968,
        n_training_fold_choices=20,n_selected_sample_rows=872,n_patient_pairs=1824,n_contrasts=32,
        all12_metric_fields_independently_recomputed=True,all_terminal_csv_npz_arrays_and_cell_orders_match=True,
        all_input_hashes_verified=True,no_invalid_or_silently_excluded_conditions=True,reference_lfine_recomposition_exact=True,
        selection_training_patients_only=True,selection_tie_rule_verified=True,
        patient_bootstrap_2000_seed42_recomputed=True,wilcoxon_pratt_recomputed=True,holm_definition_recomputed=True,
        diagnostic_training_scores_json_absolute_tolerance=5.1e-11,
        diagnostic_precision_reason='pandas to_json default double_precision=10; all source CSV/selected mean/statistics checks remain 1e-12',
        no_model_fitting=True,no_writes_to_original_L1_archive=True,whole_work_package_A_complete=False,
        source_manifest_sha256=c.sha(OUT/'source_v2/SOURCE_MANIFEST.json'),protocol_sha256=c.sha(OUT/'protocol_v2.json'),
        summary_manifest_sha256=c.sha(OUT/'summary/manifest.json'),validation_source_manifest_sha256=c.sha(source_manifest),
        approval_sha256=c.sha(OUT/'APPROVED.json'),sample_manifests={r['sample']:r['manifest_sha256'] for r in receipts},
        outputs={p.name:c.sha(p) for p in dest.iterdir() if p.is_file()},
        job_id=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),completed_at=c.utc())
    c.write(dest/'validation.json',report)
    (dest/'COMPLETE').write_text(c.sha(dest/'validation.json')+'\n')
    print('A2_LFINE_FULL_VALIDATION_PASSED',json.dumps({k:v for k,v in report.items() if k not in ['sample_manifests','outputs']}),flush=True)

if __name__=='__main__':main()
