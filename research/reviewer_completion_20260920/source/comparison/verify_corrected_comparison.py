"""Independent accounting from corrected native calls, with all cells retained."""
from pathlib import Path
import os,sys,json
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
sys.path.insert(0,str(ROOT/'handoff/paper_claim_validation_20260917'))
from common import OUT as OLD,L1,require_slurm,sha,write_json,utc
BASE=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/comparison'
CORRECTED=BASE/'scDeepSort_LogNormalize'
VIEW=BASE/'corrected_reference_view'

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    from scipy.stats import wilcoxon
    from build_markers import semantic
    top=json.loads((BASE/'corrected_comparison_manifest.json').read_text())
    assert top['status']=='completed' and top['n_samples']==121
    manifest=json.loads((VIEW/'comparison_summary/manifest.json').read_text())
    for name,digest in manifest['files'].items():assert sha(VIEW/'comparison_summary'/name)==digest,name
    for name,digest in top['input_repair_audit_files'].items():assert sha(BASE/name)==digest,name
    assert sha(VIEW/'comparison_summary/manifest.json')==top['comparison_manifest_sha256']
    assert sha(BASE/'evidence_corrected_scDeepSort/verification.json')==top['evidence_verification_sha256']
    assert json.loads((BASE/'evidence_corrected_scDeepSort/verification.json').read_text())['status']=='passed'
    cohort=pd.read_csv(OLD/'protocol/cohort.csv');assert len(cohort)==121 and cohort.patient.nunique()==59
    fresh=[]
    measures=['macroF1_present','macroF1_fixed11','accuracy','coverage','unknown_rate','mapped_coverage','off_vocabulary_rate']
    for row in cohort.itertuples():
        sample=row.sample;dest=CORRECTED/sample
        predictions=pd.read_csv(dest/'predictions.csv.gz',dtype=str,keep_default_na=False)
        truth=pd.read_csv(OLD/'evaluation_inputs'/sample/'truth.csv.gz',dtype=str,keep_default_na=False)
        assert predictions.cell_id.is_unique and predictions.cell_id.tolist()==truth.cell_id.tolist()
        m=json.loads((dest/'prediction_manifest.json').read_text());assert sha(dest/'predictions.csv.gz')==m['predictions_sha256']
        mapped=predictions.default.map(semantic).to_numpy();y=truth.L1.to_numpy();assert len(y)==row.n_cells
        assert np.isin(y,L1).all()
        counts=pd.crosstab(pd.Series(y,name='author'),pd.Series(mapped,name='prediction')).reindex(index=L1,fill_value=0)
        support=counts.sum(axis=1)
        # Count false positives from every author class; native labels outside L1
        # still consume cells and false negatives in the complete denominator.
        called=pd.Series(mapped).value_counts().reindex(L1,fill_value=0)
        tp=np.array([counts.at[lab,lab] if lab in counts else 0 for lab in L1])
        denom=support.to_numpy()+called.to_numpy()
        f1=np.divide(2*tp,denom,out=np.zeros(len(L1),dtype=float),where=denom>0)
        native_unknown=predictions.default.isin(['unknown','Unknown','Undecided','','NA']).to_numpy()
        supported=np.isin(mapped,L1)
        values=dict(sample=sample,patient=row.patient,primary=row.primary,n_cells=len(y),
            macroF1_present=float(f1[support.to_numpy()>0].mean()),macroF1_fixed11=float(f1.mean()),
            accuracy=float(tp.sum()/len(y)),coverage=float((~native_unknown).mean()),unknown_rate=float(native_unknown.mean()),
            mapped_coverage=float(supported.mean()),off_vocabulary_rate=float((~native_unknown&~supported).mean()))
        saved=pd.read_csv(VIEW/'comparators/scDeepSort'/sample/'evaluation/metrics.csv');assert len(saved)==1
        for key in measures:np.testing.assert_allclose(values[key],saved.iloc[0][key],atol=1e-14,rtol=0,err_msg=f'{sample}/{key}')
        fresh.append(values)
    fresh=pd.DataFrame(fresh);assert fresh.n_cells.sum()==429305
    assert len(fresh[fresh.primary])==97 and fresh[fresh.primary].n_cells.sum()==370751
    current=pd.read_csv(VIEW/'comparison_summary/patient_heldout_results.csv')
    summary=pd.read_csv(VIEW/'comparison_summary/patient_heldout_summary.csv')
    contrasts=pd.read_csv(VIEW/'comparison_summary/paired_patient_comparisons.csv')
    old=pd.read_csv(OLD/'comparison_summary/patient_heldout_results.csv')
    keys=['cohort','patient','method']
    pd.testing.assert_frame_equal(old[old.method!='scDeepSort'].sort_values(keys).reset_index(drop=True),
        current[current.method!='scDeepSort'].sort_values(keys).reset_index(drop=True),check_exact=False,atol=1e-14,rtol=0)
    for name,subset,n in [('all121',fresh,59),('primary97',fresh[fresh.primary],55)]:
        patient=subset.groupby('patient')[measures].mean().sort_index();assert len(patient)==n
        saved=current[(current.cohort==name)&(current.method=='scDeepSort')].set_index('patient').sort_index()
        assert patient.index.tolist()==saved.index.tolist()
        for key in measures:
            np.testing.assert_allclose(patient[key],saved[key],atol=1e-14,rtol=0,err_msg=f'{name}/{key}')
            record=summary[(summary.cohort==name)&(summary.method=='scDeepSort')];assert len(record)==1
            np.testing.assert_allclose(patient[key].mean(),record.iloc[0][key+'_mean'],atol=1e-14,rtol=0)
            assert record.iloc[0][key+'_count']==n
        all_methods=current[current.cohort==name].pivot(index='patient',columns='method',values='macroF1_present')
        pvalues={}
        for method in all_methods.columns.drop('DG-scRNA'):
            delta=(all_methods[method]-all_methods['DG-scRNA']).to_numpy()
            r=contrasts[(contrasts.cohort==name)&(contrasts.method==method)].iloc[0]
            np.testing.assert_allclose(delta.mean(),r.mean_delta,atol=1e-14,rtol=0)
            p=wilcoxon(delta).pvalue if np.any(np.abs(delta)>1e-14) else 1.
            np.testing.assert_allclose(p,r.p,atol=1e-14,rtol=0);pvalues[method]=p
        ordered=sorted(pvalues,key=pvalues.get);floor=0.
        for rank,method in enumerate(ordered):
            floor=max(floor,min(1.,pvalues[method]*(len(ordered)-rank)))
            r=contrasts[(contrasts.cohort==name)&(contrasts.method==method)].iloc[0]
            np.testing.assert_allclose(floor,r.p_Holm,atol=1e-14,rtol=0)
    for name in ['NL022_display_gene_cell_expression.csv.gz','NL022_frozen_endpoint_labels.csv.gz',
                 'NL022_all_selected_marker_expression_by_endpoint.csv.gz','NL022_display_gene_selection.csv']:
        original=pd.read_csv(BASE.parent/'evidence'/name);new=pd.read_csv(BASE/'evidence_corrected_scDeepSort'/name)
        pd.testing.assert_frame_equal(original,new,check_exact=False,atol=1e-14,rtol=0)
    proof=dict(status='passed',n_samples=121,n_cells=429305,primary_samples=97,primary_cells=370751,
        patient_denominators={'all121':59,'primary97':55},
        independent_native_call_confusion_accounting=True,all_cell_and_unknown_denominators=True,
        corrected_patient_means_and_method_summary_recomputed=True,paired_means_wilcoxon_Holm_rechecked=True,
        five_other_methods_unchanged=True,NL022_expression_and_endpoint_tables_unchanged=True,
        comparison_manifest_sha256=sha(VIEW/'comparison_summary/manifest.json'),
        corrected_comparison_manifest_sha256=sha(BASE/'corrected_comparison_manifest.json'),
        source_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],verified_at=utc())
    write_json(BASE/'corrected_comparison_verification.json',proof)
    print('INDEPENDENT_CORRECTED_COMPARISON_VERIFIED',json.dumps(proof),flush=True)

if __name__=='__main__':run()
