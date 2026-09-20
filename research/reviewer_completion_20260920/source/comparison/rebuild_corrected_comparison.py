"""Versioned comparator/consensus replay; five methods unchanged, GNN input fixed.

The original audited evaluators run against an explicit read-only source view.
Only scDeepSort prediction/evaluation and generated summaries are new.
"""
from pathlib import Path
import sys,json,os,importlib.util
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OLD_CODE=ROOT/'handoff/paper_claim_validation_20260917'
sys.path.insert(0,str(OLD_CODE))
from common import OUT as OLD,L1,require_slurm,sha,checked,write_json,complete,utc
BASE=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/comparison'
CORRECTED=BASE/'scDeepSort_LogNormalize'
VIEW=BASE/'corrected_reference_view'
EVIDENCE=BASE/'evidence_corrected_scDeepSort'

def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def repair_effects(cohort):
    """Paired descriptive audit of changing input, keeping patient weighting."""
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    frames=[]
    for sample in cohort['sample']:
        frame=pd.read_csv(CORRECTED/sample/'metrics.csv')
        assert len(frame)==2
        frames.append(frame)
    metrics=pd.concat(frames,ignore_index=True).merge(cohort[['sample','patient','primary']],on='sample',validate='many_to_one')
    metrics.to_csv(BASE/'input_repair_sample_metrics.csv',index=False)
    measures=['macroF1_present','accuracy','coverage','mapped_coverage','unknown_rate','off_vocabulary_rate']
    rows=[]
    for name,subset in [('all121',metrics),('primary97',metrics[metrics.primary])]:
        patient=subset.groupby(['patient','condition'])[measures].mean().reset_index()
        patient['cohort']=name;rows.append(patient)
    patient=pd.concat(rows,ignore_index=True)
    patient.to_csv(BASE/'input_repair_patient_metrics.csv',index=False)
    summary=patient.groupby(['cohort','condition'])[measures].agg(['mean','std','count'])
    summary.columns=['_'.join(c) for c in summary.columns]
    summary.reset_index().to_csv(BASE/'input_repair_patient_summary.csv',index=False)
    fig,axs=plt.subplots(1,2,figsize=(10,4.4),layout='constrained')
    for ax,(name,subset) in zip(axs,patient.groupby('cohort')):
        p=subset.pivot(index='patient',columns='condition',values='macroF1_present')
        for _,row in p.iterrows():
            ax.plot([0,1],[row.historical_raw_counts,row.corrected_LogNormalize],color='#9e9e9e',alpha=.42,lw=.7)
        means=p[['historical_raw_counts','corrected_LogNormalize']].mean()
        ax.plot([0,1],means.to_numpy(),color='#0072B2',marker='o',lw=2.5)
        ax.set(xticks=[0,1],xticklabels=['Historical raw counts','Corrected LogNormalize'],ylim=(0,1),ylabel='Patient mean macro-F1',
            title=f'{name}: {len(p)} patients; mean {means.iloc[0]:.3f} → {means.iloc[1]:.3f}')
    fig.suptitle('scDeepSort input repair: same pretrained Brain weights and native threshold\nDescriptive paired audit; historical input is not a compliant comparator')
    for ext in ['png','pdf','svg']:fig.savefig(BASE/f'GBM_scDeepSort_input_repair.{ext}',dpi=250,bbox_inches='tight')
    plt.close(fig)
    return {str(p.relative_to(BASE)):sha(p) for p in BASE.iterdir() if p.name.startswith(('input_repair_','GBM_scDeepSort_input_repair.'))}

def run():
    require_slurm()
    import pandas as pd
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    pilot_proof={}
    for name in ['verification.json','verification_SN040.json','R_full_input_parity.json','R_full_input_parity_SN040.json']:
        p=CORRECTED/name;assert json.loads(p.read_text())['status']=='passed',name
        pilot_proof[name]=sha(p)
    VIEW.mkdir(parents=True,exist_ok=True)
    reuse={}
    for name in ['inputs','protocol','markers','GBM','summary','DG_fixed_partition_selection','evaluation_inputs','marker_evidence_summary']:
        target=VIEW/name
        if target.is_symlink():assert target.resolve()==(OLD/name).resolve()
        else:assert not target.exists();target.symlink_to(OLD/name,target_is_directory=True)
        reuse[str(target)]=str(OLD/name)
    (VIEW/'comparators').mkdir(exist_ok=True)
    for method in ['scType','scCATCH','SCINA','SingleR']:
        target=VIEW/'comparators'/method
        if target.is_symlink():assert target.resolve()==(OLD/'comparators'/method).resolve()
        else:assert not target.exists();target.symlink_to(OLD/'comparators'/method,target_is_directory=True)
        reuse[str(target)]=str(OLD/'comparators'/method)
    cohort=pd.read_csv(OLD/'protocol/cohort.csv');assert len(cohort)==121
    proof=[];model_hashes=None
    import evaluate_comparator
    evaluate_comparator.OUT=VIEW
    for sample in cohort['sample']:
        corrected=CORRECTED/sample;assert checked(corrected),sample
        cm=json.loads((corrected/'manifest.json').read_text())
        pm=json.loads((corrected/'prediction_manifest.json').read_text())
        im=json.loads((corrected/'input_manifest.json').read_text())
        if model_hashes is None:
            model_hashes=im['pretrained_sources']
            for path,digest in model_hashes.items():assert sha(path)==digest
        else:assert im['pretrained_sources']==model_hashes
        old=OLD/'comparators/scDeepSort'/sample
        oldm=json.loads((old/'manifest.json').read_text())
        assert sha(old/'predictions.csv.gz')==oldm['predictions_sha256']==cm['old_predictions_untouched_sha256']
        assert sha(corrected/'predictions.csv.gz')==pm['predictions_sha256']==cm['files']['predictions.csv.gz']
        assert im['normalized_exactly_once'] and not im['model_input_integers_cast']
        assert im['R_parity_n_cells']==cm['n_cells'] and im['R_parity_n_genes']==2000 and im['R_parity_max_abs_difference']<1e-6
        assert im['input_csv_sha256']==sha(corrected/'lognorm_model_input.csv')
        assert im['pretrained_sources'][next(p for p in im['pretrained_sources'] if p.endswith('human-Brain.pt'))]==oldm['model_sha256']
        dest=VIEW/'comparators/scDeepSort'/sample;dest.mkdir(parents=True,exist_ok=True)
        target=dest/'predictions.csv.gz'
        if not target.exists():target.symlink_to(corrected/'predictions.csv.gz')
        else:assert target.resolve()==(corrected/'predictions.csv.gz').resolve()
        if checked(dest):
            previous_adapter=json.loads((dest/'manifest.json').read_text())
            assert previous_adapter['predictions_sha256']==pm['predictions_sha256']
            assert previous_adapter['corrected_source_manifest_sha256']==sha(corrected/'manifest.json')
        else:
            write_json(dest/'manifest.json',dict(status='completed',method='scDeepSort',sample=sample,n_cells=cm['n_cells'],
                input='Corrected Seurat LogNormalize before gene intersection; no int cast',predictions_sha256=pm['predictions_sha256'],
                model_sha256=oldm['model_sha256'],corrected_source_manifest_sha256=sha(corrected/'manifest.json'),
                official_predictor_sha256=pm['official_predictor_sha256'],new_training=False,job=os.environ['SLURM_JOB_ID']))
            complete(dest)
        evaluate_comparator.run('scDeepSort',sample)
        proof.append(dict(sample=sample,n_cells=cm['n_cells'],changed_calls=cm['changed_native_predictions_vs_raw_counts'],
            R_parity_max_abs_difference=im['R_parity_max_abs_difference'],corrected_source_manifest_sha256=sha(corrected/'manifest.json'),
            previous_predictions_unchanged=True,checkpoint_unchanged=True))
    pd.DataFrame(proof).to_csv(BASE/'corrected_121_source_verification.csv',index=False)
    repair_files=repair_effects(cohort)
    write_json(VIEW/'reuse_sources.json',reuse)
    import aggregate_comparators
    aggregate_comparators.OUT=VIEW
    aggregate_comparators.run()
    comparison=json.loads((VIEW/'comparison_summary/manifest.json').read_text())
    comparison['scDeepSort_input_contract']='Corrected Seurat LogNormalize, prior to model-gene intersection; 121 samples verified'
    comparison['scDeepSort_endpoint']='Unchanged published Brain checkpoint and native unsure_rate=2; no retraining'
    comparison['supersedes']='Historical raw-count scDeepSort comparison, which did not meet the official input contract'
    comparison['other_five_methods']='Frozen previously fitted predictions, labels, patients, marker-selection opportunities and scoring unchanged'
    write_json(VIEW/'comparison_summary/manifest.json',comparison)
    complete(VIEW/'comparison_summary')
    # Confirm the other five methods' held-out selections/metrics are unchanged.
    previous=pd.read_csv(OLD/'comparison_summary/patient_heldout_results.csv',dtype={'cutoff':str})
    current=pd.read_csv(VIEW/'comparison_summary/patient_heldout_results.csv',dtype={'cutoff':str})
    keys=['cohort','patient','method']
    a=previous[previous.method!='scDeepSort'].sort_values(keys).reset_index(drop=True)
    b=current[current.method!='scDeepSort'].sort_values(keys).reset_index(drop=True)
    pd.testing.assert_frame_equal(a,b,check_exact=False,rtol=0,atol=1e-14)
    import prediction_helpers
    prediction_helpers.OUT=VIEW
    import comparison_diagnostics
    comparison_diagnostics.OUT=VIEW
    comparison_diagnostics.run()
    evidence=module('corrected_saved_evidence',ROOT/'handoff/reviewer_completion_20260920/evidence/gbm_evidence.py')
    evidence.OUT=VIEW;evidence.DEST=EVIDENCE
    evidence.run()
    verifier=module('corrected_evidence_verifier',ROOT/'handoff/reviewer_completion_20260920/evidence/verify_evidence.py')
    verifier.OUT=EVIDENCE;verifier.run()
    # Source supplement is preserved; the new endpoint contract is explicit here.
    note='''# Corrected-input scDeepSort comparison version

All 121 GBM predictions use the official pretrained Brain checkpoint and native
threshold, with Seurat-compatible LogNormalize input. The prior raw-count results
remain preserved. Every sample's 2000-HVG expression was checked against native R;
three pilots additionally check all model-input genes and cells directly against
the original R RNA assay. No model weights were retrained. The Brain vocabulary
still lacks malignant labels; unsupported cells remain in the all-cell denominator.

The other five tools, patient folds, truth labels, marker selection opportunities
and metric definitions are unchanged. Their held-out table is explicitly checked
against the previous version. Only corrected scDeepSort changes; pooled summaries,
patient-paired comparisons, Holm correction, all-tool consensus and visualizations
are rebuilt. The previous D/F evidence is historical saved-output evidence.

The corrected_reference_view directory is a declared local source view: symlinks
point to unchanged inputs and five existing tools. It is not a second raw dataset.
reuse_sources.json describes every reused root. The generated comparison_summary
and evidence_corrected_scDeepSort directories contain new derived artifacts.
'''
    (BASE/'CORRECTED_COMPARISON.md').write_text(note)
    write_json(BASE/'corrected_comparison_manifest.json',dict(status='completed',n_samples=121,
        old_other_five_heldout_rows_unchanged=True,corrected_source_verification_sha256=sha(BASE/'corrected_121_source_verification.csv'),
        input_repair_audit_files=repair_files,
        three_pilot_independent_full_R_input_proofs=pilot_proof,current_pretrained_files=model_hashes,
        official_source_version='scDeepSort 1.0',official_setup_sha256=sha('/fs/scratch/PCON0080/yimin/tools/deepsort-1.0/setup.py'),
        comparison_manifest_sha256=sha(VIEW/'comparison_summary/manifest.json'),
        evidence_manifest_sha256=sha(EVIDENCE/'manifest.json'),evidence_verification_sha256=sha(EVIDENCE/'verification.json'),
        source_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    import verify_corrected_comparison
    verify_corrected_comparison.run()
    state=json.loads((CORRECTED/'status.json').read_text());state.update(status='completed',verified=True,updated_at=utc(),
        completed=['121/121 corrected-input predictions','Patient-paired comparator metrics and corrected D/F consensus verified'],remaining=[])
    state['jobs']=list(dict.fromkeys(state['jobs']+[os.environ['SLURM_JOB_ID']]))
    state['evidence'].append(str(BASE/'corrected_comparison_manifest.json'))
    state['evidence'].append(str(BASE/'corrected_comparison_verification.json'))
    write_json(CORRECTED/'status.json',state)
    audit=json.loads((BASE/'status.json').read_text())
    audit['updated_at']=utc()
    audit['jobs']=list(dict.fromkeys(audit['jobs']+state['jobs']))
    audit['completed']=list(dict.fromkeys(audit['completed']+['GBM 121-sample corrected scDeepSort comparison and D/F independent verification']))
    audit['remaining']=[v for v in audit['remaining'] if not v.startswith('GBM corrected-input scDeepSort')]
    audit['evidence']=list(dict.fromkeys(audit['evidence']+[str(BASE/'corrected_comparison_manifest.json'),str(BASE/'corrected_comparison_verification.json')]))
    write_json(BASE/'status.json',audit)
    print('CORRECTED_COHORT_COMPARISON_COMPLETE',flush=True)

if __name__=='__main__':run()
