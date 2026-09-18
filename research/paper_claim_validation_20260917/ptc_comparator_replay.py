"""Re-evaluate saved PTC competitors without new fitting or outcome-based mapping."""
import json
import os
from pathlib import Path
from common import ROOT,OUT,sha,checked,write_json,complete,utc
from ptc_followup_common import PTC,REFERENCE,PAPER,require_ptc,read_reference
from ptc_label_rules import broad_lineage
from ptc_evaluate import binary_metrics
from ptc_aggregate import paired

# SignacX2.2.5's published GenerateLabels hierarchy and model vocabulary.
# The finer CellStates output is saved alongside the coarser CellTypes output.
# These identities are defined by the method, without inspecting TCR outcomes.
SIGNAC_T={'T.CD4.naive','T.CD4.memory','T.regs','T.CD8.naive','T.CD8.memory',
          'T.CD8.cm','T.CD8.em','T.gd'}
SIGNAC_SOURCE='https://github.com/cran/SignacX/blob/2.2.5/R/helper_functions.R'
UNKNOWN={'','unknown','unclassified','undecided','no_annotation','none','nan','na','unassigned','unlabeled'}

def run():
    require_ptc()
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    dest=OUT/'PTC_comparator_replay';dest.mkdir(exist_ok=True)
    ref=read_reference();source=ROOT/'results/competitors/signacx/signacx_all.csv'
    data=pd.read_csv(source,keep_default_na=False)
    assert data.Sample_ID.is_unique and len(data)==107545
    assert set(ref.index)<=set(data.Sample_ID)
    parts=sorted(source.parent.glob('signacx_*.csv'))
    parts=[p for p in parts if p.name!='signacx_all.csv']
    assert len(parts)==8
    assembled=pd.concat([pd.read_csv(p,keep_default_na=False) for p in parts],ignore_index=True)
    assert assembled.Sample_ID.is_unique
    pd.testing.assert_frame_equal(data.set_index('Sample_ID').sort_index(),
        assembled.set_index('Sample_ID').sort_index(),check_like=True)
    data=data.set_index('Sample_ID');joined=data.loc[ref.index]
    assert np.array_equal(joined['sample'].to_numpy(),ref['sample'].to_numpy())
    assert set(joined['sample'])==set(ref['sample'])
    assert set(joined.CellStates[joined.CellStates.str.startswith('T.')])<=SIGNAC_T
    old=REFERENCE.parent/'archived_comparators'
    assert checked(old)
    ontology=pd.read_csv(old/'format_normalized_name_ontology.csv',keep_default_na=False)
    predictions=[];mapping_rows=[]
    def add(method,native,broad,endpoint='strict_T_name_rule',information='archived original input/reference conditions'):
        native=np.asarray(native,dtype=str);broad=np.asarray(broad,dtype=str)
        assert len(native)==len(ref) and len(broad)==len(ref)
        unknown=np.asarray([v.lower() in UNKNOWN for v in native]) | (broad=='Unknown')
        pred=broad=='T' if endpoint=='strict_T_name_rule' else broad=='TNK'
        assert not (pred&unknown).any()
        predictions.append((method,endpoint,native,broad,pred,unknown,information))
        for value in sorted(set(native)):
            mask=native==value;parents=set(broad[mask]);assert len(parents)==1
            mapping_rows.append(dict(method=method,endpoint=endpoint,native=value,broad=next(iter(parents)),
                n_cells=int(mask.sum()),unknown=bool(unknown[mask].all()),TCR_used_for_mapping=False))
    for method in ['SCINA_archived','scCATCH_archived','scType_archived','SignacX_archived']:
        vocab=ontology[ontology.method==method].set_index('native').broad
        assert vocab.index.is_unique
        native=ref[method+'_native'];assert set(native)<=set(vocab.index)
        add(method,native,native.map(vocab))
    for method in ['S2_terminal','S3_terminal']:
        native=ref[method+'_native']
        add('DG-scRNA Sup '+method,native,native.map({v:broad_lineage(v) for v in native.unique()}),
            information='saved original final annotation; source-label endpoints are not independent truth')
    states=joined.CellStates
    fine=states.map(lambda v:'T' if v in SIGNAC_T else 'Unknown' if v.lower() in UNKNOWN or v in ['TNK','Immune'] else broad_lineage(v))
    add('SignacX corrected CellStates',states,fine,
        information='saved per-sample SignacX2.2.5 atlas-trained predictions; no new fit')
    coarse=joined.CellTypes
    add('SignacX corrected CellTypes TNK',coarse,coarse.map(lambda v:'TNK' if v=='TNK' else 'Unknown' if v.lower() in UNKNOWN else 'nonTNK'),
        endpoint='coarse_TNK_sensitivity',information='coarse T+NK output; not a strict-T prediction')
    pd.DataFrame(mapping_rows).to_csv(dest/'frozen_native_label_mapping.csv',index=False)
    scopes=[('ALL','ALL','ALL',np.ones(len(ref),dtype=bool))]
    for group,g in ref.groupby('group',sort=True):
        scopes.append((group,'ALL',group,ref.index.isin(g.index)))
        for patient,p in g.groupby('patient',sort=True):
            scopes.append((group,patient,p['sample'].iloc[0],ref.index.isin(p.index)))
    truths={'productive_TCR':ref.TCR_cell_high_confidence_productive_TCR.to_numpy(dtype=bool),
            'paper_original':ref.paper_truth.to_numpy(dtype=bool),
            'any_filtered_contig':ref.TCR_any_filtered_contig.to_numpy(dtype=bool),
            'S3_supplied':ref.TCR_S3_supplied.to_numpy(dtype=bool)}
    rows=[];labels=ref[['sample','patient','group']].copy()
    for method,endpoint,native,broad,pred,unknown,information in predictions:
        labels[method+'_native']=native;labels[method+'_mapped']=broad
        for truth,y in truths.items():
            for group,patient,scope,mask in scopes:
                rows.append(dict(method=method,endpoint=endpoint,truth_definition=truth,group=group,
                    patient=patient,scope=scope,information=information,**binary_metrics(y[mask],pred[mask],unknown[mask])))
    labels.to_csv(dest/'saved_predictions_joined.csv.gz')
    metrics=pd.DataFrame(rows);metrics.to_csv(dest/'all_assay_metrics.csv',index=False)
    # Audit earlier numerical reports without forcing their endpoints to agree.
    memo=metrics[(metrics.method=='SignacX corrected CellTypes TNK')&(metrics.group=='ALL')].copy()
    expected={'TP':33378,'FP':3612,'FN':2349,'TN':53065}
    for key,value in expected.items():memo['old_memo_'+key]=value
    memo['old_memo_confusion_reproduced']=np.logical_and.reduce([memo[k].eq(v).to_numpy() for k,v in expected.items()])
    memo.to_csv(dest/'SignacX_previous_memo_replay.csv',index=False)
    # The old compact-name normalization must reproduce its saved TCR count audit.
    previous=pd.read_csv(old/'comparator_TCR_metrics.csv',keep_default_na=False)
    fieldmap={'TCR_cell_high_confidence_productive_TCR':'productive_TCR',
        'TCR_any_filtered_contig':'any_filtered_contig','TCR_S3_supplied':'S3_supplied'}
    prior=previous[(previous.mapping=='strict_T')&previous.TCR_definition.isin(fieldmap)].copy()
    prior['truth_definition']=prior.TCR_definition.map(fieldmap)
    current=metrics[metrics.method.isin(prior.method)&metrics.endpoint.eq('strict_T_name_rule')]
    current=current[(current.patient!='ALL')|(current.group=='ALL')].copy();current['sample']=current.scope
    parity=current.merge(prior,on=['method','truth_definition','sample'],validate='one_to_one')
    assert len(parity)==len(prior)
    for a,b in [('n_cells','n'),('TP','detected_and_predicted'),('FP','predicted_no_detection'),
                ('FN','detected_not_predicted'),('TN','neither')]:
        assert np.array_equal(parity[a],parity[b]),(a,b)
    assert checked(OUT/'PTC_summary')
    primary=metrics[(metrics.truth_definition=='productive_TCR')&(metrics.endpoint=='strict_T_name_rule')&
        (metrics.patient!='ALL')].copy()
    # Add the already evaluated original-R anchor and label-heldout workflow,
    # with identical productive-TCR truth and the same strict name rule.
    frames=[primary]
    for file,method in [('original_anchor_patient_metrics.csv','DG-scRNA original R anchor'),
                        ('heldout_workflow_patient_metrics.csv','DG-scRNA selected fresh workflow')]:
        d=pd.read_csv(OUT/'PTC_summary'/file,keep_default_na=False)
        assert set(d.truth_definition)=={'productive_TCR'} and set(d.endpoint)=={'strict_T_name_rule'}
        assert len(d)==8
        d['method']=method;d['information']='original-R terminal DL; group-transductive fit; selection scope explicitly named'
        frames.append(d)
    primary=pd.concat(frames,ignore_index=True)
    primary.to_csv(dest/'strict_productive_TCR_patient_comparison.csv',index=False)
    measures=['F1_T','F1_nonT_unknown_as_error','macro_F1_unknown_as_error','accuracy_unknown_as_error',
              'TCR_positive_recall','TCR_detection_yield_in_predicted_T','coverage']
    summary=primary.groupby(['group','method'])[measures].agg(['mean','std','count'])
    summary.columns=['_'.join(k) for k in summary.columns]
    summary.reset_index().to_csv(dest/'strict_productive_TCR_patient_summary.csv',index=False)
    contrasts=[];differences=[]
    for group,d in primary.groupby('group',sort=True):
        for baseline in ['DG-scRNA original R anchor','DG-scRNA selected fresh workflow']:
            b=d[d.method==baseline][['patient','F1_T']].rename(columns={'F1_T':'DG_F1'})
            for method,c in d[~d.method.str.startswith('DG-scRNA')].groupby('method',sort=True):
                p=c.merge(b,on='patient',validate='one_to_one').sort_values('patient')
                delta=p.DG_F1-p.F1_T
                contrasts.append(dict(group=group,DG_configuration=baseline,competitor=method,
                    contrast='DG minus competitor',DG_mean_F1=float(p.DG_F1.mean()),competitor_mean_F1=float(p.F1_T.mean()),
                    absolute_gain_percentage_points=float(100*delta.mean()),
                    relative_gain_percent=float(100*delta.mean()/p.F1_T.mean()) if p.F1_T.mean()>0 else None,
                    information_is_matched=False,**paired(delta)))
                for row,v in zip(p.itertuples(),delta):
                    differences.append(dict(group=group,DG_configuration=baseline,competitor=method,patient=row.patient,delta=float(v)))
    contrasts=pd.DataFrame(contrasts);contrasts['p_Holm_within_group_and_DG_configuration']=1.
    for _,d in contrasts.groupby(['group','DG_configuration']):
        ix=d.sort_values('exact_sign_flip_p',kind='stable').index
        contrasts.loc[ix,'p_Holm_within_group_and_DG_configuration']=np.minimum(1,np.maximum.accumulate(
            contrasts.loc[ix,'exact_sign_flip_p'].to_numpy()*np.arange(len(ix),0,-1)))
    contrasts.to_csv(dest/'paired_patient_gains.csv',index=False)
    pd.DataFrame(differences).to_csv(dest/'paired_patient_differences.csv',index=False)
    fig,axes=plt.subplots(1,2,figsize=(13,6),layout='constrained')
    methods=list(primary.method.drop_duplicates())
    for ax,group in zip(axes,['NMT','TTU']):
        d=primary[primary.group==group]
        for i,method in enumerate(methods):
            y=d[d.method==method].sort_values('patient').F1_T.to_numpy();assert len(y)==4
            ax.scatter(np.arange(4)*.07+i-.1,y,s=14,color='#4477AA',alpha=.7)
            ax.plot([i-.25,i+.25],[y.mean()]*2,color='#CC6677',linewidth=2)
        ax.set_xticks(range(len(methods)),[m.replace('DG-scRNA ','DG ').replace('SignacX ','SignacX\n').replace('_archived','\narchived') for m in methods],rotation=60,ha='right',fontsize=7)
        ax.set(title=group,xlabel='Saved outputs: original conditions differ',ylabel='Strict-T / productive-TCR detection F1',ylim=(0,1))
    fig.suptitle('Four patient values and their mean; equal endpoint does not imply equal prior information')
    for ext in ['png','pdf']:fig.savefig(dest/f'PTC_cached_comparators.{ext}',dpi=210,bbox_inches='tight')
    plt.close(fig)
    note='''# PTC cached competitors: corrected endpoints and SignacX interpretation

No model is fitted in this audit. The corrected SignacX concatenated output is
verified against all eight saved per-sample files and joined by sample plus barcode
to all92,404 evaluation cells. Cells excluded by the historical PTC cohort remain
in the source and are not silently introduced into its denominator.

The prior0.9180 memo used CellTypes==TNK. TNK combines T and NK. The saved finer
CellStates output supports a separate explicit-T analysis using the method's own
T-prefixed state identities. NK is not a strict-T call. Unclassified and unresolved
TNK/Immune states remain Unknown for the strict endpoint. These mappings use names,
not TCR outcomes; both native columns and every mapping decision are exported.

The original S3 SignacX column was populated but contained no lymphoid calls. It
was not an absent or crashed prediction vector. Its cause has not been isolated.
The archived wrapper's comments attributing this to an Assay5 incompatibility or
missing graph are hypotheses, not established explanations. The old figure's
comparison of SignacX T-positive F1 with published DG non-T-positive F1, and its
description of the S2 flag as productive TCR, must not be used in a revised paper.
This audit supersedes those interpretations and preserves the historical files.

All methods here use the same full-cell denominator and truth definition within
each table; Unknown-aware metrics and coverage are explicit. TCR non-detection is
not biological proof of non-T identity. Historical paper broad-T compatibility,
coarse TNK, strict T, productive detection and any-contig detection remain separate.
Prior numerical memo reproduction is checked, not assumed.

Primary descriptive comparisons retain four paired patients per NMT/TTU group,
both DG original anchors and patient-label-selected fresh workflows, absolute and
relative gains, all patient differences and conditional uncertainty. Inputs,
pooling and labelled reference information differ across these historical fits.
They are not a matched-information causal benchmark. Four patients give a minimum
two-sided sign-flip p of0.125; no stronger replication is invented from cells.
'''
    (dest/'ENDPOINT_CORRECTIONS.md').write_text(note)
    write_json(dest/'manifest.json',dict(status='completed',new_fits=0,n_input_SignacX=len(data),n_joined=len(ref),
        concatenated_output_matches_all_eight_files=True,archived_count_parity_rows=len(parity),
        SignacX_version='2.2.5',SignacX_label_source=SIGNAC_SOURCE,explicit_T_states=sorted(SIGNAC_T),
        sources={str(p):sha(p) for p in [source,*parts,old/'format_normalized_name_ontology.csv',
            REFERENCE/'manifest.json',PAPER/'original_DG_binary_pairs_for_R.csv.gz',
            ROOT/'handoff/ptc/run_signacx.R',ROOT/'handoff/ptc/fig_signacx.py']},
        old_SignacX_failure_cause='unestablished',strict_T_and_coarse_TNK_separate=True,
        all_cells_denominator=True,original_conditions_differ=True,
        source_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest)

if __name__=='__main__':run()
