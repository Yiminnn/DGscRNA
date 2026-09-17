"""Complete census and patient-level analyses of independently verified PTC outputs."""
import itertools
import json
import os
import time
from pathlib import Path
from ptc_common import BASE,RECOVERY,GROUPS,require_slurm,sha,utc,write_json

FEATURES=['all','hvg500','hvg1000','hvg2000','hvg3000','hvg5000']
ENDPOINTS=['cluster_S2_native_ARI','S2_macro_F1_reference_present',
           'productive_strict_recall','productive_strict_detection_yield',
           'productive_strict_apparent_F1','coverage']
EXTRA=['S3_macro_F1_reference_present','S2_broad_concordance','S3_broad_concordance',
       'any_contig_strict_recall','any_contig_strict_apparent_F1',
       'S3_supplied_strict_apparent_F1','paired_TRA_TRB_strict_recall',
       'productive_permissive_recall','productive_permissive_apparent_F1',
       'undetected_predicted_T_RNA_T_support','Unknown_RNA_T_support',
       'Unknown_median_nCount','Unknown_median_nFeature','Unknown_median_percent_mt',
       'n_NK','n_NKT','n_ambiguous_lymphoid','noise_fraction']
KEYS=['mode','feature','input_space','dr','dim','seed','correction','clusterer','space',
      'scoring_assay','library','cutoff','stage','family']

def effect_values(values):
    import numpy as np
    v=np.asarray(values,dtype=float);v=v[np.isfinite(v)];n=len(v)
    if not n:return dict(n_paired_patients=0,delta_mean=np.nan,ci_lower=np.nan,ci_upper=np.nan,p_exact=np.nan)
    assert n<=4
    delta=float(v.mean())
    flips=np.asarray(list(itertools.product([-1,1],repeat=n)))
    p=float((np.abs((flips*v).mean(1))>=abs(delta)-1e-12).mean())
    draws=np.asarray(list(itertools.product(range(n),repeat=n)))
    lo,hi=np.percentile(v[draws].mean(1),[2.5,97.5])
    return dict(n_paired_patients=n,delta_mean=delta,ci_lower=float(lo),ci_upper=float(hi),p_exact=p,
        n_positive_patients=int((v>0).sum()),patient_differences=json.dumps(v.tolist()))

def holm(frame,groups):
    import numpy as np
    frame=frame.copy();frame['p_holm']=np.nan
    if frame.empty:return frame
    for _,g in frame.groupby(groups,dropna=False):
        ids=g.index[g.p_exact.notna()];ids=frame.loc[ids,'p_exact'].sort_values().index
        p=frame.loc[ids,'p_exact'].to_numpy();n=len(p)
        frame.loc[ids,'p_holm']=np.minimum(1,np.maximum.accumulate(p*(n-np.arange(n))))
    return frame

def paired(candidate,baseline,name,family,extra=None,metrics=ENDPOINTS):
    import pandas as pd
    keys=['sample','patient']
    assert not candidate.duplicated(keys).any(),(name,'candidate duplicate')
    assert not baseline.duplicated(keys).any(),(name,'baseline duplicate')
    a=candidate[keys+metrics].merge(baseline[keys+metrics],on=keys,suffixes=('_candidate','_baseline'),validate='one_to_one')
    result=[]
    for metric in metrics:
        q=a[keys].copy();q['candidate']=a[metric+'_candidate'];q['baseline']=a[metric+'_baseline'];q=q.dropna()
        q['delta']=q.candidate-q.baseline;p=q.groupby('patient')[['candidate','baseline','delta']].mean()
        result.append(dict(contrast=name,family=family,metric=metric,n_paired_samples=len(q),
           candidate_patient_mean=p.candidate.mean(),baseline_patient_mean=p.baseline.mean(),
           patient_ids=json.dumps(p.index.to_list()),**effect_values(p.delta),**(extra or {})))
    return result

def summaries(df,metrics=ENDPOINTS+EXTRA):
    import pandas as pd
    keys=KEYS+['analysis_scope']
    scoped=pd.concat([df.assign(analysis_scope='ALL_8_SAMPLES'),df.assign(analysis_scope=df.group)],ignore_index=True)
    patient=scoped.groupby(keys+['patient'],dropna=False)[metrics].mean().reset_index()
    means=patient.groupby(keys,dropna=False)[metrics].mean().reset_index()
    count=scoped.groupby(keys,dropna=False).agg(n_samples=('sample','nunique'),n_patients=('patient','nunique'),
              n_sample_rows=('sample','size')).reset_index()
    return means.merge(count,on=keys,validate='one_to_one')

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    smoke=os.environ.get('PTC_ANALYSIS_SMOKE')=='1'
    dest=BASE/('summary_smoke' if smoke else 'summary');dest.mkdir(exist_ok=True)
    while not smoke and any(not (BASE/'evaluations'/f'task_{i:03d}'/'COMPLETE').exists() for i in range(510)):
        time.sleep(30)
    if not smoke:assert (BASE/'verification/geometry_census/COMPLETE').exists(),'Run after geometry census completes'
    manifests=[];tables={name:[] for name in ['sample_metrics','condition_verification','cluster_metrics','lineage_composition']}
    for i in range(510):
        if smoke and not (BASE/'evaluations'/f'task_{i:03d}'/'COMPLETE').exists():continue
        d=BASE/'evaluations'/f'task_{i:03d}';m=json.loads((d/'manifest.json').read_text())
        assert m['task_index']==i and (d/'COMPLETE').read_text().strip()==sha(d/'manifest.json')
        for name,h in m['outputs'].items():assert sha(d/name)==h
        manifests.append(dict(task=i,manifest_sha256=sha(d/'manifest.json'),**m))
        for name in tables:tables[name].append(pd.read_csv(d/(name+'.csv.gz'),keep_default_na=False))
    metrics,conditions,cl,composition=[pd.concat(tables[n],ignore_index=True) for n in tables]
    assert conditions.condition_id.is_unique
    assert not metrics.duplicated(['condition_id','sample','stage']).any()
    if not smoke:assert len(conditions)==11696 and len(metrics)==49776
    assert set(metrics.stage)=={'marker_only_ablation','terminal_DL090','terminal_DL070_sensitivity'}
    if not smoke:assert conditions.groupby('mode').size().to_dict()=={'matched':1632,'pooled':1632,'single':8432}
    # CCA RNA/integrated scoring share the same physical clustering; their
    # context-family labels differ but their cell partitions must be identical.
    ck=['partition_path','sample']
    cv=['cluster_S2_native_ARI','cluster_S2_native_V','noise_fraction','n_clusters_in_sample','n_cells']
    assert cl.groupby(ck)[cv].nunique(dropna=False).le(1).all().all()
    cl=cl.drop(columns=['scoring_assay']).drop_duplicates(ck)
    assert not cl.duplicated(['partition_path','sample']).any()
    if not smoke:assert cl.partition_path.nunique()==1288
    metrics=metrics.merge(cl[['partition_path','sample','cluster_S2_native_ARI','cluster_S2_native_V',
                             'noise_fraction','n_clusters_in_sample']],on=['partition_path','sample'],validate='many_to_one')
    for col in ENDPOINTS+EXTRA:metrics[col]=pd.to_numeric(metrics[col],errors='coerce')
    for col in ['dim','seed']:metrics[col]=pd.to_numeric(metrics[col],errors='coerce').fillna(0).astype(int)
    for frame in [metrics,conditions]:frame['precision_recovery']=frame.precision_recovery.astype(str).str.lower().eq('true')
    assert set(metrics['sample'])==set(sum(GROUPS.values(),[])) and metrics.patient.nunique()==4
    cohort=pd.read_csv(BASE/'evaluation_reference/cohort_summary.csv').set_index('sample')
    assert metrics.n_cells.eq(metrics['sample'].map(cohort.n_cells)).all()
    assert metrics.patient.eq(metrics['sample'].map(cohort.patient)).all()
    assert metrics.group.eq(metrics['sample'].map(cohort.group)).all()
    for key in conditions.cache_key.unique():assert (BASE/'verification/cache'/f'{key}.json').exists(),key
    metrics.to_csv(dest/'all_sample_stage_metrics.csv.gz',index=False)
    conditions.to_csv(dest/'condition_verification.csv.gz',index=False)
    cl.to_csv(dest/'cluster_metrics.csv.gz',index=False);composition.to_csv(dest/'lineage_composition.csv.gz',index=False)
    summary=summaries(metrics);summary.to_csv(dest/'condition_summary.csv.gz',index=False)
    terminal=metrics[metrics.stage.eq('terminal_DL090')]
    primary=terminal[terminal['mode'].eq('single')&terminal.library.eq('CellMarker_AllTissues')&
                     terminal.cutoff.eq('mean')&terminal.seed.eq(42)]
    factorial=primary[primary.family.eq('factorial')]
    if not smoke:assert len(factorial)==1008
    factorial.to_csv(dest/'primary_factorial_by_sample.csv',index=False)
    summaries(factorial).to_csv(dest/'primary_factorial_patient_means.csv',index=False)
    effects=[]
    for path,dr,dim,inp in [('direct_UMAP2','UMAP',2,'direct'),('PCA30_UMAP2','UMAP',2,'pca30'),('PCA30','PCA',30,'direct')]:
        z=primary[primary.dr.eq(dr)&primary.dim.eq(dim)&primary.input_space.eq(inp)&primary.clusterer.eq('HDBSCAN')]
        for f in FEATURES[1:]:
            if z.feature.eq(f).any():
                effects+=paired(z[z.feature.eq(f)],z[z.feature.eq('all')],f+' minus all','HVG',dict(path=path,feature=f))
    hvg=holm(pd.DataFrame(effects),['path','metric'])
    hvg['p_holm_2_primary_endpoints']=np.nan
    ix=hvg.path.eq('direct_UMAP2')&hvg.feature.eq('hvg2000')&hvg.metric.isin(ENDPOINTS[:2])
    adjusted=holm(hvg[ix].assign(primary_family='two_primary_endpoints'),['primary_family'])
    hvg.loc[adjusted.index,'p_holm_2_primary_endpoints']=adjusted.p_holm
    hvg.to_csv(dest/'HVG_paired_patient_effects.csv',index=False)
    group_effects=[]
    for group in GROUPS:
      for path,dr,dim,inp in [('direct_UMAP2','UMAP',2,'direct'),('PCA30_UMAP2','UMAP',2,'pca30'),('PCA30','PCA',30,'direct')]:
        z=primary[primary.group.eq(group)&primary.dr.eq(dr)&primary.dim.eq(dim)&primary.input_space.eq(inp)&primary.clusterer.eq('HDBSCAN')]
        for feature in FEATURES[1:]:
            if z.feature.eq(feature).any():
                group_effects+=paired(z[z.feature.eq(feature)],z[z.feature.eq('all')],feature+' minus all','HVG_by_group',dict(group=group,path=path,feature=feature))
    holm(pd.DataFrame(group_effects),['group','path','metric']).to_csv(dest/'HVG_by_group_patient_effects.csv',index=False)
    interaction=[]
    u=primary[primary.dr.eq('UMAP')&primary.clusterer.eq('HDBSCAN')&primary.feature.isin(['all','hvg2000'])]
    for metric in ENDPOINTS:
        q=u.pivot(index=['sample','patient'],columns=['input_space','feature'],values=metric).dropna()
        dif=(q['direct','hvg2000']-q['direct','all'])-(q['pca30','hvg2000']-q['pca30','all'])
        p=dif.groupby(level='patient').mean()
        interaction.append(dict(metric=metric,contrast='HVG effect direct minus after PCA30',n_paired_samples=len(q),
                                **effect_values(p),patient_ids=json.dumps(p.index.to_list())))
    holm(pd.DataFrame(interaction).assign(family='HVG_by_PCA'),['family']).to_csv(dest/'HVG_by_PCA_interaction.csv',index=False)
    methods=[]
    for feature in ['all','hvg2000']:
        z=factorial[factorial.feature.eq(feature)]
        a=z[z.dr.eq('UMAP')&z.clusterer.eq('HDBSCAN')]
        for (dr,clusterer),b in z.groupby(['dr','clusterer']):
            if (dr,clusterer)==('UMAP','HDBSCAN'):continue
            methods+=paired(a,b,'UMAP/HDBSCAN minus '+dr+'/'+clusterer,'method',dict(feature=feature,comparator=dr+'/'+clusterer))
    holm(pd.DataFrame(methods),['feature','metric']).to_csv(dest/'paired_method_comparisons.csv',index=False)
    seed=terminal[terminal['mode'].eq('single')&terminal.dr.eq('UMAP')&terminal.clusterer.eq('HDBSCAN')&
                  terminal.library.eq('CellMarker_AllTissues')&terminal.cutoff.eq('mean')&terminal.feature.isin(['all','hvg2000'])]
    seed.groupby(['sample','patient','feature','input_space'])[ENDPOINTS].agg(['mean','std','min','max','count']).to_csv(dest/'seed_variability.csv')
    recovery=[]
    affected=primary[primary.precision_recovery]
    for a in affected.itertuples():
        z=primary[primary.feature.eq(a.feature)&primary.dr.eq(a.dr)&primary.clusterer.eq(a.clusterer)&primary.dim.eq(a.dim)]
        for metric in ENDPOINTS:
            recovery.append(dict(feature=a.feature,dr=a.dr,clusterer=a.clusterer,metric=metric,recovered_sample=a.sample,
               complete_patient_mean=z.groupby('patient')[metric].mean().mean(),n_samples=len(z),
               without_recovered_sample_patient_mean=z[~z['sample'].eq(a.sample)].groupby('patient')[metric].mean().mean(),
               n_samples_without_recovery=int((~z['sample'].eq(a.sample)).sum())))
    pd.DataFrame(recovery).to_csv(dest/'numerical_recovery_sensitivity.csv',index=False)
    batch=[];legacy=[]
    z=terminal[terminal['mode'].isin(['pooled','matched'])&terminal.library.eq('CellMarker_AllTissues')&terminal.cutoff.eq('mean')]
    for group in GROUPS:
      for space in ['PCA30','UMAP2']:
       for method in ['SNN','HDBSCAN_R']:
        q=z[z.group.eq(group)&z.space.eq(space)&z.clusterer.eq(method)]
        rna=q[q.scoring_assay.eq('RNA')];none=rna[rna.correction.eq('NONE')];single=rna[rna['mode'].eq('matched')]
        ex=dict(group=group,space=space,clusterer=method)
        batch+=paired(none,single,'pooled NONE minus matched single','pooling',ex)
        for correction in ['CCA','HARMONY']:
            batch+=paired(rna[rna.correction.eq(correction)],none,correction+' minus NONE','correction',ex)
            batch+=paired(rna[rna.correction.eq(correction)],single,correction+' minus matched single','pooling_plus_correction',ex)
        legacy+=paired(q[q.scoring_assay.eq('integrated')],rna[rna.correction.eq('CCA')],
                       'CCA integrated scoring/DL minus full RNA scoring/DL','assay',ex)
    holm(pd.DataFrame(batch),['group','space','clusterer','metric']).to_csv(dest/'matched_batch_effects.csv',index=False)
    holm(pd.DataFrame(legacy),['group','metric']).to_csv(dest/'legacy_integrated_effects.csv',index=False)
    # Refinement and threshold interventions use the very same initial condition.
    delta=[]
    for stage,name in [('marker_only_ablation','terminal090 minus marker-only'),('terminal_DL070_sensitivity','terminal070 minus terminal090')]:
        other=metrics[metrics.stage.eq(stage)].set_index(['condition_id','sample'])
        main=terminal.set_index(['condition_id','sample']);assert set(main.index)==set(other.index)
        for metric in ENDPOINTS[1:]+['undetected_predicted_T_RNA_T_support']:
            q=main[KEYS[:-2]+['patient','group','family']].copy()
            q['metric']=metric
            q['delta']=(main[metric]-other[metric]) if stage=='marker_only_ablation' else (other[metric]-main[metric])
            q['intervention']=name;delta.append(q.reset_index())
    delta=pd.concat(delta,ignore_index=True)
    dkeys=[x for x in KEYS if x!='stage']+['group','metric','intervention']
    dr=delta.groupby(dkeys+['patient'],dropna=False).delta.mean().reset_index()
    dr.groupby(dkeys,dropna=False).delta.agg(['mean','min','max','count']).reset_index().to_csv(dest/'refinement_and_threshold_effects.csv.gz',index=False)
    # Select annotation context on three patients, evaluate only the held-out patient.
    selection=[]
    panels=terminal[terminal['mode'].eq('pooled')&terminal.scoring_assay.eq('RNA')]
    endpoints=ENDPOINTS+['undetected_predicted_T_RNA_T_support']
    tie={'mean':0,'none':1,'0.5':2}
    for keys,q in panels.groupby(['group','correction','space','clusterer']):
        for patient in sorted(q.patient.unique()):
            train=q[~q.patient.eq(patient)].groupby(['library','cutoff','patient']).productive_strict_apparent_F1.mean().groupby(['library','cutoff']).mean().reset_index()
            train['cutoff_order']=train.cutoff.map(tie)
            chosen=train.sort_values(['productive_strict_apparent_F1','library','cutoff_order'],ascending=[False,True,True]).iloc[0]
            policies=[('heldout_selected',chosen.library,chosen.cutoff)]
            for lib in ['CellMarker_AllTissues','CellMarker_Thyroid','Pubmed_34663816','HPA_allThyroid']:
                policies.append((lib+'_mean',lib,'mean'))
            for policy,lib,cutoff in policies:
                held=q[q.patient.eq(patient)&q.library.eq(lib)&q.cutoff.eq(cutoff)]
                assert len(held)==1,(keys,patient,lib,cutoff,len(held))
                selection.append(dict(zip(['group','correction','space','clusterer'],keys),patient=patient,policy=policy,
                    library=lib,cutoff=cutoff,training_patient_count=3,
                    selected_training_apparent_F1=chosen.productive_strict_apparent_F1,
                    **held[endpoints].iloc[0].to_dict()))
    selection=pd.DataFrame(selection);selection.to_csv(dest/'patient_heldout_panel_selection.csv',index=False)
    selection.groupby(['group','correction','space','clusterer','policy'])[endpoints].mean().reset_index().to_csv(dest/'patient_heldout_panel_summary.csv',index=False)
    # Sample/QC evidence is reconstructed independently from raw 10x files.
    audits=[]
    for sample in sum(GROUPS.values(),[]):
        p=BASE/'raw_QC_audit_symbol_normalized'/f'{sample}.json';a=json.loads(p.read_text())
        assert a['all_archived_genes_compared']
        for k in ['n_retained_cells_raw_count_exact','n_retained_cells_raw_nCount_exact','n_retained_cells_raw_nFeature_exact','n_retained_cells_raw_mt_exact','raw_candidate_expected_after_doublets']:
            assert a[k]==a['n_archived_S2']
        a['gene3_first_QC_shortfall_vs_raw_QC']=a['n_after_raw_reference_QC']-a['n_after_reference_QC']
        a['raw_expected_doublets']=a['n_after_raw_reference_QC']-a['raw_candidate_expected_after_doublets']
        audits.append({k:v for k,v in a.items() if isinstance(v,(str,int,float,bool))})
    pd.DataFrame(audits).to_csv(dest/'raw_counts_QC_parity.csv',index=False)
    # Training-condition counts differ from unique trained cache/model counts.
    conditions.groupby(['mode','dl_status','actual_model_trained']).size().rename('conditions').reset_index().to_csv(dest/'terminal_execution_census.csv',index=False)
    cache=conditions.drop_duplicates('cache_key')
    cache.groupby(['dl_status','actual_model_trained']).size().rename('unique_caches').reset_index().to_csv(dest/'unique_model_census.csv',index=False)
    for key in cache.cache_key:
        p=BASE/'verification/cache'/f'{key}.json';assert p.exists() and json.loads(p.read_text())
    write_json(dest/'evaluation_manifest_inventory.json',manifests)
    ntrained=int(cache.actual_model_trained.astype(str).str.lower().eq('true').sum())
    manifest=dict(status='preliminary_partial_diagnostic' if smoke else 'complete',n_evaluation_units=len(manifests),n_terminal_conditions=len(conditions),n_sample_stage_rows=len(metrics),
       n_samples=8,n_patients=4,n_cells=92404,n_physical_partitions=cl.partition_path.nunique(),n_unique_caches=len(cache),n_unique_trained_models=ntrained,
       n_precision_recovered_conditions=int(conditions.precision_recovery.sum()),
       statistics='Sample differences averaged within patient; exact paired sign flips and exhaustive size-4 patient bootstrap; CIs descriptive',
       panel_selection='Patient-held-out panel selection on transductively fitted predictions; not fully inductive validation',
       source_sha256=sha(Path(__file__)),analysis_plan_sha256=sha(Path(__file__).with_name('FINAL_ANALYSIS_PLAN.md')),
       geometry_census_sha256=sha(BASE/'verification/geometry_census/manifest.json') if (BASE/'verification/geometry_census/manifest.json').exists() else None,
       outputs={p.name:sha(p) for p in dest.iterdir() if p.is_file() and p.name not in ['manifest.json','COMPLETE']},
       completed_at=utc(),job=os.environ['SLURM_JOB_ID'])
    write_json(dest/'manifest.json',manifest);(dest/('SMOKE_COMPLETE' if smoke else 'COMPLETE')).write_text(sha(dest/'manifest.json')+'\n')
    print(json.dumps({k:v for k,v in manifest.items() if k!='outputs'},indent=2),flush=True)

if __name__=='__main__':run()
