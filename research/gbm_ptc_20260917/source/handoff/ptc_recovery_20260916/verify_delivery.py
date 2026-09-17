"""Final artifact, independently recomputed statistics and notebook delivery gate."""
import concurrent.futures
import itertools
import json
import os
from pathlib import Path
import statistics
import subprocess
import shutil
from ptc_common import BASE,ROOT,RECOVERY,task_list,require_slurm,sha,utc,write_json

def run():
    require_slurm()
    import numpy as np
    import pandas as pd
    import nbformat
    dest=BASE/'verification/final_delivery';dest.mkdir(exist_ok=True)
    sm=json.loads((BASE/'summary/manifest.json').read_text())
    assert sm['status']=='complete' and sm['n_terminal_conditions']==11696 and sm['n_sample_stage_rows']==49776
    for name,digest in sm['outputs'].items():assert sha(BASE/'summary'/name)==digest,name
    conditions=pd.read_csv(BASE/'summary/condition_verification.csv.gz')
    assert conditions.condition_id.nunique()==11696 and len(conditions)==11696
    # Recheck every unique saved model artifact against its already passed
    # independent full-forward verification record, without reusing fit code.
    def cache_check(key):
        d=BASE/'refinement_cache'/key;m=json.loads((d/'training_manifest.json').read_text())
        v=json.loads((BASE/'verification/cache'/f'{key}.json').read_text())
        assert v['status']=='passed' and v['cache_manifest_sha256']==sha(d/'training_manifest.json')
        assert v['model_forward_checked']==m['training_executed']
        for name,h in m['outputs'].items():assert sha(d/name)==h,(key,name)
        return dict(cache_key=key,training_executed=m['training_executed'],n_cells=m['n_cells'],
                    verification_sha256=sha(BASE/'verification/cache'/f'{key}.json'))
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:
        caches=list(ex.map(cache_check,conditions.cache_key.unique()))
    assert len(caches)==sm['n_unique_caches'] and sum(c['training_executed'] for c in caches)==sm['n_unique_trained_models']
    pd.DataFrame(caches).to_csv(dest/'all_cache_integrity.csv.gz',index=False)
    # Independently read the fixed primary evaluation files, not the aggregate.
    pieces=[]
    for task in task_list():
        if task['seed']!=42 or task['feature'] not in ['all','hvg2000'] or task['dr']!='UMAP':continue
        d=BASE/'evaluations'/f'task_{task["task_id"]:03d}'
        m=json.loads((d/'manifest.json').read_text())
        assert sha(d/'sample_metrics.csv.gz')==m['outputs']['sample_metrics.csv.gz']
        q=pd.read_csv(d/'sample_metrics.csv.gz')
        q=q[q.stage.eq('terminal_DL090')&q.library.eq('CellMarker_AllTissues')&q.cutoff.eq('mean')&q.clusterer.eq('HDBSCAN')]
        assert len(q)==1
        c=pd.read_csv(d/'cluster_metrics.csv.gz');c=c[c.clusterer.eq('HDBSCAN')];assert len(c)==1
        q=q.copy();q['cluster_S2_native_ARI']=c.cluster_S2_native_ARI.iloc[0];pieces.append(q)
    raw=pd.concat(pieces,ignore_index=True)
    declared=pd.read_csv(BASE/'summary/HVG_paired_patient_effects.csv');checks=[]
    for path,inp in [('direct_UMAP2','direct'),('PCA30_UMAP2','pca30')]:
      z=raw[raw.input_space.eq(inp)]
      for metric in ['cluster_S2_native_ARI','S2_macro_F1_reference_present','productive_strict_recall','productive_strict_detection_yield','productive_strict_apparent_F1','coverage']:
        wide=z.pivot(index=['sample','patient'],columns='feature',values=metric).dropna()
        sample_diffs=wide.hvg2000-wide['all'];patient=sample_diffs.groupby(level='patient').mean().sort_index()
        if metric!='productive_strict_detection_yield':assert len(wide)==8 and len(patient)==4
        values=patient.to_list();n=len(values)
        if n:
            mean=statistics.fmean(values)
            statistics_under_null=[statistics.fmean([s*x for s,x in zip(signs,values)]) for signs in itertools.product([-1,1],repeat=n)]
            p=sum(abs(x)>=abs(mean)-1e-12 for x in statistics_under_null)/(2**n)
            draws=[statistics.fmean(x) for x in itertools.product(values,repeat=n)]
            ci=np.quantile(draws,[.025,.975])
        else:mean=p=np.nan;ci=[np.nan,np.nan]
        r=declared[declared.path.eq(path)&declared.feature.eq('hvg2000')&declared.metric.eq(metric)].iloc[0]
        np.testing.assert_allclose([r.delta_mean,r.ci_lower,r.ci_upper,r.p_exact],[mean,*ci,p],rtol=0,atol=1e-12)
        assert r.n_paired_samples==len(wide) and r.n_paired_patients==len(patient)
        checks.append(dict(path=path,metric=metric,n_paired_samples=len(wide),n_paired_patients=n,delta_mean=mean,ci_lower=ci[0],ci_upper=ci[1],p_exact=p,passed=True))
    pd.DataFrame(checks).to_csv(dest/'independent_primary_statistics.csv',index=False)
    # Saved cell compositions and RNA/QC tables must agree with evaluated counts.
    detail=json.loads((BASE/'biology_detail/manifest.json').read_text())
    assert detail['status']=='complete'
    for name,h in detail['outputs'].items():assert sha(BASE/'biology_detail'/name)==h
    vocabulary=json.loads((BASE/'vocabulary_audit/manifest.json').read_text())
    assert vocabulary['status']=='complete' and not vocabulary['prediction_labels_changed']
    for name,h in vocabulary['outputs'].items():assert sha(BASE/'vocabulary_audit'/name)==h
    comp=pd.read_csv(BASE/'biology_detail/lineage_composition_by_sample.csv')
    expected=pd.read_csv(BASE/'evaluation_reference/cohort_summary.csv')
    allmetrics=pd.read_csv(BASE/'summary/all_sample_stage_metrics.csv.gz',usecols=['sample','patient','group','n_cells','library','S2_macro_F1_reference_present','S3_macro_F1_reference_present'])
    expected=expected.set_index('sample')
    for name in ['patient','group','n_cells']:
        assert allmetrics[name].eq(allmetrics['sample'].map(expected[name])).all()
    capacity=pd.read_csv(BASE/'vocabulary_audit/reference_vocabulary_capacity_by_sample.csv')
    for reference in ['S2','S3']:
        ceiling=capacity[capacity.reference.eq(reference)][['sample','library','theoretical_macro_F1_present_ceiling']]
        check=allmetrics.merge(ceiling,on=['sample','library'],validate='many_to_one')
        assert len(check)==len(allmetrics)
        assert check[reference+'_macro_F1_reference_present'].le(check.theoretical_macro_F1_present_ceiling+1e-12).all()
    actual=comp.groupby(['sample','correction','stage']).n.sum()
    assert actual.groupby(level='sample').nunique().eq(1).all()
    assert actual.groupby(level='sample').first().sum()==92404
    assert np.allclose(comp.groupby(['sample','correction','stage']).fraction.sum(),1)
    # All figure formats and the explicit visual review refer to these same files.
    fm=json.loads((BASE/'figures/figure_manifest.json').read_text())
    review=json.loads((BASE/'figures/VISUAL_REVIEW_COMPLETED.json').read_text())
    assert fm['status']=='complete' and review['figure_manifest_sha256']==sha(BASE/'figures/figure_manifest.json')
    assert set(review['reviewed_figures'])=={x['name'] for x in fm['figures']}
    assert fm['n_figures']==len(fm['figures'])
    for f in fm['figures']:
        for suffix,h in f['files'].items():assert sha(BASE/'figures'/(f['name']+'.'+suffix))==h
    for name in ['analysis_preflight','geometry_census','optimization_diagnostics']:
        assert (BASE/'verification'/name/'COMPLETE').exists()
    bookdir=BASE.parent/'notebooks';nm=json.loads((bookdir/'ptc_execution_manifest.json').read_text())
    assert nm['status']=='complete' and nm['preserved_GBM_archive_cells']==98
    assert sha(nm['notebook'])==nm['notebook_sha256']==sha(nm['workspace_notebook'])
    book=nbformat.read(nm['notebook'],as_version=4)
    assert len(book.cells)==nm['n_cells']
    assert not any(o.output_type=='error' for c in book.cells if c.cell_type=='code' for o in c.get('outputs',[]))
    old=nbformat.read(bookdir/'dgscrna_v6_results.ipynb',as_version=4)
    original_ids={c.id for c in old.cells};preserved=[c for c in book.cells if c.id in original_ids]
    assert len(preserved)==98 and preserved==list(old.cells)
    rm=json.loads((BASE/'report_manifest.json').read_text())
    for name,h in rm['outputs'].items():assert sha(BASE/name)==h
    # Stage GBM remains sealed; v7 supersedes only the central entry-point copy.
    g=BASE.parent;gm=json.loads((g/'GBM_DELIVERY.json').read_text())
    for key,relative in [('report_sha256','GBM_REPORT.md'),('pi_brief_sha256','PI_BRIEF_ZH.md'),
                         ('summary_sha256','summary/manifest.json'),('verification_sha256','verification/independent_manifest.json'),
                         ('notebook_manifest_sha256','notebooks/execution_manifest.json')]:assert sha(g/relative)==gm[key]
    # SLURM accounting is strictly limited to the recorded jobs for this task.
    ledger=json.loads((RECOVERY/'job_ledger.json').read_text())
    ids=sorted({str(j) for v in ledger['jobs'].values() for j in [v.get('job',''),*v.get('prior_jobs',[])] if str(j).isdigit()})
    result=subprocess.run(['sacct','-j',','.join(ids),'--starttime','2026-09-16T00:00:00','-n','-P','-o',
        'JobIDRaw,JobID%40,JobName%40,State,ExitCode,ElapsedRaw,AllocCPUS,CPUTimeRAW,MaxRSS,ReqMem'],capture_output=True,text=True,check=True)
    (dest/'SLURM_accounting.psv').write_text('JobIDRaw|JobID|JobName|State|ExitCode|ElapsedRaw|AllocCPUS|CPUTimeRAW|MaxRSS|ReqMem|\n'+result.stdout)
    rows=[line.split('|')[:10] for line in result.stdout.splitlines() if line.strip()]
    def gib(value):
        if not value:return 0.0
        suffix=value[-1];scale={'K':2**10,'M':2**20,'G':2**30,'T':2**40}.get(suffix,1)
        return float(value[:-1] if suffix in 'KMGT' else value)*scale/2**30
    jobs=[r for r in rows if '.' not in r[0]]
    resources=dict(allocated_CPU_hours_at_snapshot=sum(float(r[7] or 0) for r in jobs)/3600,
       max_reported_RSS_GiB=max(gib(r[8]) for r in rows),
       OOM_job_records=[r[1] for r in jobs if r[3].startswith('OUT_OF_MEMORY')],
       state_counts=pd.Series([r[3] for r in jobs]).value_counts().to_dict(),
       note='Includes preserved failed/canceled attempts and controllers; failures with validated retries are not missing scientific conditions. Final verifier itself may be running in this accounting snapshot.')
    write_json(dest/'resource_summary.json',resources)
    source=Path(__file__).parent;snapshot=BASE/'source_snapshot/delivery_v1';snapshot.mkdir(exist_ok=True)
    source_files={}
    for original in sorted(source.iterdir()):
        if original.is_file() and original.suffix in ['.py','.R','.sbatch','.md']:
            target=snapshot/original.name
            if target.exists():assert sha(target)==sha(original),'Delivery source snapshot already exists with different bytes'
            else:shutil.copy2(original,target)
            source_files[original.name]=sha(target)
    write_json(snapshot/'manifest.json',dict(status='complete',files=source_files,
        note='Final source snapshot; earlier frozen fitting/refinement/scoring implementations remain in their original versioned directories.',
        completed_at=utc()))
    manifest=dict(status='passed',n_terminal_conditions=11696,n_sample_stage_rows=49776,n_unique_caches=len(caches),
        n_independently_recomputed_primary_contrasts=len(checks),n_preserved_GBM_cells=98,
        n_executed_notebook_cells=len(book.cells),n_figures=len(fm['figures']),
        summary_sha256=sha(BASE/'summary/manifest.json'),figure_manifest_sha256=sha(BASE/'figures/figure_manifest.json'),
        notebook_manifest_sha256=sha(bookdir/'ptc_execution_manifest.json'),report_manifest_sha256=sha(BASE/'report_manifest.json'),
        final_source_snapshot_sha256=sha(snapshot/'manifest.json'),
        resources=resources,outputs={p.name:sha(p) for p in dest.iterdir() if p.is_file() and p.name not in ['manifest.json','COMPLETE']},
        source_sha256=sha(Path(__file__)),job=os.environ['SLURM_JOB_ID'],completed_at=utc())
    write_json(dest/'manifest.json',manifest);(dest/'COMPLETE').write_text(sha(dest/'manifest.json')+'\n')
    write_json(BASE/'PTC_DELIVERY.json',dict(status='complete',final_verification_sha256=sha(dest/'manifest.json'),
        source_replay_sha256=sha(RECOVERY/'R_baseline_replay/replay.json'),**{k:v for k,v in manifest.items() if k.endswith('_sha256')},
        terminal_conditions=11696,GBM_priority_respected=True,independent_Pu_cohort=False,completed_at=utc()))
    write_json(g/'EXPERIMENT_DELIVERY.json',dict(status='complete',GBM_delivery_sha256=sha(g/'GBM_DELIVERY.json'),
        PTC_delivery_sha256=sha(BASE/'PTC_DELIVERY.json'),notebook=str(nm['workspace_notebook']),notebook_sha256=nm['notebook_sha256'],
        decision_tree=str(BASE/'figures/decision_tree_complete.pdf'),decision_tree_sha256=sha(BASE/'figures/decision_tree_complete.pdf'),
        GBM_conditions_accounted=39930,PTC_terminal_conditions=11696,completed_at=utc()))
    print(json.dumps(manifest,indent=2),flush=True)

if __name__=='__main__':run()
