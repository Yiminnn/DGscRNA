"""Independently compare completed outputs with frozen libraries and experimental design."""
import hashlib,itertools,json,os,sys
from pathlib import Path

ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
CODE=ROOT/'handoff/r_reference_campaign_20260917'
OUT=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'
ROUTES=('PCA30_SNN','PCA30_HDBSCAN_R','UMAP2_SNN','UMAP2_HDBSCAN_R')
ARCHIVED=('seurat_clusters','hdbscan_clusters','seurat.UMAP_clusters','hdbscan.UMAP_clusters')
CUTS=('none','mean','0.5')
STAGES=('initial','final090','final070')

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def expected_units():
    roster=json.loads((OUT/'markers/marker_roster.json').read_text())['units']
    names=[r['unit'] for r in roster]
    assert len(names)==len(set(names))==70
    assert sum(r['dataset']=='HCL' for r in roster)==59
    assert {r['dataset'] for r in roster}=={
        'brain_GBM','breast_TNBC','colorectal','kidney_ccRCC','blood_DLBCL',
        'baron_human','muraro','segerstolpe','xin','immune_ALL_human','HCL'}
    assert set(names)==set((OUT/'benchmark_units.txt').read_text().splitlines())
    units=[(r['unit'],r['dataset'],OUT/'benchmark'/r['unit']/'reference_CCA2000') for r in roster]
    combined=['PTC_ALL8_CCA2000']
    combined += [f'PTC_{g}_CCA{h}' for g in ('NMT','TTU') for h in ('500','1000','2000','3000','5000','all')]
    combined += [f'PTC_{g}_{method}_RNA2000' for g in ('NMT','TTU') for method in ('NONE','HARMONY')]
    assert set(combined)=={s['unit'] for s in json.loads((CODE/'ptc_ablation_conditions.json').read_text())}
    ptc=combined+[f'PTC_{g}_GEOMETRY{h}_FIXED_CCAall_DL2000' for g in ('NMT','TTU') for h in ('500','1000','2000','3000','5000','all')]
    units += [(u,'PTC',OUT/'PTC_ablation'/u) for u in ptc]
    units += [('PTC_archived_CCA2000','PTC',OUT/'PTC_archived_CCA2000')]
    return units,{r['unit']:r for r in roster}

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import pandas as pd
    units,roster=expected_units()
    reference=pd.read_csv(ROOT/'results/hvg_ptc_20260916_v1/ptc_paper_baseline/paper_baseline_reference.csv.gz',
        usecols=['group','sample'],dtype=str,keep_default_na=False).drop_duplicates()
    records=[];incomplete=[]
    for unit,dataset,prep in units:
        marker=OUT/'markers'/('PTC_original17.json' if dataset=='PTC' else unit+'.json')
        libraries=json.loads(marker.read_text())
        if dataset=='PTC':assert len(libraries)==17
        else:assert {k:len(v) for k,v in libraries.items()}==roster[unit]['libraries']
        arms={f'L{i:02d}_{"p050" if cut=="0.5" else cut}':(name,cut)
              for i,name in enumerate(libraries) for cut in CUTS}
        routes=ARCHIVED if unit=='PTC_archived_CCA2000' else ROUTES
        found={p.parent.name for p in prep.glob('*/score_manifest.json')}
        assert found<=set(routes),(unit,'unexpected routes',found-set(routes))
        scored=[];source_hashes={}
        for route in routes:
            source=prep/route
            if not (source/'SCORE_COMPLETE').exists():continue
            manifest=source/'score_manifest.json';m=json.loads(manifest.read_text())
            assert (source/'SCORE_COMPLETE').read_text().strip()==sha(manifest)
            actual={aid:(a['library'],str(a['cutoff'])) for aid,a in m['arms'].items()}
            assert actual==arms,(unit,route,'frozen library/cutoff identities differ')
            assert all(a['arm_id']==aid and a['seed_column']==aid for aid,a in m['arms'].items())
            if 'marker_source_sha256' in m:assert m['marker_source_sha256']==sha(marker)
            source_hashes[str(manifest.resolve())]=sha(manifest)
            cols=pd.read_csv(source/'initial_calls.csv.gz',nrows=0).columns.tolist()
            assert cols==['cell_id',*arms],(unit,route,'initial columns differ')
            scored.append(route)
        if len(scored)!=4 or not (prep/'evaluation/COMPLETE').exists():
            incomplete.append(unit);continue
        terminal={(r,aid) for r in routes for aid in arms}
        found_terminal={(p.parent.parent.parent.name,p.parent.name) for p in prep.glob('*/terminal/*/TERMINAL_COMPLETE')}
        assert found_terminal==terminal,(unit,'terminal grid differs')
        statuses=pd.read_csv(prep/'evaluation/terminal_statuses.csv',dtype=str,keep_default_na=False)
        assert len(statuses)==len(terminal) and set(zip(statuses.route,statuses.arm_id))==terminal
        cl=pd.read_csv(prep/'evaluation/clustering_metrics.csv',dtype=str,keep_default_na=False)
        assert len(cl)==4 and set(cl.route)==set(routes)
        metrics=pd.read_csv(prep/'evaluation/metrics.csv.gz',dtype=str,keep_default_na=False)
        if dataset=='PTC':
            refs=reference if unit in ('PTC_archived_CCA2000','PTC_ALL8_CCA2000') else reference[reference.group.eq('NMT' if '_NMT_' in unit else 'TTU')]
            scopes={'all',*refs.group,*refs['sample']}
            endpoints=('paper_broad_T_compatibility','strict_T_name_rule')
        else:
            scopes={''};metrics['scope']=''
            endpoints=('curated_semantic','common_lineage')
        expected=set(itertools.product(routes,arms,STAGES,scopes,endpoints))
        actual=list(metrics[['route','arm_id','stage','scope','endpoint']].itertuples(index=False,name=None))
        assert len(actual)==len(expected) and set(actual)==expected,(unit,'evaluation rows differ from frozen grid')
        assert set(metrics.unit)=={unit} and set(metrics.dataset)=={dataset}
        for row in metrics[['arm_id','library','cutoff']].drop_duplicates().itertuples(index=False):
            assert (row.library,row.cutoff)==arms[row.arm_id]
        em=prep/'evaluation/evaluation_manifest.json'
        assert (prep/'evaluation/COMPLETE').read_text().strip()==sha(em)
        evaluation_sources={str((ROOT/Path(p)).resolve()):h for p,h in json.loads(em.read_text())['sources'].items()}
        assert evaluation_sources==source_hashes,(unit,'evaluation sources differ')
        records.append(dict(unit=unit,dataset=dataset,status='passed',libraries=len(libraries),
            routes=list(routes),expected_terminal_conditions=len(terminal),expected_evaluation_rows=len(expected),
            marker_sha256=sha(marker),evaluation_manifest_sha256=sha(em)))
    report=dict(status='passed' if len(records)==len(units) else 'in_progress',job=os.environ['SLURM_JOB_ID'],
        expected_units=len(units),units=len(records),incomplete_units=incomplete,
        expected_grid_source='Frozen marker JSON × three cutoffs × four prescribed routes; three endpoints stages; explicit evaluation resolutions/scopes',
        records=records,script_sha256=sha(__file__))
    path=OUT/'verification/frozen_grid_audit.json';temp=path.with_suffix('.part')
    temp.write_text(json.dumps(report,indent=2)+'\n');temp.replace(path)
    print(json.dumps({k:v for k,v in report.items() if k!='records'},indent=2),flush=True)
    if '--require-complete' in sys.argv:assert report['status']=='passed','Frozen grid is not complete'

if __name__=='__main__':run()
