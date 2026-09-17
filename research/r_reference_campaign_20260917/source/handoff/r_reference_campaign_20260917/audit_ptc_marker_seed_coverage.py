"""Audit original selected PTC libraries: retained T markers and available DL classes."""
import hashlib,json,os,sys
from pathlib import Path
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
OUT=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'
sys.path.insert(0,str(ROOT/'handoff/ptc_recovery_20260916'))
from label_rules import strict_T

def run():
    assert os.environ.get('SLURM_JOB_ID')
    import numpy as np
    import pandas as pd
    libraries=json.loads((OUT/'markers/PTC_original17.json').read_text())
    reference=pd.read_csv(ROOT/'results/hvg_ptc_20260916_v1/ptc_paper_baseline/paper_baseline_reference.csv.gz',keep_default_na=False).set_index('cell_id')
    records=[];panels=[];missing=[]
    for prep in [OUT/'PTC_archived_CCA2000',*sorted((OUT/'PTC_ablation').iterdir())]:
        if not prep.is_dir():continue
        groups=['NMT','TTU'] if prep.name in ['PTC_archived_CCA2000','PTC_ALL8_CCA2000'] else [prep.name.split('_')[1]]
        for group in groups:
            library,route,cut=('CellMarker_Thyroid','PCA30_SNN','none') if group=='NMT' else ('Pubmed_34663816','UMAP2_HDBSCAN_R','mean')
            if prep.name=='PTC_archived_CCA2000':route='seurat_clusters' if group=='NMT' else 'hdbscan.UMAP_clusters'
            source=prep/route
            if not (source/'SCORE_COMPLETE').exists():missing.append([prep.name,group]);continue
            m=json.loads((source/'score_manifest.json').read_text())
            aid=next(a for a,v in m['arms'].items() if v['library']==library and v['cutoff']==cut)
            terminal=source/'terminal'/aid
            if not (terminal/'TERMINAL_COMPLETE').exists():missing.append([prep.name,group]);continue
            ids=pd.read_csv(source/'cells.csv',dtype=str).cell_id
            mask=reference.loc[ids,'group'].eq(group).to_numpy()
            retention=source/'marker_retention.csv.gz'
            if not retention.exists():retention=source/'marker_retention.csv'
            retained=pd.read_csv(retention,keep_default_na=False)
            tp={name:genes for name,genes in libraries[library].items() if strict_T(name)}
            tt=retained[retained.library.eq(library)&retained.panel.isin(tp)]
            assert set(tt.panel)==set(tp)
            for row in tt.itertuples(index=False):
                panels.append(dict(unit=prep.name,scope=group,library=library,panel=row.panel,
                    denominator=int(row.denominator),retained=int(row.retained),genes=';'.join(tp[row.panel])))
            with np.load(terminal/'terminal.npz',allow_pickle=False) as z:
                initial=np.asarray([strict_T(str(v)) for v in z['initial']],dtype=bool)
                final=np.asarray([strict_T(str(v)) for v in z['final090']],dtype=bool)
            if not initial.any():assert not final.any(),'DL predicted a T class absent from its input class vocabulary'
            records.append(dict(unit=prep.name,scope=group,route=route,library=library,cutoff=cut,
                n_cells_fit=len(ids),n_cells_scope=int(mask.sum()),n_initial_T_fit=int(initial.sum()),
                n_initial_T_scope=int(initial[mask].sum()),n_final_T_fit=int(final.sum()),n_final_T_scope=int(final[mask].sum()),
                n_T_panels=len(tt),n_T_panels_with_retained_genes=int(tt.retained.gt(0).sum()),
                T_class_absent_from_fit_vocabulary=not bool(initial.any())))
    dest=OUT/'verification';dest.mkdir(exist_ok=True)
    pd.DataFrame(records).to_csv(dest/'PTC_selected_marker_seed_coverage.csv',index=False)
    pd.DataFrame(panels).to_csv(dest/'PTC_selected_T_marker_retention.csv',index=False)
    report=dict(status='complete' if len(records)==32 and not missing else 'in_progress',job=os.environ['SLURM_JOB_ID'],
        selected_group_routes=len(records),expected_group_routes=32,missing=missing,
        zero_T_seed_fit_conditions=[r['unit'] for r in records if r['T_class_absent_from_fit_vocabulary']],
        verified='A fitted classifier never creates a T class absent from its initial class vocabulary.',
        interpretation='Marker loss, competitive density scoring and cutoff effects must be distinguished; retained markers alone do not guarantee a winning T seed.',
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    (dest/'PTC_marker_seed_coverage_audit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
    if '--require-complete' in sys.argv:assert report['status']=='complete',report

if __name__=='__main__':run()
