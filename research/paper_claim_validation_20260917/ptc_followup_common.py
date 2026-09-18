"""PTC follow-up paths and a strict GBM-before-PTC scientific execution gate."""
import json
from pathlib import Path
from common import ROOT, OUT, OLD, require_slurm, sha, checked, write_json, complete, utc

PTC = OUT/'PTC_followups'
REFERENCE = ROOT/'results/hvg_ptc_20260916_v1/ptc_experiments/evaluation_reference'
PAPER = ROOT/'results/hvg_ptc_20260916_v1/ptc_paper_baseline'
MARKERS = OLD/'markers/PTC_original17.json'
GROUPS = {'NMT':['MT-1','MT-2','N-1','N-2'], 'TTU':['TU-1','TU-2','T-1','T-2']}
PATIENTS = {'N-1':'Patient1','T-1':'Patient1','N-2':'Patient2','T-2':'Patient2',
            'MT-1':'Patient3','TU-1':'Patient3','MT-2':'Patient4','TU-2':'Patient4'}
ANCHORS = {'NMT':dict(library='CellMarker_Thyroid',cutoff='none',route='PCA30_SNN',archive='seurat_clusters'),
           'TTU':dict(library='Pubmed_34663816',cutoff='mean',route='UMAP2_HDBSCAN_R',archive='hdbscan.UMAP_clusters')}
ROUTE_MAP = dict(seurat_clusters='PCA30_SNN',hdbscan_clusters='PCA30_HDBSCAN_R',
                 **{'seurat.UMAP_clusters':'UMAP2_SNN','hdbscan.UMAP_clusters':'UMAP2_HDBSCAN_R'})
SEEDS = [0,1,2,3,42]
CONTROL_ROUTES = ['PCA30_SNN','UMAP2_HDBSCAN_R']

def require_ptc():
    require_slurm()
    gate=OUT/'GBM_full_summary/GBM_FULL_DELIVERED.json'
    assert gate.exists(), 'New PTC science must wait for completed, verified GBM delivery'
    record=json.loads(gate.read_text())
    assert record['PTC_compute_may_start'] is True
    assert record['receipt_sha256']==sha(OUT/'GBM_full_summary/DELIVERY_RECEIPT.json')

def old_units():
    units=[OLD/'PTC_archived_CCA2000']+sorted((OLD/'PTC_ablation').glob('PTC_*'))
    units=[p for p in units if p.is_dir() and (p/'evaluation/COMPLETE').exists()]
    assert len(units)==30, len(units)
    return units

def original_arm(source,group):
    arms=json.loads((Path(source)/'score_manifest.json').read_text())['arms']
    a=ANCHORS[group]
    found=[key for key,v in arms.items() if v['library']==a['library'] and v['cutoff']==a['cutoff']]
    assert len(found)==1
    return found[0]

def read_reference():
    require_ptc()
    import numpy as np
    import pandas as pd
    assert checked(REFERENCE)
    ref=pd.read_csv(REFERENCE/'reference_cells.csv.gz',keep_default_na=False).set_index('cell_id')
    paper=pd.read_csv(PAPER/'original_DG_binary_pairs_for_R.csv.gz',keep_default_na=False).set_index('cell_id')
    native=pd.read_csv(PAPER/'paper_baseline_reference.csv.gz',keep_default_na=False).set_index('cell_id')
    assert ref.index.is_unique and len(ref)==92404 and set(ref.index)==set(paper.index)==set(native.index)
    expected=ref['sample'].map(PATIENTS)
    assert expected.notna().all()
    assert set(ref['sample'])==set(PATIENTS)
    # Verify the archived patient relation even if its labels use another spelling.
    assert ref.groupby('patient')['sample'].apply(lambda s:len(set(s.map(PATIENTS)))).eq(1).all()
    assert expected.groupby(ref.patient).nunique().eq(1).all() and ref.patient.nunique()==4
    ref['patient']=expected
    ref['group']=ref['sample'].map({s:g for g,ss in GROUPS.items() for s in ss})
    assert ref.groupby('group').size().to_dict()=={'NMT':48255,'TTU':44149}
    ref['paper_truth']=paper.loc[ref.index,'truth'].astype(int)
    ref['paper_native']=native.loc[ref.index,'paper_final_native']
    for field in ['TCR_cell_high_confidence_productive_TCR','TCR_any_filtered_contig','TCR_S3_supplied']:
        values=ref[field].astype(str).str.lower()
        assert values.isin(['true','false','1','0']).all()
        ref[field]=values.isin(['true','1'])
    assert set(np.unique(ref.paper_truth))=={0,1}
    return ref
