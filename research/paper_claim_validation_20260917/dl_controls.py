"""Bounded MLP-only seed/width/epoch sensitivity on the three prespecified pilots."""
import json
import os
from pathlib import Path
import shutil
import sys
from common import OUT, ROUTES, L1, checked, complete, require_slurm, sha, write_json, utc
CONFIGS={'original':dict(architecture=[256,128],epochs=10,model_seed=42,split_seed=42)}
CONFIGS.update({f'model_seed{s}':dict(architecture=[256,128],epochs=10,model_seed=s,split_seed=42) for s in [0,1,2,3]})
CONFIGS.update(width_small=dict(architecture=[128,64],epochs=10,model_seed=42,split_seed=42),
               width_large=dict(architecture=[512,256],epochs=10,model_seed=42,split_seed=42),
               epochs5=dict(architecture=[256,128],epochs=5,model_seed=42,split_seed=42),
               epochs20=dict(architecture=[256,128],epochs=20,model_seed=42,split_seed=42))
LIBRARIES=['CM2_glioma_other','CM2_primary_all_context']

def run(sample,budget,config,parity=False):
    require_slurm()
    import numpy as np
    import pandas as pd
    import torch
    from sklearn.metrics import precision_recall_fscore_support
    import terminal
    import refine_controls
    torch.set_num_threads(4)
    assert config in CONFIGS
    pilots=(OUT/'protocol/pilot_samples.txt').read_text().split()
    assert sample in pilots and budget in ['hvg2000','hvg5000']
    dest=OUT/'GBM_DL_controls'/sample/budget/config;dest.mkdir(parents=True,exist_ok=True)
    if checked(dest):return
    if not parity:
        proof=json.loads((OUT/'verification/DL_control_default_parity.json').read_text())
        assert proof['status']=='passed' and proof['refinement_source_sha256']==sha(refine_controls.__file__)
    refine_controls.PARAMS.update(CONFIGS[config]);terminal.refine=refine_controls
    verified=[]
    for route in ROUTES:
        original=OUT/'GBM'/sample/budget/route
        assert checked(original,'score_manifest.json','SCORE_COMPLETE')
        src=dest/route;src.mkdir(exist_ok=True)
        for name in ['score_manifest.json','SCORE_COMPLETE','cells.csv','initial_calls.csv.gz']:
            target=src/name
            if not target.exists():os.link(original/name,target)
            assert sha(target)==sha(original/name)
        sm=json.loads((src/'score_manifest.json').read_text())
        aids=[aid for aid,a in sm['arms'].items() if a['library'] in LIBRARIES and a['cutoff']=='mean']
        assert len(aids)==2
        for aid in aids:
            terminal.finish_route(src,only_arm=aid)
            if parity:
                ref=np.load(original/'terminal'/aid/'terminal.npz',allow_pickle=False)
                got=np.load(src/'terminal'/aid/'terminal.npz',allow_pickle=False)
                for key in ['initial','final090','final070','probabilities']:
                    assert np.array_equal(ref[key],got[key]),(route,aid,key)
                verified.append(route+'/'+aid)
    # Only after fitting has ended may evaluation labels be opened.
    truth=pd.read_csv(OUT/'evaluation_inputs'/sample/'truth.csv.gz',dtype=str,keep_default_na=False)
    y=truth.L1.to_numpy();mapping=pd.read_csv(OUT/'markers/panel_L1_mapping.csv',dtype=str,keep_default_na=False)
    maps={lib:dict(zip(g.panel,g.L1)) for lib,g in mapping.groupby('library')}
    rows=[]
    for route in ROUTES:
        sm=json.loads((dest/route/'score_manifest.json').read_text())
        for aid,a in sm['arms'].items():
            if a['library'] not in LIBRARIES or a['cutoff']!='mean':continue
            td=dest/route/'terminal'/aid;tm=json.loads((td/'terminal_manifest.json').read_text())
            pred=pd.read_csv(td/'predictions.csv.gz',dtype=str,keep_default_na=False)
            assert np.array_equal(pred.cell_id,truth.cell_id)
            for stage in ['final090','final070']:
                p=np.asarray([maps[a['library']].get(v,'Unknown' if v in ['Unknown','Undecided',''] else 'UNMAPPABLE') for v in pred[stage]])
                _,_,f1,support=precision_recall_fscore_support(y,p,labels=L1,zero_division=0)
                rows.append(dict(sample=sample,budget=budget,route=route,library=a['library'],cutoff='mean',
                    control=config,stage=stage,macroF1_present=float(f1[support>0].mean()),accuracy=float((y==p).mean()),
                    coverage=float((~pred[stage].isin(['Unknown','Undecided',''])).mean()),dl_status=tm['dl_status'],
                    model_seed=CONFIGS[config]['model_seed'],split_seed=42,epochs=CONFIGS[config]['epochs'],
                    architecture=str(CONFIGS[config]['architecture'])))
    pd.DataFrame(rows).to_csv(dest/'metrics.csv',index=False)
    write_json(dest/'manifest.json',dict(status='completed',sample=sample,budget=budget,control=config,
        params=CONFIGS[config],n_conditions=8,n_metric_rows=len(rows),
        scope='MLP-only control; fixed RNA features, clusters, marker seeds and90/10training split. This is not an embedding/whole-pipeline seed replicate.',
        files={'metrics.csv':sha(dest/'metrics.csv')},source_sha256=sha(__file__),
        refinement_source_sha256=sha(refine_controls.__file__),job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(dest)
    if parity:
        write_json(OUT/'verification/DL_control_default_parity.json',dict(status='passed',conditions=verified,
            exactly_equal=['initial','final090','final070','all class probabilities'],
            refinement_source_sha256=sha(refine_controls.__file__),job=os.environ['SLURM_JOB_ID']))

if __name__=='__main__':run(*sys.argv[1:4],parity='--parity' in sys.argv[4:])
