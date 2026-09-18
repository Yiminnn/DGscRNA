"""Fix RNA scoring and normalized HVG2000 DL inputs while varying geometry."""
import json
import os
import sys
from pathlib import Path
from common import OUT, ROUTES, checked, complete, require_slurm, sha, write_json, utc

LIBRARIES=['CM2_glioma_other','CM2_primary_all_context']

def run(sample,budget):
    require_slurm()
    import torch
    import evaluate
    import terminal
    torch.set_num_threads(4)
    assert budget in ['hvg2000','hvg5000','all']
    prep=OUT/'GBM'/sample/budget
    fixed=OUT/'GBM'/sample/'hvg2000'
    dest=prep/'evaluation_geometry_only_DL2000'
    if checked(dest):return
    assert checked(fixed,'prepare_manifest.json','PREPARED')
    dm=json.loads((fixed/'prepare_manifest.json').read_text())
    armids=None
    for route in ROUTES:
        source=prep/route
        assert checked(source,'score_manifest.json','SCORE_COMPLETE')
        sm=json.loads((source/'score_manifest.json').read_text())
        selected=[key for key,a in sm['arms'].items() if a['library'] in LIBRARIES and a['cutoff']=='mean']
        assert len(selected)==2
        if armids is None:armids=selected
        assert selected==armids
        for aid in selected:terminal.finish_route(source,only_arm=aid,dl_prep=fixed)
    evaluate.run(prep,family='geometry_only_fixed_DL2000',only_arms=armids)
    write_json(dest/'design.json',dict(status='completed',sample=sample,geometry_budget=budget,
        fixed_DL_features=dm['features']['DL'],fixed_DL_sha256=dm['DL_binary_sha256'],
        scoring='Full RNA assay is already identical across the single-sample native-R budget arms',
        libraries=LIBRARIES,cutoff='mean',no_performance_selected_libraries=True,
        statement='Controlled geometry comparison; no new DEG calculation and no truth labels used for fitting.',
        job=os.environ['SLURM_JOB_ID'],source_sha256=sha(__file__),completed_at=utc()))

if __name__=='__main__':run(*sys.argv[1:])
