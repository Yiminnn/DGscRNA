"""Check real accepted packaged pilot metrics against old inferential inputs."""
from pathlib import Path
import os
import json
assert os.environ.get('SLURM_JOB_ID')
import run
from revalidate import compare_rows, MEASURES
import pandas as pd
freshcore = run.CAMP / 'evaluation/pilots/v3_TKU4163_hvg2000_full'
fresha2 = run.CAMP / 'pilots/no_cluster_TKU4163_hvg2000/evaluation'
inputs = {}
for directory in [freshcore, fresha2]:
    run.checked(directory, inputs, 'completed')
oldroot = run.OLD / 'no_clustering_lfine_v1/summary'
old = run.read_table(oldroot / 'all_candidate_metrics.csv', True)
old = old[old['sample'].eq('TKU4163') & old.budget.eq('hvg2000') & old.stage.isin(run.STAGES)]
fresh = run.read_table(fresha2 / 'metrics.csv.gz')
fresh['lambda_value'] = fresh.cutoff.str.removeprefix('cellwise_lambda_').astype(float)
_, a2 = compare_rows(old, fresh, ['sample','budget','lambda_value','stage'], ['lfine_macroF1','coverage'], exact=True)
old = run.read_table(oldroot / 'all_original_route_metrics.csv', True)
old = old[old['sample'].eq('TKU4163') & old.budget.eq('hvg2000') & old.stage.isin(run.STAGES)]
fresh = run.read_table(freshcore / 'metrics.csv.gz')
fresh = fresh[fresh.library.eq('CM2_glioma_other') & fresh.cutoff.eq('mean')]
_, anchors = compare_rows(old, fresh, ['sample','budget','route','stage'], ['lfine_macroF1','coverage'], exact=True)
old = run.read_table(run.OLD / 'comparison_lfine_v1/summary/all_eligible_candidate_metrics.csv.gz', True)
old = old[old['sample'].eq('TKU4163') & old.method.eq('DG-scRNA') & old.primary_candidate]
fresh = run.read_table(freshcore / 'metrics.csv.gz', True)
fresh = fresh[fresh.stage.eq('terminal090') & fresh.route.eq('UMAP2_HDBSCAN_R')].assign(method='DG-scRNA')
_, comp = compare_rows(old, fresh, ['sample','method','budget','route','library','cutoff'], MEASURES, exact=True)
result = dict(status='passed_fresh_pilot_metric_inputs_exact', sample='TKU4163',budget='hvg2000',
    A2_candidate_rows=10,A2_anchor_rows=8,C_DG_candidate_rows=39,
    max_abs_delta=dict(A2=a2,A2_anchors=anchors,C=comp),input_hashes=inputs,
    full_campaign_or_inference_claimed=False,job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'))
out = run.CAMP / 'final_report_validation/fresh_pilot_metrics.json'
run.write(out,result)
print(json.dumps(result,indent=2))
