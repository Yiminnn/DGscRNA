"""Explain canonical-cache fits versus condition-local training claims."""
from collections import Counter
import csv
import json
import os
from collect_resources import OUT, sha, write_json, utc
assert os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_STEP_ID')
dest = OUT/'GBM_core_cache_census'
assert json.loads((dest/'verification.json').read_text())['status'] == 'passed'
with (dest/'unique_cache_entries.csv').open() as f:
    trained = [r for r in csv.DictReader(f) if r['training_executed'] == 'True']
claims = Counter(int(r['n_fresh_training_claims']) for r in trained)
v = dict(status='verified_counting_units', at=utc(), canonical_trained_models=len(trained),
    core_condition_fresh_training_claims=sum(int(r['n_fresh_training_claims']) for r in trained),
    canonical_entries_by_n_core_fresh_training_claims=dict(claims),
    no_core_fresh_claim_first_conditions=[{k:r[k] for k in ['cache_key','first_condition','producer_job','n_requested_conditions']} for r in trained if int(r['n_fresh_training_claims']) == 0],
    explanation='A canonical trained model can be reused by all core conditions when the cache was produced by an earlier condition. Actual unique models and current core condition-local training flags have different denominators. Preserved duplicate fits are reported separately.',
    source_sha256=sha(__file__), cache_csv_sha256=sha(dest/'unique_cache_entries.csv'),
    job_step=os.environ['SLURM_JOB_ID']+'.'+os.environ['SLURM_STEP_ID'])
write_json(dest/'counting_units.json',v)
print(json.dumps({k:value for k,value in v.items() if k!='no_core_fresh_claim_first_conditions'}))
