"""Independent CSV/manifest reconciliation of the bounded GBM cache census."""
from pathlib import Path
from collections import Counter
import csv
import json
import os
from collect_resources import OUT, sha, write_json, utc

assert os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_STEP_ID')
dest = OUT/'GBM_core_cache_census'
m = json.loads((dest/'manifest.json').read_text())
for name, expected in m['files'].items():
    assert sha(dest/name) == expected, name
cache = {}
with (dest/'unique_cache_entries.csv').open() as f:
    for r in csv.DictReader(f):
        assert r['cache_key'] not in cache
        cache[r['cache_key']] = r
assert len(cache) == m['n_unique_cache_entries']
statuses = Counter(); actions = Counter(); references = Counter(); fresh = Counter()
n = 0
with (dest/'terminal_conditions.csv').open() as f:
    for r in csv.DictReader(f):
        n += 1
        statuses[r['dl_status']] += 1
        action = 'cached_terminal_reuse' if r['identical_result_reused'] == 'True' else ('fresh_training' if r['training_executed'] == 'True' else 'fresh_no_training_terminal')
        actions[action] += 1
        references[r['cache_key']] += 1
        fresh[r['cache_key']] += action == 'fresh_training'
        assert cache[r['cache_key']]['training_manifest_sha256'] == r['training_manifest_sha256']
assert n == m['n_requested_conditions'] == 139392
assert dict(statuses) == m['condition_status_counts']
assert dict(actions) == m['condition_action_counts']
for key, r in cache.items():
    assert references[key] == int(r['n_requested_conditions'])
    assert fresh[key] == int(r['n_fresh_training_claims'])
trained = [r for r in cache.values() if r['training_executed'] == 'True']
assert len(trained) == m['n_unique_trained_models']
assert abs(sum(float(r['elapsed_seconds']) for r in trained)-m['canonical_training_wall_seconds_sum']) < 1e-5
with (dest/'preserved_duplicate_training_attempts.csv').open() as f:
    dup = list(csv.DictReader(f))
assert len(dup) == m['n_preserved_duplicate_entries']
assert sum(r['training_executed'] == 'True' for r in dup) == m['n_preserved_duplicate_actual_fits']
v = dict(status='passed', completed_at=utc(), census_manifest_sha256=sha(dest/'manifest.json'),
         source_sha256=sha(__file__), job_step=os.environ['SLURM_JOB_ID']+'.'+os.environ['SLURM_STEP_ID'],
         n_conditions=n, n_unique_cache_entries=len(cache), n_unique_actual_fits=len(trained),
         preserved_duplicate_actual_fits=m['n_preserved_duplicate_actual_fits'],
         all_terminal_statuses_reconciled=True, all_cache_references_reconciled=True,
         scope=m['scope'], resource_workpackage_complete=False)
write_json(dest/'verification.json', v)
status = json.loads((OUT/'status.json').read_text())
status['status'] = 'partial'
status['updated_at'] = utc()
status['completed'] = 'inventory_first_ledger_and_exact_GBM_core_cache_census'
status['jobs'] = list(dict.fromkeys(status['jobs'] + [m['job_step'], v['job_step']]))
status['remaining'] = [
    'Public/PTC and non-core GBM cache uniqueness/fresh-fit coverage is not fully censused; the verified GBM core census does not imply global closure.'
    if 'dedicated cache inventory remains outstanding' in text else text for text in status['remaining']]
status['evidence'] = list(dict.fromkeys(status['evidence'] + [str(dest/'manifest.json'), str(dest/'verification.json')]))
write_json(OUT/'status.json', status)
print(json.dumps(v), flush=True)
