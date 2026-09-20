"""Exact finite native-R GBM core terminal/cache census; no recursive result scan."""
from pathlib import Path
from collections import Counter
import csv,hashlib,json,os,sys,time
sys.path.insert(0,str(Path(__file__).resolve().parent))
from collect_resources import BASE,OUT,sha,write_json,save_csv,utc
OLD=BASE/'paper_claim_validation_20260917'
assert os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_STEP_ID')
source=OLD/'protocol/core_tasks.json';tasks=json.loads(source.read_text());assert len(tasks)==726
routes=['PCA30_SNN','PCA30_HDBSCAN_R','UMAP2_SNN','UMAP2_HDBSCAN_R']
dest=OUT/'GBM_core_cache_census';dest.mkdir(exist_ok=True);cache={};statuses=Counter();actions=Counter();conditions=0
fields=['sample','budget','route','arm_id','library','cutoff','cache_key','identical_result_reused','training_executed','dl_status','job','training_manifest_sha256','source_manifest']
with (dest/'terminal_conditions.csv').open('w') as f:
 writer=csv.DictWriter(f,fields);writer.writeheader()
 for i,t in enumerate(tasks):
  for route in routes:
   p=OLD/'GBM'/t['sample']/t['budget']/route
   sm=json.loads((p/'score_manifest.json').read_text());assert len(sm['arms'])==48
   for aid,arm in sm['arms'].items():
    td=p/'terminal'/aid;mp=td/'terminal_manifest.json';raw=mp.read_bytes();assert (td/'TERMINAL_COMPLETE').read_text().strip()==hashlib.sha256(raw).hexdigest()
    m=json.loads(raw);assert m['terminal_valid'];key=m['cache_key'];assert Path(m['cache_directory']).name==key
    statuses[m['dl_status']]+=1
    action='cached_terminal_reuse' if m['identical_result_reused'] else ('fresh_training' if m['training_executed'] else 'fresh_no_training_terminal')
    actions[action]+=1;conditions+=1
    r=dict(sample=t['sample'],budget=t['budget'],route=route,arm_id=aid,library=arm['library'],cutoff=arm['cutoff'],source_manifest=str(mp))
    r.update({k:m[k] for k in ['cache_key','identical_result_reused','training_executed','dl_status','job','training_manifest_sha256']});writer.writerow(r)
    if key not in cache:
     cp=Path(m['cache_directory']);tp=cp/'training_manifest.json';data=tp.read_bytes();digest=hashlib.sha256(data).hexdigest();assert digest==m['training_manifest_sha256'];assert (cp/'COMPLETE').read_text().strip()==digest
     tr=json.loads(data);cache[key]=dict(cache_key=key,cache_directory=str(cp),training_executed=tr['training_executed'],dl_status=tr['dl_status'],elapsed_seconds=tr['elapsed_seconds'],producer_job=tr['job'],first_condition=tr.get('provenance',{}).get('first_condition',''),training_manifest_sha256=digest,n_requested_conditions=0,n_fresh_training_claims=0)
    assert cache[key]['training_manifest_sha256']==m['training_manifest_sha256']
    cache[key]['n_requested_conditions']+=1
    cache[key]['n_fresh_training_claims']+=int(action=='fresh_training')
  if i%50==0:print('CACHE_CENSUS',i+1,'/726',conditions,'conditions',len(cache),'unique cache entries',flush=True)
assert conditions==726*4*48==139392
save_csv(dest/'unique_cache_entries.csv',cache.values())
# Only this explicitly known collision directory is enumerated, not the cache tree.
dups=[];dd=OLD/'DL_cache/duplicate_publications'
if dd.exists():
 for p in dd.iterdir():
  proof=p/'DUPLICATE_TERMINAL_EQUIVALENCE.json'
  if not proof.is_file():continue
  j=json.loads(proof.read_text());published=Path(j['published_cache']).name
  if published not in cache:continue
  tp=p/'training_manifest.json';tr=json.loads(tp.read_text())
  dups.append(dict(cache_key=published,duplicate_directory=str(p),training_executed=tr['training_executed'],dl_status=tr['dl_status'],producer_job=tr['job'],elapsed_seconds=tr['elapsed_seconds'],proof_sha256=sha(proof)))
save_csv(dest/'preserved_duplicate_training_attempts.csv',dups,['cache_key','duplicate_directory','training_executed','dl_status','producer_job','elapsed_seconds','proof_sha256'])
trained=[r for r in cache.values() if r['training_executed']]
original=dict((r['dl_status'],int(r['n'])) for r in csv.DictReader((OLD/'summary/terminal_status_counts.csv').open()));assert dict(statuses)==original
m=dict(status='verified_core_condition_and_cache_census',scope='Exact726x4x48 original native-R GBM core; excludes geometry controls, MLP/representation controls and public/PTC',n_requested_conditions=conditions,condition_status_counts=dict(statuses),condition_action_counts=dict(actions),n_unique_cache_entries=len(cache),n_unique_trained_models=len(trained),n_unique_no_training_entries=len(cache)-len(trained),canonical_training_wall_seconds_sum=sum(r['elapsed_seconds'] for r in trained),n_preserved_duplicate_entries=len(dups),n_preserved_duplicate_actual_fits=sum(r['training_executed'] for r in dups),preserved_duplicate_training_wall_seconds_sum=sum(r['elapsed_seconds'] for r in dups if r['training_executed']),physical_fit_count_scope='Canonical cache fits plus explicitly preserved verified duplicates only; unregistered failed/partial attempts remain represented by accounting, not inferred here.',orphan_partial_cache_attempts_censused=False,source_tasks_sha256=sha(source),original_status_summary_exact=True,source_sha256=sha(__file__),job_step=os.environ['SLURM_JOB_ID']+'.'+os.environ['SLURM_STEP_ID'],completed_at=utc(),files={p.name:sha(p) for p in dest.iterdir() if p.is_file()})
write_json(dest/'manifest.json',m);print(json.dumps({k:v for k,v in m.items() if k!='files'}),flush=True)
