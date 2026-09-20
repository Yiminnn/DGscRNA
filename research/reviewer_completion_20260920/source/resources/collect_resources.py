"""Bounded, explicit-registry-only own-user SLURM cost inventory; no matrix/model access."""
from pathlib import Path
from collections import defaultdict,Counter
from datetime import datetime,timezone
import csv,hashlib,json,os,pwd,re,subprocess,sys,time
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna');BASE=ROOT/'results/hvg_ptc_20260916_v1'
HERE=Path(__file__).resolve().parent;OUT=BASE/'reviewer_completion_20260920/resources'
PAT=re.compile(r'^[1-9][0-9]{5,8}(?:_[0-9]+)?(?:\.[0-9]+)?$')
FIELDS='JobID,JobIDRaw,DBIndex,SLUID,Restarts,JobName,User,State,ExitCode,Submit,Start,End,ElapsedRaw,AllocCPUS,CPUTimeRAW,TotalCPU,UserCPU,SystemCPU,ReqMem,MaxRSS,MaxRSSNode,MaxRSSTask,Partition,NodeList'

def utc():return datetime.now(timezone.utc).isoformat()
def sha(p):
 with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def write_json(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);t=p.with_name(p.name+'.tmp');t.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n');t.replace(p)
def save_csv(p,rows,fields=None):
 rows=list(rows);fields=fields or list(rows[0] if rows else {})
 with Path(p).open('w') as f:
  w=csv.DictWriter(f,fields);w.writeheader();w.writerows(rows)
def stage(text):
 t=text.lower()
 if any(s in t for s in ['scalab','scaling','cold_run']):return 'single_workflow_scaling_separate'
 if any(s in t for s in ['dispatch','controller','launcher','monitor']):return 'orchestration'
 if any(s in t for s in ['prepare','export_geometry','preprocess']):return 'prepare_integration_geometry'
 if any(s in t for s in ['representation_controls','run_controls','ptc_control_job','pipeline.py','full_array']):return 'combined_geometry_DEG_marker_DL'
 if any(s in t for s in ['scina','sccatch','singler','celltypist','deepsort','comparator']):return 'comparator'
 if any(s in t for s in ['terminal','refine','learning','mlp','dl_control','geometry_only']):return 'DL_refinement'
 if any(s in t for s in ['score','cluster','representation']):return 'clustering_DEG_marker_scoring'
 if any(s in t for s in ['aggregate','finaliz','deliver','report','plot','audit','verify','evaluat','evidence','resource']):return 'evaluation_figures_validation_delivery'
 return 'unresolved_or_combined'

def registry_paths():
 refs=[]
 for family,folder in [('native_R_public_PTC','r_reference_campaign_20260917'),('native_R_GBM_followups','paper_claim_validation_20260917')]:
  p=BASE/folder
  for f in p.iterdir():
   if f.is_file() and f.suffix in ['.json','.jsonl'] and any(v in f.name for v in ['dispatch','submission','controller','resource_resubmissions','allocation_resources']):refs.append((family,f))
 for family,folder in [('native_R_GBM_followups','paper_claim_validation_20260917/PTC_followups'),('PTC_earlier_reconciliation','ptc_experiments/operational')]:
  p=BASE/folder
  for f in p.iterdir():
   if f.is_file() and f.suffix=='.json':refs.append((family,f))
 old=BASE/'r_reference_campaign_20260917/resources/job_origin_index.json';refs.append(('native_R_public_PTC',old))
 for part in ['embedding','no_clustering','controls','evidence','comparison','manuscript','protocol']:
  p=BASE/'reviewer_completion_20260920'/part
  for name in ['jobs.json','status.json','manifest.json','verification.json','parity.json','resource_jobs.json']:
   if (p/name).is_file():refs.append(('reviewer_20260920_'+part,p/name))
 for folder in ['paper_claim_validation_20260917/resources','paper_claim_validation_20260917/summary','paper_claim_validation_20260917/GBM_full_summary','paper_claim_validation_20260917/PTC_summary','paper_claim_validation_20260917/controls_summary','paper_claim_validation_20260917/scalability_summary','r_reference_campaign_20260917/summary']:
  p=BASE/folder
  for name in ['manifest.json','aggregate_manifest.json','campaign_summary.json','accounting_manifest.json']:
   if (p/name).exists():refs.append(('native_R_public_PTC' if folder.startswith('r_reference') else 'native_R_GBM_followups',p/name))
 p=BASE/'reviewer_completion_20260920/embedding/dispatch/state.json'
 if p.exists():refs.append(('reviewer_20260920_embedding',p))
 p=OUT/'additional_explicit_jobs.json'
 if p.exists():refs.append(('reviewer_20260920_orchestration',p))
 return list(dict.fromkeys(refs))

def inventory(dest):
 prov=[];files=[]
 def extract(x,loc,context=''):
  found=[]
  if isinstance(x,dict):
   local=' '.join(str(x.get(k,'')) for k in ['purpose','command','name','stage','script'])
   for k,v in x.items():
    if PAT.fullmatch(str(k)):found.append((str(k),loc+'/'+str(k),local or context))
    if (('job' in k.lower() or k in ['previous_master','main_array','array_job']) and isinstance(v,(str,int)) and PAT.fullmatch(str(v))):found.append((str(v),loc+'/'+k,local or context))
    if ('job' in k.lower() or k in ['previous_master','main_array']) and isinstance(v,list):
     found.extend((str(a),loc+'/'+k,local or context) for a in v if isinstance(a,(str,int)) and PAT.fullmatch(str(a)))
    found.extend(extract(v,loc+'/'+str(k),local or context))
  elif isinstance(x,list):
   for i,v in enumerate(x):found.extend(extract(v,loc+'/'+str(i),context))
  return found
 for family,p in registry_paths():
  rows=[]
  if p.name=='job_origin_index.json':
   doc=json.loads(p.read_text())
   for jid,origins in doc.items():
    if PAT.fullmatch(jid):rows.append((jid,jid,'; '.join(origins[:3])))
  elif p.suffix=='.jsonl':
   for i,line in enumerate(p.open()):
    if line.strip():rows.extend(extract(json.loads(line),f'line:{i+1}'))
  elif p.name=='state.json' and p.parent.name=='dispatch':rows=extract(json.loads(p.read_text()).get('waves',{}),'root/waves')
  else:rows=extract(json.loads(p.read_text()),'root')
  files.append(dict(family=family,path=str(p),sha256=sha(p),bytes=p.stat().st_size,n_references=len(rows),n_unique_job_ids=len({r[0] for r in rows})))
  for jid,loc,context in rows:prov.append(dict(job_id=jid,family=family,stage=stage(context+' '+p.name),source=str(p),locator=loc,context=context))
 # Self-attribution is to this tiny step only, never the entire interactive allocation.
 self_step=os.environ['SLURM_JOB_ID']+'.'+os.environ['SLURM_STEP_ID']
 prov.append(dict(job_id=self_step,family='reviewer_20260920_resources',stage='evaluation_figures_validation_delivery',source=str(HERE/'collect_resources.py'),locator='SLURM_JOB_ID.SLURM_STEP_ID',context='This metadata-only resource collection step'))
 save_csv(dest/'registry_inventory.csv',files);save_csv(dest/'job_provenance.csv',prov)
 return prov,files

def seconds(s):
 if not s or s in ['Unknown','N/A']:return None
 days=0
 if '-' in s:a,s=s.split('-',1);days=int(a)
 vals=[float(v) for v in s.split(':')];n=0
 for v in vals:n=n*60+v
 return days*86400+n

def mem_kib(s):
 if not s:return None
 m=re.fullmatch(r'([0-9.]+)([KMGTPE]?)',s)
 if not m:return None
 value=float(m[1]);unit=m[2];return value*({'':1/1024,'K':1,'M':1024,'G':1024**2,'T':1024**3,'P':1024**4,'E':1024**5}[unit])

def collect(dest,prov):
 refs=defaultdict(list)
 for r in prov:refs[r['job_id']].append(r)
 query=set(refs)
 # Recorded prior aliases allow removing repeated array-element queries only when their recorded master is also queried.
 old=BASE/'r_reference_campaign_20260917/resources/slurm_job_ledger.csv'
 with old.open() as f:
  for row in csv.DictReader(f):
   visible=row['JobID'];raw=row['JobIDRaw'];base=visible.split('_')[0]
   if '_' in visible and base in query:
    if raw!=base:query.discard(raw)
    query.discard(visible)
 # Bare shared allocation is not a cost attribution unit. Specific recorded steps remain queryable.
 excluded=[]
 if '7204040' in query:query.remove('7204040');excluded.append('7204040')
 query=sorted(query);write_json(dest/'query_job_ids.json',query)
 print('RESOURCE_REGISTRY',len(refs),'references IDs;',len(query),'unique accounting queries',flush=True)
 raw={};queries=[]
 for i in range(0,len(query),180):
  ids=query[i:i+180];cmd=['sacct','--user=yimin','--jobs='+','.join(ids),'--duplicates','--array','--parsable2','--units=K','--starttime=2026-09-01','--format='+FIELDS]
  p=subprocess.run(cmd,text=True,capture_output=True,timeout=180)
  queries.append(dict(job_ids=ids,command=cmd,returncode=p.returncode,stderr=p.stderr,queried_at=utc()))
  if p.returncode:write_json(dest/'accounting_queries.json',queries);raise RuntimeError(p.stderr)
  for row in csv.DictReader(p.stdout.splitlines(),delimiter='|'):
   assert row['User'] in ['','yimin'],row['User']
   key=(row['DBIndex'],row['JobID'],row['Start'],row['End']);raw[key]=row
  if i//180%10==0:print('RESOURCE_QUERY_PROGRESS',i+len(ids),'/',len(query),'records',len(raw),flush=True)
 save_csv(dest/'accounting_raw.csv',raw.values(),FIELDS.split(','));write_json(dest/'accounting_queries.json',queries)
 save_csv(dest/'excluded_shared_allocations.csv',[dict(job_id=j,reason='Shared interactive allocation; only explicitly recorded steps can be attributed') for j in excluded],['job_id','reason'])
 return list(raw.values()),refs

def ledger(dest,records,refs):
 steps=defaultdict(list);jobs=[];seen_requested=set()
 for r in records:
  if '.' in r['JobID']:steps[(r['DBIndex'],r['JobID'].split('.')[0])].append(r)
 for r in records:
  jid=r['JobID'];raw=r['JobIDRaw'];base=jid.split('_')[0].split('.')[0]
  matches=[v for key in [jid,raw,base] for v in refs.get(key,[])]
  is_step='.' in jid
  if is_step and jid not in refs and raw not in refs:continue
  if not is_step and base=='7204040':continue
  if not matches:continue
  seen_requested.update(v['job_id'] for v in matches)
  # Explicit submission commands carry stronger stage attribution than generic summary job IDs.
  ranked=sorted(matches,key=lambda v:(0 if 'sbatch' in v['context'] else 1,0 if v['stage']!='unresolved_or_combined' else 1,v['source']))
  chosen=ranked[0];matching_steps=steps.get((r['DBIndex'],jid),[]) if not is_step else [r]
  batch=[s for s in matching_steps if s['JobID'].endswith('.batch')];rss=[(mem_kib(s['MaxRSS']),s['JobID']) for s in matching_steps if mem_kib(s['MaxRSS']) is not None]
  batchrss=[mem_kib(s['MaxRSS']) for s in batch if mem_kib(s['MaxRSS']) is not None]
  state=r['State'].split()[0];wall=float(r['ElapsedRaw'] or 0);cpus=int(r['AllocCPUS'] or 0);cpu=seconds(r['TotalCPU'])
  group=chosen['family'];phase=chosen['stage']
  if any(s in r['JobName'].lower() for s in ['scalab','scaling']):phase='single_workflow_scaling_separate'
  shared=r['JobName'] in ['claude-science','codex','jupyter'] and not is_step
  jobs.append(dict(job_id=jid,job_id_raw=raw,db_index=r['DBIndex'],sluid=r['SLUID'],restarts=r['Restarts'],family=group,phase=phase,job_name=r['JobName'],user='yimin',state=state,exit_code=r['ExitCode'],submit=r['Submit'],start=r['Start'],end=r['End'],elapsed_seconds=wall,alloc_cpus=cpus,allocated_cpu_hours=wall*cpus/3600,measured_cpu_hours=cpu/3600 if cpu is not None else '',measured_cpu_raw=r['TotalCPU'],user_cpu_raw=r['UserCPU'],system_cpu_raw=r['SystemCPU'],requested_memory=r['ReqMem'],batch_maxrss_kib=max(batchrss) if batchrss else '',max_step_maxrss_kib=max(rss)[0] if rss else '',maxrss_source_step=max(rss)[1] if rss else '',aggregate_fork_peak_memory='not_measured',gpu_peak_memory='N/A_CPU_implementation',accounting_final=state not in ['RUNNING','PENDING','COMPLETING','CONFIGURING'],failed_or_requeued=state in ['FAILED','OUT_OF_MEMORY','TIMEOUT','CANCELLED','NODE_FAIL','PREEMPTED','REQUEUED'],shared_allocation_excluded_from_totals=shared,attribution_source=chosen['source'],attribution_locator=chosen['locator'],attribution_context=chosen['context'],source_reference_count=len(matches)))
 save_csv(dest/'job_attempt_ledger.csv',jobs)
 missing=sorted(set(refs)-seen_requested-{'7204040'});write_json(dest/'missing_accounting_job_ids.json',missing)
 # A parent query can legitimately cover raw references; map exact raw IDs as observed too.
 observed_raw={r['JobIDRaw'] for r in records}|{r['JobID'] for r in records}
 missing=[v for v in missing if v not in observed_raw];write_json(dest/'missing_accounting_job_ids.json',missing)
 groups=defaultdict(list)
 for r in jobs:
  if not r['shared_allocation_excluded_from_totals']:groups[(r['family'],r['phase'])].append(r)
 summary=[]
 for (family,phase),rs in sorted(groups.items()):
  starts=[r['start'] for r in rs if r['start'] not in ['Unknown','None','']];ends=[r['end'] for r in rs if r['end'] not in ['Unknown','None','']]
  first=min(starts) if starts else '';last=max(ends) if ends else '';span=(datetime.fromisoformat(last)-datetime.fromisoformat(first)).total_seconds()/3600 if first and last else ''
  summary.append(dict(family=family,phase=phase,n_accounting_attempt_records=len(rs),n_final=sum(r['accounting_final'] for r in rs),n_live_or_pending=sum(not r['accounting_final'] for r in rs),n_failed_or_requeued=sum(r['failed_or_requeued'] for r in rs),allocated_cpu_hours=sum(r['allocated_cpu_hours'] for r in rs),measured_cpu_hours=sum(r['measured_cpu_hours'] for r in rs if r['measured_cpu_hours']!=''),sum_job_wall_hours=sum(r['elapsed_seconds'] for r in rs)/3600,first_start=first,last_known_end=last,calendar_span_hours=span,max_reported_step_rss_kib=max([r['max_step_maxrss_kib'] for r in rs if r['max_step_maxrss_kib']!=''],default=''),aggregate_peak_memory='unavailable',scope='Separate timing experiment' if phase=='single_workflow_scaling_separate' else 'Recorded search/analysis jobs; mixed-stage jobs unsplit'))
 save_csv(dest/'phase_cost_summary.csv',summary)
 return jobs,missing

def terminal_evidence(dest):
 evidence=[]
 def add(scope,path,counting_unit,counts,fresh=None,reuse=None,unique=None,overlap=''):
  evidence.append(dict(scope=scope,source=str(path),sha256=sha(path),counting_unit=counting_unit,status_counts=counts,fresh_training_count=fresh,cached_terminal_reuse_count=reuse,unique_physical_models=unique,overlap_warning=overlap))
 p=BASE/'paper_claim_validation_20260917/summary/terminal_status_counts.csv';rows=list(csv.DictReader(p.open()));add('native_R_GBM_core',p,'requested terminal condition; cached results retain original training flag',{r['dl_status']:int(r['n']) for r in rows})
 p=BASE/'r_reference_campaign_20260917/summary/all_terminal_statuses.csv';counts=Counter(r['dl_status'] for r in csv.DictReader(p.open()));add('native_R_public_and_PTC_core',p,'requested terminal condition; includes30 PTC units',dict(counts),overlap='PTC reused grid must not be counted as fresh fitting in later followups')
 p=BASE/'paper_claim_validation_20260917/PTC_summary/terminal_execution_overview.json';j=json.loads(p.read_text());add('PTC_new_followups_only',p,j['counting_unit'],j['dl_status_counts'],fresh=j['n_fresh_training'],reuse=j['n_cached_terminal_reuse'],overlap='Excludes reused original grid; copied flag training_executed alone is not fresh training')
 p=BASE/'paper_claim_validation_20260917/controls_summary/MLP_learning_training_status.csv';counts=Counter(r['dl_status'] for r in csv.DictReader(p.open()));add('GBM_old_MLP_controls',p,'requested condition; cached execution flag',dict(counts))
 write_json(dest/'terminal_execution_evidence.json',evidence)
 scaling=BASE/'paper_claim_validation_20260917/scalability_summary/manifest.json';j=json.loads(scaling.read_text());write_json(dest/'separate_scaling_evidence.json',dict(source=str(scaling),sha256=sha(scaling),n_cold_runs=j['n_cold_runs'],largest_single_run=j['largest_single_run'],timing_replicates=j['timing_replicates'],included_in_full_search_total=False,gpu_peak_memory='N/A_CPU_implementation'))
 return evidence

def main():
 assert os.environ.get('SLURM_JOB_ID') and os.environ.get('SLURM_STEP_ID'),'Metadata/scientific aggregation must run in explicit SLURM step'
 assert pwd.getpwuid(os.getuid()).pw_name=='yimin'
 OUT.mkdir(exist_ok=True,parents=True);stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ');dest=OUT/('snapshot_'+stamp);dest.mkdir()
 write_json(OUT/'status.json',dict(stage='F_resources',status='inventory_running',updated_at=utc(),jobs=[os.environ['SLURM_JOB_ID']+'.'+os.environ['SLURM_STEP_ID']],completed=0,remaining=['cost_coverage','live_campaign_final_accounting','unique_DL_cache_counts'],evidence=[str(dest)]))
 prov,files=inventory(dest);raw,refs=collect(dest,prov);jobs,missing=ledger(dest,raw,refs);terminal_evidence(dest)
 limitations=['Current reviewer campaign is still running: costs are an observed snapshot, not final total.','Stage timing cannot be split where a single job combined prepare/cluster/DEG/marker/DL and no separate timers were recorded.','Batch/step MaxRSS is not summed fork/process peak memory; aggregate peak is unavailable.','GBM/public inherited training flags do not prove fresh training or unique physical model counts; dedicated cache inventory remains outstanding.','Earlier hand-run/original-paper computations lacking explicit registries remain outside proven cost coverage.','Allocated CPU-hours use actual scheduler AllocCPUS, which can exceed the application threads because of memory allocation.','sum_job_wall_hours is the sum of allocation/step durations, not wall-clock turnaround; calendar_span is separate and includes gaps.','No CPU-only accounting record supplies GPU peak memory.45 cold-run scaling measurements are a separate experiment.']
 manifest=dict(status='partial_inventory_and_first_ledger',identity='yimin',scope='Explicit local registry job IDs only; no recursive result-tree traversal and no model/matrix loading',n_registry_files=len(files),n_job_provenance_rows=len(prov),n_unique_requested_job_ids=len(refs),n_accounting_rows=len(raw),n_accounting_attempt_records=len(jobs),n_missing_accounting_ids=len(missing),n_live_or_pending=sum(not r['accounting_final'] for r in jobs),limitations=limitations,job_step=os.environ['SLURM_JOB_ID']+'.'+os.environ['SLURM_STEP_ID'],source_sha256=sha(__file__),completed_at=utc(),files={p.name:sha(p) for p in dest.iterdir() if p.is_file()})
 write_json(dest/'manifest.json',manifest);write_json(OUT/'latest.json',dict(snapshot=str(dest),manifest_sha256=sha(dest/'manifest.json')))
 write_json(OUT/'status.json',dict(stage='F_resources',status='partial',updated_at=utc(),jobs=[manifest['job_step']],completed='inventory_and_first_ledger',remaining=limitations,evidence=[str(dest/'manifest.json'),str(dest/'phase_cost_summary.csv'),str(dest/'terminal_execution_evidence.json')]))
 print(json.dumps({k:v for k,v in manifest.items() if k not in ['files','limitations']}),flush=True)
if __name__=='__main__':main()
