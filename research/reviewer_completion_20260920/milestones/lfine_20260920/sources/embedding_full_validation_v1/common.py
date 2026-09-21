"""Metadata-only IO; scientific entry points require SLURM explicitly."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os
from functools import lru_cache
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
CAMP=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
REFERENCE=ROOT/'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917'
EMBED=CAMP/'embedding';OUT=CAMP/'embedding_full_validation_v1';CODE=Path(__file__).resolve().parent
L1=['Malignant','TAM','Lymphocyte','Oligodendrocyte','Astrocyte','OPC','Excitatory neuron','Inhibitory neuron','Endothel','Pericyte','Other']
SPACES=['noDR','PCA2','FA2','ICA2','Isomap2','UMAP2','TSNE2'];KS=[5,10,15,20,30,40]
ROUTES=[f'{method}_K{k:02d}' for method in ['KMeans','GMM'] for k in KS]+['HDBSCAN_R']
STAGES={'marker_only':'initial','terminal090':'final090','terminal070':'final070'}
PILOT_SAMPLES=['TKU4163','NL022','SN040']
MEASURES=['macroF1_present','macroF1_fixed11','weightedF1','accuracy','coverage','unknown_rate','mapped_coverage','off_vocabulary_rate','coarse10_macroF1_present']
POLICY=CAMP/'protocol/embedding_convergence_repair_20260920_v2.json'
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def js(p):return json.loads(Path(p).read_text())
def utc():return datetime.now(timezone.utc).isoformat()
def write(p,value):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
    tmp=p.with_name(p.name+f'.part.{os.getpid()}');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def checked(p,name='manifest.json',flag='COMPLETE'):
    p=Path(p);return (p/name).is_file() and (p/flag).is_file() and (p/flag).read_text().strip()==sha(p/name)
def complete(p):
    p=Path(p);(p/'COMPLETE').write_text(sha(p/'manifest.json')+'\n')
def require_slurm():
    assert os.environ.get('SLURM_JOB_ID'),'All scientific validation must use SLURM'
    assert os.environ.get('PYTHONOPTIMIZE','0')=='0'
def taskkey(t):return tuple(t[k] for k in ['sample','budget','space'])
def validate_tasks(tasks):
    roster={taskkey(t) for t in js(OUT/'tasks.json')}
    assert tasks and all(set(t)=={'sample','budget','space'} for t in tasks)
    keys=[taskkey(t) for t in tasks]
    assert len(keys)==len(set(keys)) and set(keys)<=roster
def pilot_tasks():
    return [dict(sample=s,budget='hvg2000',space='PCA2') for s in PILOT_SAMPLES]+[dict(sample='TKU4163',budget='hvg2000',space='ICA2')]
def pilot_accepted():
    proof=js(OUT/'PILOT_ACCEPTANCE.json')
    assert proof['status']=='passed' and proof['validation_contract_sha256']==sha(OUT/'contract.json')
    expected={str(reportpath(t)/'manifest.json') for t in pilot_tasks()}
    assert set(proof['proofs'])==expected
    for task in pilot_tasks():
        path=reportpath(task)/'manifest.json'
        assert valid_report(task) and sha(path)==proof['proofs'][str(path)]
    return proof
def contract():
    c=js(OUT/'contract.json');assert c['status']=='prepared_for_gated_validation'
    for path,digest in c['source_hashes'].items():assert sha(path)==digest,path
    assert sha(OUT/'tasks.json')==c['tasks_sha256']
    return c
def unitpath(t):return EMBED/t['sample']/t['budget']/t['space']
def reportpath(t):return OUT/'units'/t['sample']/t['budget']/t['space']
def ready(t):
    p=unitpath(t)
    return checked(p,'fit_manifest.json','FIT_COMPLETE') and checked(p/'evaluation') and checked(p/'figures') and (t['space']!='ICA2' or checked(p/'figures_adaptive_v1'))
def signature(t):
    p=unitpath(t);names=['config.json','representation.json','fit_manifest.json','evaluation/manifest.json','figures/manifest.json']
    if t['space']=='ICA2':names.append('figures_adaptive_v1/manifest.json')
    return {name:sha(p/name) for name in names}
@lru_cache(maxsize=None)
def _scientific_artifact_sha(path,size,mtime_ns):return sha(path)
def valid_report(t,deep=False):
    d=reportpath(t)
    if not checked(d):return False
    m=js(d/'manifest.json')
    assert m['status']=='passed' and m['input_manifests']==signature(t),'Existing proof has stale inputs; preserve and create a new validation version'
    assert m['validation_contract_sha256']==sha(OUT/'contract.json')
    for name,digest in m['outputs'].items():assert sha(d/name)==digest,name
    if deep:
        require_slurm()
        artifacts=js(d/'source_artifacts.json');assert artifacts
        for name,digest in artifacts.items():
            path=Path(name);stat=path.stat()
            assert _scientific_artifact_sha(str(path),stat.st_size,stat.st_mtime_ns)==digest,name
    return True
