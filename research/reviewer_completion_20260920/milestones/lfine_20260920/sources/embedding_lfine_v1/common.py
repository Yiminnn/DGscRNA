"""A1 Lfine addendum: immutable source identity and explicit archival acceptance gate."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,json,os
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
CAMP=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
NATIVE=ROOT/'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917'
ARCHIVE=CAMP/'embedding_full_validation_v1';EMBED=CAMP/'embedding'
OUT=CAMP/'embedding_lfine_v1';CODE=Path(__file__).resolve().parent
REFERENCE=CODE/'source_reference';FIXED='CM2_glioma_other'
STAGES={'marker_only':'initial','terminal090':'final090','terminal070':'final070'}
FIELDS=['lfine_macroF1','coverage','legacy_called_coverage','abstain_rate','offvocab_rate','acc_on_called',
    'n_distinct_calls','n_classes_hit','lfine_n_classes','lfine_scored_class_cell_fraction','marker_vocab_oracle_upper_bound','n_reference_lfine_disagreements']
MEASURES=['lfine_macroF1','coverage','legacy_called_coverage','abstain_rate','offvocab_rate','acc_on_called',
    'lfine_scored_class_cell_fraction','marker_vocab_oracle_upper_bound']
ABSTAIN={'Unknown','Undecided','Noise','nan','None',''}
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def js(p):return json.loads(Path(p).read_text())
def utc():return datetime.now(timezone.utc).isoformat()
def write(p,value):
    p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_name(p.name+f'.part.{os.getpid()}')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');tmp.replace(p)
def complete(p):
    p=Path(p);(p/'COMPLETE').write_text(sha(p/'manifest.json')+'\n')
def checked(p,name='manifest.json',flag='COMPLETE'):
    p=Path(p);return (p/name).is_file() and (p/flag).is_file() and (p/flag).read_text().strip()==sha(p/name)
def require_slurm():
    assert os.environ.get('SLURM_JOB_ID'),'All Lfine evaluation and numerical checks require SLURM'
    assert os.environ.get('PYTHONOPTIMIZE','0')=='0'
def contract():
    m=js(OUT/'contract.json');assert m['status']=='frozen_source_preparation_only'
    for p,h in m['source_hashes'].items():assert sha(p)==h,p
    assert sha(OUT/'shards.json')==m['shards_sha256']
    return m
def archive_gate(m,full=False):
    assert sha(ARCHIVE/'contract.json')==m['archive_contract_sha256']
    if not full:return None
    assert checked(ARCHIVE/'summary'),'Complete independent archival acceptance is mandatory before Lfine execution'
    proof=js(ARCHIVE/'summary/manifest.json')
    assert proof['status']=='passed' and proof['n_representations']==1694 and proof['n_candidates']==22022
    assert proof['validation_contract_sha256']==m['archive_contract_sha256']
    for name,digest in proof['outputs'].items():assert sha(ARCHIVE/'summary'/name)==digest
    links=js(ARCHIVE/'summary/unit_proof_hashes.json');assert len(links)==1694
    return links
def unit_acceptance(sample,budget,space,links):
    path=ARCHIVE/'units'/sample/budget/space/'manifest.json'
    assert checked(path.parent)
    if links is not None:assert sha(path)==links[str(path)]
    proof=js(path);assert proof['status']=='passed'
    assert proof['validation_contract_sha256']==sha(ARCHIVE/'contract.json') and proof['n_conditions']==13 and proof['n_metric_rows']==39
    for name,digest in proof['outputs'].items():assert sha(path.parent/name)==digest,name
    unit=EMBED/sample/budget/space
    for name,digest in proof['input_manifests'].items():assert sha(unit/name)==digest,name
    return path
def checked_outputs(directory):
    assert checked(directory);m=js(Path(directory)/'manifest.json')
    assert m['contract_sha256']==sha(OUT/'contract.json')
    for name,digest in m['outputs'].items():assert sha(Path(directory)/name)==digest,name
    for path,digest in m['inputs'].items():assert sha(path)==digest,path
    return m
