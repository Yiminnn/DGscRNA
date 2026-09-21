"""Metadata-only v8 freeze; v2 numerical policy and v7 regression stay immutable."""
from pathlib import Path
import ast,hashlib,json,os,shutil
assert os.environ.get('SLURM_JOB_ID')
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
CAMP=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
CODE=Path(__file__).resolve().parent
V7=CAMP/'source_snapshots/embedding_v7'
FROZEN=CAMP/'source_snapshots/embedding_v8'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def write_once(path,data):
    assert not path.exists(),f'Immutable output exists: {path}'
    path.write_text(json.dumps(data,indent=2)+'\n')
def functions(path):
    source=path.read_text()
    return {node.name:(node,ast.get_source_segment(source,node)) for node in ast.parse(source).body if isinstance(node,ast.FunctionDef)}
old=functions(V7/'run.py');new=functions(CODE/'run.py')
unchanged_functions=['call_r','get_geometry','model_metadata','fit_model','make_embedding','make_partitions','finish','parity']
for name in unchanged_functions:assert old[name][1]==new[name][1],name
# run() changes only the completed-fit branch; fresh fitting code remains exact.
assert old['run'][1][old['run'][1].index('    x, cells, gm = get_geometry(cfg)'):]==new['run'][1][new['run'][1].index('    x, cells, gm = get_geometry(cfg)'):]
assert sha(CODE/'adaptive_ica.py')==sha(V7/'adaptive_ica.py')
catalog=json.loads((V7/'CACHE_COMPATIBILITY.json').read_text())
for version,item in catalog['approved_previous_sources'].items():
    base=ROOT/item['directory'];manifest=json.loads((base/'SOURCE_MANIFEST.json').read_text())
    assert sha(base/'SOURCE_MANIFEST.json')==item['manifest_sha256']
    for name in catalog['unchanged_dependency_files']:assert sha(CODE/name)==manifest[name]==sha(base/name)
shutil.copy2(V7/'CACHE_COMPATIBILITY.json',CODE/'CACHE_COMPATIBILITY.json')
reg=CAMP/'embedding_v7_regression/TKU4163/hvg2000/validation.json'
r=json.loads(reg.read_text());assert r['status']=='passed' and r['embedding_csv_bytes_exact'] and r['all_13_partition_labels_and_cells_exact']
assert r['source_manifest_sha256']==sha(V7/'SOURCE_MANIFEST.json')
policy=CAMP/'protocol/embedding_convergence_repair_20260920_v2.json'
assert sha(policy)=='ad69004d945eecbe850cdf0c5c3c1f3fc8167d601e4d40092625c8ec3355e458'
proof=dict(status='passed',scope='V8 completed-cache hardening only; all fresh numerical paths byte-exact v7',
    unchanged_functions={name:hashlib.sha256(new[name][1].encode()).hexdigest() for name in unchanged_functions},
    unchanged_fresh_run_tail=True,adaptive_helper_sha256=sha(CODE/'adaptive_ica.py'),
    inherited_numerical_regression=dict(path=str(reg.relative_to(ROOT)),sha256=sha(reg),
        source_manifest_sha256=sha(V7/'SOURCE_MANIFEST.json'),embedding_csv_bytes_exact=True,all_13_partitions_exact=True),
    policy_sha256=sha(policy),compatibility_catalog_sha256=sha(CODE/'CACHE_COMPATIBILITY.json'),
    completed_cache_tamper_regression_required=True,activation='Not activated; root gate review required',
    job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'))
write_once(CODE/'SOURCE_COMPATIBILITY_PROOF.json',proof)
assert not FROZEN.exists();FROZEN.mkdir()
for path in CODE.rglob('*'):
    if path.is_file() and '__pycache__' not in path.parts and path.name!='SOURCE_MANIFEST.json':
        target=FROZEN/path.relative_to(CODE);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target)
write_once(FROZEN/'SOURCE_MANIFEST.json',{str(p.relative_to(FROZEN)):sha(p) for p in FROZEN.rglob('*') if p.is_file()})
print(json.dumps(dict(status='prepared_frozen_not_activated',path=str(FROZEN),manifest_sha256=sha(FROZEN/'SOURCE_MANIFEST.json'),policy_sha256=sha(policy))),flush=True)
