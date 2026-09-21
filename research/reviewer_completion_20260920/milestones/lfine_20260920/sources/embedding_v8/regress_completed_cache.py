"""SLURM-only isolated completed-cache validation/tamper tests; no fits/truth."""
from pathlib import Path
import copy,json,os,shutil,sys
sys.path.insert(0,str(Path(__file__).resolve().parent))
import run as core
core.require_slurm();core.cache_compatibility.own_source(core)
reference=core.OUTPUT/'TKU4163/hvg2000/PCA2'
assert core.checked(reference,'fit_manifest.json','FIT_COMPLETE')
base=core.CAMPAIGN/'embedding_v8_regression/completed_cache'
assert not base.exists(),'Preserve prior regression outputs'
base.mkdir(parents=True)
dest=base/'fixture';dest.mkdir()
for name in ('config.json','representation.json','REPRESENTATION_COMPLETE','embedding.csv','fit_manifest.json','FIT_COMPLETE'):
    shutil.copy2(reference/name,dest/name)
original_cfg=json.loads((dest/'config.json').read_text());cfg=copy.deepcopy(original_cfg)
cfg.update(dest=str(dest),embedding=str(dest/'embedding.csv'),hdbscan_dest=str(dest/'HDBSCAN_R'))
for condition in cfg['conditions']:
    old=Path(condition['dest']);new=dest/old.name;new.mkdir()
    for name in ('clusters.csv','partition_manifest.json','PARTITION_COMPLETE','score_input_fingerprint.json','score_manifest.json','SCORE_COMPLETE','initial_calls.csv.gz'):
        shutil.copy2(old/name,new/name)
    terminal=new/'terminal/L00_mean';terminal.mkdir(parents=True)
    for name in ('terminal_manifest.json','TERMINAL_COMPLETE','terminal.npz','predictions.csv.gz'):
        shutil.copy2(old/'terminal/L00_mean'/name,terminal/name)
    condition['dest']=str(new)
previous=json.loads((dest/'fit_manifest.json').read_text());previous['config']=cfg
core.write_json(dest/'config.json',cfg);core.write_json(dest/'fit_manifest.json',previous)
core.complete(dest,'fit_manifest.json','FIT_COMPLETE')
# This is a path-relocated TEST FIXTURE, not a new scientific producer.
core.write_json(base/'FIXTURE_PROVENANCE.json',dict(purpose='Isolated path-relocated validation fixture; not a scientific result',
    reference=str(reference),reference_fit_manifest_sha256=core.sha(reference/'fit_manifest.json'),
    config_path_fields_relocated=True,representation_predictions_and_score_bytes_unchanged=True))
original_config_fn=core.config
core.config=lambda sample,budget,space:cfg
# If a shortcut ever falls through into fitting, fail before any model is run.
core.make_embedding=lambda *args:(_ for _ in ()).throw(AssertionError('Unexpected fitting'))
core.finish=lambda *args:(_ for _ in ()).throw(AssertionError('Unexpected DL'))
checks=[]
def invoke():core.run('TKU4163','hvg2000','PCA2')
def test(name,mutate,expected):
    restore=mutate()
    try:
        invoke()
    except AssertionError as exc:
        assert expected in str(exc),(name,str(exc))
        checks.append(dict(name=name,status='rejected',error=str(exc)))
    else:raise AssertionError('Tamper accepted: '+name)
    finally:restore()
invoke();checks.append(dict(name='valid_complete_cache',status='accepted'))
def mutate_bytes(path):
    data=path.read_bytes();path.write_bytes(data+b'\n')
    return lambda:path.write_bytes(data)
test('embedding_csv_changed',lambda:mutate_bytes(dest/'embedding.csv'),'Changed embedding CSV')
# Re-sign only representation hash to prove cell-order guard acts independently.
def reorder_embedding():
    e=dest/'embedding.csv';m=dest/'representation.json';flag=dest/'REPRESENTATION_COMPLETE'
    originals={p:p.read_bytes() for p in (e,m,flag)}
    lines=e.read_text().splitlines();lines[1],lines[2]=lines[2],lines[1];e.write_text('\n'.join(lines)+'\n')
    meta=json.loads(m.read_text());meta['embedding_sha256']=core.sha(e);core.write_json(m,meta);core.complete(dest,'representation.json','REPRESENTATION_COMPLETE')
    return lambda:[p.write_bytes(data) for p,data in originals.items()]
test('embedding_cell_order_resigned',reorder_embedding,'Embedding cell order differs')
route=dest/'KMeans_K05'
test('terminal_prediction_csv_changed',lambda:mutate_bytes(route/'terminal/L00_mean/predictions.csv.gz'),'Changed terminal prediction CSV')
for field in ('clusters_sha256','prepare_manifest_sha256','cells_sha256','expression_sha256','marker_sha256','scorer_sha256','protocol_sha256'):
    def mutate_fingerprint(field=field):
        p=route/'score_input_fingerprint.json';data=p.read_bytes();m=json.loads(data);m[field]='0'*64;core.write_json(p,m)
        return lambda:p.write_bytes(data)
    test('score_fingerprint_'+field,mutate_fingerprint,'Changed complete score input fingerprint')
# Valid reuse again proves all fixture mutations restored.
invoke();checks.append(dict(name='restored_complete_cache',status='accepted'))
core.config=original_config_fn
report=dict(status='passed',sample='TKU4163',budget='hvg2000',space='PCA2',n_conditions=13,
    checks=checks,source_manifest_sha256=core.sha(core.CODE/'SOURCE_MANIFEST.json'),
    source_sha256=core.sha(core.CODE/'run.py'),cache_helper_sha256=core.sha(core.CODE/'cache_compatibility.py'),
    reference=str(reference),reference_fit_manifest_sha256=core.sha(reference/'fit_manifest.json'),
    scientific_model_fits=0,author_labels_read=False,production_caches_modified=False,
    job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),completed_at=core.utc())
core.write_json(base/'validation.json',report);print(json.dumps(report,indent=2),flush=True)
