"""Pinned saved-prediction inputs for the requested A2 Lfine endpoint revision."""
from pathlib import Path
from datetime import datetime,timezone
import hashlib,importlib.util,json,os,sys
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
CAMP=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
CODE=Path(__file__).resolve().parent
OUT=CAMP/'no_clustering_lfine_v1'
OLD=CAMP/'no_clustering'
REFERENCE=ROOT/'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917'
COMPACT=ROOT/'results/hvg_ptc_20260916_v1/lfine_compact_20260920'
FIXED='CM2_glioma_other'
ROUTES=['PCA30_SNN','PCA30_HDBSCAN_R','UMAP2_SNN','UMAP2_HDBSCAN_R']
STAGES={'initial':'initial','terminal090':'final090','terminal070':'final070'}
FIELDS=['lfine_macroF1','coverage','legacy_called_coverage','abstain_rate','offvocab_rate','acc_on_called',
        'n_distinct_calls','n_classes_hit','lfine_n_classes','lfine_scored_class_cell_fraction',
        'marker_vocab_oracle_upper_bound','n_reference_lfine_disagreements']
def sha(path):
    with Path(path).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(path):return json.loads(Path(path).read_text())
def utc():return datetime.now(timezone.utc).isoformat()
def write(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_name(path.name+f'.part.{os.getpid()}');tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n');tmp.replace(path)
def checked(path,name='manifest.json',flag='COMPLETE'):
    path=Path(path);return (path/flag).is_file() and (path/name).is_file() and (path/flag).read_text().strip()==sha(path/name)
def complete(path):
    path=Path(path);(path/'COMPLETE').write_text(sha(path/'manifest.json')+'\n')
def verify(require_approval=True):
    assert os.environ.get('SLURM_JOB_ID') and not sys.flags.optimize
    sm=read(CODE/'SOURCE_MANIFEST.json')
    for name,digest in sm['files'].items():assert sha(CODE/name)==digest,('source',name)
    p=read(OUT/'protocol_v2.json');assert p['source_manifest_sha256']==sha(CODE/'SOURCE_MANIFEST.json')
    for name,digest in p['pinned_sources'].items():assert sha(ROOT/name)==digest,('pinned source',name)
    assert sha(OUT/'inputs.json')==p['input_snapshot_sha256']
    if require_approval:
        a=read(OUT/'APPROVED.json')
        assert a['approved'] is True and a['protocol_sha256']==sha(OUT/'protocol_v2.json')
        assert a['source_manifest_sha256']==sha(CODE/'SOURCE_MANIFEST.json')
    return p

def load_provider():
    # Import the exact compact provider copy; initialization performs no output
    # writes or expression loading. Its original v5/grid helpers are pinned too.
    spec=importlib.util.spec_from_file_location('a2_pinned_lfine_provider',CODE/'vendor/evaluate_lfine.py')
    provider=importlib.util.module_from_spec(spec);spec.loader.exec_module(provider)
    provider.initialize()
    assert sha(provider.v5.__file__)==sha(CODE/'vendor/v5_final_annotations.py')
    return provider

def verify_files(files):
    for name,digest in files.items():assert sha(ROOT/name)==digest,('saved input',name)

def bool_column(frame):
    import pandas as pd
    assert pd.api.types.is_bool_dtype(frame.primary.dtype) and frame.primary.notna().all(), 'Primary flags must be parsed booleans'
