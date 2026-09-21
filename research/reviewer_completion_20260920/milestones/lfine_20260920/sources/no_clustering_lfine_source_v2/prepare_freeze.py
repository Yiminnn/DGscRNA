"""SLURM metadata-only source/protocol/input freeze; never evaluate or fit."""
from pathlib import Path
from datetime import datetime,timezone
import ast,hashlib,json,os,shutil
assert os.environ.get('SLURM_JOB_ID')
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
CAMP=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
OLD=CAMP/'no_clustering';REF=ROOT/'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917'
COMPACT=ROOT/'results/hvg_ptc_20260916_v1/lfine_compact_20260920'
OUT=CAMP/'no_clustering_lfine_v1';CODE=Path(__file__).resolve().parent;FROZEN=OUT/'source_v1'
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def read(p):return json.loads(Path(p).read_text())
def write_once(p,m):
    p=Path(p);assert not p.exists(),p;p.write_text(json.dumps(m,indent=2)+'\n')
def checked(p,manifest='manifest.json',flag='COMPLETE'):
    return (p/flag).read_text().strip()==sha(p/manifest)
for p in CODE.glob('*.py'):ast.parse(p.read_text(),filename=str(p))
original=read(OLD/'protocol.json');validation=read(OLD/'validation.json')
assert validation['status']=='passed' and validation['n_candidates']==1210
assert original['source_bundle_sha256']==sha(OLD/'source_v2/SOURCE_MANIFEST.json')
assert sha(OLD/'patient_folds.csv')==original['patient_folds_sha256']
compact=read(COMPACT/'manifest.json');assert compact['n_valid']==4719 and compact['historical_v5_endpoint_checks_passed']==24
# Preserve exact provider bytes; runtime verifies all original dependencies too.
vendor=CODE/'vendor';vendor.mkdir(exist_ok=False)
provider_paths={'evaluate_lfine.py':ROOT/'handoff/lfine_compact_20260920/evaluate_lfine.py',
    'v5_final_annotations.py':ROOT/'handoff/g274_table4/v5_final_annotations.py',
    'grid_sets.py':ROOT/'handoff/g274/grid_sets.py'}
for name,p in provider_paths.items():shutil.copy2(p,vendor/name)
pins={str(p.relative_to(ROOT)):sha(p) for p in provider_paths.values()}
for relative,digest in compact['source_hashes'].items():assert sha(ROOT/relative)==digest;pins[relative]=digest
for p in [ROOT/'handoff/NOTEBOOK_PRESENTATION_POLICY.json',ROOT/'handoff/lfine_compact_20260920/LFINE_EVALUATION.md',
    COMPACT/'manifest.json',COMPACT/'metrics_hvg24.csv.gz',COMPACT/'lfine_target_prefixes.json',
    OLD/'protocol.json',OLD/'source_v2/SOURCE_MANIFEST.json',OLD/'patient_folds.csv',OLD/'validation.json',OLD/'summary/manifest.json']:
    pins[str(p.relative_to(ROOT))]=sha(p)
assert sha(COMPACT/'metrics_hvg24.csv.gz')==compact['outputs']['metrics_hvg24.csv.gz']
inputs={'scope':'Immutable saved calls/manifests/reference-label files; no expression matrices loaded','samples':{}}
cache={}
def bind(p,files):
    p=Path(p);key=str(p.relative_to(ROOT))
    if key not in cache:cache[key]=sha(p)
    files[key]=cache[key]
    return cache[key]
for sample in original['samples']:
    files={};units={}
    for p in [REF/'evaluation_inputs'/sample/'truth.csv.gz',REF/'inputs'/sample/'input_manifest.json']:bind(p,files)
    for budget in original['budgets']:
        fit=OLD/'GBM'/sample/budget;assert checked(fit,'fit_manifest.json','FIT_COMPLETE')
        for name in ('fit_manifest.json','FIT_COMPLETE','config.json'):bind(fit/name,files)
        cfg=read(fit/'config.json');assert cfg['source_bundle_sha256']==original['source_bundle_sha256']
        prep=REF/'GBM'/sample/budget
        for name in ('prepare_manifest.json','PREPARED','cells.csv'):bind(prep/name,files)
        item={'cellwise':{'directory':str(fit/'cellwise_seed')},'reference':{}}
        routes=[(fit/'cellwise_seed',None)]
        for route in ('PCA30_SNN','PCA30_HDBSCAN_R','UMAP2_SNN','UMAP2_HDBSCAN_R'):
            directory=prep/route;score=read(directory/'score_manifest.json')
            aid=[aid for aid,a in score['arms'].items() if a['library']=='CM2_glioma_other' and a['cutoff']=='mean'];assert len(aid)==1
            item['reference'][route]={'directory':str(directory),'arm_id':aid[0]};routes.append((directory,aid[0]))
        for directory,aid in routes:
            assert checked(directory,'score_manifest.json','SCORE_COMPLETE')
            score=read(directory/'score_manifest.json')
            for name in ('score_manifest.json','SCORE_COMPLETE','initial_calls.csv.gz'):bind(directory/name,files)
            assert bind(directory/'initial_calls.csv.gz',files)==score['initial_sha256']
            arms=[aid] if aid else [a['id'] for a in original['arms']]
            for arm in arms:
                terminal=directory/'terminal'/arm
                for name in ('terminal_manifest.json','TERMINAL_COMPLETE','predictions.csv.gz','terminal.npz'):
                    if (terminal/name).exists():bind(terminal/name,files)
                if (terminal/'terminal_manifest.json').exists():
                    tm=read(terminal/'terminal_manifest.json');assert checked(terminal,'terminal_manifest.json','TERMINAL_COMPLETE')
                    assert tm['score_manifest_sha256']==sha(directory/'score_manifest.json')
                    if tm['terminal_valid']:
                        assert bind(terminal/'predictions.csv.gz',files)==tm['predictions_sha256']
                        assert bind(terminal/'terminal.npz',files)==tm['terminal_sha256']
        units[budget]=item
    inputs['samples'][sample]={'files':files,'units':units}
OUT.mkdir(parents=True,exist_ok=True);write_once(OUT/'inputs.json',inputs)
assert not FROZEN.exists();FROZEN.mkdir()
for p in CODE.rglob('*'):
    if p.is_file() and '__pycache__' not in p.parts:
        target=FROZEN/p.relative_to(CODE);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
source={'files':{str(p.relative_to(FROZEN)):sha(p) for p in FROZEN.rglob('*') if p.is_file()}}
write_once(FROZEN/'SOURCE_MANIFEST.json',source)
protocol=dict(frozen_at=datetime.now(timezone.utc).isoformat(),status='prepared_awaiting_root_review',
    scope='A2 saved-prediction endpoint revision; GBM HVG2000/5000 only; not all-gene main baseline',
    trigger='Latest user Lfine-only notebook request; post-hoc endpoint revision explicitly separated from original frozen L1 analysis',
    endpoint=compact['endpoint'],compatibility_limit=compact['compatibility_limit'],
    samples=original['samples'],budgets=original['budgets'],arms=original['arms'],marker='CM2_glioma_other',primary_lambda=1,
    primary_n_samples=97,primary_n_patients=55,all_n_samples=121,all_n_patients=59,
    primary_stage='terminal090',sensitivity_stage='terminal070',initial_stage='Diagnostic only; never final DG-scRNA result',
    original_A2_source_bundle_sha256=original['source_bundle_sha256'],original_L1_outputs_preserved=True,
    patient_folds_sha256=original['patient_folds_sha256'],
    selection='Same existing patientfolds; training-patient mean terminal090 Lfine score; samples averaged within patient; exact ties nearest1 then smallerlambda; no held-out labels choose lambda',
    inference=dict(bootstrap_replicates=2000,bootstrap_seed=42,unit='patient',wilcoxon_zero_method='pratt',
        alternative='two-sided',method='auto',holm_contrasts_per_cohort=16,total_contrasts=32,
        interval_scope='Conditional on saved predictions and selected configurations; no refitting or selection uncertainty'),
    invalid_policy='Retain unavailable/invalid/no-class rows and no-op statuses; never replace terminal with seeds; block complete ranking if endpoint unavailable rather than silently exclude',
    no_new_fitting=True,no_expression_matrix_loading=True,no_new_marker_mapping=True,
    source_manifest_sha256=sha(FROZEN/'SOURCE_MANIFEST.json'),input_snapshot_sha256=sha(OUT/'inputs.json'),
    pinned_sources=pins,metadata_freeze_job=os.environ['SLURM_JOB_ID'],metadata_freeze_step=os.environ.get('SLURM_STEP_ID'))
write_once(OUT/'protocol.json',protocol)
write_once(OUT/'status.json',dict(work_package='A',stage='A2_LFINE',status='preparing',
    completed=0,remaining=121,unit='saved-prediction sample re-evaluations',whole_work_package_A_complete=False,
    updated_at=datetime.now(timezone.utc).isoformat(),jobs=[],
    summary='Frozen source/input/endpoint revision prepared; awaiting root code review before any numerical evaluation',
    evidence=[str(OUT/'protocol.json'),str(FROZEN/'SOURCE_MANIFEST.json')]))
print(json.dumps(dict(status='prepared_no_numerical_evaluation',source=str(FROZEN),source_manifest_sha256=sha(FROZEN/'SOURCE_MANIFEST.json'),
    protocol_sha256=sha(OUT/'protocol.json'),input_snapshot_sha256=sha(OUT/'inputs.json'),n_input_files=len(cache))),flush=True)
