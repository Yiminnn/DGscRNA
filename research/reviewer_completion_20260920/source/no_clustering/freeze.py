"""Freeze allowed source, input and selection metadata before the first fit."""
import csv
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import shutil

assert os.environ.get('SLURM_JOB_ID'),'Freeze input provenance through SLURM'
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
CODE=ROOT/'handoff/reviewer_completion_20260920/no_clustering'
OUT=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/no_clustering'
REFERENCE=ROOT/'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917'
OUT.mkdir(parents=True,exist_ok=True)
def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()
def write(path,data):path.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
revision=os.environ.get('A2_FREEZE_REVISION','v1')
if (OUT/'protocol.json').exists():
    assert revision!='v1' and not any((OUT/'GBM').glob('*/*/FIT_COMPLETE')), 'No protocol replacement after candidate fits'
    prior=OUT/('protocol_before_'+revision+'.json')
    assert not prior.exists();shutil.copy2(OUT/'protocol.json',prior)
samples=(ROOT/'handoff/g274/cohort_samples.txt').read_text().split()
assert len(samples)==len(set(samples))==121
foldsource=REFERENCE/'protocol/patient_folds.csv'
folds=list(csv.DictReader(foldsource.open()))
assert {row['sample'] for row in folds}==set(samples)
assert all(len({r['fold'] for r in folds if r['patient']==row['patient']})==1 for row in folds)
sourcefiles={str(p.relative_to(CODE)):sha(p) for p in sorted(CODE.rglob('*')) if p.is_file() and p.suffix in ['.py','.R','.md','.sbatch']}
write(CODE/'SOURCE_MANIFEST.json',dict(frozen_at=datetime.now(timezone.utc).isoformat(),files=sourcefiles))
snapshot=OUT/('source' if revision=='v1' else 'source_'+revision)
shutil.copytree(CODE,snapshot,ignore=shutil.ignore_patterns('__pycache__'))
inputs={}
for sample in samples:
    im=json.loads((REFERENCE/'inputs'/sample/'input_manifest.json').read_text())
    for budget in ['hvg2000','hvg5000']:
        prep=REFERENCE/'GBM'/sample/budget;pm=json.loads((prep/'prepare_manifest.json').read_text())
        assert (prep/'PREPARED').read_text().strip()==sha(prep/'prepare_manifest.json')
        inputs[sample+'/'+budget]=dict(prepare_manifest_sha256=sha(prep/'prepare_manifest.json'),
            expression_sha256=pm['expression_sha256'],cells_sha256=sha(prep/'cells.csv'),
            DL_binary_sha256=pm['DL_binary_sha256'],n_cells=pm['n_cells'],DL_features=pm['features']['DL'],
            truth_sha256=im['evaluation_files']['truth.csv.gz'],
            display_sha256=sha(REFERENCE/'GBM'/sample/'hvg2000/UMAP2.csv'))
shutil.copy2(foldsource,OUT/'patient_folds.csv')
arms=[dict(id='lambda_'+str(value).replace('.','p'),lambda_value=value) for value in [0.,.5,1.,1.5,2.]]
for arm in arms:arm['lambda']=arm.pop('lambda_value')
protocol=dict(frozen_at=datetime.now(timezone.utc).isoformat(),scope='GBM-only A2 explicit seed-mechanism replacement',
    samples=samples,budgets=['hvg2000','hvg5000'],arms=arms,primary_lambda=1,
    primary_n_samples=sum(row['primary']=='True' for row in folds),
    primary_n_patients=len({row['patient'] for row in folds if row['primary']=='True'}),
    expected_candidate_conditions=1210,expected_sample_budget_units=242,
    scoring='positive normalized RNA-allgene marker mean / full panel length, singleton factor0.8; unique positive max; threshold lambda*sample mean max',
    reference_labels_used_for_fit=False,marker='CM2_glioma_other',
    marker_sha256=sha(REFERENCE/'markers/libraries.json'),mapping_sha256=sha(REFERENCE/'markers/panel_L1_mapping.csv'),
    patient_folds_sha256=sha(OUT/'patient_folds.csv'),source_bundle_sha256=sha(snapshot/'SOURCE_MANIFEST.json'),
    selection='training patient mean terminal090 macroF1; exact ties nearest1 then smallerlambda; no test labels',
    all_cells_retained=True,inputs=inputs)
write(OUT/'protocol.json',protocol)
write(OUT/'regression_tasks.json',[dict(sample='TKU4163',budgets=['hvg2000'])])
write(OUT/'pilot_tasks.json',[dict(sample='TKU4163',budgets=['hvg2000','hvg5000'])])
write(OUT/'full_tasks.json',[dict(sample=sample,budgets=['hvg2000','hvg5000']) for sample in samples])
print(json.dumps(dict(status='frozen',source=str(snapshot),samples=len(samples),units=len(inputs),candidates=1210)))
