"""Validate and apply a process-local scCATCH singleton-dimension repair.

The original audited caller and installed package stay unchanged. Every generated
caller records its exact source and validation proof before publishing its flag.
"""
import json
import os
from pathlib import Path
import subprocess
import sys

from common import OUT, RSCRIPT, require_slurm, checked, complete, sha, write_json, utc

SOURCE = Path(__file__).resolve().parent
BASE_SHA = 'a86b044d83f12e126cda246fb4d52cf8b422f9e90b2eee541af5222af6be29ed'
AUDIT = OUT/'verification/scCATCH_dimension_guard_audit'
NATIVE = OUT/'comparators/scCATCH/TKU4163/hvg2000/PCA30_SNN'


def generated_caller(attempt, proof_path, output=None):
    require_slurm()
    assert sha(SOURCE/'sccatch_R.R') == BASE_SHA
    proof = json.loads(proof_path.read_text())
    guard_path = SOURCE/'sccatch_dimension_guard.R'
    assert proof['status'] == 'passed' and proof['guard_source_sha256'] == sha(guard_path)
    text = (SOURCE/'sccatch_R.R').read_text()
    anchor = 'suppressPackageStartupMessages(library(scCATCH))\n'
    assert text.count(anchor) == 1
    quote = lambda p: json.dumps(str(p))
    setup = (
        f'guard_proof_path <- {quote(proof_path)}\n'
        f'stopifnot(digest(file=guard_proof_path,algo="sha256")=="{sha(proof_path)}")\n'
        'guard_proof <- fromJSON(guard_proof_path,simplifyVector=FALSE)\n'
        f'stopifnot(digest(file={quote(guard_path)},algo="sha256")=="{sha(guard_path)}")\n'
        f'source({quote(guard_path)})\n'
        'dimension_guard <- install_sccatch_dimension_guard(lapply(guard_proof$guard$functions,function(f)f$original_sha256))\n'
        f'dimension_guard$source_sha256 <- "{sha(guard_path)}"\n'
        'dimension_guard$validation_manifest <- guard_proof_path\n'
        f'dimension_guard$validation_manifest_sha256 <- "{sha(proof_path)}"\n'
    )
    text = text.replace(anchor, anchor + setup)
    if output is not None:
        line = "dest<-file.path(base,'comparators/scCATCH',sample,budget,route)"
        assert text.count(line) == 1
        text = text.replace(line, f'dest<-{quote(output)}')
    anchor = "write_json(m,file.path(dest,paste0(mode,'_manifest.json')),pretty=TRUE,auto_unbox=TRUE)"
    assert text.count(anchor) == 1
    text = text.replace(anchor, f'm$compatibility_guard <- dimension_guard\nm$unmodified_caller_sha256 <- "{BASE_SHA}"\n' + anchor)
    attempt.mkdir(parents=True, exist_ok=True)
    caller = attempt/'sccatch_guarded_caller.R'
    if caller.exists():
        assert caller.read_text() == text
    else:
        caller.write_text(text)
    write_json(attempt/'source_manifest.json', dict(
        caller=str(caller), caller_sha256=sha(caller), base_caller_sha256=BASE_SHA,
        guard=str(guard_path), guard_sha256=sha(guard_path),
        validation_manifest=str(proof_path), validation_sha256=sha(proof_path),
        driver_sha256=sha(__file__), job=os.environ['SLURM_JOB_ID']))
    return caller


def audit(unit_dir):
    require_slurm()
    unit_dir = Path(unit_dir)
    assert checked(unit_dir, 'unit_manifest.json', 'UNIT_COMPLETE')
    assert checked(NATIVE, 'audit_manifest.json', 'AUDIT_COMPLETE')
    original = json.loads((NATIVE/'audit_manifest.json').read_text())
    assert original['source_sha256'] == BASE_SHA
    proof_path = unit_dir/'unit_manifest.json'
    unit = json.loads(proof_path.read_text())
    assert unit['source_sha256'] == sha(SOURCE/'audit_sccatch_dimensions.R')
    if checked(AUDIT):
        m = json.loads((AUDIT/'manifest.json').read_text())
        assert m['guard_source_sha256'] == sha(SOURCE/'sccatch_dimension_guard.R')
        return
    attempt = AUDIT/'attempts'/os.environ['SLURM_JOB_ID']
    output = attempt/'real_control'
    caller = generated_caller(attempt, proof_path, output)
    subprocess.run([RSCRIPT, str(caller), 'TKU4163', 'cohort', 'hvg2000', 'PCA30_SNN'], check=True)
    assert checked(output, 'cohort_manifest.json', 'COHORT_COMPLETE')
    subprocess.run([RSCRIPT, str(SOURCE/'audit_sccatch_guard_real.R'), str(NATIVE), str(output)], check=True)
    parity_path = output/'real_parity_manifest.json'
    parity = json.loads(parity_path.read_text())
    assert parity['status'] == 'passed'
    write_json(AUDIT/'manifest.json', dict(status='passed', completed_at=utc(),
        job=os.environ['SLURM_JOB_ID'], guard=unit['guard'],
        guard_source_sha256=sha(SOURCE/'sccatch_dimension_guard.R'), base_caller_sha256=BASE_SHA,
        unit_manifest=str(proof_path), unit_manifest_sha256=sha(proof_path),
        real_parity_manifest=str(parity_path), real_parity_manifest_sha256=sha(parity_path),
        real_control_manifest_sha256=sha(output/'cohort_manifest.json'),
        original_audit_manifest_sha256=sha(NATIVE/'audit_manifest.json'),
        real_parity=parity, source_sha256=sha(__file__),
        scope='Two column-subset drop=FALSE arguments; original tests and scoring unchanged'))
    complete(AUDIT)
    print('SCCATCH_DIMENSION_GUARD_AUDIT_PASSED', flush=True)


def recover(sample):
    require_slurm()
    assert checked(AUDIT)
    proof = json.loads((AUDIT/'manifest.json').read_text())
    for name in ['unit_manifest', 'real_parity_manifest']:
        assert sha(proof[name]) == proof[name + '_sha256']
    assert proof['original_audit_manifest_sha256'] == sha(NATIVE/'audit_manifest.json')
    dest = OUT/'comparators/scCATCH'/sample/'hvg2000/UMAP2_HDBSCAN_R'
    failure = dest/'attempts/singleton_failure_receipt.json'
    assert failure.exists(), 'A reviewed, preserved original failure is required'
    fm = json.loads(failure.read_text())
    assert fm['sample'] == sample and fm['native_caller_sha256'] == BASE_SHA
    assert fm['exact_error'] == 'dim(X) must have a positive length'
    assert sha(fm['log']) == fm['log_sha256']
    assert fm['cluster_sha256'] == sha(OUT/'GBM'/sample/'hvg2000/UMAP2_HDBSCAN_R/clusters.csv')
    if not checked(dest, 'cohort_manifest.json', 'COHORT_COMPLETE'):
        caller = generated_caller(dest/'attempts'/os.environ['SLURM_JOB_ID'], AUDIT/'manifest.json')
        subprocess.run([RSCRIPT, str(caller), sample, 'cohort', 'hvg2000', 'UMAP2_HDBSCAN_R'], check=True)
    assert checked(dest, 'cohort_manifest.json', 'COHORT_COMPLETE')
    result = json.loads((dest/'cohort_manifest.json').read_text())
    assert result['unmodified_caller_sha256'] == BASE_SHA
    assert result['compatibility_guard']['validation_manifest_sha256'] == sha(AUDIT/'manifest.json')
    import evaluate_comparator
    evaluate_comparator.run('scCATCH', sample)


if __name__ == '__main__':
    mode, arg = sys.argv[1:]
    {'audit': audit, 'recover': recover}[mode](arg)
