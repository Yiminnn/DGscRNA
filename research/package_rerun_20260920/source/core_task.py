"""One frozen installed-package GBM task; no scheduling or legacy-cache import."""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import uuid
import zipfile

HERE = Path(__file__).resolve().parent
CONFIG = dict(preset='gbm-reference',dataset='GSE274546',route='all',library='all',
              cutoff='all',seed=42,deg_workers=4,require_slurm=True)
PILOTS = {
    'anchor': ('TKU4163','hvg2000',1),
    'full_roster': ('TKU4163','hvg2000',192),
    'medium_all_genes': ('NL022','all',1),
    'large_hvg2000': ('SN040','hvg2000',1),
}
REQUIRED_LAUNCHERS = {'manage_core.py','core_task.py','core_array.sbatch',
                     'verify_packaged_parity.py','verify_packaged_r_artifacts.R',
                     'evaluate_core_lfine.py','evaluate_core_lfine_sources.json'}


def need(value, message):
    if not value: raise RuntimeError(message)


def load(path):
    return json.loads(Path(path).read_text())


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def check_file(record):
    path=Path(record['path']).resolve()
    need(path.is_file() and sha(path)==record['sha256'],f'Frozen artifact mismatch: {path}')
    return path


def write_new(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as handle:
        json.dump(value,handle,indent=2,allow_nan=False);handle.write('\n')


def installed_files(package):
    return {str(p.relative_to(package)):sha(p) for p in sorted(package.rglob('*'))
            if p.is_file() and '__pycache__' not in p.parts and p.suffix not in {'.pyc','.pyo'}}


def validate_gate(path):
    gate=load(path)
    need(gate.get('status')=='release_and_pilots_verified','Release and four-pilot gate is not approved')
    need(gate['config']==CONFIG,'Frozen core scientific configuration differs from required protocol')
    need(gate['publication']['status']=='published','GitHub publication has not completed')
    need(gate['publication']['branch']=='align-r-reference','Publication branch differs from authorized branch')
    commit=gate['publication']['commit']
    need(len(commit)==40 and all(c in '0123456789abcdef' for c in commit),'Missing publication commit')
    check_file(gate['publication']['receipt'])
    wheel=check_file(gate['wheel'])
    tasks_path=check_file(gate['task_manifest']);check_file(gate['markers'])
    need(REQUIRED_LAUNCHERS.issubset(gate['launcher_files']), 'Gate is missing a required frozen launcher/evaluator')
    for entry in gate['launcher_files'].values():check_file(entry)
    need(Path(gate['launcher_files']['core_task.py']['path']).resolve()==Path(__file__).resolve(),
         'Task launcher is not the frozen reviewed file')
    need(Path(sys.prefix).resolve()==Path(gate['runtime']['python']).parent.parent.resolve(),
         'This process is not the frozen installed runtime prefix')
    need(Path(sys.executable).absolute()==Path(gate['runtime']['python']).absolute(),
         'Run this worker using the frozen installed Python executable')
    spec=importlib.util.find_spec('dgscrna')
    need(spec is not None and spec.origin is not None,'Installed package cannot be located')
    package=Path(spec.origin).resolve().parent
    need(package==Path(gate['installed_package']['root']).resolve(),'Wrong package import location')
    need(package.is_relative_to(Path(sys.prefix).resolve()),'Package was imported outside installed prefix')
    need(installed_files(package)==gate['installed_package']['files_sha256'],
         'Installed package contents differ from frozen wheel installation')
    with zipfile.ZipFile(wheel) as archive:
        payload={name.removeprefix('dgscrna/'):hashlib.sha256(archive.read(name)).hexdigest()
                 for name in archive.namelist() if name.startswith('dgscrna/') and not name.endswith('/')}
    need(payload==gate['installed_package']['files_sha256'],'Installed package differs from the published wheel payload')
    proofs={item['role']:item for item in gate['pilots']}
    need(len(proofs)==len(gate['pilots']) and set(proofs)==set(PILOTS),'Missing/duplicate required pilot proofs')
    for role,(sample,budget,n) in PILOTS.items():
        proof_path=check_file(proofs[role]);proof=load(proof_path)
        need(proof.get('status')=='passed_exact' and proof['sample']==sample
             and proof['budget']==budget and proof['terminal_conditions']==n,f'Invalid pilot proof: {role}')
        plan_path=proof_path.parent/'run_config.json'
        need(sha(plan_path)==proof['run_config_sha256'],f'Pilot configuration changed: {role}')
        plan=load(plan_path)
        need(plan['runtime']==gate['runtime']['fingerprint'],f'Pilot runtime differs: {role}')
        for name,digest in plan['sources'].items():
            need(gate['installed_package']['files_sha256'].get('reference/'+name)==digest,
                 f'Pilot package source differs: {role}/{name}')
        if role=='anchor':need(proof['trained_conditions']==1,'Anchor did not exercise real DL training')
        if role=='full_roster':
            need(len(proof['route_checks'])==4 and proof['valid_no_training_conditions']>0,
                 'Full pilot did not cover all routes and valid no-training states')
    tasks=load(tasks_path)
    need(len(tasks)==726 and len({(t['sample'],t['budget']) for t in tasks})==726,'Incomplete/duplicate core task manifest')
    return gate,tasks,package


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gate',type=Path,required=True)
    parser.add_argument('--index',type=int,default=None)
    parser.add_argument('--validate-only',action='store_true',help='Verify gate, package, runtime and input hashes without fitting or creating a unit')
    args=parser.parse_args()
    need(os.environ.get('SLURM_JOB_ID'),'All scientific task checks and fitting require SLURM')
    need(sys.flags.no_user_site and not sys.flags.optimize,'Use installed Python -s without -O')
    need(not os.environ.get('PYTHONPATH'),'Remove PYTHONPATH before running the installed package')
    os.chdir('/tmp')
    index=args.index if args.index is not None else int(os.environ['SLURM_ARRAY_TASK_ID'])
    gate,tasks,package=validate_gate(args.gate.resolve())
    need(0<=index<len(tasks),'Array index outside frozen task manifest')
    task=tasks[index];need(task['index']==index,'Task manifest index mismatch')
    need(set(task['routes'])=={'PCA30_SNN','PCA30_HDBSCAN_R','UMAP2_SNN','UMAP2_HDBSCAN_R'},'Missing core routes')
    need(task['cutoffs']==['none','mean','0.5'],'Unexpected cutoff roster')
    need(Path(task['markers']).resolve()==Path(gate['markers']['path']).resolve(),'Task markers differ from gate')
    for name,digest in task['input_files_sha256'].items():
        need(sha(Path(task['counts_dir'])/name)==digest,f'Count input changed: {task["sample"]}/{name}')
    from dgscrna.reference.runner import doctor, run_reference
    runtime=doctor(gate['runtime']['rscript'],gate['runtime'].get('reference_r_lib'))
    need(runtime==gate['runtime']['fingerprint'],'Installed R/Python runtime differs from accepted pilots')
    if args.validate_only:
        print(json.dumps(dict(status='gate_package_runtime_inputs_verified_no_fit',index=index,
                              sample=task['sample'],budget=task['budget'],gate_sha256=sha(args.gate))),flush=True)
        return
    root=Path(gate['output_root']).resolve();output=(root/'core'/task['sample']/task['budget']).resolve()
    need(output.is_relative_to(root/'core'),'Unsafe unit output path')
    attempt=uuid.uuid4().hex
    audit=root/'control/attempts'/f'{index:04d}_{os.environ["SLURM_JOB_ID"]}_{attempt}.json'
    lock=root/'control/locks'/f'{task["sample"]}__{task["budget"]}'
    lock.parent.mkdir(parents=True,exist_ok=True)
    try:lock.mkdir()
    except FileExistsError as error:
        raise RuntimeError(f'Existing worker lock: {lock}; inspect its job before retry') from error
    write_new(lock/'owner.json',dict(job=os.environ['SLURM_JOB_ID'],array_task=index,host=os.uname().nodename,pid=os.getpid()))
    base=dict(sample=task['sample'],budget=task['budget'],index=index,output=str(output),
              gate_sha256=sha(args.gate),wheel_sha256=gate['wheel']['sha256'],
              package_root=str(package),job=os.environ['SLURM_JOB_ID'],
              array_job=os.environ.get('SLURM_ARRAY_JOB_ID'),array_task=os.environ.get('SLURM_ARRAY_TASK_ID'),
              requested_cpus=os.environ.get('SLURM_CPUS_PER_TASK'),memory_per_node=os.environ.get('SLURM_MEM_PER_NODE'),
              started_at=datetime.now(timezone.utc).isoformat())
    try:
        if (output/'PACKAGE_VERIFIED_COMPLETE').exists():
            acceptance=load(output/'package_acceptance.json')
            need((output/'PACKAGE_VERIFIED_COMPLETE').read_text().strip()==sha(output/'package_acceptance.json'),
                 'Existing package acceptance marker is invalid')
            need(acceptance['gate_sha256']==sha(args.gate),'Existing run used a different gate')
            need(acceptance['run_manifest_sha256']==sha(output/'run_manifest.json')
                 and acceptance['parity_sha256']==sha(output/'packaged_parity.json'),'Verified outputs have changed')
            completed=load(output/'run_manifest.json')
            need((output/'COMPLETE').read_text().strip()==sha(output/'run_manifest.json'),'Package completion receipt changed')
            for name,digest in completed['exports'].items():
                need(sha(output/name)==digest,f'Previously verified export changed: {name}')
            evaluation=Path(acceptance['lfine_manifest'])
            need(sha(evaluation)==acceptance['lfine_manifest_sha256']
                 and (evaluation.parent/'COMPLETE').read_text().strip()==acceptance['lfine_manifest_sha256'],
                 'Previously verified Lfine receipt changed')
            need(sha(evaluation.parent/'metrics.csv.gz')==acceptance['lfine_metrics_sha256'],
                 'Previously verified Lfine metrics changed')
            print(json.dumps(dict(status='already_verified',**base)),flush=True)
            return
        need(not output.exists(),f'Unverified output already exists: {output}; preserve it and have the root review retry')
        manifest=run_reference(counts=task['counts_dir'],markers=gate['markers']['path'],out=str(output),
                   sample=task['sample'],features=task['budget'].removeprefix('hvg'),
                   rscript=gate['runtime']['rscript'],reference_r_lib=gate['runtime'].get('reference_r_lib'),**CONFIG)
        need(manifest['terminal_condition_count']==192,'Packaged core unit is missing terminal conditions')
        # Only exact inputs inside this new unit may share a fresh training cache.
        caches=set()
        for record in manifest['conditions']:
            terminal=output/Path(record['predictions']).parent
            tm=load(terminal/'terminal_manifest.json')
            cache=Path(tm['cache_directory']).resolve()
            need(cache.is_relative_to(output/'DL_cache'),'Terminal borrowed an old/external cache')
            training=load(terminal/'training_manifest.json')
            first=Path(training['provenance']['first_condition']).resolve()
            need(first.is_relative_to(output),'Training cache first condition belongs to another run')
            caches.add(str(cache))
        verifier=check_file(gate['launcher_files']['verify_packaged_parity.py'])
        helper=check_file(gate['launcher_files']['verify_packaged_r_artifacts.R'])
        cmd=[gate['runtime']['python'],'-s',str(verifier),'--output',str(output),
             '--reference-root',gate['reference_root'],'--r-helper',str(helper),'--allow-fresh-cache-reuse']
        with (output/'independent_verification.log').open('x') as log:
            result=subprocess.run(cmd,cwd='/tmp',stdout=log,stderr=subprocess.STDOUT)
        need(result.returncode==0,'Independent exact parity failed; inspect independent_verification.log')
        proof=load(output/'packaged_parity.json')
        need(proof['status']=='passed_exact' and proof['terminal_conditions']==192,'Incomplete parity proof')
        evaluator=check_file(gate['launcher_files']['evaluate_core_lfine.py'])
        semantics=check_file(gate['launcher_files']['evaluate_core_lfine_sources.json'])
        need(semantics.parent==evaluator.parent,'Evaluator semantic lock must accompany its frozen script')
        evaluation=root/'evaluation/units'/task['sample']/task['budget']
        need(not evaluation.exists(),'Existing unverified Lfine outputs require review before overwrite')
        command=[gate['runtime']['python'],'-s',str(evaluator),'unit','--task-index',str(index),
                 '--campaign-root',str(root)]
        with (output/'lfine_evaluation.log').open('x') as log:
            result=subprocess.run(command,cwd='/tmp',stdout=log,stderr=subprocess.STDOUT)
        need(result.returncode==0,'Frozen Lfine evaluation failed; invalid states preserved without acceptance')
        evaluation_hash=sha(evaluation/'manifest.json');em=load(evaluation/'manifest.json')
        need((evaluation/'COMPLETE').read_text().strip()==evaluation_hash,'Invalid Lfine completion checksum')
        need(em['status']=='completed' and em['n_conditions']==384 and em['n_valid']==384,
             'Lfine evaluation is missing valid terminal threshold rows')
        need(em['task_index']==index and em['sample']==task['sample'] and em['budget']==task['budget'],
             'Lfine receipt identifies a different unit')
        need(em['script_sha256']==sha(evaluator) and em['frozen_semantics']==load(semantics),
             'Lfine evaluation used different frozen semantics')
        need(sha(evaluation/'metrics.csv.gz')==em['outputs']['metrics.csv.gz'],'Lfine metric checksum mismatch')
        accepted=dict(status='package_run_independently_verified',**base,
                      trained_conditions=proof['trained_conditions'],valid_no_training_conditions=proof['valid_no_training_conditions'],
                      terminal_conditions=192,cache_paths_confined_to_unit=True,n_unique_local_caches=len(caches),
                      run_manifest_sha256=sha(output/'run_manifest.json'),parity_sha256=sha(output/'packaged_parity.json'),
                      lfine_threshold_rows=384,lfine_valid_threshold_rows=384,
                      lfine_manifest=str(evaluation/'manifest.json'),lfine_manifest_sha256=evaluation_hash,
                      lfine_metrics_sha256=em['outputs']['metrics.csv.gz'],
                      completed_at=datetime.now(timezone.utc).isoformat())
        write_new(output/'package_acceptance.json',accepted)
        with (output/'PACKAGE_VERIFIED_COMPLETE').open('x') as flag:flag.write(sha(output/'package_acceptance.json')+'\n')
        write_new(audit,accepted)
        print(json.dumps({k:accepted[k] for k in ['status','index','sample','budget','terminal_conditions','trained_conditions','valid_no_training_conditions']}),flush=True)
    except Exception as error:
        write_new(audit,dict(status='failed_preserved',**base,error=str(error),failed_at=datetime.now(timezone.utc).isoformat()))
        raise
    finally:
        (lock/'owner.json').unlink();lock.rmdir()


if __name__=='__main__':main()
