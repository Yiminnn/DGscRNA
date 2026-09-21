"""Source/manifest-only checks; does not import scientific validator modules."""
import ast,json,subprocess
from collections import Counter
import common as c
contract=c.contract();tasks=c.js(c.OUT/'tasks.json');c.validate_tasks(tasks)
assert len(tasks)==1694
groups=Counter((t['sample'],t['budget']) for t in tasks)
assert len(groups)==242 and set(groups.values())=={7}
assert contract['n_candidates']==1694*13 and contract['n_metric_rows']==1694*13*3
spec=c.js(c.CAMP/'protocol/embedding.json')
assert spec['spaces']==c.SPACES and spec['K']==c.KS and set(spec['budgets'])=={'hvg2000','hvg5000'}
old=c.CAMP/'no_clustering/a1_first_review/verifier_v2/validate_unit.py'
assert (c.CODE/'pilot_baseline.py').read_text()==old.read_text().replace(
    "OUT=CAMPAIGN/'no_clustering/a1_first_review'/sample/budget/space",
    "OUT=CAMPAIGN/'embedding_full_validation_v1/pilot_baseline'/sample/budget/space")
parsed=[];forbidden=[]
for path in sorted(c.CODE.glob('*.py')):
    tree=ast.parse(path.read_text(),filename=str(path));parsed.append(path.name)
    for node in ast.walk(tree):
        if isinstance(node,ast.Call):
            name=node.func.attr if isinstance(node.func,ast.Attribute) else node.func.id if isinstance(node.func,ast.Name) else ''
            if name in {'fit','fit_transform','fit_predict','predict','train_cache','build_model','memmap','read_h5ad','readRDS'}:
                forbidden.append((path.name,node.lineno,name))
assert not forbidden,forbidden
for name in ['job.sbatch','summary.sbatch']:
    subprocess.run(['bash','-n',str(c.CODE/name)],check=True)
    text=(c.CODE/name).read_text()
    assert '#SBATCH --cpus-per-task=2' in text and '#SBATCH --mem=8G' in text
for name in ['worker.py','verify_summary.py']:
    assert 'c.require_slurm()' in (c.CODE/name).read_text()
assert len(c.pilot_tasks())==4 and set(t['sample'] for t in c.pilot_tasks())==set(c.PILOT_SAMPLES)
assert sum(t['space']=='ICA2' for t in c.pilot_tasks())==1
lfine=c.js(c.OUT/'lfine_followup_contract.json')
assert lfine['status']=='design_only_waiting_for_complete_L1_acceptance'
assert lfine['scientific_execution_or_Lfine_implementation_started'] is False
assert str(c.ROOT/'handoff/NOTEBOOK_PRESENTATION_POLICY.json') not in contract['source_hashes']
assert str(c.ROOT/'handoff/NOTEBOOK_PRESENTATION_POLICY.json') not in lfine['source_hashes']
for path,digest in lfine['source_hashes'].items():assert c.sha(path)==digest,path
proof=dict(status='passed_preparation_only',updated_at=c.utc(),n_source_hashes=len(contract['source_hashes']),
    parsed_python_files=parsed,bash_syntax_checked=['job.sbatch','summary.sbatch'],
    n_representations=1694,n_sample_budget_shards=242,n_candidates=22022,n_metric_rows=66066,
    baseline_verifier_only_output_path_changed=True,no_model_fitting_or_expression_loading_calls=True,
    scientific_source_hashes_valid=True,Lfine_design_only=True,mutable_display_policy_not_scientific_gate=True,
    validation_contract_sha256=c.sha(c.OUT/'contract.json'),tasks_sha256=c.sha(c.OUT/'tasks.json'),
    no_scientific_data_loaded=True,no_numerical_checks_run=True,jobs_submitted=[],
    remaining=['Real SLURM four-unit validator pilot','Root review and measured resource acceptance','Full sharded unit acceptance','Full independent summary acceptance','Later separate Lfine implementation and verification'])
c.write(c.OUT/'preparation_verification.json',proof)
status=c.js(c.OUT/'status.json');status['evidence'].append(str(c.OUT/'preparation_verification.json'));c.write(c.OUT/'status.json',status)
print(json.dumps(proof,indent=2))
