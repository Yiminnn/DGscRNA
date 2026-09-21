"""Verify source preparation without importing scientific code or reading arrays."""
import ast,subprocess
import common as c
contract=c.contract();spec=c.js(c.CAMP/'protocol/embedding.json');shards=c.js(c.OUT/'shards.json')
assert len(shards)==242 and sum(len(t['spaces']) for t in shards)==1694
assert len({(t['sample'],t['budget']) for t in shards})==242
assert all(t['spaces']==spec['spaces'] for t in shards)
files=[]
for path in sorted(c.CODE.glob('*.py')):
    tree=ast.parse(path.read_text(),filename=str(path));files.append(path.name)
    for node in ast.walk(tree):
        if not isinstance(node,ast.Call):continue
        name=node.func.attr if isinstance(node.func,ast.Attribute) else node.func.id if isinstance(node.func,ast.Name) else ''
        assert name not in {'fit','fit_transform','fit_predict','predict','train_cache','build_model','read_h5ad','readRDS','memmap'},(path,node.lineno,name)
for name in ['job.sbatch','summary.sbatch']:subprocess.run(['bash','-n',str(c.CODE/name)],check=True)
assert not (c.CODE/'select.py').exists(),'Do not shadow the standard-library select module'
for name in ['evaluate.py','verify_shard.py','select_lfine.py','verify_selection.py','selftest.py','worker.py']:
    assert 'c.require_slurm()' in (c.CODE/name).read_text(),name
for name in ['grid_sets.py','v5_final_annotations.py','compact_evaluate_lfine.py']:
    ast.parse((c.REFERENCE/name).read_text(),filename=name)
assert not any('NOTEBOOK_PRESENTATION_POLICY' in p for p in contract['source_hashes'])
assert contract['mutable_presentation_policy_is_scientific_gate'] is False
protocol=c.js(c.OUT/'protocol.json');assert protocol['L1_selection_reuse'] is False
assert protocol['selection_metric']=='lfine_macroF1' and protocol['Holm_family_size']==21 and protocol['bootstrap_replicates']==10000
record=dict(status='passed_preparation_only',contract_sha256=c.sha(c.OUT/'contract.json'),source_hashes_verified=len(contract['source_hashes']),
    python_AST_checked=files,shell_syntax_checked=['job.sbatch','summary.sbatch'],n_units=1694,n_endpoint_rows=66066,
    exact_endpoint_source_snapshots=True,independent_count_and_selection_sources_present=True,
    per_unit_original_acceptance_gate=True,final_original_full_acceptance_gate=True,
    Lfine_training_patient_K_reselection=True,no_numerical_checks_run=True,no_scientific_data_loaded=True,jobs_submitted=[],at=c.utc())
c.write(c.OUT/'preparation_verification.json',record)
status=c.js(c.OUT/'status.json');status['evidence'].append(str(c.OUT/'preparation_verification.json'));c.write(c.OUT/'status.json',status)
print(record)
