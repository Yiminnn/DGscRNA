"""Freeze metadata and source identity only. Never load arrays or submit jobs."""
from pathlib import Path
import ast,json
import common as c
assert not (c.OUT/'contract.json').exists(),'Do not rewrite an active validation contract'
c.OUT.mkdir(exist_ok=True)
spec=c.js(c.CAMP/'protocol/embedding.json')
tasks=[dict(sample=s,budget=b,space=p) for s in spec['samples'] for b in spec['budgets'] for p in spec['spaces']]
assert len(tasks)==1694 and len({(x['sample'],x['budget'],x['space']) for x in tasks})==1694
old=c.CAMP/'no_clustering/a1_first_review/verifier_v2/validate_unit.py'
text=old.read_text();replaced=text.replace("OUT=CAMPAIGN/'no_clustering/a1_first_review'/sample/budget/space","OUT=CAMPAIGN/'embedding_full_validation_v1/pilot_baseline'/sample/budget/space")
assert replaced==(c.CODE/'pilot_baseline.py').read_text(),'Baseline numerical checks changed beyond isolated output path'
source_paths=[c.CAMP/'protocol/embedding.json',c.CAMP/'protocol/embedding_selection.json',c.CAMP/'protocol/embedding_patient_folds.csv',c.POLICY,
    old,old.parent/'SOURCE_MANIFEST.json',c.REFERENCE/'markers/panel_L1_mapping.csv',c.REFERENCE/'markers/libraries.json',c.REFERENCE/'protocol/cohort.csv']
policy=c.js(c.POLICY)
for name in ['trigger_evidence','native_deflation_diagnostic','producer_compatibility']:
    item=policy[name];path=c.ROOT/item['path'];assert c.sha(path)==item['sha256'];source_paths.append(path)
old_failures={str(c.unitpath(t)/'CONVERGENCE_FAILURE.json'):c.sha(c.unitpath(t)/'CONVERGENCE_FAILURE.json')
    for t in tasks if t['space']=='ICA2' and (c.unitpath(t)/'CONVERGENCE_FAILURE.json').is_file()}
source_paths.extend(Path(p) for p in old_failures)
allowed={}
for version in ['embedding_v5','embedding_v6','embedding_v8']:
    directory=c.CAMP/'source_snapshots'/version;m=c.js(directory/'SOURCE_MANIFEST.json');source_paths.append(directory/'SOURCE_MANIFEST.json')
    for name,digest in m.items():assert c.sha(directory/name)==digest,name;source_paths.append(directory/name)
    allowed[c.sha(directory/'run.py')]=version
presentation=c.CAMP/'source_snapshots/embedding_presentation_v1'
source_paths.append(presentation/'SOURCE_MANIFEST.json')
for name,digest in c.js(presentation/'SOURCE_MANIFEST.json').items():assert c.sha(presentation/name)==digest;source_paths.append(presentation/name)
for p in c.CODE.iterdir():
    if p.is_file() and p.suffix in {'.py','.md','.sbatch'}:
        source_paths.append(p)
        if p.suffix=='.py':ast.parse(p.read_text(),filename=str(p))
c.write(c.OUT/'tasks.json',tasks)
refine_ast=ast.parse((c.CAMP/'source_snapshots/embedding_v8/legacy/legacy_refine.py').read_text())
params_node=next(n.value for n in refine_ast.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='PARAMS' for t in n.targets))
assert isinstance(params_node,ast.Call) and isinstance(params_node.func,ast.Name) and params_node.func.id=='dict'
dl_params={kw.arg:ast.literal_eval(kw.value) for kw in params_node.keywords}
dl_params['input']='Native R normalized RNA, selected genes in recorded order; not scaled expression'
contract=dict(status='prepared_for_gated_validation',endpoint='archived original L1 protocol; not canonical notebook display',created_at=c.utc(),
    source_hashes={str(p):c.sha(p) for p in source_paths},tasks_sha256=c.sha(c.OUT/'tasks.json'),allowed_run_producers=allowed,
    v8_run_sha256=c.sha(c.CAMP/'source_snapshots/embedding_v8/run.py'),
    adaptive_helper_sha256=c.sha(c.CAMP/'source_snapshots/embedding_v8/adaptive_ica.py'),
    refinement_source_sha256=c.sha(c.CAMP/'source_snapshots/embedding_v8/legacy/legacy_refine.py'),
    terminal_driver_sha256=c.sha(c.CAMP/'source_snapshots/embedding_v8/legacy/terminal.py'),DL_params=dl_params,
    baseline_reuse='Three-pilot verifier retained byte-for-byte except independent output directory',
    original_failure_records=old_failures,canonical_trigger_failure=policy['trigger_evidence'],
    n_representations=1694,n_candidates=22022,n_metric_rows=66066,default_shards=242,max_representations_per_shard=7,
    validation_resources=dict(cpus=2,memory_GB=8,time='04:00:00',suggested_total_concurrency=8),
    summary_resources=dict(cpus=2,memory_GB=8,time='02:00:00'),
    no_model_fitting=True,new_scientific_jobs_submitted=0,scientific_validator_status='Requires real SLURM pilot before broad validation release')
c.write(c.OUT/'contract.json',contract)
lfine_paths=[c.ROOT/'handoff/lfine_compact_20260920/evaluate_lfine.py',c.ROOT/'handoff/lfine_compact_20260920/LFINE_EVALUATION.md',
    c.ROOT/'handoff/g274/grid_sets.py',c.ROOT/'handoff/g274_table4/v5_final_annotations.py',
    c.ROOT/'results/hvg_ptc_20260916_v1/lfine_compact_20260920/lfine_target_prefixes.json',
    c.ROOT/'results/hvg_ptc_20260916_v1/lfine_compact_20260920/mapping_audit.csv.gz',c.REFERENCE/'markers/panel_L1_mapping.csv',
    c.ROOT/'results/hvg_ptc_20260916_v1/lfine_compact_20260920/manifest.json',
    c.ROOT/'results/hvg_ptc_20260916_v1/lfine_compact_20260920/historical_v5_endpoint_parity.csv',
    c.ROOT/'results/hvg_ptc_20260916_v1/lfine_compact_20260920/sample_validation.csv']
lfine=dict(status='design_only_waiting_for_complete_L1_acceptance',source_hashes={str(p):c.sha(p) for p in lfine_paths},
    delivery_policy_observation=dict(path=str(c.ROOT/'handoff/NOTEBOOK_PRESENTATION_POLICY.json'),
        observed_sha256=c.sha(c.ROOT/'handoff/NOTEBOOK_PRESENTATION_POLICY.json'),
        role='Mutable display policy, informative only; not a scientific hash prerequisite'),
    prerequisite_manifest=str(c.OUT/'summary/manifest.json'),prerequisite_flag=str(c.OUT/'summary/COMPLETE'),
    fits='Reuse identical native terminal initial/final090/final070 outputs; zero new model fits',
    endpoint='v5 compatible-target-set Lfine macro-F1; not strict one-to-one fine-label accuracy',
    scope=dict(representations=1694,candidates=22022,budgets=spec['budgets'],cohorts=['primary97','all121']),
    semantic_parent='Frozen CM2_glioma_other native panel semantic mapping; exact v5 mapping for this library',
    target_sets='Use all observed Lfine classes for fixed prefix compatibility, including classes below support20',
    class_average='Sample-present classes with support>=20 excluding Other/nan; all cells remain in TP/FP/FN',
    confusion_rule='correct iff reference fine label in predicted parent target set; TP=class&correct, FN=class&incorrect, FP=incorrect&claimed_class&not_reference_class',
    abstention='Unknown/Undecided/unmappable receive no compatible-class credit; all cells retained',
    missing='Invalid/unavailable terminals remain NA; classes absent from averaging are not discarded from cell counts; freeze truth-derived metric eligibility and report sample/patient denominators',
    K_selection='Recompute training-patient Lfine macro-F1 choices separately for each cohort/budget/space/method/fold, exact tie smallerK; same saved folds, full candidates, heldout labels excluded from K choice',
    secondary='Use the newly selected Lfine K for initial and terminal070 matched comparisons; never reuse L1-selected K as Lfine optimum',
    paired_inference='Recompute same-budget original R anchor with the same Lfine rule; patient means, bootstrap/Wilcoxon and per-cohort/budget Holm21; explicit metric-specific denominators',
    output_namespace='New endpoint/addendum and separate evaluation/selection/verification paths; do not alter frozen L1 evidence or notebook sources',
    notebook='Latest canonical display Lfine only; actual all-gene main comparison is separate native-R evidence; this A1 two-budget grid does not become an all-gene fit',
    next_implementation=['Pin immutable Lfine addendum and snapshots','Implement scientific evaluator and separate independent per-class validator','Re-select K from Lfine training patients','Independently verify full Lfine patient selection/inference'],
    scientific_execution_or_Lfine_implementation_started=False)
c.write(c.OUT/'lfine_followup_contract.json',lfine)
c.write(c.OUT/'status.json',dict(work_package='A',stage='A1_INDEPENDENT',status='prepared_not_executed',updated_at=c.utc(),jobs=[],
    completed=['1694-unit manifest','Pilot audit reuse with isolated outputs','Extended metrics/status/figures/solver checks','Independent patient selection and paired-statistic verifier','Separate gated Lfine follow-up design'],
    remaining=['Representative SLURM validator pilot','Full sharded read-only acceptance','Independent summary acceptance','Separate Lfine evaluator/addendum after L1 acceptance'],
    evidence=[str(c.OUT/'contract.json'),str(c.OUT/'tasks.json'),str(c.OUT/'lfine_followup_contract.json')]))
print(json.dumps(dict(contract=str(c.OUT/'contract.json'),n_tasks=len(tasks),baseline_numeric_source_unchanged=True,new_jobs=0),indent=2))
