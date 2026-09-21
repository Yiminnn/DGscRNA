"""Freeze Lfine source and label semantics only; no scientific reads or submission."""
import ast,shutil
import common as c
assert not (c.OUT/'contract.json').exists(),'Do not rewrite an active Lfine source freeze'
c.OUT.mkdir(parents=True,exist_ok=True);c.REFERENCE.mkdir(exist_ok=True)
spec=c.js(c.CAMP/'protocol/embedding.json');design=c.js(c.ARCHIVE/'lfine_followup_contract.json')
for p,h in design['source_hashes'].items():assert c.sha(p)==h,p
copies={
    'grid_sets.py':c.ROOT/'handoff/g274/grid_sets.py',
    'v5_final_annotations.py':c.ROOT/'handoff/g274_table4/v5_final_annotations.py',
    'compact_evaluate_lfine.py':c.ROOT/'handoff/lfine_compact_20260920/evaluate_lfine.py',
    'panel_L1_mapping.csv':c.NATIVE/'markers/panel_L1_mapping.csv',
    'lfine_target_prefixes.json':c.ROOT/'results/hvg_ptc_20260916_v1/lfine_compact_20260920/lfine_target_prefixes.json'}
for name,source in copies.items():
    target=c.REFERENCE/name
    if target.exists():assert c.sha(target)==c.sha(source)
    else:shutil.copyfile(source,target)
    assert c.sha(source)==c.sha(target)
shards=[dict(sample=s,budget=b,spaces=spec['spaces']) for s in spec['samples'] for b in spec['budgets']]
assert len(shards)==242
c.write(c.OUT/'shards.json',shards)
protocol=dict(endpoint='Frozen v5 compatible-target-set Lfine macro-F1',endpoint_source='Exact AST-extracted grid/v5 helpers used by compact evaluator',
    prediction_semantics='Native CM2_glioma_other panels retain their frozen ontology parent; no reference-derived one-to-one relabelling',
    class_average='Support >=20, excluding Other/nan; all observed classes define compatible targets',
    all_cells='All cells remain in TP/FP/FN; Unknown/Undecided and off-vocabulary predictions receive no compatible credit',
    stages=c.STAGES,selection_stage='terminal090',selection_metric='lfine_macroF1',
    candidates=dict(samples=121,budgets=spec['budgets'],spaces=spec['spaces'],K=spec['K'],methods=spec['clusterers'],partitions=22022),
    patient_folds_sha256=c.sha(c.CAMP/'protocol/embedding_patient_folds.csv'),cohorts=dict(primary97=dict(samples=97,patients=55),all121=dict(samples=121,patients=59)),
    selection='Within-patient mean across samples; equal training-patient mean by K; exact tie smaller K. Select afresh using Lfine only, separately per cohort/budget/space/method/fold',
    L1_selection_reuse=False,secondary='Use Lfine-selected K for marker-only and terminal070 comparisons',
    missing='Invalid/unavailable native terminals remain NA; final selector rejects any incomplete candidate. Truth-only zero-scored-class eligibility is explicit and must match across all candidates and anchor',
    denominators='Report all cells and observed classes, scored-class fraction, sample and patient totals, Lfine eligibility and per-metric nonmissing counts',
    anchor='Same-budget native R PCA30 to UMAP2/HDBSCAN, CM2_glioma_other/mean, final090, exact same Lfine rule',
    bootstrap_seed=20260920,bootstrap_replicates=10000,paired_test='Wilcoxon',Holm_family_size=21,
    inference='Retrospective patient-label-held-out K choice; paired inference conditional on selected predictions',
    original_acceptance_gate='Each Lfine candidate can be evaluated after matching independent original unit acceptance; final selection requires the complete original archival acceptance',
    final_acceptance='Every Lfine unit has an independent per-class/count proof, then an independent K/denominator/statistics proof',
    no_model_fits=True,notebook_changes=False,all_gene_claim=False)
c.write(c.OUT/'protocol.json',protocol)
sources=[c.OUT/'protocol.json',c.CAMP/'protocol/embedding.json',c.CAMP/'protocol/embedding_selection.json',c.CAMP/'protocol/embedding_patient_folds.csv',
    c.ARCHIVE/'contract.json',c.ARCHIVE/'lfine_followup_contract.json',c.NATIVE/'protocol/cohort.csv']
sources.extend(copies.values());sources.extend(c.REFERENCE/name for name in copies)
compact=c.ROOT/'results/hvg_ptc_20260916_v1/lfine_compact_20260920'
sources.extend(compact/name for name in ['manifest.json','mapping_audit.csv.gz','historical_v5_endpoint_parity.csv','sample_validation.csv'])
sources.extend(compact/'sample_metrics'/f'{s}.csv' for s in spec['samples'])
for path in c.CODE.iterdir():
    if path.is_file() and path.suffix in ['.py','.md','.sbatch']:
        sources.append(path)
        if path.suffix=='.py':ast.parse(path.read_text(),filename=str(path))
pilot=[dict(sample=s,budget='hvg2000',space='PCA2') for s in ['TKU4163','NL022','SN040']]+[dict(sample='TKU4163',budget='hvg2000',space='ICA2')]
contract=dict(status='frozen_source_preparation_only',endpoint=protocol['endpoint'],created_at=c.utc(),source_hashes={str(p):c.sha(p) for p in sources},
    archive_contract_sha256=c.sha(c.ARCHIVE/'contract.json'),shards_sha256=c.sha(c.OUT/'shards.json'),pilot_units=pilot,
    scope=dict(samples=121,shards=242,representations=1694,partitions=22022,endpoint_rows=66066),
    changed_gate_from_design='Root authorized 2026-09-21: per-unit Lfine replay after matching original independent unit acceptance; full original acceptance still required for final selection',
    mutable_presentation_policy_is_scientific_gate=False,source_only=True,jobs_submitted=[],
    resources=dict(unit_shards='2 CPUs, 8 GB, 2 hours, suggested array concurrency 8 after measured pilot',summary='2 CPUs, 8 GB, 1 hour'),
    execution_status='No Lfine numerical test or evaluation has run in this namespace')
c.write(c.OUT/'contract.json',contract)
c.write(c.OUT/'status.json',dict(work_package='A',stage='A1_LFINE',status='source_prepared_not_executed',updated_at=c.utc(),jobs=[],
    completed=['Separate endpoint and immutable source contract','Per-unit original-acceptance gate','Exact endpoint evaluator plus independent per-class checker','New Lfine train-patient K selector plus independent verifier'],
    remaining=['Root source review','SLURM numerical and saved-unit pilots','Full candidate replay and per-class verification','Full selection and independent inference verification'],
    evidence=[str(c.OUT/'contract.json'),str(c.OUT/'protocol.json')]))
print(c.OUT/'contract.json')
