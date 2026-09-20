"""Protocol inventory from existing metadata only; no expression or predictions loaded."""
from pathlib import Path
import csv, json, gzip, os, hashlib
from datetime import datetime, timezone
ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
CODE=ROOT/'handoff/reviewer_completion_20260920/comparison'
OUT=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/comparison'
OLD=ROOT/'results/hvg_ptc_20260916_v1/r_reference_campaign_20260917'
CLAIM=ROOT/'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917'
MODEL=Path('/fs/scratch/PCON0080/yimin/tools/deepsort-1.0/deepsort-pretrained')
METHODS=['DG-scRNA','scType','scCATCH','SCINA','SingleR','CellTypist','CHETAH','scmap','scDeepSort','SignacX']
SOURCES={}
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def track(p):
    p=Path(p);SOURCES[str(p)]=sha(p);return p
def js(p):return json.loads(track(p).read_text())
def rows(p):
    p=track(p)
    with (gzip.open(p,'rt') if p.suffix=='.gz' else p.open()) as f:return list(csv.DictReader(f))
def outcsv(name,records):
    keys=list(dict.fromkeys(k for r in records for k in r))
    with (OUT/name).open('w') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(records)
def writejson(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')
def utc():return datetime.now(timezone.utc).isoformat()

def run():
    assert os.environ.get('SLURM_JOB_ID')
    OUT.mkdir(parents=True,exist_ok=True)
    inv=js(OLD/'dataset_inventory.json')
    assert len(inv)==11
    roster=rows(OLD/'summary/marker_context_roster.csv')
    completion=rows(OLD/'summary/completion_inventory.csv')
    compby={r['unit']:r for r in completion}
    truth=rows(OLD/'summary/curated_truth_label_mapping.csv')
    markerrows=rows(OLD/'markers/native_panel_metadata.csv')
    for path in [OLD/'summary/campaign_summary.json',CLAIM/'comparison_summary/manifest.json',
        CLAIM/'PTC_comparator_replay/manifest.json',CLAIM/'legacy_coverage_audit/manifest.json']:
        js(path)
    legacy=rows(ROOT/'results/revision/annotation_metrics.csv')+rows(ROOT/'results/bench_summary.csv')
    libraries={r['unit']:r for r in roster}
    unitinfo=[]
    for dataset in inv:
        units=js(OLD/'inputs'/dataset/'units.json')
        for u in units:
            names=set(u['batch_sizes'])
            known=names-{'unknown','Unknown','nan','NA','None',''}
            n=len(known)
            actual=compby.get(u['unit'])
            if actual is None:raise ValueError(u['unit'])
            unitinfo.append(dict(dataset=dataset,unit=u['unit'],scope=u['scope'],n_cells=u['n_cells'],
                batch_column=u['batch_column'],n_declared_batches=len(names),n_nonmissing_batch_ids=n,
                declared_batch_ids=';'.join(sorted(names)),input_semantics=u['input_semantics'],expression_layer=u['expression_layer'],
                donor_identity='metadata_declares_donor' if 'donor' in u['batch_column'].lower() or u['batch_column']=='Patient' else 'sample_batch_ID_requires_patient_identity_confirmation',
                split_rule='Freeze donor-grouped folds across all units; training labels only for selection' if n>=2 else
                    'Single nonmissing unit donor: fixed marker/config or independently sourced reference; no cell-level pseudo-donor split',
                n_curated_truth_terms=len({r['truth'] for r in truth if r['unit']==u['unit']}),
                n_common_lineages=len({r['common_lineage'] for r in truth if r['unit']==u['unit']}),
                marker_libraries=libraries[u['unit']]['libraries'],n_marker_libraries=libraries[u['unit']]['n_libraries'],
                correction=actual['correction'],features_scoring=actual['features_scoring'],features_DL=actual['features_DL'],
                existing_DG_terminal_conditions=actual['terminal_complete'],existing_DG_prediction_directory=actual['directory'],
                input_manifest=str(Path(u['path'])/'input_manifest.json'),
                marker_file=str(OLD/'markers'/f"{u['unit']}.json")))
    assert len(unitinfo)==69 and sum(u['dataset']=='HCL' for u in unitinfo)==59
    gbm=js(CLAIM/'protocol/input_audit.json')
    unitinfo.append(dict(dataset='GBM_GSE274546',unit='GBM_121_samples',scope='121 samples / primary97',n_cells=gbm['n_cells'],
        batch_column='patient',n_declared_batches=59,n_nonmissing_batch_ids=59,declared_batch_ids='protocol/cohort.csv',
        input_semantics='integer_counts',expression_layer='frozen eligible-gene counts',donor_identity='59 patient IDs already grouped in frozen 5 folds',
        split_rule='Reuse frozen patient folds; no sample of a patient crosses folds',n_curated_truth_terms=11,n_common_lineages=11,
        marker_libraries='16 frozen libraries; 13 eligible for primary selection',n_marker_libraries=16,
        correction='per-sample RNA, no CCA',features_scoring='all eligible RNA genes',features_DL='HVG2000 matched comparison',
        existing_DG_terminal_conditions='121 samples x 24 geometry/budget configurations x 48 marker/cutoff arms',
        existing_DG_prediction_directory=str(CLAIM/'GBM'),input_manifest=str(CLAIM/'protocol/input_audit.json'),
        marker_file=str(CLAIM/'markers/libraries.json')))
    ptcref=ROOT/'results/hvg_ptc_20260916_v1/ptc_experiments/evaluation_reference/manifest.json'
    js(ptcref)
    for group,scope in [('NMT','MT-1;MT-2;N-1;N-2'),('TTU','TU-1;TU-2;T-1;T-2')]:
        anchors=[r for r in completion if r['unit']==f'PTC_{group}_CCA2000']
        if not anchors:anchors=[r for r in completion if group in r['unit'] and 'CCA2000' in r['unit']]
        u=anchors[0] if anchors else {}
        unitinfo.append(dict(dataset='PTC',unit=group,scope=scope,n_cells=u.get('n_cells','group-count-in-existing-PTC-roster'),
            batch_column='patient plus sample',n_declared_batches=4,n_nonmissing_batch_ids=4,declared_batch_ids='Patient1;Patient2;Patient3;Patient4',
            input_semantics='raw GEX counts available; historical shared cell set 92404',expression_layer='counts / unintegrated normalized RNA per tool contract',
            donor_identity='N-1/T-1=P1;N-2/T-2=P2;MT-1/TU-1=P3;MT-2/TU-2=P4',
            split_rule='Same patient remains in one fold across both groups; integrated DG is transductive unless independently refit',
            n_curated_truth_terms='S2/S3 labels separate; assay endpoint binary TCR',n_common_lineages='strict T versus nonT with abstention separate',
            marker_libraries='Original17 including Thyroid; Pubmed_34663816; HPA; relevant extra-thyroid CellMarker; All',n_marker_libraries=17,
            correction='CCA/NONE/Harmony kept distinct',features_scoring='depends on original assay branch',features_DL='original2000 or named controls',
            existing_DG_terminal_conditions='30 old units plus 22 follow-up control units; terminal090',
            existing_DG_prediction_directory=str(CLAIM/'PTC_summary'),input_manifest=str(ptcref),marker_file=str(OLD/'markers/PTC_original17.json')))
    outcsv('analysis_units.csv',unitinfo)
    # Checkpoint metadata only: never open a weights tensor or expression graph.
    checkpoints=[]
    for p in sorted((MODEL/'human/models').glob('human-*.pt')):
        tissue=p.stem.removeprefix('human-');label=MODEL/'human/statistics'/f'{tissue}_cell_type.txt'
        gene=MODEL/'human/statistics'/f'{tissue}_genes.txt';graph=MODEL/'human/graphs'/f'human_{tissue}_data.npz'
        labs=track(label).read_text().splitlines() if label.exists() else []
        checkpoints.append(dict(tissue=tissue,checkpoint=str(p),checkpoint_bytes=p.stat().st_size,
            labels_file=str(label),n_labels=len(labs),labels=';'.join(labs),genes_file=str(gene),graph_file=str(graph),
            complete_bundle=label.exists() and gene.exists() and graph.exists(),checkpoint_load_currently_tested=tissue=='Brain',
            training_overlap_policy='HCL is source training atlas; Pancreas additionally trained with Baron GSE84133'))
    outcsv('scDeepSort_local_checkpoints.csv',checkpoints)
    impl={
      'DG-scRNA':('native_R','handoff/r_reference_campaign_20260917/reference_prepare.R;handoff/r_reference_campaign_20260917/reference_score.R;handoff/r_reference_campaign_20260917/terminal.py','OriginalR environment already executed','terminal090'),
      'scType':('matched_marker','handoff/paper_claim_validation_20260917/sctype_R.R','Official source parity tested; adapt GBM-hardcoded caller','official cluster score/cell thresholds'),
      'scCATCH':('matched_marker','handoff/paper_claim_validation_20260917/sccatch_R.R','scCATCH3.2.2 vendor_R; singleton guard independently tested','official marker labels; ties preserved'),
      'SCINA':('matched_marker','handoff/paper_claim_validation_20260917/scina_R.R','SCINA1.2.0 vendor_R; audited boundary guard','per-cell label; rm_overlap1 anchor /0 sensitivity'),
      'SingleR':('supervised_reference','handoff/paper_claim_validation_20260917/singler_R.R','annobench R; GBM frozen references available','pruned primary; unpruned sensitivity'),
      'CellTypist':('supervised_reference','handoff/methods/run_py_methods.py','Historical CellTypist1.7.1 local donor-trained run; current import smoke pending','no-majority-voting native predictions'),
      'CHETAH':('supervised_reference','handoff/methods/run_r_methods.R','annobench R package DESCRIPTION exists','native celltype_CHETAH incl abstention'),
      'scmap':('supervised_reference','handoff/methods/run_r_methods.R','annobench R package DESCRIPTION exists; audit counts slot/logcounts assignment before reuse','scmapCluster incl unassigned'),
      'scDeepSort':('pretrained_atlas','handoff/paper_claim_validation_20260917/deepsort_predict.py','Python3.8 CPU env installed; input correction required','native unsure_rate=2 / number_of_classes; not DG0.90'),
      'SignacX':('pretrained_atlas','handoff/ptc/run_signacx.R','deconv_r2 SignacX2.2.5; old saved outputs available','CellStates strict T vs CellTypes TNK separate')}
    inputreq={
      'DG-scRNA':'NativeR feature/assay contract; scoring and DL universes distinct',
      'scType':'RNA lognorm then gene scaling; same partition; identical candidate panels; gs2=None because shared roster lacks negatives',
      'scCATCH':'Unintegrated log-normalized RNA + matched cluster partition, original gene-level testing, valid marker metadata',
      'SCINA':'Unintegrated log-normalized RNA + same positive marker signatures; preserve all cells and numerical guard statuses',
      'SingleR':'Unintegrated log-normalized query and labelled training-reference expression; exact shared genes',
      'CellTypist':'Log1p library-size normalized expression and labelled training-reference cells; do not reuse unavailable immune model by default',
      'CHETAH':'Compatible SingleCellExperiment expression and independent labelled reference; native Unknown retained',
      'scmap':'Compatible SCE query/reference, native feature-selection inputs must be audited; no target labels in index construction',
      'scDeepSort':'Default Seurat LogNormalize BEFORE model gene intersection; real-valued RNA, no int32 cast; model-native graph normalization is additional',
      'SignacX':'Raw count-bearing compatible Seurat assay with normalized data and graph as documented; preserve native label hierarchy'}
    implementation=[]
    for method,(info,paths,environment,endpoint) in impl.items():
        for path in paths.split(';'):track(ROOT/path)
        implementation.append(dict(method=method,information_condition=info,implementation=paths,
            environment_evidence=environment,input_contract=inputreq[method],endpoint=endpoint,
            portability='Existing executable is GBM/brain/PTC-specific; new dataset adapter and validation required' if method!='DG-scRNA' else 'public11 already executed'))
    outcsv('implementation_inventory.csv',implementation)
    unitmatrix=[]
    for u in unitinfo:
      ds=u['dataset']
      for method in METHODS:
        info,paths,environment,endpoint=impl[method]
        r=dict(dataset=ds,unit=u['unit'],method=method,n_target_cells=u['n_cells'],
            information_condition=info,existing_status='missing_in_audited_summary_scope',existing_fit_or_cache='no confirmed run',
            existing_cells='unverified',missing_cells_for_new_contract='all target cells unless existing-run eligibility is established',
            eligible_for_primary_new_table=False,reference_availability='',vocabulary_support='',input_appropriateness=inputreq[method],
            split_rule=u['split_rule'],batch_column=u['batch_column'],n_nonmissing_batch_ids=u['n_nonmissing_batch_ids'],
            implementation=paths,needed_work='',evidence=u['input_manifest'],scope_policy='after_GBM_gate')
        if method=='DG-scRNA':
            r.update(existing_status='native_R_terminal_complete',existing_fit_or_cache='existing fits and valid no-op terminal statuses',
                existing_cells=u['n_cells'],missing_cells_for_new_contract=0,eligible_for_primary_new_table=True,
                reference_availability=u['marker_file'],vocabulary_support='Frozen native marker names + curated/common-lineage mapping; absent classes remain errors',
                needed_work='Reuse current final090 and frozen candidate roster; retain transductive label-heldout qualifier',evidence=u['existing_DG_prediction_directory'])
        elif info=='matched_marker':
            r.update(reference_availability=u['marker_file'],vocabulary_support='Panel universe exists; per-unit actual gene support/zero-signature conditions need exported audit',
                needed_work='Adapt official audited caller to same retained cells, frozen original-R partition and all prespecified libraries; run training-donor selection only')
        elif info=='supervised_reference':
            r.update(reference_availability='Curated donor-labelled cells exist; fold reference bundles not built for this unit',
                vocabulary_support='Class support must be tabulated from TRAIN donors only; retain unsupported test classes as errors',
                needed_work='Freeze patient/study folds, export unintegrated expression/reference bundle, audit class support and run native tool')
            if int(u['n_nonmissing_batch_ids'])<2:
                r.update(reference_availability='No within-unit donor-heldout reference: one declared donor; select independent compatible existing study before fit',
                    needed_work='Use externally sourced compatible reference or retain explicit unresolved-reference status; never split cells as pseudo-patients')
            if ds=='PTC':
                r.update(reference_availability='Existing HCL AdultThyroid candidate reference; S2/S3 donor labels are alternative concordance-only reference, not independent truth',
                    needed_work='Freeze external thyroid/immune reference and broad vocabulary; preserve four-patient boundary; TCR only a detection endpoint')
        elif method=='scDeepSort':
            models={'GBM_GSE274546':'Brain','brain_GBM':'Brain','PTC':'Thyroid','baron_human':'Pancreas','muraro':'Pancreas','segerstolpe':'Pancreas','xin':'Pancreas','kidney_ccRCC':'Kidney','colorectal':'Colorectum','blood_DLBCL':'Blood;Bone_marrow;Spleen','immune_ALL_human':'Blood;Bone_marrow'}
            r.update(reference_availability='Local candidate checkpoint(s): '+models.get(ds,'See checkpoint inventory; exact anatomy/age match unresolved'),
                vocabulary_support='Checkpoint-native labels; malignant and unsupported classes not removed or renamed using target truth',
                needed_work='Audit reference provenance, select anatomy before metrics, export correct lognorm float input; native threshold unchanged')
            if ds in ['HCL','baron_human']:
                r.update(existing_status='pretrained_training_source_overlap',reference_availability='Local weights use HCL; Pancreas additionally used Baron GSE84133',
                    needed_work='Do not count pretrained result as independent evaluation. Independent donor-excluded GNN fit would be a separate future protocol; no new fit authorized here')
            if ds=='xin':r['needed_work']+='; RPKM is not raw counts: validate input-specific handling, never round RPKM to integer'
        elif method=='SignacX':
            r.update(reference_availability='Installed SignacX2.2.5 pretrained hierarchy; exact reference support requires method-native vocabulary audit',
                vocabulary_support='Blood/immune hierarchy; strict T states distinct from coarse TNK; no assumption of tumor/organ coverage',
                needed_work='Audit native input and class coverage before declaring applicability; unsupported full-class comparison remains visible')
        if ds=='GBM_GSE274546' and method in ['scType','scCATCH','SCINA','SingleR']:
            r.update(existing_status='current_comparator_complete',existing_fit_or_cache='native/audited per-sample predictions',
                existing_cells=429305,missing_cells_for_new_contract=0,eligible_for_primary_new_table=True,
                evidence=str(CLAIM/'comparison_summary/manifest.json'),needed_work='Reuse frozen predictions and heldout selection; no full rerun')
            if method=='SingleR':r['reference_availability']=str(CLAIM/'reference_inputs')+' / 5 training-patient folds'
        if ds=='GBM_GSE274546' and method=='scDeepSort':
            r.update(existing_status='existing_input_contract_mismatch',existing_fit_or_cache='121 native predictions on raw int32 counts',existing_cells=429305,
                missing_cells_for_new_contract=429305,needed_work='Corrected GBM-only LogNormalize replay first TKU4163/NL022; preserve old predictions; then scale based on pilot resources',
                evidence=str(CLAIM/'comparators/scDeepSort')+'; '+str(CODE/'SCDEEPSORT_INPUT_AUDIT.md'))
        if ds=='GBM_GSE274546' and method=='CellTypist':
            r['existing_fit_or_cache']='Old TKU3186 pretrained-model download failed; no cohort predictions demonstrated'
            r['evidence']=str(ROOT/'results/competitors/TKU3186/status.json')
        if ds=='PTC' and method in ['scType','scCATCH','SCINA','SignacX']:
            r.update(existing_status='historical_cache_reevaluated_not_matched_refit',existing_fit_or_cache='Original S3 labels; SignacX also later per-sample cached fit',
                existing_cells='92404 shared cells across both groups',missing_cells_for_new_contract='new matched group/method contract not demonstrated',
                evidence=str(CLAIM/'PTC_comparator_replay/manifest.json'),
                needed_work='Keep original paper table; add separate matched-input analysis, candidate/reference audit and patient split; never overwrite historical labels')
            if method=='SignacX':r['needed_work']+='; HISTORICAL ZERO T CALLS IMMUTABLE, later real predictions not forced to zero'
        if ds in inv and method=='scType' and any(x.get('dataset')==ds and x.get('method')=='scType' for x in legacy):
            r.update(existing_status='legacy_sctype_style_not_current_official_contract',existing_fit_or_cache='Historical Python scType-style scores on old workflow',
                existing_cells='see saved old per-config n_eval; not assumed full current cell set',evidence=str(ROOT/'results/revision/annotation_metrics.csv')+'; handoff/bench_lib.py',
                needed_work='Retain legacy result; run official-audited R scType against current original-R input/partition and marker roster')
        if ds=='brain_GBM':
            r['scope_policy']='user_excluded_new_Darmanis_tuning_or_refitting; preserve/audit existing only'
            if method in ['SingleR','CellTypist','CHETAH','scmap']:
                r.update(existing_status='historical_saved_reference_predictions',existing_fit_or_cache='existing leave-one-donor-out fits + full-denominator replay',
                    existing_cells=3589,evidence=str(CLAIM/'legacy_coverage_audit/manifest.json'),needed_work='Use existing audit with old label/denominator limits; no new fit')
            elif method=='scDeepSort':
                r.update(existing_status='historical_saved_GNN_3567_cells',existing_fit_or_cache='old published Brain model predictions; current normalization eligibility unverified',
                    existing_cells=3567,evidence=str(ROOT/'results/b7/scdeepsort_meta.json'),needed_work='Preserve with 22-cell coverage gap and input caveat; no new Darmanis fit')
            else:r['needed_work']='Preserve existing outputs or absent condition transparently; user excludes new Darmanis tuning/fits'
        if ds=='xin':r['input_appropriateness']+='; SOURCE IS RPKM, not verified raw counts: explicit published-expression condition required'
        unitmatrix.append(r)
    outcsv('frozen_unit_method_matrix.csv',unitmatrix)
    datasetmatrix=[]
    for ds in list(inv)+['GBM_GSE274546','PTC']:
      for method in METHODS:
        subset=[r for r in unitmatrix if r['dataset']==ds and r['method']==method]
        base={k:' | '.join(sorted({str(r[k]) for r in subset})) for k in subset[0] if k not in ['unit','n_target_cells','n_nonmissing_batch_ids']}
        base.update(dataset=ds,method=method,n_analysis_units=len(subset),units=';'.join(r['unit'] for r in subset),
            n_target_cells=599926 if ds=='HCL' else 92404 if ds=='PTC' else subset[0]['n_target_cells'])
        datasetmatrix.append(base)
    outcsv('frozen_dataset_method_matrix.csv',datasetmatrix)
    assert len(datasetmatrix)==130 and len(unitmatrix)==720
    policy='''# Frozen comparison inventory — audit complete, work package C remains open

The inventory covers 13 dataset families (11 public datasets, GBM GSE274546 and
PTC), 72 primary analysis units (including 59 HCL tissues and two PTC groups),
and ten named methods: 130 dataset-method / 720 unit-method entries. It reads
existing manifests, summary tables, label dictionaries and source code only.
No expression matrix, cell predictions or model weights were loaded. No new
PTC/public fit, re-evaluation or plot was run. Absence means no qualifying result
in the explicitly inspected summary scope, not a recursive proof of file absence.

## Reuse versus work still needed

GBM current DG/scType/scCATCH/SCINA/SingleR outputs have the matching frozen
comparison. scDeepSort's saved raw-int32-count input violates the published
LogNormalize contract; preserve it as historical and replay correctly. Remaining
reference tools have portable implementations but no qualifying GBM cohort
output was located. Their labelled reference condition must remain separate.

Public11 DG terminal outputs exist for all 69 public analysis units. Historical
Python scType-style results do not establish the official method on the current
R geometry/input/marker roster. Other public comparator gaps need real work,
not a renamed result. Public marker and cell-class dictionaries already exist;
the unit table records the actual retained counts, input semantics and donors.

PTC original comparator caches and later SignacX per-sample predictions have
been re-evaluated on the 92,404-cell intersection. They are not new matched-input
group refits. Keep original S2/S3 and paper SignacX zero-T outputs immutable.
Use the four linked patient IDs across both groups. TCR supports a detection
endpoint, not independent multiclass truth. Reference-trained PTC results need
explicit external-reference or source-label concordance interpretation.

## Frozen execution rules

All new PTC/public computation waits for the current GBM gate. No new Darmanis
fits/tuning are allowed. Missing models/references remain unresolved until
actual provenance/support is audited; tissue name alone does not establish
applicability or inapplicability. HCL's 59 units are not 59 independent cohorts;
the same donor must stay in one fold across all units. Single-donor tissue units
need a fixed configuration or independently sourced reference, never a random
cell split labelled patient validation. Existing integrated DG is transductive.

Xin is published RPKM, not raw counts. Immune_ALL_human retains mixed UMI and
full-length quantification provenance. Algorithms must respect these inputs;
do not round RPKM, rebrand it as counts, or use corrected/integrated expression
where an annotator expects original RNA. Every method reports all-cell coverage,
Unknown and unsupported truth classes, with training-only library selection.

scDeepSort's actual available bundles are in scDeepSort_local_checkpoints.csv.
HCL and Baron are training-source-overlap conditions for published checkpoints;
they cannot be called independent external GNN tests. Checkpoint availability
does not prove a tissue model's environment load or class support. Only Brain
has a recorded loaded checkpoint here. See SCDEEPSORT_INPUT_AUDIT.md for official
evidence, precise code mismatch and the authorized GBM-only replay protocol.
'''
    (OUT/'README.md').write_text(policy)
    (OUT/'SCDEEPSORT_INPUT_AUDIT.md').write_text((CODE/'SCDEEPSORT_INPUT_AUDIT.md').read_text())
    track(CODE/'SCDEEPSORT_INPUT_AUDIT.md');track(__file__)
    writejson(OUT/'manifest.json',dict(status='audit_complete_not_comparison_complete',job=os.environ['SLURM_JOB_ID'],
        updated_at=utc(),n_datasets=13,n_public_datasets=11,n_HCL_tissues=59,n_analysis_units=72,
        n_dataset_method_entries=130,n_unit_method_entries=720,
        new_fits=0,expression_matrices_read=0,PTC_public_scientific_computation=0,
        sources=SOURCES,files={p.name:sha(p) for p in OUT.iterdir() if p.suffix in ['.csv','.md']}))
    writejson(OUT/'status.json',dict(stage='C_AUDIT',status='audit_complete_package_open',updated_at=utc(),jobs=[os.environ['SLURM_JOB_ID']],
        completed=['130 dataset-method and 720 unit-method information/missing-work entries','59 HCL tissue/donor metadata roster','Input and training-source overlap audit'],
        remaining=['GBM corrected-input scDeepSort pilot/replay','GBM gate before PTC/public comparator fits','Reference support/folds and eligible remaining runs; whole C is not closed'],
        evidence=[str(OUT/'frozen_dataset_method_matrix.csv'),str(OUT/'frozen_unit_method_matrix.csv'),str(OUT/'SCDEEPSORT_INPUT_AUDIT.md')]))
    print(json.dumps(dict(datasets=13,units=72,matrix_rows=len(unitmatrix),audit='complete'),indent=2))

if __name__=='__main__':run()
