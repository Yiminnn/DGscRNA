"""SLURM-only exact parity of packaged outputs against the native R/DL archive.

The verifier never opens expression-truth tables or imports the fitting package.
Marker-arm identities are (route, library, cutoff), not positional LXX IDs.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import os
import subprocess
import sys
import traceback

HERE = Path(__file__).resolve().parent
ROOT = Path('/fs/scratch/PCON0080/yimin/dgscrna')
DEFAULT_REFERENCE = ROOT/'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917'
ROUTES = ['PCA30_SNN', 'PCA30_HDBSCAN_R', 'UMAP2_SNN', 'UMAP2_HDBSCAN_R']


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load(path):
    return json.loads(Path(path).read_text())


def need(value, message):
    if not value:
        raise AssertionError(message)


def receipt(directory, manifest, flag):
    directory = Path(directory)
    need((directory/flag).read_text().strip() == sha(directory/manifest),
         f'Invalid completion receipt: {directory}/{flag}')


def write_new(path, value):
    path = Path(path)
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def byte_equal(left, right, label, checks):
    a, b = sha(left), sha(right)
    need(a == b, f'Byte mismatch: {label}')
    checks.append(dict(artifact=label, comparison='sha256_exact', sha256=a))


def frame(path):
    import pandas as pd
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def frame_equal(left, right, label, checks):
    import pandas as pd
    pd.testing.assert_frame_equal(left, right, check_exact=True)
    checks.append(dict(artifact=label, comparison='ordered_table_exact', rows=len(left)))


def arm_map(manifest):
    result = {}
    for aid, arm in manifest['arms'].items():
        key = (arm['library'], str(arm['cutoff']))
        need(key not in result, f'Duplicate semantic arm {key}')
        need(arm['arm_id'] == aid, f'Inconsistent arm ID: {aid}')
        result[key] = (aid, arm)
    return result


def compare_terminal(actual, expected, record, aarm, earm, checks, fresh):
    import numpy as np
    receipt(actual, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
    receipt(expected, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
    am, em = load(actual/'terminal_manifest.json'), load(expected/'terminal_manifest.json')
    for manifest in [am, em]:
        need(manifest['terminal_valid'] is True, 'Invalid terminal condition cannot pass parity')
        need(manifest['reference_labels_used_for_fit'] is False, 'Labels declared used for fit')
    need((am['arm']['library'], str(am['arm']['cutoff'])) == (record['library'], record['cutoff']),
         'Run condition and terminal manifest disagree')
    need((em['arm']['library'], str(em['arm']['cutoff'])) == (record['library'], record['cutoff']),
         'Reference semantic terminal arm mismatch')
    need(am['arm']['arm_id'] == aarm and em['arm']['arm_id'] == earm, 'Terminal internal arm ID mismatch')
    for key in ['dl_status', 'training_executed', 'n_known', 'n_pool', 'n_training_classes',
                'n_final_called090', 'n_final_called070', 'DL_sha256', 'DL_features',
                'known_labels_unchanged', 'threshold_reconstructed']:
        need(am[key] == em[key], f'Terminal metadata mismatch: {key}')
    need(record['dl_status'] == am['dl_status'] and record['training_executed'] == am['training_executed'],
         'Run condition reports wrong training status')
    if fresh:
        need(am['identical_result_reused'] is False, 'Requested fresh parity reused a terminal cache')
    for directory, manifest in [(actual, am), (expected, em)]:
        for name, field in [('terminal.npz','terminal_sha256'),
                            ('predictions.csv.gz','predictions_sha256'),
                            ('training_manifest.json','training_manifest_sha256')]:
            need(sha(directory/name) == manifest[field], f'Invalid terminal artifact: {directory/name}')
    arrays = []
    with np.load(actual/'terminal.npz', allow_pickle=False) as a, np.load(expected/'terminal.npz', allow_pickle=False) as b:
        need(set(a.files) == set(b.files), 'Different terminal NPZ fields')
        for key in a.files:
            need(a[key].dtype == b[key].dtype, f'Dtype mismatch: {key}')
            if np.issubdtype(a[key].dtype, np.inexact):
                np.testing.assert_allclose(a[key], b[key], rtol=0, atol=0, equal_nan=True, err_msg=key)
            else:
                np.testing.assert_array_equal(a[key], b[key], err_msg=key)
            arrays.append(key)
        initial, pool = a['initial'], a['pool_indices']
        known = np.flatnonzero(initial != 'Undecided')
        need(np.array_equal(pool, np.flatnonzero(initial == 'Undecided')), 'Incorrect unknown pool')
        for stage in ['final090', 'final070']:
            need(np.array_equal(a[stage][known], initial[known]), 'Known calls changed')
        if am['training_executed']:
            probs, classes = a['probabilities'], a['classes']
            need(probs.shape == (len(pool), len(classes)) and np.isfinite(probs).all(), 'Invalid probabilities')
            np.testing.assert_allclose(probs.sum(axis=1), 1, rtol=1e-5, atol=1e-6)
            rounded = np.asarray([round(v, 4) for v in probs.max(axis=1)], dtype=np.float32)
            for stage, threshold in [('final090', .9), ('final070', .7)]:
                expected_calls = np.where(rounded >= threshold, classes[probs.argmax(axis=1)], 'Unknown')
                need(np.array_equal(a[stage][pool], expected_calls), 'Threshold reconstruction failed')
            train, valid = a['train_indices'], a['validation_indices']
            need(not set(train).intersection(valid), 'Training and validation overlap')
            need(np.array_equal(np.sort(np.concatenate([train, valid])), known), 'Known split incomplete')
        else:
            need(am['dl_status'] in ['no_op_all_initially_known',
                 'no_known_labels_archived_Undecided_terminal', 'structural_insufficient_known_split'],
                 'Unexpected no-training status')
            need(a['probabilities'].shape[0] == 0 and len(a['train_indices']) == 0
                 and len(a['validation_indices']) == 0, 'No-op contains fabricated training')
    label = '/'.join([record['route'], record['library'], record['cutoff']])
    frame_equal(frame(actual/'predictions.csv.gz'), frame(expected/'predictions.csv.gz'), label+'/predictions', checks)
    ah, eh = actual/'training_history.json', expected/'training_history.json'
    need(ah.exists() == eh.exists() == bool(am['training_executed']), 'Training-history presence disagrees with state')
    if ah.exists():
        need(load(ah) == load(eh), f'Training history mismatch: {label}')
    at, et = load(actual/'training_manifest.json'), load(expected/'training_manifest.json')
    for key in ['n_cells','n_known','n_pool','n_training_classes','classes','input_width',
                'training_executed','prediction_executed','terminal_valid','dl_status',
                'n_train','n_validation','known_seed_validation_accuracy','n_final_called090','n_final_called070']:
        need(at.get(key) == et.get(key), f'Training metadata mismatch: {label}/{key}')
    # Portable descriptive wording is not part of the numerical training protocol.
    ignored = {'input', 'num_workers_note'}
    need({k:v for k,v in at['params'].items() if k not in ignored} ==
         {k:v for k,v in et['params'].items() if k not in ignored}, 'Numerical training parameters changed')
    for key in ['DL_sha256','n_cells','n_features','cell_order_sha256','torch_version','torch_threads']:
        need(at['provenance']['input_signature'][key] == et['provenance']['input_signature'][key],
             f'Training input/runtime signature mismatch: {key}')
    ap, ep = actual/'model_state.pt', expected/'model_state.pt'
    need(ap.exists() == ep.exists() == bool(am['training_executed']), 'Weight presence disagrees with training state')
    if ap.exists():
        import torch
        left, right = torch.load(ap,map_location='cpu',weights_only=True), torch.load(ep,map_location='cpu',weights_only=True)
        need(set(left) == set(right), 'Model parameter names changed')
        for key in left:
            need(left[key].dtype == right[key].dtype and torch.equal(left[key],right[key]), f'Model weight mismatch: {key}')
    return dict(route=record['route'], library=record['library'], cutoff=record['cutoff'],
                actual_arm=aarm, reference_arm=earm, status='passed_exact',
                dl_status=am['dl_status'], training_executed=am['training_executed'],
                reused=am['identical_result_reused'], exact_npz_fields=arrays,
                history_exact=ah.exists(), weights_exact=ap.exists(),
                descriptive_input_actual=at['params'].get('input'),
                descriptive_input_reference=et['params'].get('input'),
                actual_manifest_sha256=sha(actual/'terminal_manifest.json'),
                reference_manifest_sha256=sha(expected/'terminal_manifest.json'))


def verify(args):
    import pandas as pd
    output, reference = args.output.resolve(), args.reference_root.resolve()
    receipt(output, 'run_manifest.json', 'COMPLETE')
    manifest, plan = load(output/'run_manifest.json'), load(output/'run_config.json')
    need(manifest['status'] == 'completed', 'Package run is incomplete')
    need(manifest['reference_labels_used_for_fit'] is False, 'Run declares reference-label fitting')
    need(manifest['plan_sha256'] == sha(output/'run_config.json'), 'Run configuration hash mismatch')
    cfg = plan['config']; sample = cfg['sample']
    budget = 'all' if str(cfg['features']) == 'all' else 'hvg'+str(cfg['features']).removeprefix('hvg')
    condition = budget if cfg['seed'] == 42 else f'{budget}_seed{cfg["seed"]}'
    actual, expected = output/'GBM'/sample/condition, reference/'GBM'/sample/condition
    need(Path(manifest['prepare']).resolve() == actual, 'Run preparation path differs from configuration')
    need(expected.is_dir(), f'No existing reference: {expected}')
    for directory in [actual, expected]:
        receipt(directory, 'prepare_manifest.json', 'PREPARED')
    receipt(output/'inputs'/sample, 'input_manifest.json', 'INPUT_COMPLETE')
    checks = []
    # Only fitting files are opened: no labels.csv, obs.csv or evaluation inputs.
    for name in ['x.bin','i.bin','p.bin','cells_fit.csv','genes.csv']:
        byte_equal(output/'inputs'/sample/name, reference/'inputs'/sample/name, 'inputs/'+name, checks)
    for name in ['selected_features.txt','geometry_features.txt','scoring_features.txt','DL_features.txt','DL.float32.bin','cells.csv']:
        byte_equal(actual/name, expected/name, name, checks)
    for name in ['Seurat_gene_names.csv','native_vst_statistics.csv','PCA30.csv','UMAP2.csv']:
        need((actual/name).exists() == (expected/name).exists(), f'Preparation artifact presence: {name}')
        if (actual/name).exists():
            frame_equal(frame(actual/name), frame(expected/name), name, checks)
    ap, ep = load(actual/'prepare_manifest.json'), load(expected/'prepare_manifest.json')
    for key in ['n_cells','seed','requested_features','features','correction','assay','n_batches','DL_binary_sha256']:
        need(ap[key] == ep[key], f'Preparation semantics changed: {key}')
    need(ap['reference_labels_used_for_fitting'] is False, 'Preparation declares reference-label fitting')
    new_libraries, old_libraries = load(output/'markers/libraries.json'), load(reference/'markers/libraries.json')
    for library, panels in new_libraries.items():
        need(library in old_libraries and panels == old_libraries[library], f'Marker content differs: {library}')
    records = manifest['conditions']
    need(manifest['terminal_condition_count'] == len(records) and len(records)>0, 'Bad terminal condition count')
    keys = [(r['route'],r['library'],str(r['cutoff'])) for r in records]
    need(len(keys) == len(set(keys)), 'Duplicate terminal conditions')
    routes = ROUTES if cfg['route'] == 'all' else [cfg['route']]
    conditions, wanted, route_checks = [], set(), []
    for route in routes:
        a, e = actual/route, expected/route
        receipt(a,'score_manifest.json','SCORE_COMPLETE'); receipt(e,'score_manifest.json','SCORE_COMPLETE')
        am, em = load(a/'score_manifest.json'), load(e/'score_manifest.json')
        need(am['reference_labels_used_for_fit'] is False, 'Scorer declares reference-label fitting')
        aa, ea = arm_map(am), arm_map(em)
        need(set(aa) == {(library, cutoff) for library in new_libraries
                        for cutoff in ['none', 'mean', '0.5']}, 'Incomplete scoring arm roster')
        need(set(aa).issubset(ea), 'Package generated a scoring arm absent from the reference')
        for key in ['cells.csv','clusters.csv','density_diagnostics.csv']:
            need((a/key).exists() == (e/key).exists(), f'Clustering artifact presence: {route}/{key}')
            if (a/key).exists(): frame_equal(frame(a/key), frame(e/key), route+'/'+key, checks)
        ai, ei = frame(a/'initial_calls.csv.gz'), frame(e/'initial_calls.csv.gz')
        frame_equal(ai[['cell_id']],ei[['cell_id']],route+'/initial_cell_order',checks)
        ac, ec = frame(a/'cluster_calls.csv.gz'), frame(e/'cluster_calls.csv.gz')
        ar, er = frame(a/'marker_retention.csv.gz'), frame(e/'marker_retention.csv.gz')
        for key, (aid, arm) in aa.items():
            eid, old = ea[key]
            frame_equal(ai[[arm['seed_column']]].set_axis(['initial'],axis=1),
                        ei[[old['seed_column']]].set_axis(['initial'],axis=1),
                        route+'/'+str(key)+'/initial',checks)
            left = ac[(ac.library==key[0])&(ac.cutoff==key[1])].drop(columns='arm_id').reset_index(drop=True)
            right = ec[(ec.library==key[0])&(ec.cutoff==key[1])].drop(columns='arm_id').reset_index(drop=True)
            frame_equal(left,right,route+'/'+str(key)+'/cluster_calls',checks)
            if (cfg['library']=='all' or key[0]==cfg['library']) and (cfg['cutoff']=='all' or key[1]==cfg['cutoff']):
                wanted.add((route,*key))
        for library in new_libraries:
            frame_equal(ar[ar.library==library].reset_index(drop=True),er[er.library==library].reset_index(drop=True),
                        route+'/'+library+'/marker_retention',checks)
        # The original helper's identical(DEG,density) check is generalized only
        # to compare named libraries, so harmless library-order changes are legal.
        libraries_path = output/(args.proof_name+'.'+route+'.libraries.json')
        write_new(libraries_path,list(new_libraries))
        runtime=plan['runtime']; env=os.environ.copy()
        for name in ['R_LIBS','R_LIBS_USER','R_LIBS_SITE','DGSCRNA_REFERENCE_R_LIB']:
            env.pop(name,None)
        if runtime.get('reference_r_lib'): env['DGSCRNA_REFERENCE_R_LIB']=runtime['reference_r_lib']
        command=[runtime['rscript'],'--vanilla',str(args.r_helper),str(a),str(e),str(libraries_path)]
        completed=subprocess.run(command,env=env,text=True,capture_output=True)
        log_path=output/(args.proof_name+'.'+route+'.R.log')
        with log_path.open('x') as log:log.write(completed.stdout+'\n'+completed.stderr)
        need(completed.returncode==0 and 'R_DEG_AND_NAMED_DENSITY_OBJECTS_EXACT' in completed.stdout,
             f'R object parity failed; see {log_path}')
        route_checks.append(dict(route=route,scored_arms=len(aa),DEG_exact=True,density_exact_libraries=list(new_libraries),
                                 R_log_sha256=sha(log_path)))
        for record in records:
            if record['route'] != route: continue
            key=(record['library'],str(record['cutoff']))
            need(key in aa and key in ea, 'Unscored run condition')
            aid, _ = aa[key]; eid, _ = ea[key]
            need(record['arm']==aid, 'Run condition arm ID inconsistent')
            td=a/'terminal'/aid
            need((output/record['predictions']).resolve()==(td/'predictions.csv.gz').resolve(),'Prediction path not selected terminal')
            need(sha(td/'predictions.csv.gz')==record['predictions_sha256'],'Run prediction hash mismatch')
            conditions.append(compare_terminal(td,e/'terminal'/eid,record,aid,eid,checks,not args.allow_fresh_cache_reuse))
    need(set(keys)==wanted, 'Run manifest omits or adds a requested terminal condition')
    need(manifest['trained_condition_count']==sum(r['training_executed'] for r in conditions),'Wrong trained condition count')
    for name, expected_hash in manifest['exports'].items():
        need(sha(output/name)==expected_hash,f'Export hash mismatch: {name}')
    return dict(status='passed_exact',sample=sample,budget=budget,condition=condition,
                output=str(output),reference=str(expected),terminal_conditions=len(conditions),
                trained_conditions=sum(r['training_executed'] for r in conditions),
                valid_no_training_conditions=sum(not r['training_executed'] for r in conditions),
                require_fresh_terminal=not args.allow_fresh_cache_reuse,
                no_reference_truth_read=True,route_checks=route_checks,conditions=conditions,checks=checks,
                run_manifest_sha256=sha(output/'run_manifest.json'),run_config_sha256=sha(output/'run_config.json'),
                reference_prepare_manifest_sha256=sha(expected/'prepare_manifest.json'),
                verifier_sha256=sha(__file__),R_helper_sha256=sha(args.r_helper),
                numerical_tolerance=dict(rtol=0,atol=0),
                job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),
                completed_at=datetime.now(timezone.utc).isoformat())


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--reference-root',type=Path,default=DEFAULT_REFERENCE)
    parser.add_argument('--proof-name',default='packaged_parity.json')
    parser.add_argument('--r-helper',type=Path,default=HERE/'verify_packaged_r_artifacts.R')
    parser.add_argument('--allow-fresh-cache-reuse',action='store_true',help='Permit within-new-run identical-input reuse; never required for valid no-op')
    args=parser.parse_args()
    need(os.environ.get('SLURM_JOB_ID'), 'Scientific verification requires SLURM')
    need(not sys.flags.optimize,'Do not disable assertions')
    need(Path(args.proof_name).name==args.proof_name,'Proof name must be a filename')
    proof=args.output/args.proof_name
    need(not proof.exists(),'Proof already exists; use a new proof name')
    try:
        result=verify(args)
    except Exception as error:
        result=dict(status='failed',error=str(error),traceback=traceback.format_exc(),
                    output=str(args.output),job=os.environ['SLURM_JOB_ID'],verifier_sha256=sha(__file__),
                    timestamp=datetime.now(timezone.utc).isoformat())
        write_new(proof,result)
        raise
    write_new(proof,result)
    print(json.dumps({k:result[k] for k in ['status','sample','budget','terminal_conditions','trained_conditions','valid_no_training_conditions']}))


if __name__=='__main__':
    main()
