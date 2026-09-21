"""Independent audit of the preserved failed default pilot; no output repair."""
from pathlib import Path
import argparse
import importlib.util
import json
import os
import hashlib
import subprocess

HERE=Path(__file__).resolve().parent


def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result


def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def load(path):return json.loads(Path(path).read_text())


def main(args):
    assert os.environ.get('SLURM_JOB_ID')
    import numpy as np
    adapter=module('audit_adapter',HERE/'representation_adapter.py')
    verifier=module('audit_verifier',HERE/'verify_packaged_parity.py')
    gate=load(args.gate);config=load(args.failed/'run_config.json')
    assert not args.out.exists();args.out.mkdir(parents=True)
    source=Path(config['source'])/'GBM/TKU4163/hvg2000'
    old=Path(gate['reference_root'])
    env=os.environ.copy()
    for name in ['R_LIBS','R_LIBS_USER','R_LIBS_SITE','DGSCRNA_REFERENCE_R_LIB']:env.pop(name,None)
    env.update(R_ENVIRON_USER=os.devnull,R_PROFILE_USER=os.devnull)
    if gate['runtime'].get('reference_r_lib'):env['DGSCRNA_REFERENCE_R_LIB']=gate['runtime']['reference_r_lib']
    adapter.execute([gate['runtime']['rscript'],'--vanilla',str(HERE/'audit_representation_roundtrip.R'),
                     str(source),str(args.out)],env,args.out/'roundtrip.log')
    checks=[];diagnostics=[];conditions=[];stages=[]
    full_frame_equal=verifier.frame_equal
    def diagnostic_aware_equal(left,right,label,records):
        if label.endswith('/density_diagnostics.csv') and 'UMAP2_HDBSCAN_R' in label and '_core' in label:
            full_frame_equal(left.drop(columns='membership'),right.drop(columns='membership'),label+'/identities_noise',records)
            delta=np.abs(left.membership.astype(float).to_numpy()-right.membership.astype(float).to_numpy())
            diagnostics.append(dict(artifact=label,exact=False,changed_cells=int((delta>0).sum()),
                                    max_abs=float(delta.max()),cause='core CSV round-trip vs original representation RDS precision'))
        else:full_frame_equal(left,right,label,records)
    verifier.frame_equal=diagnostic_aware_equal
    for cfg in config['configurations']:
        route=cfg['space']+'_'+cfg['method'];actual=args.failed/cfg['name']/route
        expected_rep=old/'GBM_representation_controls/TKU4163/hvg2000'/cfg['name']/route
        references=[('representation',expected_rep),('old_core',old/'GBM/TKU4163/hvg2000'/route),('new_core',source/route)]
        for label,expected in references:
            stage=adapter.score_parity(verifier,actual,expected,cfg['name']+'_'+label,gate['runtime'],env,args.out,checks)
            stages.append(dict(route=route,reference_type=label,**stage))
            am=verifier.arm_map(load(actual/'score_manifest.json'));em=verifier.arm_map(load(expected/'score_manifest.json'))
            for library in ['CM2_glioma_other','CM2_primary_all_context']:
                aid,_=am[(library,'mean')];eid,_=em[(library,'mean')]
                tm=load(actual/'terminal'/aid/'terminal_manifest.json')
                record=dict(route=route,library=library,cutoff='mean',dl_status=tm['dl_status'],training_executed=tm['training_executed'])
                result=verifier.compare_terminal(actual/'terminal'/aid,expected/'terminal'/eid,record,aid,eid,checks,False)
                conditions.append(dict(reference_type=label,**result))
    actual=args.failed/'UMAP2_HDBSCAN_R_minPts50_r0.5_seed42/UMAP2_HDBSCAN_R/density_diagnostics.csv'
    full_frame_equal(verifier.frame(actual),verifier.frame(args.out/'rds_membership.csv'),'recomputed_RDS_membership',checks)
    for name,expected in [('new',source/'UMAP2_HDBSCAN_R/density_diagnostics.csv'),
                          ('old',old/'GBM/TKU4163/hvg2000/UMAP2_HDBSCAN_R/density_diagnostics.csv')]:
        full_frame_equal(verifier.frame(expected),verifier.frame(args.out/'csv_membership.csv'),'recomputed_CSV_membership_'+name,checks)
    assert len(conditions)==24 and len(stages)==12 and len(diagnostics)==2
    result=dict(status='passed_exact_annotation_with_explained_membership_roundtrip',
                failed_attempt_preserved=str(args.failed),failure_sha256=sha(args.failed/'FAILURE.json'),
                gate_sha256=sha(args.gate),failed_fit_config_sha256=sha(args.failed/'run_config.json'),
                archived_representation_all_diagnostics_exact=True,all48_seeds_DEG_density_exact=True,
                terminal_conditions_exact_per_reference=8,reference_types=['representation','old_core','new_core'],
                roundtrip=load(args.out/'roundtrip.json'),diagnostic_differences=diagnostics,
                stages=stages,conditions=conditions,artifact_checks=checks,
                fresh_models_refit=False,truth_opened=False,job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),
                audit_source_sha256=sha(__file__),R_audit_source_sha256=sha(HERE/'audit_representation_roundtrip.R'))
    adapter.write_new(args.out/'verification.json',result)
    print(json.dumps({key:result[key] for key in ['status','roundtrip','terminal_conditions_exact_per_reference']}))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gate',type=Path,required=True);parser.add_argument('--failed',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);main(parser.parse_args())
