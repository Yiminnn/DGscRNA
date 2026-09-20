"""Independent read-only audit of the first A1 representation output."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os

ROOT=Path('/fs/scratch/PCON0080/yimin/dgscrna')
CAMPAIGN=ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920'
UNIT=CAMPAIGN/'embedding/TKU4163/hvg2000/PCA2'
OUT=CAMPAIGN/'no_clustering/a1_first_review'
REFERENCE=ROOT/'results/hvg_ptc_20260916_v1/paper_claim_validation_20260917'
L1=['Malignant','TAM','Lymphocyte','Oligodendrocyte','Astrocyte','OPC','Excitatory neuron','Inhibitory neuron','Endothel','Pericyte','Other']


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def complete(directory, name='manifest.json', marker='COMPLETE'):
    assert (directory/marker).read_text().strip()==sha(directory/name)
    return read(directory/name)


def run():
    assert os.environ.get('SLURM_JOB_ID')
    import numpy as np
    import pandas as pd
    from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, fowlkes_mallows_score
    cfg=read(UNIT/'config.json');prep=Path(cfg['prep'])
    fit=complete(UNIT,'fit_manifest.json','FIT_COMPLETE')
    representation=read(UNIT/'representation.json')
    assert representation['embedding_sha256']==sha(UNIT/'embedding.csv')
    assert representation['protocol_sha256']==sha(CAMPAIGN/'protocol/embedding.json')
    assert fit['reference_labels_used_for_fit'] is False
    em=complete(UNIT/'evaluation');fm=complete(UNIT/'figures')
    for name,digest in em['outputs'].items():assert sha(UNIT/'evaluation'/name)==digest
    for name,digest in fm['files'].items():assert sha(UNIT/'figures'/name)==digest
    assert len(fm['files'])==26 and fm['every_candidate_plotted'] is True
    cells=pd.read_csv(prep/'cells.csv',dtype=str,keep_default_na=False)
    truth=pd.read_csv(REFERENCE/'evaluation_inputs/TKU4163/truth.csv.gz',dtype=str,keep_default_na=False)
    assert sha(UNIT/'cells.csv')==sha(prep/'cells.csv')
    assert len(cells)==178 and cells.cell_id.equals(truth.cell_id)
    assert em['truth_sha256']==sha(REFERENCE/'evaluation_inputs/TKU4163/truth.csv.gz')
    assert fm['display_sha256']==sha(Path(fm['display']))
    mpath=REFERENCE/'markers/panel_L1_mapping.csv'
    assert em['mapping_sha256']==sha(mpath)
    mapping=pd.read_csv(mpath,dtype=str,keep_default_na=False)
    mapping=mapping[mapping.library=='CM2_glioma_other'];lookup=dict(zip(mapping.panel,mapping.L1))
    metrics=pd.read_csv(UNIT/'evaluation/metrics.csv')
    cluster_metrics=pd.read_csv(UNIT/'evaluation/clustering.csv')
    expected={f'{method}_K{k:02}' for method in ['KMeans','GMM'] for k in [5,10,15,20,30,40]}|{'HDBSCAN_R'}
    assert len(cfg['conditions'])==13 and {Path(c['dest']).name for c in cfg['conditions']}==expected
    assert len(metrics)==39 and len(cluster_metrics)==13
    assert set(metrics.route)==set(cluster_metrics.route)==expected
    assert metrics.n_cells.eq(178).all() and cluster_metrics.n_cells.eq(178).all()
    checks=[];y=truth.L1.to_numpy()
    for condition in cfg['conditions']:
        route=Path(condition['dest']);name=route.name
        partition=complete(route,'partition_manifest.json','PARTITION_COMPLETE')
        assert partition['clusters_sha256']==sha(route/'clusters.csv')
        assert partition['representation_manifest_sha256']==sha(UNIT/'representation.json')
        assert partition['protocol_sha256']==cfg['frozen_protocol_sha256']
        assert partition['reference_labels_used_for_fit'] is False
        fingerprint=read(route/'score_input_fingerprint.json')
        for field,path in [('clusters_sha256',route/'clusters.csv'),('prepare_manifest_sha256',prep/'prepare_manifest.json'),
                           ('cells_sha256',prep/'cells.csv'),('expression_sha256',prep/'expression_PCA30.rds'),
                           ('marker_sha256',REFERENCE/'markers/libraries.json'),('protocol_sha256',CAMPAIGN/'protocol/embedding.json')]:
            assert fingerprint[field]==sha(path),(name,field)
        score=complete(route,'score_manifest.json','SCORE_COMPLETE')
        assert fingerprint['scorer_sha256']==score['source_sha256']==sha(Path(score['execution_source']))
        assert score['input_fingerprint_sha256']==sha(route/'score_input_fingerprint.json')
        assert score['initial_sha256']==sha(route/'initial_calls.csv.gz')
        assert score['DL_binary_sha256']==sha(Path(score['DL_binary']))
        terminal=route/'terminal/L00_mean';tm=complete(terminal,'terminal_manifest.json','TERMINAL_COMPLETE')
        assert tm['score_manifest_sha256']==sha(route/'score_manifest.json')
        assert tm['terminal_sha256']==sha(terminal/'terminal.npz')
        assert tm['predictions_sha256']==sha(terminal/'predictions.csv.gz')
        assert tm['DL_sha256']==score['DL_binary_sha256']
        assert tm['reference_labels_used_for_fit'] is False
        assert fit['terminal_manifest_hashes'][condition['route']]==sha(terminal/'terminal_manifest.json')
        assert em['terminal_manifests'][name+'/terminal/L00_mean']==sha(terminal/'terminal_manifest.json')
        cluster=pd.read_csv(route/'clusters.csv',dtype=str,keep_default_na=False)
        pred=pd.read_csv(terminal/'predictions.csv.gz',dtype=str,keep_default_na=False)
        initial=pd.read_csv(route/'initial_calls.csv.gz',dtype=str,keep_default_na=False)
        assert cluster.cell_id.equals(cells.cell_id) and pred.cell_id.equals(cells.cell_id) and initial.cell_id.equals(cells.cell_id)
        assert len(pred)==len(cluster)==178
        with np.load(terminal/'terminal.npz') as z:
            assert np.array_equal(z['initial'],initial.L00_mean)
            for key in ['initial','final090','final070']:assert np.array_equal(z[key],pred[key])
            known=z['initial']!='Undecided'
            assert np.array_equal(z['initial'][known],z['final090'][known])
        max_delta=0.0
        for stage,key in [('marker_only','initial'),('terminal090','final090'),('terminal070','final070')]:
            native=pred[key].to_numpy()
            mapped=np.asarray([lookup.get(v,'Unknown' if v in ['Unknown','Undecided','Noise',''] else 'UNMAPPABLE') for v in native])
            f1=[];supports=[]
            for label in L1:
                tp=int(((y==label)&(mapped==label)).sum())
                fp=int(((y!=label)&(mapped==label)).sum())
                fn=int(((y==label)&(mapped!=label)).sum())
                denom=2*tp+fp+fn
                f1.append(2*tp/denom if denom else 0.0);supports.append(int((y==label).sum()))
            f1=np.asarray(f1);supports=np.asarray(supports)
            abstain=np.isin(native,['Unknown','Undecided','Noise',''])
            recalculated=dict(macroF1_present=f1[supports>0].mean(),macroF1_fixed11=f1.mean(),
                              weightedF1=(f1*supports).sum()/supports.sum(),accuracy=(mapped==y).mean(),
                              coverage=(~abstain).mean(),unknown_rate=abstain.mean(),mapped_coverage=np.isin(mapped,L1).mean())
            recorded=metrics[(metrics.route==name)&(metrics.stage==stage)]
            assert len(recorded)==1
            for field,value in recalculated.items():
                delta=abs(float(recorded.iloc[0][field])-float(value));max_delta=max(max_delta,delta)
                assert delta<1e-12,(name,stage,field,delta)
        cr=cluster_metrics[cluster_metrics.route==name].iloc[0]
        for field,value in [('L1_ARI',adjusted_rand_score(y,cluster.cluster)),
                            ('L1_NMI',normalized_mutual_info_score(y,cluster.cluster)),
                            ('L1_FMI',fowlkes_mallows_score(y,cluster.cluster))]:
            assert abs(float(cr[field])-value)<1e-12,(name,field)
        checks.append(dict(route=name,n_cells=178,dl_status=tm['dl_status'],training_executed=tm['training_executed'],
                           n_known=tm['n_known'],n_pool=tm['n_pool'],max_metric_absolute_difference=max_delta))
    report=dict(status='passed',unit=str(UNIT),n_conditions=13,n_metric_rows=39,n_cells_each=178,
                all_partition_score_terminal_hash_chains_valid=True,all_cells_and_order_preserved=True,
                independent_manual_F1_accuracy_coverage_checks=True,independent_ARI_NMI_FMI_checks=True,
                conditions=checks,figure_review=dict(trained='KMeans_K10.png',noop='HDBSCAN_R.png',
                    result='Both inspected. Truth/partition/terminal panels legible; shared display and fitted representation distinguished; no-op explicitly labelled DL not executed, all Unknown, coverage 0.'),
                fit_manifest_sha256=sha(UNIT/'fit_manifest.json'),evaluation_manifest_sha256=sha(UNIT/'evaluation/manifest.json'),
                figures_manifest_sha256=sha(UNIT/'figures/manifest.json'),validation_source_sha256=sha(Path(__file__)),
                scope='First TKU4163/HVG2000/PCA2 output only; no other A1 representation or large sample validated here',
                job=os.environ['SLURM_JOB_ID'],step=os.environ.get('SLURM_STEP_ID'),completed_at=datetime.now(timezone.utc).isoformat())
    OUT.mkdir(parents=True,exist_ok=True);(OUT/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


if __name__=='__main__':run()
