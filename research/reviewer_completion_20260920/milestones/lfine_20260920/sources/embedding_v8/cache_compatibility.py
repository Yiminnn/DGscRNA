"""Explicit producer whitelist; consume old success without relabelling it."""
from pathlib import Path
from datetime import datetime,timezone
import json
import os


def catalog(core):
    path=core.CODE/'CACHE_COMPATIBILITY.json'
    data=json.loads(path.read_text())
    for version,item in data['approved_previous_sources'].items():
        directory=core.ROOT/item['directory']
        assert core.sha(directory/'SOURCE_MANIFEST.json')==item['manifest_sha256']
        manifest=json.loads((directory/'SOURCE_MANIFEST.json').read_text())
        assert manifest['run.py']==item['run_sha256']
        for name in data['unchanged_dependency_files']:
            assert core.sha(core.CODE/name)==manifest[name]==core.sha(directory/name), (version,name)
    return data,path


def accepted_producer(core,producer):
    data,path=catalog(core)
    allowed={v['run_sha256'] for v in data['approved_previous_sources'].values()}
    assert producer in allowed|{core.sha(core.CODE/'run.py')}, 'Unapproved cache producer'
    return data,path


def record_consumer(core,cfg,kind,source_path,producer):
    _,catalog_path=accepted_producer(core,producer)
    from adaptive_ica import policy
    _,policy_path=policy(core)
    directory=Path(cfg['dest'])/'cache_consumers_v8'
    directory.mkdir(exist_ok=True)
    snapshots=directory/'producer_snapshots';snapshots.mkdir(exist_ok=True)
    original_digest=core.sha(source_path)
    snapshot=snapshots/(original_digest+'.json')
    if not snapshot.exists():snapshot.write_bytes(Path(source_path).read_bytes())
    assert core.sha(snapshot)==original_digest
    suffix=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%f')
    path=directory/(str(os.environ.get('SLURM_JOB_ID','unknown'))+'_'+kind+'_'+suffix+'.json')
    core.write_json(path,dict(status='compatible_reuse',kind=kind,
        producer_source_sha256=producer,producer_manifest=str(snapshot),
        original_producer_manifest_path=str(source_path),
        producer_manifest_sha256=original_digest,consumer_source_sha256=core.sha(core.CODE/'run.py'),
        compatibility_catalog_sha256=core.sha(catalog_path),policy_sha256=core.sha(policy_path),
        producer_metadata_unchanged=True,reference_labels_used=False))


def representation(core,cfg,metadata):
    from adaptive_ica import policy,POLICY_ID
    producer=metadata['source_sha256']
    assert core.checked(Path(cfg['dest']),'representation.json','REPRESENTATION_COMPLETE')
    assert metadata['protocol_sha256']==cfg['frozen_protocol_sha256']
    assert all(metadata[k]==cfg[k] for k in ('sample','budget','space'))
    geometry=Path(cfg['geometry'])
    assert core.checked(geometry)
    gm=json.loads((geometry/'manifest.json').read_text())
    assert metadata['geometry_binary_sha256']==gm['binary_sha256']
    assert metadata['n_cells']==gm['n_cells']
    assert core.sha(geometry/'cells.csv')==gm['cells_sha256'], 'Changed geometry cell order file'
    # FIT_COMPLETE reuse does not call make_embedding: validate the actual
    # representation here, not only its signed metadata.
    import csv
    with (geometry/'cells.csv').open(newline='') as stream:
        cells=[row['cell_id'] for row in csv.DictReader(stream)]
    assert len(cells)==gm['n_cells'] and len(set(cells))==len(cells)
    if cfg['space']!='noDR':
        embedding=Path(cfg['embedding'])
        assert embedding.is_file() and core.sha(embedding)==metadata['embedding_sha256'], 'Changed embedding CSV'
        with embedding.open(newline='') as stream:
            reader=csv.DictReader(stream)
            assert reader.fieldnames==['cell_id','x','y']
            embedded_cells=[row['cell_id'] for row in reader]
        assert embedded_cells==cells, 'Embedding cell order differs from geometry'
    else:
        assert metadata['embedding_sha256'] is None
    accepted_producer(core,producer)
    if cfg['space']=='ICA2':
        rules,policy_path=policy(core)
        if producer==core.sha(core.CODE/'run.py'):
            assert metadata['comparator_policy_id']==POLICY_ID
            assert metadata['policy_sha256']==core.sha(policy_path)
            assert metadata['params']['algorithm']==metadata['actual_solver']
            assert {'algorithm':metadata['actual_solver'],'max_iter':metadata['actual_iteration_cap']} in rules['attempts']
            for attempt in metadata['attempts']:
                assert core.sha(attempt['path'])==attempt['sha256']
            assert metadata['attempts'][-1]['accepted'] is True
        else:
            expected=dict(n_components=2,max_iter=5000,tol=1e-4,whiten='unit-variance',
                random_state=42,algorithm='parallel',fun='logcosh',whiten_solver='svd',w_init=None)
            assert all(metadata['params'].get(k)==v for k,v in expected.items())
            assert 0<metadata['n_iter_']<=5000
            assert metadata['sklearn_version']=='1.9.0'
            assert not any(w.startswith('ConvergenceWarning:') for w in metadata['warnings'])
    record_consumer(core,cfg,'representation',Path(cfg['dest'])/'representation.json',producer)


def own_source(core):
    manifest=json.loads((core.CODE/'SOURCE_MANIFEST.json').read_text())
    for name,digest in manifest.items():
        assert core.sha(core.CODE/name)==digest, f'Changed frozen v8 source: {name}'
    catalog(core)


def validate_completed_fit(core,cfg,previous):
    """Validate the complete consumed chain; never infer success from flags alone."""
    dest=Path(cfg['dest']);prep=Path(cfg['prep'])
    accepted_producer(core,previous['source_sha256'])
    representation(core,cfg,json.loads((dest/'representation.json').read_text()))
    core.get_geometry(cfg)
    for key in ('sample','budget','space','prep','dest','geometry','embedding',
                'frozen_protocol_sha256','conditions'):
        assert previous['config'][key]==cfg[key], f'Changed fit configuration: {key}'
    expected_conditions=[dict(route=cfg['space']+'_'+(method if k is None else f'{method}_K{k:02d}'),
        dest=str(dest/(method if k is None else f'{method}_K{k:02d}')),method=method,k=k)
        for method in ('KMeans','GMM','HDBSCAN_R') for k in (core.KS if method!='HDBSCAN_R' else [None])]
    assert cfg['conditions']==expected_conditions, 'Incomplete or changed candidate roster'
    assert set(previous['terminal_manifest_hashes'])=={c['route'] for c in expected_conditions}
    assert core.checked(prep,'prepare_manifest.json','PREPARED')
    prepared=json.loads((prep/'prepare_manifest.json').read_text())
    assert core.sha(prep/'expression_PCA30.rds')==prepared['expression_sha256'], 'Changed scoring expression'
    assert core.sha(prepared['DL_binary'])==prepared['DL_binary_sha256'], 'Changed DL input'
    marker=core.REFERENCE/'markers/libraries.json'
    expected_shared=dict(prepare_manifest_sha256=core.sha(prep/'prepare_manifest.json'),
        cells_sha256=core.sha(prep/'cells.csv'),expression_sha256=prepared['expression_sha256'],
        marker_sha256=core.sha(marker),scorer_sha256=core.sha(core.CODE/'score_candidates.R'),
        protocol_sha256=cfg['frozen_protocol_sha256'])
    import csv
    with (prep/'cells.csv').open(newline='') as stream:
        prepared_cells=[row['cell_id'] for row in csv.DictReader(stream)]
    with (Path(cfg['geometry'])/'cells.csv').open(newline='') as stream:
        geometry_cells=[row['cell_id'] for row in csv.DictReader(stream)]
    assert prepared_cells==geometry_cells, 'Scorer cell order differs from geometry'
    for condition in expected_conditions:
        path=Path(condition['dest'])
        assert core.checked(path,'partition_manifest.json','PARTITION_COMPLETE')
        pm=json.loads((path/'partition_manifest.json').read_text())
        assert pm['protocol_sha256']==cfg['frozen_protocol_sha256']
        assert pm['representation_manifest_sha256']==core.sha(dest/'representation.json')
        assert pm['clusters_sha256']==core.sha(path/'clusters.csv')
        with (path/'clusters.csv').open(newline='') as stream:
            cluster_cells=[row['cell_id'] for row in csv.DictReader(stream)]
        assert cluster_cells==prepared_cells, 'Partition cell order differs from scorer'
        fingerprint=json.loads((path/'score_input_fingerprint.json').read_text())
        expected=dict(clusters_sha256=pm['clusters_sha256'],**expected_shared)
        assert fingerprint==expected, 'Changed complete score input fingerprint'
        assert core.checked(path,'score_manifest.json','SCORE_COMPLETE')
        sm=json.loads((path/'score_manifest.json').read_text())
        assert sm['input_fingerprint_sha256']==core.sha(path/'score_input_fingerprint.json')
        assert sm['initial_sha256']==core.sha(path/'initial_calls.csv.gz')
        assert sm['marker_source_sha256']==expected_shared['marker_sha256']
        assert sm['source_sha256']==expected_shared['scorer_sha256']
        assert sm['DL_binary_sha256']==prepared['DL_binary_sha256']
        td=path/'terminal/L00_mean'
        assert core.checked(td,'terminal_manifest.json','TERMINAL_COMPLETE')
        assert core.sha(td/'terminal_manifest.json')==previous['terminal_manifest_hashes'][condition['route']]
        tm=json.loads((td/'terminal_manifest.json').read_text())
        assert tm['terminal_valid'] is True
        assert tm['score_manifest_sha256']==core.sha(path/'score_manifest.json')
        assert tm['DL_sha256']==prepared['DL_binary_sha256']
        assert tm['terminal_sha256']==core.sha(td/'terminal.npz')
        assert tm['predictions_sha256']==core.sha(td/'predictions.csv.gz'), 'Changed terminal prediction CSV'
    record_consumer(core,cfg,'fit',dest/'fit_manifest.json',previous['source_sha256'])
