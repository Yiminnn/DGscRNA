"""Verify the delivered core archives contain every terminal result and cluster plot.

This checks saved bytes and archive links, without fitting or recomputing metrics.
Large archive reads and notebook inspection run exclusively inside SLURM.
"""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
import json
import os
from pathlib import PurePosixPath
import tarfile

from common import ROOT, OUT, FEATURES, ROUTES, require_slurm, sha, write_json, complete, utc
from delivery import STAGE


def safe_name(name):
    path = PurePosixPath(name)
    assert name and not path.is_absolute() and '..' not in path.parts, name
    assert str(path) == name, name
    return name


def inspect_archive(task):
    sample, archive_record, n_cells = task
    path = STAGE / archive_record['path']
    assert path.stat().st_size == archive_record['bytes'], str(path)
    assert sha(path) == archive_record['sha256'], str(path)
    digests, objects, flags = {}, {}, {}
    links = 0
    # Streaming extraction into hashes/small metadata only; no files are extracted.
    with tarfile.open(path, 'r|gz') as archive:
        for member in archive:
            name = safe_name(member.name)
            assert name not in digests, (sample, 'duplicate member', name)
            if member.islnk():
                target = safe_name(member.linkname)
                assert target in digests, (sample, 'missing earlier hardlink target', name, target)
                digests[name] = digests[target]
                if target in objects:
                    objects[name] = objects[target]
                if target in flags:
                    flags[name] = flags[target]
                links += 1
                continue
            assert member.isfile(), (sample, 'nonregular/external archive dependency', name, member.type)
            payload = archive.extractfile(member)
            assert payload is not None
            digest = hashlib.sha256()
            small = name.endswith('.json') or PurePosixPath(name).name in {
                'COMPLETE', 'SCORE_COMPLETE', 'TERMINAL_COMPLETE'}
            chunks = [] if small else None
            for block in iter(lambda: payload.read(1024 * 1024), b''):
                digest.update(block)
                if small:
                    chunks.append(block)
            digests[name] = digest.hexdigest()
            if name.endswith('.json'):
                objects[name] = json.loads(b''.join(chunks))
            elif small:
                flags[name] = b''.join(chunks).decode().strip()

    def verified_manifest(prefix, filename='manifest.json', flag='COMPLETE'):
        key = prefix + '/' + filename
        assert flags[prefix + '/' + flag] == digests[key], (sample, key)
        return objects[key]

    statuses = Counter()
    terminal_count = 0
    expected_arms = {f'L{i:02d}_{c}' for i in range(16) for c in ['none', 'mean', 'p050']}
    for budget in FEATURES:
        evaluation = verified_manifest(budget + '/evaluation')
        assert evaluation['status'] == 'completed' and not evaluation['missing']
        assert evaluation['sample'] == sample and evaluation['budget'] == budget
        assert evaluation['n_metric_rows'] == 576
        for filename, expected in evaluation['outputs'].items():
            assert digests[f'{budget}/evaluation/{filename}'] == expected
        expected_terminal_paths = {
            f'{route}/terminal/{arm}' for route in ROUTES for arm in expected_arms}
        assert set(evaluation['terminal_manifests']) == expected_terminal_paths
        for route in ROUTES:
            prefix = budget + '/' + route
            score = verified_manifest(prefix, 'score_manifest.json', 'SCORE_COMPLETE')
            assert score['unit'] == sample and score['n_cells'] == n_cells
            assert set(score['arms']) == expected_arms
            assert digests[prefix + '/initial_calls.csv.gz'] == score['initial_sha256']
            for filename in ['cells.csv', 'clusters.csv']:
                assert prefix + '/' + filename in digests
            for arm, definition in score['arms'].items():
                terminal_prefix = prefix + '/terminal/' + arm
                terminal = verified_manifest(terminal_prefix, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
                terminal_hash = digests[terminal_prefix + '/terminal_manifest.json']
                assert evaluation['terminal_manifests'][route + '/terminal/' + arm] == terminal_hash
                assert terminal['status'] == 'completed' and terminal['terminal_valid'] is True
                assert terminal['arm'] == definition
                assert terminal['reference_labels_used_for_fit'] is False
                assert terminal['known_labels_unchanged'] is True
                assert terminal['threshold_reconstructed'] is True
                assert terminal['score_manifest_sha256'] == digests[prefix + '/score_manifest.json']
                assert terminal['predictions_sha256'] == digests[terminal_prefix + '/predictions.csv.gz']
                assert terminal['training_manifest_sha256'] == digests[terminal_prefix + '/training_manifest.json']
                training = objects[terminal_prefix + '/training_manifest.json']
                assert training['n_cells'] == n_cells and training['terminal_valid'] is True
                assert training['training_executed'] == terminal['training_executed']
                assert training['dl_status'] == terminal['dl_status']
                if 'training_history.json' in training['outputs']:
                    assert training['outputs']['training_history.json'] == digests[terminal_prefix + '/training_history.json']
                statuses[terminal['dl_status']] += 1
                terminal_count += 1
    assert terminal_count == 1152
    return dict(sample=sample, archive_sha256=archive_record['sha256'],
                archive_members=len(digests), internal_hardlinks_resolved=links,
                terminal_predictions_verified=terminal_count, dl_status_counts=dict(statuses))


def run():
    require_slurm()
    manifest_path = OUT / 'summary/core_delivery_manifest.json'
    delivery = json.loads(manifest_path.read_text())
    receipt_path = OUT / 'summary/DELIVERY_RECEIPT.json'
    receipt = json.loads(receipt_path.read_text())
    remote_receipt = json.loads((OUT / 'summary/REMOTE_RECEIPT_UPLOADED.json').read_text())
    assert receipt['status'] == 'delivered_and_verified'
    assert receipt['manifest_sha256'] == sha(manifest_path)
    assert remote_receipt['receipt_sha256'] == sha(receipt_path)
    assert remote_receipt['status'] == 'remote_receipt_uploaded_and_verified'
    records = {r['path']: r for r in delivery['files']}
    assert len(records) == delivery['n_files']
    base = OUT.relative_to(ROOT).as_posix()
    with (STAGE / base / 'protocol/cohort.csv').open() as stream:
        cohort = list(csv.DictReader(stream))
    assert len(cohort) == 121 and len({r['sample'] for r in cohort}) == 121
    tasks = [(r['sample'], records[f"{base}/artifacts/{r['sample']}_all_core_labels_metrics.tar.gz"],
              int(r['n_cells'])) for r in cohort]
    results = []
    workers = min(4, int(os.environ.get('SLURM_CPUS_PER_TASK', '1')))
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for future in as_completed([pool.submit(inspect_archive, task) for task in tasks]):
            result = future.result()
            results.append(result)
            print('VERIFIED_ARCHIVE', len(results), result['sample'], flush=True)

    figures = 0
    for row in cohort:
        for budget in FEATURES:
            prefix = f"{base}/GBM/{row['sample']}/{budget}/figures"
            folder = STAGE / prefix
            figure_manifest = json.loads((folder / 'manifest.json').read_text())
            assert sha(folder / 'manifest.json') == (folder / 'FIGURES_COMPLETE').read_text().strip()
            assert figure_manifest['status'] == 'completed'
            assert figure_manifest['sample'] == row['sample'] and figure_manifest['budget'] == budget
            assert figure_manifest['clustering_routes'] == ROUTES
            assert figure_manifest['all_cells_plotted'] == int(row['n_cells'])
            assert {'all_cluster_routes.pdf', 'all_cluster_routes.png'} <= set(figure_manifest['files'])
            for filename, expected in figure_manifest['files'].items():
                assert expected == records[prefix + '/' + filename]['sha256']
                assert expected == sha(folder / filename)
            figures += 1

    import nbformat
    notebook_path = STAGE / 'notebooks/dgscrna_results.ipynb'
    assert sha(notebook_path) == receipt['notebook_sha256']
    assert sha(notebook_path) == sha(ROOT / 'notebooks/dgscrna_results.ipynb')
    notebook = nbformat.read(notebook_path, as_version=4)
    original_path = OUT / 'notebook_before_claim_validation.ipynb'
    update = json.loads((OUT / 'summary/notebook_update_manifest.json').read_text())
    assert sha(original_path) == update['previous_sha256']
    original = nbformat.read(original_path, as_version=4)
    assert notebook.cells[1:1 + len(original.cells)] == original.cells
    assert len(original.cells) == update['original_cells_preserved']
    assert len(notebook.cells) == len(original.cells) + update['added_cells']
    assert not any(o.get('output_type') == 'error' for c in notebook.cells for o in c.get('outputs', []))
    assert sum(r['terminal_predictions_verified'] for r in results) == 139392
    assert figures == 726
    dest = OUT / 'verification/core_archive_content_audit'
    write_json(dest / 'manifest.json', dict(
        status='completed', scope='Delivered GBM core contents only; full GBM and PTC gates remain separate',
        core_delivery_manifest_sha256=sha(manifest_path), receipt_sha256=sha(receipt_path),
        remote_check='Existing verified full-download receipt; no duplicate network transfer',
        samples=sorted(results, key=lambda r: r['sample']),
        n_archives=121, terminal_predictions_verified=139392,
        figure_sets_verified=figures, cluster_routes_shown=figures * len(ROUTES),
        original_notebook_cells_preserved=len(original.cells), notebook_cells=len(notebook.cells),
        notebook_sha256=sha(notebook_path),
        boundary='Results archives contain no external links. Absolute paths in provenance are historical metadata; raw fitting inputs and model weights are not part of this delivery audit.',
        computation='Hashes, metadata and notebook-cell equality only; predictions and metrics not recomputed',
        job=os.environ['SLURM_JOB_ID'], source_sha256=sha(__file__), completed_at=utc()))
    complete(dest)
    print('CORE_ARCHIVE_CONTENT_AUDIT_COMPLETE', flush=True)


if __name__ == '__main__':
    run()
