"""Check new GBM delivery contents without refitting or recomputing metrics."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
import json
import os
from pathlib import PurePosixPath
import tarfile

from common import ROOT, OUT, require_slurm, sha, write_json, complete, utc
from delivery import STAGE


def inspect_archive(task):
    sample, record = task
    archive_path = STAGE / record['path']
    assert sha(archive_path) == record['sha256']
    digests = {}
    links = 0
    with tarfile.open(archive_path, 'r|gz') as archive:
        for member in archive:
            path = PurePosixPath(member.name)
            assert not path.is_absolute() and '..' not in path.parts
            assert str(path) == member.name and member.name not in digests
            if member.islnk():
                assert member.linkname in digests
                digest = digests[member.linkname]
                links += 1
            else:
                assert member.isfile(), (sample, member.name, member.type)
                payload = archive.extractfile(member)
                assert payload is not None
                digest = hashlib.file_digest(payload, 'sha256').hexdigest()
            assert sha(OUT / member.name) == digest, (sample, member.name)
            digests[member.name] = digest

    # Require every saved prediction in the named families, independently of
    # the delivery writer's selected-filename list.
    expected = set()
    for method in ['scType', 'scCATCH', 'SCINA', 'SingleR', 'scDeepSort']:
        folder = OUT / 'comparators' / method / sample
        expected.update(str(p.relative_to(OUT)) for p in folder.rglob('predictions.csv.gz'))
        expected.update(str(p.relative_to(OUT)) for p in folder.rglob('cohort_predictions.csv.gz'))
        evaluation = folder / 'evaluation'
        manifest = json.loads((evaluation / 'manifest.json').read_text())
        assert manifest['status'] == 'completed' and manifest['sample'] == sample
        assert (evaluation / 'COMPLETE').read_text().strip() == sha(evaluation / 'manifest.json')
        expected.update(str(p.relative_to(OUT)) for p in evaluation.iterdir() if p.is_file())
        assert str(evaluation.relative_to(OUT)) + '/manifest.json' in digests
    for family in ['GBM_DL_controls', 'GBM_representation_controls']:
        folder = OUT / family / sample
        expected.update(str(p.relative_to(OUT)) for p in folder.rglob('predictions.csv.gz'))
    for budget in ['hvg2000', 'hvg5000', 'all']:
        folder = OUT / 'GBM' / sample / budget
        expected.update(str(p.relative_to(OUT)) for p in folder.glob('*/terminal_geometry_only_DL2000/*/predictions.csv.gz'))
        evaluation = folder / 'evaluation_geometry_only_DL2000'
        expected.update(str(p.relative_to(OUT)) for p in evaluation.iterdir() if p.is_file())
    assert expected <= set(digests), (sample, sorted(expected - set(digests))[:20])
    return dict(sample=sample, members=len(digests), internal_hardlinks=links,
                required_saved_outputs=len(expected), archive_sha256=record['sha256'])


def run():
    require_slurm()
    delivery_path = OUT / 'GBM_full_summary/delivery_manifest.json'
    delivery = json.loads(delivery_path.read_text())
    records = {r['path']: r for r in delivery['files']}
    assert len(records) == delivery['n_files']
    prefix = OUT.relative_to(ROOT).as_posix()
    cohort = list(csv.DictReader((OUT / 'protocol/cohort.csv').open()))
    assert len(cohort) == 121
    tasks = [(r['sample'], records[f"{prefix}/artifacts/{r['sample']}_controls_and_comparators.tar.gz"])
             for r in cohort]
    results = []
    with ProcessPoolExecutor(max_workers=min(4, int(os.environ['SLURM_CPUS_PER_TASK']))) as pool:
        for future in as_completed([pool.submit(inspect_archive, task) for task in tasks]):
            result = future.result()
            results.append(result)
            print('VERIFIED_FOLLOWUP_ARCHIVE', len(results), result['sample'], flush=True)
    archive_paths = {record['path'] for _, record in tasks}
    for name, record in records.items():
        path = STAGE / name
        assert path.stat().st_size == record['bytes'], name
        if name not in archive_paths:
            assert sha(path) == record['sha256'], name

    figure_counts = {}
    for index in ['controls_summary/representation_figure_index.csv',
                  'scalability_summary/all_resource_clustering_figures.csv']:
        rows = list(csv.DictReader((OUT / index).open()))
        for row in rows:
            for ext in ['png', 'pdf']:
                name = prefix + '/' + row[ext]
                assert name in records and (STAGE / name).is_file(), name
        figure_counts[index] = len(rows)
    assert sorted(figure_counts.values()) == [15, 180]

    import nbformat
    notebook_path = STAGE / 'notebooks/dgscrna_results.ipynb'
    update = json.loads((OUT / 'GBM_full_summary/notebook_manifest.json').read_text())
    assert sha(notebook_path) == update['current_sha256']
    backup = OUT / 'notebook_before_GBM_followups.ipynb'
    assert sha(backup) == update['previous_sha256']
    notebook = nbformat.read(notebook_path, as_version=4)
    previous = nbformat.read(backup, as_version=4)
    assert len(previous.cells) == update['old_cells_preserved'] == 647
    assert notebook.cells[1:1 + len(previous.cells)] == previous.cells
    assert len(notebook.cells) == len(previous.cells) + update['added_cells'] == 694
    assert not any(o.get('output_type') == 'error' for c in notebook.cells for o in c.get('outputs', []))
    dest = OUT / 'verification/GBM_followup_archive_content_audit'
    write_json(dest / 'manifest.json', dict(
        status='completed', delivery_manifest_sha256=sha(delivery_path),
        n_delivery_files=len(records), n_sample_archives=len(results),
        samples=sorted(results, key=lambda r: r['sample']), figure_index_counts=figure_counts,
        previous_notebook_cells_preserved=647, notebook_cells=694,
        notebook_sha256=sha(notebook_path),
        scope='New GBM follow-up package only; archived bytes match preserved source files, '
              'required saved predictions/evaluation files and indexed figures are present. '
              'No refitting, metric recomputation, repeat core audit or duplicate network transfer. '
              'OneDrive full-download receipt is verified separately.',
        job=os.environ['SLURM_JOB_ID'], source_sha256=sha(__file__), completed_at=utc()))
    complete(dest)
    print('GBM_FOLLOWUP_ARCHIVE_CONTENT_AUDIT_COMPLETE', flush=True)


if __name__ == '__main__':
    run()
