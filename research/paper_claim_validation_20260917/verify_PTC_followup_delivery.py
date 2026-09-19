"""Validate the first complete PTC delivery package against preserved results."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import tarfile

from common import ROOT, OUT, sha, checked, write_json, complete, utc
from delivery import STAGE
from ptc_followup_common import PTC, CONTROL_ROUTES, require_ptc


def inspect_archive(task):
    record, base, required = task
    base = Path(base)
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
                assert member.isfile(), (record['path'], member.name, member.type)
                payload = archive.extractfile(member)
                assert payload is not None
                digest = hashlib.file_digest(payload, 'sha256').hexdigest()
            assert sha(base / member.name) == digest, (record['path'], member.name)
            digests[member.name] = digest
    assert set(required) <= set(digests), (record['path'], sorted(set(required) - set(digests))[:20])
    return dict(archive=record['path'], archive_sha256=record['sha256'], members=len(digests),
                internal_hardlinks=links, required_saved_outputs=len(required))


def run():
    require_ptc()
    import nbformat
    import pandas as pd
    from ptc_control_job import context_arms
    summary = OUT / 'PTC_summary'
    assert checked(summary)
    receipt_path = summary / 'DELIVERY_RECEIPT.json'
    receipt = json.loads(receipt_path.read_text())
    assert receipt['status'] == 'delivered_and_verified'
    delivery_path = summary / 'delivery_manifest.json'
    assert receipt['manifest_sha256'] == sha(delivery_path)
    delivery = json.loads(delivery_path.read_text())
    records = {r['path']: r for r in delivery['files']}
    assert len(records) == delivery['n_files']
    prefix = OUT.relative_to(ROOT).as_posix()
    mlp = json.loads((PTC / 'selection/MLP_tasks.json').read_text())
    preparations = json.loads((PTC / 'selection/preparation_tasks.json').read_text())
    assert len(mlp) == 50 and len(preparations) == 22
    terminals = {Path(c['dest']) / 'terminal' / a for c in mlp for a in c['arm_ids']}
    for cfg in preparations:
        assert checked(Path(cfg['dest']) / 'evaluation')
        for route in CONTROL_ROUTES:
            terminals.update(Path(cfg['dest']) / route / 'terminal' / a for a in context_arms(cfg, route))
    ledger = pd.read_csv(summary / 'terminal_execution_ledger.csv.gz', keep_default_na=False)
    assert ledger.terminal_directory.is_unique
    assert set(ledger.terminal_directory) == {str(p) for p in terminals}
    required_by_unit = {}
    for terminal in terminals:
        assert checked(terminal, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
        unit = PTC.joinpath(*terminal.relative_to(PTC).parts[:2])
        required = required_by_unit.setdefault(unit, set())
        for name in ['predictions.csv.gz', 'terminal_manifest.json', 'training_manifest.json', 'TERMINAL_COMPLETE']:
            path = terminal / name
            assert path.is_file()
            required.add(str(path.relative_to(PTC)))
        if (terminal / 'training_history.json').exists():
            required.add(str((terminal / 'training_history.json').relative_to(PTC)))
    tasks = []
    family_counts = {}
    mandatory_saved_names = {'metrics.csv.gz', 'terminal_statuses.csv.gz', 'initial_calls.csv.gz',
        'clusters.csv', 'cells.csv', 'score_manifest.json', 'prepare_manifest.json', 'manifest.json',
        'COMPLETE', 'PREPARED', 'SCORE_COMPLETE', 'default_parity.json',
        'retention_invariants.json', 'retention_DEG_audit.json'}
    for family, expected_count in [('representation', 16), ('marker_retention', 6),
                                    ('MLP_seeds', 50), ('default_parity', 2), ('existing_grid_evaluation', 30)]:
        units = sorted(p for p in (PTC / family).iterdir() if p.is_dir())
        assert len(units) == expected_count, (family, len(units))
        family_counts[family] = len(units)
        for unit in units:
            required = required_by_unit.get(unit, set()).copy()
            required.update(str(p.relative_to(PTC)) for p in unit.rglob('*')
                            if p.is_file() and p.name in mandatory_saved_names)
            assert required
            record = records[f'{prefix}/PTC_followups/artifacts/{unit.name}.tar.gz']
            tasks.append((record, str(PTC), sorted(required)))
    dispatch = json.loads((PTC / 'dispatch_state.json').read_text())
    source_root = OUT / 'source_snapshots'
    source_directories = {Path(dispatch['science_source']), Path(dispatch['finalizer_source'])}
    required_source = [str(p.relative_to(source_root)) for s in source_directories for p in s.iterdir()
                       if p.suffix in ['.py', '.R', '.sbatch', '.json']]
    tasks.append((records[f'{prefix}/PTC_followups/artifacts/PTC_immutable_sources.tar.gz'],
                  str(source_root), required_source))
    archive_paths = {r['path'] for r, _, _ in tasks}
    assert archive_paths == {n for n in records if n.endswith('.tar.gz')}
    results = []
    with ProcessPoolExecutor(max_workers=min(4, int(os.environ['SLURM_CPUS_PER_TASK']))) as pool:
        for future in as_completed([pool.submit(inspect_archive, task) for task in tasks]):
            results.append(future.result())
            print('VERIFIED_PTC_ARCHIVE', len(results), results[-1]['archive'], flush=True)
    for name, record in records.items():
        path = STAGE / name
        assert path.stat().st_size == record['bytes'], name
        if name not in archive_paths:
            assert sha(path) == record['sha256'], name
    figures = list(csv.DictReader((summary / 'all_new_clustering_figures.csv').open()))
    assert len(figures) == 44
    for row in figures:
        for ext in ['png', 'pdf']:
            assert prefix + '/' + row[ext] in records
    for ext in ['png', 'pdf']:
        assert f'{prefix}/summary/workflow_decision_tree.{ext}' in records
    notebook_path = STAGE / 'notebooks/dgscrna_results.ipynb'
    update = json.loads((summary / 'notebook_manifest.json').read_text())
    assert sha(notebook_path) == update['after_sha256'] == receipt['notebook_sha256']
    backup = PTC / 'notebook_before_PTC_followups.ipynb'
    backup_relative = str(backup.relative_to(ROOT))
    assert backup_relative in records and sha(STAGE / backup_relative) == sha(backup) == update['before_sha256']
    notebook = nbformat.read(notebook_path, as_version=4)
    previous = nbformat.read(backup, as_version=4)
    assert len(previous.cells) == update['old_cells_preserved'] == 694
    assert notebook.cells[1:1 + len(previous.cells)] == previous.cells
    assert len(notebook.cells) == len(previous.cells) + update['added_cells']
    assert not any(o.get('output_type') == 'error' for c in notebook.cells for o in c.get('outputs', []))
    dest = PTC / 'verification/PTC_followup_archive_content_audit'
    write_json(dest / 'manifest.json', dict(status='completed',
        delivery_manifest_sha256=sha(delivery_path), delivery_receipt_sha256=sha(receipt_path),
        n_delivery_files=len(records), n_archives=len(results), archive_family_counts=family_counts,
        archives=sorted(results, key=lambda r: r['archive']), n_terminal_conditions=len(terminals),
        terminal_ledger_exactly_covers_requested_conditions=True, n_new_clustering_figure_pairs=len(figures),
        previous_notebook_cells_preserved=len(previous.cells), notebook_cells=len(notebook.cells),
        notebook_sha256=sha(notebook_path), backup_notebook_staged_and_hashed=True,
        scope='First complete PTC follow-up package content check. Archived bytes match preserved sources; '
              'requested terminal predictions, evaluations, figures and both notebook versions are present. '
              'No refitting, metric recomputation, repeat GBM audit or duplicate full-package transfer.',
        job=os.environ['SLURM_JOB_ID'], source_sha256=sha(__file__), completed_at=utc()))
    complete(dest)
    print('PTC_FOLLOWUP_ARCHIVE_CONTENT_AUDIT_COMPLETE', flush=True)


if __name__ == '__main__':
    run()
