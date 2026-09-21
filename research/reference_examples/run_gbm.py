"""Fresh single-sample native R reference run; requires SLURM and explicit inputs."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def source_audit():
    manifest = json.loads((HERE/'SOURCE_DERIVATION.json').read_text())
    for name, item in manifest['files'].items():
        original = HERE/'reference_sources'/name
        backend = HERE/'backend'/name
        assert sha(original) == item['source_sha256'], name
        assert sha(backend) == item['backend_sha256'], name
        text = original.read_text()
        for replacement in item['path_only_replacements']:
            assert text.count(replacement['before']) == 1
            text = text.replace(replacement['before'], replacement['after'])
        assert text == backend.read_text(), name
    return sha(HERE/'SOURCE_DERIVATION.json')


def prepare_counts(args, fixture):
    import numpy as np
    import pandas as pd
    from scipy import sparse
    from scipy.io import mmread
    raw = args.counts_dir.resolve()
    for name, expected in fixture['raw_counts_files'].items():
        assert sha(raw/name) == expected, (name, 'fixture input changed')
    cells = [v.split('\t')[0].strip() for v in (raw/'barcodes.tsv').read_text().splitlines() if v.strip()]
    genes = [v.split('\t')[-1].strip() for v in (raw/'genes.tsv').read_text().splitlines() if v.strip()]
    assert len(set(cells)) == len(cells) and len(set(genes)) == len(genes)
    counts = sparse.csr_matrix(mmread(raw/'matrix.mtx'))
    assert counts.shape == (len(genes), len(cells)), 'Fixture is genes x cells'
    counts = counts.T.tocsr()
    assert np.isfinite(counts.data).all() and (counts.data >= 0).all()
    assert np.array_equal(counts.data, np.round(counts.data))
    counts.sum_duplicates(); counts.eliminate_zeros()
    keep = np.asarray((counts > 0).sum(axis=0)).ravel() >= 3
    # Reproduce the historical float32 filtered-count export, then R float64 input.
    x = counts[:, keep].astype(np.float32).astype(np.float64).tocsr()
    x.sum_duplicates(); x.eliminate_zeros(); x.sort_indices()
    assert (np.asarray(x.sum(axis=1)).ravel() > 0).all()
    assert x.shape == (fixture['n_cells'], fixture['n_genes']) and x.nnz == fixture['nnz']
    dest = args.output/'inputs'/fixture['sample']
    dest.mkdir(parents=True)
    x.data.astype('<f8').tofile(dest/'x.bin')
    x.indices.astype('<i4').tofile(dest/'i.bin')
    x.indptr.astype('<i4').tofile(dest/'p.bin')
    pd.DataFrame({'cell_id': cells, 'batch': fixture['sample']}).to_csv(dest/'cells_fit.csv', index=False)
    pd.DataFrame({'gene': np.asarray(genes)[keep]}).to_csv(dest/'genes.csv', index=False)
    hashes = {name: sha(dest/name) for name in fixture['fitting_files']}
    assert hashes == fixture['fitting_files'], 'Raw-count adapter differs from frozen reference'
    manifest = dict(status='completed', sample=fixture['sample'], n_cells=x.shape[0], n_genes=x.shape[1],
                    nnz=x.nnz, input_semantics='author-filtered nonnegative integer RNA counts; genes detected in >=3 cells',
                    fitting_files=hashes, truth_labels_excluded_from_fit=True,
                    raw_input_files={name: sha(raw/name) for name in fixture['raw_counts_files']},
                    source_sha256=sha(__file__), job=os.environ['SLURM_JOB_ID'])
    write_json(dest/'input_manifest.json', manifest)
    (dest/'INPUT_COMPLETE').write_text(sha(dest/'input_manifest.json')+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--counts-dir', type=Path, required=True)
    parser.add_argument('--markers', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='New nonexistent directory; old predictions are never reused')
    parser.add_argument('--rscript', type=Path, required=True)
    parser.add_argument('--reference-r-lib', type=Path, help='Explicit optional patched dbscan library; omit for clean prefix')
    args = parser.parse_args()
    assert os.environ.get('SLURM_JOB_ID'), 'Scientific computation requires SLURM'
    assert not sys.flags.optimize
    args.output = args.output.resolve()
    assert not args.output.exists(), 'Use a fresh output directory'
    args.output.mkdir(parents=True)
    derivation_hash = source_audit()
    fixture = json.loads((HERE/'fixtures/TKU4163.json').read_text())
    assert sha(args.markers) == fixture['marker_libraries_sha256']
    prepare_counts(args, fixture)
    (args.output/'markers').mkdir()
    shutil.copyfile(args.markers, args.output/'markers/libraries.json')
    environment = os.environ.copy()
    environment.update(DGSCRNA_EXAMPLE_OUT=str(args.output), DGSCRNA_ONLY_ROUTE=fixture['route'],
                       DGSCRNA_REFERENCE_R_LIB=str(args.reference_r_lib.resolve()) if args.reference_r_lib else '',
                       DGSCRNA_PYTHON=sys.executable, DGSCRNA_RSCRIPT=str(args.rscript.resolve()),
                       DGSCRNA_DEG_WORKERS='2', PYTHONPATH=str(HERE/'backend'),
                       PYTHONNOUSERSITE='1', OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1',
                       MKL_NUM_THREADS='1', NUMBA_NUM_THREADS='1')
    environment.pop('DGSCRNA_DL_CACHE_ROOT', None)
    prep = args.output/'GBM'/fixture['sample']/fixture['budget']
    commands = [
        [str(args.rscript), '--vanilla', str(HERE/'backend/prepare_R.R'), fixture['sample'], fixture['budget']],
        [str(args.rscript), '--vanilla', str(HERE/'backend/score_R.R'), fixture['sample'], str(prep)],
        [sys.executable, '-s', str(HERE/'backend/terminal.py'), str(prep/fixture['route']), fixture['arm']],
    ]
    write_json(args.output/'RUNNING.json', dict(job=os.environ['SLURM_JOB_ID'], commands=commands,
               source_derivation_sha256=derivation_hash, fixture_sha256=sha(HERE/'fixtures/TKU4163.json')))
    for command in commands:
        print('RUN', command, flush=True)
        subprocess.run(command, env=environment, check=True)
    terminal = prep/fixture['route']/'terminal'/fixture['arm']
    result = json.loads((terminal/'terminal_manifest.json').read_text())
    assert result['terminal_valid'] and result['training_executed'] and not result['identical_result_reused']
    # Golden fingerprints are opened only after terminal predictions are frozen.
    from output_fingerprints import fingerprints
    expected_path = HERE/'fixtures/TKU4163_expected_outputs.json'
    expected = json.loads(expected_path.read_text())
    assert fingerprints(prep) == expected['fingerprints'], 'Fresh fixture output differs from frozen historical reference'
    write_json(args.output/'run_manifest.json', dict(status='completed_fresh_fit', sample=fixture['sample'],
               counts_adapter_exact=True, no_prediction_cache_reuse=True, terminal=str(terminal),
               dl_status=result['dl_status'], training_executed=result['training_executed'],
               python=sys.executable, rscript=str(args.rscript.resolve()),
               explicit_R_library=environment['DGSCRNA_REFERENCE_R_LIB'],
               source_derivation_sha256=derivation_hash, fixture_sha256=sha(HERE/'fixtures/TKU4163.json'),
               expected_output_sha256=sha(expected_path), expected_output_fingerprints_exact=True,
               runner_sha256=sha(__file__),
               terminal_manifest_sha256=sha(terminal/'terminal_manifest.json'), job=os.environ['SLURM_JOB_ID']))
    (args.output/'COMPLETE').write_text(sha(args.output/'run_manifest.json')+'\n')
    print('FRESH_REFERENCE_COMPLETE', terminal, flush=True)


if __name__ == '__main__':
    main()
