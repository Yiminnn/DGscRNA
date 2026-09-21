"""Package the exact small raw-count fixture and frozen markers outside the Git source tree."""
import argparse
import json
import os
from pathlib import Path
import tarfile
from run_gbm import HERE, sha, write_json


def main():
    assert os.environ.get('SLURM_JOB_ID'), 'Fixture packaging requires SLURM'
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--counts-dir', type=Path, required=True)
    p.add_argument('--markers', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    fixture = json.loads((HERE/'fixtures/TKU4163.json').read_text())
    files = {f'TKU4163/{name}': args.counts_dir/name for name in fixture['raw_counts_files']}
    for name, digest in fixture['raw_counts_files'].items():
        assert sha(args.counts_dir/name) == digest
    assert sha(args.markers) == fixture['marker_libraries_sha256']
    files['markers/libraries.json'] = args.markers
    files['fixture_manifest.json'] = HERE/'fixtures/TKU4163.json'
    assert not args.output.exists(), 'Do not overwrite prior bundles'
    args.output.mkdir(parents=True)
    bundle = args.output/'TKU4163_reference_counts_markers.tar.gz'
    with tarfile.open(bundle, 'w:gz') as archive:
        for name, source in files.items():
            archive.add(source, arcname=name, recursive=False)
    write_json(args.output/'fixture_bundle_manifest.json', dict(
        status='verified_source_fixture_packaged', bundle_sha256=sha(bundle), bytes=bundle.stat().st_size,
        members={name: sha(source) for name, source in files.items()},
        source_provenance='Public GSE274546 author-retained TKU4163 count fixture and unchanged frozen original marker libraries; no predictions, labels, model, environment or cached preprocessing.',
        source_sha256=sha(__file__), job=os.environ['SLURM_JOB_ID'], step=os.environ.get('SLURM_STEP_ID')))
    print('FIXTURE_BUNDLE_COMPLETE', bundle, flush=True)


if __name__ == '__main__':
    main()
