#!/usr/bin/env python3
"""Parallel early check of all seven new-input geometries; no DEG or DL fitting."""
from pathlib import Path
import argparse
import json
import os
import run as core


def main():
    assert os.environ.get('SLURM_JOB_ID')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--package-run', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--expected', type=Path, required=True)
    parser.add_argument('--rscript', required=True)
    args = parser.parse_args()
    assert not args.out.exists(), 'Fresh geometry smoke output required'
    core.OUTPUT, core.RSCRIPT = args.out.resolve(), str(Path(args.rscript).resolve())
    os.environ['DGSCRNA_REFERENCE_OUT'] = str(core.OUTPUT)
    os.environ['DGSCRNA_REQUIRE_SLURM'] = '1'
    manifest = core.verify_source()
    import adaptive_ica
    core.adaptive_ica = adaptive_ica
    core.load_frozen_functions(manifest)
    inputs = core.input_contract(args.package_run.resolve())
    runtime = core.installed_contract(inputs['installed_reference_sources'])
    # Match the fitting entry point's native BLAS/OpenMP initialization exactly.
    from dgscrna.reference.backend import terminal
    terminal.threads()
    cases = []
    # Finish new geometries/partitions before any historical coordinates are read.
    for space in core.SPACES:
        cfg = core.configuration(inputs, runtime, space)
        x, cells, geometry = core.get_geometry(cfg)
        z = core.make_embedding(cfg, x, cells, geometry)
        conditions = core.make_partitions(cfg, z, cells)
        assert len(conditions) == 13
        cases.append((space, cfg, conditions))
    import pandas as pd
    comparisons = []
    for space, cfg, conditions in cases:
        new = Path(cfg['dest'])
        old = args.expected.resolve() / space
        if space != 'noDR':
            pd.testing.assert_frame_equal(pd.read_csv(new / 'embedding.csv', dtype={'cell_id': str}, keep_default_na=False),
                pd.read_csv(old / 'embedding.csv', dtype={'cell_id': str}, keep_default_na=False), check_exact=True)
        for condition in conditions:
            partition = Path(condition['dest']).name
            pd.testing.assert_frame_equal(pd.read_csv(new / partition / 'clusters.csv', dtype=str, keep_default_na=False),
                pd.read_csv(old / partition / 'clusters.csv', dtype=str, keep_default_na=False), check_exact=True)
            comparisons.append(dict(space=space, partition=partition, exact=True))
    result = dict(status='passed_exact', sample=inputs['sample'], budget=inputs['budget'],
        n_spaces=7, n_partitions=len(comparisons), comparisons=comparisons,
        no_DEG_or_DL_fitting=True, no_old_geometry_input=True, reference_labels_read=False,
        adapter_sha256=core.sha(core.CODE / 'run.py'), source_manifest_sha256=core.sha(core.CODE / 'SOURCE_PROTOCOL_MANIFEST.json'),
        job=os.environ['SLURM_JOB_ID'])
    core.write_json(args.out / 'GEOMETRY_PROOF.json', result)
    print(json.dumps({k: result[k] for k in ['status', 'n_spaces', 'n_partitions']}, indent=2))


if __name__ == '__main__':
    main()
