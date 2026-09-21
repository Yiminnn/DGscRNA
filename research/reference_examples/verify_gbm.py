"""Numerical parity against an explicit frozen reference, after a fresh fit."""
from pathlib import Path
import argparse
import json
import os
import subprocess
from run_gbm import sha, write_json, source_audit, HERE


def main():
    assert os.environ.get('SLURM_JOB_ID'), 'Scientific validation requires SLURM'
    import numpy as np
    import pandas as pd
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--reference', type=Path, required=True, help='Frozen TKU4163/hvg2000 preparation directory')
    parser.add_argument('--proof-name', default='parity.json', help='New proof filename; prior proofs are preserved')
    args = parser.parse_args()
    assert Path(args.proof_name).name == args.proof_name
    assert not (args.output/args.proof_name).exists(), 'Use a new proof name'
    assert (args.output/'COMPLETE').read_text().strip() == sha(args.output/'run_manifest.json')
    source_audit()
    fixture = json.loads((HERE/'fixtures/TKU4163.json').read_text())
    actual = args.output/'GBM'/fixture['sample']/fixture['budget']
    expected = args.reference
    exact = []
    for name in ['selected_features.txt', 'geometry_features.txt', 'scoring_features.txt', 'DL_features.txt', 'DL.float32.bin', 'cells.csv']:
        assert sha(actual/name) == sha(expected/name), name
        exact.append(name)
    for name in ['PCA30.csv', 'UMAP2.csv']:
        pd.testing.assert_frame_equal(pd.read_csv(actual/name), pd.read_csv(expected/name), check_exact=True)
        exact.append(name)
    route = fixture['route']
    for name in ['clusters.csv', 'density_diagnostics.csv', 'initial_calls.csv.gz', 'cluster_calls.csv.gz', 'marker_retention.csv.gz']:
        pd.testing.assert_frame_equal(pd.read_csv(actual/route/name), pd.read_csv(expected/route/name), check_exact=True)
        exact.append(route+'/'+name)
    terminal = Path(route)/'terminal'/fixture['arm']
    with np.load(actual/terminal/'terminal.npz', allow_pickle=False) as a, np.load(expected/terminal/'terminal.npz', allow_pickle=False) as b:
        assert set(a.files) == set(b.files)
        for name in a.files:
            if np.issubdtype(a[name].dtype, np.floating):
                np.testing.assert_allclose(a[name], b[name], rtol=0, atol=0, equal_nan=True, err_msg=name)
            else:
                assert np.array_equal(a[name], b[name]), name
        exact.extend('terminal.npz:'+name for name in a.files)
    assert json.loads((actual/terminal/'training_history.json').read_text()) == json.loads((expected/terminal/'training_history.json').read_text())
    tm = json.loads((actual/terminal/'terminal_manifest.json').read_text())
    assert tm['training_executed'] and not tm['identical_result_reused']
    run_manifest = json.loads((args.output/'run_manifest.json').read_text())
    subprocess.run([run_manifest['rscript'], '--vanilla', str(HERE/'verify_r_artifacts.R'),
                    str(actual/route), str(expected/route)], check=True)
    write_json(args.output/args.proof_name, dict(status='passed_exact', scope='Fixed TKU4163 HVG2000 native R UMAP2 HDBSCAN CM2_glioma_other mean terminal090/070',
               exact_artifacts=exact, history_exact=True, R_DEG_and_density_objects_exact=True, actual_fresh_DL=True,
               reference=str(expected.resolve()), reference_prepare_manifest_sha256=sha(expected/'prepare_manifest.json'),
               reference_terminal_manifest_sha256=sha(expected/terminal/'terminal_manifest.json'),
               run_manifest_sha256=sha(args.output/'run_manifest.json'), source_sha256=sha(__file__),
               job=os.environ['SLURM_JOB_ID']))
    print('EXACT_PARITY_PASSED', args.output, flush=True)


if __name__ == '__main__':
    main()
