"""Portable semantic fingerprints of the fixed fixture; no biological truth is loaded."""
from pathlib import Path
import hashlib
import os


def fingerprints(prep):
    assert os.environ.get('SLURM_JOB_ID')
    import numpy as np
    import pandas as pd
    prep = Path(prep)
    route = prep/'UMAP2_HDBSCAN_R'
    terminal = route/'terminal/L00_mean'
    result = {}
    for name in ['selected_features.txt', 'geometry_features.txt', 'scoring_features.txt', 'DL_features.txt', 'DL.float32.bin', 'cells.csv']:
        with (prep/name).open('rb') as stream:
            result[name] = hashlib.file_digest(stream, 'sha256').hexdigest()
    for path in [prep/'PCA30.csv', prep/'UMAP2.csv', *[route/name for name in ['clusters.csv', 'density_diagnostics.csv', 'initial_calls.csv.gz', 'cluster_calls.csv.gz', 'marker_retention.csv.gz']]]:
        # Re-serialization removes gzip headers and preserves fixed parsed values/order.
        frame = pd.read_csv(path)
        result[str(path.relative_to(prep))+':parsed_csv'] = hashlib.sha256(frame.to_csv(index=False).encode()).hexdigest()
    with np.load(terminal/'terminal.npz', allow_pickle=False) as arrays:
        for name in arrays.files:
            a = arrays[name]
            result['terminal:'+name] = dict(dtype=str(a.dtype), shape=list(a.shape), sha256=hashlib.sha256(a.tobytes(order='C')).hexdigest())
    return result
