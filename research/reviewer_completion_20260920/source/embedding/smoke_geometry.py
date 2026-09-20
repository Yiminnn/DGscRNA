"""Exercise all seven reducer/clusterer interfaces and resources without reading truth.

Uses a separate output tree so running production pilots cannot race its caches.
The same frozen numerical functions and protocol are used; no DEG/DL or ranking.
"""
from pathlib import Path
import json
import os
import sys
import time

source = Path(sys.argv[1]).resolve()
sys.path.insert(0, str(source))
import run as core
core.require_slurm()
tasks = json.loads(Path(sys.argv[2]).read_text())
task = tasks[int(os.environ['SLURM_ARRAY_TASK_ID'])]
sample, budget = task['sample'], task['budget']
base = core.CAMPAIGN / 'embedding_geometry_smoke' / sample / budget
shared = core.OUTPUT / sample / budget / 'geometry'
geometry = shared if core.checked(shared) else base / 'geometry'
records = []
start = time.monotonic()
for space in core.SPACES:
    dest = base / space;dest.mkdir(parents=True, exist_ok=True)
    cfg = dict(sample=sample, budget=budget, space=space,
        prep=str(core.REFERENCE / 'GBM' / sample / budget), dest=str(dest),
        geometry=str(geometry), embedding=str(dest / 'embedding.csv'),
        hdbscan_dest=str(dest / 'HDBSCAN_R'),
        frozen_protocol_sha256=core.sha(core.CAMPAIGN / 'protocol/embedding.json'),
        source_sha256=core.sha(source / 'run.py'), reference_labels_used_for_fit=False,
        purpose='isolated interface/resource gate; no scientific ranking')
    core.write_json(dest / 'config.json', cfg)
    x, cells, gm = core.get_geometry(cfg)
    z = core.make_embedding(cfg, x, cells, gm)
    conditions = core.make_partitions(cfg, z, cells)
    records.append(dict(space=space, n_conditions=len(conditions),
        representation_sha256=core.sha(dest / 'representation.json'),
        partition_manifests={c['route']:core.sha(Path(c['dest']) / 'partition_manifest.json') for c in conditions}))
    core.log('geometry_smoke_space_complete', sample=sample, budget=budget, space=space)
core.write_json(base / 'manifest.json', dict(status='passed', sample=sample, budget=budget,
    n_cells=len(cells), n_features=x.shape[1], n_partitions=sum(r['n_conditions'] for r in records),
    spaces=records, numeric_source=str(source), numeric_source_sha256=core.sha(source / 'run.py'),
    protocol_sha256=cfg['frozen_protocol_sha256'], input_binary_sha256=gm['binary_sha256'],
    elapsed_seconds=time.monotonic()-start, job=os.environ['SLURM_JOB_ID'],
    reference_labels_used_for_fit=False, completed_at=core.utc()))
core.complete(base)
