from pathlib import Path
import hashlib
import json
import os
import sys

CODE = Path(__file__).resolve().parent
sys.path.insert(0, str(CODE/'legacy'))
from common import ROOT, OUT as REFERENCE, RSCRIPT, require_slurm, sha, utc, write_json, checked, complete, L1
OUT = ROOT/'results/hvg_ptc_20260916_v1/reviewer_completion_20260920/no_clustering'


def verify_source():
    source = json.loads((CODE/'SOURCE_MANIFEST.json').read_text())
    for relative, expected in source['files'].items():
        assert sha(CODE/relative) == expected, f'Frozen source changed: {relative}'
    protocol = json.loads((OUT/'protocol.json').read_text())
    assert protocol['source_bundle_sha256'] == sha(CODE/'SOURCE_MANIFEST.json')
    return protocol


def configure(sample, budget):
    protocol = verify_source()
    assert sample in protocol['samples'] and budget in protocol['budgets']
    dest = OUT/'GBM'/sample/budget
    dest.mkdir(parents=True, exist_ok=True)
    cfg = dict(sample=sample, budget=budget, dest=str(dest), route_dir=str(dest/'cellwise_seed'),
               prep=str(REFERENCE/'GBM'/sample/budget), marker_file=str(REFERENCE/'markers/libraries.json'),
               marker_sha256=protocol['marker_sha256'], protocol_sha256=sha(OUT/'protocol.json'),
               source_bundle_sha256=sha(CODE/'SOURCE_MANIFEST.json'),
               input=protocol['inputs'][sample+'/'+budget], arms=protocol['arms'])
    cfg['input_signature'] = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()
    if (dest/'config.json').exists():
        assert json.loads((dest/'config.json').read_text()) == cfg, 'Changed A2 input; preserve old outputs before rebuilding'
    else:
        write_json(dest/'config.json', cfg)
    return cfg


def verify_terminal(route, aid):
    route = Path(route);td = route/'terminal'/aid
    assert checked(route, 'score_manifest.json', 'SCORE_COMPLETE')
    score = json.loads((route/'score_manifest.json').read_text())
    assert score['initial_sha256'] == sha(route/'initial_calls.csv.gz')
    assert checked(td, 'terminal_manifest.json', 'TERMINAL_COMPLETE')
    tm = json.loads((td/'terminal_manifest.json').read_text())
    assert tm['score_manifest_sha256'] == sha(route/'score_manifest.json')
    assert tm['terminal_sha256'] == sha(td/'terminal.npz')
    assert tm['predictions_sha256'] == sha(td/'predictions.csv.gz')
    assert tm['DL_sha256'] == score['DL_binary_sha256']
    return tm
