"""Scheduling profiles only: no expression matrices or annotation labels are read."""
import csv
import gzip
import json
import os
from pathlib import Path
import re
from common import ROOT, OUT, require_slurm, checked, complete, sha, write_json, utc

DIMENSIONS=('n_cells','n_shared_genes','n_clusters')

def covering_pilot(profile,pilots):
    for name,pilot in sorted(pilots.items()):
        if all(profile[key]<=pilot[key] for key in DIMENSIONS):return name
    return None

def run():
    require_slurm()
    source=Path(__file__).resolve().parent
    roster_path=OUT/'markers/scCATCH/manifest.json'
    roster=json.loads(roster_path.read_text());genes=set()
    for entry in roster['libraries']:
        p=Path(entry['path']);assert sha(p)==entry['sha256']
        with gzip.open(p,'rt') as f:genes.update(r['gene'] for r in csv.DictReader(f))
    profiles={}
    for row in csv.DictReader((OUT/'protocol/sample_order.csv').open()):
        sample=row['sample'];base=OUT/'GBM'/sample/'hvg2000';route=base/'UMAP2_HDBSCAN_R'
        assert checked(base,'prepare_manifest.json','PREPARED')
        assert checked(route,'score_manifest.json','SCORE_COMPLETE')
        score=json.loads((route/'score_manifest.json').read_text())
        features=base/'scoring_features.txt'
        shared=genes.intersection(features.read_text().splitlines())
        assert score['n_cells']==int(row['n_cells']) and score['n_clusters']>0
        profiles[sample]=dict(n_cells=score['n_cells'],n_clusters=score['n_clusters'],
            n_shared_genes=len(shared),scoring_features_sha256=sha(features),
            prepare_manifest_sha256=sha(base/'prepare_manifest.json'),
            score_manifest_sha256=sha(route/'score_manifest.json'))
    # Verify the inferred dimensions against the completed native method's own log.
    middle=OUT/'comparators/scCATCH/NL022/hvg2000/UMAP2_HDBSCAN_R'
    assert checked(middle,'cohort_manifest.json','COHORT_COMPLETE')
    manifest=json.loads((middle/'cohort_manifest.json').read_text())
    assert manifest['source_sha256']==sha(source/'sccatch_R.R')
    logs=list((ROOT/'logs').glob(f'*_{manifest["job"]}_*.out'));assert len(logs)==1,logs
    observed=re.search(r'SCCATCH_INPUT NL022 hvg2000 UMAP2_HDBSCAN_R (\d+) (\d+)',logs[0].read_text())
    assert observed
    assert int(observed[1])==profiles['NL022']['n_shared_genes']
    assert int(observed[2])==profiles['NL022']['n_cells']
    # Bound each dimension separately; a small cell count never excuses more clusters.
    anchor={k:10 for k in DIMENSIONS};cases=[]
    for changed in [None,*DIMENSIONS]:
        candidate=dict(anchor)
        if changed:candidate[changed]=11
        assert covering_pilot(candidate,{'anchor':anchor})==('anchor' if changed is None else None)
        cases.append('equal_bounds' if changed is None else 'reject_larger_'+changed)
    assert covering_pilot({k:9 for k in DIMENSIONS},{'anchor':anchor})=='anchor'
    assert covering_pilot(anchor,{}) is None
    cases.extend(['smaller_all_dimensions','no_completed_pilot'])
    eligible=[s for s,p in profiles.items() if covering_pilot(p,{'NL022':profiles['NL022']})]
    target=OUT/'verification/scCATCH_resource_profiles';target.mkdir(exist_ok=True)
    write_json(target/'manifest.json',dict(status='verified_resource_profiles',profiles=profiles,
        dimensions=list(DIMENSIONS),eligible_under_completed_NL022=eligible,
        policy='Before all three pilots pass, require a completed matched-partition pilot no smaller in cells, shared marker genes, or clusters. This is conservative scheduling, not a runtime bound or selection by annotation performance.',
        marker_manifest_sha256=sha(roster_path),r_source_sha256=sha(source/'sccatch_R.R'),
        native_log_sha256=sha(logs[0]),native_input_line=observed[0],boundary_checks=cases,
        source_sha256=sha(__file__),job=os.environ['SLURM_JOB_ID'],completed_at=utc()))
    complete(target)
    print('SCCATCH_PROFILES',len(profiles),'ELIGIBLE_UNDER_NL022',len(eligible),flush=True)

if __name__=='__main__':run()
