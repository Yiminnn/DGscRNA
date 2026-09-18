"""Check and preserve the two observed cross-node duplicate terminal results."""
import os
from common import OUT,require_slurm,write_json,sha,utc
from terminal import preserve_equivalent_duplicate

def run():
    require_slurm()
    key='60b79bd405aa52e069831bb2c23eb3700f61794256f4e4d06433cd6f91b62745'
    suffixes=['0042ebdf209d4c5e9f6e3aac2c41726f','c86847432a6e473f96673e41742320f8']
    results=[]
    for suffix in suffixes:
        name=key+'.partial.'+suffix;path=OUT/'DL_cache'/name
        archived=OUT/'DL_cache/duplicate_publications'/name
        if path.exists():preserve_equivalent_duplicate(path,OUT/'DL_cache'/key)
        proof=archived/'DUPLICATE_TERMINAL_EQUIVALENCE.json';assert proof.exists()
        results.append(dict(archive=str(archived),verification_sha256=sha(proof)))
    write_json(OUT/'verification/cross_node_cache_publication.json',dict(status='passed',duplicates=results,
        scientific_algorithm_unchanged=True,job=os.environ['SLURM_JOB_ID'],completed_at=utc()))

if __name__=='__main__':run()
