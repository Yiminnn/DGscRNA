"""All core units must pass before notebook, source and authorized directory delivery."""
import os
import subprocess
import sys
from pathlib import Path
from common import OUT, ROOT, require_slurm, checked, sha, write_json, utc

def run():
    require_slurm()
    import aggregate
    import report
    import delivery
    aggregate.run()
    report.run()
    delivery.run()
    write_json(OUT/'summary/GBM_CORE_DELIVERED.json',dict(status='GBM_core_delivered',
        gate='726 units evaluated/plotted; original notebook extended; complete OneDrive payload downloaded and checked',
        followups='Geometry-only, fairness, parameter, PTC and resource gates remain separate; not whole-campaign completion.',
        notebook_sha256=sha(ROOT/'notebooks/dgscrna_results.ipynb'),
        job=os.environ['SLURM_JOB_ID'],completed_at=utc()))

if __name__=='__main__':run()
