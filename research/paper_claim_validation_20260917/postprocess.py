"""SLURM-only full-condition evaluation, validation and every-route plot."""
import sys
from common import require_slurm
import evaluate
import plot_unit

if __name__=='__main__':
    require_slurm()
    evaluate.run(sys.argv[1])
    plot_unit.run(sys.argv[1])
