"""Command line interface for the shared reference implementation."""
import argparse
import json
import sys

from .. import __version__


def main(argv=None):
    parser = argparse.ArgumentParser(prog="dgscrna", description="Original R workflow + Python DL/refinement")
    parser.add_argument("--version", action="version", version=__version__)
    sub = parser.add_subparsers(dest="command", required=True)
    doctor = sub.add_parser("doctor", help="Check R and Python dependencies without fitting")
    doctor.add_argument("--rscript", default="Rscript")
    doctor.add_argument("--reference-r-lib")
    run = sub.add_parser("run", help="Annotate already-QC single-sample counts")
    run.add_argument("--counts", required=True)
    run.add_argument("--markers", required=True)
    run.add_argument("--out", required=True)
    run.add_argument("--sample")
    run.add_argument("--preset", default="gbm-reference", choices=["gbm-reference"])
    run.add_argument("--features", default="2000", choices=["500", "1000", "2000", "3000", "5000", "all"])
    run.add_argument("--route", default="UMAP2_HDBSCAN_R", choices=["PCA30_SNN", "PCA30_HDBSCAN_R", "UMAP2_SNN", "UMAP2_HDBSCAN_R", "all"])
    run.add_argument("--library", default="all")
    run.add_argument("--cutoff", default="mean", choices=["none", "mean", "0.5", "all"])
    run.add_argument("--seed", type=int, default=42, help="R geometry seed; original DL model and split seeds remain 42")
    run.add_argument("--counts-layer", help="Required for h5ad: raw-count layer name, or X")
    run.add_argument("--dataset", default="user_supplied")
    run.add_argument("--rscript", default="Rscript")
    run.add_argument("--reference-r-lib")
    run.add_argument("--deg-workers", type=int, default=2)
    run.add_argument("--resume", action="store_true")
    run.add_argument("--stop-after", default="complete", choices=["prepare", "score", "complete"])
    run.add_argument("--require-slurm", action="store_true")
    args = vars(parser.parse_args(argv))
    command = args.pop("command")
    try:
        from . import runner
        result = runner.doctor(**args) if command == "doctor" else runner.run_reference(**args)
        if command == "doctor":
            print(json.dumps(result, indent=2))
        else:
            print(json.dumps({key: result[key] for key in ["status", "output"]}))
        return 0
    except (ValueError, RuntimeError, OSError, ImportError) as error:
        print(f"dgscrna: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
