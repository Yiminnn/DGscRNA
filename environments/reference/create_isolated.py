"""Create the two pinned reference prefixes without user/site Conda settings."""
from pathlib import Path
import argparse
import json
import os
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prefix_root", type=Path)
    args = parser.parse_args()
    root = args.prefix_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    if any((root / name).exists() for name in ("r", "python", "package_cache")):
        parser.error("Use a fresh prefix root: r, python and package_cache must not exist.")
    here = Path(__file__).resolve().parent
    cache = root / "package_cache"
    os.environ["CONDA_PKGS_DIRS"] = str(cache)
    os.environ["CONDA_SAFETY_CHECKS"] = "enabled"
    os.environ["CONDA_ALWAYS_COPY"] = "true"
    # Run with Conda base's Python, whose environment provides the conda API.
    import conda
    from conda.base.context import context, reset_context
    from conda.cli.python_api import Commands, run_command

    # Environment variables alone do not exclude merged user/site cache lists.
    reset_context(search_path=())
    if tuple(context.pkgs_dirs) != (str(cache),) or not context.always_copy:
        raise RuntimeError("Conda did not accept the isolated cache configuration")
    record = dict(
        conda_version=conda.__version__, conda_python=sys.executable,
        configuration_search_path=[], package_cache=str(cache),
        copy_install=True, safety_checks=str(context.safety_checks),
        job=os.environ.get("SLURM_JOB_ID", "local"),
    )
    (root / "installer_configuration.json").write_text(json.dumps(record, indent=2) + "\n")
    for runtime in ("r", "python"):
        _, _, status = run_command(
            Commands.CREATE, "--copy", "--prefix", str(root / runtime),
            "--file", str(here / f"{runtime}-linux-64.explicit.txt"),
            search_path=(), stdout=None, stderr=None, use_exception_handler=False,
        )
        if status != 0:
            raise RuntimeError(f"Conda installation failed for {runtime}: {status}")
        if tuple(context.pkgs_dirs) != (str(cache),):
            raise RuntimeError("Conda changed the effective package cache")
    print(f"ISOLATED_PREFIXES_INSTALLED {root}", flush=True)


if __name__ == "__main__":
    main()
