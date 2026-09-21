"""Portable backend IO; importing this module performs no scientific computation."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import sys

_out = os.environ.get("DGSCRNA_REFERENCE_OUT", os.environ.get("DGSCRNA_EXAMPLE_OUT"))
OUT = Path(_out).resolve() if _out else None
FEATURES = ["hvg500", "hvg1000", "hvg2000", "hvg3000", "hvg5000", "all"]
ROUTES = ["PCA30_SNN", "PCA30_HDBSCAN_R", "UMAP2_SNN", "UMAP2_HDBSCAN_R"]


def execution_id():
    """Record the real scheduler job, or identify an explicitly local process."""
    return os.environ.get("SLURM_JOB_ID") or f"local-{os.getpid()}"


def require_slurm():
    """Enforce a site's scheduler policy when requested by the caller."""
    if os.environ.get("DGSCRNA_REQUIRE_SLURM", "0") == "1" and not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("This execution requires a SLURM allocation")
    if sys.flags.optimize:
        raise RuntimeError("Reference validation requires Python assertions; do not use -O")
    if OUT is None:
        raise RuntimeError("Set DGSCRNA_REFERENCE_OUT to the run's output directory")


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc():
    return datetime.now(timezone.utc).isoformat()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + f".part.{os.getpid()}")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    temporary.replace(path)


def complete(dest, manifest="manifest.json", flag="COMPLETE"):
    (Path(dest) / flag).write_text(sha(Path(dest) / manifest) + "\n")


def checked(dest, manifest="manifest.json", flag="COMPLETE"):
    path = Path(dest)
    return ((path / flag).is_file() and (path / manifest).is_file()
            and (path / flag).read_text().strip() == sha(path / manifest))
