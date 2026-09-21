"""One public runner around the original numerical implementation."""
from datetime import datetime, timezone
from importlib import metadata
from pathlib import Path
import json
import os
import shutil
import subprocess
import sys
import uuid

from .. import __version__
from .io import counts_files, identifier, load_markers, marker_coverage, prepare_counts, sha, write_json

HERE = Path(__file__).resolve().parent
BACKEND = HERE / "backend"
ROUTES = ("PCA30_SNN", "PCA30_HDBSCAN_R", "UMAP2_SNN", "UMAP2_HDBSCAN_R")
BUDGETS = ("500", "1000", "2000", "3000", "5000", "all")


def source_hashes():
    return {str(p.relative_to(HERE)): sha(p) for p in sorted(HERE.rglob("*")) if p.is_file() and p.suffix in {".py", ".R", ".json"}}


def runtime_versions():
    return {name: metadata.version(name) for name in ["numpy", "scipy", "pandas", "torch"]}


def doctor(rscript="Rscript", reference_r_lib=None):
    executable = shutil.which(str(rscript))
    if not executable:
        raise ValueError("Rscript was not found. Install the pinned reference environment and pass --rscript /path/to/Rscript")
    env = os.environ.copy()
    env.pop("R_LIBS", None); env.pop("R_LIBS_USER", None); env.pop("R_LIBS_SITE", None)
    env.update(R_LIBS_USER="", R_LIBS_SITE="", R_ENVIRON_USER=os.devnull, R_PROFILE_USER=os.devnull)
    if reference_r_lib:
        env["DGSCRNA_REFERENCE_R_LIB"] = str(Path(reference_r_lib).resolve())
    else:
        env.pop("DGSCRNA_REFERENCE_R_LIB", None)
    expression = """
      x<-Sys.getenv('DGSCRNA_REFERENCE_R_LIB');.libPaths(c(if(nzchar(x))x,.Library),include.site=FALSE);
      pkgs<-c('Seurat','SeuratObject','Matrix','jsonlite','digest','future','dbscan','limma','uwot');
      missing<-pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)];
      if(length(missing))stop(paste('Missing R packages:',paste(missing,collapse=', ')));
      cat(jsonlite::toJSON(list(R=as.character(getRversion()),libraries=.libPaths(),
        packages=setNames(lapply(pkgs,function(x)as.character(packageVersion(x))),pkgs),
        package_paths=setNames(lapply(pkgs,function(x)normalizePath(find.package(x))),pkgs)),auto_unbox=TRUE));
    """
    result = subprocess.run([executable, "--vanilla", "-e", expression], env=env, text=True, capture_output=True)
    if result.returncode:
        raise ValueError("R environment check failed:\n" + result.stderr[-4000:])
    return {"python": sys.version.split()[0], "python_packages": runtime_versions(),
            "DL_threads": min(4, int(os.environ.get("SLURM_CPUS_PER_TASK", "4"))),
            "rscript": str(Path(executable).resolve()),
            "R": json.loads(result.stdout), "reference_r_lib": str(Path(reference_r_lib).resolve()) if reference_r_lib else None}


def _inventory(directory, root):
    return {str(p.relative_to(root)): sha(p) for p in sorted(directory.rglob("*"))
            if p.is_file() and not any(part.startswith(".part") for part in p.parts)}


def _checked_receipt(path, root, plan_hash):
    if not path.exists():
        return False
    receipt = json.loads(path.read_text())
    if receipt["plan_sha256"] != plan_hash:
        raise ValueError(f"Checkpoint configuration changed: {path}")
    for name, expected in receipt["artifacts"].items():
        target = root / name
        if not target.is_file() or sha(target) != expected:
            raise ValueError(f"Checkpoint is missing or modified: {target}; preserve this run and use a new output directory")
    return True


def run_reference(*, counts, markers, out, sample=None, features="2000", route="UMAP2_HDBSCAN_R",
                  library="all", cutoff="mean", preset="gbm-reference", seed=42, counts_layer=None,
                  rscript="Rscript", reference_r_lib=None, deg_workers=2, resume=False,
                  stop_after="complete", require_slurm=False, dataset="user_supplied"):
    """Run the GBM single-sample reference, returning its run manifest.

    Inputs are already-QC nonnegative integer counts (10x directory or explicit
    h5ad count layer) and nested JSON/long TSV markers. No truth labels enter fit.
    ``features`` is 500/1000/2000/3000/5000/all; ``route``, ``library`` and
    ``cutoff`` accept 'all'. ``resume`` verifies all completed stage artifacts.
    Use ``stop_after='prepare'`` or 'score' to split HPC work into stages.
    """
    if require_slurm and not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("--require-slurm needs an active SLURM allocation")
    if sys.flags.optimize:
        raise ValueError("The reference backend requires Python assertions; do not use python -O")
    if preset != "gbm-reference":
        raise ValueError("The public counts runner currently supports gbm-reference; PTC requires its separately validated grouped-input workflow")
    features = str(features).removeprefix("hvg")
    if features not in BUDGETS or route not in (*ROUTES, "all") or cutoff not in ("none", "mean", "0.5", "all"):
        raise ValueError("Unsupported feature budget, route or cutoff")
    if stop_after not in ("prepare", "score", "complete") or not 1 <= int(deg_workers) <= 4:
        raise ValueError("stop_after must be prepare/score/complete and deg_workers must be 1–4")
    if isinstance(seed, bool) or not isinstance(seed, int) or not 0 <= seed < 2**31:
        raise ValueError("seed must be in [0, 2^31)")
    counts, markers, out = Path(counts).resolve(), Path(markers).resolve(), Path(out).resolve()
    sample = identifier(sample or counts.stem)
    marker_hash = sha(markers)
    libraries = load_markers(markers)
    if sha(markers) != marker_hash:
        raise ValueError("Marker input changed while it was being read")
    if library != "all" and library not in libraries:
        raise ValueError(f"Unknown library {library!r}; available: {', '.join(libraries)}")
    raw_files = counts_files(counts)
    if counts.is_relative_to(out) or markers.is_relative_to(out):
        raise ValueError("Output must not contain input counts or marker files")
    runtime = doctor(rscript, reference_r_lib)
    budget = "all" if features == "all" else "hvg" + features
    condition = budget if seed == 42 else f"{budget}_seed{seed}"
    config = dict(preset=preset, sample=sample, dataset=dataset, features=features, route=route,
                  library=library, cutoff=cutoff, seed=int(seed), counts_layer=counts_layer,
                  deg_workers=int(deg_workers), counts=str(counts), markers=str(markers))
    plan = dict(package_version=__version__, config=config, sources=source_hashes(), runtime=runtime,
                input_sha256={**{str(p): sha(p) for p in raw_files.values()}, str(markers): marker_hash},
                scientific_contract="R reference preparation/scoring + historical Python DL; already-QC single sample")
    if out.exists():
        if not resume or not (out / "run_config.json").is_file():
            raise ValueError("Output already exists; use a new directory, or --resume for a matching packaged run")
        if json.loads((out / "run_config.json").read_text()) != plan:
            raise ValueError("Input, configuration, source or environment changed; cannot resume this output")
    else:
        out.mkdir(parents=True)
        write_json(out / "run_config.json", plan)
    lock = out / ".run.lock"
    try:
        lock.mkdir()
    except FileExistsError as error:
        raise ValueError(f"Run lock exists: {lock}. Check the recorded process/job before recovering an interrupted run") from error
    write_json(lock / "owner.json", {"pid": os.getpid(), "job": os.environ.get("SLURM_JOB_ID"), "host": os.uname().nodename})
    plan_hash = sha(out / "run_config.json")
    started = datetime.now(timezone.utc).isoformat()
    prep = out / "GBM" / sample / condition
    env = os.environ.copy()
    for key in ["PYTHONPATH", "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE", "DGSCRNA_DL_CACHE_ROOT", "DGSCRNA_SITE_ROOT", "DGSCRNA_REFERENCE_R_LIB"]:
        env.pop(key, None)
    env.update(DGSCRNA_REFERENCE_OUT=str(out), DGSCRNA_EXAMPLE_OUT=str(out),
               DGSCRNA_RSCRIPT=runtime["rscript"], DGSCRNA_PYTHON=sys.executable,
               DGSCRNA_REQUIRE_SLURM="1" if require_slurm else "0", DGSCRNA_DEG_WORKERS=str(deg_workers),
               R_LIBS_USER="", R_LIBS_SITE="", R_ENVIRON_USER=os.devnull, R_PROFILE_USER=os.devnull,
               PYTHONNOUSERSITE="1", OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1", NUMBA_NUM_THREADS="1")
    if reference_r_lib:
        env["DGSCRNA_REFERENCE_R_LIB"] = str(Path(reference_r_lib).resolve())
    env["DGSCRNA_ONLY_ROUTE"] = "" if route == "all" else route
    (out / "logs").mkdir(exist_ok=True)
    (out / "checkpoints").mkdir(exist_ok=True)

    def command(args, name):
        log = out / "logs" / f"{name}_{uuid.uuid4().hex[:8]}.log"
        print(f"DG-scRNA: {name}; log={log}", flush=True)
        with log.open("w") as handle:
            result = subprocess.run(args, env=env, stdout=handle, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"{name} failed (exit {result.returncode}); see {log}")

    def checkpoint(name, directories):
        artifacts = {}
        for directory in directories:
            artifacts.update(_inventory(directory, out))
        write_json(out / "checkpoints" / (name + ".json"), {"stage": name, "plan_sha256": plan_hash, "artifacts": artifacts})

    def checked(name):
        return _checked_receipt(out / "checkpoints" / (name + ".json"), out, plan_hash)

    try:
        if not checked("input"):
            inp = out / "inputs" / sample
            if inp.exists():
                inp.rename(inp.with_name(sample + ".interrupted." + uuid.uuid4().hex[:8]))
            im = prepare_counts(counts, inp, sample, counts_layer)
            if any(plan["input_sha256"][name] != value for name, value in im["raw_input_files"].items()):
                raise ValueError("Count inputs changed while being read; start a fresh run after inputs stabilize")
            im["dataset"] = dataset
            write_json(inp / "input_manifest.json", im)
            (inp / "INPUT_COMPLETE").write_text(sha(inp / "input_manifest.json") + "\n")
            write_json(out / "markers" / "libraries.json", libraries)
            import csv
            with (inp / "genes.csv").open() as stream:
                genes = [row["gene"] for row in csv.DictReader(stream)]
            write_json(out / "markers/coverage.json", marker_coverage(libraries, genes, library))
            checkpoint("input", [inp, out / "markers"])
        if not checked("prepare"):
            if prep.exists():
                prep.rename(prep.with_name(condition + ".interrupted." + uuid.uuid4().hex[:8]))
            command([runtime["rscript"], "--vanilla", str(BACKEND / "prepare_R.R"), sample, budget, str(seed)], "prepare")
            import csv
            with (prep / "Seurat_gene_names.csv").open() as stream:
                names = list(csv.DictReader(stream))
            coverage = marker_coverage(libraries, [row["Seurat"] for row in names], library)
            write_json(prep / "marker_coverage.json", coverage)
            changed = sum(row["source"] != row["Seurat"] for row in names)
            if changed:
                print(f"DG-scRNA: Seurat renamed {changed} gene identifiers; inspect {prep / 'Seurat_gene_names.csv'} and marker_coverage.json", flush=True)
            checkpoint("prepare", [prep])
        if stop_after == "prepare":
            return {"status": "prepared", "output": str(out), "prepare": str(prep)}
        routes = ROUTES if route == "all" else [route]
        if not checked("score"):
            for r in routes:
                directory = prep / r
                if directory.exists():
                    directory.rename(directory.with_name(r + ".interrupted." + uuid.uuid4().hex[:8]))
            for r in routes:
                # Independent processes match the original campaign and prevent
                # clustering/RNG state from one route affecting another route.
                env["DGSCRNA_ONLY_ROUTE"] = r
                command([runtime["rscript"], "--vanilla", str(BACKEND / "score_R.R"), sample, str(prep)], f"score_{r}")
            checkpoint("score", [prep / r for r in routes])
        if stop_after == "score":
            return {"status": "scored", "output": str(out), "prepare": str(prep)}
        if checked("terminal"):
            expected = sha(out / "run_manifest.json")
            flag = out / "COMPLETE"
            if flag.exists() and flag.read_text().strip() != expected:
                raise ValueError("Completion marker does not match the verified run manifest")
            flag.write_text(expected + "\n")
            return json.loads((out / "run_manifest.json").read_text())
        records = []
        # Share only exact-input fits within this fresh attempt; never import an
        # earlier campaign's predictions or a partially verified training cache.
        env["DGSCRNA_DL_CACHE_ROOT"] = str(out / "DL_cache" / uuid.uuid4().hex)
        for r in routes:
            scores = json.loads((prep / r / "score_manifest.json").read_text())
            arms = [(aid, arm) for aid, arm in scores["arms"].items()
                    if (library == "all" or arm["library"] == library) and (cutoff == "all" or arm["cutoff"] == cutoff)]
            if not arms:
                raise ValueError("No scoring arms matched the requested marker/cutoff")
            for aid, arm in arms:
                terminal = prep / r / "terminal" / aid
                # Partial terminal work is not certified; retain it and refit it.
                if terminal.exists():
                    terminal.rename(terminal.with_name(aid + ".interrupted." + uuid.uuid4().hex[:8]))
            if len(arms) == len(scores["arms"]):
                command([sys.executable, "-s", str(BACKEND / "terminal.py"), str(prep / r), "all"], f"refine_{r}_all")
            else:
                for aid, arm in arms:
                    command([sys.executable, "-s", str(BACKEND / "terminal.py"), str(prep / r), aid], f"refine_{r}_{aid}")
            for aid, arm in arms:
                terminal = prep / r / "terminal" / aid
                tm = json.loads((terminal / "terminal_manifest.json").read_text())
                if not tm["terminal_valid"] or sha(terminal / "predictions.csv.gz") != tm["predictions_sha256"]:
                    raise ValueError(f"Invalid terminal result: {terminal}")
                records.append({"route": r, "arm": aid, "library": arm["library"], "cutoff": arm["cutoff"],
                                "dl_status": tm["dl_status"], "training_executed": tm["training_executed"],
                                "predictions": str((terminal / "predictions.csv.gz").relative_to(out)),
                                "predictions_sha256": tm["predictions_sha256"]})
        _export(out, prep, records)
        manifest = {"status": "completed", "package_version": __version__, "plan_sha256": plan_hash,
                    "started_at": started, "completed_at": datetime.now(timezone.utc).isoformat(),
                    "output": str(out), "prepare": str(prep), "conditions": records,
                    "terminal_condition_count": len(records), "trained_condition_count": sum(r["training_executed"] for r in records),
                    "reference_labels_used_for_fit": False, "job": os.environ.get("SLURM_JOB_ID", "local"),
                    "exports": {name: sha(out / name) for name in ["annotations.csv.gz", "embedding.csv", "condition_summary.csv"]}}
        write_json(out / "run_manifest.json", manifest)
        receipt = {"stage": "terminal", "plan_sha256": plan_hash, "artifacts": {}}
        for r in routes:
            receipt["artifacts"].update(_inventory(prep / r / "terminal", out))
        receipt["artifacts"].update({name: sha(out / name) for name in ["run_manifest.json", *manifest["exports"]]})
        write_json(out / "checkpoints/terminal.json", receipt)
        (out / "COMPLETE").write_text(sha(out / "run_manifest.json") + "\n")
        print(f"DG-scRNA complete: {out / 'annotations.csv.gz'}", flush=True)
        return manifest
    except Exception as error:
        write_json(out / f"failure_{uuid.uuid4().hex[:8]}.json", {"error": str(error), "job": os.environ.get("SLURM_JOB_ID", "local"), "plan_sha256": plan_hash})
        raise
    finally:
        shutil.rmtree(lock)


def _export(out, prep, records):
    import gzip
    import pandas as pd
    first = True
    with gzip.open(out / "annotations.csv.gz", "wt") as stream:
        for record in records:
            frame = pd.read_csv(out / record["predictions"], keep_default_na=False, dtype=str)
            clusters = pd.read_csv(prep / record["route"] / "clusters.csv", keep_default_na=False, dtype=str)
            if frame.cell_id.tolist() != clusters.cell_id.tolist():
                raise ValueError("Cell order changed between clustering and terminal annotation")
            frame["cluster"] = clusters.cluster
            for name in ["route", "arm", "library", "cutoff", "dl_status", "training_executed"]:
                frame[name] = record[name]
            frame.to_csv(stream, index=False, header=first)
            first = False
    coordinates = pd.read_csv(prep / "UMAP2.csv", index_col=0, dtype=str, keep_default_na=False)
    coordinates.index.name = "cell_id"
    coordinates.to_csv(out / "embedding.csv")
    pd.DataFrame(records).to_csv(out / "condition_summary.csv", index=False)
