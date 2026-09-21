"""Explicit count and marker adapters; no reference labels are used in fitting."""
from pathlib import Path
import csv
import gzip
import hashlib
import json
import os
import re


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_name(path.name + f".part.{os.getpid()}")
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    temp.replace(path)


def identifier(value, label="sample"):
    value = str(value)
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value) or value in {".", ".."}:
        raise ValueError(f"{label} must contain letters, digits, dots, hyphens or underscores, starting with a letter/digit")
    return value


def text_open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else Path(path).open()


def counts_files(path):
    path = Path(path)
    if path.is_file() and path.suffix == ".h5ad":
        return {"h5ad": path}
    if not path.is_dir():
        raise ValueError("--counts must be a 10x MatrixMarket directory or .h5ad file")
    found = {}
    for key, stems in {"matrix": ["matrix.mtx"], "cells": ["barcodes.tsv"], "genes": ["genes.tsv", "features.tsv"]}.items():
        candidates = [path / (stem + ext) for stem in stems for ext in ["", ".gz"] if (path / (stem + ext)).is_file()]
        if len(candidates) != 1:
            raise ValueError(f"Expected exactly one {key} file in {path}; found {candidates}")
        found[key] = candidates[0]
    return found


def load_markers(path):
    """Read nested JSON or long TSV/CSV without changing panel denominators."""
    path = Path(path)
    if path.suffix.lower() == ".json":
        def no_duplicates(pairs):
            value = {}
            for key, item in pairs:
                if key in value:
                    raise ValueError(f"Duplicate marker JSON key: {key}")
                value[key] = item
            return value
        libraries = json.loads(path.read_text(), object_pairs_hook=no_duplicates)
    else:
        libraries = {}
        with text_open(path) as handle:
            rows = csv.DictReader(handle, delimiter="," if path.name.endswith((".csv", ".csv.gz")) else "\t")
            if not {"library", "cell_type", "gene"}.issubset(rows.fieldnames or []):
                raise ValueError("Marker table needs library, cell_type and gene columns")
            for row in rows:
                libraries.setdefault(row["library"], {}).setdefault(row["cell_type"], []).append(row["gene"])
    if not isinstance(libraries, dict) or not libraries:
        raise ValueError("Markers must contain at least one named library")
    for library, panels in libraries.items():
        if not isinstance(library, str) or not library.strip() or not isinstance(panels, dict) or not panels:
            raise ValueError("Each marker library needs a nonempty name and panel dictionary")
        for panel, genes in panels.items():
            if not isinstance(panel, str) or not panel.strip() or panel in {"Unknown", "Undecided"}:
                raise ValueError(f"Invalid or reserved marker panel name: {panel!r}")
            if not isinstance(genes, list) or not genes or any(not isinstance(g, str) or not g.strip() or g != g.strip() for g in genes):
                raise ValueError(f"{library}/{panel} needs a nonempty list of gene names without surrounding whitespace")
    return libraries


def prepare_counts(counts, destination, sample, counts_layer=None):
    import numpy as np
    import pandas as pd
    from scipy import sparse
    from scipy.io import mmread
    files = counts_files(counts)

    def validate_entries(matrix):
        values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)
        if not np.isfinite(values).all() or (values < 0).any() or not np.array_equal(values, np.round(values)):
            raise ValueError("Input must be finite nonnegative integer counts, not normalized/scaled expression")

    if "h5ad" in files:
        if not counts_layer:
            raise ValueError("For h5ad explicitly set --counts-layer counts (or X for raw counts in X)")
        try:
            import anndata
        except ImportError as error:
            raise ValueError("h5ad input requires pip install 'dgscrna[anndata]'") from error
        adata = anndata.read_h5ad(files["h5ad"])
        if counts_layer != "X" and counts_layer not in adata.layers:
            raise ValueError(f"Missing count layer {counts_layer!r}")
        matrix = adata.X if counts_layer == "X" else adata.layers[counts_layer]
        validate_entries(matrix)
        x = sparse.csr_matrix(matrix)
        cells, genes = list(map(str, adata.obs_names)), list(map(str, adata.var_names))
    else:
        if counts_layer is not None:
            raise ValueError("--counts-layer applies only to h5ad input")
        with text_open(files["cells"]) as handle:
            cells = [row[0] for row in csv.reader(handle, delimiter="\t") if row]
        with text_open(files["genes"]) as handle:
            features = [row for row in csv.reader(handle, delimiter="\t") if row]
        genes = [row[1] if len(row) >= 2 else row[0] for row in features]
        if any(len(row) >= 3 and row[2] != "Gene Expression" for row in features):
            raise ValueError("Multi-assay 10x inputs must first be restricted to Gene Expression features")
        matrix = mmread(files["matrix"])
        validate_entries(matrix)
        x = sparse.csr_matrix(matrix).T.tocsr()
    if x.shape != (len(cells), len(genes)):
        raise ValueError(f"Count dimensions {x.shape} disagree with cell/gene identifiers")
    for label, values in [("cell", cells), ("gene", genes)]:
        if len(values) != len(set(values)) or any(not v or v != v.strip() or "\n" in v or "\r" in v for v in values):
            raise ValueError(f"{label} identifiers must be nonempty and unique, without whitespace at the edges")
    if not np.isfinite(x.data).all() or (x.data < 0).any() or not np.array_equal(x.data, np.round(x.data)):
        raise ValueError("Input must be finite nonnegative integer counts, not normalized/scaled expression")
    x.sum_duplicates(); x.eliminate_zeros()
    keep = np.asarray((x > 0).sum(axis=0)).ravel() >= 3
    if (x.data > 2**24).any():
        raise ValueError("Counts exceed the exact integer precision limit of the historical float32 export (16777216)")
    x = x[:, keep].astype(np.float32).astype(np.float64).tocsr()
    x.sum_duplicates(); x.eliminate_zeros(); x.sort_indices()
    if min(x.shape) <= 31:
        raise ValueError("The fixed PCA30 reference requires more than 31 cells and retained genes")
    if x.nnz > np.iinfo(np.int32).max:
        raise ValueError("Count matrix exceeds the R sparse input index limit; split the analysis explicitly")
    if not np.isfinite(x.data).all() or (np.asarray(x.sum(axis=1)).ravel() <= 0).any():
        raise ValueError("Every cell must have finite positive counts after the >=3-cell gene filter")
    dest = Path(destination)
    dest.mkdir(parents=True, exist_ok=False)
    x.data.astype("<f8").tofile(dest / "x.bin")
    x.indices.astype("<i4").tofile(dest / "i.bin")
    x.indptr.astype("<i4").tofile(dest / "p.bin")
    pd.DataFrame({"cell_id": cells, "batch": sample}).to_csv(dest / "cells_fit.csv", index=False)
    pd.DataFrame({"gene": np.asarray(genes)[keep]}).to_csv(dest / "genes.csv", index=False)
    manifest = {
        "status": "completed", "sample": sample, "n_cells": x.shape[0], "n_genes": x.shape[1], "nnz": x.nnz,
        "n_genes_before_filter": len(genes), "input_semantics": "already-QC nonnegative integer RNA counts; genes detected in >=3 cells",
        "fitting_files": {name: sha(dest / name) for name in ["x.bin", "i.bin", "p.bin", "cells_fit.csv", "genes.csv"]},
        "raw_input_files": {str(p): sha(p) for p in files.values()},
        "counts_layer": counts_layer, "truth_labels_excluded_from_fit": True,
        "source_sha256": sha(__file__), "job": os.environ.get("SLURM_JOB_ID", "local"),
    }
    write_json(dest / "input_manifest.json", manifest)
    (dest / "INPUT_COMPLETE").write_text(sha(dest / "input_manifest.json") + "\n")
    return manifest


def marker_coverage(libraries, genes, requested="all"):
    """Reject incompatible input namespaces without changing density arithmetic."""
    genes = set(genes)
    coverage = {}
    for name, panels in libraries.items():
        union = {gene for panel in panels.values() for gene in panel}
        hit = union & genes
        coverage[name] = {"unique_marker_genes": len(union), "retained_marker_genes": len(hit),
                          "panels": {panel: {"denominator": len(markers), "retained_unique": len(set(markers) & genes)}
                                     for panel, markers in panels.items()}}
        if (requested == "all" or name == requested) and not hit:
            raise ValueError(f"Marker library {name!r} has zero overlap with retained gene identifiers; check species and gene namespace")
    return coverage
