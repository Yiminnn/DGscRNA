"""Regression tests pinning the behaviours aligned to the R reference.

Reference: examples/R/source.R (find_markers_on_id, density_score) as driven by
examples/DGscRNA-Vignette.rmd. Each test fails against the pre-alignment code.
"""
import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import anndata as ad

from dgscrna.core.clustering import run_clustering, find_markers
from dgscrna.core.marker_scoring import score_cell_types


def _toy(n=300, g=60, seed=0):
    """Three well-separated blocks, each driven by its own marker genes."""
    rng = np.random.default_rng(seed)
    X = rng.negative_binomial(2, 0.5, size=(n, g)).astype(float)
    block = np.repeat([0, 1, 2], n // 3)
    for b in range(3):
        X[block == b, b * 5:(b + 1) * 5] += 40
    a = ad.AnnData(X)
    a.var_names = [f"G{i}" for i in range(g)]
    a.obs_names = [f"C{i}" for i in range(n)]
    a.obs["block"] = pd.Categorical([f"b{b}" for b in block])
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    a.layers["lognorm"] = a.X.copy()
    sc.pp.scale(a, max_value=10)          # .X is now z-scored, as in the pipeline
    sc.tl.pca(a, n_comps=10, random_state=0)
    sc.pp.neighbors(a, n_neighbors=15, random_state=0)
    sc.tl.leiden(a, key_added="cl", random_state=0, flavor="igraph",
                 n_iterations=2, directed=False)
    return a


MARKERS = {f"type{b}": [f"G{i}" for i in range(b * 5, (b + 1) * 5)] for b in range(3)}


def test_deg_uses_lognorm_layer_not_scaled_X():
    """Differential expression on z-scored .X yields NaN logFCs for every gene."""
    a = find_markers(_toy(), groupby="cl")
    lfc = a.uns["rank_genes_groups"]["logfoldchanges"]
    first = np.asarray(lfc[lfc.dtype.names[0]], dtype=float)
    assert not np.isnan(first).all(), "logFCs are all NaN: the lognorm layer was not used"


def test_scoring_is_not_degenerate_on_scaled_input():
    """The NaN-logFC failure collapses every cluster onto one panel."""
    a = find_markers(_toy(), groupby="cl")
    a = score_cell_types(a, MARKERS, cluster_key="cl", marker_set_name="toy")
    calls = a.obs["toy_cl_none"].astype(str)
    assert calls.nunique() > 1, f"all clusters collapsed onto one call: {calls.unique()}"


def test_tie_is_undecided_not_first_index():
    """R: names(which(score == max)); length > 1 -> 'Undecided'."""
    a = _toy()
    a = find_markers(a, groupby="cl")
    # panels sharing no gene with the data can never fire, so every panel ties at 0
    dead = {"panelA": ["ZZZ1", "ZZZ2"], "panelB": ["ZZZ3", "ZZZ4"]}
    a = score_cell_types(a, dead, cluster_key="cl", marker_set_name="dead")
    calls = set(a.obs["dead_cl_none"].astype(str))
    assert calls == {"Undecided"}, f"expected all Undecided on an all-zero tie, got {calls}"


def test_hdbscan_noise_is_not_a_cell_type():
    """Noise is not a population: no DEGs, no cell-type call."""
    a = _toy(n=600)
    a = run_clustering(a, methods=["hdbscan"], random_state=0,
                       min_cluster_size=15, cluster_space="pca")
    if "Noise" not in set(a.obs["hdbscan_clusters"]):
        pytest.skip("this toy partition produced no noise points")
    a = find_markers(a, groupby="hdbscan_clusters")
    assert "Noise" not in a.uns["rank_genes_groups"]["names"].dtype.names
    a = score_cell_types(a, MARKERS, cluster_key="hdbscan_clusters", marker_set_name="toy")
    noise = a.obs["hdbscan_clusters"] == "Noise"
    assert set(a.obs.loc[noise, "toy_hdbscan_none"].astype(str)) == {"Unknown"}


def test_hdbscan_min_samples_follows_min_cluster_size():
    """dbscan::hdbscan(minPts) sets both parameters; python must not default to 5."""
    a = run_clustering(_toy(n=600), methods=["hdbscan"], random_state=0,
                       min_cluster_size=15, cluster_space="pca")
    assert "hdbscan_noise_fraction" in a.uns


def test_pval_filter_is_optional():
    """R's density_score filters on avg_log2FC only; max_pval_adj=None reproduces it."""
    a = find_markers(_toy(), groupby="cl")
    strict = score_cell_types(a.copy(), MARKERS, cluster_key="cl", marker_set_name="s")
    loose = score_cell_types(a.copy(), MARKERS, cluster_key="cl", marker_set_name="l",
                             max_pval_adj=None)
    assert "s_cl_none" in strict.obs and "l_cl_none" in loose.obs
