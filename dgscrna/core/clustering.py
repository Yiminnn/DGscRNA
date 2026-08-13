"""
Clustering module for DGscRNA package
"""

import scanpy as sc
import numpy as np
import pandas as pd
from typing import List, Optional, Union, Dict
from sklearn.cluster import KMeans
import hdbscan
import warnings
warnings.filterwarnings('ignore')

def run_clustering(
    adata,
    methods: List[str] = ['leiden', 'hdbscan', 'kmeans'],
    resolution: float = 0.5,
    n_neighbors: int = 15,
    n_clusters: Optional[int] = None,
    random_state: int = 42,
    **kwargs
):
    """
    Run multiple clustering algorithms on the data
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    methods : List[str], default=['leiden', 'hdbscan', 'kmeans']
        List of clustering methods to run
    resolution : float, default=0.5
        Resolution parameter for Leiden clustering
    n_neighbors : int, default=15
        Number of neighbors for neighborhood graph
    n_clusters : int, optional
        Number of clusters for K-means (if None, estimated from data)
    **kwargs : cluster_space {'umap','pca'} selects the HDBSCAN embedding
        (default 'umap'); min_cluster_size (default 50) and min_samples
        (defaults to min_cluster_size, matching dbscan::hdbscan's minPts)
    random_state : int, default=42
        Random state for reproducibility
    **kwargs
        Additional arguments for clustering methods
        
    Returns
    -------
    AnnData
        AnnData object with clustering results added to obs
    """
    
    # Ensure neighborhood graph exists
    if 'neighbors' not in adata.uns:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, random_state=random_state)
    
    for method in methods:
        print(f"Running {method} clustering...")
        
        if method == 'leiden':
            sc.tl.leiden(adata, resolution=resolution, random_state=random_state, key_added='leiden_clusters')
            
        elif method == 'louvain':
            sc.tl.louvain(adata, resolution=resolution, random_state=random_state, key_added='louvain_clusters')

        elif method == 'hdbscan':
            # R reference (examples/R/source.R, clustering_strategies) runs
            # dbscan::hdbscan(minPts = 50) on the PCA embedding and again on the
            # UMAP embedding. dbscan's single minPts sets BOTH the minimum
            # cluster size and the core-distance neighbour count, so min_samples
            # defaults to min_cluster_size here rather than to 5.
            if 'X_pca' not in adata.obsm:
                sc.tl.pca(adata, random_state=random_state)
            space = kwargs.get('cluster_space', 'umap')
            if space == 'umap' and 'X_umap' not in adata.obsm:
                sc.tl.umap(adata, random_state=random_state)
            rep = 'X_umap' if space == 'umap' else 'X_pca'
            
            min_cluster_size = kwargs.get('min_cluster_size', 50)
            clusterer = hdbscan.HDBSCAN(
                min_cluster_size=min_cluster_size,
                min_samples=kwargs.get('min_samples', min_cluster_size)
            )
            clusters = clusterer.fit_predict(np.asarray(adata.obsm[rep], dtype=np.float64))
            
            # HDBSCAN labels noise -1. It is kept as an explicit 'Noise' group so
            # it is never silently dropped, and downstream it is excluded from
            # differential expression and from cell-type assignment.
            adata.obs[f'{method}_clusters'] = [f'Cluster_{i}' if i >= 0 else 'Noise' for i in clusters]
            adata.uns[f'{method}_noise_fraction'] = float(np.mean(np.asarray(clusters) == -1))
            
        elif method == 'kmeans':
            # Use PCA embeddings for K-means
            if 'X_pca' not in adata.obsm:
                sc.tl.pca(adata, random_state=random_state)
            
            # Estimate number of clusters if not provided
            if n_clusters is None:
                # Simple heuristic: sqrt of number of cells
                n_clusters = int(np.sqrt(adata.n_obs))
                n_clusters = max(2, min(n_clusters, 20))  # Between 2 and 20
            
            # Run K-means (filter out non-KMeans kwargs)
            kmeans_kwargs = {k: v for k, v in kwargs.items() 
                           if k in ['init', 'n_init', 'max_iter', 'tol', 'algorithm']}
            kmeans = KMeans(n_clusters=n_clusters, random_state=random_state, **kmeans_kwargs)
            clusters = kmeans.fit_predict(adata.obsm['X_pca'])
            
            adata.obs[f'{method}_clusters'] = [f'Cluster_{i}' for i in clusters]
            
        else:
            print(f"Warning: Clustering method '{method}' not supported, skipping...")
    
    return adata

def find_markers(
    adata,
    groupby: str,
    method: str = 'wilcoxon',
    key_added: str = 'rank_genes_groups',
    n_genes: int = 100,
    layer: Optional[str] = 'lognorm',
    exclude_groups: Optional[List[str]] = ('Noise',),
    **kwargs
):
    """
    Find differentially expressed genes for each cluster
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    groupby : str
        Key in obs for grouping
    method : str, default='wilcoxon'
        Method for differential expression testing
    key_added : str, default='rank_genes_groups'
        Key to store results in uns
    n_genes : int, default=100
        Number of top genes to return per cluster
    layer : str, optional, default='lognorm'
        Layer holding log-normalised counts. The R reference implementation
        (examples/R/source.R, find_markers_on_id) runs FindAllMarkers on the
        log-normalised assay. If `adata.X` has been z-scored by
        `sc.pp.scale`, differential expression on `.X` yields NaN log
        fold-changes for every gene, which silently drives every marker score
        to zero. Falls back to `.X` when the layer is absent.
    exclude_groups : list of str, optional, default=('Noise',)
        Groups to drop before testing. HDBSCAN noise is not a cell population,
        so it gets neither differential expression nor a cell-type call.
    **kwargs
        Additional arguments for sc.tl.rank_genes_groups
        
    Returns
    -------
    AnnData
        AnnData object with marker genes results
    """
    
    # Check if groupby exists in obs
    if groupby not in adata.obs.columns:
        raise ValueError(f"Groupby key '{groupby}' not found in adata.obs")
    
    # Set the grouping
    adata.obs[groupby] = adata.obs[groupby].astype('category')
    
    # Resolve the expression layer. Testing on z-scored .X produces NaN logFCs.
    use_layer = layer if (layer is not None and layer in adata.layers) else None
    if layer is not None and use_layer is None:
        warnings.warn(
            f"layer '{layer}' not found; running differential expression on .X. "
            "If .X has been scaled, log fold-changes will be NaN and all marker "
            "scores will be zero.",
            RuntimeWarning,
        )
    
    # Groups that are not cell populations (HDBSCAN noise) are not tested
    groups = 'all'
    if exclude_groups:
        present = list(adata.obs[groupby].cat.categories)
        keep = [g for g in present if g not in set(exclude_groups)]
        if len(keep) < len(present):
            groups = keep
    
    # Find marker genes
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        key_added=key_added,
        n_genes=n_genes,
        layer=use_layer,
        use_raw=False,
        groups=groups,
        **kwargs
    )
    
    return adata

def get_marker_genes(
    adata,
    groupby: str,
    key: str = 'rank_genes_groups',
    n_genes: int = 100,
    pval_cutoff: float = 0.05,
    logfc_cutoff: float = 0.25
):
    """
    Extract marker genes from rank_genes_groups results
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    groupby : str
        Key in obs for grouping
    key : str, default='rank_genes_groups'
        Key in uns containing rank_genes_groups results
    n_genes : int, default=100
        Number of top genes per cluster
    pval_cutoff : float, default=0.05
        P-value cutoff for significance
    logfc_cutoff : float, default=0.25
        Log fold change cutoff
        
    Returns
    -------
    Dict
        Dictionary with cluster names as keys and marker gene lists as values
    """
    
    if key not in adata.uns:
        raise ValueError(f"Key '{key}' not found in adata.uns. Run find_markers first.")
    
    results = adata.uns[key]
    marker_genes = {}
    
    # Get cluster names
    cluster_names = results['names'].dtype.names
    
    for cluster in cluster_names:
        # Get genes, scores, pvals, and logfoldchanges
        genes = results['names'][cluster][:n_genes]
        scores = results['scores'][cluster][:n_genes]
        pvals = results['pvals_adj'][cluster][:n_genes]
        logfcs = results['logfoldchanges'][cluster][:n_genes]
        
        # Filter by significance
        significant = (pvals < pval_cutoff) & (logfcs > logfc_cutoff)
        marker_genes[cluster] = genes[significant].tolist()
    
    return marker_genes 