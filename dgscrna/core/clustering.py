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
            sc.tl.leiden(adata, resolution=resolution, random_state=random_state)
            
        elif method == 'louvain':
            sc.tl.louvain(adata, resolution=resolution, random_state=random_state)
            
        elif method == 'hdbscan':
            # Use UMAP embeddings for HDBSCAN clustering on all genes
            if 'X_umap' not in adata.obsm:
                # First compute PCA if not available
                if 'X_pca' not in adata.obsm:
                    sc.tl.pca(adata, random_state=random_state)
                # Then compute UMAP using all genes via PCA
                sc.tl.umap(adata, random_state=random_state)
            
            # Run HDBSCAN on UMAP coordinates
            clusterer = hdbscan.HDBSCAN(
                min_cluster_size=kwargs.get('min_cluster_size', 50),
                min_samples=kwargs.get('min_samples', 5)
            )
            clusters = clusterer.fit_predict(adata.obsm['X_umap'])
            
            # Add to obs (HDBSCAN uses -1 for noise points)
            adata.obs[f'{method}_clusters'] = [f'Cluster_{i}' if i >= 0 else 'Noise' for i in clusters]
            
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
    
    # Find marker genes
    sc.tl.rank_genes_groups(
        adata,
        groupby=groupby,
        method=method,
        key_added=key_added,
        n_genes=n_genes,
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