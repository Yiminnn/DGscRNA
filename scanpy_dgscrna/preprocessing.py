"""
Preprocessing module for scanpy-dgscrna package.

Contains utility functions for data preprocessing and quality control.
"""

import numpy as np
import pandas as pd
import scanpy as sc
from typing import Optional, List, Union


def prepare_for_annotation(
    adata: sc.AnnData,
    min_genes: int = 200,
    normalize: bool = True,
    log_transform: bool = True,
    scale: bool = False,
    copy: bool = False
) -> Optional[sc.AnnData]:
    """
    Prepare AnnData object for cell type annotation.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    min_genes : int, default 200
        Minimum number of genes expressed per cell
    normalize : bool, default True
        Whether to normalize to 10,000 reads per cell
    log_transform : bool, default True
        Whether to log-transform the data
    scale : bool, default False
        Whether to scale the data to unit variance
    copy : bool, default False
        Return a copy instead of writing to adata
        
    Returns
    -------
    AnnData or None
        Returns AnnData if copy=True, otherwise modifies adata in place
    """
    if copy:
        adata = adata.copy()
    
    print(f"Starting with {adata.n_obs} cells and {adata.n_vars} genes")
    
    # Filter cells by minimum genes
    sc.pp.filter_cells(adata, min_genes=min_genes)
    
    print(f"After filtering: {adata.n_obs} cells and {adata.n_vars} genes")
    
    # Normalization
    if normalize:
        sc.pp.normalize_total(adata, target_sum=1e4)
    
    if log_transform:
        sc.pp.log1p(adata)
    
    if scale:
        sc.pp.scale(adata, max_value=10)
    
    return adata if copy else None


def select_hvgs(
    adata: sc.AnnData,
    n_top_genes: int = 2000,
    flavor: str = 'seurat_v3',
    copy: bool = False
) -> Optional[sc.AnnData]:
    """
    Select highly variable genes for downstream analysis.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    n_top_genes : int, default 2000
        Number of highly variable genes to select
    flavor : str, default 'seurat_v3'
        Method for HVG selection
    copy : bool, default False
        Return a copy instead of writing to adata
        
    Returns
    -------
    AnnData or None
        Returns AnnData if copy=True, otherwise modifies adata in place
    """
    if copy:
        adata = adata.copy()
    
    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, flavor=flavor)
    
    print(f"Selected {adata.var.highly_variable.sum()} highly variable genes")
    
    return adata if copy else None


def add_undecided_labels(
    adata: sc.AnnData,
    cluster_column: str,
    clusters_to_mark: Optional[List] = None,
    new_column_suffix: str = "_with_undecided",
    copy: bool = False
) -> Optional[sc.AnnData]:
    """
    Add 'Undecided' labels to specific clusters for annotation workflow.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    cluster_column : str
        Column in adata.obs containing cluster assignments
    clusters_to_mark : list, optional
        Specific clusters to mark as 'Undecided'. If None, user will be prompted
    new_column_suffix : str, default "_with_undecided"
        Suffix for the new column name
    copy : bool, default False
        Return a copy instead of writing to adata
        
    Returns
    -------
    AnnData or None
        Returns AnnData if copy=True, otherwise modifies adata in place
    """
    if copy:
        adata = adata.copy()
    
    if cluster_column not in adata.obs.columns:
        raise ValueError(f"Cluster column '{cluster_column}' not found in adata.obs")
    
    # Create new annotation column
    new_column = cluster_column + new_column_suffix
    adata.obs[new_column] = adata.obs[cluster_column].astype(str)
    
    # Mark specified clusters as undecided
    if clusters_to_mark is not None:
        for cluster in clusters_to_mark:
            mask = adata.obs[cluster_column] == cluster
            adata.obs.loc[mask, new_column] = 'Undecided'
        
        print(f"Marked clusters {clusters_to_mark} as 'Undecided' in column '{new_column}'")
    else:
        print(f"Created column '{new_column}' with original cluster labels")
        print("Available clusters:", adata.obs[cluster_column].unique().tolist())
        print("You can manually set specific clusters to 'Undecided' using:")
        print(f"adata.obs.loc[adata.obs['{cluster_column}'] == 'cluster_name', '{new_column}'] = 'Undecided'")
    
    return adata if copy else None
