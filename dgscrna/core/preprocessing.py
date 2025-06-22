"""
Preprocessing module for DGscRNA package
"""

import scanpy as sc
import numpy as np
import pandas as pd
from typing import List, Optional, Union
import warnings
warnings.filterwarnings('ignore')

def preprocess_adata(
    adata,
    min_genes: int = 200,
    max_genes: Optional[int] = None,
    min_cells: int = 3,
    max_counts: Optional[int] = None,
    max_mito_pct: float = 15.0,
    normalize: bool = True,
    scale: bool = True,
    n_pcs: int = 30,
    n_neighbors: int = 15,
    random_state: int = 42
):
    """
    Preprocess AnnData object with quality control and normalization
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    min_genes : int, default=200
        Minimum number of genes expressed per cell
    max_genes : int, optional
        Maximum number of genes expressed per cell
    min_cells : int, default=3
        Minimum number of cells expressing a gene
    max_counts : int, optional
        Maximum number of counts per cell
    max_mito_pct : float, default=15.0
        Maximum percentage of mitochondrial genes
    normalize : bool, default=True
        Whether to normalize data
    scale : bool, default=True
        Whether to scale data
    n_pcs : int, default=30
        Number of principal components
    n_neighbors : int, default=15
        Number of neighbors for neighborhood graph
    random_state : int, default=42
        Random state for reproducibility
        
    Returns
    -------
    AnnData
        Preprocessed AnnData object
    """
    
    # Calculate QC metrics
    sc.pp.calculate_qc_metrics(adata)
    
    # Filter cells
    sc.pp.filter_cells(adata, min_genes=min_genes)
    if max_genes is not None:
        sc.pp.filter_cells(adata, max_genes=max_genes)
    if max_counts is not None:
        sc.pp.filter_cells(adata, max_counts=max_counts)
    
    # Filter genes
    sc.pp.filter_genes(adata, min_cells=min_cells)
    
    # Filter by mitochondrial percentage
    if 'pct_counts_mt' in adata.obs.columns:
        adata = adata[adata.obs['pct_counts_mt'] < max_mito_pct, :]
    
    # Normalize data
    if normalize:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    
    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
    # Scale data
    if scale:
        sc.pp.scale(adata, max_value=10)
    
    # PCA
    sc.tl.pca(adata, n_comps=n_pcs, random_state=random_state)
    
    # Compute neighborhood graph
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, random_state=random_state)
    
    # UMAP
    sc.tl.umap(adata, random_state=random_state)
    
    return adata

def integrate_datasets(
    adata_list: List,
    batch_key: str = 'batch',
    method: str = 'harmony',
    **kwargs
):
    """
    Integrate multiple datasets to remove batch effects
    
    Parameters
    ----------
    adata_list : List[AnnData]
        List of AnnData objects to integrate
    batch_key : str, default='batch'
        Key for batch information in obs
    method : str, default='harmony'
        Integration method ('harmony', 'bbknn', 'scvi')
    **kwargs
        Additional arguments for integration method
        
    Returns
    -------
    AnnData
        Integrated AnnData object
    """
    
    if method == 'harmony':
        try:
            import harmonypy as hp
        except ImportError:
            raise ImportError("harmonypy is required for harmony integration")
        
        # Concatenate datasets
        adata = adata_list[0].concatenate(adata_list[1:], join='outer')
        
        # Run harmony
        sc.external.pp.harmony_integrate(adata, batch_key, **kwargs)
        
    elif method == 'bbknn':
        try:
            import bbknn
        except ImportError:
            raise ImportError("bbknn is required for BBKNN integration")
        
        # Concatenate datasets
        adata = adata_list[0].concatenate(adata_list[1:], join='outer')
        
        # Run BBKNN
        sc.external.pp.bbknn(adata, batch_key=batch_key, **kwargs)
        
    else:
        raise ValueError(f"Integration method {method} not supported")
    
    return adata 