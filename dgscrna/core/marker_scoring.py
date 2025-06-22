"""
Marker scoring module for DGscRNA package
"""

import scanpy as sc
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Union
import os
import glob
import warnings
warnings.filterwarnings('ignore')

def load_marker_sets(marker_folder: str, file_format: str = 'csv') -> Dict[str, Dict[str, List[str]]]:
    """
    Load marker sets from CSV files in a folder
    
    Parameters
    ----------
    marker_folder : str
        Path to folder containing marker CSV files
    file_format : str, default='csv'
        File format of marker files
        
    Returns
    -------
    Dict[str, Dict[str, List[str]]]
        Dictionary with marker set names as keys and cell type markers as values
    """
    
    marker_sets = {}
    
    # Find all CSV files in the folder
    pattern = os.path.join(marker_folder, f"*.{file_format}")
    files = glob.glob(pattern)
    
    for file_path in files:
        # Get marker set name from filename
        marker_set_name = os.path.splitext(os.path.basename(file_path))[0]
        
        # Read CSV file
        df = pd.read_csv(file_path, index_col=0)
        
        # Convert to dictionary format
        cell_type_markers = {}
        for col in df.columns:
            # Get non-null marker genes
            markers = df[col].dropna().tolist()
            # Filter out empty strings
            markers = [marker for marker in markers if marker and marker.strip()]
            if markers:  # Only add if there are markers
                cell_type_markers[col] = markers
        
        marker_sets[marker_set_name] = cell_type_markers
    
    return marker_sets

def score_cell_types(
    adata,
    markers_dict: Dict[str, List[str]],
    cluster_key: str,
    marker_set_name: str = None,
    deg_key: str = 'rank_genes_groups',
    cutoffs: List[str] = ['0.5', 'mean', 'none'],
    min_logfc: float = 1.0,
    random_state: int = 42
):
    """
    Score cell types based on marker gene expression
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    markers_dict : Dict[str, List[str]]
        Dictionary with cell type names as keys and marker gene lists as values
    cluster_key : str
        Key in obs for cluster assignments
    marker_set_name : str, optional
        Name of the marker set (used for creating annotation column names)
    deg_key : str, default='rank_genes_groups'
        Key in uns containing differential expression results
    cutoffs : List[str], default=['0.5', 'mean', 'none']
        List of cutoff methods for scoring
    min_logfc : float, default=1.0
        Minimum log fold change for marker genes
    random_state : int, default=42
        Random state for reproducibility
        
    Returns
    -------
    AnnData
        AnnData object with cell type scores added to obs
    """
    
    # Check if cluster key exists
    if cluster_key not in adata.obs.columns:
        raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")
    
    # Check if DEG results exist
    if deg_key not in adata.uns:
        raise ValueError(f"DEG key '{deg_key}' not found in adata.uns. Run find_markers first.")
    
    deg_results = adata.uns[deg_key]
    cluster_names = deg_results['names'].dtype.names
    
    # Get expression matrix (use scaled data if available)
    if 'X_scaled' in adata.obsm:
        expr_matrix = adata.obsm['X_scaled']
    else:
        expr_matrix = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    
    # Get gene names
    gene_names = adata.var_names.tolist()
    
    for cutoff in cutoffs:
        print(f"Scoring with cutoff: {cutoff}")
        
        # Create scoring matrix
        n_cell_types = len(markers_dict)
        n_clusters = len(cluster_names)
        scoring_matrix = np.zeros((n_cell_types, n_clusters))
        
        cell_type_names = list(markers_dict.keys())
        
        # Score each cluster
        for i, cluster in enumerate(cluster_names):
            # Get DEGs for this cluster
            cluster_genes = deg_results['names'][cluster]
            cluster_logfcs = deg_results['logfoldchanges'][cluster]
            cluster_pvals = deg_results['pvals_adj'][cluster]
            
            # Filter by significance and log fold change
            significant = (cluster_pvals < 0.05) & (cluster_logfcs > min_logfc)
            significant_genes = cluster_genes[significant]
            significant_logfcs = cluster_logfcs[significant]
            
            # Score each cell type
            for j, (cell_type, marker_genes) in enumerate(markers_dict.items()):
                # Find intersection with significant genes
                intersection = set(marker_genes) & set(significant_genes)
                
                if intersection:
                    # Calculate score based on log fold changes
                    intersection_indices = [np.where(significant_genes == gene)[0][0] 
                                          for gene in intersection if gene in significant_genes]
                    intersection_logfcs = significant_logfcs[intersection_indices]
                    
                    # Average score normalized by number of markers
                    score = np.sum(intersection_logfcs) / len(marker_genes)
                    
                    # Penalize small marker sets
                    if len(marker_genes) <= 1:
                        score *= 0.8
                    
                    scoring_matrix[j, i] = score
        
        # Assign cell types based on maximum scores
        max_scores = np.max(scoring_matrix, axis=0)
        best_indices = np.argmax(scoring_matrix, axis=0)
        best_cell_types = [cell_type_names[i] for i in best_indices]
        
        # Apply cutoffs
        if cutoff == '0.5':
            threshold = 0.5
        elif cutoff == 'mean':
            threshold = np.mean(max_scores)
        else:  # 'none'
            threshold = 0
        
        # Create cluster-to-cell-type mapping
        cluster_assignments = {}
        for i, cluster in enumerate(cluster_names):
            cell_type = best_cell_types[i]
            score = max_scores[i]
            if score >= threshold:
                cluster_assignments[cluster] = cell_type
            else:
                cluster_assignments[cluster] = 'Undecided'
        
        # Map cluster assignments to individual cells
        final_assignments = []
        for cell_cluster in adata.obs[cluster_key]:
            final_assignments.append(cluster_assignments.get(cell_cluster, 'Unknown'))
        
        # Create annotation key based on marker set name and clustering method
        if marker_set_name:
            # Extract method name from cluster_key (e.g., 'hdbscan_clusters' -> 'hdbscan')
            method_name = cluster_key.replace('_clusters', '') if cluster_key.endswith('_clusters') else cluster_key
            annotation_key = f"{marker_set_name}_{method_name}_{cutoff}"
        else:
            annotation_key = f"{cluster_key}_{cutoff}"
        
        # Add to obs
        adata.obs[annotation_key] = final_assignments
    
    return adata

def calculate_marker_enrichment(
    adata,
    markers_dict: Dict[str, List[str]],
    cluster_key: str,
    method: str = 'mean',
    use_raw: bool = False
):
    """
    Calculate marker gene enrichment scores for each cluster
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    markers_dict : Dict[str, List[str]]
        Dictionary with cell type names as keys and marker gene lists as values
    cluster_key : str
        Key in obs for cluster assignments
    method : str, default='mean'
        Method for calculating enrichment ('mean', 'median', 'sum')
    use_raw : bool, default=False
        Whether to use raw counts instead of normalized data
        
    Returns
    -------
    pd.DataFrame
        DataFrame with enrichment scores for each cell type and cluster
    """
    
    # Get expression matrix
    if use_raw and adata.raw is not None:
        expr_matrix = adata.raw.X.toarray() if hasattr(adata.raw.X, 'toarray') else adata.raw.X
        gene_names = adata.raw.var_names.tolist()
    else:
        expr_matrix = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
        gene_names = adata.var_names.tolist()
    
    # Get cluster assignments
    clusters = adata.obs[cluster_key].unique()
    
    # Calculate enrichment scores
    enrichment_scores = {}
    
    for cell_type, marker_genes in markers_dict.items():
        # Find marker genes in the dataset
        available_markers = [gene for gene in marker_genes if gene in gene_names]
        
        if not available_markers:
            enrichment_scores[cell_type] = {cluster: 0 for cluster in clusters}
            continue
        
        # Get indices of marker genes
        marker_indices = [gene_names.index(gene) for gene in available_markers]
        
        # Calculate scores for each cluster
        cluster_scores = {}
        for cluster in clusters:
            # Get cells in this cluster
            cluster_mask = adata.obs[cluster_key] == cluster
            cluster_expr = expr_matrix[cluster_mask, :]
            
            # Calculate marker expression
            marker_expr = cluster_expr[:, marker_indices]
            
            if method == 'mean':
                score = np.mean(marker_expr)
            elif method == 'median':
                score = np.median(marker_expr)
            elif method == 'sum':
                score = np.sum(marker_expr)
            else:
                raise ValueError(f"Method '{method}' not supported")
            
            cluster_scores[cluster] = score
        
        enrichment_scores[cell_type] = cluster_scores
    
    return pd.DataFrame(enrichment_scores) 