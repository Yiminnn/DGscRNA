#!/usr/bin/env python3
"""
Test script for the new differential expression-based density scoring.
"""

import numpy as np
import pandas as pd
import scanpy as sc
import scanpy_dgscrna as dgsc

# Set up scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

def test_deg_density_scoring():
    """Test the new density scoring with differential expression."""
    
    print("Creating synthetic test data...")
    
    # Create synthetic data
    np.random.seed(42)
    n_cells = 500
    n_genes = 2000
    
    # Create expression matrix
    X = np.random.negative_binomial(5, 0.3, size=(n_cells, n_genes)).astype(np.float32)
    
    # Create gene names
    gene_names = [f'GENE_{i}' for i in range(n_genes)]
    
    # Add some known marker genes
    marker_genes = ['CD3D', 'CD3E', 'CD8A', 'CD4', 'CD19', 'MS4A1', 'CD14', 'LYZ']
    gene_names[:len(marker_genes)] = marker_genes
    
    # Create AnnData object
    adata = sc.AnnData(X=X)
    adata.var_names = gene_names
    adata.obs_names = [f'Cell_{i}' for i in range(n_cells)]
    
    # Add synthetic clustering
    cluster_labels = np.random.choice(['0', '1', '2', '3'], size=n_cells)
    adata.obs['leiden'] = cluster_labels
    
    # Make some marker genes more highly expressed in specific clusters
    # Cluster 0: T cell markers
    t_markers = ['CD3D', 'CD3E', 'CD8A', 'CD4']
    for marker in t_markers:
        if marker in adata.var_names:
            marker_idx = list(adata.var_names).index(marker)
            cluster_0_cells = adata.obs['leiden'] == '0'
            adata.X[cluster_0_cells, marker_idx] += np.random.negative_binomial(10, 0.2, size=cluster_0_cells.sum())
    
    # Cluster 1: B cell markers
    b_markers = ['CD19', 'MS4A1']
    for marker in b_markers:
        if marker in adata.var_names:
            marker_idx = list(adata.var_names).index(marker)
            cluster_1_cells = adata.obs['leiden'] == '1'
            adata.X[cluster_1_cells, marker_idx] += np.random.negative_binomial(10, 0.2, size=cluster_1_cells.sum())
    
    # Normalize data
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    
    print(f"Created test data: {adata.n_obs} cells × {adata.n_vars} genes")
    
    # Test find_all_markers function
    print("Testing differential expression analysis...")
    deg_markers = dgsc.find_all_markers(adata, 'leiden')
    
    print(f"Found markers for {len(deg_markers)} clusters:")
    for cluster, markers_df in deg_markers.items():
        print(f"  Cluster {cluster}: {len(markers_df)} markers")
        if len(markers_df) > 0:
            top_markers = markers_df.head(3)['gene'].tolist()
            print(f"    Top markers: {top_markers}")
    
    # Create example marker sources
    marker_sources = {
        'Example_markers': {
            'T_cell': ['CD3D', 'CD3E', 'CD8A', 'CD4'],
            'B_cell': ['CD19', 'MS4A1'],
            'Monocyte': ['CD14', 'LYZ']
        }
    }
    
    print("\nTesting density scoring with DEG...")
    
    # Test density scoring
    dgsc.density_score(
        adata,
        marker_sources,
        cluster_columns=['leiden'],
        cutoffs=['0.5', 'mean', 'none'],
        use_deg=True
    )
    
    print("\nDensity scoring completed!")
    
    # Check results
    annotation_columns = [col for col in adata.obs.columns if 'leiden_Example_markers' in col]
    print(f"Added {len(annotation_columns)} annotation columns:")
    
    for col in annotation_columns:
        if not col.endswith('_score'):
            print(f"  {col}:")
            value_counts = adata.obs[col].value_counts()
            print(f"    {dict(value_counts)}")
    
    print("\nTest completed successfully!")
    return adata

if __name__ == "__main__":
    test_adata = test_deg_density_scoring()
