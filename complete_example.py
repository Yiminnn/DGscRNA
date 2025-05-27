"""
Complete example workflow using scanpy-dgscrna with real marker files.

This script demonstrates the full workflow from data preprocessing 
to cell type annotation using the package.
"""

import scanpy as sc
import pandas as pd
import numpy as np
import sys
import os

# Add the package to the path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import scanpy_dgscrna as dgsc

# Set scanpy settings
sc.settings.verbosity = 2
sc.settings.set_figure_params(dpi=80, facecolor='white')

def create_realistic_data():
    """Create more realistic single-cell data for demonstration."""
    print("Creating realistic example data...")
    
    # Parameters for realistic data
    n_cells = 2000
    n_genes = 2000
    
    np.random.seed(42)
    
    # Create cell type structure
    cell_types = ['T_cells'] * 500 + ['B_cells'] * 400 + ['Monocytes'] * 300 + \
                 ['NK_cells'] * 200 + ['Dendritic_cells'] * 150 + ['Unknown'] * 450
    
    # Create expression matrix with cell type specific patterns
    X = np.random.negative_binomial(3, 0.3, size=(n_cells, n_genes)).astype(float)
    
    # Add cell type specific expression patterns
    # T cells: higher expression in genes 0-199
    t_mask = np.array([ct == 'T_cells' for ct in cell_types])
    X[t_mask, :200] += np.random.poisson(5, size=(t_mask.sum(), 200))
    
    # B cells: higher expression in genes 200-399  
    b_mask = np.array([ct == 'B_cells' for ct in cell_types])
    X[b_mask, 200:400] += np.random.poisson(4, size=(b_mask.sum(), 200))
    
    # Monocytes: higher expression in genes 400-599
    m_mask = np.array([ct == 'Monocytes' for ct in cell_types])
    X[m_mask, 400:600] += np.random.poisson(6, size=(m_mask.sum(), 200))
    
    # NK cells: higher expression in genes 600-799
    nk_mask = np.array([ct == 'NK_cells' for ct in cell_types])
    X[nk_mask, 600:800] += np.random.poisson(3, size=(nk_mask.sum(), 200))
    
    # Dendritic cells: higher expression in genes 800-999
    dc_mask = np.array([ct == 'Dendritic_cells' for ct in cell_types])
    X[dc_mask, 800:1000] += np.random.poisson(4, size=(dc_mask.sum(), 200))
    
    # Create gene names (some matching common marker names)
    gene_names = []
    for i in range(n_genes):
        if i < 50:
            gene_names.append(f"CD3{chr(65 + i % 26)}")  # T cell markers
        elif i < 100:
            gene_names.append(f"CD19{chr(65 + i % 26)}")  # B cell markers  
        elif i < 150:
            gene_names.append(f"CD14{chr(65 + i % 26)}")  # Monocyte markers
        elif i < 200:
            gene_names.append(f"CD56{chr(65 + i % 26)}")  # NK markers
        else:
            gene_names.append(f"Gene_{i}")
    
    # Create cell barcodes
    cell_barcodes = [f"Cell_{i}" for i in range(n_cells)]
    
    # Create AnnData object
    adata = sc.AnnData(X)
    adata.var_names = gene_names
    adata.obs_names = cell_barcodes
    
    # Add metadata
    adata.obs['true_celltype'] = cell_types
    adata.obs['total_counts'] = np.array(X.sum(axis=1)).flatten()
    adata.obs['n_genes'] = (X > 0).sum(axis=1)
    
    # Add mitochondrial genes
    mt_genes = np.random.choice(range(n_genes), 50, replace=False)
    adata.var['mt'] = False
    adata.var.iloc[mt_genes, adata.var.columns.get_loc('mt')] = True
    
    print(f"Created dataset: {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"Cell type distribution: {pd.Series(cell_types).value_counts().to_dict()}")
    
    return adata

def main():
    """Run the complete workflow."""
    print("=== Complete DGscRNA Workflow Example ===\n")
    
    # Step 1: Create or load data
    adata = create_realistic_data()
    
    # Step 2: Preprocessing
    print("\n=== Step 2: Preprocessing ===")
    
    # Basic preprocessing
    dgsc.pp.prepare_for_annotation(adata, min_genes=100, mt_threshold=25)
    print(f"After QC: {adata.n_obs} cells, {adata.n_vars} genes")
    
    # Highly variable genes
    dgsc.pp.select_hvgs(adata, n_top_genes=1500)
    
    # Basic analysis for clustering
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=15)
    sc.tl.leiden(adata, resolution=0.5)
    sc.tl.umap(adata)
    
    print(f"Found {len(adata.obs['leiden'].unique())} leiden clusters")
    
    # Step 3: Load marker gene sets from actual files
    print("\n=== Step 3: Loading Marker Gene Sets ===")
    
    try:
        marker_sources = dgsc.load_markers('./data')
        print(f"Successfully loaded {len(marker_sources)} marker sources:")
        for source_name, markers_df in marker_sources.items():
            n_celltypes = len(markers_df.columns)
            n_markers = markers_df.notna().sum().sum()
            print(f"  - {source_name}: {n_celltypes} cell types, {n_markers} total markers")
    except Exception as e:
        print(f"Could not load markers from files: {e}")
        print("Creating example markers instead...")
        marker_sources = create_example_markers()
    
    # Step 4: Density scoring
    print("\n=== Step 4: Density Scoring ===")
    
    dgsc.density_score(
        adata,
        marker_sources,
        cluster_columns=['leiden'],
        cutoffs=['0.5', 'mean', 'none']
    )
    
    # Count density score columns
    score_cols = [col for col in adata.obs.columns 
                 if any(source in col for source in marker_sources.keys())]
    print(f"Created {len(score_cols)} density score columns")
    
    # Show some example scores
    if score_cols:
        example_scores = adata.obs[score_cols[:3]].describe()
        print("Example density scores:")
        print(example_scores)
    
    # Step 5: Create annotation columns for deep learning
    print("\n=== Step 5: Preparing for Deep Learning ===")
    
    # Find the best scoring combinations for each cell type
    celltype_cols = []
    for source_name, markers_df in marker_sources.items():
        for celltype in markers_df.columns:
            for cutoff in ['0.5', 'mean']:
                col_name = f"leiden_{source_name}_{celltype}_{cutoff}"
                if col_name in adata.obs.columns:
                    celltype_cols.append(col_name)
    
    print(f"Found {len(celltype_cols)} potential annotation columns")
    
    # For each promising combination, create annotation labels
    annotation_columns = []
    for col in celltype_cols[:5]:  # Limit to first 5 for demo
        # Create binary annotation based on top scoring cells
        scores = adata.obs[col]
        if scores.max() > 0:  # Only if there are non-zero scores
            threshold = scores.quantile(0.8)  # Top 20% of cells
            
            # Create annotation column
            annotation_col = f"{col}_annotation"
            adata.obs[annotation_col] = 'Undecided'
            
            # Mark high-scoring cells as positive for this cell type
            celltype_name = col.split('_')[2]  # Extract cell type name
            high_score_mask = scores >= threshold
            adata.obs.loc[high_score_mask, annotation_col] = celltype_name
            
            annotation_columns.append(annotation_col)
            
            positive_cells = high_score_mask.sum()
            print(f"  {annotation_col}: {positive_cells} {celltype_name}, {len(scores) - positive_cells} Undecided")
    
    # Step 6: Deep learning annotation
    print("\n=== Step 6: Deep Learning Annotation ===")
    
    for annotation_col in annotation_columns[:3]:  # Limit for demo
        print(f"\nProcessing {annotation_col}...")
        
        try:
            dgsc.dgscrna_annotate(
                adata,
                annotation_column=annotation_col,
                epochs=5,  # Reduced for demo
                confidence_threshold=0.7,
                batch_size=64
            )
            
            # Check results
            result_col = f"{annotation_col}_DGscRNA"
            if result_col in adata.obs.columns:
                result_counts = adata.obs[result_col].value_counts()
                print(f"  Results: {result_counts.to_dict()}")
        
        except Exception as e:
            print(f"  Error in deep learning annotation: {e}")
            continue
    
    # Step 7: Evaluation and visualization
    print("\n=== Step 7: Results Summary ===")
    
    # Show all DGscRNA results
    dgscrna_cols = [col for col in adata.obs.columns if 'DGscRNA' in col]
    if dgscrna_cols:
        print(f"\nCreated {len(dgscrna_cols)} DGscRNA annotation columns:")
        for col in dgscrna_cols:
            counts = adata.obs[col].value_counts()
            print(f"  {col}: {dict(counts)}")
    
    # Compare with true cell types (if available)
    if 'true_celltype' in adata.obs.columns:
        print(f"\nTrue cell type distribution:")
        true_counts = adata.obs['true_celltype'].value_counts()
        print(f"  {dict(true_counts)}")
    
    # Save results
    print("\n=== Saving Results ===")
    adata.write('annotated_data.h5ad')
    print("Saved annotated data to 'annotated_data.h5ad'")
    
    # Export annotation summary
    annotation_summary = adata.obs[['leiden', 'true_celltype'] + dgscrna_cols].copy()
    annotation_summary.to_csv('annotation_results.csv')
    print("Saved annotation summary to 'annotation_results.csv'")
    
    print("\n✓ Complete workflow finished successfully!")
    return adata

def create_example_markers():
    """Create example marker gene sets as fallback."""
    marker_sources = {}
    
    # Example markers that match our synthetic data
    markers1 = pd.DataFrame({
        'T_cells': ['CD3A', 'CD3B', 'CD3C', 'CD3D', None],
        'B_cells': ['CD19A', 'CD19B', 'CD19C', None, None],
        'Monocytes': ['CD14A', 'CD14B', 'CD14C', 'CD14D', None],
        'NK_cells': ['CD56A', 'CD56B', None, None, None]
    })
    marker_sources['Example_Markers'] = markers1
    
    return marker_sources

if __name__ == "__main__":
    try:
        adata = main()
        print("\nWorkflow completed successfully!")
    except Exception as e:
        print(f"\nWorkflow failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
