# DGscRNA Tutorial

This tutorial will guide you through using the DGscRNA package for single-cell RNA-seq cell type annotation.

## Prerequisites

Before starting this tutorial, make sure you have:

1. Python 3.8 or higher installed
2. The DGscRNA package installed: `pip install dgscrna`
3. Basic knowledge of single-cell RNA-seq analysis
4. Familiarity with scanpy and pandas

## Installation

```bash
# Install from PyPI
pip install dgscrna

# Or install from source
git clone https://github.com/yourusername/DGscRNA.git
cd DGscRNA
pip install -e .
```

## Quick Start

### 1. Import Required Packages

```python
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import dgscrna as dg

# Set scanpy settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Set random seed for reproducibility
np.random.seed(42)
```

### 2. Load Your Data

```python
# Load your single-cell data
adata = sc.read_h5ad('your_data.h5ad')

print(f"Data shape: {adata.shape}")
print(f"Number of cells: {adata.n_obs}")
print(f"Number of genes: {adata.n_vars}")
```

### 3. Prepare Marker Sets

Create a folder with CSV files containing your marker genes. Each CSV file should have:
- Columns: Cell type names
- Rows: Marker genes for each cell type

Example marker file (`markers.csv`):
```csv
,T_Cells,B_Cells,Epithelial_Cells
0,CD3D,CD79A,EPCAM
1,CD3E,CD19,KRT18
2,CD3G,MS4A1,KRT19
```

### 4. Run the Complete Pipeline

```python
# Run the complete DGscRNA pipeline
results = dg.run_dgscrna_pipeline(
    adata=adata,
    marker_folder='path/to/marker/sets/',
    clustering_methods=['leiden', 'hdbscan'],
    deep_learning=True,
    epochs=10,
    batch_size=256,
    probability_threshold=0.9
)

# Get the annotated data
annotated_adata = results['adata']
```

### 5. Visualize Results

```python
# Plot UMAP with annotations
sc.pl.umap(annotated_adata, color=['leiden_clusters', 'CellMarker_Thyroid_mean_DGscRNA'])

# Get summary of results
summary = dg.summarize_results(results)
print(summary)
```

## Step-by-Step Workflow

### Step 1: Data Preprocessing

If your data is not already preprocessed, you can use the DGscRNA preprocessing function:

```python
# Preprocess the data
adata = dg.preprocess_adata(
    adata,
    min_genes=200,        # Minimum genes per cell
    max_genes=None,       # Maximum genes per cell (no limit)
    min_cells=3,          # Minimum cells expressing a gene
    max_counts=None,      # Maximum counts per cell (no limit)
    max_mito_pct=15.0,    # Maximum mitochondrial percentage
    normalize=True,       # Normalize data
    scale=True,           # Scale data
    n_pcs=30,            # Number of principal components
    n_neighbors=15,      # Number of neighbors for graph
    random_state=42
)
```

### Step 2: Clustering

Run multiple clustering algorithms:

```python
# Run clustering
adata = dg.run_clustering(
    adata,
    methods=['leiden', 'hdbscan', 'kmeans'],
    resolution=0.5,
    n_neighbors=15,
    random_state=42
)

# Check clustering results
for method in ['leiden', 'hdbscan', 'kmeans']:
    cluster_key = f"{method}_clusters"
    if cluster_key in adata.obs.columns:
        print(f"{method}: {adata.obs[cluster_key].nunique()} clusters")
```

### Step 3: Find Marker Genes

Find differentially expressed genes for each cluster:

```python
# Find markers for each clustering method
for method in ['leiden', 'hdbscan']:
    cluster_key = f"{method}_clusters"
    if cluster_key in adata.obs.columns:
        adata = dg.find_markers(
            adata,
            groupby=cluster_key,
            method='wilcoxon',
            n_genes=100
        )
```

### Step 4: Load and Score Marker Sets

```python
# Load marker sets
marker_sets = dg.load_marker_sets('path/to/marker/folder/')

# Score cell types for each clustering method and marker set
for method in ['leiden', 'hdbscan']:
    cluster_key = f"{method}_clusters"
    if cluster_key not in adata.obs.columns:
        continue
        
    for marker_set_name, markers_dict in marker_sets.items():
        adata = dg.score_cell_types(
            adata,
            markers_dict,
            cluster_key,
            cutoffs=['0.5', 'mean', 'none'],
            min_logfc=1.0
        )
```

### Step 5: Deep Learning Refinement

```python
# Train deep learning models and refine annotations
for method in ['leiden', 'hdbscan']:
    cluster_key = f"{method}_clusters"
    if cluster_key not in adata.obs.columns:
        continue
        
    for marker_set_name, markers_dict in marker_sets.items():
        # Find annotation keys
        annotation_keys = [col for col in adata.obs.columns 
                          if col.startswith(f"{marker_set_name}_{method}_")]
        
        for annotation_key in annotation_keys:
            # Train model
            model_results = dg.train_deep_model(
                adata,
                annotation_key,
                epochs=10,
                batch_size=256,
                random_state=42
            )
            
            # Predict refined annotations
            refined_annotations = dg.predict_cell_types(
                model_results,
                adata,
                annotation_key,
                probability_threshold=0.9
            )
            
            # Add to adata
            refined_key = f"{annotation_key}_DGscRNA"
            adata.obs[refined_key] = refined_annotations
```

## Advanced Usage

### Custom Deep Learning Model

You can customize the deep learning model architecture:

```python
from dgscrna.models.deep_model import DeepModel

# Create custom model
model = DeepModel(
    input_dim=2000,           # Number of genes
    num_classes=5,            # Number of cell types
    hidden_dims=[512, 256, 128]  # Custom architecture
)
```

### Batch Integration

If you have multiple datasets, you can integrate them:

```python
# Integrate multiple datasets
integrated_adata = dg.integrate_datasets(
    [adata1, adata2, adata3],
    batch_key='batch',
    method='harmony'
)
```

### Custom Marker Scoring

You can implement custom scoring methods:

```python
# Calculate marker enrichment
enrichment_scores = dg.calculate_marker_enrichment(
    adata,
    markers_dict,
    cluster_key='leiden_clusters',
    method='mean'
)
```

## Troubleshooting

### Common Issues

1. **Memory Issues**: Reduce batch size or use fewer genes
   ```python
   # Use fewer highly variable genes
   adata = dg.preprocess_adata(adata, n_pcs=20)
   
   # Reduce batch size
   model_results = dg.train_deep_model(adata, annotation_key, batch_size=128)
   ```

2. **No 'Undecided' Cells**: If all cells are already annotated, deep learning won't run
   ```python
   # Check annotation distribution
   print(adata.obs['your_annotation'].value_counts())
   ```

3. **Poor Clustering**: Adjust clustering parameters
   ```python
   # Try different resolution values
   adata = dg.run_clustering(adata, methods=['leiden'], resolution=1.0)
   ```

### Performance Tips

1. **Use GPU**: If available, the deep learning will automatically use CUDA
2. **Reduce Epochs**: For quick testing, use fewer epochs
3. **Filter Genes**: Use only highly variable genes for faster processing

## Output Interpretation

### Annotation Columns

The pipeline creates several types of annotation columns:

- `{method}_clusters`: Raw clustering results
- `{marker_set}_{method}_{cutoff}`: Marker-based annotations
- `{marker_set}_{method}_{cutoff}_DGscRNA`: Deep learning refined annotations

### Quality Metrics

Check the quality of annotations:

```python
# Get summary statistics
summary = dg.summarize_results(results)
print(summary)

# Look for high decision rates and reasonable number of cell types
good_annotations = summary[summary['Decision_Rate'] > 0.8]
print("High-quality annotations:")
print(good_annotations)
```

### Visualization

```python
# Plot multiple annotation methods
annotation_cols = ['leiden_clusters', 'CellMarker_Thyroid_mean_DGscRNA']
dg.plot_results(adata, annotation_cols, n_cols=2)

# Compare before and after refinement
sc.pl.umap(adata, color=['CellMarker_Thyroid_mean', 'CellMarker_Thyroid_mean_DGscRNA'])
```

## Best Practices

1. **Data Quality**: Ensure your data is well-preprocessed before running DGscRNA
2. **Marker Sets**: Use high-quality, tissue-specific marker sets
3. **Validation**: Always validate results with known markers or external data
4. **Parameters**: Start with default parameters and adjust based on your data
5. **Reproducibility**: Set random seeds for reproducible results

## Next Steps

After running DGscRNA, you might want to:

1. **Validate Results**: Compare with known cell type markers
2. **Downstream Analysis**: Use annotations for differential expression, trajectory analysis, etc.
3. **Custom Analysis**: Implement custom scoring methods for your specific use case
4. **Integration**: Combine with other cell type annotation tools

For more information, see the [API Documentation](api.md) and [Examples](../examples/). 