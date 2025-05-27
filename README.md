# scanpy-dgscrna

A scanpy extension for automated cell type annotation using density scoring and deep learning.

## Overview

`scanpy-dgscrna` integrates marker-based density scoring with deep learning to provide robust cell type annotation for single-cell RNA-seq data. The package combines:

- **Density scoring**: Uses marker gene sets from multiple sources with different thresholds
- **Deep learning**: Employs a neural network to classify cells based on expression patterns
- **Scanpy integration**: Seamlessly works with scanpy workflows and AnnData objects

## Installation

```bash
# Install from source
pip install -e .

# Or install with optional dependencies
pip install -e ".[all]"
```

## Quick Start

```python
import scanpy as sc
import scanpy_dgscrna as dgsc

# Load your data
adata = sc.read_h5ad('your_data.h5ad')

# Run the complete workflow
dgsc.run_dgscrna_workflow(
    adata, 
    marker_folder='path/to/markers',
    cluster_columns=['leiden', 'seurat_clusters']
)

# Check results
print(adata.obs.columns)  # See new annotation columns
```

## Detailed Usage

### 1. Prepare your data

```python
# Basic preprocessing
dgsc.pp.prepare_for_annotation(adata, min_genes=200, mt_threshold=20)

# Select highly variable genes
dgsc.pp.select_hvgs(adata, n_top_genes=2000)

# Add undecided labels for annotation
dgsc.pp.add_undecided_labels(adata, 'leiden', clusters_to_mark=[5, 8, 12])
```

### 2. Load marker gene sets

Organize your marker files in a folder structure:
```
markers/
├── CellMarker.xlsx
├── PanglaoDB.tsv
└── HumanPrimaryCellAtlas.csv
```

Each file should have:
- A gene column (default: 'gene')
- Cell type columns with marker genes

```python
# Load markers
marker_sources = dgsc.load_markers('path/to/markers')
print(f"Loaded {len(marker_sources)} marker sources")
```

### 3. Calculate density scores

```python
# Calculate density scores with different thresholds
dgsc.density_score(
    adata, 
    marker_sources, 
    cluster_columns=['leiden', 'seurat_clusters'],
    cutoffs=['0.5', 'mean', 'none']
)
```

### 4. Deep learning annotation

```python
# Annotate using deep learning
dgsc.dgscrna_annotate(
    adata, 
    annotation_column='leiden_CellMarker_Tcell_0.5',
    epochs=10,
    confidence_threshold=0.9
)
```

### 5. Visualization

```python
# Plot annotation results
dgsc.pl.plot_annotation_comparison(
    adata, 
    'leiden', 
    'leiden_CellMarker_Tcell_0.5_DGscRNA'
)

# Plot UMAP with annotations
dgsc.pl.plot_umap_annotations(
    adata, 
    ['leiden', 'leiden_CellMarker_Tcell_0.5_DGscRNA']
)

# Plot density scores
score_cols = [col for col in adata.obs.columns if 'CellMarker' in col and '0.5' in col]
dgsc.pl.plot_density_scores(adata, score_cols, 'leiden')
```

## API Reference

### Tools (`dgsc.tl`)

- `load_markers()`: Load marker gene sets from files
- `density_score()`: Calculate density scores for cell types
- `dgscrna_annotate()`: Deep learning-based annotation
- `run_dgscrna_workflow()`: Complete workflow

### Preprocessing (`dgsc.pp`)

- `prepare_for_annotation()`: Basic data preprocessing
- `select_hvgs()`: Select highly variable genes
- `add_undecided_labels()`: Add undecided labels for annotation

### Plotting (`dgsc.pl`)

- `plot_annotation_comparison()`: Compare annotations
- `plot_umap_annotations()`: UMAP plots with annotations
- `plot_density_scores()`: Visualize density scores
- `plot_marker_expression()`: Marker gene expression
- `plot_confusion_matrix()`: Confusion matrix for validation

## Input Data Format

### AnnData object
- `adata.X`: Gene expression matrix (cells × genes)
- `adata.obs`: Cell metadata including cluster assignments
- `adata.var`: Gene metadata

### Marker files
Supported formats: `.xlsx`, `.csv`, `.tsv`

Example structure:
```
gene        | Tcell      | Bcell      | Monocyte
------------|------------|------------|----------
CD3D        | 1          | 0          | 0
CD3E        | 1          | 0          | 0
CD19        | 0          | 1          | 0
CD14        | 0          | 0          | 1
```

## Parameters

### Density scoring thresholds
- `'none'`: Use raw expression scores
- `'mean'`: Use cluster mean as threshold
- `'0.5'`: Use 0.5 as threshold (for normalized data)
- Custom: Any float value

### Deep learning parameters
- `epochs`: Number of training epochs (default: 10)
- `batch_size`: Batch size for training (default: 256)
- `learning_rate`: Learning rate (default: 1e-3)
- `confidence_threshold`: Minimum confidence for predictions (default: 0.9)

## Examples

See the `examples/` directory for complete workflows:
- Basic annotation workflow
- Multi-source marker integration
- Custom threshold optimization
- Visualization gallery

## Citation

If you use this package in your research, please cite:

```
Your paper citation here
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please see CONTRIBUTING.md for guidelines.

## Support

For questions and support:
- GitHub Issues: https://github.com/yourusername/scanpy-dgscrna/issues
- Email: your.email@example.com