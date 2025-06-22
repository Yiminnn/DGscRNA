# DGscRNA

A Python package for single-cell RNA-seq cell type annotation using marker-based scoring and deep learning refinement.

## Overview

DGscRNA combines traditional marker-based cell type scoring with deep learning to resolve ambiguous cell type assignments in single-cell RNA-seq data. The workflow includes:

1. **Preprocessing**: Quality control, normalization, and dimensionality reduction
2. **Clustering**: Multiple clustering algorithms (Leiden, HDBSCAN, K-means)
3. **Marker Scoring**: Density-based scoring using known cell type markers
4. **Deep Learning**: Neural network refinement of ambiguous annotations

## Installation

```bash
pip install dgscrna
```

Or install from source:

```bash
git clone https://github.com/yourusername/DGscRNA.git
cd DGscRNA
pip install -e .
```

## Quick Start

```python
import scanpy as sc
import dgscrna as dg

# Load your data
adata = sc.read_h5ad('your_data.h5ad')

# Run the complete pipeline
results = dg.run_dgscrna_pipeline(
    adata=adata,
    marker_folder='path/to/marker/sets/',
    clustering_methods=['leiden', 'hdbscan'],
    deep_learning=True
)

# View results
sc.pl.umap(adata, color=['leiden', 'CellMarker_Thyroid_mean_DGscRNA'])
```

## Input Data Format

### Single-cell Data
- **Format**: AnnData object (scanpy/anndata)
- **Requirements**: Preprocessed and normalized gene expression matrix

### Marker Sets
- **Format**: CSV files in a folder
- **Structure**: Columns are cell type names, rows are marker genes
- **Example**:
```csv
,CellType1,CellType2,CellType3
0,Gene1,Gene4,Gene7
1,Gene2,Gene5,Gene8
2,Gene3,Gene6,Gene9
```

## Output

- **AnnData object**: With added annotation columns
- **Results dictionary**: Training scores and metrics
- **Visualization**: UMAP plots with annotations

## Documentation

- [API Reference](docs/api.md)
- [Installation Guide](docs/installation.md)
- [Tutorial](docs/tutorial.md)
- [Examples](examples/)

## Citation

If you use DGscRNA in your research, please cite:

```bibtex
@article{dgscrna2023,
  title={DGscRNA: Deep learning-guided single-cell RNA-seq cell type annotation},
  author={Your Name},
  journal={Bioinformatics},
  year={2023}
}
```

## License

MIT License - see LICENSE file for details.

## Contributing

Contributions are welcome! Please read our contributing guidelines and submit pull requests.

## Support

For questions and support, please open an issue on GitHub or contact the maintainers. 