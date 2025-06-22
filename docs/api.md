# DGscRNA API Documentation

## Core Modules

### Preprocessing Module (`dgscrna.core.preprocessing`)

#### `preprocess_adata(adata, min_genes=200, max_genes=None, min_cells=3, max_counts=None, max_mito_pct=15.0, normalize=True, scale=True, n_pcs=30, n_neighbors=15, random_state=42)`

Preprocess AnnData object with quality control and normalization.

**Parameters:**
- `adata` (AnnData): Annotated data matrix
- `min_genes` (int, default=200): Minimum number of genes expressed per cell
- `max_genes` (int, optional): Maximum number of genes expressed per cell
- `min_cells` (int, default=3): Minimum number of cells expressing a gene
- `max_counts` (int, optional): Maximum number of counts per cell
- `max_mito_pct` (float, default=15.0): Maximum percentage of mitochondrial genes
- `normalize` (bool, default=True): Whether to normalize data
- `scale` (bool, default=True): Whether to scale data
- `n_pcs` (int, default=30): Number of principal components
- `n_neighbors` (int, default=15): Number of neighbors for neighborhood graph
- `random_state` (int, default=42): Random state for reproducibility

**Returns:**
- `AnnData`: Preprocessed AnnData object

#### `integrate_datasets(adata_list, batch_key='batch', method='harmony', **kwargs)`

Integrate multiple datasets to remove batch effects.

**Parameters:**
- `adata_list` (List[AnnData]): List of AnnData objects to integrate
- `batch_key` (str, default='batch'): Key for batch information in obs
- `method` (str, default='harmony'): Integration method ('harmony', 'bbknn', 'scvi')
- `**kwargs`: Additional arguments for integration method

**Returns:**
- `AnnData`: Integrated AnnData object

### Clustering Module (`dgscrna.core.clustering`)

#### `run_clustering(adata, methods=['leiden', 'hdbscan', 'kmeans'], resolution=0.5, n_neighbors=15, n_clusters=None, random_state=42, **kwargs)`

Run multiple clustering algorithms on the data.

**Parameters:**
- `adata` (AnnData): Annotated data matrix
- `methods` (List[str], default=['leiden', 'hdbscan', 'kmeans']): List of clustering methods to run
- `resolution` (float, default=0.5): Resolution parameter for Leiden clustering
- `n_neighbors` (int, default=15): Number of neighbors for neighborhood graph
- `n_clusters` (int, optional): Number of clusters for K-means
- `random_state` (int, default=42): Random state for reproducibility
- `**kwargs`: Additional arguments for clustering methods

**Returns:**
- `AnnData`: AnnData object with clustering results added to obs

#### `find_markers(adata, groupby, method='wilcoxon', key_added='rank_genes_groups', n_genes=100, **kwargs)`

Find differentially expressed genes for each cluster.

**Parameters:**
- `adata` (AnnData): Annotated data matrix
- `groupby` (str): Key in obs for grouping
- `method` (str, default='wilcoxon'): Method for differential expression testing
- `key_added` (str, default='rank_genes_groups'): Key to store results in uns
- `n_genes` (int, default=100): Number of top genes to return per cluster
- `**kwargs`: Additional arguments for sc.tl.rank_genes_groups

**Returns:**
- `AnnData`: AnnData object with marker genes results

### Marker Scoring Module (`dgscrna.core.marker_scoring`)

#### `load_marker_sets(marker_folder, file_format='csv')`

Load marker sets from CSV files in a folder.

**Parameters:**
- `marker_folder` (str): Path to folder containing marker CSV files
- `file_format` (str, default='csv'): File format of marker files

**Returns:**
- `Dict[str, Dict[str, List[str]]]`: Dictionary with marker set names as keys and cell type markers as values

#### `score_cell_types(adata, markers_dict, cluster_key, deg_key='rank_genes_groups', cutoffs=['0.5', 'mean', 'none'], min_logfc=1.0, random_state=42)`

Score cell types based on marker gene expression.

**Parameters:**
- `adata` (AnnData): Annotated data matrix
- `markers_dict` (Dict[str, List[str]]): Dictionary with cell type names as keys and marker gene lists as values
- `cluster_key` (str): Key in obs for cluster assignments
- `deg_key` (str, default='rank_genes_groups'): Key in uns containing differential expression results
- `cutoffs` (List[str], default=['0.5', 'mean', 'none']): List of cutoff methods for scoring
- `min_logfc` (float, default=1.0): Minimum log fold change for marker genes
- `random_state` (int, default=42): Random state for reproducibility

**Returns:**
- `AnnData`: AnnData object with cell type scores added to obs

### Deep Learning Module (`dgscrna.core.deep_learning`)

#### `train_deep_model(adata, annotation_key, epochs=10, batch_size=256, lr=1e-3, use_highly_variable=True, hidden_dims=None, device=None, random_state=42)`

Train deep learning model to refine cell type annotations.

**Parameters:**
- `adata` (AnnData): Annotated data matrix
- `annotation_key` (str): Key in obs containing cell type annotations
- `epochs` (int, default=10): Number of training epochs
- `batch_size` (int, default=256): Batch size for training
- `lr` (float, default=1e-3): Learning rate
- `use_highly_variable` (bool, default=True): Whether to use only highly variable genes
- `hidden_dims` (List[int], optional): Hidden layer dimensions
- `device` (str, optional): Device to use ('cpu' or 'cuda')
- `random_state` (int, default=42): Random state for reproducibility

**Returns:**
- `Dict`: Dictionary containing training results and model

#### `predict_cell_types(model_results, adata, annotation_key, probability_threshold=0.9, use_highly_variable=True)`

Predict cell types using trained model.

**Parameters:**
- `model_results` (Dict): Results from train_deep_model
- `adata` (AnnData): Annotated data matrix
- `annotation_key` (str): Key in obs containing original annotations
- `probability_threshold` (float, default=0.9): Probability threshold for confident predictions
- `use_highly_variable` (bool, default=True): Whether to use only highly variable genes

**Returns:**
- `pd.Series`: Series with refined cell type annotations

### Utility Module (`dgscrna.core.utils`)

#### `run_dgscrna_pipeline(adata, marker_folder, clustering_methods=['leiden', 'hdbscan'], deep_learning=True, epochs=10, batch_size=256, probability_threshold=0.9, random_state=42, **kwargs)`

Complete DGscRNA pipeline for cell type annotation.

**Parameters:**
- `adata` (AnnData): Annotated data matrix (should be preprocessed and normalized)
- `marker_folder` (str): Path to folder containing marker CSV files
- `clustering_methods` (List[str], default=['leiden', 'hdbscan']): List of clustering methods to use
- `deep_learning` (bool, default=True): Whether to use deep learning refinement
- `epochs` (int, default=10): Number of training epochs for deep learning
- `batch_size` (int, default=256): Batch size for deep learning
- `probability_threshold` (float, default=0.9): Probability threshold for confident predictions
- `random_state` (int, default=42): Random state for reproducibility
- `**kwargs`: Additional arguments for preprocessing and clustering

**Returns:**
- `Dict`: Dictionary containing results and trained models

#### `summarize_results(pipeline_results)`

Summarize pipeline results.

**Parameters:**
- `pipeline_results` (Dict): Results from run_dgscrna_pipeline

**Returns:**
- `pd.DataFrame`: Summary of results

#### `plot_results(adata, annotation_cols, n_cols=2, figsize=(15, 10))`

Plot annotation results.

**Parameters:**
- `adata` (AnnData): Annotated data matrix
- `annotation_cols` (List[str]): List of annotation columns to plot
- `n_cols` (int, default=2): Number of columns in the plot grid
- `figsize` (tuple, default=(15, 10)): Figure size

**Returns:**
- None (displays plot)

## Models Module

### Deep Model (`dgscrna.models.deep_model`)

#### `DeepModel(input_dim, num_classes, hidden_dims=None)`

Deep neural network for cell type annotation refinement.

**Parameters:**
- `input_dim` (int): Input dimension (number of genes)
- `num_classes` (int): Number of cell type classes
- `hidden_dims` (list, optional): List of hidden layer dimensions. Default: [256, 128]

**Methods:**
- `forward(x)`: Forward pass
- `predict(x)`: Make predictions
- `predict_proba(x)`: Get prediction probabilities

## Data Format Specifications

### Input Data
- **AnnData object**: Preprocessed and normalized single-cell data
- **Marker folder**: CSV files where columns are cell type names and rows are marker genes

### Marker File Format
```csv
,CellType1,CellType2,CellType3
0,Gene1,Gene4,Gene7
1,Gene2,Gene5,Gene8
2,Gene3,Gene6,Gene9
```

### Output Data
- **AnnData object**: With added annotation columns
- **Results dictionary**: Training scores and metrics 