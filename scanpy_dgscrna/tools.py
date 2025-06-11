"""
Tools module for scanpy-dgscrna package.

Contains main functionality for marker loading, density scoring, and deep learning annotation.
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.utils.data as data_utils
import torchmetrics
from torch.autograd import Variable
from typing import Union, List, Dict, Optional, Tuple
import warnings
from pathlib import Path

# Suppress numba warnings
from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
warnings.simplefilter('ignore', category=NumbaDeprecationWarning)
warnings.simplefilter('ignore', category=NumbaPendingDeprecationWarning)


def load_markers(
    marker_folder: Union[str, Path],
    file_pattern: str = "*.xlsx",
    gene_column: str = "gene",
    celltype_columns: Optional[List[str]] = None
) -> Dict[str, pd.DataFrame]:
    """
    Load marker gene sets from files in a folder.
    
    Parameters
    ----------
    marker_folder : str or Path
        Path to folder containing marker files
    file_pattern : str, default "*.xlsx"
        Pattern to match marker files
    gene_column : str, default "gene"
        Column name containing gene names
    celltype_columns : list of str, optional
        Specific celltype columns to load. If None, loads all columns except gene_column
        
    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary with source names as keys and marker DataFrames as values
    """
    marker_folder = Path(marker_folder)
    marker_sources = {}
    
    # Find all marker files
    if file_pattern.endswith('.xlsx'):
        marker_files = list(marker_folder.glob("*.xlsx"))
    elif file_pattern.endswith('.tsv'):
        marker_files = list(marker_folder.glob("*.tsv"))
    elif file_pattern.endswith('.csv'):
        marker_files = list(marker_folder.glob("*.csv"))
    else:
        marker_files = list(marker_folder.glob(file_pattern))
    
    for marker_file in marker_files:
        source_name = marker_file.stem
        
        # Load file based on extension
        if marker_file.suffix == '.xlsx':
            df = pd.read_excel(marker_file)
        elif marker_file.suffix == '.tsv':
            df = pd.read_csv(marker_file, sep='\t')
        elif marker_file.suffix == '.csv':
            df = pd.read_csv(marker_file)
        else:
            print(f"Unsupported file format: {marker_file}")
            continue
            
        # Process marker data
        if celltype_columns is None:
            celltype_columns = [col for col in df.columns if col != gene_column]
        
        # Create clean marker dictionary for this source
        marker_dict = {}
        for celltype in celltype_columns:
            if celltype in df.columns:
                markers = df[df[celltype].notna()][gene_column].tolist()
                if markers:
                    marker_dict[celltype] = markers
        
        if marker_dict:
            marker_sources[source_name] = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in marker_dict.items()]))
            
    return marker_sources


def density_score(
    adata: sc.AnnData,
    marker_sources: Dict[str, pd.DataFrame],
    cluster_columns: List[str],
    cutoffs: List[str] = ['0.5', 'mean', 'none'],
    deg_markers: Optional[Dict[str, Dict[str, pd.DataFrame]]] = None,
    use_deg: bool = True,
    min_logfc_threshold: float = 1.0,
    copy: bool = False
) -> Optional[sc.AnnData]:
    """
    Calculate density scores for cell type annotation using differential expression analysis.
    
    This function implements the density scoring strategy similar to the R version,
    using differential expression markers to calculate cluster-specific scores.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    marker_sources : dict
        Dictionary of marker gene sets from different sources
    cluster_columns : list of str
        Column names in adata.obs containing cluster assignments
    cutoffs : list of str, default ['0.5', 'mean', 'none']
        Thresholds for density scoring
    deg_markers : dict, optional
        Pre-computed differential expression markers. If None, will compute them.
    use_deg : bool, default True
        Whether to use differential expression analysis or simple mean expression
    min_logfc_threshold : float, default 1.0
        Minimum log fold change threshold for considering genes as upregulated
    copy : bool, default False
        Return a copy instead of writing to adata
        
    Returns
    -------
    AnnData or None
        Returns AnnData if copy=True, otherwise modifies adata in place
    """
    if copy:
        adata = adata.copy()
    
    print("Starting density scoring with differential expression analysis...")
    
    # Compute differential expression markers if not provided
    if use_deg and deg_markers is None:
        print("Computing differential expression markers...")
        deg_markers = {}
        for clustering in cluster_columns:
            if clustering in adata.obs.columns:
                deg_markers[clustering] = find_all_markers(adata, clustering)
    
    # Process each clustering method
    for clustering in cluster_columns:
        if clustering not in adata.obs.columns:
            print(f"Clustering '{clustering}' not found in adata.obs, skipping...")
            continue
        
        print(f"Processing clustering: {clustering}")
        clusters = sorted(adata.obs[clustering].unique())
        
        # Get DEG markers for this clustering if using DEG approach
        if use_deg and deg_markers and clustering in deg_markers:
            cluster_deg_markers = deg_markers[clustering]
        else:
            cluster_deg_markers = None
        
        # Process each marker source
        for source_name, marker_dict in marker_sources.items():
            print(f"  Processing marker source: {source_name}")
            
            # Convert marker dict to list of cell types and their markers
            if isinstance(marker_dict, dict):
                cell_types = list(marker_dict.keys())
                
                # Create satisfaction matrix: rows = cell types, cols = clusters
                satisfaction_matrix = np.zeros((len(cell_types), len(clusters)))
                
                for cluster_idx, cluster in enumerate(clusters):
                    cluster_str = str(cluster)
                    
                    # Get upregulated genes for this cluster if using DEG
                    if use_deg and cluster_deg_markers and cluster_str in cluster_deg_markers:
                        cluster_markers_df = cluster_deg_markers[cluster_str]
                        upregulated_genes = cluster_markers_df[
                            cluster_markers_df['avg_log2FC'] > min_logfc_threshold
                        ]['gene'].tolist()
                    else:
                        # Fallback: use all genes (for non-DEG approach)
                        upregulated_genes = list(adata.var_names)
                    
                    # Calculate satisfaction scores for each cell type
                    for ct_idx, cell_type in enumerate(cell_types):
                        cell_type_markers = marker_dict[cell_type]
                        if not isinstance(cell_type_markers, list):
                            continue
                        
                        # Find intersection between cell type markers and upregulated genes
                        intersecting_genes = list(set(cell_type_markers) & set(upregulated_genes))
                        
                        if len(intersecting_genes) > 0 and use_deg and cluster_deg_markers and cluster_str in cluster_deg_markers:
                            # Sum log fold changes of intersecting genes
                            cluster_markers_df = cluster_deg_markers[cluster_str]
                            intersecting_logfc = cluster_markers_df[
                                cluster_markers_df['gene'].isin(intersecting_genes)
                            ]['avg_log2FC'].sum()
                            
                            # Normalize by the number of markers in the cell type
                            satisfaction_score = intersecting_logfc / len(cell_type_markers)
                            
                            # Penalize small marker sets (similar to R version)
                            if len(cell_type_markers) <= 1:
                                satisfaction_score *= 0.8
                                
                        elif len(intersecting_genes) > 0:
                            # Fallback: use simple proportion if DEG not available
                            satisfaction_score = len(intersecting_genes) / len(cell_type_markers)
                        else:
                            satisfaction_score = 0.0
                        
                        satisfaction_matrix[ct_idx, cluster_idx] = satisfaction_score
                
                # Find maximum scores and assign cell types
                max_scores = np.max(satisfaction_matrix, axis=0)
                
                # Process each cutoff strategy
                for cutoff in cutoffs:
                    cluster_annotations = []
                    cluster_scores = []
                    
                    for cluster_idx, cluster in enumerate(clusters):
                        max_score = max_scores[cluster_idx]
                        
                        # Find cell type(s) with maximum score
                        max_indices = np.where(satisfaction_matrix[:, cluster_idx] == max_score)[0]
                        
                        # Handle ties
                        if len(max_indices) > 1:
                            annotation = 'Undecided'
                        else:
                            annotation = cell_types[max_indices[0]]
                        
                        # Apply cutoff thresholds
                        if cutoff != 'none':
                            if cutoff == 'mean':
                                threshold = np.mean(max_scores)
                            elif cutoff == '0.5':
                                threshold = 0.5
                            else:
                                try:
                                    threshold = float(cutoff)
                                except ValueError:
                                    print(f"Invalid cutoff: {cutoff}, using 0.5")
                                    threshold = 0.5
                            
                            if max_score < threshold:
                                annotation = 'Undecided'
                        
                        cluster_annotations.append(annotation)
                        cluster_scores.append(max_score)
                    
                    # Map cluster annotations to cells
                    cell_annotations = []
                    cell_scores = []
                    
                    for cell_cluster in adata.obs[clustering]:
                        try:
                            cluster_idx = clusters.index(cell_cluster)
                            cell_annotations.append(cluster_annotations[cluster_idx])
                            cell_scores.append(cluster_scores[cluster_idx])
                        except ValueError:
                            cell_annotations.append('Undecided')
                            cell_scores.append(0.0)
                    
                    # Add annotations to adata
                    annotation_name = f"{clustering}_{source_name}_{cutoff}"
                    score_name = f"{clustering}_{source_name}_{cutoff}_score"
                    
                    # Clean column names
                    annotation_name = annotation_name.replace(" ", ".").replace("/", ".")
                    score_name = score_name.replace(" ", ".").replace("/", ".")
                    
                    adata.obs[annotation_name] = cell_annotations
                    adata.obs[score_name] = cell_scores
                    
                    print(f"    Added annotation column: {annotation_name}")
                    
                    # Print summary
                    unique_annotations = pd.Series(cell_annotations).value_counts()
                    print(f"      Annotation distribution: {dict(unique_annotations)}")
    
    print("Density scoring completed!")
    return adata if copy else None


class DeepModel(nn.Module):
    """Deep neural network model for cell type classification."""
    
    def __init__(self, input_shape: int, class_shape: int):
        super().__init__()
        self.activation = nn.LeakyReLU()
        self.output_activation = nn.Softmax(dim=1)
        self.fc1 = nn.Linear(input_shape, 256)
        self.fc2 = nn.Linear(256, 128)
        self.fc3 = nn.Linear(128, class_shape)

    def forward(self, x):
        x = self.activation(self.fc1(x))
        x = self.activation(self.fc2(x))
        x = self.output_activation(self.fc3(x))
        return x


def get_stats(preds: torch.Tensor, labels: torch.Tensor, num_possible_labels: int) -> Dict[str, float]:
    """Calculate classification statistics."""
    all_stats = {}
    
    # Accuracy
    acc = torchmetrics.Accuracy(task='multiclass', num_classes=num_possible_labels)
    all_stats['acc'] = acc(preds, labels).item()
    
    # F1Score
    f1 = torchmetrics.F1Score(task='multiclass', num_classes=num_possible_labels)
    all_stats['f1'] = f1(preds, labels).item()
    
    # Precision and Recall
    all_stats['p'] = torchmetrics.Precision(task='multiclass', num_classes=num_possible_labels)(preds, labels).item()
    all_stats['r'] = torchmetrics.Recall(task='multiclass', num_classes=num_possible_labels)(preds, labels).item()
    
    return all_stats


def dgscrna_annotate(
    adata: sc.AnnData,
    annotation_column: str,
    epochs: int = 10,
    batch_size: int = 256,
    learning_rate: float = 1e-3,
    train_split: float = 0.9,
    confidence_threshold: float = 0.9,
    copy: bool = False
) -> Optional[sc.AnnData]:
    """
    Perform deep learning-based cell type annotation.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    annotation_column : str
        Column in adata.obs containing initial annotations (with 'Undecided' for unknown cells)
    epochs : int, default 10
        Number of training epochs
    batch_size : int, default 256
        Batch size for training
    learning_rate : float, default 1e-3
        Learning rate for optimizer
    train_split : float, default 0.9
        Fraction of labeled cells to use for training (rest for validation)
    confidence_threshold : float, default 0.9
        Minimum confidence for accepting predictions
    copy : bool, default False
        Return a copy instead of writing to adata
        
    Returns
    -------
    AnnData or None
        Returns AnnData if copy=True, otherwise modifies adata in place
    """
    if copy:
        adata = adata.copy()
    
    if annotation_column not in adata.obs.columns:
        raise ValueError(f"Annotation column '{annotation_column}' not found in adata.obs")
    
    # Get cell types and handle undecided cells
    cell_types = adata.obs[annotation_column]
    
    if (cell_types == 'Undecided').sum() == len(cell_types):
        print('All cells are undecided/unknown')
        adata.obs[f'{annotation_column}_DGscRNA'] = ['Undecided'] * len(cell_types)
        return adata if copy else None
        
    if (cell_types == 'Undecided').sum() == 0:
        print('No undecided/unknown cells')
        adata.obs[f'{annotation_column}_DGscRNA'] = cell_types
        return adata if copy else None
    
    # Split cells into decided and undecided
    decided_mask = cell_types != 'Undecided'
    undecided_mask = cell_types == 'Undecided'
    
    decided_cells = adata[decided_mask]
    undecided_cells = adata[undecided_mask]
    
    print(f'Num decided: {decided_cells.n_obs}')
    print(f'Num undecided: {undecided_cells.n_obs}')
    
    # Create label mapping
    unique_types = sorted(list(set(cell_types) - {'Undecided'}))
    label_to_idx = {label: idx for idx, label in enumerate(unique_types)}
    idx_to_label = {idx: label for label, idx in label_to_idx.items()}
    
    # Prepare training data
    X_labeled = torch.tensor(decided_cells.X, dtype=torch.float)
    if hasattr(X_labeled, 'toarray'):
        X_labeled = torch.tensor(X_labeled.toarray(), dtype=torch.float)
    
    y_labeled = torch.tensor([label_to_idx[label] for label in decided_cells.obs[annotation_column]], dtype=torch.long)
    
    # Prepare test data (undecided cells)
    X_unlabeled = torch.tensor(undecided_cells.X, dtype=torch.float)
    if hasattr(X_unlabeled, 'toarray'):
        X_unlabeled = torch.tensor(X_unlabeled.toarray(), dtype=torch.float)
    
    # Create datasets
    labeled_dataset = data_utils.TensorDataset(X_labeled, y_labeled)
    
    # Split labeled data into train/validation
    train_size = int(len(labeled_dataset) * train_split)
    val_size = len(labeled_dataset) - train_size
    
    train_dataset, val_dataset = torch.utils.data.random_split(
        labeled_dataset, [train_size, val_size], 
        generator=torch.Generator().manual_seed(42)
    )
    
    # Create data loaders
    train_loader = data_utils.DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    val_loader = data_utils.DataLoader(val_dataset, batch_size=batch_size, shuffle=False)
    
    # Initialize model
    input_dim = X_labeled.shape[1]
    num_classes = len(unique_types)
    model = DeepModel(input_dim, num_classes)
    
    # Setup training
    criterion = nn.CrossEntropyLoss()
    optimizer = torch.optim.Adamax(model.parameters(), lr=learning_rate)
    
    # Training loop
    print("Starting training...")
    for epoch in range(epochs):
        model.train()
        total_loss = 0
        all_preds = []
        all_labels = []
        
        for batch_x, batch_y in train_loader:
            # Forward pass
            outputs = model(batch_x)
            loss = criterion(outputs, batch_y)
            
            # Backward pass
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
            preds = outputs.argmax(dim=1)
            all_preds.append(preds)
            all_labels.append(batch_y)
        
        # Calculate training stats
        all_preds = torch.cat(all_preds)
        all_labels = torch.cat(all_labels)
        train_stats = get_stats(all_preds, all_labels, num_classes)
        
        print(f"Epoch {epoch+1}/{epochs} | Train - ACC: {train_stats['acc']:.4f}, "
              f"F1: {train_stats['f1']:.4f}, Loss: {total_loss:.4f}")
    
    # Validation
    print("Validating...")
    model.eval()
    val_preds = []
    val_labels = []
    val_probs = []
    
    with torch.no_grad():
        for batch_x, batch_y in val_loader:
            outputs = model(batch_x)
            preds = outputs.argmax(dim=1)
            val_preds.append(preds)
            val_labels.append(batch_y)
            val_probs.append(outputs)
    
    val_preds = torch.cat(val_preds)
    val_labels = torch.cat(val_labels)
    val_probs = torch.cat(val_probs)
    
    val_stats = get_stats(val_preds, val_labels, num_classes)
    print(f"Validation - ACC: {val_stats['acc']:.4f}, F1: {val_stats['f1']:.4f}")
    
    # Predict on undecided cells
    print("Predicting undecided cells...")
    model.eval()
    predictions = []
    confidences = []
    
    with torch.no_grad():
        # Process in batches
        for i in range(0, len(X_unlabeled), batch_size):
            batch_x = X_unlabeled[i:i+batch_size]
            outputs = model(batch_x)
            batch_preds = outputs.argmax(dim=1)
            batch_conf = torch.max(outputs, dim=1)[0]
            
            predictions.extend(batch_preds.numpy())
            confidences.extend(batch_conf.numpy())
    
    # Apply confidence threshold
    final_annotations = cell_types.copy()
    
    for i, (pred_idx, conf) in enumerate(zip(predictions, confidences)):
        if conf >= confidence_threshold:
            predicted_label = idx_to_label[pred_idx]
            final_annotations.iloc[undecided_mask.sum() if i < undecided_mask.sum() else len(final_annotations)] = predicted_label
        else:
            final_annotations.iloc[undecided_mask.sum() if i < undecided_mask.sum() else len(final_annotations)] = 'Unknown'
    
    # Store results
    result_col = f'{annotation_column}_DGscRNA'
    adata.obs[result_col] = final_annotations
    
    print(f"Final annotation distribution:")
    print(adata.obs[result_col].value_counts())
    
    return adata if copy else None


def find_all_markers(
    adata: sc.AnnData,
    cluster_column: str,
    method: str = 'wilcoxon',
    min_logfc: float = 0.25,
    min_pct: float = 0.1,
    max_pval: float = 0.05,
    copy: bool = False
) -> Dict[str, pd.DataFrame]:
    """
    Find differentially expressed genes for each cluster (equivalent to Seurat's FindAllMarkers).
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    cluster_column : str
        Column name in adata.obs containing cluster assignments
    method : str, default 'wilcoxon'
        Test method for differential expression ('wilcoxon', 't-test', 't-test_overestim_var')
    min_logfc : float, default 0.25
        Minimum log fold change threshold
    min_pct : float, default 0.1
        Minimum percentage of cells expressing the gene in either group
    max_pval : float, default 0.05
        Maximum adjusted p-value threshold
    copy : bool, default False
        Return a copy instead of writing to adata
        
    Returns
    -------
    Dict[str, pd.DataFrame]
        Dictionary with cluster IDs as keys and marker DataFrames as values
    """
    if copy:
        adata = adata.copy()
    
    if cluster_column not in adata.obs.columns:
        raise ValueError(f"Cluster column '{cluster_column}' not found in adata.obs")
    
    print(f"Finding markers for clustering: {cluster_column}")
    
    # Get unique clusters
    clusters = sorted(adata.obs[cluster_column].unique())
    all_markers = {}
    
    for cluster in clusters:
        print(f"  Processing cluster {cluster}...")
        
        # Create mask for current cluster vs all others
        cluster_mask = adata.obs[cluster_column] == cluster
        
        if cluster_mask.sum() < 3:  # Skip clusters with too few cells
            print(f"    Skipping cluster {cluster} (only {cluster_mask.sum()} cells)")
            continue
        
        # Perform differential expression test
        try:
            # Use scanpy's rank_genes_groups for one cluster vs rest
            adata_temp = adata.copy()
            adata_temp.obs['temp_group'] = 'other'
            adata_temp.obs.loc[cluster_mask, 'temp_group'] = 'target'
            
            # Run differential expression
            sc.tl.rank_genes_groups(
                adata_temp,
                groupby='temp_group',
                groups=['target'],
                reference='other',
                method=method,
                use_raw=False,
                corr_method='benjamini-hochberg'
            )
            
            # Extract results
            result = sc.get.rank_genes_groups_df(
                adata_temp,
                group='target'
            )
            
            # Filter results based on thresholds
            filtered_result = result[
                (result['logfoldchanges'] >= min_logfc) &
                (result['pvals_adj'] <= max_pval)
            ].copy()
            
            # Add additional statistics
            if len(filtered_result) > 0:
                # Calculate percentage of cells expressing each gene
                cluster_cells = adata[cluster_mask]
                other_cells = adata[~cluster_mask]
                
                pct_1_list = []
                pct_2_list = []
                avg_log2fc_list = []
                
                for gene in filtered_result['names']:
                    if gene in adata.var_names:
                        # Get expression data
                        cluster_expr = cluster_cells[:, gene].X
                        other_expr = other_cells[:, gene].X
                        
                        if hasattr(cluster_expr, 'toarray'):
                            cluster_expr = cluster_expr.toarray().flatten()
                            other_expr = other_expr.toarray().flatten()
                        
                        # Calculate percentages
                        pct_1 = (cluster_expr > 0).mean()
                        pct_2 = (other_expr > 0).mean()
                        
                        # Calculate average log2 fold change
                        mean_cluster = np.mean(cluster_expr)
                        mean_other = np.mean(other_expr)
                        
                        # Add small pseudocount to avoid log(0)
                        avg_log2fc = np.log2((mean_cluster + 1e-9) / (mean_other + 1e-9))
                        
                        pct_1_list.append(pct_1)
                        pct_2_list.append(pct_2)
                        avg_log2fc_list.append(avg_log2fc)
                    else:
                        pct_1_list.append(0)
                        pct_2_list.append(0)
                        avg_log2fc_list.append(0)
                
                # Add to result dataframe
                filtered_result['pct_1'] = pct_1_list
                filtered_result['pct_2'] = pct_2_list
                filtered_result['avg_log2FC'] = avg_log2fc_list
                filtered_result['cluster'] = cluster
                
                # Filter by percentage threshold
                filtered_result = filtered_result[
                    (filtered_result['pct_1'] >= min_pct) |
                    (filtered_result['pct_2'] >= min_pct)
                ]
                
                # Rename columns to match Seurat output (keep avg_log2FC, rename logfoldchanges to match)
                filtered_result = filtered_result.rename(columns={
                    'names': 'gene',
                    'pvals': 'p_val',
                    'pvals_adj': 'p_val_adj',
                    'logfoldchanges': 'logfoldchanges_orig'  # Keep original but rename to avoid confusion
                })
                
                # Sort by average log fold change
                filtered_result = filtered_result.sort_values('avg_log2FC', ascending=False)
                
                all_markers[str(cluster)] = filtered_result
                print(f"    Found {len(filtered_result)} markers for cluster {cluster}")
            else:
                print(f"    No significant markers found for cluster {cluster}")
                all_markers[str(cluster)] = pd.DataFrame()
                
        except Exception as e:
            print(f"    Error processing cluster {cluster}: {e}")
            all_markers[str(cluster)] = pd.DataFrame()
    
    print(f"Completed marker finding for {len(all_markers)} clusters")
    return all_markers


def density_score_with_deg(
    adata: sc.AnnData,
    marker_sources: Dict[str, pd.DataFrame],
    cluster_columns: List[str],
    cutoffs: List[str] = ['0.5', 'mean', 'none'],
    epochs: int = 10,
    confidence_threshold: float = 0.9,
    copy: bool = False
) -> Optional[sc.AnnData]:
    """
    Calculate density scores and perform differential expression analysis.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    marker_sources : dict
        Dictionary of marker gene sets from different sources
    cluster_columns : list of str
        Column names in adata.obs containing cluster assignments
    cutoffs : list of str, default ['0.5', 'mean', 'none']
        Thresholds for density scoring
    epochs : int, default 10
        Number of training epochs for deep learning
    confidence_threshold : float, default 0.9
        Minimum confidence for accepting predictions
    copy : bool, default False
        Return a copy instead of writing to adata
        
    Returns
    -------
    AnnData or None
        Returns AnnData if copy=True, otherwise modifies adata in place
    """
    if copy:
        adata = adata.copy()
    
    print("Calculating density scores...")
    density_score(adata, marker_sources, cluster_columns, cutoffs)
    
    # Find annotation columns (those with density scores)
    annotation_columns = [col for col in adata.obs.columns 
                         if any(cutoff in col for cutoff in cutoffs) 
                         and any(source in col for source in marker_sources.keys())]
    
    print(f"Running DGscRNA annotation on {len(annotation_columns)} annotation combinations...")
    
    # Run deep learning annotation on each combination
    for annotation_col in annotation_columns:
        if 'DGscRNA' not in annotation_col:  # Avoid re-annotating already processed columns
            print(f"Processing {annotation_col}...")
            try:
                dgscrna_annotate(adata, annotation_col, epochs=epochs, 
                               confidence_threshold=confidence_threshold)
            except Exception as e:
                print(f"Error processing {annotation_col}: {e}")
                continue
    
    print("Performing differential expression analysis...")
    for annotation_col in annotation_columns:
        if 'DGscRNA' in annotation_col:  # Only process DGscRNA annotations
            cluster_column = annotation_col
            break
    else:
        print("No DGscRNA annotations found, skipping differential expression analysis")
        return adata if copy else None
    
    # Run differential expression analysis
    all_markers = find_all_markers(adata, cluster_column, copy=True)
    
    # Store results in adata
    for cluster, markers_df in all_markers.items():
        if not markers_df.empty:
            # Create a unique key for this cluster's markers
            key = f"DE_{cluster}"
            adata.uns[key] = markers_df
    
    print("Density scoring and differential expression analysis completed!")
    return adata if copy else None


def run_dgscrna_workflow(
    adata: sc.AnnData,
    marker_folder: Union[str, Path],
    cluster_columns: List[str],
    cutoffs: List[str] = ['0.5', 'mean', 'none'],
    epochs: int = 10,
    confidence_threshold: float = 0.9,
    copy: bool = False
) -> Optional[sc.AnnData]:
    """
    Run the complete DGscRNA workflow.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    marker_folder : str or Path
        Path to folder containing marker files
    cluster_columns : list of str
        Column names in adata.obs containing cluster assignments
    cutoffs : list of str, default ['0.5', 'mean', 'none']
        Thresholds for density scoring
    epochs : int, default 10
        Number of training epochs for deep learning
    confidence_threshold : float, default 0.9
        Minimum confidence for accepting predictions
    copy : bool, default False
        Return a copy instead of writing to adata
        
    Returns
    -------
    AnnData or None
        Returns AnnData if copy=True, otherwise modifies adata in place
    """
    if copy:
        adata = adata.copy()
    
    print("Loading marker gene sets...")
    marker_sources = load_markers(marker_folder)
    print(f"Loaded {len(marker_sources)} marker sources")
    
    print("Calculating density scores...")
    density_score(adata, marker_sources, cluster_columns, cutoffs)
    
    # Find annotation columns (those with density scores)
    annotation_columns = [col for col in adata.obs.columns 
                         if any(cutoff in col for cutoff in cutoffs) 
                         and any(source in col for source in marker_sources.keys())]
    
    print(f"Running DGscRNA annotation on {len(annotation_columns)} annotation combinations...")
    
    # Run deep learning annotation on each combination
    for annotation_col in annotation_columns:
        if 'DGscRNA' not in annotation_col:  # Avoid re-annotating already processed columns
            print(f"Processing {annotation_col}...")
            try:
                dgscrna_annotate(adata, annotation_col, epochs=epochs, 
                               confidence_threshold=confidence_threshold)
            except Exception as e:
                print(f"Error processing {annotation_col}: {e}")
                continue
    
    print("Performing differential expression analysis...")
    for annotation_col in annotation_columns:
        if 'DGscRNA' in annotation_col:  # Only process DGscRNA annotations
            cluster_column = annotation_col
            break
    else:
        print("No DGscRNA annotations found, skipping differential expression analysis")
        return adata if copy else None
    
    # Run differential expression analysis
    all_markers = find_all_markers(adata, cluster_column, copy=True)
    
    # Store results in adata
    for cluster, markers_df in all_markers.items():
        if not markers_df.empty:
            # Create a unique key for this cluster's markers
            key = f"DE_{cluster}"
            adata.uns[key] = markers_df
    
    print("DGscRNA workflow completed!")
    return adata if copy else None
