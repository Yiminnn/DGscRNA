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
    copy: bool = False
) -> Optional[sc.AnnData]:
    """
    Calculate density scores for cell type annotation.
    
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
    copy : bool, default False
        Return a copy instead of writing to adata
        
    Returns
    -------
    AnnData or None
        Returns AnnData if copy=True, otherwise modifies adata in place
    """
    if copy:
        adata = adata.copy()
    
    # Ensure we have the required genes
    all_markers = set()
    for source_df in marker_sources.values():
        for celltype in source_df.columns:
            markers = source_df[celltype].dropna().tolist()
            all_markers.update(markers)
    
    # Filter to genes present in the data
    available_markers = list(set(all_markers) & set(adata.var_names))
    print(f"Found {len(available_markers)} marker genes in data out of {len(all_markers)} total markers")
    
    # Calculate scores for each clustering method and marker source
    for clustering in cluster_columns:
        if clustering not in adata.obs.columns:
            print(f"Clustering '{clustering}' not found in adata.obs, skipping...")
            continue
            
        clusters = adata.obs[clustering].unique()
        
        for source_name, source_df in marker_sources.items():
            print(f"Processing {clustering} with {source_name} markers...")
            
            for celltype in source_df.columns:
                markers = source_df[celltype].dropna().tolist()
                markers = [m for m in markers if m in adata.var_names]
                
                if not markers:
                    continue
                
                # Get marker expression for all cells
                marker_expr = adata[:, markers].X
                if hasattr(marker_expr, 'toarray'):
                    marker_expr = marker_expr.toarray()
                
                # Calculate mean expression per cell for these markers
                cell_scores = np.mean(marker_expr, axis=1)
                
                for cutoff in cutoffs:
                    if cutoff == 'none':
                        # Use raw scores
                        final_scores = cell_scores
                    elif cutoff == 'mean':
                        # Use global mean as threshold
                        threshold = np.mean(cell_scores)
                        final_scores = (cell_scores >= threshold).astype(float)
                    elif cutoff == '0.5':
                        # Use 0.5 as threshold (assuming normalized data)
                        final_scores = (cell_scores >= 0.5).astype(float)
                    else:
                        # Try to parse as float
                        try:
                            threshold = float(cutoff)
                            final_scores = (cell_scores >= threshold).astype(float)
                        except ValueError:
                            print(f"Invalid cutoff: {cutoff}, skipping...")
                            continue
                    
                    # Store scores
                    score_name = f"{clustering}_{source_name}_{celltype}_{cutoff}"
                    adata.obs[score_name] = final_scores
    
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
    
    print("DGscRNA workflow completed!")
    return adata if copy else None
