"""
Deep learning module for DGscRNA package
"""

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import numpy as np
import pandas as pd
from typing import Dict, List, Optional, Tuple
import scanpy as sc
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score, precision_score, recall_score
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')

from ..models.deep_model import DeepModel

def get_stats(preds: torch.Tensor, labels: torch.Tensor, num_classes: int) -> Dict[str, float]:
    """
    Calculate evaluation statistics
    
    Parameters
    ----------
    preds : torch.Tensor
        Predicted class indices
    labels : torch.Tensor
        True labels
    num_classes : int
        Number of classes
        
    Returns
    -------
    Dict[str, float]
        Dictionary with evaluation metrics
    """
    # Convert to numpy for sklearn metrics
    preds_np = preds.cpu().numpy()
    labels_np = labels.cpu().numpy()
    
    # Calculate metrics
    acc = accuracy_score(labels_np, preds_np)
    f1 = f1_score(labels_np, preds_np, average='weighted')
    precision = precision_score(labels_np, preds_np, average='weighted', zero_division=0)
    recall = recall_score(labels_np, preds_np, average='weighted', zero_division=0)
    
    return {
        'accuracy': acc,
        'f1_score': f1,
        'precision': precision,
        'recall': recall
    }

def prepare_training_data(
    adata,
    annotation_key: str,
    use_highly_variable: bool = True,
    test_size: float = 0.1,
    random_state: int = 42
) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, Dict]:
    """
    Prepare data for deep learning training
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    annotation_key : str
        Key in obs containing cell type annotations
    use_highly_variable : bool, default=True
        Whether to use only highly variable genes
    test_size : float, default=0.1
        Fraction of data to use for testing
    random_state : int, default=42
        Random state for reproducibility
        
    Returns
    -------
    Tuple
        (X_train, y_train, X_test, y_test, label_mapping)
    """
    
    # Check if annotation key exists
    if annotation_key not in adata.obs.columns:
        raise ValueError(f"Annotation key '{annotation_key}' not found in adata.obs")
    
    # Get expression data
    if use_highly_variable and 'highly_variable' in adata.var.columns:
        # Use only highly variable genes
        adata_subset = adata[:, adata.var['highly_variable']]
    else:
        adata_subset = adata
    
    # Get expression matrix
    if 'X_scaled' in adata_subset.obsm:
        X = adata_subset.obsm['X_scaled']
    else:
        X = adata_subset.X.toarray() if hasattr(adata_subset.X, 'toarray') else adata_subset.X
    
    # Ensure X is numpy array
    if not isinstance(X, np.ndarray):
        X = np.array(X)
    
    # Get annotations
    annotations = adata_subset.obs[annotation_key].values
    
    # Remove 'Undecided' cells
    valid_mask = annotations != 'Undecided'
    X = X[valid_mask]
    annotations = annotations[valid_mask]
    
    if len(X) == 0:
        raise ValueError("No valid cells found for training")
    
    # Create label mapping
    unique_labels = sorted(set(annotations))
    label_mapping = {label: i for i, label in enumerate(unique_labels)}
    reverse_mapping = {i: label for label, i in label_mapping.items()}
    
    # Convert labels to indices
    y = np.array([label_mapping[label] for label in annotations])
    
    # Split data
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state, stratify=y
    )
    
    # Convert to tensors
    X_train = torch.FloatTensor(X_train)
    y_train = torch.LongTensor(y_train)
    X_test = torch.FloatTensor(X_test)
    y_test = torch.LongTensor(y_test)
    
    return X_train, y_train, X_test, y_test, reverse_mapping

def train_deep_model(
    adata,
    annotation_key: str,
    epochs: int = 10,
    batch_size: int = 256,
    lr: float = 1e-3,
    use_highly_variable: bool = True,
    hidden_dims: Optional[List[int]] = None,
    device: Optional[str] = None,
    random_state: int = 42
) -> Dict:
    """
    Train deep learning model to refine cell type annotations
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    annotation_key : str
        Key in obs containing cell type annotations
    epochs : int, default=10
        Number of training epochs
    batch_size : int, default=256
        Batch size for training
    lr : float, default=1e-3
        Learning rate
    use_highly_variable : bool, default=True
        Whether to use only highly variable genes
    hidden_dims : List[int], optional
        Hidden layer dimensions
    device : str, optional
        Device to use ('cpu' or 'cuda')
    random_state : int, default=42
        Random state for reproducibility
        
    Returns
    -------
    Dict
        Dictionary containing training results and model
    """
    
    # Set device
    if device is None:
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
    
    # Set random seeds
    torch.manual_seed(random_state)
    np.random.seed(random_state)
    
    # Prepare data
    X_train, y_train, X_test, y_test, label_mapping = prepare_training_data(
        adata, annotation_key, use_highly_variable, random_state=random_state
    )
    
    # Create datasets
    train_dataset = TensorDataset(X_train, y_train)
    test_dataset = TensorDataset(X_test, y_test)
    
    # Create data loaders
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    test_loader = DataLoader(test_dataset, batch_size=batch_size, shuffle=False)
    
    # Initialize model
    input_dim = X_train.shape[1]
    num_classes = len(label_mapping)
    
    model = DeepModel(input_dim, num_classes, hidden_dims)
    model.to(device)
    
    # Loss and optimizer
    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adamax(model.parameters(), lr=lr)
    
    # Training loop
    train_losses = []
    train_metrics = []
    test_metrics = []
    
    print(f"Training model with {num_classes} classes on {device}")
    print(f"Input dimension: {input_dim}")
    print(f"Training samples: {len(X_train)}, Test samples: {len(X_test)}")
    
    for epoch in range(epochs):
        # Training phase
        model.train()
        total_loss = 0
        all_preds = []
        all_labels = []
        
        for batch_X, batch_y in train_loader:
            batch_X, batch_y = batch_X.to(device), batch_y.to(device)
            
            optimizer.zero_grad()
            outputs = model(batch_X)
            loss = criterion(outputs, batch_y)
            loss.backward()
            optimizer.step()
            
            total_loss += loss.item()
            preds = torch.argmax(outputs, dim=1)
            all_preds.append(preds.cpu())
            all_labels.append(batch_y.cpu())
        
        # Calculate training metrics
        all_preds = torch.cat(all_preds)
        all_labels = torch.cat(all_labels)
        train_stats = get_stats(all_preds, all_labels, num_classes)
        train_losses.append(total_loss / len(train_loader))
        train_metrics.append(train_stats)
        
        # Validation phase
        model.eval()
        all_test_preds = []
        all_test_labels = []
        all_test_probs = []
        
        with torch.no_grad():
            for batch_X, batch_y in test_loader:
                batch_X, batch_y = batch_X.to(device), batch_y.to(device)
                outputs = model(batch_X)
                preds = torch.argmax(outputs, dim=1)
                
                all_test_preds.append(preds.cpu())
                all_test_labels.append(batch_y.cpu())
                all_test_probs.append(outputs.cpu())
        
        all_test_preds = torch.cat(all_test_preds)
        all_test_labels = torch.cat(all_test_labels)
        all_test_probs = torch.cat(all_test_probs)
        
        test_stats = get_stats(all_test_preds, all_test_labels, num_classes)
        test_metrics.append(test_stats)
        
        # Print progress
        print(f"Epoch {epoch+1}/{epochs}")
        print(f"Train - Loss: {train_losses[-1]:.4f}, Acc: {train_stats['accuracy']:.4f}, F1: {train_stats['f1_score']:.4f}")
        print(f"Test  - Acc: {test_stats['accuracy']:.4f}, F1: {test_stats['f1_score']:.4f}")
        print("-" * 50)
    
    # Prepare results
    results = {
        'model': model,
        'label_mapping': label_mapping,
        'train_losses': train_losses,
        'train_metrics': train_metrics,
        'test_metrics': test_metrics,
        'input_dim': input_dim,
        'num_classes': num_classes,
        'device': device
    }
    
    return results

def predict_cell_types(
    model_results: Dict,
    adata,
    annotation_key: str,
    probability_threshold: float = 0.9,
    use_highly_variable: bool = True
) -> pd.Series:
    """
    Predict cell types using trained model
    
    Parameters
    ----------
    model_results : Dict
        Results from train_deep_model
    adata : AnnData
        Annotated data matrix
    annotation_key : str
        Key in obs containing original annotations
    probability_threshold : float, default=0.9
        Probability threshold for confident predictions
    use_highly_variable : bool, default=True
        Whether to use only highly variable genes
        
    Returns
    -------
    pd.Series
        Series with refined cell type annotations
    """
    
    model = model_results['model']
    label_mapping = model_results['label_mapping']
    device = model_results['device']
    
    # Get expression data
    if use_highly_variable and 'highly_variable' in adata.var.columns:
        adata_subset = adata[:, adata.var['highly_variable']]
    else:
        adata_subset = adata
    
    # Get expression matrix
    if 'X_scaled' in adata_subset.obsm:
        X = adata_subset.obsm['X_scaled']
    else:
        X = adata_subset.X.toarray() if hasattr(adata_subset.X, 'toarray') else adata_subset.X
    
    # Ensure X is numpy array
    if not isinstance(X, np.ndarray):
        X = np.array(X)
    
    # Get original annotations
    original_annotations = adata_subset.obs[annotation_key].values
    
    # Initialize predictions
    refined_annotations = original_annotations.copy()
    
    # Find cells with 'Undecided' annotations
    undecided_mask = original_annotations == 'Undecided'
    
    if undecided_mask.sum() == 0:
        print("No 'Undecided' cells found for refinement")
        return pd.Series(refined_annotations, index=adata_subset.obs.index)
    
    # Predict for undecided cells
    X_undecided = X[undecided_mask]
    X_tensor = torch.FloatTensor(X_undecided).to(device)
    
    model.eval()
    with torch.no_grad():
        probabilities = model.predict_proba(X_tensor)
        predictions = torch.argmax(probabilities, dim=1)
        max_probs = torch.max(probabilities, dim=1)[0]
    
    # Convert predictions back to labels
    predicted_labels = [label_mapping[pred.item()] for pred in predictions]
    max_probs = max_probs.cpu().numpy()
    
    # Apply probability threshold
    confident_mask = max_probs >= probability_threshold
    
    # Update annotations
    undecided_indices = np.where(undecided_mask)[0]
    for i, (idx, label, prob, confident) in enumerate(zip(undecided_indices, predicted_labels, max_probs, confident_mask)):
        if confident:
            refined_annotations[idx] = label
        else:
            refined_annotations[idx] = 'Unknown'
    
    return pd.Series(refined_annotations, index=adata_subset.obs.index) 