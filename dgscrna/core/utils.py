"""
Utility functions for DGscRNA package
"""

import pandas as pd
import numpy as np
import scanpy as sc
from typing import Dict, List, Optional, Union, Tuple
import warnings
warnings.filterwarnings('ignore')

from .preprocessing import preprocess_adata
from .clustering import run_clustering, find_markers
from .marker_scoring import load_marker_sets, score_cell_types
from .deep_learning import train_deep_model, predict_cell_types

def run_dgscrna_pipeline(
    adata,
    marker_folder: str,
    clustering_methods: List[str] = ['hdbscan'],
    cutoff_strategy: Union[str, List[str]] = 'mean',
    marker_sets_to_use: Optional[List[str]] = None,
    use_deep_learning: bool = True,
    dl_epochs: int = 10,
    dl_batch_size: int = 256,
    dl_learning_rate: float = 0.001,
    probability_threshold: float = 0.9,
    random_state: int = 42,
    **kwargs
):
    """
    Complete DGscRNA pipeline for cell type annotation
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix (should be preprocessed and normalized)
    marker_folder : str
        Path to folder containing marker CSV files
    clustering_methods : List[str], default=['hdbscan'] 
        List of clustering methods to use
    cutoff_strategy : Union[str, List[str]], default='mean'
        Cutoff strategy for marker scoring. Options: 'mean', '0.5', 'none', or list of these
    marker_sets_to_use : Optional[List[str]]
        Specific marker sets to use. If None, use all CSV files in marker_folder
    use_deep_learning : bool, default=True
        Whether to use deep learning refinement
    dl_epochs : int, default=10
        Number of training epochs for deep learning
    dl_batch_size : int, default=256
        Batch size for deep learning
    dl_learning_rate : float, default=0.001
        Learning rate for deep learning training
    probability_threshold : float, default=0.9
        Probability threshold for confident predictions
    random_state : int, default=42
        Random state for reproducibility
    **kwargs
        Additional arguments for preprocessing and clustering
        
    Returns
    -------
    AnnData
        Processed data with annotations
    """
    
    print("Starting DGscRNA pipeline...")
    print(f"Input data shape: {adata.shape}")
    
    # Step 1: Preprocessing (if not already done)
    if 'X_pca' not in adata.obsm:
        print("Running preprocessing...")
        adata = preprocess_adata(adata, random_state=random_state, **kwargs)
    
    # Step 2: Clustering
    print("Running clustering...")
    adata = run_clustering(adata, methods=clustering_methods, random_state=random_state, **kwargs)
    
    # Step 3: Find marker genes for each clustering method
    print("Finding marker genes...")
    for method in clustering_methods:
        cluster_key = f"{method}_clusters"
        if cluster_key in adata.obs.columns:
            adata = find_markers(adata, groupby=cluster_key, **kwargs)
        
    # Step 4: Load marker sets
    print("Loading marker sets...")
    marker_sets = load_marker_sets(marker_folder)
    print(f"Loaded {len(marker_sets)} marker sets")
    
    # Step 5: Score cell types for each clustering method and marker set
    print("Scoring cell types...")
    results = {}
    
    # Filter marker sets if specified
    if marker_sets_to_use:
        marker_sets = {k: v for k, v in marker_sets.items() if k in marker_sets_to_use}
    
    # Ensure cutoff_strategy is a list
    if isinstance(cutoff_strategy, str):
        cutoff_strategies = [cutoff_strategy]
    else:
        cutoff_strategies = cutoff_strategy
    
    for method in clustering_methods:
        cluster_key = f"{method}_clusters" if not method.endswith('_clusters') else method
        if cluster_key not in adata.obs.columns:
            continue
            
        for marker_set_name, markers_dict in marker_sets.items():
            print(f"Scoring {marker_set_name} with {method} clustering...")
            
            # Score cell types
            adata = score_cell_types(
                adata, 
                markers_dict, 
                cluster_key,
                marker_set_name=marker_set_name,
                cutoffs=cutoff_strategy if isinstance(cutoff_strategy, list) else [cutoff_strategy],
                random_state=random_state
            )
            
            # Store results
            results[f"{marker_set_name}_{method}"] = {
                'marker_set': marker_set_name,
                'clustering_method': method,
                'markers_dict': markers_dict
            }
    
    # Step 6: Deep learning refinement (if enabled)
    if use_deep_learning:
        print("Running deep learning refinement...")
        dl_results = {}
        
        for method in clustering_methods:
            cluster_key = f"{method}_clusters"
            if cluster_key not in adata.obs.columns:
                continue
                
            for marker_set_name, markers_dict in marker_sets.items():
                # Find the annotation key with cutoff
                annotation_keys = [col for col in adata.obs.columns 
                                 if col.startswith(f"{marker_set_name}_{method}_")]
                
                for annotation_key in annotation_keys:
                    print(f"Training deep learning model for {annotation_key}...")
                    
                    try:
                        # Train model
                        model_results = train_deep_model(
                            adata,
                            annotation_key=annotation_key,
                            batch_size=dl_batch_size,
                            epochs=dl_epochs,
                            lr=dl_learning_rate,
                            random_state=random_state
                        )
                        
                        # Predict refined annotations
                        refined_annotations = predict_cell_types(
                            model_results,
                            adata,
                            annotation_key=annotation_key,
                            probability_threshold=probability_threshold
                        )
                        
                        # Create new annotation key for deep learning results
                        dl_annotation_key = f"{annotation_key}_DL"
                        adata.obs[dl_annotation_key] = refined_annotations
                        
                        # Store results
                        dl_results[annotation_key] = {
                            'model_results': model_results,
                            'original_key': annotation_key,
                            'refined_key': dl_annotation_key
                        }
                        
                    except Exception as e:
                        print(f"Error training model for {annotation_key}: {e}")
                        continue
        
        # Add optimal annotation based on combined F1 scores
        if dl_results:
            print("Determining optimal annotation based on combined F1 scores...")
            
            # Calculate combined scores for each annotation method
            scores_data = []
            combined_scores = {}
            
            for annotation_key, results in dl_results.items():
                model_results = results['model_results']
                train_metrics = model_results['train_metrics']
                test_metrics = model_results['test_metrics']
                
                # Get final epoch scores
                final_train_f1 = train_metrics[-1]['f1_score']
                final_test_f1 = test_metrics[-1]['f1_score']
                final_train_acc = train_metrics[-1]['accuracy']
                final_test_acc = test_metrics[-1]['accuracy']
                
                # Calculate combined F1 score
                combined_f1 = (final_test_f1 + final_train_f1) / 2
                combined_scores[annotation_key] = combined_f1
                
                # Store detailed scores for summary table
                scores_data.append({
                    'annotation_method': annotation_key,
                    'train_f1': final_train_f1,
                    'test_f1': final_test_f1,
                    'combined_f1': combined_f1,
                    'train_accuracy': final_train_acc,
                    'test_accuracy': final_test_acc,
                    'train_precision': train_metrics[-1]['precision'],
                    'test_precision': test_metrics[-1]['precision'],
                    'train_recall': train_metrics[-1]['recall'],
                    'test_recall': test_metrics[-1]['recall']
                })
            
            # Find the annotation method with highest combined F1 score
            best_annotation_key = max(combined_scores, key=combined_scores.get)
            best_combined_f1 = combined_scores[best_annotation_key]
            
            print(f"Best annotation method: {best_annotation_key}")
            print(f"Best combined F1 score: {best_combined_f1:.4f}")
            
            # Get the optimal annotations from the best method
            best_refined_key = dl_results[best_annotation_key]['refined_key']
            optimal_annotations = adata.obs[best_refined_key].copy()
            
            # Add optimal annotation to adata
            adata.obs['optimal_annotation'] = optimal_annotations
            
            # Create and store training scores summary
            training_scores_summary = pd.DataFrame(scores_data)
            training_scores_summary = training_scores_summary.sort_values('combined_f1', ascending=False)
            training_scores_summary['rank'] = range(1, len(training_scores_summary) + 1)
            
            # Reorder columns for better readability
            column_order = [
                'rank', 'annotation_method', 'combined_f1', 'train_f1', 'test_f1',
                'train_accuracy', 'test_accuracy', 'train_precision', 'test_precision',
                'train_recall', 'test_recall'
            ]
            training_scores_summary = training_scores_summary[column_order]
            
            # Store the training scores in adata.uns for easy access
            adata.uns['training_scores_summary'] = training_scores_summary
            
            print("Added 'optimal_annotation' column to adata.obs")
            print("Stored training scores summary in adata.uns['training_scores_summary']")
            
            # Display summary
            print("\nTraining Scores Summary (top 3):")
            print(training_scores_summary.head(3).round(4))
    
    print("DGscRNA pipeline completed successfully!")
    return adata

def summarize_results(pipeline_results: Dict) -> pd.DataFrame:
    """
    Summarize pipeline results
    
    Parameters
    ----------
    pipeline_results : Dict
        Results from run_dgscrna_pipeline
        
    Returns
    -------
    pd.DataFrame
        Summary of results
    """
    
    adata = pipeline_results.copy()
    results = []
    
    # Get all annotation columns
    annotation_cols = [col for col in adata.obs.columns if '_DL' in col or any(method in col for method in ['leiden', 'hdbscan', 'kmeans']) or 'optimal_annotation' in col]
    
    for col in annotation_cols:
        # Count cell types
        cell_type_counts = adata.obs[col].value_counts()
        
        # Calculate statistics
        total_cells = len(adata.obs[col])
        undecided_cells = cell_type_counts.get('Undecided', 0)
        unknown_cells = cell_type_counts.get('Unknown', 0)
        decided_cells = total_cells - undecided_cells - unknown_cells
        
        results.append({
            'annotation_method': col,
            'total_cells': total_cells,
            'decided_cells': decided_cells,
            'undecided_cells': undecided_cells,
            'unknown_cells': unknown_cells,
            'decision_rate': decided_cells / total_cells,
            'num_cell_types': len(cell_type_counts) - (1 if 'Undecided' in cell_type_counts else 0) - (1 if 'Unknown' in cell_type_counts else 0)
        })
    
    return pd.DataFrame(results)

def plot_results(adata, annotation_cols: List[str] = None, n_cols: int = 2, figsize: tuple = (15, 10)):
    """
    Plot annotation results
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    annotation_cols : List[str], optional
        List of annotation columns to plot. If None, automatically detect annotation columns.
    n_cols : int, default=2
        Number of columns in the plot grid
    figsize : tuple, default=(15, 10)
        Figure size
    """
    
    import matplotlib.pyplot as plt
    
    # Auto-detect annotation columns if not provided
    if annotation_cols is None:
        annotation_cols = [col for col in adata.obs.columns 
                          if any(keyword in col.lower() for keyword in 
                               ['_DL', 'annotation', 'leiden', 'hdbscan', 'kmeans', 'cluster'])]
        if not annotation_cols:
            print("No annotation columns found. Please specify annotation_cols parameter.")
            return
    
    # Filter to only existing columns
    existing_cols = [col for col in annotation_cols if col in adata.obs.columns]
    if not existing_cols:
        print(f"None of the specified columns {annotation_cols} exist in adata.obs")
        return
    
    if len(existing_cols) < len(annotation_cols):
        missing_cols = set(annotation_cols) - set(existing_cols)
        print(f"Warning: Some columns not found: {missing_cols}")
    
    annotation_cols = existing_cols
    
    # Check if UMAP coordinates exist
    if 'X_umap' not in adata.obsm:
        print("UMAP coordinates not found. Computing UMAP...")
        try:
            sc.tl.umap(adata)
        except Exception as e:
            print(f"Error computing UMAP: {e}")
            print("Please ensure your data has been preprocessed with PCA.")
            return
    
    n_plots = len(annotation_cols)
    n_rows = (n_plots + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    
    # Handle different subplot layouts
    if n_plots == 1:
        axes = [axes]
    elif n_rows == 1:
        axes = [axes] if n_cols == 1 else axes
    else:
        axes = axes.flatten()
    
    # Ensure axes is always a list
    if not isinstance(axes, (list, np.ndarray)):
        axes = [axes]
    
    # Plot each annotation
    for i, col in enumerate(annotation_cols):
        if i < len(axes):
            try:
                sc.pl.umap(adata, color=col, ax=axes[i], show=False, title=col, frameon=False)
            except Exception as e:
                print(f"Error plotting {col}: {e}")
                axes[i].text(0.5, 0.5, f"Error plotting\n{col}", 
                           ha='center', va='center', transform=axes[i].transAxes)
                axes[i].set_title(col)
    
    # Hide empty subplots
    for i in range(n_plots, len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    plt.show() 