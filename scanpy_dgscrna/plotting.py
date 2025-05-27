"""
Plotting module for scanpy-dgscrna package.

Contains visualization functions for annotation results and quality control.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from typing import Optional, List, Union, Tuple


def plot_annotation_comparison(
    adata: sc.AnnData,
    original_column: str,
    annotated_column: str,
    figsize: Tuple[int, int] = (12, 5),
    save: Optional[str] = None
) -> None:
    """
    Plot comparison between original and DGscRNA annotations.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    original_column : str
        Column name for original annotations
    annotated_column : str
        Column name for DGscRNA annotations
    figsize : tuple, default (12, 5)
        Figure size
    save : str, optional
        Path to save the plot
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    
    # Original annotations
    original_counts = adata.obs[original_column].value_counts()
    axes[0].pie(original_counts.values, labels=original_counts.index, autopct='%1.1f%%')
    axes[0].set_title(f'Original: {original_column}')
    
    # DGscRNA annotations
    annotated_counts = adata.obs[annotated_column].value_counts()
    axes[1].pie(annotated_counts.values, labels=annotated_counts.index, autopct='%1.1f%%')
    axes[1].set_title(f'DGscRNA: {annotated_column}')
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    plt.show()


def plot_density_scores(
    adata: sc.AnnData,
    score_columns: List[str],
    cluster_column: str,
    figsize: Tuple[int, int] = (15, 10),
    save: Optional[str] = None
) -> None:
    """
    Plot density scores across clusters.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    score_columns : list of str
        Column names containing density scores
    cluster_column : str
        Column name for cluster assignments
    figsize : tuple, default (15, 10)
        Figure size
    save : str, optional
        Path to save the plot
    """
    n_scores = len(score_columns)
    n_cols = min(3, n_scores)
    n_rows = (n_scores + n_cols - 1) // n_cols
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    if n_scores == 1:
        axes = [axes]
    elif n_rows == 1:
        axes = axes.reshape(1, -1)
    
    for i, score_col in enumerate(score_columns):
        row = i // n_cols
        col = i % n_cols
        ax = axes[row, col] if n_rows > 1 else axes[col]
        
        # Create boxplot
        score_data = []
        cluster_labels = []
        
        for cluster in adata.obs[cluster_column].unique():
            cluster_mask = adata.obs[cluster_column] == cluster
            scores = adata.obs.loc[cluster_mask, score_col]
            score_data.extend(scores.tolist())
            cluster_labels.extend([str(cluster)] * len(scores))
        
        df_plot = pd.DataFrame({
            'Score': score_data,
            'Cluster': cluster_labels
        })
        
        sns.boxplot(data=df_plot, x='Cluster', y='Score', ax=ax)
        ax.set_title(score_col)
        ax.tick_params(axis='x', rotation=45)
    
    # Hide empty subplots
    for i in range(n_scores, n_rows * n_cols):
        row = i // n_cols
        col = i % n_cols
        if n_rows > 1:
            axes[row, col].set_visible(False)
        else:
            axes[col].set_visible(False)
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    plt.show()


def plot_umap_annotations(
    adata: sc.AnnData,
    annotation_columns: List[str],
    basis: str = 'umap',
    figsize: Tuple[int, int] = (15, 5),
    save: Optional[str] = None
) -> None:
    """
    Plot UMAP colored by different annotations.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    annotation_columns : list of str
        Column names for annotations to plot
    basis : str, default 'umap'
        Basis for plotting (umap, tsne, etc.)
    figsize : tuple, default (15, 5)
        Figure size
    save : str, optional
        Path to save the plot
    """
    n_plots = len(annotation_columns)
    fig, axes = plt.subplots(1, n_plots, figsize=figsize)
    if n_plots == 1:
        axes = [axes]
    
    for i, col in enumerate(annotation_columns):
        sc.pl.embedding(adata, basis=basis, color=col, ax=axes[i], 
                       show=False, frameon=False, title=col)
    
    plt.tight_layout()
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    plt.show()


def plot_marker_expression(
    adata: sc.AnnData,
    markers: List[str],
    groupby: str,
    figsize: Tuple[int, int] = (12, 8),
    save: Optional[str] = None
) -> None:
    """
    Plot marker gene expression across groups.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    markers : list of str
        List of marker genes to plot
    groupby : str
        Column name to group by
    figsize : tuple, default (12, 8)
        Figure size
    save : str, optional
        Path to save the plot
    """
    # Filter markers present in data
    available_markers = [m for m in markers if m in adata.var_names]
    
    if not available_markers:
        print("No markers found in data")
        return
    
    # Create dotplot
    sc.pl.dotplot(adata, available_markers, groupby=groupby, 
                  figsize=figsize, show=False)
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    plt.show()


def plot_annotation_sankey(
    adata: sc.AnnData,
    original_column: str,
    annotated_column: str,
    save: Optional[str] = None
) -> None:
    """
    Plot Sankey diagram showing annotation transitions.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    original_column : str
        Column name for original annotations
    annotated_column : str
        Column name for DGscRNA annotations
    save : str, optional
        Path to save the plot
    """
    try:
        import plotly.graph_objects as go
        from plotly.offline import plot
    except ImportError:
        print("plotly is required for Sankey diagrams. Install with: pip install plotly")
        return
    
    # Create transition matrix
    transition_df = adata.obs.groupby([original_column, annotated_column]).size().reset_index(name='count')
    
    # Create nodes
    original_labels = adata.obs[original_column].unique()
    annotated_labels = adata.obs[annotated_column].unique()
    
    all_labels = list(original_labels) + list(annotated_labels)
    label_to_idx = {label: idx for idx, label in enumerate(all_labels)}
    
    # Create links
    source = []
    target = []
    value = []
    
    for _, row in transition_df.iterrows():
        source.append(label_to_idx[row[original_column]])
        target.append(label_to_idx[row[annotated_column]] + len(original_labels))
        value.append(row['count'])
    
    # Create Sankey diagram
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=all_labels,
            color="blue"
        ),
        link=dict(
            source=source,
            target=target,
            value=value
        )
    )])
    
    fig.update_layout(title_text="Annotation Transitions", font_size=10)
    
    if save:
        fig.write_html(save)
    else:
        plot(fig)


def plot_confusion_matrix(
    adata: sc.AnnData,
    true_column: str,
    pred_column: str,
    figsize: Tuple[int, int] = (10, 8),
    save: Optional[str] = None
) -> None:
    """
    Plot confusion matrix for annotation comparison.
    
    Parameters
    ----------
    adata : AnnData
        Annotated data matrix
    true_column : str
        Column name for true annotations
    pred_column : str
        Column name for predicted annotations
    figsize : tuple, default (10, 8)
        Figure size
    save : str, optional
        Path to save the plot
    """
    from sklearn.metrics import confusion_matrix
    
    # Filter out undecided/unknown cells for comparison
    mask = (adata.obs[true_column] != 'Undecided') & (adata.obs[true_column] != 'Unknown')
    
    y_true = adata.obs.loc[mask, true_column]
    y_pred = adata.obs.loc[mask, pred_column]
    
    # Get unique labels
    labels = sorted(list(set(y_true.unique()) | set(y_pred.unique())))
    
    # Compute confusion matrix
    cm = confusion_matrix(y_true, y_pred, labels=labels)
    
    # Plot
    plt.figure(figsize=figsize)
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', 
                xticklabels=labels, yticklabels=labels)
    plt.title(f'Confusion Matrix: {true_column} vs {pred_column}')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)
    
    if save:
        plt.savefig(save, dpi=300, bbox_inches='tight')
    plt.show()
