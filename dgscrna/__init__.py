"""
DGscRNA: Deep learning-guided single-cell RNA-seq cell type annotation
"""

__version__ = "1.0"
__author__ = "DGscRNA Team"

from .core.preprocessing import preprocess_adata, integrate_datasets
from .core.clustering import run_clustering, find_markers
from .core.marker_scoring import score_cell_types, load_marker_sets
from .core.deep_learning import train_deep_model, predict_cell_types
from .core.utils import run_dgscrna_pipeline

__all__ = [
    "preprocess_adata",
    "integrate_datasets", 
    "run_clustering",
    "find_markers",
    "score_cell_types",
    "load_marker_sets",
    "train_deep_model",
    "predict_cell_types",
    "run_dgscrna_pipeline",
] 