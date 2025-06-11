"""
scanpy-dgscrna: A scanpy extension for automated cell type annotation using density scoring and deep learning.

This package provides tools for:
- Loading marker gene sets from various sources
- Computing density scores with different thresholds
- Deep learning-based cell type annotation
- Integration with scanpy workflows
"""

from . import tools as tl
from . import plotting as pl
from . import preprocessing as pp

__version__ = "0.1.0"
__author__ = "DGscRNA Team"
__email__ = "your.email@example.com"

# Make main functions available at package level
from .tools import (
    load_markers,
    find_all_markers,
    density_score,
    dgscrna_annotate,
    run_dgscrna_workflow
)

__all__ = [
    "tl",
    "pl", 
    "pp",
    "load_markers",
    "find_all_markers",
    "density_score",
    "dgscrna_annotate",
    "run_dgscrna_workflow"
]
