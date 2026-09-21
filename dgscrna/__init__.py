"""
DGscRNA: Deep learning-guided single-cell RNA-seq cell type annotation
"""

__version__ = "2.0.0rc1"
__author__ = "DGscRNA Team"

from importlib import import_module
import warnings

_LEGACY = {
    "preprocess_adata": "preprocessing", "integrate_datasets": "preprocessing",
    "run_clustering": "clustering", "find_markers": "clustering",
    "score_cell_types": "marker_scoring", "load_marker_sets": "marker_scoring",
    "train_deep_model": "deep_learning", "predict_cell_types": "deep_learning",
    "run_dgscrna_pipeline": "utils",
}


def __getattr__(name):
    if name == "run_reference":
        return import_module(".reference.runner", __name__).run_reference
    if name in _LEGACY:
        warnings.warn(
            f"{name} uses the legacy simplified Python implementation. "
            "Use run_reference for the original R + Python DL workflow.",
            FutureWarning, stacklevel=2,
        )
        return getattr(import_module(f".core.{_LEGACY[name]}", __name__), name)
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = ["run_reference", "__version__"]
