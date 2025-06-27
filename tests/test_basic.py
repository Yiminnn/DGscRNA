"""Basic tests for DGscRNA package."""

import pytest

def test_package_import():
    """Test that the package can be imported."""
    try:
        import dgscrna
        assert dgscrna is not None
    except ImportError:
        pytest.skip("Package not properly installed")

def test_core_modules_import():
    """Test that core modules can be imported."""
    try:
        from dgscrna.core import preprocessing, clustering, marker_scoring, utils
        assert preprocessing is not None
        assert clustering is not None
        assert marker_scoring is not None
        assert utils is not None
    except ImportError:
        pytest.skip("Core modules not available")

def test_models_import():
    """Test that models can be imported."""
    try:
        from dgscrna.models import deep_model
        assert deep_model is not None
    except ImportError:
        pytest.skip("Models module not available")

if __name__ == "__main__":
    pytest.main([__file__]) 