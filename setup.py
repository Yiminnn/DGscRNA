#!/usr/bin/env python3
"""
Setup script for scanpy-dgscrna package.
"""

from setuptools import setup, find_packages

setup(
    name="scanpy-dgscrna",
    version="0.1.0",
    description="A scanpy extension for automated cell type annotation using density scoring and deep learning",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="DGscRNA Team",
    author_email="your.email@example.com",
    url="https://github.com/yourusername/scanpy-dgscrna",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=[
        "scanpy>=1.9.0",
        "pandas>=1.3.0",
        "numpy>=1.20.0",
        "torch>=1.10.0",
        "torchmetrics>=0.6.0",
        "scikit-learn>=1.0.0",
        "matplotlib>=3.5.0",
        "seaborn>=0.11.0",
        "openpyxl>=3.0.0",
        "numba>=0.56.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov",
            "black",
            "isort",
            "flake8",
        ],
        "plotting": [
            "plotly>=5.0.0",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="single-cell genomics scanpy annotation deep-learning",
    project_urls={
        "Bug Reports": "https://github.com/yourusername/scanpy-dgscrna/issues",
        "Source": "https://github.com/yourusername/scanpy-dgscrna",
        "Documentation": "https://scanpy-dgscrna.readthedocs.io",
    },
)
