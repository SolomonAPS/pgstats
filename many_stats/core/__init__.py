"""
Core data structures and analysis framework for population genetics statistics.

This module provides the main classes for managing genomic data and running
windowed and non-windowed population genetics statistics while accounting for
missing data and callable sites.
"""

# Import the main classes from dataset.py
from .dataset import GenomicDataset, WindowConfig, CallableSitesConfig

# Make these available when importing from many_stats.core
__all__ = ['GenomicDataset', 'WindowConfig', 'CallableSitesConfig']