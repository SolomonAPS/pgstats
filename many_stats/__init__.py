"""
Many-Stats: A comprehensive statistical genomics toolkit.

Built on top of sgkit with custom statistical functions implemented using Numba.
"""

__version__ = "0.1.0"
__author__ = "Your Name"
__email__ = "your.email@example.com"

from . import stats
from . import utils
from . import io

__all__ = ["stats", "utils", "io"]
