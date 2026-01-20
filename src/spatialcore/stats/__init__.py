"""Statistical utilities for expression-based cell classification.

This module provides functions for classifying cell populations based on
expression thresholding (univariate or multivariate).

Functions
---------
classify_by_threshold
    Classify cells by expression thresholding (univariate or multivariate).
"""

from spatialcore.stats.classify import classify_by_threshold

__all__ = [
    "classify_by_threshold",
]
