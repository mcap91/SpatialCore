"""
SpatialCore: Standardized spatial statistics tools for computational biology.

A thin, robust wrapper around standard libraries to ensure Python and R users
get the exact same result for the same biological question.
"""

__version__ = "0.1.0"

from spatialcore import (
    annotation,
    clustering,
    core,
    diffusion,
    nmf,
    ontology,
    spatial,
)

__all__ = [
    "annotation",
    "clustering",
    "core",
    "diffusion",
    "nmf",
    "ontology",
    "spatial",
]
