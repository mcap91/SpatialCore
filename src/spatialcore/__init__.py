"""
SpatialCore: Standardized spatial statistics tools for computational biology.

A thin, robust wrapper around standard libraries to ensure Python and R users
get the exact same result for the same biological question.
"""

__version__ = "0.1.0"

from spatialcore import (
    annotation,
    core,
    diffusion,
    nmf,
    plotting,
    spatial,
    stats,
)

__all__ = [
    "annotation",
    "core",
    "diffusion",
    "nmf",
    "plotting",
    "spatial",
    "stats",
]

# Try to import internal modules (only available on dev branch)
try:
    from spatialcore import internal
    __all__.append("internal")
except ImportError:
    internal = None  # Not available in published package
