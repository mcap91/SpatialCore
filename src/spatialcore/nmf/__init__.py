"""Spatial non-negative matrix factorization."""

from spatialcore.nmf.spanmf import (
    calculate_neighbor_expression,
    run_spanmf,
)

__all__ = [
    "calculate_neighbor_expression",
    "run_spanmf",
]
