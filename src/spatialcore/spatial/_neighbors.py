"""Internal utilities for spatial neighbor graph construction.

This module provides the core neighbor-finding logic used by multiple
spatial analysis functions. It is an internal module (_neighbors.py)
and should not be imported directly by users.

Used by:
    - spatial.neighborhoods.compute_neighborhood_profile()
    - nmf.spanmf.calculate_neighbor_expression()
"""

from typing import List, Literal, Optional, Tuple

import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix

from spatialcore.core.logging import get_logger

logger = get_logger(__name__)


def build_neighbor_graph(
    coords: np.ndarray,
    method: Literal["knn", "radius"] = "knn",
    k: int = 10,
    radius: Optional[float] = None,
    include_self: bool = False,
) -> Tuple[List[np.ndarray], csr_matrix]:
    """
    Build spatial neighbor graph with row-normalized weight matrix.

    Constructs a neighbor graph from spatial coordinates using either
    k-nearest neighbors or radius-based queries. Returns both the
    neighbor indices and a sparse weight matrix suitable for matrix
    multiplication operations.

    Parameters
    ----------
    coords
        Spatial coordinates array of shape (n_cells, 2) or (n_cells, 3).
    method
        Method for defining neighborhoods:

        - "knn": k-nearest neighbors (fixed number of neighbors)
        - "radius": all cells within distance threshold
    k
        Number of neighbors for knn method. If include_self=False,
        returns k neighbors excluding self. If include_self=True,
        returns k neighbors including self. Default: 10.
    radius
        Distance threshold for radius method. Required if method="radius".
        Units should match coordinate units (microns, pixels, etc.).
    include_self
        If True, include the cell itself in its neighbor list.
        Default: False.

        - False: neighbors are OTHER cells only (matches R knearneigh)
        - True: neighbors include self (explicit opt-in)

    Returns
    -------
    neighbor_indices
        List of length n_cells where neighbor_indices[i] is a numpy array
        of neighbor indices for cell i.
    weight_matrix
        Sparse CSR matrix of shape (n_cells, n_cells) where
        weight_matrix[i, j] = 1 / n_neighbors[i] if j is a neighbor of i,
        else 0. Each row sums to 1 (row-normalized).

    Raises
    ------
    ValueError
        If coords is not 2D array, method is invalid, radius not provided
        for radius method, k < 1, or any cell has no neighbors.

    Notes
    -----
    Uses scipy.spatial.cKDTree for O(n log n) neighbor queries.

    The weight matrix enables efficient sparse matrix operations:

    - Expression averaging: X_avg = weight_matrix @ X
    - Cell-type profiles: Use neighbor_indices directly for counting

    Examples
    --------
    >>> coords = adata.obsm["spatial"]
    >>> indices, W = build_neighbor_graph(coords, method="knn", k=10)
    >>> # For expression averaging (sparse @ dense):
    >>> X_smoothed = W @ adata.X
    >>> # For cell-type counting:
    >>> for i, neighbors in enumerate(indices):
    ...     types = cell_types[neighbors]
    """
    # Input validation
    if coords.ndim != 2:
        raise ValueError(f"coords must be 2D array, got shape {coords.shape}")

    n_cells = coords.shape[0]

    if method not in ("knn", "radius"):
        raise ValueError(f"method must be 'knn' or 'radius', got '{method}'")

    if method == "knn":
        if k < 1:
            raise ValueError(f"k must be >= 1, got {k}")
        if k >= n_cells:
            raise ValueError(f"k ({k}) must be < n_cells ({n_cells})")

    if method == "radius":
        if radius is None:
            raise ValueError("radius is required when method='radius'")
        if radius <= 0:
            raise ValueError(f"radius must be > 0, got {radius}")

    # Build KD-Tree
    logger.debug(f"Building KD-Tree for {n_cells:,} cells")
    tree = cKDTree(coords)

    # Query neighbors
    if method == "knn":
        # Query k+1 to include self, then optionally remove
        query_k = k + 1
        _, indices_array = tree.query(coords, k=query_k)
        indices_array = np.atleast_2d(indices_array)  # handle k=1 shape

        neighbor_indices = []
        for i in range(n_cells):
            neighbors = indices_array[i]
            if not include_self:
                neighbors = neighbors[neighbors != i]
            else:
                # Ensure exactly k neighbors including self
                # (self is always closest, so take first k)
                neighbors = neighbors[:k]
            neighbor_indices.append(neighbors)

    else:  # radius method
        logger.debug(f"Querying neighbors within radius={radius}")
        raw_indices = tree.query_ball_point(coords, r=radius)

        neighbor_indices = []
        for i, neighbors in enumerate(raw_indices):
            neighbors = np.array(neighbors, dtype=np.int64)
            if not include_self:
                neighbors = neighbors[neighbors != i]
            if k is not None and len(neighbors) > k:
                # Cap to k nearest within radius
                dists = np.linalg.norm(coords[neighbors] - coords[i], axis=1)
                neighbors = neighbors[np.argsort(dists)[:k]]
            neighbor_indices.append(neighbors)

    # Validate no empty neighborhoods
    empty_cells = [i for i, nb in enumerate(neighbor_indices) if len(nb) == 0]
    if empty_cells:
        n_empty = len(empty_cells)
        example_cells = empty_cells[:5]
        raise ValueError(
            f"{n_empty} cells have no neighbors (e.g., cells {example_cells}). "
            f"Increase {'k' if method == 'knn' else 'radius'} or filter isolated cells."
        )

    # Build sparse weight matrix (row-normalized)
    logger.debug("Building sparse weight matrix")
    rows = []
    cols = []
    data = []

    for i, neighbors in enumerate(neighbor_indices):
        n_neighbors = len(neighbors)
        weight = 1.0 / n_neighbors
        for j in neighbors:
            rows.append(i)
            cols.append(j)
            data.append(weight)

    weight_matrix = csr_matrix(
        (data, (rows, cols)),
        shape=(n_cells, n_cells),
        dtype=np.float32,
    )

    logger.debug(
        f"Built neighbor graph: {n_cells:,} cells, "
        f"~{np.mean([len(nb) for nb in neighbor_indices]):.1f} neighbors/cell"
    )

    return neighbor_indices, weight_matrix
