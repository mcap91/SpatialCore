"""Neighborhood profile computation and niche identification.

This module provides functionality to compute cell-type composition profiles
for spatial neighborhoods and cluster them into discrete niche types.

Terminology (see docs_refs/vocab.md):
    - Neighborhood: Local spatial vicinity (k-NN or radius) - the INPUT
    - Niche: Compositional archetype from clustering neighborhoods - location-independent
    - Workflow: Neighborhood -> Niche -> Domain (local -> what kind -> where)

Key Features:
    - KD-Tree optimized k-NN and radius-based neighbor queries
    - Cell-type composition profiles (normalized or raw counts)
    - K-Means++ clustering for niche identification
    - MiniBatchKMeans for large datasets

Examples
--------
>>> import scanpy as sc
>>> from spatialcore.spatial import compute_neighborhood_profile, identify_niches
>>> adata = sc.read_h5ad("annotated.h5ad")
>>> # Step 1: Compute neighborhood profiles
>>> adata = compute_neighborhood_profile(
...     adata,
...     celltype_column="cell_type",
...     method="knn",
...     k=15,
... )
>>> # Step 2: Identify niches
>>> adata = identify_niches(adata, n_niches=8)
>>> print(adata.obs["niche"].value_counts())
"""

from typing import List, Literal, Optional

import numpy as np
import pandas as pd
import anndata as ad
from scipy.spatial import cKDTree
from sklearn.cluster import KMeans, MiniBatchKMeans

from spatialcore.core.logging import get_logger
from spatialcore.core.metadata import update_metadata

logger = get_logger(__name__)


def compute_neighborhood_profile(
    adata: ad.AnnData,
    celltype_column: str,
    method: Literal["knn", "radius"] = "knn",
    k: int = 15,
    radius: Optional[float] = None,
    normalize: bool = True,
    spatial_key: str = "spatial",
    key_added: str = "neighborhood_profile",
    copy: bool = False,
) -> ad.AnnData:
    """
    Compute cell-type composition profiles for each cell's spatial neighborhood.

    For each cell, identifies its spatial neighbors (via k-NN or fixed radius)
    and computes the cell-type composition of that neighborhood as a vector.
    This is the first step in the Neighborhood -> Niche workflow.

    Parameters
    ----------
    adata
        AnnData object with spatial coordinates in adata.obsm[spatial_key].
    celltype_column
        Column in adata.obs containing cell type labels.
        Missing labels are not allowed; fill or filter them first.
    method
        Method to define neighborhood:

        - "knn": k-nearest neighbors (default)
        - "radius": all cells within fixed radius
    k
        Number of neighbors for knn method (excluding the center cell).
        Default: 15. Ignored if method="radius".
    radius
        Distance threshold for radius method. Required if method="radius".
        Units should match spatial coordinates (microns for Xenium, pixels for CosMx).
    normalize
        If True, normalize profiles to proportions (sum to 1).
        If False, return raw counts. Default: True.

        - normalize=True: [0.6, 0.4] (proportions)
        - normalize=False: [6, 4] (raw counts)
    spatial_key
        Key in adata.obsm containing spatial coordinates. Default: "spatial".
    key_added
        Key to store the neighborhood profile matrix in adata.obsm.
        Default: "neighborhood_profile".
    copy
        If True, operate on a copy of adata. Default: False.

    Returns
    -------
    AnnData
        AnnData with neighborhood profiles in adata.obsm[key_added].
        Matrix shape: (n_cells, n_celltypes).
        Column order stored in adata.uns[f"{key_added}_celltypes"].

    Raises
    ------
    ValueError
        If celltype_column not found, spatial coordinates missing,
        missing labels, empty neighborhoods, or invalid parameters.

    Notes
    -----
    Uses scipy.spatial.cKDTree for efficient O(n log n) neighbor queries.
    For radius method, cells at tissue edges may have fewer neighbors.
    If any cell has no neighbors, this function raises a ValueError to
    avoid silent fallbacks. Use a larger radius, switch to knn, or
    pre-filter isolated cells before profiling.

    Examples
    --------
    >>> import scanpy as sc
    >>> from spatialcore.spatial import compute_neighborhood_profile
    >>> adata = sc.read_h5ad("annotated.h5ad")
    >>> # k-NN method (default)
    >>> adata = compute_neighborhood_profile(
    ...     adata,
    ...     celltype_column="cell_type",
    ...     method="knn",
    ...     k=15,
    ... )
    >>> print(adata.obsm["neighborhood_profile"].shape)
    (10000, 12)  # 10000 cells, 12 cell types
    >>> # Radius method
    >>> adata = compute_neighborhood_profile(
    ...     adata,
    ...     celltype_column="cell_type",
    ...     method="radius",
    ...     radius=50.0,  # 50 microns
    ... )

    See Also
    --------
    identify_niches : Cluster neighborhood profiles into niches.
    """
    # Input validation
    if spatial_key not in adata.obsm:
        raise ValueError(
            f"adata.obsm['{spatial_key}'] not found. "
            "Spatial coordinates are required for neighborhood computation."
        )

    if celltype_column not in adata.obs.columns:
        raise ValueError(
            f"Column '{celltype_column}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)[:10]}..."
        )

    if method not in ["knn", "radius"]:
        raise ValueError(
            f"Invalid method: '{method}'. Must be 'knn' or 'radius'."
        )

    n_cells = adata.n_obs

    if method == "knn" and k < 1:
        raise ValueError(f"k must be >= 1, got {k}")
    if method == "knn" and k >= n_cells:
        raise ValueError(
            f"k must be < number of cells ({n_cells}), got {k}"
        )

    if method == "radius":
        if radius is None:
            raise ValueError(
                "'radius' must be provided when method='radius'."
            )
        if radius <= 0:
            raise ValueError(f"radius must be > 0, got {radius}")

    # Handle copy
    adata = adata.copy() if copy else adata

    # Get spatial coordinates
    spatial_coords = adata.obsm[spatial_key]

    # Get cell types and create index mapping
    celltype_series = adata.obs[celltype_column]
    if celltype_series.isna().any():
        n_missing = int(celltype_series.isna().sum())
        raise ValueError(
            f"{n_missing} cells have missing labels in '{celltype_column}'. "
            "Fill or remove missing labels before computing neighborhoods."
        )

    cell_types = celltype_series.values
    unique_celltypes = sorted(celltype_series.unique())
    n_celltypes = len(unique_celltypes)
    celltype_to_idx = {ct: i for i, ct in enumerate(unique_celltypes)}

    if n_celltypes < 2:
        raise ValueError(
            f"At least 2 unique cell types required, found {n_celltypes}. "
            f"Check column '{celltype_column}'."
        )

    logger.info(
        f"Computing neighborhood profiles: {n_cells:,} cells, "
        f"{n_celltypes} cell types, method={method}"
    )

    # Build spatial index
    logger.debug("Building KD-Tree spatial index")
    tree = cKDTree(spatial_coords)

    # Initialize profile matrix
    neighborhood_profile = np.zeros((n_cells, n_celltypes), dtype=np.float32)

    # Query neighbors based on method
    if method == "knn":
        logger.debug(f"Querying {k} nearest neighbors per cell")
        # Query k+1 to include self, then remove
        query_k = k + 1
        _, neighbor_indices = tree.query(spatial_coords, k=query_k)

        # Count cell types for each cell (vectorized where possible)
        for i in range(n_cells):
            neighbors = neighbor_indices[i]
            neighbors = neighbors[neighbors != i]
            if neighbors.size != k:
                raise ValueError(
                    f"Expected {k} neighbors excluding self for cell {i}, "
                    f"got {neighbors.size}."
                )
            neighbor_types = cell_types[neighbors]

            for ct in neighbor_types:
                neighborhood_profile[i, celltype_to_idx[ct]] += 1

    else:  # radius method
        logger.debug(f"Querying neighbors within radius={radius}")
        neighbor_lists = tree.query_ball_point(spatial_coords, r=radius)

        for i, neighbors in enumerate(neighbor_lists):
            neighbors = [n for n in neighbors if n != i]
            if len(neighbors) == 0:
                continue

            neighbor_types = cell_types[neighbors]

            for ct in neighbor_types:
                neighborhood_profile[i, celltype_to_idx[ct]] += 1

    row_sums = neighborhood_profile.sum(axis=1)
    n_empty = int((row_sums == 0).sum())
    if n_empty > 0:
        raise ValueError(
            f"{n_empty} cells have empty neighborhood profiles. "
            "Increase radius, switch to knn, or pre-filter isolated cells "
            "before profiling."
        )

    # Normalize to proportions if requested
    if normalize:
        neighborhood_profile = neighborhood_profile / row_sums[:, None]
        logger.debug("Normalized profiles to proportions")

    # Store results
    adata.obsm[key_added] = neighborhood_profile
    adata.uns[f"{key_added}_celltypes"] = list(unique_celltypes)

    logger.info(
        f"Stored neighborhood profiles in adata.obsm['{key_added}'] "
        f"(shape: {neighborhood_profile.shape})"
    )

    # Update metadata
    update_metadata(
        adata,
        function_name="compute_neighborhood_profile",
        parameters={
            "celltype_column": celltype_column,
            "method": method,
            "k": k if method == "knn" else None,
            "radius": radius if method == "radius" else None,
            "normalize": normalize,
            "spatial_key": spatial_key,
        },
        outputs={
            "obsm": key_added,
            "uns": f"{key_added}_celltypes",
            "n_celltypes": n_celltypes,
            "n_cells": n_cells,
        },
    )

    return adata


def identify_niches(
    adata: ad.AnnData,
    n_niches: int,
    method: Literal["kmeans", "minibatch_kmeans"] = "kmeans",
    neighborhood_key: str = "neighborhood_profile",
    key_added: str = "niche",
    random_state: int = 0,
    n_init: int = 10,
    max_iter: int = 300,
    copy: bool = False,
) -> ad.AnnData:
    """
    Cluster neighborhood profiles into discrete niche types.

    Takes precomputed neighborhood profiles (cell-type compositions) and
    clusters them to identify recurring microenvironment archetypes (niches).
    Niches are location-independent - the same niche can appear in multiple
    non-contiguous regions of the tissue.

    This is the second step in the Neighborhood -> Niche workflow.

    Parameters
    ----------
    adata
        AnnData object with neighborhood profiles in adata.obsm[neighborhood_key].
        Run compute_neighborhood_profile() first.
    n_niches
        Number of niche clusters to identify.
    method
        Clustering method:

        - "kmeans": Standard K-Means with kmeans++ initialization (default)
        - "minibatch_kmeans": Mini-batch K-Means for large datasets (>100k cells)
    neighborhood_key
        Key in adata.obsm containing neighborhood profiles.
        Default: "neighborhood_profile".
    key_added
        Column name for niche labels in adata.obs. Default: "niche".
    random_state
        Random seed for reproducibility. Default: 0.
    n_init
        Number of K-Means initializations (runs with different centroid seeds).
        Default: 10.
    max_iter
        Maximum iterations per K-Means run. Default: 300.
    copy
        If True, operate on a copy of adata. Default: False.

    Returns
    -------
    AnnData
        AnnData with:

        - adata.obs[key_added]: Niche labels (categorical)
        - adata.uns['niche_centroids']: (n_niches, n_celltypes) cluster centers
        - adata.uns['niche_params']: Clustering parameters and metadata

    Raises
    ------
    ValueError
        If neighborhood profiles not found, contain empty profiles,
        or invalid parameters.

    Notes
    -----
    Uses kmeans++ initialization for smarter centroid selection, avoiding
    poor local minima. For datasets >100k cells, consider using
    method="minibatch_kmeans" for faster clustering.
    Empty neighborhood profiles are not allowed; fix neighborhood
    computation before clustering.

    The number of niches is a hyperparameter that may require tuning.
    Consider using silhouette scores or biological interpretation to
    select the optimal number.

    Examples
    --------
    >>> import scanpy as sc
    >>> from spatialcore.spatial import compute_neighborhood_profile, identify_niches
    >>> adata = sc.read_h5ad("annotated.h5ad")
    >>> adata = compute_neighborhood_profile(adata, celltype_column="cell_type")
    >>> adata = identify_niches(adata, n_niches=8)
    >>> print(adata.obs["niche"].value_counts())
    niche_1    2500
    niche_2    1800
    ...
    >>> # View niche composition (centroids)
    >>> import pandas as pd
    >>> celltypes = adata.uns["neighborhood_profile_celltypes"]
    >>> centroids = pd.DataFrame(
    ...     adata.uns["niche_centroids"],
    ...     index=[f"niche_{i+1}" for i in range(8)],
    ...     columns=celltypes,
    ... )
    >>> print(centroids)

    See Also
    --------
    compute_neighborhood_profile : Compute neighborhood profiles (run first).
    make_spatial_domains : Create spatially contiguous domains.
    """
    # Input validation
    if neighborhood_key not in adata.obsm:
        raise ValueError(
            f"adata.obsm['{neighborhood_key}'] not found. "
            "Run compute_neighborhood_profile() first."
        )

    if method not in ["kmeans", "minibatch_kmeans"]:
        raise ValueError(
            f"Invalid method: '{method}'. Must be 'kmeans' or 'minibatch_kmeans'."
        )

    n_cells = adata.n_obs
    if n_niches < 2:
        raise ValueError(f"n_niches must be >= 2, got {n_niches}")
    if n_niches > n_cells:
        raise ValueError(
            f"n_niches ({n_niches}) cannot exceed number of cells ({n_cells})"
        )

    # Handle copy
    adata = adata.copy() if copy else adata

    # Get neighborhood profiles
    profiles = adata.obsm[neighborhood_key]
    logger.info(
        f"Identifying {n_niches} niches from {n_cells:,} cells "
        f"(method={method}, random_state={random_state})"
    )

    # Check for cells with empty profiles (all zeros)
    empty_mask = profiles.sum(axis=1) == 0
    n_empty = int(empty_mask.sum())
    if n_empty > 0:
        raise ValueError(
            f"{n_empty} cells have empty neighborhood profiles. "
            "Increase radius, switch to knn, or pre-filter isolated cells "
            "before profiling."
        )

    # Run clustering with kmeans++ initialization
    if method == "kmeans":
        logger.debug(f"Running KMeans (n_init={n_init}, max_iter={max_iter})")
        clusterer = KMeans(
            n_clusters=n_niches,
            init="k-means++",
            n_init=n_init,
            max_iter=max_iter,
            random_state=random_state,
        )
    else:  # minibatch_kmeans
        batch_size = min(1024, n_cells)
        logger.debug(
            f"Running MiniBatchKMeans (batch_size={batch_size}, "
            f"n_init={n_init}, max_iter={max_iter})"
        )
        clusterer = MiniBatchKMeans(
            n_clusters=n_niches,
            init="k-means++",
            n_init=n_init,
            max_iter=max_iter,
            random_state=random_state,
            batch_size=batch_size,
        )

    # Fit and predict
    labels = clusterer.fit_predict(profiles)
    centroids = clusterer.cluster_centers_
    inertia = clusterer.inertia_

    # Create categorical labels (1-indexed for user-friendliness)
    niche_names = [f"niche_{i + 1}" for i in range(n_niches)]
    niche_labels = pd.Categorical(
        [f"niche_{label + 1}" for label in labels],
        categories=niche_names,
    )

    # Store results
    adata.obs[key_added] = niche_labels

    adata.uns["niche_centroids"] = centroids
    adata.uns["niche_params"] = {
        "n_niches": n_niches,
        "method": method,
        "neighborhood_key": neighborhood_key,
        "random_state": random_state,
        "n_init": n_init,
        "max_iter": max_iter,
        "inertia": float(inertia),
    }

    # Log cluster sizes
    cluster_sizes = pd.Series(labels).value_counts().sort_index()
    logger.info(
        f"Niche sizes: min={cluster_sizes.min()}, "
        f"max={cluster_sizes.max()}, mean={cluster_sizes.mean():.0f}"
    )
    logger.info(
        f"Stored niche labels in adata.obs['{key_added}'] "
        f"and centroids in adata.uns['niche_centroids']"
    )

    # Update metadata
    update_metadata(
        adata,
        function_name="identify_niches",
        parameters={
            "n_niches": n_niches,
            "method": method,
            "neighborhood_key": neighborhood_key,
            "random_state": random_state,
            "n_init": n_init,
            "max_iter": max_iter,
        },
        outputs={
            "obs": key_added,
            "uns_centroids": "niche_centroids",
            "uns_params": "niche_params",
            "inertia": float(inertia),
        },
    )

    return adata
