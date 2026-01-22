"""Domain distance computation for spatial transcriptomics.

This module provides functionality to compute spatial distances between
domains created by make_spatial_domains(). Supports multiple distance
metrics and output modes.

Key Features:
    - KD-Tree optimized O(n log n) minimum distance computation
    - Per-cell distance annotations
    - Domain-to-domain distance matrix
    - Multiple distance metrics (minimum, centroid, mean)

Examples
--------
>>> import scanpy as sc
>>> from spatialcore.spatial import calculate_domain_distances
>>> adata = sc.read_h5ad("domains.h5ad")
>>> # Calculate B cell to tumor distances
>>> adata = calculate_domain_distances(
...     adata,
...     source_domain_column="bcell_domain",
...     target_domain_column="tumor_domain",
...     distance_metric="minimum",
...     output_mode="both",
... )
>>> # Access per-cell distances
>>> print(adata.obs["distance_to_target"].describe())
>>> # Access distance matrix
>>> print(adata.uns["domain_distances"]["distance_matrix"])
"""

from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd
import anndata as ad
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist

from spatialcore.core.logging import get_logger
from spatialcore.core.metadata import update_metadata

logger = get_logger(__name__)


def calculate_domain_distances(
    adata: ad.AnnData,
    source_domain_column: str,
    target_domain_column: str,
    source_domain_subset: Optional[List[str]] = None,
    target_domain_subset: Optional[List[str]] = None,
    distance_metric: str = "minimum",
    output_mode: str = "both",
    output_distance_column: str = "distance_to_target",
    output_nearest_column: str = "nearest_target_domain",
    copy: bool = False,
) -> ad.AnnData:
    """
    Calculate spatial distances between domains.

    Computes distances from source domains to target domains using
    spatial coordinates. Supports per-cell distance annotation and
    domain-to-domain distance matrices.

    Parameters
    ----------
    adata
        AnnData with spatial coordinates in adata.obsm['spatial'] and
        domain labels in adata.obs.
    source_domain_column
        Column name containing source domain labels.
    target_domain_column
        Column name containing target domain labels.
        Can be the same as source_domain_column.
    source_domain_subset
        Specific source domains to measure from. If None, uses all
        non-null domains in source_domain_column.
    target_domain_subset
        Specific target domains to measure to. If None, uses all
        non-null domains in target_domain_column.
    distance_metric
        Distance metric to use:
        - "minimum": Shortest cell-to-cell distance (KD-Tree optimized)
        - "centroid": Center-to-center distance
        - "mean": Average of all pairwise distances
    output_mode
        What to output:
        - "cell": Per-cell distance columns only
        - "matrix": Distance matrix in adata.uns only
        - "both": Both per-cell and matrix
    output_distance_column
        Column name for per-cell distances (default 'distance_to_target').
        Only used if output_mode includes 'cell'.
    output_nearest_column
        Column name for nearest target domain (default 'nearest_target_domain').
        Only used if output_mode includes 'cell'.
    copy
        If True, operate on a copy of adata.

    Returns
    -------
    ad.AnnData
        AnnData with distance annotations:
        - adata.obs[output_distance_column]: distance to nearest target (if cell mode)
        - adata.obs[output_nearest_column]: name of nearest target domain (if cell mode)
        - adata.uns['domain_distances']: distance matrix and metadata (if matrix mode)

    Raises
    ------
    ValueError
        If columns not found or invalid parameters.

    Notes
    -----
    For minimum distance with per-cell annotation, uses scipy.spatial.cKDTree
    for O(n log n) performance instead of O(n^2) pairwise distances.

    Examples
    --------
    >>> import scanpy as sc
    >>> from spatialcore.spatial import calculate_domain_distances
    >>> adata = sc.read_h5ad("domains.h5ad")
    >>> # Distance from B cell domains to tumor domains
    >>> adata = calculate_domain_distances(
    ...     adata,
    ...     source_domain_column="bcell_domain",
    ...     target_domain_column="tumor_domain",
    ...     source_domain_subset=["Bcell_1", "Bcell_2"],
    ...     distance_metric="minimum",
    ...     output_mode="both",
    ... )
    >>> # Check per-cell results
    >>> bcell_mask = adata.obs["bcell_domain"].notna()
    >>> print(adata.obs.loc[bcell_mask, "distance_to_target"].describe())
    >>> # Access distance matrix
    >>> matrix = adata.uns["domain_distances"]["distance_matrix"]
    >>> print(f"Bcell_1 to Tumor_1: {matrix['Bcell_1']['Tumor_1']}")

    See Also
    --------
    make_spatial_domains : Create spatial domains.
    """
    # Validate inputs
    if "spatial" not in adata.obsm:
        raise ValueError(
            "adata.obsm['spatial'] not found. "
            f"Available keys: {list(adata.obsm.keys())}"
        )

    if source_domain_column not in adata.obs.columns:
        raise ValueError(
            f"Source column '{source_domain_column}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    if target_domain_column not in adata.obs.columns:
        raise ValueError(
            f"Target column '{target_domain_column}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    if distance_metric not in ["minimum", "centroid", "mean"]:
        raise ValueError(
            f"Invalid distance_metric: '{distance_metric}'. "
            "Must be 'minimum', 'centroid', or 'mean'."
        )

    if output_mode not in ["cell", "matrix", "both"]:
        raise ValueError(
            f"Invalid output_mode: '{output_mode}'. "
            "Must be 'cell', 'matrix', or 'both'."
        )

    # Handle copy
    adata = adata.copy() if copy else adata

    logger.info(
        f"Calculating domain distances: {source_domain_column} → {target_domain_column} "
        f"(metric={distance_metric}, mode={output_mode})"
    )

    # Get domain lists
    source_domains = adata.obs[source_domain_column].dropna().unique().tolist()
    target_domains = adata.obs[target_domain_column].dropna().unique().tolist()

    if source_domain_subset:
        source_domains = [d for d in source_domains if d in source_domain_subset]
    if target_domain_subset:
        target_domains = [d for d in target_domains if d in target_domain_subset]

    if not source_domains:
        raise ValueError(f"No valid source domains found in '{source_domain_column}'")
    if not target_domains:
        raise ValueError(f"No valid target domains found in '{target_domain_column}'")

    logger.debug(f"Source domains ({len(source_domains)}): {source_domains[:5]}...")
    logger.debug(f"Target domains ({len(target_domains)}): {target_domains[:5]}...")

    # Initialize results
    distance_matrix = pd.DataFrame(
        index=source_domains,
        columns=target_domains,
        dtype=float,
    )

    spatial = adata.obsm["spatial"]

    # Initialize per-cell columns if needed
    if output_mode in ["cell", "both"]:
        adata.obs[output_distance_column] = np.nan
        adata.obs[output_nearest_column] = None

    # Calculate distances based on metric
    if distance_metric == "minimum" and output_mode in ["cell", "both"]:
        # KD-Tree optimized approach for per-cell minimum distances
        logger.debug("Using KD-Tree optimization for minimum distances")

        # Build KD-tree for all target cells
        target_mask = adata.obs[target_domain_column].isin(target_domains)
        target_indices = np.where(target_mask.values)[0]
        target_coords = spatial[target_indices]
        target_tree = cKDTree(target_coords)
        target_domains_arr = adata.obs[target_domain_column].iloc[target_indices].values

        # Query for all source cells
        source_mask = adata.obs[source_domain_column].isin(source_domains)
        source_indices = np.where(source_mask.values)[0]
        source_coords = spatial[source_indices]

        if len(source_coords) > 0 and len(target_coords) > 0:
            # Single vectorized query - O(n log n)
            distances, nearest_idx = target_tree.query(source_coords, k=1)
            nearest_domains = target_domains_arr[nearest_idx]

            # Assign to adata.obs
            dist_col_idx = adata.obs.columns.get_loc(output_distance_column)
            nearest_col_idx = adata.obs.columns.get_loc(output_nearest_column)

            adata.obs.iloc[source_indices, dist_col_idx] = distances
            adata.obs.iloc[source_indices, nearest_col_idx] = nearest_domains

            # Build distance matrix from per-cell results
            source_domains_arr = adata.obs[source_domain_column].iloc[source_indices].values

            for src in source_domains:
                src_mask = source_domains_arr == src
                if not src_mask.any():
                    continue

                for tgt in target_domains:
                    # Self-distance check
                    if src == tgt and source_domain_column == target_domain_column:
                        distance_matrix.loc[src, tgt] = 0.0
                        continue

                    # Find minimum distance from src cells to tgt
                    src_distances = distances[src_mask]
                    src_nearest = nearest_domains[src_mask]
                    tgt_mask_local = src_nearest == tgt

                    if tgt_mask_local.any():
                        distance_matrix.loc[src, tgt] = src_distances[tgt_mask_local].min()
                    else:
                        # No direct minimum to this target, compute explicitly
                        src_cells = source_coords[src_mask]
                        tgt_cell_mask = target_domains_arr == tgt
                        if tgt_cell_mask.any():
                            tgt_cells = target_coords[tgt_cell_mask]
                            pairwise = cdist(src_cells, tgt_cells)
                            distance_matrix.loc[src, tgt] = pairwise.min()

    elif distance_metric == "centroid":
        # Centroid-to-centroid distances
        logger.debug("Computing centroid distances")

        source_centroids = {}
        target_centroids = {}

        for src in source_domains:
            mask = adata.obs[source_domain_column] == src
            coords = spatial[mask.values]
            if len(coords) > 0:
                source_centroids[src] = coords.mean(axis=0)

        for tgt in target_domains:
            mask = adata.obs[target_domain_column] == tgt
            coords = spatial[mask.values]
            if len(coords) > 0:
                target_centroids[tgt] = coords.mean(axis=0)

        for src in source_domains:
            if src not in source_centroids:
                continue
            for tgt in target_domains:
                if src == tgt and source_domain_column == target_domain_column:
                    distance_matrix.loc[src, tgt] = 0.0
                    continue
                if tgt not in target_centroids:
                    continue
                dist = np.linalg.norm(source_centroids[src] - target_centroids[tgt])
                distance_matrix.loc[src, tgt] = dist

        # Per-cell annotation for centroid metric
        if output_mode in ["cell", "both"]:
            source_mask = adata.obs[source_domain_column].isin(source_domains)
            for i, (idx, row) in enumerate(adata.obs[source_mask].iterrows()):
                src = row[source_domain_column]
                if src not in source_centroids:
                    continue

                min_dist = np.inf
                nearest = None
                for tgt, centroid in target_centroids.items():
                    if src == tgt and source_domain_column == target_domain_column:
                        continue
                    dist = np.linalg.norm(
                        spatial[adata.obs.index.get_loc(idx)] - centroid
                    )
                    if dist < min_dist:
                        min_dist = dist
                        nearest = tgt

                adata.obs.at[idx, output_distance_column] = min_dist
                adata.obs.at[idx, output_nearest_column] = nearest

    elif distance_metric == "mean":
        # Mean pairwise distances
        logger.debug("Computing mean pairwise distances")

        for src in source_domains:
            src_mask = adata.obs[source_domain_column] == src
            src_coords = spatial[src_mask.values]

            if len(src_coords) == 0:
                continue

            for tgt in target_domains:
                if src == tgt and source_domain_column == target_domain_column:
                    distance_matrix.loc[src, tgt] = 0.0
                    continue

                tgt_mask = adata.obs[target_domain_column] == tgt
                tgt_coords = spatial[tgt_mask.values]

                if len(tgt_coords) == 0:
                    continue

                # Compute all pairwise distances
                pairwise = cdist(src_coords, tgt_coords)
                distance_matrix.loc[src, tgt] = pairwise.mean()

        # Per-cell annotation for mean metric (use minimum for per-cell)
        if output_mode in ["cell", "both"]:
            # Fall back to minimum for per-cell (mean doesn't make sense per-cell)
            logger.debug("Using minimum distance for per-cell annotation with mean metric")
            target_mask = adata.obs[target_domain_column].isin(target_domains)
            target_indices = np.where(target_mask.values)[0]
            target_coords = spatial[target_indices]
            target_tree = cKDTree(target_coords)
            target_domains_arr = adata.obs[target_domain_column].iloc[target_indices].values

            source_mask = adata.obs[source_domain_column].isin(source_domains)
            source_indices = np.where(source_mask.values)[0]
            source_coords = spatial[source_indices]

            if len(source_coords) > 0 and len(target_coords) > 0:
                distances, nearest_idx = target_tree.query(source_coords, k=1)
                nearest_domains = target_domains_arr[nearest_idx]

                dist_col_idx = adata.obs.columns.get_loc(output_distance_column)
                nearest_col_idx = adata.obs.columns.get_loc(output_nearest_column)
                adata.obs.iloc[source_indices, dist_col_idx] = distances
                adata.obs.iloc[source_indices, nearest_col_idx] = nearest_domains

    else:
        # Minimum metric without per-cell (matrix only)
        logger.debug("Computing minimum distances (matrix only)")

        for src in source_domains:
            src_mask = adata.obs[source_domain_column] == src
            src_coords = spatial[src_mask.values]

            if len(src_coords) == 0:
                continue

            for tgt in target_domains:
                if src == tgt and source_domain_column == target_domain_column:
                    distance_matrix.loc[src, tgt] = 0.0
                    continue

                tgt_mask = adata.obs[target_domain_column] == tgt
                tgt_coords = spatial[tgt_mask.values]

                if len(tgt_coords) == 0:
                    continue

                pairwise = cdist(src_coords, tgt_coords)
                distance_matrix.loc[src, tgt] = pairwise.min()

    # Summary statistics
    valid_distances = distance_matrix.values[~np.isnan(distance_matrix.values)]
    summary = {
        "min_distance": float(valid_distances.min()) if len(valid_distances) > 0 else None,
        "max_distance": float(valid_distances.max()) if len(valid_distances) > 0 else None,
        "mean_distance": float(valid_distances.mean()) if len(valid_distances) > 0 else None,
        "median_distance": float(np.median(valid_distances)) if len(valid_distances) > 0 else None,
    }

    logger.info(
        f"Distance statistics: min={summary['min_distance']:.1f}, "
        f"max={summary['max_distance']:.1f}, mean={summary['mean_distance']:.1f}"
    )

    # Store matrix in adata.uns
    if output_mode in ["matrix", "both"]:
        result_metadata = {
            "source_domain_column": source_domain_column,
            "target_domain_column": target_domain_column,
            "distance_metric": distance_metric,
            "source_domains": source_domains,
            "target_domains": target_domains,
            "summary_statistics": summary,
            "distance_matrix": distance_matrix.to_dict(orient="index"),
        }
        adata.uns["domain_distances"] = result_metadata

    # Update metadata
    outputs = {"summary_statistics": summary}
    if output_mode in ["cell", "both"]:
        outputs["obs_distance"] = output_distance_column
        outputs["obs_nearest"] = output_nearest_column
    if output_mode in ["matrix", "both"]:
        outputs["uns"] = "domain_distances"

    update_metadata(
        adata,
        function_name="calculate_domain_distances",
        parameters={
            "source_domain_column": source_domain_column,
            "target_domain_column": target_domain_column,
            "source_domain_subset": source_domain_subset,
            "target_domain_subset": target_domain_subset,
            "distance_metric": distance_metric,
            "output_mode": output_mode,
        },
        outputs=outputs,
    )

    return adata


def get_distance_matrix(
    adata: ad.AnnData,
    key: str = "domain_distances",
) -> pd.DataFrame:
    """
    Extract domain distance matrix as a pandas DataFrame.

    Parameters
    ----------
    adata
        AnnData with domain distances computed.
    key
        Key in adata.uns containing distance data.

    Returns
    -------
    pd.DataFrame
        Distance matrix with source domains as rows and target domains as columns.

    Raises
    ------
    KeyError
        If distance matrix not found in adata.uns.

    Examples
    --------
    >>> from spatialcore.spatial import get_distance_matrix
    >>> matrix = get_distance_matrix(adata)
    >>> print(matrix)
              Tumor_1  Tumor_2
    Bcell_1     123.4    234.5
    Bcell_2     345.6    456.7
    """
    if key not in adata.uns:
        raise KeyError(
            f"'{key}' not found in adata.uns. "
            "Run calculate_domain_distances() with output_mode='matrix' or 'both' first."
        )

    data = adata.uns[key]
    if "distance_matrix" not in data:
        raise KeyError(f"'distance_matrix' not found in adata.uns['{key}']")

    return pd.DataFrame(data["distance_matrix"]).T
