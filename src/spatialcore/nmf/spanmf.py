"""Spatial non-negative matrix factorization (spaNMF).

This module implements spaNMF, which identifies gene expression programs
that capture spatial tissue organization. The method:

1. Computes neighborhood-averaged expression (spatial smoothing)
2. Runs NMF to identify spatial gene modules

Key advantages over standard NMF:
    - Spatial awareness: Captures local microenvironment signals
    - Noise reduction: Averaging reduces single-cell noise
    - Interpretable: Each component represents a spatial gene program

Examples
--------
>>> import scanpy as sc
>>> from spatialcore.nmf import calculate_neighbor_expression, run_spanmf
>>> adata = sc.read_h5ad("spatial_data.h5ad")
>>> adata = calculate_neighbor_expression(adata, k=10)
>>> adata = run_spanmf(adata, n_components=10)
>>> # Cell loadings in adata.obsm["X_spanmf"]
>>> # Gene loadings in adata.varm["spanmf_loadings"]
"""

from typing import List, Literal, Optional

import numpy as np
import anndata as ad
from scipy import sparse
from sklearn.decomposition import NMF

from spatialcore.core.logging import get_logger
from spatialcore.core.metadata import update_metadata
from spatialcore.spatial._neighbors import build_neighbor_graph

logger = get_logger(__name__)


def calculate_neighbor_expression(
    adata: ad.AnnData,
    method: Literal["knn", "radius"] = "knn",
    k: int = 10,
    radius: Optional[float] = None,
    layer: Optional[str] = None,
    scale: bool = True,
    include_self: bool = True,
    spatial_key: str = "spatial",
    key_added: str = "neighbor_expression",
    copy: bool = False,
) -> ad.AnnData:
    """
    Compute neighborhood-averaged gene expression.

    For each cell, computes the mean expression across its spatial
    neighbors using efficient sparse matrix multiplication. This
    creates a smoothed expression matrix that captures local
    microenvironment signals.

    Parameters
    ----------
    adata
        AnnData with spatial coordinates in obsm[spatial_key].
    method
        Method for defining neighborhoods:

        - "knn": k-nearest neighbors (default, recommended)
        - "radius": all cells within distance threshold
    k
        Number of neighbors for knn method (including self if
        include_self=True). Default: 10.
    radius
        Distance threshold for radius method. Required if method="radius".
        Units should match spatial coordinates.
    layer
        Expression layer to use. If None, uses adata.X. Default: None.
    scale
        If True, apply min-max scaling to each gene BEFORE averaging.
        Scales to [0, 1] range. Default: True.
    include_self
        If True, include the cell's own expression in the average.
        Default: True.

        - True: average of self + neighbors (recommended for spaNMF)
        - False: average of neighbors only
    spatial_key
        Key in adata.obsm for spatial coordinates. Default: "spatial".
    key_added
        Key to store result in adata.layers. Default: "neighbor_expression".
    copy
        If True, operate on a copy of adata. Default: False.

    Returns
    -------
    AnnData
        AnnData with neighborhood-averaged expression stored in
        adata.layers[key_added].

    Raises
    ------
    ValueError
        If spatial coordinates missing, layer not found, invalid parameters,
        or any cell has no neighbors.

    Notes
    -----
    Uses sparse matrix multiplication for efficiency:
    ``X_avg = W @ X`` where W is a row-normalized weight matrix.

    For 100k cells with k=10 neighbors, this completes in seconds
    rather than minutes (compared to row-wise loops).

    Examples
    --------
    >>> adata = calculate_neighbor_expression(adata, k=15)
    >>> # Result in adata.layers["neighbor_expression"]

    See Also
    --------
    run_spanmf : Run NMF on the neighbor-averaged expression.
    """
    # Input validation
    if spatial_key not in adata.obsm:
        raise ValueError(
            f"Spatial coordinates not found in adata.obsm['{spatial_key}']. "
            f"Available keys: {list(adata.obsm.keys())}"
        )

    if layer is not None and layer not in adata.layers:
        raise ValueError(
            f"Layer '{layer}' not found. Available: {list(adata.layers.keys())}"
        )

    if method not in ("knn", "radius"):
        raise ValueError(f"method must be 'knn' or 'radius', got '{method}'")

    if method == "radius" and radius is None:
        raise ValueError("radius is required when method='radius'")

    adata = adata.copy() if copy else adata

    n_cells = adata.n_obs
    n_genes = adata.n_vars

    logger.info(
        f"Computing neighbor expression: {n_cells:,} cells, {n_genes:,} genes, "
        f"method={method}, k={k if method == 'knn' else 'N/A'}, include_self={include_self}"
    )

    # Get spatial coordinates
    coords = adata.obsm[spatial_key]

    # Build neighbor graph with weight matrix
    _, weight_matrix = build_neighbor_graph(
        coords=coords,
        method=method,
        k=k,
        radius=radius,
        include_self=include_self,
    )

    # Get expression matrix
    if layer is None:
        X = adata.X
    else:
        X = adata.layers[layer]

    # Optional min-max scaling per gene BEFORE averaging
    if scale:
        logger.debug("Applying min-max scaling before averaging")
        if sparse.issparse(X):
            X = X.toarray()
        X_min = X.min(axis=0)
        X_max = X.max(axis=0)
        X_range = X_max - X_min
        X_range[X_range == 0] = 1.0
        X = (X - X_min) / X_range

    # Sparse matrix multiplication: W @ X
    # W is (n_cells, n_cells) sparse, X is (n_cells, n_genes)
    logger.debug("Computing neighbor-averaged expression (sparse @ matrix)")

    X_avg = weight_matrix @ X
    if sparse.issparse(X_avg) and scale:
        X_avg = X_avg.toarray()

    # Store result (keep sparse if unscaled and sparse to save memory)
    if sparse.issparse(X_avg):
        adata.layers[key_added] = X_avg
    else:
        adata.layers[key_added] = X_avg.astype(np.float32)

    logger.info(f"Stored neighbor expression in adata.layers['{key_added}']")

    # Update metadata
    update_metadata(
        adata,
        function_name="calculate_neighbor_expression",
        parameters={
            "method": method,
            "k": k if method == "knn" else None,
            "radius": radius if method == "radius" else None,
            "layer": layer,
            "scale": scale,
            "include_self": include_self,
            "spatial_key": spatial_key,
        },
        outputs={
            "layers": key_added,
            "n_cells": n_cells,
            "n_genes": n_genes,
        },
    )

    return adata


def run_spanmf(
    adata: ad.AnnData,
    n_components: int = 10,
    neighbor_key: Optional[str] = None,
    genes_to_use: Optional[List[str]] = None,
    genes_to_exclude: Optional[List[str]] = None,
    random_state: int = 42,
    max_iter: int = 100,
    key_added: str = "spanmf",
    copy: bool = False,
) -> ad.AnnData:
    """
    Run spatial NMF on neighborhood-averaged expression.

    Decomposes the neighbor-averaged expression matrix into two
    non-negative matrices: cell loadings (W) and gene loadings (H).
    Each component represents a spatial gene expression program.

    Parameters
    ----------
    adata
        AnnData with spatial coordinates. Must have
        adata.layers[neighbor_key] from calculate_neighbor_expression().
    n_components
        Number of NMF components (factors) to compute. Default: 10.
    neighbor_key
        Key in adata.layers with pre-computed neighbor expression.
        If None, defaults to "neighbor_expression".
    genes_to_use
        Subset of genes to include in NMF. If None, uses all genes.
    genes_to_exclude
        Genes to exclude from NMF (e.g., low variance or problematic genes).
    random_state
        Random seed for reproducibility. Default: 42.
    max_iter
        Maximum iterations for NMF optimization. Default: 200.
    key_added
        Prefix for stored results. Default: "spanmf".

        - obsm[f"X_{key_added}"]: Cell loadings
        - varm[f"{key_added}_loadings"]: Gene loadings
        - uns[key_added]: Parameters and metadata
    copy
        If True, operate on a copy of adata. Default: False.

    Returns
    -------
    AnnData
        AnnData with NMF results:

        - obsm["X_spanmf"]: Cell loadings (n_cells, n_components)
        - varm["spanmf_loadings"]: Gene loadings (n_genes, n_components)
        - uns["spanmf"]: Dict with parameters, genes_used, reconstruction_error

    Raises
    ------
    ValueError
        If neighbor_key not found, no genes remain after filtering,
        or n_components invalid.

    Notes
    -----
    Uses sklearn.decomposition.NMF with 'nndsvda' initialization,
    which is robust for sparse or near-sparse input.

    Gene loadings are stored for ALL genes in adata.var, with zeros
    for genes that were excluded from NMF.

    Examples
    --------
    >>> # Full pipeline
    >>> adata = calculate_neighbor_expression(adata, k=10)
    >>> adata = run_spanmf(adata, n_components=10)

    See Also
    --------
    calculate_neighbor_expression : Compute neighbor-averaged expression.
    """
    # Input validation
    if n_components < 2:
        raise ValueError(f"n_components must be >= 2, got {n_components}")

    if n_components >= adata.n_obs:
        raise ValueError(
            f"n_components ({n_components}) must be < n_cells ({adata.n_obs})"
        )

    if n_components >= adata.n_vars:
        raise ValueError(
            f"n_components ({n_components}) must be < n_genes ({adata.n_vars})"
        )

    adata = adata.copy() if copy else adata

    # Get pre-computed neighbor expression
    if neighbor_key is None:
        neighbor_key = "neighbor_expression"
    if neighbor_key not in adata.layers:
        raise ValueError(
            f"Layer '{neighbor_key}' not found. "
            "Run calculate_neighbor_expression() first."
        )

    X = adata.layers[neighbor_key]
    n_cells, n_genes_total = X.shape

    # Gene filtering
    gene_mask = np.ones(n_genes_total, dtype=bool)

    if genes_to_use is not None:
        genes_to_use_set = set(genes_to_use)
        gene_mask &= np.array([g in genes_to_use_set for g in adata.var_names])
        n_requested = len(genes_to_use)
        n_found = gene_mask.sum()
        if n_found < n_requested:
            missing = genes_to_use_set - set(adata.var_names[gene_mask])
            logger.warning(
                f"Requested {n_requested} genes, found {n_found}. "
                f"Missing: {list(missing)[:5]}..."
            )

    if genes_to_exclude is not None:
        genes_to_exclude_set = set(genes_to_exclude)
        gene_mask &= np.array(
            [g not in genes_to_exclude_set for g in adata.var_names]
        )

    n_genes = int(gene_mask.sum())
    if n_genes == 0:
        raise ValueError(
            "No genes remaining after filtering. "
            "Check genes_to_use and genes_to_exclude parameters."
        )

    if n_genes < n_components:
        raise ValueError(
            f"n_genes ({n_genes}) must be >= n_components ({n_components})"
        )

    gene_names = adata.var_names[gene_mask].tolist()
    X_filtered = X[:, gene_mask]

    logger.info(
        f"Running NMF: {n_cells:,} cells, {n_genes:,} genes, "
        f"{n_components} components"
    )

    # Convert sparse to dense if needed
    if sparse.issparse(X_filtered):
        n_cells, n_genes = X_filtered.shape
        est_bytes = n_cells * n_genes * 4
        if est_bytes > 2e9:
            raise ValueError(
                "Filtered matrix is too large to densify safely. "
                "Reduce genes_to_use, subset cells, or disable scaling upstream."
            )
        X_filtered = X_filtered.toarray()

    # Ensure non-negative
    if X_filtered.min() < -1e-6:
        raise ValueError(
            "Negative values found in input matrix. "
            "Check scaling and preprocessing; NMF requires non-negative data."
        )
    if X_filtered.min() < 0:
        n_negative = (X_filtered < 0).sum()
        logger.warning(
            f"Clipping {n_negative} small negative values to 0 "
            "(floating-point artifacts)"
        )
        X_filtered = np.clip(X_filtered, 0, None)

    # Run NMF
    nmf = NMF(
        n_components=n_components,
        init="nndsvda",
        random_state=random_state,
        max_iter=max_iter,
    )

    W = nmf.fit_transform(X_filtered)  # (n_cells, n_components)
    H = nmf.components_  # (n_components, n_genes)

    logger.info(
        f"NMF converged in {nmf.n_iter_} iterations, "
        f"reconstruction error: {nmf.reconstruction_err_:.4f}"
    )

    # Store cell loadings
    adata.obsm[f"X_{key_added}"] = W.astype(np.float32)

    # Store gene loadings (full size, zeros for filtered genes)
    loadings_full = np.zeros((n_genes_total, n_components), dtype=np.float32)
    loadings_full[gene_mask, :] = H.T
    adata.varm[f"{key_added}_loadings"] = loadings_full

    # Store metadata
    adata.uns[key_added] = {
        "n_components": n_components,
        "genes_used": gene_names,
        "n_genes": n_genes,
        "reconstruction_error": float(nmf.reconstruction_err_),
        "n_iter": nmf.n_iter_,
        "random_state": random_state,
        "max_iter": max_iter,
        "neighbor_key": neighbor_key,
    }

    logger.info(
        f"Stored results: obsm['X_{key_added}'], "
        f"varm['{key_added}_loadings'], uns['{key_added}']"
    )

    # Update metadata
    update_metadata(
        adata,
        function_name="run_spanmf",
        parameters={
            "n_components": n_components,
            "neighbor_key": neighbor_key,
            "genes_to_use": genes_to_use,
            "genes_to_exclude": genes_to_exclude,
            "random_state": random_state,
            "max_iter": max_iter,
        },
        outputs={
            "obsm": f"X_{key_added}",
            "varm": f"{key_added}_loadings",
            "uns": key_added,
            "reconstruction_error": float(nmf.reconstruction_err_),
        },
    )

    return adata
