"""Spatial autocorrelation statistics: Moran's I and Lee's L (Global and Local).

This module provides spatial autocorrelation analysis for identifying
spatially variable genes and spatial clustering patterns in expression data.

Key Features:
    - Sparse spatial weights matrix construction (O(nk) memory)
    - Vectorized gene batching for memory-efficient computation
    - Permutation-based significance testing
    - FDR correction for multiple testing
    - LISA quadrant classification (HH/LL/HL/LH/NS)

Statistics:
    - **Moran's I (univariate)**: Does gene X cluster spatially?
    - **Lee's L (bivariate)**: Do genes X and Y co-localize spatially?

References
----------
Moran, P.A.P. (1950). "Notes on Continuous Stochastic Phenomena".
    Biometrika. 37(1): 17-23.

Anselin, L. (1995). "Local Indicators of Spatial Association - LISA".
    Geographical Analysis. 27(2): 93-115.

Lee, S.I. (2001). "Developing a bivariate spatial association measure".
    Journal of Geographical Systems. 3(4): 369-385.

Examples
--------
>>> import scanpy as sc
>>> from spatialcore.spatial import morans_i, local_morans_i, lees_l, lees_l_local
>>> adata = sc.read_h5ad("xenium.h5ad")
>>> # Global Moran's I for spatially variable genes
>>> adata = morans_i(adata, genes=["CD8A", "FOXP3"], n_permutations=100)
>>> # Local Lee's L for gene co-localization
>>> adata = lees_l_local(adata, gene_pairs=("CD8A", "GZMB"))
>>> print(adata.obs["CD8A_GZMB_quadrant"])
"""

from itertools import combinations
from typing import List, Literal, Optional, Tuple, Union
import time

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.sparse import csr_matrix
from sklearn.neighbors import NearestNeighbors
import squidpy as sq

from spatialcore.core.logging import get_logger
from spatialcore.core.metadata import update_metadata

logger = get_logger(__name__)

# Quadrant encoding: 0=NS, 1=HH, 2=LL, 3=HL, 4=LH
QUADRANT_LABELS = {0: "NS", 1: "HH", 2: "LL", 3: "HL", 4: "LH"}


# =============================================================================
# SPARSE MATRIX UTILITIES
# =============================================================================


def _compute_mean_sparse(X: Union[np.ndarray, csr_matrix]) -> np.ndarray:
    """Compute mean per gene efficiently on sparse matrices.

    Parameters
    ----------
    X
        Expression matrix, shape (n_cells, n_genes). Can be sparse or dense.

    Returns
    -------
    np.ndarray
        Mean per gene, shape (n_genes,).
    """
    if sparse.issparse(X):
        return np.asarray(X.mean(axis=0)).flatten()
    return np.mean(X, axis=0)


def _compute_variance_sparse(X: Union[np.ndarray, csr_matrix]) -> np.ndarray:
    """Compute variance per gene without dense conversion.

    Uses the identity: var = E[X^2] - E[X]^2

    This avoids densifying sparse matrices, which would cause memory explosion
    for large datasets.

    Parameters
    ----------
    X
        Expression matrix, shape (n_cells, n_genes). Can be sparse or dense.

    Returns
    -------
    np.ndarray
        Variance per gene, shape (n_genes,).
    """
    if sparse.issparse(X):
        mean = np.asarray(X.mean(axis=0)).flatten()
        # X.power(2) keeps result sparse
        X_squared = X.power(2)
        sq_mean = np.asarray(X_squared.mean(axis=0)).flatten()
        return sq_mean - mean**2
    return np.var(X, axis=0)


def _compute_std_sparse(X: Union[np.ndarray, csr_matrix]) -> np.ndarray:
    """Compute standard deviation per gene efficiently.

    Parameters
    ----------
    X
        Expression matrix, shape (n_cells, n_genes). Can be sparse or dense.

    Returns
    -------
    np.ndarray
        Standard deviation per gene, shape (n_genes,).
    """
    return np.sqrt(_compute_variance_sparse(X))


# =============================================================================
# FDR CORRECTION (scipy-only implementation)
# =============================================================================


def _fdr_correction_bh(p_values: np.ndarray) -> np.ndarray:
    """Apply Benjamini-Hochberg FDR correction.

    Parameters
    ----------
    p_values
        Raw p-values, shape (n,).

    Returns
    -------
    np.ndarray
        Adjusted p-values, shape (n,).
    """
    n = len(p_values)
    if n == 0:
        return p_values.copy()

    # Sort p-values
    order = np.argsort(p_values)
    sorted_p = p_values[order]

    # Compute adjusted p-values: p_adj = p * n / rank
    ranks = np.arange(1, n + 1)
    p_adj_sorted = sorted_p * n / ranks

    # Ensure monotonicity via cumulative minimum from the end
    p_adj_sorted = np.minimum.accumulate(p_adj_sorted[::-1])[::-1]

    # Restore original order
    p_adj = np.empty(n)
    p_adj[order] = p_adj_sorted

    return np.clip(p_adj, 0, 1)


def _fdr_correction_bonferroni(p_values: np.ndarray) -> np.ndarray:
    """Apply Bonferroni correction.

    Parameters
    ----------
    p_values
        Raw p-values, shape (n,).

    Returns
    -------
    np.ndarray
        Adjusted p-values, shape (n,).
    """
    n = len(p_values)
    if n == 0:
        return p_values.copy()
    return np.clip(p_values * n, 0, 1)


def _apply_fdr_correction(
    p_values: np.ndarray,
    method: Literal["bonferroni", "fdr_bh", "none"],
) -> np.ndarray:
    """Apply FDR correction based on method.

    Parameters
    ----------
    p_values
        Raw p-values, shape (n,).
    method
        Correction method: "bonferroni", "fdr_bh", or "none".

    Returns
    -------
    np.ndarray
        Adjusted p-values, shape (n,).
    """
    if method == "none":
        return p_values.copy()
    elif method == "bonferroni":
        return _fdr_correction_bonferroni(p_values)
    elif method == "fdr_bh":
        return _fdr_correction_bh(p_values)
    else:
        raise ValueError(f"Unknown FDR method: {method}")


# =============================================================================
# QUADRANT CLASSIFICATION
# =============================================================================


def _classify_quadrants(
    z_values: np.ndarray,
    lag_values: np.ndarray,
    p_values: Optional[np.ndarray] = None,
    alpha: float = 0.05,
) -> np.ndarray:
    """Classify cells into LISA quadrants.

    Classification:
        - HH (1): High value, high lag (hotspot)
        - LL (2): Low value, low lag (coldspot)
        - HL (3): High value, low lag (high outlier)
        - LH (4): Low value, high lag (low outlier)
        - NS (0): Not significant

    Parameters
    ----------
    z_values
        Standardized values, shape (n_cells,) or (n_cells, n_genes).
    lag_values
        Spatial lag values, shape (n_cells,) or (n_cells, n_genes).
    p_values
        Optional p-values for significance filtering.
    alpha
        Significance threshold for p-values.

    Returns
    -------
    np.ndarray
        Quadrant labels as int8: 0=NS, 1=HH, 2=LL, 3=HL, 4=LH.
    """
    quadrants = np.zeros(z_values.shape, dtype=np.int8)

    # HH: high value, high lag
    quadrants[(z_values > 0) & (lag_values > 0)] = 1
    # LL: low value, low lag
    quadrants[(z_values < 0) & (lag_values < 0)] = 2
    # HL: high value, low lag
    quadrants[(z_values > 0) & (lag_values < 0)] = 3
    # LH: low value, high lag
    quadrants[(z_values < 0) & (lag_values > 0)] = 4

    # Filter by significance if p-values provided
    if p_values is not None:
        quadrants[p_values >= alpha] = 0

    return quadrants


# =============================================================================
# LEE'S L CORE COMPUTATION
# =============================================================================


def _compute_lees_l_core(
    z_x: np.ndarray,
    z_y: np.ndarray,
    W: csr_matrix,
    n_permutations: int,
    rng: np.random.Generator,
) -> tuple:
    """Core Lee's L computation on pre-standardized values.

    Computes bivariate spatial correlation between two variables.

    Parameters
    ----------
    z_x
        Standardized values of first variable, shape (n_cells,).
    z_y
        Standardized values of second variable, shape (n_cells,).
    W
        Row-normalized spatial weights matrix, shape (n_cells, n_cells).
    n_permutations
        Number of permutations for p-value computation.
    rng
        NumPy random generator for reproducibility.

    Returns
    -------
    tuple
        (L_local, L_global, lag_z_y, p_value) where:
        - L_local: Per-cell local Lee's L values, shape (n_cells,)
        - L_global: Global Lee's L statistic (float)
        - lag_z_y: Spatial lag of z_y, shape (n_cells,)
        - p_value: Permutation-based p-value for global L (float)
    """
    # Spatial lag of Y
    lag_z_y = W @ z_y
    if sparse.issparse(lag_z_y):
        lag_z_y = np.asarray(lag_z_y).flatten()

    # Local Lee's L: product of z_x and spatial lag of z_y
    L_local = z_x * lag_z_y

    # Global Lee's L: sum of local values
    L_global = float(L_local.sum())

    # Permutation test for global L
    p_value = 1.0
    if n_permutations > 0:
        L_permuted = np.zeros(n_permutations)

        for p in range(n_permutations):
            # Shuffle z_y to break spatial relationship
            z_y_perm = rng.permutation(z_y)
            lag_perm = W @ z_y_perm
            if sparse.issparse(lag_perm):
                lag_perm = np.asarray(lag_perm).flatten()
            L_permuted[p] = (z_x * lag_perm).sum()

        # Two-tailed p-value
        extreme = np.sum(np.abs(L_permuted) >= np.abs(L_global))
        p_value = float((extreme + 1) / (n_permutations + 1))

    return L_local, L_global, lag_z_y, p_value


# =============================================================================
# SPATIAL WEIGHTS CONSTRUCTION
# =============================================================================


def build_spatial_weights(
    adata: ad.AnnData,
    n_neighbors: int = 6,
    spatial_key: str = "spatial",
    include_self: bool = False,
) -> csr_matrix:
    """Build row-normalized sparse spatial weights matrix.

    Uses k-nearest neighbors to define spatial relationships between cells.
    Returns a sparse CSR matrix for memory efficiency.

    Parameters
    ----------
    adata
        AnnData object with spatial coordinates in adata.obsm[spatial_key].
    n_neighbors
        Number of nearest neighbors per cell. Default: 6.
    spatial_key
        Key in adata.obsm containing spatial coordinates. Default: "spatial".
    include_self
        If True, include self-connections (diagonal = 1). Default: False.

    Returns
    -------
    csr_matrix
        Row-normalized sparse weights matrix, shape (n_cells, n_cells).
        Each row sums to 1.

    Raises
    ------
    ValueError
        If spatial_key not found in adata.obsm.

    Examples
    --------
    >>> W = build_spatial_weights(adata, n_neighbors=6)
    >>> print(f"Weights: {W.shape}, nnz={W.nnz}")
    """
    if spatial_key not in adata.obsm:
        raise ValueError(
            f"adata.obsm['{spatial_key}'] not found. "
            "Spatial coordinates are required."
        )

    coordinates = adata.obsm[spatial_key]
    n_cells = coordinates.shape[0]

    logger.debug(f"Building spatial weights: {n_cells:,} cells, k={n_neighbors}")

    # Find k nearest neighbors using Ball Tree
    # Add 1 because query includes self
    nn = NearestNeighbors(n_neighbors=n_neighbors + 1, algorithm="ball_tree")
    nn.fit(coordinates)
    distances, indices = nn.kneighbors(coordinates)

    # Skip first column (self) unless include_self=True
    start_col = 0 if include_self else 1

    row_indices = np.repeat(np.arange(n_cells), n_neighbors + 1 - start_col)
    col_indices = indices[:, start_col:].flatten()
    data = np.ones(len(row_indices), dtype=np.float32)

    W = csr_matrix((data, (row_indices, col_indices)), shape=(n_cells, n_cells))

    # Row-normalize: each row sums to 1
    row_sums = np.array(W.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    W = W.multiply(1.0 / row_sums[:, np.newaxis])

    logger.debug(f"Spatial weights: nnz={W.nnz:,}")

    return W.tocsr()


# =============================================================================
# GLOBAL MORAN'S I
# =============================================================================


def morans_i(
    adata: ad.AnnData,
    genes: Optional[Union[str, List[str]]] = None,
    layer: Optional[str] = None,
    spatial_key: str = "spatial",
    n_neighbors: int = 6,
    n_permutations: int = 10,
    seed: int = 0,
    key_added: str = "morans_i",
    copy: bool = False,
    use_existing_graph: bool = False,
) -> ad.AnnData:
    """Compute Global Moran's I for spatial autocorrelation analysis.

    Global Moran's I measures the overall spatial autocorrelation of gene
    expression across all cells. Positive values indicate spatial clustering
    (similar values near each other), negative values indicate dispersion.

    Parameters
    ----------
    adata
        AnnData object with spatial coordinates in adata.obsm[spatial_key].
    genes
        Gene(s) to analyze. Can be:
        - None: analyze all genes (not recommended for large datasets)
        - str: single gene name
        - List[str]: list of gene names
    layer
        Layer to use for expression values. If None, uses adata.X.
    spatial_key
        Key in adata.obsm containing spatial coordinates. Default: "spatial".
    n_neighbors
        Number of neighbors for spatial weights. Default: 6.
    n_permutations
        Number of permutations for p-value computation. Default: 10.
        For publication-quality results, increase to 100+ (ideally 999).
        More permutations = more accurate p-values but slower computation.
    seed
        Random seed for permutation testing. Default: 0.
    key_added
        Key to store results in adata.uns. Default: "morans_i".
    copy
        If True, operate on a copy of adata. Default: False.
    use_existing_graph
        If True, reuse existing spatial connectivity graph in adata.obsp
        instead of rebuilding. Default: False (always rebuild to ensure
        n_neighbors matches).

    Returns
    -------
    AnnData
        AnnData with results in adata.uns[key_added] as DataFrame with columns:
        gene, I, expected_I, z_score, p_value

    Raises
    ------
    ValueError
        If spatial_key not found, or genes not in var_names.

    Notes
    -----
    This function uses squidpy's ``spatial_autocorr`` as the computational
    backend for vectorized, memory-efficient computation. SpatialCore adds:

    - Consistent API with other spatialcore.spatial functions
    - Integrated logging via spatialcore.core.logging
    - Metadata tracking via spatialcore.core.metadata
    - Standardized output format in adata.uns

    For local (per-cell) Moran's I with LISA quadrant classification,
    use ``local_morans_i`` which provides additional features not available
    in squidpy (FDR correction, spatial lag storage, batch processing).

    The Global Moran's I statistic is:

        I = (N / S0) * [sum_ij w_ij (x_i - x_mean)(x_j - x_mean)] / [sum_i (x_i - x_mean)^2]

    Where N is the number of cells, S0 is the sum of weights, and w_ij is
    the spatial weight between cells i and j.

    Interpretation:
        - I > 0: Positive autocorrelation (clustering)
        - I ~ 0: Random spatial pattern
        - I < 0: Negative autocorrelation (dispersion)

    Examples
    --------
    >>> import scanpy as sc
    >>> from spatialcore.spatial import morans_i
    >>> adata = sc.read_h5ad("xenium.h5ad")
    >>> adata = morans_i(adata, genes=["CD8A", "FOXP3"], n_permutations=100)
    >>> print(adata.uns["morans_i"])
    """
    start_time = time.time()

    # Input validation
    if spatial_key not in adata.obsm:
        raise ValueError(
            f"adata.obsm['{spatial_key}'] not found. "
            "Spatial coordinates are required."
        )

    if n_neighbors < 1:
        raise ValueError(f"n_neighbors must be >= 1, got {n_neighbors}")

    if n_permutations < 0:
        raise ValueError(f"n_permutations must be >= 0, got {n_permutations}")

    # Handle copy
    adata = adata.copy() if copy else adata

    # Resolve genes
    if genes is None:
        gene_names = list(adata.var_names)
        logger.warning(
            f"No genes specified, analyzing all {len(gene_names)} genes. "
            "This may be slow for large datasets."
        )
    elif isinstance(genes, str):
        gene_names = [genes]
    else:
        gene_names = list(genes)

    # Validate genes exist
    missing = set(gene_names) - set(adata.var_names)
    if missing:
        raise ValueError(f"Genes not found in adata.var_names: {list(missing)[:10]}")

    n_cells = adata.n_obs
    n_genes = len(gene_names)

    logger.info(
        f"Computing Global Moran's I: {n_cells:,} cells, {n_genes} genes, "
        f"k={n_neighbors}, permutations={n_permutations}"
    )

    # Build spatial neighbors using squidpy
    rebuild_graph = True
    if use_existing_graph and "spatial_connectivities" in adata.obsp:
        rebuild_graph = False
        logger.info("Using existing spatial connectivity graph (use_existing_graph=True)")

    if rebuild_graph:
        logger.debug(f"Building spatial neighbors graph (k={n_neighbors})")
        sq.gr.spatial_neighbors(
            adata,
            n_neighs=n_neighbors,
            coord_type="generic",
            spatial_key=spatial_key,
        )

    # Subset to genes of interest for squidpy computation
    adata_sub = adata[:, gene_names].copy()

    # Use squidpy for fast, vectorized computation
    sq.gr.spatial_autocorr(
        adata_sub,
        mode="moran",
        n_perms=n_permutations if n_permutations > 0 else None,
        n_jobs=1,
        seed=seed,
        layer=layer,
    )

    # Extract and reformat results to our standard format
    squidpy_results = adata_sub.uns["moranI"].copy()

    # Expected I under null hypothesis
    expected_I = -1 / (n_cells - 1)

    # Build results DataFrame in our format
    results = []
    for gene_name in gene_names:
        if gene_name in squidpy_results.index:
            row = squidpy_results.loc[gene_name]
            I_value = float(row["I"])

            # Get p-value (prefer permutation-based if available)
            if n_permutations > 0 and "pval_sim" in row:
                p_value = float(row["pval_sim"])
            else:
                p_value = float(row.get("pval_norm", 1.0))

            # Get z-score
            if "var_norm" in row and row["var_norm"] > 0:
                z_score = float((I_value - expected_I) / np.sqrt(row["var_norm"]))
            else:
                z_score = 0.0

            results.append({
                "gene": gene_name,
                "I": I_value,
                "expected_I": expected_I,
                "z_score": z_score,
                "p_value": p_value,
            })
        else:
            raise RuntimeError(
                f"Gene '{gene_name}' was passed to squidpy but not found in results. "
                "This indicates an internal error in squidpy or data corruption."
            )

    # Store results
    results_df = pd.DataFrame(results)
    adata.uns[key_added] = results_df

    elapsed = time.time() - start_time
    logger.info(f"Global Moran's I completed in {elapsed:.1f}s")

    # Metadata tracking
    update_metadata(
        adata,
        function_name="morans_i",
        parameters={
            "genes": gene_names[:10] if len(gene_names) > 10 else gene_names,
            "n_genes": n_genes,
            "n_neighbors": n_neighbors,
            "n_permutations": n_permutations,
            "use_existing_graph": use_existing_graph,
            "seed": seed,
            "backend": "squidpy",
        },
        outputs={
            "uns": key_added,
        },
    )

    return adata


# =============================================================================
# LOCAL MORAN'S I (LISA)
# =============================================================================


def local_morans_i(
    adata: ad.AnnData,
    genes: Optional[Union[str, List[str]]] = None,
    layer: Optional[str] = None,
    spatial_key: str = "spatial",
    n_neighbors: int = 6,
    n_permutations: int = 10,
    fdr_correction: Literal["bonferroni", "fdr_bh", "none"] = "fdr_bh",
    alpha: float = 0.05,
    seed: int = 0,
    batch_size: int = 100,
    key_added: str = "local_morans",
    copy: bool = False,
) -> ad.AnnData:
    """Compute Local Moran's I (LISA) for spatial hotspot detection.

    Local Moran's I identifies spatial clusters and outliers for each cell.
    Cells are classified into quadrants based on their value and the
    values of their neighbors.

    Parameters
    ----------
    adata
        AnnData object with spatial coordinates in adata.obsm[spatial_key].
    genes
        Gene(s) to analyze. Can be:
        - None: analyze all genes (not recommended for large datasets)
        - str: single gene name
        - List[str]: list of gene names
    layer
        Layer to use for expression values. If None, uses adata.X.
    spatial_key
        Key in adata.obsm containing spatial coordinates. Default: "spatial".
    n_neighbors
        Number of neighbors for spatial weights. Default: 6.
    n_permutations
        Number of permutations for p-value computation. Default: 10.
        For publication-quality results, increase to 100+ (ideally 999).
        More permutations = more accurate p-values but slower computation.
    fdr_correction
        Multiple testing correction method:
        - "bonferroni": Bonferroni correction (conservative)
        - "fdr_bh": Benjamini-Hochberg FDR (default, recommended)
        - "none": No correction
    alpha
        Significance threshold for quadrant classification. Default: 0.05.
    seed
        Random seed for permutation testing. Default: 0.
    batch_size
        Number of genes to process per batch. Default: 100.
        Smaller batches use less memory.
    key_added
        Prefix for result keys. Default: "local_morans".
    copy
        If True, operate on a copy of adata. Default: False.

    Returns
    -------
    AnnData
        AnnData with results in:
        - adata.obsm[f"{key_added}_I"]: Local I values (n_cells, n_genes)
        - adata.obsm[f"{key_added}_z"]: Z-scores (n_cells, n_genes)
        - adata.obsm[f"{key_added}_lag"]: Spatial lag values (n_cells, n_genes)
        - adata.obsm[f"{key_added}_p"]: P-values (n_cells, n_genes)
        - adata.obsm[f"{key_added}_p_adj"]: Adjusted p-values (n_cells, n_genes)
        - adata.obsm[f"{key_added}_quadrant"]: Quadrant labels (n_cells, n_genes)
        - adata.uns[f"{key_added}_params"]: Parameters and gene names

    Raises
    ------
    ValueError
        If spatial_key not found, or genes not in var_names.

    Notes
    -----
    The Local Moran's I statistic for cell i is:

        I_i = z_i * sum_j(w_ij * z_j)

    Where z_i = (x_i - mean) / std and w_ij is the spatial weight.

    Quadrant classification:
        - HH (1): High value, high neighbors (hotspot)
        - LL (2): Low value, low neighbors (coldspot)
        - HL (3): High value, low neighbors (high outlier)
        - LH (4): Low value, high neighbors (low outlier)
        - NS (0): Not significant

    Examples
    --------
    >>> import scanpy as sc
    >>> from spatialcore.spatial import local_morans_i
    >>> adata = sc.read_h5ad("xenium.h5ad")
    >>> adata = local_morans_i(adata, genes=["CD8A"], n_permutations=100)
    >>> # Get hotspots for CD8A
    >>> cd8a_hotspots = adata.obsm["local_morans_quadrant"][:, 0] == 1
    """
    start_time = time.time()

    # Input validation
    if spatial_key not in adata.obsm:
        raise ValueError(
            f"adata.obsm['{spatial_key}'] not found. "
            "Spatial coordinates are required."
        )

    if n_neighbors < 1:
        raise ValueError(f"n_neighbors must be >= 1, got {n_neighbors}")

    if n_permutations < 0:
        raise ValueError(f"n_permutations must be >= 0, got {n_permutations}")

    if fdr_correction not in ["bonferroni", "fdr_bh", "none"]:
        raise ValueError(
            f"Invalid fdr_correction: '{fdr_correction}'. "
            "Must be 'bonferroni', 'fdr_bh', or 'none'."
        )

    # Handle copy
    adata = adata.copy() if copy else adata

    # Resolve genes
    if genes is None:
        gene_names = list(adata.var_names)
        logger.warning(
            f"No genes specified, analyzing all {len(gene_names)} genes. "
            "This may be slow and memory-intensive."
        )
    elif isinstance(genes, str):
        gene_names = [genes]
    else:
        gene_names = list(genes)

    # Validate genes exist
    missing = set(gene_names) - set(adata.var_names)
    if missing:
        raise ValueError(f"Genes not found in adata.var_names: {list(missing)[:10]}")

    n_cells = adata.n_obs
    n_genes = len(gene_names)
    gene_indices = np.array([adata.var_names.get_loc(g) for g in gene_names])

    logger.info(
        f"Computing Local Moran's I: {n_cells:,} cells, {n_genes} genes, "
        f"k={n_neighbors}, permutations={n_permutations}"
    )

    # Build spatial weights
    W = build_spatial_weights(adata, n_neighbors=n_neighbors, spatial_key=spatial_key)

    # Get expression matrix
    if layer is not None:
        X = adata.layers[layer]
    else:
        X = adata.X

    # Ensure CSC format for efficient column slicing (CSR is O(nnz) for columns)
    if not sparse.issparse(X):
        X_sparse = sparse.csc_matrix(X)
    else:
        X_sparse = X.tocsc()

    # Pre-compute gene statistics using sparse-friendly operations (no densification)
    logger.debug("Pre-computing gene statistics (sparse)")
    X_genes = X_sparse[:, gene_indices]
    gene_means = _compute_mean_sparse(X_genes).astype(np.float32)
    gene_stds = _compute_std_sparse(X_genes).astype(np.float32)

    # Handle zero variance
    zero_var_mask = gene_stds == 0
    n_zero_var = zero_var_mask.sum()
    zero_variance_genes = [gene_names[i] for i in np.where(zero_var_mask)[0]]
    if n_zero_var > 0:
        logger.warning(f"{n_zero_var} genes have zero variance and will be skipped: {zero_variance_genes[:5]}")
        gene_stds[zero_var_mask] = 1.0  # Avoid division by zero

    # Initialize result arrays
    local_I = np.zeros((n_cells, n_genes), dtype=np.float32)
    z_values = np.zeros((n_cells, n_genes), dtype=np.float32)
    lag_values = np.zeros((n_cells, n_genes), dtype=np.float32)
    p_values = np.ones((n_cells, n_genes), dtype=np.float32)

    # Set random seed
    rng = np.random.default_rng(seed)

    # Process genes in batches
    n_batches = (n_genes + batch_size - 1) // batch_size
    logger.info(f"Processing {n_genes} genes in {n_batches} batches")

    for batch_idx in range(n_batches):
        batch_start = batch_idx * batch_size
        batch_end = min(batch_start + batch_size, n_genes)
        batch_gene_indices = gene_indices[batch_start:batch_end]

        logger.debug(f"Processing batch {batch_idx + 1}/{n_batches}")

        # Extract and densify batch
        X_batch = X_sparse[:, batch_gene_indices].toarray().astype(np.float32)

        # Standardize using pre-computed stats
        means_batch = gene_means[batch_start:batch_end]
        stds_batch = gene_stds[batch_start:batch_end]
        Z_batch = (X_batch - means_batch) / stds_batch

        # Store z-values
        z_values[:, batch_start:batch_end] = Z_batch

        # Compute spatial lag
        lag_batch = W @ Z_batch
        if sparse.issparse(lag_batch):
            lag_batch = lag_batch.toarray()
        lag_values[:, batch_start:batch_end] = lag_batch

        # Compute local I
        local_I[:, batch_start:batch_end] = Z_batch * lag_batch

        # Permutation test for p-values
        if n_permutations > 0:
            batch_size_actual = batch_end - batch_start
            I_perm = np.zeros((n_permutations, n_cells, batch_size_actual), dtype=np.float32)

            for p in range(n_permutations):
                # Shuffle z-values (same permutation for all genes in batch for efficiency)
                perm_idx = rng.permutation(n_cells)
                Z_shuffled = Z_batch[perm_idx, :]
                lag_shuffled = W @ Z_shuffled
                if sparse.issparse(lag_shuffled):
                    lag_shuffled = lag_shuffled.toarray()
                I_perm[p, :, :] = Z_shuffled * lag_shuffled

            # Compute p-values (two-tailed)
            local_I_batch = local_I[:, batch_start:batch_end]
            for gene_local_idx in range(batch_size_actual):
                for cell_idx in range(n_cells):
                    extreme = np.sum(
                        np.abs(I_perm[:, cell_idx, gene_local_idx])
                        >= np.abs(local_I_batch[cell_idx, gene_local_idx])
                    )
                    p_values[cell_idx, batch_start + gene_local_idx] = (extreme + 1) / (
                        n_permutations + 1
                    )

        # Free batch memory
        del X_batch, Z_batch, lag_batch

    # Set zero-variance genes to NaN/NS
    if zero_var_mask.any():
        local_I[:, zero_var_mask] = 0.0
        z_values[:, zero_var_mask] = 0.0
        lag_values[:, zero_var_mask] = 0.0
        p_values[:, zero_var_mask] = 1.0

    # Apply FDR correction and classify quadrants
    if n_permutations > 0:
        # Apply FDR correction (per gene)
        logger.debug(f"Applying {fdr_correction} correction")
        p_adj = np.ones_like(p_values)
        for gene_idx in range(n_genes):
            p_adj[:, gene_idx] = _apply_fdr_correction(p_values[:, gene_idx], fdr_correction)

        # Classify quadrants with significance filtering
        logger.debug("Classifying LISA quadrants (with significance filtering)")
        quadrants = _classify_quadrants(z_values, lag_values, p_adj, alpha)
    else:
        # No permutations: classify by z/lag signs only, no significance filtering
        logger.warning(
            "n_permutations=0: Quadrants classified by z/lag signs only, "
            "without significance filtering. Consider n_permutations>=99 for p-values."
        )
        p_adj = p_values  # All 1.0
        quadrants = _classify_quadrants(z_values, lag_values, p_values=None, alpha=alpha)

    # Store results
    adata.obsm[f"{key_added}_I"] = local_I
    adata.obsm[f"{key_added}_z"] = z_values
    adata.obsm[f"{key_added}_lag"] = lag_values
    adata.obsm[f"{key_added}_p"] = p_values
    adata.obsm[f"{key_added}_p_adj"] = p_adj
    adata.obsm[f"{key_added}_quadrant"] = quadrants

    # Store parameters
    elapsed = time.time() - start_time

    adata.uns[f"{key_added}_params"] = {
        "genes": gene_names,
        "n_neighbors": n_neighbors,
        "n_permutations": n_permutations,
        "fdr_correction": fdr_correction,
        "alpha": alpha,
        "n_cells": n_cells,
        "n_genes": n_genes,
        "seed": seed,
        "computation_time_seconds": elapsed,
        "zero_variance_genes": zero_variance_genes,
    }

    # Log summary
    n_significant = (quadrants != 0).sum(axis=0)
    logger.info(
        f"Local Moran's I completed in {elapsed:.1f}s. "
        f"Significant cells per gene: min={n_significant.min()}, max={n_significant.max()}"
    )

    # Metadata tracking
    update_metadata(
        adata,
        function_name="local_morans_i",
        parameters={
            "genes": gene_names[:10] if len(gene_names) > 10 else gene_names,
            "n_genes": n_genes,
            "n_neighbors": n_neighbors,
            "n_permutations": n_permutations,
            "fdr_correction": fdr_correction,
            "alpha": alpha,
            "seed": seed,
        },
        outputs={
            "obsm_I": f"{key_added}_I",
            "obsm_z": f"{key_added}_z",
            "obsm_lag": f"{key_added}_lag",
            "obsm_p": f"{key_added}_p",
            "obsm_p_adj": f"{key_added}_p_adj",
            "obsm_quadrant": f"{key_added}_quadrant",
            "uns_params": f"{key_added}_params",
        },
    )

    return adata


# =============================================================================
# GLOBAL LEE'S L
# =============================================================================


def lees_l(
    adata: ad.AnnData,
    gene_pairs: Union[Tuple[str, str], List[Tuple[str, str]]],
    layer: Optional[str] = None,
    spatial_key: str = "spatial",
    n_neighbors: int = 6,
    n_permutations: int = 199,
    seed: int = 0,
) -> Union[dict, List[dict]]:
    """Compute Global Lee's L bivariate spatial autocorrelation.

    Lee's L measures whether two genes co-vary spatially across all cells.
    Positive values indicate spatial co-localization, negative values
    indicate spatial anti-correlation.

    Parameters
    ----------
    adata
        AnnData object with spatial coordinates in adata.obsm[spatial_key].
    gene_pairs
        Gene pair(s) to analyze. Can be:
        - Tuple[str, str]: single pair like ("CD8A", "GZMB")
        - List[Tuple[str, str]]: multiple pairs like [("CD8A", "GZMB"), ("EPCAM", "KRT8")]
    layer
        Layer to use for expression values. If None, uses adata.X.
    spatial_key
        Key in adata.obsm containing spatial coordinates. Default: "spatial".
    n_neighbors
        Number of neighbors for spatial weights. Default: 6.
    n_permutations
        Number of permutations for p-value computation. Default: 199.
    seed
        Random seed for permutation testing. Default: 0.

    Returns
    -------
    Union[dict, List[dict]]
        For single pair: dict with keys "L", "p_value", "gene_x", "gene_y"
        For multiple pairs: list of dicts, one per pair

    Raises
    ------
    ValueError
        If spatial_key not found, or genes not in var_names.

    Notes
    -----
    The Global Lee's L statistic is:

        L = sum_i(z_x_i * lag(z_y)_i)

    Where z_x and z_y are standardized gene expression values and
    lag(z_y) is the spatial lag (weighted average of neighbors' z_y values).

    Interpretation:
        - L > 0: Positive spatial correlation (co-localization)
        - L ~ 0: No spatial relationship
        - L < 0: Negative spatial correlation (anti-correlation)

    Examples
    --------
    >>> from spatialcore.spatial import lees_l
    >>> result = lees_l(adata, gene_pairs=("CD8A", "GZMB"))
    >>> print(f"Lee's L = {result['L']:.4f}, p = {result['p_value']:.4f}")
    """
    start_time = time.time()

    # Input validation
    if spatial_key not in adata.obsm:
        raise ValueError(
            f"adata.obsm['{spatial_key}'] not found. "
            "Spatial coordinates are required."
        )

    if n_neighbors < 1:
        raise ValueError(f"n_neighbors must be >= 1, got {n_neighbors}")

    if n_permutations < 0:
        raise ValueError(f"n_permutations must be >= 0, got {n_permutations}")

    # Normalize gene_pairs to list
    if isinstance(gene_pairs, tuple) and len(gene_pairs) == 2:
        if isinstance(gene_pairs[0], str):
            gene_pairs = [gene_pairs]
            single_pair = True
        else:
            single_pair = False
    else:
        single_pair = False

    # Validate genes exist
    all_genes = set(g for pair in gene_pairs for g in pair)
    missing = all_genes - set(adata.var_names)
    if missing:
        raise ValueError(f"Genes not found in adata.var_names: {list(missing)}")

    n_cells = adata.n_obs
    n_pairs = len(gene_pairs)

    logger.info(
        f"Computing Global Lee's L: {n_cells:,} cells, {n_pairs} pair(s), "
        f"k={n_neighbors}, permutations={n_permutations}"
    )

    # Build spatial weights
    W = build_spatial_weights(adata, n_neighbors=n_neighbors, spatial_key=spatial_key)

    # Get expression matrix
    if layer is not None:
        X = adata.layers[layer]
    else:
        X = adata.X

    # Convert to CSC for efficient column slicing
    if sparse.issparse(X):
        X = X.tocsc()

    # Set random seed
    rng = np.random.default_rng(seed)

    # Process each pair
    results = []
    for gene_x, gene_y in gene_pairs:
        # Extract gene expression
        gene_x_idx = adata.var_names.get_loc(gene_x)
        gene_y_idx = adata.var_names.get_loc(gene_y)

        if sparse.issparse(X):
            x_vals = np.asarray(X[:, gene_x_idx].toarray()).flatten()
            y_vals = np.asarray(X[:, gene_y_idx].toarray()).flatten()
        else:
            x_vals = np.asarray(X[:, gene_x_idx]).flatten()
            y_vals = np.asarray(X[:, gene_y_idx]).flatten()

        # Standardize
        x_std = x_vals.std()
        y_std = y_vals.std()

        if x_std == 0 or y_std == 0:
            logger.warning(
                f"Gene pair ({gene_x}, {gene_y}) has zero variance gene - "
                "setting L to 0"
            )
            results.append({
                "gene_x": gene_x,
                "gene_y": gene_y,
                "L": 0.0,
                "p_value": 1.0,
            })
            continue

        z_x = (x_vals - x_vals.mean()) / x_std
        z_y = (y_vals - y_vals.mean()) / y_std

        # Compute Lee's L
        L_local, L_global, lag_z_y, p_value = _compute_lees_l_core(
            z_x, z_y, W, n_permutations, rng
        )

        results.append({
            "gene_x": gene_x,
            "gene_y": gene_y,
            "L": L_global,
            "p_value": p_value,
        })

    elapsed = time.time() - start_time
    logger.info(f"Global Lee's L completed in {elapsed:.1f}s")

    # Return single dict for single pair, list for multiple
    if single_pair:
        return results[0]
    return results


# =============================================================================
# LOCAL LEE'S L
# =============================================================================


def lees_l_local(
    adata: ad.AnnData,
    gene_pairs: Optional[Union[Tuple[str, str], List[Tuple[str, str]]]] = None,
    genes: Optional[List[str]] = None,
    layer: Optional[str] = None,
    spatial_key: str = "spatial",
    n_neighbors: int = 6,
    n_permutations: int = 199,
    compute_cell_pvalues: bool = False,
    significance_filter: bool = False,
    alpha: float = 0.05,
    seed: int = 0,
    copy: bool = False,
) -> ad.AnnData:
    """Compute Local Lee's L for bivariate spatial hotspot detection.

    Local Lee's L identifies spatial clusters and outliers for gene pairs.
    Cells are classified into quadrants based on their expression values
    and neighbors' expression.

    Parameters
    ----------
    adata
        AnnData object with spatial coordinates in adata.obsm[spatial_key].
    gene_pairs
        Gene pair(s) to analyze. Can be:
        - Tuple[str, str]: single pair like ("CD8A", "GZMB")
        - List[Tuple[str, str]]: multiple pairs
    genes
        Alternative: list of genes to compute all pairwise combinations.
        WARNING: For N genes, this computes N*(N-1)/2 pairs which can be
        very slow for large gene sets. Not recommended for > 10 genes.
    layer
        Layer to use for expression values. If None, uses adata.X.
    spatial_key
        Key in adata.obsm containing spatial coordinates. Default: "spatial".
    n_neighbors
        Number of neighbors for spatial weights. Default: 6.
    n_permutations
        Number of permutations for global p-value. Default: 199.
        Set to 0 to skip permutation testing entirely.
    compute_cell_pvalues
        If True, compute per-cell p-values via permutation testing.
        This is expensive (O(n_permutations * n_cells)) and usually
        unnecessary unless significance_filter=True. Default: False.
    significance_filter
        If True, filter quadrant classifications by p-value significance.
        Cells with p >= alpha are classified as NS. Requires
        compute_cell_pvalues=True. Default: False.
    alpha
        Significance threshold for quadrant classification when
        significance_filter=True. Default: 0.05.
    seed
        Random seed for permutation testing. Default: 0.
    copy
        If True, operate on a copy of adata. Default: False.

    Returns
    -------
    AnnData
        AnnData with results in:
        - adata.obs["{gene_x}_{gene_y}_lees_l"]: Local L values per cell
        - adata.obs["{gene_x}_{gene_y}_quadrant"]: HH/HL/LH/LL/NS classification
        - adata.obs["{gene_x}_{gene_y}_pvalue"]: Per-cell p-values (if compute_cell_pvalues=True)
        - adata.uns["{gene_x}_{gene_y}_lees_l_params"]: Parameters and global L

    Raises
    ------
    ValueError
        If neither gene_pairs nor genes provided, or genes not found.

    Notes
    -----
    The Local Lee's L statistic for cell i is:

        L_i = z_x_i * lag(z_y)_i

    Where z_x_i is the standardized expression of gene X at cell i, and
    lag(z_y)_i is the weighted average of neighbors' standardized gene Y values.

    Quadrant classification uses z-score signs only (not p-value significance),
    following Lee (2001) and the standard LISA approach:
        - HH: z_x > 0 AND lag(z_y) > 0 (co-expression hotspot)
        - LL: z_x < 0 AND lag(z_y) < 0 (co-expression coldspot)
        - HL: z_x > 0 AND lag(z_y) < 0 (spatial outlier)
        - LH: z_x < 0 AND lag(z_y) > 0 (spatial outlier)
        - NS: z_x = 0 OR lag(z_y) = 0 exactly (rare)

    Examples
    --------
    >>> from spatialcore.spatial import lees_l_local
    >>> # Single pair
    >>> adata = lees_l_local(adata, gene_pairs=("CD8A", "GZMB"))
    >>> print(adata.obs["CD8A_GZMB_quadrant"].value_counts())
    >>> # Multiple pairs
    >>> adata = lees_l_local(adata, gene_pairs=[("CD8A", "GZMB"), ("EPCAM", "KRT8")])
    >>> # All pairs from gene list (use with caution!)
    >>> adata = lees_l_local(adata, genes=["CD8A", "GZMB", "FOXP3"])
    """
    start_time = time.time()

    # Input validation
    if gene_pairs is None and genes is None:
        raise ValueError(
            "Must provide either 'gene_pairs' or 'genes' parameter. "
            "Example: gene_pairs=('CD8A', 'GZMB') or genes=['CD8A', 'GZMB', 'FOXP3']"
        )

    if spatial_key not in adata.obsm:
        raise ValueError(
            f"adata.obsm['{spatial_key}'] not found. "
            "Spatial coordinates are required."
        )

    if n_neighbors < 1:
        raise ValueError(f"n_neighbors must be >= 1, got {n_neighbors}")

    if n_permutations < 0:
        raise ValueError(f"n_permutations must be >= 0, got {n_permutations}")

    if significance_filter and not compute_cell_pvalues:
        raise ValueError(
            "significance_filter=True requires compute_cell_pvalues=True"
        )

    # Handle all-pairs mode
    if genes is not None:
        n_pairs = len(genes) * (len(genes) - 1) // 2
        logger.warning(
            f"All-pairs mode: {len(genes)} genes = {n_pairs} pairs. "
            "This may take a very long time for large gene sets. "
            "Consider using explicit gene_pairs for better performance."
        )
        gene_pairs = list(combinations(genes, 2))
    else:
        # Normalize gene_pairs to list
        if isinstance(gene_pairs, tuple) and len(gene_pairs) == 2:
            if isinstance(gene_pairs[0], str):
                gene_pairs = [gene_pairs]

    # Validate genes exist
    all_genes = list(set(g for pair in gene_pairs for g in pair))
    missing = set(all_genes) - set(adata.var_names)
    if missing:
        raise ValueError(f"Genes not found in adata.var_names: {list(missing)}")

    # Handle copy
    adata = adata.copy() if copy else adata

    n_cells = adata.n_obs
    n_pairs = len(gene_pairs)

    logger.info(
        f"Computing Local Lee's L: {n_cells:,} cells, {n_pairs} pair(s), "
        f"k={n_neighbors}, permutations={n_permutations}"
    )

    # Build spatial weights (once)
    W = build_spatial_weights(adata, n_neighbors=n_neighbors, spatial_key=spatial_key)

    # Get expression matrix
    if layer is not None:
        X = adata.layers[layer]
    else:
        X = adata.X

    # Convert to CSC for efficient column slicing
    if sparse.issparse(X):
        X = X.tocsc()

    # Extract and standardize all unique genes (batch optimization)
    logger.debug(f"Extracting {len(all_genes)} unique genes")
    gene_data = {}
    for gene in all_genes:
        gene_idx = adata.var_names.get_loc(gene)
        if sparse.issparse(X):
            vals = np.asarray(X[:, gene_idx].toarray()).flatten()
        else:
            vals = np.asarray(X[:, gene_idx]).flatten()
        gene_data[gene] = vals

    # Standardize genes
    standardized = {}
    zero_var_genes = set()
    for gene, vals in gene_data.items():
        std = vals.std()
        if std == 0:
            zero_var_genes.add(gene)
            standardized[gene] = np.zeros_like(vals)
        else:
            standardized[gene] = (vals - vals.mean()) / std

    if zero_var_genes:
        logger.warning(f"Genes with zero variance: {zero_var_genes}")

    # Set random seed
    rng = np.random.default_rng(seed)

    # Process each pair
    for pair_idx, (gene_x, gene_y) in enumerate(gene_pairs):
        logger.debug(f"Processing pair {pair_idx + 1}/{n_pairs}: {gene_x} vs {gene_y}")

        z_x = standardized[gene_x]
        z_y = standardized[gene_y]

        # Skip if either gene has zero variance
        if gene_x in zero_var_genes or gene_y in zero_var_genes:
            key = f"{gene_x}_{gene_y}"
            adata.obs[f"{key}_lees_l"] = np.zeros(n_cells, dtype=np.float32)
            adata.obs[f"{key}_quadrant"] = pd.Categorical(
                ["NS"] * n_cells, categories=["NS", "HH", "LL", "HL", "LH"]
            )
            adata.uns[f"{key}_lees_l_params"] = {
                "gene_x": gene_x,
                "gene_y": gene_y,
                "global_L": 0.0,
                "global_pvalue": 1.0,
                "n_neighbors": n_neighbors,
                "zero_variance": True,
            }
            continue

        # Compute Lee's L
        L_local, L_global, lag_z_y, global_pvalue = _compute_lees_l_core(
            z_x, z_y, W, n_permutations, rng
        )

        # Compute per-cell p-values via permutation (optional, expensive)
        p_values = np.ones(n_cells, dtype=np.float32)
        if compute_cell_pvalues and n_permutations > 0:
            logger.debug(f"Computing per-cell p-values ({n_permutations} permutations)")
            L_perm = np.zeros((n_permutations, n_cells), dtype=np.float32)
            for p in range(n_permutations):
                z_y_perm = rng.permutation(z_y)
                lag_perm = W @ z_y_perm
                if sparse.issparse(lag_perm):
                    lag_perm = np.asarray(lag_perm).flatten()
                L_perm[p, :] = z_x * lag_perm

            # Two-tailed p-value per cell
            for i in range(n_cells):
                extreme = np.sum(np.abs(L_perm[:, i]) >= np.abs(L_local[i]))
                p_values[i] = (extreme + 1) / (n_permutations + 1)
        elif compute_cell_pvalues and n_permutations == 0:
            logger.warning("compute_cell_pvalues=True but n_permutations=0; p-values will be 1.0")

        # Classify quadrants
        if significance_filter:
            # Filter by p-value significance (requires compute_cell_pvalues=True)
            quadrants = _classify_quadrants(z_x, lag_z_y, p_values=p_values, alpha=alpha)
        else:
            # Classify by z-score signs only (no p-value filtering)
            # Following Lee (2001) and standard LISA approach
            quadrants = _classify_quadrants(z_x, lag_z_y, p_values=None, alpha=alpha)

        # Convert quadrant codes to labels
        quadrant_labels = [QUADRANT_LABELS[q] for q in quadrants]

        # Store results
        key = f"{gene_x}_{gene_y}"
        adata.obs[f"{key}_lees_l"] = L_local.astype(np.float32)
        adata.obs[f"{key}_quadrant"] = pd.Categorical(
            quadrant_labels, categories=["NS", "HH", "LL", "HL", "LH"]
        )
        adata.obs[f"{key}_pvalue"] = p_values.astype(np.float32)

        # Summarize quadrant counts
        quadrant_counts = {label: 0 for label in ["NS", "HH", "LL", "HL", "LH"]}
        for label in quadrant_labels:
            quadrant_counts[label] += 1

        adata.uns[f"{key}_lees_l_params"] = {
            "gene_x": gene_x,
            "gene_y": gene_y,
            "global_L": L_global,
            "global_pvalue": global_pvalue,
            "n_neighbors": n_neighbors,
            "n_permutations": n_permutations,
            "compute_cell_pvalues": compute_cell_pvalues,
            "significance_filter": significance_filter,
            "alpha": alpha,
            "quadrant_counts": quadrant_counts,
        }

    elapsed = time.time() - start_time
    logger.info(f"Local Lee's L completed in {elapsed:.1f}s for {n_pairs} pair(s)")

    # Metadata tracking
    pair_keys = [f"{gx}_{gy}" for gx, gy in gene_pairs]
    update_metadata(
        adata,
        function_name="lees_l_local",
        parameters={
            "gene_pairs": [(gx, gy) for gx, gy in gene_pairs[:10]],
            "n_pairs": n_pairs,
            "n_neighbors": n_neighbors,
            "n_permutations": n_permutations,
            "compute_cell_pvalues": compute_cell_pvalues,
            "significance_filter": significance_filter,
            "alpha": alpha,
            "seed": seed,
        },
        outputs={
            "obs_keys": [f"{k}_lees_l" for k in pair_keys[:5]],
            "uns_keys": [f"{k}_lees_l_params" for k in pair_keys[:5]],
        },
    )

    return adata
