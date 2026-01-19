"""
Memory-efficient AnnData loading and subsampling utilities.

This module provides utilities for:
1. Loading large h5ad files using backed mode to minimize memory usage
2. Stratified subsampling to maintain cell type proportions
3. Normalization validation and application

For large CellxGene datasets (>2GB), backed mode loads only metadata initially,
allowing efficient subsampling before loading expression data.

References:
    - Scanpy backed mode: https://scanpy.readthedocs.io/en/stable/generated/scanpy.read_h5ad.html
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any
import gc

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from spatialcore.core.logging import get_logger

logger = get_logger(__name__)


# ============================================================================
# Memory Utilities
# ============================================================================

def get_available_memory_gb() -> float:
    """
    Get available system memory in GB.

    Returns
    -------
    float
        Available memory in gigabytes.

    Notes
    -----
    Requires psutil. Returns 0 if psutil is not available.
    """
    try:
        import psutil
        return psutil.virtual_memory().available / (1024**3)
    except ImportError:
        logger.warning("psutil not installed, cannot check available memory")
        return 0.0


def estimate_adata_memory_gb(n_cells: int, n_genes: int, dtype_bytes: int = 4) -> float:
    """
    Estimate memory required for a dense expression matrix.

    Parameters
    ----------
    n_cells : int
        Number of cells.
    n_genes : int
        Number of genes.
    dtype_bytes : int, default 4
        Bytes per value (4 for float32, 8 for float64).

    Returns
    -------
    float
        Estimated memory in gigabytes.
    """
    return (n_cells * n_genes * dtype_bytes) / (1024**3)


# ============================================================================
# Stratified Sampling
# ============================================================================

def _stratified_sample_indices(
    obs_df: pd.DataFrame,
    label_column: str,
    max_cells: int,
    random_state: int = 42,
) -> np.ndarray:
    """
    Get stratified sample indices maintaining cell type proportions.

    Works with backed AnnData by only accessing .obs metadata.

    Parameters
    ----------
    obs_df : pd.DataFrame
        Observation dataframe (adata.obs).
    label_column : str
        Column containing cell type labels.
    max_cells : int
        Maximum number of cells to sample.
    random_state : int, default 42
        Random seed for reproducibility.

    Returns
    -------
    np.ndarray
        Sorted array of selected cell indices.
    """
    np.random.seed(random_state)

    labels = obs_df[label_column].values
    unique_labels = np.unique(labels)

    # Calculate cells per type (proportional)
    label_counts = pd.Series(labels).value_counts()
    total_cells = len(labels)
    indices = []

    for label in unique_labels:
        label_indices = np.where(labels == label)[0]
        # Proportional allocation
        n_sample = int(np.ceil(max_cells * len(label_indices) / total_cells))
        n_sample = min(n_sample, len(label_indices))

        sampled = np.random.choice(label_indices, size=n_sample, replace=False)
        indices.extend(sampled)

    # Trim to exact max_cells if we oversampled due to ceiling
    indices = np.array(indices)
    if len(indices) > max_cells:
        indices = np.random.choice(indices, size=max_cells, replace=False)

    return np.sort(indices)


def subsample_adata(
    adata: ad.AnnData,
    max_cells: int,
    stratify_by: Optional[str] = None,
    random_state: int = 42,
    copy: bool = True,
) -> ad.AnnData:
    """
    Subsample AnnData to max_cells, optionally maintaining cell type proportions.

    Parameters
    ----------
    adata : AnnData
        AnnData object to subsample.
    max_cells : int
        Maximum number of cells to keep.
    stratify_by : str, optional
        Column in adata.obs to use for stratified sampling.
        If None, random sampling is used.
    random_state : int, default 42
        Random seed for reproducibility.
    copy : bool, default True
        If True, return a copy. Otherwise modifies in place.

    Returns
    -------
    AnnData
        Subsampled AnnData object.

    Examples
    --------
    >>> from spatialcore.annotation import subsample_adata
    >>> adata_small = subsample_adata(adata, max_cells=10000, stratify_by="cell_type")
    """
    if adata.n_obs <= max_cells:
        logger.info(f"AnnData has {adata.n_obs:,} cells, no subsampling needed")
        return adata.copy() if copy else adata

    np.random.seed(random_state)

    if stratify_by and stratify_by in adata.obs.columns:
        indices = _stratified_sample_indices(
            adata.obs, stratify_by, max_cells, random_state
        )
        logger.info(
            f"Stratified subsampling by '{stratify_by}': "
            f"{adata.n_obs:,} → {len(indices):,} cells"
        )
    else:
        indices = np.random.choice(adata.n_obs, size=max_cells, replace=False)
        logger.info(f"Random subsampling: {adata.n_obs:,} → {max_cells:,} cells")

    return adata[indices].copy() if copy else adata[indices]


# ============================================================================
# Memory-Efficient Loading
# ============================================================================

def load_adata_backed(
    path: Union[str, Path],
    max_cells: Optional[int] = None,
    label_column: Optional[str] = None,
    large_file_threshold_gb: float = 2.0,
    random_state: int = 42,
) -> ad.AnnData:
    """
    Load AnnData with memory-efficient strategies for large files.

    Strategy:
    - Small files (<2GB): Load fully into memory
    - Large files (>=2GB): Use backed mode, subsample indices, then load subset

    Parameters
    ----------
    path : str or Path
        Path to h5ad file.
    max_cells : int, optional
        Maximum cells to load. If None, loads all cells.
    label_column : str, optional
        Column for stratified sampling. If None, random sampling.
    large_file_threshold_gb : float, default 2.0
        File size threshold for using backed mode.
    random_state : int, default 42
        Random seed for reproducibility.

    Returns
    -------
    AnnData
        AnnData loaded into memory (not backed).

    Notes
    -----
    For large files, this function:
    1. Opens in backed mode (only loads metadata)
    2. Determines subsample indices from metadata only
    3. Loads ONLY the selected cells into memory

    This can reduce memory usage significantly. For example, HLCA (584k cells)
    can be subsampled to 100k cells without ever loading the full matrix.

    Examples
    --------
    >>> from spatialcore.annotation import load_adata_backed
    >>> # Load up to 100k cells with stratified sampling
    >>> adata = load_adata_backed(
    ...     "large_dataset.h5ad",
    ...     max_cells=100000,
    ...     label_column="cell_type"
    ... )
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {path}")

    file_size_gb = path.stat().st_size / (1024**3)
    available_memory = get_available_memory_gb()

    logger.info(f"Loading: {path.name} ({file_size_gb:.2f} GB)")
    if available_memory > 0:
        logger.info(f"  Available memory: {available_memory:.1f} GB")

    if file_size_gb >= large_file_threshold_gb:
        logger.info(f"  Large file - using backed mode for memory efficiency")

        # Step 1: Open in backed mode (only loads metadata, not expression data)
        adata_backed = sc.read_h5ad(str(path), backed='r')
        logger.info(
            f"  Opened: {adata_backed.n_obs:,} cells × {adata_backed.n_vars:,} genes (backed)"
        )

        # Step 2: Determine subsample indices BEFORE loading data
        n_cells = adata_backed.n_obs
        if max_cells is not None and n_cells > max_cells:
            if label_column and label_column in adata_backed.obs.columns:
                indices = _stratified_sample_indices(
                    adata_backed.obs, label_column, max_cells, random_state
                )
                logger.info(
                    f"  Stratified subsampling: {n_cells:,} → {len(indices):,} cells"
                )
            else:
                np.random.seed(random_state)
                indices = np.random.choice(n_cells, size=max_cells, replace=False)
                logger.info(f"  Random subsampling: {n_cells:,} → {max_cells:,} cells")
        else:
            indices = np.arange(n_cells)

        # Step 3: Load only the selected cells into memory
        adata = adata_backed[indices].to_memory()
        logger.info(
            f"  Loaded into memory: {adata.n_obs:,} cells × {adata.n_vars:,} genes"
        )

        # Clean up backed reference
        adata_backed.file.close()
        gc.collect()

    else:
        # Small file - load directly
        logger.info(f"  Loading full file into memory")
        adata = sc.read_h5ad(str(path))

        # Subsample if needed
        if max_cells is not None and adata.n_obs > max_cells:
            adata = subsample_adata(
                adata, max_cells, stratify_by=label_column, random_state=random_state
            )

    return adata


# ============================================================================
# Normalization Utilities
# ============================================================================

def _copy_raw_to_x(adata: ad.AnnData, raw_source: str) -> None:
    """
    Copy raw counts from source location to adata.X (in-place).

    Parameters
    ----------
    adata : AnnData
        AnnData object to modify.
    raw_source : str
        Source location: "layers/{name}", "raw.X", or "X".
    """
    from scipy.sparse import issparse, csr_matrix

    if raw_source == "X":
        # Already in X, nothing to do
        return

    if raw_source == "raw.X":
        source_matrix = adata.raw.X
    elif raw_source.startswith("layers/"):
        layer_name = raw_source.split("/", 1)[1]
        source_matrix = adata.layers[layer_name]
    else:
        raise ValueError(f"Unknown raw_source: {raw_source}")

    # Copy to X, preserving sparsity
    if issparse(source_matrix):
        adata.X = source_matrix.copy()
    else:
        adata.X = np.array(source_matrix, copy=True)

    logger.info(f"Copied raw counts from {raw_source} to X")


def ensure_normalized(
    adata: ad.AnnData,
    target_sum: float = 1e4,
    unsafe_force: bool = False,
    copy: bool = False,
) -> ad.AnnData:
    """
    Ensure data is log1p normalized to target counts per cell.

    This function robustly detects the normalization state by:
    1. Searching for raw counts in layers["counts"], layers["raw_counts"],
       layers["raw"], adata.raw.X, and adata.X
    2. Verifying raw counts via integer test with floating-point tolerance
    3. Verifying log1p_10k via expm1 row sum estimation

    CellTypist REQUIRES: log1p(10k) with exclude_highly_expressed=False

    Parameters
    ----------
    adata : AnnData
        AnnData object to normalize.
    target_sum : float, default 1e4
        Target sum for normalization (10000 for CellTypist).
    unsafe_force : bool, default False
        **DANGEROUS**: If True, applies normalization even when data state
        cannot be verified. This may produce INCORRECT results if:

        - Data is already log-transformed (double-logging destroys signal)
        - Data uses a different target sum (e.g., CPM vs 10k)
        - Data contains negative values (z-scored/batch-corrected)
        - Data is latent space embeddings (not expression)

        Only use this if you have manually verified your data's state
        through other means. Incorrect normalization will produce
        systematically wrong cell type predictions.

        When enabled, logs a WARNING with the detected (unverified) state.
    copy : bool, default False
        If True, return a copy. Otherwise modifies in place.

    Returns
    -------
    AnnData
        Normalized AnnData object with log1p(10k) in X.

    Raises
    ------
    ValueError
        If no raw counts found and adata.X is not verified as log1p_10k,
        unless ``unsafe_force=True``.

    Notes
    -----
    **Safe normalization paths:**

    1. Raw counts found (in layers, raw.X, or X): Copy to X, normalize, log1p
    2. X verified as log1p_10k: No action needed

    **Unsafe states (require unsafe_force=True):**

    - log1p_cpm: Normalized to 1M instead of 10k
    - log1p_other: Unknown target sum
    - linear: Normalized but not log-transformed
    - negative: Contains negative values (z-scored?)
    - unknown: Cannot determine state

    Examples
    --------
    >>> from spatialcore.annotation import ensure_normalized
    >>> # Normal usage - will error if data state cannot be verified
    >>> adata = ensure_normalized(adata, target_sum=10000)

    >>> # Dangerous: force normalization on unverified data
    >>> adata = ensure_normalized(adata, unsafe_force=True)  # NOT RECOMMENDED
    """
    from spatialcore.core.utils import check_normalization_status

    if copy:
        adata = adata.copy()

    status = check_normalization_status(adata)

    logger.info(
        f"Normalization status: x_state={status['x_state']}, "
        f"raw_source={status['raw_source']}"
    )

    # Path 1: X is already log1p_10k - nothing to do
    if status["x_state"] == "log1p_10k":
        logger.info("Data already log1p normalized to 10k")
        return adata

    # Path 2: Raw counts available - normalize from raw
    if status["raw_source"] is not None:
        logger.info(f"Normalizing from raw counts ({status['raw_source']})")

        # Copy raw to X if not already there
        _copy_raw_to_x(adata, status["raw_source"])

        # Apply normalization
        # CRITICAL: exclude_highly_expressed=False for CellTypist compatibility
        sc.pp.normalize_total(
            adata, target_sum=target_sum, exclude_highly_expressed=False
        )
        sc.pp.log1p(adata)

        logger.info(f"Applied normalize_total({target_sum:.0f}) + log1p")
        return adata

    # Path 3: No raw counts and X is not log1p_10k - unsafe territory
    if not status["is_usable"]:
        error_msg = (
            f"Cannot safely normalize data.\n"
            f"  Detected X state: {status['x_state']}\n"
            f"  Estimated target_sum: {status.get('x_target_sum', 'N/A')}\n"
            f"  Raw counts found: None\n"
            f"\n"
            f"To resolve this:\n"
            f"  1. Provide raw counts in adata.layers['counts'] or adata.raw.X\n"
            f"  2. Ensure adata.X contains log1p(10k) normalized data\n"
            f"  3. Use unsafe_force=True if you have manually verified your data\n"
        )

        if unsafe_force:
            logger.warning("=" * 60)
            logger.warning("UNSAFE NORMALIZATION FORCED")
            logger.warning(f"Detected X state: {status['x_state']}")
            logger.warning(f"Estimated target_sum: {status.get('x_target_sum', 'N/A')}")
            logger.warning("This may produce INCORRECT downstream results.")
            logger.warning("You have been warned.")
            logger.warning("=" * 60)

            # Apply full pipeline anyway
            sc.pp.normalize_total(
                adata, target_sum=target_sum, exclude_highly_expressed=False
            )
            sc.pp.log1p(adata)
            logger.warning("Applied normalize_total + log1p on unverified data")
            return adata
        else:
            raise ValueError(error_msg)

    # Should not reach here, but handle gracefully
    logger.warning(f"Unexpected state: {status}")
    return adata


def get_loading_summary(adata: ad.AnnData) -> Dict[str, Any]:
    """
    Get summary statistics for loaded AnnData.

    Parameters
    ----------
    adata : AnnData
        Loaded AnnData object.

    Returns
    -------
    Dict[str, Any]
        Summary statistics including cell/gene counts, memory usage, etc.
    """
    from scipy.sparse import issparse

    summary = {
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
        "is_sparse": issparse(adata.X),
        "dtype": str(adata.X.dtype),
    }

    # Memory estimate
    if issparse(adata.X):
        summary["matrix_memory_mb"] = adata.X.data.nbytes / (1024**2)
        summary["sparsity"] = 1 - (adata.X.nnz / (adata.n_obs * adata.n_vars))
    else:
        summary["matrix_memory_mb"] = adata.X.nbytes / (1024**2)
        summary["sparsity"] = 0.0

    # Cell type info if available
    for col in ["cell_type", "celltype", "CellType"]:
        if col in adata.obs.columns:
            summary["cell_type_column"] = col
            summary["n_cell_types"] = adata.obs[col].nunique()
            break

    return summary
