"""
General utilities for SpatialCore.

This module provides cross-cutting utilities used throughout the codebase:
- Gene ID mapping (Ensembl → HUGO/HGNC symbols)
- Expression normalization status detection

These functions are not specific to any data source and can be used
with any AnnData object.
"""

from pathlib import Path
from typing import Dict, Optional, Tuple, Union, Any

import numpy as np
import pandas as pd
import anndata as ad

from spatialcore.core.logging import get_logger

logger = get_logger(__name__)


# ============================================================================
# Ensembl to HUGO Gene Mapping
# ============================================================================

BIOMART_URL = "http://www.ensembl.org/biomart/martservice"

BIOMART_QUERY_TEMPLATE = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1" count="" datasetConfigVersion="0.6">
    <Dataset name="hsapiens_gene_ensembl" interface="default">
        <Attribute name="ensembl_gene_id"/>
        <Attribute name="hgnc_symbol"/>
        <Attribute name="external_gene_name"/>
    </Dataset>
</Query>"""


def download_ensembl_mapping(
    output_path: Union[str, Path],
    force: bool = False,
    timeout: float = 30.0,
) -> Path:
    """
    Download Ensembl-to-HUGO gene mapping from BioMart.

    Downloads a TSV file mapping Ensembl gene IDs (ENSG...) to HUGO/HGNC gene symbols.

    Parameters
    ----------
    output_path : str or Path
        Path to save the mapping TSV file.
    force : bool, default False
        If True, re-download even if file exists.
    timeout : float, default 30.0
        Timeout in seconds for the BioMart download request.

    Returns
    -------
    Path
        Path to the downloaded TSV file.

    Notes
    -----
    The downloaded file contains columns:
    - Gene stable ID (Ensembl ID)
    - HGNC symbol
    - Gene name (external gene name)

    Examples
    --------
    >>> from spatialcore.core.utils import download_ensembl_mapping
    >>> path = download_ensembl_mapping("./cache/ensembl_to_hugo.tsv")
    """
    import shutil
    import socket
    import urllib.request
    import urllib.parse
    import urllib.error

    output_path = Path(output_path)

    if output_path.exists() and not force:
        logger.info(f"Mapping file already exists: {output_path}")
        return output_path

    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info("Downloading Ensembl-to-HUGO gene mapping from BioMart...")

    # Encode the query
    query = urllib.parse.quote(BIOMART_QUERY_TEMPLATE)
    url = f"{BIOMART_URL}?query={query}"

    try:
        with urllib.request.urlopen(url, timeout=timeout) as response, open(output_path, "wb") as f:
            shutil.copyfileobj(response, f)

        # Verify the download
        df = pd.read_csv(output_path, sep="\t")
        n_mappings = len(df[df["HGNC symbol"].notna() & (df["HGNC symbol"] != "")])
        logger.info(f"Downloaded {n_mappings:,} gene mappings to: {output_path}")

        return output_path

    except (urllib.error.URLError, socket.timeout, TimeoutError) as e:
        logger.error(f"Failed to download from BioMart (timeout={timeout}s): {e}")
        raise
    except Exception as e:
        logger.error(f"Failed to download from BioMart: {e}")
        raise


def load_ensembl_to_hugo_mapping(
    cache_path: Optional[Union[str, Path]] = None,
    auto_download: bool = True,
) -> Dict[str, str]:
    """
    Load Ensembl ID to HUGO gene symbol mapping.

    Parameters
    ----------
    cache_path : str or Path, optional
        Path to cached TSV file. If None, uses default cache location.
    auto_download : bool, default True
        If True and cache doesn't exist, download from BioMart.

    Returns
    -------
    Dict[str, str]
        Mapping from Ensembl ID (e.g., "ENSG00000141510") to HUGO symbol (e.g., "TP53").

    Examples
    --------
    >>> from spatialcore.core.utils import load_ensembl_to_hugo_mapping
    >>> mapping = load_ensembl_to_hugo_mapping()
    >>> mapping["ENSG00000141510"]
    'TP53'
    """
    if cache_path is None:
        # Default cache location
        cache_path = Path.home() / ".cache" / "spatialcore" / "ensembl_to_hugo.tsv"
    else:
        cache_path = Path(cache_path)

    if not cache_path.exists():
        if auto_download:
            download_ensembl_mapping(cache_path)
        else:
            raise FileNotFoundError(
                f"Mapping file not found: {cache_path}. "
                "Set auto_download=True to download from BioMart."
            )

    df = pd.read_csv(cache_path, sep="\t")

    # Filter for valid HGNC symbols
    df = df.dropna(subset=["HGNC symbol"])
    df = df[df["HGNC symbol"].str.len() > 0]

    # Create mapping (Ensembl ID -> HGNC symbol)
    mapping = dict(zip(df["Gene stable ID"], df["HGNC symbol"]))

    logger.info(f"Loaded {len(mapping):,} Ensembl to HUGO gene mappings")
    return mapping


def is_ensembl_id(gene_name: str) -> bool:
    """
    Check if a gene name looks like an Ensembl ID.

    Parameters
    ----------
    gene_name : str
        Gene name to check.

    Returns
    -------
    bool
        True if the name matches Ensembl ID patterns.

    Examples
    --------
    >>> from spatialcore.core.utils import is_ensembl_id
    >>> is_ensembl_id("ENSG00000141510")
    True
    >>> is_ensembl_id("TP53")
    False
    >>> is_ensembl_id("ENSMUSG00000059552")
    True
    """
    if not gene_name or not isinstance(gene_name, str):
        return False
    return (
        gene_name.startswith("ENSG") or      # Human gene
        gene_name.startswith("ENST") or      # Human transcript
        gene_name.startswith("ENSMUSG") or   # Mouse gene
        gene_name.startswith("ENSMUS")       # Mouse
    )


def _convert_ensembl_to_hugo(
    gene_names: np.ndarray,
    ensembl_to_hugo: Dict[str, str],
) -> Tuple[np.ndarray, Dict[str, int]]:
    """
    Convert Ensembl IDs to HUGO gene symbols where possible.

    Safe for HUGO symbols: passes them through unchanged.

    Parameters
    ----------
    gene_names : np.ndarray
        Array of gene names.
    ensembl_to_hugo : Dict[str, str]
        Mapping from Ensembl ID to HUGO symbol.

    Returns
    -------
    Tuple[np.ndarray, Dict[str, int]]
        (converted_names, stats_dict)
    """
    converted = []
    n_converted = 0
    n_already_hugo = 0
    n_unmapped = 0

    for gene in gene_names:
        gene_str = str(gene)
        if is_ensembl_id(gene_str):
            if gene_str in ensembl_to_hugo:
                converted.append(ensembl_to_hugo[gene_str])
                n_converted += 1
            else:
                # Keep unmapped Ensembl ID (will be filtered during panel subsetting)
                converted.append(gene_str)
                n_unmapped += 1
        else:
            # Already HUGO symbol - pass through
            converted.append(gene_str)
            n_already_hugo += 1

    stats = {
        "total_genes": len(gene_names),
        "converted_ensembl": n_converted,
        "already_hugo": n_already_hugo,
        "unmapped_ensembl": n_unmapped,
    }
    return np.array(converted), stats


def _normalize_var_names(
    var_names: pd.Index,
    var_df: pd.DataFrame,
    ensembl_to_hugo: Dict[str, str],
) -> Tuple[np.ndarray, Dict[str, int], bool, bool]:
    """
    Normalize var_names using feature_name and Ensembl -> HUGO mapping.

    Returns converted names, conversion stats, and flags indicating
    whether non-symbol IDs were detected and feature_name was used.
    """
    first_gene = str(var_names[0])
    uses_non_symbol_ids = (
        first_gene.isdigit() or
        first_gene.startswith("ENSG") or
        first_gene.startswith("ENST")
    )

    base_names = var_names.values
    used_feature_name = False
    if uses_non_symbol_ids and "feature_name" in var_df.columns:
        base_names = var_df["feature_name"].values.astype(str)
        used_feature_name = True

    converted_names, stats = _convert_ensembl_to_hugo(
        np.asarray(base_names), ensembl_to_hugo
    )
    return converted_names, stats, uses_non_symbol_ids, used_feature_name


def normalize_gene_names(
    adata: ad.AnnData,
    ensembl_to_hugo: Optional[Dict[str, str]] = None,
    copy: bool = False,
) -> ad.AnnData:
    """
    Normalize gene names in AnnData to use HUGO gene symbols.

    Two-step process:
    1. If var_names are Ensembl IDs/indices, use feature_name column as starting point
    2. Apply Ensembl→HUGO mapping for any remaining Ensembl IDs

    Safe to call on data that already uses HUGO symbols.

    Parameters
    ----------
    adata : AnnData
        AnnData object with genes in var.
    ensembl_to_hugo : Dict[str, str], optional
        Mapping from Ensembl ID to HUGO symbol. If None, loads from cache.
    copy : bool, default False
        If True, return a copy. Otherwise modify in place.

    Returns
    -------
    AnnData
        AnnData with normalized gene names in var_names.
        If adata.raw is present, its var_names are updated to stay aligned.

    Notes
    -----
    CellxGene Census data commonly stores gene identifiers as:
    - Numeric indices with symbols in var['feature_name']
    - Ensembl IDs with symbols in var['feature_name']
    - Mixed content in feature_name (some Ensembl, some HUGO)

    This function handles all these cases.

    Examples
    --------
    >>> from spatialcore.core.utils import normalize_gene_names, load_ensembl_to_hugo_mapping
    >>> mapping = load_ensembl_to_hugo_mapping()
    >>> adata = normalize_gene_names(adata, mapping)
    """
    if copy:
        adata = adata.copy()

    if ensembl_to_hugo is None:
        ensembl_to_hugo = load_ensembl_to_hugo_mapping()

    converted_names, stats, uses_non_symbol_ids, used_feature_name = _normalize_var_names(
        adata.var_names, adata.var, ensembl_to_hugo
    )

    if not uses_non_symbol_ids:
        logger.info("Gene names already appear to be HUGO symbols")
        if stats["converted_ensembl"] > 0:
            adata.var_names = pd.Index(converted_names)
            adata.var_names_make_unique()
            logger.info(
                f"Converted {stats['converted_ensembl']:,} remaining Ensembl IDs to HUGO"
            )
        if stats["unmapped_ensembl"] > 0:
            logger.warning(
                f"{stats['unmapped_ensembl']:,} Ensembl IDs not found in mapping; "
                "leaving them unchanged"
            )
    else:
        if used_feature_name:
            logger.info("Using 'feature_name' column as gene names")

        adata.var_names = pd.Index(converted_names)

        if stats["converted_ensembl"] > 0 or stats["unmapped_ensembl"] > 0:
            logger.info(
                f"Gene mapping: {stats['converted_ensembl']:,} converted, "
                f"{stats['already_hugo']:,} already HUGO, "
                f"{stats['unmapped_ensembl']:,} unmapped"
            )
            if stats["unmapped_ensembl"] > 0:
                logger.warning(
                    f"{stats['unmapped_ensembl']:,} Ensembl IDs not found in mapping; "
                    "leaving them unchanged"
                )
        else:
            logger.info(f"All {stats['already_hugo']:,} genes already HUGO symbols")

        adata.var_names_make_unique()

    if adata.raw is not None:
        raw_converted, raw_stats, _, raw_used_feature = _normalize_var_names(
            adata.raw.var_names, adata.raw.var, ensembl_to_hugo
        )
        raw_converted_index = pd.Index(raw_converted)

        if raw_used_feature or not raw_converted_index.equals(adata.raw.var_names):
            raw_adata = adata.raw.to_adata()
            raw_adata.var_names = raw_converted_index
            raw_adata.var_names_make_unique()
            adata.raw = raw_adata
            logger.info("Updated adata.raw.var_names to normalized HUGO symbols")
            if raw_stats["unmapped_ensembl"] > 0:
                logger.warning(
                    f"{raw_stats['unmapped_ensembl']:,} raw Ensembl IDs not found in mapping; "
                    "leaving them unchanged"
                )

    return adata


# ============================================================================
# Expression Normalization Status Detection
# ============================================================================

# Layer names to search for raw counts (in priority order)
RAW_COUNT_LAYERS = ["counts", "raw_counts", "raw"]


def _is_integer_like(
    values: np.ndarray,
    tolerance: float = 1e-6,
    threshold: float = 0.95,
) -> bool:
    """
    Check if array values are integer-like within floating point tolerance.

    Parameters
    ----------
    values : np.ndarray
        Array of values to check (should be non-zero values only).
    tolerance : float, default 1e-6
        Tolerance for floating point comparison. Handles values like
        1.0000000000000002 or 2.9999999999999996.
    threshold : float, default 0.95
        Fraction of values that must be integer-like to pass.

    Returns
    -------
    bool
        True if >= threshold fraction of values are integer-like.
    """
    if len(values) == 0:
        return False

    remainder = np.abs(np.mod(values, 1))
    is_integer = (remainder < tolerance) | (remainder > 1 - tolerance)
    fraction_integer = np.mean(is_integer)

    return fraction_integer >= threshold


def _get_matrix_sample(
    matrix,
    sample_size: int = 10000,
) -> np.ndarray:
    """
    Get a dense sample from a matrix (sparse or dense).

    Parameters
    ----------
    matrix : array-like or sparse matrix
        Expression matrix.
    sample_size : int, default 10000
        Maximum number of cells to sample.

    Returns
    -------
    np.ndarray
        Dense 2D array sample.
    """
    from scipy.sparse import issparse

    n_cells = matrix.shape[0]
    n_sample = min(sample_size, n_cells)

    if issparse(matrix):
        return matrix[:n_sample].toarray()
    else:
        return np.asarray(matrix[:n_sample])


def _check_raw_counts(
    matrix,
    sample_size: int = 10000,
    integer_tolerance: float = 1e-6,
    integer_threshold: float = 0.95,
) -> Dict[str, Any]:
    """
    Check if a matrix contains raw counts.

    Parameters
    ----------
    matrix : array-like or sparse matrix
        Expression matrix to check.
    sample_size : int, default 10000
        Number of cells to sample for checking.
    integer_tolerance : float, default 1e-6
        Tolerance for integer comparison.
    integer_threshold : float, default 0.95
        Fraction of values that must be integers.

    Returns
    -------
    Dict[str, Any]
        Dictionary with:
        - is_raw: bool
        - fraction_integer: float
        - min_val: float
        - max_val: float
    """
    sample_data = _get_matrix_sample(matrix, sample_size)

    # Get non-zero values for integer check
    flat_sample = sample_data.flatten()
    non_zero = flat_sample[flat_sample != 0]

    if len(non_zero) == 0:
        return {
            "is_raw": False,
            "fraction_integer": 0.0,
            "min_val": 0.0,
            "max_val": 0.0,
            "reason": "all_zeros",
        }

    min_val = float(np.min(sample_data))
    max_val = float(np.max(sample_data))

    # Raw counts cannot be negative
    if min_val < 0:
        return {
            "is_raw": False,
            "fraction_integer": 0.0,
            "min_val": min_val,
            "max_val": max_val,
            "reason": "negative_values",
        }

    # Check integer-like property
    is_integer = _is_integer_like(non_zero, integer_tolerance, integer_threshold)

    # Calculate actual fraction for reporting
    remainder = np.abs(np.mod(non_zero, 1))
    fraction_integer = float(np.mean(
        (remainder < integer_tolerance) | (remainder > 1 - integer_tolerance)
    ))

    return {
        "is_raw": is_integer,
        "fraction_integer": fraction_integer,
        "min_val": min_val,
        "max_val": max_val,
        "reason": "integer_check",
    }


def _estimate_target_sum(
    matrix,
    sample_size: int = 1000,
) -> Dict[str, Any]:
    """
    Estimate the target sum used for normalization by reversing log1p.

    If data is log1p(counts / total * target_sum), then:
    expm1(X).sum(axis=1) should equal target_sum for each cell.

    Parameters
    ----------
    matrix : array-like or sparse matrix
        Expression matrix (assumed to be log1p transformed).
    sample_size : int, default 1000
        Number of cells to sample.

    Returns
    -------
    Dict[str, Any]
        Dictionary with:
        - estimated_target_sum: float (median of row sums)
        - target_sum_std: float (std of row sums)
        - is_log1p_10k: bool (True if target_sum ~ 10,000)
        - is_log1p_cpm: bool (True if target_sum ~ 1,000,000)
    """
    sample_data = _get_matrix_sample(matrix, sample_size)

    # Reverse log1p transformation
    reversed_data = np.expm1(sample_data)

    # Compute row sums (should equal target_sum)
    row_sums = reversed_data.sum(axis=1)

    # Exclude empty cells (sum = 0)
    non_empty_sums = row_sums[row_sums > 0]

    if len(non_empty_sums) == 0:
        return {
            "estimated_target_sum": 0.0,
            "target_sum_std": 0.0,
            "is_log1p_10k": False,
            "is_log1p_cpm": False,
        }

    median_sum = float(np.median(non_empty_sums))
    std_sum = float(np.std(non_empty_sums))

    # Check if close to 10,000 (allow 20% tolerance)
    is_log1p_10k = 8_000 < median_sum < 12_000

    # Check if close to 1,000,000 (CPM)
    is_log1p_cpm = 800_000 < median_sum < 1_200_000

    return {
        "estimated_target_sum": median_sum,
        "target_sum_std": std_sum,
        "is_log1p_10k": is_log1p_10k,
        "is_log1p_cpm": is_log1p_cpm,
    }


def _find_raw_counts_source(
    adata: ad.AnnData,
    sample_size: int = 10000,
    integer_tolerance: float = 1e-6,
    integer_threshold: float = 0.95,
) -> Optional[str]:
    """
    Search for raw counts in layers, adata.raw, and adata.X.

    Parameters
    ----------
    adata : AnnData
        AnnData object to search.
    sample_size : int, default 10000
        Number of cells to sample for checking.
    integer_tolerance : float, default 1e-6
        Tolerance for integer comparison.
    integer_threshold : float, default 0.95
        Fraction of values that must be integers.

    Returns
    -------
    Optional[str]
        Source location if found:
        - "layers/{layer_name}" for layers
        - "raw.X" for adata.raw
        - "X" for adata.X
        - None if no raw counts found
    """
    # Check layers first (in priority order)
    for layer_name in RAW_COUNT_LAYERS:
        if layer_name in adata.layers:
            result = _check_raw_counts(
                adata.layers[layer_name],
                sample_size,
                integer_tolerance,
                integer_threshold,
            )
            if result["is_raw"]:
                logger.debug(
                    f"Found raw counts in layers['{layer_name}'] "
                    f"(fraction_integer={result['fraction_integer']:.3f})"
                )
                return f"layers/{layer_name}"

    # Check adata.raw.X
    if adata.raw is not None:
        result = _check_raw_counts(
            adata.raw.X,
            sample_size,
            integer_tolerance,
            integer_threshold,
        )
        if result["is_raw"]:
            logger.debug(
                f"Found raw counts in raw.X "
                f"(fraction_integer={result['fraction_integer']:.3f})"
            )
            return "raw.X"

    # Check adata.X as last resort
    result = _check_raw_counts(
        adata.X,
        sample_size,
        integer_tolerance,
        integer_threshold,
    )
    if result["is_raw"]:
        logger.debug(
            f"Found raw counts in X "
            f"(fraction_integer={result['fraction_integer']:.3f})"
        )
        return "X"

    return None


def check_normalization_status(
    adata: ad.AnnData,
    sample_size: int = 1000,
    integer_tolerance: float = 1e-6,
    integer_threshold: float = 0.95,
) -> Dict[str, Any]:
    """
    Detect the normalization state of expression data with robust validation.

    This function searches for raw counts in layers and adata.raw, and verifies
    log1p normalization by checking the target sum via expm1 reversal.

    Parameters
    ----------
    adata : AnnData
        AnnData object to check.
    sample_size : int, default 1000
        Number of cells to sample for detection.
    integer_tolerance : float, default 1e-6
        Tolerance for integer comparison (handles float precision issues
        like 1.0000000000000002).
    integer_threshold : float, default 0.95
        Fraction of non-zero values that must be integer-like to classify
        as raw counts.

    Returns
    -------
    Dict[str, Any]
        Dictionary with keys:

        - raw_source: str or None
            Location of raw counts ("layers/counts", "raw.X", "X", or None)
        - x_state: str
            State of adata.X: "raw", "log1p_10k", "log1p_cpm", "log1p_other",
            "linear", "negative", "unknown"
        - x_target_sum: float or None
            Estimated target sum if X appears log-transformed
        - is_usable: bool
            True if data can be safely normalized (raw available OR X is log1p_10k)
        - stats: dict
            Diagnostic statistics (mean, max, min, fraction_integer)

    Notes
    -----
    **Detection Logic:**

    1. Search for raw counts in: layers["counts"], layers["raw_counts"],
       layers["raw"], adata.raw.X, adata.X (in order)
    2. For each candidate, check if >95% of non-zero values are integers
       within floating point tolerance
    3. For adata.X, if not raw, reverse log1p and check row sums to
       verify target_sum is ~10,000

    **Usability Criteria:**

    Data is considered usable (is_usable=True) if:
    - Raw counts are found anywhere, OR
    - adata.X is verified as log1p normalized to 10,000

    Examples
    --------
    >>> from spatialcore.core.utils import check_normalization_status
    >>> status = check_normalization_status(adata)
    >>> if status["is_usable"]:
    ...     if status["raw_source"]:
    ...         print(f"Will normalize from {status['raw_source']}")
    ...     else:
    ...         print("X is already log1p_10k")
    >>> else:
    ...     print(f"Cannot use: X is {status['x_state']}")
    """
    from scipy.sparse import issparse

    # Step 1: Search for raw counts
    raw_source = _find_raw_counts_source(
        adata, sample_size * 10, integer_tolerance, integer_threshold
    )

    # Step 2: Analyze adata.X
    sample_data = _get_matrix_sample(adata.X, sample_size)

    # Empty data should fail explicitly
    if sample_data.size == 0:
        raise ValueError("Cannot check normalization status of empty AnnData (0 cells or 0 genes)")

    mean_val = float(np.mean(sample_data))
    max_val = float(np.max(sample_data))
    min_val = float(np.min(sample_data))

    # Check if X contains raw counts
    x_raw_check = _check_raw_counts(
        adata.X, sample_size * 10, integer_tolerance, integer_threshold
    )

    stats = {
        "mean": mean_val,
        "max": max_val,
        "min": min_val,
        "fraction_integer": x_raw_check["fraction_integer"],
    }

    # Determine X state
    if x_raw_check["is_raw"]:
        x_state = "raw"
        x_target_sum = None
    elif min_val < 0:
        x_state = "negative"
        x_target_sum = None
    elif max_val < 25 and mean_val < 10 and min_val >= 0:
        # Likely log-transformed, verify target sum via expm1 reversal
        # log1p(10k) typically has: max ~6-9, mean ~3-6 for typical scRNA-seq
        # log1p(1M) typically has: max ~13-15, mean ~8-12
        target_info = _estimate_target_sum(adata.X, sample_size)
        x_target_sum = target_info["estimated_target_sum"]

        if target_info["is_log1p_10k"]:
            x_state = "log1p_10k"
        elif target_info["is_log1p_cpm"]:
            x_state = "log1p_cpm"
        elif x_target_sum > 0:
            x_state = "log1p_other"
        else:
            x_state = "unknown"

        stats["estimated_target_sum"] = x_target_sum
    elif max_val > 25 and x_raw_check["fraction_integer"] < 0.5:
        x_state = "linear"
        x_target_sum = None
    else:
        x_state = "unknown"
        x_target_sum = None

    # Determine if data is usable
    is_usable = (raw_source is not None) or (x_state == "log1p_10k")

    # Check for scanpy log1p metadata as additional confirmation
    has_log1p_uns = "log1p" in adata.uns

    return {
        "raw_source": raw_source,
        "x_state": x_state,
        "x_target_sum": x_target_sum,
        "is_usable": is_usable,
        "has_log1p_uns": has_log1p_uns,
        "stats": stats,
    }
