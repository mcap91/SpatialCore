"""
Confidence score transformation, filtering, and decision score storage.

This module provides utilities for:
1. Transforming CellTypist decision scores to meaningful confidence values
2. Storing decision score matrices in AnnData for downstream analysis
3. Filtering cells by confidence or cell type count thresholds

For spatial transcriptomics, raw CellTypist confidence values may be less
informative than z-score transformed values, which capture how confident
a prediction is relative to other cell types.

References:
    - CellTypist: https://www.celltypist.org/
    - Domínguez Conde et al., Science (2022)
"""

from typing import Dict, List, Literal, Optional, Tuple

import numpy as np
import pandas as pd
import anndata as ad
from scipy.special import softmax as scipy_softmax

from spatialcore.core.logging import get_logger

logger = get_logger(__name__)


# ============================================================================
# Confidence Transformation
# ============================================================================

ConfidenceMethod = Literal["raw", "zscore", "softmax", "minmax"]


def transform_confidence(
    decision_scores: np.ndarray,
    method: ConfidenceMethod = "zscore",
) -> np.ndarray:
    """
    Transform CellTypist decision scores to meaningful confidence values.

    CellTypist produces logistic regression decision scores which can be
    negative and unbounded. This function transforms them to interpretable
    [0, 1] confidence values using different strategies.

    Parameters
    ----------
    decision_scores : np.ndarray
        Decision score matrix of shape (n_cells, n_types) from CellTypist.
        Each row contains scores for all cell types for one cell.
    method : {"raw", "zscore", "softmax", "minmax"}, default "zscore"
        Transformation method:

        - "raw": Return winning score directly (may be negative/unbounded)
        - "zscore": Sigmoid of z-score of winning type vs all types.
          Recommended for spatial data. Captures how "distinct" the prediction is.
        - "softmax": Softmax probability of winning type.
          Sums to 1 across types, good for comparing type probabilities.
        - "minmax": Min-max scaling to [0, 1] per cell.
          Simple but may not be well-calibrated.

    Returns
    -------
    np.ndarray
        Array of shape (n_cells,) with transformed confidence values.
        For "zscore", "softmax", "minmax": values in [0, 1].
        For "raw": unbounded values.

    Notes
    -----
    **Why z-score for spatial data?**

    CellTypist decision scores from logistic regression can be negative and
    don't have a natural scale. The z-score method computes how many standard
    deviations the winning type's score is above the mean, then applies
    sigmoid to get a [0, 1] value:

        confidence = sigmoid((winning_score - mean) / std)

    This is more informative than raw confidence because:
    1. A cell with scores [5.0, 0.1, 0.1, 0.1] has high z-score (clear winner)
    2. A cell with scores [5.0, 4.8, 4.7, 4.6] has low z-score (ambiguous)

    **Method comparison:**

    +----------+------------------+------------------+
    | Method   | Best For         | Output Range     |
    +==========+==================+==================+
    | zscore   | Spatial data     | [0, 1] (sigmoid) |
    | softmax  | Probability est. | [0, 1] (prob)    |
    | minmax   | Simple scaling   | [0, 1] (linear)  |
    | raw      | Debug/analysis   | unbounded        |
    +----------+------------------+------------------+

    Examples
    --------
    >>> from spatialcore.annotation.confidence import transform_confidence
    >>> import numpy as np
    >>> # Decision scores from CellTypist (n_cells=3, n_types=4)
    >>> scores = np.array([
    ...     [5.0, 0.1, 0.1, 0.1],  # Clear winner
    ...     [2.0, 1.8, 1.9, 1.7],  # Ambiguous
    ...     [0.5, -1.0, -0.5, 0.3],  # Negative scores
    ... ])
    >>> conf = transform_confidence(scores, method="zscore")
    >>> print(f"Clear winner: {conf[0]:.3f}")  # High confidence
    >>> print(f"Ambiguous: {conf[1]:.3f}")      # Low confidence
    """
    if decision_scores.ndim != 2:
        raise ValueError(
            f"Expected 2D array of shape (n_cells, n_types), "
            f"got shape {decision_scores.shape}"
        )

    n_cells, n_types = decision_scores.shape

    if n_types < 2:
        raise ValueError(
            f"Expected at least 2 cell types, got {n_types}"
        )

    # Get winning type index for each cell
    winning_idx = np.argmax(decision_scores, axis=1)
    winning_scores = decision_scores[np.arange(n_cells), winning_idx]

    if method == "raw":
        return winning_scores

    elif method == "zscore":
        # Z-score: how many std above mean is the winning score?
        mean_scores = np.mean(decision_scores, axis=1)
        std_scores = np.std(decision_scores, axis=1)

        # Avoid division by zero (all types have same score)
        std_scores = np.where(std_scores < 1e-10, 1.0, std_scores)

        z_scores = (winning_scores - mean_scores) / std_scores

        # Sigmoid to [0, 1]
        confidence = 1 / (1 + np.exp(-z_scores))
        return confidence

    elif method == "softmax":
        # Softmax probability of winning type
        probs = scipy_softmax(decision_scores, axis=1)
        confidence = probs[np.arange(n_cells), winning_idx]
        return confidence

    elif method == "minmax":
        # Min-max scaling per cell
        min_scores = np.min(decision_scores, axis=1, keepdims=True)
        max_scores = np.max(decision_scores, axis=1, keepdims=True)

        # Avoid division by zero
        score_range = max_scores - min_scores
        score_range = np.where(score_range < 1e-10, 1.0, score_range)

        scaled = (decision_scores - min_scores) / score_range
        confidence = scaled[np.arange(n_cells), winning_idx]
        return confidence

    else:
        raise ValueError(
            f"Unknown confidence method: {method}. "
            f"Expected one of: raw, zscore, softmax, minmax"
        )


# ============================================================================
# Decision Score Storage
# ============================================================================

def extract_decision_scores(
    adata: ad.AnnData,
    celltypist_result,
    key_added: str = "celltypist",
    copy: bool = False,
) -> ad.AnnData:
    """
    Store CellTypist decision scores matrix in AnnData.

    Extracts the full decision score matrix from CellTypist annotation
    results and stores it in adata.obsm for downstream analysis
    (e.g., confidence transforms, plotting, uncertainty analysis).

    Parameters
    ----------
    adata : AnnData
        AnnData object that was annotated with CellTypist.
    celltypist_result
        CellTypist AnnotationResult object from celltypist.annotate().
        Must have .decision_matrix and .cell_types attributes.
    key_added : str, default "celltypist"
        Key prefix for stored results:

        - adata.obsm[f"{key_added}_decision_scores"]: Decision matrix
        - adata.uns[f"{key_added}_cell_types"]: Cell type names
    copy : bool, default False
        Whether to return a copy or modify in-place.

    Returns
    -------
    AnnData
        AnnData with decision scores stored in obsm.

    Notes
    -----
    The decision score matrix has shape (n_cells, n_types) where each
    row contains the logistic regression decision scores for all cell
    types. Higher scores indicate stronger evidence for that cell type.

    These scores can be used for:
    - Custom confidence calculations (transform_confidence)
    - Uncertainty visualization (plotting multiple high-scoring types)
    - Ensemble methods (combining multiple model predictions)

    Examples
    --------
    >>> import celltypist
    >>> from spatialcore.annotation.confidence import extract_decision_scores
    >>> # Run CellTypist annotation
    >>> result = celltypist.annotate(adata, model=model)
    >>> # Store decision scores
    >>> adata = extract_decision_scores(adata, result, key_added="ct")
    >>> # Access stored scores
    >>> scores = adata.obsm["ct_decision_scores"]
    >>> cell_types = adata.uns["ct_cell_types"]
    """
    if copy:
        adata = adata.copy()

    # Extract decision matrix from CellTypist result
    if not hasattr(celltypist_result, "decision_matrix"):
        raise ValueError(
            "CellTypist result does not have decision_matrix. "
            "Ensure you passed the result from celltypist.annotate()."
        )

    decision_matrix = celltypist_result.decision_matrix

    # Handle DataFrame (from CellTypist) or numpy array
    if isinstance(decision_matrix, pd.DataFrame):
        cell_types = list(decision_matrix.columns)
        decision_array = decision_matrix.values
    else:
        # Numpy array - get cell types from result
        decision_array = decision_matrix
        if hasattr(celltypist_result, "cell_types"):
            cell_types = list(celltypist_result.cell_types)
        else:
            # Fallback: generate generic names
            n_types = decision_array.shape[1]
            cell_types = [f"type_{i}" for i in range(n_types)]
            logger.warning(
                f"Could not extract cell type names from result. "
                f"Using generic names: type_0, type_1, ..."
            )

    # Validate shape
    if decision_array.shape[0] != adata.n_obs:
        raise ValueError(
            f"Decision matrix has {decision_array.shape[0]} cells, "
            f"but AnnData has {adata.n_obs} cells. "
            f"Ensure the CellTypist result matches the AnnData object."
        )

    # Store in AnnData
    adata.obsm[f"{key_added}_decision_scores"] = decision_array.astype(np.float32)
    adata.uns[f"{key_added}_cell_types"] = cell_types

    logger.info(
        f"Stored decision scores: {decision_array.shape[0]:,} cells × "
        f"{decision_array.shape[1]} types in adata.obsm['{key_added}_decision_scores']"
    )

    return adata


# ============================================================================
# Confidence Filtering
# ============================================================================

def filter_low_confidence(
    adata: ad.AnnData,
    label_column: str,
    confidence_column: str,
    threshold: float = 0.5,
    unassigned_label: str = "Unassigned",
    copy: bool = False,
) -> ad.AnnData:
    """
    Mark cells below confidence threshold as Unassigned.

    Cells with confidence values below the threshold have their cell type
    label replaced with an unassigned label. This is useful for quality
    control to flag uncertain predictions.

    Parameters
    ----------
    adata : AnnData
        AnnData object with cell type labels and confidence values.
    label_column : str
        Column in adata.obs containing cell type labels.
    confidence_column : str
        Column in adata.obs containing confidence values.
    threshold : float, default 0.5
        Confidence threshold. Cells with confidence < threshold are marked
        as unassigned.
    unassigned_label : str, default "Unassigned"
        Label to assign to low-confidence cells.
    copy : bool, default False
        Whether to return a copy or modify in-place.

    Returns
    -------
    AnnData
        AnnData with low-confidence cells marked as unassigned.

    Notes
    -----
    The original labels are preserved if you use copy=True. For tracking
    purposes, consider storing the original labels in a separate column
    before filtering.

    Examples
    --------
    >>> from spatialcore.annotation.confidence import filter_low_confidence
    >>> # Filter cells with confidence < 0.6
    >>> adata = filter_low_confidence(
    ...     adata,
    ...     label_column="celltypist_prediction",
    ...     confidence_column="celltypist_confidence",
    ...     threshold=0.6,
    ... )
    >>> # Check how many were marked as unassigned
    >>> n_unassigned = (adata.obs["celltypist_prediction"] == "Unassigned").sum()
    >>> print(f"Marked {n_unassigned:,} cells as Unassigned")
    """
    if copy:
        adata = adata.copy()

    # Validate columns exist
    if label_column not in adata.obs.columns:
        raise ValueError(
            f"Label column '{label_column}' not found in adata.obs. "
            f"Available: {list(adata.obs.columns)}"
        )
    if confidence_column not in adata.obs.columns:
        raise ValueError(
            f"Confidence column '{confidence_column}' not found in adata.obs. "
            f"Available: {list(adata.obs.columns)}"
        )

    # Get confidence values
    confidence = adata.obs[confidence_column].values

    # Mark low confidence cells
    low_conf_mask = confidence < threshold
    n_low_conf = low_conf_mask.sum()

    if n_low_conf > 0:
        # Ensure label column is string type for modification
        adata.obs[label_column] = adata.obs[label_column].astype(str)
        adata.obs.loc[low_conf_mask, label_column] = unassigned_label

        pct = 100 * n_low_conf / adata.n_obs
        logger.info(
            f"Marked {n_low_conf:,} cells ({pct:.1f}%) as '{unassigned_label}' "
            f"(confidence < {threshold})"
        )
    else:
        logger.info(f"No cells below confidence threshold {threshold}")

    return adata


def filter_low_count_types(
    adata: ad.AnnData,
    label_column: str,
    min_cells: int = 15,
    unassigned_label: str = "Low_count",
    copy: bool = False,
) -> ad.AnnData:
    """
    Mark cell types with fewer than min_cells as Low_count.

    Cell types with very few cells may be unreliable annotations or
    artifacts. This function marks cells of rare types with a special
    label for downstream filtering or analysis.

    Parameters
    ----------
    adata : AnnData
        AnnData object with cell type labels.
    label_column : str
        Column in adata.obs containing cell type labels.
    min_cells : int, default 15
        Minimum cells required for a cell type. Types with fewer cells
        are marked with unassigned_label.
    unassigned_label : str, default "Low_count"
        Label to assign to cells of rare types.
    copy : bool, default False
        Whether to return a copy or modify in-place.

    Returns
    -------
    AnnData
        AnnData with rare cell types marked.

    Notes
    -----
    This is different from filter_low_confidence:

    - filter_low_confidence: Marks individual cells with low prediction confidence
    - filter_low_count_types: Marks entire cell types that have too few members

    A cell type might have high individual confidence but still be rare
    in the dataset (e.g., 5 cells all with 0.9 confidence).

    Examples
    --------
    >>> from spatialcore.annotation.confidence import filter_low_count_types
    >>> # Mark cell types with fewer than 20 cells
    >>> adata = filter_low_count_types(
    ...     adata,
    ...     label_column="celltypist_prediction",
    ...     min_cells=20,
    ...     unassigned_label="Rare_type",
    ... )
    >>> # Check which types were affected
    >>> print(adata.obs["celltypist_prediction"].value_counts())
    """
    if copy:
        adata = adata.copy()

    # Validate column exists
    if label_column not in adata.obs.columns:
        raise ValueError(
            f"Label column '{label_column}' not found in adata.obs. "
            f"Available: {list(adata.obs.columns)}"
        )

    # Count cells per type
    type_counts = adata.obs[label_column].value_counts()
    low_count_types = type_counts[type_counts < min_cells].index.tolist()

    if len(low_count_types) > 0:
        # Create mask for cells of low-count types
        low_count_mask = adata.obs[label_column].isin(low_count_types)
        n_affected = low_count_mask.sum()

        # Ensure label column is string type
        adata.obs[label_column] = adata.obs[label_column].astype(str)
        adata.obs.loc[low_count_mask, label_column] = unassigned_label

        pct = 100 * n_affected / adata.n_obs
        logger.info(
            f"Marked {n_affected:,} cells ({pct:.1f}%) from "
            f"{len(low_count_types)} rare types as '{unassigned_label}' "
            f"(types with < {min_cells} cells)"
        )
        logger.info(f"  Affected types: {low_count_types[:5]}{'...' if len(low_count_types) > 5 else ''}")
    else:
        logger.info(f"No cell types with fewer than {min_cells} cells")

    return adata


def compute_confidence_from_obsm(
    adata: ad.AnnData,
    decision_scores_key: str = "cell_type_decision_scores",
    method: ConfidenceMethod = "zscore",
    confidence_column: str = "confidence_transformed",
    copy: bool = False,
) -> ad.AnnData:
    """
    Compute transformed confidence from stored decision scores.

    Convenience function that reads decision scores from adata.obsm
    and applies transform_confidence, storing the result in adata.obs.

    Parameters
    ----------
    adata : AnnData
        AnnData with decision scores in obsm.
    decision_scores_key : str, default "cell_type_decision_scores"
        Key in adata.obsm containing decision score matrix.
    method : {"raw", "zscore", "softmax", "minmax"}, default "zscore"
        Transformation method (see transform_confidence).
    confidence_column : str, default "confidence_transformed"
        Column name in adata.obs for output confidence values.
    copy : bool, default False
        Whether to return a copy or modify in-place.

    Returns
    -------
    AnnData
        AnnData with confidence values in adata.obs[confidence_column].

    Examples
    --------
    >>> from spatialcore.annotation.confidence import compute_confidence_from_obsm
    >>> # Assuming decision scores are already stored
    >>> adata = compute_confidence_from_obsm(adata, method="zscore")
    >>> print(adata.obs["confidence_transformed"].describe())
    """
    if copy:
        adata = adata.copy()

    if decision_scores_key not in adata.obsm:
        raise ValueError(
            f"Decision scores key '{decision_scores_key}' not found in adata.obsm. "
            f"Available: {list(adata.obsm.keys())}. "
            f"Run extract_decision_scores() first."
        )

    decision_scores = adata.obsm[decision_scores_key]
    confidence = transform_confidence(decision_scores, method=method)
    adata.obs[confidence_column] = confidence

    logger.info(
        f"Computed {method} confidence in adata.obs['{confidence_column}'] "
        f"(mean={confidence.mean():.3f}, std={confidence.std():.3f})"
    )

    return adata


# ============================================================================
# Dual-Threshold Marker Validation
# ============================================================================

def filter_by_marker_validation(
    adata: ad.AnnData,
    label_column: str,
    confidence_column: str,
    canonical_markers: Optional[Dict[str, List[str]]] = None,
    confidence_threshold: float = 0.5,
    n_components: int = 3,
    min_cells_per_type: int = 15,
    unassigned_label: str = "Unassigned",
    copy: bool = False,
) -> Tuple[ad.AnnData, pd.DataFrame]:
    """
    Filter cells using BOTH confidence threshold AND GMM-3 marker validation.

    Implements dual-threshold QC from spec:
    - X-axis: Confidence (z-score transformed)
    - Y-axis: Marker expression score (GMM-3 threshold)

    Cells must pass BOTH thresholds to retain their cell type label.
    Uses GMM-3 thresholding internally for marker expression.

    Parameters
    ----------
    adata : AnnData
        AnnData object with cell type labels and confidence values.
    label_column : str
        Column in adata.obs containing cell type labels.
    confidence_column : str
        Column in adata.obs containing confidence values.
    canonical_markers : Dict[str, List[str]], optional
        Dictionary mapping cell types to marker gene lists.
        If None, uses default CANONICAL_MARKERS from markers module.
    confidence_threshold : float, default 0.5
        Minimum confidence threshold. Cells below this are marked unassigned.
    n_components : int, default 3
        Number of GMM components for marker thresholding.
        3 = trimodal (dropout/moderate/high expression).
    min_cells_per_type : int, default 15
        Minimum cells required to validate a cell type. Types with fewer
        cells are marked unassigned.
    unassigned_label : str, default "Unassigned"
        Label to assign to cells that fail validation.
    copy : bool, default False
        Whether to return a copy or modify in-place.

    Returns
    -------
    Tuple[AnnData, pd.DataFrame]
        - AnnData with filtered labels and validation columns added:
          - `{label_column}_validated`: Final validated labels
          - `marker_score`: Mean expression of canonical markers
          - `marker_passes_gmm`: Whether cell passes GMM marker threshold
          - `confidence_passes`: Whether cell passes confidence threshold
          - `validation_pass`: Whether cell passes both thresholds
        - Summary DataFrame with validation statistics per cell type.

    Notes
    -----
    **Dual-threshold rationale:**

    1. Confidence alone may miss cells that are assigned to wrong types
       with deceptively high confidence.
    2. Marker expression alone may miss cells where marker genes are
       not expressed due to dropout or technical noise.
    3. Combining both axes provides more robust QC.

    **GMM-3 thresholding:**

    For spatial data with dropouts, marker expression is often trimodal:
    - Component 0: Zero/dropout (no expression)
    - Component 1: Moderate expression
    - Component 2: High expression

    The threshold is set at the boundary between component 0 and component 1,
    identifying cells with biologically meaningful marker expression.

    Examples
    --------
    >>> from spatialcore.annotation.confidence import filter_by_marker_validation
    >>> # Filter with default markers
    >>> adata, summary = filter_by_marker_validation(
    ...     adata,
    ...     label_column="celltypist",
    ...     confidence_column="celltypist_confidence_transformed",
    ...     confidence_threshold=0.5,
    ... )
    >>> # Check validation summary
    >>> print(summary[["cell_type", "n_cells", "pct_pass"]])

    See Also
    --------
    spatialcore.stats.classify_by_threshold : GMM-3 thresholding function.
    spatialcore.plotting.validation.plot_2d_validation : 2D validation plot.
    """
    if copy:
        adata = adata.copy()

    # Validate columns exist
    if label_column not in adata.obs.columns:
        raise ValueError(
            f"Label column '{label_column}' not found in adata.obs. "
            f"Available: {list(adata.obs.columns)}"
        )
    if confidence_column not in adata.obs.columns:
        raise ValueError(
            f"Confidence column '{confidence_column}' not found in adata.obs. "
            f"Available: {list(adata.obs.columns)}"
        )

    # Load canonical markers if not provided
    if canonical_markers is None:
        from spatialcore.annotation.markers import load_canonical_markers
        canonical_markers = load_canonical_markers()
        logger.info(f"Loaded canonical markers for {len(canonical_markers)} cell types")

    if not canonical_markers:
        raise ValueError(
            "Canonical markers are required for marker validation but none were provided or loaded."
        )

    # Use internal GMM thresholding on 1D marker scores
    try:
        from spatialcore.stats._thresholding import threshold_gmm
    except ImportError:
        raise ImportError(
            "spatialcore.stats._thresholding is required for GMM-3 validation. "
            "Ensure the stats module is properly installed."
        )

    # Initialize result columns
    n_cells = adata.n_obs
    marker_scores = np.zeros(n_cells)
    marker_passes = np.zeros(n_cells, dtype=bool)
    confidence_passes = adata.obs[confidence_column].values >= confidence_threshold

    # Get unique cell types (excluding already unassigned)
    cell_types = adata.obs[label_column].astype(str).unique()
    cell_types = [ct for ct in cell_types if ct.lower() not in ["unassigned", "unknown", "low_count"]]

    # Validation summary
    summary_rows = []

    for cell_type in cell_types:
        # Get cells of this type
        type_mask = adata.obs[label_column].astype(str) == cell_type
        n_type_cells = type_mask.sum()

        # Skip if too few cells
        if n_type_cells < min_cells_per_type:
            logger.debug(f"Skipping {cell_type}: only {n_type_cells} cells (< {min_cells_per_type})")
            summary_rows.append({
                "cell_type": cell_type,
                "n_cells": n_type_cells,
                "has_markers": False,
                "gmm_threshold": np.nan,
                "n_pass_confidence": (type_mask & confidence_passes).sum(),
                "n_pass_marker": 0,
                "n_pass_both": 0,
                "pct_pass": 0.0,
            })
            continue

        # Find matching markers
        markers_for_type = canonical_markers.get(cell_type, [])

        # Try case-insensitive match if exact match fails
        if not markers_for_type:
            for marker_type, markers in canonical_markers.items():
                if marker_type.lower() == cell_type.lower():
                    markers_for_type = markers
                    break

        # Get marker genes that exist in data
        available_markers = [m for m in markers_for_type if m in adata.var_names]

        if not available_markers:
            logger.warning(f"Skipping marker validation for {cell_type}: no markers available in data")
            marker_scores[type_mask] = np.nan
            marker_passes[type_mask] = True
            n_pass_conf = (type_mask & confidence_passes).sum()
            summary_rows.append({
                "cell_type": cell_type,
                "n_cells": n_type_cells,
                "has_markers": False,
                "gmm_threshold": np.nan,
                "n_pass_confidence": n_pass_conf,
                "n_pass_marker": n_type_cells,
                "n_pass_both": n_pass_conf,
                "pct_pass": 100 * n_pass_conf / n_type_cells if n_type_cells > 0 else 0,
            })
            continue

        # Calculate mean marker expression for cells of this type
        # Use integer indexing to avoid anndata 0.12.x boolean mask bug
        mask_indices = np.where(type_mask)[0]
        gene_indices = [adata.var_names.get_loc(g) for g in available_markers]
        expr_matrix = adata.X[mask_indices][:, gene_indices]
        if hasattr(expr_matrix, "toarray"):
            expr_matrix = expr_matrix.toarray()

        mean_marker_expr = np.mean(expr_matrix, axis=1)
        marker_scores[type_mask] = mean_marker_expr

        # Fit GMM-3 threshold on marker expression
        try:
            gmm_threshold, _, _, _ = threshold_gmm(
                mean_marker_expr,
                probability_cutoff=0.3,
                n_components=n_components,
                random_state=42,
            )
            type_marker_passes = mean_marker_expr >= gmm_threshold
            marker_passes[type_mask] = type_marker_passes
        except Exception as e:
            logger.warning(f"GMM fitting failed for {cell_type}: {e}. Skipping marker validation.")
            marker_scores[type_mask] = np.nan
            marker_passes[type_mask] = True
            gmm_threshold = np.nan

        # Calculate summary stats
        n_pass_conf = (type_mask & confidence_passes).sum()
        n_pass_marker = marker_passes[type_mask].sum()
        n_pass_both = (type_mask & confidence_passes & marker_passes).sum()

        summary_rows.append({
            "cell_type": cell_type,
            "n_cells": n_type_cells,
            "has_markers": True,
            "n_markers": len(available_markers),
            "markers": available_markers[:3],  # First 3 for display
            "gmm_threshold": gmm_threshold,
            "n_pass_confidence": n_pass_conf,
            "n_pass_marker": n_pass_marker,
            "n_pass_both": n_pass_both,
            "pct_pass": 100 * n_pass_both / n_type_cells if n_type_cells > 0 else 0,
        })

    # Store validation results in adata
    adata.obs["marker_score"] = marker_scores
    adata.obs["marker_passes_gmm"] = marker_passes
    adata.obs["confidence_passes"] = confidence_passes
    adata.obs["validation_pass"] = confidence_passes & marker_passes

    # Create validated labels column
    validated_labels = adata.obs[label_column].astype(str).copy()
    fail_mask = ~adata.obs["validation_pass"]
    n_failed = fail_mask.sum()

    if n_failed > 0:
        validated_labels[fail_mask] = unassigned_label
        pct_failed = 100 * n_failed / n_cells
        logger.info(
            f"Dual-threshold validation: {n_failed:,} cells ({pct_failed:.1f}%) "
            f"marked as '{unassigned_label}'"
        )

    adata.obs[f"{label_column}_validated"] = pd.Categorical(validated_labels)

    # Create summary DataFrame
    summary_df = pd.DataFrame(summary_rows)
    if len(summary_df) > 0:
        summary_df = summary_df.sort_values("n_cells", ascending=False).reset_index(drop=True)

    # Log overall summary
    n_pass_total = adata.obs["validation_pass"].sum()
    logger.info(
        f"Validation complete: {n_pass_total:,}/{n_cells:,} cells "
        f"({100*n_pass_total/n_cells:.1f}%) passed dual-threshold QC"
    )

    return adata, summary_df
