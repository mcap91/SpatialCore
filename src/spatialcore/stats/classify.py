"""Expression-based cell classification via thresholding.

This module provides functionality to classify cell populations based on
expression of one or more features using statistical threshold detection.

**Univariate mode** (single feature): Directly thresholds expression values.

**Multivariate mode** (multiple features): Computes a metagene score that
captures joint elevation across all features, then thresholds this score.
The metagene score is high ONLY when ALL markers are elevated, distinguishing
true multi-positive cells from single-positive cells.

Functions
---------
classify_by_threshold
    Classify cells by expression thresholding (univariate or multivariate).

References
----------
.. [1] Statistical thresholding approaches for single-cell data analysis.
"""

from pathlib import Path
from typing import List, Literal, Optional, Dict, Any, Union

import numpy as np
import anndata as ad
import matplotlib
matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

from spatialcore.core.logging import get_logger
from spatialcore.core.metadata import update_metadata
from spatialcore.stats._thresholding import (
    compute_metagene_score,
    threshold_ks,
    threshold_gmm,
)

logger = get_logger(__name__)


MetageneMethod = Literal[
    "shifted_geometric_mean",
    "geometric_mean",
    "arithmetic_mean",
    "median",
    "minimum",
]

ThresholdMethod = Literal["ks", "gmm"]


def _extract_features(
    adata: ad.AnnData,
    feature_columns: List[str],
) -> np.ndarray:
    """
    Extract feature values from AnnData.

    Searches for features in order: adata.obs, adata.var_names, adata.obsm.

    Supports colon syntax for obsm column selection:
    - "local_morans_I:MS4A1" - extract column by gene name lookup
    - "local_morans_I:0" - extract column by index

    Parameters
    ----------
    adata
        AnnData object.
    feature_columns
        List of feature names to extract. Supports:

        - Simple names: searches obs → var_names → obsm (first column)
        - Colon syntax: "obsm_key:column" for explicit obsm column selection

    Returns
    -------
    np.ndarray
        Array of shape (n_cells, n_features).

    Raises
    ------
    ValueError
        If any feature is not found.
    """
    features = []

    for col in feature_columns:
        # Check for colon syntax: "obsm_key:column_spec"
        if ":" in col:
            obsm_key, col_spec = col.split(":", 1)

            if obsm_key not in adata.obsm:
                raise ValueError(
                    f"obsm key '{obsm_key}' not found in adata.obsm. "
                    f"Available keys: {list(adata.obsm.keys())}"
                )

            arr = adata.obsm[obsm_key]

            if col_spec.isdigit():
                # Numeric index
                col_idx = int(col_spec)
                if col_idx >= arr.shape[1]:
                    raise ValueError(
                        f"Column index {col_idx} out of range for obsm['{obsm_key}'] "
                        f"with {arr.shape[1]} columns"
                    )
            else:
                # Gene name lookup from uns params
                # For "local_morans_I", look up "local_morans_params"
                base_key = obsm_key.rsplit("_", 1)[0]
                params_key = f"{base_key}_params"

                if params_key not in adata.uns:
                    raise ValueError(
                        f"Cannot look up column '{col_spec}' by name: "
                        f"'{params_key}' not found in adata.uns. "
                        f"Use numeric index instead (e.g., '{obsm_key}:0')."
                    )

                gene_list = adata.uns[params_key].get("genes", [])
                if col_spec not in gene_list:
                    raise ValueError(
                        f"Column '{col_spec}' not found in {params_key}['genes']. "
                        f"Available: {gene_list[:10]}{'...' if len(gene_list) > 10 else ''}"
                    )
                col_idx = gene_list.index(col_spec)

            features.append(arr[:, col_idx])
            logger.debug(f"Extracted column {col_idx} from obsm['{obsm_key}']")

        elif col in adata.obs.columns:
            features.append(adata.obs[col].values.astype(float))
            logger.debug(f"Found '{col}' in adata.obs")

        elif col in adata.var_names:
            gene_idx = adata.var_names.get_loc(col)
            if hasattr(adata.X, "toarray"):
                features.append(adata.X[:, gene_idx].toarray().flatten())
            else:
                features.append(np.asarray(adata.X[:, gene_idx]).flatten())
            logger.debug(f"Found '{col}' in adata.var_names (gene expression)")

        elif col in adata.obsm:
            arr = adata.obsm[col]
            if arr.ndim == 1:
                features.append(arr)
            else:
                features.append(arr[:, 0])  # Take first column
            logger.debug(f"Found '{col}' in adata.obsm")

        else:
            available_obs = list(adata.obs.columns)[:10]
            available_var = list(adata.var_names)[:10]
            raise ValueError(
                f"Feature '{col}' not found in adata.obs, adata.var_names, "
                f"or adata.obsm.\n"
                f"Available obs columns (first 10): {available_obs}\n"
                f"Available genes (first 10): {available_var}"
            )

    return np.column_stack(features)


def _plot_gpairs(
    scores: np.ndarray,
    feature_values: np.ndarray,
    feature_columns: List[str],
    cluster_labels: np.ndarray,
    threshold: float,
    metagene_method: str,
    threshold_method: str,
    probability_cutoff: float,
    output_path: Path,
    n_sample: int = 20000,
    seed: int = 42,
) -> None:
    """
    Generate gpairs-style histogram + scatter matrix plot.

    Creates a combined figure with:
    - Top row: Histogram of metagene scores with density curves
    - Bottom rows: Triangular matrix of gene-gene scatter plots

    Parameters
    ----------
    scores
        Metagene scores for valid cells.
    feature_values
        Feature values for valid cells (n_cells, n_features).
    feature_columns
        List of feature names.
    cluster_labels
        Cluster assignments (0=low, 1=high).
    threshold
        Threshold value used for classification.
    metagene_method
        Method used for metagene calculation.
    threshold_method
        Method used for threshold detection.
    probability_cutoff
        Probability cutoff used for GMM cluster assignment.
    output_path
        Path to save the figure.
    n_sample
        Maximum number of points to plot (for performance).
    seed
        Random seed for reproducible subsampling.
    """
    n_genes = len(feature_columns)
    n_pairs = n_genes * (n_genes - 1) // 2

    # Subsample if needed
    n_cells = len(scores)
    if n_cells > n_sample:
        rng = np.random.default_rng(seed)
        idx = rng.choice(n_cells, size=n_sample, replace=False)
        scores_plot = scores[idx]
        features_plot = feature_values[idx]
        clusters_plot = cluster_labels[idx]
    else:
        scores_plot = scores
        features_plot = feature_values
        clusters_plot = cluster_labels
        n_sample = n_cells

    # Calculate statistics
    n_low = int((cluster_labels == 0).sum())
    n_high = int((cluster_labels == 1).sum())
    n_total = n_low + n_high
    pct_low = 100 * n_low / n_total
    pct_high = 100 * n_high / n_total

    # Colors: Orange for Low, Blue for High (matching reference image)
    color_low = "#ff7f0e"   # Orange
    color_high = "#1f77b4"  # Blue

    # Determine layout
    # 1 gene: 0 pairs -> histogram only
    # 2 genes: 1 pair -> 1 row of scatters
    # 3 genes: 3 pairs -> 1 row of 3 scatters
    # 4 genes: 6 pairs -> 2 rows of 3 scatters
    # 5 genes: 10 pairs -> 2 rows (3+3+4) or 4 rows
    # 6 genes: 15 pairs -> 3 rows of 5 or similar
    if n_pairs == 0:
        # Single gene: histogram only
        scatter_rows = 0
        scatter_cols = 1  # Need at least 1 column for histogram
    elif n_pairs <= 3:
        scatter_rows = 1
        scatter_cols = n_pairs
    elif n_pairs <= 6:
        scatter_rows = 2
        scatter_cols = 3
    else:
        scatter_cols = min(n_pairs, 5)
        scatter_rows = (n_pairs + scatter_cols - 1) // scatter_cols

    # Figure size
    fig_width = max(10, scatter_cols * 4)
    fig_height = 4 + scatter_rows * 3.5

    fig = plt.figure(figsize=(fig_width, fig_height))

    # Create GridSpec: 1 row for histogram, scatter_rows for scatters
    if scatter_rows == 0:
        gs = gridspec.GridSpec(1, 1)  # Just histogram
    else:
        gs = gridspec.GridSpec(
            1 + scatter_rows, scatter_cols,
            height_ratios=[1.2] + [1] * scatter_rows,
            hspace=0.35, wspace=0.3
        )

    # === Top row: Histogram spanning all columns ===
    ax_hist = fig.add_subplot(gs[0, :])

    # Split data by cluster for histogram
    scores_low = scores_plot[clusters_plot == 0]
    scores_high = scores_plot[clusters_plot == 1]

    # Compute common bins - guard against degenerate case
    all_scores = np.concatenate([scores_low, scores_high])
    score_range = all_scores.max() - all_scores.min()
    if score_range < 1e-10:
        # Degenerate case: all scores identical
        logger.warning(
            f"All metagene scores are identical ({all_scores[0]:.4f}). "
            "Classification may not be meaningful."
        )
        bins = np.array([all_scores.min() - 0.5, all_scores.max() + 0.5])
    else:
        bins = np.linspace(all_scores.min(), all_scores.max(), 50)

    # Plot histograms with density
    ax_hist.hist(
        scores_low, bins=bins, alpha=0.6, color=color_low,
        label=f"Low: n={n_low:,} ({pct_low:.1f}%)", density=True
    )
    ax_hist.hist(
        scores_high, bins=bins, alpha=0.6, color=color_high,
        label=f"High: n={n_high:,} ({pct_high:.1f}%)", density=True
    )

    # Add density curves (may fail for degenerate data, e.g., sparse markers)
    x_kde = np.linspace(all_scores.min(), all_scores.max(), 200)

    if len(scores_low) > 10:
        try:
            kde_low = stats.gaussian_kde(scores_low)
            ax_hist.plot(x_kde, kde_low(x_kde), color=color_low, lw=2)
        except np.linalg.LinAlgError:
            logger.debug("KDE skipped for low cluster (singular covariance)")

    if len(scores_high) > 10:
        try:
            kde_high = stats.gaussian_kde(scores_high)
            ax_hist.plot(x_kde, kde_high(x_kde), color=color_high, lw=2)
        except np.linalg.LinAlgError:
            logger.debug("KDE skipped for high cluster (singular covariance)")

    # Threshold line
    if threshold_method == "gmm":
        threshold_label = f"Threshold: {threshold:.4f} (P cutoff: {probability_cutoff})"
    else:
        threshold_label = f"Threshold: {threshold:.4f}"
    ax_hist.axvline(threshold, color="#2ca02c", lw=2.5, ls="--",
                    label=threshold_label)

    # Labels and title
    gene_str = " + ".join(feature_columns)
    ax_hist.set_title(
        f"Metagene Score Distribution: {gene_str}",
        fontsize=14, fontweight="bold"
    )
    ax_hist.set_xlabel(f"Metagene Score ({gene_str})", fontsize=11)
    ax_hist.set_ylabel("Density", fontsize=11)

    # Subtitle with method info (placed below title via figure suptitle approach)
    if threshold_method == "gmm":
        subtitle = f"Method: {metagene_method} + GMM | P(high) cutoff: {probability_cutoff}"
    else:
        subtitle = f"Method: {metagene_method} + KS"
    fig.text(0.5, 0.98, subtitle, ha="center", fontsize=9, color="gray")

    ax_hist.legend(loc="upper right", fontsize=9)

    # === Bottom rows: Gene-gene scatter plots ===
    # Generate all pairs
    pairs = []
    for i in range(n_genes):
        for j in range(i + 1, n_genes):
            pairs.append((i, j))

    for pair_idx, (i, j) in enumerate(pairs):
        row = 1 + pair_idx // scatter_cols
        col = pair_idx % scatter_cols

        ax = fig.add_subplot(gs[row, col])

        # Extract gene values
        x_vals = features_plot[:, i]
        y_vals = features_plot[:, j]

        # Plot Low cluster (orange) first, then High (blue) on top
        mask_low = clusters_plot == 0
        mask_high = clusters_plot == 1

        ax.scatter(
            x_vals[mask_low], y_vals[mask_low],
            c=color_low, s=3, alpha=0.4, rasterized=True, label="Low"
        )
        ax.scatter(
            x_vals[mask_high], y_vals[mask_high],
            c=color_high, s=3, alpha=0.6, rasterized=True, label="High"
        )

        # Diagonal reference line
        lims = [
            min(ax.get_xlim()[0], ax.get_ylim()[0]),
            max(ax.get_xlim()[1], ax.get_ylim()[1]),
        ]
        ax.plot(lims, lims, 'k--', alpha=0.3, lw=1)

        # Labels
        gene_x = feature_columns[i]
        gene_y = feature_columns[j]
        ax.set_xlabel(gene_x, fontsize=10)
        ax.set_ylabel(gene_y, fontsize=10)
        ax.set_title(f"{gene_x} vs {gene_y} (n={n_sample:,} shown)", fontsize=10)

        # Only add legend to first scatter
        if pair_idx == 0:
            ax.legend(loc="upper left", fontsize=8, markerscale=2)

    # Hide unused subplots
    total_subplots = scatter_rows * scatter_cols
    for idx in range(n_pairs, total_subplots):
        row = 1 + idx // scatter_cols
        col = idx % scatter_cols
        ax = fig.add_subplot(gs[row, col])
        ax.axis("off")

    plt.tight_layout()

    # Save figure
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    logger.info(f"Saved gpairs plot: {output_path}")


def classify_by_threshold(
    adata: ad.AnnData,
    feature_columns: List[str],
    metagene_method: MetageneMethod = "shifted_geometric_mean",
    threshold_method: ThresholdMethod = "gmm",
    pseudocount: float = 0.1,
    background_quantile: float = 0.5,
    probability_cutoff: float = 0.3,
    n_components: int = 2,
    max_cells: int = 20000,
    column_prefix: str = "threshold",
    seed: int = 42,
    plot: bool = True,
    output_dir: Optional[Union[str, Path]] = None,
    n_sample_plot: int = 20000,
    copy: bool = False,
) -> ad.AnnData:
    """
    Classify cells by expression thresholding.

    Works with single features (univariate) or multiple features (multivariate).
    For a single feature, directly thresholds expression values. For multiple
    features, computes a metagene score that captures joint elevation, then
    finds a threshold separating high vs low populations.

    Designed for sparse single-cell data with skewed distributions. The metagene
    score is high ONLY when ALL markers are elevated, unlike Euclidean clustering
    which treats single-positive cells the same as multi-positive cells.

    Parameters
    ----------
    adata
        AnnData object with features in adata.obs, adata.var_names, or adata.obsm.
    feature_columns
        List of feature names to threshold. Supports:

        - Column names in adata.obs (e.g., ["CD38_score", "MS4A1_score"])
        - Gene names in adata.var_names (e.g., ["SCGB1A1", "SCGB3A1"])
        - Keys in adata.obsm (takes first column by default)
        - Colon syntax for obsm column selection:
          - "local_morans_I:MS4A1" - extract by gene name (from uns params)
          - "local_morans_I:0" - extract by column index

        For univariate thresholding, provide a single feature: ["MS4A1"].
        For multivariate, provide multiple: ["SCGB1A1", "SCGB3A1", "SCGB3A2"].
    metagene_method
        Method for combining features into a metagene score (multivariate only):

        - "shifted_geometric_mean" (default): exp(mean(log(x + c))) - c.
          Tolerates dropout while requiring joint expression.
        - "geometric_mean": (x1 * x2 * ... * xn)^(1/n). Very strict, any
          zero collapses score to zero.
        - "arithmetic_mean": mean of features. Permissive, allows single-positive.
        - "median": middle value. Requires majority of markers elevated.
        - "minimum": strictest, score equals weakest marker.

        For univariate mode (single feature), this parameter is ignored.
    threshold_method
        Method for threshold detection:

        - "gmm" (default): Gaussian Mixture Model. Fits two Gaussians and finds
          intersection. Returns posterior probabilities. Works well when both
          populations are substantial (>10% each). May struggle with sparse
          markers where most cells have zero expression.
        - "ks": KS-inspired method. Finds where empirical distribution
          deviates maximally from estimated background. **Recommended for
          sparse markers** (<30% positive, e.g., MS4A1, PDCD1). Handles
          zero-inflated distributions well.
    pseudocount
        Pseudocount for shifted_geometric_mean method. Default: 0.1.
        Larger values increase dropout tolerance.
    background_quantile
        KS method only. Fraction of data to estimate background distribution.
        Default: 0.5 (lower 50%).
    probability_cutoff
        GMM method only. P(high) threshold for cluster assignment.
        Default: 0.3.
    n_components
        GMM method only. Number of Gaussian components. Default: 2.

        - 2: Standard bimodal separation (low vs high populations)
        - 3: Trimodal for spatial data where marker expression often has
          three peaks: dropout (~0), moderate (~1-2), high (~4-5).
          The threshold is set between background and first signal peak.
          Use this for marker validation on spatial transcriptomics data.
    max_cells
        Maximum cells for threshold fitting. Default: 20000.
        Reduces computation time for large datasets while ensuring
        reliable threshold estimation.
    column_prefix
        Prefix for result columns and keys. Results stored in:

        - adata.obs[f"{column_prefix}_score"]: metagene/expression scores
        - adata.obs[f"{column_prefix}_probability"]: P(high|score) for GMM,
          deviation score for KS
        - adata.obs[f"{column_prefix}_cluster"]: 0 (low), 1 (high), -1 (invalid)
        - adata.uns[f"{column_prefix}_params"]: dict containing threshold,
          n_high, n_low, n_invalid, feature_columns, and method-specific
          params. GMM: gmm_means, gmm_stds, gmm_weights, probability_cutoff.
          KS: background_mean, background_std, background_quantile.
    seed
        Random seed for reproducibility. Default: 42.
    plot
        If True (default), generate a gpairs-style diagnostic plot with
        histogram and gene-gene scatter matrix. Requires output_dir.
    output_dir
        Directory to save the gpairs plot. Required when plot=True.
        The plot will be saved as "{output_dir}/{column_prefix}_gpairs.png".
    n_sample_plot
        Maximum number of cells to include in the plot. Default: 20000.
        Reduces file size and rendering time for large datasets.
    copy
        If True, operate on a copy of adata. Default: False.

    Returns
    -------
    AnnData
        AnnData with classification results stored in obs and uns.
        If plot=True, also saves a gpairs plot to output_dir.

    Raises
    ------
    ValueError
        If feature_columns not found, invalid method, or insufficient data.

    Notes
    -----
    **Univariate vs Multivariate Mode:**

    - Univariate (1 feature): Directly thresholds the single feature.
    - Multivariate (2+ features): Computes metagene score first, then thresholds.

    **Metagene methods comparison** (multivariate mode):

    +------------------------+------------+------------------+------------------+
    | Method                 | Strictness | Dropout Tolerance| Joint Expression |
    +========================+============+==================+==================+
    | geometric_mean         | Very strict| None             | Yes              |
    | shifted_geometric_mean | Moderate   | Good             | Yes              |
    | arithmetic_mean        | Permissive | High             | No               |
    | median                 | Moderate   | Good             | Partial          |
    | minimum                | Very strict| None             | Yes              |
    +------------------------+------------+------------------+------------------+

    **Threshold method comparison:**

    - KS: Non-parametric, best for sparse (<5% to 30%), skewed distributions.
    - GMM: Provides posterior probabilities, best when both populations substantial.

    Examples
    --------
    Univariate thresholding (single gene):

    >>> import scanpy as sc
    >>> from spatialcore.stats import classify_by_threshold
    >>> adata = sc.read_h5ad("annotated.h5ad")
    >>> # Threshold a single marker
    >>> adata = classify_by_threshold(
    ...     adata,
    ...     feature_columns=["MS4A1"],
    ...     column_prefix="b_cell_marker",
    ...     plot=False,
    ... )

    Multivariate thresholding (joint expression of multiple markers):

    >>> # Identify cells with joint expression of epithelial markers
    >>> adata = classify_by_threshold(
    ...     adata,
    ...     feature_columns=["SCGB1A1", "SCGB3A1", "SCGB3A2"],
    ...     column_prefix="epithelial",
    ...     output_dir="./plots",
    ... )
    >>> # Check results
    >>> print(adata.obs["epithelial_cluster"].value_counts())
    0    9000
    1    1000
    >>> # Get high-expressing cells
    >>> epithelial_cells = adata[adata.obs["epithelial_cluster"] == 1]

    Using KS method for sparse markers:

    >>> adata = classify_by_threshold(
    ...     adata,
    ...     feature_columns=["PDCD1"],  # Sparse marker
    ...     threshold_method="ks",
    ...     column_prefix="pd1_high",
    ...     plot=False,
    ... )

    See Also
    --------
    spatialcore.spatial.identify_niches : Cluster cells by neighborhood composition.
    """
    # -------------------------------------------------------------------------
    # Input validation
    # -------------------------------------------------------------------------
    if not isinstance(feature_columns, list) or len(feature_columns) == 0:
        raise ValueError(
            "feature_columns must be a non-empty list of feature names."
        )

    valid_metagene_methods = [
        "shifted_geometric_mean",
        "geometric_mean",
        "arithmetic_mean",
        "median",
        "minimum",
    ]
    if metagene_method not in valid_metagene_methods:
        raise ValueError(
            f"Invalid metagene_method '{metagene_method}'. "
            f"Must be one of: {valid_metagene_methods}"
        )

    valid_threshold_methods = ["ks", "gmm"]
    if threshold_method not in valid_threshold_methods:
        raise ValueError(
            f"Invalid threshold_method '{threshold_method}'. "
            f"Must be one of: {valid_threshold_methods}"
        )

    if pseudocount <= 0:
        raise ValueError(f"pseudocount must be > 0, got {pseudocount}")

    if not 0 < background_quantile < 1:
        raise ValueError(
            f"background_quantile must be in (0, 1), got {background_quantile}"
        )

    if not 0 < probability_cutoff < 1:
        raise ValueError(
            f"probability_cutoff must be in (0, 1), got {probability_cutoff}"
        )

    if plot and output_dir is None:
        raise ValueError(
            "output_dir is required when plot=True. "
            "Provide a directory path or set plot=False."
        )

    # -------------------------------------------------------------------------
    # Handle copy
    # -------------------------------------------------------------------------
    adata = adata.copy() if copy else adata

    # -------------------------------------------------------------------------
    # Extract features
    # -------------------------------------------------------------------------
    logger.info(
        f"Classifying by threshold: {len(feature_columns)} feature(s), "
        f"metagene={metagene_method}, threshold={threshold_method}"
    )
    logger.info(f"Features: {feature_columns}")

    feature_values = _extract_features(adata, feature_columns)
    n_cells = adata.n_obs

    # -------------------------------------------------------------------------
    # Handle invalid values (NaN, Inf)
    # -------------------------------------------------------------------------
    valid_mask = np.all(np.isfinite(feature_values), axis=1)
    n_invalid = int((~valid_mask).sum())

    if n_invalid > 0:
        logger.warning(
            f"{n_invalid} cells have NaN/Inf values and will be marked "
            f"as cluster=-1"
        )

    if valid_mask.sum() < 100:
        raise ValueError(
            f"Only {valid_mask.sum()} valid cells (non-NaN/Inf). "
            "Need at least 100 cells for threshold detection."
        )

    feature_values_clean = feature_values[valid_mask]
    n_valid = int(valid_mask.sum())
    logger.info(f"Valid cells: {n_valid:,} / {n_cells:,}")

    # -------------------------------------------------------------------------
    # Check for negative values (incompatible with geometric mean methods)
    # -------------------------------------------------------------------------
    has_negative = np.any(feature_values_clean < 0)
    if has_negative and metagene_method in ("shifted_geometric_mean", "geometric_mean"):
        raise ValueError(
            f"Feature values contain negative numbers, which are incompatible with "
            f"metagene_method='{metagene_method}' (log of negative values is undefined). "
            f"Use metagene_method='arithmetic_mean' or 'median' instead.\n\n"
            f"Common cases with negative values:\n"
            f"  - Local Moran's I (negative = spatial outlier)\n"
            f"  - Z-scores or scaled expression data\n"
            f"  - Differential expression log-fold changes"
        )

    # -------------------------------------------------------------------------
    # Compute metagene scores
    # -------------------------------------------------------------------------
    metagene_scores = compute_metagene_score(
        feature_values_clean,
        method=metagene_method,
        pseudocount=pseudocount,
    )

    logger.info(
        f"Metagene scores: min={metagene_scores.min():.4f}, "
        f"max={metagene_scores.max():.4f}, mean={metagene_scores.mean():.4f}"
    )

    # -------------------------------------------------------------------------
    # Check for sparse marker expression (warn if GMM may be inappropriate)
    # -------------------------------------------------------------------------
    # Count cells with all markers at zero (metagene score ~ 0)
    zero_threshold = 1e-6
    n_all_zero = int((metagene_scores < zero_threshold).sum())
    pct_all_zero = 100 * n_all_zero / n_valid

    if pct_all_zero >= 50.0 and threshold_method == "gmm":
        import warnings
        warnings.warn(
            f"{pct_all_zero:.1f}% of cells have zero expression for all markers. "
            f"GMM will likely separate zeros from non-zeros rather than finding "
            f"a meaningful biological threshold. Consider using threshold_method='ks' "
            f"which is designed for sparse marker detection.",
            UserWarning,
            stacklevel=2,
        )

    # -------------------------------------------------------------------------
    # Downsample for threshold fitting (GMM only)
    # -------------------------------------------------------------------------
    np.random.seed(seed)
    sample_size = min(n_valid, max_cells)  # Use all cells up to max_cells

    if sample_size < n_valid and threshold_method == "gmm":
        sample_indices = np.random.choice(
            n_valid, size=sample_size, replace=False
        )
        sampled_scores = metagene_scores[sample_indices]
        logger.debug(f"Downsampled to {sample_size:,} cells for GMM fitting")
    else:
        sampled_scores = metagene_scores

    # -------------------------------------------------------------------------
    # Find threshold
    # -------------------------------------------------------------------------
    if threshold_method == "ks":
        threshold, probability_scores, method_params = threshold_ks(
            metagene_scores,
            background_quantile=background_quantile,
        )
        cluster_labels = (metagene_scores >= threshold).astype(int)

    else:  # gmm
        threshold, cluster_labels, probability_scores, method_params = threshold_gmm(
            sampled_scores if sample_size < n_valid else metagene_scores,
            probability_cutoff=probability_cutoff,
            n_components=n_components,
            random_state=seed,
        )
        # Recompute probabilities for all cells if downsampled
        if sample_size < n_valid:
            # REUSE the fitted GMM model instead of fitting a new one
            # This ensures threshold and probabilities are consistent
            gmm = method_params["gmm_model"]
            high_component = method_params["high_component_idx"]
            posteriors_all = gmm.predict_proba(metagene_scores.reshape(-1, 1))

            # For trimodal (n_components >= 3), sum probabilities of non-background
            if n_components >= 3:
                sorted_indices = method_params["sorted_component_indices"]
                signal_components = sorted_indices[1:]  # All except background
                probability_scores = posteriors_all[:, signal_components].sum(axis=1)
            else:
                probability_scores = posteriors_all[:, high_component]

            cluster_labels = (probability_scores > probability_cutoff).astype(int)

    logger.info(f"Threshold: {threshold:.4f}")

    # -------------------------------------------------------------------------
    # Store results in adata.obs
    # -------------------------------------------------------------------------
    score_col = f"{column_prefix}_score"
    prob_col = f"{column_prefix}_probability"
    cluster_col = f"{column_prefix}_cluster"

    # Initialize with NaN/-1 for invalid cells
    adata.obs[score_col] = np.nan
    adata.obs[prob_col] = np.nan
    adata.obs[cluster_col] = -1

    # Fill in valid cells
    valid_indices = adata.obs.index[valid_mask]
    adata.obs.loc[valid_indices, score_col] = metagene_scores
    adata.obs.loc[valid_indices, prob_col] = probability_scores
    adata.obs.loc[valid_indices, cluster_col] = cluster_labels

    # Convert cluster to int (was float due to NaN initialization)
    adata.obs[cluster_col] = adata.obs[cluster_col].astype(int)

    # -------------------------------------------------------------------------
    # Calculate statistics
    # -------------------------------------------------------------------------
    n_high = int((cluster_labels == 1).sum())
    n_low = int((cluster_labels == 0).sum())
    pct_high = 100 * n_high / len(cluster_labels)
    pct_low = 100 * n_low / len(cluster_labels)

    logger.info(f"Cluster 0 (low): {n_low:,} cells ({pct_low:.1f}%)")
    logger.info(f"Cluster 1 (high): {n_high:,} cells ({pct_high:.1f}%)")

    # -------------------------------------------------------------------------
    # Store parameters in adata.uns
    # -------------------------------------------------------------------------
    uns_key = f"{column_prefix}_params"
    adata.uns[uns_key] = {
        "feature_columns": feature_columns,
        "metagene_method": metagene_method,
        "threshold_method": threshold_method,
        "threshold": threshold,
        "pseudocount": pseudocount,
        "n_high": n_high,
        "n_low": n_low,
        "n_invalid": n_invalid,
        "n_total": n_cells,
        "seed": seed,
        **method_params,
    }

    # -------------------------------------------------------------------------
    # Update metadata
    # -------------------------------------------------------------------------
    update_metadata(
        adata,
        function_name="classify_by_threshold",
        parameters={
            "feature_columns": feature_columns,
            "metagene_method": metagene_method,
            "threshold_method": threshold_method,
            "pseudocount": pseudocount,
            "background_quantile": background_quantile,
            "probability_cutoff": probability_cutoff,
            "column_prefix": column_prefix,
        },
        outputs={
            "obs_score": score_col,
            "obs_probability": prob_col,
            "obs_cluster": cluster_col,
            "uns_params": uns_key,
            "threshold": threshold,
            "n_high": n_high,
            "n_low": n_low,
        },
    )

    # -------------------------------------------------------------------------
    # Generate gpairs plot
    # -------------------------------------------------------------------------
    if plot:
        output_path = Path(output_dir) / f"{column_prefix}_gpairs.png"
        _plot_gpairs(
            scores=metagene_scores,
            feature_values=feature_values_clean,
            feature_columns=feature_columns,
            cluster_labels=cluster_labels,
            threshold=threshold,
            metagene_method=metagene_method,
            threshold_method=threshold_method,
            probability_cutoff=probability_cutoff,
            output_path=output_path,
            n_sample=n_sample_plot,
            seed=seed,
        )

    return adata
