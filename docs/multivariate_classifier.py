#!/usr/bin/env python3
"""
Multivariate Metagene Classifier - Standalone Implementation

Identifies cell populations based on joint expression of multiple markers.
No dependencies on Temporal, process_activity, or workflow orchestration.

Example usage:
    python multivariate_classifier.py --input data.h5ad --markers GeneA,GeneB,GeneC

Author: Spatial Biology Team
Date: December 2025
"""

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.stats import norm
from sklearn.mixture import GaussianMixture
from typing import Literal, Tuple, List, Optional, Dict, Any
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


# =============================================================================
# Type Definitions
# =============================================================================

MetageneMethod = Literal[
    "geometric_mean",
    "shifted_geometric_mean",
    "arithmetic_mean",
    "median",
    "minimum"
]

ThresholdMethod = Literal["ks", "gmm"]


# =============================================================================
# Core Functions
# =============================================================================

def extract_features(
    adata: sc.AnnData,
    feature_columns: List[str],
) -> Tuple[np.ndarray, List[str]]:
    """
    Extract feature values for multivariate analysis.

    Searches for features in:
    1. adata.obs columns (metadata)
    2. adata.var_names (gene expression)
    3. adata.obsm columns (embeddings)

    Args:
        adata: AnnData object
        feature_columns: List of feature names to extract

    Returns:
        Tuple of (feature_values array [n_cells, n_features], valid_column_names)
    """
    features = []
    valid_columns = []

    for col in feature_columns:
        if col in adata.obs.columns:
            # Found in obs (metadata)
            features.append(adata.obs[col].values.astype(float))
            valid_columns.append(col)
        elif col in adata.var_names:
            # Found in var_names (gene expression)
            gene_idx = adata.var_names.get_loc(col)
            if hasattr(adata.X, 'toarray'):
                # Sparse matrix
                features.append(adata.X[:, gene_idx].toarray().flatten())
            else:
                features.append(adata.X[:, gene_idx].flatten())
            valid_columns.append(col)
        elif col in adata.obsm:
            # Found in obsm (embeddings)
            features.append(adata.obsm[col].flatten())
            valid_columns.append(col)
        else:
            raise ValueError(
                f"Feature column '{col}' not found in adata.obs, adata.var_names, or adata.obsm. "
                f"Available obs columns: {list(adata.obs.columns)[:10]}..."
            )

    return np.column_stack(features), valid_columns


def compute_metagene_score(
    feature_values: np.ndarray,
    method: MetageneMethod = "shifted_geometric_mean",
    pseudocount: float = 0.1,
) -> np.ndarray:
    """
    Compute a metagene score that captures joint elevation across multiple features.

    Args:
        feature_values: (n_samples, n_features) array of expression values
        method: Aggregation method:
            - "shifted_geometric_mean" (default): exp(mean(log(x + c))) - c
            - "geometric_mean": (x1 * x2 * ... * xn)^(1/n)
            - "arithmetic_mean": (x1 + x2 + ... + xn) / n
            - "median": median(x1, x2, ..., xn)
            - "minimum": min(x1, x2, ..., xn)
        pseudocount: Pseudocount for shifted_geometric_mean (default: 0.1)

    Returns:
        (n_samples,) array of metagene scores
    """
    if method == "geometric_mean":
        eps = 1e-10
        return np.exp(np.mean(np.log(feature_values + eps), axis=1))

    elif method == "shifted_geometric_mean":
        # Add pseudocount to tolerate dropout, then subtract to restore scale
        shifted = feature_values + pseudocount
        return np.exp(np.mean(np.log(shifted), axis=1)) - pseudocount

    elif method == "arithmetic_mean":
        return np.mean(feature_values, axis=1)

    elif method == "median":
        return np.median(feature_values, axis=1)

    elif method == "minimum":
        return np.min(feature_values, axis=1)

    else:
        raise ValueError(f"Unknown metagene method: {method}")


def threshold_ks(
    metagene_scores: np.ndarray,
    background_quantile: float = 0.5,
) -> Tuple[float, np.ndarray, float, float]:
    """
    Find threshold where empirical distribution deviates maximally from background.

    Uses KS-inspired approach: estimates background from lower quantile,
    then finds where data "breaks away" from expected background CDF.

    Args:
        metagene_scores: 1D array of metagene scores
        background_quantile: Fraction of data to estimate background (default: 0.5)

    Returns:
        Tuple of (threshold, deviation_scores, bg_mean, bg_std)
    """
    sorted_scores = np.sort(metagene_scores)
    n = len(sorted_scores)

    # Estimate background from lower quantile
    bg_cutoff_idx = int(n * background_quantile)
    background_scores = sorted_scores[:bg_cutoff_idx]
    bg_mean = float(np.mean(background_scores))
    bg_std = float(np.std(background_scores))

    # Handle edge case: if background has no variance
    if bg_std < 1e-10:
        q25, q75 = np.percentile(sorted_scores, [25, 75])
        iqr = q75 - q25
        if iqr > 1e-10:
            bg_std = float(iqr / 1.35)
        else:
            data_range = sorted_scores[-1] - sorted_scores[0]
            bg_std = float(max(data_range * 0.1, 1e-6))
        logger.info(f"Background std was ~0, using fallback estimate: {bg_std:.6f}")

    # Compute empirical CDF
    empirical_cdf = np.arange(1, n + 1) / n

    # Compute expected CDF if all data were background
    expected_cdf = norm.cdf(sorted_scores, loc=bg_mean, scale=bg_std)

    # D = excess over expected (positive where we have "extra" high values)
    D = empirical_cdf - expected_cdf

    # Threshold at maximum deviation
    max_idx = np.argmax(D)
    threshold = float(sorted_scores[max_idx])

    # Sanity check: threshold should be above background mean
    if threshold <= bg_mean:
        threshold = float(np.percentile(sorted_scores, 90))
        logger.info(f"KS threshold below background mean, using 90th percentile: {threshold:.4f}")

    # Compute pseudo-probability for all cells (deviation score)
    max_score = sorted_scores[-1]
    score_range = max_score - threshold
    if score_range < 1e-10:
        score_range = 1e-10

    deviation_scores = np.clip(
        (metagene_scores - threshold) / score_range,
        0, 1
    )

    return threshold, deviation_scores, bg_mean, bg_std


def threshold_gmm(
    metagene_scores: np.ndarray,
    n_components: int = 2,
    probability_cutoff: float = 0.3,
    seed: int = 42,
) -> Tuple[float, np.ndarray, Dict[str, Any]]:
    """
    Find threshold using Gaussian Mixture Model.

    Args:
        metagene_scores: 1D array of metagene scores
        n_components: Number of GMM components (default: 2)
        probability_cutoff: P(high) threshold for classification (default: 0.3)
        seed: Random seed

    Returns:
        Tuple of (threshold, probability_high, gmm_params_dict)
    """
    gmm = GaussianMixture(
        n_components=n_components,
        random_state=seed,
        n_init=10,
    )
    gmm.fit(metagene_scores.reshape(-1, 1))

    # Identify high component (higher mean)
    component_means = gmm.means_.flatten()
    component_stds = np.sqrt(gmm.covariances_.flatten())
    high_component = int(np.argmax(component_means))
    low_component = 1 - high_component

    low_mean, high_mean = component_means[low_component], component_means[high_component]

    # Find threshold at P(high) = 0.5
    x_range = np.linspace(low_mean, high_mean, 1000)
    posteriors = gmm.predict_proba(x_range.reshape(-1, 1))
    diff = posteriors[:, high_component] - 0.5
    cross_idx = np.where(np.diff(np.sign(diff)))[0]

    if len(cross_idx) > 0:
        threshold = float(x_range[cross_idx[0]])
    else:
        threshold = float((low_mean + high_mean) / 2)

    # Compute posterior P(high|score) for all cells
    posteriors_all = gmm.predict_proba(metagene_scores.reshape(-1, 1))
    probability_high = posteriors_all[:, high_component]

    gmm_params = {
        "gmm_means": component_means.tolist(),
        "gmm_stds": component_stds.tolist(),
        "gmm_weights": gmm.weights_.tolist(),
        "high_component_idx": high_component,
        "probability_cutoff": probability_cutoff,
    }

    return threshold, probability_high, gmm_params


# =============================================================================
# Main Classifier Class
# =============================================================================

class MultivariateMetageneClassifier:
    """
    Multivariate classifier for identifying cell populations based on joint marker expression.

    Example:
        >>> classifier = MultivariateMetageneClassifier()
        >>> results = classifier.fit_transform(
        ...     adata,
        ...     markers=["GeneA", "GeneB", "GeneC"],
        ...     metagene_method="shifted_geometric_mean",
        ...     threshold_method="ks"
        ... )
        >>> print(f"Found {results['n_high']} high-expression cells")
    """

    def __init__(self):
        self.threshold_ = None
        self.metagene_method_ = None
        self.threshold_method_ = None
        self.feature_columns_ = None
        self.params_ = {}

    def fit_transform(
        self,
        adata: sc.AnnData,
        markers: List[str],
        metagene_method: MetageneMethod = "shifted_geometric_mean",
        threshold_method: ThresholdMethod = "ks",
        pseudocount: float = 0.1,
        background_quantile: float = 0.5,  # KS-specific
        probability_cutoff: float = 0.3,   # GMM-specific
        column_prefix: str = "multivariate",
        seed: int = 42,
        downsample_prop: float = 0.1,
        max_cells: int = 50000,
    ) -> Dict[str, Any]:
        """
        Fit classifier and transform data.

        Args:
            adata: AnnData object with expression data
            markers: List of marker gene names
            metagene_method: Method for combining markers
            threshold_method: Method for finding threshold ("ks" or "gmm")
            pseudocount: Pseudocount for shifted_geometric_mean
            background_quantile: KS-specific: fraction for background estimation
            probability_cutoff: GMM-specific: P(high) threshold
            column_prefix: Prefix for new columns added to adata.obs
            seed: Random seed
            downsample_prop: Proportion of cells for threshold fitting
            max_cells: Maximum cells for threshold fitting

        Returns:
            Dictionary with results and metadata
        """
        self.metagene_method_ = metagene_method
        self.threshold_method_ = threshold_method
        self.feature_columns_ = markers

        logger.info(f"Fitting multivariate classifier")
        logger.info(f"  Markers: {markers}")
        logger.info(f"  Metagene method: {metagene_method}")
        logger.info(f"  Threshold method: {threshold_method}")

        # Extract features
        feature_values, valid_columns = extract_features(adata, markers)
        self.feature_columns_ = valid_columns

        # Track valid cells (no NaN/Inf)
        valid_mask = np.all(np.isfinite(feature_values), axis=1)
        feature_values_clean = feature_values[valid_mask]
        logger.info(f"  Valid cells: {np.sum(valid_mask)} / {len(valid_mask)}")

        # Compute metagene score for ALL valid cells
        metagene_scores = compute_metagene_score(
            feature_values_clean,
            method=metagene_method,
            pseudocount=pseudocount
        )
        logger.info(f"  Metagene scores: min={metagene_scores.min():.3f}, "
                   f"max={metagene_scores.max():.3f}, mean={metagene_scores.mean():.3f}")

        # Downsample for threshold fitting
        np.random.seed(seed)
        sample_size = min(int(len(metagene_scores) * downsample_prop), max_cells)
        sample_indices = np.random.choice(len(metagene_scores), size=sample_size, replace=False)
        sampled_scores = metagene_scores[sample_indices]
        logger.info(f"  Downsampled to {sample_size} cells for fitting")

        # Find threshold
        if threshold_method == "ks":
            threshold, probability_high, bg_mean, bg_std = threshold_ks(
                metagene_scores,
                background_quantile=background_quantile
            )
            cluster_labels = (metagene_scores >= threshold).astype(int)

            self.params_ = {
                "threshold_method": "ks",
                "background_quantile": background_quantile,
                "background_mean": bg_mean,
                "background_std": bg_std,
                "pseudocount": pseudocount,
            }

        elif threshold_method == "gmm":
            threshold, probability_high, gmm_params = threshold_gmm(
                sampled_scores,  # Fit on sample
                probability_cutoff=probability_cutoff,
                seed=seed
            )
            # Re-compute probabilities for ALL cells using fitted GMM
            _, probability_high, _ = threshold_gmm(
                metagene_scores,
                probability_cutoff=probability_cutoff,
                seed=seed
            )
            cluster_labels = (probability_high > probability_cutoff).astype(int)

            self.params_ = {
                "threshold_method": "gmm",
                "pseudocount": pseudocount,
                **gmm_params
            }

        self.threshold_ = threshold

        # Store results in adata.obs
        score_col = f"{column_prefix}_metagene_score"
        prob_col = f"{column_prefix}_probability"
        cluster_col = f"{column_prefix}_cluster"

        adata.obs[score_col] = np.nan
        adata.obs[prob_col] = np.nan
        adata.obs[cluster_col] = -1

        valid_indices = adata.obs.index[valid_mask]
        adata.obs.loc[valid_indices, score_col] = metagene_scores
        adata.obs.loc[valid_indices, prob_col] = probability_high
        adata.obs.loc[valid_indices, cluster_col] = cluster_labels

        n_high = int(np.sum(cluster_labels == 1))
        n_low = int(np.sum(cluster_labels == 0))

        logger.info(f"  Threshold: {threshold:.4f}")
        logger.info(f"  Cluster 0 (low): {n_low} cells ({100*n_low/len(cluster_labels):.1f}%)")
        logger.info(f"  Cluster 1 (high): {n_high} cells ({100*n_high/len(cluster_labels):.1f}%)")

        # Store metadata in adata.uns
        adata.uns[f"{column_prefix}_cutoff"] = {
            "metagene_method": metagene_method,
            "feature_columns": valid_columns,
            "threshold": threshold,
            "n_high": n_high,
            "n_low": n_low,
            **self.params_
        }

        return {
            "threshold": threshold,
            "n_high": n_high,
            "n_low": n_low,
            "n_valid": len(metagene_scores),
            "n_total": len(adata),
            "columns_created": [score_col, prob_col, cluster_col],
            **self.params_
        }

    def plot_results(
        self,
        adata: sc.AnnData,
        column_prefix: str = "multivariate",
        figsize: Tuple[int, int] = (14, 5),
        save_path: Optional[str] = None,
    ) -> plt.Figure:
        """
        Plot diagnostic visualization.

        Args:
            adata: AnnData with classifier results
            column_prefix: Prefix used when fitting
            figsize: Figure size
            save_path: Path to save figure (optional)

        Returns:
            matplotlib Figure
        """
        score_col = f"{column_prefix}_metagene_score"
        cluster_col = f"{column_prefix}_cluster"

        if score_col not in adata.obs.columns:
            raise ValueError(f"Column {score_col} not found. Run fit_transform first.")

        valid_mask = adata.obs[cluster_col] >= 0
        scores = adata.obs.loc[valid_mask, score_col].values
        clusters = adata.obs.loc[valid_mask, cluster_col].values

        fig, axes = plt.subplots(1, 3, figsize=figsize)

        # 1. Score distribution with threshold
        ax = axes[0]
        ax.hist(scores, bins=100, alpha=0.7, color='steelblue', edgecolor='white')
        ax.axvline(self.threshold_, color='red', linestyle='--', linewidth=2,
                   label=f'Threshold: {self.threshold_:.3f}')
        ax.set_xlabel('Metagene Score')
        ax.set_ylabel('Count')
        ax.set_title(f'Score Distribution ({self.metagene_method_})')
        ax.legend()

        # 2. Cluster pie chart
        ax = axes[1]
        n_high = np.sum(clusters == 1)
        n_low = np.sum(clusters == 0)
        sizes = [n_low, n_high]
        labels = [f'Low\n{n_low} ({100*n_low/len(clusters):.1f}%)',
                  f'High\n{n_high} ({100*n_high/len(clusters):.1f}%)']
        colors = ['#3498db', '#e74c3c']
        ax.pie(sizes, labels=labels, colors=colors, autopct='', startangle=90)
        ax.set_title('Cluster Distribution')

        # 3. Score by cluster (box plot)
        ax = axes[2]
        low_scores = scores[clusters == 0]
        high_scores = scores[clusters == 1]
        bp = ax.boxplot([low_scores, high_scores], labels=['Low (0)', 'High (1)'],
                       patch_artist=True)
        bp['boxes'][0].set_facecolor('#3498db')
        bp['boxes'][1].set_facecolor('#e74c3c')
        ax.axhline(self.threshold_, color='green', linestyle='--', alpha=0.7,
                   label=f'Threshold')
        ax.set_xlabel('Cluster')
        ax.set_ylabel('Metagene Score')
        ax.set_title('Score by Cluster')
        ax.legend()

        plt.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info(f"Figure saved to {save_path}")

        return fig


# =============================================================================
# Convenience Functions
# =============================================================================

def classify_cells(
    adata: sc.AnnData,
    markers: List[str],
    metagene_method: MetageneMethod = "shifted_geometric_mean",
    threshold_method: ThresholdMethod = "ks",
    pseudocount: float = 0.1,
    column_prefix: str = "multivariate",
    **kwargs
) -> Dict[str, Any]:
    """
    Convenience function to classify cells in one call.

    Args:
        adata: AnnData object (modified in place)
        markers: List of marker genes
        metagene_method: Method for combining markers
        threshold_method: Method for threshold detection
        pseudocount: Pseudocount for shifted_geometric_mean
        column_prefix: Prefix for new columns
        **kwargs: Additional arguments passed to fit_transform

    Returns:
        Dictionary with results

    Example:
        >>> results = classify_cells(
        ...     adata,
        ...     markers=["GeneA", "GeneB", "GeneC"],
        ...     column_prefix="my_population"
        ... )
    """
    classifier = MultivariateMetageneClassifier()
    return classifier.fit_transform(
        adata,
        markers=markers,
        metagene_method=metagene_method,
        threshold_method=threshold_method,
        pseudocount=pseudocount,
        column_prefix=column_prefix,
        **kwargs
    )


# =============================================================================
# Command Line Interface
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Multivariate Metagene Classifier",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    # Basic usage with default settings (shifted_geometric_mean + KS)
    python multivariate_classifier.py --input data.h5ad --markers GeneA,GeneB,GeneC

    # Use GMM threshold method
    python multivariate_classifier.py --input data.h5ad --markers GeneA,GeneB,GeneC \\
        --threshold-method gmm --probability-cutoff 0.3

    # Use arithmetic mean (more permissive)
    python multivariate_classifier.py --input data.h5ad --markers GeneA,GeneB,GeneC \\
        --metagene-method arithmetic_mean

    # Save results to new file
    python multivariate_classifier.py --input data.h5ad --markers GeneA,GeneB,GeneC \\
        --output classified.h5ad --plot results.png
        """
    )

    parser.add_argument("--input", "-i", required=True, help="Input H5AD file")
    parser.add_argument("--output", "-o", help="Output H5AD file (default: overwrites input)")
    parser.add_argument("--markers", "-m", required=True,
                       help="Comma-separated list of marker genes")
    parser.add_argument("--metagene-method", default="shifted_geometric_mean",
                       choices=["geometric_mean", "shifted_geometric_mean",
                               "arithmetic_mean", "median", "minimum"],
                       help="Method for combining markers (default: shifted_geometric_mean)")
    parser.add_argument("--threshold-method", default="ks",
                       choices=["ks", "gmm"],
                       help="Method for threshold detection (default: ks)")
    parser.add_argument("--pseudocount", type=float, default=0.1,
                       help="Pseudocount for shifted_geometric_mean (default: 0.1)")
    parser.add_argument("--background-quantile", type=float, default=0.5,
                       help="KS: fraction for background estimation (default: 0.5)")
    parser.add_argument("--probability-cutoff", type=float, default=0.3,
                       help="GMM: P(high) threshold for classification (default: 0.3)")
    parser.add_argument("--column-prefix", default="multivariate",
                       help="Prefix for new columns (default: multivariate)")
    parser.add_argument("--plot", help="Path to save diagnostic plot")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")

    args = parser.parse_args()

    # Load data
    logger.info(f"Loading {args.input}")
    adata = sc.read_h5ad(args.input)
    logger.info(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")

    # Parse markers
    markers = [m.strip() for m in args.markers.split(",")]
    logger.info(f"Markers: {markers}")

    # Run classifier
    classifier = MultivariateMetageneClassifier()
    results = classifier.fit_transform(
        adata,
        markers=markers,
        metagene_method=args.metagene_method,
        threshold_method=args.threshold_method,
        pseudocount=args.pseudocount,
        background_quantile=args.background_quantile,
        probability_cutoff=args.probability_cutoff,
        column_prefix=args.column_prefix,
        seed=args.seed,
    )

    # Print results
    print("\n" + "="*60)
    print("CLASSIFICATION RESULTS")
    print("="*60)
    print(f"Metagene method:   {args.metagene_method}")
    print(f"Threshold method:  {args.threshold_method}")
    print(f"Threshold:         {results['threshold']:.4f}")
    print(f"High-expression:   {results['n_high']} cells ({100*results['n_high']/results['n_valid']:.1f}%)")
    print(f"Low-expression:    {results['n_low']} cells ({100*results['n_low']/results['n_valid']:.1f}%)")
    print("="*60)

    # Save plot
    if args.plot:
        classifier.plot_results(adata, args.column_prefix, save_path=args.plot)

    # Save output
    output_path = args.output or args.input
    adata.write_h5ad(output_path)
    logger.info(f"Saved results to {output_path}")

    # Show new columns
    print(f"\nNew columns added to adata.obs:")
    for col in results['columns_created']:
        print(f"  - {col}")
