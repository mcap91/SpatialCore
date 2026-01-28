"""Internal threshold detection algorithms for multivariate cutoff.

This module contains pure NumPy/SciPy implementations without AnnData
dependencies, used by the public compute_multivariate_cutoff function.

These functions are internal and not part of the public API.
"""

from typing import Dict, Any, Literal, Tuple

import numpy as np
from scipy.stats import norm
from sklearn.mixture import GaussianMixture


MetageneMethod = Literal[
    "shifted_geometric_mean",
    "geometric_mean",
    "arithmetic_mean",
    "median",
    "minimum",
]

ThresholdMethod = Literal["ks", "gmm"]


def compute_metagene_score(
    feature_values: np.ndarray,
    method: MetageneMethod,
    pseudocount: float = 0.1,
) -> np.ndarray:
    """
    Compute metagene score capturing joint elevation across features.

    The metagene score aggregates multiple features into a single value
    that represents "joint expression" - high only when ALL markers are
    elevated, unlike simple sums which treat single-positive cells the
    same as multi-positive cells.

    Parameters
    ----------
    feature_values
        Array of shape (n_samples, n_features) with expression values.
        Should be non-negative.
    method
        Aggregation method:

        - "shifted_geometric_mean": exp(mean(log(x + c))) - c.
          Default choice. Tolerates dropout while requiring joint expression.
        - "geometric_mean": (x1 * x2 * ... * xn)^(1/n).
          Very strict - any zero collapses score to zero.
        - "arithmetic_mean": Simple mean of features.
          Permissive - allows single-positive cells.
        - "median": Middle value across features.
          Requires majority of markers elevated.
        - "minimum": Minimum across features.
          Strictest - score equals weakest marker.
    pseudocount
        Pseudocount for shifted_geometric_mean method. Default: 0.1.
        Larger values increase tolerance to zeros/dropout.

    Returns
    -------
    np.ndarray
        Array of shape (n_samples,) with metagene scores.

    Notes
    -----
    Method comparison for cell with values (10, 0.1, 10):

    +------------------------+-------+----------------------------------+
    | Method                 | Score | Interpretation                   |
    +========================+=======+==================================+
    | geometric_mean         | ~1.0  | Low (one low marker penalizes)   |
    | shifted_geometric_mean | ~2.1  | Moderate (tolerates dropout)     |
    | arithmetic_mean        | ~6.7  | High (averages out low marker)   |
    | median                 | 10.0  | High (ignores outlier)           |
    | minimum                | 0.1   | Very low (limited by weakest)    |
    +------------------------+-------+----------------------------------+
    """
    if method == "geometric_mean":
        eps = 1e-10
        return np.exp(np.mean(np.log(feature_values + eps), axis=1))

    elif method == "shifted_geometric_mean":
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
    scores: np.ndarray,
    background_quantile: float = 0.5,
) -> Tuple[float, np.ndarray, Dict[str, Any]]:
    """
    KS-inspired threshold detection for sparse populations.

    Finds where the empirical distribution deviates maximally from the
    estimated background distribution. Designed for sparse populations
    where most cells are "background" and a small fraction are "signal".

    Parameters
    ----------
    scores
        1D array of metagene scores.
    background_quantile
        Fraction of data (from the low end) used to estimate background
        distribution. Default: 0.5 (lower 50%).

    Returns
    -------
    threshold
        Score at maximum deviation from background expectation.
    deviation_scores
        Normalized deviation scores [0, 1] for all cells.
        0 = at or below threshold, 1 = at maximum score.
    params
        Dictionary with background_mean, background_std, background_quantile.

    Notes
    -----
    Algorithm:

    1. Sort scores and estimate background from lower quantile
    2. Fit normal distribution to background: N(mu_bg, sigma_bg)
    3. Compute empirical CDF of all data
    4. Compute expected CDF if all data were background
    5. D = empirical_CDF - expected_CDF (excess over expected)
    6. Threshold = score at max(D)

    This is inspired by the Kolmogorov-Smirnov test but adapted for
    one-sided "signal detection" where we expect an excess of high values.

    Best for sparse populations (<5% to 30% of cells).
    """
    sorted_scores = np.sort(scores)
    n = len(sorted_scores)

    # Estimate background from lower quantile
    bg_cutoff_idx = int(n * background_quantile)
    bg_cutoff_idx = max(bg_cutoff_idx, 10)  # At least 10 samples
    background_scores = sorted_scores[:bg_cutoff_idx]
    bg_mean = float(np.mean(background_scores))
    bg_std = float(np.std(background_scores))

    # Handle zero variance edge case
    if bg_std < 1e-10:
        q25, q75 = np.percentile(sorted_scores, [25, 75])
        iqr = q75 - q25
        if iqr > 1e-10:
            bg_std = float(iqr / 1.35)
        else:
            data_range = sorted_scores[-1] - sorted_scores[0]
            bg_std = float(max(data_range * 0.1, 1e-6))

    # Compute empirical CDF
    empirical_cdf = np.arange(1, n + 1) / n

    # Compute expected CDF (if all data were background)
    expected_cdf = norm.cdf(sorted_scores, loc=bg_mean, scale=bg_std)

    # D = excess over expected (positive where we have "extra" high values)
    D = empirical_cdf - expected_cdf

    # Threshold at maximum deviation
    max_idx = np.argmax(D)
    threshold = float(sorted_scores[max_idx])

    # Sanity check: threshold should be above background mean
    if threshold <= bg_mean:
        threshold = float(np.percentile(sorted_scores, 90))

    # Compute deviation scores for all cells
    max_score = sorted_scores[-1]
    score_range = max_score - threshold
    if score_range < 1e-10:
        score_range = 1e-10

    deviation_scores = np.clip((scores - threshold) / score_range, 0, 1)

    params = {
        "background_mean": bg_mean,
        "background_std": bg_std,
        "background_quantile": background_quantile,
    }

    return threshold, deviation_scores, params


def threshold_gmm(
    scores: np.ndarray,
    probability_cutoff: float = 0.3,
    n_components: int = 2,
    random_state: int = 42,
) -> Tuple[float, np.ndarray, np.ndarray, Dict[str, Any]]:
    """
    GMM-based threshold detection.

    Fits an N-component Gaussian Mixture Model to find the threshold
    separating low and high populations. Returns posterior probabilities
    for soft assignment.

    Parameters
    ----------
    scores
        1D array of metagene scores.
    probability_cutoff
        P(high) threshold for cluster assignment. Default: 0.3.
        Cells with P(high|score) > probability_cutoff are assigned to
        cluster 1 (high).
    n_components
        Number of GMM components. Default: 2.

        - 2: Standard bimodal separation (low vs high)
        - 3: Trimodal for spatial data (dropout/moderate/high).
          Threshold is set between component 0 (background) and
          component 1 (first signal peak). This handles spatial data
          where marker expression is often trimodal.
    random_state
        Random seed for GMM initialization.

    Returns
    -------
    threshold
        Score at the threshold between low and high populations.
        For n_components=2: P(high) = P(low) = 0.5.
        For n_components=3: Midpoint between component 0 and 1 means.
    cluster_labels
        Binary cluster assignments: 0 (low) or 1 (high).
    probability_high
        Posterior probability of being in "high" component(s) for all cells.
    params
        Dictionary with GMM parameters (means, stds, weights).

    Notes
    -----
    Algorithm for n_components=2:

    1. Fit 2-component GMM via EM algorithm
    2. Identify high component (higher mean)
    3. Find threshold where P(high|x) = P(low|x) = 0.5
    4. Compute posterior P(high|score) for all cells
    5. Assign clusters based on probability_cutoff

    Algorithm for n_components=3 (trimodal spatial data):

    1. Fit 3-component GMM via EM algorithm
    2. Sort components by mean: [background, moderate, high]
    3. Threshold = midpoint between component 0 and 1 means
    4. "High" = cells above threshold (in components 1 or 2)
    5. P(high) = P(component 1) + P(component 2)

    The 3-component mode handles spatial data where marker expression
    is often trimodal: peak 0 (~0) = dropout, peak 1 (~1-2) = moderate
    expression, peak 2 (~4-5) = high expression.

    Best when both populations are substantial (>10% each).
    Provides soft assignments (probabilities) unlike KS.
    """
    gmm = GaussianMixture(
        n_components=n_components,
        random_state=random_state,
        n_init=10,
        covariance_type="full",
    )
    gmm.fit(scores.reshape(-1, 1))

    # Get component statistics
    component_means = gmm.means_.flatten()
    component_stds = np.sqrt(gmm.covariances_.flatten())

    # Sort components by mean (ascending)
    sorted_indices = np.argsort(component_means)
    sorted_means = component_means[sorted_indices]

    if n_components == 2:
        # Standard bimodal: threshold at P(high) = 0.5
        high_component = int(np.argmax(component_means))
        low_component = 1 - high_component

        low_mean = component_means[low_component]
        high_mean = component_means[high_component]

        # Find threshold at P(high) = 0.5
        x_range = np.linspace(low_mean, high_mean, 1000)
        posteriors = gmm.predict_proba(x_range.reshape(-1, 1))
        diff = posteriors[:, high_component] - 0.5
        cross_idx = np.where(np.diff(np.sign(diff)))[0]

        if len(cross_idx) > 0:
            threshold = float(x_range[cross_idx[0]])
        else:
            # Fallback to midpoint if no crossing found
            threshold = float((low_mean + high_mean) / 2)

        # Compute posteriors for all cells
        posteriors_all = gmm.predict_proba(scores.reshape(-1, 1))
        probability_high = posteriors_all[:, high_component]

    else:  # n_components >= 3 (trimodal mode)
        # For trimodal: threshold between component 0 (background) and 1 (first signal)
        # Components are sorted by mean: [background, moderate, high]
        background_idx = sorted_indices[0]
        signal_idx = sorted_indices[1]  # First signal component

        # Threshold = midpoint between background and first signal
        threshold = float((sorted_means[0] + sorted_means[1]) / 2)

        # Compute posteriors for all cells
        posteriors_all = gmm.predict_proba(scores.reshape(-1, 1))

        # P(high) = sum of probabilities for all non-background components
        # (i.e., all components except the lowest-mean one)
        signal_components = sorted_indices[1:]  # All except background
        probability_high = posteriors_all[:, signal_components].sum(axis=1)

        high_component = int(signal_idx)  # Primary signal component

    # Assign clusters based on probability cutoff
    cluster_labels = (probability_high > probability_cutoff).astype(int)

    params = {
        "gmm_means": component_means.tolist(),
        "gmm_stds": component_stds.tolist(),
        "gmm_weights": gmm.weights_.tolist(),
        "n_components": n_components,
        "sorted_component_indices": sorted_indices.tolist(),
        "high_component_idx": high_component,
        "probability_cutoff": probability_cutoff,
        "gmm_model": gmm,  # Fitted model for reuse when downsampling
    }

    return threshold, cluster_labels, probability_high, params
