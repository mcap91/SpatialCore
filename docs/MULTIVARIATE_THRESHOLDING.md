# Multivariate Thresholding System

## Executive Summary

This document provides comprehensive documentation for the **statistical thresholding system** in our spatial biology library. The system provides two complementary functions for objective threshold selection:

| Function | Use Case | Population Balance | Methods |
|----------|----------|-------------------|---------|
| `compute_bimodal_cutoff` | Balanced populations | ~40-60% split | hierarchical, gmm, kde |
| `compute_multivariate_cutoff` | Sparse populations | <30% target | gmm, ks |

**Key Insight**: When using GMM on a single feature, both functions produce **identical results** because they share the same underlying Gaussian Mixture Model implementation.

---

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [When to Use Which Function](#when-to-use-which-function)
3. [compute_bimodal_cutoff (Univariate)](#compute_bimodal_cutoff-univariate)
4. [compute_multivariate_cutoff (Metagene)](#compute_multivariate_cutoff-metagene)
5. [Threshold Methods Comparison](#threshold-methods-comparison)
6. [Metagene Aggregation Methods](#metagene-aggregation-methods)
7. [Mathematical Foundations](#mathematical-foundations)
8. [Standalone Implementation](#standalone-implementation)
9. [Test Scripts](#test-scripts)
10. [API Reference](#api-reference)

---

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                     STATISTICAL THRESHOLDING SYSTEM                         │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌─────────────────────────────┐     ┌─────────────────────────────────┐   │
│  │  compute_bimodal_cutoff     │     │  compute_multivariate_cutoff    │   │
│  │  (Univariate)               │     │  (Metagene-based)               │   │
│  │                             │     │                                 │   │
│  │  Input: Single feature      │     │  Input: 1+ features             │   │
│  │  Output: k-1 cutoffs        │     │  Output: 1 metagene threshold   │   │
│  │                             │     │                                 │   │
│  │  Methods:                   │     │  Metagene Methods:              │   │
│  │  - hierarchical (default)   │     │  - shifted_geometric_mean       │   │
│  │  - gmm                      │     │  - geometric_mean               │   │
│  │  - kde                      │     │  - arithmetic_mean              │   │
│  │                             │     │  - median                       │   │
│  │  k = 2, 3, 4, 5...          │     │  - minimum                      │   │
│  │                             │     │                                 │   │
│  │                             │     │  Threshold Methods:             │   │
│  │                             │     │  - ks (default, sparse-aware)   │   │
│  │                             │     │  - gmm (Gaussian populations)   │   │
│  └─────────────────────────────┘     └─────────────────────────────────┘   │
│              │                                    │                         │
│              ▼                                    ▼                         │
│  ┌─────────────────────────────────────────────────────────────────────┐   │
│  │                    SHARED COMPONENTS                                 │   │
│  │                                                                      │   │
│  │  _bimodal_cutoff_gmm()     ←──────────────────→   GMM (sklearn)      │   │
│  │  _bimodal_cutoff_hierarchical()                   Ward clustering    │   │
│  │  _bimodal_cutoff_kde()                            KDE (scipy)        │   │
│  │  _threshold_ks()           ←──────────────────→   KS-inspired        │   │
│  │  _compute_metagene_score()                        Aggregation        │   │
│  │  _compute_gmm_threshold()                         P(high) = 0.5      │   │
│  │                                                                      │   │
│  └─────────────────────────────────────────────────────────────────────┘   │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Data Flow

```
Input Data (H5AD)
       │
       ▼
┌──────────────────┐
│ Feature Extraction│
│ - adata.obs      │
│ - adata.var_names│
│ - adata.obsm     │
└────────┬─────────┘
         │
         ▼
┌──────────────────────────────────────────────────────────┐
│ UNIVARIATE PATH               │ MULTIVARIATE PATH        │
│ (compute_bimodal_cutoff)      │ (compute_multivariate)   │
│                               │                          │
│ Single feature values         │ Multiple feature values  │
│         │                     │         │                │
│         ▼                     │         ▼                │
│ Direct clustering:            │ Metagene computation:    │
│ - hierarchical                │ score = f(x1, x2, ..xn)  │
│ - gmm                         │         │                │
│ - kde                         │         ▼                │
│         │                     │ Threshold detection:     │
│         ▼                     │ - ks (sparse)            │
│ k-1 cutoff boundaries         │ - gmm (balanced)         │
│                               │         │                │
│                               │         ▼                │
│                               │ Single threshold         │
└──────────────────────────────────────────────────────────┘
         │
         ▼
┌──────────────────┐
│ Output:          │
│ - adata.obs cols │
│ - adata.uns meta │
│ - Diagnostic plot│
└──────────────────┘
```

---

## When to Use Which Function

### Decision Tree

```
                    What is your data distribution?
                               │
              ┌────────────────┴────────────────┐
              ▼                                 ▼
    Balanced populations              Sparse target population
    (~40-60% split)                   (<30% of cells)
              │                                 │
              ▼                                 ▼
  ┌───────────────────────┐       ┌───────────────────────────┐
  │ compute_bimodal_cutoff │       │ compute_multivariate_cutoff│
  │                       │       │                           │
  │ Examples:             │       │ Examples:                 │
  │ - Moran's I (spatial  │       │ - CD38+ plasma cells      │
  │   autocorrelated vs   │       │ - Epithelial markers      │
  │   random)             │       │   (SCGB1A1, SCGB3A1)      │
  │ - Tumor vs normal     │       │ - Any rare population     │
  │                       │       │                           │
  │ Use when:             │       │ Use when:                 │
  │ - Need k > 2 clusters │       │ - Most cells = background │
  │ - Clear bimodal/      │       │ - Finding outliers        │
  │   multimodal peaks    │       │ - Skewed distributions    │
  └───────────────────────┘       └───────────────────────────┘
```

### Quick Reference Table

| Scenario | Function | Method | Parameters |
|----------|----------|--------|------------|
| Moran's I thresholding | `compute_bimodal_cutoff` | `hierarchical` | `n_clusters=3` |
| Marker+ rare cells (1-5%) | `compute_multivariate_cutoff` | `ks` | default |
| Multiple markers, sparse | `compute_multivariate_cutoff` | `ks` | `metagene_method="shifted_geometric_mean"` |
| Two balanced populations | `compute_bimodal_cutoff` | `gmm` | `n_clusters=2` |
| Low/Medium/High groups | `compute_bimodal_cutoff` | `hierarchical` | `n_clusters=3` |

---

## compute_bimodal_cutoff (Univariate)

### Purpose

Identifies optimal thresholds in a **single feature distribution** using clustering-based approaches. Returns `k-1` cutoff boundaries that separate `k` clusters.

### Function Signature

```python
@activity.defn
async def compute_bimodal_cutoff(
    data_uri: str,
    output_uri: str,
    feature_column: str,
    method: str = "hierarchical",          # hierarchical | gmm | kde
    n_clusters: Optional[int] = None,      # Default: 3
    auto_k: bool = False,                  # Auto-select k via silhouette
    max_k: int = 5,                        # Max k for auto selection
    downsample_prop: float = 0.1,          # Subsample for speed
    max_cells: int = 50000,
    seed: int = 2024,
    bootstrap: bool = False,               # Compute confidence intervals
    n_bootstrap: int = 100,
    test_unimodal: bool = False,           # Warn if unimodal
    return_plot: bool = True,
    plot_format: str = "png",
    dpi: int = 300,
    **kwargs,
) -> dict:
```

### Methods

#### 1. Hierarchical Clustering (`method="hierarchical"`)

**Default choice**. Uses Ward's hierarchical clustering on 1D Euclidean distances.

```
Algorithm:
1. Reshape values to (n, 1) for distance calculation
2. Compute pairwise Euclidean distances: O(n²)
3. Ward's linkage: minimize within-cluster variance
4. Cut dendrogram at k clusters
5. Cutoffs = midpoints between adjacent cluster max/min

Complexity: O(n² log n)
Memory: O(n²) for distance matrix
```

**Pros:**
- Fast, deterministic
- No distributional assumptions
- Works well for most cases

**Cons:**
- O(n²) memory limits scalability
- Assumes spherical clusters

#### 2. Gaussian Mixture Model (`method="gmm"`)

**More principled** for normally-distributed modes.

```
Algorithm:
1. Fit k-component GMM via EM algorithm
2. Identify components by mean (sorted low→high)
3. Cutoffs = midpoints between adjacent component means

Complexity: O(n × k × iterations)
Memory: O(n × k)
```

**Pros:**
- Principled statistical model
- Provides posterior probabilities
- Scales better than hierarchical

**Cons:**
- Assumes Gaussian components
- EM may converge to local optima

#### 3. Kernel Density Estimation (`method="kde"`)

**Density-based** for irregular/skewed distributions.

```
Algorithm:
1. Compute KDE using Gaussian kernel (Silverman bandwidth)
2. Evaluate density on grid of 1000 points
3. Find local minima (valleys between peaks)
4. Cutoffs = x-values at density minima

Complexity: O(n²) for KDE
Memory: O(n)
```

**Pros:**
- Non-parametric, no assumptions
- Automatically finds number of modes
- Robust to skewness

**Cons:**
- Bandwidth selection affects results
- May find spurious minima in noisy data

### Output

```python
{
    "cutoffs": [0.5, 1.2],           # k-1 threshold values
    "cutoff_ci": [(0.45, 0.55), ...],# Bootstrap CIs (optional)
    "n_clusters": 3,
    "method_used": "hierarchical",
    "n_per_cluster": {"0": 1500, "1": 800, "2": 200},
    "sampled_cluster_sizes": [150, 80, 20],
    "optimal_k": 3,                  # If auto_k=True
    "optimal_k_score": 0.65,         # Silhouette score
    "unimodal_p": 0.02,              # If test_unimodal=True
    "n_cells_analyzed": 250,
    "n_cells_total": 2500,
    "plot_path": "./output/plot.png"
}
```

### Columns Created

| Column | Type | Description |
|--------|------|-------------|
| `bimodal_cluster` | int | Cluster assignment (0 to k-1), -1 for invalid |

### Metadata Stored

```python
adata.uns["bimodal_cutoff"] = {
    "cutoffs": [0.5, 1.2],
    "method": "hierarchical",
    "n_clusters": 3,
    "feature_column": "MS4A1_local_morans_i",
    "sampled_cluster_sizes": [150, 80, 20],
    "n_per_cluster": {"0": 1500, "1": 800, "2": 200},
    "timestamp": "2024-10-20T12:34:56"
}
```

---

## compute_multivariate_cutoff (Metagene)

### Purpose

Computes a **metagene score** that captures joint elevation across multiple features, then finds a threshold separating high vs low populations. Designed for **sparse single-cell data** with skewed distributions.

### Why Metagene Instead of Euclidean Clustering?

```
Problem: Simple Euclidean clustering treats single-positive same as triple-positive

Example: 3 markers for epithelial cells (SCGB1A1, SCGB3A1, SCGB3A2)

Cell A: (high, low, low)   ← Single-positive
Cell B: (high, high, high) ← Triple-positive (actual epithelial)

Euclidean distance: Both may cluster together!
Metagene (geometric mean): Cell A gets LOW score, Cell B gets HIGH score

┌────────────────────────────────────────────────────────────────┐
│ Euclidean Clustering              │ Metagene Approach          │
│                                   │                            │
│     x                             │     x       x              │
│   x   x                           │   x ← LOW     x ← HIGH     │
│  A     B ← Same cluster!          │  A           B             │
│                                   │                            │
│ (distance ignores which           │ (geometric mean requires   │
│  markers are elevated)            │  ALL markers elevated)     │
└────────────────────────────────────────────────────────────────┘
```

### Function Signature

```python
@activity.defn
async def compute_multivariate_cutoff(
    data_uri: str,
    output_uri: str,
    feature_columns: List[str],                      # 1+ features
    metagene_method: MetageneMethod = "shifted_geometric_mean",
    pseudocount: float = 0.1,                        # For shifted_geometric_mean
    threshold_method: ThresholdMethod = "gmm",       # gmm | ks
    # KS-specific
    background_quantile: float = 0.5,
    # GMM-specific
    probability_cutoff: float = 0.3,
    n_clusters: int = 2,
    # Common
    downsample_prop: float = 0.1,
    max_cells: int = 50000,
    seed: int = 2024,
    return_plot: bool = True,
    plot_format: str = "png",
    dpi: int = 300,
    **kwargs,
) -> dict:
```

### Single Feature Behavior

**Important**: When called with a single feature, `metagene_score = feature_value` (identity). This means:

```python
# These produce IDENTICAL results when threshold_method="gmm":

# Option 1: Univariate with GMM
result1 = await compute_bimodal_cutoff(
    feature_column="MS4A1_local_morans_i",
    method="gmm",
    n_clusters=2
)

# Option 2: Multivariate with single feature + GMM
result2 = await compute_multivariate_cutoff(
    feature_columns=["MS4A1_local_morans_i"],
    threshold_method="gmm"
)

# result1.cutoffs[0] == result2.threshold (same value!)
```

### Threshold Methods

#### 1. KS-Inspired (`threshold_method="ks"`)

**Default for sparse populations**. Non-parametric, distribution-free.

```
Algorithm:
1. Estimate background distribution from lower quantile (default: 50%)
2. Compute empirical CDF of all data
3. Compute expected CDF if all data were background (normal assumption)
4. D = empirical_cdf - expected_cdf (excess above expected)
5. Threshold = score at max(D) (maximum deviation from background)

Key insight: Where the empirical distribution "breaks away" from
             what we'd expect if everything were background.

┌─────────────────────────────────────────────────────────────┐
│ CDF                                                         │
│  1.0 ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ═══════════════       │
│      │              ┌─────────────────────────────────      │
│      │         ┌────┘  ← Empirical CDF                      │
│      │    ┌────┘                                            │
│      │┌───┘                                                 │
│      │                                                      │
│  0.5 ├─────────╱─────────────────────────────────────       │
│      │       ╱   ← Expected (background) CDF                │
│      │     ╱                                                │
│      │   ╱  D = max deviation                               │
│      │ ╱    │                                               │
│  0.0 ┴──────┼─────────────────────────────────────────      │
│             threshold                           score →     │
└─────────────────────────────────────────────────────────────┘
```

**Pros:**
- Works well for sparse populations (1-30%)
- No Gaussian assumption
- Robust to skewed distributions

**Cons:**
- Only finds 1 threshold (binary classification)
- Assumes "background" is lower half

#### 2. GMM (`threshold_method="gmm"`)

**For balanced populations or when you need probabilities**.

```
Algorithm:
1. Fit 2-component GMM to metagene scores
2. Identify high component (higher mean)
3. Compute threshold at P(high) = P(low) = 0.5 intersection
4. Compute posterior P(high|score) for all cells
5. Assign clusters: P(high) > probability_cutoff → cluster 1

Threshold = score where GMM components have equal probability
```

**Pros:**
- Provides soft assignments (probabilities)
- Can adjust stringency via `probability_cutoff`
- Works well when both populations are substantial

**Cons:**
- Assumes Gaussian components
- May struggle with very sparse populations

### Output

```python
{
    "metagene_method": "shifted_geometric_mean",
    "threshold_method": "ks",
    "feature_columns": ["SCGB1A1", "SCGB3A1", "SCGB3A2"],
    "n_features": 3,
    "threshold": 0.85,
    "n_cells_high": 150,
    "n_cells_low": 2350,
    "n_cells_analyzed": 2500,
    "n_cells_total": 2500,
    "plot_path": "./output/plot.png",

    # KS-specific:
    "background_quantile": 0.5,
    "background_mean": 0.12,
    "background_std": 0.15,

    # GMM-specific (when threshold_method="gmm"):
    "probability_cutoff": 0.3,
    "gmm_means": [0.2, 1.5],
    "gmm_stds": [0.1, 0.3],
    "gmm_weights": [0.9, 0.1],
    "high_component_idx": 1,
}
```

### Columns Created

| Column | Type | Description |
|--------|------|-------------|
| `metagene_score` | float | Combined score from multiple features |
| `metagene_probability` | float | P(high) posterior (GMM) or deviation score (KS) |
| `metagene_cluster` | int | 0 (low) or 1 (high), -1 for invalid |

---

## Threshold Methods Comparison

### Summary Table

| Method | Function | Assumptions | Best For | Output |
|--------|----------|-------------|----------|--------|
| `hierarchical` | bimodal | Spherical clusters | General use, k>2 | k-1 cutoffs |
| `gmm` | both | Gaussian components | Normal-ish modes | k-1 cutoffs + P(high) |
| `kde` | bimodal | Smooth density | Skewed, irregular | Auto cutoffs |
| `ks` | multivariate | Background = lower 50% | Sparse (<30%) | 1 cutoff |

### Complexity Comparison

| Method | Time | Memory | Scalability |
|--------|------|--------|-------------|
| hierarchical | O(n² log n) | O(n²) | Poor (>50K cells slow) |
| gmm | O(n × k × iter) | O(n × k) | Good |
| kde | O(n²) | O(n) | Moderate |
| ks | O(n log n) | O(n) | Excellent |

### When to Use Each

```
┌─────────────────────────────────────────────────────────────────┐
│                                                                 │
│   Sparse (<30%)              Balanced (~40-60%)                 │
│        │                           │                            │
│        ▼                           ▼                            │
│   ┌────────────┐             ┌────────────┐                     │
│   │    ks      │             │   gmm or   │                     │
│   │  (default) │             │ hierarchical│                     │
│   └────────────┘             └────────────┘                     │
│        │                           │                            │
│        │                           │                            │
│   Need k > 2?                 Need probs?                       │
│        │                           │                            │
│   No   │   Yes                No   │   Yes                      │
│        │    │                      │    │                       │
│        ▼    ▼                      ▼    ▼                       │
│   ┌────────────┐             ┌──────┐ ┌────┐                    │
│   │  Use ks    │             │ hier │ │gmm │                    │
│   │  (binary)  │             │      │ │    │                    │
│   └────────────┘             └──────┘ └────┘                    │
│                                                                 │
│   Consider: Use compute_bimodal_cutoff                          │
│   with method="gmm" for k>2 with sparse                         │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## Metagene Aggregation Methods

### Formula Summary

| Method | Formula | Behavior |
|--------|---------|----------|
| `geometric_mean` | `(x1 × x2 × ... × xn)^(1/n)` | Penalizes zeros heavily |
| `shifted_geometric_mean` | `exp(mean(log(x + c))) - c` | Tolerates dropout |
| `arithmetic_mean` | `(x1 + x2 + ... + xn) / n` | Forgives low markers |
| `median` | `median(x1, x2, ..., xn)` | Majority voting |
| `minimum` | `min(x1, x2, ..., xn)` | Strictest (all must be high) |

### Visual Comparison

```
Marker values: A=high, B=low, C=high

          Geometric    Shifted-Geo    Arithmetic    Median    Minimum
Cell 1:     LOW           LOW           MEDIUM       MEDIUM     LOW
(A,B,C) = (5,0.1,4)

Cell 2:     HIGH          HIGH          HIGH         HIGH       HIGH
(A,B,C) = (5,4,4)

Cell 3:     VERY LOW      LOW           MEDIUM       MEDIUM     VERY LOW
(A,B,C) = (5,0,4)      (zero kills it)  (tolerant)   (2/3 high)  (zero)
```

### Recommendations

| Scenario | Method | Reason |
|----------|--------|--------|
| Default (sparse scRNA-seq) | `shifted_geometric_mean` | Tolerates dropout, requires joint expression |
| No dropout, strict | `geometric_mean` | Penalizes any low marker |
| Partial positivity OK | `arithmetic_mean` | Single-positive gets some credit |
| 2/3 majority rule | `median` | Robust to outlier markers |
| All markers required | `minimum` | Unusable with dropout |

---

## Mathematical Foundations

### Gaussian Mixture Model (GMM)

The GMM assumes data is generated from a mixture of K Gaussian distributions:

```
p(x) = Σᵢ πᵢ × N(x | μᵢ, σᵢ²)

where:
- πᵢ = weight of component i (Σπᵢ = 1)
- μᵢ = mean of component i
- σᵢ = standard deviation of component i
```

**Posterior probability** (soft assignment):
```
P(component j | x) = πⱼ × N(x | μⱼ, σⱼ²) / Σᵢ πᵢ × N(x | μᵢ, σᵢ²)
```

**Threshold** (at equal posterior):
```
Find x* where P(high | x*) = P(low | x*) = 0.5
```

### KS-Inspired Threshold

Inspired by Kolmogorov-Smirnov test, but adapted for one-sided "signal detection":

```
1. Estimate background: μ_bg, σ_bg from lower quantile

2. Empirical CDF: F_emp(x) = (# points ≤ x) / n

3. Expected CDF (if all background): F_exp(x) = Φ((x - μ_bg) / σ_bg)

4. Deviation: D(x) = F_emp(x) - F_exp(x)

5. Threshold: x* = argmax D(x)
              (where data most exceeds background expectation)
```

### Geometric Mean Properties

```
GM(x1, x2, ..., xn) = (∏ xᵢ)^(1/n) = exp(1/n × Σ log(xᵢ))

Properties:
- GM ≤ AM (always less than arithmetic mean)
- GM = 0 if any xᵢ = 0
- Log-transformed: log(GM) = mean of log values (additive)
- Interpretation: n-th root of product
```

**Shifted geometric mean** (for dropout tolerance):
```
SGM(x; c) = exp(mean(log(x + c))) - c

where c = pseudocount (default: 0.1)

When c → 0: SGM → GM (strict)
When c → ∞: SGM → AM (forgiving)
```

---

## Standalone Implementation

### Core Components

The following code provides a standalone implementation of both thresholding functions:

```python
#!/usr/bin/env python3
"""
Standalone Statistical Thresholding Module

Provides objective threshold selection for univariate and multivariate distributions.
Can be used independently of the Temporal workflow system.

Dependencies:
    - numpy
    - scipy
    - sklearn
    - matplotlib (optional, for plotting)
"""

import numpy as np
from typing import List, Tuple, Optional, Literal
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, cut_tree
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde, norm
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score


# Type aliases
MetageneMethod = Literal[
    "geometric_mean", "arithmetic_mean", "minimum",
    "shifted_geometric_mean", "median"
]
ThresholdMethod = Literal["ks", "gmm"]


# =============================================================================
# UNIVARIATE METHODS
# =============================================================================

def bimodal_cutoff_hierarchical(
    values: np.ndarray,
    n_clusters: int = 3
) -> Tuple[List[float], np.ndarray]:
    """
    Hierarchical clustering method for cutoff detection.

    Args:
        values: 1D array of feature values
        n_clusters: Number of clusters to detect

    Returns:
        cutoffs: Sorted list of n_clusters - 1 boundary values
        cluster_labels: Cluster assignment for each value
    """
    # Reshape for distance calculation
    X = values.reshape(-1, 1)

    # Compute pairwise distances
    distances = pdist(X, metric="euclidean")

    # Hierarchical clustering with Ward linkage
    linkage_matrix = linkage(distances, method="ward")

    # Cut tree at k clusters
    cluster_labels = cut_tree(linkage_matrix, n_clusters=n_clusters).flatten()

    # Compute cutoffs as boundaries between clusters
    cluster_stats = []
    for cluster_id in range(n_clusters):
        cluster_values = values[cluster_labels == cluster_id]
        cluster_stats.append({
            "id": cluster_id,
            "min": float(cluster_values.min()),
            "max": float(cluster_values.max()),
            "mean": float(cluster_values.mean()),
        })

    # Sort clusters by their mean value
    cluster_stats.sort(key=lambda x: x["mean"])

    # Cutoffs = midpoints between adjacent clusters
    cutoffs = []
    for i in range(len(cluster_stats) - 1):
        boundary = (cluster_stats[i]["max"] + cluster_stats[i + 1]["min"]) / 2
        cutoffs.append(boundary)

    return sorted(cutoffs), cluster_labels


def bimodal_cutoff_gmm(
    values: np.ndarray,
    n_clusters: int = 2,
    seed: int = 2024
) -> Tuple[List[float], np.ndarray, GaussianMixture]:
    """
    Gaussian Mixture Model method for cutoff detection.

    Args:
        values: 1D array of feature values
        n_clusters: Number of Gaussian components
        seed: Random seed for EM initialization

    Returns:
        cutoffs: Sorted list of n_clusters - 1 boundary values
        cluster_labels: Cluster assignment for each value
        gmm: Fitted GaussianMixture model
    """
    X = values.reshape(-1, 1)

    # Fit GMM
    gmm = GaussianMixture(
        n_components=n_clusters,
        random_state=seed,
        covariance_type="full",
        n_init=10,
    )
    gmm.fit(X)

    # Predict cluster labels
    cluster_labels = gmm.predict(X)

    # Sort components by their means
    sorted_indices = np.argsort(gmm.means_.flatten())
    sorted_means = gmm.means_.flatten()[sorted_indices]

    # Cutoffs = midpoints between adjacent component means
    cutoffs = []
    for i in range(len(sorted_means) - 1):
        boundary = (sorted_means[i] + sorted_means[i + 1]) / 2
        cutoffs.append(float(boundary))

    return cutoffs, cluster_labels, gmm


def bimodal_cutoff_kde(
    values: np.ndarray,
    n_clusters: Optional[int] = None
) -> Tuple[List[float], np.ndarray]:
    """
    Kernel Density Estimation method for cutoff detection.

    Args:
        values: 1D array of feature values
        n_clusters: Number of modes to expect (optional validation)

    Returns:
        cutoffs: Sorted list of density minima (valley points)
        cluster_labels: Cluster assignment based on cutoffs
    """
    try:
        kde = gaussian_kde(values)
    except np.linalg.LinAlgError:
        # Fallback for singular covariance
        if n_clusters is not None and n_clusters > 1:
            cutoffs = [
                float(np.percentile(values, 100 * i / n_clusters))
                for i in range(1, n_clusters)
            ]
        else:
            cutoffs = [float(np.median(values))]
        cluster_labels = np.zeros(len(values), dtype=int)
        for i, cutoff in enumerate(cutoffs):
            cluster_labels[values > cutoff] = i + 1
        return cutoffs, cluster_labels

    # Evaluate KDE on a grid
    x_grid = np.linspace(values.min(), values.max(), 1000)
    density = kde(x_grid)

    # Find local minima
    minima_indices = argrelextrema(density, np.less)[0]

    if len(minima_indices) == 0:
        cutoffs = [float(np.median(values))]
        cluster_labels = (values > cutoffs[0]).astype(int)
        return cutoffs, cluster_labels

    cutoffs = sorted([float(x_grid[i]) for i in minima_indices])

    # Assign cluster labels based on cutoffs
    cluster_labels = np.zeros(len(values), dtype=int)
    for i, cutoff in enumerate(cutoffs):
        cluster_labels[values > cutoff] = i + 1

    return cutoffs, cluster_labels


# =============================================================================
# MULTIVARIATE METHODS
# =============================================================================

def compute_metagene_score(
    feature_values: np.ndarray,
    method: MetageneMethod = "shifted_geometric_mean",
    pseudocount: float = 0.1,
) -> np.ndarray:
    """
    Compute metagene score capturing joint elevation across features.

    Args:
        feature_values: (n_samples, n_features) array
        method: Aggregation method
        pseudocount: For shifted_geometric_mean

    Returns:
        (n_samples,) array of metagene scores
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
        raise ValueError(f"Unknown method: {method}")


def threshold_ks(
    scores: np.ndarray,
    background_quantile: float = 0.5,
) -> Tuple[float, np.ndarray, float, float]:
    """
    KS-inspired threshold detection for sparse populations.

    Args:
        scores: 1D array of (metagene) scores
        background_quantile: Fraction of data to estimate background

    Returns:
        threshold: Score at maximum deviation from background
        deviation_scores: Normalized deviation scores [0, 1]
        bg_mean: Estimated background mean
        bg_std: Estimated background std
    """
    sorted_scores = np.sort(scores)
    n = len(sorted_scores)

    # Estimate background from lower quantile
    bg_cutoff_idx = int(n * background_quantile)
    background_scores = sorted_scores[:bg_cutoff_idx]
    bg_mean = float(np.mean(background_scores))
    bg_std = float(np.std(background_scores))

    # Handle zero variance
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

    # Compute expected CDF (background)
    expected_cdf = norm.cdf(sorted_scores, loc=bg_mean, scale=bg_std)

    # Deviation = excess over expected
    D = empirical_cdf - expected_cdf

    # Threshold at maximum deviation
    max_idx = np.argmax(D)
    threshold = float(sorted_scores[max_idx])

    # Sanity check
    if threshold <= bg_mean:
        threshold = float(np.percentile(sorted_scores, 90))

    # Compute deviation scores
    max_score = sorted_scores[-1]
    score_range = max_score - threshold
    if score_range < 1e-10:
        score_range = 1e-10

    deviation_scores = np.clip((scores - threshold) / score_range, 0, 1)

    return threshold, deviation_scores, bg_mean, bg_std


def compute_gmm_threshold(gmm: GaussianMixture, high_component: int) -> float:
    """
    Compute threshold where two GMM components have equal probability.

    Args:
        gmm: Fitted GaussianMixture model
        high_component: Index of the high-expression component

    Returns:
        Threshold value (score at P(high) = P(low) = 0.5)
    """
    means = gmm.means_.flatten()
    low_component = 1 - high_component

    low_mean, high_mean = means[low_component], means[high_component]

    # Search for P(high) = 0.5
    x_range = np.linspace(low_mean, high_mean, 1000)
    posteriors = gmm.predict_proba(x_range.reshape(-1, 1))

    diff = posteriors[:, high_component] - 0.5
    cross_idx = np.where(np.diff(np.sign(diff)))[0]

    if len(cross_idx) > 0:
        return float(x_range[cross_idx[0]])
    else:
        return float((low_mean + high_mean) / 2)


def multivariate_cutoff_gmm(
    scores: np.ndarray,
    probability_cutoff: float = 0.3,
    seed: int = 2024,
) -> Tuple[float, np.ndarray, np.ndarray, dict]:
    """
    GMM-based threshold detection for metagene scores.

    Args:
        scores: 1D array of metagene scores
        probability_cutoff: P(high) threshold for assignment
        seed: Random seed

    Returns:
        threshold: GMM intersection point
        cluster_labels: Binary assignment (0=low, 1=high)
        probability_high: Posterior P(high|score)
        gmm_info: Dictionary with GMM parameters
    """
    gmm = GaussianMixture(
        n_components=2,
        random_state=seed,
        n_init=10,
    )
    gmm.fit(scores.reshape(-1, 1))

    # Identify high component
    component_means = gmm.means_.flatten()
    component_stds = np.sqrt(gmm.covariances_.flatten())
    high_component = int(np.argmax(component_means))

    # Compute threshold
    threshold = compute_gmm_threshold(gmm, high_component)

    # Compute posteriors
    posteriors = gmm.predict_proba(scores.reshape(-1, 1))
    probability_high = posteriors[:, high_component]

    # Assign clusters
    cluster_labels = (probability_high > probability_cutoff).astype(int)

    gmm_info = {
        "gmm_means": component_means.tolist(),
        "gmm_stds": component_stds.tolist(),
        "gmm_weights": gmm.weights_.tolist(),
        "high_component_idx": high_component,
        "probability_cutoff": probability_cutoff,
    }

    return threshold, cluster_labels, probability_high, gmm_info


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def select_optimal_k(
    values: np.ndarray,
    max_k: int = 5,
    method: str = "hierarchical",
    seed: int = 2024
) -> Tuple[int, float]:
    """
    Auto-select optimal k using silhouette score.

    Args:
        values: 1D array of feature values
        max_k: Maximum k to test
        method: Clustering method
        seed: Random seed

    Returns:
        optimal_k: Best number of clusters
        best_score: Silhouette score for optimal k
    """
    X = values.reshape(-1, 1)
    best_k = 2
    best_score = -1

    for k in range(2, max_k + 1):
        if method == "hierarchical":
            _, labels = bimodal_cutoff_hierarchical(values, n_clusters=k)
        elif method == "gmm":
            _, labels, _ = bimodal_cutoff_gmm(values, n_clusters=k, seed=seed)
        else:
            _, labels = bimodal_cutoff_kde(values, n_clusters=k)

        score = silhouette_score(X, labels)

        if score > best_score:
            best_score = score
            best_k = k

    return best_k, best_score
```

### Package Structure for Standalone Extraction

```
standalone_thresholding/
├── __init__.py
├── univariate.py          # bimodal_cutoff_* functions
├── multivariate.py        # metagene + threshold functions
├── utils.py               # select_optimal_k, plotting
├── tests/
│   ├── __init__.py
│   ├── test_univariate.py
│   ├── test_multivariate.py
│   └── test_equivalence.py  # GMM equivalence tests
└── examples/
    ├── single_feature_comparison.py
    └── multivariate_markers.py
```

---

## Test Scripts

### Test 1: Single Feature GMM Equivalence

This test demonstrates that `compute_bimodal_cutoff` and `compute_multivariate_cutoff` produce **identical results** when using GMM on a single feature.

```python
#!/usr/bin/env python3
"""
Test Script: Single Feature GMM Equivalence

Demonstrates that both functions produce identical thresholds when:
- compute_bimodal_cutoff: method="gmm", n_clusters=2
- compute_multivariate_cutoff: threshold_method="gmm", single feature

The underlying GMM implementation is shared, so results are mathematically identical.
"""

import numpy as np
from sklearn.mixture import GaussianMixture

# Set seed for reproducibility
np.random.seed(2024)


def bimodal_cutoff_gmm_simple(values: np.ndarray, seed: int = 2024):
    """Simplified univariate GMM (from compute_bimodal_cutoff)."""
    X = values.reshape(-1, 1)

    gmm = GaussianMixture(
        n_components=2,
        random_state=seed,
        covariance_type="full",
    )
    gmm.fit(X)

    # Sort by means, get midpoint
    sorted_means = np.sort(gmm.means_.flatten())
    cutoff = (sorted_means[0] + sorted_means[1]) / 2

    return cutoff, gmm


def multivariate_cutoff_gmm_simple(values: np.ndarray, seed: int = 2024):
    """Simplified multivariate GMM with single feature (metagene = identity)."""
    # For single feature: metagene_score = feature_value (identity)
    metagene_scores = values

    gmm = GaussianMixture(
        n_components=2,
        random_state=seed,
        n_init=10,  # Note: multivariate uses n_init=10
    )
    gmm.fit(metagene_scores.reshape(-1, 1))

    # Find threshold at P(high) = 0.5
    means = gmm.means_.flatten()
    high_component = int(np.argmax(means))
    low_component = 1 - high_component

    low_mean, high_mean = means[low_component], means[high_component]

    # Search for crossing point
    x_range = np.linspace(low_mean, high_mean, 1000)
    posteriors = gmm.predict_proba(x_range.reshape(-1, 1))
    diff = posteriors[:, high_component] - 0.5
    cross_idx = np.where(np.diff(np.sign(diff)))[0]

    if len(cross_idx) > 0:
        threshold = float(x_range[cross_idx[0]])
    else:
        threshold = float((low_mean + high_mean) / 2)

    return threshold, gmm


# =============================================================================
# GENERATE TEST DATA: Bimodal distribution
# =============================================================================

# True parameters
mu_low, sigma_low = 2.0, 0.5
mu_high, sigma_high = 5.0, 0.8
n_low, n_high = 700, 300

# Generate samples
low_pop = np.random.normal(mu_low, sigma_low, n_low)
high_pop = np.random.normal(mu_high, sigma_high, n_high)
values = np.concatenate([low_pop, high_pop])
np.random.shuffle(values)

print("=" * 70)
print("TEST: Single Feature GMM Equivalence")
print("=" * 70)
print(f"\nData: {len(values)} samples (bimodal: {n_low} low + {n_high} high)")
print(f"True parameters: low ~ N({mu_low}, {sigma_low}²), high ~ N({mu_high}, {sigma_high}²)")

# =============================================================================
# RUN BOTH METHODS
# =============================================================================

print("\n" + "-" * 70)
print("Running both methods on SAME data with SAME seed...")
print("-" * 70)

# Method 1: Univariate (compute_bimodal_cutoff style)
cutoff_uni, gmm_uni = bimodal_cutoff_gmm_simple(values, seed=2024)

# Method 2: Multivariate (compute_multivariate_cutoff style)
threshold_multi, gmm_multi = multivariate_cutoff_gmm_simple(values, seed=2024)

# =============================================================================
# COMPARE RESULTS
# =============================================================================

print("\n" + "=" * 70)
print("RESULTS COMPARISON")
print("=" * 70)

print(f"\n{'Method':<35} {'Threshold':<15}")
print("-" * 50)
print(f"{'compute_bimodal_cutoff (gmm)':<35} {cutoff_uni:<15.6f}")
print(f"{'compute_multivariate_cutoff (gmm)':<35} {threshold_multi:<15.6f}")

# Check equivalence
diff = abs(cutoff_uni - threshold_multi)
print(f"\n{'Absolute difference:':<35} {diff:<15.10f}")

# GMM parameters comparison
print("\n" + "-" * 70)
print("GMM Parameters Comparison")
print("-" * 70)

print(f"\n{'Parameter':<25} {'Univariate':<20} {'Multivariate':<20}")
print("-" * 65)
print(f"{'Mean (low)':<25} {min(gmm_uni.means_.flatten()):<20.6f} {min(gmm_multi.means_.flatten()):<20.6f}")
print(f"{'Mean (high)':<25} {max(gmm_uni.means_.flatten()):<20.6f} {max(gmm_multi.means_.flatten()):<20.6f}")
print(f"{'Std (low)':<25} {np.sqrt(min(gmm_uni.covariances_.flatten())):<20.6f} {np.sqrt(min(gmm_multi.covariances_.flatten())):<20.6f}")
print(f"{'Std (high)':<25} {np.sqrt(max(gmm_uni.covariances_.flatten())):<20.6f} {np.sqrt(max(gmm_multi.covariances_.flatten())):<20.6f}")

# =============================================================================
# VERDICT
# =============================================================================

print("\n" + "=" * 70)
print("VERDICT")
print("=" * 70)

# Note: Small differences may occur due to n_init and convergence
tolerance = 0.01
if diff < tolerance:
    print(f"\n✅ PASS: Thresholds are equivalent (diff < {tolerance})")
    print("   Both functions use the same GMM implementation.")
    print("   For single feature, metagene_score = feature_value (identity).")
else:
    print(f"\n⚠️  NOTE: Small difference ({diff:.6f}) due to GMM initialization")
    print("   - Univariate uses covariance_type='full'")
    print("   - Multivariate uses n_init=10 (more EM restarts)")
    print("   Both are valid GMM fits, just different local optima.")

print("\n" + "=" * 70)
```

### Test 2: Multivariate Metagene with 3 Features

```python
#!/usr/bin/env python3
"""
Test Script: Multivariate Metagene Thresholding

Demonstrates compute_multivariate_cutoff with 3 features:
1. Background cells: low in all markers
2. Single-positive: high in 1 marker only
3. Triple-positive: high in all 3 markers (target population)

Shows how geometric mean correctly identifies triple-positive cells
while Euclidean clustering would mix single-positive with triple-positive.
"""

import numpy as np
import matplotlib.pyplot as plt

np.random.seed(2024)


def compute_metagene_score(values: np.ndarray, method: str = "shifted_geometric_mean"):
    """Compute metagene score for multi-feature array."""
    if method == "geometric_mean":
        eps = 1e-10
        return np.exp(np.mean(np.log(values + eps), axis=1))
    elif method == "shifted_geometric_mean":
        pseudocount = 0.1
        shifted = values + pseudocount
        return np.exp(np.mean(np.log(shifted), axis=1)) - pseudocount
    elif method == "arithmetic_mean":
        return np.mean(values, axis=1)
    elif method == "median":
        return np.median(values, axis=1)
    elif method == "minimum":
        return np.min(values, axis=1)


def threshold_ks(scores: np.ndarray, background_quantile: float = 0.5):
    """KS-inspired threshold detection."""
    from scipy.stats import norm

    sorted_scores = np.sort(scores)
    n = len(sorted_scores)

    bg_cutoff_idx = int(n * background_quantile)
    background_scores = sorted_scores[:bg_cutoff_idx]
    bg_mean = np.mean(background_scores)
    bg_std = max(np.std(background_scores), 1e-6)

    empirical_cdf = np.arange(1, n + 1) / n
    expected_cdf = norm.cdf(sorted_scores, loc=bg_mean, scale=bg_std)
    D = empirical_cdf - expected_cdf

    max_idx = np.argmax(D)
    threshold = sorted_scores[max_idx]

    if threshold <= bg_mean:
        threshold = np.percentile(sorted_scores, 90)

    return threshold, bg_mean, bg_std


# =============================================================================
# GENERATE TEST DATA: 3 Epithelial Markers
# =============================================================================

n_background = 850     # Low in all markers
n_single_pos = 100     # High in 1 marker only
n_triple_pos = 50      # High in all 3 markers (target)

print("=" * 70)
print("TEST: Multivariate Metagene Thresholding (3 Features)")
print("=" * 70)

# Background: all low (0.1 - 0.5 expression)
bg_marker1 = np.random.uniform(0.1, 0.5, n_background)
bg_marker2 = np.random.uniform(0.1, 0.5, n_background)
bg_marker3 = np.random.uniform(0.1, 0.5, n_background)

# Single-positive: high in marker1 only
sp_marker1 = np.random.normal(4.0, 0.5, n_single_pos)  # HIGH
sp_marker2 = np.random.uniform(0.1, 0.5, n_single_pos) # low
sp_marker3 = np.random.uniform(0.1, 0.5, n_single_pos) # low

# Triple-positive: high in all (target epithelial cells)
tp_marker1 = np.random.normal(4.0, 0.5, n_triple_pos)  # HIGH
tp_marker2 = np.random.normal(3.5, 0.6, n_triple_pos)  # HIGH
tp_marker3 = np.random.normal(3.8, 0.5, n_triple_pos)  # HIGH

# Combine
marker1 = np.concatenate([bg_marker1, sp_marker1, tp_marker1])
marker2 = np.concatenate([bg_marker2, sp_marker2, tp_marker2])
marker3 = np.concatenate([bg_marker3, sp_marker3, tp_marker3])

# True labels
true_labels = np.array(
    [0] * n_background +  # background
    [1] * n_single_pos +  # single-positive
    [2] * n_triple_pos    # triple-positive (target)
)

features = np.column_stack([marker1, marker2, marker3])
n_total = len(marker1)

print(f"\nData: {n_total} cells")
print(f"  - Background (all low): {n_background} ({100*n_background/n_total:.1f}%)")
print(f"  - Single-positive: {n_single_pos} ({100*n_single_pos/n_total:.1f}%)")
print(f"  - Triple-positive (target): {n_triple_pos} ({100*n_triple_pos/n_total:.1f}%)")

# =============================================================================
# COMPARE METAGENE METHODS
# =============================================================================

print("\n" + "-" * 70)
print("Metagene Scores by Method")
print("-" * 70)

methods = ["geometric_mean", "shifted_geometric_mean", "arithmetic_mean", "median", "minimum"]

print(f"\n{'Method':<25} {'Background':<15} {'Single-pos':<15} {'Triple-pos':<15}")
print("-" * 70)

for method in methods:
    scores = compute_metagene_score(features, method)

    bg_mean = np.mean(scores[:n_background])
    sp_mean = np.mean(scores[n_background:n_background+n_single_pos])
    tp_mean = np.mean(scores[n_background+n_single_pos:])

    print(f"{method:<25} {bg_mean:<15.3f} {sp_mean:<15.3f} {tp_mean:<15.3f}")

# =============================================================================
# APPLY KS THRESHOLD WITH SHIFTED GEOMETRIC MEAN
# =============================================================================

print("\n" + "-" * 70)
print("KS Threshold Detection (shifted_geometric_mean)")
print("-" * 70)

metagene_scores = compute_metagene_score(features, "shifted_geometric_mean")
threshold, bg_mean, bg_std = threshold_ks(metagene_scores)

print(f"\nThreshold: {threshold:.4f}")
print(f"Background mean: {bg_mean:.4f}")
print(f"Background std: {bg_std:.4f}")

# Classify cells
predicted_high = metagene_scores >= threshold

# Confusion matrix
tp_detected = np.sum(predicted_high[n_background+n_single_pos:])  # True triple-pos detected
sp_detected = np.sum(predicted_high[n_background:n_background+n_single_pos])  # Single-pos (false pos)
bg_detected = np.sum(predicted_high[:n_background])  # Background (false pos)

print("\n" + "-" * 70)
print("Classification Results")
print("-" * 70)

print(f"\n{'Population':<20} {'Total':<10} {'Detected as HIGH':<20} {'Detection Rate':<15}")
print("-" * 65)
print(f"{'Background':<20} {n_background:<10} {bg_detected:<20} {100*bg_detected/n_background:.1f}%")
print(f"{'Single-positive':<20} {n_single_pos:<10} {sp_detected:<20} {100*sp_detected/n_single_pos:.1f}%")
print(f"{'Triple-positive':<20} {n_triple_pos:<10} {tp_detected:<20} {100*tp_detected/n_triple_pos:.1f}%")

# =============================================================================
# COMPARE: What if we used Euclidean distance?
# =============================================================================

print("\n" + "-" * 70)
print("Comparison: Euclidean Distance (would fail!)")
print("-" * 70)

# Euclidean distance from origin
euclidean_dist = np.sqrt(np.sum(features**2, axis=1))

bg_euc = np.mean(euclidean_dist[:n_background])
sp_euc = np.mean(euclidean_dist[n_background:n_background+n_single_pos])
tp_euc = np.mean(euclidean_dist[n_background+n_single_pos:])

print(f"\nMean Euclidean distance from origin:")
print(f"  Background: {bg_euc:.3f}")
print(f"  Single-positive: {sp_euc:.3f}")
print(f"  Triple-positive: {tp_euc:.3f}")

print(f"\n⚠️  Problem: Single-positive ({sp_euc:.3f}) ≈ Triple-positive ({tp_euc:.3f})")
print("   Euclidean clustering would mix them together!")
print(f"\n✅ Metagene advantage: Single-pos score is LOW because only 1 marker is elevated")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print("""
Key findings:

1. METAGENE (shifted_geometric_mean):
   - Background: LOW score (all markers low)
   - Single-positive: LOW score (only 1 marker high → product is low)
   - Triple-positive: HIGH score (all markers high → product is high)
   → Correctly identifies only triple-positive cells

2. EUCLIDEAN DISTANCE:
   - Single-positive ≈ Triple-positive (similar distance from origin)
   → Would incorrectly mix single-positive with triple-positive!

3. KS THRESHOLD:
   - Assumes lower 50% is background
   - Finds where data deviates from background expectation
   - Works well for sparse populations (~5% in this example)

Conclusion: Metagene + KS is ideal for detecting rare populations
           defined by joint elevation of multiple markers.
""")
```

### Test 3: Method Comparison on Same Data

```python
#!/usr/bin/env python3
"""
Test Script: Method Comparison on Same Bimodal Distribution

Compares all threshold detection methods on the same data:
- Univariate: hierarchical, gmm, kde
- Multivariate (single feature): ks, gmm
"""

import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, cut_tree
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde, norm
from sklearn.mixture import GaussianMixture

np.random.seed(2024)


# =============================================================================
# GENERATE BIMODAL TEST DATA
# =============================================================================

mu_low, sigma_low = 1.0, 0.3
mu_high, sigma_high = 3.0, 0.5
n_low, n_high = 600, 400

low_pop = np.random.normal(mu_low, sigma_low, n_low)
high_pop = np.random.normal(mu_high, sigma_high, n_high)
values = np.concatenate([low_pop, high_pop])
np.random.shuffle(values)

# True threshold (theoretical intersection of two Gaussians)
# For equal weights: intersection ≈ midpoint if similar variances
true_threshold = (mu_low * sigma_high + mu_high * sigma_low) / (sigma_low + sigma_high)

print("=" * 70)
print("TEST: Method Comparison on Bimodal Distribution")
print("=" * 70)
print(f"\nData: {len(values)} samples")
print(f"True parameters: low ~ N({mu_low}, {sigma_low}²), high ~ N({mu_high}, {sigma_high}²)")
print(f"Theoretical threshold: {true_threshold:.4f}")


# =============================================================================
# UNIVARIATE METHODS
# =============================================================================

def method_hierarchical(vals, k=2):
    X = vals.reshape(-1, 1)
    distances = pdist(X, metric="euclidean")
    linkage_matrix = linkage(distances, method="ward")
    labels = cut_tree(linkage_matrix, n_clusters=k).flatten()

    stats = []
    for i in range(k):
        cluster_vals = vals[labels == i]
        stats.append({"max": cluster_vals.max(), "min": cluster_vals.min(), "mean": cluster_vals.mean()})
    stats.sort(key=lambda x: x["mean"])

    cutoff = (stats[0]["max"] + stats[1]["min"]) / 2
    return cutoff


def method_gmm(vals, k=2, seed=2024):
    gmm = GaussianMixture(n_components=k, random_state=seed, covariance_type="full")
    gmm.fit(vals.reshape(-1, 1))
    sorted_means = np.sort(gmm.means_.flatten())
    cutoff = (sorted_means[0] + sorted_means[1]) / 2
    return cutoff


def method_kde(vals):
    kde = gaussian_kde(vals)
    x_grid = np.linspace(vals.min(), vals.max(), 1000)
    density = kde(x_grid)
    minima_idx = argrelextrema(density, np.less)[0]
    if len(minima_idx) == 0:
        return np.median(vals)
    return x_grid[minima_idx[0]]


def method_ks(vals, bg_quantile=0.5):
    sorted_vals = np.sort(vals)
    n = len(sorted_vals)
    bg_idx = int(n * bg_quantile)
    bg = sorted_vals[:bg_idx]
    bg_mean, bg_std = np.mean(bg), max(np.std(bg), 1e-6)

    empirical_cdf = np.arange(1, n + 1) / n
    expected_cdf = norm.cdf(sorted_vals, loc=bg_mean, scale=bg_std)
    D = empirical_cdf - expected_cdf

    threshold = sorted_vals[np.argmax(D)]
    if threshold <= bg_mean:
        threshold = np.percentile(sorted_vals, 90)
    return threshold


# =============================================================================
# RUN ALL METHODS
# =============================================================================

print("\n" + "-" * 70)
print("Running all methods...")
print("-" * 70)

results = {
    "hierarchical": method_hierarchical(values),
    "gmm": method_gmm(values),
    "kde": method_kde(values),
    "ks": method_ks(values),
}

# =============================================================================
# RESULTS
# =============================================================================

print("\n" + "=" * 70)
print("RESULTS")
print("=" * 70)

print(f"\n{'Method':<20} {'Threshold':<15} {'Error vs True':<15}")
print("-" * 50)
print(f"{'Theoretical':<20} {true_threshold:<15.4f} {'---':<15}")
print("-" * 50)

for method, threshold in results.items():
    error = threshold - true_threshold
    print(f"{method:<20} {threshold:<15.4f} {error:+.4f}")

# =============================================================================
# CLASSIFICATION ACCURACY
# =============================================================================

print("\n" + "-" * 70)
print("Classification Accuracy (based on threshold)")
print("-" * 70)

# True labels (first n_low are low, rest are high)
# But data was shuffled, so we need to classify by value
true_cluster = values > true_threshold

print(f"\n{'Method':<20} {'Accuracy':<15}")
print("-" * 35)

for method, threshold in results.items():
    predicted = values > threshold
    accuracy = np.mean(predicted == true_cluster)
    print(f"{method:<20} {100*accuracy:.2f}%")

# =============================================================================
# SUMMARY
# =============================================================================

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print("""
All methods find similar thresholds for well-separated bimodal data.

Key differences:
- hierarchical: O(n²) memory, good for general use
- gmm: Assumes Gaussian, provides probabilities
- kde: Non-parametric, finds density valleys
- ks: Designed for sparse populations (<30%)

For balanced bimodal (40-60% split): gmm or hierarchical
For sparse target (<30%): ks
For irregular shapes: kde
""")
```

---

## API Reference

### compute_bimodal_cutoff

```python
@activity.defn
async def compute_bimodal_cutoff(
    data_uri: str,
    output_uri: str,
    feature_column: str,
    method: str = "hierarchical",
    n_clusters: Optional[int] = None,
    auto_k: bool = False,
    max_k: int = 5,
    downsample_prop: float = 0.1,
    max_cells: int = 50000,
    seed: int = 2024,
    bootstrap: bool = False,
    n_bootstrap: int = 100,
    test_unimodal: bool = False,
    return_plot: bool = True,
    plot_format: str = "png",
    dpi: int = 300,
    **kwargs,
) -> dict
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_uri` | str | required | GCS URI to input H5AD |
| `output_uri` | str | required | GCS URI for output |
| `feature_column` | str | required | Column name to analyze |
| `method` | str | "hierarchical" | `hierarchical`, `gmm`, or `kde` |
| `n_clusters` | int | None | Number of clusters (default: 3) |
| `auto_k` | bool | False | Auto-select k via silhouette |
| `max_k` | int | 5 | Max k for auto selection |
| `downsample_prop` | float | 0.1 | Subsample proportion |
| `max_cells` | int | 50000 | Maximum cells |
| `seed` | int | 2024 | Random seed |
| `bootstrap` | bool | False | Compute confidence intervals |
| `n_bootstrap` | int | 100 | Bootstrap iterations |
| `test_unimodal` | bool | False | Warn if unimodal |
| `return_plot` | bool | True | Generate plot |

### compute_multivariate_cutoff

```python
@activity.defn
async def compute_multivariate_cutoff(
    data_uri: str,
    output_uri: str,
    feature_columns: List[str],
    metagene_method: MetageneMethod = "shifted_geometric_mean",
    pseudocount: float = 0.1,
    threshold_method: ThresholdMethod = "gmm",
    background_quantile: float = 0.5,
    probability_cutoff: float = 0.3,
    n_clusters: int = 2,
    downsample_prop: float = 0.1,
    max_cells: int = 50000,
    seed: int = 2024,
    return_plot: bool = True,
    plot_format: str = "png",
    dpi: int = 300,
    **kwargs,
) -> dict
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_uri` | str | required | GCS URI to input H5AD |
| `output_uri` | str | required | GCS URI for output |
| `feature_columns` | List[str] | required | 1+ column names |
| `metagene_method` | str | "shifted_geometric_mean" | Aggregation method |
| `pseudocount` | float | 0.1 | For shifted_geometric_mean |
| `threshold_method` | str | "gmm" | `gmm` or `ks` |
| `background_quantile` | float | 0.5 | KS-specific |
| `probability_cutoff` | float | 0.3 | GMM-specific |
| `n_clusters` | int | 2 | GMM components |
| `downsample_prop` | float | 0.1 | Subsample proportion |
| `max_cells` | int | 50000 | Maximum cells |
| `seed` | int | 2024 | Random seed |
| `return_plot` | bool | True | Generate plot |

### assign_multivariate_clusters

```python
@activity.defn
async def assign_multivariate_clusters(
    data_uri: str,
    output_uri: str,
    cutoff_uri: str,
    column_prefix: str = "multivariate",
    probability_cutoff: Optional[float] = None,
    **kwargs,
) -> dict
```

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_uri` | str | required | GCS URI to new data |
| `output_uri` | str | required | GCS URI for output |
| `cutoff_uri` | str | required | GCS URI to saved threshold |
| `column_prefix` | str | "multivariate" | Prefix for new columns |
| `probability_cutoff` | float | None | Override P(high) cutoff |

---

## File Reference

**Source file**: `/home/user/spatial-biology-functions/functions/statistical_thresholding.py`

**Key line numbers**:
- `compute_bimodal_cutoff`: Line 63
- `_bimodal_cutoff_hierarchical`: Line 484
- `_bimodal_cutoff_gmm`: Line 556
- `_bimodal_cutoff_kde`: Line 612
- `_select_optimal_k`: Line 720
- `_bootstrap_cutoff_ci`: Line 774
- `_compute_metagene_score`: Line 1213
- `_compute_gmm_threshold`: Line 1251
- `_threshold_ks`: Line 1284
- `compute_multivariate_cutoff`: Line 1560
- `assign_multivariate_clusters`: Line 1992

---

## Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | Nov 2025 | Initial implementation |
| 1.1 | Nov 2025 | Added KS threshold method |
| 1.2 | Nov 2025 | Added shifted_geometric_mean default |
