# Spatial Autocorrelation Guide: Moran's I and Lee's L

This guide provides comprehensive documentation for spatial autocorrelation analysis using Moran's I and Lee's L statistics. It includes standalone implementations that can be used without the `process_activity` framework or Temporal orchestration.

## Table of Contents

1. [Overview](#1-overview)
2. [Moran's I](#2-morans-i)
3. [Lee's L](#3-lees-l)
4. [Visualization Functions](#4-visualization-functions)
5. [Standalone Scripts](#5-standalone-scripts)
6. [Test Examples](#6-test-examples)

---

## 1. Overview

### What is Spatial Autocorrelation?

Spatial autocorrelation measures whether values of a variable are clustered, dispersed, or random in space. In spatial transcriptomics, this helps identify:
- Genes with spatially organized expression patterns
- Tissue regions where genes co-localize
- Cell-cell interactions and microenvironments

### When to Use Each Statistic

| Statistic | Use Case | Output |
|-----------|----------|--------|
| **Moran's I** | Univariate - Is gene X spatially clustered? | Single I value per gene |
| **Lee's L** | Bivariate - Are genes X and Y spatially co-localized? | Single L value per gene pair |
| **Local Moran's I** | Where is gene X clustered in the tissue? | Per-cell I values |
| **Local Lee's L** | Where do genes X and Y co-localize? | Per-cell L values + quadrants |

### Global vs Local Statistics

- **Global**: One summary statistic for the entire tissue (e.g., "CD38 has I=0.45")
- **Local**: Per-cell values showing which regions contribute to the pattern (e.g., "These 5,000 cells show CD38 clustering")

---

## 2. Moran's I

### 2.1 Mathematical Background

**Global Moran's I Formula:**

```
I = (N / W) × Σᵢⱼ wᵢⱼ(xᵢ - x̄)(xⱼ - x̄) / Σᵢ(xᵢ - x̄)²
```

Where:
- `N` = number of spatial units (cells)
- `W` = sum of all spatial weights
- `wᵢⱼ` = spatial weight between cells i and j
- `xᵢ` = value at cell i
- `x̄` = mean value

**Local Moran's I Formula:**

```
Iᵢ = zᵢ × Σⱼ wᵢⱼ × zⱼ
```

Where:
- `zᵢ = (xᵢ - x̄) / σ` (standardized value)
- The sum is over neighbors j of cell i

**Interpretation:**
- `I > 0`: Positive autocorrelation (similar values cluster together)
- `I ≈ 0`: Random spatial distribution
- `I < 0`: Negative autocorrelation (checkerboard pattern, dissimilar neighbors)

### 2.2 Standalone Implementation

```python
import numpy as np
import scanpy as sc
import squidpy as sq
from scipy import sparse


def compute_spatial_neighbors(adata, n_neighbors=6):
    """
    Compute spatial neighbor graph using squidpy.

    Args:
        adata: AnnData with spatial coordinates in adata.obsm['spatial']
        n_neighbors: Number of nearest neighbors (default: 6)

    Returns:
        adata with spatial_connectivities in adata.obsp
    """
    sq.gr.spatial_neighbors(adata, n_neighs=n_neighbors, coord_type="generic")
    return adata


def compute_global_morans_i(adata, genes=None, permutations=200):
    """
    Compute global Moran's I for specified genes.

    Args:
        adata: AnnData with spatial_connectivities computed
        genes: List of genes (None = use highly variable or top 50)
        permutations: Number of permutations for p-value (default: 200)

    Returns:
        DataFrame with columns: I, pval_sim, pval_norm, var_norm
        Stored in adata.uns['moranI']
    """
    if genes is None:
        if 'highly_variable' in adata.var.columns:
            genes = adata.var_names[adata.var['highly_variable']].tolist()[:50]
        else:
            # Use top 50 genes by mean expression
            gene_means = np.array(adata.X.mean(axis=0)).flatten()
            top_idx = np.argsort(gene_means)[-50:]
            genes = adata.var_names[top_idx].tolist()

    # Subset to genes of interest
    adata_sub = adata[:, genes].copy()

    # Compute using squidpy
    sq.gr.spatial_autocorr(
        adata_sub,
        mode="moran",
        n_perms=permutations,
        n_jobs=1
    )

    # Copy results back
    adata.uns['moranI'] = adata_sub.uns['moranI']

    return adata.uns['moranI']


def compute_local_morans_i(adata, gene):
    """
    Compute local Moran's I for a single gene.

    Args:
        adata: AnnData with spatial_connectivities computed
        gene: Gene name to analyze

    Returns:
        Array of local Moran's I values (one per cell)
        Stored in adata.obs['{gene}_local_morans_i']
    """
    # Get spatial weights matrix
    W = adata.obsp['spatial_connectivities']

    # Row-normalize weights
    if sparse.issparse(W):
        row_sums = np.array(W.sum(axis=1)).flatten()
        row_sums[row_sums == 0] = 1
        W_norm = W.multiply(1.0 / row_sums[:, np.newaxis])
    else:
        row_sums = W.sum(axis=1)
        row_sums[row_sums == 0] = 1
        W_norm = W / row_sums[:, np.newaxis]

    # Get gene expression
    gene_idx = adata.var_names.get_loc(gene)
    if hasattr(adata.X, 'toarray'):
        x = adata.X[:, gene_idx].toarray().flatten()
    else:
        x = adata.X[:, gene_idx].flatten()

    # Standardize
    x_mean, x_std = x.mean(), x.std()
    if x_std == 0:
        return np.zeros(len(x))
    z = (x - x_mean) / x_std

    # Compute spatial lag
    if sparse.issparse(W_norm):
        lag_z = np.array(W_norm @ z).flatten()
    else:
        lag_z = W_norm @ z

    # Local Moran's I
    local_I = z * lag_z

    # Store in adata
    adata.obs[f'{gene}_local_morans_i'] = local_I

    return local_I
```

### 2.3 Usage Examples

**Discovery Workflow (Find Spatially Variable Genes):**

```python
import scanpy as sc

# Load data
adata = sc.read_h5ad("spatial_data.h5ad")

# Compute spatial neighbors
compute_spatial_neighbors(adata, n_neighbors=6)

# Compute global Moran's I for top genes
morans_df = compute_global_morans_i(adata, permutations=200)

# Find significant spatially variable genes
significant = morans_df[morans_df['pval_sim'] < 0.05].sort_values('I', ascending=False)
print(f"Found {len(significant)} spatially variable genes")
print(significant.head(10))
```

**Cell-Level Analysis (Local Moran's I):**

```python
# Compute local scores for a specific gene
gene = "MS4A1"
local_I = compute_local_morans_i(adata, gene)

# Identify hotspots (high local I with high expression)
hotspot_mask = (local_I > 0.5) & (adata[:, gene].X.toarray().flatten() > 1)
print(f"Found {hotspot_mask.sum()} hotspot cells for {gene}")
```

---

## 3. Lee's L

### 3.1 Mathematical Background

**Global Lee's L Formula:**

```
L = N × Σᵢⱼ wᵢⱼ(xᵢ - x̄)(yⱼ - ȳ) / √[Σᵢ(xᵢ - x̄)² × Σⱼ(yⱼ - ȳ)²]
```

Where:
- `xᵢ`, `yⱼ` = values of two different variables
- This measures bivariate spatial correlation

**Local Lee's L Formula:**

```
Lᵢ = zₓᵢ × lag(zᵧ)ᵢ
```

Where:
- `zₓᵢ = (xᵢ - x̄) / σₓ` (standardized gene1 at cell i)
- `lag(zᵧ)ᵢ = Σⱼ wᵢⱼ × zᵧⱼ` (spatially-lagged standardized gene2)

**Interpretation:**
- `L > 0`: Positive spatial correlation (genes co-localize)
- `L ≈ 0`: No spatial relationship
- `L < 0`: Negative spatial correlation (genes anti-correlate spatially)

### 3.2 Quadrant Classification

The quadrant classification uses standardized values to distinguish four spatial co-occurrence patterns:

| Quadrant | Condition | Interpretation | Local L Sign |
|----------|-----------|----------------|--------------|
| **HH** | `zₓ > 0` AND `lag(zᵧ) > 0` | High gene1, neighbors have high gene2 | Positive |
| **LL** | `zₓ < 0` AND `lag(zᵧ) < 0` | Low gene1, neighbors have low gene2 | Positive |
| **HL** | `zₓ > 0` AND `lag(zᵧ) < 0` | High gene1, neighbors have low gene2 | Negative |
| **LH** | `zₓ < 0` AND `lag(zᵧ) > 0` | Low gene1, neighbors have high gene2 | Negative |

**Key Insight:** Both HH and LL give positive local L values, but represent fundamentally different biological patterns:
- **HH**: Active co-localization (e.g., disease regions)
- **LL**: Co-absence (e.g., healthy baseline tissue)

### 3.3 Standalone Implementation

```python
import numpy as np
import scanpy as sc
import squidpy as sq
from scipy import sparse


def compute_local_lees_l(adata, gene1, gene2, n_neighbors=6):
    """
    Compute local Lee's L bivariate spatial autocorrelation with quadrant classification.

    Args:
        adata: AnnData with spatial coordinates in adata.obsm['spatial']
        gene1: First gene name (e.g., "CD38")
        gene2: Second gene name (e.g., "COL1A1")
        n_neighbors: Number of spatial neighbors (default: 6)

    Returns:
        dict with:
            - local_L: Array of per-cell local Lee's L values
            - quadrant: Array of quadrant classifications (HH, LL, HL, LH)
            - global_L: Sum of local L values
            - mean_L: Mean of local L values
            - n_HH, n_LL, n_HL, n_LH: Quadrant counts

    Side Effects:
        Adds columns to adata.obs:
            - {gene1}_{gene2}_local_lees_l
            - {gene1}_{gene2}_quadrant
            - {gene1}_{gene2}_z_gene1
            - {gene1}_{gene2}_lag_z_gene2

        Adds summary to adata.uns['local_lees_l']
    """
    # Validate genes exist
    for gene in [gene1, gene2]:
        if gene not in adata.var_names:
            raise ValueError(f"Gene '{gene}' not found in adata.var_names")

    # Compute spatial neighbors if needed
    if 'spatial_connectivities' not in adata.obsp:
        print(f"Computing spatial neighbors (k={n_neighbors})...")
        sq.gr.spatial_neighbors(adata, n_neighs=n_neighbors, coord_type="generic")

    # Get spatial weights matrix
    W = adata.obsp['spatial_connectivities']

    # Row-normalize weights (each row sums to 1)
    if sparse.issparse(W):
        row_sums = np.array(W.sum(axis=1)).flatten()
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        W_normalized = W.multiply(1.0 / row_sums[:, np.newaxis])
    else:
        row_sums = W.sum(axis=1)
        row_sums[row_sums == 0] = 1
        W_normalized = W / row_sums[:, np.newaxis]

    # Extract gene expression
    gene1_idx = adata.var_names.get_loc(gene1)
    gene2_idx = adata.var_names.get_loc(gene2)

    if hasattr(adata.X, 'toarray'):
        x = adata.X[:, gene1_idx].toarray().flatten()
        y = adata.X[:, gene2_idx].toarray().flatten()
    else:
        x = adata.X[:, gene1_idx].flatten()
        y = adata.X[:, gene2_idx].flatten()

    # Compute local Lee's L
    x_std, y_std = x.std(), y.std()

    if x_std == 0 or y_std == 0:
        print("WARNING: Zero variance in gene expression")
        n = len(x)
        local_L = np.zeros(n)
        quadrant = np.array(["NS"] * n)
        z_x = np.zeros(n)
        lag_z_y = np.zeros(n)
    else:
        # Standardize
        z_x = (x - x.mean()) / x_std
        z_y = (y - y.mean()) / y_std

        # Spatial lag of gene2
        if sparse.issparse(W_normalized):
            lag_z_y = np.array(W_normalized @ z_y).flatten()
        else:
            lag_z_y = W_normalized @ z_y

        # Local Lee's L: product of standardized gene1 and lagged gene2
        local_L = z_x * lag_z_y

        # Quadrant classification
        quadrant = np.empty(len(x), dtype='U2')
        quadrant[:] = "NS"  # Default (not significant / zero)
        quadrant[(z_x > 0) & (lag_z_y > 0)] = "HH"
        quadrant[(z_x < 0) & (lag_z_y < 0)] = "LL"
        quadrant[(z_x > 0) & (lag_z_y < 0)] = "HL"
        quadrant[(z_x < 0) & (lag_z_y > 0)] = "LH"

    # Store results in adata
    col_prefix = f"{gene1}_{gene2}"
    adata.obs[f"{col_prefix}_local_lees_l"] = local_L
    adata.obs[f"{col_prefix}_quadrant"] = quadrant
    adata.obs[f"{col_prefix}_z_gene1"] = z_x
    adata.obs[f"{col_prefix}_lag_z_gene2"] = lag_z_y

    # Summary statistics
    global_L = local_L.sum()
    n_HH = (quadrant == "HH").sum()
    n_LL = (quadrant == "LL").sum()
    n_HL = (quadrant == "HL").sum()
    n_LH = (quadrant == "LH").sum()

    summary = {
        "gene1": gene1,
        "gene2": gene2,
        "global_L": float(global_L),
        "mean_local_L": float(local_L.mean()),
        "n_HH": int(n_HH),
        "n_LL": int(n_LL),
        "n_HL": int(n_HL),
        "n_LH": int(n_LH),
        "n_cells": int(len(x)),
    }

    adata.uns['local_lees_l'] = summary

    print(f"\nLocal Lee's L Results for {gene1} vs {gene2}:")
    print(f"  Global L (sum): {global_L:.4f}")
    print(f"  Mean local L: {local_L.mean():.4f}")
    print(f"  Quadrant counts:")
    print(f"    HH (high-high): {n_HH} ({100*n_HH/len(x):.1f}%)")
    print(f"    LL (low-low): {n_LL} ({100*n_LL/len(x):.1f}%)")
    print(f"    HL (high-low): {n_HL} ({100*n_HL/len(x):.1f}%)")
    print(f"    LH (low-high): {n_LH} ({100*n_LH/len(x):.1f}%)")

    return {
        "local_L": local_L,
        "quadrant": quadrant,
        **summary
    }


def compute_hh_proportion_by_celltype(adata, gene_pairs, cell_type_column="cell_type", top_n=10):
    """
    Compute HH proportion for each gene pair by cell type.

    Args:
        adata: AnnData with Lee's L quadrant columns
        gene_pairs: List of (gene1, gene2) tuples
        cell_type_column: Column with cell type annotations
        top_n: Number of top cell types to include

    Returns:
        DataFrame with columns: gene_pair, cell_type, n_cells, n_HH, HH_proportion
    """
    import pandas as pd

    # Get top cell types
    cell_type_counts = adata.obs[cell_type_column].value_counts()
    cell_type_counts = cell_type_counts[
        ~cell_type_counts.index.str.lower().str.contains('unassigned|unknown')
    ]
    top_cell_types = cell_type_counts.head(top_n).index.tolist()

    results = []
    for gene1, gene2 in gene_pairs:
        quadrant_col = f"{gene1}_{gene2}_quadrant"

        if quadrant_col not in adata.obs.columns:
            print(f"Warning: {quadrant_col} not found, skipping")
            continue

        pair_name = f"{gene1}:{gene2}"

        for cell_type in top_cell_types:
            mask = adata.obs[cell_type_column] == cell_type
            n_cells = mask.sum()

            if n_cells > 0:
                n_hh = ((adata.obs[quadrant_col] == "HH") & mask).sum()
                hh_prop = n_hh / n_cells
            else:
                n_hh = 0
                hh_prop = 0

            results.append({
                'gene_pair': pair_name,
                'cell_type': cell_type,
                'n_cells': int(n_cells),
                'n_HH': int(n_hh),
                'HH_proportion': float(hh_prop),
            })

    return pd.DataFrame(results)
```

### 3.4 Usage Examples

**Single Gene Pair Analysis:**

```python
import scanpy as sc

# Load data
adata = sc.read_h5ad("spatial_data.h5ad")

# Compute local Lee's L for CD38 vs CDKN1A
result = compute_local_lees_l(adata, gene1="CD38", gene2="CDKN1A")

# Find HH cells (co-localized high expression)
hh_mask = adata.obs["CD38_CDKN1A_quadrant"] == "HH"
print(f"Found {hh_mask.sum()} cells with CD38-CDKN1A co-localization")

# Examine cell types in HH quadrant
print(adata.obs.loc[hh_mask, "cell_type"].value_counts())
```

**Multiple Gene Pairs Analysis:**

```python
# Define gene pairs to analyze
gene_pairs = [
    ("CD38", "CDKN1A"),  # CD38 + p21 (senescence)
    ("CD38", "CCL3"),    # CD38 + MIP-1α (SASP)
    ("CD38", "CCL2"),    # CD38 + MCP-1 (SASP)
    ("CD38", "IL6"),     # CD38 + IL-6 (SASP)
    ("CD38", "VWF"),     # Negative control (endothelial)
]

# Compute Lee's L for each pair
for gene1, gene2 in gene_pairs:
    compute_local_lees_l(adata, gene1, gene2)

# Compute HH proportion by cell type
hh_df = compute_hh_proportion_by_celltype(adata, gene_pairs, cell_type_column="cell_type")
print(hh_df.pivot(index='gene_pair', columns='cell_type', values='HH_proportion'))
```

---

## 4. Visualization Functions

### 4.1 Moran's I Visualizations

#### Plot Global Moran's I (4-Panel)

```python
import matplotlib.pyplot as plt
import numpy as np


def plot_morans_i_global(adata, output_path="morans_i_global.png", dpi=300):
    """
    Create 4-panel global Moran's I visualization.

    Panels:
        1. Histogram of Moran's I distribution
        2. Moran's I vs mean expression (scatter)
        3. Spatial plot of top gene
        4. Summary statistics table
    """
    morans_df = adata.uns.get('moranI')
    if morans_df is None:
        raise ValueError("No Moran's I results found. Run compute_global_morans_i first.")

    fig = plt.figure(figsize=(16, 12))

    # Panel 1: Histogram
    ax1 = fig.add_subplot(2, 2, 1)
    ax1.hist(morans_df['I'], bins=30, edgecolor='black', alpha=0.7)
    ax1.axvline(0, color='red', linestyle='--', label='Random (I=0)')
    ax1.set_xlabel("Moran's I")
    ax1.set_ylabel("Frequency")
    ax1.set_title("Distribution of Moran's I Values")
    ax1.legend()

    # Panel 2: I vs Expression
    ax2 = fig.add_subplot(2, 2, 2)
    gene_means = np.array(adata[:, morans_df.index].X.mean(axis=0)).flatten()
    significant = morans_df['pval_sim'] < 0.05
    ax2.scatter(gene_means[~significant], morans_df['I'][~significant],
                alpha=0.3, c='gray', label='Not significant')
    ax2.scatter(gene_means[significant], morans_df['I'][significant],
                alpha=0.7, c='red', label='p < 0.05')
    ax2.set_xlabel("Mean Expression")
    ax2.set_ylabel("Moran's I")
    ax2.set_title("Moran's I vs Expression Level")
    ax2.legend()

    # Panel 3: Top gene spatial plot
    ax3 = fig.add_subplot(2, 2, 3)
    top_gene = morans_df['I'].idxmax()
    coords = adata.obsm['spatial']
    expr = adata[:, top_gene].X.toarray().flatten() if hasattr(adata.X, 'toarray') else adata[:, top_gene].X.flatten()
    scatter = ax3.scatter(coords[:, 0], coords[:, 1], c=expr, s=1, cmap='viridis')
    ax3.set_title(f"Top Gene: {top_gene} (I={morans_df.loc[top_gene, 'I']:.3f})")
    ax3.set_aspect('equal')
    ax3.invert_yaxis()
    plt.colorbar(scatter, ax=ax3, label='Expression')

    # Panel 4: Summary stats
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.axis('off')
    n_sig = (morans_df['pval_sim'] < 0.05).sum()
    stats_text = f"""
    Summary Statistics
    ──────────────────
    Genes tested: {len(morans_df)}
    Significant (p<0.05): {n_sig} ({100*n_sig/len(morans_df):.1f}%)

    Median Moran's I: {morans_df['I'].median():.4f}
    Max Moran's I: {morans_df['I'].max():.4f}

    Top 5 Genes:
    """
    for gene in morans_df.nlargest(5, 'I').index:
        stats_text += f"\n    {gene}: I={morans_df.loc[gene, 'I']:.3f}"

    ax4.text(0.1, 0.9, stats_text, transform=ax4.transAxes, fontsize=11,
             verticalalignment='top', fontfamily='monospace')

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved global Moran's I plot to {output_path}")


def plot_morans_i_local(adata, gene, output_path="morans_i_local.png", dpi=300):
    """
    Create 2-panel local Moran's I visualization.

    Panels:
        1. Gene expression in spatial coordinates
        2. Local Moran's I scores
    """
    local_col = f"{gene}_local_morans_i"
    if local_col not in adata.obs.columns:
        raise ValueError(f"Local Moran's I not computed for {gene}. Run compute_local_morans_i first.")

    coords = adata.obsm['spatial']
    expr = adata[:, gene].X.toarray().flatten() if hasattr(adata.X, 'toarray') else adata[:, gene].X.flatten()
    local_I = adata.obs[local_col].values

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel 1: Expression
    sc1 = axes[0].scatter(coords[:, 0], coords[:, 1], c=expr, s=1, cmap='viridis')
    axes[0].set_title(f"{gene} Expression")
    axes[0].set_aspect('equal')
    axes[0].invert_yaxis()
    plt.colorbar(sc1, ax=axes[0], label='Expression')

    # Panel 2: Local I (symmetric colormap)
    vmax = np.abs(local_I).max()
    sc2 = axes[1].scatter(coords[:, 0], coords[:, 1], c=local_I, s=1,
                          cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    axes[1].set_title(f"Local Moran's I")
    axes[1].set_aspect('equal')
    axes[1].invert_yaxis()
    plt.colorbar(sc2, ax=axes[1], label="Local I")

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved local Moran's I plot to {output_path}")
```

### 4.2 Lee's L Visualizations

#### Plot Local Lee's L with Quadrant Coloring

```python
def plot_local_lees_l(adata, gene1, gene2, output_path="lees_l_local.png",
                      dpi=300, point_size=1.0):
    """
    Create 4-panel local Lee's L visualization with quadrant coloring.

    Panels:
        1. Gene1 expression
        2. Gene2 expression
        3. Quadrant-colored spatial plot
        4. Summary statistics
    """
    col_prefix = f"{gene1}_{gene2}"
    local_l_col = f"{col_prefix}_local_lees_l"
    quadrant_col = f"{col_prefix}_quadrant"

    if local_l_col not in adata.obs.columns:
        raise ValueError(f"Lee's L not computed for {gene1} vs {gene2}")

    coords = adata.obsm['spatial']
    local_L = adata.obs[local_l_col].values
    quadrant = adata.obs[quadrant_col].values

    # Get expression
    expr1 = adata[:, gene1].X.toarray().flatten() if hasattr(adata.X, 'toarray') else adata[:, gene1].X.flatten()
    expr2 = adata[:, gene2].X.toarray().flatten() if hasattr(adata.X, 'toarray') else adata[:, gene2].X.flatten()

    fig = plt.figure(figsize=(16, 14))

    # Panel 1: Gene1 expression
    ax1 = fig.add_subplot(2, 2, 1)
    sc1 = ax1.scatter(coords[:, 0], coords[:, 1], c=expr1, s=point_size, cmap='viridis')
    ax1.set_title(f"{gene1} Expression", fontweight='bold')
    ax1.set_aspect('equal')
    ax1.invert_yaxis()
    plt.colorbar(sc1, ax=ax1, label='Expression')

    # Panel 2: Gene2 expression
    ax2 = fig.add_subplot(2, 2, 2)
    sc2 = ax2.scatter(coords[:, 0], coords[:, 1], c=expr2, s=point_size, cmap='magma')
    ax2.set_title(f"{gene2} Expression", fontweight='bold')
    ax2.set_aspect('equal')
    ax2.invert_yaxis()
    plt.colorbar(sc2, ax=ax2, label='Expression')

    # Panel 3: Quadrant coloring
    ax3 = fig.add_subplot(2, 2, 3)

    # Define colors for each quadrant
    colors = np.zeros((len(local_L), 4))  # RGBA
    max_L = max(np.abs(local_L).max(), 0.01)

    for i, (q, L_val) in enumerate(zip(quadrant, local_L)):
        intensity = min(abs(L_val) / max_L, 1.0)

        if q == "HH":
            colors[i] = [0.8 + 0.2*intensity, 0.2*(1-intensity), 0.2*(1-intensity), 0.8]
        elif q == "LL":
            colors[i] = [0.2*(1-intensity), 0.2*(1-intensity), 0.8 + 0.2*intensity, 0.8]
        elif q == "HL":
            colors[i] = [1.0, 0.5, 0.0, 0.5]  # Orange
        elif q == "LH":
            colors[i] = [0.5, 0.0, 1.0, 0.5]  # Purple
        else:
            colors[i] = [0.8, 0.8, 0.8, 0.3]  # Gray

    ax3.scatter(coords[:, 0], coords[:, 1], c=colors, s=point_size)
    ax3.set_title("Local Lee's L Quadrants", fontweight='bold')
    ax3.set_aspect('equal')
    ax3.invert_yaxis()

    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='red', label='HH (co-localization)'),
        Patch(facecolor='blue', label='LL (co-absence)'),
        Patch(facecolor='orange', label='HL'),
        Patch(facecolor='purple', label='LH'),
    ]
    ax3.legend(handles=legend_elements, loc='upper right')

    # Panel 4: Summary statistics
    ax4 = fig.add_subplot(2, 2, 4)
    ax4.axis('off')

    n_HH = (quadrant == "HH").sum()
    n_LL = (quadrant == "LL").sum()
    n_HL = (quadrant == "HL").sum()
    n_LH = (quadrant == "LH").sum()
    n = len(quadrant)

    stats_text = f"""
    Local Lee's L Summary
    ─────────────────────
    Gene Pair: {gene1} vs {gene2}
    Total Cells: {n:,}

    Global L (sum): {local_L.sum():.4f}
    Mean Local L: {local_L.mean():.6f}
    Std Local L: {local_L.std():.6f}

    Quadrant Distribution:
        HH (both high): {n_HH:,} ({100*n_HH/n:.1f}%)
        LL (both low): {n_LL:,} ({100*n_LL/n:.1f}%)
        HL: {n_HL:,} ({100*n_HL/n:.1f}%)
        LH: {n_LH:,} ({100*n_LH/n:.1f}%)

    Interpretation:
        HH = Active co-localization
        LL = Both genes absent
    """

    ax4.text(0.1, 0.95, stats_text, transform=ax4.transAxes, fontsize=11,
             verticalalignment='top', fontfamily='monospace')

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved local Lee's L plot to {output_path}")
```

#### Plot HH Proportion by Cell Type

```python
def plot_lees_l_hh_by_celltype(adata, gene_pairs, cell_type_column="cell_type",
                                top_n_cell_types=10, output_path="lees_l_hh_by_celltype.png",
                                dpi=300):
    """
    Create grouped bar chart of HH proportion by cell type.
    """
    import pandas as pd

    # Compute HH proportions
    hh_df = compute_hh_proportion_by_celltype(adata, gene_pairs, cell_type_column, top_n_cell_types)

    # Pivot for plotting
    pivot_df = hh_df.pivot(index='gene_pair', columns='cell_type', values='HH_proportion')

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 8))

    x = np.arange(len(pivot_df.index))
    width = 0.8 / len(pivot_df.columns)
    n_cell_types = len(pivot_df.columns)

    colors = plt.cm.tab10(np.linspace(0, 1, min(n_cell_types, 10)))

    for i, (cell_type, color) in enumerate(zip(pivot_df.columns, colors)):
        offset = (i - n_cell_types/2 + 0.5) * width
        ax.bar(x + offset, pivot_df[cell_type].fillna(0), width,
               label=cell_type, color=color, edgecolor='white', linewidth=0.5)

    ax.set_xlabel('Gene Pair', fontsize=12, fontweight='bold')
    ax.set_ylabel('Proportion of Cells in HH Quadrant', fontsize=12, fontweight='bold')
    ax.set_title('Bivariate Spatial Co-localization by Cell Type\n(Lee\'s L High-High Quadrant)',
                 fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(pivot_df.index, rotation=45, ha='right')
    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9, title='Cell Type')
    ax.axhline(y=0.25, color='gray', linestyle='--', alpha=0.5, label='Random (25%)')
    ax.set_ylim(0, None)
    ax.grid(axis='y', alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close()
    print(f"Saved HH by cell type plot to {output_path}")


def plot_lees_l_heatmap(adata, gene_pairs, cell_type_column="cell_type",
                        top_n_cell_types=10, output_path="lees_l_heatmap.png",
                        dpi=300, annotate=True):
    """
    Create heatmap of HH proportions (cell types × gene pairs).
    """
    import pandas as pd

    # Compute HH proportions
    hh_df = compute_hh_proportion_by_celltype(adata, gene_pairs, cell_type_column, top_n_cell_types)

    # Pivot for heatmap
    pivot_df = hh_df.pivot(index='cell_type', columns='gene_pair', values='HH_proportion')

    fig, ax = plt.subplots(figsize=(12, 8))

    im = ax.imshow(pivot_df.values, cmap='YlOrRd', aspect='auto')

    cbar = plt.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label('HH Proportion', fontsize=10)

    ax.set_xticks(np.arange(len(pivot_df.columns)))
    ax.set_yticks(np.arange(len(pivot_df.index)))
    ax.set_xticklabels(pivot_df.columns, rotation=45, ha='right')
    ax.set_yticklabels(pivot_df.index)
    ax.set_title('Bivariate Spatial Co-localization Heatmap\n(Lee\'s L HH Proportion)',
                 fontsize=12, fontweight='bold')

    if annotate:
        for i in range(len(pivot_df.index)):
            for j in range(len(pivot_df.columns)):
                val = pivot_df.values[i, j]
                if not np.isnan(val):
                    text_color = 'white' if val > 0.15 else 'black'
                    ax.text(j, i, f'{val:.2f}', ha='center', va='center',
                            color=text_color, fontsize=8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved heatmap to {output_path}")


def plot_lees_l_summary_table(adata, gene_pairs, output_path="lees_l_summary_table.png",
                               dpi=300, title="Lee's L Spatial Co-localization Summary"):
    """
    Create formatted summary table of Lee's L statistics.
    """
    import re

    # Auto-detect gene pairs if not provided
    if gene_pairs is None:
        pattern = re.compile(r"^(.+?)_(.+?)_quadrant$")
        gene_pairs = []
        for col in adata.obs.columns:
            match = pattern.match(col)
            if match:
                gene_pairs.append((match.group(1), match.group(2)))

    # Collect statistics
    results = []
    for gene1, gene2 in gene_pairs:
        col_prefix = f"{gene1}_{gene2}"
        quadrant_col = f"{col_prefix}_quadrant"
        local_l_col = f"{col_prefix}_local_lees_l"

        if quadrant_col not in adata.obs.columns:
            continue

        quadrants = adata.obs[quadrant_col]
        n_HH = (quadrants == "HH").sum()
        n_LL = (quadrants == "LL").sum()
        n_HL = (quadrants == "HL").sum()
        n_LH = (quadrants == "LH").sum()
        n_total = len(quadrants)

        if local_l_col in adata.obs.columns:
            local_L = adata.obs[local_l_col].values
            global_L = float(local_L.sum())
            mean_L = float(local_L.mean())
        else:
            global_L = 0.0
            mean_L = 0.0

        results.append({
            "Gene Pair": f"{gene1}:{gene2}",
            "Global L": global_L,
            "Mean L": mean_L,
            "HH": n_HH,
            "LL": n_LL,
            "HL": n_HL,
            "LH": n_LH,
            "HH %": 100 * n_HH / n_total if n_total > 0 else 0,
        })

    # Sort by Global L descending
    results = sorted(results, key=lambda x: x["Global L"], reverse=True)

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 2 + len(results) * 0.5))
    ax.axis('off')

    columns = ["Gene Pair", "Global L", "Mean L", "HH", "LL", "HL", "LH", "HH %"]
    cell_text = []
    for r in results:
        cell_text.append([
            r["Gene Pair"],
            f"{r['Global L']:,.2f}",
            f"{r['Mean L']:.6f}",
            f"{r['HH']:,}",
            f"{r['LL']:,}",
            f"{r['HL']:,}",
            f"{r['LH']:,}",
            f"{r['HH %']:.1f}%",
        ])

    table = ax.table(
        cellText=cell_text,
        colLabels=columns,
        cellLoc='center',
        loc='center',
        colColours=['#4472C4'] * len(columns),
    )

    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.8)

    # Style header
    for j in range(len(columns)):
        cell = table[(0, j)]
        cell.set_text_props(weight='bold', color='white')
        cell.set_facecolor('#4472C4')

    # Color code rows
    for i, r in enumerate(results):
        global_L = r["Global L"]
        if global_L > 0:
            intensity = min(abs(global_L) / 5000, 0.3)
            color = (0.9 - intensity, 1.0, 0.9 - intensity)
        elif global_L < 0:
            intensity = min(abs(global_L) / 5000, 0.3)
            color = (1.0, 0.9 - intensity, 0.9 - intensity)
        else:
            color = 'white'

        for j in range(len(columns)):
            table[(i + 1, j)].set_facecolor(color)

    ax.set_title(f"{title}\n(n = {adata.n_obs:,} cells)", fontsize=14, fontweight='bold', pad=20)

    fig.text(0.5, 0.02,
             "Interpretation: Positive Global L = spatial co-localization | "
             "Negative Global L = spatial anti-correlation | "
             "HH = High-High quadrant (both genes high)",
             ha='center', fontsize=9, style='italic', color='gray')

    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved summary table to {output_path}")
```

---

## 5. Standalone Scripts

### 5.1 Complete Moran's I Analysis Script

```python
#!/usr/bin/env python3
"""
Standalone Moran's I Analysis Script

This script performs spatial autocorrelation analysis using Moran's I
without any Temporal or process_activity dependencies.

Usage:
    python morans_i_analysis.py input.h5ad output_dir/
"""

import sys
import numpy as np
import scanpy as sc
import squidpy as sq
from scipy import sparse
from pathlib import Path
import matplotlib.pyplot as plt


def main(input_path, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    # Load data
    print(f"Loading {input_path}...")
    adata = sc.read_h5ad(input_path)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # Compute spatial neighbors
    print("Computing spatial neighbors...")
    sq.gr.spatial_neighbors(adata, n_neighs=6, coord_type="generic")

    # Compute global Moran's I
    print("Computing global Moran's I for top 50 genes...")
    gene_means = np.array(adata.X.mean(axis=0)).flatten()
    top_genes = adata.var_names[np.argsort(gene_means)[-50:]].tolist()

    adata_sub = adata[:, top_genes].copy()
    sq.gr.spatial_autocorr(adata_sub, mode="moran", n_perms=200)
    adata.uns['moranI'] = adata_sub.uns['moranI']

    # Print results
    morans_df = adata.uns['moranI']
    significant = morans_df[morans_df['pval_sim'] < 0.05].sort_values('I', ascending=False)
    print(f"\nFound {len(significant)} significant spatially variable genes")
    print("\nTop 10 genes:")
    print(significant.head(10)[['I', 'pval_sim']])

    # Compute local Moran's I for top gene
    top_gene = morans_df['I'].idxmax()
    print(f"\nComputing local Moran's I for top gene: {top_gene}")

    W = adata.obsp['spatial_connectivities']
    if sparse.issparse(W):
        row_sums = np.array(W.sum(axis=1)).flatten()
        row_sums[row_sums == 0] = 1
        W_norm = W.multiply(1.0 / row_sums[:, np.newaxis])
    else:
        row_sums = W.sum(axis=1)
        row_sums[row_sums == 0] = 1
        W_norm = W / row_sums[:, np.newaxis]

    gene_idx = adata.var_names.get_loc(top_gene)
    x = adata.X[:, gene_idx].toarray().flatten() if hasattr(adata.X, 'toarray') else adata.X[:, gene_idx].flatten()
    z = (x - x.mean()) / x.std()
    lag_z = np.array(W_norm @ z).flatten() if sparse.issparse(W_norm) else W_norm @ z
    local_I = z * lag_z
    adata.obs[f'{top_gene}_local_morans_i'] = local_I

    # Create visualizations
    print("\nCreating visualizations...")

    # Global plot
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    axes[0, 0].hist(morans_df['I'], bins=30, edgecolor='black', alpha=0.7)
    axes[0, 0].axvline(0, color='red', linestyle='--')
    axes[0, 0].set_xlabel("Moran's I")
    axes[0, 0].set_title("Distribution of Moran's I")

    sig_mask = morans_df['pval_sim'] < 0.05
    axes[0, 1].scatter(gene_means[~sig_mask], morans_df['I'][~sig_mask], alpha=0.3, c='gray')
    axes[0, 1].scatter(gene_means[sig_mask], morans_df['I'][sig_mask], alpha=0.7, c='red')
    axes[0, 1].set_xlabel("Mean Expression")
    axes[0, 1].set_ylabel("Moran's I")
    axes[0, 1].set_title("Moran's I vs Expression")

    coords = adata.obsm['spatial']
    sc1 = axes[1, 0].scatter(coords[:, 0], coords[:, 1], c=x, s=1, cmap='viridis')
    axes[1, 0].set_title(f"{top_gene} Expression")
    axes[1, 0].set_aspect('equal')
    axes[1, 0].invert_yaxis()
    plt.colorbar(sc1, ax=axes[1, 0])

    vmax = np.abs(local_I).max()
    sc2 = axes[1, 1].scatter(coords[:, 0], coords[:, 1], c=local_I, s=1, cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    axes[1, 1].set_title(f"{top_gene} Local Moran's I")
    axes[1, 1].set_aspect('equal')
    axes[1, 1].invert_yaxis()
    plt.colorbar(sc2, ax=axes[1, 1])

    plt.tight_layout()
    plt.savefig(output_dir / "morans_i_analysis.png", dpi=300, bbox_inches='tight')
    plt.close()

    # Save results
    adata.write_h5ad(output_dir / "output.h5ad")
    morans_df.to_csv(output_dir / "morans_i_results.csv")

    print(f"\nSaved results to {output_dir}")
    print("  - morans_i_analysis.png")
    print("  - morans_i_results.csv")
    print("  - output.h5ad")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python morans_i_analysis.py input.h5ad output_dir/")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
```

### 5.2 Complete Lee's L Analysis Script

```python
#!/usr/bin/env python3
"""
Standalone Lee's L Analysis Script

This script performs bivariate spatial autocorrelation analysis using Lee's L
with quadrant classification, without any Temporal or process_activity dependencies.

Usage:
    python lees_l_analysis.py input.h5ad output_dir/ --gene1 CD38 --gene2 CDKN1A,CCL3,CCL2
"""

import sys
import argparse
import numpy as np
import scanpy as sc
import squidpy as sq
from scipy import sparse
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd


def compute_local_lees_l(adata, gene1, gene2, n_neighbors=6):
    """Compute local Lee's L with quadrant classification."""

    if 'spatial_connectivities' not in adata.obsp:
        sq.gr.spatial_neighbors(adata, n_neighs=n_neighbors, coord_type="generic")

    W = adata.obsp['spatial_connectivities']
    if sparse.issparse(W):
        row_sums = np.array(W.sum(axis=1)).flatten()
        row_sums[row_sums == 0] = 1
        W_norm = W.multiply(1.0 / row_sums[:, np.newaxis])
    else:
        row_sums = W.sum(axis=1)
        row_sums[row_sums == 0] = 1
        W_norm = W / row_sums[:, np.newaxis]

    gene1_idx = adata.var_names.get_loc(gene1)
    gene2_idx = adata.var_names.get_loc(gene2)

    if hasattr(adata.X, 'toarray'):
        x = adata.X[:, gene1_idx].toarray().flatten()
        y = adata.X[:, gene2_idx].toarray().flatten()
    else:
        x = adata.X[:, gene1_idx].flatten()
        y = adata.X[:, gene2_idx].flatten()

    x_std, y_std = x.std(), y.std()

    if x_std == 0 or y_std == 0:
        return np.zeros(len(x)), np.array(["NS"] * len(x)), {}

    z_x = (x - x.mean()) / x_std
    z_y = (y - y.mean()) / y_std

    if sparse.issparse(W_norm):
        lag_z_y = np.array(W_norm @ z_y).flatten()
    else:
        lag_z_y = W_norm @ z_y

    local_L = z_x * lag_z_y

    quadrant = np.empty(len(x), dtype='U2')
    quadrant[:] = "NS"
    quadrant[(z_x > 0) & (lag_z_y > 0)] = "HH"
    quadrant[(z_x < 0) & (lag_z_y < 0)] = "LL"
    quadrant[(z_x > 0) & (lag_z_y < 0)] = "HL"
    quadrant[(z_x < 0) & (lag_z_y > 0)] = "LH"

    # Store in adata
    col_prefix = f"{gene1}_{gene2}"
    adata.obs[f"{col_prefix}_local_lees_l"] = local_L
    adata.obs[f"{col_prefix}_quadrant"] = quadrant

    summary = {
        "gene1": gene1,
        "gene2": gene2,
        "global_L": float(local_L.sum()),
        "mean_L": float(local_L.mean()),
        "n_HH": int((quadrant == "HH").sum()),
        "n_LL": int((quadrant == "LL").sum()),
        "n_HL": int((quadrant == "HL").sum()),
        "n_LH": int((quadrant == "LH").sum()),
    }

    return local_L, quadrant, summary


def main():
    parser = argparse.ArgumentParser(description="Lee's L Spatial Analysis")
    parser.add_argument("input", help="Input h5ad file")
    parser.add_argument("output_dir", help="Output directory")
    parser.add_argument("--gene1", required=True, help="Primary gene (e.g., CD38)")
    parser.add_argument("--gene2", required=True, help="Comma-separated list of secondary genes")
    parser.add_argument("--n_neighbors", type=int, default=6, help="Number of spatial neighbors")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True)

    gene2_list = [g.strip() for g in args.gene2.split(",")]

    # Load data
    print(f"Loading {args.input}...")
    adata = sc.read_h5ad(args.input)
    print(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # Validate genes
    missing = [g for g in [args.gene1] + gene2_list if g not in adata.var_names]
    if missing:
        print(f"ERROR: Genes not found: {missing}")
        sys.exit(1)

    # Compute Lee's L for each pair
    results = []
    for gene2 in gene2_list:
        print(f"\nComputing Lee's L: {args.gene1} vs {gene2}")
        local_L, quadrant, summary = compute_local_lees_l(adata, args.gene1, gene2, args.n_neighbors)
        results.append(summary)

        print(f"  Global L: {summary['global_L']:.4f}")
        print(f"  HH: {summary['n_HH']} ({100*summary['n_HH']/len(quadrant):.1f}%)")

    # Create summary table
    results_df = pd.DataFrame(results)
    results_df['gene_pair'] = results_df['gene1'] + ':' + results_df['gene2']
    results_df['HH_pct'] = 100 * results_df['n_HH'] / adata.n_obs
    results_df = results_df.sort_values('global_L', ascending=False)

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(results_df[['gene_pair', 'global_L', 'mean_L', 'n_HH', 'HH_pct']].to_string(index=False))

    # Create visualizations
    gene_pairs = [(args.gene1, g2) for g2 in gene2_list]

    # Summary table plot
    fig, ax = plt.subplots(figsize=(14, 2 + len(results) * 0.5))
    ax.axis('off')

    columns = ["Gene Pair", "Global L", "Mean L", "HH", "LL", "HL", "LH", "HH %"]
    cell_text = []
    for r in results:
        cell_text.append([
            f"{r['gene1']}:{r['gene2']}",
            f"{r['global_L']:,.2f}",
            f"{r['mean_L']:.6f}",
            f"{r['n_HH']:,}",
            f"{r['n_LL']:,}",
            f"{r['n_HL']:,}",
            f"{r['n_LH']:,}",
            f"{100*r['n_HH']/adata.n_obs:.1f}%",
        ])

    table = ax.table(cellText=cell_text, colLabels=columns, cellLoc='center', loc='center',
                     colColours=['#4472C4'] * len(columns))
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1.2, 1.8)

    for j in range(len(columns)):
        table[(0, j)].set_text_props(weight='bold', color='white')

    ax.set_title(f"Lee's L Summary\n(n = {adata.n_obs:,} cells)", fontsize=14, fontweight='bold', pad=20)

    plt.tight_layout()
    plt.savefig(output_dir / "lees_l_summary_table.png", dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()

    # Save results
    adata.write_h5ad(output_dir / "output.h5ad")
    results_df.to_csv(output_dir / "lees_l_results.csv", index=False)

    print(f"\nSaved results to {output_dir}")


if __name__ == "__main__":
    main()
```

---

## 6. Test Examples

### 6.1 Test Moran's I Computation

```python
import numpy as np
import pytest


def test_local_morans_i_basic():
    """Test that local Moran's I computes correctly."""
    from scipy import sparse

    # Create simple test data
    np.random.seed(42)
    n_cells = 100

    # Create spatial weights (ring topology)
    W = np.zeros((n_cells, n_cells))
    for i in range(n_cells):
        W[i, (i-1) % n_cells] = 1
        W[i, (i+1) % n_cells] = 1
    W = sparse.csr_matrix(W)

    # Row normalize
    row_sums = np.array(W.sum(axis=1)).flatten()
    W_norm = W.multiply(1.0 / row_sums[:, np.newaxis])

    # Create clustered expression (first half high, second half low)
    x = np.concatenate([np.ones(50) * 2, np.zeros(50)])
    z = (x - x.mean()) / x.std()

    # Compute local I
    lag_z = np.array(W_norm @ z).flatten()
    local_I = z * lag_z

    # Cells in the middle of clusters should have positive local I
    assert local_I[25] > 0  # Middle of high cluster
    assert local_I[75] > 0  # Middle of low cluster

    # Cells at boundaries should have negative local I
    assert local_I[0] < 0  # Boundary between clusters


def test_global_morans_i_random():
    """Test that random data gives Moran's I close to 0."""
    import scanpy as sc
    import squidpy as sq

    # Create random spatial data
    np.random.seed(42)
    n_cells = 500

    adata = sc.AnnData(X=np.random.randn(n_cells, 10))
    adata.obsm['spatial'] = np.random.rand(n_cells, 2) * 100
    adata.var_names = [f"gene_{i}" for i in range(10)]

    sq.gr.spatial_neighbors(adata, n_neighs=6, coord_type="generic")
    sq.gr.spatial_autocorr(adata, mode="moran", n_perms=50)

    # All Moran's I should be close to 0 for random data
    morans_I = adata.uns['moranI']['I'].values
    assert np.abs(morans_I).mean() < 0.1
```

### 6.2 Test Lee's L Computation

```python
def test_lees_l_perfect_correlation():
    """Test Lee's L with perfectly correlated genes."""
    import scanpy as sc
    import squidpy as sq
    from scipy import sparse

    np.random.seed(42)
    n_cells = 100

    # Create spatial coordinates
    coords = np.random.rand(n_cells, 2) * 100

    # Create perfectly correlated expression
    gene1 = np.random.rand(n_cells)
    gene2 = gene1.copy()  # Perfect correlation

    adata = sc.AnnData(X=np.column_stack([gene1, gene2]))
    adata.obsm['spatial'] = coords
    adata.var_names = ['gene1', 'gene2']

    sq.gr.spatial_neighbors(adata, n_neighs=6, coord_type="generic")

    # Compute Lee's L
    W = adata.obsp['spatial_connectivities']
    row_sums = np.array(W.sum(axis=1)).flatten()
    row_sums[row_sums == 0] = 1
    W_norm = W.multiply(1.0 / row_sums[:, np.newaxis])

    z_x = (gene1 - gene1.mean()) / gene1.std()
    z_y = (gene2 - gene2.mean()) / gene2.std()
    lag_z_y = np.array(W_norm @ z_y).flatten()

    local_L = z_x * lag_z_y
    global_L = local_L.sum()

    # Global L should be strongly positive for perfect correlation
    assert global_L > 0


def test_quadrant_classification():
    """Test quadrant classification logic."""

    # Test cases: (z_x, lag_z_y) -> expected quadrant
    test_cases = [
        (1.0, 1.0, "HH"),   # High-High
        (-1.0, -1.0, "LL"), # Low-Low
        (1.0, -1.0, "HL"),  # High-Low
        (-1.0, 1.0, "LH"),  # Low-High
        (0.0, 0.0, "NS"),   # Not significant (at zero)
    ]

    for z_x, lag_z_y, expected in test_cases:
        if z_x > 0 and lag_z_y > 0:
            quadrant = "HH"
        elif z_x < 0 and lag_z_y < 0:
            quadrant = "LL"
        elif z_x > 0 and lag_z_y < 0:
            quadrant = "HL"
        elif z_x < 0 and lag_z_y > 0:
            quadrant = "LH"
        else:
            quadrant = "NS"

        assert quadrant == expected, f"Expected {expected}, got {quadrant} for ({z_x}, {lag_z_y})"


def test_hh_proportion_calculation():
    """Test HH proportion by cell type calculation."""
    import pandas as pd

    # Create mock data
    n = 100
    quadrants = np.array(["HH"] * 30 + ["LL"] * 40 + ["HL"] * 15 + ["LH"] * 15)
    cell_types = np.array(["TypeA"] * 50 + ["TypeB"] * 50)

    # Calculate HH proportion
    results = []
    for ct in ["TypeA", "TypeB"]:
        mask = cell_types == ct
        n_cells = mask.sum()
        n_hh = ((quadrants == "HH") & mask).sum()
        hh_prop = n_hh / n_cells
        results.append({"cell_type": ct, "HH_proportion": hh_prop})

    df = pd.DataFrame(results)

    # TypeA: first 50 cells, 30 HH -> 60%
    assert abs(df[df['cell_type'] == 'TypeA']['HH_proportion'].values[0] - 0.6) < 0.01

    # TypeB: last 50 cells, 0 HH -> 0%
    assert df[df['cell_type'] == 'TypeB']['HH_proportion'].values[0] == 0.0
```

---

## Dependencies

```
# Core
numpy>=1.20.0
scipy>=1.7.0
pandas>=1.3.0
scanpy>=1.9.0
squidpy>=1.2.0
anndata>=0.8.0

# Visualization
matplotlib>=3.5.0
seaborn>=0.11.0

# Optional (for enhanced visualizations)
adjustText>=0.7.3
```

---

## References

1. **Moran's I**: Moran, P. A. P. (1950). "Notes on Continuous Stochastic Phenomena". Biometrika. 37 (1): 17–23.

2. **Lee's L**: Lee, S. I. (2001). "Developing a bivariate spatial association measure: An integration of Pearson's r and Moran's I". Journal of Geographical Systems. 3 (4): 369–385.

3. **Squidpy**: Palla, G., Spitzer, H., Klein, M. et al. (2022). "Squidpy: a scalable framework for spatial omics analysis". Nat Methods 19, 171–178.

4. **Local Indicators of Spatial Association (LISA)**: Anselin, L. (1995). "Local Indicators of Spatial Association—LISA". Geographical Analysis, 27(2), 93-115.
