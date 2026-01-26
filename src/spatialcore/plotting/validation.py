"""
Marker validation visualization.

This module provides functions for validating cell type annotations
using canonical marker genes and GMM-3 thresholding.

For marker validation, we use classify_by_threshold with n_components=3
(trimodal GMM) which handles spatial data's dropout/moderate/high
expression patterns better than bimodal GMM.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import gc

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import anndata as ad

from spatialcore.core.logging import get_logger
from spatialcore.plotting.utils import (
    generate_celltype_palette,
    setup_figure,
    setup_multi_figure,
    save_figure,
    despine,
)

logger = get_logger(__name__)


def plot_marker_heatmap(
    adata: ad.AnnData,
    label_column: str,
    markers: Optional[Dict[str, List[str]]] = None,
    cluster: bool = True,
    layer: Optional[str] = None,
    figsize: Optional[tuple] = None,
    cmap: str = "RdBu_r",
    center: float = 0,
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot marker gene expression heatmap by cell type.

    Parameters
    ----------
    adata : AnnData
        Annotated data with cell type labels.
    label_column : str
        Column in adata.obs containing cell type labels.
    markers : Dict[str, List[str]], optional
        Marker genes per cell type. If None, uses canonical markers.
    cluster : bool, default True
        Hierarchically cluster cell types.
    layer : str, optional
        Layer to use. If None, uses adata.X.
    figsize : tuple, optional
        Figure size. Auto-calculated if None.
    cmap : str, default "RdBu_r"
        Colormap.
    center : float, default 0
        Value to center colormap on.
    title : str, optional
        Plot title.
    save : str or Path, optional
        Path to save figure.

    Returns
    -------
    Figure
        Matplotlib figure.

    Examples
    --------
    >>> from spatialcore.plotting.validation import plot_marker_heatmap
    >>> from spatialcore.annotation.markers import CANONICAL_MARKERS
    >>> fig = plot_marker_heatmap(
    ...     adata,
    ...     label_column="cell_type",
    ...     markers=CANONICAL_MARKERS,
    ... )
    """
    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("seaborn is required for heatmaps")

    if label_column not in adata.obs.columns:
        raise ValueError(f"Label column '{label_column}' not found.")

    # Load canonical markers if not provided
    if markers is None:
        from spatialcore.annotation.markers import load_canonical_markers
        markers = load_canonical_markers()

    # Collect all marker genes that exist in the data
    all_genes = []
    cell_types_with_markers = []

    for cell_type in adata.obs[label_column].unique():
        ct_lower = str(cell_type).lower()
        ct_markers = markers.get(ct_lower, [])
        available = [g for g in ct_markers if g in adata.var_names]
        if available:
            all_genes.extend(available)
            cell_types_with_markers.append(cell_type)

    all_genes = list(dict.fromkeys(all_genes))  # Remove duplicates, preserve order

    if not all_genes:
        raise ValueError("No marker genes found in data.")

    # Calculate mean expression per cell type
    cell_types = adata.obs[label_column].unique()
    mean_expr = pd.DataFrame(index=cell_types, columns=all_genes, dtype=float)

    for ct in cell_types:
        mask = adata.obs[label_column] == ct
        # Use integer indexing to avoid anndata 0.12.x boolean mask bug
        mask_indices = np.where(mask)[0]
        gene_indices = [adata.var_names.get_loc(g) for g in all_genes]
        if layer:
            X = adata.layers[layer][mask_indices][:, gene_indices]
        else:
            X = adata.X[mask_indices][:, gene_indices]
        if hasattr(X, "toarray"):
            X = X.toarray()
        mean_expr.loc[ct] = np.mean(X, axis=0)

    # Z-score normalize columns
    mean_expr_z = (mean_expr - mean_expr.mean()) / mean_expr.std()
    mean_expr_z = mean_expr_z.fillna(0)

    # Figure size
    if figsize is None:
        n_types = len(cell_types)
        n_genes = len(all_genes)
        figsize = (max(10, n_genes * 0.3), max(6, n_types * 0.4))

    # Cluster if requested
    if cluster:
        g = sns.clustermap(
            mean_expr_z,
            cmap=cmap,
            center=center,
            figsize=figsize,
            row_cluster=True,
            col_cluster=True,
        )
        if title:
            g.fig.suptitle(title, y=1.02)
        if save:
            g.savefig(save)
        return g.fig
    else:
        fig, ax = setup_figure(figsize=figsize)
        sns.heatmap(
            mean_expr_z,
            cmap=cmap,
            center=center,
            ax=ax,
            xticklabels=True,
            yticklabels=True,
        )
        ax.set_xlabel("Marker Genes")
        ax.set_ylabel("Cell Types")

        if title is None:
            title = "Marker Expression Heatmap"
        ax.set_title(title)

        plt.tight_layout()

        if save:
            save_figure(fig, save)

        return fig


def plot_2d_validation(
    adata: ad.AnnData,
    label_column: str,
    confidence_column: str,
    markers: Optional[Dict[str, List[str]]] = None,
    confidence_threshold: float = 0.8,
    min_cells_per_type: int = 15,
    n_components: int = 3,
    ncols: int = 4,
    figsize_per_panel: Tuple[float, float] = (3, 3),
    save: Optional[Union[str, Path]] = None,
) -> Tuple[Figure, pd.DataFrame]:
    """
    2D marker validation plot per cell type.

    For each cell type, plots confidence (x-axis) vs marker metagene score
    (y-axis). Cells are colored green if above both thresholds, red otherwise.

    Uses GMM-3 via classify_by_threshold() to find marker threshold.

    Parameters
    ----------
    adata : AnnData
        Annotated data with cell type labels and confidence.
    label_column : str
        Column in adata.obs containing cell type labels.
    confidence_column : str
        Column in adata.obs containing confidence values.
    markers : Dict[str, List[str]], optional
        Marker genes per cell type. If None, uses canonical markers.
    confidence_threshold : float, default 0.8
        Confidence threshold for validation.
    min_cells_per_type : int, default 15
        Minimum cells required to plot a cell type.
    n_components : int, default 3
        Number of GMM components (3 for trimodal spatial data).
    ncols : int, default 4
        Number of columns in subplot grid.
    figsize_per_panel : Tuple[float, float], default (3, 3)
        Size per panel.
    save : str or Path, optional
        Path to save figure.

    Returns
    -------
    Tuple[Figure, pd.DataFrame]
        Figure and validation summary DataFrame.

    Examples
    --------
    >>> from spatialcore.plotting.validation import plot_2d_validation
    >>> fig, summary = plot_2d_validation(
    ...     adata,
    ...     label_column="cell_type",
    ...     confidence_column="confidence",
    ...     confidence_threshold=0.7,
    ... )
    >>> print(summary)
    """
    from spatialcore.stats.classify import classify_by_threshold

    if label_column not in adata.obs.columns:
        raise ValueError(f"Label column '{label_column}' not found.")
    if confidence_column not in adata.obs.columns:
        raise ValueError(f"Confidence column '{confidence_column}' not found.")

    # Load canonical markers if not provided
    if markers is None:
        from spatialcore.annotation.markers import load_canonical_markers
        markers = load_canonical_markers()

    # Find cell types with enough cells and available markers
    cell_types = adata.obs[label_column].value_counts()
    cell_types = cell_types[cell_types >= min_cells_per_type].index.tolist()

    types_to_plot = []
    for ct in cell_types:
        ct_lower = str(ct).lower()
        ct_markers = markers.get(ct_lower, [])
        available = [g for g in ct_markers if g in adata.var_names]
        if len(available) >= 2:  # Need at least 2 markers
            types_to_plot.append((ct, available))

    if not types_to_plot:
        logger.warning("No cell types with sufficient markers found for 2D validation.")
        fig, ax = setup_figure(figsize=figsize_per_panel)
        ax.axis("off")
        ax.text(
            0.5, 0.5,
            "No cell types with sufficient markers",
            ha="center", va="center",
        )
        if save:
            save_figure(fig, save)
        summary_df = pd.DataFrame(columns=[
            "cell_type",
            "n_cells",
            "n_low_conf",
            "n_high_conf",
            "n_validated",
            "pct_validated",
            "pct_high_conf",
            "marker_threshold",
            "n_markers",
        ])
        return fig, summary_df

    # Pre-compute GMM results for all cell types to know which will succeed
    successful_types = []
    for cell_type, ct_markers in types_to_plot:
        mask = adata.obs[label_column] == cell_type
        # Manual AnnData construction to avoid anndata 0.12.x bug in adata[mask].copy()
        mask_indices = np.where(mask)[0]
        subset = ad.AnnData(
            X=adata.X[mask_indices].copy() if hasattr(adata.X[mask_indices], 'copy') else adata.X[mask_indices],
            obs=adata.obs.iloc[mask_indices].copy(),
            var=adata.var.copy(),
        )
        if 'norm' in adata.layers:
            subset.layers['norm'] = adata.layers['norm'][mask_indices]
        confidence = subset.obs[confidence_column].values

        try:
            subset_result = classify_by_threshold(
                subset,
                feature_columns=ct_markers,
                threshold_method="gmm",
                n_components=n_components,
                metagene_method="shifted_geometric_mean",
                plot=False,
                copy=True,
            )
        except Exception as e:
            logger.warning(f"Skipping {cell_type}: GMM failed - {e}")
            continue

        metagene_scores = subset_result.obs["threshold_score"].values
        marker_threshold = subset_result.uns.get(
            "threshold_params", {}
        ).get("threshold", np.median(metagene_scores))

        successful_types.append({
            "cell_type": cell_type,
            "markers": ct_markers,
            "subset": subset,
            "confidence": confidence,
            "metagene_scores": metagene_scores,
            "marker_threshold": marker_threshold,
        })

    if not successful_types:
        logger.warning("No cell types passed GMM validation for 2D plot.")
        fig, ax = setup_figure(figsize=figsize_per_panel)
        ax.axis("off")
        ax.text(
            0.5, 0.5,
            "No cell types passed GMM validation",
            ha="center", va="center",
        )
        if save:
            save_figure(fig, save)
        summary_df = pd.DataFrame(columns=[
            "cell_type",
            "n_cells",
            "n_low_conf",
            "n_high_conf",
            "n_validated",
            "pct_validated",
            "pct_high_conf",
            "marker_threshold",
            "n_markers",
        ])
        return fig, summary_df

    # Create grid with only successful cell types
    n_types = len(successful_types)
    nrows = int(np.ceil(n_types / ncols))
    figsize = (figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = np.atleast_2d(axes).flatten()

    summary_rows = []

    for i, data in enumerate(successful_types):
        ax = axes[i]
        cell_type = data["cell_type"]
        ct_markers = data["markers"]
        subset = data["subset"]
        confidence = data["confidence"]
        metagene_scores = data["metagene_scores"]
        marker_threshold = data["marker_threshold"]

        # Classify cells into three groups per spec section 3.1:
        # - Red: Low confidence (uncertain)
        # - Green: High confidence only (validated)
        # - Yellow: High conf + High marker (strongly validated)
        high_conf = confidence >= confidence_threshold
        high_marker = metagene_scores >= marker_threshold

        low_conf = ~high_conf  # Red
        high_conf_low_marker = high_conf & ~high_marker  # Green
        high_conf_high_marker = high_conf & high_marker  # Yellow

        # Plot in order: red (back), green, yellow (front)
        ax.scatter(
            confidence[low_conf],
            metagene_scores[low_conf],
            c="red",
            s=5,
            alpha=0.5,
            label="Low Conf",
        )
        ax.scatter(
            confidence[high_conf_low_marker],
            metagene_scores[high_conf_low_marker],
            c="green",
            s=5,
            alpha=0.5,
            label="High Conf",
        )
        ax.scatter(
            confidence[high_conf_high_marker],
            metagene_scores[high_conf_high_marker],
            c="gold",
            s=5,
            alpha=0.5,
            label="High Conf + Marker",
        )

        # Threshold lines
        ax.axvline(confidence_threshold, color="gray", linestyle="--", alpha=0.5)
        ax.axhline(marker_threshold, color="gray", linestyle="--", alpha=0.5)

        ax.set_xlabel("Confidence")
        ax.set_ylabel("Marker Score")
        ax.set_xlim(0, 1)  # Fixed x-axis for consistent comparison
        ax.set_xticks([0, 0.5, 1.0])
        ax.set_title(f"{cell_type}\n(n={len(subset)})", fontsize=9)

        # Summary stats
        n_low_conf = low_conf.sum()
        n_high_conf = high_conf_low_marker.sum()
        n_validated = high_conf_high_marker.sum()
        pct_validated = 100 * n_validated / len(subset) if len(subset) > 0 else 0
        pct_high_conf = 100 * (n_high_conf + n_validated) / len(subset) if len(subset) > 0 else 0

        summary_rows.append({
            "cell_type": cell_type,
            "n_cells": len(subset),
            "n_low_conf": n_low_conf,
            "n_high_conf": n_high_conf,
            "n_validated": n_validated,
            "pct_validated": pct_validated,
            "pct_high_conf": pct_high_conf,
            "marker_threshold": marker_threshold,
            "n_markers": len(ct_markers),
        })

    # Hide unused panels (when n_types doesn't fill full grid)
    for i in range(len(successful_types), len(axes)):
        axes[i].set_visible(False)

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    summary_df = pd.DataFrame(summary_rows)
    return fig, summary_df


def plot_marker_dotplot(
    adata: ad.AnnData,
    label_column: str,
    markers: Optional[Dict[str, List[str]]] = None,
    layer: Optional[str] = None,
    figsize: Optional[tuple] = None,
    cmap: str = "Reds",
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot marker expression as dot plot.

    Dot size represents fraction of cells expressing the marker,
    color intensity represents mean expression.

    Parameters
    ----------
    adata : AnnData
        Annotated data with cell type labels.
    label_column : str
        Column in adata.obs containing cell type labels.
    markers : Dict[str, List[str]], optional
        Marker genes per cell type.
    layer : str, optional
        Layer to use.
    figsize : tuple, optional
        Figure size.
    cmap : str, default "Reds"
        Colormap.
    title : str, optional
        Plot title.
    save : str or Path, optional
        Path to save figure.

    Returns
    -------
    Figure
        Matplotlib figure.
    """
    if label_column not in adata.obs.columns:
        raise ValueError(f"Label column '{label_column}' not found.")

    # Load canonical markers if not provided
    if markers is None:
        from spatialcore.annotation.markers import load_canonical_markers
        markers = load_canonical_markers()

    # Collect all marker genes
    all_genes = []
    for ct in adata.obs[label_column].unique():
        ct_lower = str(ct).lower()
        ct_markers = markers.get(ct_lower, [])
        available = [g for g in ct_markers if g in adata.var_names]
        all_genes.extend(available)

    all_genes = list(dict.fromkeys(all_genes))

    if not all_genes:
        raise ValueError("No marker genes found in data.")

    cell_types = sorted(adata.obs[label_column].unique())

    # Calculate expression fraction and mean
    n_types = len(cell_types)
    n_genes = len(all_genes)

    frac_expr = np.zeros((n_types, n_genes))
    mean_expr = np.zeros((n_types, n_genes))

    # Pre-compute gene indices once (outside loop)
    gene_indices = [adata.var_names.get_loc(g) for g in all_genes]

    for i, ct in enumerate(cell_types):
        mask = adata.obs[label_column] == ct
        # Use integer indexing to avoid anndata 0.12.x boolean mask bug
        mask_indices = np.where(mask)[0]
        if layer:
            X = adata.layers[layer][mask_indices][:, gene_indices]
        else:
            X = adata.X[mask_indices][:, gene_indices]
        if hasattr(X, "toarray"):
            X = X.toarray()

        frac_expr[i] = (X > 0).mean(axis=0)
        mean_expr[i] = X.mean(axis=0)

    # Normalize mean expression per gene
    mean_expr_norm = (mean_expr - mean_expr.min(axis=0)) / (
        mean_expr.max(axis=0) - mean_expr.min(axis=0) + 1e-10
    )

    if figsize is None:
        figsize = (max(10, n_genes * 0.4), max(6, n_types * 0.4))

    fig, ax = setup_figure(figsize=figsize)

    # Create dot plot
    for i, ct in enumerate(cell_types):
        for j, gene in enumerate(all_genes):
            size = frac_expr[i, j] * 300  # Scale for visibility
            color = mean_expr_norm[i, j]
            ax.scatter(j, i, s=size, c=[color], cmap=cmap, vmin=0, vmax=1)

    ax.set_xticks(range(n_genes))
    ax.set_xticklabels(all_genes, rotation=90)
    ax.set_yticks(range(n_types))
    ax.set_yticklabels(cell_types)

    ax.set_xlabel("Marker Genes")
    ax.set_ylabel("Cell Types")

    if title is None:
        title = "Marker Expression Dotplot"
    ax.set_title(title)

    # Add legend for size
    size_legend = [0.25, 0.5, 0.75, 1.0]
    legend_x = n_genes + 0.5
    for k, frac in enumerate(size_legend):
        ax.scatter(legend_x, k, s=frac * 300, c="gray")
        ax.text(legend_x + 0.3, k, f"{int(frac*100)}%", va="center")
    ax.text(legend_x, len(size_legend) + 0.5, "% Expressing", fontsize=9)

    ax.set_xlim(-0.5, n_genes + 2)

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    return fig


def plot_celltype_confidence(
    adata: ad.AnnData,
    label_column: str,
    confidence_column: str,
    spatial_key: str = "spatial",
    threshold: float = 0.8,
    max_cell_types: int = 20,
    figsize: Tuple[float, float] = (14, 6),
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Cell type confidence visualization with spatial and jitter plots.

    Creates a two-panel figure:
    - Left panel: Spatial scatter plot colored by cell type
    - Right panel: Jitter plot (x=cell type, y=confidence) with threshold line

    Parameters
    ----------
    adata : AnnData
        Annotated data with cell type labels and confidence scores.
    label_column : str
        Column in adata.obs containing cell type labels.
    confidence_column : str
        Column in adata.obs containing confidence values.
    spatial_key : str, default "spatial"
        Key in adata.obsm for spatial coordinates.
    threshold : float, default 0.8
        Confidence threshold line to display.
    max_cell_types : int, default 20
        Maximum number of cell types to show in jitter plot.
    figsize : Tuple[float, float], default (14, 6)
        Figure size (width, height).
    save : str or Path, optional
        Path to save figure.

    Returns
    -------
    Figure
        Matplotlib figure.

    Examples
    --------
    >>> from spatialcore.plotting.validation import plot_celltype_confidence
    >>> fig = plot_celltype_confidence(
    ...     adata,
    ...     label_column="cell_type",
    ...     confidence_column="cell_type_confidence",
    ...     threshold=0.8,
    ... )
    """
    if label_column not in adata.obs.columns:
        raise ValueError(f"Label column '{label_column}' not found.")
    if confidence_column not in adata.obs.columns:
        raise ValueError(f"Confidence column '{confidence_column}' not found.")

    fig, (ax_spatial, ax_jitter) = plt.subplots(1, 2, figsize=figsize)

    # Get cell types and palette
    cell_types = adata.obs[label_column].unique()
    colors = generate_celltype_palette(cell_types)

    # Left: Spatial plot colored by confidence z-score
    if spatial_key in adata.obsm:
        coords = adata.obsm[spatial_key]
        conf_values = adata.obs[confidence_column].values

        # Create scatter colored by confidence
        scatter = ax_spatial.scatter(
            coords[:, 0],
            coords[:, 1],
            s=1,
            alpha=0.7,
            c=conf_values,
            cmap="RdYlGn",  # Red (low) -> Yellow -> Green (high)
            vmin=np.nanpercentile(conf_values, 5),
            vmax=np.nanpercentile(conf_values, 95),
            rasterized=True,
        )
        cbar = plt.colorbar(scatter, ax=ax_spatial, shrink=0.7, pad=0.02)
        cbar.set_label("Confidence (z-score)", fontsize=9)
        ax_spatial.set_title("Spatial Confidence")
        ax_spatial.set_xlabel("X")
        ax_spatial.set_ylabel("Y")
        ax_spatial.axis("equal")
        despine(ax_spatial)
    else:
        ax_spatial.text(
            0.5, 0.5,
            f"No spatial coordinates found\n(key: '{spatial_key}')",
            ha="center", va="center", transform=ax_spatial.transAxes
        )
        ax_spatial.set_title("Spatial Confidence (N/A)")

    # Right: Jitter plot - sort by median confidence (flipped: cell types on Y-axis)
    ct_median = adata.obs.groupby(label_column, observed=True)[confidence_column].median()
    ct_order = ct_median.sort_values(ascending=True).index.tolist()  # Ascending so highest at top

    # Limit number of cell types shown
    if len(ct_order) > max_cell_types:
        # Keep top by median confidence (will be at top of plot)
        ct_order = ct_median.sort_values(ascending=False).index.tolist()[:max_cell_types]
        ct_order = ct_order[::-1]  # Reverse so highest at top
        logger.info(f"Showing top {max_cell_types} cell types by median confidence")

    for i, ct in enumerate(ct_order):
        mask = adata.obs[label_column] == ct
        conf = adata.obs.loc[mask, confidence_column].values
        # Add jitter on Y-axis (cell type position)
        y = np.random.normal(i, 0.15, len(conf))
        color = colors.get(str(ct), "gray")
        ax_jitter.scatter(conf, y, s=3, alpha=0.3, c=[color], rasterized=True)

    # Threshold line (vertical now)
    ax_jitter.axvline(
        threshold, color="red", linestyle="--", lw=2, label=f"Threshold={threshold}"
    )

    ax_jitter.set_yticks(range(len(ct_order)))
    ax_jitter.set_yticklabels(
        [str(ct)[:25] for ct in ct_order],  # Truncate long names
        fontsize=8
    )

    # Set x-axis to fixed 0-1 range with 0.5 increments for consistent comparison
    ax_jitter.set_xlim(0, 1)
    ax_jitter.set_xticks([0, 0.5, 1.0])

    ax_jitter.set_xlabel("Confidence (z-score)")
    ax_jitter.set_ylabel("")
    ax_jitter.set_title("Confidence by Cell Type")
    ax_jitter.legend(loc="lower right")
    despine(ax_jitter)

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    return fig


def plot_deg_heatmap(
    adata: ad.AnnData,
    label_column: str,
    n_genes: int = 5,
    method: str = "wilcoxon",
    layer: Optional[str] = None,
    figsize: Optional[Tuple[float, float]] = None,
    cmap: str = "viridis",
    save: Optional[Union[str, Path]] = None,
    title: Optional[str] = None,
) -> Figure:
    """
    DEG heatmap showing top marker genes per cell type.

    Runs differential expression analysis and displays a heatmap with genes
    on rows, cell types on columns, row color annotation bar, and cell type legend.
    Uses exact plotting code from celltypist_demonstration_plots.py.

    Parameters
    ----------
    adata : AnnData
        Annotated data with cell type labels.
    label_column : str
        Column in adata.obs containing cell type labels.
    n_genes : int, default 5
        Number of top DEGs per cell type to include.
    method : str, default "wilcoxon"
        DEG method ("wilcoxon", "t-test", "t-test_overestim_var", "logreg").
    layer : str, optional
        Expression layer to use. If None, uses adata.X.
    figsize : Tuple[float, float], optional
        Figure size. Auto-calculated if None.
    cmap : str, default "viridis"
        Colormap for heatmap.
    save : str or Path, optional
        Path to save figure.
    title : str, optional
        Plot title. Defaults to "Marker Genes ({label_column})".

    Returns
    -------
    Figure
        Matplotlib figure.
    """
    import scanpy as sc
    import seaborn as sns

    if label_column not in adata.obs.columns:
        raise ValueError(f"Label column '{label_column}' not found.")

    # Build filter mask on original data (no copies)
    mask_assigned = adata.obs[label_column] != "Unassigned"

    # Get valid cell types from view (no copy)
    ct_counts = adata.obs.loc[mask_assigned, label_column].value_counts()
    valid_cts = ct_counts[ct_counts >= 10].index.tolist()
    if len(valid_cts) < 2:
        raise ValueError("Need at least 2 cell types with >= 10 cells each for DEG analysis.")

    # Combined mask - applied once
    mask = mask_assigned & adata.obs[label_column].isin(valid_cts)

    # Manual AnnData construction to avoid anndata 0.12.x memory bug in adata[mask].copy()
    indices = np.where(mask)[0]
    X_sub = adata.X[indices]
    if layer and layer in adata.layers:
        X_sub = adata.layers[layer][indices]
    obs_sub = adata.obs.iloc[indices].copy()
    obs_sub[label_column] = obs_sub[label_column].astype(str).astype("category")
    adata_deg = ad.AnnData(X=X_sub, obs=obs_sub, var=adata.var.copy())

    logger.info(f"Running rank_genes_groups with method={method} on {len(valid_cts)} cell types...")
    sc.tl.rank_genes_groups(adata_deg, label_column, method=method, n_genes=50)

    # Get marker genes per cell type
    results = adata_deg.uns["rank_genes_groups"]
    cell_types = sorted(results["names"].dtype.names)

    # Collect genes grouped by cell type (preserve order)
    all_genes = []
    gene_to_celltype = {}
    seen = set()
    for ct in cell_types:
        genes = results["names"][ct][:n_genes]
        for gene in genes:
            if gene not in seen and gene in adata_deg.var_names:
                all_genes.append(gene)
                gene_to_celltype[gene] = ct
                seen.add(gene)

    if len(all_genes) == 0:
        raise ValueError("No valid marker genes found")

    # Calculate mean expression per cell type
    expr_matrix = np.zeros((len(all_genes), len(cell_types)))
    for j, ct in enumerate(cell_types):
        ct_mask = adata_deg.obs[label_column] == ct
        if ct_mask.sum() == 0:
            continue
        adata_ct = adata_deg[ct_mask]
        for i, gene in enumerate(all_genes):
            gene_idx = adata_deg.var_names.get_loc(gene)
            if hasattr(adata_ct.X, "toarray"):
                expr = adata_ct.X[:, gene_idx].toarray().flatten()
            else:
                expr = np.asarray(adata_ct.X[:, gene_idx]).flatten()
            expr_matrix[i, j] = np.mean(expr)

    # Release the copy - no longer needed after expression matrix is built
    del adata_deg
    gc.collect()

    # Z-score normalize across cell types (rows)
    expr_scaled = np.zeros_like(expr_matrix)
    for i in range(expr_matrix.shape[0]):
        row = expr_matrix[i, :]
        if row.std() > 0:
            expr_scaled[i, :] = (row - row.mean()) / row.std()
        else:
            expr_scaled[i, :] = 0

    # Create color palette for cell types
    n_cts = len(cell_types)
    palette = sns.color_palette("tab20", n_cts)
    ct_to_color = {ct: palette[i] for i, ct in enumerate(cell_types)}

    # Create row colors array
    row_colors = [ct_to_color[gene_to_celltype[gene]] for gene in all_genes]

    # Create figure with gridspec for custom layout
    fig_height = max(10, len(all_genes) * 0.11)
    fig = plt.figure(figsize=(12, fig_height))

    # GridSpec with annotation bar flush against heatmap
    gs = fig.add_gridspec(
        1, 5,
        width_ratios=[0.06, 0.012, 1, 0.02, 0.015],
        wspace=0.0,
        left=0.01,
        right=0.78,
        top=0.95,
        bottom=0.12,
    )

    # Gene labels axis (far left)
    ax_labels = fig.add_subplot(gs[0, 0])
    ax_labels.set_ylim(0, len(all_genes))
    ax_labels.set_xlim(0, 1)
    ax_labels.invert_yaxis()
    for i, gene in enumerate(all_genes):
        ax_labels.text(
            0.98, i + 0.5, gene,
            ha="right", va="center",
            fontsize=5,
            color="black",
        )
    ax_labels.axis("off")

    # Row colors axis (annotation bar)
    ax_rowcolors = fig.add_subplot(gs[0, 1])
    for i, color in enumerate(row_colors):
        ax_rowcolors.add_patch(plt.Rectangle(
            (0, i), 1, 1,
            facecolor=color,
            edgecolor="none",
        ))
    ax_rowcolors.set_xlim(0, 1)
    ax_rowcolors.set_ylim(0, len(all_genes))
    ax_rowcolors.invert_yaxis()
    ax_rowcolors.axis("off")

    # Main heatmap axis
    ax_heatmap = fig.add_subplot(gs[0, 2])
    im = ax_heatmap.imshow(
        expr_scaled,
        aspect="auto",
        cmap=cmap,
        vmin=-2.5,
        vmax=2.5,
    )

    ax_heatmap.set_yticks([])
    ax_heatmap.set_xticks(range(len(cell_types)))
    ax_heatmap.set_xticklabels(cell_types, rotation=45, ha="right", fontsize=8)

    # Gap between heatmap and colorbar
    ax_gap = fig.add_subplot(gs[0, 3])
    ax_gap.axis("off")

    # Colorbar axis
    ax_cbar = fig.add_subplot(gs[0, 4])
    cbar = plt.colorbar(im, cax=ax_cbar)
    cbar.set_label("Scaled expression", fontsize=9)
    cbar.ax.tick_params(labelsize=7)

    # Cell type legend
    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, facecolor=ct_to_color[ct], label=ct)
        for ct in cell_types
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(0.80, 0.95),
        fontsize=6,
        title="Cell type",
        title_fontsize=7,
        frameon=True,
        fancybox=True,
    )

    # Title
    if title is None:
        title = f"Marker Genes ({label_column})"
    fig.suptitle(title, fontsize=14, fontweight="bold")

    if save:
        save_figure(fig, save)

    return fig


def plot_ontology_mapping(
    adata: ad.AnnData,
    source_label_column: str,
    ontology_name_column: str,
    ontology_id_column: str,
    mapping_table: Optional[pd.DataFrame] = None,
    title: Optional[str] = None,
    figsize: Tuple[float, float] = (14, 8),
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot ontology mapping table showing original labels mapped to Cell Ontology.

    Creates a table visualization with:
    - Original model labels
    - Mapped ontology names and CL IDs
    - Match tier (tier0=pattern, tier1=exact, etc.)
    - Matching score (actual score from matching, not hardcoded)
    - Cell counts

    Parameters
    ----------
    adata : AnnData
        Annotated data with ontology mapping columns.
    source_label_column : str
        Column with original model labels (e.g., "cell_type_predicted").
    ontology_name_column : str
        Column with mapped ontology names.
    ontology_id_column : str
        Column with Cell Ontology IDs (CL:XXXXXXX).
    mapping_table : pd.DataFrame, optional
        Pre-computed mapping table from OntologyMappingResult.table.
        If provided, uses this directly instead of building from adata.
        Should have columns: input_label, ontology_name, ontology_id,
        match_tier, score, n_cells.
    title : str, optional
        Plot title. Auto-generated if None.
    figsize : tuple
        Figure size.
    save : Path, optional
        Path to save figure.

    Returns
    -------
    Figure
        Matplotlib figure with ontology mapping table.

    Examples
    --------
    >>> # Using adata columns (scores read from {id_column}_score)
    >>> fig = plot_ontology_mapping(
    ...     adata,
    ...     source_label_column="cell_type_predicted",
    ...     ontology_name_column="cell_type_ontology_label",
    ...     ontology_id_column="cell_type_ontology_term_id",
    ... )

    >>> # Using pre-computed mapping table
    >>> from spatialcore.annotation import add_ontology_ids
    >>> adata, mappings, result = add_ontology_ids(adata, "cell_type_predicted", save_mapping="./")
    >>> fig = plot_ontology_mapping(
    ...     adata,
    ...     source_label_column="cell_type_predicted",
    ...     ontology_name_column="cell_type_ontology_label",
    ...     ontology_id_column="cell_type_ontology_term_id",
    ...     mapping_table=result.table,
    ... )
    """
    # If mapping_table is provided, use it directly
    if mapping_table is not None:
        summary = mapping_table.copy()
        # Rename columns to display names
        col_map = {
            "input_label": "Source Label",
            "ontology_name": "Ontology Name",
            "ontology_id": "CL ID",
            "match_tier": "Match Tier",
            "score": "Score",
            "n_cells": "Cells",
        }
        summary = summary.rename(columns=col_map)

        # Format score column
        summary["Score"] = summary["Score"].apply(
            lambda x: f"{x:.2f}" if pd.notna(x) and x > 0 else "-"
        )
    else:
        # Build from adata columns
        for col in [source_label_column, ontology_name_column, ontology_id_column]:
            if col not in adata.obs.columns:
                raise ValueError(
                    f"Column '{col}' not found in adata.obs. "
                    f"Available columns: {list(adata.obs.columns)}"
                )

        # Check for tier and score columns (CellxGene standard)
        tier_column = None
        score_column = None
        if ontology_id_column.endswith("_term_id"):
            tier_column = ontology_id_column.replace("_term_id", "_tier")
            score_column = ontology_id_column.replace("_term_id", "_score")
        has_tier_column = tier_column in adata.obs.columns if tier_column else False
        has_score_column = score_column in adata.obs.columns if score_column else False

        # Build mapping summary - convert to string to avoid categorical issues
        cols_to_use = [source_label_column, ontology_name_column, ontology_id_column]
        if has_tier_column:
            cols_to_use.append(tier_column)
        if has_score_column:
            cols_to_use.append(score_column)

        df = adata.obs[cols_to_use].copy()
        df[source_label_column] = df[source_label_column].astype(str)
        df[ontology_name_column] = df[ontology_name_column].astype(str)
        df[ontology_id_column] = df[ontology_id_column].astype(str)
        if has_tier_column:
            df[tier_column] = df[tier_column].astype(str)
        df["count"] = 1

        # Replace 'nan' strings with empty string
        df = df.replace("nan", "")

        # Aggregate by source label
        agg_dict = {
            ontology_name_column: "first",
            ontology_id_column: "first",
            "count": "sum",
        }
        if has_tier_column:
            agg_dict[tier_column] = "first"
        if has_score_column:
            agg_dict[score_column] = "first"

        summary = df.groupby(source_label_column).agg(agg_dict).reset_index()

        # Rename columns based on what we have
        if has_tier_column and has_score_column:
            summary.columns = ["Source Label", "Ontology Name", "CL ID", "Cells", "Match Tier", "Score"]
        elif has_tier_column:
            summary.columns = ["Source Label", "Ontology Name", "CL ID", "Cells", "Match Tier"]
        elif has_score_column:
            summary.columns = ["Source Label", "Ontology Name", "CL ID", "Cells", "Score"]
        else:
            summary.columns = ["Source Label", "Ontology Name", "CL ID", "Cells"]

        # Add tier column if missing
        if "Match Tier" not in summary.columns:
            def get_tier(row):
                cl_id = str(row["CL ID"]).strip()
                if not cl_id or cl_id == "" or cl_id == "-" or cl_id == "nan" or cl_id == "unknown" or cl_id == "skipped":
                    return "unmapped"
                return "tier0_pattern"  # Simplified - no tier info available

            summary["Match Tier"] = summary.apply(get_tier, axis=1)

        # Format score column - use actual scores if available
        if "Score" in summary.columns:
            summary["Score"] = summary["Score"].apply(
                lambda x: f"{float(x):.2f}" if pd.notna(x) and str(x) not in ["", "nan", "None"] and float(x) > 0 else "-"
            )
        else:
            # Fallback: estimate score based on tier (less accurate)
            tier_scores = {
                "tier0_pattern": "0.95",
                "tier1_exact": "1.00",
                "tier2_token": "0.75",
                "tier3_overlap": "0.60",
                "unmapped": "-",
                "skipped": "-",
            }
            summary["Score"] = summary["Match Tier"].map(lambda t: tier_scores.get(t, "-"))

    # Fill empty values with placeholder
    summary["CL ID"] = summary["CL ID"].replace("", "-")
    summary.loc[summary["Ontology Name"] == "", "Ontology Name"] = summary.loc[summary["Ontology Name"] == "", "Source Label"]

    # Sort by cell count descending
    summary = summary.sort_values("Cells", ascending=False).reset_index(drop=True)

    # Calculate stats
    n_labels = len(summary)
    n_mapped = (~summary["Match Tier"].isin(["unmapped", "skipped"])).sum()
    total_cells = summary["Cells"].sum()
    mapped_cells = summary[~summary["Match Tier"].isin(["unmapped", "skipped"])]["Cells"].sum()

    # Tier colors
    tier_colors = {
        "tier0_pattern": "#d4edda",  # Green
        "tier1_exact": "#cce5ff",    # Blue
        "tier2_token": "#fff3cd",    # Orange/Yellow
        "tier3_overlap": "#f8d7da",  # Red
        "unmapped": "#e9ecef",       # Gray
        "skipped": "#e9ecef",        # Gray (for Unassigned, Unknown, etc.)
    }

    # Calculate figure height based on number of rows
    n_rows = len(summary)
    fig_height = max(6, 2 + n_rows * 0.4)  # Dynamic height
    figsize = (figsize[0], fig_height)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    ax.axis("off")

    # Title and stats at top
    if title is None:
        title = "Source Label to Cell Ontology Mapping"

    fig.text(
        0.5, 0.97,
        f"Labels: {n_mapped}/{n_labels} mapped ({100*n_mapped/n_labels:.1f}%) | "
        f"Cells: {mapped_cells:,}/{total_cells:,} mapped ({100*mapped_cells/total_cells:.1f}%)",
        ha="center", fontsize=10, color="green"
    )
    fig.text(0.5, 0.93, title, ha="center", fontsize=14, fontweight="bold")

    # Reorder columns for display
    display_cols = ["Source Label", "Ontology Name", "CL ID", "Match Tier", "Score", "Cells"]
    summary = summary[display_cols]

    # Create table
    table = ax.table(
        cellText=summary.values,
        colLabels=summary.columns,
        cellLoc="left",
        loc="upper center",
        colColours=["#2c3e50"] * len(display_cols),
    )

    # Style table
    table.auto_set_font_size(False)
    table.set_fontsize(8)

    # Scale table to fit - adjust row height based on number of rows
    row_height = min(0.08, 0.8 / (n_rows + 1))
    table.scale(1.0, 1.5)

    # Set column widths
    col_widths = [0.22, 0.22, 0.12, 0.14, 0.08, 0.10]
    for j, width in enumerate(col_widths):
        for i in range(n_rows + 1):
            table[(i, j)].set_width(width)

    # Color header text white
    for j in range(len(display_cols)):
        table[(0, j)].get_text().set_color("white")
        table[(0, j)].get_text().set_fontweight("bold")

    # Color rows by tier
    for i in range(len(summary)):
        tier = summary.iloc[i]["Match Tier"]
        color = tier_colors.get(tier, "#ffffff")
        for j in range(len(display_cols)):
            table[(i + 1, j)].set_facecolor(color)

    # Add legend at bottom
    legend_text = (
        "Tier Colors:  Green = Pattern Match (tier0)  |  Blue = Exact Match (tier1)  |  "
        "Orange = Token Match (tier2)  |  Red = Overlap (tier3)  |  Gray = Unmapped"
    )
    fig.text(0.5, 0.02, legend_text, ha="center", fontsize=8, style="italic")

    if save:
        save_figure(fig, save)

    return fig


def generate_annotation_plots(
    adata: ad.AnnData,
    label_column: str = "cell_type",
    confidence_column: str = "cell_type_confidence",
    output_dir: Optional[Union[str, Path]] = None,
    prefix: str = "celltyping",
    confidence_threshold: float = 0.8,
    markers: Optional[Dict[str, List[str]]] = None,
    n_deg_genes: int = 10,
    spatial_key: str = "spatial",
    source_label_column: Optional[str] = None,
    ontology_name_column: Optional[str] = None,
    ontology_id_column: Optional[str] = None,
) -> Dict:
    """
    Generate all cell typing validation plots.

    Produces four standard validation outputs per the spec:
    1. Ontology mapping table - original labels to Cell Ontology
    2. 2D marker validation (GMM-3) - faceted by cell type
    3. Cell type confidence - spatial + jitter plot
    4. DEG heatmap - top genes per cell type

    This function should be called after annotation to validate results.
    Per spec: "Validation: 2D Multivariate QC" (Step H in workflow).

    Parameters
    ----------
    adata : AnnData
        Annotated data with cell type labels and confidence scores.
    label_column : str, default "cell_type"
        Column containing cell type labels (CellxGene standard).
        If ontology_name_column is provided and exists, plots use that
        column to display mapped ontology names.
    confidence_column : str, default "cell_type_confidence"
        Column containing confidence values (z-score transformed, CellxGene standard).
    output_dir : str or Path, optional
        Directory to save plots. If None, plots are returned but not saved.
    prefix : str, default "celltyping"
        Filename prefix for saved plots.
    confidence_threshold : float, default 0.8
        Threshold for confidence validation.
    markers : Dict[str, List[str]], optional
        Custom marker genes per cell type. If None, uses canonical markers
        from C:/SpatialCore/Data/markers/canonical_markers.json.
    n_deg_genes : int, default 10
        Number of top DEGs per cell type for heatmap.
    source_label_column : str, optional
        Original model label column (for ontology mapping table).
        If None, defaults to "cell_type_predicted" when present, else "cell_type".
    ontology_name_column : str, optional
        Ontology name column. If None, uses "cell_type_ontology_label" if present.
    ontology_id_column : str, optional
        Ontology ID column. If None, uses "cell_type_ontology_term_id" if present.
    spatial_key : str, default "spatial"
        Key for spatial coordinates in adata.obsm.

    Returns
    -------
    Dict
        Dictionary with keys:
        - "figures": dict of matplotlib Figure objects
        - "summary": validation summary DataFrame (from 2D validation)
        - "paths": dict of saved file paths (if output_dir provided)

    Examples
    --------
    >>> from spatialcore.plotting import generate_annotation_plots
    >>> results = generate_annotation_plots(
    ...     adata,
    ...     output_dir="./qc_plots",
    ...     prefix="lung_cancer",
    ... )
    >>> print(results["summary"])
    """
    output_dir = Path(output_dir) if output_dir else None
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)

    results = {"figures": {}, "summary": None, "paths": {}}

    # Infer column names if not provided (CellxGene standard only)
    if source_label_column is None:
        source_label_column = (
            "cell_type_predicted" if "cell_type_predicted" in adata.obs.columns else "cell_type"
        )

    if ontology_name_column is None and "cell_type_ontology_label" in adata.obs.columns:
        ontology_name_column = "cell_type_ontology_label"

    if ontology_id_column is None and "cell_type_ontology_term_id" in adata.obs.columns:
        ontology_id_column = "cell_type_ontology_term_id"

    # Prefer mapped ontology names for plots when available
    plot_label_column = (
        ontology_name_column
        if ontology_name_column and ontology_name_column in adata.obs.columns
        else label_column
    )

    # Guard against duplicate mapping columns
    mapping_cols = [c for c in [source_label_column, ontology_name_column, ontology_id_column] if c]
    if len(mapping_cols) != len(set(mapping_cols)):
        raise ValueError(
            "plot_ontology_mapping requires distinct source, ontology name, and ontology ID columns."
        )

    # 0. Ontology Mapping Table
    logger.info("Generating ontology mapping table...")
    if source_label_column and ontology_name_column and ontology_id_column:
        path_ontology = output_dir / f"{prefix}_ontology_mapping.png" if output_dir else None
        fig_ontology = plot_ontology_mapping(
            adata,
            source_label_column=source_label_column,
            ontology_name_column=ontology_name_column,
            ontology_id_column=ontology_id_column,
            save=path_ontology,
        )
        results["figures"]["ontology_mapping"] = fig_ontology
        results["paths"]["ontology_mapping"] = path_ontology
        logger.info("  Ontology mapping table generated")
    else:
        logger.info(f"  Skipping ontology table - columns not found")
        logger.info(f"    source_label_column: {source_label_column}")
        logger.info(f"    ontology_name_column: {ontology_name_column}")
        logger.info(f"    ontology_id_column: {ontology_id_column}")

    # 1. 2D Marker Validation (GMM-3)
    logger.info("Generating 2D marker validation plot...")
    path_2d = output_dir / f"{prefix}_2d_validation.png" if output_dir else None
    fig_2d, summary = plot_2d_validation(
        adata,
        label_column=plot_label_column,
        confidence_column=confidence_column,
        markers=markers,
        confidence_threshold=confidence_threshold,
        n_components=3,  # GMM-3 per spec
        save=path_2d,
    )
    results["figures"]["2d_validation"] = fig_2d
    results["summary"] = summary
    results["paths"]["2d_validation"] = path_2d
    logger.info(f"  2D validation: {len(summary)} cell types analyzed")

    # 2. Cell Type Confidence
    logger.info("Generating confidence plot...")
    path_conf = output_dir / f"{prefix}_confidence.png" if output_dir else None
    fig_conf = plot_celltype_confidence(
        adata,
        label_column=plot_label_column,
        confidence_column=confidence_column,
        spatial_key=spatial_key,
        threshold=confidence_threshold,
        save=path_conf,
    )
    results["figures"]["confidence"] = fig_conf
    results["paths"]["confidence"] = path_conf
    logger.info("  Confidence plot generated")

    # 3. DEG Heatmap
    logger.info("Generating DEG heatmap...")
    path_deg = output_dir / f"{prefix}_deg_heatmap.png" if output_dir else None
    try:
        fig_deg = plot_deg_heatmap(
            adata,
            label_column=plot_label_column,
            n_genes=n_deg_genes,
            save=path_deg,
        )
    except ValueError as exc:
        if "Need at least 2 cell types" in str(exc):
            logger.warning("Skipping DEG heatmap: %s", exc)
            fig_deg = None
            path_deg = None
        else:
            raise
    results["figures"]["deg_heatmap"] = fig_deg
    results["paths"]["deg_heatmap"] = path_deg
    if fig_deg is not None:
        logger.info("  DEG heatmap generated")

    if output_dir:
        logger.info(f"Annotation plots saved to: {output_dir}")

    return results
