"""Visualization functions for spatial NMF results.

This module provides plotting functions for interpreting spaNMF
components, including spatial distributions, gene loadings, and
component correlations.
"""

from pathlib import Path
from typing import List, Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.figure import Figure

import anndata as ad

from spatialcore.core.logging import get_logger
from spatialcore.plotting.utils import save_figure

logger = get_logger(__name__)


def plot_top_genes(
    adata: ad.AnnData,
    n_genes: int = 15,
    components: Optional[List[int]] = None,
    key: str = "spanmf",
    ncols: int = 4,
    figsize_per_panel: Tuple[float, float] = (3, 4),
    color: str = "#3498db",
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot top genes for each NMF component.

    Creates a grid of horizontal bar plots showing the highest-loading
    genes for each component.

    Parameters
    ----------
    adata
        AnnData with NMF results from run_spanmf().
    n_genes
        Number of top genes to show per component. Default: 15.
    components
        Which components to plot (0-indexed). If None, plots all.
    key
        Key prefix used in run_spanmf(). Default: "spanmf".
    ncols
        Number of columns in subplot grid. Default: 4.
    figsize_per_panel
        Size per panel (width, height). Default: (3, 4).
    color
        Bar color. Default: "#3498db" (blue).
    save
        Path to save figure. If None, does not save.

    Returns
    -------
    Figure
        Matplotlib figure.

    Raises
    ------
    KeyError
        If NMF results not found in adata.

    Examples
    --------
    >>> from spatialcore.plotting.nmf import plot_top_genes
    >>> fig = plot_top_genes(adata, n_genes=12, ncols=5)
    """
    # Validate inputs
    loadings_key = f"{key}_loadings"
    if loadings_key not in adata.varm:
        raise KeyError(
            f"Gene loadings not found in adata.varm['{loadings_key}']. "
            "Run run_spanmf() first."
        )

    loadings = adata.varm[loadings_key]  # (n_genes, n_components)
    gene_names = adata.var_names
    n_components = loadings.shape[1]

    # Determine which components to plot
    if components is None:
        components = list(range(n_components))
    else:
        for c in components:
            if c < 0 or c >= n_components:
                raise ValueError(
                    f"Component {c} out of range [0, {n_components - 1}]"
                )

    n_plots = len(components)
    nrows = int(np.ceil(n_plots / ncols))

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows),
    )
    axes = np.atleast_2d(axes).flatten()

    for idx, comp in enumerate(components):
        ax = axes[idx]

        # Get top genes for this component
        gene_loadings = loadings[:, comp]
        top_idx = np.argsort(gene_loadings)[-n_genes:][::-1]
        top_genes = gene_names[top_idx]
        top_values = gene_loadings[top_idx]

        # Plot horizontal bar chart (solid color)
        ax.barh(range(n_genes), top_values, color=color)
        ax.set_yticks(range(n_genes))
        ax.set_yticklabels(top_genes, fontsize=8)
        ax.invert_yaxis()
        ax.set_xlabel("Loading", fontsize=9)
        ax.set_title(f"Component {comp + 1}", fontsize=10, fontweight="bold")

    # Hide empty panels
    for idx in range(n_plots, len(axes)):
        axes[idx].set_visible(False)

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    return fig


def plot_gene_loading(
    adata: ad.AnnData,
    gene: str,
    n_genes: int = 10,
    key: str = "spanmf",
    figsize: Tuple[float, float] = (6, 5),
    highlight_color: str = "#e74c3c",
    other_color: str = "#3498db",
    save: Optional[Union[str, Path]] = None,
) -> Tuple[Figure, int]:
    """
    Show which component a gene loads highest on.

    Displays bar plot of top genes in that component,
    highlighting the query gene.

    Parameters
    ----------
    adata
        AnnData with NMF results.
    gene
        Gene name to query.
    n_genes
        Number of genes to show in bar plot. Default: 10.
    key
        Key prefix from run_spanmf(). Default: "spanmf".
    figsize
        Figure size. Default: (6, 5).
    highlight_color
        Color for the query gene bar. Default: red.
    other_color
        Color for other gene bars. Default: blue.
    save
        Path to save figure.

    Returns
    -------
    fig
        Matplotlib figure.
    component
        0-indexed component where the gene has highest loading.

    Raises
    ------
    ValueError
        If gene not found in dataset.
    KeyError
        If NMF results not found.

    Examples
    --------
    >>> fig, comp = plot_gene_loading(adata, gene="CD8A")
    >>> print(f"CD8A loads highest on component {comp + 1}")
    """
    if gene not in adata.var_names:
        raise ValueError(
            f"Gene '{gene}' not found. "
            f"Available genes: {list(adata.var_names[:5])}..."
        )

    loadings_key = f"{key}_loadings"
    if loadings_key not in adata.varm:
        raise KeyError(
            f"Gene loadings not found in adata.varm['{loadings_key}']. "
            "Run run_spanmf() first."
        )

    loadings = adata.varm[loadings_key]
    gene_names = adata.var_names

    # Find gene's index and its max component
    gene_idx = adata.var_names.get_loc(gene)
    gene_loadings = loadings[gene_idx, :]
    max_component = int(np.argmax(gene_loadings))

    # Get top genes in that component
    component_loadings = loadings[:, max_component]
    top_idx = np.argsort(component_loadings)[-n_genes:][::-1]

    # Ensure query gene is included
    if gene_idx not in top_idx:
        top_idx = np.append(top_idx, gene_idx)

    top_genes = gene_names[top_idx]
    top_values = component_loadings[top_idx]

    # Sort by loading value
    sort_order = np.argsort(top_values)[::-1]
    top_genes = top_genes[sort_order]
    top_values = top_values[sort_order]

    # Create colors
    colors = [highlight_color if g == gene else other_color for g in top_genes]

    # Plot
    fig, ax = plt.subplots(figsize=figsize)

    ax.barh(range(len(top_genes)), top_values, color=colors)
    ax.set_yticks(range(len(top_genes)))
    ax.set_yticklabels(top_genes, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("Loading", fontsize=10)
    ax.set_title(
        f"Component {max_component + 1} (highest for {gene})",
        fontsize=11,
        fontweight="bold",
    )

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    return fig, max_component


def plot_component_spatial(
    adata: ad.AnnData,
    components: Union[int, List[int], str] = "all",
    key: str = "spanmf",
    spatial_key: str = "spatial",
    ncols: int = 4,
    spot_size: float = 1.0,
    cmap: str = "viridis",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    figsize_per_panel: Tuple[float, float] = (4, 4),
    flip_xy: bool = True,
    dark_background: bool = False,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot NMF component loadings on spatial coordinates.

    Creates a grid of spatial plots, one per component, showing
    how each NMF factor is distributed across the tissue.

    Parameters
    ----------
    adata
        AnnData with NMF results and spatial coordinates.
    components
        Which components to plot:

        - int: single component (0-indexed)
        - List[int]: specific components
        - "all": all components (default)
    key
        Key prefix from run_spanmf(). Default: "spanmf".
    spatial_key
        Key in adata.obsm for spatial coordinates. Default: "spatial".
    ncols
        Number of columns in grid. Default: 4.
    spot_size
        Point size for scatter. Default: 1.0.
    cmap
        Colormap for component values. Default: "viridis".
    vmin, vmax
        Color scale limits. If None, uses per-component min/max.
    figsize_per_panel
        Size per panel. Default: (4, 4).
    flip_xy
        If True, swap X and Y coordinates. Default: True.
    dark_background
        Use dark background style. Default: False.
    save
        Path to save figure.

    Returns
    -------
    Figure
        Matplotlib figure.

    Raises
    ------
    ValueError
        If spatial coordinates or NMF results not found.

    Examples
    --------
    >>> from spatialcore.plotting.nmf import plot_component_spatial
    >>> plot_component_spatial(adata, components=[0, 2, 5])
    """
    # Validate inputs
    cell_loadings_key = f"X_{key}"
    if cell_loadings_key not in adata.obsm:
        raise KeyError(
            f"Cell loadings not found in adata.obsm['{cell_loadings_key}']. "
            "Run run_spanmf() first."
        )

    if spatial_key not in adata.obsm:
        raise ValueError(
            f"Spatial coordinates not found in adata.obsm['{spatial_key}']."
        )

    W = adata.obsm[cell_loadings_key]  # (n_cells, n_components)
    coords = adata.obsm[spatial_key]
    n_components = W.shape[1]

    # Parse components argument
    if components == "all":
        components = list(range(n_components))
    elif isinstance(components, int):
        components = [components]

    for c in components:
        if c < 0 or c >= n_components:
            raise ValueError(
                f"Component {c} out of range [0, {n_components - 1}]"
            )

    n_plots = len(components)
    nrows = int(np.ceil(n_plots / ncols))

    if dark_background:
        plt.style.use("dark_background")

    fig, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows),
    )
    axes = np.atleast_2d(axes).flatten()

    for idx, comp in enumerate(components):
        ax = axes[idx]
        values = W[:, comp]

        # Optionally flip XY coordinates
        x_coord = coords[:, 1] if flip_xy else coords[:, 0]
        y_coord = coords[:, 0] if flip_xy else coords[:, 1]

        scatter = ax.scatter(
            x_coord,
            y_coord,
            c=values,
            cmap=cmap,
            s=spot_size,
            vmin=vmin,
            vmax=vmax,
            rasterized=True,
        )

        ax.set_aspect("equal")
        ax.set_title(f"Component {comp + 1}", fontsize=10, fontweight="bold")
        ax.set_xticks([])
        ax.set_yticks([])

        plt.colorbar(scatter, ax=ax, shrink=0.6)

    # Hide empty panels
    for idx in range(n_plots, len(axes)):
        axes[idx].set_visible(False)

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    if dark_background:
        plt.style.use("default")

    return fig


def plot_component_correlation(
    adata: ad.AnnData,
    key: str = "spanmf",
    method: Literal["pearson", "spearman"] = "pearson",
    cluster: bool = True,
    annot: bool = True,
    cmap: str = "RdBu_r",
    figsize: Tuple[float, float] = (8, 7),
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot correlation matrix between NMF components.

    Uses cell loadings to compute pairwise correlations between
    components. Helps identify redundant or anti-correlated factors.

    Parameters
    ----------
    adata
        AnnData with NMF results.
    key
        Key prefix from run_spanmf(). Default: "spanmf".
    method
        Correlation method: "pearson" or "spearman". Default: "pearson".
    cluster
        If True, hierarchically cluster rows/columns. Default: True.
    annot
        If True, show correlation values in cells. Default: True.
    cmap
        Colormap (diverging recommended). Default: "RdBu_r".
    figsize
        Figure size. Default: (8, 7).
    save
        Path to save figure.

    Returns
    -------
    Figure
        Matplotlib figure.

    Raises
    ------
    KeyError
        If NMF results not found.

    Notes
    -----
    High correlation (|r| > 0.7) between components suggests
    n_components may be too high. Anti-correlated components
    often represent opposing biological programs.

    Examples
    --------
    >>> fig = plot_component_correlation(adata, method="spearman")
    """
    cell_loadings_key = f"X_{key}"
    if cell_loadings_key not in adata.obsm:
        raise KeyError(
            f"Cell loadings not found in adata.obsm['{cell_loadings_key}']. "
            "Run run_spanmf() first."
        )

    W = adata.obsm[cell_loadings_key]
    n_components = W.shape[1]

    # Compute correlation matrix
    df = pd.DataFrame(
        W,
        columns=[f"C{i + 1}" for i in range(n_components)],
    )
    corr = df.corr(method=method)

    # Plot
    fig, ax = plt.subplots(figsize=figsize)

    if cluster:
        # Use seaborn clustermap then extract order
        g = sns.clustermap(
            corr,
            cmap=cmap,
            vmin=-1,
            vmax=1,
            annot=annot,
            fmt=".2f",
            figsize=figsize,
        )
        plt.close(g.fig)  # Close clustermap figure

        # Recreate with clustering order
        order = g.dendrogram_row.reordered_ind
        corr_ordered = corr.iloc[order, order]

        sns.heatmap(
            corr_ordered,
            cmap=cmap,
            vmin=-1,
            vmax=1,
            annot=annot,
            fmt=".2f",
            ax=ax,
            square=True,
        )
    else:
        sns.heatmap(
            corr,
            cmap=cmap,
            vmin=-1,
            vmax=1,
            annot=annot,
            fmt=".2f",
            ax=ax,
            square=True,
        )

    ax.set_title(
        f"Component Correlation ({method.capitalize()})",
        fontsize=12,
        fontweight="bold",
    )

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    return fig
