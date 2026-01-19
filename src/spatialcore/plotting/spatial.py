"""
Spatial visualization of cell types and confidence.

This module provides functions for visualizing cell annotations
on spatial coordinates.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import anndata as ad

from spatialcore.core.logging import get_logger
from spatialcore.plotting.utils import (
    generate_celltype_palette,
    setup_figure,
    save_figure,
    format_axis_labels,
)

logger = get_logger(__name__)


def plot_spatial_celltype(
    adata: ad.AnnData,
    label_column: str,
    spatial_key: str = "spatial",
    colors: Optional[Dict[str, str]] = None,
    point_size: float = 1.0,
    alpha: float = 0.8,
    figsize: tuple = (10, 10),
    dark_background: bool = True,
    legend_loc: str = "right margin",
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot cell types on spatial coordinates.

    Parameters
    ----------
    adata : AnnData
        Annotated data with spatial coordinates.
    label_column : str
        Column in adata.obs containing cell type labels.
    spatial_key : str, default "spatial"
        Key in adata.obsm for spatial coordinates.
    colors : Dict[str, str], optional
        Color mapping for cell types.
    point_size : float, default 1.0
        Size of points.
    alpha : float, default 0.8
        Point transparency.
    figsize : tuple, default (10, 10)
        Figure size.
    dark_background : bool, default True
        Use dark background (better for spatial data).
    legend_loc : str, default "right margin"
        Legend location: "right margin", "on data", "none".
    xlim : Tuple[float, float], optional
        X-axis limits.
    ylim : Tuple[float, float], optional
        Y-axis limits.
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
    >>> from spatialcore.plotting.spatial import plot_spatial_celltype
    >>> fig = plot_spatial_celltype(
    ...     adata,
    ...     label_column="cell_type",
    ...     point_size=0.5,
    ... )
    """
    if label_column not in adata.obs.columns:
        raise ValueError(f"Label column '{label_column}' not found.")

    if spatial_key not in adata.obsm:
        raise ValueError(
            f"Spatial key '{spatial_key}' not found. "
            f"Available: {list(adata.obsm.keys())}"
        )

    coords = adata.obsm[spatial_key]
    cell_types = adata.obs[label_column].astype(str)
    unique_types = sorted(cell_types.unique())

    if colors is None:
        colors = generate_celltype_palette(unique_types)

    # Adjust figure size for legend
    if legend_loc == "right margin":
        figsize = (figsize[0] + 3, figsize[1])

    if dark_background:
        plt.style.use("dark_background")

    fig, ax = setup_figure(figsize=figsize)

    # Plot each cell type
    for cell_type in unique_types:
        mask = cell_types == cell_type
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            c=colors.get(cell_type, "#888888"),
            label=cell_type,
            s=point_size,
            alpha=alpha,
            rasterized=True,
        )

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect("equal")

    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    # Legend
    if legend_loc == "right margin":
        ax.legend(
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            markerscale=5,
            frameon=False,
        )
    elif legend_loc == "on data":
        ax.legend(loc="best", markerscale=5)

    if title is None:
        title = f"Spatial Cell Types ({label_column})"
    ax.set_title(title, color="white" if dark_background else "black")

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    # Reset style
    if dark_background:
        plt.style.use("default")

    return fig


def plot_spatial_confidence(
    adata: ad.AnnData,
    confidence_column: str,
    spatial_key: str = "spatial",
    cmap: str = "viridis",
    point_size: float = 1.0,
    alpha: float = 0.8,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    figsize: tuple = (10, 10),
    dark_background: bool = True,
    colorbar: bool = True,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot confidence scores on spatial coordinates.

    Parameters
    ----------
    adata : AnnData
        Annotated data with spatial coordinates and confidence.
    confidence_column : str
        Column in adata.obs containing confidence values.
    spatial_key : str, default "spatial"
        Key in adata.obsm for spatial coordinates.
    cmap : str, default "viridis"
        Colormap name.
    point_size : float, default 1.0
        Size of points.
    alpha : float, default 0.8
        Point transparency.
    vmin : float, optional
        Minimum value for color scale.
    vmax : float, optional
        Maximum value for color scale.
    figsize : tuple, default (10, 10)
        Figure size.
    dark_background : bool, default True
        Use dark background.
    colorbar : bool, default True
        Show colorbar.
    xlim : Tuple[float, float], optional
        X-axis limits.
    ylim : Tuple[float, float], optional
        Y-axis limits.
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
    >>> from spatialcore.plotting.spatial import plot_spatial_confidence
    >>> fig = plot_spatial_confidence(
    ...     adata,
    ...     confidence_column="confidence",
    ...     cmap="plasma",
    ... )
    """
    if confidence_column not in adata.obs.columns:
        raise ValueError(f"Confidence column '{confidence_column}' not found.")

    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial key '{spatial_key}' not found.")

    coords = adata.obsm[spatial_key]
    confidence = adata.obs[confidence_column].values

    if dark_background:
        plt.style.use("dark_background")

    fig, ax = setup_figure(figsize=figsize)

    scatter = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=confidence,
        cmap=cmap,
        s=point_size,
        alpha=alpha,
        vmin=vmin,
        vmax=vmax,
        rasterized=True,
    )

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect("equal")

    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    if colorbar:
        cbar = plt.colorbar(scatter, ax=ax, shrink=0.6)
        cbar.set_label("Confidence", fontsize=10)

    if title is None:
        title = f"Spatial Confidence ({confidence_column})"
    ax.set_title(title, color="white" if dark_background else "black")

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    if dark_background:
        plt.style.use("default")

    return fig


def plot_spatial_gene(
    adata: ad.AnnData,
    gene: str,
    spatial_key: str = "spatial",
    layer: Optional[str] = None,
    cmap: str = "Reds",
    point_size: float = 1.0,
    alpha: float = 0.8,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    figsize: tuple = (10, 10),
    dark_background: bool = True,
    colorbar: bool = True,
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot gene expression on spatial coordinates.

    Parameters
    ----------
    adata : AnnData
        Annotated data with spatial coordinates.
    gene : str
        Gene name to plot.
    spatial_key : str, default "spatial"
        Key in adata.obsm for spatial coordinates.
    layer : str, optional
        Layer to use. If None, uses adata.X.
    cmap : str, default "Reds"
        Colormap name.
    point_size : float, default 1.0
        Size of points.
    alpha : float, default 0.8
        Point transparency.
    vmin : float, optional
        Minimum value for color scale.
    vmax : float, optional
        Maximum value for color scale.
    figsize : tuple, default (10, 10)
        Figure size.
    dark_background : bool, default True
        Use dark background.
    colorbar : bool, default True
        Show colorbar.
    title : str, optional
        Plot title.
    save : str or Path, optional
        Path to save figure.

    Returns
    -------
    Figure
        Matplotlib figure.
    """
    if gene not in adata.var_names:
        raise ValueError(
            f"Gene '{gene}' not found. "
            f"Available genes: {list(adata.var_names[:10])}..."
        )

    if spatial_key not in adata.obsm:
        raise ValueError(f"Spatial key '{spatial_key}' not found.")

    coords = adata.obsm[spatial_key]

    # Get expression
    gene_idx = adata.var_names.get_loc(gene)
    if layer is not None:
        expression = adata.layers[layer][:, gene_idx]
    else:
        X = adata.X
        if hasattr(X, "toarray"):
            expression = X[:, gene_idx].toarray().flatten()
        else:
            expression = X[:, gene_idx].flatten()

    if dark_background:
        plt.style.use("dark_background")

    fig, ax = setup_figure(figsize=figsize)

    scatter = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=expression,
        cmap=cmap,
        s=point_size,
        alpha=alpha,
        vmin=vmin,
        vmax=vmax,
        rasterized=True,
    )

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_aspect("equal")

    if colorbar:
        cbar = plt.colorbar(scatter, ax=ax, shrink=0.6)
        cbar.set_label("Expression", fontsize=10)

    if title is None:
        title = f"{gene} Expression"
    ax.set_title(title, color="white" if dark_background else "black")

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    if dark_background:
        plt.style.use("default")

    return fig


def plot_spatial_multi_gene(
    adata: ad.AnnData,
    genes: List[str],
    spatial_key: str = "spatial",
    layer: Optional[str] = None,
    cmap: str = "Reds",
    point_size: float = 0.5,
    ncols: int = 3,
    figsize_per_panel: Tuple[float, float] = (4, 4),
    dark_background: bool = True,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot multiple genes on spatial coordinates.

    Parameters
    ----------
    adata : AnnData
        Annotated data with spatial coordinates.
    genes : List[str]
        List of gene names to plot.
    spatial_key : str, default "spatial"
        Key in adata.obsm for spatial coordinates.
    layer : str, optional
        Layer to use.
    cmap : str, default "Reds"
        Colormap name.
    point_size : float, default 0.5
        Size of points.
    ncols : int, default 3
        Number of columns in subplot grid.
    figsize_per_panel : Tuple[float, float], default (4, 4)
        Size per panel.
    dark_background : bool, default True
        Use dark background.
    save : str or Path, optional
        Path to save figure.

    Returns
    -------
    Figure
        Matplotlib figure.
    """
    # Filter to available genes
    available_genes = [g for g in genes if g in adata.var_names]
    missing = set(genes) - set(available_genes)
    if missing:
        logger.warning(f"Genes not found: {missing}")

    if not available_genes:
        raise ValueError("No valid genes found.")

    n_genes = len(available_genes)
    nrows = int(np.ceil(n_genes / ncols))

    figsize = (figsize_per_panel[0] * ncols, figsize_per_panel[1] * nrows)

    if dark_background:
        plt.style.use("dark_background")

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    axes = np.atleast_2d(axes).flatten()

    coords = adata.obsm[spatial_key]

    for i, gene in enumerate(available_genes):
        ax = axes[i]

        gene_idx = adata.var_names.get_loc(gene)
        if layer is not None:
            expression = adata.layers[layer][:, gene_idx]
        else:
            X = adata.X
            if hasattr(X, "toarray"):
                expression = X[:, gene_idx].toarray().flatten()
            else:
                expression = X[:, gene_idx].flatten()

        ax.scatter(
            coords[:, 0],
            coords[:, 1],
            c=expression,
            cmap=cmap,
            s=point_size,
            rasterized=True,
        )
        ax.set_aspect("equal")
        ax.set_title(gene, fontsize=10)
        ax.set_xticks([])
        ax.set_yticks([])

    # Hide empty panels
    for i in range(n_genes, len(axes)):
        axes[i].set_visible(False)

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    if dark_background:
        plt.style.use("default")

    return fig
