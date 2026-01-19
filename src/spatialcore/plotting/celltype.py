"""
Cell type distribution and UMAP visualization.

This module provides functions for visualizing cell type distributions
and embeddings.
"""

from pathlib import Path
from typing import Dict, List, Optional, Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
import anndata as ad

from spatialcore.core.logging import get_logger
from spatialcore.plotting.utils import (
    generate_celltype_palette,
    setup_figure,
    save_figure,
    despine,
    format_axis_labels,
)

logger = get_logger(__name__)


def plot_celltype_distribution(
    adata: ad.AnnData,
    label_column: str,
    colors: Optional[Dict[str, str]] = None,
    horizontal: bool = False,
    top_n: Optional[int] = None,
    figsize: Optional[tuple] = None,
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot cell type distribution as bar chart.

    Parameters
    ----------
    adata : AnnData
        Annotated data with cell type labels.
    label_column : str
        Column in adata.obs containing cell type labels.
    colors : Dict[str, str], optional
        Color mapping for cell types. If None, auto-generated.
    horizontal : bool, default False
        Plot horizontal bars (easier to read with many types).
    top_n : int, optional
        Only show top N most frequent cell types.
    figsize : tuple, optional
        Figure size. Auto-calculated if None.
    title : str, optional
        Plot title. Default: "Cell Type Distribution".
    save : str or Path, optional
        Path to save figure (without extension).

    Returns
    -------
    Figure
        Matplotlib figure.

    Examples
    --------
    >>> from spatialcore.plotting.celltype import plot_celltype_distribution
    >>> fig = plot_celltype_distribution(
    ...     adata,
    ...     label_column="cell_type",
    ...     horizontal=True,
    ...     top_n=20,
    ... )
    """
    if label_column not in adata.obs.columns:
        raise ValueError(
            f"Label column '{label_column}' not found. "
            f"Available: {list(adata.obs.columns)}"
        )

    # Get counts
    counts = adata.obs[label_column].value_counts()

    if top_n is not None:
        counts = counts.head(top_n)

    cell_types = counts.index.tolist()
    values = counts.values

    # Generate colors
    if colors is None:
        colors = generate_celltype_palette(cell_types)

    bar_colors = [colors.get(ct, "#888888") for ct in cell_types]

    # Calculate figure size
    n_types = len(cell_types)
    if figsize is None:
        if horizontal:
            figsize = (8, max(4, 0.3 * n_types))
        else:
            figsize = (max(8, 0.5 * n_types), 6)

    fig, ax = setup_figure(figsize=figsize)

    if horizontal:
        y_pos = np.arange(len(cell_types))
        ax.barh(y_pos, values, color=bar_colors)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(cell_types)
        ax.invert_yaxis()  # Top to bottom
        format_axis_labels(ax, xlabel="Number of Cells")
    else:
        x_pos = np.arange(len(cell_types))
        ax.bar(x_pos, values, color=bar_colors)
        ax.set_xticks(x_pos)
        ax.set_xticklabels(cell_types, rotation=45, ha="right")
        format_axis_labels(ax, ylabel="Number of Cells")

    despine(ax)

    if title is None:
        title = "Cell Type Distribution"
    ax.set_title(title)

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    return fig


def plot_celltype_pie(
    adata: ad.AnnData,
    label_column: str,
    colors: Optional[Dict[str, str]] = None,
    min_pct: float = 2.0,
    other_label: str = "Other",
    figsize: tuple = (8, 8),
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot cell type distribution as pie chart.

    Parameters
    ----------
    adata : AnnData
        Annotated data with cell type labels.
    label_column : str
        Column in adata.obs containing cell type labels.
    colors : Dict[str, str], optional
        Color mapping for cell types.
    min_pct : float, default 2.0
        Minimum percentage to show as separate slice.
        Smaller types are grouped into "Other".
    other_label : str, default "Other"
        Label for grouped small cell types.
    figsize : tuple, default (8, 8)
        Figure size.
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

    counts = adata.obs[label_column].value_counts()
    pcts = 100 * counts / counts.sum()

    # Group small types
    main_types = pcts[pcts >= min_pct]
    other_pct = pcts[pcts < min_pct].sum()

    if other_pct > 0:
        main_types[other_label] = other_pct

    cell_types = main_types.index.tolist()
    values = main_types.values

    # Generate colors
    if colors is None:
        colors = generate_celltype_palette(cell_types)
        colors[other_label] = "#888888"

    pie_colors = [colors.get(ct, "#888888") for ct in cell_types]

    fig, ax = setup_figure(figsize=figsize)

    ax.pie(
        values,
        labels=cell_types,
        colors=pie_colors,
        autopct="%1.1f%%",
        pctdistance=0.8,
    )

    if title is None:
        title = "Cell Type Distribution"
    ax.set_title(title)

    if save:
        save_figure(fig, save)

    return fig


def plot_celltype_umap(
    adata: ad.AnnData,
    label_column: str,
    colors: Optional[Dict[str, str]] = None,
    obsm_key: str = "X_umap",
    point_size: float = 1.0,
    alpha: float = 0.7,
    legend_loc: str = "right margin",
    figsize: tuple = (10, 8),
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot cell types on UMAP embedding.

    Parameters
    ----------
    adata : AnnData
        Annotated data with UMAP coordinates.
    label_column : str
        Column in adata.obs containing cell type labels.
    colors : Dict[str, str], optional
        Color mapping for cell types.
    obsm_key : str, default "X_umap"
        Key in adata.obsm for embedding coordinates.
    point_size : float, default 1.0
        Size of points.
    alpha : float, default 0.7
        Point transparency.
    legend_loc : str, default "right margin"
        Legend location: "right margin", "on data", "none".
    figsize : tuple, default (10, 8)
        Figure size.
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
    >>> from spatialcore.plotting.celltype import plot_celltype_umap
    >>> fig = plot_celltype_umap(
    ...     adata,
    ...     label_column="cell_type",
    ...     point_size=0.5,
    ... )
    """
    if label_column not in adata.obs.columns:
        raise ValueError(f"Label column '{label_column}' not found.")

    if obsm_key not in adata.obsm:
        raise ValueError(
            f"Embedding '{obsm_key}' not found. "
            f"Available: {list(adata.obsm.keys())}"
        )

    coords = adata.obsm[obsm_key]
    cell_types = adata.obs[label_column].astype(str)
    unique_types = sorted(cell_types.unique())

    # Generate colors
    if colors is None:
        colors = generate_celltype_palette(unique_types)

    # Adjust figure size for legend
    if legend_loc == "right margin":
        figsize = (figsize[0] + 3, figsize[1])

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
            rasterized=True,  # Faster for large datasets
        )

    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.set_aspect("equal")

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
    # else: no legend

    if title is None:
        title = f"Cell Types ({label_column})"
    ax.set_title(title)

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    return fig
