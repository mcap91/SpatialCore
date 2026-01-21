"""
Confidence score visualization.

This module provides functions for visualizing prediction confidence
distributions and comparing confidence across cell types.
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
    setup_multi_figure,
    save_figure,
    despine,
    format_axis_labels,
)

logger = get_logger(__name__)


def plot_confidence_histogram(
    adata: ad.AnnData,
    confidence_column: str,
    bins: int = 50,
    threshold: Optional[float] = None,
    threshold_color: str = "#FF0000",
    figsize: tuple = (8, 5),
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot histogram of confidence scores.

    Parameters
    ----------
    adata : AnnData
        Annotated data with confidence values.
    confidence_column : str
        Column in adata.obs containing confidence values.
    bins : int, default 50
        Number of histogram bins.
    threshold : float, optional
        Confidence threshold to highlight with vertical line.
    threshold_color : str, default "#FF0000"
        Color for threshold line.
    figsize : tuple, default (8, 5)
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
    >>> from spatialcore.plotting.confidence import plot_confidence_histogram
    >>> fig = plot_confidence_histogram(
    ...     adata,
    ...     confidence_column="celltypist_confidence",
    ...     threshold=0.5,
    ... )
    """
    if confidence_column not in adata.obs.columns:
        raise ValueError(
            f"Confidence column '{confidence_column}' not found. "
            f"Available: {list(adata.obs.columns)}"
        )

    values = adata.obs[confidence_column].values

    fig, ax = setup_figure(figsize=figsize)

    ax.hist(values, bins=bins, color="#3784FE", edgecolor="white", alpha=0.8)

    if threshold is not None:
        ax.axvline(threshold, color=threshold_color, linestyle="--", linewidth=2)
        below = (values < threshold).sum()
        pct = 100 * below / len(values)
        ax.text(
            threshold + 0.02,
            ax.get_ylim()[1] * 0.9,
            f"{pct:.1f}% below\nthreshold",
            color=threshold_color,
            fontsize=10,
        )

    format_axis_labels(
        ax,
        xlabel="Confidence Score",
        ylabel="Number of Cells",
    )
    despine(ax)

    if title is None:
        title = f"Confidence Distribution ({confidence_column})"
    ax.set_title(title)

    # Add statistics
    stats_text = (
        f"Mean: {np.mean(values):.3f}\n"
        f"Median: {np.median(values):.3f}\n"
        f"Std: {np.std(values):.3f}"
    )
    ax.text(
        0.02,
        0.98,
        stats_text,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
    )

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    return fig


def plot_confidence_by_celltype(
    adata: ad.AnnData,
    label_column: str,
    confidence_column: str,
    colors: Optional[Dict[str, str]] = None,
    top_n: Optional[int] = 20,
    figsize: Optional[tuple] = None,
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot confidence distribution per cell type as box plots.

    Parameters
    ----------
    adata : AnnData
        Annotated data with cell type labels and confidence.
    label_column : str
        Column in adata.obs containing cell type labels.
    confidence_column : str
        Column in adata.obs containing confidence values.
    colors : Dict[str, str], optional
        Color mapping for cell types.
    top_n : int, optional, default 20
        Only show top N most frequent cell types.
    figsize : tuple, optional
        Figure size. Auto-calculated if None.
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
    >>> from spatialcore.plotting.confidence import plot_confidence_by_celltype
    >>> fig = plot_confidence_by_celltype(
    ...     adata,
    ...     label_column="cell_type",
    ...     confidence_column="confidence",
    ... )
    """
    if label_column not in adata.obs.columns:
        raise ValueError(f"Label column '{label_column}' not found.")
    if confidence_column not in adata.obs.columns:
        raise ValueError(f"Confidence column '{confidence_column}' not found.")

    # Get data
    df = pd.DataFrame({
        "cell_type": adata.obs[label_column].values,
        "confidence": adata.obs[confidence_column].values,
    })

    # Get top N types by count
    type_order = df["cell_type"].value_counts().index.tolist()
    if top_n is not None:
        type_order = type_order[:top_n]
        df = df[df["cell_type"].isin(type_order)]

    n_types = len(type_order)

    # Generate colors
    if colors is None:
        colors = generate_celltype_palette(type_order)

    # Calculate figure size
    if figsize is None:
        figsize = (max(8, 0.5 * n_types), 6)

    fig, ax = setup_figure(figsize=figsize)

    # Create box plot
    positions = range(len(type_order))
    box_data = [
        df[df["cell_type"] == ct]["confidence"].values
        for ct in type_order
    ]

    bp = ax.boxplot(
        box_data,
        positions=positions,
        patch_artist=True,
        widths=0.6,
    )

    # Color boxes
    for patch, ct in zip(bp["boxes"], type_order):
        patch.set_facecolor(colors.get(ct, "#888888"))
        patch.set_alpha(0.7)

    ax.set_xticks(positions)
    ax.set_xticklabels(type_order, rotation=45, ha="right")

    format_axis_labels(ax, ylabel="Confidence Score")
    despine(ax)

    if title is None:
        title = "Confidence by Cell Type"
    ax.set_title(title)

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    return fig


def plot_confidence_violin(
    adata: ad.AnnData,
    label_column: str,
    confidence_column: str,
    colors: Optional[Dict[str, str]] = None,
    top_n: Optional[int] = 15,
    figsize: Optional[tuple] = None,
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot confidence distribution per cell type as violin plots.

    Parameters
    ----------
    adata : AnnData
        Annotated data with cell type labels and confidence.
    label_column : str
        Column in adata.obs containing cell type labels.
    confidence_column : str
        Column in adata.obs containing confidence values.
    colors : Dict[str, str], optional
        Color mapping for cell types.
    top_n : int, optional, default 15
        Only show top N most frequent cell types.
    figsize : tuple, optional
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
    try:
        import seaborn as sns
    except ImportError:
        raise ImportError("seaborn is required for violin plots")

    if label_column not in adata.obs.columns:
        raise ValueError(f"Label column '{label_column}' not found.")
    if confidence_column not in adata.obs.columns:
        raise ValueError(f"Confidence column '{confidence_column}' not found.")

    df = pd.DataFrame({
        "cell_type": adata.obs[label_column].values,
        "confidence": adata.obs[confidence_column].values,
    })

    type_order = df["cell_type"].value_counts().index.tolist()
    if top_n is not None:
        type_order = type_order[:top_n]
        df = df[df["cell_type"].isin(type_order)]

    n_types = len(type_order)

    if colors is None:
        colors = generate_celltype_palette(type_order)

    if figsize is None:
        figsize = (max(8, 0.6 * n_types), 6)

    fig, ax = setup_figure(figsize=figsize)

    palette = {ct: colors.get(ct, "#888888") for ct in type_order}

    sns.violinplot(
        data=df,
        x="cell_type",
        y="confidence",
        order=type_order,
        palette=palette,
        ax=ax,
    )

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    format_axis_labels(ax, xlabel="", ylabel="Confidence Score")
    despine(ax)

    if title is None:
        title = "Confidence by Cell Type"
    ax.set_title(title)

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    return fig


def plot_model_contribution(
    adata: ad.AnnData,
    model_column: str = "cell_type_model",
    figsize: tuple = (8, 6),
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot which model contributed each cell's prediction.

    Useful for hierarchical or multi-model annotation pipelines
    to see model coverage.

    Parameters
    ----------
    adata : AnnData
        Annotated data with model source column.
    model_column : str, default "cell_type_model"
        Column in adata.obs indicating which model made the prediction.
    figsize : tuple, default (8, 6)
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
    if model_column not in adata.obs.columns:
        raise ValueError(
            f"Model column '{model_column}' not found. "
            f"Available: {list(adata.obs.columns)}"
        )

    counts = adata.obs[model_column].value_counts()

    fig, ax = setup_figure(figsize=figsize)

    colors = generate_celltype_palette(counts.index.tolist())
    bar_colors = [colors.get(m, "#888888") for m in counts.index]

    ax.bar(range(len(counts)), counts.values, color=bar_colors)
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index, rotation=45, ha="right")

    format_axis_labels(ax, ylabel="Number of Cells")
    despine(ax)

    if title is None:
        title = "Model Contribution"
    ax.set_title(title)

    # Add percentages on bars
    total = counts.sum()
    for i, (count, model) in enumerate(zip(counts.values, counts.index)):
        pct = 100 * count / total
        ax.text(
            i,
            count + total * 0.01,
            f"{pct:.1f}%",
            ha="center",
            fontsize=9,
        )

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    return fig
