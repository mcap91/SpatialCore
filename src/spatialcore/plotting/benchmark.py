"""
Benchmark visualization for annotation method comparison.

This module provides functions for comparing annotation methods
and visualizing classification performance.
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


def plot_method_comparison(
    df: pd.DataFrame,
    metrics: List[str] = None,
    method_column: str = "method",
    figsize: tuple = (10, 6),
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot comparison of annotation methods across metrics.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with methods as rows and metrics as columns.
        Should have a column identifying the method.
    metrics : List[str], optional
        Metrics to compare. Default: all numeric columns.
    method_column : str, default "method"
        Column containing method names.
    figsize : tuple, default (10, 6)
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
    >>> from spatialcore.plotting.benchmark import plot_method_comparison
    >>> df = pd.DataFrame({
    ...     "method": ["CellTypist", "HieraType", "Manual"],
    ...     "Accuracy": [0.85, 0.88, 0.92],
    ...     "Silhouette": [0.45, 0.52, 0.48],
    ... })
    >>> fig = plot_method_comparison(df, metrics=["Accuracy", "Silhouette"])
    """
    if method_column not in df.columns:
        raise ValueError(f"Method column '{method_column}' not found.")

    if metrics is None:
        metrics = [c for c in df.columns if c != method_column and np.issubdtype(df[c].dtype, np.number)]

    if not metrics:
        raise ValueError("No numeric metrics found.")

    methods = df[method_column].tolist()
    n_methods = len(methods)
    n_metrics = len(metrics)

    fig, ax = setup_figure(figsize=figsize)

    x = np.arange(n_metrics)
    width = 0.8 / n_methods

    colors = generate_celltype_palette(methods)

    for i, method in enumerate(methods):
        values = df[df[method_column] == method][metrics].values.flatten()
        offset = (i - n_methods / 2 + 0.5) * width
        bars = ax.bar(
            x + offset,
            values,
            width,
            label=method,
            color=colors.get(method, "#888888"),
        )

        # Add value labels on bars
        for bar, val in zip(bars, values):
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.01,
                f"{val:.2f}",
                ha="center",
                va="bottom",
                fontsize=8,
            )

    ax.set_xticks(x)
    ax.set_xticklabels(metrics)
    ax.set_ylabel("Score")
    ax.legend()

    despine(ax)

    if title is None:
        title = "Method Comparison"
    ax.set_title(title)

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    return fig


def plot_confusion_matrix(
    true_labels: np.ndarray,
    pred_labels: np.ndarray,
    labels: Optional[List[str]] = None,
    normalize: bool = True,
    cmap: str = "Blues",
    figsize: Optional[tuple] = None,
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot confusion matrix.

    Parameters
    ----------
    true_labels : np.ndarray
        True class labels.
    pred_labels : np.ndarray
        Predicted class labels.
    labels : List[str], optional
        Class labels. If None, inferred from data.
    normalize : bool, default True
        Normalize by true class (row sums to 1).
    cmap : str, default "Blues"
        Colormap.
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

    Examples
    --------
    >>> from spatialcore.plotting.benchmark import plot_confusion_matrix
    >>> true = adata.obs["true_label"].values
    >>> pred = adata.obs["predicted_label"].values
    >>> fig = plot_confusion_matrix(true, pred, normalize=True)
    """
    from sklearn.metrics import confusion_matrix

    if labels is None:
        labels = sorted(set(true_labels) | set(pred_labels))

    cm = confusion_matrix(true_labels, pred_labels, labels=labels)

    if normalize:
        cm = cm.astype(float) / cm.sum(axis=1, keepdims=True)
        cm = np.nan_to_num(cm)

    n_labels = len(labels)

    if figsize is None:
        figsize = (max(8, n_labels * 0.5), max(6, n_labels * 0.5))

    fig, ax = setup_figure(figsize=figsize)

    im = ax.imshow(cm, cmap=cmap, aspect="auto")

    # Add colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label("Fraction" if normalize else "Count")

    # Add text annotations
    thresh = cm.max() / 2
    for i in range(n_labels):
        for j in range(n_labels):
            val = cm[i, j]
            if normalize:
                text = f"{val:.2f}"
            else:
                text = f"{int(val)}"
            ax.text(
                j,
                i,
                text,
                ha="center",
                va="center",
                color="white" if val > thresh else "black",
                fontsize=8,
            )

    ax.set_xticks(range(n_labels))
    ax.set_yticks(range(n_labels))
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticklabels(labels)

    ax.set_xlabel("Predicted")
    ax.set_ylabel("True")

    if title is None:
        title = "Confusion Matrix"
    ax.set_title(title)

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    return fig


def plot_classification_report(
    true_labels: np.ndarray,
    pred_labels: np.ndarray,
    labels: Optional[List[str]] = None,
    figsize: Optional[tuple] = None,
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot classification metrics (precision, recall, F1) per class.

    Parameters
    ----------
    true_labels : np.ndarray
        True class labels.
    pred_labels : np.ndarray
        Predicted class labels.
    labels : List[str], optional
        Class labels. If None, inferred from data.
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
    from sklearn.metrics import precision_recall_fscore_support

    if labels is None:
        labels = sorted(set(true_labels) | set(pred_labels))

    precision, recall, f1, support = precision_recall_fscore_support(
        true_labels, pred_labels, labels=labels, zero_division=0
    )

    n_labels = len(labels)

    if figsize is None:
        figsize = (max(10, n_labels * 0.5), 6)

    fig, ax = setup_figure(figsize=figsize)

    x = np.arange(n_labels)
    width = 0.25

    ax.bar(x - width, precision, width, label="Precision", color="#3784FE")
    ax.bar(x, recall, width, label="Recall", color="#33FF33")
    ax.bar(x + width, f1, width, label="F1", color="#FF6B6B")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_ylabel("Score")
    ax.set_ylim(0, 1.1)
    ax.legend()

    despine(ax)

    if title is None:
        title = "Classification Metrics by Class"
    ax.set_title(title)

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    return fig


def plot_agreement_heatmap(
    adata: ad.AnnData,
    columns: List[str],
    figsize: Optional[tuple] = None,
    cmap: str = "Greens",
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot agreement matrix between annotation methods.

    Parameters
    ----------
    adata : AnnData
        Annotated data with multiple annotation columns.
    columns : List[str]
        Columns in adata.obs to compare.
    figsize : tuple, optional
        Figure size.
    cmap : str, default "Greens"
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
    for col in columns:
        if col not in adata.obs.columns:
            raise ValueError(f"Column '{col}' not found.")

    n_methods = len(columns)

    # Calculate agreement matrix
    agreement = np.zeros((n_methods, n_methods))
    for i, col1 in enumerate(columns):
        for j, col2 in enumerate(columns):
            agreement[i, j] = (adata.obs[col1] == adata.obs[col2]).mean()

    if figsize is None:
        figsize = (max(6, n_methods * 1.2), max(5, n_methods))

    fig, ax = setup_figure(figsize=figsize)

    im = ax.imshow(agreement, cmap=cmap, vmin=0, vmax=1)
    plt.colorbar(im, ax=ax, label="Agreement")

    # Add text
    for i in range(n_methods):
        for j in range(n_methods):
            ax.text(
                j,
                i,
                f"{agreement[i, j]:.2f}",
                ha="center",
                va="center",
                color="white" if agreement[i, j] > 0.5 else "black",
            )

    ax.set_xticks(range(n_methods))
    ax.set_yticks(range(n_methods))
    ax.set_xticklabels(columns, rotation=45, ha="right")
    ax.set_yticklabels(columns)

    if title is None:
        title = "Method Agreement"
    ax.set_title(title)

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    return fig


def plot_silhouette_by_type(
    adata: ad.AnnData,
    label_column: str,
    embedding_key: str = "X_pca",
    sample_size: int = 5000,
    random_state: int = 42,
    figsize: tuple = (10, 6),
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot silhouette scores by cell type.

    Parameters
    ----------
    adata : AnnData
        Annotated data.
    label_column : str
        Column in adata.obs containing cell type labels.
    embedding_key : str, default "X_pca"
        Key in adata.obsm for embedding.
    sample_size : int, default 5000
        Number of cells to sample (for speed).
    random_state : int, default 42
        Random seed.
    figsize : tuple, default (10, 6)
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
    from sklearn.metrics import silhouette_samples

    if label_column not in adata.obs.columns:
        raise ValueError(f"Label column '{label_column}' not found.")
    if embedding_key not in adata.obsm:
        raise ValueError(f"Embedding '{embedding_key}' not found.")

    # Sample if too large
    if adata.n_obs > sample_size:
        np.random.seed(random_state)
        idx = np.random.choice(adata.n_obs, sample_size, replace=False)
        X = adata.obsm[embedding_key][idx]
        labels = adata.obs[label_column].values[idx]
    else:
        X = adata.obsm[embedding_key]
        labels = adata.obs[label_column].values

    # Calculate silhouette scores
    sil_scores = silhouette_samples(X, labels)

    # Get mean per type
    df = pd.DataFrame({"label": labels, "silhouette": sil_scores})
    type_scores = df.groupby("label")["silhouette"].mean().sort_values()

    fig, ax = setup_figure(figsize=figsize)

    colors = generate_celltype_palette(type_scores.index.tolist())
    bar_colors = [colors.get(ct, "#888888") for ct in type_scores.index]

    y_pos = np.arange(len(type_scores))
    ax.barh(y_pos, type_scores.values, color=bar_colors)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(type_scores.index)
    ax.axvline(0, color="gray", linestyle="--")

    format_axis_labels(ax, xlabel="Silhouette Score")
    despine(ax)

    if title is None:
        title = f"Silhouette Scores by Cell Type\n(mean={sil_scores.mean():.3f})"
    ax.set_title(title)

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    return fig
