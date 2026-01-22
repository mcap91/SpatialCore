"""
Spatial visualization of cell types and confidence.

This module provides functions for visualizing cell annotations
on spatial coordinates.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.colors import Normalize
from matplotlib.collections import LineCollection
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


def plot_domain_distances(
    adata: ad.AnnData,
    source_domain_column: str,
    target_domain_column: Optional[str] = None,
    spatial_key: str = "spatial",
    distance_key: str = "domain_distances",
    top_n_connections: int = 1,
    line_cmap: str = "coolwarm_r",
    line_width: float = 2.0,
    point_size: float = 0.5,
    point_alpha: float = 0.3,
    domain_point_size: float = 3.0,
    domain_point_alpha: float = 0.7,
    figsize: Tuple[float, float] = (14, 12),
    title: Optional[str] = None,
    save: Optional[Union[str, Path]] = None,
) -> Figure:
    """
    Plot spatial domains with lines showing inter-domain distances.

    Shows all cells (non-domain cells in grey), with domain cells colored.
    Draws lines between domain centroids colored by distance.
    Works with any distance metric (minimum, centroid, mean).

    Parameters
    ----------
    adata
        AnnData with spatial coordinates and domain distances computed.
    source_domain_column
        Column in adata.obs containing source domain labels.
    target_domain_column
        Column in adata.obs containing target domain labels.
        If None, uses source_domain_column (intra-domain distances).
    spatial_key
        Key in adata.obsm for spatial coordinates.
    distance_key
        Key in adata.uns containing distance matrix from calculate_domain_distances().
    top_n_connections
        Number of closest connections to show per domain. Default 1 shows only
        the nearest neighbor for each domain. Set to 0 or None to show all.
    line_cmap
        Colormap for distance lines (default coolwarm_r: blue=close, red=far).
    line_width
        Width of distance lines.
    point_size
        Size of background (non-domain) cell scatter points.
    point_alpha
        Transparency of background cell scatter points.
    domain_point_size
        Size of domain cell scatter points.
    domain_point_alpha
        Transparency of domain cell scatter points.
    figsize
        Figure size.
    title
        Plot title.
    save
        Path to save figure.

    Returns
    -------
    Figure
        Matplotlib figure.

    Raises
    ------
    KeyError
        If distance_key not found in adata.uns.
    ValueError
        If domain columns not found in adata.obs.

    Examples
    --------
    >>> from spatialcore.plotting import plot_domain_distances
    >>> fig = plot_domain_distances(
    ...     adata,
    ...     source_domain_column="bcell_domain",
    ...     top_n_connections=1,  # Show only nearest neighbor per domain
    ... )
    """
    if target_domain_column is None:
        target_domain_column = source_domain_column

    if spatial_key not in adata.obsm:
        raise ValueError(
            f"Spatial key '{spatial_key}' not found in adata.obsm. "
            f"Available: {list(adata.obsm.keys())}"
        )

    if source_domain_column not in adata.obs.columns:
        raise ValueError(
            f"Source domain column '{source_domain_column}' not found in adata.obs."
        )

    if target_domain_column not in adata.obs.columns:
        raise ValueError(
            f"Target domain column '{target_domain_column}' not found in adata.obs."
        )

    if distance_key not in adata.uns:
        raise KeyError(
            f"Distance key '{distance_key}' not found in adata.uns. "
            "Run calculate_domain_distances() first."
        )

    dist_data = adata.uns[distance_key]
    if "distance_matrix" not in dist_data:
        raise KeyError(f"'distance_matrix' not found in adata.uns['{distance_key}'].")

    distance_matrix = pd.DataFrame(dist_data["distance_matrix"]).T
    coords = adata.obsm[spatial_key]

    # Compute centroids for all domains
    all_domains = set(distance_matrix.index) | set(distance_matrix.columns)
    centroids = {}

    for domain in all_domains:
        mask_src = adata.obs[source_domain_column] == domain
        mask_tgt = adata.obs[target_domain_column] == domain
        mask = mask_src | mask_tgt

        if mask.any():
            centroids[domain] = coords[mask.values].mean(axis=0)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Get domain masks
    source_mask = adata.obs[source_domain_column].notna()
    target_mask = adata.obs[target_domain_column].notna()
    domain_mask = source_mask | target_mask
    background_mask = ~domain_mask

    # Plot background cells (non-domain) in grey
    if background_mask.any():
        ax.scatter(
            coords[background_mask.values, 0],
            coords[background_mask.values, 1],
            c="lightgrey",
            s=point_size,
            alpha=point_alpha,
            rasterized=True,
            zorder=0,
        )

    # Get unique domains for coloring
    source_domains = adata.obs[source_domain_column].dropna().unique()
    target_domains = adata.obs[target_domain_column].dropna().unique()
    unique_domains = sorted(set(source_domains) | set(target_domains))
    n_colors = min(20, len(unique_domains))
    domain_colors = dict(zip(
        unique_domains,
        plt.cm.tab20(np.linspace(0, 1, n_colors))
    ))

    # Plot domain cells colored by domain
    for domain in unique_domains:
        mask = (adata.obs[source_domain_column] == domain) | (
            adata.obs[target_domain_column] == domain
        )
        if mask.any():
            ax.scatter(
                coords[mask.values, 0],
                coords[mask.values, 1],
                c=[domain_colors[domain]],
                s=domain_point_size,
                alpha=domain_point_alpha,
                rasterized=True,
                zorder=1,
            )

    # Collect line segments - filter to top_n closest per domain
    segments = []
    distances = []
    seen_pairs = set()  # Track (src, tgt) pairs to avoid duplicates

    if top_n_connections and top_n_connections > 0:
        # For each source domain, get top_n closest targets
        for src in distance_matrix.index:
            if src not in centroids:
                continue

            # Get distances from this source to all targets (excluding self)
            row = distance_matrix.loc[src].drop(src, errors='ignore').dropna()
            if row.empty:
                continue

            # Get top_n closest
            closest = row.nsmallest(top_n_connections)

            for tgt, dist in closest.items():
                if tgt not in centroids:
                    continue
                # Avoid duplicate lines (A->B and B->A)
                pair_key = tuple(sorted([src, tgt]))
                if pair_key in seen_pairs:
                    continue
                seen_pairs.add(pair_key)
                segments.append([centroids[src], centroids[tgt]])
                distances.append(dist)
    else:
        # Show all connections
        for src in distance_matrix.index:
            if src not in centroids:
                continue
            for tgt in distance_matrix.columns:
                if tgt not in centroids or src == tgt:
                    continue
                if src >= tgt:  # Avoid duplicates
                    continue

                dist = distance_matrix.loc[src, tgt]
                if pd.isna(dist):
                    continue

                segments.append([centroids[src], centroids[tgt]])
                distances.append(dist)

    if not segments:
        raise ValueError("No valid distance segments to plot.")

    # Normalize distances for coloring
    distances = np.array(distances)
    norm = Normalize(vmin=distances.min(), vmax=distances.max())

    # Create line collection
    lc = LineCollection(
        segments,
        cmap=line_cmap,
        norm=norm,
        linewidths=line_width,
        zorder=2,
    )
    lc.set_array(distances)
    ax.add_collection(lc)

    # Add colorbar
    cbar = plt.colorbar(lc, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label("Distance", fontsize=11)

    # Plot centroids as larger markers with labels
    for domain, centroid in centroids.items():
        ax.scatter(
            centroid[0],
            centroid[1],
            c=[domain_colors.get(domain, "#888888")],
            s=100,
            edgecolors="white",
            linewidths=1.5,
            zorder=3,
        )
        ax.annotate(
            domain.split("_")[-1],
            (centroid[0], centroid[1]),
            fontsize=9,
            fontweight="bold",
            ha="center",
            va="center",
            color="white",
            zorder=4,
        )

    ax.set_xlabel("X (pixels)", fontsize=12)
    ax.set_ylabel("Y (pixels)", fontsize=12)
    ax.set_aspect("equal")
    ax.invert_yaxis()

    if title is None:
        metric = dist_data.get("distance_metric", "unknown")
        title = f"Domain Distances ({metric})"
    ax.set_title(title, fontsize=14, fontweight="bold")

    plt.tight_layout()

    if save:
        save_figure(fig, save)

    return fig
