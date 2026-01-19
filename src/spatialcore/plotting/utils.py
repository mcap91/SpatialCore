"""
Plotting utilities for SpatialCore visualization.

This module provides utility functions for:
1. Color palette generation and loading
2. Figure setup and configuration
3. Saving figures in multiple formats

These utilities are used by all other plotting modules.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union
import json

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from spatialcore.core.logging import get_logger

logger = get_logger(__name__)


# ============================================================================
# Color Palettes
# ============================================================================

# Default high-contrast palette for dark backgrounds
DEFAULT_PALETTE = [
    "#F5252F",  # Red
    "#FB3FFC",  # Magenta
    "#00FFFF",  # Cyan
    "#33FF33",  # Green
    "#FFB300",  # Amber
    "#9966FF",  # Purple
    "#FF6B6B",  # Coral
    "#3784FE",  # Blue
    "#FF8000",  # Orange
    "#66CCCC",  # Teal
    "#CC66FF",  # Lavender
    "#99FF99",  # Light green
    "#FF6699",  # Pink
    "#E3E1E3",  # Light gray
    "#FFB3CC",  # Light pink
    "#E6E680",  # Yellow-green
    "#CC9966",  # Tan
    "#8080FF",  # Periwinkle
    "#FF9999",  # Salmon
    "#66FF99",  # Mint
]

# Colorblind-safe palette (from ColorBrewer)
COLORBLIND_PALETTE = [
    "#E69F00",  # Orange
    "#56B4E9",  # Sky blue
    "#009E73",  # Bluish green
    "#F0E442",  # Yellow
    "#0072B2",  # Blue
    "#D55E00",  # Vermillion
    "#CC79A7",  # Reddish purple
    "#000000",  # Black
]


def generate_celltype_palette(
    cell_types: List[str],
    palette: str = "default",
    custom_colors: Optional[Dict[str, str]] = None,
) -> Dict[str, str]:
    """
    Generate color palette for cell types.

    Parameters
    ----------
    cell_types : List[str]
        List of cell type names.
    palette : str, default "default"
        Palette name: "default", "colorblind", "tab20", "Set1", "Set2", "Set3".
    custom_colors : Dict[str, str], optional
        Custom color overrides for specific cell types.

    Returns
    -------
    Dict[str, str]
        Mapping from cell type name to hex color.

    Examples
    --------
    >>> from spatialcore.plotting.utils import generate_celltype_palette
    >>> colors = generate_celltype_palette(
    ...     ["T cell", "B cell", "Macrophage"],
    ...     custom_colors={"T cell": "#FF0000"},
    ... )
    """
    # Get base palette
    if palette == "default":
        base_colors = DEFAULT_PALETTE
    elif palette == "colorblind":
        base_colors = COLORBLIND_PALETTE
    elif palette in ["tab20", "Set1", "Set2", "Set3"]:
        cmap = plt.get_cmap(palette)
        n_colors = cmap.N if hasattr(cmap, "N") else 20
        base_colors = [plt.colors.rgb2hex(cmap(i)) for i in range(n_colors)]
    else:
        logger.warning(f"Unknown palette '{palette}', using default")
        base_colors = DEFAULT_PALETTE

    if custom_colors is None:
        custom_colors = {}

    color_map = {}
    color_idx = 0

    for cell_type in sorted(cell_types):
        if cell_type in custom_colors:
            color_map[cell_type] = custom_colors[cell_type]
        else:
            color_map[cell_type] = base_colors[color_idx % len(base_colors)]
            color_idx += 1

    return color_map


def load_celltype_palette(path: Path) -> Dict[str, str]:
    """
    Load cell type color palette from JSON file.

    Parameters
    ----------
    path : Path
        Path to JSON file with color mapping.

    Returns
    -------
    Dict[str, str]
        Mapping from cell type name to hex color.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Palette file not found: {path}")

    with open(path, "r") as f:
        colors = json.load(f)

    logger.info(f"Loaded {len(colors)} colors from {path}")
    return colors


def save_celltype_palette(
    colors: Dict[str, str],
    path: Path,
) -> None:
    """
    Save cell type color palette to JSON file.

    Parameters
    ----------
    colors : Dict[str, str]
        Mapping from cell type name to hex color.
    path : Path
        Output path for JSON file.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w") as f:
        json.dump(colors, f, indent=2)

    logger.info(f"Saved {len(colors)} colors to {path}")


# ============================================================================
# Figure Setup
# ============================================================================

def setup_figure(
    figsize: Tuple[float, float] = (8, 6),
    dpi: int = 150,
    style: str = "ticks",
    context: str = "notebook",
    dark_background: bool = False,
) -> Tuple[Figure, Axes]:
    """
    Create figure with consistent styling.

    Parameters
    ----------
    figsize : Tuple[float, float], default (8, 6)
        Figure size in inches (width, height).
    dpi : int, default 150
        Dots per inch for rendering.
    style : str, default "ticks"
        Seaborn style: "ticks", "whitegrid", "darkgrid", "white", "dark".
    context : str, default "notebook"
        Seaborn context: "paper", "notebook", "talk", "poster".
    dark_background : bool, default False
        Use dark background for spatial plots.

    Returns
    -------
    Tuple[Figure, Axes]
        Matplotlib figure and axes.

    Examples
    --------
    >>> from spatialcore.plotting.utils import setup_figure
    >>> fig, ax = setup_figure(figsize=(10, 8), dpi=300)
    >>> ax.scatter(x, y)
    >>> fig.savefig("plot.png")
    """
    try:
        import seaborn as sns
        sns.set_style(style)
        sns.set_context(context)
    except ImportError:
        pass  # Seaborn optional

    if dark_background:
        plt.style.use("dark_background")

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    return fig, ax


def setup_multi_figure(
    nrows: int = 1,
    ncols: int = 1,
    figsize: Optional[Tuple[float, float]] = None,
    dpi: int = 150,
    sharex: bool = False,
    sharey: bool = False,
) -> Tuple[Figure, np.ndarray]:
    """
    Create multi-panel figure.

    Parameters
    ----------
    nrows : int, default 1
        Number of subplot rows.
    ncols : int, default 1
        Number of subplot columns.
    figsize : Tuple[float, float], optional
        Figure size. If None, auto-calculated from panel count.
    dpi : int, default 150
        Dots per inch.
    sharex : bool, default False
        Share x-axis across subplots.
    sharey : bool, default False
        Share y-axis across subplots.

    Returns
    -------
    Tuple[Figure, np.ndarray]
        Figure and array of axes.
    """
    if figsize is None:
        figsize = (4 * ncols, 4 * nrows)

    fig, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=figsize,
        dpi=dpi,
        sharex=sharex,
        sharey=sharey,
    )

    return fig, np.atleast_2d(axes) if nrows == 1 or ncols == 1 else axes


# ============================================================================
# Figure Saving
# ============================================================================

def save_figure(
    fig: Figure,
    path: Union[str, Path],
    formats: List[str] = None,
    dpi: int = 300,
    bbox_inches: str = "tight",
    transparent: bool = False,
) -> List[Path]:
    """
    Save figure in multiple formats.

    Parameters
    ----------
    fig : Figure
        Matplotlib figure to save.
    path : str or Path
        Base output path (without extension).
    formats : List[str], optional
        Output formats. Default: ["png"].
    dpi : int, default 300
        Resolution for raster formats.
    bbox_inches : str, default "tight"
        Bounding box setting.
    transparent : bool, default False
        Transparent background.

    Returns
    -------
    List[Path]
        Paths to saved files.

    Examples
    --------
    >>> from spatialcore.plotting.utils import save_figure
    >>> paths = save_figure(fig, "output/my_plot", formats=["png", "pdf", "svg"])
    """
    if formats is None:
        formats = ["png"]

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    saved_paths = []
    for fmt in formats:
        output_path = path.with_suffix(f".{fmt}")
        fig.savefig(
            output_path,
            format=fmt,
            dpi=dpi,
            bbox_inches=bbox_inches,
            transparent=transparent,
        )
        saved_paths.append(output_path)
        logger.debug(f"Saved figure to: {output_path}")

    logger.info(f"Saved figure to {len(saved_paths)} formats: {path}")
    return saved_paths


def close_figure(fig: Figure) -> None:
    """
    Close figure to free memory.

    Parameters
    ----------
    fig : Figure
        Matplotlib figure to close.
    """
    plt.close(fig)


# ============================================================================
# Axis Formatting
# ============================================================================

def format_axis_labels(
    ax: Axes,
    xlabel: Optional[str] = None,
    ylabel: Optional[str] = None,
    title: Optional[str] = None,
    fontsize: int = 12,
) -> Axes:
    """
    Format axis labels and title.

    Parameters
    ----------
    ax : Axes
        Matplotlib axes.
    xlabel : str, optional
        X-axis label.
    ylabel : str, optional
        Y-axis label.
    title : str, optional
        Plot title.
    fontsize : int, default 12
        Font size for labels.

    Returns
    -------
    Axes
        Modified axes.
    """
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=fontsize)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=fontsize)
    if title:
        ax.set_title(title, fontsize=fontsize + 2)
    return ax


def despine(ax: Axes, top: bool = True, right: bool = True) -> Axes:
    """
    Remove spines from axes.

    Parameters
    ----------
    ax : Axes
        Matplotlib axes.
    top : bool, default True
        Remove top spine.
    right : bool, default True
        Remove right spine.

    Returns
    -------
    Axes
        Modified axes.
    """
    if top:
        ax.spines["top"].set_visible(False)
    if right:
        ax.spines["right"].set_visible(False)
    return ax
