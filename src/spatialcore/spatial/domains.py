"""Spatial domain creation using R's sf package.

This module provides functionality to create spatial domains (contiguous regions)
from groups of cells in spatial transcriptomics data using the Buffer-Union-Shrink
algorithm implemented in R's sf package.

Key Features:
    - Creates smooth, contiguous domain polygons
    - Unified filter_expression (ontology IDs, column equality, boolean columns)
    - Dual cell filtering (min_target_cells_domain, min_total_cells_domain)
    - Heterogeneity analysis by assigning ALL cells to domains
    - Deterministic domain numbering (largest = 1) for reproducibility
    - Platform-aware defaults for CosMx, Xenium, and Visium
    - Domain expansion warning when assigned/target ratio is high

Notes
-----
This module requires R with the following packages:
    - sf (Simple Features for geometry operations)
    - concaveman (Concave hull algorithm)
    - dplyr, purrr (Data manipulation)
    - jsonlite (JSON output)

Install R packages:
    install.packages(c("sf", "concaveman", "dplyr", "purrr", "jsonlite"))

References
----------
- sf package: https://r-spatial.github.io/sf/
- concaveman: https://github.com/joelgombin/concaveman
- Buffer-Union-Shrink: Common GIS technique for region delineation
- DOMAINS.md specification: C:\\SpatialCore\\Data\\ref_docs\\DOMAINS.md

Examples
--------
>>> import scanpy as sc
>>> from spatialcore.spatial import make_spatial_domains
>>> adata = sc.read_h5ad("annotated.h5ad")
>>> # Create B cell domains using ontology ID
>>> adata = make_spatial_domains(
...     adata,
...     filter_expression="CL:0000236",  # B cell
...     domain_prefix="Bcell",
... )
>>> # Compound filter expression
>>> adata = make_spatial_domains(
...     adata,
...     filter_expression="CL:0000236 & NCIT:C4349",  # B cell AND tumor
...     domain_prefix="Bcell_Tumor",
...     platform="xenium",
... )
>>> # Column-based filtering with dual cell thresholds
>>> adata = make_spatial_domains(
...     adata,
...     filter_expression="metagene_cluster == 1",
...     domain_prefix="Tumor",
...     min_target_cells_domain=10,
...     min_total_cells_domain=15000,
... )
"""

import json
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd
import anndata as ad

from spatialcore.core.logging import get_logger
from spatialcore.core.metadata import update_metadata
from spatialcore.r_bridge import check_r_available, run_r_script
from spatialcore.r_bridge.subprocess_runner import RNotFoundError, RExecutionError

logger = get_logger(__name__)

# Path to R functions
R_FUNCTIONS_PATH = Path(__file__).parent / "r_functions.R"

# Platform-specific defaults for cell_dist parameter
# These values represent typical cell spacing for each platform
PLATFORM_DEFAULTS: Dict[str, float] = {
    "cosmx": 400.0,   # CosMx uses pixel coordinates (typically 0-120000)
    "xenium": 50.0,   # Xenium uses microns (typically 0-10000)
    "visium": 200.0,  # Visium uses spot array units (typically 0-50000)
}

# Coordinate range thresholds for platform auto-detection
# Based on typical coordinate ranges for each platform
PLATFORM_COORD_RANGES: Dict[str, Tuple[float, float]] = {
    "cosmx": (50000.0, float("inf")),   # >50000 suggests CosMx pixels
    "xenium": (0.0, 15000.0),            # 0-15000 suggests Xenium microns
    "visium": (15000.0, 50000.0),        # 15000-50000 suggests Visium
}


def _detect_platform(adata: ad.AnnData) -> Optional[str]:
    """
    Auto-detect spatial platform based on coordinate ranges.

    Detection is based on the maximum coordinate value in the spatial data:
    - CosMx: coordinates typically 0-120000 (pixels)
    - Xenium: coordinates typically 0-10000 (microns)
    - Visium: coordinates typically 0-50000 (spot array units)

    Parameters
    ----------
    adata
        AnnData object with spatial coordinates in adata.obsm['spatial'].

    Returns
    -------
    Optional[str]
        Detected platform name ('cosmx', 'xenium', 'visium') or None if
        detection is ambiguous or coordinates are not available.

    Notes
    -----
    This heuristic is based on typical coordinate ranges observed in practice.
    For reliable results, explicitly specify the platform parameter.
    """
    if "spatial" not in adata.obsm:
        return None

    spatial_coords = adata.obsm["spatial"]
    if spatial_coords.shape[0] == 0:
        return None

    # Use max coordinate value for detection
    max_coord = np.max(np.abs(spatial_coords))

    # Check against platform ranges
    if max_coord > PLATFORM_COORD_RANGES["cosmx"][0]:
        return "cosmx"
    elif max_coord <= PLATFORM_COORD_RANGES["xenium"][1]:
        return "xenium"
    elif max_coord <= PLATFORM_COORD_RANGES["visium"][1]:
        return "visium"

    return None


def _get_platform_defaults(platform: str) -> float:
    """
    Get default cell_dist value for a given platform.

    Parameters
    ----------
    platform
        Platform name ('cosmx', 'xenium', 'visium').

    Returns
    -------
    float
        Default cell_dist value for the platform.

    Raises
    ------
    ValueError
        If platform is not recognized.
    """
    platform_lower = platform.lower()
    if platform_lower not in PLATFORM_DEFAULTS:
        valid_platforms = list(PLATFORM_DEFAULTS.keys())
        raise ValueError(
            f"Unknown platform '{platform}'. "
            f"Valid platforms are: {valid_platforms}"
        )
    return PLATFORM_DEFAULTS[platform_lower]


def _evaluate_filter_expression(
    filter_expression: str,
    adata: ad.AnnData,
) -> pd.Series:
    """
    Evaluate a filter expression to produce a boolean mask.

    Supports multiple filter formats:
    - Ontology IDs: "CL:0000236" (matches cells with this ID in ontology columns)
    - Boolean expressions: "CL:0000236 & NCIT:C4349"
    - Column equality: "cell_type == 'B cell'" or "metagene_cluster == 1"
    - Boolean column: "is_tumor" (column containing True/False)
    - Mixed: "CL:0000236 & is_tumor"

    Parameters
    ----------
    filter_expression
        Filter expression string.
    adata
        AnnData object to filter.

    Returns
    -------
    pd.Series
        Boolean mask with True for cells matching the filter.
    """
    import re

    expr = filter_expression.strip()

    # Check if this looks like an ontology expression (contains CL:, NCIT:, UBERON:, etc.)
    ontology_pattern = r'[A-Z]+:[0-9]+'
    has_ontology_ids = bool(re.search(ontology_pattern, expr))

    if has_ontology_ids:
        # Use the ontology expression evaluator
        from spatialcore.ontology.expression import evaluate_ontology_expression
        return evaluate_ontology_expression(expr, adata)

    # Check if this is a simple column equality (e.g., "cell_type == 'B cell'")
    equality_match = re.match(r"^(\w+)\s*==\s*['\"]?(.+?)['\"]?$", expr)
    if equality_match:
        col_name = equality_match.group(1)
        col_value = equality_match.group(2)
        if col_name not in adata.obs.columns:
            raise ValueError(
                f"Column '{col_name}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)[:10]}..."
            )
        return adata.obs[col_name] == col_value

    # Check if this is a simple boolean column reference
    if expr in adata.obs.columns:
        col = adata.obs[expr]
        if col.dtype == bool or set(col.dropna().unique()).issubset({True, False, "True", "False"}):
            # Convert to boolean
            if col.dtype == object:
                return col.map({"True": True, "False": False, True: True, False: False}).fillna(False)
            return col.fillna(False).astype(bool)
        else:
            raise ValueError(
                f"Column '{expr}' exists but is not boolean. "
                f"Use equality syntax like \"{expr} == 'value'\" instead."
            )

    # Try to evaluate as a complex pandas expression
    # This handles expressions like "cell_type == 'B cell' & is_tumor"
    try:
        result = adata.obs.eval(expr)
        return result.astype(bool)
    except Exception as e:
        raise ValueError(
            f"Could not evaluate filter expression: '{filter_expression}'. "
            f"Error: {e}\n"
            "Supported formats:\n"
            "  - Ontology ID: 'CL:0000236'\n"
            "  - Column equality: \"cell_type == 'B cell'\"\n"
            "  - Boolean column: 'is_tumor'\n"
            "  - Compound: 'CL:0000236 & is_tumor'"
        ) from e


def _generate_domain_prefix(filter_expression: str) -> str:
    """
    Generate a domain prefix from a filter expression.

    Parameters
    ----------
    filter_expression
        Filter expression string.

    Returns
    -------
    str
        A sanitized prefix for domain names.
    """
    import re

    # Check for ontology IDs first
    ontology_match = re.search(r'([A-Z]+):([0-9]+)', filter_expression)
    if ontology_match:
        # Use "CL_0000236" format
        return f"{ontology_match.group(1)}_{ontology_match.group(2)}"

    # Check for column equality
    equality_match = re.match(r"^(\w+)\s*==\s*['\"]?(.+?)['\"]?$", filter_expression)
    if equality_match:
        # Use the value, sanitized
        value = equality_match.group(2)
        return re.sub(r'[^a-zA-Z0-9_]', '_', value)

    # Use first part of expression, sanitized
    prefix = re.sub(r'[^a-zA-Z0-9_]', '_', filter_expression[:20])
    return prefix if prefix else "domain"


def make_spatial_domains(
    adata: ad.AnnData,
    filter_expression: Optional[str] = None,
    cell_dist_um: Optional[float] = None,
    shrink_margin_um: float = 25.0,
    domain_prefix: Optional[str] = None,
    min_target_cells_domain: int = 10,
    min_total_cells_domain: Optional[int] = None,
    output_column: str = "spatial_domain",
    assign_all_cells: bool = True,
    domain_expansion_warn_ratio: float = 10.0,
    r_functions_path: Optional[Union[str, Path]] = None,
    copy: bool = False,
    platform: Optional[Literal["cosmx", "xenium", "visium"]] = None,
) -> ad.AnnData:
    """
    Create spatial domains from cell groupings using R's sf package.

    Uses the Buffer-Union-Shrink algorithm:

    1. Buffer each target cell by cell_dist_um
    2. Union overlapping buffers into polygons
    3. Shrink by (cell_dist_um - shrink_margin_um) to tighten boundaries
    4. Create concave hulls for natural shapes

    This function supports platform-aware defaults for the ``cell_dist_um``
    parameter. When ``cell_dist_um`` is not explicitly provided, the function
    will auto-detect the spatial platform based on coordinate ranges and
    apply appropriate defaults.

    Parameters
    ----------
    adata
        AnnData object with spatial coordinates in ``adata.obsm['spatial']``.
        Coordinates should be in the native units of the platform (pixels
        for CosMx, microns for Xenium, array units for Visium).
    filter_expression
        Expression to filter target cells. Supports multiple formats:

        - Ontology IDs: ``"CL:0000236"`` (B cell)
        - Boolean expressions: ``"CL:0000236 & NCIT:C4349"`` (B cell AND tumor)
        - Column equality: ``"cell_type == 'B cell'"``
        - Boolean column: ``"is_tumor"``
        - Mixed: ``"CL:0000236 & is_tumor"``

        Operators: ``&`` (AND), ``|`` (OR), ``~`` (NOT), parentheses for grouping.
    cell_dist_um
        Buffer distance in coordinate units. If None (default), uses
        platform-specific defaults:

        - CosMx: 400 pixels
        - Xenium: 50 microns
        - Visium: 200 array units

        If a value is explicitly provided, it overrides platform defaults.
        Larger values create more connected domains.
    shrink_margin_um
        Margin to keep when shrinking. Default: 25.
        Actual shrink = cell_dist_um - shrink_margin_um.
        Larger values create tighter polygons.
    domain_prefix
        Prefix for domain names. Auto-generated from filter_expression if None.
    min_target_cells_domain
        Minimum TARGET cells per domain. Default: 10. Counts only cells
        matching the filter_expression. Smaller domains are merged with
        neighbors or removed.
    min_total_cells_domain
        Minimum TOTAL cells per domain. Default: None (no filter).
        Counts ALL cells after assign_all_cells expansion. Use to keep
        only large spatial regions.
    output_column
        Column name for domain assignments. Default: 'spatial_domain'.
    assign_all_cells
        If True, assign ALL cells to domains (for heterogeneity analysis).
        If False, only assign cells matching the filter. Default: True.
    domain_expansion_warn_ratio
        Warn if (assigned cells / target cells) exceeds this ratio.
        Default: 10.0. Helps catch unexpected domain expansion.
    r_functions_path
        Path to R functions file. Uses bundled r_functions.R if None.
    copy
        If True, operate on a copy of adata. Default: False.
    platform
        Spatial platform type. One of 'cosmx', 'xenium', 'visium', or None.
        If None (default), platform is auto-detected based on coordinate
        ranges:

        - CosMx: max coordinate > 50,000 (pixel coordinates)
        - Xenium: max coordinate <= 15,000 (micron coordinates)
        - Visium: max coordinate 15,000-50,000 (array units)

        Explicitly specifying the platform ensures correct defaults are used.

    Returns
    -------
    AnnData
        AnnData with domain assignments in ``adata.obs[output_column]``.
        Cells not within any domain have NaN. Domains are numbered by size
        (largest domain = 1) for reproducible workflows.

    Raises
    ------
    ValueError
        If spatial coordinates are missing, filter_expression is not provided,
        no cells match the filter, or an unknown platform is specified.
    RNotFoundError
        If R is not installed or not in PATH.
    RExecutionError
        If R script fails.

    Notes
    -----
    Requires R with sf, concaveman, dplyr, purrr, jsonlite packages.

    R's sf package is used instead of Python's shapely because sf produces
    correct polygon geometries where shapely fails with buffer operations
    on complex point clouds.

    Platform Detection Heuristics
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The auto-detection uses the maximum absolute coordinate value:

    - **CosMx** uses pixel coordinates, typically ranging 0-120,000
    - **Xenium** uses micron coordinates, typically ranging 0-10,000
    - **Visium** uses spot array units, typically ranging 0-50,000

    For datasets with unusual coordinate ranges, explicitly specify the
    ``platform`` parameter to ensure correct defaults.

    Deterministic Domain Numbering
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Domains are renumbered by cell count after filtering, with the largest
    domain always numbered 1. This ensures reproducible workflows where
    ``{prefix}_1`` always refers to the largest domain.

    Examples
    --------
    Basic usage with ontology ID:

    >>> import scanpy as sc
    >>> from spatialcore.spatial import make_spatial_domains
    >>> adata = sc.read_h5ad("annotated.h5ad")
    >>> # Create B cell domains (CL:0000236 = B cell)
    >>> adata = make_spatial_domains(
    ...     adata,
    ...     filter_expression="CL:0000236",
    ...     domain_prefix="Bcell",
    ... )
    >>> print(adata.obs["spatial_domain"].value_counts())
    Bcell_1    150
    Bcell_2     89
    Name: spatial_domain, dtype: int64

    Compound filter expression:

    >>> # B cells in tumor regions
    >>> adata = make_spatial_domains(
    ...     adata,
    ...     filter_expression="CL:0000236 & NCIT:C4349",
    ...     domain_prefix="Bcell_Tumor",
    ... )

    Column-based filtering:

    >>> # Filter by column value
    >>> adata = make_spatial_domains(
    ...     adata,
    ...     filter_expression="cell_type == 'B cell'",
    ...     domain_prefix="Bcell",
    ... )

    Dual cell filtering:

    >>> # Keep only large tumor regions (>15000 total cells)
    >>> adata = make_spatial_domains(
    ...     adata,
    ...     filter_expression="metagene_cluster == 1",
    ...     domain_prefix="Tumor",
    ...     min_target_cells_domain=10,
    ...     min_total_cells_domain=15000,
    ... )

    See Also
    --------
    calculate_domain_distances : Compute distances between domains.
    get_domain_summary : Get summary statistics for spatial domains.
    """
    # Validate inputs
    if "spatial" not in adata.obsm:
        raise ValueError(
            "adata.obsm['spatial'] not found. "
            "Spatial coordinates are required for domain creation."
        )

    if filter_expression is None:
        raise ValueError(
            "'filter_expression' must be provided. Examples:\n"
            "  - Ontology ID: 'CL:0000236' (B cell)\n"
            "  - Boolean expression: 'CL:0000236 & NCIT:C4349'\n"
            "  - Column equality: \"cell_type == 'B cell'\"\n"
            "  - Boolean column: 'is_tumor'"
        )

    # Check R availability
    if not check_r_available():
        raise RNotFoundError(
            "R is not installed or not in PATH. "
            "Spatial domain creation requires R with packages: "
            "sf, concaveman, dplyr, purrr, jsonlite.\n"
            "Install R via: mamba install -c conda-forge r-base r-sf r-concaveman r-dplyr r-purrr r-jsonlite"
        )

    # Validate platform parameter if provided
    if platform is not None:
        platform_lower = platform.lower()
        if platform_lower not in PLATFORM_DEFAULTS:
            valid_platforms = list(PLATFORM_DEFAULTS.keys())
            raise ValueError(
                f"Unknown platform '{platform}'. "
                f"Valid platforms are: {valid_platforms}"
            )

    # Resolve platform and cell_dist_um
    # If cell_dist_um is explicitly provided, use it regardless of platform
    # Otherwise, detect platform and use platform-specific defaults
    effective_cell_dist_um: float
    effective_platform: Optional[str] = platform

    if cell_dist_um is not None:
        # User provided explicit cell_dist_um - use it directly
        effective_cell_dist_um = cell_dist_um
        logger.debug(f"Using user-provided cell_dist_um={cell_dist_um}")
    else:
        # No explicit cell_dist_um - detect platform and apply defaults
        if platform is None:
            # Auto-detect platform
            detected_platform = _detect_platform(adata)
            if detected_platform is not None:
                effective_platform = detected_platform
                effective_cell_dist_um = _get_platform_defaults(detected_platform)
                max_coord = np.max(np.abs(adata.obsm["spatial"]))
                logger.info(
                    f"Auto-detected platform '{detected_platform}' "
                    f"(max coordinate: {max_coord:.1f}), "
                    f"using cell_dist_um={effective_cell_dist_um}"
                )
            else:
                raise ValueError(
                    "Could not auto-detect platform from coordinate ranges. "
                    "Provide 'platform' or 'cell_dist_um' explicitly."
                )
        else:
            # Platform specified, use its defaults
            effective_cell_dist_um = _get_platform_defaults(platform)
            logger.info(
                f"Using platform '{platform}' defaults: cell_dist_um={effective_cell_dist_um}"
            )

    # Handle copy
    adata = adata.copy() if copy else adata

    # Resolve R functions path
    r_path = Path(r_functions_path) if r_functions_path else R_FUNCTIONS_PATH
    if not r_path.exists():
        raise FileNotFoundError(f"R functions file not found: {r_path}")

    logger.info("Creating spatial domains using R's sf package")

    # Parse and evaluate filter_expression
    logger.info(f"Evaluating filter expression: {filter_expression}")
    mask = _evaluate_filter_expression(filter_expression, adata)
    n_target_cells = mask.sum()

    if n_target_cells == 0:
        raise ValueError(
            f"No cells match filter expression: '{filter_expression}'. "
            "Check that column names and values are correct."
        )

    # Create temporary filter column for R
    effective_group = "_filter"
    effective_group_subset = "True"
    adata.obs[effective_group] = mask.astype(str)
    logger.info(f"Filter expression matched {n_target_cells:,} cells")

    # Determine domain prefix from filter_expression if not provided
    if domain_prefix is None:
        domain_prefix = _generate_domain_prefix(filter_expression)

    # Prepare data for R
    with tempfile.TemporaryDirectory() as tmpdir:
        input_csv = Path(tmpdir) / "input.csv"
        output_csv = Path(tmpdir) / "output.csv"

        # Extract spatial coordinates
        spatial_coords = adata.obsm["spatial"]
        if spatial_coords.shape[1] < 2:
            raise ValueError(
                f"Spatial coordinates must have at least 2 columns, "
                f"got shape {spatial_coords.shape}"
            )

        # Build input dataframe
        input_df = pd.DataFrame({
            "cell": adata.obs.index,
            "x": spatial_coords[:, 0],
            "y": spatial_coords[:, 1],
            effective_group: adata.obs[effective_group].values,
        })

        # Write input CSV
        input_df.to_csv(input_csv, index=False)
        logger.debug(f"Wrote {len(input_df):,} cells to {input_csv}")

        # Build R code to call make_spatial_domains
        group_subset_str = str(effective_group_subset).replace("'", "\\'")
        min_total_cells_r = "NULL" if min_total_cells_domain is None else str(min_total_cells_domain)
        r_code = f"""
source('{str(r_path).replace(chr(92), "/")}')

make_spatial_domains(
    input_csv = '{str(input_csv).replace(chr(92), "/")}',
    output_csv = '{str(output_csv).replace(chr(92), "/")}',
    group = '{effective_group}',
    group_subset = '{group_subset_str}',
    cell_dist = {effective_cell_dist_um},
    shrink_margin = {shrink_margin_um},
    domain_prefix = '{domain_prefix}',
    min_target_cells_domain = {min_target_cells_domain},
    min_total_cells_domain = {min_total_cells_r},
    assign_all_cells = {'TRUE' if assign_all_cells else 'FALSE'}
)
"""

        # Execute R via subprocess
        logger.debug("Executing R domain creation...")
        from spatialcore.r_bridge import run_r_code
        result = run_r_code(r_code, timeout=1200)

        if "stdout" in result and "parse_error" in result:
            logger.warning(f"R output parsing issue: {result.get('parse_error')}")

        # Read output CSV
        if not output_csv.exists():
            raise RExecutionError(
                "R domain creation did not produce output file. "
                f"R output: {result.get('stdout', '')}"
            )

        output_df = pd.read_csv(output_csv)
        logger.debug(f"Read {len(output_df):,} cells from R output")

        # Merge domain assignments back to adata
        output_df["cell"] = output_df["cell"].astype(str)
        domain_map = dict(zip(output_df["cell"], output_df["domain"]))

        if assign_all_cells is False:
            assigned_count = output_df["domain"].notna().sum()
            if assigned_count == 0:
                raise ValueError(
                    "No cells were assigned to any domain. "
                    "Try relaxing filters, adjusting cell_dist_um, or setting "
                    "assign_all_cells=True."
                )

        # Map domains to adata
        adata.obs[output_column] = adata.obs.index.map(domain_map)

        # Convert "NA" strings to actual NaN
        adata.obs[output_column] = adata.obs[output_column].replace(
            {"NA": np.nan, "domain_NA": np.nan}
        )

        # Also handle domain_prefix_NA pattern
        na_pattern = f"{domain_prefix}_NA"
        adata.obs[output_column] = adata.obs[output_column].replace({na_pattern: np.nan})

        # Renumber domains sequentially (1 to m) by cell count
        domain_counts = adata.obs[output_column].value_counts()
        if len(domain_counts) > 0:
            # Create mapping: old_name -> domain_prefix_1, domain_prefix_2, ...
            renumber_map = {}
            for i, old_name in enumerate(domain_counts.index, start=1):
                renumber_map[old_name] = f"{domain_prefix}_{i}"

            # Apply renumbering
            adata.obs[output_column] = adata.obs[output_column].map(
                lambda x: renumber_map.get(x, x) if pd.notna(x) else x
            )
            logger.debug(f"Renumbered {len(renumber_map)} domains to sequential 1-{len(renumber_map)}")

        # Clean up temporary filter column
        if "_filter" in adata.obs.columns:
            del adata.obs["_filter"]

    # Calculate summary statistics
    n_domains = adata.obs[output_column].nunique()
    n_assigned = adata.obs[output_column].notna().sum()
    domains_list = adata.obs[output_column].dropna().unique().tolist()

    logger.info(
        f"Created {n_domains} domains, assigned {n_assigned:,}/{adata.n_obs:,} cells "
        f"({100 * n_assigned / adata.n_obs:.1f}%)"
    )

    # Domain expansion warning
    if n_target_cells > 0:
        expansion_ratio = n_assigned / n_target_cells
        if expansion_ratio > domain_expansion_warn_ratio:
            logger.warning(
                f"Domain expansion ratio {expansion_ratio:.1f}x exceeds threshold "
                f"({domain_expansion_warn_ratio}x). This means {n_assigned:,} cells were "
                f"assigned to domains defined by only {n_target_cells:,} target cells. "
                "Review assign_all_cells setting if this is unexpected."
            )

    # Update metadata
    update_metadata(
        adata,
        function_name="make_spatial_domains",
        parameters={
            "filter_expression": filter_expression,
            "cell_dist_um": effective_cell_dist_um,
            "cell_dist_um_user_provided": cell_dist_um is not None,
            "platform": effective_platform,
            "platform_user_provided": platform is not None,
            "shrink_margin_um": shrink_margin_um,
            "domain_prefix": domain_prefix,
            "min_target_cells_domain": min_target_cells_domain,
            "min_total_cells_domain": min_total_cells_domain,
            "output_column": output_column,
            "assign_all_cells": assign_all_cells,
            "domain_expansion_warn_ratio": domain_expansion_warn_ratio,
        },
        outputs={
            "obs": output_column,
            "n_domains": n_domains,
            "n_cells_assigned": int(n_assigned),
            "n_target_cells": int(n_target_cells),
            "domains": domains_list,
        },
    )

    return adata


def get_domain_summary(
    adata: ad.AnnData,
    domain_column: str = "spatial_domain",
) -> pd.DataFrame:
    """
    Get summary statistics for spatial domains.

    Parameters
    ----------
    adata
        AnnData with domain assignments.
    domain_column
        Column containing domain labels.

    Returns
    -------
    pd.DataFrame
        Summary with columns: domain, n_cells, percent, centroid_x, centroid_y.

    Examples
    --------
    >>> from spatialcore.spatial import get_domain_summary
    >>> summary = get_domain_summary(adata, "spatial_domain")
    >>> print(summary)
             domain  n_cells  percent  centroid_x  centroid_y
    0       Bcell_1      150    15.0      1234.5      5678.9
    1       Bcell_2       89     8.9      2345.6      6789.0
    """
    if domain_column not in adata.obs.columns:
        raise ValueError(
            f"Column '{domain_column}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    if "spatial" not in adata.obsm:
        raise ValueError(
            "adata.obsm['spatial'] not found. "
            f"Available keys: {list(adata.obsm.keys())}"
        )

    spatial = adata.obsm["spatial"]
    domains = adata.obs[domain_column]

    summaries = []
    for domain in domains.dropna().unique():
        mask = domains == domain
        n_cells = mask.sum()
        coords = spatial[mask.values]

        summaries.append({
            "domain": domain,
            "n_cells": n_cells,
            "percent": 100 * n_cells / len(domains),
            "centroid_x": coords[:, 0].mean(),
            "centroid_y": coords[:, 1].mean(),
        })

    return pd.DataFrame(summaries).sort_values("n_cells", ascending=False)
