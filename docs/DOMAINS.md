# Spatial Domains: Standalone Implementation Guide

**Date**: 2025-12-16
**Status**: Production-ready, documented for standalone extraction
**Context**: Spatial domain creation and distance computation for spatial biology workflows

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Architecture Overview](#architecture-overview)
3. [Why R's sf Package (Not Python Shapely)](#why-rs-sf-package-not-python-shapely)
4. [Core Components](#core-components)
5. [Algorithm Details](#algorithm-details)
6. [Python-R Integration](#python-r-integration)
7. [Domain Distance Computation](#domain-distance-computation)
8. [Standalone Extraction Guide](#standalone-extraction-guide)
9. [API Reference](#api-reference)
10. [Usage Examples](#usage-examples)
11. [Troubleshooting](#troubleshooting)
12. [References](#references)

---

## Executive Summary

### What This Does

Creates **spatial domains** (contiguous regions) from groups of cells in spatial transcriptomics data, then computes **distances between domains** for spatial relationship analysis.

**Key Use Cases**:
- Tumor microenvironment analysis (immune infiltration distance to tumor)
- Cell niche identification (B cell regions, T cell zones)
- Spatial heterogeneity quantification
- Therapeutic target identification by spatial proximity

### Key Features

| Feature | Description |
|---------|-------------|
| **Buffer-Union-Shrink Algorithm** | Creates smooth, contiguous domain polygons |
| **R's sf Package** | Superior geometry handling vs Python shapely |
| **Ontology-Aware Filtering** | Filter by CL/NCIT/UBERON codes or boolean expressions |
| **Heterogeneity Analysis** | Assign ALL cells to domains (not just target cells) |
| **KD-Tree Distance Computation** | O(n log n) optimized distance calculations |
| **Multi-Domain Distance Matrix** | Source → Target domain relationships |

### Why R Over Python

After extensive testing, **R's sf package produces correct polygon geometries** where Python's shapely fails with buffer operations on complex point clouds. See [Why R's sf Package](#why-rs-sf-package-not-python-shapely) for details.

---

## Architecture Overview

### System Diagram

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                        SPATIAL DOMAINS SYSTEM                                │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                              │
│  ┌────────────────┐     ┌────────────────┐     ┌────────────────────────┐   │
│  │  Python Caller │────▶│   R Bridge     │────▶│  R sf Functions        │   │
│  │                │     │                │     │                        │   │
│  │  make_spatial_ │     │  subprocess    │     │  MakeDomains()         │   │
│  │  domains()     │     │  + Rscript     │     │  ReduceDomains()       │   │
│  │                │     │                │     │                        │   │
│  │  ontology      │     │  h5ad ↔ CSV    │     │  st_buffer()           │   │
│  │  expression    │     │  serialization │     │  st_union()            │   │
│  │  evaluation    │     │                │     │  concaveman()          │   │
│  └────────────────┘     └────────────────┘     └────────────────────────┘   │
│         │                                                │                   │
│         ▼                                                ▼                   │
│  ┌────────────────────────────────────────────────────────────────────────┐ │
│  │                        OUTPUT: Cell → Domain Mapping                    │ │
│  │                                                                         │ │
│  │  adata.obs['spatial_domain'] = ['MS4A1_1', 'MS4A1_2', ..., NaN, ...]  │ │
│  │  adata.uns['available_domain_columns'] = ['spatial_domain', ...]       │ │
│  └────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
│  ┌────────────────────────────────────────────────────────────────────────┐ │
│  │                    DOMAIN DISTANCES (Python)                            │ │
│  │                                                                         │ │
│  │  calculate_domain_distances()                                           │ │
│  │    └── KD-Tree optimization                                             │ │
│  │    └── Per-cell: distance_to_target, nearest_target_domain             │ │
│  │    └── Matrix: domain_distances[source][target]                         │ │
│  └────────────────────────────────────────────────────────────────────────┘ │
│                                                                              │
└─────────────────────────────────────────────────────────────────────────────┘
```

### Data Flow

```
1. Input: h5ad with spatial coordinates + cell annotations
         ↓
2. Python: Evaluate ontology_expression filter (if provided)
         ↓
3. R: MakeDomains() - Buffer-Union-Shrink polygon creation
         ↓
4. R: ReduceDomains() - Merge small domains
         ↓
5. Output: CSV with cell → domain mapping
         ↓
6. Python: Merge mapping back to h5ad
         ↓
7. Python: calculate_domain_distances() (optional)
```

---

## Why R's sf Package (Not Python Shapely)

### The Problem with Shapely

When applying `buffer()` → `union()` → `negative buffer()` operations on complex cell point clouds, **Python's shapely library produces incorrect geometries**:

1. **Self-intersecting polygons** after shrink operations
2. **Lost interior holes** in complex shapes
3. **Invalid geometries** that break downstream operations
4. **Unpredictable failures** on different cell configurations

### Why sf Works

R's **sf package** (Simple Features) uses **GEOS** with more robust algorithms:

1. **Proper topology handling** during buffer/union operations
2. **Automatic validity repair** after operations
3. **Consistent results** across complex geometries
4. **Better memory management** for large point clouds

### Benchmark Comparison

| Operation | shapely (Python) | sf (R) |
|-----------|------------------|--------|
| Buffer 50k points | Often fails | Succeeds |
| Union complex shapes | Self-intersections | Clean polygons |
| Shrink polygons | Invalid geometries | Valid geometries |
| Memory efficiency | Poor on large data | Better |

### Code Example (The Difference)

**Python shapely (FAILS on complex data)**:
```python
from shapely.ops import unary_union
from shapely.geometry import Point

# Create buffered points
buffered = [Point(x, y).buffer(50) for x, y in coords]
# Union - often creates invalid geometries
merged = unary_union(buffered)
# Shrink - frequently fails or produces self-intersections
shrunk = merged.buffer(-25)  # <-- PROBLEM: self-intersecting result
```

**R sf (WORKS correctly)**:
```r
library(sf)

# Create buffered points
points_sf <- st_as_sf(coords, coords = c("x", "y"))
buffered <- st_buffer(points_sf, dist = 50)

# Union - handles complex topologies correctly
merged <- st_union(buffered)

# Shrink - produces valid geometries
shrunk <- st_buffer(merged, dist = -25)  # <-- WORKS: clean result
```

---

## Core Components

### File Structure

```
spatial_domains/
├── functions/
│   ├── spatial_analysis.py      # Python entry point (make_spatial_domains)
│   ├── r_functions.R            # R implementation (MakeDomains, ReduceDomains)
│   └── r_bridge.py              # Generic R bridge utility
├── docs/
│   └── DOMAINS.md               # This file
```

### Key Functions

| Function | Language | Purpose | Location |
|----------|----------|---------|----------|
| `make_spatial_domains()` | Python | Entry point, orchestration | `spatial_analysis.py:4793` |
| `MakeDomains()` | R | Buffer-Union-Shrink algorithm | `r_functions.R:1228` |
| `ReduceDomains()` | R | Merge small domains | `r_functions.R:1425` |
| `calculate_domain_distances()` | Python | Distance computation | `spatial_analysis.py:5117` |
| `RBridge.call_r_function()` | Python | Generic R call utility | `r_bridge.py:39` |

### Dependencies

**Python**:
```
scanpy          # AnnData I/O
numpy           # Numerical operations
pandas          # DataFrame operations
scipy.spatial   # KD-Tree for distances
```

**R** (required):
```r
sf              # Simple Features for geometry operations
dplyr           # Data manipulation
purrr           # Functional programming
concaveman      # Concave hull algorithm
anndata         # Read h5ad files from R
reticulate      # Python interop (for anndata)
```

**Install R packages**:
```r
install.packages(c("sf", "concaveman", "dplyr", "purrr"))
# For anndata (requires Python environment):
install.packages("reticulate")
remotes::install_github("scverse/anndata")
```

---

## Algorithm Details

### Buffer-Union-Shrink Algorithm

The core algorithm creates smooth domain polygons from discrete cell positions:

```
┌─────────────────────────────────────────────────────────────────┐
│  STEP 1: BUFFER                                                 │
│                                                                 │
│  Each target cell → circle with radius = cell_dist (50 pixels)  │
│                                                                 │
│     ○  ○     ○        becomes        ◉  ◉     ◉                │
│       ○   ○                            ◉   ◉                    │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 2: UNION                                                  │
│                                                                 │
│  Merge overlapping circles into single polygon                  │
│                                                                 │
│     ◉◉◉                             ┌─────────────┐             │
│    ◉◉◉◉                             │             │             │
│       ◉◉◉     becomes               │    BLOB    │             │
│                                      │             │             │
│                                      └─────────────┘             │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 3: SHRINK                                                 │
│                                                                 │
│  Negative buffer by (cell_dist - shrink_margin) = 25 pixels    │
│                                                                 │
│     ┌─────────────┐                ┌───────────┐                │
│     │             │                │           │                │
│     │    BLOB    │   becomes      │  TIGHT   │                │
│     │             │                │           │                │
│     └─────────────┘                └───────────┘                │
│                                                                 │
│  This removes bridges and creates tighter polygons              │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│  STEP 4: CONCAVE HULL                                           │
│                                                                 │
│  Convert shrunk polygon to concave hull for cleaner boundary    │
│                                                                 │
│     Uses concaveman algorithm for natural-looking shapes        │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

### R Implementation (MakeDomains)

```r
MakeDomains <- function(input, group, group_subset, cell_dist = 225,
                        shrink_margin = 25, domain_prefix = NULL,
                        assign_all_cells = FALSE) {

  # 1. Filter to target cells
  target_cells <- object %>%
    filter(.[[group]] == group_subset)

  # 2. Buffer-Union-Shrink pipeline
  shrunken_polygons <- target_cells %>%
    st_as_sf(coords = c("x", "y")) %>%     # Convert to sf points
    st_buffer(dist = cell_dist) %>%         # Buffer each point
    st_union() %>%                          # Union all buffers
    st_cast("POLYGON") %>%                  # Extract individual polygons
    st_buffer(-(cell_dist - shrink_margin)) %>%  # Shrink
    st_cast("MULTIPOINT") %>%               # Get boundary points
    st_coordinates()                        # Extract coordinates

  # 3. Create concave hulls for each polygon
  for (i in unique(shrunken_polygons$L1)) {
    tmp <- shrunken_polygons %>% filter(L1 == i)
    conc <- concaveman::concaveman(tmp, length_threshold = 0, concavity = 0.999999)
    poly_list[[i]] <- conc
  }

  # 4. Assign cells to polygons
  if (assign_all_cells) {
    # Assign ALL cells (for heterogeneity analysis)
    cells_use <- object %>% select(x, y, cell, all_of(group))
  } else {
    # Only target cells
    cells_use <- target_cells %>% select(x, y, cell, all_of(group))
  }

  # 5. Spatial join: which cells fall within which polygons
  joined_data <- st_join(cells_use, polygon_data, join = st_intersects)

  return(list(cell_data = joined_data, polygon_data = polygon_data))
}
```

### ReduceDomains Algorithm

Merges small domains into neighboring larger domains:

```r
ReduceDomains <- function(result, min.cells.domain = 10, group = NULL,
                          group_subset = NULL) {

  # 1. Count TARGET cells per domain (not all cells)
  if (!is.null(group) && !is.null(group_subset)) {
    target_counts <- result$cell_data %>%
      filter(.data[[group]] == group_subset) %>%
      group_by(domain) %>%
      summarise(n_target = n())

    domains_merge <- target_counts$domain[target_counts$n_target <= min.cells.domain]
  }

  # 2. Find neighboring polygons for small domains
  neighbors <- st_intersects(small_polygons, polygon_data, sparse = FALSE)

  # 3. Merge small polygons into neighbors
  for (i in 1:nrow(small_polygons)) {
    neighbor_indices <- which(neighbors[i, ])
    if (!is_empty(neighbor_indices)) {
      # Merge geometry with first neighbor
      merged_polygon <- st_union(neighbor_polygon, small_polygon)
      # Update cell assignments
      result$cell_data$domain <- ifelse(domain == sp$domain, merged_polygon$domain, domain)
    } else {
      # No neighbors - drop this domain
      result$polygon_data <- result$polygon_data %>% filter(domain != sp$domain)
    }
  }

  return(result)
}
```

### Parameter Tuning Guide

| Parameter | Default | Effect | When to Adjust |
|-----------|---------|--------|----------------|
| `cell_dist` | 50 | Buffer radius (pixels) | Larger = more merging, smaller = distinct regions |
| `shrink_margin` | 25 | How much to subtract from cell_dist | Larger = tighter polygons, less bridging |
| `min_cells_domain` | 10 | Minimum target cells per domain | Increase for cleaner domains |
| `assign_all_cells` | True | Assign ALL cells or just targets | True for heterogeneity analysis |

---

## Python-R Integration

### How make_spatial_domains() Calls R

```python
@activity.defn
async def make_spatial_domains(data_uri, output_uri, group, group_subset, ...):
    """Python entry point that orchestrates R execution."""

    async def process_fn(input_data, params):
        # 1. Handle ontology_expression (Python-side)
        if params.get("ontology_expression"):
            adata = sc.read_h5ad(input_path)
            mask = evaluate_ontology_expression(params["ontology_expression"], adata)
            adata.obs["_ontology_filter"] = mask.astype(str)
            adata.write_h5ad(modified_input_path)
            effective_group = "_ontology_filter"
            effective_group_subset = "True"

        # 2. Generate R code dynamically
        r_code = f"""
        library(sf)
        library(anndata)
        source('{r_functions_path}')

        result <- MakeDomains(
            input = '{modified_input_path}',
            group = '{effective_group}',
            group_subset = '{group_subset_str}',
            cell_dist = {params['cell_dist']},
            shrink_margin = {params['shrink_margin']},
            domain_prefix = '{params['domain_prefix']}',
            assign_all_cells = {'TRUE' if params['assign_all_cells'] else 'FALSE'}
        )

        result <- ReduceDomains(
            result = result,
            min.cells.domain = {params['min_cells_domain']},
            group = '{effective_group}',
            group_subset = '{group_subset_str}'
        )

        write.csv(result$cell_data, '{result_csv_path}', row.names = FALSE)
        """

        # 3. Execute R via subprocess
        result = subprocess.run(
            ["Rscript", "-e", r_code],
            capture_output=True,
            text=True,
            timeout=1200
        )

        # 4. Read results and merge back to AnnData
        result_df = pd.read_csv(result_csv_path)
        domain_map = dict(zip(result_df["cell"], result_df["domain"]))
        adata.obs[params["output_column"]] = adata.obs.index.map(domain_map)

        return {"n_domains": adata.obs[output_column].nunique(), ...}
```

### RBridge (Generic R Bridge)

For more complex R integrations, use the generic RBridge class:

```python
# r_bridge.py
class RBridge:
    """Bridge for executing R functions from Python."""

    def __init__(self, r_functions_path: Optional[str] = None):
        self.r_functions_path = Path(r_functions_path or "r_functions.R")

    async def call_r_function(self, function_name: str, params: Dict,
                              timeout: int = 3600) -> Dict:
        """
        Call an R function with parameters and return the result.

        1. Write params to JSON file
        2. Generate R wrapper script
        3. Execute via Rscript subprocess
        4. Read JSON output
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            params_file = Path(tmpdir) / "params.json"
            output_file = Path(tmpdir) / "output.json"

            # Write params
            with open(params_file, "w") as f:
                json.dump(params, f)

            # Generate wrapper script
            wrapper_script = f"""
            library(jsonlite)
            source("{self.r_functions_path}")
            params <- fromJSON("{params_file}")
            result <- do.call("{function_name}", params)
            write_json(list(status="success", result=result), "{output_file}")
            """

            # Execute
            result = subprocess.run(
                ["Rscript", "-e", wrapper_script],
                capture_output=True,
                timeout=timeout
            )

            # Read output
            with open(output_file) as f:
                return json.load(f)
```

---

## Domain Distance Computation

### calculate_domain_distances() Overview

After creating domains, compute spatial distances between them:

```python
@activity.defn
async def calculate_domain_distances(
    data_uri: str,
    output_uri: str,
    source_domain_column: str,
    target_domain_column: str,
    source_domain_subset: Optional[List[str]] = None,
    target_domain_subset: Optional[List[str]] = None,
    distance_metric: str = "minimum",  # minimum | centroid | mean
    output_mode: str = "cell",         # cell | matrix | both
    output_distance_column: str = "distance_to_target",
    output_nearest_column: str = "nearest_target_domain",
) -> str:
```

### Distance Metrics

| Metric | Description | Use Case |
|--------|-------------|----------|
| `minimum` | Shortest cell-to-cell distance | Boundary-to-boundary proximity |
| `centroid` | Center-to-center distance | Domain-level relationships |
| `mean` | Average of all pairwise distances | Statistical analysis |

### KD-Tree Optimization

For `minimum` distance with per-cell annotation, we use KD-Tree for O(n log n) performance:

```python
from scipy.spatial import cKDTree

# Build KD-tree for all target cells (once)
target_mask = adata.obs[target_domain_column].isin(target_domains)
target_coords = adata.obsm["spatial"][target_mask.values]
target_tree = cKDTree(target_coords)

# Query all source cells at once
source_mask = adata.obs[source_domain_column].isin(source_domains)
source_coords = adata.obsm["spatial"][source_mask.values]

# Single vectorized query - O(n log n) vs O(n^2)
distances, nearest_idx = target_tree.query(source_coords, k=1)

# Map to domain names
target_domains_array = adata.obs[target_domain_column].iloc[target_indices].values
nearest_domains = target_domains_array[nearest_idx]

# Assign to adata.obs
adata.obs.iloc[source_indices, dist_col] = distances
adata.obs.iloc[source_indices, nearest_col] = nearest_domains
```

### Output Modes

**Cell Mode** (`output_mode="cell"`):
```python
# Per-cell columns added to adata.obs
adata.obs["distance_to_target"]      # Float: distance in coordinate units
adata.obs["nearest_target_domain"]   # String: name of nearest domain
```

**Matrix Mode** (`output_mode="matrix"`):
```python
# Distance matrix stored in adata.uns
adata.uns["domain_distances"] = {
    "distance_matrix": {
        "MS4A1_1": {"CMS2_5": 123.4, "CMS2_6": 456.7},
        "MS4A1_2": {"CMS2_5": 234.5, "CMS2_6": 567.8},
    },
    "summary_statistics": {
        "min_distance": 123.4,
        "max_distance": 567.8,
        "mean_distance": 345.6,
        "median_distance": 345.1,
    }
}
```

---

## Standalone Extraction Guide

### Minimal Standalone Package

```
spatial_domains/
├── __init__.py
├── make_domains.py          # Python entry point
├── distance.py              # Domain distance computation
├── r_functions.R            # R implementation
├── r_bridge.py              # R bridge utility
└── setup.py
```

### __init__.py

```python
from .make_domains import make_spatial_domains
from .distance import calculate_domain_distances
from .r_bridge import RBridge, call_r_function

__all__ = [
    'make_spatial_domains',
    'calculate_domain_distances',
    'RBridge',
    'call_r_function',
]
```

### make_domains.py (Standalone)

```python
"""
Standalone spatial domain creation using R's sf package.
"""
import subprocess
import pandas as pd
from pathlib import Path
from typing import Optional, Any
import tempfile


def make_spatial_domains(
    h5ad_path: str,
    output_path: str,
    group: str,
    group_subset: Any,
    cell_dist: float = 50.0,
    shrink_margin: float = 25.0,
    domain_prefix: str = "domain",
    min_cells_domain: int = 10,
    output_column: str = "spatial_domain",
    assign_all_cells: bool = True,
    r_functions_path: Optional[str] = None,
) -> dict:
    """
    Create spatial domains from cell groupings using R's sf package.

    Args:
        h5ad_path: Path to input h5ad file with spatial coordinates
        output_path: Path for output h5ad with domain assignments
        group: Column name to filter by (e.g., 'cell_type', 'leiden')
        group_subset: Value to filter in the group column
        cell_dist: Buffer distance in pixels (default 50)
        shrink_margin: Margin to subtract when shrinking (default 25)
        domain_prefix: Prefix for domain names
        min_cells_domain: Minimum target cells per domain (default 10)
        output_column: Output column name (default 'spatial_domain')
        assign_all_cells: Assign ALL cells to domains (default True)
        r_functions_path: Path to r_functions.R

    Returns:
        Dict with n_domains, n_cells_assigned, output_column
    """
    import scanpy as sc
    import shutil

    # Check R availability
    if not shutil.which("Rscript"):
        raise RuntimeError(
            "Rscript not found. Install R and packages: sf, concaveman, dplyr, purrr, anndata"
        )

    # Determine R functions path
    if r_functions_path is None:
        r_functions_path = Path(__file__).parent / "r_functions.R"

    # Result file
    result_csv = tempfile.mktemp(suffix=".csv")

    # Handle group_subset type
    group_subset_str = str(group_subset).replace("'", "\\'")

    # Generate R code
    r_code = f"""
library(dplyr)
library(sf)
library(purrr)
library(concaveman)

# Configure reticulate
library(reticulate)
python_path <- Sys.which("python3")
if (nchar(python_path) > 0) {{
    use_python(python_path, required = TRUE)
}}

library(anndata)
source('{r_functions_path}')

result <- MakeDomains(
    input = '{h5ad_path}',
    group = '{group}',
    group_subset = '{group_subset_str}',
    cell_dist = {cell_dist},
    shrink_margin = {shrink_margin},
    domain_prefix = '{domain_prefix}',
    assign_all_cells = {'TRUE' if assign_all_cells else 'FALSE'}
)

result <- ReduceDomains(
    result = result,
    min.cells.domain = {min_cells_domain},
    group = '{group}',
    group_subset = '{group_subset_str}'
)

write.csv(result$cell_data, '{result_csv}', row.names = FALSE)
"""

    try:
        # Run R
        result = subprocess.run(
            ["Rscript", "-e", r_code],
            capture_output=True,
            text=True,
            timeout=1200
        )

        if result.returncode != 0:
            raise RuntimeError(f"R failed: {result.stderr}")

        # Read results
        result_df = pd.read_csv(result_csv)

        # Load adata and add domain column
        adata = sc.read_h5ad(h5ad_path)
        result_df["cell"] = result_df["cell"].astype(str)
        domain_map = dict(zip(result_df["cell"], result_df["domain"]))
        adata.obs[output_column] = adata.obs.index.map(domain_map)

        # Save
        adata.write_h5ad(output_path)

        return {
            "n_domains": int(adata.obs[output_column].nunique()),
            "n_cells_assigned": int(adata.obs[output_column].notna().sum()),
            "output_column": output_column,
        }

    finally:
        import os
        if os.path.exists(result_csv):
            os.remove(result_csv)
```

### distance.py (Standalone)

```python
"""
Standalone domain distance computation.
"""
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from typing import List, Optional, Dict, Any
import scanpy as sc


def calculate_domain_distances(
    h5ad_path: str,
    output_path: str,
    source_domain_column: str,
    target_domain_column: str,
    source_domain_subset: Optional[List[str]] = None,
    target_domain_subset: Optional[List[str]] = None,
    distance_metric: str = "minimum",
    output_mode: str = "both",
    output_distance_column: str = "distance_to_target",
    output_nearest_column: str = "nearest_target_domain",
) -> Dict[str, Any]:
    """
    Calculate spatial distances between domains.

    Args:
        h5ad_path: Path to h5ad with spatial domains
        output_path: Path for output h5ad with distances
        source_domain_column: Column with source domain labels
        target_domain_column: Column with target domain labels
        source_domain_subset: Specific source domains to measure from
        target_domain_subset: Specific target domains to measure to
        distance_metric: 'minimum', 'centroid', or 'mean'
        output_mode: 'cell', 'matrix', or 'both'
        output_distance_column: Name for per-cell distance column
        output_nearest_column: Name for nearest domain column

    Returns:
        Dict with distance matrix and summary statistics
    """
    adata = sc.read_h5ad(h5ad_path)

    # Validate columns
    if source_domain_column not in adata.obs.columns:
        raise ValueError(f"Source column '{source_domain_column}' not found")
    if "spatial" not in adata.obsm:
        raise ValueError("Spatial coordinates not found in adata.obsm['spatial']")

    # Get domain lists
    source_domains = adata.obs[source_domain_column].dropna().unique().tolist()
    target_domains = adata.obs[target_domain_column].dropna().unique().tolist()

    if source_domain_subset:
        source_domains = [d for d in source_domains if d in source_domain_subset]
    if target_domain_subset:
        target_domains = [d for d in target_domains if d in target_domain_subset]

    # Initialize distance matrix
    distance_matrix = pd.DataFrame(
        index=source_domains, columns=target_domains, dtype=float
    )

    # Initialize cell columns
    if output_mode in ["cell", "both"]:
        adata.obs[output_distance_column] = np.nan
        adata.obs[output_nearest_column] = None

    # Calculate distances
    if distance_metric == "minimum" and output_mode in ["cell", "both"]:
        # KD-Tree optimized approach
        target_mask = adata.obs[target_domain_column].isin(target_domains)
        target_indices = np.where(target_mask.values)[0]
        target_coords = adata.obsm["spatial"][target_indices]
        target_tree = cKDTree(target_coords)
        target_domains_arr = adata.obs[target_domain_column].iloc[target_indices].values

        source_mask = adata.obs[source_domain_column].isin(source_domains)
        source_indices = np.where(source_mask.values)[0]
        source_coords = adata.obsm["spatial"][source_indices]

        distances, nearest_idx = target_tree.query(source_coords, k=1)
        nearest_domains = target_domains_arr[nearest_idx]

        # Assign to adata
        adata.obs.iloc[source_indices, adata.obs.columns.get_loc(output_distance_column)] = distances
        adata.obs.iloc[source_indices, adata.obs.columns.get_loc(output_nearest_column)] = nearest_domains

        # Build matrix from per-cell results
        source_domains_arr = adata.obs[source_domain_column].iloc[source_indices].values
        for src in source_domains:
            src_mask = source_domains_arr == src
            if not src_mask.any():
                continue
            for tgt in target_domains:
                if src == tgt and source_domain_column == target_domain_column:
                    distance_matrix.loc[src, tgt] = 0.0
                    continue
                tgt_mask = nearest_domains[src_mask] == tgt
                if tgt_mask.any():
                    distance_matrix.loc[src, tgt] = distances[src_mask][tgt_mask].min()

    elif distance_metric == "centroid":
        # Centroid distances
        source_centroids = {}
        target_centroids = {}

        for src in source_domains:
            mask = adata.obs[source_domain_column] == src
            coords = adata.obsm["spatial"][mask.values]
            if len(coords) > 0:
                source_centroids[src] = coords.mean(axis=0)

        for tgt in target_domains:
            mask = adata.obs[target_domain_column] == tgt
            coords = adata.obsm["spatial"][mask.values]
            if len(coords) > 0:
                target_centroids[tgt] = coords.mean(axis=0)

        for src in source_domains:
            if src not in source_centroids:
                continue
            for tgt in target_domains:
                if src == tgt:
                    distance_matrix.loc[src, tgt] = 0.0
                    continue
                if tgt not in target_centroids:
                    continue
                dist = np.linalg.norm(source_centroids[src] - target_centroids[tgt])
                distance_matrix.loc[src, tgt] = dist

    # Summary statistics
    valid = distance_matrix.values[~np.isnan(distance_matrix.values)]
    summary = {
        "min_distance": float(valid.min()) if len(valid) > 0 else None,
        "max_distance": float(valid.max()) if len(valid) > 0 else None,
        "mean_distance": float(valid.mean()) if len(valid) > 0 else None,
        "median_distance": float(np.median(valid)) if len(valid) > 0 else None,
    }

    # Store in adata.uns
    result_metadata = {
        "source_domain_column": source_domain_column,
        "target_domain_column": target_domain_column,
        "distance_metric": distance_metric,
        "source_domains": source_domains,
        "target_domains": target_domains,
        "summary_statistics": summary,
    }

    if output_mode in ["matrix", "both"]:
        result_metadata["distance_matrix"] = distance_matrix.to_dict(orient="index")
        adata.uns["domain_distances"] = result_metadata

    # Save
    adata.write_h5ad(output_path)

    return result_metadata
```

### setup.py

```python
from setuptools import setup, find_packages

setup(
    name="spatial-domains",
    version="1.0.0",
    description="Spatial domain creation and distance computation for spatial biology",
    packages=find_packages(),
    package_data={
        "spatial_domains": ["r_functions.R"],
    },
    install_requires=[
        "scanpy>=1.9.0",
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
    ],
    extras_require={
        "full": [
            "anndata>=0.8.0",
        ],
    },
    python_requires=">=3.8",
)
```

---

## API Reference

### make_spatial_domains()

```python
async def make_spatial_domains(
    data_uri: str,
    output_uri: str,
    group: Optional[str] = None,
    group_subset: Optional[Any] = None,
    ontology_expression: Optional[str] = None,
    cell_dist: float = 50.0,
    shrink_margin: float = 25.0,
    domain_prefix: Optional[str] = None,
    min_cells_domain: int = 10,
    output_column: str = "spatial_domain",
    assign_all_cells: bool = True,
    workflow_focus: Optional[str] = None,
    workflow_stage: Optional[str] = None,
) -> str:
```

**Parameters**:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_uri` | str | required | GCS URI to spatial data |
| `output_uri` | str | required | GCS URI for output |
| `group` | str | None | Column name to filter by |
| `group_subset` | Any | None | Value to filter in group column |
| `ontology_expression` | str | None | Boolean expression (e.g., "CL:0000236 & NCIT:C4349") |
| `cell_dist` | float | 50.0 | Buffer distance in pixels |
| `shrink_margin` | float | 25.0 | Margin for shrinking |
| `domain_prefix` | str | None | Prefix for domain names |
| `min_cells_domain` | int | 10 | Minimum target cells per domain |
| `output_column` | str | "spatial_domain" | Output column name |
| `assign_all_cells` | bool | True | Assign ALL cells to domains |

**Returns**: Output URI with domain assignments

### calculate_domain_distances()

```python
async def calculate_domain_distances(
    data_uri: str,
    output_uri: str,
    source_domain_column: str,
    target_domain_column: str,
    source_domain_subset: Optional[List[str]] = None,
    target_domain_subset: Optional[List[str]] = None,
    distance_metric: str = "minimum",
    output_mode: str = "cell",
    annotate_cells: bool = True,
    output_distance_column: str = "distance_to_target",
    output_nearest_column: str = "nearest_target_domain",
) -> str:
```

**Parameters**:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `source_domain_column` | str | required | Column with source domains |
| `target_domain_column` | str | required | Column with target domains |
| `source_domain_subset` | List[str] | None | Filter source domains |
| `target_domain_subset` | List[str] | None | Filter target domains |
| `distance_metric` | str | "minimum" | minimum, centroid, or mean |
| `output_mode` | str | "cell" | cell, matrix, or both |
| `output_distance_column` | str | "distance_to_target" | Per-cell distance column |
| `output_nearest_column` | str | "nearest_target_domain" | Per-cell nearest domain column |

**Returns**: Output URI with distance annotations

---

## Usage Examples

### Basic Domain Creation

```python
# Create B cell domains
result = await make_spatial_domains(
    data_uri="gs://bucket/annotated.h5ad",
    output_uri="gs://bucket/bcell_domains.h5ad",
    group="cell_type_ontology_id",
    group_subset="CL:0000236",  # B cell
    cell_dist=50,
    shrink_margin=25,
    domain_prefix="Bcell",
    min_cells_domain=10,
    assign_all_cells=True
)
print(f"Created {result['n_domains']} B cell domains")
```

### Using Ontology Expressions

```python
# B cells in tumor regions (compound filter)
result = await make_spatial_domains(
    data_uri="gs://bucket/annotated.h5ad",
    output_uri="gs://bucket/bcell_tumor.h5ad",
    ontology_expression="CL:0000236 & NCIT:C4349",  # B cell AND tumor
    cell_dist=75,
    domain_prefix="Bcell_Tumor"
)
```

### Cross-Domain Distance Calculation

```python
# Distance from B cell domains to tumor domains
result = await calculate_domain_distances(
    data_uri="gs://bucket/both_domains.h5ad",
    output_uri="gs://bucket/distances.h5ad",
    source_domain_column="spatial_domain",
    target_domain_column="spatial_domain",
    source_domain_subset=["Bcell_1", "Bcell_2", "Bcell_3"],
    target_domain_subset=["Tumor_1"],
    distance_metric="minimum",
    output_mode="both",
    output_distance_column="dist_to_tumor",
    output_nearest_column="nearest_tumor"
)

print(f"Distance summary: {result['summary_statistics']}")
```

### Full Workflow Example

```python
# Step 1: Create B cell domains
bcell_result = await make_spatial_domains(
    data_uri="gs://bucket/step12.h5ad",
    output_uri="gs://bucket/step15_bcell.h5ad",
    group="cell_type_ontology_id",
    group_subset="CL:0000236",
    domain_prefix="MS4A1",
    output_column="bcell_domain"
)

# Step 2: Create tumor domains
tumor_result = await make_spatial_domains(
    data_uri="gs://bucket/step15_bcell.h5ad",
    output_uri="gs://bucket/step16_tumor.h5ad",
    ontology_expression="NCIT:C4910",  # CMS2 tumor
    domain_prefix="CMS2",
    output_column="tumor_domain"
)

# Step 3: Calculate B cell → Tumor distances
dist_result = await calculate_domain_distances(
    data_uri="gs://bucket/step16_tumor.h5ad",
    output_uri="gs://bucket/step17_distances.h5ad",
    source_domain_column="bcell_domain",
    target_domain_column="tumor_domain",
    distance_metric="minimum",
    output_mode="both"
)

# Step 4: Find markers in closest B cell domains
# (Use select_markers_by_proximity)
```

---

## Troubleshooting

### Common Issues

#### "Rscript not found"

**Cause**: R not installed or not in PATH

**Solution**:
```bash
# Install R
sudo apt-get install r-base

# Or on macOS
brew install r

# Verify
which Rscript
```

#### "R package 'sf' not found"

**Cause**: Missing R dependencies

**Solution**:
```r
install.packages(c("sf", "concaveman", "dplyr", "purrr"))
remotes::install_github("scverse/anndata")
```

#### "No cells match group_subset"

**Cause**: Type mismatch between column values and filter

**Solution**:
- Check column type: `adata.obs[group].dtype`
- Ensure group_subset matches (string vs int vs bool)
- Use ontology_expression for complex filters

#### "Invalid geometry" warnings in R

**Cause**: Complex point configurations

**Solution**:
- Increase `cell_dist` (larger buffer)
- Increase `shrink_margin` (less aggressive shrink)
- Check for outlier cells far from main clusters

#### KD-Tree memory error

**Cause**: Too many cells for available memory

**Solution**:
- Subset to region of interest first
- Use `output_mode="matrix"` (less memory)
- Increase available memory

### Performance Tips

1. **Large datasets**: Use `output_mode="matrix"` to skip per-cell annotation
2. **Many domains**: Filter to specific domain subsets
3. **Memory**: Process samples individually, not merged
4. **Speed**: Pre-compute domains in parallel, then merge for distance

---

## References

### R Packages

- **sf (Simple Features)**: https://r-spatial.github.io/sf/
- **concaveman**: https://github.com/joelgombin/concaveman
- **anndata (R)**: https://github.com/scverse/anndata

### Python Libraries

- **scanpy**: https://scanpy.readthedocs.io/
- **scipy.spatial.cKDTree**: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.cKDTree.html

### Algorithms

- **Buffer-Union-Shrink**: Common GIS technique for region delineation
- **Concave Hull**: Park & Oh (2012) "A New Concave Hull Algorithm"
- **KD-Tree**: Bentley (1975) "Multidimensional Binary Search Trees"

### Internal References

- Temporal activity pattern: `docs/PROCESS_ACTIVITY_MIGRATION_GUIDE.md`
- Ontology expressions: `functions/annotation_utils.py`

---

## Changelog

### 2025-12-16
- Initial comprehensive documentation
- Standalone extraction guide
- Algorithm details with diagrams
- KD-Tree optimization documentation
- R vs Python comparison

### Prior Work
- Original R implementation from lab protocols
- Python-R bridge development
- KD-Tree optimization for large datasets
