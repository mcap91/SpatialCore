# Using SpatialCore from R

SpatialCore is a Python library for spatial biology analysis. This guide shows
how R users can leverage the full SpatialCore library from within R using the
reticulate package.

## Why Use SpatialCore from R?

- **Access all SpatialCore functions** without learning Python syntax
- **In-memory conversion** - no file I/O overhead
- **Round-trip workflows** - analyze in Python, visualize in R
- **Works with existing Seurat pipelines**

## Quick Start

```r
library(Seurat)
library(reticulate)

# Setup (once per session)
use_condaenv("spatialcore")

# Source conversion functions
source("https://raw.githubusercontent.com/mcap91/SpatialCore/main/docs/r_bridge/spatialcore_convert.R")

# Load your Seurat object
seurat_obj <- readRDS("my_spatial_data.rds")

# Convert to AnnData
adata <- seurat_to_adata(seurat_obj)

# Use ANY SpatialCore function
spc <- import("spatialcore")
spc$spatial$compute_morans_i(adata, genes = c("EPCAM", "CD3E"))
spc$clustering$leiden(adata, resolution = 0.5)

# Convert back with new annotations
seurat_obj <- adata_to_seurat(adata, seurat_obj)

# Continue with Seurat visualization
DimPlot(seurat_obj, group.by = "leiden")
```

## Installation

### Step 1: Install R Packages

```r
install.packages(c("reticulate", "Matrix", "Seurat"))
```

### Step 2: Create Python Environment

In your terminal (not R):

```bash
# Using mamba (recommended)
mamba create -n spatialcore python=3.10
mamba activate spatialcore
pip install spatialcore scanpy anndata

# Or using conda
conda create -n spatialcore python=3.10
conda activate spatialcore
pip install spatialcore scanpy anndata
```

### Step 3: Configure R

```r
library(reticulate)

# Point to your conda environment
use_condaenv("spatialcore", required = TRUE)

# Verify
py_config()
```

### Step 4: Get Conversion Functions

Option A: Source from GitHub (recommended)

```r
source("https://raw.githubusercontent.com/mcap91/SpatialCore/main/docs/r_bridge/spatialcore_convert.R")
```

Option B: Download locally

Download `spatialcore_convert.R` from [GitHub](https://github.com/mcap91/SpatialCore/blob/main/docs/r_bridge/spatialcore_convert.R) and source it:

```r
source("path/to/spatialcore_convert.R")
```

## Conversion Functions

### seurat_to_adata()

Converts a Seurat object to AnnData.

```r
adata <- seurat_to_adata(
    seurat_obj,
    assay = "RNA",           # Which assay to convert
    x_slot = "data",         # "data" (normalized) or "counts" (raw) for X
    include_counts = TRUE,   # Include raw counts in layers["counts"]
    include_scale = FALSE,   # Include scaled data (large!)
    include_reductions = TRUE,
    include_graphs = TRUE,
    spatial_key = NULL       # Auto-detect spatial coordinates
)
```

**Default behavior:**

- `X` contains **normalized data** (from `data` slot)
- Raw counts are stored in `layers["counts"]`
- `uns["X_slot"]` records what X contains for safe round-trip

### adata_to_seurat()

Converts AnnData back to Seurat.

```r
# Create new Seurat object
seurat_new <- adata_to_seurat(adata)

# Or merge into existing object (preserves additional Seurat data)
seurat_updated <- adata_to_seurat(adata, seurat_obj = original_seurat)
```

**Merge behavior:**

- New metadata columns are added
- Existing columns are updated
- New reductions are added
- Existing reductions are replaced

## Data Mapping

### Seurat to AnnData

| Seurat | AnnData | Notes |
|--------|---------|-------|
| `GetAssayData(slot="data")` | `X` | Normalized (default) |
| `GetAssayData(slot="counts")` | `layers["counts"]` | Raw counts |
| `seurat_obj[[]]` | `obs` | Cell metadata |
| `GetAssay()[[]]` | `var` | Gene metadata |
| `@reductions$pca` | `obsm["X_pca"]` | PCA |
| `@reductions$umap` | `obsm["X_umap"]` | UMAP |
| `@reductions$spatial` | `obsm["spatial"]` | Coordinates |
| `@reductions$pca@feature.loadings` | `varm["PCs"]` | Gene loadings |
| `@graphs` | `obsp` | Neighbor graphs |

### AnnData to Seurat

| AnnData | Seurat | Notes |
|---------|--------|-------|
| `X` | `data` or `counts` | Based on `uns["X_slot"]` |
| `layers["counts"]` | `counts` slot | Preferred |
| `obs` | `seurat_obj[[]]` | Cell metadata |
| `obsm["X_pca"]` | `@reductions$pca` | With `PC_` key |
| `obsm["X_umap"]` | `@reductions$umap` | With `UMAP_` key |
| `obsm["spatial"]` | `@reductions$spatial` | With `Spatial_` key |

## Example Workflows

### Spatial Autocorrelation Analysis

```r
library(Seurat)
library(reticulate)

use_condaenv("spatialcore")
source("spatialcore_convert.R")

# Load spatial data
seurat_obj <- readRDS("xenium_sample.rds")

# Convert
adata <- seurat_to_adata(seurat_obj)

# Compute Moran's I for spatial autocorrelation
spc <- import("spatialcore")
genes <- c("EPCAM", "CD3E", "CD68", "COL1A1", "ACTA2")
spc$spatial$compute_morans_i(adata, genes = genes)

# Get results
morans_df <- py_to_r(adata$uns["moranI"])
print(morans_df)

# High Moran's I = spatially clustered expression
```

### Leiden Clustering on Spatial Graph

```r
# Convert
adata <- seurat_to_adata(seurat_obj)

# Build spatial neighbors (if not already present)
spc <- import("spatialcore")
# spc$spatial$compute_spatial_neighbors(adata, n_neighbors = 15L)

# Leiden clustering
spc$clustering$leiden(adata, resolution = 0.8)

# Get results back
seurat_obj <- adata_to_seurat(adata, seurat_obj)

# Visualize
DimPlot(seurat_obj, reduction = "spatial", group.by = "leiden", pt.size = 0.5)
```

### Cell Type Annotation with CellTypist

```r
# Convert
adata <- seurat_to_adata(seurat_obj)

# Run CellTypist
spc <- import("spatialcore")
spc$annotation$run_celltypist(
    adata,
    model = "Immune_All_Low.pkl",
    majority_voting = TRUE
)

# Results in adata$obs["predicted_labels"]
seurat_obj <- adata_to_seurat(adata, seurat_obj)

# Plot cell types
DimPlot(seurat_obj, reduction = "spatial", group.by = "predicted_labels")
```

### Batch Processing Multiple Samples

```r
library(purrr)

# List of Seurat objects
samples <- list(
    sample1 = readRDS("sample1.rds"),
    sample2 = readRDS("sample2.rds"),
    sample3 = readRDS("sample3.rds")
)

# Process each sample
processed <- map(samples, function(seurat_obj) {
    # Convert
    adata <- seurat_to_adata(seurat_obj)

    # Run analysis
    spc$spatial$compute_morans_i(adata)
    spc$clustering$leiden(adata, resolution = 0.5)

    # Convert back
    adata_to_seurat(adata, seurat_obj)
})

# Combine results
combined_metadata <- map_dfr(processed, ~ .x[[]], .id = "sample")
```

## Troubleshooting

### "Python is not available"

```r
# Check Python configuration
py_available()

# List conda environments
conda_list()

# Set the correct environment
use_condaenv("spatialcore", required = TRUE)

# If still failing, try specifying the Python path directly
use_python("/path/to/miniforge3/envs/spatialcore/bin/python")
```

### "Module 'spatialcore' not found"

```bash
# Activate your conda environment
mamba activate spatialcore

# Install spatialcore
pip install spatialcore
```

### "BPCells on-disk matrices are not yet supported"

Seurat v5 can use disk-backed matrices via BPCells. These are not yet supported.
Convert to in-memory first:

```r
# Convert BPCells layer to standard matrix
seurat_obj[["RNA"]]@layers$counts <- as.matrix(seurat_obj[["RNA"]]@layers$counts)
```

### "Multiple graphs match _nn or _snn pattern"

Your Seurat object has multiple neighbor graphs. Specify which to use:

```r
adata <- seurat_to_adata(
    seurat_obj,
    graph_key_map = c(
        "RNA_nn" = "connectivities",
        "RNA_snn" = "distances"
    )
)
```

### "Cell names do not match"

When merging back to Seurat, cell names must match exactly. If cells were filtered
or reordered in Python, create a new Seurat object instead:

```r
# Don't pass seurat_obj - create new one
seurat_new <- adata_to_seurat(adata)
```

### Slow Conversion

For large objects:

1. Disable scaled data: `include_scale = FALSE`
2. Disable graphs if not needed: `include_graphs = FALSE`
3. Use sparse matrices (default for most Seurat objects)

### Memory Issues

For very large datasets:

```r
# Convert only essential data
adata <- seurat_to_adata(
    seurat_obj,
    include_counts = FALSE,  # Skip if X has what you need
    include_scale = FALSE,
    include_graphs = FALSE,
    include_reductions = FALSE
)
```

## Seurat Version Compatibility

Both Seurat v4 and v5 are supported:

| Feature | Seurat v4 | Seurat v5 |
|---------|-----------|-----------|
| Assay access | `slot=` | `layer=` |
| BPCells | Not available | **Not supported** (convert to in-memory) |
| Standard assays | Supported | Supported |

The conversion functions automatically detect and handle version differences.

## Limitations

- **BPCells**: On-disk matrices not yet supported (future SpatialCore feature)
- **Multi-assay**: Only one assay converted at a time
- **Images**: Seurat image slots not converted (spatial coords are)
- **Large scale.data**: Must cover all features to be included

## API Reference

### seurat_to_adata()

```r
seurat_to_adata(
    seurat_obj,
    assay = NULL,           # Default: DefaultAssay(seurat_obj)
    x_slot = "data",        # "data" or "counts"
    include_counts = TRUE,
    include_scale = FALSE,
    include_reductions = TRUE,
    include_graphs = TRUE,
    graph_key_map = NULL,   # Named vector for explicit graph mapping
    spatial_key = NULL      # Override spatial reduction name
)
```

### adata_to_seurat()

```r
adata_to_seurat(
    adata,
    seurat_obj = NULL,      # NULL = create new, otherwise merge
    assay = "RNA",
    project = "SpatialCore"
)
```

### setup_spatialcore()

```r
setup_spatialcore(
    conda_env = "spatialcore",
    verbose = TRUE
)
```

### is_spatialcore_ready()
```r
is_spatialcore_ready(conda_env = "spatialcore")
# Returns TRUE/FALSE
```

### get_spatialcore()

```r
spc <- get_spatialcore()
# Returns imported spatialcore module
```
