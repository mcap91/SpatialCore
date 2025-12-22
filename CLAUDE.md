# SpatialCore - Claude Instructions

## Project Overview

SpatialCore is a Python/R package for spatial biology analysis. It provides standardized, robust wrappers around existing statistical tools for computational biologists working with spatial transcriptomics data (Xenium, CosMx, Visium).

**Philosophy:** Not novel math—just better, standardized engineering for the community.

## Branch Strategy

- **`main` branch (public):** Published package via pip - spatial statistics, clustering, NMF, annotation tools
- **`dev` branch (internal):** Preprocessing pipelines, data ingestion, QC workflows - never published

## Environment

- Use the `spatialcore` mamba environment for all Python/R operations
- Primary data format: AnnData (.h5ad)
- Secondary data format: SeuratObject (.Rds)
- Supported platforms: CosMx, Xenium, Visium

## Coding Conventions

### Naming
- Use `snake_case` for functions, variables, and files (NOT dot separators)
- Each function should do ONE specific thing with clear input/output

### Required for All Functions
1. **Docstrings** (Python) or **roxygen2** (R) with references where possible
2. **Type hints** (Python) or parameter typing (R)
3. **Input validation** with descriptive error messages
4. **Logging** - use clear logging statements via `spatialcore.core.logging`
5. **Metadata tracking** - update metadata via `spatialcore.core.metadata`

### scanpy Conventions
- Modify AnnData in-place by default, offer `copy=True` option
- Store results in appropriate slots:
  - `adata.obs` - cell-level annotations
  - `adata.var` - gene-level annotations
  - `adata.obsm` - cell embeddings (PCA, UMAP, etc.)
  - `adata.uns` - unstructured data (parameters, etc.)

## When Creating Functions

1. **Always ask** which platform(s) the function should support
2. **Prefer vectorized operations** over loops
3. **Use existing packages** as building blocks:
   - Python: scanpy, squidpy, spatialdata
   - R: Seurat
4. **Cache intermediate results** to `.cache/` directory as `.h5ad` files

## Parallelization Strategy

- **Default:** Sequential sample processing (memory-safe for large spatial datasets)
- **Within-function:** Use numba/parallel operations where available (scanpy defaults)
- **Optional:** Parallel samples via multiprocessing if memory permits

## R Integration

For analyses requiring R/Seurat:

```python
# rpy2 bridge - call R functions directly from Python
from spatialcore.r_bridge import run_seurat_analysis
adata = run_seurat_analysis(adata, method="sctransform")
```

Alternative: h5ad interchange (save AnnData, load in R, process, save back)

## Project Structure (main branch)

```
SpatialCore/
├── src/
│   └── spatialcore/
│       ├── __init__.py
│       ├── core/              # Logging, metadata, caching utilities
│       ├── spatial/           # Moran's I, Lee's L, spatial autocorrelation
│       ├── clustering/        # Domains, neighborhoods, communities
│       ├── nmf/               # Spatial NMF
│       ├── diffusion/         # Diffusion maps, pseudotime
│       ├── ontology/          # Ontology conversion tools
│       ├── annotation/        # CellTypist wrappers, benchmarking
│       └── r_bridge/          # R/Python integration
├── R/                         # R package source
├── tests/
│   ├── tests.yaml             # YAML-based test definitions
│   ├── fixtures/              # Test data files
│   ├── conftest.py            # pytest fixtures
│   └── test_runner.py         # Cross-language test runner
├── benchmarks/                # Benchmarking scripts
├── examples/
├── docs/
└── .cache/                    # Intermediate results (gitignored)
```

## Module Responsibilities (main branch)

| Module | Purpose |
|--------|---------|
| `core` | Logging, metadata tracking, caching utilities |
| `spatial` | Moran's I, Lee's L (HH/LL/HL/LH), spatial autocorrelation |
| `clustering` | Neighborhoods, domains, niches, community detection |
| `nmf` | Spatial non-negative matrix factorization |
| `diffusion` | Diffusion maps, pseudotime analysis |
| `ontology` | Cell ontology conversion and mapping |
| `annotation` | CellTypist wrappers, custom model training, benchmarking |
| `r_bridge` | rpy2 wrappers, Seurat integration |

## Testing

Tests use a YAML-based framework for cross-language validation:

```yaml
tests:
  - language: python
    module: spatial
    function: morans_i
    inputs:
      - adata: "fixtures/test_spatial.h5ad"
        gene: "CD8A"
    expected_outputs:
      - I: 0.45
        p_value: 0.001
    tolerance: 0.01
```

Run tests:
- Python: `pytest tests/`
- R: `testthat::test_dir("tests/")`
- Cross-language: `python tests/test_runner.py`

## Key Dependencies

### Python
- scanpy, squidpy, spatialdata
- anndata, pandas, numpy, scipy
- celltypist (optional, for annotation)
- rpy2 (optional, for R integration)
- numba (optional, for parallelization)

### R
- Seurat
- anndata (R package)
- testthat

## Terminology Standards

Use these definitions consistently:
- **Neighborhood:** Immediate spatial context of a cell
- **Domain:** Larger tissue region with coherent cell composition
- **Niche:** Functional microenvironment
- **Community:** Graph-based cell grouping
