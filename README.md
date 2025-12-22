# SpatialCore

**Standardized spatial statistics for computational biology.**

License & Commercial Use SpatialCore is licensed under the Apache License 2.0. This means you are free to use, modify, and distribute this software, including for commercial products. However, the SpatialCore name and trademarks are reserved to ensure the community can rely on the "Standardized" quality of the core library.

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

---

## The Problem - NEEDS UPDATE

We spend too much time re-implementing basic spatial statistics.

## The Solution - NEEDS UPDATE

A thin, robust wrapper around standard libraries to ensure a Python user and an R user get the exact same result for the same biological question.

## The Goal - NEEDS UPDATE

Not novel math—just better, standardized engineering for the community.

---

## Installation

```bash
pip install spatialcore
```

With optional dependencies:
```bash
# CellTypist integration
pip install spatialcore[celltypist]

# R integration via rpy2
pip install spatialcore[r]

# All optional dependencies
pip install spatialcore[all]
```

## Quick Start

```python
import scanpy as sc
import spatialcore as sp

# Load your spatial data (Xenium, CosMx, Visium, etc.)
adata = sc.read_h5ad("spatial_data.h5ad")

# Spatial autocorrelation
sp.spatial.morans_i(adata, gene="CD8A", spatial_key="spatial")
sp.spatial.lees_l(adata, gene_x="CD8A", gene_y="CD4", spatial_key="spatial")

# Domain identification
sp.clustering.compute_domains(adata, method="leiden", resolution=0.5)
sp.clustering.domain_distance(adata, domain_key="domain")

# Spatial NMF
sp.nmf.spatial_nmf(adata, n_components=10)
```

## Features

| Module | Description |
|--------|-------------|
| `spatialcore.spatial` | Moran's I, Lee's L with HH/LL/HL/LH classification |
| `spatialcore.clustering` | Domains, neighborhoods, niches, community detection |
| `spatialcore.nmf` | Spatial non-negative matrix factorization |
| `spatialcore.diffusion` | Diffusion maps, pseudotime analysis |
| `spatialcore.ontology` | Cell ontology conversion and mapping |
| `spatialcore.annotation` | CellTypist wrappers, custom model training, benchmarking |
| `spatialcore.r_bridge` | Seurat integration via rpy2 |

## Supported Platforms

- 10x Genomics Xenium
- NanoString CosMx
- 10x Genomics Visium
- Any AnnData-compatible spatial data

## Terminology

We use consistent terminology throughout the package:

| Term | Definition |
|------|------------|
| **Neighborhood** | Immediate spatial context of a cell and similar groups of cells across the tissue|
| **Niche** | Functional microenvironment | - NEEDS UPDATE

## Ecosystem

SpatialCore integrates seamlessly with:

- [scanpy](https://scanpy.readthedocs.io/) - Single-cell analysis
- [squidpy](https://squidpy.readthedocs.io/) - Spatial omics analysis
- [spatialdata](https://spatialdata.scverse.org/) - Spatial data structures
- [Seurat](https://satijalab.org/seurat/) - R integration

## Package Philosophy

This package is **for computational biologists, by computational biologists**.

We prioritize:
- **Reproducibility** - Same inputs = same outputs across Python and R
- **Scalability** - Designed for large spatial datasets (millions of cells)
- **Simplicity** - Thin wrappers, not black boxes
- **Documentation** - Clear docstrings with references to original methods

We don't aim to:
- Invent new statistical methods
- Replace existing tools
- Abstract away the underlying math

## Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

## Citation

If you use SpatialCore in your research, please cite:

```bibtex
@software{spatialcore,
  title = {SpatialCore: Standardized spatial statistics for computational biology},
  url = {https://github.com/spatialcore/spatialcore},
  license = {Apache-2.0}
}
```

## License

Apache License 2.0. See [LICENSE](LICENSE) for details.

The SpatialCore name and trademarks are reserved to ensure the community can rely on the "Standardized" quality of the core library.
