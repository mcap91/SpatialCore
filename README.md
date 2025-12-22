# SpatialCore

**Standardized spatial statistics for computational biology.**

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![PyPI](https://img.shields.io/pypi/v/spatialcore.svg)](https://pypi.org/project/spatialcore/)
[![Code Style](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

---

## 🎯 The Mission

### The Problem
Spatial biology analysis is fragmented. Implementations of basic statistics often differ between languages (R vs Python) and even between packages, making reproducibility difficult and benchmarking impossible. We spend too much time wondering if a result is biological or an artifact of the implementation.

### The Solution
SpatialCore serves as a **"ground truth" engineering layer**. It provides robust, standardized implementations of core spatial statistics that ensure identical results across platforms, wrapping high-performance libraries where available.

### The Goal
To make spatial analysis engineering boring, so you can focus on the exciting biology. **Standardized. Scalable. Reproducible.**

---

## 📦 Installation

```bash
pip install spatialcore
```

**Optional Dependencies:**

```bash
# For CellTypist automated annotation
pip install spatialcore[celltypist]

# For R/Seurat integration
pip install spatialcore[r]

# Install everything
pip install spatialcore[all]
```

---

## 🚀 Quick Start

```python
import scanpy as sc
import spatialcore as sp

# 1. Load your spatial data (Xenium, CosMx, Visium, etc.)
adata = sc.read_h5ad("spatial_data.h5ad")

# 2. Spatial Autocorrelation
# Calculate Moran's I for a specific gene
sp.spatial.morans_i(adata, gene="CD8A", spatial_key="spatial")

# Calculate Lee's L bivariate autocorrelation (e.g., co-localization)
sp.spatial.lees_l(adata, gene_x="CD8A", gene_y="CD4", spatial_key="spatial")

# 3. Define Neighborhoods & Domains
# Compute spatial domains using Leiden clustering on spatial graph
sp.clustering.compute_domains(adata, method="leiden", resolution=0.5)

# Calculate physical distances between different tissue domains
sp.clustering.domain_distance(adata, domain_key="domain")

# 4. Spatial Factorization
# Decompose expression into spatially coherent programs
sp.nmf.spatial_nmf(adata, n_components=10)
```

---

## 🧩 Modules & Features

| Module | Features |
|--------|-------------|
| **`spatialcore.spatial`** | • Global & Local Moran's I<br>• Bivariate Lee's L<br>• HH/LL/HL/LH Local Classification |
| **`spatialcore.clustering`** | • Spatial Domain Identification<br>• Niche Analysis<br>• Neighborhood definition |
| **`spatialcore.nmf`** | • Spatially-aware Non-negative Matrix Factorization<br>• Pattern extraction |
| **`spatialcore.diffusion`** | • Diffusion Maps<br>• Spatial Pseudotime Analysis |
| **`spatialcore.ontology`** | • Cell Ontology (CL) Mapping<br>• Standardization of labels |
| **`spatialcore.annotation`** | • Automated CellTypist Wrappers<br>• Custom Model Training<br>• Benchmarking |
| **`spatialcore.r_bridge`** | • **Seurat Integration** (via rpy2)<br>• Cross-language object conversion |

---

## 📚 Terminology

We strictly define our spatial units to ensure clarity:

| Term | Definition |
|------|------------|
| **Neighborhood** | The immediate spatial vicinity of a cell (e.g., k-Nearest Neighbors or fixed radius). |
| **Niche** | A functional microenvironment defined by a specific composition of cell types (e.g., "Tumor-immune border"). |
| **Domain** | A macroscopic, continuous tissue region with shared structural characteristics (e.g., "Cortex", "Medulla"). |

---

## 🤝 Ecosystem Integration

SpatialCore is designed to play nice with others. It fits seamlessly into the existing Python spatial biology stack:

*   **[Scanpy](https://scanpy.readthedocs.io/)**: The backbone for single-cell analysis.
*   **[Squidpy](https://squidpy.readthedocs.io/)**: Advanced spatial omics analysis.
*   **[SpatialData](https://spatialdata.scverse.org/)**: The OME-NGFF standard for spatial data storage.
*   **[Seurat](https://satijalab.org/seurat/)**: Direct R interoperability for teams working across languages.

---

## ⚖️ Philosophy

This package is **for computational biologists, by computational biologists**.

*   **Reproducibility**: Same inputs = Same outputs. Period.
*   **Scalability**: Built for the era of millions of cells (Xenium/CosMx).
*   **Transparency**: Thin wrappers, not black boxes. We verify, we don't obfuscate.
*   **Documentation**: Clear docstrings with academic references.

**What we are NOT:**
*   Inventing new, unproven math.
*   Replacing Scanpy or Seurat.

---

## 📝 Citation

If SpatialCore aids your research, please cite:

```bibtex
@software{spatialcore,
  title = {SpatialCore: Standardized spatial statistics for computational biology},
  url = {https://github.com/spatialcore/spatialcore},
  license = {Apache-2.0}
}
```

## License

**Apache License 2.0**

The SpatialCore name and trademarks are reserved to ensure the community can rely on the "Standardized" quality of the core library. You are free to use, modify, and distribute the code, including for commercial use.