# <img src="assets/logo.png" width="150" align="absmiddle" style="margin-right: 30px;"> SpatialCore

**Standardized spatial statistics for computational biology.**

---

## 🎯 The Mission

Spatial biology analysis is fragmented. Implementations of basic statistics often differ between languages (R vs Python) and even between packages, making reproducibility difficult. We spend too much time wondering which method or package is best, and how to interpet them.

**SpatialCore serves as a foundational engineering layer.**

We provide robust, intuiative, and standardized implementations of core spatial statistics that ensure identical results across platforms. Our goal is to make spatial analysis engineering boring, so you can focus on the exciting biology.

**GitHub:** [mcap91/SpatialCore](https://github.com/mcap91/SpatialCore)

**Core Principles**{: .section-label }

This package is **for computational biologists, by computational biologists**.

*   **Reproducibility**: Same inputs = Same outputs. Period.
*   **Scalability**: Built for the era of millions of cells (Xenium/CosMx).
*   **Transparency**: Thin wrappers, not black boxes. We verify, we don't obfuscate.
*   **Scope**: We do **not** invent new math; we make existing math work reliably.

**Ecosystem Integration**{: .section-label }

SpatialCore is designed to play nice with others. It fits seamlessly into the existing Python spatial biology stack:

*   **[Scanpy](https://scanpy.readthedocs.io/)**: The backbone for single-cell analysis.
*   **[Squidpy](https://squidpy.readthedocs.io/)**: Advanced spatial omics analysis.
*   **[Seurat](https://satijalab.org/seurat/)**: Direct R interoperability for teams working across languages.

Given the acceleration and widespread adoption of agentic coding assistants, this package is primarily written in python and is designed for both manual implementation, as well as discovery by agents. There are some paralell R function implementations, for those who wish to use them. Some functinality is eaither exclusive or best suited for R, and these functions are called via the r_bridge functionality. More information can be found in the doc strings.

---

## 📚 Terminology

We strictly define our spatial units to ensure clarity and alignment with established literature:

*   **Neighborhood**: The immediate spatial vicinity of a cell, typically defined by k-Nearest Neighbors (k-NN) or a fixed radius. Neighborhoods are the **input** for downstream niche and domain analysis.
    > Citations: [COVET, Nat. Biotech. (2024)](https://www.nature.com/articles/s41587-024-02193-4) | [NicheFlow (2025)](https://www.alphaxiv.org/overview/2511.00977v1)

*   **Niche**: A cellular microenvironment archetype defined by cell-type composition. Niches represent clusters of neighborhoods with similar profiles, **independent of their physical location**.
    > Citations: [Hu et al., Nat. Gen. (2025)](https://www.nature.com/articles/s41588-025-02120-6) | [ONTraC, Genome Biology (2025)](https://link.springer.com/article/10.1186/s13059-025-03588-5) | [BANKSY, Nat. Gen. (2024)](https://www.nature.com/articles/s41588-024-01664-3)

*   **Domain**: A spatially contiguous tissue region with coherent expression patterns and/or cell-type composition. Domains partition the tissue into discrete, bounded areas.
    > Citations: [BANKSY, Nat. Gen. (2024)](https://www.nature.com/articles/s41588-024-01664-3) | [Benchmark study, NAR (2024)](https://academic.oup.com/nar/article/53/7/gkaf303/8114322)

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

| Module | Status | Features |
|--------|--------|----------|
| **`spatialcore.core`** | ✅ Available | Logging, metadata tracking, caching utilities |
| **[`spatialcore.annotation`](celltyping/index.md)** | ✅ Available | CellTypist wrappers, custom model training, benchmarking |
| **`spatialcore.spatial`** | 🔜 Coming soon | Moran's I, Lee's L, neighborhoods, niches, domains |
| **`spatialcore.nmf`** | 🔜 Coming soon | Spatial non-negative matrix factorization |
| **`spatialcore.diffusion`** | 🔜 Coming soon | Diffusion maps, pseudotime analysis |
| **`spatialcore.r_bridge`** | 🔜 Coming soon | Seurat integration via subprocess |

---

## Development

*Details for developers and contributors will be added here.*

---

## 📝 Citation

If SpatialCore aids your research, please cite:

```bibtex
@software{spatialcore,
  title = {SpatialCore: Standardized spatial statistics for computational biology},
  url = {https://github.com/mcap91/SpatialCore},
  license = {Apache-2.0}
}
```

## License

**Apache License 2.0**

The SpatialCore name and trademarks are reserved to ensure the community can rely on the "Standardized" quality of the core library. You are free to use, modify, and distribute the code, including for commercial use.
