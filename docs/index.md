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
*   **Scope**: We do **not** invent new methods, alogithims, or complex solutions; we make existing math work reliably.

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

### Recommended: Conda/Mamba Environment

For the best experience, we recommend using a conda or mamba environment:

```bash
# Create environment with Python 3.11
mamba create -n spatialcore python=3.11
mamba activate spatialcore

# Install SpatialCore
pip install spatialcore
```

This installs all core Python dependencies including CellTypist for custom model training for cell type annotation.


Run upgrade to get the latest modules, features, and fixes
```bash
# Activate environment
mamba activate spatialcore

# Upgrade
pip install --upgrade spatialcore
```

### R Requirements

SpatialCore uses R for certain operations that are statistically optimized or perform better in R. The `r_bridge` module handles R integration via subprocess (no rpy2 required).

**Install R packages in your environment:**

```bash
# If using conda/mamba (recommended)
mamba install -c conda-forge r-base r-sf r-concaveman r-dplyr r-purrr r-jsonlite

# If using system R (Linux/macOS)
sudo apt-get install r-base  # Ubuntu/Debian
R -e "install.packages(c('sf', 'concaveman', 'dplyr', 'purrr', 'jsonlite'), repos='https://cloud.r-project.org/')"
```

**Verify R is configured correctly:**

```python
from spatialcore.r_bridge import check_r_available, get_r_version
print(check_r_available())  # True
print(get_r_version())      # R version 4.x.x
```

### How r_bridge Works

The `r_bridge` automatically detects your environment:

| Environment | R Execution Method |
|-------------|-------------------|
| Conda/Mamba | `mamba run -n env_name Rscript ...` |
| System R | `Rscript` directly |

No manual configuration needed - it just works.

---


## 🚀 Quick Start

```python
import scanpy as sc
import spatialcore

# Check what's available in your installation
spatialcore.print_info()
# SpatialCore v0.4.5
# Available modules: core, annotation, nmf, r_bridge, spatial, ...

# Load your spatial data
adata = sc.read_h5ad("spatial_data.h5ad")

# Spatial domain detection on B cells
from spatialcore.spatial import make_spatial_domains
from spatialcore.plotting import plot_domains

adata = make_spatial_domains(
    adata,
    filter_expression="hieratype_ontology_name == 'B cell'",
    output_column="bcell_domain",
    domain_prefix="Bcell",
    platform="cosmx",
)

# Visualize spatial domains
plot_domains(adata, domain_col="bcell_domain", title="B Cell Domains")
```

See the [Domain Detection documentation](domains/domain_detection.md) for detailed tutorials and API reference.

---

## 🧩 Modules & Features

| Module | Status | Features |
|--------|--------|----------|
| **[`spatialcore.annotation`](celltyping/index.md)** | ✅ Available | CellTypist wrappers, custom model training, benchmarking |
| **[`spatialcore.nmf`](nmf/spanmf.md)** | ✅ Available | Spatial non-negative matrix factorization (spaNMF) |
| **[`spatialcore.spatial`](spatial/index.md)** | ✅ Available | Moran's I, Lee's L, spatial autocorrelation |
| **[`spatialcore.domains`](domains/index.md)** | ✅ Available | Neighborhood profiling, niche identification, domain detection |
| **[`spatialcore.thresholding`](thresholding/index.md)** | ✅ Available | Cell classification, oncogene thresholding |
| **[`spatialcore.r_bridge`](r_bridge/index.md)** | ✅ Available | Seurat integration, R interoperability via subprocess |
| **`spatialcore.diffusion`** | 🔜 Coming soon | Diffusion maps, pseudotime analysis |

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
