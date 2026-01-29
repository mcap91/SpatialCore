# SpatialCore

**Standardized spatial statistics for computational biology.**

[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![PyPI](https://img.shields.io/pypi/v/spatialcore.svg)](https://pypi.org/project/spatialcore/)
[![Code Style](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

---

## 🎯 The Mission

### The Problem
Tools for spatial biology analysis are fragmented. Implementations packages and complex functions often differ between languages (R vs Python) and even between packages, making reproducibility difficult and benchmarking impossible. We believe simple statistical tools solve many of the core problems of spatial biology, and we have desinged them to be intuative, easy to use, and scaleable for millions of cells. 

### The Solution
SpatialCore serves as a package for computational biologists, by computational biologists. It provides robust, standardized implementations of core spatial statistics that ensure identical results across platforms, wrapping high-performance libraries where available.

### The Goal
To make spatial analysis engineering boring, so you can focus on the exciting biology. **Standardized. Scalable. Reproducible.**

📖 **See the full documentation:** for modules, examples, and benchmarks [mcap91.github.io/SpatialCore](https://mcap91.github.io/SpatialCore)

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
import spatialcore

# Check what's available in your installation
spatialcore.print_info()
# SpatialCore v0.1.3
# Available modules: core, annotation

```


---

## 🧩 Modules & Features

| Module | Status | Features |
|--------|--------|----------|
| **`spatialcore.core`** | ✅ Available | Logging, metadata tracking, caching utilities |
| **`spatialcore.annotation`** | ✅ Available | CellTypist wrappers, custom model training, benchmarking |
| **`spatialcore.spatial`** | 🔜 Coming soon | Moran's I, Lee's L, neighborhoods, niches, domains |
| **`spatialcore.nmf`** | ✅ Available | Spatial non-negative matrix factorization (spaNMF) |
| **`spatialcore.diffusion`** | 🔜 Coming soon | Diffusion maps, pseudotime analysis |

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
*   Replacing Scanpy, Seurat, or other methods.

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