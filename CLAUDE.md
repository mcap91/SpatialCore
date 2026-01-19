# SpatialCore - Claude Instructions

## Project Overview

SpatialCore is a Python/R package for spatial biology analysis. It provides standardized, robust wrappers around existing statistical tools for computational biologists working with spatial transcriptomics data (Xenium, CosMx, Visium).

**Philosophy:** Not novel math—just better, standardized engineering for the community.

## Repository Structure

SpatialCore uses a **two-repository model** for public/private separation:

| Repository | Visibility | Purpose |
|------------|------------|---------|
| `SpatialCore-Dev` | Private | Development repo with all code (main + dev branches) |
| `SpatialCore` | Public | Public-facing repo, auto-synced from main branch |

### GitHub Actions Sync

A GitHub Action (`.github/workflows/sync-public.yml`) automatically pushes the `main` branch to the public `SpatialCore` repo on every push. The sync uses an SSH deploy key.

**Secrets required in SpatialCore-Dev:**
- `DEPLOY_SSH_KEY` - Private key for pushing to public repo

**Deploy key in SpatialCore (public):**
- Public key with write access

## Branch Strategy

- **`dev` branch:** All development happens here - public code, preprocessing, tests, everything
- **`main` branch:** Public-ready code only - selectively copied from dev

### Environment
- **`spatialcore`** - single mamba environment for all development
- **Mamba path:** `C:\SpatialCore\miniforge3\envs\spatialcore\`
- **Python executable:** `C:\SpatialCore\miniforge3\envs\spatialcore\python.exe`

### Git Workflow
- **All development** happens on `dev` branch
- When ready to publish, use `scripts/publish.ps1` to selectively copy files to main
- **Never merge** dev → main (use selective file copying instead)
- Main branch only contains public-ready code

### Publishing Workflow

The publish script (`scripts/publish.ps1`) allows **selective** copying of modules and docs to main.

**Interactive mode** (recommended):
```powershell
./scripts/publish.ps1 -Interactive
```

**Direct mode** (specify what to include):
```powershell
./scripts/publish.ps1 -SrcModules "annotation","core" -DocsSections "celltyping"
```

**Dry run** (preview without changes):
```powershell
./scripts/publish.ps1 -DryRun -SrcModules "annotation" -DocsSections "celltyping"
```

**What gets copied:**
- Selected `src/spatialcore/` modules (e.g., annotation, spatial, core)
- Selected `docs/` sections (e.g., celltyping, spatial)
- `.github/workflows/` (sync-public.yml, claude.yml, docs.yml)
- Core files: `pyproject.toml`, `mkdocs.yml`, `README.md`, `LICENSE`, `CLAUDE.md`

**What the script does:**
1. Stashes uncommitted changes on dev
2. Switches to main branch
3. Copies only selected files/folders from dev
4. Commits and pushes to origin/main
5. Returns to dev and restores stash
6. GitHub Action automatically syncs main → public repo

## Cross-Platform Development (Windows → Linux)

SpatialCore is developed on **Windows** and deployed on **Linux** (production/cluster). This creates specific challenges that the codebase addresses.

### Platform Differences

| Aspect | Windows Development | Linux Production |
|--------|--------------------|--------------------|
| Shell | PowerShell | Bash |
| Path separator | Backslash `\` | Forward slash `/` |
| Environment activation | `mamba activate env` | `mamba activate env` |
| R execution | `mamba run` (automatic via r_bridge) | `mamba run` or direct `Rscript` |
| Python executable | `.exe` suffix | No suffix |

### Path Handling Rules (CRITICAL)

**All paths passed to R must use forward slashes**, even on Windows. R's `source()` function fails with backslashes.

```python
# CORRECT - Use pathlib and convert to posix for R
from pathlib import Path
r_script_path = Path(__file__).parent / "r_functions.R"
r_code = f"source('{r_script_path.as_posix()}')"

# CORRECT - Replace backslashes in string paths
csv_path = str(input_csv).replace("\\", "/")

# WRONG - Backslashes break R's source()
r_code = f"source('{r_script_path}')"  # Fails on Windows!
```

See `src/spatialcore/spatial/domains.py` for the canonical pattern.

### Environment Variables

| Variable | Purpose | Set By |
|----------|---------|--------|
| `RETICULATE_PYTHON` | Tells R's reticulate which Python to use | r_bridge automatically |
| `CONDA_DEFAULT_ENV` | Current mamba environment name | Shell activation |
| `CONDA_EXE` | Path to mamba executable | Mamba installation |

### PowerShell Configuration (Windows)

**One-time setup** - Initialize mamba for PowerShell:
```powershell
C:\SpatialCore\miniforge3\condabin\mamba.bat init powershell
# Restart PowerShell after this
```

**If scripts are blocked**, set execution policy:
```powershell
Set-ExecutionPolicy -ExecutionPolicy RemoteSigned -Scope CurrentUser
```

**Activate environment and verify**:
```powershell
mamba activate spatialcore
$env:CONDA_DEFAULT_ENV          # Should show "spatialcore"
python --version                # Should show Python from mamba env
Get-Command python              # Should point to mamba env
```

**Running Python scripts**:
```powershell
mamba activate spatialcore
python C:\SpatialCore\SpatialCore\scripts\vignettes\generate_domain_detection.py
```

### WSL Warning

Do **NOT** mix WSL mamba with Windows mamba - they are completely separate installations:
- Windows paths (`C:\...`) won't work in WSL
- WSL uses `/mnt/c/...` for Windows drives
- Keep Windows and WSL environments separate

## Environment

- Use the `spatialcore` mamba environment for all development
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
5. **Check upstream documentation** before implementing - use web search to verify current API signatures and best practices:
   - scanpy: https://scanpy.readthedocs.io/
   - squidpy: https://squidpy.readthedocs.io/
   - Seurat: https://satijalab.org/seurat/

## Parallelization Strategy

- **Default:** Sequential sample processing (memory-safe for large spatial datasets)
- **Within-function:** Use numba/parallel operations where available (scanpy defaults)
- **Optional:** Parallel samples via multiprocessing if memory permits

## R Integration

SpatialCore uses **subprocess-based R execution** (NOT rpy2) for maximum portability and to avoid DLL/library conflicts. The `r_bridge` module provides:

- `run_r_script()` - Execute R scripts with proper environment configuration
- `run_r_code()` - Execute inline R code snippets
- `check_r_available()` - Verify R installation
- `get_r_version()` - Get installed R version

### How It Works

1. Python calls R via subprocess using `mamba run -n env_name Rscript` (Windows) or direct `Rscript` (Linux)
2. `RETICULATE_PYTHON` is set automatically to `sys.executable`
3. Data exchange uses CSV files for tabular data, JSON for results
4. R scripts output JSON on their last stdout line for parsing

### Required R Packages

For spatial domain creation and related operations:
```
sf              # Simple Features for geometry operations
concaveman      # Concave hull algorithm
dplyr           # Data manipulation
purrr           # Functional programming
jsonlite        # JSON output for Python parsing
```

### Example Usage

```python
from spatialcore.r_bridge import check_r_available, run_r_script

if check_r_available():
    result = run_r_script(
        "my_analysis.R",
        args=["input.csv", "output.csv"],
        timeout=600,
    )
    print(result)  # Parsed JSON from R stdout
```

### Why Subprocess Instead of rpy2?

1. **DLL isolation** - rpy2 links R DLLs into Python process, causing conflicts
2. **Environment flexibility** - Works with mamba, system R, or mixed setups
3. **Debugging** - Easier to debug R scripts in isolation
4. **Cross-platform** - Same approach works on Windows, Linux, macOS

### Windows R Configuration (Development)

On Windows, R must be installed within the mamba environment to ensure DLL paths are correctly configured.

**Install R and packages**:
```powershell
mamba activate spatialcore
mamba install -c conda-forge r-base r-sf r-concaveman r-dplyr r-purrr r-jsonlite
```

**Why `mamba run` is required on Windows**: When Python launches a subprocess, the mamba environment's DLLs are not in PATH. R fails with errors like `DLL 'proj.dll' not found`. The r_bridge automatically uses `mamba run -n env_name Rscript` on Windows (see `subprocess_runner.py` lines 339-359).

**R location in mamba environment**:
```
C:\SpatialCore\miniforge3\envs\spatialcore\
├── Lib\R\
│   ├── bin\
│   │   ├── Rscript.exe
│   │   └── R.exe
│   └── library\
│       ├── sf\
│       ├── dplyr\
│       └── ...
└── Library\mingw-w64\bin\  # DLLs: proj.dll, gdal.dll, etc.
```

**Verify R installation**:
```powershell
mamba run -n spatialcore Rscript -e "library(sf); print('sf loaded')"

# Or from Python
python -c "from spatialcore.r_bridge import check_r_available, get_r_version; print(check_r_available(), get_r_version())"
```

### Linux R Configuration (pip Users)

Users who `pip install spatialcore` on Linux need to configure R separately.

**Option 1: System R (Recommended for production)**:
```bash
# Ubuntu/Debian
sudo apt-get install r-base

# Install R packages
R -e "install.packages(c('sf', 'concaveman', 'dplyr', 'purrr', 'jsonlite'), repos='https://cloud.r-project.org/')"
```

**Option 2: Mamba R (Recommended for isolation)**:
```bash
mamba create -n spatialcore python=3.11
mamba activate spatialcore
pip install spatialcore
mamba install -c conda-forge r-base r-sf r-concaveman r-dplyr r-purrr r-jsonlite
```

**Verify R works with spatialcore**:
```bash
python -c "from spatialcore.r_bridge import check_r_available; print(check_r_available())"
```

### R Troubleshooting

**"Rscript not found in conda environment or PATH"**
- Windows: `mamba install -c conda-forge r-base`
- Linux: `apt-get install r-base` or `mamba install r-base`

**"DLL not found" / "sf package load failed" (Windows)**
- Verify environment activation: `$env:CONDA_DEFAULT_ENV` should show environment name
- Verify mamba is found: `python -c "from spatialcore.r_bridge.subprocess_runner import _find_mamba_or_conda; print(_find_mamba_or_conda())"`
- Test mamba run directly: `mamba run -n spatialcore Rscript -e "library(sf)"`

**"Error in source('path\\to\\script.R')"**
- Backslash path issue - use `Path.as_posix()` or `str(path).replace("\\", "/")`

**"reticulate could not find Python"**
- r_bridge sets `RETICULATE_PYTHON` automatically
- For manual R sessions: `Sys.setenv(RETICULATE_PYTHON = "path/to/python")`

**R script timeout**
- Default: 60s for `run_r_code()`, 3600s for `run_r_script()`
- Increase for large datasets: `run_r_script("script.R", timeout=7200)`

**Debugging R subprocess calls**:
```python
import logging
logging.getLogger("spatialcore.r_bridge").setLevel(logging.DEBUG)

from spatialcore.r_bridge import run_r_code
result = run_r_code('cat(jsonlite::toJSON(list(x=1)))')
# Logs will show: Command, Working directory, RETICULATE_PYTHON
```

**PowerShell-specific issues**:
- "mamba: command not found" → run `mamba init powershell` and restart shell
- "Execution policy" errors → `Set-ExecutionPolicy RemoteSigned -Scope CurrentUser`
- Environment not activated → check `$env:CONDA_DEFAULT_ENV` is set
- Wrong Python → check `Get-Command python` points to mamba env

## Project Structure

### Main Branch (synced to public SpatialCore repo)

```
SpatialCore/
├── .github/workflows/         # GitHub Actions (sync-public.yml, claude.yml, docs.yml)
├── src/
│   └── spatialcore/
│       ├── __init__.py
│       ├── core/              # Logging, metadata, caching utilities
│       ├── spatial/           # Spatial stats, neighborhoods, niches, domains
│       ├── clustering/        # Graph-based clustering (Leiden, Louvain)
│       ├── nmf/               # Spatial NMF
│       ├── diffusion/         # Diffusion maps, pseudotime
│       ├── ontology/          # Ontology conversion tools
│       ├── annotation/        # CellTypist wrappers, benchmarking
│       └── r_bridge/          # R/Python integration
├── docs/                      # Documentation (mkdocs)
├── .gitignore
├── LICENSE
├── README.md
├── pyproject.toml
└── mkdocs.yml
```

### Dev Branch (private, all development)

```
SpatialCore/
├── (everything from main, plus:)
├── preprocessing/             # Internal preprocessing pipelines
│   ├── ingestion/             # Data loading from CosMx, Xenium, Visium
│   ├── qc/                    # Quality control metrics and filtering
│   ├── normalization/         # Expression normalization pipelines
│   ├── celltyping/            # Cell type assignment workflows
│   └── pipelines/             # Config-driven preprocessing pipelines
├── tests/                     # Test suite
├── scripts/                   # Vignette generation, utilities
├── R/                         # R package source
└── .cache/                    # Intermediate results (gitignored)
```

## Module Responsibilities (main branch)

| Module | Purpose |
|--------|---------|
| `core` | Logging, metadata tracking, caching utilities |
| `spatial` | Moran's I, Lee's L, neighborhoods, niches, domains, distances |
| `clustering` | Graph-based clustering (Leiden, Louvain on spatial graphs) |
| `nmf` | Spatial non-negative matrix factorization |
| `diffusion` | Diffusion maps, pseudotime analysis |
| `ontology` | Cell ontology conversion and mapping |
| `annotation` | CellTypist wrappers, custom model training, benchmarking |
| `r_bridge` | Subprocess-based R execution, RETICULATE_PYTHON config, JSON parsing |

## Preprocessing Modules (dev branch only)

The `preprocessing/` directory contains internal data processing code that is never published:

| Module | Purpose |
|--------|---------|
| `ingestion` | Data loading from CosMx, Xenium, Visium platforms |
| `qc` | Quality control metrics and filtering |
| `normalization` | Expression normalization pipelines |
| `celltyping` | Cell type assignment workflows |
| `pipelines` | Config-driven preprocessing runner |

### Using Preprocessing Modules

```python
# In spatialcore environment:
pip install -e .

# Import preprocessing functions:
from preprocessing.ingestion import load_xenium, load_cosmx
from preprocessing.pipelines import run_pipeline

# Use core utilities:
from spatialcore.core.logging import get_logger
```

## Documentation & Vignettes

### Strategy: Static Markdown (No Code Execution)

All public-facing vignettes are **static markdown files** with pre-generated images. Reports do NOT execute code at render time.

**Why static?**
- Reproducibility: Results are fixed at generation time
- Hosting compatibility: Works with ReadTheDocs, GitHub Pages, Sphinx
- No environment dependencies: Users can view without Python setup

### Vignette Structure

Vignettes live directly in `docs/[topic]/` alongside API documentation. Each vignette has an adjacent `_images/` folder for plots.

```
docs/
├── spatial/
│   ├── morans_i.md                # Vignette markdown
│   └── morans_i_images/           # Pre-generated plots
├── thresholding/
│   ├── multivariate_cutoff.md
│   └── multivariate_cutoff_images/
└── ...
```

### Vignette Generation Scripts

Scripts that generate vignette plots live in `scripts/vignettes/`. They save plots directly to `docs/[topic]/[name]_images/`.

```
scripts/
└── vignettes/
    ├── _template.py               # Copy this for new vignettes
    ├── generate_morans_i.py
    └── generate_lees_l.py
```

### Vignette Workflow

1. **Create/update script** in `scripts/vignettes/` that loads data and generates plots
2. **Script saves plots directly** to `docs/[topic]/[name]_images/`
3. **Write/update markdown** in `docs/[topic]/[name].md`
4. **Commit both** markdown and images together

### Data References

Data stays outside `docs/`. Each vignette should include a data block:

```markdown
!!! info "Data"
    This analysis uses the CosMx Colon dataset.
    **Public access:** [CELLxGENE Census](https://cellxgene.cziscience.com/...)
```

### Legacy Quarto (Internal Only)

The `reports/` directory contains Quarto templates (`.qmd`) for internal preprocessing reports on the dev branch. These execute code and are not published.

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
- numba (optional, for parallelization)

### R (for spatial domain creation)
- sf (Simple Features for geometry)
- concaveman (concave hull algorithm)
- dplyr, purrr (data manipulation)
- jsonlite (JSON output for Python parsing)

### R (for testing - optional)
- testthat

## Terminology Standards

Use these definitions consistently (see `docs_refs/vocab.md` for full rationale with literature citations):

| Term | Definition | Key Property |
|------|------------|--------------|
| **Neighborhood** | Immediate spatial vicinity of a cell (k-NN or radius). The *input* for niche/domain analysis. | Local spatial context |
| **Niche** | Cellular microenvironment archetype defined by cell-type composition. Can appear in multiple non-contiguous locations. | Compositional similarity (location-independent) |
| **Domain** | Spatially contiguous tissue region. Partitions tissue into discrete, bounded areas. | Spatial contiguity |

**Conceptual workflow:** `Neighborhood → Niche → Domain` (local → what kind → where)

**Deprecated:** "Community" — use "Niche" (compositional) or "Domain" (spatial) instead.
