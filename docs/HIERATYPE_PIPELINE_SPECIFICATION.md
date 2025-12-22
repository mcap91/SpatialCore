# HieraType Hierarchical Cell Typing Pipeline Specification

## Purpose

This document provides a complete specification for rebuilding the HieraType hierarchical cell typing pipeline as standalone scripts without Temporal/process_activity dependencies. Use this to create a stepwise data processing pipeline.

---

## Implementation Notes (Updated 2025-12-07)

### Key Findings from Production Implementation

This section documents critical implementation details discovered during development that differ from or extend the original specification.

#### 1. R Bridge Configuration with Reticulate

When calling R from Python via subprocess, reticulate must be configured to use the correct Python environment:

```python
# In r_bridge.py - CRITICAL: Pass RETICULATE_PYTHON to subprocess
import os
import sys

env = os.environ.copy()
env["RETICULATE_PYTHON"] = sys.executable  # Use current Python interpreter

result = subprocess.run(
    ["Rscript", str(wrapper_script)],
    capture_output=True,
    text=True,
    timeout=timeout,
    cwd=r_functions_path.parent,
    env=env,  # Pass environment with RETICULATE_PYTHON
)
```

**Why this matters:** Without `RETICULATE_PYTHON`, reticulate will download its own Python installation to `~/.cache/R/reticulate/`, which won't have anndata, scipy, or other required packages.

#### 2. Sparse Matrix Handling with Reticulate

Reticulate **automatically converts** scipy sparse matrices to R Matrix classes:
- `scipy.sparse.csr_matrix` → `dgRMatrix`
- `scipy.sparse.csc_matrix` → `dgCMatrix`

**Do NOT manually extract indices/indptr/data.** Instead:

```r
# CORRECT - reticulate auto-converts
X <- adata$X
if (inherits(X, "sparseMatrix")) {
    counts_matrix <- as(X, "CsparseMatrix")
    rownames(counts_matrix) <- cell_names
    colnames(counts_matrix) <- gene_names
}
```

```r
# WRONG - will fail with "$ operator not defined for this S4 class"
X_csc <- scipy_sparse$csc_matrix(X)
indices <- X_csc$indices  # ERROR!
```

#### 3. Accessing AnnData obsp Keys

The `adata$obsp` object is a Python dictionary-like object. To get keys:

```r
# CORRECT - use iterate() to convert Python iterator
obsp_keys <- unlist(iterate(adata$obsp$keys()))

# WRONG - returns method names, not keys
obsp_keys <- names(py_to_r(adata$obsp))  # Returns: "as_dict", "keys", "values", etc.
```

#### 4. Matrix Orientation for HieraType

HieraType's `fit_metagene_scores` function expects genes in **columns** (cells × genes format):

```r
# CORRECT - pass cells × genes, genes in colnames
hieratype_result <- run_pipeline(
    pipeline = pipeline,
    counts_matrix = counts_matrix,  # cells × genes
    ...
)
```

The function checks `colnames(counts_matrix)` for gene names. If genes are not found in colnames, it checks rownames and auto-transposes if needed. However, if not ALL marker genes are present (which is common), the fallback logic may fail.

#### 5. HieraType Marker Validation

HieraType's internal `validate_markerslist` function:
- **Errors** if ALL index markers for a cell type are missing
- **Warns** (but continues) if only SOME index markers are missing

The function automatically filters missing markers from the list. Our R wrapper adds additional validation to provide clearer logging.

#### 6. Required R Packages

Beyond the core packages, these are required at runtime:

```bash
# Install via mamba in spatial-bio environment
mamba install -n spatial-bio -y r-mvtnorm r-matrixstats r-fnn r-reticulate
```

| Package | Purpose |
|---------|---------|
| `r-mvtnorm` | Multivariate normal distribution (used by HieraType clustering) |
| `r-matrixstats` | Fast matrix statistics |
| `r-fnn` | Fast nearest neighbors (for spatial kNN option) |
| `r-reticulate` | Python-R interoperability |

#### 7. Output Visualization

The test script generates 5 standard plots:

| File | Description |
|------|-------------|
| `spatial_celltype.png` | XY scatter plot colored by cell type |
| `celltype_pie.png` | Pie chart of cell type composition |
| `celltype_umi_counts.png` | Bar plot: mean UMI counts per cell type |
| `celltype_feature_counts.png` | Bar plot: mean genes detected per cell type |
| `celltype_summary_table.png` | Summary table with counts, %, scores, UMI, genes |

---

## Complete Data Flow

```
┌─────────────────────────────────────────────────────────────────┐
│ STEP 1: PREPARE INPUT DATA                                       │
│ Input:  Normalized h5ad with SNN graph OR spatial coordinates   │
│ Process: Validate log-normalization, ensure neighbors computed  │
│ Output: h5ad ready for HieraType                                │
└─────────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 2: LOAD & VALIDATE MARKER CONFIGURATION                    │
│ Input:  JSON marker hierarchy file                              │
│ Process: Parse levels → Validate markers against genes          │
│ Output: HieraType markerlist objects (L1, L2, LT)               │
└─────────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 3: PREPARE ADJACENCY MATRIX                                │
│ Input:  adata.obsp['connectivities'] OR adata.obsm['spatial']  │
│ Process: Load SNN graph OR compute spatial kNN                  │
│ Output: Sparse cells × cells adjacency matrix                   │
└─────────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 4: RUN HIERATYPE PIPELINE                                  │
│ Input:  Expression matrix + adjacency + marker pipeline         │
│ Process: NNLS scoring → neighborhood averaging → hierarchical   │
│ Output: Cell type assignments at L1, L2, LT levels              │
└─────────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 5: EXTRACT & SAVE RESULTS                                  │
│ Input:  HieraType result object                                 │
│ Process: Extract per-level assignments → Create final labels   │
│ Output: adata.obs with 8 new annotation columns                 │
└─────────────────────────────────────────────────────────────────┘
```

---

## Marker Configuration Format

### JSON Structure (3-Level Hierarchy)

```json
{
  "l1": {
    "cell_type_name": {
      "index_marker": ["GENE1", "GENE2"],
      "predictors": ["GENE1", "GENE2", "GENE3", ...],
      "use_offclass_markers_as_negative_predictors": true
    }
  },
  "immunemajor": { ... },
  "tcellmajor": { ... }
}
```

### Marker Properties

- **index_marker**: Essential markers (must have at least one present)
  - Can be string or array of strings
  - Used to filter out cell types with no detectable markers

- **predictors**: Comprehensive gene list for scoring (50-200 genes)
  - Includes markers, co-factors, downstream targets

- **use_offclass_markers_as_negative_predictors**: Use genes from OTHER cell types as negative signal (default: true)

### Complete Default Configuration

```json
{
  "l1": {
    "epithelial": {
      "index_marker": ["EPCAM", "KRT18", "KRT19", "CDH1"],
      "predictors": ["EPCAM", "KRT18", "KRT8", "KRT19", "KRT7", "KRT5", "KRT14", "CDH1", "CDH3", "CD24", "CLDN1", "CLDN3", "CLDN4", "CLDN7", "OCLN", "TJP1", "DSP", "DSG1", "DSG2", "DSG3", "MUC1", "MUC4", "MUC16", "CEACAM5", "CEACAM6", "SFN", "SPINT1", "ANPEP", "SLC2A1", "GJB2", "GJB6", "TFF1", "TFF3", "CD9", "ITGA6", "ITGB4", "CD151", "F11R", "JUP", "GATA3", "FOXA1", "TACSTD2", "FXYD3", "ST14", "AGR2", "AGR3", "CLCA2", "CXADR", "SLPI", "CSTA", "SERPINB5", "SERPINB3", "S100A14", "S100A16", "S100A2", "S100P"],
      "use_offclass_markers_as_negative_predictors": true
    },
    "endothelial": {
      "index_marker": ["PECAM1", "VWF", "CDH5"],
      "predictors": ["PECAM1", "VWF", "CDH5", "CLDN5", "ESAM", "ENG", "FLT1", "KDR", "TIE1", "TEK", "NOTCH4", "ECSCR", "ICAM2", "MCAM", "CAV1", "CD34", "SELE", "VCAM1", "ICAM1", "PLVAP", "PDGFB", "EFNB2", "ESM1", "ANGPT2", "ANGPT1", "CD36", "LYVE1", "PROX1", "ACKR1", "S100A4", "NR2F2", "ERG"],
      "use_offclass_markers_as_negative_predictors": true
    },
    "fibroblast": {
      "index_marker": ["COL1A1", "DCN", "PDGFRA"],
      "predictors": ["COL1A2", "VIM", "S100A4", "PDGFRB", "PDGFRA", "ACTA2", "COL1A1", "COL3A1", "COL5A1", "FN1", "THY1", "DCN", "LUM", "MMP2", "MMP9", "MMP1", "MMP3", "TIMP1", "TIMP2", "SPARC", "FAP", "DDR2", "POSTN", "SDC1", "SDC4", "PDPN", "CD44", "ITGB1", "MMP14"],
      "use_offclass_markers_as_negative_predictors": true
    },
    "plasma": {
      "index_marker": ["IGHG1", "IGHG1/2", "IGHA1", "MZB1", "JCHAIN"],
      "predictors": ["IGHG1", "IGHG1/2", "IGHA1", "MZB1", "PRDM1", "XBP1", "JCHAIN", "CD38", "IRF4", "PPIB", "SDC1", "TNFRSF17", "HSP90B1", "HYOU1", "SSR4", "DERL3", "FKBP11", "SEL1L", "CST3", "SEC61B", "SEC61A1", "GPRC5D", "HSPA5", "CXCR4"],
      "use_offclass_markers_as_negative_predictors": true
    },
    "immune": {
      "index_marker": ["PTPRC", "ITGB2", "SPI1", "CD52", "CORO1A", "HLA-DRA"],
      "predictors": ["PTPRC", "CD3E", "CD3D", "CD4", "CD8A", "CD8B", "CD19", "CD20", "CD14", "CD68", "ITGAM", "ITGAX", "CD56", "FCGR3A", "CD86", "CD80", "HLA-DRA", "CSF1R", "CD163", "MRC1", "CD33", "CEACAM8", "FCER1A", "IL3RA", "NCF1", "LYZ", "CCR7", "CXCR3", "CXCR4", "SELL", "ITGAL", "CD27", "CD28", "CD38", "CD40", "CD44", "CD69", "CD83", "GZMB", "PRF1", "NKG7", "KLRD1"],
      "use_offclass_markers_as_negative_predictors": true
    },
    "smooth_muscle": {
      "index_marker": ["MYH11", "ACTA2", "TAGLN", "DES"],
      "predictors": ["MYH11", "ACTA2", "TAGLN", "DES", "CNN1", "SYNPO2", "LMOD1", "MYLK", "MYL9", "TPM1", "TPM2", "TPM4", "CALD1", "PLN", "PDGFRB", "NOTCH3", "ENG", "COL1A1", "COL3A1", "COL5A1", "MMP2", "TIMP1", "FN1", "SPARC", "VIM"],
      "use_offclass_markers_as_negative_predictors": true
    }
  },
  "immunemajor": {
    "tcell": {
      "index_marker": ["CD3D", "CD3G", "CD3E", "IL7R", "TRBC1", "CD2"],
      "predictors": ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "IL7R", "CD27", "SELL", "CCR7", "TCF7", "LEF1", "FOXP1", "BACH2", "KLF2", "CD28", "ICOS", "STAT3", "STAT5", "RUNX3", "BATF", "EOMES", "TBX21", "PRDM1", "PDCD1", "LAG3", "TIGIT", "HAVCR2", "GZMK", "GZMA", "GZMB", "PRF1", "GNLY", "KLRG1", "KLRD1", "IFNG", "CCL5", "CX3CR1", "NKG7", "FAS", "FASLG", "IL2", "CD244", "CXCR3", "KLRB1", "TOX", "NR4A2", "BCL6", "CXCL13", "CD4", "CD40LG", "RUNX1", "GATA3", "RORC", "FOXP3", "IL2RA", "CTLA4"],
      "use_offclass_markers_as_negative_predictors": true
    },
    "bcell": {
      "index_marker": ["CD19", "MS4A1", "CD79A", "CD79B"],
      "predictors": ["CD19", "MS4A1", "CD79A", "CD79B", "PAX5", "BLK", "VPREB1", "BANK1", "SPIB", "CD22", "CR2", "CXCR5", "CD40", "TNFRSF13B", "TNFRSF13C", "FCMR", "CD72", "IL4R", "CD180", "CD24", "CD38", "IGHM", "IGHD", "AICDA", "BACH2", "IRF4", "IRF8", "BCL6", "PRDM1"],
      "use_offclass_markers_as_negative_predictors": true
    },
    "macrophage": {
      "index_marker": ["CD163", "CD68", "MARCO", "CSF1R"],
      "predictors": ["CD68", "CD163", "MARCO", "CSF1R", "MERTK", "SPP1", "C1QA", "C1QB", "C1QC", "SIRPA", "AIF1", "TGFB1", "NOS2", "ARG1", "IL10", "CD206", "MRC1", "IL4R", "CD80", "CD86", "HLA-DRA", "HLA-DRB1", "PTPRC", "FCGR1A", "ITGAM", "ITGAX"],
      "use_offclass_markers_as_negative_predictors": true
    },
    "neutrophil": {
      "index_marker": ["MPO", "FCGR3B", "S100A8", "S100A9"],
      "predictors": ["MPO", "FCGR3B", "S100A8", "S100A9", "CXCR2", "CSF3R", "IL1B", "NLRP3", "CD177", "ELANE", "PRTN3", "RETN", "LTF", "ARG1", "OLFM4", "TLR4", "CEBPB", "CXCL8", "G0S2", "ITGAM", "FPR1", "TLR2", "CEACAM8"],
      "use_offclass_markers_as_negative_predictors": true
    },
    "mast": {
      "index_marker": ["KIT", "TPSAB1", "TPSB2", "TPSAB1/2", "CPA3"],
      "predictors": ["KIT", "TPSAB1", "CPA3", "TPSB2", "TPSAB1/2", "HDC", "FCER1A", "MS4A2", "FUT3", "SIGLEC6", "SIGLEC8", "CD123", "CD9", "IL4", "IL13", "IL3RA", "PTGDR2", "CCL2", "CCL3", "CCL4", "XBP1", "GATA2"],
      "use_offclass_markers_as_negative_predictors": true
    },
    "monocyte": {
      "index_marker": ["CD14", "FCGR3A", "CCR2", "LYZ"],
      "predictors": ["CD14", "FCGR3A", "CCR2", "LYZ", "CSF1R", "HLA-DRA", "HLA-DRB1", "S100A8", "S100A9", "S100A12", "IL1B", "TNF", "CX3CR1", "CD36", "VCAN", "CD86", "FCGR1A", "TLR2", "TLR4", "FPR1", "NLRP3", "PTPRC", "IRF5", "IRF7", "C1QA", "C1QB", "C1QC"],
      "use_offclass_markers_as_negative_predictors": true
    },
    "dendritic": {
      "index_marker": ["ITGAX", "CLEC9A", "FLT3", "XCR1", "CD1C"],
      "predictors": ["ITGAX", "CLEC9A", "FLT3", "XCR1", "LAMP3", "IL3RA", "HLA-DRA", "HLA-DRB1", "IRF4", "IRF8", "CD80", "CD86", "CD274", "CCR7", "BATF3", "CD1C", "CD141", "THBD", "SIRPA", "NLRP3", "IL12B", "IL15", "CCL17", "CCL19", "CCL22"],
      "use_offclass_markers_as_negative_predictors": true
    },
    "natural_killer_cell": {
      "index_marker": ["CD69", "GNLY", "KLRD1", "KLRF1", "NKG7", "FCGR3A", "FCGR3A/B"],
      "predictors": ["CD69", "GNLY", "KLRD1", "KLRF1", "NKG7", "GZMA", "GZMB", "GZMH", "PRF1", "FCGR3A", "FCGR3A/B", "NCAM1", "KLRC1", "KLRC2", "IL2RB", "KLRK1", "NCR1", "NCR3", "XCL1", "XCL2", "XCL1/2", "TYROBP", "IFNG", "TNF"],
      "use_offclass_markers_as_negative_predictors": true
    },
    "eosinophil": {
      "index_marker": ["CLC", "PRG2", "SIGLEC8", "EPX", "JAML"],
      "predictors": ["CLC", "PRG2", "SIGLEC8", "EPX", "JAML", "ALOX5AP", "CCR3", "HRH4", "IL9", "IL5RA", "IL5", "RNASE2", "RNASE3", "TNFRSF12A", "ITGA4", "CEACAM8"],
      "use_offclass_markers_as_negative_predictors": true
    }
  },
  "tcellmajor": {
    "cd8t": {
      "index_marker": ["CD8B", "CD8A"],
      "predictors": ["CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "IL7R", "CD27", "SELL", "CCR7", "TCF7", "LEF1", "FOXP1", "BACH2", "KLF2", "CD28", "ICOS", "STAT3", "STAT5", "RUNX3", "BATF", "EOMES", "TBX21", "PRDM1", "PDCD1", "LAG3", "TIGIT", "HAVCR2", "GZMK", "GZMA", "GZMB", "PRF1", "GNLY", "KLRG1", "KLRD1", "IFNG", "CCL5", "CX3CR1", "NKG7", "FAS", "FASLG", "IL2", "CD244", "CXCR3", "KLRB1", "TOX", "NR4A2", "BCL6", "CXCL13"],
      "use_offclass_markers_as_negative_predictors": true
    },
    "cd4t": {
      "index_marker": "CD4",
      "predictors": ["CD3D", "CD3E", "CD3G", "CD4", "IL7R", "CD27", "SELL", "CCR7", "TCF7", "LEF1", "FOXP1", "BACH2", "KLF2", "CD28", "ICOS", "CD40LG", "STAT3", "STAT5", "RUNX1", "RUNX3", "BATF", "GATA3", "TBX21", "RORC", "FOXP3", "IL2RA", "CTLA4", "PDCD1", "LAG3", "TIGIT", "HAVCR2", "IL2", "IFNG", "IL4", "IL5", "IL13", "IL17A", "IL17F", "IL10", "CCL5", "CXCR3", "CCR6", "CXCR5", "PRDM1", "BCL6", "EOMES", "CXCL13", "TOX"],
      "use_offclass_markers_as_negative_predictors": true
    }
  }
}
```

### Hierarchy Levels

| Level | Name | Cell Types | Applies To |
|-------|------|------------|------------|
| L1 | Broad | epithelial, endothelial, fibroblast, plasma, immune, smooth_muscle | All cells |
| L2 | immunemajor | tcell, bcell, macrophage, neutrophil, mast, monocyte, dendritic, natural_killer_cell, eosinophil | L1 = "immune" only |
| LT | tcellmajor | cd8t, cd4t | L2 = "tcell" only |

---

## R Implementation

### Main Function (Production Version)

This is the corrected implementation that handles reticulate's automatic type conversions properly.

```r
annotate_hieratype_r <- function(
    h5ad_path,
    output_h5ad_path,
    marker_config_path,
    celltype_call_threshold = 0.5,
    use_spatial_graph = FALSE,
    k_neighbors = 20,
    remove_missing_marker_types = TRUE,
    key_added = "hieratype"
) {
    library(HieraType)
    library(Matrix)
    library(jsonlite)

    # ================================================================
    # STEP 1: LOAD DATA USING reticulate/anndata
    # ================================================================

    cat("Loading h5ad file...\n")

    # Use reticulate to load h5ad via Python anndata
    # Note: RETICULATE_PYTHON env var must be set by caller (r_bridge.py)
    library(reticulate)

    # Import Python modules
    ad <- import("anndata")
    scipy_sparse <- import("scipy.sparse")
    np <- import("numpy")

    # Read h5ad
    adata <- ad$read_h5ad(h5ad_path)

    n_cells <- as.integer(adata$n_obs)
    n_genes <- as.integer(adata$n_vars)
    cat(sprintf("  Loaded: %d cells x %d genes\n", n_cells, n_genes))

    # Get gene and cell names
    gene_names <- as.character(adata$var_names$to_list())
    cell_names <- as.character(adata$obs_names$to_list())

    # Extract expression matrix (cells x genes)
    # IMPORTANT: reticulate auto-converts scipy sparse to R Matrix classes
    X <- adata$X
    cat(sprintf("  Expression matrix class: %s\n", paste(class(X), collapse=", ")))

    if (inherits(X, "sparseMatrix")) {
        # X is already an R sparse matrix (dgRMatrix or dgCMatrix)
        counts_matrix <- as(X, "CsparseMatrix")
        rownames(counts_matrix) <- cell_names
        colnames(counts_matrix) <- gene_names
    } else if (scipy_sparse$issparse(X)) {
        # Fallback: manually extract if not auto-converted
        X_csc <- scipy_sparse$csc_matrix(X)
        indices <- as.integer(py_to_r(np$asarray(py_get_attr(X_csc, "indices"))))
        indptr <- as.integer(py_to_r(np$asarray(py_get_attr(X_csc, "indptr"))))
        data_vals <- as.numeric(py_to_r(np$asarray(py_get_attr(X_csc, "data"))))

        counts_matrix <- sparseMatrix(
            i = indices + 1L,
            p = indptr,
            x = data_vals,
            dims = c(n_cells, n_genes),
            dimnames = list(cell_names, gene_names)
        )
    } else {
        # Dense matrix
        counts_matrix <- as.matrix(py_to_r(np$asarray(X)))
        rownames(counts_matrix) <- cell_names
        colnames(counts_matrix) <- gene_names
    }

    # Validate log-normalization
    max_val <- max(counts_matrix@x)
    if (max_val > 20) {
        cat(sprintf("  WARNING: Data may not be log-normalized (max: %.2f)\n", max_val))
    }

  # ================================================================
  # STEP 2: LOAD & CONVERT MARKER CONFIGURATION
  # ================================================================

  cat("\nLoading marker configuration...\n")
  marker_json <- fromJSON(marker_config_path)

  # Convert JSON format to HieraType markerlist format
  convert_marker_level <- function(level_data) {
    markerlist <- lapply(names(level_data), function(ct_name) {
      ct_data <- level_data[[ct_name]]
      list(
        index_marker = ct_data$index_marker,
        predictors = ct_data$predictors,
        use_offclass_markers_as_negative_predictors =
          ifelse(is.null(ct_data$use_offclass_markers_as_negative_predictors),
                 TRUE, ct_data$use_offclass_markers_as_negative_predictors)
      )
    })
    names(markerlist) <- names(level_data)
    class(markerlist) <- c("list", "markerslist")
    return(markerlist)
  }

  l1_markers <- convert_marker_level(marker_json$l1)
  l2_markers <- convert_marker_level(marker_json$immunemajor)
  lt_markers <- convert_marker_level(marker_json$tcellmajor)

  cat(sprintf("  L1 types: %d\n", length(l1_markers)))
  cat(sprintf("  L2 types: %d\n", length(l2_markers)))
  cat(sprintf("  LT types: %d\n", length(lt_markers)))

    # ================================================================
    # STEP 3: VALIDATE AND FILTER MARKERS (Required for HieraType)
    # ================================================================

    cat("\nValidating and filtering marker genes...\n")

    # HieraType's validate_markerslist requires ALL specified markers to be present.
    # Filter marker lists to only include genes that exist in the data.
    filter_markers <- function(markerlist, level_name) {
        filtered_list <- list()

        for (ct_name in names(markerlist)) {
            ct_data <- markerlist[[ct_name]]
            index_markers <- ct_data$index_marker
            predictors <- ct_data$predictors

            # Filter to only markers present in data
            valid_index <- index_markers[index_markers %in% gene_names]
            valid_predictors <- predictors[predictors %in% gene_names]

            if (length(valid_index) == 0) {
                # No valid index markers - remove this cell type
                cat(sprintf("  REMOVING %s (%s): No index markers found in data\n",
                           ct_name, level_name))
            } else {
                # Update with filtered markers
                if (length(valid_index) < length(index_markers)) {
                    cat(sprintf("  %s (%s): Using %d/%d index markers\n",
                               ct_name, level_name, length(valid_index), length(index_markers)))
                }
                ct_data$index_marker <- valid_index
                ct_data$predictors <- valid_predictors
                filtered_list[[ct_name]] <- ct_data
            }
        }

        # Apply markerslist class
        class(filtered_list) <- c("list", "markerslist")
        return(filtered_list)
    }

    l1_markers <- filter_markers(l1_markers, "L1")
    l2_markers <- filter_markers(l2_markers, "L2")
    lt_markers <- filter_markers(lt_markers, "LT")

    cat(sprintf("  Filtered L1 types: %d\n", length(l1_markers)))
    cat(sprintf("  Filtered L2 types: %d\n", length(l2_markers)))
    cat(sprintf("  Filtered LT types: %d\n", length(lt_markers)))

  # ================================================================
  # STEP 4: BUILD HIERATYPE PIPELINE
  # ================================================================

  cat("\nBuilding HieraType pipeline...\n")
  pipeline <- make_pipeline(
    markerslists = list(
      "l1" = l1_markers,
      "l2" = l2_markers,
      "lt" = lt_markers
    ),
    priors = list(
      "lt" = "l2",      # LT refines L2 results
      "l2" = "l1"       # L2 refines L1 results
    ),
    priors_category = list(
      "lt" = "tcell",   # LT applies only to cells classified as "tcell" in L2
      "l2" = "immune"   # L2 applies only to cells classified as "immune" in L1
    )
  )

    # ================================================================
    # STEP 5: PREPARE ADJACENCY MATRIX
    # ================================================================

    cat("\nPreparing neighborhood graph...\n")

    if (use_spatial_graph) {
        # Option A: Compute spatial kNN from coordinates
        cat(sprintf("  Computing spatial kNN graph (k=%d)...\n", k_neighbors))
        library(FNN)

        # Extract spatial coordinates from adata.obsm
        if ("spatial" %in% names(py_to_r(adata$obsm))) {
            coords <- as.matrix(np$array(adata$obsm[["spatial"]]))
        } else if ("X_spatial" %in% names(py_to_r(adata$obsm))) {
            coords <- as.matrix(np$array(adata$obsm[["X_spatial"]]))
        } else {
            stop("No spatial coordinates found in adata.obsm['spatial'] or adata.obsm['X_spatial']")
        }

        # Compute kNN
        knn_result <- get.knn(coords, k = k_neighbors)

        # Build sparse symmetric adjacency matrix
        n <- nrow(coords)
        i_indices <- rep(1:n, each = k_neighbors)
        j_indices <- as.vector(t(knn_result$nn.index))
        adjacency_matrix <- sparseMatrix(
            i = i_indices, j = j_indices, x = 1, dims = c(n, n)
        )
        adjacency_matrix <- adjacency_matrix + t(adjacency_matrix)
        adjacency_matrix@x[adjacency_matrix@x > 0] <- 1
        adjacency_matrix <- as(adjacency_matrix, "CsparseMatrix")

    } else {
        # Option B: Use pre-computed SNN graph from clustering
        cat("  Using SNN graph from adata.obsp...\n")

        # Get obsp keys properly (iterate() converts Python iterator to R list)
        obsp_keys <- unlist(iterate(adata$obsp$keys()))
        cat(sprintf("  Available obsp keys: %s\n", paste(obsp_keys, collapse=", ")))

        if ("connectivities" %in% obsp_keys) {
            conn <- adata$obsp[["connectivities"]]
            # reticulate auto-converts scipy sparse to R Matrix classes
            if (inherits(conn, "sparseMatrix")) {
                adjacency_matrix <- as(conn, "CsparseMatrix")
            } else {
                # Fallback for non-auto-converted matrices
                conn_csc <- scipy_sparse$csc_matrix(conn)
                indices <- as.integer(py_to_r(np$asarray(py_get_attr(conn_csc, "indices"))))
                indptr <- as.integer(py_to_r(np$asarray(py_get_attr(conn_csc, "indptr"))))
                data_vals <- as.numeric(py_to_r(np$asarray(py_get_attr(conn_csc, "data"))))
                adjacency_matrix <- sparseMatrix(
                    i = indices + 1L,
                    p = indptr,
                    x = data_vals,
                    dims = c(n_cells, n_cells)
                )
            }
        } else if ("distances" %in% obsp_keys) {
            # Convert distances to similarities
            dist <- adata$obsp[["distances"]]
            if (inherits(dist, "sparseMatrix")) {
                dist_csc <- as(dist, "CsparseMatrix")
                # Convert distances to similarities: 1 / (1 + distance)
                dist_csc@x <- 1 / (1 + dist_csc@x)
                adjacency_matrix <- dist_csc
            } else {
                dist_csc <- scipy_sparse$csc_matrix(dist)
                indices <- as.integer(py_to_r(np$asarray(py_get_attr(dist_csc, "indices"))))
                indptr <- as.integer(py_to_r(np$asarray(py_get_attr(dist_csc, "indptr"))))
                data_vals <- as.numeric(py_to_r(np$asarray(py_get_attr(dist_csc, "data"))))
                adjacency_matrix <- sparseMatrix(
                    i = indices + 1L,
                    p = indptr,
                    x = 1 / (1 + data_vals),
                    dims = c(n_cells, n_cells)
                )
            }
        } else {
            stop("No neighborhood graph found in adata.obsp. Run sc.pp.neighbors() first.")
        }

        adjacency_matrix <- as(adjacency_matrix, "CsparseMatrix")
    }

    cat(sprintf("  Adjacency matrix: %d x %d, %d non-zeros\n",
                nrow(adjacency_matrix), ncol(adjacency_matrix),
                length(adjacency_matrix@x)))

    # ================================================================
    # STEP 6: RUN HIERATYPE PIPELINE
    # ================================================================

    cat("\nRunning HieraType annotation...\n")
    cat(sprintf("  Confidence threshold: %.2f\n", celltype_call_threshold))

    start_time <- Sys.time()

    # IMPORTANT: Pass counts_matrix as cells x genes (NOT transposed)
    # HieraType's fit_metagene_scores checks colnames for genes, and will
    # auto-transpose if genes are found in rownames instead
    hieratype_result <- run_pipeline(
        pipeline = pipeline,
        counts_matrix = counts_matrix,  # cells x genes, genes in colnames
        adjacency_matrix = adjacency_matrix,
        celltype_call_threshold = celltype_call_threshold,
        return_all_columns_postprobs = FALSE
    )

    elapsed_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat(sprintf("  Completed in %.1f seconds\n", elapsed_time))

  # ================================================================
  # STEP 7: EXTRACT RESULTS
  # ================================================================

  cat("\nExtracting annotation results...\n")

  extract_level_results <- function(level_name) {
    if (level_name %in% names(hieratype_result$post_probs)) {
      level_results <- as.data.frame(hieratype_result$post_probs[[level_name]])
      return(list(
        celltype = level_results$celltype_thresh,
        score = level_results$best_score_thresh
      ))
    } else {
      return(list(
        celltype = rep("Unknown", n_cells),
        score = rep(0, n_cells)
      ))
    }
  }

  l1_results <- extract_level_results("l1")
  l2_results <- extract_level_results("l2")
  lt_results <- extract_level_results("lt")

  # Add to adata.obs
  adata$obs$hieratype_l1 <- l1_results$celltype
  adata$obs$hieratype_l1_score <- l1_results$score
  adata$obs$hieratype_l2 <- l2_results$celltype
  adata$obs$hieratype_l2_score <- l2_results$score
  adata$obs$hieratype_lt <- lt_results$celltype
  adata$obs$hieratype_lt_score <- lt_results$score

  # Create final annotation (most specific level)
  adata$obs$hieratype_final <- ifelse(
    lt_results$celltype != "Unknown" & lt_results$celltype != "lt",
    lt_results$celltype,
    ifelse(
      l2_results$celltype != "Unknown" & l2_results$celltype != "l2",
      l2_results$celltype,
      l1_results$celltype
    )
  )
  adata$obs$hieratype_final_score <- pmax(
    l1_results$score, l2_results$score, lt_results$score
  )

  # ================================================================
  # STEP 8: SAVE AND RETURN
  # ================================================================

  cat("\nSaving annotated h5ad...\n")
  write_h5ad(adata, output_h5ad_path)

  # Return summary
  return(list(
    n_cells = n_cells,
    n_genes = n_genes,
    elapsed_time_seconds = elapsed_time,
    n_cell_types = length(unique(adata$obs$hieratype_final)),
    l1_counts = as.list(table(adata$obs$hieratype_l1)),
    l2_counts = as.list(table(adata$obs$hieratype_l2)),
    threshold = celltype_call_threshold
  ))
}
```

---

## Python Wrapper (Optional - for Calling R)

```python
import subprocess
import json
import tempfile
from pathlib import Path

def run_hieratype(
    h5ad_path: str,
    output_path: str,
    marker_config_path: str,
    celltype_call_threshold: float = 0.5,
    use_spatial_graph: bool = False,
    k_neighbors: int = 20,
) -> dict:
    """
    Run HieraType R function via subprocess.

    Args:
        h5ad_path: Path to input h5ad (normalized, with neighbors)
        output_path: Path for output h5ad
        marker_config_path: Path to marker JSON config
        celltype_call_threshold: Minimum confidence (0-1)
        use_spatial_graph: Use spatial kNN instead of SNN
        k_neighbors: k for spatial kNN

    Returns:
        dict with annotation summary
    """
    # Create R script
    r_script = f'''
    source("r_functions.R")
    result <- annotate_hieratype_r(
        h5ad_path = "{h5ad_path}",
        output_h5ad_path = "{output_path}",
        marker_config_path = "{marker_config_path}",
        celltype_call_threshold = {celltype_call_threshold},
        use_spatial_graph = {"TRUE" if use_spatial_graph else "FALSE"},
        nn_method = "{"spatial_knn" if use_spatial_graph else "snn"}",
        k_neighbors = {k_neighbors},
        remove_missing_marker_types = TRUE
    )
    cat(jsonlite::toJSON(result, auto_unbox = TRUE))
    '''

    with tempfile.NamedTemporaryFile(mode='w', suffix='.R', delete=False) as f:
        f.write(r_script)
        script_path = f.name

    try:
        result = subprocess.run(
            ["Rscript", script_path],
            capture_output=True,
            text=True,
            timeout=3600
        )

        if result.returncode != 0:
            raise RuntimeError(f"R script failed: {result.stderr}")

        # Parse JSON output from last line
        output_lines = result.stdout.strip().split('\n')
        json_output = output_lines[-1]
        return json.loads(json_output)

    finally:
        Path(script_path).unlink()
```

---

## Output Columns Created

| Column | Type | Description | Values |
|--------|------|-------------|--------|
| `hieratype_l1` | string | Broad cell type | epithelial, endothelial, fibroblast, plasma, immune, smooth_muscle, Unknown |
| `hieratype_l1_score` | float | L1 confidence | 0.0-1.0 |
| `hieratype_l2` | string | Immune subtype | tcell, bcell, macrophage, dendritic, neutrophil, mast, monocyte, natural_killer_cell, eosinophil, Unknown |
| `hieratype_l2_score` | float | L2 confidence | 0.0-1.0 |
| `hieratype_lt` | string | T cell subtype | cd8t, cd4t, Unknown |
| `hieratype_lt_score` | float | LT confidence | 0.0-1.0 |
| `hieratype_final` | string | Most specific | Best of LT > L2 > L1 |
| `hieratype_final_score` | float | Max score | 0.0-1.0 |

---

## Input Requirements

### Expression Matrix (`adata.X`)
- Sparse or dense matrix (cells × genes)
- **MUST be log-normalized** (max value < 20)
- Gene symbols as column names

### Neighborhood Graph (one of)

**Option A: SNN Graph** (default, from clustering)
- `adata.obsp['connectivities']`: Sparse adjacency matrix
- Created by `scanpy.pp.neighbors()` or `compute_neighbors()`

**Option B: Spatial kNN** (for spatial data)
- `adata.obsm['spatial']`: Dense array (cells × 2)
- X, Y coordinates in microns
- kNN computed by function using FNN package

---

## Configuration & Defaults

```python
# Default parameters
DEFAULT_THRESHOLD = 0.5        # Confidence cutoff
DEFAULT_K_NEIGHBORS = 20       # For spatial kNN
DEFAULT_USE_SPATIAL = False    # Use SNN by default

# Marker file location
DEFAULT_MARKER_CONFIG = "config/hieratype_markers.json"
```

---

## Dependencies

### R Packages
```r
# Core
HieraType           # Main classification package (install from GitHub)
Matrix              # Sparse matrix operations
anndata             # h5ad file reading

# Utilities
jsonlite            # JSON parsing
FNN                 # Fast nearest neighbor (for spatial kNN)
```

### Install HieraType
```r
# HieraType is not on CRAN - install from GitHub
devtools::install_github("path/to/HieraType")
```

### Python (for wrapper only)
```
subprocess          # Built-in
json                # Built-in
pathlib             # Built-in
```

---

## Usage Example

```r
# Standalone R usage
result <- annotate_hieratype_r(
    h5ad_path = "./data/normalized_clustered.h5ad",
    output_h5ad_path = "./data/hieratype_annotated.h5ad",
    marker_config_path = "./config/hieratype_markers.json",
    celltype_call_threshold = 0.5,
    use_spatial_graph = FALSE,  # Use SNN from clustering
    remove_missing_marker_types = TRUE
)

print(result$l1_counts)
# $epithelial: 5000
# $immune: 12000
# $fibroblast: 3000
# ...
```

```python
# Python wrapper usage
result = run_hieratype(
    h5ad_path="./data/normalized_clustered.h5ad",
    output_path="./data/hieratype_annotated.h5ad",
    marker_config_path="./config/hieratype_markers.json",
    celltype_call_threshold=0.5,
)

print(f"Annotated {result['n_cells']} cells into {result['n_cell_types']} types")
```

---

## Algorithm Details

### HieraType Scoring Method

1. **Metagene Construction**: For each cell type, compute metagene scores using non-negative least squares (NNLS) on predictor genes

2. **Neighborhood Averaging**: Average metagene scores across neighbors defined by adjacency matrix for robustness

3. **Cell Type Assignment**:
   - Compare scores across all cell types at current level
   - Assign cell type with highest score if score ≥ threshold
   - Otherwise assign "Unknown"

4. **Hierarchical Refinement**:
   - L1 runs on ALL cells
   - L2 runs ONLY on cells where L1 = "immune"
   - LT runs ONLY on cells where L2 = "tcell"

### Computational Complexity
- O(n_cells × n_neighbors × n_marker_genes)
- Memory: ~16GB for 100k cells

---

## When to Use HieraType

**Good For:**
- Spatial transcriptomics with tissue coordinates
- Immune cell classification (strong built-in markers)
- Hierarchical annotations (broad → fine-grained)
- Leveraging spatial/SNN neighborhood context
- Custom marker hierarchies

**Not Recommended For:**
- Non-spatial scRNA-seq without clustering (no meaningful neighbors)
- Very fast annotation (use CellTypist instead)
- No marker knowledge (run find_markers first)
- Tissues not covered by default markers (need custom config)
