# CellTypist Cell Typing Pipeline Specification

## Purpose

This document provides a complete specification for rebuilding the CellTypist cell typing pipeline as standalone Python scripts without Temporal/process_activity dependencies. Use this to create a stepwise data processing pipeline.

---

## Complete Data Flow

```
┌─────────────────────────────────────────────────────────────────┐
│ STEP 1: DOWNLOAD CELLXGENE REFERENCES                           │
│ Input:  Dataset key (e.g., "colon_immune_niches")              │
│ Output: CellxGene h5ad (30k+ cells × 15k+ genes)               │
│ Key:    Gene symbols in var['feature_name'] → normalize        │
└─────────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 2: COMBINE & VALIDATE REFERENCES                           │
│ Input:  Multiple reference h5ad files + label columns          │
│ Process: Subsample → Normalize → Harmonize labels → Intersect  │
│ Output: Combined reference (100k cells × shared genes)         │
└─────────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 3: TRAIN CUSTOM CELLTYPIST MODEL                           │
│ Input:  Combined reference + panel gene list                   │
│ Process: Subset to panel genes → Train logistic regression     │
│ Output: .pkl model file + metadata                             │
└─────────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 4: ANNOTATE QUERY DATA                                     │
│ Input:  Query h5ad (normalized) + tissue type                  │
│ Process: Multi-model ensemble → Batched prediction → Voting    │
│ Output: adata.obs['cell_type', 'cell_type_confidence', ...]    │
└─────────────────────────────────────────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────────────┐
│ STEP 5: CONVERT TO ONTOLOGY CODES                               │
│ Input:  cell_type labels (text)                                │
│ Process: Pattern matching → Multi-tier ontology search         │
│ Output: CL IDs (e.g., "CL:0000236")                            │
└─────────────────────────────────────────────────────────────────┘
```

---

## Reference Data Sources

CellTypist training supports three reference data sources, each with different access patterns and data formats:

### Overview

| Source | Format | Access | Use Case |
|--------|--------|--------|----------|
| **CellxGene Census** | h5ad | Public API | Published atlases, broad cell types |
| **Synapse** | h5ad/RDS | Token auth | Controlled-access data (AMP, HTAN) |
| **Local scRNAseq** | Seurat RDS | File system | Pre-curated institutional data |

### Environment Configuration

Configure data paths via environment variables:

```bash
# Optional: CellxGene download cache
export CELLXGENE_CACHE_DIR="/path/to/cellxgene/cache"

# Required for Synapse downloads
export SYNAPSE_AUTH_TOKEN="your_synapse_personal_access_token"

# Required for local scRNAseq data
export SCRNASEQ_DATA_DIR="/path/to/scrnaseq/data"
```

---

### Source 1: CellxGene Census (Public)

CellxGene Census provides public single-cell datasets downloadable via API. No authentication required.

**Registry structure:**
```python
from spacebio.celltypist_utils import CELLXGENE_DATASETS

# Available datasets with metadata
CELLXGENE_DATASETS = {
    "colon_immune_niches": {
        "dataset_id": "2872f4b0-b171-46e2-abc6-befcf6de6306",
        "description": "Distinct microbial and immune niches of the human colon",
        "tissue": "colon",
        "cell_type_column": "cell_type",
    },
    "human_lung_cell_atlas": {...},
    # ... more datasets
}
```

**Download example:**
```python
from spacebio.celltypist_utils import download_cellxgene_reference, load_reference

# Direct download
h5ad_path = download_cellxgene_reference("colon_immune_niches", output_dir=Path("./references"))

# Or via unified interface
adata = load_reference("cellxgene", "colon_immune_niches", output_dir=Path("./references"))
```

---

### Source 2: Synapse (Controlled Access)

Synapse provides controlled-access datasets requiring authentication. Some datasets require explicit access approval.

**Setup:**
1. Create Synapse account at https://synapse.org
2. Generate personal access token: Account Settings → Personal Access Tokens
3. Set environment variable: `export SYNAPSE_AUTH_TOKEN="your_token"`
4. For restricted datasets, request access via Synapse web interface

**Registry structure:**
```python
from spacebio.celltypist_utils import SYNAPSE_DATASETS

SYNAPSE_DATASETS = {
    "amp_ra_synovial_atlas": {
        "synapse_id": "syn47579773",
        "description": "AMP Phase II RA Synovial Atlas (~314k cells)",
        "tissue": "synovium",
        "condition": "rheumatoid arthritis",
        "cell_type_column": "cell_type",
        "requires_approval": True,
    },
}
```

**Download example:**
```python
from spacebio.celltypist_utils import download_synapse_reference, load_reference

# Direct download (uses SYNAPSE_AUTH_TOKEN env var)
h5ad_path = download_synapse_reference("syn47579773", output_dir=Path("./references"))

# Or via unified interface
adata = load_reference("synapse", "syn47579773", output_dir=Path("./references"))
```

---

### Source 3: Local scRNAseq Data

Local scRNAseq data consists of pre-curated Seurat objects with standardized metadata CSVs for discovery.

**Directory structure:**
```
$SCRNASEQ_DATA_DIR/
├── metadata_summary.csv      # Study-level metadata
├── sample_summary.csv        # Sample-level with cell type breakdown
└── seurat_objects/
    ├── 21.curated.seurat.rds
    ├── 42.curated.seurat.rds
    └── ...
```

**Metadata schema:**

`metadata_summary.csv` columns:
| Column | Description |
|--------|-------------|
| `study_id` | Unique ID, maps to `{study_id}.curated.seurat.rds` |
| `condition_standardized` | Disease/condition (e.g., "rheumatoid arthritis") |
| `tissue_standardized` | Tissue type (e.g., "synovial membrane of synovial joint") |
| `sample_site_standardized` | Location (e.g., "lesional", "normal") |
| `species_standardized` | Organism (e.g., "homo sapiens") |
| `n_subjects` | Number of donors |
| `n_samples` | Number of samples |
| `n_cells` | Total cell count |

`sample_summary.csv` columns:
| Column | Description |
|--------|-------------|
| `study_id` | Links to metadata_summary |
| `sample_id` | Individual sample identifier |
| `celltype` | Cell type annotation |
| `n` | Cell count for this sample/celltype |

**Query presets:**
```python
from spacebio.celltypist_utils import LOCAL_PRESETS

LOCAL_PRESETS = {
    "ra_synovial": {
        "description": "RA synovial tissue/fluid studies",
        "condition": "rheumatoid arthritis",
        "tissue_pattern": "synovial",
        "min_cells": 10000,
    },
    "arthritis_all": {
        "description": "All arthritis types (RA, PsA, AS) synovial",
        "condition_pattern": "arthritis|spondylitis",
        "tissue_pattern": "synovial",
        "min_cells": 5000,
    },
}
```

**Query local data:**
```python
from spacebio.celltypist_utils import (
    load_local_metadata,
    query_local_references,
    query_local_preset,
    load_reference,
)
from pathlib import Path

# Load metadata CSVs
metadata_df, sample_df = load_local_metadata(
    metadata_csv=Path("metadata_summary.csv"),
    sample_csv=Path("sample_summary.csv"),
)

# Query by condition/tissue
ra_studies = query_local_references(
    metadata_df,
    condition="rheumatoid arthritis",
    tissue="synovial",
    min_cells=10000,
)
print(f"Found {len(ra_studies)} RA synovial studies")

# Or use preset
ra_studies = query_local_preset("ra_synovial", metadata_df)

# Load Seurat as AnnData (auto-converts and caches)
adata = load_reference(
    "local",
    study_id="42",
    seurat_dir=Path("seurat_objects"),
    metadata_csv=Path("metadata_summary.csv"),
)
```

**Cell type column auto-detection:**

When loading local data, the system auto-detects cell type columns in this priority order:
1. `celltype_lvl_1_standardized`
2. `cell_type`
3. `celltype`
4. `CellType`
5. `cell_type_ontology_term_id`
6. `cluster`
7. `seurat_clusters`

If no column is detected, available columns are printed for manual selection:
```python
adata = load_reference("local", "42", celltype_column="my_celltype_annotation", ...)
```

---

### Unified Reference Loading

The `load_reference()` function provides a unified interface to load data from any source:

```python
from spacebio.celltypist_utils import load_reference

# CellxGene
adata = load_reference("cellxgene", "colon_immune_niches")

# Synapse
adata = load_reference("synapse", "syn47579773")

# Local
adata = load_reference(
    "local",
    "42",
    seurat_dir=Path("seurat_objects"),
    metadata_csv=Path("metadata_summary.csv"),
)
```

---

## STEP 1: Download CellxGene References

### Dataset Registry

```python
RECOMMENDED_DATASETS = {
    # Colon/GI datasets
    "colon_immune_niches": {
        "description": "Distinct microbial and immune niches of the human colon",
        "tissue": "colon",
        "dataset_id": "2872f4b0-b171-46e2-abc6-befcf6de6306",
        "cell_type_column": "cell_type",
        "expected_cells": "~41,650",
        "manual_url": "https://cellxgene.cziscience.com/collections/dde06e0f-ab3b-46be-96a2-a8082383c4a1",
    },
    "colon_ulcerative_colitis": {
        "description": "Human Colon during Ulcerative Colitis (Smillie et al.)",
        "tissue": "colon",
        "dataset_id": "4dd00779-7f73-4f50-89bb-e2d3c6b71b18",
        "cell_type_column": "cell_type",
        "expected_cells": "~34,772",
    },
    "colon_crohns_immune": {
        "description": "Crohn's disease colon immune cells",
        "tissue": "colon",
        "dataset_id": "518d9049-2a76-44f8-8abc-1e2b59ab5ba1",
        "cell_type_column": "cell_type",
        "expected_cells": "~152,509",
    },
    # Lung datasets
    "human_lung_cell_atlas": {
        "description": "Human Lung Cell Atlas (HLCA) - Azimuth",
        "tissue": "lung",
        "dataset_id": "f72958f5-7f42-4ebb-98da-445b0c6de516",
        "cell_type_column": "ann_finest_level",  # HLCA hierarchical
        "expected_cells": "~584,884",
    },
    "lung_covid": {
        "description": "Molecular single-cell lung atlas of lethal COVID-19",
        "tissue": "lung",
        "dataset_id": "d8da613f-e681-4c69-b463-e94f5e66847f",
        "cell_type_column": "cell_type",
        "expected_cells": "~116,313",
    },
    # Immune datasets
    "tabula_sapiens_blood": {
        "description": "Tabula Sapiens - Blood (immune cells)",
        "tissue": "immune",
        "dataset_id": None,  # Search by terms
        "cell_type_column": "cell_type",
        "expected_cells": "~100,000+",
    },
    # CRC datasets
    "crc_htan_epithelial_discovery": {
        "description": "HTAN VUMC CRC Polyps - Epithelial (Discovery)",
        "tissue": "crc",
        "dataset_id": "e40c6272-af77-4a10-9385-62a398884f27",
        "cell_type_column": "cell_type",
        "expected_cells": "~65,088",
    },
    # ... additional datasets available
}
```

### Gene Name Normalization (CRITICAL)

CellxGene Census stores gene identifiers in various formats. The `var_names` may contain:
- **Numeric indices** (0, 1, 2, ...) with gene symbols in `var['feature_name']`
- **Ensembl IDs** (ENSG00000000003, ...) with symbols in `var['feature_name']`
- **Mixed content** in `feature_name` (some Ensembl IDs, some HUGO symbols)

**Critical Finding:** Some CellxGene datasets (e.g., COVID lung atlas) have Ensembl IDs in `var_names` AND mixed Ensembl/HUGO content in `feature_name`. A robust gene mapping step is essential.

#### Ensembl to HUGO Gene Mapping

```python
import pandas as pd
import numpy as np

# Gene mapping file location (pre-built from Ensembl BioMart)
GENE_MAPPING_FILE = Path("data_for_docker_dep/gene_mappings/ensembl_to_hugo_human.tsv")

def load_ensembl_to_hugo_mapping() -> dict:
    """Load Ensembl ID to HUGO gene symbol mapping from TSV file."""
    if not GENE_MAPPING_FILE.exists():
        logger.warning(f"Gene mapping file not found: {GENE_MAPPING_FILE}")
        return {}

    df = pd.read_csv(GENE_MAPPING_FILE, sep="\t")
    df = df.dropna(subset=["HGNC symbol"])
    df = df[df["HGNC symbol"].str.len() > 0]
    mapping = dict(zip(df["Gene stable ID"], df["HGNC symbol"]))
    logger.info(f"Loaded {len(mapping):,} Ensembl to HUGO gene mappings")
    return mapping


def is_ensembl_id(gene_name: str) -> bool:
    """Check if a gene name looks like an Ensembl ID."""
    if not gene_name or not isinstance(gene_name, str):
        return False
    return (
        gene_name.startswith("ENSG") or      # Human gene
        gene_name.startswith("ENST") or      # Human transcript
        gene_name.startswith("ENSMUSG") or   # Mouse gene
        gene_name.startswith("ENSMUS")       # Mouse
    )


def convert_ensembl_to_hugo(
    gene_names: np.ndarray,
    ensembl_to_hugo: dict
) -> tuple[np.ndarray, dict]:
    """
    Convert Ensembl IDs to HUGO gene symbols where possible.

    SAFE for HUGO symbols: passes them through unchanged.

    Returns:
        (converted_names, stats_dict)
    """
    converted = []
    n_converted = 0
    n_already_hugo = 0
    n_unmapped = 0

    for gene in gene_names:
        gene_str = str(gene)
        if is_ensembl_id(gene_str):
            if gene_str in ensembl_to_hugo:
                converted.append(ensembl_to_hugo[gene_str])
                n_converted += 1
            else:
                # Keep unmapped Ensembl ID (will be filtered during panel subsetting)
                converted.append(gene_str)
                n_unmapped += 1
        else:
            # Already HUGO symbol - pass through
            converted.append(gene_str)
            n_already_hugo += 1

    stats = {
        "total_genes": len(gene_names),
        "converted_ensembl": n_converted,
        "already_hugo": n_already_hugo,
        "unmapped_ensembl": n_unmapped,
    }
    return np.array(converted), stats
```

#### Gene Name Normalization Function

```python
def normalize_gene_names(adata: sc.AnnData, ensembl_to_hugo: dict = None) -> sc.AnnData:
    """
    Normalize gene names in AnnData to use HUGO gene symbols.

    Two-step process:
    1. If var_names are Ensembl IDs, use feature_name column as starting point
    2. Apply Ensembl→HUGO mapping for any remaining Ensembl IDs

    Safe to call on data that already uses HUGO symbols.
    """
    first_gene = str(adata.var_names[0])
    uses_non_symbol_ids = (
        first_gene.isdigit() or
        first_gene.startswith("ENSG") or
        first_gene.startswith("ENST")
    )

    if not uses_non_symbol_ids:
        return adata  # Already uses symbols

    # Step 1: Use feature_name column if available
    if "feature_name" in adata.var.columns:
        feature_names = adata.var["feature_name"].values.astype(str)
        adata.var_names = pd.Index(feature_names)
        print(f"  Using 'feature_name' column as gene names")

    # Step 2: Apply Ensembl to HUGO mapping for any remaining Ensembl IDs
    if ensembl_to_hugo:
        converted_names, stats = convert_ensembl_to_hugo(
            adata.var_names.values, ensembl_to_hugo
        )
        adata.var_names = pd.Index(converted_names)

        if stats['converted_ensembl'] > 0:
            print(f"  Applied gene mapping: {stats['converted_ensembl']:,} converted, "
                  f"{stats['already_hugo']:,} already HUGO, {stats['unmapped_ensembl']:,} unmapped")
        else:
            print(f"  All {stats['already_hugo']:,} genes already HUGO symbols")

    adata.var_names_make_unique()
    return adata
```

**Key Points:**
- **Always apply mapping**: HUGO symbols pass through unchanged, so it's safe to always apply
- **Unmapped Ensembl IDs don't hurt training**: They're filtered out during panel gene subsetting since they won't match HUGO panel gene names
- **Load mapping once**: Load the Ensembl→HUGO mapping at the start of `combine_references()` and reuse for all datasets

### Download Function

```python
def download_via_census_api(dataset_key: str, output_dir: Path) -> Optional[Path]:
    """Download dataset using CellxGene Census API."""
    import cellxgene_census

    dataset_info = RECOMMENDED_DATASETS[dataset_key]
    dataset_id = dataset_info.get("dataset_id")

    if not dataset_id:
        # Search by terms
        dataset_id = search_dataset_id(dataset_info["search_terms"])

    output_file = output_dir / f"{dataset_key}.h5ad"

    # Download original h5ad
    cellxgene_census.download_source_h5ad(dataset_id, to_path=str(output_file))

    # Normalize gene names
    normalize_gene_names(output_file)

    return output_file
```

---

## STEP 2: Combine & Validate References

### Memory-Efficient Loading with Backed Mode

Large CellxGene reference datasets (e.g., HLCA with 584k cells, 5.6GB) can exceed available memory when loaded fully. Use **backed mode** for memory-efficient processing.

```python
import scanpy as sc
import psutil
from pathlib import Path

def get_available_memory_gb() -> float:
    """Get available system memory in GB."""
    return psutil.virtual_memory().available / (1024**3)


def load_reference_memory_efficient(
    path: str,
    max_cells: int,
    label_column: str,
    ensembl_to_hugo: dict = None,
    large_file_threshold_gb: float = 2.0,
) -> sc.AnnData:
    """
    Load reference data with memory-efficient strategies for large files.

    Strategy:
    - Small files (<2GB): Load fully into memory
    - Large files (>=2GB): Use backed mode, subsample indices, then load subset

    Returns:
        AnnData loaded into memory (not backed)
    """
    file_path = Path(path.replace("file://", ""))
    file_size_gb = file_path.stat().st_size / (1024**3)
    available_memory = get_available_memory_gb()

    print(f"  Available memory: {available_memory:.1f} GB")

    if file_size_gb >= large_file_threshold_gb:
        print(f"  Large file ({file_size_gb:.1f}GB) - using backed mode for memory efficiency")

        # Step 1: Open in backed mode (only loads metadata, not expression data)
        adata_backed = sc.read_h5ad(path, backed='r')
        print(f"  Opened: {adata_backed.n_obs:,} cells × {adata_backed.n_vars:,} genes (backed)")

        # Step 2: Determine subsample indices BEFORE loading data
        n_cells = adata_backed.n_obs
        if n_cells > max_cells:
            # Stratified sampling by cell type
            if label_column in adata_backed.obs.columns:
                indices = _stratified_sample_indices(
                    adata_backed.obs, label_column, max_cells
                )
            else:
                indices = np.random.choice(n_cells, size=max_cells, replace=False)
            print(f"  Subsampling {n_cells:,} → {len(indices):,} cells")
        else:
            indices = np.arange(n_cells)

        # Step 3: Load only the selected cells into memory
        # This is the key memory optimization - we never load the full matrix
        adata = adata_backed[indices].to_memory()
        print(f"  Loaded into memory: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

        # Clean up backed reference
        adata_backed.file.close()

    else:
        # Small file - load directly
        print(f"  Loading: {file_size_gb:.2f}GB file")
        adata = sc.read_h5ad(path)

        # Subsample if needed
        if adata.n_obs > max_cells:
            adata = subsample_adata(adata, max_cells, stratify_by=label_column)

    # Normalize gene names (Ensembl → HUGO)
    adata = normalize_gene_names(adata, ensembl_to_hugo)

    return adata


def _stratified_sample_indices(
    obs_df: pd.DataFrame,
    label_column: str,
    max_cells: int,
    random_state: int = 42,
) -> np.ndarray:
    """
    Get stratified sample indices maintaining cell type proportions.

    Works with backed AnnData by only accessing .obs metadata.
    """
    np.random.seed(random_state)

    labels = obs_df[label_column].values
    unique_labels = np.unique(labels)

    # Calculate cells per type (proportional)
    label_counts = pd.Series(labels).value_counts()
    total_cells = len(labels)
    indices = []

    for label in unique_labels:
        label_indices = np.where(labels == label)[0]
        # Proportional allocation
        n_sample = int(np.ceil(max_cells * len(label_indices) / total_cells))
        n_sample = min(n_sample, len(label_indices))

        sampled = np.random.choice(label_indices, size=n_sample, replace=False)
        indices.extend(sampled)

    # Trim to exact max_cells if we oversampled due to ceiling
    indices = np.array(indices)
    if len(indices) > max_cells:
        indices = np.random.choice(indices, size=max_cells, replace=False)

    return np.sort(indices)
```

**Key Memory Optimization:**
- `sc.read_h5ad(path, backed='r')` opens the file without loading the expression matrix
- `adata_backed[indices].to_memory()` loads ONLY the selected cells
- For HLCA (584k cells × 27k genes = ~60GB dense), we load only 100k cells (~10GB)

**Real-world Example:**
```
[1/2] Loading: human_lung_cell_atlas.h5ad
  Available memory: 30.5 GB
  Large file (5.6GB) - using backed mode for memory efficiency
  Opened: 584,884 cells × 27,402 genes (backed)
  Subsampling 584,884 → 100,000 cells
  Loaded into memory: 100,000 cells × 27,402 genes
```

### Normalization Detection

```python
def check_normalization_status(adata: sc.AnnData, sample_size: int = 1000) -> dict:
    """
    Detect normalization type.

    Returns:
        {
            "is_normalized": bool,
            "is_log_transformed": bool,
            "normalization_type": "raw" | "log1p_10k" | "normalized_not_logged" | "unknown",
            "stats": {"mean": float, "max": float, "min": float, "fraction_integer": float}
        }
    """
    from scipy.sparse import issparse

    n_sample = min(sample_size, adata.n_obs)
    if issparse(adata.X):
        sample_data = adata.X[:n_sample].toarray()
    else:
        sample_data = adata.X[:n_sample]

    mean_val = float(np.mean(sample_data))
    max_val = float(np.max(sample_data))
    min_val = float(np.min(sample_data))

    # Check for integer counts (raw data indicator)
    flat_sample = sample_data.flatten()[:10000]
    non_zero = flat_sample[flat_sample != 0]
    fraction_integer = np.mean(np.equal(np.mod(non_zero, 1), 0)) if len(non_zero) > 0 else 1.0

    stats = {"mean": mean_val, "max": max_val, "min": min_val, "fraction_integer": fraction_integer}

    # Decision logic
    if fraction_integer > 0.9 and max_val > 100:
        return {"is_normalized": False, "is_log_transformed": False,
                "normalization_type": "raw", "stats": stats}

    if min_val >= 0 and mean_val < 5 and max_val < 20:
        return {"is_normalized": True, "is_log_transformed": True,
                "normalization_type": "log1p_10k", "stats": stats}

    if mean_val < 10 and (max_val > 20 or min_val < 0):
        return {"is_normalized": True, "is_log_transformed": False,
                "normalization_type": "normalized_not_logged", "stats": stats}

    return {"is_normalized": False, "is_log_transformed": False,
            "normalization_type": "unknown", "stats": stats}
```

### Ensure Normalization (CellTypist Requirement)

```python
def ensure_normalized(adata: sc.AnnData, target_sum: float = 1e4) -> sc.AnnData:
    """
    Ensure data is log-normalized to 10k counts/cell.
    CellTypist REQUIRES: log1p(CPM to 10k) with exclude_highly_expressed=False
    """
    status = check_normalization_status(adata)

    if status["normalization_type"] == "log1p_10k":
        return adata  # Already normalized

    if status["normalization_type"] == "raw":
        # CRITICAL: exclude_highly_expressed=False for CellTypist compatibility
        sc.pp.normalize_total(adata, target_sum=target_sum, exclude_highly_expressed=False)
        sc.pp.log1p(adata)
        return adata

    if status["normalization_type"] == "normalized_not_logged":
        sc.pp.log1p(adata)
        return adata

    raise ValueError(f"Unknown normalization state: {status}")
```

### Subsample with Stratification

```python
def subsample_adata(
    adata: sc.AnnData,
    max_cells: int,
    stratify_by: Optional[str] = None,
    random_state: int = 42,
) -> sc.AnnData:
    """Subsample AnnData to max_cells, maintaining cell type proportions."""
    if adata.n_obs <= max_cells:
        return adata

    np.random.seed(random_state)

    if stratify_by and stratify_by in adata.obs.columns:
        indices = []
        groups = adata.obs[stratify_by].unique()
        n_per_group = max_cells // len(groups)

        for group in groups:
            group_idx = np.where(adata.obs[stratify_by] == group)[0]
            n_sample = min(len(group_idx), n_per_group)
            sampled = np.random.choice(group_idx, size=n_sample, replace=False)
            indices.extend(sampled)

        # Fill remaining with random cells
        remaining = max_cells - len(indices)
        if remaining > 0:
            available = list(set(range(adata.n_obs)) - set(indices))
            if available:
                extra = np.random.choice(available, size=min(remaining, len(available)), replace=False)
                indices.extend(extra)

        indices = np.array(indices)
    else:
        indices = np.random.choice(adata.n_obs, size=max_cells, replace=False)

    return adata[indices].copy()
```

### Combine Multiple References

```python
def combine_references(
    reference_paths: List[str],
    label_columns: List[str],
    output_column: str = "unified_cell_type",
    max_cells_per_ref: int = 500000,
    target_genes: Optional[List[str]] = None,
    harmonize_labels: bool = True,
    normalize_data: bool = True,
) -> sc.AnnData:
    """
    Combine multiple reference datasets for CellTypist training.

    Key steps:
    1. Load Ensembl→HUGO gene mapping (once, reused for all references)
    2. Load each reference with memory-efficient backed mode
    3. Normalize gene names (Ensembl → HUGO via feature_name + mapping)
    4. Subsample if > max_cells_per_ref (stratified by cell type)
    5. Check/apply normalization
    6. Harmonize cell type labels via Cell Ontology
    7. Subset to target panel genes
    8. Find shared genes (inner join)
    9. Concatenate with memory checks
    """
    import gc

    # Load gene mapping once at start (reused for all references)
    ensembl_to_hugo = load_ensembl_to_hugo_mapping()
    if ensembl_to_hugo:
        print(f"Loaded {len(ensembl_to_hugo):,} Ensembl to HUGO gene mappings")

    adatas = []

    for i, (path, label_col) in enumerate(zip(reference_paths, label_columns)):
        print(f"\n[{i+1}/{len(reference_paths)}] Loading: {Path(path).name}")

        # Memory-efficient loading with backed mode for large files
        adata = load_reference_memory_efficient(
            path=path,
            max_cells=max_cells_per_ref,
            label_column=label_col,
            ensembl_to_hugo=ensembl_to_hugo,  # Pass mapping for gene normalization
        )

        # Check normalization status
        if normalize_data:
            status = check_normalization_status(adata)
            if status["normalization_type"] == "log1p_10k":
                print(f"  ✓ adata.X already log-normalized")
            else:
                adata = ensure_normalized(adata)
                print(f"  Applied log1p(10k) normalization")

        # Harmonize labels via Cell Ontology
        if harmonize_labels:
            unique_labels = adata.obs[label_col].nunique()
            print(f"  Harmonizing {unique_labels} unique cell type labels...")
            adata = harmonize_cell_type_labels(adata, label_col, output_column)
            unified_labels = adata.obs[output_column].nunique()
            print(f"  Unified to {unified_labels} harmonized cell types")
        else:
            adata.obs[output_column] = adata.obs[label_col]

        adata.obs["reference_source"] = Path(path).stem
        adatas.append(adata)
        gc.collect()

    # Subset to target genes if provided (BEFORE finding shared genes)
    if target_genes:
        print(f"\nSubsetting to {len(target_genes)} target genes...")
        target_set = set(target_genes)
        for i, adata in enumerate(adatas):
            overlap = list(set(adata.var_names) & target_set)
            adatas[i] = adata[:, overlap].copy()
            print(f"  Reference {i}: {len(overlap)} genes after subset")

    # Find shared genes (inner join)
    print("\nFinding shared genes across all references...")
    shared_genes = set(adatas[0].var_names)
    for adata in adatas[1:]:
        shared_genes &= set(adata.var_names)

    if len(shared_genes) == 0:
        raise ValueError("No shared genes found across references!")

    print(f"  Shared genes: {len(shared_genes)}")

    # Subset all to shared genes (sorted for consistency)
    shared_genes_sorted = sorted(shared_genes)
    for i in range(len(adatas)):
        adatas[i] = adatas[i][:, shared_genes_sorted].copy()

    # Memory check before concatenation
    available_gb = get_available_memory_gb()
    total_cells = sum(a.n_obs for a in adatas)
    estimated_gb = (total_cells * len(shared_genes) * 4) / (1024**3)  # float32
    print(f"\nMemory check before concat:")
    print(f"  Available: {available_gb:.1f} GB")
    print(f"  Estimated need: ~{estimated_gb:.1f} GB")

    # Concatenate
    print("\nConcatenating references...")
    combined = sc.concat(adatas, join="inner", label="batch", index_unique="-")
    print(f"  ✓ Combined: {combined.n_obs:,} cells × {combined.n_vars:,} genes")

    # Print cell type distribution
    print(f"\n  Cell type distribution:")
    ct_counts = combined.obs[output_column].value_counts()
    for ct, count in ct_counts.head(10).items():
        print(f"    {ct}: {count:,} cells")
    if len(ct_counts) > 10:
        print(f"    ... and {len(ct_counts) - 10} more types")

    return combined
```

**Key Improvements:**
1. **Gene mapping loaded once**: Ensembl→HUGO mapping loaded at start and passed to all `load_reference_memory_efficient()` calls
2. **Memory-efficient loading**: Uses backed mode for large files (>2GB)
3. **Panel gene subsetting FIRST**: Subset to target genes before finding shared genes (reduces memory)
4. **Memory checks**: Estimates memory needed before concatenation
5. **Cell type harmonization**: Reports unified cell type counts

---

## STEP 3: Train Custom CellTypist Model

### Training Function

```python
def train_celltypist_model(
    adata: sc.AnnData,
    label_column: str = "unified_cell_type",
    model_name: str = "custom_model",
    output_path: str = "./model.pkl",
    use_SGD: bool = True,
    feature_selection: bool = False,
    n_jobs: int = -1,
    max_iter: int = 100,
) -> dict:
    """
    Train a custom CellTypist logistic regression model.

    Args:
        adata: Reference data (already subset to panel genes, normalized)
        label_column: Cell type label column
        model_name: Name for the model
        output_path: Where to save .pkl file
        use_SGD: Use stochastic gradient descent (faster for large data)
        feature_selection: Perform feature selection (False = use all genes)
        n_jobs: Parallel jobs (-1 = all cores)
        max_iter: Maximum training iterations

    Returns:
        dict with model metadata
    """
    import celltypist

    # Ensure normalized (double-check)
    status = check_normalization_status(adata)
    if status["normalization_type"] != "log1p_10k":
        adata = ensure_normalized(adata)

    # Train model
    model = celltypist.train(
        adata,
        labels=label_column,
        check_expression=False,  # Already validated
        use_SGD=use_SGD,
        feature_selection=feature_selection,
        n_jobs=n_jobs,
        max_iter=max_iter,
    )

    # Save model
    model.write(output_path)

    return {
        "model_path": output_path,
        "n_cells_trained": adata.n_obs,
        "n_genes": len(model.features),
        "n_cell_types": len(model.cell_types),
        "cell_types": model.cell_types.tolist(),
    }
```

### Training Data Balancing Strategies

CellTypist training can suffer from class imbalance (e.g., 192k T cells vs 277 mast cells). Two approaches are available:

#### Option 1: CellTypist Built-in Balancing (Recommended)

CellTypist's `train()` function supports mini-batch SGD with cell type balancing:

```python
import celltypist

model = celltypist.train(
    adata,
    labels="celltype_lvl_1_standardized",
    check_expression=False,
    use_SGD=True,
    mini_batch=True,           # Enable mini-batch training
    balance_cell_type=True,    # Balance rare cell types in batches
    epochs=10,                 # Multiple training passes
    batch_number=200,          # Batches per epoch
    batch_size=1000,           # Cells per batch
)
```

**How it works:** `balance_cell_type=True` samples rare cell types with higher probability, ensuring close-to-even distributions in mini-batches. This is the recommended approach for large datasets.

#### Option 2: Pre-training Subsampling

For more control, subsample before training:

```python
from spacebio.celltypist_utils import subsample_balanced

# Cap abundant types, keep all rare types
adata_balanced = subsample_balanced(
    adata,
    label_column="celltype_lvl_1_standardized",
    max_cells_per_type=50000,   # Cap at 50k cells per type
    min_cells_per_type=100,     # Exclude types with <100 cells
    random_state=42,
)

# Train on balanced data
model = celltypist.train(adata_balanced, labels="celltype_lvl_1_standardized")
```

#### Combining Multiple Presets

For broader models, combine multiple data presets:

```python
from spacebio.celltypist_utils import combine_presets, load_local_metadata

metadata_df, _ = load_local_metadata()

# Combine RA-only and all-arthritis presets
combined_df = combine_presets(["ra_synovial", "arthritis_all"], metadata_df)
# Returns deduplicated studies from both presets
```

**CLI Example:**
```bash
python scripts/celltypist_training/train_celltypist_RA.py \
    --presets ra_synovial arthritis_all \
    --mini-batch \
    --balance-cell-types \
    --epochs 10
```

### Cell Type Hierarchy

Local scRNAseq data uses multiple annotation levels:

| Column | Description | Example Types |
|--------|-------------|---------------|
| `celltype_lvl_1_standardized` | Major lineages (broad) | T cell, fibroblast, mononuclear phagocyte |
| `celltype_lvl_2_standardized` | Subtypes (if available) | CD4+ T cell, CD8+ T cell, Treg |
| `author_celltype` | Study-specific (not harmonized) | Varies by study |

**Choosing annotation level:**
- **Level 1**: Best for spatial biology (matches spatial resolution)
- **Level 2**: Use for finer subtyping if available
- **Author**: Use for study-specific fine types (requires harmonization)

---

## STEP 4: Annotate Query Data (CellTypist Ensemble)

### Tissue Model Presets

```python
TISSUE_MODEL_PRESETS = {
    "colon": [
        "Immune_All_Low.pkl",
        "Pan_Fetal_Human.pkl",
        "Cells_Intestinal_Tract.pkl",
    ],
    "lung": ["Immune_All_Low.pkl", "Pan_Fetal_Human.pkl", "Human_Lung_Atlas.pkl"],
    "lung_airway": ["Immune_All_Low.pkl", "Pan_Fetal_Human.pkl", "Cells_Lung_Airway.pkl"],
    "liver": ["Immune_All_Low.pkl", "Pan_Fetal_Human.pkl", "Healthy_Human_Liver.pkl"],
    "heart": ["Immune_All_Low.pkl", "Pan_Fetal_Human.pkl", "Healthy_Adult_Heart.pkl"],
    "breast": ["Immune_All_Low.pkl", "Pan_Fetal_Human.pkl", "Cells_Adult_Breast.pkl"],
    "skin": ["Immune_All_Low.pkl", "Pan_Fetal_Human.pkl", "Adult_Human_Skin.pkl"],
    "pancreas": ["Immune_All_Low.pkl", "Pan_Fetal_Human.pkl", "Adult_Human_PancreaticIslet.pkl"],
    "brain": ["Immune_All_Low.pkl", "Pan_Fetal_Human.pkl", "Adult_Human_MTG.pkl"],
    "tonsil": ["Immune_All_Low.pkl", "Pan_Fetal_Human.pkl", "Cells_Human_Tonsil.pkl"],
    "unknown": ["Immune_All_Low.pkl", "Pan_Fetal_Human.pkl"],
}
```

### Native CellTypist Prediction

**IMPORTANT:** CellTypist's `celltypist.annotate()` function does NOT support internal batching.
It processes all cells at once. For large datasets, ensure sufficient memory is available.

When using native majority voting (`majority_voting=True`), CellTypist requires access to all cells
simultaneously to compute cluster-level statistics. Custom batching would break this functionality.

### Main Annotation Function

```python
def annotate_celltypist(
    adata: sc.AnnData,
    tissue: str = "unknown",
    ensemble_mode: bool = True,
    custom_model_path: Optional[str] = None,
    majority_voting: bool = True,
    over_clustering: Optional[str] = None,  # "leiden" or custom column
    min_prop: float = 0.0,
) -> sc.AnnData:
    """
    Annotate cells using CellTypist with native majority voting.

    Algorithm:
    1. Load tissue-specific model preset (2-3 models)
    2. Validate gene overlap for each model (skip if <25%)
    3. Re-normalize from raw for CellTypist compatibility
    4. Run prediction with native celltypist.annotate() including majority voting
    5. Ensemble: take highest confidence per cell across models

    Parameters:
    - majority_voting: Use CellTypist's native majority voting (recommended)
    - over_clustering: Column in adata.obs for cluster-based voting (e.g., "leiden")
    - min_prop: Minimum proportion for subcluster assignment (0.0 = no threshold)

    Output columns:
    - celltypist: Final cell type labels
    - celltypist_original: Per-cell predictions (before voting)
    - celltypist_confidence: Confidence scores (0-1)
    - celltypist_model: Which model contributed each cell
    """
    import celltypist
    from celltypist import models
    import gc

    # Determine models to run
    if custom_model_path:
        models_to_run = ["custom"]
    elif ensemble_mode:
        models_to_run = TISSUE_MODEL_PRESETS.get(tissue.lower(), TISSUE_MODEL_PRESETS["unknown"]).copy()
    else:
        models_to_run = ["Immune_All_Low.pkl"]

    # Load models and validate gene overlap
    loaded_models = {}
    all_overlap_genes = set()

    for model_name in models_to_run:
        if model_name == "custom":
            loaded_model = models.Model.load(custom_model_path)
        else:
            try:
                loaded_model = models.Model.load(model=model_name)
            except Exception:
                models.download_models(model=model_name)
                loaded_model = models.Model.load(model=model_name)

        # Calculate gene overlap
        model_genes = set(loaded_model.features)
        data_genes = set(adata.var_names)
        overlap_genes = model_genes & data_genes
        overlap_pct = 100 * len(overlap_genes) / len(model_genes)

        if overlap_pct < 25:
            print(f"Skipping {model_name}: only {overlap_pct:.1f}% gene overlap")
            continue

        loaded_models[model_name] = loaded_model
        all_overlap_genes.update(overlap_genes)

    if not loaded_models:
        raise ValueError("No models passed gene overlap threshold")

    # Prepare data for CellTypist (re-normalize from raw)
    original_X = adata.X.copy()
    original_raw = adata.raw.to_adata() if adata.raw is not None else None

    if adata.raw is not None:
        raw_adata = adata.raw.to_adata()
        sc.pp.normalize_total(raw_adata, target_sum=10000, exclude_highly_expressed=False)
        sc.pp.log1p(raw_adata)
        adata.X = raw_adata.X
        adata.raw = None  # Remove to pass CellTypist validation

    # Subset to overlapping genes
    genes_mask = adata.var_names.isin(all_overlap_genes)
    adata_subset = adata[:, genes_mask].copy()
    gc.collect()

    # Determine cluster column for native voting
    cluster_col = over_clustering or (
        "leiden" if "leiden" in adata.obs.columns else None
    )

    # Copy cluster info to subset if using native voting
    if majority_voting and cluster_col and cluster_col in adata.obs.columns:
        adata_subset.obs[cluster_col] = adata.obs[cluster_col].values

    # Run predictions for each model using native CellTypist
    all_model_predictions = {}

    for model_name, loaded_model in loaded_models.items():
        # Use native CellTypist annotate with majority voting
        prediction = celltypist.annotate(
            adata_subset,
            model=loaded_model,
            mode="best match",
            majority_voting=majority_voting,
            over_clustering=cluster_col if majority_voting else None,
            min_prop=min_prop,
        )

        # Get labels based on whether voting was enabled
        if majority_voting and "majority_voting" in prediction.predicted_labels.columns:
            labels = prediction.predicted_labels["majority_voting"]
        else:
            labels = prediction.predicted_labels["predicted_labels"]

        confidence = prediction.probability_matrix.max(axis=1).values
        all_model_predictions[model_name] = (labels, confidence)

    # Restore original data
    adata.X = original_X
    adata.raw = original_raw

    # Combine predictions (ensemble: highest confidence wins)
    if len(loaded_models) == 1:
        model_name = list(all_model_predictions.keys())[0]
        labels, confidence = all_model_predictions[model_name]
        per_cell_predictions = labels
        per_cell_confidence = confidence
        per_cell_source_model = pd.Series([model_name] * len(labels), index=labels.index)
    else:
        # Multi-model ensemble: pick highest confidence per cell
        cell_indices = list(all_model_predictions.values())[0][0].index
        final_labels = []
        final_confidence = []
        final_source_model = []

        for i, cell_idx in enumerate(cell_indices):
            best_conf = -1.0
            best_label = "Unknown"
            best_model = "none"

            for model_name, (labels, confidence) in all_model_predictions.items():
                cell_conf = confidence[i]
                if cell_conf > best_conf:
                    best_conf = cell_conf
                    best_label = labels.iloc[i]
                    best_model = model_name

            final_labels.append(best_label)
            final_confidence.append(best_conf)
            final_source_model.append(best_model)

        per_cell_predictions = pd.Series(final_labels, index=cell_indices)
        per_cell_confidence = np.array(final_confidence)
        per_cell_source_model = pd.Series(final_source_model, index=cell_indices)

    # Store results (native voting already applied if enabled)
    adata.obs["celltypist_original"] = per_cell_predictions.values
    adata.obs["celltypist_confidence"] = per_cell_confidence
    adata.obs["celltypist_model"] = per_cell_source_model.values
    adata.obs["celltypist"] = pd.Categorical(per_cell_predictions.values)

    return adata
```

---

## STEP 5: Convert to Ontology Codes

### Cell Type Patterns (Regex → Canonical Terms)

```python
CELL_TYPE_PATTERNS = {
    # Lymphoid lineage
    r"t\s*cells?,?\s*cd4\+?|cd4\+?\s*t": "cd4-positive, alpha-beta t cell",
    r"t\s*cells?,?\s*cd8\+?|cd8\+?\s*t": "cd8-positive, alpha-beta t cell",
    r"cd19.*cd20.*b\b|cd20.*cd19.*b\b": "b cell",
    r"cd19.*b\b": "b cell",
    r"cd20.*b\b": "b cell",
    r"t.*helper.*17|th17": "t-helper 17 cell",
    r"regulatory.*t|t.*regulatory|treg": "regulatory t cell",
    r"gamma.*delta.*t|gammadelta.*t|gdt": "gamma-delta t cell",
    r"mait": "mucosal invariant t cell",
    r"nkt|natural.*killer.*t|inkt": "invariant natural killer t-cell",
    r"^t\s*cell|^t\s+cells": "t cell",
    r"^b\s*cell|^b\s+cells?$": "b cell",
    r"\bnk\s*cell|\bnatural\s*killer": "natural killer cell",
    r"iga\+?\s*plasma": "iga plasma cell",
    r"igg\+?\s*plasma": "igg plasma cell",
    r"^plasma\s*$|plasmablast": "plasma cell",
    r"innate.*lymphoid|^ilc\b|^ilc\d": "innate lymphoid cell",

    # Myeloid lineage
    r"monocyte": "monocyte",
    r"macrophage|^m\d+": "macrophage",
    r"langerhans": "langerhans cell",
    r"^pdc\b|plasmacytoid\s*dc|plasmacytoid\s*dendritic": "plasmacytoid dendritic cell",
    r"dendritic\s*cells?|^dc\d?\b|^cdc\d?\b": "dendritic cell",
    r"mast\s*cell": "mast cell",
    r"neutrophils?": "neutrophil",
    r"basophils?": "basophil",
    r"eosinophils?": "eosinophil",

    # Stromal
    r"myofibroblast": "myofibroblast cell",
    r"fibroblast|^s\d+|reticular": "fibroblast",
    r"smooth\s*muscle": "smooth muscle cell",
    r"pericyte": "pericyte",
    r"stromal": "stromal cell",

    # Endothelial
    r"endotheli|^ve\b|^ec\b|ecs\b": "endothelial cell",
    r"lymphatic.*ec|lymphatic.*endothel": "lymphatic endothelial cell",

    # Epithelial
    r"epitheli": "epithelial cell",
    r"keratinocytes?": "keratinocyte",
    r"basal": "basal cell",
    r"goblet": "goblet cell",
    r"enterocyte|entero": "enterocyte",
    r"colonocyte|colono": "colon glandular cell",
    r"paneth": "paneth cell",
    r"tuft": "intestinal tuft cell",
    r"ciliated": "ciliated epithelial cell",

    # Neural
    r"astrocytes?": "astrocyte",
    r"oligodendrocyte": "oligodendrocyte",
    r"glia": "glial cell",
    r"neuron|neural\s*cell": "neuron",

    # Stem/Progenitor
    r"^hsc\b|hematopoietic.*stem": "hematopoietic stem cell",
    r"^msc\b|mesenchymal.*stem": "mesenchymal stem cell",
    r"stem": "stem cell",
    r"progenitor|precursor": "progenitor cell",
    r"transit.*amplifying|^ta\b": "transit amplifying cell",

    # Other
    r"hepatocyte": "hepatocyte",
    r"adipocyte|adipose": "adipocyte",
    r"platelets?|thrombocyte": "platelet",
    r"^rbc\b|red\s*blood\s*cell": "erythrocyte",
    # ... ~100+ more patterns in full implementation
}
```

### Token Extraction

```python
def extract_biological_tokens(label: str) -> Dict[str, List[str]]:
    """
    Extract key biological identifiers for token-based matching.

    Returns:
        {
            'markers': ['cd4', 'cd8', 'cd19', ...],  # CD markers
            'proteins': ['igg', 'iga', 'spp1', ...],  # Immunoglobulins, genes
            'core_words': ['b', 't', 'plasma', 'endothelial', ...],
            'modifiers': ['positive', 'mature', 'activated', ...]
        }
    """
    import re

    label_lower = label.lower().strip()
    tokens = {"markers": [], "proteins": [], "core_words": [], "modifiers": []}

    # Extract CD markers: CD4, CD8, CD19, etc.
    cd_markers = re.findall(r"cd\d+", label_lower)
    tokens["markers"].extend(cd_markers)

    # Extract immunoglobulin types: IgG, IgA, IgM
    ig_types = re.findall(r"ig[gamedGAMED]", label_lower)
    tokens["proteins"].extend([ig.lower() for ig in ig_types])

    # Extract gene names (uppercase + plus sign)
    gene_markers = re.findall(r"\b[A-Z0-9]{3,}\+", label)
    tokens["proteins"].extend([g.replace("+", "").lower() for g in gene_markers])

    # Clean and extract core words
    cleaned = re.sub(r"cd\d+", "", label_lower)
    cleaned = re.sub(r"ig[gamed]", "", cleaned)
    cleaned = re.sub(r"[+\-]", " ", cleaned)
    cleaned = re.sub(r"\d+", "", cleaned)
    cleaned = re.sub(r"\s+", " ", cleaned).strip()

    MODIFIERS = {"positive", "negative", "like", "type", "mature", "immature",
                 "activated", "resting", "proliferating", "pro", "inflammatory"}
    GENERIC = {"cell", "cells"}
    MEANINGFUL_SHORT = {"b", "t", "nk", "dc", "ec", "ve", "ta", "m1", "m2"}

    for word in cleaned.split():
        if word in MODIFIERS:
            tokens["modifiers"].append(word)
        elif word in GENERIC:
            pass
        elif word in MEANINGFUL_SHORT:
            tokens["core_words"].append(word)
        elif len(word) > 1:
            tokens["core_words"].append(word)

    return tokens
```

### Multi-Tier Search Function

```python
def search_ontology_index(
    labels: List[str],
    ontology_index: Dict[str, Dict[str, Dict[str, str]]],
    annotation_type: str = "cell_type",
    min_score: float = 0.7,
) -> Dict[str, List[Dict]]:
    """
    Fast ontology search with sequential tier matching.

    Tiers:
    - Tier 0: Pattern-based canonicalization (score 0.95)
    - Tier 1: Exact match (score 1.0)
    - Tier 2: Token-based (score 0.70-0.85)
    - Tier 3: Word overlap (score 0.5-0.7)

    Args:
        labels: Cell type labels to search
        ontology_index: Pre-loaded index {onto: {label_lower: {id, name}}}
        annotation_type: "cell_type" (CL), "pathology" (NCIT), "anatomy" (UBERON)
        min_score: Minimum match score to accept

    Returns:
        {label: [{id, name, ontology, score, match_type}, ...]}
    """
    import re

    # Select ontologies based on type
    if annotation_type == "cell_type":
        ontologies = ["cl"]
    elif annotation_type == "pathology":
        ontologies = ["ncit", "cl"]
    elif annotation_type == "anatomy":
        ontologies = ["uberon", "cl"]
    else:
        ontologies = ["cl", "ncit", "uberon"]

    results = {label: [] for label in labels}

    for label in labels:
        label_lower = label.lower().strip().replace("_", " ")
        label_normalized = label.replace("_", " ")

        # Extract tokens for Tier 2
        tokens = extract_biological_tokens(label_normalized)

        # TIER 0: Pattern-based canonicalization
        canonical_term = None
        for pattern, search_term in CELL_TYPE_PATTERNS.items():
            if re.search(pattern, label_lower):
                canonical_term = search_term
                break

        search_label = canonical_term if canonical_term else label_lower
        is_pattern_match = canonical_term is not None

        if is_pattern_match:
            tokens = extract_biological_tokens(canonical_term)

        # Search across ontologies
        tier_results = []

        for onto_prefix in ontologies:
            onto_dict = ontology_index.get(onto_prefix.lower(), {})
            ontology_matches = []

            # TIER 1: Exact lookup
            if search_label in onto_dict:
                term = onto_dict[search_label]
                ontology_matches.append({
                    "id": term["id"],
                    "name": term["name"],
                    "ontology": onto_prefix,
                    "score": 0.95 if is_pattern_match else 1.0,
                    "match_type": "tier0_pattern" if is_pattern_match else "tier1_exact",
                })
            else:
                # TIER 2-3: Fuzzy matching
                for term_label_lower, term_data in onto_dict.items():
                    # Skip imported terms
                    term_id_prefix = term_data["id"].split(":")[0].upper()
                    if term_id_prefix != onto_prefix.upper():
                        continue

                    # Skip obsolete
                    if "obsolete" in term_data["name"].lower():
                        continue

                    score = score_match(search_label, term_label_lower, tokens, is_pattern_match)

                    if score >= min_score:
                        ontology_matches.append({
                            "id": term_data["id"],
                            "name": term_data["name"],
                            "ontology": onto_prefix,
                            "score": score,
                            "match_type": "tier2_token" if score >= 0.7 else "tier3_overlap",
                        })

            if ontology_matches:
                tier_results.extend(ontology_matches)
                if onto_prefix == "cl" and any(m["score"] >= 0.8 for m in ontology_matches):
                    break  # Good CL match found

        # Sort by score, deduplicate
        seen_ids = set()
        unique_results = []
        for result in sorted(tier_results, key=lambda x: x["score"], reverse=True):
            if result["id"] not in seen_ids:
                seen_ids.add(result["id"])
                unique_results.append(result)

        results[label] = unique_results

    return results
```

---

## Configuration & Defaults

```python
# Training defaults
DEFAULT_MAX_CELLS_PER_REF = 500000
DEFAULT_TARGET_SUM = 10000  # CellTypist normalization target
MIN_CELLS_PER_TYPE = 50    # Minimum cells for training

# Annotation defaults
DEFAULT_CONFIDENCE_THRESHOLD = 0.3
DEFAULT_BATCH_SIZE = 5000
MIN_GENE_OVERLAP_PCT = 25  # Skip model if below

# Ontology matching defaults
MIN_MATCH_SCORE = 0.7
```

---

## Dependencies

```
# Core
scanpy>=1.9.0
anndata>=0.8.0
pandas>=1.3.0
numpy>=1.20.0
scipy>=1.7.0

# CellTypist
celltypist>=1.5.0
# OR celltypist-SO for GPL-free (no majority voting)

# CellxGene download
cellxgene-census>=1.0.0

# Ontology (optional, for hierarchy operations)
owlready2>=0.40

# Memory monitoring
psutil>=5.8.0
```

---

## Critical Findings for Spatial Data

The following findings emerged from investigation of CellTypist performance on Xenium spatial transcriptomics data.

### Training Custom Models (Findings 1-4)

### 1. Custom Models Dramatically Improve Gene Overlap

**Problem:** Pre-trained CellTypist models have poor gene overlap with targeted spatial panels:
- Human_IPF_Lung.pkl: Only 7.8% overlap with 430-gene Xenium panel (34/430 genes)
- CellTypist requires minimum 25% overlap to produce reliable predictions
- Low overlap causes annotation failures or highly uncertain predictions

**Solution:** Train custom models on CellxGene reference data subsetted to panel genes:

| Metric | Pre-trained Model | Custom Panel Model |
|--------|-------------------|-------------------|
| Gene overlap | 7.8% (34/430) | **100%** (430/430) |
| Usable features | 34 genes | 430 genes |
| Training data | Unknown | 200k cells from HLCA + COVID lung |
| Cell types | 44 | 44 |

**Result:** 100% gene utilization enables CellTypist to leverage the full panel for classification.

### 2. Ensembl-HUGO Gene Mapping is Essential

**Problem:** CellxGene datasets use inconsistent gene identifiers:
- HLCA: Ensembl IDs in var_names (ENSG...) with HUGO symbols in `feature_name`
- COVID lung: Mixed Ensembl and HUGO in `feature_name` column
- Xenium panels: Use HUGO gene symbols

Without mapping, gene name mismatches cause 0% overlap between reference and query data.

**Solution:** Two-step gene name normalization:
1. Use `feature_name` column as primary gene names
2. Apply Ensembl→HUGO mapping for any remaining Ensembl IDs

```python
# Load mapping once (46,920 mappings from Ensembl BioMart)
ensembl_to_hugo = load_ensembl_to_hugo_mapping()

# Apply to each reference dataset
adata = normalize_gene_names(adata, ensembl_to_hugo)
```

**Key insight:** HUGO symbols pass through unchanged, so it's safe to ALWAYS apply the mapping.

### 3. Memory-Efficient Backed Mode for Large References

**Problem:** Large CellxGene datasets exceed available memory:
- HLCA: 584,884 cells × 27,402 genes = 5.6GB compressed, ~60GB dense
- Loading fully causes OOM on 32GB machines

**Solution:** Use scanpy's backed mode for memory-efficient loading:

```python
# Open without loading expression matrix
adata_backed = sc.read_h5ad(path, backed='r')

# Determine sample indices from metadata only
indices = _stratified_sample_indices(adata_backed.obs, label_column, max_cells)

# Load ONLY the selected cells
adata = adata_backed[indices].to_memory()
```

**Memory savings:** For HLCA, load 100k of 584k cells = ~17% of data loaded.

### 4. Panel Gene Subsetting Enables Full Gene Utilization

**Problem:** Pre-trained models contain thousands of genes, but spatial panels are targeted (300-500 genes). When subsetting a pre-trained model to panel genes, most discriminative features are lost.

**Solution:** Train custom models on reference data already subsetted to panel genes:

```python
# Get panel genes from Xenium workflow
panel_genes = list(xenium_adata.var_names)  # 430 genes

# Combine references with panel subsetting
combined = combine_references(
    reference_paths=["hlca.h5ad", "covid_lung.h5ad"],
    label_columns=["ann_finest_level", "cell_type"],
    target_genes=panel_genes,  # Subset BEFORE training
)

# Train on panel-subsetted data
model = celltypist.train(combined, labels="unified_cell_type", ...)
```

**Workflow output:**
```
[1/2] Loading: human_lung_cell_atlas.h5ad
  Applied gene mapping: 228 converted, 22,279 already HUGO
  Gene names now: TSPAN6, TNMD, ...

Subsetting to 430 target genes...
  Reference 0: 430 genes after subset
  Reference 1: 430 genes after subset

Training CellTypist model...
  🧬 430 features are selected
  ✅ Model training done!
```

### Annotation Considerations (Findings 5-8)

### 5. Majority Voting Collapses Cell Types (Spatial-Specific)

**Problem:** CellTypist's native majority voting assigns the dominant cell type within each cluster to all cells in that cluster. With coarse clustering (e.g., 9 leiden clusters), this can dramatically reduce cell type diversity.

**Example:** A 122k cell dataset with 13 CellTypist-predicted cell types was collapsed to only 2 types when majority voting was enabled:
- 8 of 9 clusters had Fibroblast as dominant type → all cells labeled "Fibroblast"
- 1 cluster had Macrophage as dominant → all cells labeled "Macrophage"

**Recommendation:** Set `majority_voting=False` by default for spatial data. Only enable when:
- Fine-grained over-clustering is available (many small clusters)
- Spatial consistency is more important than cell type diversity

### 6. Re-normalization from Raw is Essential

**Problem:** CellTypist expects log1p(10k CPM) normalization where cells sum to exactly 10,000. However, scanpy's default `normalize_total()` uses `exclude_highly_expressed=True`, which produces variable cell sums.

**Solution:** Always re-normalize from raw counts before CellTypist prediction:
```python
if adata.raw is not None:
    raw_adata = adata.raw.to_adata()
    sc.pp.normalize_total(raw_adata, target_sum=10000, exclude_highly_expressed=False)
    sc.pp.log1p(raw_adata)
    adata.X = raw_adata.X
```

**Key parameter:** `exclude_highly_expressed=False` ensures all cells sum to exactly 10k.

### 7. Domain Shift Between scRNAseq and Spatial

**Inherent challenge:** Spatial transcriptomics (e.g., Xenium) captures fewer transcripts per cell than scRNAseq, resulting in sparser expression matrices:

| Metric | scRNAseq (Study 619) | Spatial (Xenium) |
|--------|---------------------|------------------|
| Non-zero gene fraction | ~10.6% | ~1.6% |

**Impact:** Lower confidence scores on spatial data (mean ~0.24 vs. ~0.70 for scRNAseq). This is expected behavior, not a pipeline bug.

**Recommendation:** Use confidence thresholds appropriate for spatial data (e.g., 0.1-0.2 instead of 0.3).

### 8. Pipeline Workflow for Annotation

The recommended two-step workflow:

1. **CellTypist annotation** → Raw predictions with confidence scores
   ```python
   adata = annotate_celltypist(adata, tissue="synovium", majority_voting=False)
   # Creates: celltypist, celltypist_confidence, celltypist_model columns
   ```

2. **Ontology standardization** → Map to Cell Ontology IDs and canonical names
   ```python
   adata, mappings = add_ontology_ids(
       adata,
       source_col="celltypist",
       target_col="celltypist_ontology_id",
       name_col="celltypist_ontology_name"
   )
   # Creates: celltypist_ontology_id (CL:0000xxx), celltypist_ontology_name columns
   ```

---

## Usage Example

```python
# 1. Download reference
download_cellxgene("colon_immune_niches", output_dir="./references")

# 2. Combine references
combined = combine_references(
    reference_paths=["./references/colon_immune_niches.h5ad",
                     "./references/colon_ulcerative_colitis.h5ad"],
    label_columns=["cell_type", "cell_type"],
    max_cells_per_ref=100000,
)

# 3. Train custom model (optional - for panel-specific)
# First subset to panel genes
panel_genes = list(xenium_adata.var_names)
combined_subset = combined[:, combined.var_names.isin(panel_genes)].copy()

train_celltypist_model(
    combined_subset,
    label_column="unified_cell_type",
    output_path="./custom_colon_model.pkl",
)

# 4. Annotate query data
query_adata = sc.read_h5ad("./query_data.h5ad")
annotated = annotate_celltypist(
    query_adata,
    tissue="colon",
    ensemble_mode=True,
    confidence_threshold=0.3,
)

# 5. Convert to ontology codes
unique_labels = annotated.obs["cell_type"].unique().tolist()
ontology_index = load_ontology_index()
mappings = search_ontology_index(unique_labels, ontology_index)

# Apply mappings
annotated.obs["cell_type_ontology_id"] = annotated.obs["cell_type"].map(
    {label: mappings[label][0]["id"] if mappings[label] else "Unknown"
     for label in unique_labels}
)
```
