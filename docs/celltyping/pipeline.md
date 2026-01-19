# Pipeline & API Reference

**Complete technical specification for the SpatialCore cell typing workflow.**

---

## Overview

The pipeline trains custom CellTypist models on scRNA-seq references, produces calibrated confidence scores via z-score transformation, and maps predictions to Cell Ontology (CL) IDs.

![SpatialCore Cell Typing Pipeline](images/spatialcore_celltypist_pipeline.png)

---

## Column Naming Convention (CellxGene Standard)

All outputs use the [CellxGene schema](https://github.com/chanzuckerberg/single-cell-curation/blob/main/schema/3.0.0/schema.md):

| Column | Type | Description |
|--------|------|-------------|
| `cell_type` | `str` | Final predicted cell type |
| `cell_type_confidence` | `float` | Z-score transformed confidence [0, 1] |
| `cell_type_confidence_raw` | `float` | Raw CellTypist decision score |
| `cell_type_ontology_term_id` | `str` | Cell Ontology ID (e.g., `CL:0000624`) |
| `cell_type_ontology_label` | `str` | Canonical ontology name |
| `original_label` | `str` | Raw label from source reference |
| `reference_source` | `str` | Which reference file the cell came from |

---

## Phase 1: Data Acquisition

### `acquire_reference()`

Download reference data from a source and store to a destination.

```python
def acquire_reference(
    source: str,
    output: Union[str, Path],
    force: bool = False,
    **kwargs,
) -> str:
```

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `source` | `str` | Source URI (see supported schemes below) |
| `output` | `str` or `Path` | Destination path or URI |
| `force` | `bool` | Re-download even if output exists |
| `**kwargs` | | Source-specific options (e.g., `max_cells`, `auth_token`) |

**Supported Source Schemes:**

| Scheme | Format | Example |
|--------|--------|---------|
| CellxGene dataset | `cellxgene://dataset_key` | `cellxgene://human_lung_cell_atlas` |
| CellxGene query | `cellxgene://?tissue=X&disease=Y` | `cellxgene://?tissue=liver&disease=normal` |
| Synapse | `synapse://synXXXXXXXX` | `synapse://syn12345678` |

**Supported Destination Schemes:**

| Scheme | Format | Auth Required |
|--------|--------|---------------|
| Local | `/path/to/file.h5ad` | No |
| GCS | `gs://bucket/path/file.h5ad` | `GOOGLE_APPLICATION_CREDENTIALS` |
| S3 | `s3://bucket/path/file.h5ad` | `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY` |

**Example:**

```python
from spatialcore.annotation import acquire_reference

# Download from CellxGene → store locally
path = acquire_reference(
    source="cellxgene://human_lung_cell_atlas",
    output="/data/references/hlca.h5ad",
)

# CellxGene query with filters → store to GCS
gcs_path = acquire_reference(
    source="cellxgene://?tissue=liver&disease=normal",
    output="gs://my-bucket/references/healthy_liver.h5ad",
    max_cells=100000,
)
```

### Available CellxGene Datasets

```python
from spatialcore.annotation import list_available_datasets

datasets = list_available_datasets()
print(datasets)
```

| Dataset Key | Tissue | Description |
|-------------|--------|-------------|
| `healthy_human_liver` | liver | Healthy human liver scRNA-seq |
| `colon_immune_niches` | colon | Colon immune microenvironment |
| `human_lung_cell_atlas` | lung | HLCA reference atlas |
| `lung_covid` | lung | COVID-19 lung atlas |

---

## Phase 2: Training & Annotation

### High-Level API: `train_and_annotate()`

Complete pipeline in a single call. Recommended for most users.

```python
def train_and_annotate(
    adata: AnnData,
    references: List[Union[str, Path]],
    tissue: str = "unknown",
    label_columns: Optional[List[str]] = None,
    balance_strategy: Literal["proportional", "equal"] = "proportional",
    max_cells_per_type: int = 10000,
    max_cells_per_ref: int = 100000,
    confidence_threshold: float = 0.8,
    model_output: Optional[Union[str, Path]] = None,
    plot_output: Optional[Union[str, Path]] = None,
    add_ontology: bool = True,
    generate_plots: bool = True,
    copy: bool = False,
) -> AnnData:
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `adata` | `AnnData` | required | Spatial data to annotate |
| `references` | `List[str]` | required | Paths/URIs to reference files |
| `tissue` | `str` | `"unknown"` | Tissue type for model naming |
| `label_columns` | `List[str]` | `None` | Cell type column per reference (auto-detect if None) |
| `balance_strategy` | `str` | `"proportional"` | Source balancing strategy |
| `max_cells_per_type` | `int` | `10000` | Max cells per type after balancing |
| `max_cells_per_ref` | `int` | `100000` | Max cells to load per reference |
| `confidence_threshold` | `float` | `0.8` | Below this → "Unassigned" |
| `model_output` | `Path` | `None` | Save model to this path |
| `plot_output` | `Path` | `None` | Save plots to this directory |
| `add_ontology` | `bool` | `True` | Map predictions to CL IDs |
| `generate_plots` | `bool` | `True` | Generate validation plots |

**Example:**

```python
from spatialcore.annotation import train_and_annotate
import scanpy as sc

adata = sc.read_h5ad("xenium_lung.h5ad")

adata = train_and_annotate(
    adata,
    references=[
        "gs://my-bucket/references/hlca.h5ad",
        "/local/data/lung.h5ad",
    ],
    tissue="lung",
    balance_strategy="proportional",
    confidence_threshold=0.8,
    model_output="./models/lung_custom_v1.pkl",
    plot_output="./qc_plots/",
)

print(adata.obs["cell_type"].value_counts())
```

### Config-Driven API: `TrainingConfig`

For reproducible workflows, use YAML configuration.

```yaml
# training_config.yaml
tissue: lung
references:
  - gs://my-bucket/references/hlca.h5ad
  - /local/data/lung.h5ad
balance_strategy: proportional
max_cells_per_type: 10000
max_cells_per_ref: 100000
confidence_threshold: 0.8
add_ontology: true
generate_plots: true
```

```python
from spatialcore.annotation import TrainingConfig, train_and_annotate_config

config = TrainingConfig.from_yaml("training_config.yaml")
adata = train_and_annotate_config(adata, config, plot_output="./qc/")
```

---

## Low-Level Functions

For users who need fine-grained control over each stage.

### `combine_references()`

Combine multiple reference datasets with memory-efficient loading.

```python
def combine_references(
    reference_paths: List[Union[str, Path]],
    label_columns: List[str],
    output_column: str = "original_label",
    max_cells_per_ref: int = 100000,
    target_genes: Optional[List[str]] = None,
    normalize_data: bool = True,
    min_cells_per_type: int = 10,
    exclude_labels: Optional[List[str]] = None,
    filter_min_cells: bool = True,
) -> AnnData:
```

**Parameters:**

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `reference_paths` | `List[str]` | required | Paths/URIs to reference h5ad files |
| `label_columns` | `List[str]` | required | Cell type column for each reference |
| `output_column` | `str` | `"original_label"` | Column name for unified labels |
| `max_cells_per_ref` | `int` | `100000` | Max cells to load per reference |
| `target_genes` | `List[str]` | `None` | Panel genes to subset to |
| `normalize_data` | `bool` | `True` | Apply log1p(10k) normalization |
| `min_cells_per_type` | `int` | `10` | Minimum cells per type |
| `exclude_labels` | `List[str]` | `None` | Labels to exclude |
| `filter_min_cells` | `bool` | `True` | Remove types below min_cells_per_type |

**Label Filtering (`exclude_labels`):**

By default, ambiguous labels are removed:

```python
DEFAULT_EXCLUDE_LABELS = [
    "unknown", "Unknown", "UNKNOWN",
    "unassigned", "Unassigned",
    "na", "NA", "N/A", "n/a",
    "none", "None", "null",
    "doublet", "Doublet",
    "low quality", "Low quality",
]
```

**Supported reference_paths:**
- Local: `/data/references/lung.h5ad`
- GCS: `gs://bucket/references/lung.h5ad`
- S3: `s3://bucket/references/lung.h5ad`

### `subsample_balanced()`

Source-aware balanced subsampling with semantic grouping.

```python
def subsample_balanced(
    adata: AnnData,
    label_column: str,
    max_cells_per_type: int = 10000,
    min_cells_per_type: int = 50,
    source_column: Optional[str] = "reference_source",
    source_balance: Literal["proportional", "equal"] = "proportional",
    group_by_column: Optional[str] = None,
    target_proportions: Optional[Union[Dict[str, float], str, Path]] = None,
    random_state: int = 42,
) -> AnnData:
```

**The `group_by_column` parameter:**

Different references may use different names for the same cell type:

| Reference A | Reference B | CL ID |
|-------------|-------------|-------|
| "CD4-positive, alpha-beta T cell" | "CD4+ T cells" | CL:0000624 |

By setting `group_by_column="cell_type_ontology_term_id"`, cells are grouped by CL ID:

```python
balanced = subsample_balanced(
    combined,
    label_column="original_label",
    group_by_column="cell_type_ontology_term_id",
    source_balance="proportional",
)
```

**"Cap & Fill" Algorithm:**

```
FOR each cell_type (or CL ID group):
  1. IDENTIFY SOURCES that have this type
  2. CALCULATE per-source targets:
     IF source_balance == "proportional":
       target[src] = total × (src_count / total_count)
     ELSE (equal):
       target[src] = total / n_sources
  3. ENFORCE minimums and maximums
  4. FILL SHORTFALL from sources with unused capacity
  5. SAMPLE from each source
```

**Example:**

```
Macrophage: 35K total (Study1: 30K, Study2: 5K)
Target: 10K cells

PROPORTIONAL BALANCE:
  Study1: 10K × (30K/35K) = 8,571 cells
  Study2: 10K × (5K/35K)  = 1,429 cells

EQUAL BALANCE:
  Study1: 5,000 cells
  Study2: 5,000 cells
```

**The `target_proportions` parameter:**

For FACS-sorted or enriched references where a cell type only exists in the enriched source:

| Data Source | NK Cells | Problem |
|-------------|----------|---------|
| Tissue atlas (500K cells) | 0 | NK rare/absent |
| FACS pure NK (5K cells) | 5,000 (100%) | Artificially enriched |
| Naive combination | 5,000 | Model thinks NK = 10% |
| Biological reality | ~0.25% | NK should be rare |

```python
balanced = subsample_balanced(
    combined,
    label_column="original_label",
    max_cells_per_type=10000,
    target_proportions={"NK cell": 0.0025},  # 0.25%
)
```

**Supported formats for target_proportions:**

```python
# Dict (inline)
target_proportions={"NK cell": 0.0025, "plasma cell": 0.001}

# JSON file
target_proportions="proportions.json"

# CSV file
target_proportions="proportions.csv"
```

### `train_celltypist_model()`

Train a custom CellTypist logistic regression model.

```python
def train_celltypist_model(
    adata: AnnData,
    label_column: str = "unified_cell_type",
    model_name: str = "custom_model",
    output_path: Optional[Union[str, Path]] = None,
    use_SGD: bool = True,
    mini_batch: bool = True,
    balance_cell_type: bool = True,
) -> Dict[str, Any]:
```

**Artifacts saved:**

| File | Description |
|------|-------------|
| `{name}.pkl` | CellTypist model weights |
| `{name}_celltypist.json` | Training metadata |
| `{name}_colors.json` | Color palette for visualization |

### `annotate_celltypist()`

Apply model to spatial data with z-score confidence transformation.

```python
def annotate_celltypist(
    adata: AnnData,
    custom_model_path: Optional[Union[str, Path]] = None,
    majority_voting: bool = False,
    confidence_transform: Optional[str] = "zscore",
) -> AnnData:
```

**IMPORTANT: `majority_voting=False` for spatial data**

```
SCRNASEQ (voting OK):
• Clusters are fine-grained
• Voting improves consistency

SPATIAL (voting DANGEROUS):
• Spatial clustering may be coarse
• A single "immune" cluster might contain:
  - 1000 macrophages
  - 50 T cells
  - 30 B cells
• Voting assigns dominant type to ALL cells
• Result: All 1080 become macrophages (WRONG!)

SOLUTION: Use majority_voting=False
```

### Complete Low-Level Example

```python
from spatialcore.annotation import (
    get_panel_genes,
    combine_references,
    add_ontology_ids,
    subsample_balanced,
    train_celltypist_model,
    annotate_celltypist,
    filter_low_confidence,
)
import scanpy as sc

# STAGE 1: Load spatial data and extract panel genes
xenium = sc.read_h5ad("xenium_lung.h5ad")
panel_genes = get_panel_genes(xenium)

# STAGE 2: Combine references
combined = combine_references(
    reference_paths=[
        "gs://my-bucket/references/hlca.h5ad",
        "/local/data/inhouse_lung.h5ad",
    ],
    label_columns=["cell_type", "cell_type"],
    output_column="original_label",
    max_cells_per_ref=100000,
    target_genes=panel_genes,
)

# STAGE 3: Fill missing ontology IDs
combined, _, _ = add_ontology_ids(
    combined,
    source_col="original_label",
    target_col="cell_type_ontology_term_id",
    skip_if_exists=True,
)

# STAGE 4: Source-aware balanced subsampling
balanced = subsample_balanced(
    combined,
    label_column="original_label",
    group_by_column="cell_type_ontology_term_id",
    source_column="reference_source",
    source_balance="proportional",
    max_cells_per_type=10000,
)

# STAGE 5: Train CellTypist model
result = train_celltypist_model(
    balanced,
    label_column="cell_type_ontology_term_id",
    output_path="./models/lung_custom_v1.pkl",
)

# STAGE 6: Annotate spatial data
xenium = annotate_celltypist(
    xenium,
    custom_model_path="./models/lung_custom_v1.pkl",
    majority_voting=False,
    confidence_transform="zscore",
)

# STAGE 7: Apply confidence threshold
xenium = filter_low_confidence(
    xenium,
    label_column="cell_type",
    confidence_column="cell_type_confidence",
    threshold=0.8,
)

# STAGE 8: Add ontology IDs to predictions
xenium, _, _ = add_ontology_ids(
    xenium,
    source_col="cell_type",
    target_col="cell_type_ontology_term_id",
)
```

---

## Phase 3: Plotting & Validation

### `generate_annotation_plots()`

Generate all validation plots in one call.

```python
def generate_annotation_plots(
    adata: AnnData,
    label_column: str = "cell_type",
    confidence_column: str = "cell_type_confidence",
    output_dir: Optional[Union[str, Path]] = None,
    prefix: str = "celltyping",
    confidence_threshold: float = 0.8,
    markers: Optional[Dict[str, List[str]]] = None,
    n_deg_genes: int = 10,
) -> Dict:
```

**Output files:**

| Plot | Filename | Description |
|------|----------|-------------|
| DEG Heatmap | `{prefix}_deg_heatmap.png` | Top N DEGs per cell type |
| 2D Validation | `{prefix}_2d_validation.png` | Confidence vs marker (GMM-3) |
| Confidence | `{prefix}_confidence.png` | Spatial + jitter with threshold |
| Ontology Table | `{prefix}_ontology_mapping.png` | Mapping statistics |

---

## Expression Normalization Detection

CellTypist requires log1p(10k) normalized data. SpatialCore uses robust detection:

```
PHASE 1: Search for Raw Counts
──────────────────────────────────────────────────────────────
Check in priority order:
  1. adata.layers["counts"]
  2. adata.layers["raw_counts"]
  3. adata.layers["raw"]
  4. adata.raw.X
  5. adata.X

Integer Test (with floating-point tolerance):
  • Sample 10,000 non-zero values
  • Check: |value - round(value)| < 1e-6
  • Pass: >95% of values are integer-like

If raw counts found → SAFE PATH: normalize from raw

PHASE 2: Verify X is log1p(10k) via expm1 Reversal
──────────────────────────────────────────────────────────────
Key insight: If X = log1p(counts / total * target_sum), then:
             expm1(X).sum(axis=1) ≈ target_sum

Verification:
  • Reverse log1p: reversed = expm1(X_sample)
  • Compute row sums: row_sums = reversed.sum(axis=1)
  • Check median: 8,000 < median(row_sums) < 12,000 → log1p_10k
                 800,000 < median < 1,200,000 → log1p_cpm

If verified log1p_10k → SAFE PATH: use as-is
If log1p_cpm or unknown → ERROR (unless unsafe_force=True)
```

**Decision Matrix:**

| `raw_source` | `x_state` | Action |
|--------------|-----------|--------|
| Found (any location) | Any | Normalize from raw → **SAFE** |
| None | `log1p_10k` (verified) | Use X as-is → **SAFE** |
| None | `raw` | Normalize X directly → **SAFE** |
| None | `log1p_cpm` | **ERROR** (wrong scale) |
| None | `unknown` | **ERROR** |

---

## Ontology Mapping

### 4-Tier Matching System

![4-Tier Ontology Matching](images/matching_tiers.png)

```
INPUT: "CD4+ T cells"
           │
           ▼
┌──────────────────────────────────────────────────────────┐
│ TIER 0: Pattern Canonicalization                          │
│ "CD4+ T cells" → "cd4-positive, alpha-beta t cell"       │
│ Score: 0.95                                               │
└──────────────────────────────────────────────────────────┘
           │
           ▼
┌──────────────────────────────────────────────────────────┐
│ TIER 1: Exact Match                                       │
│ "cd4-positive, alpha-beta t cell" → CL:0000624           │
│ Score: 1.0 (exact) or 0.95 (via pattern)                 │
└──────────────────────────────────────────────────────────┘
           │
           ▼ (if no exact match)
┌──────────────────────────────────────────────────────────┐
│ TIER 2: Token-Based Match                                 │
│ Extracts biological tokens (CD markers, core words)       │
│ Score: 0.60-0.85                                          │
└──────────────────────────────────────────────────────────┘
           │
           ▼ (if no token match)
┌──────────────────────────────────────────────────────────┐
│ TIER 3: Word Overlap (Jaccard Similarity)                 │
│ Score: 0.50-0.70                                          │
└──────────────────────────────────────────────────────────┘
           │
           ▼
OUTPUT: CL:0000624, "cd4-positive, alpha-beta t cell"
```

### Files

| File | Location | Purpose |
|------|----------|---------|
| `ontology_index.json` | `data/ontology_mappings/` | Pre-built CL term lookup |
| `patterns.py` | `annotation/` | Regex patterns for Tier 0 |
| `canonical_markers.json` | `data/markers/` | Marker genes for 50+ cell types |

### Customizing Pattern Matching

Add patterns in `src/spatialcore/annotation/patterns.py`:

```python
CELL_TYPE_PATTERNS = {
    r"regex_pattern": "canonical cl term name",
    ...
}
```

**Common regex patterns:**

| Pattern | Meaning | Example |
|---------|---------|---------|
| `\b` | Word boundary | `\bnk\b` matches "NK" but not "unk" |
| `.*` | Any characters | `cd4.*t` matches "CD4+ T cell" |
| `?` | Optional character | `cells?` matches "cell" or "cells" |
| `\|` | OR | `nk\|natural killer` matches either |

**Testing your pattern:**

```python
from spatialcore.annotation.patterns import get_canonical_term

print(get_canonical_term("Club (nasal)"))  # → "club cell"
print(get_canonical_term("Migratory DCs")) # → "migratory dendritic cell"
```

---

## API Summary

### Phase 1: Data Acquisition

| Function | Purpose |
|----------|---------|
| `acquire_reference()` | Download from CellxGene/Synapse → store to local/cloud |
| `list_available_datasets()` | List predefined CellxGene datasets |

### Phase 2: Training & Annotation

| Function | Purpose |
|----------|---------|
| `train_and_annotate()` | **Complete pipeline in one call** |
| `train_and_annotate_config()` | Config-driven version |
| `TrainingConfig` | YAML-serializable configuration |
| `get_panel_genes()` | Extract gene list from spatial data |
| `combine_references()` | Load + normalize + filter + concatenate |
| `add_ontology_ids()` | Map labels to CL IDs |
| `subsample_balanced()` | Source-aware balanced subsampling |
| `train_celltypist_model()` | Train custom CellTypist model |
| `annotate_celltypist()` | Apply model with z-score confidence |
| `filter_low_confidence()` | Mark low-confidence as Unassigned |

### Phase 3: Plotting

| Function | Purpose |
|----------|---------|
| `generate_annotation_plots()` | **All validation plots in one call** |
| `plot_deg_heatmap()` | DEG heatmap with top markers |
| `plot_2d_validation()` | GMM-3 marker validation |
| `plot_celltype_confidence()` | Spatial + jitter confidence |
| `plot_ontology_mapping()` | Ontology mapping table |

---

## Error Handling

| Error | Cause | Solution |
|-------|-------|----------|
| `Cannot safely normalize data` | No raw counts and X not verified | Provide raw counts or use `unsafe_force=True` |
| `No shared genes found` | Gene format mismatch | Check Ensembl vs HUGO format |
| `Label column not found` | Wrong column name | List columns with `list(adata.obs.columns)` |
| `ImportError: cellxgene-census` | Missing dependency | `pip install cellxgene-census` |

---

## Cloud Authentication

```python
import os

# Google Cloud Storage
os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = "/path/to/service-account.json"

# Amazon S3
os.environ["AWS_ACCESS_KEY_ID"] = "your-key"
os.environ["AWS_SECRET_ACCESS_KEY"] = "your-secret"

# Synapse
os.environ["SYNAPSE_AUTH_TOKEN"] = "your-token"
```
