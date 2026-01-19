# Benchmark: Custom vs Pre-trained Cell Type Annotation

This benchmark demonstrates SpatialCore's approach to cell typing on spatial transcriptomics data, comparing it against out-of-box CellTypist with a pre-trained model.

## The Problem

Pre-trained cell type classifiers face two fundamental challenges on spatial data:

1. **Gene overlap mismatch** - Pre-trained models learn from RNA-seq (~20,000 genes). Spatial panels have 300-1,000 targeted genes. Only 5-10% of learned features are available at inference time.

2. **Reference bias** - Single-source training data introduces dataset-specific biases. Cell types from dominant sources are over-represented.

SpatialCore addresses both problems with:

- **Panel-specific training** - Train on exactly the genes in your spatial panel (100% overlap)
- **CellxGene integration** - Download tissue-matched references from 60M+ cells
- **Source-aware balancing** - `subsample_balanced()` ensures fair representation across sources

---

## Dataset

| Property | Value |
|----------|-------|
| **Platform** | 10x Genomics Xenium |
| **Tissue** | Human lung (NSCLC) |
| **Cells** | 93,162 |
| **Panel genes** | 518 |

---

## Step 1: Acquiring Reference Data

SpatialCore integrates with CZ CELLxGENE Discover Census to download tissue-matched scRNA-seq references. This provides high-quality, ontology-standardized training data.

```python
from spatialcore.annotation import acquire_reference
from pathlib import Path

REFERENCE_DIR = Path("references/cellxgene/lung")
REFERENCE_DIR.mkdir(parents=True, exist_ok=True)

# Download healthy lung tissue (~100k cells)
acquire_reference(
    source="cellxgene://?tissue=lung&disease=normal",
    output=REFERENCE_DIR / "healthy_lung.h5ad",
    max_cells=100000,
)

# Download NSCLC tumor samples (~100k cells)
acquire_reference(
    source="cellxgene://?tissue=lung&disease=non-small cell lung carcinoma",
    output=REFERENCE_DIR / "nsclc.h5ad",
    max_cells=100000,
)
```

**Why CellxGene?**

- 60M+ cells from 700+ curated datasets
- Standardized Cell Ontology (CL) labels
- Tissue and disease filtering
- Programmatic access via Census API

For validation of the CellxGene download and subsampling approach, see [validation.md](validation.md).

---

## Step 2: The Baseline (Standalone CellTypist)

First, we establish a baseline using out-of-box CellTypist with the pre-trained Human Lung Atlas model.

```python
import scanpy as sc
import celltypist
from celltypist import models

# Load spatial data
adata = sc.read_h5ad("xenium_lung_cancer_clustered.h5ad")

# Download and load pre-trained model
models.download_models(model="Human_Lung_Atlas.pkl")
model = models.Model.load(model="Human_Lung_Atlas.pkl")

# Check gene overlap
model_genes = set(model.features)
query_genes = set(adata.var_names)
overlap = model_genes & query_genes
overlap_pct = 100 * len(overlap) / len(model_genes)

print(f"Model genes: {len(model_genes):,}")
print(f"Query genes: {len(query_genes):,}")
print(f"Overlap: {len(overlap):,} ({overlap_pct:.1f}%)")
# Output: Overlap: 356 (7.1%)

# Run annotation
predictions = celltypist.annotate(adata, model=model, majority_voting=False)

# Check confidence
confidence = predictions.probability_matrix.max(axis=1).values
low_conf = (confidence < 0.5).mean()
print(f"Below 0.5 threshold: {low_conf:.1%}")
# Output: Below 0.5 threshold: 98.0%
```

**Result: 98% of cells are unassigned.**

With only 7% gene overlap, the model lacks the features it needs to make confident predictions.

---

## Step 3: The SpatialCore Solution

SpatialCore's `train_and_annotate()` solves both problems in a single API call.

### The Full Pipeline

```python
from spatialcore.annotation import train_and_annotate, discover_training_data
import scanpy as sc

# Load spatial data
adata = sc.read_h5ad("xenium_lung_cancer_clustered.h5ad")

# Discover available references
datasets = discover_training_data("references/cellxgene/lung")
reference_paths = [ds.path for ds in datasets]

# Train custom model and annotate
adata = train_and_annotate(
    adata,
    references=reference_paths,
    tissue="lung",
    balance_strategy="proportional",
    max_cells_per_type=10000,
    max_cells_per_ref=100000,
    confidence_threshold=0.8,
    model_output="models/lung_nsclc_custom_v1.pkl",
    plot_output="plots/",
    add_ontology=True,
    generate_plots=True,
)

# Check results
print(f"Cell types: {adata.obs['cell_type'].nunique()}")
print(f"Mean confidence: {adata.obs['cell_type_confidence'].mean():.3f}")
print(f"Unassigned: {(adata.obs['cell_type'] == 'Unassigned').mean():.1%}")
# Output: Unassigned: 0.5%
```

### What `train_and_annotate()` Does

The pipeline executes 9 stages:

1. **Extract panel genes** - Gets gene names from spatial data
2. **Load references** - Combines multiple h5ad files with Ensembl-to-HUGO normalization
3. **Fill ontology IDs** - Maps cell type labels to Cell Ontology (CL) terms
4. **Source-aware balancing** - `subsample_balanced()` with "Cap & Fill" strategy
5. **Train CellTypist model** - SGD classifier on balanced, panel-subset data
6. **Annotate spatial data** - Apply model with z-score confidence
7. **Apply threshold** - Mark low-confidence cells as Unassigned
8. **Map to ontology** - Add CL IDs to predictions
9. **Generate plots** - DEG heatmap, 2D validation, confidence plots

### Source-Aware Balancing

The `subsample_balanced()` function prevents reference bias through "Cap & Fill" balancing:

```python
from spatialcore.annotation import subsample_balanced

# Balance training data across sources and cell types
balanced = subsample_balanced(
    combined_references,
    label_column="cell_type_ontology_label",
    group_by_column="cell_type_ontology_term_id",  # Group by CL ID
    source_column="reference_source",
    source_balance="proportional",
    max_cells_per_type=10000,
    copy=True,
)
```

**Why this matters:**

- **Source balance** - Each reference contributes proportionally to each cell type
- **CL ID grouping** - Semantic synonyms (e.g., "CD4+ T cell" and "CD4-positive, alpha-beta T cell") are grouped together
- **Cell type balance** - Rare types get adequate representation

For detailed scenarios and validation, see [validation.md](validation.md).

For the full API reference, see [pipeline.md](pipeline.md).

---

## Results

### Summary Table

| Metric | Standalone | SpatialCore | Improvement |
|--------|-----------|-------------|-------------|
| Gene Overlap (%) | 7.1% | 100% | 14x |
| Unknown Cells (%) | 98.0% | 0.5% | 196x |
| Marker CV | 1.77 | 1.43 | 19% lower |
| Marker log2FC | 1.50 | 1.95 | 30% higher |
| DEG log2FC | 3.93 | 4.04 | 3% higher |
| Marker Purity (%) | 39.0% | 46.3% | 19% higher |
| Contamination | 0.85 | 0.72 | 15% lower |

**SpatialCore wins on all 7 metrics.**

![Benchmark Summary](images/benchmark_summary_table.png)

### Gene Overlap

The Human Lung Atlas model has 5,017 features. Our Xenium panel has 518 genes. Only 356 genes overlap (7.1%).

SpatialCore trains on the exact panel genes, achieving 100% overlap by design.

![Gene Overlap Comparison](images/gene_overlap_comparison.png)

### Unassigned Rate

With 7% gene overlap, standalone CellTypist marks 98% of cells as unassigned (below 0.5 confidence threshold).

SpatialCore, with 100% overlap and z-score normalization, marks only 0.5% as unassigned (below 0.8 threshold).

![Unknown Cell Rates](images/unknown_celltype_calls.png)

### Confidence Distribution

The two methods produce fundamentally different confidence distributions:

![Confidence Distributions](images/confidence_distribution_comparison.png)

- **Standalone**: Raw probabilities cluster near 0 (lacks informative features)
- **SpatialCore**: Z-score normalized values spread across [0,1]

### Biological Validation

All biological metrics are evaluated on **all cells** without threshold filtering for fair comparison.

**Marker Expression Consistency (CV)** - Lower is better:

![Marker CV Comparison](images/marker_cv_comparison.png)

**Marker Specificity (log2FC)** - Higher is better:

![Marker Fold Change](images/marker_foldchange_comparison.png)

**Canonical Marker Recovery** - Higher is better:

![Canonical Marker Recovery](images/canonical_marker_presence.png)

**DEG Effect Size** - Higher is better:

![DEG Effect Size](images/deg_effect_size_comparison.png)

**Marker Purity** - Higher is better:

![Marker Purity](images/marker_purity_comparison.png)

**Cross-Type Contamination** - Lower is better:

![Contamination Score](images/contamination_comparison.png)

### Cell Type Inventory

![Cell Type Inventory](images/cell_type_inventory.png)

---

## Validation Plots

Both methods generate identical validation plot suites for fair comparison.

### SpatialCore

| Plot | Description |
|------|-------------|
| ![DEG Heatmap](images/lung_celltyping_deg_heatmap.png) | Top DEGs per cell type |
| ![2D Validation](images/lung_celltyping_2d_validation.png) | Canonical marker expression |
| ![Confidence](images/lung_celltyping_confidence.png) | Confidence by cell type |
| ![Ontology](images/lung_celltyping_ontology_mapping.png) | Ontology mapping table |

### Standalone CellTypist

| Plot | Description |
|------|-------------|
| ![DEG Heatmap](images/standalone_deg_heatmap.png) | Top DEGs per cell type |
| ![2D Validation](images/standalone_2d_validation.png) | Canonical marker expression |
| ![Confidence](images/standalone_confidence.png) | Confidence by cell type |
| ![Ontology](images/standalone_ontology_mapping.png) | Ontology mapping table |

---

## Conclusion

The gene overlap problem is real and measurable. Pre-trained models with 7% feature overlap produce 98% unassigned cells.

SpatialCore solves this with three innovations:

1. **CellxGene integration** - Download tissue-matched references from 60M+ cells
2. **Source-aware balancing** - Fair representation across sources and cell types
3. **Panel-specific training** - 100% gene overlap by design

The result: 0.5% unassigned rate and improved biological coherence across all metrics.

---

## References

**Spatial Data**

- 10x Genomics (2023). FFPE Human Lung Cancer with Immuno-Oncology Panel. [10xgenomics.com/datasets](https://www.10xgenomics.com/datasets/ffpe-human-lung-cancer-data-with-human-immuno-oncology-profiling-panel-and-custom-add-on-1-standard)

**CellTypist**

- Dominguez Conde C, et al. (2022). Cross-tissue immune cell analysis reveals tissue-specific features in humans. *Science*. [DOI: 10.1126/science.abl5197](https://doi.org/10.1126/science.abl5197)
- GitHub: [github.com/Teichlab/celltypist](https://github.com/Teichlab/celltypist) | License: Apache 2.0

**CellxGene Census**

- CZI Single-Cell Biology, et al. (2023). CZ CELLxGENE Discover. *bioRxiv*. [DOI: 10.1101/2023.10.30.563174](https://doi.org/10.1101/2023.10.30.563174)
- Docs: [chanzuckerberg.github.io/cellxgene-census](https://chanzuckerberg.github.io/cellxgene-census/) | License: CC-BY 4.0
