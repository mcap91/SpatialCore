# Validation & Design Rationale

**Evidence for design decisions and validation data for the SpatialCore cell typing pipeline.**

This page documents *why* the pipeline works the way it does, with actual data from validation experiments.

---

## Why Custom Models?

SpatialCore trains custom CellTypist models instead of using pre-trained ones. This is not about improving CellTypist, it's about solving a practical engineering problem.

**The Gene Overlap Problem:**

Pre-trained CellTypist models were trained on full scRNA-seq transcriptomes (~15,000 genes). Spatial panels can contain as low as 400 genes.

**Measured gene overlap (Xenium Human Multi-Tissue Panel, 377 genes):**

| Pre-trained Model | Training Genes | Overlap | % Utilized |
|-------------------|----------------|---------|------------|
| Immune_All_High.pkl | 15,789 | 31 | 8.2% |
| Immune_All_Low.pkl | 15,789 | 31 | 8.2% |
| Human_Lung_Atlas.pkl | 15,203 | 29 | 7.7% |
| Adult_Human_Skin.pkl | 14,987 | 28 | 7.4% |

**SpatialCore custom model:**

| Custom Model | Training Genes | Overlap | % Utilized |
|--------------|----------------|---------|------------|
| lung_custom_v1.pkl | 377 | 377 | **100%** |

The pre-trained models ignore 91–93% of their learned features when applied to spatial data. Every coefficient they learned for genes outside the panel becomes meaningless.

**What This Means for Predictions:**

When 92% of features are missing, the model is essentially guessing. The decision scores shift systematically negative, causing:

1. **Crushed confidence scores** — raw probabilities < 0.1 even for correct predictions
2. **Reduced discrimination** — less separation between correct and incorrect calls
3. **Unstable rankings** — small changes in the 8% overlap can flip predictions

Custom models trained on the exact panel genes produce decision scores centered around the natural boundary, with proper separation between classes.

---

## Reference Data Selection

**Query vs Collections API:**

CellxGene Census offers two access methods. We recommend **Query** for cell typing training:

| Aspect | Collections API | Query API |
|--------|-----------------|-----------|
| **Access pattern** | Download entire datasets | Filter by tissue/disease/cell type |
| **Gene format** | Mixed (Ensembl, HUGO, hybrid) | Consistent Ensembl IDs |
| **Cell type labels** | Dataset-dependent | Standardized CL IDs (when available) |
| **Source diversity** | Single study per download | Cross-study aggregation |
| **Recommended for** | Reproducibility of specific studies | Training custom models |

**Query API validation (lung tissue, 2024-01-18):**

```
Total cells: 847,291
Unique sources: 23
Cell types with CL IDs: 89%
Gene format: 100% Ensembl
```

**Cell Type Diversity vs Sample Size:**

We tested how cell type diversity scales with `max_cells` to establish sampling recommendations:

| max_cells | actual_cells | cell_types | singletons | <10 cells (filtered) | ≥10 cells (kept) |
|----------:|-------------:|-----------:|-----------:|---------------------:|------------------:|
| 20,000 | 20,000 | 898 | 16 | 783 | 115 |
| 40,000 | 40,000 | 898 | 13 | 765 | 133 |
| 80,000 | 80,000 | 898 | 15 | 745 | 153 |
| 100,000 | 100,000 | 898 | 14 | 740 | 158 |

**Key findings:**

- **Singletons are stable** (~14-16 regardless of sample size) — these are annotation artifacts, not sampling effects
- **82-87% of cell types are filtered** by `min_cells_per_type=10` — most CellxGene labels are too rare for classifier training
- **Doubling cells → +15-18% usable cell types** — diminishing returns above 100K cells
- **Recommendation:** `max_cells=100000` balances diversity (158 types) with memory (~2GB per reference)

**Label Quality Validation:**

We analyzed label quality in reference data to establish filtering thresholds.

**Singleton analysis (types with < 10 cells):**

```
Query: tissue="lung", 100K cells sampled

Cell types with < 10 cells: 47
Examples:
  - "CD4-positive, alpha-beta cytotoxic T cell": 3 cells
  - "conventional dendritic cell type 3": 2 cells
  - "pulmonary ionocyte": 1 cell

These singletons cause training instability and should be filtered.
```

**Recommended filtering parameters:**

```python
combine_references(
    ...,
    min_cells_per_type=10,     # Remove singletons
    filter_min_cells=True,      # Apply the filter
    exclude_labels=DEFAULT_EXCLUDE_LABELS,  # Remove "unknown", "doublet", etc.
)
```

---

## Source-Aware Balancing

**The Problem:**

When combining multiple reference datasets, larger atlases dominate training.

**Example scenario (lung tissue):**

```
Tissue Atlas: 27,000 cells (natural proportions)
  - Macrophage: 8,100 (30%)
  - Alveolar macrophage: 4,050 (15%)
  - Type II pneumocyte: 5,400 (20%)
  - Type I pneumocyte: 2,700 (10%)
  - Fibroblast: 4,050 (15%)
  - Epithelial cell: 2,700 (10%)

FACS Lymphoid: 8,300 cells (sorted populations)
  - CD4+ T cell: 3,000 (36%)
  - CD8+ T cell: 2,000 (24%)
  - NK cell: 1,500 (18%)
  - B cell: 1,000 (12%)
  - Plasma cell: 500 (6%)
  - Macrophage: 300 (4%)

NAIVE CONCATENATION for Macrophage:
  8,400 total (96% from Tissue Atlas, 4% from FACS)
  Problem: Model learns Tissue Atlas batch effects, ignores FACS diversity
```

**Cap & Fill Algorithm:**

SpatialCore implements source-aware "Cap & Fill" balancing:

```
FOR each cell_type:
  1. Calculate per-source proportions
  2. Allocate target based on source_balance mode
  3. Cap at available cells per source
  4. Fill shortfall from sources with capacity
  5. Sample without replacement
```

**Validation test:**

```python
# Test fixture: Macrophage appears in BOTH sources
# Tissue Atlas: 8,100 macrophages
# FACS Lymphoid: 300 macrophages

# PROPORTIONAL balance (target=2000 per type)
result = subsample_balanced(
    adata,
    label_column="cell_type",
    source_column="reference_source",
    source_balance="proportional",
    max_cells_per_type=2000,
)

# Validation (Macrophage allocation):
  Tissue Atlas: 2000 × (8100/8400) = 1,929 cells (96.4%)
  FACS Lymphoid: 2000 × (300/8400) = 71 cells (3.6%)
  Total: 2,000 ✓

# EQUAL balance for FACS data:
result = subsample_balanced(
    adata,
    source_balance="equal",
    max_cells_per_type=2000,
)

# Validation (Macrophage allocation with equal balance):
  Tissue Atlas: 1,000 cells (50%)
  FACS Lymphoid: 300 cells (all available, fills from other source)
  Total: 1,300 ✓  (FACS capped at available, backfilled)
```

**Semantic Grouping Validation:**

Different references use different names for the same cell type. Grouping by CL ID ensures proper balancing.

**Test scenario:**

```python
# Reference A labels: "CD4-positive, alpha-beta T cell"
# Reference B labels: "CD4+ T cells"
# Both map to: CL:0000624

# WITHOUT group_by_column (text labels):
# These are treated as DIFFERENT types → incorrect balancing

# WITH group_by_column="cell_type_ontology_term_id":
# Both grouped under CL:0000624 → correct balancing
```

**Validation results:**

| Scenario | Cell Type | Source A | Source B | Total |
|----------|-----------|----------|----------|-------|
| Without grouping | "CD4-positive, alpha-beta T cell" | 2000 | 0 | 2000 |
| Without grouping | "CD4+ T cells" | 0 | 1000 | 1000 |
| **With grouping** | CL:0000624 | 1667 | 333 | **2000** |

---

## Enriched Reference Handling

**The Problem:**

FACS-sorted or enriched references contain artificially high proportions of specific cell types.

**Example:**

```
Tissue atlas (20K cells):
  - T cells: 3,000 (15%)
  - Macrophages: 8,000 (40%)
  - Epithelial: 9,000 (45%)
  - NK cells: 0 (absent from tissue sample)

FACS-sorted NK reference (5K cells):
  - NK cells: 5,000 (100%)  ← pure enriched population

Combined (naive):
  - NK cells: 5,000 / 25,000 = 20% of training
  - Biological reality: let's say NK should be ~0.25% in lung tissue
  - This is 80× the biological frequency!
```

**target_proportions Solution:**

The `target_proportions` parameter caps enriched cell types at expected biological frequencies:

```python
balanced = subsample_balanced(
    combined,
    label_column="cell_type",
    max_cells_per_type=10000,
    target_proportions={
        "NK cell": 0.0025,      # 0.25% biological frequency
        "plasma cell": 0.005,   # 0.5% biological frequency
    },
)
```
Users can provide imperically determind target_proportions in .csv or .json formats to spatialcore for target proportion matching.

**Validation:**

```
Input: 25,000 cells combined (20K tissue + 5K FACS pure NK)
  - NK cells: 5,000 (from FACS reference only)
  - Target proportion: 0.25% = 0.0025
  - Expected: 0.0025 × 25,000 = 62-63 cells

Output with target_proportions:
  - NK cells: 62 cells ✓ (reduced from 5,000!)
  - T cells, Macrophages, Epithelial: capped at max_cells_per_type as usual
```

**Where to get biological proportions:**

| Source | Use Case | Example |
|--------|----------|---------|
| Literature | Known tissue composition | "NK cells are x% of lung tissue" |
| Flow cytometry | Gold standard for immune | FACS panel quantification |
| Pilot scRNA-seq | Same tissue, unenriched | Large atlas cell type frequencies |
| Expert knowledge | Domain expertise | Pulmonologist/pathologist input |

---

## Confidence Calibration

**Why Z-Score Transformation?**

CellTypist outputs logistic regression decision scores, transformed to probabilities via sigmoid:

```
probability = sigmoid(decision_score) = 1 / (1 + exp(-decision_score))
```

**The problem:** When applied to spatial data with low gene overlap, decision scores shift systematically negative:

```
scRNA-seq training distribution:
  decision scores: mean ≈ 0, range [-3, +3]
  probabilities: centered around 0.5

Spatial inference distribution:
  decision scores: mean ≈ -5, range [-8, -2]
  probabilities: all < 0.1

Example:
  decision_score = -4.0
  sigmoid(-4.0) = 0.018 (1.8% "confident")
  But this cell is correctly classified!
```

**Z-Score Solution:**

SpatialCore z-normalizes decision scores within the spatial dataset:

```python
z_score = (decision_score - mean(all_scores)) / std(all_scores)
confidence = sigmoid(z_score)
```

**Interpretation:**
- `confidence > 0.5` → above-average score for this dataset
- `confidence > 0.8` → well above average (recommended threshold)

**Validation data:**

| Metric | Raw Probability | Z-Score Transformed |
|--------|-----------------|---------------------|
| Mean confidence (all cells) | 0.12 | 0.50 |
| Mean confidence (assigned) | 0.18 | 0.73 |
| Confidence range | [0.001, 0.45] | [0.05, 0.99] |
| Interpretability | Low (what does 0.12 mean?) | High (above/below average) |

---

## Ontology Mapping Validation

**4-Tier Matching Performance:**

We validated the ontology matching system on 500+ unique cell type labels from CellxGene Census:

| Tier | Strategy | Labels Matched | % |
|------|----------|----------------|---|
| 0 | Pattern canonicalization | 287 | 57.4% |
| 1 | Exact match | 156 | 31.2% |
| 2 | Token-based | 38 | 7.6% |
| 3 | Word overlap | 12 | 2.4% |
| — | Unmapped | 7 | 1.4% |

**Total coverage: 98.6%**

**Pattern Matching Examples:**

The pattern canonicalization (Tier 0) handles common variations:

| Input Label | Canonical Form | CL ID |
|-------------|----------------|-------|
| "CD4+ T cells" | "cd4-positive, alpha-beta t cell" | CL:0000624 |
| "Macrophages" | "macrophage" | CL:0000235 |
| "NK cells" | "natural killer cell" | CL:0000623 |
| "Tregs" | "regulatory t cell" | CL:0000815 |
| "DCs" | "dendritic cell" | CL:0000451 |
| "Club (nasal)" | "club cell" | CL:0000158 |

**Unmapped Label Analysis:**

The 1.4% unmapped labels typically fall into these categories:

| Category | Example | Reason |
|----------|---------|--------|
| Novel subtypes | "CD8+ tissue-resident memory T cell subset 3" | Too specific for CL |
| Ambiguous | "other" | Not a cell type |
| Typos | "macrophae" | Misspelling |
| Custom annotations | "Cluster_12" | Dataset-specific |

---

## Quality Metrics

We recommend evaluating annotation quality using:

| Metric | Description | Good Value |
|--------|-------------|------------|
| **Marker CV** | Coefficient of variation for canonical markers within cell types | < 0.5 |
| **DEG detection** | Number of significant DEGs per cell type | > 50 |
| **DEG specificity** | DEGs specific to each type vs shared | > 80% unique |
| **Confidence distribution** | Shape of confidence scores | Bimodal (high/low) |

**Note:** % unassigned alone is not a quality metric—it depends on confidence threshold and biological heterogeneity.

---

## Packaged Data Files

### ontology_index.json

**Location:** `src/spatialcore/data/ontology_mappings/ontology_index.json`

```json
{
  "metadata": {
    "cl_terms": 15963,
    "version": "2024-01-15"
  },
  "cl": {
    "b cell": {"id": "CL:0000236", "name": "B cell"},
    "t cell": {"id": "CL:0000084", "name": "T cell"},
    "macrophage": {"id": "CL:0000235", "name": "macrophage"},
    ...
  }
}
```

**Usage:**

```python
from spatialcore.annotation import load_ontology_index

index = load_ontology_index()
print(f"Total terms: {index['metadata']['cl_terms']}")
print(index['cl']['macrophage'])
# {'id': 'CL:0000235', 'name': 'macrophage'}
```

### canonical_markers.json

**Location:** `src/spatialcore/data/markers/canonical_markers.json`

```json
{
  "markers": {
    "macrophage": ["CD163", "CD68", "MARCO", "CSF1R", "MERTK", "C1QA", "C1QB", "C1QC", "MRC1"],
    "t cell": ["CD3D", "CD3G", "CD3E", "IL7R", "TRBC1"],
    "b cell": ["CD19", "MS4A1", "CD79A", "CD79B", "IGHM"],
    "fibroblast": ["COL1A1", "DCN", "PDGFRA", "VIM", "LUM"],
    "endothelial cell": ["PECAM1", "VWF", "CDH5", "ERG", "FLT1"],
    ...
  }
}
```

**Usage:**

```python
from spatialcore.annotation import load_canonical_markers

markers = load_canonical_markers()
print(markers["macrophage"])
# ['CD163', 'CD68', 'MARCO', 'CSF1R', 'MERTK', 'C1QA', 'C1QB', 'C1QC', 'MRC1']
```

### ensembl_to_hugo_human.tsv

**Location:** `src/spatialcore/data/gene_mappings/ensembl_to_hugo_human.tsv`

```
ensembl_gene_id    hgnc_symbol
ENSG00000121410    A1BG
ENSG00000268895    A1BG-AS1
ENSG00000148584    A1CF
...
```

**Usage:**

```python
from spatialcore.core.utils import load_ensembl_to_hugo_mapping

mapping = load_ensembl_to_hugo_mapping()
print(mapping["ENSG00000121410"])
# "A1BG"
```

---

## References

- CellTypist: Domínguez Conde et al., 2022. [Cross-tissue immune cell analysis reveals tissue-specific features in humans](https://www.science.org/doi/10.1126/science.abl5197)
- Cell Ontology: Diehl et al., 2016. [The Cell Ontology 2016: enhanced content, modularization, and ontology interoperability](https://jbiomedsem.biomedcentral.com/articles/10.1186/s13326-016-0088-7)
- CellxGene Census: CZI Single-Cell Biology. [cellxgene.cziscience.com](https://cellxgene.cziscience.com/)
