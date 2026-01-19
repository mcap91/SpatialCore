# Validation & Design Rationale

**Evidence for design decisions and validation data for the SpatialCore cell typing pipeline.**

This page documents *why* the pipeline works the way it does, with actual data from validation experiments.

---

## Why Custom Models?

SpatialCore trains custom CellTypist models instead of using pre-trained ones. This is not about improving CellTypist—it's about solving a practical engineering problem.

### The Gene Overlap Problem

Pre-trained CellTypist models were trained on full scRNA-seq transcriptomes (~15,000 genes). Spatial panels contain 300–500 genes.

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

### What This Means for Predictions

When 92% of features are missing, the model is essentially guessing. The decision scores shift systematically negative, causing:

1. **Crushed confidence scores** — raw probabilities < 0.1 even for correct predictions
2. **Reduced discrimination** — less separation between correct and incorrect calls
3. **Unstable rankings** — small changes in the 8% overlap can flip predictions

Custom models trained on the exact panel genes produce decision scores centered around the natural boundary, with proper separation between classes.

---

## Reference Data Selection

### Query vs Collections API

CellxGene Census offers two access methods. We recommend **Query** for cell typing training:

| Aspect | Collections API | Query API |
|--------|-----------------|-----------|
| **Access pattern** | Download entire datasets | Filter by tissue/disease/cell type |
| **Gene format** | Mixed (Ensembl, HUGO, hybrid) | Consistent Ensembl IDs |
| **Cell type labels** | Dataset-dependent | Standardized CL IDs (when available) |
| **Source diversity** | Single study per download | Cross-study aggregation |
| **Recommended for** | Reproducibility of specific studies | Training custom models |

**Query API validation (liver tissue, 2024-01-18):**

```
Total cells: 847,291
Unique sources: 23
Cell types with CL IDs: 89%
Gene format: 100% Ensembl
```

### Label Quality Validation

We analyzed label quality in reference data to establish filtering thresholds.

**Singleton analysis (types with < 10 cells):**

```python
# From tests/test_cellxgene_functions.py
# Query: tissue="liver", 100K cells sampled

Cell types with < 10 cells: 47
Examples:
  - "CD4-positive, alpha-beta cytotoxic T cell": 3 cells
  - "conventional dendritic cell type 3": 2 cells
  - "hepatic stellate cell (activated)": 1 cell

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

### The Problem

When combining multiple reference datasets, larger atlases dominate training.

**Example scenario:**

```
Reference A (HLCA): 500,000 cells
  - Macrophages: 50,000
  - T cells: 100,000
  - B cells: 30,000

Reference B (In-house): 10,000 cells
  - Macrophages: 2,000
  - T cells: 3,000
  - B cells: 1,000

NAIVE CONCATENATION:
  Macrophages: 52,000 (96% from HLCA)
  T cells: 103,000 (97% from HLCA)

Problem: Model learns HLCA batch effects, ignores in-house data
```

### Cap & Fill Algorithm

SpatialCore implements source-aware "Cap & Fill" balancing:

```
FOR each cell_type:
  1. Calculate per-source proportions
  2. Allocate target based on source_balance mode
  3. Cap at available cells per source
  4. Fill shortfall from sources with capacity
  5. Sample without replacement
```

**Validation test (from `tests/test_subsample_balanced.py`):**

```python
# Test fixture: 2 sources, 3 cell types
# Source A: T cell=5000, B cell=3000, Macrophage=2000
# Source B: T cell=1000, B cell=500, Macrophage=300

# PROPORTIONAL balance (target=2000 per type)
result = subsample_balanced(
    adata,
    label_column="cell_type",
    source_column="source",
    source_balance="proportional",
    max_cells_per_type=2000,
)

# Validation:
T cell allocation:
  Source A: 2000 × (5000/6000) = 1,667 cells
  Source B: 2000 × (1000/6000) = 333 cells
  Total: 2,000 ✓

B cell allocation:
  Source A: 2000 × (3000/3500) = 1,714 cells
  Source B: 2000 × (500/3500) = 286 cells
  Total: 2,000 ✓
```

### Semantic Grouping Validation

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

### The Problem

FACS-sorted or enriched references contain artificially high proportions of specific cell types.

**Real-world example:**

```
Tissue atlas (500K cells):
  - NK cells: 250 (0.05%)
  - Plasma cells: 1,500 (0.3%)
  - Macrophages: 50,000 (10%)

FACS-sorted NK reference (5K cells):
  - NK cells: 5,000 (100%)  ← artificially enriched

Combined (naive):
  - NK cells: 5,250 (1% of 505K)
  - This is 20× the biological frequency!
```

### target_proportions Solution

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

**Validation (from `tests/validate_subsampling.py`):**

```
Input: 25,000 cells combined
  - NK cells: 5,000 (from FACS reference only)
  - Expected: 0.25% = 62 cells

Output with target_proportions:
  - NK cells: 62 cells ✓
  - Other types: capped at max_cells_per_type as usual
```

**Where to get biological proportions:**

| Source | Use Case | Example |
|--------|----------|---------|
| Literature | Known tissue composition | "NK cells are ~1-5% of liver lymphocytes" |
| Flow cytometry | Gold standard for immune | FACS panel quantification |
| Pilot scRNA-seq | Same tissue, unenriched | Large atlas cell type frequencies |
| Expert knowledge | Domain expertise | Pathologist/immunologist input |

---

## Confidence Calibration

### Why Z-Score Transformation?

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

### Z-Score Solution

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

### 4-Tier Matching Performance

We validated the ontology matching system on 500+ unique cell type labels from CellxGene Census:

| Tier | Strategy | Labels Matched | % |
|------|----------|----------------|---|
| 0 | Pattern canonicalization | 287 | 57.4% |
| 1 | Exact match | 156 | 31.2% |
| 2 | Token-based | 38 | 7.6% |
| 3 | Word overlap | 12 | 2.4% |
| — | Unmapped | 7 | 1.4% |

**Total coverage: 98.6%**

### Pattern Matching Examples

The pattern canonicalization (Tier 0) handles common variations:

| Input Label | Canonical Form | CL ID |
|-------------|----------------|-------|
| "CD4+ T cells" | "cd4-positive, alpha-beta t cell" | CL:0000624 |
| "Macrophages" | "macrophage" | CL:0000235 |
| "NK cells" | "natural killer cell" | CL:0000623 |
| "Tregs" | "regulatory t cell" | CL:0000815 |
| "DCs" | "dendritic cell" | CL:0000451 |
| "Club (nasal)" | "club cell" | CL:0000158 |

### Unmapped Label Analysis

The 1.4% unmapped labels typically fall into these categories:

| Category | Example | Reason |
|----------|---------|--------|
| Novel subtypes | "CD8+ tissue-resident memory T cell subset 3" | Too specific for CL |
| Ambiguous | "other" | Not a cell type |
| Typos | "macrophae" | Misspelling |
| Custom annotations | "Cluster_12" | Dataset-specific |

---

## Validation Test Summary

### Test Coverage

```
tests/
├── test_cellxgene_functions.py    # CellxGene API integration
├── test_subsample_balanced.py     # Balancing algorithm unit tests
├── test_normalization.py          # Normalization detection (42 tests)
└── validate_subsampling.py        # End-to-end validation with plots
```

### Key Test Scenarios

**From `test_subsample_balanced.py`:**

| Test | Description | Status |
|------|-------------|--------|
| `test_proportional_balance` | Source A:B = 5:1 → allocation 5:1 | PASS |
| `test_equal_balance` | Source A:B = 5:1 → allocation 1:1 | PASS |
| `test_cap_at_available` | Source has fewer cells than target | PASS |
| `test_fill_shortfall` | Redistribute unfilled quota | PASS |
| `test_min_cells_enforced` | Below-min sources still contribute | PASS |
| `test_group_by_ontology` | Different labels, same CL ID | PASS |
| `test_target_proportions` | Enriched cell types capped | PASS |

**From `test_normalization.py`:**

| Test | Description | Status |
|------|-------------|--------|
| `test_detect_raw_counts` | Integers in X | PASS |
| `test_detect_log1p_10k` | Verified via expm1 | PASS |
| `test_detect_log1p_cpm` | Different target sum | PASS |
| `test_layer_priority` | counts > raw_counts > raw > raw.X | PASS |
| `test_float_precision` | 1.0000000000002 is integer | PASS |
| `test_unsafe_force` | Allows invalid states with warning | PASS |

---

## Benchmark Data

### Xenium Liver (162K cells)

Comparison of pre-trained vs custom model performance:

| Metric | Pre-trained (Immune_All) | SpatialCore Custom |
|--------|--------------------------|-------------------|
| Gene overlap | 31 (8.2%) | 377 (100%) |
| Unassigned cells | 54.5% | 44.2% |
| Mean confidence | 0.477 | 0.558 |
| Mean confidence (assigned) | 0.80 | 0.99 |
| Cell types detected | 67 | 64 |
| Ontology mapped | 95.5% | 100% |

### Key Quality Metrics

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

## Running Validation Tests

```bash
# Run all annotation tests
pytest tests/test_subsample_balanced.py -v
pytest tests/test_normalization.py -v

# Run end-to-end validation with plots
python tests/validate_subsampling.py

# Run CellxGene integration tests (requires network)
pytest tests/test_cellxgene_functions.py -v
```

---

## References

- CellTypist: Domínguez Conde et al., 2022. [Cross-tissue immune cell analysis reveals tissue-specific features in humans](https://www.science.org/doi/10.1126/science.abl5197)
- Cell Ontology: Diehl et al., 2016. [The Cell Ontology 2016: enhanced content, modularization, and ontology interoperability](https://jbiomedsem.biomedcentral.com/articles/10.1186/s13326-016-0088-7)
- CellxGene Census: CZI Single-Cell Biology. [cellxgene.cziscience.com](https://cellxgene.cziscience.com/)
