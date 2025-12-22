# Ontology Mapping: Standalone Implementation Guide

**Date**: 2025-11-21 (Updated 2025-12-16)
**Status**: Production-ready, documented for standalone extraction
**Context**: Cell type label → ontology code conversion

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Architecture Overview](#architecture-overview)
3. [Core Components](#core-components)
4. [Implementation Details](#implementation-details)
5. [Standalone Extraction Guide](#standalone-extraction-guide)
6. [API Reference](#api-reference)
7. [Extension Guide](#extension-guide)
8. [Troubleshooting](#troubleshooting)
9. [Performance Benchmarks](#performance-benchmarks)
10. [References](#references)

---

## Executive Summary

### What This Does

Converts human-readable biological labels (e.g., "CD4+ T cells", "B cell", "Macrophage") to standardized ontology codes (e.g., `CL:0000624`, `CL:0000236`, `CL:0000235`) using a multi-tier matching system.

### Key Features

| Feature | Description |
|---------|-------------|
| **93%+ match rate** | High accuracy with preserved specificity |
| **Zero external APIs** | Fully offline operation |
| **Fast** | ~2s index load, milliseconds per label |
| **Multi-ontology** | CL (cells), NCIT (pathology), UBERON (anatomy) |
| **Hierarchy-aware** | Voting across multiple annotation sources |
| **Extensible** | Pattern-based system for custom mappings |

### Decision: Custom Hybrid Approach

We evaluated and rejected:
- **CellOntologyMapper** (7.5GB Docker bloat, requires internet)
- **text2term** (general-purpose, not cell-type optimized)
- **CellO** (gene expression-based, different use case)

Our custom implementation provides better performance with zero dependencies beyond standard Python libraries.

---

## Architecture Overview

### System Diagram

```
┌─────────────────────────────────────────────────────────────────────┐
│                        ONTOLOGY MAPPING SYSTEM                       │
├─────────────────────────────────────────────────────────────────────┤
│                                                                      │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────────────┐   │
│  │ Input Labels │───▶│ Tier Matcher │───▶│ Ontology Codes       │   │
│  │              │    │              │    │                      │   │
│  │ "CD4+ T"     │    │ T0: Pattern  │    │ CL:0000624           │   │
│  │ "B cell"     │    │ T1: Exact    │    │ CL:0000236           │   │
│  │ "NK cells"   │    │ T2: Token    │    │ CL:0000623           │   │
│  └──────────────┘    │ T3: Overlap  │    └──────────────────────┘   │
│                      └──────────────┘                                │
│                             │                                        │
│                             ▼                                        │
│                    ┌──────────────────┐                              │
│                    │ Ontology Index   │                              │
│                    │ (JSON, 28.8 MB)  │                              │
│                    │                  │                              │
│                    │ CL: 2,847 terms  │                              │
│                    │ NCIT: 172,000+   │                              │
│                    │ UBERON: 15,000+  │                              │
│                    └──────────────────┘                              │
│                                                                      │
└─────────────────────────────────────────────────────────────────────┘
```

### Three-Layer System

1. **Pre-built JSON Index** (`ontology_index.json`)
   - 245,372 terms from CL, NCIT, UBERON
   - 28.8 MB, loads in ~2 seconds
   - O(1) lookups vs 90+ minutes parsing OWL files

2. **Annotation Type Routing**
   - `"cell_type"` → Uses **CL only** (most specific)
   - `"pathology"` → Uses **NCIT first**, then CL
   - `"anatomy"` → Uses **UBERON first**, then CL
   - `"all"` → Searches CL, NCIT, UBERON in priority order

3. **Sequential Tier Matching**
   - **Tier 0**: Pattern canonicalization (known abbreviations)
   - **Tier 1**: Exact/partial match (0.9-1.0)
   - **Tier 2**: Token-based match (0.60-0.80)
   - **Tier 3**: Word overlap fallback (0.5-0.7)

---

## Core Components

### File Structure

**Primary Implementation** (this repo):
```
spatial-biology/
├── spacebio/
│   └── ontology.py            # Core matching logic (~1,000 lines)
├── docs/
│   ├── ontology_index.json    # Pre-built term index (28.8 MB)
│   └── ONTOLOGY_MAPPING_RESEARCH.md  # This file
└── configs/
    └── canonical_markers.json # Curated marker panel for benchmarking
```

**Standalone Package** (for extraction):
```
ontology_mapping/
├── ontology_utils.py          # Core matching logic
├── convert_to_ontology_codes.py   # Temporal activity wrapper
├── data/
│   └── ontology_index.json    # Pre-built term index (28.8 MB)
└── docs/
    └── ONTOLOGY_MAPPING_RESEARCH.md
```

### Key Functions

| Function | Purpose | Location |
|----------|---------|----------|
| `search_ontology_index()` | Main search API | `spacebio/ontology.py` |
| `load_ontology_index()` | Load JSON index | `spacebio/ontology.py` |
| `map_celltype_to_ontology()` | Convenience wrapper | `spacebio/ontology.py` |
| `add_ontology_ids()` | AnnData integration | `spacebio/ontology.py` |
| `hierarchy_aware_vote()` | Multi-reference voting | `spacebio/ontology.py` |
| `extract_biological_tokens()` | Token extraction | `spacebio/ontology.py` |
| `_apply_cell_type_patterns()` | Pattern matching (Tier 0) | `spacebio/ontology.py` |
| `_score_match()` | Multi-tier scoring | `spacebio/ontology.py` |
| `validate_cl_term()` | Check term exists in CL | `spacebio/ontology.py` |

### Dependencies

**Required (Core)**:
```
re          # Standard library - regex
json        # Standard library - JSON parsing
pathlib     # Standard library - path handling
typing      # Standard library - type hints
```

**Required (Index Building Only)**:
```
owlready2   # OWL file parsing (pip install owlready2)
```

**Optional (Hierarchy Expansion)**:
```
owlready2   # For expand_term(), get_ancestors(), is_ancestor_of()
```

---

## Implementation Details

### Tier 0: Pattern Canonicalization

**Purpose**: Transform known abbreviations/variants to canonical search terms BEFORE fuzzy matching.

**Why First**: Prevents fuzzy matching from giving wrong results for known patterns.

```python
CELL_TYPE_PATTERNS = {
    # CD marker-based cell types
    r"t\s*cells?,?\s*cd4\+?|cd4\+?\s*t": "cd4-positive, alpha-beta t cell",
    r"t\s*cells?,?\s*cd8\+?|cd8\+?\s*t": "cd8-positive, alpha-beta t cell",

    # T cell subtypes
    r"t.*helper.*17|th17": "t-helper 17 cell",
    r"regulatory.*t|t.*regulatory|treg": "regulatory t cell",
    r"gamma.*delta.*t|gammadelta.*t|gdt": "gamma-delta t cell",

    # B cell variants
    r"^b\s*cell|^b\s+cells?$": "b cell",
    r"germinal.*center.*b|gc.*b\s*cell": "germinal center B cell",

    # NK cells
    r"\bnk\s*cell|\bnatural\s*killer": "natural killer cell",

    # Myeloid lineage
    r"^pdc\b|plasmacytoid\s*dc": "plasmacytoid dendritic cell",
    r"dendritic\s*cells?|^dc\d?\b|^cdc\d?\b": "dendritic cell",

    # ... 100+ patterns total
}
```

**Adding New Patterns**:
1. Identify failing label
2. Find correct CL term in index
3. Add pattern: `r'pattern': 'canonical term'`
4. Test with `search_ontology_index(["test label"])`

### Tier 1: Exact & Direct Matches (0.9-1.0)

```python
# Exact match
"b cell" == "b cell" → 1.0

# Cleaned match (symbols removed)
"cd4+ t cells" → "cd4 t cells" → 0.95

# Contains match (substring, min 4 chars)
"helper t" in "cd4-positive helper t cell" → 0.9
```

### Tier 2: Smart Token Extraction (0.60-0.80)

**Key Innovation**: Filter out generic terms like "cell/cells" to prevent false matches.

```python
def _extract_biological_tokens(label: str) -> Dict[str, List[str]]:
    """
    Extract key biological identifiers from label.

    Returns:
        {
            'markers': ['cd4', 'cd8', 'cd19'],     # CD markers
            'proteins': ['igg', 'iga', 'spp1'],    # Immunoglobulins, genes
            'core_words': ['helper', 't', 'plasma'], # Main bio terms
            'modifiers': ['positive', 'mature']    # Descriptors (ignored)
        }
    """
    # CRITICAL: Filter out generic terms
    GENERIC_TERMS = {'cell', 'cells'}  # Too generic for matching

    # Keep meaningful short identifiers
    MEANINGFUL_SHORT_WORDS = {'b', 't', 'nk', 'dc', 'ec', 've', 'ta', 'm1', 'm2'}
```

**Why Filter "cell/cells"?**
- Old: "CD4+ T cells" → tokens: `['cd4', 't', 'cells']` → matched "bronchioalveolar stem **cells**"
- New: "CD4+ T cells" → tokens: `['cd4', 't']` → correct match

**Scoring Logic**:
```python
# All core words present
if all(word in term_words for word in tokens["core_words"]):
    base_score = 0.70

    # Penalty for ambiguous short tokens
    if len(tokens["core_words"]) == 1 and len(tokens["core_words"][0]) <= 2:
        base_score -= 0.15  # Drop to 0.55 (below threshold)

    # Penalty for unwanted developmental prefixes
    unwanted_prefixes = ["pro", "pre", "post", "immature", "ecto", "endo"]
    if has_unwanted_prefix and not label_has_prefix:
        base_score -= 0.15

    # Bonus if markers also match
    if tokens["markers"] and any(m in term_clean for m in tokens["markers"]):
        base_score = max(base_score, 0.75)

    # Specificity boosts for number matching, TCR type, word count
    # ...
```

### Tier 3: Word Overlap Fallback (0.5-0.7)

Jaccard similarity on word sets:
```python
label_words = {"cd4", "t"}
term_words = {"cd4", "positive", "alpha", "beta", "t", "cell"}
common = label_words & term_words
jaccard = len(common) / len(label_words | term_words)
score = 0.5 + (0.2 * jaccard)  # Range: 0.5-0.7
```

### Hierarchy-Aware Voting

For multi-reference annotation (e.g., SingleR with multiple atlases):

```python
def hierarchy_aware_vote(
    cl_codes: List[str],
    scores: Optional[List[float]] = None,
    use_hierarchy: bool = True,
) -> Tuple[str, float]:
    """
    Algorithm:
    1. Group codes by lineage (codes in same ancestor chain)
    2. Within a lineage, prefer most specific (deepest) code
    3. A parent code counts as a vote for any child in that lineage
    4. Confidence = fraction of refs that agree (directly or via hierarchy)

    Example:
        Input: ["CL:0000084", "CL:0000624", "CL:0000624"]
               (T cell,       CD4+ T,       CD4+ T)

        Lineage grouping: {CL:0000084, CL:0000624} - same lineage
        T cell is ancestor of CD4+ T cell

        Result: CL:0000624 (most specific), confidence=1.0 (all agree via hierarchy)
    """
```

---

## Standalone Extraction Guide

### Minimal Standalone Package

To extract this as a standalone package, you need:

#### 1. Core Files (Required)

```python
# ontology_mapper/__init__.py
from .matcher import search_ontology_index, load_ontology_index
from .patterns import CELL_TYPE_PATTERNS
from .voting import hierarchy_aware_vote

__all__ = [
    'search_ontology_index',
    'load_ontology_index',
    'hierarchy_aware_vote',
    'CELL_TYPE_PATTERNS',
]
```

#### 2. Matcher Module

```python
# ontology_mapper/matcher.py
"""
Core ontology matching functionality.
Extract from ontology_utils.py lines 1317-1495 (search_ontology_index)
Plus supporting functions:
- load_ontology_index (lines 793-893)
- _extract_biological_tokens (lines 905-988)
- _apply_cell_type_patterns (lines 1130-1150)
- _score_match (lines 1153-1314)
"""

from typing import List, Dict, Optional, Any
from pathlib import Path
import json
import re

def load_ontology_index(index_uri: str = None) -> Dict[str, Dict[str, Dict[str, str]]]:
    """
    Load pre-built ontology index from JSON file.

    Args:
        index_uri: Path to ontology_index.json (default: bundled)

    Returns:
        {
            "cl": {"b cell": {"id": "CL:0000236", "name": "B cell"}, ...},
            "ncit": {...},
            "uberon": {...}
        }
    """
    if index_uri is None:
        # Use bundled index
        module_dir = Path(__file__).parent
        possible_paths = [
            module_dir / "data" / "ontology_index.json",
            module_dir.parent / "data" / "ontology_index.json",
        ]
        for path in possible_paths:
            if path.exists():
                index_uri = str(path)
                break
        else:
            raise FileNotFoundError("Ontology index not found")

    with open(index_uri, "r") as f:
        index_data = json.load(f)

    return {
        "cl": index_data.get("cl", {}),
        "ncit": index_data.get("ncit", {}),
        "uberon": index_data.get("uberon", {}),
    }


def search_ontology_index(
    labels: List[str],
    ontology_index: Dict = None,
    index_uri: str = None,
    ontologies: List[str] = None,
    annotation_type: str = "cell_type",
    min_score: float = 0.7,
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Search ontology index for matching terms.

    Args:
        labels: List of human-readable labels to search
        ontology_index: Pre-loaded index (optional)
        index_uri: Path to index JSON (optional)
        ontologies: Which ontologies to search (overrides annotation_type)
        annotation_type: "cell_type" | "pathology" | "anatomy" | "all"
        min_score: Minimum match score (0.0-1.0)

    Returns:
        {
            "B cell": [{"id": "CL:0000236", "name": "B cell", "score": 1.0, ...}],
            ...
        }
    """
    # Implementation from ontology_utils.py lines 1364-1493
    ...
```

#### 3. Patterns Module

```python
# ontology_mapper/patterns.py
"""
Cell type pattern definitions.
Extract CELL_TYPE_PATTERNS from ontology_utils.py lines 991-1127
"""

CELL_TYPE_PATTERNS = {
    # === Lymphoid lineage ===
    r"t\s*cells?,?\s*cd4\+?|cd4\+?\s*t": "cd4-positive, alpha-beta t cell",
    r"t\s*cells?,?\s*cd8\+?|cd8\+?\s*t": "cd8-positive, alpha-beta t cell",
    # ... full pattern list
}
```

#### 4. Voting Module

```python
# ontology_mapper/voting.py
"""
Hierarchy-aware voting for multi-reference annotation.
Extract from ontology_utils.py lines 360-457
"""

from typing import List, Tuple, Optional

def hierarchy_aware_vote(
    cl_codes: List[str],
    scores: Optional[List[float]] = None,
    use_hierarchy: bool = True,
) -> Tuple[str, float]:
    """
    Vote on CL codes with hierarchy awareness.
    """
    # Implementation from ontology_utils.py
    ...
```

#### 5. Index Builder (Optional)

```python
# ontology_mapper/builder.py
"""
Build ontology index from OWL files.
Only needed if you want to rebuild the index.
Requires: pip install owlready2
"""

def build_ontology_index(output_path: str = "ontology_index.json") -> Dict:
    """
    Build complete ontology index from CL, NCIT, and UBERON.

    Downloads OWL files, extracts terms, saves to JSON.
    This is a ONE-TIME BUILD step.
    """
    # Implementation from ontology_utils.py lines 650-790
    ...
```

### Package Structure

```
ontology_mapper/
├── __init__.py
├── matcher.py           # Core search functionality
├── patterns.py          # CELL_TYPE_PATTERNS dict
├── voting.py            # hierarchy_aware_vote()
├── builder.py           # Index builder (optional)
├── data/
│   └── ontology_index.json  # Pre-built index (28.8 MB)
├── tests/
│   ├── test_matcher.py
│   ├── test_patterns.py
│   └── test_voting.py
└── setup.py
```

### setup.py

```python
from setuptools import setup, find_packages

setup(
    name="ontology-mapper",
    version="1.0.0",
    description="Map biological labels to ontology codes (CL, NCIT, UBERON)",
    packages=find_packages(),
    package_data={
        "ontology_mapper": ["data/ontology_index.json"],
    },
    install_requires=[],  # No dependencies for core functionality
    extras_require={
        "builder": ["owlready2>=0.40"],  # For rebuilding index
        "hierarchy": ["owlready2>=0.40"],  # For hierarchy operations
    },
    python_requires=">=3.8",
)
```

### Usage Example (spacebio)

```python
from spacebio.ontology import search_ontology_index, add_ontology_ids

# Search for CL terms
results = search_ontology_index(
    labels=["CD4+ T cells", "B cell", "Macrophage", "NK cells"],
    annotation_type="cell_type"
)

for label, matches in results.items():
    if matches:
        best = matches[0]
        print(f"{label:20} → {best['id']:15} {best['name']} (score: {best['score']:.2f})")

# Add ontology IDs to AnnData
adata, mappings = add_ontology_ids(
    adata,
    source_col="celltypist",
    target_col="cell_type_ontology_id",
    name_col="cell_type_ontology_name"
)
```

### Usage Example (Standalone)

```python
from ontology_mapper import search_ontology_index

# Basic usage
results = search_ontology_index(
    labels=["CD4+ T cells", "B cell", "Macrophage", "NK cells"],
    annotation_type="cell_type"
)

for label, matches in results.items():
    if matches:
        best = matches[0]
        print(f"{label:20} → {best['id']:15} {best['name']} (score: {best['score']:.2f})")
    else:
        print(f"{label:20} → No match")

# Output:
# CD4+ T cells         → CL:0000624      CD4-positive, alpha-beta T cell (score: 0.95)
# B cell               → CL:0000236      B cell (score: 1.00)
# Macrophage           → CL:0000235      macrophage (score: 1.00)
# NK cells             → CL:0000623      natural killer cell (score: 0.95)
```

---

## API Reference

### search_ontology_index()

```python
def search_ontology_index(
    labels: List[str],
    ontology_index: Dict[str, Dict[str, Dict[str, str]]] = None,
    index_uri: str = None,
    ontologies: List[str] = None,
    annotation_type: str = "cell_type",
    min_score: float = 0.7,
) -> Dict[str, List[Dict[str, Any]]]:
```

**Parameters**:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `labels` | `List[str]` | required | Labels to search |
| `ontology_index` | `Dict` | `None` | Pre-loaded index |
| `index_uri` | `str` | `None` | Path to index JSON |
| `ontologies` | `List[str]` | `None` | Override ontologies to search |
| `annotation_type` | `str` | `"cell_type"` | Type routing |
| `min_score` | `float` | `0.7` | Minimum match score |

**Returns**:

```python
{
    "label1": [
        {
            "id": "CL:0000236",
            "name": "B cell",
            "ontology": "cl",
            "score": 1.0,
            "match_type": "tier1_exact"
        },
        # ... more matches sorted by score
    ],
    "label2": [...],
}
```

### hierarchy_aware_vote()

```python
def hierarchy_aware_vote(
    cl_codes: List[str],
    scores: Optional[List[float]] = None,
    use_hierarchy: bool = True,
) -> Tuple[str, float]:
```

**Parameters**:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cl_codes` | `List[str]` | required | CL codes from references |
| `scores` | `List[float]` | `None` | Confidence scores per code |
| `use_hierarchy` | `bool` | `True` | Use hierarchy-aware voting |

**Returns**: `(winning_code, confidence)`

### load_ontology_index()

```python
def load_ontology_index(
    index_uri: str = None
) -> Dict[str, Dict[str, Dict[str, str]]]:
```

**Parameters**:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `index_uri` | `str` | `None` | Path to JSON index |

**Returns**: Nested dict `{ontology: {label: {id, name}}}`

---

## Extension Guide

### Adding New Cell Type Patterns

**When to Add**: Label consistently fails to match or matches wrong term.

**Steps**:

1. **Identify the problem**:
```python
results = search_ontology_index(["Failing Label"], annotation_type="cell_type")
print(results)  # Shows no match or wrong match
```

2. **Find correct CL term**:
```python
index = load_ontology_index()
for term, data in index['cl'].items():
    if 'keyword' in term:
        print(f"{data['id']:15} {data['name']}")
```

3. **Add pattern to CELL_TYPE_PATTERNS**:
```python
# In ontology_utils.py or patterns.py
CELL_TYPE_PATTERNS = {
    # ... existing patterns ...

    # Your new pattern
    r'your.*pattern|variant': 'canonical cl term',
}
```

4. **Test**:
```python
results = search_ontology_index(
    ["Failing Label", "Variant 1", "Variant 2"],
    annotation_type="cell_type"
)
for label, matches in results.items():
    if matches:
        print(f"✓ {label} → {matches[0]['id']}")
    else:
        print(f"✗ {label} → No match")
```

### Pattern Tips

| Pattern | Matches | Notes |
|---------|---------|-------|
| `r'\bnk\b'` | "NK" but not "PINK" | Word boundary |
| `r'cd4.*t'` | "CD4+ T", "CD4 T cells" | Flexible |
| `r'^b\s*cell'` | "B cell", "B cells" | Anchored |
| `r'th1\b\|type\s*1\s*helper'` | "Th1", "Type 1 helper" | Alternatives |

### Adding Manual Overrides

For one-off cases, use `manual_overrides` parameter:

```python
await convert_to_ontology_codes(
    data_uri="...",
    output_uri="...",
    source_column="cell_type",
    manual_overrides={
        "Weird Custom Label": "CL:0000236",
        "Unknown": None,  # Skip this label
        "CMS2": "NCIT:C4910",  # Pathology code
    }
)
```

### Rebuilding the Index

If ontology versions update:

```bash
# Install builder dependency
pip install owlready2

# Rebuild index
python -c "from ontology_utils import build_ontology_index; build_ontology_index('data/ontology_index.json')"
```

The builder downloads OWL files from:
- CL: `gs://spatial-bio-data/ontology/ontology_data/raw/cl.owl`
- NCIT: `gs://spatial-bio-data/ontology/ontology_data/raw/ncit.owl`
- UBERON: `gs://spatial-bio-data/ontology/ontology_data/raw/uberon.owl`

---

## Troubleshooting

### Common Issues

#### "My label matches but to the wrong term"

**Cause**: Token matching too permissive OR wrong ontology

**Solution 1**: Use more specific `annotation_type`
```python
annotation_type="cell_type"  # Uses CL only, not NCIT
```

**Solution 2**: Add specific pattern
```python
CELL_TYPE_PATTERNS[r'your_specific_label'] = 'exact_cl_term'
```

#### "Label with CD markers doesn't match"

**Cause**: Core words filtered out, no pattern exists

**Solution**: Add CD marker pattern
```python
CELL_TYPE_PATTERNS[r'cd25.*cd4'] = 'cd4-positive, cd25-positive, alpha-beta regulatory t cell'
```

#### "Too many NCIT matches for cell types"

**Cause**: Using `annotation_type="all"` instead of `"cell_type"`

**Solution**:
```python
search_ontology_index(labels, annotation_type="cell_type")  # CL only
```

#### "Generic term like 'T cells' matches too specifically"

**Expected behavior**: The system preserves specificity when available.

If you want generic "T cells" → CL:0000084 (generic T cell), that's correct. The pattern system only applies when more specific info is available.

#### "Index not found"

**Cause**: `ontology_index.json` not in expected location

**Solution**: Check paths:
```python
# Priority order checked:
# 1. data_for_docker_dep/ontology_index.json
# 2. data_for_dev/ontology_index.json
# 3. data/ontology_index.json
```

Or specify explicitly:
```python
load_ontology_index("path/to/ontology_index.json")
```

---

## Performance Benchmarks

### Match Quality (vs Previous Implementation)

| Metric | Before | After | Change |
|--------|--------|-------|--------|
| Overall match rate | 93% | 93%+ | Same |
| T cell specificity | Generic | Specific | **FIXED** |
| NCIT pollution | Yes | None | **FIXED** |
| Match quality | Mixed | High | **IMPROVED** |

### Speed Benchmarks

| Operation | Time | Notes |
|-----------|------|-------|
| Load index | ~2s | One-time per session |
| Search 1 label | <1ms | After index loaded |
| Search 100 labels | ~50ms | Batch operation |
| Search 1000 labels | ~500ms | Batch operation |
| Build index | ~90min | One-time, offline |

### Memory Usage

| Component | Size |
|-----------|------|
| Index JSON file | 28.8 MB |
| Index in memory | ~100 MB |
| OWL files (for build) | ~2 GB |

---

## Test Cases

Use these to verify any changes:

```python
from ontology_utils import search_ontology_index

test_cases = {
    # CD marker combinations
    "CD19+CD20+ B": "CL:0000236",           # B cell
    "CD4+ T cells": "CL:0000624",           # CD4-positive, alpha-beta T cell
    "CD8+ T cells": "CL:0000625",           # CD8-positive, alpha-beta T cell

    # Immunoglobulin subtypes (maps to SPECIFIC subtypes when available)
    "IgG+ Plasma": "CL:0000985",            # IgG plasma cell (specific)
    "IgA+ Plasma": "CL:0000987",            # IgA plasma cell (specific)

    # T cell subtypes (preserve specificity!)
    "T helper 17 cells": "CL:0000899",      # T-helper 17 cell
    "Regulatory T cells": "CL:0000815",     # regulatory T cell
    "Treg": "CL:0000815",                   # regulatory T cell

    # Tissue-specific
    "Mature Enterocytes type 2": "CL:0000584",  # enterocyte
    "Goblet cells": "CL:0000160",           # goblet cell

    # Functional states
    "Stalk-like ECs": "CL:0000115",         # endothelial cell
    "Tip-like ECs": "CL:0000115",           # endothelial cell

    # Abbreviations
    "NK cells": "CL:0000623",               # natural killer cell
    "DC": "CL:0000451",                     # dendritic cell
    "pDC": "CL:0000784",                    # plasmacytoid dendritic cell

    # Standard terms
    "Smooth muscle cells": "CL:0000192",    # smooth muscle cell
    "Mast cells": "CL:0000097",             # mast cell
    "Macrophage": "CL:0000235",             # macrophage
    "Fibroblast": "CL:0000057",             # fibroblast
}

# Run test
results = search_ontology_index(
    list(test_cases.keys()),
    annotation_type="cell_type"
)

# Verify
passed = 0
failed = 0
for label, expected_id in test_cases.items():
    matches = results.get(label, [])
    actual_id = matches[0]['id'] if matches else None
    if actual_id == expected_id:
        print(f"✓ {label:30} → {actual_id}")
        passed += 1
    else:
        print(f"✗ {label:30} → {actual_id} (expected: {expected_id})")
        failed += 1

print(f"\nResults: {passed}/{passed+failed} passed ({100*passed/(passed+failed):.1f}%)")
```

**Success criteria**: ≥ 90% correct matches with preserved specificity

---

## References

### Ontology Sources

- **Cell Ontology (CL)**: https://github.com/obophenotype/cell-ontology
- **NCI Thesaurus (NCIt)**: https://ncithesaurus.nci.nih.gov/
- **UBERON**: http://uberon.github.io/
- **OBO Foundry**: http://www.obofoundry.org/

### Tools Evaluated (Not Used)

- **CellOntologyMapper**: https://github.com/Starlitnightly/CellOntologyMapper
  - Pros: High accuracy (83-91%), cell-type specialized
  - Cons: 7.5 GB Docker bloat, requires internet for LLM

- **text2term**: https://github.com/rsgoncalves/text2term
  - Pros: Academic validation, well-tested, offline
  - Cons: General-purpose (not cell-type optimized)

- **CellO**: https://github.com/deweylab/CellO
  - Pros: Uses gene expression (more accurate)
  - Cons: Different use case (classification vs label mapping)

### Internal References

- R reference implementation: `test_scripts/ontology_converter_example.r`
- SingleR voting logic: `test_scripts/singleR_voting.r`

---

## Changelog

### 2025-12-16 (v2)
- Implemented full module in `spacebio/ontology.py`
- 135 cell type patterns (expanded from ~50)
- Fixed critical mappings: MAIT→mucosal invariant T cell, NK→natural killer cell, ILC3→group 3 ILC
- Added `add_ontology_ids()` for AnnData integration
- Added `validate_cl_term()`, `get_cl_id()`, `list_available_terms()` utilities
- Updated test cases: IgG/IgA plasma cells now map to specific CL IDs (CL:0000985, CL:0000987)
- Updated documentation with spacebio usage examples

### 2025-12-16 (v1)
- Added comprehensive standalone extraction guide
- Documented all tiers with code examples
- Added API reference section
- Added extension guide with patterns
- Added troubleshooting section
- Added performance benchmarks

### 2025-11-21
- Initial implementation
- Sequential tier matching (fixed from parallel)
- Annotation type routing
- CELL_TYPE_PATTERNS expanded to 100+ patterns
- Hierarchy-aware voting for multi-reference
