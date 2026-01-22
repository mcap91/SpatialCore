# Domain Detection

Identifying spatially contiguous tissue regions.

!!! note "Content Coming Soon"
    Full documentation will be added.

---

## Overview

**Domains** are spatially contiguous tissue regions with coherent expression patterns and/or cell-type composition. Unlike [niches](neighborhood_analysis.md) (which are compositional archetypes that can appear anywhere), domains partition the tissue into discrete, bounded spatial areas.

### Conceptual Workflow

| Step | Operation | Output |
|------|-----------|--------|
| **1. Neighborhoods** | Define local context (k-NN or radius) | Per-cell context features |
| **2. Niche Typing** | Cluster neighborhoods by composition | Niche labels (recurring archetypes) |
| **3. Domain Detection** | Spatially segment tissue (with smoothing) | Domain labels (contiguous regions) |

The key distinction: **Step 2** (niche typing) groups cells by *what kind* of microenvironment they're in (location-independent). **Step 3** (domain detection) identifies *where* in the tissue those cells are located (spatially contiguous regions).

---

## Example Domain Types

| Domain | Characteristics |
|--------|-----------------|
| **Tumor core** | High malignant cell density, low immune infiltration |
| **Tumor margin** | Mixed tumor/immune, active interface |
| **Immune aggregates** | Tertiary lymphoid structures, B/T cell zones |
| **Stromal regions** | Fibroblasts, vasculature, ECM-rich |

---

## Use Cases

- Tumor microenvironment characterization
- Mapping tissue architecture in healthy vs diseased tissue
- Identifying anatomical regions (cortex, medulla, etc.)
- Spatial segmentation for downstream analysis

---

## Related

- [Neighborhood & Niche Analysis](neighborhood_analysis.md) — Compositional clustering (location-independent)
- [Terminology](../index.md#terminology) — Canonical definitions with literature citations
