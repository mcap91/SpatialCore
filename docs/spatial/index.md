# Spatial Statistics

Spatial autocorrelation and bivariate correlation analysis for spatial transcriptomics.

---

## Overview

SpatialCore provides implementations of classical spatial statistics adapted for single-cell resolution spatial transcriptomics data, including both global (tissue-wide) and local (per-cell) variants.

## Methods

### Moran's I

Measures spatial autocorrelation of a single variable (gene expression).

| Variant | Description |
|---------|-------------|
| Global Moran's I | Tissue-wide autocorrelation statistic |
| Local Moran's I (LISA) | Per-cell hotspot/coldspot identification |

### Lee's L

Measures bivariate spatial correlation between two variables.

| Variant | Description |
|---------|-------------|
| Global Lee's L | Tissue-wide bivariate correlation |
| Local Lee's L | Per-cell HH/LL/HL/LH classification |
