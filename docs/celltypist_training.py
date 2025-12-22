"""
CellTypist custom model training for targeted gene panels.

This module provides functions to train custom CellTypist models using CellxGene
reference data, specifically for targeted gene panels like Xenium (~400 genes).

The problem: Pre-trained CellTypist models are trained on full transcriptomes
(15-30k genes), resulting in only 10-25% gene overlap with targeted panels.
This causes unreliable predictions and many "Unknown" assignments.

The solution: Train custom models on CellxGene atlas data subsetted to exactly
the genes available in the target panel.

Functions:
- prepare_celltypist_model: Composite orchestrator for model preparation
- validate_reference_data: Validate CellxGene reference data
- subset_reference_to_panel: Subset reference to panel genes
- train_celltypist_model: Train CellTypist model on subset
- list_registered_models: Query the model registry

Model Registry:
Models are stored in GCS at gs://spatial-bio-output/models/celltypist/
with versioning and metadata for reuse across samples.

Example Usage:
    # Prepare custom model for Xenium colon panel
    result = await prepare_celltypist_model(
        reference_uri="file:///path/to/cellxgene_colon_atlas.h5ad",
        output_uri="file:///output/metadata.json",
        xenium_data_uri="file:///path/to/xenium_sample.h5ad",
        tissue_name="colon",
        model_version="v1",
    )

    # Use the trained model for annotation
    await annotate_celltypist(
        data_uri=xenium_data_uri,
        output_uri=annotation_output,
        model_uri=result["model_uri"],
    )
"""

import json
import logging
import os
from datetime import datetime
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
from temporalio import activity

from universal_io import process_activity, universal_write
from gcs_utils import download_from_gcs, upload_to_gcs

logger = logging.getLogger(__name__)

# Default registry location
DEFAULT_REGISTRY_PATH = os.getenv(
    "CELLTYPIST_REGISTRY_PATH", "gs://spatial-bio-output/models/celltypist"
)

# Default max cells per reference for memory safety
# CellTypist uses mini-batch SGD (1000 cells × 100 batches), so training handles
# large datasets well. The constraint is loading h5ad files before gene subsetting.
# - 200k: Safe for 16GB workers
# - 500k: Good for 32GB workers (default)
# - None: No limit for 64GB+ workers
# Cell types need ≥100 cells for good predictions (CellTypist docs)
DEFAULT_MAX_CELLS_PER_REF = int(os.getenv("CELLTYPIST_MAX_CELLS_PER_REF", "500000"))


# ============================================================================
# UTILITY FUNCTIONS: Normalization detection, subsampling, label harmonization
# ============================================================================


def check_normalization_status(adata: sc.AnnData, sample_size: int = 1000) -> dict:
    """
    Check if data is normalized and what type of normalization was applied.

    Detects:
    - Raw counts (integers, high values)
    - Log-normalized (log1p to 10k, typical mean 0-3)
    - Already normalized but not logged

    Args:
        adata: AnnData object to check
        sample_size: Number of cells to sample for speed

    Returns:
        Dict with keys:
        - is_normalized: bool
        - is_log_transformed: bool
        - normalization_type: str ("raw", "log1p_10k", "normalized_not_logged", "unknown")
        - stats: dict with mean, max, min, fraction_integer
    """
    # Sample data for speed
    n_sample = min(sample_size, adata.n_obs)
    if issparse(adata.X):
        sample_data = adata.X[:n_sample].toarray()
    else:
        sample_data = adata.X[:n_sample]

    # Compute statistics
    mean_val = float(np.mean(sample_data))
    max_val = float(np.max(sample_data))
    min_val = float(np.min(sample_data))

    # Check for integer counts (raw data indicator)
    # Sample a subset to check if values are integers
    flat_sample = sample_data.flatten()[:10000]
    non_zero = flat_sample[flat_sample != 0]
    if len(non_zero) > 0:
        fraction_integer = np.mean(np.equal(np.mod(non_zero, 1), 0))
    else:
        fraction_integer = 1.0

    stats = {
        "mean": round(mean_val, 4),
        "max": round(max_val, 4),
        "min": round(min_val, 4),
        "fraction_integer": round(fraction_integer, 4),
    }

    # Decision logic
    # Raw counts: mostly integers with high max (sparse data may have low mean)
    # Key indicator: fraction_integer > 0.9 AND max > 100 (much higher than log range)
    if fraction_integer > 0.9 and max_val > 100:
        return {
            "is_normalized": False,
            "is_log_transformed": False,
            "normalization_type": "raw",
            "stats": stats,
        }

    # Log-normalized: mean typically 0-3, max typically 5-15, non-negative
    if min_val >= 0 and mean_val < 5 and max_val < 20:
        return {
            "is_normalized": True,
            "is_log_transformed": True,
            "normalization_type": "log1p_10k",
            "stats": stats,
        }

    # Normalized but maybe not logged: low mean but high max, or negative values
    if mean_val < 10 and (max_val > 20 or min_val < 0):
        return {
            "is_normalized": True,
            "is_log_transformed": False,
            "normalization_type": "normalized_not_logged",
            "stats": stats,
        }

    # Unknown state
    return {
        "is_normalized": False,
        "is_log_transformed": False,
        "normalization_type": "unknown",
        "stats": stats,
    }


def ensure_normalized(
    adata: sc.AnnData,
    target_sum: float = 1e4,
    force_renormalize: bool = False,
) -> sc.AnnData:
    """
    Ensure data is log-normalized (log1p to 10k counts).

    Checks current normalization status and applies normalization if needed.
    CellTypist requires log-normalized data for proper model training.

    Args:
        adata: AnnData object (modified in place)
        target_sum: Target library size (default 1e4)
        force_renormalize: If True, renormalize even if already normalized

    Returns:
        AnnData with log-normalized expression in X
    """
    status = check_normalization_status(adata)

    print(f"  Normalization status: {status['normalization_type']}")
    print(f"    Stats: mean={status['stats']['mean']}, max={status['stats']['max']}")

    if status["normalization_type"] == "log1p_10k" and not force_renormalize:
        print("  ✓ Data already log-normalized, skipping")
        return adata

    if status["normalization_type"] == "raw" or force_renormalize:
        print(f"  Normalizing: total counts to {target_sum:.0f}, then log1p...")
        # exclude_highly_expressed=False required for CellTypist compatibility
        sc.pp.normalize_total(
            adata, target_sum=target_sum, exclude_highly_expressed=False
        )
        sc.pp.log1p(adata)
        print("  ✓ Normalization complete")
        return adata

    if status["normalization_type"] == "normalized_not_logged":
        print("  Data normalized but not log-transformed, applying log1p...")
        sc.pp.log1p(adata)
        print("  ✓ Log transformation complete")
        return adata

    # Unknown state - raise error rather than silently normalizing
    raise ValueError(
        f"Unknown normalization state: {status['normalization_type']}. "
        f"Stats: mean={status['stats']['mean']}, max={status['stats']['max']}, min={status['stats']['min']}. "
        "Please provide raw counts or log-normalized data."
    )


def get_log_normalized_data(adata: sc.AnnData) -> sc.AnnData:
    """
    Extract log-normalized expression with priority-based data selection.

    CellTypist requires log-normalized data for training. This function
    checks multiple locations for the appropriate data:

    Priority order:
    1. adata.layers["log"] - if available, use directly
    2. adata.raw - if X appears scaled (z-score), use raw
    3. adata.X - normalize if needed

    Args:
        adata: AnnData object (may have scaled X, raw, or layers)

    Returns:
        AnnData with log-normalized expression in X
    """
    # Priority 1: Check for "log" layer
    if adata.layers is not None and "log" in adata.layers:
        print("  Found adata.layers['log'] - using log-normalized layer")
        result = adata.copy()
        result.X = result.layers["log"]
        return result

    # Priority 2: Check normalization status
    status = check_normalization_status(adata)

    # If already log-normalized, return as-is
    if status["normalization_type"] == "log1p_10k":
        print("  ✓ adata.X already log-normalized")
        return adata

    # If raw counts detected, normalize
    if status["normalization_type"] == "raw":
        print("  adata.X is raw counts, normalizing...")
        return ensure_normalized(adata)

    # Scaled data indicators: mean ~0, negative values (but NOT raw sparse data)
    is_scaled = status["stats"]["min"] < -1 or (
        abs(status["stats"]["mean"]) < 0.5 and status["stats"]["max"] < 20
    )

    if is_scaled and adata.raw is not None:
        print("  Data appears scaled (z-score), checking adata.raw...")
        raw_adata = adata.raw.to_adata()
        raw_status = check_normalization_status(raw_adata)

        if raw_status["normalization_type"] == "log1p_10k":
            print("  ✓ Using adata.raw (log-normalized)")
            # Preserve obs from original (may have more columns)
            raw_adata.obs = adata.obs.loc[raw_adata.obs_names]
            return raw_adata
        else:
            print(f"  adata.raw is {raw_status['normalization_type']}, normalizing...")
            return ensure_normalized(raw_adata)

    if is_scaled:
        raise ValueError(
            "Data appears scaled (z-scored) but no adata.raw available. "
            "CellTypist requires log1p normalized data (not scaled). "
            "Please provide data with raw counts or log-normalized expression."
        )

    # Unknown state - try to normalize
    print(f"  adata.X is {status['normalization_type']}, normalizing...")
    return ensure_normalized(adata)


def subsample_adata(
    adata: sc.AnnData,
    max_cells: int,
    stratify_by: Optional[str] = None,
    random_state: int = 42,
) -> sc.AnnData:
    """
    Subsample AnnData to max_cells for memory management.

    Optionally stratifies by a column to maintain cell type proportions.

    Args:
        adata: AnnData to subsample
        max_cells: Maximum number of cells to keep
        stratify_by: Column to stratify by (e.g., 'cell_type')
        random_state: Random seed for reproducibility

    Returns:
        Subsampled AnnData
    """
    if adata.n_obs <= max_cells:
        return adata

    np.random.seed(random_state)

    if stratify_by and stratify_by in adata.obs.columns:
        # Stratified sampling to maintain proportions
        indices = []
        groups = adata.obs[stratify_by].unique()
        n_per_group = max_cells // len(groups)

        for group in groups:
            group_idx = np.where(adata.obs[stratify_by] == group)[0]
            n_sample = min(len(group_idx), n_per_group)
            sampled = np.random.choice(group_idx, size=n_sample, replace=False)
            indices.extend(sampled)

        # If we have room for more, fill with random cells
        remaining = max_cells - len(indices)
        if remaining > 0:
            all_idx = set(range(adata.n_obs))
            available = list(all_idx - set(indices))
            if available:
                extra = np.random.choice(
                    available, size=min(remaining, len(available)), replace=False
                )
                indices.extend(extra)

        indices = np.array(indices)
    else:
        # Random sampling
        indices = np.random.choice(adata.n_obs, size=max_cells, replace=False)

    print(f"  Subsampled {adata.n_obs:,} → {len(indices):,} cells")
    return adata[indices].copy()


def harmonize_cell_type_labels(
    adata: sc.AnnData,
    label_column: str,
    output_column: str = "unified_cell_type",
    use_ontology: bool = True,
    min_match_score: float = 0.7,
) -> sc.AnnData:
    """
    Harmonize cell type labels using Cell Ontology.

    Maps varied labels to canonical CL-based names for consistency
    across multiple reference datasets:
    - "B cell", "B cells", "B_cell" → "B cell" (via CL:0000236)
    - "CD4+ T", "CD4-positive T cell" → "CD4-positive, alpha-beta T cell"

    Args:
        adata: AnnData with cell type labels
        label_column: Column containing cell type labels
        output_column: Column to store harmonized labels
        use_ontology: If True, use ontology for harmonization
        min_match_score: Minimum score for ontology match (0-1)

    Returns:
        AnnData with harmonized labels in output_column
    """
    if label_column not in adata.obs.columns:
        raise ValueError(f"Label column '{label_column}' not found in adata.obs")

    unique_labels = adata.obs[label_column].dropna().unique().tolist()
    print(f"  Harmonizing {len(unique_labels)} unique cell type labels...")

    label_mapping = {}

    if use_ontology:
        try:
            from ontology_utils import search_ontology_index, get_term_name

            # Search all labels at once
            results = search_ontology_index(
                unique_labels,
                annotation_type="cell_type",
                min_score=min_match_score,
            )

            matched = 0
            for label in unique_labels:
                matches = results.get(label, [])
                if matches and matches[0]["score"] >= min_match_score:
                    # Use the canonical name from ontology
                    cl_code = matches[0]["id"]
                    canonical = matches[0]["name"]
                    # Try to get official term name if available
                    try:
                        official_name = get_term_name(cl_code)
                        if official_name:
                            canonical = official_name
                    except Exception:
                        pass
                    label_mapping[label] = canonical
                    matched += 1
                else:
                    # Keep original if no good match
                    label_mapping[label] = label

            print(f"  ✓ Matched {matched}/{len(unique_labels)} labels to Cell Ontology")

        except ImportError:
            print("  ⚠️  ontology_utils not available, using simple normalization")
            use_ontology = False

    if not use_ontology:
        # Simple normalization: lowercase, strip, replace underscores
        for label in unique_labels:
            normalized = str(label).strip()
            normalized = normalized.replace("_", " ")
            # Capitalize first letter of each word
            normalized = " ".join(word.capitalize() for word in normalized.split())
            label_mapping[label] = normalized

    # Apply mapping
    adata.obs[output_column] = adata.obs[label_column].map(label_mapping)

    # Handle any NaN values - convert to string first to avoid Categorical issues
    if adata.obs[output_column].isna().any():
        adata.obs[output_column] = (
            adata.obs[output_column].astype(str).replace("nan", "Unknown")
        )

    # Show mapping summary
    n_unique_output = adata.obs[output_column].nunique()
    print(f"  Unified to {n_unique_output} harmonized cell types")

    return adata


def combine_references(
    reference_paths: List[str],
    label_columns: List[str],
    output_column: str = "unified_cell_type",
    max_cells_per_ref: Optional[int] = None,
    target_genes: Optional[List[str]] = None,
    harmonize_labels: bool = True,
    normalize_data: bool = True,
) -> sc.AnnData:
    """
    Combine multiple reference datasets for CellTypist training.

    Handles:
    1. Loading references from local paths or GCS
    2. Subsampling large references for memory management
    3. Checking/applying normalization consistently
    4. Harmonizing cell type labels across datasets
    5. Concatenating with shared gene space

    Args:
        reference_paths: List of paths to h5ad files (local or GCS)
        label_columns: Cell type column name for each reference
        output_column: Unified cell type column name
        max_cells_per_ref: Max cells to sample per reference (None = no limit)
        target_genes: Optional list of genes to subset to
        harmonize_labels: Use ontology to harmonize labels
        normalize_data: Ensure all data is log-normalized

    Returns:
        Combined AnnData with unified cell types

    Raises:
        MemoryError: If combined data is too large (caught and reported)
        ValueError: If no shared genes found
    """
    import gc
    import psutil

    if len(reference_paths) != len(label_columns):
        raise ValueError("Must provide one label_column per reference_path")

    if max_cells_per_ref is None:
        max_cells_per_ref = DEFAULT_MAX_CELLS_PER_REF

    print(f"\n{'='*60}", flush=True)
    print(f"Combining {len(reference_paths)} reference datasets", flush=True)
    print(f"{'='*60}", flush=True)
    print(f"Max cells per reference: {max_cells_per_ref:,}", flush=True)

    adatas = []

    for i, (path, label_col) in enumerate(zip(reference_paths, label_columns)):
        print(
            f"\n[{i+1}/{len(reference_paths)}] Loading: {Path(path).name}", flush=True
        )

        # Check available memory
        mem = psutil.virtual_memory()
        print(f"  Available memory: {mem.available / 1e9:.1f} GB", flush=True)

        # Load reference with comprehensive error handling
        try:
            if path.startswith("gs://"):
                local_path = Path(f"/tmp/ref_{i}.h5ad")
                try:
                    download_from_gcs(path, local_path)
                    adata = sc.read_h5ad(local_path)
                finally:
                    # Always clean up temp file
                    if local_path.exists():
                        local_path.unlink()
            else:
                local_path = (
                    path.replace("file://", "") if path.startswith("file://") else path
                )
                adata = sc.read_h5ad(local_path)

            print(
                f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes", flush=True
            )

            # Fix gene names if they're numeric indices (CellxGene Census format)
            # Gene symbols are often stored in var['feature_name'] column
            if adata.var_names[0].isdigit() and "feature_name" in adata.var.columns:
                print(
                    "  ⚠️  Detected numeric var_names, recovering from 'feature_name' column",
                    flush=True,
                )
                adata.var_names = adata.var["feature_name"].values
                # Handle duplicates by making unique
                adata.var_names_make_unique()

            # Subsample if needed
            if adata.n_obs > max_cells_per_ref:
                adata = subsample_adata(
                    adata,
                    max_cells=max_cells_per_ref,
                    stratify_by=label_col if label_col in adata.obs.columns else None,
                )

            # Check label column exists
            if label_col not in adata.obs.columns:
                # Try common alternatives
                alternatives = [
                    "cell_type",
                    "celltype",
                    "CellType",
                    "cell_ontology_class",
                ]
                found = None
                for alt in alternatives:
                    if alt in adata.obs.columns:
                        found = alt
                        break
                if found:
                    print(
                        f"  ⚠️  Column '{label_col}' not found, using '{found}'",
                        flush=True,
                    )
                    label_col = found
                else:
                    raise ValueError(
                        f"Label column '{label_col}' not found. Available: {list(adata.obs.columns)[:10]}"
                    )

            # Get log-normalized data (checks layers["log"], raw, then X)
            if normalize_data:
                adata = get_log_normalized_data(adata)

            # Harmonize labels
            if harmonize_labels:
                adata = harmonize_cell_type_labels(
                    adata,
                    label_column=label_col,
                    output_column=output_column,
                )
            else:
                adata.obs[output_column] = adata.obs[label_col]

            # Add source tracking
            adata.obs["reference_source"] = Path(path).stem

            # Clean up to save memory
            # Keep only essential columns
            essential_cols = [output_column, "reference_source", label_col]
            cols_to_drop = [c for c in adata.obs.columns if c not in essential_cols]
            adata.obs = adata.obs.drop(columns=cols_to_drop, errors="ignore")

            adatas.append(adata)

        except MemoryError:
            print(
                f"  ❌ MemoryError loading reference {i}: {Path(path).name}", flush=True
            )
            print(
                "     Try reducing max_cells_per_ref or using fewer references",
                flush=True,
            )
            raise
        except Exception as e:
            print(
                f"  ❌ Error loading reference {i} ({Path(path).name}): {e}", flush=True
            )
            import traceback

            traceback.print_exc()
            raise RuntimeError(f"Failed to load reference {path}") from e

        gc.collect()

    # Subset to target genes if provided
    if target_genes:
        print(f"\nSubsetting to {len(target_genes)} target genes...", flush=True)
        target_set = set(target_genes)
        for i, adata in enumerate(adatas):
            overlap = set(adata.var_names) & target_set
            if len(overlap) < len(target_set) * 0.1:
                print(
                    f"  ⚠️  Reference {i} has only {len(overlap)} of {len(target_set)} target genes",
                    flush=True,
                )
            adatas[i] = adata[:, list(overlap)].copy()
            print(f"  Reference {i}: {adatas[i].n_vars} genes after subset", flush=True)

    # Find shared genes (inner join)
    print("\nFinding shared genes across all references...", flush=True)
    shared_genes = set(adatas[0].var_names)
    for adata in adatas[1:]:
        shared_genes &= set(adata.var_names)

    if len(shared_genes) == 0:
        raise ValueError("No shared genes found across references!")

    print(f"  Shared genes: {len(shared_genes)}", flush=True)

    # Subset all to shared genes
    for i in range(len(adatas)):
        adatas[i] = adatas[i][:, list(shared_genes)].copy()

    # Check memory before concatenation
    mem = psutil.virtual_memory()
    total_cells = sum(a.n_obs for a in adatas)
    total_genes = adatas[0].n_vars if adatas else 0
    estimated_mb = (total_cells * total_genes * 8 * 2) / 1e6  # 2x for concat overhead
    print("\nMemory check before concat:", flush=True)
    print(f"  Available: {mem.available/1e9:.1f} GB", flush=True)
    print(f"  Estimated need: ~{estimated_mb/1e3:.1f} GB", flush=True)

    # Concatenate
    print("\nConcatenating references...", flush=True)
    try:
        combined = sc.concat(adatas, join="inner", label="batch", index_unique="-")
        print(
            f"  ✓ Combined: {combined.n_obs:,} cells × {combined.n_vars:,} genes",
            flush=True,
        )

        # Show cell type distribution
        print("\n  Cell type distribution:", flush=True)
        for ct, count in combined.obs[output_column].value_counts().head(10).items():
            print(f"    {ct}: {count:,} cells", flush=True)
        if combined.obs[output_column].nunique() > 10:
            print(
                f"    ... and {combined.obs[output_column].nunique() - 10} more types",
                flush=True,
            )

    except MemoryError:
        print("  ❌ MemoryError during concatenation", flush=True)
        print(f"     Total cells: {sum(a.n_obs for a in adatas):,}", flush=True)
        print("     Try reducing max_cells_per_ref", flush=True)
        raise

    # Final cleanup
    del adatas
    gc.collect()

    return combined


# ============================================================================
# ATOMIC FUNCTION: combine_reference_datasets (Temporal activity wrapper)
# ============================================================================


@activity.defn
async def combine_reference_datasets(
    reference_uris: List[str],
    output_uri: str,
    label_columns: Optional[List[str]] = None,
    xenium_data_uri: Optional[str] = None,
    max_cells_per_ref: int = 100000,
    harmonize_labels: bool = True,
    normalize_data: bool = True,
) -> str:
    """
    Combine multiple reference datasets for CellTypist training.

    Temporal activity wrapper for combine_references(). Loads multiple
    reference h5ad files, harmonizes labels using Cell Ontology, ensures
    consistent normalization, and outputs a combined reference suitable
    for training a comprehensive CellTypist model.

    Features:
    - Automatic normalization detection and log1p normalization if needed
    - Stratified subsampling to prevent memory issues with large references
    - Cell type label harmonization via Cell Ontology
    - Gene intersection with optional subsetting to Xenium panel genes

    Args:
        reference_uris: List of GCS/local URIs to reference h5ad files
        output_uri: GCS URI for combined reference output
        label_columns: Cell type column name for each reference (auto-detected if None)
        xenium_data_uri: Optional URI to Xenium data for panel gene subsetting
        max_cells_per_ref: Max cells to sample per reference (default 100k)
        harmonize_labels: Use ontology to harmonize labels (default True)
        normalize_data: Ensure all data is log-normalized (default True)

    Returns:
        GCS URI to combined reference h5ad file

    Example:
        # Combine epithelial and immune references for CRC
        result = await combine_reference_datasets(
            reference_uris=[
                "gs://spatial-bio-output/references/cellxgene/crc_progressive_plasticity/Epithelial.h5ad",
                "gs://spatial-bio-output/references/cellxgene/colon_immune_niches/colon_immune_niches.h5ad",
            ],
            output_uri="gs://spatial-bio-output/references/combined/crc_comprehensive.h5ad",
            label_columns=["celltype", "cell_type"],
            xenium_data_uri="gs://my-bucket/xenium_sample.h5ad",  # Optional: subset to panel genes
            max_cells_per_ref=50000,  # Keep 50k cells per reference
        )

    @metadata
    is_composite: false
    ready_for_production: false
    leads_to: [prepare_celltypist_model, train_celltypist_model]
    follows_from: []
    category: preprocessing
    runtime: "10-60min"
    memory: "32GB"
    algorithm_type: data_preparation
    method_type: reference_combination
    prerequisites:
      - reference_h5ads_available: true
      - cell_type_labels_required: true
    output_guarantees:
      - combined_h5ad: true
      - unified_cell_type_column: true
      - log_normalized_data: true
    best_practices: |
      Combines multiple CellxGene reference datasets for comprehensive model training.

      When to use:
      - Training models that need both epithelial and immune cell coverage
      - Combining tissue-specific references from different studies
      - Creating comprehensive CRC models from Progressive Plasticity + immune data

      Memory management:
      - Set max_cells_per_ref to limit memory usage
      - Default 100k cells per reference = ~400k total for 4 references
      - Uses stratified sampling to maintain cell type proportions

      Label harmonization:
      - Enabled by default, uses Cell Ontology for standardization
      - Maps varied labels like "B cell", "B cells", "B_cell" → "B cell"
      - Creates "unified_cell_type" column for training
    workflow_position: preprocessing
    recommended_sequence: [combine_reference_datasets, prepare_celltypist_model, annotate_celltypist]
    suppress_from_mcp: true
    """
    import gc
    from pathlib import Path

    print(f"\n{'='*80}")
    print("combine_reference_datasets - Multi-Reference Combination")
    print(f"{'='*80}")
    print(f"References: {len(reference_uris)}")
    print(f"Max cells per ref: {max_cells_per_ref:,}")

    # Auto-detect label columns if not provided
    if label_columns is None:
        label_columns = ["cell_type"] * len(reference_uris)
        print("Using default label column: 'cell_type'")

    if len(label_columns) != len(reference_uris):
        raise ValueError(
            f"Must provide one label_column per reference. "
            f"Got {len(label_columns)} columns for {len(reference_uris)} references."
        )

    # Get panel genes if xenium_data_uri provided
    target_genes = None
    if xenium_data_uri:
        print(f"\nExtracting panel genes from: {xenium_data_uri}")
        temp_path = Path("/tmp/xenium_panel.h5ad")
        if xenium_data_uri.startswith("gs://"):
            download_from_gcs(xenium_data_uri, temp_path)
        else:
            import shutil

            local_path = xenium_data_uri.replace("file://", "")
            shutil.copy(local_path, temp_path)

        xenium_adata = sc.read_h5ad(temp_path)
        target_genes = list(xenium_adata.var_names)
        print(f"  Panel genes: {len(target_genes)}")
        del xenium_adata
        temp_path.unlink()
        gc.collect()

    # Combine references
    combined = combine_references(
        reference_paths=reference_uris,
        label_columns=label_columns,
        output_column="unified_cell_type",
        max_cells_per_ref=max_cells_per_ref,
        target_genes=target_genes,
        harmonize_labels=harmonize_labels,
        normalize_data=normalize_data,
    )

    # Save combined reference
    print("\nSaving combined reference...")
    output_path = Path("/tmp/combined_reference.h5ad")
    combined.write_h5ad(output_path)

    if output_uri.startswith("gs://"):
        upload_to_gcs(output_path, output_uri)
        output_path.unlink()
    else:
        import shutil

        local_out = output_uri.replace("file://", "")
        Path(local_out).parent.mkdir(parents=True, exist_ok=True)
        shutil.move(output_path, local_out)

    print(f"  ✓ Saved to: {output_uri}")

    # Create metadata
    result = {
        "output_uri": output_uri,
        "n_cells": combined.n_obs,
        "n_genes": combined.n_vars,
        "n_cell_types": combined.obs["unified_cell_type"].nunique(),
        "cell_types": combined.obs["unified_cell_type"]
        .value_counts()
        .head(20)
        .to_dict(),
        "references_combined": len(reference_uris),
        "max_cells_per_ref": max_cells_per_ref,
        "harmonized_labels": harmonize_labels,
        "panel_genes_used": len(target_genes) if target_genes else None,
    }

    # Write metadata
    metadata_uri = output_uri.replace(".h5ad", "_metadata.json")
    universal_write(metadata_uri, result)

    del combined
    gc.collect()

    return output_uri


# ============================================================================
# ATOMIC FUNCTION: validate_reference_data
# ============================================================================


@activity.defn
async def validate_reference_data(
    data_uri: str,
    output_uri: str,
    label_column: str = "cell_type",
    min_cells_per_type: int = 50,
    check_normalization: bool = True,
    expected_normalization: str = "log1p_10k",
) -> str:
    """
    Validate CellxGene reference data before custom model training.

    Performs comprehensive validation of reference data to ensure it's suitable
    for CellTypist model training. Checks data integrity, normalization status,
    cell type distribution, and available annotation columns.

    This is a critical first step before training custom models - catching
    data issues early prevents wasted compute and failed training runs.

    Args:
        data_uri: GCS/local URI to CellxGene h5ad file
        output_uri: GCS URI for validation report
        label_column: Column containing cell type labels (default: "cell_type")
        min_cells_per_type: Minimum cells required per cell type (default: 50)
        check_normalization: Whether to validate normalization (default: True)
        expected_normalization: Expected normalization type (default: "log1p_10k")

    Returns:
        GCS URI to validation report with detailed diagnostics

    @metadata
    is_composite: false
    ready_for_production: false
    leads_to: [subset_reference_to_panel, train_celltypist_model]
    follows_from: []
    category: utility
    runtime: "2-5min"
    memory: "8GB"
    algorithm_type: data_validation
    method_type: reference_validation
    best_for: [cellxgene_data, atlas_data, reference_preparation]
    prerequisites:
      - h5ad_file_required: true
      - cell_type_labels_required: true
    output_guarantees:
      - validation_report: true
      - normalization_check: true
      - cell_type_distribution: true
    best_practices: |
      Validates reference data before CellTypist model training:
      1. Checks data can be loaded (valid h5ad format)
      2. Verifies cell type label column exists
      3. Checks cell count per type (warns if below threshold)
      4. Validates normalization (CellTypist expects log1p to 10k)
      5. Lists available label columns for selection

      When to use:
      - Before any custom model training
      - When using new CellxGene reference datasets
      - To diagnose training failures

      When NOT to use:
      - For already-validated reference data
      - If validation was performed in a previous run
    workflow_position: utility
    recommended_sequence: [validate_reference_data, subset_reference_to_panel, train_celltypist_model]
    suppress_from_mcp: true
    """
    params = {
        "data_uri": data_uri,
        "output_uri": output_uri,
        "label_column": label_column,
        "min_cells_per_type": min_cells_per_type,
        "check_normalization": check_normalization,
        "expected_normalization": expected_normalization,
    }

    async def process_fn(input_data: dict, params: dict) -> dict:
        """Validate reference data for CellTypist training."""
        print(f"\n{'='*80}")
        print("validate_reference_data - Validating CellxGene Reference")
        print(f"{'='*80}")

        input_path = Path("./input/data.h5ad")
        label_column = params["label_column"]
        min_cells = params["min_cells_per_type"]

        validation_result = {
            "valid": True,
            "warnings": [],
            "errors": [],
            "n_cells": 0,
            "n_genes": 0,
            "n_cell_types": 0,
            "cell_type_counts": {},
            "cell_types_below_threshold": [],
            "normalization_check": {},
            "available_label_columns": [],
        }

        # Load data
        print("\n[1/4] Loading reference data...")
        try:
            adata = sc.read_h5ad(input_path)
            validation_result["n_cells"] = adata.n_obs
            validation_result["n_genes"] = adata.n_vars
            print(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
        except Exception as e:
            validation_result["valid"] = False
            validation_result["errors"].append(f"Failed to load h5ad: {str(e)}")
            return validation_result

        # Check available label columns
        print("\n[2/4] Checking available annotation columns...")
        obs_columns = list(adata.obs.columns)
        validation_result["available_label_columns"] = obs_columns

        # Look for common cell type columns
        common_cell_type_cols = [
            "cell_type",
            "cell_ontology_class",
            "celltype",
            "CellType",
            "cell_type_ontology_term_id",
            "author_cell_type",
        ]
        found_cols = [c for c in common_cell_type_cols if c in obs_columns]
        print(f"  Found cell type columns: {found_cols}")
        print(
            f"  All obs columns: {obs_columns[:20]}{'...' if len(obs_columns) > 20 else ''}"
        )

        # Check requested label column
        if label_column not in adata.obs.columns:
            validation_result["valid"] = False
            validation_result["errors"].append(
                f"Label column '{label_column}' not found. "
                f"Available: {found_cols or obs_columns[:10]}"
            )
            return validation_result

        # Cell type distribution
        print(f"\n[3/4] Analyzing cell type distribution (column: {label_column})...")
        cell_type_counts = adata.obs[label_column].value_counts().to_dict()
        validation_result["cell_type_counts"] = {
            str(k): int(v) for k, v in cell_type_counts.items()
        }
        validation_result["n_cell_types"] = len(cell_type_counts)

        print(f"  Found {len(cell_type_counts)} unique cell types")

        # Check for types below threshold
        below_threshold = [
            ct for ct, count in cell_type_counts.items() if count < min_cells
        ]
        validation_result["cell_types_below_threshold"] = below_threshold

        if below_threshold:
            validation_result["warnings"].append(
                f"{len(below_threshold)} cell types have fewer than {min_cells} cells"
            )
            print(
                f"  ⚠️  {len(below_threshold)} cell types below threshold ({min_cells} cells):"
            )
            for ct in below_threshold[:5]:
                print(f"      - {ct}: {cell_type_counts[ct]} cells")
            if len(below_threshold) > 5:
                print(f"      ... and {len(below_threshold) - 5} more")

        # Show top cell types
        sorted_types = sorted(cell_type_counts.items(), key=lambda x: -x[1])
        print("\n  Top 10 cell types:")
        for ct, count in sorted_types[:10]:
            print(f"      {ct}: {count:,} cells")

        # Normalization check
        print("\n[4/4] Checking normalization...")
        if params["check_normalization"]:
            # Sample some values
            if issparse(adata.X):
                sample_data = (
                    adata.X[:1000].toarray()
                    if adata.n_obs > 1000
                    else adata.X.toarray()
                )
            else:
                sample_data = adata.X[:1000] if adata.n_obs > 1000 else adata.X

            # Check characteristics
            mean_expr = float(np.mean(sample_data))
            max_expr = float(np.max(sample_data))
            min_expr = float(np.min(sample_data))

            # CellTypist expects log1p normalized data (typically mean 0-2, max 5-15)
            appears_normalized = mean_expr < 10 and max_expr < 50 and min_expr >= 0

            validation_result["normalization_check"] = {
                "appears_normalized": appears_normalized,
                "mean_expression": round(mean_expr, 4),
                "max_expression": round(max_expr, 4),
                "min_expression": round(min_expr, 4),
            }

            if appears_normalized:
                print("  ✓ Data appears log-normalized")
                print(
                    f"    Mean: {mean_expr:.4f}, Max: {max_expr:.4f}, Min: {min_expr:.4f}"
                )
            else:
                validation_result["warnings"].append(
                    f"Data may not be log-normalized (mean={mean_expr:.2f}, max={max_expr:.2f})"
                )
                print("  ⚠️  Data may not be properly normalized")
                print(f"    Mean: {mean_expr:.4f}, Max: {max_expr:.4f}")
                print("    CellTypist expects log1p(CPM to 10k) normalized data")

        # Summary
        print(f"\n{'='*80}")
        print("VALIDATION SUMMARY")
        print(f"{'='*80}")
        print(f"  Valid: {validation_result['valid']}")
        print(f"  Cells: {validation_result['n_cells']:,}")
        print(f"  Genes: {validation_result['n_genes']:,}")
        print(f"  Cell types: {validation_result['n_cell_types']}")
        print(f"  Warnings: {len(validation_result['warnings'])}")
        print(f"  Errors: {len(validation_result['errors'])}")

        if validation_result["errors"]:
            print("\nERRORS:")
            for err in validation_result["errors"]:
                print(f"  ❌ {err}")

        if validation_result["warnings"]:
            print("\nWARNINGS:")
            for warn in validation_result["warnings"]:
                print(f"  ⚠️  {warn}")

        return validation_result

    return await process_activity("validate_reference_data", params, process_fn)


# ============================================================================
# ATOMIC FUNCTION: subset_reference_to_panel
# ============================================================================


@activity.defn
async def subset_reference_to_panel(
    data_uri: str,
    output_uri: str,
    panel_genes_uri: Optional[str] = None,
    panel_genes: Optional[List[str]] = None,
    xenium_data_uri: Optional[str] = None,
    min_overlap_fraction: float = 0.2,
    gene_name_column: Optional[str] = None,
) -> str:
    """
    Subset reference data to only the genes in the target panel.

    Takes CellxGene reference data and subsets it to only include genes present
    in the target panel (e.g., Xenium ~400 genes). This is essential for training
    panel-specific CellTypist models that will have 100% gene utilization.

    Panel genes can be provided in three ways (in priority order):
    1. xenium_data_uri: Extract from existing Xenium data (recommended)
    2. panel_genes_uri: Load from JSON/CSV file
    3. panel_genes: Provide directly as list

    Args:
        data_uri: GCS URI to reference h5ad (from validate_reference_data)
        output_uri: GCS URI for output metadata
        panel_genes_uri: Optional URI to JSON/CSV with panel genes
        panel_genes: Optional list of panel gene names
        xenium_data_uri: Optional URI to Xenium data to extract panel genes
        min_overlap_fraction: Minimum acceptable gene overlap (default: 0.2 = 20%)
        gene_name_column: Column in var for gene names (default: var_names)

    Returns:
        GCS URI to subset reference data

    @metadata
    is_composite: false
    ready_for_production: false
    leads_to: [train_celltypist_model]
    follows_from: [validate_reference_data]
    category: preprocessing
    runtime: "3-10min"
    memory: "12GB"
    algorithm_type: gene_subsetting
    method_type: panel_intersection
    prerequisites:
      - validated_reference: true
      - panel_genes_required: true
    output_guarantees:
      - subset_h5ad: true
      - gene_overlap_stats: true
      - missing_genes_list: true
    best_practices: |
      Subsets reference data to panel genes for custom model training:
      1. Loads reference h5ad from previous step
      2. Gets panel genes from xenium data, file, or direct list
      3. Computes intersection of reference and panel genes
      4. Validates minimum overlap (fails if too low)
      5. Saves subset with only overlapping genes

      When to use:
      - After validating reference data
      - Before training custom CellTypist model
      - When reference has full transcriptome but panel is targeted

      Panel Gene Sources (priority order):
      1. xenium_data_uri: Extract from Xenium h5ad (most reliable)
      2. panel_genes_uri: Load from JSON/CSV file
      3. panel_genes: Direct list (for testing/debugging)
    workflow_position: preprocessing
    recommended_sequence: [validate_reference_data, subset_reference_to_panel, train_celltypist_model]
    suppress_from_mcp: true
    """
    params = {
        "data_uri": data_uri,
        "output_uri": output_uri,
        "panel_genes_uri": panel_genes_uri,
        "panel_genes": panel_genes,
        "xenium_data_uri": xenium_data_uri,
        "min_overlap_fraction": min_overlap_fraction,
        "gene_name_column": gene_name_column,
    }

    async def process_fn(input_data: dict, params: dict) -> dict:
        """Subset reference to panel genes."""
        print(f"\n{'='*80}")
        print("subset_reference_to_panel - Gene Panel Subsetting")
        print(f"{'='*80}")

        input_path = Path("./input/data.h5ad")
        output_path = Path("./output/data.h5ad")

        # Load reference data
        print("\n[1/4] Loading reference data...")
        adata_ref = sc.read_h5ad(input_path)
        print(f"  Reference: {adata_ref.n_obs:,} cells × {adata_ref.n_vars:,} genes")

        # Get panel genes
        print("\n[2/4] Getting panel genes...")
        target_genes = None

        if params["xenium_data_uri"]:
            print(f"  Extracting from Xenium data: {params['xenium_data_uri']}")
            # Download Xenium data to get panel
            xenium_path = Path("./input/xenium_panel.h5ad")
            if params["xenium_data_uri"].startswith("gs://"):
                download_from_gcs(params["xenium_data_uri"], xenium_path)
            else:
                import shutil

                local_path = params["xenium_data_uri"].replace("file://", "")
                shutil.copy(local_path, xenium_path)

            adata_xenium = sc.read_h5ad(xenium_path)
            target_genes = list(adata_xenium.var_names)
            print(f"  Extracted {len(target_genes)} genes from Xenium panel")

        elif params["panel_genes_uri"]:
            print(f"  Loading from file: {params['panel_genes_uri']}")
            panel_path = Path("./input/panel_genes")
            if params["panel_genes_uri"].startswith("gs://"):
                download_from_gcs(params["panel_genes_uri"], panel_path)
            else:
                import shutil

                local_path = params["panel_genes_uri"].replace("file://", "")
                shutil.copy(local_path, panel_path)

            # Try JSON first, then CSV
            try:
                with open(panel_path) as f:
                    target_genes = json.load(f)
            except json.JSONDecodeError:
                df = pd.read_csv(panel_path)
                target_genes = df.iloc[:, 0].tolist()

            print(f"  Loaded {len(target_genes)} genes from file")

        elif params["panel_genes"]:
            target_genes = params["panel_genes"]
            print(f"  Using provided gene list: {len(target_genes)} genes")

        else:
            raise ValueError(
                "Must provide one of: xenium_data_uri, panel_genes_uri, or panel_genes"
            )

        # Compute gene overlap
        print("\n[3/4] Computing gene overlap...")
        ref_genes = set(adata_ref.var_names)
        panel_genes_set = set(target_genes)

        overlap_genes = ref_genes & panel_genes_set
        missing_from_ref = panel_genes_set - ref_genes

        overlap_fraction = len(overlap_genes) / len(panel_genes_set)

        print(f"  Panel genes: {len(panel_genes_set)}")
        print(f"  Reference genes: {len(ref_genes)}")
        print(f"  Overlap genes: {len(overlap_genes)}")
        print(f"  Overlap fraction: {overlap_fraction:.1%}")

        if missing_from_ref:
            print(f"  Panel genes missing from reference: {len(missing_from_ref)}")
            if len(missing_from_ref) <= 10:
                print(f"    {list(missing_from_ref)}")

        # Validate overlap
        if overlap_fraction < params["min_overlap_fraction"]:
            raise ValueError(
                f"Gene overlap {overlap_fraction:.1%} is below minimum "
                f"{params['min_overlap_fraction']:.1%}. "
                f"Consider using a different reference dataset."
            )

        # Subset reference to overlapping genes
        print("\n[4/4] Subsetting reference data...")
        overlap_genes_list = sorted(list(overlap_genes))
        adata_subset = adata_ref[:, overlap_genes_list].copy()

        print(f"  Subset: {adata_subset.n_obs:,} cells × {adata_subset.n_vars:,} genes")

        # Save subset
        adata_subset.write_h5ad(output_path)
        print(f"  Saved to: {output_path}")

        return {
            "n_panel_genes": len(panel_genes_set),
            "n_reference_genes": len(ref_genes),
            "n_overlap_genes": len(overlap_genes),
            "overlap_fraction": round(overlap_fraction, 4),
            "missing_panel_genes": list(missing_from_ref)[:50],  # First 50
            "n_missing_panel_genes": len(missing_from_ref),
            "subset_shape": [adata_subset.n_obs, adata_subset.n_vars],
        }

    return await process_activity("subset_reference_to_panel", params, process_fn)


# ============================================================================
# ATOMIC FUNCTION: train_celltypist_model
# ============================================================================


@activity.defn
async def train_celltypist_model(
    data_uri: str,
    output_uri: str,
    label_column: str = "cell_type",
    model_name: str = "custom_xenium_model",
    tissue_name: str = "custom",
    model_version: str = "v1",
    use_SGD: bool = True,
    feature_selection: bool = False,
    n_jobs: int = -1,
    check_expression: bool = False,
    max_iter: int = 100,
    cv_folds: int = 5,
    report_cv_metrics: bool = True,
    register_model: bool = True,
    registry_path: str = DEFAULT_REGISTRY_PATH,
) -> str:
    """
    Train a custom CellTypist logistic regression model.

    Trains a CellTypist model on the prepared (subsetted) reference data.
    The model is a multi-class logistic regression classifier that predicts
    cell types from gene expression profiles.

    After training, the model can optionally be registered in a central GCS
    registry for reuse across samples and workflows.

    Args:
        data_uri: GCS URI to subset reference h5ad (from subset_reference_to_panel)
        output_uri: GCS URI for output metadata
        label_column: Column with cell type labels (default: "cell_type")
        model_name: Name for the trained model (default: "custom_xenium_model")
        tissue_name: Tissue type for registry organization (default: "custom")
        model_version: Version tag for the model (default: "v1")
        use_SGD: Use stochastic gradient descent - faster for large data (default: True)
        feature_selection: Perform feature selection (default: False - use all genes)
        n_jobs: Number of parallel jobs (-1 = all cores)
        check_expression: Skip expression validation (default: False)
        max_iter: Maximum training iterations (default: 100)
        cv_folds: Cross-validation folds for metrics (default: 5)
        report_cv_metrics: Compute and report CV metrics (default: True)
        register_model: Register model in central registry (default: True)
        registry_path: GCS path for model registry

    Returns:
        GCS URI to training results including model location

    @metadata
    is_composite: false
    ready_for_production: false
    leads_to: [annotate_celltypist]
    follows_from: [subset_reference_to_panel]
    category: analysis
    runtime: "5-30min"
    memory: "16GB"
    algorithm_type: model_training
    method_type: logistic_regression
    prerequisites:
      - subset_reference_required: true
      - normalized_data_required: true
      - cell_type_labels_required: true
    output_guarantees:
      - model_pkl_file: true
      - training_metrics: true
      - model_metadata: true
    best_practices: |
      Trains CellTypist logistic regression model on panel-specific reference:
      1. Loads subset reference (already filtered to panel genes)
      2. Validates data normalization (re-normalizes if needed)
      3. Trains multi-class logistic regression classifier
      4. Computes cross-validation metrics
      5. Saves model and registers in GCS registry

      Training Parameters:
      - use_SGD=True: Faster for large datasets (>10k cells), uses mini-batches
      - use_SGD=False: More stable for small datasets, full gradient descent
      - feature_selection=False: Use all genes (panel already selected)
      - check_expression=False: Skip validation (already done)

      Model Quality Indicators:
      - High CV accuracy (>80%): Model generalizes well
      - Low per-class F1: Some cell types hard to distinguish
      - Training time: Should be minutes, not hours

      When to use:
      - After subsetting reference to panel genes
      - When training tissue-specific custom models
      - For creating reusable models for a panel type

      When NOT to use:
      - If pre-trained model has >50% gene overlap (use that instead)
      - For one-off annotations (pre-trained faster)
    workflow_position: analysis
    recommended_sequence: [subset_reference_to_panel, train_celltypist_model, annotate_celltypist]
    suppress_from_mcp: true
    """
    params = {
        "data_uri": data_uri,
        "output_uri": output_uri,
        "label_column": label_column,
        "model_name": model_name,
        "tissue_name": tissue_name,
        "model_version": model_version,
        "use_SGD": use_SGD,
        "feature_selection": feature_selection,
        "n_jobs": n_jobs,
        "check_expression": check_expression,
        "max_iter": max_iter,
        "cv_folds": cv_folds,
        "report_cv_metrics": report_cv_metrics,
        "register_model": register_model,
        "registry_path": registry_path,
    }

    async def process_fn(input_data: dict, params: dict) -> dict:
        """Train CellTypist model on subset reference."""
        import celltypist
        from sklearn.model_selection import cross_val_score
        import gc

        print(f"\n{'='*80}")
        print("train_celltypist_model - Custom Model Training")
        print(f"{'='*80}")

        input_path = Path("./input/data.h5ad")
        output_dir = Path("./output")
        output_dir.mkdir(exist_ok=True)

        # Load subset reference
        print("\n[1/5] Loading subset reference data...")
        adata = sc.read_h5ad(input_path)
        print(f"  Data: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

        label_column = params["label_column"]
        if label_column not in adata.obs.columns:
            raise ValueError(f"Label column '{label_column}' not found in data")

        n_cell_types = adata.obs[label_column].nunique()
        print(f"  Cell types: {n_cell_types} (column: {label_column})")

        # Ensure proper normalization for CellTypist
        print("\n[2/5] Validating/fixing normalization...")

        # Check if data looks normalized
        if issparse(adata.X):
            sample_data = (
                adata.X[:1000].toarray() if adata.n_obs > 1000 else adata.X.toarray()
            )
        else:
            sample_data = adata.X[:1000] if adata.n_obs > 1000 else adata.X

        mean_val = float(np.mean(sample_data))
        max_val = float(np.max(sample_data))

        if mean_val > 10 or max_val > 50:
            print("  Data appears to be raw counts - normalizing...")
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            print("  ✓ Normalized to log1p(10k CPM)")
        else:
            print(
                f"  ✓ Data appears already normalized (mean={mean_val:.2f}, max={max_val:.2f})"
            )

        # Train model
        print("\n[3/5] Training CellTypist model...")
        print(
            f"  Method: {'SGD (mini-batch)' if params['use_SGD'] else 'Full gradient descent'}"
        )
        print(f"  Feature selection: {params['feature_selection']}")
        print(f"  Max iterations: {params['max_iter']}")

        t0 = datetime.now()

        model = celltypist.train(
            adata,
            labels=label_column,
            check_expression=params["check_expression"],
            use_SGD=params["use_SGD"],
            feature_selection=params["feature_selection"],
            n_jobs=params["n_jobs"],
            max_iter=params["max_iter"],
        )

        training_time = (datetime.now() - t0).total_seconds()
        print(f"  ✓ Training complete in {training_time:.1f}s")

        # Model info
        print("\n  Model Statistics:")
        print(f"    Cell types: {len(model.cell_types)}")
        print(f"    Features (genes): {len(model.features)}")

        # Save model locally
        model_filename = f"{params['model_name']}.pkl"
        local_model_path = output_dir / model_filename
        model.write(str(local_model_path))
        print(f"  Saved model to: {local_model_path}")

        # Cross-validation metrics (optional)
        cv_metrics = {}
        if params["report_cv_metrics"]:
            print("\n[4/5] Computing cross-validation metrics...")
            try:
                # Get classifier from model
                classifier = model.classifier
                X = adata.X.toarray() if issparse(adata.X) else adata.X
                y = adata.obs[label_column].values

                # Accuracy via CV
                cv_scores = cross_val_score(
                    classifier,
                    X,
                    y,
                    cv=min(params["cv_folds"], n_cell_types),
                    scoring="accuracy",
                    n_jobs=params["n_jobs"],
                )

                cv_metrics = {
                    "cv_accuracy_mean": float(np.mean(cv_scores)),
                    "cv_accuracy_std": float(np.std(cv_scores)),
                    "cv_folds": len(cv_scores),
                }
                print(
                    f"  CV Accuracy: {cv_metrics['cv_accuracy_mean']:.3f} ± {cv_metrics['cv_accuracy_std']:.3f}"
                )

            except Exception as e:
                print(f"  ⚠️  CV metrics failed: {e}")
                cv_metrics = {"error": str(e)}
        else:
            print("\n[4/5] Skipping CV metrics (report_cv_metrics=False)")

        # Register model in GCS
        print("\n[5/5] Registering model...")
        model_uri = None
        metadata_uri = None

        if params["register_model"]:
            registry_path = params["registry_path"]
            model_dir = f"{params['tissue_name']}_{params['model_name']}_{params['model_version']}"

            if registry_path.startswith("gs://"):
                # GCS registration
                model_uri = f"{registry_path}/{model_dir}/model.pkl"
                metadata_uri = f"{registry_path}/{model_dir}/metadata.json"

                print(f"  Uploading to: {model_uri}")
                upload_to_gcs(local_model_path, model_uri)

                # Create and upload metadata
                metadata = {
                    "model_name": params["model_name"],
                    "tissue": params["tissue_name"],
                    "version": params["model_version"],
                    "created": datetime.now().isoformat(),
                    "n_cells_trained": adata.n_obs,
                    "n_genes": adata.n_vars,
                    "n_cell_types": len(model.cell_types),
                    "cell_types": model.cell_types.tolist(),
                    "features": model.features.tolist()[:100],  # First 100 genes
                    "training_time_seconds": training_time,
                    "cv_metrics": cv_metrics,
                    "training_params": {
                        "use_SGD": params["use_SGD"],
                        "feature_selection": params["feature_selection"],
                        "max_iter": params["max_iter"],
                        "label_column": label_column,
                    },
                    "celltypist_version": celltypist.__version__,
                }

                # Upload metadata
                metadata_path = output_dir / "metadata.json"
                with open(metadata_path, "w") as f:
                    json.dump(metadata, f, indent=2)
                upload_to_gcs(metadata_path, metadata_uri)

                print(f"  ✓ Registered: {model_uri}")
            else:
                # Local registration
                local_registry = Path(registry_path.replace("file://", ""))
                model_dir_path = local_registry / model_dir
                model_dir_path.mkdir(parents=True, exist_ok=True)

                import shutil

                model_uri = f"file://{model_dir_path}/model.pkl"
                shutil.copy(local_model_path, model_dir_path / "model.pkl")

                metadata = {
                    "model_name": params["model_name"],
                    "tissue": params["tissue_name"],
                    "version": params["model_version"],
                    "created": datetime.now().isoformat(),
                    "n_cells_trained": adata.n_obs,
                    "n_genes": adata.n_vars,
                    "n_cell_types": len(model.cell_types),
                    "cell_types": model.cell_types.tolist(),
                    "training_time_seconds": training_time,
                    "cv_metrics": cv_metrics,
                    "celltypist_version": celltypist.__version__,
                }
                with open(model_dir_path / "metadata.json", "w") as f:
                    json.dump(metadata, f, indent=2)

                print(f"  ✓ Registered locally: {model_uri}")
        else:
            model_uri = f"file://{local_model_path.absolute()}"
            print(f"  Model saved (not registered): {model_uri}")

        # Cleanup
        gc.collect()

        return {
            "model_uri": model_uri,
            "model_name": params["model_name"],
            "tissue": params["tissue_name"],
            "version": params["model_version"],
            "training_stats": {
                "n_cells_trained": adata.n_obs,
                "n_genes": adata.n_vars,
                "n_cell_types": len(model.cell_types),
                "cell_types": model.cell_types.tolist()[:20],  # First 20
                "training_time_seconds": training_time,
            },
            "cv_metrics": cv_metrics,
            "registered": params["register_model"],
        }

    return await process_activity("train_celltypist_model", params, process_fn)


# ============================================================================
# COMPOSITE FUNCTION: prepare_celltypist_model
# ============================================================================


@activity.defn
async def prepare_celltypist_model(
    reference_uri: str,
    output_uri: str,
    panel_genes_uri: Optional[str] = None,
    panel_genes: Optional[List[str]] = None,
    xenium_data_uri: Optional[str] = None,
    label_column: str = "cell_type",
    tissue_name: str = "custom",
    model_version: str = "v1",
    use_SGD: bool = True,
    min_cells_per_type: int = 50,
    min_panel_overlap: float = 0.2,
    register_model: bool = True,
    registry_path: str = DEFAULT_REGISTRY_PATH,
) -> str:
    """
    Prepare a custom CellTypist model from CellxGene reference data.

    This composite function orchestrates the complete model preparation workflow:
    1. Validate reference data (check format, labels, normalization)
    2. Subset reference to panel genes (Xenium ~400 genes)
    3. Train CellTypist logistic regression model
    4. Register model in central GCS registry

    The trained model can then be used for cell type annotation with
    annotate_celltypist(model_uri=...).

    Panel genes can be provided in three ways (in priority order):
    1. xenium_data_uri: Extract from existing Xenium data (recommended)
    2. panel_genes_uri: Load from JSON/CSV file
    3. panel_genes: Provide directly as list

    Args:
        reference_uri: GCS/local URI to CellxGene h5ad file
        output_uri: GCS URI for workflow output metadata
        panel_genes_uri: Optional URI to JSON/CSV with panel genes
        panel_genes: Optional list of panel gene names
        xenium_data_uri: Optional URI to Xenium data to extract panel genes
        label_column: Column with cell type labels (default: "cell_type")
        tissue_name: Tissue type for naming and organization (default: "custom")
        model_version: Version tag for the model (default: "v1")
        use_SGD: Use stochastic gradient descent (default: True)
        min_cells_per_type: Minimum cells per type for training (default: 50)
        min_panel_overlap: Minimum gene overlap fraction (default: 0.2)
        register_model: Register in central registry (default: True)
        registry_path: GCS path for model registry

    Returns:
        GCS URI to output metadata containing model_uri

    @metadata
    is_composite: true
    ready_for_production: false
    leads_to: [annotate_celltypist]
    follows_from: []
    category: analysis
    runtime: "10-60min"
    memory: "16GB"
    available_methods: [validate, subset, train, full]
    default_method: full
    recommended_for: [xenium, targeted_panels, low_gene_overlap]
    prerequisites:
      - reference_h5ad_available: true
      - panel_genes_defined: true
      - cell_type_labels_present: true
    output_guarantees:
      - trained_model_uri: GCS URI to .pkl model file
      - model_metadata: training stats, gene overlap, cell types
      - ready_for_annotation: true
    best_practices: |
      Orchestrates custom CellTypist model training for targeted gene panels.

      Reference Data Selection:
      - Use tissue-matched CellxGene data (e.g., colon atlas for colon Xenium)
      - Ensure adequate cell type diversity (50+ cells per type)
      - Prefer normalized data (log1p to 10k counts per cell)

      Panel Gene Sources (priority order):
      1. xenium_data_uri: Extract from existing Xenium data (recommended)
      2. panel_genes_uri: Load from JSON/CSV file
      3. panel_genes: Provide directly as list

      Expected Gene Overlap:
      - Good overlap: >60% panel genes in reference
      - Acceptable overlap: >40% panel genes in reference
      - Poor overlap: <30% (may need different reference)

      When to use:
      - Pre-trained models have <50% gene overlap with your panel
      - Training tissue-specific models for reuse
      - Xenium or other targeted panel data

      When NOT to use:
      - Pre-trained model has good gene overlap (>50%)
      - One-off annotations (pre-trained is faster)
      - No suitable CellxGene reference available
    workflow_position: model_preparation
    recommended_sequence: [prepare_celltypist_model, annotate_celltypist]
    suppress_from_mcp: true
    """
    params = {
        "data_uri": reference_uri,
        "output_uri": output_uri,
        "panel_genes_uri": panel_genes_uri,
        "panel_genes": panel_genes,
        "xenium_data_uri": xenium_data_uri,
        "label_column": label_column,
        "tissue_name": tissue_name,
        "model_version": model_version,
        "use_SGD": use_SGD,
        "min_cells_per_type": min_cells_per_type,
        "min_panel_overlap": min_panel_overlap,
        "register_model": register_model,
        "registry_path": registry_path,
    }

    async def process_fn(input_data: dict, params: dict) -> dict:
        """Orchestrate model preparation workflow."""
        print(f"\n{'='*80}")
        print("prepare_celltypist_model - Custom Model Preparation")
        print(f"{'='*80}")

        input_path = Path("./input/data.h5ad")
        output_dir = Path("./output")
        output_dir.mkdir(exist_ok=True)

        timings = {}

        # ================================================================
        # STEP 1: VALIDATE REFERENCE DATA
        # ================================================================
        print("\n" + "=" * 60)
        print("STEP 1: Validating Reference Data")
        print("=" * 60)

        t0 = datetime.now()

        adata = sc.read_h5ad(input_path)
        print(f"  Loaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

        # Check label column
        label_column = params["label_column"]
        if label_column not in adata.obs.columns:
            # Try common alternatives
            alternatives = ["cell_type", "cell_ontology_class", "celltype", "CellType"]
            found = None
            for alt in alternatives:
                if alt in adata.obs.columns:
                    found = alt
                    break
            if found:
                print(f"  ⚠️  Column '{label_column}' not found, using '{found}'")
                label_column = found
            else:
                raise ValueError(
                    f"Label column '{label_column}' not found. "
                    f"Available: {list(adata.obs.columns)[:20]}"
                )

        # Cell type stats
        cell_type_counts = adata.obs[label_column].value_counts()
        n_cell_types = len(cell_type_counts)
        below_threshold = (cell_type_counts < params["min_cells_per_type"]).sum()

        print(f"  Cell types: {n_cell_types}")
        print(
            f"  Cell types below {params['min_cells_per_type']} cells: {below_threshold}"
        )

        # Normalization check
        if issparse(adata.X):
            sample = (
                adata.X[:1000].toarray() if adata.n_obs > 1000 else adata.X.toarray()
            )
        else:
            sample = adata.X[:1000] if adata.n_obs > 1000 else adata.X

        mean_val = float(np.mean(sample))
        needs_normalization = mean_val > 10

        if needs_normalization:
            print(f"  ⚠️  Data appears to be raw counts (mean={mean_val:.2f})")
            print("      Will normalize before training")
        else:
            print(f"  ✓ Data appears normalized (mean={mean_val:.2f})")

        timings["validation"] = (datetime.now() - t0).total_seconds()
        print(f"  Validation time: {timings['validation']:.1f}s")

        # ================================================================
        # STEP 2: GET PANEL GENES AND SUBSET
        # ================================================================
        print("\n" + "=" * 60)
        print("STEP 2: Subsetting to Panel Genes")
        print("=" * 60)

        t0 = datetime.now()

        # Get panel genes
        target_genes = None

        if params["xenium_data_uri"]:
            print(f"  Source: Xenium data ({params['xenium_data_uri']})")
            xenium_path = Path("./input/xenium_temp.h5ad")

            if params["xenium_data_uri"].startswith("gs://"):
                download_from_gcs(params["xenium_data_uri"], xenium_path)
            else:
                import shutil

                local_path = params["xenium_data_uri"].replace("file://", "")
                shutil.copy(local_path, xenium_path)

            adata_xenium = sc.read_h5ad(xenium_path)
            target_genes = list(adata_xenium.var_names)
            print(f"  Extracted {len(target_genes)} genes from Xenium panel")

        elif params["panel_genes_uri"]:
            print(f"  Source: Gene list file ({params['panel_genes_uri']})")
            panel_path = Path("./input/panel_temp")

            if params["panel_genes_uri"].startswith("gs://"):
                download_from_gcs(params["panel_genes_uri"], panel_path)
            else:
                import shutil

                local_path = params["panel_genes_uri"].replace("file://", "")
                shutil.copy(local_path, panel_path)

            try:
                with open(panel_path) as f:
                    target_genes = json.load(f)
            except json.JSONDecodeError:
                df = pd.read_csv(panel_path)
                target_genes = df.iloc[:, 0].tolist()

            print(f"  Loaded {len(target_genes)} genes from file")

        elif params["panel_genes"]:
            target_genes = params["panel_genes"]
            print(f"  Source: Direct list ({len(target_genes)} genes)")

        else:
            raise ValueError(
                "Must provide one of: xenium_data_uri, panel_genes_uri, or panel_genes"
            )

        # Compute overlap
        ref_genes = set(adata.var_names)
        panel_genes_set = set(target_genes)
        overlap_genes = ref_genes & panel_genes_set
        overlap_fraction = len(overlap_genes) / len(panel_genes_set)

        print(f"  Panel genes: {len(panel_genes_set)}")
        print(f"  Reference genes: {len(ref_genes)}")
        print(f"  Overlap: {len(overlap_genes)} ({overlap_fraction:.1%})")

        if overlap_fraction < params["min_panel_overlap"]:
            raise ValueError(
                f"Gene overlap {overlap_fraction:.1%} is below minimum "
                f"{params['min_panel_overlap']:.1%}. "
                f"Consider using a different reference dataset."
            )

        # Subset
        overlap_genes_list = sorted(list(overlap_genes))
        adata_subset = adata[:, overlap_genes_list].copy()
        print(f"  Subset: {adata_subset.n_obs:,} cells × {adata_subset.n_vars:,} genes")

        timings["subsetting"] = (datetime.now() - t0).total_seconds()
        print(f"  Subsetting time: {timings['subsetting']:.1f}s")

        # ================================================================
        # STEP 3: NORMALIZE IF NEEDED
        # ================================================================
        if needs_normalization:
            print("\n" + "=" * 60)
            print("STEP 3: Normalizing Data")
            print("=" * 60)

            t0 = datetime.now()
            sc.pp.normalize_total(adata_subset, target_sum=1e4)
            sc.pp.log1p(adata_subset)
            timings["normalization"] = (datetime.now() - t0).total_seconds()
            print(
                f"  ✓ Normalized to log1p(10k CPM) in {timings['normalization']:.1f}s"
            )

        # ================================================================
        # STEP 4: TRAIN MODEL
        # ================================================================
        print("\n" + "=" * 60)
        print("STEP 4: Training CellTypist Model")
        print("=" * 60)

        import celltypist

        t0 = datetime.now()

        print(
            f"  Method: {'SGD (mini-batch)' if params['use_SGD'] else 'Full gradient descent'}"
        )
        print(f"  Label column: {label_column}")

        model = celltypist.train(
            adata_subset,
            labels=label_column,
            check_expression=False,
            use_SGD=params["use_SGD"],
            feature_selection=False,
            n_jobs=-1,
        )

        timings["training"] = (datetime.now() - t0).total_seconds()
        print(f"  ✓ Training complete in {timings['training']:.1f}s")
        print(f"    Cell types: {len(model.cell_types)}")
        print(f"    Features: {len(model.features)}")

        # ================================================================
        # STEP 5: SAVE AND REGISTER MODEL
        # ================================================================
        print("\n" + "=" * 60)
        print("STEP 5: Saving and Registering Model")
        print("=" * 60)

        # Construct model name
        full_model_name = f"{params['tissue_name']}_{params['model_version']}"

        # Save locally
        local_model_path = output_dir / f"{full_model_name}.pkl"
        model.write(str(local_model_path))
        print(f"  Saved locally: {local_model_path}")

        # Register if requested
        model_uri = None
        if params["register_model"]:
            registry = params["registry_path"]
            model_dir = f"{params['tissue_name']}_{params['model_version']}"

            if registry.startswith("gs://"):
                model_uri = f"{registry}/{model_dir}/model.pkl"
                metadata_uri = f"{registry}/{model_dir}/metadata.json"

                # Upload model
                upload_to_gcs(local_model_path, model_uri)
                print(f"  Uploaded: {model_uri}")

                # Create metadata
                metadata = {
                    "model_name": full_model_name,
                    "tissue": params["tissue_name"],
                    "version": params["model_version"],
                    "created": datetime.now().isoformat(),
                    "n_cells_trained": adata_subset.n_obs,
                    "n_genes": adata_subset.n_vars,
                    "n_cell_types": len(model.cell_types),
                    "cell_types": model.cell_types.tolist(),
                    "gene_overlap_pct": round(overlap_fraction * 100, 1),
                    "panel_genes_total": len(panel_genes_set),
                    "training_time_seconds": timings["training"],
                    "celltypist_version": celltypist.__version__,
                }

                metadata_path = output_dir / "metadata.json"
                with open(metadata_path, "w") as f:
                    json.dump(metadata, f, indent=2)
                upload_to_gcs(metadata_path, metadata_uri)

            else:
                # Local registry
                local_registry = Path(registry.replace("file://", ""))
                model_dir_path = local_registry / model_dir
                model_dir_path.mkdir(parents=True, exist_ok=True)

                import shutil

                model_uri = f"file://{model_dir_path}/model.pkl"
                shutil.copy(local_model_path, model_dir_path / "model.pkl")

                metadata = {
                    "model_name": full_model_name,
                    "tissue": params["tissue_name"],
                    "version": params["model_version"],
                    "created": datetime.now().isoformat(),
                    "n_cells_trained": adata_subset.n_obs,
                    "n_genes": adata_subset.n_vars,
                    "n_cell_types": len(model.cell_types),
                    "cell_types": model.cell_types.tolist(),
                    "gene_overlap_pct": round(overlap_fraction * 100, 1),
                }
                with open(model_dir_path / "metadata.json", "w") as f:
                    json.dump(metadata, f, indent=2)

            print(f"  ✓ Registered: {model_uri}")
        else:
            model_uri = f"file://{local_model_path.absolute()}"
            print(f"  Model saved (not registered): {model_uri}")

        # Summary
        print("\n" + "=" * 60)
        print("MODEL PREPARATION COMPLETE")
        print("=" * 60)
        print(f"  Model URI: {model_uri}")
        print(f"  Cell types: {len(model.cell_types)}")
        print(f"  Genes: {len(model.features)}")
        print(f"  Gene overlap: {overlap_fraction:.1%}")
        print(f"  Total time: {sum(timings.values()):.1f}s")

        return {
            "model_uri": model_uri,
            "model_name": full_model_name,
            "tissue": params["tissue_name"],
            "version": params["model_version"],
            "n_cells_trained": adata_subset.n_obs,
            "n_genes": adata_subset.n_vars,
            "n_cell_types": len(model.cell_types),
            "cell_types": model.cell_types.tolist()[:30],  # First 30
            "gene_overlap_pct": round(overlap_fraction * 100, 1),
            "panel_genes_total": len(panel_genes_set),
            "timings": timings,
            "registered": params["register_model"],
        }

    return await process_activity("prepare_celltypist_model", params, process_fn)


# ============================================================================
# UTILITY FUNCTION: list_registered_models
# ============================================================================


@activity.defn
async def list_registered_models(
    output_uri: str,
    registry_path: str = DEFAULT_REGISTRY_PATH,
    tissue_filter: Optional[str] = None,
) -> str:
    """
    List all registered custom CellTypist models.

    Queries the model registry to list available custom models with their
    metadata. Useful for discovering pre-trained models before annotation.

    Args:
        output_uri: GCS URI for output
        registry_path: GCS path to model registry
        tissue_filter: Optional filter by tissue type

    Returns:
        GCS URI to JSON listing of available models

    @metadata
    is_composite: false
    ready_for_production: false
    leads_to: [annotate_celltypist]
    follows_from: []
    category: utility
    runtime: "1min"
    memory: "2GB"
    algorithm_type: registry_query
    method_type: list_models
    prerequisites: []
    output_guarantees:
      - model_list: true
      - model_metadata: true
    best_practices: |
      Lists available custom CellTypist models from the registry.

      When to use:
      - Before annotation to discover available models
      - To check if a tissue-specific model exists
      - To get model metadata (cell types, gene count)

      When NOT to use:
      - Already know which model to use
      - Using pre-trained CellTypist models (not custom)
    workflow_position: utility
    suppress_from_mcp: true
    """
    from google.cloud import storage

    print(f"\n{'='*80}")
    print("list_registered_models - Registry Query")
    print(f"{'='*80}")

    models = []

    if registry_path.startswith("gs://"):
        # Parse GCS path
        bucket_name = registry_path.replace("gs://", "").split("/")[0]
        prefix = "/".join(registry_path.replace("gs://", "").split("/")[1:])

        print(f"  Registry: {registry_path}")

        try:
            client = storage.Client()
            bucket = client.bucket(bucket_name)

            # List all model directories
            blobs = bucket.list_blobs(prefix=prefix + "/")

            seen_dirs = set()
            for blob in blobs:
                # Extract model directory name
                parts = blob.name.replace(prefix + "/", "").split("/")
                if len(parts) >= 2 and parts[0] not in seen_dirs:
                    model_dir = parts[0]
                    seen_dirs.add(model_dir)

                    # Try to read metadata
                    metadata_path = f"{prefix}/{model_dir}/metadata.json"
                    metadata_blob = bucket.blob(metadata_path)

                    if metadata_blob.exists():
                        try:
                            content = metadata_blob.download_as_text()
                            metadata = json.loads(content)
                            metadata["model_dir"] = model_dir
                            metadata["model_uri"] = (
                                f"gs://{bucket_name}/{prefix}/{model_dir}/model.pkl"
                            )

                            # Apply filter
                            if tissue_filter:
                                if (
                                    tissue_filter.lower()
                                    in metadata.get("tissue", "").lower()
                                ):
                                    models.append(metadata)
                            else:
                                models.append(metadata)

                        except Exception as e:
                            print(
                                f"    ⚠️  Could not read metadata for {model_dir}: {e}"
                            )

            print(f"  Found {len(models)} models")

        except Exception as e:
            print(f"  ❌ Error accessing registry: {e}")

    else:
        # Local registry
        local_path = Path(registry_path.replace("file://", ""))
        if local_path.exists():
            for model_dir in local_path.iterdir():
                if model_dir.is_dir():
                    metadata_file = model_dir / "metadata.json"
                    if metadata_file.exists():
                        with open(metadata_file) as f:
                            metadata = json.load(f)
                        metadata["model_dir"] = model_dir.name
                        metadata["model_uri"] = f"file://{model_dir}/model.pkl"

                        if tissue_filter:
                            if (
                                tissue_filter.lower()
                                in metadata.get("tissue", "").lower()
                            ):
                                models.append(metadata)
                        else:
                            models.append(metadata)

            print(f"  Found {len(models)} models in local registry")

    # Sort by creation date (newest first)
    models.sort(key=lambda x: x.get("created", ""), reverse=True)

    # Print summary
    if models:
        print("\n  Available Models:")
        for m in models[:10]:  # Show first 10
            print(f"    - {m.get('model_name', 'unknown')}")
            print(f"      Tissue: {m.get('tissue', 'unknown')}")
            print(f"      Cell types: {m.get('n_cell_types', '?')}")
            print(f"      Genes: {m.get('n_genes', '?')}")

    result = {
        "registry_path": registry_path,
        "tissue_filter": tissue_filter,
        "n_models": len(models),
        "models": models,
    }

    # Write to output
    universal_write(output_uri, result)

    return output_uri
