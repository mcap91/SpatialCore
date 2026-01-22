"""
CellTypist model training utilities.

This module provides utilities for:
1. Combining multiple reference datasets for training
2. Training custom CellTypist models
3. Panel gene subsetting for spatial transcriptomics

For spatial data (e.g., Xenium), custom models trained on panel-specific genes
achieve ~100% gene utilization vs ~8% with pre-trained models.

References:
    - CellTypist: https://www.celltypist.org/
    - Domínguez Conde et al., Science (2022)
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union, Any
from datetime import datetime
import json
import gc

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from spatialcore.core.logging import get_logger
from spatialcore.core.utils import (
    load_ensembl_to_hugo_mapping,
    normalize_gene_names,
    check_normalization_status,
)
from spatialcore.annotation.loading import (
    load_adata_backed,
    ensure_normalized,
    get_available_memory_gb,
)
from spatialcore.annotation.validation import validate_cell_type_column
from spatialcore.annotation.acquisition import resolve_uri_to_local

logger = get_logger(__name__)

# Default cache directory for downloaded references
DEFAULT_CACHE_DIR = Path.home() / ".spatialcore" / "cache" / "references"

# Default labels to exclude from training (ambiguous/uninformative)
# These are matched EXACTLY (case-sensitive, no partial matching)
# "unknown cells" would NOT be filtered (not an exact match to "unknown")
DEFAULT_EXCLUDE_LABELS = [
    "unknown",
    "Unknown",
    "UNKNOWN",
    "unassigned",
    "Unassigned",
    "na",
    "NA",
    "N/A",
    "n/a",
    "nan",  # Python str(np.nan) produces "nan"
    "NaN",
    "NAN",
    "none",
    "None",
    "null",
    "doublet",
    "Doublet",
    "multiplet",
    "Multiplet",
    "low quality",
    "Low quality",
    "low_count",
    "Low_count",
    "LOW_COUNT",
    "low count",
    "Low count",
]


# ============================================================================
# Reference Combination
# ============================================================================

def combine_references(
    reference_paths: List[Union[str, Path]],
    label_columns: List[str],
    output_column: str = "original_label",
    max_cells_per_ref: int = 100000,
    target_genes: Optional[List[str]] = None,
    normalize_data: bool = True,
    random_state: int = 42,
    validate_labels: bool = True,
    min_cells_per_type: int = 10,
    strict_validation: bool = False,
    cache_dir: Optional[Path] = None,
    exclude_labels: Optional[List[str]] = None,
    filter_min_cells: bool = True,
) -> ad.AnnData:
    """
    Combine multiple reference datasets for CellTypist training.

    This function handles:
    1. Memory-efficient loading with per-reference cell caps (stratified)
    2. Gene name normalization (Ensembl → HUGO)
    3. Expression normalization (log1p to 10k)
    4. Gene intersection (with optional panel gene subsetting)
    5. Concatenation with source tracking

    IMPORTANT: This function does NOT perform post-combine balancing.
    Call subsample_balanced() on the output for source-aware balancing.
    For semantic grouping by ontology ID, first run add_ontology_ids() then use
    subsample_balanced(group_by_column="cell_type_ontology_term_id").

    Parameters
    ----------
    reference_paths : List[str or Path]
        Paths to reference h5ad files. Supports:

        - Local paths: ``/data/references/lung.h5ad``
        - GCS URIs: ``gs://bucket/references/lung.h5ad``
        - S3 URIs: ``s3://bucket/references/lung.h5ad``

        Cloud files are automatically downloaded to cache_dir and loaded
        in memory-efficient backed mode.
    label_columns : List[str]
        Cell type label column for each reference.
    output_column : str, default "original_label"
        Column name for unified cell type labels in output.
        Use "original_label" (new default) for clarity that this is the raw
        source label before any harmonization.
    max_cells_per_ref : int, default 100000
        Maximum cells to load per reference. Uses stratified sampling
        to preserve natural cell type proportions within each reference.
        This is for MEMORY MANAGEMENT during loading, not training balance.
    target_genes : List[str], optional
        Panel genes to subset to (e.g., from spatial data via get_panel_genes()).
        If provided, each reference is subset to these genes before
        finding the intersection. This ensures maximum gene utilization.
    normalize_data : bool, default True
        Whether to ensure log1p(10k) normalization.
    random_state : int, default 42
        Random seed for reproducibility.
    validate_labels : bool, default True
        Run cell type column validation before combining. Checks for
        null values, suspicious patterns, and cardinality issues.
    min_cells_per_type : int, default 10
        Minimum cells required per cell type. Used for validation warnings
        and, when ``filter_min_cells=True``, to remove low-count types
        after concatenation.
    strict_validation : bool, default False
        If True, fail on validation warnings (not just errors).
    cache_dir : Path, optional
        Directory for caching downloaded cloud files. Defaults to
        ``~/.spatialcore/cache/references/``. Only used for gs:// and s3:// URIs.
    exclude_labels : List[str], optional
        Cell type labels to exclude from the combined output. Uses exact
        case-sensitive matching (no partial matches). Cells with these labels
        are removed after concatenation. If None (default), uses
        DEFAULT_EXCLUDE_LABELS which includes common ambiguous labels like
        "unknown", "Unknown", "unassigned", "NA", "doublet", etc.
        Pass an empty list ``[]`` to disable label filtering entirely.
    filter_min_cells : bool, default True
        If True, remove cell types with fewer than ``min_cells_per_type``
        cells after concatenation. If False, only warn about low-count
        types (original behavior). Filtering is recommended for training
        as singleton or very-low-count cell types can destabilize models.

    Returns
    -------
    AnnData
        Combined reference data with:
        - Unified cell type labels in output_column
        - Source tracking in .obs["reference_source"]
        - Gene intersection applied
        - Ready for subsample_balanced()

    Notes
    -----
    **Two-stage workflow for multi-reference training:**

    1. **combine_references()** (this function): Load, normalize, intersect genes
    2. **subsample_balanced()**: Source-aware balancing for training quality

    **For spatial transcriptomics**, providing `target_genes` (panel genes)
    before training is critical:
    - Pre-trained models: ~8% gene overlap with 400-gene panel
    - Custom panel models: 100% gene overlap

    Examples
    --------
    >>> from spatialcore.annotation import (
    ...     combine_references, subsample_balanced, get_panel_genes
    ... )
    >>> # Step 1: Get panel genes from spatial data
    >>> panel_genes = get_panel_genes(xenium_adata)
    >>>
    >>> # Step 2: Combine references (loading + gene intersection)
    >>> # Supports local paths and cloud URIs (GCS, S3)
    >>> combined = combine_references(
    ...     reference_paths=[
    ...         "gs://my-bucket/references/hlca.h5ad",  # GCS
    ...         "s3://my-bucket/references/liver.h5ad",  # S3
    ...         "/local/data/study3.h5ad",  # Local
    ...     ],
    ...     label_columns=["cell_type", "cell_type", "cell_type"],
    ...     max_cells_per_ref=100000,  # Memory cap during loading
    ...     target_genes=panel_genes,   # Intersect with spatial panel
    ... )
    >>> # Output: combined AnnData with reference_source tracked in .obs
    >>>
    >>> # Step 3: Source-aware balancing (SEPARATE STEP)
    >>> balanced = subsample_balanced(
    ...     combined,
    ...     label_column="original_label",
    ...     max_cells_per_type=5000,
    ...     source_balance="proportional",
    ... )
    """
    if len(reference_paths) != len(label_columns):
        raise ValueError(
            f"Number of paths ({len(reference_paths)}) must match "
            f"number of label columns ({len(label_columns)})"
        )

    # Initialize cache directory for cloud downloads
    if cache_dir is None:
        cache_dir = DEFAULT_CACHE_DIR
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Load gene mapping once at start (reused for all references)
    logger.info("Loading Ensembl to HUGO gene mapping...")
    ensembl_to_hugo = load_ensembl_to_hugo_mapping()
    logger.info(f"Loaded {len(ensembl_to_hugo):,} gene mappings")

    adatas = []
    reference_names = []
    normalization_statuses = []
    validation_results = []  # Track validation results for each reference

    for i, (ref_path, label_col) in enumerate(zip(reference_paths, label_columns)):
        ref_path_str = str(ref_path)

        # Extract source name for logging (handle URIs and local paths)
        if ref_path_str.startswith(("gs://", "s3://")):
            # Cloud URI: extract filename from URI
            source_name = Path(ref_path_str.split("/")[-1]).stem
            logger.info(f"\n[{i+1}/{len(reference_paths)}] Loading: {ref_path_str}")
        else:
            source_name = Path(ref_path_str).stem
            logger.info(f"\n[{i+1}/{len(reference_paths)}] Loading: {Path(ref_path_str).name}")

        # Resolve URI to local path (downloads cloud files if needed)
        local_path = resolve_uri_to_local(ref_path_str, cache_dir)

        # Memory-efficient loading with backed mode for large files
        adata = load_adata_backed(
            path=local_path,
            max_cells=max_cells_per_ref,
            label_column=label_col,
            random_state=random_state,
        )

        # Validate cell type column BEFORE gene normalization
        if validate_labels:
            logger.info(f"  Validating cell type column: {label_col}")
            val_result = validate_cell_type_column(
                adata,
                label_col,
                min_cells_per_type=min_cells_per_type,
            )
            validation_results.append({
                "path": ref_path_str,
                "column": label_col,
                "is_valid": val_result.is_valid,
                "n_cell_types": val_result.n_cell_types,
                "n_cells": val_result.n_cells,
                "errors": [str(e) for e in val_result.errors],
                "warnings": [str(w) for w in val_result.warnings],
            })

            if not val_result.is_valid:
                raise ValueError(
                    f"Validation failed for {source_name}:\n{val_result.summary()}"
                )

            if strict_validation and val_result.warnings:
                raise ValueError(
                    f"Validation warnings (strict mode) for {source_name}:\n"
                    f"{val_result.summary()}"
                )

        # Normalize gene names (Ensembl → HUGO)
        adata = normalize_gene_names(adata, ensembl_to_hugo)

        # Record normalization status before any gene subsetting
        status = None
        if normalize_data:
            status = check_normalization_status(adata)
        normalization_statuses.append(status)

        # Copy cell type labels to unified column
        if label_col not in adata.obs.columns:
            available = list(adata.obs.columns)
            raise ValueError(
                f"Label column '{label_col}' not found in {source_name}. "
                f"Available columns: {available}"
            )
        adata.obs[output_column] = adata.obs[label_col].astype(str)

        # Add source reference info (use source name from URI or local path)
        adata.obs["reference_source"] = source_name

        reference_names.append(source_name)
        adatas.append(adata)
        gc.collect()

    # Subset to target genes if provided (BEFORE finding shared genes)
    if target_genes:
        logger.info(f"\nSubsetting to {len(target_genes)} target genes...")
        target_set = set(target_genes)
        for i, adata in enumerate(adatas):
            overlap = list(set(adata.var_names) & target_set)
            if len(overlap) == 0:
                raise ValueError(
                    f"No overlap between reference {i} and target genes. "
                    f"Check gene name format (HUGO symbols expected)."
                )
            adatas[i] = adata[:, overlap].copy()
            logger.info(f"  Reference {i}: {len(overlap)} genes after subset")

    # Find shared genes (inner join)
    logger.info("\nFinding shared genes across all references...")
    shared_genes = set(adatas[0].var_names)
    for adata in adatas[1:]:
        shared_genes &= set(adata.var_names)

    if len(shared_genes) == 0:
        raise ValueError(
            "No shared genes found across references! "
            "Check that gene names are in the same format (HUGO symbols)."
        )

    logger.info(f"  Shared genes: {len(shared_genes):,}")

    # Subset all to shared genes (sorted for consistency)
    shared_genes_sorted = sorted(shared_genes)
    for i in range(len(adatas)):
        adatas[i] = adatas[i][:, shared_genes_sorted].copy()

    # Normalize after final gene subsetting (panel + shared genes)
    if normalize_data:
        logger.info("\nNormalizing references after gene subsetting...")
        for i, adata in enumerate(adatas):
            status = normalization_statuses[i]
            source_name = reference_names[i]
            if status is None:
                raise ValueError(
                    f"Missing normalization status for {source_name}. "
                    "This is an internal error."
                )

            if status["raw_source"] is not None:
                if (
                    status["raw_source"] == "raw.X"
                    and adata.raw is not None
                    and not adata.raw.var_names.equals(adata.var_names)
                ):
                    raw_adata = adata.raw.to_adata()[:, adata.var_names].copy()
                    adata.raw = raw_adata

                adata = ensure_normalized(adata, copy=False)
                logger.info(f"  {source_name}: applied log1p(10k) normalization")
            else:
                if status["x_state"] == "log1p_10k":
                    logger.warning(
                        f"  {source_name}: no raw counts found; using existing "
                        "log1p(10k) values after subsetting. Totals may be <10k."
                    )
                else:
                    raise ValueError(
                        f"Cannot safely normalize {source_name}. "
                        f"Detected X state: {status['x_state']}. "
                        "Raw counts not found. Provide raw counts or log1p(10k) data."
                    )

            adatas[i] = adata

    # Memory check before concatenation
    available_gb = get_available_memory_gb()
    total_cells = sum(a.n_obs for a in adatas)
    estimated_gb = (total_cells * len(shared_genes) * 4) / (1024**3)

    logger.info(f"\nMemory check before concatenation:")
    logger.info(f"  Total cells: {total_cells:,}")
    if available_gb > 0:
        logger.info(f"  Available: {available_gb:.1f} GB")
        logger.info(f"  Estimated need: ~{estimated_gb:.1f} GB")

    # Concatenate
    logger.info("\nConcatenating references...")
    combined = sc.concat(adatas, join="inner", label="batch", index_unique="-")
    logger.info(f"  Combined: {combined.n_obs:,} cells × {combined.n_vars:,} genes")

    # Filter excluded labels (exact match only, no partial matching)
    if exclude_labels is None:
        exclude_labels = DEFAULT_EXCLUDE_LABELS

    labels_to_exclude = set(exclude_labels)
    label_values = combined.obs[output_column].astype(str)
    exclude_mask = label_values.isin(labels_to_exclude)
    n_excluded_labels = exclude_mask.sum()

    if n_excluded_labels > 0:
        excluded_counts = combined.obs.loc[exclude_mask, output_column].value_counts()
        logger.info(f"\nFiltering excluded labels:")
        for label, count in excluded_counts.items():
            logger.info(f"  Removing '{label}': {count:,} cells")
        combined = combined[~exclude_mask].copy()
        logger.info(f"  Remaining: {combined.n_obs:,} cells")

    # Filter low-count cell types
    if filter_min_cells and min_cells_per_type > 0:
        type_counts = combined.obs[output_column].value_counts()
        low_count_types = type_counts[type_counts < min_cells_per_type].index.tolist()

        if low_count_types:
            low_count_mask = combined.obs[output_column].isin(low_count_types)
            n_removed = low_count_mask.sum()
            logger.info(f"\nFiltering low-count cell types (<{min_cells_per_type} cells):")
            logger.info(f"  Removing {len(low_count_types)} types, {n_removed:,} cells")
            for ct in low_count_types[:10]:
                logger.info(f"    {ct}: {type_counts[ct]} cells")
            if len(low_count_types) > 10:
                logger.info(f"    ... and {len(low_count_types) - 10} more types")
            combined = combined[~low_count_mask].copy()
            logger.info(f"  Remaining: {combined.n_obs:,} cells, {combined.obs[output_column].nunique()} types")

    # Print cell type distribution
    logger.info(f"\n  Cell type distribution:")
    ct_counts = combined.obs[output_column].value_counts()
    for ct, count in ct_counts.head(10).items():
        logger.info(f"    {ct}: {count:,} cells")
    if len(ct_counts) > 10:
        logger.info(f"    ... and {len(ct_counts) - 10} more types")

    # Store validation results in uns if validation was performed
    # Convert to JSON string for h5ad serialization (lists of dicts not supported)
    if validate_labels and validation_results:
        combined.uns["validation_results"] = json.dumps(validation_results)

    logger.info(
        f"\nCombined reference ready. Call subsample_balanced() for "
        f"source-aware balancing before training."
    )

    return combined


def get_panel_genes(adata: ad.AnnData) -> List[str]:
    """
    Extract panel genes from AnnData (e.g., Xenium spatial data).

    Parameters
    ----------
    adata : AnnData
        AnnData object (typically from spatial platform).

    Returns
    -------
    List[str]
        List of gene names (panel genes).
    """
    return list(adata.var_names)


# ============================================================================
# CellTypist Training
# ============================================================================

def _save_model_metadata(
    metadata_path: Path,
    model,
    adata: ad.AnnData,
    label_column: str,
    training_params: Dict[str, Any],
) -> None:
    """
    Save model training metadata for reproducibility.

    Parameters
    ----------
    metadata_path : Path
        Path to save JSON metadata.
    model
        Trained CellTypist model.
    adata : AnnData
        Training data.
    label_column : str
        Cell type label column.
    training_params : Dict
        Training parameters used.
    """
    try:
        import celltypist
        celltypist_version = celltypist.__version__
    except Exception:
        celltypist_version = "unknown"

    try:
        import spatialcore
        spatialcore_version = spatialcore.__version__
    except Exception:
        spatialcore_version = "unknown"

    # Cell type counts
    cell_type_counts = adata.obs[label_column].value_counts().to_dict()

    # Reference sources if tracked
    reference_sources = []
    if "reference_source" in adata.obs.columns:
        for source in adata.obs["reference_source"].unique():
            source_mask = adata.obs["reference_source"] == source
            reference_sources.append({
                "name": source,
                "n_cells_used": int(source_mask.sum()),
            })

    metadata = {
        "model_name": metadata_path.stem.replace("_metadata", ""),
        "created_at": datetime.now().isoformat(),
        "spatialcore_version": spatialcore_version,
        "celltypist_version": celltypist_version,
        "training": {
            "n_cells": int(adata.n_obs),
            "n_genes": int(len(model.features)),
            "n_cell_types": int(len(model.cell_types)),
            "label_column": label_column,
            **training_params,
        },
        "references": reference_sources,
        "panel_genes": {
            "n_genes": int(len(model.features)),
            "genes": list(model.features)[:50],  # First 50 for preview
            "genes_truncated": len(model.features) > 50,
        },
        "cell_type_summary": {
            str(k): int(v) for k, v in cell_type_counts.items()
        },
    }

    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2)


def train_celltypist_model(
    adata: ad.AnnData,
    label_column: str = "unified_cell_type",
    model_name: str = "custom_model",
    output_path: Optional[Union[str, Path]] = None,
    use_SGD: bool = True,
    mini_batch: bool = True,
    balance_cell_type: bool = True,
    feature_selection: bool = False,
    n_jobs: int = -1,
    max_iter: int = 100,
    epochs: int = 10,
    batch_size: int = 1000,
    batch_number: int = 200,
) -> Dict[str, Any]:
    """
    Train a custom CellTypist logistic regression model.

    Parameters
    ----------
    adata : AnnData
        Reference data (already subset to panel genes, normalized).
    label_column : str, default "unified_cell_type"
        Cell type label column in adata.obs.
    model_name : str, default "custom_model"
        Name for the trained model.
    output_path : str or Path, optional
        Where to save .pkl file. If None, doesn't save.
    use_SGD : bool, default True
        Use stochastic gradient descent (faster for large data).
    mini_batch : bool, default True
        Use mini-batch training (recommended for large datasets).
    balance_cell_type : bool, default True
        Balance rare cell types in mini-batches.
    feature_selection : bool, default False
        Perform feature selection (False = use all genes).
    n_jobs : int, default -1
        Parallel jobs (-1 = all cores).
    max_iter : int, default 100
        Maximum training iterations (for non-mini-batch).
    epochs : int, default 10
        Training epochs (for mini-batch).
    batch_size : int, default 1000
        Cells per batch (for mini-batch).
    batch_number : int, default 200
        Batches per epoch (for mini-batch).

    Returns
    -------
    Dict[str, Any]
        Model metadata including:
        - model_path: Path to saved model (if output_path provided)
        - n_cells_trained: Number of cells used for training
        - n_genes: Number of features/genes
        - n_cell_types: Number of cell types
        - cell_types: List of cell type names
        - model: The trained CellTypist model object

    Notes
    -----
    For imbalanced data (e.g., 192k T cells vs 277 mast cells), use
    `balance_cell_type=True` with `mini_batch=True` to ensure rare
    cell types are adequately represented in training batches.

    Examples
    --------
    >>> from spatialcore.annotation import train_celltypist_model
    >>> result = train_celltypist_model(
    ...     combined_adata,
    ...     label_column="unified_cell_type",
    ...     output_path="./models/liver_xenium_v1.pkl",
    ...     mini_batch=True,
    ...     balance_cell_type=True,
    ... )
    >>> print(f"Trained on {result['n_cells_trained']:,} cells")
    >>> print(f"Cell types: {result['n_cell_types']}")
    """
    try:
        import celltypist
    except ImportError:
        raise ImportError(
            "celltypist is required for model training. "
            "Install with: pip install celltypist"
        )

    # Validate label column
    if label_column not in adata.obs.columns:
        available = list(adata.obs.columns)
        raise ValueError(
            f"Label column '{label_column}' not found. Available: {available}"
        )

    # Validate batch_size for mini-batch training
    n_cells = adata.n_obs
    if mini_batch and n_cells <= batch_size:
        raise ValueError(
            f"Dataset has {n_cells:,} cells but batch_size={batch_size}. "
            f"Either reduce batch_size (e.g., batch_size={max(50, n_cells // 2)}) "
            f"or set mini_batch=False for full-batch training."
        )

    # Log training parameters
    logger.info("Training CellTypist model...")
    logger.info(f"  Cells: {n_cells:,}")
    logger.info(f"  Genes: {adata.n_vars:,}")
    logger.info(f"  Cell types: {adata.obs[label_column].nunique()}")
    logger.info(f"  Mini-batch: {mini_batch}")
    logger.info(f"  Balance cell types: {balance_cell_type}")

    # Train model
    if mini_batch:
        model = celltypist.train(
            adata,
            labels=label_column,
            check_expression=False,  # Already validated
            use_SGD=use_SGD,
            mini_batch=True,
            balance_cell_type=balance_cell_type,
            feature_selection=feature_selection,
            n_jobs=n_jobs,
            epochs=epochs,
            batch_size=batch_size,
            batch_number=batch_number,
        )
    else:
        model = celltypist.train(
            adata,
            labels=label_column,
            check_expression=False,
            use_SGD=use_SGD,
            feature_selection=feature_selection,
            n_jobs=n_jobs,
            max_iter=max_iter,
        )

    # Save model if output path provided
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        model.write(str(output_path))
        logger.info(f"  Saved model to: {output_path}")

        # Save metadata JSON for reproducibility
        metadata_path = output_path.with_suffix("").with_name(
            output_path.stem + "_celltypist.json"
        )
        _save_model_metadata(
            metadata_path=metadata_path,
            model=model,
            adata=adata,
            label_column=label_column,
            training_params={
                "use_SGD": use_SGD,
                "mini_batch": mini_batch,
                "balance_cell_type": balance_cell_type,
                "feature_selection": feature_selection,
                "n_jobs": n_jobs,
                "max_iter": max_iter,
                "epochs": epochs,
                "batch_size": batch_size,
                "batch_number": batch_number,
            },
        )
        logger.info(f"  Saved metadata to: {metadata_path}")

    return {
        "model_path": str(output_path) if output_path else None,
        "metadata_path": str(metadata_path) if output_path else None,
        "n_cells_trained": adata.n_obs,
        "n_genes": len(model.features),
        "n_cell_types": len(model.cell_types),
        "cell_types": list(model.cell_types),
        "model": model,
    }


def get_model_gene_overlap(
    model_path: Union[str, Path],
    query_genes: List[str],
) -> Dict[str, Any]:
    """
    Calculate gene overlap between a CellTypist model and query data.

    Parameters
    ----------
    model_path : str or Path
        Path to CellTypist model (.pkl file).
    query_genes : List[str]
        Gene names from query data (e.g., Xenium panel).

    Returns
    -------
    Dict[str, Any]
        - n_model_genes: Total genes in model
        - n_query_genes: Total genes in query
        - n_overlap: Number of overlapping genes
        - overlap_pct: Percentage of model genes in query
        - overlapping_genes: List of overlapping gene names
        - missing_genes: List of model genes missing from query

    Examples
    --------
    >>> from spatialcore.annotation import get_model_gene_overlap
    >>> overlap = get_model_gene_overlap(
    ...     "Healthy_Human_Liver.pkl",
    ...     list(xenium_adata.var_names)
    ... )
    >>> print(f"Gene overlap: {overlap['overlap_pct']:.1f}%")
    """
    try:
        from celltypist import models
    except ImportError:
        raise ImportError("celltypist is required. Install with: pip install celltypist")

    model = models.Model.load(str(model_path))
    model_genes = set(model.features)
    query_genes_set = set(query_genes)

    overlap = model_genes & query_genes_set
    missing = model_genes - query_genes_set

    return {
        "n_model_genes": len(model_genes),
        "n_query_genes": len(query_genes_set),
        "n_overlap": len(overlap),
        "overlap_pct": 100 * len(overlap) / len(model_genes) if model_genes else 0,
        "overlapping_genes": sorted(overlap),
        "missing_genes": sorted(missing),
    }


def get_training_summary(combined_adata: ad.AnnData, label_column: str) -> pd.DataFrame:
    """
    Get summary of cell type distribution for training data.

    Parameters
    ----------
    combined_adata : AnnData
        Combined reference data.
    label_column : str
        Cell type label column.

    Returns
    -------
    pd.DataFrame
        Summary with columns: cell_type, n_cells, pct_total.
    """
    counts = combined_adata.obs[label_column].value_counts()
    df = pd.DataFrame({
        "cell_type": counts.index,
        "n_cells": counts.values,
        "pct_total": 100 * counts.values / combined_adata.n_obs,
    })
    return df


# ============================================================================
# Color Palettes
# ============================================================================

# High-contrast palette optimized for dark backgrounds
# Designed for spatial maps and UMAP visualizations
HIGH_CONTRAST_PALETTE = [
    "#F5252F",  # Red
    "#FB3FFC",  # Magenta
    "#00FFFF",  # Cyan
    "#33FF33",  # Green
    "#FFB300",  # Amber
    "#9966FF",  # Purple
    "#FF6B6B",  # Coral
    "#3784FE",  # Blue
    "#FF8000",  # Orange
    "#66CCCC",  # Teal
    "#CC66FF",  # Lavender
    "#99FF99",  # Light green
    "#FF6699",  # Pink
    "#E3E1E3",  # Light gray
    "#FFB3CC",  # Light pink
    "#E6E680",  # Yellow-green
    "#CC9966",  # Tan
    "#8080FF",  # Periwinkle
    "#FF9999",  # Salmon
    "#66FF99",  # Mint
]


def generate_color_scheme(
    cell_types: List[str],
    custom_colors: Optional[Dict[str, str]] = None,
    palette: Optional[List[str]] = None,
) -> Dict[str, str]:
    """
    Generate deterministic color mapping for cell types.

    Creates a mapping from cell type names to hex colors, using
    custom overrides if provided, then filling remaining types
    from the palette.

    Parameters
    ----------
    cell_types : List[str]
        List of cell type names.
    custom_colors : Dict[str, str], optional
        Custom color overrides. Keys are cell type names, values
        are hex color codes (e.g., "#FF0000").
    palette : List[str], optional
        Color palette to use. If None, uses HIGH_CONTRAST_PALETTE.

    Returns
    -------
    Dict[str, str]
        Mapping from cell type to hex color.

    Notes
    -----
    Colors are assigned deterministically based on sorted cell type names,
    ensuring consistent colors across runs.

    Examples
    --------
    >>> from spatialcore.annotation.training import generate_color_scheme
    >>> colors = generate_color_scheme(
    ...     ["T cell", "B cell", "Macrophage"],
    ...     custom_colors={"T cell": "#FF0000"},
    ... )
    >>> print(colors)
    """
    if palette is None:
        palette = HIGH_CONTRAST_PALETTE

    if custom_colors is None:
        custom_colors = {}

    color_scheme = {}
    palette_idx = 0

    # Sort for deterministic ordering
    for cell_type in sorted(cell_types):
        if cell_type in custom_colors:
            color_scheme[cell_type] = custom_colors[cell_type]
        else:
            color_scheme[cell_type] = palette[palette_idx % len(palette)]
            palette_idx += 1

    return color_scheme


# ============================================================================
# Model Artifact Management
# ============================================================================

def save_model_artifacts(
    model,
    output_dir: Union[str, Path],
    model_name: str,
    training_metadata: Optional[Dict[str, Any]] = None,
    custom_colors: Optional[Dict[str, str]] = None,
) -> Dict[str, Path]:
    """
    Save model with metadata and color scheme.

    Creates a complete model artifact package with:
    - Model pickle file (.pkl)
    - Training metadata JSON
    - Color scheme JSON for visualization

    Parameters
    ----------
    model
        Trained CellTypist model.
    output_dir : str or Path
        Directory to save artifacts.
    model_name : str
        Base name for output files.
    training_metadata : Dict[str, Any], optional
        Additional metadata to include (e.g., training parameters,
        reference sources).
    custom_colors : Dict[str, str], optional
        Custom color overrides for specific cell types.

    Returns
    -------
    Dict[str, Path]
        Paths to saved files:
        - model_path: Path to .pkl model file
        - metadata_path: Path to metadata JSON
        - colors_path: Path to color scheme JSON

    Notes
    -----
    The color scheme is saved separately so it can be loaded by
    visualization functions without loading the full model.

    File naming:
    - {model_name}.pkl
    - {model_name}_celltypist.json
    - {model_name}_colors.json

    Examples
    --------
    >>> from spatialcore.annotation.training import save_model_artifacts
    >>> result = train_celltypist_model(adata, label_column="cell_type")
    >>> paths = save_model_artifacts(
    ...     result["model"],
    ...     output_dir="./models",
    ...     model_name="liver_v1",
    ...     training_metadata={"reference": "cellxgene"},
    ... )
    >>> print(f"Model saved to: {paths['model_path']}")
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Save model
    model_path = output_dir / f"{model_name}.pkl"
    model.write(str(model_path))
    logger.info(f"Saved model to: {model_path}")

    # Generate color scheme
    cell_types = list(model.cell_types)
    colors = generate_color_scheme(cell_types, custom_colors)

    # Save color scheme
    colors_path = output_dir / f"{model_name}_colors.json"
    with open(colors_path, "w") as f:
        json.dump(colors, f, indent=2)
    logger.info(f"Saved color scheme to: {colors_path}")

    # Build metadata
    try:
        import celltypist
        celltypist_version = celltypist.__version__
    except Exception:
        celltypist_version = "unknown"

    try:
        import spatialcore
        spatialcore_version = spatialcore.__version__
    except Exception:
        spatialcore_version = "unknown"

    metadata = {
        "model_name": model_name,
        "created_at": datetime.now().isoformat(),
        "spatialcore_version": spatialcore_version,
        "celltypist_version": celltypist_version,
        "n_genes": len(model.features),
        "n_cell_types": len(model.cell_types),
        "cell_types": cell_types,
        "genes": list(model.features),
    }

    if training_metadata:
        metadata["training"] = training_metadata

    # Save metadata
    metadata_path = output_dir / f"{model_name}_celltypist.json"
    with open(metadata_path, "w") as f:
        json.dump(metadata, f, indent=2)
    logger.info(f"Saved metadata to: {metadata_path}")

    return {
        "model_path": model_path,
        "metadata_path": metadata_path,
        "colors_path": colors_path,
    }


# ============================================================================
# Balanced Subsampling
# ============================================================================

def _load_target_proportions(
    target_proportions: Union[Dict[str, float], str, Path, None]
) -> Optional[Dict[str, float]]:
    """
    Load and validate target proportions from dict, JSON, or CSV.

    Parameters
    ----------
    target_proportions : dict, str, Path, or None
        - Dict: Used directly
        - str/Path ending in .json: Load as JSON
        - str/Path ending in .csv: Load as CSV with columns (cell_type, proportion)
        - None: Return None

    Returns
    -------
    dict or None
        Validated proportions dict, or None if input was None.

    Raises
    ------
    ValueError
        If proportions are invalid (negative, >1.0, or wrong format).
    """
    if target_proportions is None:
        return None

    # Load from file if path provided
    if isinstance(target_proportions, (str, Path)):
        path = Path(target_proportions)

        if not path.exists():
            raise ValueError(f"Target proportions file not found: {path}")

        if path.suffix.lower() == ".json":
            with open(path, "r") as f:
                props = json.load(f)
        elif path.suffix.lower() == ".csv":
            df = pd.read_csv(path)
            if "cell_type" not in df.columns or "proportion" not in df.columns:
                raise ValueError(
                    f"CSV must have 'cell_type' and 'proportion' columns. "
                    f"Found: {list(df.columns)}"
                )
            props = dict(zip(df["cell_type"], df["proportion"]))
        else:
            raise ValueError(
                f"Unsupported file format: {path.suffix}. Use .json or .csv"
            )
    else:
        props = target_proportions

    # Validate proportions
    if not isinstance(props, dict):
        raise ValueError(
            f"target_proportions must be a dict, got {type(props).__name__}"
        )

    for cell_type, proportion in props.items():
        if not isinstance(proportion, (int, float)):
            raise ValueError(
                f"Invalid proportion for '{cell_type}': {proportion}. "
                f"Must be a number."
            )
        if proportion < 0 or proportion > 1.0:
            raise ValueError(
                f"Invalid proportion for '{cell_type}': {proportion}. "
                f"Must be between 0.0 and 1.0."
            )

    return props


def _resolve_target_totals(
    type_counts: pd.Series,
    min_cells_per_type: int,
    max_cells_per_type: int,
    props: Optional[Dict[str, float]],
) -> Dict[str, int]:
    """
    Resolve per-type target counts, honoring target proportions against the
    final output size.
    """
    if not props:
        targets: Dict[str, int] = {}
        for cell_type, n_available in type_counts.items():
            if n_available <= min_cells_per_type:
                targets[cell_type] = int(n_available)
            else:
                targets[cell_type] = int(min(max_cells_per_type, n_available))
        return targets

    missing = sorted(set(props) - set(type_counts.index))
    if missing:
        raise ValueError(
            "target_proportions include cell types not found in data: "
            + ", ".join(missing)
        )

    sum_props = float(sum(props.values()))
    eps = 1e-6
    non_prop_types = [ct for ct in type_counts.index if ct not in props]

    if sum_props > 1.0 + eps:
        raise ValueError(
            f"target_proportions sum to {sum_props:.4f}, must be <= 1.0"
        )
    if sum_props >= 1.0 - eps and non_prop_types:
        raise ValueError(
            "target_proportions sum to 1.0 but there are cell types without "
            "target proportions. Provide proportions for all types or "
            "reduce the total."
        )

    fixed_counts: Dict[str, int] = {}
    variable_props: Dict[str, int] = {}

    for cell_type, n_available in type_counts.items():
        if cell_type in props and n_available > min_cells_per_type:
            variable_props[cell_type] = int(n_available)
        else:
            if n_available <= min_cells_per_type:
                fixed_counts[cell_type] = int(n_available)
            else:
                fixed_counts[cell_type] = int(min(max_cells_per_type, n_available))

    fixed_total = sum(fixed_counts.values())

    if not variable_props:
        return fixed_counts

    sum_props_variable = sum(props[ct] for ct in variable_props)
    if sum_props_variable >= 1.0 - eps and fixed_total > 0:
        raise ValueError(
            "target_proportions leave no room for fixed counts. "
            "Reduce target_proportions or min_cells_per_type."
        )

    total = fixed_total
    if sum_props_variable > eps and (1.0 - sum_props_variable) > eps:
        total = int(round(fixed_total / (1.0 - sum_props_variable)))
    total = max(total, fixed_total)

    targets_var: Dict[str, int] = {}
    for _ in range(50):
        targets_var = {}
        for cell_type, n_available in variable_props.items():
            desired = int(props[cell_type] * total)
            target = max(min_cells_per_type, desired)
            target = min(target, n_available)
            targets_var[cell_type] = int(target)
        new_total = fixed_total + sum(targets_var.values())
        if new_total == total:
            break
        total = new_total
    else:
        raise RuntimeError(
            "Failed to resolve target_proportions. "
            "Check target_proportions and cell counts."
        )

    targets = dict(fixed_counts)
    targets.update(targets_var)

    total = sum(targets.values())
    for cell_type, prop in props.items():
        n_available = int(type_counts[cell_type])
        desired = int(prop * total)
        target = targets[cell_type]
        if n_available <= min_cells_per_type:
            logger.warning(
                f"Target proportion for '{cell_type}' cannot be met: "
                f"only {n_available} cells available (min_cells_per_type="
                f"{min_cells_per_type})."
            )
            continue
        if desired < min_cells_per_type and n_available >= min_cells_per_type:
            logger.warning(
                f"Target proportion for '{cell_type}' below min_cells_per_type; "
                f"using floor {min_cells_per_type} instead of {desired}."
            )
        elif desired > n_available:
            logger.warning(
                f"Target proportion for '{cell_type}' exceeds availability; "
                f"capping at {n_available} instead of {desired}."
            )

    return targets


def subsample_balanced(
    adata: ad.AnnData,
    label_column: str,
    max_cells_per_type: int = 5000,
    min_cells_per_type: int = 50,
    source_column: Optional[str] = "reference_source",
    source_balance: str = "proportional",
    min_cells_per_source: int = 50,
    group_by_column: Optional[str] = None,
    target_proportions: Optional[Union[Dict[str, float], str, Path]] = None,
    random_state: int = 42,
    copy: bool = True,
) -> ad.AnnData:
    """
    Source-aware balanced subsampling for multi-reference training data.

    When combining multiple scRNA-seq references, this function ensures each
    cell type is sampled proportionally from all sources that contain it.
    This prevents the model from learning source-specific artifacts.

    Parameters
    ----------
    adata : AnnData
        Combined reference data (output of combine_references()).
    label_column : str
        Column in adata.obs containing cell type labels (for logging/display).
    max_cells_per_type : int, default 5000
        Maximum cells per cell type in output.
    min_cells_per_type : int, default 50
        Minimum cells required to keep a cell type. Types with fewer
        cells are removed before balancing.
    source_column : str, optional, default "reference_source"
        Column identifying which reference each cell came from.
        Set to None to disable source-aware balancing (simple capping).
    source_balance : {"proportional", "equal"}, default "proportional"
        How to distribute sampling across sources:
        - "proportional": Draw from each source proportionally to its
          contribution for that cell type. (RECOMMENDED)
        - "equal": Draw equally from each source (up to available).
    min_cells_per_source : int, default 50
        Minimum cells to draw from a source that has the cell type.
        Ensures rare sources still contribute to the model.
    group_by_column : str, optional
        If provided, group cells by values in this column instead of label_column.
        This enables semantic grouping where different text labels map to the
        same identity. For example, using "cell_type_ontology_term_id":
        - "CD4+ T cells" and "CD4-positive, alpha-beta T cell" -> CL:0000624
        - Both are grouped together for balancing purposes.
        If None, groups by label_column text (current behavior).
    target_proportions : dict, str, or Path, optional
        Target proportions for specific cell types. Accepts:
        - Dict mapping cell type names to proportions (0.0-1.0)
        - Path to JSON file: ``{"NK cell": 0.0025, "plasma cell": 0.001}``
        - Path to CSV file with columns: ``cell_type``, ``proportion``

        For cell types in this mapping, the target count is calculated as:
        ``proportion x total_output_cells`` instead of using max_cells_per_type.

        Cell types NOT in the mapping use normal max_cells_per_type capping.

        This is essential for handling pure/enriched references (e.g., FACS-sorted
        cells) where a cell type exists only in the enriched source and would
        otherwise dominate training.
    random_state : int, default 42
        Random seed for reproducibility.
    copy : bool, default True
        Return a copy or modify in-place.

    Returns
    -------
    AnnData
        Subsampled data with source-balanced cell type representation.

    Raises
    ------
    ValueError
        If label_column not found in adata.obs.
        If source_column specified but not found in adata.obs.
        If source_balance is not "proportional" or "equal".

    Examples
    --------
    >>> from spatialcore.annotation import combine_references, subsample_balanced
    >>> # Step 1: Combine (no balancing)
    >>> combined = combine_references(
    ...     reference_paths=["study1.h5ad", "study2.h5ad"],
    ...     label_columns=["cell_type", "cell_type"],
    ...     target_genes=panel_genes,
    ... )
    >>> # Step 2: Source-aware balancing
    >>> balanced = subsample_balanced(
    ...     combined,
    ...     label_column="unified_cell_type",
    ...     max_cells_per_type=5000,
    ...     source_balance="proportional",
    ... )

    Notes
    -----
    **Why source-aware balancing?**

    When combining Study 1 (30K Macrophages) + Study 2 (5K Macrophages) and
    capping at 10K, naive capping takes mostly from Study 1. The model learns
    Study 1's version of Macrophage, not a consensus.

    Source-aware balancing draws proportionally: ~8.5K from Study 1, ~1.5K
    from Study 2, ensuring both studies contribute to each shared cell type.

    **Handling FACS-Enriched / Pure Cell Type References**

    When combining tissue-derived scRNA-seq (natural cell type proportions)
    with FACS-sorted or enriched populations (e.g., 100% pure T cells, sorted
    NK cells), special considerations apply:

    - **source_balance="proportional"** (default): Best when all references
      are tissue-derived with natural proportions. Each source contributes
      proportionally to its cell count for each type.

    - **source_balance="equal"**: Recommended when combining tissue-derived
      data with FACS-enriched pure populations. Forces equal contribution
      from each source, preventing pure populations from overwhelming the
      model's concept of that cell type.

    Example scenario:
    - Reference A (tissue): 50K cells with natural proportions (5% T cells = 2.5K)
    - Reference B (FACS): 100K pure T cells (100% T cells)

    With "proportional": T cells would be 97.5% from FACS, 2.5% from tissue.
    With "equal": T cells would be 50% from each, learning a consensus.

    **Future Improvements**

    - Automatic detection of enriched/pure cell type sources based on
      cell type distribution entropy
    - Per-cell-type source_balance overrides (e.g., use "equal" only for
      T cells that appear in FACS source)
    - Integration with Harmony batch correction for multi-source training via PCAs and predicition on spatial data with >5000 genes
    """
    if copy:
        adata = adata.copy()

    rng = np.random.default_rng(random_state)

    # Load and validate target proportions
    props = _load_target_proportions(target_proportions)
    if props:
        logger.info(f"Using target proportions for {len(props)} cell type(s)")
        for ct, prop in props.items():
            logger.debug(f"  {ct}: {prop:.4f} ({prop*100:.2f}%)")

    # Validate label_column
    if label_column not in adata.obs.columns:
        raise ValueError(
            f"Label column '{label_column}' not found. "
            f"Available: {list(adata.obs.columns)}"
        )

    # Determine grouping column: use group_by_column if provided, else label_column
    if group_by_column is not None:
        if group_by_column not in adata.obs.columns:
            raise ValueError(
                f"Group-by column '{group_by_column}' not found. "
                f"Available: {list(adata.obs.columns)}"
            )
        # Use group_by_column for grouping (e.g., CL IDs)
        cell_types = adata.obs[group_by_column].astype(str)
        logger.info(
            f"Grouping by '{group_by_column}' instead of '{label_column}' "
            f"(semantic grouping enabled)"
        )
    else:
        cell_types = adata.obs[label_column].astype(str)

    type_counts = cell_types.value_counts()
    if min_cells_per_type > 0:
        low_count_types = type_counts[type_counts < min_cells_per_type].index.tolist()
        if low_count_types:
            n_removed = int(type_counts[type_counts < min_cells_per_type].sum())
            logger.info(
                f"\nFiltering low-count cell types (<{min_cells_per_type} cells) before balancing:"
            )
            logger.info(f"  Removing {len(low_count_types)} types, {n_removed:,} cells")
            for ct in low_count_types[:10]:
                logger.info(f"    {ct}: {type_counts[ct]} cells")
            if len(low_count_types) > 10:
                logger.info(f"    ... and {len(low_count_types) - 10} more types")

            keep_mask = ~cell_types.isin(low_count_types)
            adata = adata[keep_mask].copy()

            if group_by_column is not None:
                cell_types = adata.obs[group_by_column].astype(str)
            else:
                cell_types = adata.obs[label_column].astype(str)

            if props:
                dropped = sorted(set(props) & set(low_count_types))
                if dropped:
                    for ct in dropped:
                        props.pop(ct, None)
                    logger.warning(
                        "Dropping target_proportions for low-count types: %s",
                        ", ".join(dropped),
                    )

            type_counts = cell_types.value_counts()

    unique_types = cell_types.unique()
    target_totals = _resolve_target_totals(
        type_counts=type_counts,
        min_cells_per_type=min_cells_per_type,
        max_cells_per_type=max_cells_per_type,
        props=props,
    )

    # =========================================================================
    # Source-unaware mode (simple capping)
    # =========================================================================
    if source_column is None:
        logger.info(
            f"Subsampling {adata.n_obs:,} cells "
            f"(cap={max_cells_per_type:,}, no source balancing)"
        )
        return _subsample_simple_cap(
            adata, cell_types, unique_types,
            max_cells_per_type, min_cells_per_type, rng,
            target_totals=target_totals,
        )

    # =========================================================================
    # Source-aware mode (nested balancing)
    # =========================================================================
    if source_column not in adata.obs.columns:
        raise ValueError(
            f"Source column '{source_column}' not found. "
            f"Available: {list(adata.obs.columns)}. "
            f"Set source_column=None to disable source-aware balancing."
        )

    if source_balance not in ("proportional", "equal"):
        raise ValueError(
            f"Invalid source_balance: '{source_balance}'. "
            f"Must be 'proportional' or 'equal'."
        )

    sources = adata.obs[source_column].astype(str)
    unique_sources = sources.unique()
    n_sources = len(unique_sources)

    logger.info(
        f"Subsampling {adata.n_obs:,} cells "
        f"(source-aware, {n_sources} sources, {source_balance} balance, "
        f"max={max_cells_per_type:,}/type)"
    )

    selected_indices = []

    for cell_type in unique_types:
        type_mask = cell_types == cell_type
        type_indices = np.where(type_mask)[0]
        n_available = len(type_indices)

        # Skip empty types
        if n_available == 0:
            continue

        target_total = target_totals[cell_type]
        if target_total >= n_available:
            selected_indices.extend(type_indices)
            logger.debug(
                f"  {cell_type}: keeping all {n_available} (target {target_total})"
            )
            continue

        # Identify which sources have this cell type
        type_sources = sources.iloc[type_indices]
        sources_with_type = type_sources.unique()
        n_sources_with_type = len(sources_with_type)

        # Warn if only one source has this type (when multiple sources exist)
        if n_sources_with_type == 1 and n_sources > 1:
            logger.warning(
                f"  {cell_type}: only in '{sources_with_type[0]}' "
                f"(no cross-source balancing)"
            )

        # Calculate per-source targets
        source_targets = _calculate_source_targets(
            type_sources=type_sources,
            sources_with_type=sources_with_type,
            target_total=target_total,
            source_balance=source_balance,
            min_cells_per_source=min_cells_per_source,
        )

        # Sample from each source
        for source_name, (target, available) in source_targets.items():
            source_type_mask = (cell_types == cell_type) & (sources == source_name)
            source_indices = np.where(source_type_mask)[0]

            if target >= len(source_indices):
                selected_indices.extend(source_indices)
            else:
                sampled = rng.choice(source_indices, size=target, replace=False)
                selected_indices.extend(sampled)

        # Debug log
        total_sampled = sum(t for t, _ in source_targets.values())
        source_summary = ", ".join(f"{s}:{t}" for s, (t, _) in source_targets.items())
        logger.debug(f"  {cell_type}: {n_available} -> {total_sampled} [{source_summary}]")

    # Sort and subset
    selected_indices = sorted(set(selected_indices))
    adata_sub = adata[selected_indices].copy()

    # Log final summary
    new_counts = adata_sub.obs[label_column].value_counts()
    logger.info(
        f"Subsampled: {adata.n_obs:,} -> {adata_sub.n_obs:,} cells "
        f"({len(new_counts)} types)"
    )

    return adata_sub


def _calculate_source_targets(
    type_sources: pd.Series,
    sources_with_type: np.ndarray,
    target_total: int,
    source_balance: str,
    min_cells_per_source: int,
) -> Dict[str, Tuple[int, int]]:
    """
    Calculate how many cells to sample from each source for one cell type.

    Parameters
    ----------
    type_sources : pd.Series
        Series of source labels for cells of this type.
    sources_with_type : np.ndarray
        Unique sources that have this cell type.
    target_total : int
        Total cells to sample for this cell type.
    source_balance : str
        "proportional" or "equal".
    min_cells_per_source : int
        Minimum cells to draw from each source.

    Returns
    -------
    Dict[str, Tuple[int, int]]
        Mapping: source_name -> (target_count, available_count)
    """
    source_counts = type_sources.value_counts().to_dict()
    total_available = sum(source_counts.values())
    n_sources = len(sources_with_type)

    targets = {}

    if source_balance == "proportional":
        # Draw proportionally to each source's contribution
        for source_name in sources_with_type:
            available = source_counts[source_name]
            proportion = available / total_available
            target = int(np.ceil(target_total * proportion))

            # Enforce minimum (if source has enough cells)
            if available >= min_cells_per_source:
                target = max(target, min_cells_per_source)

            # Can't exceed available
            target = min(target, available)
            targets[source_name] = (target, available)

    elif source_balance == "equal":
        # Draw equally from each source
        per_source = target_total // n_sources
        remainder = target_total % n_sources

        for i, source_name in enumerate(sorted(sources_with_type)):
            available = source_counts[source_name]
            target = per_source + (1 if i < remainder else 0)
            target = min(target, available)
            targets[source_name] = (target, available)

    # Redistribute shortfall (when some sources can't provide enough)
    total_targeted = sum(t for t, _ in targets.values())
    shortfall = target_total - total_targeted

    if shortfall > 0:
        for source_name in sources_with_type:
            if shortfall <= 0:
                break
            target, available = targets[source_name]
            capacity = available - target
            if capacity > 0:
                additional = min(capacity, shortfall)
                targets[source_name] = (target + additional, available)
                shortfall -= additional

    return targets


def _subsample_simple_cap(
    adata: ad.AnnData,
    cell_types: pd.Series,
    unique_types: np.ndarray,
    max_cells_per_type: int,
    min_cells_per_type: int,
    rng: np.random.Generator,
    target_totals: Optional[Dict[str, int]] = None,
) -> ad.AnnData:
    """
    Simple per-type capping without source awareness.

    Used when source_column=None is passed to subsample_balanced().
    """
    selected_indices = []

    for cell_type in unique_types:
        type_mask = cell_types == cell_type
        type_indices = np.where(type_mask)[0]
        n_available = len(type_indices)

        if n_available == 0:
            continue

        if target_totals is None:
            if n_available <= min_cells_per_type:
                selected_indices.extend(type_indices)
            elif n_available <= max_cells_per_type:
                selected_indices.extend(type_indices)
            else:
                sampled = rng.choice(type_indices, size=max_cells_per_type, replace=False)
                selected_indices.extend(sampled)
            continue

        target_total = target_totals[cell_type]
        if target_total >= n_available:
            selected_indices.extend(type_indices)
        elif target_total > 0:
            sampled = rng.choice(type_indices, size=target_total, replace=False)
            selected_indices.extend(sampled)

    selected_indices = sorted(selected_indices)
    return adata[selected_indices].copy()



# NOTE: harmonize_labels() function was removed in favor of using:
# 1. add_ontology_ids() to fill missing CL IDs
# 2. subsample_balanced(group_by_column="cell_type_ontology_term_id") for semantic grouping
# See pipeline.py for the new recommended workflow.
