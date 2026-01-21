"""
High-level cell typing pipeline for spatial transcriptomics.

This module provides consolidated entry points for the full annotation workflow:
- train_and_annotate(): Train custom model + annotate in one call
- TrainingConfig: YAML-serializable configuration for reproducibility

The pipeline integrates:
1. Reference loading and combination
2. Ontology ID filling (for semantic grouping)
3. Source-aware "Cap & Fill" balancing
4. CellTypist model training
5. Annotation with z-score confidence transformation
6. Ontology mapping
7. Validation plot generation

Column Naming (CellxGene Standard):
- cell_type: Final predicted cell type
- cell_type_confidence: Z-score transformed confidence
- cell_type_ontology_term_id: CL:XXXXX ontology ID
- cell_type_ontology_label: Canonical ontology name
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Literal, Optional, Union
import json
import gc

import anndata as ad
import numpy as np
import pandas as pd

from spatialcore.core.logging import get_logger

logger = get_logger(__name__)


# ============================================================================
# Configuration
# ============================================================================


@dataclass
class TrainingConfig:
    """
    YAML-serializable configuration for reproducible training + annotation.

    This configuration captures all parameters needed to reproduce a training
    run, from reference selection to confidence thresholds.

    Parameters
    ----------
    tissue : str, default "unknown"
        Tissue type for the model (used in model naming and selection).
    references : List[str]
        Paths to reference h5ad files. Required.
    label_columns : List[str]
        Cell type label column for each reference. Required: must be provided
        explicitly to avoid misinterpreting CL IDs as labels.
    balance_strategy : {"proportional", "equal"}, default "proportional"
        How to distribute sampling across sources when balancing.
    max_cells_per_type : int, default 10000
        Maximum cells per cell type after balancing.
    max_cells_per_ref : int, default 100000
        Maximum cells to load per reference (memory management).
    target_proportions : Dict[str, float], str, or Path, optional
        Expected biological proportions for specific cell types.
    confidence_threshold : float, default 0.8
        Threshold for marking low-confidence predictions as Unassigned.
    add_ontology : bool, default True
        Whether to add ontology IDs to predictions.
    generate_plots : bool, default True
        Whether to generate validation plots.

    Examples
    --------
    >>> from spatialcore.annotation.pipeline import TrainingConfig
    >>> config = TrainingConfig(
    ...     tissue="lung",
    ...     references=["ref1.h5ad", "ref2.h5ad"],
    ...     balance_strategy="proportional",
    ... )
    >>> config.to_yaml("training_config.yaml")
    """

    tissue: str = "unknown"
    references: List[str] = field(default_factory=list)
    label_columns: Optional[List[str]] = None
    balance_strategy: Literal["proportional", "equal"] = "proportional"
    max_cells_per_type: int = 5000
    max_cells_per_ref: int = 100000
    target_proportions: Optional[Union[Dict[str, float], str, Path]] = None
    confidence_threshold: float = 0.8
    add_ontology: bool = True
    generate_plots: bool = True

    @classmethod
    def from_yaml(cls, path: Union[str, Path]) -> "TrainingConfig":
        """Load configuration from YAML file."""
        try:
            import yaml
        except ImportError:
            raise ImportError("PyYAML is required for YAML config loading")

        path = Path(path)
        with open(path) as f:
            data = yaml.safe_load(f)

        return cls(**data)

    def to_yaml(self, path: Union[str, Path]) -> None:
        """Save configuration to YAML file."""
        try:
            import yaml
        except ImportError:
            raise ImportError("PyYAML is required for YAML config saving")

        path = Path(path)
        with open(path, "w") as f:
            yaml.dump(self.__dict__, f, default_flow_style=False)

        logger.info(f"Saved training config to: {path}")

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return self.__dict__.copy()

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "TrainingConfig":
        """Create from dictionary."""
        return cls(**data)

    def validate(self) -> None:
        """
        Validate required fields.

        Ensures label_columns is provided and non-empty to avoid
        misinterpreting ontology ID columns as raw labels.
        """
        if not self.label_columns:
            raise ValueError(
                "TrainingConfig.label_columns is required and cannot be empty. "
                "Provide one label column per reference."
            )


# ============================================================================
# High-Level Pipeline
# ============================================================================


def train_and_annotate(
    adata: ad.AnnData,
    references: List[Union[str, Path]],
    tissue: str = "unknown",
    label_columns: Optional[List[str]] = None,
    balance_strategy: Literal["proportional", "equal"] = "proportional",
    max_cells_per_type: int = 5000,
    max_cells_per_ref: int = 100000,
    target_proportions: Optional[Union[Dict[str, float], str, Path]] = None,
    confidence_threshold: float = 0.8,
    model_output: Optional[Union[str, Path]] = None,
    plot_output: Optional[Union[str, Path]] = None,
    add_ontology: bool = True,
    generate_plots: bool = True,
    copy: bool = False,
) -> ad.AnnData:
    """
    Full workflow: train custom model on references, then annotate spatial data.

    This is the core SpatialCore value proposition - NOT a thin wrapper around
    CellTypist. The function provides significant added value:

    1. **Panel-specific training** - Subsets references to spatial panel genes,
       achieving ~100% gene overlap vs ~5-9% with pre-trained models.
    2. **Source-aware balancing** - "Cap & Fill" strategy ensures all references
       contribute to each cell type, preventing source-specific biases.
    3. **CL ID-based grouping** - Groups semantically equivalent labels
       (e.g., "CD4+ T cells" and "CD4-positive, alpha-beta T cell") by
       their Cell Ontology ID for proper balancing.
    4. **Z-score confidence** - Transforms raw logistic regression scores to
       interpretable [0,1] confidence values that handle domain shift.
    5. **Multi-tier ontology mapping** - Maps predictions to Cell Ontology
       using pattern matching, exact match, token matching, and fuzzy overlap.
    6. **Automatic validation plots** - Generates DEG heatmap, 2D marker
       validation, confidence plots, and ontology mapping table.

    Pipeline Stages
    ---------------
    1. Get panel genes from spatial data (var_names)
    2. Load + combine references (Ensembl→HUGO normalization, log1p 10k)
    3. Fill missing ontology IDs (add_ontology_ids with skip_if_exists)
    4. Balance by CL ID (source-aware "Cap & Fill")
    5. Train CellTypist model
    6. Annotate spatial data
    7. Transform confidence (z-score)
    8. Add ontology IDs to predictions
    9. Generate validation plots

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics data to annotate. Can contain:
        - Raw counts in .X (will be normalized automatically)
        - Raw counts in .layers['counts'] or .raw.X
        - Already log1p(10k) normalized data in .X
        Gene names should be HUGO symbols.
    references : List[str or Path]
        Paths to reference h5ad files for training.
    tissue : str, default "unknown"
        Tissue type for model naming.
    label_columns : List[str]
        Cell type label column for each reference. Required: must be provided
        explicitly to avoid mis-mapping CL IDs as labels.
    balance_strategy : {"proportional", "equal"}, default "proportional"
        How to distribute sampling across sources:
        - "proportional": Sample proportionally to source contribution
        - "equal": Sample equally from each source
    max_cells_per_type : int, default 10000
        Maximum cells per cell type after balancing.
    max_cells_per_ref : int, default 100000
        Maximum cells to load per reference (memory management).
    target_proportions : Dict[str, float], str, or Path, optional
        Expected biological proportions for specific cell types. Use when
        combining tissue references with FACS-sorted or enriched populations.
        Accepts dict, JSON file path, or CSV file path. Keys are cell type
        labels, values are target proportions (0-1). Cell types not in the
        dict use default max_cells_per_type behavior.
    confidence_threshold : float, default 0.8
        Cells with confidence below this threshold are marked "Unassigned".
    model_output : str or Path, optional
        Path to save trained CellTypist model (.pkl). If None, model is
        not saved to disk.
    plot_output : str or Path, optional
        Directory to save validation plots. If None, uses current directory
        when generate_plots=True.
    add_ontology : bool, default True
        Whether to map predictions to Cell Ontology IDs.
    generate_plots : bool, default True
        Whether to generate validation plots (DEG heatmap, 2D validation,
        confidence plots, ontology mapping table).
    copy : bool, default False
        If True, return a copy of adata; otherwise modify in-place.

    Returns
    -------
    AnnData
        Annotated data with new columns (CellxGene standard names):
        - cell_type: Final cell type labels (ontology-mapped, confidence-filtered).
          Low-confidence cells are marked "Unassigned".
        - cell_type_confidence: Z-score transformed confidence [0, 1]
        - cell_type_ontology_term_id: CL:XXXXX (if add_ontology=True)
        - cell_type_ontology_label: Canonical ontology name for ALL cells
          (unfiltered, preserves predictions for low-confidence cells)

        And metadata in uns:
        - spatialcore_annotation: Dict with training parameters and stats

    Examples
    --------
    >>> from spatialcore.annotation.pipeline import train_and_annotate
    >>> import scanpy as sc
    >>>
    >>> # Load spatial data
    >>> adata = sc.read_h5ad("xenium_lung.h5ad")
    >>>
    >>> # Train and annotate
    >>> adata = train_and_annotate(
    ...     adata,
    ...     references=["hlca_core.h5ad", "tabula_sapiens_lung.h5ad"],
    ...     tissue="lung",
    ...     balance_strategy="proportional",
    ...     confidence_threshold=0.8,
    ...     plot_output="./qc_plots/",
    ... )
    >>>
    >>> # Check results
    >>> print(adata.obs["cell_type"].value_counts())

    See Also
    --------
    train_and_annotate_config : Config-driven version for reproducibility.
    annotate_celltypist : Lower-level annotation function.
    combine_references : Reference combination without annotation.
    """
    from spatialcore.annotation.training import (
        combine_references,
        get_panel_genes,
        subsample_balanced,
        train_celltypist_model,
        save_model_artifacts,
    )
    from spatialcore.annotation.annotate import annotate_celltypist
    from spatialcore.annotation.ontology import add_ontology_ids
    from spatialcore.plotting import generate_annotation_plots

    if copy:
        if adata.isbacked:
            raise ValueError(
                "Cannot use copy=True with backed AnnData. "
                "Either use copy=False or load data into memory first with adata.to_memory()"
            )
        adata = adata.copy()

    logger.info("=" * 60)
    logger.info("SpatialCore Cell Typing Pipeline")
    logger.info("=" * 60)

    # -------------------------------------------------------------------------
    # Stage 1: Get panel genes from spatial data
    # -------------------------------------------------------------------------
    logger.info("Stage 1: Extracting panel genes from spatial data...")
    panel_genes = get_panel_genes(adata)
    logger.info(f"  Panel genes: {len(panel_genes)}")

    # -------------------------------------------------------------------------
    # Stage 2: Load and combine references
    # -------------------------------------------------------------------------
    logger.info("Stage 2: Loading and combining references...")

    if label_columns is None:
        raise ValueError(
            "label_columns must be provided (one per reference). "
            "Auto-detection was removed to prevent misinterpreting CL ID "
            "columns as raw labels."
        )

    combined = combine_references(
        reference_paths=references,
        label_columns=label_columns,
        output_column="original_label",
        max_cells_per_ref=max_cells_per_ref,
        target_genes=panel_genes,
        normalize_data=True,
    )
    logger.info(f"  Combined: {combined.n_obs:,} cells, {combined.n_vars:,} genes")

    # -------------------------------------------------------------------------
    # Stage 3: Fill missing ontology IDs
    # -------------------------------------------------------------------------
    logger.info("Stage 3: Filling missing ontology IDs...")

    combined, _, ontology_result = add_ontology_ids(
        combined,
        source_col="original_label",
        target_col="cell_type_ontology_term_id",
        name_col="cell_type_ontology_label",
        skip_if_exists=True,  # Preserve CellxGene's native IDs
        copy=False,
    )

    # -------------------------------------------------------------------------
    # Stage 4: Balance by CL ID (source-aware)
    # -------------------------------------------------------------------------
    logger.info("Stage 4: Balancing training data (source-aware)...")

    balanced = subsample_balanced(
        combined,
        label_column="original_label",
        group_by_column="cell_type_ontology_term_id",  # Group by CL ID!
        source_column="reference_source",
        source_balance=balance_strategy,
        max_cells_per_type=max_cells_per_type,
        target_proportions=target_proportions,
        copy=True,
    )
    logger.info(f"  Balanced: {balanced.n_obs:,} cells")

    # Release combined reference data - no longer needed after balancing
    del combined
    gc.collect()

    # -------------------------------------------------------------------------
    # Stage 5: Train CellTypist model
    # -------------------------------------------------------------------------
    logger.info("Stage 5: Training CellTypist model...")

    # Train on canonical ontology names (human-readable, semantically grouped)
    # Balancing grouped by CL ID, so synonyms are consolidated under the same
    # canonical name in cell_type_ontology_label (e.g., "CD4+ T cell" and
    # "CD4-positive, alpha-beta T cell" both become the canonical CL name)
    training_result = train_celltypist_model(
        balanced,
        label_column="cell_type_ontology_label",  # Train on canonical names
        feature_selection=False,  # Use panel genes as-is
        n_jobs=-1,
    )
    model = training_result["model"]
    n_training_cells = balanced.n_obs

    # Save model to user-defined path (use temp dir if not specified)
    if model_output is None:
        import tempfile
        temp_dir = tempfile.mkdtemp(prefix="spatialcore_model_")
        model_path = Path(temp_dir) / "model"
        logger.info("  No model_output specified, using temporary directory")
    else:
        model_path = Path(model_output)
        model_path.parent.mkdir(parents=True, exist_ok=True)
    artifacts = save_model_artifacts(
        model,
        output_dir=model_path.parent,
        model_name=model_path.stem,
        training_metadata={
            "references": [str(p) for p in references],
            "n_cells": training_result["n_cells_trained"],
            "n_genes": training_result["n_genes"],
            "n_cell_types": training_result["n_cell_types"],
            "cell_types": training_result["cell_types"],
        },
    )
    # Use actual model path returned by save_model_artifacts (includes .pkl extension)
    actual_model_path = artifacts["model_path"]
    logger.info(f"  Model saved to: {actual_model_path}")

    # Release training data
    del balanced
    gc.collect()

    # -------------------------------------------------------------------------
    # Stage 6: Annotate spatial data
    # -------------------------------------------------------------------------
    logger.info("Stage 6: Annotating spatial data...")

    # Ensure spatial data is normalized (handles raw counts, pre-normalized, etc.)
    from spatialcore.core.utils import check_normalization_status
    from spatialcore.annotation.loading import ensure_normalized

    status = check_normalization_status(adata)
    if status["x_state"] != "log1p_10k":
        logger.info(f"  Normalizing spatial data (detected: {status['x_state']})")
        adata = ensure_normalized(adata, target_sum=1e4, copy=False)

    adata = annotate_celltypist(
        adata,
        custom_model_path=actual_model_path,
        confidence_transform="zscore",
        store_decision_scores=True,
        min_confidence=0.0,
        copy=False,
    )

    # -------------------------------------------------------------------------
    # Stage 7: Add ontology IDs to predictions
    # -------------------------------------------------------------------------
    if add_ontology:
        logger.info("Stage 7: Mapping predictions to Cell Ontology...")

        adata, _, _ = add_ontology_ids(
            adata,
            source_col="cell_type",
            target_col="cell_type_ontology_term_id",
            name_col="cell_type_ontology_label",
            skip_if_exists=False,  # Map all predictions
            copy=False,
        )

    # -------------------------------------------------------------------------
    # Stage 8: Generate validation plots
    # -------------------------------------------------------------------------
    if generate_plots:
        logger.info("Stage 8: Generating validation plots...")

        # Release training artifacts before memory-intensive plot generation
        del model
        del training_result
        gc.collect()

        try:
            output_dir = Path(plot_output) if plot_output else Path(".")
            output_dir.mkdir(parents=True, exist_ok=True)

            generate_annotation_plots(
                adata,
                label_column="cell_type_ontology_label",
                confidence_column="cell_type_confidence",
                output_dir=output_dir,
                prefix=f"{tissue}_celltyping",
                confidence_threshold=confidence_threshold,
                # Ontology columns for mapping table
                source_label_column="cell_type",
                ontology_name_column="cell_type_ontology_label",
                ontology_id_column="cell_type_ontology_term_id",
            )
        except Exception as exc:
            logger.warning(
                "Plot generation failed; continuing without plots: %s",
                exc,
                exc_info=True,
            )

    # -------------------------------------------------------------------------
    # Stage 9: Apply confidence threshold (after plots, so plots show all cells)
    # -------------------------------------------------------------------------
    if confidence_threshold > 0:
        conf = adata.obs["cell_type_confidence"].values
        low_conf_mask = conf < confidence_threshold
        n_low = low_conf_mask.sum()

        if n_low > 0:
            # Mark low-confidence cells as Unassigned in cell_type
            labels = adata.obs["cell_type"].astype(str).copy()
            labels[low_conf_mask] = "Unassigned"
            adata.obs["cell_type"] = pd.Categorical(labels)

            pct = 100 * n_low / adata.n_obs
            logger.info(
                f"Stage 9: Marked {n_low:,} cells ({pct:.1f}%) as Unassigned "
                f"(confidence < {confidence_threshold})"
            )

    # -------------------------------------------------------------------------
    # Store metadata
    # -------------------------------------------------------------------------
    adata.uns["spatialcore_annotation"] = {
        "tissue": tissue,
        "n_references": len(references),
        "references": [str(r) for r in references],
        "panel_genes": len(panel_genes),
        "training_cells": n_training_cells,
        "balance_strategy": balance_strategy,
        "max_cells_per_type": max_cells_per_type,
        "confidence_threshold": confidence_threshold,
        "n_cell_types": adata.obs["cell_type"].nunique(),
    }

    logger.info("=" * 60)
    logger.info("Pipeline complete!")
    logger.info(f"  Cell types: {adata.obs['cell_type'].nunique()}")
    logger.info(f"  Mean confidence: {adata.obs['cell_type_confidence'].mean():.3f}")
    logger.info("=" * 60)

    return adata


def train_and_annotate_config(
    adata: ad.AnnData,
    config: TrainingConfig,
    model_output: Optional[Union[str, Path]] = None,
    plot_output: Optional[Union[str, Path]] = None,
    copy: bool = False,
) -> ad.AnnData:
    """
    Config-driven training + annotation for reproducible workflows.

    This is a convenience wrapper around train_and_annotate() that accepts
    a TrainingConfig object instead of individual parameters.

    Parameters
    ----------
    adata : AnnData
        Spatial transcriptomics data to annotate.
    config : TrainingConfig
        Configuration object with training parameters.
    model_output : str or Path, optional
        Path to save trained model.
    plot_output : str or Path, optional
        Directory to save validation plots.
    copy : bool, default False
        If True, return a copy of adata.

    Returns
    -------
    AnnData
        Annotated data with cell type predictions.

    Examples
    --------
    >>> from spatialcore.annotation.pipeline import (
    ...     TrainingConfig,
    ...     train_and_annotate_config,
    ... )
    >>>
    >>> # Load config from YAML
    >>> config = TrainingConfig.from_yaml("training_config.yaml")
    >>>
    >>> # Run pipeline
    >>> adata = train_and_annotate_config(adata, config)
    """
    config.validate()
    return train_and_annotate(
        adata=adata,
        references=config.references,
        tissue=config.tissue,
        label_columns=config.label_columns,
        balance_strategy=config.balance_strategy,
        max_cells_per_type=config.max_cells_per_type,
        max_cells_per_ref=config.max_cells_per_ref,
        target_proportions=config.target_proportions,
        confidence_threshold=config.confidence_threshold,
        model_output=model_output,
        plot_output=plot_output,
        add_ontology=config.add_ontology,
        generate_plots=config.generate_plots,
        copy=copy,
    )


# ============================================================================
# Helper Functions
# ============================================================================


def _detect_label_columns(
    references: List[Union[str, Path]],
) -> List[str]:
    """
    Auto-detect cell type label columns in reference files.

    Searches for common column names: cell_type, celltype, cell_type_ontology_term_id,
    Cell_type, CellType, etc.

    Parameters
    ----------
    references : List[str or Path]
        Paths to reference h5ad files.

    Returns
    -------
    List[str]
        Detected label column for each reference.
    """
    common_columns = [
        "cell_type",
        "celltype",
        "cell_type_ontology_term_id",
        "Cell_type",
        "CellType",
        "cell_type_label",
        "annotation",
        "cluster",
        "leiden",
    ]

    label_columns = []

    for ref_path in references:
        ref_path = Path(ref_path)

        # Read just the obs columns (don't load full data)
        try:
            import h5py

            with h5py.File(ref_path, "r") as f:
                if "obs" in f:
                    if "__categories" in f["obs"]:
                        # Categorical columns
                        obs_cols = list(f["obs"]["__categories"].keys())
                    else:
                        obs_cols = [
                            k for k in f["obs"].keys()
                            if not k.startswith("_")
                        ]
                else:
                    obs_cols = []
        except Exception:
            # Fallback: load full file
            import anndata
            adata_ref = anndata.read_h5ad(ref_path, backed="r")
            obs_cols = list(adata_ref.obs.columns)
            adata_ref.file.close()

        # Find first matching column
        found = None
        for col in common_columns:
            if col in obs_cols:
                found = col
                break

        if found is None:
            raise ValueError(
                f"Could not auto-detect label column in {ref_path.name}. "
                f"Available columns: {obs_cols[:10]}... "
                f"Please provide label_columns explicitly."
            )

        label_columns.append(found)
        logger.debug(f"  {ref_path.name}: detected label column '{found}'")

    return label_columns
