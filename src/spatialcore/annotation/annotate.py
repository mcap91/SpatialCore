"""
CellTypist annotation wrapper for cell type annotation.

This module provides a convenience wrapper around CellTypist for:
1. Tissue-specific model selection
2. Ensemble annotation across multiple models
3. Gene overlap validation
4. Proper re-normalization for CellTypist compatibility

References:
    - CellTypist: https://www.celltypist.org/
    - Model documentation: See docs/CELLTYPIST_MODELS.md
"""

from pathlib import Path
from typing import Dict, List, Literal, Optional, Union, Any
from types import SimpleNamespace
import gc

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from spatialcore.core.logging import get_logger
from spatialcore.core.utils import check_normalization_status
from spatialcore.annotation.confidence import (
    extract_decision_scores,
    transform_confidence,
    ConfidenceMethod,
)

logger = get_logger(__name__)


# ============================================================================
# Tissue-Specific Model Presets
# ============================================================================

TISSUE_MODEL_PRESETS: Dict[str, List[str]] = {
    # General (use for unknown tissues)
    "unknown": [
        "Immune_All_Low.pkl",
        "Pan_Fetal_Human.pkl",
    ],
    # Digestive system
    "colon": [
        "Immune_All_Low.pkl",
        "Pan_Fetal_Human.pkl",
        "Cells_Intestinal_Tract.pkl",
    ],
    "intestine": [
        "Immune_All_Low.pkl",
        "Pan_Fetal_Human.pkl",
        "Cells_Intestinal_Tract.pkl",
    ],
    # Liver
    "liver": [
        "Immune_All_Low.pkl",
        "Pan_Fetal_Human.pkl",
        "Healthy_Human_Liver.pkl",
    ],
    # Lung
    "lung": [
        "Immune_All_Low.pkl",
        "Pan_Fetal_Human.pkl",
        "Human_Lung_Atlas.pkl",
    ],
    "lung_airway": [
        "Immune_All_Low.pkl",
        "Pan_Fetal_Human.pkl",
        "Cells_Lung_Airway.pkl",
    ],
    "lung_cancer": [
        "Immune_All_Low.pkl",
        "Human_Lung_Atlas.pkl",
    ],
    # Heart
    "heart": [
        "Immune_All_Low.pkl",
        "Pan_Fetal_Human.pkl",
        "Healthy_Adult_Heart.pkl",
    ],
    # Breast
    "breast": [
        "Immune_All_Low.pkl",
        "Pan_Fetal_Human.pkl",
        "Cells_Adult_Breast.pkl",
    ],
    # Skin
    "skin": [
        "Immune_All_Low.pkl",
        "Pan_Fetal_Human.pkl",
        "Adult_Human_Skin.pkl",
    ],
    # Pancreas
    "pancreas": [
        "Immune_All_Low.pkl",
        "Pan_Fetal_Human.pkl",
        "Adult_Human_PancreaticIslet.pkl",
    ],
    # Brain
    "brain": [
        "Immune_All_Low.pkl",
        "Pan_Fetal_Human.pkl",
        "Adult_Human_MTG.pkl",
    ],
    # Tonsil
    "tonsil": [
        "Immune_All_Low.pkl",
        "Pan_Fetal_Human.pkl",
        "Cells_Human_Tonsil.pkl",
    ],
    # Blood/Immune
    "blood": [
        "Immune_All_Low.pkl",
        "Pan_Fetal_Human.pkl",
    ],
    "pbmc": [
        "Immune_All_Low.pkl",
        "Pan_Fetal_Human.pkl",
    ],
}


def get_models_for_tissue(tissue: str) -> List[str]:
    """
    Get recommended CellTypist models for a tissue type.

    Parameters
    ----------
    tissue : str
        Tissue name (e.g., "liver", "lung", "colon").

    Returns
    -------
    List[str]
        List of model names/paths.

    Examples
    --------
    >>> from spatialcore.annotation import get_models_for_tissue
    >>> models = get_models_for_tissue("liver")
    >>> print(models)
    ['Immune_All_Low.pkl', 'Pan_Fetal_Human.pkl', 'Healthy_Human_Liver.pkl']
    """
    tissue_lower = tissue.lower().strip()
    return TISSUE_MODEL_PRESETS.get(tissue_lower, TISSUE_MODEL_PRESETS["unknown"])


# ============================================================================
# Model Validation
# ============================================================================

def _validate_gene_overlap(
    model,
    data_genes: set,
    min_overlap_pct: float = 25.0,
) -> Dict[str, Any]:
    """
    Validate gene overlap between model and data.

    Parameters
    ----------
    model
        Loaded CellTypist model.
    data_genes : set
        Gene names from query data.
    min_overlap_pct : float, default 25.0
        Minimum required overlap percentage.

    Returns
    -------
    Dict[str, Any]
        Overlap statistics and pass/fail status.
    """
    model_genes = set(model.features)
    overlap_genes = model_genes & data_genes
    overlap_pct = 100 * len(overlap_genes) / len(model_genes) if model_genes else 0

    return {
        "n_model_genes": len(model_genes),
        "n_data_genes": len(data_genes),
        "n_overlap": len(overlap_genes),
        "overlap_pct": overlap_pct,
        "passes_threshold": overlap_pct >= min_overlap_pct,
    }


def _prepare_for_celltypist(
    adata: ad.AnnData,
    copy: bool = True,
    status: Optional[Dict[str, Any]] = None,
) -> ad.AnnData:
    """
    Prepare AnnData for CellTypist prediction.

    Uses check_normalization_status() to detect data state and
    ensure_normalized() to normalize if needed.

    Parameters
    ----------
    adata : AnnData
        Input data (raw counts or normalized).
    copy : bool, default True
        Return a copy (CellTypist prediction is destructive).
    status : dict, optional
        Normalization status from check_normalization_status() on the full data.

    Returns
    -------
    AnnData
        Copy with log1p(10k) normalized data in X.

    Raises
    ------
    ValueError
        If data cannot be safely normalized.
    """
    from spatialcore.core.utils import check_normalization_status
    from spatialcore.annotation.loading import ensure_normalized

    # Copy if requested (CellTypist modifies X during prediction)
    adata_ct = adata.copy() if copy else adata

    # Check current state (allow caller to pass full-data status)
    if status is None:
        status = check_normalization_status(adata_ct)
    logger.info(f"Input data state: x_state={status['x_state']}, raw_source={status['raw_source']}")

    # Use ensure_normalized to handle all cases
    if status["x_state"] != "log1p_10k":
        if status["raw_source"] is None and not status["is_usable"]:
            raise ValueError(
                f"Cannot prepare data for CellTypist.\n"
                f"  Detected X state: {status['x_state']}\n"
                f"  Estimated target_sum: {status.get('x_target_sum', 'N/A')}\n"
                f"  Raw counts found: None\n\n"
                f"Provide data with:\n"
                f"  1. Raw counts in .X, .layers['counts'], or .raw.X\n"
                f"  2. Or log1p(10k) normalized data in .X"
            )
        adata_ct = ensure_normalized(adata_ct, target_sum=1e4, copy=False)

    logger.info(f"Prepared: {adata_ct.n_obs:,} cells x {adata_ct.n_vars:,} genes")
    return adata_ct


# ============================================================================
# Main Annotation Function
# ============================================================================

def annotate_celltypist(
    adata: ad.AnnData,
    tissue: str = "unknown",
    ensemble_mode: bool = True,
    custom_model_path: Optional[Union[str, Path]] = None,
    majority_voting: bool = False,
    over_clustering: Optional[str] = None,
    min_prop: float = 0.0,
    min_gene_overlap_pct: float = 25.0,
    min_confidence: float = 0.5,
    store_decision_scores: bool = True,
    confidence_transform: Optional[ConfidenceMethod] = "zscore",
    batch_size: Optional[int] = None,
    copy: bool = False,
) -> ad.AnnData:
    """
    Annotate cells using CellTypist with tissue-specific models.

    Algorithm:
    1. Load tissue-specific model preset (or custom model)
    2. Validate gene overlap for each model (skip if <25%)
    3. Validate normalization in specified layer (log1p, ~10k sum)
    4. Run prediction with native celltypist.annotate()
    5. Ensemble: take highest confidence per cell across models

    Parameters
    ----------
    adata : AnnData
        AnnData object to annotate.
    tissue : str, default "unknown"
        Tissue type for model selection (e.g., "liver", "lung", "colon").
    ensemble_mode : bool, default True
        Use multiple tissue-specific models and ensemble results.
    custom_model_path : str or Path, optional
        Path to custom .pkl model (overrides tissue preset).
    majority_voting : bool, default False
        Use CellTypist's native majority voting within clusters.
        **Default False for spatial data** - voting can collapse cell types.
    over_clustering : str, optional
        Column in adata.obs for cluster-based voting (e.g., "leiden").
    min_prop : float, default 0.0
        Minimum proportion for subcluster assignment (0.0 = no threshold).
    min_gene_overlap_pct : float, default 25.0
        Skip models with less than this gene overlap.
    min_confidence : float, default 0.5
        Minimum confidence threshold for cell type assignment.
        Cells below this threshold are labeled "Unassigned".
        Set to 0.0 to disable filtering (assign all cells).
    store_decision_scores : bool, default True
        Store full decision score matrix in adata.obsm for downstream analysis.
        Stores in adata.obsm["cell_type_decision_scores"].
    confidence_transform : {"raw", "zscore", "softmax", "minmax"} or None, default "zscore"
        Transform method for confidence scores. "zscore" is recommended for
        spatial data. Set to None to skip transformation.
    batch_size : int, optional
        If provided, run CellTypist prediction in batches of this many cells to
        reduce peak memory usage. Requires majority_voting=False.
    copy : bool, default False
        If True, return a copy.

    Returns
    -------
    AnnData
        AnnData with new columns in obs (CellxGene standard names):
        - cell_type: Final cell type labels
        - cell_type_confidence: Transformed confidence (z-score by default)
        - cell_type_confidence_raw: Winning-model probability (CellTypist default).
          Decision scores are stored separately in obsm when available.
        - cell_type_model: Which model contributed each prediction
        - cell_type_predicted: Per-cell predictions (before any voting)

        And optionally in obsm (if store_decision_scores=True):
        - cell_type_decision_scores: Full decision score matrix (n_cells x n_types)

    Notes
    -----
    **Input data handling:**

    Input data can be:
    - Raw counts in .X (will be normalized automatically)
    - Raw counts in .layers['counts'] or .raw.X (will be found and normalized)
    - Already log1p(10k) normalized in .X (will be used directly)

    The function auto-detects data state using check_normalization_status()
    and normalizes using ensure_normalized() if needed.

    **Majority voting:**

    For spatial data, majority_voting=False is recommended because:
    1. Spatial clustering may be coarse (few clusters)
    2. Voting assigns dominant type to ALL cells in cluster
    3. This can collapse 13 cell types to 2 types

    **Ensemble mode confidence (LIMITATION):**

    When using ensemble_mode=True with pre-trained models, the z-score confidence
    transformation is NOT applied. This is because each model has different cell
    types in its decision matrix, making it complex to combine them for z-score
    computation. In this case, raw confidence values are used instead.

    This limitation does NOT affect:
    - Single custom model: ``annotate_celltypist(custom_model_path="model.pkl")``
    - ``train_and_annotate()`` pipeline (always uses single custom model)

    Future versions will compute z-score during the ensemble loop using each
    winning model's probability row.

    Examples
    --------
    >>> from spatialcore.annotation import annotate_celltypist
    >>> adata = annotate_celltypist(
    ...     adata,
    ...     tissue="liver",
    ...     ensemble_mode=True,
    ...     majority_voting=False,  # Default for spatial
    ... )
    >>> adata.obs[["cell_type", "cell_type_confidence"]].head()
    """
    try:
        import celltypist
        from celltypist import models
    except ImportError:
        raise ImportError(
            "celltypist is required. Install with: pip install celltypist"
        )

    if batch_size is not None:
        if isinstance(batch_size, bool) or not isinstance(batch_size, int) or batch_size <= 0:
            raise ValueError("batch_size must be a positive integer or None.")
        if majority_voting:
            raise ValueError(
                "batch_size requires majority_voting=False because voting "
                "cannot be computed across batches."
            )

    if copy:
        adata = adata.copy()

    # Determine models to run
    if custom_model_path:
        models_to_run = [str(custom_model_path)]
        logger.info(f"Using custom model: {custom_model_path}")
    elif ensemble_mode:
        models_to_run = get_models_for_tissue(tissue)
        logger.info(f"Using {len(models_to_run)} models for tissue '{tissue}'")
    else:
        models_to_run = ["Immune_All_Low.pkl"]
        logger.info("Using single model: Immune_All_Low.pkl")

    # Load models and validate gene overlap
    loaded_models = {}
    all_overlap_genes = set()
    data_genes = set(adata.var_names)

    for model_name in models_to_run:
        try:
            if Path(model_name).exists():
                loaded_model = models.Model.load(model_name)
            else:
                # Try to load from CellTypist's model collection
                try:
                    loaded_model = models.Model.load(model=model_name)
                except Exception:
                    logger.info(f"Downloading model: {model_name}")
                    models.download_models(model=model_name)
                    loaded_model = models.Model.load(model=model_name)
        except Exception as e:
            logger.warning(f"Failed to load model {model_name}: {e}")
            continue

        # Validate gene overlap
        overlap_info = _validate_gene_overlap(
            loaded_model, data_genes, min_gene_overlap_pct
        )

        if not overlap_info["passes_threshold"]:
            logger.warning(
                f"Skipping {model_name}: only {overlap_info['overlap_pct']:.1f}% gene overlap"
            )
            continue

        logger.info(
            f"  {model_name}: {overlap_info['overlap_pct']:.1f}% overlap "
            f"({overlap_info['n_overlap']}/{overlap_info['n_model_genes']} genes)"
        )

        loaded_models[model_name] = loaded_model
        all_overlap_genes.update(set(loaded_model.features) & data_genes)

    if not loaded_models:
        raise ValueError(
            "No models passed gene overlap threshold. "
            "Consider training a custom model for your panel genes."
        )

    # Determine cluster column for voting (validate early before gene subsetting)
    cluster_col = over_clustering or (
        "leiden" if "leiden" in adata.obs.columns else None
    )

    if majority_voting:
        if cluster_col is None or cluster_col not in adata.obs.columns:
            raise ValueError(
                "majority_voting=True requires a valid cluster column. "
                "Provide over_clustering or add a cluster column (e.g., 'leiden') "
                f"to adata.obs. Available columns: {list(adata.obs.columns)}"
            )

    # Validate and prepare data for CellTypist (avoid full-data copy)
    status = check_normalization_status(adata)

    # Subset to overlapping genes (copy only the needed genes)
    genes_mask = adata.var_names.isin(all_overlap_genes)
    adata_subset = adata[:, genes_mask].copy()

    # If gene set changed, re-normalize after subsetting to match training basis
    if adata_subset.n_vars != adata.n_vars:
        subset_status = check_normalization_status(adata_subset)
        if subset_status["raw_source"] is not None:
            from spatialcore.annotation.loading import _copy_raw_to_x

            logger.info(
                f"Re-normalizing after gene subset ({adata_subset.n_vars}/{adata.n_vars} genes) "
                f"using raw counts ({subset_status['raw_source']})."
            )
            if (
                subset_status["raw_source"] == "raw.X"
                and adata_subset.raw is not None
                and adata_subset.raw.n_vars != adata_subset.n_vars
            ):
                adata_subset.X = adata_subset.raw[:, adata_subset.var_names].X.copy()
            else:
                _copy_raw_to_x(adata_subset, subset_status["raw_source"])
            sc.pp.normalize_total(
                adata_subset, target_sum=1e4, exclude_highly_expressed=False
            )
            sc.pp.log1p(adata_subset)
            status = None
        elif subset_status["x_state"] == "log1p_10k":
            logger.info(
                f"Re-normalizing after gene subset ({adata_subset.n_vars}/{adata.n_vars} genes) "
                "from log1p(10k) data."
            )
            from scipy.sparse import issparse

            if issparse(adata_subset.X):
                adata_subset.X = adata_subset.X.copy()
                adata_subset.X.data = np.expm1(adata_subset.X.data)
            else:
                adata_subset.X = np.expm1(adata_subset.X)
            sc.pp.normalize_total(
                adata_subset, target_sum=1e4, exclude_highly_expressed=False
            )
            sc.pp.log1p(adata_subset)
            status = None
        else:
            raise ValueError(
                "Gene subset requires re-normalization after subsetting, but no raw "
                "counts or log1p(10k) data are available. Provide raw counts in "
                "adata.layers['counts'] or adata.raw.X, or ensure adata.X is "
                "log1p(10k) before annotation."
            )

    adata_subset = _prepare_for_celltypist(adata_subset, copy=False, status=status)
    logger.info(f"Predicting on {adata_subset.n_obs:,} cells × {adata_subset.n_vars:,} genes")

    # Copy cluster info to subset if using voting
    if majority_voting and cluster_col and cluster_col in adata.obs.columns:
        adata_subset.obs[cluster_col] = adata.obs[cluster_col].values

    # Run predictions for each model
    all_model_predictions = {}
    decision_matrix_for_scores = None
    collect_decision_scores = store_decision_scores and len(loaded_models) == 1

    for model_name, loaded_model in loaded_models.items():
        logger.info(f"  Running {model_name}...")

        if batch_size is None or batch_size >= adata_subset.n_obs:
            prediction = celltypist.annotate(
                adata_subset,
                model=loaded_model,
                mode="best match",
                majority_voting=majority_voting,
                over_clustering=cluster_col if majority_voting else None,
                min_prop=min_prop,
            )

            # Get labels based on whether voting was enabled
            if majority_voting and "majority_voting" in prediction.predicted_labels.columns:
                labels = prediction.predicted_labels["majority_voting"]
            else:
                labels = prediction.predicted_labels["predicted_labels"]

            confidence = prediction.probability_matrix.max(axis=1).values
            if collect_decision_scores:
                decision_matrix_for_scores = prediction.decision_matrix
        else:
            n_cells = adata_subset.n_obs
            total_batches = (n_cells + batch_size - 1) // batch_size
            logger.info(f"    Batching enabled (batch_size={batch_size}, n_batches={total_batches})")

            all_labels = []
            all_confidences = []
            all_decision_scores = []

            for start_idx in range(0, n_cells, batch_size):
                end_idx = min(start_idx + batch_size, n_cells)
                adata_batch = adata_subset[start_idx:end_idx].copy()

                prediction = celltypist.annotate(
                    adata_batch,
                    model=loaded_model,
                    mode="best match",
                    majority_voting=False,
                    over_clustering=None,
                    min_prop=min_prop,
                )

                all_labels.append(prediction.predicted_labels["predicted_labels"])
                all_confidences.append(prediction.probability_matrix.max(axis=1).values)
                if collect_decision_scores:
                    all_decision_scores.append(prediction.decision_matrix)

                del adata_batch
                gc.collect()

            labels = pd.concat(all_labels)
            confidence = np.concatenate(all_confidences)
            if collect_decision_scores:
                decision_matrix_for_scores = pd.concat(all_decision_scores)

        all_model_predictions[model_name] = (labels, confidence)

    gc.collect()

    # Combine predictions (ensemble: highest confidence wins)
    if len(loaded_models) == 1:
        model_name = list(all_model_predictions.keys())[0]
        labels, confidence = all_model_predictions[model_name]
        per_cell_predictions = labels
        per_cell_confidence = confidence
        per_cell_source_model = pd.Series([model_name] * len(labels), index=labels.index)
    else:
        # Multi-model ensemble
        cell_indices = list(all_model_predictions.values())[0][0].index
        final_labels = []
        final_confidence = []
        final_source_model = []

        for i, cell_idx in enumerate(cell_indices):
            best_conf = -1.0
            best_label = "Unknown"
            best_model = "none"

            for model_name, (labels, confidence) in all_model_predictions.items():
                cell_conf = confidence[i]
                if cell_conf > best_conf:
                    best_conf = cell_conf
                    best_label = labels.iloc[i]
                    best_model = model_name

            final_labels.append(best_label)
            final_confidence.append(best_conf)
            final_source_model.append(best_model)

        per_cell_predictions = pd.Series(final_labels, index=cell_indices)
        per_cell_confidence = np.array(final_confidence)
        per_cell_source_model = pd.Series(final_source_model, index=cell_indices)

    # Store results (CellxGene standard + required provenance)
    adata.obs["cell_type_predicted"] = per_cell_predictions.values
    adata.obs["cell_type_confidence_raw"] = per_cell_confidence
    adata.obs["cell_type_model"] = per_cell_source_model.values

    # Apply confidence threshold (post-hoc filter)
    # Convert to numpy array of strings to allow "Unassigned" assignment
    final_labels = np.array(per_cell_predictions.values, dtype=object)
    if min_confidence > 0.0:
        low_conf_mask = per_cell_confidence < min_confidence
        n_unassigned = low_conf_mask.sum()
        if n_unassigned > 0:
            final_labels[low_conf_mask] = "Unassigned"
            logger.info(
                f"Confidence filter: {n_unassigned:,} cells ({100*n_unassigned/len(final_labels):.1f}%) "
                f"below {min_confidence} threshold -> 'Unassigned'"
            )
    adata.obs["cell_type"] = pd.Categorical(final_labels)

    # Store decision scores if requested (uses last model's prediction for now)
    # In ensemble mode, only stores the winning model's scores
    if store_decision_scores:
        if len(loaded_models) == 1:
            if decision_matrix_for_scores is not None:
                adata = extract_decision_scores(
                    adata,
                    SimpleNamespace(decision_matrix=decision_matrix_for_scores),
                    key_added="cell_type",
                )
                logger.info(
                    "Stored decision scores in adata.obsm['cell_type_decision_scores']"
                )
            else:
                # For ensemble, we'd need to combine scores - for now, store best model's scores
                # Get the last prediction result for decision matrix access
                model_name = list(loaded_models.keys())[0]  # Use first model for scores
                loaded_model = loaded_models[model_name]

                # Re-run prediction to get decision matrix (this is a limitation for ensemble)
                # Future: Store during prediction loop
                prediction_for_scores = celltypist.annotate(
                    adata_subset,
                    model=loaded_model,
                    mode="best match",
                    majority_voting=False,
                )
                adata = extract_decision_scores(
                    adata,
                    prediction_for_scores,
                    key_added="cell_type",
                )
                logger.info(
                    "Stored decision scores in adata.obsm['cell_type_decision_scores']"
                )
        else:
            logger.warning(
                "store_decision_scores=True requested, but decision scores are not "
                "stored in ensemble mode (multiple models)."
            )

    # Apply confidence transformation if requested
    # Store as main confidence column (cell_type_confidence) per CellxGene standard
    if confidence_transform is not None and "cell_type_decision_scores" in adata.obsm:
        decision_scores = adata.obsm["cell_type_decision_scores"]
        transformed_conf = transform_confidence(decision_scores, method=confidence_transform)
        adata.obs["cell_type_confidence"] = transformed_conf
        logger.info(
            f"Applied {confidence_transform} confidence transform "
            f"(mean={transformed_conf.mean():.3f})"
        )
    else:
        if confidence_transform is not None and "cell_type_decision_scores" not in adata.obsm:
            logger.warning(
                "Confidence transform requested, but decision scores are missing. "
                "Using raw confidence values instead."
            )
        # Use raw confidence if no transform available (ensemble mode limitation)
        # TODO: Compute z-score during ensemble loop using winning model's probability row
        adata.obs["cell_type_confidence"] = per_cell_confidence

    # Log summary
    n_types = adata.obs["cell_type"].nunique()
    mean_conf = np.mean(per_cell_confidence)
    logger.info(f"Annotation complete: {n_types} cell types, mean confidence: {mean_conf:.3f}")

    return adata


def get_annotation_summary(adata: ad.AnnData) -> pd.DataFrame:
    """
    Get summary of CellTypist annotations.

    Parameters
    ----------
    adata : AnnData
        Annotated AnnData object.

    Returns
    -------
    pd.DataFrame
        Summary with columns: cell_type, n_cells, pct_total, mean_confidence.

    Examples
    --------
    >>> from spatialcore.annotation import get_annotation_summary
    >>> summary = get_annotation_summary(adata)
    >>> print(summary.head())
    """
    label_col = "cell_type"
    conf_col = "cell_type_confidence"

    if label_col not in adata.obs.columns or conf_col not in adata.obs.columns:
        raise ValueError("No cell type annotations found. Run annotate_celltypist first.")

    summary = adata.obs.groupby(label_col).agg({
        conf_col: ["count", "mean"],
    })
    summary.columns = ["n_cells", "mean_confidence"]
    summary["pct_total"] = 100 * summary["n_cells"] / adata.n_obs
    summary = summary.sort_values("n_cells", ascending=False).reset_index()
    summary.columns = ["cell_type", "n_cells", "mean_confidence", "pct_total"]

    return summary
