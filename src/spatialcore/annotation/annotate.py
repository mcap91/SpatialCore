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
import gc

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from spatialcore.core.logging import get_logger
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


def _validate_celltypist_input(
    adata: ad.AnnData,
    norm_layer: str = "norm",
) -> ad.AnnData:
    """
    Validate and prepare AnnData for CellTypist.

    Validates that data is properly normalized (log1p, ~10k sum).
    Does NOT modify or normalize data - errors if validation fails.

    Parameters
    ----------
    adata : AnnData
        Input data with normalized layer.
    norm_layer : str, default "norm"
        Layer containing log1p(10k) normalized data.

    Returns
    -------
    AnnData
        Copy with validated layer in X.

    Raises
    ------
    ValueError
        If layer doesn't exist or data is not properly normalized.
    """
    # Step 1: Check layer exists (NO fallback)
    if norm_layer not in adata.layers:
        available = list(adata.layers.keys())
        raise ValueError(
            f"Layer '{norm_layer}' not found in adata.layers.\n"
            f"Available layers: {available}\n"
            f"Ensure normalization has been run before CellTypist annotation."
        )

    logger.info(f"Validating adata.layers['{norm_layer}'] for CellTypist...")
    data = adata.layers[norm_layer]

    # Handle sparse matrices
    if hasattr(data, "toarray"):
        sample = data[:1000].toarray() if data.shape[0] > 1000 else data.toarray()
    else:
        sample = data[:1000] if data.shape[0] > 1000 else data

    # Step 2: Check log-transformed (NO fallback)
    data_max = float(np.max(sample))
    if data_max > 50:
        raise ValueError(
            f"Data in layer '{norm_layer}' does not appear to be log-transformed.\n"
            f"  Max value: {data_max:.2f} (expected < 15 for log1p data)\n"
            f"  Run normalization pipeline before CellTypist annotation."
        )
    logger.info(f"  [OK] Log-transformed (max={data_max:.2f})")

    # Step 3: Check sum ~10000 (NO fallback)
    original_sum = np.expm1(sample).sum(axis=1)
    mean_sum = float(np.mean(original_sum))
    if abs(mean_sum - 10000) / 10000 > 0.1:  # 10% tolerance
        raise ValueError(
            f"Data in layer '{norm_layer}' not normalized to 10000 counts.\n"
            f"  Observed mean sum: {mean_sum:.0f} (expected ~10000)\n"
            f"  Normalize with: sc.pp.normalize_total(adata, target_sum=10000)"
        )
    logger.info(f"  [OK] Sum validation passed (mean={mean_sum:.0f})")

    # Step 4: Create copy with validated layer in X
    adata_ct = adata.copy()
    adata_ct.X = adata.layers[norm_layer].copy()
    logger.info(f"  [OK] Prepared: {adata_ct.n_obs:,} cells x {adata_ct.n_vars:,} genes")

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
    norm_layer: str = "norm",
    store_decision_scores: bool = True,
    confidence_transform: Optional[ConfidenceMethod] = "zscore",
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
    norm_layer : str, default "norm"
        Layer containing log1p(10k) normalized data. Must exist in adata.layers.
    store_decision_scores : bool, default True
        Store full decision score matrix in adata.obsm for downstream analysis.
        Stores in adata.obsm["celltypist_decision_scores"].
    confidence_transform : {"raw", "zscore", "softmax", "minmax"} or None, default "zscore"
        Transform method for confidence scores. "zscore" is recommended for
        spatial data. Set to None to skip transformation.
    copy : bool, default False
        If True, return a copy.

    Returns
    -------
    AnnData
        AnnData with new columns in obs (CellxGene standard names):
        - cell_type: Final cell type labels
        - cell_type_confidence: Transformed confidence (z-score by default)
        - cell_type_confidence_raw: Raw confidence scores from CellTypist
        - cell_type_model: Which model contributed each prediction
        - cell_type_original: Per-cell predictions (before any voting)

        And optionally in obsm (if store_decision_scores=True):
        - cell_type_decision_scores: Full decision score matrix (n_cells x n_types)

    Notes
    -----
    For spatial data, majority_voting=False is recommended because:
    1. Spatial clustering may be coarse (few clusters)
    2. Voting assigns dominant type to ALL cells in cluster
    3. This can collapse 13 cell types to 2 types

    Examples
    --------
    >>> from spatialcore.annotation import annotate_celltypist
    >>> adata = annotate_celltypist(
    ...     adata,
    ...     tissue="liver",
    ...     ensemble_mode=True,
    ...     majority_voting=False,  # Default for spatial
    ... )
    >>> adata.obs[["celltypist", "celltypist_confidence"]].head()
    """
    try:
        import celltypist
        from celltypist import models
    except ImportError:
        raise ImportError(
            "celltypist is required. Install with: pip install celltypist"
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

    # Validate and prepare data for CellTypist
    adata_for_prediction = _validate_celltypist_input(adata, norm_layer=norm_layer)

    # Subset to overlapping genes
    genes_mask = adata_for_prediction.var_names.isin(all_overlap_genes)
    adata_subset = adata_for_prediction[:, genes_mask].copy()
    logger.info(f"Predicting on {adata_subset.n_obs:,} cells × {adata_subset.n_vars:,} genes")

    # Determine cluster column for voting
    cluster_col = over_clustering or (
        "leiden" if "leiden" in adata.obs.columns else None
    )

    # Copy cluster info to subset if using voting
    if majority_voting and cluster_col and cluster_col in adata.obs.columns:
        adata_subset.obs[cluster_col] = adata.obs[cluster_col].values

    # Run predictions for each model
    all_model_predictions = {}

    for model_name, loaded_model in loaded_models.items():
        logger.info(f"  Running {model_name}...")

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

    # Store results (CellxGene standard column names)
    adata.obs["cell_type_original"] = per_cell_predictions.values
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
        # For ensemble, we'd need to combine scores - for now, store best model's scores
        # Get the last prediction result for decision matrix access
        model_name = list(loaded_models.keys())[0]  # Use first model for scores
        loaded_model = loaded_models[model_name]

        # Re-run prediction to get decision matrix (this is a limitation for ensemble)
        # Future: Store during prediction loop
        if len(loaded_models) == 1:
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
                f"Stored decision scores in adata.obsm['cell_type_decision_scores']"
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
        # Use raw confidence if no transform available
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
    # Support both old (celltypist) and new (cell_type) column names
    label_col = "cell_type" if "cell_type" in adata.obs.columns else "celltypist"
    conf_col = "cell_type_confidence" if "cell_type_confidence" in adata.obs.columns else "celltypist_confidence"

    if label_col not in adata.obs.columns:
        raise ValueError("No cell type annotations found. Run annotate_celltypist first.")

    summary = adata.obs.groupby(label_col).agg({
        conf_col: ["count", "mean"],
    })
    summary.columns = ["n_cells", "mean_confidence"]
    summary["pct_total"] = 100 * summary["n_cells"] / adata.n_obs
    summary = summary.sort_values("n_cells", ascending=False).reset_index()
    summary.columns = ["cell_type", "n_cells", "mean_confidence", "pct_total"]

    return summary
