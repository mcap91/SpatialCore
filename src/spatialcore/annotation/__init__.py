"""
Cell type annotation: Unified hub for CellTypist, ontology mapping, and reference data.

This module consolidates all cell typing functionality:
- CellTypist model training and annotation
- Confidence scoring and filtering
- Cell Ontology (CL) ID mapping
- Reference data loading from CellxGene, Synapse, GCS
- Marker gene management
- Validation utilities
"""

# Training & model management
from spatialcore.annotation.training import (
    combine_references,
    train_celltypist_model,
    get_panel_genes,
    get_model_gene_overlap,
    get_training_summary,
    HIGH_CONTRAST_PALETTE,
    generate_color_scheme,
    save_model_artifacts,
    subsample_balanced,
    DEFAULT_EXCLUDE_LABELS,
)

# High-level pipeline
from spatialcore.annotation.pipeline import (
    train_and_annotate,
    train_and_annotate_config,
    TrainingConfig,
)

# Annotation (CellTypist wrapper)
from spatialcore.annotation.annotate import (
    annotate_celltypist,
    get_models_for_tissue,
    get_annotation_summary,
    TISSUE_MODEL_PRESETS,
)

# Validation
from spatialcore.annotation.validation import (
    validate_cell_type_column,
    validate_multiple_columns,
    CellTypeValidationResult,
    ValidationIssue,
)

# Confidence scoring & filtering
from spatialcore.annotation.confidence import (
    transform_confidence,
    extract_decision_scores,
    filter_low_confidence,
    filter_low_count_types,
    compute_confidence_from_obsm,
    filter_by_marker_validation,
    ConfidenceMethod,
)

# Markers
from spatialcore.annotation.markers import (
    load_canonical_markers,
    match_to_canonical,
    get_markers_for_type,
    list_available_cell_types,
)

# Ontology mapping (consolidated from ontology/)
from spatialcore.annotation.ontology import (
    add_ontology_ids,
    has_ontology_ids,
    search_ontology_index,
    load_ontology_index,
    create_mapping_table,
    OntologyMappingResult,
    UNKNOWN_CELL_TYPE_ID,
    UNKNOWN_CELL_TYPE_NAME,
)

from spatialcore.annotation.patterns import (
    CELL_TYPE_PATTERNS,
    get_canonical_term,
)

from spatialcore.annotation.expression import (
    evaluate_ontology_expression,
    get_ontology_ids_in_expression,
)

# Reference data loading (consolidated from reference/)
from spatialcore.annotation.cellxgene import (
    download_cellxgene_reference,
    list_available_datasets,
    load_ensembl_to_hugo_mapping,
    normalize_gene_names,
    check_normalization_status,
    query_cellxgene_census,
)

from spatialcore.annotation.loading import (
    load_adata_backed,
    subsample_adata,
    ensure_normalized,
    get_available_memory_gb,
    estimate_adata_memory_gb,
    get_loading_summary,
)

from spatialcore.annotation.discovery import (
    discover_training_data,
    DiscoveredDataset,
    print_discovery_summary,
    load_local_metadata,
    query_local_references,
    create_metadata_template,
)

# Data acquisition (Phase 1: download from sources → store to local/cloud)
from spatialcore.annotation.acquisition import (
    acquire_reference,
    resolve_uri_to_local,
)

from spatialcore.annotation.synapse import (
    download_synapse_reference,
    authenticate_synapse,
    get_synapse_entity_info,
    list_synapse_folder,
    download_synapse_folder,
)

__all__ = [
    # =====================================================================
    # Training & Model Management
    # =====================================================================
    "combine_references",
    "train_celltypist_model",
    "get_panel_genes",
    "get_model_gene_overlap",
    "get_training_summary",
    "HIGH_CONTRAST_PALETTE",
    "generate_color_scheme",
    "save_model_artifacts",
    "subsample_balanced",
    "DEFAULT_EXCLUDE_LABELS",

    # =====================================================================
    # High-Level Pipeline
    # =====================================================================
    "train_and_annotate",
    "train_and_annotate_config",
    "TrainingConfig",

    # =====================================================================
    # Annotation (CellTypist)
    # =====================================================================
    "annotate_celltypist",
    "get_models_for_tissue",
    "get_annotation_summary",
    "TISSUE_MODEL_PRESETS",

    # =====================================================================
    # Validation
    # =====================================================================
    "validate_cell_type_column",
    "validate_multiple_columns",
    "CellTypeValidationResult",
    "ValidationIssue",

    # =====================================================================
    # Confidence Scoring & Filtering
    # =====================================================================
    "transform_confidence",
    "extract_decision_scores",
    "filter_low_confidence",
    "filter_low_count_types",
    "compute_confidence_from_obsm",
    "filter_by_marker_validation",
    "ConfidenceMethod",

    # =====================================================================
    # Markers
    # =====================================================================
    "load_canonical_markers",
    "match_to_canonical",
    "get_markers_for_type",
    "list_available_cell_types",

    # =====================================================================
    # Ontology Mapping
    # =====================================================================
    "add_ontology_ids",
    "has_ontology_ids",
    "search_ontology_index",
    "load_ontology_index",
    "create_mapping_table",
    "OntologyMappingResult",
    "UNKNOWN_CELL_TYPE_ID",
    "UNKNOWN_CELL_TYPE_NAME",
    "CELL_TYPE_PATTERNS",
    "get_canonical_term",
    "evaluate_ontology_expression",
    "get_ontology_ids_in_expression",

    # =====================================================================
    # Reference Data: CellxGene
    # =====================================================================
    "download_cellxgene_reference",
    "list_available_datasets",
    "load_ensembl_to_hugo_mapping",
    "normalize_gene_names",
    "check_normalization_status",
    "query_cellxgene_census",

    # =====================================================================
    # Reference Data: Loading & Memory
    # =====================================================================
    "load_adata_backed",
    "subsample_adata",
    "ensure_normalized",
    "get_available_memory_gb",
    "estimate_adata_memory_gb",
    "get_loading_summary",

    # =====================================================================
    # Reference Data: Discovery
    # =====================================================================
    "discover_training_data",
    "DiscoveredDataset",
    "print_discovery_summary",
    "load_local_metadata",
    "query_local_references",
    "create_metadata_template",

    # =====================================================================
    # Data Acquisition (Phase 1)
    # =====================================================================
    # Use acquire_reference() to download from CellxGene/Synapse and store
    # to local filesystem or cloud storage (GCS, S3). Then use the stored
    # paths with train_and_annotate() or combine_references().
    "acquire_reference",
    "resolve_uri_to_local",

    # =====================================================================
    # Reference Data: Synapse
    # =====================================================================
    "download_synapse_reference",
    "authenticate_synapse",
    "get_synapse_entity_info",
    "list_synapse_folder",
    "download_synapse_folder",
]
