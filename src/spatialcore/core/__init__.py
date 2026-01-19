"""Core utilities for SpatialCore: logging, metadata tracking, caching, gene mapping."""

from spatialcore.core.logging import get_logger, setup_logging
from spatialcore.core.metadata import MetadataTracker, update_metadata
from spatialcore.core.cache import cache_result, clear_cache, get_cache_path
from spatialcore.core.utils import (
    # Gene ID mapping (Ensembl → HUGO)
    load_ensembl_to_hugo_mapping,
    normalize_gene_names,
    download_ensembl_mapping,
    is_ensembl_id,
    # Expression normalization detection
    check_normalization_status,
)

__all__ = [
    # Logging
    "get_logger",
    "setup_logging",
    # Metadata
    "MetadataTracker",
    "update_metadata",
    # Caching
    "cache_result",
    "clear_cache",
    "get_cache_path",
    # Gene ID mapping
    "load_ensembl_to_hugo_mapping",
    "normalize_gene_names",
    "download_ensembl_mapping",
    "is_ensembl_id",
    # Expression normalization
    "check_normalization_status",
]
