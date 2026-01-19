"""
Training data discovery for local and cloud storage.

This module provides unified data discovery across:
1. Local filesystem paths
2. Google Cloud Storage (gs://) paths

The discovery function auto-detects path type and returns
consistent metadata regardless of storage backend.

Examples
--------
>>> from spatialcore.annotation import discover_training_data
>>> # Local discovery
>>> datasets = discover_training_data("/data/references/liver/")
>>> for ds in datasets:
...     print(f"{ds.name}: {ds.size_human}")

>>> # GCS discovery
>>> datasets = discover_training_data("gs://my-bucket/cellxgene/")
>>> len(datasets)
5
"""

from pathlib import Path
from typing import List, Optional, Union, Literal
from dataclasses import dataclass
from datetime import datetime
import re
import time

from spatialcore.core.logging import get_logger

logger = get_logger(__name__)


@dataclass
class DiscoveredDataset:
    """Metadata for a discovered training dataset.

    Attributes
    ----------
    path : str
        Full path to the dataset (local path or gs:// URL).
    name : str
        Dataset name (filename stem without extension).
    size_bytes : int, optional
        File size in bytes (None if unavailable).
    size_human : str
        Human-readable size (e.g., "2.3 GB").
    storage_type : str
        Storage backend: "local" or "gcs".
    last_modified : str, optional
        ISO timestamp of last modification (None if unavailable).
    """

    path: str
    name: str
    size_bytes: Optional[int]
    size_human: str
    storage_type: Literal["local", "gcs"]
    last_modified: Optional[str]

    def __repr__(self) -> str:
        return (
            f"DiscoveredDataset(name='{self.name}', "
            f"size='{self.size_human}', type='{self.storage_type}')"
        )


def discover_training_data(
    path: Union[str, Path],
    pattern: str = "*.h5ad",
    recursive: bool = False,
) -> List[DiscoveredDataset]:
    """
    Discover available training data files at a path.

    Supports both local filesystem paths and Google Cloud Storage (gs://) paths.
    Auto-detects path type based on prefix.

    Parameters
    ----------
    path : str or Path
        Directory to search. Supports:
        - Local paths: "/data/references/", "C:/Data/references/"
        - GCS paths: "gs://bucket-name/references/"
    pattern : str, default "*.h5ad"
        Glob pattern for matching files (e.g., "*.h5ad", "*.parquet").
    recursive : bool, default False
        If True, search subdirectories recursively.

    Returns
    -------
    List[DiscoveredDataset]
        List of discovered datasets with metadata, sorted by name.

    Raises
    ------
    FileNotFoundError
        If local path does not exist.
    ValueError
        If GCS path is malformed.
    ImportError
        If google-cloud-storage is required but not installed.
    PermissionError
        If access to GCS path is denied.

    Examples
    --------
    >>> from spatialcore.annotation import discover_training_data
    >>> # Local discovery
    >>> datasets = discover_training_data("./references/")
    >>> for ds in datasets:
    ...     print(f"{ds.name}: {ds.size_human}")
    healthy_liver: 2.3 GB
    hcc_liver: 1.8 GB

    >>> # GCS discovery (requires google-cloud-storage)
    >>> datasets = discover_training_data("gs://my-bucket/cellxgene/")
    >>> len(datasets)
    5

    >>> # Recursive search
    >>> datasets = discover_training_data("./references/", recursive=True)

    See Also
    --------
    download_cellxgene_reference : Download pre-configured CellxGene datasets.
    list_available_datasets : List pre-configured CellxGene dataset IDs.
    """
    path_str = str(path)

    if path_str.startswith("gs://"):
        return _discover_gcs_with_retry(path_str, pattern, recursive)
    else:
        return _discover_local(Path(path_str), pattern, recursive)


def _discover_local(
    path: Path,
    pattern: str,
    recursive: bool,
) -> List[DiscoveredDataset]:
    """Discover datasets on local filesystem."""
    if not path.exists():
        raise FileNotFoundError(f"Path does not exist: {path}")

    if not path.is_dir():
        raise ValueError(f"Path is not a directory: {path}")

    glob_method = path.rglob if recursive else path.glob
    files = sorted(glob_method(pattern))

    datasets = []
    for f in files:
        try:
            stat = f.stat()
            size_bytes = stat.st_size
            last_modified = _format_timestamp(stat.st_mtime)
        except OSError:
            size_bytes = None
            last_modified = None

        datasets.append(
            DiscoveredDataset(
                path=str(f.resolve()),
                name=f.stem,
                size_bytes=size_bytes,
                size_human=_format_size(size_bytes) if size_bytes else "unknown",
                storage_type="local",
                last_modified=last_modified,
            )
        )

    logger.info(f"Discovered {len(datasets)} datasets at {path}")
    return datasets


def _discover_gcs(
    path: str,
    pattern: str,
    recursive: bool,
) -> List[DiscoveredDataset]:
    """Discover datasets on Google Cloud Storage."""
    try:
        from google.cloud import storage
    except ImportError:
        raise ImportError(
            "google-cloud-storage is required for GCS paths. "
            "Install with: pip install google-cloud-storage"
        )

    # Parse gs://bucket/prefix
    match = re.match(r"gs://([^/]+)(?:/(.*))?", path)
    if not match:
        raise ValueError(f"Invalid GCS path format: {path}. Expected gs://bucket/prefix")

    bucket_name = match.group(1)
    prefix = match.group(2) or ""
    if prefix and not prefix.endswith("/"):
        prefix += "/"

    # Convert glob pattern to regex for matching
    pattern_regex = _glob_to_regex(pattern)

    client = storage.Client()
    bucket = client.bucket(bucket_name)

    # List blobs - use delimiter only if not recursive
    delimiter = None if recursive else "/"
    blobs = bucket.list_blobs(prefix=prefix, delimiter=delimiter)

    datasets = []
    for blob in blobs:
        # Skip "directory" markers
        if blob.name.endswith("/"):
            continue

        # Check if matches pattern
        filename = blob.name.split("/")[-1]
        if not re.match(pattern_regex, filename):
            continue

        datasets.append(
            DiscoveredDataset(
                path=f"gs://{bucket_name}/{blob.name}",
                name=Path(filename).stem,
                size_bytes=blob.size,
                size_human=_format_size(blob.size) if blob.size else "unknown",
                storage_type="gcs",
                last_modified=blob.updated.isoformat() if blob.updated else None,
            )
        )

    datasets.sort(key=lambda d: d.name)
    logger.info(f"Discovered {len(datasets)} datasets at {path}")
    return datasets


def _discover_gcs_with_retry(
    path: str,
    pattern: str,
    recursive: bool,
    max_retries: int = 3,
    backoff_seconds: float = 1.0,
) -> List[DiscoveredDataset]:
    """GCS discovery with exponential backoff retry."""
    try:
        from google.api_core import exceptions as gcs_exceptions
    except ImportError:
        # If we can't import google.api_core, just try once without retry
        return _discover_gcs(path, pattern, recursive)

    last_error = None
    for attempt in range(max_retries):
        try:
            return _discover_gcs(path, pattern, recursive)
        except gcs_exceptions.ServiceUnavailable as e:
            last_error = e
            if attempt < max_retries - 1:
                wait_time = backoff_seconds * (2**attempt)
                logger.warning(f"GCS unavailable, retrying in {wait_time}s...")
                time.sleep(wait_time)
        except gcs_exceptions.Forbidden as e:
            raise PermissionError(f"Access denied to GCS path: {path}. {e}")
        except gcs_exceptions.NotFound as e:
            raise FileNotFoundError(f"GCS path not found: {path}. {e}")

    raise ConnectionError(
        f"GCS unavailable after {max_retries} retries: {last_error}"
    )


def _format_size(size_bytes: Optional[int]) -> str:
    """Format bytes as human-readable string."""
    if size_bytes is None:
        return "unknown"

    size = float(size_bytes)
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size < 1024:
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} PB"


def _format_timestamp(timestamp: float) -> str:
    """Format Unix timestamp as ISO string."""
    return datetime.fromtimestamp(timestamp).isoformat()


def _glob_to_regex(pattern: str) -> str:
    """Convert glob pattern to regex.

    Examples
    --------
    >>> _glob_to_regex("*.h5ad")
    '^.*\\.h5ad$'
    >>> _glob_to_regex("data_?.csv")
    '^data_.\\.csv$'
    """
    # Escape special regex chars except * and ?
    pattern = re.escape(pattern)
    # Convert glob wildcards to regex
    pattern = pattern.replace(r"\*", ".*")
    pattern = pattern.replace(r"\?", ".")
    return f"^{pattern}$"


def print_discovery_summary(datasets: List[DiscoveredDataset]) -> None:
    """Print a formatted summary of discovered datasets.

    Parameters
    ----------
    datasets : List[DiscoveredDataset]
        List of discovered datasets from discover_training_data.

    Examples
    --------
    >>> datasets = discover_training_data("./references/")
    >>> print_discovery_summary(datasets)
    Found 3 datasets:
      1. healthy_liver     2.3 GB   (local)
      2. hcc_liver         1.8 GB   (local)
      3. colon_atlas       5.1 GB   (local)
    Total: 9.2 GB
    """
    if not datasets:
        print("No datasets found.")
        return

    print(f"Found {len(datasets)} datasets:")
    total_bytes = 0
    for i, ds in enumerate(datasets, 1):
        storage_tag = f"({ds.storage_type})"
        print(f"  {i}. {ds.name:<20} {ds.size_human:>10}   {storage_tag}")
        if ds.size_bytes:
            total_bytes += ds.size_bytes

    if total_bytes > 0:
        print(f"Total: {_format_size(total_bytes)}")


# ============================================================================
# Local Metadata CSV Support
# ============================================================================

def load_local_metadata(
    metadata_csv: Union[str, Path],
    sample_csv: Optional[Union[str, Path]] = None,
) -> "tuple[pd.DataFrame, Optional[pd.DataFrame]]":
    """
    Load local scRNAseq metadata and sample summaries.

    Metadata CSV should contain columns like:
    - file_path: Path to h5ad file
    - tissue: Tissue type
    - condition: Disease/healthy status
    - n_cells: Number of cells
    - label_column: Cell type column name

    Parameters
    ----------
    metadata_csv : str or Path
        Path to metadata CSV file.
    sample_csv : str or Path, optional
        Path to sample-level summary CSV.

    Returns
    -------
    Tuple[pd.DataFrame, Optional[pd.DataFrame]]
        (metadata DataFrame, sample DataFrame or None)

    Examples
    --------
    >>> from spatialcore.annotation.discovery import load_local_metadata
    >>> meta, samples = load_local_metadata("references_metadata.csv")
    >>> print(meta.columns.tolist())
    ['file_path', 'tissue', 'condition', 'n_cells', 'label_column']
    """
    import pandas as pd

    metadata_csv = Path(metadata_csv)
    if not metadata_csv.exists():
        raise FileNotFoundError(f"Metadata CSV not found: {metadata_csv}")

    metadata_df = pd.read_csv(metadata_csv)
    logger.info(f"Loaded metadata: {len(metadata_df)} entries from {metadata_csv}")

    sample_df = None
    if sample_csv is not None:
        sample_csv = Path(sample_csv)
        if sample_csv.exists():
            sample_df = pd.read_csv(sample_csv)
            logger.info(f"Loaded sample summary: {len(sample_df)} entries")
        else:
            logger.warning(f"Sample CSV not found: {sample_csv}")

    return metadata_df, sample_df


def query_local_references(
    metadata_df: "pd.DataFrame",
    tissue: Optional[str] = None,
    condition: Optional[str] = None,
    min_cells: int = 1000,
    file_column: str = "file_path",
    tissue_column: str = "tissue",
    condition_column: str = "condition",
    cells_column: str = "n_cells",
) -> "pd.DataFrame":
    """
    Query local references by tissue/condition filters.

    Parameters
    ----------
    metadata_df : pd.DataFrame
        Metadata DataFrame from load_local_metadata.
    tissue : str, optional
        Filter by tissue type (case-insensitive substring match).
    condition : str, optional
        Filter by condition (e.g., "healthy", "cancer").
    min_cells : int, default 1000
        Minimum cells required.
    file_column : str, default "file_path"
        Column containing file paths.
    tissue_column : str, default "tissue"
        Column containing tissue types.
    condition_column : str, default "condition"
        Column containing conditions.
    cells_column : str, default "n_cells"
        Column containing cell counts.

    Returns
    -------
    pd.DataFrame
        Filtered metadata DataFrame.

    Examples
    --------
    >>> from spatialcore.annotation.discovery import load_local_metadata, query_local_references
    >>> meta, _ = load_local_metadata("references_metadata.csv")
    >>> liver_refs = query_local_references(meta, tissue="liver", min_cells=5000)
    >>> print(liver_refs[[file_column, tissue_column, cells_column]])
    """
    result = metadata_df.copy()

    # Filter by tissue
    if tissue is not None and tissue_column in result.columns:
        tissue_lower = tissue.lower()
        mask = result[tissue_column].astype(str).str.lower().str.contains(tissue_lower, na=False)
        result = result[mask]
        logger.info(f"  Filtered by tissue '{tissue}': {len(result)} remaining")

    # Filter by condition
    if condition is not None and condition_column in result.columns:
        condition_lower = condition.lower()
        mask = result[condition_column].astype(str).str.lower().str.contains(condition_lower, na=False)
        result = result[mask]
        logger.info(f"  Filtered by condition '{condition}': {len(result)} remaining")

    # Filter by min cells
    if cells_column in result.columns:
        result = result[result[cells_column] >= min_cells]
        logger.info(f"  Filtered by min_cells={min_cells}: {len(result)} remaining")

    return result


def create_metadata_template(
    output_path: Union[str, Path],
    discovered_datasets: Optional[List[DiscoveredDataset]] = None,
) -> Path:
    """
    Create a metadata CSV template for local references.

    Optionally pre-populate with discovered datasets.

    Parameters
    ----------
    output_path : str or Path
        Path to save template CSV.
    discovered_datasets : List[DiscoveredDataset], optional
        Pre-populate with discovered datasets.

    Returns
    -------
    Path
        Path to created template.

    Examples
    --------
    >>> from spatialcore.annotation.discovery import discover_training_data, create_metadata_template
    >>> datasets = discover_training_data("./references/")
    >>> create_metadata_template("metadata_template.csv", datasets)
    """
    import pandas as pd

    output_path = Path(output_path)

    if discovered_datasets:
        data = [
            {
                "file_path": ds.path,
                "name": ds.name,
                "tissue": "",
                "condition": "",
                "n_cells": None,
                "label_column": "cell_type",
                "notes": "",
            }
            for ds in discovered_datasets
        ]
        df = pd.DataFrame(data)
    else:
        df = pd.DataFrame(columns=[
            "file_path",
            "name",
            "tissue",
            "condition",
            "n_cells",
            "label_column",
            "notes",
        ])

    df.to_csv(output_path, index=False)
    logger.info(f"Created metadata template: {output_path}")
    return output_path
