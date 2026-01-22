"""
Data acquisition utilities for reference datasets.

This module provides a unified interface for downloading reference data
from various sources (CellxGene, Synapse) and storing to local or cloud
destinations (GCS, S3).

This is Phase 1 of the training workflow - run once upstream to acquire
and store reference data before training.

Workflow
--------
1. Use acquire_reference() to download from source → store to destination
2. Use train_and_annotate() with the stored paths for training

Example
-------
>>> from spatialcore.annotation import acquire_reference
>>> # Download from CellxGene → store to GCS
>>> path = acquire_reference(
...     source="cellxgene://human_lung_cell_atlas",
...     output="gs://my-bucket/references/hlca.h5ad",
... )
>>> # Later, use in training pipeline
>>> adata = train_and_annotate(
...     spatial_adata,
...     references=["gs://my-bucket/references/hlca.h5ad"],
...     tissue="lung",
... )
"""

import os
import tempfile
from pathlib import Path
from typing import Dict, Optional, Union, Any

import anndata as ad

from spatialcore.core.logging import get_logger

logger = get_logger(__name__)


# ============================================================================
# Cloud I/O Utilities
# ============================================================================

def _upload_to_gcs(local_path: Path, gcs_uri: str) -> str:
    """
    Upload a local file to Google Cloud Storage.

    Parameters
    ----------
    local_path : Path
        Local file to upload.
    gcs_uri : str
        Destination URI (gs://bucket/path/file.h5ad).

    Returns
    -------
    str
        The GCS URI.
    """
    try:
        from google.cloud import storage
    except ImportError:
        raise ImportError(
            "google-cloud-storage is required for GCS uploads. "
            "Install with: pip install google-cloud-storage"
        )

    # Parse GCS URI: gs://bucket/path/to/file.h5ad
    if not gcs_uri.startswith("gs://"):
        raise ValueError(f"Invalid GCS URI: {gcs_uri}. Must start with gs://")

    uri_path = gcs_uri[5:]  # Remove "gs://"
    parts = uri_path.split("/", 1)
    bucket_name = parts[0]
    blob_path = parts[1] if len(parts) > 1 else Path(local_path).name

    logger.info(f"Uploading to GCS: {gcs_uri}")
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_path)
    blob.upload_from_filename(str(local_path))
    logger.info(f"Upload complete: {gcs_uri}")

    return gcs_uri


def _upload_to_s3(local_path: Path, s3_uri: str) -> str:
    """
    Upload a local file to Amazon S3.

    Parameters
    ----------
    local_path : Path
        Local file to upload.
    s3_uri : str
        Destination URI (s3://bucket/path/file.h5ad).

    Returns
    -------
    str
        The S3 URI.
    """
    try:
        import boto3
    except ImportError:
        raise ImportError(
            "boto3 is required for S3 uploads. "
            "Install with: pip install boto3"
        )

    # Parse S3 URI: s3://bucket/path/to/file.h5ad
    if not s3_uri.startswith("s3://"):
        raise ValueError(f"Invalid S3 URI: {s3_uri}. Must start with s3://")

    uri_path = s3_uri[5:]  # Remove "s3://"
    parts = uri_path.split("/", 1)
    bucket_name = parts[0]
    object_key = parts[1] if len(parts) > 1 else Path(local_path).name

    logger.info(f"Uploading to S3: {s3_uri}")
    s3_client = boto3.client("s3")
    s3_client.upload_file(str(local_path), bucket_name, object_key)
    logger.info(f"Upload complete: {s3_uri}")

    return s3_uri


def _write_output(adata: ad.AnnData, output: Union[str, Path]) -> str:
    """
    Write AnnData to local path or cloud storage.

    Parameters
    ----------
    adata : AnnData
        Data to write.
    output : str or Path
        Output location. Supports:
        - Local path: /data/refs/lung.h5ad
        - GCS: gs://bucket/refs/lung.h5ad
        - S3: s3://bucket/refs/lung.h5ad

    Returns
    -------
    str
        The output path/URI.
    """
    output_str = str(output)

    if output_str.startswith("gs://"):
        # Write to temp file, upload to GCS
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
            tmp_path = Path(tmp.name)
        try:
            logger.info(f"Writing to temporary file: {tmp_path}")
            adata.write_h5ad(tmp_path)
            return _upload_to_gcs(tmp_path, output_str)
        finally:
            if tmp_path.exists():
                tmp_path.unlink()

    elif output_str.startswith("s3://"):
        # Write to temp file, upload to S3
        with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
            tmp_path = Path(tmp.name)
        try:
            logger.info(f"Writing to temporary file: {tmp_path}")
            adata.write_h5ad(tmp_path)
            return _upload_to_s3(tmp_path, output_str)
        finally:
            if tmp_path.exists():
                tmp_path.unlink()

    else:
        # Local path
        local_path = Path(output_str)
        local_path.parent.mkdir(parents=True, exist_ok=True)
        logger.info(f"Writing to local file: {local_path}")
        adata.write_h5ad(local_path)
        return str(local_path)


def _download_from_gcs(gcs_uri: str, local_path: Path) -> Path:
    """Download file from GCS to local path."""
    try:
        from google.cloud import storage
    except ImportError:
        raise ImportError(
            "google-cloud-storage is required for GCS. "
            "Install with: pip install google-cloud-storage"
        )

    if not gcs_uri.startswith("gs://"):
        raise ValueError(f"Invalid GCS URI: {gcs_uri}")

    uri_path = gcs_uri[5:]
    parts = uri_path.split("/", 1)
    bucket_name = parts[0]
    blob_path = parts[1] if len(parts) > 1 else ""

    local_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Downloading from GCS: {gcs_uri}")
    client = storage.Client()
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_path)
    blob.download_to_filename(str(local_path))
    logger.info(f"Downloaded to: {local_path}")

    return local_path


def _download_from_s3(s3_uri: str, local_path: Path) -> Path:
    """Download file from S3 to local path."""
    try:
        import boto3
        from botocore.exceptions import ClientError
    except ImportError:
        raise ImportError(
            "boto3 is required for S3. "
            "Install with: pip install boto3"
        )

    if not s3_uri.startswith("s3://"):
        raise ValueError(f"Invalid S3 URI: {s3_uri}")

    uri_path = s3_uri[5:]
    parts = uri_path.split("/", 1)
    bucket_name = parts[0]
    object_key = parts[1] if len(parts) > 1 else ""

    local_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Downloading from S3: {s3_uri}")
    try:
        s3_client = boto3.client("s3")
        s3_client.download_file(bucket_name, object_key, str(local_path))
        logger.info(f"Downloaded to: {local_path}")
    except ClientError as e:
        error_code = e.response.get("Error", {}).get("Code", "Unknown")
        if error_code == "404":
            raise FileNotFoundError(f"S3 object not found: {s3_uri}")
        elif error_code == "403":
            raise PermissionError(f"Access denied to S3 object: {s3_uri}")
        raise

    return local_path


def resolve_uri_to_local(
    uri: str,
    cache_dir: Path,
    force: bool = False,
) -> Path:
    """
    Resolve a URI to a local file path, downloading if necessary.

    Parameters
    ----------
    uri : str
        Source URI. Supports:
        - Local path: /data/refs/lung.h5ad
        - GCS: gs://bucket/refs/lung.h5ad
        - S3: s3://bucket/refs/lung.h5ad
    cache_dir : Path
        Directory for downloaded files.
    force : bool, default False
        Re-download even if cached.

    Returns
    -------
    Path
        Local file path.
    """
    if uri.startswith("gs://"):
        # Extract filename from URI
        filename = Path(uri[5:].split("/", 1)[1]).name if "/" in uri[5:] else "data.h5ad"
        local_path = cache_dir / filename
        if not local_path.exists() or force:
            return _download_from_gcs(uri, local_path)
        else:
            logger.info(f"Using cached file: {local_path}")
            return local_path

    elif uri.startswith("s3://"):
        filename = Path(uri[5:].split("/", 1)[1]).name if "/" in uri[5:] else "data.h5ad"
        local_path = cache_dir / filename
        if not local_path.exists() or force:
            return _download_from_s3(uri, local_path)
        else:
            logger.info(f"Using cached file: {local_path}")
            return local_path

    else:
        # Local path
        local_path = Path(uri)
        if not local_path.exists():
            raise FileNotFoundError(f"File not found: {local_path}")
        return local_path


# ============================================================================
# Main Acquisition Function
# ============================================================================

def acquire_reference(
    source: str,
    output: Union[str, Path],
    force: bool = False,
    **kwargs,
) -> str:
    """
    Download reference data from a source and store to a destination.

    This is the unified data acquisition function for Phase 1 of the
    training workflow. It handles downloading from public databases
    and storing to local filesystem or cloud storage.

    Parameters
    ----------
    source : str
        Source to download from. Supported schemes:

        - ``cellxgene://dataset_key`` - CellxGene Census dataset
          (e.g., "cellxgene://human_lung_cell_atlas")
        - ``cellxgene://?tissue=lung&disease=normal`` - CellxGene query
        - ``synapse://syn12345678`` - Synapse entity ID

    output : str or Path
        Destination to store the data. Supports:

        - Local path: ``/data/refs/lung.h5ad``
        - GCS: ``gs://bucket/refs/lung.h5ad``
        - S3: ``s3://bucket/refs/lung.h5ad``

    force : bool, default False
        Re-download and overwrite even if output exists.

    **kwargs
        Source-specific options:

        - ``max_cells`` (int): Maximum cells to download (for CellxGene query)
        - ``resolve_hierarchy`` (str): "remove_parents" to drop parent labels
        - ``auth_token`` (str): Synapse authentication token
        - ``tissue``, ``disease``, ``cell_type`` (str): CellxGene query filters

    Returns
    -------
    str
        The output path/URI (same as input output parameter).
        Use this path in train_and_annotate() references list.

    Raises
    ------
    ValueError
        If source scheme is not recognized.
    ImportError
        If required cloud SDK is not installed.

    Examples
    --------
    >>> from spatialcore.annotation import acquire_reference
    >>> # Download from CellxGene → store locally
    >>> path = acquire_reference(
    ...     source="cellxgene://human_lung_cell_atlas",
    ...     output="/data/references/hlca.h5ad",
    ... )
    >>> print(path)
    /data/references/hlca.h5ad

    >>> # Download from CellxGene → store to GCS
    >>> gcs_path = acquire_reference(
    ...     source="cellxgene://human_lung_cell_atlas",
    ...     output="gs://my-bucket/references/hlca.h5ad",
    ... )

    >>> # CellxGene query with filters → store to S3
    >>> s3_path = acquire_reference(
    ...     source="cellxgene://?tissue=liver&disease=normal",
    ...     output="s3://my-bucket/references/healthy_liver.h5ad",
    ...     max_cells=100000,
    ... )

    >>> # Download from Synapse → store locally
    >>> path = acquire_reference(
    ...     source="synapse://syn12345678",
    ...     output="/data/references/lung_ref.h5ad",
    ...     auth_token=os.environ.get("SYNAPSE_AUTH_TOKEN"),
    ... )

    See Also
    --------
    train_and_annotate : Use acquired references for training.
    download_cellxgene_reference : Direct CellxGene download (low-level).
    query_cellxgene_census : CellxGene query with filters (low-level).
    download_synapse_reference : Direct Synapse download (low-level).

    Notes
    -----
    **Recommended Workflow:**

    1. Run ``acquire_reference()`` once to download and store data
    2. Use the returned path in ``train_and_annotate()`` references list
    3. Cloud storage (GCS/S3) enables reproducible training across team members

    **Cloud Storage Benefits:**

    - Versioned datasets for reproducibility
    - Team collaboration with shared references
    - No local disk requirements for large datasets
    """
    output_str = str(output)

    # Check if output exists (for cloud, we'd need to check remotely)
    if not force and not output_str.startswith(("gs://", "s3://")):
        local_output = Path(output_str)
        if local_output.exists():
            logger.info(f"Output already exists: {output_str}")
            return output_str

    # Parse source scheme
    if source.startswith("cellxgene://"):
        adata = _acquire_from_cellxgene(source, **kwargs)

    elif source.startswith("synapse://"):
        adata = _acquire_from_synapse(source, **kwargs)

    else:
        raise ValueError(
            f"Unknown source scheme: {source}. "
            "Supported: cellxgene://, synapse://"
        )

    # Write to output (local or cloud)
    result_path = _write_output(adata, output)

    logger.info(f"Acquisition complete: {source} → {result_path}")
    return result_path


def _acquire_from_cellxgene(source: str, **kwargs) -> ad.AnnData:
    """
    Acquire data from CellxGene Census.

    Supports two formats:
    - cellxgene://dataset_key - Download predefined dataset
    - cellxgene://?tissue=lung&disease=normal - Query with filters
    """
    from urllib.parse import parse_qs, urlparse

    # Remove scheme
    rest = source[12:]  # len("cellxgene://") = 12

    if rest.startswith("?"):
        # Query format: cellxgene://?tissue=lung&disease=normal
        from spatialcore.annotation.cellxgene import query_cellxgene_census

        parsed = urlparse(source)
        params = parse_qs(parsed.query)

        # Extract query parameters
        tissue = params.get("tissue", [None])[0]
        disease = params.get("disease", [None])[0]
        cell_type = params.get("cell_type", [None])[0]
        assay = params.get("assay", [None])[0]

        logger.info(f"Querying CellxGene Census with filters...")
        logger.info(f"  tissue={tissue}, disease={disease}, cell_type={cell_type}")

        adata = query_cellxgene_census(
            tissue=tissue,
            disease=disease,
            cell_type=cell_type,
            assay=assay,
            max_cells=kwargs.get("max_cells"),
            random_state=kwargs.get("random_state", 42),
            resolve_hierarchy=kwargs.get("resolve_hierarchy", "none"),
            validate_labels=kwargs.get("validate_labels", True),
        )

    else:
        # Dataset key format: cellxgene://human_lung_cell_atlas
        from spatialcore.annotation.cellxgene import (
            download_cellxgene_reference,
            CELLXGENE_DATASETS,
        )

        dataset_key = rest

        if dataset_key not in CELLXGENE_DATASETS:
            available = ", ".join(CELLXGENE_DATASETS.keys())
            raise ValueError(
                f"Unknown CellxGene dataset: '{dataset_key}'. "
                f"Available: {available}. "
                "For custom queries, use: cellxgene://?tissue=lung&disease=normal"
            )

        # Download to temp location, then load
        with tempfile.TemporaryDirectory() as tmp_dir:
            tmp_path = download_cellxgene_reference(
                dataset_key=dataset_key,
                output_dir=tmp_dir,
                force=True,
            )
            adata = ad.read_h5ad(tmp_path)

    logger.info(f"Acquired from CellxGene: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    return adata


def _acquire_from_synapse(source: str, **kwargs) -> ad.AnnData:
    """Acquire data from Synapse."""
    from spatialcore.annotation.synapse import download_synapse_reference

    # Remove scheme: synapse://syn12345678 -> syn12345678
    synapse_id = source[10:]  # len("synapse://") = 10

    # Download to temp location
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_path = download_synapse_reference(
            synapse_id=synapse_id,
            output_dir=Path(tmp_dir),
            auth_token=kwargs.get("auth_token"),
            force=True,
        )
        adata = ad.read_h5ad(tmp_path)

    logger.info(f"Acquired from Synapse: {adata.n_obs:,} cells × {adata.n_vars:,} genes")
    return adata
