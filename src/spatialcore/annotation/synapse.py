"""
Synapse data download utilities.

This module provides functions for downloading reference datasets from
Synapse (synapse.org), a platform commonly used for sharing biomedical
research data.

Authentication requires a Synapse account and auth token. Set the
SYNAPSE_AUTH_TOKEN environment variable or pass the token directly.

References:
    - Synapse: https://www.synapse.org/
    - synapseclient docs: https://python-docs.synapse.org/
"""

import os
from pathlib import Path
from typing import Dict, List, Optional, Any

from spatialcore.core.logging import get_logger

logger = get_logger(__name__)


def authenticate_synapse(auth_token: Optional[str] = None) -> bool:
    """
    Authenticate with Synapse.

    Parameters
    ----------
    auth_token : str, optional
        Synapse authentication token. If None, uses SYNAPSE_AUTH_TOKEN
        environment variable.

    Returns
    -------
    bool
        True if authentication successful.

    Notes
    -----
    To get a Synapse auth token:
    1. Create account at synapse.org
    2. Go to Account Settings > Personal Access Tokens
    3. Generate token with "Download" scope
    4. Set SYNAPSE_AUTH_TOKEN environment variable

    Examples
    --------
    >>> from spatialcore.annotation.synapse import authenticate_synapse
    >>> # Using environment variable
    >>> success = authenticate_synapse()
    >>> # Using direct token
    >>> success = authenticate_synapse("my_token")
    """
    try:
        import synapseclient
    except ImportError:
        raise ImportError(
            "synapseclient is required for Synapse downloads. "
            "Install with: pip install synapseclient"
        )

    token = auth_token or os.environ.get("SYNAPSE_AUTH_TOKEN")

    if not token:
        logger.warning(
            "No Synapse auth token provided. "
            "Set SYNAPSE_AUTH_TOKEN environment variable or pass auth_token."
        )
        return False

    try:
        syn = synapseclient.Synapse()
        syn.login(authToken=token, silent=True)
        logger.info("Successfully authenticated with Synapse")
        return True
    except Exception as e:
        logger.error(f"Synapse authentication failed: {e}")
        return False


def download_synapse_reference(
    synapse_id: str,
    output_dir: Path,
    auth_token: Optional[str] = None,
    force: bool = False,
) -> Path:
    """
    Download reference dataset from Synapse.

    Parameters
    ----------
    synapse_id : str
        Synapse entity ID (e.g., "syn12345678").
    output_dir : Path
        Directory to save downloaded file.
    auth_token : str, optional
        Synapse authentication token. If None, uses SYNAPSE_AUTH_TOKEN
        environment variable.
    force : bool, default False
        Force re-download even if file exists.

    Returns
    -------
    Path
        Path to downloaded file.

    Raises
    ------
    ImportError
        If synapseclient is not installed.
    ValueError
        If authentication fails or entity not found.

    Examples
    --------
    >>> from spatialcore.annotation.synapse import download_synapse_reference
    >>> path = download_synapse_reference(
    ...     "syn12345678",
    ...     output_dir=Path("./data"),
    ... )
    >>> print(f"Downloaded to: {path}")
    """
    try:
        import synapseclient
    except ImportError:
        raise ImportError(
            "synapseclient is required for Synapse downloads. "
            "Install with: pip install synapseclient"
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Get auth token
    token = auth_token or os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        raise ValueError(
            "Synapse authentication required. "
            "Set SYNAPSE_AUTH_TOKEN environment variable or pass auth_token."
        )

    # Initialize client
    syn = synapseclient.Synapse()
    try:
        syn.login(authToken=token, silent=True)
    except Exception as e:
        raise ValueError(f"Synapse authentication failed: {e}")

    # Get entity info
    try:
        entity = syn.get(synapse_id, downloadFile=False)
        filename = entity.name
        expected_path = output_dir / filename
    except Exception as e:
        raise ValueError(f"Failed to get Synapse entity '{synapse_id}': {e}")

    # Check if already downloaded
    if expected_path.exists() and not force:
        logger.info(f"File already exists: {expected_path}")
        return expected_path

    # Download
    logger.info(f"Downloading Synapse entity: {synapse_id}")
    logger.info(f"  Name: {filename}")

    try:
        entity = syn.get(synapse_id, downloadLocation=str(output_dir))
        downloaded_path = Path(entity.path)
        logger.info(f"  Downloaded to: {downloaded_path}")
        return downloaded_path
    except Exception as e:
        raise ValueError(f"Failed to download '{synapse_id}': {e}")


def get_synapse_entity_info(
    synapse_id: str,
    auth_token: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Get metadata for a Synapse entity.

    Parameters
    ----------
    synapse_id : str
        Synapse entity ID.
    auth_token : str, optional
        Synapse authentication token.

    Returns
    -------
    Dict[str, Any]
        Entity metadata including:
        - name: File name
        - id: Synapse ID
        - content_type: MIME type
        - size_mb: File size in MB
        - md5: MD5 checksum
        - created_on: Creation date
        - modified_on: Last modified date

    Examples
    --------
    >>> from spatialcore.annotation.synapse import get_synapse_entity_info
    >>> info = get_synapse_entity_info("syn12345678")
    >>> print(f"Size: {info['size_mb']:.1f} MB")
    """
    try:
        import synapseclient
    except ImportError:
        raise ImportError(
            "synapseclient is required. "
            "Install with: pip install synapseclient"
        )

    token = auth_token or os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        raise ValueError("Synapse authentication required.")

    syn = synapseclient.Synapse()
    syn.login(authToken=token, silent=True)

    entity = syn.get(synapse_id, downloadFile=False)

    # Extract metadata
    info = {
        "name": entity.name,
        "id": entity.id,
        "content_type": getattr(entity, "contentType", "unknown"),
        "created_on": getattr(entity, "createdOn", None),
        "modified_on": getattr(entity, "modifiedOn", None),
    }

    # Get file-specific info
    if hasattr(entity, "contentSize"):
        info["size_mb"] = entity.contentSize / (1024 * 1024)
    if hasattr(entity, "md5"):
        info["md5"] = entity.md5

    return info


def list_synapse_folder(
    folder_id: str,
    auth_token: Optional[str] = None,
    file_types: Optional[List[str]] = None,
) -> List[Dict[str, Any]]:
    """
    List contents of a Synapse folder.

    Parameters
    ----------
    folder_id : str
        Synapse folder ID.
    auth_token : str, optional
        Synapse authentication token.
    file_types : List[str], optional
        Filter by file extensions (e.g., [".h5ad", ".h5"]).

    Returns
    -------
    List[Dict[str, Any]]
        List of entity metadata dictionaries.

    Examples
    --------
    >>> from spatialcore.annotation.synapse import list_synapse_folder
    >>> files = list_synapse_folder("syn12345", file_types=[".h5ad"])
    >>> for f in files:
    ...     print(f"{f['name']}: {f['size_mb']:.1f} MB")
    """
    try:
        import synapseclient
    except ImportError:
        raise ImportError(
            "synapseclient is required. "
            "Install with: pip install synapseclient"
        )

    token = auth_token or os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        raise ValueError("Synapse authentication required.")

    syn = synapseclient.Synapse()
    syn.login(authToken=token, silent=True)

    # Get folder children
    children = list(syn.getChildren(folder_id))

    results = []
    for child in children:
        # Filter by type if specified
        if file_types is not None:
            name = child.get("name", "")
            if not any(name.endswith(ext) for ext in file_types):
                continue

        info = {
            "name": child.get("name"),
            "id": child.get("id"),
            "type": child.get("type"),
        }

        # Add size if available
        if "dataFileHandleId" in child:
            try:
                entity = syn.get(child["id"], downloadFile=False)
                if hasattr(entity, "contentSize"):
                    info["size_mb"] = entity.contentSize / (1024 * 1024)
            except Exception:
                pass

        results.append(info)

    return results


def download_synapse_folder(
    folder_id: str,
    output_dir: Path,
    auth_token: Optional[str] = None,
    file_types: Optional[List[str]] = None,
    force: bool = False,
) -> List[Path]:
    """
    Download all files from a Synapse folder.

    Parameters
    ----------
    folder_id : str
        Synapse folder ID.
    output_dir : Path
        Directory to save downloaded files.
    auth_token : str, optional
        Synapse authentication token.
    file_types : List[str], optional
        Only download files with these extensions.
    force : bool, default False
        Force re-download even if files exist.

    Returns
    -------
    List[Path]
        Paths to downloaded files.

    Examples
    --------
    >>> from spatialcore.annotation.synapse import download_synapse_folder
    >>> paths = download_synapse_folder(
    ...     "syn12345",
    ...     output_dir=Path("./data"),
    ...     file_types=[".h5ad"],
    ... )
    """
    # List folder contents
    files = list_synapse_folder(folder_id, auth_token, file_types)

    if not files:
        logger.warning(f"No matching files found in folder {folder_id}")
        return []

    logger.info(f"Found {len(files)} files to download")

    downloaded = []
    for file_info in files:
        if file_info.get("type") == "org.sagebionetworks.repo.model.FileEntity":
            try:
                path = download_synapse_reference(
                    file_info["id"],
                    output_dir,
                    auth_token,
                    force,
                )
                downloaded.append(path)
            except Exception as e:
                logger.warning(f"Failed to download {file_info['name']}: {e}")

    logger.info(f"Downloaded {len(downloaded)} files")
    return downloaded
