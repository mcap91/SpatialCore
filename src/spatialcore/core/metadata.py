"""Metadata tracking for AnnData operations."""

import json
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Optional, Union

import anndata as ad


class MetadataTracker:
    """
    Track operations performed on AnnData objects.

    Stores metadata in adata.uns['spatialcore_metadata'] and optionally
    writes to a JSON file.

    Parameters
    ----------
    adata
        AnnData object to track.
    json_path
        Optional path to write metadata JSON file.
    """

    def __init__(
        self,
        adata: ad.AnnData,
        json_path: Optional[Union[str, Path]] = None,
    ) -> None:
        self.adata = adata
        self.json_path = Path(json_path) if json_path else None

        if "spatialcore_metadata" not in adata.uns:
            adata.uns["spatialcore_metadata"] = {
                "created": datetime.now().isoformat(),
                "operations": [],
            }
        else:
            # Ensure operations is a list (may be numpy array after h5ad reload)
            meta = adata.uns["spatialcore_metadata"]
            if "operations" not in meta:
                meta["operations"] = []
            elif not isinstance(meta["operations"], list):
                meta["operations"] = list(meta["operations"])

    def log_operation(
        self,
        function_name: str,
        parameters: Dict[str, Any],
        outputs: Optional[Dict[str, str]] = None,
    ) -> None:
        """
        Log an operation to the metadata.

        Parameters
        ----------
        function_name
            Name of the function that was called.
        parameters
            Dictionary of parameters passed to the function.
        outputs
            Dictionary describing where outputs were stored
            (e.g., {'obs': 'cell_type', 'obsm': 'X_spatial_nmf'}).
        """
        entry = {
            "timestamp": datetime.now().isoformat(),
            "function": function_name,
            "parameters": _serialize_params(parameters),
        }
        if outputs:
            entry["outputs"] = outputs

        self.adata.uns["spatialcore_metadata"]["operations"].append(entry)

        if self.json_path:
            self._write_json()

    def _write_json(self) -> None:
        """Write metadata to JSON file."""
        with open(self.json_path, "w") as f:
            json.dump(self.adata.uns["spatialcore_metadata"], f, indent=2)

    def get_history(self) -> list:
        """Return list of all operations performed."""
        return self.adata.uns["spatialcore_metadata"]["operations"]


def update_metadata(
    adata: ad.AnnData,
    function_name: str,
    parameters: Dict[str, Any],
    outputs: Optional[Dict[str, str]] = None,
) -> None:
    """
    Convenience function to update metadata without creating a tracker.

    Parameters
    ----------
    adata
        AnnData object to update.
    function_name
        Name of the function that was called.
    parameters
        Dictionary of parameters passed to the function.
    outputs
        Dictionary describing where outputs were stored.
    """
    tracker = MetadataTracker(adata)
    tracker.log_operation(function_name, parameters, outputs)


def prepare_metadata_for_h5ad(adata: ad.AnnData) -> None:
    """
    Convert metadata to h5ad-compatible format.

    The operations list contains nested dicts that can't be serialized
    directly to h5ad. This function converts it to a JSON string.

    Parameters
    ----------
    adata
        AnnData object with spatialcore_metadata.
    """
    if "spatialcore_metadata" in adata.uns:
        meta = adata.uns["spatialcore_metadata"]
        if "operations" in meta and isinstance(meta["operations"], list):
            # Convert operations list to JSON string for h5ad compatibility
            meta["operations_json"] = json.dumps(meta["operations"])
            meta["operations"] = []  # Clear the list to avoid serialization issues


def _serialize_params(params: Dict[str, Any]) -> Dict[str, Any]:
    """Convert parameters to JSON-serializable format."""
    serialized = {}
    for key, value in params.items():
        if isinstance(value, (str, int, float, bool, type(None))):
            serialized[key] = value
        elif isinstance(value, (list, tuple)):
            serialized[key] = list(value)
        elif isinstance(value, dict):
            serialized[key] = _serialize_params(value)
        elif isinstance(value, Path):
            serialized[key] = str(value)
        elif isinstance(value, ad.AnnData):
            serialized[key] = f"<AnnData: {value.shape[0]} obs x {value.shape[1]} var>"
        else:
            serialized[key] = str(type(value).__name__)
    return serialized
