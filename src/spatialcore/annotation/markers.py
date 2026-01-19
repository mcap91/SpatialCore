"""
Canonical markers for cell type validation.

This module provides canonical marker gene definitions for common cell types.
These markers are used for validation of cell type annotations, typically
via GMM-3 thresholding (see spatialcore.stats.classify).

Marker genes are curated from literature and validated on spatial
transcriptomics platforms (Xenium, CosMx).

References:
    - Domínguez Conde et al., Science (2022) - Immune cell markers
    - Tabula Sapiens Consortium (2022) - Pan-tissue markers
    - Human Cell Atlas marker databases
"""

from pathlib import Path
from typing import Dict, List, Optional
import json

from spatialcore.core.logging import get_logger

logger = get_logger(__name__)

# Default path for canonical markers (package data directory)
DEFAULT_MARKERS_PATH = Path(__file__).parent.parent / "data" / "markers" / "canonical_markers.json"


# ============================================================================
# Canonical Marker Definitions
# ============================================================================
# NOTE: All canonical markers are now defined in a single source of truth:
#       src/spatialcore/data/markers/canonical_markers.json
#
# Use load_canonical_markers() to access them. The JSON file contains 75+
# cell types with curated marker genes from literature.
# ============================================================================


# ============================================================================
# Marker Loading and Lookup
# ============================================================================

def load_canonical_markers(
    config_path: Optional[Path] = None,
) -> Dict[str, List[str]]:
    """
    Load canonical markers from JSON config.

    Parameters
    ----------
    config_path : Path, optional
        Path to JSON file with custom marker definitions. If None,
        loads from the default canonical_markers.json in the package
        data directory.

        JSON can be either:
        - Simple format: ``{"cell_type": ["GENE1", "GENE2"]}``
        - Extended format: ``{"cell_type": {"index_marker": ["GENE1"], "description": "..."}}``

    Returns
    -------
    Dict[str, List[str]]
        Dictionary mapping cell type names to marker gene lists.

    Raises
    ------
    FileNotFoundError
        If the markers JSON file does not exist.

    Notes
    -----
    Cell type names should be in lowercase and match Cell Ontology (CL)
    naming conventions where possible (e.g., "cd4-positive, alpha-beta t cell").

    All canonical markers are defined in a single source of truth:
    ``src/spatialcore/data/markers/canonical_markers.json``

    The function supports two JSON formats:

    Simple format (list of genes):

    .. code-block:: json

        {
            "my custom type": ["GENE1", "GENE2", "GENE3"],
            "another type": ["GENE4", "GENE5"]
        }

    Extended format (with metadata):

    .. code-block:: json

        {
            "my custom type": {
                "index_marker": ["GENE1", "GENE2"],
                "description": "Description text"
            }
        }

    Examples
    --------
    >>> from spatialcore.annotation.markers import load_canonical_markers
    >>> markers = load_canonical_markers()
    >>> print(markers["macrophage"])
    ['CD163', 'CD68', 'MARCO', 'CSF1R', 'MERTK', 'C1QA', 'C1QB', 'C1QC', 'MRC1']
    >>> # Load custom markers from a different file
    >>> markers = load_canonical_markers(Path("custom_markers.json"))
    """
    # Determine which file to load
    markers_path = config_path if config_path is not None else DEFAULT_MARKERS_PATH

    if not markers_path.exists():
        raise FileNotFoundError(
            f"Canonical markers file not found: {markers_path}. "
            "Ensure the package data directory contains canonical_markers.json"
        )

    markers = _load_markers_from_json(markers_path)
    logger.debug(f"Loaded {len(markers)} markers from {markers_path}")

    return markers


def _load_markers_from_json(path: Path) -> Dict[str, List[str]]:
    """
    Load markers from a JSON file, supporting multiple formats.

    Supported formats:
    1. Simple: ``{"cell_type": ["GENE1", "GENE2"]}``
    2. Extended: ``{"cell_type": {"index_marker": ["GENE1"], ...}}``
    3. Wrapped: ``{"metadata": {...}, "markers": {"cell_type": [...]}}``

    Parameters
    ----------
    path : Path
        Path to JSON file.

    Returns
    -------
    Dict[str, List[str]]
        Dictionary mapping cell type names to marker gene lists.
    """
    with open(path, "r", encoding="utf-8") as f:
        data = json.load(f)

    # Handle wrapped format with "markers" key
    if "markers" in data and isinstance(data["markers"], dict):
        data = data["markers"]

    markers = {}
    for key, value in data.items():
        # Skip metadata keys
        if key.startswith("_") or key == "metadata":
            continue

        if isinstance(value, list):
            # Simple format: {"cell_type": ["GENE1", "GENE2"]}
            if not all(isinstance(g, str) for g in value):
                raise ValueError(f"All marker genes must be strings for '{key}'")
            markers[key] = value
        elif isinstance(value, dict):
            # Extended format: {"cell_type": {"index_marker": [...], "description": "..."}}
            if "index_marker" in value:
                gene_list = value["index_marker"]
                if not isinstance(gene_list, list):
                    raise ValueError(
                        f"'index_marker' must be a list for '{key}'"
                    )
                if not all(isinstance(g, str) for g in gene_list):
                    raise ValueError(
                        f"All marker genes must be strings for '{key}'"
                    )
                markers[key] = gene_list
            else:
                logger.warning(
                    f"Skipping '{key}': dict format requires 'index_marker' key"
                )
        else:
            raise ValueError(
                f"Marker config values must be lists or dicts, got {type(value)} for '{key}'"
            )

    return markers


def match_to_canonical(
    cell_type: str,
    markers: Optional[Dict[str, List[str]]] = None,
) -> Optional[str]:
    """
    Match a cell type name to a canonical cell type in the markers dictionary.

    Uses exact case-insensitive matching only. No fuzzy/substring matching.

    Parameters
    ----------
    cell_type : str
        Cell type name to match.
    markers : Dict[str, List[str]], optional
        Marker dictionary. If None, loads from canonical_markers.json.

    Returns
    -------
    str or None
        Canonical cell type name if matched, None otherwise.

    Examples
    --------
    >>> from spatialcore.annotation.markers import match_to_canonical
    >>> match_to_canonical("Macrophage")
    'macrophage'
    >>> match_to_canonical("B cell")
    'b cell'
    >>> match_to_canonical("some unknown type")
    None
    """
    if cell_type in ["Unassigned", "Unknown", "unknown", "cell", None, ""]:
        return None

    if markers is None:
        markers = load_canonical_markers()

    cell_type_lower = cell_type.lower().strip()

    for canonical_name in markers.keys():
        if canonical_name.lower().strip() == cell_type_lower:
            return canonical_name

    return None


def get_markers_for_type(
    cell_type: str,
    markers: Optional[Dict[str, List[str]]] = None,
) -> List[str]:
    """
    Get marker genes for a specific cell type.

    Uses exact case-insensitive matching only.

    Parameters
    ----------
    cell_type : str
        Cell type name to look up.
    markers : Dict[str, List[str]], optional
        Marker dictionary. If None, loads from canonical_markers.json.

    Returns
    -------
    List[str]
        Marker gene list. Empty list if no match found.

    Examples
    --------
    >>> from spatialcore.annotation.markers import get_markers_for_type
    >>> markers = get_markers_for_type("Macrophage")
    >>> print(markers)
    ['CD68', 'CD163', 'CD14', 'FCGR3A', 'CSF1R']
    """
    if markers is None:
        markers = load_canonical_markers()

    matched = match_to_canonical(cell_type, markers)
    if matched is not None:
        return markers[matched]

    return []


def list_available_cell_types(
    markers: Optional[Dict[str, List[str]]] = None,
) -> List[str]:
    """
    List all cell types with defined markers.

    Parameters
    ----------
    markers : Dict[str, List[str]], optional
        Marker dictionary. If None, loads from canonical_markers.json.

    Returns
    -------
    List[str]
        Sorted list of cell type names.

    Examples
    --------
    >>> from spatialcore.annotation.markers import list_available_cell_types
    >>> types = list_available_cell_types()
    >>> print(f"Available types: {len(types)}")
    >>> print(types[:5])
    """
    if markers is None:
        markers = load_canonical_markers()

    return sorted(markers.keys())
