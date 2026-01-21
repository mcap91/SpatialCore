"""
Ontology mapping: Convert cell type labels to Cell Ontology (CL) codes.

This module provides a 4-tier matching system:
1. Tier 0: Pattern canonicalization (known abbreviations → CL term names)
2. Tier 1: Exact match (score 1.0)
3. Tier 2: Token-based match (score 0.60-0.85)
4. Tier 3: Word overlap fallback (score 0.5-0.7)

The system uses a pre-built JSON index for fast, offline operation.

References:
    - Cell Ontology (CL): https://github.com/obophenotype/cell-ontology
    - NCI Thesaurus (NCIT): https://ncithesaurus.nci.nih.gov/
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any, Union
from datetime import datetime
import re
import json

import numpy as np
import pandas as pd
import anndata as ad

from spatialcore.core.logging import get_logger
from spatialcore.annotation.patterns import CELL_TYPE_PATTERNS, get_canonical_term

logger = get_logger(__name__)


# ============================================================================
# Constants: Unknown Cell Type Definition
# ============================================================================

# For unmapped labels that don't match any ontology term, we use a
# standardized "Unknown" category. CL has no native "unknown cell" term.
UNKNOWN_CELL_TYPE_ID = "unknown"
UNKNOWN_CELL_TYPE_NAME = "Unknown"


# ============================================================================
# Data Classes for Structured Results
# ============================================================================

@dataclass
class OntologyMappingResult:
    """
    Complete result of an ontology mapping operation.

    Contains the mapping table, metadata, and any errors encountered.
    This is the primary output structure from create_mapping_table().

    Attributes
    ----------
    table : pd.DataFrame
        Mapping table with columns: input_label, ontology_name, ontology_id,
        match_tier, score, n_cells, canonical_term
    metadata : Dict[str, Any]
        Full metadata dictionary for JSON serialization
    errors : List[Dict[str, Any]]
        List of errors/warnings encountered during mapping
    """

    table: pd.DataFrame
    metadata: Dict[str, Any] = field(default_factory=dict)
    errors: List[Dict[str, Any]] = field(default_factory=list)

    def to_json(self, path: Union[str, Path]) -> Path:
        """
        Save metadata to JSON file.

        Parameters
        ----------
        path : str or Path
            Output path for JSON file.

        Returns
        -------
        Path
            Path to saved file.
        """
        path = Path(path)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(self.metadata, f, indent=2, default=str)
        return path

    def to_csv(self, path: Union[str, Path]) -> Path:
        """
        Save mapping table to CSV.

        Parameters
        ----------
        path : str or Path
            Output path for CSV file.

        Returns
        -------
        Path
            Path to saved file.
        """
        path = Path(path)
        self.table.to_csv(path, index=False)
        return path


# ============================================================================
# Index Loading
# ============================================================================

_ONTOLOGY_INDEX_CACHE: Optional[Dict] = None


def load_ontology_index(
    index_path: Optional[Union[str, Path]] = None,
    use_cache: bool = True,
) -> Dict[str, Dict[str, Dict[str, str]]]:
    """
    Load pre-built ontology index from JSON file.

    Parameters
    ----------
    index_path : str or Path, optional
        Path to ontology_index.json. If None, uses default location.
    use_cache : bool, default True
        Cache the loaded index for faster subsequent calls.

    Returns
    -------
    Dict[str, Dict[str, Dict[str, str]]]
        Nested dictionary: {ontology: {label_lower: {id, name}}}
        - ontology: "cl", "ncit", or "uberon"
        - label_lower: lowercase term name
        - id: ontology ID (e.g., "CL:0000624")
        - name: canonical term name

    Examples
    --------
    >>> from spatialcore.annotation import load_ontology_index
    >>> index = load_ontology_index()
    >>> index["cl"]["b cell"]
    {'id': 'CL:0000236', 'name': 'B cell'}
    """
    global _ONTOLOGY_INDEX_CACHE

    if use_cache and _ONTOLOGY_INDEX_CACHE is not None:
        return _ONTOLOGY_INDEX_CACHE

    if index_path is None:
        # Default: look in package data directory first, then fallback locations
        possible_paths = []

        # Primary: Package data directory
        package_data_path = Path(__file__).parent.parent / "data" / "ontology_mappings" / "ontology_index.json"
        possible_paths.append(package_data_path)

        # Fallback: User cache
        possible_paths.append(
            Path.home() / ".cache" / "spatialcore" / "ontology_index.json"
        )

        for path in possible_paths:
            if path.exists():
                index_path = path
                break
        else:
            raise FileNotFoundError(
                f"Ontology index not found. Searched: {[str(p) for p in possible_paths]}. "
                "Ensure spatialcore is installed correctly with data files."
            )
    else:
        index_path = Path(index_path)

    logger.info(f"Loading ontology index from: {index_path}")

    with open(index_path, "r", encoding="utf-8") as f:
        raw_index = json.load(f)

    # Extract term dictionaries (skip metadata)
    index = {
        "cl": raw_index.get("cl", {}),
        "ncit": raw_index.get("ncit", {}),
        "uberon": raw_index.get("uberon", {}),
    }

    # Log stats
    if "metadata" in raw_index:
        meta = raw_index["metadata"]
        logger.info(
            f"  CL: {meta.get('cl_terms', len(index['cl'])):,} terms, "
            f"NCIT: {meta.get('ncit_terms', len(index['ncit'])):,} terms, "
            f"UBERON: {meta.get('uberon_terms', len(index['uberon'])):,} terms"
        )

    if use_cache:
        _ONTOLOGY_INDEX_CACHE = index

    return index


# ============================================================================
# Token Extraction
# ============================================================================

# Words that are too generic for meaningful matching
GENERIC_TERMS = {"cell", "cells", "type", "like"}

# CL IDs that are too generic to be valid match results
# These terms exist in the ontology but should never be returned as matches
# because they provide no useful information about cell type identity
BLACKLISTED_CL_IDS = {
    "CL:0000000",  # "cell" - root term, too generic
    "CL:0000003",  # "native cell" - too generic
    "CL:0000255",  # "eukaryotic cell" - too generic
}

# Modifiers that describe state, not identity
MODIFIER_TERMS = {
    "positive", "negative", "high", "low", "like", "type",
    "mature", "immature", "activated", "resting", "proliferating",
    "pro", "pre", "post", "inflammatory", "naive", "memory",
    "effector", "resident", "circulating",
}

# Short tokens that ARE meaningful for cell types
MEANINGFUL_SHORT_TOKENS = {
    "b", "t", "nk", "dc", "ec", "ve", "ta",
    "m1", "m2", "cd", "th", "ilc",
}


def extract_biological_tokens(label: str) -> Dict[str, List[str]]:
    """
    Extract key biological identifiers from a cell type label.

    Parameters
    ----------
    label : str
        Cell type label to tokenize.

    Returns
    -------
    Dict[str, List[str]]
        Dictionary with keys:
        - markers: CD markers (cd4, cd8, cd19, ...)
        - proteins: Immunoglobulins, gene names (igg, iga, spp1, ...)
        - core_words: Main biological terms (helper, plasma, ...)
        - modifiers: Descriptors (positive, mature, ...)

    Examples
    --------
    >>> tokens = extract_biological_tokens("CD4+ T cells")
    >>> tokens["markers"]
    ['cd4']
    >>> tokens["core_words"]
    ['t']
    """
    label_lower = label.lower().strip()
    tokens = {
        "markers": [],
        "proteins": [],
        "core_words": [],
        "modifiers": [],
    }

    # Extract CD markers: CD4, CD8, CD19, etc.
    cd_markers = re.findall(r"cd\d+", label_lower)
    tokens["markers"].extend(cd_markers)

    # Extract immunoglobulin types: IgG, IgA, IgM
    ig_types = re.findall(r"ig[gamedGAMED]", label_lower)
    tokens["proteins"].extend([ig.lower() for ig in ig_types])

    # Extract gene names (uppercase + plus sign pattern)
    gene_markers = re.findall(r"\b[A-Z0-9]{3,}\+", label)
    tokens["proteins"].extend([g.replace("+", "").lower() for g in gene_markers])

    # Clean and extract core words
    cleaned = re.sub(r"cd\d+", "", label_lower)  # Remove CD markers
    cleaned = re.sub(r"ig[gamed]", "", cleaned)   # Remove Ig types
    cleaned = re.sub(r"[+\-]", " ", cleaned)      # Replace +/- with space
    cleaned = re.sub(r"\d+", "", cleaned)         # Remove numbers
    cleaned = re.sub(r"\s+", " ", cleaned).strip()

    for word in cleaned.split():
        if word in MODIFIER_TERMS:
            tokens["modifiers"].append(word)
        elif word in GENERIC_TERMS:
            pass  # Skip generic terms
        elif word in MEANINGFUL_SHORT_TOKENS:
            tokens["core_words"].append(word)
        elif len(word) > 1:
            tokens["core_words"].append(word)

    return tokens


# ============================================================================
# Scoring Functions
# ============================================================================

def _score_match(
    search_label: str,
    term_label: str,
    tokens: Dict[str, List[str]],
    is_pattern_match: bool,
) -> float:
    """
    Calculate match score between search label and ontology term.

    Parameters
    ----------
    search_label : str
        Canonicalized search label.
    term_label : str
        Ontology term label (lowercase).
    tokens : Dict
        Extracted tokens from search label.
    is_pattern_match : bool
        Whether search_label came from pattern canonicalization.

    Returns
    -------
    float
        Match score (0.0 to 1.0).
    """
    search_clean = search_label.lower().strip()
    term_clean = term_label.lower().strip()

    # Tier 1: Exact match
    if search_clean == term_clean:
        return 0.95 if is_pattern_match else 1.0

    # Tier 1b: Clean match (remove symbols)
    search_no_symbols = re.sub(r"[+\-,]", " ", search_clean)
    search_no_symbols = re.sub(r"\s+", " ", search_no_symbols).strip()
    term_no_symbols = re.sub(r"[+\-,]", " ", term_clean)
    term_no_symbols = re.sub(r"\s+", " ", term_no_symbols).strip()

    if search_no_symbols == term_no_symbols:
        return 0.92 if is_pattern_match else 0.95

    # Tier 1c: Word boundary contains match (avoid false positives)
    if len(search_clean) >= 4:
        # Only match if it's a word boundary match, not arbitrary substring
        if re.search(rf'\b{re.escape(search_clean)}\b', term_clean):
            return 0.88 if is_pattern_match else 0.90
        elif re.search(rf'\b{re.escape(term_clean)}\b', search_clean):
            return 0.86 if is_pattern_match else 0.88

    # Tier 2: Token-based matching
    term_words = set(term_clean.replace("-", " ").replace(",", " ").split())
    core_words = tokens.get("core_words", [])
    markers = tokens.get("markers", [])

    if core_words:
        # Check if all core words present (exact word match only)
        matches_all_core = all(
            any(word == tw for tw in term_words)
            for word in core_words
        )

        if matches_all_core:
            base_score = 0.70

            # Penalty for single short token (too ambiguous)
            if len(core_words) == 1 and len(core_words[0]) <= 2:
                base_score -= 0.15

            # Penalty for unwanted prefixes in term but not in label
            unwanted_prefixes = ["pro", "pre", "post", "immature", "ecto", "endo"]
            label_has_prefix = any(p in search_clean for p in unwanted_prefixes)
            term_has_prefix = any(p in term_clean for p in unwanted_prefixes)
            if term_has_prefix and not label_has_prefix:
                base_score -= 0.15

            # Bonus if markers also match
            if markers and any(m in term_clean for m in markers):
                base_score = max(base_score, 0.75)

            # Bonus for word count similarity
            if len(core_words) >= 2:
                base_score = min(base_score + 0.05, 0.85)

            return max(base_score, 0.0)

    # Tier 3: Word overlap (Jaccard similarity)
    label_words = set(search_clean.replace("-", " ").replace(",", " ").split())
    label_words -= GENERIC_TERMS

    if label_words and term_words:
        common = label_words & term_words
        jaccard = len(common) / len(label_words | term_words)
        score = 0.5 + (0.4 * jaccard)  # Range: 0.5-0.9
        return score

    return 0.0


# ============================================================================
# Main Search Function
# ============================================================================

def search_ontology_index(
    labels: List[str],
    ontology_index: Optional[Dict] = None,
    index_path: Optional[Union[str, Path]] = None,
    annotation_type: str = "cell_type",
    min_score: float = 0.7,
) -> Dict[str, List[Dict[str, Any]]]:
    """
    Search ontology index for matching terms.

    Uses a 4-tier matching system:
    1. Tier 0: Pattern canonicalization
    2. Tier 1: Exact/partial match
    3. Tier 2: Token-based match
    4. Tier 3: Word overlap fallback

    Parameters
    ----------
    labels : List[str]
        Cell type labels to search.
    ontology_index : Dict, optional
        Pre-loaded ontology index. If None, loads from file.
    index_path : str or Path, optional
        Path to ontology index JSON.
    annotation_type : str, default "cell_type"
        Type of annotation:
        - "cell_type": Search CL only (most specific)
        - "pathology": Search NCIT first, then CL
        - "anatomy": Search UBERON first, then CL
        - "all": Search CL, NCIT, UBERON
    min_score : float, default 0.7
        Minimum match score to accept.

    Returns
    -------
    Dict[str, List[Dict[str, Any]]]
        {label: [{id, name, ontology, score, match_type}, ...]}

    Examples
    --------
    >>> from spatialcore.annotation import search_ontology_index
    >>> results = search_ontology_index(["CD4+ T cells", "B cell", "NK cells"])
    >>> results["CD4+ T cells"][0]
    {'id': 'CL:0000624', 'name': 'CD4-positive, alpha-beta T cell', 'score': 0.95, ...}
    """
    if ontology_index is None:
        ontology_index = load_ontology_index(index_path)

    # Select ontologies based on annotation type
    if annotation_type == "cell_type":
        ontologies = ["cl"]
    elif annotation_type == "pathology":
        ontologies = ["ncit", "cl"]
    elif annotation_type == "anatomy":
        ontologies = ["uberon", "cl"]
    else:
        ontologies = ["cl", "ncit", "uberon"]

    results = {label: [] for label in labels}

    for label in labels:
        label_lower = label.lower().strip().replace("_", " ")
        label_normalized = label.replace("_", " ")

        # Tier 0: Pattern canonicalization
        canonical_term = get_canonical_term(label_normalized)
        search_label = canonical_term if canonical_term else label_lower
        is_pattern_match = canonical_term is not None

        # Extract tokens for fuzzy matching
        tokens = extract_biological_tokens(
            canonical_term if canonical_term else label_normalized
        )

        # Search across ontologies
        tier_results = []

        for onto_prefix in ontologies:
            onto_dict = ontology_index.get(onto_prefix.lower(), {})
            ontology_matches = []

            # Tier 1: Exact lookup
            if search_label in onto_dict:
                term = onto_dict[search_label]
                ontology_matches.append({
                    "id": term["id"],
                    "name": term["name"],
                    "ontology": onto_prefix,
                    "score": 0.95 if is_pattern_match else 1.0,
                    "match_type": "tier0_pattern" if is_pattern_match else "tier1_exact",
                })
            else:
                # Tier 2-3: Fuzzy matching
                for term_label_lower, term_data in onto_dict.items():
                    # Skip imported terms (e.g., UBERON term in CL)
                    term_id_prefix = term_data["id"].split(":")[0].upper()
                    if term_id_prefix != onto_prefix.upper():
                        continue

                    # Skip obsolete terms
                    if "obsolete" in term_data["name"].lower():
                        continue

                    # Skip blacklisted generic terms (e.g., CL:0000000 "cell")
                    if term_data["id"] in BLACKLISTED_CL_IDS:
                        continue

                    score = _score_match(search_label, term_label_lower, tokens, is_pattern_match)

                    if score >= min_score:
                        ontology_matches.append({
                            "id": term_data["id"],
                            "name": term_data["name"],
                            "ontology": onto_prefix,
                            "score": score,
                            "match_type": "tier2_token" if score >= 0.7 else "tier3_overlap",
                        })

            if ontology_matches:
                tier_results.extend(ontology_matches)
                # If good CL match found, don't search other ontologies
                if onto_prefix == "cl" and any(m["score"] >= 0.8 for m in ontology_matches):
                    break

        # Sort by score, deduplicate by ID
        seen_ids = set()
        unique_results = []
        for result in sorted(tier_results, key=lambda x: x["score"], reverse=True):
            if result["id"] not in seen_ids:
                seen_ids.add(result["id"])
                unique_results.append(result)

        results[label] = unique_results

    return results


# ============================================================================
# Mapping Result Saving
# ============================================================================

def _save_ontology_mapping_results(
    save_dir: Path,
    mappings: Dict[str, List[Dict]],
    adata: ad.AnnData,
    source_col: str,
    dataset_name: str,
) -> Tuple[Path, Optional[Path]]:
    """
    Save ontology mapping results for reproducibility.

    Creates two files:
    - {dataset}_{timestamp}_mapping.json: Full mapping results
    - {dataset}_{timestamp}_missed.json: Unmapped terms (if any)

    Parameters
    ----------
    save_dir : Path
        Directory to save results.
    mappings : Dict
        Mapping results from search_ontology_index.
    adata : AnnData
        Source data for cell counts.
    source_col : str
        Source column name.
    dataset_name : str
        Name for output files.

    Returns
    -------
    Tuple[Path, Optional[Path]]
        Paths to mapping.json and missed.json (or None if no missed).
    """
    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base_name = f"{dataset_name}_{timestamp}"

    try:
        import spatialcore
        spatialcore_version = spatialcore.__version__
    except Exception:
        spatialcore_version = "unknown"

    # Count cells per label
    cell_counts = adata.obs[source_col].value_counts().to_dict()

    # Build mapping results
    mapped_entries = []
    missed_entries = []
    tier_counts = {"pattern": 0, "exact": 0, "token": 0, "overlap": 0, "unmapped": 0}

    for label, matches in mappings.items():
        n_cells = int(cell_counts.get(label, 0))

        if matches:
            best = matches[0]
            tier = best.get("match_type", "unknown")

            # Simplify tier names for counting
            if "pattern" in tier:
                tier_counts["pattern"] += 1
            elif "exact" in tier:
                tier_counts["exact"] += 1
            elif "token" in tier:
                tier_counts["token"] += 1
            else:
                tier_counts["overlap"] += 1

            canonical = get_canonical_term(label)
            mapped_entries.append({
                "input_label": label,
                "canonical_term": canonical,
                "ontology_id": best["id"],
                "ontology_name": best["name"],
                "match_tier": tier,
                "score": round(best["score"], 3),
                "n_cells": n_cells,
            })
        else:
            tier_counts["unmapped"] += 1
            # Get closest matches for missed terms
            closest = search_ontology_index([label], min_score=0.4).get(label, [])[:3]
            missed_entries.append({
                "input_label": label,
                "n_cells": n_cells,
                "closest_matches": [
                    {"term": m["name"], "id": m["id"], "score": round(m["score"], 2)}
                    for m in closest
                ],
                "suggested_action": "manual_review",
            })

    # Save main mapping file
    mapping_data = {
        "dataset": dataset_name,
        "created_at": datetime.now().isoformat(),
        "spatialcore_version": spatialcore_version,
        "source_column": source_col,
        "summary": {
            "total_labels": len(mappings),
            "mapped": len(mapped_entries),
            "unmapped": len(missed_entries),
            "match_rate": round(len(mapped_entries) / len(mappings), 3) if mappings else 0,
        },
        "tier_breakdown": tier_counts,
        "mappings": mapped_entries,
    }

    mapping_path = save_dir / f"{base_name}_mapping.json"
    with open(mapping_path, "w") as f:
        json.dump(mapping_data, f, indent=2)
    logger.info(f"Saved mapping results to: {mapping_path}")

    # Save missed terms file (if any)
    missed_path = None
    if missed_entries:
        missed_data = {
            "dataset": dataset_name,
            "created_at": datetime.now().isoformat(),
            "source_column": source_col,
            "missed_terms": missed_entries,
            "recommendations": [
                f"Review '{m['input_label']}' ({m['n_cells']} cells) - may need manual mapping"
                for m in missed_entries[:5]
            ],
        }
        missed_path = save_dir / f"{base_name}_missed.json"
        with open(missed_path, "w") as f:
            json.dump(missed_data, f, indent=2)
        logger.info(f"Saved missed terms to: {missed_path}")

    return mapping_path, missed_path


# ============================================================================
# Mapping Table and Metadata Generation
# ============================================================================


def create_mapping_table(
    mappings: Dict[str, List[Dict]],
    cell_counts: Dict[str, int],
    skipped_labels: Optional[List[str]] = None,
    index_source: Optional[str] = None,
    min_score: float = 0.7,
    dataset_name: str = "ontology_mapping",
) -> OntologyMappingResult:
    """
    Create a structured mapping table from search results.

    Produces a DataFrame and metadata suitable for visualization and
    JSON export. This is the primary interface for extracting mapping
    results in a structured format.

    Parameters
    ----------
    mappings : Dict[str, List[Dict]]
        Search results from search_ontology_index().
    cell_counts : Dict[str, int]
        Number of cells per input label.
    skipped_labels : List[str], optional
        Labels that were skipped (e.g., "Unassigned").
    index_source : str, optional
        Path or description of the ontology index used.
    min_score : float, default 0.7
        Minimum score threshold used for matching.
    dataset_name : str, default "ontology_mapping"
        Name for the dataset (used in metadata).

    Returns
    -------
    OntologyMappingResult
        Contains:
        - table: DataFrame with columns [input_label, ontology_name,
          ontology_id, match_tier, score, n_cells, canonical_term]
        - metadata: Full metadata dict for JSON serialization
        - errors: List of mapping errors/warnings

    Examples
    --------
    >>> from spatialcore.annotation import search_ontology_index, create_mapping_table
    >>> labels = ["CD4+ T cells", "B cells", "Unknown_cluster"]
    >>> mappings = search_ontology_index(labels)
    >>> cell_counts = {"CD4+ T cells": 1000, "B cells": 500, "Unknown_cluster": 50}
    >>> result = create_mapping_table(mappings, cell_counts)
    >>> print(result.table)
    >>> result.to_json("mapping_metadata.json")
    """
    try:
        import spatialcore
        spatialcore_version = spatialcore.__version__
    except Exception:
        spatialcore_version = "unknown"

    skipped_labels = skipped_labels or []
    errors = []

    # Build rows for the table
    rows = []
    tier_counts = {
        "tier0_pattern": 0,
        "tier1_exact": 0,
        "tier2_token": 0,
        "tier3_overlap": 0,
        "unmapped": 0,
        "skipped": 0,
    }

    # Process mapped labels
    for label, matches in mappings.items():
        n_cells = cell_counts.get(label, 0)
        canonical = get_canonical_term(label)

        if matches:
            best = matches[0]
            tier = best.get("match_type", "unknown")
            score = best.get("score", 0.0)

            # Normalize tier name for counting
            if "pattern" in tier:
                tier_counts["tier0_pattern"] += 1
            elif "exact" in tier:
                tier_counts["tier1_exact"] += 1
            elif "token" in tier:
                tier_counts["tier2_token"] += 1
            elif "overlap" in tier:
                tier_counts["tier3_overlap"] += 1
            else:
                tier_counts["unmapped"] += 1

            rows.append({
                "input_label": label,
                "ontology_name": best["name"],
                "ontology_id": best["id"],
                "match_tier": tier,
                "score": round(score, 3),
                "n_cells": n_cells,
                "canonical_term": canonical,
            })
        else:
            # Unmapped - use Unknown
            tier_counts["unmapped"] += 1
            rows.append({
                "input_label": label,
                "ontology_name": UNKNOWN_CELL_TYPE_NAME,
                "ontology_id": UNKNOWN_CELL_TYPE_ID,
                "match_tier": "unmapped",
                "score": 0.0,
                "n_cells": n_cells,
                "canonical_term": canonical,
            })

            # Record as error for review
            errors.append({
                "type": "unmapped",
                "label": label,
                "n_cells": n_cells,
                "message": f"No ontology match found for '{label}'",
            })

    # Process skipped labels (Unassigned, Unknown, etc.)
    for label in skipped_labels:
        n_cells = cell_counts.get(label, 0)
        tier_counts["skipped"] += 1
        rows.append({
            "input_label": label,
            "ontology_name": label,  # Keep original
            "ontology_id": "skipped",
            "match_tier": "skipped",
            "score": None,
            "n_cells": n_cells,
            "canonical_term": None,
        })

    # Create DataFrame
    table = pd.DataFrame(rows)

    # Sort by tier (best matches first), then by cell count within tier
    if len(table) > 0:
        # Define tier order (best to worst)
        tier_order = {
            "tier1_exact": 0,
            "tier0_pattern": 1,
            "tier2_token": 2,
            "tier3_overlap": 3,
            "unmapped": 4,
            "skipped": 5,
        }
        table["_tier_order"] = table["match_tier"].map(tier_order).fillna(6)
        table = table.sort_values(
            ["_tier_order", "n_cells"],
            ascending=[True, False]
        ).drop(columns=["_tier_order"]).reset_index(drop=True)

    # Calculate summary statistics
    total_labels = len(mappings) + len(skipped_labels)
    mapped_labels = sum(1 for r in rows if r["match_tier"] not in ["unmapped", "skipped"])
    unmapped_labels = tier_counts["unmapped"]
    skipped_count = tier_counts["skipped"]

    total_cells = sum(cell_counts.values())
    mapped_cells = sum(r["n_cells"] for r in rows if r["match_tier"] not in ["unmapped", "skipped"])
    unmapped_cells = sum(r["n_cells"] for r in rows if r["match_tier"] == "unmapped")

    # Build metadata
    metadata = {
        "dataset_name": dataset_name,
        "created_at": datetime.now().isoformat(),
        "spatialcore_version": spatialcore_version,
        "index_source": str(index_source) if index_source else "package_default",
        "min_score": min_score,
        "summary": {
            "total_labels": total_labels,
            "mapped_labels": mapped_labels,
            "unmapped_labels": unmapped_labels,
            "skipped_labels": skipped_count,
            "mapping_rate": round(mapped_labels / max(total_labels - skipped_count, 1), 3),
            "total_cells": total_cells,
            "mapped_cells": mapped_cells,
            "unmapped_cells": unmapped_cells,
            "cell_mapping_rate": round(mapped_cells / max(total_cells, 1), 3),
        },
        "tier_breakdown": tier_counts,
        "mappings": rows,
        "errors": errors,
    }

    return OntologyMappingResult(
        table=table,
        metadata=metadata,
        errors=errors,
    )


# ============================================================================
# AnnData Integration
# ============================================================================

# Labels that should not be mapped to ontology (placeholders, not cell types)
SKIP_LABELS = {
    "Unassigned", "unassigned", "Unknown", "unknown", "NA", "N/A", "nan",
    "Other", "other", "Doublet", "doublet", "Doublets", "doublets",
    "Low quality", "low quality", "Filtered", "filtered",
}


def has_ontology_ids(
    adata: ad.AnnData,
    id_col: str = "cell_type_ontology_term_id",
    label_col: str = "cell_type",
) -> Dict[str, Any]:
    """
    Check if AnnData has existing ontology IDs and their coverage.

    Use this to decide whether label harmonization is needed before
    subsample_balanced() and training.

    Parameters
    ----------
    adata : AnnData
        Data to check (typically after combine_references()).
    id_col : str, default "cell_type_ontology_term_id"
        Column that may contain existing CL IDs.

    Returns
    -------
    Dict[str, Any]
        Dictionary with keys:
        - has_column: bool - whether id_col exists in adata.obs
        - coverage: float - fraction of cells with valid CL IDs (0.0-1.0)
        - n_with_ids: int - count of cells with valid CL IDs
        - n_without_ids: int - count of cells without valid CL IDs
        - unique_ids: List[str] - unique CL IDs found
        - by_source: Dict[str, float] - coverage by reference_source (if present)

    Examples
    --------
    >>> from spatialcore.annotation import has_ontology_ids
    >>> status = has_ontology_ids(combined)
    >>> print(f"Coverage: {status['coverage']:.1%}")
    Coverage: 65.0%

    >>> # Check per-source coverage
    >>> status['by_source']
    {'cellxgene_lung': 1.0, 'inhouse_batch1': 0.0}

    >>> # Decision logic
    >>> if status['coverage'] < 1.0:
    ...     print("Harmonization recommended")
    """
    result = {
        "has_column": False,
        "coverage": 0.0,
        "n_with_ids": 0,
        "n_without_ids": adata.n_obs,
        "unique_ids": [],
        "by_source": {},
    }

    if id_col not in adata.obs.columns:
        logger.info(f"Column '{id_col}' not found in adata.obs")
        return result

    result["has_column"] = True

    # Check for valid CL IDs (not null and starts with "CL:")
    ids = adata.obs[id_col]
    valid_mask = ids.notna() & ids.astype(str).str.startswith("CL:")

    result["n_with_ids"] = int(valid_mask.sum())
    result["n_without_ids"] = int((~valid_mask).sum())
    result["coverage"] = result["n_with_ids"] / adata.n_obs if adata.n_obs > 0 else 0.0
    result["unique_ids"] = ids[valid_mask].unique().tolist()

    # Calculate per-source coverage if reference_source exists
    if "reference_source" in adata.obs.columns:
        by_source = {}
        for source in adata.obs["reference_source"].unique():
            source_mask = adata.obs["reference_source"] == source
            source_valid = valid_mask & source_mask
            source_total = source_mask.sum()
            by_source[source] = source_valid.sum() / source_total if source_total > 0 else 0.0
        result["by_source"] = by_source

    logger.info(
        f"Ontology ID coverage: {result['coverage']:.1%} "
        f"({result['n_with_ids']:,}/{adata.n_obs:,} cells)"
    )

    return result


def add_ontology_ids(
    adata: ad.AnnData,
    source_col: str,
    target_col: str = "cell_type_ontology_term_id",
    name_col: Optional[str] = "cell_type_ontology_label",
    min_score: float = 0.7,
    index_path: Optional[Union[str, Path]] = None,
    save_mapping: Optional[Union[str, Path]] = None,
    dataset_name: Optional[str] = None,
    skip_labels: Optional[set] = None,
    skip_if_exists: bool = True,
    copy: bool = False,
) -> Tuple[ad.AnnData, Dict[str, List[Dict]], Optional[OntologyMappingResult]]:
    """
    Add ontology IDs to AnnData based on cell type labels.

    Uses CellxGene standard column naming conventions by default.

    Parameters
    ----------
    adata : AnnData
        AnnData object with cell type labels.
    source_col : str
        Column in adata.obs containing cell type labels.
    target_col : str, default "cell_type_ontology_term_id"
        Column to store ontology IDs (CellxGene standard).
    name_col : str, optional, default "cell_type_ontology_label"
        Column to store canonical ontology names. If None, skip.
    min_score : float, default 0.7
        Minimum match score.
    index_path : str or Path, optional
        Path to ontology index JSON file. If None, uses package default.
    save_mapping : str or Path, optional
        Directory to save mapping results. Creates:
        - {dataset}_ontology_mapping.csv: Mapping table
        - {dataset}_ontology_mapping_metadata.json: Full metadata
    dataset_name : str, optional
        Name for output files. If None, uses 'ontology_mapping'.
    skip_labels : set, optional
        Labels to skip (mark as unmapped). If None, uses SKIP_LABELS default
        which includes "Unassigned", "Unknown", "Doublet", etc.
    skip_if_exists : bool, default True
        If True and target_col already exists with valid CL IDs, preserve
        existing IDs and only map cells with missing/invalid IDs. This is
        useful when combining CellxGene data (which has native CL IDs) with
        other references that need mapping.
    copy : bool, default False
        If True, return a copy.

    Returns
    -------
    Tuple[AnnData, Dict, Optional[OntologyMappingResult]]
        - Updated AnnData with new columns (including _tier and _score)
        - Mapping dictionary with all matches per label
        - OntologyMappingResult with table and metadata (if save_mapping provided)

    Notes
    -----
    The function adds the following columns to adata.obs:
    - {target_col}: Ontology ID (e.g., "CL:0000624") or "unknown"/"skipped"
    - {name_col}: Canonical ontology name (if name_col is not None)
    - {target_col}_tier: Match tier (tier0_pattern, tier1_exact, etc.)
    - {target_col}_score: Match score (0.0-1.0)

    Examples
    --------
    >>> from spatialcore.annotation import add_ontology_ids
    >>> adata, mappings, result = add_ontology_ids(
    ...     adata,
    ...     source_col="celltypist",
    ...     save_mapping="./output/",
    ... )
    >>> # View the mapping table
    >>> print(result.table)
    >>> # Save metadata JSON
    >>> result.to_json("ontology_mapping_metadata.json")
    """
    if copy:
        adata = adata.copy()

    if source_col not in adata.obs.columns:
        raise ValueError(
            f"Source column '{source_col}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    # Use default skip labels if not provided
    labels_to_skip = skip_labels if skip_labels is not None else SKIP_LABELS

    # Check for existing valid CL IDs
    labels_with_existing_ids = set()
    existing_id_map = {}  # label -> existing CL ID
    existing_name_map = {}  # label -> existing ontology name

    if skip_if_exists and target_col in adata.obs.columns:
        # Load ontology index for name lookups
        ontology_index = load_ontology_index(index_path)
        cl_index = ontology_index.get("cl", {})

        # Build ID -> name lookup
        id_to_name = {}
        for term_lower, term_data in cl_index.items():
            id_to_name[term_data["id"]] = term_data["name"]

        # Check which labels already have valid CL IDs
        for label in adata.obs[source_col].dropna().unique():
            label_mask = adata.obs[source_col] == label
            existing_ids = adata.obs.loc[label_mask, target_col].dropna()

            # Find valid CL IDs
            valid_ids = [
                eid for eid in existing_ids.unique()
                if isinstance(eid, str) and eid.startswith("CL:")
            ]

            if valid_ids:
                # Use the most common valid ID
                id_counts = existing_ids.value_counts()
                for top_id in id_counts.index:
                    if isinstance(top_id, str) and top_id.startswith("CL:"):
                        existing_id_map[label] = top_id
                        existing_name_map[label] = id_to_name.get(top_id, str(label))
                        labels_with_existing_ids.add(label)
                        break

        if labels_with_existing_ids:
            logger.info(
                f"Preserving existing CL IDs for {len(labels_with_existing_ids)} labels "
                f"(skip_if_exists=True)"
            )

    # Get unique labels
    unique_labels = adata.obs[source_col].dropna().unique().tolist()

    # Filter out skip labels and labels with existing IDs before searching
    labels_to_map = [
        l for l in unique_labels
        if l not in labels_to_skip and l not in labels_with_existing_ids
    ]
    skipped = [l for l in unique_labels if l in labels_to_skip]

    if skipped:
        logger.info(f"Skipping {len(skipped)} non-cell-type labels: {skipped}")

    logger.info(f"Mapping {len(labels_to_map)} unique cell types to ontology...")

    # Search for matches (only labels_to_map, not skipped ones)
    mappings = search_ontology_index(
        labels_to_map,
        index_path=index_path,
        annotation_type="cell_type",  # CL-only
        min_score=min_score,
    )

    # Create label → ID, name, tier, and score mappings
    label_to_id = {}
    label_to_name = {}
    label_to_tier = {}
    label_to_score = {}

    n_matched = 0
    n_unmatched = 0
    n_preserved = 0

    # Handle labels with existing IDs (preserved from input data)
    for label in labels_with_existing_ids:
        label_to_id[label] = existing_id_map[label]
        label_to_name[label] = existing_name_map[label]
        label_to_tier[label] = "existing"
        label_to_score[label] = 1.0  # Full confidence for existing IDs
        n_preserved += 1

    # Handle skipped labels (mark as skipped)
    for label in skipped:
        label_to_id[label] = "skipped"
        label_to_name[label] = label  # Keep original
        label_to_tier[label] = "skipped"
        label_to_score[label] = None

    # Process search results
    for label, matches in mappings.items():
        if matches:
            best_match = matches[0]
            label_to_id[label] = best_match["id"]
            label_to_name[label] = best_match["name"]
            label_to_tier[label] = best_match.get("match_type", "unknown")
            label_to_score[label] = round(best_match.get("score", 0.0), 3)
            n_matched += 1
        else:
            label_to_id[label] = UNKNOWN_CELL_TYPE_ID
            label_to_name[label] = UNKNOWN_CELL_TYPE_NAME
            label_to_tier[label] = "unmapped"
            label_to_score[label] = 0.0
            n_unmatched += 1

    total_to_map = len(labels_to_map)
    if total_to_map > 0:
        logger.info(f"  Matched: {n_matched}/{total_to_map} ({100*n_matched/total_to_map:.1f}%)")
    else:
        logger.info("  No labels to map")
    if n_preserved > 0:
        logger.info(f"  Preserved (existing CL IDs): {n_preserved}")
    if n_unmatched > 0:
        unmatched = [l for l, m in mappings.items() if not m]
        logger.warning(f"  Unmatched labels: {unmatched[:5]}{'...' if len(unmatched) > 5 else ''}")

    # Apply mappings to adata
    adata.obs[target_col] = adata.obs[source_col].map(label_to_id)
    if name_col:
        adata.obs[name_col] = adata.obs[source_col].map(label_to_name)

    # Store tier and score information (derive column names from target_col)
    # Handle both old format (cell_type_ontology_id) and new CellxGene format (cell_type_ontology_term_id)
    if "_term_id" in target_col:
        tier_col = target_col.replace("_term_id", "_tier")
        score_col = target_col.replace("_term_id", "_score")
    else:
        tier_col = target_col.replace("_id", "_tier")
        score_col = target_col.replace("_id", "_score")
    adata.obs[tier_col] = adata.obs[source_col].map(label_to_tier)
    adata.obs[score_col] = adata.obs[source_col].map(label_to_score)

    # Create mapping result with table and metadata
    mapping_result = None
    if save_mapping:
        save_dir = Path(save_mapping)
        save_dir.mkdir(parents=True, exist_ok=True)

        # Count cells per label
        cell_counts = adata.obs[source_col].value_counts().to_dict()

        # Create structured mapping result
        mapping_result = create_mapping_table(
            mappings=mappings,
            cell_counts=cell_counts,
            skipped_labels=list(skipped),
            index_source=str(index_path) if index_path else None,
            min_score=min_score,
            dataset_name=dataset_name or "ontology_mapping",
        )

        # Save artifacts
        name = dataset_name or "ontology_mapping"
        csv_path = save_dir / f"{name}_ontology_mapping.csv"
        json_path = save_dir / f"{name}_ontology_mapping_metadata.json"

        mapping_result.to_csv(csv_path)
        mapping_result.to_json(json_path)

        logger.info(f"Saved mapping table to: {csv_path}")
        logger.info(f"Saved mapping metadata to: {json_path}")

    return adata, mappings, mapping_result


def validate_cl_term(term_id: str, ontology_index: Optional[Dict] = None) -> bool:
    """
    Check if a CL term ID exists in the ontology.

    Parameters
    ----------
    term_id : str
        Cell Ontology ID (e.g., "CL:0000624").
    ontology_index : Dict, optional
        Pre-loaded ontology index.

    Returns
    -------
    bool
        True if term exists in CL ontology.
    """
    if ontology_index is None:
        ontology_index = load_ontology_index()

    cl_index = ontology_index.get("cl", {})
    for term_data in cl_index.values():
        if term_data.get("id") == term_id:
            return True
    return False


def get_cl_id(term_name: str, ontology_index: Optional[Dict] = None) -> Optional[str]:
    """
    Get CL ID for a term name (exact match).

    Parameters
    ----------
    term_name : str
        Term name to look up.
    ontology_index : Dict, optional
        Pre-loaded ontology index.

    Returns
    -------
    str or None
        CL ID if found, None otherwise.
    """
    if ontology_index is None:
        ontology_index = load_ontology_index()

    term_lower = term_name.lower().strip()
    cl_index = ontology_index.get("cl", {})

    if term_lower in cl_index:
        return cl_index[term_lower]["id"]
    return None
