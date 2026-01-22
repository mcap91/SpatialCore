"""
CellxGene reference data download utilities.

This module provides utilities for downloading reference datasets from
CellxGene Census, including:
- Downloading predefined datasets by key
- Querying Census with flexible filters (tissue, disease, cell type)
- Listing available datasets

Gene mapping utilities (Ensembl → HUGO) have been moved to spatialcore.core.utils
and are re-exported here for backward compatibility.

References:
    - CellxGene Census: https://chanzuckerberg.github.io/cellxgene-census/
"""

from pathlib import Path
from typing import Dict, List, Optional, Any, Union, Set
import re

import numpy as np
import pandas as pd
import anndata as ad

from spatialcore.core.logging import get_logger

# Re-export gene mapping utilities from core/utils for backward compatibility
from spatialcore.core.utils import (
    load_ensembl_to_hugo_mapping,
    normalize_gene_names,
    check_normalization_status,
    download_ensembl_mapping,
    is_ensembl_id,
)

logger = get_logger(__name__)

_LABEL_TOKEN_PATTERN = re.compile(r"[^a-z0-9]+")


def _label_tokens(label: str) -> Set[str]:
    """Normalize a label into lowercase alphanumeric tokens."""
    if label is None:
        return set()
    normalized = _LABEL_TOKEN_PATTERN.sub(" ", str(label).lower()).strip()
    if not normalized:
        return set()
    return set(token for token in normalized.split() if token)


def _detect_parent_child_conflicts(
    adata: ad.AnnData,
    label_column: str,
    ontology_column: str,
    min_parent_tokens: int = 2,
) -> Dict[str, List[str]]:
    """
    Detect parent/child conflicts based on label token containment.

    Returns dict mapping parent CL IDs to list of child CL IDs present in data.
    """
    pairs = adata.obs[[label_column, ontology_column]].dropna()
    if pairs.empty:
        return {}

    pairs = pairs.astype(str)
    valid_mask = pairs[ontology_column].str.startswith("CL:")
    pairs = pairs.loc[valid_mask, [label_column, ontology_column]]
    if pairs.empty:
        return {}

    # Map CL ID -> most common label for that ID
    id_to_label = (
        pairs.groupby(ontology_column)[label_column]
        .agg(lambda values: values.value_counts().idxmax())
        .to_dict()
    )
    id_to_tokens = {cl_id: _label_tokens(label) for cl_id, label in id_to_label.items()}

    conflicts: Dict[str, List[str]] = {}
    for parent_id, parent_tokens in id_to_tokens.items():
        if len(parent_tokens) < min_parent_tokens:
            continue
        for child_id, child_tokens in id_to_tokens.items():
            if parent_id == child_id:
                continue
            if parent_tokens == child_tokens:
                continue
            if len(child_tokens) <= len(parent_tokens):
                continue
            if parent_tokens.issubset(child_tokens):
                conflicts.setdefault(parent_id, []).append(child_id)

    return conflicts

# ============================================================================
# CellxGene Dataset Registry
# ============================================================================

CELLXGENE_DATASETS: Dict[str, Dict[str, Any]] = {
    # Liver datasets
    "healthy_human_liver": {
        "dataset_id": "4f88c1be-5156-463d-b64d-a3a3a8e0da6d",
        "description": "Cell types from scRNA-seq and snRNA-seq of healthy human liver",
        "tissue": "liver",
        "cell_type_column": "cell_type",
        "expected_cells": "~100,000",
    },
    # Colon / GI datasets
    "colon_immune_niches": {
        "dataset_id": "2872f4b0-b171-46e2-abc6-befcf6de6306",
        "description": "Distinct microbial and immune niches of the human colon",
        "tissue": "colon",
        "cell_type_column": "cell_type",
        "expected_cells": "~41,650",
    },
    "colon_ulcerative_colitis": {
        "dataset_id": "4dd00779-7f73-4f50-89bb-e2d3c6b71b18",
        "description": "Human Colon during Ulcerative Colitis (Smillie et al.)",
        "tissue": "colon",
        "cell_type_column": "cell_type",
        "expected_cells": "~34,772",
    },
    "colon_crohns_immune": {
        "dataset_id": "518d9049-2a76-44f8-8abc-1e2b59ab5ba1",
        "description": "Crohn's disease colon immune cells",
        "tissue": "colon",
        "cell_type_column": "cell_type",
        "expected_cells": "~152,509",
    },
    # Lung datasets
    "human_lung_cell_atlas": {
        "dataset_id": "f72958f5-7f42-4ebb-98da-445b0c6de516",
        "description": "Human Lung Cell Atlas (HLCA) - Azimuth",
        "tissue": "lung",
        "cell_type_column": "ann_finest_level",
        "expected_cells": "~584,884",
    },
    "lung_covid": {
        "dataset_id": "d8da613f-e681-4c69-b463-e94f5e66847f",
        "description": "Molecular single-cell lung atlas of lethal COVID-19",
        "tissue": "lung",
        "cell_type_column": "cell_type",
        "expected_cells": "~116,313",
    },
    # CRC datasets
    "crc_htan_epithelial_discovery": {
        "dataset_id": "e40c6272-af77-4a10-9385-62a398884f27",
        "description": "HTAN VUMC CRC Polyps - Epithelial (Discovery)",
        "tissue": "colon",
        "cell_type_column": "cell_type",
        "expected_cells": "~65,088",
    },
}


def list_available_datasets() -> pd.DataFrame:
    """
    List all available CellxGene datasets with metadata.

    Returns
    -------
    pd.DataFrame
        DataFrame with dataset keys, descriptions, tissues, and expected cell counts.
    """
    records = []
    for key, info in CELLXGENE_DATASETS.items():
        records.append({
            "dataset_key": key,
            "description": info["description"],
            "tissue": info["tissue"],
            "cell_type_column": info["cell_type_column"],
            "expected_cells": info.get("expected_cells", "unknown"),
        })
    return pd.DataFrame(records)


def download_cellxgene_reference(
    dataset_key: str,
    output_dir: Union[str, Path],
    force: bool = False,
) -> Path:
    """
    Download a reference dataset from CellxGene Census.

    Parameters
    ----------
    dataset_key : str
        Key from CELLXGENE_DATASETS registry (e.g., "healthy_human_liver").
    output_dir : str or Path
        Directory to save the downloaded h5ad file.
    force : bool, default False
        If True, re-download even if file exists.

    Returns
    -------
    Path
        Path to the downloaded h5ad file.

    Raises
    ------
    ValueError
        If dataset_key is not in CELLXGENE_DATASETS.
    ImportError
        If cellxgene-census is not installed.

    Examples
    --------
    >>> from spatialcore.annotation import download_cellxgene_reference
    >>> path = download_cellxgene_reference("healthy_human_liver", "./references")
    >>> print(path)
    references/healthy_human_liver.h5ad
    """
    if dataset_key not in CELLXGENE_DATASETS:
        available = ", ".join(CELLXGENE_DATASETS.keys())
        raise ValueError(
            f"Unknown dataset: '{dataset_key}'. Available: {available}"
        )

    try:
        import cellxgene_census
    except ImportError:
        raise ImportError(
            "cellxgene-census is required for downloading CellxGene data. "
            "Install with: pip install cellxgene-census"
        )

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / f"{dataset_key}.h5ad"

    if output_file.exists() and not force:
        logger.info(f"Dataset already exists: {output_file}")
        return output_file

    dataset_info = CELLXGENE_DATASETS[dataset_key]
    dataset_id = dataset_info["dataset_id"]

    logger.info(f"Downloading {dataset_key} (ID: {dataset_id})...")
    logger.info(f"  Description: {dataset_info['description']}")
    logger.info(f"  Expected cells: {dataset_info.get('expected_cells', 'unknown')}")

    # Download using Census API
    cellxgene_census.download_source_h5ad(
        dataset_id,
        to_path=str(output_file),
    )

    logger.info(f"Downloaded to: {output_file}")
    return output_file


def query_cellxgene_census(
    tissue: Optional[str] = None,
    disease: Optional[str] = None,
    cell_type: Optional[str] = None,
    assay: Optional[str] = None,
    organism: str = "Homo sapiens",
    obs_columns: Optional[List[str]] = None,
    max_cells: Optional[int] = None,
    output_path: Optional[Union[str, Path]] = None,
    random_state: int = 42,
    validate_labels: bool = True,
    resolve_hierarchy: str = "none",
) -> ad.AnnData:
    """
    Query cells from CellxGene Census with flexible filters.

    This provides more flexibility than download_cellxgene_reference() by
    allowing arbitrary tissue/disease/cell_type combinations.

    Parameters
    ----------
    tissue : str, optional
        Tissue filter (e.g., "liver", "lung", "colon").
    disease : str, optional
        Disease filter (e.g., "normal", "hepatocellular carcinoma").
    cell_type : str, optional
        Cell type filter (e.g., "T cell", "hepatocyte").
    assay : str, optional
        Assay filter (e.g., "10x 3' v3", "Smart-seq2").
    organism : str, default "Homo sapiens"
        Organism to query.
    obs_columns : List[str], optional
        Columns to include in obs. Default: cell_type, disease, assay, tissue.
    max_cells : int, optional
        Maximum cells to return. Default None downloads ALL matching cells
        (recommended for production). If specified, uses memory-efficient
        sampling: queries cell IDs first, samples in memory, then downloads
        only the sampled cells. Use this for testing/development to avoid
        OOM errors on memory-constrained systems.
    output_path : str or Path, optional
        If provided, save result to this h5ad file.
    random_state : int, default 42
        Random seed for subsampling (only used when max_cells is specified).
    validate_labels : bool, default True
        If True, check for label-to-ontology inconsistencies in CellxGene
        columns (cell_type vs cell_type_ontology_term_id) and log warnings.
    resolve_hierarchy : str, default "none"
        If "remove_parents", drop cells labeled with parent terms when any
        child terms are present (based on label token containment). Use "none"
        to keep current behavior.

    Returns
    -------
    AnnData
        AnnData object with queried cells.

    Raises
    ------
    ImportError
        If cellxgene-census is not installed (Linux only, no Windows support).
    ValueError
        If no filter criteria provided.

    Examples
    --------
    >>> from spatialcore.annotation import query_cellxgene_census
    >>> # Production: Download ALL healthy liver cells
    >>> adata = query_cellxgene_census(
    ...     tissue="liver",
    ...     disease="normal",
    ...     output_path="./references/healthy_liver.h5ad"
    ... )
    >>> # Testing: Sample 5000 cells (memory-efficient for development)
    >>> sample = query_cellxgene_census(
    ...     tissue="liver",
    ...     disease="hepatocellular carcinoma",
    ...     max_cells=5000,  # Only for testing
    ... )
    """
    try:
        import cellxgene_census
    except ImportError:
        raise ImportError(
            "cellxgene-census is required for querying CellxGene data. "
            "Install with: pip install cellxgene-census"
        )

    # Build filter string
    filters = ["is_primary_data == True"]
    if tissue:
        filters.append(f"tissue == '{tissue}'")
    if disease:
        filters.append(f"disease == '{disease}'")
    if cell_type:
        filters.append(f"cell_type == '{cell_type}'")
    if assay:
        filters.append(f"assay == '{assay}'")

    if len(filters) == 1:
        raise ValueError(
            "At least one filter (tissue, disease, cell_type, or assay) is required"
        )

    filter_string = " and ".join(filters)

    # Default obs columns - includes ontology ID if available in Census
    if obs_columns is None:
        obs_columns = [
            "cell_type",
            "cell_type_ontology_term_id",  # CL ID from CellxGene curators
            "disease",
            "assay",
            "dataset_id",
            "tissue",
        ]

    logger.info("Querying CellxGene Census...")
    logger.info(f"  Organism: {organism}")
    logger.info(f"  Filter: {filter_string}")

    with cellxgene_census.open_soma() as census:
        # Memory-efficient approach: sample cell IDs BEFORE downloading expression data
        # This prevents OOM by only fetching the cells we actually need

        # Convert organism name to Census key format (e.g., "Homo sapiens" -> "homo_sapiens")
        organism_key = organism.lower().replace(" ", "_")

        if max_cells:
            # Step 1: Get cell IDs matching filter (lightweight - no expression data)
            logger.info("  Step 1: Counting cells matching filter...")
            human = census["census_data"][organism_key]
            obs_df = human.obs.read(
                value_filter=filter_string,
                column_names=["soma_joinid"],  # Only get IDs, very lightweight
            ).concat().to_pandas()

            total_cells = len(obs_df)
            logger.info(f"  Found {total_cells:,} cells matching filter")

            # Step 2: Sample cell IDs if needed
            if total_cells > max_cells:
                logger.info(f"  Step 2: Sampling {max_cells:,} cell IDs (memory-efficient)...")
                np.random.seed(random_state)
                sampled_ids = np.random.choice(
                    obs_df["soma_joinid"].values,
                    size=max_cells,
                    replace=False,
                )
            else:
                sampled_ids = obs_df["soma_joinid"].values
                logger.info(f"  Step 2: Using all {len(sampled_ids):,} cells (under max_cells limit)")

            # Step 3: Download only sampled cells (key memory optimization!)
            logger.info(f"  Step 3: Downloading expression data for {len(sampled_ids):,} cells...")
            adata = cellxgene_census.get_anndata(
                census=census,
                organism=organism,
                obs_coords=sampled_ids,  # Only fetch these specific cells!
                obs_column_names=obs_columns,
            )
        else:
            # No max_cells limit - download everything (use with caution!)
            logger.warning("  No max_cells limit set - downloading ALL matching cells!")
            logger.warning("  This may use significant memory. Consider setting max_cells.")
            adata = cellxgene_census.get_anndata(
                census=census,
                organism=organism,
                obs_value_filter=filter_string,
                obs_column_names=obs_columns,
            )

    logger.info(f"  Downloaded: {adata.n_obs:,} cells × {adata.n_vars:,} genes")

    if resolve_hierarchy not in {"none", "remove_parents"}:
        raise ValueError("resolve_hierarchy must be 'none' or 'remove_parents'")

    if resolve_hierarchy == "remove_parents":
        if (
            "cell_type" not in adata.obs.columns
            or "cell_type_ontology_term_id" not in adata.obs.columns
        ):
            raise ValueError(
                "resolve_hierarchy='remove_parents' requires "
                "cell_type and cell_type_ontology_term_id in adata.obs"
            )

        conflicts = _detect_parent_child_conflicts(
            adata,
            label_column="cell_type",
            ontology_column="cell_type_ontology_term_id",
        )

        if conflicts:
            parent_ids = set(conflicts.keys())
            parent_mask = adata.obs["cell_type_ontology_term_id"].isin(parent_ids)
            removed = int(parent_mask.sum())
            adata = adata[~parent_mask].copy()
            logger.info(
                "Removed %d parent-labeled cells due to hierarchy conflicts",
                removed,
            )

    if validate_labels:
        if (
            "cell_type" in adata.obs.columns
            and "cell_type_ontology_term_id" in adata.obs.columns
        ):
            from spatialcore.annotation.validation import (
                check_label_ontology_consistency,
            )

            consistency = check_label_ontology_consistency(
                adata,
                label_column="cell_type",
                ontology_column="cell_type_ontology_term_id",
            )

            if consistency.n_labels_with_multiple_ids > 0:
                examples = []
                for label in sorted(consistency.labels_with_multiple_ids.keys())[:5]:
                    ids = ", ".join(consistency.labels_with_multiple_ids[label])
                    examples.append(f"{label} -> {ids}")
                logger.warning(
                    "CellxGene label/ontology mismatch: %d labels map to multiple CL IDs. "
                    "Examples: %s",
                    consistency.n_labels_with_multiple_ids,
                    "; ".join(examples),
                )

            if consistency.n_hierarchical_labels > 0:
                logger.warning(
                    "CellxGene labels look hierarchical (parent/child in one label): %s",
                    ", ".join(sorted(consistency.hierarchical_labels)[:5]),
                )

    # Save if output path provided
    if output_path:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(output_path)
        logger.info(f"  Saved to: {output_path}")

    return adata
