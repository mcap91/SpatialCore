"""
Cell type column validation for CellTypist training.

This module provides comprehensive validation of cell type annotation columns
before training to catch data quality issues early:

1. Column existence
2. Null/empty value detection
3. Cardinality checks (1-500 cell types is reasonable)
4. Minimum cells per type
5. Suspicious pattern detection (e.g., numeric-only labels, placeholder values)
6. Class imbalance detection

References:
    - CellTypist: https://www.celltypist.org/
    - Best practices for cell type annotation in single-cell genomics
"""

from typing import Dict, List, Optional, Any, Literal
from dataclasses import dataclass, field
import re

import pandas as pd
import anndata as ad

from spatialcore.core.logging import get_logger

logger = get_logger(__name__)


# Known placeholder/suspicious patterns
SUSPICIOUS_PATTERNS = [
    (r"^[0-9]+$", "Numeric-only cell type labels"),
    (r"^(unknown|unassigned|na|n/a|none|null)$", "Placeholder values"),
    (r"^cluster_?[0-9]+$", "Cluster-based labels (not biological)"),
    (r"^cell_?[0-9]+$", "Generic cell labels"),
    (r"^type_?[0-9]+$", "Generic type labels"),
    (r"^leiden_?[0-9]*$", "Leiden clustering labels"),
    (r"^louvain_?[0-9]*$", "Louvain clustering labels"),
]


@dataclass
class ValidationIssue:
    """A single validation issue.

    Attributes
    ----------
    severity : str
        One of "error", "warning", or "info".
    code : str
        Machine-readable code (e.g., "NULL_VALUES", "LOW_CELL_COUNTS").
    message : str
        Human-readable description of the issue.
    details : dict, optional
        Additional details about the issue.
    """

    severity: Literal["error", "warning", "info"]
    code: str
    message: str
    details: Optional[Dict[str, Any]] = None

    def __repr__(self) -> str:
        return f"[{self.severity.upper()}] {self.code}: {self.message}"


@dataclass
class CellTypeValidationResult:
    """Complete validation result for a cell type column.

    Attributes
    ----------
    column_name : str
        Name of the validated column.
    is_valid : bool
        True if no errors were found (warnings are OK).
    n_cells : int
        Total number of cells.
    n_cell_types : int
        Number of unique cell types.
    issues : list
        List of ValidationIssue objects.
    cell_type_counts : pd.Series, optional
        Counts per cell type, sorted descending.
    """

    column_name: str
    is_valid: bool
    n_cells: int
    n_cell_types: int
    issues: List[ValidationIssue] = field(default_factory=list)
    cell_type_counts: Optional[pd.Series] = None

    @property
    def errors(self) -> List[ValidationIssue]:
        """Get all error-level issues."""
        return [i for i in self.issues if i.severity == "error"]

    @property
    def warnings(self) -> List[ValidationIssue]:
        """Get all warning-level issues."""
        return [i for i in self.issues if i.severity == "warning"]

    @property
    def infos(self) -> List[ValidationIssue]:
        """Get all info-level issues."""
        return [i for i in self.issues if i.severity == "info"]

    def summary(self) -> str:
        """Return human-readable summary of validation results."""
        status = "PASSED" if self.is_valid else "FAILED"
        lines = [
            f"Validation: {status}",
            f"  Column: {self.column_name}",
            f"  Cells: {self.n_cells:,}",
            f"  Cell types: {self.n_cell_types}",
        ]
        if self.errors:
            lines.append(f"  Errors: {len(self.errors)}")
            for e in self.errors:
                lines.append(f"    - {e.message}")
        if self.warnings:
            lines.append(f"  Warnings: {len(self.warnings)}")
            for w in self.warnings:
                lines.append(f"    - {w.message}")
        return "\n".join(lines)


def validate_cell_type_column(
    adata: ad.AnnData,
    column: str,
    min_cells_per_type: int = 10,
    max_cell_types: int = 500,
    min_cell_types: int = 1,
    allow_nulls: bool = False,
    max_null_fraction: float = 0.05,
    check_suspicious_patterns: bool = True,
) -> CellTypeValidationResult:
    """
    Validate a cell type annotation column for CellTypist training.

    Performs comprehensive checks:
    1. Column existence
    2. Null/empty values (error if >5% by default)
    3. Cardinality (1-500 cell types)
    4. Minimum cells per type (10 by default)
    5. Suspicious patterns (numeric-only, placeholders)
    6. Class imbalance (warn if >1000x ratio)

    Parameters
    ----------
    adata : AnnData
        AnnData object to validate.
    column : str
        Name of cell type column in adata.obs.
    min_cells_per_type : int, default 10
        Minimum cells required per cell type. Types with fewer cells
        will generate a warning.
    max_cell_types : int, default 500
        Maximum number of cell types (warn if exceeded).
    min_cell_types : int, default 1
        Minimum number of cell types required (error if not met).
    allow_nulls : bool, default False
        If True, allow null values (still warns if >max_null_fraction).
    max_null_fraction : float, default 0.05
        Maximum fraction of null values before error (if allow_nulls=False).
    check_suspicious_patterns : bool, default True
        Check for suspicious label patterns like numeric-only or cluster labels.

    Returns
    -------
    CellTypeValidationResult
        Complete validation result with issues and statistics.

    Examples
    --------
    >>> from spatialcore.annotation import validate_cell_type_column
    >>> result = validate_cell_type_column(adata, "cell_type")
    >>> if not result.is_valid:
    ...     print(result.summary())
    ...     raise ValueError("Validation failed")
    >>> print(f"Found {result.n_cell_types} cell types")

    >>> # Inspect low-count types
    >>> low_count = result.cell_type_counts[result.cell_type_counts < 100]
    >>> print(low_count)
    """
    issues: List[ValidationIssue] = []

    # Check 1: Column exists
    if column not in adata.obs.columns:
        available = list(adata.obs.columns)
        return CellTypeValidationResult(
            column_name=column,
            is_valid=False,
            n_cells=adata.n_obs,
            n_cell_types=0,
            issues=[
                ValidationIssue(
                    severity="error",
                    code="COLUMN_NOT_FOUND",
                    message=f"Column '{column}' not found in adata.obs",
                    details={"available_columns": available[:20]},
                )
            ],
        )

    labels = adata.obs[column]
    n_cells = len(labels)

    # Check 2: Null values
    null_mask = labels.isna() | (labels.astype(str).str.strip() == "")
    n_nulls = int(null_mask.sum())
    null_fraction = n_nulls / n_cells if n_cells > 0 else 0.0

    if n_nulls > 0:
        if null_fraction > max_null_fraction and not allow_nulls:
            issues.append(
                ValidationIssue(
                    severity="error",
                    code="EXCESSIVE_NULLS",
                    message=f"{n_nulls:,} null values ({null_fraction:.1%} of cells)",
                    details={"n_nulls": n_nulls, "fraction": float(null_fraction)},
                )
            )
        else:
            issues.append(
                ValidationIssue(
                    severity="warning",
                    code="NULL_VALUES",
                    message=f"{n_nulls:,} null values ({null_fraction:.1%} of cells)",
                    details={"n_nulls": n_nulls, "fraction": float(null_fraction)},
                )
            )

    # Filter out nulls for remaining checks
    valid_labels = labels[~null_mask].astype(str)
    cell_type_counts = valid_labels.value_counts()
    n_cell_types = len(cell_type_counts)

    # Check 3: Cardinality - too few
    if n_cell_types < min_cell_types:
        issues.append(
            ValidationIssue(
                severity="error",
                code="TOO_FEW_TYPES",
                message=f"Only {n_cell_types} cell types found (minimum: {min_cell_types})",
                details={"n_cell_types": n_cell_types, "min_required": min_cell_types},
            )
        )

    # Check 4: Cardinality - too many
    if n_cell_types > max_cell_types:
        issues.append(
            ValidationIssue(
                severity="warning",
                code="MANY_CELL_TYPES",
                message=f"{n_cell_types} cell types found (may indicate over-annotation)",
                details={"n_cell_types": n_cell_types, "max_typical": max_cell_types},
            )
        )

    # Check 5: Minimum cells per type
    low_count_types = cell_type_counts[cell_type_counts < min_cells_per_type]
    if len(low_count_types) > 0:
        issues.append(
            ValidationIssue(
                severity="warning",
                code="LOW_CELL_COUNTS",
                message=f"{len(low_count_types)} cell types have <{min_cells_per_type} cells",
                details={
                    "affected_types": dict(low_count_types.head(10)),
                    "n_affected": len(low_count_types),
                },
            )
        )

    # Check 6: Suspicious patterns
    if check_suspicious_patterns and n_cell_types > 0:
        suspicious_types = []
        for cell_type in cell_type_counts.index:
            for pattern, description in SUSPICIOUS_PATTERNS:
                if re.match(pattern, str(cell_type), re.IGNORECASE):
                    suspicious_types.append((str(cell_type), description))
                    break

        if suspicious_types:
            issues.append(
                ValidationIssue(
                    severity="warning",
                    code="SUSPICIOUS_LABELS",
                    message=f"{len(suspicious_types)} cell types have suspicious labels",
                    details={
                        "suspicious_types": [
                            {"label": t, "reason": r} for t, r in suspicious_types[:10]
                        ],
                        "n_total": len(suspicious_types),
                    },
                )
            )

    # Check 7: Highly imbalanced data
    if n_cell_types >= 2:
        max_count = int(cell_type_counts.max())
        min_count = int(cell_type_counts.min())
        imbalance_ratio = max_count / max(min_count, 1)

        if imbalance_ratio > 1000:
            issues.append(
                ValidationIssue(
                    severity="warning",
                    code="HIGHLY_IMBALANCED",
                    message=f"Extreme class imbalance: {imbalance_ratio:.0f}x ratio",
                    details={
                        "largest_type": str(cell_type_counts.idxmax()),
                        "largest_count": max_count,
                        "smallest_type": str(cell_type_counts.idxmin()),
                        "smallest_count": min_count,
                        "imbalance_ratio": float(imbalance_ratio),
                    },
                )
            )

    # Determine overall validity (only errors make it invalid)
    is_valid = not any(i.severity == "error" for i in issues)

    result = CellTypeValidationResult(
        column_name=column,
        is_valid=is_valid,
        n_cells=n_cells,
        n_cell_types=n_cell_types,
        issues=issues,
        cell_type_counts=cell_type_counts,
    )

    # Log summary
    if is_valid:
        logger.info(
            f"Validation passed: {column} ({n_cell_types} cell types, {n_cells:,} cells)"
        )
        for issue in issues:
            logger.warning(f"  {issue.code}: {issue.message}")
    else:
        logger.error(f"Validation failed: {column}")
        for issue in issues:
            log_fn = logger.error if issue.severity == "error" else logger.warning
            log_fn(f"  {issue.code}: {issue.message}")

        return result

    return result


@dataclass
class LabelOntologyConsistencyResult:
    """Result of checking label to ontology ID consistency."""

    label_column: str
    ontology_column: str
    n_labels: int
    n_labels_with_multiple_ids: int
    labels_with_multiple_ids: Dict[str, List[str]]
    n_hierarchical_labels: int
    hierarchical_labels: List[str]


_HIERARCHY_PATTERN = re.compile(r"(?:\s>\s|\s->\s|;|\|)")


def check_label_ontology_consistency(
    adata: ad.AnnData,
    label_column: str,
    ontology_column: str,
    detect_hierarchy: bool = True,
) -> LabelOntologyConsistencyResult:
    """
    Check whether each label maps to a single ontology ID.

    Flags labels that map to multiple valid CL IDs, which can cause label
    collapsing when IDs are inferred from labels.
    """
    if label_column not in adata.obs.columns:
        raise ValueError(
            f"Label column '{label_column}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )
    if ontology_column not in adata.obs.columns:
        raise ValueError(
            f"Ontology column '{ontology_column}' not found in adata.obs. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    labels = adata.obs[label_column].dropna().astype(str)
    n_labels = int(labels.nunique())

    pairs = adata.obs[[label_column, ontology_column]].dropna()
    pairs = pairs.drop_duplicates().astype(str)
    valid_mask = pairs[ontology_column].str.startswith("CL:")
    unique_pairs = pairs.loc[valid_mask, [label_column, ontology_column]]

    labels_with_multiple_ids: Dict[str, List[str]] = {}
    if not unique_pairs.empty:
        grouped = unique_pairs.groupby(label_column)[ontology_column].unique()
        for label, ids in grouped.items():
            unique_ids = sorted(set(ids))
            if len(unique_ids) > 1:
                labels_with_multiple_ids[str(label)] = unique_ids

    hierarchical_labels: List[str] = []
    if detect_hierarchy:
        for label in labels.unique():
            if _HIERARCHY_PATTERN.search(str(label)):
                hierarchical_labels.append(str(label))

    return LabelOntologyConsistencyResult(
        label_column=label_column,
        ontology_column=ontology_column,
        n_labels=n_labels,
        n_labels_with_multiple_ids=len(labels_with_multiple_ids),
        labels_with_multiple_ids=labels_with_multiple_ids,
        n_hierarchical_labels=len(hierarchical_labels),
        hierarchical_labels=hierarchical_labels,
    )


def validate_multiple_columns(
    adatas: List[ad.AnnData],
    columns: List[str],
    raise_on_error: bool = True,
    **kwargs,
) -> List[CellTypeValidationResult]:
    """
    Validate cell type columns across multiple AnnData objects.

    Convenience function for validating multiple references before combining.

    Parameters
    ----------
    adatas : List[AnnData]
        List of AnnData objects to validate.
    columns : List[str]
        Cell type column name for each AnnData. Must have same length as adatas.
    raise_on_error : bool, default True
        If True, raise ValueError when any validation fails with errors.
    **kwargs
        Additional arguments passed to validate_cell_type_column.

    Returns
    -------
    List[CellTypeValidationResult]
        Validation results for each AnnData.

    Raises
    ------
    ValueError
        If lengths of adatas and columns don't match.
        If any validation fails with errors and raise_on_error=True.

    Examples
    --------
    >>> from spatialcore.annotation import validate_multiple_columns
    >>> results = validate_multiple_columns(
    ...     [adata1, adata2],
    ...     ["cell_type", "annotation"],
    ... )
    >>> for r in results:
    ...     print(r.summary())
    """
    if len(adatas) != len(columns):
        raise ValueError(
            f"Number of adatas ({len(adatas)}) must match columns ({len(columns)})"
        )

    results = []
    all_valid = True

    for i, (adata, col) in enumerate(zip(adatas, columns)):
        logger.info(f"Validating reference {i + 1}/{len(adatas)}: column '{col}'")
        result = validate_cell_type_column(adata, col, **kwargs)
        results.append(result)
        if not result.is_valid:
            all_valid = False

    if not all_valid and raise_on_error:
        failed = [r for r in results if not r.is_valid]
        error_msgs = []
        for r in failed:
            for e in r.errors:
                error_msgs.append(f"{r.column_name}: {e.message}")
        raise ValueError(
            f"Validation failed for {len(failed)} reference(s):\n"
            + "\n".join(error_msgs)
        )

    return results
