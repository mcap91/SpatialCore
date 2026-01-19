"""Boolean expression evaluation for ontology IDs.

This module provides functionality to evaluate boolean expressions
on ontology IDs (CL, NCIT, UBERON) for filtering cells in spatial data.

Examples
--------
>>> from spatialcore.annotation.expression import evaluate_ontology_expression
>>> import anndata as ad
>>> adata = ad.read_h5ad("annotated.h5ad")
>>> # Single ontology filter
>>> mask = evaluate_ontology_expression("CL:0000236", adata)  # B cells
>>> # AND expression
>>> mask = evaluate_ontology_expression("CL:0000236 & NCIT:C4349", adata)  # B cells AND tumor
>>> # OR expression
>>> mask = evaluate_ontology_expression("CL:0000236 | CL:0000624", adata)  # B cells OR CD4+ T cells
>>> # NOT expression
>>> mask = evaluate_ontology_expression("~NCIT:C4349", adata)  # NOT tumor
>>> # Complex expression with parentheses
>>> mask = evaluate_ontology_expression("(CL:0000236 | CL:0000624) & ~NCIT:C4349", adata)
"""

import re
from typing import List, Optional, Set

import numpy as np
import anndata as ad

from spatialcore.core.logging import get_logger

logger = get_logger(__name__)

# Default columns to search for ontology IDs
DEFAULT_ONTOLOGY_COLUMNS = [
    "cell_type_ontology_id",
    "disease_ontology_id",
    "tissue_ontology_id",
    "cell_type_ontology_term_id",  # Alternative naming
    "celltype_ontology_id",  # Alternative naming
]

# Regex pattern for ontology IDs (e.g., CL:0000236, NCIT:C4349, UBERON:0002107)
ONTOLOGY_ID_PATTERN = re.compile(r"([A-Z]+):([A-Z0-9]+)", re.IGNORECASE)


def _tokenize_expression(expression: str) -> List[str]:
    """
    Tokenize a boolean expression into operators and ontology IDs.

    Parameters
    ----------
    expression
        Boolean expression string.

    Returns
    -------
    List of tokens: ontology IDs, operators (&, |, ~), and parentheses.
    """
    # Normalize whitespace
    expression = expression.strip()

    tokens = []
    i = 0
    while i < len(expression):
        char = expression[i]

        # Skip whitespace
        if char.isspace():
            i += 1
            continue

        # Operators and parentheses
        if char in "&|~()":
            tokens.append(char)
            i += 1
            continue

        # Ontology ID (starts with letters, contains colon)
        if char.isalpha():
            # Capture the full ontology ID
            match = ONTOLOGY_ID_PATTERN.match(expression[i:])
            if match:
                tokens.append(match.group(0).upper())
                i += len(match.group(0))
            else:
                raise ValueError(
                    f"Invalid token at position {i}: '{expression[i:i+10]}...'. "
                    "Expected ontology ID (e.g., CL:0000236, NCIT:C4349)."
                )
            continue

        raise ValueError(
            f"Unexpected character '{char}' at position {i} in expression: '{expression}'"
        )

    return tokens


def _find_ontology_columns(adata: ad.AnnData, ontology_columns: Optional[List[str]]) -> List[str]:
    """
    Find available ontology columns in AnnData.

    Parameters
    ----------
    adata
        AnnData object.
    ontology_columns
        Explicit list of columns to search, or None to use defaults.

    Returns
    -------
    List of available column names.
    """
    if ontology_columns is not None:
        # Validate provided columns exist
        available = [col for col in ontology_columns if col in adata.obs.columns]
        if not available:
            raise ValueError(
                f"None of the specified ontology columns found in adata.obs. "
                f"Requested: {ontology_columns}. "
                f"Available: {list(adata.obs.columns)}"
            )
        return available

    # Use defaults
    available = [col for col in DEFAULT_ONTOLOGY_COLUMNS if col in adata.obs.columns]
    if not available:
        raise ValueError(
            f"No ontology columns found in adata.obs. "
            f"Expected one of: {DEFAULT_ONTOLOGY_COLUMNS}. "
            f"Run add_ontology_ids() first to add ontology annotations."
        )
    return available


def _evaluate_single_id(
    ontology_id: str,
    adata: ad.AnnData,
    ontology_columns: List[str],
) -> np.ndarray:
    """
    Evaluate a single ontology ID against all cells.

    Parameters
    ----------
    ontology_id
        Ontology ID to match (e.g., "CL:0000236").
    adata
        AnnData object.
    ontology_columns
        Columns to search for the ID.

    Returns
    -------
    Boolean mask where True indicates cells matching the ontology ID.
    """
    # Normalize ID to uppercase
    ontology_id = ontology_id.upper()

    # Search across all ontology columns
    mask = np.zeros(adata.n_obs, dtype=bool)

    for col in ontology_columns:
        col_values = adata.obs[col].astype(str).str.upper()
        # Exact match
        mask |= col_values == ontology_id

    return mask


def _parse_and_evaluate(
    tokens: List[str],
    adata: ad.AnnData,
    ontology_columns: List[str],
) -> np.ndarray:
    """
    Parse and evaluate tokenized boolean expression using recursive descent.

    Grammar:
        expression := term (('|') term)*
        term := factor (('&') factor)*
        factor := '~' factor | '(' expression ')' | ONTOLOGY_ID

    Parameters
    ----------
    tokens
        List of tokens from _tokenize_expression.
    adata
        AnnData object.
    ontology_columns
        Columns to search.

    Returns
    -------
    Boolean mask for cells matching the expression.
    """
    pos = [0]  # Use list to allow mutation in nested functions

    def peek() -> Optional[str]:
        if pos[0] < len(tokens):
            return tokens[pos[0]]
        return None

    def consume() -> str:
        token = tokens[pos[0]]
        pos[0] += 1
        return token

    def parse_expression() -> np.ndarray:
        """Parse OR expressions."""
        result = parse_term()
        while peek() == "|":
            consume()  # consume '|'
            right = parse_term()
            result = result | right
        return result

    def parse_term() -> np.ndarray:
        """Parse AND expressions."""
        result = parse_factor()
        while peek() == "&":
            consume()  # consume '&'
            right = parse_factor()
            result = result & right
        return result

    def parse_factor() -> np.ndarray:
        """Parse NOT, parentheses, or ontology IDs."""
        token = peek()

        if token == "~":
            consume()  # consume '~'
            operand = parse_factor()
            return ~operand

        if token == "(":
            consume()  # consume '('
            result = parse_expression()
            if peek() != ")":
                raise ValueError("Missing closing parenthesis in expression")
            consume()  # consume ')'
            return result

        if token is None:
            raise ValueError("Unexpected end of expression")

        # Must be an ontology ID
        if ONTOLOGY_ID_PATTERN.match(token):
            consume()
            return _evaluate_single_id(token, adata, ontology_columns)

        raise ValueError(f"Unexpected token: '{token}'")

    result = parse_expression()

    if pos[0] < len(tokens):
        raise ValueError(
            f"Unexpected token after expression: '{tokens[pos[0]]}'"
        )

    return result


def evaluate_ontology_expression(
    expression: str,
    adata: ad.AnnData,
    ontology_columns: Optional[List[str]] = None,
) -> np.ndarray:
    """
    Evaluate boolean expression on ontology IDs.

    Supports boolean operators: & (AND), | (OR), ~ (NOT), and parentheses.

    Parameters
    ----------
    expression
        Boolean expression using ontology IDs.
        Examples:
        - "CL:0000236" - B cells
        - "CL:0000236 & NCIT:C4349" - B cells AND tumor
        - "CL:0000236 | CL:0000624" - B cells OR CD4+ T cells
        - "~NCIT:C4349" - NOT tumor
        - "(CL:0000236 | CL:0000624) & ~NCIT:C4349" - (B OR CD4+) AND NOT tumor
    adata
        AnnData object with ontology ID columns in .obs.
    ontology_columns
        Specific columns to search for ontology IDs. If None, searches
        default columns: cell_type_ontology_id, disease_ontology_id,
        tissue_ontology_id, etc.

    Returns
    -------
    np.ndarray
        Boolean mask (shape: n_obs) where True indicates cells matching
        the expression.

    Raises
    ------
    ValueError
        If expression syntax is invalid or no ontology columns found.

    Notes
    -----
    Ontology IDs are matched case-insensitively. The function searches
    across all specified ontology columns and returns True if any column
    matches for a cell.

    Examples
    --------
    >>> import anndata as ad
    >>> from spatialcore.annotation.expression import evaluate_ontology_expression
    >>> adata = ad.read_h5ad("annotated.h5ad")
    >>> # Find B cells
    >>> b_cell_mask = evaluate_ontology_expression("CL:0000236", adata)
    >>> adata_bcells = adata[b_cell_mask]
    >>> # Find B cells in tumor regions
    >>> mask = evaluate_ontology_expression("CL:0000236 & NCIT:C4349", adata)
    """
    if not expression or not expression.strip():
        raise ValueError("Expression cannot be empty")

    logger.debug(f"Evaluating ontology expression: {expression}")

    # Find available ontology columns
    available_columns = _find_ontology_columns(adata, ontology_columns)
    logger.debug(f"Using ontology columns: {available_columns}")

    # Tokenize expression
    tokens = _tokenize_expression(expression)
    logger.debug(f"Tokenized expression: {tokens}")

    # Parse and evaluate
    mask = _parse_and_evaluate(tokens, adata, available_columns)

    n_matching = mask.sum()
    logger.info(
        f"Ontology expression '{expression}' matched {n_matching:,} / {adata.n_obs:,} cells "
        f"({100 * n_matching / adata.n_obs:.1f}%)"
    )

    return mask


def get_ontology_ids_in_expression(expression: str) -> Set[str]:
    """
    Extract all ontology IDs from an expression.

    Parameters
    ----------
    expression
        Boolean expression string.

    Returns
    -------
    Set of ontology IDs found in the expression.

    Examples
    --------
    >>> get_ontology_ids_in_expression("CL:0000236 & NCIT:C4349")
    {'CL:0000236', 'NCIT:C4349'}
    """
    tokens = _tokenize_expression(expression)
    return {token.upper() for token in tokens if ONTOLOGY_ID_PATTERN.match(token)}
