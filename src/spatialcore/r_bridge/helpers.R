# ==============================================================================
# helpers.R - Internal Helper Functions
# ==============================================================================
#
# Internal functions for Seurat version compatibility and data conversion.
# These are not exported; they are used by convert.R.
#
# Note: seurat_to_adata() now uses anndataR for conversion, so most R->Python
# matrix conversion helpers are no longer needed. This file retains helpers for:
# - Seurat version detection (.get_seurat_version, .is_assay5)
# - BPCells detection (.check_bpcells)
# - Seurat data setters (.set_assay_data) for adata_to_seurat
# - Python->R matrix conversion (.from_python_matrix) for Python AnnData input
# - Spatial detection (.detect_spatial_reduction)
# - Validation utilities (.validate_names, .factors_to_character)
# - Python dict key conversion (.py_keys_to_r)
#
# ==============================================================================

# ------------------------------------------------------------------------------
# Seurat Version Detection
# ------------------------------------------------------------------------------

#' Get Seurat Major Version
#'
#' @return Integer major version (3, 4, or 5)
#' @keywords internal
.get_seurat_version <- function() {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop(
            "Package 'Seurat' is required but not installed.\n",
            "Install with: install.packages('Seurat')"
        )
    }
    version <- packageVersion("Seurat")
    as.integer(version$major)
}


#' Check if Assay is Seurat v5 Assay5
#'
#' @param assay_obj Seurat assay object
#' @return Logical
#' @keywords internal
.is_assay5 <- function(assay_obj) {
    inherits(assay_obj, "Assay5")
}


#' Set Assay Data (v3/v4/v5 Compatible)
#'
#' Handles the slot= (v3/v4) vs layer= (v5) API difference.
#'
#' @param seurat_obj Seurat object
#' @param assay Character. Assay name.
#' @param slot_name Character. One of "counts", "data", "scale.data"
#' @param new_data Matrix to set
#'
#' @return Modified Seurat object
#' @keywords internal
.set_assay_data <- function(seurat_obj, assay, slot_name, new_data) {

    assay_obj <- Seurat::GetAssay(seurat_obj, assay = assay)

    if (.is_assay5(assay_obj)) {
        # Seurat v5: use layer=
        layer_name <- switch(
            slot_name,
            "counts" = "counts",
            "data" = "data",
            "scale.data" = "scale.data",
            slot_name
        )
        return(Seurat::SetAssayData(seurat_obj, assay = assay, layer = layer_name, new.data = new_data))

    } else {
        # Seurat v3/v4: use slot=
        return(Seurat::SetAssayData(seurat_obj, assay = assay, slot = slot_name, new.data = new_data))
    }
}


# ------------------------------------------------------------------------------
# BPCells Detection
# ------------------------------------------------------------------------------

#' Check for BPCells On-Disk Matrices
#'
#' BPCells matrices are not supported. This function checks and errors if found.
#'
#' @param seurat_obj Seurat object
#' @param assay Character. Assay name.
#'
#' @return NULL (invisible). Raises error if BPCells detected.
#' @keywords internal
.check_bpcells <- function(seurat_obj, assay) {

    assay_obj <- Seurat::GetAssay(seurat_obj, assay = assay)

    if (!.is_assay5(assay_obj)) {
        # v3/v4 Assay objects don't have BPCells
        return(invisible(NULL))
    }

    # Check each layer for IterableMatrix (BPCells)
    layer_names <- SeuratObject::Layers(assay_obj)

    for (layer_name in layer_names) {
        # Access layer data - let errors propagate (fail loud)
        layer_data <- SeuratObject::LayerData(assay_obj, layer = layer_name)

        if (inherits(layer_data, "IterableMatrix")) {
            stop(
                "BPCells on-disk matrices (Seurat v5) are not yet supported.\n",
                "Layer '", layer_name, "' in assay '", assay, "' is an IterableMatrix.\n\n",
                "To convert to in-memory matrix:\n",
                "  seurat_obj[['", assay, "']]@layers$", layer_name, " <- ",
                "as.matrix(seurat_obj[['", assay, "']]@layers$", layer_name, ")\n\n",
                "Note: Future versions of SpatialCore will support BPCells directly."
            )
        }
    }

    invisible(NULL)
}


# ------------------------------------------------------------------------------
# Matrix Conversion: Python to R
# ------------------------------------------------------------------------------
# Note: R->Python matrix conversion (.to_python_matrix) has been removed
# since seurat_to_adata() now uses anndataR which handles this internally.

#' Convert Python Matrix to R (genes x cells)
#'
#' Converts Python matrix (cells x genes) to R format (genes x cells).
#' Handles both sparse and dense matrices.
#'
#' @param mat_py Python matrix (cells x genes)
#' @param gene_names Character vector of gene names
#' @param cell_names Character vector of cell names
#'
#' @return R matrix (genes x cells)
#' @keywords internal
.from_python_matrix <- function(mat_py, gene_names, cell_names) {

    # Check if sparse - check ALL classes, not just the first one
    all_classes <- class(mat_py)
    is_sparse <- any(grepl("sparse", all_classes, ignore.case = TRUE))

    if (is_sparse) {
        # Convert to CSC for efficient column operations in R
        mat_csc <- mat_py$tocsc()

        # Extract components
        data <- as.numeric(reticulate::py_to_r(mat_csc$data))
        indices <- as.integer(reticulate::py_to_r(mat_csc$indices))
        indptr <- as.integer(reticulate::py_to_r(mat_csc$indptr))
        shape <- as.integer(reticulate::py_to_r(mat_csc$shape))

        # Create dgCMatrix (cells x genes in R terms)
        mat_r <- Matrix::sparseMatrix(
            i = indices + 1L,  # 0-indexed to 1-indexed
            p = indptr,
            x = data,
            dims = shape,
            dimnames = list(cell_names, gene_names)
        )

        # Transpose to genes x cells
        return(Matrix::t(mat_r))

    } else {
        # Dense matrix
        mat_r <- reticulate::py_to_r(mat_py)

        # Transpose: cells x genes -> genes x cells
        mat_r <- t(mat_r)

        rownames(mat_r) <- gene_names
        colnames(mat_r) <- cell_names

        return(mat_r)
    }
}


# ------------------------------------------------------------------------------
# Spatial Reduction Detection
# ------------------------------------------------------------------------------

#' Detect Spatial Reduction in Seurat Object
#'
#' Auto-detects spatial coordinates from reduction names.
#' Checks for: "spatial", "SP", or names containing "spatial" (case-insensitive).
#'
#' @param seurat_obj Seurat object
#'
#' @return Character name of spatial reduction, or NULL if not found
#' @keywords internal
.detect_spatial_reduction <- function(seurat_obj) {

    red_names <- names(seurat_obj@reductions)

    if (length(red_names) == 0) {
        return(NULL)
    }

    # Priority 1: exact match "spatial"
    if ("spatial" %in% red_names) {
        return("spatial")
    }

    # Priority 2: exact match "SP" (CosMx convention)
    if ("SP" %in% red_names) {
        return("SP")
    }

    # Priority 3: case-insensitive "spatial"
    spatial_idx <- grep("^spatial$", red_names, ignore.case = TRUE)
    if (length(spatial_idx) > 0) {
        return(red_names[spatial_idx[1]])
    }

    # Priority 4: contains "spatial"
    spatial_idx <- grep("spatial", red_names, ignore.case = TRUE)
    if (length(spatial_idx) > 0) {
        return(red_names[spatial_idx[1]])
    }

    # Priority 5: starts with "SP"
    sp_idx <- grep("^SP", red_names)
    if (length(sp_idx) > 0) {
        return(red_names[sp_idx[1]])
    }

    return(NULL)
}


# ------------------------------------------------------------------------------
# Data Validation Helpers
# ------------------------------------------------------------------------------
# Note: .has_data() has been removed since seurat_to_adata() now uses anndataR
# which handles data validation internally.

#' Validate Name Alignment
#'
#' Checks that row/column names match expected values.
#' Raises error on mismatch.
#'
#' @param actual Character vector of actual names
#' @param expected Character vector of expected names
#' @param context Character describing what is being checked (for error message)
#'
#' @return NULL (invisible). Raises error on mismatch.
#' @keywords internal
.validate_names <- function(actual, expected, context) {

    if (!identical(actual, expected)) {
        # Find first difference for helpful error
        if (length(actual) != length(expected)) {
            stop(
                context, ": length mismatch.\n",
                "Expected ", length(expected), " names, got ", length(actual), "."
            )
        }

        diff_idx <- which(actual != expected)[1]
        stop(
            context, ": names do not match.\n",
            "First mismatch at position ", diff_idx, ":\n",
            "  Expected: '", expected[diff_idx], "'\n",
            "  Got: '", actual[diff_idx], "'"
        )
    }

    invisible(NULL)
}


#' Convert Factors to Character
#'
#' Converts all factor columns in a data.frame to character.
#' This prevents loss of factor labels when crossing R/Python boundary.
#'
#' @param df data.frame
#'
#' @return data.frame with factors converted to character
#' @keywords internal
.factors_to_character <- function(df) {

    if (ncol(df) == 0) {
        return(df)
    }

    for (col in names(df)) {
        if (is.factor(df[[col]])) {
            df[[col]] <- as.character(df[[col]])
        }
    }

    return(df)
}


# ------------------------------------------------------------------------------
# Python Object Conversion Helpers
# ------------------------------------------------------------------------------

#' Convert Python Dict Keys to R Character Vector
#'
#' Python dict_keys objects don't convert properly with py_to_r().
#' This function properly converts them to R character vectors.
#'
#' @param py_dict_or_mapping Python dict or mapping object (e.g., adata$layers)
#'
#' @return Character vector of keys, or empty character vector if NULL/empty
#' @keywords internal
.py_keys_to_r <- function(py_dict_or_mapping) {
    if (is.null(py_dict_or_mapping)) {
        return(character(0))
    }

    # Check if it's a valid Python object with keys method
    if (!inherits(py_dict_or_mapping, "python.builtin.object")) {
        return(character(0))
    }

    # Get keys and convert via Python list()
    keys_obj <- py_dict_or_mapping$keys()
    if (is.null(keys_obj)) {
        return(character(0))
    }

    # Use Python's list() to convert dict_keys to list, then to R
    # import_builtins() gives us access to Python built-in functions
    builtins <- reticulate::import_builtins()
    py_list <- builtins$list(keys_obj)
    result <- reticulate::py_to_r(py_list)

    if (is.null(result) || length(result) == 0) {
        return(character(0))
    }

    return(as.character(result))
}
