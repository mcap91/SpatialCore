# ==============================================================================
# spatialcore_convert.R - Seurat <-> AnnData Conversion for SpatialCore
# ==============================================================================
#
# Standalone script for converting between Seurat and AnnData objects.
# Source this file to use SpatialCore functions from R.
#
# Usage:
#   source("spatialcore_convert.R")
#   adata <- seurat_to_adata(seurat_obj)
#   spc <- reticulate::import("spatialcore")
#   spc$spatial$compute_morans_i(adata)
#   seurat_obj <- adata_to_seurat(adata, seurat_obj)
#
# Requirements:
#   R packages: reticulate, Matrix, Seurat
#   Python packages: anndata, numpy, scipy, spatialcore (in conda env)
#
# Documentation: https://mcap91.github.io/SpatialCore/r_bridge/r_integration/
#
# ==============================================================================


# ------------------------------------------------------------------------------
# Internal Helpers (with all Phase 1 & Phase 2 fixes applied)
# ------------------------------------------------------------------------------

.get_seurat_version <- function() {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' is required but not installed.\n",
             "Install with: install.packages('Seurat')")
    }
    as.integer(packageVersion("Seurat")$major)
}

.is_assay5 <- function(assay_obj) {
    inherits(assay_obj, "Assay5")
}

.get_assay_data <- function(seurat_obj, assay, slot_name) {
    assay_obj <- Seurat::GetAssay(seurat_obj, assay = assay)
    if (.is_assay5(assay_obj)) {
        layer_name <- switch(slot_name,
            "counts" = "counts", "data" = "data", "scale.data" = "scale.data",
            slot_name)
        available_layers <- SeuratObject::Layers(assay_obj)
        if (!layer_name %in% available_layers) {
            stop("Layer '", layer_name, "' not found in assay '", assay, "'.\n",
                 "Available layers: ", paste(available_layers, collapse = ", "), "\n",
                 "Ensure the required data slot exists before conversion.")
        }
        return(Seurat::GetAssayData(seurat_obj, assay = assay, layer = layer_name))
    }
    Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot_name)
}

.set_assay_data <- function(seurat_obj, assay, slot_name, new_data) {
    assay_obj <- Seurat::GetAssay(seurat_obj, assay = assay)
    if (.is_assay5(assay_obj)) {
        layer_name <- switch(slot_name,
            "counts" = "counts", "data" = "data", "scale.data" = "scale.data",
            slot_name)
        return(Seurat::SetAssayData(seurat_obj, assay = assay, layer = layer_name, new.data = new_data))
    }
    Seurat::SetAssayData(seurat_obj, assay = assay, slot = slot_name, new.data = new_data)
}

.check_bpcells <- function(seurat_obj, assay) {
    assay_obj <- Seurat::GetAssay(seurat_obj, assay = assay)
    if (!.is_assay5(assay_obj)) return(invisible(NULL))
    for (layer_name in SeuratObject::Layers(assay_obj)) {
        # Let errors propagate (fail loud)
        layer_data <- SeuratObject::LayerData(assay_obj, layer = layer_name)
        if (inherits(layer_data, "IterableMatrix")) {
            stop("BPCells on-disk matrices (Seurat v5) are not yet supported.\n",
                 "Layer '", layer_name, "' is an IterableMatrix.\n",
                 "Convert to in-memory first or wait for future SpatialCore support.")
        }
    }
    invisible(NULL)
}

.to_python_matrix <- function(mat_r, np, scipy_sparse) {
    mat_r <- Matrix::t(mat_r)
    if (inherits(mat_r, "dgCMatrix")) {
        if (length(mat_r@x) == 0) {
            return(scipy_sparse$csr_matrix(
                reticulate::tuple(as.integer(nrow(mat_r)), as.integer(ncol(mat_r))),
                dtype = "float64"))
        }
        # Convert CSC (dgCMatrix) to CSR (RsparseMatrix) for correct component extraction
        mat_csr <- as(mat_r, "RsparseMatrix")
        return(scipy_sparse$csr_matrix(
            reticulate::tuple(
                np$array(as.numeric(mat_csr@x), dtype = "float64"),
                np$array(as.integer(mat_csr@j), dtype = "int32"),
                np$array(as.integer(mat_csr@p), dtype = "int32")),
            shape = reticulate::tuple(as.integer(nrow(mat_csr)), as.integer(ncol(mat_csr)))))
    }
    np$array(as.matrix(mat_r), dtype = "float64")
}

.from_python_matrix <- function(mat_py, gene_names, cell_names) {
    # Check ALL classes for sparse, not just the first one
    all_classes <- class(mat_py)
    is_sparse <- any(grepl("sparse", all_classes, ignore.case = TRUE))
    if (is_sparse) {
        mat_csc <- mat_py$tocsc()
        mat_r <- Matrix::sparseMatrix(
            i = as.integer(reticulate::py_to_r(mat_csc$indices)) + 1L,
            p = as.integer(reticulate::py_to_r(mat_csc$indptr)),
            x = as.numeric(reticulate::py_to_r(mat_csc$data)),
            dims = as.integer(reticulate::py_to_r(mat_csc$shape)),
            dimnames = list(cell_names, gene_names))
        return(Matrix::t(mat_r))
    }
    mat_r <- t(reticulate::py_to_r(mat_py))
    rownames(mat_r) <- gene_names
    colnames(mat_r) <- cell_names
    mat_r
}

.detect_spatial_reduction <- function(seurat_obj) {
    red_names <- names(seurat_obj@reductions)
    if (length(red_names) == 0) return(NULL)
    if ("spatial" %in% red_names) return("spatial")
    if ("SP" %in% red_names) return("SP")
    idx <- grep("spatial", red_names, ignore.case = TRUE)
    if (length(idx) > 0) return(red_names[idx[1]])
    idx <- grep("^SP", red_names)
    if (length(idx) > 0) return(red_names[idx[1]])
    NULL
}

.has_data <- function(mat) {
    if (is.null(mat) || nrow(mat) == 0 || ncol(mat) == 0) return(FALSE)
    # Check all sparse matrix types, not just dgCMatrix
    if (inherits(mat, "sparseMatrix")) return(length(mat@x) > 0)
    # For dense: use sum(abs()) which handles NA and is efficient
    sum(abs(mat), na.rm = TRUE) > 0
}

.validate_names <- function(actual, expected, context) {
    if (!identical(actual, expected)) {
        if (length(actual) != length(expected)) {
            stop(context, ": length mismatch. Expected ", length(expected), ", got ", length(actual), ".")
        }
        diff_idx <- which(actual != expected)[1]
        stop(context, ": mismatch at position ", diff_idx, ". Expected '", expected[diff_idx], "', got '", actual[diff_idx], "'.")
    }
    invisible(NULL)
}

.factors_to_character <- function(df) {
    if (ncol(df) == 0) return(df)
    for (col in names(df)) {
        if (is.factor(df[[col]])) df[[col]] <- as.character(df[[col]])
    }
    df
}

# Helper for graph conversion (cells x cells, no transpose needed)
.graph_to_python <- function(graph_mat, np, scipy_sparse, n_cells) {
    if (length(graph_mat@x) == 0) {
        return(scipy_sparse$csr_matrix(
            reticulate::tuple(as.integer(n_cells), as.integer(n_cells)),
            dtype = "float64"))
    }
    mat_csr <- as(graph_mat, "RsparseMatrix")
    scipy_sparse$csr_matrix(
        reticulate::tuple(
            np$array(as.numeric(mat_csr@x), dtype = "float64"),
            np$array(as.integer(mat_csr@j), dtype = "int32"),
            np$array(as.integer(mat_csr@p), dtype = "int32")),
        shape = reticulate::tuple(as.integer(n_cells), as.integer(n_cells)))
}

# Helper for graph conversion from Python (cells x cells, no transpose)
.graph_from_python <- function(graph_py, obs_names, n_cells) {
    graph_csc <- graph_py$tocsc()
    data <- as.numeric(reticulate::py_to_r(graph_csc$data))
    indices <- as.integer(reticulate::py_to_r(graph_csc$indices))
    indptr <- as.integer(reticulate::py_to_r(graph_csc$indptr))
    graph_r <- Matrix::sparseMatrix(
        i = indices + 1L, p = indptr, x = data,
        dims = c(n_cells, n_cells),
        dimnames = list(obs_names, obs_names))
    as(graph_r, "dgCMatrix")
}

# Helper to convert Python dict_keys to R vector
.py_keys_to_r <- function(keys_obj) {
    if (is.null(keys_obj)) return(character(0))
    builtins <- reticulate::import_builtins()
    reticulate::py_to_r(builtins$list(keys_obj))
}


# ------------------------------------------------------------------------------
# Setup Functions
# ------------------------------------------------------------------------------

#' Setup SpatialCore Environment
#' @export
setup_spatialcore <- function(conda_env = "spatialcore", verbose = TRUE) {
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("Package 'reticulate' required. Install with: install.packages('reticulate')")
    }
    envs <- tryCatch(reticulate::conda_list(), error = function(e) {
        stop("Could not list conda environments. Ensure conda/mamba is installed.")
    })
    if (!conda_env %in% envs$name) {
        stop("Conda environment '", conda_env, "' not found.\n",
             "Available: ", paste(envs$name, collapse = ", "))
    }
    reticulate::use_condaenv(conda_env, required = TRUE)
    if (!reticulate::py_available(initialize = TRUE)) {
        stop("Python not available after activating '", conda_env, "'.")
    }
    required <- c("anndata", "numpy", "scipy")
    missing <- required[!sapply(required, reticulate::py_module_available)]
    if (length(missing) > 0) {
        stop("Required Python packages not found: ", paste(missing, collapse = ", "))
    }
    if (verbose) {
        message("SpatialCore ready. Python: ", reticulate::py_config()$python)
    }
    invisible(TRUE)
}

#' Check if SpatialCore Environment is Ready
#' @export
is_spatialcore_ready <- function(conda_env = "spatialcore") {
    if (!requireNamespace("reticulate", quietly = TRUE)) return(FALSE)
    if (!reticulate::py_available()) return(FALSE)
    reticulate::py_module_available("anndata")
}

#' Get SpatialCore Python Module
#' @export
get_spatialcore <- function() {
    if (!reticulate::py_module_available("spatialcore")) {
        stop("Python package 'spatialcore' not found. Install with: pip install spatialcore")
    }
    reticulate::import("spatialcore")
}


# ------------------------------------------------------------------------------
# Conversion Functions
# ------------------------------------------------------------------------------

#' Convert Seurat Object to AnnData
#'
#' @param seurat_obj A Seurat object (v3, v4, or v5)
#' @param assay Assay name. Default: DefaultAssay(seurat_obj)
#' @param x_slot What to store in X: "data" (normalized) or "counts" (raw)
#' @param include_counts Include raw counts in layers["counts"]
#' @param include_scale Include scaled data in layers["scale"]
#' @param include_reductions Include PCA, UMAP, spatial, etc.
#' @param include_graphs Include neighbor graphs
#' @param graph_key_map Named vector mapping Seurat graph names to obsp keys
#' @param spatial_key Override spatial reduction name (auto-detected if NULL)
#'
#' @return AnnData object (Python object via reticulate)
#' @export
seurat_to_adata <- function(
    seurat_obj,
    assay = NULL,
    x_slot = c("data", "counts"),
    include_counts = TRUE,
    include_scale = FALSE,
    include_reductions = TRUE,
    include_graphs = TRUE,
    graph_key_map = NULL,
    spatial_key = NULL
) {
    if (!requireNamespace("reticulate", quietly = TRUE)) stop("Package 'reticulate' required.")
    if (!reticulate::py_available()) stop("Python not available. Run setup_spatialcore() first.")
    if (!inherits(seurat_obj, "Seurat")) stop("Expected Seurat object, got: ", class(seurat_obj)[1])

    if (is.null(assay)) assay <- Seurat::DefaultAssay(seurat_obj)
    if (!assay %in% Seurat::Assays(seurat_obj)) {
        stop("Assay '", assay, "' not found. Available: ", paste(Seurat::Assays(seurat_obj), collapse = ", "))
    }
    x_slot <- match.arg(x_slot)
    .check_bpcells(seurat_obj, assay)

    anndata <- reticulate::import("anndata", convert = FALSE)
    np <- reticulate::import("numpy", convert = FALSE)
    scipy_sparse <- reticulate::import("scipy.sparse", convert = FALSE)

    x_data_r <- .get_assay_data(seurat_obj, assay, x_slot)
    if (!.has_data(x_data_r)) stop("No data in '", x_slot, "' slot for assay '", assay, "'.")
    X_py <- .to_python_matrix(x_data_r, np, scipy_sparse)

    layers <- list()
    if (include_counts && x_slot == "data") {
        counts_r <- tryCatch(
            .get_assay_data(seurat_obj, assay, "counts"),
            error = function(e) NULL
        )
        if (!is.null(counts_r) && .has_data(counts_r)) {
            layers[["counts"]] <- .to_python_matrix(counts_r, np, scipy_sparse)
        } else {
            warning("Raw counts not available in assay '", assay, "'.\n",
                    "layers['counts'] will not be created.")
        }
    }
    if (include_scale) {
        scale_r <- tryCatch(
            .get_assay_data(seurat_obj, assay, "scale.data"),
            error = function(e) NULL
        )
        if (!is.null(scale_r) && .has_data(scale_r)) {
            full_genes <- rownames(Seurat::GetAssay(seurat_obj, assay = assay))
            if (identical(rownames(scale_r), full_genes)) {
                layers[["scale"]] <- .to_python_matrix(scale_r, np, scipy_sparse)
            } else {
                warning("scale.data does not cover all features (", nrow(scale_r),
                        " vs ", length(full_genes), "). layers['scale'] not created.")
            }
        } else {
            warning("include_scale=TRUE but scale.data is empty in assay '", assay, "'.")
        }
    }

    obs_df <- .factors_to_character(seurat_obj[[]])
    cell_names <- colnames(seurat_obj)
    rownames(obs_df) <- cell_names

    assay_obj <- Seurat::GetAssay(seurat_obj, assay = assay)
    var_df <- .factors_to_character(assay_obj[[]])
    gene_names <- rownames(assay_obj)
    rownames(var_df) <- gene_names

    obsm <- list()
    varm <- list()
    n_cells <- length(cell_names)
    if (include_reductions && length(seurat_obj@reductions) > 0) {
        for (red_name in names(seurat_obj@reductions)) {
            red <- seurat_obj@reductions[[red_name]]
            emb <- red@cell.embeddings
            if (nrow(emb) == 0) {
                warning("Reduction '", red_name, "' has no cell embeddings. Skipping.")
                next
            }
            obsm_key <- switch(tolower(red_name),
                "pca" = "X_pca", "umap" = "X_umap", "tsne" = "X_tsne",
                "spatial" = "spatial", "sp" = "spatial", paste0("X_", red_name))
            if (obsm_key %in% names(obsm)) {
                stop("obsm key collision: '", obsm_key, "' already exists.")
            }
            obsm[[obsm_key]] <- np$array(emb, dtype = "float64")
            if (!is.null(red@feature.loadings) && nrow(red@feature.loadings) > 0) {
                if (identical(rownames(red@feature.loadings), gene_names)) {
                    varm_key <- if (tolower(red_name) == "pca") "PCs" else paste0(red_name, "_loadings")
                    if (!varm_key %in% names(varm)) varm[[varm_key]] <- np$array(as.matrix(red@feature.loadings), dtype = "float64")
                } else {
                    warning("Reduction '", red_name, "' loadings do not cover all genes. Skipping varm.")
                }
            }
        }
    }
    if (!is.null(spatial_key) && spatial_key %in% names(seurat_obj@reductions)) {
        obsm[["spatial"]] <- np$array(seurat_obj@reductions[[spatial_key]]@cell.embeddings, dtype = "float64")
    } else if (!"spatial" %in% names(obsm)) {
        detected <- .detect_spatial_reduction(seurat_obj)
        if (!is.null(detected)) obsm[["spatial"]] <- np$array(seurat_obj@reductions[[detected]]@cell.embeddings, dtype = "float64")
    }

    obsp <- list()
    if (include_graphs && length(seurat_obj@graphs) > 0) {
        graph_names <- names(seurat_obj@graphs)
        nn_matches <- grep("_nn$", graph_names, value = TRUE)
        snn_matches <- grep("_snn$", graph_names, value = TRUE)
        if ((length(nn_matches) > 1 || length(snn_matches) > 1) && is.null(graph_key_map)) {
            stop("Multiple graphs match _nn/_snn. Provide graph_key_map to disambiguate.\n",
                 "  _nn graphs: ", paste(nn_matches, collapse = ", "), "\n",
                 "  _snn graphs: ", paste(snn_matches, collapse = ", "))
        }
        for (gname in graph_names) {
            graph <- seurat_obj@graphs[[gname]]
            if (!inherits(graph, "Graph")) {
                warning("Graph '", gname, "' is not a Seurat Graph object. Skipping.")
                next
            }
            if (!is.null(graph_key_map) && gname %in% names(graph_key_map)) {
                obsp_key <- graph_key_map[[gname]]
            } else if (gname %in% nn_matches && length(nn_matches) == 1) {
                obsp_key <- "connectivities"
            } else if (gname %in% snn_matches && length(snn_matches) == 1) {
                obsp_key <- "distances"
            } else {
                obsp_key <- gname
            }
            if (obsp_key %in% names(obsp)) {
                stop("obsp key collision: '", obsp_key, "' already exists.")
            }
            graph_mat <- as(graph, "dgCMatrix")
            if (nrow(graph_mat) != n_cells || ncol(graph_mat) != n_cells) {
                stop("Graph '", gname, "' dimensions do not match cell count.")
            }
            # Use graph-specific converter (no transpose for square cell x cell)
            obsp[[obsp_key]] <- .graph_to_python(graph_mat, np, scipy_sparse, n_cells)
        }
    }

    uns <- list(
        "X_slot" = x_slot,
        "seurat_assay" = assay,
        "seurat_project" = seurat_obj@project.name,
        "conversion" = list(source = "seurat_to_adata", timestamp = format(Sys.time()))
    )
    if (!is.null(graph_key_map)) uns[["graph_key_map"]] <- as.list(graph_key_map)
    if (length(seurat_obj@misc) > 0) uns[["seurat_misc"]] <- seurat_obj@misc

    adata <- anndata$AnnData(
        X = X_py, obs = obs_df, var = var_df,
        layers = if (length(layers) > 0) reticulate::r_to_py(layers) else NULL,
        obsm = if (length(obsm) > 0) reticulate::r_to_py(obsm) else NULL,
        obsp = if (length(obsp) > 0) reticulate::r_to_py(obsp) else NULL,
        varm = if (length(varm) > 0) reticulate::r_to_py(varm) else NULL,
        uns = reticulate::r_to_py(uns)
    )

    message("Converted: ", length(cell_names), " cells, ", length(gene_names), " genes")
    adata
}


#' Convert AnnData to Seurat Object
#'
#' @param adata AnnData object (Python object via reticulate)
#' @param seurat_obj Optional existing Seurat object to merge into
#' @param assay Assay name for new Seurat object
#' @param project Project name for new Seurat object
#'
#' @return Seurat object
#' @export
adata_to_seurat <- function(
    adata,
    seurat_obj = NULL,
    assay = "RNA",
    project = "SpatialCore"
) {
    if (!requireNamespace("reticulate", quietly = TRUE)) stop("Package 'reticulate' required.")
    if (!reticulate::py_available()) stop("Python not available.")
    if (!inherits(adata, "python.builtin.object")) stop("Expected Python AnnData object.")

    obs_names <- reticulate::py_to_r(adata$obs_names$tolist())
    var_names <- reticulate::py_to_r(adata$var_names$tolist())
    n_cells <- length(obs_names)
    n_genes <- length(var_names)

    # Get X slot semantics - REQUIRED for safe conversion
    x_slot <- NULL
    if (!is.null(adata$uns)) {
        x_slot_val <- reticulate::py_to_r(adata$uns$get("X_slot"))
        if (!is.null(x_slot_val)) x_slot <- x_slot_val
    }
    if (is.null(x_slot)) {
        stop("adata.uns['X_slot'] not found.\n",
             "Cannot safely determine if X contains 'data' or 'counts'.\n",
             "Set adata.uns['X_slot'] = 'data' or 'counts' before conversion.")
    }
    if (!x_slot %in% c("data", "counts")) {
        stop("Invalid adata.uns['X_slot'] value: '", x_slot, "'. Expected 'data' or 'counts'.")
    }

    if (is.null(adata$X)) stop("adata.X is empty.")
    X_r <- .from_python_matrix(adata$X, var_names, obs_names)

    # Get layers safely using .py_keys_to_r helper
    layer_keys <- character(0)
    if (!is.null(adata$layers) && length(reticulate::py_to_r(adata$layers)) > 0) {
        layer_keys <- .py_keys_to_r(adata$layers$keys())
    }
    counts_r <- NULL
    data_r <- NULL
    if (length(layer_keys) > 0 && "counts" %in% layer_keys) {
        counts_r <- .from_python_matrix(adata$layers[["counts"]], var_names, obs_names)
    }

    if (x_slot == "counts") {
        if (is.null(counts_r)) counts_r <- X_r
    } else {
        data_r <- X_r
    }

    # Ensure we have counts - REQUIRED for Seurat
    if (is.null(counts_r)) {
        stop("No raw counts found for Seurat object creation.\n",
             "  - adata.uns['X_slot'] = '", x_slot, "' (so X is not counts)\n",
             "  - adata.layers['counts'] is missing\n\n",
             "To fix: set adata.layers['counts'] or adata.uns['X_slot'] = 'counts'")
    }

    obs_df <- reticulate::py_to_r(adata$obs)
    rownames(obs_df) <- obs_names
    var_df <- reticulate::py_to_r(adata$var)
    rownames(var_df) <- var_names

    if (is.null(seurat_obj)) {
        seurat_obj <- Seurat::CreateSeuratObject(counts = counts_r, project = project, assay = assay, meta.data = obs_df)
        if (!is.null(data_r)) seurat_obj <- .set_assay_data(seurat_obj, assay, "data", data_r)
    } else {
        if (!inherits(seurat_obj, "Seurat")) stop("seurat_obj must be a Seurat object.")
        .validate_names(obs_names, colnames(seurat_obj), "Cell names")
        .validate_names(var_names, rownames(Seurat::GetAssay(seurat_obj, assay = assay)), "Gene names")
        for (col in colnames(obs_df)) seurat_obj[[col]] <- obs_df[[col]]
        if (!is.null(data_r)) seurat_obj <- .set_assay_data(seurat_obj, assay, "data", data_r)
    }

    # Get obsm keys safely using .py_keys_to_r helper
    obsm_keys <- character(0)
    if (!is.null(adata$obsm) && length(reticulate::py_to_r(adata$obsm)) > 0) {
        obsm_keys <- .py_keys_to_r(adata$obsm$keys())
    }
    if (length(obsm_keys) > 0) {
        for (obsm_key in obsm_keys) {
            emb <- reticulate::py_to_r(adata$obsm[[obsm_key]])
            if (nrow(emb) != n_cells) {
                warning("obsm['", obsm_key, "'] has ", nrow(emb), " rows but ", n_cells, " expected. Skipping.")
                next
            }
            red_name <- switch(obsm_key, "X_pca" = "pca", "X_umap" = "umap", "X_tsne" = "tsne", "spatial" = "spatial", gsub("^X_", "", obsm_key))
            key <- switch(red_name, "pca" = "PC_", "umap" = "UMAP_", "tsne" = "tSNE_", "spatial" = "Spatial_", paste0(toupper(substr(red_name, 1, 1)), "_"))
            rownames(emb) <- colnames(seurat_obj)
            colnames(emb) <- paste0(key, seq_len(ncol(emb)))
            seurat_obj[[red_name]] <- Seurat::CreateDimReducObject(embeddings = emb, key = key, assay = assay)
        }
    }

    # Get varm keys safely using .py_keys_to_r helper
    varm_keys <- character(0)
    if (!is.null(adata$varm) && length(reticulate::py_to_r(adata$varm)) > 0) {
        varm_keys <- .py_keys_to_r(adata$varm$keys())
    }
    if (length(varm_keys) > 0) {
        for (varm_key in varm_keys) {
            loadings <- reticulate::py_to_r(adata$varm[[varm_key]])
            red_name <- switch(varm_key, "PCs" = "pca", gsub("_loadings$", "", varm_key))
            if (!red_name %in% names(seurat_obj@reductions)) {
                warning("varm['", varm_key, "'] maps to reduction '", red_name, "' which does not exist. Skipping.")
                next
            }
            if (nrow(loadings) != n_genes) {
                warning("varm['", varm_key, "'] has ", nrow(loadings), " rows but ", n_genes, " expected. Skipping.")
                next
            }
            rownames(loadings) <- rownames(Seurat::GetAssay(seurat_obj, assay = assay))
            seurat_obj[[red_name]]@feature.loadings <- loadings
        }
    }

    # Get obsp keys safely using .py_keys_to_r helper
    obsp_keys <- character(0)
    if (!is.null(adata$obsp) && length(reticulate::py_to_r(adata$obsp)) > 0) {
        obsp_keys <- .py_keys_to_r(adata$obsp$keys())
    }
    if (length(obsp_keys) > 0) {
        for (obsp_key in obsp_keys) {
            graph_py <- adata$obsp[[obsp_key]]
            graph_shape <- reticulate::py_to_r(graph_py$shape)
            if (graph_shape[1] != n_cells || graph_shape[2] != n_cells) {
                warning("obsp['", obsp_key, "'] shape does not match n_cells. Skipping.")
                next
            }
            # Use graph-specific converter (no transpose for square cell x cell)
            graph_r <- .graph_from_python(graph_py, obs_names, n_cells)
            graph_name <- switch(obsp_key, "connectivities" = paste0(assay, "_nn"), "distances" = paste0(assay, "_snn"), obsp_key)
            seurat_obj@graphs[[graph_name]] <- as(graph_r, "Graph")
        }
    }

    message("Converted: ", ncol(seurat_obj), " cells, ", nrow(seurat_obj), " genes")
    seurat_obj
}


# ==============================================================================
# End of spatialcore_convert.R
# ==============================================================================
