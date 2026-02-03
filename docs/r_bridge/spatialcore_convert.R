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
#   R packages: anndataR (Bioconductor), reticulate, Matrix, Seurat
#   Python packages: anndata, numpy, scipy, spatialcore (in conda env)
#
# Documentation: https://mcap91.github.io/SpatialCore/r_bridge/r_integration/
#
# ==============================================================================


# ------------------------------------------------------------------------------
# Internal Helpers
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
        layer_data <- SeuratObject::LayerData(assay_obj, layer = layer_name)
        if (inherits(layer_data, "IterableMatrix")) {
            stop("BPCells on-disk matrices (Seurat v5) are not yet supported.\n",
                 "Layer '", layer_name, "' is an IterableMatrix.\n",
                 "Convert to in-memory first or wait for future SpatialCore support.")
        }
    }
    invisible(NULL)
}

.from_python_matrix <- function(mat_py, gene_names, cell_names) {
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

.py_keys_to_r <- function(py_dict_or_mapping) {
    if (is.null(py_dict_or_mapping)) return(character(0))
    if (!inherits(py_dict_or_mapping, "python.builtin.object")) return(character(0))
    keys_obj <- py_dict_or_mapping$keys()
    if (is.null(keys_obj)) return(character(0))
    builtins <- reticulate::import_builtins()
    result <- reticulate::py_to_r(builtins$list(keys_obj))
    if (is.null(result) || length(result) == 0) return(character(0))
    as.character(result)
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

#' Convert Seurat Object to AnnData (using anndataR)
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
#' @return AnnData object (R-based from anndataR)
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
    # Check anndataR
    if (!requireNamespace("anndataR", quietly = TRUE)) {
        stop(
            "Package 'anndataR' is required but not installed.\n",
            "Install from Bioconductor with:\n",
            "  if (!require('BiocManager', quietly = TRUE))\n",
            "      install.packages('BiocManager')\n",
            "  BiocManager::install('anndataR')"
        )
    }

    if (!inherits(seurat_obj, "Seurat")) stop("Expected Seurat object, got: ", class(seurat_obj)[1])

    if (is.null(assay)) assay <- Seurat::DefaultAssay(seurat_obj)
    if (!assay %in% Seurat::Assays(seurat_obj)) {
        stop("Assay '", assay, "' not found. Available: ", paste(Seurat::Assays(seurat_obj), collapse = ", "))
    }
    x_slot <- match.arg(x_slot)
    .check_bpcells(seurat_obj, assay)

    # Build layers_mapping for anndataR
    layers_mapping <- NULL
    if (include_counts && x_slot == "data") {
        layers_mapping <- c("counts")
    }
    if (include_scale) {
        layers_mapping <- c(layers_mapping, "scale.data")
    }

    # Core conversion via anndataR
    adata <- anndataR::as_AnnData(
        seurat_obj,
        assay_name = assay,
        x_mapping = x_slot,
        layers_mapping = layers_mapping,
        obs_mapping = TRUE,
        var_mapping = TRUE,
        obsm_mapping = include_reductions,
        obsp_mapping = include_graphs,
        uns_mapping = TRUE
    )

    # Post-processing: Add SpatialCore-specific uns fields
    adata$uns[["X_slot"]] <- x_slot
    adata$uns[["seurat_assay"]] <- assay
    adata$uns[["seurat_project"]] <- seurat_obj@project.name
    adata$uns[["conversion"]] <- list(
        source = "seurat_to_adata",
        backend = "anndataR",
        seurat_version = as.character(packageVersion("Seurat")),
        anndataR_version = as.character(packageVersion("anndataR")),
        r_version = as.character(getRversion()),
        timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    )

    if (!is.null(graph_key_map)) {
        adata$uns[["graph_key_map"]] <- as.list(graph_key_map)
    }
    if (length(seurat_obj@misc) > 0) {
        adata$uns[["seurat_misc"]] <- seurat_obj@misc
    }

    # Handle spatial_key override
    if (!is.null(spatial_key)) {
        if (!spatial_key %in% names(seurat_obj@reductions)) {
            stop("Spatial reduction '", spatial_key, "' not found.\n",
                 "Available reductions: ", paste(names(seurat_obj@reductions), collapse = ", "))
        }
        red <- seurat_obj@reductions[[spatial_key]]
        adata$obsm[["spatial"]] <- red@cell.embeddings
    }

    # Auto-detect spatial if not already present
    obsm_keys <- names(adata$obsm)
    if (!"spatial" %in% obsm_keys) {
        detected <- .detect_spatial_reduction(seurat_obj)
        if (!is.null(detected) && detected %in% names(seurat_obj@reductions)) {
            red <- seurat_obj@reductions[[detected]]
            adata$obsm[["spatial"]] <- red@cell.embeddings
        }
    }

    # Handle graph_key_map renaming
    if (!is.null(graph_key_map) && include_graphs) {
        obsp_keys <- names(adata$obsp)
        for (old_key in names(graph_key_map)) {
            new_key <- graph_key_map[[old_key]]
            if (old_key %in% obsp_keys && old_key != new_key) {
                adata$obsp[[new_key]] <- adata$obsp[[old_key]]
                adata$obsp[[old_key]] <- NULL
            }
        }
    }

    # Get dimensions - anndataR uses methods, not properties
    n_cells <- nrow(adata$obs)
    n_genes <- nrow(adata$var)
    layer_keys <- names(adata$layers)
    obsm_keys <- names(adata$obsm)
    obsp_keys <- names(adata$obsp)

    message("Converted Seurat to AnnData (via anndataR):\n",
            "  Cells: ", n_cells, "\n",
            "  Genes: ", n_genes, "\n",
            "  X contains: ", x_slot, "\n",
            "  Layers: ", ifelse(length(layer_keys) > 0, paste(layer_keys, collapse = ", "), "none"), "\n",
            "  Reductions: ", ifelse(length(obsm_keys) > 0, paste(obsm_keys, collapse = ", "), "none"), "\n",
            "  Graphs: ", ifelse(length(obsp_keys) > 0, paste(obsp_keys, collapse = ", "), "none"))

    adata
}


#' Convert AnnData to Seurat Object
#'
#' @param adata AnnData object (R-based from anndataR or Python via reticulate)
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
    # Detect AnnData type
    is_r_anndata <- inherits(adata, "AnnData") || inherits(adata, "AbstractAnnData")
    is_py_anndata <- inherits(adata, "python.builtin.object")

    if (!is_r_anndata && !is_py_anndata) {
        stop("Expected an AnnData object, got R class: ", class(adata)[1], "\n",
             "Provide an anndataR AnnData or Python anndata.AnnData object.")
    }

    # Extract data based on AnnData type
    if (is_r_anndata) {
        obs_names <- rownames(adata$obs)
        var_names <- rownames(adata$var)
        n_cells <- length(obs_names)
        n_genes <- length(var_names)

        x_slot <- adata$uns[["X_slot"]]
        if (is.null(x_slot)) {
            stop("adata$uns['X_slot'] not found.\n",
                 "Cannot safely determine if X contains 'data' or 'counts'.\n",
                 "Set adata$uns['X_slot'] <- 'data' or 'counts' before conversion.")
        }
        if (!x_slot %in% c("data", "counts")) {
            stop("Invalid adata$uns['X_slot'] value: '", x_slot, "'. Expected 'data' or 'counts'.")
        }

        X_r <- adata$X
        if (is.null(X_r)) stop("adata$X is NULL or empty.")
        X_r <- Matrix::t(X_r)
        rownames(X_r) <- var_names
        colnames(X_r) <- obs_names

        layer_keys <- names(adata$layers)
        counts_r <- NULL
        data_r <- NULL
        scale_r <- NULL

        if (!is.null(layer_keys) && length(layer_keys) > 0) {
            if ("counts" %in% layer_keys) {
                counts_r <- Matrix::t(adata$layers[["counts"]])
                rownames(counts_r) <- var_names
                colnames(counts_r) <- obs_names
            }
            if ("data" %in% layer_keys) {
                data_r <- Matrix::t(adata$layers[["data"]])
                rownames(data_r) <- var_names
                colnames(data_r) <- obs_names
            }
            if ("scale" %in% layer_keys || "scale.data" %in% layer_keys) {
                scale_key <- if ("scale" %in% layer_keys) "scale" else "scale.data"
                scale_r <- Matrix::t(adata$layers[[scale_key]])
                rownames(scale_r) <- var_names
                colnames(scale_r) <- obs_names
            }
        }

        obs_df <- adata$obs
        var_df <- adata$var
        obsm_keys <- names(adata$obsm)
        obsm_data <- adata$obsm
        varm_keys <- names(adata$varm)
        varm_data <- adata$varm
        obsp_keys <- names(adata$obsp)
        obsp_data <- adata$obsp

    } else {
        # Python-based AnnData
        if (!requireNamespace("reticulate", quietly = TRUE)) {
            stop("Package 'reticulate' is required for Python AnnData objects.")
        }
        if (!reticulate::py_available()) {
            stop("Python is not available. Configure reticulate first.")
        }

        adata_class <- class(adata)[1]
        if (!grepl("AnnData", adata_class, ignore.case = TRUE)) {
            stop("Expected an AnnData object, got: ", adata_class)
        }

        obs_names <- reticulate::py_to_r(adata$obs_names$tolist())
        var_names <- reticulate::py_to_r(adata$var_names$tolist())
        n_cells <- length(obs_names)
        n_genes <- length(var_names)

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
        x_shape <- reticulate::py_to_r(adata$X$shape)
        if (x_shape[1] != n_cells || x_shape[2] != n_genes) {
            stop("adata.X shape does not match obs/var dimensions.")
        }

        X_r <- .from_python_matrix(adata$X, var_names, obs_names)

        layer_keys <- .py_keys_to_r(adata$layers)
        counts_r <- NULL
        data_r <- NULL
        scale_r <- NULL

        if (!is.null(layer_keys) && length(layer_keys) > 0) {
            if ("counts" %in% layer_keys) {
                counts_r <- .from_python_matrix(adata$layers[["counts"]], var_names, obs_names)
            }
            if ("data" %in% layer_keys) {
                data_r <- .from_python_matrix(adata$layers[["data"]], var_names, obs_names)
            }
            if ("scale" %in% layer_keys) {
                scale_r <- .from_python_matrix(adata$layers[["scale"]], var_names, obs_names)
            }
        }

        obs_df <- reticulate::py_to_r(adata$obs)
        if (is.null(rownames(obs_df))) rownames(obs_df) <- obs_names
        var_df <- reticulate::py_to_r(adata$var)
        if (is.null(rownames(var_df))) rownames(var_df) <- var_names

        obsm_keys <- .py_keys_to_r(adata$obsm)
        obsm_data <- adata$obsm
        varm_keys <- .py_keys_to_r(adata$varm)
        varm_data <- adata$varm
        obsp_keys <- .py_keys_to_r(adata$obsp)
        obsp_data <- adata$obsp
    }

    # Assign X based on x_slot
    if (x_slot == "counts") {
        if (is.null(counts_r)) counts_r <- X_r
    } else {
        if (is.null(data_r)) data_r <- X_r
    }

    if (is.null(counts_r)) {
        stop("No raw counts found for Seurat object creation.\n",
             "  - adata.uns['X_slot'] = '", x_slot, "' (so X is not counts)\n",
             "  - adata.layers['counts'] is missing\n\n",
             "To fix: set adata.layers['counts'] or adata.uns['X_slot'] = 'counts'")
    }

    if (is.null(rownames(obs_df))) rownames(obs_df) <- obs_names
    if (is.null(rownames(var_df))) rownames(var_df) <- var_names

    # Create or merge Seurat object
    if (is.null(seurat_obj)) {
        seurat_obj <- Seurat::CreateSeuratObject(counts = counts_r, project = project, assay = assay, meta.data = obs_df)
        if (!is.null(data_r)) seurat_obj <- .set_assay_data(seurat_obj, assay, "data", data_r)
        if (!is.null(scale_r)) {
            if (identical(rownames(scale_r), var_names) && identical(colnames(scale_r), obs_names)) {
                seurat_obj <- .set_assay_data(seurat_obj, assay, "scale.data", scale_r)
            }
        }
        if (ncol(var_df) > 0) {
            assay_obj <- Seurat::GetAssay(seurat_obj, assay = assay)
            if (identical(rownames(var_df), rownames(assay_obj))) {
                assay_obj[[names(var_df)]] <- var_df
                seurat_obj[[assay]] <- assay_obj
            }
        }
    } else {
        if (!inherits(seurat_obj, "Seurat")) stop("seurat_obj must be a Seurat object.")
        .validate_names(obs_names, colnames(seurat_obj), "Cell names")
        .validate_names(var_names, rownames(Seurat::GetAssay(seurat_obj, assay = assay)), "Gene names")
        for (col in colnames(obs_df)) seurat_obj[[col]] <- obs_df[[col]]
        if (!is.null(data_r)) seurat_obj <- .set_assay_data(seurat_obj, assay, "data", data_r)
        if (!is.null(scale_r)) {
            seurat_genes <- rownames(Seurat::GetAssay(seurat_obj, assay = assay))
            seurat_cells <- colnames(seurat_obj)
            if (identical(rownames(scale_r), seurat_genes) && identical(colnames(scale_r), seurat_cells)) {
                seurat_obj <- .set_assay_data(seurat_obj, assay, "scale.data", scale_r)
            }
        }
    }

    # Add reductions
    if (!is.null(obsm_keys) && length(obsm_keys) > 0) {
        for (obsm_key in obsm_keys) {
            if (is_r_anndata) {
                emb <- obsm_data[[obsm_key]]
            } else {
                emb <- reticulate::py_to_r(obsm_data[[obsm_key]])
            }
            if (nrow(emb) != n_cells) next
            red_name <- switch(obsm_key, "X_pca" = "pca", "X_umap" = "umap", "X_tsne" = "tsne", "spatial" = "spatial", gsub("^X_", "", obsm_key))
            key <- switch(red_name, "pca" = "PC_", "umap" = "UMAP_", "tsne" = "tSNE_", "spatial" = "Spatial_", paste0(toupper(substr(red_name, 1, 1)), "_"))
            rownames(emb) <- colnames(seurat_obj)
            colnames(emb) <- paste0(key, seq_len(ncol(emb)))
            seurat_obj[[red_name]] <- Seurat::CreateDimReducObject(embeddings = emb, key = key, assay = assay)
        }
    }

    # Add feature loadings
    if (!is.null(varm_keys) && length(varm_keys) > 0) {
        for (varm_key in varm_keys) {
            if (is_r_anndata) {
                loadings <- varm_data[[varm_key]]
            } else {
                loadings <- reticulate::py_to_r(varm_data[[varm_key]])
            }
            red_name <- switch(varm_key, "PCs" = "pca", gsub("_loadings$", "", varm_key))
            if (!red_name %in% names(seurat_obj@reductions)) next
            if (nrow(loadings) != n_genes) next
            rownames(loadings) <- rownames(Seurat::GetAssay(seurat_obj, assay = assay))
            seurat_obj[[red_name]]@feature.loadings <- loadings
        }
    }

    # Add graphs
    if (!is.null(obsp_keys) && length(obsp_keys) > 0) {
        for (obsp_key in obsp_keys) {
            if (is_r_anndata) {
                graph_r <- obsp_data[[obsp_key]]
                if (nrow(graph_r) != n_cells || ncol(graph_r) != n_cells) next
                if (!inherits(graph_r, "dgCMatrix")) graph_r <- as(graph_r, "dgCMatrix")
                rownames(graph_r) <- obs_names
                colnames(graph_r) <- obs_names
            } else {
                graph_py <- obsp_data[[obsp_key]]
                graph_shape <- reticulate::py_to_r(graph_py$shape)
                if (graph_shape[1] != n_cells || graph_shape[2] != n_cells) next
                graph_csc <- graph_py$tocsc()
                graph_r <- Matrix::sparseMatrix(
                    i = as.integer(reticulate::py_to_r(graph_csc$indices)) + 1L,
                    p = as.integer(reticulate::py_to_r(graph_csc$indptr)),
                    x = as.numeric(reticulate::py_to_r(graph_csc$data)),
                    dims = c(n_cells, n_cells),
                    dimnames = list(obs_names, obs_names))
                graph_r <- as(graph_r, "dgCMatrix")
            }
            graph_name <- switch(obsp_key, "connectivities" = paste0(assay, "_nn"), "distances" = paste0(assay, "_snn"), obsp_key)
            seurat_obj@graphs[[graph_name]] <- as(graph_r, "Graph")
        }
    }

    message("Converted AnnData to Seurat:\n",
            "  Cells: ", ncol(seurat_obj), "\n",
            "  Genes: ", nrow(seurat_obj), "\n",
            "  Metadata columns: ", ncol(seurat_obj[[]]), "\n",
            "  Reductions: ", paste(names(seurat_obj@reductions), collapse = ", "), "\n",
            "  Graphs: ", paste(names(seurat_obj@graphs), collapse = ", "))

    seurat_obj
}


# ==============================================================================
# End of spatialcore_convert.R
# ==============================================================================
