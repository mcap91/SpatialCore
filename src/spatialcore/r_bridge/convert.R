# ==============================================================================
# convert.R - Seurat <-> AnnData Conversion Functions
# ==============================================================================
#
# Bidirectional conversion between Seurat (R) and AnnData (Python) objects.
# Enables R users to use any SpatialCore function via reticulate.
#
# Dependencies:
#   - helpers.R (internal functions)
#   - reticulate, Matrix, Seurat packages
#   - Python: anndata, numpy, scipy
#
# ==============================================================================

# Source helpers if not already loaded
# (In package context, this is handled by NAMESPACE/Collate)


# ==============================================================================
# seurat_to_adata: Seurat -> AnnData
# ==============================================================================

#' Convert Seurat Object to AnnData
#'
#' Converts a Seurat object to a Python AnnData object in memory via reticulate.
#' The resulting AnnData can be used with any SpatialCore or scanpy function.
#'
#' @param seurat_obj A Seurat object (v3, v4, or v5)
#' @param assay Character. Name of the assay to convert. Default: DefaultAssay(seurat_obj)
#' @param x_slot Character. What to store in X: "data" (normalized, default) or "counts" (raw).
#' @param include_counts Logical. Include raw counts in layers["counts"]. Default: TRUE
#' @param include_scale Logical. Include scaled data in layers["scale"]. Default: FALSE
#' @param include_reductions Logical. Include dimensionality reductions. Default: TRUE
#' @param include_graphs Logical. Include neighbor graphs. Default: TRUE
#' @param graph_key_map Named character vector. Explicit mapping of Seurat graph names to obsp keys.
#'   Required when multiple graphs match _nn or _snn pattern.
#' @param spatial_key Character or NULL. Name of spatial reduction. NULL = auto-detect.
#'
#' @return An AnnData object (Python object via reticulate)
#'
#' @details
#' ## Data Mapping
#' By default, normalized data goes to X and raw counts go to layers["counts"].
#' This matches scanpy/SpatialCore conventions where X contains processed values.
#'
#' The `uns["X_slot"]` field records whether X contains "data" (normalized) or
#' "counts" (raw), enabling safe round-trip conversion.
#'
#' ## Seurat Version Compatibility
#' This function supports Seurat v3, v4, and v5 objects. BPCells on-disk
#' matrices (Seurat v5) are NOT supported and will raise an error.
#'
#' ## Graph Handling
#' Graphs are auto-mapped when unambiguous:
#' - Single `*_nn` graph -> `obsp["connectivities"]`
#' - Single `*_snn` graph -> `obsp["distances"]`
#'
#' If multiple graphs match these patterns, you must provide `graph_key_map`
#' to explicitly specify the mapping.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(reticulate)
#' use_condaenv("spatialcore")
#'
#' # Basic conversion
#' adata <- seurat_to_adata(seurat_obj)
#'
#' # Use raw counts in X instead of normalized
#' adata <- seurat_to_adata(seurat_obj, x_slot = "counts")
#'
#' # Include scaled data
#' adata <- seurat_to_adata(seurat_obj, include_scale = TRUE)
#'
#' # Explicit graph mapping
#' adata <- seurat_to_adata(
#'     seurat_obj,
#'     graph_key_map = c("RNA_nn" = "connectivities", "RNA_snn" = "distances")
#' )
#'
#' # Use with SpatialCore
#' sc <- import("spatialcore")
#' sc$spatial$compute_morans_i(adata)
#' }
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

    # =========================================================================
    # Input Validation
    # =========================================================================

    # Check reticulate
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop(
            "Package 'reticulate' is required but not installed.\n",
            "Install with: install.packages('reticulate')"
        )
    }

    # Check Python
    if (!reticulate::py_available()) {
        stop(
            "Python is not available.\n",
            "Configure reticulate first:\n",
            "  reticulate::use_condaenv('spatialcore')\n",
            "Or run: setup_spatialcore()"
        )
    }

    # Check Seurat object
    if (!inherits(seurat_obj, "Seurat")) {
        stop(
            "Expected a Seurat object, got: ", class(seurat_obj)[1], "\n",
            "Provide a valid Seurat v3, v4, or v5 object."
        )
    }

    # Default assay
    if (is.null(assay)) {
        assay <- Seurat::DefaultAssay(seurat_obj)
    }

    # Validate assay exists
    available_assays <- Seurat::Assays(seurat_obj)
    if (!assay %in% available_assays) {
        stop(
            "Assay '", assay, "' not found in Seurat object.\n",
            "Available assays: ", paste(available_assays, collapse = ", ")
        )
    }

    # Validate x_slot
    x_slot <- match.arg(x_slot)

    # Check for BPCells (hard error)
    .check_bpcells(seurat_obj, assay)

    # =========================================================================
    # Import Python Modules
    # =========================================================================

    anndata <- reticulate::import("anndata", convert = FALSE)
    np <- reticulate::import("numpy", convert = FALSE)
    scipy_sparse <- reticulate::import("scipy.sparse", convert = FALSE)

    # =========================================================================
    # Extract Expression Data
    # =========================================================================

    # Get data for X
    if (x_slot == "data") {
        x_data_r <- .get_assay_data(seurat_obj, assay, "data")
        x_slot_label <- "data"
    } else {
        x_data_r <- .get_assay_data(seurat_obj, assay, "counts")
        x_slot_label <- "counts"
    }

    # Validate X has data
    if (!.has_data(x_data_r)) {
        stop(
            "No data found in '", x_slot, "' slot for assay '", assay, "'.\n",
            "Ensure the Seurat object has ", x_slot, " data, or use x_slot='",
            ifelse(x_slot == "data", "counts", "data"), "'."
        )
    }

    # Convert X to Python
    X_py <- .to_python_matrix(x_data_r, np, scipy_sparse)

    # =========================================================================
    # Build Layers
    # =========================================================================

    layers <- list()

    # Include counts in layers (if X is normalized and counts exist)
    if (include_counts && x_slot == "data") {
        counts_r <- .get_assay_data(seurat_obj, assay, "counts")
        if (.has_data(counts_r)) {
            layers[["counts"]] <- .to_python_matrix(counts_r, np, scipy_sparse)
        } else {
            warning(
                "Raw counts not available in assay '", assay, "'.\n",
                "layers['counts'] will not be created."
            )
        }
    }

    # Include scaled data in layers
    if (include_scale) {
        scale_r <- .get_assay_data(seurat_obj, assay, "scale.data")
        if (.has_data(scale_r)) {
            # Verify scale.data covers all features
            full_genes <- rownames(Seurat::GetAssay(seurat_obj, assay = assay))
            scale_genes <- rownames(scale_r)

            if (!identical(scale_genes, full_genes)) {
                warning(
                    "scale.data does not cover all features (", length(scale_genes),
                    " vs ", length(full_genes), ").\n",
                    "layers['scale'] will not be created to avoid shape mismatch.\n",
                    "Run ScaleData() on all features if you need scale.data in AnnData."
                )
            } else {
                layers[["scale"]] <- .to_python_matrix(scale_r, np, scipy_sparse)
            }
        } else {
            warning(
                "include_scale=TRUE but scale.data is empty in assay '", assay, "'.\n",
                "Run ScaleData() first if you need scaled data in AnnData."
            )
        }
    }

    # =========================================================================
    # Cell Metadata (obs)
    # =========================================================================

    obs_df <- seurat_obj[[]]
    cell_names <- colnames(seurat_obj)

    # Convert factors to character
    obs_df <- .factors_to_character(obs_df)

    # Ensure rownames are set
    if (is.null(rownames(obs_df)) || !identical(rownames(obs_df), cell_names)) {
        rownames(obs_df) <- cell_names
    }

    # =========================================================================
    # Gene Metadata (var)
    # =========================================================================

    assay_obj <- Seurat::GetAssay(seurat_obj, assay = assay)
    var_df <- assay_obj[[]]
    gene_names <- rownames(assay_obj)

    # Convert factors to character
    var_df <- .factors_to_character(var_df)

    # Ensure rownames are set
    if (is.null(rownames(var_df)) || !identical(rownames(var_df), gene_names)) {
        rownames(var_df) <- gene_names
    }

    # =========================================================================
    # Validate Dimensions
    # =========================================================================

    n_cells <- length(cell_names)
    n_genes <- length(gene_names)

    # Check X dimensions (Python format: cells x genes)
    x_shape <- reticulate::py_to_r(X_py$shape)
    if (x_shape[1] != n_cells || x_shape[2] != n_genes) {
        stop(
            "X matrix dimensions do not match obs/var.\n",
            "X shape: ", x_shape[1], " x ", x_shape[2], "\n",
            "Expected: ", n_cells, " cells x ", n_genes, " genes"
        )
    }

    # Check layers dimensions
    for (layer_name in names(layers)) {
        layer_shape <- reticulate::py_to_r(layers[[layer_name]]$shape)
        if (layer_shape[1] != n_cells || layer_shape[2] != n_genes) {
            stop(
                "Layer '", layer_name, "' dimensions do not match X.\n",
                "Layer shape: ", layer_shape[1], " x ", layer_shape[2], "\n",
                "Expected: ", n_cells, " x ", n_genes
            )
        }
    }

    # =========================================================================
    # Reductions (obsm) and Feature Loadings (varm)
    # =========================================================================

    obsm <- list()
    varm <- list()

    if (include_reductions && length(seurat_obj@reductions) > 0) {
        for (red_name in names(seurat_obj@reductions)) {
            red <- seurat_obj@reductions[[red_name]]
            embeddings <- red@cell.embeddings

            if (nrow(embeddings) == 0) {
                warning(
                    "Reduction '", red_name, "' has no cell embeddings. Skipping."
                )
                next
            }

            # Validate cell alignment
            .validate_names(
                rownames(embeddings),
                cell_names,
                paste0("Reduction '", red_name, "' cell names")
            )

            # Map to AnnData obsm key
            obsm_key <- switch(
                tolower(red_name),
                "pca" = "X_pca",
                "umap" = "X_umap",
                "tsne" = "X_tsne",
                "spatial" = "spatial",
                "sp" = "spatial",
                paste0("X_", red_name)
            )

            # Check for key collision
            if (obsm_key %in% names(obsm)) {
                stop(
                    "obsm key collision: '", obsm_key, "' already exists.\n",
                    "This can happen if you have both 'spatial' and 'SP' reductions.\n",
                    "Rename one of the reductions in your Seurat object."
                )
            }

            obsm[[obsm_key]] <- np$array(embeddings, dtype = "float64")

            # Feature loadings (varm)
            if (!is.null(red@feature.loadings) && nrow(red@feature.loadings) > 0) {
                loadings <- red@feature.loadings

                # Only include if loadings match all genes
                if (identical(rownames(loadings), gene_names)) {
                    varm_key <- switch(
                        tolower(red_name),
                        "pca" = "PCs",
                        paste0(red_name, "_loadings")
                    )

                    if (!varm_key %in% names(varm)) {
                        varm[[varm_key]] <- np$array(as.matrix(loadings), dtype = "float64")
                    }
                } else {
                    warning(
                        "Reduction '", red_name, "' feature loadings do not cover all genes (",
                        nrow(loadings), " vs ", length(gene_names), "). Skipping varm."
                    )
                }
            }
        }
    }

    # Handle explicit spatial_key override
    if (!is.null(spatial_key)) {
        if (!spatial_key %in% names(seurat_obj@reductions)) {
            stop(
                "Spatial reduction '", spatial_key, "' not found.\n",
                "Available reductions: ", paste(names(seurat_obj@reductions), collapse = ", ")
            )
        }
        red <- seurat_obj@reductions[[spatial_key]]
        obsm[["spatial"]] <- np$array(red@cell.embeddings, dtype = "float64")
    }

    # Auto-detect spatial if not already present
    if (!"spatial" %in% names(obsm)) {
        detected <- .detect_spatial_reduction(seurat_obj)
        if (!is.null(detected) && detected %in% names(seurat_obj@reductions)) {
            red <- seurat_obj@reductions[[detected]]
            obsm[["spatial"]] <- np$array(red@cell.embeddings, dtype = "float64")
        }
    }

    # =========================================================================
    # Graphs (obsp)
    # =========================================================================

    obsp <- list()

    if (include_graphs && length(seurat_obj@graphs) > 0) {
        graph_names <- names(seurat_obj@graphs)

        # Check for ambiguous graph mapping
        nn_matches <- grep("_nn$", graph_names, value = TRUE)
        snn_matches <- grep("_snn$", graph_names, value = TRUE)

        if ((length(nn_matches) > 1 || length(snn_matches) > 1) && is.null(graph_key_map)) {
            stop(
                "Multiple graphs match _nn or _snn pattern.\n",
                "  _nn graphs: ", paste(nn_matches, collapse = ", "), "\n",
                "  _snn graphs: ", paste(snn_matches, collapse = ", "), "\n\n",
                "Provide graph_key_map to specify which graphs to use:\n",
                "  graph_key_map = c('", nn_matches[1], "' = 'connectivities', '",
                snn_matches[1], "' = 'distances')"
            )
        }

        for (graph_name in graph_names) {
            graph <- seurat_obj@graphs[[graph_name]]

            # Determine obsp key
            if (!is.null(graph_key_map) && graph_name %in% names(graph_key_map)) {
                obsp_key <- graph_key_map[[graph_name]]
            } else if (graph_name %in% nn_matches && length(nn_matches) == 1) {
                obsp_key <- "connectivities"
            } else if (graph_name %in% snn_matches && length(snn_matches) == 1) {
                obsp_key <- "distances"
            } else {
                obsp_key <- graph_name
            }

            # Check for collision
            if (obsp_key %in% names(obsp)) {
                stop(
                    "obsp key collision: '", obsp_key, "' already exists.\n",
                    "Provide graph_key_map to disambiguate graph names."
                )
            }

            # Convert graph to sparse matrix
            if (inherits(graph, "Graph")) {
                graph_mat <- as(graph, "dgCMatrix")

                # Validate dimensions
                if (nrow(graph_mat) != n_cells || ncol(graph_mat) != n_cells) {
                    stop(
                        "Graph '", graph_name, "' dimensions (", nrow(graph_mat),
                        " x ", ncol(graph_mat), ") do not match cell count (", n_cells, ")."
                    )
                }

                # Convert graph directly to Python CSR (no transpose needed for square cell x cell)
                # Note: .to_python_matrix() is for genes x cells -> cells x genes conversion
                # For graphs, we convert directly to avoid double transpose
                if (length(graph_mat@x) == 0) {
                    shape <- reticulate::tuple(as.integer(n_cells), as.integer(n_cells))
                    obsp[[obsp_key]] <- scipy_sparse$csr_matrix(shape, dtype = "float64")
                } else {
                    mat_csr <- as(graph_mat, "RsparseMatrix")
                    data <- np$array(as.numeric(mat_csr@x), dtype = "float64")
                    indices <- np$array(as.integer(mat_csr@j), dtype = "int32")
                    indptr <- np$array(as.integer(mat_csr@p), dtype = "int32")
                    shape <- reticulate::tuple(as.integer(n_cells), as.integer(n_cells))
                    obsp[[obsp_key]] <- scipy_sparse$csr_matrix(
                        reticulate::tuple(data, indices, indptr),
                        shape = shape
                    )
                }
            } else {
                warning(
                    "Graph '", graph_name, "' is not a Seurat Graph object. Skipping."
                )
            }
        }
    }

    # =========================================================================
    # Unstructured Data (uns)
    # =========================================================================

    uns <- list(
        "X_slot" = x_slot_label,
        "seurat_assay" = assay,
        "seurat_project" = seurat_obj@project.name,
        "conversion" = list(
            "source" = "seurat_to_adata",
            "seurat_version" = as.character(packageVersion("Seurat")),
            "r_version" = as.character(getRversion()),
            "timestamp" = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        )
    )

    # Store graph_key_map for round-trip
    if (!is.null(graph_key_map)) {
        uns[["graph_key_map"]] <- as.list(graph_key_map)
    }

    # Store misc if present
    if (length(seurat_obj@misc) > 0) {
        uns[["seurat_misc"]] <- seurat_obj@misc
    }

    # =========================================================================
    # Create AnnData Object
    # =========================================================================

    # Convert R lists to Python dicts
    layers_py <- if (length(layers) > 0) reticulate::r_to_py(layers) else NULL
    obsm_py <- if (length(obsm) > 0) reticulate::r_to_py(obsm) else NULL
    obsp_py <- if (length(obsp) > 0) reticulate::r_to_py(obsp) else NULL
    varm_py <- if (length(varm) > 0) reticulate::r_to_py(varm) else NULL
    uns_py <- reticulate::r_to_py(uns)

    adata <- anndata$AnnData(
        X = X_py,
        obs = obs_df,
        var = var_df,
        layers = layers_py,
        obsm = obsm_py,
        obsp = obsp_py,
        varm = varm_py,
        uns = uns_py
    )

    # =========================================================================
    # Summary Message
    # =========================================================================

    message(
        "Converted Seurat to AnnData:\n",
        "  Cells: ", n_cells, "\n",
        "  Genes: ", n_genes, "\n",
        "  X contains: ", x_slot_label, "\n",
        "  Layers: ", ifelse(length(layers) > 0, paste(names(layers), collapse = ", "), "none"), "\n",
        "  Reductions: ", ifelse(length(obsm) > 0, paste(names(obsm), collapse = ", "), "none"), "\n",
        "  Graphs: ", ifelse(length(obsp) > 0, paste(names(obsp), collapse = ", "), "none")
    )

    return(adata)
}


# ==============================================================================
# adata_to_seurat: AnnData -> Seurat
# ==============================================================================

#' Convert AnnData to Seurat Object
#'
#' Converts a Python AnnData object to a Seurat object. Can either create a
#' new Seurat object or merge results into an existing one.
#'
#' @param adata An AnnData object (Python object via reticulate)
#' @param seurat_obj Optional. Existing Seurat object to merge into.
#'   If NULL, creates a new Seurat object.
#' @param assay Character. Name for the assay in Seurat. Default: "RNA"
#' @param project Character. Project name for new Seurat object. Default: "SpatialCore"
#'
#' @return A Seurat object
#'
#' @details
#' ## X Slot Handling
#' The function reads `adata.uns["X_slot"]` to determine how X was stored:
#' - "data": X contains normalized values -> goes to `data` slot
#' - "counts": X contains raw counts -> goes to `counts` slot
#'
#' If `uns["X_slot"]` is missing, X is assumed to be normalized data.
#'
#' ## Merge Behavior
#' When `seurat_obj` is provided:
#' - New columns in `adata.obs` are ADDED to Seurat metadata
#' - Existing columns are UPDATED with `adata.obs` values
#' - New reductions are ADDED
#' - Existing reductions are REPLACED
#' - Cell and gene names must match exactly
#'
#' ## Reduction Key Conventions
#' Seurat expects specific key prefixes for reductions:
#' - PCA: `PC_1`, `PC_2`, ...
#' - UMAP: `UMAP_1`, `UMAP_2`, ...
#' - t-SNE: `tSNE_1`, `tSNE_2`, ...
#' - Spatial: `Spatial_1`, `Spatial_2`, ...
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Create new Seurat object from AnnData
#' seurat_obj <- adata_to_seurat(adata)
#'
#' # Merge results back into existing Seurat object
#' seurat_obj <- adata_to_seurat(adata, seurat_obj = original_seurat)
#'
#' # Specify assay and project names
#' seurat_obj <- adata_to_seurat(adata, assay = "Spatial", project = "MyProject")
#' }
adata_to_seurat <- function(
    adata,
    seurat_obj = NULL,
    assay = "RNA",
    project = "SpatialCore"
) {

    # =========================================================================
    # Input Validation
    # =========================================================================

    # Check reticulate
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop("Package 'reticulate' is required but not installed.")
    }

    # Check Python
    if (!reticulate::py_available()) {
        stop("Python is not available. Configure reticulate first.")
    }

    # Check it's a Python object
    if (!inherits(adata, "python.builtin.object")) {
        stop(
            "Expected a Python AnnData object, got R class: ", class(adata)[1], "\n",
            "Use reticulate to import AnnData objects from Python."
        )
    }

    # Basic AnnData check
    adata_class <- class(adata)[1]
    if (!grepl("AnnData", adata_class, ignore.case = TRUE)) {
        stop(
            "Expected an AnnData object, got: ", adata_class, "\n",
            "Ensure you're passing an anndata.AnnData object."
        )
    }

    # =========================================================================
    # Extract Basic Info
    # =========================================================================

    # Get obs/var names
    obs_names <- reticulate::py_to_r(adata$obs_names$tolist())
    var_names <- reticulate::py_to_r(adata$var_names$tolist())
    n_cells <- length(obs_names)
    n_genes <- length(var_names)

    # Get X slot semantics - REQUIRED for safe conversion
    x_slot <- NULL
    if (!is.null(adata$uns)) {
        x_slot_val <- reticulate::py_to_r(adata$uns$get("X_slot"))
        if (!is.null(x_slot_val)) {
            x_slot <- x_slot_val
        }
    }

    if (is.null(x_slot)) {
        stop(
            "adata.uns['X_slot'] not found.\n",
            "Cannot safely determine if X contains 'data' (normalized) or 'counts' (raw).\n",
            "Set adata.uns['X_slot'] = 'data' or 'counts' before conversion."
        )
    }

    if (!x_slot %in% c("data", "counts")) {
        stop(
            "Invalid adata.uns['X_slot'] value: '", x_slot, "'.\n",
            "Expected 'data' or 'counts'."
        )
    }

    # =========================================================================
    # Validate X Matrix
    # =========================================================================

    if (is.null(adata$X)) {
        stop("adata.X is NULL or empty. Cannot create Seurat object without expression data.")
    }

    x_shape <- reticulate::py_to_r(adata$X$shape)
    if (x_shape[1] != n_cells || x_shape[2] != n_genes) {
        stop(
            "adata.X shape (", x_shape[1], " x ", x_shape[2], ") does not match ",
            "obs (", n_cells, ") and var (", n_genes, ") dimensions."
        )
    }

    # =========================================================================
    # Extract Expression Matrices
    # =========================================================================

    # X matrix
    X_r <- .from_python_matrix(adata$X, var_names, obs_names)

    # Layers - check if layers exists and has keys
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

    # Assign X based on x_slot
    if (x_slot == "counts") {
        if (is.null(counts_r)) {
            counts_r <- X_r
        }
    } else {
        if (is.null(data_r)) {
            data_r <- X_r
        }
    }

    # Ensure we have counts for CreateSeuratObject - REQUIRED
    if (is.null(counts_r)) {
        stop(
            "No raw counts found for Seurat object creation.\n",
            "  - adata.uns['X_slot'] = '", x_slot, "' (so X is not counts)\n",
            "  - adata.layers['counts'] is missing\n\n",
            "To fix, either:\n",
            "  1. Set adata.layers['counts'] = <raw count matrix>\n",
            "  2. Or if X contains counts, set adata.uns['X_slot'] = 'counts'"
        )
    }

    # =========================================================================
    # Cell Metadata (obs)
    # =========================================================================

    obs_df <- reticulate::py_to_r(adata$obs)

    # Ensure rownames
    if (is.null(rownames(obs_df))) {
        rownames(obs_df) <- obs_names
    } else {
        .validate_names(rownames(obs_df), obs_names, "obs rownames")
    }

    # =========================================================================
    # Gene Metadata (var)
    # =========================================================================

    var_df <- reticulate::py_to_r(adata$var)

    # Ensure rownames
    if (is.null(rownames(var_df))) {
        rownames(var_df) <- var_names
    } else {
        .validate_names(rownames(var_df), var_names, "var rownames")
    }

    # =========================================================================
    # Create or Update Seurat Object
    # =========================================================================

    if (is.null(seurat_obj)) {
        # ----- Create New Seurat Object -----

        seurat_obj <- Seurat::CreateSeuratObject(
            counts = counts_r,
            project = project,
            assay = assay,
            meta.data = obs_df
        )

        # Add normalized data if available
        if (!is.null(data_r)) {
            seurat_obj <- .set_assay_data(seurat_obj, assay, "data", data_r)
        }

        # Add scaled data if available
        if (!is.null(scale_r)) {
            # Verify alignment
            if (identical(rownames(scale_r), var_names) && identical(colnames(scale_r), obs_names)) {
                seurat_obj <- .set_assay_data(seurat_obj, assay, "scale.data", scale_r)
            } else {
                warning("layers['scale'] does not align with features/cells. Skipping.")
            }
        }

        # Add gene metadata
        if (ncol(var_df) > 0) {
            assay_obj <- Seurat::GetAssay(seurat_obj, assay = assay)
            # Only add if rownames match
            if (identical(rownames(var_df), rownames(assay_obj))) {
                assay_obj[[names(var_df)]] <- var_df
                seurat_obj[[assay]] <- assay_obj
            } else {
                warning(
                    "adata.var rownames do not match assay gene names.\n",
                    "Gene metadata (var) will not be added to Seurat object."
                )
            }
        }

    } else {
        # ----- Merge into Existing Seurat Object -----

        if (!inherits(seurat_obj, "Seurat")) {
            stop("seurat_obj must be a Seurat object.")
        }

        # Validate cell alignment
        seurat_cells <- colnames(seurat_obj)
        .validate_names(obs_names, seurat_cells, "Cell names (AnnData vs Seurat)")

        # Validate gene alignment
        seurat_genes <- rownames(Seurat::GetAssay(seurat_obj, assay = assay))
        .validate_names(var_names, seurat_genes, "Gene names (AnnData vs Seurat)")

        # Merge cell metadata (add/update columns)
        for (col in colnames(obs_df)) {
            seurat_obj[[col]] <- obs_df[[col]]
        }

        # Update expression data if provided
        if (!is.null(data_r)) {
            seurat_obj <- .set_assay_data(seurat_obj, assay, "data", data_r)
        }
        if (!is.null(counts_r) && x_slot == "counts") {
            seurat_obj <- .set_assay_data(seurat_obj, assay, "counts", counts_r)
        }
        if (!is.null(scale_r)) {
            if (identical(rownames(scale_r), seurat_genes) && identical(colnames(scale_r), seurat_cells)) {
                seurat_obj <- .set_assay_data(seurat_obj, assay, "scale.data", scale_r)
            } else {
                warning(
                    "layers['scale'] dimensions do not match existing Seurat object.\n",
                    "scale.data will not be updated during merge."
                )
            }
        }
    }

    # =========================================================================
    # Add Reductions (obsm)
    # =========================================================================

    # Get obsm keys safely
    obsm_keys <- .py_keys_to_r(adata$obsm)

    if (!is.null(obsm_keys) && length(obsm_keys) > 0) {
        for (obsm_key in obsm_keys) {
            embeddings <- reticulate::py_to_r(adata$obsm[[obsm_key]])

            # Validate dimensions
            if (nrow(embeddings) != n_cells) {
                warning(
                    "obsm['", obsm_key, "'] has ", nrow(embeddings), " rows but ",
                    n_cells, " cells expected. Skipping."
                )
                next
            }

            # Map to Seurat reduction name
            red_name <- switch(
                obsm_key,
                "X_pca" = "pca",
                "X_umap" = "umap",
                "X_tsne" = "tsne",
                "spatial" = "spatial",
                gsub("^X_", "", obsm_key)
            )

            # Set key convention
            key <- switch(
                red_name,
                "pca" = "PC_",
                "umap" = "UMAP_",
                "tsne" = "tSNE_",
                "spatial" = "Spatial_",
                paste0(toupper(substr(red_name, 1, 1)), tolower(substr(red_name, 2, nchar(red_name))), "_")
            )

            # Set rownames and colnames
            rownames(embeddings) <- colnames(seurat_obj)
            colnames(embeddings) <- paste0(key, seq_len(ncol(embeddings)))

            # Create and add reduction
            seurat_obj[[red_name]] <- Seurat::CreateDimReducObject(
                embeddings = embeddings,
                key = key,
                assay = assay
            )
        }
    }

    # =========================================================================
    # Add Feature Loadings (varm)
    # =========================================================================

    # Get varm keys safely
    varm_keys <- .py_keys_to_r(adata$varm)

    if (!is.null(varm_keys) && length(varm_keys) > 0) {
        for (varm_key in varm_keys) {
            loadings <- reticulate::py_to_r(adata$varm[[varm_key]])

            # Map to reduction name
            red_name <- switch(
                varm_key,
                "PCs" = "pca",
                gsub("_loadings$", "", varm_key)
            )

            # Only add if reduction exists
            if (!red_name %in% names(seurat_obj@reductions)) {
                warning(
                    "varm['", varm_key, "'] maps to reduction '", red_name,
                    "' which does not exist. Skipping feature loadings."
                )
                next
            }

            # Validate dimensions
            if (nrow(loadings) != n_genes) {
                warning(
                    "varm['", varm_key, "'] has ", nrow(loadings), " rows but ",
                    n_genes, " genes expected. Skipping."
                )
                next
            }

            # Set rownames
            rownames(loadings) <- rownames(Seurat::GetAssay(seurat_obj, assay = assay))

            # Add to reduction
            seurat_obj[[red_name]]@feature.loadings <- loadings
        }
    }

    # =========================================================================
    # Add Graphs (obsp)
    # =========================================================================

    # Get obsp keys safely
    obsp_keys <- .py_keys_to_r(adata$obsp)

    if (!is.null(obsp_keys) && length(obsp_keys) > 0) {
        for (obsp_key in obsp_keys) {
            graph_py <- adata$obsp[[obsp_key]]

            # Validate dimensions
            graph_shape <- reticulate::py_to_r(graph_py$shape)
            if (graph_shape[1] != n_cells || graph_shape[2] != n_cells) {
                warning(
                    "obsp['", obsp_key, "'] shape (", graph_shape[1], " x ", graph_shape[2],
                    ") is not n_cells x n_cells. Skipping."
                )
                next
            }

            # Convert graph from Python to R (cells x cells, no transpose needed)
            # Note: .from_python_matrix() is for cells x genes -> genes x cells
            # For square graphs, we convert directly without transpose
            graph_csc <- graph_py$tocsc()
            data <- as.numeric(reticulate::py_to_r(graph_csc$data))
            indices <- as.integer(reticulate::py_to_r(graph_csc$indices))
            indptr <- as.integer(reticulate::py_to_r(graph_csc$indptr))

            graph_r <- Matrix::sparseMatrix(
                i = indices + 1L,
                p = indptr,
                x = data,
                dims = c(n_cells, n_cells),
                dimnames = list(obs_names, obs_names)
            )
            # Ensure dgCMatrix format for Seurat
            graph_r <- as(graph_r, "dgCMatrix")

            # Map to Seurat graph name
            graph_name <- switch(
                obsp_key,
                "connectivities" = paste0(assay, "_nn"),
                "distances" = paste0(assay, "_snn"),
                obsp_key
            )

            # Convert to Seurat Graph and add
            seurat_obj@graphs[[graph_name]] <- as(graph_r, "Graph")
        }
    }

    # =========================================================================
    # Summary Message
    # =========================================================================

    message(
        "Converted AnnData to Seurat:\n",
        "  Cells: ", ncol(seurat_obj), "\n",
        "  Genes: ", nrow(seurat_obj), "\n",
        "  Metadata columns: ", ncol(seurat_obj[[]]), "\n",
        "  Reductions: ", paste(names(seurat_obj@reductions), collapse = ", "), "\n",
        "  Graphs: ", paste(names(seurat_obj@graphs), collapse = ", ")
    )

    return(seurat_obj)
}
