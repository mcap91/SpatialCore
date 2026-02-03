# ==============================================================================
# convert.R - Seurat <-> AnnData Conversion Functions
# ==============================================================================
#
# Bidirectional conversion between Seurat (R) and AnnData (Python) objects.
# Enables R users to use any SpatialCore function via reticulate.
#
# Dependencies:
#   - helpers.R (internal functions)
#   - anndataR package (Bioconductor)
#   - reticulate, Matrix, Seurat packages
#
# ==============================================================================

# Source helpers if not already loaded
# (In package context, this is handled by NAMESPACE/Collate)


# ==============================================================================
# seurat_to_adata: Seurat -> AnnData (using anndataR)
# ==============================================================================

#' Convert Seurat Object to AnnData
#'
#' Converts a Seurat object to a Python AnnData object in memory using anndataR.
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
#' This function supports Seurat v3, v4, and v5 objects via anndataR.
#' BPCells on-disk matrices (Seurat v5) are NOT supported and will raise an error.
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

    # Check for BPCells (hard error - keep our helpful message)
    .check_bpcells(seurat_obj, assay)

    # =========================================================================
    # Build layers_mapping for anndataR
    # =========================================================================

    layers_mapping <- NULL
    if (include_counts && x_slot == "data") {
        # Include counts in layers when X contains normalized data
        layers_mapping <- c("counts")
    }
    if (include_scale) {
        layers_mapping <- c(layers_mapping, "scale.data")
    }

    # =========================================================================
    # Core Conversion via anndataR
    # =========================================================================

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

    # =========================================================================
    # Post-processing: Add SpatialCore-specific uns fields
    # =========================================================================

    # Set uns["X_slot"] for round-trip safety (REQUIRED for adata_to_seurat)
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

    # Store graph_key_map for round-trip
    if (!is.null(graph_key_map)) {
        adata$uns[["graph_key_map"]] <- as.list(graph_key_map)
    }

    # Store misc if present
    if (length(seurat_obj@misc) > 0) {
        adata$uns[["seurat_misc"]] <- seurat_obj@misc
    }

    # =========================================================================
    # Post-processing: Handle spatial_key override
    # =========================================================================

    if (!is.null(spatial_key)) {
        if (!spatial_key %in% names(seurat_obj@reductions)) {
            stop(
                "Spatial reduction '", spatial_key, "' not found.\n",
                "Available reductions: ", paste(names(seurat_obj@reductions), collapse = ", ")
            )
        }
        # Override spatial coordinates
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

    # =========================================================================
    # Post-processing: Handle graph_key_map renaming
    # =========================================================================

    if (!is.null(graph_key_map) && include_graphs) {
        obsp_keys <- names(adata$obsp)
        for (old_key in names(graph_key_map)) {
            new_key <- graph_key_map[[old_key]]
            if (old_key %in% obsp_keys && old_key != new_key) {
                # Rename the graph key
                adata$obsp[[new_key]] <- adata$obsp[[old_key]]
                adata$obsp[[old_key]] <- NULL
            }
        }
    }

    # =========================================================================
    # Summary Message
    # =========================================================================

    # Get dimensions - anndataR uses methods, not properties
    n_cells <- nrow(adata$obs)
    n_genes <- nrow(adata$var)
    layer_keys <- names(adata$layers)
    obsm_keys <- names(adata$obsm)
    obsp_keys <- names(adata$obsp)

    message(
        "Converted Seurat to AnnData (via anndataR):\n",
        "  Cells: ", n_cells, "\n",
        "  Genes: ", n_genes, "\n",
        "  X contains: ", x_slot, "\n",
        "  Layers: ", ifelse(length(layer_keys) > 0, paste(layer_keys, collapse = ", "), "none"), "\n",
        "  Reductions: ", ifelse(length(obsm_keys) > 0, paste(obsm_keys, collapse = ", "), "none"), "\n",
        "  Graphs: ", ifelse(length(obsp_keys) > 0, paste(obsp_keys, collapse = ", "), "none")
    )

    return(adata)
}


# ==============================================================================
# adata_to_seurat: AnnData -> Seurat
# ==============================================================================

#' Convert AnnData to Seurat Object
#'
#' Converts an AnnData object (R-based from anndataR or Python-based via reticulate)
#' to a Seurat object. Can either create a new Seurat object or merge results
#' into an existing one.
#'
#' @param adata An AnnData object. Can be:
#'   - R-based AnnData from anndataR (output of seurat_to_adata)
#'   - Python AnnData object via reticulate
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
#' If `uns["X_slot"]` is missing, conversion fails with an error.
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
    # Detect AnnData Type (R-based anndataR or Python-based)
    # =========================================================================

    is_r_anndata <- inherits(adata, "AnnData") || inherits(adata, "AbstractAnnData")
    is_py_anndata <- inherits(adata, "python.builtin.object")

    if (!is_r_anndata && !is_py_anndata) {
        stop(
            "Expected an AnnData object, got R class: ", class(adata)[1], "\n",
            "Provide an anndataR AnnData or Python anndata.AnnData object."
        )
    }

    # =========================================================================
    # Extract Data Based on AnnData Type
    # =========================================================================

    if (is_r_anndata) {
        # ----- R-based AnnData (from anndataR) -----

        obs_names <- rownames(adata$obs)
        var_names <- rownames(adata$var)
        n_cells <- length(obs_names)
        n_genes <- length(var_names)

        # Get X slot semantics
        x_slot <- adata$uns[["X_slot"]]
        if (is.null(x_slot)) {
            stop(
                "adata$uns['X_slot'] not found.\n",
                "Cannot safely determine if X contains 'data' (normalized) or 'counts' (raw).\n",
                "Set adata$uns['X_slot'] <- 'data' or 'counts' before conversion."
            )
        }
        if (!x_slot %in% c("data", "counts")) {
            stop("Invalid adata$uns['X_slot'] value: '", x_slot, "'. Expected 'data' or 'counts'.")
        }

        # Extract X matrix (already in genes x cells format for R)
        X_r <- adata$X
        if (is.null(X_r)) {
            stop("adata$X is NULL or empty. Cannot create Seurat object without expression data.")
        }
        # anndataR stores as cells x genes, transpose to genes x cells
        X_r <- Matrix::t(X_r)
        rownames(X_r) <- var_names
        colnames(X_r) <- obs_names

        # Extract layers
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

        # Cell and gene metadata
        obs_df <- adata$obs
        var_df <- adata$var

        # Reductions (obsm)
        obsm_keys <- names(adata$obsm)
        obsm_data <- adata$obsm

        # Feature loadings (varm)
        varm_keys <- names(adata$varm)
        varm_data <- adata$varm

        # Graphs (obsp)
        obsp_keys <- names(adata$obsp)
        obsp_data <- adata$obsp

    } else {
        # ----- Python-based AnnData (via reticulate) -----

        if (!requireNamespace("reticulate", quietly = TRUE)) {
            stop("Package 'reticulate' is required for Python AnnData objects.")
        }
        if (!reticulate::py_available()) {
            stop("Python is not available. Configure reticulate first.")
        }

        # Basic AnnData check
        adata_class <- class(adata)[1]
        if (!grepl("AnnData", adata_class, ignore.case = TRUE)) {
            stop("Expected an AnnData object, got: ", adata_class)
        }

        obs_names <- reticulate::py_to_r(adata$obs_names$tolist())
        var_names <- reticulate::py_to_r(adata$var_names$tolist())
        n_cells <- length(obs_names)
        n_genes <- length(var_names)

        # Get X slot semantics
        x_slot <- NULL
        if (!is.null(adata$uns)) {
            x_slot_val <- reticulate::py_to_r(adata$uns$get("X_slot"))
            if (!is.null(x_slot_val)) x_slot <- x_slot_val
        }
        if (is.null(x_slot)) {
            stop(
                "adata.uns['X_slot'] not found.\n",
                "Cannot safely determine if X contains 'data' (normalized) or 'counts' (raw).\n",
                "Set adata.uns['X_slot'] = 'data' or 'counts' before conversion."
            )
        }
        if (!x_slot %in% c("data", "counts")) {
            stop("Invalid adata.uns['X_slot'] value: '", x_slot, "'. Expected 'data' or 'counts'.")
        }

        # Validate X
        if (is.null(adata$X)) {
            stop("adata.X is NULL or empty. Cannot create Seurat object without expression data.")
        }
        x_shape <- reticulate::py_to_r(adata$X$shape)
        if (x_shape[1] != n_cells || x_shape[2] != n_genes) {
            stop("adata.X shape (", x_shape[1], " x ", x_shape[2], ") does not match obs/var dimensions.")
        }

        # Extract X matrix
        X_r <- .from_python_matrix(adata$X, var_names, obs_names)

        # Extract layers
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

        # Cell and gene metadata
        obs_df <- reticulate::py_to_r(adata$obs)
        if (is.null(rownames(obs_df))) rownames(obs_df) <- obs_names
        var_df <- reticulate::py_to_r(adata$var)
        if (is.null(rownames(var_df))) rownames(var_df) <- var_names

        # Store keys and data references for later
        obsm_keys <- .py_keys_to_r(adata$obsm)
        obsm_data <- adata$obsm
        varm_keys <- .py_keys_to_r(adata$varm)
        varm_data <- adata$varm
        obsp_keys <- .py_keys_to_r(adata$obsp)
        obsp_data <- adata$obsp
    }

    # =========================================================================
    # Assign X based on x_slot
    # =========================================================================

    if (x_slot == "counts") {
        if (is.null(counts_r)) counts_r <- X_r
    } else {
        if (is.null(data_r)) data_r <- X_r
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

    # Ensure rownames on metadata
    if (is.null(rownames(obs_df))) rownames(obs_df) <- obs_names
    if (is.null(rownames(var_df))) rownames(var_df) <- var_names

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

    if (!is.null(obsm_keys) && length(obsm_keys) > 0) {
        for (obsm_key in obsm_keys) {
            # Get embeddings based on AnnData type
            if (is_r_anndata) {
                embeddings <- obsm_data[[obsm_key]]
            } else {
                embeddings <- reticulate::py_to_r(obsm_data[[obsm_key]])
            }

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

    if (!is.null(varm_keys) && length(varm_keys) > 0) {
        for (varm_key in varm_keys) {
            # Get loadings based on AnnData type
            if (is_r_anndata) {
                loadings <- varm_data[[varm_key]]
            } else {
                loadings <- reticulate::py_to_r(varm_data[[varm_key]])
            }

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

    if (!is.null(obsp_keys) && length(obsp_keys) > 0) {
        for (obsp_key in obsp_keys) {
            # Get graph based on AnnData type
            if (is_r_anndata) {
                graph_r <- obsp_data[[obsp_key]]
                # Validate dimensions
                if (nrow(graph_r) != n_cells || ncol(graph_r) != n_cells) {
                    warning(
                        "obsp['", obsp_key, "'] shape (", nrow(graph_r), " x ", ncol(graph_r),
                        ") is not n_cells x n_cells. Skipping."
                    )
                    next
                }
                # Ensure proper format
                if (!inherits(graph_r, "dgCMatrix")) {
                    graph_r <- as(graph_r, "dgCMatrix")
                }
                rownames(graph_r) <- obs_names
                colnames(graph_r) <- obs_names
            } else {
                graph_py <- obsp_data[[obsp_key]]

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
                graph_r <- as(graph_r, "dgCMatrix")
            }

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
