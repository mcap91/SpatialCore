# ==============================================================================
# setup.R - SpatialCore R Environment Setup
# ==============================================================================
#
# Functions for validating and configuring the R/Python environment
# for SpatialCore usage via reticulate.
#
# ==============================================================================

#' Check and Setup SpatialCore Environment
#'
#' Validates that reticulate, Python, and required packages are available.
#' Call this once per R session before using conversion functions.
#'
#' @param conda_env Character. Name of conda environment. Default: "spatialcore"
#' @param verbose Logical. Print status messages. Default: TRUE
#'
#' @return Invisibly returns TRUE if setup successful.
#' @export
#'
#' @examples
#' \dontrun{
#' setup_spatialcore()
#' setup_spatialcore(conda_env = "my_env", verbose = FALSE)
#' }
setup_spatialcore <- function(conda_env = "spatialcore", verbose = TRUE) {

    # -------------------------------------------------------------------------
    # Check reticulate is installed
    # -------------------------------------------------------------------------
    if (!requireNamespace("reticulate", quietly = TRUE)) {
        stop(
            "Package 'reticulate' is required but not installed.\n",
            "Install with: install.packages('reticulate')"
        )
    }

    # -------------------------------------------------------------------------
    # Check conda environment exists
    # -------------------------------------------------------------------------
    envs <- tryCatch(
        reticulate::conda_list(),
        error = function(e) {
            stop(
                "Could not list conda environments.\n",
                "Ensure conda/mamba is installed and accessible.\n",
                "Original error: ", conditionMessage(e)
            )
        }
    )

    if (!conda_env %in% envs$name) {
        stop(
            "Conda environment '", conda_env, "' not found.\n",
            "Available environments: ", paste(envs$name, collapse = ", "), "\n\n",
            "Create the environment with:\n",
            "  mamba create -n ", conda_env, " python=3.10\n",
            "  mamba activate ", conda_env, "\n",
            "  pip install spatialcore scanpy anndata"
        )
    }

    # -------------------------------------------------------------------------
    # Activate the conda environment
    # -------------------------------------------------------------------------
    tryCatch(
        reticulate::use_condaenv(conda_env, required = TRUE),
        error = function(e) {
            stop(
                "Could not activate conda environment '", conda_env, "'.\n",
                "Original error: ", conditionMessage(e)
            )
        }
    )

    # -------------------------------------------------------------------------
    # Verify Python is available
    # -------------------------------------------------------------------------
    if (!reticulate::py_available(initialize = TRUE)) {
        stop(
            "Python is not available after activating environment '", conda_env, "'.\n",
            "Check your conda/reticulate configuration."
        )
    }

    # -------------------------------------------------------------------------
    # Check required Python packages
    # -------------------------------------------------------------------------
    required_packages <- c("anndata", "numpy", "scipy", "pandas")
    optional_packages <- c("scanpy", "spatialcore")

    missing_required <- character(0)
    missing_optional <- character(0)

    for (pkg in required_packages) {
        if (!reticulate::py_module_available(pkg)) {
            missing_required <- c(missing_required, pkg)
        }
    }

    for (pkg in optional_packages) {
        if (!reticulate::py_module_available(pkg)) {
            missing_optional <- c(missing_optional, pkg)
        }
    }

    if (length(missing_required) > 0) {
        stop(
            "Required Python packages not found: ", paste(missing_required, collapse = ", "), "\n",
            "Install with: pip install ", paste(missing_required, collapse = " ")
        )
    }

    # Always warn about missing optional packages (regardless of verbose)
    if (length(missing_optional) > 0) {
        warning(
            "Optional Python packages not found: ", paste(missing_optional, collapse = ", "), "\n",
            "Some functionality may be limited.\n",
            "Install with: pip install ", paste(missing_optional, collapse = " ")
        )
    }

    # -------------------------------------------------------------------------
    # Print status if verbose
    # -------------------------------------------------------------------------
    if (verbose) {
        py_config <- reticulate::py_config()
        message("SpatialCore environment ready:")
        message("  Python: ", py_config$python)
        message("  Version: ", py_config$version)

        # Get package versions
        anndata <- reticulate::import("anndata")
        message("  anndata: ", anndata$`__version__`)

        if (reticulate::py_module_available("spatialcore")) {
            sc <- reticulate::import("spatialcore")
            message("  spatialcore: ", sc$`__version__`)
        }
    }

    invisible(TRUE)
}


#' Check if SpatialCore Environment is Ready
#'
#' Quick check without activation. Returns TRUE/FALSE.
#'
#' @param conda_env Character. Name of conda environment. Default: "spatialcore"
#'
#' @return Logical. TRUE if environment is ready, FALSE otherwise.
#' @export
#'
#' @examples
#' \dontrun{
#' if (!is_spatialcore_ready()) {
#'     setup_spatialcore()
#' }
#' }
is_spatialcore_ready <- function(conda_env = "spatialcore") {

    if (!requireNamespace("reticulate", quietly = TRUE)) {
        return(FALSE)
    }

    if (!reticulate::py_available()) {
        return(FALSE)
    }

    # Check if anndata is available (minimum requirement)
    if (!reticulate::py_module_available("anndata")) {
        return(FALSE)
    }

    return(TRUE)
}


#' Get SpatialCore Python Module
#'
#' Import and return the spatialcore Python module.
#' Raises error if not available.
#'
#' @return Python module object for spatialcore
#' @export
#'
#' @examples
#' \dontrun{
#' sc <- get_spatialcore()
#' sc$spatial$compute_morans_i(adata)
#' }
get_spatialcore <- function() {

    if (!reticulate::py_module_available("spatialcore")) {
        stop(
            "Python package 'spatialcore' not found.\n",
            "Install with: pip install spatialcore"
        )
    }

    reticulate::import("spatialcore")
}
