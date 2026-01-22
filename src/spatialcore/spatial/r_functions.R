# Spatial domain creation using sf package
# Called via spatialcore.r_bridge.run_r_code()
#
# Based on original implementation from docs/domains_r_code_old.md
# Modified for Python integration per docs/DOMAINS.md
#
# Required packages: sf, dplyr, purrr, concaveman, jsonlite
# Install: install.packages(c("sf", "concaveman", "dplyr", "purrr", "jsonlite"))

suppressPackageStartupMessages({
    library(sf)
    library(dplyr)
    library(purrr)
    library(concaveman)
    library(jsonlite)
})

if (Sys.info()[["sysname"]] == "Windows") {
    sf::sf_use_s2(FALSE)
    message("Windows detected: disabling s2 for sf operations")
}

#' Internal: Create spatial domains using Buffer-Union-Shrink algorithm
#'
#' @param input_csv Path to CSV with columns: cell, x, y, <group>
#' @param group Column name containing cell groupings
#' @param group_subset Value to filter for target cells
#' @param cell_dist Buffer distance in coordinate units (default 225)
#' @param shrink_margin Margin to keep when shrinking (default 25). Shrink = cell_dist - shrink_margin
#' @param domain_prefix Prefix for domain names
#' @param assign_all_cells If TRUE, assign ALL cells to domains (for heterogeneity analysis)
#' @return List with cell_data (cell-domain mapping), polygon_data (domain geometries), and original_data
#' @keywords internal
.MakeDomains <- function(input_csv, group, group_subset,
                        cell_dist = 225, shrink_margin = 25,
                        domain_prefix = "domain", assign_all_cells = TRUE) {

    if (is.null(group_subset)) {
        stop("Please provide a factor level or subset from the cell grouping 'group' input.")
    }

    # Read input data from CSV
    # check.names = FALSE preserves Python column names like "_filter"
    object <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)

    # Validate required columns
    required_cols <- c("cell", "x", "y", group)
    missing_cols <- setdiff(required_cols, names(object))
    if (length(missing_cols) > 0) {
        stop(paste("Missing columns:", paste(missing_cols, collapse = ", ")))
    }

    # Use explicit x, y column names
    x <- "x"
    y <- "y"

    # Buffer-Union-Shrink pipeline (preserving original logic)
    shrunken_polygons <- object %>%
        select(all_of(c(x, y, group))) %>%
        filter(.[[group]] == group_subset) %>%
        st_as_sf(coords = c(x, y)) %>%
        st_buffer(dist = cell_dist) %>%
        st_union() %>%
        st_cast("POLYGON") %>%
        st_buffer(-(cell_dist - shrink_margin)) %>%  # Now configurable
        st_cast("MULTIPOINT") %>%
        st_coordinates() %>%
        as.data.frame() %>%
        st_as_sf(coords = c("X", "Y"))

    if (nrow(shrunken_polygons) == 0) {
        stop("No polygons created. Try larger cell_dist or check input data.")
    }

    # Create concave hulls for each polygon (L1 identifies polygon ID)
    poly_list <- list()
    for (i in unique(shrunken_polygons$L1)) {
        d <- unique(shrunken_polygons$L1)[i]
        tmp <- shrunken_polygons %>% filter(L1 == d)
        conc <- concaveman::concaveman(tmp, length_threshold = 0, concavity = 0.999999)
        conc$domain <- d
        poly_list[[d]] <- conc
    }
    polygon_data <- do.call(rbind, poly_list)

    # Assign cells to polygons - PRESERVE GROUP COLUMN for target cell counting
    if (assign_all_cells) {
        # Assign ALL cells (for heterogeneity analysis) but keep group column
        cells_use <- object %>%
            select(all_of(c(x, y)), cell, all_of(group)) %>%
            st_as_sf(coords = c(x, y))
    } else {
        # Only target cells (original behavior)
        cells_use <- object %>%
            filter(.data[[group]] == group_subset) %>%
            select(all_of(c(x, y)), cell, all_of(group)) %>%
            st_as_sf(coords = c(x, y))
    }

    # Spatial join with largest = TRUE to handle overlaps
    joined_data <- st_join(cells_use, polygon_data, join = st_intersects, largest = TRUE) %>%
        st_drop_geometry()

    # Apply domain prefix
    joined_data$domain <- paste0(domain_prefix, "_", joined_data$domain)
    polygon_data$domain <- paste0(domain_prefix, "_", polygon_data$domain)

    # Filter polygon_data to only include domains with cells
    polygon_data <- polygon_data %>%
        filter(domain %in% joined_data$domain)

    message(paste("Created", length(unique(polygon_data$domain)), "domains"))
    message(paste("Assigned", sum(!is.na(joined_data$domain)), "cells to domains"))

    result <- list(
        cell_data = joined_data,
        polygon_data = polygon_data,
        group = group,
        group_subset = group_subset
    )

    gc()
    return(result)
}


#' Internal: Merge small domains into neighboring larger domains
#'
#' Must contain cell_data and polygon_data in a list (output from .MakeDomains)
#' Supports dual filtering by target cells and total cells.
#'
#' @param result Output from .MakeDomains
#' @param min.target.cells Minimum TARGET cells per domain (default 10).
#'        Counts only cells matching group/group_subset filter.
#' @param min.total.cells Minimum TOTAL cells per domain (default NULL = no filter).
#'        Counts ALL cells in domain after assign_all_cells expansion.
#' @param group Column name for cell groupings (for target cell counting)
#' @param group_subset Value identifying target cells
#' @return Updated result list with merged domains
#' @keywords internal
.ReduceDomains <- function(result,
                          min.target.cells = 10,
                          min.total.cells = NULL,
                          group = NULL,
                          group_subset = NULL) {

    metadata <- result$cell_data
    polygon_data <- result$polygon_data

    # Use group info from result if not provided
    if (is.null(group)) group <- result$group
    if (is.null(group_subset)) group_subset <- result$group_subset
    if (is.null(group) || is.null(group_subset)) {
        stop("group and group_subset are required for domain reduction.")
    }

    # Calculate TOTAL cells per domain (all cells)
    metadata <- metadata %>%
        group_by(domain) %>%
        mutate(ncell_total = n()) %>%
        ungroup()

    # Calculate TARGET cells per domain (cells matching filter)
    if (!(group %in% names(metadata))) {
        stop(paste("Missing group column in metadata:", group))
    }
    target_counts <- metadata %>%
        filter(!is.na(domain) & .data[[group]] == group_subset) %>%
        group_by(domain) %>%
        summarise(ncell_target = n(), .groups = "drop")

    metadata <- metadata %>%
        left_join(target_counts, by = "domain") %>%
        mutate(ncell_target = ifelse(is.na(ncell_target), 0, ncell_target))

    # Identify domains to merge based on DUAL filtering
    domains_merge <- c()

    # Filter 1: By TARGET cells
    domains_below_target <- unique(
        metadata$domain[metadata$ncell_target <= min.target.cells & !is.na(metadata$domain)]
    )
    domains_merge <- c(domains_merge, domains_below_target)
    message(paste("Domains below target cell threshold:", length(domains_below_target)))

    # Filter 2: By TOTAL cells (if specified)
    if (!is.null(min.total.cells)) {
        domains_below_total <- unique(
            metadata$domain[metadata$ncell_total <= min.total.cells & !is.na(metadata$domain)]
        )
        domains_merge <- c(domains_merge, domains_below_total)
        message(paste("Domains below total cell threshold:", length(domains_below_total)))
    }

    # Deduplicate
    domains_merge <- unique(domains_merge)

    if (length(domains_merge) == 0) {
        message("No domains below threshold, nothing to merge")
        result$cell_data <- metadata
        return(result)
    }

    message(paste("Merging", length(domains_merge), "small domains"))

    small_polygons <- polygon_data %>%
        filter(domain %in% domains_merge)
    polygon_data <- polygon_data %>%
        filter(!domain %in% domains_merge)

    # Find neighbors via intersection
    neighbors <- st_intersects(small_polygons, polygon_data, sparse = FALSE)

    if (nrow(small_polygons) > 0) {
        for (i in 1:nrow(small_polygons)) {
            sp <- small_polygons[i, ]

            neighbor_indices <- which(neighbors[i, ])
            neighbor_indices <- neighbor_indices[neighbor_indices != i]  # Exclude self

            if (!is_empty(neighbor_indices)) {
                # Merge with the first neighbor
                neighbor_index <- neighbor_indices[1]
                neighbor_polygon <- polygon_data[neighbor_index, ]

                # Merge geometries
                merged_polygon <- st_union(neighbor_polygon, sp) %>%
                    select(-domain.1)
                polygon_data[merged_polygon$domain, ] <- merged_polygon

                # Update domain in metadata
                metadata <- metadata %>%
                    mutate(domain = ifelse(domain == sp$domain, merged_polygon$domain, domain))

                message(paste("Merged", sp$domain, "into", merged_polygon$domain))
            } else {
                # No neighbors - drop domain entirely
                polygon_data <- polygon_data %>%
                    filter(domain != sp$domain)
                metadata <- metadata %>%
                    mutate(domain = ifelse(domain == sp$domain, NA_character_, domain))

                message(paste("Removed isolated domain:", sp$domain))
            }
        }
    }

    # Recalculate cell counts after merging
    metadata <- metadata %>%
        group_by(domain) %>%
        mutate(ncell_total = n()) %>%
        ungroup()

    message(paste("Final domain count:", length(unique(na.omit(metadata$domain)))))

    result$cell_data <- metadata
    result$polygon_data <- polygon_data
    return(result)
}


#' Main entry point - called from Python via r_bridge
#' Named to match Python function: make_spatial_domains()
#'
#' @param input_csv Path to input CSV
#' @param output_csv Path to output CSV for cell-domain mapping
#' @param group Column name for cell groupings
#' @param group_subset Value to filter for target cells
#' @param cell_dist Buffer distance
#' @param shrink_margin Shrink margin
#' @param domain_prefix Domain name prefix
#' @param min_target_cells_domain Minimum TARGET cells per domain
#' @param min_total_cells_domain Minimum TOTAL cells per domain (NULL = no filter)
#' @param assign_all_cells Whether to assign all cells
make_spatial_domains <- function(input_csv, output_csv, group, group_subset,
                                 cell_dist = 225, shrink_margin = 25,
                                 domain_prefix = "domain",
                                 min_target_cells_domain = 10,
                                 min_total_cells_domain = NULL,
                                 assign_all_cells = TRUE) {

    # Create domains (internal helper)
    result <- .MakeDomains(
        input_csv = input_csv,
        group = group,
        group_subset = group_subset,
        cell_dist = cell_dist,
        shrink_margin = shrink_margin,
        domain_prefix = domain_prefix,
        assign_all_cells = assign_all_cells
    )

    # Reduce small domains (internal helper) with dual filtering
    result <- .ReduceDomains(
        result = result,
        min.target.cells = min_target_cells_domain,
        min.total.cells = min_total_cells_domain,
        group = group,
        group_subset = group_subset
    )

    # Write output CSV
    write.csv(result$cell_data, output_csv, row.names = FALSE)

    # Return summary as JSON (printed to stdout, parsed by Python)
    summary <- list(
        n_domains = length(unique(na.omit(result$cell_data$domain))),
        n_cells_assigned = sum(!is.na(result$cell_data$domain)),
        n_cells_total = nrow(result$cell_data),
        domains = unique(result$cell_data$domain[!is.na(result$cell_data$domain)])
    )

    cat(toJSON(summary, auto_unbox = TRUE))
}
