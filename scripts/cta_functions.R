# Project Information ----
# SSECR: Diversity-Stability Relationships
# Script Objective: Compositional Stability CTA
# Co-Leads: James Sturges & Dr. Junna Wang
# Last Updated: February 9th, 2026
# 0. Libraries / Saving Outputs / Diagnostic Files ----
library(tidyverse)
library(vegan)
library(ecotraj)
library(RColorBrewer)
library(here)

# helper: save diagnostic CSV (kept separate from CTA metrics)
save_diag <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  write_csv(df, path)
}

# helper: standardized path builder for CTA outputs
make_save_path <- function(site_id, group, method, filename) {
  base_dir <- here("outputs", "cta", site_id, group, method)
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(base_dir, filename)
}
# 1. Function: match_site_years cleans producer/consumer data sets by matching site-years ----
# Returns filtered versions and a diagnostics tibble of removed combos
match_site_years <- function(prod_df, con_df, by_plot = TRUE) {
  # expects both tables have columns: plot, year (can be numeric or character)
  prod_df <- prod_df %>% mutate(year = as.character(year))
  con_df  <- con_df  %>% mutate(year = as.character(year))
  
  # plots that appear in both datasets
  common_plots <- intersect(unique(prod_df$plot), unique(con_df$plot))
  
  # For each common plot, keep only years present in both
  keep_pairs <- map_df(common_plots, function(p) {
    years_prod <- prod_df %>% filter(plot == p) %>% pull(year) %>% unique()
    years_con  <- con_df  %>% filter(plot == p) %>% pull(year) %>% unique()
    tibble(plot = p, keep_year = intersect(years_prod, years_con),
           prod_years = paste(sort(years_prod), collapse = ";"),
           con_years  = paste(sort(years_con), collapse = ";"))
  })
  
  # build diagnostics: removed years per plot
  diag_removed <- keep_pairs %>%
    mutate(removed_prod = map2(prod_years, keep_year, ~ setdiff(str_split(.x, ";")[[1]], .y) %>% paste(collapse = ";")),
           removed_con  = map2(con_years,  keep_year, ~ setdiff(str_split(.x, ";")[[1]], .y) %>% paste(collapse = ";"))) %>%
    select(plot, keep_year, removed_prod, removed_con)
  
  # filter the original dfs
  prod_filt <- prod_df %>% filter(plot %in% common_plots) %>% semi_join(keep_pairs %>% select(plot, keep_year) %>% unnest(cols = c(keep_year)) %>% rename(year = keep_year), by = c("plot", "year"))
  con_filt  <- con_df  %>% filter(plot %in% common_plots) %>% semi_join(keep_pairs %>% select(plot, keep_year) %>% unnest(cols = c(keep_year)) %>% rename(year = keep_year), by = c("plot", "year"))
  
  list(prod = prod_filt, con = con_filt, diagnostics = diag_removed)
}

# 2. Function: Prepare_cta_matrix creates community assemblage matrix ----
# - method: "bray", "jaccard", "hellinger", "chord"
# - returns list(community_matrix, metadata, removed_zero_rows_tbl)
prepare_cta_matrix <- function(df, method = "bray", nmds_trymax = 100) {
  # Ensure year is numeric where possible for ordering
  df <- df %>% mutate(year = if_else(is.na(as.numeric(year)), year, as.character(as.numeric(year))))
  df <- df %>% mutate(site.year = paste(plot, year, sep = "-"))
  
  # Metadata columns we want to drop from the community matrix if present
  metadata_cols <- c("site", "ecosystem", "habitat_broad", "habitat_fine",
                     "guild", "scale_abundance", "biome", "plot", "year", "unit_abundance")
  
  # Build community matrix (rows = site.year)
  comm_matrix <- df %>%
    column_to_rownames(var = "site.year") %>%
    select(where(is.numeric)) %>%
    select(-any_of(metadata_cols))  # safe drop
  
  # Remove taxa that are all zeros across all rows (no information)
  nonzero_taxa <- which(colSums(comm_matrix, na.rm = TRUE) > 0)
  comm_matrix <- comm_matrix[, nonzero_taxa, drop = FALSE]
  
  # Identify rows (site.year) that are all zero → can't compute distances/ordination
  zero_row_idx <- which(rowSums(comm_matrix, na.rm = TRUE) == 0)
  removed_zero_rows_tbl <- tibble(
    site.year = rownames(comm_matrix)[zero_row_idx],
    reason = "all-zero row (removed for distance/ordination)"
  )
  
  if (length(zero_row_idx) > 0) {
    comm_matrix <- comm_matrix[-zero_row_idx, , drop = FALSE]
  }
  
  # Method-specific transforms
  if (tolower(method) %in% c("jaccard", "jaccardian", "jaccardian_distance")) {
    # Convert to presence / absence (1/0)
    comm_matrix_pa <- (comm_matrix > 0) * 1
    comm_matrix <- comm_matrix_pa
  } else if (tolower(method) == "hellinger") {
    # Hellinger: decostand(..., method="hellinger") in vegan
    comm_matrix <- vegan::decostand(comm_matrix, method = "hellinger")
  } else if (tolower(method) == "chord") {
    # Chord: normalized to unit length per row (decostand "normalize")
    comm_matrix <- vegan::decostand(comm_matrix, method = "normalize")
  } else {
    # For "bray" and others, leave raw numeric (counts or biomass)
    # Consider optionally standardizing or log-transforming upstream if desired
    comm_matrix <- comm_matrix
  }
  
  # NMDS: ensure using a distance that makes sense — for jaccard use distance="jaccard"; for hellinger/chord we can use "euclidean"
  nmds_distance <- case_when(
    tolower(method) == "jaccard" ~ "jaccard",
    tolower(method) %in% c("hellinger", "chord") ~ "euclidean",
    TRUE ~ "bray"
  )
  
  # Run NMDS; make vector_nums numeric and ordered.
  # Use tryCatch to avoid breaking whole pipeline for failures - return NA coords if fails.
  nmds <- tryCatch(
    vegan::metaMDS(comm_matrix, distance = nmds_distance, trymax = nmds_trymax, autotransform = FALSE),
    error = function(e) {
      warning("metaMDS failed: ", e$message)
      return(NULL)
    }
  )
  
  if (!is.null(nmds)) {
    nmds_coords <- data.frame(nmds$points, site.year = rownames(comm_matrix))
  } else {
    # create an empty coords frame so later join won't fail; PCoA plotting may fail if no coords
    nmds_coords <- tibble(site.year = rownames(comm_matrix))
  }
  
  # Attach coords back to original metadata (but only for rows that remain in comm_matrix)
  df_joined <- df %>%
    mutate(site.year = paste(plot, year, sep = "-")) %>%
    semi_join(tibble(site.year = rownames(comm_matrix)), by = "site.year") %>%
    left_join(nmds_coords, by = "site.year") %>%
    mutate(
      year = as.numeric(year),
      vector_nums = as.integer(year - min(year, na.rm = TRUE) + 1)
    ) %>%
    arrange(plot, year) # ensure ordering by numeric year within plot
  
  list(
    community_matrix = comm_matrix,
    metadata = df_joined,
    removed_zero_rows = removed_zero_rows_tbl
  )
}

# 3. Function: generate_custom_palette creates unique trajectory colors per plot ----
generate_custom_palette <- function(sites) {
  uniq <- unique(sites)
  num_sites <- length(uniq)
  n_palette <- min(max(3, num_sites), 12) # get palette baseline size
  base_colors <- brewer.pal(n_palette, "Set3")
  colors <- colorRampPalette(base_colors)(num_sites)
  names(colors) <- uniq
  # return named vector mapping plot -> hex color
  colors
}
# 4. Function: run_cta_plot creates plot structure ----
run_cta_plot <- function(D, metadata, palette, output_path) {
  # D: distance object or matrix
  # metadata must contain plot and vector_nums (numeric) and be ordered sensibly
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  png(filename = output_path, width = 12, height = 8, res = 300, units = "in")
  # Ensure unique plot ordering is consistent
  unique_plots <- unique(metadata$plot)
  traj_colors <- unname(palette[unique_plots]) # one color per unique plot, in plot order
  
  trajectoryPCoA(
    D,
    sites   = metadata$plot,
    surveys = metadata$vector_nums,
    traj.colors = traj_colors,
    lwd = 1,
    survey.labels = TRUE
  )
  legend("topright", legend = unique_plots, col = traj_colors, lwd = 2, bty = "n", cex = 0.8)
  dev.off()
}

# 5. Function: extract_all_cta_metrics + run full pipeline ----
extract_all_cta_metrics <- function(D, metadata, site_id, output_dir, method = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # --- Lengths ---
  lengths <- trajectoryLengths(D, sites = metadata$plot, surveys = metadata$vector_nums)
  lengths_df <- as.data.frame(lengths) %>% tibble::rownames_to_column("plot")
  lengths_df$total_trajectory <- if ("Trajectory" %in% names(lengths_df)) lengths_df$Trajectory else lengths_df[[ncol(lengths_df)]]
  segment_cols <- grep("^S\\d+$", names(lengths_df), value = TRUE)
  lengths_df$mean_segment_length <- if (length(segment_cols) > 0) rowMeans(lengths_df[, segment_cols, drop = FALSE], na.rm = TRUE) else NA_real_
  lengths_df$site <- site_id
  write_csv(lengths_df, file.path(output_dir, paste0(site_id, "_lengths", if(!is.null(method)) paste0("_", method), ".csv")))
  write_csv(lengths_df %>% select(site, plot, mean_segment_length, total_trajectory),
            file.path(output_dir, paste0(site_id, "_lengths_summary", if(!is.null(method)) paste0("_", method), ".csv")))
  
  # --- Angles ---
  angles <- trajectoryAngles(D, sites = metadata$plot, surveys = metadata$vector_nums)
  angles_df <- as.data.frame(angles) %>% tibble::rownames_to_column("plot")
  numeric_angle_cols <- names(angles_df)[sapply(angles_df, is.numeric)]
  angles_df$mean_angle <- if (length(numeric_angle_cols) > 0) rowMeans(angles_df[, numeric_angle_cols, drop = FALSE], na.rm = TRUE) else NA_real_
  angles_df$site <- site_id
  write_csv(angles_df, file.path(output_dir, paste0(site_id, "_angles", if(!is.null(method)) paste0("_", method), ".csv")))
  write_csv(angles_df %>% select(site, plot, mean_angle),
            file.path(output_dir, paste0(site_id, "_angles_summary", if(!is.null(method)) paste0("_", method), ".csv")))
  
  # --- Directionality ---
  dirct <- trajectoryDirectionality(D, sites = metadata$plot, surveys = metadata$vector_nums)
  dirct_df <- tibble::tibble(plot = names(dirct),
                             directionality = as.numeric(dirct),
                             site = site_id)
  write_csv(dirct_df, file.path(output_dir, paste0(site_id, "_directionality", if(!is.null(method)) paste0("_", method), ".csv")))
  
  # --- Convergence ---
  conv <- trajectoryConvergence(D, sites = metadata$plot, surveys = metadata$vector_nums)
  conv_df <- as.data.frame(conv) %>%
    tibble::rownames_to_column("from_plot") %>%
    tidyr::pivot_longer(-from_plot, names_to = "to_plot", values_to = "convergence") %>%
    dplyr::mutate(site = site_id, method = ifelse(is.null(method), NA_character_, as.character(method)))
  write_csv(conv_df, file.path(output_dir, paste0(site_id, "_convergence", if(!is.null(method)) paste0("_", method), ".csv")))
  
  invisible(list(lengths_full = lengths_df,
                 lengths_summary = lengths_df %>% select(site, plot, mean_segment_length, total_trajectory),
                 angles_full = angles_df,
                 angles_summary = angles_df %>% select(site, plot, mean_angle),
                 directionality = dirct_df,
                 convergence = conv_df))
}

run_cta_pipeline <- function(df, site_id, fig_path, output_dir,
                             method = "bray", nmds_trymax = 100) {
  prep <- prepare_cta_matrix(df, method = method, nmds_trymax = nmds_trymax)
  if (is.null(prep$community_matrix) || nrow(prep$community_matrix) == 0) {
    warning("No rows in community matrix after filtering. Skipping site: ", site_id, " (method: ", method, ").")
    return(invisible(NULL))
  }
  
  if (!is.null(prep$removed_zero_rows) && nrow(prep$removed_zero_rows) > 0) {
    diag_path <- file.path(output_dir, paste0("removed_zero_rows_", site_id, "_", method, ".csv"))
    save_diag(prep$removed_zero_rows, diag_path)
  }
  
  prep$metadata <- prep$metadata %>%
    dplyr::group_by(plot) %>%
    dplyr::arrange(year, .by_group = TRUE) %>%
    dplyr::mutate(vector_nums = as.integer(dplyr::dense_rank(year))) %>%
    dplyr::ungroup()
  
  veg_method <- switch(tolower(method),
                       "jaccard" = "jaccard",
                       "hellinger" = "euclidean",
                       "chord" = "euclidean",
                       "bray" = "bray",
                       tolower(method))
  
  D <- vegan::vegdist(prep$community_matrix, method = veg_method)
  
  pal <- generate_custom_palette(prep$metadata$plot)
  run_cta_plot(D, prep$metadata, pal, fig_path)
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  metrics <- extract_all_cta_metrics(D, prep$metadata, site_id, output_dir, method = method)
  
  invisible(metrics)
}

# 6. Function: Helper function to load model output length files ----
load_lengths <- function(path, trophic, method) {
  read_csv(path) %>%
    mutate(
      trophic = trophic,
      method = method
    )
}