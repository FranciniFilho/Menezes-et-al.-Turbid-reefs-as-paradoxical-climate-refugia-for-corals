# ============================================================================ #
# SCRIPT: JSDM_Data_Integration_V4.R
# PURPOSE: Prepare unified dataset for JSDM hypothesis-testing pipeline
# VERSION: 4.0 (Corrected from V3)
# ============================================================================ #

rm(list = ls())
gc()

# Required packages
libs <- c("dplyr", "tidyr", "readxl", "janitor", "digest")
invisible(lapply(libs, library, character.only = TRUE))

# ============================================================================ #
# CONFIGURATION
# ============================================================================ #

cat("\n", rep("=", 70), "\n", sep = "")
cat("=== JSDM V4: DATA INTEGRATION ===\n")
cat(rep("=", 70), "\n", sep = "")

# --- Input Paths ---
path_benthos <- "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######22.04.23/DATA/BENTHOS TEMPORAL CLEAN_v3.xlsx"

# --- MULTI-CV DATASET PATHS ---
dataset_paths <- list(
  CV_02 = "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_RESULTS/#####output_local_PCA_CV_2_FINAL/dados_abundancia_integrados_long_format.csv",
  CV_30 = "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_RESULTS/#####output_local_PCA_CV_30_FINAL/dados_abundancia_integrados_long_format.csv",
  CV_ALL = "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_RESULTS/#####output_local_PCA_CV_all_FINAL/dados_abundancia_integrados_long_format.csv"
)

# --- Output Path ---
output_dir <- "C:/Users/rbfra/OneDrive/JSDM_Dirichlet_HypothesisTesting_V4"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================ #
# FUNCTION: prepare_jsdm_data_for_cv
# PURPOSE: Integrate benthos + environmental data for a specific CV scenario
# CRITICAL FIX: Uses explicit path parameter instead of undefined variable
# ============================================================================ #

prepare_jsdm_data_for_cv <- function(cv_name, cv_path, benthos_path, output_dir) {
  #' @param cv_name Name of CV scenario (CV_02, CV_30, CV_ALL)
  #' @param cv_path Full path to environmental data CSV
  #' @param benthos_path Path to benthos community Excel file
  #' @param output_dir Directory to save integrated datasets
  #' @return Integrated data frame

  cat(sprintf("\n=== Processing CV Dataset: %s ===\n", cv_name))

  # --- Step 1: Load integrated environmental data (with ARCH) ---
  cat("  [1/5] Loading environmental data...\n")
  env_data <- tryCatch(
    {
      read.csv2(cv_path, stringsAsFactors = TRUE)
    },
    error = function(e) {
      read.csv(cv_path, stringsAsFactors = TRUE)
    }
  )
  colnames(env_data) <- toupper(colnames(env_data))

  # Verify required columns
  required_cols <- c(
    "SITE", "HAB", "ARCH", "PC1_MAGNITUDE", "PC2_MAGNITUDE",
    "PC1_VARIABILITY", "PC2_VARIABILITY", "DEPTH_M"
  )
  missing_cols <- setdiff(required_cols, colnames(env_data))
  if (length(missing_cols) > 0) {
    warning(sprintf(
      "Missing columns in environmental data: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  cat(sprintf("     -> Loaded %d rows from environmental data\n", nrow(env_data)))

  # --- Step 2: Load benthos community data ---
  cat("  [2/5] Loading benthos community data...\n")
  benthos_raw <- read_excel(benthos_path, guess_max = 5000) %>%
    janitor::clean_names()

  cat(sprintf("     -> Loaded %d rows from benthos data\n", nrow(benthos_raw)))

  # Identify organism columns (exclude metadata)
  metadata_cols <- c("site", "reef", "hab", "year", "transect", "point", "photo_code")
  organism_cols <- setdiff(names(benthos_raw), metadata_cols)

  cat(sprintf("     -> Identified %d organism columns\n", length(organism_cols)))

  # --- Step 3: Extract mappings from environmental data ---
  cat("  [3/5] Extracting SITE->ARCH and SITE->DEPTH mappings...\n")

  site_arch_map <- env_data %>%
    select(SITE, ARCH) %>%
    distinct()

  site_depth_map <- env_data %>%
    select(SITE, DEPTH_M) %>%
    distinct() %>%
    group_by(SITE) %>%
    summarise(DEPTH_M = mean(DEPTH_M, na.rm = TRUE), .groups = "drop")

  # PCA scores per SITE-HAB combination
  pca_scores <- env_data %>%
    select(SITE, HAB, PC1_MAGNITUDE, PC2_MAGNITUDE, PC1_VARIABILITY, PC2_VARIABILITY) %>%
    distinct() %>%
    group_by(SITE, HAB) %>%
    summarise(across(everything(), ~ mean(., na.rm = TRUE)), .groups = "drop")

  cat(sprintf("     -> ARCH mappings: %d sites\n", nrow(site_arch_map)))
  cat(sprintf("     -> PCA scores: %d SITE-HAB combinations\n", nrow(pca_scores)))

  # --- Step 4: Merge all data ---
  cat("  [4/5] Merging datasets...\n")

  # Standardize column names
  benthos_raw <- benthos_raw %>%
    rename(SITE = site, REEF = reef, HAB = hab)

  n_start <- nrow(benthos_raw)

  benthos_merged <- benthos_raw %>%
    left_join(site_arch_map, by = "SITE") %>%
    left_join(site_depth_map, by = "SITE") %>%
    left_join(pca_scores, by = c("SITE", "HAB")) %>%
    filter(!is.na(ARCH))

  n_after_filter <- nrow(benthos_merged)

  cat(sprintf("     -> Records before ARCH filter: %d\n", n_start))
  cat(sprintf(
    "     -> Records after ARCH filter: %d (removed %d)\n",
    n_after_filter, n_start - n_after_filter
  ))

  # --- Step 5: Clean and normalize organism data ---
  cat("  [5/5] Processing organism data...\n")

  # Select only organism columns with metadata
  final_cols <- c(
    "SITE", "REEF", "HAB", "ARCH", "DEPTH_M",
    "PC1_MAGNITUDE", "PC2_MAGNITUDE", "PC1_VARIABILITY", "PC2_VARIABILITY"
  )

  # Check which organism columns exist in merged data
  available_organisms <- intersect(organism_cols, names(benthos_merged))
  missing_organisms <- setdiff(organism_cols, names(benthos_merged))

  if (length(missing_organisms) > 0) {
    warning(sprintf(
      "Missing organism columns after merge: %s",
      paste(missing_organisms, collapse = ", ")
    ))
  }

  # Create final dataset
  benthos_final <- benthos_merged %>%
    select(all_of(c(final_cols, available_organisms)))

  # Report summary
  cat(sprintf("\n  Summary for %s:\n", cv_name))
  cat(sprintf("    - Total observations: %d\n", nrow(benthos_final)))
  cat(sprintf("    - Unique sites: %d\n", n_distinct(benthos_final$SITE)))
  cat(sprintf("    - Unique habitats: %d\n", n_distinct(benthos_final$HAB)))
  cat(sprintf("    - Organism columns: %d\n", length(available_organisms)))

  # Save CV-specific integrated dataset
  output_file <- file.path(output_dir, sprintf("JSDM_integrated_%s.rds", cv_name))
  saveRDS(list(
    data = benthos_final,
    organism_cols = available_organisms,
    cv_name = cv_name,
    data_hash = digest::digest(benthos_final),
    metadata = list(
      n_observations = nrow(benthos_final),
      n_sites = n_distinct(benthos_final$SITE),
      n_habitats = n_distinct(benthos_final$HAB),
      n_organisms = length(available_organisms)
    )
  ), output_file)

  cat(sprintf("  -> Saved to: %s\n", output_file))

  return(benthos_final)
}

# ============================================================================ #
# MAIN EXECUTION
# ============================================================================ #

cat("\n### PROCESSING ALL CV DATASETS ###\n")

# Verify input files exist
cat("\nVerifying input files...\n")

if (!file.exists(path_benthos)) {
  stop(sprintf("Benthos file not found: %s", path_benthos))
}
cat(sprintf("  Benthos file: OK\n"))

for (cv_name in names(dataset_paths)) {
  cv_path <- dataset_paths[[cv_name]]
  if (!file.exists(cv_path)) {
    stop(sprintf("Dataset file not found for %s: %s", cv_name, cv_path))
  }
  cat(sprintf("  %s: OK\n", cv_name))
}

cat("\n--- Starting data integration ---\n")

# Process all CV datasets
all_datasets <- list()
summaries <- list()

for (cv_name in names(dataset_paths)) {
  result <- prepare_jsdm_data_for_cv(
    cv_name, dataset_paths[[cv_name]],
    path_benthos, output_dir
  )
  all_datasets[[cv_name]] <- result

  # Store summary
  summaries[[cv_name]] <- list(
    n_obs = nrow(result),
    n_sites = n_distinct(result$SITE),
    n_habs = n_distinct(result$HAB)
  )
}

# ============================================================================ #
# GENERATE SUMMARY REPORT
# ============================================================================ #

cat("\n", rep("=", 70), "\n", sep = "")
cat("=== DATA INTEGRATION SUMMARY ===\n")
cat(rep("=", 70), "\n", sep = "")

summary_df <- do.call(rbind, lapply(names(summaries), function(cv_name) {
  data.frame(
    Dataset = cv_name,
    Observations = summaries[[cv_name]]$n_obs,
    Sites = summaries[[cv_name]]$n_sites,
    Habitats = summaries[[cv_name]]$n_habs,
    stringsAsFactors = FALSE
  )
}))

print(summary_df)

# Save summary
write.csv(summary_df,
  file.path(output_dir, "data_integration_summary.csv"),
  row.names = FALSE
)

# Verify consistency across datasets
common_sites <- Reduce(intersect, lapply(all_datasets, function(d) unique(d$SITE)))
common_habs <- Reduce(intersect, lapply(all_datasets, function(d) unique(d$HAB)))

cat(sprintf("\nConsistency check:\n"))
cat(sprintf("  - Common sites across all CVs: %d\n", length(common_sites)))
cat(sprintf("  - Common habitats across all CVs: %d\n", length(common_habs)))

# Get common organism columns
all_organism_cols <- lapply(all_datasets, function(d) {
  organism_cols <- setdiff(names(d), c(
    "SITE", "REEF", "HAB", "ARCH", "DEPTH_M",
    "PC1_MAGNITUDE", "PC2_MAGNITUDE",
    "PC1_VARIABILITY", "PC2_VARIABILITY"
  ))
  organism_cols
})
common_organisms <- Reduce(intersect, all_organism_cols)

cat(sprintf("  - Common organisms across all CVs: %d\n", length(common_organisms)))

# Save common organisms list
saveRDS(
  common_organisms,
  file.path(output_dir, "common_organisms_all_cv.rds")
)

cat(sprintf("\nAll integrated datasets saved to: %s\n", output_dir))
cat("\nNext step: Run JSDM_Dirichlet_HypothesisTesting_Engine_V4.R\n")
