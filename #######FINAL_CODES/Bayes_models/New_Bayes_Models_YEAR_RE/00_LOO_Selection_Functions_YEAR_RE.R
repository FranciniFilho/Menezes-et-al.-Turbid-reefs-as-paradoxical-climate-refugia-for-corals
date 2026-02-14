# ============================================================================
# 00_LOO_Selection_Functions_YEAR_RE.R
# ============================================================================
# Extended functions for YEAR Random Effect models
# Based on: 00_LOO_Selection_Functions.R
# Adds: YEAR handling for ZOIB, Health, JSDM models
# Key Change: CROSSED random effects (1|REEF) + (1|YEAR) - NOT nested!
# ============================================================================

# First, load brms (needed for prior() function)
suppressMessages(library(brms))

# Then source the original functions
source("C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_CODES/Bayes_models/00_LOO_Selection_Functions.R")

# ============================================================================
# SECTION 1: PATH DEFINITIONS FOR YEAR RE MODELS
# ============================================================================

# Base output directory for YEAR RE models (separate from original outputs)
BASE_OUTPUT_DIR_YEAR_RE <- "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/"

# Create base directory
dir.create(BASE_OUTPUT_DIR_YEAR_RE, showWarnings = FALSE, recursive = TRUE)

# Subdirectories for YEAR RE models
output_dirs_year_re <- list(
  ZOIB_YEAR_RE    = file.path(BASE_OUTPUT_DIR_YEAR_RE, "ZOIB_Abundance_YEAR_RE/"),
  RGR             = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_RGR/"),
  HEALTH_YEAR_RE  = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_Health_YEAR_RE/"),
  JSDM_YEAR_RE    = file.path(BASE_OUTPUT_DIR_YEAR_RE, "JSDM_Dirichlet_YEAR_RE/")
)

# Raw health data path (with yearly observations)
# Using list.files to avoid # character issues
RAW_HEALTH_EXCEL <- list.files(
  path = "C:/Users/rbfra/OneDrive/",
  pattern = "Vitality and size_new\\.xlsx$",
  recursive = TRUE,
  full.names = TRUE
)[1]

# PCA scores path (for loading existing loadings - NOT recalculating!)
# Using list.files to avoid # character issues
PCA_SCORES_PATH <- list.files(
  path = "C:/Users/rbfra/OneDrive/",
  pattern = "scores_PCA_Saude\\.xlsx$",
  recursive = TRUE,
  full.names = TRUE
)[1]

# ============================================================================
# SECTION 2: UTILITY FUNCTIONS FOR YEAR VALIDATION
# ============================================================================

#' Check YEAR column validity for random effect
#' @param data data.frame to check
#' @param min_levels Minimum number of year levels required (default: 2)
#' @return list with valid (logical), n_years (int), and summary info
check_year_levels <- function(data, min_levels = 2) {

  if (!"YEAR" %in% names(data)) {
    warning("YEAR column not found. Models will not have YEAR random effect.")
    return(list(valid = FALSE, n_years = 0, summary = "No YEAR column"))
  }

  # Convert YEAR to character for consistent handling
  data$YEAR <- as.character(data$YEAR)
  n_years <- length(unique(data$YEAR))
  obs_per_year <- table(data$YEAR)

  cat(sprintf("\n=== YEAR SUMMARY ===\n"))
  cat(sprintf("  Number of years: %d\n", n_years))
  cat(sprintf("  Years present: %s\n", paste(sort(unique(data$YEAR)), collapse = ", ")))
  cat("  Observations per year:\n")
  for (y in names(obs_per_year)) {
    cat(sprintf("    Year %s: %d\n", y, obs_per_year[y]))
  }

  if (n_years < min_levels) {
    warning(sprintf("Only %d year(s) in data. Minimum recommended: %d",
                   n_years, min_levels))
    return(list(valid = FALSE, n_years = n_years,
                summary = sprintf("Insufficient years: %d < %d", n_years, min_levels)))
  }

  # Check for highly unbalanced years
  min_n <- min(obs_per_year)
  max_n <- max(obs_per_year)
  if (max_n / min_n > 10) {
    warning(sprintf("YEAR levels are highly unbalanced (min: %d, max: %d)",
                   min_n, max_n))
  }

  return(list(valid = TRUE, n_years = n_years,
              summary = sprintf("Valid: %d years, range [%d, %d] obs/year",
                               n_years, min_n, max_n)))
}

# ============================================================================
# SECTION 3: DATA PREPARATION FUNCTIONS WITH YEAR
# ============================================================================

#' Prepare ZOIB data with YEAR as crossed random effect
#' @param data_path Path to dados_abundancia_integrados_long_format.csv
#' @param cv_label CV scenario label
#' @return data.frame prepared for brms with YEAR column
prepare_zoib_data_with_year <- function(data_path, cv_label) {

  cat(sprintf("\n=== Preparing ZOIB data with YEAR: %s ===\n", cv_label))

  # Use original prepare_zoib_data() as base
  data <- prepare_zoib_data(data_path, cv_label)

  # Add YEAR handling - YEAR column already exists in the data
  if ("YEAR" %in% names(data)) {
    data$YEAR <- factor(data$YEAR)
    n_years <- length(unique(data$YEAR))
    cat(sprintf("  YEARS in ZOIB data: %s (N=%d)\n",
                paste(unique(data$YEAR), collapse = ", "), n_years))

    # Check for sufficient year levels
    if (n_years < 2) {
      warning(sprintf("Only %d year(s) in ZOIB data. YEAR random effect may be poorly identified.", n_years))
    }
  } else {
    warning("YEAR column not found in ZOIB data. Using (1 | REEF) only.")
    data$YEAR <- factor(NA)
  }

  return(data)
}


#' Prepare health data with YEAR from raw Excel file
#' @param excel_path Path to ######Vitality and size_new.xlsx
#' @param pca_scores_path Path to scores_PCA_Saude.xlsx (existing PCA loadings)
#' @param env_data_path Path to CV-specific environmental data CSV
#' @param final_data_path Path to CV-specific final data CSV with interaction variables (NEW)
#' @param cv_label CV scenario label (CV_02, CV_30, CV_ALL)
#' @return data.frame prepared for brms with YEAR column
prepare_health_data_with_year <- function(excel_path, pca_scores_path,
                                          env_data_path, final_data_path = NULL, cv_label) {

  cat(sprintf("\n=== Preparing HEALTH data with YEAR: %s ===\n", cv_label))

  library(readxl)

  # --- STEP 1: Read raw Excel file ---
  raw_data <- tryCatch({
    read_excel(excel_path)
  }, error = function(e) {
    stop(paste("Cannot read Excel file:", e$message))
  })

  # Standardize column names
  colnames(raw_data) <- toupper(trimws(colnames(raw_data)))

  # Verify required columns
  required_cols <- c("YEAR", "HEALTH %", "BLEACHING %", "DEAD %", "SITE", "COL")
  missing <- setdiff(required_cols, names(raw_data))
  if (length(missing) > 0) {
    # Try alternative column names
    alt_cols <- list(
      "HEALTH %" = c("HEALTH", "HEALTH_PERCENT", "HEALTH_PCT"),
      "BLEACHING %" = c("BLEACHING", "BLEACH", "BLEACHING_PERCENT"),
      "DEAD %" = c("DEAD", "MORTALITY", "DEAD_PERCENT")
    )
    for (orig in names(alt_cols)) {
      if (orig %in% missing) {
        for (alt in alt_cols[[orig]]) {
          if (alt %in% names(raw_data)) {
            raw_data[[orig]] <- raw_data[[alt]]
            missing <- setdiff(missing, orig)
            cat(sprintf("  Using alternative column: %s -> %s\n", alt, orig))
            break
          }
        }
      }
    }
    if (length(missing) > 0) {
      stop(paste("Missing columns in raw Excel:", paste(missing, collapse = ", ")))
    }
  }

  # --- STEP 2: Create unique colony identifier ---
  raw_data$SITE_COL <- paste(raw_data$SITE, raw_data$COL, sep = "_")
  raw_data$YEAR <- as.character(raw_data$YEAR)

  # --- STEP 3: Merge with existing PCA scores ---
  # This ensures HEALTH_PC1 and HEALTH_PC2 are calculated with the SAME
  # PCA loadings as the original analysis (comparability!)
  cat("  Loading existing PCA scores for comparability...\n")

  pca_scores <- tryCatch({
    read_excel(pca_scores_path)
  }, error = function(e) {
    stop(paste("Cannot read PCA scores file:", e$message))
  })

  colnames(pca_scores) <- toupper(colnames(pca_scores))
  pca_scores$SITE_COL <- toupper(trimws(pca_scores$SITE_COL))

  # Verify PCA scores have required columns
  if (!all(c("SITE_COL", "HEALTH_PC1", "HEALTH_PC2") %in% names(pca_scores))) {
    stop("PCA scores file missing required columns: SITE_COL, HEALTH_PC1, HEALTH_PC2")
  }

  # Merge PCA scores by SITE_COL
  health_data <- merge(raw_data,
                       pca_scores[, c("SITE_COL", "HEALTH_PC1", "HEALTH_PC2")],
                       by = "SITE_COL", all.x = TRUE)

  # Check for NA PCA scores
  if (any(is.na(health_data$HEALTH_PC1)) || any(is.na(health_data$HEALTH_PC2))) {
    warning("Some observations have NA PCA scores after merge")
  }

  # --- STEP 4: Read and aggregate environmental data BEFORE merge ---
  cat("  Reading and aggregating environmental data...\n")

  env_data <- tryCatch({
    read.csv2(env_data_path, stringsAsFactors = TRUE)
  }, error = function(e) {
    read.csv(env_data_path, stringsAsFactors = TRUE)
  })

  colnames(env_data) <- toupper(colnames(env_data))
  env_data$YEAR <- as.character(env_data$YEAR)

  # CRITICAL FIX: Aggregate env_data to ONE ROW per SITE-HAB-YEAR
  # The env_data has multiple organisms (MUSSISMILIA, TURF, CCA, etc.)
  # We need to take the mean of environmental variables per SITE-HAB-YEAR
  # to avoid creating a Cartesian product in the merge
  cat("  Aggregating environmental data to one row per SITE-HAB-YEAR...\n")

  # Select only the environmental columns we need
  env_cols_to_keep <- c("SITE", "HAB", "YEAR", "PC1_MAGNITUDE", "PC2_MAGNITUDE",
                        "PC1_VARIABILITY", "PC2_VARIABILITY", "DEPTH_M", "ARCH")

  # Keep only columns that exist
  env_cols_to_keep <- intersect(env_cols_to_keep, names(env_data))

  # Aggregate: take the first (or mean) of each env variable per SITE-HAB-YEAR
  # Since environmental data is constant across organisms at the same site-hab-year
  env_agg <- env_data[, env_cols_to_keep]

  # Remove duplicates based on SITE-HAB-YEAR, keeping first occurrence
  env_agg <- env_agg[!duplicated(env_agg[, c("SITE", "HAB", "YEAR")]), ]

  cat(sprintf("    Environmental data: %d rows -> %d unique SITE-HAB-YEAR combinations\n",
              nrow(env_data), nrow(env_agg)))

  # Merge by SITE, HAB, and YEAR
  health_data <- merge(health_data, env_agg,
                       by = c("SITE", "HAB", "YEAR"), all.x = TRUE)

  # --- STEP 5: Create derived variables (same as original) ---
  # Create REEF from first 3 characters of SITE
  health_data$REEF <- factor(substr(as.character(health_data$SITE), 1, 3))

  # Create HABMERGED (PA vs RR_TP)
  health_data$HAB <- as.character(health_data$HAB)
  health_data$HABMERGED <- ifelse(health_data$HAB %in% c("RR", "TP"), "RR_TP", health_data$HAB)
  health_data$HABMERGED <- factor(health_data$HABMERGED, levels = c("PA", "RR_TP"))

  # --- STEP 5.5: Read and merge interaction data from final data file (NEW) ---
  if (!is.null(final_data_path) && file.exists(final_data_path)) {
    cat("  Reading interaction data from final data file...\n")

    final_data <- tryCatch({
      read.csv2(final_data_path, stringsAsFactors = TRUE)
    }, error = function(e) {
      read.csv(final_data_path, stringsAsFactors = TRUE)
    })

    colnames(final_data) <- toupper(colnames(final_data))

    # Select only interaction variables and key columns for merging
    # Note: final data doesn't have YEAR - merge by SITE_COL or SITE+HAB
    interaction_cols <- c("SITE_COL", "SITE", "HAB", "PC1_INTERACAO", "PC2_INTERACAO")
    interaction_cols <- intersect(interaction_cols, names(final_data))

    if (all(c("PC1_INTERACAO", "PC2_INTERACAO") %in% interaction_cols)) {
      # Use SITE_COL for merging (unique colony identifier in final data)
      if ("SITE_COL" %in% interaction_cols) {
        # Create SITE_COL in health_data if it doesn't exist
        if (!"SITE_COL" %in% names(health_data)) {
          health_data$SITE_COL <- paste(health_data$SITE, health_data$COL, sep = "_")
        }

        # Remove duplicates based on SITE_COL
        final_data_agg <- final_data[, interaction_cols]
        final_data_agg <- final_data_agg[!duplicated(final_data_agg$SITE_COL), ]

        # Merge interaction data by SITE_COL
        health_data <- merge(health_data, final_data_agg,
                           by = "SITE_COL", all.x = TRUE)

        cat(sprintf("    Interaction data merged: %d rows\n", nrow(health_data)))
      } else {
        # Fallback: merge by SITE and HAB
        final_data_agg <- final_data[, interaction_cols]
        final_data_agg <- final_data_agg[!duplicated(final_data_agg[, c("SITE", "HAB")]), ]

        # Merge interaction data
        health_data <- merge(health_data, final_data_agg,
                           by = c("SITE", "HAB"), all.x = TRUE)

        cat(sprintf("    Interaction data merged: %d rows\n", nrow(health_data)))
      }
    } else {
      cat("    WARNING: Interaction variables not found in final data file\n")
    }
  } else {
    cat("    WARNING: No final data path provided for interaction variables\n")
  }

  # --- STEP 6: Scale predictors (remove underscores for brms) ---
  # Scale existing environmental variables
  for (col in c("PC1_MAGNITUDE", "PC2_MAGNITUDE",
                "PC1_VARIABILITY", "PC2_VARIABILITY", "DEPTH_M")) {
    if (col %in% names(health_data)) {
      scaled_name <- gsub("_", "", col)
      health_data[[scaled_name]] <- scale(as.numeric(health_data[[col]]))[, 1]
    }
  }

  # NEW: Scale biological interaction variables (PC1_INTERACAO, PC2_INTERACAO)
  for (col in c("PC1_INTERACAO", "PC2_INTERACAO")) {
    if (col %in% names(health_data)) {
      # Remove underscore for brms compatibility (s() doesn't accept underscores)
      scaled_name <- gsub("_", "", col)  # PC1_INTERACAO -> PC1INTERACAO
      health_data[[scaled_name]] <- scale(as.numeric(health_data[[col]]))[, 1]
      cat(sprintf("  Scaled %s -> %s\n", col, scaled_name))
    } else {
      cat(sprintf("  WARNING: %s not found in data\n", col))
    }
  }

  # --- STEP 7: Convert YEAR to factor for random effect ---
  health_data$YEAR <- factor(health_data$YEAR)

  # --- STEP 8: Add metadata ---
  health_data$CVLABEL <- cv_label

  # --- STEP 9: Validate YEAR ---
  year_check <- check_year_levels(health_data)

  cat(sprintf("  Health data prepared: N=%d, %d colonies, %d REEFs, %d years\n",
              nrow(health_data), length(unique(health_data$SITE_COL)),
              length(unique(health_data$REEF)), length(unique(health_data$YEAR))))

  return(as.data.frame(health_data))
}


#' Prepare JSDM data with YEAR as crossed random effect
#' @param data_path Path to dados_abundancia_integrados_long_format.csv
#' @param cv_label CV scenario label
#' @return data.frame with community composition proportions and YEAR
prepare_jsdm_data_with_year <- function(data_path, cv_label) {

  cat(sprintf("\n=== Preparing JSDM data with YEAR: %s ===\n", cv_label))

  # Read and prepare data (same as original but we keep YEAR)
  data <- tryCatch({
    read.csv2(data_path, stringsAsFactors = TRUE)
  }, error = function(e) {
    read.csv(data_path, stringsAsFactors = TRUE)
  })

  colnames(data) <- toupper(colnames(data))

  # Convert coverage
  data$COBERTURA <- as.numeric(gsub(",", ".", as.character(data$COBERTURA)))
  data <- data[!is.na(data$COBERTURA), ]
  data$COVER_PROP <- data$COBERTURA / 100

  # Define organism categories
  data$CATEGORY <- case_when(
    data$ORGANISMO == "MUSSISMILIA_HISPIDA" ~ "MUSSISMILIA",
    grepl("TURF", data$ORGANISMO, ignore.case = TRUE) ~ "TURF",
    grepl("CCA|CORALLINE", data$ORGANISMO, ignore.case = TRUE) ~ "CCA",
    grepl("CYANO", data$ORGANISMO, ignore.case = TRUE) ~ "CYANO",
    grepl("MACRO", data$ORGANISMO, ignore.case = TRUE) ~ "MACROALGAE",
    TRUE ~ "OTHER"
  )

  # Create SAMPLE_UNIT with YEAR INCLUDED (this is the key fix!)
  data$SAMPLE_UNIT <- paste(data$SITE, data$HAB, data$YEAR, sep = "_")

  # Aggregate coverage by category and sampling unit (including YEAR!)
  composition <- data %>%
    group_by(SAMPLE_UNIT, SITE, HAB, YEAR, ARCH,
             PC1_MAGNITUDE, PC2_MAGNITUDE,
             PC1_VARIABILITY, PC2_VARIABILITY,
             DEPTH_M) %>%
    summarise(
      MUSSISMILIA_raw = sum(COVER_PROP[CATEGORY == "MUSSISMILIA"], na.rm = TRUE),
      TURF_raw = sum(COVER_PROP[CATEGORY == "TURF"], na.rm = TRUE),
      CCA_raw = sum(COVER_PROP[CATEGORY == "CCA"], na.rm = TRUE),
      CYANO_raw = sum(COVER_PROP[CATEGORY == "CYANO"], na.rm = TRUE),
      MACROALGAE_raw = sum(COVER_PROP[CATEGORY == "MACROALGAE"], na.rm = TRUE),
      OTHER_raw = sum(COVER_PROP[CATEGORY == "OTHER"], na.rm = TRUE),
      .groups = "drop"
    )

  # Normalize to sum = 1
  total_cover <- rowSums(composition[, c("MUSSISMILIA_raw", "TURF_raw", "CCA_raw",
                                        "CYANO_raw", "MACROALGAE_raw", "OTHER_raw")])
  total_cover[total_cover == 0] <- 1

  composition$MUSSISMILIA_prop <- composition$MUSSISMILIA_raw / total_cover
  composition$TURF_prop <- composition$TURF_raw / total_cover
  composition$CCA_prop <- composition$CCA_raw / total_cover
  composition$CYANO_prop <- composition$CYANO_raw / total_cover
  composition$MACROALGAE_prop <- composition$MACROALGAE_raw / total_cover
  composition$OTHER_prop <- composition$OTHER_raw / total_cover

  # Adjust zeros (Dirichlet doesn't allow exactly 0 or 1)
  epsilon <- 0.001
  prop_cols <- c("MUSSISMILIA_prop", "TURF_prop", "CCA_prop",
                 "CYANO_prop", "MACROALGAE_prop", "OTHER_prop")

  for (col in prop_cols) {
    composition[[col]] <- pmax(epsilon, pmin(1 - epsilon, composition[[col]]))
  }

  # Re-normalize after epsilon adjustment
  row_sums <- rowSums(composition[, prop_cols])
  for (col in prop_cols) {
    composition[[col]] <- composition[[col]] / row_sums
  }

  # Remove proportion columns that are all zeros
  # But first, identify which columns to keep
  prop_cols_to_keep <- prop_cols[sapply(prop_cols, function(col) {
    col %in% names(composition) && !all(composition[[col]] <= epsilon * 2)
  })]

  # Remove zero columns
  for (col in prop_cols) {
    if (col %in% names(composition) && all(composition[[col]] <= epsilon * 2)) {
      composition[[col]] <- NULL
    }
  }

  # CRITICAL FIX: Re-normalize remaining columns to sum to 1
  # After removing zero columns, proportions no longer sum to 1
  if (length(prop_cols_to_keep) > 0) {
    final_row_sums <- rowSums(composition[, prop_cols_to_keep, drop = FALSE])
    # Avoid division by zero
    final_row_sums[final_row_sums == 0] <- 1
    for (col in prop_cols_to_keep) {
      if (col %in% names(composition)) {
        composition[[col]] <- composition[[col]] / final_row_sums
      }
    }
  }

  # Verify proportions sum to 1 (within tolerance)
  if (length(prop_cols_to_keep) > 0) {
    test_sums <- rowSums(composition[, prop_cols_to_keep, drop = FALSE])
    if (any(abs(test_sums - 1) > 0.01)) {
      warning(sprintf("Some rows don't sum to 1 after normalization. Range: [%.4f, %.4f]",
                      min(test_sums), max(test_sums)))
    }
  }

  # Create REEF and HABMERGED
  composition$REEF <- factor(substr(as.character(composition$SITE), 1, 3))
  composition$HABMERGED <- ifelse(composition$HAB %in% c("RR", "TP"), "RR_TP", composition$HAB)
  composition$HABMERGED <- factor(composition$HABMERGED)
  composition$HABMERGED <- droplevels(composition$HABMERGED)

  if (length(unique(composition$HABMERGED)) < 2) {
    composition$HABMERGED <- NULL
  }

  composition$ARCH <- factor(composition$ARCH)
  if (length(unique(composition$ARCH)) < 2) {
    composition$ARCH <- NULL
  }

  # Scale predictors (remove underscores for brms)
  for (col in c("PC1_MAGNITUDE", "PC2_MAGNITUDE",
                "PC1_VARIABILITY", "PC2_VARIABILITY", "DEPTH_M")) {
    if (col %in% names(composition)) {
      scaled_name <- gsub("_", "", col)
      composition[[scaled_name]] <- scale(as.numeric(composition[[col]]))[, 1]
    }
  }

  # Convert YEAR to factor for random effect
  composition$YEAR <- factor(composition$YEAR)

  # Metadata
  composition$CVLABEL <- cv_label

  # Validate YEAR
  n_years <- length(unique(composition$YEAR))
  cat(sprintf("  JSDM data has %d YEAR levels: %s\n",
              n_years, paste(unique(composition$YEAR), collapse = ", ")))

  cat(sprintf("✓ JSDM data prepared for %s: N=%d sampling units, %d REEFs\n",
              cv_label, nrow(composition), length(unique(composition$REEF))))

  return(as.data.frame(composition))
}

# ============================================================================
# SECTION 4: FORMULA BUILDERS WITH CROSSED RANDOM EFFECTS
# ============================================================================

#' Build formula for Gaussian Health models with CROSSED YEAR random effect
#' @param response_var Response variable ("HEALTH_PC1", "HEALTH_PC2", "RGR")
#' @param include_hab Include HABMERGED?
#' @param include_depth Include s(DEPTHM)?
#' @param include_interaction_pca Include interaction PCA terms? (NEW)
#' @return formula object for brms
make_gaussian_formula_year_re <- function(response_var, include_hab = TRUE,
                                         include_depth = TRUE,
                                         include_interaction_pca = TRUE) {
  k_val <- 3

  # Fixed effects (same as original)
  main_part <- paste0(response_var, " ~ PC1MAGNITUDE + PC2MAGNITUDE")
  main_part <- paste0(main_part, " + s(PC1MAGNITUDE, k = ", k_val, ")")
  main_part <- paste0(main_part, " + s(PC2MAGNITUDE, k = ", k_val, ")")
  main_part <- paste0(main_part, " + s(PC1VARIABILITY, k = ", k_val, ")")
  main_part <- paste0(main_part, " + s(PC2VARIABILITY, k = ", k_val, ")")

  # NEW: Add interaction PCA terms as splines
  if (include_interaction_pca) {
    main_part <- paste0(main_part, " + s(PC1INTERACAO, k = ", k_val, ")")
    main_part <- paste0(main_part, " + s(PC2INTERACAO, k = ", k_val, ")")
  }

  if (include_hab) {
    main_part <- paste0(main_part, " + HABMERGED")
  }

  if (include_depth) {
    main_part <- paste0(main_part, " + s(DEPTHM, k = ", k_val, ")")
  }

  # CROSSED random effects (not nested!)
  # (1|REEF) + (1|YEAR) means YEAR is crossed with REEF
  # NOT (1|REEF/YEAR) which would be nested
  main_part <- paste0(main_part, " + (1 | REEF) + (1 | YEAR)")

  stats::as.formula(main_part)
}


#' Build formula for ZOIB models with CROSSED YEAR random effect
#' @param include_hab Include HABMERGED?
#' @param include_depth Include s(DEPTHM)?
#' @param include_arch_interaction Include ARCH interactions?
#' @return brmsformula object for brms
make_zoib_formula_year_re <- function(include_hab = TRUE, include_depth = TRUE,
                                     include_arch_interaction = TRUE) {
  k_val <- 3

  main_part <- "COVER ~ PC1MAGNITUDE + PC2MAGNITUDE"
  main_part <- paste0(main_part, " + s(PC1VARIABILITY, k = ", k_val, ")")
  main_part <- paste0(main_part, " + s(PC2VARIABILITY, k = ", k_val, ")")

  if (include_arch_interaction) {
    main_part <- paste0(main_part, " + PC1VARIABILITY:ARCH")
    main_part <- paste0(main_part, " + PC2VARIABILITY:ARCH")
  }

  if (include_hab) {
    main_part <- paste0(main_part, " + HABMERGED")
  }

  if (include_depth) {
    main_part <- paste0(main_part, " + s(DEPTHM, k = ", k_val, ")")
  }

  # CROSSED random effects (not nested!)
  main_part <- paste0(main_part, " + (1 | REEF) + (1 | YEAR)")

  main_formula <- stats::as.formula(main_part)

  return(main_formula)
}


#' Build formula for JSDM Dirichlet models with CROSSED YEAR random effect
#' @param include_hab Include HABMERGED?
#' @param include_depth Include DEPTHM?
#' @param include_arch_interaction Include ARCH interactions?
#' @param data Data frame to check for column availability
#' @return bf() object for brms with Dirichlet
make_jsdm_formula_year_re <- function(include_hab = TRUE, include_depth = TRUE,
                                     include_arch_interaction = TRUE, data = NULL) {

  # Check column availability
  has_habmerged <- if (!is.null(data) && "HABMERGED" %in% names(data)) {
    length(unique(data$HABMERGED)) >= 2
  } else {
    FALSE
  }
  has_arch <- if (!is.null(data) && "ARCH" %in% names(data)) {
    length(unique(data$ARCH)) >= 2
  } else {
    FALSE
  }

  # Fixed effects (linear for JSDM)
  mag_terms <- "PC1MAGNITUDE + PC2MAGNITUDE"
  var_terms <- "PC1VARIABILITY + PC2VARIABILITY"

  arch_terms <- if (include_arch_interaction && has_arch) {
    " + PC1VARIABILITY:ARCH + PC2VARIABILITY:ARCH"
  } else ""

  hab_term <- if (include_hab && has_habmerged) " + HABMERGED" else ""
  depth_term <- if (include_depth) " + DEPTHM" else ""

  # Response columns
  response_cols <- c("MUSSISMILIA_prop", "TURF_prop", "CCA_prop",
                     "CYANO_prop", "MACROALGAE_prop")
  if (!is.null(data)) {
    response_cols <- intersect(response_cols, names(data))
  }
  response_str <- paste0("cbind(", paste(response_cols, collapse = ", "), ")")

  # CROSSED random effects (not nested!)
  formula_str <- paste0(
    response_str, " ~ ",
    mag_terms, " + ", var_terms, arch_terms, hab_term, depth_term,
    " + (1 | REEF) + (1 | YEAR)"
  )

  bf(stats::as.formula(formula_str))
}

# ============================================================================
# SECTION 5: MODEL COMBINATIONS DEFINITION (YEAR RE SPECIFIC)
# ============================================================================

# Model combinations for YEAR RE models
model_combinations_year_re <- list(
  # Existing models (without interaction PCA) - KEEP AS IS
  list(
    name = "FULL",
    include_hab = TRUE,
    include_depth = TRUE,
    include_arch_interaction = TRUE,
    include_interaction_pca = FALSE,  # NEW: without interaction
    description = "Full model with HAB, DEPTH (no Interaction PCA)"
  ),
  list(
    name = "NOHABITAT",
    include_hab = FALSE,
    include_depth = TRUE,
    include_arch_interaction = TRUE,
    include_interaction_pca = FALSE,
    description = "No HAB, with DEPTH (no Interaction PCA)"
  ),
  list(
    name = "NODEPTH",
    include_hab = TRUE,
    include_depth = FALSE,
    include_arch_interaction = TRUE,
    include_interaction_pca = FALSE,
    description = "With HAB, no DEPTH (no Interaction PCA)"
  ),
  list(
    name = "MINIMAL",
    include_hab = FALSE,
    include_depth = FALSE,
    include_arch_interaction = FALSE,
    include_interaction_pca = FALSE,
    description = "Minimal model (PCA only, no Interaction PCA)"
  ),

  # NEW: Models with interaction PCA
  list(
    name = "FULL_INTERACTION",
    include_hab = TRUE,
    include_depth = TRUE,
    include_arch_interaction = TRUE,
    include_interaction_pca = TRUE,  # NEW: with interaction
    description = "Full model with HAB, DEPTH, and Interaction PCA"
  ),
  list(
    name = "NOHABITAT_INTERACTION",
    include_hab = FALSE,
    include_depth = TRUE,
    include_arch_interaction = TRUE,
    include_interaction_pca = TRUE,
    description = "No HAB, with DEPTH and Interaction PCA"
  ),
  list(
    name = "NODEPTH_INTERACTION",
    include_hab = TRUE,
    include_depth = FALSE,
    include_arch_interaction = TRUE,
    include_interaction_pca = TRUE,
    description = "With HAB and Interaction PCA, no DEPTH"
  )
)

# ============================================================================
# SECTION 6: MODEL FITTING WITH CACHE (YEAR RE)
# ============================================================================

#' Load model from cache or fit new model - YEAR RE version
#' @param model_name Model name (e.g., "FULL")
#' @param cv_name CV name (e.g., "CV_02")
#' @param model_type Model type (e.g., "ZOIB_YEAR_RE", "HEALTH_PC1")
#' @param formula Model formula
#' @param data Prepared data
#' @param family Model family
#' @param prior Priors
#' @param out_dir Output directory
#' @param brms_args Additional arguments for brm()
#' @return List with fit and cache flag
load_or_fit_model_year_re <- function(model_name, cv_name, model_type,
                                      formula, data, family, prior,
                                      out_dir, brms_args) {

  # File name - YEAR RE models include YEAR_RE in type name
  model_filename <- sprintf("%s_%s_%s.rds",
                            tolower(model_type), tolower(model_name), cv_name)
  model_path <- file.path(out_dir, model_filename)

  # Check cache
  if (file.exists(model_path)) {
    cat(sprintf("  📦 Cache found: %s\n", model_filename))
    fit <- readRDS(model_path)
    return(list(fit = fit, cached = TRUE, path = model_path))
  }

  # Fit new model
  cat(sprintf("  🔧 Fitting model: %s\n", model_filename))
  start_time <- Sys.time()

  fit <- do.call(brm, c(
    list(
      formula = formula,
      data = data,
      family = family,
      prior = prior
    ),
    brms_args
  ))

  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")
  cat(sprintf("  ✓ Model fitted in %.1f minutes\n", as.numeric(elapsed)))

  # Save immediately
  saveRDS(fit, model_path)
  cat(sprintf("  💾 Model saved: %s\n", model_path))

  return(list(fit = fit, cached = FALSE, path = model_path,
              elapsed_mins = as.numeric(elapsed)))
}


# ============================================================================
# SECTION 7: CONVERGENCE CHECKING
# ============================================================================

#' Check convergence criteria for a brms model - YEAR RE compatible
#' @param fit brmsfit object
#' @param model_name Model name for logging
#' @param cv_label CV label for logging
#' @return List with convergence status and metrics
check_convergence <- function(fit, model_name, cv_label) {

  # Extract summary
  summary_fit <- summary(fit)

  # --- R-hat ---
  rhat_vals <- c(summary_fit$fixed[, "Rhat"],
                 summary_fit$random$REEF[, "Rhat"])
  # Add YEAR random effect R-hat if present
  if ("YEAR" %in% names(summary_fit$random)) {
    rhat_vals <- c(rhat_vals, summary_fit$random$YEAR[, "Rhat"])
  }
  rhat_vals <- rhat_vals[!is.na(rhat_vals)]
  rhat_max <- max(rhat_vals)
  rhat_ok <- rhat_max < 1.01

  # --- ESS ---
  ess_bulk_vals <- c(summary_fit$fixed[, "Bulk_ESS"],
                     summary_fit$random$REEF[, "Bulk_ESS"])
  ess_tail_vals <- c(summary_fit$fixed[, "Tail_ESS"],
                     summary_fit$random$REEF[, "Tail_ESS"])

  # Add YEAR random effect ESS if present
  if ("YEAR" %in% names(summary_fit$random)) {
    ess_bulk_vals <- c(ess_bulk_vals, summary_fit$random$YEAR[, "Bulk_ESS"])
    ess_tail_vals <- c(ess_tail_vals, summary_fit$random$YEAR[, "Tail_ESS"])
  }

  ess_bulk_vals <- ess_bulk_vals[!is.na(ess_bulk_vals)]
  ess_tail_vals <- ess_tail_vals[!is.na(ess_tail_vals)]

  ess_bulk_min <- min(ess_bulk_vals)
  ess_tail_min <- min(ess_tail_vals)
  ess_ok <- (ess_bulk_min > 400) && (ess_tail_min > 400)

  # --- Divergences ---
  n_div <- 0
  tryCatch({
    np <- nuts_params(fit)
    n_div <- sum(np$Value[np$Parameter == "divergent__"])
  }, error = function(e) {
    n_div <- NA
  })
  div_ok <- is.na(n_div) || (n_div == 0)

  # Overall result
  all_ok <- rhat_ok && ess_ok && div_ok

  # Logging
  status <- if (all_ok) "  ✓ CONVERGED" else "  ⚠ PROBLEMS"
  cat(sprintf("%s %s_%s: Rhat=%.4f, ESS_bulk=%d, ESS_tail=%d, Div=%s\n",
              status, model_name, cv_label, rhat_max,
              round(ess_bulk_min), round(ess_tail_min),
              ifelse(is.na(n_div), "NA", as.character(n_div))))

  return(list(
    passed = all_ok,
    rhat_max = rhat_max,
    ess_bulk_min = ess_bulk_min,
    ess_tail_min = ess_tail_min,
    n_divergent = n_div,
    details = list(rhat_ok = rhat_ok, ess_ok = ess_ok, div_ok = div_ok)
  ))
}


# ============================================================================
# SECTION 8: ENHANCED CONVERGENCE DIAGNOSTICS WITH ADAPTIVE DELTA RECOMMENDATIONS
# ============================================================================

#' Comprehensive convergence diagnostics with adaptive delta recommendations
#' @param fit brmsfit object
#' @param model_name Model name for logging
#' @param cv_label CV label for logging
#' @param current_adapt_delta Current adapt_delta value used
#' @param current_max_treedepth Current max_treedepth value used
#' @return List with detailed diagnostics and recommendations
diagnose_convergence_with_recommendations <- function(fit, model_name, cv_label,
                                                      current_adapt_delta = 0.97,
                                                      current_max_treedepth = 12) {

  cat(sprintf("\n=== DETAILED CONVERGENCE DIAGNOSTICS: %s_%s ===\n", model_name, cv_label))

  # Extract summary and diagnostics
  summary_fit <- summary(fit)

  # Initialize results list
  diagnostics <- list(
    model_name = model_name,
    cv_label = cv_label,
    current_settings = list(adapt_delta = current_adapt_delta, max_treedepth = current_max_treedepth),
    issues = list(),
    warnings = list(),
    recommendations = list()
  )

  # --- 1. R-hat Diagnostics ---
  rhat_all <- c()
  if (!is.null(summary_fit$fixed)) {
    rhat_all <- c(rhat_all, summary_fit$fixed[, "Rhat"])
  }
  if (!is.null(summary_fit$random) && "REEF" %in% names(summary_fit$random)) {
    rhat_all <- c(rhat_all, summary_fit$random$REEF[, "Rhat"])
  }
  if (!is.null(summary_fit$random) && "YEAR" %in% names(summary_fit$random)) {
    rhat_all <- c(rhat_all, summary_fit$random$YEAR[, "Rhat"])
  }
  rhat_all <- rhat_all[!is.na(rhat_all)]
  rhat_max <- max(rhat_all)
  rhat_n_fail <- sum(rhat_all > 1.01)

  diagnostics$metrics$rhat <- list(max = rhat_max, n_above_threshold = rhat_n_fail)

  if (rhat_max > 1.05) {
    diagnostics$issues$rhat <- sprintf("CRITICAL: Rhat max = %.4f (> 1.05). %d parameters above 1.01",
                                       rhat_max, rhat_n_fail)
  } else if (rhat_max > 1.01) {
    diagnostics$warnings$rhat <- sprintf("WARNING: Rhat max = %.4f (> 1.01). %d parameters affected",
                                         rhat_max, rhat_n_fail)
  } else {
    diagnostics$metrics$rhat$status <- "PASS"
  }

  # --- 2. ESS Diagnostics ---
  ess_bulk_all <- c()
  ess_tail_all <- c()
  if (!is.null(summary_fit$fixed)) {
    ess_bulk_all <- c(ess_bulk_all, summary_fit$fixed[, "Bulk_ESS"])
    ess_tail_all <- c(ess_tail_all, summary_fit$fixed[, "Tail_ESS"])
  }
  if (!is.null(summary_fit$random) && "REEF" %in% names(summary_fit$random)) {
    ess_bulk_all <- c(ess_bulk_all, summary_fit$random$REEF[, "Bulk_ESS"])
    ess_tail_all <- c(ess_tail_all, summary_fit$random$REEF[, "Tail_ESS"])
  }
  if (!is.null(summary_fit$random) && "YEAR" %in% names(summary_fit$random)) {
    ess_bulk_all <- c(ess_bulk_all, summary_fit$random$YEAR[, "Bulk_ESS"])
    ess_tail_all <- c(ess_tail_all, summary_fit$random$YEAR[, "Tail_ESS"])
  }

  ess_bulk_all <- ess_bulk_all[!is.na(ess_bulk_all)]
  ess_tail_all <- ess_tail_all[!is.na(ess_tail_all)]

  ess_bulk_min <- min(ess_bulk_all)
  ess_tail_min <- min(ess_tail_all)
  ess_bulk_mean <- mean(ess_bulk_all)
  ess_tail_mean <- mean(ess_tail_all)

  diagnostics$metrics$ess <- list(
    bulk_min = ess_bulk_min, bulk_mean = ess_bulk_mean,
    tail_min = ess_tail_min, tail_mean = ess_tail_mean
  )

  if (ess_bulk_min < 100 || ess_tail_min < 100) {
    diagnostics$issues$ess <- sprintf("CRITICAL: ESS very low (Bulk: %d, Tail: %d)",
                                      ess_bulk_min, ess_tail_min)
  } else if (ess_bulk_min < 400 || ess_tail_min < 400) {
    diagnostics$warnings$ess <- sprintf("WARNING: ESS below recommended (Bulk: %d, Tail: %d < 400)",
                                       ess_bulk_min, ess_tail_min)
  } else {
    diagnostics$metrics$ess$status <- "PASS"
  }

  # --- 3. Divergent Transitions ---
  n_div <- 0
  div_per_chain <- NULL
  tryCatch({
    np <- nuts_params(fit)
    n_div <- sum(np$Value[np$Parameter == "divergent__"])

    # Divergences per chain
    div_by_chain <- tapply(np$Value[np$Parameter == "divergent__"],
                          np$Chain[np$Parameter == "divergent__"], sum)
    div_per_chain <- as.numeric(div_by_chain)
    names(div_per_chain) <- paste("Chain", 1:length(div_per_chain))
  }, error = function(e) {
    n_div <- NA
  })

  diagnostics$metrics$divergences <- list(n = n_div, per_chain = div_per_chain)

  if (!is.na(n_div) && n_div > 100) {
    diagnostics$issues$divergences <- sprintf("CRITICAL: %d divergent transitions detected", n_div)
  } else if (!is.na(n_div) && n_div > 0) {
    diagnostics$warnings$divergences <- sprintf("WARNING: %d divergent transitions detected", n_div)
  } else if (is.na(n_div)) {
    diagnostics$warnings$divergences <- "Could not assess divergences"
  } else {
    diagnostics$metrics$divergences$status <- "PASS"
  }

  # --- 4. BFMI (Bayesian Fraction of Missing Information) ---
  bfmi_values <- NULL
  tryCatch({
    np <- nuts_params(fit)
    bfmi_vals <- unique(np$Value[np$Parameter == "energy__"])
    # Calculate BFMI from energy - use per-chain values
    chain_ids <- unique(np$Chain[np$Parameter == "energy__"])
    bfmi_per_chain <- sapply(chain_ids, function(ch) {
      energies <- np$Value[np$Parameter == "energy__" & np$Chain == ch]
      # BFMI approximation: var(energies) / mean((diff(energies))^2)
      n <- length(energies)
      if (n > 1) {
        variance <- var(energies)
        num_integrals <- sum((diff(energies))^2) / (n - 1)
        return(variance / num_integrals)
      }
      return(NA)
    })
    bfmi_values <- bfmi_per_chain
  }, error = function(e) {
    bfmi_values <- NULL
  })

  diagnostics$metrics$bfmi <- list(values = bfmi_values)

  if (!is.null(bfmi_values)) {
    bfmi_min <- min(bfmi_values)
    if (bfmi_min < 0.2) {
      diagnostics$warnings$bfmi <- sprintf("WARNING: Low BFMI detected (min: %.3f < 0.2)", bfmi_min)
    } else {
      diagnostics$metrics$bfmi$status <- "PASS"
    }
  }

  # --- 5. Tree Depth ---
  max_td_hit <- FALSE
  tryCatch({
    np <- nuts_params(fit)
    max_td_values <- unique(np$Value[np$Parameter == "treedepth__"])
    if (any(max_td_values >= current_max_treedepth - 1)) {
      max_td_hit <- TRUE
      diagnostics$warnings$treedepth <- sprintf("WARNING: Some chains hit max_treedepth (%d)", current_max_treedepth)
    }
  }, error = function(e) {
    # Cannot assess treedepth
  })

  # --- GENERATE RECOMMENDATIONS ---
  recommendations <- list()

  # Priority 1: Divergences -> Increase adapt_delta
  if (!is.na(n_div) && n_div > 0) {
    # Calculate recommended adapt_delta
    if (n_div > 100) {
      # Many divergences - big jump
      new_delta <- min(0.999, current_adapt_delta + 0.02)
    } else if (n_div > 10) {
      # Some divergences - moderate increase
      new_delta <- min(0.995, current_adapt_delta + 0.01)
    } else {
      # Few divergences - small increase
      new_delta <- min(0.99, current_adapt_delta + 0.005)
    }
    recommendations$adapt_delta <- list(
      current = current_adapt_delta,
      recommended = new_delta,
      reason = sprintf("%d divergent transitions detected", n_div),
      priority = "HIGH"
    )
  }

  # Priority 2: Max treedepth hit
  if (max_td_hit) {
    recommendations$max_treedepth <- list(
      current = current_max_treedepth,
      recommended = current_max_treedepth + 1,
      reason = "Chains hitting maximum tree depth",
      priority = "MEDIUM"
    )
  }

  # Priority 3: Low ESS -> Increase iterations
  if (ess_bulk_min < 400 || ess_tail_min < 400) {
    # Calculate multiplier needed
    multiplier <- ceiling(400 / min(ess_bulk_min, ess_tail_min))
    recommendations$iterations <- list(
      current_iter = 4000,
      recommended_iter = 4000 * multiplier,
      reason = sprintf("Low ESS (Bulk: %d, Tail: %d)", ess_bulk_min, ess_tail_min),
      priority = "MEDIUM"
    )
  }

  # Priority 4: High Rhat -> Usually requires model reformulation
  if (rhat_max > 1.05) {
    recommendations$model_reformulation <- list(
      reason = "Very high Rhat suggests fundamental convergence issues",
      suggestions = c(
        "Check for collinearity among predictors",
        "Simplify model structure (remove interactions)",
        "Check data separation for ZOIB models",
        "Consider stronger priors"
      ),
      priority = "HIGH"
    )
  }

  diagnostics$recommendations <- recommendations

  # --- SUMMARY AND STATUS ---
  n_issues <- length(diagnostics$issues)
  n_warnings <- length(diagnostics$warnings)

  if (n_issues == 0 && n_warnings == 0) {
    diagnostics$status <- "EXCELLENT"
    diagnostics$passed <- TRUE
    cat(sprintf("  ✓✓✓ EXCELLENT CONVERGENCE: Rhat=%.4f, ESS_bulk=%d, ESS_tail=%d, Div=0\n",
                rhat_max, ess_bulk_min, ess_tail_min))
  } else if (n_issues == 0) {
    diagnostics$status <- "GOOD"
    diagnostics$passed <- TRUE
    cat(sprintf("  ✓✓ GOOD CONVERGENCE: Rhat=%.4f, ESS_bulk=%d, ESS_tail=%d, Div=%d\n",
                rhat_max, ess_bulk_min, ess_tail_min, n_div))
    if (n_warnings > 0) {
      cat(sprintf("     Minor warnings: %s\n", paste(names(diagnostics$warnings), collapse = ", ")))
    }
  } else {
    diagnostics$status <- "PROBLEMS"
    diagnostics$passed <- FALSE
    cat(sprintf("  ⚠⚠⚠ CONVERGENCE PROBLEMS DETECTED:\n"))
    for (issue_name in names(diagnostics$issues)) {
      cat(sprintf("     - %s\n", diagnostics$issues[[issue_name]]))
    }
  }

  # --- PRINT RECOMMENDATIONS ---
  if (length(recommendations) > 0) {
    cat(sprintf("\n  🔧 RECOMMENDED ACTIONS:\n"))
    for (rec_name in names(recommendations)) {
      rec <- recommendations[[rec_name]]
      if (rec_name == "adapt_delta") {
        cat(sprintf("     [%s] Set adapt_delta = %.3f (current: %.3f) - %s\n",
                    rec$priority, rec$recommended, rec$current, rec$reason))
      } else if (rec_name == "max_treedepth") {
        cat(sprintf("     [%s] Set max_treedepth = %d (current: %d) - %s\n",
                    rec$priority, rec$recommended, rec$current, rec$reason))
      } else if (rec_name == "iterations") {
        cat(sprintf("     [%s] Increase iter to %d (current: %d) - %s\n",
                    rec$priority, rec$recommended_iter, rec$current_iter, rec$reason))
      } else if (rec_name == "model_reformulation") {
        cat(sprintf("     [%s] Model reformulation needed:\n", rec$priority))
        for (s in rec$suggestions) {
          cat(sprintf("         • %s\n", s))
        }
      }
    }
  }

  return(diagnostics)
}


#' Generate updated brms_args based on convergence diagnostics
#' @param base_args Base brms_args list (e.g., brms_args_year_re)
#' @param diagnostics Output from diagnose_convergence_with_recommendations()
#' @return Updated brms_args list
generate_updated_brms_args <- function(base_args, diagnostics) {

  new_args <- base_args

  if (is.null(diagnostics$recommendations) || length(diagnostics$recommendations) == 0) {
    return(new_args)
  }

  # Update adapt_delta if recommended
  if (!is.null(diagnostics$recommendations$adapt_delta)) {
    rec <- diagnostics$recommendations$adapt_delta
    new_args$control$adapt_delta <- rec$recommended
    cat(sprintf("  Updated adapt_delta: %.3f -> %.3f\n",
                base_args$control$adapt_delta, rec$recommended))
  }

  # Update max_treedepth if recommended
  if (!is.null(diagnostics$recommendations$max_treedepth)) {
    rec <- diagnostics$recommendations$max_treedepth
    new_args$control$max_treedepth <- rec$recommended
    cat(sprintf("  Updated max_treedepth: %d -> %d\n",
                base_args$control$max_treedepth, rec$recommended))
  }

  # Update iterations if recommended
  if (!is.null(diagnostics$recommendations$iterations)) {
    rec <- diagnostics$recommendations$iterations
    new_args$iter <- rec$recommended_iter
    new_args$warmup <- rec$recommended_iter / 2
    cat(sprintf("  Updated iterations: %d -> %d\n",
                base_args$iter, rec$recommended_iter))
  }

  return(new_args)
}


#' Auto-refit model with improved settings if convergence failed
#' @param model_name Model name
#' @param cv_name CV name
#' @param model_type Model type
#' @param formula Model formula
#' @param data Prepared data
#' @param family Model family
#' @param prior Priors
#' @param out_dir Output directory
#' @param base_brms_args Base brms args
#' @param max_refits Maximum number of refit attempts (default: 2)
#' @return List with fit, diagnostics, and refit history
auto_fit_with_convergence_check <- function(model_name, cv_name, model_type,
                                            formula, data, family, prior,
                                            out_dir, base_brms_args,
                                            max_refits = 2) {

  refit_history <- list()
  current_args <- base_brms_args

  for (attempt in 0:max_refits) {
    attempt_suffix <- if (attempt == 0) "" else sprintf("_attempt%d", attempt)

    # File name for this attempt
    model_filename <- sprintf("%s_%s_%s%s.rds",
                              tolower(model_type), tolower(model_name),
                              cv_name, attempt_suffix)
    model_path <- file.path(out_dir, model_filename)

    # Check cache first
    if (file.exists(model_path) && attempt == 0) {
      cat(sprintf("  📦 Cache found: %s\n", model_filename))
      fit <- readRDS(model_path)
      refit_history$cached <- TRUE
    } else {
      # Fit with current settings
      cat(sprintf("\n  🔧 Fitting attempt %d for %s_%s\n", attempt, model_name, cv_name))
      if (attempt > 0) {
        cat(sprintf("     adapt_delta=%.3f, max_treedepth=%d\n",
                    current_args$control$adapt_delta,
                    current_args$control$max_treedepth))
      }

      start_time <- Sys.time()

      fit <- tryCatch({
        do.call(brm, c(
          list(formula = formula, data = data, family = family, prior = prior),
          current_args
        ))
      }, error = function(e) {
        cat(sprintf("  ❌ Fitting error: %s\n", e$message))
        return(NULL)
      })

      if (is.null(fit)) {
        refit_history[[paste0("attempt_", attempt)]] <- list(
          success = FALSE,
          error = "Fitting failed"
        )
        next
      }

      end_time <- Sys.time()
      elapsed <- as.numeric(difftime(end_time, start_time, units = "mins"))

      # Save model
      saveRDS(fit, model_path)
      cat(sprintf("  💾 Model saved: %s (%.1f min)\n", model_filename, elapsed))
    }

    # Run diagnostics
    current_delta <- current_args$control$adapt_delta
    current_td <- current_args$control$max_treedepth

    diags <- diagnose_convergence_with_recommendations(
      fit, model_name, cv_name,
      current_adapt_delta = current_delta,
      current_max_treedepth = current_td
    )

    refit_history[[paste0("attempt_", attempt)]] <- list(
      success = TRUE,
      diagnostics = diags,
      settings = list(adapt_delta = current_delta, max_treedepth = current_td)
    )

    # Check if convergence achieved
    if (diags$passed) {
      cat(sprintf("\n  ✅ Convergence achieved on attempt %d\n", attempt))
      return(list(
        fit = fit,
        diagnostics = diags,
        refit_history = refit_history,
        final_attempt = attempt
      ))
    }

    # Prepare for refit if not at max attempts
    if (attempt < max_refits) {
      cat(sprintf("\n  🔄 Preparing refit attempt %d...\n", attempt + 1))
      current_args <- generate_updated_brms_args(current_args, diags)
    } else {
      cat(sprintf("\n  ⚠ Maximum refit attempts (%d) reached. Model may have convergence issues.\n",
                  max_refits))
    }
  }

  # Return final fit even if convergence issues remain
  return(list(
    fit = fit,
    diagnostics = diags,
    refit_history = refit_history,
    final_attempt = max_refits
  ))
}

# ============================================================================
# SECTION 9: LOO CALCULATION AND COMPARISON
# ============================================================================

#' Calculate LOO with error handling and moment matching
#' @param fit brmsfit object
#' @param model_name Model name
#' @param cv_label CV label
#' @param use_moment_match Try moment_match for high Pareto k?
#' @return loo object or NULL on error
safe_loo <- function(fit, model_name, cv_label, use_moment_match = TRUE) {

  cat(sprintf("  📊 Calculating LOO for %s_%s...\n", model_name, cv_label))

  loo_result <- tryCatch({
    loo_obj <- loo(fit, cores = 4)

    # Check Pareto k
    k_vals <- loo_obj$diagnostics$pareto_k
    n_high_k <- sum(k_vals > 0.7, na.rm = TRUE)

    if (n_high_k > 0 && use_moment_match) {
      cat(sprintf("  ⚠ %d observations with Pareto k > 0.7. Trying moment_match...\n", n_high_k))

      loo_obj <- tryCatch({
        loo(fit, cores = 4, moment_match = TRUE, k_threshold = 0.7)
      }, error = function(e) {
        cat("  ⚠ moment_match failed. Using original LOO.\n")
        loo_obj
      })
    }

    loo_obj

  }, error = function(e) {
    cat(sprintf("  ❌ LOO error for %s_%s: %s\n", model_name, cv_label, e$message))
    NULL
  })

  return(loo_result)
}


#' Compare models LOO within a CV
#' @param loo_list Named list of loo objects
#' @param cv_label CV label
#' @return data.frame with comparison and winner model
compare_loo_within_cv <- function(loo_list, cv_label) {

  # Remove NULLs
  loo_list <- Filter(Negate(is.null), loo_list)

  if (length(loo_list) < 2) {
    cat(sprintf("  ⚠ Less than 2 valid models for comparison in %s\n", cv_label))
    return(NULL)
  }

  # Compare
  comparison <- loo_compare(loo_list)

  # Extract results
  comp_df <- as.data.frame(comparison)
  comp_df$Model <- rownames(comp_df)
  comp_df$CV <- cv_label
  comp_df$Rank <- 1:nrow(comp_df)

  # Determine winner
  winner <- rownames(comparison)[1]
  delta_second <- abs(comparison[2, "elpd_diff"]) * 2  # ΔLOOIC = 2 * Δelpd
  se_second <- comparison[2, "se_diff"]

  # Significance
  significance <- if (requireNamespace("dplyr", quietly = TRUE)) {
    dplyr::case_when(
      delta_second > 2 * se_second ~ "Significant",
      delta_second > se_second ~ "Moderate",
      TRUE ~ "Not significant"
    )
  } else {
    if (delta_second > 2 * se_second) "Significant"
    else if (delta_second > se_second) "Moderate"
    else "Not significant"
  }

  cat(sprintf("  🏆 Winner for %s: %s (Δ=%.1f, SE=%.1f, %s)\n",
              cv_label, winner, delta_second, se_second, significance))

  return(list(
    comparison = comp_df,
    winner = winner,
    delta_looic_to_second = delta_second,
    se_diff = se_second,
    significance = significance
  ))
}


# ============================================================================
# SECTION 10: PRIORS (same as original)
# ============================================================================

# Priors for ZOIB (Abundance) - YEAR RE
# Note: Removed dpar priors that cause issues with auto-generated formulas
priors_zoib_year_re <- c(
  prior(normal(0, 1), class = "sds"),
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(2), class = "sd"),
  prior(normal(0, 2), class = "Intercept")
)

# Priors for Gaussian (Health/RGR) - YEAR RE
priors_gaussian_year_re <- c(
  prior(normal(0, 1), class = "sds"),
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(2), class = "sd"),
  prior(normal(0, 2), class = "Intercept")
)

# Priors for JSDM Dirichlet - YEAR RE
# Using NULL to use brms defaults (recommended for Dirichlet)
priors_jsdm_year_re <- NULL

# Standard brms configuration for YEAR RE
# IMPROVED PARAMETERS (Task 1): Better convergence for Health models
# iter: 8000 (was 4000), warmup: 4000 (was 2000), adapt_delta: 0.99 (was 0.97)
brms_args_year_re <- list(
  backend = "cmdstanr",
  cores = 4,
  iter = 8000,
  warmup = 4000,
  chains = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 12),
  refresh = 500,
  save_pars = save_pars(all = TRUE)
)

# More conservative configuration for JSDM YEAR RE
brms_args_jsdm_year_re <- list(
  backend = "cmdstanr",
  cores = 4,
  iter = 6000,
  warmup = 3000,
  chains = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),
  refresh = 500,
  save_pars = save_pars(all = TRUE)
)

# ============================================================================
# SECTION 11: PROJECT PATHS (for YEAR RE)
# ============================================================================

# Project root (same as original)
PROJECT_ROOT <- "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/"
RESULTS_DIR <- file.path(PROJECT_ROOT, "#######FINAL_RESULTS/")

# Abundance dataset paths (for ZOIB and JSDM) - contain YEAR
dataset_paths_abundance <- list(
  CV_02 = file.path(RESULTS_DIR, "#####output_local_PCA_CV_2_FINAL/dados_abundancia_integrados_long_format.csv"),
  CV_30 = file.path(RESULTS_DIR, "#####output_local_PCA_CV_30_FINAL/dados_abundancia_integrados_long_format.csv"),
  CV_ALL = file.path(RESULTS_DIR, "#####output_local_PCA_CV_all_FINAL/dados_abundancia_integrados_long_format.csv")
)

# ============================================================================
# END OF 00_LOO_Selection_Functions_YEAR_RE.R
# Added:
#   - load_or_fit_model_year_re(), safe_loo(), check_convergence()
#   - diagnose_convergence_with_recommendations()
#   - generate_updated_brms_args()
#   - auto_fit_with_convergence_check()
# ============================================================================
