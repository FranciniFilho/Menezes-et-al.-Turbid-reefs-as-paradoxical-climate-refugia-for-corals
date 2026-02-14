# ============================================================================
# 00_LOO_Selection_Functions.R
# ============================================================================
# Common helper functions for Full LOO Model Selection
# Unified implementation for ZOIB, Gaussian, and JSDM models
# ============================================================================
# Author: Auto-generated from UNIFIED_MODEL_SELECTION_LOO_IMPLEMENTATION_v1.md
# Date: 2026-02-04
# ============================================================================

# ============================================================================
# SECTION 1: DATA PREPARATION FUNCTIONS
# ============================================================================

#' Prepare abundance data for ZOIB models
#' @param data_path Path to CSV file
#' @param cv_label CV scenario label ("CV_02", "CV_30", "CV_ALL")
#' @return data.frame prepared for brms
prepare_zoib_data <- function(data_path, cv_label) {

  # --- 1. READ WITH DELIMITER DETECTION ---
  data <- tryCatch({
    read.csv2(data_path, stringsAsFactors = TRUE)
  }, error = function(e) {
    read.csv(data_path, stringsAsFactors = TRUE)
  })

  # --- 2. STANDARDIZE COLUMN NAMES ---
  colnames(data) <- toupper(colnames(data))

  # --- 3. FILTER FOR M. HISPIDA ---
  data <- data[data$ORGANISMO == "MUSSISMILIA_HISPIDA", ]

  if (nrow(data) == 0) {
    stop(paste("No M. hispida observations in", cv_label))
  }

  # --- 4. CONVERT COVERAGE TO PROPORTION ---
  data$COBERTURA <- as.numeric(gsub(",", ".", as.character(data$COBERTURA)))
  data <- data[!is.na(data$COBERTURA), ]

  # Check range [0,1]
  cover_prop <- data$COBERTURA / 100
  if (any(cover_prop < 0 | cover_prop > 1, na.rm = TRUE)) {
    warning("Coverage outside [0,1]. Correcting...")
    cover_prop <- pmax(0, pmin(1, cover_prop))
  }

  # --- 5. CREATE REEF FROM SITE ---
  data$REEF <- factor(substr(as.character(data$SITE), 1, 3))

  # --- 6. CREATE HABMERGED (PA vs RR_TP) ---
  data$HAB <- as.character(data$HAB)
  data$HABMERGED <- ifelse(data$HAB %in% c("RR", "TP"), "RR_TP", data$HAB)
  data$HABMERGED <- factor(data$HABMERGED, levels = c("PA", "RR_TP"))

  # Check ARCH
  if ("ARCH" %in% names(data)) {
    data$ARCH <- factor(data$ARCH)
  } else {
    warning("ARCH column not found. ARCH interactions will not be included.")
    data$ARCH <- NA
  }

  # --- 7. SCALE PREDICTORS ---
  # NOTE: brms splines (s()) cannot use variable names with underscores or dots
  # We use period (.) as separator which is valid for data but not for brms parameter names
  # So we'll use no separator for spline-compatible names
  pca_predictors <- c("PC1_MAGNITUDE", "PC2_MAGNITUDE",
                      "PC1_VARIABILITY", "PC2_VARIABILITY",
                      "DEPTH_M")

  for (col in pca_predictors) {
    if (col %in% names(data)) {
      # Create scaled version with name compatible with brms splines
      # Replace underscores with nothing for compatibility
      scaled_name <- gsub("_", "", col)
      data[[scaled_name]] <- scale(as.numeric(data[[col]]))[, 1]
    } else {
      warning(paste("Predictor", col, "not found"))
    }
  }

  # --- 8. CREATE RESPONSE VARIABLE WITHOUT UNDERSCORE ---
  data$COVER <- cover_prop

  # --- 9. METADATA ---
  data$CVLABEL <- cv_label

  # --- 10. FINAL VALIDATION ---
  required_cols <- c("COVER", "REEF", "HABMERGED",
                     "PC1MAGNITUDE", "PC2MAGNITUDE",
                     "PC1VARIABILITY", "PC2VARIABILITY")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop(paste("Required columns missing:", paste(missing_cols, collapse = ", ")))
  }

  cat(sprintf("✓ ZOIB data prepared for %s: N=%d observations, %d REEFs\n",
              cv_label, nrow(data), length(unique(data$REEF))))

  return(data)
}


#' Prepare health/growth data for Gaussian models
#' @param data_path Path to CSV file
#' @param cv_label CV scenario label
#' @return data.frame prepared for brms
prepare_gaussian_data <- function(data_path, cv_label) {

  # --- 1. READ WITH DELIMITER DETECTION ---
  data <- tryCatch({
    read.csv2(data_path, stringsAsFactors = TRUE)
  }, error = function(e) {
    read.csv(data_path, stringsAsFactors = TRUE)
  })

  # --- 2. STANDARDIZE COLUMN NAMES ---
  colnames(data) <- toupper(colnames(data))

  # --- 3. VERIFY RESPONSE VARIABLES ---
  response_vars <- c("HEALTH_PC1", "HEALTH_PC2", "RGR")
  for (rv in response_vars) {
    if (!(rv %in% names(data))) {
      stop(paste("Response variable", rv, "not found"))
    }
    data[[rv]] <- as.numeric(gsub(",", ".", as.character(data[[rv]])))
  }

  # --- 4. CREATE REEF FROM SITE ---
  if ("SITE" %in% names(data)) {
    data$REEF <- factor(substr(as.character(data$SITE), 1, 3))
  } else if ("REEF" %in% names(data)) {
    data$REEF <- factor(data$REEF)
  } else {
    stop("Neither SITE nor REEF found in data")
  }

  # --- 5. CREATE HABMERGED (PA vs RR_TP) ---
  # NOTE: Renamed from HAB_MERGED to HABMERGED (no underscores) for brms compatibility
  if ("HAB" %in% names(data)) {
    data$HAB <- as.character(data$HAB)
    data$HABMERGED <- ifelse(data$HAB %in% c("RR", "TP"), "RR_TP", data$HAB)
    data$HABMERGED <- factor(data$HABMERGED, levels = c("PA", "RR_TP"))
  } else {
    stop("HAB column not found")
  }

  # Check ARCH
  if ("ARCH" %in% names(data)) {
    data$ARCH <- factor(data$ARCH)
  }

  # --- 6. SCALE PREDICTORS ---
  # Names may vary: PC1_ENV_MAG or PC1_MAGNITUDE
  # NOTE: Renamed scaled variables to remove underscores for brms spline compatibility
  pca_magnitude <- c("PC1_ENV_MAG", "PC2_ENV_MAG", "PC1_MAGNITUDE", "PC2_MAGNITUDE")
  pca_variability <- c("PC1_ENV_VAR", "PC2_ENV_VAR", "PC1_VARIABILITY", "PC2_VARIABILITY")
  depth_cols <- c("DEPTH", "DEPTH_M", "DEPTH_SCALED")

  # Magnitude - rename to PC1MAGNITUDE, PC2MAGNITUDE (no underscores)
  for (col in pca_magnitude) {
    if (col %in% names(data)) {
      target_name <- gsub("_ENV_MAG", "_MAGNITUDE", col)
      # Remove underscores for brms compatibility
      scaled_name <- gsub("_", "", target_name)
      data[[scaled_name]] <- scale(as.numeric(data[[col]]))[, 1]
    }
  }

  # Variability - rename to PC1VARIABILITY, PC2VARIABILITY (no underscores)
  for (col in pca_variability) {
    if (col %in% names(data)) {
      target_name <- gsub("_ENV_VAR", "_VARIABILITY", col)
      # Remove underscores for brms compatibility
      scaled_name <- gsub("_", "", target_name)
      data[[scaled_name]] <- scale(as.numeric(data[[col]]))[, 1]
    }
  }

  # Depth - rename to DEPTHM (no underscores)
  for (col in depth_cols) {
    if (col %in% names(data)) {
      data$DEPTHM <- scale(as.numeric(data[[col]]))[, 1]
      break
    }
  }

  # NEW: Scale biological interaction variables (PC1_INTERACAO, PC2_INTERACAO)
  # Remove underscore for brms compatibility (s() doesn't accept underscores)
  for (col in c("PC1_INTERACAO", "PC2_INTERACAO")) {
    if (col %in% names(data)) {
      scaled_name <- gsub("_", "", col)  # PC1_INTERACAO -> PC1INTERACAO
      data[[scaled_name]] <- scale(as.numeric(data[[col]]))[, 1]
      cat(sprintf("  Scaled %s -> %s in prepare_gaussian_data\n", col, scaled_name))
    }
  }

  # --- 7. METADATA ---
  data$CVLABEL <- cv_label

  cat(sprintf("✓ Gaussian data prepared for %s: N=%d colonies, %d REEFs\n",
              cv_label, nrow(data), length(unique(data$REEF))))

  return(data)
}


#' Prepare community composition data for JSDM Dirichlet models
#' @param data_path Path to abundance CSV
#' @param cv_label CV scenario label
#' @return data.frame with community composition proportions
prepare_jsdm_data <- function(data_path, cv_label) {

  # --- 1. READ WITH DELIMITER DETECTION ---
  data <- tryCatch({
    read.csv2(data_path, stringsAsFactors = TRUE)
  }, error = function(e) {
    read.csv(data_path, stringsAsFactors = TRUE)
  })

  colnames(data) <- toupper(colnames(data))

  # --- 2. CONVERT COVERAGE ---
  data$COBERTURA <- as.numeric(gsub(",", ".", as.character(data$COBERTURA)))
  data <- data[!is.na(data$COBERTURA), ]
  data$COVER_PROP <- data$COBERTURA / 100

  # --- 3. DEFINE ORGANISM CATEGORIES ---
  data$CATEGORY <- case_when(
    data$ORGANISMO == "MUSSISMILIA_HISPIDA" ~ "MUSSISMILIA",
    grepl("TURF", data$ORGANISMO, ignore.case = TRUE) ~ "TURF",
    grepl("CCA|CORALLINE", data$ORGANISMO, ignore.case = TRUE) ~ "CCA",
    grepl("CYANO", data$ORGANISMO, ignore.case = TRUE) ~ "CYANO",
    grepl("MACRO", data$ORGANISMO, ignore.case = TRUE) ~ "MACROALGAE",
    TRUE ~ "OTHER"
  )

  # --- 4. AGGREGATE BY SAMPLING UNIT ---
  data$SAMPLE_UNIT <- paste(data$SITE, data$HAB, sep = "_")

  # Aggregate coverage by category and sampling unit
  composition <- data %>%
    group_by(SAMPLE_UNIT, SITE, HAB, ARCH,
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

  # --- 5. NORMALIZE TO SUM = 1 ---
  total_cover <- rowSums(composition[, c("MUSSISMILIA_raw", "TURF_raw", "CCA_raw",
                                          "CYANO_raw", "MACROALGAE_raw", "OTHER_raw")])

  # Avoid division by zero
  total_cover[total_cover == 0] <- 1

  composition$MUSSISMILIA_prop <- composition$MUSSISMILIA_raw / total_cover
  composition$TURF_prop <- composition$TURF_raw / total_cover
  composition$CCA_prop <- composition$CCA_raw / total_cover
  composition$CYANO_prop <- composition$CYANO_raw / total_cover
  composition$MACROALGAE_prop <- composition$MACROALGAE_raw / total_cover
  composition$OTHER_prop <- composition$OTHER_raw / total_cover

  # --- 6. ADJUST ZEROS (Dirichlet doesn't allow exactly 0 or 1) ---
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

  # Remove proportion columns that are all zeros (no data for that category)
  # This prevents Dirichlet model from failing due to missing categories
  prop_cols_to_remove <- c()
  for (col in prop_cols) {
    if (all(composition[[col]] <= epsilon * 2)) {  # Essentially all zeros
      cat(sprintf("Removing %s: all values are near zero\n", col))
      composition[[col]] <- NULL
      prop_cols_to_remove <- c(prop_cols_to_remove, col)
    }
  }

  if (length(prop_cols_to_remove) > 0) {
    cat(sprintf("Removed %d zero-value proportion columns\n", length(prop_cols_to_remove)))

    # Update prop_cols to only include remaining columns
    prop_cols <- setdiff(prop_cols, prop_cols_to_remove)

    # Re-normalize remaining proportions to sum to 1
    row_sums <- rowSums(composition[, prop_cols])
    for (col in prop_cols) {
      composition[[col]] <- composition[[col]] / row_sums
    }
    cat("Re-normalized remaining proportions to sum to 1\n")
  }

  # --- 7. CREATE REEF AND HABMERGED ---
  # NOTE: Renamed from HAB_MERGED to HABMERGED (no underscores)
  composition$REEF <- factor(substr(as.character(composition$SITE), 1, 3))
  composition$HABMERGED <- ifelse(composition$HAB %in% c("RR", "TP"), "RR_TP", composition$HAB)
  # Create factor without pre-specifying levels, then drop unused levels
  composition$HABMERGED <- factor(composition$HABMERGED)
  composition$HABMERGED <- droplevels(composition$HABMERGED)

  # Check HABMERGED has sufficient levels - remove if not
  # brms requires factors with 2+ levels for treatment contrasts
  if (length(unique(composition$HABMERGED)) < 2) {
    cat(sprintf("⚠ HABMERGED has only %d level(s) in %s. Removing HABMERGED from data.\n",
                length(unique(composition$HABMERGED)), cv_label))
    composition$HABMERGED <- NULL
  }

  composition$ARCH <- factor(composition$ARCH)
  # Check ARCH has sufficient levels
  if (length(unique(composition$ARCH)) < 2) {
    cat(sprintf("⚠ ARCH has only %d level(s) in %s. Removing ARCH from data.\n",
                length(unique(composition$ARCH)), cv_label))
    composition$ARCH <- NULL
  }

  # --- 8. SCALE PREDICTORS ---
  # NOTE: Renamed to remove underscores for brms spline compatibility
  for (col in c("PC1_MAGNITUDE", "PC2_MAGNITUDE",
                "PC1_VARIABILITY", "PC2_VARIABILITY", "DEPTH_M")) {
    if (col %in% names(composition)) {
      # Remove underscores for brms compatibility
      scaled_name <- gsub("_", "", col)
      composition[[scaled_name]] <- scale(as.numeric(composition[[col]]))[, 1]
    }
  }

  # --- 9. METADATA ---
  composition$CVLABEL <- cv_label

  cat(sprintf("✓ JSDM data prepared for %s: N=%d sampling units, %d REEFs\n",
              cv_label, nrow(composition), length(unique(composition$REEF))))

  return(as.data.frame(composition))
}


# ============================================================================
# SECTION 2: FORMULA BUILDING FUNCTIONS
# ============================================================================

#' Build formula for ZOIB models
#' @param include_hab Include HABMERGED?
#' @param include_depth Include s(DEPTHM)?
#' @param include_arch_interaction Include ARCH interactions?
#' @return brmsformula object for brms
make_zoib_formula <- function(include_hab = TRUE, include_depth = TRUE,
                               include_arch_interaction = TRUE) {
  k_val <- 3

  # NOTE: Variable names must not contain underscores or dots for brms splines
  # Using names: PC1MAGNITUDE, PC2MAGNITUDE, PC1VARIABILITY, PC2VARIABILITY, DEPTHM

  # Build the formula string (main part)
  main_part <- "COVER ~ PC1MAGNITUDE + PC2MAGNITUDE"

  # Add variability splines
  main_part <- paste0(main_part, " + s(PC1VARIABILITY, k = ", k_val, ")")
  main_part <- paste0(main_part, " + s(PC2VARIABILITY, k = ", k_val, ")")

  # Add ARCH interactions if requested
  if (include_arch_interaction) {
    main_part <- paste0(main_part, " + PC1VARIABILITY:ARCH")
    main_part <- paste0(main_part, " + PC2VARIABILITY:ARCH")
  }

  # Add habitat if requested
  if (include_hab) {
    main_part <- paste0(main_part, " + HABMERGED")
  }

  # Add depth spline if requested
  if (include_depth) {
    main_part <- paste0(main_part, " + s(DEPTHM, k = ", k_val, ")")
  }

  # Add random effect
  main_part <- paste0(main_part, " + (1 | REEF)")

  # Convert to formula object using stats::as.formula()
  # This avoids environment variable lookup issues
  main_formula <- stats::as.formula(main_part)

  # Convert to brmsformula
  brms_form <- brmsformula(main_formula)

  # Add ZOIB submodels using lf() (linear formula)
  brms_form <- brms_form + lf(zoi ~ 1) + lf(phi ~ 1)

  return(brms_form)
}


#' Build formula for Gaussian models (Health/RGR)
#' @param response_var Response variable name ("HEALTH_PC1", "HEALTH_PC2", "RGR")
#' @param include_hab Include HABMERGED?
#' @param include_depth Include s(DEPTHM)?
#' @param include_interaction_pca Include interaction PCA terms? (NEW)
#' @return formula object for brms
make_gaussian_formula <- function(response_var, include_hab = TRUE,
                                   include_depth = TRUE,
                                   include_interaction_pca = TRUE) {
  k_val <- 3

  # NOTE: Variable names must not contain underscores for brms splines
  # Using names: PC1MAGNITUDE, PC2MAGNITUDE, PC1VARIABILITY, PC2VARIABILITY, DEPTHM, HABMERGED

  # Build formula string
  main_part <- paste0(response_var, " ~ PC1MAGNITUDE + PC2MAGNITUDE")

  # Add magnitude splines
  main_part <- paste0(main_part, " + s(PC1MAGNITUDE, k = ", k_val, ")")
  main_part <- paste0(main_part, " + s(PC2MAGNITUDE, k = ", k_val, ")")

  # Add variability splines
  main_part <- paste0(main_part, " + s(PC1VARIABILITY, k = ", k_val, ")")
  main_part <- paste0(main_part, " + s(PC2VARIABILITY, k = ", k_val, ")")

  # NEW: Add interaction PCA splines if requested
  if (include_interaction_pca) {
    main_part <- paste0(main_part, " + s(PC1INTERACAO, k = ", k_val, ")")
    main_part <- paste0(main_part, " + s(PC2INTERACAO, k = ", k_val, ")")
  }

  # Add habitat if requested
  if (include_hab) {
    main_part <- paste0(main_part, " + HABMERGED")
  }

  # Add depth spline if requested
  if (include_depth) {
    main_part <- paste0(main_part, " + s(DEPTHM, k = ", k_val, ")")
  }

  # Add random effect
  main_part <- paste0(main_part, " + (1 | REEF)")

  # Convert to formula object using stats::as.formula()
  stats::as.formula(main_part)
}


#' Build formula for JSDM Dirichlet models
#' @param include_hab Include HABMERGED?
#' @param include_depth Include DEPTHM?
#' @param include_arch_interaction Include ARCH interactions?
#' @param data Data frame to check for column availability
#' @return bf() object for brms with Dirichlet
make_jsdm_formula <- function(include_hab = TRUE, include_depth = TRUE,
                               include_arch_interaction = TRUE, data = NULL) {

  # NOTE: Variable names must not contain underscores
  # Using names: PC1MAGNITUDE, PC2MAGNITUDE, PC1VARIABILITY, PC2VARIABILITY, DEPTHM, HABMERGED
  # JSDM uses linear terms (not splines) due to small sample size

  # Check if data is provided and if columns exist
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

  # Magnitude terms (linear)
  mag_terms <- "PC1MAGNITUDE + PC2MAGNITUDE"

  # Variability terms (linear)
  var_terms <- "PC1VARIABILITY + PC2VARIABILITY"

  # ARCH interactions - only if ARCH exists in data
  arch_terms <- if (include_arch_interaction && has_arch) {
    " + PC1VARIABILITY:ARCH + PC2VARIABILITY:ARCH"
  } else ""

  # Optional terms - only if columns exist in data
  hab_term <- if (include_hab && has_habmerged) " + HABMERGED" else ""
  depth_term <- if (include_depth) " + DEPTHM" else ""

  # Formula with multivariate response
  # Note: OTHER_prop might be missing if no OTHER category in data
  # We'll use only the categories that exist
  response_cols <- c("MUSSISMILIA_prop", "TURF_prop", "CCA_prop", "CYANO_prop", "MACROALGAE_prop")
  if (!is.null(data)) {
    response_cols <- intersect(response_cols, names(data))
  }
  response_str <- paste0("cbind(", paste(response_cols, collapse = ", "), ")")

  formula_str <- paste0(
    response_str, " ~ ",
    mag_terms, " + ", var_terms, arch_terms, hab_term, depth_term, " + (1 | REEF)"
  )

  bf(stats::as.formula(formula_str))
}


# ============================================================================
# SECTION 3: MODEL FITTING WITH CACHE
# ============================================================================

#' Load model from cache or fit new model
#' @param model_name Model name (e.g., "FULL")
#' @param cv_name CV name (e.g., "CV_02")
#' @param model_type Model type (e.g., "ZOIB", "HEALTH_PC1")
#' @param formula Model formula
#' @param data Prepared data
#' @param family Model family
#' @param prior Priors
#' @param out_dir Output directory
#' @param brms_args Additional arguments for brm()
#' @return List with fit and cache flag
load_or_fit_model <- function(model_name, cv_name, model_type,
                               formula, data, family, prior,
                               out_dir, brms_args) {

  # File name
  model_filename <- sprintf("%s_%s_%s.rds",
                            tolower(model_type), tolower(model_name), cv_name)
  model_path <- file.path(out_dir, model_filename)

  # Check cache
  if (file.exists(model_path)) {
    cat(sprintf("📦 Cache found: %s\n", model_filename))
    fit <- readRDS(model_path)
    return(list(fit = fit, cached = TRUE, path = model_path))
  }

  # Fit new model
  cat(sprintf("🔧 Fitting model: %s\n", model_filename))
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
  cat(sprintf("✓ Model fitted in %.1f minutes\n", as.numeric(elapsed)))

  # Save immediately
  saveRDS(fit, model_path)
  cat(sprintf("💾 Model saved: %s\n", model_path))

  return(list(fit = fit, cached = FALSE, path = model_path,
              elapsed_mins = as.numeric(elapsed)))
}


# ============================================================================
# SECTION 4: CONVERGENCE CHECKING
# ============================================================================

#' Check convergence criteria for a brms model
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
  rhat_vals <- rhat_vals[!is.na(rhat_vals)]
  rhat_max <- max(rhat_vals)
  rhat_ok <- rhat_max < 1.01

  # --- ESS ---
  ess_bulk_vals <- c(summary_fit$fixed[, "Bulk_ESS"],
                     summary_fit$random$REEF[, "Bulk_ESS"])
  ess_tail_vals <- c(summary_fit$fixed[, "Tail_ESS"],
                     summary_fit$random$REEF[, "Tail_ESS"])

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
  status <- if (all_ok) "✓ CONVERGED" else "⚠ PROBLEMS"
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
# SECTION 5: LOO CALCULATION AND COMPARISON
# ============================================================================

#' Calculate LOO with error handling
#' @param fit brmsfit object
#' @param model_name Model name
#' @param cv_label CV label
#' @param use_moment_match Try moment_match for high Pareto k?
#' @return loo object or NULL on error
safe_loo <- function(fit, model_name, cv_label, use_moment_match = TRUE) {

  cat(sprintf("📊 Calculating LOO for %s_%s...\n", model_name, cv_label))

  loo_result <- tryCatch({
    loo_obj <- loo(fit, cores = 4)

    # Check Pareto k
    k_vals <- loo_obj$diagnostics$pareto_k
    n_high_k <- sum(k_vals > 0.7, na.rm = TRUE)

    if (n_high_k > 0 && use_moment_match) {
      cat(sprintf("⚠ %d observations with Pareto k > 0.7. Trying moment_match...\n", n_high_k))

      loo_obj <- tryCatch({
        loo(fit, cores = 4, moment_match = TRUE)
      }, error = function(e) {
        cat("⚠ moment_match failed. Using original LOO.\n")
        loo_obj
      })
    }

    loo_obj

  }, error = function(e) {
    cat(sprintf("❌ LOO error for %s_%s: %s\n", model_name, cv_label, e$message))
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
    cat(sprintf("⚠ Less than 2 valid models for comparison in %s\n", cv_label))
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
  significance <- case_when(
    delta_second > 2 * se_second ~ "Significant",
    delta_second > se_second ~ "Moderate",
    TRUE ~ "Not significant"
  )

  cat(sprintf("🏆 Winner for %s: %s (Δ=%.1f, SE=%.1f, %s)\n",
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
# SECTION 6: MODEL COMBINATIONS DEFINITION
# ============================================================================

# Definition of model combinations (applicable to all types)
model_combinations <- list(
  # Existing models (without interaction PCA) - KEEP AS IS
  list(name = "FULL",      include_hab = TRUE,  include_depth = TRUE,  include_interaction_pca = FALSE, description = "Full model with HAB, DEPTH (no Interaction PCA)"),
  list(name = "NOHABITAT", include_hab = FALSE, include_depth = TRUE,  include_interaction_pca = FALSE, description = "No HAB, with DEPTH (no Interaction PCA)"),
  list(name = "NODEPTH",   include_hab = TRUE,  include_depth = FALSE, include_interaction_pca = FALSE, description = "With HAB, no DEPTH (no Interaction PCA)"),
  list(name = "MINIMAL",   include_hab = FALSE, include_depth = FALSE, include_interaction_pca = FALSE, description = "Minimal model (PCA only, no Interaction PCA)"),

  # NEW: Models with interaction PCA
  list(name = "FULL_INTERACTION",      include_hab = TRUE,  include_depth = TRUE,  include_interaction_pca = TRUE, description = "Full model with HAB, DEPTH, and Interaction PCA"),
  list(name = "NOHABITAT_INTERACTION", include_hab = FALSE, include_depth = TRUE,  include_interaction_pca = TRUE, description = "No HAB, with DEPTH and Interaction PCA"),
  list(name = "NODEPTH_INTERACTION",   include_hab = TRUE,  include_depth = FALSE, include_interaction_pca = TRUE, description = "With HAB and Interaction PCA, no DEPTH")
)


# ============================================================================
# SECTION 7: PRIORS AND CONFIGURATIONS
# ============================================================================

# Priors for ZOIB (Abundance)
priors_zoib <- c(
  prior(normal(0, 1), class = "sds"),           # Smooths
  prior(normal(0, 0.5), class = "b"),           # Fixed effects
  prior(exponential(2), class = "sd"),          # Random effects SD
  prior(normal(0, 2), class = "Intercept")      # Intercept
)

# Priors for Gaussian (Health/RGR)
priors_gaussian <- c(
  prior(normal(0, 1), class = "sds"),
  prior(normal(0, 0.5), class = "b"),
  prior(exponential(2), class = "sd"),
  prior(normal(0, 1), class = "Intercept")
)

# Priors for JSDM Dirichlet
# NOTE: Dirichlet models have different parameter structure
# Using NULL to use brms defaults (recommended for Dirichlet)
priors_jsdm <- NULL  # Use brms defaults for Dirichlet models

# Standard brms configuration
brms_args_standard <- list(
  backend = "cmdstanr",
  cores = 4,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  control = list(adapt_delta = 0.97, max_treedepth = 12),
  refresh = 500,
  save_pars = save_pars(all = TRUE)
)

# More conservative configuration for JSDM (more complex models)
brms_args_jsdm <- list(
  backend = "cmdstanr",
  cores = 4,
  iter = 6000,       # More iterations
  warmup = 3000,     # More warmup
  chains = 4,
  control = list(adapt_delta = 0.99, max_treedepth = 15),  # More conservative
  refresh = 500,
  save_pars = save_pars(all = TRUE)
)


# ============================================================================
# SECTION 8: PATH DEFINITIONS
# ============================================================================

# Base output directory (hardcoded)
BASE_OUTPUT_DIR <- "C:/Users/rbfra/OneDrive/Bayesian_Full_LOO_Selection_v5/"

# Subdirectories by response type
output_dirs <- list(
  ZOIB       = file.path(BASE_OUTPUT_DIR, "ZOIB_Abundance/"),
  HEALTH_PC1 = file.path(BASE_OUTPUT_DIR, "Gaussian_HEALTH_PC1/"),
  HEALTH_PC2 = file.path(BASE_OUTPUT_DIR, "Gaussian_HEALTH_PC2/"),
  RGR        = file.path(BASE_OUTPUT_DIR, "Gaussian_RGR/"),
  JSDM       = file.path(BASE_OUTPUT_DIR, "JSDM_Dirichlet/")
)

# Project paths
PROJECT_ROOT <- "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/"
RESULTS_DIR <- file.path(PROJECT_ROOT, "#######FINAL_RESULTS/")

# Abundance dataset paths (for ZOIB and JSDM)
dataset_paths_abundance <- list(
  CV_02 = file.path(RESULTS_DIR, "#####output_local_PCA_CV_2_FINAL/dados_abundancia_integrados_long_format.csv"),
  CV_30 = file.path(RESULTS_DIR, "#####output_local_PCA_CV_30_FINAL/dados_abundancia_integrados_long_format.csv"),
  CV_ALL = file.path(RESULTS_DIR, "#####output_local_PCA_CV_all_FINAL/dados_abundancia_integrados_long_format.csv")
)

# Health/RGR dataset paths (for Gaussian models)
dataset_paths_health <- list(
  CV_02 = file.path(RESULTS_DIR, "output_DADOS_FINAIS_PARA_MODELAGEM_cv_2/dados_finais_para_modelagem_com_ARCH.csv"),
  CV_30 = file.path(RESULTS_DIR, "output_DADOS_FINAIS_PARA_MODELAGEM_cv_30/dados_finais_para_modelagem_com_ARCH.csv"),
  CV_ALL = file.path(RESULTS_DIR, "output_DADOS_FINAIS_PARA_MODELAGEM_cv_all/dados_finais_para_modelagem_com_ARCH.csv")
)


# ============================================================================
# END OF 00_LOO_Selection_Functions.R
# ============================================================================
