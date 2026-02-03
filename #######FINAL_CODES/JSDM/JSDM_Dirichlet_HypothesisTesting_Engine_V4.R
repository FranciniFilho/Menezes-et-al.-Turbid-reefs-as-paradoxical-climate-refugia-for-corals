# ============================================================================ #
# SCRIPT: JSDM_Dirichlet_HypothesisTesting_Engine_V4.R
# PURPOSE: Fit 10 candidate models across 3 CV datasets with CORRECT LOO comparison
# VERSION: 4.0 (All V3 issues resolved)
# ============================================================================ #

rm(list = ls())
gc()

# ============================================================================ #
# PACKAGES AND CONFIGURATION
# ============================================================================ #

libs <- c(
  "brms", "dplyr", "loo", "cmdstanr", "bayesplot", "tidyr",
  "readxl", "janitor", "digest", "ggplot2"
)
invisible(lapply(libs, library, character.only = TRUE))

options(mc.cores = 4)
set.seed(42)

cat("\n", rep("=", 70), "\n", sep = "")
cat("=== JSDM V4: CORRECTED MULTI-CV ANALYSIS ===\n")
cat(rep("=", 70), "\n\n", sep = "")

# ============================================================================ #
# CONFIGURATION
# ============================================================================ #

# Convergence thresholds
CONVERGENCE_THRESHOLDS <- list(
  rhat_max = 1.01,
  ess_bulk_min = 400,
  ess_tail_min = 400,
  divergent_max_pct = 1
)

# Model fitting parameters - OPTIMIZED FOR STABILITY AND EFFICIENCY
# Based on successful ZOIB configuration and user feedback on hardware
MODEL_CONFIG <- list(
  chains = 4,
  cores = 4,
  threads = 2, # Intra-chain parallelization
  iter = 4000, # Lower iterations for faster initial feedback
  warmup = 1500,
  adapt_delta = 0.95, # High quality but faster than 0.99
  max_treedepth = 12, # Balanced depth
  use_opencl = FALSE # Disabled to avoid initialization hangs
)

# Data preparation parameters
DATA_CONFIG <- list(
  n_top_species = 10,
  focal_species = "mussismilia_hispida",
  min_prevalence = 0.10,
  zero_method = "half_min"
)

# ============================================================================ #
# PATHS
# ============================================================================ #

# --- Paths ---
main_out_dir <- "C:/Users/rbfra/OneDrive/JSDM_Dirichlet_HypothesisTesting_V4"
dir.create(main_out_dir, showWarnings = FALSE, recursive = TRUE)
brms_cache_dir <- file.path(main_out_dir, "brms_cache")
dir.create(brms_cache_dir, showWarnings = FALSE, recursive = TRUE)

# Source utility functions
source("C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_CODES/JSDM/JSDM_utils_V4.R")

# ============================================================================ #
# FUNCTION: prepare_data_for_modeling
# PURPOSE: Transform integrated data into format ready for brms
# ============================================================================ #

prepare_data_for_modeling <- function(integrated_data, data_config) {
  #' @param integrated_data List from JSDM_Data_Integration_V4.R
  #' @param data_config List of data preparation parameters
  #' @return List with prepared data and metadata

  cat("\n  Preparing data for modeling...\n")

  data_raw <- integrated_data$data
  organism_cols <- integrated_data$organism_cols

  # Select organism columns only
  organism_data <- data_raw[, organism_cols, drop = FALSE]

  # --- Step 1: Species selection with "Others" ---
  cat("    [1/4] Selecting top species with 'Others' category...\n")
  selection_result <- select_top_species_with_others(
    community_matrix = as.matrix(organism_data),
    n_top = data_config$n_top_species,
    focal_species = data_config$focal_species,
    min_prevalence = data_config$min_prevalence
  )

  final_organism_vars <- selection_result$species
  composition_matrix <- selection_result$matrix

  # --- Step 2: Zero substitution ---
  cat("    [2/4] Handling zeros for Dirichlet model...\n")
  dirichlet_ready <- prepare_dirichlet_data(
    composition_matrix = composition_matrix,
    method = data_config$zero_method
  )

  # --- Step 3: Prepare predictors ---
  cat("    [3/4] Scaling predictors...\n")

  # Get predictor columns
  predictor_cols <- c(
    "SITE", "HAB", "PC1_MAGNITUDE", "PC2_MAGNITUDE",
    "PC1_VARIABILITY", "PC2_VARIABILITY"
  )

  # Extract and scale predictors
  predictor_data <- data_raw[, predictor_cols, drop = FALSE]

  # Scale continuous predictors
  continuous_vars <- c(
    "PC1_MAGNITUDE", "PC2_MAGNITUDE",
    "PC1_VARIABILITY", "PC2_VARIABILITY"
  )
  predictor_data[, continuous_vars] <- scale(predictor_data[, continuous_vars])

  # Convert factors
  predictor_data$SITE <- as.factor(predictor_data$SITE)
  predictor_data$HAB <- as.factor(predictor_data$HAB)

  # --- Step 4: Combine and verify ---
  cat("    [4/4] Combining and verifying data...\n")

  data_for_brms <- cbind(predictor_data, as.data.frame(dirichlet_ready))

  # Verify integrity
  integrity_check <- verify_data_integrity(data_for_brms, final_organism_vars)

  return(list(
    data = data_for_brms,
    organism_vars = final_organism_vars,
    metadata = list(
      cv_name = integrated_data$cv_name,
      n_obs = nrow(data_for_brms),
      n_sites = n_distinct(data_for_brms$SITE),
      n_habs = n_distinct(data_for_brms$HAB),
      n_species = length(final_organism_vars),
      pct_in_others = selection_result$pct_in_others,
      integrity_passed = integrity_check$passed
    )
  ))
}

# ============================================================================ #
# FUNCTION: define_candidate_models
# PURPOSE: Define 10 candidate model formulas
# ============================================================================ #

define_candidate_models <- function(organism_vars) {
  #' @param organism_vars Character vector of organism (species) column names
  #' @return Named list of brms formulas

  # Build response part: cbind(species1, species2, ..., others)
  response_part <- sprintf("cbind(%s)", paste(organism_vars, collapse = ", "))

  models <- list(
    # Model 1: Null model (intercept only)
    model_null = bf(sprintf("%s ~ 1 + (1 | SITE / HAB)", response_part),
      family = dirichlet()
    ),

    # Model 2: Magnitude only (linear)
    model_magnitude_linear = bf(
      sprintf("%s ~ PC1_MAGNITUDE + PC2_MAGNITUDE + (1 | SITE / HAB)", response_part),
      family = dirichlet()
    ),

    # Model 3: Variability only (linear)
    model_variability_linear = bf(
      sprintf("%s ~ PC1_VARIABILITY + PC2_VARIABILITY + (1 | SITE / HAB)", response_part),
      family = dirichlet()
    ),

    # Model 4: Full linear (all predictors)
    model_full_linear = bf(
      sprintf(
        "%s ~ PC1_MAGNITUDE + PC2_MAGNITUDE + PC1_VARIABILITY + PC2_VARIABILITY + (1 | SITE / HAB)",
        response_part
      ),
      family = dirichlet()
    ),

    # Model 5: Magnitude splines
    model_magnitude_spline = bf(
      sprintf(
        "%s ~ s(PC1_MAGNITUDE, k = 3) + s(PC2_MAGNITUDE, k = 3) + (1 | SITE / HAB)",
        response_part
      ),
      family = dirichlet()
    ),

    # Model 6: Variability splines
    model_variability_spline = bf(
      sprintf(
        "%s ~ s(PC1_VARIABILITY, k = 3) + s(PC2_VARIABILITY, k = 3) + (1 | SITE / HAB)",
        response_part
      ),
      family = dirichlet()
    ),

    # Model 7: Full splines
    model_full_spline = bf(
      sprintf(
        "%s ~ s(PC1_MAGNITUDE, k = 3) + s(PC2_MAGNITUDE, k = 3) + s(PC1_VARIABILITY, k = 3) + s(PC2_VARIABILITY, k = 3) + (1 | SITE / HAB)",
        response_part
      ),
      family = dirichlet()
    ),

    # Model 8: Magnitude linear + habitat
    model_magnitude_habitat = bf(
      sprintf("%s ~ PC1_MAGNITUDE + PC2_MAGNITUDE + HAB + (1 | SITE / HAB)", response_part),
      family = dirichlet()
    ),

    # Model 9: Variability linear + habitat
    model_variability_habitat = bf(
      sprintf("%s ~ PC1_VARIABILITY + PC2_VARIABILITY + HAB + (1 | SITE / HAB)", response_part),
      family = dirichlet()
    ),

    # Model 10: Full linear + habitat
    model_full_habitat = bf(
      sprintf(
        "%s ~ PC1_MAGNITUDE + PC2_MAGNITUDE + PC1_VARIABILITY + PC2_VARIABILITY + HAB + (1 | SITE / HAB)",
        response_part
      ),
      family = dirichlet()
    )
  )

  return(models)
}

# ============================================================================ #
# FUNCTION: define_priors
# PURPOSE: Define priors for models
# ============================================================================ #

define_priors <- function() {
  #' @return List of priors for brms

  priors_weakly <- c(
    prior(student_t(3, 0, 2.5), class = "Intercept"),
    prior(normal(0, 1.5), class = "b"),
    prior(normal(0, 1), class = "sds"),
    prior(exponential(1), class = "sd")
  )

  return(priors_weakly)
}

# ============================================================================ #
# FUNCTION: fit_single_model
# PURPOSE: Fit a single model with caching
# ============================================================================ #

fit_single_model <- function(model_name, model_formula, data_for_brms,
                             priors, config, cache_dir, dataset_name) {
  #' @param model_name Name of the model
  #' @param model_formula brms formula object
  #' @param data_for_brms Prepared data frame
  #' @param priors Prior specification (unused if we switch to defaults)
  #' @param config Model fitting configuration
  #' @param cache_dir Directory for cached models
  #' @param dataset_name Name of CV dataset
  #' @return Fitted brms object or NULL on failure

  # NOTE: We are using DEFAULT priors because Dirichlet models have complex
  # parameter structures that conflict with global prior definitions.
  priors_to_use <- NULL

  cache_file <- file.path(cache_dir, sprintf("JSDM_%s_%s", dataset_name, model_name))

  # Generate hash for deterministic caching
  data_hash <- generate_cache_hash(data_for_brms, model_formula, priors_to_use)
  hash_file <- paste0(cache_file, "_hash.rds")

  needs_refit <- TRUE
  if (file.exists(hash_file) && file.exists(paste0(cache_file, ".rds"))) {
    stored_hash <- readRDS(hash_file)
    if (stored_hash == data_hash) {
      cat("     -> Using cached model\n")
      fit <- readRDS(paste0(cache_file, ".rds"))
      needs_refit <- FALSE
    }
  }

  if (needs_refit) {
    cat(sprintf("     -> Fitting model %s (defaults priors)...\n", model_name))

    fit <- tryCatch(
      {
        brm(
          formula = model_formula,
          data = data_for_brms,
          # Using NULL to let brms use its own (usually student-t) defaults
          prior = NULL,
          family = dirichlet(),
          chains = config$chains,
          cores = config$cores,
          threads = threading(config$threads), # Intra-chain parallelization
          opencl = if (config$use_opencl) opencl(ids = c(0, 0)) else NULL, # GPU acceleration
          iter = config$iter,
          warmup = config$warmup,
          control = list(
            adapt_delta = config$adapt_delta,
            max_treedepth = config$max_treedepth
          ),
          backend = "cmdstanr",
          save_pars = save_pars(all = FALSE),
          init = 0,
          silent = 2,
          refresh = 500
        )
      },
      error = function(e) {
        cat(sprintf("     ERROR fitting %s: %s\n", model_name, e$message))
        return(NULL)
      }
    )

    if (!is.null(fit)) {
      # Save with hash
      saveRDS(fit, paste0(cache_file, ".rds"))
      saveRDS(data_hash, hash_file)
    }
  }

  return(fit)
}

# ============================================================================ #
# MAIN EXECUTION
# ============================================================================ #

# --- Global Storage ---
all_loo_by_dataset <- list()
all_winners_by_dataset <- list()
convergence_issues <- list()
model_metadata <- list()

# --- Load integrated datasets ---
cat("\n### LOADING INTEGRATED DATASETS ###\n")

dataset_names <- c("CV_02", "CV_30", "CV_ALL")

for (dataset_name in dataset_names) {
  input_file <- file.path(main_out_dir, sprintf("JSDM_integrated_%s.rds", dataset_name))

  if (!file.exists(input_file)) {
    stop(sprintf("Integrated dataset not found: %s\nRun JSDM_Data_Integration_V4.R first.", input_file))
  }

  cat(sprintf("Loading %s...\n", dataset_name))
  integrated_data <- readRDS(input_file)
  assign(sprintf("integrated_%s", dataset_name), integrated_data)
}

# ============================================================================ #
# MAIN LOOP: PROCESS EACH CV DATASET INDEPENDENTLY
# ============================================================================ #

for (dataset_name in dataset_names) {
  cat(sprintf("\n### PROCESSING DATASET: %s ###\n", dataset_name))
  cat(rep("-", 70), "\n", sep = "")

  # Load integrated data
  integrated_data <- get(sprintf("integrated_%s", dataset_name))

  # --- Data Preparation ---
  cat("\n[Step 1] Data Preparation\n")
  prep_result <- prepare_data_for_modeling(integrated_data, DATA_CONFIG)
  data_for_brms <- prep_result$data
  organism_vars <- prep_result$organism_vars

  if (!prep_result$metadata$integrity_passed) {
    warning(sprintf("Data integrity check failed for %s", dataset_name))
  }

  cat(sprintf(
    "  Final data dimensions: %d obs x %d variables\n",
    nrow(data_for_brms), ncol(data_for_brms)
  ))
  cat(sprintf("  Species included: %s\n", paste(organism_vars, collapse = ", ")))

  # --- Define Candidate Models ---
  cat("\n[Step 2] Defining Candidate Models\n")
  candidate_models <- define_candidate_models(organism_vars)
  cat(sprintf("  Defined %d candidate models\n", length(candidate_models)))

  # --- Define Priors ---
  priors <- define_priors()

  # --- FIT ALL MODELS FOR THIS DATASET ---
  cat("\n[Step 3] Fitting Models\n")
  cat(sprintf(
    "  Configuration: %d chains, %d iterations, %d warmup\n",
    MODEL_CONFIG$chains, MODEL_CONFIG$iter, MODEL_CONFIG$warmup
  ))
  cat("  Note: This will take MANY hours. Be patient.\n\n")

  fitted_models <- list()
  loo_results <- list()
  dataset_metadata <- list()

  for (model_name in names(candidate_models)) {
    cat(sprintf("  >> Fitting: %s (%s)\n", model_name, dataset_name))

    model_formula <- candidate_models[[model_name]]

    # Fit model
    fit <- fit_single_model(
      model_name = model_name,
      model_formula = model_formula,
      data_for_brms = data_for_brms,
      priors = priors,
      config = MODEL_CONFIG,
      cache_dir = brms_cache_dir,
      dataset_name = dataset_name
    )

    if (is.null(fit)) {
      cat(sprintf("  !! SKIPPED: %s failed to fit\n", model_name))
      next
    }

    # --- Check Convergence ---
    conv_check <- check_convergence(fit, CONVERGENCE_THRESHOLDS)
    dataset_metadata[[model_name]] <- list(
      converged = conv_check$passed,
      convergence_metrics = conv_check$metrics
    )

    if (!conv_check$passed) {
      convergence_issues[[paste(dataset_name, model_name, sep = "_")]] <- conv_check$metrics
      cat(sprintf("  !! WARNING: Convergence issues detected\n"))
    }

    # --- Compute LOO ---
    cat("     -> Computing LOO...\n")
    loo_result <- tryCatch(
      {
        loo(fit, moment_match = FALSE, cores = MODEL_CONFIG$cores)
      },
      error = function(e) {
        cat(sprintf("     ERROR computing LOO: %s\n", e$message))
        NULL
      }
    )

    if (!is.null(loo_result)) {
      loo_results[[model_name]] <- loo_result
    }

    fitted_models[[model_name]] <- list(
      fit = fit,
      loo = loo_result,
      converged = conv_check$passed
    )

    # --- Memory Cleanup ---
    rm(fit)
    cleanup_after_model(model_name)
  }

  # --- WITHIN-DATASET LOO COMPARISON (CORRECT APPROACH) ---
  cat(sprintf("\n[Step 4] LOO Comparison for %s\n", dataset_name))

  if (length(loo_results) > 1) {
    comp <- loo_compare(loo_results)
    all_loo_by_dataset[[dataset_name]] <- list(
      results = loo_results,
      comparison = comp
    )

    winner <- rownames(comp)[1]
    all_winners_by_dataset[[dataset_name]] <- winner

    cat(sprintf("\n  LOO Comparison Table:\n"))
    print(comp)

    cat(sprintf("\n  Winner: %s\n", winner))

    # Save per-dataset comparison
    capture.output(print(comp),
      file = file.path(
        main_out_dir,
        sprintf("LOO_Comparison_%s.txt", dataset_name)
      )
    )
  } else {
    cat(sprintf(
      "  WARNING: Only %d model(s) successfully fitted for LOO comparison\n",
      length(loo_results)
    ))
    if (length(loo_results) == 1) {
      all_winners_by_dataset[[dataset_name]] <- names(loo_results)[1]
    }
  }

  # Store metadata
  model_metadata[[dataset_name]] <- dataset_metadata

  cat(sprintf("\n=== Completed %s ===\n", dataset_name))
}

# ============================================================================ #
# CROSS-DATASET CONSISTENCY ANALYSIS (CORRECT APPROACH)
# ============================================================================ #

cat("\n", rep("=", 70), "\n", sep = "")
cat("=== CROSS-DATASET CONSISTENCY ANALYSIS ===\n")
cat(rep("=", 70), "\n\n", sep = "")

cat("Winners by dataset:\n")
for (ds in names(all_winners_by_dataset)) {
  cat(sprintf("  - %s: %s\n", ds, all_winners_by_dataset[[ds]]))
}

# Check consistency
unique_winners <- unique(unlist(all_winners_by_dataset))

if (length(unique_winners) == 1) {
  cat(sprintf("\n CONSISTENT WINNER: %s wins across all CV scenarios\n", unique_winners[1]))
  cat("   -> Model selection is robust to temporal scale choice\n")
  final_winner <- unique_winners[1]
} else {
  cat(sprintf("\n VARIED RESULTS: %d different winners\n", length(unique_winners)))
  cat("   -> Model performance depends on temporal scale\n")
  cat("   -> Report this as a finding (temporal scale sensitivity)\n")

  # Use majority vote
  winner_freq <- table(unlist(all_winners_by_dataset))
  final_winner <- names(sort(winner_freq, decreasing = TRUE))[1]
  cat(sprintf("   -> Majority winner: %s\n", final_winner))

  print(winner_freq)
}

# ============================================================================ #
# SAVE FINAL OUTPUTS
# ============================================================================ #

cat("\n", rep("=", 70), "\n", sep = "")
cat("=== SAVING FINAL OUTPUTS ===\n")
cat(rep("=", 70), "\n\n", sep = "")

winner_dir <- file.path(main_out_dir, "FINAL_WINNER_OUTPUTS")
dir.create(winner_dir, showWarnings = FALSE)

# Save summary
summary_df <- data.frame(
  dataset = names(all_winners_by_dataset),
  winner = unlist(all_winners_by_dataset)
)
write.csv(summary_df, file.path(winner_dir, "winners_by_dataset.csv"), row.names = FALSE)

# Save model metadata
saveRDS(model_metadata, file.path(winner_dir, "model_metadata.rds"))

# Save all LOO results
saveRDS(all_loo_by_dataset, file.path(winner_dir, "all_loo_results.rds"))

# Report convergence issues
if (length(convergence_issues) > 0) {
  cat("\n Models with convergence issues:\n")
  for (model_id in names(convergence_issues)) {
    cat(sprintf("  - %s\n", model_id))
  }
  saveRDS(convergence_issues, file.path(winner_dir, "convergence_issues.rds"))
} else {
  cat("\n All models converged successfully!\n")
}

# ============================================================================ #
# FINAL REPORT
# ============================================================================ #

cat("\n", rep("=", 70), "\n", sep = "")
cat("=== ANALYSIS COMPLETE ===\n")
cat(rep("=", 70), "\n\n", sep = "")

cat(sprintf("Results directory: %s\n", main_out_dir))
cat(sprintf("Winner models: %s\n", file.path(winner_dir, "winners_by_dataset.csv")))
cat(sprintf("All LOO results: %s\n", file.path(winner_dir, "all_loo_results.rds")))

cat("\nNext steps:\n")
cat("  1. Review LOO comparison files in main output directory\n")
cat("  2. Check convergence_issues.rds if any models failed\n")
cat("  3. Run JSDM_Generate_Winner_Viz_V4.R for publication figures\n")
cat(sprintf("  4. Final winner: %s\n", final_winner))

cat("\n")
