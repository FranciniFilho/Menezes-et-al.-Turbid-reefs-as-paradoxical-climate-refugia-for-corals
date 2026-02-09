# ============================================================================
# 01_ZOIB_YEAR_RE_Full_LOO_Selection.R
# ============================================================================
# ZOIB Abundance models with YEAR Random Effect
# Full LOO-based model selection across HAB × DEPTH × CV scenarios
# YEAR Random Effect: YES - (1|REEF) + (1|YEAR) - CROSSED (not nested!)
# ============================================================================
# Source: Based on 01_ZOIB_Abundance_Full_LOO_Selection.R
# ============================================================================
# CHANGES FROM ORIGINAL:
# - Now uses load_or_fit_model_year_re() for caching
# - Uses safe_loo() for robust LOO calculation
# - Uses check_convergence() for validation
# - Fixed loo_compare() to use list of loo objects (not brmsfit)
# ============================================================================

# Set working directory to script location
setwd("C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_CODES/Bayes_models/New_Bayes_Models_YEAR_RE/")

# Source extended functions (includes original + new cache/safe_loo functions)
source("00_LOO_Selection_Functions_YEAR_RE.R")

# Load required libraries
library(brms)
library(cmdstanr)
library(loo)
library(dplyr)

# Set options
options(mc.cores = 4)
set.seed(42)

# ============================================================================
# MODEL CONFIGURATIONS
# ============================================================================

# CV scenarios
cv_scenarios <- c("CV_02", "CV_30", "CV_ALL")

# Abundance data paths (contain YEAR column)
dataset_paths_abundance <- list(
  CV_02 = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/#####output_local_PCA_CV_2_FINAL/dados_abundancia_integrados_long_format.csv"),
  CV_30 = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/#####output_local_PCA_CV_30_FINAL/dados_abundancia_integrados_long_format.csv"),
  CV_ALL = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/#####output_local_PCA_CV_all_FINAL/dados_abundancia_integrados_long_format.csv")
)

# Output directory
base_output_dir <- output_dirs_year_re$ZOIB_YEAR_RE
dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# MAIN MODELING LOOP
# ============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   ZOIB ABUNDANCE MODELS WITH YEAR RANDOM EFFECT (CROSSED)     ║\n")
cat("║   Full LOO-based Model Selection                               ║\n")
cat("║   Using: load_or_fit_model, safe_loo, check_convergence       ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

# Store all results for summary
all_results <- data.frame()

for (cv_label in cv_scenarios) {

  cat(sprintf("\n═══════════════════════════════════════════════════════════════\n"))
  cat(sprintf("  PROCESSING: %s\n", cv_label))
  cat(sprintf("═══════════════════════════════════════════════════════════════\n\n"))

  # Create output directory for this CV
  cv_output_dir <- file.path(base_output_dir, cv_label)
  dir.create(cv_output_dir, showWarnings = FALSE, recursive = TRUE)

  # Prepare data with YEAR
  prepared_data <- tryCatch({
    prepare_zoib_data_with_year(
      data_path = dataset_paths_abundance[[cv_label]],
      cv_label = cv_label
    )
  }, error = function(e) {
    cat(sprintf("  ❌ Data preparation failed: %s\n", e$message))
    NULL
  })

  if (is.null(prepared_data)) {
    cat(sprintf("    ⚠ Skipping %s due to data preparation error\n", cv_label))
    next
  }

  # Check YEAR validity
  year_check <- check_year_levels(prepared_data)
  if (!year_check$valid) {
    warning(sprintf("Skipping %s: %s", cv_label, year_check$summary))
    next
  }

  # Store LOO objects for comparison
  loo_list <- list()
  model_info_list <- list()

  for (model_config in model_combinations_year_re) {

    model_name <- sprintf("ZOIB_%s_%s", model_config$name, cv_label)

    cat(sprintf("\n─── Model: %s ───\n", model_config$name))
    cat(sprintf("    %s\n", model_config$description))

    # Build formula with CROSSED YEAR random effect
    formula_obj <- make_zoib_formula_year_re(
      include_hab = model_config$include_hab,
      include_depth = model_config$include_depth,
      include_arch_interaction = model_config$include_arch_interaction
    )

    cat(sprintf("    Formula: %s\n\n", deparse(formula_obj)[1]))

    # Fit or load model from cache
    result <- tryCatch({
      load_or_fit_model_year_re(
        model_name = model_config$name,
        cv_name = cv_label,
        model_type = "ZOIB_YEAR_RE",
        formula = formula_obj,
        data = prepared_data,
        family = zero_one_inflated_beta(),
        prior = priors_zoib_year_re,
        out_dir = cv_output_dir,
        brms_args = brms_args_year_re
      )
    }, error = function(e) {
      cat(sprintf("  ❌ Error fitting model: %s\n", e$message))
      list(fit = NULL, cached = FALSE, error = e$message)
    })

    if (is.null(result$fit)) {
      next
    }

    # Check convergence
    conv_check <- tryCatch({
      check_convergence(result$fit, model_config$name, cv_label)
    }, error = function(e) {
      cat(sprintf("  ⚠ Convergence check error: %s\n", e$message))
      list(passed = FALSE, rhat_max = NA, ess_bulk_min = NA,
           ess_tail_min = NA, n_divergent = NA)
    })

    # Calculate LOO with safe_loo (handles Pareto k issues)
    loo_result <- safe_loo(result$fit, model_config$name, cv_label)

    if (!is.null(loo_result)) {
      loo_list[[model_config$name]] <- loo_result

      # Store results
      result_row <- data.frame(
        Model_Type = "ZOIB_YEAR_RE",
        Response = "COVER",
        Model = model_config$name,
        CV = cv_label,
        LOOIC = loo_result$estimates["looic", "Estimate"],
        SE_LOOIC = loo_result$estimates["looic", "SE"],
        Converged = conv_check$passed,
        Rhat_Max = conv_check$rhat_max,
        ESS_Bulk_Min = conv_check$ess_bulk_min,
        ESS_Tail_Min = conv_check$ess_tail_min,
        N_Divergent = conv_check$n_divergent,
        Cached = result$cached,
        Elapsed_Mins = ifelse(is.null(result$elapsed_mins), NA, result$elapsed_mins),
        stringsAsFactors = FALSE
      )
      all_results <- rbind(all_results, result_row)
    }

    # Store model info
    model_info_list[[model_name]] <- list(
      model_id = model_name,
      cv = cv_label,
      combination = model_config$name,
      fit = result$fit,
      loo = loo_result,
      convergence = conv_check,
      cached = result$cached,
      path = result$path
    )

    # Garbage collection between models
    gc()
  }

  # ======================================================================
  # COMPARE MODELS WITHIN CV
  # ======================================================================

  if (length(loo_list) > 1) {
    cat(sprintf("\n─── LOO Comparison for ZOIB_%s ───\n", cv_label))
    cv_comparison <- tryCatch({
      compare_loo_within_cv(loo_list, cv_label)
    }, error = function(e) {
      cat(sprintf("  ⚠ Comparison error: %s\n", e$message))
      NULL
    })

    # Save winner model copy
    if (!is.null(cv_comparison)) {
      winner_path <- file.path(cv_output_dir,
                               sprintf("WINNER_ZOIB_%s.rds", cv_label))
      source_name <- cv_comparison$winner
      source_path <- file.path(cv_output_dir,
                               sprintf("zoib_year_re_%s_%s.rds",
                                       tolower(source_name), cv_label))
      if (file.exists(source_path) && !file.exists(winner_path)) {
        file.copy(source_path, winner_path)
        cat(sprintf("  🏆 Winner saved: %s\n", basename(winner_path)))
      }
    }
  }
}

# ============================================================================
# GLOBAL SUMMARY
# ============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   ALL ZOIB MODELS COMPLETED                                      ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")

if (nrow(all_results) > 0) {
  cat(sprintf("Total models fitted: %d\n", nrow(all_results)))
  cat(sprintf("Models with successful LOO: %d\n", sum(!is.na(all_results$LOOIC))))
  cat(sprintf("Models converged: %d/%d\n",
              sum(all_results$Converged, na.rm = TRUE),
              sum(!is.na(all_results$Converged))))

  cat("\nResults summary:\n")
  print(all_results[, c("Model", "CV", "LOOIC", "SE_LOOIC", "Converged")])

  # Save combined results
  results_path <- file.path(base_output_dir, "ZOIB_YEAR_RE_combined_results.csv")
  write.csv(all_results, results_path, row.names = FALSE)
  cat(sprintf("\n✓ Combined results saved: %s\n", results_path))
} else {
  cat("⚠ No results to save\n")
}

cat(sprintf("\nOutput directory: %s\n", base_output_dir))
cat("\n✓ Script 01 (ZOIB YEAR RE) completed!\n")
