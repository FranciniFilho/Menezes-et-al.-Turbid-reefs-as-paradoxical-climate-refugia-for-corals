# ============================================================================
# 02a_Gaussian_RGR_Full_LOO_Selection.R
# ============================================================================
# Gaussian models for Relative Growth Rate (RGR)
# Full LOO-based model selection across HAB × DEPTH × CV scenarios
# YEAR Random Effect: NO (RGR is time-integrated metric)
# ============================================================================
# Rationale: RGR = (ln(AR_final) - ln(AR_initial)) / (year_final - year_initial)
# RGR is a time-integrated metric - YEAR variation is already incorporated
# Adding YEAR as random effect would create circularity
# ============================================================================
# Source: Based on 02_Gaussian_Health_RGR_Full_LOO_Selection.R
# ============================================================================
# CHANGES FROM ORIGINAL:
# - Now uses load_or_fit_model_year_re() for caching
# - Uses safe_loo() for robust LOO calculation
# - Uses check_convergence() for validation
# - Fixed loo_compare() to use list of loo objects (not brmsfit)
# ============================================================================

# Set working directory
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

# Response variables (RGR ONLY - NO HEALTH)
response_variables <- c("RGR")

# CV scenarios
cv_scenarios <- c("CV_02", "CV_30", "CV_ALL")

# Gaussian data paths (same as original - NO YEAR in RGR data)
dataset_paths_gaussian <- list(
  CV_02 = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/output_DADOS_FINAIS_PARA_MODELAGEM_cv_2/dados_finais_para_modelagem_com_ARCH.csv"),
  CV_30 = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/output_DADOS_FINAIS_PARA_MODELAGEM_cv_30/dados_finais_para_modelagem_com_ARCH.csv"),
  CV_ALL = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/output_DADOS_FINAIS_PARA_MODELAGEM_cv_all/dados_finais_para_modelagem_com_ARCH.csv")
)

# Output directory
base_output_dir <- output_dirs_year_re$RGR
dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# MAIN MODELING LOOP
# ============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   GAUSSIAN RGR MODELS (NO YEAR RANDOM EFFECT)                  ║\n")
cat("║   Full LOO-based Model Selection                               ║\n")
cat("║   Using: load_or_fit_model, safe_loo, check_convergence       ║\n")
cat("╚════════════════════════════════════════════════════════════════╝\n")
cat("\n")
cat("NOTE: RGR is a time-integrated metric. YEAR variation is already\n")
cat("      incorporated in its calculation. Using (1|REEF) only.\n\n")

# Store all results for summary
all_results <- data.frame()

for (cv_label in cv_scenarios) {
  for (response_var in response_variables) {

    cat(sprintf("\n═══════════════════════════════════════════════════════════════\n"))
    cat(sprintf("  PROCESSING: %s - %s\n", cv_label, response_var))
    cat(sprintf("═══════════════════════════════════════════════════════════════\n\n"))

    # Create output directory for this response
    response_output_dir <- file.path(base_output_dir, response_var, cv_label)
    dir.create(response_output_dir, showWarnings = FALSE, recursive = TRUE)

    # Prepare data (using original function - NO YEAR handling)
    prepared_data <- tryCatch({
      prepare_gaussian_data(
        data_path = dataset_paths_gaussian[[cv_label]],
        cv_label = cv_label
      )
    }, error = function(e) {
      cat(sprintf("  ❌ Data preparation failed: %s\n", e$message))
      NULL
    })

    if (is.null(prepared_data)) {
      cat(sprintf("    ⚠ Skipping %s - %s due to data preparation error\n", cv_label, response_var))
      next
    }

    # Verify response variable exists
    if (!(response_var %in% names(prepared_data))) {
      cat(sprintf("    ⚠ Response variable %s not found in data\n", response_var))
      next
    }

    # Remove NA values for this response
    prepared_data <- prepared_data[!is.na(prepared_data[[response_var]]), ]
    cat(sprintf("    Data after NAs removed: N=%d\n", nrow(prepared_data)))

    # Data summary
    cat(sprintf("    Response range: [%.2f, %.2f]\n",
                min(prepared_data[[response_var]], na.rm = TRUE),
                max(prepared_data[[response_var]], na.rm = TRUE)))
    cat(sprintf("    REEFs: %d, HABMERGED levels: %d\n",
                length(unique(prepared_data$REEF)),
                length(unique(prepared_data$HABMERGED))))

    # Store LOO objects for comparison
    loo_list <- list()
    model_info_list <- list()

    for (model_config in model_combinations_year_re) {

      model_name <- sprintf("%s_%s_%s", response_var, model_config$name, cv_label)

      cat(sprintf("\n─── Model: %s ───\n", model_config$name))
      cat(sprintf("    %s\n", model_config$description))

      # Build formula (using original function - NO YEAR)
      formula_obj <- make_gaussian_formula(
        response_var = response_var,
        include_hab = model_config$include_hab,
        include_depth = model_config$include_depth,
        include_interaction_pca = model_config$include_interaction_pca  # NEW
      )

      cat(sprintf("    Formula: %s\n\n", deparse(formula_obj)[1]))

      # Fit or load model from cache
      result <- tryCatch({
        load_or_fit_model_year_re(
          model_name = model_config$name,
          cv_name = cv_label,
          model_type = sprintf("Gaussian_%s", response_var),
          formula = formula_obj,
          data = prepared_data,
          family = gaussian(),
          prior = priors_gaussian_year_re,
          out_dir = response_output_dir,
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
          Model_Type = sprintf("Gaussian_%s", response_var),
          Response = response_var,
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
        response = response_var,
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
      cat(sprintf("\n─── LOO Comparison for %s_%s ───\n", response_var, cv_label))
      cv_comparison <- tryCatch({
        compare_loo_within_cv(loo_list, cv_label)
      }, error = function(e) {
        cat(sprintf("  ⚠ Comparison error: %s\n", e$message))
        NULL
      })

      # Save winner model copy
      if (!is.null(cv_comparison)) {
        winner_path <- file.path(response_output_dir,
                                 sprintf("WINNER_%s_%s.rds",
                                         tolower(response_var),
                                         tolower(cv_label)))
        source_name <- cv_comparison$winner
        source_path <- file.path(response_output_dir,
                                 sprintf("gaussian_%s_%s_%s.rds",
                                         tolower(response_var),
                                         tolower(source_name), cv_label))
        if (file.exists(source_path) && !file.exists(winner_path)) {
          file.copy(source_path, winner_path)
          cat(sprintf("  🏆 Winner saved: %s\n", basename(winner_path)))
        }
      }
    }
  }
}

# ============================================================================
# GLOBAL SUMMARY
# ============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════╗\n")
cat("║   RGR MODELS COMPLETED                                          ║\n")
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
  results_path <- file.path(base_output_dir, "RGR_combined_results.csv")
  write.csv(all_results, results_path, row.names = FALSE)
  cat(sprintf("\n✓ Combined results saved: %s\n", results_path))
} else {
  cat("⚠ No results to save\n")
}

cat(sprintf("\nOutput directory: %s\n", base_output_dir))
cat("\n✓ Script 02a (Gaussian RGR) completed!\n")
