# ============================================================================
# 02b_Gaussian_Health_YEAR_RE_Full_LOO_Selection.R
# ============================================================================
# Gaussian models for Health (PC1, PC2) with YEAR Random Effect
# Full LOO-based model selection across HAB ?? DEPTH ?? CV scenarios
# YEAR Random Effect: YES - (1|REEF) + (1|YEAR) - CROSSED (not nested!)
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

# Response variables (Health ONLY - NO RGR)
response_variables <- c("HEALTH_PC1", "HEALTH_PC2")

# CV scenarios
cv_scenarios <- c("CV_02", "CV_30", "CV_ALL")

# Abundance data paths (contain YEAR column - needed for merging)
dataset_paths_abundance <- list(
  CV_02 = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/#####output_local_PCA_CV_2_FINAL/dados_abundancia_integrados_long_format.csv"),
  CV_30 = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/#####output_local_PCA_CV_30_FINAL/dados_abundancia_integrados_long_format.csv"),
  CV_ALL = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/#####output_local_PCA_CV_all_FINAL/dados_abundancia_integrados_long_format.csv")
)

# Final data paths (contain PC1_INTERACAO, PC2_INTERACAO - NEW)
dataset_paths_final <- list(
  CV_02 = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/output_DADOS_FINAIS_PARA_MODELAGEM_cv_2/dados_finais_para_modelagem_com_ARCH.csv"),
  CV_30 = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/output_DADOS_FINAIS_PARA_MODELAGEM_cv_30/dados_finais_para_modelagem_com_ARCH.csv"),
  CV_ALL = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/output_DADOS_FINAIS_PARA_MODELAGEM_cv_all/dados_finais_para_modelagem_com_ARCH.csv")
)

# Output directory
base_output_dir <- output_dirs_year_re$HEALTH_YEAR_RE
dir.create(base_output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================================
# MAIN MODELING LOOP
# ============================================================================

cat("\n")
cat("??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????\n")
cat("???   GAUSSIAN HEALTH MODELS WITH YEAR RANDOM EFFECT (CROSSED)     ???\n")
cat("???   Full LOO-based Model Selection                               ???\n")
cat("???   Using: load_or_fit_model, safe_loo, check_convergence       ???\n")
cat("??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????\n")
cat("\n")

# Store all results for summary
all_results <- data.frame()
cv_window_loo_registry <- list()
cv_window_convergence_registry <- list()

for (cv_label in cv_scenarios) {
  for (response_var in response_variables) {

    cat(sprintf("\n?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????\n"))
    cat(sprintf("  PROCESSING: %s - %s\n", cv_label, response_var))
    cat(sprintf("?????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????\n\n"))

    # Create output directory for this response
    response_output_dir <- file.path(base_output_dir, response_var, cv_label)
    dir.create(response_output_dir, showWarnings = FALSE, recursive = TRUE)

    # Prepare data with YEAR (from raw Excel)
    prepared_data <- tryCatch({
      prepare_health_data_with_year(
        excel_path = RAW_HEALTH_EXCEL,
        pca_scores_path = PCA_SCORES_PATH,
        env_data_path = dataset_paths_abundance[[cv_label]],
        final_data_path = dataset_paths_final[[cv_label]],  # NEW: for interaction variables
        cv_label = cv_label
      )
    }, error = function(e) {
      cat(sprintf("  ??? Data preparation failed: %s\n", e$message))
      NULL
    })

    if (is.null(prepared_data)) {
      cat(sprintf("    ??? Skipping %s - %s due to data preparation error\n", cv_label, response_var))
      next
    }

    # Check YEAR validity
    year_check <- check_year_levels(prepared_data)
    if (!year_check$valid) {
      warning(sprintf("Skipping %s - %s: %s", cv_label, response_var, year_check$summary))
      next
    }

    # Verify response variable exists
    if (!(response_var %in% names(prepared_data))) {
      cat(sprintf("    ??? Response variable %s not found in data\n", response_var))
      next
    }

    # Remove NA values for this response
    prepared_data <- prepared_data[!is.na(prepared_data[[response_var]]), ]
    cat(sprintf("    Data after NAs removed: N=%d\n", nrow(prepared_data)))

    # Data summary
    cat(sprintf("    Response range: [%.2f, %.2f]\n",
                min(prepared_data[[response_var]], na.rm = TRUE),
                max(prepared_data[[response_var]], na.rm = TRUE)))
    cat(sprintf("    REEFs: %d, HABMERGED levels: %d, YEARS: %d\n",
                length(unique(prepared_data$REEF)),
                length(unique(prepared_data$HABMERGED)),
                length(unique(prepared_data$YEAR))))

    # Store LOO objects for comparison
    loo_list <- list()
    convergence_map <- list()
    model_info_list <- list()

    for (model_config in model_combinations_year_re) {

      model_name <- sprintf("%s_%s_%s", response_var, model_config$name, cv_label)

      cat(sprintf("\n????????? Model: %s ?????????\n", model_config$name))
      cat(sprintf("    %s\n", model_config$description))

      # Build formula with CROSSED YEAR random effect
      formula_obj <- make_gaussian_formula_year_re(
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
          model_type = sprintf("Gaussian_%s_YEAR_RE", response_var),
          formula = formula_obj,
          data = prepared_data,
          family = gaussian(),
          prior = priors_gaussian_year_re,
          out_dir = response_output_dir,
          brms_args = brms_args_year_re
        )
      }, error = function(e) {
        cat(sprintf("  ??? Error fitting model: %s\n", e$message))
        list(fit = NULL, cached = FALSE, error = e$message)
      })

      if (is.null(result$fit)) {
        next
      }

      # Check convergence
      conv_check <- tryCatch({
        check_convergence(result$fit, model_config$name, cv_label)
      }, error = function(e) {
        cat(sprintf("  ??? Convergence check error: %s\n", e$message))
        list(passed = FALSE, rhat_max = NA, ess_bulk_min = NA,
             ess_tail_min = NA, n_divergent = NA)
      })

      # Calculate LOO with safe_loo (handles Pareto k issues)
      loo_result <- safe_loo(result$fit, model_config$name, cv_label)
      convergence_map[[model_config$name]] <- isTRUE(conv_check$passed)

      if (!is.null(loo_result)) {
        loo_list[[model_config$name]] <- loo_result
        key_name <- sprintf("%s__%s", response_var, model_config$name)
        if (is.null(cv_window_loo_registry[[key_name]])) {
          cv_window_loo_registry[[key_name]] <- list()
          cv_window_convergence_registry[[key_name]] <- list()
        }
        cv_window_loo_registry[[key_name]][[cv_label]] <- loo_result
        cv_window_convergence_registry[[key_name]][[cv_label]] <- isTRUE(conv_check$passed)

        # Store results
        result_row <- data.frame(
          Model_Type = sprintf("Gaussian_%s_YEAR_RE", response_var),
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
      cat(sprintf("\n????????? LOO Comparison for %s_%s ?????????\n", response_var, cv_label))
      cv_comparison <- tryCatch({
        compare_loo_within_cv(loo_list, cv_label, convergence_map)
      }, error = function(e) {
        cat(sprintf("  ??? Comparison error: %s\n", e$message))
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
                                 sprintf("gaussian_%s_year_re_%s_%s.rds",
                                         tolower(response_var),
                                         tolower(source_name), cv_label))
        if (file.exists(source_path) && !file.exists(winner_path)) {
          file.copy(source_path, winner_path)
          cat(sprintf("  ???? Winner saved: %s\n", basename(winner_path)))
        }
      }
    }
  }
}

# ============================================================================
# FORMAL LOO COMPARISON ACROSS CV WINDOWS (same model specification)
# ============================================================================
cv_window_summary <- data.frame()
for (key_name in names(cv_window_loo_registry)) {
  comparison_path <- file.path(
    base_output_dir,
    sprintf("CV_WINDOW_LOO_HEALTH_%s.csv", tolower(gsub("__", "_", key_name)))
  )
  cv_comp <- compare_cv_windows_formal_loo(
    loo_by_cv = cv_window_loo_registry[[key_name]],
    context_label = sprintf("HEALTH_%s", key_name),
    converged_by_cv = cv_window_convergence_registry[[key_name]],
    out_csv_path = comparison_path
  )
  cv_window_summary <- rbind(
    cv_window_summary,
    data.frame(
      Model_Key = key_name,
      Status = cv_comp$status,
      Winner_CV = ifelse(is.null(cv_comp$winner_cv), NA, cv_comp$winner_cv),
      Message = cv_comp$message,
      stringsAsFactors = FALSE
    )
  )
}
if (nrow(cv_window_summary) > 0) {
  cv_summary_path <- file.path(base_output_dir, "CV_WINDOW_LOO_SUMMARY_HEALTH.csv")
  write.csv(cv_window_summary, cv_summary_path, row.names = FALSE)
  cat(sprintf("\nCV-window LOO summary saved: %s\n", cv_summary_path))
}
# ============================================================================
# GLOBAL SUMMARY
# ============================================================================

cat("\n")
cat("??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????\n")
cat("???   ALL HEALTH MODELS COMPLETED                                   ???\n")
cat("??????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????????\n")
cat("\n")

if (nrow(all_results) > 0) {
  cat(sprintf("Total models fitted: %d\n", nrow(all_results)))
  cat(sprintf("Models with successful LOO: %d\n", sum(!is.na(all_results$LOOIC))))
  cat(sprintf("Models converged: %d/%d\n",
              sum(all_results$Converged, na.rm = TRUE),
              sum(!is.na(all_results$Converged))))

  cat("\nResults summary:\n")
  print(all_results[, c("Response", "Model", "CV", "LOOIC", "SE_LOOIC", "Converged")])

  # Save combined results
  results_path <- file.path(base_output_dir, "Health_YEAR_RE_combined_results.csv")
  write.csv(all_results, results_path, row.names = FALSE)
  cat(sprintf("\n??? Combined results saved: %s\n", results_path))
} else {
  cat("??? No results to save\n")
}

cat(sprintf("\nOutput directory: %s\n", base_output_dir))
cat("\n??? Script 02b (Gaussian Health YEAR RE) completed!\n")


