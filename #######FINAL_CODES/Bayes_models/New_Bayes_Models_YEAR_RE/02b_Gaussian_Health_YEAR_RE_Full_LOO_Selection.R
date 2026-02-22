# ============================================================================
# 02b_Gaussian_Health_YEAR_RE_Full_LOO_Selection.R
# ============================================================================
# Full-grid prior sensitivity for Gaussian Health YEAR_RE models
# Runs WeaklyInformative and Informative priors for all candidate models and CV
# Writes outputs only to PS_FULLGRID namespace
# ============================================================================

setwd("C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_CODES/Bayes_models/New_Bayes_Models_YEAR_RE/")
source("00_LOO_Selection_Functions_YEAR_RE.R")

suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(loo)
  library(dplyr)
})

options(mc.cores = BRMS_CORES_DEFAULT)
set.seed(42)

response_variables <- c("HEALTH_PC1", "HEALTH_PC2")
cv_scenarios <- c("CV_02", "CV_30", "CV_ALL")
prior_scenarios <- prior_scenarios_year_re

dataset_paths_abundance <- list(
  CV_02 = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/#####output_local_PCA_CV_2_FINAL/dados_abundancia_integrados_long_format.csv"),
  CV_30 = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/#####output_local_PCA_CV_30_FINAL/dados_abundancia_integrados_long_format.csv"),
  CV_ALL = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/#####output_local_PCA_CV_all_FINAL/dados_abundancia_integrados_long_format.csv")
)

dataset_paths_final <- list(
  CV_02 = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/output_DADOS_FINAIS_PARA_MODELAGEM_cv_2/dados_finais_para_modelagem_com_ARCH.csv"),
  CV_30 = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/output_DADOS_FINAIS_PARA_MODELAGEM_cv_30/dados_finais_para_modelagem_com_ARCH.csv"),
  CV_ALL = file.path(PROJECT_ROOT, "#######FINAL_RESULTS/output_DADOS_FINAIS_PARA_MODELAGEM_cv_all/dados_finais_para_modelagem_com_ARCH.csv")
)

base_output_dir <- output_dirs_year_re_prior_sens$HEALTH_YEAR_RE
dir.create(base_output_dir, recursive = TRUE, showWarnings = FALSE)

all_results <- data.frame()
cv_window_loo_registry <- list()
cv_window_convergence_registry <- list()

# PROGRESS TRACKING
total_models <- length(cv_scenarios) * length(response_variables) * length(prior_scenarios) * length(model_combinations_year_re)
model_counter <- 0
start_global_time <- Sys.time()

for (cv_label in cv_scenarios) {

  prepared_data <- tryCatch({
    prepare_health_data_with_year(
      excel_path = RAW_HEALTH_EXCEL,
      pca_scores_path = PCA_SCORES_PATH,
      env_data_path = dataset_paths_abundance[[cv_label]],
      final_data_path = dataset_paths_final[[cv_label]],
      cv_label = cv_label
    )
  }, error = function(e) {
    cat(sprintf("Data preparation failed for %s: %s
", cv_label, e$message))
    NULL
  })
  if (is.null(prepared_data)) next

  year_check <- check_year_levels(prepared_data)
  if (!year_check$valid) next

  hab_rrtp_mean <- compute_habmerged_prior_mean(prepared_data, fallback = -1.15)

  for (response_var in response_variables) {
    data_resp <- prepared_data[!is.na(prepared_data[[response_var]]), , drop = FALSE]
    if (nrow(data_resp) == 0) next

    for (scenario_name in prior_scenarios) {
      scenario_name <- normalize_prior_scenario(scenario_name)
      scenario_tag <- tolower(gsub("[^a-z0-9]", "", scenario_name))

      out_dir <- file.path(base_output_dir, response_var, cv_label, scenario_name)
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

      loo_list <- list()
      convergence_map <- list()

      for (model_config in model_combinations_year_re) {
        model_counter <- model_counter + 1
        
        cat(sprintf("\n%s\n", paste(rep("=", 60), collapse="")))
        cat(sprintf("MODEL %d of %d (%.1f%%) | Start: %s\n", 
                    model_counter, total_models, 100 * (model_counter-1) / total_models,
                    format(Sys.time(), "%H:%M:%S")))
        cat(sprintf("Response: %s | Model: %s | CV: %s | Prior: %s\n", 
                    response_var, model_config$name, cv_label, scenario_name))
        cat(sprintf("%s\n\n", paste(rep("-", 60), collapse="")))

        formula_obj <- make_gaussian_formula_year_re(

          response_var = response_var,
          include_hab = model_config$include_hab,
          include_depth = model_config$include_depth,
          include_interaction_pca = model_config$include_interaction_pca
        )

        priors_to_use <- build_prior_set_gaussian(
          scenario_name = scenario_name,
          cv_label = cv_label,
          response_var = response_var,
          formula_obj = formula_obj,
          data = data_resp,
          hab_rrtp_mean = hab_rrtp_mean
        )

        result <- load_or_fit_model_year_re_with_prior(
          model_name = model_config$name,
          cv_name = cv_label,
          model_type = sprintf("Gaussian_%s_YEAR_RE", response_var),
          formula = formula_obj,
          data = data_resp,
          family = gaussian(),
          prior = priors_to_use,
          out_dir = out_dir,
          brms_args = brms_args_year_re,
          prior_tag = scenario_name
        )
        if (is.null(result$fit)) next

        conv_check <- check_convergence(result$fit, model_config$name, cv_label)
        loo_result <- safe_loo(result$fit, model_config$name, cv_label)
        convergence_map[[model_config$name]] <- isTRUE(conv_check$passed)
        if (is.null(loo_result)) next

        loo_list[[model_config$name]] <- loo_result
        reg_key <- sprintf("%s__%s__%s", response_var, model_config$name, scenario_name)
        if (is.null(cv_window_loo_registry[[reg_key]])) {
          cv_window_loo_registry[[reg_key]] <- list()
          cv_window_convergence_registry[[reg_key]] <- list()
        }
        cv_window_loo_registry[[reg_key]][[cv_label]] <- loo_result
        cv_window_convergence_registry[[reg_key]][[cv_label]] <- isTRUE(conv_check$passed)

        all_results <- rbind(all_results, data.frame(
          Model_Type = sprintf("Gaussian_%s_YEAR_RE", response_var),
          Response = response_var,
          Model = model_config$name,
          CV = cv_label,
          Prior_Scenario = scenario_name,
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
      ))
      
      # PROGRESS LOG - CONCLUSION
      elapsed_global <- as.numeric(difftime(Sys.time(), start_global_time, units = "mins"))
      avg_time <- elapsed_global / model_counter
      rem_models <- total_models - model_counter
      eta_h <- (avg_time * rem_models) / 60
      
      cat(sprintf("\n[PROGRESS] Model %d completed. Time: %.1f min | Total Elapsed: %.1f min\n",
                  model_counter, ifelse(is.null(result$elapsed_mins), 0, result$elapsed_mins), elapsed_global))
      cat(sprintf("[PROGRESS] Avg Time: %.1f min/model | Estimated remaining: %.1f hours\n",
                  avg_time, eta_h))
      cat(sprintf("%s\n", paste(rep("=", 60), collapse="")))
    }


      if (length(loo_list) > 1) {
        cv_comp <- compare_loo_within_cv(loo_list, paste(cv_label, response_var, scenario_name, sep = "_"), convergence_map)
        if (!is.null(cv_comp)) {
          winner_path <- file.path(out_dir, sprintf("WINNER_%s_%s_%s.rds", tolower(response_var), cv_label, scenario_name))
          src_path <- file.path(
            out_dir,
            sprintf("gaussian_%s_year_re_%s_%s_prior_%s.rds", tolower(response_var), tolower(cv_comp$winner), cv_label, scenario_tag)
          )
          if (file.exists(src_path)) {
            file.copy(src_path, winner_path, overwrite = TRUE)
          }
        }
      }
    }
  }
}

for (scenario_name in prior_scenarios) {
  scenario_name <- normalize_prior_scenario(scenario_name)
  scenario_tag_upper <- toupper(gsub("[^A-Za-z0-9]", "", scenario_name))
  scenario_keys <- names(cv_window_loo_registry)[grepl(sprintf("__%s$", scenario_name), names(cv_window_loo_registry))]

  cv_window_summary <- data.frame()
  for (reg_key in scenario_keys) {
    model_key <- sub(sprintf("__%s$", scenario_name), "", reg_key)
    comparison_path <- file.path(
      base_output_dir,
      sprintf("CV_WINDOW_LOO_HEALTH_%s_%s_PS_FULLGRID.csv", tolower(gsub("__", "_", model_key)), scenario_tag_upper)
    )
    cv_comp <- compare_cv_windows_formal_loo(
      loo_by_cv = cv_window_loo_registry[[reg_key]],
      context_label = sprintf("HEALTH_%s_%s", model_key, scenario_name),
      converged_by_cv = cv_window_convergence_registry[[reg_key]],
      out_csv_path = comparison_path
    )
    cv_window_summary <- rbind(cv_window_summary, data.frame(
      Model_Key = model_key,
      Prior_Scenario = scenario_name,
      Status = cv_comp$status,
      Winner_CV = ifelse(is.null(cv_comp$winner_cv), NA, cv_comp$winner_cv),
      Message = cv_comp$message,
      stringsAsFactors = FALSE
    ))
  }

  if (nrow(cv_window_summary) > 0) {
    cv_summary_path <- file.path(
      base_output_dir,
      sprintf("CV_WINDOW_LOO_SUMMARY_HEALTH_%s_PS_FULLGRID.csv", scenario_tag_upper)
    )
    write.csv(cv_window_summary, cv_summary_path, row.names = FALSE)
  }
}

if (nrow(all_results) > 0) {
  results_path <- file.path(base_output_dir, "Health_YEAR_RE_combined_results_all_priors_PS_FULLGRID.csv")
  write.csv(all_results, results_path, row.names = FALSE)

  prior_sensitivity <- compute_prior_sensitivity_summary(all_results, model_family_label = "Gaussian_HEALTH_YEAR_RE")
  if (nrow(prior_sensitivity) > 0) {
    prior_sensitivity_path <- file.path(base_output_dir, "PRIOR_SENSITIVITY_HEALTH_PS_FULLGRID.csv")
    write.csv(prior_sensitivity, prior_sensitivity_path, row.names = FALSE)
  }
}

cat(sprintf("Output directory: %s\\n", base_output_dir))
cat("Script 02b completed (Gaussian Health YEAR_RE prior sensitivity full-grid).\\n")
