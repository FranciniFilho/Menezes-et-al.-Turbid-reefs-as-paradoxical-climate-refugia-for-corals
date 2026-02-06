# ============================================================================
# 02_Gaussian_Health_RGR_Full_LOO_Selection.R
# ============================================================================
# Gaussian models for coral health (HEALTH_PC1, HEALTH_PC2) and growth (RGR)
# Full LOO-based model selection across HAB × DEPTH × CV scenarios
# ============================================================================
# Model Matrix: 4 combinations × 3 CVs × 3 response types = 36 models
# ============================================================================

# --- SETUP ---
rm(list = ls())
gc()

libs <- c("brms", "dplyr", "ggplot2", "cmdstanr", "bayesplot",
          "modelsummary", "loo", "tidyr", "stringr")
invisible(lapply(libs, library, character.only = TRUE))

options(mc.cores = 4)
set.seed(42)

# Source common functions
source("00_LOO_Selection_Functions.R")

# ============================================================================
# PATHS AND DIRECTORIES
# ============================================================================

cat("\n========== GAUSSIAN HEALTH/RGR MODEL SELECTION ==========\n")
cat("Initializing directories and paths...\n")

# Create output directories
for (resp in c("HEALTH_PC1", "HEALTH_PC2", "RGR")) {
  dir_name <- sprintf("Gaussian_%s", resp)
  dir.create(file.path(BASE_OUTPUT_DIR, dir_name),
             showWarnings = FALSE, recursive = TRUE)
}
cat(sprintf("Base output directory: %s\n", BASE_OUTPUT_DIR))

# ============================================================================
# MODEL DEFINITIONS
# ============================================================================

response_variables <- c("HEALTH_PC1", "HEALTH_PC2", "RGR")

cat("\nResponse variables:\n")
for (rv in response_variables) {
  cat(sprintf("  - %s\n", rv))
}

cat("\nModel combinations to fit per response:\n")
for (combo in model_combinations) {
  cat(sprintf("  - %s: HAB=%s, DEPTH=%s\n",
              combo$name,
              ifelse(combo$include_hab, "Y", "N"),
              ifelse(combo$include_depth, "Y", "N")))
}

cat(sprintf("\nTotal models: 4 combinations × 3 CVs × 3 responses = 36\n"))

# ============================================================================
# MAIN LOOP: RESPONSE VARIABLES
# ============================================================================

# Initialize results storage
all_results <- data.frame()
all_loo <- list()
model_info_list <- list()

for (response_var in response_variables) {

  cat(sprintf("\n========== RESPONSE: %s ==========\n", response_var))

  # Set output directory for this response
  output_dir <- file.path(BASE_OUTPUT_DIR, sprintf("Gaussian_%s", response_var))
  cat(sprintf("Output directory: %s\n", output_dir))

  # ========================================================================
  # INNER LOOP: CV SCENARIOS
  # ========================================================================

  for (cv_name in names(dataset_paths_health)) {

    cat(sprintf("\n--- CV: %s ---\n", cv_name))

    # Verify data path exists
    if (!file.exists(dataset_paths_health[[cv_name]])) {
      cat(sprintf("⚠ Data file not found: %s\n", dataset_paths_health[[cv_name]]))
      next
    }

    # Prepare data
    cat("Preparing Gaussian data...\n")
    data_cv <- tryCatch({
      prepare_gaussian_data(dataset_paths_health[[cv_name]], cv_name)
    }, error = function(e) {
      cat(sprintf("❌ Error preparing data: %s\n", e$message))
      NULL
    })

    if (is.null(data_cv)) {
      next
    }

    # Verify response variable exists
    if (!(response_var %in% names(data_cv))) {
      cat(sprintf("⚠ Response variable %s not found in data\n", response_var))
      next
    }

    # Remove NA values for this response
    data_cv <- data_cv[!is.na(data_cv[[response_var]]), ]
    cat(sprintf("Data after NAs removed: N=%d\n", nrow(data_cv)))

    # Data summary
    cat(sprintf("Response range: [%.2f, %.2f]\n",
                min(data_cv[[response_var]], na.rm = TRUE),
                max(data_cv[[response_var]], na.rm = TRUE)))
    cat(sprintf("REEFs: %d, HABMERGED levels: %d\n",
                length(unique(data_cv$REEF)),
                length(unique(data_cv$HABMERGED))))

    # LOO list for this CV
    loo_list_cv <- list()

    # ======================================================================
    # INNER LOOP: MODEL COMBINATIONS
    # ======================================================================

    for (combo in model_combinations) {

      model_id <- sprintf("%s_%s_%s", response_var, combo$name, cv_name)
      cat(sprintf("\n--- Model: %s ---\n", model_id))

      # Build formula
      formula <- tryCatch({
        make_gaussian_formula(
          response_var = response_var,
          include_hab = combo$include_hab,
          include_depth = combo$include_depth
        )
      }, error = function(e) {
        cat(sprintf("❌ Error building formula: %s\n", e$message))
        NULL
      })

      if (is.null(formula)) {
        next
      }

      # Fit or load model
      result <- tryCatch({
        load_or_fit_model(
          model_name = combo$name,
          cv_name = cv_name,
          model_type = sprintf("Gaussian_%s", response_var),
          formula = formula,
          data = data_cv,
          family = gaussian(),
          prior = priors_gaussian,
          out_dir = output_dir,
          brms_args = brms_args_standard
        )
      }, error = function(e) {
        cat(sprintf("❌ Error fitting model: %s\n", e$message))
        list(fit = NULL, cached = FALSE, error = e$message)
      })

      if (is.null(result$fit)) {
        next
      }

      # Check convergence
      conv_check <- tryCatch({
        check_convergence(result$fit, combo$name, cv_name)
      }, error = function(e) {
        cat(sprintf("⚠ Convergence check error: %s\n", e$message))
        list(passed = FALSE, rhat_max = NA, ess_bulk_min = NA,
             ess_tail_min = NA, n_divergent = NA)
      })

      # Calculate LOO
      loo_result <- safe_loo(result$fit, combo$name, cv_name)
      loo_list_cv[[combo$name]] <- loo_result

      # Store results
      if (!is.null(loo_result)) {
        result_row <- data.frame(
          Model_Type = sprintf("Gaussian_%s", response_var),
          Response = response_var,
          Model = combo$name,
          CV = cv_name,
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
      model_info_list[[model_id]] <- list(
        model_id = model_id,
        response = response_var,
        cv = cv_name,
        combination = combo$name,
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

    if (length(loo_list_cv) > 1) {
      cat(sprintf("\n--- Comparing models for %s_%s ---\n", response_var, cv_name))
      cv_comparison <- tryCatch({
        compare_loo_within_cv(loo_list_cv, cv_name)
      }, error = function(e) {
        cat(sprintf("⚠ Comparison error: %s\n", e$message))
        NULL
      })

      all_loo[[sprintf("%s_%s", response_var, cv_name)]] <- list(
        loo_objects = loo_list_cv,
        comparison = cv_comparison
      )

      # Save winner model copy
      if (!is.null(cv_comparison)) {
        winner_path <- file.path(output_dir,
                                 sprintf("WINNER_Gaussian_%s_%s_%s.rds",
                                         tolower(response_var),
                                         tolower(cv_comparison$winner), cv_name))
        source_path <- file.path(output_dir,
                                 sprintf("gaussian_%s_%s_%s.rds",
                                         tolower(response_var),
                                         tolower(cv_comparison$winner), cv_name))
        if (file.exists(source_path) && !file.exists(winner_path)) {
          file.copy(source_path, winner_path)
          cat(sprintf("🏆 Winner saved: %s\n", basename(winner_path)))
        }
      }
    }
  }
}

# ============================================================================
# GLOBAL COMPARISON (ACROSS CVs AND RESPONSES)
# ============================================================================

cat("\n========== GLOBAL COMPARISON ==========\n")

if (nrow(all_results) > 0) {
  # Find best model per response and CV
  summary_by_response_cv <- all_results %>%
    group_by(Response, CV) %>%
    filter(LOOIC == min(LOOIC, na.rm = TRUE)) %>%
    ungroup()

  cat("\nBest models by response and CV:\n")
  print(summary_by_response_cv[, c("Response", "CV", "Model", "LOOIC", "SE_LOOIC")])

  # Best global
  best_global <- all_results[which.min(all_results$LOOIC), ]
  cat(sprintf("\n🏆 GLOBAL WINNER: %s, %s_%s (LOOIC=%.1f, SE=%.1f)\n",
              best_global$Response, best_global$Model, best_global$CV,
              best_global$LOOIC, best_global$SE_LOOIC))

  # Best per response
  cat("\nBest models by response:\n")
  for (rv in response_variables) {
    rv_results <- all_results[all_results$Response == rv, ]
    if (nrow(rv_results) > 0) {
      best_rv <- rv_results[which.min(rv_results$LOOIC), ]
      cat(sprintf("  %s: %s_%s (LOOIC=%.1f)\n",
                  rv, best_rv$Model, best_rv$CV, best_rv$LOOIC))
    }
  }
} else {
  cat("⚠ No results to compare\n")
}

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\n========== SAVING RESULTS ==========\n")

# Save combined results table
if (nrow(all_results) > 0) {
  results_path <- file.path(BASE_OUTPUT_DIR, "Gaussian_combined_results.csv")
  write.csv(all_results, results_path, row.names = FALSE)
  cat(sprintf("✓ Combined results saved: %s\n", results_path))

  # Save LOO objects
  loo_path <- file.path(BASE_OUTPUT_DIR, "Gaussian_all_loo_objects.rds")
  saveRDS(all_loo, loo_path)
  cat(sprintf("✓ LOO objects saved: %s\n", loo_path))

  # Save model info
  info_path <- file.path(BASE_OUTPUT_DIR, "Gaussian_model_info.rds")
  saveRDS(model_info_list, info_path)
  cat(sprintf("✓ Model info saved: %s\n", info_path))
}

# Save per-response results
for (rv in response_variables) {
  rv_results <- all_results[all_results$Response == rv, ]
  if (nrow(rv_results) > 0) {
    rv_path <- file.path(BASE_OUTPUT_DIR, sprintf("Gaussian_%s_results.csv", rv))
    write.csv(rv_results, rv_path, row.names = FALSE)
    cat(sprintf("✓ %s results saved: %s\n", rv, rv_path))
  }
}

# ============================================================================
# GENERATE SUMMARY REPORT
# ============================================================================

cat("\n========== GENERATING REPORT ==========\n")

report_path <- file.path(BASE_OUTPUT_DIR, "Gaussian_summary_report.md")
sink(report_path)

cat("# Gaussian Health/RGR Model Selection Report\n\n")
cat(sprintf("**Generated:** %s\n\n", Sys.time()))
cat("**Model Family:** Gaussian\n\n")
cat("**Response Variables:**\n")
cat("- HEALTH_PC1: Health component 1 from PCA\n")
cat("- HEALTH_PC2: Health component 2 from PCA\n")
cat("- RGR: Relative Growth Rate\n\n")

cat("## Model Matrix\n\n")
cat("| Combination | HAB | DEPTH | Responses | N Models |\n")
cat("|-------------|-----|-------|-----------|----------|\n")
cat("| FULL        | Y   | Y     | 3         | 9 |\n")
cat("| NOHABITAT   | N   | Y     | 3         | 9 |\n")
cat("| NODEPTH     | Y   | N     | 3         | 9 |\n")
cat("| MINIMAL     | N   | N     | 3         | 9 |\n")
cat("| **Total**   |     |       | **3**     | **36** |\n\n")

cat("## Results Summary\n\n")
if (nrow(all_results) > 0) {
  # Table ordered by LOOIC
  top_models <- all_results[order(all_results$LOOIC), ]
  print(knitr::kable(top_models[, c("Response", "Model", "CV", "LOOIC", "SE_LOOIC", "Converged")]))
}

cat("\n## Best Model by Response and CV\n\n")
if (nrow(all_results) > 0 && exists("summary_by_response_cv")) {
  for (rv in response_variables) {
    cat(sprintf("### %s\n", rv))
    rv_summary <- summary_by_response_cv[summary_by_response_cv$Response == rv, ]
    if (nrow(rv_summary) > 0) {
      for (i in 1:nrow(rv_summary)) {
        cat(sprintf("\n#### %s\n", rv_summary$CV[i]))
        cat(sprintf("- **Model:** %s\n", rv_summary$Model[i]))
        cat(sprintf("- **LOOIC:** %.1f (SE: %.1f)\n",
                    rv_summary$LOOIC[i], rv_summary$SE_LOOIC[i]))
        cat(sprintf("- **Converged:** %s\n", ifelse(rv_summary$Converged[i], "Yes", "No")))
      }
    }
    cat("\n")
  }
}

cat("\n## Global Winner\n\n")
if (exists("best_global")) {
  cat(sprintf("- **Response:** %s\n", best_global$Response))
  cat(sprintf("- **Model:** %s\n", best_global$Model))
  cat(sprintf("- **CV:** %s\n", best_global$CV))
  cat(sprintf("- **LOOIC:** %.1f (SE: %.1f)\n",
              best_global$LOOIC, best_global$SE_LOOIC))
  cat(sprintf("- **Converged:** %s\n", ifelse(best_global$Converged, "Yes", "No")))
}

cat("\n---\n")
cat("*Generated by 02_Gaussian_Health_RGR_Full_LOO_Selection.R*\n")

sink()

cat(sprintf("✓ Report saved: %s\n", report_path))

# ============================================================================
# VISUALIZATION: LOOIC COMPARISON
# ============================================================================

if (nrow(all_results) > 0) {
  cat("\n========== GENERATING PLOTS ==========\n")

  # LOOIC comparison plot (facet by response)
  p1 <- ggplot(all_results, aes(x = Model, y = LOOIC, fill = CV)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = LOOIC - SE_LOOIC, ymax = LOOIC + SE_LOOIC),
                  position = position_dodge(0.9), width = 0.25) +
    facet_wrap(~ Response, scales = "free_y") +
    labs(title = "Gaussian Model LOOIC Comparison",
         subtitle = "Health (PC1, PC2) and Growth (RGR)",
         x = "Model Combination", y = "LOOIC (lower is better)") +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1))

  plot_path <- file.path(BASE_OUTPUT_DIR, "Gaussian_looic_comparison.png")
  ggsave(plot_path, p1, width = 12, height = 8, dpi = 300)
  cat(sprintf("✓ Plot saved: %s\n", plot_path))

  # Convergence status heatmap (facet by response)
  if (any(!is.na(all_results$Converged))) {
    conv_summary <- all_results %>%
      group_by(Response, Model, CV) %>%
      summarise(Converged = first(Converged), .groups = "drop")

    p2 <- ggplot(conv_summary, aes(x = Model, y = CV, fill = Converged)) +
      geom_tile() +
      scale_fill_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
      facet_wrap(~ Response) +
      labs(title = "Gaussian Model Convergence Status",
           x = "Model Combination", y = "CV Scenario") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    conv_plot_path <- file.path(BASE_OUTPUT_DIR, "Gaussian_convergence_status.png")
    ggsave(conv_plot_path, p2, width = 12, height = 8, dpi = 300)
    cat(sprintf("✓ Convergence plot saved: %s\n", conv_plot_path))
  }
}

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n========== GAUSSIAN SCRIPT COMPLETE ==========\n")
cat(sprintf("Total models fitted/loaded: %d\n", length(model_info_list)))
cat(sprintf("Models with successful LOO: %d\n", sum(!is.na(all_results$LOOIC))))
cat(sprintf("Models converged: %d/%d\n",
            sum(all_results$Converged, na.rm = TRUE),
            sum(!is.na(all_results$Converged))))

cat("\nBreakdown by response:\n")
for (rv in response_variables) {
  rv_count <- sum(all_results$Response == rv, na.rm = TRUE)
  rv_conv <- sum(all_results$Converged[all_results$Response == rv], na.rm = TRUE)
  cat(sprintf("  %s: %d models, %d converged\n", rv, rv_count, rv_conv))
}

cat(sprintf("\nOutput directory: %s\n", BASE_OUTPUT_DIR))

cat("\n✓ Script 02 (Gaussian Health/RGR) completed!\n")
