# ============================================================================
# 01_ZOIB_Abundance_Full_LOO_Selection.R
# ============================================================================
# Zero-One Inflated Beta (ZOIB) models for coral abundance (M. hispida)
# Full LOO-based model selection across HAB × DEPTH × CV scenarios
# ============================================================================
# Model Matrix: 4 combinations × 3 CVs = 12 models
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

cat("\n========== ZOIB ABUNDANCE MODEL SELECTION ==========\n")
cat("Initializing directories and paths...\n")

# Create output directory
OUTPUT_DIR_ZOIB <- output_dirs$ZOIB
dir.create(OUTPUT_DIR_ZOIB, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Output directory: %s\n", OUTPUT_DIR_ZOIB))

# ============================================================================
# MODEL DEFINITIONS
# ============================================================================

cat("\nModel combinations to fit:\n")
for (combo in model_combinations) {
  cat(sprintf("  - %s: HAB=%s, DEPTH=%s\n",
              combo$name,
              ifelse(combo$include_hab, "Y", "N"),
              ifelse(combo$include_depth, "Y", "N")))
}

# ============================================================================
# MAIN LOOP: CV SCENARIOS
# ============================================================================

# Initialize results storage
all_results <- data.frame()
all_loo <- list()
model_info_list <- list()  # Store detailed model info

for (cv_name in names(dataset_paths_abundance)) {

  cat(sprintf("\n========== CV: %s ==========\n", cv_name))

  # Verify data path exists
  if (!file.exists(dataset_paths_abundance[[cv_name]])) {
    cat(sprintf("⚠ Data file not found: %s\n", dataset_paths_abundance[[cv_name]]))
    next
  }

  # Prepare data
  cat("Preparing ZOIB data...\n")
  data_cv <- tryCatch({
    prepare_zoib_data(dataset_paths_abundance[[cv_name]], cv_name)
  }, error = function(e) {
    cat(sprintf("❌ Error preparing data: %s\n", e$message))
    NULL
  })

  if (is.null(data_cv)) {
    next
  }

  # Data summary
  cat(sprintf("Data summary: N=%d, REEFs=%d, HABMERGED levels=%d\n",
              nrow(data_cv),
              length(unique(data_cv$REEF)),
              length(unique(data_cv$HABMERGED))))

  # LOO list for this CV
  loo_list_cv <- list()

  # ========================================================================
  # INNER LOOP: MODEL COMBINATIONS
  # ========================================================================

  for (combo in model_combinations) {

    model_id <- paste0(combo$name, "_", cv_name)
    cat(sprintf("\n--- Modelo: %s ---\n", model_id))

    # Build formula
    formula <- tryCatch({
      make_zoib_formula(
        include_hab = combo$include_hab,
        include_depth = combo$include_depth,
        include_arch_interaction = TRUE
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
        model_type = "ZOIB",
        formula = formula,
        data = data_cv,
        family = zero_one_inflated_beta(),
        prior = priors_zoib,
        out_dir = OUTPUT_DIR_ZOIB,
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
        Model_Type = "ZOIB",
        Response = "COVER_PROP",
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

  # ========================================================================
  # COMPARE MODELS WITHIN CV
  # ========================================================================

  if (length(loo_list_cv) > 1) {
    cat(sprintf("\n--- Comparing models for %s ---\n", cv_name))
    cv_comparison <- tryCatch({
      compare_loo_within_cv(loo_list_cv, cv_name)
    }, error = function(e) {
      cat(sprintf("⚠ Comparison error: %s\n", e$message))
      NULL
    })

    all_loo[[cv_name]] <- list(
      loo_objects = loo_list_cv,
      comparison = cv_comparison
    )

    # Save winner model copy
    if (!is.null(cv_comparison)) {
      winner_path <- file.path(OUTPUT_DIR_ZOIB,
                               sprintf("WINNER_ZOIB_%s_%s.rds",
                                       tolower(cv_comparison$winner), cv_name))
      source_path <- file.path(OUTPUT_DIR_ZOIB,
                               sprintf("zoib_%s_%s.rds",
                                       tolower(cv_comparison$winner), cv_name))
      if (file.exists(source_path) && !file.exists(winner_path)) {
        file.copy(source_path, winner_path)
        cat(sprintf("🏆 Winner saved: %s\n", basename(winner_path)))
      }
    }
  }
}

# ============================================================================
# GLOBAL COMPARISON (ACROSS CVs)
# ============================================================================

cat("\n========== GLOBAL COMPARISON ==========\n")

if (nrow(all_results) > 0) {
  # Find best model per CV
  summary_by_cv <- all_results %>%
    group_by(CV) %>%
    filter(LOOIC == min(LOOIC, na.rm = TRUE)) %>%
    ungroup()

  cat("\nBest models by CV:\n")
  print(summary_by_cv[, c("CV", "Model", "LOOIC", "SE_LOOIC")])

  # Best global
  best_global <- all_results[which.min(all_results$LOOIC), ]
  cat(sprintf("\n🏆 GLOBAL WINNER: %s_%s (LOOIC=%.1f, SE=%.1f)\n",
              best_global$Model, best_global$CV,
              best_global$LOOIC, best_global$SE_LOOIC))
} else {
  cat("⚠ No results to compare\n")
}

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\n========== SAVING RESULTS ==========\n")

# Save comparison table
if (nrow(all_results) > 0) {
  results_path <- file.path(OUTPUT_DIR_ZOIB, "zoib_full_comparison_results.csv")
  write.csv(all_results, results_path, row.names = FALSE)
  cat(sprintf("✓ Results saved: %s\n", results_path))

  # Save LOO objects
  loo_path <- file.path(OUTPUT_DIR_ZOIB, "zoib_all_loo_objects.rds")
  saveRDS(all_loo, loo_path)
  cat(sprintf("✓ LOO objects saved: %s\n", loo_path))

  # Save model info
  info_path <- file.path(OUTPUT_DIR_ZOIB, "zoib_model_info.rds")
  saveRDS(model_info_list, info_path)
  cat(sprintf("✓ Model info saved: %s\n", info_path))
}

# ============================================================================
# GENERATE SUMMARY REPORT
# ============================================================================

cat("\n========== GENERATING REPORT ==========\n")

report_path <- file.path(OUTPUT_DIR_ZOIB, "zoib_summary_report.md")
sink(report_path)

cat("# ZOIB Abundance Model Selection Report\n\n")
cat(sprintf("**Generated:** %s\n\n", Sys.time()))
cat("**Model Family:** Zero-One Inflated Beta (ZOIB)\n\n")
cat("**Response Variable:** COVER_PROP (Mussismilia hispida coverage proportion)\n\n")

cat("## Model Matrix\n\n")
cat("| Combination | HAB | DEPTH | N Models |\n")
cat("|-------------|-----|-------|----------|\n")
cat("| FULL        | Y   | Y     | 3 |\n")
cat("| NOHABITAT   | N   | Y     | 3 |\n")
cat("| NODEPTH     | Y   | N     | 3 |\n")
cat("| MINIMAL     | N   | N     | 3 |\n")
cat("| **Total**   |     |       | **12** |\n\n")

cat("## Results Summary\n\n")
if (nrow(all_results) > 0) {
  print(knitr::kable(all_results[order(all_results$LOOIC),
                                  c("Model", "CV", "LOOIC", "SE_LOOIC",
                                    "Converged", "Rhat_Max", "ESS_Bulk_Min")]))
}

cat("\n## Best Model by CV\n\n")
if (nrow(all_results) > 0 && exists("summary_by_cv")) {
  for (i in 1:nrow(summary_by_cv)) {
    cat(sprintf("### %s\n", summary_by_cv$CV[i]))
    cat(sprintf("- **Model:** %s\n", summary_by_cv$Model[i]))
    cat(sprintf("- **LOOIC:** %.1f (SE: %.1f)\n",
                summary_by_cv$LOOIC[i], summary_by_cv$SE_LOOIC[i]))
    cat(sprintf("- **Converged:** %s\n", ifelse(summary_by_cv$Converged[i], "Yes", "No")))
    cat("\n")
  }
}

cat("\n## Global Winner\n\n")
if (exists("best_global")) {
  cat(sprintf("- **Model:** %s\n", best_global$Model))
  cat(sprintf("- **CV:** %s\n", best_global$CV))
  cat(sprintf("- **LOOIC:** %.1f (SE: %.1f)\n",
              best_global$LOOIC, best_global$SE_LOOIC))
  cat(sprintf("- **Converged:** %s\n", ifelse(best_global$Converged, "Yes", "No")))
}

cat("\n---\n")
cat("*Generated by 01_ZOIB_Abundance_Full_LOO_Selection.R*\n")

sink()

cat(sprintf("✓ Report saved: %s\n", report_path))

# ============================================================================
# VISUALIZATION: LOOIC COMPARISON
# ============================================================================

if (nrow(all_results) > 0) {
  cat("\n========== GENERATING PLOTS ==========\n")

  # LOOIC comparison plot
  p1 <- ggplot(all_results, aes(x = Model, y = LOOIC, fill = CV)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = LOOIC - SE_LOOIC, ymax = LOOIC + SE_LOOIC),
                  position = position_dodge(0.9), width = 0.25) +
    labs(title = "ZOIB Model LOOIC Comparison",
         subtitle = "Mussismilia hispida abundance (COVER_PROP)",
         x = "Model Combination", y = "LOOIC (lower is better)") +
    theme_minimal() +
    theme(legend.position = "bottom")

  plot_path <- file.path(OUTPUT_DIR_ZOIB, "zoib_looic_comparison.png")
  ggsave(plot_path, p1, width = 10, height = 6, dpi = 300)
  cat(sprintf("✓ Plot saved: %s\n", plot_path))

  # Convergence status plot
  if (any(!is.na(all_results$Converged))) {
    conv_summary <- all_results %>%
      group_by(Model, CV) %>%
      summarise(Converged = first(Converged), .groups = "drop")

    p2 <- ggplot(conv_summary, aes(x = Model, y = CV, fill = Converged)) +
      geom_tile() +
      scale_fill_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
      labs(title = "ZOIB Model Convergence Status",
           x = "Model Combination", y = "CV Scenario") +
      theme_minimal()

    conv_plot_path <- file.path(OUTPUT_DIR_ZOIB, "zoib_convergence_status.png")
    ggsave(conv_plot_path, p2, width = 8, height = 4, dpi = 300)
    cat(sprintf("✓ Convergence plot saved: %s\n", conv_plot_path))
  }
}

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n========== ZOIB SCRIPT COMPLETE ==========\n")
cat(sprintf("Total models fitted/loaded: %d\n", length(model_info_list)))
cat(sprintf("Models with successful LOO: %d\n", sum(!is.na(all_results$LOOIC))))
cat(sprintf("Models converged: %d/%d\n",
            sum(all_results$Converged, na.rm = TRUE),
            sum(!is.na(all_results$Converged))))
cat(sprintf("\nOutput directory: %s\n", OUTPUT_DIR_ZOIB))

cat("\n✓ Script 01 (ZOIB Abundance) completed!\n")
