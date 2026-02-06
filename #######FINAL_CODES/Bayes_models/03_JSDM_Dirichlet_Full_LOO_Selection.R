# ============================================================================
# 03_JSDM_Dirichlet_Full_LOO_Selection.R
# ============================================================================
# Joint Species Distribution Model (JSDM) with Dirichlet family
# Full LOO-based model selection across HAB × DEPTH × CV scenarios
# Community composition: 6 categories (MUSSISMILIA, TURF, CCA, CYANO, MACROALGAE, OTHER)
# ============================================================================
# Model Matrix: 4 combinations × 3 CVs = 12 models
# NOTE: LOO may fail for Dirichlet; k-fold CV fallback implemented
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

cat("\n========== JSDM DIRICHLET MODEL SELECTION ==========\n")
cat("Initializing directories and paths...\n")

# Create output directory
OUTPUT_DIR_JSDM <- output_dirs$JSDM
dir.create(OUTPUT_DIR_JSDM, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Output directory: %s\n", OUTPUT_DIR_JSDM))

# ============================================================================
# JSDM-SPECIFIC FUNCTIONS
# ============================================================================

#' Calculate k-fold CV as fallback when LOO fails
#' @param fit brmsfit object
#' @param model_name Model name
#' @param cv_label CV label
#' @param K Number of folds (default 5)
#' @return kfold object or NULL
safe_kfold <- function(fit, model_name, cv_label, K = 5) {

  cat(sprintf("📊 Calculating %d-fold CV for %s_%s...\n", K, model_name, cv_label))

  kfold_result <- tryCatch({
    kfold(fit, K = K, cores = 4)
  }, error = function(e) {
    cat(sprintf("❌ k-fold CV error for %s_%s: %s\n", model_name, cv_label, e$message))
    NULL
  })

  return(kfold_result)
}


#' Compare models using either LOOIC or K-fold IC
#' @param ic_list Named list of loo or kfold objects
#' @param cv_label CV label
#' @param ic_type "LOO" or "Kfold"
#' @return data.frame with comparison
compare_ic_within_cv <- function(ic_list, cv_label, ic_type = "LOO") {

  # Remove NULLs
  ic_list <- Filter(Negate(is.null), ic_list)

  if (length(ic_list) < 2) {
    cat(sprintf("⚠ Less than 2 valid models for comparison in %s\n", cv_label))
    return(NULL)
  }

  # Compare
  comparison <- loo_compare(ic_list)

  # Extract results
  comp_df <- as.data.frame(comparison)
  comp_df$Model <- rownames(comp_df)
  comp_df$CV <- cv_label
  comp_df$Rank <- 1:nrow(comp_df)
  comp_df$IC_Type <- ic_type

  # Determine winner
  winner <- rownames(comparison)[1]
  delta_second <- abs(comparison[2, "elpd_diff"]) * 2
  se_second <- comparison[2, "se_diff"]

  # Significance
  significance <- case_when(
    delta_second > 2 * se_second ~ "Significant",
    delta_second > se_second ~ "Moderate",
    TRUE ~ "Not significant"
  )

  cat(sprintf("🏆 Winner for %s (%s): %s (Δ=%.1f, SE=%.1f, %s)\n",
              cv_label, ic_type, winner, delta_second, se_second, significance))

  return(list(
    comparison = comp_df,
    winner = winner,
    delta_ic_to_second = delta_second,
    se_diff = se_second,
    significance = significance
  ))
}


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

cat("\nCommunity composition categories:\n")
cat("  1. MUSSISMILIA (target species)\n")
cat("  2. TURF\n")
cat("  3. CCA (crustose coralline algae)\n")
cat("  4. CYANO (cyanobacteria)\n")
cat("  5. MACROALGAE\n")
cat("  6. OTHER\n")

# ============================================================================
# MAIN LOOP: CV SCENARIOS
# ============================================================================

# Initialize results storage
all_results <- data.frame()
all_ic <- list()
model_info_list <- list()
use_kfold <- FALSE  # Flag to track if we need to use k-fold

for (cv_name in names(dataset_paths_abundance)) {

  cat(sprintf("\n========== CV: %s ==========\n", cv_name))

  # Verify data path exists
  if (!file.exists(dataset_paths_abundance[[cv_name]])) {
    cat(sprintf("⚠ Data file not found: %s\n", dataset_paths_abundance[[cv_name]]))
    next
  }

  # Prepare data
  cat("Preparing JSDM data...\n")
  data_cv <- tryCatch({
    prepare_jsdm_data(dataset_paths_abundance[[cv_name]], cv_name)
  }, error = function(e) {
    cat(sprintf("❌ Error preparing data: %s\n", e$message))
    NULL
  })

  if (is.null(data_cv)) {
    next
  }

  # Data summary
  cat(sprintf("Data summary: N=%d sampling units, %d REEFs\n",
              nrow(data_cv),
              length(unique(data_cv$REEF))))

  # Check for sufficient sample size
  if (nrow(data_cv) < 30) {
    cat(sprintf("⚠ Small sample size (N=%d). Consider more conservative priors.\n", nrow(data_cv)))
  }

  # IC list for this CV
  ic_list_cv <- list()

  # ========================================================================
  # INNER LOOP: MODEL COMBINATIONS
  # ========================================================================

  for (combo in model_combinations) {

    model_id <- paste0(combo$name, "_", cv_name)
    cat(sprintf("\n--- Model: %s ---\n", model_id))

    # Build formula (pass data to check if HABMERGED and ARCH exist)
    formula <- tryCatch({
      make_jsdm_formula(
        include_hab = combo$include_hab,
        include_depth = combo$include_depth,
        include_arch_interaction = TRUE,
        data = data_cv
      )
    }, error = function(e) {
      cat(sprintf("❌ Error building formula: %s\n", e$message))
      NULL
    })

    if (is.null(formula)) {
      next
    }

    # Fit or load model (using JSDM-specific conservative settings)
    result <- tryCatch({
      load_or_fit_model(
        model_name = combo$name,
        cv_name = cv_name,
        model_type = "JSDM",
        formula = formula,
        data = data_cv,
        family = dirichlet(),
        prior = priors_jsdm,
        out_dir = OUTPUT_DIR_JSDM,
        brms_args = brms_args_jsdm  # More conservative settings
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

    # Try LOO first, fall back to k-fold if it fails
    ic_result <- safe_loo(result$fit, combo$name, cv_name, use_moment_match = TRUE)

    ic_type <- "LOO"
    if (is.null(ic_result)) {
      cat("⚠ LOO failed, trying k-fold CV...\n")
      ic_result <- safe_kfold(result$fit, combo$name, cv_name, K = 5)
      ic_type <- "Kfold"
      if (!is.null(ic_result)) {
        use_kfold <<- TRUE  # Switch global flag
      }
    }

    ic_list_cv[[combo$name]] <- ic_result

    # Store results
    if (!is.null(ic_result)) {
      # Extract IC value - handle different row names (ic, looic, etc.)
      ic_row_name <- intersect(c("ic", "looic", "kfoldic"), rownames(ic_result$estimates))[1]
      if (is.na(ic_row_name)) {
        # Fallback: try to get the first row that contains "IC" in the name
        ic_row_name <- grep("IC", rownames(ic_result$estimates), value = TRUE, ignore.case = TRUE)[1]
        if (is.na(ic_row_name)) {
          # Last resort: use the first row
          ic_row_name <- rownames(ic_result$estimates)[1]
        }
      }

      result_row <- data.frame(
        Model_Type = "JSDM_Dirichlet",
        Response = "Community_Composition",
        Model = combo$name,
        CV = cv_name,
        IC_Type = ic_type,
        IC = ic_result$estimates[ic_row_name, "Estimate"],
        SE_IC = ic_result$estimates[ic_row_name, "SE"],
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
      ic = ic_result,
      ic_type = ic_type,
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

  if (length(ic_list_cv) > 1) {
    # Determine which IC type to use for comparison
    cv_ic_type <- ifelse(use_kfold, "Kfold", "LOO")

    cat(sprintf("\n--- Comparing models for %s (%s) ---\n", cv_name, cv_ic_type))
    cv_comparison <- tryCatch({
      compare_ic_within_cv(ic_list_cv, cv_name, ic_type = cv_ic_type)
    }, error = function(e) {
      cat(sprintf("⚠ Comparison error: %s\n", e$message))
      NULL
    })

    all_ic[[cv_name]] <- list(
      ic_objects = ic_list_cv,
      comparison = cv_comparison,
      ic_type = cv_ic_type
    )

    # Save winner model copy
    if (!is.null(cv_comparison)) {
      winner_path <- file.path(OUTPUT_DIR_JSDM,
                               sprintf("WINNER_JSDM_%s_%s.rds",
                                       tolower(cv_comparison$winner), cv_name))
      source_path <- file.path(OUTPUT_DIR_JSDM,
                               sprintf("jsdm_%s_%s.rds",
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
  # Separate LOO and Kfold results
  loo_results <- all_results[all_results$IC_Type == "LOO", ]
  kfold_results <- all_results[all_results$IC_Type == "Kfold", ]

  if (nrow(loo_results) > 0) {
    # Find best model per CV (LOO)
    summary_by_cv_loo <- loo_results %>%
      group_by(CV) %>%
      filter(IC == min(IC, na.rm = TRUE)) %>%
      ungroup()

    cat("\nBest models by CV (LOOIC):\n")
    print(summary_by_cv_loo[, c("CV", "Model", "IC", "SE_IC")])
  }

  if (nrow(kfold_results) > 0) {
    cat("\n⚠ Note: Some models used K-fold CV instead of LOO\n")
    cat("K-fold results:\n")
    print(kfold_results[, c("CV", "Model", "IC", "SE_IC")])
  }

  # Best global (prefer LOO, fallback to Kfold)
  if (nrow(loo_results) > 0) {
    best_global <- loo_results[which.min(loo_results$IC), ]
  } else {
    best_global <- kfold_results[which.min(kfold_results$IC), ]
  }

  cat(sprintf("\n🏆 GLOBAL WINNER: %s_%s (%s=%.1f, SE=%.1f)\n",
              best_global$Model, best_global$CV,
              best_global$IC_Type, best_global$IC, best_global$SE_IC))
} else {
  cat("⚠ No results to compare\n")
}

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\n========== SAVING RESULTS ==========\n")

# Save comparison table
if (nrow(all_results) > 0) {
  results_path <- file.path(OUTPUT_DIR_JSDM, "jsdm_full_comparison_results.csv")
  write.csv(all_results, results_path, row.names = FALSE)
  cat(sprintf("✓ Results saved: %s\n", results_path))

  # Save IC objects
  ic_path <- file.path(OUTPUT_DIR_JSDM, "jsdm_all_ic_objects.rds")
  saveRDS(all_ic, ic_path)
  cat(sprintf("✓ IC objects saved: %s\n", ic_path))

  # Save model info
  info_path <- file.path(OUTPUT_DIR_JSDM, "jsdm_model_info.rds")
  saveRDS(model_info_list, info_path)
  cat(sprintf("✓ Model info saved: %s\n", info_path))
}

# ============================================================================
# GENERATE SUMMARY REPORT
# ============================================================================

cat("\n========== GENERATING REPORT ==========\n")

report_path <- file.path(OUTPUT_DIR_JSDM, "jsdm_summary_report.md")
sink(report_path)

cat("# JSDM Dirichlet Model Selection Report\n\n")
cat(sprintf("**Generated:** %s\n\n", Sys.time()))
cat("**Model Family:** Dirichlet (Joint Species Distribution Model)\n\n")
cat("**Response:** Community composition proportions (6 categories)\n\n")
cat("**Categories:**\n")
cat("1. MUSSISMILIA (target species: *Mussismilia hispida*)\n")
cat("2. TURF (turf algae)\n")
cat("3. CCA (crustose coralline algae)\n")
cat("4. CYANO (cyanobacteria)\n")
cat("5. MACROALGAE\n")
cat("6. OTHER\n\n")

cat("## Model Matrix\n\n")
cat("| Combination | HAB | DEPTH | N Models |\n")
cat("|-------------|-----|-------|----------|\n")
cat("| FULL        | Y   | Y     | 3 |\n")
cat("| NOHABITAT   | N   | Y     | 3 |\n")
cat("| NODEPTH     | Y   | N     | 3 |\n")
cat("| MINIMAL     | N   | N     | 3 |\n")
cat("| **Total**   |     |       | **12** |\n\n")

cat("## Notes\n\n")
cat("- **Sample Size:** JSDM models use aggregated sampling units (SITE × HAB)\n")
cat("- **LOO Limitation:** LOO may fail for Dirichlet models due to small N\n")
cat("- **K-fold Fallback:** 5-fold CV used when LOO fails\n")
cat("- **Conservative Settings:** iter=6000, adapt_delta=0.99\n\n")

cat("## Results Summary\n\n")
if (nrow(all_results) > 0) {
  # Table ordered by IC
  top_models <- all_results[order(all_results$IC), ]
  print(knitr::kable(top_models[, c("Model", "CV", "IC_Type", "IC", "SE_IC", "Converged")]))
}

cat("\n## Best Model by CV\n\n")
if (nrow(all_results) > 0) {
  for (cv in unique(all_results$CV)) {
    cv_results <- all_results[all_results$CV == cv, ]
    if (nrow(cv_results) > 0) {
      best_cv <- cv_results[which.min(cv_results$IC), ]
      cat(sprintf("### %s\n", cv))
      cat(sprintf("- **Model:** %s\n", best_cv$Model))
      cat(sprintf("- **IC Type:** %s\n", best_cv$IC_Type))
      cat(sprintf("- **IC:** %.1f (SE: %.1f)\n",
                  best_cv$IC, best_cv$SE_IC))
      cat(sprintf("- **Converged:** %s\n", ifelse(best_cv$Converged, "Yes", "No")))
      cat("\n")
    }
  }
}

cat("\n## Global Winner\n\n")
if (exists("best_global")) {
  cat(sprintf("- **Model:** %s\n", best_global$Model))
  cat(sprintf("- **CV:** %s\n", best_global$CV))
  cat(sprintf("- **IC Type:** %s\n", best_global$IC_Type))
  cat(sprintf("- **IC:** %.1f (SE: %.1f)\n",
              best_global$IC, best_global$SE_IC))
  cat(sprintf("- **Converged:** %s\n", ifelse(best_global$Converged, "Yes", "No")))
}

cat("\n---\n")
cat("*Generated by 03_JSDM_Dirichlet_Full_LOO_Selection.R*\n")

sink()

cat(sprintf("✓ Report saved: %s\n", report_path))

# ============================================================================
# VISUALIZATION: IC COMPARISON
# ============================================================================

if (nrow(all_results) > 0) {
  cat("\n========== GENERATING PLOTS ==========\n")

  # IC comparison plot
  ic_label <- ifelse(use_kfold, "IC (LOOIC/KfoldIC)", "LOOIC")

  p1 <- ggplot(all_results, aes(x = Model, y = IC, fill = CV)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = IC - SE_IC, ymax = IC + SE_IC),
                  position = position_dodge(0.9), width = 0.25) +
    labs(title = "JSDM Dirichlet Model IC Comparison",
         subtitle = "Community composition (6 categories)",
         x = "Model Combination",
         y = sprintf("%s (lower is better)", ic_label)) +
    theme_minimal() +
    theme(legend.position = "bottom")

  plot_path <- file.path(OUTPUT_DIR_JSDM, "jsdm_ic_comparison.png")
  ggsave(plot_path, p1, width = 10, height = 6, dpi = 300)
  cat(sprintf("✓ Plot saved: %s\n", plot_path))

  # IC type indicator
  if (length(unique(all_results$IC_Type)) > 1) {
    p2 <- ggplot(all_results, aes(x = CV, fill = IC_Type)) +
      geom_bar(position = "fill") +
      labs(title = "Proportion of LOO vs K-fold by CV",
           x = "CV Scenario", y = "Proportion") +
      theme_minimal()

    ic_type_path <- file.path(OUTPUT_DIR_JSDM, "jsdm_ic_type_distribution.png")
    ggsave(ic_type_path, p2, width = 6, height = 4, dpi = 300)
    cat(sprintf("✓ IC type plot saved: %s\n", ic_type_path))
  }

  # Convergence status plot
  if (any(!is.na(all_results$Converged))) {
    conv_summary <- all_results %>%
      group_by(Model, CV) %>%
      summarise(Converged = first(Converged), .groups = "drop")

    p3 <- ggplot(conv_summary, aes(x = Model, y = CV, fill = Converged)) +
      geom_tile() +
      scale_fill_manual(values = c("TRUE" = "green", "FALSE" = "red")) +
      labs(title = "JSDM Model Convergence Status",
           x = "Model Combination", y = "CV Scenario") +
      theme_minimal()

    conv_plot_path <- file.path(OUTPUT_DIR_JSDM, "jsdm_convergence_status.png")
    ggsave(conv_plot_path, p3, width = 8, height = 4, dpi = 300)
    cat(sprintf("✓ Convergence plot saved: %s\n", conv_plot_path))
  }
}

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n========== JSDM SCRIPT COMPLETE ==========\n")
cat(sprintf("Total models fitted/loaded: %d\n", length(model_info_list)))
cat(sprintf("Models with successful IC calculation: %d\n", nrow(all_results)))
cat(sprintf("Models using LOO: %d\n", sum(all_results$IC_Type == "LOO")))
cat(sprintf("Models using K-fold: %d\n", sum(all_results$IC_Type == "Kfold")))
cat(sprintf("Models converged: %d/%d\n",
            sum(all_results$Converged, na.rm = TRUE),
            sum(!is.na(all_results$Converged))))
cat(sprintf("\nOutput directory: %s\n", OUTPUT_DIR_JSDM))

cat("\n✓ Script 03 (JSDM Dirichlet) completed!\n")
