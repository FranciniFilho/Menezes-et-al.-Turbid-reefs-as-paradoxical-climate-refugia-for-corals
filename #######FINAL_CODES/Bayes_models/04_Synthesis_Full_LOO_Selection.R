# ============================================================================
# 04_Synthesis_Full_LOO_Selection.R
# ============================================================================
# Synthesis of all LOO model selection results
# Creates master comparison, visualizations, and final report
# ============================================================================
# Loads results from:
#   - 01_ZOIB_Abundance_Full_LOO_Selection.R (12 models)
#   - 02_Gaussian_Health_RGR_Full_LOO_Selection.R (36 models)
#   - 03_JSDM_Dirichlet_Full_LOO_Selection.R (12 models)
# Total: 60 models
# ============================================================================

# --- SETUP ---
rm(list = ls())
gc()

libs <- c("brms", "dplyr", "ggplot2", "tidyr", "readr", "stringr",
          "viridis", "RColorBrewer", "gridExtra", "knitr")
invisible(lapply(libs, library, character.only = TRUE))

# Source common functions for paths
source("00_LOO_Selection_Functions.R")

# ============================================================================
# PATHS AND DIRECTORIES
# ============================================================================

cat("\n========== FULL LOO SELECTION SYNTHESIS ==========\n")
cat("Loading results from all model selection scripts...\n")

# Base output directory
SYNTHESIS_DIR <- file.path(BASE_OUTPUT_DIR, "synthesis/")
dir.create(SYNTHESIS_DIR, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Synthesis output directory: %s\n", SYNTHESIS_DIR))

# ============================================================================
# LOAD RESULTS FROM ALL SCRIPTS
# ============================================================================

cat("\n========== LOADING RESULTS ==========\n")

# File paths
zoib_results_path <- file.path(output_dirs$ZOIB, "zoib_full_comparison_results.csv")
zoib_info_path <- file.path(output_dirs$ZOIB, "zoib_model_info.rds")

gaussian_results_path <- file.path(BASE_OUTPUT_DIR, "Gaussian_combined_results.csv")
gaussian_info_path <- file.path(BASE_OUTPUT_DIR, "Gaussian_model_info.rds")

jsdm_results_path <- file.path(output_dirs$JSDM, "jsdm_full_comparison_results.csv")
jsdm_info_path <- file.path(output_dirs$JSDM, "jsdm_model_info.rds")

# Load ZOIB results
zoib_results <- NULL
zoib_info <- NULL
if (file.exists(zoib_results_path)) {
  zoib_results <- read.csv(zoib_results_path, stringsAsFactors = FALSE)
  cat(sprintf("✓ ZOIB results loaded: %d models\n", nrow(zoib_results)))
} else {
  cat("⚠ ZOIB results not found\n")
}

if (file.exists(zoib_info_path)) {
  zoib_info <- readRDS(zoib_info_path)
}

# Load Gaussian results
gaussian_results <- NULL
gaussian_info <- NULL
if (file.exists(gaussian_results_path)) {
  gaussian_results <- read.csv(gaussian_results_path, stringsAsFactors = FALSE)
  cat(sprintf("✓ Gaussian results loaded: %d models\n", nrow(gaussian_results)))
} else {
  cat("⚠ Gaussian results not found\n")
}

if (file.exists(gaussian_info_path)) {
  gaussian_info <- readRDS(gaussian_info_path)
}

# Load JSDM results
jsdm_results <- NULL
jsdm_info <- NULL
if (file.exists(jsdm_results_path)) {
  jsdm_results <- read.csv(jsdm_results_path, stringsAsFactors = FALSE)
  cat(sprintf("✓ JSDM results loaded: %d models\n", nrow(jsdm_results)))
} else {
  cat("⚠ JSDM results not found\n")
}

if (file.exists(jsdm_info_path)) {
  jsdm_info <- readRDS(jsdm_info_path)
}

# ============================================================================
# CREATE MASTER COMPARISON TABLE
# ============================================================================

cat("\n========== CREATING MASTER COMPARISON ==========\n")

# Combine all results
master_table <- bind_rows(
  zoib_results %>% mutate(IC = LOOIC, SE_IC = SE_LOOIC) %>% select(-LOOIC, -SE_LOOIC),
  gaussian_results %>% mutate(IC = LOOIC, SE_IC = SE_LOOIC) %>% select(-LOOIC, -SE_LOOIC),
  jsdm_results
)

if (nrow(master_table) > 0) {
  cat(sprintf("✓ Master table created: %d total models\n", nrow(master_table)))

  # Add rankings within each model type and CV
  master_table <- master_table %>%
    group_by(Model_Type, CV) %>%
    arrange(IC) %>%
    mutate(Rank = row_number()) %>%
    ungroup()

  # Add overall rank
  master_table <- master_table %>%
    arrange(IC) %>%
    mutate(Overall_Rank = row_number())

  # Add difference from best within group
  master_table <- master_table %>%
    group_by(Model_Type, CV) %>%
    mutate(Delta_IC = IC - min(IC, na.rm = TRUE)) %>%
    ungroup()

  # Save master table
  master_path <- file.path(SYNTHESIS_DIR, "master_comparison_all_models.csv")
  write.csv(master_table, master_path, row.names = FALSE)
  cat(sprintf("✓ Master table saved: %s\n", master_path))
} else {
  cat("⚠ No results to synthesize\n")
}

# ============================================================================
# CONSISTENCY ANALYSIS
# ============================================================================

cat("\n========== CONSISTENCY ANALYSIS ==========\n")

if (nrow(master_table) > 0) {
  # Analyze which model types tend to win
  winner_by_cv <- master_table %>%
    group_by(Model_Type, CV) %>%
    filter(Rank == 1) %>%
    ungroup()

  # Count wins by combination
  combination_wins <- master_table %>%
    group_by(Model, Model_Type) %>%
    summarise(
      Wins = sum(Rank == 1, na.rm = TRUE),
      Total = n(),
      Win_Rate = Wins / Total,
      .groups = "drop"
    ) %>%
    arrange(desc(Wins))

  cat("\nModel combination win rates:\n")
  print(combination_wins)

  # HAB importance: Do models with HAB win more often?
  hab_importance <- master_table %>%
    mutate(Has_HAB = ifelse(Model %in% c("FULL", "NODEPTH"), "Yes", "No")) %>%
    group_by(Model_Type, Has_HAB) %>%
    summarise(
      Mean_Rank = mean(Rank, na.rm = TRUE),
      Median_Rank = median(Rank, na.rm = TRUE),
      Wins = sum(Rank == 1, na.rm = TRUE),
      .groups = "drop"
    )

  cat("\nHAB importance analysis:\n")
  print(hab_importance)

  # DEPTH importance: Do models with DEPTH win more often?
  depth_importance <- master_table %>%
    mutate(Has_DEPTH = ifelse(Model %in% c("FULL", "NOHABITAT"), "Yes", "No")) %>%
    group_by(Model_Type, Has_DEPTH) %>%
    summarise(
      Mean_Rank = mean(Rank, na.rm = TRUE),
      Median_Rank = median(Rank, na.rm = TRUE),
      Wins = sum(Rank == 1, na.rm = TRUE),
      .groups = "drop"
    )

  cat("\nDEPTH importance analysis:\n")
  print(depth_importance)

  # Save consistency analysis
  consistency_path <- file.path(SYNTHESIS_DIR, "consistency_analysis.csv")
  consistency_results <- list(
    combination_wins = combination_wins,
    hab_importance = hab_importance,
    depth_importance = depth_importance
  )
  saveRDS(consistency_results, consistency_path)
  cat(sprintf("✓ Consistency analysis saved: %s\n", consistency_path))
}

# ============================================================================
# GENERATE VISUALIZATIONS
# ============================================================================

cat("\n========== GENERATING VISUALIZATIONS ==========\n")

if (nrow(master_table) > 0) {

  # -----------------------------------------------------------------------
  # 1. LOOIC/IC Heatmap (Model Type × Model Combination × CV)
  # -----------------------------------------------------------------------

  # Normalize IC within model type for better visualization
  heatmap_data <- master_table %>%
    group_by(Model_Type) %>%
    mutate(Normalized_IC = scale(IC)[,1]) %>%
    ungroup() %>%
    mutate(Model_Label = paste0(Model, "\n(", CV, ")"))

  p_heatmap <- ggplot(heatmap_data, aes(x = Model, y = CV, fill = Normalized_IC)) +
    geom_tile(color = "white") +
    scale_fill_viridis_c(option = "inferno", name = "Normalized IC\n(lower = better)") +
    facet_wrap(~ Model_Type, scales = "free_x") +
    labs(title = "Model Selection Heatmap",
         subtitle = "IC values normalized within each model type",
         x = "Model Combination",
         y = "CV Scenario") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    )

  heatmap_path <- file.path(SYNTHESIS_DIR, "fig_ic_heatmap.png")
  ggsave(heatmap_path, p_heatmap, width = 14, height = 8, dpi = 300)
  cat(sprintf("✓ Heatmap saved: %s\n", heatmap_path))

  # -----------------------------------------------------------------------
  # 2. Winner Barplot by Model Type
  # -----------------------------------------------------------------------

  winners <- master_table %>%
    group_by(Model_Type, CV) %>%
    filter(Rank == 1) %>%
    ungroup() %>%
    mutate(Winner_Label = paste0(Model, " (", CV, ")"))

  p_winners <- ggplot(winners, aes(x = reorder(Winner_Label, -IC), y = IC, fill = Model_Type)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = IC - SE_IC, ymax = IC + SE_IC), width = 0.25) +
    coord_flip() +
    labs(title = "Winning Models by Type and CV",
         subtitle = "Best model (lowest IC) for each model type × CV combination",
         x = "Model (CV)",
         y = "IC (lower is better)") +
    theme_minimal() +
    theme(legend.position = "bottom")

  winners_path <- file.path(SYNTHESIS_DIR, "fig_winner_barplot.png")
  ggsave(winners_path, p_winners, width = 12, height = 8, dpi = 300)
  cat(sprintf("✓ Winner barplot saved: %s\n", winners_path))

  # -----------------------------------------------------------------------
  # 3. IC Distribution by Model Combination
  # -----------------------------------------------------------------------

  p_boxplot <- ggplot(master_table, aes(x = Model, y = IC, fill = Model_Type)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitterdodge(0.2), alpha = 0.5) +
    facet_wrap(~ Model_Type, scales = "free_y") +
    labs(title = "IC Distribution by Model Combination",
         subtitle = "Spread of IC values across CVs for each model type",
         x = "Model Combination",
         y = "IC") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  boxplot_path <- file.path(SYNTHESIS_DIR, "fig_ic_distribution.png")
  ggsave(boxplot_path, p_boxplot, width = 14, height = 10, dpi = 300)
  cat(sprintf("✓ Distribution plot saved: %s\n", boxplot_path))

  # -----------------------------------------------------------------------
  # 4. Convergence Status Summary
  # -----------------------------------------------------------------------

  conv_summary <- master_table %>%
    group_by(Model_Type, Model) %>%
    summarise(
      Converged_Rate = mean(Converged == TRUE, na.rm = TRUE),
      Total = n(),
      .groups = "drop"
    ) %>%
    mutate(Converged_Status = ifelse(Converged_Rate == 1, "All",
                                     ifelse(Converged_Rate > 0.5, "Partial", "None")))

  p_conv <- ggplot(conv_summary, aes(x = Model, y = Converged_Rate, fill = Model_Type)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    facet_wrap(~ Model_Type) +
    labs(title = "Model Convergence Rate",
         subtitle = "Proportion of models that converged (Rhat < 1.01, ESS > 400)",
         x = "Model Combination",
         y = "Convergence Rate") +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  conv_path <- file.path(SYNTHESIS_DIR, "fig_convergence_summary.png")
  ggsave(conv_path, p_conv, width = 12, height = 8, dpi = 300)
  cat(sprintf("✓ Convergence summary saved: %s\n", conv_path))

  # -----------------------------------------------------------------------
  # 5. CV Comparison: Which CV scenario performs best?
  # -----------------------------------------------------------------------

  cv_summary <- master_table %>%
    group_by(Model_Type, CV) %>%
    summarise(
      Mean_IC = mean(IC, na.rm = TRUE),
      Median_IC = median(IC, na.rm = TRUE),
      Best_IC = min(IC, na.rm = TRUE),
      .groups = "drop"
    )

  p_cv <- ggplot(cv_summary, aes(x = CV, y = Mean_IC, fill = Model_Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = Best_IC, ymax = Best_IC),
                  position = position_dodge(0.9), width = 0.25) +
    labs(title = "IC Summary by CV Scenario",
         subtitle = "Bars: Mean IC | Lines: Best (lowest) IC",
         x = "CV Scenario",
         y = "Mean IC") +
    theme_minimal() +
    theme(legend.position = "bottom")

  cv_path <- file.path(SYNTHESIS_DIR, "fig_cv_comparison.png")
  ggsave(cv_path, p_cv, width = 10, height = 6, dpi = 300)
  cat(sprintf("✓ CV comparison saved: %s\n", cv_path))
}

# ============================================================================
# FINAL REPORT GENERATION
# ============================================================================

cat("\n========== GENERATING FINAL REPORT ==========\n")

report_path <- file.path(SYNTHESIS_DIR, "final_report.md")
sink(report_path)

cat("# Full LOO Model Selection - Final Report\n\n")
cat(sprintf("**Generated:** %s\n\n", Sys.time()))
cat("**Project:** Assessing Environmental Refugia for Coral Health\n")
cat("**Species:** *Mussismilia hispida*\n")
cat("**Location:** Abrolhos Bank, Brazil\n\n")

cat("## Executive Summary\n\n")
cat(sprintf("This report synthesizes results from **%d** Bayesian models fitted across:\n\n", nrow(master_table)))
cat("- **Model Types:** ZOIB (abundance), Gaussian (health/growth), JSDM Dirichlet (community)\n")
cat("- **Covariate Combinations:** FULL, NOHABITAT, NODEPTH, MINIMAL\n")
cat("- **CV Scenarios:** CV_02 (short-term), CV_30 (long-term), CV_ALL (complete)\n\n")

cat("## Model Selection Summary\n\n")
if (nrow(master_table) > 0) {
  cat("### Overall Statistics\n\n")
  cat(sprintf("- **Total Models:** %d\n", nrow(master_table)))
  cat(sprintf("- **Models Converged:** %d (%.1f%%)\n",
              sum(master_table$Converged, na.rm = TRUE),
              100 * mean(master_table$Converged == TRUE, na.rm = TRUE)))

  cat("\n### Best Models by Type\n\n")

  for (mt in unique(master_table$Model_Type)) {
    mt_results <- master_table[master_table$Model_Type == mt, ]
    if (nrow(mt_results) > 0) {
      best <- mt_results[which.min(mt_results$IC), ]
      cat(sprintf("#### %s\n", mt))
      cat(sprintf("- **Best Model:** %s\n", best$Model))
      cat(sprintf("- **CV:** %s\n", best$CV))
      cat(sprintf("- **IC:** %.1f (SE: %.1f)\n", best$IC, best$SE_IC))
      cat(sprintf("- **Converged:** %s\n\n", ifelse(best$Converged, "Yes", "No")))
    }
  }
}

cat("## Key Findings\n\n")

if (exists("hab_importance")) {
  cat("### Habitat (HAB) Importance\n\n")
  for (mt in unique(hab_importance$Model_Type)) {
    mt_hab <- hab_importance[hab_importance$Model_Type == mt, ]
    hab_yes <- mt_hab[mt_hab$Has_HAB == "Yes", ]
    hab_no <- mt_hab[mt_hab$Has_HAB == "No", ]
    if (nrow(hab_yes) > 0 && nrow(hab_no) > 0) {
      cat(sprintf("**%s:** ", mt))
      if (hab_yes$Mean_Rank < hab_no$Mean_Rank) {
        cat(sprintf("Models WITH HAB have better mean rank (%.1f vs %.1f)\n",
                    hab_yes$Mean_Rank, hab_no$Mean_Rank))
      } else {
        cat(sprintf("Models WITHOUT HAB have better mean rank (%.1f vs %.1f)\n",
                    hab_no$Mean_Rank, hab_yes$Mean_Rank))
      }
    }
  }
  cat("\n")
}

if (exists("depth_importance")) {
  cat("### Depth Importance\n\n")
  for (mt in unique(depth_importance$Model_Type)) {
    mt_depth <- depth_importance[depth_importance$Model_Type == mt, ]
    depth_yes <- mt_depth[mt_depth$Has_DEPTH == "Yes", ]
    depth_no <- mt_depth[mt_depth$Has_DEPTH == "No", ]
    if (nrow(depth_yes) > 0 && nrow(depth_no) > 0) {
      cat(sprintf("**%s:** ", mt))
      if (depth_yes$Mean_Rank < depth_no$Mean_Rank) {
        cat(sprintf("Models WITH DEPTH have better mean rank (%.1f vs %.1f)\n",
                    depth_yes$Mean_Rank, depth_no$Mean_Rank))
      } else {
        cat(sprintf("Models WITHOUT DEPTH have better mean rank (%.1f vs %.1f)\n",
                    depth_no$Mean_Rank, depth_yes$Mean_Rank))
      }
    }
  }
  cat("\n")
}

cat("## CV Scenario Comparison\n\n")
if (exists("cv_summary")) {
  for (mt in unique(cv_summary$Model_Type)) {
    mt_cv <- cv_summary[cv_summary$Model_Type == mt, ]
    best_cv <- mt_cv[which.min(mt_cv$Best_IC), ]
    cat(sprintf("- **%s:** Best CV scenario is **%s** (IC: %.1f)\n",
                mt, best_cv$CV, best_cv$Best_IC))
  }
}

cat("\n## Recommendations\n\n")
cat("Based on the model selection results:\n\n")
cat("1. **Winner Models:** Use the WINNER_*.rds files for downstream analysis and visualization\n")
cat("2. **Covariate Selection: ")
if (exists("combination_wins") && nrow(combination_wins) > 0) {
  top_combo <- combination_wins[which.max(combination_wins$Wins), ]
  cat(sprintf("The %s combination has the most wins (%d/%d)\n", top_combo$Model, top_combo$Wins, top_combo$Total))
} else {
  cat("\n")
}
cat("3. **CV Selection: Use the CV scenario with the best performance for each response type\n")
cat("4. **Visualization:** Use MASTER_Viz_Pipeline.R with the WINNER models\n\n")

cat("## Files Generated\n\n")
cat("### Master Tables\n")
cat("- `master_comparison_all_models.csv`: Combined results from all model types\n")
cat("- `consistency_analysis.rds`: Win rates and importance analysis\n\n")

cat("### Visualizations\n")
cat("- `fig_ic_heatmap.png`: Heatmap of IC values across all models\n")
cat("- `fig_winner_barplot.png`: Barplot of winning models\n")
cat("- `fig_ic_distribution.png`: Distribution of IC by combination\n")
cat("- `fig_convergence_summary.png`: Convergence rates\n")
cat("- `fig_cv_comparison.png`: Comparison across CV scenarios\n\n")

cat("## Next Steps\n\n")
cat("1. Review the WINNER models in each output directory\n")
cat("2. Run `05_MASTER_Viz_Pipeline_v5_FINAL.R` for publication-ready figures\n")
cat("3. Extract posterior summaries for manuscript\n\n")

cat("---\n")
cat("*Generated by 04_Synthesis_Full_LOO_Selection.R*\n")
cat(sprintf("*%s*\n", Sys.time()))

sink()

cat(sprintf("✓ Final report saved: %s\n", report_path))

# ============================================================================
# CREATE EXECUTIVE SUMMARY CSV
# ============================================================================

if (nrow(master_table) > 0) {
  exec_summary <- master_table %>%
    group_by(Model_Type) %>%
    filter(Rank == 1) %>%
    summarise(
      Best_Model = first(Model),
      Best_CV = first(CV),
      Best_IC = first(IC),
      SE_IC = first(SE_IC),
      Converged = first(Converged),
      .groups = "drop"
    )

  exec_path <- file.path(SYNTHESIS_DIR, "executive_summary.csv")
  write.csv(exec_summary, exec_path, row.names = FALSE)
  cat(sprintf("✓ Executive summary saved: %s\n", exec_path))
}

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n========== SYNTHESIS COMPLETE ==========\n")
cat(sprintf("\nTotal models analyzed: %d\n", nrow(master_table)))
cat(sprintf("Models converged: %d (%.1f%%)\n",
            sum(master_table$Converged, na.rm = TRUE),
            100 * mean(master_table$Converged == TRUE, na.rm = TRUE)))
cat(sprintf("Visualizations generated: 5\n"))
cat(sprintf("Reports generated: 2\n"))
cat(sprintf("\nAll outputs saved to: %s\n", SYNTHESIS_DIR))

cat("\n✓ Script 04 (Synthesis) completed!\n")
cat("\nNext: Run MASTER_Viz_Pipeline.R for publication-ready figures\n")
