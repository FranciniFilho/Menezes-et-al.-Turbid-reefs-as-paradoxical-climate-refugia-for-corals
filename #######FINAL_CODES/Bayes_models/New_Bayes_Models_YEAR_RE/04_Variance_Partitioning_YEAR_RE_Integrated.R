# ============================================================================
# 04_Variance_Partitioning_YEAR_RE_Integrated.R
# ============================================================================
# Variance Partitioning for YEAR_RE Models - HYBRID APPROACH
#
# Unlike the legacy VP pipeline (which fits models with combined CV predictors),
# this script analyzes models that were ALREADY FITTED separately for each CV scale.
#
# Key features:
#   1. Loads winner models from YEAR_RE output directories
#   2. Extracts R² and YEAR random effect variance for each CV scale
#   3. Compares variance explained by CV_02, CV_30, CV_ALL
#   4. Generates publication-ready plots
#
# Date: 2026-02-09
# ============================================================================

# ============================================================================
# 1. CONFIGURATION AND LIBRARIES
# ============================================================================

suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
  library(readr)
  library(scales)
  library(patchwork)
})

# Configure options
options(mc.cores = parallel::detectCores() - 2)
set.seed(42)

# ============================================================================
# 2. PATH DEFINITIONS
# ============================================================================

BASE_OUTPUT_DIR_YEAR_RE <- "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/"
BASE_OUTPUT_DIR_YEAR_RE_PRIOR_SENS <- "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output_PRIOR_SENSITIVITY_FULLGRID_v1/"

normalize_prior_scenario <- function(prior_tag) {
  tag <- tolower(trimws(as.character(prior_tag)))
  if (tag %in% c("weaklyinformative", "weakly_informative", "wi")) return("WeaklyInformative")
  if (tag %in% c("informative", "inf")) return("Informative")
  stop(sprintf("Unknown prior scenario: %s", prior_tag))
}

resolve_vp_namespace <- function(run_namespace = c("canonical", "prior_sens_fullgrid"),
                                 prior_scenario_target = "WeaklyInformative") {
  run_namespace <- match.arg(run_namespace)
  prior_scenario_target <- normalize_prior_scenario(prior_scenario_target)

  if (run_namespace == "canonical") {
    return(list(
      run_namespace = run_namespace,
      prior_scenario_target = prior_scenario_target,
      scenario_tag_upper = toupper(gsub("[^A-Za-z0-9]", "", prior_scenario_target)),
      base_dir = BASE_OUTPUT_DIR_YEAR_RE,
      vp_output_dir = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Variance_Partitioning_YEAR_RE/"),
      model_dirs = list(
        ZOIB_Abundance = file.path(BASE_OUTPUT_DIR_YEAR_RE, "ZOIB_Abundance_YEAR_RE/"),
        RGR = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_RGR/RGR/"),
        Health_PC1 = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_Health_YEAR_RE/HEALTH_PC1/"),
        Health_PC2 = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_Health_YEAR_RE/HEALTH_PC2/"),
        JSDM = file.path(BASE_OUTPUT_DIR_YEAR_RE, "JSDM_Dirichlet_YEAR_RE/")
      )
    ))
  }

  list(
    run_namespace = run_namespace,
    prior_scenario_target = prior_scenario_target,
    scenario_tag_upper = toupper(gsub("[^A-Za-z0-9]", "", prior_scenario_target)),
    base_dir = BASE_OUTPUT_DIR_YEAR_RE_PRIOR_SENS,
    vp_output_dir = file.path(BASE_OUTPUT_DIR_YEAR_RE_PRIOR_SENS, "PS_FULLGRID_Variance_Partitioning_YEAR_RE/"),
    model_dirs = list(
      ZOIB_Abundance = file.path(BASE_OUTPUT_DIR_YEAR_RE_PRIOR_SENS, "PS_FULLGRID_ZOIB_YEAR_RE/"),
      RGR = file.path(BASE_OUTPUT_DIR_YEAR_RE_PRIOR_SENS, "PS_FULLGRID_GAUSSIAN_RGR_YEAR_RE/RGR/"),
      Health_PC1 = file.path(BASE_OUTPUT_DIR_YEAR_RE_PRIOR_SENS, "PS_FULLGRID_GAUSSIAN_HEALTH_YEAR_RE/HEALTH_PC1/"),
      Health_PC2 = file.path(BASE_OUTPUT_DIR_YEAR_RE_PRIOR_SENS, "PS_FULLGRID_GAUSSIAN_HEALTH_YEAR_RE/HEALTH_PC2/"),
      JSDM = file.path(BASE_OUTPUT_DIR_YEAR_RE_PRIOR_SENS, "PS_FULLGRID_JSDM_DIRICHLET_YEAR_RE/")
    )
  )
}

RUN_NAMESPACE <- Sys.getenv("RUN_NAMESPACE", "canonical")
PRIOR_SCENARIO_TARGET <- Sys.getenv("PRIOR_SCENARIO_TARGET", "WeaklyInformative")
VP_NAMESPACE_CFG <- resolve_vp_namespace(RUN_NAMESPACE, PRIOR_SCENARIO_TARGET)

VP_OUTPUT_DIR <- VP_NAMESPACE_CFG$vp_output_dir
dir.create(VP_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
model_dirs <- VP_NAMESPACE_CFG$model_dirs

# CV scales to compare
cv_scales <- c("CV_02", "CV_30", "CV_ALL")

# ============================================================================
# 3. DATA LOADING FUNCTIONS
# ============================================================================

#' Find winner model files for a given model type and CV
#' @param model_dir Directory containing model files
#' @param cv_label CV label (e.g., "CV_02")
#' @param model_type Type identifier (e.g., "zoib_abundance_year_re")
#' @return Path to winner model file or NULL
find_winner_model <- function(model_dir, cv_label, model_type,
                              run_namespace = VP_NAMESPACE_CFG$run_namespace,
                              prior_scenario_target = VP_NAMESPACE_CFG$prior_scenario_target) {

  # Look for WINNER files with the CV label
  all_files <- list.files(model_dir, pattern = "\\.rds$", full.names = TRUE, recursive = TRUE)
  all_files <- all_files[grepl("WINNER", all_files, ignore.case = TRUE)]

  # Filter by CV label (e.g., "CV_02", "CV_30", "CV_ALL")
  cv_pattern <- gsub("_", "", cv_label)  # Remove underscores for matching
  cv_pattern <- gsub("0", "", cv_pattern) # Remove zeros for looser matching
  matching_files <- all_files[grepl(cv_label, all_files, ignore.case = TRUE)]

  if (run_namespace == "prior_sens_fullgrid") {
    matching_files <- matching_files[grepl(prior_scenario_target, matching_files, ignore.case = TRUE)]
  }

  if (length(matching_files) > 0) {
    matching_files <- sort(matching_files)
    return(matching_files[1])
  }

  return(NULL)
}

#' Load all winner models for a response type across CV scales
#' @param response_type "ZOIB_Abundance", "Health_PC1", "Health_PC2", or "JSDM"
#' @return List with models for each CV scale
load_winner_models_by_cv <- function(response_type,
                                     run_namespace = VP_NAMESPACE_CFG$run_namespace,
                                     prior_scenario_target = VP_NAMESPACE_CFG$prior_scenario_target) {

  cat(sprintf("\n=== Loading %s models ===\n", response_type))

  models <- list()

  # Define model_type identifier based on response_type
  model_type_map <- list(
    "ZOIB_Abundance" = "zoib_abundance_year_re",
    "RGR" = "rgr",
    "Health_PC1" = "gaussian_health_pc1_year_re",
    "Health_PC2" = "gaussian_health_pc2_year_re",
    "JSDM" = "jsdm_year_re"
  )

  model_type <- model_type_map[[response_type]]
  model_dir <- model_dirs[[response_type]]

  if (is.null(model_dir) || !dir.exists(model_dir)) {
    warning(sprintf("Directory not found for %s: %s", response_type, model_dir))
    return(NULL)
  }

  for (cv in cv_scales) {
    model_path <- find_winner_model(
      model_dir,
      cv,
      model_type,
      run_namespace = run_namespace,
      prior_scenario_target = prior_scenario_target
    )

    if (!is.null(model_path)) {
      cat(sprintf("  ✓ %s: %s\n", cv, basename(model_path)))
      models[[cv]] <- list(
        fit = readRDS(model_path),
        path = model_path,
        cv = cv
      )
    } else {
      cat(sprintf("  ⚠ %s: No winner model found\n", cv))
      models[[cv]] <- NULL
    }
  }

  # Count loaded models
  n_loaded <- sum(sapply(models, function(x) !is.null(x)))
  cat(sprintf("  Loaded %d of %d models\n", n_loaded, length(cv_scales)))

  return(models)
}

#' Load combined results CSV if available
#' @param response_type Response type
#' @return Dataframe with combined results
load_combined_results <- function(response_type,
                                  run_namespace = VP_NAMESPACE_CFG$run_namespace,
                                  prior_scenario_target = VP_NAMESPACE_CFG$prior_scenario_target) {

  # Map response types to result files
  if (run_namespace == "prior_sens_fullgrid") {
    result_files <- list(
      "ZOIB_Abundance" = file.path(VP_NAMESPACE_CFG$base_dir, "PS_FULLGRID_ZOIB_YEAR_RE", "ZOIB_YEAR_RE_combined_results_all_priors_PS_FULLGRID.csv"),
      "RGR" = file.path(VP_NAMESPACE_CFG$base_dir, "PS_FULLGRID_GAUSSIAN_RGR_YEAR_RE", "RGR_combined_results_all_priors_PS_FULLGRID.csv"),
      "Health_PC1" = file.path(VP_NAMESPACE_CFG$base_dir, "PS_FULLGRID_GAUSSIAN_HEALTH_YEAR_RE", "Health_YEAR_RE_combined_results_all_priors_PS_FULLGRID.csv"),
      "Health_PC2" = file.path(VP_NAMESPACE_CFG$base_dir, "PS_FULLGRID_GAUSSIAN_HEALTH_YEAR_RE", "Health_YEAR_RE_combined_results_all_priors_PS_FULLGRID.csv"),
      "JSDM" = file.path(VP_NAMESPACE_CFG$base_dir, "PS_FULLGRID_JSDM_DIRICHLET_YEAR_RE", "JSDM_YEAR_RE_combined_results_all_priors_PS_FULLGRID.csv")
    )
  } else {
    result_files <- list(
      "ZOIB_Abundance" = file.path(model_dirs$ZOIB_Abundance, "ZOIB_Abundance_YEAR_RE_combined_results.csv"),
      "RGR" = file.path(VP_NAMESPACE_CFG$base_dir, "Gaussian_RGR", "RGR_combined_results.csv"),
      "Health_PC1" = file.path(VP_NAMESPACE_CFG$base_dir, "Gaussian_Health_YEAR_RE", "Health_YEAR_RE_combined_results.csv"),
      "Health_PC2" = file.path(VP_NAMESPACE_CFG$base_dir, "Gaussian_Health_YEAR_RE", "Health_YEAR_RE_combined_results.csv"),
      "JSDM" = file.path(model_dirs$JSDM, "JSDM_YEAR_RE_combined_results.csv")
    )
  }

  result_file <- result_files[[response_type]]

  if (!is.null(result_file) && file.exists(result_file)) {
    df <- read.csv(result_file, stringsAsFactors = FALSE)
    if (run_namespace == "prior_sens_fullgrid" && ("Prior_Scenario" %in% names(df))) {
      df <- df[df$Prior_Scenario == prior_scenario_target, , drop = FALSE]
    }
    return(df)
  }

  return(NULL)
}

# ============================================================================
# 4. VARIANCE EXTRACTION FUNCTIONS
# ============================================================================

#' Extract R² from a brmsfit object
#' @param fit brmsfit object
#' @param resp Response component (for ZOIB: "mu", NULL for Gaussian)
#' @return Vector of posterior R² samples
extract_r2 <- function(fit, resp = NULL) {

  r2_samples <- tryCatch({
    if (!is.null(resp)) {
      as.vector(bayes_R2(fit, resp = resp, summary = FALSE))
    } else {
      as.vector(bayes_R2(fit, summary = FALSE))
    }
  }, error = function(e) {
    warning(sprintf("Could not extract R²: %s", e$message))
    return(NA)
  })

  return(r2_samples)
}

#' Extract YEAR random effect variance from a brmsfit object
#' @param fit brmsfit object
#' @return Vector of posterior YEAR SD samples
extract_year_variance <- function(fit) {

  # Get VarCorr
  vc <- VarCorr(fit)

  # Extract YEAR random effect
  if ("YEAR" %in% names(vc)) {
    year_samples <- as.vector(vc$YEAR$sd[1, ])  # SD samples
    year_variance <- year_samples^2  # Convert to variance
    return(year_variance)
  }

  warning("YEAR random effect not found in model")
  return(NULL)
}

#' Extract REEF random effect variance from a brmsfit object
#' @param fit brmsfit object
#' @return Vector of posterior REEF SD samples
extract_reef_variance <- function(fit) {

  vc <- VarCorr(fit)

  if ("REEF" %in% names(vc)) {
    reef_samples <- as.vector(vc$REEF$sd[1, ])
    reef_variance <- reef_samples^2
    return(reef_variance)
  }

  return(NULL)
}

#' Calculate total variance components for a model
#' @param fit brmsfit object
#' @param response_type "ZOIB_Abundance", "Health_PC1", "Health_PC2", or "JSDM"
#' @return List with R² and variance components
calculate_variance_components <- function(fit, response_type) {

  components <- list()

  # R² extraction
  if (response_type == "ZOIB_Abundance") {
    components$r2 <- extract_r2(fit, resp = "mu")
  } else if (response_type == "JSDM") {
    # For JSDM, use the first category's R² or overall
    components$r2 <- extract_r2(fit, resp = NULL)
  } else {
    components$r2 <- extract_r2(fit, resp = NULL)
  }

  # Random effect variances
  components$year_variance <- extract_year_variance(fit)
  components$reef_variance <- extract_reef_variance(fit)

  # Calculate proportion of variance due to YEAR
  if (!is.null(components$year_variance) && !is.null(components$r2)) {
    # Approximate: YEAR variance / Total variance
    # Total variance ≈ YEAR variance + REEF variance + residual
    # For Gaussian models, we can get residual variance from the fit

    sigma <- tryCatch({
      as.vector(fit$sigma$sigma)
    }, error = function(e) NULL)

    if (!is.null(sigma)) {
      total_variance <- components$year_variance + sigma
      components$year_proportion <- components$year_variance / total_variance
    }
  }

  return(components)
}

#' Compare variance components across CV scales
#' @param models List of models for each CV scale
#' @param response_type Response type
#' @return Dataframe with comparison
compare_cv_scales <- function(models, response_type) {

  cat(sprintf("\n=== Comparing CV scales for %s ===\n", response_type))

  results <- list()

  for (cv in names(models)) {
    if (is.null(models[[cv]])) next

    fit <- models[[cv]]$fit
    components <- calculate_variance_components(fit, response_type)

    results[[cv]] <- data.frame(
      CV = cv,
      R2_Mean = mean(components$r2, na.rm = TRUE),
      R2_SD = sd(components$r2, na.rm = TRUE),
      R2_Lower = quantile(components$r2, 0.025, na.rm = TRUE),
      R2_Upper = quantile(components$r2, 0.975, na.rm = TRUE),
      YEAR_Var_Mean = if (!is.null(components$year_variance)) mean(components$year_variance) else NA,
      YEAR_Var_SD = if (!is.null(components$year_variance)) sd(components$year_variance) else NA,
      REEF_Var_Mean = if (!is.null(components$reef_variance)) mean(components$reef_variance) else NA,
      REEF_Var_SD = if (!is.null(components$reef_variance)) sd(components$reef_variance) else NA,
      stringsAsFactors = FALSE
    )
  }

  if (length(results) == 0) {
    warning("No valid results for comparison")
    return(NULL)
  }

  combined <- do.call(rbind, results)

  # Rank by R²
  combined$R2_Rank <- rank(-combined$R2_Mean)

  cat(sprintf("\n  R² Ranking:\n"))
  for (i in 1:nrow(combined)) {
    cat(sprintf("    %d. %s: R² = %.3f [%.3f, %.3f]\n",
                combined$R2_Rank[i], combined$CV[i],
                combined$R2_Mean[i], combined$R2_Lower[i], combined$R2_Upper[i]))
  }

  return(combined)
}

# ============================================================================
# 5. PLOTTING FUNCTIONS
# ============================================================================

#' Publication-ready theme for VP plots
vp_plot_theme <- function() {

  theme_bw(base_size = 12) +
    theme(
      panel.grid.major = element_line(color = "gray85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.text = element_text(color = "black", size = 11),
      axis.title = element_text(size = 13, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10)
    )
}

#' Plot R² comparison across CV scales
#' @param comparison_data Dataframe from compare_cv_scales()
#' @param response_type Response type for title
#' @return ggplot object
plot_r2_comparison <- function(comparison_data, response_type) {

  if (is.null(comparison_data) || nrow(comparison_data) == 0) {
    warning("No data to plot")
    return(NULL)
  }

  # Check if all R2 values are NA (e.g., Dirichlet models)
  if (all(is.na(comparison_data$R2_Mean))) {
    # Return NULL and let the calling function handle this case
    return(NULL)
  }

  # Define all CV levels explicitly
  all_cv_levels <- c("CV_02", "CV_30", "CV_ALL")
  comparison_data$CV <- factor(comparison_data$CV, levels = all_cv_levels)

  # Filter out rows with NA R2 for plotting
  plot_data <- comparison_data[!is.na(comparison_data$R2_Mean), ]

  # Order by R² for display (only present levels)
  present_levels <- unique(plot_data$CV)
  plot_data$CV <- factor(plot_data$CV, levels = present_levels[order(plot_data$R2_Mean)])

  # Color palette - ensure all levels have colors
  cv_colors <- c("CV_02" = "#E69F00", "CV_30" = "#56B4E9", "CV_ALL" = "#009E73")

  p <- ggplot(plot_data, aes(x = CV, y = R2_Mean, fill = CV)) +
    geom_col(color = "black", linewidth = 0.3, width = 0.7) +
    geom_errorbar(aes(ymin = R2_Lower, ymax = R2_Upper),
                  linewidth = 0.8, width = 0.25) +
    scale_fill_manual(values = cv_colors, drop = FALSE) +
    labs(
      title = sprintf("Variance Explained by Temporal Scale\n%s", response_type),
      x = "Temporal Scale (Coefficient of Variation Window)",
      y = expression(R^2),
      fill = "CV Scale"
    ) +
    vp_plot_theme() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, hjust = 0.5))

  return(p)
}

#' Plot YEAR random effect variance comparison
#' @param comparison_data Dataframe from compare_cv_scales()
#' @param response_type Response type for title
#' @return ggplot object
plot_year_variance_comparison <- function(comparison_data, response_type) {

  if (is.null(comparison_data) || nrow(comparison_data) == 0) {
    return(NULL)
  }

  # Define all CV levels explicitly
  all_cv_levels <- c("CV_02", "CV_30", "CV_ALL")
  comparison_data$CV <- factor(comparison_data$CV, levels = all_cv_levels)

  # Filter out rows with NA YEAR variance for plotting
  plot_data <- comparison_data[!is.na(comparison_data$YEAR_Var_Mean), ]

  if (nrow(plot_data) == 0) {
    return(NULL)
  }

  # Order by YEAR variance for display (only present levels)
  present_levels <- unique(plot_data$CV)
  plot_data$CV <- factor(plot_data$CV, levels = present_levels[order(plot_data$YEAR_Var_Mean)])

  cv_colors <- c("CV_02" = "#E69F00", "CV_30" = "#56B4E9", "CV_ALL" = "#009E73")

  p <- ggplot(plot_data, aes(x = CV, y = YEAR_Var_Mean, fill = CV)) +
    geom_col(color = "black", linewidth = 0.3, width = 0.7) +
    scale_fill_manual(values = cv_colors, drop = FALSE) +
    labs(
      title = sprintf("Inter-annual Variance by Temporal Scale\n%s", response_type),
      x = "Temporal Scale (Coefficient of Variation Window)",
      y = expression(paste("YEAR Variance ", sigma[YEAR]^2)),
      fill = "CV Scale"
    ) +
    vp_plot_theme() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 0, hjust = 0.5))

  return(p)
}

#' Combined plot for a response type
#' @param models List of models for each CV scale
#' @param response_type Response type
#' @param output_dir Output directory
#' @return List with plots
create_combined_vp_plot <- function(models, response_type, output_dir) {

  comparison_data <- compare_cv_scales(models, response_type)

  if (is.null(comparison_data)) {
    warning(sprintf("Could not create plots for %s", response_type))
    return(NULL)
  }

  # Create individual plots
  p1 <- plot_r2_comparison(comparison_data, response_type)
  p2 <- plot_year_variance_comparison(comparison_data, response_type)

  # Handle case where R2 is not available (e.g., JSDM/Dirichlet models)
  if (is.null(p1) && !is.null(p2)) {
    # Add note that R2 is not available
    p2 <- p2 + labs(subtitle = "Note: R2 not available for this model type")
    combined <- p2
  } else if (!is.null(p1) && !is.null(p2)) {
    combined <- p1 / p2 +
      plot_layout(heights = c(1, 1)) +
      plot_annotation(
        title = sprintf("Variance Partitioning: %s", response_type),
        theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
      )
  } else if (!is.null(p1)) {
    combined <- p1
  } else {
    combined <- NULL
  }

  # Save
  if (!is.null(combined)) {
    # Create filename
    safe_name <- gsub(" ", "_", response_type)
    output_png <- file.path(output_dir, sprintf("VP_Comparison_%s.png", safe_name))
    output_pdf <- file.path(output_dir, sprintf("VP_Comparison_%s.pdf", safe_name))

    ggsave(output_png, combined, width = 10, height = ifelse(is.null(p1), 6, 8), dpi = 300)
    ggsave(output_pdf, combined, width = 10, height = ifelse(is.null(p1), 6, 8), dpi = 300)

    cat(sprintf("  ✓ Saved: %s\n", basename(output_png)))
  }

  # Save comparison data
  data_output <- file.path(output_dir, sprintf("VP_Data_%s.csv", safe_name))
  write.csv(comparison_data, data_output, row.names = FALSE)
  cat(sprintf("  ✓ Saved data: %s\n", basename(data_output)))

  return(list(
    comparison = combined,
    r2_plot = p1,
    variance_plot = p2,
    data = comparison_data
  ))
}

#' Master plot comparing all response types
#' @param all_results List of results for all response types
#' @param output_dir Output directory
#' @return Combined plot
create_master_vp_comparison <- function(all_results, output_dir) {

  cat("\n=== Creating master VP comparison ===\n")

  # Combine all data
  combined_data <- list()
  for (resp in names(all_results)) {
    if (!is.null(all_results[[resp]]$data)) {
      df <- all_results[[resp]]$data
      df$Response <- resp
      combined_data[[resp]] <- df
    }
  }

  if (length(combined_data) == 0) {
    warning("No data to combine for master plot")
    return(NULL)
  }

  master_df <- do.call(rbind, combined_data)

  # Filter out NaN R2 values (for JSDM/Dirichlet models)
  plot_df <- master_df[!is.na(master_df$R2_Mean), ]

  # Add annotation for responses without R2
  note_df <- master_df[is.na(master_df$R2_Mean) & !duplicated(master_df$Response), ]
  if (nrow(note_df) > 0) {
    for (i in 1:nrow(note_df)) {
      cat(sprintf("  Note: %s - R2 not available (Dirichlet model)\n", note_df$Response[i]))
    }
  }

  # Create multi-panel plot (only for responses with valid R2)
  p_master <- ggplot(plot_df, aes(x = CV, y = R2_Mean, fill = CV)) +
    geom_col(color = "black", linewidth = 0.3) +
    geom_errorbar(aes(ymin = R2_Lower, ymax = R2_Upper), linewidth = 0.6, width = 0.2) +
    facet_wrap(~Response, scales = "free_y") +
    scale_fill_manual(values = c("CV_02" = "#E69F00", "CV_30" = "#56B4E9", "CV_ALL" = "#009E73")) +
    labs(
      title = "Variance Explained by Temporal Scale Across All Responses",
      subtitle = "Note: R2 not available for Dirichlet (JSDM) models",
      x = "Coefficient of Variation Window",
      y = expression(R^2),
      fill = "CV Scale"
    ) +
    vp_plot_theme()

  # Save
  output_png <- file.path(output_dir, "VP_Master_Comparison_All_Responses.png")
  output_pdf <- file.path(output_dir, "VP_Master_Comparison_All_Responses.pdf")

  ggsave(output_png, p_master, width = 12, height = 10, dpi = 300)
  ggsave(output_pdf, p_master, width = 12, height = 10, dpi = 300)

  cat(sprintf("  ✓ Saved master plot: %s\n", basename(output_png)))

  # Save master data
  write.csv(master_df, file.path(output_dir, "VP_Master_Data.csv"), row.names = FALSE)

  return(p_master)
}

# ============================================================================
# 6. SUMMARY STATISTICS
# ============================================================================

#' Generate summary statistics for all response types
#' @param all_results List of results for all response types
#' @return Summary dataframe
generate_vp_summary <- function(all_results) {

  cat("\n=== GENERATING VP SUMMARY ===\n")

  summary_list <- list()

  for (resp in names(all_results)) {
    data <- all_results[[resp]]$data

    if (!is.null(data) && nrow(data) > 0) {
      # Check if R2 values are valid (not all NA)
      valid_r2 <- !is.na(data$R2_Mean)

      if (sum(valid_r2) == 0) {
        # No valid R2 values (e.g., Dirichlet models) - use YEAR variance instead
        best_cv_idx <- which.min(data$YEAR_Var_Mean)  # Lower variance = more stable
        best_cv <- data$CV[best_cv_idx]
        year_var <- data$YEAR_Var_Mean[best_cv_idx]

        summary_list[[resp]] <- data.frame(
          Response = resp,
          Best_CV = best_cv,
          Best_R2 = NA,
          R2_Range = NA,
          N_models = nrow(data),
          YEAR_Var_Mean = year_var,
          Note = "R2 not available - used YEAR variance",
          stringsAsFactors = FALSE
        )
      } else {
        # Find best CV among valid R2 values
        valid_data <- data[valid_r2, ]
        best_cv_idx <- which.max(valid_data$R2_Mean)
        best_cv <- valid_data$CV[best_cv_idx]
        best_r2 <- valid_data$R2_Mean[best_cv_idx]

        # Calculate R² range (using only valid values)
        r2_range <- max(valid_data$R2_Mean, na.rm = TRUE) - min(valid_data$R2_Mean, na.rm = TRUE)

        summary_list[[resp]] <- data.frame(
          Response = resp,
          Best_CV = best_cv,
          Best_R2 = best_r2,
          R2_Range = r2_range,
          N_models = nrow(valid_data),
          YEAR_Var_Mean = data$YEAR_Var_Mean[best_cv_idx],
          Note = NA,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(summary_list) > 0) {
    summary_df <- do.call(rbind, summary_list)
    rownames(summary_df) <- NULL

    cat("\nVARIANCE PARTITIONING SUMMARY:\n")
    print(summary_df)

    return(summary_df)
  }

  return(NULL)
}

#' Interpret VP results in plain text
#' @param summary_df Summary dataframe from generate_vp_summary()
#' @return Character vector with interpretations
interpret_vp_results <- function(summary_df) {

  if (is.null(summary_df)) {
    return(NULL)
  }

  interpretations <- character(nrow(summary_df))

  for (i in 1:nrow(summary_df)) {
    resp <- summary_df$Response[i]
    best_cv <- summary_df$Best_CV[i]
    r2_range <- summary_df$R2_Range[i]

    # Handle JSDM case where R2 is not available
    if (is.na(r2_range)) {
      year_var <- summary_df$YEAR_Var_Mean[i]
      interpretations[i] <- sprintf(
        "%s: R2 not available (Dirichlet model). YEAR variance: %.3f (lower = more stable). Best scale: %s",
        resp, year_var, best_cv
      )
    } else if (r2_range < 0.05) {
      interpretations[i] <- sprintf(
        "%s: Similar variance explained across all temporal scales (range < 5%%).",
        resp
      )
    } else if (best_cv == "CV_02") {
      interpretations[i] <- sprintf(
        "%s: Short-term variability (CV_02) explains most variance, suggesting rapid environmental responses.",
        resp
      )
    } else if (best_cv == "CV_30") {
      interpretations[i] <- sprintf(
        "%s: Medium-term variability (CV_30) explains most variance, suggesting intermediate temporal responses.",
        resp
      )
    } else if (best_cv == "CV_ALL") {
      interpretations[i] <- sprintf(
        "%s: Long-term variability (CV_ALL) explains most variance, suggesting cumulative or integrated effects.",
        resp
      )
    }
  }

  return(interpretations)
}

# ============================================================================
# 7. MAIN EXECUTION FUNCTION
# ============================================================================

#' Run complete variance partitioning analysis for YEAR_RE models
#' @param response_types Vector of response types to analyze
#' @param output_dir Output directory for results
#' @return List with all results
run_year_re_variance_partitioning <- function(
    response_types = c("ZOIB_Abundance", "RGR", "Health_PC1", "Health_PC2", "JSDM"),
    output_dir = VP_OUTPUT_DIR,
    run_namespace = VP_NAMESPACE_CFG$run_namespace,
    prior_scenario_target = VP_NAMESPACE_CFG$prior_scenario_target) {

  cat("\n")
  cat("============================================================================\n")
  cat("VARIANCE PARTITIONING: YEAR_RE MODELS (HYBRID APPROACH)\n")
  cat(sprintf("Namespace: %s\n", run_namespace))
  if (run_namespace == "prior_sens_fullgrid") {
    cat(sprintf("Prior scenario target: %s\n", prior_scenario_target))
  }
  cat("============================================================================\n")

  # Create output subdirectories
  dir.create(file.path(output_dir, "plots"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "data"), recursive = TRUE, showWarnings = FALSE)

  all_results <- list()

  for (resp_type in response_types) {
    cat(sprintf("\n>>> Processing: %s <<<\n", resp_type))

    # Load models
    models <- load_winner_models_by_cv(
      resp_type,
      run_namespace = run_namespace,
      prior_scenario_target = prior_scenario_target
    )

    if (is.null(models) || all(sapply(models, is.null))) {
      cat(sprintf("  ⚠ Skipping %s: no models found\n", resp_type))
      next
    }

    # Create plots
    plot_dir <- file.path(output_dir, "plots")
    result <- create_combined_vp_plot(models, resp_type, plot_dir)

    all_results[[resp_type]] <- result
  }

  # Create master comparison
  master_plot <- create_master_vp_comparison(all_results, file.path(output_dir, "plots"))

  # Generate summary
  summary_df <- generate_vp_summary(all_results)

  if (!is.null(summary_df)) {
    interpretations <- interpret_vp_results(summary_df)

    cat("\n=== INTERPRETATION ===\n")
    for (i in seq_along(interpretations)) {
      cat(sprintf("%s\n", interpretations[i]))
    }

    # Save summary
    write.csv(summary_df, file.path(output_dir, "VP_Summary.csv"), row.names = FALSE)
  }

  cat("\n")
  cat("============================================================================\n")
  cat("VARIANCE PARTITIONING COMPLETE\n")
  cat(sprintf("Results saved to: %s\n", output_dir))
  cat("============================================================================\n")

  return(list(
    results = all_results,
    summary = summary_df,
    interpretations = interpretations,
    master_plot = master_plot
  ))
}

# ============================================================================
# END OF 04_Variance_Partitioning_YEAR_RE_Integrated.R
# ============================================================================

# Execute the variance partitioning analysis
if (!interactive() || TRUE) {
  results <- run_year_re_variance_partitioning(
    response_types = c("ZOIB_Abundance", "RGR", "Health_PC1", "Health_PC2", "JSDM"),
    output_dir = VP_OUTPUT_DIR,
    run_namespace = VP_NAMESPACE_CFG$run_namespace,
    prior_scenario_target = VP_NAMESPACE_CFG$prior_scenario_target
  )
}
