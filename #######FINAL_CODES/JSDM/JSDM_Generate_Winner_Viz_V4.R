# ============================================================================ #
# SCRIPT: JSDM_Generate_Winner_Viz_V4.R
# PURPOSE: Publication-quality figures for JSDM winner model
# VERSION: 4.0 (Aligned with ZOIB visualization pipeline)
# ============================================================================ #

rm(list = ls())
gc()

# Required packages
libs <- c("brms", "dplyr", "ggplot2", "patchwork", "tidybayes",
          "bayesplot", "scales", "ggdist")
invisible(lapply(libs, library, character.only = TRUE))

# ============================================================================ #
# CONFIGURATION
# ============================================================================ #

cat("\n", rep("=", 70), "\n", sep = "")
cat("=== JSDM V4: VISUALIZATION PIPELINE ===\n")
cat(rep("=", 70), "\n\n", sep = "")

# Paths
# Paths
main_out_dir <- "C:/Users/rbfra/ONEDRIVE_NOVO/OneDrive/JSDM_Dirichlet_HypothesisTesting_V4"
winner_dir <- file.path(main_out_dir, "FINAL_WINNER_OUTPUTS")
viz_dir <- file.path(main_out_dir, "FIGURES")
dir.create(viz_dir, showWarnings = FALSE, recursive = TRUE)

# Load utility functions
source("C:/Users/rbfra/ONEDRIVE_NOVO/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_CODES/JSDM/JSDM_utils_V4.R")

# ============================================================================ #
# THEME: Publication-Quality (Aligned with ZOIB)
# ============================================================================ #

theme_publication <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(face = "bold", size = 14),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      plot.title = element_text(face = "bold", size = 16, hjust = 0),
      plot.subtitle = element_text(size = 11, color = "gray30"),
      plot.caption = element_text(size = 9, color = "gray50", hjust = 1),
      legend.position = "bottom",
      legend.text = element_text(size = 11),
      legend.title = element_text(face = "bold", size = 12),
      legend.background = element_rect(fill = "white", color = "gray50"),
      strip.background = element_rect(fill = "gray95", color = NA),
      strip.text = element_text(face = "bold", size = 12),
      panel.grid = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# Color palette (Okabe-Ito: color-blind friendly)
okabe_ito <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#000000",
  "#5694c7", "#c45d4f", "#8c8c8c", "#a65628",
  "#f781bf", "#98fb98", "#d9d9d9", "#a6cee3"
)

# ============================================================================ #
# FUNCTION: extract_fixed_effects
# PURPOSE: Safely extract fixed effects from fitted model
# ============================================================================ #

extract_fixed_effects <- function(brms_fit) {
  #' @param brms_fit Fitted brms model
  #' @return Data frame with posterior draws for fixed effects

  # Get summary
  model_summary <- summary(brms_fit)

  # Extract fixed effects
  fixed_effects <- model_summary$fixed

  if (is.null(fixed_effects) || nrow(fixed_effects) == 0) {
    stop("No fixed effects found in model")
  }

  # Get posterior draws for fixed effects
  posterior_draws <- as_draws_df(brms_fit, variable = "^b_", regex = TRUE)

  return(list(
    summary = fixed_effects,
    draws = posterior_draws
  ))
}

# ============================================================================ #
# FUNCTION: plot_fixed_effects_forest
# PURPOSE: Create forest plot of fixed effects (half-eye)
# ============================================================================ #

plot_fixed_effects_forest <- function(brms_fit, species_names,
                                       viz_dir, dataset_name) {
  #' @param brms_fit Fitted brms model
  #' @param species_names Character vector of species names
  #' @param viz_dir Directory to save figures
  #' @param dataset_name Name of CV dataset

  cat("\n  Generating fixed effects forest plot...\n")

  # Extract draws
  draws_data <- extract_fixed_effects(brms_fit)
  summary_data <- draws_data$summary

  # Remove intercept for cleaner visualization
  summary_plot <- summary_data[!grepl("^Intercept", rownames(summary_data)), , drop = FALSE]
  summary_plot$Parameter <- rownames(summary_plot)

  # Parse parameter names (format: b_species_predictor or b_predictor)
  summary_plot <- summary_plot %>%
    mutate(
      param_clean = gsub("^b_", "", Parameter),
      # Handle species-specific parameters (Dirichlet)
      species = ifelse(grepl("_", param_clean),
                       sapply(strsplit(param_clean, "_"), `[`, 1),
                       "Overall"),
      predictor = ifelse(grepl("_", param_clean),
                         sapply(strsplit(param_clean, "_"), `[`, 2),
                         param_clean)
    )

  # Filter out reference levels (usually ending in 1)
  summary_plot <- summary_plot %>%
    filter(!predictor %in% c("1", "Intercept"))

  # Assign colors by species
  n_species <- length(unique(summary_plot$species))
  species_colors <- setNames(okabe_ito[1:n_species], unique(summary_plot$species))

  # Create forest plot
  p <- ggplot(summary_plot,
              aes(x = Estimate, y = reorder(Parameter, Estimate))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = `l-95% CI`, xmax = `u+95% CI`),
                   height = 0.2, color = "gray50") +
    geom_point(aes(color = species), size = 3) +
    scale_color_manual(values = species_colors, name = "Species") +
    labs(
      title = sprintf("Fixed Effects: %s", dataset_name),
      subtitle = "Posterior means with 95% credible intervals",
      x = "Effect size (log-odds)",
      y = "Parameter"
    ) +
    theme_publication() +
    theme(
      legend.position = "right",
      axis.text.y = element_text(size = 10)
    )

  # Save
  output_file <- file.path(viz_dir, sprintf("forest_fixed_effects_%s.png", dataset_name))
  ggsave(output_file, p, width = 12, height = 10, dpi = 300)

  cat(sprintf("    -> Saved: %s\n", output_file))

  return(p)
}

# ============================================================================ #
# FUNCTION: plot_conditional_effects
# PURPOSE: Plot conditional effects for continuous predictors
# ============================================================================ #

plot_conditional_effects <- function(brms_fit, predictor, species_names,
                                      viz_dir, dataset_name) {
  #' @param brms_fit Fitted brms model
  #' @param predictor Name of predictor variable
  #' @param species_names Character vector of species names
  #' @param viz_dir Directory to save figures
  #' @param dataset_name Name of CV dataset

  cat(sprintf("  Generating conditional effects for %s...\n", predictor))

  # Get conditional effects
  cond_eff <- conditional_effects(brms_fit, effects = predictor)

  # Plot
  p <- plot(cond_eff)[[1]] +
    labs(
      title = sprintf("Conditional Effect: %s", predictor),
      subtitle = sprintf("Dataset: %s", dataset_name)
    ) +
    theme_publication()

  # Save
  safe_predictor <- gsub("[^A-Za-z0-9]", "_", predictor)
  output_file <- file.path(viz_dir, sprintf("conditional_%s_%s.png",
                                             safe_predictor, dataset_name))
  ggsave(output_file, p, width = 10, height = 8, dpi = 300)

  cat(sprintf("    -> Saved: %s\n", output_file))

  return(p)
}

# ============================================================================ #
# FUNCTION: plot_posterior_predictive_check
# PURPOSE: Create posterior predictive check plots
# ============================================================================ #

plot_posterior_predictive_check <- function(brms_fit, species_names,
                                             viz_dir, dataset_name) {
  #' @param brms_fit Fitted brms model
  #' @param species_names Character vector of species names
  #' @param viz_dir Directory to save figures
  #' @param dataset_name Name of CV dataset

  cat("  Generating posterior predictive checks...\n")

  # PPC plots
  y_rep <- posterior_predict(brms_fit, draws = 100)

  # Density overlay for first few species
  n_species_plot <- min(4, length(species_names))

  for (i in 1:n_species_plot) {
    species <- species_names[i]

    # Get observed data
    y_obs <- brms_fit$data[[species]]

    # Create PPC density overlay
    p <- ppc_dens_overlay(y_obs, y_rep[,, i]) +
      labs(
        title = sprintf("Posterior Predictive Check: %s", species),
        subtitle = sprintf("Dataset: %s", dataset_name)
      ) +
      theme_publication()

    safe_species <- gsub("[^A-Za-z0-9]", "_", species)
    output_file <- file.path(viz_dir, sprintf("ppc_density_%s_%s.png",
                                               safe_species, dataset_name))
    ggsave(output_file, p, width = 8, height = 6, dpi = 300)
  }

  cat("    -> Saved PPC plots\n")

  return(invisible(NULL))
}

# ============================================================================ #
# FUNCTION: plot_species_responses
# PURPOSE: Plot species-specific responses to environmental gradients
# ============================================================================ #

plot_species_responses <- function(brms_fit, species_names, predictors,
                                    viz_dir, dataset_name) {
  #' @param brms_fit Fitted brms model
  #' @param species_names Character vector of species names
  #' @param predictors Character vector of predictor names
  #' @param viz_dir Directory to save figures
  #' @param dataset_name Name of CV dataset

  cat("  Generating species response plots...\n")

  # Get posterior draws for mu (linear predictor)
  draws <- as_draws_df(brms_fit)

  # Find mu coefficients
  mu_cols <- grep("^mu_", names(draws), value = TRUE)

  if (length(mu_cols) == 0) {
    cat("    -> No mu parameters found, skipping species response plots\n")
    return(invisible(NULL))
  }

  # Parse coefficients to extract species-predictor combinations
  coef_info <- data.frame(
    param = mu_cols,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      species = gsub("^mu_(.+)_.*$", "\\1", param),
      predictor = gsub("^mu_.+_(.+)$", "\\1", param)
    )

  # Filter for selected predictors
  coef_info <- coef_info %>%
    filter(predictor %in% predictors)

  if (nrow(coef_info) == 0) {
    cat("    -> No matching predictor-species combinations found\n")
    return(invisible(NULL))
  }

  # Create plot for each predictor
  for (pred in predictors) {
    pred_data <- coef_info %>%
      filter(predictor == pred)

    if (nrow(pred_data) == 0) next

    # Extract posterior samples for this predictor
    plot_data <- data.frame(draws[, pred_data$param, drop = FALSE])
    names(plot_data) <- pred_data$species

    # Convert to long format
    plot_long <- plot_data %>%
      mutate(draw = 1:n()) %>%
      pivot_longer(cols = -draw, names_to = "species", values_to = "value")

    # Create ridge plot or violin plot
    p <- ggplot(plot_long, aes(x = value, y = species, fill = species)) +
      geom_violin(scale = "width", alpha = 0.7) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
      scale_fill_manual(values = okabe_ito[1:length(unique(pred_data$species))]) +
      labs(
        title = sprintf("Species Responses: %s", pred),
        subtitle = sprintf("Dataset: %s", dataset_name),
        x = "Effect size (linear predictor)",
        y = "Species"
      ) +
      theme_publication() +
      theme(legend.position = "none")

    safe_pred <- gsub("[^A-Za-z0-9]", "_", pred)
    output_file <- file.path(viz_dir, sprintf("species_response_%s_%s.png",
                                               safe_pred, dataset_name))
    ggsave(output_file, p, width = 10, height = 8, dpi = 300)

    cat(sprintf("    -> Saved: %s\n", output_file))
  }

  return(invisible(NULL))
}

# ============================================================================ #
# FUNCTION: create_summary_figure
# PURPOSE: Create multi-panel summary figure
# ============================================================================ #

create_summary_figure <- function(all_plots, viz_dir, dataset_name) {
  #' @param all_plots List of ggplot objects
  #' @param viz_dir Directory to save figures
  #' @param dataset_name Name of CV dataset

  cat("  Creating summary composite figure...\n")

  # Combine plots using patchwork
  n_plots <- min(length(all_plots), 4)

  if (n_plots == 1) {
    composite <- all_plots[[1]] + plot_annotation(
      title = sprintf("JSDM Results Summary: %s", dataset_name)
    )
  } else if (n_plots == 2) {
    composite <- all_plots[[1]] | all_plots[[2]] +
      plot_annotation(title = sprintf("JSDM Results Summary: %s", dataset_name))
  } else if (n_plots == 3) {
    composite <- (all_plots[[1]] | all_plots[[2]]) / all_plots[[3]] +
      plot_annotation(title = sprintf("JSDM Results Summary: %s", dataset_name))
  } else {
    composite <- (all_plots[[1]] | all_plots[[2]]) /
                 (all_plots[[3]] | all_plots[[4]]) +
      plot_annotation(title = sprintf("JSDM Results Summary: %s", dataset_name))
  }

  composite <- composite + theme_publication()

  # Save
  output_file <- file.path(viz_dir, sprintf("summary_composite_%s.png", dataset_name))
  ggsave(output_file, composite, width = 16, height = 12, dpi = 300)

  cat(sprintf("    -> Saved: %s\n", output_file))

  return(composite)
}

# ============================================================================ #
# MAIN EXECUTION
# ============================================================================ #

cat("\n### LOADING WINNER MODELS ###\n")

# Load winners summary
winners_file <- file.path(winner_dir, "winners_by_dataset.csv")

if (!file.exists(winners_file)) {
  stop(sprintf("Winners file not found: %s\nRun the main engine first.", winners_file))
}

winners_df <- read.csv(winners_file)
cat("\nWinner models:\n")
print(winners_df)

# Process each dataset's winner
all_viz_results <- list()

for (i in 1:nrow(winners_df)) {
  dataset_name <- winners_df$dataset[i]
  winner_model <- winners_df$winner[i]

  cat(sprintf("\n### PROCESSING: %s (Winner: %s) ###\n", dataset_name, winner_model))

  # Construct model file path
  model_file <- file.path(main_out_dir, "brms_cache",
                          sprintf("JSDM_%s_%s.rds", dataset_name, winner_model))

  if (!file.exists(model_file)) {
    warning(sprintf("Model file not found: %s\nSkipping this dataset.", model_file))
    next
  }

  cat(sprintf("Loading model from: %s\n", model_file))

  # Load model
  winner_fit <- readRDS(model_file)

  # Get species names from model data
  response_cols <- attr(winner_fit$data, "response")
  if (is.null(response_cols)) {
    # Try to extract from formula
    formula_str <- deparse(winner_fit$formula$formula)
    # Parse cbind() to get species names
    species_match <- regmatches(formula_str, regexpr("cbind\\(([^)]+)\\)", formula_str))
    if (length(species_match) > 0) {
      species_str <- gsub("cbind\\(([^)]+)\\)", "\\1", species_match)
      species_names <- trimws(strsplit(species_str, ",")[[1]])
    } else {
      species_names <- colnames(winner_fit$data)[grep("^mussismilia|^others|^porites|^agaricia",
                                                       colnames(winner_fit$data), ignore.case = TRUE)]
    }
  } else {
    species_names <- response_cols
  }

  cat(sprintf("Species in model: %s\n", paste(species_names, collapse = ", ")))

  # Create dataset-specific viz directory
  dataset_viz_dir <- file.path(viz_dir, dataset_name)
  dir.create(dataset_viz_dir, showWarnings = FALSE, recursive = TRUE)

  # --- Generate Plots ---
  all_plots <- list()

  # 1. Forest plot of fixed effects
  tryCatch({
    p1 <- plot_fixed_effects_forest(winner_fit, species_names, dataset_viz_dir, dataset_name)
    all_plots[["forest"]] <- p1
  }, error = function(e) {
    cat(sprintf("  ERROR creating forest plot: %s\n", e$message))
  })

  # 2. Conditional effects for continuous predictors
  predictors <- c("PC1_MAGNITUDE", "PC2_MAGNITUDE", "PC1_VARIABILITY", "PC2_VARIABILITY")

  for (pred in predictors) {
    tryCatch({
      p <- plot_conditional_effects(winner_fit, pred, species_names,
                                     dataset_viz_dir, dataset_name)
      all_plots[[paste0("conditional_", pred)]] <- p
    }, error = function(e) {
      cat(sprintf("  ERROR creating conditional plot for %s: %s\n", pred, e$message))
    })
  }

  # 3. Posterior predictive checks
  tryCatch({
    plot_posterior_predictive_check(winner_fit, species_names, dataset_viz_dir, dataset_name)
  }, error = function(e) {
    cat(sprintf("  ERROR creating PPC: %s\n", e$message))
  })

  # 4. Species-specific responses
  tryCatch({
    plot_species_responses(winner_fit, species_names, predictors, dataset_viz_dir, dataset_name)
  }, error = function(e) {
    cat(sprintf("  ERROR creating species response plots: %s\n", e$message))
  })

  # 5. Summary composite
  tryCatch({
    create_summary_figure(all_plots, dataset_viz_dir, dataset_name)
  }, error = function(e) {
    cat(sprintf("  ERROR creating summary figure: %s\n", e$message))
  })

  all_viz_results[[dataset_name]] <- list(
    winner = winner_model,
    plots = all_plots,
    viz_dir = dataset_viz_dir
  )

  cat(sprintf("\n=== Completed visualization for %s ===\n", dataset_name))
}

# ============================================================================ #
# FINAL REPORT
# ============================================================================ #

cat("\n", rep("=", 70), "\n", sep = "")
cat("=== VISUALIZATION COMPLETE ===\n")
cat(rep("=", 70), "\n\n", sep = "")

cat(sprintf("All figures saved to: %s\n", viz_dir))

for (ds in names(all_viz_results)) {
  cat(sprintf("\n  %s:\n", ds))
  cat(sprintf("    - Winner model: %s\n", all_viz_results[[ds]]$winner))
  cat(sprintf("    - Plots generated: %d\n", length(all_viz_results[[ds]]$plots)))
  cat(sprintf("    - Location: %s\n", all_viz_results[[ds]]$viz_dir))
}

cat("\nNext steps:\n")
cat("  1. Review generated figures\n")
cat("  2. Select publication-ready figures\n")
cat("  3. Integrate into manuscript\n")

cat("\n")
