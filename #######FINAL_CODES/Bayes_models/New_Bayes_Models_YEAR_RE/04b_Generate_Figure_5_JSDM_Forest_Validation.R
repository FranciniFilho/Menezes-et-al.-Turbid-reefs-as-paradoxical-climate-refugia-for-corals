### ========================================================================= ###
### === VALIDATION FIGURE 5: JSDM FOREST WITH ABBREVIATED TAXA LABELS     === ###
### ========================================================================= ###

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(brms)
  library(ggplot2)
  library(dplyr)
})

theme_publication_validation <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(face = "bold", size = 14),
      axis.line = element_line(color = "black", linewidth = 0.5),
      plot.title = element_text(face = "bold", size = 16, hjust = 0),
      plot.subtitle = element_text(size = 11, color = "gray30"),
      plot.caption = element_text(size = 9, color = "gray50", hjust = 1),
      legend.position = "right",
      legend.text = element_text(size = 11),
      legend.title = element_text(face = "bold", size = 12),
      panel.grid = element_blank()
    )
}

get_current_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  getwd()
}

quote_plotmath_text <- function(x) {
  sprintf("'%s'", gsub("'", "\\'", x, fixed = TRUE))
}

clean_predictor_name_validation <- function(x) {
  x <- gsub("MAGNITUDE", "Mag.", x, fixed = TRUE)
  x <- gsub("VARIABILITY", "Var.", x, fixed = TRUE)
  x <- gsub("HABMERGED", "Habitat", x, fixed = TRUE)
  x <- gsub("DEPTHM", "Depth", x, fixed = TRUE)
  x <- gsub("ARCH", "xARC", x, fixed = TRUE)
  x
}

main <- function() {
  script_dir <- get_current_script_dir()
  repo_root <- normalizePath(file.path(script_dir, "..", "..", ".."), mustWork = TRUE)

  model_path <- file.path(
    repo_root, "#######FINAL_RESULTS", "BAYES_MODELS_FINAL_RES",
    "PS_FULLGRID_Global_Winners_YEAR_RE", "WINNER_GLOBAL_jsdm_weaklyinformative.rds"
  )
  output_dir <- file.path(script_dir, "validation_outputs")

  if (!file.exists(model_path)) stop(sprintf("Model file not found: %s", model_path))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  message("Reading model: ", model_path)
  model <- readRDS(model_path)
  if (!inherits(model, "brmsfit")) stop("The supplied RDS is not a brmsfit object.")
  if (!grepl("dirichlet", family(model)[["family"]], ignore.case = TRUE)) {
    stop("The supplied model is not a Dirichlet JSDM.")
  }

  fixed_effects <- summary(model)$fixed
  plot_data <- data.frame(
    parameter = rownames(fixed_effects),
    estimate = as.numeric(fixed_effects[, "Estimate"]),
    l95 = as.numeric(fixed_effects[, "l-95% CI"]),
    u95 = as.numeric(fixed_effects[, "u-95% CI"]),
    stringsAsFactors = FALSE
  )

  plot_data <- plot_data[!grepl("Intercept$", plot_data$parameter), , drop = FALSE]
  plot_data$species <- gsub("^mu([A-Z]+)prop.*$", "\\1", plot_data$parameter)
  plot_data$predictor <- gsub("^mu[A-Z]+prop_(.*)$", "\\1", plot_data$parameter)
  plot_data <- plot_data[plot_data$predictor != plot_data$parameter, , drop = FALSE]

  all_predictors <- unique(plot_data$predictor)
  mussismilia_rows <- data.frame(
    parameter = paste0("muMUSSISMILIAprop_", all_predictors),
    estimate = 0,
    l95 = 0,
    u95 = 0,
    species = "MUSSISMILIA",
    predictor = all_predictors,
    stringsAsFactors = FALSE
  )
  plot_data <- bind_rows(plot_data, mussismilia_rows)

  species_levels <- c("MUSSISMILIA", "TURF", "CCA", "CYANO", "MACROALGAE")
  axis_codes <- c(MUSSISMILIA = "M. hispida", TURF = "TA", CCA = "CCA", CYANO = "CYANO", MACROALGAE = "MACRO")
  legend_expr <- c(
    MUSSISMILIA = "italic('M. hispida')~'(ref)'",
    TURF = "'Turf Algae (TA)'",
    CCA = "'Crustose Coralline Algae (CCA)'",
    CYANO = "'Cyanobacteria (CYANO)'",
    MACROALGAE = "'Macroalgae (MACRO)'"
  )
  species_colors <- c(MUSSISMILIA = "#D55E00", TURF = "#009E73", CCA = "#0072B2", CYANO = "#CC79A7", MACROALGAE = "#E69F00")

  plot_data$predictor_clean <- clean_predictor_name_validation(plot_data$predictor)
  plot_data$species <- factor(plot_data$species, levels = species_levels)
  plot_data$is_reference <- plot_data$species == "MUSSISMILIA"
  plot_data <- plot_data |>
    arrange(species, abs(estimate)) |>
    mutate(label_id = paste0("row_", seq_len(n())))

  plot_data$axis_label_expr <- ifelse(
    plot_data$is_reference,
    paste0("italic('M. hispida')~'(ref):'~", quote_plotmath_text(plot_data$predictor_clean)),
    paste0(quote_plotmath_text(axis_codes[as.character(plot_data$species)]), "~':'~", quote_plotmath_text(plot_data$predictor_clean))
  )
  axis_label_map <- setNames(plot_data$axis_label_expr, plot_data$label_id)

  p <- ggplot(plot_data, aes(x = estimate, y = factor(label_id, levels = label_id))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.7) +
    geom_segment(aes(x = l95, xend = u95, yend = factor(label_id, levels = label_id), color = species), linewidth = 1) +
    geom_point(aes(color = species, size = is_reference, shape = is_reference)) +
    scale_y_discrete(labels = function(x) parse(text = unname(axis_label_map[x]))) +
    scale_color_manual(
      values = species_colors,
      breaks = names(legend_expr),
      labels = function(x) parse(text = unname(legend_expr[x])),
      name = "Species",
      guide = guide_legend(override.aes = list(shape = c(18, 16, 16, 16, 16), size = c(5, 3.5, 3.5, 3.5, 3.5)))
    ) +
    scale_size_manual(values = c("TRUE" = 5, "FALSE" = 3.5), guide = "none") +
    scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16), guide = "none") +
    labs(
      title = "JSDM: Species-Specific Environmental Responses",
      subtitle = "Model: GLOBAL_jsdm_weaklyinformative",
      x = expression("Effect size (log-ratio relative to " * italic("M. hispida") * ")"),
      y = NULL,
      caption = "Points: Posterior median | Error bars: 95% CI | Reference at 0"
    ) +
    theme_publication_validation() +
    theme(axis.text.y = element_text(size = 9), plot.subtitle = element_text(size = 10))

  output_png <- file.path(output_dir, "FIGURE_5_JSDM_Forest_abbrev.png")
  output_pdf <- file.path(output_dir, "FIGURE_5_JSDM_Forest_abbrev.pdf")
  ggsave(output_png, p, width = 12, height = 8, dpi = 600, bg = "white")
  ggsave(output_pdf, p, width = 12, height = 8, bg = "white")

  message("Saved PNG: ", output_png)
  message("Saved PDF: ", output_pdf)
}

main()
