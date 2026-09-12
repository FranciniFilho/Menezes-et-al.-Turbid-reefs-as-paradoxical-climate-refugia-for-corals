### ==================================================================================== ###
### === FIGURE: SUMMARY OF COLONY-LEVEL BAYESIAN MODELS                              === ###
### ===                                                                              === ###
### === Panels:                                                                       === ###
### ===   A - HABMERGED marginal effects: HEALTH_PC1 & HEALTH_PC2                   === ###
### ===   B - PC1 Variability smooth PDP: HEALTH_PC1 only                           === ###
### ===   C - RGR fixed effects forest plot                                          === ###
### ===                                                                              === ###
### === Models: WINNER_GLOBAL_health_pc1.rds / health_pc2.rds / rgr.rds             === ###
### === Output: Summary_Colony_Models.png + .pdf  (300 Dpi, 10 × 8 in)             === ###
### === JSDM is EXCLUDED from this figure.                                           === ###
### ==================================================================================== ###

rm(list = ls())
gc()

# ---------------------------------------------------------------------------
# 1. LIBRARIES
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(brms)
  library(tidybayes)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ---------------------------------------------------------------------------
# 2. SHARED SETUP (exact copy from 04_ pipeline)
# ---------------------------------------------------------------------------

## Publication theme — base_size = 14 suits the wider 14-inch canvas
theme_publication <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      text            = element_text(family = "sans", color = "black"),
      axis.text       = element_text(color = "black", size = base_size),
      axis.title      = element_text(face = "bold", size = base_size + 2),
      axis.line       = element_line(color = "black", linewidth = 0.5),
      plot.title      = element_text(face = "bold", size = base_size + 4, hjust = 0),
      plot.subtitle   = element_text(size = base_size, color = "gray30"),
      plot.caption    = element_text(size = base_size - 2, color = "gray50", hjust = 1),
      legend.position = "bottom",
      legend.text     = element_text(size = base_size),
      legend.title    = element_text(face = "bold", size = base_size + 1),
      strip.background = element_rect(fill = "gray95", color = NA),
      strip.text      = element_text(face = "bold", size = base_size + 1),
      panel.grid      = element_blank()
    )
}
theme_set(theme_publication())

## Okabe-Ito colorblind-safe palette
color_primary   <- "#0072B2"
color_secondary <- "#E69F00"
hab_colors      <- c("PA" = "#009E73", "RR_TP" = "#D55E00")

## Directories
WINNER_DIR <- "#######FINAL_RESULTS/BAYES_MODELS_YEAR_RE/Global_Winners_YEAR_RE"
OUT_DIR    <- "#######FINAL_RESULTS/Bayesian_Figures_YEAR_RE_v7_CATEGORICAL/Summary_Colony_Models"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

NDRAWS <- 2000

cat("\n", rep("=", 70), "\n", sep = "")
cat("SUMMARY COLONY-LEVEL MODELS FIGURE\n")
cat("Output dir:", OUT_DIR, "\n")
cat(rep("=", 70), "\n\n", sep = "")

# ---------------------------------------------------------------------------
# 3. LOAD MODELS
# ---------------------------------------------------------------------------
cat("Loading models...\n")
model_h1  <- readRDS(file.path(WINNER_DIR, "WINNER_GLOBAL_health_pc1.rds"))
model_h2  <- readRDS(file.path(WINNER_DIR, "WINNER_GLOBAL_health_pc2.rds"))
model_rgr <- readRDS(file.path(WINNER_DIR, "WINNER_GLOBAL_rgr.rds"))
cat("  ✓ WINNER_GLOBAL_health_pc1 loaded\n")
cat("  ✓ WINNER_GLOBAL_health_pc2 loaded\n")
cat("  ✓ WINNER_GLOBAL_rgr       loaded\n\n")

# ---------------------------------------------------------------------------
# 4. HELPER: newdata at mean predictors (all continuous = 0, i.e. scaled mean)
# ---------------------------------------------------------------------------
make_newdata_means <- function(model, habmerged_level) {
  d         <- model$data
  resp_col  <- as.character(model$formula$resp)
  num_vars  <- names(d)[sapply(d, is.numeric) & names(d) != resp_col]
  nd        <- as.data.frame(setNames(as.list(rep(0, length(num_vars))), num_vars))
  hab_levels <- levels(d$HABMERGED)
  nd$HABMERGED <- factor(habmerged_level, levels = hab_levels)
  nd
}

# ---------------------------------------------------------------------------
# 5. PANEL A — HABMERGED marginal predictions: HEALTH_PC1 & HEALTH_PC2
# ---------------------------------------------------------------------------
cat("Building Panel A: HABMERGED marginal effects...\n")

compute_hab_predictions <- function(model, response_label) {
  d          <- model$data
  hab_levels <- levels(d$HABMERGED)

  rows <- lapply(hab_levels, function(h) {
    nd   <- make_newdata_means(model, h)
    eprd <- posterior_epred(model, newdata = nd, ndraws = NDRAWS, re.form = NA)
    vals <- as.numeric(eprd)
    data.frame(
      response_label = response_label,
      hab            = h,
      y              = median(vals),
      lower89        = quantile(vals, 0.055),
      upper89        = quantile(vals, 0.945),
      lower95        = quantile(vals, 0.025),
      upper95        = quantile(vals, 0.975),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

hab_h1 <- compute_hab_predictions(model_h1, "Health PC1\n(FULL, CV_30)")
hab_h2 <- compute_hab_predictions(model_h2, "Health PC2\n(FULL, CV_02)")
hab_df <- rbind(hab_h1, hab_h2)
hab_df$response_label <- factor(hab_df$response_label,
                                 levels = c("Health PC1\n(FULL, CV_30)",
                                            "Health PC2\n(FULL, CV_02)"))
hab_df$hab <- factor(hab_df$hab, levels = c("PA", "RR_TP"))

panel_a <- ggplot(hab_df, aes(x = hab, y = y, color = hab, fill = hab)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.7) +
  geom_linerange(aes(ymin = lower95, ymax = upper95), linewidth = 1.2, alpha = 0.4) +
  geom_linerange(aes(ymin = lower89, ymax = upper89), linewidth = 2.0, alpha = 0.7) +
  geom_point(size = 3.5, shape = 21, color = "white", stroke = 1.5) +
  facet_wrap(~ response_label, scales = "free_y", nrow = 1) +
  scale_color_manual(values = hab_colors, guide = "none") +
  scale_fill_manual(values = hab_colors, guide = "none") +
  scale_x_discrete(labels = c("PA" = "Reef Wall\n(PA)", "RR_TP" = "Rocky Reef\n+ Reef Top")) +
  labs(
    title   = "Habitat effect on coral health",
    x       = NULL,
    y       = "Predicted value (population-level)",
    caption = "Points: median  |  Bars: 89% & 95% CrI  |  Continuous predictors at mean"
  ) +
  theme_publication() +
  theme(axis.text.x = element_text(size = 12, lineheight = 1.2))

cat("  ✓ Panel A built\n")

# ---------------------------------------------------------------------------
# 6. PANEL B — PC1 Variability smooth PDP (HEALTH_PC1 only)
# ---------------------------------------------------------------------------
cat("Building Panel B: PC1 Variability smooth...\n")

# Detect the exact column name for PC1 variability
d_h1 <- model_h1$data
pc1var_col <- grep("^PC1.VARIABILITY$|^PC1_VARIABILITY$", names(d_h1),
                   ignore.case = TRUE, value = TRUE)
if (length(pc1var_col) == 0) {
  # Broader fallback search
  pc1var_col <- grep("PC1.*VAR|VAR.*PC1", names(d_h1),
                     ignore.case = TRUE, value = TRUE)
}
if (length(pc1var_col) == 0) {
  stop("Cannot locate PC1 Variability column in model_h1$data. ",
       "Available columns: ", paste(names(d_h1), collapse = ", "))
}
pc1var_col <- pc1var_col[1]
cat(sprintf("  PC1 Variability column detected: '%s'\n", pc1var_col))

## Build newdata: PC1VARIABILITY varies over full range; all others at 0
pc1var_range <- range(d_h1[[pc1var_col]], na.rm = TRUE)
pc1var_seq   <- seq(pc1var_range[1], pc1var_range[2], length.out = 200)

nd_smooth <- make_newdata_means(model_h1, habmerged_level = levels(d_h1$HABMERGED)[1])
nd_smooth <- nd_smooth[rep(1, 200), ]
nd_smooth[[pc1var_col]] <- pc1var_seq
rownames(nd_smooth) <- NULL

## posterior_epred → NDRAWS × 200 matrix
eprd_smooth <- posterior_epred(model_h1, newdata = nd_smooth, ndraws = NDRAWS, re.form = NA)

smooth_df <- data.frame(
  x       = pc1var_seq,
  estimate = apply(eprd_smooth, 2, median),
  lower89  = apply(eprd_smooth, 2, quantile, probs = 0.055),
  upper89  = apply(eprd_smooth, 2, quantile, probs = 0.945),
  lower95  = apply(eprd_smooth, 2, quantile, probs = 0.025),
  upper95  = apply(eprd_smooth, 2, quantile, probs = 0.975)
)

panel_b <- ggplot(smooth_df, aes(x = x, y = estimate)) +
  geom_ribbon(aes(ymin = lower95, ymax = upper95),
              fill = color_primary, alpha = 0.15) +
  geom_ribbon(aes(ymin = lower89, ymax = upper89),
              fill = color_primary, alpha = 0.25) +
  geom_line(color = color_primary, linewidth = 1.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.7) +
  labs(
    title   = "PC1 Variability smooth — Health PC1",
    x       = "PC1 Variability (scaled)",
    y       = "Predicted Health PC1",
    caption = "Habitat fixed at PA  |  Others at mean  |  Inner: 89% CI  |  Outer: 95% CI"
  ) +
  theme_publication()

cat("  ✓ Panel B built\n")

# ---------------------------------------------------------------------------
# 7. PANEL C — RGR forest plot (all fixed effects)
# ---------------------------------------------------------------------------
cat("Building Panel C: RGR forest plot...\n")

## Human-readable label lookup
label_map <- c(
  "b_Intercept"            = "Intercept",
  "b_HABMERGEDRR_TP"       = "Habitat: RR+TP vs PA",
  "b_PC1MAGNITUDE"         = "PC1 Magnitude",
  "b_PC2MAGNITUDE"         = "PC2 Magnitude",
  "b_PC1VARIABILITY"       = "PC1 Variability",
  "b_PC2VARIABILITY"       = "PC2 Variability",
  "b_PC1_MAGNITUDE"        = "PC1 Magnitude",
  "b_PC2_MAGNITUDE"        = "PC2 Magnitude",
  "b_PC1_VARIABILITY"      = "PC1 Variability",
  "b_PC2_VARIABILITY"      = "PC2 Variability",
  "b_DEPTHM"               = "Depth",
  "b_PC1INTERACAO"         = "PC1 Biotic (local)",
  "b_PC2INTERACAO"         = "PC2 Biotic (local)"
)

clean_param_label <- function(nm) {
  # Try exact lookup first
  if (nm %in% names(label_map)) return(label_map[[nm]])
  # Smooth basis columns: bs_sPC1VARIABILITY_1 → "s(PC1 Variability)"
  if (grepl("^bs_s", nm)) {
    inner <- sub("^bs_s", "", nm)
    inner <- sub("_[0-9]+$", "", inner)
    inner <- gsub("_", " ", inner)
    inner <- gsub("VARIABILITY", "Variability", inner)
    inner <- gsub("MAGNITUDE", "Magnitude", inner)
    return(sprintf("s(%s)", inner))
  }
  # Generic cleanup
  nm <- sub("^b_", "", nm)
  nm <- gsub("_", " ", nm)
  nm
}

## Extract posterior draws — parametric (b_) and smooth basis (bs_s)
draws_df <- as_draws_df(model_rgr)
param_cols <- grep("^b_|^bs_s", names(draws_df), value = TRUE)
param_cols <- setdiff(param_cols, "b_Intercept")  # exclude intercept from main panel

if (length(param_cols) == 0) {
  # Graceful fallback: show intercept only with a note
  warning("No non-intercept parametric draws found for RGR model. Showing all b_ params.")
  param_cols <- grep("^b_", names(draws_df), value = TRUE)
}

## Long format for plotting
rgr_long <- draws_df %>%
  dplyr::select(dplyr::all_of(param_cols)) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "param_raw", values_to = ".value") %>%
  dplyr::mutate(param_label = sapply(param_raw, clean_param_label))

## Compute medians for ordering
param_order <- rgr_long %>%
  dplyr::group_by(param_label) %>%
  dplyr::summarise(abs_med = abs(median(.value)), .groups = "drop") %>%
  dplyr::arrange(abs_med) %>%
  dplyr::pull(param_label)

rgr_long$param_label <- factor(rgr_long$param_label, levels = param_order)

panel_c <- ggplot(rgr_long, aes(y = param_label, x = .value)) +
  stat_halfeye(
    .width         = c(0.89, 0.95),
    alpha          = 0.8,
    point_interval = median_qi,
    fill           = "gray70",
    color          = "gray30",
    point_size     = 1.8
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray30", linewidth = 0.8) +
  labs(
    title    = "RGR \u2014 Fixed effects (no credible driver)",
    subtitle = "Bayesian R\u00B2 = 0.294 [0.175, 0.404]",
    x        = "Effect estimate",
    y        = NULL,
    caption  = "All 95% CrI include zero  |  Points: median  |  Bars: 89% & 95% CrI"
  ) +
  coord_cartesian(xlim = c(-1.5, 1.5)) +   # clip tails; densities remain correct
  theme_publication() +
  theme(axis.text.y = element_text(size = 12))

cat("  ✓ Panel C built\n")

# ---------------------------------------------------------------------------
# 8. ASSEMBLE & SAVE
# ---------------------------------------------------------------------------
cat("\nAssembling composite figure...\n")

## Top row: panel A (2 facets) gets 55% of width; panel B gets 45%
## At 14-in canvas Panel B gets ~6.4 in — enough for its full title
top_row <- wrap_plots(panel_a, panel_b, widths = c(1.2, 1))

## Landscape-ish layout: rows roughly equal in height
fig <- top_row / panel_c +
  plot_layout(heights = c(1.1, 1)) +
  plot_annotation(
    title    = "Colony-level Bayesian model results",
    subtitle = "HEALTH_PC1 (R\u00B2 = 0.53) \u00B7 HEALTH_PC2 (R\u00B2 = 0.25) \u00B7 RGR (R\u00B2 = 0.29)",
    tag_levels = "A",
    theme    = theme_publication(base_size = 15)
  )

out_png <- file.path(OUT_DIR, "Summary_Colony_Models.png")
out_pdf <- file.path(OUT_DIR, "Summary_Colony_Models.pdf")

## 14 × 11 in (wider than tall — good for landscape/journal double-column)
ggsave(out_png, fig, width = 14, height = 11, dpi = 300, bg = "white")
ggsave(out_pdf, fig, width = 14, height = 11, bg = "white")

cat(sprintf("\n\u2713 Figure saved to: %s\n", OUT_DIR))
cat(sprintf("  PNG: %s\n", basename(out_png)))
cat(sprintf("  PDF: %s\n", basename(out_pdf)))
cat(rep("=", 70), "\n", sep = "")
