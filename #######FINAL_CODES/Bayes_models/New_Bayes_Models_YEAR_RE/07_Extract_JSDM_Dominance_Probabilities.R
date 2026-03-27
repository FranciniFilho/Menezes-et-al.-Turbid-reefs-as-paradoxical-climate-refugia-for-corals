# ============================================================================
# 07_Extract_JSDM_Dominance_Probabilities.R
# ============================================================================
# Extract posterior predictions and calculate dominance probabilities
# from the JSDM Dirichlet winner model.
#
# BUGS FIXED vs. previous version:
#   BUG 1 (CRITICAL) - Scaling inversion:
#     model$data holds ALREADY-SCALED predictors (scale() was applied inside
#     prepare_jsdm_data_with_year()). Extracting scale_params from model$data
#     yields center≈0, scale≈1, so raw PCA scores passed through unchanged
#     while the model expected properly re-centred values. Result: all
#     predictions driven to one extreme (100 % macroalgae everywhere).
#     FIX: use fitted(model, summary = FALSE) WITHOUT newdata → the model's
#          own correctly-scaled data is used; no external scaling needed.
#
#   BUG 2 (CRITICAL) - Missing HABMERGED column:
#     pred_data was built with HAB, but the winner model formula contains
#     HABMERGED (PA stays; RR + TP merged). posterior_predict(newdata = .)
#     failed silently (or used wrong columns), producing garbage output.
#     FIX: same as above – fitted() without newdata bypasses the issue.
#
#   BUG 3 - Response-column index assumed, not verified:
#     resp_names was a hardcoded vector; its order might not match the order
#     of the 3rd dimension of the prediction array.
#     FIX: response indices are now derived from dimnames(fitted_vals)[[3]].
#
# OUTPUT DIRECTORY changed to C:/Users/rbfra/OneDrive/JSDM_Dominance_Figures/
# to avoid Windows MAX_PATH issues caused by deeply-nested long filenames.
# ============================================================================

setwd(paste0(
  "C:/Users/rbfra/OneDrive/########PUBLICACOES/",
  "############Menezes et al. Mus his distribution and abundance Abrolhos/",
  "######FINAL/#######FINAL_CODES/Bayes_models/New_Bayes_Models_YEAR_RE/"
))

suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(dplyr)
  library(tidyr)
  library(tidybayes)
  library(ggplot2)
})

options(mc.cores = 4)
set.seed(42)

# ============================================================================
# PATH CONFIGURATION
# ============================================================================

winner_paths <- c(
  "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/Global_Winners_YEAR_RE/WINNER_GLOBAL_jsdm.rds",
  "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output_PRIOR_SENSITIVITY_FULLGRID_v1/PS_FULLGRID_Global_Winners_YEAR_RE/WINNER_GLOBAL_jsdm_weaklyinformative.rds",
  "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/JSDM_Dirichlet_YEAR_RE/CV_ALL/jsdm_year_re_nohabitat_interaction_CV_ALL.rds"
)

data_path <- paste0(
  "C:/Users/rbfra/OneDrive/########PUBLICACOES/",
  "############Menezes et al. Mus his distribution and abundance Abrolhos/",
  "######FINAL/#######FINAL_RESULTS/#####output_local_PCA_CV_all_FINAL/",
  "dados_abundancia_integrados_long_format.csv"
)

# Short output path avoids Windows MAX_PATH corruption
output_dir <- "C:/Users/rbfra/OneDrive/JSDM_Dominance_Figures/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# LOAD MODEL
# ============================================================================

cat("\n=== Loading JSDM Winner Model ===\n")

model <- NULL
for (p in winner_paths) {
  if (file.exists(p)) {
    cat(sprintf("  Found model at: %s\n", p))
    model <- tryCatch(readRDS(p), error = function(e) {
      cat(sprintf("  Failed to load: %s\n", e$message)); NULL
    })
    if (!is.null(model)) {
      cat("  Model loaded successfully!\n")
      break
    }
  }
}

if (is.null(model)) stop("Could not load JSDM winner model from any known location.")

cat("\n=== Model Summary ===\n")
print(model$formula)
cat(sprintf("\nFamily: %s\n", model$family$family))

# ============================================================================
# LOAD METADATA
# Two purposes:
#   1. site_geo  – geographic coordinates for spatial figures (from raw CSV)
#   2. obs_meta  – per-observation SITE/HAB/YEAR/ARCH/REEF identifiers that
#                  brms dropped from model$data (it only keeps formula columns)
#
# WHY re-prepare with prepare_jsdm_data_with_year():
#   brms stores only variables referenced in the formula (REEF, YEAR,
#   HABMERGED, predictor columns, responses). SITE, HAB, ARCH, SAMPLE_UNIT
#   are not in the formula → they are absent from model$data. The function
#   is deterministic: same CSV → same row order as model$data, so row i in
#   obs_meta corresponds exactly to row i in model$data.
# ============================================================================

cat("\n=== Loading geographic metadata ===\n")

raw_meta <- tryCatch(
  read.csv2(data_path, sep = ";", header = TRUE, stringsAsFactors = FALSE),
  error = function(e) read.csv(data_path, stringsAsFactors = FALSE)
)
names(raw_meta) <- toupper(names(raw_meta))

# One row per SITE x HAB (geographic info repeats across organisms)
site_geo <- raw_meta %>%
  mutate(SITE = as.character(SITE), HAB = as.character(HAB)) %>%
  select(SITE, HAB, LATITUDE, LONGITUDE, REEF_NAME) %>%
  distinct()

cat(sprintf("  Geographic metadata: %d unique SITE x HAB combinations\n", nrow(site_geo)))

# ── Re-prepare observation metadata (matches model$data row order exactly) ──
cat("\n=== Re-preparing observation metadata (for SITE/HAB/ARCH labels) ===\n")
cat("  Sourcing prepare_jsdm_data_with_year() ... ")

obs_meta <- tryCatch({
  # source() loads the function definitions into the current R session;
  # setwd() at top of this script ensures the file is found.
  source("00_LOO_Selection_Functions_YEAR_RE.R")
  cat("OK\n")
  prepare_jsdm_data_with_year(data_path, "CV_ALL")
}, error = function(e) {
  cat(sprintf("\n  Warning: could not source functions (%s)\n  Falling back to manual preparation.\n", e$message))
  NULL
})

if (is.null(obs_meta)) {
  # Manual lightweight fallback: replicate core logic of prepare_jsdm_data_with_year
  organisms_jsdm <- c("MUSSISMILIA_HISPIDA", "MACROALGAE", "CCA", "TURF", "CYANO")
  obs_meta <- raw_meta %>%
    mutate(
      COBERTURA = as.numeric(gsub(",", ".", as.character(COBERTURA))),
      COVER_PROP = COBERTURA / 100,
      CATEGORY = dplyr::case_when(
        ORGANISMO == "MUSSISMILIA_HISPIDA"              ~ "MUSSISMILIA",
        grepl("TURF",  ORGANISMO, ignore.case = TRUE)  ~ "TURF",
        grepl("CCA",   ORGANISMO, ignore.case = TRUE)  ~ "CCA",
        grepl("CYANO", ORGANISMO, ignore.case = TRUE)  ~ "CYANO",
        grepl("MACRO", ORGANISMO, ignore.case = TRUE)  ~ "MACROALGAE",
        TRUE                                            ~ "OTHER"
      ),
      SAMPLE_UNIT = paste(SITE, HAB, YEAR, sep = "_")
    ) %>%
    filter(!is.na(COBERTURA), CATEGORY != "OTHER") %>%
    group_by(SAMPLE_UNIT, SITE, HAB, YEAR, ARCH,
             PC1_MAGNITUDE, PC2_MAGNITUDE, PC1_VARIABILITY, PC2_VARIABILITY, DEPTH_M) %>%
    summarise(
      MUSSISMILIA_raw = sum(COVER_PROP[CATEGORY == "MUSSISMILIA"], na.rm = TRUE),
      MACROALGAE_raw  = sum(COVER_PROP[CATEGORY == "MACROALGAE"],  na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      REEF      = substr(as.character(SITE), 1, 3),
      HABMERGED = ifelse(HAB %in% c("RR", "TP"), "RR_TP", HAB),
      YEAR      = as.character(YEAR)
    )
}

obs_meta <- obs_meta %>%
  mutate(
    obs_idx  = row_number(),
    SITE     = as.character(SITE),
    HAB      = as.character(HAB),
    YEAR     = as.character(YEAR),
    ARCH     = as.character(ARCH),
    REEF     = as.character(REEF)
  )

cat(sprintf("  obs_meta rows: %d (model$data rows will be checked after fitted())\n",
            nrow(obs_meta)))

# ============================================================================
# EXTRACT FITTED VALUES — IN-SAMPLE, NO newdata
# ============================================================================
# fitted(summary = FALSE) returns [S, N, K]:
#   S = posterior draws
#   N = number of observations in model$data
#   K = number of Dirichlet categories
#
# Using the model's own data guarantees:
#   * correct predictor scaling
#   * correct column names (HABMERGED, REEF, DEPTHM, …)
#   * no "allow_new_levels" complications
# ============================================================================

cat("\n=== Extracting fitted values (in-sample, summary = FALSE) ===\n")
cat("  This may take several minutes...\n")

fitted_vals <- tryCatch({
  fitted(model, summary = FALSE, ndraws = 2000)
}, error = function(e) {
  cat(sprintf("  fitted() with ndraws failed: %s\n  Trying without ndraws...\n", e$message))
  tryCatch(
    fitted(model, summary = FALSE),
    error = function(e2) {
      cat(sprintf("  fitted() failed entirely: %s\n", e2$message))
      NULL
    }
  )
})

if (is.null(fitted_vals)) stop("Could not extract fitted values from JSDM model.")

cat(sprintf("  Fitted array dimensions: %s\n", paste(dim(fitted_vals), collapse = " x ")))

# ============================================================================
# IDENTIFY RESPONSE COLUMNS BY NAME (dimnames), NOT by hardcoded index
# ============================================================================

resp_dimnames <- dimnames(fitted_vals)[[3]]

cat(sprintf("  Response categories in model: %s\n",
            if (is.null(resp_dimnames)) "unnamed (using positional heuristic)"
            else paste(resp_dimnames, collapse = ", ")))

if (!is.null(resp_dimnames)) {
  muss_idx  <- which(grepl("MUSSISMILIA", resp_dimnames, ignore.case = TRUE))
  macro_idx <- which(grepl("MACROALGAE",  resp_dimnames, ignore.case = TRUE))
  cca_idx   <- which(grepl("CCA",           resp_dimnames, ignore.case = TRUE))
  # Filter out false positives (e.g. if a category name contained "CCA" as substring)
  cca_idx   <- setdiff(cca_idx, macro_idx)
} else {
  # Fallback: derive from formula cbind() order
  # make_jsdm_formula_year_re orders: MUSSISMILIA, TURF, CCA, CYANO, MACROALGAE
  cat("  WARNING: dimnames not available – using formula-order fallback\n")
  cat("  Assumed order: MUSSISMILIA=1, TURF=2, CCA=3, CYANO=4, MACROALGAE=5\n")
  muss_idx  <- 1L
  macro_idx <- 5L
  cca_idx   <- 3L
}

if (length(muss_idx) == 0 || length(macro_idx) == 0) {
  stop(sprintf(
    "Could not locate MUSSISMILIA (%d) or MACROALGAE (%d) in fitted output. dimnames: %s",
    length(muss_idx), length(macro_idx),
    paste(resp_dimnames, collapse = ", ")
  ))
}

cat(sprintf("  Index mapping → MUSSISMILIA: %d, MACROALGAE: %d, CCA: %d\n",
            muss_idx, macro_idx,
            if (length(cca_idx) > 0) cca_idx else NA_integer_))

# ============================================================================
# CALCULATE DOMINANCE PROBABILITIES (vectorised over observations)
# ============================================================================

cat("\n=== Calculating dominance probabilities ===\n")

n_obs   <- dim(fitted_vals)[2]
n_draws <- dim(fitted_vals)[1]

cat(sprintf("  %d observations × %d posterior draws\n", n_obs, n_draws))

# Vectorised: extract matrices [draws × obs] for each category
muss_mat  <- fitted_vals[, , muss_idx,  drop = TRUE]  # [S, N]
macro_mat <- fitted_vals[, , macro_idx, drop = TRUE]
cca_mat   <- if (length(cca_idx) > 0) fitted_vals[, , cca_idx, drop = TRUE] else matrix(NA, n_draws, n_obs)

# Per-observation summaries
dominance_results <- data.frame(
  obs_idx          = seq_len(n_obs),
  P_macro_coral    = colMeans(macro_mat > muss_mat,  na.rm = TRUE),
  P_cca_coral      = colMeans(cca_mat   > muss_mat,  na.rm = TRUE),
  mean_MUSSISMILIA = colMeans(muss_mat,  na.rm = TRUE),
  mean_MACROALGAE  = colMeans(macro_mat, na.rm = TRUE),
  mean_CCA         = colMeans(cca_mat,   na.rm = TRUE),
  lwr_MUSSISMILIA  = apply(muss_mat,  2, quantile, 0.025, na.rm = TRUE),
  upr_MUSSISMILIA  = apply(muss_mat,  2, quantile, 0.975, na.rm = TRUE),
  lwr_MACROALGAE   = apply(macro_mat, 2, quantile, 0.025, na.rm = TRUE),
  upr_MACROALGAE   = apply(macro_mat, 2, quantile, 0.975, na.rm = TRUE)
)

# ============================================================================
# MARGINAL PREDICTIONS (re_formula = NA) — for general ARC/HAB claims
# ============================================================================
# Conditional predictions (above) include REEF random effects, which for
# MACROALGAE span ±3.2 on the log-ratio scale (= 25× multiplier between
# extreme reefs). This makes conditional P(Macro > Coral) extreme at
# individual reefs but unrepresentative of the general environmental gradient.
#
# Marginal predictions integrate over the RE distribution, giving the
# population-average response to the environmental predictors only.
# Use MARGINAL for general ARC/HAB/gradient claims in the manuscript.
# Use CONDITIONAL for reef-specific claims only.
# ============================================================================

cat("\n=== Extracting MARGINAL fitted values (re_formula = NA) ===\n")

marginal_vals <- tryCatch({
  fitted(model, summary = FALSE, re_formula = NA, ndraws = 2000)
}, error = function(e) {
  cat(sprintf("  Marginal fitted() failed: %s\n", e$message))
  NULL
})

if (!is.null(marginal_vals)) {
  cat(sprintf("  Marginal array dimensions: %s\n", paste(dim(marginal_vals), collapse = " x ")))

  # Use same indices as conditional (verified above)
  marg_muss_mat  <- marginal_vals[, , muss_idx,  drop = TRUE]
  marg_macro_mat <- marginal_vals[, , macro_idx, drop = TRUE]
  marg_cca_mat   <- if (length(cca_idx) > 0) marginal_vals[, , cca_idx, drop = TRUE] else matrix(NA, dim(marginal_vals)[1], n_obs)

  dominance_results$P_macro_coral_marginal <- colMeans(marg_macro_mat > marg_muss_mat, na.rm = TRUE)
  dominance_results$P_cca_coral_marginal   <- colMeans(marg_cca_mat   > marg_muss_mat, na.rm = TRUE)
  dominance_results$mean_MUSSISMILIA_marg  <- colMeans(marg_muss_mat,  na.rm = TRUE)
  dominance_results$mean_MACROALGAE_marg   <- colMeans(marg_macro_mat, na.rm = TRUE)
  dominance_results$mean_CCA_marg          <- colMeans(marg_cca_mat,   na.rm = TRUE)

  cat("  Marginal dominance computed successfully.\n")
} else {
  cat("  WARNING: Marginal extraction failed — only conditional results available.\n")
}

# ============================================================================
# MERGE: dominance probabilities ← obs_meta ← geographic metadata
# ============================================================================

cat("\n=== Merging with observation metadata ===\n")

# Sanity check: row counts must agree
if (nrow(obs_meta) != n_obs) {
  warning(sprintf(
    "obs_meta has %d rows but fitted_vals has %d observations. ",
    nrow(obs_meta), n_obs,
    "Metadata may be misaligned – check CV_ALL data path."
  ))
  # Trim or pad obs_meta to n_obs rows as best-effort
  if (nrow(obs_meta) > n_obs) obs_meta <- obs_meta[seq_len(n_obs), ]
}

# Make sure obs_idx is correct (1:n_obs)
obs_meta <- obs_meta %>% mutate(obs_idx = row_number())

# Join dominance results to obs_meta (row-aligned, obs_idx is the key)
dominance_full <- dominance_results %>%
  left_join(obs_meta %>%
              select(obs_idx, SITE, HAB, YEAR, ARCH, REEF,
                     any_of("HABMERGED")),
            by = "obs_idx")

cat(sprintf("  After obs_meta join: %d rows, cols: %s\n",
            nrow(dominance_full),
            paste(names(dominance_full), collapse = ", ")))

# Join geographic coordinates via SITE + HAB
dominance_full <- dominance_full %>%
  left_join(site_geo, by = c("SITE", "HAB"))

n_geo_matched <- sum(!is.na(dominance_full$LATITUDE))
cat(sprintf("  Geographic join: %d/%d rows matched lat/lon\n",
            n_geo_matched, nrow(dominance_full)))

# Status classification
dominance_full <- dominance_full %>%
  mutate(
    status = case_when(
      P_macro_coral >= 0.95 ~ "Macroalgae-dominated",
      P_macro_coral <= 0.05 ~ "Coral-favored",
      P_macro_coral >= 0.50 ~ "Macroalgae-advantaged",
      TRUE                  ~ "Uncertain"
    ),
    arch_status = paste(ARCH, status, sep = ": ")
  )

cat(sprintf("  Total observations in output: %d\n", nrow(dominance_full)))

# Quick sanity check
status_table <- table(dominance_full$status)
cat("\n  Status distribution:\n")
print(status_table)

if (all(dominance_full$status == "Macroalgae-dominated")) {
  warning(paste(
    "ALL observations classified as Macroalgae-dominated.",
    "Check that the model loaded correctly and that fitted() returned",
    "sensible values. mean_MUSSISMILIA range:",
    round(min(dominance_full$mean_MUSSISMILIA), 4), "to",
    round(max(dominance_full$mean_MUSSISMILIA), 4)
  ))
}

# ============================================================================
# SUMMARIES BY REEF x HAB AND ARCH x HAB
# ============================================================================

cat("\n=== Creating summaries ===\n")

summary_reef_hab <- dominance_full %>%
  group_by(REEF, HAB) %>%
  summarise(
    n_obs              = n(),
    mean_P_macro_coral = mean(P_macro_coral, na.rm = TRUE),
    min_P_macro_coral  = min(P_macro_coral,  na.rm = TRUE),
    max_P_macro_coral  = max(P_macro_coral,  na.rm = TRUE),
    sd_P_macro_coral   = sd(P_macro_coral,   na.rm = TRUE),
    mean_MUSSISMILIA   = mean(mean_MUSSISMILIA, na.rm = TRUE),
    mean_MACROALGAE    = mean(mean_MACROALGAE,  na.rm = TRUE),
    mean_CCA           = mean(mean_CCA,         na.rm = TRUE),
    dominant_status    = names(which.max(table(status))),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_P_macro_coral))

summary_arch_hab <- dominance_full %>%
  group_by(ARCH, HAB) %>%
  summarise(
    n_sites            = n_distinct(SITE),
    n_obs              = n(),
    mean_P_macro_coral = mean(P_macro_coral, na.rm = TRUE),
    min_P_macro_coral  = min(P_macro_coral,  na.rm = TRUE),
    max_P_macro_coral  = max(P_macro_coral,  na.rm = TRUE),
    mean_MUSSISMILIA   = mean(mean_MUSSISMILIA, na.rm = TRUE),
    mean_MACROALGAE    = mean(mean_MACROALGAE,  na.rm = TRUE),
    mean_CCA           = mean(mean_CCA,         na.rm = TRUE),
    # Marginal (population-average) dominance — preferred for general claims
    mean_P_macro_marginal = if ("P_macro_coral_marginal" %in% names(pick(everything())))
      mean(P_macro_coral_marginal, na.rm = TRUE) else NA_real_,
    mean_P_cca_marginal = if ("P_cca_coral_marginal" %in% names(pick(everything())))
      mean(P_cca_coral_marginal, na.rm = TRUE) else NA_real_,
    .groups = "drop"
  )

# Extremes
refugia_sites <- dominance_full %>%
  arrange(P_macro_coral) %>%
  head(10) %>%
  select(REEF, SITE, HAB, ARCH, P_macro_coral,
         mean_MUSSISMILIA, mean_MACROALGAE, mean_CCA, status)

stress_sites <- dominance_full %>%
  arrange(desc(P_macro_coral)) %>%
  head(10) %>%
  select(REEF, SITE, HAB, ARCH, P_macro_coral,
         mean_MUSSISMILIA, mean_MACROALGAE, mean_CCA, status)

# ============================================================================
# SAVE RESULTS
# ============================================================================

cat("\n=== Saving results to:", output_dir, "===\n")

write.csv(dominance_full,
          file.path(output_dir, "jsdm_dominance_probabilities_full.csv"),
          row.names = FALSE)

write.csv(summary_reef_hab,
          file.path(output_dir, "jsdm_dominance_summary_reef_hab.csv"),
          row.names = FALSE)

write.csv(summary_arch_hab,
          file.path(output_dir, "jsdm_dominance_summary_arch_hab.csv"),
          row.names = FALSE)

write.csv(refugia_sites,
          file.path(output_dir, "jsdm_refugia_sites.csv"),
          row.names = FALSE)

write.csv(stress_sites,
          file.path(output_dir, "jsdm_stress_sites.csv"),
          row.names = FALSE)

cat("  All CSVs saved.\n")

# ============================================================================
# PRINT KEY FINDINGS
# ============================================================================

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("KEY FINDINGS: DOMINANCE PROBABILITIES\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("TOP 10 CORAL REFUGIA (lowest P_macro > coral):\n")
print(refugia_sites[, c("REEF", "SITE", "HAB", "ARCH", "P_macro_coral",
                         "mean_MUSSISMILIA", "mean_MACROALGAE")])

cat("\nTOP 10 MACROALGAE-DOMINATED ZONES (highest P_macro > coral):\n")
print(stress_sites[, c("REEF", "SITE", "HAB", "ARCH", "P_macro_coral",
                        "mean_MUSSISMILIA", "mean_MACROALGAE")])

cat("\nSUMMARY BY ARCH x HAB (CONDITIONAL — includes REEF random effects):\n")
print(summary_arch_hab[, c("ARCH", "HAB", "n_sites", "n_obs",
                             "mean_P_macro_coral", "mean_MUSSISMILIA",
                             "mean_MACROALGAE")])

if ("mean_P_macro_marginal" %in% names(summary_arch_hab)) {
  cat("\nSUMMARY BY ARCH x HAB (MARGINAL — population-average, preferred for general claims):\n")
  print(summary_arch_hab[, c("ARCH", "HAB", "n_sites", "n_obs",
                               "mean_P_macro_marginal", "mean_P_cca_marginal")])
}

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("Output directory: %s\n", output_dir))
cat("Script 07 completed successfully!\n")
