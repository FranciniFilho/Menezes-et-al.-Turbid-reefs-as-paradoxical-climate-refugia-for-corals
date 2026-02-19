# ============================================================================
# 06_WINNER_DETAILED_REPORT_YEAR_RE.R
# ============================================================================
# Build detailed winner-level outputs for manuscript support:
#  - CSV tables per winner (metrics, fixed effects, random effects, key effects)
#  - Combined CSV tables across winners
#  - Markdown report with interpretative text and tables
# ============================================================================

suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(readr)
  library(stringr)
})

BASE_OUTPUT_DIR_YEAR_RE <- "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/"
GLOBAL_WINNER_DIR <- file.path(BASE_OUTPUT_DIR_YEAR_RE, "Global_Winners_YEAR_RE")
DETAIL_DIR <- file.path(BASE_OUTPUT_DIR_YEAR_RE, "Winner_Detailed_Reports_YEAR_RE")
SUPPLEMENT_DIR <- file.path(DETAIL_DIR, "Supplement_Ready")

dir.create(DETAIL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPPLEMENT_DIR, recursive = TRUE, showWarnings = FALSE)

combined_file_map <- list(
  COVER = file.path(BASE_OUTPUT_DIR_YEAR_RE, "ZOIB_Abundance_YEAR_RE", "ZOIB_YEAR_RE_combined_results.csv"),
  RGR = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_RGR", "RGR_combined_results.csv"),
  HEALTH_PC1 = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_Health_YEAR_RE", "Health_YEAR_RE_combined_results.csv"),
  HEALTH_PC2 = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_Health_YEAR_RE", "Health_YEAR_RE_combined_results.csv"),
  JSDM = file.path(BASE_OUTPUT_DIR_YEAR_RE, "JSDM_Dirichlet_YEAR_RE", "JSDM_YEAR_RE_combined_results.csv")
)

response_lookup_map <- c(
  COVER = "COVER",
  RGR = "RGR",
  HEALTH_PC1 = "HEALTH_PC1",
  HEALTH_PC2 = "HEALTH_PC2",
  JSDM = "COMMUNITY"
)

response_order <- c("COVER", "RGR", "HEALTH_PC1", "HEALTH_PC2", "JSDM")

response_label_map <- c(
  COVER = "Abundance (COVER)",
  RGR = "Growth (RGR)",
  HEALTH_PC1 = "Health PC1",
  HEALTH_PC2 = "Health PC2",
  JSDM = "Community Composition (JSDM)"
)

apply_response_order <- function(df) {
  if (is.null(df) || nrow(df) == 0 || !("Response" %in% names(df))) return(df)
  df %>%
    mutate(Response = factor(Response, levels = response_order)) %>%
    arrange(Response)
}

safe_read_csv <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readr::read_csv(path, show_col_types = FALSE), error = function(e) NULL)
}

safe_num <- function(x, digits = 4) {
  ifelse(is.na(x), NA_character_, format(round(as.numeric(x), digits), nsmall = digits, trim = TRUE))
}

md_table <- function(df, digits = 4, max_rows = 20) {
  if (is.null(df) || nrow(df) == 0) return("_No rows available._")

  if (nrow(df) > max_rows) {
    df <- df[seq_len(max_rows), , drop = FALSE]
  }

  fmt <- df
  for (j in seq_along(fmt)) {
    if (is.numeric(fmt[[j]])) {
      fmt[[j]] <- safe_num(fmt[[j]], digits = digits)
    } else {
      fmt[[j]] <- as.character(fmt[[j]])
    }
  }

  header <- paste0("| ", paste(names(fmt), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(fmt)), collapse = " | "), " |")
  rows <- apply(fmt, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  paste(c(header, sep, rows), collapse = "\n")
}

pretty_parameter <- function(x) {
  out <- x
  out <- str_replace(out, "^muMUSSISMILIAprop_", "MUSSISMILIA: ")
  out <- str_replace(out, "^muTURFprop_", "TURF: ")
  out <- str_replace(out, "^muCCAprop_", "CCA: ")
  out <- str_replace(out, "^muCYANOprop_", "CYANO: ")
  out <- str_replace(out, "^muMACROALGAEprop_", "MACROALGAE: ")
  out <- str_replace_all(out, "_", " ")
  out <- str_replace_all(out, ":", " x ")
  out
}

extract_fixed_effects <- function(model, response_label) {
  fx <- as.data.frame(brms::fixef(model, summary = TRUE))
  sm <- as.data.frame(summary(model)$fixed)

  out <- data.frame(
    Response = response_label,
    Parameter = rownames(fx),
    Parameter_Label = pretty_parameter(rownames(fx)),
    Estimate = fx$Estimate,
    Est_Error = fx$Est.Error,
    CI_low = fx$Q2.5,
    CI_high = fx$Q97.5,
    stringsAsFactors = FALSE
  )

  if (all(c("Rhat", "Bulk_ESS", "Tail_ESS") %in% names(sm))) {
    out$Rhat <- sm$Rhat[match(out$Parameter, rownames(sm))]
    out$ESS_Bulk <- sm$Bulk_ESS[match(out$Parameter, rownames(sm))]
    out$ESS_Tail <- sm$Tail_ESS[match(out$Parameter, rownames(sm))]
  } else {
    out$Rhat <- NA_real_
    out$ESS_Bulk <- NA_real_
    out$ESS_Tail <- NA_real_
  }

  out <- out %>%
    mutate(
      Credible_NonZero = CI_low * CI_high > 0,
      Direction = case_when(
        Estimate > 0 ~ "Positive",
        Estimate < 0 ~ "Negative",
        TRUE ~ "Neutral"
      ),
      Abs_Estimate = abs(Estimate)
    )

  out
}

extract_random_effects <- function(model, response_label) {
  vc <- brms::VarCorr(model, summary = TRUE)
  groups <- names(vc)
  if (length(groups) == 0) return(data.frame())

  bind_rows(lapply(groups, function(g) {
    sd_tbl <- as.data.frame(vc[[g]]$sd)
    if (nrow(sd_tbl) == 0) return(data.frame())

    data.frame(
      Response = response_label,
      Group = g,
      Term = rownames(sd_tbl),
      Estimate_SD = sd_tbl$Estimate,
      Est_Error_SD = sd_tbl$Est.Error,
      CI_low_SD = sd_tbl$Q2.5,
      CI_high_SD = sd_tbl$Q97.5,
      stringsAsFactors = FALSE
    )
  }))
}

extract_divergences <- function(model) {
  np <- tryCatch(brms::nuts_params(model), error = function(e) NULL)
  if (is.null(np) || nrow(np) == 0) return(NA_integer_)
  sum(np$Parameter == "divergent__" & np$Value == 1, na.rm = TRUE)
}

extract_r2 <- function(model) {
  fam <- tolower(model$family$family)
  if (fam == "dirichlet") return(c(NA_real_, NA_real_, NA_real_))

  r2 <- tryCatch(brms::bayes_R2(model, summary = TRUE), error = function(e) NULL)
  if (is.null(r2)) return(c(NA_real_, NA_real_, NA_real_))

  if (is.matrix(r2) && ncol(r2) >= 4) {
    return(c(r2[1, 1], r2[1, 3], r2[1, 4]))
  }

  c(NA_real_, NA_real_, NA_real_)
}

build_interpretation_text <- function(metrics_row, fixed_df, random_df) {
  response <- metrics_row$Response
  family <- metrics_row$Family

  non_intercept <- fixed_df %>% filter(!str_detect(Parameter, "Intercept"))
  sig <- non_intercept %>% filter(Credible_NonZero)
  top_sig <- sig %>% arrange(desc(Abs_Estimate)) %>% head(5)

  n_terms <- nrow(non_intercept)
  n_sig <- nrow(sig)
  n_pos <- sum(sig$Estimate > 0, na.rm = TRUE)
  n_neg <- sum(sig$Estimate < 0, na.rm = TRUE)

  year_re <- random_df %>% filter(Group == "YEAR")
  reef_re <- random_df %>% filter(Group == "REEF")

  lines <- c(
    sprintf("The selected winner for **%s** is **%s** at **%s** (family: `%s`).",
            response, metrics_row$Model, metrics_row$Winner_CV, family),
    sprintf("Model fit metrics indicate LOOIC = **%.3f** and %s divergences, with max R-hat %.3f and minimum bulk ESS %.0f.",
            metrics_row$LOOIC,
            ifelse(is.na(metrics_row$N_Divergent), "unknown", as.character(metrics_row$N_Divergent)),
            metrics_row$Rhat_Max,
            metrics_row$ESS_Bulk_Min),
    sprintf("Across %d non-intercept fixed terms, %d show credible non-zero effects (%d positive, %d negative).",
            n_terms, n_sig, n_pos, n_neg)
  )

  if (nrow(top_sig) > 0) {
    top_lines <- apply(top_sig, 1, function(r) {
      sprintf("- %s: estimate %.3f (95%% CrI %.3f to %.3f)",
              r[["Parameter_Label"]], as.numeric(r[["Estimate"]]),
              as.numeric(r[["CI_low"]]), as.numeric(r[["CI_high"]]))
    })
    lines <- c(lines, "Key strongest credible effects:", top_lines)
  } else {
    lines <- c(lines, "No non-intercept fixed effect showed a 95% credible interval fully away from zero.")
  }

  if (nrow(year_re) > 0) {
    lines <- c(lines,
               sprintf("YEAR random-effect SD range: %.3f to %.3f.",
                       min(year_re$Estimate_SD, na.rm = TRUE),
                       max(year_re$Estimate_SD, na.rm = TRUE)))
  }

  if (nrow(reef_re) > 0) {
    lines <- c(lines,
               sprintf("REEF random-effect SD range: %.3f to %.3f.",
                       min(reef_re$Estimate_SD, na.rm = TRUE),
                       max(reef_re$Estimate_SD, na.rm = TRUE)))
  }

  paste(lines, collapse = "\n")
}

build_supplement_outputs <- function(metrics_all, fixed_all, random_all, key_all) {
  metrics_all <- apply_response_order(metrics_all)
  fixed_all <- apply_response_order(fixed_all)
  random_all <- apply_response_order(random_all)
  key_all <- apply_response_order(key_all)

  s1 <- metrics_all %>%
    mutate(Outcome = unname(response_label_map[as.character(Response)])) %>%
    transmute(
      Outcome,
      Winner_Model = Model,
      Temporal_Window = Winner_CV,
      Family,
      N_Obs,
      LOOIC,
      SE_LOOIC,
      Converged,
      Max_Rhat = Rhat_Max,
      Min_ESS_Bulk = ESS_Bulk_Min,
      N_Divergent,
      R2_Mean,
      R2_CI_low,
      R2_CI_high
    )

  s2 <- key_all %>%
    mutate(Outcome = unname(response_label_map[as.character(Response)])) %>%
    group_by(Outcome) %>%
    arrange(desc(Abs_Estimate), .by_group = TRUE) %>%
    slice_head(n = 12) %>%
    ungroup() %>%
    transmute(
      Outcome,
      Parameter = Parameter_Label,
      Estimate,
      CrI_95_low = CI_low,
      CrI_95_high = CI_high,
      Credible_NonZero,
      Direction,
      Rhat,
      ESS_Bulk
    )

  s3 <- random_all %>%
    mutate(Outcome = unname(response_label_map[as.character(Response)])) %>%
    group_by(Outcome, Group) %>%
    summarise(
      N_Terms = n(),
      Mean_SD = mean(Estimate_SD, na.rm = TRUE),
      Min_SD = min(Estimate_SD, na.rm = TRUE),
      Max_SD = max(Estimate_SD, na.rm = TRUE),
      Mean_CI_Width = mean(CI_high_SD - CI_low_SD, na.rm = TRUE),
      .groups = "drop"
    )

  s4 <- fixed_all %>%
    mutate(Outcome = unname(response_label_map[as.character(Response)])) %>%
    arrange(Response, desc(Abs_Estimate)) %>%
    transmute(
      Outcome,
      Parameter = Parameter_Label,
      Estimate,
      Est_Error,
      CrI_95_low = CI_low,
      CrI_95_high = CI_high,
      Credible_NonZero,
      Direction,
      Rhat,
      ESS_Bulk,
      ESS_Tail
    )

  s1_path <- file.path(SUPPLEMENT_DIR, "Table_S1_Winner_Overview.csv")
  s2_path <- file.path(SUPPLEMENT_DIR, "Table_S2_Key_Fixed_Effects.csv")
  s3_path <- file.path(SUPPLEMENT_DIR, "Table_S3_Random_Effect_Summary.csv")
  s4_path <- file.path(SUPPLEMENT_DIR, "Table_S4_All_Fixed_Effects.csv")
  md_path <- file.path(SUPPLEMENT_DIR, "SUPPLEMENT_WINNER_TABLES_YEAR_RE.md")

  write.csv(s1, s1_path, row.names = FALSE)
  write.csv(s2, s2_path, row.names = FALSE)
  write.csv(s3, s3_path, row.names = FALSE)
  write.csv(s4, s4_path, row.names = FALSE)

  md_lines <- c(
    "# Supplementary Tables: YEAR_RE Winner Models",
    "",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "This supplement-ready package presents the final winner models in publication appendix format.",
    "",
    "## Table S1. Winner overview across outcomes",
    "",
    md_table(s1, digits = 4, max_rows = 10),
    "",
    sprintf("Full CSV: `%s`", s1_path),
    "",
    "## Table S2. Key fixed effects (top effects per outcome)",
    "",
    md_table(s2, digits = 4, max_rows = 30),
    "",
    sprintf("Full CSV: `%s`", s2_path),
    "",
    "## Table S3. Random-effect summary",
    "",
    md_table(s3, digits = 4, max_rows = 20),
    "",
    sprintf("Full CSV: `%s`", s3_path),
    "",
    "## Table S4. Complete fixed-effects table",
    "",
    "Table S4 can be large and is provided primarily as CSV for direct appendix import.",
    "",
    sprintf("Full CSV: `%s`", s4_path),
    ""
  )

  writeLines(md_lines, md_path, useBytes = TRUE)

  list(
    table_s1 = s1_path,
    table_s2 = s2_path,
    table_s3 = s3_path,
    table_s4 = s4_path,
    supplement_md = md_path
  )
}

build_winner_detailed_reports_year_re <- function() {
  cat("\n============================================================================\n")
  cat("BUILDING DETAILED WINNER REPORTS (YEAR_RE)\n")
  cat("============================================================================\n")

  summary_path <- file.path(GLOBAL_WINNER_DIR, "GLOBAL_WINNER_SUMMARY_YEAR_RE.csv")
  if (!file.exists(summary_path)) {
    stop(sprintf("Global summary not found: %s", summary_path))
  }

  global_summary <- readr::read_csv(summary_path, show_col_types = FALSE)
  winners <- global_summary %>% filter(Status == "OK")
  if (nrow(winners) == 0) stop("No global winners with Status == OK")

  all_metrics <- list()
  all_fixed <- list()
  all_random <- list()
  all_key <- list()

  md_lines <- c(
    "# YEAR_RE Winner Detailed Report",
    "",
    sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    "This report consolidates detailed posterior outputs for each global winner model, including diagnostics, fixed effects, random effects, and interpretation-ready summaries.",
    ""
  )

  for (i in seq_len(nrow(winners))) {
    row <- winners[i, ]
    response <- as.character(row$Response)
    model_path <- as.character(row$Global_Path)
    model_name <- as.character(row$Model)
    winner_cv <- as.character(row$Winner_CV)

    cat(sprintf("\n>>> Processing winner: %s (%s, %s)\n", response, model_name, winner_cv))

    if (!file.exists(model_path)) {
      cat(sprintf("  WARNING: Winner file missing: %s\n", model_path))
      next
    }

    model <- readRDS(model_path)
    sm <- summary(model)

    combined_tbl <- safe_read_csv(combined_file_map[[response]])
    lookup_response <- response_lookup_map[[response]]

    combined_row <- NULL
    if (!is.null(combined_tbl)) {
      combined_row <- combined_tbl %>%
        mutate(
          Response = toupper(as.character(Response)),
          Model = toupper(as.character(Model)),
          CV = toupper(as.character(CV))
        ) %>%
        filter(Response == toupper(lookup_response),
               Model == toupper(model_name),
               CV == toupper(winner_cv)) %>%
        slice(1)
    }

    fixed_df <- extract_fixed_effects(model, response)
    random_df <- extract_random_effects(model, response)

    key_df <- fixed_df %>%
      filter(!str_detect(Parameter, "Intercept"), Credible_NonZero) %>%
      arrange(desc(Abs_Estimate))
    if (nrow(key_df) == 0) {
      key_df <- fixed_df %>%
        filter(!str_detect(Parameter, "Intercept")) %>%
        arrange(desc(Abs_Estimate))
    }
    key_df <- head(key_df, 12)

    r2_vals <- extract_r2(model)
    rhat_max <- suppressWarnings(max(c(sm$fixed[, "Rhat"], unlist(lapply(sm$random, function(x) x[, "Rhat"]))), na.rm = TRUE))
    ess_bulk_min <- suppressWarnings(min(c(sm$fixed[, "Bulk_ESS"], unlist(lapply(sm$random, function(x) x[, "Bulk_ESS"]))), na.rm = TRUE))
    ess_tail_min <- suppressWarnings(min(c(sm$fixed[, "Tail_ESS"], unlist(lapply(sm$random, function(x) x[, "Tail_ESS"]))), na.rm = TRUE))

    metrics_row <- data.frame(
      Response = response,
      Family = model$family$family,
      Model = model_name,
      Winner_CV = winner_cv,
      Formula = paste(deparse(formula(model)$formula), collapse = " "),
      N_Obs = nrow(model$data),
      LOOIC = as.numeric(row$LOOIC),
      SE_LOOIC = if (!is.null(combined_row) && "SE_LOOIC" %in% names(combined_row)) as.numeric(combined_row$SE_LOOIC[1]) else NA_real_,
      Converged = if (!is.null(combined_row) && "Converged" %in% names(combined_row)) as.character(combined_row$Converged[1]) else NA_character_,
      Rhat_Max = rhat_max,
      ESS_Bulk_Min = ess_bulk_min,
      ESS_Tail_Min = ess_tail_min,
      N_Divergent = extract_divergences(model),
      R2_Mean = r2_vals[1],
      R2_CI_low = r2_vals[2],
      R2_CI_high = r2_vals[3],
      Source_Path = as.character(row$Source_Path),
      Global_Path = model_path,
      stringsAsFactors = FALSE
    )

    resp_dir <- file.path(DETAIL_DIR, response)
    dir.create(resp_dir, recursive = TRUE, showWarnings = FALSE)

    metrics_csv <- file.path(resp_dir, sprintf("%s_winner_metrics.csv", tolower(response)))
    fixed_csv <- file.path(resp_dir, sprintf("%s_fixed_effects.csv", tolower(response)))
    random_csv <- file.path(resp_dir, sprintf("%s_random_effects.csv", tolower(response)))
    key_csv <- file.path(resp_dir, sprintf("%s_key_effects.csv", tolower(response)))

    write.csv(metrics_row, metrics_csv, row.names = FALSE)
    write.csv(fixed_df, fixed_csv, row.names = FALSE)
    write.csv(random_df, random_csv, row.names = FALSE)
    write.csv(key_df, key_csv, row.names = FALSE)

    all_metrics[[response]] <- metrics_row
    all_fixed[[response]] <- fixed_df
    all_random[[response]] <- random_df
    all_key[[response]] <- key_df

    interpretation <- build_interpretation_text(metrics_row, fixed_df, random_df)

    md_lines <- c(
      md_lines,
      sprintf("## %s", response),
      "",
      interpretation,
      "",
      "### Winner Metrics",
      "",
      md_table(metrics_row %>%
                 select(Response, Family, Model, Winner_CV, N_Obs, LOOIC, SE_LOOIC, Converged, Rhat_Max, ESS_Bulk_Min, N_Divergent, R2_Mean, R2_CI_low, R2_CI_high),
               digits = 4, max_rows = 5),
      "",
      "### Top Effects Table",
      "",
      md_table(key_df %>%
                 select(Parameter_Label, Estimate, CI_low, CI_high, Credible_NonZero, Rhat, ESS_Bulk),
               digits = 4, max_rows = 12),
      "",
      "### Random Effects (SD)",
      "",
      md_table(random_df %>%
                 select(Group, Term, Estimate_SD, CI_low_SD, CI_high_SD),
               digits = 4, max_rows = 20),
      "",
      sprintf("CSV files: `%s`, `%s`, `%s`, `%s`", metrics_csv, fixed_csv, random_csv, key_csv),
      ""
    )
  }

  metrics_all <- bind_rows(all_metrics)
  fixed_all <- bind_rows(all_fixed)
  random_all <- bind_rows(all_random)
  key_all <- bind_rows(all_key)

  metrics_all_path <- file.path(DETAIL_DIR, "WINNER_METRICS_ALL.csv")
  fixed_all_path <- file.path(DETAIL_DIR, "WINNER_FIXED_EFFECTS_ALL.csv")
  random_all_path <- file.path(DETAIL_DIR, "WINNER_RANDOM_EFFECTS_ALL.csv")
  key_all_path <- file.path(DETAIL_DIR, "WINNER_KEY_EFFECTS_ALL.csv")
  report_path <- file.path(DETAIL_DIR, "WINNER_DETAILED_REPORT_YEAR_RE.md")

  write.csv(metrics_all, metrics_all_path, row.names = FALSE)
  write.csv(fixed_all, fixed_all_path, row.names = FALSE)
  write.csv(random_all, random_all_path, row.names = FALSE)
  write.csv(key_all, key_all_path, row.names = FALSE)

  supplement_outputs <- build_supplement_outputs(metrics_all, fixed_all, random_all, key_all)

  md_lines <- c(
    md_lines,
    "## Consolidated CSV Outputs",
    "",
    sprintf("- `%s`", metrics_all_path),
    sprintf("- `%s`", fixed_all_path),
    sprintf("- `%s`", random_all_path),
    sprintf("- `%s`", key_all_path),
    ""
  )

  writeLines(md_lines, report_path, useBytes = TRUE)

  cat("\nDetailed winner outputs saved:\n")
  cat(sprintf("  - %s\n", metrics_all_path))
  cat(sprintf("  - %s\n", fixed_all_path))
  cat(sprintf("  - %s\n", random_all_path))
  cat(sprintf("  - %s\n", key_all_path))
  cat(sprintf("  - %s\n", report_path))
  cat("Supplement-ready outputs saved:\n")
  cat(sprintf("  - %s\n", supplement_outputs$table_s1))
  cat(sprintf("  - %s\n", supplement_outputs$table_s2))
  cat(sprintf("  - %s\n", supplement_outputs$table_s3))
  cat(sprintf("  - %s\n", supplement_outputs$table_s4))
  cat(sprintf("  - %s\n", supplement_outputs$supplement_md))

  invisible(list(
    metrics_csv = metrics_all_path,
    fixed_csv = fixed_all_path,
    random_csv = random_all_path,
    key_csv = key_all_path,
    report_md = report_path,
    supplement_s1_csv = supplement_outputs$table_s1,
    supplement_s2_csv = supplement_outputs$table_s2,
    supplement_s3_csv = supplement_outputs$table_s3,
    supplement_s4_csv = supplement_outputs$table_s4,
    supplement_md = supplement_outputs$supplement_md
  ))
}

called_args <- commandArgs(trailingOnly = FALSE)
is_direct_call <- any(grepl("^--file=.*06_WINNER_DETAILED_REPORT_YEAR_RE\\.R$", called_args))

if (!interactive() && is_direct_call) {
  build_winner_detailed_reports_year_re()
}
