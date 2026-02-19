# ============================================================================
# 05_CONSOLIDATE_GLOBAL_WINNERS_YEAR_RE.R
# ============================================================================
# Build one global winner per response for final publication figures.
#
# Strategy:
#   1) Use formal CV-window LOO summaries (CV_WINDOW_LOO_SUMMARY_*.csv) to get
#      the best CV window per model specification.
#   2) For each response, evaluate each model at its CV winner using the
#      combined results table.
#   3) Select the best global model (minimum LOOIC) and copy the matching .rds
#      to WINNER_GLOBAL_<response>.rds.
# ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

BASE_OUTPUT_DIR_YEAR_RE <- "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/"
GLOBAL_WINNER_DIR <- file.path(BASE_OUTPUT_DIR_YEAR_RE, "Global_Winners_YEAR_RE")
dir.create(GLOBAL_WINNER_DIR, recursive = TRUE, showWarnings = FALSE)

read_csv_flexible <- function(path) {
  if (!file.exists(path)) return(NULL)
  out <- tryCatch(
    readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)
  tryCatch(read.csv(path, stringsAsFactors = FALSE), error = function(e) NULL)
}

to_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  as.character(x) %in% c("TRUE", "True", "true", "1")
}

build_response_config <- function() {
  list(
    COVER = list(
      label = "COVER",
      summary_file = file.path(BASE_OUTPUT_DIR_YEAR_RE, "ZOIB_Abundance_YEAR_RE", "CV_WINDOW_LOO_SUMMARY_ZOIB.csv"),
      combined_file = file.path(BASE_OUTPUT_DIR_YEAR_RE, "ZOIB_Abundance_YEAR_RE", "ZOIB_YEAR_RE_combined_results.csv"),
      model_dir = file.path(BASE_OUTPUT_DIR_YEAR_RE, "ZOIB_Abundance_YEAR_RE"),
      model_prefix = "zoib_year_re",
      source_subdir = function(cv) cv
    ),
    RGR = list(
      label = "RGR",
      summary_file = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_RGR", "CV_WINDOW_LOO_SUMMARY_RGR.csv"),
      combined_file = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_RGR", "RGR_combined_results.csv"),
      model_dir = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_RGR"),
      model_prefix = "gaussian_rgr",
      source_subdir = function(cv) file.path("RGR", cv)
    ),
    HEALTH_PC1 = list(
      label = "HEALTH_PC1",
      summary_file = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_Health_YEAR_RE", "CV_WINDOW_LOO_SUMMARY_HEALTH.csv"),
      combined_file = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_Health_YEAR_RE", "Health_YEAR_RE_combined_results.csv"),
      model_dir = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_Health_YEAR_RE"),
      model_prefix = "gaussian_health_pc1_year_re",
      source_subdir = function(cv) file.path("HEALTH_PC1", cv)
    ),
    HEALTH_PC2 = list(
      label = "HEALTH_PC2",
      summary_file = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_Health_YEAR_RE", "CV_WINDOW_LOO_SUMMARY_HEALTH.csv"),
      combined_file = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_Health_YEAR_RE", "Health_YEAR_RE_combined_results.csv"),
      model_dir = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_Health_YEAR_RE"),
      model_prefix = "gaussian_health_pc2_year_re",
      source_subdir = function(cv) file.path("HEALTH_PC2", cv)
    ),
    JSDM = list(
      label = "JSDM",
      summary_file = file.path(BASE_OUTPUT_DIR_YEAR_RE, "JSDM_Dirichlet_YEAR_RE", "CV_WINDOW_LOO_SUMMARY_JSDM.csv"),
      combined_file = file.path(BASE_OUTPUT_DIR_YEAR_RE, "JSDM_Dirichlet_YEAR_RE", "JSDM_YEAR_RE_combined_results.csv"),
      model_dir = file.path(BASE_OUTPUT_DIR_YEAR_RE, "JSDM_Dirichlet_YEAR_RE"),
      model_prefix = "jsdm_year_re",
      source_subdir = function(cv) cv
    )
  )
}

normalize_summary_rows <- function(summary_df, response_label) {
  if (is.null(summary_df) || nrow(summary_df) == 0) return(NULL)

  key_col <- if ("Model_Key" %in% names(summary_df)) "Model_Key" else if ("Model" %in% names(summary_df)) "Model" else NULL
  if (is.null(key_col) || !("Winner_CV" %in% names(summary_df)) || !("Status" %in% names(summary_df))) return(NULL)

  out <- summary_df %>%
    mutate(
      model_key = as.character(.data[[key_col]]),
      model_name = ifelse(grepl("__", model_key), sub("^.*__", "", model_key), model_key),
      response_in_key = ifelse(grepl("__", model_key), sub("__.*$", "", model_key), response_label),
      response_in_key = toupper(response_in_key),
      Winner_CV = as.character(Winner_CV),
      Status = as.character(Status)
    )

  out <- out %>%
    filter(Status == "OK", !is.na(Winner_CV), Winner_CV != "", Winner_CV != "NA")

  if (response_label == "COVER") {
    return(out)
  }

  out %>% filter(response_in_key == response_label)
}

select_global_winner <- function(cfg) {
  summary_df <- read_csv_flexible(cfg$summary_file)
  combined_df <- read_csv_flexible(cfg$combined_file)
  if (is.null(combined_df)) return(NULL)

  candidate_rows <- normalize_summary_rows(summary_df, cfg$label)

  if (!("Response" %in% names(combined_df)) || !("Model" %in% names(combined_df)) ||
      !("CV" %in% names(combined_df)) || !("LOOIC" %in% names(combined_df))) {
    return(NULL)
  }

  if (!("Converged" %in% names(combined_df))) {
    combined_df$Converged <- TRUE
  }

  response_filter <- ifelse(cfg$label == "JSDM", "COMMUNITY", cfg$label)

  converged_tbl <- combined_df %>%
    mutate(
      Response = toupper(as.character(Response)),
      Model = toupper(as.character(Model)),
      CV = toupper(as.character(CV)),
      Converged = to_logical_safe(Converged),
      LOOIC = suppressWarnings(as.numeric(LOOIC))
    ) %>%
    filter(Response == response_filter, Converged, !is.na(LOOIC)) %>%
    select(Model, CV, LOOIC, SE_LOOIC)

  if (nrow(converged_tbl) == 0) return(NULL)

  eval_tbl <- NULL

  if (!is.null(candidate_rows) && nrow(candidate_rows) > 0) {
    eval_tbl <- candidate_rows %>%
      mutate(
        model_name = toupper(model_name),
        Winner_CV = toupper(Winner_CV)
      ) %>%
      left_join(
        converged_tbl,
        by = c("model_name" = "Model", "Winner_CV" = "CV")
      ) %>%
      filter(!is.na(LOOIC)) %>%
      arrange(LOOIC)
  }

  if (is.null(eval_tbl) || nrow(eval_tbl) == 0) {
    cat(sprintf("  [Fallback] %s: selecting global winner from converged combined results.\n", cfg$label))
    eval_tbl <- converged_tbl %>%
      transmute(
        model_name = Model,
        Winner_CV = CV,
        LOOIC = LOOIC,
        SE_LOOIC = SE_LOOIC
      ) %>%
      arrange(LOOIC)
  }

  if (nrow(eval_tbl) == 0) return(NULL)

  best <- eval_tbl[1, , drop = FALSE]
  model_name <- tolower(best$model_name[[1]])
  winner_cv <- toupper(best$Winner_CV[[1]])

  src_file <- sprintf("%s_%s_%s.rds", cfg$model_prefix, model_name, winner_cv)
  src_path <- file.path(cfg$model_dir, cfg$source_subdir(winner_cv), src_file)

  if (!file.exists(src_path)) {
    return(list(
      status = "SOURCE_NOT_FOUND",
      response = cfg$label,
      model = best$model_name[[1]],
      cv = winner_cv,
      looic = best$LOOIC[[1]],
      source_path = src_path,
      global_path = NA_character_
    ))
  }

  global_name <- sprintf("WINNER_GLOBAL_%s.rds", tolower(cfg$label))
  global_path <- file.path(GLOBAL_WINNER_DIR, global_name)
  ok <- file.copy(src_path, global_path, overwrite = TRUE)

  list(
    status = ifelse(ok, "OK", "COPY_FAILED"),
    response = cfg$label,
    model = best$model_name[[1]],
    cv = winner_cv,
    looic = best$LOOIC[[1]],
    source_path = src_path,
    global_path = global_path
  )
}

build_global_winners_year_re <- function() {
  cat("\n============================================================================\n")
  cat("CONSOLIDATING GLOBAL WINNERS (YEAR_RE)\n")
  cat("============================================================================\n")

  cfg_list <- build_response_config()
  results <- lapply(cfg_list, select_global_winner)

  result_tbl <- bind_rows(lapply(names(results), function(resp_name) {
    x <- results[[resp_name]]
    if (is.null(x)) {
      data.frame(
        Status = "MISSING_INPUTS",
        Response = resp_name,
        Model = NA_character_,
        Winner_CV = NA_character_,
        LOOIC = NA_real_,
        Source_Path = NA_character_,
        Global_Path = NA_character_,
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        Status = x$status,
        Response = x$response,
        Model = x$model,
        Winner_CV = x$cv,
        LOOIC = x$looic,
        Source_Path = x$source_path,
        Global_Path = x$global_path,
        stringsAsFactors = FALSE
      )
    }
  }))

  out_csv <- file.path(GLOBAL_WINNER_DIR, "GLOBAL_WINNER_SUMMARY_YEAR_RE.csv")
  write.csv(result_tbl, out_csv, row.names = FALSE)

  n_ok <- sum(result_tbl$Status == "OK", na.rm = TRUE)
  cat(sprintf("Global winners generated: %d\n", n_ok))
  cat(sprintf("Summary: %s\n", out_csv))

  invisible(result_tbl)
}

if (!interactive()) {
  build_global_winners_year_re()
}
