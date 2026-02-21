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
BASE_OUTPUT_DIR_YEAR_RE_PRIOR_SENS <- "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output_PRIOR_SENSITIVITY_FULLGRID_v1/"

normalize_prior_scenario <- function(prior_tag) {
  tag <- tolower(trimws(as.character(prior_tag)))
  if (tag %in% c("weaklyinformative", "weakly_informative", "wi")) return("WeaklyInformative")
  if (tag %in% c("informative", "inf")) return("Informative")
  stop(sprintf("Unknown prior scenario: %s", prior_tag))
}

resolve_namespace <- function(run_namespace = c("canonical", "prior_sens_fullgrid"),
                              prior_scenario_target = "WeaklyInformative") {
  run_namespace <- match.arg(run_namespace)
  prior_scenario_target <- normalize_prior_scenario(prior_scenario_target)

  if (run_namespace == "canonical") {
    return(list(
      run_namespace = run_namespace,
      prior_scenario_target = prior_scenario_target,
      scenario_tag_upper = toupper(gsub("[^A-Za-z0-9]", "", prior_scenario_target)),
      scenario_tag_lower = tolower(gsub("[^A-Za-z0-9]", "", prior_scenario_target)),
      base_dir = BASE_OUTPUT_DIR_YEAR_RE,
      global_winner_dir = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Global_Winners_YEAR_RE")
    ))
  }

  list(
    run_namespace = run_namespace,
    prior_scenario_target = prior_scenario_target,
    scenario_tag_upper = toupper(gsub("[^A-Za-z0-9]", "", prior_scenario_target)),
    scenario_tag_lower = tolower(gsub("[^A-Za-z0-9]", "", prior_scenario_target)),
    base_dir = BASE_OUTPUT_DIR_YEAR_RE_PRIOR_SENS,
    global_winner_dir = file.path(BASE_OUTPUT_DIR_YEAR_RE_PRIOR_SENS, "PS_FULLGRID_Global_Winners_YEAR_RE")
  )
}

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

build_response_config <- function(ns_cfg) {
  if (ns_cfg$run_namespace == "canonical") {
    return(list(
      COVER = list(
        label = "COVER",
        summary_file = file.path(ns_cfg$base_dir, "ZOIB_Abundance_YEAR_RE", "CV_WINDOW_LOO_SUMMARY_ZOIB.csv"),
        combined_file = file.path(ns_cfg$base_dir, "ZOIB_Abundance_YEAR_RE", "ZOIB_YEAR_RE_combined_results.csv"),
        model_dir = file.path(ns_cfg$base_dir, "ZOIB_Abundance_YEAR_RE"),
        model_prefix = "zoib_year_re",
        source_subdir = function(cv, prior_scenario) cv
      ),
      RGR = list(
        label = "RGR",
        summary_file = file.path(ns_cfg$base_dir, "Gaussian_RGR", "CV_WINDOW_LOO_SUMMARY_RGR.csv"),
        combined_file = file.path(ns_cfg$base_dir, "Gaussian_RGR", "RGR_combined_results.csv"),
        model_dir = file.path(ns_cfg$base_dir, "Gaussian_RGR"),
        model_prefix = "gaussian_rgr",
        source_subdir = function(cv, prior_scenario) file.path("RGR", cv)
      ),
      HEALTH_PC1 = list(
        label = "HEALTH_PC1",
        summary_file = file.path(ns_cfg$base_dir, "Gaussian_Health_YEAR_RE", "CV_WINDOW_LOO_SUMMARY_HEALTH.csv"),
        combined_file = file.path(ns_cfg$base_dir, "Gaussian_Health_YEAR_RE", "Health_YEAR_RE_combined_results.csv"),
        model_dir = file.path(ns_cfg$base_dir, "Gaussian_Health_YEAR_RE"),
        model_prefix = "gaussian_health_pc1_year_re",
        source_subdir = function(cv, prior_scenario) file.path("HEALTH_PC1", cv)
      ),
      HEALTH_PC2 = list(
        label = "HEALTH_PC2",
        summary_file = file.path(ns_cfg$base_dir, "Gaussian_Health_YEAR_RE", "CV_WINDOW_LOO_SUMMARY_HEALTH.csv"),
        combined_file = file.path(ns_cfg$base_dir, "Gaussian_Health_YEAR_RE", "Health_YEAR_RE_combined_results.csv"),
        model_dir = file.path(ns_cfg$base_dir, "Gaussian_Health_YEAR_RE"),
        model_prefix = "gaussian_health_pc2_year_re",
        source_subdir = function(cv, prior_scenario) file.path("HEALTH_PC2", cv)
      ),
      JSDM = list(
        label = "JSDM",
        summary_file = file.path(ns_cfg$base_dir, "JSDM_Dirichlet_YEAR_RE", "CV_WINDOW_LOO_SUMMARY_JSDM.csv"),
        combined_file = file.path(ns_cfg$base_dir, "JSDM_Dirichlet_YEAR_RE", "JSDM_YEAR_RE_combined_results.csv"),
        model_dir = file.path(ns_cfg$base_dir, "JSDM_Dirichlet_YEAR_RE"),
        model_prefix = "jsdm_year_re",
        source_subdir = function(cv, prior_scenario) cv
      )
    ))
  }

  list(
    COVER = list(
      label = "COVER",
      summary_file = file.path(ns_cfg$base_dir, "PS_FULLGRID_ZOIB_YEAR_RE", sprintf("CV_WINDOW_LOO_SUMMARY_ZOIB_%s_PS_FULLGRID.csv", ns_cfg$scenario_tag_upper)),
      combined_file = file.path(ns_cfg$base_dir, "PS_FULLGRID_ZOIB_YEAR_RE", "ZOIB_YEAR_RE_combined_results_all_priors_PS_FULLGRID.csv"),
      model_dir = file.path(ns_cfg$base_dir, "PS_FULLGRID_ZOIB_YEAR_RE"),
      model_prefix = "zoib_year_re",
      source_subdir = function(cv, prior_scenario) file.path(cv, prior_scenario)
    ),
    RGR = list(
      label = "RGR",
      summary_file = file.path(ns_cfg$base_dir, "PS_FULLGRID_GAUSSIAN_RGR_YEAR_RE", sprintf("CV_WINDOW_LOO_SUMMARY_RGR_%s_PS_FULLGRID.csv", ns_cfg$scenario_tag_upper)),
      combined_file = file.path(ns_cfg$base_dir, "PS_FULLGRID_GAUSSIAN_RGR_YEAR_RE", "RGR_combined_results_all_priors_PS_FULLGRID.csv"),
      model_dir = file.path(ns_cfg$base_dir, "PS_FULLGRID_GAUSSIAN_RGR_YEAR_RE"),
      model_prefix = "gaussian_rgr",
      source_subdir = function(cv, prior_scenario) file.path("RGR", cv, prior_scenario)
    ),
    HEALTH_PC1 = list(
      label = "HEALTH_PC1",
      summary_file = file.path(ns_cfg$base_dir, "PS_FULLGRID_GAUSSIAN_HEALTH_YEAR_RE", sprintf("CV_WINDOW_LOO_SUMMARY_HEALTH_%s_PS_FULLGRID.csv", ns_cfg$scenario_tag_upper)),
      combined_file = file.path(ns_cfg$base_dir, "PS_FULLGRID_GAUSSIAN_HEALTH_YEAR_RE", "Health_YEAR_RE_combined_results_all_priors_PS_FULLGRID.csv"),
      model_dir = file.path(ns_cfg$base_dir, "PS_FULLGRID_GAUSSIAN_HEALTH_YEAR_RE"),
      model_prefix = "gaussian_health_pc1_year_re",
      source_subdir = function(cv, prior_scenario) file.path("HEALTH_PC1", cv, prior_scenario)
    ),
    HEALTH_PC2 = list(
      label = "HEALTH_PC2",
      summary_file = file.path(ns_cfg$base_dir, "PS_FULLGRID_GAUSSIAN_HEALTH_YEAR_RE", sprintf("CV_WINDOW_LOO_SUMMARY_HEALTH_%s_PS_FULLGRID.csv", ns_cfg$scenario_tag_upper)),
      combined_file = file.path(ns_cfg$base_dir, "PS_FULLGRID_GAUSSIAN_HEALTH_YEAR_RE", "Health_YEAR_RE_combined_results_all_priors_PS_FULLGRID.csv"),
      model_dir = file.path(ns_cfg$base_dir, "PS_FULLGRID_GAUSSIAN_HEALTH_YEAR_RE"),
      model_prefix = "gaussian_health_pc2_year_re",
      source_subdir = function(cv, prior_scenario) file.path("HEALTH_PC2", cv, prior_scenario)
    ),
    JSDM = list(
      label = "JSDM",
      summary_file = file.path(ns_cfg$base_dir, "PS_FULLGRID_JSDM_DIRICHLET_YEAR_RE", sprintf("CV_WINDOW_LOO_SUMMARY_JSDM_%s_PS_FULLGRID.csv", ns_cfg$scenario_tag_upper)),
      combined_file = file.path(ns_cfg$base_dir, "PS_FULLGRID_JSDM_DIRICHLET_YEAR_RE", "JSDM_YEAR_RE_combined_results_all_priors_PS_FULLGRID.csv"),
      model_dir = file.path(ns_cfg$base_dir, "PS_FULLGRID_JSDM_DIRICHLET_YEAR_RE"),
      model_prefix = "jsdm_year_re",
      source_subdir = function(cv, prior_scenario) file.path(cv, prior_scenario)
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

select_global_winner <- function(cfg, ns_cfg) {
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

  if (ns_cfg$run_namespace == "prior_sens_fullgrid" && ("Prior_Scenario" %in% names(combined_df))) {
    combined_df <- combined_df %>%
      mutate(Prior_Scenario = as.character(Prior_Scenario)) %>%
      filter(Prior_Scenario == ns_cfg$prior_scenario_target)
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

  if (ns_cfg$run_namespace == "prior_sens_fullgrid") {
    src_file <- sprintf("%s_%s_%s_prior_%s.rds", cfg$model_prefix, model_name, winner_cv, ns_cfg$scenario_tag_lower)
  } else {
    src_file <- sprintf("%s_%s_%s.rds", cfg$model_prefix, model_name, winner_cv)
  }
  src_path <- file.path(cfg$model_dir, cfg$source_subdir(winner_cv, ns_cfg$prior_scenario_target), src_file)

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

  if (ns_cfg$run_namespace == "prior_sens_fullgrid") {
    global_name <- sprintf("WINNER_GLOBAL_%s_%s.rds", tolower(cfg$label), ns_cfg$scenario_tag_lower)
  } else {
    global_name <- sprintf("WINNER_GLOBAL_%s.rds", tolower(cfg$label))
  }
  global_path <- file.path(ns_cfg$global_winner_dir, global_name)
  ok <- file.copy(src_path, global_path, overwrite = TRUE)

  list(
    status = ifelse(ok, "OK", "COPY_FAILED"),
    response = cfg$label,
    model = best$model_name[[1]],
    cv = winner_cv,
    looic = best$LOOIC[[1]],
    prior_scenario = ns_cfg$prior_scenario_target,
    source_path = src_path,
    global_path = global_path
  )
}

build_global_winners_year_re <- function(run_namespace = "canonical",
                                         prior_scenario_target = "WeaklyInformative") {
  ns_cfg <- resolve_namespace(run_namespace, prior_scenario_target)
  dir.create(ns_cfg$global_winner_dir, recursive = TRUE, showWarnings = FALSE)

  cat("\n============================================================================\n")
  cat(sprintf("CONSOLIDATING GLOBAL WINNERS (YEAR_RE) - namespace: %s\n", ns_cfg$run_namespace))
  if (ns_cfg$run_namespace == "prior_sens_fullgrid") {
    cat(sprintf("Prior scenario target: %s\n", ns_cfg$prior_scenario_target))
  }
  cat("============================================================================\n")

  cfg_list <- build_response_config(ns_cfg)
  results <- lapply(cfg_list, function(cfg) select_global_winner(cfg, ns_cfg))

  result_tbl <- bind_rows(lapply(names(results), function(resp_name) {
    x <- results[[resp_name]]
    if (is.null(x)) {
      data.frame(
        Status = "MISSING_INPUTS",
        Response = resp_name,
        Model = NA_character_,
        Winner_CV = NA_character_,
        LOOIC = NA_real_,
        Prior_Scenario = ns_cfg$prior_scenario_target,
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
        Prior_Scenario = x$prior_scenario,
        Source_Path = x$source_path,
        Global_Path = x$global_path,
        stringsAsFactors = FALSE
      )
    }
  }))

  if (ns_cfg$run_namespace == "prior_sens_fullgrid") {
    out_csv <- file.path(
      ns_cfg$global_winner_dir,
      sprintf("GLOBAL_WINNER_SUMMARY_YEAR_RE_%s.csv", ns_cfg$scenario_tag_upper)
    )
  } else {
    out_csv <- file.path(ns_cfg$global_winner_dir, "GLOBAL_WINNER_SUMMARY_YEAR_RE.csv")
  }
  write.csv(result_tbl, out_csv, row.names = FALSE)

  n_ok <- sum(result_tbl$Status == "OK", na.rm = TRUE)
  cat(sprintf("Global winners generated: %d\n", n_ok))
  cat(sprintf("Summary: %s\n", out_csv))

  invisible(result_tbl)
}

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  run_namespace <- if (length(args) >= 1) args[[1]] else "canonical"
  prior_scenario_target <- if (length(args) >= 2) args[[2]] else "WeaklyInformative"
  build_global_winners_year_re(run_namespace = run_namespace, prior_scenario_target = prior_scenario_target)
}
