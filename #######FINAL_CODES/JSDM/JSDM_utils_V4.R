# ============================================================================ #
# SCRIPT: JSDM_utils_V4.R
# PURPOSE: Utility functions for JSDM Dirichlet hypothesis-testing pipeline
# VERSION: 4.0 (Corrected from V3)
# ============================================================================ #




# Required packages
libs <- c("dplyr", "tidyr", "brms", "loo", "bayesplot", "digest")
invisible(lapply(libs, library, character.only = TRUE))

# ============================================================================ #
# CONVERGENCE CRITERIA FOR JSDM DIRICHLET MODELS
# ============================================================================ #

CONVERGENCE_THRESHOLDS <- list(
  rhat_max = 1.01, # Gelman-Rubin statistic (ideal: < 1.01)
  ess_bulk_min = 400, # Effective sample size (bulk)
  ess_tail_min = 400, # Effective sample size (tail)
  divergent_max_pct = 1, # Maximum % of divergent transitions
  max_treedepth_pct = 5 # Maximum % hitting max_treedepth
)

# ============================================================================ #
# FUNCTION: check_convergence
# PURPOSE: Evaluate model convergence against defined criteria
# ============================================================================ #

check_convergence <- function(brms_fit, thresholds = CONVERGENCE_THRESHOLDS) {
  #' @param brms_fit A fitted brms model object
  #' @param thresholds Named list of convergence thresholds
  #' @return List with passed status, checks, and metrics

  # Extract diagnostics
  np <- nuts_params(brms_fit)
  rhat_vals <- rhat(brms_fit)
  neff <- neff_ratio(brms_fit)

  # Calculate metrics
  n_divergent <- sum(np$Value[np$Parameter == "divergent__"])
  n_chains <- length(unique(np$Chain))
  n_samples_per_chain <- nrow(np) / n_chains
  pct_divergent <- 100 * n_divergent / (n_samples_per_chain * n_chains)

  max_rhat <- max(rhat_vals, na.rm = TRUE)

  # Approximate absolute ESS from ratio
  post_draws <- posterior::as_draws_df(brms_fit$fit)
  n_total_draws <- nrow(post_draws)
  min_ess_bulk <- min(neff$Bulk_ESS, na.rm = TRUE) * n_total_draws / n_chains
  min_ess_tail <- min(neff$Tail_ESS, na.rm = TRUE) * n_total_draws / n_chains

  # Evaluate criteria
  checks <- list(
    rhat_ok = max_rhat < thresholds$rhat_max,
    ess_bulk_ok = min_ess_bulk > thresholds$ess_bulk_min,
    ess_tail_ok = min_ess_tail > thresholds$ess_tail_min,
    divergent_ok = pct_divergent < thresholds$divergent_max_pct
  )

  all_passed <- all(unlist(checks))

  # Report
  cat(sprintf("Convergence diagnostics:\n"))
  cat(sprintf("  - Max Rhat: %.4f (%s)\n", max_rhat, ifelse(checks$rhat_ok, "OK", "FAIL")))
  cat(sprintf("  - Min Bulk ESS: %.0f (%s)\n", min_ess_bulk, ifelse(checks$ess_bulk_ok, "OK", "FAIL")))
  cat(sprintf("  - Min Tail ESS: %.0f (%s)\n", min_ess_tail, ifelse(checks$ess_tail_ok, "OK", "FAIL")))
  cat(sprintf("  - Divergent: %.2f%% (%s)\n", pct_divergent, ifelse(checks$divergent_ok, "OK", "FAIL")))
  cat(sprintf("  Overall: %s\n", ifelse(all_passed, "CONVERGED", "ISSUES DETECTED")))

  return(list(
    passed = all_passed,
    checks = checks,
    metrics = list(
      max_rhat = max_rhat,
      min_ess_bulk = min_ess_bulk,
      min_ess_tail = min_ess_tail,
      pct_divergent = pct_divergent
    )
  ))
}

# ============================================================================ #
# FUNCTION: select_top_species_with_others
# PURPOSE: Deterministic species selection with "Others" aggregation
# CRITICAL FIX: Ensures compositional validity (proportions sum to 1)
# ============================================================================ #

select_top_species_with_others <- function(community_matrix, n_top = 10,
                                           focal_species = "mussismilia_hispida",
                                           min_prevalence = 0.10) {
  #' @param community_matrix Matrix with species as columns, observations as rows
  #' @param n_top Number of top species to select
  #' @param focal_species Name of focal species to force include
  #' @param min_prevalence Minimum prevalence (proportion of non-zero observations)
  #' @return List with matrix, species names, and percentage in "Others"

  # Calculate mean abundance and prevalence per species
  species_stats <- data.frame(
    species = colnames(community_matrix),
    mean_abundance = colMeans(community_matrix, na.rm = TRUE),
    prevalence = colMeans(community_matrix > 0, na.rm = TRUE)
  ) %>%
    dplyr::filter(prevalence >= min_prevalence) %>%
    dplyr::arrange(dplyr::desc(mean_abundance))

  # Select top N, ensuring focal species is included
  selected <- head(species_stats$species, n_top)

  if (!(focal_species %in% selected) && focal_species %in% species_stats$species) {
    # Replace last selected with focal species
    selected <- c(head(selected, n_top - 1), focal_species)
  }

  # --- CRITICAL: Calculate "Others" category ---
  excluded_species <- setdiff(colnames(community_matrix), selected)

  if (length(excluded_species) > 0) {
    others_abundance <- rowSums(community_matrix[, excluded_species, drop = FALSE], na.rm = TRUE)
    abundance_in_others <- sum(colMeans(community_matrix[, excluded_species, drop = FALSE], na.rm = TRUE))
    total_abundance <- sum(colMeans(community_matrix, na.rm = TRUE))
    pct_in_others <- 100 * abundance_in_others / total_abundance

    cat(sprintf("Species selection summary:\n"))
    cat(sprintf("  - Selected: %d species\n", length(selected)))
    cat(sprintf("  - Excluded: %d species -> aggregated as 'Others'\n", length(excluded_species)))
    cat(sprintf("  - Abundance in 'Others': %.1f%% of total\n", pct_in_others))

    # Build final matrix with Others
    selected_matrix <- community_matrix[, selected, drop = FALSE]
    final_matrix <- cbind(selected_matrix, others = others_abundance)

    # Normalize to ensure sum = 1 (compositional validity)
    final_matrix <- final_matrix / rowSums(final_matrix)

    # Verify normalization
    max_deviation <- max(abs(rowSums(final_matrix) - 1))
    if (max_deviation > 1e-10) {
      warning(sprintf("Normalization check: max deviation from 1.0 = %.2e", max_deviation))
    }

    return(list(
      matrix = final_matrix,
      species = c(selected, "others"),
      pct_in_others = pct_in_others,
      excluded_count = length(excluded_species)
    ))
  } else {
    # All species selected (rare case)
    final_matrix <- community_matrix[, selected, drop = FALSE] / rowSums(community_matrix[, selected, drop = FALSE])
    return(list(
      matrix = final_matrix,
      species = selected,
      pct_in_others = 0,
      excluded_count = 0
    ))
  }
}

# ============================================================================ #
# FUNCTION: prepare_dirichlet_data
# PURPOSE: Handle zeros for Dirichlet model with documented strategy
# ============================================================================ #

prepare_dirichlet_data <- function(composition_matrix, method = "half_min") {
  #' @param composition_matrix Normalized matrix with rows summing to 1
  #' @param method One of: "half_min", "uniform"
  #'   - half_min: Replace zeros with min(non-zero)/2 (V3 approach)
  #'   - uniform: Replace zeros with 1/(2*n_species)
  #' @return Matrix with zeros replaced, re-normalized to sum = 1

  n_zeros_original <- sum(composition_matrix == 0)
  n_total <- prod(dim(composition_matrix))
  pct_zeros <- 100 * n_zeros_original / n_total

  cat(sprintf("Zero handling:\n"))
  cat(sprintf("  - Original zeros: %d (%.1f%% of cells)\n", n_zeros_original, pct_zeros))
  cat(sprintf("  - Method: %s\n", method))

  if (method == "half_min") {
    min_non_zero <- min(composition_matrix[composition_matrix > 0], na.rm = TRUE)
    replacement_value <- min_non_zero / 2
    cat(sprintf("  - Replacement value: %.6f (half of min non-zero)\n", replacement_value))
  } else if (method == "uniform") {
    n_species <- ncol(composition_matrix)
    replacement_value <- 1 / (2 * n_species)
    cat(sprintf("  - Replacement value: %.6f (uniform prior)\n", replacement_value))
  } else {
    stop("Unknown method. Use 'half_min' or 'uniform'.")
  }

  # Apply replacement
  result <- composition_matrix
  result[result == 0] <- replacement_value

  # Re-normalize rows to sum = 1
  result <- result / rowSums(result)

  # Verify
  max_deviation <- max(abs(rowSums(result) - 1))
  if (max_deviation > 1e-10) {
    stop("Post-normalization check failed: rows do not sum to 1.")
  }

  cat("  Zero substitution complete. Rows normalized.\n")

  return(result)
}

# ============================================================================ #
# FUNCTION: cleanup_after_model
# PURPOSE: Memory management after model fitting
# ============================================================================ #

cleanup_after_model <- function(model_name) {
  #' @param model_name Name of the model for logging
  cat(sprintf("  - Cleaning memory after %s...\n", model_name))
  gc(verbose = FALSE)
}

# ============================================================================ #
# FUNCTION: analyze_loo_per_dataset
# PURPOSE: Compare LOO within each dataset (CORRECTED APPROACH)
# ============================================================================ #

analyze_loo_per_dataset <- function(loo_results_by_dataset) {
  #' @param loo_results_by_dataset Named list: dataset -> list of loo objects
  #' @return Summary table with winner per dataset

  results_summary <- data.frame(
    dataset = character(),
    winner = character(),
    delta_elpd_2nd = numeric(),
    n_models = integer(),
    stringsAsFactors = FALSE
  )

  for (dataset_name in names(loo_results_by_dataset)) {
    loo_list <- loo_results_by_dataset[[dataset_name]]

    if (length(loo_list) > 1) {
      comp <- loo_compare(loo_list)
      winner <- rownames(comp)[1]
      delta_2nd <- if (nrow(comp) > 1) abs(comp[2, "elpd_diff"]) else NA

      results_summary <- rbind(results_summary, data.frame(
        dataset = dataset_name,
        winner = winner,
        delta_elpd_2nd = delta_2nd,
        n_models = length(loo_list),
        stringsAsFactors = FALSE
      ))

      cat(sprintf("\n=== LOO Comparison: %s ===\n", dataset_name))
      print(comp)
    }
  }

  # Analyze consistency
  unique_winners <- unique(results_summary$winner)

  if (length(unique_winners) == 1) {
    cat(sprintf(
      "\n CONSISTENT WINNER: %s (same across all %d datasets)\n",
      unique_winners[1], nrow(results_summary)
    ))
    cat("   -> High confidence in model selection\n")
  } else {
    cat(sprintf("\n INCONSISTENT WINNERS: %s\n", paste(unique_winners, collapse = ", ")))
    cat("   -> Further investigation required (sensitivity to temporal scale)\n")

    most_common <- names(sort(table(results_summary$winner), decreasing = TRUE))[1]
    cat(sprintf("   -> Majority winner: %s\n", most_common))
  }

  return(results_summary)
}

# ============================================================================ #
# FUNCTION: verify_data_integrity
# PURPOSE: Check data quality before modeling
# ============================================================================ #

verify_data_integrity <- function(data_for_brms, organism_vars) {
  #' @param data_for_brms Prepared data frame for brms
  #' @param organism_vars Character vector of organism column names
  #' @return List of integrity check results

  cat("\n=== DATA INTEGRITY VERIFICATION ===\n")

  checks <- list()

  # Check 1: No NA values in response variables
  na_in_response <- sapply(data_for_brms[, organism_vars], function(x) sum(is.na(x)))
  checks$na_in_response <- all(na_in_response == 0)
  cat(sprintf("  [1] NA in response: %s\n", ifelse(checks$na_in_response, "PASS", "FAIL")))

  # Check 2: Proportions sum to approximately 1
  row_sums <- rowSums(data_for_brms[, organism_vars])
  checks$proportions_sum_ok <- max(abs(row_sums - 1)) < 0.01
  cat(sprintf(
    "  [2] Proportions sum to 1: %s (max deviation: %.4f)\n",
    ifelse(checks$proportions_sum_ok, "PASS", "FAIL"),
    max(abs(row_sums - 1))
  ))

  # Check 3: All proportions are positive
  checks$all_positive <- all(data_for_brms[, organism_vars] > 0)
  cat(sprintf("  [3] All proportions positive: %s\n", ifelse(checks$all_positive, "PASS", "FAIL")))

  # Check 4: Sufficient observations per site
  obs_per_site <- table(data_for_brms$SITE)
  checks$sufficient_obs_per_site <- all(obs_per_site >= 5)
  cat(sprintf(
    "  [4] Min obs per site >= 5: %s (min: %d)\n",
    ifelse(checks$sufficient_obs_per_site, "PASS", "FAIL"),
    min(obs_per_site)
  ))

  # Check 5: Sufficient sites for random effects
  checks$sufficient_sites <- length(unique(data_for_brms$SITE)) >= 5
  cat(sprintf(
    "  [5] Sufficient sites (>=5): %s (n=%d)\n",
    ifelse(checks$sufficient_sites, "PASS", "FAIL"),
    length(unique(data_for_brms$SITE))
  ))

  all_passed <- all(unlist(checks))

  cat(sprintf("\nOverall: %s\n", ifelse(all_passed, "ALL CHECKS PASSED", "SOME CHECKS FAILED")))

  return(list(
    passed = all_passed,
    checks = checks,
    n_observations = nrow(data_for_brms),
    n_sites = length(unique(data_for_brms$SITE)),
    n_organisms = length(organism_vars)
  ))
}

# ============================================================================ #
# FUNCTION: generate_cache_hash
# PURPOSE: Create deterministic hash for caching
# ============================================================================ #

generate_cache_hash <- function(data_for_brms, model_formula, prior_list) {
  #' @param data_for_brms Data frame for modeling
  #' @param model_formula Model formula object
  #' @param prior_list List of priors (or NULL)
  #' @return Character hash string

  data_hash <- digest::digest(data_for_brms)
  formula_hash <- digest::digest(deparse(model_formula))
  prior_hash <- if (is.null(prior_list)) "none" else digest::digest(prior_list)

  combined_hash <- digest::digest(paste(data_hash, formula_hash, prior_hash, sep = "_"))

  return(combined_hash)
}

cat("\n=== JSDM V4 UTILITIES LOADED ===\n")
cat("Available functions:\n")
cat("  - check_convergence()\n")
cat("  - select_top_species_with_others()\n")
cat("  - prepare_dirichlet_data()\n")
cat("  - cleanup_after_model()\n")
cat("  - analyze_loo_per_dataset()\n")
cat("  - verify_data_integrity()\n")
cat("  - generate_cache_hash()\n")
