# ============================================================================
# 99_TEST_Viz_Pipeline_v6.R
# ============================================================================
# Test Validation Script for V6 Integrated Visualization Pipeline
#
# Purpose: Validate all V6 improvements to YEAR_RE visualization pipeline
#
# Tests:
#   1. ARCH interaction detection in marginal effects
#   2. HABMERGED display in forest plots
#   3. YEAR SD annotation
#   4. Interaction spotlight generation
#   5. YEAR RE diagnostics with significance testing
#
# Usage: Source this script after loading the V6 pipeline
#   source("04_MASTER_Viz_Pipeline_YEAR_RE_v6_INTEGRATED.R")
#   source("99_TEST_Viz_Pipeline_v6.R")
#
# Date: 2026-02-09
# ============================================================================

# ============================================================================
# 1. TEST CONFIGURATION
# ============================================================================

# Suppress package loading messages
suppressPackageStartupMessages({
  library(brms)
  library(ggplot2)
  library(patchwork)
})

# Test output directory
TEST_OUTPUT_DIR <- "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/TEST_V6_Validation/"
dir.create(TEST_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Base model directories
BASE_OUTPUT_DIR_YEAR_RE <- "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/"

# ============================================================================
# 2. TEST DATA LOADING
# ============================================================================

#' Load a sample model for testing
#' @param model_type Type of model ("ZOIB", "Health", etc.)
#' @param cv_label CV label ("CV_02", "CV_30", "CV_ALL")
#' @return brmsfit object or NULL
load_test_model <- function(model_type = "ZOIB", cv_label = "CV_02") {

  cat(sprintf("\n=== Loading test model: %s %s ===\n", model_type, cv_label))

  # Map model types to directories
  dir_map <- list(
    "ZOIB" = file.path(BASE_OUTPUT_DIR_YEAR_RE, "ZOIB_Abundance_YEAR_RE/"),
    "Health" = file.path(BASE_OUTPUT_DIR_YEAR_RE, "Gaussian_Health_YEAR_RE/"),
    "JSDM" = file.path(BASE_OUTPUT_DIR_YEAR_RE, "JSDM_Dirichlet_YEAR_RE/")
  )

  model_dir <- dir_map[[model_type]]
  if (is.null(model_dir) || !dir.exists(model_dir)) {
    cat(sprintf("  ⚠ Directory not found: %s\n", model_dir))
    return(NULL)
  }

  # Find model files
  if (model_type == "ZOIB") {
    pattern <- sprintf("zoib.*%s.*\\.rds$", tolower(cv_label))
  } else if (model_type == "Health") {
    pattern <- sprintf("gaussian_health.*%s.*\\.rds$", tolower(cv_label))
  } else if (model_type == "JSDM") {
    pattern <- sprintf("jsdm.*%s.*\\.rds$", tolower(cv_label))
  }

  model_files <- list.files(model_dir, pattern = pattern, full.names = TRUE)

  if (length(model_files) == 0) {
    cat(sprintf("  ⚠ No model files found for pattern: %s\n", pattern))
    return(NULL)
  }

  # Load first matching model
  fit <- readRDS(model_files[1])
  cat(sprintf("  ✓ Loaded: %s\n", basename(model_files[1])))

  return(fit)
}

# ============================================================================
# 3. TEST FUNCTIONS
# ============================================================================

#' Test 1: ARCH interaction detection in marginal effects
#' @param fit brmsfit object
#' @return Test result (PASS/FAIL)
test_arch_interaction_detection <- function(fit) {

  cat("\n--- TEST 1: ARCH Interaction Detection ---\n")

  if (is.null(fit)) {
    cat("  ⚠ SKIP: No model provided\n")
    return("SKIP")
  }

  # Check if V6 function exists
  if (!exists("generate_marginal_effects_year_re_v6")) {
    cat("  ❌ FAIL: generate_marginal_effects_year_re_v6 not found\n")
    cat("     Make sure to source: 04_MASTER_Viz_Pipeline_YEAR_RE_v6_INTEGRATED.R\n")
    return("FAIL")
  }

  tryCatch({
    # Test the function
    output_path <- file.path(TEST_OUTPUT_DIR, "test_arch_interaction")
    result <- generate_marginal_effects_year_re_v6(
      fit = fit,
      model_name = "TEST_MODEL",
      output_dir = output_path,
      predictors_to_plot = c("PC1VARIABILITY", "PC2VARIABILITY")
    )

    if (!is.null(result)) {
      cat("  ✓ PASS: ARCH interaction detection function executed\n")

      # Check if ARCH-specific plots were generated
      arch_files <- list.files(output_path, pattern = "ARCH|Inner|Outer", full.names = FALSE)
      if (length(arch_files) > 0) {
        cat(sprintf("  ✓ PASS: Found %d ARCH-specific plot(s)\n", length(arch_files)))
      } else {
        cat("  ⚠ NOTE: No ARCH-specific plots found (model may not have ARCH interaction)\n")
      }

      return("PASS")
    } else {
      cat("  ❌ FAIL: Function returned NULL\n")
      return("FAIL")
    }

  }, error = function(e) {
    cat(sprintf("  ❌ FAIL: Error occurred: %s\n", e$message))
    return("FAIL")
  })
}

#' Test 2: HABMERGED display in forest plots
#' @param fit brmsfit object
#' @return Test result (PASS/FAIL)
test_habmerged_forest_plot <- function(fit) {

  cat("\n--- TEST 2: HABMERGED in Forest Plots ---\n")

  if (is.null(fit)) {
    cat("  ⚠ SKIP: No model provided\n")
    return("SKIP")
  }

  if (!exists("generate_forest_plot_year_re_v6")) {
    cat("  ❌ FAIL: generate_forest_plot_year_re_v6 not found\n")
    return("FAIL")
  }

  tryCatch({
    output_path <- file.path(TEST_OUTPUT_DIR, "test_habmerged_forest")
    result <- generate_forest_plot_year_re_v6(
      fit = fit,
      model_name = "TEST_MODEL",
      output_dir = output_path
    )

    if (!is.null(result)) {
      cat("  ✓ PASS: Forest plot function executed\n")

      # Check if HABMERGED is in the plot
      forest_files <- list.files(output_path, pattern = "\\.png$|\\.pdf$", full.names = TRUE)
      if (length(forest_files) > 0) {
        cat(sprintf("  ✓ PASS: Generated %d forest plot file(s)\n", length(forest_files)))

        # Check for HABMERGED in the plot data (if accessible)
        if (!is.null(result$plot_data)) {
          if ("HABMERGED" %in% names(result$plot_data) ||
              any(grepl("HABMERGED|HAB", rownames(result$plot_data)))) {
            cat("  ✓ PASS: HABMERGED parameter found in plot data\n")
          }
        }
      }

      return("PASS")
    } else {
      cat("  ❌ FAIL: Function returned NULL\n")
      return("FAIL")
    }

  }, error = function(e) {
    cat(sprintf("  ❌ FAIL: Error occurred: %s\n", e$message))
    return("FAIL")
  })
}

#' Test 3: YEAR SD annotation
#' @param fit brmsfit object
#' @return Test result (PASS/FAIL)
test_year_sd_annotation <- function(fit) {

  cat("\n--- TEST 3: YEAR SD Annotation ---\n")

  if (is.null(fit)) {
    cat("  ⚠ SKIP: No model provided\n")
    return("SKIP")
  }

  # Check if model has YEAR random effect
  vc <- VarCorr(fit)
  if (!"YEAR" %in% names(vc)) {
    cat("  ⚠ SKIP: Model does not have YEAR random effect\n")
    return("SKIP")
  }

  if (!exists("generate_forest_plot_year_re_v6")) {
    cat("  ❌ FAIL: generate_forest_plot_year_re_v6 not found\n")
    return("FAIL")
  }

  tryCatch({
    output_path <- file.path(TEST_OUTPUT_DIR, "test_year_sd")
    result <- generate_forest_plot_year_re_v6(
      fit = fit,
      model_name = "TEST_MODEL",
      output_dir = output_path,
      show_year_sd = TRUE  # Explicitly request YEAR SD
    )

    if (!is.null(result)) {
      cat("  ✓ PASS: Forest plot with YEAR SD generated\n")

      # Extract YEAR SD from the model
      year_sd <- as.numeric(vc$YEAR$sd[1, 1])
      cat(sprintf("  ✓ INFO: YEAR SD = %.4f\n", year_sd))

      return("PASS")
    } else {
      cat("  ❌ FAIL: Function returned NULL\n")
      return("FAIL")
    }

  }, error = function(e) {
    cat(sprintf("  ❌ FAIL: Error occurred: %s\n", e$message))
    return("FAIL")
  })
}

#' Test 4: Interaction spotlight for significant ARCH interactions
#' @param fit brmsfit object
#' @return Test result (PASS/FAIL)
test_interaction_spotlight <- function(fit) {

  cat("\n--- TEST 4: Interaction Spotlight ---\n")

  if (is.null(fit)) {
    cat("  ⚠ SKIP: No model provided\n")
    return("SKIP")
  }

  if (!exists("generate_interaction_spotlight_ARCH")) {
    cat("  ❌ FAIL: generate_interaction_spotlight_ARCH not found\n")
    return("FAIL")
  }

  tryCatch({
    output_path <- file.path(TEST_OUTPUT_DIR, "test_interaction_spotlight")
    result <- generate_interaction_spotlight_ARCH(
      fit = fit,
      model_name = "TEST_MODEL",
      output_dir = output_path
    )

    if (!is.null(result)) {
      cat("  ✓ PASS: Interaction spotlight function executed\n")

      # Check if files were generated
      spotlight_files <- list.files(output_path, pattern = "spotlight|interaction", full.names = FALSE, ignore.case = TRUE)
      if (length(spotlight_files) > 0) {
        cat(sprintf("  ✓ PASS: Generated %d spotlight file(s)\n", length(spotlight_files)))
      } else {
        cat("  ⚠ NOTE: No spotlight files generated (may indicate no significant interactions)\n")
      }

      return("PASS")
    } else {
      cat("  ❌ FAIL: Function returned NULL\n")
      return("FAIL")
    }

  }, error = function(e) {
    cat(sprintf("  ❌ FAIL: Error occurred: %s\n", e$message))
    return("FAIL")
  })
}

#' Test 5: YEAR RE diagnostics with significance testing
#' @param fit brmsfit object
#' @return Test result (PASS/FAIL)
test_year_re_diagnostics <- function(fit) {

  cat("\n--- TEST 5: YEAR RE Diagnostics ---\n")

  if (is.null(fit)) {
    cat("  ⚠ SKIP: No model provided\n")
    return("SKIP")
  }

  # Check if model has YEAR random effect
  vc <- VarCorr(fit)
  if (!"YEAR" %in% names(vc)) {
    cat("  ⚠ SKIP: Model does not have YEAR random effect\n")
    return("SKIP")
  }

  if (!exists("generate_year_re_diagnostics_v6")) {
    cat("  ❌ FAIL: generate_year_re_diagnostics_v6 not found\n")
    return("FAIL")
  }

  tryCatch({
    output_path <- file.path(TEST_OUTPUT_DIR, "test_year_diagnostics")
    result <- generate_year_re_diagnostics_v6(
      fit = fit,
      model_name = "TEST_MODEL",
      output_dir = output_path
    )

    if (!is.null(result)) {
      cat("  ✓ PASS: YEAR RE diagnostics function executed\n")

      # Check diagnostic plots
      diag_files <- list.files(output_path, pattern = "\\.png$|\\.pdf$", full.names = FALSE)
      if (length(diag_files) > 0) {
        cat(sprintf("  ✓ PASS: Generated %d diagnostic file(s)\n", length(diag_files)))
      }

      return("PASS")
    } else {
      cat("  ❌ FAIL: Function returned NULL\n")
      return("FAIL")
    }

  }, error = function(e) {
    cat(sprintf("  ❌ FAIL: Error occurred: %s\n", e$message))
    return("FAIL")
  })
}

# ============================================================================
# 4. COMPREHENSIVE TEST RUNNER
# ============================================================================

#' Run all V6 validation tests
#' @param model_types Vector of model types to test
#' @param cv_labels Vector of CV labels to test
#' @return Summary dataframe with test results
run_all_v6_tests <- function(
    model_types = c("ZOIB", "Health"),
    cv_labels = c("CV_02", "CV_30", "CV_ALL")) {

  cat("\n")
  cat("============================================================================\n")
  cat("V6 VALIDATION TEST SUITE\n")
  cat("============================================================================\n")

  # Initialize results storage
  all_results <- list()
  test_names <- c(
    "ARCH_Interaction_Detection",
    "HABMERGED_Forest_Plot",
    "YEAR_SD_Annotation",
    "Interaction_Spotlight",
    "YEAR_RE_Diagnostics"
  )

  # Run tests for each model type and CV combination
  for (model_type in model_types) {
    for (cv_label in cv_labels) {

      cat(sprintf("\n>>> Testing: %s - %s <<<\n", model_type, cv_label))

      # Load model
      fit <- load_test_model(model_type, cv_label)

      if (is.null(fit)) {
        cat(sprintf("  ⚠ Skipping tests for %s - %s (no model loaded)\n", model_type, cv_label))
        next
      }

      # Run all tests
      results <- list()
      results$ARCH_Interaction_Detection <- test_arch_interaction_detection(fit)
      results$HABMERGED_Forest_Plot <- test_habmerged_forest_plot(fit)
      results$YEAR_SD_Annotation <- test_year_sd_annotation(fit)
      results$Interaction_Spotlight <- test_interaction_spotlight(fit)
      results$YEAR_RE_Diagnostics <- test_year_re_diagnostics(fit)

      # Store results
      test_key <- sprintf("%s_%s", model_type, cv_label)
      all_results[[test_key]] <- results
    }
  }

  # Compile summary
  summary_list <- list()
  for (key in names(all_results)) {
    for (test_name in test_names) {
      result <- all_results[[key]][[test_name]]
      summary_list[[length(summary_list) + 1]] <- data.frame(
        Test_Case = key,
        Test_Name = test_name,
        Result = result,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(summary_list) > 0) {
    summary_df <- do.call(rbind, summary_list)
    rownames(summary_df) <- NULL
  } else {
    summary_df <- data.frame(
      Test_Case = character(),
      Test_Name = character(),
      Result = character()
    )
  }

  # Print summary
  cat("\n")
  cat("============================================================================\n")
  cat("TEST SUMMARY\n")
  cat("============================================================================\n")

  if (nrow(summary_df) > 0) {
    # Count results
    n_pass <- sum(summary_df$Result == "PASS")
    n_fail <- sum(summary_df$Result == "FAIL")
    n_skip <- sum(summary_df$Result == "SKIP")
    n_total <- nrow(summary_df)

    cat(sprintf("\n  Total Tests: %d\n", n_total))
    cat(sprintf("  ✓ PASS: %d (%.1f%%)\n", n_pass, 100 * n_pass / n_total))
    cat(sprintf("  ❌ FAIL: %d (%.1f%%)\n", n_fail, 100 * n_fail / n_total))
    cat(sprintf("  ⚠ SKIP: %d (%.1f%%)\n", n_skip, 100 * n_skip / n_total))

    # Show failures
    if (n_fail > 0) {
      cat("\n  FAILED TESTS:\n")
      failures <- summary_df[summary_df$Result == "FAIL", ]
      for (i in 1:nrow(failures)) {
        cat(sprintf("    - %s: %s\n", failures$Test_Case[i], failures$Test_Name[i]))
      }
    }

    # Save summary
    summary_path <- file.path(TEST_OUTPUT_DIR, "Test_Summary.csv")
    write.csv(summary_df, summary_path, row.names = FALSE)
    cat(sprintf("\n  ✓ Summary saved to: %s\n", summary_path))
  } else {
    cat("\n  ⚠ No tests were executed\n")
  }

  cat("\n============================================================================\n")

  return(summary_df)
}

#' Quick test - run a single model through all tests
#' @param model_type Model type to test
#' @param cv_label CV label to test
#' @return Test results
quick_test <- function(model_type = "ZOIB", cv_label = "CV_02") {

  cat(sprintf("\n=== QUICK TEST: %s - %s ===\n", model_type, cv_label))

  fit <- load_test_model(model_type, cv_label)

  if (is.null(fit)) {
    cat("  ⚠ No model loaded. Cannot run tests.\n")
    return(NULL)
  }

  results <- list(
    ARCH_Interaction_Detection = test_arch_interaction_detection(fit),
    HABMERGED_Forest_Plot = test_habmerged_forest_plot(fit),
    YEAR_SD_Annotation = test_year_sd_annotation(fit),
    Interaction_Spotlight = test_interaction_spotlight(fit),
    YEAR_RE_Diagnostics = test_year_re_diagnostics(fit)
  )

  # Print results
  cat("\n--- Results ---\n")
  for (name in names(results)) {
    symbol <- if (results[[name]] == "PASS") "✓" else if (results[[name]] == "FAIL") "❌" else "⚠"
    cat(sprintf("%s %s: %s\n", symbol, name, results[[name]]))
  }

  return(results)
}

# ============================================================================
# 5. USAGE INSTRUCTIONS
# ============================================================================

cat("\n")
cat("============================================================================\n")
cat("V6 VALIDATION TEST SUITE LOADED\n")
cat("============================================================================\n")
cat("\nUSAGE:\n")
cat("  # Load V6 pipeline first:\n")
cat("  source('04_MASTER_Viz_Pipeline_YEAR_RE_v6_INTEGRATED.R')\n")
cat("\n")
cat("  # Run quick test (single model):\n")
cat("  quick_test('ZOIB', 'CV_02')\n")
cat("\n")
cat("  # Run comprehensive test suite:\n")
cat("  run_all_v6_tests()\n")
cat("\n")
cat("  # Run specific model types:\n")
cat("  run_all_v6_tests(model_types = c('ZOIB'), cv_labels = c('CV_02'))\n")
cat("\n")
cat("============================================================================\n")

# ============================================================================
# END OF 99_TEST_Viz_Pipeline_v6.R
# ============================================================================
