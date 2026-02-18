# ============================================================================
# 00_RUN_ALL_YEAR_RE_MODELS.R
# ============================================================================
# Master script to run all YEAR RE models with checkpointing
# Order of execution: ZOIB -> RGR -> Health -> JSDM
# ============================================================================
# Total models: 60
# - ZOIB: 12 (3 CV × 4 combos)
# - RGR: 12 (3 CV × 4 combos)
# - Health: 24 (3 CV × 4 combos × 2 responses)
# - JSDM: 12 (3 CV × 4 combos)
# ============================================================================

# Set working directory
setwd("C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_CODES/Bayes_models/New_Bayes_Models_YEAR_RE/")

# ============================================================================
# CONFIGURATION
# ============================================================================

# Scripts to run (in order)
scripts_to_run <- list(
  list(
    name = "ZOIB Abundance (YEAR RE)",
    script = "01_ZOIB_YEAR_RE_Full_LOO_Selection.R",
    n_models = 12,
    description = "ZOIB models with (1|REEF) + (1|YEAR)"
  ),
  list(
    name = "Gaussian RGR (NO YEAR RE)",
    script = "02a_Gaussian_RGR_Full_LOO_Selection.R",
    n_models = 12,
    description = "RGR models with (1|REEF) only (time-integrated metric)"
  ),
  list(
    name = "Gaussian Health (YEAR RE)",
    script = "02b_Gaussian_Health_YEAR_RE_Full_LOO_Selection.R",
    n_models = 24,
    description = "Health models with (1|REEF) + (1|YEAR)"
  ),
  list(
    name = "JSDM Dirichlet (YEAR RE)",
    script = "03_JSDM_YEAR_RE_Full_LOO_Selection.R",
    n_models = 12,
    description = "JSDM models with (1|REEF) + (1|YEAR)"
  )
)

# Total models
total_models <- sum(sapply(scripts_to_run, function(x) x$n_models))

# ============================================================================
# EXECUTION
# ============================================================================

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════════╗\n")
cat("║         MASTER SCRIPT: ALL YEAR RE MODELS                         ║\n")
cat("╚═══════════════════════════════════════════════════════════════════╝\n")
cat("\n")
cat(sprintf("Total scripts to run: %d\n", length(scripts_to_run)))
cat(sprintf("Total models to fit: %d\n", total_models))
cat(sprintf("Output directory: %s\n", "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/"))
cat("\n")

# Track results
results <- data.frame(
  Script = character(),
  Models = integer(),
  Status = character(),
  Time_Mins = numeric(),
  stringsAsFactors = FALSE
)

# Run each script
for (i in seq_along(scripts_to_run)) {
  script_info <- scripts_to_run[[i]]

  cat(sprintf("\n"))
  cat("════════════════════════════════════════════════════════════════\n")
  cat(sprintf("  SCRIPT %d of %d: %s\n", i, length(scripts_to_run), script_info$name))
  cat(sprintf("  %s\n", script_info$description))
  cat(sprintf("  Expected models: %d\n", script_info$n_models))
  cat("════════════════════════════════════════════════════════════════\n")
  cat("\n")

  start_time <- Sys.time()

  # Run the script
  status <- tryCatch({
    source(script_info$script)
    "SUCCESS"
  }, error = function(e) {
    cat(sprintf("\n❌ Script failed: %s\n", e$message))
    "FAILED"
  })

  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")

  # Store results
  results <- rbind(results, data.frame(
    Script = script_info$name,
    Models = script_info$n_models,
    Status = status,
    Time_Mins = as.numeric(elapsed),
    stringsAsFactors = FALSE
  ))

  cat(sprintf("\n✓ Script %s completed in %.1f minutes\n", script_info$name, as.numeric(elapsed)))
}

# ============================================================================
# FINAL SUMMARY
# ============================================================================

cat("\n")
cat("╔═══════════════════════════════════════════════════════════════════╗\n")
cat("║                    ALL SCRIPTS COMPLETED                          ║\n")
cat("╚═══════════════════════════════════════════════════════════════════╝\n")
cat("\n")

cat("Results Summary:\n")
print(results)

cat(sprintf("\nTotal time: %.1f minutes (%.2f hours)\n",
            sum(results$Time_Mins), sum(results$Time_Mins) / 60))

# Count successes
n_success <- sum(results$Status == "SUCCESS")
n_failed <- sum(results$Status == "FAILED")

cat(sprintf("\nScripts succeeded: %d/%d\n", n_success, nrow(results)))
cat(sprintf("Scripts failed: %d/%d\n", n_failed, nrow(results)))

if (n_success > 0) {
  cat("\n✓ All scripts completed successfully!\n")
  cat("\nTo view results, check:\n")
  cat("  - C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/\n")
  cat("\nWinner models are saved as WINNER_*.rds files\n")
  cat("\nRunning global winner consolidation...\n")
  tryCatch({
    source("05_CONSOLIDATE_GLOBAL_WINNERS_YEAR_RE.R")
    global_summary <- build_global_winners_year_re()
    n_global <- sum(global_summary$Status == "OK", na.rm = TRUE)
    cat(sprintf("Global winners created: %d\n", n_global))
    cat("Global directory: C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/Global_Winners_YEAR_RE/\n")
  }, error = function(e) {
    cat(sprintf("\n⚠ Global winner consolidation failed: %s\n", e$message))
  })
} else {
  cat("\n⚠ Some scripts failed. Check error messages above.\n")
}
