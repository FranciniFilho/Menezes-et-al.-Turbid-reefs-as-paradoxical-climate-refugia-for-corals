# ============================================================================
# 00_LOO_Selection_Functions.R (minimal stub)
# ============================================================================
# Original file was removed; this stub provides the minimal interface
# expected by 00_LOO_Selection_Functions_YEAR_RE.R.
# All substantive functions are defined in the YEAR_RE companion file.
# ============================================================================

suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(loo)
  library(dplyr)
  library(tidyr)
})

options(mc.cores = parallel::detectCores())
set.seed(42)
