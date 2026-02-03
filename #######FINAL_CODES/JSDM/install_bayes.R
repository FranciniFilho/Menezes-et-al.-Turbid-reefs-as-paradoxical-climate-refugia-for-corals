
# Install Bayesian packages
pkgs <- c("brms", "loo", "cmdstanr", "bayesplot", "tidybayes", "ggplot2", "patchwork", "ggdist")
missing <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]

if(length(missing) > 0) {
  message("Installing: ", paste(missing, collapse=", "))
  install.packages(missing, repos = "http://cran.us.r-project.org")
} else {
  message("All packages installed.")
}

# Check cmdstanr
library(cmdstanr)
tryCatch({
  check_cmdstan_toolchain(fix = TRUE)
  # We might need to install cmdstan itself if not present, but let's check first
  # install_cmdstan() 
}, error = function(e) {
  message("CmdStan check warning: ", e$message)
})
