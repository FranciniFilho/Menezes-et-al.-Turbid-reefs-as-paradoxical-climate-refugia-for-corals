# Test only prior construction syntax
library(brms)

cat("Testing prior construction with extracted values...\n")

# Simulate UNIFIED_PRIORS_HYPOTHESIS
UNIFIED_PRIORS_HYPOTHESIS <- list(
  PC1_MAGNITUDE_scaled = c(-0.5, 0.5),
  PC2_MAGNITUDE_scaled = c(0.5, 0.5),
  PC1_VARIABILITY_scaled = c(-0.5, 0.5),
  PC2_VARIABILITY_scaled = c(-0.5, 0.5),
  DEPTH_M_scaled = c(0.2, 0.4),
  HABRR = c(-1.5, 0.4),
  HABTP = c(-0.8, 0.4)
)

# Extract values
pc1_mag_mean <- UNIFIED_PRIORS_HYPOTHESIS$PC1_MAGNITUDE_scaled[1]
pc1_mag_sd   <- UNIFIED_PRIORS_HYPOTHESIS$PC1_MAGNITUDE_scaled[2]
pc2_mag_mean <- UNIFIED_PRIORS_HYPOTHESIS$PC2_MAGNITUDE_scaled[1]
pc2_mag_sd   <- UNIFIED_PRIORS_HYPOTHESIS$PC2_MAGNITUDE_scaled[2]
pc1_var_mean <- UNIFIED_PRIORS_HYPOTHESIS$PC1_VARIABILITY_scaled[1]
pc1_var_sd   <- UNIFIED_PRIORS_HYPOTHESIS$PC1_VARIABILITY_scaled[2]
pc2_var_mean <- UNIFIED_PRIORS_HYPOTHESIS$PC2_VARIABILITY_scaled[1]
pc2_var_sd   <- UNIFIED_PRIORS_HYPOTHESIS$PC2_VARIABILITY_scaled[2]
depth_mean   <- UNIFIED_PRIORS_HYPOTHESIS$DEPTH_M_scaled[1]
depth_sd     <- UNIFIED_PRIORS_HYPOTHESIS$DEPTH_M_scaled[2]
habrr_mean   <- UNIFIED_PRIORS_HYPOTHESIS$HABRR[1]
habrr_sd     <- UNIFIED_PRIORS_HYPOTHESIS$HABRR[2]
habtp_mean   <- UNIFIED_PRIORS_HYPOTHESIS$HABTP[1]
habtp_sd     <- UNIFIED_PRIORS_HYPOTHESIS$HABTP[2]

# Build priors
priors <- c(
  prior(student_t(3, 0, 2.5), class = "Intercept"),
  prior(logistic(0, 1), class = "Intercept", dpar = "zoi"),
  prior(exponential(1), class = "sd"),
  prior(normal(pc1_mag_mean, pc1_mag_sd), class = "b", coef = "PC1_MAGNITUDE_scaled"),
  prior(normal(pc2_mag_mean, pc2_mag_sd), class = "b", coef = "PC2_MAGNITUDE_scaled"),
  prior(normal(pc1_var_mean, pc1_var_sd), class = "b", coef = "PC1_VARIABILITY_scaled"),
  prior(normal(pc2_var_mean, pc2_var_sd), class = "b", coef = "PC2_VARIABILITY_scaled"),
  prior(normal(depth_mean, depth_sd), class = "b", coef = "DEPTH_M_scaled"),
  prior(normal(habrr_mean, habrr_sd), class = "b", coef = "HABRR"),
  prior(normal(habtp_mean, habtp_sd), class = "b", coef = "HABTP")
)

cat("Successfully created priors:\n")
print(priors)
cat("\nTest passed: No syntax errors in prior construction.\n")