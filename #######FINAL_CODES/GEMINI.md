# Project Overview

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Marine ecology research project analyzing **Mussismilia hispida** (coral) distribution, abundance, health, and growth in the Abrolhos region, Brazil. Uses advanced Bayesian statistics (ZOIB models via brms) to understand coral responses to environmental gradients (temperature, light, chlorophyll, depth) across different habitats and reef architectures.

## Directory Structure

```
#######FINAL_CODES/          # All computational scripts (current location)
#######FINAL_RESULTS/        # Model outputs and analysis results
#######LAST ROUND MS/        # Scientific manuscript (Main text.docx, References.docx)
#######PHOTOS AND MAPS/      # Visualizations and figures
```

## Complete Analysis Workflow

### Step 1: Environmental Data Summarization (Python)
```bash
python PCA_LOCAL_bubbleplot.py
```
- Performs PCA on environmental variables (SST, DLI, CHL, DHW)
- Generates bubble plots for site-level environmental characterization
- Output: PCA scores and visualizations in `#####output_local_PCA_CV_*_FINAL/`

### Step 2: PERMANOVA Analysis (Python)
```bash
python PERMANOVA_Local_PCA.py
```
- Evaluates if REEF factor captures environmental variability
- **Justification for removing REEF as a categorical variable** from subsequent models
- If REEF explains significant environmental variation, it is excluded to avoid multicollinearity

### Step 3: Health & Growth Calculations (Python)
```bash
python Calculate_RGR_&_health_PCA.py
```
- Calculates Relative Growth Rate (RGR) from coral measurements
- Performs PCA on health metrics (proportion of healthy, dead and bleached coral tissues per colony) to derive HEALTH_PC1 and HEALTH_PC2
- Output: `resultados_biologicos_por_colonia.csv`

### Step 4: Data Integration (Python)
```bash
python FINAL_DATA_INTEGRATION.py
```
- Integrates biological (health/growth) + environmental (PCA) + frequency data per reef, site and habitat types (two habitats: PA - reef wall, RR - rocky reef, TP - reef top)
- Adds ARCH (inner/outer reef arcs) and lat/lon coordinates
- Output: `dados_abundancia_integrados_long_format.csv` in `#####output_local_PCA_CV_*_FINAL/`
- **Three CV scenarios (based on the Coefficient of Variability, measured as sliding window in python PCA_LOCAL_bubbleplot.py)**: CV_02, CV_30, CV_ALL

### Step 5: Modeling (R)

#### 5a. Bayesian Models (Primary Analysis)
```bash
Rscript Abundance_models_ZOIB_Bayes_DEPTH_&_HABITAT_REFACTORED.R
```
- Zero-One Inflated Beta (ZOIB) models for coral coverage (proportion [0,1])
- 3 datasets (CV_02, CV_30, CV_ALL) x 2 prior scenarios (WeaklyInformative, Informative)
- Multiple candidate models (null, linear, spline, ARCH interaction)
- Output: `C:/Users/rbfra/OneDrive/Bayesian_Analyses_ZOIB_Abundance_MULTI_CV/`

```bash
Rscript Health_&_Growth_models_Bayes_DEPTH_&_HABITAT_REFACTORED.R
```
- Gaussian models for HEALTH_PC1, HEALTH_PC2, RGR
- Same multi-dataset x multi-prior structure
- Output: `C:/Users/rbfra/OneDrive/Bayesian_Analyses_Health_Growth_MultiDataset/`

#### 5b. BRT Models (Complementary Interpretation)
```bash
Rscript Abundance_models_LM_BRT.R
Rscript Health_&_Growth_models_LM_BRT_NEW.R
```
- Boosted Regression Trees for abundance and health/growth
- Provides complementary interpretation to Bayesian models
- Useful for capturing non-linear relationships and variable importance

### Step 6: Visualization (R)

#### 6a. Master Visualization Pipeline
```bash
Rscript MASTER_Viz_Pipeline.R
```
- Auto-discovers winner models (files with "WINNER" in filename)
- Generates publication-quality forest plots (half-eye) for fixed effects
- Includes categorical variables (HAB levels shown as separate coefficients)
- Output: `C:/Users/rbfra/OneDrive/Bayesian_Figures_Publication/`

#### 6b. Bayesian-Specific Visualizations
```bash
Rscript Bayes_Viz_REFACTORED.R
```
- Posterior predictive checks
- Conditional effects plots
- Diagnostic plots (R-hat, ESS, trace plots)

#### 6c. BRT Visualizations
```bash
Rscript BRT_PDP_PLOTS.R
```
- Partial Dependence Plots for BRT models
- Shows marginal effects of predictors on response variables

## Key R Dependencies

```r
# Core Bayesian modeling
libs <- c("brms", "cmdstanr", "bayesplot", "tidybayes", "ggdist")

# BRT models
libs <- c("gbm", "dismo", "randomForest")

# Visualization
libs <- c("ggplot2", "patchwork", "cowplot")

# Data manipulation
libs <- c("dplyr", "tidyverse", "readxl")
```

**Critical**: Requires CmdStanR backend for brms. Set cores before running:
```r
options(mc.cores = 4)  # Adjust based on available CPU
set.seed(42)
```

## Python Dependencies

```python
# Data processing
pandas, numpy, xarray

# Analysis
sklearn (PCA), scipy

# Remote sensing downloads
requests, PyPDF2, tqdm
```

## Model Architecture

### ZOIB Models (Abundance)
- **Response**: `COVER_PROP` (proportion [0,1], includes zeros and ones)
- **Family**: `zero_one_inflated_beta()`
- **Predictors**: Scaled PCA components (PC1_MAGNITUDE, PC2_MAGNITUDE, PC1_VARIABILITY, PC2_VARIABILITY), DEPTH_M, HAB (habitat type), ARCH (architecture)
- **Random Effects**: `(1 | SITE)` for hierarchical structure
- **Priors**: Two scenarios - WeaklyInformative and Informative (dataset-specific)

### Gaussian Models (Health/Growth)
- **Response**: HEALTH_PC1, HEALTH_PC2, RGR (continuous)
- **Family**: `gaussian()`
- Similar predictor structure to ZOIB models

### Model Types
- `model_null`: Intercept only
- `model_linear_*`: Linear terms + optional depth/habitat
- `model_spline_*`: GAM splines (k=5) + optional depth/habitat
- `model_interaction_arch_*`: Splines with ARCH interaction

## Code Patterns

### Data Preparation Pattern (R)
All modeling scripts follow this pattern:
1. Load data with automatic delimiter detection (`read.csv2` fallback to `read.csv`)
2. Filter for `ORGANISMO == "MUSSISMILIA_HISPIDA"`
3. Convert coverage to proportion: `COVER_PROP = COBERTURA / 100`
4. Scale predictors: `scale()` for all continuous variables
5. Verify data integrity (range checks, NA handling)

### Model Fitting Pattern
```r
brm(
  formula = model_formula,
  data = prepared_data,
  family = zero_one_inflated_beta(),
  prior = priors,
  backend = "cmdstanr",
  cores = 4,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  control = list(adapt_delta = 0.95)
)
```

### Winner Model Selection
Models saved with `WINNER_` prefix based on LOOIC comparison via `loo_compare()`

## Important File Locations

### Input Data Paths (Hardcoded - Update if Moving Project)
```r
# PCA outputs for 3 CV scenarios
"../#######FINAL_RESULTS/#####output_local_PCA_CV_2_FINAL/dados_abundancia_integrados_long_format.csv"
"../#######FINAL_RESULTS/#####output_local_PCA_CV_30_FINAL/dados_abundancia_integrados_long_format.csv"
"../#######FINAL_RESULTS/#####output_local_PCA_CV_all_FINAL/dados_abundancia_integrados_long_format.csv"
```

### Output Directories
- **PCA outputs**: `../#######FINAL_RESULTS/#####output_local_PCA_CV_*_FINAL/`
- **ZOIB Abundance**: `C:/Users/rbfra/OneDrive/Bayesian_Analyses_ZOIB_Abundance_MULTI_CV/`
- **Health/Growth Bayesian**: `C:/Users/rbfra/OneDrive/Bayesian_Analyses_Health_Growth_MultiDataset/`
- **BRT Abundance**: `../#######FINAL_RESULTS/BRT_Abundance_CV_*/`
- **BRT Health/Growth**: `../#######FINAL_RESULTS/BRT_Health_Growth_CV_*/`
- **Figures**: `C:/Users/rbfra/OneDrive/Bayesian_Figures_Publication/`

## Common Tasks

### Testing Prior Syntax
Use `test_prior_syntax.R` to verify prior construction before running full models.

### Running Single Model
Edit the candidate_models loop to fit only the model of interest, or extract the specific formula and prior for manual fitting.

### Regenerating Figures
- **Bayesian models**: Run `MASTER_Viz_Pipeline.R` (auto-discovers winner models) and `Bayes_Viz_REFACTORED.R` (diagnostics and conditional effects)
- **BRT models**: Run `BRT_PDP_PLOTS.R` for partial dependence plots

### Adding New Predictor Variables
1. Add to `pca_predictors` or `predictor_lists` in data preparation function
2. Add scaled version to data preparation loop
3. Add prior specification in `dataset_specific_priors`
4. Update model formulas in `candidate_models`

## Environmental Variables (PCA Components)

- **PC1_MAGNITUDE/PC2_MAGNITUDE**: Mean environmental conditions (SST, DLI, CHL)
- **PC1_VARIABILITY/PC2_VARIABILITY**: Temporal variability (CV, frequency metrics)
- Derived from Lomb-Scargle periodogram analysis via `TIME_SERIES_LAG_temporal_environ.py`

## Habitat Types (HAB)
- **RR**: Reef rim
- **TP**: Tide pool
- (See sites metadata for full list)

## Architecture Types (ARCH)
- **Inner reef**: Protected environments
- **Outer reef**: Exposed environments

## Performance Notes

- Bayesian models can take **hours to days** depending on data size and model complexity
- Use `mc.cores` for parallel processing
- Monitor convergence via R-hat statistics (should be < 1.01)
- Check effective sample size (ESS) for reliable estimates
- Brms cache in output directories can be cleared for fresh runs: `unlink(brms_cache_dir, recursive = TRUE)`
