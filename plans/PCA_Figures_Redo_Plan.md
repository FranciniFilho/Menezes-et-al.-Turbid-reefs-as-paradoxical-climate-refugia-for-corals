# Plan: Refactoring PCA Bubble Plot Figures

## Overview
Create a standalone script to generate professional composite PCA figures from pre-computed results, eliminating the need to re-calculate environmental metrics from raw satellite data.

---

## Current State Analysis

### Existing Script: `PCA_LOCAL_bubbleplot.py`
- **Purpose**: Calculates environmental metrics from satellite data AND generates PCA plots
- **Input Paths**: Raw satellite data (SST, MODIS, DHW, etc.)
- **Output** (3 CV scenarios):
  - `#####output_local_PCA_CV_2_FINAL/`
  - `#####output_local_PCA_CV_30_FINAL/`
  - `#####output_local_PCA_CV_all_FINAL/`

### Available Output Files (per CV scenario)
| File | Description |
|------|-------------|
| `dados_consolidados_com_scores_das_duas_PCAs.xlsx` | Site data with PCA scores |
| `loadings_PCA_Magnitude_CV_*.csv` | Loadings for Magnitude PCA |
| `loadings_PCA_Variability_CV_*.csv` | Loadings for Variability PCA |
| `resumo_segmentos_por_site.xlsx` | Segment summary |

### Current Figure Outputs
1. **Magnitude**: 2 figures (bubble plots composite + PCA+loadings+importance)
2. **Variability**: 2 figures (bubble plots composite + PCA+loadings+importance)

### Magnitude PCA Variables
- `sst_mean` (log-transformed)
- `mean_DLI_local` (log-transformed)
- `chl_mean` (log-transformed)
- `prop_DHW_gt4` (sqrt-transformed)

### Variability PCA Variables
- `sst_cv_*` (CV_2, CV_30, or CV_all)
- `dli_cv_*`
- `chl_cv_*`

---

## Target Specification

### New Figure Requirements

#### Figure 1: Magnitude (5 panels, 2×3 grid)
| Panel | Content |
|-------|---------|
| (1,1) | Bubble plot: `sst_mean` |
| (1,2) | Bubble plot: `mean_DLI_local` |
| (1,3) | Bubble plot: `chl_mean` |
| (2,1) | Bubble plot: `prop_DHW_gt4` |
| (2,2) | **Loadings plot** |
| (2,3) | **Legend panel** (symbols, colors, arch styles) |

#### Figure 2: Variability (4 panels, 2×2 grid)
| Panel | Content |
|-------|---------|
| (1,1) | Bubble plot: `sst_cv_*` |
| (1,2) | Bubble plot: `dli_cv_*` |
| (2,1) | Bubble plot: `chl_cv_*` |
| (2,2) | **Loadings plot** |
| **Legend** | Positioned outside (right side) |

---

## Visual Style Reference

### From `Health_&_Growth_models_Bayes_DEPTH_&_HABITAT_REFACTORED.R`

```r
theme_publication <- theme_classic(base_size = 14) +
  theme(
    text = element_text(color = "black"),
    plot.title = element_text(face = "bold", size = 16)
  )
```

**Key style elements to replicate:**
- `theme_classic()` base
- Clean black text
- No grid lines (minimalist)
- Professional publication quality
- Bold titles
- Base font size: 14pt

### Current Bubble Plot Style (from Python script)
- **Colors**: `tab10` colormap (10 distinct colors for reefs)
- **Shapes** (habitats): `o`, `s`, `^`, `D`, `v`, `<`, `>`
- **Arch styles**:
  - `inner`: black edge, linewidth=2.0
  - `outer`: darkgrey edge, linewidth=0.75
- **Bubble size**: Scaled proportionally to variable values

---

## Implementation Language Decision: R vs Python

### Recommendation: **R**

#### Rationale for R:

| Factor | R | Python |
|--------|---|--------|
| **Style consistency** | ✅ Matches `theme_publication` from existing R scripts | ❌ Would need manual theme recreation |
| **ggplot2 features** | ✅ Excellent for faceted composite figures | ⚠️ Matplotlib requires more boilerplate |
| **Data integration** | ✅ Native Excel read (readxl) | ⚠️ Requires openpyxl |
| **Patchwork** | ✅ Easy panel composition with legends | ⚠️ More complex subplot management |
| **Publication quality** | ✅ Excellent vector graphics (PDF/SVG) | ⚠️ Good but requires more tuning |
| **Code reusability** | ✅ Can leverage existing R plotting functions | ❌ No existing Python plotting utilities in project |

#### R Stack to Use:
```r
libs <- c(
  "readxl",      # Read Excel files
  "ggplot2",     # Main plotting
  "patchwork",   # Composite figures
  "cowplot",     # Legend extraction
  "dplyr",       # Data manipulation
  "scales"       # Color scales
)
```

---

## Data Structure Required

### From `dados_consolidados_com_scores_das_duas_PCAs.xlsx`:
| Column | Description |
|--------|-------------|
| `Site_name` | Site identifier |
| `Reef_name` | Reef (for coloring) |
| `HAB` | Habitat type (for shape) |
| `Arch` | Inner/Outer arc (for border) |
| `PC1_Magnitude_CV_*` | PC1 score |
| `PC2_Magnitude_CV_*` | PC2 score |
| `PC1_Variability_CV_*` | PC1 score |
| `PC2_Variability_CV_*` | PC2 score |
| `sst_mean`, `chl_mean`, `mean_DLI_local`, `prop_DHW_gt4` | Magnitude vars |
| `sst_cv_*`, `dli_cv_*`, `chl_cv_*` | Variability vars |

### From `loadings_PCA_*.csv`:
```csv
,PC1,PC2
variable_name,loading1,loading2
```

---

## Implementation Plan

### Step 1: Input Configuration
```r
# Define paths to the 3 CV scenarios
scenarios <- list(
  CV_2  = "#####FINAL_RESULTS/#####output_local_PCA_CV_2_FINAL/",
  CV_30 = "#####FINAL_RESULTS/#####output_local_PCA_CV_30_FINAL/",
  CV_ALL = "#####FINAL_RESULTS/#####output_local_PCA_CV_all_FINAL/"
)
```

### Step 2: Data Loading Function
```r
load_pca_results <- function(scenario_path) {
  # Load Excel with scores
  df_scores <- read_excel(file.path(scenario_path, "dados_consolidados_com_scores_das_duas_PCAs.xlsx"))

  # Load loadings
  loadings_mag <- read_csv(file.path(scenario_path, "loadings_PCA_Magnitude_CV_*.csv"))
  loadings_var <- read_csv(file.path(scenario_path, "loadings_PCA_Variability_CV_*.csv"))

  # Return list
}
```

### Step 3: Bubble Plot Function
```r
create_bubble_plot <- function(df, pc1_col, pc2_col, size_var,
                               reef_colors, habitat_shapes, arch_styles) {
  # Create ggplot with:
  # - x = PC1, y = PC2
  # - color = Reef_name
  # - shape = HAB
  # - size = size_var
  # - aes(linewidth, color) for Arch border
}
```

### Step 4: Loadings Plot Function
```r
create_loadings_plot <- function(loadings_df, explained_var) {
  # Biplot-style arrows
  # Variables labeled with transformed names
  # Circle for unit variance
}
```

### Step 5: Legend Function
```r
create_legend_panel <- function(reef_colors, habitat_shapes, arch_styles) {
  # Create standalone legend ggplot
  # Three sub-legends: Reef (color), Habitat (shape), Arch (border)
}
```

### Step 6: Composite Figure Assembly
```r
# Magnitude (2×3 with legend in last cell)
p_mag <- (
  create_bubble_plot(sst)    | create_bubble_plot(dli) | create_bubble_plot(chl)
) / (
  create_bubble_plot(dhw)    | create_loadings_plot()  | create_legend_panel()
)

# Variability (2×2, legend external)
p_var <- (
  create_bubble_plot(sst_cv) | create_bubble_plot(dli_cv)
) / (
  create_bubble_plot(chl_cv) | create_loadings_plot()
)
p_var_with_legend <- p_var | plot_layout(guides = "collect")
```

### Step 7: Output
```r
# Save figures
ggsave("PCA_Magnitude_Composite.png", p_mag, width = 14, height = 10, dpi = 300)
ggsave("PCA_Variability_Composite.png", p_var_with_legend, width = 10, height = 10, dpi = 300)
```

---

## Detailed Visual Specifications

### Color Scheme (Reefs)
```r
reef_palette <- c(
  "Abrolhos" = "#D55E00",
  "Parcel das Paredes" = "#0072B2",
  "Sebastiao Gomes" = "#009E73",
  "Santa Barbara" = "#CC79A7",
  "California" = "#F0E442",
  # ... use tab10 equivalent from RColorBrewer or viridis
)
```

### Shape Scheme (Habitats)
```r
habitat_shapes <- c(
  "PA" = 16,  # Circle
  "RR" = 15,  # Square
  "TP" = 17   # Triangle
)
```

### Border Scheme (Arch)
```r
arch_linetypes <- list(
  "inner"  = list(color = "black", size = 1.5),
  "outer"  = list(color = "grey50", size = 0.5)
)
```

### Axis Labels
- X: "PC1 (XX.X%)"
- Y: "PC2 (YY.Y%)"
- Percentage from explained variance

### Titles
- **NO main figure title** (as requested)
- Panel subtitles optional but not required

### English Labels for Variables
| Original | English Label |
|----------|---------------|
| `sst_mean` | SST Mean |
| `mean_DLI_local` | Mean DLI |
| `chl_mean` | Chlorophyll Mean |
| `prop_DHW_gt4` | DHW > 4°C |
| `sst_cv_*` | SST CV |
| `dli_cv_*` | DLI CV |
| `chl_cv_*` | Chl CV |

---

## File Structure

### New Script: `PCA_Visualization_Only.R`
```
#######FINAL_CODES/
├── PCA_Visualization_Only.R          # NEW standalone script
├── PCA_LOCAL_bubbleplot.py           # Original (keep unchanged)
└── ...
```

### Output Files (generated per CV scenario)
```
#######FINAL_RESULTS/#####output_local_PCA_CV_*/FINAL/
├── PCA_Magnitude_Composite.png
├── PCA_Magnitude_Composite.pdf
├── PCA_Variability_Composite.png
├── PCA_Variability_Composite.pdf
└── PCA_Figures_Log.txt
```

---

## Validation Checklist

- [ ] All 3 CV scenarios (CV_2, CV_30, CV_all) produce figures
- [ ] Magnitude figure has 5 panels (2×3 grid)
- [ ] Variability figure has 4 panels (2×2 grid)
- [ ] Loadings arrows correctly oriented
- [ ] Bubble sizes proportional to variable values
- [ ] Colors match reefs consistently across panels
- [ ] Shapes match habitats
- [ ] Arch borders (inner/outer) clearly distinct
- [ ] All labels in English
- [ ] No main figure title
- [ ] Style matches `theme_publication` from R script
- [ ] PC1/PC2 axis labels show explained variance percentage

---

## Open Questions for User

1. **CV Scenario Selection**: Should the script generate figures for all 3 scenarios automatically, or should there be a selector parameter?

2. **Output Format**: PNG only, or also PDF/SVG for publication?

3. **Variable Labels**: Confirm the exact English labels for the transformed variables (e.g., "log(SST+1)" vs "SST Mean")?

4. **Color Palette**: Should we use the same `tab10` colors from Python, or switch to a colorblind-friendly palette?

5. **Legend Position**: For the Variability figure, should the legend be positioned to the right (wider figure) or bottom (taller figure)?

---

## Next Steps (Upon User Approval)

1. Create `#######FINAL_CODES/PCA_Visualization_Only.R`
2. Implement data loading functions
3. Implement bubble plot function with proper styling
4. Implement loadings plot function
5. Implement legend panel function
6. Implement composite figure assembly with patchwork
7. Test on one CV scenario
8. Validate output against current figures
9. Generate all 3 CV scenarios
10. Document usage

---

*Document created: 2026-01-12*
*Project: Mussismilia hispida distribution and abundance - Abrolhos*
