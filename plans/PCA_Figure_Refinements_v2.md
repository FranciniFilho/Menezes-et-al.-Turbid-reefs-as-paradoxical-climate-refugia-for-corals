# Plan: PCA Figure Refinements (v2)

## Issues Identified

### Issue 1: Overlapping Legends in Magnitude Figures
- **Problem**: Reef, HAB, and Arc legends are stacked vertically and overlapping.
- **Solution**: Arrange Reef and HAB legends side-by-side horizontally, with Arc legend below.

### Issue 2: Text Sizes Need 20% Increase
- **Affected elements**:
  - Axis labels (X and Y): Currently `size = 11` → New: `size = 13.2` (~13)
  - Panel titles (e.g., "SST (°C)"): Currently `size = 12` → New: `size = 14.4` (~14)
  - Panel tags (a, b, c...): Currently `size = 12` → New: `size = 14.4` (~14)
  - Loading arrow labels: Currently `size = 3` → New: `size = 3.6` (~4)

### Issue 3: Variability Loadings Show CV IDs Instead of Variable Names
- **Problem**: In `create_loadings_panel`, the `loadings_labels` lookup fails for variability variables because they include CV suffixes (e.g., `sst_cv_2`, `dli_cv_30`).
- **Root Cause**: The `loadings_labels` dictionary only has generic keys like `sst_cv`, `dli_cv`, not scenario-specific ones.
- **Solution**: Expand `loadings_labels` to include all scenario-specific variable names OR strip the CV suffix before lookup.

---

## Proposed Changes

### 1. Refactor `create_legend_panel`

```r
# New layout: Reef | HAB side by side, Arc below
create_legend_panel <- function(reefs, habs) {
  # ... existing legend creation ...
  
  # Combine Reef and HAB horizontally
  top_row <- cowplot::plot_grid(leg_reef, leg_hab, ncol = 2, rel_widths = c(1, 1))
  
  # Stack with Arc below
  res <- cowplot::plot_grid(top_row, p_arch, ncol = 1, rel_heights = c(1.5, 1))
  
  return(res)
}
```

### 2. Update `theme_publication` for 20% Size Increase

```r
theme_publication <- theme_classic(base_size = 14.4) +  # Was 12
  theme(
    axis.title = element_text(size = 13, face = "bold"),  # Was 11
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Was 12
    plot.tag = element_text(face = "bold", size = 14),  # Was 12
    # ... rest unchanged ...
  )
```

### 3. Update Loading Labels in `create_loadings_panel`

```r
# Expand loadings_labels to include scenario-specific names
loadings_labels <- c(
  # Magnitude (existing)
  "log(sst_mean+1)" = "log(SST+1)",
  "log(mean_DLI_local+1)" = "log(DLI+1)",
  "log(chl_mean+1)" = "log(Chl+1)",
  "sqrt(DHW>4)" = "√(DHW>4)",
  
  # Variability CV_02
  "sst_cv_2" = "SST CV",
  "dli_cv_2" = "DLI CV",
  "chl_cv_2" = "Chl CV",
  
  # Variability CV_30
  "sst_cv_30" = "SST CV",
  "dli_cv_30" = "DLI CV",
  "chl_cv_30" = "Chl CV",
  
  # Variability CV_ALL
  "sst_cv_all" = "SST CV",
  "cv_DLI_local" = "DLI CV",
  "chl_cv_all" = "Chl CV"
)
```

Also increase label size in `geom_text`:
```r
geom_text(aes(...), size = 4, fontface = "bold")  # Was size = 3
```

---

## Files to Modify

| File | Changes |
|------|---------|
| [PCA_Visualization_Publication.R](file:///c:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes%20et%20al.%20Mus%20his%20distribution%20and%20abundance%20Abrolhos/######FINAL/#######FINAL_CODES/PCA_Visualization_Publication.R) | All changes above |

---

## Verification Checklist

- [ ] Magnitude figures: Legends side-by-side without overlap
- [ ] All text sizes increased by ~20%
- [ ] Variability figures: Loading labels show "SST CV", "DLI CV", "Chl CV" instead of "cv 2", "cv 30", etc.
- [ ] Panel tags (a, b, c...) visible and larger
- [ ] 12 figures regenerated successfully
