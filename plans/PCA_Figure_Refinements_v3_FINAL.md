# Plan: PCA Figure Refinements (v3 FINAL)

> **Status**: Plano definitivo para implementação
> **Baseado em**: v2 + análise de lacunas
> **Data**: 2026-01-12

---

## Issues Identified

### Issue 1: Overlapping Legends in Magnitude Figures
- **Problem**: Reef, HAB, and Arc legends are stacked vertically and overlapping.
- **Root Cause**: Current layout uses vertical stacking with `wrap_elements()` which compresses items.
- **Solution**: Arrange Reef and HAB legends side-by-side horizontally, with Arc legend below.

### Issue 2: Text Sizes Need 20% Increase
- **Affected elements**:
  - Base font size: Currently `12` → New: `14.4` (~14.5)
  - Axis labels (X and Y): Currently `size = 11` → New: `size = 13.2` (~13)
  - Panel titles (e.g., "SST (°C)"): Currently `size = 12` → New: `size = 14.4` (~14)
  - Panel tags (a, b, c...): Currently `size = 12` → New: `size = 14.4` (~14)
  - Loading arrow labels: Currently `size = 3` → New: `size = 3.6` (~4)
  - Legend text: Will inherit from base_size increase

### Issue 3: Variability Loadings Show CV IDs Instead of Variable Names
- **Problem**: In `create_loadings_panel`, the `loadings_labels` lookup fails for variability variables because they include CV suffixes (e.g., `sst_cv_2`, `dli_cv_30`, `cv_DLI_local`).
- **Root Cause**: The `loadings_labels` dictionary only has generic keys like `sst_cv`, `dli_cv`, not scenario-specific ones.
- **Solution**: Expand `loadings_labels` to include all scenario-specific variable names.

---

## Verification: Actual Variable Names in Loadings CSVs

| Scenario | Variable Names (as they appear) |
|----------|----------------------------------|
| CV_2 | `sst_cv_2`, `dli_cv_2`, `chl_cv_2` |
| CV_30 | `sst_cv_30`, `dli_cv_30`, `chl_cv_30` |
| CV_ALL | `sst_cv_all`, `cv_DLI_local`, `chl_cv_all` |

**Conclusion**: All nine (9) scenario-specific combinations must be mapped.

---

## Proposed Changes

### Change 1: Refactor `create_legend_panel` for Horizontal Layout

**Current problematic pattern:**
```r
# Vertical stacking causes overlap
res <- wrap_elements(full = leg_reef) /
    wrap_elements(full = leg_hab) /
    wrap_elements(full = p_arch)
```

**New layout:**
```r
create_legend_panel <- function(reefs, habs) {
    # ... [criação de p_reef, p_hab, p_arch como antes] ...

    # Função interna robusta de extração
    get_leg <- function(p) {
        leg <- cowplot::get_plot_component(p, "guide-box", return_all = TRUE)
        if (is.list(leg) && length(leg) > 0) leg <- leg[[1]]
        if (is.null(leg)) leg <- cowplot::get_legend(p)
        return(leg)
    }

    cat("      - Extraindo legendas Reef e HAB...\n")
    leg_reef <- get_leg(p_reef)
    leg_hab <- get_leg(p_hab)

    # NOVO: Dispor Reef e HAB lado a lado (horizontalmente)
    cat("      - Combinando Reef | HAB (lado a lado)...\n")
    top_row <- cowplot::plot_grid(
        leg_reef,
        leg_hab,
        ncol = 2,                    # 2 colunas = lado a lado
        rel_widths = c(1, 1),        # Larguras iguais
        align = "v"                  # Alinhar verticalmente
    )

    # NOVO: Arc logo abaixo do par Reef|HAB
    cat("      - Posicionando Arc abaixo...\n")
    res <- cowplot::plot_grid(
        top_row,                      # Reef | HAB (linha superior)
        p_arch,                       # Arc (linha inferior)
        ncol = 1,                     # 1 coluna = empilhar verticalmente
        rel_heights = c(1.5, 1),      # Top row 50% mais alto
        align = "h"                   # Alinhar horizontalmente
    )

    return(res)
}
```

**Visual result:**
```
┌─────────────────────┐
│  Reef    |    HAB   │  ← top_row (ncol=2)
├─────────────────────┤
│       Arc           │  ← p_arch (abaixo)
└─────────────────────┘
```

---

### Change 2: Update `theme_publication` for 20% Size Increase

**Location:** Lines 63-78 in current script

**Current:**
```r
theme_publication <- theme_classic(base_size = 12) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major = element_line(color = "grey90", linetype = "dashed", linewidth = 0.3),
    strip.background = element_blank(),
    plot.tag = element_text(face = "bold", size = 12)
  )
```

**New (20% increase applied):**
```r
theme_publication <- theme_classic(base_size = 14.4) +  # 12 * 1.2 = 14.4
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(size = 12),           # 10 * 1.2 = 12
    axis.title = element_text(size = 13, face = "bold"),  # 11 * 1.2 = 13.2 → 13
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # 12 * 1.2 = 14.4 → 14
    legend.title = element_text(size = 12, face = "bold"),  # 10 * 1.2 = 12
    legend.text = element_text(size = 11),          # 9 * 1.2 = 10.8 → 11
    panel.grid.major = element_line(color = "grey90", linetype = "dashed", linewidth = 0.3),
    strip.background = element_blank(),
    plot.tag = element_text(face = "bold", size = 14)  # 12 * 1.2 = 14.4 → 14
  )
```

**Side effects of base_size increase:**
- All text elements inheriting from base_size will automatically increase
- This is DESIRED for consistent appearance

---

### Change 3: Expand `loadings_labels` for All Scenarios

**Location:** Lines 51-60 in current script

**Current (insufficient):**
```r
loadings_labels <- c(
  "log(sst_mean+1)" = "log(SST+1)",
  "log(mean_DLI_local+1)" = "log(DLI+1)",
  "log(chl_mean+1)" = "log(Chl+1)",
  "sqrt(DHW>4)" = "√(DHW>4)",
  "sst_cv" = "SST CV",
  "dli_cv" = "DLI CV",
  "chl_cv" = "Chl CV",
  "cv_DLI_local" = "DLI CV"
)
```

**New (complete coverage):**
```r
loadings_labels <- c(
  # ===== Magnitude (inalterado) =====
  "log(sst_mean+1)" = "log(SST+1)",
  "log(mean_DLI_local+1)" = "log(DLI+1)",
  "log(chl_mean+1)" = "log(Chl+1)",
  "sqrt(DHW>4)" = "√(DHW>4)",

  # ===== Variability - CV_02 =====
  "sst_cv_2" = "SST CV",
  "dli_cv_2" = "DLI CV",
  "chl_cv_2" = "Chl CV",

  # ===== Variability - CV_30 =====
  "sst_cv_30" = "SST CV",
  "dli_cv_30" = "DLI CV",
  "chl_cv_30" = "Chl CV",

  # ===== Variability - CV_ALL =====
  "sst_cv_all" = "SST CV",
  "cv_DLI_local" = "DLI CV",
  "chl_cv_all" = "Chl CV"
)
```

---

### Change 4: Increase Loading Label Size in `create_loadings_panel`

**Location:** Line 219 in current script

**Current:**
```r
geom_text(aes(x = xend * 1.2, y = yend * 1.2, label = label), size = 3, fontface = "bold")
```

**New:**
```r
geom_text(aes(x = xend * 1.2, y = yend * 1.2, label = label), size = 4, fontface = "bold")
```

**Calculation:** 3 * 1.2 = 3.6 → round to 4

---

### Change 5: (OPTIONAL) Increase Legend Point Size for Consistency

**Location:** Lines 229, 236 in current script

**If legend points look too small after text increase:**

**Current:**
```r
geom_point(shape = 21, size = 3.5)
```

**New (20% increase):**
```r
geom_point(shape = 21, size = 4.2)  # 3.5 * 1.2 = 4.2
```

**Note:** This is OPTIONAL - evaluate visually after other changes.

---

## Summary of All Changes

| Line # | Element | Current | New | Change |
|---------|---------|---------|-----|--------|
| 63 | `base_size` | 12 | 14.4 | +20% |
| 67 | `axis.text` | 10 | 12 | +20% |
| 68 | `axis.title` | 11 | 13 | +20% |
| 69 | `plot.title` | 12 | 14 | +20% |
| 70 | `legend.title` | 10 | 12 | +20% |
| 71 | `legend.text` | 9 | 11 | +20% |
| 77 | `plot.tag` | 12 | 14 | +20% |
| 51-60 | `loadings_labels` | 8 entries | 17 entries | +9 mappings |
| 219 | `geom_text` size | 3 | 4 | +33% |
| 226-282 | `create_legend_panel` | Vertical stack | Horizontal+below | Layout refactor |
| 229/236 | legend `size` | 3.5 | 4.2 (optional) | +20% |

---

## Implementation Order

1. **Update `theme_publication`** (Lines 63-78) - Affects all text globally
2. **Expand `loadings_labels`** (Lines 51-60) - Fixes variability labels
3. **Update `geom_text` size** (Line 219) - Increases loading label size
4. **Refactor `create_legend_panel`** (Lines 226-282) - Fixes overlapping legends
5. **OPTIONAL: Update legend point sizes** (Lines 229, 236) - If needed visually

---

## Code Diff (Ready to Apply)

```diff
--- a/PCA_Visualization_Publication.R
+++ b/PCA_Visualization_Publication.R
@@ -48,19 +48,28 @@
 var_labels_var_default <- c("SST CV (%)", "DLI CV (%)", "Chl-a CV (%)")

 # Labels para loadings (abreviados e limpos)
 loadings_labels <- c(
-    "log(sst_mean+1)" = "log(SST+1)",
-    "log(mean_DLI_local+1)" = "log(DLI+1)",
-    "log(chl_mean+1)" = "log(Chl+1)",
-    "sqrt(DHW>4)" = "√(DHW>4)",
-    "sst_cv" = "SST CV",
-    "dli_cv" = "DLI CV",
-    "chl_cv" = "Chl CV",
-    "cv_DLI_local" = "DLI CV"
+    # Magnitude
+    "log(sst_mean+1)" = "log(SST+1)",
+    "log(mean_DLI_local+1)" = "log(DLI+1)",
+    "log(chl_mean+1)" = "log(Chl+1)",
+    "sqrt(DHW>4)" = "√(DHW>4)",
+    # Variability - CV_02
+    "sst_cv_2" = "SST CV",
+    "dli_cv_2" = "DLI CV",
+    "chl_cv_2" = "Chl CV",
+    # Variability - CV_30
+    "sst_cv_30" = "SST CV",
+    "dli_cv_30" = "DLI CV",
+    "chl_cv_30" = "Chl CV",
+    # Variability - CV_ALL
+    "sst_cv_all" = "SST CV",
+    "cv_DLI_local" = "DLI CV",
+    "chl_cv_all" = "Chl CV"
 )

 # Tema profissional
 theme_publication <- theme_classic(base_size = 14.4) +
     theme(
         text = element_text(color = "black"),
-        axis.text = element_text(size = 10),
-        axis.title = element_text(size = 11, face = "bold"),
-        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
-        legend.title = element_text(size = 10, face = "bold"),
-        legend.text = element_text(size = 9),
+        axis.text = element_text(size = 12),
+        axis.title = element_text(size = 13, face = "bold"),
+        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
+        legend.title = element_text(size = 12, face = "bold"),
+        legend.text = element_text(size = 11),
         panel.grid.major = element_line(
             color = "grey90",
             linetype = "dashed",
             linewidth = 0.3
         ),
         strip.background = element_blank(),
-        plot.tag = element_text(face = "bold", size = 12)
+        plot.tag = element_text(face = "bold", size = 14)
     )

@@ -216,7 +225,7 @@
         geom_segment(aes(x = 0, y = 0, xend = xend, yend = yend),
             arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
             color = "#d62728", linewidth = 0.7
         ) +
-        geom_text(aes(x = xend * 1.2, y = yend * 1.2, label = label), size = 3, fontface = "bold") +
+        geom_text(aes(x = xend * 1.2, y = yend * 1.2, label = label), size = 4, fontface = "bold") +
         labs(title = "Loadings", x = sprintf("PC1 (%.1f%%)", exp_var[1]), y = sprintf("PC2 (%.1f%%)", exp_var[2])) +
         theme_publication +
         coord_fixed(ratio = 1, xlim = c(-2.2, 2.2), ylim = c(-2.2, 2.2))
@@ -226,31 +235,33 @@

 # Painel de Legenda (Refatorado para patchwork)
 create_legend_panel <- function(reefs, habs) {
     # Legenda de recifes (cores)
-    p_reef <- ggplot(data.frame(Reef = factor(reefs, levels = reefs)), aes(x = 1, y = seq_along(Reef), fill = Reef)) +
-        geom_point(shape = 21, size = 3.5) +
+    p_reef <- ggplot(data.frame(Reef = factor(reefs, levels = reefs)), aes(x = 1, y = seq_along(Reef), fill = Reef)) +
+        geom_point(shape = 21, size = 4.2) +
         scale_fill_manual(values = reef_colors, name = "Reef") +
         theme_publication +
         theme(legend.position = "right")

     # Legenda de habitats (formas)
-    p_hab <- ggplot(data.frame(HAB = factor(habs, levels = habs)), aes(x = 1, y = seq_along(HAB), shape = HAB)) +
-        geom_point(size = 3.5, fill = "grey60") +
+    p_hab <- ggplot(data.frame(HAB = factor(habs, levels = habs)), aes(x = 1, y = seq_along(HAB), shape = HAB)) +
+        geom_point(size = 4.2, fill = "grey60") +
         scale_shape_manual(values = habitat_shapes, name = "Habitat") +
         theme_publication +
         theme(legend.position = "right")

     # Legenda manual para Arch
-    p_arch <- ggplot(data.frame(x = 1, y = 2:1, label = c("Inner Arc", "Outer Arc"), s = c(1.2, 0.4), c = c("black", "grey60")), aes(x, y)) +
-        geom_point(shape = 21, size = 3.5, fill = "grey70", color = c("black", "grey60"), stroke = c(1.2, 0.4)) +
-        geom_text(aes(x + 1.0, y, label = label), hjust = 0, size = 3) +
+    p_arch <- ggplot(data.frame(x = 1, y = 2:1, label = c("Inner Arc", "Outer Arc"), s = c(1.2, 0.4), c = c("black", "grey60")), aes(x, y)) +
+        geom_point(shape = 21, size = 4.2, fill = "grey70", color = c("black", "grey60"), stroke = c(1.2, 0.4)) +
+        geom_text(aes(x + 1.0, y, label = label), hjust = 0, size = 4) +
         xlim(0.5, 4.0) +
         ylim(0.5, 2.5) +
         labs(title = "Arc") +
@@ -268,16 +279,23 @@

     cat("      - Preparando leg_hab...\n")
     leg_hab <- get_leg(p_hab)

-    cat("      - Combinando legendas...\n")
-    # Combinar usando patchwork
-    res <- wrap_elements(full = leg_reef) /
-        wrap_elements(full = leg_hab) /
-        wrap_elements(full = p_arch) +
-        plot_layout(heights = c(1.2, 0.8, 1))
+    cat("      - Combinando Reef | HAB (lado a lado)...\n")
+    # NOVO: Reef e HAB lado a lado
+    top_row <- cowplot::plot_grid(
+        leg_reef,
+        leg_hab,
+        ncol = 2,
+        rel_widths = c(1, 1),
+        align = "v"
+    )
+
+    cat("      - Posicionando Arc abaixo...\n")
+    # NOVO: Arc abaixo do par Reef|HAB
+    res <- cowplot::plot_grid(
+        top_row,
+        p_arch,
+        ncol = 1,
+        rel_heights = c(1.5, 1),
+        align = "h"
+    )

     return(res)
 }
```

---

## Verification Checklist

### After Running Script

- [ ] Magnitude figures: Reef and HAB legends are side-by-side (not overlapping)
- [ ] Magnitude figures: Arc legend is positioned below Reef|HAB row
- [ ] All text is ~20% larger than before
- [ ] Axis labels (PC1, PC2) are readable at size 13
- [ ] Panel titles (e.g., "SST (°C)") are readable at size 14
- [ ] Panel tags (a, b, c...) are clearly visible at size 14
- [ ] Loading arrow labels show "SST CV", "DLI CV", "Chl CV" (not CV IDs)
- [ ] Legend text is readable at size 11
- [ ] All 12 figures (6 PNG + 6 PDF) regenerated successfully

### Visual Checks

- [ ] No text overlaps anywhere
- [ ] Legend boxes don't overlap each other
- [ ] Loading arrows and labels are clearly visible
- [ ] Color contrast remains good with increased text
- [ ] Overall figure composition is balanced

---

## Files to Modify

| File | Lines | Changes |
|------|-------|---------|
| `PCA_Visualization_Publication.R` | 51-60 | Expand `loadings_labels` |
| `PCA_Visualization_Publication.R` | 63-78 | Update `theme_publication` |
| `PCA_Visualization_Publication.R` | 219 | Increase `geom_text` size to 4 |
| `PCA_Visualization_Publication.R` | 226-282 | Refactor `create_legend_panel` |
| `PCA_Visualization_Publication.R` | 229, 236 | (OPTIONAL) Increase point size to 4.2 |

---

## Rollback Plan (If Needed)

If the new layout creates issues, the old `create_legend_panel` can be restored:

```r
# OLD VERSION (backup)
res <- wrap_elements(full = leg_reef) /
    wrap_elements(full = leg_hab) /
    wrap_elements(full = p_arch) +
    plot_layout(heights = c(1.2, 0.8, 1))
```

---

*Plan v3 FINAL - Ready for Implementation*
*Date: 2026-01-12*
