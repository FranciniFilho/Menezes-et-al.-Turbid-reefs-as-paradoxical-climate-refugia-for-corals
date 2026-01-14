# Implementation Plan: Publication-Quality Boxplot Figures for Cover and RGR

## Goal Description

Create a new R script that generates publication-quality versions of two existing Python-generated boxplot figures:
1. **`boxplots_facet_reef_hab_auto.png`** (coral cover by reef and habitat)
2. **`boxplot_RGR_por_Reef_e_Hab.png`** (relative growth rate by reef and habitat)

The new figures must match the aesthetic standards of Nature/Science publications and maintain visual consistency with the existing R scripts: `PCA_Health_Interactions_Composite_GLM.R` and `PCA_Environ_vGeminiPro.R`.

---

## Current State Analysis

### Current Python Figures

| Figure | Script | Current Issues |
|--------|--------|----------------|
| Cover Boxplot | `BOXPLOT_cover_REEF_HAB.py` | Basic seaborn styling, Portuguese labels, no ARCH distinction, default colors |
| RGR Boxplot | `Calculate_RGR_&_health_PCA.py` | Similar issues, inconsistent with R publication figures |

### Current Python Aesthetic Limitations

```python
# Current simple styling (from Python scripts)
sns.set_theme(style="whitegrid", palette="muted")
ax.set_title(f"Cobertura (%) de {org} por Reef e HAB", fontsize=14)
```

**Problems:**
- Uses seaborn's default muted palette (not matching `reef_colors` scheme)
- Portuguese labels instead of English
- No inner/outer arc distinction (critical for the publication)
- Default serif fonts and grid styling (doesn't match R figures)
- No individual data points overlaid (reduces transparency)
- Standard boxplot geometry (less elegant than alternatives)

---

## Reference Aesthetic Standards

From the existing R scripts, these are the key design elements to replicate:

### Color Palette (`reef_colors`)
```r
reef_colors <- c(
    "ARC" = "#1f77b4", # Blue
    "ITA" = "#ff7f0e", # Orange
    "PAB" = "#2ca02c", # Green
    "UCR" = "#d62728", # Red
    "TIM" = "#9467bd"  # Purple
)
```

### Shape Scheme (`habitat_shapes`)
```r
habitat_shapes <- c(
    "PA" = 21, # Filled Circle
    "RR" = 22, # Filled Square
    "TP" = 24  # Filled Triangle
)
```

### Theme Publication
```r
theme_publication <- theme_classic(base_size = 14.4) +
    theme(
        text = element_text(color = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        panel.grid.major = element_line(
            color = "grey90", linetype = "dashed", linewidth = 0.3
        ),
        strip.background = element_blank(),
        plot.tag = element_text(face = "bold", size = 14)
    )
```

### Arc Distinction
Inner Arc reefs (ARC, ITA, PAB) are distinguished by **thick black borders**, while Outer Arc reefs (UCR, TIM) have **thin grey borders**.

---

## Proposed Changes

### [NEW] Boxplot_Cover_RGR_Publication.R

**Path**: `#######FINAL_CODES/Boxplot_Cover_RGR_Publication.R`

A new R script that generates both figures with consistent aesthetics.

---

## User Review Required

> [!IMPORTANT]
> **Design Decisions Requiring Approval**

The following design alternatives were evaluated. Please confirm the recommended approach or indicate preferences:

### 1. Plot Type: Boxplot vs. Violin Plot vs. Hybrid

| Option | Description | Recommendation |
|--------|-------------|----------------|
| **A. Classic Boxplot** | Traditional quartile representation | Less informative for small samples |
| **B. Violin Plot** | Shows full distribution density | Can be misleading for small n |
| **C. Half-Violin + Boxplot (Raincloud)** | Combines distribution with summary stats | **RECOMMENDED** |
| **D. Jittered Points + Mean ± SE** | Shows all data points with error bars | Good alternative |

**Recommendation**: **Option C (Raincloud Plot)** using `ggdist::stat_halfeye()` combined with `geom_boxplot()` and jittered raw data points. This approach:
- Shows the full distribution shape (unlike boxplots)
- Avoids the "violin symmetry" that can misrepresent data
- Overlays individual data points for maximum transparency
- Includes summary statistics (median, quartiles)

### 2. Layout Architecture

| Option | Description | Recommendation |
|--------|-------------|----------------|
| **A. Faceted by Reef** | One panel per reef, HAB as color/shape | **RECOMMENDED** - cleaner |
| **B. Reef on X-axis** | HAB as dodged groups (current Python) | Crowded with 5 reefs x 3 habitats |
| **C. Faceted by HAB** | One panel per habitat, Reef as color | Alternative |

**Recommendation**: **Option A (Faceted by Reef)** with a `1x5` horizontal strip layout for Cover, and a separate simplified layout for RGR.

### 3. Data Point Overlay

| Option | Description | Recommendation |
|--------|-------------|----------------|
| **A. No overlay** | Clean boxplots/violins only | Hides sample size variation |
| **B. Jittered points** | Shows all raw data | **RECOMMENDED** |
| **C. Sina (beeswarm)** | Points follow distribution contour | Alternative |

**Recommendation**: **Option B (Jittered Points)** with transparency (`alpha = 0.6`) and small size (`size = 2`), positioned outside the boxplot for clarity.

### 4. Color Encoding

For the boxplot/violin fill, two approaches are possible:

| Option | Description | Recommendation |
|--------|-------------|----------------|
| **A. Fill by REEF** | Uses `reef_colors` palette | Redundant if faceted by Reef |
| **B. Fill by HAB** | Uses habitat-specific colors | **RECOMMENDED** |
| **C. Fill by ARCH** | Inner=solid, Outer=hatched | Alternative |

**Recommendation**: **Option B (Fill by HAB)** with a new complementary palette that doesn't clash with `reef_colors`:

```r
habitat_fill_colors <- c(
    "PA" = "#56B4E9",  # Light blue (Reef Wall)
    "RR" = "#E69F00",  # Orange (Rocky Reef)
    "TP" = "#009E73"   # Teal (Reef Top)
)
```

### 5. Organism Scope (UPDATED per User Feedback)

> [!NOTE]
> **Cover/Abundance Figure**: Will include **5 key organisms** to show community composition:
> - CCA (Crustose Coralline Algae)
> - Cyano (Cyanobacteria)
> - Macroalgae
> - M. hispida (Mussismilia hispida)
> - Turf
>
> **RGR Figure**: Will focus **only on M. hispida** (the focal species for growth analysis).

---

## Detailed Implementation

### Script Structure

```
Boxplot_Cover_RGR_Publication.R
├── 1. PACKAGES
├── 2. PATHS
├── 3. VISUAL SCHEME (reef_colors, habitat_shapes, organism_colors, theme_publication)
├── 4. HELPER FUNCTIONS
│   ├── load_cover_data()
│   ├── load_rgr_data()
│   ├── create_raincloud_panel()
│   └── create_legend_strip()
├── 5. FIGURE GENERATION
│   ├── create_cover_figure()  # Multi-organism
│   └── create_rgr_figure()    # M. hispida only
├── 6. MAIN EXECUTION
│   └── main()
```

### Figure 1: Cover Raincloud Plot (Multi-Organism)

**Layout**: 5-row x 5-column grid (organisms as rows, reefs as columns)

```
              ARC        ITA        PAB        UCR        TIM
           (inner)    (inner)    (inner)    (outer)    (outer)
         ┌─────────┬─────────┬─────────┬─────────┬─────────┐
   CCA   │         │         │         │         │         │
         ├─────────┼─────────┼─────────┼─────────┼─────────┤
  Cyano  │         │         │         │         │         │
         ├─────────┼─────────┼─────────┼─────────┼─────────┤
Macroalg │         │         │         │         │         │
         ├─────────┼─────────┼─────────┼─────────┼─────────┤
M.hispida│         │         │         │         │         │
         ├─────────┼─────────┼─────────┼─────────┼─────────┤
  Turf   │         │         │         │         │         │
         └─────────┴─────────┴─────────┴─────────┴─────────┘
                    ↓ Shared Legend Below ↓
              [Habitat: PA | RR | TP]  [Arc: Inner | Outer]
```

**Alternative Layout** (if 5x5 is too crowded): Facet by Organism (rows) with Reef on X-axis, HAB as fill color

```
         ┌────────────────────────────────────────────────────┐
   CCA   │  [ARC] [ITA] [PAB] [UCR] [TIM]  (HAB as colors)    │
         ├────────────────────────────────────────────────────┤
  Cyano  │  [ARC] [ITA] [PAB] [UCR] [TIM]                     │
         ├────────────────────────────────────────────────────┤
Macroalg │  [ARC] [ITA] [PAB] [UCR] [TIM]                     │
         ├────────────────────────────────────────────────────┤
M.hispida│  [ARC] [ITA] [PAB] [UCR] [TIM]                     │
         ├────────────────────────────────────────────────────┤
  Turf   │  [ARC] [ITA] [PAB] [UCR] [TIM]                     │
         └────────────────────────────────────────────────────┘
```

**Key Features**:
- Rows: ORGANISMO (5 organisms)
- X-axis: REEF (ARC, ITA, PAB, UCR, TIM) or HAB depending on layout
- Y-axis: Cover (%)
- HAB encoded as fill color (PA=blue, RR=orange, TP=teal)
- Raincloud geometry (half-violin + boxplot + jittered points)
- Arc distinction visible in REEF ordering (Inner: ARC, ITA, PAB first; Outer: UCR, TIM last)
- English labels with proper organism formatting (e.g., *M. hispida* in italics)

### Figure 2: RGR Raincloud Plot (M. hispida Only)

**Layout**: 1x5 horizontal strip (one per reef) or 2x3 grid

```
┌─────────┬─────────┬─────────┬─────────┬─────────┐
│   ARC   │   ITA   │   PAB   │   UCR   │   TIM   │
│ (inner) │ (inner) │ (inner) │ (outer) │ (outer) │
└─────────┴─────────┴─────────┴─────────┴─────────┘
              ↓ Shared Legend Below ↓
        [Habitat: PA | RR | TP]  [Arc: Inner | Outer]
```

**Key Features**:
- X-axis: HAB (PA, RR, TP)
- Y-axis: Relative Growth Rate (year^-1)
- Reference line at y = 0 (zero growth, red dashed)
- Raincloud geometry (half-violin + boxplot + jittered points)
- HAB encoded as fill color
- Title: "Relative Growth Rate of *Mussismilia hispida*"

---

## Verification Plan

### Automated Tests

1. **Script Execution**:
   ```r
   source("Boxplot_Cover_RGR_Publication.R")
   ```
   Expected: No errors, all figures saved successfully.

2. **Output File Existence**:
   - `Cover_MultiOrganism_Raincloud_Publication.png` (300 DPI)
   - `Cover_MultiOrganism_Raincloud_Publication.pdf` (vector)
   - `RGR_Mussismilia_hispida_Raincloud_Publication.png` (300 DPI)
   - `RGR_Mussismilia_hispida_Raincloud_Publication.pdf` (vector)

3. **Data Integrity Checks** (within script):
   - Verify all 5 reefs are present
   - Verify all 3 habitats are present
   - Log sample sizes per group

### Manual Verification

1. **Visual Comparison**:
   - Compare new figures side-by-side with existing PCA figures
   - Verify color palette consistency
   - Verify font sizes and weights match
   - Verify Arc distinction is clear

2. **Publication Checklist**:
   - [ ] All labels in English
   - [ ] Axis labels with proper units (% for cover, year^-1 for RGR)
   - [ ] Consistent color scheme with PCA figures
   - [ ] Clear legend (not overlapping with data)
   - [ ] No clipped elements
   - [ ] High resolution (300+ DPI for PNG)
   - [ ] Vector export (PDF) available

---

## Dependencies

### Required R Packages

```r
libs <- c(
    "readxl",     # Read Excel files
    "readr",      # Read CSV files
    "dplyr",      # Data manipulation
    "ggplot2",    # Base plotting
    "ggdist",     # stat_halfeye for raincloud plots
    "patchwork",  # Panel composition
    "cowplot",    # Legend extraction and assembly
    "scales",     # Rescale functions
    "ggbeeswarm", # Optional: sina plots (alternative to jitter)
    "grid"        # For graphical parameters
)
```

### Data Sources

| Figure | Data File | Key Columns |
|--------|-----------|-------------|
| Cover | `dados_integrados_long_format.csv` | REEF, HAB, ORGANISMO, COBERTURA |
| RGR | `resultados_biologicos_por_colonia.csv` | REEF, HAB, RGR |

**Data Paths** (from existing Python scripts):
- Cover: `C:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\########NEW RESULTS\BOX_PLOTS\dados_integrados_long_format.csv`
- RGR: `C:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\########NEW RESULTS\output_ANALISE_BIOLOGICA_boxplot_PCAnew\resultados_biologicos_por_colonia.csv`

---

## Output Location

All figures will be saved to:

```
C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_RESULTS/Boxplot_Publication_Figures/
```

---

## Timeline Estimate

| Phase | Estimated Time |
|-------|----------------|
| Script development | 1-2 hours |
| Visual testing and iteration | 30 min |
| Documentation | 15 min |
| **Total** | **~2-3 hours** |

---

## Alternative Approaches Considered

### 1. Python with `plotnine` (ggplot2 for Python)

**Pros**: Keep everything in Python
**Cons**: Less mature than ggplot2, harder to replicate exact R aesthetics

**Decision**: Use R for full aesthetic parity with existing figures.

### 2. Simple Boxplot Enhancement

**Pros**: Minimal code changes
**Cons**: Doesn't show full data distribution, less informative

**Decision**: Adopt raincloud/hybrid approach for maximum information density.

### 3. Separate Scripts

**Pros**: Modular
**Cons**: Redundant theme/function definitions

**Decision**: Single unified script with shared visual scheme.

---

## Summary

This plan proposes a new R script that will:

1. **Replicate the exact visual standards** from `PCA_Health_Interactions_Composite_GLM.R`
2. **Use raincloud plots** (half-violin + boxplot + jittered points) for maximum information
3. **Show 5 organisms** (CCA, Cyano, Macroalgae, M. hispida, Turf) in the Cover figure
4. **Focus on M. hispida only** for the RGR (growth) figure
5. **Facet by REEF and ORGANISMO** for cleaner presentation of Cover data
6. **Include Arc distinction** via panel styling or ordering
7. **Generate both PNG (300 DPI) and PDF (vector)** outputs

The remaining design decisions (layout options, color encoding) can be refined during implementation based on data visualization best practices.

