---
name: map_figure_layouts
description: Provides publication-quality map and figure creation following Nature/Science standards. Useful when generating forest plots, multi-panel figures, or maps. Keywords: ggplot2, theme, Nature, DPI, patchwork.
---

# Map & Figure Layouts

This skill guides the creation of publication-quality visualizations following Nature/Science journal standards. It applies to maps, multi-panel figures, forest plots, and marginal effects plots, ensuring consistency in aesthetics and technical specifications.

## When to use this skill

- Use this when generating figures for publication in high-impact scientific journals.
- Use this to apply a standardized, color-blind friendly theme to `ggplot2` or `matplotlib` plots.
- Use this when creating multi-panel figures with systematic tagging (A, B, C).
- Use this to export figures in specific formats (PNG, PDF, TIFF) with required DPI settings.

## How to use it

### Step 1: Define the Theme
Apply a standard publication theme (e.g., `theme_classic` base) with black axis text, bold titles, and no grid lines.

### Step 2: Select Color Palettes
Use muted, professional color palettes (e.g., Viridis, Magma, or custom muted blue/orange) that are color-blind friendly.

### Step 3: Generate Specific Plot Types
Follow patterns for:
- **Forest Plots**: For Bayesian coefficients using `ggdist`.
- **Ribbon Plots**: For marginal effects of continuous predictors.
- **Point-range Plots**: For categorical effects.

### Step 4: Assemble Multi-panel Figures
Use `patchwork` (R) or `subplots` (Python) to combine individual plots. Ensure uppercase panel tags (A, B, ...) are present.

### Step 5: Export with Precision
Export at the required resolution:
- **PNG**: 300 DPI for review.
- **TIFF**: 600 DPI (LZW compression) for journal submission.
- **PDF**: Vector format for final publication.

## Common pitfalls

1. **Gray text**: Ensure all axis text and titles are explicitly set to black (`color = "black"`).
2. **Default colors**: Never use default `ggplot2` or `matplotlib` colors; always apply a curated palette.
3. **Low resolution**: Avoid exporting at default DPI; explicitly set 300 or 600 DPI.
4. **Crowded legends**: Default legend position is the bottom; move or simplify if the plot area is compromised.
