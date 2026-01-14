# Diagnostic Report: PCA Figure Composition Solutions

This document analyzes the differences between the current implementation (`PCA_Health_Interactions_Composite_vGeminiPro.R`) and the GLM implementation (`PCA_Health_Interactions_Composite_GLM.R`) regarding layout, spacing, and visual density.

## 1. Architectural Comparison

| Feature | GeminiPro Solution (Current) | GLM Solution (Reference) |
| :--- | :--- | :--- |
| **Assembly Engine** | `patchwork` (Plot Arithmetic) | `cowplot::plot_grid` (Grid Matrix) |
| **Alignment Logic** | Implicit theme-based margins | Explicit `align = "vh"` and `axis = "tblr"` |
| **Axis Limits** | Symmetric Automated (max abs * 1.15) | Manual Custom (Specific ranges per PCA) |
| **Coord System** | `coord_fixed` (Square) | `coord_fixed` (Square) |
| **Fig Dimensions** | 14x15 (Vertical) | 11x12 (Compact) |

---

## 2. Pros and Cons

### solution: `PCA_Health_Interactions_Composite_vGeminiPro.R`

**Pros:**
- **Full Automation**: Dynamically calculates symmetric limits, making it "plug-and-play" for any new dataset.
- **Visual Consistency**: Guaranteed square panels regardless of data distribution.
- **Modern Syntax**: Uses `patchwork` which is the standard for modern R plot composition.

**Cons:**
- **White Space Issue**: `patchwork` can be "greedy" with margins. When combined with `coord_fixed`, it creates vast empty gutters if the output file aspect ratio doesn't perfectly match the grid.
- **Symmetry Penalty**: Forcing X and Y ranges to be identical (symmetric) when data is naturally oblong (e.g., -3 to 3 on X but -2 to 5 on Y) inevitably creates "empty air" inside the charts.
- **Dimension bloat**: Requires very high `ggsave` height (15) to minimize gutters.

### solution: `PCA_Health_Interactions_Composite_GLM.R`

**Pros:**
- **Superior Space Optimization**: Using `cowplot::plot_grid` with `align = "vh"` and `axis = "tblr"` physically "crushes" the plots together by aligning their axis lines, removing almost all inter-plot gutters.
- **Data Density**: By using specific limits (`-3 to 4` and `-2.5 to 5`), the data fills the available square area much more efficiently without "dead zones".
- **Precise Control**: Manual limits allow the user to dictate exactly how much "breathing room" is given to outliers.

**Cons:**
- **Less Flexible**: If a new dataset has values outside these hardcoded limits, the samples will be clipped.
- **Hardcoded**: Requires manual adjustment if the PCA scores change significantly in future analyses.

---

## 3. Conclusion

The **GLM Solution** is technically superior for this specific use case because it prioritizes **visual density** and **mathematical alignment** over automation. 

In scientific publication, the "gutters" (space between plots) are often restricted to the absolute minimum allowed by the publisher. The `cowplot` implementation achieves this by essentially ignoring the individual plot margins and forcing a shared coordinate grid for the entire composite.

### Why GLM won:
1. **The Grid**: `plot_grid` is a more "brute force" approach to layout that is better suited for scientific posters and papers where whitespace is expensive.
2. **The Limits**: GEMINI tried to be mathematically "neat" with symmetric squares; GLM was "practically" neat by providing enough room exactly where the data needed it, resulting in a much "fuller" and more professional look.

**Recommendation**: The architectural logic of the GLM script (`plot_grid` + `manual limits`) should be the standard for final publication versions of these figures.
