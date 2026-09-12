# ============================================================================
# 08_Visualize_JSDM_Dominance_Probabilities.R  [v2 – corrected figures]
# ============================================================================
# Fixes vs. previous version:
#
#   FIG 8  – Added P(CCA > Coral) alongside P(Macro > Coral) as grouped bars,
#             so both competitors are compared per ARC × HAB.
#
#   FIG 9  – Map was squished horizontally due to coord_fixed(ratio=1.5).
#             Fixed with coord_quickmap(). Legend was all-black because
#             shapes 21/22 require fill= but the legend didn't inherit it.
#             Fixed by switching to solid shapes (16/17) + color= for status.
#             Legend keys are now explicitly overridden for clarity.
#
#   FIG 10 – Added panel B showing community composition of stress sites
#             (previously only refugia composition was shown). Layout is now
#             3-panel: (A) Refugia | (B) Stress / (C) Mean comparison.
#             Added TURF+CYANO as "Other benthic" for full composition view.
#
#   FIG 11 – Changed from categorical fill to grouped bars showing both
#             P(Macro > Coral) and P(CCA > Coral) per reef, ordered by
#             P(Macro), so the dominance gradient is unambiguous.
# ============================================================================

rm(list = ls()); gc()

# ── 1. PACKAGES ────────────────────────────────────────────────────────────

libs <- c("ggplot2", "dplyr", "tidyr", "patchwork", "cowplot",
          "scales", "ggrepel")
missing_pkgs <- libs[!libs %in% installed.packages()[, "Package"]]
if (length(missing_pkgs) > 0) install.packages(missing_pkgs, dependencies = TRUE)
invisible(lapply(libs, library, character.only = TRUE))

# ── 2. PUBLICATION THEME ────────────────────────────────────────────────────

theme_pub <- function(base_size = 13) {
  theme_classic(base_size = base_size) +
    theme(
      text            = element_text(color = "black"),
      axis.text       = element_text(color = "black", size = 11),
      axis.title      = element_text(face = "bold", size = 13),
      axis.line       = element_line(color = "black", linewidth = 0.5),
      plot.title      = element_text(face = "bold", size = 14, hjust = 0),
      plot.subtitle   = element_text(size = 10, color = "gray30"),
      plot.caption    = element_text(size = 8, color = "gray50", hjust = 1),
      legend.position = "bottom",
      legend.text     = element_text(size = 10),
      legend.title    = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "gray93", color = NA),
      strip.text      = element_text(face = "bold", size = 12),
      panel.grid      = element_blank()
    )
}
theme_set(theme_pub())

# Colour palette (Okabe-Ito, colour-blind friendly)
COMP_COLS  <- c("Macroalgae > Coral" = "#D55E00", "CCA > Coral" = "#009E73")
HAB_COLS   <- c("Reef Wall (PA)"   = "#D55E00",
                "Reef Top (TP)"    = "#009E73",
                "Rocky Reef (RR)"  = "#CC79A7")
STATUS_COLS <- c("Coral-favored"         = "#0072B2",
                 "Uncertain"             = "#999999",
                 "Macroalgae-advantaged" = "#E69F00",
                 "Macroalgae-dominated"  = "#D55E00")
COMP_COLS_STACK <- c("Coral (M. hispida)" = "#0072B2",
                     "CCA"               = "#009E73",
                     "Macroalgae"        = "#D55E00",
                     "Other benthic"     = "#CCCCCC")

# ── 3. I/O DIRECTORIES ──────────────────────────────────────────────────────

input_dir  <- "#######FINAL_RESULTS/JSDM_Dominance/"
output_dir <- "#######FINAL_RESULTS/JSDM_Dominance/"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

save_fig <- function(fig, name, w, h) {
  ggsave(file.path(output_dir, paste0(name, ".png")), fig,
         width = w, height = h, dpi = 600, bg = "white")
  ggsave(file.path(output_dir, paste0(name, ".pdf")), fig,
         width = w, height = h)
  cat(sprintf("  Saved: %s\n", name))
}

# ── 4. LOAD & ENRICH DATA ───────────────────────────────────────────────────

cat("\n=== Loading data ===\n")
dominance_full   <- read.csv(file.path(input_dir, "jsdm_dominance_probabilities_full.csv"),
                              stringsAsFactors = FALSE)
summary_reef_hab <- read.csv(file.path(input_dir, "jsdm_dominance_summary_reef_hab.csv"),
                              stringsAsFactors = FALSE)
summary_arch_hab <- read.csv(file.path(input_dir, "jsdm_dominance_summary_arch_hab.csv"),
                              stringsAsFactors = FALSE)

cat(sprintf("  Full data: %d obs, %d cols\n", nrow(dominance_full), ncol(dominance_full)))

# Compute "Other benthic" = TURF + CYANO (remainder from Dirichlet sum ≈ 1)
dominance_full <- dominance_full %>%
  mutate(
    mean_other = pmax(0, 1 - mean_MUSSISMILIA - mean_MACROALGAE - mean_CCA),
    ARC = factor(ARCH, levels = c("inner", "outer"),
                 labels = c("Inner Arc", "Outer Arc")),
    HAB_label = factor(HAB,
                       levels = c("PA", "TP", "RR"),
                       labels = c("Reef Wall\n(PA)", "Reef Top\n(TP)", "Rocky Reef\n(RR)"))
  )

# Site-level summary (collapse across years)
site_level <- dominance_full %>%
  group_by(REEF_NAME, SITE, HAB, HAB_label, ARCH, ARC, LATITUDE, LONGITUDE) %>%
  summarise(
    P_macro_coral    = mean(P_macro_coral,    na.rm = TRUE),
    P_cca_coral      = mean(P_cca_coral,      na.rm = TRUE),
    mean_MUSSISMILIA = mean(mean_MUSSISMILIA, na.rm = TRUE),
    mean_MACROALGAE  = mean(mean_MACROALGAE,  na.rm = TRUE),
    mean_CCA         = mean(mean_CCA,         na.rm = TRUE),
    mean_other       = mean(mean_other,        na.rm = TRUE),
    status = names(which.max(table(status))),
    .groups = "drop"
  ) %>%
  mutate(
    status = factor(status,
                    levels = c("Coral-favored", "Uncertain",
                               "Macroalgae-advantaged", "Macroalgae-dominated")),
    site_label = paste(SITE, HAB, sep = "_")
  )

cat(sprintf("  Site-level: %d unique SITE x HAB combinations\n", nrow(site_level)))
cat("  Status distribution:\n")
print(table(site_level$status))

# ── Env. variability data from CV_ALL window (JSDM model predictor) ──────────
cv_all_path <- "#######FINAL_DATA/02_PCA_environmental/CV_ALL/dados_abundancia_integrados_long_format.csv"

first_line_cv <- readLines(cv_all_path, n = 1L, encoding = "UTF-8")
delim_cv      <- if (grepl(";", first_line_cv, fixed = TRUE)) ";" else ","
dec_cv        <- if (delim_cv == ";") "," else "."        # BR locale: ";" sep + "," decimal
env_raw       <- read.csv(cv_all_path, sep = delim_cv, dec = dec_cv,
                          stringsAsFactors = FALSE, fileEncoding = "UTF-8-BOM")

# Site × HAB means – PCA scores AND raw inter-annual CVs (collapse across years)
# Raw CVs are extracted alongside PCA scores to avoid CV_ALL axis-orientation
# artefacts: in CV_ALL the SST loading on PC1 flips sign relative to CV_02/CV_30,
# so the raw metrics are used directly for Panel D (unambiguous interpretation).
env_site <- env_raw %>%
  group_by(SITE, HAB) %>%
  summarise(
    PC1_VAR = mean(PC1_VARIABILITY, na.rm = TRUE),
    PC2_VAR = mean(PC2_VARIABILITY, na.rm = TRUE),
    SST_CV  = mean(SST_CV_ALL,  na.rm = TRUE),   # raw thermal variability (%)
    CHL_CV  = mean(CHL_CV_ALL,  na.rm = TRUE),   # raw productivity variability (%)
    DLI_CV  = mean(CV_DLI_LOCAL, na.rm = TRUE),  # raw light variability (%)
    .groups = "drop"
  )

# Left-join into site_level (all downstream subsets inherit these columns)
site_level <- site_level %>%
  left_join(env_site, by = c("SITE", "HAB"))

cat(sprintf("  Env. variability joined: %d/%d sites matched\n",
            sum(!is.na(site_level$PC1_VAR)), nrow(site_level)))

# ── 5. FIGURE 8: P(Competitor > Coral) by ARC × HAB ────────────────────────
# Two competitors shown as grouped bars within each ARC × HAB cell.
# ─────────────────────────────────────────────────────────────────────────────

cat("\n=== Figure 8: Competitor probabilities by ARC × HAB ===\n")

arch_hab_df <- dominance_full %>%
  group_by(ARC, HAB_label) %>%
  summarise(
    n_obs        = n(),
    macro_mean   = mean(P_macro_coral, na.rm = TRUE),
    macro_min    = min(P_macro_coral,  na.rm = TRUE),
    macro_max    = max(P_macro_coral,  na.rm = TRUE),
    cca_mean     = mean(P_cca_coral,   na.rm = TRUE),
    cca_min      = min(P_cca_coral,    na.rm = TRUE),
    cca_max      = max(P_cca_coral,    na.rm = TRUE),
    .groups = "drop"
  )

# Pivot to long for grouped bars
arch_hab_long <- bind_rows(
  arch_hab_df %>% transmute(ARC, HAB_label, n_obs,
                             Competitor = "Macroalgae > Coral",
                             P_mean = macro_mean, P_min = macro_min, P_max = macro_max),
  arch_hab_df %>% transmute(ARC, HAB_label, n_obs,
                             Competitor = "CCA > Coral",
                             P_mean = cca_mean,   P_min = cca_min,   P_max = cca_max)
) %>%
  mutate(Competitor = factor(Competitor,
                             levels = c("Macroalgae > Coral", "CCA > Coral")))

fig8 <- ggplot(arch_hab_long,
               aes(x = HAB_label, y = P_mean, fill = Competitor,
                   group = Competitor)) +
  geom_bar(stat = "identity", position = position_dodge(0.75),
           width = 0.65, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = P_min, ymax = P_max),
                position = position_dodge(0.75), width = 0.25, linewidth = 0.7) +
  geom_hline(yintercept = c(0.05, 0.50, 0.95),
             linetype = "dashed", color = "gray40", linewidth = 0.45,
             alpha = 0.8) +
  facet_wrap(~ ARC, nrow = 1) +
  scale_fill_manual(values = COMP_COLS,
                    name   = expression(paste(P, "(Competitor > ", italic("M. hispida"), ")"))) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2),
                     labels = number_format(accuracy = 0.1)) +
  labs(
    title    = "Probability of Competitor Dominance over Coral by Arc and Habitat",
    subtitle = expression(paste("Posterior ", P(Competitor > italic("M. hispida")),
                                " from JSDM Dirichlet model")),
    x        = "Habitat",
    y        = expression(P(Competitor > Coral)),
    caption  = "Error bars: range across sites within each ARC × HAB combination.\nDashed lines: 0.05 (refugia), 0.50 (uncertain) and 0.95 (dominance) thresholds."
  ) +
  theme(
    strip.text      = element_text(face = "bold", size = 13),
    axis.text.x     = element_text(size = 10),
    legend.key.size = unit(0.5, "cm")
  )

save_fig(fig8, "FIGURE_8_Dominance_by_ARC_HAB", 11, 7)

# ── 6. FIGURE 9: Spatial map (fixed legend + aspect ratio) ──────────────────
# FIX: switched from fill= (broken with shapes 21/22) to color= + solid shapes
# FIX: coord_quickmap() replaces coord_fixed(1.5) which squished the map
# ─────────────────────────────────────────────────────────────────────────────

cat("\n=== Figure 9: Spatial dominance map ===\n")

# Label only the most extreme sites
p05 <- quantile(site_level$P_macro_coral, 0.10, na.rm = TRUE)
p95 <- quantile(site_level$P_macro_coral, 0.90, na.rm = TRUE)

# Compute data extents and expand longitude to ~double for readability.
# At ~18°S the bank is taller than wide in degrees (lat_span > lon_span),
# so coord_quickmap() produces a compressed x-axis.
lon_rng  <- range(site_level$LONGITUDE, na.rm = TRUE)
lat_rng  <- range(site_level$LATITUDE,  na.rm = TRUE)
lon_span <- diff(lon_rng)
lat_span <- diff(lat_rng)

# The data has lat_span >> lon_span (~1.1° vs ~0.5°), so coord_quickmap()
# produces a map taller than wide. To ensure x >= 2× original we need
# lon_span_final × cos(18°) > lat_span_padded, i.e. lon_span > lat/cos(18°).
# Padding lon by 1.2× each side gives lon_span × 3.4 which is sufficient.
lon_xlim <- c(lon_rng[1] - lon_span * 1.2,  lon_rng[2] + lon_span * 1.2)
lat_ylim <- c(lat_rng[1] - lat_span * 0.12, lat_rng[2] + lat_span * 0.12)

cat(sprintf("  Lon range (data): %.3f to %.3f  →  xlim: %.3f to %.3f\n",
            lon_rng[1], lon_rng[2], lon_xlim[1], lon_xlim[2]))

fig9 <- ggplot(site_level, aes(x = LONGITUDE, y = LATITUDE)) +
  # Background points (all sites, semi-transparent)
  geom_point(aes(color = status, shape = ARC, size = P_macro_coral),
             alpha = 0.88, stroke = 1.2) +
  geom_text_repel(
    data = site_level %>% filter(P_macro_coral <= p05 | P_macro_coral >= p95),
    aes(label = site_label),
    size = 2.8, max.overlaps = 20, box.padding = 0.4,
    segment.color = "gray50", segment.size = 0.4,
    min.segment.length = 0.2
  ) +
  # Aesthetics
  scale_color_manual(
    values = STATUS_COLS,
    name   = "Dominance status",
    guide  = guide_legend(
      override.aes = list(size = 4, stroke = 1.5, shape = 16),
      order = 1
    )
  ) +
  scale_shape_manual(
    values = c("Inner Arc" = 16, "Outer Arc" = 17),
    name   = "Arc position",
    guide  = guide_legend(
      override.aes = list(size = 4, color = "black"),
      order = 2
    )
  ) +
  scale_size_continuous(
    range  = c(2, 7),
    name   = expression(P(Macroalgae > Coral)),
    breaks = c(0.05, 0.25, 0.50, 0.75, 0.95),
    labels = c("0.05\n(Refugia)", "0.25", "0.50", "0.75", "0.95\n(Stress)"),
    guide  = guide_legend(order = 3,
                          override.aes = list(color = "gray40", shape = 16))
  ) +
  coord_quickmap(xlim = lon_xlim, ylim = lat_ylim) +  # expanded x for readability
  labs(
    title    = "Spatial Distribution of Coral vs. Macroalgae Competitive Dominance",
    subtitle = "Each point: unique SITE × HAB combination (P averaged across years). Shape = Arc position.",
    x        = "Longitude (°W)",
    y        = "Latitude (°S)",
    caption  = paste("Point size: P(macroalgae > coral).",
                     "Labels shown for top/bottom 10 % extreme sites.",
                     sep = "\n")
  ) +
  theme(
    legend.position  = "right",
    legend.box       = "vertical",
    legend.spacing.y = unit(0.3, "cm")
  )

save_fig(fig9, "FIGURE_9_Dominance_Map", 16, 8)

# ── 7. FIGURE 10: Refugia vs Stress – 3-panel composition ──────────────────
# Panel A: community composition at REFUGIA sites (P < 0.05)
# Panel B: community composition at STRESS sites (P > 0.95)
# Panel C: mean comparison summary (refugia vs stress)
# ─────────────────────────────────────────────────────────────────────────────

cat("\n=== Figure 10: Refugia vs Stress composition ===\n")

refugia_sites <- site_level %>%
  filter(P_macro_coral < 0.05) %>%
  arrange(P_macro_coral)

stress_sites <- site_level %>%
  filter(P_macro_coral > 0.95) %>%
  arrange(desc(P_macro_coral))

cat(sprintf("  CCA-dominated sites (panel A): %d\n",    nrow(refugia_sites)))
cat(sprintf("  Macroalgae-dominated sites (panel B): %d\n", nrow(stress_sites)))

# Fallback if thresholds yield no sites
if (nrow(refugia_sites) == 0)
  refugia_sites <- site_level %>% arrange(P_macro_coral) %>% head(8)
if (nrow(stress_sites) == 0)
  stress_sites  <- site_level %>% arrange(desc(P_macro_coral)) %>% head(8)

# Tidy helper: pivot composition to long
make_comp_long <- function(df, cat_label) {
  df %>%
    mutate(site_hab = paste0(SITE, "\n", HAB)) %>%
    select(site_hab, ARC,
           `Coral (M. hispida)` = mean_MUSSISMILIA,
           `CCA`                = mean_CCA,
           `Macroalgae`         = mean_MACROALGAE,
           `Other benthic`      = mean_other) %>%
    pivot_longer(cols = c(`Coral (M. hispida)`, CCA, Macroalgae, `Other benthic`),
                 names_to = "Functional_group", values_to = "Proportion") %>%
    mutate(
      Functional_group = factor(Functional_group,
                                levels = c("Coral (M. hispida)", "CCA",
                                           "Macroalgae", "Other benthic")),
      Category = cat_label
    )
}

refugia_long <- make_comp_long(refugia_sites, "CCA-dominated Zones")
stress_long  <- make_comp_long(stress_sites,  "Macroalgae-dominated Zones")

# Preserve site ordering within each panel
refugia_order <- unique(refugia_long$site_hab)
stress_order  <- unique(stress_long$site_hab)

refugia_long <- refugia_long %>%
  mutate(site_hab = factor(site_hab, levels = refugia_order))
stress_long <- stress_long %>%
  mutate(site_hab = factor(site_hab, levels = stress_order))

# Zone colour palette for Panel D violin plot
ZONE_COLS <- c(
  "CCA-dominated\n(Refugia)"       = "#0072B2",   # blue (coral-favored)
  "Macroalgae-dominated\n(Stress)" = "#D55E00"    # orange (stress)
)

# Panel A – Refugia
fig10a <- ggplot(refugia_long,
                 aes(x = site_hab, y = Proportion, fill = Functional_group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75,
           color = "black", linewidth = 0.2) +
  scale_fill_manual(values = COMP_COLS_STACK, name = "Functional group") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.01)) +
  labs(
    title    = "(a) CCA-dominated Zones",
    subtitle = expression(paste(P(Macro > Coral), " < 0.05")),
    x        = NULL, y = "Posterior mean proportion"
  ) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 8.5),
        legend.position = "none")

# Panel B – Stress
fig10b <- ggplot(stress_long,
                 aes(x = site_hab, y = Proportion, fill = Functional_group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.75,
           color = "black", linewidth = 0.2) +
  scale_fill_manual(values = COMP_COLS_STACK, name = "Functional group") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.01)) +
  labs(
    title    = "(b) Macroalgae-dominated Zones",
    subtitle = expression(paste(P(Macro > Coral), " > 0.95")),
    x        = NULL, y = "Posterior mean proportion"
  ) +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 8.5),
        legend.position = "none")

# Panel C – Mean comparison
summary_comp <- bind_rows(
  refugia_sites %>% summarise(
    Category              = "CCA-dominated Zones",
    `Coral (M. hispida)`  = mean(mean_MUSSISMILIA),
    `CCA`                 = mean(mean_CCA),
    `Macroalgae`          = mean(mean_MACROALGAE),
    `Other benthic`       = mean(mean_other)),
  stress_sites %>% summarise(
    Category              = "Macroalgae-dominated Zones",
    `Coral (M. hispida)`  = mean(mean_MUSSISMILIA),
    `CCA`                 = mean(mean_CCA),
    `Macroalgae`          = mean(mean_MACROALGAE),
    `Other benthic`       = mean(mean_other))
) %>%
  pivot_longer(cols = -Category, names_to = "Functional_group", values_to = "Proportion") %>%
  mutate(
    Functional_group = factor(Functional_group,
                              levels = c("Coral (M. hispida)", "CCA",
                                         "Macroalgae", "Other benthic")),
    Category = factor(Category, levels = c("CCA-dominated Zones", "Macroalgae-dominated Zones"))
  )

fig10c <- ggplot(summary_comp,
                 aes(x = Category, y = Proportion, fill = Functional_group)) +
  geom_bar(stat = "identity", position = "stack", width = 0.55,
           color = "black", linewidth = 0.3) +
  geom_text(
    aes(label = ifelse(Proportion > 0.03,
                       sprintf("%.1f%%", Proportion * 100), "")),
    position = position_stack(vjust = 0.5),
    size = 3.5, fontface = "bold", color = "white"
  ) +
  scale_fill_manual(
    values = COMP_COLS_STACK,
    name   = "Functional group",
    labels = c(
      expression(paste("Coral (", italic("M. hispida"), ")")),
      "CCA",
      "Macroalgae",
      "Other benthic"
    )
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title    = "(c) Mean community composition: CCA-dominated vs Macroalgae-dominated Zones",
    subtitle = "Posterior means averaged across sites within each category",
    x        = NULL, y = "Posterior mean proportion"
  ) +
  theme(axis.text.x   = element_text(face = "bold", size = 12),
        legend.position = "right")

# ── Panel D: Environmental variability – violin + jitter (raw CVs) ────────────
# PCA scores from CV_ALL are NOT used here because the SST loading on PC1 flips
# sign between CV_02/CV_30 and CV_ALL (artefact of PCA axis orientation), making
# the direction of the PC1 gradient ambiguous for interpretation. Instead, the
# three most informative raw inter-annual CVs: SST_CV (thermal stability),
# CHL_CV (productivity / nutrient-pulse hypothesis), DLI_CV (light variability).
env_violin_df <- bind_rows(
  refugia_sites %>% mutate(Zone = "CCA-dominated\n(Refugia)"),
  stress_sites  %>% mutate(Zone = "Macroalgae-dominated\n(Stress)")
) %>%
  select(Zone, SITE, HAB, SST_CV, CHL_CV, DLI_CV) %>%
  pivot_longer(cols = c(SST_CV, CHL_CV, DLI_CV),
               names_to  = "Variable",
               values_to = "Value") %>%
  mutate(
    Zone = factor(Zone, levels = c("CCA-dominated\n(Refugia)",
                                   "Macroalgae-dominated\n(Stress)")),
    Variable = factor(
      dplyr::recode(Variable,
        "SST_CV" = "SST inter-annual CV\n(Thermal, %)",
        "CHL_CV" = "Chlorophyll inter-annual CV\n(Productivity, %)",
        "DLI_CV" = "DLI inter-annual CV\n(Light, %)"
      ),
      levels = c("SST inter-annual CV\n(Thermal, %)",
                 "Chlorophyll inter-annual CV\n(Productivity, %)",
                 "DLI inter-annual CV\n(Light, %)")
    )
  )

set.seed(42)
fig10d <- ggplot(env_violin_df,
                 aes(x = Zone, y = Value, fill = Zone, color = Zone)) +
  geom_violin(trim = FALSE, alpha = 0.35, linewidth = 0.5, color = NA) +
  geom_jitter(width = 0.08, size = 2.4, alpha = 0.88, shape = 16) +
  facet_wrap(~ Variable, scales = "free_y", nrow = 1) +
  scale_fill_manual(values  = ZONE_COLS, guide = "none") +
  scale_color_manual(values = ZONE_COLS, guide = "none") +
  labs(
    title    = "(d) Environmental drivers: refugia vs. stress zones",
    subtitle = "CV_ALL window \u00b7 Raw inter-annual coefficients of variation (JSDM predictor window)",
    x        = NULL,
    y        = "Coefficient of variation (%)",
    caption  = paste(
      "Raw inter-annual CVs used instead of PCA scores to avoid CV_ALL axis-orientation artefact",
      "(SST loading on PC1 flips sign between CV_02/CV_30 and CV_ALL).",
      "Higher CHL and DLI CVs in refugia reflect dynamic nutrient/light regimes; lower CVs in stress zones = chronic oligotrophy.",
      sep = "\n"
    )
  ) +
  theme(
    axis.text.x = element_text(size = 9.5, lineheight = 1.1),
    strip.text  = element_text(face = "bold", size = 11)
  )

# Assemble with patchwork only (avoids cowplot/ggsave incompatibility)
fig10 <- (
  (fig10a + theme(legend.position = "none")) |
  (fig10b + theme(legend.position = "none"))
) /
  (fig10c | fig10d) +
  plot_layout(heights = c(1.2, 1.1), guides = "collect") +
  plot_annotation(
    title   = "Community Composition and Environmental Context: CCA-dominated vs Macroalgae-dominated Zones",
    caption = paste(
      "Stacked bars: posterior mean proportions from JSDM Dirichlet model (CV_ALL, WeaklyInformative prior).",
      "A: CCA-dominated zones [P(Macro > Coral) < 0.05].  B: Macroalgae-dominated zones [P(Macro > Coral) > 0.95].",
      "C: Mean composition comparison across categories.",
      "D: Raw inter-annual CVs (SST and Chlorophyll) from CV_ALL window — shown instead of PCA scores to avoid axis-orientation artefact.",
      sep = "\n"
    )
  ) &
  theme(legend.position = "bottom")

save_fig(fig10, "FIGURE_10_Refugia_vs_Stress", 22, 14)   # width: 18 → 22 (3 facets in panel D)

# ── 8. FIGURE 11: Reef-level gradient – both competitors ────────────────────
# Grouped bars: P(Macro > Coral) vs P(CCA > Coral) per reef,
# ordered by P(Macro > Coral) descending so the gradient is explicit.
# ─────────────────────────────────────────────────────────────────────────────

cat("\n=== Figure 11: Reef-level dominance gradient ===\n")

reef_summary <- dominance_full %>%
  group_by(REEF_NAME, ARC) %>%
  summarise(
    P_macro     = mean(P_macro_coral, na.rm = TRUE),
    P_macro_min = min(P_macro_coral,  na.rm = TRUE),
    P_macro_max = max(P_macro_coral,  na.rm = TRUE),
    P_cca       = mean(P_cca_coral,   na.rm = TRUE),
    P_cca_min   = min(P_cca_coral,    na.rm = TRUE),
    P_cca_max   = max(P_cca_coral,    na.rm = TRUE),
    n_site_hab  = n_distinct(paste(SITE, HAB)),
    .groups = "drop"
  ) %>%
  # Order within each ARC group by descending P_macro
  group_by(ARC) %>%
  arrange(ARC, desc(P_macro)) %>%
  ungroup() %>%
  mutate(REEF_NAME = factor(REEF_NAME, levels = unique(REEF_NAME)))

reef_long <- bind_rows(
  reef_summary %>% transmute(REEF_NAME, ARC,
                              Competitor = "Macroalgae > Coral",
                              P = P_macro, P_min = P_macro_min, P_max = P_macro_max),
  reef_summary %>% transmute(REEF_NAME, ARC,
                              Competitor = "CCA > Coral",
                              P = P_cca,   P_min = P_cca_min,   P_max = P_cca_max)
) %>%
  mutate(Competitor = factor(Competitor, levels = c("Macroalgae > Coral", "CCA > Coral")))

fig11 <- ggplot(reef_long,
                aes(x = REEF_NAME, y = P, fill = Competitor, group = Competitor)) +
  geom_bar(stat = "identity", position = position_dodge(0.75),
           width = 0.65, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = P_min, ymax = P_max),
                position = position_dodge(0.75),
                width = 0.25, linewidth = 0.65) +
  geom_hline(yintercept = c(0.05, 0.50, 0.95),
             linetype = "dashed", color = "gray40", linewidth = 0.4) +
  scale_fill_manual(values = COMP_COLS,
                    name   = expression(paste(P, "(Competitor > ",
                                              italic("M. hispida"), ")"))) +
  # ARC shown via x-axis label colour using axis.text customisation
  # (geom_tile removed; colour is conveyed by facet via REEF_NAME label colouring)
  facet_grid(. ~ ARC, scales = "free_x", space = "free_x") +
  scale_y_continuous(limits = c(0, 1.05),
                     breaks = seq(0, 1, 0.2),
                     labels = number_format(accuracy = 0.1)) +
  labs(
    title    = "Competitor Dominance Gradient Across Reefs",
    subtitle = expression(paste("Mean posterior ", P(Competitor > italic("M. hispida")),
                                " by reef (ordered highest to lowest macroalgae probability)")),
    x        = NULL,
    y        = expression(P(Competitor > Coral)),
    caption  = paste(
      "Error bars: range of site × habitat means within each reef.",
      "Dashed lines: 0.05 (refugia), 0.50 (uncertain), 0.95 (dominance) thresholds.",
      "Facets separate Inner Arc from Outer Arc reefs.",
      sep = "\n"
    )
  ) +
  theme(
    axis.text.x     = element_text(angle = 40, hjust = 1, size = 10),
    legend.key.size = unit(0.5, "cm")
  )

save_fig(fig11, "FIGURE_11_Dominance_by_Reef", 12, 7)

# ── 9. COMBINED PANEL ────────────────────────────────────────────────────────

cat("\n=== Combined panel ===\n")

combined <- (fig8 / fig9) +
  plot_annotation(
    tag_levels = "a", tag_prefix = "(", tag_suffix = ")",
    title    = "Benthic Community Competitive Dynamics – JSDM Posterior Analysis",
    subtitle = "Dirichlet JSDM with YEAR random effects (CV_ALL, WeaklyInformative prior)",
    caption  = "Model: cbind(MUSSISMILIA, TURF, CCA, CYANO, MACROALGAE) ~ PCA + ARCH interactions + (1|REEF) + (1|YEAR)"
  ) &
  theme(plot.tag = element_text(face = "bold", size = 14))

save_fig(combined, "FIGURE_COMBINED_Dominance_Analysis", 12, 15)

# ── 10. MANUSCRIPT TABLE ─────────────────────────────────────────────────────

cat("\n=== Manuscript table ===\n")

# Recompute with CCA included
tbl_df <- dominance_full %>%
  group_by(ARC, HAB_label) %>%
  summarise(
    `N sites`          = n_distinct(SITE),
    `N obs`            = n(),
    `P(Macro > Coral)` = sprintf("%.3f (%.3f–%.3f)",
                                  mean(P_macro_coral), min(P_macro_coral), max(P_macro_coral)),
    `P(CCA > Coral)`   = sprintf("%.3f (%.3f–%.3f)",
                                  mean(P_cca_coral),   min(P_cca_coral),   max(P_cca_coral)),
    `Coral (%)`        = sprintf("%.1f", mean(mean_MUSSISMILIA) * 100),
    `CCA (%)`          = sprintf("%.1f", mean(mean_CCA)         * 100),
    `Macroalgae (%)`   = sprintf("%.1f", mean(mean_MACROALGAE)  * 100),
    `Other benthic (%)` = sprintf("%.1f", mean(mean_other)       * 100),
    .groups = "drop"
  ) %>%
  rename(Arc = ARC, Habitat = HAB_label)

write.csv(tbl_df, file.path(output_dir, "TABLE_Dominance_Summary.csv"), row.names = FALSE)
cat("  Saved: TABLE_Dominance_Summary.csv\n")

# ── 11. FINAL SUMMARY ────────────────────────────────────────────────────────

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("All figures saved to: %s\n\n", output_dir))
cat("  FIGURE_8_Dominance_by_ARC_HAB     – P(Macro) and P(CCA) per ARC×HAB\n")
cat("  FIGURE_9_Dominance_Map             – Spatial map (fixed aspect + legend)\n")
cat("  FIGURE_10_Refugia_vs_Stress        – 4-panel composition (A:refugia B:stress C:mean D:env variability)\n")
cat("  FIGURE_11_Dominance_by_Reef        – Both competitors per reef\n")
cat("  FIGURE_COMBINED_Dominance_Analysis – Combined Fig8 + Fig9\n")
cat("  TABLE_Dominance_Summary.csv\n\n")

cat("CORAL REFUGIA (site-level, P_macro < 0.05):\n")
print(refugia_sites %>%
        select(REEF_NAME, SITE, HAB, ARC, P_macro_coral, P_cca_coral,
               mean_MUSSISMILIA, mean_MACROALGAE) %>%
        arrange(P_macro_coral) %>% head(8))

cat("\nSTRESS ZONES (site-level, P_macro > 0.95):\n")
print(stress_sites %>%
        select(REEF_NAME, SITE, HAB, ARC, P_macro_coral, P_cca_coral,
               mean_MUSSISMILIA, mean_MACROALGAE) %>%
        arrange(desc(P_macro_coral)) %>% head(8))

cat("\nScript 08 completed successfully!\n")
