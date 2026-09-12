# ============================================================================
# Boxplot_Cover_RGR_Publication.R
# Publication-quality raincloud plots for Cover (multi-organism) and RGR (M. hispida)
# Version: Opus 1.0 - Nature/Science Standards
#
# Generates:
#   1. Cover figure: 5 organisms (CCA, Cyano, Macroalgae, M. hispida, Turf)
#   2. RGR figure: M. hispida only
#
# Aesthetic standards match: PCA_Health_Interactions_Composite_GLM.R
# ============================================================================

# --- 1. PACKAGES ---
libs <- c(
    "readr", # Read CSV files
    "dplyr", # Data manipulation
    "tidyr", # Data reshaping
    "ggplot2", # Base plotting
    "ggdist", # stat_halfeye for raincloud plots
    "patchwork", # Panel composition
    "cowplot", # Legend extraction and assembly
    "scales", # Rescale functions
    "grid" # For graphical parameters
)

message("Loading packages...")
invisible(lapply(libs, function(lib) {
    if (!require(lib, character.only = TRUE)) {
        install.packages(lib)
        library(lib, character.only = TRUE)
    }
}))

# --- 2. PATHS ---
base_dir <- "."

# Input paths (corrected based on actual file locations)
cover_data_path <- "#######FINAL_DATA/05_benthic_cover/dados_integrados_long_format.csv"
rgr_data_path <- "#######FINAL_DATA/04_health_growth_PCA/resultados_biologicos_por_colonia.csv"

# Output directory with _Opus suffix
output_dir <- "#######FINAL_RESULTS/Boxplot_Publication_Figures_Opus"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

message(paste("Output directory:", output_dir))

# --- 3. VISUAL SCHEME ---

# Reef colors (consistent with PCA scripts)
reef_colors <- c(
    "ARC" = "#1f77b4", # Blue
    "ITA" = "#ff7f0e", # Orange
    "PAB" = "#2ca02c", # Green
    "UCR" = "#d62728", # Red
    "TIM" = "#9467bd" # Purple
)

# Habitat fill colors (complementary palette)
habitat_fill_colors <- c(
    "PA" = "#56B4E9", # Light blue (Reef Wall)
    "RR" = "#E69F00", # Orange (Rocky Reef)
    "TP" = "#009E73" # Teal (Reef Top)
)

# Habitat shapes (consistent with PCA scripts)
habitat_shapes <- c(
    "PA" = 21, # Filled Circle
    "RR" = 22, # Filled Square
    "TP" = 24 # Filled Triangle
)

# Organism display names (for plot labels)
organism_labels <- c(
    "CCA" = "CCA",
    "CYANO" = "Cyanobacteria",
    "MACROALGAE" = "Macroalgae",
    "MUSSISMILIA_HISPIDA" = expression(italic("M. hispida")),
    "TURF" = "Turf"
)

# Organism order for consistent display
organism_order <- c("CCA", "CYANO", "MACROALGAE", "MUSSISMILIA_HISPIDA", "TURF")

# Reef order (Inner Arc first, then Outer Arc)
reef_order <- c("ARC", "ITA", "PAB", "UCR", "TIM")

# Inner Arc reefs for styling distinction
inner_arc_reefs <- c("ARC", "ITA", "PAB")

# Theme publication (consistent with PCA scripts)
theme_publication <- theme_classic(base_size = 14.4) +
    theme(
        text = element_text(color = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 5)),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        panel.grid.major.y = element_line(
            color = "grey90",
            linetype = "dashed",
            linewidth = 0.3
        ),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(size = 12, face = "bold"),
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(5, 5, 5, 5)
    )

# --- 4. HELPER FUNCTIONS ---

load_cover_data <- function(file_path) {
    message("Loading cover data...")

    # Read first line to detect delimiter
    first_line <- readLines(file_path, n = 1)

    if (grepl(";", first_line)) {
        df <- read.csv2(file_path, stringsAsFactors = FALSE)
    } else {
        df <- read.csv(file_path, stringsAsFactors = FALSE)
    }

    # Standardize column names
    names(df) <- toupper(trimws(names(df)))

    # Debug: print column names
    message(sprintf("Columns found: %s", paste(names(df), collapse = ", ")))

    # Normalize HAB values (RO -> RR as per user clarification)
    if ("HAB" %in% names(df)) {
        df$HAB <- toupper(trimws(df$HAB))
        df$HAB[df$HAB == "RO"] <- "RR"
    } else {
        stop("HAB column not found in cover data!")
    }

    # Filter for target organisms
    target_organisms <- c("CCA", "CYANO", "MACROALGAE", "MUSSISMILIA_HISPIDA", "TURF")
    df <- df %>%
        filter(ORGANISMO %in% target_organisms) %>%
        mutate(
            ORGANISMO = factor(ORGANISMO, levels = organism_order),
            REEF = factor(REEF, levels = reef_order),
            HAB = factor(HAB),
            ARCH = ifelse(REEF %in% inner_arc_reefs, "Inner", "Outer"),
            COVER_PROP = COBERTURA # Keep as percentage for display
        )

    message(sprintf("Loaded %d observations for %d organisms", nrow(df), length(target_organisms)))
    message(sprintf("Reefs: %s", paste(unique(df$REEF), collapse = ", ")))
    message(sprintf("Habitats: %s", paste(unique(df$HAB), collapse = ", ")))

    return(df)
}

load_rgr_data <- function(file_path) {
    message("Loading RGR data...")

    # Read first line to detect delimiter
    first_line <- readLines(file_path, n = 1)

    if (grepl(";", first_line)) {
        df <- read.csv2(file_path, stringsAsFactors = FALSE)
    } else {
        df <- read.csv(file_path, stringsAsFactors = FALSE)
    }

    # Standardize column names
    names(df) <- toupper(trimws(names(df)))

    # Debug: print column names
    message(sprintf("RGR columns found: %s", paste(names(df), collapse = ", ")))

    # Normalize HAB values (RO -> RR as per user clarification)
    if ("HAB" %in% names(df)) {
        df$HAB <- toupper(trimws(df$HAB))
        df$HAB[df$HAB == "RO"] <- "RR"
    } else {
        stop("HAB column not found in RGR data!")
    }

    # Filter for M. hispida only (RGR is only for this species)
    df <- df %>%
        filter(!is.na(RGR)) %>%
        mutate(
            REEF = factor(REEF, levels = reef_order),
            HAB = factor(HAB),
            ARCH = ifelse(REEF %in% inner_arc_reefs, "Inner", "Outer")
        )

    message(sprintf("Loaded %d RGR observations", nrow(df)))
    message(sprintf("Reefs: %s", paste(unique(df$REEF), collapse = ", ")))
    message(sprintf("Habitats: %s", paste(unique(df$HAB), collapse = ", ")))

    return(df)
}

create_legend_strip <- function(habs, size_scale = 1.0) {
    legend_theme <- theme(
        legend.title = element_text(size = 14 * size_scale, face = "bold"),
        legend.text = element_text(size = 12 * size_scale),
        legend.key.size = unit(0.8 * size_scale, "cm"),
        legend.box.margin = margin(5, 0, 5, 0)
    )

    # Habitat Legend
    p_hab <- ggplot(data.frame(HAB = factor(habs, levels = habs)), aes(x = 1, y = 1, fill = HAB)) +
        geom_point(shape = 22, size = 6 * size_scale) +
        scale_fill_manual(values = habitat_fill_colors, name = "Habitat") +
        theme_publication +
        legend_theme
    leg_hab <- cowplot::get_legend(p_hab)

    # Arc Legend (Custom manual)
    p_arc <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
        annotate("rect",
            xmin = 0.8, xmax = 1.2, ymin = 0.8, ymax = 1.2,
            fill = "grey70", color = "black", linewidth = 1.2
        ) +
        annotate("text", x = 1.5, y = 1, label = "Inner Arc", hjust = 0, size = 4.5 * size_scale) +
        annotate("rect",
            xmin = 3.3, xmax = 3.7, ymin = 0.8, ymax = 1.2,
            fill = "grey70", color = "grey60", linewidth = 0.4
        ) +
        annotate("text", x = 4.0, y = 1, label = "Outer Arc", hjust = 0, size = 4.5 * size_scale) +
        xlim(0.5, 6) +
        ylim(0, 2) +
        theme_void() +
        labs(title = "Arc") +
        theme(
            plot.title = element_text(face = "bold", size = 14 * size_scale, hjust = 0.1, margin = margin(b = 5)),
            plot.margin = margin(10, 0, 0, 0)
        )

    legend_strip <- cowplot::plot_grid(
        leg_hab, p_arc,
        ncol = 2,
        rel_widths = c(1, 1.2),
        align = "vh"
    )

    return(legend_strip)
}

# --- 5. FIGURE GENERATION ---

create_cover_figure <- function(df, output_dir) {
    message("\n--- Creating Multi-Organism Cover Figure ---")

    # Custom labeller for organisms (with italics for M. hispida)
    organism_labeller <- function(labels) {
        lapply(labels, function(x) {
            if (x == "MUSSISMILIA_HISPIDA") {
                expression(italic("M. hispida"))
            } else if (x == "CYANO") {
                "Cyanobacteria"
            } else if (x == "MACROALGAE") {
                "Macroalgae"
            } else {
                as.character(x)
            }
        })
    }

    # Create the main plot with facet_grid
    p <- ggplot(df, aes(x = REEF, y = COVER_PROP, fill = HAB)) +
        # Half-eye (half violin) distribution
        stat_halfeye(
            aes(group = interaction(REEF, HAB)),
            adjust = 0.8,
            width = 0.5,
            .width = 0,
            justification = -0.2,
            point_colour = NA,
            alpha = 0.7,
            position = position_dodge(width = 0.8)
        ) +
        # Boxplot
        geom_boxplot(
            width = 0.15,
            outlier.shape = NA,
            alpha = 0.8,
            position = position_dodge(width = 0.8),
            color = "grey30"
        ) +
        # Jittered points
        geom_point(
            aes(shape = HAB),
            position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8),
            size = 1.5,
            alpha = 0.5,
            color = "grey30"
        ) +
        # Scales
        scale_fill_manual(values = habitat_fill_colors, name = "Habitat") +
        scale_shape_manual(values = c("PA" = 21, "RR" = 22, "TP" = 24), name = "Habitat") +
        # Facet by organism
        facet_wrap(
            ~ORGANISMO,
            ncol = 1,
            scales = "free_y",
            labeller = labeller(ORGANISMO = c(
                "CCA" = "CCA",
                "CYANO" = "Cyanobacteria",
                "MACROALGAE" = "Macroalgae",
                "MUSSISMILIA_HISPIDA" = "M. hispida",
                "TURF" = "Turf"
            ))
        ) +
        # Labels
        labs(
            x = "Reef",
            y = "Cover (%)"
        ) +
        theme_publication +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 0, hjust = 0.5),
            strip.text = element_text(face = "bold.italic"),
            panel.spacing = unit(0.8, "lines")
        )

    # Add Arc distinction via panel border styling
    # Inner Arc reefs will have thicker left border indication (visual grouping)

    # Create legend strip
    leg_strip <- create_legend_strip(levels(df$HAB))

    # Combine plot and legend
    final_combined <- cowplot::plot_grid(
        p,
        NULL,
        leg_strip,
        ncol = 1,
        rel_heights = c(15, 0.3, 1.5)
    )

    # Add title
    title <- cowplot::ggdraw() +
        cowplot::draw_label(
            "Benthic Cover by Reef and Habitat",
            fontface = "bold",
            size = 16,
            x = 0.5,
            hjust = 0.5
        )

    final_with_title <- cowplot::plot_grid(
        title,
        final_combined,
        ncol = 1,
        rel_heights = c(0.04, 1)
    )

    # Save PNG
    out_path_png <- file.path(output_dir, "Cover_MultiOrganism_Raincloud_Publication.png")
    ggsave(out_path_png, final_with_title, width = 10, height = 14, dpi = 300, bg = "white")
    message(paste("Saved:", out_path_png))

    # Save PDF (vector)
    out_path_pdf <- file.path(output_dir, "Cover_MultiOrganism_Raincloud_Publication.pdf")
    ggsave(out_path_pdf, final_with_title, width = 10, height = 14, device = cairo_pdf)
    message(paste("Saved:", out_path_pdf))

    return(invisible(NULL))
}

create_rgr_figure <- function(df, output_dir) {
    message("\n--- Creating RGR Figure (M. hispida Only) ---")

    # Create the main plot
    p <- ggplot(df, aes(x = HAB, y = RGR, fill = HAB)) +
        # Reference line at zero growth
        geom_hline(yintercept = 0, color = "#d62728", linetype = "dashed", linewidth = 0.8) +
        # Half-eye (half violin) distribution
        stat_halfeye(
            adjust = 0.8,
            width = 0.6,
            .width = 0,
            justification = -0.15,
            point_colour = NA,
            alpha = 0.7
        ) +
        # Boxplot
        geom_boxplot(
            width = 0.12,
            outlier.shape = NA,
            alpha = 0.8,
            color = "grey30"
        ) +
        # Jittered points
        geom_point(
            aes(shape = HAB),
            position = position_jitter(width = 0.05),
            size = 2,
            alpha = 0.6,
            color = "grey30"
        ) +
        # Scales
        scale_fill_manual(values = habitat_fill_colors, name = "Habitat") +
        scale_shape_manual(values = habitat_shapes, name = "Habitat") +
        # Facet by Reef
        facet_wrap(
            ~REEF,
            ncol = 5,
            labeller = labeller(REEF = function(x) {
                paste0(x, "\n(", ifelse(x %in% inner_arc_reefs, "Inner", "Outer"), ")")
            })
        ) +
        # Labels
        labs(
            x = "Habitat",
            y = expression(paste("Relative Growth Rate (", year^-1, ")"))
        ) +
        theme_publication +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 0, hjust = 0.5),
            strip.text = element_text(face = "bold"),
            panel.spacing = unit(1, "lines")
        )

    # Create legend strip
    leg_strip <- create_legend_strip(levels(df$HAB))

    # Combine plot and legend
    final_combined <- cowplot::plot_grid(
        p,
        NULL,
        leg_strip,
        ncol = 1,
        rel_heights = c(8, 0.3, 1.5)
    )

    # Add title
    title <- cowplot::ggdraw() +
        cowplot::draw_label(
            expression(paste("Relative Growth Rate of ", italic("Mussismilia hispida"), " by Reef and Habitat")),
            fontface = "bold",
            size = 16,
            x = 0.5,
            hjust = 0.5
        )

    final_with_title <- cowplot::plot_grid(
        title,
        final_combined,
        ncol = 1,
        rel_heights = c(0.06, 1)
    )

    # Save PNG
    out_path_png <- file.path(output_dir, "RGR_Mussismilia_hispida_Raincloud_Publication.png")
    ggsave(out_path_png, final_with_title, width = 14, height = 7, dpi = 300, bg = "white")
    message(paste("Saved:", out_path_png))

    # Save PDF (vector)
    out_path_pdf <- file.path(output_dir, "RGR_Mussismilia_hispida_Raincloud_Publication.pdf")
    ggsave(out_path_pdf, final_with_title, width = 14, height = 7, device = cairo_pdf)
    message(paste("Saved:", out_path_pdf))

    return(invisible(NULL))
}

# --- 6. MAIN EXECUTION ---

main <- function() {
    message("\n========================================")
    message("Boxplot Publication Figures Generator")
    message("Version: Opus 1.0")
    message("Standards: Nature / Science")
    message("========================================\n")

    # Verify input files exist
    if (!file.exists(cover_data_path)) {
        stop(paste("Error: Cover data file not found:", cover_data_path))
    }
    if (!file.exists(rgr_data_path)) {
        stop(paste("Error: RGR data file not found:", rgr_data_path))
    }

    # Load data
    cover_df <- load_cover_data(cover_data_path)
    rgr_df <- load_rgr_data(rgr_data_path)

    # Generate figures
    create_cover_figure(cover_df, output_dir)
    create_rgr_figure(rgr_df, output_dir)

    message("\n=== ALL FIGURES GENERATED SUCCESSFULLY ===")
    message(paste("Output directory:", output_dir))
    message("\nGenerated files:")
    message("  - Cover_MultiOrganism_Raincloud_Publication.png")
    message("  - Cover_MultiOrganism_Raincloud_Publication.pdf")
    message("  - RGR_Mussismilia_hispida_Raincloud_Publication.png")
    message("  - RGR_Mussismilia_hispida_Raincloud_Publication.pdf")
}

# Run main
main()
