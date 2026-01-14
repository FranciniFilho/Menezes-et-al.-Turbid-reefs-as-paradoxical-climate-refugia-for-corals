# ============================================================================
# Boxplot_Publication_Grade_GLM.R
# Generation of publication-quality boxplots (Nature/Science standards)
# Version: 1.0 - Consistent with PCA figures
#
# FEATURES:
#   - Boxplot + stripplot with jitter (data transparency)
#   - Identical color scheme to PCA figures (reef_colors, habitat_shapes)
#   - theme_publication consistent across all figures
#   - Cowplot assembly for tight spacing
#   - PNG (300 dpi) + PDF (vector) export
# ============================================================================

# --- 1. PACKAGES ---
libs <- c(
    "readxl",       # Read Excel
    "readr",        # Read CSV
    "dplyr",        # Data manipulation
    "ggplot2",      # Base plotting
    "cowplot",      # Legend extraction and assembly
    "scales",       # Rescale functions
    "grid",         # unit() for dimensions
    "tools"         # toTitleCase()
    # ggbeeswarm is optional - we use geom_point with position_jitterdodge
)

message("Loading packages...")
invisible(lapply(libs, function(lib) {
    if (!require(lib, character.only = TRUE)) {
        install.packages(lib)
        library(lib, character.only = TRUE)
    }
}))

# --- 2. PATHS ---
base_dir <- "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_RESULTS"
output_dir <- file.path(base_dir, "BOXPLOTS_Publication_GLM")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Input data paths
coverage_data_path <- file.path(base_dir, "%%BOX_PLOTS_benthic_cover/dados_integrados_long_format.csv")
rgr_data_path <- file.path(base_dir, "#output_ANALISE_BIOLOGICA/resultados_biologicos_por_colonia.csv")

# --- 3. VISUAL SCHEME (Identical to PCA scripts) ---

# Reef colors (EXACT match with PCA scripts)
reef_colors <- c(
    "ARC" = "#1f77b4",  # Blue
    "ITA" = "#ff7f0e",  # Orange
    "PAB" = "#2ca02c",  # Green
    "UCR" = "#d62728",  # Red
    "TIM" = "#9467bd"   # Purple
)

# Habitat fill colors (distinct, publication-quality palette)
habitat_colors <- c(
    "PA" = "#E69F00",  # Orange/Gold
    "RR" = "#56B4E9",  # Sky Blue
    "TP" = "#009E73"   # Greenish
)

# Habitat shapes (EXACT match with PCA scripts)
habitat_shapes <- c(
    "PA" = 21,  # Filled Circle
    "RR" = 22,  # Filled Square
    "TP" = 24   # Filled Triangle
)

# Publication theme (EXACT match with PCA scripts)
theme_publication <- theme_classic(base_size = 14.4) +
    theme(
        text = element_text(color = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 3)),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        panel.grid.major = element_line(
            color = "grey90",
            linetype = "dashed",
            linewidth = 0.3
        ),
        strip.background = element_blank(),
        strip.text = element_text(size = 13, face = "bold"),
        plot.tag = element_text(face = "bold", size = 14),
        plot.margin = margin(5, 5, 5, 5)
    )

# --- 4. DATA LOADING FUNCTIONS ---

load_coverage_data <- function(file_path) {
    message("  Loading coverage data...")
    df <- read.csv(file_path, stringsAsFactors = FALSE)

    # Standardize column names
    df <- df %>%
        rename(
            ORGANISMO = ORGANISMO,
            REEF = REEF,
            HAB = HAB,
            COBERTURA = COBERTURA
        )

    # Standardize REEF and HAB
    df$REEF <- toupper(trimws(df$REEF))
    df$HAB <- toupper(trimws(df$HAB))

    # Standardize RO -> RR (Rocky Reef)
    df$HAB[df$HAB == "RO"] <- "RR"

    # Filter for complete cases
    df <- df[!is.na(df$COBERTURA), ]

    message(sprintf("    Loaded %d observations", nrow(df)))
    return(df)
}

load_rgr_data <- function(file_path) {
    message("  Loading RGR data...")
    df <- read.csv(file_path, stringsAsFactors = FALSE)

    # Standardize column names
    df <- df %>%
        rename(
            REEF = REEF,
            HAB = HAB,
            RGR = RGR
        )

    # Standardize REEF and HAB
    df$REEF <- toupper(trimws(df$REEF))
    df$HAB <- toupper(trimws(df$HAB))

    # Remove NA RGR values
    df <- df[!is.na(df$RGR), ]

    message(sprintf("    Loaded %d observations", nrow(df)))
    return(df)
}

# --- 5. PLOTTING FUNCTIONS ---

#' Create a boxplot + stripplot panel with publication quality
#'
#' @param df Data frame with the data
#' @param x_var Variable name for x-axis (categorical)
#' @param y_var Variable name for y-axis (continuous)
#' @param fill_var Variable name for fill color (categorical)
#' @param title Plot title
#' @param y_label Y-axis label
#' @param show_ref_line Whether to show horizontal reference line at y=0
#' @param y_limits Custom y-axis limits (NULL for auto)
#' @return ggplot object
create_boxplot_stripplot_panel <- function(df, x_var, y_var, fill_var,
                                           title = NULL, y_label = NULL,
                                           show_ref_line = FALSE,
                                           y_limits = NULL) {

    # Convert column names to symbols for aes()
    x_sym <- sym(x_var)
    y_sym <- sym(y_var)
    fill_sym <- sym(fill_var)

    # Determine y-axis limits
    if (is.null(y_limits)) {
        y_data <- df[[y_var]]
        y_range <- range(y_data, na.rm = TRUE)
        y_padding <- diff(y_range) * 0.1
        y_limits <- c(y_range[1] - y_padding, y_range[2] + y_padding)
    }

    # Build the plot
    p <- ggplot(df, aes(x = !!x_sym, y = !!y_sym, fill = !!fill_sym)) +
        # Reference line if requested
        (if (show_ref_line) geom_hline(yintercept = 0, color = "grey60", linetype = "dashed", linewidth = 0.5) else NULL) +
        # Boxplot (without outliers)
        geom_boxplot(
            outlier.shape = NA,
            width = 0.6,
            linewidth = 0.5,
            alpha = 0.9,
            color = "black"
        ) +
        # Stripplot with jitter
        geom_point(
            position = position_jitterdodge(
                jitter.width = 0.15,
                dodge.width = 0.6
            ),
            size = 1.5,
            alpha = 0.5,
            shape = 21,
            color = "black",
            stroke = 0.3
        ) +
        # Scales - use habitat_colors for fill
        scale_fill_manual(values = habitat_colors, name = "Habitat") +
        # Labels
        labs(
            title = title,
            x = "Reef",
            y = if (!is.null(y_label)) y_label else y_var
        ) +
        # Theme
        theme_publication +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1)
        ) +
        coord_cartesian(ylim = y_limits)

    return(p)
}

#' Create a faceted boxplot for multiple groups
#'
#' @param df Data frame with the data
#' @param facet_var Variable to facet by
#' @param x_var Variable for x-axis
#' @param y_var Variable for y-axis
#' @param fill_var Variable for fill color
#' @param y_label Y-axis label
#' @param ncol Number of columns in facet layout
#' @return ggplot object
create_faceted_boxplot <- function(df, facet_var, x_var, y_var, fill_var,
                                   y_label = NULL, ncol = 3) {

    facet_sym <- sym(facet_var)
    x_sym <- sym(x_var)
    y_sym <- sym(y_var)
    fill_sym <- sym(fill_var)

    # Order the facet variable
    df[[facet_var]] <- factor(df[[facet_var]], levels = sort(unique(df[[facet_var]])))

    # Build the plot
    p <- ggplot(df, aes(x = !!x_sym, y = !!y_sym, fill = !!fill_sym)) +
        # Boxplot (without outliers)
        geom_boxplot(
            outlier.shape = NA,
            width = 0.6,
            linewidth = 0.5,
            alpha = 0.9,
            color = "black"
        ) +
        # Stripplot with jitter
        geom_point(
            position = position_jitterdodge(
                jitter.width = 0.15,
                dodge.width = 0.6
            ),
            size = 1.2,
            alpha = 0.5,
            shape = 21,
            color = "black",
            stroke = 0.3
        ) +
        # Scales - use habitat_colors for fill
        scale_fill_manual(values = habitat_colors, name = "Habitat") +
        # Faceting
        facet_wrap(as.formula(paste("~", facet_var)), ncol = ncol, scales = "free_y") +
        # Labels
        labs(
            x = "Reef",
            y = if (!is.null(y_label)) y_label else y_var
        ) +
        # Theme
        theme_publication +
        theme(
            legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            strip.text = element_text(size = 11, face = "bold")
        )

    return(p)
}

#' Create unified legend strip for reef and habitat
#'
#' @param reefs Vector of reef names
#' @param habitats Vector of habitat names
#' @param size_scale Scaling factor for legend elements
#' @return Cowplot legend object
create_legend_strip <- function(reefs, habitats, size_scale = 1.0) {
    legend_theme <- theme(
        legend.title = element_text(size = 14 * size_scale, face = "bold"),
        legend.text = element_text(size = 12 * size_scale),
        legend.key.size = unit(1.0 * size_scale, "cm"),
        legend.box.margin = margin(10, 0, 10, 0)
    )

    # Reef Legend
    p_reef <- ggplot(data.frame(Reef = factor(reefs, levels = reefs)), aes(x = 1, y = 1, fill = Reef)) +
        geom_point(shape = 21, size = 5 * size_scale) +
        scale_fill_manual(values = reef_colors, name = "Reef") +
        theme_publication +
        legend_theme
    leg_reef <- cowplot::get_legend(p_reef)

    # Habitat Legend
    p_hab <- ggplot(data.frame(HAB = factor(habitats, levels = habitats)), aes(x = 1, y = 1, fill = HAB)) +
        geom_bar(stat = "identity", width = 0.5) +
        scale_fill_manual(values = c("PA" = "#E69F00", "RR" = "#56B4E9", "TP" = "#009E73"), name = "Habitat") +
        theme_publication +
        legend_theme
    leg_hab <- cowplot::get_legend(p_hab)

    # Combine legends horizontally
    legend_strip <- cowplot::plot_grid(
        leg_reef, leg_hab,
        ncol = 2,
        rel_widths = c(1, 0.8),
        align = "vh"
    )

    return(legend_strip)
}

# --- 6. FIGURE GENERATION FUNCTIONS ---

#' Generate coverage boxplot figure (faceted by organism)
#'
#' @param coverage_df Data frame with coverage data
#' @param output_dir Output directory path
generate_coverage_figure <- function(coverage_df, output_dir) {
    message("\n=== Generating Coverage Boxplot Figure ===")

    # Get unique organisms
    organisms <- unique(coverage_df$ORGANISMO)
    message(sprintf("  Found %d unique organisms", length(organisms)))

    # Create faceted plot
    p <- create_faceted_boxplot(
        df = coverage_df,
        facet_var = "ORGANISMO",
        x_var = "REEF",
        y_var = "COBERTURA",
        fill_var = "HAB",
        y_label = "Coverage (%)",
        ncol = 3
    )

    # Create legend
    leg_strip <- create_legend_strip(
        reefs = unique(coverage_df$REEF),
        habitats = unique(coverage_df$HAB)
    )

    # Combine plot with legend
    final_combined <- cowplot::plot_grid(
        p,
        NULL,
        leg_strip,
        ncol = 1,
        rel_heights = c(10, 0.5, 1.5),
        align = "v"
    )

    # Export PNG
    out_path_png <- file.path(output_dir, "Boxplot_Coverage_By_Reef_Habitat_GLM.png")
    ggsave(out_path_png, final_combined, width = 14, height = 12, dpi = 300, bg = "white")
    message(sprintf("  Saved PNG: %s", out_path_png))

    # Export PDF (vector)
    out_path_pdf <- file.path(output_dir, "Boxplot_Coverage_By_Reef_Habitat_GLM.pdf")
    ggsave(out_path_pdf, final_combined, width = 14, height = 12, device = cairo_pdf)
    message(sprintf("  Saved PDF: %s", out_path_pdf))

    message("  Coverage figure generation complete!")
    return(invisible(final_combined))
}

#' Generate RGR boxplot figure (single panel with reference line)
#'
#' @param rgr_df Data frame with RGR data
#' @param output_dir Output directory path
generate_rgr_figure <- function(rgr_df, output_dir) {
    message("\n=== Generating RGR Boxplot Figure ===")

    message(sprintf("  Found %d observations", nrow(rgr_df)))

    # Create panel
    p <- create_boxplot_stripplot_panel(
        df = rgr_df,
        x_var = "REEF",
        y_var = "RGR",
        fill_var = "HAB",
        title = "Relative Growth Rate (RGR)",
        y_label = "Relative Growth Rate (RGR)",
        show_ref_line = TRUE
    )

    # Create legend
    leg_strip <- create_legend_strip(
        reefs = unique(rgr_df$REEF),
        habitats = unique(rgr_df$HAB)
    )

    # Combine plot with legend
    final_combined <- cowplot::plot_grid(
        p,
        NULL,
        leg_strip,
        ncol = 1,
        rel_heights = c(10, 0.5, 1.5),
        align = "v"
    )

    # Export PNG
    out_path_png <- file.path(output_dir, "Boxplot_RGR_By_Reef_Habitat_GLM.png")
    ggsave(out_path_png, final_combined, width = 10, height = 8, dpi = 300, bg = "white")
    message(sprintf("  Saved PNG: %s", out_path_png))

    # Export PDF (vector)
    out_path_pdf <- file.path(output_dir, "Boxplot_RGR_By_Reef_Habitat_GLM.pdf")
    ggsave(out_path_pdf, final_combined, width = 10, height = 8, device = cairo_pdf)
    message(sprintf("  Saved PDF: %s", out_path_pdf))

    message("  RGR figure generation complete!")
    return(invisible(final_combined))
}

# --- 7. MAIN EXECUTION ---

main <- function() {
    message("\n========================================")
    message("Boxplot Publication Figures Generator")
    message("Standards: Nature / Science")
    message("========================================\n")

    # Load data
    message("Loading data...")
    coverage_df <- load_coverage_data(coverage_data_path)
    rgr_df <- load_rgr_data(rgr_data_path)

    # Generate coverage figure
    generate_coverage_figure(coverage_df, output_dir)

    # Generate RGR figure
    generate_rgr_figure(rgr_df, output_dir)

    message("\n=== ALL PROCESSES COMPLETED ===")
    message(sprintf("Output Directory: %s", output_dir))
}

# Run main function
main()
