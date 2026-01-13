# ============================================================================
# PCA_Visualization_Publication_vGeminiPro.R
# Geração de figuras PCA de alta qualidade para publicação
# Versão: Gemini Pro Refined V5 Hybrid (Footer fix: Increased Heights)
# ============================================================================

# --- 1. PACOTES E CONFIGURAÇÃO ---
libs <- c(
    "readxl", # Ler Excel
    "readr", # Ler CSV com controle de separador
    "dplyr", # Manipulação de dados
    "ggplot2", # Plotagem base
    "patchwork", # Composição de painéis ( CRUCIAL para alinhamento)
    "cowplot", # Extração de legendas e montagem final
    "scales", # Rescale de tamanhos
    "grid" # arrow() e unit() para loadings
)
invisible(lapply(libs, library, character.only = TRUE))

# --- Caminhos ---
base_dir <- "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_RESULTS"

output_dir <- file.path(base_dir, "PCA_Publication_Figures_vGemini_Pro")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 2. ESQUEMA VISUAL ---

reef_colors <- c(
    "ARC" = "#1f77b4", # Azul
    "ITA" = "#ff7f0e", # Laranja
    "PAB" = "#2ca02c", # Verde
    "UCR" = "#d62728", # Vermelho
    "TIM" = "#9467bd" # Roxo
)

habitat_shapes <- c(
    "PA" = 21, # Círculo
    "RR" = 22, # Quadrado
    "TP" = 24 # Triângulo
)

var_labels_var_default <- c("SST CV (%)", "DLI CV (%)", "Chl-a CV (%)")

loadings_labels <- c(
    # Magnitude
    "log(sst_mean+1)" = "log(SST+1)",
    "log(mean_DLI_local+1)" = "log(DLI+1)",
    "log(chl_mean+1)" = "log(Chl+1)",
    "sqrt(DHW>4)" = "sqrt(DHW>4)",
    # Variability - CV_02
    "sst_cv_2" = "SST CV",
    "dli_cv_2" = "DLI CV",
    "chl_cv_2" = "Chl CV",
    # Variability - CV_30
    "sst_cv_30" = "SST CV",
    "dli_cv_30" = "DLI CV",
    "chl_cv_30" = "Chl CV",
    # Variability - CV_ALL
    "sst_cv_all" = "SST CV",
    "cv_DLI_local" = "DLI CV",
    "chl_cv_all" = "Chl CV"
)

theme_publication <- theme_classic(base_size = 14.4) +
    theme(
        text = element_text(color = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 13, face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 11),
        panel.grid.major = element_line(
            color = "grey90",
            linetype = "dashed",
            linewidth = 0.3
        ),
        strip.background = element_blank(),
        plot.tag = element_text(face = "bold", size = 14)
    )

# --- 3. CONFIGURAÇÃO DOS CENÁRIOS ---

scenarios <- list(
    CV_02 = list(
        dir = file.path(base_dir, "#####output_local_PCA_CV_2_FINAL"),
        loadings_mag = "loadings_PCA_Magnitude_CV_2.csv",
        loadings_var = "loadings_PCA_Variability_CV_2.csv",
        csv_sep = ",", csv_dec = ".",
        pc1_mag = "PC1_Magnitude", pc2_mag = "PC2_Magnitude",
        pc1_var = "PC1_Variability", pc2_var = "PC2_Variability",
        var_mag_cols = c("sst_mean", "mean_DLI_local", "chl_mean", "prop_DHW_gt4"),
        var_mag_labels = c("SST (deg C)", "DLI (mol m^-2 d^-1)", "Chl-a (mg m^-3)", "DHW>4 (freq)"),
        var_var_cols = c("sst_cv_2", "dli_cv_2", "chl_cv_2"),
        var_var_labels = var_labels_var_default
    ),
    CV_30 = list(
        dir = file.path(base_dir, "#####output_local_PCA_CV_30_FINAL"),
        loadings_mag = "loadings_PCA_Magnitude_CV_30.csv",
        loadings_var = "loadings_PCA_Variability_CV_30.csv",
        csv_sep = ",", csv_dec = ".",
        pc1_mag = "PC1_Magnitude", pc2_mag = "PC2_Magnitude",
        pc1_var = "PC1_Variability", pc2_var = "PC2_Variability",
        var_mag_cols = c("sst_mean", "mean_DLI_local", "chl_mean", "prop_DHW_gt4"),
        var_mag_labels = c("SST (deg C)", "DLI (mol m^-2 d^-1)", "Chl-a (mg m^-3)", "DHW>4 (freq)"),
        var_var_cols = c("sst_cv_30", "dli_cv_30", "chl_cv_30"),
        var_var_labels = var_labels_var_default
    ),
    CV_ALL = list(
        dir = file.path(base_dir, "#####output_local_PCA_CV_all_FINAL"),
        loadings_mag = "loadings_PCA_Magnitude_CV_all.csv",
        loadings_var = "loadings_PCA_Variability_CV_all.csv",
        csv_sep = ";", csv_dec = ",",
        pc1_mag = "PC1_Magnitude", pc2_mag = "PC2_Magnitude",
        pc1_var = "PC1_Variability", pc2_var = "PC2_Variability",
        var_mag_cols = c("sst_mean", "mean_DLI_local", "chl_mean", "prop_DHW_gt4"),
        var_mag_labels = c("SST (deg C)", "DLI (mol m^-2 d^-1)", "Chl-a (mg m^-3)", "DHW>4 (freq)"),
        var_var_cols = c("sst_cv_all", "cv_DLI_local", "chl_cv_all"),
        var_var_labels = var_labels_var_default
    )
)

# --- 4. FUNÇÕES AUXILIARES ---

load_loadings <- function(file_path, sep, dec) {
    if (sep == ";") {
        df <- read.csv2(file_path, stringsAsFactors = FALSE)
    } else {
        df <- read.csv(file_path, stringsAsFactors = FALSE)
    }
    colnames(df)[1] <- "variable"
    df$PC1 <- as.numeric(gsub(",", ".", df$PC1))
    df$PC2 <- as.numeric(gsub(",", ".", df$PC2))
    return(df)
}

load_scenario_data <- function(scenario) {
    df <- readxl::read_excel(file.path(scenario$dir, "dados_consolidados_com_scores_das_duas_PCAs.xlsx"))
    df$pc1_mag <- df[[scenario$pc1_mag]]
    df$pc2_mag <- df[[scenario$pc2_mag]]
    df$pc1_var <- df[[scenario$pc1_var]]
    df$pc2_var <- df[[scenario$pc2_var]]
    if (is.null(df$pc1_mag)) {
        stop(sprintf("ERRO: Coluna %s nao encontrada no Excel!", scenario$pc1_mag))
    }
    return(df)
}

calculate_explained_variance <- function(df, var_cols) {
    data_clean <- na.omit(df[, var_cols, drop = FALSE])
    if (nrow(data_clean) < 3) {
        return(c(NA, NA))
    }
    pca_result <- prcomp(data_clean, scale. = TRUE)
    exp_var <- summary(pca_result)$importance[2, 1:2] * 100
    return(exp_var)
}

# --- 5. FUNÇÕES DE PLOTAGEM ---

create_bubble_panel <- function(df, var_col, pc1_col, pc2_col, exp_var, var_label) {
    var_values <- df[[var_col]]
    df$size_scaled <- scales::rescale(var_values, to = c(2, 8), from = range(var_values, na.rm = TRUE))

    border_color <- ifelse(df$Arch == "inner", "black", "grey60")
    border_width <- ifelse(df$Arch == "inner", 1.2, 0.4)

    vals <- c(df[[pc1_col]], df[[pc2_col]])
    limit <- max(abs(vals), na.rm = TRUE) * 1.1

    ggplot(df, aes(x = .data[[pc1_col]], y = .data[[pc2_col]])) +
        geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
        geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
        geom_point(aes(fill = Reef_name, shape = HAB, size = size_scaled),
            color = border_color, stroke = border_width, alpha = 0.8
        ) +
        scale_fill_manual(values = reef_colors, name = "Reef") +
        scale_shape_manual(values = habitat_shapes, name = "Habitat") +
        scale_size_identity() +
        labs(
            title = var_label,
            x = sprintf("PC1 (%.1f%%)", exp_var[1]),
            y = sprintf("PC2 (%.1f%%)", exp_var[2])
        ) +
        xlim(-limit, limit) +
        ylim(-limit, limit) +
        theme_publication +
        theme(legend.position = "none") +
        coord_fixed(ratio = 1)
}

create_loadings_panel <- function(loadings_df, exp_var, arrow_scale = 1.5) {
    loadings_df$xend <- loadings_df$PC1 * arrow_scale
    loadings_df$yend <- loadings_df$PC2 * arrow_scale
    loadings_df$label <- sapply(loadings_df$variable, function(v) {
        v_clean <- trimws(v)
        v_clean <- iconv(v_clean, to = "ASCII//TRANSLIT")
        if (v_clean %in% names(loadings_labels)) {
            return(loadings_labels[[v_clean]])
        }
        v_f <- gsub("^(sst_|dli_|chl_|cv_)", "", v_clean)
        v_f <- gsub("_", " ", v_f)
        return(v_f)
    })

    limit <- 2.2

    ggplot(loadings_df) +
        geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
        geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
        annotate("path",
            x = cos(seq(0, 2 * pi, length.out = 100)), y = sin(seq(0, 2 * pi, length.out = 100)),
            color = "grey80", linetype = "dashed", linewidth = 0.3
        ) +
        geom_segment(aes(x = 0, y = 0, xend = xend, yend = yend),
            arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
            color = "#d62728", linewidth = 0.7
        ) +
        geom_text(aes(x = xend * 1.2, y = yend * 1.2, label = label), size = 4, fontface = "bold") +
        labs(title = "Loadings", x = sprintf("PC1 (%.1f%%)", exp_var[1]), y = sprintf("PC2 (%.1f%%)", exp_var[2])) +
        xlim(-limit, limit) +
        ylim(-limit, limit) +
        theme_publication +
        coord_fixed(ratio = 1)
}

# --- 6. FUNÇÃO DE LEGENDA UNIFICADA (Horizontal) ---

create_full_legend_strip <- function(reefs, habs, size_scale = 1.0) {
    legend_theme <- theme(
        legend.title = element_text(size = 14 * size_scale, face = "bold"),
        legend.text = element_text(size = 12 * size_scale),
        legend.key.size = unit(1.0 * size_scale, "cm"),
        legend.box.margin = margin(0, 0, 0, 0)
    )

    p_reef <- ggplot(data.frame(Reef = factor(reefs, levels = reefs)), aes(x = 1, y = 1, fill = Reef)) +
        geom_point(shape = 21, size = 5 * size_scale) +
        scale_fill_manual(values = reef_colors, name = "Reef") +
        theme_publication +
        legend_theme
    leg_reef <- cowplot::get_legend(p_reef)

    p_hab <- ggplot(data.frame(HAB = factor(habs, levels = habs)), aes(x = 1, y = 1, shape = HAB)) +
        geom_point(size = 5 * size_scale, fill = "grey60") +
        scale_shape_manual(values = habitat_shapes, name = "Habitat") +
        theme_publication +
        legend_theme
    leg_hab <- cowplot::get_legend(p_hab)

    p_arch_manual <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
        annotate("point", x = 1, y = 1, shape = 21, size = 5 * size_scale, fill = "grey70", color = "black", stroke = 1.2) +
        annotate("text", x = 1.2, y = 1, label = "Inner Arc", hjust = 0, size = 4.5 * size_scale) +
        annotate("point", x = 3.5, y = 1, shape = 21, size = 5 * size_scale, fill = "grey70", color = "grey60", stroke = 0.4) +
        annotate("text", x = 3.7, y = 1, label = "Outer Arc", hjust = 0, size = 4.5 * size_scale) +
        xlim(0.8, 6) +
        ylim(0, 2) +
        theme_void() +
        labs(title = "Arc") +
        theme(
            plot.title = element_text(face = "bold", size = 14 * size_scale, hjust = 0.1, margin = margin(b = 5)),
            plot.margin = margin(10, 0, 0, 0)
        )

    legend_strip <- plot_grid(
        leg_reef,
        leg_hab,
        p_arch_manual,
        ncol = 3,
        rel_widths = c(1, 0.8, 1.2),
        align = "vh"
    )

    return(legend_strip)
}


# --- 7. ARQUITETURA DE MONTAGEM (Híbrida: Patchwork Plots + Cowplot Assembly) ---

create_magnitude_figure_v2 <- function(df, loadings_df, exp_var, scenario) {
    cat("    - Criando panels de Magnitude...\n")

    b1 <- create_bubble_panel(df, scenario$var_mag_cols[1], "pc1_mag", "pc2_mag", exp_var, scenario$var_mag_labels[1])
    b2 <- create_bubble_panel(df, scenario$var_mag_cols[2], "pc1_mag", "pc2_mag", exp_var, scenario$var_mag_labels[2])
    b3 <- create_bubble_panel(df, scenario$var_mag_cols[3], "pc1_mag", "pc2_mag", exp_var, scenario$var_mag_labels[3])
    b4 <- create_bubble_panel(df, scenario$var_mag_cols[4], "pc1_mag", "pc2_mag", exp_var, scenario$var_mag_labels[4])
    load_plot <- create_loadings_panel(loadings_df, exp_var)
    leg_strip <- create_full_legend_strip(unique(df$Reef_name), unique(df$HAB))

    # Grid de Plots via Patchwork (Alinhamento perfeito)
    layout_design <- "
    ABC
    DE#
    "

    # Adicionando um spacer (plot_spacer()) para preencher o slot #
    plots_grid <- b1 + b2 + b3 + b4 + load_plot + plot_spacer() +
        plot_layout(design = layout_design, widths = c(1, 1, 1)) +
        plot_annotation(tag_levels = "a", tag_suffix = ")")

    # Montagem Final Híbrida via Cowplot
    cat("    - Montando layout final (Híbrido)...\n")

    final_combined <- cowplot::plot_grid(
        plots_grid,
        NULL, # Espaçador Rígido Vertical
        leg_strip,
        ncol = 1,
        rel_heights = c(10, 1, 2.5), # Aumentado ratio da legenda
        align = "v"
    )

    return(final_combined)
}

create_variability_figure_v2 <- function(df, loadings_df, exp_var, scenario) {
    cat("    - Criando panels de Variabilidade...\n")

    b1 <- create_bubble_panel(df, scenario$var_var_cols[1], "pc1_var", "pc2_var", exp_var, scenario$var_var_labels[1])
    b2 <- create_bubble_panel(df, scenario$var_var_cols[2], "pc1_var", "pc2_var", exp_var, scenario$var_var_labels[2])
    b3 <- create_bubble_panel(df, scenario$var_var_cols[3], "pc1_var", "pc2_var", exp_var, scenario$var_var_labels[3])
    load_plot <- create_loadings_panel(loadings_df, exp_var)
    leg_strip <- create_full_legend_strip(unique(df$Reef_name), unique(df$HAB))

    layout_design <- "
    AB
    CD
    "

    plots_grid <- b1 + b2 + b3 + load_plot +
        plot_layout(design = layout_design, widths = c(1, 1)) +
        plot_annotation(tag_levels = "a", tag_suffix = ")")

    cat("    - Montando layout final (Híbrido)...\n")

    final_combined <- cowplot::plot_grid(
        plots_grid,
        NULL, # Espaçador Rígido Vertical
        leg_strip,
        ncol = 1,
        rel_heights = c(10, 1, 2.5), # Aumentado ratio da legenda
        align = "v"
    )

    return(final_combined)
}


# --- 8. LOOP PRINCIPAL ---

for (sn in names(scenarios)) {
    cat(sprintf("\n=== Processing scenario: %s ===\n", sn))
    sc <- scenarios[[sn]]

    cat("  - Loading data...\n")
    df <- load_scenario_data(sc)

    cat("  - Loading loadings...\n")
    loadings_mag <- load_loadings(file.path(sc$dir, sc$loadings_mag), sc$csv_sep, sc$csv_dec)
    loadings_var <- load_loadings(file.path(sc$dir, sc$loadings_var), sc$csv_sep, sc$csv_dec)

    cat("  - Calculating variance...\n")
    ev_mag <- calculate_explained_variance(df, sc$var_mag_cols)
    ev_var <- calculate_explained_variance(df, sc$var_var_cols)

    cat("  - Generating Magnitude figure...\n")
    fig_mag <- create_magnitude_figure_v2(df, loadings_mag, ev_mag, sc)
    mag_path <- file.path(output_dir, sprintf("Figure_PCA_Magnitude_%s", sn))
    ggsave(paste0(mag_path, ".png"), fig_mag, width = 14, height = 13, dpi = 300) # Fix: increased height

    cat("  - Generating Variability figure...\n")
    fig_var <- create_variability_figure_v2(df, loadings_var, ev_var, sc)
    var_path <- file.path(output_dir, sprintf("Figure_PCA_Variability_%s", sn))
    ggsave(paste0(var_path, ".png"), fig_var, width = 12, height = 13, dpi = 300) # Fix: increased height

    cat(sprintf("  ✓ %s finished\n", sn))
}

cat("\n=== PROCESSAMENTO CONCLUÍDO ===\n")
cat("Figuras salvas em:", output_dir, "\n")
