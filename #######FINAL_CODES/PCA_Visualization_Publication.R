# ============================================================================
# PCA_Visualization_Publication.R
# Geração de figuras PCA de alta qualidade para publicação
# Projeto: Mussismilia hispida - Abrolhos
# Versão: 3.0 (Definitiva - Baseada no Plano v3)
# ============================================================================

# --- 1. PACOTES E CONFIGURAÇÃO ---
libs <- c(
    "readxl", # Ler Excel
    "readr", # Ler CSV com controle de separador
    "dplyr", # Manipulação de dados
    "ggplot2", # Plotagem base
    "patchwork", # Composição de painéis
    "cowplot", # Extração de legendas e plot_grid
    "scales", # Rescale de tamanhos
    "grid" # arrow() e unit() para loadings
)
invisible(lapply(libs, library, character.only = TRUE))

# --- Caminhos ---
# Ajuste o caminho base se necessário.
# Assume-se que o script está em #######FINAL_CODES/
base_dir <- "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_RESULTS"

output_dir <- file.path(base_dir, "PCA_Publication_Figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 2. ESQUEMA VISUAL ---

# Cores tab10 (ARC, ITA, PAB, UCR, TIM)
reef_colors <- c(
    "ARC" = "#1f77b4", # Azul
    "ITA" = "#ff7f0e", # Laranja
    "PAB" = "#2ca02c", # Verde
    "UCR" = "#d62728", # Vermelho
    "TIM" = "#9467bd" # Roxo
)

# Formas preenchíveis (21=círculo, 22=quadrado, 24=triângulo)
habitat_shapes <- c(
    "PA" = 21, # Círculo
    "RR" = 22, # Quadrado
    "TP" = 24 # Triângulo
)

# Labels para variabilidade (usados nos títulos dos painéis)
var_labels_var_default <- c("SST CV (%)", "DLI CV (%)", "Chl-a CV (%)")

# Labels para loadings (abreviados e limpos)
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

# Tema profissional
theme_publication <- theme_classic(base_size = 12) +
    theme(
        text = element_text(color = "black"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11, face = "bold"),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 9),
        panel.grid.major = element_line(
            color = "grey90",
            linetype = "dashed",
            linewidth = 0.3
        ),
        strip.background = element_blank(),
        plot.tag = element_text(face = "bold", size = 12)
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
        var_mag_labels = c("SST (°C)", "DLI (mol m⁻² d⁻¹)", "Chl-a (mg m⁻³)", "DHW>4 (freq)"),
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
        var_mag_labels = c("SST (°C)", "DLI (mol m⁻² d⁻¹)", "Chl-a (mg m⁻³)", "DHW>4 (freq)"),
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
        var_mag_labels = c("SST (°C)", "DLI (mol m⁻² d⁻¹)", "Chl-a (mg m⁻³)", "DHW>4 (freq)"),
        var_var_cols = c("sst_cv_all", "cv_DLI_local", "chl_cv_all"),
        var_var_labels = var_labels_var_default
    )
)

# --- 4. FUNÇÕES AUXILIARES ---

# Carregar loadings
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

# Carregar dados do cenário
load_scenario_data <- function(scenario) {
    df <- readxl::read_excel(file.path(scenario$dir, "dados_consolidados_com_scores_das_duas_PCAs.xlsx"))

    # Mapear scores PCA para nomes padronizados
    df$pc1_mag <- df[[scenario$pc1_mag]]
    df$pc2_mag <- df[[scenario$pc2_mag]]
    df$pc1_var <- df[[scenario$pc1_var]]
    df$pc2_var <- df[[scenario$pc2_var]]

    if (is.null(df$pc1_mag)) {
        stop(sprintf("ERRO: Coluna %s nao encontrada no Excel!", scenario$pc1_mag))
    }

    return(df)
}

# Calcular variância explicada
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

# Bubble Plot Individual
create_bubble_panel <- function(df, var_col, pc1_col, pc2_col, exp_var, var_label) {
    var_values <- df[[var_col]]
    df$size_scaled <- scales::rescale(var_values, to = c(2, 8), from = range(var_values, na.rm = TRUE))

    border_color <- ifelse(df$Arch == "inner", "black", "grey60")
    border_width <- ifelse(df$Arch == "inner", 1.2, 0.4) # Slightly adjusted for clarity

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
        theme_publication +
        theme(legend.position = "none") +
        coord_fixed(ratio = 1)
}

# Loadings Biplot
create_loadings_panel <- function(loadings_df, exp_var, arrow_scale = 1.5) {
    loadings_df$xend <- loadings_df$PC1 * arrow_scale
    loadings_df$yend <- loadings_df$PC2 * arrow_scale
    loadings_df$label <- sapply(loadings_df$variable, function(v) {
        v_clean <- trimws(v)
        if (v_clean %in% names(loadings_labels)) {
            return(loadings_labels[[v_clean]])
        }
        # Fallback names
        v_f <- gsub("^(sst_|dli_|chl_|cv_)", "", v_clean)
        v_f <- gsub("_", " ", v_f)
        return(v_f)
    })

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
        geom_text(aes(x = xend * 1.2, y = yend * 1.2, label = label), size = 3, fontface = "bold") +
        labs(title = "Loadings", x = sprintf("PC1 (%.1f%%)", exp_var[1]), y = sprintf("PC2 (%.1f%%)", exp_var[2])) +
        theme_publication +
        coord_fixed(ratio = 1, xlim = c(-2.2, 2.2), ylim = c(-2.2, 2.2))
}

# Painel de Legenda (Refatorado para patchwork)
create_legend_panel <- function(reefs, habs) {
    # Legenda de recifes (cores)
    p_reef <- ggplot(data.frame(Reef = factor(reefs, levels = reefs)), aes(x = 1, y = seq_along(Reef), fill = Reef)) +
        geom_point(shape = 21, size = 3.5) +
        scale_fill_manual(values = reef_colors, name = "Reef") +
        theme_publication +
        theme(legend.position = "right")

    # Legenda de habitats (formas)
    p_hab <- ggplot(data.frame(HAB = factor(habs, levels = habs)), aes(x = 1, y = seq_along(HAB), shape = HAB)) +
        geom_point(size = 3.5, fill = "grey60") +
        scale_shape_manual(values = habitat_shapes, name = "Habitat") +
        theme_publication +
        theme(legend.position = "right")

    # Legenda manual para Arch
    p_arch <- ggplot(data.frame(x = 1, y = 2:1, label = c("Inner Arc", "Outer Arc"), s = c(1.2, 0.4), c = c("black", "grey60")), aes(x, y)) +
        geom_point(shape = 21, size = 3.5, fill = "grey70", color = c("black", "grey60"), stroke = c(1.2, 0.4)) +
        geom_text(aes(x + 1.0, y, label = label), hjust = 0, size = 3) +
        xlim(0.5, 4.0) +
        ylim(0.5, 2.5) +
        labs(title = "Arc") +
        theme_void() +
        theme(
            plot.title = element_text(face = "bold", size = 10, hjust = 0),
            plot.margin = margin(5, 5, 5, 5)
        )

    # Função interna para extrair legenda de forma extremamente robusta
    get_leg <- function(p) {
        cat("      - Extraindo componente de legenda...\n")
        # Tentar via cowplot::get_plot_component com return_all=TRUE para evitar erro 3.5.0+
        leg <- cowplot::get_plot_component(p, "guide-box", return_all = TRUE)
        if (is.list(leg) && length(leg) > 0) {
            leg <- leg[[1]]
        }
        if (is.null(leg)) {
            # Tentar cowplot::get_legend antigo
            leg <- cowplot::get_legend(p)
        }
        return(leg)
    }

    cat("      - Preparando leg_reef...\n")
    leg_reef <- get_leg(p_reef)
    cat("      - Preparando leg_hab...\n")
    leg_hab <- get_leg(p_hab)

    cat("      - Combinando legendas...\n")
    # Combinar usando patchwork
    res <- wrap_elements(full = leg_reef) /
        wrap_elements(full = leg_hab) /
        wrap_elements(full = p_arch) +
        plot_layout(heights = c(1.2, 0.8, 1))

    return(res)
}

# --- 6. FUNÇÕES DE COMPOSIÇÃO ---

create_magnitude_figure <- function(df, loadings_df, exp_var, scenario) {
    cat("    - Criando bubbles...\n")
    bubbles <- lapply(seq_along(scenario$var_mag_cols), function(i) {
        create_bubble_panel(df, scenario$var_mag_cols[i], "pc1_mag", "pc2_mag", exp_var, scenario$var_mag_labels[i])
    })

    cat("    - Criando loadings...\n")
    p_loadings <- create_loadings_panel(loadings_df, exp_var)

    cat("    - Criando legenda...\n")
    p_legend <- create_legend_panel(unique(df$Reef_name), unique(df$HAB))

    cat("    - Agrupando...\n")
    # Grid 2x3
    fig <- (bubbles[[1]] | bubbles[[2]] | bubbles[[3]]) /
        (bubbles[[4]] | p_loadings | p_legend) +
        plot_annotation(tag_levels = "a", tag_suffix = ")")

    return(fig)
}

create_variability_figure <- function(df, loadings_df, exp_var, scenario) {
    cat("    - Criando bubbles...\n")
    bubbles <- lapply(seq_along(scenario$var_var_cols), function(i) {
        create_bubble_panel(df, scenario$var_var_cols[i], "pc1_var", "pc2_var", exp_var, scenario$var_var_labels[i])
    })

    cat("    - Criando loadings...\n")
    p_loadings <- create_loadings_panel(loadings_df, exp_var)

    cat("    - Criando legenda...\n")
    p_legend <- create_legend_panel(unique(df$Reef_name), unique(df$HAB))

    cat("    - Agrupando...\n")
    # Grid principal 2x2
    main_grid <- (bubbles[[1]] | bubbles[[2]]) / (bubbles[[3]] | p_loadings) +
        plot_annotation(tag_levels = "a", tag_suffix = ")")

    # Compor com legenda à direita
    final_fig <- main_grid | wrap_elements(full = p_legend)
    final_fig <- final_fig + plot_layout(widths = c(4, 1))

    return(final_fig)
}

# --- 7. LOOP PRINCIPAL ---

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
    fig_mag <- create_magnitude_figure(df, loadings_mag, ev_mag, sc)
    mag_path <- file.path(output_dir, sprintf("Figure_PCA_Magnitude_%s", sn))
    ggsave(paste0(mag_path, ".png"), fig_mag, width = 14, height = 10, dpi = 300)
    ggsave(paste0(mag_path, ".pdf"), fig_mag, width = 14, height = 10)

    cat("  - Generating Variability figure...\n")
    fig_var <- create_variability_figure(df, loadings_var, ev_var, sc)
    var_path <- file.path(output_dir, sprintf("Figure_PCA_Variability_%s", sn))
    ggsave(paste0(var_path, ".png"), fig_var, width = 12, height = 10, dpi = 300)
    ggsave(paste0(var_path, ".pdf"), fig_var, width = 12, height = 10)

    cat(sprintf("  ✓ %s finished\n", sn))
}

cat("\n=== PROCESSAMENTO CONCLUÍDO ===\n")
cat("Figuras salvas em:", output_dir, "\n")
