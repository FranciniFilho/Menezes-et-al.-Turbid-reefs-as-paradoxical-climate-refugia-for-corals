# ============================================================================
# PCA_Environ_v2.R
# Geração de figura PCA de alta qualidade para publicação
# Versão: v2 - Figura única combinada (Magnitude + Variabilidade)
# - 4 linhas x 2 colunas
# - Linha 1: Magnitude (PCA | Loadings)
# - Linhas 2-4: Variabilidade CV_2, CV_30, CV_all (PCA | Loadings)
# - Abandona bubble plots: tamanhos de amostra constantes
# ============================================================================

# --- 1. PACOTES E CONFIGURAÇÃO ---
libs <- c(
    "readxl",
    "readr",
    "dplyr",
    "ggplot2",
    "patchwork",
    "cowplot",
    "grid"
)
invisible(lapply(libs, library, character.only = TRUE))

# --- Caminhos ---
base_dir <- "#######FINAL_RESULTS"

output_dir <- file.path(base_dir, "PCA_Publication_Figures_v2")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 2. ESQUEMA VISUAL ---

reef_colors <- c(
    "ARC" = "#1f77b4",
    "ITA" = "#ff7f0e",
    "PAB" = "#2ca02c",
    "UCR" = "#d62728",
    "TIM" = "#9467bd"
)

habitat_shapes <- c(
    "PA" = 21,
    "RR" = 22,
    "TP" = 24
)

var_labels_var <- c("SST CV (%)", "DLI CV (%)", "Chl-a CV (%)")

loadings_labels <- c(
    "log(sst_mean+1)" = "log(SST+1)",
    "log(mean_DLI_local+1)" = "log(DLI+1)",
    "log(chl_mean+1)" = "log(Chl+1)",
    "sqrt(DHW>4)" = "sqrt(DHW>4)",
    "sst_cv_2" = "SST CV",
    "dli_cv_2" = "DLI CV",
    "chl_cv_2" = "Chl CV",
    "sst_cv_30" = "SST CV",
    "dli_cv_30" = "DLI CV",
    "chl_cv_30" = "Chl CV",
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

labels_mag_formatted <- list(
    expression(paste("SST (", degree, "C)")),
    expression(paste("DLI (mol ", m^-2, " ", d^-1, ")")),
    expression(paste("Chl-a (mg ", m^-3, ")")),
    "DHW>4 (freq)"
)

scenarios <- list(
    CV_02 = list(
        dir = file.path(base_dir, "#####output_local_PCA_CV_2_FINAL"),
        loadings_mag = "loadings_PCA_Magnitude_CV_2.csv",
        loadings_var = "loadings_PCA_Variability_CV_2.csv",
        csv_sep = ",", csv_dec = ".",
        pc1_mag = "PC1_Magnitude", pc2_mag = "PC2_Magnitude",
        pc1_var = "PC1_Variability", pc2_var = "PC2_Variability",
        var_mag_cols = c("sst_mean", "mean_DLI_local", "chl_mean", "prop_DHW_gt4"),
        var_mag_labels = labels_mag_formatted,
        var_var_cols = c("sst_cv_2", "dli_cv_2", "chl_cv_2"),
        var_var_labels = var_labels_var,
        label = "CV_2 window"
    ),
    CV_30 = list(
        dir = file.path(base_dir, "#####output_local_PCA_CV_30_FINAL"),
        loadings_mag = "loadings_PCA_Magnitude_CV_30.csv",
        loadings_var = "loadings_PCA_Variability_CV_30.csv",
        csv_sep = ",", csv_dec = ".",
        pc1_mag = "PC1_Magnitude", pc2_mag = "PC2_Magnitude",
        pc1_var = "PC1_Variability", pc2_var = "PC2_Variability",
        var_mag_cols = c("sst_mean", "mean_DLI_local", "chl_mean", "prop_DHW_gt4"),
        var_mag_labels = labels_mag_formatted,
        var_var_cols = c("sst_cv_30", "dli_cv_30", "chl_cv_30"),
        var_var_labels = var_labels_var,
        label = "CV_30 window"
    ),
    CV_ALL = list(
        dir = file.path(base_dir, "#####output_local_PCA_CV_all_FINAL"),
        loadings_mag = "loadings_PCA_Magnitude_CV_all.csv",
        loadings_var = "loadings_PCA_Variability_CV_all.csv",
        csv_sep = ";", csv_dec = ",",
        pc1_mag = "PC1_Magnitude", pc2_mag = "PC2_Magnitude",
        pc1_var = "PC1_Variability", pc2_var = "PC2_Variability",
        var_mag_cols = c("sst_mean", "mean_DLI_local", "chl_mean", "prop_DHW_gt4"),
        var_mag_labels = labels_mag_formatted,
        var_var_cols = c("sst_cv_all", "cv_DLI_local", "chl_cv_all"),
        var_var_labels = var_labels_var,
        label = "CV_all window"
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

# --- 5. FUNÇÕES DE PLOTAGEM (SEM BUBBLE) ---

create_pca_panel <- function(df, pc1_col, pc2_col, exp_var, title = "", point_size = 4) {
    
    border_color <- ifelse(df$Arch == "inner", "black", "grey60")
    border_width <- ifelse(df$Arch == "inner", 1.2, 0.4)
    
    vals <- c(df[[pc1_col]], df[[pc2_col]])
    limit <- max(abs(vals), na.rm = TRUE) * 1.15
    
    ggplot(df, aes(x = .data[[pc1_col]], y = .data[[pc2_col]])) +
        geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
        geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
        geom_point(
            aes(fill = Reef_name, shape = HAB),
            color = border_color, 
            stroke = border_width, 
            size = point_size,
            alpha = 0.85
        ) +
        scale_fill_manual(values = reef_colors, name = "Reef") +
        scale_shape_manual(values = habitat_shapes, name = "Habitat") +
        labs(
            title = title,
            x = sprintf("PC1 (%.1f%%)", exp_var[1]),
            y = sprintf("PC2 (%.1f%%)", exp_var[2])
        ) +
        xlim(-limit, limit) +
        ylim(-limit, limit) +
        theme_publication +
        theme(legend.position = "none") +
        coord_fixed(ratio = 1)
}

create_loadings_panel <- function(loadings_df, exp_var, arrow_scale = 1.5, title = "Loadings") {
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
            x = cos(seq(0, 2 * pi, length.out = 100)), 
            y = sin(seq(0, 2 * pi, length.out = 100)),
            color = "grey80", linetype = "dashed", linewidth = 0.3
        ) +
        geom_segment(
            aes(x = 0, y = 0, xend = xend, yend = yend),
            arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
            color = "#d62728", linewidth = 0.7
        ) +
        geom_text(aes(x = xend * 1.2, y = yend * 1.2, label = label), 
                  size = 4, fontface = "bold") +
        labs(
            title = title, 
            x = sprintf("PC1 (%.1f%%)", exp_var[1]), 
            y = sprintf("PC2 (%.1f%%)", exp_var[2])
        ) +
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
    
    p_reef <- ggplot(data.frame(Reef = factor(reefs, levels = reefs)), 
                     aes(x = 1, y = 1, fill = Reef)) +
        geom_point(shape = 21, size = 5 * size_scale) +
        scale_fill_manual(values = reef_colors, name = "Reef") +
        theme_publication +
        legend_theme
    leg_reef <- cowplot::get_legend(p_reef)
    
    p_hab <- ggplot(data.frame(HAB = factor(habs, levels = habs)), 
                    aes(x = 1, y = 1, shape = HAB)) +
        geom_point(size = 5 * size_scale, fill = "grey60") +
        scale_shape_manual(values = habitat_shapes, name = "Habitat") +
        theme_publication +
        legend_theme
    leg_hab <- cowplot::get_legend(p_hab)
    
    p_arch_manual <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
        annotate("point", x = 1, y = 1, shape = 21, size = 5 * size_scale, 
                 fill = "grey70", color = "black", stroke = 1.2) +
        annotate("text", x = 1.2, y = 1, label = "Inner Arc", hjust = 0, 
                 size = 4.5 * size_scale) +
        annotate("point", x = 3.5, y = 1, shape = 21, size = 5 * size_scale, 
                 fill = "grey70", color = "grey60", stroke = 0.4) +
        annotate("text", x = 3.7, y = 1, label = "Outer Arc", hjust = 0, 
                 size = 4.5 * size_scale) +
        xlim(0.8, 6) +
        ylim(0, 2) +
        theme_void() +
        labs(title = "Arc") +
        theme(
            plot.title = element_text(face = "bold", size = 14 * size_scale, 
                                      hjust = 0.1, margin = margin(b = 5)),
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

# --- 7. FIGURA COMBINADA (MAGNITUDE + VARIABILIDADE) ---

create_combined_figure <- function(df_mag, loadings_mag, ev_mag,
                                    all_scenario_data, all_loadings_var, all_exp_var_var) {
    cat("    - Criando figura combinada (4 linhas x 2 colunas)...\n")
    
    panels_list <- list()
    
    # --- Linha 1: Magnitude ---
    pca_mag <- create_pca_panel(df_mag, "pc1_mag", "pc2_mag", ev_mag, 
                                title = "Environmental Magnitude", point_size = 4)
    loadings_mag_plot <- create_loadings_panel(loadings_mag, ev_mag, 
                                                title = "Loadings (Magnitude)")
    panels_list[[1]] <- pca_mag
    panels_list[[2]] <- loadings_mag_plot
    
    # --- Linhas 2-4: Variabilidade ---
    tag_counter <- 3
    for (sn in names(all_scenario_data)) {
        df <- all_scenario_data[[sn]]$data
        loadings <- all_loadings_var[[sn]]
        exp_var <- all_exp_var_var[[sn]]
        label <- all_scenario_data[[sn]]$label
        pc1_col <- all_scenario_data[[sn]]$pc1_col
        pc2_col <- all_scenario_data[[sn]]$pc2_col
        
        pca_plot <- create_pca_panel(df, pc1_col, pc2_col, exp_var, 
                                      title = label, point_size = 3.5)
        
        loadings_plot <- create_loadings_panel(loadings, exp_var, 
                                                title = sprintf("Loadings (%s)", label))
        
        panels_list[[tag_counter]] <- pca_plot
        panels_list[[tag_counter + 1]] <- loadings_plot
        tag_counter <- tag_counter + 2
    }
    
    # --- Legenda ---
    leg_strip <- create_full_legend_strip(unique(df_mag$Reef_name), 
                                           unique(df_mag$HAB),
                                           size_scale = 0.9)
    
    # --- Montar grid ---
    plots_combined <- wrap_plots(panels_list, ncol = 2, byrow = TRUE) +
        plot_annotation(tag_levels = "a", tag_suffix = ")") &
        theme(plot.tag = element_text(face = "bold", size = 12))
    
    cat("    - Montando layout final...\n")
    
    final_combined <- cowplot::plot_grid(
        plots_combined,
        NULL,
        leg_strip,
        ncol = 1,
        rel_heights = c(20, 0.3, 2),
        align = "v"
    )
    
    return(final_combined)
}

# --- 8. EXECUÇÃO PRINCIPAL ---

cat("\n============================================================================\n")
cat("PCA_Environ_v2.R - Processamento Iniciado\n")
cat("============================================================================\n")

# --- 8.1 Carregar dados de todos os cenários ---

all_scenario_data <- list()
all_loadings_var <- list()
all_exp_var_var <- list()

for (sn in names(scenarios)) {
    cat(sprintf("\n  Carregando dados de %s...\n", sn))
    sc <- scenarios[[sn]]
    
    df <- load_scenario_data(sc)
    loadings_var <- load_loadings(file.path(sc$dir, sc$loadings_var), sc$csv_sep, sc$csv_dec)
    ev_var <- calculate_explained_variance(df, sc$var_var_cols)
    
    all_scenario_data[[sn]] <- list(
        data = df,
        label = sc$label,
        pc1_col = "pc1_var",
        pc2_col = "pc2_var"
    )
    all_loadings_var[[sn]] <- loadings_var
    all_exp_var_var[[sn]] <- ev_var
    
    cat(sprintf("  ✓ %s carregado (PC1: %.1f%%, PC2: %.1f%%)\n", sn, ev_var[1], ev_var[2]))
}

# --- 8.2 Carregar dados de Magnitude (usando CV_02 como referência) ---
# NOTA: Os scores e loadings de Magnitude de CV_02 têm sinais invertidos em
# relação a CV_ALL (causa: force_positive_pc1="SST" aplicado inconsistentemente).
# O JSDM vencedor (CV_ALL) usa a convenção CV_ALL, portanto ambos os eixos
# (PC1 e PC2) são explicitamente invertidos aqui para harmonizar a figura com
# o texto do manuscrito e com os coeficientes do modelo.

cat("\n------------------------------------------------------------------------\n")
cat("Carregando dados de Magnitude (CV_02, com inversão de sinal para convenção CV_ALL)...\n")
cat("------------------------------------------------------------------------\n")

sc_mag <- scenarios[["CV_02"]]
df_mag <- all_scenario_data[["CV_02"]]$data
loadings_mag <- load_loadings(file.path(sc_mag$dir, sc_mag$loadings_mag), 
                               sc_mag$csv_sep, sc_mag$csv_dec)
ev_mag <- calculate_explained_variance(df_mag, sc_mag$var_mag_cols)

# Inverter ambos os eixos de Magnitude para convenção CV_ALL
# (PC1 CV_ALL: SST negativo, Chl positivo; PC2 CV_ALL: DLI positivo)
df_mag$pc1_mag <- -df_mag$pc1_mag
df_mag$pc2_mag <- -df_mag$pc2_mag
loadings_mag$PC1 <- -loadings_mag$PC1
loadings_mag$PC2 <- -loadings_mag$PC2

cat(sprintf("  Variância explicada: PC1 = %.1f%%, PC2 = %.1f%%\n", ev_mag[1], ev_mag[2]))

# --- 8.3 Gerar Figura Combinada ---

cat("\n------------------------------------------------------------------------\n")
cat("Gerando Figura Combinada (Magnitude + Variabilidade)...\n")
cat("------------------------------------------------------------------------\n")

fig_combined <- create_combined_figure(df_mag, loadings_mag, ev_mag,
                                        all_scenario_data, all_loadings_var, all_exp_var_var)

combined_path <- file.path(output_dir, "Figure_PCA_Combined_v2")
ggsave(paste0(combined_path, ".png"), fig_combined, width = 12, height = 22, dpi = 300)
ggsave(paste0(combined_path, ".pdf"), fig_combined, width = 12, height = 22)

cat(sprintf("  ✓ Figura salva: %s\n", paste0(combined_path, ".png")))

# --- 9. RESUMO FINAL ---

cat("\n============================================================================\n")
cat("PROCESSAMENTO CONCLUÍDO\n")
cat("============================================================================\n")
cat(sprintf("Diretório de saída: %s\n", output_dir))
cat("Arquivo gerado:\n")
cat("  - Figure_PCA_Combined_v2.png/.pdf\n")
cat("============================================================================\n")