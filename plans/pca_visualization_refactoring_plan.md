# Plano Detalhado: Código de Visualização PCA

## Resumo das Decisões

| Item | Decisão |
|------|---------|
| **Cenários** | Gerar figuras para CV_02, CV_30 e CV_all (Opção A) |
| **Linguagem** | R (ggplot2 + patchwork + cowplot) |
| **Nomes de variáveis** | Abreviações técnicas profissionais |
| **Variância explicada** | Calculada a partir dos loadings |

---

## Estrutura dos Dados de Entrada

### Excel: `dados_consolidados_com_scores_das_duas_PCAs.xlsx`

| Coluna | Descrição | Uso |
|--------|-----------|-----|
| `Reef_name` | Nome do recife (ARC, ITA, PAB, UCR, TIM) | Cor dos pontos |
| `HAB` | Habitat (RR, TP, PA) | Forma dos pontos |
| `Arch` | Arco (inner, outer) | Estilo da borda |
| `PC1_Variability`, `PC2_Variability` | Scores PCA Variabilidade | Eixos X e Y |
| `PC1_Magnitude`, `PC2_Magnitude` | Scores PCA Magnitude | Eixos X e Y |
| `sst_mean`, `chl_mean`, `mean_DLI_local`, `prop_DHW_gt4` | Variáveis de Magnitude | Tamanho dos pontos |
| `sst_cv_*`, `chl_cv_*`, `dli_cv_*` | Variáveis de Variabilidade | Tamanho dos pontos |

### CSVs de Loadings

> **ATENÇÃO**: Separadores inconsistentes entre cenários!
> - CV_02 e CV_30: `,` (vírgula) como separador, `.` decimal
> - CV_all: `;` (ponto e vírgula) como separador, `,` decimal

| Arquivo | Variáveis |
|---------|-----------|
| `loadings_PCA_Magnitude_CV_*.csv` | `log(sst_mean+1)`, `log(mean_DLI_local+1)`, `log(chl_mean+1)`, `sqrt(DHW>4)` |
| `loadings_PCA_Variability_CV_*.csv` | `sst_cv_*`, `dli_cv_*`, `chl_cv_*` |

---

## Nomenclatura Técnica para Labels

| Variável Original | Label Técnico (Inglês) |
|-------------------|------------------------|
| `sst_mean` | SST (°C) |
| `mean_DLI_local` | DLI (mol m⁻² d⁻¹) |
| `chl_mean` | Chl-a (mg m⁻³) |
| `prop_DHW_gt4` | DHW>4 (freq) |
| `sst_cv_*` | SST CV (%) |
| `dli_cv_*` / `cv_DLI_local` | DLI CV (%) |
| `chl_cv_*` | Chl-a CV (%) |

---

## Cálculo da Variância Explicada

A variância explicada não está salva nos arquivos. Pode ser:

1. **Recalculada via PCA** (mais preciso - recomendado)
2. **Hardcoded** baseado nos valores originais do script Python

Para a opção 1, usamos os dados brutos do Excel:

```r
# Recalcular variância explicada
recalculate_explained_variance <- function(df, vars, n_components = 2) {
  data_clean <- df[, vars] %>% na.omit()
  pca_result <- prcomp(data_clean, scale. = TRUE)
  var_exp <- summary(pca_result)$importance[2, 1:n_components] * 100
  return(var_exp)
}
```

---

## Estrutura do Script R

```r
# ============================================================================
# PCA_Visualization_Publication.R
# Script independente para gerar figuras de PCA de alta qualidade
# Autor: [Seu nome]
# Data: 2026-01-12
# ============================================================================

# --- 1. PACOTES E CONFIGURAÇÃO ---
libs <- c("readxl", "readr", "dplyr", "tidyr", "ggplot2", 
          "patchwork", "cowplot", "RColorBrewer")
invisible(lapply(libs, library, character.only = TRUE))

# --- 2. CAMINHOS E CONFIGURAÇÕES ---
base_dir <- "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_RESULTS"

scenarios <- list(
  CV_02 = list(
    dir = file.path(base_dir, "#####output_local_PCA_CV_2_FINAL"),
    scores_file = "dados_consolidados_com_scores_das_duas_PCAs.xlsx",
    loadings_mag = "loadings_PCA_Magnitude_CV_2.csv",
    loadings_var = "loadings_PCA_Variability_CV_2.csv",
    csv_sep = ",",
    csv_dec = ".",
    var_cv_cols = c("sst_cv_2", "dli_cv_2", "chl_cv_2")
  ),
  CV_30 = list(
    dir = file.path(base_dir, "#####output_local_PCA_CV_30_FINAL"),
    scores_file = "dados_consolidados_com_scores_das_duas_PCAs.xlsx",
    loadings_mag = "loadings_PCA_Magnitude_CV_30.csv",
    loadings_var = "loadings_PCA_Variability_CV_30.csv",
    csv_sep = ",",
    csv_dec = ".",
    var_cv_cols = c("sst_cv_30", "dli_cv_30", "chl_cv_30")
  ),
  CV_ALL = list(
    dir = file.path(base_dir, "#####output_local_PCA_CV_all_FINAL"),
    scores_file = "dados_consolidados_com_scores_das_duas_PCAs.xlsx",
    loadings_mag = "loadings_PCA_Magnitude_CV_all.csv",
    loadings_var = "loadings_PCA_Variability_CV_all.csv",
    csv_sep = ";",
    csv_dec = ",",
    var_cv_cols = c("sst_cv_all", "cv_DLI_local", "chl_cv_all")
  )
)

output_dir <- file.path(base_dir, "PCA_Publication_Figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 3. ESQUEMA VISUAL ---
# Cores para recifes (tab10 colormap equivalente)
reef_colors <- c(
  "ARC" = "#1f77b4",  # azul
  "ITA" = "#ff7f0e",  # laranja
  "PAB" = "#2ca02c",  # verde
  "UCR" = "#d62728",  # vermelho
  "TIM" = "#9467bd"   # roxo
)

# Formas para habitats
habitat_shapes <- c("PA" = 21, "RR" = 22, "TP" = 24)  # círculo, quadrado, triângulo (preenchíveis)

# Mapeamento de nomes de variáveis para labels técnicos
var_labels <- list(
  # Magnitude
  "sst_mean" = "SST (°C)",
  "mean_DLI_local" = "DLI (mol m⁻² d⁻¹)",
  "chl_mean" = "Chl-a (mg m⁻³)",
  "prop_DHW_gt4" = "DHW>4 (freq)",
  # Variabilidade (genérico)
  "sst_cv" = "SST CV (%)",
  "dli_cv" = "DLI CV (%)",
  "chl_cv" = "Chl-a CV (%)"
)

# Mapeamento de loadings para labels
loadings_labels <- list(
  "log(sst_mean+1)" = "log(SST+1)",
  "log(mean_DLI_local+1)" = "log(DLI+1)",
  "log(chl_mean+1)" = "log(Chl-a+1)",
  "sqrt(DHW>4)" = "√(DHW>4)",
  "sst_cv_2" = "SST CV",
  "sst_cv_30" = "SST CV",
  "sst_cv_all" = "SST CV",
  "dli_cv_2" = "DLI CV",
  "dli_cv_30" = "DLI CV",
  "cv_DLI_local" = "DLI CV",
  "chl_cv_2" = "Chl-a CV",
  "chl_cv_30" = "Chl-a CV",
  "chl_cv_all" = "Chl-a CV"
)

# Tema profissional (inspirado no script Bayesiano)
theme_publication <- theme_classic(base_size = 12) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major = element_line(color = "grey90", linetype = "dashed", linewidth = 0.3),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 11)
  )
```

---

## Funções de Plotagem

### 4.1 Função: Bubble Plot Individual

```r
create_bubble_panel <- function(df, var_col, pc1_col, pc2_col, 
                                 exp_var, var_label, show_legend = FALSE) {
  # Escalar tamanhos (min 2, max 8)
  size_vals <- df[[var_col]]
  size_scaled <- scales::rescale(size_vals, to = c(2, 8), 
                                  from = range(size_vals, na.rm = TRUE))
  df$size_scaled <- size_scaled
  
  p <- ggplot(df, aes(x = .data[[pc1_col]], y = .data[[pc2_col]])) +
    # Linha de referência em 0
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
    # Pontos com estética completa
    geom_point(aes(fill = Reef_name, 
                   shape = HAB, 
                   size = size_scaled,
                   stroke = ifelse(Arch == "inner", 1.5, 0.5)),
               color = ifelse(df$Arch == "inner", "black", "grey60"),
               alpha = 0.8) +
    # Escalas
    scale_fill_manual(values = reef_colors, name = "Reef") +
    scale_shape_manual(values = habitat_shapes, name = "Habitat") +
    scale_size_identity() +  # Usar tamanhos já calculados
    # Labels
    labs(
      title = var_label,
      x = sprintf("PC1 (%.1f%%)", exp_var[1]),
      y = sprintf("PC2 (%.1f%%)", exp_var[2])
    ) +
    theme_publication +
    coord_fixed(ratio = 1)  # Manter proporção dos eixos
  
  if (!show_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  return(p)
}
```

### 4.2 Função: Painel de Loadings (Biplot)

```r
create_loadings_panel <- function(loadings_df, exp_var, arrow_scale = 1.5) {
  # Preparar dados para setas
  loadings_df$xend <- loadings_df$PC1 * arrow_scale
  loadings_df$yend <- loadings_df$PC2 * arrow_scale
  
  # Criar labels limpos
  loadings_df$label <- sapply(loadings_df$variable, function(v) {
    if (v %in% names(loadings_labels)) loadings_labels[[v]] else v
  })
  
  p <- ggplot(loadings_df) +
    # Linhas de referência
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
    # Círculo unitário (opcional, para referência)
    annotate("path",
             x = cos(seq(0, 2*pi, length.out = 100)),
             y = sin(seq(0, 2*pi, length.out = 100)),
             color = "grey80", linetype = "dashed", linewidth = 0.3) +
    # Setas de loadings
    geom_segment(aes(x = 0, y = 0, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
                 color = "#d62728", linewidth = 0.8) +
    # Labels das variáveis
    geom_text(aes(x = xend * 1.15, y = yend * 1.15, label = label),
              size = 3.5, fontface = "bold", color = "black") +
    # Labels
    labs(
      title = "Loadings",
      x = sprintf("PC1 (%.1f%%)", exp_var[1]),
      y = sprintf("PC2 (%.1f%%)", exp_var[2])
    ) +
    theme_publication +
    coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))
  
  return(p)
}
```

### 4.3 Função: Legenda Composta

```r
create_legend_panel <- function(unique_reefs, unique_habs) {
  # Criar dados dummy para gerar legenda
  legend_data <- expand.grid(Reef_name = unique_reefs, HAB = unique_habs[1])
  
  # Plot dummy para extrair legendas
  p_reef <- ggplot(legend_data, aes(x = 1, y = 1, fill = Reef_name)) +
    geom_point(size = 5, shape = 21) +
    scale_fill_manual(values = reef_colors, name = "Reef") +
    theme_void() +
    theme(legend.position = "right")
  
  p_hab <- ggplot(data.frame(HAB = unique_habs), aes(x = 1, y = 1, shape = HAB)) +
    geom_point(size = 5, fill = "grey50") +
    scale_shape_manual(values = habitat_shapes, name = "Habitat") +
    theme_void() +
    theme(legend.position = "right")
  
  p_arch <- ggplot(data.frame(Arch = c("inner", "outer")), 
                   aes(x = 1, y = 1:2)) +
    geom_point(size = 5, shape = 21, fill = "grey70",
               color = c("black", "grey60"), stroke = c(1.5, 0.5)) +
    annotate("text", x = 1.5, y = 1:2, label = c("Inner Arc", "Outer Arc"), 
             hjust = 0, size = 3.5) +
    xlim(0.5, 3) +
    labs(title = "Arc") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 10))
  
  # Extrair grobs de legenda
  legend_reef <- cowplot::get_legend(p_reef)
  legend_hab <- cowplot::get_legend(p_hab)
  legend_arch <- cowplot::get_legend(p_arch)
  
  # Combinar legendas verticalmente
  combined_legend <- cowplot::plot_grid(
    legend_reef, legend_hab, legend_arch,
    ncol = 1, align = "v", rel_heights = c(1.2, 0.8, 0.6)
  )
  
  return(combined_legend)
}
```

---

## Funções de Composição

### 5.1 Figura de Magnitude (2×3)

```r
create_magnitude_figure <- function(df, loadings_df, exp_var, scenario_name) {
  # Variáveis de magnitude
  mag_vars <- c("sst_mean", "mean_DLI_local", "chl_mean", "prop_DHW_gt4")
  mag_labels <- c("SST (°C)", "DLI (mol m⁻² d⁻¹)", "Chl-a (mg m⁻³)", "DHW>4 (freq)")
  
  # Criar bubble plots
  bubbles <- lapply(seq_along(mag_vars), function(i) {
    create_bubble_panel(df, mag_vars[i], "PC1_Magnitude", "PC2_Magnitude",
                        exp_var, mag_labels[i], show_legend = FALSE)
  })
  
  # Loadings panel
  p_loadings <- create_loadings_panel(loadings_df, exp_var)
  
  # Legend panel
  p_legend <- create_legend_panel(unique(df$Reef_name), unique(df$HAB))
  
  # Compor figura: 2 linhas × 3 colunas
  # Linha 1: 3 bubbles (SST, DLI, CHL)
  # Linha 2: DHW, Loadings, Legenda
  
  row1 <- bubbles[[1]] | bubbles[[2]] | bubbles[[3]]
  row2 <- bubbles[[4]] | p_loadings | p_legend
  
  final_figure <- row1 / row2 +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(
      tag_levels = 'a',
      tag_suffix = ')',
      theme = theme(plot.tag = element_text(face = "bold", size = 12))
    )
  
  return(final_figure)
}
```

### 5.2 Figura de Variabilidade (2×2)

```r
create_variability_figure <- function(df, loadings_df, exp_var, 
                                       cv_cols, scenario_name) {
  # Labels para variáveis de CV
  cv_labels <- c("SST CV (%)", "DLI CV (%)", "Chl-a CV (%)")
  
  # Criar bubble plots
  bubbles <- lapply(seq_along(cv_cols), function(i) {
    create_bubble_panel(df, cv_cols[i], "PC1_Variability", "PC2_Variability",
                        exp_var, cv_labels[i], show_legend = FALSE)
  })
  
  # Loadings panel
  p_loadings <- create_loadings_panel(loadings_df, exp_var)
  
  # Compor figura: 2 linhas × 2 colunas
  # Com legenda ao lado direito da figura inteira
  
  main_grid <- (bubbles[[1]] | bubbles[[2]]) / (bubbles[[3]] | p_loadings) +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(
      tag_levels = 'a',
      tag_suffix = ')',
      theme = theme(plot.tag = element_text(face = "bold", size = 12))
    )
  
  # Adicionar legenda externa à direita
  p_legend <- create_legend_panel(unique(df$Reef_name), unique(df$HAB))
  
  final_figure <- cowplot::plot_grid(
    main_grid, p_legend,
    ncol = 2, rel_widths = c(4, 1)
  )
  
  return(final_figure)
}
```

---

## Loop Principal

```r
# --- 6. PROCESSAMENTO PARA CADA CENÁRIO ---
for (scenario_name in names(scenarios)) {
  cat(sprintf("\n=== Processando cenário: %s ===\n", scenario_name))
  
  sc <- scenarios[[scenario_name]]
  
  # 6.1 Carregar dados de scores
  scores_path <- file.path(sc$dir, sc$scores_file)
  df <- readxl::read_excel(scores_path)
  cat(sprintf("  Scores carregados: %d sítios\n", nrow(df)))
  
  # 6.2 Carregar loadings (com tratamento de separador)
  load_loadings <- function(file_path, sep, dec) {
    if (sep == ";") {
      ld <- read.csv2(file_path, row.names = 1)
    } else {
      ld <- read.csv(file_path, row.names = 1)
    }
    ld$variable <- rownames(ld)
    return(ld)
  }
  
  loadings_mag <- load_loadings(file.path(sc$dir, sc$loadings_mag), sc$csv_sep, sc$csv_dec)
  loadings_var <- load_loadings(file.path(sc$dir, sc$loadings_var), sc$csv_sep, sc$csv_dec)
  
  # 6.3 Calcular variância explicada (PCA de magnitude)
  mag_vars <- c("sst_mean", "mean_DLI_local", "chl_mean", "prop_DHW_gt4")
  df_mag_clean <- df[complete.cases(df[, mag_vars]), mag_vars]
  pca_mag <- prcomp(df_mag_clean, scale. = TRUE)
  exp_var_mag <- summary(pca_mag)$importance[2, 1:2] * 100
  
  # 6.4 Calcular variância explicada (PCA de variabilidade)
  var_cols <- sc$var_cv_cols
  df_var_clean <- df[complete.cases(df[, var_cols]), var_cols]
  pca_var <- prcomp(df_var_clean, scale. = TRUE)
  exp_var_var <- summary(pca_var)$importance[2, 1:2] * 100
  
  cat(sprintf("  Var. explicada (Magnitude): PC1=%.1f%%, PC2=%.1f%%\n", 
              exp_var_mag[1], exp_var_mag[2]))
  cat(sprintf("  Var. explicada (Variability): PC1=%.1f%%, PC2=%.1f%%\n", 
              exp_var_var[1], exp_var_var[2]))
  
  # 6.5 Gerar figura de Magnitude
  fig_mag <- create_magnitude_figure(df, loadings_mag, exp_var_mag, scenario_name)
  ggsave(
    file.path(output_dir, sprintf("Figure_PCA_Magnitude_%s.png", scenario_name)),
    fig_mag, width = 14, height = 10, dpi = 300
  )
  cat(sprintf("  ✓ Figura Magnitude salva\n"))
  
  # 6.6 Gerar figura de Variabilidade
  fig_var <- create_variability_figure(df, loadings_var, exp_var_var, var_cols, scenario_name)
  ggsave(
    file.path(output_dir, sprintf("Figure_PCA_Variability_%s.png", scenario_name)),
    fig_var, width = 12, height = 10, dpi = 300
  )
  cat(sprintf("  ✓ Figura Variability salva\n"))
}

cat("\n=== PROCESSAMENTO CONCLUÍDO ===\n")
cat(sprintf("Figuras salvas em: %s\n", output_dir))
```

---

## Saídas Esperadas

| Arquivo | Dimensões | Cenário |
|---------|-----------|---------|
| `Figure_PCA_Magnitude_CV_02.png` | 14×10 in, 300 DPI | CV 2 dias |
| `Figure_PCA_Variability_CV_02.png` | 12×10 in, 300 DPI | CV 2 dias |
| `Figure_PCA_Magnitude_CV_30.png` | 14×10 in, 300 DPI | CV 30 dias |
| `Figure_PCA_Variability_CV_30.png` | 12×10 in, 300 DPI | CV 30 dias |
| `Figure_PCA_Magnitude_CV_ALL.png` | 14×10 in, 300 DPI | CV total |
| `Figure_PCA_Variability_CV_ALL.png` | 12×10 in, 300 DPI | CV total |

---

## Próximo Passo

Após aprovação deste plano, o script completo será criado em:
```
#######FINAL_CODES/PCA_Visualization_Publication.R
```
