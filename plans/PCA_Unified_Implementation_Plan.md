# Plano Unificado: Visualização de PCA para Publicação

## Diagnóstico Comparativo

### Convergências ✅

| Aspecto | Ambos os Planos Concordam |
|---------|---------------------------|
| **Linguagem** | R com ggplot2 + patchwork + cowplot |
| **Pacotes** | readxl, ggplot2, patchwork, cowplot, dplyr, scales |
| **3 cenários** | CV_02, CV_30, CV_ALL processados automaticamente |
| **Tema** | `theme_classic()` com base_size 12-14 |
| **Símbolos HAB** | shapes 16/21 (círculo), 15/22 (quadrado), 17/24 (triângulo) |
| **Bordas Arch** | Inner: preto + grosso (1.5-2.0); Outer: cinza + fino (0.5-0.75) |
| **Layout Magnitude** | 2×3 grid (4 bubbles + loadings + legenda) |
| **Layout Variability** | 2×2 grid (3 bubbles + loadings) + legenda externa |
| **Sem título principal** | Confirmado: sem `fig.suptitle` no topo |
| **Labels em inglês** | Todas as legendas e eixos em inglês |
| **Variância explicada** | Nos rótulos dos eixos: "PC1 (XX.X%)" |
| **Loadings biplot** | Setas com círculo unitário de referência |
| **Saída PNG** | 300 DPI, ~14×10 in (Magnitude), ~12×10 in (Variability) |

### Diferenças e Resolução 🔄

| Aspecto | Plano 1 (Redo) | Plano 2 (Refactoring) | **Resolução Unificada** |
|---------|----------------|----------------------|-------------------------|
| **Cores recifes** | Paleta custom ("Abrolhos"="#D55E00", etc.) | tab10 (ARC="#1f77b4", ITA="#ff7f0e", etc.) | **Usar tab10** (consistente com script Python original) |
| **Nomes recifes** | Nomes longos (Abrolhos, Parcel das Paredes) | Códigos curtos (ARC, ITA, PAB, UCR, TIM) | **Códigos curtos** (como estão nos dados) |
| **Variância explicada** | Não especifica fonte | Recalcula via `prcomp()` | **Recalcular** (mais preciso e robusto) |
| **Separadores CSV** | Assume read.csv genérico | Detecta `,` vs `;` por cenário | **Detectar separador** (CV_all usa `;`) |
| **Labels variáveis** | "SST Mean", "DLI CV" simples | "SST (°C)", "DLI (mol m⁻² d⁻¹)" com unidades | **Com unidades** para publicação |
| **Labels loadings** | Nomes transformados longos | Abreviados: log(SST+1), √(DHW>4) | **Abreviados limpos** |
| **Saída PDF** | Sugere PNG + PDF/SVG | Apenas PNG | **PNG + PDF** (para submissão) |
| **Arquivo de script** | `PCA_Visualization_Only.R` | `PCA_Visualization_Publication.R` | **`PCA_Visualization_Publication.R`** |

---

## Especificações Finais Unificadas

### Dados de Entrada

```
#######FINAL_RESULTS/
├── #####output_local_PCA_CV_2_FINAL/
│   ├── dados_consolidados_com_scores_das_duas_PCAs.xlsx
│   ├── loadings_PCA_Magnitude_CV_2.csv      (sep=",", dec=".")
│   └── loadings_PCA_Variability_CV_2.csv    (sep=",", dec=".")
├── #####output_local_PCA_CV_30_FINAL/
│   ├── dados_consolidados_com_scores_das_duas_PCAs.xlsx
│   ├── loadings_PCA_Magnitude_CV_30.csv     (sep=",", dec=".")
│   └── loadings_PCA_Variability_CV_30.csv   (sep=",", dec=".")
└── #####output_local_PCA_CV_all_FINAL/
    ├── dados_consolidados_com_scores_das_duas_PCAs.xlsx
    ├── loadings_PCA_Magnitude_CV_all.csv    (sep=";", dec=",")  # ⚠️ Diferente!
    └── loadings_PCA_Variability_CV_all.csv  (sep=";", dec=",")
```

### Colunas do Excel (23 colunas confirmadas)

| # | Coluna | Tipo | Uso |
|---|--------|------|-----|
| 0 | `Reef_name` | Factor (5 níveis) | Cor dos pontos |
| 2 | `HAB` | Factor (3 níveis) | Forma dos pontos |
| 6 | `Arch` | Factor (inner/outer) | Estilo da borda |
| 11 | `sst_mean` | Numeric | Bubble size (Magnitude) |
| 8 | `mean_DLI_local` | Numeric | Bubble size (Magnitude) |
| 16 | `chl_mean` | Numeric | Bubble size (Magnitude) |
| 14 | `prop_DHW_gt4` | Numeric | Bubble size (Magnitude) |
| 12/9 | `sst_cv_*` | Numeric | Bubble size (Variability) |
| 9 | `dli_cv_*` / `cv_DLI_local` | Numeric | Bubble size (Variability) |
| 17 | `chl_cv_*` | Numeric | Bubble size (Variability) |
| 21 | `PC1_Magnitude` | Numeric | Eixo X (Magnitude) |
| 22 | `PC2_Magnitude` | Numeric | Eixo Y (Magnitude) |
| 19 | `PC1_Variability` | Numeric | Eixo X (Variability) |
| 20 | `PC2_Variability` | Numeric | Eixo Y (Variability) |

---

## Esquema Visual Final

### Cores (Recifes) - tab10

```r
reef_colors <- c(
  "ARC" = "#1f77b4",
  "ITA" = "#ff7f0e",
  "PAB" = "#2ca02c",
  "UCR" = "#d62728",
  "TIM" = "#9467bd"
)
```

### Formas (Habitats)

```r
habitat_shapes <- c("PA" = 21, "RR" = 22, "TP" = 24)  # preenchíveis
```

### Bordas (Arco)

```r
arch_styles <- list(
  inner = list(color = "black",   stroke = 1.5),
  outer = list(color = "grey60",  stroke = 0.5)
)
```

### Labels Técnicos

| Variável | Bubble Title | Loading Label |
|----------|--------------|---------------|
| `sst_mean` | SST (°C) | log(SST+1) |
| `mean_DLI_local` | DLI (mol m⁻² d⁻¹) | log(DLI+1) |
| `chl_mean` | Chl-a (mg m⁻³) | log(Chl+1) |
| `prop_DHW_gt4` | DHW>4 (freq) | √(DHW>4) |
| `sst_cv_*` | SST CV (%) | SST CV |
| `dli_cv_*` | DLI CV (%) | DLI CV |
| `chl_cv_*` | Chl-a CV (%) | Chl CV |

### Tema Profissional

```r
theme_publication <- theme_classic(base_size = 12) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major = element_line(color = "grey90", linetype = "dashed", linewidth = 0.3),
    strip.background = element_blank(),
    plot.tag = element_text(face = "bold", size = 12)
  )
```

---

## Layout das Figuras

### Figura 1: Magnitude (2×3 = 6 células)

```
┌──────────────┬──────────────┬──────────────┐
│  a) SST (°C) │ b) DLI (mol) │ c) Chl-a     │
│   [bubble]   │   [bubble]   │   [bubble]   │
├──────────────┼──────────────┼──────────────┤
│ d) DHW>4     │ e) Loadings  │  [LEGENDA]   │
│   [bubble]   │   [biplot]   │  Reef/HAB/   │
│              │              │  Arc         │
└──────────────┴──────────────┴──────────────┘
        ↓ PC1 (XX.X%) / PC2 (YY.Y%) ↓
```

### Figura 2: Variability (2×2 + legenda externa)

```
┌──────────────┬──────────────┐     ┌────────┐
│ a) SST CV(%) │ b) DLI CV(%) │     │ Reef   │
│   [bubble]   │   [bubble]   │     │ ────── │
├──────────────┼──────────────┤     │ HAB    │
│ c) Chl-a CV  │ d) Loadings  │     │ ────── │
│   [bubble]   │   [biplot]   │     │ Arc    │
└──────────────┴──────────────┘     └────────┘
```

---

## Código Completo (Estrutura Modular)

### 1. Configuração e Pacotes

```r
# ============================================================================
# PCA_Visualization_Publication.R
# Geração de figuras PCA de alta qualidade para publicação
# Projeto: Mussismilia hispida - Abrolhos
# ============================================================================

libs <- c("readxl", "readr", "dplyr", "ggplot2", "patchwork", "cowplot", "scales")
invisible(lapply(libs, library, character.only = TRUE))

# --- Caminhos ---
base_dir <- "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_RESULTS"

output_dir <- file.path(base_dir, "PCA_Publication_Figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
```

### 2. Configuração dos Cenários

```r
scenarios <- list(
  CV_02 = list(
    dir = file.path(base_dir, "#####output_local_PCA_CV_2_FINAL"),
    loadings_mag = "loadings_PCA_Magnitude_CV_2.csv",
    loadings_var = "loadings_PCA_Variability_CV_2.csv",
    csv_sep = ",", csv_dec = ".",
    var_cv_cols = c("sst_cv_2", "dli_cv_2", "chl_cv_2")
  ),
  CV_30 = list(
    dir = file.path(base_dir, "#####output_local_PCA_CV_30_FINAL"),
    loadings_mag = "loadings_PCA_Magnitude_CV_30.csv",
    loadings_var = "loadings_PCA_Variability_CV_30.csv",
    csv_sep = ",", csv_dec = ".",
    var_cv_cols = c("sst_cv_30", "dli_cv_30", "chl_cv_30")
  ),
  CV_ALL = list(
    dir = file.path(base_dir, "#####output_local_PCA_CV_all_FINAL"),
    loadings_mag = "loadings_PCA_Magnitude_CV_all.csv",
    loadings_var = "loadings_PCA_Variability_CV_all.csv",
    csv_sep = ";", csv_dec = ",",
    var_cv_cols = c("sst_cv_all", "cv_DLI_local", "chl_cv_all")
  )
)
```

### 3. Esquema Visual

```r
# Cores tab10 (consistente com Python original)
reef_colors <- c(
  "ARC" = "#1f77b4", "ITA" = "#ff7f0e", "PAB" = "#2ca02c",
  "UCR" = "#d62728", "TIM" = "#9467bd"
)

# Formas preenchíveis
habitat_shapes <- c("PA" = 21, "RR" = 22, "TP" = 24)

# Labels técnicos com unidades
var_labels_mag <- c(
  "sst_mean" = "SST (°C)",
  "mean_DLI_local" = "DLI (mol m⁻² d⁻¹)",
  "chl_mean" = "Chl-a (mg m⁻³)",
  "prop_DHW_gt4" = "DHW>4 (freq)"
)

var_labels_var <- c(
  "SST CV (%)", "DLI CV (%)", "Chl-a CV (%)"
)

# Labels para loadings (abreviados)
loadings_labels <- c(
  "log(sst_mean+1)" = "log(SST+1)",
  "log(mean_DLI_local+1)" = "log(DLI+1)",
  "log(chl_mean+1)" = "log(Chl+1)",
  "sqrt(DHW>4)" = "√(DHW>4)"
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
    panel.grid.major = element_line(color = "grey90", linetype = "dashed", linewidth = 0.3),
    plot.tag = element_text(face = "bold", size = 12)
  )
```

### 4. Funções de Plotagem

```r
# --- 4.1 Bubble Plot ---
create_bubble_panel <- function(df, var_col, pc1_col, pc2_col, exp_var, var_label) {
  # Escalar tamanhos
  df$size_scaled <- scales::rescale(df[[var_col]], to = c(2, 8), 
                                     from = range(df[[var_col]], na.rm = TRUE))
  # Estilo de borda baseado em Arch
  df$border_color <- ifelse(df$Arch == "inner", "black", "grey60")
  df$border_width <- ifelse(df$Arch == "inner", 1.5, 0.5)
  
  ggplot(df, aes(x = .data[[pc1_col]], y = .data[[pc2_col]])) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
    geom_point(aes(fill = Reef_name, shape = HAB, size = size_scaled),
               color = df$border_color, stroke = df$border_width, alpha = 0.85) +
    scale_fill_manual(values = reef_colors, name = "Reef") +
    scale_shape_manual(values = habitat_shapes, name = "Habitat") +
    scale_size_identity() +
    labs(title = var_label,
         x = sprintf("PC1 (%.1f%%)", exp_var[1]),
         y = sprintf("PC2 (%.1f%%)", exp_var[2])) +
    theme_publication +
    theme(legend.position = "none") +
    coord_fixed(ratio = 1)
}

# --- 4.2 Loadings Biplot ---
create_loadings_panel <- function(loadings_df, exp_var, arrow_scale = 1.5) {
  loadings_df$xend <- loadings_df$PC1 * arrow_scale
  loadings_df$yend <- loadings_df$PC2 * arrow_scale
  loadings_df$label <- sapply(loadings_df$variable, function(v) {
    if (v %in% names(loadings_labels)) loadings_labels[[v]] else gsub("_", " ", v)
  })
  
  ggplot(loadings_df) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
    annotate("path", x = cos(seq(0, 2*pi, length.out = 100)),
             y = sin(seq(0, 2*pi, length.out = 100)),
             color = "grey80", linetype = "dashed", linewidth = 0.3) +
    geom_segment(aes(x = 0, y = 0, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                 color = "#d62728", linewidth = 0.8) +
    geom_text(aes(x = xend * 1.15, y = yend * 1.15, label = label),
              size = 3.5, fontface = "bold") +
    labs(title = "Loadings",
         x = sprintf("PC1 (%.1f%%)", exp_var[1]),
         y = sprintf("PC2 (%.1f%%)", exp_var[2])) +
    theme_publication +
    coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))
}

# --- 4.3 Painel de Legenda ---
create_legend_panel <- function(reefs, habs) {
  # Legendas extraídas de plots dummy
  p_reef <- ggplot(data.frame(Reef = reefs), aes(x = 1, y = seq_along(Reef), fill = Reef)) +
    geom_point(shape = 21, size = 4) +
    scale_fill_manual(values = reef_colors, name = "Reef") +
    theme_void() + theme(legend.position = "right")
  
  p_hab <- ggplot(data.frame(HAB = habs), aes(x = 1, y = seq_along(HAB), shape = HAB)) +
    geom_point(size = 4, fill = "grey60") +
    scale_shape_manual(values = habitat_shapes, name = "Habitat") +
    theme_void() + theme(legend.position = "right")
  
  # Legenda manual para Arch
  arch_grob <- cowplot::ggdraw() +
    cowplot::draw_label("Arc", x = 0.5, y = 0.9, fontface = "bold", size = 10) +
    cowplot::draw_plot(
      ggplot(data.frame(x = 1, y = 1:2, 
                        label = c("Inner Arc", "Outer Arc"),
                        stroke = c(1.5, 0.5),
                        color = c("black", "grey60"))) +
        geom_point(aes(x, y), shape = 21, size = 4, fill = "grey70",
                   color = c("black", "grey60"), stroke = c(1.5, 0.5)) +
        geom_text(aes(x + 0.4, y, label = label), hjust = 0, size = 3) +
        xlim(0.5, 2.5) + theme_void(),
      x = 0, y = 0, width = 1, height = 0.8
    )
  
  cowplot::plot_grid(
    cowplot::get_legend(p_reef),
    cowplot::get_legend(p_hab),
    arch_grob,
    ncol = 1, rel_heights = c(1.2, 0.8, 0.8)
  )
}
```

### 5. Funções de Composição

```r
# --- 5.1 Figura Magnitude (2×3) ---
create_magnitude_figure <- function(df, loadings_df, exp_var) {
  mag_vars <- c("sst_mean", "mean_DLI_local", "chl_mean", "prop_DHW_gt4")
  
  bubbles <- lapply(seq_along(mag_vars), function(i) {
    create_bubble_panel(df, mag_vars[i], "PC1_Magnitude", "PC2_Magnitude",
                        exp_var, var_labels_mag[mag_vars[i]])
  })
  
  p_loadings <- create_loadings_panel(loadings_df, exp_var)
  p_legend <- create_legend_panel(unique(df$Reef_name), unique(df$HAB))
  
  # Composição patchwork
  (bubbles[[1]] | bubbles[[2]] | bubbles[[3]]) /
  (bubbles[[4]] | p_loadings | p_legend) +
    plot_annotation(tag_levels = 'a', tag_suffix = ')')
}

# --- 5.2 Figura Variability (2×2 + legenda externa) ---
create_variability_figure <- function(df, loadings_df, exp_var, cv_cols) {
  bubbles <- lapply(seq_along(cv_cols), function(i) {
    create_bubble_panel(df, cv_cols[i], "PC1_Variability", "PC2_Variability",
                        exp_var, var_labels_var[i])
  })
  
  p_loadings <- create_loadings_panel(loadings_df, exp_var)
  p_legend <- create_legend_panel(unique(df$Reef_name), unique(df$HAB))
  
  main_grid <- (bubbles[[1]] | bubbles[[2]]) / (bubbles[[3]] | p_loadings) +
    plot_annotation(tag_levels = 'a', tag_suffix = ')')
  
  cowplot::plot_grid(main_grid, p_legend, ncol = 2, rel_widths = c(4, 1))
}
```

### 6. Loop Principal

```r
# --- 6. PROCESSAMENTO ---
for (scenario_name in names(scenarios)) {
  cat(sprintf("\n=== Processando: %s ===\n", scenario_name))
  sc <- scenarios[[scenario_name]]
  
  # Carregar dados
  df <- readxl::read_excel(file.path(sc$dir, "dados_consolidados_com_scores_das_duas_PCAs.xlsx"))
  
  load_loadings <- function(path, sep) {
    if (sep == ";") read.csv2(path, row.names = 1) else read.csv(path, row.names = 1)
  }
  loadings_mag <- load_loadings(file.path(sc$dir, sc$loadings_mag), sc$csv_sep)
  loadings_var <- load_loadings(file.path(sc$dir, sc$loadings_var), sc$csv_sep)
  loadings_mag$variable <- rownames(loadings_mag)
  loadings_var$variable <- rownames(loadings_var)
  
  # Calcular variância explicada
  mag_vars <- c("sst_mean", "mean_DLI_local", "chl_mean", "prop_DHW_gt4")
  exp_var_mag <- summary(prcomp(na.omit(df[, mag_vars]), scale. = TRUE))$importance[2, 1:2] * 100
  exp_var_var <- summary(prcomp(na.omit(df[, sc$var_cv_cols]), scale. = TRUE))$importance[2, 1:2] * 100
  
  # Gerar e salvar figuras
  fig_mag <- create_magnitude_figure(df, loadings_mag, exp_var_mag)
  ggsave(file.path(output_dir, sprintf("Figure_PCA_Magnitude_%s.png", scenario_name)),
         fig_mag, width = 14, height = 10, dpi = 300)
  ggsave(file.path(output_dir, sprintf("Figure_PCA_Magnitude_%s.pdf", scenario_name)),
         fig_mag, width = 14, height = 10)
  
  fig_var <- create_variability_figure(df, loadings_var, exp_var_var, sc$var_cv_cols)
  ggsave(file.path(output_dir, sprintf("Figure_PCA_Variability_%s.png", scenario_name)),
         fig_var, width = 12, height = 10, dpi = 300)
  ggsave(file.path(output_dir, sprintf("Figure_PCA_Variability_%s.pdf", scenario_name)),
         fig_var, width = 12, height = 10)
  
  cat(sprintf("  ✓ %s completo\n", scenario_name))
}

cat("\n=== PROCESSAMENTO CONCLUÍDO ===\n")
cat(sprintf("Figuras em: %s\n", output_dir))
```

---

## Saídas Esperadas

| Arquivo | Formato | Dimensões |
|---------|---------|-----------|
| `Figure_PCA_Magnitude_CV_02.png` | PNG 300 DPI | 14×10 in |
| `Figure_PCA_Magnitude_CV_02.pdf` | PDF vetorial | 14×10 in |
| `Figure_PCA_Variability_CV_02.png` | PNG 300 DPI | 12×10 in |
| `Figure_PCA_Variability_CV_02.pdf` | PDF vetorial | 12×10 in |
| `Figure_PCA_Magnitude_CV_30.png` | PNG 300 DPI | 14×10 in |
| `Figure_PCA_Magnitude_CV_30.pdf` | PDF vetorial | 14×10 in |
| `Figure_PCA_Variability_CV_30.png` | PNG 300 DPI | 12×10 in |
| `Figure_PCA_Variability_CV_30.pdf` | PDF vetorial | 12×10 in |
| `Figure_PCA_Magnitude_CV_ALL.png` | PNG 300 DPI | 14×10 in |
| `Figure_PCA_Magnitude_CV_ALL.pdf` | PDF vetorial | 14×10 in |
| `Figure_PCA_Variability_CV_ALL.png` | PNG 300 DPI | 12×10 in |
| `Figure_PCA_Variability_CV_ALL.pdf` | PDF vetorial | 12×10 in |

**Total: 12 arquivos (6 PNG + 6 PDF)**

---

## Checklist de Validação

- [ ] Todas as 3 cenários (CV_02, CV_30, CV_ALL) processados
- [ ] Magnitude: 6 painéis (4 bubbles + loadings + legenda)
- [ ] Variability: 4 painéis + legenda externa
- [ ] Cores correspondem aos recifes (tab10)
- [ ] Formas correspondem aos habitats (PA/RR/TP)
- [ ] Bordas distinguem arcos (inner preto grosso, outer cinza fino)
- [ ] Tamanhos proporcionais aos valores das variáveis
- [ ] Variância explicada nos rótulos dos eixos
- [ ] Loadings com setas e círculo unitário
- [ ] Labels em inglês
- [ ] Sem título principal no topo
- [ ] Tags de painel (a, b, c...)
- [ ] PNG 300 DPI + PDF vetorial

---

## Próximo Passo

Script será criado em:
```
#######FINAL_CODES/PCA_Visualization_Publication.R
```

Após execução, figuras serão salvas em:
```
#######FINAL_RESULTS/PCA_Publication_Figures/
```

---

*Plano unificado criado: 2026-01-12*
*Origem: Fusão de `PCA_Figures_Redo_Plan.md` + `pca_visualization_refactoring_plan.md`*
