# Plano Unificado v2: Visualização de PCA para Publicação

## Diagnóstico de Lacunas (v1 → v2)

### Problemas Críticos Identificados e Correções

| # | Problema (v1) | Impacto | Solução (v2) |
|---|---------------|---------|--------------|
| 1 | **Nome da coluna PCA inconsistente**: `PC1_Magnitude` vs `PC1_Magnitude_CV_2` | Erro: coluna não encontrada | Detectar automaticamente o sufixo correto |
| 2 | **Nomes de colunas CV inconsistentes**: `dli_cv_2` vs `cv_DLI_local` | Erro: coluna não encontrada no CV_all | Mapeamento explícito por cenário |
| 3 | **Loadings sem coluna 'variable'**: `rownames()` usaria índices numéricos | Labels incorretos nos loadings | Adicionar coluna `variable` ao ler CSV |
| 4 | **`border_color/width` no `aes()`**: Não funciona corretamente | Bordas não aplicadas | Mover para fora de `aes()` |
| 5 | **`exp_var` como vetor numérico**: Usado como `exp_var[1]`, `exp_var[2]` | OK mas frágil | Validar comprimento e adicionar check |
| 6 | **`var_labels_var[i]`**: Indexação numérica em vetor nomeado | Labels errados se ordem mudar | Usar `names(cv_cols)` para mapeamento |
| 7 | **`patchwork` importado mas não usado**: Uso exclusivo de `cowplot` | Confuso no código | Usar `patchwork` para composição principal |
| 8 | **Falta `library(grid)`**: Funções `arrow()` e `unit()` | Erro ao executar | Adicionar às dependências |
| 9 | **Nomenclatura de saída**: `CV_02` vs `CV_2` inconsistente | Confusão nos nomes de arquivos | Padronizar para `CV_02`, `CV_30`, `CV_ALL` |
| 10 | **Variância explicada recalculada**: Pode diferir do Python | Inconsistência com PCA original | Recalcular É válido, mas documentar |
| 11 | **Separador CSV**: `read.csv2()` vs `read.csv()` com parâmetros | Inconsistente | Padronizar com `read_delim()` |

---

## Especificações Finais Unificadas (v2)

### Dados de Entrada - Estrutura Validada

```
#######FINAL_RESULTS/
├── #####output_local_PCA_CV_2_FINAL/
│   ├── dados_consolidados_com_scores_das_duas_PCAs.xlsx
│   ├── loadings_PCA_Magnitude_CV_2.csv      (sep=",", rownames como 1ª coluna)
│   └── loadings_PCA_Variability_CV_2.csv
├── #####output_local_PCA_CV_30_FINAL/
│   ├── dados_consolidados_com_scores_das_duas_PCAs.xlsx
│   ├── loadings_PCA_Magnitude_CV_30.csv
│   └── loadings_PCA_Variability_CV_30.csv
└── #####output_local_PCA_CV_all_FINAL/
    ├── dados_consolidados_com_scores_das_duas_PCAs.xlsx
    ├── loadings_PCA_Magnitude_CV_all.csv    (sep=";") ⚠️
    └── loadings_PCA_Variability_CV_all.csv
```

### Colunas do Excel - Mapeamento Corrigido

| Categoria | Colunas CV_2/CV_30 | Colunas CV_ALL | Nota |
|-----------|-------------------|----------------|------|
| **PCA Magnitude** | `PC1_Magnitude_CV_2` / `PC2_Magnitude_CV_2` | `PC1_Magnitude` / `PC2_Magnitude` | ⚠️ **Nomes diferentes!** |
| **PCA Variability** | `PC1_Variability_CV_2` / `PC2_Variability_CV_2` | `PC1_Variability` / `PC2_Variability` | ⚠️ **Nomes diferentes!** |
| **Vars Magnitude** | `sst_mean`, `mean_DLI_local`, `chl_mean`, `prop_DHW_gt4` | (mesmo) | ✓ |
| **Vars Variability** | `sst_cv_2`, `dli_cv_2`, `chl_cv_2` | `sst_cv_all`, `cv_DLI_local`, `chl_cv_all` | ⚠️ **dli_cv vs cv_DLI_local!** |

### Colunas de Metadados

| # | Coluna | Tipo | Uso |
|---|--------|------|-----|
| 0 | `Reef_name` | Factor (ARC, ITA, PAB, UCR, TIM) | Cor dos pontos |
| 2 | `HAB` | Factor (PA, RR, TP) | Forma dos pontos |
| 6 | `Arch` | Factor (inner/outer) | Estilo da borda |

---

## Configuração dos Cenários (CORRIGIDA)

```r
scenarios <- list(
  CV_02 = list(
    dir = file.path(base_dir, "#####output_local_PCA_CV_2_FINAL"),
    loadings_mag = "loadings_PCA_Magnitude_CV_2.csv",
    loadings_var = "loadings_PCA_Variability_CV_2.csv",
    csv_sep = ",", csv_dec = ".",
    # Colunas PCA scores (COM sufixo)
    pc1_mag = "PC1_Magnitude_CV_2",
    pc2_mag = "PC2_Magnitude_CV_2",
    pc1_var = "PC1_Variability_CV_2",
    pc2_var = "PC2_Variability_CV_2",
    # Colunas para bubbles (SEM sufixo nos nomes base)
    var_mag_cols = c("sst_mean", "mean_DLI_local", "chl_mean", "prop_DHW_gt4"),
    var_var_cols = c("sst_cv_2", "dli_cv_2", "chl_cv_2"),
    # Labels para bubbles (ordem corresponde a var_var_cols)
    var_var_labels = c("SST CV (%)", "DLI CV (%)", "Chl-a CV (%)")
  ),

  CV_30 = list(
    dir = file.path(base_dir, "#####output_local_PCA_CV_30_FINAL"),
    loadings_mag = "loadings_PCA_Magnitude_CV_30.csv",
    loadings_var = "loadings_PCA_Variability_CV_30.csv",
    csv_sep = ",", csv_dec = ".",
    pc1_mag = "PC1_Magnitude_CV_30",
    pc2_mag = "PC2_Magnitude_CV_30",
    pc1_var = "PC1_Variability_CV_30",
    pc2_var = "PC2_Variability_CV_30",
    var_mag_cols = c("sst_mean", "mean_DLI_local", "chl_mean", "prop_DHW_gt4"),
    var_var_cols = c("sst_cv_30", "dli_cv_30", "chl_cv_30"),
    var_var_labels = c("SST CV (%)", "DLI CV (%)", "Chl-a CV (%)")
  ),

  CV_ALL = list(
    dir = file.path(base_dir, "#####output_local_PCA_CV_all_FINAL"),
    loadings_mag = "loadings_PCA_Magnitude_CV_all.csv",
    loadings_var = "loadings_PCA_Variability_CV_all.csv",
    csv_sep = ";", csv_dec = ",",
    # ⚠️ CV_ALL tem nomes DIFERENTES para PCA scores (sem sufixo)
    pc1_mag = "PC1_Magnitude",
    pc2_mag = "PC2_Magnitude",
    pc1_var = "PC1_Variability",
    pc2_var = "PC2_Variability",
    var_mag_cols = c("sst_mean", "mean_DLI_local", "chl_mean", "prop_DHW_gt4"),
    # ⚠️ CV_ALL usa cv_DLI_local, não dli_cv_all
    var_var_cols = c("sst_cv_all", "cv_DLI_local", "chl_cv_all"),
    var_var_labels = c("SST CV (%)", "DLI CV (%)", "Chl-a CV (%)")
  )
)
```

---

## Código Completo (Estrutura Modular v2)

### 1. Configuração e Pacotes

```r
# ============================================================================
# PCA_Visualization_Publication.R
# Geração de figuras PCA de alta qualidade para publicação
# Projeto: Mussismilia hispida - Abrolhos
# Versão: 2.0 (Corrigida)
# ============================================================================

libs <- c("readxl", "readr", "dplyr", "ggplot2", "patchwork",
          "cowplot", "scales", "grid")
invisible(lapply(libs, library, character.only = TRUE))

# --- Caminhos ---
base_dir <- "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_RESULTS"

output_dir <- file.path(base_dir, "PCA_Publication_Figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Log file
log_file <- file.path(output_dir, "processing_log.txt")
sink(log_file); on.exit(sink(), add = TRUE)
```

### 2. Esquema Visual

```r
# Cores tab10 (consistente com Python original)
reef_colors <- c(
  "ARC" = "#1f77b4", "ITA" = "#ff7f0e", "PAB" = "#2ca02c",
  "UCR" = "#d62728", "TIM" = "#9467bd"
)

# Formas preenchíveis (21=círculo, 22=quadrado, 24=triângulo)
habitat_shapes <- c("PA" = 21, "RR" = 22, "TP" = 24)

# Labels técnicos com unidades para Magnitude
var_labels_mag <- c(
  "sst_mean" = "SST (°C)",
  "mean_DLI_local" = "DLI (mol m⁻² d⁻¹)",
  "chl_mean" = "Chl-a (mg m⁻³)",
  "prop_DHW_gt4" = "DHW>4 (freq)"
)

# Labels para loadings (abreviados e limpos)
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
    panel.grid.major = element_line(color = "grey90",
                                     linetype = "dashed",
                                     linewidth = 0.3),
    strip.background = element_blank(),
    plot.tag = element_text(face = "bold", size = 12)
  )
```

### 3. Funções de Carregamento de Dados

```r
# --- Carregar loadings com tratamento de separador e rownames ---
load_loadings <- function(file_path, sep) {
  # Ler com readr para melhor controle
  df <- readr::read_delim(
    file_path,
    delim = sep,
    col_names = TRUE,
    show_col_types = FALSE
  )

  # A primeira coluna contém os nomes das variáveis
  df$variable <- df[[1]]

  # Renomear PC1 e PC2 (remover espaços ou caracteres especiais)
  colnames(df) <- gsub("^\\s+|\\s+$", "", colnames(df))

  # Remover a coluna original de nomes (agora está em 'variable')
  df <- df[, c("variable", setdiff(colnames(df), "variable"))]

  return(df)
}

# --- Carregar dados do Excel e validar colunas ---
load_scenario_data <- function(scenario) {
  cat("  - Carregando Excel...\n")

  # Carregar dados
  df <- readxl::read_excel(
    file.path(scenario$dir, "dados_consolidados_com_scores_das_duas_PCAs.xlsx")
  )

  # Validar colunas PCA
  required_cols <- c(scenario$pc1_mag, scenario$pc2_mag,
                     scenario$pc1_var, scenario$pc2_var)
  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop(sprintf("Colunas PCA não encontradas: %s\nColunas disponíveis: %s",
                 paste(missing_cols, collapse = ", "),
                 paste(names(df), collapse = ", ")))
  }

  # Validar colunas de variáveis
  all_var_cols <- c(scenario$var_mag_cols, scenario$var_var_cols)
  missing_vars <- setdiff(all_var_cols, names(df))

  if (length(missing_vars) > 0) {
    warning(sprintf("Variáveis não encontradas: %s",
                    paste(missing_vars, collapse = ", ")))
  }

  return(df)
}
```

### 4. Funções de Plotagem (CORRIGIDAS)

```r
# --- 4.1 Bubble Plot (CORRIGIDO - border fora de aes) ---
create_bubble_panel <- function(df, var_col, pc1_col, pc2_col,
                                exp_var, var_label) {
  # Escalar tamanhos
  var_values <- df[[var_col]]
  size_scaled <- scales::rescale(var_values, to = c(2, 8),
                                 from = range(var_values, na.rm = TRUE))
  df$size_scaled <- size_scaled

  # Estilo de borda baseado em Arch (fora de aes!)
  border_color <- ifelse(df$Arch == "inner", "black", "grey60")
  border_width <- ifelse(df$Arch == "inner", 1.5, 0.5)

  ggplot(df, aes(x = .data[[pc1_col]], y = .data[[pc2_col]])) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
    geom_point(aes(fill = Reef_name, shape = HAB, size = size_scaled),
               color = border_color,        # Fora de aes
               stroke = border_width,       # Fora de aes
               alpha = 0.85) +
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

# --- 4.2 Loadings Biplot (CORRIGIDO - com validação) ---
create_loadings_panel <- function(loadings_df, exp_var,
                                   arrow_scale = 1.5) {
  # Validar estrutura
  required_cols <- c("variable", "PC1", "PC2")
  missing <- setdiff(required_cols, names(loadings_df))
  if (length(missing) > 0) {
    stop(sprintf("Loadings DF missing columns: %s",
                 paste(missing, collapse = ", ")))
  }

  # Calcular pontas das setas
  loadings_df$xend <- loadings_df$PC1 * arrow_scale
  loadings_df$yend <- loadings_df$PC2 * arrow_scale

  # Aplicar labels abreviados
  loadings_df$label <- sapply(loadings_df$variable, function(v) {
    if (v %in% names(loadings_labels)) {
      loadings_labels[[v]]
    } else {
      gsub("_", " ", v)
    }
  })

  ggplot(loadings_df) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
    # Círculo unitário
    annotate("path",
             x = cos(seq(0, 2*pi, length.out = 100)),
             y = sin(seq(0, 2*pi, length.out = 100)),
             color = "grey80", linetype = "dashed", linewidth = 0.3) +
    # Setas
    geom_segment(aes(x = 0, y = 0, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
                 color = "#d62728", linewidth = 0.8) +
    # Labels
    geom_text(aes(x = xend * 1.15, y = yend * 1.15, label = label),
              size = 3.5, fontface = "bold") +
    labs(
      title = "Loadings",
      x = sprintf("PC1 (%.1f%%)", exp_var[1]),
      y = sprintf("PC2 (%.1f%%)", exp_var[2])
    ) +
    theme_publication +
    coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))
}

# --- 4.3 Painel de Legenda ---
create_legend_panel <- function(reefs, habs) {
  # Legenda de recifes (cores)
  p_reef <- ggplot(data.frame(Reef = reefs),
                   aes(x = 1, y = seq_along(Reef), fill = Reef)) +
    geom_point(shape = 21, size = 4) +
    scale_fill_manual(values = reef_colors, name = "Reef") +
    theme_void() + theme(legend.position = "right")

  # Legenda de habitats (formas)
  p_hab <- ggplot(data.frame(HAB = habs),
                  aes(x = 1, y = seq_along(HAB), shape = HAB)) +
    geom_point(size = 4, fill = "grey60") +
    scale_shape_manual(values = habitat_shapes, name = "Habitat") +
    theme_void() + theme(legend.position = "right")

  # Legenda manual para Arch
  arch_legend_data <- data.frame(
    x = 1, y = 2:1,
    label = c("Inner Arc", "Outer Arc"),
    stroke = c(1.5, 0.5),
    color = c("black", "grey60")
  )

  p_arch <- ggplot(arch_legend_data, aes(x = x, y = y)) +
    geom_point(shape = 21, size = 4, fill = "grey70",
               color = arch_legend_data$color,
               stroke = arch_legend_data$stroke) +
    geom_text(aes(x = x + 0.4, y = y, label = label), hjust = 0, size = 3) +
    xlim(0.5, 2.5) + ylim(0.5, 2.5) +
    theme_void() +
    labs(title = "Arc")

  # Combinar usando cowplot
  cowplot::plot_grid(
    cowplot::get_legend(p_reef),
    cowplot::get_legend(p_hab),
    cowplot::get_legend(p_arch),
    ncol = 1, rel_heights = c(1, 1, 1.2)
  )
}
```

### 5. Funções de Composição (CORRIGIDAS)

```r
# --- 5.1 Figura Magnitude (2×3) ---
create_magnitude_figure <- function(df, loadings_df, exp_var, var_cols) {
  cat("    - Criando painéis de bubble...\n")

  bubbles <- list()
  for (i in seq_along(var_cols)) {
    bubbles[[i]] <- create_bubble_panel(
      df, var_cols[i], "pc1_mag", "pc2_mag",
      exp_var, var_labels_mag[var_cols[i]]
    )
  }

  cat("    - Criando painel de loadings...\n")
  p_loadings <- create_loadings_panel(loadings_df, exp_var)

  cat("    - Criando painel de legenda...\n")
  p_legend <- create_legend_panel(unique(df$Reef_name), unique(df$HAB))

  cat("    - Compondo figura final...\n")
  # Composição patchwork
  (bubbles[[1]] | bubbles[[2]] | bubbles[[3]]) /
  (bubbles[[4]] | p_loadings | p_legend) +
    plot_annotation(tag_levels = 'a', tag_suffix = ')')
}

# --- 5.2 Figura Variability (2×2 + legenda externa) ---
create_variability_figure <- function(df, loadings_df, exp_var,
                                      var_cols, var_labels) {
  cat("    - Criando painéis de bubble...\n")

  bubbles <- list()
  for (i in seq_along(var_cols)) {
    bubbles[[i]] <- create_bubble_panel(
      df, var_cols[i], "pc1_var", "pc2_var",
      exp_var, var_labels[i]
    )
  }

  cat("    - Criando painel de loadings...\n")
  p_loadings <- create_loadings_panel(loadings_df, exp_var)

  cat("    - Criando painel de legenda...\n")
  p_legend <- create_legend_panel(unique(df$Reef_name), unique(df$HAB))

  cat("    - Compondo figura final...\n")
  # Grid principal 2×2
  main_grid <- (bubbles[[1]] | bubbles[[2]]) /
               (bubbles[[3]] | p_loadings) +
    plot_annotation(tag_levels = 'a', tag_suffix = ')')

  # Adicionar legenda à direita
  cowplot::plot_grid(main_grid, p_legend, ncol = 2,
                     rel_widths = c(4, 1))
}
```

### 6. Loop Principal (CORRIGIDO)

```r
# ============================================================================
# --- 6. PROCESSAMENTO PRINCIPAL ---
# ============================================================================

cat("\n===============================================\n")
cat("PCA Visualization - Publication Figures\n")
cat("===============================================\n\n")

for (scenario_name in names(scenarios)) {
  cat(sprintf("\n=== Processando: %s ===\n", scenario_name))
  sc <- scenarios[[scenario_name]]

  # Carregar dados
  df <- load_scenario_data(sc)

  # Carregar loadings
  cat("  - Carregando loadings...\n")
  loadings_mag <- load_loadings(
    file.path(sc$dir, sc$loadings_mag),
    sc$csv_sep
  )
  loadings_var <- load_loadings(
    file.path(sc$dir, sc$loadings_var),
    sc$csv_sep
  )

  # Calcular variância explicada (recalculando para consistência)
  cat("  - Calculando variância explicada...\n")

  # Magnitude
  mag_data <- na.omit(df[, sc$var_mag_cols])
  pca_mag <- prcomp(mag_data, scale. = TRUE)
  exp_var_mag <- summary(pca_mag)$importance[2, 1:2] * 100

  # Variability
  var_data <- na.omit(df[, sc$var_var_cols])
  pca_var <- prcomp(var_data, scale. = TRUE)
  exp_var_var <- summary(pca_var)$importance[2, 1:2] * 100

  cat(sprintf("    Magnitude: PC1=%.1f%%, PC2=%.1f%%\n",
              exp_var_mag[1], exp_var_mag[2]))
  cat(sprintf("    Variability: PC1=%.1f%%, PC2=%.1f%%\n",
              exp_var_var[1], exp_var_var[2]))

  # Criar figura Magnitude
  cat("  - Gerando figura Magnitude...\n")
  fig_mag <- create_magnitude_figure(df, loadings_mag, exp_var_mag, sc$var_mag_cols)

  # Salvar Magnitude
  mag_basename <- sprintf("Figure_PCA_Magnitude_%s", scenario_name)
  ggsave(
    filename = file.path(output_dir, paste0(mag_basename, ".png")),
    plot = fig_mag,
    width = 14, height = 10, dpi = 300
  )
  ggsave(
    filename = file.path(output_dir, paste0(mag_basename, ".pdf")),
    plot = fig_mag,
    width = 14, height = 10
  )
  cat(sprintf("    ✓ Salvo: %s\n", mag_basename))

  # Criar figura Variability
  cat("  - Gerando figura Variability...\n")
  fig_var <- create_variability_figure(df, loadings_var, exp_var_var,
                                       sc$var_var_cols, sc$var_var_labels)

  # Salvar Variability
  var_basename <- sprintf("Figure_PCA_Variability_%s", scenario_name)
  ggsave(
    filename = file.path(output_dir, paste0(var_basename, ".png")),
    plot = fig_var,
    width = 12, height = 10, dpi = 300
  )
  ggsave(
    filename = file.path(output_dir, paste0(var_basename, ".pdf")),
    plot = fig_var,
    width = 12, height = 10
  )
  cat(sprintf("    ✓ Salvo: %s\n", var_basename))

  cat(sprintf("  ✓ %s completo\n", scenario_name))
}

cat("\n===============================================\n")
cat("PROCESSAMENTO CONCLUÍDO\n")
cat("===============================================\n")
cat(sprintf("\nFiguras salvas em: %s\n", output_dir))
cat(sprintf("Log salvo em: %s\n", log_file))
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

**Total: 12 arquivos (6 PNG + 6 PDF) + 1 log**

---

## Checklist de Validação (v2)

### Estrutura de Dados
- [ ] Colunas PCA com sufixo corretas (CV_2, CV_30) vs sem sufixo (CV_ALL)
- [ ] Colunas CV: `dli_cv_2`/`dli_cv_30` vs `cv_DLI_local` (CV_ALL)
- [ ] Loadings CSV com separador correto por cenário
- [ ] Coluna 'variable' criada ao ler loadings

### Visual
- [ ] Magnitude: 6 painéis (4 bubbles + loadings + legenda interna)
- [ ] Variability: 4 painéis + legenda externa
- [ ] Cores correspondem aos recifes (tab10: ARC, ITA, PAB, UCR, TIM)
- [ ] Formas correspondem aos habitats (PA=21, RR=22, TP=24)
- [ ] Bordas distinguem arcos (inner preto 1.5pt, outer cinza 0.5pt)
- [ ] Tamanhos proporcionais aos valores das variáveis
- [ ] Variância explicada nos rótulos dos eixos
- [ ] Loadings com setas e círculo unitário

### Estilo
- [ ] Labels em inglês
- [ ] Sem título principal no topo
- [ ] Tags de painel (a, b, c...)
- [ ] Tema `theme_publication` aplicado
- [ ] PNG 300 DPI + PDF vetorial

### Código
- [ ] `library(grid)` carregado para arrow/unit
- [ ] Bordas fora de aes()
- [ ] Labels de variabilidade mapeados corretamente
- [ ] Nomenclatura de saída consistente (CV_02, CV_30, CV_ALL)
- [ ] Log de processamento gerado

---

## Resumo das Mudanças v1 → v2

| Categoria | Mudança | Justificativa |
|-----------|---------|---------------|
| **Nomes de colunas** | `pc1_mag`, `pc2_mag`, etc. definidos por cenário | CV_all tem nomes diferentes |
| **Variáveis CV** | Mapeamento explícito `var_var_cols` + `var_var_labels` | `dli_cv_*` vs `cv_DLI_local` |
| **Loadings** | Função `load_loadings()` extrai coluna de nomes | CSV tem nomes como 1ª coluna |
| **Bordas** | `border_color/width` movidos fora de `aes()` | ggplot requer isso |
| **Dependências** | Adicionado `library(grid)` | Funções `arrow()` e `unit()` |
| **Validação** | Função `load_scenario_data()` com checks | Prevenir erros de coluna |
| **Composição** | `patchwork` para grid principal | Código mais limpo |
| **Nomenclatura** | Padronizado para `CV_02`, `CV_30`, `CV_ALL` | Consistência |
| **Logging** | Adicionado `sink()` para log file | Rastreabilidade |

---

## Próximo Passo

Script será criado em:
```
#######FINAL_CODES/PCA_Visualization_Publication.R
```

Após execução, figuras serão salvas em:
```
#######FINAL_RESULTS/PCA_Publication_Figures/
├── Figure_PCA_Magnitude_CV_02.png
├── Figure_PCA_Magnitude_CV_02.pdf
├── Figure_PCA_Variability_CV_02.png
├── Figure_PCA_Variability_CV_02.pdf
├── ... (CV_30 e CV_ALL)
└── processing_log.txt
```

---

*Plano unificado v2 criado: 2026-01-12*
*Origem: Revisão de `PCA_Unified_Implementation_Plan.md`*
*Correções: 11 problemas críticos identificados e resolvidos*
