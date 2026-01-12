# Plano Definitivo: Visualização PCA para Publicação

> **Versão**: FINAL (v3)
> **Data**: 2026-01-12
> **Objetivo**: Orientação técnica detalhada para implementação por LLM

---

## 1. CONTEXTO E OBJETIVO

### 1.1 Problema
O script `PCA_LOCAL_bubbleplot.py` calcula métricas ambientais de sensoriamento remoto E gera figuras PCA. É necessário criar um script **independente** apenas para visualização, usando os resultados já calculados.

### 1.2 Entradas
- **3 cenários de CV** (janela deslizante): CV_02, CV_30, CV_ALL
- **Dados consolidados** em Excel (scores PCA + metadados + valores das variáveis)
- **Loadings** em CSV separados para Magnitude e Variabilidade

### 1.3 Saídas
- **12 figuras** (2 por cenário × 2 formatos)
- Formato PNG (300 DPI) + PDF vetorial
- Dimensões: Magnitude 14×10 in, Variability 12×10 in

---

## 2. ESTRUTURA DE ARQUIVOS DE ENTRADA

### 2.1 Hierarquia de Diretórios

```
C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_RESULTS/
├── #####output_local_PCA_CV_2_FINAL/
│   ├── dados_consolidados_com_scores_das_duas_PCAs.xlsx
│   ├── loadings_PCA_Magnitude_CV_2.csv
│   └── loadings_PCA_Variability_CV_2.csv
├── #####output_local_PCA_CV_30_FINAL/
│   ├── dados_consolidados_com_scores_das_duas_PCAs.xlsx
│   ├── loadings_PCA_Magnitude_CV_30.csv
│   └── loadings_PCA_Variability_CV_30.csv
└── #####output_local_PCA_CV_all_FINAL/
    ├── dados_consolidados_com_scores_das_duas_PCAs.xlsx
    ├── loadings_PCA_Magnitude_CV_all.csv   ⚠️ sep=";" dec=","
    └── loadings_PCA_Variability_CV_all.csv ⚠️ sep=";" dec=","
```

### 2.2 Estrutura do Excel (23 colunas)

| # | Nome da Coluna | Tipo | Descrição |
|---|----------------|------|-----------|
| 0 | `Reef_name` | Factor | ARC, ITA, PAB, UCR, TIM → Cor |
| 1 | `Site_name` | String | Identificador do sítio |
| 2 | `HAB` | Factor | PA, RR, TP → Forma |
| 3 | `Latitude` | Numeric | Coordenada |
| 4 | `Longitude` | Numeric | Coordenada |
| 5 | `Depth_m` | Numeric | Profundidade |
| 6 | `Arch` | Factor | inner, outer → Estilo borda |
| 7 | `unique_id` | String | Site_HAB combinado |
| 8 | `mean_DLI_local` | Numeric | DLI médio (bubble size) |
| 9 | `dli_cv_*` ou `cv_DLI_local`¹ | Numeric | CV do DLI |
| 10 | `dli_cv_*_segment_count` | Numeric | Contagem de segmentos |
| 11 | `sst_mean` | Numeric | SST média (bubble size) |
| 12 | `sst_cv_*` | Numeric | CV da SST (bubble size) |
| 13 | `sst_cv_*_segment_count` | Numeric | Contagem de segmentos |
| 14 | `prop_DHW_gt4` | Numeric | Frequência DHW>4 (bubble size) |
| 15 | `prop_DHW_gt8` | Numeric | Frequência DHW>8 |
| 16 | `chl_mean` | Numeric | Clorofila média (bubble size) |
| 17 | `chl_cv_*` | Numeric | CV da clorofila (bubble size) |
| 18 | `chl_cv_*_segment_count` | Numeric | Contagem de segmentos |
| 19 | `PC1_Variability`² | Numeric | Score PC1 Variabilidade |
| 20 | `PC2_Variability`² | Numeric | Score PC2 Variabilidade |
| 21 | `PC1_Magnitude`² | Numeric | Score PC1 Magnitude |
| 22 | `PC2_Magnitude`² | Numeric | Score PC2 Magnitude |

**Notas:**
1. CV_ALL usa `cv_DLI_local`, enquanto CV_02/CV_30 usam `dli_cv_2`/`dli_cv_30`
2. CV_02/CV_30 têm sufixo (ex: `PC1_Magnitude_CV_2`), CV_ALL não tem sufixo

### 2.3 Estrutura dos CSVs de Loadings

#### Formato CV_02 e CV_30 (sep=",", dec=".")
```csv
,PC1,PC2
log(sst_mean+1),0.5236547297496962,0.17194345477470258
log(mean_DLI_local+1),-0.18678339408208877,-0.898585541754425
log(chl_mean+1),-0.6462970241768868,0.07879889803836992
sqrt(DHW>4),-0.5226833116199978,0.395942175295818
```

#### Formato CV_ALL (sep=";", dec=",")
```csv
;PC1;PC2
log(sst_mean+1);-0,5236538518;-0,1719064287
log(mean_DLI_local+1);0,1867340419;0,8986426059
sqrt(DHW>4);0,5226988547;-0,3958288268
log(chl_mean+1);0,6462994263;-0,0787983913
```

---

## 3. MAPEAMENTO DE NOMES POR CENÁRIO

### 3.1 Configuração Explícita (CRÍTICO)

```r
scenarios <- list(
  CV_02 = list(
    # Diretório
    dir = file.path(base_dir, "#####output_local_PCA_CV_2_FINAL"),
    
    # Arquivos de loadings
    loadings_mag = "loadings_PCA_Magnitude_CV_2.csv",
    loadings_var = "loadings_PCA_Variability_CV_2.csv",
    
    # Separadores CSV
    csv_sep = ",",
    csv_dec = ".",
    
    # Colunas de scores PCA (⚠️ COM sufixo _CV_2)
    pc1_mag = "PC1_Magnitude_CV_2",
    pc2_mag = "PC2_Magnitude_CV_2",
    pc1_var = "PC1_Variability_CV_2",
    pc2_var = "PC2_Variability_CV_2",
    
    # Variáveis para bubble size (Magnitude) - sempre iguais
    var_mag_cols = c("sst_mean", "mean_DLI_local", "chl_mean", "prop_DHW_gt4"),
    var_mag_labels = c("SST (°C)", "DLI (mol m⁻² d⁻¹)", "Chl-a (mg m⁻³)", "DHW>4 (freq)"),
    
    # Variáveis para bubble size (Variability) - específicas por cenário
    var_var_cols = c("sst_cv_2", "dli_cv_2", "chl_cv_2"),
    var_var_labels = c("SST CV (%)", "DLI CV (%)", "Chl-a CV (%)")
  ),
  
  CV_30 = list(
    dir = file.path(base_dir, "#####output_local_PCA_CV_30_FINAL"),
    loadings_mag = "loadings_PCA_Magnitude_CV_30.csv",
    loadings_var = "loadings_PCA_Variability_CV_30.csv",
    csv_sep = ",",
    csv_dec = ".",
    pc1_mag = "PC1_Magnitude_CV_30",
    pc2_mag = "PC2_Magnitude_CV_30",
    pc1_var = "PC1_Variability_CV_30",
    pc2_var = "PC2_Variability_CV_30",
    var_mag_cols = c("sst_mean", "mean_DLI_local", "chl_mean", "prop_DHW_gt4"),
    var_mag_labels = c("SST (°C)", "DLI (mol m⁻² d⁻¹)", "Chl-a (mg m⁻³)", "DHW>4 (freq)"),
    var_var_cols = c("sst_cv_30", "dli_cv_30", "chl_cv_30"),
    var_var_labels = c("SST CV (%)", "DLI CV (%)", "Chl-a CV (%)")
  ),
  
  CV_ALL = list(
    dir = file.path(base_dir, "#####output_local_PCA_CV_all_FINAL"),
    loadings_mag = "loadings_PCA_Magnitude_CV_all.csv",
    loadings_var = "loadings_PCA_Variability_CV_all.csv",
    csv_sep = ";",  # ⚠️ DIFERENTE!
    csv_dec = ",",  # ⚠️ DIFERENTE!
    # ⚠️ SEM sufixo nos nomes de colunas PCA
    pc1_mag = "PC1_Magnitude",
    pc2_mag = "PC2_Magnitude",
    pc1_var = "PC1_Variability",
    pc2_var = "PC2_Variability",
    var_mag_cols = c("sst_mean", "mean_DLI_local", "chl_mean", "prop_DHW_gt4"),
    var_mag_labels = c("SST (°C)", "DLI (mol m⁻² d⁻¹)", "Chl-a (mg m⁻³)", "DHW>4 (freq)"),
    # ⚠️ cv_DLI_local (não dli_cv_all!)
    var_var_cols = c("sst_cv_all", "cv_DLI_local", "chl_cv_all"),
    var_var_labels = c("SST CV (%)", "DLI CV (%)", "Chl-a CV (%)")
  )
)
```

---

## 4. ESPECIFICAÇÕES VISUAIS

### 4.1 Cores (Recifes) - Paleta tab10

```r
reef_colors <- c(
  "ARC" = "#1f77b4",  # Azul
  "ITA" = "#ff7f0e",  # Laranja
  "PAB" = "#2ca02c",  # Verde
  "UCR" = "#d62728",  # Vermelho
  "TIM" = "#9467bd"   # Roxo
)
```

### 4.2 Formas (Habitats) - Shapes Preenchíveis

```r
habitat_shapes <- c(
  "PA" = 21,  # Círculo preenchível
  "RR" = 22,  # Quadrado preenchível
  "TP" = 24   # Triângulo preenchível
)
```

### 4.3 Bordas (Arco)

| Arco | Cor | Espessura (stroke) |
|------|-----|-------------------|
| inner | `"black"` | 1.5 |
| outer | `"grey60"` | 0.5 |

**IMPORTANTE**: Bordas devem ser definidas **FORA** de `aes()`:
```r
geom_point(aes(fill = ..., shape = ..., size = ...),
           color = border_color_vector,  # FORA de aes()
           stroke = border_width_vector) # FORA de aes()
```

### 4.4 Tamanho dos Bubbles

- Escala: valor mínimo → 2pt, valor máximo → 8pt
- Função: `scales::rescale(values, to = c(2, 8), from = range(values, na.rm = TRUE))`
- Usar `scale_size_identity()` para aplicar tamanhos calculados

### 4.5 Labels para Loadings

```r
loadings_labels <- c(
  "log(sst_mean+1)" = "log(SST+1)",
  "log(mean_DLI_local+1)" = "log(DLI+1)",
  "log(chl_mean+1)" = "log(Chl+1)",
  "sqrt(DHW>4)" = "√(DHW>4)",
  # Variability (para qualquer sufixo)
  "sst_cv" = "SST CV",
  "dli_cv" = "DLI CV",
  "chl_cv" = "Chl CV",
  "cv_DLI_local" = "DLI CV"
)
```

### 4.6 Tema Profissional

```r
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
```

---

## 5. LAYOUT DAS FIGURAS

### 5.1 Figura de Magnitude (2×3 = 6 células)

```
┌─────────────────┬─────────────────┬─────────────────┐
│  a) SST (°C)    │  b) DLI         │  c) Chl-a       │
│   [bubble]      │   [bubble]      │   [bubble]      │
│                 │   (mol m⁻²d⁻¹)  │   (mg m⁻³)      │
├─────────────────┼─────────────────┼─────────────────┤
│  d) DHW>4       │  e) Loadings    │  [LEGENDA]      │
│   (freq)        │   [biplot]      │   Reef          │
│   [bubble]      │                 │   Habitat       │
│                 │                 │   Arc           │
└─────────────────┴─────────────────┴─────────────────┘
```

**Dimensões**: 14×10 in

### 5.2 Figura de Variabilidade (2×2 + legenda externa)

```
┌─────────────────┬─────────────────┐   ┌───────────┐
│  a) SST CV (%)  │  b) DLI CV (%)  │   │   Reef    │
│   [bubble]      │   [bubble]      │   │ ────────  │
│                 │                 │   │  Habitat  │
├─────────────────┼─────────────────┤   │ ────────  │
│  c) Chl-a CV    │  d) Loadings    │   │   Arc     │
│   (%)           │   [biplot]      │   │           │
│   [bubble]      │                 │   │           │
└─────────────────┴─────────────────┘   └───────────┘
         MAIN GRID (4:1)                  LEGEND
```

**Dimensões**: 12×10 in

---

## 6. FUNÇÕES DO SCRIPT

### 6.1 Dependências

```r
libs <- c(
  "readxl",     # Ler Excel
  "readr",      # Ler CSV com controle de separador
  "dplyr",      # Manipulação de dados
  "ggplot2",    # Plotagem base
  "patchwork",  # Composição de painéis
  "cowplot",    # Extração de legendas e plot_grid
  "scales",     # Rescale de tamanhos
  "grid"        # arrow() e unit() para loadings
)
```

### 6.2 Função: Carregar Loadings

```r
load_loadings <- function(file_path, sep, dec) {
  # Ler arquivo com separador apropriado
  if (sep == ";") {
    df <- read.csv2(file_path, stringsAsFactors = FALSE)
  } else {
    df <- read.csv(file_path, stringsAsFactors = FALSE)
  }
  
  # A primeira coluna contém os nomes das variáveis (sem header)
  # Renomear para 'variable'
  colnames(df)[1] <- "variable"
  
  # Garantir que PC1 e PC2 são numéricos
  df$PC1 <- as.numeric(gsub(",", ".", df$PC1))
  df$PC2 <- as.numeric(gsub(",", ".", df$PC2))
  
  return(df)
}
```

### 6.3 Função: Carregar Dados do Cenário

```r
load_scenario_data <- function(scenario) {
  # Carregar Excel
  df <- readxl::read_excel(
    file.path(scenario$dir, "dados_consolidados_com_scores_das_duas_PCAs.xlsx")
  )
  
  # Validar colunas PCA obrigatórias
  required_pca <- c(scenario$pc1_mag, scenario$pc2_mag,
                    scenario$pc1_var, scenario$pc2_var)
  missing_pca <- setdiff(required_pca, names(df))
  
  if (length(missing_pca) > 0) {
    stop(sprintf(
      "ERRO: Colunas PCA não encontradas: %s\nColunas disponíveis: %s",
      paste(missing_pca, collapse = ", "),
      paste(head(names(df), 10), collapse = ", ")
    ))
  }
  
  # Validar colunas de variáveis
  all_vars <- c(scenario$var_mag_cols, scenario$var_var_cols)
  missing_vars <- setdiff(all_vars, names(df))
  
  if (length(missing_vars) > 0) {
    warning(sprintf("AVISO: Variáveis não encontradas: %s",
                    paste(missing_vars, collapse = ", ")))
  }
  
  # Renomear colunas PCA para nomes padronizados (facilita código posterior)
  df$pc1_mag <- df[[scenario$pc1_mag]]
  df$pc2_mag <- df[[scenario$pc2_mag]]
  df$pc1_var <- df[[scenario$pc1_var]]
  df$pc2_var <- df[[scenario$pc2_var]]
  
  return(df)
}
```

### 6.4 Função: Calcular Variância Explicada

```r
calculate_explained_variance <- function(df, var_cols) {
  # Remover NAs das colunas relevantes
  data_clean <- na.omit(df[, var_cols, drop = FALSE])
  
  if (nrow(data_clean) < 3) {
    warning("Menos de 3 observações válidas para PCA")
    return(c(NA, NA))
  }
  
  # Executar PCA
  pca_result <- prcomp(data_clean, scale. = TRUE)
  
  # Extrair proporção de variância explicada (%)
  exp_var <- summary(pca_result)$importance[2, 1:2] * 100
  
  return(exp_var)
}
```

### 6.5 Função: Bubble Plot

```r
create_bubble_panel <- function(df, var_col, pc1_col, pc2_col,
                                 exp_var, var_label) {
  # Validar se a coluna existe
  if (!var_col %in% names(df)) {
    stop(sprintf("Coluna '%s' não encontrada no dataframe", var_col))
  }
  
  # Escalar tamanhos (min=2, max=8)
  var_values <- df[[var_col]]
  df$size_scaled <- scales::rescale(
    var_values,
    to = c(2, 8),
    from = range(var_values, na.rm = TRUE)
  )
  
  # Calcular bordas ANTES do plot (vetores, não aes)
  border_color <- ifelse(df$Arch == "inner", "black", "grey60")
  border_width <- ifelse(df$Arch == "inner", 1.5, 0.5)
  
  # Criar plot
  p <- ggplot(df, aes(x = .data[[pc1_col]], y = .data[[pc2_col]])) +
    # Linhas de referência
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
    # Pontos
    geom_point(
      aes(fill = Reef_name, shape = HAB, size = size_scaled),
      color = border_color,
      stroke = border_width,
      alpha = 0.85
    ) +
    # Escalas
    scale_fill_manual(values = reef_colors, name = "Reef") +
    scale_shape_manual(values = habitat_shapes, name = "Habitat") +
    scale_size_identity() +
    # Labels
    labs(
      title = var_label,
      x = sprintf("PC1 (%.1f%%)", exp_var[1]),
      y = sprintf("PC2 (%.1f%%)", exp_var[2])
    ) +
    # Tema
    theme_publication +
    theme(legend.position = "none") +
    coord_fixed(ratio = 1)
  
  return(p)
}
```

### 6.6 Função: Loadings Biplot

```r
create_loadings_panel <- function(loadings_df, exp_var, arrow_scale = 1.5) {
  # Validar estrutura
  if (!"variable" %in% names(loadings_df)) {
    stop("Loadings DF deve conter coluna 'variable'")
  }
  if (!all(c("PC1", "PC2") %in% names(loadings_df))) {
    stop("Loadings DF deve conter colunas 'PC1' e 'PC2'")
  }
  
  # Calcular pontas das setas
  loadings_df$xend <- loadings_df$PC1 * arrow_scale
  loadings_df$yend <- loadings_df$PC2 * arrow_scale
  
  # Aplicar labels limpos
  loadings_df$label <- sapply(loadings_df$variable, function(v) {
    v_clean <- trimws(v)
    if (v_clean %in% names(loadings_labels)) {
      return(loadings_labels[[v_clean]])
    }
    # Fallback: remover underscores e prefixos
    return(gsub("_", " ", gsub("^(sst_|dli_|chl_|cv_)", "", v_clean)))
  })
  
  # Criar plot
  p <- ggplot(loadings_df) +
    # Linhas de referência
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
    # Círculo unitário
    annotate(
      "path",
      x = cos(seq(0, 2*pi, length.out = 100)),
      y = sin(seq(0, 2*pi, length.out = 100)),
      color = "grey80",
      linetype = "dashed",
      linewidth = 0.3
    ) +
    # Setas de loadings
    geom_segment(
      aes(x = 0, y = 0, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.2, "cm"), type = "closed"),
      color = "#d62728",
      linewidth = 0.8
    ) +
    # Labels das variáveis
    geom_text(
      aes(x = xend * 1.15, y = yend * 1.15, label = label),
      size = 3.5,
      fontface = "bold"
    ) +
    # Labels dos eixos
    labs(
      title = "Loadings",
      x = sprintf("PC1 (%.1f%%)", exp_var[1]),
      y = sprintf("PC2 (%.1f%%)", exp_var[2])
    ) +
    # Tema
    theme_publication +
    coord_fixed(ratio = 1, xlim = c(-2, 2), ylim = c(-2, 2))
  
  return(p)
}
```

### 6.7 Função: Painel de Legenda

```r
create_legend_panel <- function(reefs, habs) {
  # --- Legenda de Recifes (cores) ---
  p_reef <- ggplot(
    data.frame(Reef = factor(reefs, levels = reefs)),
    aes(x = 1, y = seq_along(Reef), fill = Reef)
  ) +
    geom_point(shape = 21, size = 4, stroke = 0.5) +
    scale_fill_manual(values = reef_colors, name = "Reef") +
    theme_void() +
    theme(legend.position = "right")
  
  # --- Legenda de Habitats (formas) ---
  p_hab <- ggplot(
    data.frame(HAB = factor(habs, levels = habs)),
    aes(x = 1, y = seq_along(HAB), shape = HAB)
  ) +
    geom_point(size = 4, fill = "grey60", stroke = 0.5) +
    scale_shape_manual(values = habitat_shapes, name = "Habitat") +
    theme_void() +
    theme(legend.position = "right")
  
  # --- Legenda de Arco (manual) ---
  arch_data <- data.frame(
    x = 1,
    y = 2:1,
    label = c("Inner Arc", "Outer Arc"),
    stroke = c(1.5, 0.5),
    color = c("black", "grey60")
  )
  
  p_arch <- ggplot(arch_data, aes(x = x, y = y)) +
    geom_point(
      shape = 21, size = 4, fill = "grey70",
      color = arch_data$color,
      stroke = arch_data$stroke
    ) +
    geom_text(
      aes(x = x + 0.5, y = y, label = label),
      hjust = 0, size = 3
    ) +
    xlim(0.5, 3) +
    ylim(0.5, 2.5) +
    labs(title = "Arc") +
    theme_void() +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0))
  
  # --- Combinar legendas ---
  combined <- cowplot::plot_grid(
    cowplot::get_legend(p_reef),
    cowplot::get_legend(p_hab),
    p_arch,
    ncol = 1,
    rel_heights = c(1.2, 0.8, 1)
  )
  
  return(combined)
}
```

### 6.8 Função: Figura de Magnitude Completa

```r
create_magnitude_figure <- function(df, loadings_df, exp_var, scenario) {
  var_cols <- scenario$var_mag_cols
  var_labels <- scenario$var_mag_labels
  
  # Criar 4 bubble plots
  bubbles <- lapply(seq_along(var_cols), function(i) {
    create_bubble_panel(
      df = df,
      var_col = var_cols[i],
      pc1_col = "pc1_mag",  # Nome padronizado
      pc2_col = "pc2_mag",
      exp_var = exp_var,
      var_label = var_labels[i]
    )
  })
  
  # Criar painel de loadings
  p_loadings <- create_loadings_panel(loadings_df, exp_var)
  
  # Criar painel de legenda
  p_legend <- create_legend_panel(
    reefs = unique(df$Reef_name),
    habs = unique(df$HAB)
  )
  
  # Compor com patchwork (2 linhas × 3 colunas)
  final_figure <- (
    (bubbles[[1]] | bubbles[[2]] | bubbles[[3]]) /
    (bubbles[[4]] | p_loadings | p_legend)
  ) +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(
      tag_levels = 'a',
      tag_suffix = ')'
    )
  
  return(final_figure)
}
```

### 6.9 Função: Figura de Variabilidade Completa

```r
create_variability_figure <- function(df, loadings_df, exp_var, scenario) {
  var_cols <- scenario$var_var_cols
  var_labels <- scenario$var_var_labels
  
  # Criar 3 bubble plots
  bubbles <- lapply(seq_along(var_cols), function(i) {
    create_bubble_panel(
      df = df,
      var_col = var_cols[i],
      pc1_col = "pc1_var",  # Nome padronizado
      pc2_col = "pc2_var",
      exp_var = exp_var,
      var_label = var_labels[i]
    )
  })
  
  # Criar painel de loadings
  p_loadings <- create_loadings_panel(loadings_df, exp_var)
  
  # Criar painel de legenda
  p_legend <- create_legend_panel(
    reefs = unique(df$Reef_name),
    habs = unique(df$HAB)
  )
  
  # Compor grid principal (2×2)
  main_grid <- (
    (bubbles[[1]] | bubbles[[2]]) /
    (bubbles[[3]] | p_loadings)
  ) +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(
      tag_levels = 'a',
      tag_suffix = ')'
    )
  
  # Adicionar legenda à direita usando cowplot
  final_figure <- cowplot::plot_grid(
    main_grid,
    p_legend,
    ncol = 2,
    rel_widths = c(4, 1)
  )
  
  return(final_figure)
}
```

---

## 7. LOOP PRINCIPAL

```r
# --- Configuração inicial ---
base_dir <- "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_RESULTS"

output_dir <- file.path(base_dir, "PCA_Publication_Figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Processar cada cenário ---
for (scenario_name in names(scenarios)) {
  cat(sprintf("\n=== Processando: %s ===\n", scenario_name))
  sc <- scenarios[[scenario_name]]
  
  # 1. Carregar dados
  cat("  Carregando dados...\n")
  df <- load_scenario_data(sc)
  
  # 2. Carregar loadings
  cat("  Carregando loadings...\n")
  loadings_mag <- load_loadings(
    file.path(sc$dir, sc$loadings_mag),
    sc$csv_sep,
    sc$csv_dec
  )
  loadings_var <- load_loadings(
    file.path(sc$dir, sc$loadings_var),
    sc$csv_sep,
    sc$csv_dec
  )
  
  # 3. Calcular variância explicada
  cat("  Calculando variância explicada...\n")
  exp_var_mag <- calculate_explained_variance(df, sc$var_mag_cols)
  exp_var_var <- calculate_explained_variance(df, sc$var_var_cols)
  
  cat(sprintf("    Magnitude: PC1=%.1f%%, PC2=%.1f%%\n",
              exp_var_mag[1], exp_var_mag[2]))
  cat(sprintf("    Variability: PC1=%.1f%%, PC2=%.1f%%\n",
              exp_var_var[1], exp_var_var[2]))
  
  # 4. Gerar figura Magnitude
  cat("  Gerando figura Magnitude...\n")
  fig_mag <- create_magnitude_figure(df, loadings_mag, exp_var_mag, sc)
  
  # 5. Salvar Magnitude
  mag_base <- sprintf("Figure_PCA_Magnitude_%s", scenario_name)
  ggsave(file.path(output_dir, paste0(mag_base, ".png")),
         fig_mag, width = 14, height = 10, dpi = 300)
  ggsave(file.path(output_dir, paste0(mag_base, ".pdf")),
         fig_mag, width = 14, height = 10)
  cat(sprintf("    ✓ %s salvo\n", mag_base))
  
  # 6. Gerar figura Variability
  cat("  Gerando figura Variability...\n")
  fig_var <- create_variability_figure(df, loadings_var, exp_var_var, sc)
  
  # 7. Salvar Variability
  var_base <- sprintf("Figure_PCA_Variability_%s", scenario_name)
  ggsave(file.path(output_dir, paste0(var_base, ".png")),
         fig_var, width = 12, height = 10, dpi = 300)
  ggsave(file.path(output_dir, paste0(var_base, ".pdf")),
         fig_var, width = 12, height = 10)
  cat(sprintf("    ✓ %s salvo\n", var_base))
  
  cat(sprintf("  ✓ %s completo\n", scenario_name))
}

cat("\n=== PROCESSAMENTO CONCLUÍDO ===\n")
cat(sprintf("Figuras em: %s\n", output_dir))
```

---

## 8. CHECKLIST DE VALIDAÇÃO

### Antes da Execução
- [ ] Pacotes instalados: readxl, readr, dplyr, ggplot2, patchwork, cowplot, scales, grid
- [ ] Diretórios de entrada existem
- [ ] Arquivos Excel e CSV presentes em cada cenário

### Durante Execução
- [ ] Nenhum erro de coluna não encontrada
- [ ] Variância explicada calculada (não NA)
- [ ] 12 arquivos gerados (6 PNG + 6 PDF)

### Após Execução (Verificação Visual)
- [ ] 5 recifes com cores corretas (ARC=azul, ITA=laranja, PAB=verde, UCR=vermelho, TIM=roxo)
- [ ] 3 habitats com formas corretas (PA=círculo, RR=quadrado, TP=triângulo)
- [ ] Arcos distinguíveis (inner=borda preta grossa, outer=borda cinza fina)
- [ ] Tamanhos dos bubbles variam proporcionalmente
- [ ] Loadings com setas vermelhas apontando corretamente
- [ ] Variância explicada nos rótulos dos eixos
- [ ] Tags de painel (a, b, c, d, e)
- [ ] Sem título principal no topo da figura
- [ ] Labels em inglês

---

## 9. SAÍDAS ESPERADAS

```
#######FINAL_RESULTS/PCA_Publication_Figures/
├── Figure_PCA_Magnitude_CV_02.png     (14×10 in, 300 DPI)
├── Figure_PCA_Magnitude_CV_02.pdf     (14×10 in, vetorial)
├── Figure_PCA_Variability_CV_02.png   (12×10 in, 300 DPI)
├── Figure_PCA_Variability_CV_02.pdf   (12×10 in, vetorial)
├── Figure_PCA_Magnitude_CV_30.png
├── Figure_PCA_Magnitude_CV_30.pdf
├── Figure_PCA_Variability_CV_30.png
├── Figure_PCA_Variability_CV_30.pdf
├── Figure_PCA_Magnitude_CV_ALL.png
├── Figure_PCA_Magnitude_CV_ALL.pdf
├── Figure_PCA_Variability_CV_ALL.png
└── Figure_PCA_Variability_CV_ALL.pdf
```

**Total: 12 arquivos**

---

## 10. INFORMAÇÕES PARA O LLM IMPLEMENTADOR

### 10.1 Arquivo de Destino
```
C:/Users/rbfra/.../######FINAL/#######FINAL_CODES/PCA_Visualization_Publication.R
```

### 10.2 Ordem de Implementação
1. Carregar pacotes
2. Definir constantes (cores, formas, labels, tema)
3. Definir configuração de cenários
4. Implementar funções auxiliares (load_loadings, load_scenario_data, calculate_explained_variance)
5. Implementar funções de plotagem (create_bubble_panel, create_loadings_panel, create_legend_panel)
6. Implementar funções de composição (create_magnitude_figure, create_variability_figure)
7. Implementar loop principal
8. Testar com um cenário antes de rodar todos

### 10.3 Pontos de Atenção Críticos
1. **Nomes de colunas PCA variam entre cenários** - usar mapeamento explícito
2. **CV_ALL usa separador diferente** (`;`) - detectar ao carregar loadings
3. **Bordas devem estar FORA de `aes()`** - usar vetores pré-calculados
4. **Primeira coluna do CSV é o nome da variável** - renomear para 'variable'
5. **Recalcular variância explicada** - não está salva nos arquivos

### 10.4 Verificação Pós-Implementação
```r
# Verificar se as figuras foram geradas
list.files(output_dir, pattern = "*.png|*.pdf")
# Deve retornar 12 arquivos
```

---

*Plano Definitivo Versão FINAL (v3)*
*Data: 2026-01-12*
*Pronto para implementação por LLM*
