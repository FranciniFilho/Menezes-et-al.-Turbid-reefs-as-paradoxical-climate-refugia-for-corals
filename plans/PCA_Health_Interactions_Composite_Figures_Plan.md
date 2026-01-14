# Plano de Implementação: Figuras PCA Compostas (Saúde e Interações)
## Estilo Gemini Pro - Bubble Plots + Loadings em Figura Única

**Data**: 14/01/2026
**Objetivo**: Criar figuras compostas para PCA Saúde e PCA Interações com bubble plots das variáveis originais + painel de loadings com setas, seguindo o estilo visual de `PCA_Visualization_Publication_vGeminiPro.R`

---

## 1. Visão Geral da Arquitetura

### 1.1 Estrutura do Novo Script
```
#######FINAL_CODES/
├── PCA_Health_Interactions_Composite_vGeminiPro.R  # NOVO SCRIPT
```

### 1.2 Fluxo de Dados
```
Vitality Data (Python PCA output)
    ↓
scores_PCA_Saude.xlsx + loadings_PCA_Saude.csv
scores_PCA_Interacoes.xlsx + loadings_PCA_Interacoes.csv
    ↓
Novo Script R → Figuras Compostas (PNG 300dpi)
```

---

## 2. Especificação Visual Detalhada

### 2.1 Esquema de Cores (Igual ao GeminiPro.R)

```r
# Cores dos recifes (REEF)
reef_colors <- c(
    "ARC" = "#1f77b4",  # Azul
    "ITA" = "#ff7f0e",  # Laranja
    "PAB" = "#2ca02c",  # Verde
    "UCR" = "#d62728",  # Vermelho
    "TIM" = "#9467bd"   # Roxo
)

# Formas dos habitats (HAB)
habitat_shapes <- c(
    "PA" = 21,  # Círculo preenchido
    "RR" = 22,  # Quadrado preenchido
    "TP" = 24   # Triângulo preenchido
)
```

### 2.2 Tema de Publicação (theme_publication)

```r
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
```

### 2.3 Layout das Figuras Compostas

#### Figura PCA Saúde:
```
┌─────────────┬─────────────┬─────────────┐
│  Bubble     │  Bubble     │  Bubble     │
│  HEALTH %   │  BLEACHING  │  DEAD %     │
│     (a)     │     (b)     │     (c)     │
├─────────────┴─────────────┴─────────────┤
│  Loadings (d)              │  LEGENDA    │
│                             │             │
└────────────────────────────┴─────────────┘
```

#### Figura PCA Interações:
```
┌─────────────┬─────────────┬─────────────┐
│  Bubble     │  Bubble     │  Bubble     │
│  SUR_TURF   │  SUR_CCA    │  SUR_CYANO  │
│     (a)     │     (b)     │     (c)     │
├─────────────┼─────────────┼─────────────┤
│  Bubble     │  Bubble     │             │
│  SUR_MACRO  │  SUR_PALY   │  LEGENDA    │
│     (d)     │     (e)     │             │
├─────────────┴─────────────┴─────────────┤
│  Loadings (f)                          │
└────────────────────────────────────────┘
```

**NOTA**: A variância explicada (%) aparece nos labels dos eixos PC1 e PC2 de cada painel.

**Dimensões finais**:
- PCA Saúde: 12 x 7 polegadas
- PCA Interações: 14 x 8 polegadas
- DPI: 300

---

## 3. Estrutura de Dados de Entrada

### 3.1 Arquivo scores_PCA_Saude.xlsx
```
SITE_COL | HEALTH_PC1 | HEALTH_PC2 | SITE | HAB | REEF | HEALTH % | BLEACHING % | DEAD %
```

### 3.2 Arquivo loadings_PCA_Saude.csv
```csv
variable,PC1,PC2
HEALTH %,0.707,0.123
BLEACHING %,-0.234,0.912
DEAD %,-0.667,-0.345
```

### 3.3 Arquivo scores_PCA_Interacoes.xlsx
```
SITE | HAB | REEF | PC1_INTERACAO | PC2_INTERACAO | SUR_TURF % | SUR_CCA % | SUR_CYANO % | SUR_DICTYOTA % | SUR_OTHMACR % | SUR_PALYTHOA % | SUR_SAND % | SUR_NON-BIOTIC %
```

### 3.4 Arquivo loadings_PCA_Interacoes.csv
```csv
variable,PC1,PC2
SUR_TURF,0.523,-0.412
SUR_CCA,0.234,0.567
...
```

---

## 4. Funções Auxiliares - Implementação Detalhada

### 4.1 Função: `load_pca_data()` - Carregar dados da PCA

```r
load_pca_data <- function(base_dir, pca_type) {
    # pca_type: "Saude" ou "Interacoes"

    scores_file <- file.path(
        base_dir,
        paste0("scores_PCA_", pca_type, ".xlsx")
    )
    loadings_file <- file.path(
        base_dir,
        paste0("loadings_PCA_", pca_type, ".csv")
    )

    # Carregar scores
    df_scores <- readxl::read_excel(scores_file)

    # Carregar loadings
    df_loadings <- read.csv(loadings_file, stringsAsFactors = FALSE)

    return(list(
        scores = df_scores,
        loadings = df_loadings
    ))
}
```

### 4.2 Função: `create_bubble_panel()` - Criar painel de bubble plot

```r
create_bubble_panel <- function(df_scores, var_col, pc1_col, pc2_col,
                                exp_var, var_label,
                                reef_colors, habitat_shapes) {

    # Escalonar tamanho dos pontos (bubbles)
    var_values <- df_scores[[var_col]]
    df_scores$size_scaled <- scales::rescale(
        var_values,
        to = c(2, 8),
        from = range(var_values, na.rm = TRUE)
    )

    # Bordas mais grossas para inner arc (se aplicável)
    # Neste caso não temos Arch, então borda uniforme
    border_color <- "black"
    border_width <- 0.8

    # Calcular limites dos eixos
    vals <- c(df_scores[[pc1_col]], df_scores[[pc2_col]])
    limit <- max(abs(vals), na.rm = TRUE) * 1.1

    # Criar o bubble plot
    p <- ggplot(df_scores, aes(x = .data[[pc1_col]], y = .data[[pc2_col]])) +
        # Linhas de referência
        geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
        geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
        # Bolhas com cor por Reef, forma por HAB, tamanho por variável
        geom_point(aes(fill = REEF, shape = HAB, size = size_scaled),
                   color = border_color,
                   stroke = border_width,
                   alpha = 0.8) +
        # Escalas
        scale_fill_manual(values = reef_colors, name = "Reef") +
        scale_shape_manual(values = habitat_shapes, name = "Habitat") +
        scale_size_identity() +  # Usar valores de tamanho diretamente
        # Labels e títulos
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

    return(p)
}
```

### 4.3 Função: `create_loadings_panel()` - Criar painel de loadings

```r
create_loadings_panel <- function(loadings_df, exp_var, arrow_scale = 1.5) {

    # Calcular coordenadas finais das setas
    loadings_df$xend <- loadings_df$PC1 * arrow_scale
    loadings_df$yend <- loadings_df$PC2 * arrow_scale

    # Limpar labels das variáveis
    loadings_df$label <- sapply(loadings_df$variable, function(v) {
        v_clean <- trimws(v)
        # Remover sufixos de porcentagem
        v_clean <- gsub(" %$", "", v_clean)
        # Substituir underscores por espaços
        v_clean <- gsub("_", " ", v_clean)
        # Capitalizar
        v_clean <- tools::toTitleCase(v_clean)
        return(v_clean)
    })

    limit <- 2.2

    p <- ggplot(loadings_df) +
        # Linhas de referência
        geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
        geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
        # Círculo de referência unitário
        annotate("path",
                 x = cos(seq(0, 2 * pi, length.out = 100)),
                 y = sin(seq(0, 2 * pi, length.out = 100)),
                 color = "grey80",
                 linetype = "dashed",
                 linewidth = 0.3) +
        # Setas de loadings
        geom_segment(aes(x = 0, y = 0, xend = xend, yend = yend),
                     arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
                     color = "#d62728",  # Vermelho
                     linewidth = 0.7) +
        # Labels das variáveis
        geom_text(aes(x = xend * 1.2, y = yend * 1.2, label = label),
                  size = 4,
                  fontface = "bold") +
        # Títulos
        labs(title = "Loadings",
             x = sprintf("PC1 (%.1f%%)", exp_var[1]),
             y = sprintf("PC2 (%.1f%%)", exp_var[2])) +
        xlim(-limit, limit) +
        ylim(-limit, limit) +
        theme_publication +
        coord_fixed(ratio = 1)

    return(p)
}
```

### 4.4 Função: `create_legend_strip()` - Criar legenda unificada

```r
create_legend_strip <- function(reefs, habs, size_scale = 1.0) {

    legend_theme <- theme(
        legend.title = element_text(size = 14 * size_scale, face = "bold"),
        legend.text = element_text(size = 12 * size_scale),
        legend.key.size = unit(1.0 * size_scale, "cm"),
        legend.box.margin = margin(0, 0, 0, 0)
    )

    # Legenda de Reef (cores)
    p_reef <- ggplot(
        data.frame(REEF = factor(reefs, levels = reefs)),
        aes(x = 1, y = 1, fill = REEF)
    ) +
        geom_point(shape = 21, size = 5 * size_scale) +
        scale_fill_manual(values = reef_colors, name = "Reef") +
        theme_publication +
        legend_theme
    leg_reef <- cowplot::get_legend(p_reef)

    # Legenda de Habitat (formas)
    p_hab <- ggplot(
        data.frame(HAB = factor(habs, levels = habs)),
        aes(x = 1, y = 1, shape = HAB)
    ) +
        geom_point(size = 5 * size_scale, fill = "grey60") +
        scale_shape_manual(values = habitat_shapes, name = "Habitat") +
        theme_publication +
        legend_theme
    leg_hab <- cowplot::get_legend(p_hab)

    # Combinar legendas horizontalmente
    legend_strip <- plot_grid(
        leg_reef,
        leg_hab,
        ncol = 2,
        rel_widths = c(1, 0.8),
        align = "vh"
    )

    return(legend_strip)
}
```

### 4.6 Função: `calculate_explained_variance()` - Calcular variância explicada

```r
calculate_explained_variance <- function(df_scores, loadings_df) {

    # A variância explicada pode ser calculada a partir dos autovalores
    # Ou assumimos que está armazenada em algum lugar
    # Alternativa: recalcular PCA a partir dos dados originais

    # Para simplificar, podemos extrair dos eigenvalues se disponíveis
    # Ou usar a proporção da soma dos quadrados dos loadings

    # Método: Usar a soma dos quadrados dos loadings normalizados
    eigenvalues <- colSums(loadings_df[, c("PC1", "PC2")]^2)
    exp_var <- (eigenvalues / sum(eigenvalues)) * 100

    return(exp_var)
}
```

---

## 5. Funções Principais de Montagem

### 5.1 Função: `create_health_pca_figure()`

```r
create_health_pca_figure <- function(data_dir, output_dir) {

    cat("\n=== Criando Figura PCA Saúde ===\n")

    # 1. Carregar dados
    pca_data <- load_pca_data(data_dir, "Saude")
    df_scores <- pca_data$scores
    df_loadings <- pca_data$loadings

    # 2. Calcular variância explicada
    exp_var <- calculate_explained_variance(df_scores, df_loadings)

    # 3. Definir variáveis para bubble plots
    bubble_vars <- c("HEALTH %", "BLEACHING %", "DEAD %")
    bubble_labels <- c("Healthy Tissue", "Bleaching", "Dead Tissue")

    # 4. Criar painéis de bubble plots
    cat("  - Criando painéis de bubble plots...\n")
    bubble_panels <- list()
    for (i in seq_along(bubble_vars)) {
        bubble_panels[[i]] <- create_bubble_panel(
            df_scores = df_scores,
            var_col = bubble_vars[i],
            pc1_col = "HEALTH_PC1",
            pc2_col = "HEALTH_PC2",
            exp_var = exp_var,
            var_label = bubble_labels[i],
            reef_colors = reef_colors,
            habitat_shapes = habitat_shapes
        )
    }

    # 5. Criar painel de loadings
    cat("  - Criando painel de loadings...\n")
    loadings_panel <- create_loadings_panel(df_loadings, exp_var)

    # 6. Criar legenda
    cat("  - Criando legenda...\n")
    legend_strip <- create_legend_strip(
        unique(df_scores$REEF),
        unique(df_scores$HAB)
    )

    # 7. Montar figura composta
    cat("  - Montando figura composta...\n")

    # Layout usando patchwork
    # Top row: 3 bubble plots
    # Bottom row: loadings | legend
    layout_design <- "
    ABC
    D#E
    "

    figure_composite <- bubble_panels[[1]] + bubble_panels[[2]] +
        bubble_panels[[3]] + loadings_panel + plot_spacer() + legend_strip +
        plot_layout(
            design = layout_design,
            widths = c(1, 1, 1)
        ) +
        plot_annotation(tag_levels = "a", tag_suffix = ")")

    # 8. Salvar figura
    output_path <- file.path(output_dir, "PCA_Saude_Figura_Composta_GeminiPro.png")
    ggsave(output_path, figure_composite,
           width = 12, height = 7, dpi = 300)

    cat(sprintf("  ✓ Figura salva: %s\n", output_path))

    return(figure_composite)
}
```

### 5.2 Função: `create_interactions_pca_figure()`

```r
create_interactions_pca_figure <- function(data_dir, output_dir) {

    cat("\n=== Criando Figura PCA Interações ===\n")

    # 1. Carregar dados
    pca_data <- load_pca_data(data_dir, "Interacoes")
    df_scores <- pca_data$scores
    df_loadings <- pca_data$loadings

    # 2. Calcular variância explicada
    exp_var <- calculate_explained_variance(df_scores, df_loadings)

    # 3. Definir variáveis para bubble plots
    # Usar as variáveis originais mais importantes
    bubble_vars <- c("SUR_TURF %", "SUR_CCA %", "SUR_CYANO %",
                     "SUR_DICTYOTA %", "SUR_OTHMACR %")

    # Combinar macroalgae para o plot
    df_scores$SUR_MACRO_ALGAE <-
        df_scores$`SUR_DICTYOTA %` %na% 0 + df_scores$`SUR_OTHMACR %` %na% 0

    bubble_vars_plot <- c("SUR_TURF %", "SUR_CCA %", "SUR_CYANO %",
                          "SUR_MACRO_ALGAE")

    bubble_labels <- c("Turf", "Crustose Coralline Algae", "Cyanobacteria",
                       "Macroalgae")

    # Criar cópia com nomes limpos para plotagem
    df_plot <- df_scores
    names(df_plot)[names(df_plot) == "SUR_TURF %"] <- "SUR_TURF"
    names(df_plot)[names(df_plot) == "SUR_CCA %"] <- "SUR_CCA"
    names(df_plot)[names(df_plot) == "SUR_CYANO %"] <- "SUR_CYANO"

    bubble_vars_final <- c("SUR_TURF", "SUR_CCA", "SUR_CYANO", "SUR_MACRO_ALGAE")

    # 4. Criar painéis de bubble plots
    cat("  - Criando painéis de bubble plots...\n")
    bubble_panels <- list()
    for (i in seq_along(bubble_vars_final)) {
        bubble_panels[[i]] <- create_bubble_panel(
            df_scores = df_plot,
            var_col = bubble_vars_final[i],
            pc1_col = "PC1_INTERACAO",
            pc2_col = "PC2_INTERACAO",
            exp_var = exp_var,
            var_label = bubble_labels[i],
            reef_colors = reef_colors,
            habitat_shapes = habitat_shapes
        )
    }

    # 5. Criar painel de loadings
    cat("  - Criando painel de loadings...\n")
    loadings_panel <- create_loadings_panel(df_loadings, exp_var)

    # 6. Criar legenda
    cat("  - Criando legenda...\n")
    legend_strip <- create_legend_strip(
        unique(df_scores$REEF),
        unique(df_scores$HAB)
    )

    # 7. Montar figura composta
    cat("  - Montando figura composta...\n")

    # Layout: 2 linhas de bubble plots + linha com loadings
    layout_design <- "
    ABC
    DE#
    FFF
    "

    figure_composite <- bubble_panels[[1]] + bubble_panels[[2]] +
        bubble_panels[[3]] + bubble_panels[[4]] +
        plot_spacer() + legend_strip +
        loadings_panel +
        plot_layout(
            design = layout_design,
            widths = c(1, 1, 1)
        ) +
        plot_annotation(tag_levels = "a", tag_suffix = ")")

    # 8. Salvar figura
    output_path <- file.path(output_dir, "PCA_Interacoes_Figura_Composta_GeminiPro.png")
    ggsave(output_path, figure_composite,
           width = 14, height = 8, dpi = 300)

    cat(sprintf("  ✓ Figura salva: %s\n", output_path))

    return(figure_composite)
}
```

---

## 6. Script Principal - Estrutura Completa

```r
# ============================================================================
# PCA_Health_Interactions_Composite_vGeminiPro.R
# Geração de figuras PCA compostas para Saúde e Interações
# Estilo: GeminiPro (bubble plots + loadings)
#
# Layout:
#   PCA Saúde: 3 bubble plots + loadings + legenda
#   PCA Interações: 4-5 bubble plots + loadings + legenda
#   (Variância explicada aparece nos eixos de cada painel)
# ============================================================================

# --- 1. PACOTES ---
libs <- c(
    "readxl",      # Ler Excel
    "readr",       # Ler CSV
    "dplyr",       # Manipulação de dados
    "ggplot2",     # Plotagem base
    "patchwork",   # Composição de painéis
    "cowplot",     # Extração de legendas
    "scales",      # Rescale de tamanhos
    "grid",        # arrow() e unit()
    "tools"        # toTitleCase()
)

# Carregar pacotes
invisible(lapply(libs, function(lib) {
    if (!require(lib, character.only = TRUE)) {
        install.packages(lib)
        library(lib, character.only = TRUE)
    }
}))

# --- 2. CAMINHOS ---
base_dir <- "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_RESULTS/output_ANALISE_BIOLOGICA_boxplot_PCAnew"
output_dir <- file.path(base_dir, "PCA_Composite_Figures_vGeminiPro")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 3. ESQUEMA VISUAL (reef_colors, habitat_shapes, theme_publication) ---
# [Inserir código da seção 2.2 acima]

# --- 4. FUNÇÕES AUXILIARES ---
# [Inserir todas as funções das seções 4.1 a 4.6 acima]

# --- 5. EXECUÇÃO PRINCIPAL ---
main <- function() {

    cat("\n")
    cat("========================================\n")
    cat("PCA Composite Figures Generator\n")
    cat("Estilo: GeminiPro\n")
    cat("========================================\n")

    # Verificar se os dados existem
    if (!file.exists(file.path(base_dir, "scores_PCA_Saude.xlsx"))) {
        stop("ERRO: scores_PCA_Saude.xlsx não encontrado. Execute primeiro Calculate_RGR_&_health_PCA.py")
    }

    # Criar figura PCA Saúde
    create_health_pca_figure(base_dir, output_dir)

    # Criar figura PCA Interações
    create_interactions_pca_figure(base_dir, output_dir)

    cat("\n=== PROCESSAMENTO CONCLUÍDO ===\n")
    cat(sprintf("Figuras salvas em: %s\n", output_dir))
}

# Executar
main()
```

---

## 7. Checklist de Implementação

- [ ] **Passo 1**: Criar arquivo `PCA_Health_Interactions_Composite_vGeminiPro.R`
- [ ] **Passo 2**: Implementar todas as funções auxiliares
  - [ ] `load_pca_data()`
  - [ ] `create_bubble_panel()`
  - [ ] `create_loadings_panel()`
  - [ ] `create_legend_strip()`
  - [ ] `calculate_explained_variance()`
- [ ] **Passo 3**: Implementar funções principais de montagem
  - [ ] `create_health_pca_figure()`
  - [ ] `create_interactions_pca_figure()`
- [ ] **Passo 4**: Implementar script principal com `main()`
- [ ] **Passo 5**: Testar executando o script
- [ ] **Passo 6**: Verificar saída das figuras
- [ ] **Passo 7**: Ajustar dimensões/layout se necessário

---

## 8. Dependências R

```r
# Versões mínimas sugeridas
# R >= 4.2.0
# ggplot2 >= 3.4.0
# patchwork >= 1.1.2
# cowplot >= 1.1.1
# dplyr >= 1.1.0
# readxl >= 1.4.2
# scales >= 1.2.1
```

Instalar se necessário:
```r
install.packages(c("readxl", "dplyr", "ggplot2", "patchwork",
                   "cowplot", "scales", "tools"))
```

---

## 9. Notas Importantes

1. **Pré-requisito**: O script Python `Calculate_RGR_&_health_PCA.py` deve ser executado primeiro para gerar os arquivos de entrada

2. **Tratamento de valores NA**: Usar `%na% 0` ou similar para lidar com valores faltantes

3. **Coordenadas fixas**: `coord_fixed(ratio = 1)` garante que os ângulos dos loadings sejam representados corretamente

4. **Tamanho dos pontos**: O rescale para (2, 8) garante bolhas visíveis mas não excessivamente grandes

5. **Cores das setas**: Vermelho (`#d62728`) para destacar sobre os pontos coloridos

---

## 10. Exemplo de Uso

```r
# No console R:
source("#######FINAL_CODES/PCA_Health_Interactions_Composite_vGeminiPro.R")

# Ou via linha de comando:
# Rscript PCA_Health_Interactions_Composite_vGeminiPro.R
```

---

**Fim do Plano de Implementação**
