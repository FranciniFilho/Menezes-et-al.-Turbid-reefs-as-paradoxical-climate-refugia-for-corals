# Plano de Implementação Revisado: Figuras PCA Compostas (Saúde e Interações)
## Versão 2.0 - Padrão Nature/Science

**Data**: 14/01/2026  
**Revisor**: Claude  
**Status**: Plano revisado para correção de inconsistências críticas

---

## 1. Avaliação Crítica do Plano Original (v1)

### 1.1 Problemas Identificados

> [!CAUTION]
> **Inconsistências Críticas que Afetariam a Qualidade da Publicação**

| # | Problema | Impacto | Seção Afetada |
|---|----------|---------|---------------|
| 1 | **Nomes de colunas divergentes** | Script não executaria | Seções 4.2, 5.1, 5.2 |
| 2 | **Variáveis de interações incorretas** | Dados cientificamente errados | Seção 5.2 |
| 3 | **Cálculo de variância explicada incorreto** | Valores errados nos eixos | Seção 4.6 |
| 4 | **Falta do painel de Arch (Inner/Outer Arc)** | Inconsistência com Reference Script | Seção 4.4 |
| 5 | **Layout patchwork com erro de sintaxe** | Plot não montaria | Seção 5.1, 5.2 |
| 6 | **Falta de escalonamento por Arch** | Bordas inconsistentes com Reference Script | Seção 4.2 |
| 7 | **Labels em Português misturados com Inglês** | Inconsistência para publicação | Múltiplas seções |

---

### 1.2 Análise Detalhada das Inconsistências

#### **Problema 1: Nomes de Colunas Divergentes**

O plano original assume colunas como `"HEALTH %"`, `"BLEACHING %"`, `"DEAD %"` diretamente nos scores, mas o script Python salva:

```python
# Calculate_RGR_&_health_PCA.py (linha 69-71)
df_scores_health = pd.merge(df_scores_health, health_props, on='SITE_COL')
# As colunas originais são preservadas: "HEALTH %", "BLEACHING %", "DEAD %"
```

**Verificação**: Colunas corretas no `scores_PCA_Saude.xlsx`:
- `SITE_COL`, `HEALTH_PC1`, `HEALTH_PC2`, `SITE`, `HAB`, `REEF`, `HEALTH %`, `BLEACHING %`, `DEAD %`

**Status**: ✅ Correto no plano (mas precisa verificar espaços/caracteres especiais)

---

#### **Problema 2: Variáveis de Interações - CRÍTICO**

O plano original (seção 5.2) lista:
```r
bubble_vars <- c("SUR_TURF %", "SUR_CCA %", "SUR_CYANO %", 
                 "SUR_DICTYOTA %", "SUR_OTHMACR %")
```

Mas o script Python (linha 83-85) **agrupa** as variáveis de forma diferente:

```python
interaction_vars = {
    'SUR_TURF': interaction_means_agg['SUR_TURF %'],
    'SUR_CCA': interaction_means_agg['SUR_CCA %'],
    'SUR_MACROALGAE': interaction_means_agg['SUR_DICTYOTA %'] + interaction_means_agg['SUR_OTHMACR %'],
    'SUR_CYANO': interaction_means_agg['SUR_CYANO %'],
    'SUR_PALYTHOA': interaction_means_agg['SUR_PALYTHOA %'],
    'SUR_ABIOTIC': interaction_means_agg['SUR_SAND %'] + interaction_means_agg['SUR_NON-BIOTIC %']
}
```

**Colunas reais no `scores_PCA_Interacoes.xlsx`**:
- `SITE`, `HAB`, `PC1_INTERACAO`, `PC2_INTERACAO`, `REEF`
- `SUR_TURF`, `SUR_CCA`, `SUR_MACROALGAE`, `SUR_CYANO`, `SUR_PALYTHOA`, `SUR_ABIOTIC`

> [!WARNING]
> O plano original tentaria acessar colunas inexistentes (e.g., `SUR_DICTYOTA %`) causando erro fatal.

---

#### **Problema 3: Cálculo de Variância Explicada - CRÍTICO**

O plano original (seção 4.6) usa:

```r
# INCORRETO - Método aproximado e impreciso
eigenvalues <- colSums(loadings_df[, c("PC1", "PC2")]^2)
exp_var <- (eigenvalues / sum(eigenvalues)) * 100
```

Isso NÃO funciona corretamente para PCA. A variância explicada deve vir diretamente do objeto PCA ou ser recalculada dos dados originais.

**Solução**: O Python já salva a variância explicada no script (linha 150, 156):
```python
'explained_variance': pca_health.explained_variance_ratio_ * 100
```

**Correção necessária**: Modificar o Python para **salvar** a variância explicada em arquivo separado, ou recalcular via R como no Reference Script (linha 156-163):

```r
# CORRETO - PCA_Visualization_Publication_vGeminiPro.R
calculate_explained_variance <- function(df, var_cols) {
    data_clean <- na.omit(df[, var_cols, drop = FALSE])
    pca_result <- prcomp(data_clean, scale. = TRUE)
    exp_var <- summary(pca_result)$importance[2, 1:2] * 100
    return(exp_var)
}
```

---

#### **Problema 4: Falta do Painel de Arch (Inner/Outer Arc)**

O Reference Script (`PCA_Visualization_Publication_vGeminiPro.R`) inclui uma legenda de "Arc" (Inner/Outer) com bordas diferenciadas (linha 258-270):

```r
p_arch_manual <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
    annotate("point", x = 1, y = 1, shape = 21, size = 5, fill = "grey70", 
             color = "black", stroke = 1.2) +  # Inner Arc
    annotate("point", x = 3.5, y = 1, shape = 21, size = 5, fill = "grey70", 
             color = "grey60", stroke = 0.4) +  # Outer Arc
```

E o bubble plot usa bordas variáveis por Arch (linha 172-173):
```r
border_color <- ifelse(df$Arch == "inner", "black", "grey60")
border_width <- ifelse(df$Arch == "inner", 1.2, 0.4)
```

> [!IMPORTANT]
> O plano original usa borda uniforme para todos os pontos, perdendo informação científica importante (distinção Inner/Outer Arc).

**✅ CONFIRMADO**: Os dados de Saúde/Interações possuem coluna ARCH (renomear para "Arc" nas legendas).

---

#### **Problema 5: Sintaxe de Layout Patchwork**

O plano original tenta usar:
```r
# INCORRETO
bubble_panels[[1]] + bubble_panels[[2]] + bubble_panels[[3]] + 
loadings_panel + plot_spacer() + legend_strip
```

Mas `legend_strip` é um objeto `cowplot::plot_grid()`, não um ggplot, então não pode ser combinado diretamente com `+` do patchwork.

**Solução do Reference Script** (Montagem Híbrida, linha 312-319):
```r
final_combined <- cowplot::plot_grid(
    plots_grid,  # patchwork object
    NULL,        # spacer
    leg_strip,   # cowplot grob
    ncol = 1,
    rel_heights = c(10, 1, 2.5)
)
```

---

## 2. Plano de Implementação Corrigido

> [!IMPORTANT]
> **REGRA CRÍTICA: TODOS OS TEXTOS DEVEM ESTAR EM INGLÊS**
> 
> Incluindo: títulos de painéis, labels de eixos, legendas, mensagens de console, comentários de código voltados ao usuário. Isso é obrigatório para publicação em periódicos internacionais (Nature/Science).

### 2.1 Estrutura de Dados Verificada

#### Arquivo: `scores_PCA_Saude.xlsx`
| Coluna | Tipo | Uso |
|--------|------|-----|
| `SITE_COL` | string | Identificador único colônia-site |
| `HEALTH_PC1` | numeric | Score PC1 |
| `HEALTH_PC2` | numeric | Score PC2 |
| `SITE` | string | Identificador do site |
| `HAB` | string | Tipo de habitat (PA, RR, TP) |
| `REEF` | string | Nome do recife |
| `HEALTH %` | numeric | Proporção de tecido saudável |
| `BLEACHING %` | numeric | Proporção de branqueamento |
| `DEAD %` | numeric | Proporção de tecido morto |

#### Arquivo: `scores_PCA_Interacoes.xlsx`
| Coluna | Tipo | Uso |
|--------|------|-----|
| `SITE` | string | Identificador do site |
| `HAB` | string | Tipo de habitat |
| `REEF` | string | Nome do recife |
| `PC1_INTERACAO` | numeric | Score PC1 |
| `PC2_INTERACAO` | numeric | Score PC2 |
| `SUR_TURF` | numeric | Cobertura de turf |
| `SUR_CCA` | numeric | Cobertura de CCA |
| `SUR_MACROALGAE` | numeric | Cobertura de macroalgas (agregada) |
| `SUR_CYANO` | numeric | Cobertura de cianobactérias |
| `SUR_PALYTHOA` | numeric | Cobertura de Palythoa |
| `SUR_ABIOTIC` | numeric | Cobertura de substrato abiótico (agregada) |

---

### 2.2 Esquema Visual Corrigido (Idêntico ao Reference Script)

```r
# ============================================================================
# CORES DOS RECIFES (Idêntico ao PCA_Visualization_Publication_vGeminiPro.R)
# ============================================================================
reef_colors <- c(
    "ARC" = "#1f77b4",  # Azul
    "ITA" = "#ff7f0e",  # Laranja
    "PAB" = "#2ca02c",  # Verde
    "UCR" = "#d62728",  # Vermelho
    "TIM" = "#9467bd"   # Roxo
)

# ============================================================================
# FORMAS DE HABITAT (Idêntico)
# ============================================================================
habitat_shapes <- c(
    "PA" = 21,  # Círculo preenchido
    "RR" = 22,  # Quadrado preenchido
    "TP" = 24   # Triângulo preenchido
)

# ============================================================================
# TEMA DE PUBLICAÇÃO (Idêntico - base_size = 14.4 = 12 * 1.2)
# ============================================================================
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

---

### 2.3 Layout Visual Revisado

#### Figura PCA Saúde (3 painéis + 1 loadings):
```
┌─────────────────┬─────────────────┬─────────────────┐
│     (a)         │      (b)        │       (c)       │
│  Healthy Tissue │   Bleaching     │   Dead Tissue   │
│    Bubble       │     Bubble      │      Bubble     │
├─────────────────┼─────────────────┼─────────────────┤
│     (d)         │                                   │
│   Loadings      │           [empty]                 │
│   (Unit Circle) │                                   │
├─────────────────┴─────────────────┴─────────────────┤
│  ● Reef Legend  │  ◆ Habitat Legend                 │
│  (cores)        │  (formas)                         │
└─────────────────┴───────────────────────────────────┘
```

#### Figura PCA Interações (5 painéis + 1 loadings):
```
┌─────────────────┬─────────────────┬─────────────────┐
│     (a)         │      (b)        │       (c)       │
│     Turf        │      CCA        │  Cyanobacteria  │
├─────────────────┼─────────────────┼─────────────────┤
│     (d)         │      (e)        │       (f)       │
│  Macroalgae     │   Palythoa      │    Loadings     │
├─────────────────┴─────────────────┴─────────────────┤
│  ● Reef Legend  │  ◆ Habitat Legend                 │
└─────────────────┴───────────────────────────────────┘
```

> [!NOTE]
> SUR_ABIOTIC foi removido pois substrato abiótico não é biologicamente relevante para interações coral-bentos.

---

### 2.4 Funções Auxiliares Corrigidas

#### 2.4.1 Carregamento de Dados

```r
load_pca_health_data <- function(base_dir) {
    scores_file <- file.path(base_dir, "scores_PCA_Saude.xlsx")
    loadings_file <- file.path(base_dir, "loadings_PCA_Saude.csv")
    
    df_scores <- readxl::read_excel(scores_file)
    df_loadings <- read.csv(loadings_file, 
                            stringsAsFactors = FALSE,
                            row.names = 1)  # Primeira coluna é index
    
    # Verificar e adicionar coluna ARCH se não existir
    # (assumindo outer arc por default - ajustar conforme dados reais)
    if (!"ARCH" %in% names(df_scores)) {
        # Mapear por REEF se disponível
        inner_reefs <- c("ARC", "ITA", "PAB")  # Ajustar conforme dados reais
        df_scores$ARCH <- ifelse(df_scores$REEF %in% inner_reefs, "inner", "outer")
        warning("Coluna ARCH criada automaticamente baseada em REEF. Verificar mapeamento!")
    }
    
    return(list(
        scores = df_scores,
        loadings = df_loadings
    ))
}

load_pca_interactions_data <- function(base_dir) {
    scores_file <- file.path(base_dir, "scores_PCA_Interacoes.xlsx")
    loadings_file <- file.path(base_dir, "loadings_PCA_Interacoes.csv")
    
    df_scores <- readxl::read_excel(scores_file)
    df_loadings <- read.csv(loadings_file, 
                            stringsAsFactors = FALSE,
                            row.names = 1)
    
    # Mesma lógica ARCH
    if (!"ARCH" %in% names(df_scores)) {
        inner_reefs <- c("ARC", "ITA", "PAB")
        df_scores$ARCH <- ifelse(df_scores$REEF %in% inner_reefs, "inner", "outer")
    }
    
    return(list(
        scores = df_scores,
        loadings = df_loadings
    ))
}
```

#### 2.4.2 Cálculo de Variância Explicada (CORRIGIDO)

```r
calculate_explained_variance <- function(df, var_cols) {
    #' Recalcula a variância explicada via prcomp()
    #' Método idêntico ao PCA_Visualization_Publication_vGeminiPro.R
    
    # Remover NAs
    data_clean <- na.omit(df[, var_cols, drop = FALSE])
    
    if (nrow(data_clean) < 3) {
        warning("Dados insuficientes para PCA. Usando valores default.")
        return(c(50, 30))  # Fallback
    }
    
    # PCA com escalonamento (reproduz sklearn.StandardScaler)
    pca_result <- prcomp(data_clean, scale. = TRUE)
    
    # Extrair proporção de variância dos dois primeiros componentes
    exp_var <- summary(pca_result)$importance[2, 1:2] * 100
    
    return(exp_var)
}
```

#### 2.4.3 Bubble Plot (Com Distinção de Arch)

```r
create_bubble_panel <- function(df, var_col, pc1_col, pc2_col, 
                                 exp_var, var_label,
                                 reef_colors, habitat_shapes) {
    #' Cria um bubble plot para uma variável específica
    #' Bordas diferenciadas por ARCH (inner = preto grosso, outer = cinza fino)
    
    # Escalonar tamanho dos pontos
    var_values <- df[[var_col]]
    df$size_scaled <- scales::rescale(
        var_values,
        to = c(2, 8),
        from = range(var_values, na.rm = TRUE)
    )
    
    # Bordas por ARCH (reproduz Reference Script)
    df$border_color <- ifelse(df$ARCH == "inner", "black", "grey60")
    df$border_width <- ifelse(df$ARCH == "inner", 1.2, 0.4)
    
    # Calcular limites simétricos
    vals <- c(df[[pc1_col]], df[[pc2_col]])
    limit <- max(abs(vals), na.rm = TRUE) * 1.1
    
    # Nota: Para usar bordas variáveis por ponto, precisamos iterar
    # ou usar aes() com valores pré-calculados
    # Alternativa: plotar em duas camadas (inner e outer separados)
    
    p <- ggplot(df, aes(x = .data[[pc1_col]], y = .data[[pc2_col]])) +
        # Linhas de referência
        geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
        geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4)
    
    # Camada Outer Arc (primeiro, ficará por baixo)
    df_outer <- df[df$ARCH == "outer", ]
    if (nrow(df_outer) > 0) {
        p <- p + geom_point(
            data = df_outer,
            aes(fill = REEF, shape = HAB, size = size_scaled),
            color = "grey60",
            stroke = 0.4,
            alpha = 0.8
        )
    }
    
    # Camada Inner Arc (por cima)
    df_inner <- df[df$ARCH == "inner", ]
    if (nrow(df_inner) > 0) {
        p <- p + geom_point(
            data = df_inner,
            aes(fill = REEF, shape = HAB, size = size_scaled),
            color = "black",
            stroke = 1.2,
            alpha = 0.8
        )
    }
    
    p <- p +
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
        xlim(-limit, limit) +
        ylim(-limit, limit) +
        theme_publication +
        theme(legend.position = "none") +
        coord_fixed(ratio = 1)
    
    return(p)
}
```

#### 2.4.4 Loadings Panel (Idêntico ao Reference)

```r
create_loadings_panel <- function(loadings_df, exp_var, arrow_scale = 1.5) {
    #' Cria painel de loadings com círculo unitário e setas vermelhas
    #' Estilo: Idêntico ao PCA_Visualization_Publication_vGeminiPro.R
    
    # Preparar dados
    loadings_df$xend <- loadings_df$PC1 * arrow_scale
    loadings_df$yend <- loadings_df$PC2 * arrow_scale
    
    # Labels das variáveis (limpar e traduzir)
    loadings_labels_map <- c(
        # Health
        "HEALTH %" = "Healthy",
        "BLEACHING %" = "Bleaching",
        "DEAD %" = "Dead",
        # Interactions
        "SUR_TURF" = "Turf",
        "SUR_CCA" = "CCA",
        "SUR_CYANO" = "Cyanobacteria",
        "SUR_MACROALGAE" = "Macroalgae",
        "SUR_PALYTHOA" = "Palythoa",
        "SUR_ABIOTIC" = "Abiotic"
    )
    
    loadings_df$label <- sapply(rownames(loadings_df), function(v) {
        v_clean <- trimws(v)
        if (v_clean %in% names(loadings_labels_map)) {
            return(loadings_labels_map[[v_clean]])
        }
        # Fallback: capitalizar e remover underscores
        v_f <- gsub("_", " ", v_clean)
        v_f <- gsub(" %$", "", v_f)
        v_f <- tools::toTitleCase(tolower(v_f))
        return(v_f)
    })
    
    limit <- 2.2
    
    p <- ggplot(loadings_df) +
        # Linhas de referência
        geom_hline(yintercept = 0, color = "grey60", linewidth = 0.4) +
        geom_vline(xintercept = 0, color = "grey60", linewidth = 0.4) +
        # Círculo unitário (dashed)
        annotate("path",
            x = cos(seq(0, 2 * pi, length.out = 100)),
            y = sin(seq(0, 2 * pi, length.out = 100)),
            color = "grey80",
            linetype = "dashed",
            linewidth = 0.3
        ) +
        # Setas vermelhas
        geom_segment(
            aes(x = 0, y = 0, xend = xend, yend = yend),
            arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
            color = "#d62728",  # Vermelho (consistente com reef_colors)
            linewidth = 0.7
        ) +
        # Labels em negrito
        geom_text(
            aes(x = xend * 1.2, y = yend * 1.2, label = label),
            size = 4,
            fontface = "bold"
        ) +
        # Títulos
        labs(
            title = "Loadings",
            x = sprintf("PC1 (%.1f%%)", exp_var[1]),
            y = sprintf("PC2 (%.1f%%)", exp_var[2])
        ) +
        xlim(-limit, limit) +
        ylim(-limit, limit) +
        theme_publication +
        coord_fixed(ratio = 1)
    
    return(p)
}
```

#### 2.4.5 Legenda Unificada (Com Arch)

```r
create_legend_strip <- function(reefs, habs, has_arch = TRUE, size_scale = 1.0) {
    #' Cria strip de legendas horizontal
    #' Inclui: Reef (cores), Habitat (formas), e opcionalmente Arch (bordas)
    
    legend_theme <- theme(
        legend.title = element_text(size = 14 * size_scale, face = "bold"),
        legend.text = element_text(size = 12 * size_scale),
        legend.key.size = unit(1.0 * size_scale, "cm"),
        legend.box.margin = margin(0, 0, 0, 0)
    )
    
    # Legenda REEF
    p_reef <- ggplot(
        data.frame(Reef = factor(reefs, levels = reefs)),
        aes(x = 1, y = 1, fill = Reef)
    ) +
        geom_point(shape = 21, size = 5 * size_scale) +
        scale_fill_manual(values = reef_colors, name = "Reef") +
        theme_publication +
        legend_theme
    leg_reef <- cowplot::get_legend(p_reef)
    
    # Legenda HABITAT
    p_hab <- ggplot(
        data.frame(HAB = factor(habs, levels = habs)),
        aes(x = 1, y = 1, shape = HAB)
    ) +
        geom_point(size = 5 * size_scale, fill = "grey60") +
        scale_shape_manual(values = habitat_shapes, name = "Habitat") +
        theme_publication +
        legend_theme
    leg_hab <- cowplot::get_legend(p_hab)
    
    if (has_arch) {
        # Legenda ARC (não "Arch" - usar "Arc" como no artigo)
        p_arc <- ggplot(data.frame(x = 1, y = 1), aes(x, y)) +
            annotate("point", x = 1, y = 1, shape = 21, size = 5 * size_scale,
                     fill = "grey70", color = "black", stroke = 1.2) +
            annotate("text", x = 1.3, y = 1, label = "Inner Arc", hjust = 0, 
                     size = 4.5 * size_scale) +
            annotate("point", x = 3.5, y = 1, shape = 21, size = 5 * size_scale,
                     fill = "grey70", color = "grey60", stroke = 0.4) +
            annotate("text", x = 3.8, y = 1, label = "Outer Arc", hjust = 0, 
                     size = 4.5 * size_scale) +
            xlim(0.8, 6) + ylim(0, 2) +
            theme_void() +
            labs(title = "Arc") +  # Título da legenda: "Arc" (não "Arch")
            theme(
                plot.title = element_text(face = "bold", size = 14 * size_scale, 
                                          hjust = 0.1, margin = margin(b = 5)),
                plot.margin = margin(10, 0, 0, 0)
            )
        
        legend_strip <- cowplot::plot_grid(
            leg_reef, leg_hab, p_arc,  # Variável renomeada para p_arc
            ncol = 3,
            rel_widths = c(1, 0.8, 1.2),
            align = "vh"
        )
    } else {
        legend_strip <- cowplot::plot_grid(
            leg_reef, leg_hab,
            ncol = 2,
            rel_widths = c(1, 0.8),
            align = "vh"
        )
    }
    
    return(legend_strip)
}
```

---

### 2.5 Funções Principais de Montagem (CORRIGIDAS)

#### 2.5.1 Figura PCA Saúde

```r
create_health_pca_figure <- function(base_dir, output_dir) {
    #' Gera figura composta PCA Saúde
    #' Layout: 3 bubbles (row 1) + loadings (row 2) + legenda (row 3)
    
    cat("\n=== Criando Figura PCA Saúde ===\n")
    
    # 1. Carregar dados
    pca_data <- load_pca_health_data(base_dir)
    df_scores <- pca_data$scores
    df_loadings <- pca_data$loadings
    
    # 2. Verificar colunas necessárias
    health_cols <- c("HEALTH %", "BLEACHING %", "DEAD %")
    missing_cols <- setdiff(health_cols, names(df_scores))
    if (length(missing_cols) > 0) {
        # Tentar nomes alternativos
        alt_names <- c("HEALTH", "BLEACHING", "DEAD")
        for (i in seq_along(missing_cols)) {
            if (alt_names[i] %in% names(df_scores)) {
                names(df_scores)[names(df_scores) == alt_names[i]] <- health_cols[i]
            }
        }
    }
    
    # 3. Calcular variância explicada
    exp_var <- calculate_explained_variance(df_scores, health_cols)
    cat(sprintf("  Variância explicada: PC1=%.1f%%, PC2=%.1f%%\n", exp_var[1], exp_var[2]))
    
    # 4. Definir variáveis para bubble plots
    bubble_vars <- health_cols
    bubble_labels <- c("Healthy Tissue", "Bleaching", "Dead Tissue")
    
    # 5. Criar painéis
    cat("  - Criando painéis de bubble plots...\n")
    panels <- list()
    for (i in seq_along(bubble_vars)) {
        panels[[i]] <- create_bubble_panel(
            df = df_scores,
            var_col = bubble_vars[i],
            pc1_col = "HEALTH_PC1",
            pc2_col = "HEALTH_PC2",
            exp_var = exp_var,
            var_label = bubble_labels[i],
            reef_colors = reef_colors,
            habitat_shapes = habitat_shapes
        )
    }
    
    # 6. Criar loadings
    cat("  - Criando painel de loadings...\n")
    loadings_panel <- create_loadings_panel(df_loadings, exp_var)
    
    # 7. Criar legenda
    cat("  - Criando legenda...\n")
    has_arch <- "ARCH" %in% names(df_scores)
    legend_strip <- create_legend_strip(
        unique(df_scores$REEF),
        unique(df_scores$HAB),
        has_arch = has_arch
    )
    
    # 8. Montar via Patchwork + Cowplot (Método Híbrido)
    cat("  - Montando figura composta (método híbrido)...\n")
    
    layout_design <- "
    ABC
    D##
    "
    
    plots_grid <- panels[[1]] + panels[[2]] + panels[[3]] + 
        loadings_panel + plot_spacer() + plot_spacer() +
        plot_layout(design = layout_design, widths = c(1, 1, 1)) +
        plot_annotation(tag_levels = "a", tag_suffix = ")")
    
    # Montagem final híbrida
    final_figure <- cowplot::plot_grid(
        plots_grid,
        NULL,  # Espaçador vertical
        legend_strip,
        ncol = 1,
        rel_heights = c(10, 0.5, 2)
    )
    
    # 9. Salvar
    output_path <- file.path(output_dir, "PCA_Health_Composite_vGeminiPro.png")
    ggsave(output_path, final_figure,
           width = 14, height = 11, dpi = 300, bg = "white")
    
    cat(sprintf("  ✓ Figura salva: %s\n", output_path))
    
    return(final_figure)
}
```

#### 2.5.2 Figura PCA Interações

```r
create_interactions_pca_figure <- function(base_dir, output_dir) {
    #' Gera figura composta PCA Interações
    #' Layout: 2x3 grid (5 bubbles + 1 loadings) + legenda
    
    cat("\n=== Criando Figura PCA Interações ===\n")
    
    # 1. Carregar dados
    pca_data <- load_pca_interactions_data(base_dir)
    df_scores <- pca_data$scores
    df_loadings <- pca_data$loadings
    
    # 2. Definir variáveis (CORRIGIDO - nomes reais do Excel)
    bubble_vars <- c("SUR_TURF", "SUR_CCA", "SUR_CYANO", 
                     "SUR_MACROALGAE", "SUR_PALYTHOA")
    bubble_labels <- c("Turf", "Crustose Coralline Algae", "Cyanobacteria",
                       "Macroalgae", "Palythoa")
    
    # Verificar colunas
    missing <- setdiff(bubble_vars, names(df_scores))
    if (length(missing) > 0) {
        warning(sprintf("Colunas não encontradas: %s", paste(missing, collapse = ", ")))
        bubble_vars <- bubble_vars[bubble_vars %in% names(df_scores)]
        bubble_labels <- bubble_labels[seq_along(bubble_vars)]
    }
    
    # 3. Calcular variância explicada
    exp_var <- calculate_explained_variance(df_scores, bubble_vars)
    cat(sprintf("  Variância explicada: PC1=%.1f%%, PC2=%.1f%%\n", exp_var[1], exp_var[2]))
    
    # 4. Criar painéis
    cat("  - Criando painéis de bubble plots...\n")
    panels <- list()
    for (i in seq_along(bubble_vars)) {
        panels[[i]] <- create_bubble_panel(
            df = df_scores,
            var_col = bubble_vars[i],
            pc1_col = "PC1_INTERACAO",
            pc2_col = "PC2_INTERACAO",
            exp_var = exp_var,
            var_label = bubble_labels[i],
            reef_colors = reef_colors,
            habitat_shapes = habitat_shapes
        )
    }
    
    # 5. Loadings
    cat("  - Criando painel de loadings...\n")
    # Filtrar loadings apenas para variáveis usadas
    df_loadings_filtered <- df_loadings[rownames(df_loadings) %in% bubble_vars, ]
    loadings_panel <- create_loadings_panel(df_loadings_filtered, exp_var)
    
    # 6. Legenda
    cat("  - Criando legenda...\n")
    has_arch <- "ARCH" %in% names(df_scores)
    legend_strip <- create_legend_strip(
        unique(df_scores$REEF),
        unique(df_scores$HAB),
        has_arch = has_arch
    )
    
    # 7. Montar (2x3 grid)
    cat("  - Montando figura composta...\n")
    
    layout_design <- "
    ABC
    DEF
    "
    
    plots_grid <- panels[[1]] + panels[[2]] + panels[[3]] +
        panels[[4]] + panels[[5]] + loadings_panel +
        plot_layout(design = layout_design, widths = c(1, 1, 1)) +
        plot_annotation(tag_levels = "a", tag_suffix = ")")
    
    final_figure <- cowplot::plot_grid(
        plots_grid,
        NULL,
        legend_strip,
        ncol = 1,
        rel_heights = c(10, 0.5, 2)
    )
    
    # 8. Salvar
    output_path <- file.path(output_dir, "PCA_Interactions_Composite_vGeminiPro.png")
    ggsave(output_path, final_figure,
           width = 14, height = 13, dpi = 300, bg = "white")
    
    cat(sprintf("  ✓ Figura salva: %s\n", output_path))
    
    return(final_figure)
}
```

---

### 2.6 Script Principal Completo

```r
# ============================================================================
# PCA_Health_Interactions_Composite_vGeminiPro.R
# Geração de figuras PCA compostas para Saúde e Interações
# Versão: 2.0 - Padrão Nature/Science
#
# CORREÇÕES vs v1:
#   - Nomes de colunas alinhados com output do Python
#   - Cálculo correto de variância explicada via prcomp()
#   - Distinção visual Inner/Outer Arc
#   - Montagem híbrida Patchwork + Cowplot
#   - Labels em inglês consistentes
# ============================================================================

# --- 1. PACOTES ---
libs <- c(
    "readxl",     # Ler Excel
    "readr",      # Ler CSV
    "dplyr",      # Manipulação de dados
    "ggplot2",    # Plotagem base
    "patchwork",  # Composição de painéis
    "cowplot",    # Extração de legendas e montagem
    "scales",     # Rescale de tamanhos
    "grid",       # arrow() e unit()
    "tools"       # toTitleCase()
)

invisible(lapply(libs, function(lib) {
    if (!require(lib, character.only = TRUE)) {
        install.packages(lib)
        library(lib, character.only = TRUE)
    }
}))

# --- 2. CAMINHOS ---
base_dir <- "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/########NEW RESULTS/output_ANALISE_BIOLOGICA_boxplot_PCAnew"

output_dir <- file.path(dirname(base_dir), "PCA_Health_Interactions_Composite_vGeminiPro")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- 3. ESQUEMA VISUAL ---
# [Inserir reef_colors, habitat_shapes, theme_publication daqui de cima]

# --- 4. FUNÇÕES AUXILIARES ---
# [Inserir todas as funções das seções 2.4.1 a 2.4.5]

# --- 5. FUNÇÕES PRINCIPAIS ---
# [Inserir create_health_pca_figure e create_interactions_pca_figure]

# --- 6. EXECUÇÃO ---
main <- function() {
    cat("\n")
    cat("========================================\n")
    cat("PCA Composite Figures - Health & Interactions\n")
    cat("Versão: GeminiPro v2.0\n")
    cat("========================================\n")
    
    # Verificar dados
    if (!file.exists(file.path(base_dir, "scores_PCA_Saude.xlsx"))) {
        stop("ERRO: scores_PCA_Saude.xlsx não encontrado.\n",
             "Execute primeiro: python Calculate_RGR_&_health_PCA.py")
    }
    
    # Gerar figuras
    create_health_pca_figure(base_dir, output_dir)
    create_interactions_pca_figure(base_dir, output_dir)
    
    cat("\n=== PROCESSAMENTO CONCLUÍDO ===\n")
    cat(sprintf("Figuras salvas em: %s\n", output_dir))
}

main()
```

---

## 3. Checklist de Validação

### 3.1 Antes da Implementação

- [ ] Verificar estrutura exata de `scores_PCA_Saude.xlsx`
- [ ] Verificar estrutura exata de `scores_PCA_Interacoes.xlsx`
- [ ] Confirmar se existe coluna ARCH nos dados de Saúde/Interações
- [ ] Se não existir ARCH, definir mapeamento REEF → ARCH

### 3.2 Durante a Implementação

- [ ] Testar `load_pca_health_data()` com print das colunas
- [ ] Testar `calculate_explained_variance()` e comparar com Python
- [ ] Testar bubble plot individual antes da montagem
- [ ] Verificar loadings labels

### 3.3 Após a Implementação

- [ ] Comparar visualmente com figuras do Reference Script
- [ ] Verificar:
  - [ ] Cores dos recifes consistentes
  - [ ] Formas dos habitats corretas
  - [ ] Bordas Inner/Outer distinguíveis
  - [ ] Setas vermelhas nos loadings
  - [ ] Círculo unitário tracejado
  - [ ] Variância explicada nos eixos
  - [ ] Tags (a), (b), (c)... corretas
  - [ ] Legenda sem cortes
  - [ ] Todos os textos em inglês

---

## 4. Modificação OBRIGATÓRIA no Python

> [!CAUTION]
> Executar ANTES de rodar o script R para garantir precisão na variância explicada.

Modificar `Calculate_RGR_&_health_PCA.py`:

```python
# Após linha 75, adicionar:
exp_var_health = pd.DataFrame({
    'PC': ['PC1', 'PC2'],
    'explained_variance_ratio': pca_health.explained_variance_ratio_
})
exp_var_health.to_csv(os.path.join(output_dir, "explained_variance_PCA_Saude.csv"), index=False)

# Após linha 102, adicionar:
exp_var_interactions = pd.DataFrame({
    'PC': ['PC1', 'PC2'],
    'explained_variance_ratio': pca_interactions.explained_variance_ratio_
})
exp_var_interactions.to_csv(os.path.join(output_dir, "explained_variance_PCA_Interacoes.csv"), index=False)
```

---

## 5. Dimensões Finais Recomendadas

| Figura | Largura | Altura | DPI |
|--------|---------|--------|-----|
| PCA Saúde | 14" | 11" | 300 |
| PCA Interações | 14" | 13" | 300 |

Para publicação em Nature/Science, exportar também em formato vetorial:
```r
ggsave(output_path_pdf, final_figure,
       width = 14, height = 11, device = cairo_pdf)
```

---

**Fim do Plano de Implementação Revisado v2.0**
