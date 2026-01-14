# Plano de Implementação: Boxplots de Publicação Nature/Science
## Figuras de Cobertura Bêntica e Taxa de Crescimento Relativa (RGR)

**Data:** 2026-01-14
**Status:** Planejamento
**Prioridade:** Alta

---

## 1. Análise das Figuras Atuais

### 1.1 Figura Atual: Boxplot Cobertura (BOXPLOT_cover_REEF_HAB.py)

**Arquivo:** `#######FINAL_CODES/BOXPLOT_cover_REEF_HAB.py`

**Características atuais:**
- FacetGrid com múltiplos organismos em painéis separados
- Boxplot simples sem outliers
- Escala y automática por painel
- Seaborn default theme ("whitegrid", "muted")
- Legenda posicionada externamente

**Problemas identificados:**
1. **Inconsistência visual** com figuras PCA (tema, cores, fontes diferentes)
2. **Falta de informações estatísticas** (sem sobreposição de dados brutos)
3. **Densidade visual subutilizada** (boxplot apenas mostra resumos)
4. **Legenda externa** desperdiça espaço valioso
5. **Sem padrão de cor para reef/habitat** consistente com PCA
6. **Tamanho de fonte não otimizado** para publicação

### 1.2 Figura Atual: Boxplot RGR (Calculate_RGR_&_health_PCA.py)

**Arquivo:** `#######FINAL_CODES/Calculate_RGR_&_health_PCA.py` (linhas 197-238)

**Características atuais:**
- Boxplot único para RGR por REEF e HAB
- Linha horizontal em y=0 (referência de crescimento zero)
- Sem outliers
- Seaborn default theme

**Problemas identificados:**
1. Mesma inconsistência visual com figuras PCA
2. Distribuições assimétricas não são bem representadas por boxplot
3. Não mostra a variabilidade natural dos dados
4. Sem informações sobre tamanho amostral (n)

---

## 2. Padrões Estéticos de Referência (PCA Scripts)

### 2.1 Cores e Paletas (Extraído de PCA_Health_Interactions_Composite_GLM.R e PCA_Environ_vGeminiPro.R)

```r
# Cores dos Recifes (REEF)
reef_colors <- c(
    "ARC" = "#1f77b4",  # Azul
    "ITA" = "#ff7f0e",  # Laranja
    "PAB" = "#2ca02c",  # Verde
    "UCR" = "#d62728",  # Vermelho
    "TIM" = "#9467bd"   # Roxo
)

# Shapes dos Habitats (HAB)
habitat_shapes <- c(
    "PA" = 21,  # Círculo preenchido
    "RR" = 22,  # Quadrado preenchido
    "TP" = 24   # Triângulo preenchido
)

# Distinção Inner/Outer Arc
inner_arc <- "black"   # stroke = 1.2
outer_arc <- "grey60"  # stroke = 0.4
```

### 2.2 Tema Publication-Grade (theme_publication)

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

### 2.3 Padrões de Montagem (Baseado em Diagnostic_PCA_Composite_Solutions.md)

- **Engine:** `cowplot::plot_grid` para alinhamento preciso
- **Align:** `align = "vh"` e `axis = "tblr"` para eliminar gutters
- **Legenda:** Strip horizontal unificado usando cowplot
- **Output:** PNG (300 dpi) + PDF (vector) para flexibilidade

---

## 3. Decisões de Design: Tipo de Gráfico

### 3.1 Avaliação de Opções

| Opção | Vantagens | Desvantagens | Recomendação |
|-------|-----------|--------------|---------------|
| **Boxplot simples** | Simples, familiar | Oculta distribuição, perde informações | ❌ Não usar |
| **Violin plot** | Mostra distribuição completa | Pode ser confuso com poucos dados | ⚠️ Usar com cautela |
| **Boxplot + jitter points** | Mostra resumo + dados brutos | Pode ficar overcrowded | ✅ **RECOMENDADO** |
| **Raincloud plot** | Mais informativo, moderno | Requer pacotes adicionais | ⚠️ Alternativa avançada |
| **Sina plot** | Estatisticamente robusto | Menos familiar para revisores | ⚠️ Alternativa complementar |

### 3.2 Recomendação Final: **Boxplot + Stripplot com Jitter**

**Justificativa:**
1. **Padrão Nature/Science**: Combina estatística (box) com transparência (pontos)
2. **Complementaridade**: Boxplot para tendências centrais, pontos para variabilidade
3. **Familiaridade**: Revisores entendem imediatamente
4. **Flexibilidade**: Funciona bem com diferentes tamanhos amostrais
5. **Consistência**: Pode usar cores/shapes já definidos para PCA

### 3.3 Design Específico por Figura

#### Figura 1: Cobertura Bêntica (Multiple Organismos)
- **Estrutura:** Faceted plot (1 organizmo por painel)
- **Layout:** 3 colunas x n linhas (organismos)
- **Eixo X:** REEF (ordenado alfabeticamente)
- **Eixo Y:** COBERTURA (%)
- **Hue/Color:** HAB (PA, RR, TP)
- **Elementos:**
  - Boxplot (sem outliers, largura 0.6)
  - Stripplot com jitter (width=0.15, alpha=0.5)
  - Tamanho amostral (n) acima de cada box

#### Figura 2: RGR (Growth Rate)
- **Estrutura:** Painel único
- **Eixo X:** REEF (ordenado)
- **Eixo Y:** RGR (taxa de crescimento)
- **Hue/Color:** HAB
- **Elementos adicionais:**
  - Linha horizontal tracejada em y=0 (crescimento zero)
  - Destaque visual para valores positivos vs negativos
  - Notação de significância estatística (se aplicável)

---

## 4. Arquitetura do Script R Proposto

### 4.1 Estrutura do Arquivo

```
#######FINAL_CODES/Boxplot_Publication_Grade.R
├── 1. Setup e Pacotes
├── 2. Definição de Caminhos
├── 3. Definição de Tema e Cores (consistente com PCA)
├── 4. Funções Auxiliares
│   ├── load_coverage_data()
│   ├── load_rgr_data()
│   ├── create_boxplot_stripplot_panel()
│   └── create_legend_strip()
├── 5. Geração da Figura de Cobertura
├── 6. Geração da Figura de RGR
└── 7. Exportação (PNG + PDF)
```

### 4.2 Pacotes Necessários

```r
libs <- c(
    "readxl",           # Leitura de Excel/CSV
    "readr",            # Leitura de CSV
    "dplyr",            # Manipulação de dados
    "ggplot2",          # Plotagem base
    "cowplot",          # Composição de painéis e legendas
    "scales",           # Formatação de escalas
    "ggdist",           # Distribuições avançadas (opcional)
    "ggbeeswarm",       # Jitter melhorado (opcional)
    "patchwork"         # Composição alternativa (se necessário)
)
```

### 4.3 Função Principal para Painéis

```r
create_boxplot_stripplot_panel <- function(df, x_var, y_var, color_var,
                                          title, y_label, reef_colors,
                                          habitat_shapes, show_ref_line = FALSE) {
    # Camadas do ggplot2:
    # 1. Boxplot (geom_boxplot)
    # 2. Stripplot com jitter (geom_point ou geom_beeswarm)
    # 3. Linha de referência (se show_ref_line = TRUE)
    # 4. Anotações de tamanho amostral
    # 5. Theme_publication
    # 6. Escalas manuais de cor e shape
}
```

---

## 5. Especificações Técnicas Detalhadas

### 5.1 Dimensões das Figuras

#### Figura Cobertura (Faceted)
- **Width:** 14-16 polegadas (dependendo do número de organismos)
- **Height:** 10-12 polegadas (proporcional ao layout)
- **DPI:** 300 (PNG)
- **Vector:** PDF (para editor)

#### Figura RGR (Single Panel)
- **Width:** 10-12 polegadas
- **Height:** 6-8 polegadas
- **DPI:** 300 (PNG)
- **Vector:** PDF (para editor)

### 5.2 Paleta de Cores

**Manter consistência EXATA com PCA scripts:**
- Não usar cores default do ggplot2
- Usar `reef_colors` definido acima
- Para HAB como hue em cobertura, usar escala cinza ou alternativa visual

**Opção para HAB em cobertura (sem sobrepor com REEF):**
```r
hab_fill_colors <- c(
    "PA" = "#E69F00",  # Amarelo/ouro
    "RR" = "#56B4E9",  # Azul claro
    "TP" = "#009E73"   # Verde azulado
)
```

### 5.3 Tipografia

- **Font:** Arial ou Helvetica (padrão publicação)
- **Base size:** 14.4 pt (equivalente a 12pt em texto)
- **Axis text:** 12 pt
- **Axis title:** 13 pt bold
- **Legend title:** 12 pt bold
- **Legend text:** 11 pt
- **Plot title:** 14 pt bold, hjust=0.5

### 5.4 Anotações e Metadados

**Tamanho amostral (n):**
```r
# Adicionar acima de cada box
stat_summary(fun = length, geom = "text",
             label = sprintf("n=%i", ..y..),
             vjust = -0.5, size = 3)
```

**Significância estatística (opcional):**
```r
# Usar geom_signif ou anotação manual
annotate("text", x = pos_x, y = pos_y,
         label = "*", size = 8)
```

---

## 6. Implementação em Etapas

### Etapa 1: Setup e Consistência Visual (1-2 horas)
- [ ] Criar script `Boxplot_Publication_Grade.R`
- [ ] Definir `theme_publication` consistente com PCA
- [ ] Definir `reef_colors` e `habitat_shapes`
- [ ] Testar carregamento de dados

### Etapa 2: Função de Painel Base (2-3 horas)
- [ ] Implementar `create_boxplot_stripplot_panel()`
- [ ] Testar boxplot + stripplot combinados
- [ ] Ajustar parâmetros de jitter e transparência
- [ ] Adicionar anotações de tamanho amostral

### Etapa 3: Figura de Cobertura (2-3 horas)
- [ ] Filtrar dados por organismo
- [ ] Criar layout facetado
- [ ] Ajustar escalas y por painel
- [ ] Adicionar legenda unificada
- [ ] Exportar PNG + PDF

### Etapa 4: Figura de RGR (1-2 horas)
- [ ] Implementar painel único
- [ ] Adicionar linha de referência y=0
- [ ] Destacar crescimento positivo vs negativo
- [ ] Adicionar legenda unificada
- [ ] Exportar PNG + PDF

### Etapa 5: Validação e Ajustes Finais (1-2 horas)
- [ ] Comparar com figuras PCA para consistência
- [ ] Verificar resolução e tamanhos de fonte
- [ ] Testar em diferentes contextos (paper, apresentação)
- [ ] Documentar código

**Tempo total estimado:** 7-12 horas

---

## 7. Alternativas e Melhorias Futuras

### 7.1 Alternativas ao Boxplot + Stripplot

**Raincloud Plot (se aprovado):**
```r
library(ggdist)
ggplot(df, aes(x = REEF, y = COBERTURA, fill = HAB)) +
    stat_halfeye(alpha = 0.7) +           # Distribuição
    stat_boxplot(width = 0.1, outlier.alpha = 0) +  # Box
    stat_dots(alpha = 0.5) +               # Pontos
    theme_publication
```

**Sina Plot (alternativa robusta):**
```r
library(ggsignif)
# Mais estatisticamente robusto que boxplot
```

### 7.2 Melhorias de Usabilidade

- **Parâmetros configuráveis:** Permitir fácil ajuste de cores, temas, dimensões
- **Output múltiplo:** Gerar versões para papel (dpi 300) e apresentação (dpi 150)
- **Legendas modulares:** Função separada para gerar legendas consistentes
- **Documentação inline:** Comentários explicando cada parâmetro

### 7.3 Validação Estatística

- **Testes de normalidade:** Verificar se boxplot é apropriado
- **Transformações:** Aplicar log/sqrt se necessário para RGR
- **Notação de significância:** Adicionar testes post-hoc se relevante

---

## 8. Checklist de Validação

Antes de considerar as figuras finalizadas:

### Consistência Visual
- [ ] Cores de REEF idênticas às figuras PCA
- [ ] Shapes de HAB idênticos às figuras PCA
- [ ] Fontes e tamanhos consistentes
- [ ] Espessuras de linha consistentes
- [ ] Grid lines idênticas (grey90, dashed, 0.3)

### Qualidade de Publicação
- [ ] Resolução mínima 300 DPI
- [ ] Versão vector (PDF) disponível
- [ ] Sem elementos cortados nas bordas
- [ ] Legenda legível em tamanho reduzido
- [ ] Títulos autoexplicativos

### Transparência de Dados
- [ ] Todos os pontos de dados visíveis (stripplot)
- [ ] Tamanho amostral indicado
- [ ] Outliers tratados apropriadamente
- [ ] Escala y apropriada para cada variável

### Documentação
- [ ] Código comentado
- [ ] Parâmetros configuráveis documentados
- [ ] Instruções de reprodução incluídas

---

## 9. Comparação: Antes vs Depois

### Antes (Python Seaborn - BOXPLOT_cover_REEF_HAB.py)
```
- Tema: seaborn default ("whitegrid", "muted")
- Cores: Não especificadas (default matplotlib)
- Sem dados brutos sobrepostos
- Legenda externa desperdiça espaço
- Sem consistência com figuras PCA
```

### Depois (R ggplot2 Proposto)
```
- Tema: theme_publication (consistente com PCA)
- Cores: reef_colors exatamente definidas
- Boxplot + stripplot com dados brutos
- Legenda unificada otimizada
- Fontes, dimensões, layout consistentes com PCA
- Export PNG (300 dpi) + PDF (vector)
```

---

## 10. Referências e Inspirations

### Nature/Science Boxplots
- **Nature Methods:** Box + jitter plots com cores consistentes
- **Science:** Minimalist boxplots com anotações de n
- **Current Biology:** Faceted boxplots com legendas compartilhadas

### R Packages Documentation
- **ggplot2:** geom_boxplot, geom_point (jitter)
- **ggdist:** stat_halfeye (raincloud plots)
- **ggbeeswarm:** geom_beeswarm (jitter otimizado)
- **cowplot:** plot_grid, get_legend

### Paletas de Cores
- **Okabe-Ito:** Color-blind friendly palette (usada em PCA)
- **Nature Methods:** Paletas de alta distinção visual

---

## 11. Próximos Passos

1. **Aprovação do plano:** Revisar e aprovar este plano
2. **Implementação:** Criar script `Boxplot_Publication_Grade.R`
3. **Testes:** Gerar figuras e validar contra PCA figures
4. **Iteração:** Ajustar baseado em feedback
5. **Documentação:** Adicionar ao CLAUDE.md se aprovado

---

**Documento preparado por:** Claude Code (Anthropic)
**Versão:** 1.0
**Data:** 2026-01-14
