> **HISTORICAL NOTE (2026-09-12):** this guide documents the v5-era JSDM correction work.
> The active pipeline is `#######FINAL_CODES/Bayes_models/New_Bayes_Models_YEAR_RE/`.
> Scripts referenced below now live under `#######FINAL_CODES/Legacy/`; OneDrive output paths
> have been replaced by `#######FINAL_RESULTS/` subdirectories.

# Guia de Correções JSDM Dirichlet

## Modelagem Bayesiana de Distribuição de Corais (*Mussismilia hispida*) no Banco dos Abrolhos

**Data:** 4 de Fevereiro de 2026  
**Projeto:** Menezes et al. - Distribution and Abundance of *Mussismilia hispida* in Abrolhos  
**Autor:** Assistente de Pesquisa (Kimi Code CLI)

---

## 📋 Índice

1. [Contexto do Problema](#contexto-do-problema)
2. [Problemas Identificados](#problemas-identificados)
3. [Soluções Implementadas](#soluções-implementadas)
4. [Arquivos Criados](#arquivos-criados)
5. [Passo a Passo para Execução](#passo-a-passo-para-execução)
6. [Explicação Técnica Detalhada](#explicação-técnica-detalhada)
7. [Solução de Problemas](#solução-de-problemas)
8. [Referências](#referências)

---

## Contexto do Problema

Este projeto utiliza **Joint Species Distribution Models (JSDM)** com distribuição **Dirichlet** via pacote `brms` no R para modelar a composição de comunidades bentônicas no Banco dos Abrolhos, Brasil.

### Estrutura do Modelo

```r
Family: dirichlet
Link: logit
Formula: cbind(MUSSISMILIA_prop, TURF_prop, CCA_prop, CYANO_prop, MACROALGAE_prop) ~ 
    PC1MAGNITUDE + PC2MAGNITUDE + PC1VARIABILITY + PC2VARIABILITY +
    PC1VARIABILITY:ARCH + PC2VARIABILITY:ARCH + HABMERGED + DEPTHM + (1 | REEF)
```

---

## Problemas Identificados

### ❌ Problema 1: Ausência de *Mussismilia hispida* nas Figuras JSDM

**Situação:**
- *M. hispida* é a **categoria de referência** no modelo Dirichlet
- Não aparece nos Forest Plots porque não tem coeficientes estimados
- Os coeficientes são nomeados como: `muTURFprop_PC1MAGNITUDE`, `muCCAprop_PC2MAGNITUDE`, etc.
- **NÃO existe** `muMUSSISMILIAprop_*` porque é o baseline

**Resultado Esperado:**
- *M. hispida* deve ser incluída nas figuras com efeito = 0 (log-ratio relativo a si mesmo)

---

### ❌ Problema 2: PDP/Marginal Effects com Linha Azul Única

**Situação:**
- `conditional_effects()` com `categorical=TRUE` combina todas as espécies em uma única curva
- Tentativa de usar `categorical=FALSE` retorna erro: *"Please set 'categorical' to TRUE"*
- Resultado: uma curva azul única representando a "resposta média da comunidade"

**Resultado Esperado:**
- Curvas **SEPARADAS** para cada uma das 5 espécies taxonômicas

---

## Soluções Implementadas

### ✅ Solução 1: Incluir *M. hispida* no Forest Plot

**Abordagem:** Adicionar manualmente entradas dummy para *M. hispida* com:
- `Estimate = 0` (log-ratio relativo a si mesmo)
- `l95 = 0, u95 = 0` (sem incerteza, pois é a referência)

**Código:**
```r
# Identificar todos os preditores únicos
all_predictors <- unique(plot_data$predictor)

# Criar entradas para Mussismilia (referência = 0)
mussismilia_rows <- data.frame(
    Parameter = paste0("muMUSSISMILIAprop_", all_predictors),
    Estimate = 0,      # Referência: log-ratio = 0
    l95 = 0,           # Sem incerteza (baseline)
    u95 = 0,
    species = "MUSSISMILIA",
    predictor = all_predictors,
    stringsAsFactors = FALSE
)

# Combinar dados originais com M. hispida
plot_data <- bind_rows(plot_data, mussismilia_rows)
```

**Visual:** *M. hispida* aparece como um **diamante** na posição 0, destacada das outras espécies.

---

### ✅ Solução 2: PDPs com Curvas Separadas

**Abordagem:** Usar `posterior_epred()` diretamente em vez de `conditional_effects()`

**Por que funciona:**
- `posterior_epred()` retorna array `[ndraws, n_observations, n_species]`
- Cada "slice" do array representa uma espécie
- Permite plotar cada espécie como linha separada

**Código:**
```r
# Criar grade de predição
newdata <- data.frame(
    PC1MAGNITUDE = seq(min, max, length.out = 100),
    PC2MAGNITUDE = mean(data$PC2MAGNITUDE),  # Manter no valor médio
    ARCH = "inner",                          # Valor de referência
    ...
)

# Obter predições para TODAS as espécies
pred_array <- posterior_epred(model, newdata = newdata, ndraws = 500)
# Dimensões: [500 draws, 100 observations, 5 species]

# Extrair cada espécie
for (sp_idx in 1:n_species) {
    sp_draws <- pred_array[, , sp_idx]  # Todas as amostras para esta espécie
    
    # Calcular estatísticas
    sp_summary <- data.frame(
        x = pred_seq,
        estimate = apply(sp_draws, 2, median),
        lower = apply(sp_draws, 2, quantile, 0.025),
        upper = apply(sp_draws, 2, quantile, 0.975)
    )
}
```

**Visual:** **5 curvas coloridas**, cada uma representando uma espécie, com áreas sombreadas para o 95% CI.

---

## Arquivos Criados

### Localização
```
#######FINAL_CODES/
├── JSDM_Fix_Mhispida_and_PDPs.R          # Funções completas e documentadas
├── JSDM_Correction_Module.R               # Módulo otimizado para integração
├── MASTER_Viz_Pipeline_v5_JSDM_FIXED.R   # Pipeline completo corrigido
├── JSDM_Demo_and_Validation.R             # Demonstração com dados simulados
└── README_JSDM_Corrections.md             # Documentação técnica
```

### Descrição dos Arquivos

| Arquivo | Propósito | Quando Usar |
|---------|-----------|-------------|
| `JSDM_Fix_Mhispida_and_PDPs.R` | Funções completas com múltiplas opções de visualização | Para explorar diferentes tipos de gráficos |
| `JSDM_Correction_Module.R` | Versão otimizada das funções corrigidas | Para integrar no pipeline existente |
| `MASTER_Viz_Pipeline_v5_JSDM_FIXED.R` | Pipeline V5 completo com correções integradas | **Recomendado** - Executar diretamente |
| `JSDM_Demo_and_Validation.R` | Demonstração com dados simulados | Para entender/testar as correções |

---

## Passo a Passo para Execução

### Opção 1: Pipeline Completo (RECOMENDADO) ⭐

**Executar diretamente o pipeline corrigido:**

```r
# Limpar ambiente
rm(list = ls())
gc()

# Executar pipeline completo
source("#######FINAL_CODES/Legacy/MASTER_Viz_Pipeline_v5_JSDM_FIXED.R")
```

**O que acontece:**
1. O script busca automaticamente todos os modelos WINNER nos diretórios configurados
2. Processa cada modelo e gera figuras corrigidas
3. Salva na pasta: `#######FINAL_RESULTS/Bayesian_Figures_Publication_v5_FIXED`

**Figuras geradas para modelos JSDM:**
- `FIGURE_1_Forest_Plot.png` - Forest plot padrão
- `FIGURE_3_Validation_PPC.png` - Posterior predictive checks
- `FIGURE_4_Diagnostics.png` - Diagnósticos MCMC
- `FIGURE_5_JSDM_Forest_with_Reference.png` - **Forest com M. hispida**
- `FIGURE_6_JSDM_PDP_Separate_Curves.png` - **PDPs com curvas separadas**
- `FIGURE_7_JSDM_PPC_by_Species.png` - PPC por espécie

---

### Opção 2: Funções Específicas

**Para processar um modelo específico:**

```r
# 1. Carregar o módulo de correções
source("#######FINAL_CODES/JSDM_Correction_Module.R")

# 2. Carregar o modelo
model <- readRDS("#######FINAL_RESULTS/Bayesian_Full_LOO_Selection_v5/JSDM_Dirichlet/WINNER_JSDM_full_CV_ALL.rds")

# 3. Gerar Forest Plot com M. hispida
p_forest <- generate_jsdm_forest_plot_errorbar_CORRECTED(model, "JSDM Model")
ggsave("forest_com_mhispida.png", p_forest, width = 14, height = 10, dpi = 300)

# 4. Gerar PDPs com curvas separadas
p_pdps <- generate_jsdm_pdp_community_overlaid_CORRECTED(model, "JSDM Model")
ggsave("pdps_curvas_separadas.png", p_pdps, width = 16, height = 12, dpi = 300)

# 5. Gerar PDPs em painéis separados
p_per_species <- generate_jsdm_pdp_per_species_CORRECTED(model, "JSDM Model")
ggsave("pdps_por_especie.png", p_per_species, width = 14, height = 10, dpi = 300)
```

---

### Opção 3: Demonstração com Dados Simulados

**Para entender/testar as correções:**

```r
# Executar demonstração
source("#######FINAL_CODES/Legacy/JSDM_Demo_and_Validation.R")
```

**O que acontece:**
1. Simula dados de comunidade bentônica com 5 espécies
2. Gera 3 figuras de demonstração:
   - `DEMO_Forest_Plot_with_Mhispida.png`
   - `DEMO_PDP_Separate_Curves.png`
   - `DEMO_Comparison_Before_After.png`

---

### Opção 4: Integração no Pipeline Existente

**Para integrar no `MASTER_Viz_Pipeline_v5_FINAL.R`:**

1. **Abrir** o arquivo `MASTER_Viz_Pipeline_v5_FINAL.R`

2. **Localizar** a função `generate_jsdm_forest_plot_errorbar()` (aprox. linhas 567-650)

3. **Substituir** pelo conteúdo de `generate_jsdm_forest_plot_errorbar_CORRECTED()` em `JSDM_Correction_Module.R`

4. **Localizar** a função `generate_jsdm_pdp_community_overlaid()` (aprox. linhas 653-737)

5. **Substituir** pelo conteúdo de `generate_jsdm_pdp_community_overlaid_CORRECTED()`

6. **Salvar** e executar o pipeline normalmente

---

## Explicação Técnica Detalhada

### Por que *M. hispida* não aparecia?

#### Modelo Dirichlet em `brms`

O modelo Dirichlet modela composições (proporções que somam 1) usando uma parametrização multivariada:

```
Y_i ~ Dirichlet(α_i1, α_i2, ..., α_iK)
```

Onde:
- `Y_i` = vetor de proporções para observação `i`
- `α_ik` = parâmetro de concentração para categoria `k`

#### Link Logit e Categoria de Referência

O `brms` usa link logit para modelar os log-ratios:

```
log(α_ik / α_i1) = η_ik   para k = 2, ..., K
```

Onde:
- Categoria 1 (*M. hispida*) é a **referência**
- `η_ik` é o preditor linear para categoria `k` relativo à referência
- Por definição: `log(α_i1 / α_i1) = log(1) = 0`

#### Implicação

- Não existe coeficiente para *M. hispida* (é o baseline)
- Todos os outros coeficientes são **log-ratios relativos a M. hispida**
- Portanto, *M. hispida* deve aparecer no plot com effect = 0

---

### Por que `conditional_effects()` não funcionava?

#### Comportamento de `conditional_effects()`

Para modelos categóricos/multinomiais:

```r
conditional_effects(model, categorical = TRUE)
```

- Combina todas as categorias em uma visualização única
- Mostra probabilidades condicionais médias
- Não permite separar por espécie

#### Alternativa: `posterior_epred()`

```r
posterior_epred(model, newdata = newdata)
```

- Retorna predições para **todas** as respostas simultaneamente
- Array 3D: `[draws, observations, responses]`
- Permite extrair e plotar cada resposta separadamente

---

## Solução de Problemas

### Erro: "resp argument not supported"

**Causa:** Tentando usar `resp=` em `posterior_predict()` para modelo Dirichlet

**Solução:** Usar `posterior_epred()` sem argumento `resp` e extrair dimensões do array

```r
# ERRADO
pred <- posterior_predict(model, resp = "TURF_prop")

# CERTO
pred_array <- posterior_epred(model, newdata = newdata)
turf_pred <- pred_array[, , 2]  # Extrair segunda espécie
```

---

### Erro: "dimensions do not match"

**Causa:** `newdata` não contém todas as colunas necessárias

**Solução:** Verificar que todas as variáveis do modelo estão em `newdata`

```r
# Verificar variáveis do modelo
all.vars(formula(model))

# Criar newdata completo
newdata <- data.frame(
    PC1MAGNITUDE = pred_seq,
    PC2MAGNITUDE = mean(data$PC2MAGNITUDE),
    PC1VARIABILITY = mean(data$PC1VARIABILITY),
    PC2VARIABILITY = mean(data$PC2VARIABILITY),
    DEPTHM = mean(data$DEPTHM),
    ARCH = "inner",
    HABMERGED = "Plateau"
)
```

---

### M. hispida não aparece no plot

**Causa:** Esqueceu de adicionar as linhas dummy

**Solução:** Verificar que `bind_rows()` foi chamado corretamente

```r
# Verificar
plot_data <- bind_rows(plot_data, mussismilia_rows)
nrow(plot_data)  # Deve aumentar
"MUSSISMILIA" %in% plot_data$species  # Deve ser TRUE
```

---

## Referências

### Documentação
- `brms`: https://paul-buerkner.github.io/brms/reference/brms-package.html
- Dirichlet regression: https://en.wikipedia.org/wiki/Dirichlet_distribution
- Compositional data analysis: Aitchison, J. (1986). *The Statistical Analysis of Compositional Data*

### Pacotes R
```r
library(brms)      # Modelagem Bayesiana
library(ggplot2)   # Visualização
library(dplyr)     # Manipulação de dados
library(tidyr)     # Reformatação
library(patchwork) # Combinação de plots
library(tidybayes) # Visualização Bayesiana
```

---

## Contato e Suporte

Para dúvidas ou problemas:
1. Verificar este guia na seção de Solução de Problemas
2. Consultar a documentação técnica em `#######FINAL_CODES/README_JSDM_Corrections.md`
3. Executar o script de demonstração para validar as correções

---

**Última atualização:** 4 de Fevereiro de 2026  
**Versão:** 1.0  
**Status:** ✅ Aprovado para uso em produção
