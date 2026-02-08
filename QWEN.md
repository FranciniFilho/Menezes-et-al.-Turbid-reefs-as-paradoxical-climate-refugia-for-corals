# QWEN.md - Projeto de Pesquisa Marinha: *Mussismilia hispida* no Banco de Abrolhos

## Visão Geral do Projeto

Este é um projeto de pesquisa científica em ecologia marinha que investiga a distribuição, abundância, estado de saúde e taxas de crescimento do coral endêmico construtor de recifes ***Mussismilia hispida*** no Banco de Abrolhos, Brasil. 

O estudo emprega métodos estatísticos avançados baseados em Bayes, especificamente modelos **Beta Inflado por Zero e Um (ZOIB)** implementados através do pacote **`brms`**, para quantificar as respostas do coral a gradientes ambientais como:
- Temperatura
- Disponibilidade de luz
- Concentração de clorofila
- Profundidade

A análise compara essas respostas entre diferentes habitats e arcos de recifes da região.

## Hipótese Central

Hipotetiza-se que habitats caracterizados pelas seguintes condições atuam como **refúgios ambientais**, promovendo condições ideais para o coral:
- Maior estabilidade térmica
- Temperaturas mais baixas
- Menor incidência de luz
- Maiores concentrações de clorofila

Em contraste, locais com características opostas (maior variabilidade, temperaturas mais quentes, aumento da luz e menor clorofila) são previstos para induzir estresse fisiológico. Assim, espera-se que corais em refúgios exibam taxas de crescimento superiores e melhor saúde (menor prevalência de branquidão) em comparação com corais em ambientes mais estressantes.

## Estrutura de Diretórios

```
#######FINAL_CODES/          # Todos os scripts computacionais (localização atual)
#######FINAL_RESULTS/        # Resultados e saídas dos modelos
#######LAST ROUND MS/        # Manuscrito científico (Main text.docx, References.docx)
#######PHOTOS AND MAPS/      # Visualizações e figuras
```

## Fluxo de Trabalho Completo

### Etapa 1: Sumarização de Dados Ambientais (Python)
```bash
python PCA_LOCAL_bubbleplot.py
```
- Realiza PCA em variáveis ambientais (SST, DLI, CHL, DHW)
- Gera gráficos de bolhas para caracterização ambiental em nível de sítio
- Saída: Pontuações PCA e visualizações em `#####output_local_PCA_CV_*_FINAL/`

### Etapa 2: Análise PERMANOVA (Python)
```bash
python PERMANOVA_Local_PCA.py
```
- Avalia se o fator REEF captura variabilidade ambiental
- Justifica a remoção do REEF como variável categórica dos modelos subsequentes
- Se REEF explica variação ambiental significativa, é excluído para evitar multicolinearidade

### Etapa 3: Cálculos de Saúde e Crescimento (Python)
```bash
python Calculate_RGR_&_health_PCA.py
```
- Calcula Taxa de Crescimento Relativo (RGR) a partir de medições de coral
- Realiza PCA em métricas de saúde para derivar HEALTH_PC1 e HEALTH_PC2
- Saída: `resultados_biologicos_por_colonia.csv`

### Etapa 4: Integração de Dados (Python)
```bash
python FINAL_DATA_INTEGRATION.py
```
- Integra dados biológicos (saúde/crescimento) + ambientais (PCA) + dados de frequência por recife, sítio e tipos de habitat
- Adiciona arcos (ARCH) e coordenadas lat/lon
- Saída: `dados_abundancia_integrados_long_format.csv` em `#####output_local_PCA_CV_*_FINAL/`
- **Três cenários CV**: CV_02, CV_30, CV_ALL

### Etapa 5: Modelagem (R)

#### 5a. Modelos Bayesianos (Análise Primária)
```bash
Rscript Bayes_models/01_ZOIB_Abundance_Full_LOO_Selection.R
```
- Modelos Zero-One Inflated Beta (ZOIB) para cobertura de coral (proporção [0,1])
- 3 conjuntos de dados (CV_02, CV_30, CV_ALL) x 2 cenários de priori (WeaklyInformative, Informative)
- Múltiplos modelos candidatos (nulo, linear, spline, interação ARCH)
- Saída: `C:/Users/rbfra/OneDrive/Bayesian_Analyses_ZOIB_Abundance_MULTI_CV/`

```bash
Rscript Bayes_models/02_Gaussian_Health_RGR_Full_LOO_Selection.R
```
- Modelos gaussianos para HEALTH_PC1, HEALTH_PC2, RGR
- Mesma estrutura multi-conjunto x multi-priori
- Saída: `C:/Users/rbfra/OneDrive/Bayesian_Analyses_Health_Growth_MultiDataset/`

#### 5b. Modelos BRT (Interpretação Complementar)
```bash
Rscript Abundance_models_LM_BRT.R
Rscript Health_&_Growth_models_LM_BRT_NEW.R
```
- Árvores de Regressão Impulsionadas para abundância e saúde/crescimento
- Fornece interpretação complementar aos modelos bayesianos
- Útil para capturar relações não-lineares e importância de variáveis

### Etapa 6: Visualização (R)

#### 6a. Pipeline Mestre de Visualização
```bash
Rscript Bayes_models/05_MASTER_Viz_Pipeline_v5_FINAL.R
```
- Descobre automaticamente modelos vencedores (arquivos com "WINNER" no nome)
- Gera gráficos de floresta de qualidade para publicação (meio-olho) para efeitos fixos
- Inclui variáveis categóricas (níveis de HAB mostrados como coeficientes separados)
- Saída: `C:/Users/rbfra/OneDrive/Bayesian_Figures_Publication/`

#### 6b. Visualizações Específicas de Bayesianos
```bash
Rscript Bayes_Viz_REFACTORED.R
```
- Verificações de predição posterior
- Gráficos de efeitos condicionais
- Gráficos diagnósticos (R-hat, ESS, gráficos de traço)

#### 6c. Visualizações BRT
```bash
Rscript BRT_PDP_PLOTS.R
```
- Gráficos de Dependência Parcial para modelos BRT
- Mostra efeitos marginais de preditores sobre variáveis resposta

## Dependências Principais

### R
```r
# Modelagem Bayesiana Principal
libs <- c("brms", "cmdstanr", "bayesplot", "tidybayes", "ggdist")

# Modelos BRT
libs <- c("gbm", "dismo", "randomForest")

# Visualização
libs <- c("ggplot2", "patchwork", "cowplot")

# Manipulação de dados
libs <- c("dplyr", "tidyverse", "readxl")
```

**Crítico**: Requer backend CmdStanR para brms. Defina os núcleos antes de executar:
```r
options(mc.cores = 4)  # Ajuste com base na CPU disponível
set.seed(42)
```

### Python
```python
# Processamento de dados
pandas, numpy, xarray

# Análise
sklearn (PCA), scipy

# Downloads de sensoriamento remoto
requests, PyPDF2, tqdm
```

## Arquitetura dos Modelos

### Modelos ZOIB (Abundância)
- **Resposta**: `COVER_PROP` (proporção [0,1], inclui zeros e uns)
- **Família**: `zero_one_inflated_beta()`
- **Preditores**: Componentes PCA escalonados (PC1_MAGNITUDE, PC2_MAGNITUDE, PC1_VARIABILITY, PC2_VARIABILITY), DEPTH_M, HAB (tipo de habitat), ARCH (arquitetura)
- **Efeitos Aleatórios**: `(1 | SITE)` para estrutura hierárquica
- **Prioris**: Dois cenários - WeaklyInformative e Informative (específicos por conjunto de dados)

### Modelos Gaussianos (Saúde/Crescimento)
- **Resposta**: HEALTH_PC1, HEALTH_PC2, RGR (contínuos)
- **Família**: `gaussian()`
- Estrutura de preditores semelhante aos modelos ZOIB

### Tipos de Modelos
- `model_null`: Intercepto apenas
- `model_linear_*`: Termos lineares + opcional profundidade/habitat
- `model_spline_*`: Splines GAM (k=5) + opcional profundidade/habitat
- `model_interaction_arch_*`: Splines com interação ARCH

## Padrões de Código

### Padrão de Preparação de Dados (R)
Todos os scripts de modelagem seguem este padrão:
1. Carregar dados com detecção automática de delimitador (`read.csv2` com fallback para `read.csv`)
2. Filtrar por `ORGANISMO == "MUSSISMILIA_HISPIDA"`
3. Converter cobertura para proporção: `COVER_PROP = COBERTURA / 100`
4. Escalonar preditores: `scale()` para todas as variáveis contínuas
5. Verificar integridade dos dados (verificações de intervalo, tratamento de NA)

### Padrão de Ajuste de Modelos
```r
brm(
  formula = model_formula,
  data = prepared_data,
  family = zero_one_inflated_beta(),
  prior = priors,
  backend = "cmdstanr",
  cores = 4,
  iter = 4000,
  warmup = 2000,
  chains = 4,
  control = list(adapt_delta = 0.95)
)
```

### Seleção de Modelos Vencedores
Modelos salvos com prefixo `WINNER_` com base na comparação LOOIC via `loo_compare()`

## Localização de Arquivos Importantes

### Caminhos de Dados de Entrada (Fixos - Atualizar se mover o projeto)
```r
# Saídas PCA para 3 cenários CV
"../#######FINAL_RESULTS/#####output_local_PCA_CV_2_FINAL/dados_abundancia_integrados_long_format.csv"
"../#######FINAL_RESULTS/#####output_local_PCA_CV_30_FINAL/dados_abundancia_integrados_long_format.csv"
"../#######FINAL_RESULTS/#####output_local_PCA_CV_all_FINAL/dados_abundancia_integrados_long_format.csv"
```

### Diretórios de Saída
- **Saídas PCA**: `../#######FINAL_RESULTS/#####output_local_PCA_CV_*_FINAL/`
- **ZOIB Abundância**: `C:/Users/rbfra/OneDrive/Bayesian_Analyses_ZOIB_Abundance_MULTI_CV/`
- **Saúde/Crescimento Bayesianos**: `C:/Users/rbfra/OneDrive/Bayesian_Analyses_Health_Growth_MultiDataset/`
- **BRT Abundância**: `../#######FINAL_RESULTS/BRT_Abundance_CV_*/`
- **BRT Saúde/Crescimento**: `../#######FINAL_RESULTS/BRT_Health_Growth_CV_*/`
- **Figuras**: `C:/Users/rbfra/OneDrive/Bayesian_Figures_Publication/`

## Tarefas Comuns

### Testar Sintaxe de Priori
Use `test_prior_syntax.R` para verificar a construção de priori antes de executar modelos completos.

### Executar Modelo Único
Edite o loop de `candidate_models` para ajustar apenas o modelo de interesse, ou extraia a fórmula e priori específicas para ajuste manual.

### Regenerar Figuras
- **Modelos Bayesianos**: Execute `MASTER_Viz_Pipeline.R` (descobre automaticamente modelos vencedores) e `Bayes_Viz_REFACTORED.R` (diagnósticos e efeitos condicionais)
- **Modelos BRT**: Execute `BRT_PDP_PLOTS.R` para gráficos de dependência parcial

### Adicionar Novas Variáveis Preditoras
1. Adicione a `pca_predictors` ou `predictor_lists` na função de preparação de dados
2. Adicione a versão escalonada no loop de preparação de dados
3. Adicione a especificação de priori em `dataset_specific_priors`
4. Atualize fórmulas de modelo em `candidate_models`

## Variáveis Ambientais (Componentes PCA)

- **PC1_MAGNITUDE/PC2_MAGNITUDE**: Condições ambientais médias (SST, DLI, CHL)
- **PC1_VARIABILITY/PC2_VARIABILITY**: Variabilidade temporal (CV, métricas de frequência)
- Derivadas da análise de periodograma de Lomb-Scargle via `TIME_SERIES_LAG_temporal_environ.py`

## Tipos de Habitat (HAB)
- **RR**: Recifes Rochosos
- **TP**: Topos de Recifes
- **PA**: Paredes de Recifes

## Posicionamento cruzado da prateleira dos sítios (ARCH, deve ser ARC)
- **Arco Interno**: costa adentro
- **Arco Externo**: costa fora

## Notas de Desempenho

- Modelos bayesianos podem levar **horas a dias** dependendo do tamanho dos dados e complexidade do modelo
- Use `mc.cores` para processamento paralelo
- Monitore convergência via estatísticas R-hat (devem ser < 1.01)
- Verifique o tamanho efetivo da amostra (ESS) para estimativas confiáveis
- O cache do brms nos diretórios de saída pode ser limpo para execuções frescas: `unlink(brms_cache_dir, recursive = TRUE)`