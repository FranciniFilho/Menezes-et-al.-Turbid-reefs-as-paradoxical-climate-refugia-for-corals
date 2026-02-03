### ===================================================================================== ###
### SCRIPT DE AVALIAÇÃO GERAL DO MODELO JSDM DIRICHLET
### ===================================================================================== ###
#
# OBJETIVO:
# Este script realiza uma avaliação abrangente e multifacetada do desempenho
# de um modelo JSDM Dirichlet ajustado com 'brms'. Ele gera três tipos principais
# de diagnósticos para avaliar o poder explicativo do modelo "como um todo".
#
# AUTOR: Assistente de IA
# DATA: 2024-10-27
#
### ===================================================================================== ###

### ------------------------------------------------------------------------------------- ###
### 0. CONFIGURAÇÃO DO AMBIENTE
### ------------------------------------------------------------------------------------- ###

# Limpar ambiente para garantir um início limpo
rm(list = ls())
gc()

cat("--- Iniciando Script de Avaliação Geral do JSDM ---\n")

# Pacotes necessários para a análise
libs <- c("brms", "dplyr", "tidyr", "ggplot2", "patchwork", "posterior", 
          "loo", "vegan", "viridis", "ggrepel")

# Instalar pacotes faltantes
cat("Verificando pacotes necessários...\n")
inst <- libs[!libs %in% installed.packages()[, "Package"]]
if(length(inst)) install.packages(inst, dependencies = TRUE)

# Carregar pacotes
invisible(lapply(libs, library, character.only = TRUE))

# Configurações gráficas
theme_set(theme_light(base_size = 12) +
            theme(plot.title = element_text(face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(hjust = 0.5),
                  legend.position = "bottom"))

### ------------------------------------------------------------------------------------- ###
### 1. DEFINIR CAMINHOS E CARREGAR DADOS
### ------------------------------------------------------------------------------------- ###

cat("\n--- 1. Carregando Dados e Modelo ---\n")

# Diretório principal de saída (onde os resultados do Script 1 estão)
main_output_dir <- file.path(getwd(), "JSDM_Dirichlet_Publication_Ready")

# Novo diretório para os resultados desta avaliação
evaluation_dir <- file.path(main_output_dir, "Avaliacao_Geral_JSDM")
dir.create(evaluation_dir, showWarnings = FALSE)

# Verificar se os arquivos essenciais existem
required_files <- c("modelo_principal_publicacao.rds", 
                    "data_for_brms_completo.rds", 
                    "organism_vars_completo.rds")

if (!all(file.exists(file.path(main_output_dir, required_files)))) {
  stop("ERRO: Um ou mais arquivos necessários não foram encontrados em '", main_output_dir, "'. Execute os scripts anteriores primeiro.")
}

# Carregar objetos
jsdm_model <- readRDS(file.path(main_output_dir, "modelo_principal_publicacao.rds"))
analysis_data <- readRDS(file.path(main_output_dir, "data_for_brms_completo.rds"))
organism_vars <- readRDS(file.path(main_output_dir, "organism_vars_completo.rds"))

cat("✅ Modelo e dados carregados com sucesso.\n")

### ------------------------------------------------------------------------------------- ###
### 2. MÉTODO 1: R² BAYESIANO MULTIVARIADO
### ------------------------------------------------------------------------------------- ###

cat("\n--- 2. Calculando R² Bayesiano Multivariado ---\n")

tryCatch({
  # Extrair R² para cada variável de resposta (espécie)
  r2_draws <- bayes_R2(jsdm_model, summary = FALSE)
  r2_values <- r2_draws$R2
  
  # Resumo estatístico da distribuição dos R²s
  r2_summary <- summary(r2_values)
  
  # Salvar resumo em texto
  sink(file.path(evaluation_dir, "relatorio_r2_multivariado.txt"))
  cat("=== RESUMO DO R² BAYESIANO MULTIVARIADO ===\n\n")
  cat("Esta métrica quantifica a porcentagem da variância explicada pelos\n")
  cat("preditores ambientais para cada espécie individualmente.\n\n")
  print(r2_summary)
  cat("\n\n--- R² por Espécie (Top 5 e Bottom 5) ---\n")
  r2_per_species <- data.frame(
    organismo = organism_vars,
    R2_medio = colMeans(r2_values)
  ) %>% arrange(desc(R2_medio))
  cat("\nMelhor explicados:\n")
  print(head(r2_per_species, 5))
  cat("\nPior explicados:\n")
  print(tail(r2_per_species, 5))
  sink()
  
  # Visualizar a distribuição dos R²s
  r2_plot <- ggplot(data.frame(R2 = r2_values), aes(x = R2)) +
    geom_histogram(aes(y = after_stat(density)), bins = 15, fill = "skyblue", color = "black", alpha = 0.7) +
    geom_density(color = "darkblue", linewidth = 1) +
    labs(title = "Distribuição do Poder Explicativo Ambiental entre as Espécies",
         subtitle = paste0("R² Médio: ", round(mean(r2_values) * 100, 1), "%"),
         x = "R² Bayesiano",
         y = "Densidade")
  
  ggsave(file.path(evaluation_dir, "distribuicao_r2_multivariado.png"), 
         plot = r2_plot, width = 8, height = 6, dpi = 300)
  
  cat("✅ Análise do R² Multivariado concluída.\n")
  
}, error = function(e) {
  cat("ERRO na análise de R²:", e$message, "\n")
})


### ------------------------------------------------------------------------------------- ###
### 3. MÉTODO 2: VALIDAÇÃO DA ESTRUTURA DA COMUNIDADE
### ------------------------------------------------------------------------------------- ###

cat("\n--- 3. Validando Estrutura da Comunidade Predita ---\n")

# Esta etapa pode consumir bastante memória.
cat("   Calculando predições posteriores (pode levar um momento)...\n")
predicted_draws <- posterior_epred(jsdm_model)
predicted_means <- apply(predicted_draws, c(2, 3), mean) # Média das predições para cada amostra
colnames(predicted_means) <- organism_vars

# --- 3.1 Comparação da Composição Média ---
tryCatch({
  observed_props <- colMeans(analysis_data[, organism_vars])
  predicted_props <- colMeans(predicted_means)
  
  comparison_df <- data.frame(Organismo = organism_vars, Observado = observed_props, Predito = predicted_props) %>%
    pivot_longer(cols = c(Observado, Predito), names_to = "Tipo", values_to = "Proporcao")
  
  composition_plot <- ggplot(comparison_df, aes(x = Tipo, y = Proporcao, fill = Organismo)) +
    geom_bar(stat = "identity", position = "fill", color = "white") +
    scale_fill_viridis_d(name = "Organismo") +
    guides(fill = guide_legend(ncol = 2)) +
    labs(title = "Composição Média da Comunidade",
         subtitle = "Observado vs. Média Predita pelo Modelo",
         x = "", y = "Proporção Relativa")
  
  ggsave(file.path(evaluation_dir, "comparacao_composicao_media.png"),
         plot = composition_plot, width = 10, height = 8, dpi = 300)
  
  cat("✅ Gráfico de composição média salvo.\n")
}, error = function(e) {
  cat("ERRO no gráfico de composição:", e$message, "\n")
})


# --- 3.2 Comparação por Ordenação (PCA) ---
tryCatch({
  pca_observed <- rda(analysis_data[, organism_vars])
  pca_predicted <- rda(predicted_means)
  
  # Extrair scores para ggplot
  scores_obs <- as.data.frame(scores(pca_observed, display = "sites"))
  scores_pred <- as.data.frame(scores(pca_predicted, display = "sites"))
  
  # Limites consistentes para os eixos
  x_limits <- range(c(scores_obs$PC1, scores_pred$PC1))
  y_limits <- range(c(scores_obs$PC2, scores_pred$PC2))
  
  # Plots
  plot_obs <- ggplot(scores_obs, aes(x = PC1, y = PC2)) +
    geom_point(alpha = 0.6, color = "darkgreen") +
    xlim(x_limits) + ylim(y_limits) +
    labs(title = "Comunidade Observada", x = "PC1", y = "PC2")
  
  plot_pred <- ggplot(scores_pred, aes(x = PC1, y = PC2)) +
    geom_point(alpha = 0.6, color = "purple") +
    xlim(x_limits) + ylim(y_limits) +
    labs(title = "Comunidade Predita (Média)", x = "PC1", y = "PC2")
  
  # Combinar plots
  ordination_plot <- (plot_obs | plot_pred) +
    plot_annotation(title = "Comparação da Estrutura da Comunidade via Ordenação PCA",
                    subtitle = "O modelo deve capturar padrões de variação similares aos dados observados")
  
  ggsave(file.path(evaluation_dir, "comparacao_ordenacao_pca.png"),
         plot = ordination_plot, width = 12, height = 7, dpi = 300)
  
  cat("✅ Gráfico de ordenação PCA salvo.\n")
}, error = function(e) {
  cat("ERRO na ordenação PCA:", e$message, "\n")
})


### ------------------------------------------------------------------------------------- ###
### 4. MÉTODO 3: COMPARAÇÃO COM MODELO NULO (LOO)
### ------------------------------------------------------------------------------------- ###
cat("\n--- 4. Comparando com Modelo Nulo via LOO ---\n")
cat("AVISO: Esta etapa ajustará um novo modelo e pode ser demorada.\n")

tryCatch({
  # Definir e ajustar o modelo nulo (apenas interceptos aleatórios)
  null_formula <- update(jsdm_model$formula, new = . ~ 1 + (1|SITE/HAB))
  
  null_model <- brm(
    formula = null_formula,
    data = analysis_data,
    family = dirichlet(),
    chains = 4,
    iter = 4000, # Pode ser reduzido se o tempo for um problema
    warmup = 1000,
    backend = "cmdstanr",
    file = file.path(evaluation_dir, "modelo_nulo_cache"), # Cache para não re-rodar
    file_refit = "on_change"
  )
  
  # Calcular LOO para ambos os modelos
  cat("   Calculando LOO para os modelos (pode levar um momento)...\n")
  loo_full <- loo(jsdm_model, reloo = TRUE)
  loo_null <- loo(null_model, reloo = TRUE)
  
  # Comparar os modelos
  loo_comparison <- loo_compare(loo_full, loo_null)
  
  # Salvar relatório
  sink(file.path(evaluation_dir, "relatorio_comparacao_loo.txt"))
  cat("=== COMPARAÇÃO DE MODELOS (COMPLETO vs. NULO) ===\n\n")
  cat("Esta análise quantifica o quanto os preditores ambientais (efeitos fixos)\n")
  cat("melhoram o poder preditivo do modelo como um todo, em comparação com\n")
  cat("um modelo que apenas considera a variação entre locais (efeitos aleatórios).\n\n")
  cat("Modelos:\n")
  cat(" - jsdm_model: Modelo completo com preditores ambientais\n")
  cat(" - null_model: Modelo nulo, apenas com interceptos aleatórios\n\n")
  cat("--- Resultado da Comparação ---\n")
  print(loo_comparison)
  cat("\n--- Interpretação ---\n")
  cat("A métrica chave é 'elpd_diff'. Um valor positivo para o 'jsdm_model' indica\n")
  cat("que ele tem um poder preditivo fora da amostra melhor que o modelo nulo.\n")
  cat("Por convenção, se o 'elpd_diff' for maior que 2x o seu erro padrão ('se_diff'),\n")
  cat("a diferença é considerada significativa.\n")
  sink()
  
  cat("✅ Comparação LOO concluída.\n")
  
}, error = function(e) {
  cat("ERRO na comparação LOO:", e$message, "\n")
})


### ===================================================================================== ###
cat("\n\n✅ ANÁLISE GERAL CONCLUÍDA!\n")
cat("📁 Verifique a pasta '", evaluation_dir, "' para todos os relatórios e gráficos.\n")
### ===================================================================================== ###