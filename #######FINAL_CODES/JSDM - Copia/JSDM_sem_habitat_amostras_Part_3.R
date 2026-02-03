### ------------------------------------------------------------------------------------- ###
### SCRIPT 3 - REFINAMENTO COM PRIORS INFORMATIVOS (OPCIONAL)
### ------------------------------------------------------------------------------------- ###

rm(list = ls())
gc()

cat("\n", rep("=", 70), "\n", sep = "")
cat("🎯 SCRIPT 3: REFINAMENTO COM PRIORS INFORMATIVOS (OPCIONAL)\n")
cat(rep("=", 70), "\n", sep = "")

# Pacotes necessários
libs <- c("brms", "dplyr", "tidyr", "posterior", "bayesplot")
invisible(lapply(libs, require, character.only = TRUE))

output_dir <- file.path(getwd(), "JSDM_Dirichlet_Publication_Ready")
refinement_dir <- file.path(output_dir, "Refinamento_Priors")
dir.create(refinement_dir, showWarnings = FALSE)

### ------------------------------------------------------------------------------------- ###
### 1. CARREGAR MODELO EXISTENTE E ANALISAR ESTRUTURA
### ------------------------------------------------------------------------------------- ###

cat("\n--- CARREGANDO MODELO EXISTENTE ---\n")

if (!file.exists(file.path(output_dir, "modelo_principal_publicacao.rds"))) {
  stop("❌ Execute primeiro o Script 1 para gerar o modelo base!")
}

modelo_base <- readRDS(file.path(output_dir, "modelo_principal_publicacao.rds"))
data_for_brms <- readRDS(file.path(output_dir, "data_for_brms_completo.rds"))
organism_vars <- readRDS(file.path(output_dir, "organism_vars_completo.rds"))

cat("✅ Modelo base carregado:\n")
cat("   - Fórmula:", as.character(modelo_base$formula$formula)[1], "\n")
cat("   - Família:", modelo_base$family$family, "\n")

### ------------------------------------------------------------------------------------- ###
### 2. ANALISAR ESTRUTURA REAL DOS PARÂMETROS
### ------------------------------------------------------------------------------------- ###

cat("\n🔍 ANALISANDO ESTRUTURA DOS PARÂMETROS...\n")

# Extrair parâmetros reais do modelo
posterior_samples <- as_draws_df(modelo_base)
parametros_reais <- names(posterior_samples)

cat("📊 PARÂMETROS REAIS DO MODELO:\n")

# Coeficientes de regressão
coef_params <- grep("^b_", parametros_reais, value = TRUE)
cat("\n--- COEFICIENTES DE REGRESSÃO (b_) ---\n")
print(head(coef_params, 10))
cat("Total de coeficientes:", length(coef_params), "\n")

# Efeitos aleatórios
sd_params <- grep("^sd_", parametros_reais, value = TRUE)
cat("\n--- EFEITOS ALEATÓRIOS (sd_) ---\n")
print(sd_params)

# Correlações
cor_params <- grep("^cor_", parametros_reais, value = TRUE)
cat("\n--- CORRELAÇÕES (cor_) ---\n")
print(head(cor_params, 5))

# Interceptos
intercept_params <- grep("^b_.*Intercept", parametros_reais, value = TRUE)
cat("\n--- INTERCEPTOS ---\n")
print(intercept_params)

# Salvar análise de parâmetros
param_analysis <- data.frame(
  tipo_parametro = c("Coeficientes", "Efeitos Aleatórios", "Correlações", "Interceptos"),
  n_parametros = c(length(coef_params), length(sd_params), length(cor_params), length(intercept_params)),
  exemplos = c(
    paste(head(coef_params, 3), collapse = ", "),
    paste(sd_params, collapse = ", "),
    paste(head(cor_params, 3), collapse = ", "),
    paste(intercept_params, collapse = ", ")
  )
)

write.csv(param_analysis, file.path(refinement_dir, "analise_estrutura_parametros.csv"), row.names = FALSE)

### ------------------------------------------------------------------------------------- ###
### 3. DEFINIR PRIORS INFORMATIVOS BASEADOS NA ESTRUTURA REAL
### ------------------------------------------------------------------------------------- ###

cat("\n🎯 CONFIGURANDO PRIORS INFORMATIVOS...\n")

# Primeiro verificar os priors padrão que foram usados
cat("📋 PRIORS PADRÃO UTILIZADOS:\n")
print(prior_summary(modelo_base))

# --- CORREÇÃO PRINCIPAL AQUI ---
# Para modelos multivariados como o Dirichlet, precisamos ser mais explícitos.
# Usamos `dpar` para especificar a qual parâmetro da distribuição o prior se aplica.
# Usar dpar = "" aplica o prior a TODOS os parâmetros distribucionais (ex: mucca, mucyano, etc.) daquela classe.

priors_informativos <- c(
  # Prior para TODOS os coeficientes de inclinação (slopes) de TODAS as espécies
  prior(normal(0, 0.5), class = "b", dpar = ""),
  
  # Prior para TODOS os interceptos de TODAS as espécies
  prior(student_t(3, -2, 2.5), class = "Intercept", dpar = ""),
  
  # Prior para TODOS os desvios padrão dos efeitos aleatórios
  # A classe 'sd' não precisa de 'dpar' pois se refere à estrutura de grupo
  prior(exponential(1), class = "sd")
)

cat("\n📋 PRIORS INFORMATIVOS CONFIGURADOS (CORRIGIDO):\n")

# Loop de impressão corrigido (da nossa conversa anterior)
if(nrow(priors_informativos) > 0) {
  for(i in 1:nrow(priors_informativos)) {
    prior_str <- sprintf("Prior: %s, Classe: %s", 
                         priors_informativos$prior[i], 
                         priors_informativos$class[i])
    if (!is.na(priors_informativos$dpar[i]) && priors_informativos$dpar[i] != "") {
      prior_str <- paste(prior_str, ", Dpar:", priors_informativos$dpar[i])
    }
    cat(paste0(" - ", prior_str, "\n"))
  }
}

### ------------------------------------------------------------------------------------- ###
### 4. EXECUTAR MODELO REFINADO COM PRIORS INFORMATIVOS
### ------------------------------------------------------------------------------------- ###

cat("\n🔄 EXECUTANDO MODELO REFINADO COM PRIORS INFORMATIVOS...\n")

# Usar mesma fórmula do modelo base
formula_refinada <- modelo_base$formula

# Parâmetros MCMC (podem ser reduzidos para teste)
CHAINS_REF <- 4
ITER_REF <- 4000  # Reduzido para teste rápido
WARMUP_REF <- 2000

cat("⏰ Início:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

modelo_refinado <- brm(
  formula = formula_refinada,
  data = data_for_brms,
  prior = priors_informativos,
  chains = CHAINS_REF,
  iter = ITER_REF,
  warmup = WARMUP_REF,
  control = list(
    adapt_delta = 0.95,  # Ligeiramente reduzido para melhor convergência
    max_treedepth = 12
  ),
  backend = "cmdstanr",
  file = file.path(refinement_dir, "modelo_refinado_priors"),
  silent = 2,
  refresh = 200,
  save_pars = save_pars(all = TRUE)
)

cat("⏰ Término:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

### ------------------------------------------------------------------------------------- ###
### 5. COMPARAÇÃO ENTRE MODELOS
### ------------------------------------------------------------------------------------- ###

cat("\n📊 COMPARANDO MODELO BASE VS REFINADO...\n")

tryCatch({
  # Comparação usando LOO-CV
  loo_base <- loo(modelo_base)
  loo_refinado <- loo(modelo_refinado)
  
  comparacao <- loo_compare(loo_base, loo_refinado)
  
  sink(file.path(refinement_dir, "comparacao_modelos.txt"))
  cat("=== COMPARAÇÃO: MODELO BASE vs MODELO REFINADO ===\n\n")
  cat("DATA:", format(Sys.Date(), "%Y-%m-%d"), "\n\n")
  
  cat("--- COMPARAÇÃO LOO-CV ---\n")
  print(comparacao)
  
  cat("\n--- INTERPRETAÇÃO ---\n")
  if(comparacao["modelo_refinado", "elpd_diff"] > 2) {
    cat("✅ Modelo refinado é SIGNIFICATIVAMENTE melhor\n")
  } else if(comparacao["modelo_refinado", "elpd_diff"] > 0) {
    cat("⚠️  Modelo refinado é ligeiramente melhor\n")
  } else {
    cat("🔁 Modelo base é equivalente ou melhor\n")
  }
  
  cat("\n--- RESUMO DOS PRIORS ---\n")
  cat("Priors utilizados no modelo refinado:\n")
  print(priors_informativos)
  
  sink()
  
  # Salvar modelo refinado
  saveRDS(modelo_refinado, file.path(refinement_dir, "modelo_refinado_priors.rds"))
  
  cat("✅ Comparação de modelos concluída\n")
  cat("📄 Relatório salvo: 'comparacao_modelos.txt'\n")
  
}, error = function(e) {
  cat("❌ Erro na comparação de modelos:", e$message, "\n")
  cat("💡 Salvando modelo refinado mesmo sem comparação...\n")
  saveRDS(modelo_refinado, file.path(refinement_dir, "modelo_refinado_priors.rds"))
})

### ------------------------------------------------------------------------------------- ###
### 6. ANÁLISE DE SENSIBILIDADE AOS PRIORS
### ------------------------------------------------------------------------------------- ###

cat("\n🔍 ANALISANDO SENSIBILIDADE AOS PRIORS...\n")

tryCatch({
  # Comparar distribuições prior-posterior para parâmetros-chave
  posterior_refinado <- as_draws_df(modelo_refinado)
  
  # Selecionar alguns parâmetros para análise
  params_analise <- c(
    grep("^b_pc1_magnitude", names(posterior_refinado), value = TRUE)[1],
    grep("^b_pc2_magnitude", names(posterior_refinado), value = TRUE)[1],
    grep("^sd_SITE", names(posterior_refinado), value = TRUE)[1]
  )
  
  if(length(params_analise) > 0) {
    # Criar gráficos de sensibilidade
    plot_list <- list()
    
    for(param in params_analise) {
      if(param %in% names(posterior_refinado)) {
        p <- ggplot(data.frame(value = posterior_refinado[[param]]), aes(x = value)) +
          geom_density(fill = "steelblue", alpha = 0.7, color = NA) +
          labs(title = paste("Posterior:", param),
               x = "Valor", y = "Densidade") +
          theme_minimal()
        
        # Adicionar prior se for um coeficiente
        if(grepl("^b_", param)) {
          prior_x <- seq(-2, 2, length.out = 100)
          prior_y <- dnorm(prior_x, mean = 0, sd = 0.5)
          p <- p + 
            geom_line(data = data.frame(x = prior_x, y = prior_y),
                      aes(x = x, y = y), color = "red", linewidth = 1, linetype = "dashed") +
            annotate("text", x = -1.5, y = max(prior_y)*0.8, 
                     label = "Prior: normal(0, 0.5)", color = "red", size = 3)
        }
        
        plot_list[[param]] <- p
      }
    }
    
    if(length(plot_list) > 0) {
      sensibilidade_plot <- wrap_plots(plot_list, ncol = 2) +
        plot_annotation(title = "Análise de Sensibilidade: Prior vs Posterior")
      
      ggsave(file.path(refinement_dir, "sensibilidade_priors.png"),
             sensibilidade_plot, width = 12, height = 8, dpi = 300, bg = "white")
      
      cat("✅ Análise de sensibilidade salva\n")
    }
  }
  
}, error = function(e) {
  cat("❌ Erro na análise de sensibilidade:", e$message, "\n")
})

### ------------------------------------------------------------------------------------- ###
### 7. RELATÓRIO FINAL DO REFINAMENTO
### ------------------------------------------------------------------------------------- ###

cat("\n", rep("=", 70), "\n", sep = "")
cat("📋 RELATÓRIO FINAL DO REFINAMENTO\n")
cat(rep("=", 70), "\n", sep = "")

sink(file.path(refinement_dir, "RELATORIO_REFINAMENTO.txt"))
cat("=== RELATÓRIO DE REFINAMENTO COM PRIORS INFORMATIVOS ===\n\n")
cat("DATA:", format(Sys.Date(), "%Y-%m-%d"), "\n\n")

cat("--- ESTRATÉGIA ADOTADA ---\n")
cat("1. Modelo base com priors padrão (Script 1)\n")
cat("2. Análise da estrutura real de parâmetros\n")
cat("3. Definição de priors informativos baseados na análise\n")
cat("4. Comparação entre modelos base e refinado\n\n")

cat("--- RECOMENDAÇÕES ---\n")
cat("• Se modelo refinado for significativamente melhor: USAR PARA PUBLICAÇÃO\n")
cat("• Se diferença for pequena: MANTER MODELO BASE (mais conservador)\n")
cat("• Se modelo base for melhor: REVISAR PRIORS INFORMATIVOS\n\n")

cat("--- PRÓXIMOS PASSOS ---\n")
cat("1. Verificar arquivo 'comparacao_modelos.txt'\n")
cat("2. Analisar gráficos de sensibilidade\n")
cat("3. Decidir qual modelo usar na publicação\n")
cat("4. Atualizar métodos do manuscrito conforme decisão\n")

sink()

cat("\n✅ SCRIPT 3 CONCLUÍDO!\n")
cat("📁 Resultados do refinamento salvos em:", refinement_dir, "\n")
cat("🎯 Verifique 'comparacao_modelos.txt' para decisão final\n")
cat(rep("=", 70), "\n", sep = "")