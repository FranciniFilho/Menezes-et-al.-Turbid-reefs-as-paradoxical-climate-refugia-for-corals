### ------------------------------------------------------------------------------------- ###
### ANÁLISES PÓS-MODELO - VERSÃO PUBLICAÇÃO
### ------------------------------------------------------------------------------------- ###

rm(list = ls())
gc()

### ------------------------------------------------------------------------------------- ###
### 0. CONFIGURAÇÃO - PACOTES PARA PUBLICAÇÃO
### ------------------------------------------------------------------------------------- ###

cat("\n", rep("=", 70), "\n", sep = "")
cat("📊 SCRIPT 2: ANÁLISES PROFISSIONAIS PARA PUBLICAÇÃO\n")
cat(rep("=", 70), "\n", sep = "")

# Pacotes para análise profissional
libs <- c("brms", "dplyr", "tidyr", "ggplot2", "patchwork", "bayesplot", 
          "ggpubr", "posterior", "loo", "kableExtra", "gridExtra", "viridis")
invisible(lapply(libs, require, character.only = TRUE))

# Configurações gráficas profissionais
theme_set(theme_minimal(base_size = 11) +
            theme(plot.title = element_text(face = "bold", size = 12),
                  axis.title = element_text(face = "bold"),
                  legend.position = "bottom"))

output_dir <- file.path(getwd(), "JSDM_Dirichlet_Publication_Ready")

### ------------------------------------------------------------------------------------- ###
### 1. CARREGAMENTO E VERIFICAÇÃO - ROBUSTO
### ------------------------------------------------------------------------------------- ###

cat("\n--- CARREGAMENTO DE DADOS ---\n")

required_files <- c("modelo_principal_publicacao.rds", 
                    "data_for_brms_completo.rds", 
                    "organism_vars_completo.rds")

missing_files <- required_files[!file.exists(file.path(output_dir, required_files))]
if(length(missing_files) > 0) {
  stop("❌ Arquivos faltantes: ", paste(missing_files, collapse = ", "),
       "\nExecute o Script 1 primeiro!")
}

# Carregar dados
jsdm_brms_fit_final <- readRDS(file.path(output_dir, "modelo_principal_publicacao.rds"))
data_for_brms <- readRDS(file.path(output_dir, "data_for_brms_completo.rds"))
organism_vars <- readRDS(file.path(output_dir, "organism_vars_completo.rds"))
metadata <- readRDS(file.path(output_dir, "metadata_publicacao.rds"))

cat("✅ Dados carregados:\n")
cat("   - Modelo:", class(jsdm_brms_fit_final), "\n")
cat("   - Observações:", nrow(data_for_brms), "\n")
cat("   - Organismos:", length(organism_vars), "\n")

### ------------------------------------------------------------------------------------- ###
### 2. DIAGNÓSTICOS DE CONVERGÊNCIA COMPLETOS
### ------------------------------------------------------------------------------------- ###

cat("\n", rep("=", 70), "\n", sep = "")
cat("🔍 DIAGNÓSTICOS DE CONVERGÊNCIA - ANÁLISE PROFISSIONAL\n")
cat(rep("=", 70), "\n", sep = "")

diagnostics_dir <- file.path(output_dir, "Diagnosticos_Convergencia")
dir.create(diagnostics_dir, showWarnings = FALSE)

# Sumário completo do modelo
model_summary <- summary(jsdm_brms_fit_final)

# Relatório de diagnóstico profissional
sink(file.path(diagnostics_dir, "relatorio_convergencia_completo.txt"))
cat("=== RELATÓRIO DE CONVERGÊNCIA - MODELO JSDM DIRICHLET ===\n\n")
cat("Data da análise:", format(Sys.Date(), "%Y-%m-%d"), "\n")
cat("Fórmula:", metadata$model_formula, "\n\n")

cat("--- ESTATÍSTICAS DE CONVERGÊNCIA ---\n")
cat("EFEITOS FIXOS:\n")
print(model_summary$fixed)

if(!is.null(model_summary$random)) {
  cat("\nEFEITOS ALEATÓRIOS:\n")
  print(model_summary$random)
}

# Análise de Rhat
rhat_fixed <- model_summary$fixed$Rhat
cat(sprintf("\n--- ANÁLISE RHAT ---\n"))
cat(sprintf("Mínimo: %.4f\n", min(rhat_fixed)))
cat(sprintf("Médio:  %.4f\n", mean(rhat_fixed)))
cat(sprintf("Máximo: %.4f\n", max(rhat_fixed)))
cat(sprintf("> 1.01: %d/%d (%.1f%%)\n", 
            sum(rhat_fixed > 1.01), length(rhat_fixed),
            mean(rhat_fixed > 1.01) * 100))
cat(sprintf("> 1.05: %d/%d (%.1f%%)\n",
            sum(rhat_fixed > 1.05), length(rhat_fixed),
            mean(rhat_fixed > 1.05) * 100))

# Tamanhos amostrais efetivos
cat(sprintf("\n--- TAMANHOS AMOSTRAIS EFETIVOS (ESS) ---\n"))
cat(sprintf("Bulk ESS mínimo: %.0f\n", min(model_summary$fixed$Bulk_ESS)))
cat(sprintf("Tail ESS mínimo: %.0f\n", min(model_summary$fixed$Tail_ESS)))
cat(sprintf("Bulk ESS < 400: %d/%d\n", 
            sum(model_summary$fixed$Bulk_ESS < 400), length(rhat_fixed)))
cat(sprintf("Tail ESS < 400: %d/%d\n",
            sum(model_summary$fixed$Tail_ESS < 400), length(rhat_fixed)))

sink()

### ------------------------------------------------------------------------------------- ###
### 3. VISUALIZAÇÕES DE DIAGNÓSTICO PROFISSIONAIS
### ------------------------------------------------------------------------------------- ###

cat("\n📈 GERANDO VISUALIZAÇÕES DE DIAGNÓSTICO...\n")

tryCatch({
  # 1. Trace plots para parâmetros principais
  posterior_samples <- as_draws_df(jsdm_brms_fit_final)
  
  # Selecionar parâmetros representativos
  param_patterns <- c("b_", "sd_SITE", "cor_SITE")
  selected_params <- names(posterior_samples)[
    sapply(param_patterns, function(p) grepl(p, names(posterior_samples)))
  ]
  
  if(length(selected_params) > 0) {
    # Limitar a 12 parâmetros para clareza visual
    if(length(selected_params) > 12) {
      selected_params <- selected_params[1:12]
    }
    
    trace_data <- posterior_samples[, selected_params] %>%
      mutate(iteration = 1:n()) %>%
      pivot_longer(cols = -iteration, names_to = "parameter", values_to = "value")
    
    p_trace <- ggplot(trace_data, aes(x = iteration, y = value, color = parameter)) +
      geom_line(alpha = 0.7, linewidth = 0.3) +
      facet_wrap(~ parameter, scales = "free_y", ncol = 3) +
      labs(title = "Trace Plots - Parâmetros Principais",
           subtitle = "Diagnóstico de Convergência") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(diagnostics_dir, "trace_plots_profissionais.png"), 
           p_trace, width = 14, height = 10, dpi = 300, bg = "white")
  }
  
  # 2. Distribuição de Rhat
  rhat_data <- data.frame(
    parameter = rownames(model_summary$fixed),
    rhat = model_summary$fixed$Rhat,
    type = ifelse(grepl("Intercept", rownames(model_summary$fixed)), "Intercept", "Slope")
  )
  
  p_rhat <- ggplot(rhat_data, aes(x = rhat, fill = type)) +
    geom_histogram(bins = 30, alpha = 0.7) +
    geom_vline(xintercept = 1.05, linetype = "dashed", color = "red", linewidth = 1) +
    geom_vline(xintercept = 1.01, linetype = "dotted", color = "orange", linewidth = 1) +
    labs(title = "Distribuição de Rhat - Diagnóstico de Convergência",
         subtitle = "Linhas: 1.01 (laranja) e 1.05 (vermelho)",
         x = "Rhat", y = "Frequência") +
    scale_fill_viridis_d() +
    theme_minimal()
  
  ggsave(file.path(diagnostics_dir, "rhat_distribution_profissional.png"), 
         p_rhat, width = 10, height = 6, dpi = 300, bg = "white")
  
  cat("✅ Visualizações de diagnóstico salvas\n")
  
}, error = function(e) {
  cat("❌ Erro nas visualizações de diagnóstico:", e$message, "\n")
})

### ------------------------------------------------------------------------------------- ###
### 4. ANÁLISE DETALHADA MUSSISMILIA HISPIDA - PROFISSIONAL
### ------------------------------------------------------------------------------------- ###

cat("\n", rep("=", 70), "\n", sep = "")
cat("🎯 ANÁLISE DETALHADA: MUSSISMILIA HISPIDA\n")
cat(rep("=", 70), "\n", sep = "")

mh_analysis_dir <- file.path(output_dir, "Analise_Mussismilia_Hispida")
dir.create(mh_analysis_dir, showWarnings = FALSE)

if("mussismilia_hispida" %in% organism_vars) {
  
  mh_position <- which(organism_vars == "mussismilia_hispida")
  
  # Métricas de desempenho robustas
  tryCatch({
    predictions <- posterior_epred(jsdm_brms_fit_final)
    
    if(length(dim(predictions)) == 3) {
      predicted_means <- apply(predictions[,,mh_position], 2, mean)
      predicted_intervals <- apply(predictions[,,mh_position], 2, quantile, 
                                   probs = c(0.025, 0.975))
    } else {
      predicted_means <- colMeans(predictions[, mh_position])
      predicted_intervals <- apply(predictions[, mh_position], 2, quantile,
                                   probs = c(0.025, 0.975))
    }
    
    observed_mh <- data_for_brms$mussismilia_hispida
    
    # Múltiplas métricas de avaliação
    performance_metrics <- list(
      r_squared = cor(predicted_means, observed_mh, use = "complete.obs")^2,
      mae = mean(abs(predicted_means - observed_mh), na.rm = TRUE),
      rmse = sqrt(mean((predicted_means - observed_mh)^2, na.rm = TRUE)),
      bias = mean(predicted_means - observed_mh, na.rm = TRUE),
      correlation = cor(predicted_means, observed_mh, use = "complete.obs")
    )
    
    # Gráfico de desempenho preditivo
    performance_df <- data.frame(
      Observado = observed_mh,
      Predito = predicted_means,
      LI = predicted_intervals[1,],
      LS = predicted_intervals[2,]
    )
    
    p_performance <- ggplot(performance_df, aes(x = Observado, y = Predito)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +
      geom_point(alpha = 0.6, color = "steelblue") +
      geom_errorbar(aes(ymin = LI, ymax = LS), alpha = 0.3, width = 0) +
      geom_smooth(method = "lm", se = FALSE, color = "darkblue") +
      labs(title = "Desempenho Preditivo: Mussismilia hispida",
           subtitle = sprintf("R² = %.1f%%, RMSE = %.4f", 
                              performance_metrics$r_squared * 100,
                              performance_metrics$rmse),
           x = "Abundância Observada", 
           y = "Abundância Predita") +
      theme_minimal()
    
    ggsave(file.path(mh_analysis_dir, "desempenho_preditivo_mushis.png"),
           p_performance, width = 10, height = 8, dpi = 300, bg = "white")
    
    # Relatório detalhado
    sink(file.path(mh_analysis_dir, "relatorio_detalhado_mussismilia.txt"))
    cat("=== ANÁLISE DETALHADA: MUSSISMILIA HISPIDA ===\n\n")
    cat("DATA:", format(Sys.Date(), "%Y-%m-%d"), "\n\n")
    
    cat("--- MÉTRICAS DE DESEMPENHO ---\n")
    cat(sprintf("R² (Poder Preditivo): %.1f%%\n", performance_metrics$r_squared * 100))
    cat(sprintf("Correlação: %.3f\n", performance_metrics$correlation))
    cat(sprintf("Erro Absoluto Médio (MAE): %.4f\n", performance_metrics$mae))
    cat(sprintf("Raiz do Erro Quadrático Médio (RMSE): %.4f\n", performance_metrics$rmse))
    cat(sprintf("Viés Médio: %.4f\n", performance_metrics$bias))
    
    cat("\n--- INTERPRETAÇÃO ---\n")
    if(performance_metrics$r_squared < 0.1) {
      cat("• Poder preditivo BAIXO: variáveis ambientais explicam pouca variação\n")
      cat("• Efeitos aleatórios (SITE/HAB) podem dominar a variação\n")
      cat("• Considerar outros preditores ou processos ecológicos\n")
    } else if(performance_metrics$r_squared < 0.3) {
      cat("• Poder preditivo MODERADO\n")
      cat("• Variáveis ambientais têm influência detectável\n")
      cat("• Modelo captura padrões gerais mas não toda variação\n")
    } else {
      cat("• Poder preditivo ALTO\n")
      cat("• Variáveis ambientais são bons preditores da abundância\n")
    }
    
    sink()
    
    cat("✅ Análise detalhada de M. hispida concluída\n")
    
  }, error = function(e) {
    cat("❌ Erro na análise de M. hispida:", e$message, "\n")
  })
}

### ------------------------------------------------------------------------------------- ###
### 5. VISUALIZAÇÕES PROFISSIONAIS DOS EFEITOS PCA
### ------------------------------------------------------------------------------------- ###

cat("\n", rep("=", 70), "\n", sep = "")
cat("📊 VISUALIZAÇÕES PROFISSIONAIS DOS EFEITOS AMBIENTAIS\n")
cat(rep("=", 70), "\n", sep = "")

visualizations_dir <- file.path(output_dir, "Visualizacoes_Profissionais")
dir.create(visualizations_dir, showWarnings = FALSE)

tryCatch({
  pca_predictors <- c("pc1_magnitude", "pc2_magnitude", "pc1_variability", "pc2_variability")
  cat("Gerando efeitos condicionais profissionais...\n")
  professional_plots <- list()
  
  for(pred in pca_predictors) {
    cat("  Processando:", pred, "\n")
    tryCatch({
      cond_effect <- conditional_effects(jsdm_brms_fit_final, 
                                         effects = pred, 
                                         categorical = TRUE,
                                         points = FALSE,
                                         plot = FALSE)
      
      effect_df <- cond_effect[[1]]
      
      # --- MUDANÇA 1: ATIVAR A LEGENDA E DAR UM TÍTULO A ELA ---
      p <- ggplot(effect_df, aes(x = .data[[pred]], y = estimate__, color = cats__)) +
        geom_ribbon(aes(ymin = lower__, ymax = upper__, fill = cats__), 
                    alpha = 0.2, linetype = 0) +
        geom_line(linewidth = 0.8) +
        labs(x = toupper(pred), 
             y = "Proporção Prevista",
             title = paste("Resposta da Comunidade a", toupper(pred))) +
        scale_color_viridis_d(name = "Espécie / Grupo") +
        scale_fill_viridis_d(name = "Espécie / Grupo") +
        theme_minimal() +
        theme(plot.title = element_text(face = "bold", size = 10))
      
      professional_plots[[pred]] <- p
      
    }, error = function(e) {
      cat("    ❌ Erro em", pred, ":", e$message, "\n")
    })
  }
  
  if(length(professional_plots) > 0) {
    # --- MUDANÇA 2: USAR plot_layout() PARA CRIAR UMA LEGENDA ÚNICA ---
    combined_plot <- wrap_plots(professional_plots, ncol = 2) +
      plot_layout(guides = 'collect') + 
      plot_annotation(
        title = "Respostas da Comunidade Bentônica aos Gradientes Ambientais PCA",
        subtitle = "Modelo JSDM Dirichlet com Efeitos Aleatórios Aninhados",
        theme = theme(plot.title = element_text(face = "bold", size = 14),
                      plot.subtitle = element_text(size = 11),
                      legend.position = 'bottom') # Posição da legenda compartilhada
      )
    
    ggsave(file.path(visualizations_dir, "efeitos_pca_profissionais_com_legenda.png"), 
           combined_plot, width = 16, height = 12, dpi = 300, bg = "white")
    
    cat("✅ Visualização profissional dos efeitos PCA salva\n")
  }
  
}, error = function(e) {
  cat("❌ Erro nas visualizações profissionais:", e$message, "\n")
})

### ------------------------------------------------------------------------------------- ###
### 5.1 VISUALIZAÇÃO FOCADA: MUSSISMILIA HISPIDA
### ------------------------------------------------------------------------------------- ###

cat("\n", rep("=", 70), "\n", sep = "")
cat("🎯 VISUALIZAÇÃO FOCADA: EFEITOS AMBIENTAIS EM MUSSISMILIA HISPIDA\n")
cat(rep("=", 70), "\n", sep = "")

# Diretório para os gráficos focados
mh_visualizations_dir <- file.path(output_dir, "Visualizacoes_Profissionais", "Foco_Mussismilia_Hispida")
dir.create(mh_visualizations_dir, showWarnings = FALSE)

# Os mesmos preditores PCA
pca_predictors <- c("pc1_magnitude", "pc2_magnitude", "pc1_variability", "pc2_variability")

mh_plots <- list()

for(pred in pca_predictors) {
  cat("  Processando efeito de", pred, "para M. hispida\n")
  
  tryCatch({
    # 1. Gera os efeitos para TODAS as espécies (como antes)
    cond_effect_full <- conditional_effects(jsdm_brms_fit_final, 
                                            effects = pred, 
                                            categorical = TRUE,
                                            plot = FALSE) # Importante: plot = FALSE para pegar os dados
    
    # 2. Extrai o dataframe de dados dos efeitos
    effect_df_full <- cond_effect_full[[1]]
    
    # 3. <<< A MÁGICA ACONTECE AQUI >>>
    # Filtra o dataframe para manter APENAS mussismilia_hispida
    effect_df_mh <- effect_df_full %>%
      filter(cats__ == "mussismilia_hispida")
    
    # Verifica se o filtro funcionou
    if(nrow(effect_df_mh) > 0) {
      # 4. Cria o gráfico usando APENAS os dados filtrados
      p_mh <- ggplot(effect_df_mh, aes(x = .data[[pred]], y = estimate__)) +
        geom_ribbon(aes(ymin = lower__, ymax = upper__), 
                    alpha = 0.3, fill = "skyblue") + # Cor única para clareza
        geom_line(linewidth = 1, color = "darkblue") +
        labs(x = toupper(pred), 
             y = "Proporção Prevista de M. hispida",
             title = paste("Efeito de", toupper(pred))) +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(face = "bold", size = 10))
      
      mh_plots[[pred]] <- p_mh
      
      # Salva o gráfico individual
      ggsave(file.path(mh_visualizations_dir, paste0("efeito_", pred, "_em_mh.png")),
             p_mh, width = 7, height = 5, dpi = 300, bg = "white")
      
    } else {
      cat("    AVISO: Não foram encontrados dados de M. hispida para o preditor", pred, "\n")
    }
    
  }, error = function(e) {
    cat("    ❌ Erro ao gerar gráfico focado para", pred, ":", e$message, "\n")
  })
}

# Combina os 4 gráficos focados em um único painel
if(length(mh_plots) > 0) {
  combined_mh_plot <- wrap_plots(mh_plots, ncol = 2) +
    plot_annotation(
      title = "Respostas de Mussismilia hispida aos Gradientes Ambientais",
      subtitle = "Efeitos parciais do modelo JSDM Dirichlet (R² preditivo ~4.8%)",
      theme = theme(plot.title = element_text(face = "bold", size = 16),
                    plot.subtitle = element_text(size = 12))
    )
  
  ggsave(file.path(mh_visualizations_dir, "efeitos_combinados_em_mh.png"),
         combined_mh_plot, width = 12, height = 10, dpi = 300, bg = "white")
  
  cat("✅ Gráficos focados em M. hispida salvos com sucesso!\n")
}

### ------------------------------------------------------------------------------------- ###
### 5.2 VISUALIZAÇÕES INDIVIDUAIS PARA CADA ORGANISMO (TOP 20)
### ------------------------------------------------------------------------------------- ###

cat("\n", rep("=", 70), "\n", sep = "")
cat("✨ GERANDO FIGURAS INDIVIDUAIS (4 PAINÉIS) PARA CADA ORGANISMO\n")
cat(rep("=", 70), "\n", sep = "")

# Diretório para salvar as 20 figuras individuais
all_organisms_dir <- file.path(output_dir, "Visualizacoes_Profissionais", "Figuras_Por_Organismo")
dir.create(all_organisms_dir, showWarnings = FALSE)

# --- PASSO 1: OTIMIZAÇÃO - Pré-calcular os dataframes de efeitos ---
# Isso é muito mais rápido do que chamar conditional_effects 20x4 vezes.
# Chamamos apenas 4 vezes, uma para cada preditor.
pca_predictors <- c("pc1_magnitude", "pc2_magnitude", "pc1_variability", "pc2_variability")
precomputed_effects <- list()

cat("--- Otimização: Pré-calculando dataframes de efeitos para os 4 preditores ---\n")
for (pred in pca_predictors) {
  cat("  Calculando para:", pred, "...\n")
  tryCatch({
    cond_effect <- conditional_effects(jsdm_brms_fit_final, 
                                       effects = pred, 
                                       categorical = TRUE, 
                                       plot = FALSE)
    precomputed_effects[[pred]] <- cond_effect[[1]]
  }, error = function(e) {
    cat("    AVISO: Falha ao calcular efeitos para", pred, "-", e$message, "\n")
  })
}
cat("--- Pré-cálculo concluído ---\n")


# --- PASSO 2: Loop principal para gerar uma figura para cada organismo ---
for (org_name in organism_vars) {
  
  cat(paste0("\n--- Gerando figura de 4 painéis para: ", org_name, " ---\n"))
  
  # Lista para armazenar os 4 gráficos de um organismo
  organism_plots <- list()
  
  # Loop interno para criar um painel para cada preditor PCA
  for (pred in pca_predictors) {
    
    # Pega o dataframe pré-calculado
    effect_df_full <- precomputed_effects[[pred]]
    
    # Se o dataframe existir, continue
    if (!is.null(effect_df_full)) {
      
      # Filtra para o organismo atual
      effect_df_filtered <- effect_df_full %>% 
        filter(cats__ == org_name)
      
      # Cria o gráfico se houver dados para este organismo
      if (nrow(effect_df_filtered) > 0) {
        
        # Define um título mais curto para o painel
        panel_title <- case_when(
          pred == "pc1_magnitude" ~ "PC1 (Magnitude)",
          pred == "pc2_magnitude" ~ "PC2 (Magnitude)",
          pred == "pc1_variability" ~ "PC1 (Variabilidade)",
          pred == "pc2_variability" ~ "PC2 (Variabilidade)",
          TRUE ~ pred
        )
        
        p <- ggplot(effect_df_filtered, aes(x = .data[[pred]], y = estimate__)) +
          geom_ribbon(aes(ymin = lower__, ymax = upper__), 
                      alpha = 0.3, fill = "#0D0887FF") + # Usando uma cor do viridis
          geom_line(linewidth = 1, color = "#0D0887FF") +
          labs(
            x = NULL, # Remove o rótulo X para um visual mais limpo no painel
            y = if(pred %in% c("pc1_magnitude", "pc1_variability")) "Proporção Prevista" else "",
            title = panel_title
          ) +
          theme_minimal(base_size = 11) +
          theme(
            plot.title = element_text(size = 10, hjust = 0.5),
            axis.text.x = element_text(angle = 45, hjust = 1)
          )
        
        organism_plots[[pred]] <- p
      }
    }
  } # Fim do loop de preditores
  
  # Combina os 4 painéis em uma única figura
  if (length(organism_plots) == 4) {
    
    # Limpa o nome do organismo para usar em nomes de arquivos
    safe_org_name <- gsub("[^A-Za-z0-9_]", "_", org_name)
    
    # Usa patchwork para combinar os gráficos
    combined_figure <- wrap_plots(organism_plots, ncol = 2) +
      plot_annotation(
        title = paste("Efeitos Ambientais na Abundância de:", org_name),
        subtitle = "Efeitos parciais do modelo JSDM Dirichlet",
        theme = theme(
          plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5)
        )
      )
    
    # Salva a figura combinada
    file_name <- paste0("efeitos_pca_para_", safe_org_name, ".png")
    ggsave(
      filename = file.path(all_organisms_dir, file_name),
      plot = combined_figure,
      width = 10,
      height = 8,
      dpi = 300,
      bg = "white"
    )
    cat(paste0("  ✅ Figura salva em: ", file.path(all_organisms_dir, file_name), "\n"))
    
  } else {
    cat(paste0("  AVISO: Não foi possível gerar a figura completa para ", org_name, ". Apenas ", length(organism_plots), " painéis foram criados.\n"))
  }
  
} # Fim do loop de organismos

cat("\n", rep("=", 70), "\n", sep = "")
cat("✅ Processo de geração de figuras individuais concluído!\n")
cat("📁 Verifique o diretório:", all_organisms_dir, "\n")


### ------------------------------------------------------------------------------------- ###
### 6. RELATÓRIO FINAL DE PUBLICAÇÃO
### ------------------------------------------------------------------------------------- ###

cat("\n", rep("=", 70), "\n", sep = "")
cat("📋 RELATÓRIO FINAL PARA PUBLICAÇÃO\n")
cat(rep("=", 70), "\n", sep = "")

# Gerar relatório resumido profissional
sink(file.path(output_dir, "RELATORIO_FINAL_PUBLICACAO.txt"))
cat("=== RELATÓRIO FINAL: JSDM DIRICHLET PARA PUBLICAÇÃO ===\n\n")
cat("DATA:", format(Sys.Date(), "%Y-%m-%d"), "\n\n")

cat("--- RESUMO DO MODELO ---\n")
cat("• Tipo: JSDM Dirichlet com efeitos aleatórios aninhados\n")
cat("• Organismos:", length(organism_vars), "espécies\n")
cat("• Amostras:", nrow(data_for_brms), "observações\n")
cat("• Sítios:", metadata$n_sites, "| Habitats:", metadata$n_habitats, "\n")
cat("• Preditores: 4 variáveis PCA (magnitude e variabilidade PC1/PC2)\n\n")

cat("--- DESEMPENHO DO MODELO ---\n")
if(exists("performance_metrics")) {
  cat(sprintf("• Mussismilia hispida - R²: %.1f%%\n", performance_metrics$r_squared * 100))
  cat(sprintf("• Correlação: %.3f\n", performance_metrics$correlation))
  cat(sprintf("• RMSE: %.4f\n", performance_metrics$rmse))
}

cat("\n--- CONVERGÊNCIA ---\n")
cat(sprintf("• Rhat máximo: %.4f\n", max(model_summary$fixed$Rhat)))
cat(sprintf("• Bulk ESS mínimo: %.0f\n", min(model_summary$fixed$Bulk_ESS)))
cat(sprintf("• Tail ESS mínimo: %.0f\n", min(model_summary$fixed$Tail_ESS)))

cat("\n--- PRÓXIMOS PASSOS PARA ANÁLISE ---\n")
cat("1. Verificar diagnósticos em: Diagnosticos_Convergencia/\n")
cat("2. Analisar visualizações em: Visualizacoes_Profissionais/\n")
cat("3. Revisar análise detalhada de M. hispida\n")
cat("4. Preparar tabelas e figuras para o manuscrito\n")

sink()

cat("\n", rep("=", 70), "\n", sep = "")
cat("✅ SCRIPT 2 CONCLUÍDO! ANÁLISES PROFISSIONAIS FINALIZADAS\n")
cat("📁 Todos os resultados salvos em:", output_dir, "\n")
cat("🎯 Relatório final gerado: 'RELATORIO_FINAL_PUBLICACAO.txt'\n")
cat(rep("=", 70), "\n", sep = "")