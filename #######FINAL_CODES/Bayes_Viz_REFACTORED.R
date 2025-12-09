### ==================================================================================== ###
### === SCRIPT DE VISUALIZAÇÃO PROFISSIONAL (PÓS-PROCESSAMENTO) === ###
### === Foca exclusivamente em Estética, Tidybayes e Comunicação === ###
### ==================================================================================== ###

rm(list = ls())
gc()

# --- 1. PACOTES DE ELITE PARA VISUALIZAÇÃO ---
libs <- c("brms", "ggplot2", "dplyr", "tidybayes", "ggdist", "bayesplot", "patchwork", "cowplot")
# Se não tiver algum: install.packages(c("tidybayes", "ggdist", "cowplot"))
invisible(lapply(libs, library, character.only = TRUE))

# --- 2. CONFIGURAÇÃO ESTÉTICA (TEMA DE REVISTA) ---
theme_nature <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(family = "sans", color = "black"),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(face = "bold", size = 14),
      plot.title = element_text(face = "bold", size = 16, hjust = 0),
      plot.subtitle = element_text(size = 12, color = "gray40"),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      strip.background = element_rect(fill = "gray95", color = NA),
      strip.text = element_text(face = "bold", size = 12)
    )
}
theme_set(theme_nature())

# Cores acessíveis e profissionais (Paleta Okabe-Ito)
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# ==================================================================================== #
# === PASSO 1: CARREGUE SEU MODELO VENCEDOR AQUI === #
# ==================================================================================== #

# Caminho para o arquivo .rds gerado pelo script de modelagem
# EXEMPLO: Substitua pelo caminho real do seu arquivo vencedor
model_path <- "C:/Users/rbfra/OneDrive/.../GLOBAL_WINNER_....rds" 

if (!file.exists(model_path)) stop("Arquivo do modelo não encontrado! Verifique o caminho.")

cat("--- Carregando o modelo vencedor (isso pode levar alguns segundos)... \n")
model <- readRDS(model_path)
cat("--- Modelo carregado com sucesso!\n")

# Detectar se é ZOIB ou Gaussiano para ajustar legendas
is_zoib <- family(model)$family == "zero_one_inflated_beta"
response_label <- if(is_zoib) "Cobertura (Prop)" else "Taxa/Saúde"

# ==================================================================================== #
# === VISUALIZAÇÃO A: FOREST PLOT COM "HALF-EYE" (EFEITOS FIXOS) === #
# ==================================================================================== #
# Mostra a distribuição completa da incerteza dos coeficientes, não apenas linhas.

cat("--- Gerando Forest Plot (Half-Eye)...\n")

# Extrair efeitos fixos (remove Intercepto e sds de spline que poluem)
variables_to_plot <- get_variables(model)[grep("^b_", get_variables(model))]
variables_to_plot <- variables_to_plot[!grepl("Intercept", variables_to_plot)]

if(length(variables_to_plot) > 0) {
  p_forest <- model %>%
    gather_draws(!!sym(variables_to_plot[1]):!!sym(tail(variables_to_plot, 1))) %>%
    ggplot(aes(y = .variable, x = .value, fill = after_stat(x > 0))) +
    
    # O "Half-Eye" combina densidade + intervalo
    stat_halfeye(
      .width = c(0.89, 0.95),  # Intervalos de Credibilidade Padrão Bayesiano
      alpha = 0.8,
      point_interval = median_qi
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("gray60", "#0072B2"), guide = "none") + # Azul se positivo, Cinza se negativo
    labs(
      title = "Magnitude dos Efeitos",
      subtitle = "Distribuição Posterior dos Coeficientes Fixos",
      x = "Estimativa do Efeito (Escala Padronizada)",
      y = NULL
    ) +
    theme_nature()
  
  print(p_forest)
  # ggsave("Forest_Plot_Profissional.png", p_forest, width = 10, height = 6, dpi = 300)
} else {
  cat("Modelo dominado por Splines (sem efeitos lineares fixos simples para plotar).\n")
}

# ==================================================================================== #
# === VISUALIZAÇÃO B: CONDITIONAL EFFECTS COM "RIBBON" (SPLINES/REGRESSÃO) === #
# ==================================================================================== #
# O gráfico de linhas com bandas de confiança sombreadas.

cat("--- Gerando Ribbon Plots para Efeitos Marginais...\n")

# Calcular efeitos
ce <- conditional_effects(model, robust = TRUE, method = "posterior_epred") # method posterior_epred é mais suave

plot_list <- list()

for(eff in names(ce)) {
  data_eff <- ce[[eff]]
  
  # Identificar x e y
  x_var <- names(data_eff)[1] # A primeira coluna é sempre a variável de efeito
  
  p_ribbon <- ggplot(data_eff, aes_string(x = x_var, y = "estimate__")) +
    # Banda de incerteza (95%)
    geom_ribbon(aes_string(ymin = "lower__", ymax = "upper__"), fill = "#0072B2", alpha = 0.2) +
    # Linha média
    geom_line(color = "#0072B2", size = 1.2) +
    labs(
      title = paste("Efeito de", x_var),
      y = paste("Predição:", response_label),
      x = x_var
    ) +
    theme_nature()
  
  # Se houver pontos de dados reais (rug plot), adicionar
  # (Opcional, remove se poluir muito)
  # p_ribbon <- p_ribbon + geom_rug(data = model$data, aes_string(x = x_var), alpha=0.1, sides="b")
  
  plot_list[[eff]] <- p_ribbon
}

# Juntar todos em um painel
p_effects <- wrap_plots(plot_list, ncol = 2) + 
  plot_annotation(
    title = "Efeitos Marginais Condicionais",
    subtitle = "Linha: Mediana Posterior | Banda: Intervalo de Credibilidade 95%",
    theme = theme(plot.title = element_text(face="bold", size=18))
  )

print(p_effects)
# ggsave("Effects_Ribbon_Profissional.png", p_effects, width = 12, height = 8, dpi = 300)


# ==================================================================================== #
# === VISUALIZAÇÃO C: DIAGNÓSTICOS DE CONVERGÊNCIA (ESTILO REPORT) === #
# ==================================================================================== #

cat("--- Gerando Diagnósticos Visuais...\n")

# Trace Plot (Estilo Clean)
# Pega apenas 4 parâmetros aleatórios para não poluir
params_to_trace <- sample(variables(model)[1:min(4, length(variables(model)))], 4)
p_trace <- mcmc_trace(model, pars = params_to_trace) + 
  theme_nature() + 
  scale_color_manual(values = okabe_ito) +
  labs(title = "Trace Plots (Convergência)", subtitle = "Cadeias devem parecer 'lagartas peludas' misturadas")

# Autocorrelação
p_acf <- mcmc_acf(model, pars = params_to_trace) + 
  theme_nature() +
  labs(title = "Autocorrelação", subtitle = "Deve cair para zero rapidamente")

p_diag <- (p_trace / p_acf)
print(p_diag)
# ggsave("Diagnosticos_Profissional.png", p_diag, width = 10, height = 8, dpi = 300)


# ==================================================================================== #
# === VISUALIZAÇÃO D: PPC (POSTERIOR PREDICTIVE CHECK) === #
# ==================================================================================== #

cat("--- Gerando PPC (Densidade)...\n")

p_ppc <- pp_check(model, ndraws = 100, type = "dens_overlay") +
  scale_color_manual(values = c("black", "#56B4E9")) + # Preto (Obs) vs Azul Claro (Sim)
  labs(
    title = "Validação Preditiva (PPC)",
    subtitle = "Preto: Dados Observados | Azul: Simulações do Modelo",
    x = response_label,
    y = "Densidade"
  ) +
  theme_nature()

print(p_ppc)
# ggsave("PPC_Profissional.png", p_ppc, width = 8, height = 6, dpi = 300)

cat("\n=== SCRIPT FINALIZADO ===\n")