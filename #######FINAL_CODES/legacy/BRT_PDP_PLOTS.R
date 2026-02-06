################################################################################
##### SCRIPT A POSTERIORI PARA GERAÇÃO DE GRÁFICOS DE MODELOS BRT (V. CORRIGIDA) ####
#
# DESCRIÇÃO:
# Este script carrega um objeto de modelo BRT salvo (best_brt_object.rds)
# de uma análise anterior e gera um conjunto completo de saídas visuais.
#
# INSTRUÇÕES:
# 1. Configure os 3 campos na Seção 1 para apontar para o modelo desejado.
# 2. Execute o script.
#
################################################################################

### 1. CONFIGURAÇÃO: APONTE PARA O MELHOR MODELO -----------------------------

# --- PASSO 1: Especifique o diretório principal da sua análise de comparação ---
### CORREÇÃO PRINCIPAL AQUI ###
# Certifique-se de que este é o diretório de saída do script de modelagem (BRT_Health_V4).
main_output_folder <- "C:/Users/rbfra/OneDrive/BRT_Health_Growth_CV_30"

# --- PASSO 2: Especifique a pasta da variável resposta ---
# (Ex: "RGR", "HEALTH_PC1", "HEALTH_PC2")
response_variable_folder <- "HEALTH_PC1"

# --- PASSO 3: Especifique o nome da pasta do modelo que você quer analisar ---
# (Copie o nome da pasta do melhor modelo identificado na análise de comparação)
model_name_folder <- "model_main_with_depth_with_hab"


### 2. Pacotes e Configurações Globais -----------------------------------------
libs <- c("gbm", "dismo", "dplyr", "ggplot2", "pdp", "gridExtra", "purrr", "parallel",
          "lmerTest", "performance", "spdep", "readr", "ggeffects", "car", "plotly",
          "plot3D", "forcats", "tibble", "htmlwidgets")
invisible(lapply(libs, library, character.only = TRUE))

set.seed(123)
USE_PARALLEL <- TRUE
N_CORES <- if (USE_PARALLEL) max(1, parallel::detectCores() - 1) else 1
N_BOOTSTRAP <- 100 # Número de iterações de bootstrap para os PDPs


### 3. Funções Utilitárias (CONJUNTO COMPLETO) ---------------------------------
# (Nenhuma mudança necessária aqui, as funções estão corretas)
pdp_bootstrap_fixed_refit <- function(model, data, var, optimal_iter, response_name, model_params, n_grid = 50, n_bootstrap = N_BOOTSTRAP, sample_frac = 0.7, n_cores = N_CORES) {
  grid_vals <- seq(min(data[[var]], na.rm = TRUE), max(data[[var]], na.rm = TRUE), length.out = n_grid); baseline <- mean(predict(model, newdata = data, n.trees = optimal_iter)); pdp_main_raw <- sapply(grid_vals, function(v) { X_tmp <- data; X_tmp[[var]] <- v; mean(predict(model, newdata = X_tmp, n.trees = optimal_iter)) }); pdp_main <- pdp_main_raw - baseline
  bootstrap_iteration <- function(i) {
    set.seed(12345 + i); idx <- sample(seq_len(nrow(data)), size = floor(sample_frac * nrow(data)), replace = TRUE); Xb <- data[idx, , drop = FALSE]; form <- as.formula(paste(response_name, "~", paste(model$var.names, collapse = "+"))); tc <- model_params$tc; lr <- model_params$lr; bf <- model_params$bf
    refit_model <- try(gbm(form, data = Xb, distribution = "gaussian", n.trees = model$n.trees, interaction.depth = tc, shrinkage = lr, bag.fraction = bf, n.minobsinnode = model$n.minobsinnode, keep.data = TRUE, verbose = FALSE), silent = TRUE)
    if (inherits(refit_model, "try-error")) return(rep(NA_real_, n_grid))
    it_refit <- optimal_iter; if(it_refit > length(refit_model$train.error) || it_refit < 1) { it_refit <- length(refit_model$train.error) }; if(it_refit < 1) return(rep(NA_real_, n_grid))
    b_baseline <- mean(predict(refit_model, newdata = Xb, n.trees = it_refit)); preds_on_grid <- sapply(grid_vals, function(v) { X_tmp <- Xb; X_tmp[[var]] <- v; mean(predict(refit_model, newdata = X_tmp, n.trees = it_refit)) }); return(preds_on_grid - b_baseline)
  }
  if (USE_PARALLEL && n_cores > 1) {
    cl <- parallel::makeCluster(n_cores); on.exit(parallel::stopCluster(cl), add = TRUE); parallel::clusterEvalQ(cl, { library(gbm) }); parallel::clusterExport(cl, varlist = c("data", "var", "optimal_iter", "response_name", "model_params", "sample_frac", "n_grid", "grid_vals", "bootstrap_iteration", "model"), envir = environment()); boot_curves <- parallel::parSapply(cl, 1:n_bootstrap, bootstrap_iteration)
  } else { boot_curves <- sapply(1:n_bootstrap, bootstrap_iteration) }
  if (!is.matrix(boot_curves) || nrow(boot_curves) != n_grid || ncol(boot_curves) != n_bootstrap) { lower <- upper <- rep(NA_real_, n_grid) } else { lower <- apply(boot_curves, 1, quantile, 0.025, na.rm = TRUE); upper <- apply(boot_curves, 1, quantile, 0.975, na.rm = TRUE) }
  df_plot <- data.frame(x = grid_vals, y = pdp_main, lower = lower, upper = upper); valid_y_range <- range(c(df_plot$y, df_plot$lower, df_plot$upper), na.rm = TRUE); y_max_abs <- max(abs(valid_y_range), na.rm = TRUE); y_lims <- if(is.finite(y_max_abs)) c(-y_max_abs * 1.1, y_max_abs * 1.1) else NULL
  p <- ggplot(df_plot, aes(x = x, y = y)) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey70", alpha = 0.4) + geom_line(linewidth = 1.0, colour = "blue") + geom_rug(data = data, aes(x = .data[[var]]), inherit.aes = FALSE, sides = "b", alpha = 0.1) + theme_classic(base_size = 12) + coord_cartesian(ylim = y_lims) + labs(x = var, y = "Efeito parcial") + theme(axis.title = ggplot2::element_text(size=11)); return(p)
}
pdp_categorical_refit <- function(model, data, var, optimal_iter, response_name, model_params, n_bootstrap = N_BOOTSTRAP, sample_frac = 0.7, n_cores = N_CORES, grouping_info_df = NULL, grouping_var_name = "REEF") {
  data[[var]] <- factor(data[[var]]); levels_var <- levels(data[[var]])
  bootstrap_iteration_cat <- function(i, current_level) {
    set.seed(54321 + i + match(current_level, levels_var) * n_bootstrap); idx <- sample(seq_len(nrow(data)), size = floor(sample_frac * nrow(data)), replace = TRUE); Xb <- data[idx, , drop = FALSE]; Xb[[var]] <- factor(Xb[[var]], levels = levels_var); form <- as.formula(paste(response_name, "~", paste(model$var.names, collapse = "+"))); tc <- model_params$tc; lr <- model_params$lr; bf <- model_params$bf
    refit_model <- try(gbm(form, data = Xb, distribution = "gaussian", n.trees = model$n.trees, interaction.depth = tc, shrinkage = lr, bag.fraction = bf, n.minobsinnode = model$n.minobsinnode, keep.data = TRUE, verbose = FALSE), silent = TRUE)
    if (inherits(refit_model, "try-error")) return(NA_real_)
    it_refit <- optimal_iter; if(it_refit > length(refit_model$train.error) || it_refit < 1) { it_refit <- length(refit_model$train.error) }; if(it_refit < 1) return(NA_real_)
    b_baseline <- mean(predict(refit_model, newdata = Xb, n.trees = it_refit)); d_boot <- Xb; d_boot[[var]] <- factor(current_level, levels = levels_var); pred_for_level <- mean(predict(refit_model, newdata = d_boot, n.trees = it_refit)); return(pred_for_level - b_baseline)
  }
  pdp_list <- list(); baseline <- mean(predict(model, newdata = data, n.trees = optimal_iter))
  for (l in levels_var) {
    dtmp <- data; dtmp[[var]] <- factor(l, levels = levels_var); pred_mean <- mean(predict(model, newdata = dtmp, n.trees = optimal_iter)) - baseline
    if (USE_PARALLEL && n_cores > 1) {
      cl <- parallel::makeCluster(n_cores); parallel::clusterEvalQ(cl, { library(gbm) }); parallel::clusterExport(cl, varlist = c("data", "var", "optimal_iter", "response_name", "model_params", "sample_frac", "levels_var", "bootstrap_iteration_cat", "model"), envir = environment()); boot_preds <- parallel::parLapply(cl, 1:n_bootstrap, bootstrap_iteration_cat, current_level = l); parallel::stopCluster(cl)
    } else { boot_preds <- lapply(1:n_bootstrap, bootstrap_iteration_cat, current_level = l) }
    boot_preds <- unlist(boot_preds); pdp_list[[as.character(l)]] <- data.frame(level = l, mean_pred = pred_mean, lower = quantile(boot_preds, 0.025, na.rm = TRUE), upper = quantile(boot_preds, 0.975, na.rm = TRUE))
  }
  if(length(pdp_list) == 0) return(NULL); pdp_df <- do.call(rbind, pdp_list)
  if (!is.null(grouping_info_df) && grouping_var_name %in% names(grouping_info_df) && var %in% names(grouping_info_df)) {
    grouping_info_df <- rename(grouping_info_df, "group_col" = all_of(grouping_var_name))
    pdp_df <- left_join(pdp_df, grouping_info_df, by = c("level" = var))
    pdp_df$group_col[is.na(pdp_df$group_col)] <- "Desconhecido"; pdp_df$group_col <- as.factor(pdp_df$group_col)
  } else { pdp_df$group_col <- "default" }
  pdp_df <- pdp_df %>% mutate(level = forcats::fct_reorder(level, mean_pred, .desc = TRUE))
  valid_y_range <- range(c(pdp_df$mean_pred, pdp_df$lower, pdp_df$upper), na.rm = TRUE); y_max_abs <- max(abs(valid_y_range), na.rm = TRUE); y_lims <- if(is.finite(y_max_abs)) c(-y_max_abs * 1.1, y_max_abs * 1.1) else NULL
  p <- ggplot(pdp_df, aes(x = level, y = mean_pred, color = group_col)) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, linewidth=0.6) + geom_point(size = 3) +
    scale_color_viridis_d(name = grouping_var_name) + theme_classic(base_size = 12) + coord_cartesian(ylim = y_lims) +
    labs(x = var, y = "Efeito parcial") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=8), axis.title = element_text(size=11),
          legend.position = if(length(unique(pdp_df$group_col)) <= 1) "none" else "right") 
  return(p)
}
plot_interaction_categorical_continuous <- function(model, data, n.trees, var_cont, var_cat, n_grid = 50) {
  data[[var_cat]] <- as.factor(data[[var_cat]]); levels_cat <- levels(data[[var_cat]]); grid_cont <- seq(min(data[[var_cont]], na.rm = TRUE), max(data[[var_cont]], na.rm = TRUE), length.out = n_grid); pdp_data_list <- list(); all_model_vars <- model$var.names
  for (level in levels_cat) {
    pred_grid <- data.frame(matrix(nrow = n_grid, ncol = 0)); pred_grid[[var_cont]] <- grid_cont; pred_grid[[var_cat]] <- factor(level, levels = levels_cat); other_vars <- setdiff(all_model_vars, c(var_cont, var_cat))
    for(v in other_vars) { if(is.numeric(data[[v]])) { pred_grid[[v]] <- mean(data[[v]], na.rm = TRUE) } else { pred_grid[[v]] <- factor(levels(data[[v]])[1], levels = levels(data[[v]])) } }
    pred_grid <- pred_grid[, all_model_vars, drop = FALSE]; predictions <- predict(model, newdata = pred_grid, n.trees = n.trees); pdp_data_list[[level]] <- data.frame(cont_value = grid_cont, prediction = predictions, cat_level = level)
  }
  pdp_df <- do.call(rbind, pdp_data_list)
  p <- ggplot(pdp_df, aes(x = cont_value, y = prediction, color = cat_level, group = cat_level)) + geom_line(linewidth = 1.2) + labs(title = paste("Interação:", var_cont, "vs", var_cat), x = var_cont, y = "Efeito Parcial (Predição)", color = var_cat) + theme_classic(base_size = 12)
  return(p)
}
plot_interaction_continuous_continuous <- function(model, data, n.trees, var1, var2, n_grid = 30) {
  grid1 <- seq(min(data[[var1]], na.rm = TRUE), max(data[[var1]], na.rm = TRUE), length.out = n_grid); grid2 <- seq(min(data[[var2]], na.rm = TRUE), max(data[[var2]], na.rm = TRUE), length.out = n_grid); pred_grid <- expand.grid(grid1, grid2); names(pred_grid) <- c(var1, var2); other_vars <- setdiff(model$var.names, c(var1, var2))
  for(v in other_vars) { if(is.numeric(data[[v]])) { pred_grid[[v]] <- mean(data[[v]], na.rm = TRUE) } else { pred_grid[[v]] <- factor(levels(data[[v]])[1], levels = levels(data[[v]])) } }
  predictions <- predict(model, newdata = pred_grid[, model$var.names], n.trees = n.trees); z_matrix <- matrix(predictions, nrow = n_grid, ncol = n_grid)
  p <- plotly::plot_ly(x = ~grid1, y = ~grid2, z = ~z_matrix, type = "surface") %>% plotly::layout(title = paste("Interação:", var1, "vs", var2), scene = list(xaxis = list(title = var1), yaxis = list(title = var2), zaxis = list(title = "Efeito Parcial (Predição)")))
  return(p)
}
plot_interaction_2d_heatmap <- function(model, data, n.trees, var1, var2, n_grid = 50) {
  grid1 <- seq(min(data[[var1]], na.rm = TRUE), max(data[[var1]], na.rm = TRUE), length.out = n_grid); grid2 <- seq(min(data[[var2]], na.rm = TRUE), max(data[[var2]], na.rm = TRUE), length.out = n_grid); pred_grid <- expand.grid(X1 = grid1, X2 = grid2); names(pred_grid) <- c(var1, var2); other_vars <- setdiff(model$var.names, c(var1, var2))
  for(v in other_vars) { if(is.numeric(data[[v]])) { pred_grid[[v]] <- mean(data[[v]], na.rm = TRUE) } else { pred_grid[[v]] <- factor(levels(data[[v]])[1], levels = levels(data[[v]])) } }
  pred_grid <- pred_grid[, model$var.names]; predictions <- predict(model, newdata = pred_grid, n.trees = n.trees); pred_grid$prediction <- predictions
  p <- ggplot(pred_grid, aes(x = .data[[var1]], y = .data[[var2]], fill = prediction)) + geom_tile() + scale_fill_viridis_c(option = "plasma") + labs(title = paste("Interação:", var1, "e", var2), x = var1, y = var2, fill = "Predição") + theme_classic(base_size = 14)
  return(p)
}
plot_interaction_3d_static <- function(model, data, n.trees, var1, var2, n_grid = 50, theta = 30, phi = 20) {
  grid1 <- seq(min(data[[var1]], na.rm = TRUE), max(data[[var1]], na.rm = TRUE), length.out = n_grid); grid2 <- seq(min(data[[var2]], na.rm = TRUE), max(data[[var2]], na.rm = TRUE), length.out = n_grid); pred_grid <- expand.grid(X1 = grid1, X2 = grid2); names(pred_grid) <- c(var1, var2); other_vars <- setdiff(model$var.names, c(var1, var2))
  for(v in other_vars) { if(is.numeric(data[[v]])) { pred_grid[[v]] <- mean(data[[v]], na.rm = TRUE) } else { pred_grid[[v]] <- factor(levels(data[[v]])[1], levels = levels(data[[v]])) } }
  pred_grid <- pred_grid[, model$var.names]; predictions <- predict(model, newdata = pred_grid, n.trees = n.trees); z_matrix <- matrix(predictions, nrow = n_grid, ncol = n_grid)
  plot3D::persp3D(x = grid1, y = grid2, z = z_matrix, xlab = var1, ylab = var2, zlab = "Predição", main = paste("Interação:", var1, "e", var2), theta = theta, phi = phi, ticktype = "detailed", colkey = list(side = 4, length = 0.5))
}


### 4. CARREGAMENTO E PROCESSAMENTO DO MODELO ----------------------------------

# Construir o caminho completo para o arquivo do modelo
model_file_path <- file.path(main_output_folder, "Outputs", response_variable_folder, model_name_folder, "BRT_Model_Outputs", "best_brt_object.rds")

cat("--- Tentando carregar o modelo de:", model_file_path, "---\n")

if (!file.exists(model_file_path)) {
  stop("ERRO: O arquivo do modelo não foi encontrado. Verifique os 3 campos de configuração na Seção 1.")
}

best_brt <- readRDS(model_file_path)
cat("--- Modelo carregado com sucesso! ---\n\n")

# Extrair informações essenciais do objeto carregado
model <- best_brt$model
optimal_iter <- best_brt$iter
df_model <- best_brt$data
model_params <- as.list(best_brt$params)
resp <- names(df_model)[1] # A variável resposta é a primeira coluna nos dados salvos

# Criar um diretório de saída para os gráficos
plots_out_dir <- file.path(main_output_folder, "Outputs", response_variable_folder, model_name_folder, "Best_Model_Visuals_A_Posteriori")
dir.create(plots_out_dir, showWarnings = FALSE, recursive = TRUE)
cat("--- Gráficos serão salvos em:", plots_out_dir, "---\n")

# Carregar dados de coordenadas para colorir gráficos de SITE
coords_path_new <- "C:/Users/rbfra/OneDrive/########CEBIMAR/####PROJETOS/#####Coral trade offs"
coords_file_new <- file.path(coords_path_new, "sites_list_full.csv")
coords_df <- readr::read_delim(coords_file_new, delim = ";", show_col_types = FALSE) %>%
  rename(REEF = Reef_name, SITE = Site_name) %>%
  dplyr::select(REEF, SITE) %>%
  distinct(SITE, .keep_all = TRUE)


### 5. GERAÇÃO DE GRÁFICOS DETALHADOS ------------------------------------------

# --- 5.1 Importância das Variáveis ---
cat("\n----- 1. Calculando a Importância das Variáveis -----\n")
importance_df <- try(summary(model, n.trees=optimal_iter, plotit=F), silent=T)
if (!inherits(importance_df, "try-error") && nrow(importance_df) > 0) {
  write.csv(importance_df, file.path(plots_out_dir, "importance_brt.csv"), row.names = FALSE)
  cat("  - Tabela de importância salva com sucesso.\n")
} else {
  cat("  - AVISO: Não foi possível calcular a importância das variáveis.\n")
}

# --- 5.2 Partial Dependence Plots (PDPs) com Refit ---
cat("\n----- 2. Gerando Partial Dependence Plots (PDPs) com Bootstrap -----\n")
predictors_used <- if(!is.null(importance_df) && nrow(importance_df)>0) as.character(importance_df$var) else model$var.names
plot_list <- list()
for (var in predictors_used) {
  cat(paste0("    - Processando PDP para variável: ", var, "...\n"))
  plot_obj <- NULL
  if (is.numeric(df_model[[var]])) { 
    plot_obj <- try(pdp_bootstrap_fixed_refit(model, df_model, var, optimal_iter, resp, model_params), silent = TRUE) 
  } else if (is.factor(df_model[[var]])) {
    if (var == "SITE" && !is.null(coords_df)) {
      plot_obj <- try(pdp_categorical_refit(model, df_model, var, optimal_iter, resp, model_params, grouping_info_df = coords_df, grouping_var_name = "REEF"), silent = TRUE)
    } else {
      plot_obj <- try(pdp_categorical_refit(model, df_model, var, optimal_iter, resp, model_params), silent = TRUE)
    }
  }
  if (inherits(plot_obj, "ggplot")) { 
    plot_obj <- plot_obj + labs(title = var)
    plot_list[[var]] <- plot_obj 
  } else {
    cat(paste0("    - AVISO: Falha ao gerar PDP para '", var, "'. Erro: ", attr(plot_obj, "condition")$message, "\n"))
  }
}
if(length(plot_list) > 0) { 
  ncol_combo <- min(ceiling(sqrt(length(plot_list))), 4)
  combo_width <- 5 * ncol_combo
  combo_height <- 4 * ceiling(length(plot_list) / ncol_combo)
  title_text <- paste("Efeitos Parciais (BRT) para Resposta:", resp, "| Modelo:", model_name_folder)
  combined_grob <- gridExtra::arrangeGrob(grobs = plot_list, ncol = ncol_combo, top = grid::textGrob(title_text, gp=grid::gpar(fontsize=16)))
  ggsave(filename = file.path(plots_out_dir, "PDPs_Combined_Refit.png"), plot = combined_grob, width = combo_width, height = combo_height, dpi = 300, limitsize = FALSE)
  cat("  - Figura combinada de PDPs salva com sucesso.\n")
}

# --- 5.3 Análise de Interações ---
cat("\n----- 3. Analisando e Plotando Interações -----\n")
if (model_params$tc < 2) {
  cat("  - Nenhuma interação avaliada (tc < 2 no modelo treinado).\n")
} else {
  
  ### MELHORIA AQUI: Lógica mais robusta para recriar a lista de interações ###
  interactions_to_test <- list()
  # Se 'ARCH' está nos preditores, assumimos que as interações com ele são as de interesse.
  if ("ARCH" %in% model$var.names) {
    cat("  - Modelo de interação com ARCH detectado. Testando interações relevantes...\n")
    # Pega todos os outros preditores do modelo para testar a interação com ARCH
    other_predictors <- setdiff(model$var.names, "ARCH")
    interactions_to_test <- lapply(other_predictors, function(p) c(p, "ARCH"))
  } else {
    cat("  - Modelo de efeitos principais. Testando todas as interações de 2 vias entre preditores.\n")
    if(length(model$var.names) >= 2) {
      comb_matrix <- combn(model$var.names, 2, simplify = TRUE)
      interactions_to_test <- lapply(seq_len(ncol(comb_matrix)), function(i) comb_matrix[, i])
    }
  }
  
  if (length(interactions_to_test) == 0) {
    cat("  - Nenhuma interação válida para testar foi encontrada.\n")
  } else {
    cat(paste("  - Testando a força de", length(interactions_to_test), "interações potenciais...\n"))
    interaction_strengths <- list()
    for (pair in interactions_to_test) {
      if (all(c(pair[1], pair[2]) %in% model$var.names)) {
        strength <- try(gbm::interact.gbm(model, data=df_model, i.var=pair, n.trees=optimal_iter), silent=T)
        if (!inherits(strength, "try-error")) { 
          interaction_strengths[[paste(pair, collapse = ":")]] <- data.frame(var1.names = pair[1], var2.names = pair[2], strength = strength) 
        }
      }
    }
    
    if (length(interaction_strengths) > 0) {
      full_interaction_list <- do.call(rbind, interaction_strengths) %>% arrange(desc(strength))
      write.csv(full_interaction_list, file.path(plots_out_dir, "interactions_ranked.csv"), row.names=F)
      cat("  - Tabela de interações ranqueadas salva.\n")
      print(head(full_interaction_list, 5))
      
      num_to_plot <- min(nrow(full_interaction_list), 3)
      if (num_to_plot > 0) {
        cat(paste0("  - Gerando gráficos para as ", num_to_plot, " interações mais fortes...\n"))
        for (i in 1:num_to_plot) {
          top_int <- full_interaction_list[i, ]; var1 <- as.character(top_int$var1.names); var2 <- as.character(top_int$var2.names)
          cat(paste0("    [#", i, "] Plotando '", var1, "' vs '", var2, "'...\n"))
          is_v1_num <- is.numeric(df_model[[var1]]); is_v2_num <- is.numeric(df_model[[var2]])
          var1_short <- substr(gsub("_scaled", "", var1), 1, 10); var2_short <- substr(gsub("_scaled", "", var2), 1, 10)
          file_base <- paste0("Int_BRT_R", i, "_", var1_short, "_vs_", var2_short)
          
          if (is_v1_num && is_v2_num) {
            p_hm <- plot_interaction_2d_heatmap(model, df_model, optimal_iter, var1, var2)
            ggsave(file.path(plots_out_dir, paste0(file_base, "_HM.png")), p_hm, width=8, height=6, dpi=300)
            p_int <- plot_interaction_continuous_continuous(model, df_model, optimal_iter, var1, var2)
            try(htmlwidgets::saveWidget(p_int, file.path(plots_out_dir, paste0(file_base, ".html")), selfcontained=T))
            view_angles <- list(v1_fr=list(theta=45,phi=25), v2_fl=list(theta=135,phi=25), v3_bl=list(theta=225,phi=25), v4_br=list(theta=315,phi=25))
            for (view_name in names(view_angles)) {
              png(file.path(plots_out_dir, paste0(file_base, "_3D_", view_name, ".png")), width=8, height=7, units="in", res=300)
              plot_interaction_3d_static(model, df_model, optimal_iter, var1, var2, theta=view_angles[[view_name]]$theta, phi=view_angles[[view_name]]$phi)
              dev.off()
            }
          } else if ((is_v1_num && !is_v2_num) || (!is_v1_num && is_v2_num)) {
            var_c <- if(is_v1_num) var1 else var2; var_f <- if(is_v1_num) var2 else var1
            p <- plot_interaction_categorical_continuous(model, df_model, optimal_iter, var_cont=var_c, var_cat=var_f)
            ggsave(file.path(plots_out_dir, paste0(file_base, ".png")), p, width=8, height=6, dpi=300)
          } else {
            cat("    - Pulando (ambas categóricas).\n")
          }
        }
      }
    } else {
      cat("  - Nenhuma interação pôde ser calculada.\n")
    }
  }
}

cat("\n\n===== GERAÇÃO DE GRÁFICOS A POSTERIORI FINALIZADA. Verifique o diretório de saída. =====\n")