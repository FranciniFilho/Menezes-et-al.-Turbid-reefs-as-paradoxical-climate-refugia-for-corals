################################################################################
##### BRT_Abundancia_V6_Bayesian_Comparison.R (VERSÃO CORRIGIDA E COMPLETA) ####
#
# DESCRIÇÃO:
# Versão adaptada para espelhar a estrutura de análise do pipeline Bayesiano
# (PIPELINE_BAYESIANO_V8.5) para Mussismilia hispida.
#  - Foca exclusivamente em M. hispida.
#  - Usa os mesmos dados de entrada do script Bayesiano.
#  - Testa a mesma lista de 13 modelos candidatos (com/sem DEPTH, com/sem HAB,
#    e interações com ARCH).
#  - Organiza as saídas por modelo candidato para comparação direta.
#  - Gera uma tabela final comparando a performance de todos os modelos.
#
# AUTOR: Adaptado por IA com base nos scripts de R.B. Francini-Filho
# DATA: [Data Atual]
#
################################################################################

### 1. Pacotes -------------------------------------------------------------------
libs <- c("gbm", "dismo", "dplyr", "ggplot2", "pdp", "gridExtra", "purrr", "parallel",
          "lmerTest", "performance", "spdep", "readr", "ggeffects", "car", "plotly",
          "plot3D", "forcats", "tibble", "htmlwidgets", "tidyr")
invisible(lapply(libs, library, character.only = TRUE))


### 2. DADOS E VARIÁVEIS ---------------------------------------------------------
cat("--- [1/4] Carregando dados de Abundância e Coordenadas...\n")
abundance_file <- "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/########NEW RESULTS/#####output_local_PCA_CV_30_FINAL/dados_abundancia_integrados_long_format.csv"
tryCatch({
  data_raw <- read.csv2(abundance_file, stringsAsFactors = TRUE)
}, error = function(e) {
  data_raw <- read.csv(abundance_file, stringsAsFactors = TRUE)
})

coords_path_new <- "C:/Users/rbfra/OneDrive/########CEBIMAR/####PROJETOS/#####Coral trade offs"
coords_file_new <- file.path(coords_path_new, "sites_list_full.csv")
coords_df <- readr::read_delim(coords_file_new, delim = ";", show_col_types = FALSE) %>%
  rename(REEF = Reef_name, SITE = Site_name, LAT = Latitude, LONG = Longitude) %>%
  dplyr::select(REEF, SITE, LAT, LONG) %>%
  filter(!is.na(LAT) & !is.na(LONG)) %>%
  distinct(SITE, .keep_all = TRUE) %>%
  mutate(SITE = factor(SITE), REEF = factor(REEF))
site_to_reef_map <- coords_df %>% dplyr::select(SITE, REEF) %>% distinct()

cat("--- [2/4] Filtrando e preparando dados para M. hispida...\n")
colnames(data_raw) <- toupper(colnames(data_raw))
data <- data_raw %>%
  filter(ORGANISMO == "MUSSISMILIA_HISPIDA") %>%
  mutate(COBERTURA = as.numeric(gsub(",", ".", as.character(COBERTURA)))) %>%
  drop_na(COBERTURA)

if ("DEPTH_M" %in% names(data)) {
  if ("DEPTH" %in% names(data)) data <- data %>% dplyr::select(-DEPTH)
  data <- data %>% rename(DEPTH = DEPTH_M)
}
data$HAB <- factor(data$HAB, levels = c("PA", "RR", "TP"))

cat("--- [3/4] Escalonando preditores contínuos (PCs e Profundidade)...\n")
pca_predictors <- c("PC1_MAGNITUDE", "PC2_MAGNITUDE", "PC1_VARIABILITY", "PC2_VARIABILITY")
continuous_predictors_to_scale <- c(pca_predictors, "DEPTH")

data_scaled <- data
for(col in continuous_predictors_to_scale) {
  if(col %in% names(data_scaled)) {
    data_scaled[[paste0(col, "_scaled")]] <- scale(as.numeric(data_scaled[[col]]))[,1]
  }
}
data <- data_scaled
cat(sprintf("    -> Total de observações para M. hispida: %d\n", nrow(data)))

cat("--- [4/4] Definindo modelos candidatos para a comparação...\n")
pca_predictors_scaled <- paste0(pca_predictors, "_scaled")
depth_term_scaled <- "DEPTH_scaled"
hab_term <- "HAB"
arch_term <- "ARCH"
site_term <- "SITE" 
variability_predictors_scaled <- grep("_VARIABILITY_scaled", pca_predictors_scaled, value = TRUE)
interactions_arch <- lapply(variability_predictors_scaled, function(p) c(p, arch_term))

candidate_models <- list(
  "model_null" = list(predictors = c(site_term), interactions = list()),
  "model_linear_no_depth_no_hab" = list(predictors = c(pca_predictors_scaled, site_term), interactions = list()),
  "model_linear_with_depth_no_hab" = list(predictors = c(depth_term_scaled, pca_predictors_scaled, site_term), interactions = list()),
  "model_linear_no_depth_with_hab" = list(predictors = c(hab_term, pca_predictors_scaled, site_term), interactions = list()),
  "model_linear_with_depth_with_hab" = list(predictors = c(depth_term_scaled, hab_term, pca_predictors_scaled, site_term), interactions = list()),
  "model_spline_no_depth_no_hab" = list(predictors = c(pca_predictors_scaled, site_term), interactions = list()),
  "model_spline_with_depth_no_hab" = list(predictors = c(depth_term_scaled, pca_predictors_scaled, site_term), interactions = list()),
  "model_spline_no_depth_with_hab" = list(predictors = c(hab_term, pca_predictors_scaled, site_term), interactions = list()),
  "model_spline_with_depth_with_hab" = list(predictors = c(depth_term_scaled, hab_term, pca_predictors_scaled, site_term), interactions = list()),
  "model_interaction_arch_no_depth_no_hab" = list(predictors = c(pca_predictors_scaled, arch_term, site_term), interactions = interactions_arch),
  "model_interaction_arch_with_depth_no_hab" = list(predictors = c(depth_term_scaled, pca_predictors_scaled, arch_term, site_term), interactions = interactions_arch),
  "model_interaction_arch_no_depth_with_hab" = list(predictors = c(hab_term, pca_predictors_scaled, arch_term, site_term), interactions = interactions_arch),
  "model_interaction_arch_with_depth_with_hab" = list(predictors = c(depth_term_scaled, hab_term, pca_predictors_scaled, arch_term, site_term), interactions = interactions_arch)
)
candidate_models <- candidate_models[!duplicated(candidate_models)]
cat(sprintf("    -> Total de modelos únicos a serem testados: %d\n", length(candidate_models)))


### 3. Configurações OTIMIZADAS -----------------------------------------------
set.seed(42)
RUN_MODE <- "COMPLETO" 
min_obs_in_node <- 5
SAMPLE_FRAC <- 0.7
USE_PARALLEL <- TRUE
N_CORES <- if (USE_PARALLEL) max(1, parallel::detectCores() - 1) else 1
if (RUN_MODE == "TESTE") {
  n.folds <- 3; max_trees <- 500; tree_complexities <- c(1, 2); learning_rates <- c(0.01); bag_fractions <- c(0.75); N_BOOTSTRAP <- 10
} else {
  n.folds <- 5; max_trees <- 3000; tree_complexities <- c(1, 2, 3, 4); learning_rates <- c(0.001, 0.005, 0.01); bag_fractions <- c(0.5, 0.75, 0.9); N_BOOTSTRAP <- 100
}

main_out_dir <- "C:/Users/rbfra/OneDrive/BRT_Abundance_CV_30"
dir.create(main_out_dir, showWarnings = FALSE, recursive = TRUE)


### 4. Funções Utilitárias (SEÇÃO COMPLETA E CORRIGIDA) -----------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

run_linear_comparison <- function(data, response_var, predictor_vars, interactions = NULL, coords_df, random_effect_var = "SITE", p_threshold = 0.05) {
  all_vars_needed <- unique(c(response_var, predictor_vars, random_effect_var, unlist(interactions)))
  df_complete <- data[, intersect(names(data), all_vars_needed)] %>% filter(complete.cases(.))
  if (nrow(df_complete) < 20) { return(list(model_type = "N/A", message = "Dados insuficientes.", model_object_full = NULL, model_object_final = NULL)) }
  main_effects_in_interaction <- unique(unlist(interactions)); main_effects_only <- setdiff(predictor_vars, main_effects_in_interaction)
  parts <- c(); if (length(main_effects_only) > 0) { parts <- c(parts, paste(main_effects_only, collapse = " + ")) }
  if (!is.null(interactions) && length(interactions) > 0) { interaction_strings <- sapply(interactions, function(pair) paste(pair, collapse = " * ")); parts <- c(parts, interaction_strings) }
  formula_fixed_part_full <- paste(parts, collapse = " + ")
  formula_lm_str_full <- paste(response_var, "~", formula_fixed_part_full)
  model_full <- NULL; model_type <- "N/A"
  formula_lm_full <- as.formula(formula_lm_str_full)
  if (random_effect_var %in% names(df_complete) && nlevels(df_complete[[random_effect_var]]) > 1) {
    formula_lmm <- as.formula(paste(formula_lm_str_full, paste0("+ (1 | ", random_effect_var, ")")))
    suppressMessages({ lmm_model <- try(lmer(formula_lmm, data = df_complete, REML = FALSE), silent = TRUE) })
    if (!inherits(lmm_model, "try-error")) {
      if (isSingular(lmm_model)) { model_type <- "LM (RE variance ~0)"; model_full <- lm(formula_lm_full, data = df_complete)
      } else { model_type <- "LMM"; model_full <- lmm_model }
    } else { model_type <- "LM (LMM failed)"; model_full <- lm(formula_lm_full, data = df_complete) }
  } else { model_type <- "LM (No valid RE)"; model_full <- lm(formula_lm_full, data = df_complete) }
  if (is.null(model_full)) { return(list(model_type = "N/A", message = "Falha no ajuste do modelo completo.", model_object_full = NULL, model_object_final = NULL)) }
  model_final <- model_full; significant_terms <- NULL
  try({
    anova_results <- car::Anova(model_full, type = "II")
    p_value_col <- if ("Pr(>F)" %in% names(anova_results)) "Pr(>F)" else "Pr(>Chisq)"
    significant_rows <- anova_results[which(anova_results[[p_value_col]] < p_threshold), ]
    significant_terms <- rownames(significant_rows)
    if (length(significant_terms) > 0) {
      formula_fixed_part_final <- paste(significant_terms, collapse = " + ")
      formula_lm_str_final <- paste(response_var, "~", formula_fixed_part_final)
      if (grepl("LMM", model_type)) {
        formula_lmm_final <- as.formula(paste(formula_lm_str_final, paste0("+ (1 | ", random_effect_var, ")")))
        model_final <- lmer(formula_lmm_final, data = df_complete, REML = FALSE)
      } else { model_final <- lm(as.formula(formula_lm_str_final), data = df_complete) }
    } else {
      if (grepl("LMM", model_type)) {
        model_final <- lmer(as.formula(paste(response_var, "~ 1 + (1 |", random_effect_var, ")")), data = df_complete, REML = FALSE)
      } else { model_final <- lm(as.formula(paste(response_var, "~ 1")), data = df_complete) }
    }
  }, silent = TRUE)
  r2_vals_full <- try(performance::r2(model_full), silent = TRUE); r2_marg_full <- if (inherits(r2_vals_full, "try-error")) NA else r2_vals_full$R2_marginal %||% r2_vals_full$R2_adjusted
  r2_vals_final <- try(performance::r2(model_final), silent = TRUE); r2_marg_final <- if (inherits(r2_vals_final, "try-error")) NA else r2_vals_final$R2_marginal %||% r2_vals_final$R2_adjusted
  full_summary_text <- capture.output(summary(model_full))
  relative_importance <- tryCatch({
    anova_full <- car::Anova(model_full, type = "II"); anova_df <- as.data.frame(anova_full)
    p_col_name <- if ("Pr(>F)" %in% names(anova_df)) "Pr(>F)" else "Pr(>Chisq)"
    importance_stat_col_name <- if ("Sum Sq" %in% names(anova_df)) "Sum Sq" else "Chisq"
    importance_df <- anova_df %>%
      tibble::rownames_to_column("Predictor") %>%
      dplyr::mutate(Importance_Stat = as.numeric(.data[[importance_stat_col_name]]), p_value = as.numeric(.data[[p_col_name]]),
                    Relative.Importance.perc = (Importance_Stat / sum(Importance_Stat, na.rm = TRUE)) * 100) %>%
      rename(Value = Importance_Stat) %>%
      dplyr::select(Predictor, Value, Relative.Importance.perc, p_value) %>%
      dplyr::arrange(desc(Relative.Importance.perc))
    importance_df
  }, error = function(e) { NULL })
  moran_result <- list(message = "Não calculado.")
  return(list(model_type = model_type, r2_marginal_full = r2_marg_full, r2_marginal_final = r2_marg_final, model_object_full = model_full, model_object_final = model_final, significant_terms = significant_terms, moran_test = moran_result, full_summary_text = full_summary_text, relative_importance = relative_importance))
}
drop_zero_var <- function(df) {
  if (is.null(df) || ncol(df) == 0) return(df)
  keep <- vapply(df, function(x) length(unique(x[!is.na(x)])) > 1, logical(1))
  df[, keep, drop = FALSE]
}
grid_search <- function(df, resp, preds) {
  cols_to_select <- c(resp, intersect(preds, names(df))); if (length(cols_to_select) <= 1) return(NULL)
  df_sub <- df[, cols_to_select, drop = FALSE]; df_complete <- df_sub[complete.cases(df_sub), ]; if(nrow(df_complete) < 10) return(NULL)
  preds_in_complete <- intersect(preds, names(df_complete)); if(length(preds_in_complete) == 0) return(NULL)
  df_preds_only <- df_complete[, preds_in_complete, drop = FALSE]; df_preds_ok <- drop_zero_var(df_preds_only); preds_ok <- names(df_preds_ok)
  if (length(preds_ok) == 0) return(NULL)
  df_final <- cbind(df_complete[, resp, drop = FALSE], df_preds_ok); best <- list(deviance = Inf)
  grid <- expand.grid(tc = tree_complexities, lr = learning_rates, bf = bag_fractions)
  for (i in seq_len(nrow(grid))) {
    set.seed(123 + i); m <- try(gbm(as.formula(paste(resp, "~", paste(preds_ok, collapse = "+"))), data = df_final, distribution = "gaussian", n.trees = max_trees, interaction.depth = grid$tc[i], shrinkage = grid$lr[i], bag.fraction = grid$bf[i], n.minobsinnode = min_obs_in_node, cv.folds = n.folds, keep.data = TRUE, verbose = FALSE), silent = TRUE)
    if (inherits(m, "try-error") || length(m$var.names) == 0 || length(m$cv.error) == 0) next
    it <- try(gbm.perf(m, method = "cv", plot.it = FALSE), silent = TRUE)
    if (inherits(it, "try-error") || !is.numeric(it) || length(it)==0 || it < 1 || it > length(m$cv.error) || is.na(it)) { it <- which.min(m$cv.error); if(length(it) == 0 || !is.finite(m$cv.error[it])) { it <- length(m$cv.error); if (it == 0) next }}
    dev <- m$cv.error[it]; if (!is.na(dev) && dev < best$deviance) { best <- list(model = m, iter = it, deviance = dev, params = grid[i, ], data = df_final) }
  }
  if (is.infinite(best$deviance)) return(NULL); return(best)
}
pdp_bootstrap_fixed_refit <- function(model, data, var, optimal_iter, response_name, model_params, n_grid = 50, n_bootstrap = N_BOOTSTRAP, sample_frac = SAMPLE_FRAC, n_cores = N_CORES) {
  grid_vals <- seq(min(data[[var]], na.rm = TRUE), max(data[[var]], na.rm = TRUE), length.out = n_grid); baseline <- mean(predict(model, newdata = data, n.trees = optimal_iter)); pdp_main_raw <- sapply(grid_vals, function(v) { X_tmp <- data; X_tmp[[var]] <- v; mean(predict(model, newdata = X_tmp, n.trees = optimal_iter)) }); pdp_main <- pdp_main_raw - baseline
  bootstrap_iteration <- function(i) {
    set.seed(12345 + i); idx <- sample(seq_len(nrow(data)), size = floor(sample_frac * nrow(data)), replace = TRUE); Xb <- data[idx, , drop = FALSE]; form <- as.formula(paste(response_name, "~", paste(model$var.names, collapse = "+"))); tc <- model_params$tc; lr <- model_params$lr; bf <- model_params$bf
    refit_model <- try(gbm(form, data = Xb, distribution = "gaussian", n.trees = max_trees, interaction.depth = tc, shrinkage = lr, bag.fraction = bf, n.minobsinnode = min_obs_in_node, keep.data = TRUE, verbose = FALSE), silent = TRUE)
    if (inherits(refit_model, "try-error")) return(rep(NA_real_, n_grid))
    it_refit <- optimal_iter; if(it_refit > length(refit_model$train.error) || it_refit < 1) { it_refit <- length(refit_model$train.error) }; if(it_refit < 1) return(rep(NA_real_, n_grid))
    b_baseline <- mean(predict(refit_model, newdata = Xb, n.trees = it_refit)); preds_on_grid <- sapply(grid_vals, function(v) { X_tmp <- Xb; X_tmp[[var]] <- v; mean(predict(refit_model, newdata = X_tmp, n.trees = it_refit)) }); return(preds_on_grid - b_baseline)
  }
  if (USE_PARALLEL && n_cores > 1) {
    cl <- makeCluster(n_cores); on.exit(stopCluster(cl), add = TRUE); clusterEvalQ(cl, { library(gbm) }); clusterExport(cl, varlist = c("data", "var", "optimal_iter", "response_name", "model_params", "max_trees", "min_obs_in_node", "sample_frac", "n_grid", "grid_vals", "bootstrap_iteration"), envir = environment()); boot_curves <- parSapply(cl, 1:n_bootstrap, bootstrap_iteration)
  } else { boot_curves <- sapply(1:n_bootstrap, bootstrap_iteration) }
  if (!is.matrix(boot_curves) || nrow(boot_curves) != n_grid || ncol(boot_curves) != n_bootstrap) { lower <- upper <- rep(NA_real_, n_grid) } else { lower <- apply(boot_curves, 1, quantile, 0.025, na.rm = TRUE); upper <- apply(boot_curves, 1, quantile, 0.975, na.rm = TRUE) }
  df_plot <- data.frame(x = grid_vals, y = pdp_main, lower = lower, upper = upper); valid_y_range <- range(c(df_plot$y, df_plot$lower, df_plot$upper), na.rm = TRUE); y_max_abs <- max(abs(valid_y_range), na.rm = TRUE); y_lims <- if(is.finite(y_max_abs)) c(-y_max_abs * 1.1, y_max_abs * 1.1) else NULL
  p <- ggplot(df_plot, aes(x = x, y = y)) + geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.3) + geom_line(linewidth = 1.0, colour = "blue") + geom_rug(data = data, aes(x = .data[[var]]), inherit.aes = FALSE, sides = "b", alpha = 0.1) + theme_classic(base_size = 12) + coord_cartesian(ylim = y_lims) + labs(x = var, y = "Efeito parcial") + theme(axis.title = ggplot2::element_text(size=11)); return(p)
}
pdp_categorical_refit <- function(model, data, var, optimal_iter, response_name, model_params, n_bootstrap = N_BOOTSTRAP, sample_frac = SAMPLE_FRAC, n_cores = N_CORES, grouping_info_df = NULL, grouping_var_name = "REEF") {
  data[[var]] <- factor(data[[var]]); levels_var <- levels(data[[var]])
  bootstrap_iteration_cat <- function(i, current_level) {
    set.seed(54321 + i + match(current_level, levels_var) * n_bootstrap); idx <- sample(seq_len(nrow(data)), size = floor(sample_frac * nrow(data)), replace = TRUE); Xb <- data[idx, , drop = FALSE]; Xb[[var]] <- factor(Xb[[var]], levels = levels_var); form <- as.formula(paste(response_name, "~", paste(model$var.names, collapse = "+"))); tc <- model_params$tc; lr <- model_params$lr; bf <- model_params$bf
    refit_model <- try(gbm(form, data = Xb, distribution = "gaussian", n.trees = max_trees, interaction.depth = tc, shrinkage = lr, bag.fraction = bf, n.minobsinnode = min_obs_in_node, keep.data = TRUE, verbose = FALSE), silent = TRUE)
    if (inherits(refit_model, "try-error")) return(NA_real_)
    it_refit <- optimal_iter; if(it_refit > length(refit_model$train.error) || it_refit < 1) { it_refit <- length(refit_model$train.error) }; if(it_refit < 1) return(NA_real_)
    b_baseline <- mean(predict(refit_model, newdata = Xb, n.trees = it_refit)); d_boot <- Xb; d_boot[[var]] <- factor(current_level, levels = levels_var); pred_for_level <- mean(predict(refit_model, newdata = d_boot, n.trees = it_refit)); return(pred_for_level - b_baseline)
  }
  pdp_list <- list(); baseline <- mean(predict(model, newdata = data, n.trees = optimal_iter))
  for (l in levels_var) {
    dtmp <- data; dtmp[[var]] <- factor(l, levels = levels_var); pred_mean <- mean(predict(model, newdata = dtmp, n.trees = optimal_iter)) - baseline
    if (USE_PARALLEL && n_cores > 1) {
      cl <- makeCluster(n_cores); clusterEvalQ(cl, { library(gbm) }); clusterExport(cl, varlist = c("data", "var", "optimal_iter", "response_name", "model_params", "max_trees", "min_obs_in_node", "sample_frac", "levels_var", "bootstrap_iteration_cat"), envir = environment()); boot_preds <- parLapply(cl, 1:n_bootstrap, bootstrap_iteration_cat, current_level = l); stopCluster(cl)
    } else { boot_preds <- lapply(1:n_bootstrap, bootstrap_iteration_cat, current_level = l) }
    boot_preds <- unlist(boot_preds); pdp_list[[as.character(l)]] <- data.frame(level = l, mean_pred = pred_mean, lower = quantile(boot_preds, 0.025, na.rm = TRUE), upper = quantile(boot_preds, 0.975, na.rm = TRUE))
  }
  if(length(pdp_list) == 0) return(NULL); pdp_df <- do.call(rbind, pdp_list)
  if (!is.null(grouping_info_df) && grouping_var_name %in% names(grouping_info_df)) {
    grouping_info_df <- rename(grouping_info_df, "group_col" = all_of(grouping_var_name))
    pdp_df <- left_join(pdp_df, grouping_info_df, by = c("level" = var)); pdp_df$group_col[is.na(pdp_df$group_col)] <- "Desconhecido"; pdp_df$group_col <- as.factor(pdp_df$group_col)
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


### 5. LOOP DE MODELAGEM POR CANDIDATO -------------------------------------------
all_model_results <- list()
log_file <- file.path(main_out_dir, "BRT_vs_Bayesian_Analysis_Log.txt")
if (file.exists(log_file)) file.remove(log_file)

model_counter <- 0
for (model_name in names(candidate_models)) {
  model_counter <- model_counter + 1
  model_config    <- candidate_models[[model_name]]
  resp            <- "COBERTURA"
  preds           <- model_config$predictors
  interactions_to_test <- model_config$interactions
  
  cat(sprintf("\n\n#################################################################\n"))
  cat(sprintf("##### [MODELO %d/%d] INICIANDO: %s #####\n", model_counter, length(candidate_models), model_name))
  cat(sprintf("#################################################################\n"))
  cat(sprintf("\n\n===== INICIANDO ANÁLISE PARA MODELO: %s =====\n", model_name), file = log_file, append = TRUE)
  
  model_out_dir <- file.path(main_out_dir, "Outputs", model_name)
  dir.create(model_out_dir, showWarnings = FALSE, recursive = TRUE)
  
  preds_final <- intersect(preds, names(data))
  df_model_subset <- data[, unique(c(resp, preds_final)), drop = FALSE]
  preds_final <- names(drop_zero_var(df_model_subset[, preds_final, drop = FALSE]))
  if(length(preds_final) == 0 && !grepl("null", model_name)){ cat("  AVISO: Nenhum preditor com variância. Pulando.\n"); next }
  
  current_results <- list(model_name = model_name, lm_results = NULL, brt_results = NULL)
  
  # --- ETAPA 1: ANÁLISE LINEAR (LMM/LM) ---
  cat("\n----- 1. Análise Linear (LMM/LM) -----\n")
  lm_preds <- if(grepl("null", model_name)) "1" else preds_final
  lm_interactions <- if(grepl("null", model_name)) NULL else interactions_to_test
  
  linear_results <- run_linear_comparison(data, resp, lm_preds, lm_interactions, coords_df, random_effect_var = "SITE")
  
  if (!is.null(linear_results$model_object_full)){
    cat(paste("  - Tipo de Modelo:", linear_results$model_type, "| R² (Final):", round(linear_results$r2_marginal_final, 3), "\n"))
    current_results$lm_results <- list(r2_final = linear_results$r2_marginal_final)
    lm_out_dir <- file.path(model_out_dir, "Linear_Model_Outputs")
    dir.create(lm_out_dir, showWarnings = FALSE, recursive = TRUE)
    capture.output(summary(linear_results$model_object_final), file = file.path(lm_out_dir, "summary_final_model.txt"))
    if(!is.null(linear_results$relative_importance)) {
      write.csv(linear_results$relative_importance, file.path(lm_out_dir, "importance_anova.csv"), row.names = FALSE)
    }
  } else {
    cat("  - Falha ao ajustar modelo linear.\n")
  }
  
  # --- ETAPA 2: ANÁLISE BRT COMPLETA ---
  cat("\n----- 2. Análise Principal (BRT) -----\n")
  if (length(preds_final) == 0){
    cat("  - BRT pulado (modelo nulo não tem preditores para BRT).\n")
  } else {
    best_brt <- grid_search(data, resp, preds_final)
    if (is.null(best_brt)) {
      cat("  Nenhum modelo BRT válido encontrado.\n")
    } else {
      y_response <- best_brt$data[[resp]]; null_deviance <- var(y_response); model_deviance <- best_brt$deviance
      cv_explained_deviance <- 1 - (model_deviance / null_deviance)
      cat(paste("  --> Melhor BRT. CV Explained Deviance:", paste0(round(cv_explained_deviance * 100, 2), "%\n")))
      current_results$brt_results <- list(cv_explained_deviance = cv_explained_deviance)
      
      brt_out_dir <- file.path(model_out_dir, "BRT_Model_Outputs")
      dir.create(brt_out_dir, showWarnings = FALSE, recursive = TRUE)
      saveRDS(best_brt, file.path(brt_out_dir, "best_brt_object.rds"))
      
      # --- Bloco de Geração de Saídas do BRT ---
      model <- best_brt$model; optimal_iter <- best_brt$iter; df_model <- best_brt$data; model_params <- as.list(best_brt$params)
      
      # Importância
      importance_df <- try(summary(model, n.trees=optimal_iter, plotit=F), silent=T)
      if (!inherits(importance_df, "try-error") && nrow(importance_df) > 0) {
        write.csv(importance_df, file.path(brt_out_dir, "importance_brt.csv"), row.names = FALSE)
      }
      
      # PDPs
      predictors_used <- if(!is.null(importance_df) && nrow(importance_df)>0) importance_df$var else model$var.names
      plot_list <- list()
      for (var in predictors_used) { 
        plot_obj <- NULL
        if (is.numeric(df_model[[var]])) { 
          plot_obj <- try(pdp_bootstrap_fixed_refit(model, df_model, var, optimal_iter, resp, model_params), silent = TRUE) 
        } else if (is.factor(df_model[[var]])) {
          if (var == "SITE" && !is.null(site_to_reef_map)) {
            plot_obj <- try(pdp_categorical_refit(model, df_model, var, optimal_iter, resp, model_params, grouping_info_df = site_to_reef_map, grouping_var_name = "REEF"), silent = TRUE)
          } else {
            plot_obj <- try(pdp_categorical_refit(model, df_model, var, optimal_iter, resp, model_params), silent = TRUE)
          }
        }
        if (inherits(plot_obj, "ggplot")) { plot_obj <- plot_obj + labs(title = var); plot_list[[var]] <- plot_obj } 
      }
      if(length(plot_list) > 0) { 
        ncol_combo <- min(ceiling(sqrt(length(plot_list))), 3); combo_width <- 5 * ncol_combo; combo_height <- 4 * ceiling(length(plot_list) / ncol_combo)
        combined_grob <- gridExtra::arrangeGrob(grobs = plot_list, ncol = ncol_combo, top = grid::textGrob(paste("Efeitos (BRT) para Modelo:", model_name), gp=grid::gpar(fontsize=16)))
        ggsave(filename = file.path(brt_out_dir, "PDPs_Combined.png"), plot = combined_grob, width = combo_width, height = combo_height, dpi = 300, limitsize = FALSE)
      }
      
      # Interações
      if (best_brt$params$tc >= 2 && length(interactions_to_test) > 0) {
        interaction_strengths <- list()
        for (pair in interactions_to_test) {
          if (all(c(pair[1], pair[2]) %in% model$var.names)) {
            strength <- try(gbm::interact.gbm(model, data=df_model, i.var=pair, n.trees=optimal_iter), silent=T)
            if (!inherits(strength, "try-error")) { interaction_strengths[[paste(pair, collapse = ":")]] <- data.frame(var1.names = pair[1], var2.names = pair[2], strength = strength) }
          }
        }
        if (length(interaction_strengths) > 0) {
          full_interaction_list <- do.call(rbind, interaction_strengths) %>% arrange(desc(strength))
          write.csv(full_interaction_list, file.path(brt_out_dir, "interactions_ranked.csv"), row.names=F)
        }
      }
    }
  }
  all_model_results[[model_name]] <- current_results
}

### 6. COMPARAÇÃO FINAL DOS MODELOS ----------------------------------------------
cat("\n\n#################################################################\n")
cat("##### COMPARAÇÃO FINAL DE PERFORMANCE ENTRE TODOS OS MODELOS #####\n")
cat("#################################################################\n")

comparison_summary <- lapply(all_model_results, function(res) {
  data.frame(
    model_name = res$model_name,
    lm_r2_final = if(!is.null(res$lm_results)) round(res$lm_results$r2_final, 4) else NA,
    brt_cv_deviance = if(!is.null(res$brt_results)) round(res$brt_results$cv_explained_deviance, 4) else NA
  )
})
comparison_df <- do.call(rbind, comparison_summary)
comparison_df <- comparison_df %>%
  arrange(desc(brt_cv_deviance), desc(lm_r2_final))

cat("\n--- Tabela de Comparação de Performance ---\n")
print(as.data.frame(comparison_df))

write.csv(comparison_df, file.path(main_out_dir, "FINAL_MODEL_PERFORMANCE_COMPARISON.csv"), row.names = FALSE)
capture.output(print(as.data.frame(comparison_df)), file = file.path(main_out_dir, "FINAL_MODEL_PERFORMANCE_COMPARISON.txt"))

cat("\n\n===== PROCESSAMENTO FINALIZADO. Verifique a tabela de comparação e os diretórios de saída. =====\n")