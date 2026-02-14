### ==================================================================================== ###
### === MASTER VISUALIZATION PIPELINE V7 - YEAR RANDOM EFFECT MODELS (CATEGORICAL)   === ###
### ===                                                                              === ###
### === BASED ON: 04_MASTER_Viz_Pipeline_YEAR_RE_v6_INTEGRATED.R                     === ###
### ===                                                                              === ###
### === V7 IMPROVEMENTS (NEW):                                                        === ###
### === 1. PDPs for CATEGORICAL variables (HABMERGED, ARCH, REEF)                   === ###
### === 2. Categorical × ARCH interaction detection with position_dodge              === ###
### === 3. Combined continuous + categorical plots in single figure                  === ###
### ===                                                                              === ###
### === V6 IMPROVEMENTS (INHERITED):                                                 === ###
### === 1. ARCH interaction detection in Marginal Effects                           === ###
### === 2. HABMERGED display in Forest Plot                                         === ###
### === 3. YEAR SD annotation in all relevant plots                                 === ###
### === 4. Interaction Spotlight for significant ARCH interactions                   === ###
### === 5. Enhanced YEAR RE Diagnostics with significance testing                   === ###
### ==================================================================================== ###

rm(list = ls())
gc()

### 1. PACOTES NECESSARIOS --------------------------------------------------------
libs <- c("brms", "ggplot2", "dplyr", "tidybayes", "ggdist", "bayesplot",
          "patchwork", "cowplot", "loo", "purrr", "tidyverse", "scales")
missing_packages <- libs[!libs %in% installed.packages()[, "Package"]]
if (length(missing_packages > 0)) {
    cat("Instalando pacotes faltantes:", paste(missing_packages, collapse = ", "), "\n")
    install.packages(missing_packages, dependencies = TRUE)
}
invisible(lapply(libs, library, character.only = TRUE))

### 2. TEMA PROFISSIONAL ESTILO NATURE/SCIENCE ------------------------------------
### EXACT COPY FROM ORIGINAL - DO NOT MODIFY
theme_publication <- function(base_size = 14) {
    theme_classic(base_size = base_size) +
        theme(
            text = element_text(family = "sans", color = "black"),
            axis.text = element_text(color = "black", size = 12),
            axis.title = element_text(face = "bold", size = 14),
            axis.line = element_line(color = "black", linewidth = 0.5),
            plot.title = element_text(face = "bold", size = 16, hjust = 0),
            plot.subtitle = element_text(size = 11, color = "gray30"),
            plot.caption = element_text(size = 9, color = "gray50", hjust = 1),
            legend.position = "bottom",
            legend.text = element_text(size = 11),
            legend.title = element_text(face = "bold", size = 12),
            strip.background = element_rect(fill = "gray95", color = NA),
            strip.text = element_text(face = "bold", size = 12),
            panel.grid = element_blank()
        )
}
theme_set(theme_publication())

# Paleta Okabe-Ito (acessível para daltonismo) - EXACT COPY FROM ORIGINAL
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
color_primary <- "#0072B2"
color_secondary <- "#E69F00"
arc_colors <- c("Inner" = "#0072B2", "Outer" = "#E69F00")
arc_colors_full <- c("Inner Arc" = "#0072B2", "Outer Arc" = "#E69F00")

### 3. CONFIGURAÇÃO DOS DIRETÓRIOS -----------------------------------------------
search_directories_YEAR_RE_v7 <- c(
    "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/ZOIB_Abundance_YEAR_RE",
    "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/Gaussian_RGR",
    "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/Gaussian_Health_YEAR_RE",
    "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/JSDM_Dirichlet_YEAR_RE"
)

output_dir_YEAR_RE_v7 <- "C:/Users/rbfra/OneDrive/Bayesian_Figures_YEAR_RE_v7_CATEGORICAL"
dir.create(output_dir_YEAR_RE_v7, showWarnings = FALSE, recursive = TRUE)

cat("\n", rep("=", 80), "\n", sep = "")
cat("MASTER VISUALIZATION PIPELINE V7 - YEAR RANDOM EFFECT MODELS (CATEGORICAL)\n")
cat("V7 NEW: Categorical PDPs (HABMERGED, ARCH, REEF)\n")
cat("V6 Inherited: ARCH interactions, HABMERGED display, YEAR SD annotation\n")
cat(rep("=", 80), "\n", sep = "")
cat("Output directory:", output_dir_YEAR_RE_v7, "\n")

### 4. FUNÇÕES AUXILIARES --------------------------------------------------------

discover_winner_models <- function(search_dirs) {
    cat("\n--- Buscando modelos WINNER (.rds) ---\n")
    all_rds_files <- c()
    for (dir_path in search_dirs) {
        if (dir.exists(dir_path)) {
            files <- list.files(path = dir_path, pattern = "^WINNER_.*\\.rds$",
                               recursive = TRUE, full.names = TRUE, ignore.case = FALSE)
            files <- files[!grepl("_BACKUP\\.rds$", files, ignore.case = TRUE)]
            all_rds_files <- c(all_rds_files, files)
            cat(sprintf("  ✓ %s: %d modelo(s)\n", basename(dirname(dir_path)), length(files)))
        }
    }
    cat(sprintf("\nTotal: %d modelos WINNER\n", length(all_rds_files)))
    return(all_rds_files)
}

detect_model_type_year_re_v7 <- function(model, model_name) {
    if (!inherits(model, "brmsfit")) {
        return(list(type = "UNKNOWN", response_label = "Unknown",
                    uses_arch = FALSE, has_year_re = FALSE))
    }
    fam <- family(model)$family

    # Check for YEAR random effect
    formula_str <- paste(deparse(formula(model)), collapse = " ")
    has_year_re <- grepl("\\(1.*\\|.*YEAR\\)", formula_str) ||
                   grepl("\\(1.*\\|.*REEF.*YEAR\\)", formula_str)

    # Check for ARCH in model
    has_arch <- "ARCH" %in% names(model$data)

    if (grepl("zero_one_inflated_beta", fam, ignore.case = TRUE)) {
        return(list(type = "ZOIB",
                    response_label = ifelse(has_year_re, "Cover (Proportion) + Year RE", "Cover (Proportion)"),
                    uses_arch = has_arch, has_year_re = has_year_re))
    } else if (grepl("gaussian", fam, ignore.case = TRUE)) {
        resp_var <- all.vars(formula(model))[1]
        label_suffix <- ifelse(has_year_re, " + Year RE", "")
        return(list(type = "GAUSSIAN",
                    response_label = paste(resp_var, label_suffix),
                    uses_arch = has_arch, has_year_re = has_year_re))
    } else if (grepl("dirichlet", fam, ignore.case = TRUE)) {
        return(list(type = "DIRICHLET",
                    response_label = ifelse(has_year_re, "Composition + Year RE", "Composition"),
                    uses_arch = has_arch, has_year_re = has_year_re))
    }
    return(list(type = "UNKNOWN", response_label = "Unknown",
                uses_arch = FALSE, has_year_re = FALSE))
}

#' Função auxiliar: limpar nomes de preditores para display
clean_predictor_name <- function(name) {
    name <- gsub("MAGNITUDE", " Magnitude", name)
    name <- gsub("VARIABILITY", " Variability", name)
    name <- gsub("HABMERGED", "Habitat", name)
    name <- gsub("DEPTHM", "Depth", name)
    name <- gsub("HABMERGEDRR_TP", "Habitat: RR_TP", name)
    name <- gsub("HABMERGEDPA", "Habitat: PA", name)
    name <- gsub("_", " ", name)
    name <- trimws(name)
    return(name)
}

#' Função auxiliar: criar template de newdata
create_newdata_template <- function(model_data, n) {
    cols_needed <- names(model_data)
    newdata <- data.frame(matrix(NA, nrow = n, ncol = length(cols_needed)))
    names(newdata) <- cols_needed

    for (col in cols_needed) {
        if (is.numeric(model_data[[col]])) {
            newdata[[col]] <- mean(model_data[[col]], na.rm = TRUE)
        } else {
            newdata[[col]] <- names(which.max(table(model_data[[col]])))[1]
        }
    }

    for (col in cols_needed) {
        if (is.factor(model_data[[col]])) {
            newdata[[col]] <- factor(newdata[[col]], levels = levels(model_data[[col]]))
        }
    }

    return(newdata)
}

### ============================================================================
### V7 NEW: DETECÇÃO DE VARIÁVEIS CATEGÓRICAS
### ============================================================================

#' Detecta variáveis categóricas no modelo
#' @param model_data dados do modelo
#' @param model_formula fórmula do modelo
#' @return lista com variáveis categóricas e seus níveis
detect_categorical_predictors <- function(model_data, model_formula) {
    all_vars <- all.vars(model_formula)
    response_var <- all_vars[1]
    predictors <- setdiff(all_vars, response_var)
    
    categorical_preds <- list()
    
    for (pred in predictors) {
        if (pred %in% names(model_data)) {
            if (is.factor(model_data[[pred]]) || is.character(model_data[[pred]])) {
                levels_vec <- unique(model_data[[pred]])
                levels_vec <- levels_vec[!is.na(levels_vec)]
                
                categorical_preds[[pred]] <- list(
                    levels = levels_vec,
                    n_levels = length(levels_vec)
                )
            }
        }
    }
    
    return(categorical_preds)
}

### ============================================================================
### V7 NEW: PDP PARA VARIÁVEIS CATEGÓRICAS
### ============================================================================

#' Gera PDP para variáveis categóricas
#' @param model modelo brmsfit
#' @param pred nome do preditor categórico
#' @param model_data dados do modelo
#' @param model_info informações do modelo
#' @param ndraws número de draws posteriores
#' @param has_arch_interaction se há interação com ARCH
#' @return objeto ggplot2
generate_categorical_pdp <- function(model, pred, model_data, model_info,
                                     ndraws = 500, has_arch_interaction = FALSE) {
    
    cat(sprintf("      → Gerando PDP categórico para %s...\n", pred))
    
    # Obter níveis da variável
    pred_levels <- unique(model_data[[pred]])
    pred_levels <- pred_levels[!is.na(pred_levels)]
    n_levels <- length(pred_levels)
    
    if (n_levels == 0) {
        cat(sprintf("        ⚠ Nenhum nível encontrado para %s\n", pred))
        return(NULL)
    }
    
    # Verificar se há ARCH para interação
    arch_levels <- if ("ARCH" %in% names(model_data)) {
        unique(model_data$ARCH)
    } else {
        NULL
    }
    
    if (has_arch_interaction && !is.null(arch_levels) && length(arch_levels) >= 2) {
        # PDP com interação ARCH (múltiplos pontos agrupados)
        plot_data_list <- list()
        
        for (arch in arch_levels) {
            for (level in pred_levels) {
                newdata <- create_newdata_template(model_data, 1)
                newdata[[pred]] <- level
                newdata$ARCH <- arch
                
                pred_arr <- tryCatch({
                    posterior_epred(model, newdata = newdata, ndraws = ndraws)
                }, error = function(e) NULL)
                
                if (!is.null(pred_arr)) {
                    plot_data_list[[paste(arch, level, sep = "_")]] <- data.frame(
                        level = as.character(level),
                        estimate = median(pred_arr),
                        lower = quantile(pred_arr, 0.025),
                        upper = quantile(pred_arr, 0.975),
                        ARCH = as.character(arch),
                        stringsAsFactors = FALSE
                    )
                }
            }
        }
        
        if (length(plot_data_list) == 0) {
            cat(sprintf("        ⚠ Falha nas predições para %s × ARCH\n", pred))
            return(NULL)
        }
        
        plot_data <- bind_rows(plot_data_list)
        
        # Normalizar nomes de ARC
        plot_data$ARC <- ifelse(grepl("inner", plot_data$ARCH, ignore.case = TRUE),
                                "Inner Arc", "Outer Arc")
        plot_data$ARC <- factor(plot_data$ARC, levels = c("Inner Arc", "Outer Arc"))
        
        p <- ggplot(plot_data, aes(x = level, y = estimate, color = ARC, group = ARC)) +
            geom_pointrange(aes(ymin = lower, ymax = upper),
                           size = 1, linewidth = 1.2,
                           position = position_dodge(width = 0.5)) +
            scale_color_manual(values = arc_colors_full, name = "ARC") +
            labs(
                title = paste("Effect of", clean_predictor_name(pred), "× ARC"),
                x = clean_predictor_name(pred),
                y = ifelse(model_info$type == "ZOIB", "Predicted Proportion", "Predicted Value"),
                caption = "Point: Median | Bar: 95% CI"
            ) +
            theme_publication() +
            theme(legend.position = "bottom",
                  axis.text.x = element_text(angle = 45, hjust = 1))
        
    } else {
        # PDP simples (sem interação)
        plot_data_list <- list()
        
        for (level in pred_levels) {
            newdata <- create_newdata_template(model_data, 1)
            newdata[[pred]] <- level
            
            pred_arr <- tryCatch({
                posterior_epred(model, newdata = newdata, ndraws = ndraws)
            }, error = function(e) NULL)
            
            if (!is.null(pred_arr)) {
                plot_data_list[[as.character(level)]] <- data.frame(
                    level = as.character(level),
                    estimate = median(pred_arr),
                    lower = quantile(pred_arr, 0.025),
                    upper = quantile(pred_arr, 0.975),
                    stringsAsFactors = FALSE
                )
            }
        }
        
        if (length(plot_data_list) == 0) {
            cat(sprintf("        ⚠ Falha nas predições para %s\n", pred))
            return(NULL)
        }
        
        plot_data <- bind_rows(plot_data_list)
        
        p <- ggplot(plot_data, aes(x = level, y = estimate)) +
            geom_pointrange(aes(ymin = lower, ymax = upper),
                           color = color_primary, size = 1, linewidth = 1.2) +
            labs(
                title = paste("Effect of", clean_predictor_name(pred)),
                x = clean_predictor_name(pred),
                y = ifelse(model_info$type == "ZOIB", "Predicted Proportion", "Predicted Value"),
                caption = "Point: Median | Bar: 95% CI"
            ) +
            theme_publication() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
    }
    
    return(p)
}

### ============================================================================
### V6 IMPROVEMENT 1: FOREST PLOT COM HABMERGED E YEAR SD
### ============================================================================

generate_forest_plot_year_re_v7 <- function(model, model_name, model_info) {
    cat("    -> Gerando Forest Plot V7 (com HABMERGED, YEAR SD, interações)...\n")

    # Para DIRICHLET, usar função JSDM existente
    if (model_info$type == "DIRICHLET") {
        cat("    → Modelo DIRICHLET detectado: usando função JSDM\n")
        p_jsdm <- generate_jsdm_forest_plot_errorbar_CORRECTED(model, model_name)
        if (!is.null(p_jsdm)) {
            return(p_jsdm)
        }
    }

    # Standard forest plot for GAUSSIAN and ZOIB
    all_vars <- get_variables(model)
    fixed_vars <- all_vars[grepl("^b_", all_vars)]
    fixed_vars <- fixed_vars[!grepl("Intercept|sds_|^b_sigma|^b_phi|^b_zoi|^b_coi", fixed_vars)]

    if (length(fixed_vars) == 0) {
        cat("    ⚠ Nenhuma variável fixa encontrada\n")
        return(NULL)
    }

    draws_df <- as_draws_df(model)
    vars_in_draws <- intersect(fixed_vars, names(draws_df))
    if (length(vars_in_draws) == 0) return(NULL)

    draws_long <- draws_df %>%
        dplyr::select(all_of(vars_in_draws)) %>%
        tidyr::pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
        mutate(parameter = gsub("^b_", "", parameter))

    # V7: Adicionar YEAR SD annotation quando aplicável
    subtitle_text <- paste("Model:", gsub("^WINNER_", "", gsub("\\.rds$", "", model_name)))
    if (model_info$has_year_re) {
        year_sd <- tryCatch({
            summary(model)$random$YEAR["sd_(Intercept)", "Estimate"]
        }, error = function(e) NA)
        if (!is.na(year_sd)) {
            subtitle_text <- paste(subtitle_text, sprintf("| Year SD: %.3f", year_sd))
        }
    }

    # V7: Adicionar nota sobre interações ARCH se presentes
    if (model_info$uses_arch) {
        formula_str <- paste(deparse(formula(model)), collapse = " ")
        has_arch_interaction <- grepl(":ARCH", formula_str)
        if (has_arch_interaction) {
            subtitle_text <- paste(subtitle_text, "| ARCH interactions: YES")
        }
    }

    # Criar plot
    if (model_info$type == "GAUSSIAN") {
        p <- ggplot(draws_long, aes(y = reorder(parameter, abs(value)), x = value)) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.7) +
            stat_slab(aes(fill = after_stat(x > 0)), alpha = 0.5, slab_size = 0.3, scale = 0.6) +
            stat_pointinterval(aes(fill = after_stat(x > 0)), point_interval = median_qi,
                              point_size = 1.2, interval_size = 0.8, .width = c(0.89, 0.95)) +
            scale_fill_manual(values = c("TRUE" = color_primary, "FALSE" = "gray60"), guide = "none") +
            labs(title = paste("Fixed Effects -", model_info$type),
                 subtitle = subtitle_text,
                 x = "Effect Estimate", y = NULL,
                 caption = "Points: Median | Bars: 89% and 95% CI") +
            theme_publication()
    } else {
        # ZOIB
        p <- ggplot(draws_long, aes(y = reorder(parameter, abs(value)), x = value, fill = after_stat(x > 0))) +
            stat_halfeye(.width = c(0.89, 0.95), alpha = 0.85, point_interval = median_qi,
                        point_size = 1.5, slab_alpha = 0.6, slab_size = 0.3) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.7) +
            scale_fill_manual(values = c("TRUE" = color_primary, "FALSE" = "gray60"), guide = "none") +
            labs(title = paste("Fixed Effects -", model_info$type),
                 subtitle = subtitle_text,
                 x = "Effect Estimate", y = NULL,
                 caption = "Points: Median | Bars: 89% and 95% CI") +
            theme_publication()
    }
    return(p)
}

### ============================================================================
### V7 NEW: MARGINAL EFFECTS COM VARIÁVEIS CATEGÓRICAS
### ============================================================================

generate_marginal_effects_year_re_v7 <- function(model, model_name, model_info,
                                                resolution = 200, ndraws = 500) {
    cat(sprintf("    -> Gerando PDPs V7 (res=%d, ndraws=%d)...\n", resolution, ndraws))

    if (model_info$type == "DIRICHLET") {
        cat("    ⚠ Marginal effects simplificado para DIRICHLET (usar JSDM functions)\n")
        return(NULL)
    }

    model_data <- model$data
    model_formula <- if (!is.null(model$formula)) model$formula else formula(model)
    if (is.list(model_formula) && !is.null(model_formula$formula)) {
        model_formula <- model_formula$formula
    }

    all_vars <- all.vars(model_formula)
    if (length(all_vars) < 2) return(NULL)

    response_var <- all_vars[1]
    predictors <- setdiff(all_vars, response_var)

    # V7: Detectar interações ARCH na fórmula
    formula_str <- paste(deparse(model_formula), collapse = " ")
    has_arch <- "ARCH" %in% names(model_data)

    cat(sprintf("    → Formula check: has_arch=%s\n", has_arch))

    # ─────────────────────────────────────────────────────────────
    # PARTE 1: Preditores Contínuos (código existente)
    # ─────────────────────────────────────────────────────────────
    continuous_preds <- c()
    for (pred in predictors) {
        if (pred %in% names(model_data) && is.numeric(model_data[[pred]])) {
            continuous_preds <- c(continuous_preds, pred)
        }
    }

    priority <- c("PC1MAGNITUDE", "PC2MAGNITUDE", "PC1VARIABILITY", "PC2VARIABILITY", "DEPTHM")
    continuous_preds <- c(intersect(priority, continuous_preds),
                          setdiff(continuous_preds, priority))

    plot_list <- list()

    for (pred in continuous_preds) {
        cat(sprintf("      → Processando contínuo: %s...\n", pred))

        pred_range <- range(model_data[[pred]], na.rm = TRUE)
        pred_seq <- seq(pred_range[1], pred_range[2], length.out = resolution)

        # V7: Verificar se há interação com ARCH para este preditor
        interaction_pattern <- paste0("(\\b", pred, "\\b.*:ARCH|ARCH.*:", pred, "\\b)")
        has_interaction <- grepl(interaction_pattern, formula_str) && has_arch

        if (has_interaction) {
            cat(sprintf("        ✓ INTERAÇÃO ARCH detectada para %s (duas curvas)\n", pred))

            # Criar newdata para cada ARC
            newdata_inner <- create_newdata_template(model_data, resolution)
            newdata_outer <- newdata_inner
            newdata_inner$ARCH <- model_data$ARCH[1]
            newdata_outer$ARCH <- model_data$ARCH[1]
            if (length(unique(model_data$ARCH)) >= 2) {
                newdata_inner$ARCH <- unique(model_data$ARCH)[1]
                newdata_outer$ARCH <- unique(model_data$ARCH)[2]
            }
            newdata_inner[[pred]] <- pred_seq
            newdata_outer[[pred]] <- pred_seq

            # Predições separadas
            pred_inner <- tryCatch({
                posterior_epred(model, newdata = newdata_inner, ndraws = ndraws)
            }, error = function(e) {
                cat(sprintf("        ⚠ Erro predição Inner: %s\n", substr(e$message, 1, 60)))
                NULL
            })

            pred_outer <- tryCatch({
                posterior_epred(model, newdata = newdata_outer, ndraws = ndraws)
            }, error = function(e) {
                cat(sprintf("        ⚠ Erro predição Outer: %s\n", substr(e$message, 1, 60)))
                NULL
            })

            if (!is.null(pred_inner) && !is.null(pred_outer)) {
                plot_data <- bind_rows(
                    data.frame(
                        x = pred_seq,
                        estimate = apply(pred_inner, 2, median),
                        lower = apply(pred_inner, 2, quantile, 0.025),
                        upper = apply(pred_inner, 2, quantile, 0.975),
                        ARC = unique(newdata_inner$ARCH)
                    ),
                    data.frame(
                        x = pred_seq,
                        estimate = apply(pred_outer, 2, median),
                        lower = apply(pred_outer, 2, quantile, 0.025),
                        upper = apply(pred_outer, 2, quantile, 0.975),
                        ARC = unique(newdata_outer$ARCH)
                    )
                )

                plot_data$ARC <- gsub("inner", "Inner Arc", plot_data$ARC, ignore.case = TRUE)
                plot_data$ARC <- gsub("outer", "Outer Arc", plot_data$ARC, ignore.case = TRUE)
                plot_data$ARC <- factor(plot_data$ARC, levels = c("Inner Arc", "Outer Arc"))

                p <- ggplot(plot_data, aes(x = x, y = estimate, color = ARC, fill = ARC, group = ARC)) +
                    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, color = NA) +
                    geom_line(linewidth = 1.2) +
                    scale_color_manual(values = arc_colors_full, name = "ARC") +
                    scale_fill_manual(values = arc_colors_full, name = "ARC") +
                    labs(
                        title = paste("Effect of", clean_predictor_name(pred), "× ARC"),
                        x = clean_predictor_name(pred),
                        y = ifelse(model_info$type == "ZOIB", "Predicted Proportion", "Predicted Value"),
                        caption = sprintf("Line: Median | Band: 95%% CI | Resolution: %d points", resolution)
                    ) +
                    theme_publication() +
                    theme(legend.position = "bottom")
            } else {
                # Fallback para curva única
                cat(sprintf("        ⚠ Usando curva única (erro em predição separada)\n"))
                newdata <- create_newdata_template(model_data, resolution)
                newdata[[pred]] <- pred_seq
                pred_arr <- tryCatch({
                    posterior_epred(model, newdata = newdata, ndraws = ndraws)
                }, error = function(e) NULL)

                if (!is.null(pred_arr)) {
                    plot_data <- data.frame(
                        x = pred_seq,
                        estimate = apply(pred_arr, 2, median),
                        lower = apply(pred_arr, 2, quantile, 0.025),
                        upper = apply(pred_arr, 2, quantile, 0.975)
                    )
                    p <- ggplot(plot_data, aes(x = x, y = estimate)) +
                        geom_ribbon(aes(ymin = lower, ymax = upper), fill = color_primary, alpha = 0.25) +
                        geom_line(color = color_primary, linewidth = 1.2) +
                        labs(
                            title = paste("Effect of", clean_predictor_name(pred)),
                            x = clean_predictor_name(pred),
                            y = ifelse(model_info$type == "ZOIB", "Predicted Proportion", "Predicted Value")
                        ) +
                        theme_publication()
                } else {
                    p <- NULL
                }
            }

        } else {
            # Uma curva só (sem interação)
            cat(sprintf("        Sem interação ARCH para %s (curva única)\n", pred))

            newdata <- create_newdata_template(model_data, resolution)
            newdata[[pred]] <- pred_seq

            pred_arr <- tryCatch({
                posterior_epred(model, newdata = newdata, ndraws = ndraws)
            }, error = function(e) {
                cat(sprintf("        ⚠ Erro em posterior_epred: %s\n", substr(e$message, 1, 60)))
                NULL
            })

            if (!is.null(pred_arr)) {
                plot_data <- data.frame(
                    x = pred_seq,
                    estimate = apply(pred_arr, 2, median),
                    lower = apply(pred_arr, 2, quantile, 0.025),
                    upper = apply(pred_arr, 2, quantile, 0.975)
                )

                p <- ggplot(plot_data, aes(x = x, y = estimate)) +
                    geom_ribbon(aes(ymin = lower, ymax = upper), fill = color_primary, alpha = 0.25) +
                    geom_line(color = color_primary, linewidth = 1.2) +
                    labs(
                        title = paste("Effect of", clean_predictor_name(pred)),
                        x = clean_predictor_name(pred),
                        y = ifelse(model_info$type == "ZOIB", "Predicted Proportion", "Predicted Value"),
                        caption = sprintf("Line: Median | Band: 95%% CI | Resolution: %d points", resolution)
                    ) +
                    theme_publication()
            } else {
                p <- NULL
            }
        }

        if (!is.null(p)) {
            plot_list[[paste0("cont_", pred)]] <- p
        }
    }

    # ─────────────────────────────────────────────────────────────
    # PARTE 2: Preditores Categóricos (NOVO V7)
    # ─────────────────────────────────────────────────────────────
    categorical_preds <- detect_categorical_predictors(model_data, model_formula)
    
    if (length(categorical_preds) > 0) {
        cat(sprintf("    → Detectados %d preditor(es) categórico(s)\n", length(categorical_preds)))
        
        # Prioridade para HABMERGED
        priority_categorical <- c("HABMERGED", "ARCH", "REEF")
        cat_order <- intersect(priority_categorical, names(categorical_preds))
        cat_order <- c(cat_order, setdiff(names(categorical_preds), priority_categorical))
        
        for (pred in cat_order) {
            cat(sprintf("      → Processando categórico: %s (%d níveis)\n",
                        pred, categorical_preds[[pred]]$n_levels))
            
            # Verificar se há interação com ARCH
            has_arch_int <- grepl(paste0(pred, ".*:ARCH|ARCH.*:", pred), formula_str)
            
            # Não gerar PDP para ARCH isolado se já está em interação
            if (pred == "ARCH" && has_arch_int) {
                cat("        ⚠ ARCH em interação - pulando PDP isolado\n")
                next
            }
            
            p_cat <- generate_categorical_pdp(
                model = model,
                pred = pred,
                model_data = model_data,
                model_info = model_info,
                ndraws = ndraws,
                has_arch_interaction = has_arch_int
            )
            
            if (!is.null(p_cat)) {
                plot_list[[paste0("cat_", pred)]] <- p_cat
            }
        }
    }

    if (length(plot_list) == 0) return(NULL)

    # Combinar todos os plots (contínuos + categóricos)
    combined <- wrap_plots(plot_list, ncol = 2) +
        plot_annotation(
            title = "Partial Dependence Plots (V7: Continuous + Categorical)",
            subtitle = paste("Model:", gsub("^WINNER_", "", gsub("\\.rds$", "", model_name)))
        )

    return(combined)
}

### ============================================================================
### V6 IMPROVEMENT 3: YEAR RE DIAGNOSTICS COM TESTE DE SIGNIFICÂNCIA
### ============================================================================

generate_year_re_diagnostics_v7 <- function(model, model_name) {

    # Check if model has YEAR random effect
    formula_str <- paste(deparse(formula(model)), collapse = " ")
    if (!grepl("YEAR", formula_str, fixed = TRUE)) {
        return(NULL)
    }

    # Extract YEAR random effects
    year_re_list <- ranef(model)$YEAR
    if (is.null(year_re_list)) {
        return(NULL)
    }

    # Handle different structures of ranef output
    if (is.list(year_re_list) && length(year_re_list) > 0) {
        year_re <- year_re_list[[1]]
    } else {
        year_re <- year_re_list
    }

    # Extract standard errors if available
    year_se_list <- tryCatch({
        se(ranef(model)$YEAR)
    }, error = function(e) NULL)

    if (is.list(year_se_list) && length(year_se_list) > 0) {
        year_se <- year_se_list[[1]]
    } else {
        year_se <- year_se_list
    }

    # Prepare data - handle both matrix and vector cases
    if (is.matrix(year_re)) {
        year_df <- data.frame(
            YEAR = rownames(year_re),
            Estimate = year_re[, 1]
        )
        if (!is.null(year_se) && is.matrix(year_se)) {
            year_df$SE <- year_se[, 1]
        } else {
            year_df$SE <- NA
        }
    } else if (is.numeric(year_re)) {
        year_names <- names(year_re)
        if (is.null(year_names) || all(year_names == "")) {
            if ("YEAR" %in% names(model$data)) {
                year_names <- sort(unique(model$data$YEAR))
            } else {
                year_names <- seq_along(year_re)
            }
        }
        year_df <- data.frame(
            YEAR = year_names,
            Estimate = as.numeric(year_re)
        )
        if (!is.null(year_se) && is.numeric(year_se)) {
            year_df$SE <- as.numeric(year_se)
        } else {
            year_df$SE <- NA
        }
    } else {
        return(NULL)
    }

    # Clean YEAR names
    year_df$YEAR <- gsub("^.*\\)|\\(|\\)", "", year_df$YEAR)
    year_df$YEAR <- as.integer(as.character(year_df$YEAR))

    # Remove NA years and sort
    year_df <- year_df[!is.na(year_df$YEAR), ]
    year_df <- year_df[order(year_df$YEAR), ]

    if (nrow(year_df) == 0) return(NULL)

    # V7: Adicionar teste de significância
    if (any(!is.na(year_df$SE))) {
        year_df$Significant <- abs(year_df$Estimate) > (2 * year_df$SE)
        year_df$CI_Lower <- year_df$Estimate - 1.96 * year_df$SE
        year_df$CI_Upper <- year_df$Estimate + 1.96 * year_df$SE
    } else {
        year_df$Significant = FALSE
        year_df$CI_Lower = NA
        year_df$CI_Upper = NA
    }

    # Create plot
    caption_text <- "Points: Median | Line: Trend"
    if (any(!is.na(year_df$SE))) {
        caption_text <- "Points: Median | Error bars: 95% CI | Red = Significant (|Est| > 2*SE)"
    }

    p <- ggplot(year_df, aes(x = YEAR, y = Estimate)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_point(aes(color = Significant, size = Significant)) +
        geom_line(aes(group = 1), color = okabe_ito[1], alpha = 0.5) +
        scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "#0072B2")) +
        scale_size_manual(values = c("TRUE" = 4, "FALSE" = 2)) +
        labs(
            title = "Year Random Effects (V7)",
            subtitle = "Deviation from overall intercept by year",
            x = "Year",
            y = "Random Effect (deviation)",
            caption = caption_text
        ) +
        theme_publication() +
        theme(legend.position = "none")

    # Add error bars if SE is available
    if (any(!is.na(year_df$SE))) {
        p <- p + geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2, linewidth = 0.8)
    }

    return(p)
}

### ============================================================================
### V6 IMPROVEMENT 4: INTERACTION SPOTLIGHT PARA INTERAÇÕES SIGNIFICATIVAS
### ============================================================================

#' Visualização destacada para interações ARCH significativas
#' @param model Modelo brmsfit
#' @param interaction_term Nome da interação (ex: "muMACROALGAEprop_PC1VARIABILITY:ARCHouter")
#' @param species Espécie focal (ex: "MACROALGAEprop")
#' @param predictor Preditor focal (ex: "PC1VARIABILITY")
#' @return Objeto ggplot2 combinado ou NULL
generate_interaction_spotlight_ARCH <- function(model, interaction_term, species, predictor) {

    cat(sprintf("    -> Gerando Interaction Spotlight: %s\n", interaction_term))

    # Verificar se a interação existe no modelo
    model_summary <- summary(model)
    fixed_effects <- model_summary$fixed

    # Verificar se o coeficiente existe e é significativo
    if (!interaction_term %in% rownames(fixed_effects)) {
        cat(sprintf("        ⚠ Interação %s não encontrada no modelo\n", interaction_term))
        return(NULL)
    }

    coef_estimate <- fixed_effects[interaction_term, "Estimate"]
    ci_lower <- fixed_effects[interaction_term, "l-95% CI"]
    ci_upper <- fixed_effects[interaction_term, "u-95% CI"]

    is_significant <- (ci_lower > 0) || (ci_upper < 0)

    cat(sprintf("        → Estimate: %.3f [%.3f, %.3f] %s\n",
                coef_estimate, ci_lower, ci_upper,
                ifelse(is_significant, "SIGNIFICANT", "ns")))

    if (!is_significant) {
        cat("        ⚠ Interação não é significativa (95% CI inclui 0)\n")
    }

    model_data <- model$data

    # Verificar se ARCH existe nos dados
    if (!"ARCH" %in% names(model_data)) {
        cat("        ⚠ Modelo não tem variável ARCH\n")
        return(NULL)
    }

    # Verificar níveis de ARCH
    arch_levels <- unique(model_data$ARCH)
    if (length(arch_levels) < 2) {
        cat("        ⚠ Menos de 2 níveis de ARC nos dados\n")
        return(NULL)
    }

    # Criar grade de predição
    if (!predictor %in% names(model_data)) {
        cat(sprintf("        ⚠ Preditor %s não encontrado nos dados do modelo\n", predictor))
        return(NULL)
    }
    pred_range <- range(model_data[[predictor]], na.rm = TRUE)
    if (any(!is.finite(pred_range))) {
        cat(sprintf("        ⚠ Preditor %s tem valores inválidos\n", predictor))
        return(NULL)
    }
    pred_seq <- seq(pred_range[1], pred_range[2], length.out = 200)

    # Criar newdata para cada nível de ARCH
    plot_data_list <- list()

    for (arc_level in arch_levels) {
        newdata <- create_newdata_template(model_data, length(pred_seq))
        newdata$ARCH <- arc_level
        newdata[[predictor]] <- pred_seq

        pred_arr <- tryCatch({
            posterior_epred(model, newdata = newdata, ndraws = 1000)
        }, error = function(e) {
            cat(sprintf("        ⚠ Erro em posterior_epred para %s: %s\n", arc_level, substr(e$message, 1, 60)))
            NULL
        })

        if (!is.null(pred_arr)) {
            # Para JSDM, extrair espécie específica
            if (model$family$family == "dirichlet") {
                response_cols <- grep("prop$", colnames(pred_arr), value = TRUE)
                sp_col <- grep(species, response_cols, ignore.case = TRUE)
                if (length(sp_col) > 0) {
                    sp_draws <- pred_arr[, , sp_col[1]]
                } else {
                    sp_draws <- pred_arr[, , 1]
                }
            } else {
                sp_draws <- pred_arr
            }

            plot_data_list[[arc_level]] <- data.frame(
                x = pred_seq,
                estimate = apply(sp_draws, 2, median),
                lower = apply(sp_draws, 2, quantile, 0.025),
                upper = apply(sp_draws, 2, quantile, 0.975),
                ARC = arc_level
            )
        }
    }

    if (length(plot_data_list) < 2) {
        cat("        ⚠ Não foi possível gerar predições para ambos os ARC\n")
        return(NULL)
    }

    plot_data <- bind_rows(plot_data_list)

    # Calcular diferença entre níveis de ARC
    if (length(arch_levels) == 2) {
        level1 <- plot_data_list[[arch_levels[1]]]
        level2 <- plot_data_list[[arch_levels[2]]]

        diff_draws <- level1$estimate - level2$estimate

        plot_data$diff_estimate <- diff_draws
        plot_data$diff_lower <- level1$lower - level2$upper
        plot_data$diff_upper <- level1$upper - level2$lower
        plot_data$diff_significant <- (plot_data$diff_lower > 0 & plot_data$diff_upper > 0) |
                                      (plot_data$diff_lower < 0 & plot_data$diff_upper < 0)
    }

    # PLOT 1: Curvas sobrepostas
    p1 <- ggplot(plot_data, aes(x = x, y = estimate, color = ARC, fill = ARC, group = ARC)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, color = NA) +
        geom_line(linewidth = 1.5) +
        scale_color_manual(values = arc_colors_full, name = "ARC") +
        scale_fill_manual(values = arc_colors_full, name = "ARC") +
        labs(
            title = sprintf("Interaction: %s × ARCH", clean_predictor_name(predictor)),
            subtitle = sprintf("%s response | %s (Est: %.3f)",
                            clean_predictor_name(species),
                            ifelse(is_significant, "SIGNIFICANT", "ns"),
                            coef_estimate),
            x = clean_predictor_name(predictor),
            y = "Predicted Proportion",
            caption = "Bands: 95% CI"
        ) +
        theme_publication()

    # PLOT 2: Diferença entre ARC (se houver 2 níveis)
    if (length(arch_levels) == 2 && "diff_estimate" %in% names(plot_data)) {
        plot_data$sig_label <- ifelse(plot_data$diff_significant, "Sig.", "ns")

        p2 <- ggplot(plot_data, aes(x = x, y = diff_estimate)) +
            geom_ribbon(aes(ymin = diff_lower, ymax = diff_upper), alpha = 0.3, fill = color_primary) +
            geom_line(color = color_primary, linewidth = 1.2) +
            geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
            geom_point(aes(color = sig_label), alpha = 0.6) +
            scale_color_manual(values = c("Sig." = "#D55E00", "ns" = "gray70")) +
            labs(
                title = sprintf("Difference (%s - %s)", arch_levels[1], arch_levels[2]),
                subtitle = "Positive = First level > Second level",
                x = clean_predictor_name(predictor),
                y = "Difference in Predicted Proportion",
                caption = "Band = 95% CI | Red = Significant"
            ) +
            theme_publication() +
            theme(legend.position = "none")

        combined <- (p1 / p2) + plot_annotation(
            title = sprintf("Interaction Spotlight V7: %s × ARCH", predictor),
            subtitle = sprintf("Species: %s | Estimate: %.3f [%.3f, %.3f]",
                              clean_predictor_name(species), coef_estimate, ci_lower, ci_upper)
        )
    } else {
        combined <- p1
    }

    return(combined)
}

### ============================================================================
### V7: JSDM-SPECIFIC FUNCTIONS (IMPORTED FROM V6)
### ============================================================================

#' Valida a ordem das espécies no modelo JSDM Dirichlet
validate_jsdm_species_order <- function(model) {
    cat("    -> [VALIDAÇÃO] Verificando ordem das espécies JSDM...\n")

    model_formula <- if (!is.null(model$formula)) model$formula else formula(model)
    if (is.list(model_formula) && !is.null(model_formula$formula)) {
        model_formula <- model_formula$formula
    }
    formula_str <- paste(deparse(model_formula), collapse = " ")

    species_match <- regmatches(formula_str, regexpr("cbind\\(([^)]+)\\)", formula_str))

    if (length(species_match) == 0) {
        cat("    ⚠ Modelo não é JSDM Dirichlet (sem cbind na resposta)\n")
        return(NULL)
    }

    species_str <- gsub("cbind\\(([^)]+)\\)", "\\1", species_match)
    formula_species <- trimws(strsplit(species_str, ",")[[1]])

    model_data <- model$data
    data_species <- grep("_prop$", names(model_data), value = TRUE)

    validation <- list(
        is_jsdm = TRUE,
        formula_species = formula_species,
        data_species = data_species,
        n_formula = length(formula_species),
        n_data = length(data_species),
        match_formula_data = setequal(formula_species, data_species),
        mhispida_in_formula = any(grepl("MUSSISMILIA", formula_species, ignore.case = TRUE)),
        mhispida_in_data = any(grepl("MUSSISMILIA", data_species, ignore.case = TRUE)),
        formula_str = formula_str
    )

    if (validation$mhispida_in_formula) {
        validation$mhispida_idx <- which(
            grepl("MUSSISMILIA", formula_species, ignore.case = TRUE)
        )[1]
    }

    cat(sprintf("      → Fórmula (%d espécies): %s\n",
                validation$n_formula,
                paste(formula_species, collapse = ", ")))
    cat(sprintf("      → Dados (%d espécies): %s\n",
                validation$n_data,
                paste(data_species, collapse = ", ")))

    if (validation$mhispida_in_formula) {
        cat(sprintf("      → M. hispida é a espécie #%d\n", validation$mhispida_idx))
    }

    tryCatch({
        test_data <- model_data[1:2, ]
        test_pred <- posterior_epred(model, newdata = test_data, ndraws = 10)
        validation$test_dims <- dim(test_pred)
        validation$test_success <- TRUE
        cat(sprintf("      → Dimensões: %d draws × %d obs × %d espécies\n",
                    dim(test_pred)[1], dim(test_pred)[2], dim(test_pred)[3]))
    }, error = function(e) {
        cat(sprintf("      ⚠ Erro no teste: %s\n", substr(e$message, 1, 80)))
        validation$test_success <- FALSE
    })

    return(validation)
}

#' Forest Plot JSDM com M. hispida como Referência
generate_jsdm_forest_plot_errorbar_CORRECTED <- function(model, model_name) {
    cat("    -> Gerando Forest Plot JSDM (com M. hispida como referência)...\n")

    model_summary <- summary(model)
    fixed_effects <- model_summary$fixed

    if (is.null(fixed_effects) || nrow(fixed_effects) == 0) {
        cat("    ⚠ Nenhum efeito fixo encontrado\n")
        return(NULL)
    }

    plot_data <- data.frame(
        Parameter = rownames(fixed_effects),
        Estimate = as.numeric(fixed_effects[, "Estimate"]),
        Est.Error = as.numeric(fixed_effects[, "Est.Error"]),
        l95 = as.numeric(fixed_effects[, "l-95% CI"]),
        u95 = as.numeric(fixed_effects[, "u-95% CI"]),
        stringsAsFactors = FALSE
    )

    plot_data <- plot_data[!grepl("Intercept$", plot_data$Parameter), ]

    if (nrow(plot_data) == 0) return(NULL)

    plot_data$species <- gsub("^mu([A-Z]+)prop.*$", "\\1", plot_data$Parameter)
    plot_data$predictor <- gsub("^mu[A-Z]+prop_(.*)$", "\\1", plot_data$Parameter)
    plot_data <- plot_data[plot_data$predictor != plot_data$Parameter, ]

    if (nrow(plot_data) == 0) return(NULL)

    # Adicionar M. hispida como referência
    all_predictors <- unique(plot_data$predictor)
    mussismilia_rows <- data.frame(
        Parameter = paste0("muMUSSISMILIAprop_", all_predictors),
        Estimate = 0,
        Est.Error = 0,
        l95 = 0,
        u95 = 0,
        species = "MUSSISMILIA",
        predictor = all_predictors,
        stringsAsFactors = FALSE
    )

    plot_data <- bind_rows(plot_data, mussismilia_rows)

    species_labels <- c(
        "MUSSISMILIA" = "*M. hispida* (ref)",
        "TURF" = "Turf Algae",
        "CCA" = "Crustose Coralline Algae",
        "CYANO" = "Cyanobacteria",
        "MACROALGAE" = "Macroalgae"
    )

    plot_data$species_clean <- species_labels[plot_data$species]
    plot_data$predictor_clean <- gsub("MAGNITUDE", "Mag.", plot_data$predictor)
    plot_data$predictor_clean <- gsub("VARIABILITY", "Var.", plot_data$predictor_clean)
    plot_data$predictor_clean <- gsub("HABMERGED", "Habitat", plot_data$predictor_clean)
    plot_data$predictor_clean <- gsub("ARCH", "×ARC", plot_data$predictor_clean)

    plot_data$label <- paste0(plot_data$species_clean, ": ", plot_data$predictor_clean)
    plot_data$species_order <- factor(plot_data$species,
                                      levels = c("MUSSISMILIA", "TURF", "CCA", "CYANO", "MACROALGAE"))
    plot_data <- plot_data %>%
        arrange(species_order, abs(Estimate)) %>%
        mutate(label = factor(label, levels = unique(label)))

    plot_data$is_reference <- plot_data$species == "MUSSISMILIA"

    species_colors <- c(
        "*M. hispida* (ref)" = "#D55E00",
        "Turf Algae" = "#009E73",
        "Crustose Coralline Algae" = "#0072B2",
        "Cyanobacteria" = "#CC79A7",
        "Macroalgae" = "#E69F00"
    )

    p <- ggplot(plot_data, aes(x = Estimate, y = label)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.7) +
        geom_errorbarh(aes(xmin = l95, xmax = u95, color = species_clean),
                       height = 0.3, linewidth = 1) +
        geom_point(aes(color = species_clean, size = is_reference, shape = is_reference)) +
        scale_color_manual(values = species_colors, name = "Species") +
        scale_size_manual(values = c("TRUE" = 5, "FALSE" = 3.5), guide = "none") +
        scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16), guide = "none") +
        labs(
            title = "JSDM: Species-Specific Environmental Responses",
            subtitle = paste("Model:", gsub("^WINNER_", "", gsub("\\.rds$", "", model_name))),
            x = "Effect size (log-ratio relative to *M. hispida*)",
            y = NULL,
            caption = "Points: Posterior median | Error bars: 95% CI | Reference at 0"
        ) +
        theme_publication() +
        theme(
            legend.position = "right",
            axis.text.y = element_text(size = 9, family = "mono"),
            plot.subtitle = element_text(size = 10)
        )

    return(p)
}

#' PDPs JSDM: Curvas Separadas por Espécie
generate_jsdm_pdp_community_overlaid_CORRECTED <- function(model, model_name,
                                                            ndraws = 500, resolution = 100) {
    cat("    -> Gerando PDP JSDM: Curvas Separadas (corrigido)...\n")

    validation <- validate_jsdm_species_order(model)
    if (is.null(validation) || !validation$test_success) {
        cat("    ⚠ Validação falhou\n")
        return(NULL)
    }

    model_data <- model$data
    species_codes <- validation$formula_species
    mhispida_idx <- validation$mhispida_idx

    pred <- "PC1MAGNITUDE"
    if (!pred %in% names(model_data)) {
        pred <- names(model_data)[sapply(model_data, is.numeric)][1]
    }

    pred_range <- range(model_data[[pred]], na.rm = TRUE)
    pred_seq <- seq(pred_range[1], pred_range[2], length.out = resolution)

    newdata <- create_newdata_template(model_data, resolution)
    newdata[[pred]] <- pred_seq

    pred_array <- tryCatch({
        posterior_epred(model, newdata = newdata, ndraws = ndraws)
    }, error = function(e) {
        cat(sprintf("    ⚠ Erro: %s\n", substr(e$message, 1, 80)))
        return(NULL)
    })

    if (is.null(pred_array)) return(NULL)

    all_species_data <- list()
    for (sp_idx in 1:min(dim(pred_array)[3], length(species_codes))) {
        sp_code <- gsub("_prop$", "", species_codes[sp_idx])
        sp_draws <- pred_array[, , sp_idx]

        all_species_data[[sp_idx]] <- data.frame(
            x = pred_seq,
            estimate = apply(sp_draws, 2, median),
            lower = apply(sp_draws, 2, quantile, 0.025),
            upper = apply(sp_draws, 2, quantile, 0.975),
            species_code = sp_code,
            is_mhispida = (sp_idx == mhispida_idx)
        )
    }

    plot_data <- bind_rows(all_species_data)

    species_labels <- c(
        "MUSSISMILIA" = "*M. hispida*",
        "TURF" = "Turf Algae",
        "CCA" = "CCA",
        "CYANO" = "Cyanobacteria",
        "MACROALGAE" = "Macroalgae"
    )

    plot_data$species_label <- species_labels[plot_data$species_code]
    plot_data$species_label <- factor(plot_data$species_label, levels = species_labels)

    species_colors <- c(
        "*M. hispida*" = "#D55E00",
        "Turf Algae" = "#009E73",
        "CCA" = "#0072B2",
        "Cyanobacteria" = "#CC79A7",
        "Macroalgae" = "#E69F00"
    )

    p <- ggplot(plot_data, aes(x = x, y = estimate, color = species_label,
                                fill = species_label, group = species_label)) +
        geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
        geom_line(aes(linewidth = is_mhispida)) +
        scale_color_manual(values = species_colors, name = "Species") +
        scale_fill_manual(values = species_colors, name = "Species") +
        scale_linewidth_manual(values = c("TRUE" = 2, "FALSE" = 1), guide = "none") +
        labs(
            title = "JSDM: Species Response Curves",
            subtitle = paste("Response to", pred, "| Thick line = *M. hispida*"),
            x = pred,
            y = "Predicted Proportion",
            caption = "Line: Median | Band: 95% CI"
        ) +
        theme_publication() +
        theme(legend.position = "right")

    return(p)
}

#' PDPs JSDM: Painéis individuais por espécie
generate_jsdm_pdp_per_species_CORRECTED <- function(model, model_name,
                                                     ndraws = 500, resolution = 100) {
    cat("    -> Gerando PDPs JSDM: Per-Species (CORRIGIDO)...\n")

    validation <- validate_jsdm_species_order(model)
    if (is.null(validation)) return(NULL)

    model_data <- model$data
    species_codes <- validation$formula_species

    pred <- "PC1MAGNITUDE"
    if (!pred %in% names(model_data)) {
        numeric_cols <- names(model_data)[sapply(model_data, is.numeric)]
        pred <- numeric_cols[1]
    }

    pred_range <- range(model_data[[pred]], na.rm = TRUE)
    pred_seq <- seq(pred_range[1], pred_range[2], length.out = resolution)

    newdata <- create_newdata_template(model_data, resolution)
    newdata[[pred]] <- pred_seq

    pred_array <- posterior_epred(model, newdata = newdata, ndraws = ndraws)
    n_resp <- dim(pred_array)[3]

    all_data <- list()
    colors <- c("#D55E00", "#009E73", "#0072B2", "#CC79A7", "#E69F00")

    for (sp_idx in 1:min(n_resp, length(species_codes))) {
        sp_draws <- pred_array[, , sp_idx]

        all_data[[sp_idx]] <- data.frame(
            x = pred_seq,
            estimate = apply(sp_draws, 2, median),
            lower = apply(sp_draws, 2, quantile, 0.025),
            upper = apply(sp_draws, 2, quantile, 0.975),
            species = gsub("_prop$", "", species_codes[sp_idx])
        )
    }

    plot_data <- bind_rows(all_data)

    species_labels <- c(
        "MUSSISMILIA" = "*M. hispida*",
        "TURF" = "Turf Algae",
        "CCA" = "CCA",
        "CYANO" = "Cyanobacteria",
        "MACROALGAE" = "Macroalgae"
    )
    plot_data$species_clean <- species_labels[plot_data$species]

    p <- ggplot(plot_data, aes(x = x, y = estimate)) +
        geom_ribbon(aes(ymin = lower, ymax = upper, fill = species_clean), alpha = 0.3) +
        geom_line(aes(color = species_clean), linewidth = 1.2) +
        facet_wrap(~species_clean, scales = "free_y", ncol = 3) +
        scale_color_manual(values = colors[1:length(unique(plot_data$species))]) +
        scale_fill_manual(values = colors[1:length(unique(plot_data$species))]) +
        labs(
            title = paste("Species-Specific Responses to", pred),
            x = pred,
            y = "Predicted Proportion"
        ) +
        theme_publication() +
        theme(legend.position = "none", strip.text = element_text(face = "bold", size = 11))

    return(p)
}

### ============================================================================
### V7 EXPANDED: JSDM MULTI-PREDICTOR VISUALIZATIONS
### ============================================================================

#' Generate expanded JSDM community plot showing ALL predictors
#' @param model brmsfit JSDM model
#' @param model_name Model name for labeling
#' @param ndraws Number of posterior draws
#' @param resolution Number of points for prediction sequence
#' @return ggplot object with all predictors as panels
generate_jsdm_pdp_community_all_predictors <- function(model, model_name,
                                                       ndraws = 500, resolution = 100) {
    cat("    -> Gerando JSDM Community Plot: ALL PREDICTORS (EXPANDED)...\n")

    validation <- validate_jsdm_species_order(model)
    if (is.null(validation) || !validation$test_success) {
        cat("    ⚠ Validação falhou\n")
        return(NULL)
    }

    model_data <- model$data
    species_codes <- validation$formula_species
    mhispida_idx <- validation$mhispida_idx

    # ALL predictors to visualize
    all_predictors <- c("PC1MAGNITUDE", "PC2MAGNITUDE", "PC1VARIABILITY", "PC2VARIABILITY", "DEPTHM")
    available_preds <- intersect(all_predictors, names(model_data))

    cat(sprintf("    → Predictors available: %s\n", paste(available_preds, collapse = ", ")))

    if (length(available_preds) == 0) {
        cat("    ⚠ No predictors found\n")
        return(NULL)
    }

    # Species labels and colors
    species_labels <- c(
        "MUSSISMILIA" = "*M. hispida*",
        "TURF" = "Turf Algae",
        "CCA" = "CCA",
        "CYANO" = "Cyanobacteria",
        "MACROALGAE" = "Macroalgae"
    )

    species_colors <- c(
        "*M. hispida*" = "#D55E00",
        "Turf Algae" = "#009E73",
        "CCA" = "#0072B2",
        "Cyanobacteria" = "#CC79A7",
        "Macroalgae" = "#E69F00"
    )

    # Generate plot for each predictor
    plot_list <- list()

    for (pred in available_preds) {
        cat(sprintf("      → Processing %s...\n", pred))

        pred_range <- range(model_data[[pred]], na.rm = TRUE)
        if (any(!is.finite(pred_range))) {
            cat(sprintf("        ⚠ Invalid range for %s, skipping\n", pred))
            next
        }
        pred_seq <- seq(pred_range[1], pred_range[2], length.out = resolution)

        newdata <- create_newdata_template(model_data, resolution)
        newdata[[pred]] <- pred_seq

        pred_array <- tryCatch({
            posterior_epred(model, newdata = newdata, ndraws = ndraws)
        }, error = function(e) {
            cat(sprintf("        ⚠ Error: %s\n", substr(e$message, 1, 60)))
            return(NULL)
        })

        if (is.null(pred_array)) next

        # Process all species
        all_species_data <- list()
        for (sp_idx in 1:min(dim(pred_array)[3], length(species_codes))) {
            sp_code <- gsub("_prop$", "", species_codes[sp_idx])
            sp_draws <- pred_array[, , sp_idx]

            all_species_data[[sp_idx]] <- data.frame(
                x = pred_seq,
                estimate = apply(sp_draws, 2, median),
                lower = apply(sp_draws, 2, quantile, 0.025),
                upper = apply(sp_draws, 2, quantile, 0.975),
                species_code = sp_code,
                is_mhispida = (sp_idx == mhispida_idx)
            )
        }

        plot_data <- bind_rows(all_species_data)
        plot_data$species_label <- species_labels[plot_data$species_code]
        plot_data$species_label <- factor(plot_data$species_label, levels = species_labels)

        # Create individual panel plot
        p <- ggplot(plot_data, aes(x = x, y = estimate, color = species_label,
                                    fill = species_label, group = species_label)) +
            geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.15, color = NA) +
            geom_line(aes(linewidth = is_mhispida)) +
            scale_color_manual(values = species_colors, name = "Species") +
            scale_fill_manual(values = species_colors, name = "Species") +
            scale_linewidth_manual(values = c("TRUE" = 1.5, "FALSE" = 0.8), guide = "none") +
            labs(
                title = clean_predictor_name(pred),
                x = clean_predictor_name(pred),
                y = "Predicted Proportion"
            ) +
            theme_publication() +
            theme(legend.position = "none",
                  plot.title = element_text(face = "bold", size = 12),
                  axis.title = element_text(size = 10))

        plot_list[[pred]] <- p
    }

    if (length(plot_list) == 0) return(NULL)

    # Arrange in 2x3 or 3x2 grid (5 panels total)
    n_plots <- length(plot_list)
    n_cols <- if (n_plots <= 3) 2 else 3
    n_rows <- ceiling(n_plots / n_cols)

    combined <- wrap_plots(plot_list, ncol = n_cols, nrow = n_rows) +
        plot_annotation(
            title = "JSDM: Species Response Curves - All Predictors",
            subtitle = paste("Model:", gsub("^WINNER_", "", gsub("\\.rds$", "", model_name))),
            theme = theme(plot.title = element_text(face = "bold", size = 16),
                         plot.subtitle = element_text(size = 11))
        )

    return(combined)
}


#' Generate per-species JSDM plots showing ALL predictors
#' @param model brmsfit JSDM model
#' @param model_name Model name for labeling
#' @param target_species Species code (e.g., "MUSSISMILIA", "TURF")
#' @param ndraws Number of posterior draws
#' @param resolution Number of points for prediction sequence
#' @return ggplot object with all predictors as panels for one species
generate_jsdm_pdp_one_species_all_predictors <- function(model, model_name,
                                                         target_species = "MUSSISMILIA",
                                                         ndraws = 500, resolution = 100) {
    cat(sprintf("    -> Gerando JSDM Per-Species Plot: %s - ALL PREDICTORS...\n", target_species))

    validation <- validate_jsdm_species_order(model)
    if (is.null(validation)) return(NULL)

    model_data <- model$data
    species_codes <- validation$formula_species

    # Find target species index
    target_idx <- which(grepl(target_species, species_codes, ignore.case = TRUE))
    if (length(target_idx) == 0) {
        cat(sprintf("    ⚠ Species %s not found\n", target_species))
        return(NULL)
    }
    target_idx <- target_idx[1]

    # ALL predictors to visualize
    all_predictors <- c("PC1MAGNITUDE", "PC2MAGNITUDE", "PC1VARIABILITY", "PC2VARIABILITY", "DEPTHM")
    available_preds <- intersect(all_predictors, names(model_data))

    cat(sprintf("    → Species: %s (index %d)\n", species_codes[target_idx], target_idx))
    cat(sprintf("    → Predictors: %s\n", paste(available_preds, collapse = ", ")))

    # Species color mapping
    species_colors <- c(
        "MUSSISMILIA" = "#D55E00",
        "TURF" = "#009E73",
        "CCA" = "#0072B2",
        "CYANO" = "#CC79A7",
        "MACROALGAE" = "#E69F00"
    )
    species_color <- species_colors[target_species]

    species_labels <- c(
        "MUSSISMILIA" = "*M. hispida*",
        "TURF" = "Turf Algae",
        "CCA" = "Crustose Coralline Algae",
        "CYANO" = "Cyanobacteria",
        "MACROALGAE" = "Macroalgae"
    )
    species_label <- species_labels[target_species]

    # Generate plot for each predictor
    plot_list <- list()

    for (pred in available_preds) {
        pred_range <- range(model_data[[pred]], na.rm = TRUE)
        if (any(!is.finite(pred_range))) next

        pred_seq <- seq(pred_range[1], pred_range[2], length.out = resolution)

        newdata <- create_newdata_template(model_data, resolution)
        newdata[[pred]] <- pred_seq

        pred_array <- tryCatch({
            posterior_epred(model, newdata = newdata, ndraws = ndraws)
        }, error = function(e) NULL)

        if (is.null(pred_array)) next

        sp_draws <- pred_array[, , target_idx]

        plot_data <- data.frame(
            x = pred_seq,
            estimate = apply(sp_draws, 2, median),
            lower = apply(sp_draws, 2, quantile, 0.025),
            upper = apply(sp_draws, 2, quantile, 0.975),
            predictor = pred
        )

        p <- ggplot(plot_data, aes(x = x, y = estimate)) +
            geom_ribbon(aes(ymin = lower, ymax = upper), fill = species_color, alpha = 0.3) +
            geom_line(color = species_color, linewidth = 1.2) +
            labs(
                title = clean_predictor_name(pred),
                x = clean_predictor_name(pred),
                y = "Predicted Proportion"
            ) +
            theme_publication() +
            theme(plot.title = element_text(face = "bold", size = 12),
                  axis.title = element_text(size = 10))

        plot_list[[pred]] <- p
    }

    if (length(plot_list) == 0) return(NULL)

    # Arrange in grid
    n_plots <- length(plot_list)
    n_cols <- if (n_plots <= 3) 2 else 3
    n_rows <- ceiling(n_plots / n_cols)

    combined <- wrap_plots(plot_list, ncol = n_cols, nrow = n_rows) +
        plot_annotation(
            title = sprintf("JSDM: %s Response to All Predictors", species_label),
            subtitle = paste("Model:", gsub("^WINNER_", "", gsub("\\.rds$", "", model_name))),
            theme = theme(plot.title = element_text(face = "bold", size = 16),
                         plot.subtitle = element_text(size = 11))
        )

    return(combined)
}


#' Generate ALL per-species JSDM plots (for each species)
#' @param model brmsfit JSDM model
#' @param model_name Model name for labeling
#' @param output_dir Output directory to save figures
#' @param ndraws Number of posterior draws
#' @param resolution Number of points for prediction sequence
generate_all_jsdm_per_species_plots <- function(model, model_name, output_dir,
                                               ndraws = 500, resolution = 100) {
    cat("    -> Gerando JSDM Per-Species Plots: ALL SPECIES...\n")

    validation <- validate_jsdm_species_order(model)
    if (is.null(validation)) return(NULL)

    species_codes <- validation$formula_species

    # Species to plot
    target_species <- c("MUSSISMILIA", "TURF", "CCA", "CYANO", "MACROALGAE")
    available_species <- intersect(target_species,
                                   gsub("_prop$", "", species_codes))

    cat(sprintf("    → Species to plot: %s\n", paste(available_species, collapse = ", ")))

    results <- list()

    for (sp in available_species) {
        cat(sprintf("      → Generating %s plot...\n", sp))

        p <- generate_jsdm_pdp_one_species_all_predictors(
            model = model,
            model_name = model_name,
            target_species = sp,
            ndraws = ndraws,
            resolution = resolution
        )

        if (!is.null(p)) {
            # Save individual species figure
            fig_name <- sprintf("FIGURE_7%s_JSDM_Per_Species_%s.png",
                               letters[which(target_species == sp)], sp)
            fig_path <- file.path(output_dir, fig_name)

            ggsave(fig_path, p, width = 14, height = 10, dpi = 600, bg = "white")
            cat(sprintf("        ✓ Saved: %s\n", fig_name))

            # Also save PDF
            fig_path_pdf <- sub("\\.png$", ".pdf", fig_path)
            ggsave(fig_path_pdf, p, width = 14, height = 10, bg = "white")

            results[[sp]] <- list(plot = p, path = fig_path)
        }
    }

    return(results)
}

### ============================================================================
### V7: PPC PLOTS
### ============================================================================

generate_ppc_plots_year_re_v7 <- function(model, model_name, model_info) {
    cat("    -> PPC...\n")

    if (model_info$type == "DIRICHLET") {
        cat("    ⚠ PPC não disponível para modelos DIRICHLET\n")
        return(NULL)
    }

    ppc_dens <- pp_check(model, ndraws = 100, type = "dens_overlay") +
        scale_color_manual(values = c("y" = "black", "yrep" = color_primary)) +
        labs(title = "Density Validation", subtitle = "Black: Observed | Blue: Simulated") +
        theme_publication() + theme(legend.position = "none")

    if (model_info$type %in% c("ZOIB", "BETA")) {
        ppc_zeros <- tryCatch({
            pp_check(model, ndraws = 100, type = "stat", stat = function(y) mean(y == 0)) +
                labs(title = "Zero Validation", subtitle = "Black line should be under blue distribution") +
                theme_publication() + theme(legend.position = "none")
        }, error = function(e) NULL)

        if (!is.null(ppc_zeros)) {
            return((ppc_dens | ppc_zeros) +
                   plot_annotation(title = "Posterior Predictive Checks"))
        }
    }

    ppc_ecdf <- pp_check(model, ndraws = 100, type = "ecdf_overlay") +
        scale_color_manual(values = c("y" = "black", "yrep" = color_primary)) +
        labs(title = "Cumulative ECDF") +
        theme_publication() + theme(legend.position = "none")

    return((ppc_dens | ppc_ecdf) + plot_annotation(title = "Posterior Predictive Checks"))
}

### ============================================================================
### V7: DIAGNOSTIC PLOTS
### ============================================================================

generate_diagnostics_year_re_v7 <- function(model, model_name) {
    cat("    -> Diagnósticos...\n")
    all_vars <- variables(model)
    params_to_trace <- head(all_vars[grepl("^b_", all_vars)], 4)
    if (length(params_to_trace) < 2) params_to_trace <- head(all_vars, 4)

    p_trace <- mcmc_trace(model, pars = params_to_trace) + theme_publication()
    p_acf <- mcmc_acf(model, pars = params_to_trace) + theme_publication()

    max_rhat <- max(summary(model)$fixed[, "Rhat"], na.rm = TRUE)

    # V7: Adicionar YEAR Rhat se disponível
    if ("YEAR" %in% names(summary(model)$random)) {
        year_rhat <- max(summary(model)$random$YEAR[, "Rhat"], na.rm = TRUE)
        subtitle <- sprintf("Max R-hat: Fixed=%.3f | YEAR=%.3f %s",
                             max_rhat, year_rhat,
                             ifelse(max(max_rhat, year_rhat) < 1.01, "✓", "✗"))
    } else {
        subtitle <- sprintf("Max R-hat: %.3f %s", max_rhat,
                             ifelse(max_rhat < 1.01, "✓", "✗"))
    }

    return((p_trace / p_acf) +
           plot_annotation(title = "MCMC Convergence Diagnostics",
                          subtitle = subtitle))
}

### ============================================================================
### V7: PIPELINE PRINCIPAL
### ============================================================================

process_model_year_re_v7 <- function(rds_file, output_base_dir) {
    model <- readRDS(rds_file)
    model_name <- gsub("\\.rds$", "", basename(rds_file))

    model_output_dir <- file.path(output_base_dir, model_name)
    if (!dir.exists(model_output_dir)) {
        dir.create(model_output_dir, recursive = TRUE)
    }

    cat(sprintf("\n🔍 PROCESSING V7: %s\n", model_name))

    model_info <- detect_model_type_year_re_v7(model, model_name)

    # Diagnosticar estrutura
    cat(sprintf("\n  Tipo: %s | Has Year RE: %s | Uses ARCH: %s\n",
                model_info$type, model_info$has_year_re, model_info$uses_arch))

    # ─────────────────────────────────────────────────────────────
    # FIGURA 1: Forest Plot V7
    # ─────────────────────────────────────────────────────────────
    cat("\n  [1/7] Forest Plot V7 (com HABMERGED, YEAR SD)...\n")
    p_forest <- generate_forest_plot_year_re_v7(model, model_name, model_info)
    if (!is.null(p_forest)) {
        ggsave(file.path(model_output_dir, "FIGURE_1_Forest_Plot.png"),
               p_forest, width = 10, height = 6, dpi = 600, bg = "white")
        ggsave(file.path(model_output_dir, "FIGURE_1_Forest_Plot.pdf"),
               p_forest, width = 10, height = 6, bg = "white")
        cat("    ✓ FIGURE_1 salvo\n")
    }

    # ─────────────────────────────────────────────────────────────
    # FIGURA 2: Marginal Effects V7 (COM CATEGÓRICOS!)
    # ─────────────────────────────────────────────────────────────
    cat("\n  [2/7] Marginal Effects V7 (Continuous + Categorical)...\n")
    if (model_info$type %in% c("GAUSSIAN", "ZOIB")) {
        p_effects <- generate_marginal_effects_year_re_v7(model, model_name, model_info)
        if (!is.null(p_effects)) {
            ggsave(file.path(model_output_dir, "FIGURE_2_Marginal_Effects.png"),
                   p_effects, width = 14, height = 12, dpi = 600, bg = "white")
            ggsave(file.path(model_output_dir, "FIGURE_2_Marginal_Effects.pdf"),
                   p_effects, width = 14, height = 12, bg = "white")
            cat("    ✓ FIGURE_2 salvo (com categóricos!)\n")
        }
    }

    # ─────────────────────────────────────────────────────────────
    # FIGURA 3: PPC
    # ─────────────────────────────────────────────────────────────
    cat("\n  [3/7] PPC...\n")
    p_ppc <- generate_ppc_plots_year_re_v7(model, model_name, model_info)
    if (!is.null(p_ppc)) {
        ggsave(file.path(model_output_dir, "FIGURE_3_PPC.png"),
               p_ppc, width = 12, height = 10, dpi = 600, bg = "white")
        ggsave(file.path(model_output_dir, "FIGURE_3_PPC.pdf"),
               p_ppc, width = 12, height = 10, bg = "white")
        cat("    ✓ FIGURE_3 salvo\n")
    }

    # ─────────────────────────────────────────────────────────────
    # FIGURA 4: Diagnostics
    # ─────────────────────────────────────────────────────────────
    cat("\n  [4/7] Diagnostics...\n")
    p_diag <- generate_diagnostics_year_re_v7(model, model_name)
    if (!is.null(p_diag)) {
        ggsave(file.path(model_output_dir, "FIGURE_4_Diagnostics.png"),
               p_diag, width = 14, height = 10, dpi = 600, bg = "white")
        ggsave(file.path(model_output_dir, "FIGURE_4_Diagnostics.pdf"),
               p_diag, width = 14, height = 10, bg = "white")
        cat("    ✓ FIGURE_4 salvo\n")
    }

    # ─────────────────────────────────────────────────────────────
    # FIGURA 5-8: JSDM-specific (se aplicável)
    # ─────────────────────────────────────────────────────────────
    if (model_info$type == "DIRICHLET") {
        cat("\n  [5/7] JSDM: Forest Plot com M. hispida...\n")
        p_jsdm_forest <- generate_jsdm_forest_plot_errorbar_CORRECTED(model, model_name)
        if (!is.null(p_jsdm_forest)) {
            ggsave(file.path(model_output_dir, "FIGURE_5_JSDM_Forest.png"),
                   p_jsdm_forest, width = 12, height = 8, dpi = 600, bg = "white")
            ggsave(file.path(model_output_dir, "FIGURE_5_JSDM_Forest.pdf"),
                   p_jsdm_forest, width = 12, height = 8, bg = "white")
            cat("    ✓ FIGURE_5 salvo\n")
        }

        cat("\n  [6/7] JSDM: Community Overlaid (PC1MAGNITUDE only - original)...\n")
        p_jsdm_comm <- generate_jsdm_pdp_community_overlaid_CORRECTED(model, model_name)
        if (!is.null(p_jsdm_comm)) {
            ggsave(file.path(model_output_dir, "FIGURE_6_JSDM_Community.png"),
                   p_jsdm_comm, width = 14, height = 10, dpi = 600, bg = "white")
            cat("    ✓ FIGURE_6 salvo (original)\n")
        }

        # ── NEW: Expanded Community Plot with ALL predictors ────────────────
        cat("\n  [6b/7] JSDM: Community Expanded (ALL PREDICTORS)...\n")
        p_jsdm_comm_exp <- generate_jsdm_pdp_community_all_predictors(model, model_name)
        if (!is.null(p_jsdm_comm_exp)) {
            ggsave(file.path(model_output_dir, "FIGURE_6_JSDM_Community_Expanded.png"),
                   p_jsdm_comm_exp, width = 16, height = 12, dpi = 600, bg = "white")
            ggsave(file.path(model_output_dir, "FIGURE_6_JSDM_Community_Expanded.pdf"),
                   p_jsdm_comm_exp, width = 16, height = 12, bg = "white")
            cat("    ✓ FIGURE_6_EXPANDED salvo (all predictors)\n")
        }

        cat("\n  [7/7] JSDM: Per Species (original - PC1MAGNITUDE only)...\n")
        p_jsdm_species <- generate_jsdm_pdp_per_species_CORRECTED(model, model_name)
        if (!is.null(p_jsdm_species)) {
            ggsave(file.path(model_output_dir, "FIGURE_7_JSDM_Per_Species.png"),
                   p_jsdm_species, width = 14, height = 12, dpi = 600, bg = "white")
            cat("    ✓ FIGURE_7 salvo (original)\n")
        }

        # ── NEW: Per-Species Plots with ALL predictors ────────────────────────
        cat("\n  [7b/7] JSDM: Per-Species Expanded (ALL PREDICTORS)...\n")
        per_species_results <- generate_all_jsdm_per_species_plots(
            model = model,
            model_name = model_name,
            output_dir = model_output_dir
        )
        cat(sprintf("    ✓ Generated %d per-species expanded figures\n",
                    length(per_species_results)))
    }

    # ─────────────────────────────────────────────────────────────
    # FIGURA EXTRA: YEAR RE Diagnostics
    # ─────────────────────────────────────────────────────────────
    if (model_info$has_year_re) {
        cat("\n  [EXTRA] YEAR RE Diagnostics V7...\n")
        p_year <- generate_year_re_diagnostics_v7(model, model_name)
        if (!is.null(p_year)) {
            fig_num <- ifelse(model_info$type == "DIRICHET", 8, 5)
            ggsave(file.path(model_output_dir, sprintf("FIGURE_%d_YEAR_RE_Diagnostics.png", fig_num)),
                   p_year, width = 8, height = 6, dpi = 600, bg = "white")
            ggsave(file.path(model_output_dir, sprintf("FIGURE_%d_YEAR_RE_Diagnostics.pdf", fig_num)),
                   p_year, width = 8, height = 6, bg = "white")
            cat(sprintf("    ✓ FIGURA_%d salvo (YEAR RE)\n", fig_num))
        }
    }

    # ─────────────────────────────────────────────────────────────
    # V7 EXTRA: Interaction Spotlight para interações significativas
    # ─────────────────────────────────────────────────────────────
    if (model_info$uses_arch) {
        # Verificar interações significativas no sumário
        model_summary <- summary(model)
        int_terms <- rownames(model_summary$fixed)
        arch_ints <- int_terms[grepl(":ARCH", int_terms)]

        if (length(arch_ints) > 0) {
            cat(sprintf("\n  [EXTRA] %d interações ARCH encontradas. Gerando Spotlight...\n", length(arch_ints)))

            # Para cada interação significativa, gerar spotlight
            for (int_term in arch_ints) {
                tryCatch({
                    # Extrair espécie e preditor
                    # Padrão: muMACROALGAEprop_PC1VARIABILITY:ARCHouter
                    if (grepl("^mu", int_term)) {
                        parts <- strsplit(int_term, ":")[[1]]
                        species_part <- parts[1]  # muMACROALGAEprop
                        predictor_part <- if ("ARCHouter" %in% parts) {
                            gsub(":ARCHouter", "", parts[1])  # PC1VARIABILITY
                        } else {
                            gsub(":ARCHinner", "", parts[1])
                        }

                        # Limpar nomes
                        species_clean <- gsub("prop$", "", species_part)
                        predictor_clean <- predictor_part

                        p_spot <- generate_interaction_spotlight_ARCH(
                            model = model,
                            interaction_term = int_term,
                            species = species_clean,
                            predictor = predictor_clean
                        )

                        if (!is.null(p_spot)) {
                            fig_name <- gsub("[^a-zA-Z0-9]", "_", int_term)
                            ggsave(file.path(model_output_dir, paste0("FIGURE_SPOTLIGHT_", fig_name, ".png")),
                                   p_spot, width = 12, height = 10, dpi = 600, bg = "white")
                            cat(sprintf("    ✓ Spotlight salvo: %s\n", fig_name))
                        }
                    }
                }, error = function(e) {
                    cat(sprintf("    ⚠ Erro ao gerar spotlight para %s: %s\n", int_term, substr(e$message, 1, 60)))
                })
            }
        }
    }

    cat(sprintf("\n✅ MODELO CONCLUÍDO: %s\n", model_name))
}

### ============================================================================
### V7: EXECUÇÃO PRINCIPAL
### ============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("🚀 INICIANDO PIPELINE V7 - YEAR RANDOM EFFECT MODELS (CATEGORICAL)\n")
cat("V7 NEW:\n")
cat("  ✓ PDPs for CATEGORICAL variables (HABMERGED, ARCH, REEF)\n")
cat("  ✓ Categorical × ARCH interaction detection with position_dodge\n")
cat("  ✓ Combined continuous + categorical plots in single figure\n")
cat("V6 Inherited:\n")
cat("  ✓ ARCH interaction detection in Marginal Effects\n")
cat("  ✓ HABMERGED display in Forest Plot\n")
cat("  ✓ YEAR SD annotation in all relevant plots\n")
cat("  ✓ Interaction Spotlight for significant ARCH interactions\n")
cat("  ✓ Enhanced YEAR RE Diagnostics with significance testing\n")
cat(rep("=", 80), "\n", sep = "")

winner_files <- discover_winner_models(search_directories_YEAR_RE_v7)

processed_info <- list()
for (rds_file in winner_files) {
    result <- process_model_year_re_v7(rds_file, output_dir_YEAR_RE_v7)
    if (!is.null(result)) processed_info[[length(processed_info) + 1]] <- result
}

cat("\n", rep("=", 80), "\n", sep = "")
cat("📋 RELATÓRIO FINAL - PIPELINE V7 CATEGORICAL\n")
cat(rep("=", 80), "\n", sep = "")

cat(sprintf("\n✓ Modelos processados: %d\n", length(processed_info)))
cat(sprintf("📁 Saída: %s\n", output_dir_YEAR_RE_v7))

cat("\nFiguras geradas (V7):\n")
cat("  [ZOIB/GAUSSIAN MODELS:]\n")
cat("  - FIGURE_1_Forest_Plot.png/pdf (com YEAR SD annotation)\n")
cat("  - FIGURE_2_Marginal_Effects.png/pdf (V7: Continuous + Categorical!)\n")
cat("  - FIGURE_3_PPC.png/pdf\n")
cat("  - FIGURE_4_Diagnostics.png/pdf\n")
cat("  - FIGURE_5_YEAR_RE_Diagnostics.png/pdf (se aplicável)\n")
cat("\n  [JSDM/DIRICHLET MODELS:]\n")
cat("  - FIGURE_1_Forest_Plot.png/pdf (com M. hispida como referência)\n")
cat("  - FIGURE_5_JSDM_Forest.png/pdf\n")
cat("  - FIGURE_6_JSDM_Community.png/pdf (original: PC1MAGNITUDE only)\n")
cat("  - FIGURE_6_JSDM_Community_Expanded.png/pdf (NEW: ALL 5 predictors!)\n")
cat("  - FIGURE_7_JSDM_Per_Species.png/pdf (original: PC1MAGNITUDE only)\n")
cat("  - FIGURE_7A_JSDM_Per_Species_MUSSISMILIA.png/pdf (NEW: all predictors)\n")
cat("  - FIGURE_7B_JSDM_Per_Species_TURF.png/pdf (NEW: all predictors)\n")
cat("  - FIGURE_7C_JSDM_Per_Species_CCA.png/pdf (NEW: all predictors)\n")
cat("  - FIGURE_7D_JSDM_Per_Species_CYANO.png/pdf (NEW: all predictors)\n")
cat("  - FIGURE_7E_JSDM_Per_Species_MACROALGAE.png/pdf (NEW: all predictors)\n")
cat("  - FIGURE_8_YEAR_RE_Diagnostics.png/pdf (se aplicável)\n")
cat("\n  [EXTRA: Interaction Spotlight]\n")
cat("  - FIGURE_SPOTLIGHT_*.png/pdf (para interações ARCH significativas)\n")

cat("\n🎉 PIPELINE V7 CONCLUÍDO!\n\n")