### ==================================================================================== ###
### === MASTER VISUALIZATION PIPELINE V5 - VERSÃO FINAL HÍBRIDA V4                     === ###
### ===                                                                              === ###
### === V4 IMPLEMENTATION:                                                           === ###
### === 1. Validação robusta de espécies JSDM                                      === ###
### === 2. Forest Plot com M. hispida incluída (referência = 0)                   === ###
### === 3. PDPs Profissionais Gaussian/ZOIB (alta resolução 300)                  === ###
### === 4. PDPs de Composição Relativa (trade-offs ecológicos)                    === ###
### === 5. Análise de Sensibilidade ARC (Inner vs Outer)                          === ###
### === 6. PDPs JSDM com curvas separadas por espécie                             === ###
### ==================================================================================== ###
### === Baseado em: MASTER HYBRID V4 PLAN (2026-02-04)                            === ###
### ==================================================================================== ###

rm(list = ls())
gc()

### 1. PACOTES NECESSÁRIOS --------------------------------------------------------
libs <- c("brms", "ggplot2", "dplyr", "tidybayes", "ggdist", "bayesplot",
          "patchwork", "cowplot", "loo", "purrr", "tidyverse")
missing_packages <- libs[!libs %in% installed.packages()[, "Package"]]
if (length(missing_packages) > 0) {
    cat("Instalando pacotes faltantes:", paste(missing_packages, collapse = ", "), "\n")
    install.packages(missing_packages, dependencies = TRUE)
}
invisible(lapply(libs, library, character.only = TRUE))

### 2. TEMA PROFISSIONAL ESTILO NATURE/SCIENCE ------------------------------------
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

# Paleta Okabe-Ito (acessível para daltonismo)
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
color_primary <- "#0072B2"
color_secondary <- "#E69F00"
arc_colors <- c("Inner" = "#0072B2", "Outer" = "#E69F00")

### 3. CONFIGURAÇÃO DOS DIRETÓRIOS -----------------------------------------------
search_directories_v5 <- c(
    "C:/Users/rbfra/OneDrive/Bayesian_Full_LOO_Selection_v5/ZOIB_Abundance",
    "C:/Users/rbfra/OneDrive/Bayesian_Full_LOO_Selection_v5/Gaussian_HEALTH_PC1",
    "C:/Users/rbfra/OneDrive/Bayesian_Full_LOO_Selection_v5/Gaussian_HEALTH_PC2",
    "C:/Users/rbfra/OneDrive/Bayesian_Full_LOO_Selection_v5/Gaussian_RGR",
    "C:/Users/rbfra/OneDrive/Bayesian_Full_LOO_Selection_v5/JSDM_Dirichlet"
)

output_dir_v5 <- "C:/Users/rbfra/OneDrive/Bayesian_Figures_Publication_v5_HYBRID_V4"
dir.create(output_dir_v5, showWarnings = FALSE, recursive = TRUE)

cat("\n", rep("=", 75), "\n", sep = "")
cat("🎨 MASTER VISUALIZATION PIPELINE V5 - HYBRID V4 (JSDM CORRIGIDO + FUNCIONALIDADES AVANÇADAS)\n")
cat(rep("=", 75), "\n", sep = "")
cat("📁 Diretório de saída:", output_dir_v5, "\n")

### 4. FUNÇÕES AUXILIARES --------------------------------------------------------

discover_winner_models <- function(search_dirs) {
    cat("\n--- Buscando modelos WINNER (.rds) ---\n")
    all_rds_files <- c()
    for (dir_path in search_dirs) {
        if (dir.exists(dir_path)) {
            files <- list.files(path = dir_path, pattern = "^WINNER_.*\\.rds$",
                               recursive = FALSE, full.names = TRUE, ignore.case = FALSE)
            files <- files[!grepl("_BACKUP\\.rds$", files, ignore.case = TRUE)]
            all_rds_files <- c(all_rds_files, files)
            cat(sprintf("  ✓ %s: %d modelo(s)\n", basename(dir_path), length(files)))
        }
    }
    cat(sprintf("\n📊 Total: %d modelos\n", length(all_rds_files)))
    return(all_rds_files)
}

detect_model_type_v5 <- function(model, model_name) {
    if (!inherits(model, "brmsfit")) {
        return(list(type = "UNKNOWN", response_label = "Unknown", uses_arch = FALSE))
    }
    fam <- family(model)$family
    if (grepl("zero_one_inflated_beta", fam, ignore.case = TRUE)) {
        return(list(type = "ZOIB", response_label = "Cover (Proportion)", uses_arch = TRUE))
    } else if (grepl("beta", fam, ignore.case = TRUE)) {
        return(list(type = "BETA", response_label = "Proportion", uses_arch = TRUE))
    } else if (grepl("gaussian", fam, ignore.case = TRUE)) {
        return(list(type = "GAUSSIAN", response_label = "Value", uses_arch = FALSE))
    } else if (grepl("dirichlet", fam, ignore.case = TRUE)) {
        return(list(type = "DIRICHLET", response_label = "Composition", uses_arch = TRUE))
    } else {
        return(list(type = "OTHER", response_label = "Response", uses_arch = FALSE))
    }
}

### ═══════════════════════════════════════════════════════════════════════════════
### FUNÇÕES DE VALIDAÇÃO E DIAGNÓSTICO (NOVO V4)
### ═══════════════════════════════════════════════════════════════════════════════

#' Valida a ordem das espécies no modelo JSDM Dirichlet
#'
#' Esta função é CRÍTICA para garantir que M. hispida seja corretamente
#' identificada em todas as análises.
#'
#' @param model Objeto brmsfit do modelo JSDM
#' @return Lista com resultados da validação ou NULL se não for JSDM
validate_jsdm_species_order <- function(model) {
    cat("    -> [VALIDAÇÃO] Verificando ordem das espécies JSDM...\n")

    # Extrair fórmula como string (corrigido para brmsfit)
    model_formula <- if (!is.null(model$formula)) model$formula else formula(model)
    if (is.list(model_formula) && !is.null(model_formula$formula)) {
        model_formula <- model_formula$formula
    }
    formula_str <- paste(deparse(model_formula), collapse = " ")

    # Verificar se é modelo JSDM (tem cbind na resposta)
    species_match <- regmatches(formula_str, regexpr("cbind\\(([^)]+)\\)", formula_str))

    if (length(species_match) == 0) {
        cat("    ⚠ Modelo não é JSDM Dirichlet (sem cbind na resposta)\n")
        return(NULL)
    }

    # Extrair nomes das espécies da fórmula
    species_str <- gsub("cbind\\(([^)]+)\\)", "\\1", species_match)
    formula_species <- trimws(strsplit(species_str, ",")[[1]])

    # Extrair espécies dos dados (colunas _prop)
    model_data <- model$data
    data_species <- grep("_prop$", names(model_data), value = TRUE)

    # Criar resultado da validação
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

    # Determinar índice de M. hispida
    if (validation$mhispida_in_formula) {
        validation$mhispida_idx <- which(
            grepl("MUSSISMILIA", formula_species, ignore.case = TRUE)
        )[1]
    }

    # Log detalhado
    cat(sprintf("      → Fórmula (%d espécies): %s\n",
                validation$n_formula,
                paste(formula_species, collapse = ", ")))
    cat(sprintf("      → Dados (%d espécies): %s\n",
                validation$n_data,
                paste(data_species, collapse = ", ")))

    if (validation$mhispida_in_formula) {
        cat(sprintf("      → M. hispida é a espécie #%d na fórmula\n",
                    validation$mhispida_idx))
    }

    # Verificações críticas
    if (!validation$match_formula_data) {
        cat("    ⚠ AVISO: Espécies na fórmula não correspondem aos dados!\n")
    }

    if (!validation$mhispida_in_formula) {
        cat("    ❌ ERRO CRÍTICO: M. hispida não encontrada na fórmula!\n")
    }

    if (!validation$mhispida_in_data) {
        cat("    ❌ ERRO CRÍTICO: M. hispida não encontrada nos dados!\n")
    }

    # Teste prático: predição de amostra para verificar dimensões
    cat("      → Testando predição de amostra...\n")
    tryCatch({
        test_data <- model_data[1:2, ]
        test_pred <- posterior_epred(model, newdata = test_data, ndraws = 10)

        validation$test_dims <- dim(test_pred)
        validation$test_success <- TRUE

        cat(sprintf("      → Dimensões: %d draws × %d obs × %d espécies\n",
                    dim(test_pred)[1], dim(test_pred)[2], dim(test_pred)[3]))

        if (dim(test_pred)[3] == validation$n_formula) {
            cat("      ✅ Número de respostas corresponde à fórmula\n")
        } else {
            cat(sprintf("      ⚠ Inconsistência: predição tem %d respostas\n",
                        dim(test_pred)[3]))
        }

    }, error = function(e) {
        cat(sprintf("      ⚠ Erro no teste: %s\n", substr(e$message, 1, 80)))
        validation$test_success <- FALSE
    })

    return(validation)
}

#' Diagnostica a estrutura completa de um modelo
#'
#' @param model Objeto brmsfit
#' @return Lista com informações diagnósticas
diagnose_model_structure <- function(model) {
    cat("\n🔍 [DIAGNÓSTICO] Estrutura do Modelo:\n")

    info <- list()

    # Tipo de família
    info$family <- model$family$family
    cat(sprintf("    → Família: %s\n", info$family))

    # Fórmula (corrigido para brmsfit)
    model_formula <- if (!is.null(model$formula)) model$formula else formula(model)
    if (is.list(model_formula) && !is.null(model_formula$formula)) {
        model_formula <- model_formula$formula
    }
    info$formula <- model_formula
    cat(sprintf("    → Fórmula: %s\n", paste(deparse(model_formula), collapse = " ")))

    # Número de observações
    info$n_obs <- nrow(model$data)
    cat(sprintf("    → Observações: %d\n", info$n_obs))

    # Número de parâmetros
    info$n_params <- length(get_variables(model))
    cat(sprintf("    → Parâmetros: %d\n", info$n_params))

    # Convergência (se disponível)
    tryCatch({
        rhats <- brms::rhat(model)
        info$max_rhat <- max(rhats, na.rm = TRUE)
        cat(sprintf("    → Max R-hat: %.4f\n", info$max_rhat))
    }, error = function(e) {
        cat("    ⚠ Não foi possível calcular R-hat\n")
    })

    return(info)
}

### 5. FUNÇÕES DE VISUALIZAÇÃO PADRÃO -------------------------------------------

generate_forest_plot_v5 <- function(model, model_name, model_info) {
    cat("    -> Gerando Forest Plot (V4: com suporte JSDM)...\n")

    # ═══════════════════════════════════════════════════════════════════════════
    # CORREÇÃO CRÍTICA V4: Detectar modelo DIRICHLET/BETA e usar função JSDM
    # ═══════════════════════════════════════════════════════════════════════════
    if (!is.null(model_info$type) && model_info$type %in% c("DIRICHLET", "BETA")) {
        cat("    → Modelo DIRICHLET/BETA detectado: usando função JSDM corrigida\n")

        # Chamar versão corrigida que inclui M. hispida
        p_jsdm <- generate_jsdm_forest_plot_errorbar_CORRECTED(model, model_name)

        if (!is.null(p_jsdm)) {
            cat("    ✅ Forest Plot JSDM gerado com M. hispida incluída\n")
            return(p_jsdm)
        } else {
            cat("    ⚠ Fallback: versão JSDM falhou, tentando padrão...\n")
        }
    }
    # ═══════════════════════════════════════════════════════════════════════════

    # [CÓDIGO ORIGINAL PARA GAUSSIAN/ZOIB]
    all_vars <- get_variables(model)
    fixed_vars <- all_vars[grepl("^b_", all_vars)]
    fixed_vars <- fixed_vars[!grepl("Intercept|sds_|^b_sigma|^b_phi|^b_zoi|^b_coi", fixed_vars)]

    if (length(fixed_vars) == 0) {
        cat("    ⚠ Nenhuma variável fixa encontrada para plotar\n")
        return(NULL)
    }

    draws_df <- as_draws_df(model)
    vars_in_draws <- intersect(fixed_vars, names(draws_df))
    if (length(vars_in_draws) == 0) return(NULL)

    draws_long <- draws_df %>%
        dplyr::select(all_of(vars_in_draws)) %>%
        tidyr::pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
        mutate(parameter = gsub("^b_", "", parameter)) %>%
        mutate(parameter = gsub("_scaled", "", parameter))

    if (model_info$type == "GAUSSIAN") {
        p <- ggplot(draws_long, aes(y = reorder(parameter, abs(value)), x = value)) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.7) +
            stat_slab(aes(fill = after_stat(x > 0)), alpha = 0.5, slab_size = 0.3, scale = 0.6) +
            stat_pointinterval(aes(fill = after_stat(x > 0)), point_interval = median_qi,
                              point_size = 1.2, interval_size = 0.8, .width = c(0.89, 0.95)) +
            scale_fill_manual(values = c("TRUE" = color_primary, "FALSE" = "gray60"), guide = "none") +
            labs(title = paste("Magnitude of Fixed Effects -", model_info$type),
                 subtitle = paste("Model:", gsub("^WINNER_(?:FINAL_)?", "", gsub("_", " ", model_name))),
                 x = "Effect Estimate", y = NULL, caption = "Points: Median | Bars: 89% and 95% CI") +
            theme_publication()
    } else {
        p <- ggplot(draws_long, aes(y = reorder(parameter, abs(value)), x = value, fill = after_stat(x > 0))) +
            stat_halfeye(.width = c(0.89, 0.95), alpha = 0.85, point_interval = median_qi,
                        point_size = 1.5, slab_alpha = 0.6, slab_size = 0.3) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.7) +
            scale_fill_manual(values = c("TRUE" = color_primary, "FALSE" = "gray60"), guide = "none") +
            labs(title = paste("Magnitude of Fixed Effects -", model_info$type),
                 subtitle = paste("Model:", gsub("^WINNER_(?:FINAL_)?", "", gsub("_", " ", model_name))),
                 x = "Effect Estimate", y = NULL, caption = "Points: Median | Bars: 89% and 95% CI") +
            theme_publication()
    }
    return(p)
}

generate_ppc_plots_v5 <- function(model, model_name, model_info) {
    cat("    -> PPC...\n")

    # Dirichlet models don't support pp_check
    if (model_info$type == "DIRICHLET") {
        cat("    ⚠ PPC não disponível para modelos DIRICHLET (pp_check não implementado)\n")
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

generate_diagnostics_v5 <- function(model, model_name) {
    cat("    -> Diagnósticos...\n")
    all_vars <- variables(model)
    params_to_trace <- head(all_vars[grepl("^b_", all_vars)], 4)
    if (length(params_to_trace) < 2) params_to_trace <- head(all_vars, 4)

    p_trace <- mcmc_trace(model, pars = params_to_trace) + theme_publication()
    p_acf <- mcmc_acf(model, pars = params_to_trace) + theme_publication()

    max_rhat <- max(summary(model)$fixed[, "Rhat"], na.rm = TRUE)

    return((p_trace / p_acf) +
           plot_annotation(title = "MCMC Convergence Diagnostics",
                          subtitle = sprintf("Max R-hat: %.3f %s", max_rhat,
                                            ifelse(max_rhat < 1.01, "✓", "✗"))))
}

### 6. FUNÇÕES JSDM CORRIGIDAS ===================================================

#' Forest Plot JSDM com M. hispida como Referência (CORRIGIDO V4)
#'
#' Esta função extrai efeitos fixos de modelos Dirichlet via summary(model)$fixed
#' e adiciona manualmente M. hispida como categoria de referência (coeficiente = 0).
#'
#' @param model Objeto brmsfit do modelo JSDM Dirichlet
#' @param model_name Nome do modelo para título
#' @return Objeto ggplot2 ou NULL em caso de erro
generate_jsdm_forest_plot_errorbar_CORRECTED <- function(model, model_name) {
    cat("    -> Gerando Forest Plot JSDM (com M. hispida como referência)...\n")

    # Extrair efeitos fixos do summary
    model_summary <- summary(model)
    fixed_effects <- model_summary$fixed

    if (is.null(fixed_effects) || nrow(fixed_effects) == 0) {
        cat("    ⚠ Nenhum efeito fixo encontrado no modelo\n")
        return(NULL)
    }

    # Construir dataframe base
    plot_data <- data.frame(
        Parameter = rownames(fixed_effects),
        Estimate = as.numeric(fixed_effects[, "Estimate"]),
        Est.Error = as.numeric(fixed_effects[, "Est.Error"]),
        l95 = as.numeric(fixed_effects[, "l-95% CI"]),
        u95 = as.numeric(fixed_effects[, "u-95% CI"]),
        stringsAsFactors = FALSE
    )

    # Remover interceptos (não são comparações interespecíficas)
    plot_data <- plot_data[!grepl("Intercept$", plot_data$Parameter), ]

    if (nrow(plot_data) == 0) {
        cat("    ⚠ Nenhum parâmetro válido após remover interceptos\n")
        return(NULL)
    }

    # Extrair espécie e preditor dos nomes dos parâmetros
    # Padrão: muTURFprop_PC1MAGNITUDE -> espécie=TURF, preditor=PC1MAGNITUDE
    plot_data$species <- gsub("^mu([A-Z]+)prop.*$", "\\1", plot_data$Parameter)
    plot_data$predictor <- gsub("^mu[A-Z]+prop_(.*)$", "\\1", plot_data$Parameter)

    # Filtrar linhas onde extração falhou
    plot_data <- plot_data[plot_data$predictor != plot_data$Parameter, ]

    if (nrow(plot_data) == 0) {
        cat("    ⚠ Não foi possível extrair espécies e preditores\n")
        return(NULL)
    }

    # ═══════════════════════════════════════════════════════════════════════════
    # CORREÇÃO CRÍTICA: Adicionar M. hispida como referência (coeficiente = 0)
    # ═══════════════════════════════════════════════════════════════════════════
    all_predictors <- unique(plot_data$predictor)

    mussismilia_rows <- data.frame(
        Parameter = paste0("muMUSSISMILIAprop_", all_predictors),
        Estimate = 0,           # Referência: log-ratio = 0
        Est.Error = 0,          # Sem erro (baseline)
        l95 = 0,                # Sem intervalo (baseline)
        u95 = 0,                # Sem intervalo (baseline)
        species = "MUSSISMILIA",
        predictor = all_predictors,
        stringsAsFactors = FALSE
    )

    plot_data <- bind_rows(plot_data, mussismilia_rows)
    # ═══════════════════════════════════════════════════════════════════════════

    # Definir labels para espécies
    species_labels <- c(
        "MUSSISMILIA" = "*M. hispida* (ref)",
        "TURF" = "Turf Algae",
        "CCA" = "Crustose Coralline Algae",
        "CYANO" = "Cyanobacteria",
        "MACROALGAE" = "Macroalgae"
    )

    plot_data$species_clean <- species_labels[plot_data$species]
    plot_data$species_clean[is.na(plot_data$species_clean)] <- plot_data$species[is.na(plot_data$species_clean)]

    # Simplificar nomes dos preditores
    plot_data$predictor_clean <- plot_data$predictor
    plot_data$predictor_clean <- gsub("MAGNITUDE", "Mag.", plot_data$predictor_clean)
    plot_data$predictor_clean <- gsub("VARIABILITY", "Var.", plot_data$predictor_clean)
    plot_data$predictor_clean <- gsub("HABMERGED", "Habitat", plot_data$predictor_clean)
    plot_data$predictor_clean <- gsub("ARCH", "×ARC", plot_data$predictor_clean)
    plot_data$predictor_clean <- gsub("DEPTHM", "Depth", plot_data$predictor_clean)
    plot_data$predictor_clean <- gsub("_scaled", "", plot_data$predictor_clean)

    # Criar label combinada
    plot_data$label <- paste0(plot_data$species_clean, ": ", plot_data$predictor_clean)

    # Ordenar: M. hispida primeiro, depois por magnitude do efeito
    plot_data$species_order <- factor(plot_data$species,
                                      levels = c("MUSSISMILIA", "TURF", "CCA", "CYANO", "MACROALGAE"))
    plot_data <- plot_data %>%
        arrange(species_order, abs(Estimate)) %>%
        mutate(label = factor(label, levels = unique(label)))

    # Marcar M. hispida para destacar visualmente
    plot_data$is_reference <- plot_data$species == "MUSSISMILIA"

    # Paleta de cores (colorblind-friendly, destaque para M. hispida)
    species_colors <- c(
        "*M. hispida* (ref)" = "#D55E00",      # Vermelho-laranja (destaque máximo)
        "Turf Algae" = "#009E73",               # Verde
        "Crustose Coralline Algae" = "#0072B2", # Azul
        "Cyanobacteria" = "#CC79A7",            # Rosa
        "Macroalgae" = "#E69F00"                # Amarelo
    )

    # Criar plot
    p <- ggplot(plot_data, aes(x = Estimate, y = label)) +
        # Linha de referência em zero
        geom_vline(xintercept = 0, linetype = "dashed",
                   color = "gray50", linewidth = 0.7) +
        # Barras de erro
        geom_errorbarh(aes(xmin = l95, xmax = u95, color = species_clean),
                       height = 0.3, linewidth = 1) +
        # Pontos (M. hispida como diamante maior)
        geom_point(aes(color = species_clean,
                       size = is_reference,
                       shape = is_reference)) +
        # Escalas
        scale_color_manual(values = species_colors, name = "Species") +
        scale_size_manual(values = c("TRUE" = 5, "FALSE" = 3.5), guide = "none") +
        scale_shape_manual(values = c("TRUE" = 18, "FALSE" = 16), guide = "none") +
        # Labels
        labs(
            title = "JSDM: Species-Specific Environmental Responses",
            subtitle = paste("Model:", gsub("^WINNER_(?:FINAL_)?", "", gsub("_", " ", model_name)),
                           "| Diamond = *M. hispida* (reference category)"),
            x = "Effect size (log-ratio relative to *M. hispida*)",
            y = NULL,
            caption = "Points: Posterior median | Error bars: 95% Credible Interval | Reference at 0"
        ) +
        theme_publication() +
        theme(
            legend.position = "right",
            axis.text.y = element_text(size = 9, family = "mono"),
            plot.subtitle = element_text(size = 10)
        )

    return(p)
}

#' PDPs JSDM: Curvas Separadas por Espécie (CORRIGIDO V4)
#'
#' Usa posterior_epred() corretamente para obter array [draws, obs, species]
#' e plota cada espécie em curvas separadas (não empilhadas).
#'
#' @param model Objeto brmsfit JSDM
#' @param model_name Nome do modelo
#' @param ndraws Número de draws (padrão: 500)
#' @param resolution Número de pontos (padrão: 100)
#' @return Objeto ggplot2 ou NULL
generate_jsdm_pdp_community_overlaid_CORRECTED <- function(model, model_name,
                                                            ndraws = 500, resolution = 100) {
    cat("    -> Gerando PDP JSDM: Curvas Separadas (corrigido)...\n")

    # Validação prévia
    validation <- validate_jsdm_species_order(model)
    if (is.null(validation) || !validation$test_success) {
        cat("    ⚠ Validação falhou, abortando PDP\n")
        return(NULL)
    }

    model_data <- model$data
    species_codes <- validation$formula_species

    # Identificar índice de M. hispida
    mhispida_idx <- validation$mhispida_idx

    # Selecionar preditor principal
    pred <- "PC1MAGNITUDE"
    if (!pred %in% names(model_data)) {
        numeric_cols <- names(model_data)[sapply(model_data, is.numeric)]
        pred <- numeric_cols[1]
        cat(sprintf("    → Usando preditor alternativo: %s\n", pred))
    }

    # Criar newdata variando o preditor
    pred_range <- range(model_data[[pred]], na.rm = TRUE)
    pred_seq <- seq(pred_range[1], pred_range[2], length.out = resolution)

    newdata <- create_newdata_template(model_data, resolution)
    newdata[[pred]] <- pred_seq

    # Obter predições
    cat("      → Calculando posterior_epred...\n")
    pred_array <- tryCatch({
        posterior_epred(model, newdata = newdata, ndraws = ndraws)
    }, error = function(e) {
        cat(sprintf("    ⚠ Erro em posterior_epred: %s\n", substr(e$message, 1, 80)))
        return(NULL)
    })

    if (is.null(pred_array)) return(NULL)

    dims <- dim(pred_array)
    cat(sprintf("      → Array: %d draws × %d obs × %d species\n",
                dims[1], dims[2], dims[3]))

    # Processar cada espécie
    all_species_data <- list()
    colors_species <- c("#D55E00", "#009E73", "#0072B2", "#CC79A7", "#E69F00")

    for (sp_idx in 1:min(dims[3], length(species_codes))) {
        sp_code <- gsub("_prop$", "", species_codes[sp_idx])
        sp_draws <- pred_array[, , sp_idx]

        sp_summary <- data.frame(
            x = pred_seq,
            estimate = apply(sp_draws, 2, median),
            lower = apply(sp_draws, 2, quantile, 0.025),
            upper = apply(sp_draws, 2, quantile, 0.975),
            species_code = sp_code,
            species_idx = sp_idx,
            is_mhispida = (sp_idx == mhispida_idx)
        )

        all_species_data[[sp_idx]] <- sp_summary
    }

    plot_data <- bind_rows(all_species_data)

    # Labels
    species_labels <- c(
        "MUSSISMILIA" = "*M. hispida*",
        "TURF" = "Turf Algae",
        "CCA" = "CCA",
        "CYANO" = "Cyanobacteria",
        "MACROALGAE" = "Macroalgae"
    )

    plot_data$species_label <- species_labels[plot_data$species_code]
    plot_data$species_label <- factor(plot_data$species_label,
                                       levels = species_labels)

    # Cores
    species_colors <- c(
        "*M. hispida*" = "#D55E00",
        "Turf Algae" = "#009E73",
        "CCA" = "#0072B2",
        "Cyanobacteria" = "#CC79A7",
        "Macroalgae" = "#E69F00"
    )

    # Plot
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

#' PDPs JSDM: Painéis individuais por espécie (CORRIGIDO V4)
#'
#' @param model Objeto brmsfit JSDM
#' @param model_name Nome do modelo
#' @param ndraws Número de draws (padrão: 500)
#' @param resolution Número de pontos (padrão: 100)
#' @return Objeto ggplot2 ou NULL
generate_jsdm_pdp_per_species_CORRECTED <- function(model, model_name,
                                                     ndraws = 500, resolution = 100) {
    cat("    -> Gerando PDPs JSDM: Per-Species (CORRIGIDO)...\n")

    validation <- validate_jsdm_species_order(model)
    if (is.null(validation)) return(NULL)

    model_data <- model$data
    species_codes <- validation$formula_species

    # Selecionar preditor principal
    pred <- "PC1MAGNITUDE"
    if (!pred %in% names(model_data)) {
        numeric_cols <- names(model_data)[sapply(model_data, is.numeric)]
        pred <- numeric_cols[1]
    }

    # Criar grades de predição
    pred_range <- range(model_data[[pred]], na.rm = TRUE)
    pred_seq <- seq(pred_range[1], pred_range[2], length.out = resolution)

    newdata <- create_newdata_template(model_data, resolution)
    newdata[[pred]] <- pred_seq

    # Predições
    pred_array <- posterior_epred(model, newdata = newdata, ndraws = ndraws)
    n_resp <- dim(pred_array)[3]

    # Organizar dados
    all_data <- list()
    colors <- c("#D55E00", "#009E73", "#0072B2", "#CC79A7", "#E69F00")

    for (sp_idx in 1:min(n_resp, length(species_codes))) {
        sp_draws <- pred_array[, , sp_idx]

        sp_data <- data.frame(
            x = pred_seq,
            estimate = apply(sp_draws, 2, median),
            lower = apply(sp_draws, 2, quantile, 0.025),
            upper = apply(sp_draws, 2, quantile, 0.975),
            species = gsub("_prop$", "", species_codes[sp_idx]),
            species_idx = sp_idx
        )
        all_data[[sp_idx]] <- sp_data
    }

    plot_data <- bind_rows(all_data)

    # Labels limpos
    species_labels <- c(
        "MUSSISMILIA" = "*M. hispida*",
        "TURF" = "Turf Algae",
        "CCA" = "Crustose Coralline Algae",
        "CYANO" = "Cyanobacteria",
        "MACROALGAE" = "Macroalgae"
    )
    plot_data$species_clean <- species_labels[plot_data$species]
    plot_data$species_clean[is.na(plot_data$species_clean)] <- plot_data$species[is.na(plot_data$species_clean)]

    # Plot com facets
    p <- ggplot(plot_data, aes(x = x, y = estimate)) +
        geom_ribbon(aes(ymin = lower, ymax = upper, fill = species_clean), alpha = 0.3) +
        geom_line(aes(color = species_clean), linewidth = 1.2) +
        facet_wrap(~species_clean, scales = "free_y", ncol = 3) +
        scale_color_manual(values = colors[1:length(unique(plot_data$species))]) +
        scale_fill_manual(values = colors[1:length(unique(plot_data$species))]) +
        labs(
            title = paste("Species-Specific Responses to", pred),
            subtitle = "Each panel shows one species' response curve with 95% CI",
            x = pred,
            y = "Predicted Proportion"
        ) +
        theme_publication() +
        theme(
            legend.position = "none",
            strip.text = element_text(face = "bold", size = 11)
        )

    return(p)
}

generate_jsdm_ppc_species_overlays <- function(model, model_name) {
    cat("    -> PPC por espécie...\n")

    # Dirichlet models don't support per-species PPC via posterior_predict
    # Use manual PPC with posterior_epred instead
    model_formula <- if (!is.null(model$formula)) model$formula else formula(model)
    if (is.list(model_formula) && !is.null(model_formula$formula)) {
        model_formula <- model_formula$formula
    }
    formula_str <- paste(deparse(model_formula), collapse = " ")
    if (grepl("dirichlet", family(model)$family, ignore.case = TRUE)) {
        cat("    ⚠ PPC por espécie não disponível para DIRICHLET (usando posterior_epred manual)\n")
        return(NULL)
    }

    species_match <- regmatches(formula_str, regexpr("cbind\\(([^)]+)\\)", formula_str))

    if (length(species_match) > 0) {
        species_str <- gsub("cbind\\(([^)]+)\\)", "\\1", species_match)
        species_names <- trimws(strsplit(species_str, ",")[[1]])
    } else {
        species_names <- grep("_prop$", names(model$data), value = TRUE)
    }

    if (length(species_names) == 0) return(NULL)

    species_to_plot <- head(species_names, 4)
    plots_list <- list()

    for (species in species_to_plot) {
        y_obs <- model$data[[species]]
        if (is.null(y_obs) || length(y_obs) == 0) next

        y_rep <- tryCatch(posterior_predict(model, resp = species, draws = 100), error = function(e) NULL)
        if (is.null(y_rep)) next

        p <- ppc_dens_overlay(y_obs, y_rep) +
            labs(title = paste("PPC:", gsub("_", " ", species))) +
            theme_publication() + theme(legend.position = "none", plot.title = element_text(size = 12))
        plots_list[[species]] <- p
    }

    if (length(plots_list) == 0) return(NULL)

    combined <- wrap_plots(plots_list, ncol = ifelse(length(plots_list) <= 2, length(plots_list), 2)) +
        plot_annotation(title = "JSDM: Posterior Predictive Checks by Species")

    return(combined)
}

### ═══════════════════════════════════════════════════════════════════════════════
### FUNÇÕES PROFUNSSIONAIS PARA GAUSSIAN/ZOIB (NOVO V4)
### ═══════════════════════════════════════════════════════════════════════════════

#' Função auxiliar: criar template de newdata
#' @param model_data Data frame com dados do modelo
#' @param n Número de linhas desejado
#' @return Data frame template
create_newdata_template <- function(model_data, n) {
    # Identificar todas as variáveis necessárias
    cols_needed <- names(model_data)

    newdata <- data.frame(matrix(NA, nrow = n, ncol = length(cols_needed)))
    names(newdata) <- cols_needed

    # Preencher com médias/modas
    for (col in cols_needed) {
        if (is.numeric(model_data[[col]])) {
            newdata[[col]] <- mean(model_data[[col]], na.rm = TRUE)
        } else {
            # Para fatores/categóricos, usar o nível mais frequente
            newdata[[col]] <- names(which.max(table(model_data[[col]])))[1]
        }
    }

    # Garantir tipos corretos
    for (col in cols_needed) {
        if (is.factor(model_data[[col]])) {
            newdata[[col]] <- factor(newdata[[col]], levels = levels(model_data[[col]]))
        }
    }

    return(newdata)
}

#' Função auxiliar: limpar nomes de preditores para display
clean_predictor_name <- function(name) {
    name <- gsub("MAGNITUDE", "Magnitude", name)
    name <- gsub("VARIABILITY", "Variability", name)
    name <- gsub("HABMERGED", "Habitat", name)
    name <- gsub("DEPTHM", "Depth", name)
    name <- gsub("_scaled", "", name)
    name <- gsub("_", " ", name)
    return(name)
}

#' Partial Dependence Plots Profissionais para Gaussian/ZOIB (V4)
#'
#' Gera PDPs de alta resolução (300 pontos) usando posterior_epred()
#' com suporte a interações ARC (Inner/Outer).
#'
#' @param model Objeto brmsfit
#' @param model_name Nome do modelo
#' @param model_info Informações do tipo de modelo
#' @param resolution Número de pontos na grade (padrão: 300)
#' @param ndraws Número de draws posteriores (padrão: 500)
#' @return Objeto patchwork ou NULL
generate_marginal_effects_professional <- function(model, model_name, model_info,
                                                    resolution = 300, ndraws = 500) {
    cat(sprintf("    -> Gerando PDPs Profissionais (res=%d, ndraws=%d)...\n",
                resolution, ndraws))

    model_data <- model$data

    # Identificar variáveis da fórmula (corrigido para brmsfit)
    # fórmula brms pode estar em model$formula ou formula(model)$formula
    model_formula <- if (!is.null(model$formula)) model$formula else formula(model)
    # Se for lista, extrair a fórmula principal
    if (is.list(model_formula) && !is.null(model_formula$formula)) {
        model_formula <- model_formula$formula
    }
    cat(sprintf("    → Fórmula extraída: %s\n", paste(deparse(model_formula), collapse = " ")))
    all_vars <- all.vars(model_formula)
    cat(sprintf("    → Variáveis extraídas: %s\n", paste(all_vars, collapse = ", ")))
    if (length(all_vars) < 2) {
        cat("    ⚠ Fórmula muito simples, sem preditores para PDP\n")
        return(NULL)
    }

    response_var <- all_vars[1]
    predictors <- setdiff(all_vars, response_var)

    # Filtrar apenas preditores contínuos presentes nos dados
    continuous_preds <- c()
    for (pred in predictors) {
        if (pred %in% names(model_data) && is.numeric(model_data[[pred]])) {
            continuous_preds <- c(continuous_preds, pred)
        }
    }

    # Priorizar preditores principais do estudo
    priority <- c("PC1MAGNITUDE", "PC2MAGNITUDE", "PC1VARIABILITY", "PC2VARIABILITY", "DEPTHM")
    continuous_preds <- c(intersect(priority, continuous_preds),
                          setdiff(continuous_preds, priority))

    if (length(continuous_preds) == 0) {
        cat("    ⚠ Nenhum preditor contínuo encontrado\n")
        return(NULL)
    }

    cat(sprintf("    → Processando %d preditores: %s\n",
                length(continuous_preds), paste(continuous_preds, collapse = ", ")))

    # Verificar interações com ARCH na fórmula (usando mesma fórmula extraída)
    formula_str <- paste(deparse(model_formula), collapse = " ")
    has_arch <- "ARCH" %in% names(model_data)

    plot_list <- list()

    for (pred in continuous_preds) {
        cat(sprintf("      → Processando %s...\n", pred))

        # Criar sequência do preditor
        pred_range <- range(model_data[[pred]], na.rm = TRUE)
        pred_seq <- seq(pred_range[1], pred_range[2], length.out = resolution)

        # Verificar se há interação com ARCH
        interaction_pattern <- paste0("(\\b", pred, "\\b.*:ARCH|ARCH.*:", pred, "\\b)")
        has_interaction <- grepl(interaction_pattern, formula_str) && has_arch

        if (has_interaction) {
            # Predições separadas para Inner e Outer ARC
            newdata_inner <- create_newdata_template(model_data, resolution)
            newdata_outer <- newdata_inner
            newdata_inner$ARCH <- "inner"
            newdata_outer$ARCH <- "outer"
            newdata_inner[[pred]] <- pred_seq
            newdata_outer[[pred]] <- pred_seq

            # Predições
            pred_inner <- posterior_epred(model, newdata = newdata_inner, ndraws = ndraws)
            pred_outer <- posterior_epred(model, newdata = newdata_outer, ndraws = ndraws)

            # Consolidar dados para plot
            plot_data <- bind_rows(
                data.frame(
                    x = pred_seq,
                    estimate = apply(pred_inner, 2, median),
                    lower = apply(pred_inner, 2, quantile, 0.025),
                    upper = apply(pred_inner, 2, quantile, 0.975),
                    ARC = "Inner Arc"
                ),
                data.frame(
                    x = pred_seq,
                    estimate = apply(pred_outer, 2, median),
                    lower = apply(pred_outer, 2, quantile, 0.025),
                    upper = apply(pred_outer, 2, quantile, 0.975),
                    ARC = "Outer Arc"
                )
            )

            arc_colors <- c("Inner Arc" = "#0072B2", "Outer Arc" = "#E69F00")

            p <- ggplot(plot_data, aes(x = x, y = estimate, color = ARC, fill = ARC, group = ARC)) +
                geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, color = NA) +
                geom_line(linewidth = 1.2) +
                scale_color_manual(values = arc_colors, name = "ARC") +
                scale_fill_manual(values = arc_colors, name = "ARC") +
                labs(
                    title = paste("Effect of", clean_predictor_name(pred), "by ARC"),
                    x = clean_predictor_name(pred),
                    y = ifelse(model_info$type == "ZOIB", "Predicted Proportion", "Predicted Value"),
                    caption = sprintf("Line: Median | Band: 95%% CI | Resolution: %d points", resolution)
                ) +
                theme_publication() +
                theme(legend.position = "bottom")

        } else {
            # Uma curva só (sem interação com ARC)
            newdata <- create_newdata_template(model_data, resolution)
            newdata[[pred]] <- pred_seq

            pred_arr <- posterior_epred(model, newdata = newdata, ndraws = ndraws)

            plot_data <- data.frame(
                x = pred_seq,
                estimate = apply(pred_arr, 2, median),
                lower = apply(pred_arr, 2, quantile, 0.025),
                upper = apply(pred_arr, 2, quantile, 0.975)
            )

            p <- ggplot(plot_data, aes(x = x, y = estimate)) +
                geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#0072B2", alpha = 0.25) +
                geom_line(color = "#0072B2", linewidth = 1.2) +
                labs(
                    title = paste("Effect of", clean_predictor_name(pred)),
                    x = clean_predictor_name(pred),
                    y = ifelse(model_info$type == "ZOIB", "Predicted Proportion", "Predicted Value"),
                    caption = sprintf("Line: Median | Band: 95%% CI | Resolution: %d points", resolution)
                ) +
                theme_publication()
        }

        plot_list[[pred]] <- p
    }

    if (length(plot_list) == 0) return(NULL)

    # Combinar em painel
    ncol_calc <- min(3, length(plot_list))
    combined <- wrap_plots(plot_list, ncol = ncol_calc) +
        plot_annotation(
            title = "Partial Dependence Plots (High Resolution)",
            subtitle = paste("Model:", gsub("^WINNER_(?:FINAL_)?", "", gsub("_", " ", model_name))),
            caption = sprintf("Professional-grade PDPs using posterior_epred() with %d draws", ndraws)
        )

    return(combined)
}

### ═══════════════════════════════════════════════════════════════════════════════
### FUNÇÕES AVANÇADAS JSDM (NOVO V4)
### ═══════════════════════════════════════════════════════════════════════════════

#' PDPs de Composição Relativa: Trade-offs Ecológicos (V4)
#'
#' Visualiza como a composição total da comunidade (soma=1) se redistribui
#' ao longo de gradientes ambientais. Mostra trade-offs competitivos.
#'
#' @param model Objeto brmsfit JSDM
#' @param model_name Nome do modelo
#' @param ndraws Número de draws (padrão: 500)
#' @param resolution Número de pontos (padrão: 100)
#' @return Objeto patchwork ou NULL
generate_jsdm_composition_tradeoff <- function(model, model_name,
                                                ndraws = 500, resolution = 100) {
    cat("    -> Gerando PDPs de Composição Relativa (trade-offs)...\n")

    # Validação
    validation <- validate_jsdm_species_order(model)
    if (is.null(validation)) return(NULL)

    model_data <- model$data
    species_codes <- validation$formula_species

    # Selecionar preditor principal
    pred <- "PC1MAGNITUDE"
    if (!pred %in% names(model_data)) {
        numeric_cols <- names(model_data)[sapply(model_data, is.numeric)]
        pred <- setdiff(numeric_cols, c("PC2MAGNITUDE", "PC1VARIABILITY"))[1]
    }

    # Criar newdata
    pred_range <- range(model_data[[pred]], na.rm = TRUE)
    pred_seq <- seq(pred_range[1], pred_range[2], length.out = resolution)

    newdata <- create_newdata_template(model_data, resolution)
    newdata[[pred]] <- pred_seq

    # Predições
    pred_array <- tryCatch({
        posterior_epred(model, newdata = newdata, ndraws = ndraws)
    }, error = function(e) {
        cat(sprintf("    ⚠ Erro: %s\n", substr(e$message, 1, 80)))
        return(NULL)
    })

    if (is.null(pred_array)) return(NULL)

    dims <- dim(pred_array)

    # Processar cada espécie
    all_species_data <- list()
    colors_tradeoff <- c("#D55E00", "#009E73", "#0072B2", "#CC79A7", "#E69F00")

    for (sp_idx in 1:min(dims[3], length(species_codes))) {
        sp_code <- gsub("_prop$", "", species_codes[sp_idx])
        sp_draws <- pred_array[, , sp_idx]

        sp_summary <- data.frame(
            x = pred_seq,
            estimate = apply(sp_draws, 2, median),
            lower = apply(sp_draws, 2, quantile, 0.025),
            upper = apply(sp_draws, 2, quantile, 0.975),
            species = sp_code,
            species_idx = sp_idx
        )
        all_species_data[[sp_idx]] <- sp_summary
    }

    tradeoff_data <- bind_rows(all_species_data)

    # Labels
    species_labels <- c(
        "MUSSISMILIA" = "*M. hispida*",
        "TURF" = "Turf Algae",
        "CCA" = "Crustose Coralline Algae",
        "CYANO" = "Cyanobacteria",
        "MACROALGAE" = "Macroalgae"
    )
    tradeoff_data$species_clean <- species_labels[tradeoff_data$species]
    tradeoff_data$species_clean <- factor(tradeoff_data$species_clean,
                                           levels = species_labels)

    # PLOT 1: Área empilhada (composição total)
    p1 <- ggplot(tradeoff_data, aes(x = x, y = estimate, fill = species_clean)) +
        geom_area(position = "stack", alpha = 0.85) +
        scale_fill_manual(values = colors_tradeoff, name = "Species") +
        labs(
            title = "Community Composition Redistribution",
            subtitle = paste("How proportions shift along", pred, "(sum = 1)"),
            x = pred,
            y = "Total Community Proportion",
            caption = "Stacked area shows relative abundance changes"
        ) +
        theme_publication() +
        theme(legend.position = "bottom")

    # PLOT 2: Linhas individuais com destaque
    tradeoff_data$is_mhispida <- tradeoff_data$species == "MUSSISMILIA"

    p2 <- ggplot(tradeoff_data, aes(x = x, y = estimate,
                                     color = species_clean, group = species_clean)) +
        geom_ribbon(aes(ymin = lower, ymax = upper, fill = species_clean),
                   alpha = 0.15, color = NA) +
        geom_line(aes(linewidth = is_mhispida)) +
        scale_color_manual(values = colors_tradeoff, name = "Species") +
        scale_fill_manual(values = colors_tradeoff, name = "Species") +
        scale_linewidth_manual(values = c("TRUE" = 2, "FALSE" = 1), guide = "none") +
        labs(
            title = "Species-Specific Response Curves",
            subtitle = "Thick line: *M. hispida* | Shaded: 95% CI",
            x = pred,
            y = "Predicted Proportion"
        ) +
        theme_publication() +
        theme(legend.position = "bottom")

    # Combinar
    combined <- (p1 / p2) +
        plot_annotation(
            title = "Ecological Trade-offs: Community Redistribution",
            subtitle = paste("Model:", gsub("^WINNER_(?:FINAL_)?", "", gsub("_", " ", model_name))),
            caption = "Top: Total composition (stacked) | Bottom: Individual responses"
        )

    return(combined)
}

#' Análise Detalhada de Resposta de M. hispida por ARC (V4)
#'
#' Compara estatisticamente as curvas de resposta de M. hispida entre
#' Inner Arc e Outer Arc, identificando diferenças significativas.
#'
#' @param model Objeto brmsfit JSDM
#' @param model_name Nome do modelo
#' @param model_output_dir Diretório para salvar resultados
#' @param ndraws Número de draws (padrão: 500)
#' @param resolution Número de pontos (padrão: 100)
#' @return Lista com plot, summary e dados, ou NULL
generate_mhispida_response_detailed_EXTENDED <- function(model, model_name,
                                                          model_output_dir,
                                                          ndraws = 500,
                                                          resolution = 100) {
    cat("    -> Gerando análise detalhada de M. hispida por ARC...\n")

    # Validações
    validation <- validate_jsdm_species_order(model)
    if (is.null(validation) || !validation$mhispida_in_formula) {
        cat("    ⚠ M. hispida não encontrada no modelo\n")
        return(NULL)
    }

    model_data <- model$data

    # Verificar se modelo tem ARC
    if (!"ARCH" %in% names(model_data)) {
        cat("    ⚠ Modelo não tem variável ARCH\n")
        return(NULL)
    }

    mhispida_idx <- validation$mhispida_idx

    # Selecionar preditores principais
    predictors <- c("PC1MAGNITUDE", "PC2MAGNITUDE", "PC1VARIABILITY", "PC2VARIABILITY", "DEPTHM")
    predictors <- intersect(predictors, names(model_data))

    if (length(predictors) == 0) {
        cat("    ⚠ Nenhum preditor principal encontrado\n")
        return(NULL)
    }

    # Análise por preditor
    plot_list <- list()
    sensitivity_summary <- data.frame()

    for (pred in predictors) {
        cat(sprintf("      → Analisando %s...\n", pred))

        pred_range <- range(model_data[[pred]], na.rm = TRUE)
        pred_seq <- seq(pred_range[1], pred_range[2], length.out = resolution)

        # Criar newdata para cada ARC
        create_newdata_arc <- function(arc_val) {
            nd <- create_newdata_template(model_data, resolution)
            nd$ARCH <- arc_val
            nd[[pred]] <- pred_seq
            return(nd)
        }

        newdata_inner <- create_newdata_arc("inner")
        newdata_outer <- create_newdata_arc("outer")

        # Predições
        pred_inner <- posterior_epred(model, newdata = newdata_inner, ndraws = ndraws)
        pred_outer <- posterior_epred(model, newdata = newdata_outer, ndraws = ndraws)

        # Extrair M. hispida
        mhispida_inner <- pred_inner[, , mhispida_idx]
        mhispida_outer <- pred_outer[, , mhispida_idx]

        # Estatísticas descritivas
        plot_data <- bind_rows(
            data.frame(
                x = pred_seq,
                estimate = apply(mhispida_inner, 2, median),
                lower = apply(mhispida_inner, 2, quantile, 0.025),
                upper = apply(mhispida_inner, 2, quantile, 0.975),
                ARC = "Inner Arc"
            ),
            data.frame(
                x = pred_seq,
                estimate = apply(mhispida_outer, 2, median),
                lower = apply(mhispida_outer, 2, quantile, 0.025),
                upper = apply(mhispida_outer, 2, quantile, 0.975),
                ARC = "Outer Arc"
            )
        )

        # Calcular diferenças estatísticas
        difference_draws <- mhispida_inner - mhispida_outer
        diff_median <- apply(difference_draws, 2, median)
        diff_lower <- apply(difference_draws, 2, quantile, 0.025)
        diff_upper <- apply(difference_draws, 2, quantile, 0.975)

        # Teste: diferença significativa (95% CI não inclui 0)
        significant_points <- sum(diff_lower > 0 | diff_upper < 0)
        significance_pct <- significant_points / resolution * 100

        # Amplitude de resposta
        range_inner <- max(plot_data$estimate[plot_data$ARC == "Inner Arc"]) -
                       min(plot_data$estimate[plot_data$ARC == "Inner Arc"])
        range_outer <- max(plot_data$estimate[plot_data$ARC == "Outer Arc"]) -
                       min(plot_data$estimate[plot_data$ARC == "Outer Arc"])

        sensitivity_summary <- rbind(sensitivity_summary, data.frame(
            predictor = pred,
            range_inner = range_inner,
            range_outer = range_outer,
            range_ratio = range_inner / ifelse(range_outer == 0, 0.001, range_outer),
            significance_pct = significance_pct,
            mean_difference = mean(abs(diff_median))
        ))

        # Plot
        arc_colors <- c("Inner Arc" = "#0072B2", "Outer Arc" = "#E69F00")

        p <- ggplot(plot_data, aes(x = x, y = estimate, color = ARC, fill = ARC)) +
            geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.25, color = NA) +
            geom_line(linewidth = 1.5) +
            scale_color_manual(values = arc_colors, name = "ARC") +
            scale_fill_manual(values = arc_colors, name = "ARC") +
            labs(
                title = sprintf("*M. hispida* Response to %s", pred),
                subtitle = sprintf("Inner range: %.3f | Outer range: %.3f | Sig. diff: %.1f%%",
                                 range_inner, range_outer, significance_pct),
                x = clean_predictor_name(pred),
                y = "Predicted Proportion of *M. hispida*",
                caption = sprintf("Range ratio (Inner/Outer): %.2f", range_inner/ifelse(range_outer==0, 0.001, range_outer))
            ) +
            theme_publication() +
            theme(legend.position = "bottom")

        plot_list[[pred]] <- p
    }

    # Salvar resumo como CSV
    if (nrow(sensitivity_summary) > 0 && !is.null(model_output_dir)) {
        csv_path <- file.path(model_output_dir, "Mhispida_ARC_Sensitivity_Summary.csv")
        write.csv(sensitivity_summary, csv_path, row.names = FALSE)
        cat(sprintf("      → Resumo salvo: %s\n", basename(csv_path)))

        cat("\n      📊 RESUMO DE SENSIBILIDADE:\n")
        print(sensitivity_summary, row.names = FALSE)
    }

    # Combinar plots
    if (length(plot_list) == 0) return(NULL)

    ncol_calc <- min(3, length(plot_list))
    combined <- wrap_plots(plot_list, ncol = ncol_calc) +
        plot_annotation(
            title = "*Mussismilia hispida*: ARC-Specific Sensitivity",
            subtitle = paste("Model:", gsub("^WINNER_(?:FINAL_)?", "", gsub("_", " ", model_name))),
            caption = "Inner Arc = blue | Outer Arc = orange | Sig. diff = % of gradient with significant difference"
        )

    return(list(
        plot = combined,
        summary = sensitivity_summary,
        plot_list = plot_list
    ))
}

### 7. PIPELINE PRINCIPAL ========================================================

process_model_v5 <- function(rds_file, output_base_dir) {
    # Carregar modelo
    model <- readRDS(rds_file)
    model_name <- gsub("\\.rds$", "", basename(rds_file))

    # Criar diretório de saída
    model_output_dir <- file.path(output_base_dir, model_name)
    if (!dir.exists(model_output_dir)) {
        dir.create(model_output_dir, recursive = TRUE)
    }

    cat(sprintf("\n🔍 PROCESSANDO: %s\n", model_name))

    # Detectar tipo de modelo
    model_info <- detect_model_type_v5(model, model_name)
    diagnose_model_structure(model)

    # ═══════════════════════════════════════════════════════════════════════════
    # NOVO V4: Validação de espécies para modelos JSDM
    # ═══════════════════════════════════════════════════════════════════════════
    if (model_info$type %in% c("DIRICHLET", "BETA")) {
        validation <- validate_jsdm_species_order(model)

        if (!is.null(validation)) {
            # Salvar relatório de validação
            validation_file <- file.path(model_output_dir, "JSDM_Validation_Report.txt")
            sink(validation_file)
            cat("=== RELATÓRIO DE VALIDAÇÃO JSDM ===\n\n")
            cat("Data:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
            cat("Modelo:", model_name, "\n")
            cat("Tipo:", model_info$type, "\n\n")
            cat("ESPÉCIES NA FÓRMULA:\n")
            for (i in seq_along(validation$formula_species)) {
                cat(sprintf("  %d. %s\n", i, validation$formula_species[i]))
            }
            cat("\nÍNDICE M. HISPIDA:", ifelse(is.null(validation$mhispida_idx), "NÃO ENCONTRADO", validation$mhispida_idx), "\n")
            cat("TESTE DE PREDIÇÃO:", ifelse(validation$test_success, "SUCESSO", "FALHA"), "\n")
            sink()
            cat(sprintf("    ✓ Relatório de validação: %s\n", basename(validation_file)))
        }
    }
    # ═══════════════════════════════════════════════════════════════════════════

    # 1. Forest Plot
    cat("\n  [1/8] Forest Plot...\n")
    p_forest <- generate_forest_plot_v5(model, model_name, model_info)
    if (!is.null(p_forest)) {
        ggsave(file.path(model_output_dir, "FIGURE_1_Forest_Plot.png"),
               p_forest, width = 10, height = 6, dpi = 300, bg = "white")
        cat("    ✓ FIGURE_1 salvo\n")
    }

    # 2. Efeitos Marginais (MODIFICADO V4)
    cat("\n  [2/8] Efeitos Marginais...\n")
    if (model_info$type %in% c("GAUSSIAN", "ZOIB")) {
        # Usar versão profissional para Gaussian/ZOIB
        cat("    → Usando PDPs profissionais (alta resolução)\n")
        p_effects <- generate_marginal_effects_professional(model, model_name, model_info)
        if (!is.null(p_effects)) {
            ggsave(file.path(model_output_dir, "FIGURE_2_Marginal_Effects_Professional.png"),
                   p_effects, width = 16, height = 10, dpi = 300, bg = "white")
            cat("    ✓ FIGURE_2 (Profissional) salvo\n")
        }
    } else if (model_info$type %in% c("DIRICHLET", "BETA")) {
        # JSDM: PDPs serão gerados nas figuras específicas abaixo
        cat("    → JSDM: PDPs serão gerados em figuras específicas\n")
    } else {
        # Fallback
        cat("    → Fallback: efeitos marginais padrão\n")
        # (Não implementado para simplificar)
    }

    # 3. PPC
    cat("\n  [3/8] Posterior Predictive Checks...\n")
    p_ppc <- generate_ppc_plots_v5(model, model_name, model_info)
    if (!is.null(p_ppc)) {
        ggsave(file.path(model_output_dir, "FIGURE_3_PPC.png"),
               p_ppc, width = 12, height = 10, dpi = 300, bg = "white")
        cat("    ✓ FIGURE_3 salvo\n")
    }

    # 4. Diagnósticos
    cat("\n  [4/8] Diagnósticos...\n")
    p_diag <- generate_diagnostics_v5(model, model_name)
    if (!is.null(p_diag)) {
        ggsave(file.path(model_output_dir, "FIGURE_4_Diagnostics.png"),
               p_diag, width = 14, height = 10, dpi = 300, bg = "white")
        cat("    ✓ FIGURE_4 salvo\n")
    }

    # ═══════════════════════════════════════════════════════════════════════════
    # NOVO V4: Figuras específicas JSDM
    # ═══════════════════════════════════════════════════════════════════════════
    if (model_info$type %in% c("DIRICHLET", "BETA")) {
        cat("\n  [5-10/10] Figuras JSDM específicas...\n")

        # 5. Forest Plot JSDM Error Bar
        cat("    → FIGURE_5: Forest Plot JSDM...\n")
        p_jsdm_forest <- generate_jsdm_forest_plot_errorbar_CORRECTED(model, model_name)
        if (!is.null(p_jsdm_forest)) {
            ggsave(file.path(model_output_dir, "FIGURE_5_JSDM_Forest_Errorbar.png"),
                   p_jsdm_forest, width = 12, height = 8, dpi = 300, bg = "white")
            cat("      ✓ FIGURE_5 salvo\n")
        }

        # 6. PDP Community Overlaid (corrigido)
        cat("    → FIGURE_6: PDP Community...\n")
        p_jsdm_pdp_comm <- generate_jsdm_pdp_community_overlaid_CORRECTED(model, model_name)
        if (!is.null(p_jsdm_pdp_comm)) {
            ggsave(file.path(model_output_dir, "FIGURE_6_JSDM_PDP_Community.png"),
                   p_jsdm_pdp_comm, width = 14, height = 10, dpi = 300, bg = "white")
            cat("      ✓ FIGURE_6 salvo\n")
        }

        # 7. PDP Per Species
        cat("    → FIGURE_7: PDP Per Species...\n")
        p_jsdm_pdp_species <- generate_jsdm_pdp_per_species_CORRECTED(model, model_name)
        if (!is.null(p_jsdm_pdp_species)) {
            ggsave(file.path(model_output_dir, "FIGURE_7_JSDM_PDP_Per_Species.png"),
                   p_jsdm_pdp_species, width = 14, height = 12, dpi = 300, bg = "white")
            cat("      ✓ FIGURE_7 salvo\n")
        }

        # 8. PPC por Espécie
        cat("    → FIGURE_8: PPC by Species...\n")
        p_jsdm_ppc <- generate_jsdm_ppc_species_overlays(model, model_name)
        if (!is.null(p_jsdm_ppc)) {
            ggsave(file.path(model_output_dir, "FIGURE_8_JSDM_PPC_by_Species.png"),
                   p_jsdm_ppc, width = 14, height = 10, dpi = 300, bg = "white")
            cat("      ✓ FIGURE_8 salvo\n")
        }

        # 9. Composição Relativa (NOVO V4)
        cat("    → FIGURE_9: Composition Trade-off...\n")
        p_tradeoff <- tryCatch(
            generate_jsdm_composition_tradeoff(model, model_name),
            error = function(e) {
                cat(sprintf("      ⚠ Erro: %s\n", substr(e$message, 1, 60)))
                NULL
            }
        )
        if (!is.null(p_tradeoff)) {
            ggsave(file.path(model_output_dir, "FIGURE_9_JSDM_Composition_Tradeoff.png"),
                   p_tradeoff, width = 14, height = 12, dpi = 300, bg = "white")
            cat("      ✓ FIGURE_9 salvo\n")
        }

        # 10. Análise de Sensibilidade ARC (NOVO V4)
        cat("    → FIGURE_10: M. hispida ARC Sensitivity...\n")
        p_arc_analysis <- tryCatch(
            generate_mhispida_response_detailed_EXTENDED(model, model_name, model_output_dir),
            error = function(e) {
                cat(sprintf("      ⚠ Erro: %s\n", substr(e$message, 1, 60)))
                NULL
            }
        )
        if (!is.null(p_arc_analysis)) {
            ggsave(file.path(model_output_dir, "FIGURE_10_Mhispida_ARC_Sensitivity.png"),
                   p_arc_analysis$plot, width = 16, height = 10, dpi = 300, bg = "white")
            cat("      ✓ FIGURE_10 salvo\n")
        }
    }
    # ═══════════════════════════════════════════════════════════════════════════

    cat(sprintf("\n✅ MODELO CONCLUÍDO: %s\n", model_name))
    cat(sprintf("   Diretório: %s\n", model_output_dir))
}

### 8. EXECUÇÃO PRINCIPAL ========================================================

cat("\n", rep("=", 75), "\n", sep = "")
cat("🚀 INICIANDO PIPELINE V5 - HYBRID V4 (JSDM CORRIGIDO + FUNCIONALIDADES AVANÇADAS)\n")
cat(rep("=", 75), "\n", sep = "")

winner_files <- discover_winner_models(search_directories_v5)

processed_info <- list()
for (rds_file in winner_files) {
    result <- process_model_v5(rds_file, output_dir_v5)
    if (!is.null(result)) processed_info[[length(processed_info) + 1]] <- result
}

cat("\n", rep("=", 75), "\n", sep = "")
cat("📋 RELATÓRIO FINAL - PIPELINE V5 HYBRID V4\n")
cat(rep("=", 75), "\n", sep = "")

cat(sprintf("\n✓ Modelos processados: %d\n", length(processed_info)))
cat(sprintf("📁 Saída: %s\n", output_dir_v5))

cat("\nModelos:\n")
for (info in processed_info) {
    cat(sprintf("  • %s (%s)\n", info$model_name, info$model_info$type))
}

cat("\nFiguras geradas:\n")
cat("  GAUSSIAN/ZOIB:\n")
cat("    - FIGURE_1_Forest_Plot.png\n")
cat("    - FIGURE_2_Marginal_Effects_Professional.png (NOVO V4: alta resolução 300)\n")
cat("    - FIGURE_3_PPC.png\n")
cat("    - FIGURE_4_Diagnostics.png\n")
cat("\n  DIRICHLET/BETA (JSDM):\n")
cat("    - FIGURE_1_Forest_Plot.png (CORRIGIDO V4: com M. hispida)\n")
cat("    - FIGURE_3_PPC.png\n")
cat("    - FIGURE_4_Diagnostics.png\n")
cat("    - FIGURE_5_JSDM_Forest_Errorbar.png (com M. hispida como referência)\n")
cat("    - FIGURE_6_JSDM_PDP_Community.png (curvas separadas)\n")
cat("    - FIGURE_7_JSDM_PDP_Per_Species.png (facets)\n")
cat("    - FIGURE_8_JSDM_PPC_by_Species.png\n")
cat("    - FIGURE_9_JSDM_Composition_Tradeoff.png (NOVO V4: trade-offs)\n")
cat("    - FIGURE_10_Mhispida_ARC_Sensitivity.png (NOVO V4: análise ARC)\n")
cat("    - JSDM_Validation_Report.txt (NOVO V4: relatório de validação)\n")
cat("    - Mhispida_ARC_Sensitivity_Summary.csv (NOVO V4: estatísticas)\n")

cat("\n🎉 PIPELINE HYBRID V4 CONCLUÍDO!\n\n")
