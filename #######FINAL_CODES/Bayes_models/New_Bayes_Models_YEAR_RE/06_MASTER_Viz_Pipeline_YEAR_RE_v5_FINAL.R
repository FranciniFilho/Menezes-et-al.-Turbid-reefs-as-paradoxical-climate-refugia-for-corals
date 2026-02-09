### ==================================================================================== ###
### === MASTER VISUALIZATION PIPELINE V5 - YEAR RANDOM EFFECT MODELS                   === ###
### ===                                                                              === ###
### === Based on: 05_MASTER_Viz_Pipeline_v5_FINAL.R                                 === ###
### === Adds: YEAR random effect diagnostic plots                                   === ###
### === Maintains: EXACT theme_publication() and Okabe-Ito palette                  === ###
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

### 3. CONFIGURAÇÃO DOS DIRETÓRIOS -----------------------------------------------
search_directories_YEAR_RE <- c(
    "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/ZOIB_Abundance_YEAR_RE",
    "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/Gaussian_RGR",
    "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/Gaussian_Health_YEAR_RE",
    "C:/Users/rbfra/OneDrive/New_Bayes_Models_Output/JSDM_Dirichlet_YEAR_RE"
)

output_dir_YEAR_RE <- "C:/Users/rbfra/OneDrive/Bayesian_Figures_YEAR_RE_v5"
dir.create(output_dir_YEAR_RE, showWarnings = FALSE, recursive = TRUE)

cat("\n", rep("=", 75), "\n", sep = "")
cat("MASTER VISUALIZATION PIPELINE V5 - YEAR RANDOM EFFECT MODELS\n")
cat(rep("=", 75), "\n", sep = "")
cat("Output directory:", output_dir_YEAR_RE, "\n")

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

detect_model_type_year_re <- function(model, model_name) {
    if (!inherits(model, "brmsfit")) {
        return(list(type = "UNKNOWN", response_label = "Unknown",
                    uses_arch = FALSE, has_year_re = FALSE))
    }
    fam <- family(model)$family

    # Check for YEAR random effect
    formula_str <- paste(deparse(formula(model)), collapse = " ")
    has_year_re <- grepl("\\(1.*\\|.*YEAR\\)", formula_str) ||
                   grepl("\\(1.*\\|.*REEF.*YEAR\\)", formula_str)

    if (grepl("zero_one_inflated_beta", fam, ignore.case = TRUE)) {
        return(list(type = "ZOIB",
                    response_label = ifelse(has_year_re, "Cover (Proportion) + Year RE", "Cover (Proportion)"),
                    uses_arch = TRUE, has_year_re = has_year_re))
    } else if (grepl("gaussian", fam, ignore.case = TRUE)) {
        resp_var <- all.vars(formula(model))[1]
        label_suffix <- ifelse(has_year_re, " + Year RE", "")
        return(list(type = "GAUSSIAN",
                    response_label = paste(resp_var, label_suffix),
                    uses_arch = FALSE, has_year_re = has_year_re))
    } else if (grepl("dirichlet", fam, ignore.case = TRUE)) {
        return(list(type = "DIRICHLET",
                    response_label = ifelse(has_year_re, "Composition + Year RE", "Composition"),
                    uses_arch = TRUE, has_year_re = has_year_re))
    }
    return(list(type = "UNKNOWN", response_label = "Unknown",
                uses_arch = FALSE, has_year_re = FALSE))
}

### 5. YEAR RANDOM EFFECT DIAGNOSTIC PLOT (NEW) ----------------------------------

#' Generate YEAR random effect diagnostic plot
#' @param model brmsfit object
#' @return ggplot2 object or NULL
generate_year_re_diagnostics <- function(model) {

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
    # year_re_list can be a list with one element (the parameter) or a matrix
    if (is.list(year_re_list) && !is.data.frame(year_re_list)) {
        # If it's a list, get the first element (usually the intercept)
        if (length(year_re_list) > 0) {
            year_re <- year_re_list[[1]]
        } else {
            return(NULL)
        }
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
        # If year_re has no names, use sequential numbers or extract from model data
        year_names <- names(year_re)
        if (is.null(year_names) || all(year_names == "")) {
            # Try to get years from model data
            if ("YEAR" %in% names(model$data)) {
                year_names <- sort(unique(model$data$YEAR))
            } else {
                # Fall back to sequential numbers
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

    # Clean YEAR names (remove factor levels if present)
    year_df$YEAR <- gsub("^.*\\)|\\(|\\)", "", year_df$YEAR)
    year_df$YEAR <- as.integer(as.character(year_df$YEAR))

    # Remove NA years and sort
    year_df <- year_df[!is.na(year_df$YEAR), ]
    year_df <- year_df[order(year_df$YEAR), ]

    if (nrow(year_df) == 0) {
        return(NULL)
    }

    # Create plot using EXACT theme
    caption_text <- "Points: Median | Line: Trend"
    if (any(!is.na(year_df$SE))) {
        caption_text <- "Points: Median | Error bars: 95% CI | Line: Trend"
    }

    p <- ggplot(year_df, aes(x = YEAR, y = Estimate)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
        geom_point(size = 3, color = okabe_ito[1]) +
        geom_line(aes(group = 1), color = okabe_ito[1], alpha = 0.5) +
        labs(
            title = "Year Random Effects",
            subtitle = "Deviation from overall intercept by year",
            x = "Year",
            y = "Random Effect (deviation)",
            caption = caption_text
        ) +
        theme_publication()

    # Add error bars if SE is available
    if (any(!is.na(year_df$SE))) {
        p <- p + geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE),
                               width = 0.2, linewidth = 0.8)
    }

    return(p)
}

### 6. FOREST PLOT GENERATION ----------------------------------------------------

generate_forest_plot_year_re <- function(model, model_name, model_info) {
    cat("    -> Gerando Forest Plot...\n")

    # For DIRICHLET models, use the specialized JSDM function
    if (model_info$type == "DIRICHLET") {
        cat("    -> Modelo DIRICHLET detectado: usando função JSDM\n")
        # This would need the JSDM-specific function from original
        # For now, return a simplified version
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

    # Add YEAR RE variance annotation if applicable
    subtitle_text <- paste("Model:", gsub("^WINNER_", "", gsub("\\.rds$", "", model_name)))
    if (model_info$has_year_re) {
        # Extract YEAR SD
        year_sd <- tryCatch({
            summary(model)$random$YEAR["sd_(Intercept)", "Estimate"]
        }, error = function(e) NA)
        if (!is.na(year_sd)) {
            subtitle_text <- paste(subtitle_text, sprintf("| Year SD: %.3f", year_sd))
        }
    }

    if (model_info$type == "GAUSSIAN") {
        p <- ggplot(draws_long, aes(y = reorder(parameter, abs(value)), x = value)) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.7) +
            stat_slab(aes(fill = after_stat(x > 0)), alpha = 0.5, slab_size = 0.3, scale = 0.6) +
            stat_pointinterval(aes(fill = after_stat(x > 0)), point_interval = median_qi,
                              point_size = 1.2, interval_size = 0.8, .width = c(0.89, 0.95)) +
            scale_fill_manual(values = c("TRUE" = color_primary, "FALSE" = "gray60"), guide = "none") +
            labs(title = paste("Magnitude of Fixed Effects -", model_info$type),
                 subtitle = subtitle_text,
                 x = "Effect Estimate", y = NULL, caption = "Points: Median | Bars: 89% and 95% CI") +
            theme_publication()
    } else {
        p <- ggplot(draws_long, aes(y = reorder(parameter, abs(value)), x = value, fill = after_stat(x > 0))) +
            stat_halfeye(.width = c(0.89, 0.95), alpha = 0.85, point_interval = median_qi,
                        point_size = 1.5, slab_alpha = 0.6, slab_size = 0.3) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.7) +
            scale_fill_manual(values = c("TRUE" = color_primary, "FALSE" = "gray60"), guide = "none") +
            labs(title = paste("Magnitude of Fixed Effects -", model_info$type),
                 subtitle = subtitle_text,
                 x = "Effect Estimate", y = NULL, caption = "Points: Median | Bars: 89% and 95% CI") +
            theme_publication()
    }
    return(p)
}

### 7. PPC PLOTS ----------------------------------------------------------------

generate_ppc_plots_year_re <- function(model, model_name, model_info) {
    cat("    -> PPC...\n")

    # Dirichlet models don't support pp_check
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

### 8. MARGINAL EFFECTS PLOTS ----------------------------------------------------

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

clean_predictor_name <- function(name) {
    name <- gsub("MAGNITUDE", "Magnitude", name)
    name <- gsub("VARIABILITY", "Variability", name)
    name <- gsub("HABMERGED", "Habitat", name)
    name <- gsub("DEPTHM", "Depth", name)
    return(name)
}

generate_marginal_effects_year_re <- function(model, model_name, model_info,
                                             resolution = 200, ndraws = 500) {
    cat("    -> Marginal Effects...\n")

    if (model_info$type == "DIRICHLET") {
        cat("    ⚠ Marginal effects simplificado para DIRICHLET\n")
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

    continuous_preds <- c()
    for (pred in predictors) {
        if (pred %in% names(model_data) && is.numeric(model_data[[pred]])) {
            continuous_preds <- c(continuous_preds, pred)
        }
    }

    priority <- c("PC1MAGNITUDE", "PC2MAGNITUDE", "PC1VARIABILITY", "PC2VARIABILITY", "DEPTHM")
    continuous_preds <- c(intersect(priority, continuous_preds),
                          setdiff(continuous_preds, priority))

    if (length(continuous_preds) == 0) return(NULL)

    plot_list <- list()

    for (pred in head(continuous_preds, 4)) {
        pred_range <- range(model_data[[pred]], na.rm = TRUE)
        pred_seq <- seq(pred_range[1], pred_range[2], length.out = resolution)

        newdata <- create_newdata_template(model_data, resolution)
        newdata[[pred]] <- pred_seq

        pred_arr <- tryCatch({
            posterior_epred(model, newdata = newdata, ndraws = ndraws)
        }, error = function(e) NULL)

        if (is.null(pred_arr)) next

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

        plot_list[[pred]] <- p
    }

    if (length(plot_list) == 0) return(NULL)

    combined <- wrap_plots(plot_list, ncol = 2) +
        plot_annotation(
            title = "Partial Dependence Plots",
            subtitle = paste("Model:", gsub("^WINNER_", "", gsub("\\.rds$", "", model_name)))
        )

    return(combined)
}

### 9. JSDM-SPECIFIC PDP FUNCTIONS (from v5 Original) -----------------------------

#' Validate JSDM species order and structure
#' @param model brmsfit object
#' @return Validation list with species info
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

#' PDPs JSDM: Trade-offs ecológicos (composição relativa)
#'
#' @param model Objeto brmsfit JSDM
#' @param model_name Nome do modelo
#' @param ndraws Número de draws (padrão: 500)
#' @param resolution Número de pontos (padrão: 100)
#' @return Objeto ggplot2 combinado ou NULL
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
        theme(legend.position = "right")

    # Combinar os dois plots
    combined <- (p1 | p2) +
        plot_annotation(
            title = "JSDM: Ecological Trade-offs in Community Composition",
            subtitle = paste("Model:", model_name)
        )

    return(combined)
}

### 10. DIAGNOSTIC PLOTS --------------------------------------------------------

generate_diagnostics_year_re <- function(model, model_name) {
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
                                            ifelse(max_rhat < 1.01, "(OK)", "(FAIL)"))))
}

### 10. MAIN PROCESSING PIPELINE -----------------------------------------------

process_model_year_re <- function(rds_file, output_base_dir) {
    model <- readRDS(rds_file)
    model_name <- gsub("\\.rds$", "", basename(rds_file))

    model_output_dir <- file.path(output_base_dir, model_name)
    if (!dir.exists(model_output_dir)) {
        dir.create(model_output_dir, recursive = TRUE)
    }

    cat(sprintf("\n🔍 PROCESSING: %s\n", model_name))

    model_info <- detect_model_type_year_re(model, model_name)

    # Specialized pipeline for DIRICHLET (JSDM) models
    if (model_info$type == "DIRICHLET") {
        cat("\n  🔄 DIRICHLET (JSDM) MODEL DETECTED - Using specialized JSDM pipeline\n")

        # JSDM PDP 1: Community Overlaid (curvas separadas)
        cat("\n  [1/4] JSDM PDP: Community Overlaid (curvas separadas)...\n")
        p_jsdm_community <- generate_jsdm_pdp_community_overlaid_CORRECTED(model, model_name)
        if (!is.null(p_jsdm_community)) {
            ggsave(file.path(model_output_dir, "FIGURE_1_JSDM_Community_Overlaid.png"),
                   p_jsdm_community, width = 10, height = 7, dpi = 600, bg = "white")
            ggsave(file.path(model_output_dir, "FIGURE_1_JSDM_Community_Overlaid.pdf"),
                   p_jsdm_community, width = 10, height = 7, bg = "white")
            cat("    ✓ FIGURE_1_JSDM_Community_Overlaid salvo\n")
        }

        # JSDM PDP 2: Per-Species (facets)
        cat("\n  [2/4] JSDM PDP: Per-Species (facets)...\n")
        p_jsdm_species <- generate_jsdm_pdp_per_species_CORRECTED(model, model_name)
        if (!is.null(p_jsdm_species)) {
            ggsave(file.path(model_output_dir, "FIGURE_2_JSDM_Per_Species.png"),
                   p_jsdm_species, width = 14, height = 10, dpi = 600, bg = "white")
            ggsave(file.path(model_output_dir, "FIGURE_2_JSDM_Per_Species.pdf"),
                   p_jsdm_species, width = 14, height = 10, bg = "white")
            cat("    ✓ FIGURE_2_JSDM_Per_Species salvo\n")
        }

        # JSDM PDP 3: Composition Tradeoff
        cat("\n  [3/4] JSDM PDP: Composition Tradeoff (trade-offs ecológicos)...\n")
        p_jsdm_tradeoff <- generate_jsdm_composition_tradeoff(model, model_name)
        if (!is.null(p_jsdm_tradeoff)) {
            ggsave(file.path(model_output_dir, "FIGURE_3_JSDM_Composition_Tradeoff.png"),
                   p_jsdm_tradeoff, width = 16, height = 10, dpi = 600, bg = "white")
            ggsave(file.path(model_output_dir, "FIGURE_3_JSDM_Composition_Tradeoff.pdf"),
                   p_jsdm_tradeoff, width = 16, height = 10, bg = "white")
            cat("    ✓ FIGURE_3_JSDM_Composition_Tradeoff salvo\n")
        }

        # JSDM: Diagnostics (shared)
        cat("\n  [4/4] JSDM Diagnostics...\n")
        p_diag <- generate_diagnostics_year_re(model, model_name)
        if (!is.null(p_diag)) {
            ggsave(file.path(model_output_dir, "FIGURE_4_Diagnostics.png"),
                   p_diag, width = 14, height = 10, dpi = 600, bg = "white")
            ggsave(file.path(model_output_dir, "FIGURE_4_Diagnostics.pdf"),
                   p_diag, width = 14, height = 10, bg = "white")
            cat("    ✓ FIGURE_4_Diagnostics salvo\n")
        }

        # YEAR RE Diagnostic if applicable
        if (model_info$has_year_re) {
            cat("\n  [EXTRA] YEAR RE Diagnostics...\n")
            p_year <- generate_year_re_diagnostics(model)
            if (!is.null(p_year)) {
                ggsave(file.path(model_output_dir, "FIGURE_5_YEAR_RE_Diagnostics.png"),
                       p_year, width = 8, height = 6, dpi = 600, bg = "white")
                ggsave(file.path(model_output_dir, "FIGURE_5_YEAR_RE_Diagnostics.pdf"),
                       p_year, width = 8, height = 6, bg = "white")
                cat("    ✓ FIGURE_5_YEAR_RE_Diagnostics salvo\n")
            }
        }

        cat(sprintf("\n✅ MODELO JSDM CONCLUÍDO: %s\n", model_name))
        return(invisible(NULL))
    }

    # Standard pipeline for GAUSSIAN and ZOIB models
    # 1. Forest Plot
    cat("\n  [1/5] Forest Plot...\n")
    p_forest <- generate_forest_plot_year_re(model, model_name, model_info)
    if (!is.null(p_forest)) {
        ggsave(file.path(model_output_dir, "FIGURE_1_Forest_Plot.png"),
               p_forest, width = 10, height = 6, dpi = 600, bg = "white")
        ggsave(file.path(model_output_dir, "FIGURE_1_Forest_Plot.pdf"),
               p_forest, width = 10, height = 6, bg = "white")
        cat("    ✓ FIGURE_1 salvo\n")
    }

    # 2. Marginal Effects
    cat("\n  [2/5] Marginal Effects...\n")
    p_effects <- generate_marginal_effects_year_re(model, model_name, model_info)
    if (!is.null(p_effects)) {
        ggsave(file.path(model_output_dir, "FIGURE_2_Marginal_Effects.png"),
               p_effects, width = 12, height = 10, dpi = 600, bg = "white")
        ggsave(file.path(model_output_dir, "FIGURE_2_Marginal_Effects.pdf"),
               p_effects, width = 12, height = 10, bg = "white")
        cat("    ✓ FIGURE_2 salvo\n")
    }

    # 3. PPC
    cat("\n  [3/5] PPC...\n")
    p_ppc <- generate_ppc_plots_year_re(model, model_name, model_info)
    if (!is.null(p_ppc)) {
        ggsave(file.path(model_output_dir, "FIGURE_3_PPC.png"),
               p_ppc, width = 12, height = 10, dpi = 600, bg = "white")
        ggsave(file.path(model_output_dir, "FIGURE_3_PPC.pdf"),
               p_ppc, width = 12, height = 10, bg = "white")
        cat("    ✓ FIGURE_3 salvo\n")
    }

    # 4. Diagnostics
    cat("\n  [4/5] Diagnostics...\n")
    p_diag <- generate_diagnostics_year_re(model, model_name)
    if (!is.null(p_diag)) {
        ggsave(file.path(model_output_dir, "FIGURE_4_Diagnostics.png"),
               p_diag, width = 14, height = 10, dpi = 600, bg = "white")
        ggsave(file.path(model_output_dir, "FIGURE_4_Diagnostics.pdf"),
               p_diag, width = 14, height = 10, bg = "white")
        cat("    ✓ FIGURE_4 salvo\n")
    }

    # 5. YEAR RE Diagnostic (NEW!)
    cat("\n  [5/5] YEAR RE Diagnostics...\n")
    if (model_info$has_year_re) {
        p_year <- generate_year_re_diagnostics(model)
        if (!is.null(p_year)) {
            ggsave(file.path(model_output_dir, "FIGURE_5_YEAR_RE_Diagnostics.png"),
                   p_year, width = 8, height = 6, dpi = 600, bg = "white")
            ggsave(file.path(model_output_dir, "FIGURE_5_YEAR_RE_Diagnostics.pdf"),
                   p_year, width = 8, height = 6, bg = "white")
            cat("    ✓ FIGURE_5 salvo (YEAR RE)\n")
        }
    } else {
        cat("    ⚠ Modelo não possui YEAR random effect\n")
    }

    cat(sprintf("\n✅ MODELO CONCLUÍDO: %s\n", model_name))
}

### 11. EXECUÇÃO PRINCIPAL =======================================================

cat("\n", rep("=", 75), "\n", sep = "")
cat("🚀 INICIANDO PIPELINE V5 - YEAR RANDOM EFFECT MODELS\n")
cat(rep("=", 75), "\n", sep = "")

winner_files <- discover_winner_models(search_directories_YEAR_RE)

processed_info <- list()
for (rds_file in winner_files) {
    result <- process_model_year_re(rds_file, output_dir_YEAR_RE)
    if (!is.null(result)) processed_info[[length(processed_info) + 1]] <- result
}

cat("\n", rep("=", 75), "\n", sep = "")
cat("📋 RELATÓRIO FINAL - YEAR RE PIPELINE V5\n")
cat(rep("=", 75), "\n", sep = "")

cat(sprintf("\n✓ Modelos processados: %d\n", length(processed_info)))
cat(sprintf("📁 Saída: %s\n", output_dir_YEAR_RE))

cat("\nFiguras geradas:\n")
cat("  [ZOIB/GAUSSIAN MODELS:]\n")
cat("  - FIGURE_1_Forest_Plot.png/pdf (com YEAR SD annotation quando aplicável)\n")
cat("  - FIGURE_2_Marginal_Effects.png/pdf\n")
cat("  - FIGURE_3_PPC.png/pdf\n")
cat("  - FIGURE_4_Diagnostics.png/pdf\n")
cat("  - FIGURE_5_YEAR_RE_Diagnostics.png/pdf (apenas para modelos com YEAR RE)\n")
cat("\n  [JSDM/DIRICHLET MODELS:]\n")
cat("  - FIGURE_1_JSDM_Community_Overlaid.png/pdf (curvas separadas por espécie)\n")
cat("  - FIGURE_2_JSDM_Per_Species.png/pdf (facets por espécie)\n")
cat("  - FIGURE_3_JSDM_Composition_Tradeoff.png/pdf (trade-offs ecológicos)\n")
cat("  - FIGURE_4_Diagnostics.png/pdf\n")
cat("  - FIGURE_5_YEAR_RE_Diagnostics.png/pdf (apenas para modelos com YEAR RE)\n")

cat("\n🎉 YEAR RE PIPELINE V5 CONCLUÍDO!\n\n")
