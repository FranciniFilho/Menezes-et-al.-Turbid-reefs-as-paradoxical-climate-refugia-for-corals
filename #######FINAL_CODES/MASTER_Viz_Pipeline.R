### ==================================================================================== ###
### === MASTER VISUALIZATION PIPELINE - FIGURAS DE PUBLICAÇÃO (AUTOMÁTICO)           === ###
### === Varre diretórios, detecta modelos vencedores e gera visualizações Q1         === ###
### ==================================================================================== ###

rm(list = ls())
gc()

### 1. PACOTES NECESSÁRIOS --------------------------------------------------------
libs <- c("brms", "ggplot2", "dplyr", "tidybayes", "ggdist", "bayesplot", "patchwork", "cowplot")
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
            panel.grid = element_blank() # Sem gridlines
        )
}
theme_set(theme_publication())

# Paleta Okabe-Ito (acessível para daltonismo)
okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
color_primary <- "#0072B2"
color_secondary <- "#E69F00"

### 3. CONFIGURAÇÃO DOS DIRETÓRIOS DE BUSCA ---------------------------------------
# Diretórios onde os scripts de modelagem salvam os modelos vencedores
search_directories <- c(
    "C:/Users/rbfra/OneDrive/Bayesian_Analyses_ZOIB_Abundance_MULTI_CV",
    "C:/Users/rbfra/OneDrive/Bayesian_Analyses_Health_Growth_MultiDataset"
)

# Diretório de saída para figuras de publicação
output_dir <- "C:/Users/rbfra/OneDrive/Bayesian_Figures_Publication"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n", rep("=", 75), "\n", sep = "")
cat("🎨 MASTER VISUALIZATION PIPELINE - FIGURAS DE PUBLICAÇÃO\n")
cat(rep("=", 75), "\n", sep = "")
cat("📁 Diretório de saída:", output_dir, "\n")

### 4. AUTO-DESCOBERTA DE MODELOS VENCEDORES (.rds com WINNER) --------------------
discover_winner_models <- function(search_dirs) {
    cat("\n--- Buscando modelos vencedores (.rds) ---\n")

    all_rds_files <- c()
    for (dir_path in search_dirs) {
        if (dir.exists(dir_path)) {
            files <- list.files(
                path = dir_path,
                pattern = "WINNER.*\\.rds$",
                recursive = TRUE,
                full.names = TRUE,
                ignore.case = TRUE
            )
            all_rds_files <- c(all_rds_files, files)
            cat(sprintf("  ✓ %s: %d modelo(s) encontrado(s)\n", basename(dir_path), length(files)))
        } else {
            cat(sprintf("  ⚠ Diretório não existe: %s\n", dir_path))
        }
    }

    if (length(all_rds_files) == 0) {
        stop("ERRO: Nenhum modelo vencedor encontrado. Execute os scripts de modelagem primeiro.")
    }

    cat(sprintf("\n📊 Total de modelos vencedores encontrados: %d\n", length(all_rds_files)))
    return(all_rds_files)
}

### 5. DETECÇÃO AUTOMÁTICA DO TIPO DE MODELO (ZOIB vs GAUSSIANO) ------------------
detect_model_type <- function(model) {
    fam <- family(model)$family
    if (grepl("zero_one_inflated_beta", fam, ignore.case = TRUE)) {
        return("ZOIB")
    } else if (grepl("gaussian", fam, ignore.case = TRUE)) {
        return("GAUSSIAN")
    } else {
        return("OTHER")
    }
}

### 6. GERAÇÃO DE FOREST PLOT (HALF-EYE) ------------------------------------------
generate_forest_plot <- function(model, model_name) {
    cat("    -> Gerando Forest Plot (Half-Eye)...\n")

    # Extrair variáveis de efeitos fixos (sem intercepto e sds)
    all_vars <- get_variables(model)
    fixed_vars <- all_vars[grepl("^b_", all_vars)]
    fixed_vars <- fixed_vars[!grepl("Intercept|sds_|^b_sigma|^b_phi|^b_zoi|^b_coi", fixed_vars)]

    if (length(fixed_vars) == 0) {
        cat("    ⚠ Modelo sem efeitos fixos lineares (pode ser dominado por splines)\n")
        return(NULL)
    }

    # Extrair draws para todas as variáveis
    draws_df <- as_draws_df(model)

    # Filtrar apenas as variáveis de interesse
    vars_in_draws <- intersect(fixed_vars, names(draws_df))
    if (length(vars_in_draws) == 0) {
        return(NULL)
    }

    # Reformatar para tidybayes format
    draws_long <- draws_df %>%
        dplyr::select(all_of(vars_in_draws)) %>%
        tidyr::pivot_longer(
            cols = everything(),
            names_to = "parameter",
            values_to = "value"
        ) %>%
        mutate(parameter = gsub("^b_", "", parameter)) # Limpar nomes

    # Criar plot
    p <- ggplot(draws_long, aes(y = reorder(parameter, abs(value)), x = value, fill = after_stat(x > 0))) +
        stat_halfeye(
            .width = c(0.89, 0.95),
            alpha = 0.85,
            point_interval = median_qi,
            point_size = 3
        ) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.7) +
        scale_fill_manual(values = c("TRUE" = color_primary, "FALSE" = "gray60"), guide = "none") +
        labs(
            title = "Magnitude of Fixed Effects",
            subtitle = paste("Model:", gsub("_", " ", gsub("WINNER_", "", model_name))),
            x = "Effect Estimate (Standardized Scale)",
            y = NULL,
            caption = "Points: Median | Bars: 89% and 95% CI"
        ) +
        theme_publication()

    return(p)
}

### 7. GERAÇÃO DE EFEITOS MARGINAIS (RIBBON PLOTS) --------------------------------
generate_marginal_effects <- function(model, model_name) {
    cat("    -> Gerando Efeitos Marginais (Ribbon/Point)...\n")

    ce <- tryCatch(
        {
            conditional_effects(model, robust = TRUE, method = "posterior_epred")
        },
        error = function(e) {
            cat(sprintf("    ⚠ Erro ao calcular efeitos: %s\n", e$message))
            return(NULL)
        }
    )

    if (is.null(ce) || length(ce) == 0) {
        return(NULL)
    }

    # Cores padronizadas para ARC (Inner/Outer)
    arc_colors <- c("Inner" = "#0072B2", "Outer" = "#E69F00")

    plot_list <- list()

    for (eff_name in names(ce)) {
        data_eff <- ce[[eff_name]]
        x_var <- names(data_eff)[1]

        # Detectar se a variável X é categórica ou numérica
        x_is_categorical <- is.factor(data_eff[[x_var]]) || is.character(data_eff[[x_var]])

        # Detectar coluna ARCH (pode variar dependendo do tipo de modelo)
        arch_col <- NULL
        if ("ARCH" %in% names(data_eff)) {
            arch_col <- "ARCH"
        } else if ("effect2__" %in% names(data_eff)) {
            arch_col <- "effect2__"
        } else if ("cond__" %in% names(data_eff)) {
            arch_col <- "cond__"
        }

        # Verificar se este é um efeito de interação com ARCH
        has_arch_interaction <- !is.null(arch_col) && grepl(":ARCH|ARCH:", eff_name)

        if (x_is_categorical) {
            # Para variáveis categóricas: usar pointrange (pontos com barras de erro)
            if (has_arch_interaction && length(unique(data_eff[[arch_col]])) > 1) {
                # Com interação ARCH - cores separadas
                p <- ggplot(data_eff, aes(x = .data[[x_var]], y = estimate__,
                                          color = .data[[arch_col]], group = .data[[arch_col]])) +
                    geom_pointrange(aes(ymin = lower__, ymax = upper__), size = 1, linewidth = 1.2,
                                    position = position_dodge(width = 0.5)) +
                    scale_color_manual(values = arc_colors, name = "ARC") +
                    labs(
                        title = paste("Effect of", gsub("_scaled", "", x_var), "by ARC"),
                        y = "Prediction",
                        x = gsub("_scaled", "", x_var)
                    ) +
                    theme_publication() +
                    theme(legend.position = "bottom")
            } else {
                p <- ggplot(data_eff, aes(x = .data[[x_var]], y = estimate__)) +
                    geom_pointrange(
                        aes(ymin = lower__, ymax = upper__),
                        color = color_primary,
                        size = 1,
                        linewidth = 1.2
                    ) +
                    labs(
                        title = paste("Effect of", gsub("_scaled", "", x_var)),
                        y = "Prediction",
                        x = gsub("_scaled", "", x_var)
                    ) +
                    theme_publication()
            }
        } else {
            # Para variáveis contínuas: usar ribbon + line
            if (has_arch_interaction && length(unique(data_eff[[arch_col]])) > 1) {
                # Com interação ARCH - linhas sobrepostas com cores diferentes
                cat(sprintf("      → Plotando interação %s com cores por ARC (coluna: %s)\n", eff_name, arch_col))
                p <- ggplot(data_eff, aes(x = .data[[x_var]], y = estimate__,
                                          color = .data[[arch_col]], fill = .data[[arch_col]],
                                          group = .data[[arch_col]])) +
                    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, linewidth = 0) +
                    geom_line(linewidth = 1.2) +
                    scale_color_manual(values = arc_colors, name = "ARC") +
                    scale_fill_manual(values = arc_colors, name = "ARC") +
                    labs(
                        title = paste("Effect of", gsub("_scaled", "", x_var), "by ARC"),
                        y = "Prediction",
                        x = gsub("_scaled", "", x_var)
                    ) +
                    theme_publication() +
                    theme(legend.position = "bottom")
            } else {
                # Efeito principal - cor única
                p <- ggplot(data_eff, aes(x = .data[[x_var]], y = estimate__)) +
                    geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = color_primary, alpha = 0.25) +
                    geom_line(color = color_primary, linewidth = 1.2) +
                    labs(
                        title = paste("Effect of", gsub("_scaled", "", x_var)),
                        y = "Prediction",
                        x = gsub("_scaled", "", x_var)
                    ) +
                    theme_publication()
            }
        }

        plot_list[[eff_name]] <- p
    }

    # Combinar plots
    n_plots <- length(plot_list)
    ncol_calc <- min(3, n_plots)

    combined <- wrap_plots(plot_list, ncol = ncol_calc) +
        plot_annotation(
            title = "Conditional Marginal Effects",
            subtitle = paste("Model:", gsub("_", " ", gsub("WINNER_", "", model_name))),
            caption = "Line/Point: Median | Band/Bar: 95% CI | Colored lines indicate ARC interaction",
            theme = theme(
                plot.title = element_text(face = "bold", size = 18),
                plot.subtitle = element_text(size = 12)
            )
        )

    return(combined)
}

### 8. VALIDAÇÃO PREDITIVA (PPC) --------------------------------------------------
generate_ppc_plots <- function(model, model_name, model_type) {
    cat("    -> Gerando Validação Preditiva (PPC)...\n")

    response_label <- if (model_type == "ZOIB") "Cover (Proportion)" else "Value"

    # Densidade
    ppc_dens <- pp_check(model, ndraws = 100, type = "dens_overlay") +
        scale_color_manual(values = c("y" = "black", "yrep" = color_primary)) +
        labs(
            title = "Density Validation",
            subtitle = "Black: Observed | Blue: Simulated",
            x = response_label,
            y = "Density"
        ) +
        theme_publication() +
        theme(legend.position = "none")

    # Para ZOIB: adicionar validação de zeros
    if (model_type == "ZOIB") {
        ppc_zeros <- tryCatch(
            {
                pp_check(model, ndraws = 100, type = "stat", stat = function(y) mean(y == 0)) +
                    labs(
                        title = "Zero Validation",
                        subtitle = "Black line should be under blue distribution",
                        x = "Proportion of Zeros",
                        y = "Frequency"
                    ) +
                    theme_publication() +
                    theme(legend.position = "none")
            },
            error = function(e) NULL
        )

        if (!is.null(ppc_zeros)) {
            combined_ppc <- (ppc_dens | ppc_zeros) +
                plot_annotation(
                    title = "Posterior Predictive Checks (ZOIB)",
                    subtitle = paste("Modelo:", gsub("_", " ", gsub("WINNER_", "", model_name)))
                )
            return(combined_ppc)
        }
    }

    # Para Gaussiano: ECDF
    ppc_ecdf <- pp_check(model, ndraws = 100, type = "ecdf_overlay") +
        scale_color_manual(values = c("y" = "black", "yrep" = color_primary)) +
        labs(
            title = "Cumulative ECDF",
            subtitle = "Cumulative distribution fit",
            x = response_label,
            y = "Cumulative Probability"
        ) +
        theme_publication() +
        theme(legend.position = "none")

    combined_ppc <- (ppc_dens | ppc_ecdf) +
        plot_annotation(
            title = "Posterior Predictive Checks",
            subtitle = paste("Modelo:", gsub("_", " ", gsub("WINNER_", "", model_name)))
        )

    return(combined_ppc)
}

### 9. DIAGNÓSTICOS DE CONVERGÊNCIA -----------------------------------------------
generate_diagnostics <- function(model, model_name) {
    cat("    -> Gerando Diagnósticos de Convergência...\n")

    # Selecionar parâmetros para trace plots (máximo 4)
    all_vars <- variables(model)
    params_to_trace <- head(all_vars[grepl("^b_", all_vars)], 4)

    if (length(params_to_trace) < 2) {
        params_to_trace <- head(all_vars, 4)
    }

    # Trace Plot
    p_trace <- mcmc_trace(model, pars = params_to_trace) +
        theme_publication() +
        scale_color_manual(values = okabe_ito[1:4]) +
        labs(
            title = "Trace Plots (Convergence)",
            subtitle = "Chains should look like well-mixed 'hairy caterpillars'"
        )

    # Autocorrelação
    p_acf <- mcmc_acf(model, pars = params_to_trace) +
        theme_publication() +
        labs(
            title = "Autocorrelation",
            subtitle = "Should decay rapidly to zero"
        )

    combined <- (p_trace / p_acf) +
        plot_annotation(
            title = "MCMC Convergence Diagnostics",
            subtitle = paste("Model:", gsub("_", " ", gsub("WINNER_", "", model_name)))
        )

    return(combined)
}

### 10. PIPELINE PRINCIPAL --------------------------------------------------------
process_model <- function(rds_path, output_base_dir) {
    cat("\n", rep("-", 60), "\n", sep = "")

    model_name <- tools::file_path_sans_ext(basename(rds_path))
    cat(sprintf("📊 Processando: %s\n", model_name))

    # Carregar modelo
    cat("    -> Carregando modelo...\n")
    model <- tryCatch(
        {
            readRDS(rds_path)
        },
        error = function(e) {
            cat(sprintf("    ❌ Erro ao carregar: %s\n", e$message))
            return(NULL)
        }
    )

    if (is.null(model)) {
        return(NULL)
    }

    # Detectar tipo
    model_type <- detect_model_type(model)
    cat(sprintf("    -> Tipo detectado: %s\n", model_type))

    # Criar diretório de saída para este modelo
    model_output_dir <- file.path(output_base_dir, model_name)
    dir.create(model_output_dir, showWarnings = FALSE, recursive = TRUE)

    # 1. Forest Plot
    p_forest <- generate_forest_plot(model, model_name)
    if (!is.null(p_forest)) {
        ggsave(
            file.path(model_output_dir, "FIGURE_1_Forest_Plot.png"),
            p_forest,
            width = 10, height = 6, dpi = 300, bg = "white"
        )
        cat("    ✓ Forest Plot salvo\n")
    }

    # 2. Efeitos Marginais
    p_effects <- generate_marginal_effects(model, model_name)
    if (!is.null(p_effects)) {
        n_effects <- length(conditional_effects(model))
        h_calc <- ceiling(n_effects / 3) * 4.5
        ggsave(
            file.path(model_output_dir, "FIGURE_2_Marginal_Effects.png"),
            p_effects,
            width = 14, height = max(6, h_calc), dpi = 300, bg = "white"
        )
        cat("    ✓ Efeitos Marginais salvos\n")
    }

    # 3. PPC (Validação Preditiva)
    p_ppc <- generate_ppc_plots(model, model_name, model_type)
    if (!is.null(p_ppc)) {
        ggsave(
            file.path(model_output_dir, "FIGURE_3_Validation_PPC.png"),
            p_ppc,
            width = 12, height = 6, dpi = 300, bg = "white"
        )
        cat("    ✓ Validação PPC salva\n")
    }

    # 4. Diagnósticos
    p_diag <- generate_diagnostics(model, model_name)
    if (!is.null(p_diag)) {
        ggsave(
            file.path(model_output_dir, "FIGURE_4_Diagnostics.png"),
            p_diag,
            width = 10, height = 10, dpi = 300, bg = "white"
        )
        cat("    ✓ Diagnósticos salvos\n")
    }

    cat(sprintf("✅ Concluído: %s\n", model_name))
    return(model_output_dir)
}

### 11. EXECUÇÃO PRINCIPAL --------------------------------------------------------
cat("\n", rep("=", 75), "\n", sep = "")
cat("🚀 INICIANDO PIPELINE DE VISUALIZAÇÃO\n")
cat(rep("=", 75), "\n", sep = "")

# Descobrir modelos
winner_files <- discover_winner_models(search_directories)

# Processar cada modelo
processed_dirs <- c()
for (rds_file in winner_files) {
    result_dir <- process_model(rds_file, output_dir)
    if (!is.null(result_dir)) {
        processed_dirs <- c(processed_dirs, result_dir)
    }
}

### 12. RELATÓRIO FINAL -----------------------------------------------------------
cat("\n", rep("=", 75), "\n", sep = "")
cat("📋 RELATÓRIO FINAL\n")
cat(rep("=", 75), "\n")
cat(sprintf("✓ Modelos processados: %d\n", length(processed_dirs)))
cat(sprintf("📁 Figuras salvas em: %s\n", output_dir))
cat("\nDiretórios de saída:\n")
for (d in processed_dirs) {
    cat(sprintf("  • %s\n", d))
}
cat("\n🎉 PIPELINE CONCLUÍDO COM SUCESSO!\n")
