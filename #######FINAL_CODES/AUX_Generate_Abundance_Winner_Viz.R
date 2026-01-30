### ==================================================================================== ###
### === AUXILIARY SCRIPT: GERAR FIGURAS DO VENCEDOR DE ABUNDÂNCIA (FROM CACHE)      === ###
### === Recupera o modelo do cache e gera as figuras faltantes sem re-rodar tudo    === ###
### ==================================================================================== ###
### === SOLUÇÃO: Usar diretório de saída MASTER para evitar problemas de caminho  === ###
### ==================================================================================== ###

rm(list = ls())
gc()

### 1. CONFIGURAÇÃO E PACOTES --------------------------------------------------------
libs <- c("brms", "ggplot2", "dplyr", "tidybayes", "ggdist", "bayesplot", "patchwork", "cowplot")
invisible(lapply(libs, library, character.only = TRUE))

# --- CAMINHOS: Usar diretório com caminho curto (igual MASTER_Viz_Pipeline) ---
# O MASTER_Viz_Pipeline funciona porque os arquivos WINNER estão salvos em:
# C:/Users/rbfra/OneDrive/Bayesian_Analyses_ZOIB_Abundance_MULTI_CV
# que tem um caminho MUITO mais curto!

search_dir <- "C:/Users/rbfra/OneDrive/Bayesian_Analyses_ZOIB_Abundance_MULTI_CV"
output_dir <- "C:/Users/rbfra/OneDrive/Bayesian_Figures_Publication"

# --- MODELO VENCEDOR (Especificado pelo usuário) ---
winner_id <- "PCA_Informative_model_interaction_arch_no_depth_with_hab"
cached_filename <- paste0("ZOIB_", winner_id, ".rds")

cat(sprintf("🔍 Buscando modelo: %s\n", cached_filename))
cat(sprintf("   Diretório: %s\n", search_dir))

# Buscar o arquivo recursivamente (igual MASTER_Viz_Pipeline)
cached_files <- list.files(
  path = search_dir,
  pattern = cached_filename,
  full.names = TRUE,
  recursive = TRUE
)

if (length(cached_files) == 0) {
  stop("ERRO CRÍTICO: Arquivo do modelo não encontrado! Verifique se os modelos foram salvos com WINNER_ no nome.")
}

cached_path <- cached_files[1]
cat(sprintf("   ✓ Arquivo encontrado!\n"))

### 2. CARREGAR MODELO ---------------------------------------------------------------
cat("📦 Carregando modelo (pode demorar alguns segundos)...\n")
model <- readRDS(cached_path)

### 3. TEMA E FUNÇÕES DE VISUALIZAÇÃO (NATURA/SCIENCE STYLE) -------------------------
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

okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")
color_primary <- "#0072B2"

# Cores padronizadas para ARC (Inner/Outer) - usar em todos os plots
arc_colors <- c("Inner" = "#0072B2", "Outer" = "#E69F00")

### 4. GERAÇÃO DAS FIGURAS -----------------------------------------------------------

model_name <- paste0("WINNER_ZOIB_", winner_id)
model_fig_dir <- file.path(output_dir, model_name)
dir.create(model_fig_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("🎨 Gerando figuras em:\n  %s\n", model_fig_dir))

# --- A. Forest Plot (Half-Eye) ---
cat("   -> Gerando Forest Plot...\n")
all_vars <- get_variables(model)
fixed_vars <- all_vars[grepl("^b_", all_vars)]
fixed_vars <- fixed_vars[!grepl("Intercept|sds_|^b_sigma|^b_phi|^b_zoi|^b_coi", fixed_vars)]

if (length(fixed_vars) > 0) {
    draws_df <- as_draws_df(model)
    vars_in_draws <- intersect(fixed_vars, names(draws_df))

    if (length(vars_in_draws) > 0) {
        draws_long <- draws_df %>%
            dplyr::select(all_of(vars_in_draws)) %>%
            tidyr::pivot_longer(cols = everything(), names_to = "parameter", values_to = "value") %>%
            mutate(parameter = gsub("^b_", "", parameter))

        p_forest <- ggplot(draws_long, aes(y = reorder(parameter, abs(value)), x = value, fill = after_stat(x > 0))) +
            stat_halfeye(.width = c(0.89, 0.95), alpha = 0.85, point_interval = median_qi, point_size = 3) +
            geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.7) +
            scale_fill_manual(values = c("TRUE" = color_primary, "FALSE" = "gray60"), guide = "none") +
            labs(
                title = "Magnitude of Fixed Effects (Abundance)",
                subtitle = paste("Model:", winner_id),
                x = "Effect Estimate (Standardized Scale)", y = NULL
            ) + theme_publication()

        ggsave(file.path(model_fig_dir, "FIGURE_1_Forest_Plot.png"), p_forest, width = 10, height = 6, dpi = 300, bg = "white")
        cat("   ✓ Forest Plot salvo\n")
    }
} else {
    cat("   ⚠ Aviso: Modelo dominado por splines, sem efeitos fixos lineares para Forest Plot.\n")
}

# --- B. Marginal Effects (Ribbon/Point) ---
cat("   -> Gerando Efeitos Marginais...\n")
ce <- try(conditional_effects(model, robust = TRUE, method = "posterior_epred"), silent = TRUE)

if (!inherits(ce, "try-error") && !is.null(ce)) {
    plot_list <- list()

    # FILTRAR: Remover efeitos principais redundantes quando existe interação correspondente
    # Ex: Se existe "PC1_VARIABILITY:ARCH", remover "PC1_VARIABILITY" sozinho
    effect_names <- names(ce)

    # Identificar quais efeitos têm correspondente com interação ARCH
    to_remove <- c()
    for (name in effect_names) {
        # Se este efeito tem ":ARCH", adicione a versão sem ":ARCH" à lista de remoção
        if (grepl(":ARCH", name)) {
            base_name <- gsub(":ARCH", "", name)
            to_remove <- c(to_remove, base_name)
        }
    }

    # Remover efeitos redundantes
    if (length(to_remove) > 0) {
        ce <- ce[!names(ce) %in% to_remove]
        cat(sprintf("   → Removidos %d efeito(s) principal(is) redundante(s)\n", length(to_remove)))
    }

    for (eff_name in names(ce)) {
        data_eff <- ce[[eff_name]]
        x_var <- names(data_eff)[1]
        x_is_cat <- is.factor(data_eff[[x_var]]) || is.character(data_eff[[x_var]])

        # --- CORREÇÃO CRÍTICA: Detectar coluna ARCH corretamente ---
        # conditional_effects() pode usar "ARCH", "effect2__", "cond__" ou outras colunas
        # dependendo do tipo de spline e da estrutura do modelo
        arch_col <- NULL
        if ("ARCH" %in% names(data_eff)) {
            arch_col <- "ARCH"
        } else if ("effect2__" %in% names(data_eff)) {
            # Para splines com by=ARCH, a segunda variável vai para effect2__
            arch_col <- "effect2__"
        } else if ("cond__" %in% names(data_eff)) {
            arch_col <- "cond__"
        }

        # Verificar se este é um efeito de interação com ARCH
        has_arch_interaction <- !is.null(arch_col) && grepl(":ARCH|ARCH:", eff_name)

        if (x_is_cat) {
            if (has_arch_interaction && !is.null(arch_col) && length(unique(data_eff[[arch_col]])) > 1) {
                # --- CATEGÓRICA COM INTERAÇÃO ARCH: Cores separadas para Inner/Outer ---
                cat(sprintf("      → Plotando interação categórica %s com cores por ARC (coluna: %s)\n", eff_name, arch_col))
                p <- ggplot(data_eff, aes(x = .data[[x_var]], y = estimate__,
                                          color = .data[[arch_col]], group = .data[[arch_col]])) +
                    geom_pointrange(aes(ymin = lower__, ymax = upper__), size = 1, linewidth = 1.2,
                                    position = position_dodge(width = 0.5)) +
                    scale_color_manual(values = arc_colors, name = "ARC") +
                    labs(title = paste("Effect of", gsub("_scaled", "", x_var), "by ARC"),
                         y = "Predicted Coverage (Prop)", x = NULL) +
                    theme_publication() +
                    theme(legend.position = "bottom")
            } else {
                p <- ggplot(data_eff, aes(x = .data[[x_var]], y = estimate__)) +
                    geom_pointrange(aes(ymin = lower__, ymax = upper__), color = color_primary, size = 1, linewidth = 1.2) +
                    labs(title = paste("Effect of", gsub("_scaled", "", x_var)), y = "Predicted Coverage (Prop)", x = NULL) +
                    theme_publication()
            }
        } else {
            if (has_arch_interaction && !is.null(arch_col) && length(unique(data_eff[[arch_col]])) > 1) {
                # --- EFEITO DE INTERAÇÃO: Cores separadas para Inner/Outer ---
                cat(sprintf("      → Plotando interação %s com cores por ARC (coluna: %s)\n", eff_name, arch_col))

                # Normalizar valores da coluna ARCH para garantir correspondência correta das cores
                # Se a coluna contiver fatores ou caracteres, mapear para "Inner" e "Outer"
                data_eff_plot <- data_eff
                if (is.factor(data_eff_plot[[arch_col]]) || is.character(data_eff_plot[[arch_col]])) {
                    # Verificar os valores únicos e criar um mapeamento
                    unique_vals <- unique(data_eff_plot[[arch_col]])
                    # Se os valores já são "Inner" e "Outer", usar diretamente
                    if (all(unique_vals %in% c("Inner", "Outer"))) {
                        # Já está no formato correto
                    } else {
                        # Tentar mapear baseado em algum critério (ex: "inner"/"outer" case-insensitive)
                        data_eff_plot[[arch_col]] <- ifelse(
                            grepl("inner", data_eff_plot[[arch_col]], ignore.case = TRUE),
                            "Inner",
                            ifelse(grepl("outer", data_eff_plot[[arch_col]], ignore.case = TRUE),
                                   "Outer", as.character(data_eff_plot[[arch_col]]))
                        )
                    }
                }

                p <- ggplot(data_eff_plot, aes(x = .data[[x_var]], y = estimate__,
                                                color = .data[[arch_col]], fill = .data[[arch_col]],
                                                group = .data[[arch_col]])) +
                    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.2, linewidth = 0, color = NA) +
                    geom_line(aes(color = .data[[arch_col]]), linewidth = 1.2) +
                    scale_color_manual(values = arc_colors, name = "ARC", drop = FALSE) +
                    scale_fill_manual(values = arc_colors, name = "ARC", drop = FALSE) +
                    labs(
                        title = paste("Effect of", gsub("_scaled", "", x_var), "by ARC"),
                        y = "Predicted Coverage (Prop)",
                        x = gsub("_scaled", "", x_var)
                    ) +
                    theme_publication() +
                    theme(legend.position = "bottom")
            } else {
                # --- EFEITO PRINCIPAL: Cor única ---
                p <- ggplot(data_eff, aes(x = .data[[x_var]], y = estimate__)) +
                    geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = color_primary, alpha = 0.25) +
                    geom_line(color = color_primary, linewidth = 1.2) +
                    labs(
                        title = paste("Effect of", gsub("_scaled", "", x_var)),
                        y = "Predicted Coverage (Prop)",
                        x = gsub("_scaled", "", x_var)
                    ) +
                    theme_publication()
            }
        }
        plot_list[[eff_name]] <- p
    }

    if (length(plot_list) > 0) {
        combined_eff <- wrap_plots(plot_list, ncol = min(3, length(plot_list))) +
            plot_annotation(
                title = "Conditional Marginal Effects (Abundance)",
                subtitle = paste("Model:", winner_id, "- Interaction effects shown separately by ARCH"),
                caption = "Blue line: Inner Arc | Orange line: Outer Arc",
                theme = theme(plot.title = element_text(face = "bold", size = 18))
            )

        h_calc <- ceiling(length(plot_list) / 3) * 4.5
        ggsave(file.path(model_fig_dir, "FIGURE_2_Marginal_Effects.png"), combined_eff, width = 14, height = max(6, h_calc), dpi = 300, bg = "white")
        cat("   ✓ Efeitos Marginais salvos\n")
    }
}

# --- C. PPC (Posterior Predictive Checks) ---
cat("   -> Gerando PPC (Zero-Validation & Density)...\n")
ppc_dens <- pp_check(model, ndraws = 100, type = "dens_overlay") +
    scale_color_manual(values = c("y" = "black", "yrep" = color_primary)) +
    labs(title = "Density Validation", subtitle = "Black: Observed | Blue: Simulated") +
    theme_publication() + theme(legend.position = "none")

ppc_zeros <- tryCatch({
    pp_check(model, ndraws = 100, type = "stat", stat = function(y) mean(y == 0)) +
        labs(title = "Zero-Inflation Validation", subtitle = "Model should capture % of zeros", x = "Proportion of Zeros", y = "Freq") +
        theme_publication() + theme(legend.position = "none")
}, error = function(e) NULL)

combined_ppc <- if (!is.null(ppc_zeros)) (ppc_dens | ppc_zeros) else ppc_dens
ggsave(file.path(model_fig_dir, "FIGURE_3_Validation_PPC.png"), combined_ppc, width = 12, height = 6, dpi = 300, bg = "white")
cat("   ✓ PPC salvo\n")

# --- D. Diagnostics (Clean) ---
cat("   -> Gerando Diagnósticos...\n")
params_to_trace <- head(all_vars[grepl("^b_|^sds_", all_vars)], 4)
if (length(params_to_trace) == 0) params_to_trace <- head(all_vars, 4)

p_trace <- mcmc_trace(model, pars = params_to_trace) +
    theme_publication() + scale_color_manual(values = okabe_ito[1:4]) +
    labs(title = "Trace Plots", subtitle = "Convergence check")

p_acf <- mcmc_acf(model, pars = params_to_trace) +
    theme_publication() + labs(title = "Autocorrelation")

combined_diag <- (p_trace / p_acf)
ggsave(file.path(model_fig_dir, "FIGURE_4_Diagnostics.png"), combined_diag, width = 10, height = 10, dpi = 300, bg = "white")
cat("   ✓ Diagnósticos salvos\n")

cat("\n✅ SUCESSO! Todas as figuras foram geradas em:\n")
cat(sprintf("   %s\n", model_fig_dir))
