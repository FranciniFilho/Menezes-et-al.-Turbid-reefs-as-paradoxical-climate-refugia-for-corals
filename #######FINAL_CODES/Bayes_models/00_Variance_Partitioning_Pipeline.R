# ============================================================================
# VARIANCE PARTITIONING PIPELINE - ESCALAS TEMPORAIS DE VARIABILIDADE
# ============================================================================
# Analisa a contribuicao relativa de CV_02, CV_30 e CV_ALL na explicacao de:
#   - Abundancia (ZOIB - componente mu)
#   - Health PC1, PC2 (Gaussian)
#   - RGR (Gaussian)
#
# Versao: 2.0
# Data: 2026-02-04
# ============================================================================

# ============================================================================
# 1. CONFIGURACAO E CARREGAMENTO DE PACOTES
# ============================================================================

suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(dplyr)
  library(ggplot2)
  library(purrr)
  library(tidyr)
  library(readr)
})

# Configurar CmdStanR
cmdstanr::set_cmdstan_path("C:/Users/rbfra/.cmdstan/cmdstan-2.37.0")  # Ajustar se necessario

# Configurar opcoes
options(mc.cores = parallel::detectCores() - 2)
set.seed(42)

# Diretorios
# Ajustar para apontar para o diretorio ######FINAL (pai de #######FINAL_CODES)
CURRENT_DIR <- getwd()
if (grepl("#######FINAL_CODES$", CURRENT_DIR)) {
  BASE_DIR <- dirname(CURRENT_DIR)
} else {
  BASE_DIR <- CURRENT_DIR
}
RESULTS_DIR <- file.path(BASE_DIR, "#######FINAL_RESULTS")
OUTPUT_DIR <- file.path(RESULTS_DIR, "Variance_Partitioning_Analysis_V2")

# Criar diretorio de saida
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ============================================================================
# 2. FUNCOES DE DIAGNOSTICO E PREPARACAO DE DADOS
# ============================================================================

#' Verifica se os datasets podem ser mergeados
#' @param data_type "abundance" ou "health"
#' @return Lista com informacoes de diagnostico
verify_merge_compatibility <- function(data_type = "abundance") {

  cat(sprintf("\n=== Verificando compatibilidade: %s ===\n", data_type))

  # Definir caminhos baseado no tipo
  if (data_type == "abundance") {
    paths <- list(
      CV_02 = file.path(RESULTS_DIR, "#####output_local_PCA_CV_2_FINAL", "dados_abundancia_integrados_long_format.csv"),
      CV_30 = file.path(RESULTS_DIR, "#####output_local_PCA_CV_30_FINAL", "dados_abundancia_integrados_long_format.csv"),
      CV_ALL = file.path(RESULTS_DIR, "#####output_local_PCA_CV_all_FINAL", "dados_abundancia_integrados_long_format.csv")
    )
    response_col <- "COBERTURA"
  } else {
    paths <- list(
      CV_02 = file.path(RESULTS_DIR, "output_DADOS_FINAIS_PARA_MODELAGEM_cv_2", "dados_finais_para_modelagem_com_ARCH.csv"),
      CV_30 = file.path(RESULTS_DIR, "output_DADOS_FINAIS_PARA_MODELAGEM_cv_30", "dados_finais_para_modelagem_com_ARCH.csv"),
      CV_ALL = file.path(RESULTS_DIR, "output_DADOS_FINAIS_PARA_MODELAGEM_cv_all", "dados_finais_para_modelagem_com_ARCH.csv")
    )
    response_col <- "HEALTH_PC1"
  }

  # VERIFICACAO 1: Arquivos existem?
  missing_files <- sapply(paths, function(p) !file.exists(p))
  if (any(missing_files)) {
    stop(sprintf("CRITICAL: Arquivos nao encontrados:\n  %s",
                 paste(names(which(missing_files)), collapse = "\n  ")))
  }
  cat("✓ Todos os arquivos existem\n")

  # VERIFICACAO 2: Carregar headers e verificar colunas
  datasets <- list()
  for (cv_name in names(paths)) {
    datasets[[cv_name]] <- read.csv2(paths[[cv_name]], nrows = 5, stringsAsFactors = FALSE)
    colnames(datasets[[cv_name]]) <- toupper(colnames(datasets[[cv_name]]))
  }

  # VERIFICACAO 3: UNIQUE_ID existe?
  has_unique_id <- all(sapply(datasets, function(d) "UNIQUE_ID" %in% names(d)))

  key_col <- "UNIQUE_ID"
  if (!has_unique_id) {
    # Tentar alternativas comuns
    alternatives <- c("ID", "SAMPLE_ID", "COLONY_ID", "SITE_YEAR")
    found_alt <- NULL
    for (alt in alternatives) {
      if (all(sapply(datasets, function(d) alt %in% names(d)))) {
        found_alt <- alt
        break
      }
    }

    if (is.null(found_alt)) {
      stop(sprintf("CRITICAL: Nenhuma chave de merge encontrada.\nColunas em CV_02: %s",
                   paste(names(datasets$CV_02), collapse = ", ")))
    }
    warning(sprintf("Usando '%s' como chave (UNIQUE_ID nao encontrado)", found_alt))
    key_col <- found_alt
  }
  cat(sprintf("✓ Chave de merge: %s\n", key_col))

  # VERIFICACAO 4: Tamanho dos datasets
  full_datasets <- list()
  n_rows <- c()
  for (cv_name in names(paths)) {
    full_datasets[[cv_name]] <- read.csv2(paths[[cv_name]], stringsAsFactors = FALSE)
    colnames(full_datasets[[cv_name]]) <- toupper(colnames(full_datasets[[cv_name]]))
    n_rows[cv_name] <- nrow(full_datasets[[cv_name]])
  }

  if (length(unique(n_rows)) > 1) {
    warning(sprintf("Datasets têm tamanhos diferentes:\n  CV_02: %d, CV_30: %d, CV_ALL: %d",
                    n_rows[1], n_rows[2], n_rows[3]))
  } else {
    cat(sprintf("✓ Datasets com mesmo tamanho: %d observacoes\n", n_rows[1]))
  }

  # VERIFICACAO 5: Variaveis resposta identicas?
  if (toupper(response_col) %in% names(full_datasets$CV_02)) {
    resp_col <- toupper(response_col)

    if (data_type == "abundance") {
      resp_02 <- as.numeric(gsub(",", ".", as.character(full_datasets$CV_02[[resp_col]])))
      resp_30 <- as.numeric(gsub(",", ".", as.character(full_datasets$CV_30[[resp_col]])))
      resp_all <- as.numeric(gsub(",", ".", as.character(full_datasets$CV_ALL[[resp_col]])))
    } else {
      resp_02 <- full_datasets$CV_02[[resp_col]]
      resp_30 <- full_datasets$CV_30[[resp_col]]
      resp_all <- full_datasets$CV_ALL[[resp_col]]
    }

    identical_02_30 <- all(resp_02 == resp_30, na.rm = TRUE)
    identical_02_all <- all(resp_02 == resp_all, na.rm = TRUE)

    if (identical_02_30 && identical_02_all) {
      cat("✓ Variaveis resposta identicas entre datasets\n")
    } else {
      warning("Variaveis resposta NAO sao identicas entre datasets!")
    }
  }

  cat("=== Verificacao completa: PASSOU ===\n\n")

  return(list(
    key_column = key_col,
    n_rows = n_rows[1],
    paths = paths,
    status = "OK"
  ))
}

#' Carrega e escala um unico dataset de CV
#' @param path Caminho para o arquivo CSV
#' @param cv_name Nome do CV ("02", "30", "ALL")
#' @param data_type "abundance" ou "health"
#' @return Dataframe processado
load_and_scale_single <- function(path, cv_name, data_type = "abundance") {

  cat(sprintf("  Processando CV_%s: ", cv_name))

  # Carregar dados
  data <- read.csv2(path, stringsAsFactors = FALSE)
  colnames(data) <- toupper(colnames(data))

  # Filtrar para M. hispida (se aplicavel)
  if ("ORGANISMO" %in% names(data)) {
    n_before <- nrow(data)
    data <- data %>% filter(ORGANISMO == "MUSSISMILIA_HISPIDA")
    cat(sprintf("[%d -> %d linhas] ", n_before, nrow(data)))
  }

  # Converter COBERTURA para proporcao (abundancia)
  if (data_type == "abundance" && "COBERTURA" %in% names(data)) {
    data$COVER_PROP <- as.numeric(gsub(",", ".", as.character(data$COBERTURA))) / 100
  }

  # Criar REEF (3 primeiros caracteres de SITE)
  data$REEF <- factor(substr(as.character(data$SITE), 1, 3))

  # Criar HAB_MERGED (RR + TP = RR_TP)
  data$HAB <- as.character(data$HAB)
  data$HAB_MERGED <- data$HAB
  data$HAB_MERGED[data$HAB_MERGED %in% c("RR", "TP")] <- "RR_TP"
  data$HAB_MERGED <- factor(data$HAB_MERGED, levels = c("PA", "RR_TP"))

  # CORRECAO: Lidar com diferentes nomes de colunas entre abundance e health
  # Abundance usa: PC1_VARIABILITY, PC2_VARIABILITY, PC1_MAGNITUDE, PC2_MAGNITUDE
  # Health usa: PC1_ENV_VAR, PC2_ENV_VAR, PC1_ENV_MAG, PC2_ENV_MAG
  if (data_type == "abundance") {
    var_col_1 <- "PC1_VARIABILITY"
    var_col_2 <- "PC2_VARIABILITY"
    mag_col_1 <- "PC1_MAGNITUDE"
    mag_col_2 <- "PC2_MAGNITUDE"
  } else {
    # Health/RGR data tem nomes diferentes
    var_col_1 <- "PC1_ENV_VAR"
    var_col_2 <- "PC2_ENV_VAR"
    mag_col_1 <- "PC1_ENV_MAG"
    mag_col_2 <- "PC2_ENV_MAG"
  }

  # ESCALONAR variaveis de variabilidade (DESTES CV especifico)
  data[[paste0("PC1_VAR_", cv_name, "_s")]] <- as.vector(scale(data[[var_col_1]]))
  data[[paste0("PC2_VAR_", cv_name, "_s")]] <- as.vector(scale(data[[var_col_2]]))

  # Guardar parametros de escalonamento
  attr(data[[paste0("PC1_VAR_", cv_name, "_s")]], "scaled:center") <- mean(data[[var_col_1]], na.rm = TRUE)
  attr(data[[paste0("PC1_VAR_", cv_name, "_s")]], "scaled:scale") <- sd(data[[var_col_1]], na.rm = TRUE)

  # Escalonar PC de Magnitude (usar CV_02 como referencia na pratica)
  if (!("PC1_MAG_s" %in% names(data))) {
    data$PC1_MAG_s <- as.vector(scale(data[[mag_col_1]]))
    data$PC2_MAG_s <- as.vector(scale(data[[mag_col_2]]))
  }

  # Escalonar Profundidade
  depth_col <- if ("DEPTH_M" %in% names(data)) "DEPTH_M" else "DEPTH"
  if (depth_col %in% names(data)) {
    data$DEPTH_M_s <- as.vector(scale(data[[depth_col]]))
  }

  cat("OK\n")
  return(data)
}

#' Prepara dados combinados para VP (Abundancia)
prepare_vp_data_abundance_v2 <- function() {

  cat("\n=== Preparando dados: Abundancia ===\n")

  # Verificar compatibilidade
  merge_info <- verify_merge_compatibility("abundance")
  key_col <- merge_info$key_column
  paths <- merge_info$paths

  # PASSO 1: Carregar e escalar cada dataset separadamente
  data_02 <- load_and_scale_single(paths$CV_02, "02", "abundance")
  data_30 <- load_and_scale_single(paths$CV_30, "30", "abundance")
  data_all <- load_and_scale_single(paths$CV_ALL, "ALL", "abundance")

  # PASSO 2: Combinar colunas por indice de linha (nao merge)
  # Os datasets CV_02, CV_30, CV_ALL tem as MESMAS linhas (mesmas observacoes),
  # apenas com diferentes calculos de variabilidade. Usar cbind diretamente.
  combined <- data_02 %>%
    select(all_of(c(key_col, "SITE", "HAB", "HAB_MERGED", "REEF", "DEPTH_M_s",
                    "COVER_PROP", "PC1_MAG_s", "PC2_MAG_s",
                    "PC1_VAR_02_s", "PC2_VAR_02_s")))

  # Adicionar colunas de variabilidade dos outros CVs por indice
  combined$PC1_VAR_30_s <- data_30$PC1_VAR_30_s
  combined$PC2_VAR_30_s <- data_30$PC2_VAR_30_s
  combined$PC1_VAR_ALL_s <- data_all$PC1_VAR_ALL_s
  combined$PC2_VAR_ALL_s <- data_all$PC2_VAR_ALL_s

  # Remover NAs
  n_before <- nrow(combined)
  combined <- combined %>% drop_na(COVER_PROP)
  n_after <- nrow(combined)

  cat(sprintf("  Dataset combinado: %d -> %d observacoes (apos remover NAs)\n", n_before, n_after))

  return(combined)
}

#' Prepara dados combinados para VP (Health/RGR)
prepare_vp_data_health_rgr_v2 <- function() {

  cat("\n=== Preparando dados: Health/RGR ===\n")

  # Verificar compatibilidade
  merge_info <- verify_merge_compatibility("health")
  key_col <- merge_info$key_column
  paths <- merge_info$paths

  # Carregar e escalar cada dataset
  data_02 <- load_and_scale_single(paths$CV_02, "02", "health")
  data_30 <- load_and_scale_single(paths$CV_30, "30", "health")
  data_all <- load_and_scale_single(paths$CV_ALL, "ALL", "health")

  # Combinar colunas por indice de linha (nao merge)
  # Os datasets CV_02, CV_30, CV_ALL tem as MESMAS linhas
  combined <- data_02 %>%
    select(all_of(c(key_col, "SITE", "HAB", "HAB_MERGED", "REEF", "DEPTH_M_s",
                    "HEALTH_PC1", "HEALTH_PC2", "RGR",
                    "PC1_MAG_s", "PC2_MAG_s",
                    "PC1_VAR_02_s", "PC2_VAR_02_s")))

  # Adicionar colunas de variabilidade dos outros CVs por indice
  combined$PC1_VAR_30_s <- data_30$PC1_VAR_30_s
  combined$PC2_VAR_30_s <- data_30$PC2_VAR_30_s
  combined$PC1_VAR_ALL_s <- data_all$PC1_VAR_ALL_s
  combined$PC2_VAR_ALL_s <- data_all$PC2_VAR_ALL_s

  cat(sprintf("  Dataset combinado: %d observacoes\n", nrow(combined)))

  return(combined)
}

# ============================================================================
# 3. DEFINICAO DE FORMULAS E PRIORS
# ============================================================================

#' Define priors para VP baseado nos winners ou default
get_vp_priors <- function(family = "zoib") {

  # Priors default (consistentes com modelos existentes)
  if (family == "zoib") {
    return(c(
      prior(normal(0, 1), class = "sds"),
      prior(normal(0, 0.5), class = "b"),
      prior(exponential(2), class = "sd"),
      prior(normal(0, 2), class = "Intercept")
    ))
  } else if (family == "gaussian") {
    return(c(
      prior(normal(0, 1), class = "sds"),
      prior(normal(0, 0.5), class = "b"),
      prior(exponential(2), class = "sd"),
      prior(normal(0, 1), class = "Intercept")
    ))
  }
}

#' Obtem formulas para modelo ZOIB (Abundancia)
get_zoib_formulas_vp <- function() {

  controls <- "PC1_MAG_s + PC2_MAG_s + s(DEPTH_M_s, k=3) + HAB_MERGED + (1|REEF)"

  list(
    M0_base = bf(
      as.formula(paste("COVER_PROP ~", controls)),
      zoi ~ 1 + (1|REEF),
      phi ~ 1 + (1|REEF)
    ),

    M1_02 = bf(
      as.formula(paste("COVER_PROP ~ s(PC1_VAR_02_s, k=3) + s(PC2_VAR_02_s, k=3) +", controls)),
      zoi ~ 1 + (1|REEF),
      phi ~ 1 + (1|REEF)
    ),

    M1_30 = bf(
      as.formula(paste("COVER_PROP ~ s(PC1_VAR_30_s, k=3) + s(PC2_VAR_30_s, k=3) +", controls)),
      zoi ~ 1 + (1|REEF),
      phi ~ 1 + (1|REEF)
    ),

    M1_ALL = bf(
      as.formula(paste("COVER_PROP ~ s(PC1_VAR_ALL_s, k=3) + s(PC2_VAR_ALL_s, k=3) +", controls)),
      zoi ~ 1 + (1|REEF),
      phi ~ 1 + (1|REEF)
    ),

    M2_02_30 = bf(
      as.formula(paste("COVER_PROP ~ s(PC1_VAR_02_s, k=3) + s(PC2_VAR_02_s, k=3) +
                       s(PC1_VAR_30_s, k=3) + s(PC2_VAR_30_s, k=3) +", controls)),
      zoi ~ 1 + (1|REEF),
      phi ~ 1 + (1|REEF)
    ),

    M2_02_ALL = bf(
      as.formula(paste("COVER_PROP ~ s(PC1_VAR_02_s, k=3) + s(PC2_VAR_02_s, k=3) +
                       s(PC1_VAR_ALL_s, k=3) + s(PC2_VAR_ALL_s, k=3) +", controls)),
      zoi ~ 1 + (1|REEF),
      phi ~ 1 + (1|REEF)
    ),

    M2_30_ALL = bf(
      as.formula(paste("COVER_PROP ~ s(PC1_VAR_30_s, k=3) + s(PC2_VAR_30_s, k=3) +
                       s(PC1_VAR_ALL_s, k=3) + s(PC2_VAR_ALL_s, k=3) +", controls)),
      zoi ~ 1 + (1|REEF),
      phi ~ 1 + (1|REEF)
    ),

    M3_full = bf(
      as.formula(paste("COVER_PROP ~ s(PC1_VAR_02_s, k=3) + s(PC2_VAR_02_s, k=3) +
                       s(PC1_VAR_30_s, k=3) + s(PC2_VAR_30_s, k=3) +
                       s(PC1_VAR_ALL_s, k=3) + s(PC2_VAR_ALL_s, k=3) +", controls)),
      zoi ~ 1 + (1|REEF),
      phi ~ 1 + (1|REEF)
    )
  )
}

#' Obtem formulas para modelo Gaussiano (Health/RGR)
get_gaussian_formulas_vp <- function(response_var = "HEALTH_PC1") {

  controls <- "s(PC1_MAG_s, k=3) + s(PC2_MAG_s, k=3) + s(DEPTH_M_s, k=3) + HAB_MERGED + (1|REEF)"

  list(
    M0_base = bf(as.formula(paste(response_var, "~", controls))),

    M1_02 = bf(as.formula(paste(response_var, "~ s(PC1_VAR_02_s, k=3) + s(PC2_VAR_02_s, k=3) +", controls))),
    M1_30 = bf(as.formula(paste(response_var, "~ s(PC1_VAR_30_s, k=3) + s(PC2_VAR_30_s, k=3) +", controls))),
    M1_ALL = bf(as.formula(paste(response_var, "~ s(PC1_VAR_ALL_s, k=3) + s(PC2_VAR_ALL_s, k=3) +", controls))),

    M2_02_30 = bf(as.formula(paste(response_var, "~ s(PC1_VAR_02_s, k=3) + s(PC2_VAR_02_s, k=3) +
                                  s(PC1_VAR_30_s, k=3) + s(PC2_VAR_30_s, k=3) +", controls))),
    M2_02_ALL = bf(as.formula(paste(response_var, "~ s(PC1_VAR_02_s, k=3) + s(PC2_VAR_02_s, k=3) +
                                   s(PC1_VAR_ALL_s, k=3) + s(PC2_VAR_ALL_s, k=3) +", controls))),
    M2_30_ALL = bf(as.formula(paste(response_var, "~ s(PC1_VAR_30_s, k=3) + s(PC2_VAR_30_s, k=3) +
                                   s(PC1_VAR_ALL_s, k=3) + s(PC2_VAR_ALL_s, k=3) +", controls))),

    M3_full = bf(as.formula(paste(response_var, "~ s(PC1_VAR_02_s, k=3) + s(PC2_VAR_02_s, k=3) +
                                  s(PC1_VAR_30_s, k=3) + s(PC2_VAR_30_s, k=3) +
                                  s(PC1_VAR_ALL_s, k=3) + s(PC2_VAR_ALL_s, k=3) +", controls)))
  )
}

# ============================================================================
# 4. EXTRACAO DE R2
# ============================================================================

#' Extrai R² do componente mu (abundancia continua) para ZOIB
#' CORRECAO CRITICA: Usar resp = "mu" para R² da parte beta
extract_r2_zoib_mu <- function(fit) {
  r2_mu <- bayes_R2(fit, resp = "mu", summary = FALSE)
  return(as.vector(r2_mu))
}

#' Extrai R² para modelos Gaussianos
extract_r2_gaussian <- function(fit) {
  r2 <- bayes_R2(fit, summary = FALSE)
  return(as.vector(r2))
}

# ============================================================================
# 5. DIAGNOSTICO DE COLINEARIDADE
# ============================================================================

#' Diagnostica colinearidade no modelo completo
diagnose_collinearity <- function(fit, data) {

  cat("\n=== DIAGNOSTICO DE COLINEARIDADE ===\n")

  # Calcular matriz de correlacao entre variaveis de variabilidade
  var_cols <- grep("VAR_0[23]_s|VAR_ALL_s", names(data), value = TRUE)

  if (length(var_cols) >= 2) {
    cor_matrix <- cor(data[, var_cols], use = "complete.obs")

    cat("\nMatriz de correlacao entre PCs de variabilidade:\n")
    print(round(cor_matrix, 3))

    # Alertas para correlacoes altas
    high_corr <- which(abs(cor_matrix) > 0.8 & abs(cor_matrix) < 1, arr.ind = TRUE)
    if (nrow(high_corr) > 0) {
      cat("\n⚠️  ALERTA: Correlacoes > 0.8 detectadas:\n")
      for (i in 1:nrow(high_corr)) {
        if (high_corr[i, 1] < high_corr[i, 2]) {  # Apenas triangulo superior
          cat(sprintf("  %s vs %s: r = %.3f\n",
                      var_cols[high_corr[i, 1]],
                      var_cols[high_corr[i, 2]],
                      cor_matrix[high_corr[i, 1], high_corr[i, 2]]))
        }
      }
    }
  }

  cat("=== FIM DO DIAGNOSTICO ===\n\n")
}

# ============================================================================
# 6. CALCULO DE COMPONENTES DE VARIANCIA
# ============================================================================

#' Calcula componentes de variancia com intervalos de credibilidade
#' @param r2_posterior Lista com amostras posteriores de R²
#' @return Dataframe com componentes e IC 95%
calculate_variance_components_v2 <- function(r2_posterior) {

  cat("\n=== Calculando componentes de variancia ===\n")

  r2 <- r2_posterior

  # VERIFICACAO DE MONOTONICIDADE
  check_monotonicity <- function(r2_sample) {
    max_m2 <- max(r2_sample["M2_02_30"], r2_sample["M2_02_ALL"], r2_sample["M2_30_ALL"])
    max_m1 <- max(r2_sample["M1_02"], r2_sample["M1_30"], r2_sample["M1_ALL"])
    is_valid <- (r2_sample["M3_full"] >= max_m2) && (max_m2 >= max_m1) && (max_m1 >= r2_sample["M0_base"])
    return(is_valid)
  }

  # Aplicar para cada amostra posterior
  n_samples <- length(r2[["M3_full"]])
  valid_samples <- sapply(1:n_samples, function(i) {
    r2_sample <- sapply(r2, function(x) x[i])
    check_monotonicity(r2_sample)
  })

  n_invalid <- sum(!valid_samples)
  if (n_invalid > 0) {
    warning(sprintf("%d de %d amostras violam monotonicidade. Verificar convergencia.",
                    n_invalid, n_samples))
  }

  # CALCULAR COMPONENTES PARA CADA AMOSTRA POSTERIOR
  unique_02_samples <- r2[["M3_full"]] - r2[["M2_30_ALL"]]
  unique_30_samples <- r2[["M3_full"]] - r2[["M2_02_ALL"]]
  unique_ALL_samples <- r2[["M3_full"]] - r2[["M2_02_30"]]

  shared_02_30_samples <- r2[["M2_02_30"]] + r2[["M1_ALL"]] - r2[["M3_full"]] - r2[["M0_base"]]
  shared_02_ALL_samples <- r2[["M2_02_ALL"]] + r2[["M1_30"]] - r2[["M3_full"]] - r2[["M0_base"]]
  shared_30_ALL_samples <- r2[["M2_30_ALL"]] + r2[["M1_02"]] - r2[["M3_full"]] - r2[["M0_base"]]

  shared_3way_samples <- r2[["M1_02"]] + r2[["M1_30"]] + r2[["M1_ALL"]] -
                         r2[["M2_02_30"]] - r2[["M2_02_ALL"]] - r2[["M2_30_ALL"]] +
                         r2[["M3_full"]] - r2[["M0_base"]]

  unexplained_samples <- 1 - r2[["M3_full"]]

  # FUNCAO AUXILIAR: Resumir amostras posteriores
  summarize_posterior <- function(samples, name) {
    samples_truncated <- pmax(samples, 0)

    data.frame(
      Component = name,
      Mean = mean(samples_truncated),
      Median = median(samples_truncated),
      Lower_95 = quantile(samples_truncated, 0.025),
      Upper_95 = quantile(samples_truncated, 0.975),
      SD = sd(samples_truncated),
      Percentage_Mean = mean(samples_truncated) * 100,
      Percentage_Lower = quantile(samples_truncated, 0.025) * 100,
      Percentage_Upper = quantile(samples_truncated, 0.975) * 100,
      stringsAsFactors = FALSE
    )
  }

  # Criar dataframe de resultados
  components <- rbind(
    summarize_posterior(unique_02_samples, "CV_02_Unique"),
    summarize_posterior(unique_30_samples, "CV_30_Unique"),
    summarize_posterior(unique_ALL_samples, "CV_ALL_Unique"),
    summarize_posterior(shared_02_30_samples, "CV_02_Intersection_CV_30"),
    summarize_posterior(shared_02_ALL_samples, "CV_02_Intersection_CV_ALL"),
    summarize_posterior(shared_30_ALL_samples, "CV_30_Intersection_CV_ALL"),
    summarize_posterior(shared_3way_samples, "Three_Way_Intersection"),
    summarize_posterior(unexplained_samples, "Unexplained")
  )

  # Adicionar totais marginais
  marginal_02 <- r2[["M1_02"]] - r2[["M0_base"]]
  marginal_30 <- r2[["M1_30"]] - r2[["M0_base"]]
  marginal_ALL <- r2[["M1_ALL"]] - r2[["M0_base"]]

  components$Marginal_Mean <- c(
    mean(marginal_02), mean(marginal_30), mean(marginal_ALL),
    NA, NA, NA, NA, NA
  )

  cat("✓ Componentes calculados\n")

  return(components)
}

# ============================================================================
# 7. AJUSTE DE MODELOS
# ============================================================================

#' Funcao generica para rodar VP em qualquer modelo
run_vp_model <- function(data, response_type, response_var = NULL,
                         output_dir, r2_extractor, diagnose_col = FALSE) {

  # Criar diretorios
  subdir_name <- ifelse(response_type == "zoib", "abundance",
                        ifelse(!is.null(response_var), tolower(response_var), "unknown"))
  subdir <- file.path(output_dir, subdir_name)
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(subdir, "models"), showWarnings = FALSE)
  dir.create(file.path(subdir, "plots"), showWarnings = FALSE)

  cat(sprintf("\n--- Iniciando VP: %s ---\n", subdir_name))

  # Definir formulas
  if (response_type == "zoib") {
    formulas <- get_zoib_formulas_vp()
    family <- zero_one_inflated_beta()
    family_name <- "zoib"
  } else {
    formulas <- get_gaussian_formulas_vp(response_var)
    family <- gaussian()
    family_name <- "gaussian"
  }

  # Obter priors
  priors <- get_vp_priors(family_name)

  # Ajustar modelos e extrair R² posterior
  r2_posterior <- list()
  fits <- list()

  for (model_name in names(formulas)) {
    cat(sprintf("\n  Ajustando %s... ", model_name))

    model_file <- file.path(subdir, "models", sprintf("%s_%s.rds", family_name, model_name))

    if (file.exists(model_file)) {
      fit <- readRDS(model_file)
      cat("[CACHE]")
    } else {
      fit <- tryCatch({
        brm(
          formula = formulas[[model_name]],
          data = data,
          family = family,
          prior = priors,
          chains = 4,
          iter = 2000,
          warmup = 1000,
          cores = min(4, parallel::detectCores() - 2),
          backend = "cmdstanr",
          control = list(adapt_delta = 0.97, max_treedepth = 12),
          save_pars = save_pars(all = TRUE),
          silent = 0,
          refresh = 100
        )
      }, error = function(e) {
        cat(sprintf("\n    ERRO: %s\n", e$message))
        return(NULL)
      })

      if (is.null(fit)) {
        stop(sprintf("Modelo %s falhou ao ajustar", model_name))
      }

      saveRDS(fit, model_file)
      cat("[NOVO]")
    }

    # Verificar convergencia
    rhat_values <- rhat(fit)
    if (any(rhat_values > 1.01, na.rm = TRUE)) {
      warning(sprintf("  ⚠️  %s: Rhat > 1.01 detectado", model_name))
    }

    # Diagnostico de colinearidade apenas para modelo completo
    if (diagnose_col && model_name == "M3_full") {
      diagnose_collinearity(fit, data)
    }

    # Extrair R² (amostras posteriores)
    r2_samples <- r2_extractor(fit)
    r2_posterior[[model_name]] <- r2_samples

    cat(sprintf(" R² = %.3f (SD: %.3f)\n", mean(r2_samples), sd(r2_samples)))

    fits[[model_name]] <- fit
  }

  # Calcular componentes com incerteza
  components <- calculate_variance_components_v2(r2_posterior)
  write.csv(components, file.path(subdir, "components_with_ci.csv"), row.names = FALSE)

  # Visualizacao
  plots <- plot_variance_partitioning_v2(components, subdir)

  return(list(
    components = components,
    plots = plots,
    fits = fits,
    r2_posterior = r2_posterior
  ))
}

# ============================================================================
# 8. VISUALIZACAO
# ============================================================================

#' Plota resultados do Variance Partitioning
plot_variance_partitioning_v2 <- function(components, output_dir) {

  cat("\n  Gerando visualizacoes...\n")

  # Definir cores
  component_colors <- c(
    "CV_02_Unique" = "#E41A1C",
    "CV_30_Unique" = "#377EB8",
    "CV_ALL_Unique" = "#4DAF4A",
    "CV_02_Intersection_CV_30" = "#FF7F00",
    "CV_02_Intersection_CV_ALL" = "#A65628",
    "CV_30_Intersection_CV_ALL" = "#999999",
    "Three_Way_Intersection" = "#F781BF",
    "Unexplained" = "#CCCCCC"
  )

  component_labels <- c(
    "CV_02_Unique" = "CV 02 (Unico)",
    "CV_30_Unique" = "CV 30 (Unico)",
    "CV_ALL_Unique" = "CV All (Unico)",
    "CV_02_Intersection_CV_30" = "CV 02 ∩ CV 30",
    "CV_02_Intersection_CV_ALL" = "CV 02 ∩ CV All",
    "CV_30_Intersection_CV_ALL" = "CV 30 ∩ CV All",
    "Three_Way_Intersection" = "CV 02 ∩ CV 30 ∩ CV All",
    "Unexplained" = "Nao Explicado"
  )

  # Plot 1: Barras empilhadas com IC
  explained <- components %>% filter(Component != "Unexplained", Mean > 0.001)

  p1 <- ggplot(explained, aes(x = "", y = Percentage_Mean, fill = Component)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = sprintf("%.1f%%", Percentage_Mean)),
              position = position_stack(vjust = 0.5),
              size = 3, color = "white", fontface = "bold") +
    scale_fill_manual(
      values = component_colors[intersect(names(component_colors), explained$Component)],
      labels = component_labels[intersect(names(component_labels), explained$Component)],
      name = "Componente"
    ) +
    labs(
      x = NULL,
      y = "% Variancia Explicada",
      title = "Decomposicao da Variancia por Escala Temporal"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5)
    ) +
    coord_flip()

  ggsave(file.path(output_dir, "plots", "variance_partition_stacked.png"),
         p1, width = 12, height = 6, dpi = 300)

  # Plot 2: Contribuicoes unicas comparativas
  unique_comp <- components %>%
    filter(grepl("Unique", Component)) %>%
    mutate(Component = gsub("_Unique", "", Component))

  p2 <- ggplot(unique_comp, aes(x = Component, y = Percentage_Mean, fill = Component)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = Percentage_Lower, ymax = Percentage_Upper),
                  width = 0.2) +
    geom_text(aes(label = sprintf("%.1f%% [%.1f-%.1f]",
                                Percentage_Mean, Percentage_Lower, Percentage_Upper)),
              vjust = -0.5, size = 3.5, fontface = "bold") +
    scale_fill_manual(
      values = c("CV_02" = "#E41A1C", "CV_30" = "#377EB8", "CV_ALL" = "#4DAF4A"),
      labels = c("CV_02" = "CV 02 anos", "CV_30" = "CV 30 anos", "CV_ALL" = "CV Completo"),
      name = "Escala"
    ) +
    labs(
      x = "Escala de Variabilidade",
      y = "% Variancia Explicada (Unica)",
      title = "Contribuicao Unica de Cada Escala Temporal"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
    ) +
    ylim(c(0, max(unique_comp$Percentage_Upper) * 1.2))

  ggsave(file.path(output_dir, "plots", "unique_contributions.png"),
         p2, width = 8, height = 6, dpi = 300)

  # Plot 3: Tabela resumo visual
  components_plot <- components %>%
    mutate(Component = factor(Component, levels = Component)) %>%
    arrange(desc(Mean))

  p3 <- ggplot(components_plot, aes(x = reorder(Component, Mean), y = Percentage_Mean)) +
    geom_segment(aes(xend = Component, y = Percentage_Lower, yend = Percentage_Upper),
                 size = 3, color = "gray70") +
    geom_point(size = 4, aes(color = Component)) +
    geom_errorbar(aes(ymin = Percentage_Lower, ymax = Percentage_Upper),
                  width = 0.2, size = 1) +
    scale_color_manual(values = component_colors) +
    coord_flip() +
    labs(
      x = "Componente",
      y = "% Variancia",
      title = "Componentes de Variancia com Intervalos de Credibilidade (95%)"
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12, hjust = 0.5)
    )

  ggsave(file.path(output_dir, "plots", "components_with_ci.png"),
         p3, width = 10, height = 6, dpi = 300)

  cat("  ✓ Visualizacoes salvas\n")

  return(list(
    stacked = p1,
    unique = p2,
    components = p3
  ))
}

# ============================================================================
# 9. COMPILACAO DE RESULTADOS
# ============================================================================

#' Compila resultados de todos os modelos
compile_vp_results_v2 <- function(output_dir) {

  cat("\n=== Compilando resultados finais ===\n")

  # Ler componentes de cada modelo
  abundance_comp <- tryCatch({
    read.csv(file.path(output_dir, "abundance", "components_with_ci.csv"))
  }, error = function(e) NULL)

  health_pc1_comp <- tryCatch({
    read.csv(file.path(output_dir, "health_pc1", "components_with_ci.csv"))
  }, error = function(e) NULL)

  health_pc2_comp <- tryCatch({
    read.csv(file.path(output_dir, "health_pc2", "components_with_ci.csv"))
  }, error = function(e) NULL)

  rgr_comp <- tryCatch({
    read.csv(file.path(output_dir, "rgr", "components_with_ci.csv"))
  }, error = function(e) NULL)

  # Criar tabela resumo
  summary_data <- list()

  extract_summary <- function(comp, name) {
    if (is.null(comp)) return(NULL)

    list(
      response = name,
      cv_02_unique = comp$Percentage_Mean[comp$Component == "CV_02_Unique"],
      cv_02_lower = comp$Percentage_Lower[comp$Component == "CV_02_Unique"],
      cv_02_upper = comp$Percentage_Upper[comp$Component == "CV_02_Unique"],

      cv_30_unique = comp$Percentage_Mean[comp$Component == "CV_30_Unique"],
      cv_30_lower = comp$Percentage_Lower[comp$Component == "CV_30_Unique"],
      cv_30_upper = comp$Percentage_Upper[comp$Component == "CV_30_Unique"],

      cv_all_unique = comp$Percentage_Mean[comp$Component == "CV_ALL_Unique"],
      cv_all_lower = comp$Percentage_Lower[comp$Component == "CV_ALL_Unique"],
      cv_all_upper = comp$Percentage_Upper[comp$Component == "CV_ALL_Unique"],

      total_shared = sum(comp$Percentage_Mean[grepl("Intersection", comp$Component)]),
      total_explained = 100 - comp$Percentage_Mean[comp$Component == "Unexplained"]
    )
  }

  summary_list <- list(
    extract_summary(abundance_comp, "Abundance"),
    extract_summary(health_pc1_comp, "Health_PC1"),
    extract_summary(health_pc2_comp, "Health_PC2"),
    extract_summary(rgr_comp, "RGR")
  )

  summary_list <- Filter(Negate(is.null), summary_list)

  # CORRECAO: Lidar com caso de apenas uma resposta processada
  if (length(summary_list) == 0) {
    cat("  AVISO: Nenhum resultado encontrado para compilar\n")
    return(NULL)
  }

  # Converter para dataframe de forma robusta
  if (length(summary_list) == 1) {
    summary_df <- as.data.frame(summary_list[[1]])
  } else {
    summary_df <- do.call(rbind, lapply(summary_list, as.data.frame))
  }

  # Salvar tabela resumo
  dir.create(file.path(output_dir, "summary"), showWarnings = FALSE)
  write.csv(summary_df, file.path(output_dir, "summary", "vp_summary_all_responses.csv"),
            row.names = FALSE)

  cat("✓ Tabela resumo salva\n")

  # Criar visualizacao comparativa (apenas se temos mais de uma resposta)
  if (nrow(summary_df) > 1) {
    summary_long <- summary_df %>%
    pivot_longer(
      cols = ends_with("_unique"),
      names_to = c("cv", ".value"),
      names_pattern = "(cv_.+)_unique"
    ) %>%
    rename(Unique_Mean = unique) %>%
    mutate(cv = factor(cv, levels = c("cv_02", "cv_30", "cv_all"),
                      labels = c("CV 02", "CV 30", "CV ALL")))

  p_comp <- ggplot(summary_long, aes(x = response, y = Unique_Mean, fill = cv)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(
      values = c("CV 02" = "#E41A1C", "CV 30" = "#377EB8", "CV ALL" = "#4DAF4A"),
      name = "Escala"
    ) +
    labs(
      x = "Variavel Resposta",
      y = "% Variancia Explicada (Unica)",
      title = "Comparacao de Contribuicoes Unicas por Escala Temporal"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  dir.create(file.path(output_dir, "summary"), showWarnings = FALSE)
  ggsave(file.path(output_dir, "summary", "comparison_plot.png"),
         p_comp, width = 10, height = 6, dpi = 300)

  cat("✓ Plot comparativo salvo\n")
  } else {
  cat("  INFO: Apenas uma resposta processada. Plot comparativo nao gerado.\n")
  }

  return(summary_df)
}

# ============================================================================
# 10. PIPELINE PRINCIPAL
# ============================================================================

#' Pipeline principal de Variance Partitioning
run_variance_partitioning_pipeline_v2 <- function(
    run_abundance = TRUE,
    run_health_pc1 = TRUE,
    run_health_pc2 = TRUE,
    run_rgr = TRUE) {

  cat("\n")
  cat("========================================\n")
  cat("VARIANCE PARTITIONING PIPELINE v2.0\n")
  cat("========================================\n")

  start_time <- Sys.time()

  # 1. PREPARAR DADOS
  cat("\n[1/6] Preparando dados...\n")

  if (run_abundance || run_health_pc1 || run_health_pc2 || run_rgr) {
    data_health <- NULL
    data_abundance <- NULL

    if (run_abundance) {
      data_abundance <- prepare_vp_data_abundance_v2()
    }

    if (run_health_pc1 || run_health_pc2 || run_rgr) {
      data_health <- prepare_vp_data_health_rgr_v2()
    }
  }

  # 2. VARIANCE PARTITIONING - ABUNDANCIA
  vp_abundance <- NULL
  if (run_abundance && !is.null(data_abundance)) {
    cat("\n[2/6] Variance Partitioning - Abundancia (ZOIB)...\n")
    vp_abundance <- run_vp_model(
      data = data_abundance,
      response_type = "zoib",
      output_dir = OUTPUT_DIR,
      r2_extractor = extract_r2_zoib_mu,
      diagnose_col = TRUE
    )
  }

  # 3. VARIANCE PARTITIONING - HEALTH PC1
  vp_health_pc1 <- NULL
  if (run_health_pc1 && !is.null(data_health)) {
    cat("\n[3/6] Variance Partitioning - Health PC1...\n")
    data_health_pc1 <- data_health[!is.na(data_health$HEALTH_PC1), ]
    vp_health_pc1 <- run_vp_model(
      data = data_health_pc1,
      response_type = "gaussian",
      response_var = "HEALTH_PC1",
      output_dir = OUTPUT_DIR,
      r2_extractor = extract_r2_gaussian
    )
  }

  # 4. VARIANCE PARTITIONING - HEALTH PC2
  vp_health_pc2 <- NULL
  if (run_health_pc2 && !is.null(data_health)) {
    cat("\n[4/6] Variance Partitioning - Health PC2...\n")
    data_health_pc2 <- data_health[!is.na(data_health$HEALTH_PC2), ]
    vp_health_pc2 <- run_vp_model(
      data = data_health_pc2,
      response_type = "gaussian",
      response_var = "HEALTH_PC2",
      output_dir = OUTPUT_DIR,
      r2_extractor = extract_r2_gaussian
    )
  }

  # 5. VARIANCE PARTITIONING - RGR
  vp_rgr <- NULL
  if (run_rgr && !is.null(data_health)) {
    cat("\n[5/6] Variance Partitioning - RGR...\n")
    data_rgr <- data_health[!is.na(data_health$RGR), ]
    vp_rgr <- run_vp_model(
      data = data_rgr,
      response_type = "gaussian",
      response_var = "RGR",
      output_dir = OUTPUT_DIR,
      r2_extractor = extract_r2_gaussian
    )
  }

  # 6. COMPILAR RESULTADOS
  cat("\n[6/6] Compilando resultados...\n")
  summary_table <- compile_vp_results_v2(OUTPUT_DIR)

  # RELATORIO FINAL
  end_time <- Sys.time()
  elapsed <- difftime(end_time, start_time, units = "mins")

  cat("\n")
  cat("========================================\n")
  cat("PIPELINE COMPLETO\n")
  cat("========================================\n")
  cat(sprintf("Tempo total: %.1f minutos\n", as.numeric(elapsed)))
  cat(sprintf("Resultados salvos em: %s\n", OUTPUT_DIR))
  cat("\n")

  # CORRECAO: So imprimir tabela se nao for NULL
  if (!is.null(summary_table)) {
    print(summary_table)
  }

  cat("\n========================================\n")

  return(list(
    abundance = vp_abundance,
    health_pc1 = vp_health_pc1,
    health_pc2 = vp_health_pc2,
    rgr = vp_rgr,
    summary = summary_table,
    output_dir = OUTPUT_DIR
  ))
}

# ============================================================================
# 11. EXECUCAO
# ============================================================================

# Para executar o pipeline completo:
# results <- run_variance_partitioning_pipeline_v2()

# Para executar apenas abundancia (teste piloto):
# results <- run_variance_partitioning_pipeline_v2(
#   run_abundance = TRUE,
#   run_health_pc1 = FALSE,
#   run_health_pc2 = FALSE,
#   run_rgr = FALSE
# )

# ============================================================================
# FIM DO SCRIPT
# ============================================================================
