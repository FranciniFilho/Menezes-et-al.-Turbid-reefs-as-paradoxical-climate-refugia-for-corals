### ------------------------------------------------------------------------------------- ###
### 0. CONFIGURAÇÃO E PACOTES - OTIMIZADO
### ------------------------------------------------------------------------------------- ###

rm(list = ls())
gc()

# Pacotes otimizados para publicação
libs <- c("brms", "dplyr", "tidyr", "readxl", "ggplot2", "patchwork", 
          "bayesplot", "posterior", "loo", "cmdstanr", "janitor", "ggpubr")
invisible(lapply(libs, require, character.only = TRUE))

# Configurações de alta performance
options(mc.cores = parallel::detectCores(),
        brms.backend = "cmdstanr",
        mc.cores = 4)

set.seed(123)

# Diretório de output profissional
output_dir <- file.path(getwd(), "JSDM_Dirichlet_Publication_Ready")
dir.create(output_dir, showWarnings = FALSE)

### ------------------------------------------------------------------------------------- ###
### 1. CARREGAMENTO E PREPARAÇÃO DE DADOS - COMPLETO
### ------------------------------------------------------------------------------------- ###

cat("\n", rep("=", 70), "\n", sep = "")
cat("📊 CARREGAMENTO E PREPARAÇÃO DE DADOS PARA PUBLICAÇÃO\n")
cat(rep("=", 70), "\n", sep = "")

# --- Carregamento dos dados (mantido do script original) ---
path_abundance <- "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######22.04.23/DATA/BENTHOS TEMPORAL CLEAN_v3.xlsx"
path_env <- "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/########NEW RESULTS/#####output_local_PCA_CV_all_FINAL/dados_consolidados_com_scores_das_duas_PCAs.xlsx"

abundance_data <- read_excel(path_abundance, guess_max = 2000) %>%
  janitor::clean_names() %>%
  rename(SITE = site, REEF = reef, HAB = hab, YEAR = year)

env_data <- read_excel(path_env) %>%
  janitor::clean_names() %>%
  rename(SITE = site_name, REEF = reef_name, HAB = hab) %>%
  distinct(SITE, HAB, .keep_all = TRUE)

# --- Integração dos dados ---
full_data_raw <- abundance_data %>%
  select(-YEAR) %>% 
  left_join(env_data %>% select(-REEF), by = c("SITE", "HAB"))

cat("Dados integrados. Dimensões das amostras individuais:", dim(full_data_raw), "\n")

# --- Identificação e limpeza de variáveis ---
non_organism_cols <- unique(c(names(env_data), "transect", "point", "photo_code"))
organism_vars <- setdiff(names(full_data_raw), non_organism_cols)
clean_organism_vars <- janitor::make_clean_names(organism_vars)
names(full_data_raw)[names(full_data_raw) %in% organism_vars] <- clean_organism_vars
organism_vars <- clean_organism_vars

# --- Seleção, Limpeza e Escalonamento ---
combined_df <- full_data_raw %>%
  select(SITE, REEF, HAB, pc1_magnitude, pc2_magnitude, pc1_variability, pc2_variability, all_of(organism_vars))

clean_df <- combined_df %>%
  tidyr::drop_na()
cat(sprintf("\nForam removidas %d amostras individuais devido a NAs.\n", nrow(combined_df) - nrow(clean_df)))
cat("Amostras restantes para modelagem:", nrow(clean_df), "\n")

clean_df_scaled <- clean_df %>%
  mutate(across(starts_with("pc"), ~ as.numeric(scale(.))))

# --- Tratamento de zeros e redução de dimensionalidade ---
organism_data <- clean_df_scaled[, organism_vars]
min_non_zero <- min(organism_data[organism_data > 0])
correction_value <- min_non_zero / 2
organism_data_corrected <- organism_data
organism_data_corrected[organism_data_corrected == 0] <- correction_value

organism_abundance <- colMeans(organism_data_corrected)
N_TOP_ORGS <- 20
important_orgs <- names(sort(organism_abundance, decreasing = TRUE))[1:N_TOP_ORGS]
if (!"mussismilia_hispida" %in% important_orgs) {
  important_orgs <- c(important_orgs[1:(N_TOP_ORGS-1)], "mussismilia_hispida")
}
organism_vars <- important_orgs
organism_data_corrected <- organism_data_corrected[, organism_vars]

# --- Normalização ---
organism_data_normalized <- organism_data_corrected / rowSums(organism_data_corrected)

# --- CRIAÇÃO DO DATAFRAME FINAL PARA O MODELO ---
data_for_brms <- clean_df_scaled %>%
  select(SITE, REEF, HAB, pc1_magnitude, pc2_magnitude, pc1_variability, pc2_variability) %>%
  bind_cols(as.data.frame(organism_data_normalized)) %>%
  mutate(
    SITE = as.factor(SITE),
    REEF = as.factor(REEF),
    HAB = as.factor(HAB)
  )

cat("\n--- Estrutura final dos dados (amostras individuais) ---\n")
cat("Número total de amostras individuais:", nrow(data_for_brms), "\n")
cat("Número de SITES únicos:", n_distinct(data_for_brms$SITE), "\n")
cat("Número de HABITATS únicos:", n_distinct(data_for_brms$HAB), "\n")

# Verificação do aninhamento
habitat_counts_per_site <- data_for_brms %>%
  group_by(SITE) %>%
  summarise(n_unique_hab = n_distinct(HAB))

sites_com_um_hab <- habitat_counts_per_site %>%
  filter(n_unique_hab == 1)

if (nrow(sites_com_um_hab) > 0) {
  cat("INFO: Foram encontrados", nrow(sites_com_um_hab), "sítios com apenas UM habitat.\n")
}

### ------------------------------------------------------------------------------------- ###
### 2. VERIFICAÇÃO DE QUALIDADE DOS DADOS - EXPANDIDA
### ------------------------------------------------------------------------------------- ###

cat("\n🔍 VERIFICAÇÃO DE QUALIDADE PARA PUBLICAÇÃO\n")

# Verificação de completude
completeness_report <- data.frame(
  Variável = names(data_for_brms),
  Tipo = sapply(data_for_brms, class),
  N_Completo = sapply(data_for_brms, function(x) sum(!is.na(x))),
  Percentual_Completo = round(sapply(data_for_brms, function(x) mean(!is.na(x)) * 100), 1)
)

write.csv(completeness_report, file.path(output_dir, "data_completeness_report.csv"), row.names = FALSE)

# Verificação de variância zero
zero_variance_vars <- sapply(data_for_brms[, organism_vars], function(x) var(x, na.rm = TRUE) == 0)
if(any(zero_variance_vars)) {
  cat("⚠️  Variáveis com variância zero:", names(zero_variance_vars)[zero_variance_vars], "\n")
}

# Sumário estatístico dos dados
data_summary <- data_for_brms %>%
  summarise(
    n_observations = n(),
    n_sites = n_distinct(SITE),
    n_habitats = n_distinct(HAB),
    n_organisms = length(organism_vars),
    mean_abundance_per_sample = mean(rowSums(across(all_of(organism_vars)))),
    sd_abundance_per_sample = sd(rowSums(across(all_of(organism_vars))))
  )

cat("📊 RESUMO ESTATÍSTICO DOS DADOS:\n")
cat("  - Observações:", data_summary$n_observations, "\n")
cat("  - Sítios:", data_summary$n_sites, "\n")
cat("  - Habitats:", data_summary$n_habitats, "\n")
cat("  - Organismos:", data_summary$n_organisms, "\n")
cat("  - Abundância média por amostra:", round(data_summary$mean_abundance_per_sample, 4), "\n")

### ------------------------------------------------------------------------------------- ###
### 3. MODELO FINAL PARA PUBLICAÇÃO - ABORDAGEM CONSERVADORA
### ------------------------------------------------------------------------------------- ###

cat("\n", rep("=", 70), "\n", sep = "")
cat("🎯 MODELO JSDM DIRICHLET - ABORDAGEM CONSERVADORA\n")
cat(rep("=", 70), "\n", sep = "")

# PARÂMETROS FINAIS PARA PUBLICAÇÃO (MANTIDOS)
CHAINS_FINAL <- 4
ITER_FINAL <- 8000
WARMUP_FINAL <- 3000
ADAPT_DELTA_FINAL <- 0.99
MAX_TREEDEPTH_FINAL <- 15

# Fórmula final (MANTIDA)
response_vars <- paste0("cbind(", paste(organism_vars, collapse = ", "), ")")
formula_final <- paste0(response_vars, " ~ ",
                        "pc1_magnitude + pc2_magnitude + pc1_variability + pc2_variability + ",
                        "(1 | SITE / HAB)")

jsdm_formula_final <- bf(as.formula(formula_final), family = dirichlet())

cat("Fórmula Final:\n", formula_final, "\n")
cat("Parâmetros:\n")
cat(" - Cadeias:", CHAINS_FINAL, "\n")
cat(" - Iterações:", ITER_FINAL, "\n") 
cat(" - Warmup:", WARMUP_FINAL, "\n")
cat(" - adapt_delta:", ADAPT_DELTA_FINAL, "\n")
cat(" - max_treedepth:", MAX_TREEDEPTH_FINAL, "\n")

# <<< CORREÇÃO CRÍTICA: SEM PRIORS CUSTOMIZADOS >>>
cat("\n🔧 CONFIGURAÇÃO: Usando priors DEFAULT do brms para evitar problemas\n")
cat("   Motivo: Estrutura complexa do Dirichlet dificulta priors customizados\n")
cat("   Vantagem: Garante que o modelo rode sem erros de prior\n")

cat("\n⏳ Ajustando modelo final (pode levar várias horas)...\n")

model_cache_path <- file.path(output_dir, "jsdm_dirichlet_publication_model")

jsdm_brms_fit_final <- brm(
  formula = jsdm_formula_final,
  data = data_for_brms,
  # ⬅️ REMOVIDO: prior = priors_final
  chains = CHAINS_FINAL,
  iter = ITER_FINAL,
  warmup = WARMUP_FINAL,
  control = list(
    adapt_delta = ADAPT_DELTA_FINAL,
    max_treedepth = MAX_TREEDEPTH_FINAL
  ),
  backend = "cmdstanr",
  file = model_cache_path,
  silent = 2,
  refresh = 500,
  save_pars = save_pars(all = TRUE)
)

### ------------------------------------------------------------------------------------- ###
### 4. SALVAMENTO COMPLETO PARA ANÁLISES
### ------------------------------------------------------------------------------------- ###

cat("\n💾 SALVANDO DADOS COMPLETOS PARA ANÁLISES...\n")

# Salvar modelo principal
saveRDS(jsdm_brms_fit_final, file.path(output_dir, "modelo_principal_publicacao.rds"))

# Salvar dados completos para Script 2
saveRDS(data_for_brms, file.path(output_dir, "data_for_brms_completo.rds"))
saveRDS(organism_vars, file.path(output_dir, "organism_vars_completo.rds"))

# Metadados expandidos
metadata_publicacao <- list(
  data_preparation_date = Sys.time(),
  n_observations = nrow(data_for_brms),
  n_sites = n_distinct(data_for_brms$SITE),
  n_habitats = n_distinct(data_for_brms$HAB),
  n_organisms = length(organism_vars),
  model_formula = formula_final,
  model_parameters = list(
    chains = CHAINS_FINAL,
    iter = ITER_FINAL,
    warmup = WARMUP_FINAL,
    adapt_delta = ADAPT_DELTA_FINAL,
    max_treedepth = MAX_TREEDEPTH_FINAL
  ),
  organisms_included = organism_vars,
  data_summary = data_summary
)

saveRDS(metadata_publicacao, file.path(output_dir, "metadata_publicacao.rds"))

# Relatório rápido de diagnóstico
sink(file.path(output_dir, "modelo_quick_diagnostic.txt"))
cat("=== DIAGNÓSTICO RÁPIDO DO MODELO ===\n\n")
print(summary(jsdm_brms_fit_final)$fixed)
cat("\n--- CONVERGÊNCIA ---\n")
if(!is.null(summary(jsdm_brms_fit_final)$fixed)) {
  cat("Rhat máximo:", max(summary(jsdm_brms_fit_final)$fixed$Rhat), "\n")
}
sink()

cat("\n✅ SCRIPT 1 CONCLUÍDO! Modelo pronto para análises detalhadas.\n")
cat("📁 Arquivos salvos em:", output_dir, "\n")
cat("📊 Resumo dos dados:\n")
cat("   -", data_summary$n_observations, "observações\n")
cat("   -", data_summary$n_sites, "sítios\n") 
cat("   -", data_summary$n_habitats, "habitats\n")
cat("   -", data_summary$n_organisms, "organismos\n")
cat("🚀 Execute o Script 2 para análises completas\n")