### ==================================================================================== ###
### === AUXILIARY: GERAR SUMÁRIO E VALIDAÇÃO LOO PARA MODELO VENCEDOR DE ABUNDÂNCIA  === ###
### === Gera model_summary.txt e VALIDACAO_LOO_RIGOROSA.txt a partir do modelo salvo === ###
### ==================================================================================== ###

rm(list = ls())
gc()

### 1. PACOTES NECESSÁRIOS -----------------------------------------------------------
libs <- c("brms", "loo", "cmdstanr")
missing_packages <- libs[!libs %in% installed.packages()[, "Package"]]
if (length(missing_packages) > 0) {
    cat("Instalando pacotes faltantes:", paste(missing_packages, collapse = ", "), "\n")
    install.packages(missing_packages, dependencies = TRUE)
}
invisible(lapply(libs, library, character.only = TRUE))

### 2. CAMINHOS ----------------------------------------------------------------------
# O modelo vencedor já existe como arquivo .rds na raiz do diretório de análise
model_path <- "C:/Users/rbfra/OneDrive/Bayesian_Analyses_ZOIB_Abundance_MULTI_CV/ZOIB_PCA_Informative_model_interaction_arch_no_depth_with_hab.rds"
output_dir <- "C:/Users/rbfra/OneDrive/Bayesian_Analyses_ZOIB_Abundance_MULTI_CV/GLOBAL_WINNER_OUTPUTS"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n", rep("=", 75), "\n", sep = "")
cat("📊 GERADOR DE SUMÁRIO E VALIDAÇÃO - MODELO DE ABUNDÂNCIA\n")
cat(rep("=", 75), "\n", sep = "")

### 3. VERIFICAR EXISTÊNCIA DO MODELO -----------------------------------------------
if (!file.exists(model_path)) {
    stop(sprintf("ERRO: Modelo não encontrado em:\n  %s", model_path))
}

### 4. CARREGAR MODELO --------------------------------------------------------------
cat("\n📦 Carregando modelo vencedor (pode demorar alguns segundos)...\n")
model <- readRDS(model_path)
winner_id <- "PCA_Informative_model_interaction_arch_no_depth_with_hab"
file_prefix <- paste0("GLOBAL_WINNER_ABUNDANCE_", winner_id)

cat(sprintf("   ✓ Modelo carregado: %s\n", winner_id))

### 5. SALVAR CÓPIA COM PREFIXO PARA COMPATIBILIDADE COM MASTER_VIZ_PIPELINE -------
cat("\n💾 Salvando cópia com prefixo WINNER para compatibilidade...\n")
winner_rds_path <- file.path(output_dir, paste0("GLOBAL_WINNER_ABUNDANCE_", winner_id, ".rds"))
saveRDS(model, file = winner_rds_path)
cat(sprintf("   ✓ Modelo salvo em: %s\n", winner_rds_path))

### 6. GERAR SUMÁRIO DO MODELO ------------------------------------------------------
cat("\n📊 Gerando sumário do modelo...\n")
summary_file <- file.path(output_dir, sprintf("model_summary_%s.txt", file_prefix))

# Sumário básico
capture.output(summary(model), file = summary_file)

# Adicionar R² Bayesiano
cat("\n\n--- R² Bayesiano ---\n", file = summary_file, append = TRUE)
r2 <- bayes_R2(model)
capture.output(print(r2), file = summary_file, append = TRUE)

# Adicionar informações sobre convergência
cat("\n\n--- Diagnósticos de Convergência ---\n", file = summary_file, append = TRUE)
rhat_summary <- summary(model)$fixed
cat(sprintf("   Número de parâmetros fixos: %d\n", nrow(rhat_summary)), file = summary_file, append = TRUE)
cat(sprintf("   Range de Rhat: [%.3f, %.3f]\n", min(rhat_summary$Rhat), max(rhat_summary$Rhat)), file = summary_file, append = TRUE)

cat(sprintf("   ✓ Sumário salvo em:\n     %s\n", summary_file))

### 7. GERAR VALIDAÇÃO LOO RIGOROSA -------------------------------------------------
cat("\n🔬 Executando LOO Rigoroso (moment_match=TRUE)...\n")
cat("   ⚠ AVISO: Isso pode demorar 1-5 minutos dependendo do hardware.\n")

loo_rigoroso <- tryCatch(
  {
    loo(model, moment_match = TRUE, cores = 2)
  },
  error = function(e) {
    cat(sprintf("   ⚠ AVISO: moment_match falhou com erro: %s\n", e$message))
    cat("   → Calculando LOO padrão (sem moment matching)...\n")
    return(loo(model, moment_match = FALSE, cores = 2))
  }
)

# Salvar validação LOO
loo_file <- file.path(output_dir, sprintf("VALIDACAO_LOO_RIGOROSA_%s.txt", file_prefix))
capture.output(print(loo_rigoroso), file = loo_file)

# Adicionar diagnósticos sobre Pareto-k
k_values <- loo_rigoroso$diagnostics$pareto_k
n_bad_k <- sum(k_values > 0.7)
pct_bad_k <- 100 * n_bad_k / length(k_values)

cat("\n\n--- Diagnósticos Pareto-k ---\n", file = loo_file, append = TRUE)
cat(sprintf("   Observações com k > 0.7: %d / %d (%.1f%%)\n", 
            n_bad_k, length(k_values), pct_bad_k), file = loo_file, append = TRUE)
cat(sprintf("   Range de k: [%.3f, %.3f]\n", min(k_values), max(k_values)), file = loo_file, append = TRUE)

cat(sprintf("   ✓ Validação LOO salva em:\n     %s\n", loo_file))

### 8. RELATÓRIO FINAL --------------------------------------------------------------
cat("\n", rep("=", 75), "\n", sep = "")
cat("✅ PROCESSO CONCLUÍDO COM SUCESSO!\n")
cat(rep("=", 75), "\n", sep = "")
cat(sprintf("\n📁 Diretório de saída: %s\n", output_dir))
cat("\nArquivos gerados:\n")
cat(sprintf("  • %s\n", basename(winner_rds_path)))
cat(sprintf("  • %s\n", basename(summary_file)))
cat(sprintf("  • %s\n", basename(loo_file)))
cat("\n")
