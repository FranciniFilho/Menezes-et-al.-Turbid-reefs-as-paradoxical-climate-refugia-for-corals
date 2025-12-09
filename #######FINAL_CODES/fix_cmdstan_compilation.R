# Script de Diagnóstico e Correção Automática do CmdStan (Versão 2 - Robustez Melhorada)
# Execute este script para corrigir os erros de compilação (OpenCL duplicado e C++20 warnings)

library(cmdstanr)

cat("=== DIAGNÓSTICO DO CMDSTAN (V2) ===\n")

# 1. Verificar caminho e versão
cmdstan_path_dir <- tryCatch(cmdstan_path(), error = function(e) NULL)

if (is.null(cmdstan_path_dir)) {
    stop("CmdStan não encontrado. Verifique se o cmdstanr está configurado corretamente.")
}

cat(sprintf("Caminho do CmdStan: %s\n", cmdstan_path_dir))
cat(sprintf("Versão do CmdStan: %s\n", cmdstan_version()))

# 2. Verificar arquivo make/local
make_local_path <- file.path(cmdstan_path_dir, "make", "local")
cat(sprintf("Arquivo de configuração (make/local): %s\n", make_local_path))

current_content <- if (file.exists(make_local_path)) readLines(make_local_path) else character(0)

cat("\n--- Conteúdo Atual do make/local ---\n")
if (length(current_content) > 0) {
    cat(paste(current_content, collapse = "\n"))
} else {
    cat("(Arquivo vazio ou inexistente)")
}
cat("\n------------------------------------\n")

# 3. Aplicar Correções
cat("\n=== APLICANDO CORREÇÕES ===\n")
new_content <- current_content

# Correção 1: Remover STAN_OPENCL=true (ou TRUE) se existir
# O brms já adiciona essa flag quando opencl() é usado.
# A regex agora é case-insensitive e mais robusta
lines_to_remove <- grep("^STAN_OPENCL\\s*=\\s*(true|TRUE)", new_content, ignore.case = TRUE)

if (length(lines_to_remove) > 0) {
    cat("Detectado 'STAN_OPENCL=TRUE' global. Removendo para evitar conflito com brms...\n")
    new_content <- new_content[-lines_to_remove]
}

# Correção 2: Forçar C++14 para evitar warnings/erros de C++20 com TBB
# Verifica se já existe alguma flag de CXXFLAGS com std=c++14
if (!any(grepl("CXXFLAGS.*-std=c\\+\\+14", new_content))) {
    cat("Adicionando flag '-std=c++14' para compatibilidade com TBB...\n")
    new_content <- c(new_content, "CXXFLAGS += -std=c++14")
}

# Salvar alterações
if (!identical(current_content, new_content)) {
    writeLines(new_content, make_local_path)
    cat("\nArquivo make/local atualizado com sucesso.\n")

    cat("\n--- Novo Conteúdo do make/local ---\n")
    cat(paste(new_content, collapse = "\n"))
    cat("\n-----------------------------------\n")

    # 4. Recompilar CmdStan
    cat("\n=== RECOMPILANDO CMDSTAN ===\n")
    cat("Isso pode levar alguns minutos...\n")

    # Tenta recompilar. Se falhar, avisa o usuário.
    tryCatch(
        {
            cmdstan_make()
            cat("\nCmdStan reconstruído com sucesso!\n")
        },
        error = function(e) {
            cat("\nERRO ao reconstruir CmdStan: ", e$message, "\n")
            cat("Tente rodar 'cmdstan_make()' manualmente no console.\n")
        }
    )
} else {
    cat("\nNenhuma alteração necessária no make/local.\n")
    cat("Se o erro persistir, tente reconstruir o CmdStan manualmente com: cmdstan_make()\n")
}

cat("\n=== CONCLUÍDO ===\n")
cat("Tente rodar seu script de modelagem novamente.\n")
