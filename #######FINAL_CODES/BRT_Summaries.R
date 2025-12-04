################################################################################
##### SCRIPT A POSTERIORI PARA EXTRAÇÃO DE DETALHES DE MODELOS BRT         ####
#
# DESCRIÇÃO:
# Este script carrega um objeto de modelo BRT salvo (best_brt_object.rds)
# e extrai os principais detalhes de performance e hiperparâmetros.
# Ele gera um resumo em formato .csv e .txt para fácil consulta.
#
# INSTRUÇÕES:
# 1. Configure os 3 campos na Seção 1 para apontar para o modelo desejado.
# 2. Execute o script.
#
################################################################################

# --- PASSO 1: Especifique o diretório principal da sua análise de comparação ---
# Adicione \Outputs ao final do caminho principal
main_output_folder <- "C:/Users/rbfra/OneDrive/BRT_Health_Growth_CV_30/Outputs"

# --- PASSO 2: Especifique a pasta da variável resposta ---
# Mantenha como está
response_variable_folder <- "HEALTH_PC1"

# --- PASSO 3: Especifique o nome da pasta do modelo que você quer analisar ---
# Mantenha como está
model_name_folder <- "model_main_with_depth_with_hab"


### 2. Pacotes (Nenhum pacote extra é necessário para esta tarefa) -------------
# Este script utiliza apenas funções base do R.


### 3. LÓGICA DE EXTRAÇÃO E SALVAMENTO -----------------------------------------

cat("--- Iniciando extração de detalhes do modelo ---\n")

# Construir o caminho completo para o arquivo do modelo
model_file_path <- file.path(main_output_folder, 
                             response_variable_folder, # Para o script de abundância, isso será "Outputs"
                             model_name_folder, 
                             "BRT_Model_Outputs", 
                             "best_brt_object.rds")

cat("--- Tentando carregar o modelo de:", model_file_path, "---\n")

# Verificar se o arquivo do modelo existe
if (!file.exists(model_file_path)) {
  stop("ERRO: O arquivo 'best_brt_object.rds' não foi encontrado. Verifique os 3 campos de configuração na Seção 1.")
}

# Carregar o objeto do modelo
best_brt <- readRDS(model_file_path)
cat("--- Modelo carregado com sucesso! ---\n\n")

# --- Extraindo as informações do objeto 'best_brt' ---
cat("--- Extraindo métricas e hiperparâmetros... ---\n")

# Calcular a Deviance Explicada pela Validação Cruzada (CV Explained Deviance)
response_variable_name <- names(best_brt$data)[1] # A variável resposta é sempre a primeira coluna
y_response <- best_brt$data[[response_variable_name]]
null_deviance <- var(y_response, na.rm = TRUE)
model_deviance <- best_brt$deviance
cv_explained_deviance <- 1 - (model_deviance / null_deviance)

# Obter os preditores que foram efetivamente usados no modelo final
final_predictors <- best_brt$model$var.names

# Montar um data.frame com todos os detalhes
summary_df <- data.frame(
  Metrica = c(
    "Modelo Analisado",
    "Variável Resposta",
    "--------------------------------",
    "CV Explained Deviance",
    "CV Deviance (residual)",
    "Null Deviance (total)",
    "--------------------------------",
    "Número Ótimo de Árvores (iter)",
    "Complexidade da Árvore (tc)",
    "Taxa de Aprendizado (lr)",
    "Fração de Bag (bf)",
    "--------------------------------",
    "Número de Observações (n)",
    "Número de Preditores Usados",
    "Preditores Usados"
  ),
  Valor = c(
    model_name_folder,
    response_variable_name,
    "--------------------------------",
    paste0(round(cv_explained_deviance * 100, 3), " %"),
    round(model_deviance, 5),
    round(null_deviance, 5),
    "--------------------------------",
    best_brt$iter,
    best_brt$params$tc,
    best_brt$params$lr,
    best_brt$params$bf,
    "--------------------------------",
    nrow(best_brt$data),
    length(final_predictors),
    paste(final_predictors, collapse = ", ")
  )
)

# --- Exibir e Salvar os Resultados ---

# Criar um diretório de saída para o resumo
summary_out_dir <- file.path(main_output_folder, 
                             response_variable_folder, 
                             model_name_folder, 
                             "Best_Model_Summary_A_Posteriori")
dir.create(summary_out_dir, showWarnings = FALSE, recursive = TRUE)

# Salvar como um arquivo .csv de fácil leitura
csv_file_path <- file.path(summary_out_dir, "model_summary_details.csv")
write.csv(summary_df, csv_file_path, row.names = FALSE)
cat("--- Resumo salvo em formato CSV em:", csv_file_path, "---\n")

# Salvar como um arquivo de texto para visualização rápida
txt_file_path <- file.path(summary_out_dir, "model_summary_details.txt")
# Usar write.table com separador de tabulação para melhor alinhamento
write.table(summary_df, 
            file = txt_file_path, 
            sep = "\t", 
            row.names = FALSE, 
            quote = FALSE)
cat("--- Resumo salvo em formato TXT em:", txt_file_path, "---\n\n")


# Imprimir o resumo no console
cat("===== RESUMO DOS DETALHES DO MODELO =====\n")
print(summary_df, row.names = FALSE)
cat("=========================================\n")

cat("\n\n===== Extração de detalhes finalizada com sucesso! =====\n")