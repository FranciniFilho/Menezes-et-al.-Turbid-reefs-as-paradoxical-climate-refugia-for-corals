### INTEGRAÇÃO FINAL DOS DADOS (NÍVEL COLÔNIA) - VERSÃO 7 (COM LAT/LON) ###
# 1. Adiciona as colunas 'LATITUDE' e 'LONGITUDE' ao dataframe final.
# 2. Mantém a adição da coluna 'ARCH'.
# 3. Garante que as novas colunas sejam incluídas na ordem correta no arquivo de saída.

import os
import pandas as pd
import glob

# ---------------------------
# 1. CONFIGURAÇÕES E CAMINHOS
# ---------------------------
print("--- Iniciando Script de Integração Final (Nível Colônia) - Versão 7 (com Lat/Lon) ---")

# --- Caminhos de Entrada ---
sites_csv_path = r"#######FINAL_DATA/00_sites_metadata/sites_list_full.csv"
path_pca_local = r"#######FINAL_DATA/02_PCA_environmental/CV_ALL"
file_ambiental_scores = os.path.join(path_pca_local, "dados_consolidados_com_scores_das_duas_PCAs.xlsx")
path_biologico = r"#######FINAL_DATA/04_health_growth_PCA"
file_resultados_biologicos = os.path.join(path_biologico, "resultados_biologicos_por_colonia.csv")
path_frequencia = r"#######FINAL_DATA/07_lomb_frequencies"

# --- Diretório de Saída ---
output_dir = r"#######FINAL_RESULTS/output_DADOS_FINAIS_PARA_MODELAGEM_cv_all_lat"
os.makedirs(output_dir, exist_ok=True)
print(f"Diretório de saída definido para: {output_dir}")

# ---------------------------
# 2. LEITURA E PADRONIZAÇÃO IMEDIATA
# ---------------------------
print("\nLendo e padronizando arquivos de entrada...")
try:
    df_sites_meta = pd.read_csv(sites_csv_path, sep=';')
    df_sites_meta.columns = [col.strip().upper() for col in df_sites_meta.columns]
    print(">>> Arquivo de metadados dos sítios lido e padronizado.")

    df_ambiental = pd.read_excel(file_ambiental_scores)
    # Substitui vírgulas por pontos nas colunas de coordenadas e converte para numérico
    df_ambiental['Latitude'] = df_ambiental['Latitude'].astype(str).str.replace(',', '.').astype(float)
    df_ambiental['Longitude'] = df_ambiental['Longitude'].astype(str).str.replace(',', '.').astype(float)
    df_ambiental.columns = [col.strip().upper() for col in df_ambiental.columns]
    print(">>> Arquivo de dados ambientais lido e padronizado (Lat/Lon convertidas para numérico).")

    df_biologico = pd.read_csv(file_resultados_biologicos)
    df_biologico.columns = [col.strip().upper() for col in df_biologico.columns]
    print(">>> Arquivo de dados biológicos lido e padronizado.")

    frequency_files = glob.glob(os.path.join(path_frequencia, "frequencias_*.csv"))
    if frequency_files:
        df_frequencia = pd.read_csv(frequency_files[0], sep=';', decimal=',')
        df_frequencia.columns = [col.strip().upper() for col in df_frequencia.columns]
        print(f">>> Arquivo de frequência '{os.path.basename(frequency_files[0])}' lido e padronizado.")
    else:
        df_frequencia = None
        print("AVISO: Nenhum arquivo de frequência encontrado.")

except Exception as e:
    print(f"\nERRO CRÍTICO DURANTE A LEITURA DE UM DOS ARQUIVOS DE ENTRADA: {e}")
    exit()

# ---------------------------
# 3. CRIAÇÃO DE CHAVES E INTEGRAÇÃO
# ---------------------------
print("\nCriando chaves de junção e integrando os dados...")

df_ambiental.rename(columns={'SITE_NAME': 'SITE'}, inplace=True, errors='ignore')
df_biologico.rename(columns={'SITE_NAME': 'SITE'}, inplace=True, errors='ignore')
if df_frequencia is not None:
    df_frequencia.rename(columns={'SITE_NAME': 'SITE'}, inplace=True, errors='ignore')
df_sites_meta.rename(columns={'SITE_NAME': 'SITE'}, inplace=True, errors='ignore')

if df_frequencia is not None and 'HAB' not in df_frequencia.columns:
    meta_subset = df_sites_meta[['SITE', 'HAB']].drop_duplicates()
    df_frequencia = pd.merge(df_frequencia, meta_subset, on='SITE', how='left')
    print("Coluna 'HAB' adicionada ao dataframe de frequência.")

for df, name in [(df_ambiental, 'ambiental'), (df_biologico, 'biologico'), (df_frequencia, 'frequencia')]:
    if df is not None and 'SITE' in df.columns and 'HAB' in df.columns:
        df['UNIQUE_ID'] = df['SITE'].astype(str).str.strip() + '_' + df['HAB'].astype(str).str.strip()
        print(f"Chave 'UNIQUE_ID' criada/verificada para o dataframe '{name}'.")

# --- SELEÇÃO DE COLUNAS E MERGE ---
### ALTERAÇÃO 1: Adicionar 'LATITUDE' e 'LONGITUDE' à lista de colunas a manter.
cols_ambientais_a_manter = [
    'UNIQUE_ID', 'ARCH', 'LATITUDE', 'LONGITUDE', 
    'PC1_MAGNITUDE', 'PC2_MAGNITUDE', 'PC1_VARIABILITY', 'PC2_VARIABILITY',
    'SST_MEAN', 'SST_CV_ALL', 'MEAN_DLI_LOCAL', 'CV_DLI_LOCAL',
    'CHL_MEAN', 'CHL_CV_ALL', 'PROP_DHW_GT4', 'PROP_DHW_GT8', 'DEPTH_M'
]
existing_cols_ambientais = [col for col in cols_ambientais_a_manter if col in df_ambiental.columns]
df_ambiental_subset = df_ambiental[existing_cols_ambientais].drop_duplicates(subset=['UNIQUE_ID'])
df_ambiental_subset = df_ambiental_subset.rename(columns={
    'PC1_MAGNITUDE': 'PC1_ENV_MAG', 'PC2_MAGNITUDE': 'PC2_ENV_MAG',
    'PC1_VARIABILITY': 'PC1_ENV_VAR', 'PC2_VARIABILITY': 'PC2_ENV_VAR',
    'DEPTH_M': 'DEPTH'
})
print(f"\nColunas ambientais selecionadas para merge: {df_ambiental_subset.columns.tolist()}")

# --- Merge 1: Biológico + Ambiental ---
dados_finais = pd.merge(df_biologico, df_ambiental_subset, on='UNIQUE_ID', how='left')
print(f"Shape após merge com dados ambientais: {dados_finais.shape}")

# --- Merge 2: Adicionar dados de Frequência ---
if df_frequencia is not None:
    cols_frequencia_desejadas = [
        'UNIQUE_ID', 'SST_PERIOD_PRIMARY', 'DLI_PERIOD_PRIMARY', 'CHL_PERIOD_PRIMARY',
        'SST_PERIOD_SECONDARY', 'DLI_PERIOD_SECONDARY', 'CHL_PERIOD_SECONDARY'
    ]
    existing_cols_freq = [col for col in cols_frequencia_desejadas if col in df_frequencia.columns]
    df_frequencia_subset = df_frequencia[existing_cols_freq].drop_duplicates(subset=['UNIQUE_ID'])
    dados_finais = pd.merge(dados_finais, df_frequencia_subset, on='UNIQUE_ID', how='left')
    print(f"Shape após merge com dados de frequência: {dados_finais.shape}")

# ---------------------------
# 4. FORMATAÇÃO E SALVAMENTO FINAL
# ---------------------------
if 'MEAN_DEAD_TISSUE' not in dados_finais.columns and 'MEAN_HEALTH' in dados_finais.columns and 'MEAN_BLEACH' in dados_finais.columns:
    dados_finais['MEAN_DEAD_TISSUE'] = 100 - dados_finais['MEAN_HEALTH'] - dados_finais['MEAN_BLEACH']

### ALTERAÇÃO 2: Adicionar 'LATITUDE' e 'LONGITUDE' à ordem final das colunas.
final_columns_order = [
    # Identificadores da Colônia e Sítio
    'SITE_COL', 'SITE', 'HAB', 'REEF', 'ARCH', 'LATITUDE', 'LONGITUDE', 'UNIQUE_ID',
    
    # Variáveis de Resposta Biológica
    'RGR', 'AR_TOTAL_INICIAL', 'MEAN_AR_TOTAL', 'MEAN_HEALTH', 'MEAN_BLEACH', 'MEAN_DEAD_TISSUE',
    
    # Scores das PCAs Biológicas
    'HEALTH_PC1', 'HEALTH_PC2', 'PC1_INTERACAO', 'PC2_INTERACAO',
    
    # Scores das PCAs Ambientais
    'PC1_ENV_MAG', 'PC2_ENV_MAG', 'PC1_ENV_VAR', 'PC2_ENV_VAR',
    
    # Variáveis Ambientais Brutas e de Frequência
    'SST_MEAN', 'SST_CV_ALL', 'MEAN_DLI_LOCAL', 'CV_DLI_LOCAL',
    'CHL_MEAN', 'CHL_CV_ALL', 'PROP_DHW_GT4', 'PROP_DHW_GT8', 'DEPTH',
    'SST_PERIOD_PRIMARY', 'DLI_PERIOD_PRIMARY', 'CHL_PERIOD_PRIMARY',
    'SST_PERIOD_SECONDARY', 'DLI_PERIOD_SECONDARY', 'CHL_PERIOD_SECONDARY'
]
final_columns_exist = [col for col in final_columns_order if col in dados_finais.columns]
dados_finais_ordenados = dados_finais[final_columns_exist]

print("\nVerificação de valores ausentes (NaN) na planilha final:")
print(dados_finais_ordenados.isnull().sum())

# Define um nome de arquivo para esta nova versão
output_file = os.path.join(output_dir, "dados_finais_para_modelagem_com_ARCH_e_Coords.csv")
dados_finais_ordenados.to_csv(output_file, index=False, sep=';', decimal=',')

print(f"\n>>> SUCESSO! Planilha final (com ARCH, Latitude e Longitude) integrada salva em: {output_file}")