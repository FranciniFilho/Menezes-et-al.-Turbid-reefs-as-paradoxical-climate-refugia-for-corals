# Roda PERMANOVA, ANOVA e Levene para REEF e HAB
# PCA LOCAL

# =============================================================================
# ANÁLISE ESTATÍSTICA (PERMANOVA/ANOVA) DOS SCORES DA PCA LOCAL
# Versão 2.0 - Adaptado para ler a saída da PCA Local e analisar
#              separadamente a PCA de Magnitude e a de Variabilidade.
# =============================================================================

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from skbio.stats.distance import DistanceMatrix
from skbio.stats.distance import permanova
from scipy.spatial.distance import pdist, squareform
from scipy.stats import f_oneway, levene

# === CAMINHOS E CONFIGURAÇÕES ===
# Aponta para o arquivo de saída padrão do script PCA LOCAL
input_file = r"#######FINAL_DATA/02_PCA_environmental/CV_02/dados_consolidados_com_scores_das_duas_PCAs.xlsx"

# Diretório principal para os resultados desta análise
base_output_dir = r"#######FINAL_RESULTS/PERMANOVA_ANOVA_LOCAL_PCA_CV_2"
os.makedirs(base_output_dir, exist_ok=True)


# === FUNÇÃO PARA RODAR A ANÁLISE COMPLETA ===
def run_full_analysis(df, pc1_col, pc2_col, analysis_name, output_dir):
    """
    Executa PERMANOVA, ANOVA, Levene e gera gráficos para um conjunto de scores de PCA.
    
    Args:
        df (pd.DataFrame): DataFrame com os dados.
        pc1_col (str): Nome da coluna do PC1.
        pc2_col (str): Nome da coluna do PC2.
        analysis_name (str): Nome da análise (ex: "Magnitude", "Variability").
        output_dir (str): Diretório para salvar os resultados.
    """
    os.makedirs(output_dir, exist_ok=True)
    print("\n" + "="*60)
    print(f"  INICIANDO ANÁLISE PARA: PCA {analysis_name.upper()}")
    print("="*60)

    # --- 1. Preparação dos Dados e Matriz de Distância ---
    # Colunas essenciais para esta análise específica
    required_cols = [pc1_col, pc2_col, 'Reef_name', 'Site_name', 'HAB']
    df_valid = df.dropna(subset=required_cols).copy()
    
    if len(df_valid) < 3:
        print(f"AVISO: Dados insuficientes ({len(df_valid)} linhas) para a análise {analysis_name}. Pulando.")
        return

    print(f"Dados válidos para {analysis_name}: {len(df_valid)} linhas.")
    
    data_matrix = df_valid[[pc1_col, pc2_col]].values
    dist_matrix = squareform(pdist(data_matrix, metric='euclidean'))
    # Garante que os IDs sejam únicos para a matriz de distância
    unique_ids = df_valid.reset_index().apply(lambda row: f"{row['Site_name']}_{row['HAB']}_{row['index']}", axis=1).tolist()
    dm = DistanceMatrix(dist_matrix, ids=unique_ids)

    # --- 2. Análise por REEF ---
    print(f"\n--- Analisando por REEF ({analysis_name}) ---")
    grouping_reef = df_valid.reset_index().set_index(pd.Index(unique_ids))['Reef_name']
    permanova_reef = permanova(distance_matrix=dm, grouping=grouping_reef, permutations=999)

    grouped_pc1_reef = [g[pc1_col].values for _, g in df_valid.groupby('Reef_name')]
    anova_pc1_reef = f_oneway(*grouped_pc1_reef)
    levene_pc1_reef = levene(*grouped_pc1_reef)

    grouped_pc2_reef = [g[pc2_col].values for _, g in df_valid.groupby('Reef_name')]
    anova_pc2_reef = f_oneway(*grouped_pc2_reef)
    levene_pc2_reef = levene(*grouped_pc2_reef)

    fig_reef, axes = plt.subplots(1, 2, figsize=(14, 6))
    sns.boxplot(ax=axes[0], x='Reef_name', y=pc1_col, data=df_valid, palette='viridis')
    axes[0].set_title(f'{pc1_col} por REEF\nANOVA p={anova_pc1_reef.pvalue:.4f} | Levene p={levene_pc1_reef.pvalue:.4f}')
    axes[0].tick_params(axis='x', rotation=45)
    sns.boxplot(ax=axes[1], x='Reef_name', y=pc2_col, data=df_valid, palette='viridis')
    axes[1].set_title(f'{pc2_col} por REEF\nANOVA p={anova_pc2_reef.pvalue:.4f} | Levene p={levene_pc2_reef.pvalue:.4f}')
    axes[1].tick_params(axis='x', rotation=45)
    fig_reef.suptitle(f"Análise por Recife - PCA de {analysis_name}", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(os.path.join(output_dir, f"boxplots_PCA_{analysis_name}_vs_REEF.png"), dpi=300)
    plt.close()

    # --- 3. Análise por HABITAT ---
    print(f"\n--- Analisando por HAB ({analysis_name}) ---")
    grouping_hab = df_valid.reset_index().set_index(pd.Index(unique_ids))['HAB']
    permanova_hab = permanova(distance_matrix=dm, grouping=grouping_hab, permutations=999)
    
    grouped_pc1_hab = [g[pc1_col].values for _, g in df_valid.groupby('HAB')]
    anova_pc1_hab = f_oneway(*grouped_pc1_hab)
    levene_pc1_hab = levene(*grouped_pc1_hab)

    grouped_pc2_hab = [g[pc2_col].values for _, g in df_valid.groupby('HAB')]
    anova_pc2_hab = f_oneway(*grouped_pc2_hab)
    levene_pc2_hab = levene(*grouped_pc2_hab)
    
    fig_hab, axes = plt.subplots(1, 2, figsize=(12, 5))
    sns.boxplot(ax=axes[0], x='HAB', y=pc1_col, data=df_valid, palette='pastel')
    axes[0].set_title(f'{pc1_col} por HABITAT\nANOVA p={anova_pc1_hab.pvalue:.4f} | Levene p={levene_pc1_hab.pvalue:.4f}')
    sns.boxplot(ax=axes[1], x='HAB', y=pc2_col, data=df_valid, palette='pastel')
    axes[1].set_title(f'{pc2_col} por HABITAT\nANOVA p={anova_pc2_hab.pvalue:.4f} | Levene p={levene_pc2_hab.pvalue:.4f}')
    fig_hab.suptitle(f"Análise por Habitat - PCA de {analysis_name}", fontsize=16)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(os.path.join(output_dir, f"boxplots_PCA_{analysis_name}_vs_HAB.png"), dpi=300)
    plt.close()

    # --- 4. Salvar Resultados Numéricos ---
    output_txt = os.path.join(output_dir, f"resultados_estatisticos_{analysis_name}.txt")
    with open(output_txt, 'w') as f:
        f.write("="*50 + "\n")
        f.write(f"RESULTADOS ESTATÍSTICOS - PCA de {analysis_name.upper()}\n")
        f.write("="*50 + "\n\n")

        f.write("PERMANOVA - ESTRUTURAÇÃO MULTIVARIADA\n")
        f.write("-" * 50 + "\n")
        f.write("PERMANOVA baseada em REEF:\n")
        f.write(str(permanova_reef) + "\n\n")
        f.write("PERMANOVA baseada em HAB:\n")
        f.write(str(permanova_hab) + "\n\n")

        f.write("ANOVA/LEVENE - ESTRUTURAÇÃO UNIVARIADA POR EIXO\n")
        f.write("-" * 50 + "\n")
        f.write("Resultados para REEF:\n")
        f.write(f"  {pc1_col} vs REEF: ANOVA F={anova_pc1_reef.statistic:.3f}, p={anova_pc1_reef.pvalue:.4f} | Levene W={levene_pc1_reef.statistic:.3f}, p={levene_pc1_reef.pvalue:.4f}\n")
        f.write(f"  {pc2_col} vs REEF: ANOVA F={anova_pc2_reef.statistic:.3f}, p={anova_pc2_reef.pvalue:.4f} | Levene W={levene_pc2_reef.statistic:.3f}, p={levene_pc2_reef.pvalue:.4f}\n\n")
        f.write("Resultados para HAB:\n")
        f.write(f"  {pc1_col} vs HAB: ANOVA F={anova_pc1_hab.statistic:.3f}, p={anova_pc1_hab.pvalue:.4f} | Levene W={levene_pc1_hab.statistic:.3f}, p={levene_pc1_hab.pvalue:.4f}\n")
        f.write(f"  {pc2_col} vs HAB: ANOVA F={anova_pc2_hab.statistic:.3f}, p={anova_pc2_hab.pvalue:.4f} | Levene W={levene_pc2_hab.statistic:.3f}, p={levene_pc2_hab.pvalue:.4f}\n\n")
    
    print(f"Análise para PCA de {analysis_name} concluída. Resultados em: {output_dir}")


# === SCRIPT PRINCIPAL ===

# --- Carrega e prepara os dados UMA VEZ ---
try:
    df_full = pd.read_excel(input_file)
    print(f"Arquivo '{os.path.basename(input_file)}' carregado com sucesso.")
except FileNotFoundError:
    print(f"ERRO: O arquivo de entrada não foi encontrado em: {input_file}")
    exit()

# --- ANÁLISE 1: PCA DE MAGNITUDE ---
# Define as colunas de score para a PCA de Magnitude
pc1_mag_col = 'PC1_Magnitude'
pc2_mag_col = 'PC2_Magnitude'
output_mag_dir = os.path.join(base_output_dir, "analise_magnitude")

if pc1_mag_col in df_full.columns and pc2_mag_col in df_full.columns:
    run_full_analysis(df_full, pc1_mag_col, pc2_mag_col, "Magnitude", output_mag_dir)
else:
    print(f"AVISO: Colunas para PCA de Magnitude ('{pc1_mag_col}', '{pc2_mag_col}') não encontradas. Análise pulada.")


# <--- Continuação do código ---

# --- ANÁLISE 2: PCA DE VARIABILIDADE ---
# Define as colunas de score para a PCA de Variabilidade
pc1_var_col = 'PC1_Variability'
pc2_var_col = 'PC2_Variability'
output_var_dir = os.path.join(base_output_dir, "analise_variabilidade")

if pc1_var_col in df_full.columns and pc2_var_col in df_full.columns:
    run_full_analysis(df_full, pc1_var_col, pc2_var_col, "Variability", output_var_dir)
else:
    print(f"AVISO: Colunas para PCA de Variabilidade ('{pc1_var_col}', '{pc2_var_col}') não encontradas. Análise pulada.")


# --- ANÁLISE EXTRA: JUSTIFICATIVA DO HABITAT (opcional, mas útil) ---
print("\n" + "="*60)
print("  ANÁLISE EXTRA: Justificativa do Fator HABITAT")
print("="*60)

output_justificativa_dir = os.path.join(base_output_dir, "justificativa_habitat")
os.makedirs(output_justificativa_dir, exist_ok=True)

# Usa o DataFrame completo carregado no início
df_valid_justificativa = df_full.dropna(subset=['Depth_m', 'mean_DLI_local', 'HAB']).copy()

if len(df_valid_justificativa) > 3 and df_valid_justificativa['HAB'].nunique() > 1:
    # ANOVA para Depth vs HAB
    grouped_depth_hab = [g['Depth_m'].values for _, g in df_valid_justificativa.groupby('HAB')]
    anova_depth_hab = f_oneway(*grouped_depth_hab)

    # ANOVA para DLI vs HAB
    grouped_dli_hab = [g['mean_DLI_local'].values for _, g in df_valid_justificativa.groupby('HAB')]
    anova_dli_hab = f_oneway(*grouped_dli_hab)
    
    # Boxplots
    fig_direct, axes = plt.subplots(1, 2, figsize=(12, 5))
    sns.boxplot(ax=axes[0], x='HAB', y='Depth_m', data=df_valid_justificativa, palette='mako')
    axes[0].set_title(f'Profundidade por HABITAT\nANOVA p={anova_depth_hab.pvalue:.4f}')
    axes[0].set_ylabel('Profundidade (m)')
    
    sns.boxplot(ax=axes[1], x='HAB', y='mean_DLI_local', data=df_valid_justificativa, palette='rocket')
    axes[1].set_title(f'DLI Médio por HABITAT\nANOVA p={anova_dli_hab.pvalue:.4f}')
    axes[1].set_ylabel('DLI Médio (mol/m²/dia)')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_justificativa_dir, "boxplots_Fisico_vs_HAB.png"), dpi=300)
    plt.close()

    # Salvar resultados numéricos da justificativa
    output_txt_just = os.path.join(output_justificativa_dir, "resultados_justificativa_habitat.txt")
    with open(output_txt_just, 'w') as f:
        f.write("="*50 + "\n")
        f.write("ANÁLISE DO LINK DIRETO ENTRE FATORES FÍSICOS E HABITAT\n")
        f.write("="*50 + "\n\n")
        f.write(f"  Profundidade (Depth_m) vs HAB: ANOVA F={anova_depth_hab.statistic:.3f}, p={anova_depth_hab.pvalue:.4f}\n")
        f.write(f"  DLI Médio (mean_DLI_local) vs HAB: ANOVA F={anova_dli_hab.statistic:.3f}, p={anova_dli_hab.pvalue:.4f}\n")
    print(f"Análise de justificativa do habitat concluída. Resultados em: {output_justificativa_dir}")
else:
    print("Dados insuficientes ou apenas um tipo de habitat para a análise de justificativa.")


print("\nAnálise estatística completa finalizada!")