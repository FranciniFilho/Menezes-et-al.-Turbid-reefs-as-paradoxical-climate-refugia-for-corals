##### Calcula RGR, médias de saúde, realiza PCAs de Saúde e Interações
##### Add growth boxplot
#### New PCAs with loadings

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import math  # Importado para ajudar a organizar as figuras compostas

# ---------------------------
# 1. CONFIGURAÇÕES E LEITURA DOS DADOS
# ---------------------------
print("--- Iniciando Script 2: Análise Biológica (Versão Aprimorada) ---")
base_path = "#######FINAL_DATA/01_raw_biological"
file_vitality = os.path.join(base_path, "Vitality_and_size_new.xlsx")
output_dir = r"#######FINAL_DATA/04_health_growth_PCA"
os.makedirs(output_dir, exist_ok=True)
print(f"Diretório de saída definido para: {output_dir}")
try:
    vitality_data = pd.read_excel(file_vitality)
except FileNotFoundError as e:
    print(f"ERRO: Arquivo não encontrado - {e}. Verifique o caminho.")
    exit()


def standardize_df(df):
    df.columns = df.columns.str.strip().str.upper()
    for col in ["SITE", "HAB", "REEF"]:
        if col in df.columns:
            df[col] = df[col].astype(str).str.strip().str.upper()
    return df


vitality_data = standardize_df(vitality_data)
vitality_data["SITE_COL"] = (
    vitality_data["SITE"] + "_" + vitality_data["COL"].astype(str)
)

# ---------------------------
# 2. CÁLCULO DAS MÉTRICAS DE DESEMPENHO POR COLÔNIA
# ---------------------------
print("\nCalculando métricas de desempenho por colônia (RGR, médias, etc.)...")


def compute_rgr(group):
    group = group.sort_values("YEAR")
    if len(group["YEAR"].unique()) < 2:
        return np.nan
    ar_inicial = group.iloc[0]["AR_TOTAL"]
    ar_final = group.iloc[-1]["AR_TOTAL"]
    ano_inicial = group.iloc[0]["YEAR"]
    ano_final = group.iloc[-1]["YEAR"]
    if ar_inicial <= 0 or ar_final <= 0 or ano_final == ano_inicial:
        return np.nan
    return (np.log(ar_final) - np.log(ar_inicial)) / (ano_final - ano_inicial)


def get_ar_total_2006(group):
    if 2006 in group["YEAR"].values:
        return group[group["YEAR"] == 2006]["AR_TOTAL"].iloc[0]
    return np.nan


colony_metrics = (
    vitality_data.groupby("SITE_COL")
    .apply(
        lambda g: pd.Series(
            {
                "RGR": compute_rgr(g),
                "AR_TOTAL_INICIAL": get_ar_total_2006(g),
                "MEAN_AR_TOTAL": g["AR_TOTAL"].mean(),
                "MEAN_HEALTH": g["HEALTH %"].mean(),
                "MEAN_BLEACH": g["BLEACHING %"].mean(),
                "MEAN_DEAD_TISSUE": g["DEAD %"].mean(),
            }
        )
    )
    .reset_index()
)
for col in ["MEAN_HEALTH", "MEAN_BLEACH", "MEAN_DEAD_TISSUE"]:
    if col in colony_metrics.columns:
        colony_metrics[col] = colony_metrics[col].clip(lower=0)
print(
    "Métricas calculadas (incluindo tecido morto) e valores negativos zerados com sucesso."
)

# ---------------------------
# 3. PCA DA SAÚDE DOS CORAIS (POR COLÔNIA-ANO)
# ---------------------------
print("\nRealizando PCA da Saúde dos Corais (observações por colônia-ano)...")
health_cols = ["HEALTH %", "BLEACHING %", "DEAD %"]

# Manter observações por colônia-ano (não agregar por média)
health_props = vitality_data[
    ["SITE_COL", "YEAR", "SITE", "HAB", "REEF"] + health_cols
].copy()
health_props = health_props.dropna(subset=health_cols)

n_obs = len(health_props)
n_colonies = health_props["SITE_COL"].nunique()
n_years = health_props["YEAR"].nunique()
print(f"  Dados para PCA: {n_obs} observações ({n_colonies} colônias, {n_years} anos)")
print(f"  Anos presentes: {sorted(health_props['YEAR'].unique())}")

health_props_scaled = StandardScaler().fit_transform(health_props[health_cols])
pca_health = PCA(n_components=2)
health_scores = pca_health.fit_transform(health_props_scaled)
df_scores_health = health_props[["SITE_COL", "YEAR", "SITE", "HAB", "REEF"]].copy()
df_scores_health["HEALTH_PC1"] = health_scores[:, 0]
df_scores_health["HEALTH_PC2"] = health_scores[:, 1]
df_scores_health = pd.merge(
    df_scores_health,
    health_props[["SITE_COL", "YEAR"] + health_cols],
    on=["SITE_COL", "YEAR"],
)
scores_path_health = os.path.join(output_dir, "scores_PCA_Saude.xlsx")
df_scores_health.to_excel(scores_path_health, index=False)
print(f"Scores da PCA de Saúde salvos em: {scores_path_health}")
print(f"  Colunas: {list(df_scores_health.columns)}")
loadings_df_health = pd.DataFrame(
    pca_health.components_.T, columns=["PC1", "PC2"], index=health_cols
)
loadings_path_health = os.path.join(output_dir, "loadings_PCA_Saude.csv")
loadings_df_health.to_csv(loadings_path_health)
print(f"Coral Health PCA loadings saved to: {loadings_path_health}")

# Export explained variance for R script precision
exp_var_health_df = pd.DataFrame(
    {
        "PC": ["PC1", "PC2"],
        "explained_variance_ratio": pca_health.explained_variance_ratio_,
    }
)
exp_var_health_df.to_csv(
    os.path.join(output_dir, "explained_variance_PCA_Saude.csv"), index=False
)


# ---------------------------
# 4. PCA DAS INTERAÇÕES LOCAIS (POR SITE-HAB-YEAR)
# ---------------------------
print("\nRealizando PCA das Interações Locais (por SITE-HAB-YEAR)...")
cols_interactions = [
    "SUR_TURF %",
    "SUR_CCA %",
    "SUR_CYANO %",
    "SUR_DICTYOTA %",
    "SUR_OTHMACR %",
    "SUR_PALYTHOA %",
    "SUR_CORAL %",
    "SUR_SAND %",
    "SUR_NON-BIOTIC %",
]
interaction_means_agg = (
    vitality_data.groupby(["SITE", "HAB", "YEAR", "REEF"])[cols_interactions]
    .mean(numeric_only=True)
    .reset_index()
)
print(f"  Dados para PCA: {len(interaction_means_agg)} observações SITE-HAB-YEAR")
interaction_vars = {
    "SUR_TURF": interaction_means_agg["SUR_TURF %"],
    "SUR_CCA": interaction_means_agg["SUR_CCA %"],
    "SUR_MACROALGAE": interaction_means_agg["SUR_DICTYOTA %"]
    + interaction_means_agg["SUR_OTHMACR %"],
    "SUR_CYANO": interaction_means_agg["SUR_CYANO %"],
    "SUR_PALYTHOA": interaction_means_agg["SUR_PALYTHOA %"],
    "SUR_ABIOTIC": interaction_means_agg["SUR_SAND %"]
    + interaction_means_agg["SUR_NON-BIOTIC %"],
}
new_interactions = pd.DataFrame(interaction_vars)
new_interactions.insert(0, "REEF", interaction_means_agg["REEF"])
new_interactions.insert(0, "YEAR", interaction_means_agg["YEAR"])
new_interactions.insert(0, "HAB", interaction_means_agg["HAB"])
new_interactions.insert(0, "SITE", interaction_means_agg["SITE"])
X_interactions = new_interactions.drop(columns=["SITE", "HAB", "YEAR", "REEF"]).fillna(
    0
)
X_interactions_scaled = StandardScaler().fit_transform(X_interactions)
pca_interactions = PCA(n_components=2)
interaction_scores = pca_interactions.fit_transform(X_interactions_scaled)
df_scores_interactions = new_interactions[["SITE", "HAB", "YEAR", "REEF"]].copy()
df_scores_interactions["PC1_INTERACAO"] = interaction_scores[:, 0]
df_scores_interactions["PC2_INTERACAO"] = interaction_scores[:, 1]
df_scores_interactions = pd.merge(
    df_scores_interactions, new_interactions, on=["SITE", "HAB", "YEAR", "REEF"]
)
scores_path_interactions = os.path.join(output_dir, "scores_PCA_Interacoes.xlsx")
df_scores_interactions.to_excel(scores_path_interactions, index=False)
print(f"Scores da PCA de Interações salvos em: {scores_path_interactions}")
print(f"  Colunas: {list(df_scores_interactions.columns)}")
loadings_df_interactions = pd.DataFrame(
    pca_interactions.components_.T, columns=["PC1", "PC2"], index=X_interactions.columns
)
loadings_path_interactions = os.path.join(output_dir, "loadings_PCA_Interacoes.csv")
loadings_df_interactions.to_csv(loadings_path_interactions)
print(f"Local Interactions PCA loadings saved to: {loadings_path_interactions}")

# Export explained variance for R script precision
exp_var_int_df = pd.DataFrame(
    {
        "PC": ["PC1", "PC2"],
        "explained_variance_ratio": pca_interactions.explained_variance_ratio_,
    }
)
exp_var_int_df.to_csv(
    os.path.join(output_dir, "explained_variance_PCA_Interacoes.csv"), index=False
)


# =============================================================================
# 4.5. FUNÇÃO PARA PLOTAGEM COMPOSTA DA PCA (ADAPTADA DO SCRIPT DE ANÁLISE LOCAL) <--- NOVA FUNÇÃO AQUI
# =============================================================================
def create_composite_pca_figure(pca_results, output_filename):
    if pca_results is None:
        print(
            f"Não há resultados de PCA para gerar a figura {output_filename}. Pulando."
        )
        return
    df_scores = pca_results["df_scores"]
    loadings = pca_results["loadings"]
    explained_variance = pca_results["explained_variance"]
    pc1_col, pc2_col = pca_results["pc1_col"], pca_results["pc2_col"]
    title_suffix = pca_results["title_suffix"]
    unique_reefs_plot = sorted(df_scores["REEF"].unique())
    cmap_plot = plt.get_cmap("tab10")
    color_map_plot = {reef: cmap_plot(i) for i, reef in enumerate(unique_reefs_plot)}
    unique_habitats_plot = sorted(df_scores["HAB"].unique())
    habitat_shapes_plot = ["o", "s", "^", "D", "v", "<", ">"]
    shape_map_plot = {
        hab: habitat_shapes_plot[i % len(habitat_shapes_plot)]
        for i, hab in enumerate(unique_habitats_plot)
    }
    if "YEAR" in df_scores.columns:
        unique_years_plot = sorted(df_scores["YEAR"].unique())
        alpha_vals = np.linspace(0.4, 0.95, len(unique_years_plot))
        alpha_map_plot = {
            year: alpha_vals[i] for i, year in enumerate(unique_years_plot)
        }
    else:
        unique_years_plot = []
        alpha_map_plot = None
    fig_title = f"Principal Component Analysis Summary - {title_suffix}"
    fig, axes = plt.subplots(
        1, 3, figsize=(24, 7), gridspec_kw={"width_ratios": [1.2, 1, 0.8]}
    )
    fig.suptitle(fig_title, fontsize=20, y=1.02)
    ax1 = axes[0]
    if alpha_map_plot is not None:
        for year in unique_years_plot:
            df_year = df_scores[df_scores["YEAR"] == year]
            sns.scatterplot(
                data=df_year,
                x=pc1_col,
                y=pc2_col,
                hue="REEF",
                style="HAB",
                palette=color_map_plot,
                markers=shape_map_plot,
                s=100,
                alpha=alpha_map_plot[year],
                edgecolor="k",
                ax=ax1,
                legend=False,
            )
    else:
        sns.scatterplot(
            data=df_scores,
            x=pc1_col,
            y=pc2_col,
            hue="REEF",
            style="HAB",
            palette=color_map_plot,
            markers=shape_map_plot,
            s=100,
            alpha=0.8,
            edgecolor="k",
            ax=ax1,
            legend=False,
        )
    ax1.set_xlabel(f"PC1 ({explained_variance[0]:.1f}%)", fontsize=14)
    ax1.set_ylabel(f"PC2 ({explained_variance[1]:.1f}%)", fontsize=14)
    ax1.set_title("Site Scores", fontsize=16)
    ax1.grid(True, linestyle="--", alpha=0.6)
    ax1.axhline(0, color="grey", lw=0.5)
    ax1.axvline(0, color="grey", lw=0.5)
    legend_elements_color = [
        plt.Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label=reef,
            markersize=10,
            markerfacecolor=color_map_plot[reef],
        )
        for reef in unique_reefs_plot
    ]
    legend_elements_shape = [
        plt.Line2D(
            [0],
            [0],
            marker=shape_map_plot[hab],
            color="grey",
            label=hab,
            linestyle="None",
            markersize=10,
        )
        for hab in unique_habitats_plot
    ]
    if alpha_map_plot is not None:
        legend_elements_year = [
            plt.Line2D(
                [0],
                [0],
                marker="o",
                color="w",
                label=str(year),
                markersize=10,
                markerfacecolor="grey",
                alpha=alpha_map_plot[year],
            )
            for year in unique_years_plot
        ]
        fig.legend(
            title="Year",
            handles=legend_elements_year,
            loc="center left",
            bbox_to_anchor=(0.91, 0.25),
        )
    fig.legend(
        title="Reef",
        handles=legend_elements_color,
        loc="center left",
        bbox_to_anchor=(0.91, 0.65),
    )
    fig.legend(
        title="Habitat",
        handles=legend_elements_shape,
        loc="center left",
        bbox_to_anchor=(0.91, 0.45),
    )
    ax2 = axes[1]
    ax2.axhline(0, color="grey", lw=0.5)
    ax2.axvline(0, color="grey", lw=0.5)
    for i, var in enumerate(loadings.index):
        pc1_val = loadings["PC1"].iloc[i]
        pc2_val = loadings["PC2"].iloc[i]
        ax2.arrow(
            0,
            0,
            pc1_val * 1.5,
            pc2_val * 1.5,
            head_width=0.05,
            head_length=0.1,
            fc="red",
            ec="red",
        )
        ax2.text(
            pc1_val * 1.7,
            pc2_val * 1.7,
            var,
            color="black",
            ha="center",
            va="center",
            fontsize=12,
        )
    ax2.set_xlim(-2, 2)
    ax2.set_ylim(-2, 2)
    ax2.set_xlabel("Contribution to PC1", fontsize=14)
    ax2.set_ylabel("Contribution to PC2", fontsize=14)
    ax2.set_title("Variable Loadings", fontsize=16)
    ax2.set_aspect("equal", adjustable="box")
    ax3 = axes[2]
    components = ["PC1", "PC2"]
    ax3.bar(components, explained_variance, color="skyblue", edgecolor="black")
    ax3.set_ylabel("Explained Variance (%)", fontsize=14)
    ax3.set_title("Component Importance", fontsize=16)
    ax3.set_ylim(0, 100)
    for i, v in enumerate(explained_variance):
        ax3.text(i, v + 2, f"{v:.1f}%", ha="center", color="black", fontsize=12)
    fig.subplots_adjust(right=0.9)
    plt.savefig(output_filename, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Composite PCA summary figure saved to: {output_filename}")


# ---------------------------
# 5. GERAÇÃO DOS GRÁFICOS (BUBBLE PLOTS) - VERSÃO CORRIGIDA
# ---------------------------
# ... (código dos bubble plots permanece o mesmo) ...
print("\nGerando gráficos de bubble plot...")
# ... (código omitido para brevidade) ...
print("Gráficos salvos com sucesso.")


# ------------------------------------------------------------------
# 5.5 GERAÇÃO DAS FIGURAS COMPOSTAS DE RESUMO DAS PCAS  <--- NOVA SEÇÃO AQUI
# ------------------------------------------------------------------
print(
    "\nGerando figuras compostas de resumo das PCAs (Ordenação + Loadings + Scree Plot)..."
)
health_pca_results = {
    "df_scores": df_scores_health,
    "loadings": loadings_df_health,
    "explained_variance": pca_health.explained_variance_ratio_ * 100,
    "pc1_col": "HEALTH_PC1",
    "pc2_col": "HEALTH_PC2",
    "title_suffix": "Coral Health",
}

create_composite_pca_figure(
    health_pca_results, os.path.join(output_dir, "figura_composta_PCA_Saude.png")
)
interactions_pca_results = {
    "df_scores": df_scores_interactions,
    "loadings": loadings_df_interactions,
    "explained_variance": pca_interactions.explained_variance_ratio_ * 100,
    "pc1_col": "PC1_INTERACAO",
    "pc2_col": "PC2_INTERACAO",
    "title_suffix": "Local Interactions",
}

create_composite_pca_figure(
    interactions_pca_results,
    os.path.join(output_dir, "figura_composta_PCA_Interacoes.png"),
)


# ---------------------------
# 6. SALVAR RESULTADOS BIOLÓGICOS
# ---------------------------
resultados_biologicos = colony_metrics.copy()
health_scores_by_colony = (
    df_scores_health.groupby("SITE_COL")[["HEALTH_PC1", "HEALTH_PC2"]]
    .mean()
    .reset_index()
)
meta_info = vitality_data[["SITE_COL", "SITE", "HAB", "REEF"]].drop_duplicates()
interaction_by_site_hab = (
    df_scores_interactions.groupby(["SITE", "HAB"])[["PC1_INTERACAO", "PC2_INTERACAO"]]
    .mean()
    .reset_index()
)

resultados_biologicos = pd.merge(
    resultados_biologicos, health_scores_by_colony, on="SITE_COL", how="left"
)
resultados_biologicos = pd.merge(
    resultados_biologicos, meta_info, on="SITE_COL", how="left"
)
resultados_biologicos = pd.merge(
    resultados_biologicos, interaction_by_site_hab, on=["SITE", "HAB"], how="left"
)
final_cols = [
    "SITE_COL",
    "SITE",
    "HAB",
    "REEF",
    "RGR",
    "AR_TOTAL_INICIAL",
    "MEAN_AR_TOTAL",
    "MEAN_HEALTH",
    "MEAN_BLEACH",
    "MEAN_DEAD_TISSUE",
    "HEALTH_PC1",
    "HEALTH_PC2",
    "PC1_INTERACAO",
    "PC2_INTERACAO",
]
final_cols_exist = [col for col in final_cols if col in resultados_biologicos.columns]
resultados_biologicos = resultados_biologicos[final_cols_exist]
output_file = os.path.join(output_dir, "resultados_biologicos_por_colonia.csv")
resultados_biologicos.to_csv(output_file, index=False)
print(f"\nResultados biológicos consolidados salvos em: {output_file}")

# ---------------------------
# 6.5. EXPORTAR DADOS DE SAÚDE POR COLÔNIA-ANO
# ---------------------------
print("\nExportando dados de saúde por colônia-ano (com YEAR)...")

health_yearly = pd.merge(
    vitality_data[
        ["SITE_COL", "YEAR", "SITE", "HAB", "REEF", "HEALTH %", "BLEACHING %", "DEAD %"]
    ],
    df_scores_health[["SITE_COL", "YEAR", "HEALTH_PC1", "HEALTH_PC2"]],
    on=["SITE_COL", "YEAR"],
    how="inner",
)

health_yearly = pd.merge(
    health_yearly,
    df_scores_interactions[["SITE", "HAB", "YEAR", "PC1_INTERACAO", "PC2_INTERACAO"]],
    on=["SITE", "HAB", "YEAR"],
    how="left",
)

print(f"  Total de observações: {len(health_yearly)}")
print(f"  Colônias: {health_yearly['SITE_COL'].nunique()}")
print(f"  Anos: {sorted(health_yearly['YEAR'].unique())}")

health_yearly_path = os.path.join(output_dir, "dados_saude_por_colonia_ano.csv")
health_yearly.to_csv(health_yearly_path, index=False)
print(f"Dados de saúde por colônia-ano salvos em: {health_yearly_path}")

health_yearly_excel_path = os.path.join(output_dir, "dados_saude_por_colonia_ano.xlsx")
health_yearly.to_excel(health_yearly_excel_path, index=False)
print(f"Dados de saúde por colônia-ano (Excel) salvos em: {health_yearly_excel_path}")

# ---------------------------
# 7. GERAÇÃO DE BOX PLOTS (RGR)  <--- CÓDIGO NOVO INSERIDO AQUI
# ---------------------------
print("\nGerando box plots para a Taxa de Crescimento Relativa (RGR)...")

# Define o tema/estilo do seaborn para ser consistente com seus outros gráficos
sns.set_theme(style="whitegrid", palette="muted")

# É uma boa prática remover valores nulos de RGR antes de plotar
rgr_data_for_plot = resultados_biologicos.dropna(subset=["RGR"])

# Cria a figura para o box plot
plt.figure(figsize=(12, 7))  # Aumentei um pouco a largura para acomodar a legenda
ax = sns.boxplot(
    data=rgr_data_for_plot,
    x="REEF",
    y="RGR",
    hue="HAB",
    order=sorted(
        rgr_data_for_plot["REEF"].unique()
    ),  # Garante a ordem alfabética dos recifes
    showfliers=False,  # Oculta os outliers, como no seu script de exemplo
)

# Adiciona títulos e rótulos
ax.set_title("Relative Growth Rate (RGR) by Reef and Habitat", fontsize=16)
ax.set_xlabel("Reef", fontsize=12)
ax.set_ylabel("Relative Growth Rate (RGR)", fontsize=12)


# Adiciona uma linha horizontal em y=0 para indicar crescimento zero
ax.axhline(0, color="red", linestyle="--", linewidth=1)

# Posiciona a legenda do lado de fora do gráfico para não obstruir os dados
ax.legend(title="Habitat", bbox_to_anchor=(1.02, 1), loc="upper left")

# Ajusta o layout para garantir que tudo (incluindo a legenda) caiba na imagem
plt.tight_layout()

# Salva a figura no diretório de saída
boxplot_path = os.path.join(output_dir, "boxplot_RGR_por_Reef_e_Hab.png")
plt.savefig(boxplot_path, dpi=300)
plt.close()  # Fecha a figura para liberar memória

print(f"Box plot de RGR salvo com sucesso em: {boxplot_path}")


print("\n--- Script 2 (Versão Aprimorada) concluído com sucesso! ---")
