"""
PCA Figure Generation Script v2 - Uses Pre-computed Data
=========================================================

Generates publication-quality PCA figures from existing computed data:
1. Magnitude PCA: Single figure (means don't change across CV windows)
2. Variability PCA: 3-row figure (one per CV window: 2, 30, all)

Changes from v1:
- Abandoned bubble plots - uniform marker sizes
- Professional legend placement
- Consistent axis orientation (SST positive on PC1)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

output_base = r"C:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_RESULTS"
output_dir = os.path.join(output_base, "#####output_local_PCA_v2_FINAL")
os.makedirs(output_dir, exist_ok=True)

cv_windows = ["2", "30", "all"]

data_dirs = {
    "2": os.path.join(output_base, "#####output_local_PCA_CV_2_FINAL"),
    "30": os.path.join(output_base, "#####output_local_PCA_CV_30_FINAL"),
    "all": os.path.join(output_base, "#####output_local_PCA_CV_all_FINAL"),
}

print("=" * 60)
print("PCA Figure Generation v2 - Using Pre-computed Data")
print("=" * 60)


def load_data_and_loadings(cv_window):
    """Load pre-computed PCA scores and loadings for a given CV window."""
    data_dir = data_dirs[cv_window]

    if cv_window == "2":
        data_file = os.path.join(
            data_dir, "dados_consolidados_com_scores_das_duas_PCAs_CV2.csv"
        )
    elif cv_window == "30":
        data_file = os.path.join(
            data_dir, "dados_consolidados_com_scores_das_duas_PCAs_CV30.csv"
        )
    else:
        data_file = os.path.join(
            data_dir, "dados_consolidados_com_scores_das_duas_PCAs_CVall.csv"
        )

    df = pd.read_csv(data_file, sep=";", decimal=",")

    for col in df.select_dtypes(include=["object"]).columns:
        df[col] = df[col].astype(str).str.strip()

    if cv_window == "2":
        loadings_mag_file = os.path.join(data_dir, "loadings_PCA_Magnitude_CV_2.csv")
        loadings_var_file = os.path.join(data_dir, "loadings_PCA_Variability_CV_2.csv")
    elif cv_window == "30":
        loadings_mag_file = os.path.join(data_dir, "loadings_PCA_Magnitude_CV_30.csv")
        loadings_var_file = os.path.join(data_dir, "loadings_PCA_Variability_CV_30.csv")
    else:
        loadings_mag_file = os.path.join(data_dir, "loadings_PCA_Magnitude_CV_all.csv")
        loadings_var_file = os.path.join(
            data_dir, "loadings_PCA_Variability_CV_all.csv"
        )

    with open(loadings_mag_file, "r") as f:
        first_line = f.readline()
    sep = ";" if ";" in first_line else ","
    decimal = "," if sep == ";" else "."

    loadings_mag = pd.read_csv(loadings_mag_file, index_col=0, sep=sep, decimal=decimal)
    loadings_var = pd.read_csv(loadings_var_file, index_col=0, sep=sep, decimal=decimal)

    return df, loadings_mag, loadings_var


def get_visual_attributes(df):
    """Extract consistent visual attributes across all figures."""
    unique_reefs = df["Reef_name"].unique()
    cmap = plt.get_cmap("tab10")
    color_map = {r: cmap(i % 10) for i, r in enumerate(unique_reefs)}

    unique_habitats = df["HAB"].unique()
    habitat_shapes = ["o", "s", "^", "D", "v", "<", ">"]
    shape_map = {
        hab: habitat_shapes[i % len(habitat_shapes)]
        for i, hab in enumerate(unique_habitats)
    }

    arch_styles = {
        "inner": {"edgecolor": "black", "linewidth": 2.0},
        "outer": {"edgecolor": "darkgrey", "linewidth": 0.75},
    }

    return color_map, shape_map, arch_styles, unique_reefs, unique_habitats


def create_legend_elements(
    color_map, shape_map, arch_styles, unique_reefs, unique_habitats
):
    """Create legend elements for all three categories."""
    legend_elements_reef = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label=reef,
            markersize=10,
            markerfacecolor=color_map[reef],
        )
        for reef in unique_reefs
    ]

    legend_elements_habitat = [
        Line2D(
            [0],
            [0],
            marker=shape_map[hab],
            color="grey",
            label=hab,
            linestyle="None",
            markersize=10,
        )
        for hab in unique_habitats
    ]

    legend_elements_arch = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="Inner Arc",
            markersize=10,
            markeredgecolor=arch_styles["inner"]["edgecolor"],
            markeredgewidth=arch_styles["inner"]["linewidth"],
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label="Outer Arc",
            markersize=10,
            markeredgecolor=arch_styles["outer"]["edgecolor"],
            markeredgewidth=arch_styles["outer"]["linewidth"],
        ),
    ]

    return legend_elements_reef, legend_elements_habitat, legend_elements_arch


def plot_pca_scores(
    ax, df, pc1_col, pc2_col, color_map, shape_map, arch_styles, marker_size=100
):
    """Plot PCA scores with uniform marker sizes."""
    for _, row in df.iterrows():
        if pd.notna(row[pc1_col]) and pd.notna(row[pc2_col]):
            arch_val = str(row["Arch"]).lower().strip()
            style = arch_styles.get(arch_val, {"edgecolor": "grey", "linewidth": 0.5})
            ax.scatter(
                row[pc1_col],
                row[pc2_col],
                color=color_map.get(row["Reef_name"]),
                marker=shape_map.get(row["HAB"]),
                s=marker_size,
                alpha=0.8,
                **style,
            )
    ax.axhline(0, color="grey", lw=0.5, zorder=0)
    ax.axvline(0, color="grey", lw=0.5, zorder=0)
    ax.grid(True, linestyle="--", alpha=0.6, zorder=0)


def plot_loadings(ax, loadings_df, scale=2.0):
    """Plot loading vectors with arrows and labels."""
    ax.axhline(0, color="grey", lw=0.5, zorder=0)
    ax.axvline(0, color="grey", lw=0.5, zorder=0)

    for i, var in enumerate(loadings_df.index):
        x = loadings_df.loc[var, "PC1"] * scale
        y = loadings_df.loc[var, "PC2"] * scale
        ax.arrow(
            0, 0, x, y, head_width=0.08, head_length=0.12, fc="red", ec="red", zorder=3
        )
        ax.text(
            x * 1.15,
            y * 1.15,
            var,
            color="black",
            ha="center",
            va="center",
            fontsize=11,
            fontweight="bold",
            zorder=4,
        )

    ax.set_xlim(-2.5, 2.5)
    ax.set_ylim(-2.5, 2.5)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linestyle="--", alpha=0.6, zorder=0)


def compute_explained_variance(df, pc1_col, pc2_col):
    """Compute explained variance from PCA scores."""
    pc1_var = df[pc1_col].var()
    pc2_var = df[pc2_col].var()
    total_var = pc1_var + pc2_var
    return [pc1_var / total_var * 100, pc2_var / total_var * 100]


print("\n--- Loading pre-computed data ---")
all_data = {}
for cv in cv_windows:
    df, loadings_mag, loadings_var = load_data_and_loadings(cv)
    all_data[cv] = {
        "df": df,
        "loadings_mag": loadings_mag,
        "loadings_var": loadings_var,
    }
    print(f"  CV {cv}: {len(df)} sites loaded")

df_reference = all_data["2"]["df"]
color_map, shape_map, arch_styles, unique_reefs, unique_habitats = (
    get_visual_attributes(df_reference)
)
legend_elements_reef, legend_elements_habitat, legend_elements_arch = (
    create_legend_elements(
        color_map, shape_map, arch_styles, unique_reefs, unique_habitats
    )
)

print("\n--- Creating Magnitude PCA Figure ---")

fig_mag, axes_mag = plt.subplots(
    1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1, 1]}
)
fig_mag.suptitle("Environmental Magnitude PCA", fontsize=18, fontweight="bold", y=0.98)

df_cv2 = all_data["2"]["df"]
loadings_mag = all_data["2"]["loadings_mag"]

sst_loading = loadings_mag.loc["log(sst_mean+1)", "PC1"]
if sst_loading < 0:
    print("  Flipping Magnitude PC1 sign for consistency")
    for cv in cv_windows:
        all_data[cv]["df"]["PC1_Magnitude"] = -all_data[cv]["df"]["PC1_Magnitude"]
        all_data[cv]["loadings_mag"]["PC1"] = -all_data[cv]["loadings_mag"]["PC1"]

ax1 = axes_mag[0]
plot_pca_scores(
    ax1,
    df_cv2,
    "PC1_Magnitude",
    "PC2_Magnitude",
    color_map,
    shape_map,
    arch_styles,
    marker_size=120,
)

var_exp = compute_explained_variance(df_cv2, "PC1_Magnitude", "PC2_Magnitude")
ax1.set_xlabel(f"PC1 ({var_exp[0]:.1f}%)", fontsize=14, fontweight="bold")
ax1.set_ylabel(f"PC2 ({var_exp[1]:.1f}%)", fontsize=14, fontweight="bold")
ax1.set_title("Site Ordination", fontsize=15, fontweight="bold")

ax2 = axes_mag[1]
plot_loadings(ax2, loadings_mag)
ax2.set_xlabel("Contribution to PC1", fontsize=14, fontweight="bold")
ax2.set_ylabel("Contribution to PC2", fontsize=14, fontweight="bold")
ax2.set_title("Variable Loadings", fontsize=15, fontweight="bold")

leg1 = fig_mag.legend(
    handles=legend_elements_reef,
    title="Reef",
    loc="center left",
    bbox_to_anchor=(1.0, 0.72),
    fontsize=10,
    title_fontsize=11,
)
leg2 = fig_mag.legend(
    handles=legend_elements_habitat,
    title="Habitat",
    loc="center left",
    bbox_to_anchor=(1.0, 0.42),
    fontsize=10,
    title_fontsize=11,
)
leg3 = fig_mag.legend(
    handles=legend_elements_arch,
    title="Arc",
    loc="center left",
    bbox_to_anchor=(1.0, 0.18),
    fontsize=10,
    title_fontsize=11,
)

plt.tight_layout(rect=[0, 0, 0.92, 0.95])
output_mag = os.path.join(output_dir, "PCA_Magnitude_publication_v2.png")
plt.savefig(output_mag, dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig_mag)
print(f"Magnitude figure saved: {output_mag}")

print("\n--- Creating Variability PCA Figure ---")

n_windows = len(cv_windows)
fig_var, axes_var = plt.subplots(
    n_windows, 2, figsize=(14, 6 * n_windows), gridspec_kw={"width_ratios": [1, 1]}
)
fig_var.suptitle(
    "Environmental Variability PCA (Multiple CV Windows)",
    fontsize=18,
    fontweight="bold",
    y=0.995,
)

for row_idx, cv in enumerate(cv_windows):
    df_cv = all_data[cv]["df"]
    loadings_var = all_data[cv]["loadings_var"]

    sst_loading_var = (
        loadings_var.iloc[0]["PC1"]
        if "sst" in loadings_var.index[0].lower()
        else loadings_var.iloc[0]["PC1"]
    )

    ax_scores = axes_var[row_idx, 0]
    plot_pca_scores(
        ax_scores,
        df_cv,
        "PC1_Variability",
        "PC2_Variability",
        color_map,
        shape_map,
        arch_styles,
        marker_size=120,
    )

    var_exp_var = compute_explained_variance(
        df_cv, "PC1_Variability", "PC2_Variability"
    )
    ax_scores.set_xlabel(f"PC1 ({var_exp_var[0]:.1f}%)", fontsize=14, fontweight="bold")
    ax_scores.set_ylabel(f"PC2 ({var_exp_var[1]:.1f}%)", fontsize=14, fontweight="bold")
    ax_scores.set_title(
        f"CV Window: {cv} - Site Ordination", fontsize=14, fontweight="bold"
    )

    ax_loadings = axes_var[row_idx, 1]
    plot_loadings(ax_loadings, loadings_var)
    ax_loadings.set_xlabel("Contribution to PC1", fontsize=14, fontweight="bold")
    ax_loadings.set_ylabel("Contribution to PC2", fontsize=14, fontweight="bold")
    ax_loadings.set_title(
        f"CV Window: {cv} - Variable Loadings", fontsize=14, fontweight="bold"
    )

leg1 = fig_var.legend(
    handles=legend_elements_reef,
    title="Reef",
    loc="center left",
    bbox_to_anchor=(1.0, 0.72),
    fontsize=10,
    title_fontsize=11,
)
leg2 = fig_var.legend(
    handles=legend_elements_habitat,
    title="Habitat",
    loc="center left",
    bbox_to_anchor=(1.0, 0.42),
    fontsize=10,
    title_fontsize=11,
)
leg3 = fig_var.legend(
    handles=legend_elements_arch,
    title="Arc",
    loc="center left",
    bbox_to_anchor=(1.0, 0.18),
    fontsize=10,
    title_fontsize=11,
)

plt.tight_layout(rect=[0, 0, 0.92, 0.97])
output_var = os.path.join(output_dir, "PCA_Variability_publication_v2.png")
plt.savefig(output_var, dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig_var)
print(f"Variability figure saved: {output_var}")

print("\n--- Saving loadings files ---")
for cv in cv_windows:
    loadings_mag_out = os.path.join(
        output_dir, f"loadings_PCA_Magnitude_CV_{cv}_v2.csv"
    )
    loadings_var_out = os.path.join(
        output_dir, f"loadings_PCA_Variability_CV_{cv}_v2.csv"
    )
    all_data[cv]["loadings_mag"].to_csv(loadings_mag_out)
    all_data[cv]["loadings_var"].to_csv(loadings_var_out)
    print(f"  Saved: loadings_PCA_*_CV_{cv}_v2.csv")

print("\n" + "=" * 60)
print("PROCESS COMPLETE")
print("=" * 60)
print(f"Output directory: {output_dir}")
print("Generated files:")
for f in os.listdir(output_dir):
    if f.endswith(".png") or f.endswith(".csv"):
        print(f"  - {f}")
print("=" * 60)
