### PCA LOCAL v2 - Figure Generation for Publication ###
#
# CHANGES FROM v1:
# 1. Single Magnitude PCA figure (means don't change across CV windows)
# 2. Composite Variability PCA figure with 3 rows (one per CV window)
# 3. Abandoned bubble plots - uniform marker sizes
# 4. Professional legend placement
# 5. Separate loading panels for each analysis
# 6. Consistent axis orientation (avoids sign inversion)

import os

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
os.environ["OMP_NUM_THREADS"] = "1"
import glob
import re
import logging
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import dask
import psutil
import gc

logging.basicConfig(
    filename="pca_local_v2_errors.log",
    filemode="w",
    level=logging.ERROR,
    format="%(asctime)s - %(levelname)s - %(message)s",
)

output_base = r"C:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_RESULTS"
output_dir = os.path.join(output_base, "#####output_local_PCA_v2_FINAL")
os.makedirs(output_dir, exist_ok=True)

dask.config.set(
    {
        "array.slicing.split_large_chunks": True,
        "scheduler": "single-threaded",
        "temporary_directory": os.path.join(output_dir, "dask_temp"),
    }
)
os.makedirs(os.path.join(output_dir, "dask_temp"), exist_ok=True)

cv_windows = ["2", "30", "all"]

print("=" * 60)
print("PCA LOCAL v2 - Publication-Quality Figure Generation")
print(f"CV Windows: {cv_windows}")
print("=" * 60)

sst_dir = r"F:\remote sensing\CRW_SST_FULL"
modis_dir = r"F:\remote sensing\MODIS_DATA_FULL"
chl_dir = r"F:\remote sensing\MODIS_DATA_FULL"
dhw_dir = r"F:\remote sensing\CRW_DHW_FULL"

sst_period = (2020, 2024)
dhw_period = (2020, 2024)
light_period = (2020, 2024)
chl_period = (2020, 2024)

sites_csv_file = r"C:\Users\rbfra\OneDrive\########CEBIMAR\####PROJETOS\#####Coral trade offs\sites_list_full.csv"
lat_min, lat_max = -20.5, -14.5
lon_min, lon_max = -40.5, -35.5
sst_pattern = "coraltemp_v3.1_*.nc"
dhw_pattern = "*.nc"
kd490_pattern = "AQUA_MODIS.*.L3m.DAY.KD.Kd_490.4km.nc"
par_pattern = "AQUA_MODIS.*.L3m.DAY.PAR.par.4km.nc"
chl_pattern = "AQUA_MODIS.*.L3m.DAY.CHL.chlor_a.4km.nc"
KDPAR_GATTUSO_A, KDPAR_GATTUSO_B, KDPAR_GATTUSO_C = 0.0665, 0.874, 0.00121
KDPAR_KD490_MIN_THRESHOLD = 0.001

try:
    df_sites_info = pd.read_csv(sites_csv_file, sep=";")
    df_sites_info.columns = [
        col.strip().replace(" ", "_") for col in df_sites_info.columns
    ]
    for col in df_sites_info.select_dtypes(["object"]).columns:
        df_sites_info[col] = df_sites_info[col].str.strip()
    df_sites_info["unique_id"] = df_sites_info["Site_name"] + "_" + df_sites_info["HAB"]
    sites = list(zip(df_sites_info["Latitude"], df_sites_info["Longitude"]))
    print(f"Sites loaded: {len(sites)}")
except Exception as e:
    logging.critical(f"Failed to load sites file: {e}")
    raise


def verify_file_integrity(filepath):
    if not os.path.exists(filepath) or os.path.getsize(filepath) == 0:
        return False
    try:
        with xr.open_dataset(filepath) as ds:
            _ = ds.attrs
        return True
    except Exception:
        return False


def load_satellite_data(
    pattern,
    base_dir,
    period,
    var_name_options,
    chunks={"time": 30, "lat": 100, "lon": 100},
):
    full_search_path = os.path.join(base_dir, pattern)
    print(f"\n--- Loading: {pattern} for {period} ---")
    all_files = sorted(glob.glob(full_search_path))
    if not all_files:
        logging.error(f"No files found for '{pattern}' in {base_dir}")
        return None
    candidate_paths = [
        f
        for f in all_files
        if (match := re.search(r"(\d{4})\d{4}", os.path.basename(f)))
        and period[0] <= int(match.group(1)) <= period[1]
    ]
    if not candidate_paths:
        logging.warning(f"No files found for '{pattern}' in period {period}")
        return None
    valid_paths = [f for f in candidate_paths if verify_file_integrity(f)]
    if not valid_paths:
        logging.error(f"No valid files for '{pattern}' in period {period}")
        return None
    print(f"Found {len(valid_paths)} valid files")

    def preprocess_with_time(ds):
        var_name = next((v for v in var_name_options if v in ds.data_vars), None)
        if var_name is None:
            return xr.Dataset()
        ds_subset = ds[[var_name]].astype("float32")
        if "time" not in ds_subset.coords:
            match = re.search(
                r"(\d{8})", os.path.basename(ds.encoding.get("source", ""))
            )
            if match:
                dt = pd.to_datetime(match.group(1), format="%Y%m%d")
                return ds_subset.expand_dims(time=[dt])
        return ds_subset

    try:
        ds_raw = xr.open_mfdataset(
            valid_paths,
            preprocess=preprocess_with_time,
            combine="by_coords",
            parallel=True,
            chunks=chunks,
            engine="netcdf4",
        )
        if not ds_raw.data_vars:
            print(f"WARNING: No data loaded for '{pattern}'")
            return None
        lat_slice = (
            slice(lat_max, lat_min)
            if ds_raw["lat"].values[0] > ds_raw["lat"].values[-1]
            else slice(lat_min, lat_max)
        )
        ds_final = ds_raw.sel(lat=lat_slice, lon=slice(lon_min, lon_max))
        return ds_final.sortby("time")[list(ds_final.data_vars)[0]]
    except Exception as e:
        logging.critical(f"Failed to concatenate files for '{pattern}': {e}")
        return None


def calculate_kdpar_gattuso(kd490_da):
    kd490_safe = xr.where(
        (kd490_da.isnull()) | (kd490_da <= KDPAR_KD490_MIN_THRESHOLD), np.nan, kd490_da
    )
    epsilon = 1e-9
    kdpar = (
        KDPAR_GATTUSO_A
        + KDPAR_GATTUSO_B * kd490_safe
        - KDPAR_GATTUSO_C / (kd490_safe + epsilon)
    )
    kdpar_final = xr.where((kdpar.isnull()) | (kdpar <= 0), np.nan, kdpar)
    kdpar_final.name = "kdpar"
    return kdpar_final


def extract_and_compute_site_metrics(
    da, sites_coords, metrics_to_calc, cv_windows_list
):
    if da is None:
        return {}
    if isinstance(da, xr.Dataset):
        da = da[list(da.data_vars)[0]]
    lats = xr.DataArray([s[0] for s in sites_coords], dims="site")
    lons = xr.DataArray([s[1] for s in sites_coords], dims="site")
    site_timeseries = da.sel(lat=lats, lon=lons, method="nearest").load()

    results = {}
    for metric in metrics_to_calc:
        print(f"  Computing metric: {metric}...")
        if metric == "mean":
            results["mean"] = site_timeseries.mean("time", skipna=True).values
        elif metric == "cv_all":
            mean_vals = site_timeseries.mean("time", skipna=True)
            std_vals = site_timeseries.std("time", skipna=True)
            results["cv_all"] = (std_vals / (mean_vals + 1e-9) * 100).values
            results["cv_all_segment_count"] = (
                site_timeseries.notnull().sum("time").values
            )
        elif metric.startswith("cv_"):
            try:
                window = int(metric.split("_")[1])
                min_p = max(2, int(window * 0.25))
                rolling_mean = site_timeseries.rolling(
                    time=window, min_periods=min_p, center=True
                ).mean()
                rolling_std = site_timeseries.rolling(
                    time=window, min_periods=min_p, center=True
                ).std()
                cv_ts = (rolling_std / (rolling_mean + 1e-9)) * 100
                results[metric] = cv_ts.mean("time", skipna=True).values
                results[f"{metric}_segment_count"] = cv_ts.notnull().sum("time").values
            except (ValueError, IndexError):
                print(f"    WARNING: Could not process metric '{metric}'. Skipping.")
                continue
        elif metric.startswith("prop_gt_"):
            threshold = float(metric.split("_")[-1])
            results[metric] = (
                (site_timeseries > threshold).mean("time", skipna=True).values
            )
    return results


def calculate_site_specific_dli(kd490_da, par_da, sites_info_df, cv_windows_list):
    if kd490_da is None or par_da is None:
        return {}
    lats = xr.DataArray(sites_info_df["Latitude"].values, dims="site")
    lons = xr.DataArray(sites_info_df["Longitude"].values, dims="site")
    depths = xr.DataArray(sites_info_df["Depth_m"].values, dims="site")

    kd490_points = kd490_da.sel(lat=lats, lon=lons, method="nearest").load()
    par_points = par_da.sel(lat=lats, lon=lons, method="nearest").load()
    kdpar_points = calculate_kdpar_gattuso(kd490_points)
    benthic_dli_timeseries = par_points * np.exp(-kdpar_points * depths)

    results = {}
    mean_dli = benthic_dli_timeseries.mean("time", skipna=True)
    results["mean_DLI_local"] = mean_dli.values

    std_dli = benthic_dli_timeseries.std("time", skipna=True)
    results["cv_DLI_local"] = (std_dli / (mean_dli + 1e-9) * 100).values
    results["cv_DLI_local_segment_count"] = (
        benthic_dli_timeseries.notnull().sum("time").values
    )

    for window in cv_windows_list:
        if window == "all":
            continue
        w = int(window)
        min_p = max(2, int(w * 0.25))
        rolling_mean = benthic_dli_timeseries.rolling(
            time=w, min_periods=min_p, center=True
        ).mean()
        rolling_std = benthic_dli_timeseries.rolling(
            time=w, min_periods=min_p, center=True
        ).std()
        cv_ts = (rolling_std / (rolling_mean + 1e-9)) * 100
        results[f"dli_cv_{w}"] = cv_ts.mean("time", skipna=True).values
        results[f"dli_cv_{w}_segment_count"] = cv_ts.notnull().sum("time").values

    return results


print("\n" + "=" * 60)
print("STEP 1: Loading Satellite Data")
print("=" * 60)
sst_all_ds = load_satellite_data(
    sst_pattern, sst_dir, sst_period, ["analysed_sst", "sea_surface_temperature", "sst"]
)
dhw_all_ds = load_satellite_data(
    dhw_pattern, dhw_dir, dhw_period, ["degree_heating_week"]
)
chl_all_ds = load_satellite_data(chl_pattern, chl_dir, chl_period, ["chlor_a"])
kd490_all_ds = load_satellite_data(kd490_pattern, modis_dir, light_period, ["Kd_490"])
par_all_ds = load_satellite_data(par_pattern, modis_dir, light_period, ["par"])

print("\n" + "=" * 60)
print("STEP 2: Computing Site Metrics")
print("=" * 60)
df_sites_local = df_sites_info.copy()

dli_metrics = calculate_site_specific_dli(
    kd490_all_ds, par_all_ds, df_sites_info, cv_windows
)
for key, values in dli_metrics.items():
    df_sites_local[key] = values
del kd490_all_ds, par_all_ds, dli_metrics
gc.collect()

sst_metrics_to_calculate = ["mean", "cv_all"] + [
    f"cv_{w}" for w in cv_windows if w != "all"
]
sst_metrics = extract_and_compute_site_metrics(
    sst_all_ds, sites, sst_metrics_to_calculate, cv_windows
)
for key, values in sst_metrics.items():
    df_sites_local[f"sst_{key}"] = values
del sst_all_ds, sst_metrics
gc.collect()

dhw_metrics = extract_and_compute_site_metrics(
    dhw_all_ds, sites, ["prop_gt_4", "prop_gt_8"], cv_windows
)
if "prop_gt_4" in dhw_metrics:
    df_sites_local["prop_DHW_gt4"] = dhw_metrics["prop_gt_4"]
if "prop_gt_8" in dhw_metrics:
    df_sites_local["prop_DHW_gt8"] = dhw_metrics["prop_gt_8"]
del dhw_all_ds, dhw_metrics
gc.collect()

chl_metrics_to_calculate = ["mean", "cv_all"] + [
    f"cv_{w}" for w in cv_windows if w != "all"
]
chl_metrics = extract_and_compute_site_metrics(
    chl_all_ds, sites, chl_metrics_to_calculate, cv_windows
)
for key, values in chl_metrics.items():
    df_sites_local[f"chl_{key}"] = values
del chl_all_ds, chl_metrics
gc.collect()

print("\n" + "=" * 60)
print("STEP 3: Running PCA Analyses")
print("=" * 60)


def run_pca_analysis(df_input, pca_vars_list, output_suffix, force_positive_pc1=None):
    """
    Execute PCA analysis and return scores, loadings, and metadata.

    force_positive_pc1: Variable name to ensure its loading is positive on PC1.
                        If the loading is negative, PC1 sign is flipped.
    """
    print(f"\n--- Running PCA: {output_suffix} ---")

    existing_vars = [var for var in pca_vars_list if var in df_input.columns]
    if len(existing_vars) < 2:
        print(
            f"WARNING: PCA for '{output_suffix}' skipped. Variables found: {existing_vars}"
        )
        return df_input, None

    print(f"Variables used: {existing_vars}")
    df_subset = df_input.dropna(subset=existing_vars)
    if len(df_subset) < 2:
        print(f"WARNING: PCA for '{output_suffix}' skipped. Insufficient data points.")
        return df_input, None

    data_for_pca = df_subset[existing_vars].copy()
    transformed_names = list(data_for_pca.columns)

    for i, col_name in enumerate(data_for_pca.columns):
        if col_name == "prop_DHW_gt4":
            data_for_pca[col_name] = np.sqrt(data_for_pca[col_name] + 0.5)
            transformed_names[i] = "sqrt(DHW>4)"
        elif col_name in ["chl_mean", "sst_mean", "mean_DLI_local"]:
            data_for_pca[col_name] = np.log1p(data_for_pca[col_name])
            transformed_names[i] = f"log({col_name.split('_')[0]}+1)"

    scaler = StandardScaler()
    data_matrix = scaler.fit_transform(data_for_pca)
    pca = PCA(n_components=2)
    scores = pca.fit_transform(data_matrix)

    if force_positive_pc1:
        var_idx = None
        for i, vname in enumerate(transformed_names):
            if force_positive_pc1.lower() in vname.lower():
                var_idx = i
                break
        if var_idx is not None and pca.components_[0, var_idx] < 0:
            print(
                f"  Flipping PC1 sign to ensure positive {force_positive_pc1} loading"
            )
            scores[:, 0] = -scores[:, 0]
            pca.components_[0, :] = -pca.components_[0, :]

    explained_variance = pca.explained_variance_ratio_ * 100

    loadings_df = pd.DataFrame(
        pca.components_.T, columns=["PC1", "PC2"], index=transformed_names
    )
    loadings_path = os.path.join(output_dir, f"loadings_PCA{output_suffix}.csv")
    loadings_df.to_csv(loadings_path)
    print(f"Loadings saved: {loadings_path}")

    df_scores = df_input.copy()
    df_scores.loc[df_subset.index, f"PC1{output_suffix}"] = scores[:, 0]
    df_scores.loc[df_subset.index, f"PC2{output_suffix}"] = scores[:, 1]

    pca_results = {
        "df_scores": df_scores,
        "loadings": loadings_df,
        "explained_variance": explained_variance,
        "pc1_col": f"PC1{output_suffix}",
        "pc2_col": f"PC2{output_suffix}",
        "title_suffix": output_suffix.replace("_", " ").strip(),
    }

    return df_scores, pca_results


df_processed = df_sites_local.copy()

pca_vars_magnitude = ["sst_mean", "mean_DLI_local", "chl_mean", "prop_DHW_gt4"]
df_processed, magnitude_results = run_pca_analysis(
    df_processed, pca_vars_magnitude, "_Magnitude", force_positive_pc1="SST"
)

variability_results_dict = {}
for cv_window in cv_windows:
    if cv_window == "all":
        cv_sst_col = "sst_cv_all"
        cv_dli_col = "cv_DLI_local"
        cv_chl_col = "chl_cv_all"
    else:
        cv_sst_col = f"sst_cv_{cv_window}"
        cv_dli_col = f"dli_cv_{cv_window}"
        cv_chl_col = f"chl_cv_{cv_window}"

    pca_vars_variability = [cv_sst_col, cv_dli_col, cv_chl_col]
    suffix = f"_Variability_CV{cv_window}"
    df_processed, var_results = run_pca_analysis(
        df_processed, pca_vars_variability, suffix, force_positive_pc1="SST"
    )
    if var_results is not None:
        variability_results_dict[cv_window] = var_results

final_scores_path = os.path.join(
    output_dir, "dados_consolidados_com_scores_PCAs_v2.xlsx"
)
df_processed.to_excel(final_scores_path, index=False)
print(f"\nConsolidated data saved: {final_scores_path}")

print("\n" + "=" * 60)
print("STEP 4: Generating Publication Figures")
print("=" * 60)


def get_visual_attributes(df_scores):
    """Extract consistent visual attributes across all figures."""
    unique_reefs = df_scores["Reef_name"].unique()
    cmap = plt.get_cmap("tab10")
    color_map = {r: cmap(i % 10) for i, r in enumerate(unique_reefs)}

    unique_habitats = df_scores["HAB"].unique()
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
    ax, df_scores, pc1_col, pc2_col, color_map, shape_map, arch_styles, marker_size=100
):
    """Plot PCA scores with uniform marker sizes."""
    for _, row in df_scores.iterrows():
        if pd.notna(row[pc1_col]) and pd.notna(row[pc2_col]):
            style = arch_styles.get(
                row["Arch"].lower(), {"edgecolor": "grey", "linewidth": 0.5}
            )
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
        x = loadings_df["PC1"][i] * scale
        y = loadings_df["PC2"][i] * scale
        ax.arrow(
            0, 0, x, y, head_width=0.08, head_length=0.12, fc="red", ec="red", zorder=3
        )

        label_offset = 0.15
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


def create_magnitude_figure(pca_results, output_filename):
    """
    Create publication-quality figure for Magnitude PCA.
    Layout: 1 row × 2 columns (Scores | Loadings)
    """
    if pca_results is None:
        print(f"No results for Magnitude figure. Skipping.")
        return

    print(f"\n--- Creating Magnitude Figure ---")

    df_scores = pca_results["df_scores"]
    loadings = pca_results["loadings"]
    explained_variance = pca_results["explained_variance"]
    pc1_col = pca_results["pc1_col"]
    pc2_col = pca_results["pc2_col"]

    fig, axes = plt.subplots(
        1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1, 1]}
    )
    fig.suptitle("Environmental Magnitude PCA", fontsize=18, fontweight="bold", y=0.98)

    color_map, shape_map, arch_styles, unique_reefs, unique_habitats = (
        get_visual_attributes(df_scores)
    )
    legend_elements_reef, legend_elements_habitat, legend_elements_arch = (
        create_legend_elements(
            color_map, shape_map, arch_styles, unique_reefs, unique_habitats
        )
    )

    ax1 = axes[0]
    plot_pca_scores(
        ax1,
        df_scores,
        pc1_col,
        pc2_col,
        color_map,
        shape_map,
        arch_styles,
        marker_size=120,
    )
    ax1.set_xlabel(
        f"PC1 ({explained_variance[0]:.1f}%)", fontsize=14, fontweight="bold"
    )
    ax1.set_ylabel(
        f"PC2 ({explained_variance[1]:.1f}%)", fontsize=14, fontweight="bold"
    )
    ax1.set_title("Site Ordination", fontsize=15, fontweight="bold")

    ax2 = axes[1]
    plot_loadings(ax2, loadings)
    ax2.set_xlabel("Contribution to PC1", fontsize=14, fontweight="bold")
    ax2.set_ylabel("Contribution to PC2", fontsize=14, fontweight="bold")
    ax2.set_title("Variable Loadings", fontsize=15, fontweight="bold")

    leg1 = fig.legend(
        handles=legend_elements_reef,
        title="Reef",
        loc="center left",
        bbox_to_anchor=(1.0, 0.72),
        fontsize=10,
        title_fontsize=11,
    )
    leg2 = fig.legend(
        handles=legend_elements_habitat,
        title="Habitat",
        loc="center left",
        bbox_to_anchor=(1.0, 0.42),
        fontsize=10,
        title_fontsize=11,
    )
    leg3 = fig.legend(
        handles=legend_elements_arch,
        title="Arc",
        loc="center left",
        bbox_to_anchor=(1.0, 0.18),
        fontsize=10,
        title_fontsize=11,
    )

    plt.tight_layout(rect=[0, 0, 0.92, 0.95])
    plt.savefig(output_filename, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Magnitude figure saved: {output_filename}")


def create_variability_figure(variability_results_dict, output_filename):
    """
    Create publication-quality figure for Variability PCA.
    Layout: 3 rows × 2 columns (one row per CV window)
    Each row: Scores | Loadings
    """
    if not variability_results_dict:
        print("No variability results. Skipping figure.")
        return

    print(f"\n--- Creating Variability Figure ---")

    n_windows = len(variability_results_dict)
    fig, axes = plt.subplots(
        n_windows, 2, figsize=(14, 6 * n_windows), gridspec_kw={"width_ratios": [1, 1]}
    )

    if n_windows == 1:
        axes = axes.reshape(1, -1)

    fig.suptitle(
        "Environmental Variability PCA (Multiple CV Windows)",
        fontsize=18,
        fontweight="bold",
        y=0.995,
    )

    color_map = None
    shape_map = None
    arch_styles = None
    unique_reefs = None
    unique_habitats = None

    for row_idx, (cv_window, pca_results) in enumerate(
        sorted(variability_results_dict.items())
    ):
        df_scores = pca_results["df_scores"]
        loadings = pca_results["loadings"]
        explained_variance = pca_results["explained_variance"]
        pc1_col = pca_results["pc1_col"]
        pc2_col = pca_results["pc2_col"]

        if color_map is None:
            color_map, shape_map, arch_styles, unique_reefs, unique_habitats = (
                get_visual_attributes(df_scores)
            )

        ax_scores = axes[row_idx, 0]
        plot_pca_scores(
            ax_scores,
            df_scores,
            pc1_col,
            pc2_col,
            color_map,
            shape_map,
            arch_styles,
            marker_size=120,
        )
        ax_scores.set_xlabel(
            f"PC1 ({explained_variance[0]:.1f}%)", fontsize=14, fontweight="bold"
        )
        ax_scores.set_ylabel(
            f"PC2 ({explained_variance[1]:.1f}%)", fontsize=14, fontweight="bold"
        )
        ax_scores.set_title(
            f"CV Window: {cv_window} - Site Ordination", fontsize=14, fontweight="bold"
        )

        ax_loadings = axes[row_idx, 1]
        plot_loadings(ax_loadings, loadings)
        ax_loadings.set_xlabel("Contribution to PC1", fontsize=14, fontweight="bold")
        ax_loadings.set_ylabel("Contribution to PC2", fontsize=14, fontweight="bold")
        ax_loadings.set_title(
            f"CV Window: {cv_window} - Variable Loadings",
            fontsize=14,
            fontweight="bold",
        )

    legend_elements_reef, legend_elements_habitat, legend_elements_arch = (
        create_legend_elements(
            color_map, shape_map, arch_styles, unique_reefs, unique_habitats
        )
    )

    leg1 = fig.legend(
        handles=legend_elements_reef,
        title="Reef",
        loc="center left",
        bbox_to_anchor=(1.0, 0.72),
        fontsize=10,
        title_fontsize=11,
    )
    leg2 = fig.legend(
        handles=legend_elements_habitat,
        title="Habitat",
        loc="center left",
        bbox_to_anchor=(1.0, 0.42),
        fontsize=10,
        title_fontsize=11,
    )
    leg3 = fig.legend(
        handles=legend_elements_arch,
        title="Arc",
        loc="center left",
        bbox_to_anchor=(1.0, 0.18),
        fontsize=10,
        title_fontsize=11,
    )

    plt.tight_layout(rect=[0, 0, 0.92, 0.97])
    plt.savefig(output_filename, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Variability figure saved: {output_filename}")


create_magnitude_figure(
    magnitude_results, os.path.join(output_dir, "PCA_Magnitude_publication_v2.png")
)

create_variability_figure(
    variability_results_dict,
    os.path.join(output_dir, "PCA_Variability_publication_v2.png"),
)

print("\n" + "=" * 60)
print("STEP 5: Saving Segment Summary")
print("=" * 60)
segment_cols = [col for col in df_processed.columns if "segment_count" in col]
summary_cols = ["Site_name", "HAB", "Arch"] + segment_cols
if segment_cols:
    df_segments_summary = df_processed[summary_cols]
    summary_path = os.path.join(output_dir, "resumo_segmentos_por_site_v2.xlsx")
    df_segments_summary.to_excel(summary_path, index=False)
    print(f"Segment summary saved: {summary_path}")
else:
    print("No segment count columns found.")

print("\n" + "=" * 60)
print("PROCESS COMPLETE")
print("=" * 60)
print(f"Output directory: {output_dir}")
print("Generated files:")
print("  - PCA_Magnitude_publication_v2.png")
print("  - PCA_Variability_publication_v2.png")
print("  - dados_consolidados_com_scores_PCAs_v2.xlsx")
print("  - loadings_PCA_*.csv")
print("=" * 60)
