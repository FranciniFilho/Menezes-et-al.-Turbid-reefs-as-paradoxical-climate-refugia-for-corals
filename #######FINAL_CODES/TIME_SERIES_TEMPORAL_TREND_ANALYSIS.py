### TIME_SERIES_TEMPORAL_TREND_ANALYSIS.py ###
#
# PURPOSE: Test whether environmental variability and extreme event
#          frequency are increasing over time (2002-2025)
#
# FOCUS: DLI (primary), with SST and Chl-a for context
#
# OUTPUTS:
#   1. Annual_Metrics_by_Site.csv
#   2. Annual_Metrics_by_Arc.csv
#   3. Trend_Test_Results_MannKendall.csv
#   4. Fig_1_Temporal_Trends_Panel.png (study period highlighted)
#   5. Fig_2_DLI_Focus_Trends.png (enlarged DLI panels)
#   6. Fig_3_Significance_Heatmap.png
#
# =====================================================================

import os

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
import glob
import re
import logging
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy import stats
from datetime import datetime
import warnings

warnings.filterwarnings("ignore")

# ============================================
# CONFIGURATION
# ============================================

# --- Temporal Configuration ---
START_YEAR = 2002
END_YEAR = 2025
STUDY_PERIOD = (2006, 2008)  # Highlighted region in plots

# --- Data Directories ---
SST_DIR = r"H:\remote sensing\CRW_SST_FULL"
MODIS_DIR = r"H:\remote sensing\MODIS_DATA_FULL"

# --- Sites CSV ---
SITES_CSV_PATH = r"C:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######22.04.23\DATA\sites_list_full.csv"

# --- Output Directory ---
# Keep outputs inside this repository (consistent with other pipelines).
OUTPUT_DIR = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "#######FINAL_RESULTS",
        "Temporal_Trend_Analysis",
    )
)
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- File Patterns ---
SST_PATTERN = "coraltemp_v3.1_*.nc"
KD490_PATTERN = "AQUA_MODIS.*.L3m.DAY.KD.Kd_490.4km.nc"
PAR_PATTERN = "AQUA_MODIS.*.L3m.DAY.PAR.par.4km.nc"
CHL_PATTERN = "AQUA_MODIS.*.L3m.DAY.CHL.chlor_a.4km.nc"

# --- DLI Calculation Constants (Gattuso et al. 2006) ---
KDPAR_A, KDPAR_B, KDPAR_C = 0.0665, 0.874, 0.00121

# --- Color Palette (MATCHES VARIABILITY_ANALYSIS_GEMINI.py) ---
COLOR_INNER = "#CD853F"  # Light brown (Peru)
COLOR_OUTER = "#4169E1"  # Royal Blue
COLOR_SST = "#D55E00"  # Vermelho-alaranjado
COLOR_DLI = "#0072B2"  # Azul escuro
COLOR_CHL = "#009E73"  # Verde-teal
STUDY_SHADE_COLOR = "#FFD700"  # Gold for study period

# --- Publication Style Settings (MATCHES VARIABILITY_ANALYSIS_GEMINI.py) ---
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 11,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

# ============================================
# UTILITY FUNCTIONS
# ============================================


def parse_date_from_filename(filename):
    """Extract date from filename (YYYYMMDD format)."""
    basename = os.path.basename(filename)
    match = re.search(r"(\d{8})", basename)
    if match:
        try:
            return pd.to_datetime(match.group(1), format="%Y%m%d")
        except ValueError:
            return None
    return None


def filter_files_by_period(files, start_year, end_year):
    """Filter files by year range."""
    filtered = []
    for f in files:
        date = parse_date_from_filename(f)
        if date is not None and start_year <= date.year <= end_year:
            filtered.append(f)
    return sorted(filtered)


# ============================================
# DATA LOADING (Memory-Efficient Year-by-Year)
# ============================================


def load_timeseries_by_year(
    files, var_names, sites_df, var_label, start_year, end_year
):
    """
    Memory-efficient loader: processes one year at a time.
    Returns dict: {site: DataFrame with 'date', 'value', 'year' columns}
    """
    print(f"\n  Loading {var_label} ({start_year}-{end_year})...")

    # Group files by year first
    files_by_year = {}
    for f in files:
        date = parse_date_from_filename(f)
        if date and start_year <= date.year <= end_year:
            year = date.year
            if year not in files_by_year:
                files_by_year[year] = []
            files_by_year[year].append((f, date))

    # Initialize storage
    site_data = {row["Site_HAB"]: [] for _, row in sites_df.iterrows()}

    # Process year by year
    for year in sorted(files_by_year.keys()):
        print(
            f"    Processing {var_label} year {year} ({len(files_by_year[year])} files)..."
        )

        for f, date in files_by_year[year]:
            try:
                with xr.open_dataset(f) as ds:
                    var_name = next((v for v in var_names if v in ds.data_vars), None)
                    if var_name is None:
                        continue

                    data = ds[var_name]
                    if "time" in data.dims:
                        data = data.isel(time=0)

                    for _, row in sites_df.iterrows():
                        site_key = row["Site_HAB"]
                        try:
                            value = float(
                                data.sel(
                                    lat=row["Latitude"],
                                    lon=row["Longitude"],
                                    method="nearest",
                                ).values
                            )
                            if not np.isnan(value):
                                site_data[site_key].append(
                                    {"date": date, "value": value, "year": date.year}
                                )
                        except:
                            continue
            except:
                continue

    # Convert to DataFrames
    result = {}
    for site_name, records in site_data.items():
        if records:
            result[site_name] = (
                pd.DataFrame(records).sort_values("date").reset_index(drop=True)
            )

    print(f"    [OK] Loaded data for {len(result)} sites")
    return result


def calculate_dli_by_year(kd490_files, par_files, sites_df, start_year, end_year):
    """
    Calculate DLI time series for each site (memory-efficient).
    """
    print(f"\n  Calculating DLI ({start_year}-{end_year})...")

    # Create date-indexed dictionaries
    par_by_date = {}
    for f in par_files:
        date = parse_date_from_filename(f)
        if date and start_year <= date.year <= end_year:
            par_by_date[date] = f

    # Group kd490 files by year
    kd_by_year = {}
    for f in kd490_files:
        date = parse_date_from_filename(f)
        if date and start_year <= date.year <= end_year and date in par_by_date:
            year = date.year
            if year not in kd_by_year:
                kd_by_year[year] = []
            kd_by_year[year].append((f, date))

    # Initialize storage
    site_data = {row["Site_HAB"]: [] for _, row in sites_df.iterrows()}

    # Process year by year
    for year in sorted(kd_by_year.keys()):
        print(f"    Processing DLI year {year} ({len(kd_by_year[year])} file pairs)...")

        for f_kd, date in kd_by_year[year]:
            f_par = par_by_date[date]

            try:
                with xr.open_dataset(f_kd) as ds_kd, xr.open_dataset(f_par) as ds_par:
                    kd490_data = ds_kd["Kd_490"]
                    par_data = ds_par["par"]

                    if "time" in kd490_data.dims:
                        kd490_data = kd490_data.isel(time=0)
                    if "time" in par_data.dims:
                        par_data = par_data.isel(time=0)

                    for _, row in sites_df.iterrows():
                        site_key = row["Site_HAB"]
                        lat, lon, depth = (
                            row["Latitude"],
                            row["Longitude"],
                            row["Depth_m"],
                        )

                        try:
                            kd490_val = float(
                                kd490_data.sel(
                                    lat=lat, lon=lon, method="nearest"
                                ).values
                            )
                            par_val = float(
                                par_data.sel(lat=lat, lon=lon, method="nearest").values
                            )

                            if (
                                np.isnan(kd490_val)
                                or np.isnan(par_val)
                                or kd490_val <= 0
                            ):
                                continue

                            # Gattuso equation for Kd_PAR
                            kdpar = (
                                KDPAR_A
                                + KDPAR_B * kd490_val
                                - KDPAR_C / (kd490_val + 1e-9)
                            )
                            dli = par_val * np.exp(-kdpar * depth)

                            if not np.isnan(dli) and dli > 0:
                                site_data[site_key].append(
                                    {"date": date, "value": dli, "year": date.year}
                                )
                        except:
                            continue
            except:
                continue

    # Convert to DataFrames
    result = {}
    for site_name, records in site_data.items():
        if records:
            result[site_name] = (
                pd.DataFrame(records).sort_values("date").reset_index(drop=True)
            )

    print(f"    [OK] Calculated DLI for {len(result)} sites")
    return result


# ============================================
# ANNUAL METRICS COMPUTATION
# ============================================


def compute_annual_metrics(site_data_dict, sites_df, global_percentile=90):
    """
    Compute variability metrics per site per year.
    Global percentile computed from FULL time series (2002-2025).
    """
    print("\n  Computing annual metrics...")
    results = []

    for site_key, df in site_data_dict.items():
        if len(df) == 0:
            continue

        # Get metadata for this Site_HAB
        site_meta = sites_df[sites_df["Site_HAB"] == site_key]
        if len(site_meta) == 0:
            continue
        site_meta = site_meta.iloc[0]
        site_name = site_meta["Site_name"]
        hab = site_meta["HAB"]
        arc = site_meta["Arc"]

        # Compute GLOBAL thresholds from full series
        all_values = df["value"].dropna().values
        global_p90 = np.percentile(all_values, global_percentile)
        global_p95 = np.percentile(all_values, 95)

        # Compute metrics per year
        for year, year_df in df.groupby("year"):
            vals = year_df["value"]
            n = len(vals)

            if n < 10:  # Minimum observations threshold
                continue

            mean_val = vals.mean()
            std_val = vals.std()

            results.append(
                {
                    "Site_HAB": site_key,
                    "Site_name": site_name,
                    "HAB": hab,
                    "Arc": arc,
                    "year": int(year),
                    "n_obs": n,
                    "mean": mean_val,
                    "std": std_val,
                    "var": vals.var(),
                    "cv": (std_val / mean_val * 100) if mean_val != 0 else np.nan,
                    "range": vals.max() - vals.min(),
                    "iqr": vals.quantile(0.75) - vals.quantile(0.25),
                    "min": vals.min(),
                    "max": vals.max(),
                    "extreme_p90_count": (vals > global_p90).sum(),
                    "extreme_p90_freq": (vals > global_p90).sum() / n * 100,
                    "extreme_p95_count": (vals > global_p95).sum(),
                    "extreme_p95_freq": (vals > global_p95).sum() / n * 100,
                    "global_p90_threshold": global_p90,
                    "global_p95_threshold": global_p95,
                }
            )

    return pd.DataFrame(results)


def aggregate_metrics_by_arc(annual_metrics_df):
    """
    Aggregate annual metrics across sites within each Arc.
    Returns mean values per Arc per year.
    """
    print("\n  Aggregating metrics by Arc...")

    # Columns to aggregate (mean)
    agg_cols = [
        "n_obs",
        "mean",
        "std",
        "var",
        "cv",
        "range",
        "iqr",
        "extreme_p90_count",
        "extreme_p90_freq",
        "extreme_p95_count",
        "extreme_p95_freq",
    ]

    # Group by Arc and year
    arc_metrics = (
        annual_metrics_df.groupby(["Arc", "year"])[agg_cols].mean().reset_index()
    )

    # Also compute site count per year
    site_counts = (
        annual_metrics_df.groupby(["Arc", "year"]).size().reset_index(name="n_sites")
    )
    arc_metrics = arc_metrics.merge(site_counts, on=["Arc", "year"])

    return arc_metrics


# ============================================
# TREND ANALYSIS (Mann-Kendall + Sen's Slope)
# ============================================


def mann_kendall_test(data):
    """
    Perform Mann-Kendall trend test with Sen's slope estimator.

    Args:
        data: 1D array of values (ordered by time)

    Returns:
        dict with 'trend', 'p_value', 'sens_slope', 'significance'
    """
    n = len(data)
    if n < 4:
        return {
            "trend": "insufficient_data",
            "p_value": np.nan,
            "z_score": np.nan,
            "sens_slope": np.nan,
            "significance": "",
        }

    # Remove NaN values
    data = np.array([x for x in data if not np.isnan(x)])
    n = len(data)

    if n < 4:
        return {
            "trend": "insufficient_data",
            "p_value": np.nan,
            "z_score": np.nan,
            "sens_slope": np.nan,
            "significance": "",
        }

    # Mann-Kendall S statistic
    s = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            s += np.sign(data[j] - data[i])

    # Variance of S
    var_s = n * (n - 1) * (2 * n + 5) / 18

    # Z-score
    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)
    else:
        z = 0

    # Two-tailed p-value
    p_value = 2 * (1 - stats.norm.cdf(abs(z)))

    # Sen's slope estimator
    slopes = []
    for i in range(n - 1):
        for j in range(i + 1, n):
            slopes.append((data[j] - data[i]) / (j - i))
    sens_slope = np.median(slopes)

    # Significance stars
    if p_value < 0.001:
        sig = "***"
    elif p_value < 0.01:
        sig = "**"
    elif p_value < 0.05:
        sig = "*"
    else:
        sig = "ns"

    trend = "increasing" if sens_slope > 0 else "decreasing"

    return {
        "trend": trend,
        "p_value": p_value,
        "z_score": z,
        "sens_slope": sens_slope,
        "significance": sig,
    }


def run_all_trend_tests(arc_metrics_df, variable_name):
    """
    Run Mann-Kendall tests for all Arc × Metric combinations.
    """
    print(f"\n  Running Mann-Kendall trend tests for {variable_name}...")

    metrics_to_test = ["cv", "std", "range", "extreme_p90_freq", "extreme_p95_freq"]

    results = []

    for arc in ["inner", "outer"]:
        arc_data = arc_metrics_df[arc_metrics_df["Arc"] == arc].sort_values("year")
        years = arc_data["year"].values

        for metric in metrics_to_test:
            if metric not in arc_data.columns:
                continue

            values = arc_data[metric].values

            mk_result = mann_kendall_test(values)

            results.append(
                {
                    "Variable": variable_name,
                    "Arc": arc,
                    "Metric": metric,
                    "n_years": len(values),
                    "year_start": int(years[0]) if len(years) > 0 else np.nan,
                    "year_end": int(years[-1]) if len(years) > 0 else np.nan,
                    "mean_value": np.nanmean(values),
                    **mk_result,
                }
            )

    return pd.DataFrame(results)


# ============================================
# VISUALIZATION FUNCTIONS
# ============================================


def create_trend_panel_figure(arc_metrics_dict, output_path, study_period=(2006, 2008)):
    """
    Create publication-quality trend figure with study period highlighted.

    Layout: 3 rows (SST, DLI, Chl-a) × 4 cols (CV, Std, Extreme P90, Range)
    """
    print("\n  Creating trend panel figure...")

    variables = ["SST", "DLI", "Chl-a"]
    var_colors = {"SST": COLOR_SST, "DLI": COLOR_DLI, "Chl-a": COLOR_CHL}
    metrics = ["cv", "std", "extreme_p90_freq", "range"]
    metric_labels = ["CV (%)", "Std Dev", "Extreme Freq (%)", "Range"]

    fig, axes = plt.subplots(3, 4, figsize=(18, 12))
    fig.suptitle(
        f"Temporal Trends in Environmental Variability ({START_YEAR}-{END_YEAR})\n"
        f"Study Period Highlighted: {study_period[0]}-{study_period[1]}",
        fontsize=14,
        fontweight="bold",
    )

    for row, var in enumerate(variables):
        if var not in arc_metrics_dict:
            for col in range(4):
                axes[row, col].text(
                    0.5,
                    0.5,
                    f"{var} data unavailable",
                    ha="center",
                    va="center",
                    transform=axes[row, col].transAxes,
                )
            continue

        arc_metrics = arc_metrics_dict[var]

        for col, (metric, ylabel) in enumerate(zip(metrics, metric_labels)):
            ax = axes[row, col]

            # Add study period shading
            ax.axvspan(
                study_period[0] - 0.5,
                study_period[1] + 0.5,
                alpha=0.25,
                color=STUDY_SHADE_COLOR,
                label="Study Period" if row == 0 and col == 0 else "",
            )

            for arc in ["inner", "outer"]:
                arc_color = COLOR_INNER if arc == "inner" else COLOR_OUTER
                arc_data = arc_metrics[arc_metrics["Arc"] == arc].sort_values("year")

                if len(arc_data) == 0 or metric not in arc_data.columns:
                    continue

                years = arc_data["year"].values
                values = arc_data[metric].values

                # Plot annual values
                ax.plot(
                    years,
                    values,
                    "o-",
                    color=arc_color,
                    lw=2,
                    ms=7,
                    label=f"{arc.capitalize()} Arc" if row == 0 and col == 3 else "",
                    markeredgecolor="white",
                    markeredgewidth=0.5,
                )

                # Linear trend line
                if len(years) >= 4:
                    valid_mask = ~np.isnan(values)
                    if valid_mask.sum() >= 4:
                        z = np.polyfit(years[valid_mask], values[valid_mask], 1)
                        p_fit = np.poly1d(z)
                        trend_years = np.arange(years.min(), years.max() + 1)
                        ax.plot(
                            trend_years,
                            p_fit(trend_years),
                            "--",
                            color=arc_color,
                            alpha=0.6,
                            lw=1.5,
                        )

                        # Mann-Kendall annotation
                        mk = mann_kendall_test(values)
                        if mk["significance"] and mk["significance"] != "ns":
                            ax.text(
                                0.95,
                                0.95,
                                f"{mk['significance']}\nslope={mk['sens_slope']:.4f}/yr",
                                transform=ax.transAxes,
                                ha="right",
                                va="top",
                                fontsize=9,
                                bbox=dict(
                                    boxstyle="round",
                                    facecolor="white",
                                    alpha=0.9,
                                    edgecolor=arc_color,
                                    linewidth=1.5,
                                ),
                            )

            # Styling
            ax.set_xlim(START_YEAR - 0.5, END_YEAR + 0.5)
            ax.grid(True, alpha=0.3, linestyle="--")

            if row == 0:
                ax.set_title(ylabel, fontsize=12, fontweight="bold")
            if col == 0:
                ax.set_ylabel(
                    var, fontsize=12, fontweight="bold", color=var_colors[var]
                )
            if row == 2:
                ax.set_xlabel("Year", fontsize=10)

            ax.tick_params(axis="both", labelsize=9)

    # Legend
    handles = [
        Line2D([0], [0], color=COLOR_INNER, lw=2, marker="o", ms=7, label="Inner Arc"),
        Line2D([0], [0], color=COLOR_OUTER, lw=2, marker="o", ms=7, label="Outer Arc"),
        mpatches.Patch(
            color=STUDY_SHADE_COLOR,
            alpha=0.25,
            label=f"Study Period ({study_period[0]}-{study_period[1]})",
        ),
    ]
    fig.legend(
        handles=handles, loc="upper right", bbox_to_anchor=(0.98, 0.98), fontsize=11
    )

    plt.tight_layout(rect=(0, 0, 0.92, 0.95))
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"    [OK] Saved: {os.path.basename(output_path)}")


def create_dli_focus_figure(arc_metrics_dict, output_path, study_period=(2006, 2008)):
    """
    Create enlarged DLI-focused figure (2x2 panel).
    """
    print("\n  Creating DLI focus figure...")

    if "DLI" not in arc_metrics_dict:
        print("    WARNING: DLI data not available")
        return

    arc_metrics = arc_metrics_dict["DLI"]
    metrics = ["cv", "std", "extreme_p90_freq", "range"]
    metric_labels = [
        "CV (%)",
        "Std Dev (mol m-2 d-1)",
        "Extreme Freq (%)",
        "Range (mol m-2 d-1)",
    ]

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()

    fig.suptitle(
        f"DLI Temporal Variability Trends ({START_YEAR}-{END_YEAR})\n"
        f"Primary Focus Variable - Study Period: {study_period[0]}-{study_period[1]}",
        fontsize=14,
        fontweight="bold",
        color=COLOR_DLI,
    )

    for i, (metric, ylabel) in enumerate(zip(metrics, metric_labels)):
        ax = axes[i]

        # Study period shading
        ax.axvspan(
            study_period[0] - 0.5,
            study_period[1] + 0.5,
            alpha=0.25,
            color=STUDY_SHADE_COLOR,
        )

        for arc in ["inner", "outer"]:
            arc_color = COLOR_INNER if arc == "inner" else COLOR_OUTER
            arc_data = arc_metrics[arc_metrics["Arc"] == arc].sort_values("year")

            if len(arc_data) == 0 or metric not in arc_data.columns:
                continue

            years = arc_data["year"].values
            values = arc_data[metric].values

            ax.plot(
                years,
                values,
                "o-",
                color=arc_color,
                lw=2.5,
                ms=8,
                label=f"{arc.capitalize()} Arc",
                markeredgecolor="white",
                markeredgewidth=0.5,
            )

            # Trend line and test
            if len(years) >= 4:
                valid_mask = ~np.isnan(values)
                if valid_mask.sum() >= 4:
                    z = np.polyfit(years[valid_mask], values[valid_mask], 1)
                    p_fit = np.poly1d(z)
                    ax.plot(
                        np.arange(years.min(), years.max() + 1),
                        p_fit(np.arange(years.min(), years.max() + 1)),
                        "--",
                        color=arc_color,
                        alpha=0.6,
                        lw=2,
                    )

                    mk = mann_kendall_test(values)
                    ax.text(
                        0.95,
                        0.15 if arc == "outer" else 0.85,
                        f"{arc.capitalize()}: {mk['significance']}\nslope={mk['sens_slope']:.4f}/yr",
                        transform=ax.transAxes,
                        ha="right",
                        va="bottom" if arc == "outer" else "top",
                        fontsize=10,
                        color=arc_color,
                        fontweight="bold",
                        bbox=dict(boxstyle="round", facecolor="white", alpha=0.9),
                    )

        ax.set_xlim(START_YEAR - 0.5, END_YEAR + 0.5)
        ax.grid(True, alpha=0.3, linestyle="--")
        ax.set_title(metric_labels[i], fontsize=12, fontweight="bold")
        ax.set_xlabel("Year", fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.tick_params(axis="both", labelsize=9)

    # Legend
    handles = [
        Line2D(
            [0], [0], color=COLOR_INNER, lw=2.5, marker="o", ms=8, label="Inner Arc"
        ),
        Line2D(
            [0], [0], color=COLOR_OUTER, lw=2.5, marker="o", ms=8, label="Outer Arc"
        ),
        mpatches.Patch(color=STUDY_SHADE_COLOR, alpha=0.25, label="Study Period"),
    ]
    fig.legend(
        handles=handles, loc="upper right", bbox_to_anchor=(0.98, 0.98), fontsize=11
    )

    plt.tight_layout(rect=(0, 0, 0.92, 0.93))
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"    [OK] Saved: {os.path.basename(output_path)}")


def create_significance_heatmap(trend_results, output_path):
    """
    Create heatmap showing trend significance for all combinations.
    """
    print("\n  Creating significance heatmap...")

    metrics = ["cv", "std", "range", "extreme_p90_freq"]
    metric_labels = ["CV", "Std Dev", "Range", "Extreme Freq"]
    variables = ["SST", "DLI", "Chl-a"]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(
        "Trend Significance Heatmap (Mann-Kendall Test, Sen's Slope)",
        fontsize=14,
        fontweight="bold",
    )

    for i, arc in enumerate(["inner", "outer"]):
        ax = axes[i]
        arc_data = trend_results[trend_results["Arc"] == arc]

        # Create matrix
        matrix = np.zeros((len(variables), len(metrics)))
        annot = np.empty((len(variables), len(metrics)), dtype=object)

        for j, var in enumerate(variables):
            for k, metric in enumerate(metrics):
                row = arc_data[
                    (arc_data["Variable"] == var) & (arc_data["Metric"] == metric)
                ]
                if len(row) > 0:
                    slope = row["sens_slope"].values[0]
                    sig = row["significance"].values[0]
                    matrix[j, k] = slope
                    annot[j, k] = f"{slope:.4f}\n{sig}" if sig else f"{slope:.4f}"
                else:
                    matrix[j, k] = np.nan
                    annot[j, k] = "N/A"

        # Plot heatmap
        vmax = np.nanmax(np.abs(matrix)) if np.any(~np.isnan(matrix)) else 1
        cmap = plt.cm.get_cmap("RdBu_r")
        im = ax.imshow(matrix, cmap=cmap, aspect="auto", vmin=-vmax, vmax=vmax)

        # Add annotations
        for j in range(len(variables)):
            for k in range(len(metrics)):
                if not np.isnan(matrix[j, k]):
                    text_color = "white" if abs(matrix[j, k]) > vmax * 0.5 else "black"
                    ax.text(
                        k,
                        j,
                        annot[j, k],
                        ha="center",
                        va="center",
                        fontsize=9,
                        color=text_color,
                    )

        ax.set_xticks(range(len(metrics)))
        ax.set_xticklabels(metric_labels, fontsize=10)
        ax.set_yticks(range(len(variables)))
        ax.set_yticklabels(variables, fontsize=11)
        ax.set_title(f"{arc.capitalize()} Arc", fontsize=12, fontweight="bold")

        plt.colorbar(im, ax=ax, label="Sen's Slope (per year)")

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"    [OK] Saved: {os.path.basename(output_path)}")


def create_extreme_events_timeline(
    arc_metrics_dict, output_path, study_period=(2006, 2008)
):
    """
    Create figure showing extreme events frequency over time for all variables.
    """
    print("\n  Creating extreme events timeline figure...")

    variables = ["SST", "DLI", "Chl-a"]
    var_colors = {"SST": COLOR_SST, "DLI": COLOR_DLI, "Chl-a": COLOR_CHL}

    fig, axes = plt.subplots(3, 1, figsize=(14, 10), sharex=True)
    fig.suptitle(
        f"Extreme Events Frequency Over Time ({START_YEAR}-{END_YEAR})\n"
        f"Events > 90th Percentile (Global Threshold)",
        fontsize=14,
        fontweight="bold",
    )

    for i, var in enumerate(variables):
        ax = axes[i]

        if var not in arc_metrics_dict:
            ax.text(
                0.5,
                0.5,
                f"{var} data unavailable",
                ha="center",
                va="center",
                transform=ax.transAxes,
            )
            continue

        arc_metrics = arc_metrics_dict[var]

        # Study period shading
        ax.axvspan(
            study_period[0] - 0.5,
            study_period[1] + 0.5,
            alpha=0.25,
            color=STUDY_SHADE_COLOR,
        )

        for arc in ["inner", "outer"]:
            arc_color = COLOR_INNER if arc == "inner" else COLOR_OUTER
            arc_data = arc_metrics[arc_metrics["Arc"] == arc].sort_values("year")

            if len(arc_data) == 0 or "extreme_p90_freq" not in arc_data.columns:
                continue

            years = arc_data["year"].values
            values = arc_data["extreme_p90_freq"].values

            ax.bar(
                [y + (0.2 if arc == "outer" else -0.2) for y in years],
                values,
                width=0.35,
                color=arc_color,
                alpha=0.7,
                label=f"{arc.capitalize()} Arc",
                edgecolor="white",
            )

        ax.set_xlim(START_YEAR - 0.5, END_YEAR + 0.5)
        ax.grid(True, alpha=0.3, linestyle="--", axis="y")
        ax.set_ylabel(f"{var}\nExtreme Freq (%)", fontsize=11, color=var_colors[var])
        ax.tick_params(axis="y", labelsize=9)

        if i == 0:
            ax.legend(loc="upper right", fontsize=10)

    axes[-1].set_xlabel("Year", fontsize=11)
    axes[-1].tick_params(axis="x", labelsize=10)

    plt.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()
    print(f"    [OK] Saved: {os.path.basename(output_path)}")


# ============================================
# MAIN EXECUTION
# ============================================

if __name__ == "__main__":
    print("=" * 70)
    print("TEMPORAL TREND ANALYSIS: Environmental Variability 2002-2025")
    print("Testing: Is variability and extreme frequency increasing over time?")
    print("=" * 70)

    # --- STEP 1: Load Sites ---
    print("\n[STEP 1] Loading sites metadata...")
    sites_df = pd.read_csv(SITES_CSV_PATH, sep=";")
    sites_df.columns = [c.strip() for c in sites_df.columns]
    sites_df = sites_df.rename(columns={"Arch": "Arc"})
    sites_df["Arc"] = (
        sites_df["Arc"].astype(str).str.strip().str.lower().str.split("_").str[0]
    )
    sites_df = sites_df[sites_df["Arc"].isin(["inner", "outer"])].copy()
    sites_df["Site_name"] = [str(v).strip() for v in sites_df["Site_name"]]
    sites_df["HAB"] = [str(v).strip() for v in sites_df["HAB"]]
    sites_df["Site_HAB"] = [
        f"{site}_{hab}" for site, hab in zip(sites_df["Site_name"], sites_df["HAB"])
    ]
    print(f"  [OK] Loaded {len(sites_df)} sites")
    print(f"    Inner Arc: {(sites_df['Arc'] == 'inner').sum()} sites")
    print(f"    Outer Arc: {(sites_df['Arc'] == 'outer').sum()} sites")

    # --- STEP 2: Find and Filter Files ---
    print("\n[STEP 2] Finding satellite data files...")
    sst_files = filter_files_by_period(
        glob.glob(os.path.join(SST_DIR, SST_PATTERN)), START_YEAR, END_YEAR
    )
    kd490_files = filter_files_by_period(
        glob.glob(os.path.join(MODIS_DIR, KD490_PATTERN)), START_YEAR, END_YEAR
    )
    par_files = filter_files_by_period(
        glob.glob(os.path.join(MODIS_DIR, PAR_PATTERN)), START_YEAR, END_YEAR
    )
    chl_files = filter_files_by_period(
        glob.glob(os.path.join(MODIS_DIR, CHL_PATTERN)), START_YEAR, END_YEAR
    )

    print(f"  SST files: {len(sst_files)}")
    print(f"  Kd490 files: {len(kd490_files)}")
    print(f"  PAR files: {len(par_files)}")
    print(f"  Chl-a files: {len(chl_files)}")

    # --- STEP 3: Load Time Series ---
    print("\n[STEP 3] Loading time series data (year-by-year)...")
    sst_data = load_timeseries_by_year(
        sst_files, ["analysed_sst", "sst"], sites_df, "SST", START_YEAR, END_YEAR
    )
    chl_data = load_timeseries_by_year(
        chl_files, ["chlor_a"], sites_df, "Chl-a", START_YEAR, END_YEAR
    )
    dli_data = calculate_dli_by_year(
        kd490_files, par_files, sites_df, START_YEAR, END_YEAR
    )

    # --- STEP 4: Compute Annual Metrics ---
    print("\n[STEP 4] Computing annual variability metrics...")
    sst_annual = compute_annual_metrics(sst_data, sites_df)
    sst_annual["Variable"] = "SST"

    dli_annual = compute_annual_metrics(dli_data, sites_df)
    dli_annual["Variable"] = "DLI"

    chl_annual = compute_annual_metrics(chl_data, sites_df)
    chl_annual["Variable"] = "Chl-a"

    # Combine all
    all_annual = pd.concat([sst_annual, dli_annual, chl_annual], ignore_index=True)
    all_annual.to_csv(
        os.path.join(OUTPUT_DIR, "Annual_Metrics_by_Site.csv"), index=False
    )
    print(f"  [OK] Saved: Annual_Metrics_by_Site.csv ({len(all_annual)} records)")

    # --- STEP 5: Aggregate by Arc ---
    print("\n[STEP 5] Aggregating metrics by Arc...")
    sst_arc = aggregate_metrics_by_arc(sst_annual)
    sst_arc["Variable"] = "SST"

    dli_arc = aggregate_metrics_by_arc(dli_annual)
    dli_arc["Variable"] = "DLI"

    chl_arc = aggregate_metrics_by_arc(chl_annual)
    chl_arc["Variable"] = "Chl-a"

    all_arc = pd.concat([sst_arc, dli_arc, chl_arc], ignore_index=True)
    all_arc.to_csv(os.path.join(OUTPUT_DIR, "Annual_Metrics_by_Arc.csv"), index=False)
    print(f"  [OK] Saved: Annual_Metrics_by_Arc.csv ({len(all_arc)} records)")

    # --- STEP 6: Trend Tests ---
    print("\n[STEP 6] Running Mann-Kendall trend tests...")

    # Run tests for each variable
    all_trends = pd.concat(
        [
            run_all_trend_tests(sst_arc, "SST"),
            run_all_trend_tests(dli_arc, "DLI"),
            run_all_trend_tests(chl_arc, "Chl-a"),
        ],
        ignore_index=True,
    )

    all_trends.to_csv(os.path.join(OUTPUT_DIR, "Trend_Test_Results.csv"), index=False)
    print(f"  [OK] Saved: Trend_Test_Results.csv ({len(all_trends)} tests)")

    # --- STEP 7: Visualizations ---
    print("\n[STEP 7] Creating publication figures...")

    arc_metrics_dict = {"SST": sst_arc, "DLI": dli_arc, "Chl-a": chl_arc}

    create_trend_panel_figure(
        arc_metrics_dict,
        os.path.join(OUTPUT_DIR, "Fig_1_Temporal_Trends_Panel.png"),
        STUDY_PERIOD,
    )

    create_dli_focus_figure(
        arc_metrics_dict,
        os.path.join(OUTPUT_DIR, "Fig_2_DLI_Focus_Trends.png"),
        STUDY_PERIOD,
    )

    create_significance_heatmap(
        all_trends, os.path.join(OUTPUT_DIR, "Fig_3_Significance_Heatmap.png")
    )

    create_extreme_events_timeline(
        arc_metrics_dict,
        os.path.join(OUTPUT_DIR, "Fig_4_Extreme_Events_Timeline.png"),
        STUDY_PERIOD,
    )

    # --- Summary Report ---
    print("\n" + "=" * 70)
    print("SUMMARY: Significant Trends Detected")
    print("=" * 70)

    sig_trends = all_trends[all_trends["significance"].isin(["*", "**", "***"])]
    if len(sig_trends) > 0:
        print(f"\n{len(sig_trends)} significant trends found:")
        for _, row in sig_trends.iterrows():
            print(
                f"  {row['Variable']} | {str(row['Arc']).capitalize()} Arc | {row['Metric']}: "
                f"{row['trend']} (p={row['p_value']:.4f}, slope={row['sens_slope']:.4f})"
            )
    else:
        print("\nNo significant trends detected at alpha=0.05")

    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE!")
    print(f"All outputs saved to: {OUTPUT_DIR}")
    print("=" * 70)
