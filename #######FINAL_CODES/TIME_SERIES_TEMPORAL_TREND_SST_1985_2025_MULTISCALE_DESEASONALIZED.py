### TIME_SERIES_TEMPORAL_TREND_SST_1985_2025_MULTISCALE_DESEASONALIZED.py ###
#
# PURPOSE
#   Temporal trend analysis focused on SST only (CoralTemp; 1985-2025).
#   Quantifies whether thermal variability and extreme-event frequency increase
#   over time, comparing Inner vs Outer arc sites, with study period highlighted.
#
# DESIGN PRINCIPLES
#   - Efficient extraction: read only SST at site pixels from each daily NetCDF.
#   - Resume-from-checkpoint: caches per-year site tables to avoid re-reading files.
#   - Deseasonalized metrics: anomalies from a fixed baseline DOY climatology.
#   - Multi-scale variability: mean rolling STD at 2d and 30d on anomalies.
#   - Extremes: exceedance frequency and clustered event counts above P90/P95
#     thresholds computed from baseline anomalies.
#
# OUTPUTS
#   Written to: #######FINAL_RESULTS/Temporal_Trend_Analysis_SST_1985_2025
#     - Annual_SST_Metrics_by_Site.csv
#     - Annual_SST_Metrics_by_Arc.csv
#     - Trend_Test_Results_SST.csv
#     - Fig_1_SST_Mean_Trend.png
#     - Fig_2_SST_Variability_Multiscale.png
#     - Fig_3_SST_Extremes_Frequency_Events.png
#     - Fig_S1_Data_Availability_SST.png
#

import os

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import glob
import pickle
import re
import warnings
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy import stats


warnings.filterwarnings("ignore")


# =========================================================
# Configuration
# =========================================================

START_YEAR = 1985
END_YEAR = 2025

# Highlight study period (visual only)
STUDY_PERIOD = (2006, 2008)

# Fixed baseline for seasonal climatology and extremes thresholds
CLIM_BASELINE_START = 1985
CLIM_BASELINE_END = 2012

# Threshold percentiles for extremes (on anomalies)
EXTREME_P90 = 90
EXTREME_P95 = 95

# Event definition
EVENT_MIN_GAP_DAYS = 2  # new event if gap > this

# Cache behavior
USE_CACHE = True
FORCE_REBUILD_CACHE = False

# Data paths
SST_DIR = r"H:\remote sensing\CRW_SST_FULL"
SST_PATTERN = "coraltemp_v3.1_*.nc"

SITES_CSV_PATH = r"C:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######22.04.23\DATA\sites_list_full.csv"

# Output directory inside repository
OUTPUT_DIR = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "#######FINAL_RESULTS",
        "Temporal_Trend_Analysis_SST_1985_2025",
    )
)
CACHE_DIR = os.path.join(OUTPUT_DIR, "cache")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(CACHE_DIR, exist_ok=True)

# Colors & style (aligned with TIME_SERIES_VARIABILITY_ANALYSIS_GEMINI.py)
COLOR_INNER = "#CD853F"  # Peru
COLOR_OUTER = "#4169E1"  # Royal Blue
COLOR_SST = "#D55E00"  # Vermelho-alaranjado
STUDY_SHADE_COLOR = "#FFD700"  # Gold

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


# =========================================================
# Helpers
# =========================================================


def get_netcdf_engine():
    try:
        import netCDF4  # noqa: F401

        return "netcdf4"
    except Exception:
        pass
    try:
        import h5netcdf  # noqa: F401

        return "h5netcdf"
    except Exception:
        pass
    return None


NETCDF_ENGINE = get_netcdf_engine()


def parse_date_from_filename(filename: str):
    m = re.search(r"(\d{8})", os.path.basename(filename))
    if not m:
        return None
    try:
        return pd.to_datetime(m.group(1), format="%Y%m%d")
    except Exception:
        return None


def filter_files_by_period(
    files: List[str], start_year: int, end_year: int
) -> List[str]:
    out = []
    for f in files:
        d = parse_date_from_filename(f)
        if d is not None and start_year <= d.year <= end_year:
            out.append(f)
    return sorted(out)


def load_sites(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=";")
    df.columns = [c.strip() for c in df.columns]
    if "Arch" in df.columns and "Arc" not in df.columns:
        df = df.rename(columns={"Arch": "Arc"})

    required = {"Site_name", "HAB", "Latitude", "Longitude", "Depth_m", "Arc"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in sites file: {missing}")

    df["Arc"] = df["Arc"].astype(str).str.strip().str.lower().str.split("_").str[0]
    df = df[df["Arc"].isin(["inner", "outer"])].copy()
    df["Site_name"] = [str(v).strip() for v in df["Site_name"]]
    df["HAB"] = [str(v).strip() for v in df["HAB"]]
    df["Site_HAB"] = [f"{site}_{hab}" for site, hab in zip(df["Site_name"], df["HAB"])]
    return df


def cache_load(path: str):
    if USE_CACHE and (not FORCE_REBUILD_CACHE) and os.path.exists(path):
        with open(path, "rb") as f:
            return pickle.load(f)
    return None


def cache_save(path: str, obj) -> None:
    if USE_CACHE:
        with open(path, "wb") as f:
            pickle.dump(obj, f)


def sanitize_sst(val: float) -> float:
    # CoralTemp is in degC; treat absurd values and common fill values as missing.
    if np.isnan(val):
        return np.nan
    if val < -5 or val > 45:
        return np.nan
    return float(val)


def daily_index_for_year(year: int) -> pd.DatetimeIndex:
    return pd.date_range(f"{year}-01-01", f"{year}-12-31", freq="D")


# =========================================================
# Efficient SST extraction for sites (per-year cache)
# =========================================================


def compute_site_indexers_from_reference_file(
    ref_file: str, sites_df: pd.DataFrame
) -> Tuple[List[str], xr.DataArray, xr.DataArray]:
    """Compute lat/lon integer indexers (paired) for each site."""
    open_kwargs = {}
    if NETCDF_ENGINE:
        open_kwargs["engine"] = NETCDF_ENGINE

    with xr.open_dataset(ref_file, **open_kwargs) as ds:
        lat = ds["lat"].values
        lon = ds["lon"].values

    site_names = sites_df["Site_HAB"].astype(str).tolist()
    ilat = []
    ilon = []
    for _, r in sites_df.iterrows():
        ilat.append(int(np.argmin(np.abs(lat - float(r["Latitude"])))))
        ilon.append(int(np.argmin(np.abs(lon - float(r["Longitude"])))))

    ilat_da = xr.DataArray(
        np.array(ilat, dtype=int), dims=("site",), coords={"site": site_names}
    )
    ilon_da = xr.DataArray(
        np.array(ilon, dtype=int), dims=("site",), coords={"site": site_names}
    )
    return site_names, ilat_da, ilon_da


def extract_sst_sites_for_year(
    year: int,
    year_files: List[str],
    site_names: List[str],
    ilat_da: xr.DataArray,
    ilon_da: xr.DataArray,
) -> pd.DataFrame:
    """Extract analysed_sst for all sites for all available days in a year."""
    open_kwargs = {}
    if NETCDF_ENGINE:
        open_kwargs["engine"] = NETCDF_ENGINE

    dates: List[pd.Timestamp] = []
    values: List[np.ndarray] = []

    for i, f in enumerate(year_files, 1):
        if i % 200 == 0:
            print(f"  SST {year}: {i}/{len(year_files)}")
        d = parse_date_from_filename(f)
        if d is None:
            continue
        try:
            with xr.open_dataset(f, **open_kwargs) as ds:
                da = ds["analysed_sst"]
                if "time" in da.dims:
                    da = da.isel(time=0)

                # Vectorized paired indexing using same 'site' dimension
                v = da.isel(lat=ilat_da, lon=ilon_da).values
                v = np.array([sanitize_sst(float(x)) for x in v], dtype=float)
                dates.append(d)
                values.append(v)
        except Exception:
            continue

    if not dates:
        return pd.DataFrame(columns=site_names)

    df = pd.DataFrame(
        np.vstack(values), index=pd.DatetimeIndex(dates), columns=site_names
    ).sort_index()

    # Normalize to a complete daily index for the year
    idx = daily_index_for_year(year)
    df = df[~df.index.duplicated(keep="first")].reindex(idx)
    return df


def load_or_build_sst_daily_table(
    sst_files: List[str], sites_df: pd.DataFrame
) -> pd.DataFrame:
    """Builds (or loads) a daily SST table (date x site) using per-year caches."""
    if len(sst_files) == 0:
        raise ValueError("No SST files found after filtering.")

    # Reference file for grid coords
    ref_file = sst_files[0]
    site_names, ilat_da, ilon_da = compute_site_indexers_from_reference_file(
        ref_file, sites_df
    )

    # Group files by year
    by_year: Dict[int, List[str]] = {}
    for f in sst_files:
        d = parse_date_from_filename(f)
        if d is None:
            continue
        by_year.setdefault(int(d.year), []).append(f)

    parts = []
    for year in range(START_YEAR, END_YEAR + 1):
        if year not in by_year:
            continue

        cache_parquet = os.path.join(
            CACHE_DIR, f"sst_sites_daily_site_hab_{year}.parquet"
        )
        if USE_CACHE and (not FORCE_REBUILD_CACHE) and os.path.exists(cache_parquet):
            df_y = pd.read_parquet(cache_parquet)
            # Ensure columns order
            df_y = df_y.reindex(columns=site_names)
            parts.append(df_y)
            continue

        print(f"\nExtracting SST for year {year} ({len(by_year[year])} files)...")
        df_y = extract_sst_sites_for_year(
            year, by_year[year], site_names, ilat_da, ilon_da
        )
        df_y.to_parquet(cache_parquet)
        parts.append(df_y)
        print(f"  [OK] Cached: {cache_parquet}")

    if not parts:
        raise ValueError("No per-year SST site tables could be built.")

    df = pd.concat(parts, axis=0).sort_index()
    df = df[~df.index.duplicated(keep="first")]
    return df


# =========================================================
# Deseasonalization and metrics
# =========================================================


def compute_doy_climatology(
    series: pd.Series, baseline_years: Tuple[int, int], smooth_window_days: int = 31
) -> pd.Series:
    """Compute smoothed DOY climatology (median) from baseline years."""
    s = series.dropna().sort_index()
    if len(s) < 365:
        return pd.Series(dtype=float)

    s_base = s[
        (s.index.year >= baseline_years[0]) & (s.index.year <= baseline_years[1])
    ]
    if len(s_base) < 365:
        # Fallback: use all available data
        s_base = s

    doy = s_base.index.dayofyear.to_numpy()
    doy = np.where(doy == 366, 365, doy)

    clim = pd.DataFrame({"doy": doy, "v": s_base.values}).groupby("doy")["v"].median()
    clim = clim.reindex(np.arange(1, 366)).interpolate(limit_direction="both")

    pad = smooth_window_days // 2
    arr = clim.values
    ext = np.r_[arr[-pad:], arr, arr[:pad]]
    sm = (
        pd.Series(ext)
        .rolling(smooth_window_days, center=True, min_periods=1)
        .median()
        .values
    )
    return pd.Series(sm[pad : pad + 365], index=np.arange(1, 366), dtype=float)


def compute_anomalies_from_climatology(
    series: pd.Series, clim_doy: pd.Series
) -> pd.Series:
    s = series.dropna().sort_index()
    if len(s) == 0 or len(clim_doy) != 365:
        return pd.Series(dtype=float)
    doy = np.where(s.index.dayofyear == 366, 365, s.index.dayofyear)
    clim_on_dates = pd.Series(doy, index=s.index).map(clim_doy).astype(float)
    return s - clim_on_dates.values


def rolling_std_mean(series: pd.Series, window: int) -> float:
    if len(series) < max(5, window):
        return np.nan
    rs = series.rolling(
        window=window, min_periods=max(2, int(window * 0.6)), center=True
    ).std()
    return float(rs.mean())


def exceedance_frequency_percent(series: pd.Series, threshold: float) -> float:
    s = series
    denom = float(s.notna().sum())
    if denom <= 0:
        return np.nan
    num = float(((s > threshold) & s.notna()).sum())
    return num / denom * 100.0


def count_events(series: pd.Series, threshold: float, min_gap_days: int) -> int:
    s = series.dropna().sort_index()
    above = s.index[s > threshold]
    if len(above) == 0:
        return 0
    events = 1
    for i in range(1, len(above)):
        if (above[i] - above[i - 1]).days > min_gap_days:
            events += 1
    return int(events)


def compute_annual_sst_metrics_by_site(
    sst_daily_df: pd.DataFrame, sites_df: pd.DataFrame
) -> pd.DataFrame:
    """Compute annual SST metrics per site (raw + anomaly-based)."""
    results = []

    for _, r in sites_df.iterrows():
        site = str(r["Site_HAB"])
        site_name = str(r["Site_name"])
        hab = str(r["HAB"])
        arc = str(r["Arc"])
        if site not in sst_daily_df.columns:
            continue

        s_raw = sst_daily_df[site].copy()
        s_raw.name = "sst"
        if s_raw.notna().sum() < 365:
            continue

        clim = compute_doy_climatology(s_raw, (CLIM_BASELINE_START, CLIM_BASELINE_END))
        anom = compute_anomalies_from_climatology(s_raw, clim)

        # thresholds computed on baseline anomalies
        anom_base = anom[
            (anom.index.year >= CLIM_BASELINE_START)
            & (anom.index.year <= CLIM_BASELINE_END)
        ]
        if anom_base.notna().sum() < 365:
            anom_base = anom.dropna()

        thr90 = float(np.nanpercentile(anom_base.values, EXTREME_P90))
        thr95 = float(np.nanpercentile(anom_base.values, EXTREME_P95))

        for year in range(START_YEAR, END_YEAR + 1):
            idx = daily_index_for_year(year)
            s_y = s_raw.reindex(idx)
            a_y = anom.reindex(idx)

            n_obs = int(a_y.notna().sum())
            if n_obs < 60:
                continue

            # Raw metrics (context warming)
            mean_sst = float(s_y.mean(skipna=True))
            max_sst = float(s_y.max(skipna=True))

            # Variability metrics (anomalies)
            msv_2d = rolling_std_mean(a_y, 2)
            msv_30d = rolling_std_mean(a_y, 30)
            std_full = float(a_y.std(skipna=True))
            iqr_full = float(a_y.quantile(0.75) - a_y.quantile(0.25))

            # Extremes (anomalies)
            dayfreq_p90 = exceedance_frequency_percent(a_y, thr90)
            dayfreq_p95 = exceedance_frequency_percent(a_y, thr95)
            events_p90 = count_events(a_y, thr90, EVENT_MIN_GAP_DAYS)
            events_p95 = count_events(a_y, thr95, EVENT_MIN_GAP_DAYS)

            results.append(
                {
                    "Site_HAB": site,
                    "Site_name": site_name,
                    "HAB": hab,
                    "Arc": arc,
                    "year": int(year),
                    "n_obs": n_obs,
                    "mean_sst": mean_sst,
                    "max_sst": max_sst,
                    "msv_2d": msv_2d,
                    "msv_30d": msv_30d,
                    "std_full": std_full,
                    "iqr_full": iqr_full,
                    "dayfreq_p90": dayfreq_p90,
                    "dayfreq_p95": dayfreq_p95,
                    "events_p90": events_p90,
                    "events_p95": events_p95,
                    "thr90_anom": thr90,
                    "thr95_anom": thr95,
                }
            )

    return pd.DataFrame(results)


def aggregate_by_arc_year(df_site: pd.DataFrame) -> pd.DataFrame:
    agg_cols = [
        "n_obs",
        "mean_sst",
        "max_sst",
        "msv_2d",
        "msv_30d",
        "std_full",
        "iqr_full",
        "dayfreq_p90",
        "dayfreq_p95",
        "events_p90",
        "events_p95",
    ]
    out = df_site.groupby(["Arc", "year"])[agg_cols].mean().reset_index()
    out["n_sites"] = df_site.groupby(["Arc", "year"]).size().values
    return out


# =========================================================
# Trend tests
# =========================================================


def mann_kendall_sens(y: np.ndarray) -> dict:
    y = np.array([v for v in y if not np.isnan(v)], dtype=float)
    n = len(y)
    if n < 4:
        return {
            "trend": "insufficient",
            "p_value": np.nan,
            "sens_slope": np.nan,
            "significance": "ns",
        }

    s = 0
    slopes = []
    for i in range(n - 1):
        for j in range(i + 1, n):
            s += np.sign(y[j] - y[i])
            slopes.append((y[j] - y[i]) / (j - i))

    var_s = n * (n - 1) * (2 * n + 5) / 18
    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)
    else:
        z = 0.0

    p = 2 * (1 - stats.norm.cdf(abs(z)))
    slope = float(np.median(slopes))
    trend = "increasing" if slope > 0 else "decreasing"

    sig = "ns"
    if p < 0.001:
        sig = "***"
    elif p < 0.01:
        sig = "**"
    elif p < 0.05:
        sig = "*"

    return {
        "trend": trend,
        "p_value": float(p),
        "sens_slope": slope,
        "significance": sig,
    }


def run_trend_tests(df_arc: pd.DataFrame) -> pd.DataFrame:
    metrics = ["mean_sst", "msv_2d", "msv_30d", "std_full", "dayfreq_p95", "events_p95"]
    rows = []
    for arc in ["inner", "outer"]:
        sub = df_arc[df_arc["Arc"] == arc].sort_values("year")
        for m in metrics:
            if m not in sub.columns:
                continue
            res = mann_kendall_sens(sub[m].values)
            rows.append({"Arc": arc, "Metric": m, "n_years": int(len(sub)), **res})
    return pd.DataFrame(rows)


# =========================================================
# Figures
# =========================================================


def _shade_study(ax):
    ax.axvspan(
        STUDY_PERIOD[0] - 0.5,
        STUDY_PERIOD[1] + 0.5,
        color=STUDY_SHADE_COLOR,
        alpha=0.25,
    )


def _plot_series_with_trend(ax, years, values, color, label):
    ax.plot(
        years,
        values,
        "o-",
        color=color,
        lw=2,
        ms=7,
        markeredgecolor="white",
        markeredgewidth=0.5,
        label=label,
    )
    y = np.asarray(values, dtype=float)
    x = np.asarray(years, dtype=float)
    vm = ~np.isnan(y)
    if vm.sum() >= 4:
        z = np.polyfit(x[vm], y[vm], 1)
        p = np.poly1d(z)
        xx = np.arange(int(np.nanmin(x)), int(np.nanmax(x)) + 1)
        ax.plot(xx, p(xx), "--", color=color, alpha=0.6, lw=1.5)


def create_fig_mean_trend(
    df_arc: pd.DataFrame, trend_df: pd.DataFrame, out_png: str
) -> None:
    fig, ax = plt.subplots(1, 1, figsize=(14, 5))
    _shade_study(ax)

    for arc, c in [("inner", COLOR_INNER), ("outer", COLOR_OUTER)]:
        d = df_arc[df_arc["Arc"] == arc].sort_values("year")
        _plot_series_with_trend(
            ax, d["year"].values, d["mean_sst"].values, c, arc.capitalize()
        )

    ax.set_title(
        f"Annual Mean SST (CoralTemp) ({START_YEAR}-{END_YEAR})", fontweight="bold"
    )
    ax.set_xlabel("Year")
    ax.set_ylabel("Mean SST (degC)")
    ax.grid(True, alpha=0.3, linestyle="--")
    ax.set_xlim(START_YEAR - 0.5, END_YEAR + 0.5)
    ax.legend(loc="upper left")

    # annotate MK results
    for arc, y_pos in [("inner", 0.95), ("outer", 0.83)]:
        row = trend_df[(trend_df["Arc"] == arc) & (trend_df["Metric"] == "mean_sst")]
        if len(row):
            row = row.iloc[0]
            ax.text(
                0.99,
                y_pos,
                f"{arc.capitalize()}: {row['significance']} slope={row['sens_slope']:.4g}/yr (p={row['p_value']:.3g})",
                transform=ax.transAxes,
                ha="right",
                va="top",
                fontsize=9,
                color=COLOR_INNER if arc == "inner" else COLOR_OUTER,
                bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
            )

    handles = [
        Line2D([0], [0], color=COLOR_INNER, lw=2, marker="o", label="Inner Arc"),
        Line2D([0], [0], color=COLOR_OUTER, lw=2, marker="o", label="Outer Arc"),
        mpatches.Patch(
            color=STUDY_SHADE_COLOR,
            alpha=0.25,
            label=f"Study period ({STUDY_PERIOD[0]}-{STUDY_PERIOD[1]})",
        ),
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.98, 0.98))
    plt.tight_layout(rect=(0, 0, 0.92, 1))
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def create_fig_variability_multiscale(
    df_arc: pd.DataFrame, trend_df: pd.DataFrame, out_png: str
) -> None:
    metrics = [
        ("msv_2d", "MSV 2d (anom SD)"),
        ("msv_30d", "MSV 30d (anom SD)"),
        ("std_full", "STD full-year (anom)"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), sharex=True)

    for ax, (m, title) in zip(axes, metrics):
        _shade_study(ax)
        for arc, c in [("inner", COLOR_INNER), ("outer", COLOR_OUTER)]:
            d = df_arc[df_arc["Arc"] == arc].sort_values("year")
            _plot_series_with_trend(
                ax, d["year"].values, d[m].values, c, arc.capitalize()
            )

        ax.set_title(title, fontweight="bold")
        ax.grid(True, alpha=0.3, linestyle="--")
        ax.set_xlim(START_YEAR - 0.5, END_YEAR + 0.5)
        ax.set_xlabel("Year")

        for arc, y_pos in [("inner", 0.95), ("outer", 0.83)]:
            row = trend_df[(trend_df["Arc"] == arc) & (trend_df["Metric"] == m)]
            if len(row):
                row = row.iloc[0]
                ax.text(
                    0.99,
                    y_pos,
                    f"{arc.capitalize()}: {row['significance']} slope={row['sens_slope']:.3g}/yr",
                    transform=ax.transAxes,
                    ha="right",
                    va="top",
                    fontsize=9,
                    color=COLOR_INNER if arc == "inner" else COLOR_OUTER,
                    bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
                )

    axes[0].set_ylabel("Anomaly units (degC)")
    axes[0].legend(loc="upper left")

    fig.suptitle(
        f"SST Variability Trends (Deseasonalized; baseline {CLIM_BASELINE_START}-{CLIM_BASELINE_END})",
        fontweight="bold",
        y=1.02,
    )
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def create_fig_extremes(
    df_arc: pd.DataFrame, trend_df: pd.DataFrame, out_png: str
) -> None:
    panels = [
        ("dayfreq_p95", "Daily exceedance freq > P95 (%)"),
        ("events_p95", f"Events > P95 (count; gap>{EVENT_MIN_GAP_DAYS}d)"),
        ("dayfreq_p90", "Daily exceedance freq > P90 (%)"),
        ("events_p90", f"Events > P90 (count; gap>{EVENT_MIN_GAP_DAYS}d)"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), sharex=True)
    axes = axes.flatten()

    for ax, (m, title) in zip(axes, panels):
        _shade_study(ax)
        for arc, c in [("inner", COLOR_INNER), ("outer", COLOR_OUTER)]:
            d = df_arc[df_arc["Arc"] == arc].sort_values("year")
            _plot_series_with_trend(
                ax, d["year"].values, d[m].values, c, arc.capitalize()
            )
        ax.set_title(title, fontweight="bold")
        ax.grid(True, alpha=0.3, linestyle="--")
        ax.set_xlim(START_YEAR - 0.5, END_YEAR + 0.5)

        # annotate only for the main p95 metrics
        if m in ["dayfreq_p95", "events_p95"]:
            for arc, y_pos in [("inner", 0.95), ("outer", 0.83)]:
                row = trend_df[(trend_df["Arc"] == arc) & (trend_df["Metric"] == m)]
                if len(row):
                    row = row.iloc[0]
                    ax.text(
                        0.99,
                        y_pos,
                        f"{arc.capitalize()}: {row['significance']} slope={row['sens_slope']:.3g}/yr",
                        transform=ax.transAxes,
                        ha="right",
                        va="top",
                        fontsize=9,
                        color=COLOR_INNER if arc == "inner" else COLOR_OUTER,
                        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
                    )

    for ax in axes[-2:]:
        ax.set_xlabel("Year")

    axes[0].set_ylabel("% of available days")
    axes[1].set_ylabel("Events / year")
    axes[2].set_ylabel("% of available days")
    axes[3].set_ylabel("Events / year")
    axes[0].legend(loc="upper left")

    fig.suptitle(
        f"SST Extremes Trends (Anomalies; thresholds from baseline {CLIM_BASELINE_START}-{CLIM_BASELINE_END})",
        fontweight="bold",
        y=1.02,
    )
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def create_fig_availability(
    sst_daily_df: pd.DataFrame, sites_df: pd.DataFrame, out_png: str
) -> None:
    """Availability figure: scatter of dates with data by arc."""
    fig, ax = plt.subplots(1, 1, figsize=(16, 3.5))
    fig.suptitle("SST Data Availability (CoralTemp) by Arc", fontweight="bold")

    for y, arc in enumerate(["inner", "outer"]):
        arc_sites = (
            sites_df.loc[sites_df["Arc"] == arc, "Site_HAB"].astype(str).tolist()
        )
        arc_sites = [s for s in arc_sites if s in sst_daily_df.columns]
        if not arc_sites:
            continue
        has_data = sst_daily_df[arc_sites].notna().any(axis=1)
        dates = sst_daily_df.index[has_data]
        ax.scatter(
            dates,
            np.full(len(dates), y),
            s=4,
            marker="|",
            color=COLOR_INNER if arc == "inner" else COLOR_OUTER,
            alpha=0.8,
        )

    ax.set_yticks([0, 1])
    ax.set_yticklabels(["Inner", "Outer"])
    ax.set_xlabel("Date")
    ax.set_ylim(-0.5, 1.5)
    ax.grid(True, axis="x", alpha=0.2)
    plt.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


# =========================================================
# Main
# =========================================================


if __name__ == "__main__":
    print("=" * 70)
    print("SST TEMPORAL TREND ANALYSIS (1985-2025) | Multi-scale + Deseasonalized")
    print("=" * 70)

    if NETCDF_ENGINE:
        print(f"[OK] NetCDF engine: {NETCDF_ENGINE}")
    else:
        print("[WARN] NetCDF engine not detected; continuing with auto-detection.")

    sites_df = load_sites(SITES_CSV_PATH)
    print(f"[OK] Sites loaded: {len(sites_df)}")
    print(
        f"  Inner: {(sites_df['Arc'] == 'inner').sum()} | Outer: {(sites_df['Arc'] == 'outer').sum()}"
    )

    sst_files = filter_files_by_period(
        glob.glob(os.path.join(SST_DIR, SST_PATTERN)), START_YEAR, END_YEAR
    )
    print(f"[OK] SST files (filtered): {len(sst_files)}")
    if len(sst_files) == 0:
        raise SystemExit("No SST files found. Check SST_DIR and SST_PATTERN.")

    # Step 1: Extract / load SST daily table (date x site)
    sst_daily = load_or_build_sst_daily_table(sst_files, sites_df)
    print(
        f"[OK] SST daily table: {sst_daily.shape[0]} days x {sst_daily.shape[1]} sites"
    )

    # Availability figure
    create_fig_availability(
        sst_daily,
        sites_df,
        os.path.join(OUTPUT_DIR, "Fig_S1_Data_Availability_SST.png"),
    )

    # Step 2: Compute annual metrics (site-level)
    df_site = compute_annual_sst_metrics_by_site(sst_daily, sites_df)
    site_csv = os.path.join(OUTPUT_DIR, "Annual_SST_Metrics_by_Site.csv")
    df_site.to_csv(site_csv, index=False)
    print(f"[OK] Saved: {site_csv} ({len(df_site)} rows)")

    # Step 3: Aggregate by arc/year
    df_arc = aggregate_by_arc_year(df_site)
    arc_csv = os.path.join(OUTPUT_DIR, "Annual_SST_Metrics_by_Arc.csv")
    df_arc.to_csv(arc_csv, index=False)
    print(f"[OK] Saved: {arc_csv} ({len(df_arc)} rows)")

    # Step 4: Trend tests
    trend = run_trend_tests(df_arc)
    trend_csv = os.path.join(OUTPUT_DIR, "Trend_Test_Results_SST.csv")
    trend.to_csv(trend_csv, index=False)
    print(f"[OK] Saved: {trend_csv} ({len(trend)} rows)")

    # Step 5: Figures
    create_fig_mean_trend(
        df_arc, trend, os.path.join(OUTPUT_DIR, "Fig_1_SST_Mean_Trend.png")
    )
    create_fig_variability_multiscale(
        df_arc, trend, os.path.join(OUTPUT_DIR, "Fig_2_SST_Variability_Multiscale.png")
    )
    create_fig_extremes(
        df_arc,
        trend,
        os.path.join(OUTPUT_DIR, "Fig_3_SST_Extremes_Frequency_Events.png"),
    )

    print("=" * 70)
    print("DONE")
    print(f"Outputs: {OUTPUT_DIR}")
    print("=" * 70)
