### TIME_SERIES_TEMPORAL_TREND_MULTISCALE_DESEASONALIZED.py ###
#
# PURPOSE
#   Multi-scale and deseasonalized temporal trend analysis (2002-2025).
#   Tests whether environmental variability and extreme-event frequency are
#   increasing over time.
#
# CONTEXT
#   Complements:
#     - TIME_SERIES_VARIABILITY_ANALYSIS_GEMINI.py (cross-sectional Inner vs Outer)
#     - TIME_SERIES_TEMPORAL_TREND_ANALYSIS.py (annual, non-deseasonalized)
#
# CORE IDEAS
#   1) Deseasonalize daily series using a smoothed day-of-year climatology.
#   2) Compute annual metrics from anomalies:
#        - Multi-scale variability (MSV): mean rolling STD at 2d and 30d
#        - Full-year anomaly STD / IQR
#        - Extreme frequency (daily exceedance %) and event count (clustered)
#   3) Test temporal trends in Inner vs Outer arc using Mann-Kendall + Sen slope.
#
# OUTPUTS (written to #######FINAL_RESULTS/Temporal_Trend_Analysis_Multiscale_Deseasonalized)
#   - Annual_Multiscale_Deseasonalized_by_Site.csv
#   - Annual_Multiscale_Deseasonalized_by_Arc.csv
#   - Trend_Test_Results_Multiscale_Deseasonalized.csv
#   - Fig_1_Multiscale_Deseasonalized_Trends.png

import os

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"

import glob
import pickle
import re
import warnings
from typing import Dict, List

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy import stats


warnings.filterwarnings("ignore")


# =========================================================
# CONFIG
# =========================================================

START_YEAR = 2002
END_YEAR = 2025
STUDY_PERIOD = (2006, 2008)

# Cache behavior (pickled dict of site -> pd.Series)
USE_CACHE = True
FORCE_REBUILD_CACHE = False

SST_DIR = r"H:\remote sensing\CRW_SST_FULL"
MODIS_DIR = r"H:\remote sensing\MODIS_DATA_FULL"

SITES_CSV_PATH = r"C:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######22.04.23\DATA\sites_list_full.csv"

OUTPUT_DIR = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        "..",
        "#######FINAL_RESULTS",
        "Temporal_Trend_Analysis_Multiscale_Deseasonalized",
    )
)
CACHE_DIR = os.path.join(OUTPUT_DIR, "cache")
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(CACHE_DIR, exist_ok=True)

SST_PATTERN = "coraltemp_v3.1_*.nc"
KD490_PATTERN = "AQUA_MODIS.*.L3m.DAY.KD.Kd_490.4km.nc"
PAR_PATTERN = "AQUA_MODIS.*.L3m.DAY.PAR.par.4km.nc"
CHL_PATTERN = "AQUA_MODIS.*.L3m.DAY.CHL.chlor_a.4km.nc"

# DLI Calculation Constants (Gattuso et al. 2006)
KDPAR_A, KDPAR_B, KDPAR_C = 0.0665, 0.874, 0.00121

# Palette aligned with TIME_SERIES_VARIABILITY_ANALYSIS_GEMINI.py
COLOR_INNER = "#CD853F"
COLOR_OUTER = "#4169E1"
COLOR_SST = "#D55E00"
COLOR_DLI = "#0072B2"
COLOR_CHL = "#009E73"
STUDY_SHADE_COLOR = "#FFD700"

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
# NetCDF engine detection
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
if NETCDF_ENGINE:
    print(f"[OK] Using NetCDF engine: {NETCDF_ENGINE}")
else:
    print("[WARN] No NetCDF engine found; xarray will auto-detect.")


# =========================================================
# Utility
# =========================================================


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


# =========================================================
# Data loading (site time series)
# =========================================================


def load_timeseries_for_sites(
    files: List[str],
    var_candidates: List[str],
    sites_df: pd.DataFrame,
    label: str,
    cache_name: str,
) -> Dict[str, pd.Series]:
    cache_path = os.path.join(CACHE_DIR, cache_name)
    cached = cache_load(cache_path)
    if cached is not None:
        print(f"[OK] Loaded cached {label}: {cache_path}")
        return cached

    print(f"\nLoading {label} from {len(files)} files...")

    out = {r["Site_HAB"]: {"t": [], "v": []} for _, r in sites_df.iterrows()}

    open_kwargs = {}
    if NETCDF_ENGINE:
        open_kwargs["engine"] = NETCDF_ENGINE

    for i, f in enumerate(files, 1):
        if i % 500 == 0:
            print(f"  {label}: {i}/{len(files)}")
        d = parse_date_from_filename(f)
        if d is None:
            continue
        try:
            with xr.open_dataset(f, **open_kwargs) as ds:
                vname = next((v for v in var_candidates if v in ds.data_vars), None)
                if vname is None:
                    continue
                da = ds[vname]
                if "time" in da.dims:
                    da = da.isel(time=0)

                for _, r in sites_df.iterrows():
                    site_key = r["Site_HAB"]
                    try:
                        val = float(
                            da.sel(
                                lat=r["Latitude"], lon=r["Longitude"], method="nearest"
                            ).values
                        )
                        if np.isnan(val):
                            continue
                        out[site_key]["t"].append(d)
                        out[site_key]["v"].append(val)
                    except Exception:
                        continue
        except Exception:
            continue

    result: Dict[str, pd.Series] = {}
    for s, dct in out.items():
        if dct["t"]:
            result[s] = pd.Series(
                dct["v"], index=pd.DatetimeIndex(dct["t"])
            ).sort_index()

    cache_save(cache_path, result)
    print(f"[OK] {label} sites loaded: {len(result)}")
    return result


def calculate_dli_for_sites(
    kd_files: List[str],
    par_files: List[str],
    sites_df: pd.DataFrame,
    cache_name: str,
) -> Dict[str, pd.Series]:
    cache_path = os.path.join(CACHE_DIR, cache_name)
    cached = cache_load(cache_path)
    if cached is not None:
        print(f"[OK] Loaded cached DLI: {cache_path}")
        return cached

    print(f"\nCalculating DLI from {len(kd_files)} KD files...")

    par_by_date = {}
    for f in par_files:
        d = parse_date_from_filename(f)
        if d is not None:
            par_by_date[d] = f

    out = {r["Site_HAB"]: {"t": [], "v": []} for _, r in sites_df.iterrows()}

    open_kwargs = {}
    if NETCDF_ENGINE:
        open_kwargs["engine"] = NETCDF_ENGINE

    for i, f_kd in enumerate(kd_files, 1):
        if i % 500 == 0:
            print(f"  DLI: {i}/{len(kd_files)}")
        d = parse_date_from_filename(f_kd)
        if d is None or d not in par_by_date:
            continue

        f_par = par_by_date[d]
        try:
            with (
                xr.open_dataset(f_kd, **open_kwargs) as ds_kd,
                xr.open_dataset(f_par, **open_kwargs) as ds_par,
            ):
                kd = ds_kd["Kd_490"]
                par = ds_par["par"]
                if "time" in kd.dims:
                    kd = kd.isel(time=0)
                if "time" in par.dims:
                    par = par.isel(time=0)

                for _, r in sites_df.iterrows():
                    site_key = r["Site_HAB"]
                    try:
                        kdval = float(
                            kd.sel(
                                lat=r["Latitude"], lon=r["Longitude"], method="nearest"
                            ).values
                        )
                        parval = float(
                            par.sel(
                                lat=r["Latitude"], lon=r["Longitude"], method="nearest"
                            ).values
                        )
                        if np.isnan(kdval) or np.isnan(parval) or kdval <= 0:
                            continue

                        kdpar = KDPAR_A + KDPAR_B * kdval - KDPAR_C / (kdval + 1e-9)
                        dli = parval * np.exp(-kdpar * r["Depth_m"])
                        if np.isnan(dli) or dli <= 0:
                            continue

                        out[site_key]["t"].append(d)
                        out[site_key]["v"].append(dli)
                    except Exception:
                        continue
        except Exception:
            continue

    result: Dict[str, pd.Series] = {}
    for s, dct in out.items():
        if dct["t"]:
            result[s] = pd.Series(
                dct["v"], index=pd.DatetimeIndex(dct["t"])
            ).sort_index()

    cache_save(cache_path, result)
    print(f"[OK] DLI sites loaded: {len(result)}")
    return result


# =========================================================
# Deseasonalization + multi-scale metrics
# =========================================================


def deseasonalize_by_doy(series: pd.Series, smooth_window_days: int = 31) -> pd.Series:
    """Deseasonalize by subtracting a smoothed day-of-year climatology.

    Notes:
      - Uses median climatology by DOY (robust to outliers).
      - Applies rolling median smoothing on the DOY climatology.
      - Returns anomalies with the same timestamps as the input series.
    """

    s = series.dropna().sort_index()
    if len(s) < 365:
        return pd.Series(dtype=float)

    doy = s.index.dayofyear.to_numpy()
    doy = np.where(doy == 366, 365, doy)

    clim = pd.DataFrame({"doy": doy, "v": s.values}).groupby("doy")["v"].median()
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
    clim_smooth = pd.Series(sm[pad : pad + 365], index=np.arange(1, 366))

    clim_on_dates = (
        pd.Series(
            np.where(s.index.dayofyear == 366, 365, s.index.dayofyear), index=s.index
        )
        .map(clim_smooth)
        .astype(float)
    )
    return s - clim_on_dates.values


def rolling_std_mean(series: pd.Series, window: int) -> float:
    """Mean rolling standard deviation (centered)."""
    if len(series) < max(5, window):
        return np.nan
    rs = series.rolling(
        window=window, min_periods=max(2, int(window * 0.6)), center=True
    ).std()
    return float(rs.mean())


def count_events_over_threshold(
    series: pd.Series, threshold: float, min_gap_days: int = 2
) -> int:
    """Count clustered exceedances as events (separated by > min_gap_days)."""
    s = series.dropna().sort_index()
    above_idx = s.index[s > threshold]
    if len(above_idx) == 0:
        return 0
    events = 1
    for i in range(1, len(above_idx)):
        if (above_idx[i] - above_idx[i - 1]).days > min_gap_days:
            events += 1
    return int(events)


def compute_annual_multiscale_metrics(
    ts_dict: Dict[str, pd.Series],
    sites_df: pd.DataFrame,
    variable_name: str,
) -> pd.DataFrame:
    rows = []

    for site_key, s_raw in ts_dict.items():
        s_anom = deseasonalize_by_doy(s_raw)
        if len(s_anom) < 365:
            continue

        # Site-specific thresholds on anomalies (global over full record)
        thr90 = float(np.nanpercentile(s_anom.values, 90))
        thr95 = float(np.nanpercentile(s_anom.values, 95))

        site_meta = sites_df.loc[sites_df["Site_HAB"] == site_key]
        if len(site_meta) == 0:
            continue
        site_meta = site_meta.iloc[0]
        site_name = str(site_meta["Site_name"])
        hab = str(site_meta["HAB"])
        arc = str(site_meta["Arc"])

        for year, s_year in s_anom.groupby(s_anom.index.year):
            if year < START_YEAR or year > END_YEAR:
                continue
            if len(s_year) < 60:
                continue

            msv_2d = rolling_std_mean(s_year, 2)
            msv_30d = rolling_std_mean(s_year, 30)
            std_full = float(s_year.std())
            iqr_full = float(s_year.quantile(0.75) - s_year.quantile(0.25))

            dayfreq90 = float((s_year > thr90).mean() * 100.0)
            dayfreq95 = float((s_year > thr95).mean() * 100.0)
            ev90 = count_events_over_threshold(s_year, thr90, min_gap_days=2)
            ev95 = count_events_over_threshold(s_year, thr95, min_gap_days=2)

            rows.append(
                {
                    "Variable": variable_name,
                    "Site_HAB": site_key,
                    "Site_name": site_name,
                    "HAB": hab,
                    "Arc": arc,
                    "year": int(year),
                    "n_obs": int(len(s_year)),
                    "msv_2d": msv_2d,
                    "msv_30d": msv_30d,
                    "std_full": std_full,
                    "iqr_full": iqr_full,
                    "dayfreq_p90": dayfreq90,
                    "dayfreq_p95": dayfreq95,
                    "events_p90": ev90,
                    "events_p95": ev95,
                    "thr90_anom": thr90,
                    "thr95_anom": thr95,
                }
            )

    return pd.DataFrame(rows)


def aggregate_by_arc_year(df_site: pd.DataFrame) -> pd.DataFrame:
    agg_cols = [
        "n_obs",
        "msv_2d",
        "msv_30d",
        "std_full",
        "iqr_full",
        "dayfreq_p90",
        "dayfreq_p95",
        "events_p90",
        "events_p95",
    ]
    out = df_site.groupby(["Variable", "Arc", "year"])[agg_cols].mean().reset_index()
    out["n_sites"] = df_site.groupby(["Variable", "Arc", "year"]).size().values
    return out


# =========================================================
# Trend tests (Mann-Kendall + Sen slope)
# =========================================================


def mann_kendall_sens(y: np.ndarray) -> dict:
    y = np.array([v for v in y if not np.isnan(v)])
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
    metrics = ["msv_2d", "msv_30d", "std_full", "dayfreq_p95", "events_p95"]
    rows = []
    for var in ["SST", "DLI", "Chl-a"]:
        for arc in ["inner", "outer"]:
            sub = df_arc[
                (df_arc["Variable"] == var) & (df_arc["Arc"] == arc)
            ].sort_values("year")
            for m in metrics:
                if m not in sub.columns:
                    continue
                res = mann_kendall_sens(sub[m].values)
                rows.append(
                    {
                        "Variable": var,
                        "Arc": arc,
                        "Metric": m,
                        "n_years": int(len(sub)),
                        **res,
                    }
                )
    return pd.DataFrame(rows)


# =========================================================
# Figures
# =========================================================


def _plot_with_trend(ax, x, y, color, label):
    ax.plot(
        x,
        y,
        "o-",
        color=color,
        lw=2,
        ms=6,
        markeredgecolor="white",
        markeredgewidth=0.5,
        label=label,
    )
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    vm = ~np.isnan(y)
    if vm.sum() >= 4:
        z = np.polyfit(x[vm], y[vm], 1)
        p = np.poly1d(z)
        xx = np.arange(int(np.nanmin(x)), int(np.nanmax(x)) + 1)
        ax.plot(xx, p(xx), "--", color=color, alpha=0.6, lw=1.5)


def create_multiscale_panel(
    df_arc: pd.DataFrame, trend_df: pd.DataFrame, out_png: str
) -> None:
    vars_ = ["SST", "DLI", "Chl-a"]
    metrics = ["msv_2d", "msv_30d", "std_full", "events_p95"]
    labels = [
        "MSV 2d (anom SD)",
        "MSV 30d (anom SD)",
        "STD full-year (anom)",
        "Extreme events P95 (count)",
    ]

    fig, axes = plt.subplots(3, 4, figsize=(18, 12))
    fig.suptitle(
        f"Deseasonalized Multi-scale Trends ({START_YEAR}-{END_YEAR})\n"
        f"Study Period Highlighted: {STUDY_PERIOD[0]}-{STUDY_PERIOD[1]}",
        fontweight="bold",
    )

    for i, v in enumerate(vars_):
        for j, m in enumerate(metrics):
            ax = axes[i, j]
            ax.axvspan(
                STUDY_PERIOD[0] - 0.5,
                STUDY_PERIOD[1] + 0.5,
                color=STUDY_SHADE_COLOR,
                alpha=0.25,
            )

            for arc, c in [("inner", COLOR_INNER), ("outer", COLOR_OUTER)]:
                d = df_arc[
                    (df_arc["Variable"] == v) & (df_arc["Arc"] == arc)
                ].sort_values("year")
                if len(d) == 0:
                    continue
                _plot_with_trend(ax, d["year"].values, d[m].values, c, arc.capitalize())

            # Annotate trend test summaries (inner and outer)
            t = trend_df[
                (trend_df["Variable"] == v)
                & (trend_df["Arc"] == "inner")
                & (trend_df["Metric"] == m)
            ]
            if len(t):
                ax.text(
                    0.02,
                    0.98,
                    f"Inner {t.iloc[0]['significance']} slope={t.iloc[0]['sens_slope']:.3g}",
                    transform=ax.transAxes,
                    va="top",
                    fontsize=8,
                    color=COLOR_INNER,
                )

            t = trend_df[
                (trend_df["Variable"] == v)
                & (trend_df["Arc"] == "outer")
                & (trend_df["Metric"] == m)
            ]
            if len(t):
                ax.text(
                    0.02,
                    0.88,
                    f"Outer {t.iloc[0]['significance']} slope={t.iloc[0]['sens_slope']:.3g}",
                    transform=ax.transAxes,
                    va="top",
                    fontsize=8,
                    color=COLOR_OUTER,
                )

            ax.grid(True, alpha=0.3, linestyle="--")
            ax.set_xlim(START_YEAR - 0.5, END_YEAR + 0.5)
            if i == 0:
                ax.set_title(labels[j], fontweight="bold")
            if j == 0:
                ax.set_ylabel(v, fontweight="bold")

    handles = [
        Line2D([0], [0], color=COLOR_INNER, lw=2, marker="o", label="Inner Arc"),
        Line2D([0], [0], color=COLOR_OUTER, lw=2, marker="o", label="Outer Arc"),
        mpatches.Patch(
            color=STUDY_SHADE_COLOR,
            alpha=0.25,
            label=f"Study period ({STUDY_PERIOD[0]}-{STUDY_PERIOD[1]})",
        ),
    ]
    fig.legend(handles=handles, loc="upper right")
    plt.tight_layout(rect=(0, 0, 0.93, 0.95))
    fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()


# =========================================================
# MAIN
# =========================================================


if __name__ == "__main__":
    print("=" * 70)
    print("DESEASONALIZED MULTISCALE TEMPORAL TREND ANALYSIS")
    print("=" * 70)

    sites = load_sites(SITES_CSV_PATH)
    print(
        f"[OK] Sites: {len(sites)} | Inner={(sites['Arc'] == 'inner').sum()} | Outer={(sites['Arc'] == 'outer').sum()}"
    )

    sst_files = filter_files_by_period(
        glob.glob(os.path.join(SST_DIR, SST_PATTERN)), START_YEAR, END_YEAR
    )
    kd_files = filter_files_by_period(
        glob.glob(os.path.join(MODIS_DIR, KD490_PATTERN)), START_YEAR, END_YEAR
    )
    par_files = filter_files_by_period(
        glob.glob(os.path.join(MODIS_DIR, PAR_PATTERN)), START_YEAR, END_YEAR
    )
    chl_files = filter_files_by_period(
        glob.glob(os.path.join(MODIS_DIR, CHL_PATTERN)), START_YEAR, END_YEAR
    )

    print(
        f"[OK] Files | SST={len(sst_files)} KD={len(kd_files)} PAR={len(par_files)} CHL={len(chl_files)}"
    )

    sst_ts = load_timeseries_for_sites(
        sst_files,
        ["analysed_sst", "sst"],
        sites,
        "SST",
        "cache_sst_site_hab.pkl",
    )
    chl_ts = load_timeseries_for_sites(
        chl_files,
        ["chlor_a"],
        sites,
        "Chl-a",
        "cache_chl_site_hab.pkl",
    )
    dli_ts = calculate_dli_for_sites(
        kd_files, par_files, sites, "cache_dli_site_hab.pkl"
    )

    sst_site = compute_annual_multiscale_metrics(sst_ts, sites, "SST")
    dli_site = compute_annual_multiscale_metrics(dli_ts, sites, "DLI")
    chl_site = compute_annual_multiscale_metrics(chl_ts, sites, "Chl-a")

    site_all = pd.concat([sst_site, dli_site, chl_site], ignore_index=True)
    site_csv = os.path.join(OUTPUT_DIR, "Annual_Multiscale_Deseasonalized_by_Site.csv")
    site_all.to_csv(site_csv, index=False)
    print(f"[OK] Saved site metrics: {site_csv} ({len(site_all)} rows)")

    arc_all = aggregate_by_arc_year(site_all)
    arc_csv = os.path.join(OUTPUT_DIR, "Annual_Multiscale_Deseasonalized_by_Arc.csv")
    arc_all.to_csv(arc_csv, index=False)
    print(f"[OK] Saved arc metrics: {arc_csv} ({len(arc_all)} rows)")

    trend = run_trend_tests(arc_all)
    trend_csv = os.path.join(
        OUTPUT_DIR, "Trend_Test_Results_Multiscale_Deseasonalized.csv"
    )
    trend.to_csv(trend_csv, index=False)
    print(f"[OK] Saved trend tests: {trend_csv} ({len(trend)} rows)")

    fig_png = os.path.join(OUTPUT_DIR, "Fig_1_Multiscale_Deseasonalized_Trends.png")
    create_multiscale_panel(arc_all, trend, fig_png)
    print(f"[OK] Saved figure: {fig_png}")

    print("=" * 70)
    print("DONE")
    print(f"Outputs: {OUTPUT_DIR}")
    print("=" * 70)
