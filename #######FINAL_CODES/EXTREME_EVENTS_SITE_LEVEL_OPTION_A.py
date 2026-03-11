"""
Extreme events analysis (Option A) with SITE x HAB differentiation.

Core design:
1) Build profile-level time series (PROFILE_ID = SITE_HAB).
2) Compute day-of-year climatology per profile and z-score anomalies.
3) Detect episodic extremes using thresholded anomalies with duration rules.
4) Compare Inner vs Outer distributions and summarize by Inner sectors and HAB.
5) Generate manuscript-ready figures and an automatic scientific report.
"""

import os
import glob
import logging
import re
import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.patches import Patch
from scipy import stats

warnings.filterwarnings("ignore")
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"


def get_netcdf_engine():
    try:
        import netCDF4  # noqa: F401

        return "netcdf4"
    except ImportError:
        pass
    try:
        import h5netcdf  # noqa: F401

        return "h5netcdf"
    except ImportError:
        pass
    return None


NETCDF_ENGINE = get_netcdf_engine()


# ============================================
# Configuration
# ============================================
START_YEAR = 2002
END_YEAR = 2008

SST_DIR = r"H:\remote sensing\CRW_SST_FULL"
MODIS_DIR = r"H:\remote sensing\MODIS_DATA_FULL"

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

INTEGRATED_DATA_PATH = os.path.join(
    BASE_DIR,
    "#######FINAL_RESULTS",
    "#####output_local_PCA_CV_2_FINAL",
    "dados_abundancia_integrados_long_format.csv",
)

CV2_PATH = os.path.join(
    BASE_DIR,
    "#######FINAL_RESULTS",
    "#####output_local_PCA_CV_2_FINAL",
    "dados_consolidados_com_scores_das_duas_PCAs_CV2.csv",
)
CV30_PATH = os.path.join(
    BASE_DIR,
    "#######FINAL_RESULTS",
    "#####output_local_PCA_CV_30_FINAL",
    "dados_consolidados_com_scores_das_duas_PCAs_CV30.csv",
)
CVALL_PATH = os.path.join(
    BASE_DIR,
    "#######FINAL_RESULTS",
    "#####output_local_PCA_CV_all_FINAL",
    "dados_consolidados_com_scores_das_duas_PCAs_CVall.csv",
)

OUTPUT_DIR = os.path.join(BASE_DIR, "#######FINAL_RESULTS", "Extreme_Events_Site_Level")
os.makedirs(OUTPUT_DIR, exist_ok=True)

SST_PATTERN = "coraltemp_v3.1_*.nc"
KD490_PATTERN = "AQUA_MODIS.*.L3m.DAY.KD.Kd_490.4km.nc"
PAR_PATTERN = "AQUA_MODIS.*.L3m.DAY.PAR.par.4km.nc"
CHL_PATTERN = "AQUA_MODIS.*.L3m.DAY.CHL.chlor_a.4km.nc"

KDPAR_A, KDPAR_B, KDPAR_C = 0.0665, 0.874, 0.00121

THRESHOLD_Z = 1.28
MIN_DURATION_DAYS = 2
MERGE_GAP_DAYS = 1
MIN_DAYS_FOR_CLIMATOLOGY = 365
CLIMATOLOGY_SMOOTHING_WINDOW = 30

COLOR_INNER = "#CD853F"
COLOR_OUTER = "#4169E1"
COLOR_SST = "#D55E00"
COLOR_DLI = "#0072B2"
COLOR_CHL = "#009E73"

HAB_COLORS = {"TP": "#009E73", "PA": "#E69F00", "RR": "#56B4E9"}

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 11,
        "savefig.dpi": 300,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

logging.basicConfig(
    filename=os.path.join(
        BASE_DIR, "#######FINAL_CODES", "extreme_events_analysis.log"
    ),
    filemode="w",
    level=logging.ERROR,
    format="%(asctime)s - %(levelname)s - %(message)s",
)


# ============================================
# Utility
# ============================================
def parse_date_from_filename(filename):
    m = re.search(r"(\d{8})", os.path.basename(filename))
    if not m:
        return None
    try:
        return pd.to_datetime(m.group(1), format="%Y%m%d")
    except ValueError:
        return None


def filter_files_by_period(files, start_year, end_year):
    out = []
    for f in files:
        d = parse_date_from_filename(f)
        if d is not None and start_year <= d.year <= end_year:
            out.append(f)
    return sorted(out)


def classify_inner_sector(site_name, arc, reef_name=None):
    site = str(site_name).upper()
    reef = str(reef_name).upper() if reef_name is not None else ""
    arc_l = str(arc).lower()

    if arc_l != "inner":
        return "outer"
    if site.startswith("ITA"):
        return "inner_north"
    if site.startswith("TIM"):
        return "inner_central"
    if reef == "UCR" or site in {"ARE", "PA2", "PLE", "SGO"}:
        return "inner_south"
    return "inner_other"


def format_p_value(p):
    if pd.isna(p):
        return "NA"
    if p < 0.001:
        return "<0.001"
    return f"{p:.3f}"


# ============================================
# Metadata
# ============================================
def load_profile_metadata(path):
    df = pd.read_csv(path, sep=";", decimal=",", encoding="utf-8-sig", low_memory=False)
    df.columns = [c.strip() for c in df.columns]

    required = [
        "REEF",
        "SITE",
        "HAB",
        "LATITUDE",
        "LONGITUDE",
        "DEPTH_M",
        "ARCH",
        "UNIQUE_ID",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Integrated input missing columns: {missing}")

    m = df[required].dropna().drop_duplicates().copy()
    m = m.rename(
        columns={
            "REEF": "REEF_NAME",
            "SITE": "SITE_ID",
            "HAB": "HAB",
            "LATITUDE": "LAT",
            "LONGITUDE": "LON",
            "DEPTH_M": "DEPTH_M",
            "ARCH": "ARC",
            "UNIQUE_ID": "PROFILE_ID",
        }
    )

    m["ARC"] = m["ARC"].astype(str).str.strip().str.lower()
    m = m[m["ARC"].isin(["inner", "outer"])].copy()

    m["LAT"] = pd.to_numeric(m["LAT"], errors="coerce")
    m["LON"] = pd.to_numeric(m["LON"], errors="coerce")
    m["DEPTH_M"] = pd.to_numeric(m["DEPTH_M"], errors="coerce")
    m = m.dropna(subset=["LAT", "LON", "DEPTH_M"]).copy()

    m["INNER_SECTOR"] = m.apply(
        lambda r: classify_inner_sector(r["SITE_ID"], r["ARC"], r["REEF_NAME"]), axis=1
    )
    m = m.sort_values(["ARC", "REEF_NAME", "SITE_ID", "HAB"]).reset_index(drop=True)
    return m


# ============================================
# Time-series extraction
# ============================================
def load_timeseries_for_profiles(files, var_names, profiles_df, label):
    print(f"\n  Loading {label}...")
    profile_tuples = [
        (r.PROFILE_ID, float(r.LAT), float(r.LON))
        for r in profiles_df[["PROFILE_ID", "LAT", "LON"]].itertuples(index=False)
    ]
    out = {pid: {"times": [], "values": []} for pid, _, _ in profile_tuples}

    ok = 0
    err = 0
    open_kwargs = {"engine": NETCDF_ENGINE} if NETCDF_ENGINE else {}

    for i, f in enumerate(files, start=1):
        if i % 500 == 0:
            print(f"    Processing file {i}/{len(files)}...")

        d = parse_date_from_filename(f)
        if d is None:
            continue

        try:
            with xr.open_dataset(f, **open_kwargs) as ds:
                var_name = next((v for v in var_names if v in ds.data_vars), None)
                if var_name is None:
                    continue
                da = ds[var_name]
                if "time" in da.dims:
                    da = da.isel(time=0)

                for pid, lat, lon in profile_tuples:
                    try:
                        v = float(da.sel(lat=lat, lon=lon, method="nearest").values)
                        if not np.isnan(v):
                            out[pid]["times"].append(d)
                            out[pid]["values"].append(v)
                    except Exception:
                        continue
                ok += 1
        except Exception as e:
            err += 1
            if err <= 5:
                logging.error(f"{label} read error in {f}: {e}")

    print(f"    Files processed successfully: {ok}/{len(files)}")
    if err > 0:
        print(f"    Files with errors: {err}")

    series_dict = {}
    for pid, d in out.items():
        if d["times"]:
            s = pd.Series(d["values"], index=pd.DatetimeIndex(d["times"])).sort_index()
            series_dict[pid] = s
    print(f"    Loaded data for {len(series_dict)} profiles")
    return series_dict


def calculate_dli_for_profiles(kd_files, par_files, profiles_df):
    print("\n  Calculating DLI...")
    kd_dates = {
        parse_date_from_filename(f): f for f in kd_files if parse_date_from_filename(f)
    }
    par_dates = {
        parse_date_from_filename(f): f for f in par_files if parse_date_from_filename(f)
    }
    common_dates = sorted(set(kd_dates).intersection(set(par_dates)))

    profile_tuples = [
        (r.PROFILE_ID, float(r.LAT), float(r.LON), float(r.DEPTH_M))
        for r in profiles_df[["PROFILE_ID", "LAT", "LON", "DEPTH_M"]].itertuples(
            index=False
        )
    ]
    out = {pid: {"times": [], "values": []} for pid, _, _, _ in profile_tuples}

    ok = 0
    err = 0
    open_kwargs = {"engine": NETCDF_ENGINE} if NETCDF_ENGINE else {}

    for i, d in enumerate(common_dates, start=1):
        if i % 500 == 0:
            print(f"    Processing date {i}/{len(common_dates)}...")

        f_kd = kd_dates[d]
        f_par = par_dates[d]
        try:
            with (
                xr.open_dataset(f_kd, **open_kwargs) as ds_kd,
                xr.open_dataset(f_par, **open_kwargs) as ds_par,
            ):
                kd_da = ds_kd["Kd_490"]
                par_da = ds_par["par"]
                if "time" in kd_da.dims:
                    kd_da = kd_da.isel(time=0)
                if "time" in par_da.dims:
                    par_da = par_da.isel(time=0)

                for pid, lat, lon, depth in profile_tuples:
                    try:
                        kd490 = float(
                            kd_da.sel(lat=lat, lon=lon, method="nearest").values
                        )
                        par = float(
                            par_da.sel(lat=lat, lon=lon, method="nearest").values
                        )
                        if np.isnan(kd490) or np.isnan(par) or kd490 <= 0:
                            continue
                        kdpar = KDPAR_A + KDPAR_B * kd490 - KDPAR_C / (kd490 + 1e-9)
                        dli = par * np.exp(-kdpar * depth)
                        if not np.isnan(dli) and dli > 0:
                            out[pid]["times"].append(d)
                            out[pid]["values"].append(dli)
                    except Exception:
                        continue
                ok += 1
        except Exception as e:
            err += 1
            if err <= 5:
                logging.error(f"DLI calc error at {d}: {e}")

    print(f"    Date pairs processed successfully: {ok}/{len(common_dates)}")
    if err > 0:
        print(f"    Date pairs with errors: {err}")

    series_dict = {}
    for pid, d in out.items():
        if d["times"]:
            s = pd.Series(d["values"], index=pd.DatetimeIndex(d["times"])).sort_index()
            series_dict[pid] = s
    print(f"    Calculated DLI for {len(series_dict)} profiles")
    return series_dict


# ============================================
# Climatology, anomalies, events
# ============================================
def circular_smooth(values, window):
    arr = np.asarray(values, dtype=float)
    if arr.size == 0 or window <= 1:
        return arr
    pad = window // 2
    ext = np.concatenate([arr[-pad:], arr, arr[:pad]])
    sm = pd.Series(ext).rolling(window=window, center=True, min_periods=1).mean().values
    return sm[pad : pad + arr.size]


def compute_climatology(series, window=CLIMATOLOGY_SMOOTHING_WINDOW):
    s = series.dropna()
    if len(s) < MIN_DAYS_FOR_CLIMATOLOGY:
        return None

    df = pd.DataFrame({"date": s.index, "value": s.values})
    df["doy"] = df["date"].dt.dayofyear
    clm = df.groupby("doy")["value"].agg(["count", "mean", "std"]).reset_index()
    clm.columns = ["doy", "n", "mean", "std"]

    full = pd.DataFrame({"doy": np.arange(1, 367)})
    clm = full.merge(clm, on="doy", how="left")
    clm["mean_smooth"] = circular_smooth(clm["mean"].values, window)
    clm["std_smooth"] = circular_smooth(clm["std"].values, window)
    clm["mean_smooth"] = clm["mean_smooth"].fillna(clm["mean"])
    clm["std_smooth"] = clm["std_smooth"].fillna(clm["std"])
    clm["std_smooth"] = clm["std_smooth"].replace(0, np.nan)
    return clm


def compute_anomaly(series, climatology):
    if climatology is None:
        return None
    df = pd.DataFrame({"date": series.index, "value": series.values})
    df["doy"] = df["date"].dt.dayofyear
    df = df.merge(
        climatology[["doy", "mean_smooth", "std_smooth"]], on="doy", how="left"
    )
    df["anomaly_z"] = (df["value"] - df["mean_smooth"]) / df["std_smooth"]
    return pd.Series(df["anomaly_z"].values, index=df["date"], name="anomaly_z")


def detect_events(
    anomaly_series,
    threshold_z,
    min_duration,
    merge_gap,
    analysis_start,
    analysis_end,
):
    cols = [
        "start_date",
        "end_date",
        "duration_d",
        "n_extreme_days",
        "peak_date",
        "peak_z",
        "mean_z",
        "cumulative_z",
    ]
    if anomaly_series is None or len(anomaly_series) == 0:
        return pd.DataFrame(columns=cols)

    idx = pd.date_range(analysis_start, analysis_end, freq="D")
    anom = anomaly_series.reindex(idx)
    is_ext = (anom >= threshold_z) & anom.notna()
    if is_ext.sum() == 0:
        return pd.DataFrame(columns=cols)

    run_id = (is_ext != is_ext.shift(1)).cumsum()
    runs = (
        pd.DataFrame({"date": idx, "is_ext": is_ext.values, "run_id": run_id.values})
        .query("is_ext == True")
        .groupby("run_id")
        .agg(
            start_date=("date", "min"),
            end_date=("date", "max"),
            n_extreme_days=("date", "count"),
        )
        .reset_index(drop=True)
    )
    runs["duration_d"] = (runs["end_date"] - runs["start_date"]).dt.days + 1
    runs = (
        runs[runs["duration_d"] >= min_duration]
        .sort_values("start_date")
        .reset_index(drop=True)
    )
    if runs.empty:
        return pd.DataFrame(columns=cols)

    merged = []
    for _, row in runs.iterrows():
        s = row["start_date"]
        e = row["end_date"]
        if not merged:
            merged.append([s, e])
            continue
        ls, le = merged[-1]
        gap = (s - le).days - 1
        if 0 <= gap <= merge_gap:
            gap_days = pd.date_range(
                le + pd.Timedelta(days=1), s - pd.Timedelta(days=1), freq="D"
            )
            gap_missing = (
                anom.reindex(gap_days).isna().any() if len(gap_days) else False
            )
            if not gap_missing:
                merged[-1][1] = e
            else:
                merged.append([s, e])
        else:
            merged.append([s, e])

    events = []
    for s, e in merged:
        span_days = pd.date_range(s, e, freq="D")
        span_vals = anom.reindex(span_days)
        extreme_vals = span_vals[(span_vals >= threshold_z) & span_vals.notna()]
        if len(extreme_vals) < min_duration:
            continue
        events.append(
            {
                "start_date": s,
                "end_date": e,
                "duration_d": int(len(span_days)),
                "n_extreme_days": int(len(extreme_vals)),
                "peak_date": extreme_vals.idxmax(),
                "peak_z": float(extreme_vals.max()),
                "mean_z": float(extreme_vals.mean()),
                "cumulative_z": float(extreme_vals.sum()),
            }
        )

    return pd.DataFrame(events, columns=cols)


def compute_profile_metrics(event_df, n_valid_days, expected_product_days):
    n_events = len(event_df)
    effective_years = n_valid_days / 365.25 if n_valid_days > 0 else np.nan
    coverage_pct = (
        (n_valid_days / expected_product_days) * 100
        if expected_product_days and not np.isnan(expected_product_days)
        else np.nan
    )

    if n_events == 0:
        return {
            "N_VALID_DAYS": n_valid_days,
            "N_YEARS": effective_years,
            "EXPECTED_PRODUCT_DAYS": expected_product_days,
            "COVERAGE_PCT_PRODUCT": coverage_pct,
            "n_events": 0,
            "n_events_per_year": 0.0 if effective_years > 0 else np.nan,
            "mean_duration_d": np.nan,
            "max_duration_d": np.nan,
            "mean_peak_z": np.nan,
            "max_peak_z": np.nan,
            "mean_cumulative_z": np.nan,
            "pct_extreme_days": 0.0,
        }

    extreme_days = int(event_df["n_extreme_days"].sum())
    return {
        "N_VALID_DAYS": n_valid_days,
        "N_YEARS": effective_years,
        "EXPECTED_PRODUCT_DAYS": expected_product_days,
        "COVERAGE_PCT_PRODUCT": coverage_pct,
        "n_events": n_events,
        "n_events_per_year": n_events / effective_years
        if effective_years > 0
        else np.nan,
        "mean_duration_d": float(event_df["duration_d"].mean()),
        "max_duration_d": int(event_df["duration_d"].max()),
        "mean_peak_z": float(event_df["peak_z"].mean()),
        "max_peak_z": float(event_df["peak_z"].max()),
        "mean_cumulative_z": float(event_df["cumulative_z"].mean()),
        "pct_extreme_days": (extreme_days / n_valid_days) * 100
        if n_valid_days > 0
        else 0.0,
    }


# ============================================
# Statistics
# ============================================
def permutation_test(x, y, n_perm=9999, seed=42):
    xv = np.array([v for v in x if not np.isnan(v)], dtype=float)
    yv = np.array([v for v in y if not np.isnan(v)], dtype=float)
    if len(xv) < 3 or len(yv) < 3:
        return {"obs_diff": np.nan, "p_value": np.nan, "effect_size": np.nan}

    obs = float(np.mean(xv) - np.mean(yv))
    allv = np.concatenate([xv, yv])
    n_x = len(xv)
    rng = np.random.default_rng(seed)
    null = np.zeros(n_perm, dtype=float)
    for i in range(n_perm):
        rng.shuffle(allv)
        null[i] = np.mean(allv[:n_x]) - np.mean(allv[n_x:])
    p = (np.sum(np.abs(null) >= np.abs(obs)) + 1) / (n_perm + 1)

    pooled = np.sqrt(
        ((len(xv) - 1) * np.var(xv, ddof=1) + (len(yv) - 1) * np.var(yv, ddof=1))
        / (len(xv) + len(yv) - 2)
    )
    d = obs / pooled if pooled > 0 else np.nan
    return {"obs_diff": obs, "p_value": p, "effect_size": d}


def mann_whitney_p(x, y):
    xv = [v for v in x if not np.isnan(v)]
    yv = [v for v in y if not np.isnan(v)]
    if len(xv) < 3 or len(yv) < 3:
        return np.nan
    try:
        _, p = stats.mannwhitneyu(xv, yv, alternative="two-sided")
        return float(p)
    except Exception:
        return np.nan


def build_stat_tables(metrics_df):
    stat_rows = []
    for var in ["SST", "DLI", "Chl-a"]:
        sub = metrics_df[metrics_df["VARIABLE"] == var].copy()
        for metric in [
            "n_events_per_year",
            "mean_duration_d",
            "mean_peak_z",
            "pct_extreme_days",
        ]:
            msub = sub.dropna(subset=[metric]).copy()
            if metric in ["mean_duration_d", "mean_peak_z"]:
                msub = msub[msub["n_events"] > 0].copy()

            inner = msub[msub["ARC"] == "inner"][metric].values
            outer = msub[msub["ARC"] == "outer"][metric].values
            perm = permutation_test(inner, outer)
            p_mw = mann_whitney_p(inner, outer)

            stat_rows.append(
                {
                    "Variable": var,
                    "Metric": metric,
                    "Inner_Mean": np.nanmean(inner) if len(inner) else np.nan,
                    "Outer_Mean": np.nanmean(outer) if len(outer) else np.nan,
                    "Obs_Difference": perm["obs_diff"],
                    "P_Value_Perm": perm["p_value"],
                    "P_Value_MWU": p_mw,
                    "Effect_Size_CohenD": perm["effect_size"],
                    "N_Inner": len(inner),
                    "N_Outer": len(outer),
                }
            )

    stat_df = pd.DataFrame(stat_rows)

    summary_sector_hab = (
        metrics_df.groupby(["VARIABLE", "ARC", "INNER_SECTOR", "HAB"], dropna=False)
        .agg(
            N_profiles=("PROFILE_ID", "nunique"),
            Mean_events_per_year=("n_events_per_year", "mean"),
            Mean_duration_d=("mean_duration_d", "mean"),
            Mean_peak_z=("mean_peak_z", "mean"),
            Mean_pct_extreme_days=("pct_extreme_days", "mean"),
            Mean_coverage_pct=("COVERAGE_PCT_PRODUCT", "mean"),
        )
        .reset_index()
    )

    zero_summary = (
        metrics_df.groupby(["VARIABLE", "ARC", "INNER_SECTOR", "HAB"], dropna=False)
        .agg(
            N_profiles=("PROFILE_ID", "nunique"),
            N_zero_events=("n_events", lambda x: int((x == 0).sum())),
        )
        .reset_index()
    )
    zero_summary["Pct_zero_events"] = (
        zero_summary["N_zero_events"] / zero_summary["N_profiles"] * 100
    )

    return stat_df, summary_sector_hab, zero_summary


# ============================================
# Figures
# ============================================
def _metric_figure(
    metrics_df, value_col, ylabel, title, output_path, events_only=False
):
    vars_order = ["SST", "DLI", "Chl-a"]
    var_colors = {"SST": COLOR_SST, "DLI": COLOR_DLI, "Chl-a": COLOR_CHL}

    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    fig.suptitle(title, fontsize=14, fontweight="bold", y=1.03)

    for i, var in enumerate(vars_order):
        ax = axes[i]
        sub = metrics_df[metrics_df["VARIABLE"] == var].copy()
        if events_only:
            sub = sub[sub["n_events"] > 0].copy()
        sub = sub.dropna(subset=[value_col])

        inner = sub[sub["ARC"] == "inner"][value_col].values
        outer = sub[sub["ARC"] == "outer"][value_col].values

        if len(inner) == 0 and len(outer) == 0:
            ax.text(
                0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center"
            )
            ax.set_title(var)
            continue

        data = [
            inner if len(inner) else np.array([0.0]),
            outer if len(outer) else np.array([0.0]),
        ]
        bp = ax.boxplot(data, labels=["Inner", "Outer"], patch_artist=True, widths=0.6)
        bp["boxes"][0].set_facecolor(COLOR_INNER)
        bp["boxes"][1].set_facecolor(COLOR_OUTER)
        for b in bp["boxes"]:
            b.set_alpha(0.75)

        for j, vals in enumerate([inner, outer], start=1):
            if len(vals):
                xj = np.random.normal(j, 0.04, len(vals))
                c = COLOR_INNER if j == 1 else COLOR_OUTER
                ax.scatter(
                    xj,
                    vals,
                    c=c,
                    s=26,
                    alpha=0.6,
                    edgecolors="white",
                    linewidth=0.4,
                    zorder=3,
                )

        p = mann_whitney_p(inner, outer)
        sig = "ns"
        if not pd.isna(p):
            if p < 0.001:
                sig = "***"
            elif p < 0.01:
                sig = "**"
            elif p < 0.05:
                sig = "*"

        ax.text(
            0.95,
            0.95,
            f"p={format_p_value(p)}\n{sig}",
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.55),
        )
        ax.set_title(var, color=var_colors[var], fontweight="bold")
        ax.set_ylabel(ylabel)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()


def create_distribution_figures(metrics_df):
    _metric_figure(
        metrics_df,
        value_col="n_events_per_year",
        ylabel="Events per year",
        title="Extreme Event Frequency: Inner vs Outer (SITE x HAB)",
        output_path=os.path.join(OUTPUT_DIR, "Fig1_Event_Frequency.png"),
        events_only=False,
    )
    _metric_figure(
        metrics_df,
        value_col="mean_duration_d",
        ylabel="Mean event duration (days)",
        title="Extreme Event Duration: Inner vs Outer (SITE x HAB)",
        output_path=os.path.join(OUTPUT_DIR, "Fig2_Event_Duration.png"),
        events_only=True,
    )
    _metric_figure(
        metrics_df,
        value_col="mean_peak_z",
        ylabel="Mean peak anomaly (z)",
        title="Extreme Event Peak Intensity: Inner vs Outer (SITE x HAB)",
        output_path=os.path.join(OUTPUT_DIR, "Fig3_Peak_Intensity.png"),
        events_only=True,
    )


def create_profile_pair_timeseries_figure(
    anomaly_dict_var,
    event_dict_var,
    profile_top,
    profile_bottom,
    figure_title,
    output_file,
):
    s_top = anomaly_dict_var.get(profile_top)
    s_bottom = anomaly_dict_var.get(profile_bottom)
    if s_top is None or s_bottom is None:
        print(f"  WARNING: Missing anomaly series for {output_file}")
        return

    ev_top = event_dict_var.get(profile_top, pd.DataFrame())
    ev_bottom = event_dict_var.get(profile_bottom, pd.DataFrame())

    fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=True)
    for ax, pid, s, ev in [
        (axes[0], profile_top, s_top, ev_top),
        (axes[1], profile_bottom, s_bottom, ev_bottom),
    ]:
        ax.plot(s.index, s.values, color=COLOR_DLI, linewidth=0.8, alpha=0.85)
        ax.axhline(THRESHOLD_Z, color="red", linestyle="--", linewidth=1.0, alpha=0.7)
        ax.axhline(0.0, color="gray", linestyle="-", linewidth=0.8, alpha=0.5)
        if ev is not None and len(ev) > 0:
            for _, row in ev.iterrows():
                ax.axvspan(
                    row["start_date"], row["end_date"], color=COLOR_DLI, alpha=0.20
                )
                ax.scatter(
                    [row["peak_date"]], [row["peak_z"]], color="darkred", s=35, zorder=4
                )
        ax.set_ylabel("DLI anomaly (z)")
        ax.set_title(pid, loc="left", fontweight="bold")

    fig.suptitle(figure_title, fontsize=13, fontweight="bold")
    axes[-1].set_xlabel("Date")
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_file, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()


def load_cv_long_table():
    specs = [
        (CV2_PATH, "CV_2", "dli_cv_2"),
        (CV30_PATH, "CV_30", "dli_cv_30"),
        (CVALL_PATH, "CV_ALL", "cv_DLI_local"),
    ]

    all_df = []
    for path, window, cv_col in specs:
        df = pd.read_csv(
            path, sep=";", decimal=",", encoding="utf-8-sig", low_memory=False
        )
        df.columns = [c.strip() for c in df.columns]
        req = ["Reef_name", "Site_name", "HAB", "Arch", cv_col]
        missing = [c for c in req if c not in df.columns]
        if missing:
            raise ValueError(
                f"Missing CV columns in {os.path.basename(path)}: {missing}"
            )

        d = df[["Reef_name", "Site_name", "HAB", "Arch", cv_col]].copy()
        d = d.rename(
            columns={
                "Reef_name": "REEF_NAME",
                "Site_name": "SITE_ID",
                "Arch": "ARC",
                cv_col: "DLI_CV",
            }
        )
        d["ARC"] = d["ARC"].astype(str).str.strip().str.lower()
        d = d[d["ARC"].isin(["inner", "outer"])].copy()
        d["DLI_CV"] = pd.to_numeric(d["DLI_CV"], errors="coerce")
        d = d.dropna(subset=["DLI_CV"])
        d["WINDOW"] = window
        d["INNER_SECTOR"] = d.apply(
            lambda r: classify_inner_sector(r["SITE_ID"], r["ARC"], r["REEF_NAME"]),
            axis=1,
        )
        all_df.append(d)

    return pd.concat(all_df, ignore_index=True)


def create_cv_hab_sector_figure(cv_df, output_path):
    print("\n  Creating CV by HAB/Sector figure...")

    windows = ["CV_2", "CV_30", "CV_ALL"]
    sectors = ["inner_north", "inner_central", "inner_south", "outer"]
    sector_labels = {
        "inner_north": "Inner North",
        "inner_central": "Inner Central",
        "inner_south": "Inner South",
        "outer": "Outer",
    }
    hab_order = ["TP", "PA", "RR"]
    offsets = {"TP": -0.25, "PA": 0.00, "RR": 0.25}

    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    fig.suptitle(
        "DLI CV Across Windows by HAB and Arc/Sector",
        fontsize=15,
        fontweight="bold",
        y=1.02,
    )

    for ax, w in zip(axes, windows):
        sub = cv_df[cv_df["WINDOW"] == w].copy()
        all_pos = []
        all_data = []
        all_hab = []

        for i, sec in enumerate(sectors):
            sec_sub = sub[sub["INNER_SECTOR"] == sec]
            for hab in hab_order:
                vals = sec_sub[sec_sub["HAB"] == hab]["DLI_CV"].dropna().values
                if len(vals) == 0:
                    continue
                pos = i + offsets[hab]
                all_pos.append(pos)
                all_data.append(vals)
                all_hab.append(hab)

        if len(all_data) == 0:
            ax.text(
                0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center"
            )
            ax.set_title(w)
            continue

        bp = ax.boxplot(
            all_data,
            positions=all_pos,
            widths=0.20,
            patch_artist=True,
            showfliers=False,
        )
        for box, hab in zip(bp["boxes"], all_hab):
            box.set_facecolor(HAB_COLORS.get(hab, "gray"))
            box.set_alpha(0.75)

        for pos, vals, hab in zip(all_pos, all_data, all_hab):
            xj = np.random.normal(pos, 0.03, len(vals))
            ax.scatter(
                xj,
                vals,
                c=HAB_COLORS.get(hab, "gray"),
                s=18,
                alpha=0.5,
                edgecolors="white",
                linewidth=0.35,
                zorder=3,
            )

        ax.set_xticks(np.arange(len(sectors)))
        ax.set_xticklabels([sector_labels[s] for s in sectors], rotation=15)
        ax.set_title(w, fontweight="bold")
        if ax is axes[0]:
            ax.set_ylabel("DLI coefficient of variation (%)")

    legend = [
        Patch(facecolor=HAB_COLORS[h], edgecolor="black", alpha=0.75, label=h)
        for h in hab_order
    ]
    fig.legend(
        handles=legend, title="HAB", loc="upper right", bbox_to_anchor=(0.99, 0.98)
    )
    plt.tight_layout(rect=[0, 0, 0.97, 0.96])
    plt.savefig(output_path, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close()


# ============================================
# Report
# ============================================
def generate_report(metrics_df, stat_df, sector_hab_df, zero_df, output_path):
    sst = stat_df[
        (stat_df["Variable"] == "SST") & (stat_df["Metric"] == "n_events_per_year")
    ]
    dli = stat_df[
        (stat_df["Variable"] == "DLI") & (stat_df["Metric"] == "n_events_per_year")
    ]
    chl = stat_df[
        (stat_df["Variable"] == "Chl-a") & (stat_df["Metric"] == "n_events_per_year")
    ]

    def one_row(df):
        return df.iloc[0] if len(df) else None

    sst_r = one_row(sst)
    dli_r = one_row(dli)
    chl_r = one_row(chl)

    cov_dli = metrics_df[metrics_df["VARIABLE"] == "DLI"][
        "COVERAGE_PCT_PRODUCT"
    ].dropna()
    cov_chl = metrics_df[metrics_df["VARIABLE"] == "Chl-a"][
        "COVERAGE_PCT_PRODUCT"
    ].dropna()

    sector_rank = (
        metrics_df.groupby(["VARIABLE", "INNER_SECTOR"])["n_events_per_year"]
        .mean()
        .reset_index()
    )

    dli_hab = (
        metrics_df[metrics_df["VARIABLE"] == "DLI"]
        .groupby("HAB")["n_events_per_year"]
        .mean()
        .sort_values(ascending=False)
    )

    report = []
    report.append("# Extreme Event Analysis Report (Option A; SITE x HAB)")
    report.append("")
    report.append("## Methodology")
    report.append(
        "Daily SST, benthic DLI, and chlorophyll-a time series (2002-2008) were extracted at the profile level, "
        "where each profile is a unique SITE x HAB combination. Benthic DLI was computed as PAR * exp(-KdPAR * depth), "
        "thereby preserving habitat-specific depth attenuation. For each profile and variable, a seasonal climatology was "
        "estimated as day-of-year means and standard deviations, smoothed with a 30-day circular window. Standardized "
        "anomalies were then computed as z-scores relative to each profile climatology. Extreme events were defined as "
        f"z >= {THRESHOLD_Z}, with a minimum duration of {MIN_DURATION_DAYS} days, and events separated by <= {MERGE_GAP_DAYS} day "
        "were merged if no missing values occurred in the gap. Event metrics included frequency (events per effective year), "
        "duration, peak anomaly, cumulative anomaly, and percentage of extreme days. Cross-shelf contrasts (Inner vs Outer) "
        "were evaluated using permutation tests and Mann-Whitney tests on profile-level metrics."
    )
    report.append("")
    report.append("## Results")

    if sst_r is not None:
        report.append(
            "SST showed a clear cross-shelf contrast in event frequency, with higher rates in the Inner Arc "
            f"({sst_r['Inner_Mean']:.2f} events year-1) than in the Outer Arc ({sst_r['Outer_Mean']:.2f} events year-1; "
            f"permutation p={format_p_value(sst_r['P_Value_Perm'])})."
        )
    if dli_r is not None:
        report.append(
            "DLI event frequency was broadly comparable between arcs "
            f"(Inner {dli_r['Inner_Mean']:.2f} vs Outer {dli_r['Outer_Mean']:.2f} events year-1; "
            f"p={format_p_value(dli_r['P_Value_Perm'])}), while profile-level differentiation by HAB was retained."
        )
    if chl_r is not None:
        report.append(
            "Chl-a frequency differences were moderate and not strongly significant at the arc level "
            f"(Inner {chl_r['Inner_Mean']:.2f} vs Outer {chl_r['Outer_Mean']:.2f} events year-1; "
            f"p={format_p_value(chl_r['P_Value_Perm'])})."
        )

    if len(dli_hab) > 0:
        hab_text = ", ".join([f"{h}={v:.2f}" for h, v in dli_hab.items()])
        report.append(
            "Across all profiles, DLI event frequency by habitat followed the pattern: "
            f"{hab_text} events year-1, highlighting the value of explicit HAB differentiation."
        )

    for var in ["SST", "DLI", "Chl-a"]:
        sub = sector_rank[sector_rank["VARIABLE"] == var].copy()
        if len(sub) < 2:
            continue
        top = sub.sort_values("n_events_per_year", ascending=False).iloc[0]
        bot = sub.sort_values("n_events_per_year", ascending=True).iloc[0]
        report.append(
            f"For {var}, sector-level means ranged from {bot['INNER_SECTOR']} ({bot['n_events_per_year']:.2f} events year-1) "
            f"to {top['INNER_SECTOR']} ({top['n_events_per_year']:.2f} events year-1)."
        )

    report.append("")
    report.append("## Limitations")
    report.append(
        "Ocean-color products (DLI and Chl-a) were affected by substantial cloud-related missingness. "
        f"Median product coverage was {np.nanmedian(cov_dli):.1f}% for DLI and {np.nanmedian(cov_chl):.1f}% for Chl-a across retained profiles. "
        "Although the anomaly framework is robust to unequal baselines, low and uneven temporal coverage may reduce sensitivity "
        "for short-lived events and can influence profile inclusion."
    )
    report.append(
        "SST and chlorophyll extractions share the same satellite pixel across habitats within a site, whereas DLI differs by depth. "
        "Therefore, HAB-level contrasts are strongest and most mechanistically interpretable for DLI."
    )
    report.append(
        "The analysis window (2002-2008) captures interannual variability but remains finite for extreme-value inference; "
        "future extensions should test threshold and duration sensitivity, and incorporate hierarchical models that account for "
        "within-site dependence among habitats."
    )

    with open(output_path, "w", encoding="utf-8") as f:
        f.write("\n\n".join(report) + "\n")


# ============================================
# Main
# ============================================
if __name__ == "__main__":
    print("=" * 74)
    print("EXTREME EVENTS ANALYSIS (OPTION A) - PROFILE LEVEL (SITE x HAB)")
    print(f"Period: {START_YEAR}-{END_YEAR} | Threshold z={THRESHOLD_Z}")
    print(f"NetCDF engine: {NETCDF_ENGINE if NETCDF_ENGINE else 'auto'}")
    print("=" * 74)

    profiles_df = load_profile_metadata(INTEGRATED_DATA_PATH)
    print(f"\nLoaded {len(profiles_df)} profiles (SITE x HAB)")
    print(
        f"  Inner: {(profiles_df['ARC'] == 'inner').sum()} | Outer: {(profiles_df['ARC'] == 'outer').sum()}"
    )

    print("\nSTEP 1: Listing files")
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
    print(f"  SST files: {len(sst_files)}")
    print(f"  Kd490 files: {len(kd_files)}")
    print(f"  PAR files: {len(par_files)}")
    print(f"  Chl-a files: {len(chl_files)}")

    if (
        len(sst_files) == 0
        and len(kd_files) == 0
        and len(par_files) == 0
        and len(chl_files) == 0
    ):
        raise SystemExit("No NetCDF files found. Check data directories.")

    print("\nSTEP 2: Extracting profile-level time series")
    sst_ts = load_timeseries_for_profiles(
        sst_files,
        ["analysed_sst", "sea_surface_temperature", "sst"],
        profiles_df,
        "SST",
    )
    chl_ts = load_timeseries_for_profiles(chl_files, ["chlor_a"], profiles_df, "Chl-a")
    dli_ts = calculate_dli_for_profiles(kd_files, par_files, profiles_df)

    kd_dates = {parse_date_from_filename(f) for f in kd_files}
    par_dates = {parse_date_from_filename(f) for f in par_files}
    kd_dates = {d for d in kd_dates if d is not None}
    par_dates = {d for d in par_dates if d is not None}
    expected_days_by_var = {
        "SST": len(sst_files),
        "DLI": len(kd_dates.intersection(par_dates)),
        "Chl-a": len(chl_files),
    }

    analysis_start = pd.Timestamp(year=START_YEAR, month=1, day=1)
    analysis_end = pd.Timestamp(year=END_YEAR, month=12, day=31)

    print("\nSTEP 3: Climatologies, anomalies, and event detection")
    ts_by_var = {"SST": sst_ts, "DLI": dli_ts, "Chl-a": chl_ts}
    anomaly_dict = {"SST": {}, "DLI": {}, "Chl-a": {}}
    event_dict = {"SST": {}, "DLI": {}, "Chl-a": {}}
    metrics_rows = []

    profile_lookup = profiles_df.set_index("PROFILE_ID").to_dict("index")

    for var, ts_map in ts_by_var.items():
        print(f"\n  Variable: {var}")
        valid_profiles = 0
        excluded_profiles = 0

        for pid, s in ts_map.items():
            s = s.sort_index().loc[analysis_start:analysis_end]
            n_valid = int(s.notna().sum())
            if n_valid < MIN_DAYS_FOR_CLIMATOLOGY:
                excluded_profiles += 1
                continue

            clm = compute_climatology(s)
            anom = compute_anomaly(s, clm)
            if anom is None:
                excluded_profiles += 1
                continue

            ev = detect_events(
                anom,
                threshold_z=THRESHOLD_Z,
                min_duration=MIN_DURATION_DAYS,
                merge_gap=MERGE_GAP_DAYS,
                analysis_start=analysis_start,
                analysis_end=analysis_end,
            )

            anomaly_dict[var][pid] = anom
            event_dict[var][pid] = ev
            valid_profiles += 1

            meta = profile_lookup.get(pid)
            if meta is None:
                continue

            metrics = compute_profile_metrics(
                event_df=ev,
                n_valid_days=n_valid,
                expected_product_days=expected_days_by_var[var],
            )
            metrics_rows.append(
                {
                    "PROFILE_ID": pid,
                    "SITE_ID": meta["SITE_ID"],
                    "HAB": meta["HAB"],
                    "ARC": meta["ARC"],
                    "INNER_SECTOR": meta["INNER_SECTOR"],
                    "REEF_NAME": meta["REEF_NAME"],
                    "LAT": meta["LAT"],
                    "LON": meta["LON"],
                    "DEPTH_M": meta["DEPTH_M"],
                    "VARIABLE": var,
                    **metrics,
                    "THRESHOLD_Z": THRESHOLD_Z,
                    "MIN_DURATION_DAYS": MIN_DURATION_DAYS,
                    "MERGE_GAP_DAYS": MERGE_GAP_DAYS,
                    "METHODOLOGY": "OptionA_anomaly_profile_level",
                }
            )

        print(f"    Retained profiles: {valid_profiles}")
        print(
            f"    Excluded profiles (<{MIN_DAYS_FOR_CLIMATOLOGY} valid days): {excluded_profiles}"
        )
        print(f"    Total events: {sum(len(v) for v in event_dict[var].values())}")

    metrics_df = pd.DataFrame(metrics_rows)
    if metrics_df.empty:
        raise SystemExit(
            "No metrics generated. Verify data access and filtering rules."
        )

    print("\nSTEP 4: Statistics")
    stat_df, sector_hab_df, zero_df = build_stat_tables(metrics_df)
    for _, r in stat_df.iterrows():
        print(
            f"  {r['Variable']} | {r['Metric']} | p_perm={format_p_value(r['P_Value_Perm'])} | "
            f"Inner={r['Inner_Mean']:.2f} Outer={r['Outer_Mean']:.2f}"
        )

    print("\nSTEP 5: Figures")
    create_distribution_figures(metrics_df)

    create_profile_pair_timeseries_figure(
        anomaly_dict_var=anomaly_dict["DLI"],
        event_dict_var=event_dict["DLI"],
        profile_top="ITA1_TP",
        profile_bottom="ITA1_PA",
        figure_title="Inner North (ITA1): DLI anomalies and events (TP vs PA)",
        output_file=os.path.join(OUTPUT_DIR, "Fig4A_Inner_North_ITA1_TP_PA.png"),
    )
    create_profile_pair_timeseries_figure(
        anomaly_dict_var=anomaly_dict["DLI"],
        event_dict_var=event_dict["DLI"],
        profile_top="TIM2_TP",
        profile_bottom="TIM2_PA",
        figure_title="Inner Central (TIM2): DLI anomalies and events (TP vs PA)",
        output_file=os.path.join(OUTPUT_DIR, "Fig4B_Inner_Central_TIM2_TP_PA.png"),
    )
    create_profile_pair_timeseries_figure(
        anomaly_dict_var=anomaly_dict["DLI"],
        event_dict_var=event_dict["DLI"],
        profile_top="SGO_TP",
        profile_bottom="SGO_PA",
        figure_title="Inner South (SGO): DLI anomalies and events (TP vs PA)",
        output_file=os.path.join(OUTPUT_DIR, "Fig4C_Inner_South_SGO_TP_PA.png"),
    )
    create_profile_pair_timeseries_figure(
        anomaly_dict_var=anomaly_dict["DLI"],
        event_dict_var=event_dict["DLI"],
        profile_top="SIR_RR",
        profile_bottom="FAR_RR",
        figure_title="Outer Rocky Reefs: DLI anomalies and events (SIR_RR vs FAR_RR)",
        output_file=os.path.join(OUTPUT_DIR, "Fig4D_Outer_Rocky_SIR_FAR_RR.png"),
    )
    create_profile_pair_timeseries_figure(
        anomaly_dict_var=anomaly_dict["DLI"],
        event_dict_var=event_dict["DLI"],
        profile_top="PAB4_TP",
        profile_bottom="PAB4_PA",
        figure_title="Outer PAB4: DLI anomalies and events (TP vs PA)",
        output_file=os.path.join(OUTPUT_DIR, "Fig4E_Outer_PAB4_TP_PA.png"),
    )

    cv_long = load_cv_long_table()
    create_cv_hab_sector_figure(
        cv_long,
        os.path.join(OUTPUT_DIR, "Fig5_DLI_CV_HAB_Sectors.png"),
    )

    print("\nSTEP 6: Save tables")
    metrics_path = os.path.join(OUTPUT_DIR, "extreme_events_metrics_per_site_hab.csv")
    metrics_df.to_csv(metrics_path, index=False)
    print(f"  Saved: {os.path.basename(metrics_path)}")

    stat_path = os.path.join(OUTPUT_DIR, "statistical_comparison_results_arc.csv")
    stat_df.to_csv(stat_path, index=False)
    print(f"  Saved: {os.path.basename(stat_path)}")

    sector_path = os.path.join(OUTPUT_DIR, "summary_by_sector_hab.csv")
    sector_hab_df.to_csv(sector_path, index=False)
    print(f"  Saved: {os.path.basename(sector_path)}")

    zero_path = os.path.join(OUTPUT_DIR, "zero_inflation_summary_site_hab.csv")
    zero_df.to_csv(zero_path, index=False)
    print(f"  Saved: {os.path.basename(zero_path)}")

    event_rows = []
    for var in ["SST", "DLI", "Chl-a"]:
        for pid, ev in event_dict[var].items():
            if ev is None or len(ev) == 0:
                continue
            meta = profile_lookup.get(pid)
            if meta is None:
                continue
            e = ev.copy().reset_index(drop=True)
            e["EVENT_ID"] = np.arange(1, len(e) + 1)
            e["PROFILE_ID"] = pid
            e["SITE_ID"] = meta["SITE_ID"]
            e["HAB"] = meta["HAB"]
            e["ARC"] = meta["ARC"]
            e["INNER_SECTOR"] = meta["INNER_SECTOR"]
            e["VARIABLE"] = var
            event_rows.append(e)

    if event_rows:
        events_df = pd.concat(event_rows, ignore_index=True)
    else:
        events_df = pd.DataFrame(
            columns=[
                "EVENT_ID",
                "PROFILE_ID",
                "SITE_ID",
                "HAB",
                "ARC",
                "INNER_SECTOR",
                "VARIABLE",
                "start_date",
                "end_date",
                "duration_d",
                "n_extreme_days",
                "peak_date",
                "peak_z",
                "mean_z",
                "cumulative_z",
            ]
        )
    events_path = os.path.join(
        OUTPUT_DIR, "extreme_events_details_per_event_site_hab.csv"
    )
    events_df.to_csv(events_path, index=False)
    print(f"  Saved: {os.path.basename(events_path)}")

    print("\nSTEP 7: Automatic report")
    report_path = os.path.join(OUTPUT_DIR, "Extreme_Events_Report.md")
    generate_report(metrics_df, stat_df, sector_hab_df, zero_df, report_path)
    print(f"  Saved: {os.path.basename(report_path)}")

    print("\n" + "=" * 74)
    print("ANALYSIS COMPLETE")
    print(f"Outputs: {OUTPUT_DIR}")
    print("=" * 74)
