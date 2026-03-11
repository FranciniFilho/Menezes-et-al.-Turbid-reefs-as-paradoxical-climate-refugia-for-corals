"""Improved STL synchrony analysis with diagnostics and uncertainty.

This module is a non-breaking v2 companion for:
`TIME_SERIES_UNIFIED_SYNCHRONY_ANALYSIS.py`.

Key improvements:
1) Explicit daily regularization before interpolation.
2) Transparent coverage diagnostics (original vs interpolated vs missing).
3) STL quality diagnostics (seasonal/trend strength and residual tests).
4) Residual synchrony uncertainty via circular block bootstrap.
5) Robustness checks for STL robust=True vs robust=False.
"""

from __future__ import annotations

import glob
import os
import re
import warnings
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import matplotlib.dates as mdates
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr
from scipy.stats import chi2, kurtosis
from scipy import stats
from statsmodels.tsa.seasonal import STL

warnings.filterwarnings("ignore")


PAIR_DEFS: List[Tuple[str, str, str]] = [
    ("SST", "DLI", "SST_DLI"),
    ("SST", "CHL", "SST_CHL"),
    ("DLI", "CHL", "DLI_CHL"),
]


# ============================================================
# Configuration (edit these for local reruns)
# ============================================================

START_YEAR = 2002
END_YEAR = 2008

SST_DIR = r"K:\remote sensing\CRW_SST_FULL"
MODIS_DIR = r"K:\remote sensing\MODIS_DATA_FULL"

SITES_CSV_PATH = (
    r"C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES"
    r"\############Menezes et al. Mus his distribution and abundance Abrolhos"
    r"\######22.04.23\DATA\sites_list_full_clean.csv"
)

OUTPUT_ROOT_V2 = (
    r"C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES"
    r"\############Menezes et al. Mus his distribution and abundance Abrolhos"
    r"\######FINAL\#######FINAL_RESULTS\Variability_Analysis_V2"
)

STL_DIAGNOSTICS_DIRNAME = "STL_Diagnostics"

SST_PATTERN = "coraltemp_v3.1_*.nc"
KD490_PATTERN = "AQUA_MODIS.*.L3m.DAY.KD.Kd_490.4km.nc"
PAR_PATTERN = "AQUA_MODIS.*.L3m.DAY.PAR.par.4km.nc"
CHL_PATTERN = "AQUA_MODIS.*.L3m.DAY.CHL.chlor_a.4km.nc"

KDPAR_A = 0.0665
KDPAR_B = 0.874
KDPAR_C = 0.00121

STL_PERIOD = 365
STL_MIN_YEARS = 2.0
STL_INTERPOLATION_LIMIT_MAIN = 7
STL_INTERPOLATION_LIMIT_SENSITIVITY = [7, 14, 30]
STL_ROBUST = True
STL_BOOTSTRAP_BLOCK_SIZE = 14
STL_BOOTSTRAP_ITERATIONS = 1000
STL_RANDOM_STATE = 42
SMOOTHING_WINDOW_DAYS = 14

COLOR_INNER = "#CD853F"
COLOR_OUTER = "#4169E1"
COLOR_SST = "#D55E00"
COLOR_DLI = "#0072B2"
COLOR_CHL = "#009E73"

PAIR_LABEL_TO_COLUMNS = {
    "SST_DLI": (
        "Pulse_Corr_SST_DLI",
        "Pulse_CI95_Low_SST_DLI",
        "Pulse_CI95_High_SST_DLI",
    ),
    "SST_CHL": (
        "Pulse_Corr_SST_CHL",
        "Pulse_CI95_Low_SST_CHL",
        "Pulse_CI95_High_SST_CHL",
    ),
    "DLI_CHL": (
        "Pulse_Corr_DLI_CHL",
        "Pulse_CI95_Low_DLI_CHL",
        "Pulse_CI95_High_DLI_CHL",
    ),
}


@dataclass
class InterpolationDiagnostics:
    """Container for interpolation diagnostics for one variable.

    Attributes:
        n_total_days: Number of days in the regularized daily index.
        n_original_days: Number of days with observed values before interpolation.
        n_interpolated_days: Number of days filled by interpolation.
        n_missing_after_interp: Number of days still missing after interpolation.
        pct_original: Percentage of original observations.
        pct_interpolated: Percentage of interpolated values.
        pct_missing_after_interp: Percentage still missing after interpolation.
        max_gap_original_days: Maximum consecutive missing gap before interpolation.
        max_gap_after_interp_days: Maximum consecutive missing gap after interpolation.
    """

    n_total_days: int
    n_original_days: int
    n_interpolated_days: int
    n_missing_after_interp: int
    pct_original: float
    pct_interpolated: float
    pct_missing_after_interp: float
    max_gap_original_days: int
    max_gap_after_interp_days: int


def _max_consecutive_true(mask: pd.Series) -> int:
    """Return the maximum length of consecutive True values in a boolean Series."""
    if mask.empty:
        return 0
    values = mask.astype(int).to_numpy()
    max_run = 0
    run = 0
    for val in values:
        if val == 1:
            run += 1
            if run > max_run:
                max_run = run
        else:
            run = 0
    return int(max_run)


def _to_daily_series(
    series: pd.Series,
    start_date: pd.Timestamp,
    end_date: pd.Timestamp,
) -> pd.Series:
    """Convert a series to daily frequency with explicit missing days.

    Args:
        series: Input time series. Index must be datetime-like.
        start_date: Start date for reindexing.
        end_date: End date for reindexing.

    Returns:
        Daily series indexed from start_date to end_date.
    """
    clean = series.copy()
    clean.index = pd.to_datetime(clean.index)
    clean = clean.groupby(level=0).mean().sort_index()
    daily_index = pd.date_range(start_date, end_date, freq="D")
    return clean.reindex(daily_index)


def interpolate_series_with_diagnostics(
    series_daily: pd.Series,
    max_gap_days: int = 7,
    method: str = "time",
) -> Tuple[pd.Series, pd.Series, InterpolationDiagnostics]:
    """Interpolate a daily series with strict gap diagnostics.

    Args:
        series_daily: Daily series with explicit NaN for missing dates.
        max_gap_days: Maximum consecutive NaN values to fill.
        method: Pandas interpolation method.

    Returns:
        Tuple with:
            - interpolated daily series,
            - boolean mask for interpolated points,
            - diagnostics dataclass.
    """
    observed_mask = series_daily.notna()

    interpolated = series_daily.interpolate(
        method=method,
        limit=max_gap_days,
        limit_area="inside",
    )

    interpolated_mask = series_daily.isna() & interpolated.notna()
    remaining_missing_mask = interpolated.isna()

    n_total = int(len(series_daily))
    n_original = int(observed_mask.sum())
    n_interp = int(interpolated_mask.sum())
    n_missing_after = int(remaining_missing_mask.sum())

    diag = InterpolationDiagnostics(
        n_total_days=n_total,
        n_original_days=n_original,
        n_interpolated_days=n_interp,
        n_missing_after_interp=n_missing_after,
        pct_original=100.0 * n_original / n_total if n_total else np.nan,
        pct_interpolated=100.0 * n_interp / n_total if n_total else np.nan,
        pct_missing_after_interp=100.0 * n_missing_after / n_total
        if n_total
        else np.nan,
        max_gap_original_days=_max_consecutive_true(series_daily.isna()),
        max_gap_after_interp_days=_max_consecutive_true(remaining_missing_mask),
    )
    return interpolated, interpolated_mask, diag


def aggregate_arc_daily_mean(
    ts_dict: Dict[str, pd.Series],
    sites_df: pd.DataFrame,
    arc: str,
    min_sites_per_day: int = 1,
) -> Tuple[pd.Series, pd.Series]:
    """Aggregate site-level data into arc-level daily means.

    Args:
        ts_dict: Dictionary mapping site name to time series.
        sites_df: DataFrame with columns `Site_name` and `Arc`.
        arc: Arc label (e.g., 'inner', 'outer').
        min_sites_per_day: Minimum contributing sites required per day.

    Returns:
        Tuple of:
            - daily arc mean series (irregular index, pre-regularization),
            - number of contributing sites per day.
    """
    arc_sites = sites_df.loc[
        sites_df["Arc"].str.strip().str.lower() == arc.lower(), "Site_name"
    ].tolist()
    site_series = [ts_dict[s] for s in arc_sites if s in ts_dict]
    if not site_series:
        return pd.Series(dtype=float), pd.Series(dtype=float)

    concat = pd.concat(site_series, axis=1)
    n_sites = concat.notna().sum(axis=1)
    arc_mean = concat.mean(axis=1)
    arc_mean.loc[n_sites < min_sites_per_day] = np.nan
    return arc_mean.sort_index(), n_sites.sort_index()


def stl_strength_metrics(
    observed: pd.Series, trend: pd.Series, seasonal: pd.Series, resid: pd.Series
) -> Dict[str, float]:
    """Compute STL component strength metrics (Wang et al.-style).

    Definitions:
        seasonal_strength = max(0, 1 - Var(R) / Var(S + R))
        trend_strength = max(0, 1 - Var(R) / Var(T + R))

    Args:
        observed: Original (interpolated) observed series.
        trend: STL trend component.
        seasonal: STL seasonal component.
        resid: STL residual component.

    Returns:
        Dictionary with strength and variance partition metrics.
    """
    var_obs = float(np.nanvar(observed, ddof=1))
    var_resid = float(np.nanvar(resid, ddof=1))
    var_trend = float(np.nanvar(trend, ddof=1))
    var_seasonal = float(np.nanvar(seasonal, ddof=1))
    var_sr = float(np.nanvar(seasonal + resid, ddof=1))
    var_tr = float(np.nanvar(trend + resid, ddof=1))

    seasonal_strength = max(0.0, 1.0 - var_resid / var_sr) if var_sr > 0 else np.nan
    trend_strength = max(0.0, 1.0 - var_resid / var_tr) if var_tr > 0 else np.nan

    return {
        "seasonal_strength": seasonal_strength,
        "trend_strength": trend_strength,
        "var_observed": var_obs,
        "var_trend": var_trend,
        "var_seasonal": var_seasonal,
        "var_residual": var_resid,
        "residual_variance_share": var_resid / var_obs if var_obs > 0 else np.nan,
    }


def residual_pulse_diagnostics(
    resid: pd.Series, ljungbox_lag: int = 30
) -> Dict[str, float]:
    """Quantify whether residuals look pulse-like or mostly white noise.

    Args:
        resid: STL residual series.
        ljungbox_lag: Lag for Ljung-Box whiteness test.

    Returns:
        Diagnostics including autocorrelation and tail-heaviness.
    """
    clean = resid.dropna()
    if len(clean) < max(40, ljungbox_lag + 5):
        return {
            "resid_acf_lag1": np.nan,
            "resid_kurtosis_excess": np.nan,
            "resid_extreme_frac_abs_z_gt_2": np.nan,
            "ljungbox_pvalue": np.nan,
        }

    z = (clean - clean.mean()) / (clean.std(ddof=1) + 1e-12)

    x = clean.to_numpy(dtype=float)
    x = x - np.mean(x)
    denom = np.dot(x, x)
    if denom <= 0:
        ljungbox_pvalue = np.nan
    else:
        acf_vals = []
        for k in range(1, ljungbox_lag + 1):
            num = np.dot(x[k:], x[:-k])
            acf_vals.append(num / denom)
        n = len(x)
        q_stat = (
            n
            * (n + 2)
            * np.sum(
                [(acf_vals[k - 1] ** 2) / (n - k) for k in range(1, ljungbox_lag + 1)]
            )
        )
        ljungbox_pvalue = float(1.0 - chi2.cdf(q_stat, df=ljungbox_lag))

    return {
        "resid_acf_lag1": float(clean.autocorr(lag=1)),
        "resid_kurtosis_excess": float(
            kurtosis(clean.to_numpy(), fisher=True, bias=False)
        ),
        "resid_extreme_frac_abs_z_gt_2": float((np.abs(z) > 2.0).mean()),
        "ljungbox_pvalue": ljungbox_pvalue,
    }


def circular_block_bootstrap_corr(
    x: pd.Series,
    y: pd.Series,
    block_size: int = 14,
    n_boot: int = 1000,
    random_state: int = 42,
) -> Dict[str, float]:
    """Estimate uncertainty of correlation using circular block bootstrap.

    This preserves local temporal autocorrelation better than iid bootstrap.

    Args:
        x: Series 1.
        y: Series 2.
        block_size: Block length in days.
        n_boot: Number of bootstrap resamples.
        random_state: Random seed.

    Returns:
        Dictionary with observed correlation, bootstrap CI, and sign stability.
    """
    df = pd.DataFrame({"x": x, "y": y}).dropna()
    n = len(df)
    if n < max(60, block_size * 3):
        return {
            "corr": np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "p_two_sided_sign": np.nan,
            "n_effective": int(n),
        }

    x_vals = df["x"].to_numpy()
    y_vals = df["y"].to_numpy()
    obs_corr = float(np.corrcoef(x_vals, y_vals)[0, 1])

    rng = np.random.default_rng(random_state)
    boot_corr = np.empty(n_boot, dtype=float)
    starts_max = n

    for b in range(n_boot):
        idx_blocks: List[int] = []
        while len(idx_blocks) < n:
            start = int(rng.integers(0, starts_max))
            block = [(start + k) % n for k in range(block_size)]
            idx_blocks.extend(block)
        idx = np.array(idx_blocks[:n], dtype=int)
        boot_corr[b] = np.corrcoef(x_vals[idx], y_vals[idx])[0, 1]

    ci_low, ci_high = np.nanpercentile(boot_corr, [2.5, 97.5])
    p_two_sided_sign = 2.0 * min((boot_corr >= 0).mean(), (boot_corr <= 0).mean())

    return {
        "corr": obs_corr,
        "ci_low": float(ci_low),
        "ci_high": float(ci_high),
        "p_two_sided_sign": float(p_two_sided_sign),
        "n_effective": int(n),
    }


def leave_one_year_out_corr(
    x: pd.Series, y: pd.Series, min_points: int = 120
) -> Dict[str, float]:
    """Assess temporal stability of correlation using leave-one-year-out CV.

    Args:
        x: Series 1.
        y: Series 2.
        min_points: Minimum points required in training subset.

    Returns:
        Summary of LOYO correlation variability.
    """
    df = pd.DataFrame({"x": x, "y": y}).dropna()
    if len(df) < min_points:
        return {
            "loyo_n_folds": 0,
            "loyo_corr_mean": np.nan,
            "loyo_corr_sd": np.nan,
            "loyo_corr_min": np.nan,
            "loyo_corr_max": np.nan,
        }

    years = sorted(df.index.year.unique())
    fold_corr: List[float] = []
    for yr in years:
        train = df[df.index.year != yr]
        if len(train) < min_points:
            continue
        fold_corr.append(float(train["x"].corr(train["y"])))

    if not fold_corr:
        return {
            "loyo_n_folds": 0,
            "loyo_corr_mean": np.nan,
            "loyo_corr_sd": np.nan,
            "loyo_corr_min": np.nan,
            "loyo_corr_max": np.nan,
        }

    return {
        "loyo_n_folds": int(len(fold_corr)),
        "loyo_corr_mean": float(np.nanmean(fold_corr)),
        "loyo_corr_sd": float(np.nanstd(fold_corr, ddof=1))
        if len(fold_corr) > 1
        else 0.0,
        "loyo_corr_min": float(np.nanmin(fold_corr)),
        "loyo_corr_max": float(np.nanmax(fold_corr)),
    }


def calculate_pulse_synchrony_stl_v2(
    sst: pd.Series,
    dli: pd.Series,
    chl: pd.Series,
    period: int = 365,
    min_years: float = 2.0,
    interpolation_limit_days: int = 7,
    interpolation_method: str = "time",
    robust: bool = True,
    seasonal: Optional[int] = 13,
    trend: Optional[int] = None,
    low_pass: Optional[int] = None,
    chl_transform: str = "log1p",
    bootstrap_block_size: int = 14,
    bootstrap_iterations: int = 1000,
    random_state: int = 42,
) -> Dict[str, object]:
    """Run improved STL pulse synchrony analysis with diagnostics.

    Args:
        sst: Arc-level SST series.
        dli: Arc-level DLI series.
        chl: Arc-level Chl-a series.
        period: Seasonal period in days.
        min_years: Minimum number of seasonal cycles required.
        interpolation_limit_days: Max consecutive days to interpolate.
        interpolation_method: Pandas interpolation method.
        robust: STL robust option.
        seasonal: STL seasonal smoother length.
        trend: STL trend smoother length.
        low_pass: STL low-pass smoother length.
        chl_transform: Optional transform for CHL ('log1p' or 'none').
        bootstrap_block_size: Block size for bootstrap CI.
        bootstrap_iterations: Number of bootstrap iterations.
        random_state: Random seed.

    Returns:
        Dictionary with coverage diagnostics, STL components, and pairwise synchrony.
    """
    series_map: Dict[str, pd.Series] = {
        "SST": sst,
        "DLI": dli,
        "CHL": chl,
    }

    if any(len(v) == 0 for v in series_map.values()):
        return {
            "status": "failed",
            "reason": "One or more input series are empty.",
        }

    global_start = min(pd.to_datetime(v.index).min() for v in series_map.values())
    global_end = max(pd.to_datetime(v.index).max() for v in series_map.values())

    prepared: Dict[str, Dict[str, object]] = {}
    coverage_rows: List[Dict[str, object]] = []

    for var_name, raw in series_map.items():
        daily = _to_daily_series(raw, global_start, global_end)
        if var_name == "CHL" and chl_transform.lower() == "log1p":
            daily = np.log1p(daily.clip(lower=0))

        filled, interp_mask, diag = interpolate_series_with_diagnostics(
            daily,
            max_gap_days=interpolation_limit_days,
            method=interpolation_method,
        )

        prepared[var_name] = {
            "daily_original": daily,
            "daily_filled": filled,
            "interp_mask": interp_mask,
            "diag": diag,
        }

        coverage_rows.append(
            {
                "Variable": var_name,
                "N_Total_Days": diag.n_total_days,
                "N_Original_Days": diag.n_original_days,
                "N_Interpolated_Days": diag.n_interpolated_days,
                "N_Missing_After_Interp": diag.n_missing_after_interp,
                "Pct_Original": diag.pct_original,
                "Pct_Interpolated": diag.pct_interpolated,
                "Pct_Missing_After_Interp": diag.pct_missing_after_interp,
                "Max_Gap_Original_Days": diag.max_gap_original_days,
                "Max_Gap_After_Interp_Days": diag.max_gap_after_interp_days,
            }
        )

    df_original = pd.DataFrame(
        {k: prepared[k]["daily_original"] for k in ["SST", "DLI", "CHL"]}
    )
    df_filled = pd.DataFrame(
        {k: prepared[k]["daily_filled"] for k in ["SST", "DLI", "CHL"]}
    )

    n_joint_original = int(df_original.dropna().shape[0])
    n_joint_after_interp = int(df_filled.dropna().shape[0])
    min_days_required = int(period * min_years)

    status = "ok" if n_joint_after_interp >= min_days_required else "insufficient_data"
    reason = ""
    if status != "ok":
        reason = (
            f"Joint complete days after interpolation ({n_joint_after_interp}) "
            f"< minimum required ({min_days_required})."
        )

    result: Dict[str, object] = {
        "status": status,
        "reason": reason,
        "coverage_by_variable": pd.DataFrame(coverage_rows),
        "joint_days_original": n_joint_original,
        "joint_days_after_interpolation": n_joint_after_interp,
        "min_days_required": min_days_required,
        "stl_executed": status == "ok",
    }

    if status != "ok":
        return result

    use_df = df_filled.dropna().copy()
    stl_models: Dict[str, object] = {}
    component_rows: List[Dict[str, object]] = []

    for var_name in ["SST", "DLI", "CHL"]:
        fit = STL(
            use_df[var_name],
            period=period,
            robust=robust,
            seasonal=seasonal,
            trend=trend,
            low_pass=low_pass,
        ).fit()
        stl_models[var_name] = fit

        metrics = stl_strength_metrics(
            observed=use_df[var_name],
            trend=fit.trend,
            seasonal=fit.seasonal,
            resid=fit.resid,
        )
        metrics.update(residual_pulse_diagnostics(fit.resid))
        metrics["Variable"] = var_name
        component_rows.append(metrics)

    residual_df = pd.DataFrame(
        {k: stl_models[k].resid for k in ["SST", "DLI", "CHL"]}
    ).dropna()

    pair_rows: List[Dict[str, object]] = []
    for var_a, var_b, pair_label in PAIR_DEFS:
        boot = circular_block_bootstrap_corr(
            residual_df[var_a],
            residual_df[var_b],
            block_size=bootstrap_block_size,
            n_boot=bootstrap_iterations,
            random_state=random_state,
        )
        loyo = leave_one_year_out_corr(residual_df[var_a], residual_df[var_b])

        pair_rows.append(
            {
                "Pair": pair_label,
                "Pulse_Corr": boot["corr"],
                "CI95_Low": boot["ci_low"],
                "CI95_High": boot["ci_high"],
                "Bootstrap_p_two_sided_sign": boot["p_two_sided_sign"],
                "N_Effective": boot["n_effective"],
                "LOYO_N_Folds": loyo["loyo_n_folds"],
                "LOYO_Corr_Mean": loyo["loyo_corr_mean"],
                "LOYO_Corr_SD": loyo["loyo_corr_sd"],
                "LOYO_Corr_Min": loyo["loyo_corr_min"],
                "LOYO_Corr_Max": loyo["loyo_corr_max"],
            }
        )

    robust_alt = not robust
    alt_models = {
        var_name: STL(
            use_df[var_name],
            period=period,
            robust=robust_alt,
            seasonal=seasonal,
            trend=trend,
            low_pass=low_pass,
        ).fit()
        for var_name in ["SST", "DLI", "CHL"]
    }

    robust_rows: List[Dict[str, object]] = []
    for var_a, var_b, pair_label in PAIR_DEFS:
        corr_main = pd.Series(stl_models[var_a].resid).corr(
            pd.Series(stl_models[var_b].resid)
        )
        corr_alt = pd.Series(alt_models[var_a].resid).corr(
            pd.Series(alt_models[var_b].resid)
        )
        robust_rows.append(
            {
                "Pair": pair_label,
                "Corr_robust_main": float(corr_main),
                "Corr_robust_alt": float(corr_alt),
                "Delta_main_minus_alt": float(corr_main - corr_alt),
                "Main_robust_value": robust,
                "Alt_robust_value": robust_alt,
            }
        )

    result.update(
        {
            "residual_df": residual_df,
            "component_diagnostics": pd.DataFrame(component_rows),
            "pairwise_synchrony": pd.DataFrame(pair_rows),
            "robust_sensitivity": pd.DataFrame(robust_rows),
            "stl_models": stl_models,
            "prepared_series": prepared,
        }
    )
    return result


def diagnose_stl_by_arc(
    sst_ts: Dict[str, pd.Series],
    dli_ts: Dict[str, pd.Series],
    chl_ts: Dict[str, pd.Series],
    sites_df: pd.DataFrame,
    arcs: Optional[List[str]] = None,
    period: int = 365,
    interpolation_limit_days: int = 7,
    min_years: float = 2.0,
    robust: bool = True,
    bootstrap_iterations: int = 1000,
) -> Dict[str, object]:
    """Diagnose and run STL synchrony by arc with publication-ready outputs.

    Args:
        sst_ts: Site-level SST series dict.
        dli_ts: Site-level DLI series dict.
        chl_ts: Site-level CHL series dict.
        sites_df: DataFrame with Arc/Site_name metadata.
        arcs: Optional subset of arcs.
        period: STL period.
        interpolation_limit_days: Interpolation limit.
        min_years: Minimum cycles required.
        robust: STL robust argument.
        bootstrap_iterations: Number of bootstrap iterations.

    Returns:
        Dictionary with arc summaries and per-arc detailed analysis objects.
    """
    if arcs is None:
        arcs = sorted(
            sites_df["Arc"].str.strip().str.lower().dropna().unique().tolist()
        )

    status_rows: List[Dict[str, object]] = []
    coverage_rows: List[pd.DataFrame] = []
    pair_rows: List[pd.DataFrame] = []
    component_rows: List[pd.DataFrame] = []
    robust_rows: List[pd.DataFrame] = []
    details: Dict[str, Dict[str, object]] = {}

    for arc in arcs:
        sst_arc, n_sst = aggregate_arc_daily_mean(sst_ts, sites_df, arc)
        dli_arc, n_dli = aggregate_arc_daily_mean(dli_ts, sites_df, arc)
        chl_arc, n_chl = aggregate_arc_daily_mean(chl_ts, sites_df, arc)

        out = calculate_pulse_synchrony_stl_v2(
            sst=sst_arc,
            dli=dli_arc,
            chl=chl_arc,
            period=period,
            min_years=min_years,
            interpolation_limit_days=interpolation_limit_days,
            robust=robust,
            bootstrap_iterations=bootstrap_iterations,
        )

        details[arc] = {
            "analysis": out,
            "n_sites_daily_sst": n_sst,
            "n_sites_daily_dli": n_dli,
            "n_sites_daily_chl": n_chl,
        }

        status_rows.append(
            {
                "Arc": arc,
                "STL_Executed": out.get("stl_executed", False),
                "Status": out.get("status"),
                "Reason": out.get("reason", ""),
                "Joint_Days_Original": out.get("joint_days_original", np.nan),
                "Joint_Days_After_Interp": out.get(
                    "joint_days_after_interpolation", np.nan
                ),
                "Min_Days_Required": out.get("min_days_required", np.nan),
            }
        )

        cov_df = out.get("coverage_by_variable")
        if isinstance(cov_df, pd.DataFrame) and not cov_df.empty:
            cov_df = cov_df.copy()
            cov_df["Arc"] = arc
            coverage_rows.append(cov_df)

        pair_df = out.get("pairwise_synchrony")
        if isinstance(pair_df, pd.DataFrame) and not pair_df.empty:
            pair_df = pair_df.copy()
            pair_df["Arc"] = arc
            pair_rows.append(pair_df)

        comp_df = out.get("component_diagnostics")
        if isinstance(comp_df, pd.DataFrame) and not comp_df.empty:
            comp_df = comp_df.copy()
            comp_df["Arc"] = arc
            component_rows.append(comp_df)

        robust_df = out.get("robust_sensitivity")
        if isinstance(robust_df, pd.DataFrame) and not robust_df.empty:
            robust_df = robust_df.copy()
            robust_df["Arc"] = arc
            robust_rows.append(robust_df)

    return {
        "status_df": pd.DataFrame(status_rows),
        "coverage_df": pd.concat(coverage_rows, ignore_index=True)
        if coverage_rows
        else pd.DataFrame(),
        "pairwise_df": pd.concat(pair_rows, ignore_index=True)
        if pair_rows
        else pd.DataFrame(),
        "component_df": pd.concat(component_rows, ignore_index=True)
        if component_rows
        else pd.DataFrame(),
        "robustness_df": pd.concat(robust_rows, ignore_index=True)
        if robust_rows
        else pd.DataFrame(),
        "details": details,
    }


def plot_stl_components_diagnostic(
    arc_result: Dict[str, object],
    arc_label: str,
    output_path: Optional[str] = None,
) -> Optional[plt.Figure]:
    """Plot observed/trend/seasonal/residual STL components for one arc.

    Args:
        arc_result: Output from `calculate_pulse_synchrony_stl_v2`.
        arc_label: Arc label for title.
        output_path: Optional path to save figure.

    Returns:
        Matplotlib figure, or None if STL did not run.
    """
    if arc_result.get("status") != "ok":
        return None

    stl_models = arc_result["stl_models"]
    prepared = arc_result["prepared_series"]
    variables = ["SST", "DLI", "CHL"]
    rows = ["Observed", "Trend", "Seasonal", "Residual"]

    fig, axes = plt.subplots(4, 3, figsize=(18, 10), sharex=True)
    fig.suptitle(
        f"STL Diagnostic Components - {arc_label}", fontsize=15, fontweight="bold"
    )

    for col, var_name in enumerate(variables):
        model = stl_models[var_name]
        daily_orig = prepared[var_name]["daily_original"]
        daily_fill = prepared[var_name]["daily_filled"]
        interp_mask = prepared[var_name]["interp_mask"]

        ax_obs = axes[0, col]
        ax_obs.plot(
            daily_fill.index,
            daily_fill.values,
            color="#1b9e77",
            linewidth=1.0,
            alpha=0.7,
            label="Filled",
        )
        ax_obs.scatter(
            daily_orig.index[daily_orig.notna()],
            daily_orig[daily_orig.notna()].values,
            s=4,
            color="#2c3e50",
            alpha=0.6,
            label="Original",
        )
        ax_obs.scatter(
            daily_fill.index[interp_mask],
            daily_fill[interp_mask].values,
            s=6,
            color="#e67e22",
            alpha=0.8,
            label="Interpolated",
        )
        ax_obs.set_title(var_name, fontweight="bold")
        ax_obs.grid(alpha=0.2)

        axes[1, col].plot(
            model.trend.index, model.trend.values, color="#4c78a8", linewidth=1.2
        )
        axes[1, col].grid(alpha=0.2)

        axes[2, col].plot(
            model.seasonal.index, model.seasonal.values, color="#59a14f", linewidth=1.0
        )
        axes[2, col].grid(alpha=0.2)

        axes[3, col].plot(
            model.resid.index, model.resid.values, color="#e15759", linewidth=0.9
        )
        axes[3, col].axhline(0, color="black", linewidth=0.8, alpha=0.5)
        axes[3, col].grid(alpha=0.2)

    for r, name in enumerate(rows):
        axes[r, 0].set_ylabel(name)

    handles, labels = axes[0, 0].get_legend_handles_labels()
    if handles:
        fig.legend(handles, labels, loc="upper right", ncol=3, frameon=False)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
    return fig


# ============================================================
# End-to-end v2 workflow (keeps original 5x3 main figure)
# ============================================================


def parse_date_from_filename(filename: str) -> Optional[pd.Timestamp]:
    """Extract YYYYMMDD date from filename."""
    base = os.path.basename(filename)
    match = re.search(r"(\d{8})", base)
    if not match:
        return None
    try:
        return pd.to_datetime(match.group(1), format="%Y%m%d")
    except ValueError:
        return None


def filter_files_by_period(
    files: List[str], start_year: int, end_year: int
) -> List[str]:
    """Filter netCDF files to an inclusive year interval."""
    selected: List[str] = []
    for fpath in files:
        date = parse_date_from_filename(fpath)
        if date is not None and start_year <= date.year <= end_year:
            selected.append(fpath)
    return sorted(selected)


def load_timeseries_for_sites_v2(
    files: List[str],
    var_names: List[str],
    sites_df: pd.DataFrame,
    var_label: str,
) -> Dict[str, pd.Series]:
    """Load one variable time series per site from daily netCDF files."""
    print(f"  Loading {var_label} ...")
    out: Dict[str, Dict[str, List[object]]] = {
        row["Site_name"]: {"times": [], "values": []} for _, row in sites_df.iterrows()
    }
    for idx, fpath in enumerate(files):
        if (idx + 1) % 1000 == 0:
            print(f"    processed {idx + 1}/{len(files)} files")
        date = parse_date_from_filename(fpath)
        if date is None:
            continue
        try:
            with xr.open_dataset(fpath) as ds:
                var_name = next((v for v in var_names if v in ds.data_vars), None)
                if var_name is None:
                    continue
                data = ds[var_name]
                if "time" in data.dims:
                    data = data.isel(time=0)

                for _, row in sites_df.iterrows():
                    site = row["Site_name"]
                    try:
                        value = float(
                            data.sel(
                                lat=row["Latitude"],
                                lon=row["Longitude"],
                                method="nearest",
                            ).values
                        )
                        if not np.isnan(value):
                            out[site]["times"].append(date)
                            out[site]["values"].append(value)
                    except Exception:
                        continue
        except Exception:
            continue

    series_dict: Dict[str, pd.Series] = {}
    for site, payload in out.items():
        if payload["times"]:
            series_dict[site] = pd.Series(
                payload["values"], index=pd.DatetimeIndex(payload["times"])
            ).sort_index()
    return series_dict


def calculate_dli_for_sites_v2(
    kd490_files: List[str],
    par_files: List[str],
    sites_df: pd.DataFrame,
) -> Dict[str, pd.Series]:
    """Compute benthic DLI per site from KD490 and PAR files."""
    print("  Calculating DLI ...")
    par_by_date = {
        parse_date_from_filename(f): f for f in par_files if parse_date_from_filename(f)
    }
    out: Dict[str, Dict[str, List[object]]] = {
        row["Site_name"]: {"times": [], "values": []} for _, row in sites_df.iterrows()
    }

    for idx, kd_file in enumerate(kd490_files):
        if (idx + 1) % 1000 == 0:
            print(f"    processed {idx + 1}/{len(kd490_files)} files")
        date = parse_date_from_filename(kd_file)
        if date is None or date not in par_by_date:
            continue
        par_file = par_by_date[date]

        try:
            with xr.open_dataset(kd_file) as ds_kd, xr.open_dataset(par_file) as ds_par:
                kd_data = ds_kd["Kd_490"]
                par_data = ds_par["par"]
                if "time" in kd_data.dims:
                    kd_data = kd_data.isel(time=0)
                if "time" in par_data.dims:
                    par_data = par_data.isel(time=0)

                for _, row in sites_df.iterrows():
                    site = row["Site_name"]
                    try:
                        kd_val = float(
                            kd_data.sel(
                                lat=row["Latitude"],
                                lon=row["Longitude"],
                                method="nearest",
                            ).values
                        )
                        par_val = float(
                            par_data.sel(
                                lat=row["Latitude"],
                                lon=row["Longitude"],
                                method="nearest",
                            ).values
                        )
                        if np.isnan(kd_val) or np.isnan(par_val) or kd_val <= 0:
                            continue
                        kdpar = KDPAR_A + KDPAR_B * kd_val - KDPAR_C / (kd_val + 1e-9)
                        dli = par_val * np.exp(-kdpar * row["Depth_m"])
                        if not np.isnan(dli) and dli > 0:
                            out[site]["times"].append(date)
                            out[site]["values"].append(dli)
                    except Exception:
                        continue
        except Exception:
            continue

    series_dict: Dict[str, pd.Series] = {}
    for site, payload in out.items():
        if payload["times"]:
            series_dict[site] = pd.Series(
                payload["values"], index=pd.DatetimeIndex(payload["times"])
            ).sort_index()
    return series_dict


def calculate_pearson_synchrony_v2(
    sst_ts: Dict[str, pd.Series],
    dli_ts: Dict[str, pd.Series],
    chl_ts: Dict[str, pd.Series],
    sites_df: pd.DataFrame,
) -> pd.DataFrame:
    """Compute site-level Pearson synchrony for the three variable pairs."""
    rows: List[Dict[str, object]] = []
    for _, row in sites_df.iterrows():
        site = row["Site_name"]
        if site in sst_ts and site in dli_ts and site in chl_ts:
            df = pd.DataFrame(
                {"SST": sst_ts[site], "DLI": dli_ts[site], "CHL": chl_ts[site]}
            ).dropna()
            if len(df) > 30:
                corr = df.corr()
                rows.append(
                    {
                        "Site_name": site,
                        "Arc": row["Arc"],
                        "Corr_SST_DLI": corr.loc["SST", "DLI"],
                        "Corr_SST_CHL": corr.loc["SST", "CHL"],
                        "Corr_DLI_CHL": corr.loc["DLI", "CHL"],
                    }
                )
    return pd.DataFrame(rows)


def calculate_cross_correlation_full_v2(
    series1: pd.Series,
    series2: pd.Series,
    max_lag: int = 60,
) -> Optional[Dict[str, object]]:
    """Compute full lagged cross-correlation and extract optimal lag."""
    df = pd.DataFrame({"s1": series1, "s2": series2}).dropna()
    if len(df) < 60:
        return None

    lags = list(range(-max_lag, max_lag + 1))
    corr = [df["s2"].corr(df["s1"].shift(lag)) for lag in lags]
    corr_arr = np.array(corr, dtype=float)
    if np.isnan(corr_arr).all():
        return None
    best = int(np.nanargmax(np.abs(corr_arr)))
    return {
        "lags": np.array(lags, dtype=int),
        "correlations": corr_arr,
        "optimal_lag": int(lags[best]),
        "optimal_corr": float(corr_arr[best]),
    }


def mann_whitney_test_v2(
    inner_vals: pd.Series, outer_vals: pd.Series
) -> Tuple[float, str]:
    """Mann-Whitney U test with compact p-value string."""
    inner = [float(v) for v in inner_vals if not np.isnan(v)]
    outer = [float(v) for v in outer_vals if not np.isnan(v)]
    if len(inner) < 3 or len(outer) < 3:
        return np.nan, "N/A"
    pvalue = float(stats.mannwhitneyu(inner, outer, alternative="two-sided").pvalue)
    if pvalue < 0.001:
        return pvalue, "p < 0.001"
    if pvalue < 0.05:
        return pvalue, f"p = {pvalue:.3f}"
    return pvalue, f"p = {pvalue:.2f}"


def _pair_metrics_for_arc(
    pairwise_df: pd.DataFrame, arc: str, pair: str
) -> Tuple[float, float, float]:
    """Extract corr and CI bounds for one arc/pair from pairwise diagnostics table."""
    if pairwise_df.empty:
        return np.nan, np.nan, np.nan
    row = pairwise_df[(pairwise_df["Arc"] == arc) & (pairwise_df["Pair"] == pair)]
    if row.empty:
        return np.nan, np.nan, np.nan
    return (
        float(row["Pulse_Corr"].iloc[0]),
        float(row["CI95_Low"].iloc[0]),
        float(row["CI95_High"].iloc[0]),
    )


def build_stl_results_table_for_main_figure(pairwise_df: pd.DataFrame) -> pd.DataFrame:
    """Build arc-level table for panel I, including bootstrap CI for each pair."""
    rows: List[Dict[str, object]] = []
    for arc in ["inner", "outer"]:
        row: Dict[str, object] = {"Arc": arc}
        for pair in ["SST_DLI", "SST_CHL", "DLI_CHL"]:
            corr, low, high = _pair_metrics_for_arc(pairwise_df, arc, pair)
            col_corr, col_low, col_high = PAIR_LABEL_TO_COLUMNS[pair]
            row[col_corr] = corr
            row[col_low] = low
            row[col_high] = high
        rows.append(row)
    return pd.DataFrame(rows)


def create_main_figure_5x3_v2(
    arc_data_plot: Dict[str, Dict[str, pd.Series]],
    pearson_df: pd.DataFrame,
    arc_results: List[Dict[str, object]],
    stl_table: pd.DataFrame,
    output_path: str,
) -> None:
    """Create main 5x3 figure preserving the original layout.

    Rows:
        1-2: Smoothed time series (Inner/Outer).
        3: Pearson boxplots.
        4: Cross-correlation curves.
        5: STL residual synchrony bars with bootstrap CI.
    """
    fig = plt.figure(figsize=(18, 24))
    gs = fig.add_gridspec(
        5, 3, height_ratios=[1.2, 1.2, 1, 1, 0.8], hspace=0.35, wspace=0.25
    )

    for i, arc in enumerate(["inner", "outer"]):
        ax = fig.add_subplot(gs[i, :])
        sst_smooth = (
            arc_data_plot[arc]["sst"].rolling(SMOOTHING_WINDOW_DAYS, center=True).mean()
        )
        dli_smooth = (
            arc_data_plot[arc]["dli"].rolling(SMOOTHING_WINDOW_DAYS, center=True).mean()
        )
        chl_smooth = (
            arc_data_plot[arc]["chl"].rolling(SMOOTHING_WINDOW_DAYS, center=True).mean()
        )

        ax.plot(sst_smooth.index, sst_smooth.values, color=COLOR_SST, lw=2, label="SST")
        ax.set_ylabel("SST (degC)", color=COLOR_SST, fontweight="bold")
        ax.tick_params(axis="y", labelcolor=COLOR_SST)

        ax2 = ax.twinx()
        ax2.plot(
            dli_smooth.index,
            dli_smooth.values,
            color=COLOR_DLI,
            lw=1.5,
            alpha=0.8,
            label="DLI",
        )
        ax2.set_ylabel("DLI (mol m-2 d-1)", color=COLOR_DLI, fontweight="bold")
        ax2.tick_params(axis="y", labelcolor=COLOR_DLI)

        ax3 = ax.twinx()
        ax3.spines["right"].set_position(("outward", 60))
        ax3.plot(
            chl_smooth.index,
            chl_smooth.values,
            color=COLOR_CHL,
            lw=1.5,
            alpha=0.8,
            label="Chl-a",
        )
        ax3.set_ylabel("Chl-a (mg m-3)", color=COLOR_CHL, fontweight="bold")
        ax3.tick_params(axis="y", labelcolor=COLOR_CHL)

        ax.set_title(
            f"Panel {chr(65 + i)}: {arc.capitalize()} Arc Combined Dynamics ({SMOOTHING_WINDOW_DAYS}-day smooth)",
            fontweight="bold",
            loc="left",
            fontsize=11,
        )

    pearson_pairs = [
        ("Corr_SST_DLI", "SST<->DLI"),
        ("Corr_SST_CHL", "SST<->Chl"),
        ("Corr_DLI_CHL", "DLI<->Chl"),
    ]
    for j, (col, label) in enumerate(pearson_pairs):
        ax = fig.add_subplot(gs[2, j])
        if pearson_df.empty or col not in pearson_df.columns:
            ax.text(
                0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes
            )
            continue
        inner_vals = pearson_df[pearson_df["Arc"] == "inner"][col].dropna()
        outer_vals = pearson_df[pearson_df["Arc"] == "outer"][col].dropna()

        bp = ax.boxplot(
            [inner_vals, outer_vals],
            labels=["Inner", "Outer"],
            patch_artist=True,
            widths=0.75,
        )
        bp["boxes"][0].set_facecolor(COLOR_INNER)
        bp["boxes"][1].set_facecolor(COLOR_OUTER)
        for box in bp["boxes"]:
            box.set_alpha(0.7)

        _, p_str = mann_whitney_test_v2(inner_vals, outer_vals)
        ax.text(
            0.88,
            0.82,
            p_str,
            transform=ax.transAxes,
            ha="right",
            va="top",
            fontsize=9,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5),
        )
        ax.axhline(0, color="gray", alpha=0.3, ls="--")
        ax.set_title(
            f"Panel {chr(67 + j)}: {label} (Pearson)",
            fontweight="bold",
            loc="left",
            fontsize=10,
        )
        ax.set_ylabel("Correlation")

    ccf_pairs = [
        ("SST_DLI", "SST -> DLI"),
        ("SST_CHL", "SST -> Chl"),
        ("DLI_CHL", "DLI -> Chl"),
    ]
    for j, (pair_key, label) in enumerate(ccf_pairs):
        ax = fig.add_subplot(gs[3, j])
        for res in [item for item in arc_results if item["Pair"] == pair_key]:
            color = COLOR_INNER if res["Arc"] == "inner" else COLOR_OUTER
            ax.plot(
                res["lags"],
                res["correlations"],
                color=color,
                lw=2,
                label=f"{res['Arc'].capitalize()} (Lag: {res['optimal_lag']}d)",
            )
        ax.axvline(0, color="black", alpha=0.3, ls="--")
        ax.axhline(0, color="gray", alpha=0.3)
        ax.legend(fontsize=8, loc="best")
        ax.set_title(
            f"Panel {chr(70 + j)}: CCF ({label})",
            fontweight="bold",
            loc="left",
            fontsize=10,
        )
        ax.set_xlabel("Lag (days)")
        ax.set_ylabel("Correlation")

    ax_i = fig.add_subplot(gs[4, :])
    if not stl_table.empty:
        x = np.arange(len(stl_table))
        width = 0.25

        sst_dli = stl_table["Pulse_Corr_SST_DLI"].to_numpy(dtype=float)
        sst_chl = stl_table["Pulse_Corr_SST_CHL"].to_numpy(dtype=float)
        dli_chl = stl_table["Pulse_Corr_DLI_CHL"].to_numpy(dtype=float)

        yerr_sst_dli = np.vstack(
            [
                sst_dli - stl_table["Pulse_CI95_Low_SST_DLI"].to_numpy(dtype=float),
                stl_table["Pulse_CI95_High_SST_DLI"].to_numpy(dtype=float) - sst_dli,
            ]
        )
        yerr_sst_chl = np.vstack(
            [
                sst_chl - stl_table["Pulse_CI95_Low_SST_CHL"].to_numpy(dtype=float),
                stl_table["Pulse_CI95_High_SST_CHL"].to_numpy(dtype=float) - sst_chl,
            ]
        )
        yerr_dli_chl = np.vstack(
            [
                dli_chl - stl_table["Pulse_CI95_Low_DLI_CHL"].to_numpy(dtype=float),
                stl_table["Pulse_CI95_High_DLI_CHL"].to_numpy(dtype=float) - dli_chl,
            ]
        )

        bars1 = ax_i.bar(
            x - width,
            sst_dli,
            width,
            yerr=yerr_sst_dli,
            capsize=4,
            label="SST<->DLI",
            color=COLOR_SST,
            alpha=0.85,
        )
        bars2 = ax_i.bar(
            x,
            sst_chl,
            width,
            yerr=yerr_sst_chl,
            capsize=4,
            label="SST<->Chl",
            color="#888888",
            alpha=0.85,
        )
        bars3 = ax_i.bar(
            x + width,
            dli_chl,
            width,
            yerr=yerr_dli_chl,
            capsize=4,
            label="DLI<->Chl",
            color=COLOR_CHL,
            alpha=0.85,
        )

        for bars in [bars1, bars2, bars3]:
            for bar in bars:
                height = bar.get_height()
                ax_i.annotate(
                    f"{height:.2f}",
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3 if height >= 0 else -10),
                    textcoords="offset points",
                    ha="center",
                    va="bottom" if height >= 0 else "top",
                    fontsize=8,
                )

        ax_i.set_xticks(x)
        ax_i.set_xticklabels(stl_table["Arc"].str.capitalize(), fontsize=11)
        ax_i.legend(fontsize=10, loc="upper right")
    ax_i.axhline(0, color="black", alpha=0.5)
    ax_i.set_title(
        "Panel I: STL Residual (Pulse) Correlations with Bootstrap 95% CI",
        fontweight="bold",
        loc="left",
        fontsize=11,
    )
    ax_i.set_ylabel("Residual Correlation", fontsize=10)
    ax_i.set_ylim(-1, 1)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def run_interpolation_sensitivity(
    sst_ts: Dict[str, pd.Series],
    dli_ts: Dict[str, pd.Series],
    chl_ts: Dict[str, pd.Series],
    sites_df: pd.DataFrame,
    interpolation_limits: List[int],
    period: int,
    min_years: float,
    robust: bool,
    bootstrap_iterations: int,
) -> pd.DataFrame:
    """Run STL for multiple interpolation limits and compare residual correlations."""
    rows: List[pd.DataFrame] = []
    for lim in interpolation_limits:
        out = diagnose_stl_by_arc(
            sst_ts=sst_ts,
            dli_ts=dli_ts,
            chl_ts=chl_ts,
            sites_df=sites_df,
            arcs=["inner", "outer"],
            period=period,
            interpolation_limit_days=lim,
            min_years=min_years,
            robust=robust,
            bootstrap_iterations=bootstrap_iterations,
        )
        pair_df = out["pairwise_df"]
        if not pair_df.empty:
            pair_df = pair_df.copy()
            pair_df["Interpolation_Limit_Days"] = lim
            rows.append(pair_df)

    if not rows:
        return pd.DataFrame()

    all_df = pd.concat(rows, ignore_index=True)
    baseline = interpolation_limits[0]
    base_df = all_df[all_df["Interpolation_Limit_Days"] == baseline][
        ["Arc", "Pair", "Pulse_Corr"]
    ].rename(columns={"Pulse_Corr": "Pulse_Corr_Baseline"})
    merged = all_df.merge(base_df, on=["Arc", "Pair"], how="left")
    merged["Delta_vs_Baseline"] = merged["Pulse_Corr"] - merged["Pulse_Corr_Baseline"]
    return merged


def _format_ci(corr: float, low: float, high: float) -> str:
    if np.isnan(corr) or np.isnan(low) or np.isnan(high):
        return "NA"
    return f"{corr:.3f} [{low:.3f}, {high:.3f}]"


def _pulse_evidence_label(corr: float, low: float, high: float) -> str:
    if np.isnan(corr) or np.isnan(low) or np.isnan(high):
        return "indeterminate"
    if low > 0 or high < 0:
        return "evidence_of_synchrony"
    return "weak_or_uncertain"


def build_manuscript_text(
    status_df: pd.DataFrame,
    coverage_df: pd.DataFrame,
    pairwise_df: pd.DataFrame,
    component_df: pd.DataFrame,
) -> Dict[str, str]:
    """Build paper-ready methods, results, and limitations text."""
    if status_df.empty:
        return {
            "methods": "STL analysis could not be completed because no valid arc-level time series were available after preprocessing.",
            "results": "No robust STL-based synchrony estimates were generated.",
            "limitations": "Data availability constraints prevented decomposition and pulse synchrony inference.",
        }

    executed = (
        int(status_df["STL_Executed"].sum())
        if "STL_Executed" in status_df.columns
        else 0
    )
    n_arcs = len(status_df)
    joint_days = ", ".join(
        [
            f"{row['Arc']}: {int(row['Joint_Days_After_Interp'])} days"
            for _, row in status_df.iterrows()
            if not pd.isna(row.get("Joint_Days_After_Interp", np.nan))
        ]
    )

    methods = (
        "Daily SST, DLI, and Chl-a time series were first aggregated by reef arc (inner, outer) as the daily mean across sites. "
        "For each arc and variable, we then projected observations onto a complete daily grid and filled only short internal gaps by time interpolation with a capped maximum gap length. "
        "We decomposed each arc-level series using STL (period = 365 days, robust fitting), and defined pulse synchrony as the Pearson correlation between STL residuals for each variable pair (SST-DLI, SST-Chl-a, DLI-Chl-a). "
        "Uncertainty was estimated with circular block bootstrap (95% CI), preserving local autocorrelation, and sensitivity to interpolation assumptions was evaluated by repeating the analysis across multiple interpolation limits."
    )

    result_lines: List[str] = [
        f"STL ran successfully for {executed}/{n_arcs} arcs.",
        f"Joint complete days after interpolation: {joint_days if joint_days else 'NA'}.",
    ]
    if not pairwise_df.empty:
        for arc in ["inner", "outer"]:
            arc_df = pairwise_df[pairwise_df["Arc"] == arc]
            if arc_df.empty:
                continue
            parts = []
            for pair in ["SST_DLI", "SST_CHL", "DLI_CHL"]:
                row = arc_df[arc_df["Pair"] == pair]
                if row.empty:
                    continue
                corr = float(row["Pulse_Corr"].iloc[0])
                low = float(row["CI95_Low"].iloc[0])
                high = float(row["CI95_High"].iloc[0])
                parts.append(f"{pair}: {_format_ci(corr, low, high)}")
            if parts:
                result_lines.append(
                    f"{arc.capitalize()} arc residual synchrony -> "
                    + "; ".join(parts)
                    + "."
                )

    if not component_df.empty:
        dli_strength = component_df[component_df["Variable"] == "DLI"][
            ["Arc", "seasonal_strength", "residual_variance_share"]
        ]
        if not dli_strength.empty:
            dli_txt = ", ".join(
                [
                    f"{r['Arc']}: seasonal_strength={r['seasonal_strength']:.2f}, residual_share={r['residual_variance_share']:.2f}"
                    for _, r in dli_strength.iterrows()
                    if not pd.isna(r["seasonal_strength"])
                    and not pd.isna(r["residual_variance_share"])
                ]
            )
            if dli_txt:
                result_lines.append(f"DLI decomposition diagnostics -> {dli_txt}.")

    results = " ".join(result_lines)

    limitations = (
        "MODIS-derived DLI and Chl-a contain cloud-related missingness that can be spatially and temporally non-random, especially in offshore sectors. "
        "Therefore, interpolation may attenuate extreme short-lived events or, if too permissive, introduce smooth artificial structure. "
        "Residual correlations should be interpreted jointly with coverage diagnostics, bootstrap intervals, and interpolation sensitivity. "
        "When confidence intervals include zero or change materially across interpolation limits, evidence for pulse synchrony should be considered uncertain."
    )

    return {"methods": methods, "results": results, "limitations": limitations}


def write_stl_analysis_report(
    report_path: str,
    status_df: pd.DataFrame,
    coverage_df: pd.DataFrame,
    pairwise_df: pd.DataFrame,
    component_df: pd.DataFrame,
    robustness_df: pd.DataFrame,
    sensitivity_df: pd.DataFrame,
    manuscript_text: Dict[str, str],
) -> None:
    """Write a transparent text report for STL diagnostics and interpretation."""
    lines: List[str] = []
    lines.append("STL ANALYSIS REPORT - Menezes et al. (Abrolhos)")
    lines.append("=" * 72)
    lines.append("")

    lines.append("1) STL execution status by arc")
    lines.append("-" * 72)
    if status_df.empty:
        lines.append("No status information available.")
    else:
        for _, row in status_df.iterrows():
            lines.append(
                f"Arc={row['Arc']}; executed={row['STL_Executed']}; "
                f"joint_days_after_interp={row['Joint_Days_After_Interp']}; reason={row['Reason']}"
            )
    lines.append("")

    lines.append("2) Data coverage by arc and variable")
    lines.append("-" * 72)
    if coverage_df.empty:
        lines.append("No coverage table available.")
    else:
        for _, row in coverage_df.iterrows():
            lines.append(
                f"Arc={row['Arc']}, Var={row['Variable']}, original={row['Pct_Original']:.1f}%, "
                f"interpolated={row['Pct_Interpolated']:.1f}%, missing_after={row['Pct_Missing_After_Interp']:.1f}%, "
                f"max_gap_original={int(row['Max_Gap_Original_Days'])}"
            )
    lines.append("")

    lines.append("3) Residual pulse synchrony with bootstrap 95% CI")
    lines.append("-" * 72)
    if pairwise_df.empty:
        lines.append("No pairwise STL synchrony available.")
    else:
        for _, row in pairwise_df.iterrows():
            evidence = _pulse_evidence_label(
                row["Pulse_Corr"], row["CI95_Low"], row["CI95_High"]
            )
            lines.append(
                f"Arc={row['Arc']}, Pair={row['Pair']}, corr={row['Pulse_Corr']:.3f}, "
                f"CI95=[{row['CI95_Low']:.3f}, {row['CI95_High']:.3f}], evidence={evidence}"
            )
    lines.append("")

    lines.append("4) Component diagnostics (seasonal strength and residual share)")
    lines.append("-" * 72)
    if component_df.empty:
        lines.append("No component diagnostics available.")
    else:
        for _, row in component_df.iterrows():
            lines.append(
                f"Arc={row['Arc']}, Var={row['Variable']}, seasonal_strength={row['seasonal_strength']:.3f}, "
                f"trend_strength={row['trend_strength']:.3f}, residual_variance_share={row['residual_variance_share']:.3f}, "
                f"ljungbox_p={row['ljungbox_pvalue']:.3f}"
            )
    lines.append("")

    lines.append("5) Robustness: robust=True vs robust=False")
    lines.append("-" * 72)
    if robustness_df.empty:
        lines.append("No robustness table available.")
    else:
        for _, row in robustness_df.iterrows():
            lines.append(
                f"Arc={row['Arc']}, Pair={row['Pair']}, delta_main_minus_alt={row['Delta_main_minus_alt']:.3f}"
            )
    lines.append("")

    lines.append("6) Interpolation sensitivity (limits = 7, 14, 30)")
    lines.append("-" * 72)
    if sensitivity_df.empty:
        lines.append("No interpolation sensitivity table available.")
    else:
        for _, row in sensitivity_df.iterrows():
            lines.append(
                f"Arc={row['Arc']}, Pair={row['Pair']}, limit={int(row['Interpolation_Limit_Days'])}, "
                f"corr={row['Pulse_Corr']:.3f}, delta_vs_baseline={row['Delta_vs_Baseline']:.3f}"
            )
    lines.append("")

    lines.append("7) Interpretation")
    lines.append("-" * 72)
    if pairwise_df.empty:
        lines.append("No conclusion possible: pairwise synchrony was not estimated.")
    else:
        sig_count = int(
            ((pairwise_df["CI95_Low"] > 0) | (pairwise_df["CI95_High"] < 0)).sum()
        )
        total = int(len(pairwise_df))
        if sig_count == 0:
            lines.append(
                "There is no clear evidence of pulse synchrony (all CIs include zero)."
            )
        elif sig_count < total:
            lines.append(
                "There is partial evidence of pulse synchrony (some pairs/arcs exclude zero)."
            )
        else:
            lines.append(
                "There is consistent evidence of pulse synchrony (all pairs/arcs exclude zero)."
            )
    lines.append("")

    lines.append("8) Text for manuscript")
    lines.append("-" * 72)
    lines.append("Methods paragraph:")
    lines.append(manuscript_text["methods"])
    lines.append("")
    lines.append("Results paragraph:")
    lines.append(manuscript_text["results"])
    lines.append("")
    lines.append("Limitations paragraph:")
    lines.append(manuscript_text["limitations"])
    lines.append("")

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))


def run_complete_stl_analysis(
    start_year: int = START_YEAR,
    end_year: int = END_YEAR,
    output_root: str = OUTPUT_ROOT_V2,
    stl_period: int = STL_PERIOD,
    stl_min_years: float = STL_MIN_YEARS,
    interpolation_limit_main: int = STL_INTERPOLATION_LIMIT_MAIN,
    interpolation_limits_sensitivity: Optional[List[int]] = None,
    stl_robust: bool = STL_ROBUST,
    bootstrap_iterations: int = STL_BOOTSTRAP_ITERATIONS,
    max_lag_ccf: int = 60,
) -> Dict[str, object]:
    """Run complete v2 pipeline with figure + diagnostics + report.

    Processing order for STL is explicitly:
        1) Aggregate site series by arc (daily mean).
        2) Build complete daily index.
        3) Interpolate short gaps with configurable limit.
        4) Fit STL and compute residual synchrony.
    """
    if interpolation_limits_sensitivity is None:
        interpolation_limits_sensitivity = STL_INTERPOLATION_LIMIT_SENSITIVITY

    os.makedirs(output_root, exist_ok=True)
    diagnostics_dir = os.path.join(output_root, STL_DIAGNOSTICS_DIRNAME)
    os.makedirs(diagnostics_dir, exist_ok=True)

    print("=" * 70)
    print("UNIFIED SYNCHRONY ANALYSIS V2")
    print("=" * 70)

    sites_df = pd.read_csv(SITES_CSV_PATH, sep=";")
    sites_df["Arc"] = sites_df["Arc"].astype(str).str.strip().str.lower()
    sites_df = sites_df[sites_df["Arc"].isin(["inner", "outer"])].copy()

    print("Finding input files ...")
    sst_files = filter_files_by_period(
        glob.glob(os.path.join(SST_DIR, SST_PATTERN)), start_year, end_year
    )
    kd490_files = filter_files_by_period(
        glob.glob(os.path.join(MODIS_DIR, KD490_PATTERN)), start_year, end_year
    )
    par_files = filter_files_by_period(
        glob.glob(os.path.join(MODIS_DIR, PAR_PATTERN)), start_year, end_year
    )
    chl_files = filter_files_by_period(
        glob.glob(os.path.join(MODIS_DIR, CHL_PATTERN)), start_year, end_year
    )

    print(f"  SST files: {len(sst_files)}")
    print(f"  KD490 files: {len(kd490_files)}")
    print(f"  PAR files: {len(par_files)}")
    print(f"  CHL files: {len(chl_files)}")

    sst_ts = load_timeseries_for_sites_v2(
        sst_files, ["analysed_sst", "sst"], sites_df, "SST"
    )
    chl_ts = load_timeseries_for_sites_v2(chl_files, ["chlor_a"], sites_df, "Chl-a")
    dli_ts = calculate_dli_for_sites_v2(kd490_files, par_files, sites_df)

    pearson_df = calculate_pearson_synchrony_v2(sst_ts, dli_ts, chl_ts, sites_df)
    pearson_df.to_csv(
        os.path.join(output_root, "Pearson_Synchrony_by_Site.csv"), index=False
    )

    arc_data_plot: Dict[str, Dict[str, pd.Series]] = {}
    arc_ccf_results: List[Dict[str, object]] = []
    for arc in ["inner", "outer"]:
        sst_arc_raw, _ = aggregate_arc_daily_mean(sst_ts, sites_df, arc)
        dli_arc_raw, _ = aggregate_arc_daily_mean(dli_ts, sites_df, arc)
        chl_arc_raw, _ = aggregate_arc_daily_mean(chl_ts, sites_df, arc)

        bounds = [
            (series.index.min(), series.index.max())
            for series in [sst_arc_raw, dli_arc_raw, chl_arc_raw]
            if len(series) > 0
        ]
        if not bounds:
            arc_data_plot[arc] = {
                "sst": pd.Series(dtype=float),
                "dli": pd.Series(dtype=float),
                "chl": pd.Series(dtype=float),
            }
            continue

        start = min(item[0] for item in bounds)
        end = max(item[1] for item in bounds)

        sst_daily = _to_daily_series(sst_arc_raw, start, end)
        dli_daily = _to_daily_series(dli_arc_raw, start, end)
        chl_daily = _to_daily_series(chl_arc_raw, start, end)

        sst_fill, _, _ = interpolate_series_with_diagnostics(
            sst_daily, max_gap_days=interpolation_limit_main
        )
        dli_fill, _, _ = interpolate_series_with_diagnostics(
            dli_daily, max_gap_days=interpolation_limit_main
        )
        chl_fill, _, _ = interpolate_series_with_diagnostics(
            chl_daily, max_gap_days=interpolation_limit_main
        )

        arc_data_plot[arc] = {"sst": sst_fill, "dli": dli_fill, "chl": chl_fill}

        ccf_sd = calculate_cross_correlation_full_v2(
            sst_fill, dli_fill, max_lag=max_lag_ccf
        )
        if ccf_sd:
            arc_ccf_results.append({**ccf_sd, "Arc": arc, "Pair": "SST_DLI"})
        ccf_sc = calculate_cross_correlation_full_v2(
            sst_fill, chl_fill, max_lag=max_lag_ccf
        )
        if ccf_sc:
            arc_ccf_results.append({**ccf_sc, "Arc": arc, "Pair": "SST_CHL"})
        ccf_dc = calculate_cross_correlation_full_v2(
            dli_fill, chl_fill, max_lag=max_lag_ccf
        )
        if ccf_dc:
            arc_ccf_results.append({**ccf_dc, "Arc": arc, "Pair": "DLI_CHL"})

    main_diag = diagnose_stl_by_arc(
        sst_ts=sst_ts,
        dli_ts=dli_ts,
        chl_ts=chl_ts,
        sites_df=sites_df,
        arcs=["inner", "outer"],
        period=stl_period,
        interpolation_limit_days=interpolation_limit_main,
        min_years=stl_min_years,
        robust=stl_robust,
        bootstrap_iterations=bootstrap_iterations,
    )

    status_df = main_diag["status_df"]
    coverage_df = main_diag["coverage_df"]
    pairwise_df = main_diag["pairwise_df"]
    component_df = main_diag["component_df"]
    robustness_df = main_diag["robustness_df"]

    status_df.to_csv(
        os.path.join(diagnostics_dir, "STL_Status_by_Arc.csv"), index=False
    )
    coverage_df.to_csv(
        os.path.join(diagnostics_dir, "STL_Coverage_by_Arc_Variable.csv"), index=False
    )
    pairwise_df.to_csv(
        os.path.join(diagnostics_dir, "STL_Pulse_Corr_with_Bootstrap_CI.csv"),
        index=False,
    )
    component_df.to_csv(
        os.path.join(diagnostics_dir, "STL_Component_Diagnostics.csv"), index=False
    )
    robustness_df.to_csv(
        os.path.join(diagnostics_dir, "STL_Robust_vs_Nonrobust.csv"), index=False
    )

    sensitivity_df = run_interpolation_sensitivity(
        sst_ts=sst_ts,
        dli_ts=dli_ts,
        chl_ts=chl_ts,
        sites_df=sites_df,
        interpolation_limits=interpolation_limits_sensitivity,
        period=stl_period,
        min_years=stl_min_years,
        robust=stl_robust,
        bootstrap_iterations=bootstrap_iterations,
    )
    sensitivity_path = os.path.join(
        diagnostics_dir, "STL_Sensitivity_to_Interpolation_Limit.csv"
    )
    sensitivity_df.to_csv(sensitivity_path, index=False)

    for arc in ["inner", "outer"]:
        detail = main_diag["details"].get(arc, {})
        analysis = detail.get("analysis", {})
        fig_name = f"Fig_STL_Components_{arc.capitalize()}.png"
        plot_stl_components_diagnostic(
            analysis,
            arc_label=arc.capitalize(),
            output_path=os.path.join(diagnostics_dir, fig_name),
        )

    stl_table = build_stl_results_table_for_main_figure(pairwise_df)
    create_main_figure_5x3_v2(
        arc_data_plot=arc_data_plot,
        pearson_df=pearson_df,
        arc_results=arc_ccf_results,
        stl_table=stl_table,
        output_path=os.path.join(output_root, "Fig_Main_Temporal_Synchrony_v2.png"),
    )

    pd.DataFrame(arc_ccf_results).to_csv(
        os.path.join(output_root, "Cross_Correlation_Results_v2.csv"), index=False
    )
    stl_table.to_csv(
        os.path.join(output_root, "STL_Pulse_Synchrony_MainFigure_v2.csv"), index=False
    )

    manuscript_text = build_manuscript_text(
        status_df, coverage_df, pairwise_df, component_df
    )
    report_path = os.path.join(diagnostics_dir, "STL_Analysis_Report.txt")
    write_stl_analysis_report(
        report_path=report_path,
        status_df=status_df,
        coverage_df=coverage_df,
        pairwise_df=pairwise_df,
        component_df=component_df,
        robustness_df=robustness_df,
        sensitivity_df=sensitivity_df,
        manuscript_text=manuscript_text,
    )

    paper_text_path = os.path.join(diagnostics_dir, "STL_Text_for_Paper.txt")
    with open(paper_text_path, "w", encoding="utf-8") as f:
        f.write("Methods\n")
        f.write(manuscript_text["methods"] + "\n\n")
        f.write("Results\n")
        f.write(manuscript_text["results"] + "\n\n")
        f.write("Limitations\n")
        f.write(manuscript_text["limitations"] + "\n")

    print("=" * 70)
    print("V2 STL ANALYSIS COMPLETE")
    print(f"Main outputs: {output_root}")
    print(f"Diagnostics: {diagnostics_dir}")
    print("=" * 70)

    return {
        "output_root": output_root,
        "diagnostics_dir": diagnostics_dir,
        "status_df": status_df,
        "coverage_df": coverage_df,
        "pairwise_df": pairwise_df,
        "component_df": component_df,
        "robustness_df": robustness_df,
        "sensitivity_df": sensitivity_df,
        "manuscript_text": manuscript_text,
    }


if __name__ == "__main__":
    run_complete_stl_analysis()
