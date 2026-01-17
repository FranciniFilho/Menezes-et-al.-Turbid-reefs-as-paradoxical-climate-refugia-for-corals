### TIME_SERIES_UNIFIED_SYNCHRONY_ANALYSIS_GLM.py ###
# Temporal Synchrony Analysis: Environmental Dynamics in Abrolhos Reefs
#
# This script analyzes temporal synchrony between environmental variables
# (SST, DLI, Chl-a) to test the hypothesis that Inner Arc reefs exhibit
# stronger coupling between light attenuation, temperature, and chlorophyll.
#
# Hypothesis: ↓ DLI → ↓ SST (after lag) → ↑ Chl-a (nutrient pulses)
#
# Features:
# 1. Pearson Synchrony (correlation between variables by site)
# 2. Cross-Correlation with Lag (temporal coupling analysis)
# 3. STL Decomposition (pulse/residual synchrony)
# 4. Spectral Coherence (frequency-domain synchrony)
#
# All analyses include Mann-Whitney U tests comparing Inner vs Outer arcs.

import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import glob
import re
import logging
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from scipy import stats
from scipy.signal import lombscargle, coherence
from statsmodels.tsa.seasonal import STL
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Detect available NetCDF engine
def get_netcdf_engine():
    try:
        import netCDF4
        return 'netcdf4'
    except ImportError:
        pass
    try:
        import h5netcdf
        return 'h5netcdf'
    except ImportError:
        pass
    try:
        import scipy.io.netcdf
        return 'scipy'
    except ImportError:
        pass
    return None

NETCDF_ENGINE = get_netcdf_engine()
if NETCDF_ENGINE:
    print(f"Using NetCDF engine: {NETCDF_ENGINE}")
else:
    print("WARNING: No NetCDF engine found. Will attempt auto-detection.")

# ============================================
# Configuration
# ============================================
logging.basicConfig(
    filename='synchrony_analysis_errors.log',
    filemode='w',
    level=logging.ERROR,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# --- Time Period ---
START_YEAR = 2002
END_YEAR = 2008

# --- Data Directories ---
SST_DIR = r'K:\remote sensing\CRW_SST_FULL'
MODIS_DIR = r'K:\remote sensing\MODIS_DATA_FULL'

# --- Sites CSV ---
SITES_CSV_PATH = r'C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######22.04.23\DATA\sites_list_full_clean.csv'

# --- Output Directory ---
OUTPUT_DIR = r'C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_RESULTS\Variability_Analysis_GLM'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- File Patterns ---
SST_PATTERN = 'coraltemp_v3.1_*.nc'
KD490_PATTERN = 'AQUA_MODIS.*.L3m.DAY.KD.Kd_490.4km.nc'
PAR_PATTERN = 'AQUA_MODIS.*.L3m.DAY.PAR.par.4km.nc'
CHL_PATTERN = 'AQUA_MODIS.*.L3m.DAY.CHL.chlor_a.4km.nc'

# --- DLI Calculation Constants (Gattuso) ---
KDPAR_A, KDPAR_B, KDPAR_C = 0.0665, 0.874, 0.00121

# --- Color Palette (consistent with variability script) ---
COLOR_INNER = '#CD853F'  # Light brown (Peru)
COLOR_OUTER = '#4169E1'  # Royal Blue
COLORS_ARC = {'inner': COLOR_INNER, 'outer': COLOR_OUTER}
COLOR_SST = '#D55E00'    # Red-orange
COLOR_DLI = '#0072B2'    # Dark blue
COLOR_CHL = '#009E73'    # Green-teal

# --- Publication Style Settings ---
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 11,
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.spines.top': False,
    'axes.spines.right': False
})

# ============================================
# Load Sites Data
# ============================================
print("=" * 60)
print("TEMPORAL SYNCHRONY ANALYSIS: Inner Arc vs Outer Arc")
print(f"Period: {START_YEAR} - {END_YEAR}")
print("=" * 60)

try:
    sites_df = pd.read_csv(SITES_CSV_PATH, sep=';')
    sites_df.columns = [col.strip() for col in sites_df.columns]

    # Verify required columns
    required_cols = {'Site_name', 'Latitude', 'Longitude', 'Depth_m', 'Arc'}
    if not required_cols.issubset(sites_df.columns):
        missing = required_cols - set(sites_df.columns)
        raise ValueError(f"Missing columns: {missing}. Available: {list(sites_df.columns)}")

    # Clean Arc column and standardize
    sites_df['Arc'] = sites_df['Arc'].str.strip().str.lower()
    sites_df = sites_df[sites_df['Arc'].isin(['inner', 'outer'])].copy()

    print(f"\nOK Loaded {len(sites_df)} sites")
    print(f"  Inner Arc: {(sites_df['Arc'] == 'inner').sum()} sites")
    print(f"  Outer Arc: {(sites_df['Arc'] == 'outer').sum()} sites")

except FileNotFoundError:
    print(f"ERROR: Sites file not found: {SITES_CSV_PATH}")
    exit(1)
except Exception as e:
    print(f"ERROR loading sites: {e}")
    exit(1)

# ============================================
# Utility Functions (from variability script)
# ============================================

def parse_date_from_filename(filename):
    """Extract date from filename."""
    basename = os.path.basename(filename)
    match = re.search(r'(\d{8})', basename)
    if match:
        try:
            return pd.to_datetime(match.group(1), format='%Y%m%d')
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

def load_timeseries_for_sites(files, var_names, sites_df, var_label):
    """Load time series data for all sites from a list of files."""
    print(f"\n  Loading {var_label}...")

    site_data = {row['Site_name']: {'times': [], 'values': []}
                 for _, row in sites_df.iterrows()}

    total_files = len(files)
    success_count = 0
    error_count = 0

    for i, f in enumerate(files):
        if (i + 1) % 500 == 0:
            print(f"    Processing file {i+1}/{total_files}...")

        date = parse_date_from_filename(f)
        if date is None:
            continue

        try:
            open_kwargs = {}
            if NETCDF_ENGINE:
                open_kwargs['engine'] = NETCDF_ENGINE

            with xr.open_dataset(f, **open_kwargs) as ds:
                var_name = next((v for v in var_names if v in ds.data_vars), None)
                if var_name is None:
                    continue

                data = ds[var_name]

                if 'time' in data.dims:
                    data = data.isel(time=0)

                for _, row in sites_df.iterrows():
                    site_name = row['Site_name']
                    lat, lon = row['Latitude'], row['Longitude']

                    try:
                        value = float(data.sel(lat=lat, lon=lon, method='nearest').values)
                        if not np.isnan(value):
                            site_data[site_name]['times'].append(date)
                            site_data[site_name]['values'].append(value)
                    except Exception as e:
                        continue

                success_count += 1

        except Exception as e:
            error_count += 1
            if error_count <= 5:
                logging.error(f"Error reading {f}: {e}")
            continue

    print(f"    Files processed successfully: {success_count}/{total_files}")

    result = {}
    for site_name, data in site_data.items():
        if len(data['times']) > 0:
            result[site_name] = pd.Series(data['values'], index=pd.DatetimeIndex(data['times'])).sort_index()

    print(f"    OK Loaded data for {len(result)} sites")
    return result

def calculate_dli_for_sites(kd490_files, par_files, sites_df):
    """Calculate benthic DLI time series for each site."""
    print("\n  Calculating DLI...")

    par_by_date = {parse_date_from_filename(f): f for f in par_files if parse_date_from_filename(f)}

    site_data = {row['Site_name']: {'times': [], 'values': []}
                 for _, row in sites_df.iterrows()}

    total_files = len(kd490_files)
    success_count = 0
    error_count = 0

    open_kwargs = {}
    if NETCDF_ENGINE:
        open_kwargs['engine'] = NETCDF_ENGINE

    for i, f_kd in enumerate(kd490_files):
        if (i + 1) % 500 == 0:
            print(f"    Processing file {i+1}/{total_files}...")

        date = parse_date_from_filename(f_kd)
        if date is None or date not in par_by_date:
            continue

        f_par = par_by_date[date]

        try:
            with xr.open_dataset(f_kd, **open_kwargs) as ds_kd, xr.open_dataset(f_par, **open_kwargs) as ds_par:
                kd490_data = ds_kd['Kd_490']
                par_data = ds_par['par']

                if 'time' in kd490_data.dims:
                    kd490_data = kd490_data.isel(time=0)
                if 'time' in par_data.dims:
                    par_data = par_data.isel(time=0)

                for _, row in sites_df.iterrows():
                    site_name = row['Site_name']
                    lat, lon, depth = row['Latitude'], row['Longitude'], row['Depth_m']

                    try:
                        kd490_val = float(kd490_data.sel(lat=lat, lon=lon, method='nearest').values)
                        par_val = float(par_data.sel(lat=lat, lon=lon, method='nearest').values)

                        if np.isnan(kd490_val) or np.isnan(par_val) or kd490_val <= 0:
                            continue

                        kdpar = KDPAR_A + KDPAR_B * kd490_val - KDPAR_C / (kd490_val + 1e-9)
                        dli = par_val * np.exp(-kdpar * depth)

                        if not np.isnan(dli) and dli > 0:
                            site_data[site_name]['times'].append(date)
                            site_data[site_name]['values'].append(dli)
                    except Exception:
                        continue

                success_count += 1

        except Exception as e:
            error_count += 1
            if error_count <= 5:
                logging.error(f"Error calculating DLI: {e}")
            continue

    print(f"    Files processed successfully: {success_count}/{total_files}")

    result = {}
    for site_name, data in site_data.items():
        if len(data['times']) > 0:
            result[site_name] = pd.Series(data['values'], index=pd.DatetimeIndex(data['times'])).sort_index()

    print(f"    OK Calculated DLI for {len(result)} sites")
    return result

# ============================================
# Data Aggregation Functions
# ============================================

def aggregate_timeseries_by_arc(timeseries_dict, sites_df):
    """
    Aggregate site-level time series into arc-level (Inner/Outer) series.
    Returns a dict with 'inner' and 'outer' keys containing daily mean Series.
    """
    arc_series = {'inner': None, 'outer': None}

    for arc in ['inner', 'outer']:
        arc_sites = sites_df[sites_df['Arc'] == arc]['Site_name'].tolist()

        # Collect all data points for this arc
        all_data = []
        for site in arc_sites:
            if site in timeseries_dict:
                ts = timeseries_dict[site]
                for date, value in ts.items():
                    all_data.append({'date': date, 'value': value})

        if len(all_data) == 0:
            continue

        # Create DataFrame and aggregate by date (mean across sites)
        df = pd.DataFrame(all_data)
        agg = df.groupby('date')['value'].mean()

        arc_series[arc] = agg

    return arc_series

def smooth_timeseries(series, window=14):
    """Apply rolling mean smoothing to a time series."""
    if series is None or len(series) < window:
        return series
    return series.rolling(window=window, min_periods=1, center=True).mean()

def interpolate_gaps(series, max_gap=7):
    """Interpolate small gaps in time series."""
    if series is None:
        return series
    return series.interpolate(method='time', limit=max_gap)

# ============================================
# Synchrony Analysis Functions
# ============================================

def calculate_pearson_synchrony(sst_ts, dli_ts, chl_ts, sites_df):
    """
    Calculate Pearson correlation between variables for each site.
    Returns DataFrame with correlations by site and variable pair.
    """
    results = []

    for _, row in sites_df.iterrows():
        site_name = row['Site_name']
        arc = row['Arc']

        if site_name not in sst_ts or site_name not in dli_ts or site_name not in chl_ts:
            continue

        # Get series for this site
        sst = sst_ts[site_name]
        dli = dli_ts[site_name]
        chl = chl_ts[site_name]

        # Align series (common dates only)
        df = pd.DataFrame({'SST': sst, 'DLI': dli, 'CHL': chl}).dropna()

        if len(df) < 30:
            continue

        # Calculate correlations
        corr_sst_dli = df['SST'].corr(df['DLI'])
        corr_sst_chl = df['SST'].corr(df['CHL'])
        corr_dli_chl = df['DLI'].corr(df['CHL'])

        results.append({
            'Site_name': site_name,
            'Arc': arc,
            'Corr_SST_DLI': corr_sst_dli,
            'Corr_SST_CHL': corr_sst_chl,
            'Corr_DLI_CHL': corr_dli_chl
        })

    return pd.DataFrame(results)

def calculate_cross_correlation_full(series1, series2, max_lag=60):
    """
    Cross-correlation lag analysis between two variables.

    Returns:
    - lags: array of lags tested (-max_lag to +max_lag)
    - correlations: correlation at each lag
    - optimal_lag: lag with maximum absolute correlation
    - optimal_corr: correlation at optimal lag
    """
    df = pd.DataFrame({'s1': series1, 's2': series2}).dropna()

    if len(df) < 60:
        return None

    lags = range(-max_lag, max_lag + 1)
    correlations = []

    for lag in lags:
        corr = df['s2'].corr(df['s1'].shift(lag))
        correlations.append(corr)

    correlations = np.array(correlations)
    optimal_idx = np.nanargmax(np.abs(correlations))

    return {
        'lags': np.array(list(lags)),
        'correlations': correlations,
        'optimal_lag': list(lags)[optimal_idx],
        'optimal_corr': correlations[optimal_idx]
    }

def calculate_cross_correlation_by_arc(sst_arc, dli_arc, chl_arc, max_lag=60):
    """
    Calculate cross-correlation for each arc between variable pairs.
    Returns dict with results for each arc and variable pair.
    """
    results = {'inner': {}, 'outer': {}}

    for arc in ['inner', 'outer']:
        arc_data = {'SST': sst_arc[arc], 'DLI': dli_arc[arc], 'CHL': chl_arc[arc]}

        for var1, var2 in [('SST', 'DLI'), ('SST', 'CHL'), ('DLI', 'CHL')]:
            if arc_data[var1] is None or arc_data[var2] is None:
                results[arc][f'{var1}_{var2}'] = None
                continue

            # Interpolate gaps before CCF
            s1_interp = interpolate_gaps(arc_data[var1])
            s2_interp = interpolate_gaps(arc_data[var2])

            ccf_result = calculate_cross_correlation_full(s1_interp, s2_interp, max_lag)
            results[arc][f'{var1}_{var2}'] = ccf_result

    return results

def calculate_pulse_synchrony_stl(sst, dli, chl, period=365):
    """
    Decompose series using STL and calculate correlation of RESIDUALS (pulses).
    STL is more robust than seasonal_decompose for missing data and outliers.
    """
    # Interpolate gaps (max 7 days)
    sst_interp = interpolate_gaps(sst, max_gap=7)
    dli_interp = interpolate_gaps(dli, max_gap=7)
    chl_interp = interpolate_gaps(chl, max_gap=7)

    # Align series
    df = pd.DataFrame({'SST': sst_interp, 'DLI': dli_interp, 'CHL': chl_interp}).dropna()

    if len(df) < period * 2:
        return None

    try:
        # STL Decomposition
        stl_sst = STL(df['SST'], period=period, robust=True).fit()
        stl_dli = STL(df['DLI'], period=period, robust=True).fit()
        stl_chl = STL(df['CHL'], period=period, robust=True).fit()

        return {
            'Pulse_Corr_SST_DLI': stl_sst.resid.corr(stl_dli.resid),
            'Pulse_Corr_SST_CHL': stl_sst.resid.corr(stl_chl.resid),
            'Pulse_Corr_DLI_CHL': stl_dli.resid.corr(stl_chl.resid),
            'Resid_SST': stl_sst.resid,
            'Resid_DLI': stl_dli.resid,
            'Resid_CHL': stl_chl.resid
        }
    except Exception as e:
        logging.error(f"STL decomposition error: {e}")
        return None

def calculate_pulse_synchrony_by_arc(sst_arc, dli_arc, chl_arc):
    """
    Calculate pulse synchrony (STL residual correlations) for each arc.
    """
    results = {'inner': None, 'outer': None}

    for arc in ['inner', 'outer']:
        if sst_arc[arc] is None or dli_arc[arc] is None or chl_arc[arc] is None:
            results[arc] = None
            continue

        results[arc] = calculate_pulse_synchrony_stl(
            sst_arc[arc], dli_arc[arc], chl_arc[arc]
        )

    return results

def calculate_spectral_coherence_pair(series1, series2, fs=1.0):
    """
    Calculate spectral coherence between two series.
    Returns coherence values and corresponding frequencies.
    """
    # Interpolate and align
    s1 = interpolate_gaps(series1)
    s2 = interpolate_gaps(series2)

    df = pd.DataFrame({'s1': s1, 's2': s2}).dropna()

    if len(df) < 256:  # Minimum length for coherence
        return None

    try:
        # Create regular time series (fill missing dates)
        full_range = pd.date_range(df.index.min(), df.index.max(), freq='D')
        df_regular = df.reindex(full_range).interpolate(method='linear', limit=7)

        # Calculate coherence
        freqs, coh = coherence(df_regular['s1'].values, df_regular['s2'].values, fs=fs, nperseg=256)

        # Extract coherence at specific periods
        periods = 1.0 / freqs
        coh_7d = coh[(periods >= 6) & (periods <= 8)].mean() if len(coh) > 0 else np.nan
        coh_30d = coh[(periods >= 25) & (periods <= 35)].mean() if len(coh) > 0 else np.nan

        return {
            'frequencies': freqs,
            'coherence': coh,
            'coh_7d': coh_7d,
            'coh_30d': coh_30d
        }
    except Exception as e:
        logging.error(f"Spectral coherence error: {e}")
        return None

def calculate_spectral_coherence_by_arc(sst_arc, dli_arc, chl_arc):
    """
    Calculate spectral coherence for each arc between variable pairs.
    """
    results = {'inner': {}, 'outer': {}}

    for arc in ['inner', 'outer']:
        arc_data = {'SST': sst_arc[arc], 'DLI': dli_arc[arc], 'CHL': chl_arc[arc]}

        for var1, var2 in [('SST', 'DLI'), ('SST', 'CHL'), ('DLI', 'CHL')]:
            if arc_data[var1] is None or arc_data[var2] is None:
                results[arc][f'{var1}_{var2}'] = None
                continue

            coh_result = calculate_spectral_coherence_pair(arc_data[var1], arc_data[var2])
            results[arc][f'{var1}_{var2}'] = coh_result

    return results

# ============================================
# Statistical Testing Functions
# ============================================

def mann_whitney_test(inner_values, outer_values):
    """Perform Mann-Whitney U test and return formatted p-value."""
    inner_clean = [v for v in inner_values if not np.isnan(v)]
    outer_clean = [v for v in outer_values if not np.isnan(v)]

    if len(inner_clean) < 3 or len(outer_clean) < 3:
        return np.nan, "N/A"

    try:
        stat, p_value = stats.mannwhitneyu(inner_clean, outer_clean, alternative='two-sided')

        if p_value < 0.001:
            p_str = "p < 0.001"
        elif p_value < 0.01:
            p_str = f"p = {p_value:.3f}"
        else:
            p_str = f"p = {p_value:.2f}"

        return p_value, p_str
    except Exception:
        return np.nan, "N/A"

def format_significance(p_value):
    """Format significance stars."""
    if np.isnan(p_value):
        return ""
    if p_value < 0.001:
        return "***"
    elif p_value < 0.01:
        return "**"
    elif p_value < 0.05:
        return "*"
    return "ns"

# ============================================
# Visualization Functions
# ============================================

def create_main_synchrony_figure(sst_arc, dli_arc, chl_arc, pearson_df, ccf_results, stl_results,
                                 coh_results, output_path):
    """
    Create the main 6-panel figure following the plan's layout:
    - Panel A: Inner Arc Time Series (3 Y-axes, full width)
    - Panel B: Outer Arc Time Series (3 Y-axes, full width)
    - Panel C: Pearson Synchrony (boxplot)
    - Panel D: Cross-Correlation (LAG)
    - Panel E: Spectral Coherence
    - Panel F: Pulse Synchrony (STL Residuals)
    """
    print("\n  Creating Main Synchrony Figure...")

    fig = plt.figure(figsize=(16, 18))
    gs = fig.add_gridspec(4, 2,
                          height_ratios=[1.5, 1.5, 1, 1],
                          hspace=0.25, wspace=0.30)

    import matplotlib.dates as mdates

    # ========== Panel A: Inner Arc Time Series ==========
    ax_inner = fig.add_subplot(gs[0, :])

    sst_in = smooth_timeseries(sst_arc['inner'], window=14)
    dli_in = smooth_timeseries(dli_arc['inner'], window=14)
    chl_in = smooth_timeseries(chl_arc['inner'], window=14)

    if sst_in is not None and len(sst_in) > 0:
        ax_inner.plot(sst_in.index, sst_in.values, color=COLOR_SST, linewidth=1.5, label='SST', alpha=0.9)

    ax_inner.set_ylabel('SST (°C)', color=COLOR_SST, fontsize=11, fontweight='bold')
    ax_inner.tick_params(axis='y', labelcolor=COLOR_SST)

    ax_inner2 = ax_inner.twinx()
    if dli_in is not None and len(dli_in) > 0:
        ax_inner2.plot(dli_in.index, dli_in.values, color=COLOR_DLI, linewidth=1.5, label='DLI', alpha=0.9)
    ax_inner2.set_ylabel('DLI (mol m⁻² d⁻¹)', color=COLOR_DLI, fontsize=11, fontweight='bold')
    ax_inner2.tick_params(axis='y', labelcolor=COLOR_DLI)

    ax_inner3 = ax_inner.twinx()
    ax_inner3.spines['right'].set_position(('outward', 60))
    if chl_in is not None and len(chl_in) > 0:
        ax_inner3.plot(chl_in.index, chl_in.values, color=COLOR_CHL, linewidth=1.5, label='Chl-a', alpha=0.9)
    ax_inner3.set_ylabel('Chl-a (mg m⁻³)', color=COLOR_CHL, fontsize=11, fontweight='bold')
    ax_inner3.tick_params(axis='y', labelcolor=COLOR_CHL)

    ax_inner.set_title('Panel A: Inner Arc - Environmental Time Series (14-day smooth)',
                       fontsize=13, fontweight='bold', loc='left')
    ax_inner.grid(True, alpha=0.3)

    # ========== Panel B: Outer Arc Time Series ==========
    ax_outer = fig.add_subplot(gs[1, :])

    sst_out = smooth_timeseries(sst_arc['outer'], window=14)
    dli_out = smooth_timeseries(dli_arc['outer'], window=14)
    chl_out = smooth_timeseries(chl_arc['outer'], window=14)

    if sst_out is not None and len(sst_out) > 0:
        ax_outer.plot(sst_out.index, sst_out.values, color=COLOR_SST, linewidth=1.5, label='SST', alpha=0.9)

    ax_outer.set_ylabel('SST (°C)', color=COLOR_SST, fontsize=11, fontweight='bold')
    ax_outer.tick_params(axis='y', labelcolor=COLOR_SST)

    ax_outer2 = ax_outer.twinx()
    if dli_out is not None and len(dli_out) > 0:
        ax_outer2.plot(dli_out.index, dli_out.values, color=COLOR_DLI, linewidth=1.5, label='DLI', alpha=0.9)
    ax_outer2.set_ylabel('DLI (mol m⁻² d⁻¹)', color=COLOR_DLI, fontsize=11, fontweight='bold')
    ax_outer2.tick_params(axis='y', labelcolor=COLOR_DLI)

    ax_outer3 = ax_outer.twinx()
    ax_outer3.spines['right'].set_position(('outward', 60))
    if chl_out is not None and len(chl_out) > 0:
        ax_outer3.plot(chl_out.index, chl_out.values, color=COLOR_CHL, linewidth=1.5, label='Chl-a', alpha=0.9)
    ax_outer3.set_ylabel('Chl-a (mg m⁻³)', color=COLOR_CHL, fontsize=11, fontweight='bold')
    ax_outer3.tick_params(axis='y', labelcolor=COLOR_CHL)

    ax_outer.set_title('Panel B: Outer Arc - Environmental Time Series (14-day smooth)',
                       fontsize=13, fontweight='bold', loc='left')
    ax_outer.grid(True, alpha=0.3)

    # ========== Panel C: Pearson Synchrony ==========
    ax_pearson = fig.add_subplot(gs[2, 0])

    if len(pearson_df) > 0:
        pair_labels = ['SST-DLI', 'SST-CHL', 'DLI-CHL']
        pair_cols = ['Corr_SST_DLI', 'Corr_SST_CHL', 'Corr_DLI_CHL']

        plot_data = []
        plot_colors = []
        positions = []
        pos = 0

        for i, (label, col) in enumerate(zip(pair_labels, pair_cols)):
            if col not in pearson_df.columns:
                continue

            inner_vals = pearson_df[pearson_df['Arc'] == 'inner'][col].dropna()
            outer_vals = pearson_df[pearson_df['Arc'] == 'outer'][col].dropna()

            if len(inner_vals) > 0:
                plot_data.append(inner_vals.values)
                plot_colors.append(COLOR_INNER)
                positions.append(pos)
                pos += 1

            if len(outer_vals) > 0:
                plot_data.append(outer_vals.values)
                plot_colors.append(COLOR_OUTER)
                positions.append(pos)
                pos += 1

            pos += 0.3

        if len(plot_data) > 0:
            bp = ax_pearson.boxplot(plot_data, positions=positions, widths=0.6, patch_artist=True)

            for patch, color in zip(bp['boxes'], plot_colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)

            for i, (data, pos_val) in enumerate(zip(plot_data, positions)):
                x = np.random.normal(pos_val, 0.03, len(data))
                ax_pearson.scatter(x, data, alpha=0.5, s=25, c=plot_colors[i],
                                  edgecolors='white', linewidth=0.5, zorder=3)

            ax_pearson.set_xticks([positions[i] for i in range(0, len(positions), 2)])
            ax_pearson.set_xticklabels(pair_labels, fontsize=10)

            # Add arc labels
            yrange = ax_pearson.get_ylim()[1] - ax_pearson.get_ylim()[0]
            ax_pearson.text(-0.5, ax_pearson.get_ylim()[0] - 0.15 * yrange, 'Inner',
                           color=COLOR_INNER, fontsize=9, fontweight='bold', ha='center')
            ax_pearson.text(0.5, ax_pearson.get_ylim()[0] - 0.15 * yrange, 'Outer',
                           color=COLOR_OUTER, fontsize=9, fontweight='bold', ha='center')

        ax_pearson.set_ylabel('Pearson Correlation', fontsize=11)
        ax_pearson.set_ylim(-1, 1)
        ax_pearson.axhline(0, color='gray', linestyle='--', alpha=0.5)
        ax_pearson.set_title('Panel C: Pearson Synchrony (by Site)', fontsize=12, fontweight='bold')

    # ========== Panel D: Cross-Correlation (LAG) ==========
    ax_ccf = fig.add_subplot(gs[2, 1])

    # Plot CCF for Inner Arc (SST-DLI as example)
    for arc_name, arc_key, color in [('Inner', 'inner', COLOR_INNER), ('Outer', 'outer', COLOR_OUTER)]:
        if ccf_results[arc_key].get('SST_DLI') is not None:
            ccf = ccf_results[arc_key]['SST_DLI']
            ax_ccf.plot(ccf['lags'], ccf['correlations'], color=color, linewidth=2,
                       label=f'{arc_name} Arc (opt lag: {ccf["optimal_lag"]}d)')

    ax_ccf.axvline(0, color='gray', linestyle='--', alpha=0.5)
    ax_ccf.set_xlabel('Lag (days)', fontsize=10)
    ax_ccf.set_ylabel('Correlation', fontsize=10)
    ax_ccf.set_title('Panel D: Cross-Correlation SST ↔ DLI', fontsize=12, fontweight='bold')
    ax_ccf.legend(loc='upper right', fontsize=9)
    ax_ccf.set_ylim(-0.5, 0.8)
    ax_ccf.grid(True, alpha=0.3)

    # ========== Panel E: Spectral Coherence ==========
    ax_spectral = fig.add_subplot(gs[3, 0])

    if len(pearson_df) > 0:
        coh_labels = ['SST-DLI\n(7d)', 'SST-CHL\n(7d)', 'DLI-CHL\n(7d)',
                      'SST-DLI\n(30d)', 'SST-CHL\n(30d)', 'DLI-CHL\n(30d)']

        plot_data = []
        plot_colors = []
        positions = []
        pos = 0

        for arc in ['inner', 'outer']:
            for var_pair in ['SST_DLI', 'SST_CHL', 'DLI_CHL']:
                for period in ['7d', '30d']:
                    coh_key = f'coh_{period.replace("d", "")}'

                    if coh_results[arc].get(var_pair) is not None:
                        coh_val = coh_results[arc][var_pair].get(coh_key, np.nan)
                        if not np.isnan(coh_val):
                            plot_data.append([coh_val])
                            plot_colors.append(COLOR_INNER if arc == 'inner' else COLOR_OUTER)
                            positions.append(pos)
                            pos += 1

        if len(plot_data) > 0:
            bp = ax_spectral.boxplot(plot_data, positions=positions, widths=0.5, patch_artist=True)

            for patch, color in zip(bp['boxes'], plot_colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)

            ax_spectral.set_xticks(positions)
            ax_spectral.set_xticklabels(coh_labels[:len(positions)], fontsize=8)

        ax_spectral.set_ylabel('Coherence', fontsize=10)
        ax_spectral.set_ylim(0, 1)
        ax_spectral.set_title('Panel E: Spectral Coherence', fontsize=12, fontweight='bold')

    # ========== Panel F: Pulse Synchrony (STL Residuals) ==========
    ax_stl = fig.add_subplot(gs[3, 1])

    if len(pearson_df) > 0 and stl_results['inner'] is not None and stl_results['outer'] is not None:
        stl_labels = ['SST-DLI', 'SST-CHL', 'DLI-CHL']
        stl_keys = ['Pulse_Corr_SST_DLI', 'Pulse_Corr_SST_CHL', 'Pulse_Corr_DLI_CHL']

        plot_data = []
        plot_colors = []
        positions = []
        pos = 0

        for arc in ['inner', 'outer']:
            for key in stl_keys:
                val = stl_results[arc].get(key, np.nan)
                if not np.isnan(val):
                    plot_data.append([val])
                    plot_colors.append(COLOR_INNER if arc == 'inner' else COLOR_OUTER)
                    positions.append(pos)
                    pos += 1

        if len(plot_data) > 0:
            bp = ax_stl.boxplot(plot_data, positions=positions, widths=0.5, patch_artist=True)

            for patch, color in zip(bp['boxes'], plot_colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)

            ax_stl.set_xticks(positions)
            ax_stl.set_xticklabels(stl_labels[:len(positions)], fontsize=9)

        ax_stl.set_ylabel('Residual Correlation', fontsize=10)
        ax_stl.set_ylim(-1, 1)
        ax_stl.axhline(0, color='gray', linestyle='--', alpha=0.5)
        ax_stl.set_title('Panel F: Pulse Synchrony (STL Residuals)', fontsize=12, fontweight='bold')

    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"    OK Saved: {os.path.basename(output_path)}")

def create_raw_timeseries_figure(sst_arc, dli_arc, chl_arc, output_path):
    """
    Create figure with raw daily data (points) for high-resolution temporal visualization.
    Following the style of TIME_SERIES_LAG_temporal_environ.py
    """
    print("\n  Creating Raw Time Series Figure...")

    fig, axes = plt.subplots(2, 1, figsize=(18, 10), sharex=True)
    fig.suptitle('Raw Daily Environmental Data (2002-2008)', fontsize=14, fontweight='bold')

    for i, (arc, title, color_arc) in enumerate([('inner', 'Inner Arc', COLOR_INNER),
                                                  ('outer', 'Outer Arc', COLOR_OUTER)]):
        ax = axes[i]

        sst = sst_arc[arc]
        dli = dli_arc[arc]
        chl = chl_arc[arc]

        if sst is not None and len(sst) > 0:
            ax.scatter(sst.index, sst.values, s=2, alpha=0.6, c=COLOR_SST, label='SST')

        ax.set_ylabel('SST (°C)', color=COLOR_SST, fontsize=11, fontweight='bold')
        ax.tick_params(axis='y', labelcolor=COLOR_SST)
        ax.set_title(f'{title}', fontsize=12, fontweight='bold', loc='left')
        ax.grid(True, alpha=0.3)

        ax2 = ax.twinx()
        if dli is not None and len(dli) > 0:
            ax2.scatter(dli.index, dli.values, s=2, alpha=0.6, c=COLOR_DLI, label='DLI')
        ax2.set_ylabel('DLI (mol m⁻² d⁻¹)', color=COLOR_DLI, fontsize=11, fontweight='bold')
        ax2.tick_params(axis='y', labelcolor=COLOR_DLI)

        ax3 = ax.twinx()
        ax3.spines['right'].set_position(('outward', 60))
        if chl is not None and len(chl) > 0:
            ax3.scatter(chl.index, chl.values, s=2, alpha=0.6, c=COLOR_CHL, label='Chl-a')
        ax3.set_ylabel('Chl-a (mg m⁻³)', color=COLOR_CHL, fontsize=11, fontweight='bold')
        ax3.tick_params(axis='y', labelcolor=COLOR_CHL)

    import matplotlib.dates as mdates
    axes[1].xaxis.set_major_locator(mdates.YearLocator())
    axes[1].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"    OK Saved: {os.path.basename(output_path)}")

def create_ccf_detail_figure(ccf_results, output_path):
    """
    Create detailed cross-correlation figure showing all variable pairs.
    """
    print("\n  Creating Cross-Correlation Detail Figure...")

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle('Cross-Correlation Analysis by Variable Pair', fontsize=14, fontweight='bold')

    var_pairs = [('SST', 'DLI'), ('SST', 'CHL'), ('DLI', 'CHL')]
    pair_titles = ['SST ↔ DLI', 'SST ↔ Chl-a', 'DLI ↔ Chl-a']

    for col, (var1, var2, title) in enumerate(zip(['SST', 'SST', 'DLI'],
                                                    ['DLI', 'CHL', 'CHL'],
                                                    pair_titles)):
        for row, (arc, arc_name, color) in enumerate([('inner', 'Inner Arc', COLOR_INNER),
                                                       ('outer', 'Outer Arc', COLOR_OUTER)]):
            ax = axes[row, col]

            key = f'{var1}_{var2}'
            if ccf_results[arc].get(key) is not None:
                ccf = ccf_results[arc][key]

                ax.plot(ccf['lags'], ccf['correlations'], color=color, linewidth=2)

                # Mark optimal lag
                ax.axvline(ccf['optimal_lag'], color='red', linestyle='--', alpha=0.7,
                          label=f'Optimal lag: {ccf["optimal_lag"]}d')

                # Add correlation value
                ax.text(0.05, 0.95, f'r = {ccf["optimal_corr"]:.3f}',
                       transform=ax.transAxes, ha='left', va='top',
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))

            ax.axvline(0, color='gray', linestyle='-', alpha=0.3)
            ax.set_xlabel('Lag (days)', fontsize=9)
            ax.set_ylabel('Correlation', fontsize=9)
            ax.set_ylim(-0.6, 0.8)
            ax.grid(True, alpha=0.3)
            ax.legend(loc='upper right', fontsize=8)

            if row == 0:
                ax.set_title(title, fontsize=11, fontweight='bold')

            if col == 0:
                ax.text(-0.15, 0.5, arc_name, transform=ax.transAxes,
                       fontsize=10, fontweight='bold', rotation=90, va='center', color=color)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"    OK Saved: {os.path.basename(output_path)}")

# ============================================
# Main Execution
# ============================================

if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("STEP 1: Loading Satellite Data")
    print("=" * 60)

    # List and filter files
    print("\n  Finding files...")
    sst_files = filter_files_by_period(glob.glob(os.path.join(SST_DIR, SST_PATTERN)), START_YEAR, END_YEAR)
    kd490_files = filter_files_by_period(glob.glob(os.path.join(MODIS_DIR, KD490_PATTERN)), START_YEAR, END_YEAR)
    par_files = filter_files_by_period(glob.glob(os.path.join(MODIS_DIR, PAR_PATTERN)), START_YEAR, END_YEAR)
    chl_files = filter_files_by_period(glob.glob(os.path.join(MODIS_DIR, CHL_PATTERN)), START_YEAR, END_YEAR)

    print(f"    SST files: {len(sst_files)}")
    print(f"    Kd490 files: {len(kd490_files)}")
    print(f"    PAR files: {len(par_files)}")
    print(f"    Chl-a files: {len(chl_files)}")

    # Load time series for each variable
    sst_ts = load_timeseries_for_sites(sst_files, ['analysed_sst', 'sea_surface_temperature', 'sst'], sites_df, 'SST')
    chl_ts = load_timeseries_for_sites(chl_files, ['chlor_a'], sites_df, 'Chl-a')
    dli_ts = calculate_dli_for_sites(kd490_files, par_files, sites_df)

    print("\n" + "=" * 60)
    print("STEP 2: Aggregating by Arc")
    print("=" * 60)

    # Aggregate by arc
    sst_arc = aggregate_timeseries_by_arc(sst_ts, sites_df)
    dli_arc = aggregate_timeseries_by_arc(dli_ts, sites_df)
    chl_arc = aggregate_timeseries_by_arc(chl_ts, sites_df)

    print(f"  Inner Arc SST: {len(sst_arc['inner']) if sst_arc['inner'] is not None else 0} days")
    print(f"  Outer Arc SST: {len(sst_arc['outer']) if sst_arc['outer'] is not None else 0} days")

    print("\n" + "=" * 60)
    print("STEP 3: Calculating Pearson Synchrony")
    print("=" * 60)

    pearson_df = calculate_pearson_synchrony(sst_ts, dli_ts, chl_ts, sites_df)
    print(f"  OK Pearson correlations calculated for {len(pearson_df)} sites")

    print("\n" + "=" * 60)
    print("STEP 4: Calculating Cross-Correlation with Lag")
    print("=" * 60)

    ccf_results = calculate_cross_correlation_by_arc(sst_arc, dli_arc, chl_arc, max_lag=60)
    print("  OK Cross-correlation analysis complete")

    print("\n" + "=" * 60)
    print("STEP 5: Calculating Pulse Synchrony (STL Decomposition)")
    print("=" * 60)

    stl_results = calculate_pulse_synchrony_by_arc(sst_arc, dli_arc, chl_arc)
    print("  OK STL decomposition complete")

    print("\n" + "=" * 60)
    print("STEP 6: Calculating Spectral Coherence")
    print("=" * 60)

    coh_results = calculate_spectral_coherence_by_arc(sst_arc, dli_arc, chl_arc)
    print("  OK Spectral coherence analysis complete")

    print("\n" + "=" * 60)
    print("STEP 7: Generating Figures")
    print("=" * 60)

    # Main synchrony figure
    create_main_synchrony_figure(
        sst_arc, dli_arc, chl_arc,
        pearson_df, ccf_results, stl_results, coh_results,
        os.path.join(OUTPUT_DIR, 'Fig_Main_Temporal_Synchrony.png')
    )

    # Raw time series figure
    create_raw_timeseries_figure(
        sst_arc, dli_arc, chl_arc,
        os.path.join(OUTPUT_DIR, 'Fig_S2_Raw_Time_Series.png')
    )

    # Cross-correlation detail figure
    create_ccf_detail_figure(
        ccf_results,
        os.path.join(OUTPUT_DIR, 'Fig_S3_Cross_Correlation_Detail.png')
    )

    print("\n" + "=" * 60)
    print("STEP 8: Saving Results")
    print("=" * 60)

    # Save Pearson correlations
    pearson_df.to_csv(os.path.join(OUTPUT_DIR, 'Pearson_Synchrony_by_Site.csv'), index=False)
    print("  OK Saved: Pearson_Synchrony_by_Site.csv")

    # Save cross-correlation results
    ccf_summary = []
    for arc in ['inner', 'outer']:
        for pair, result in ccf_results[arc].items():
            if result is not None:
                ccf_summary.append({
                    'Arc': arc,
                    'Variable_Pair': pair,
                    'Optimal_Lag_Days': result['optimal_lag'],
                    'Optimal_Correlation': result['optimal_corr']
                })
    pd.DataFrame(ccf_summary).to_csv(os.path.join(OUTPUT_DIR, 'Cross_Correlation_Results.csv'), index=False)
    print("  OK Saved: Cross_Correlation_Results.csv")

    # Save STL pulse synchrony
    stl_summary = []
    for arc in ['inner', 'outer']:
        if stl_results[arc] is not None:
            stl_summary.append({
                'Arc': arc,
                'Pulse_Corr_SST_DLI': stl_results[arc].get('Pulse_Corr_SST_DLI', np.nan),
                'Pulse_Corr_SST_CHL': stl_results[arc].get('Pulse_Corr_SST_CHL', np.nan),
                'Pulse_Corr_DLI_CHL': stl_results[arc].get('Pulse_Corr_DLI_CHL', np.nan)
            })
    pd.DataFrame(stl_summary).to_csv(os.path.join(OUTPUT_DIR, 'STL_Pulse_Synchrony.csv'), index=False)
    print("  OK Saved: STL_Pulse_Synchrony.csv")

    print("\n" + "=" * 60)
    print("SYNCHRONY ANALYSIS COMPLETE!")
    print(f"All outputs saved to: {OUTPUT_DIR}")
    print("=" * 60)
