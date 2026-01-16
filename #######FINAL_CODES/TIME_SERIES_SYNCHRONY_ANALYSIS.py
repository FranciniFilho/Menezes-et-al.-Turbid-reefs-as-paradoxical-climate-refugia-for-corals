### TIME_SERIES_SYNCHRONY_ANALYSIS.py ###
# Temporal Synchrony Analysis: Inner Arc vs Outer Arc Environmental Coupling
#
# This script evaluates whether Inner Arc reefs exhibit higher temporal
# synchrony (coupling) between SST, DLI, and Chl-a than Outer Arc reefs.
#
# Features:
# 1. Publication-quality time series (smooth + raw data)
# 2. Pearson correlation analysis for synchrony
# 3. Spectral coherence analysis (Lomb-Scargle)
# 4. Statistical comparison (Mann-Whitney U)
#
# Period: 2002-2008
# Reference style: TIME_SERIES_VARIABILITY_ANALYSIS_GEMINI.py

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
from scipy import stats
from scipy.signal import lombscargle
from datetime import datetime
from scipy.stats import mannwhitneyu
import warnings
warnings.filterwarnings('ignore')

# ============================================
# Configuration (Following GEMINI.py style)
# ============================================

# --- Time Period ---
START_YEAR = 2002
END_YEAR = 2008

# --- Data Directories ---
SST_DIR = r'H:\remote sensing\CRW_SST_FULL'
MODIS_DIR = r'H:\remote sensing\MODIS_DATA_FULL'

# --- Sites CSV ---
SITES_CSV_PATH = r'C:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######22.04.23\DATA\sites_list_full_clean.csv'

# --- Output Directory ---
OUTPUT_DIR = r'C:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_RESULTS\Temporal_Synchrony_Analysis'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Region of Interest ---
LAT_MIN, LAT_MAX = -18.5, -16.5
LON_MIN, LON_MAX = -39.5, -38.0

# --- File Patterns ---
SST_PATTERN = 'coraltemp_v3.1_*.nc'
KD490_PATTERN = 'AQUA_MODIS.*.L3m.DAY.KD.Kd_490.4km.nc'
PAR_PATTERN = 'AQUA_MODIS.*.L3m.DAY.PAR.par.4km.nc'
CHL_PATTERN = 'AQUA_MODIS.*.L3m.DAY.CHL.chlor_a.4km.nc'

# --- DLI Calculation Constants ---
KDPAR_A, KDPAR_B, KDPAR_C = 0.0665, 0.874, 0.00121

# --- Color Palette (From GEMINI.py) ---
COLOR_INNER = '#CD853F'  # Light brown (Peru)
COLOR_OUTER = '#4169E1'  # Royal Blue
COLORS_ARC = {'inner': COLOR_INNER, 'outer': COLOR_OUTER}

# --- Variable Colors ---
COLOR_SST = '#D55E00'    # Professional red-orange
COLOR_DLI = '#0072B2'    # Professional blue
COLOR_CHL = '#009E73'    # Professional green

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
# Logging Configuration
# ============================================
logging.basicConfig(
    filename=os.path.join(OUTPUT_DIR, 'temporal_synchrony_analysis_errors.log'),
    filemode='w',
    level=logging.ERROR,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# ============================================
# Print Header
# ============================================
print("=" * 60)
print("TEMPORAL SYNCHRONY ANALYSIS: Inner vs Outer Arc")
print(f"Period: {START_YEAR} - {END_YEAR}")
print("=" * 60)

# ============================================
# Phase 1: Load Sites Data
# ============================================

def load_sites_data():
    """Load and standardize sites data."""
    print("\n" + "=" * 60)
    print("PHASE 1: Loading Sites Data")
    print("=" * 60)

    try:
        sites_df = pd.read_csv(SITES_CSV_PATH, sep=';')
        sites_df.columns = [col.strip() for col in sites_df.columns]

        # Verify required columns
        required_cols = {'Site_name', 'Latitude', 'Longitude', 'Depth_m', 'Arc'}
        if not required_cols.issubset(sites_df.columns):
            missing = required_cols - set(sites_df.columns)
            raise ValueError(f"Missing columns: {missing}. Available: {list(sites_df.columns)}")

        # Clean and standardize
        sites_df['Arc'] = sites_df['Arc'].str.strip().str.lower()
        sites_df = sites_df[sites_df['Arc'].isin(['inner', 'outer'])].copy()

        print(f"\nLoaded {len(sites_df)} sites")
        print(f"  Inner Arc: {(sites_df['Arc'] == 'inner').sum()} sites")
        print(f"  Outer Arc: {(sites_df['Arc'] == 'outer').sum()} sites")

        return sites_df

    except FileNotFoundError:
        print(f"ERROR: Sites file not found: {SITES_CSV_PATH}")
        exit(1)
    except Exception as e:
        print(f"ERROR loading sites: {e}")
        exit(1)

sites_df = load_sites_data()

# ============================================
# Phase 2: Utility Functions
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

# ============================================
# Phase 3: Load Satellite Data
# ============================================

def load_timeseries_for_sites(files, var_names, sites_df, var_label):
    """Load time series data for all sites."""
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
            with xr.open_dataset(f) as ds:
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
                    except Exception:
                        continue

                success_count += 1

        except Exception as e:
            error_count += 1
            if error_count <= 5:
                logging.error(f"Error reading {f}: {e}")
                print(f"    WARNING: Error reading {os.path.basename(f)}: {str(e)[:50]}")
            continue

    print(f"    Files processed successfully: {success_count}/{total_files}")
    if error_count > 0:
        print(f"    Files with errors: {error_count}")

    # Convert to pandas Series
    result = {}
    for site_name, data in site_data.items():
        if len(data['times']) > 0:
            result[site_name] = pd.Series(data['values'], index=pd.DatetimeIndex(data['times'])).sort_index()

    print(f"    Loaded data for {len(result)} sites")
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

    for i, f_kd in enumerate(kd490_files):
        if (i + 1) % 500 == 0:
            print(f"    Processing file {i+1}/{total_files}...")

        date = parse_date_from_filename(f_kd)
        if date is None or date not in par_by_date:
            continue

        f_par = par_by_date[date]

        try:
            with xr.open_dataset(f_kd) as ds_kd, xr.open_dataset(f_par) as ds_par:
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
                print(f"    WARNING: Error with DLI calculation: {str(e)[:50]}")
            continue

    print(f"    Files processed successfully: {success_count}/{total_files}")
    if error_count > 0:
        print(f"    Files with errors: {error_count}")

    # Convert to pandas Series
    result = {}
    for site_name, data in site_data.items():
        if len(data['times']) > 0:
            result[site_name] = pd.Series(data['values'], index=pd.DatetimeIndex(data['times'])).sort_index()

    print(f"    Calculated DLI for {len(result)} sites")
    return result

def aggregate_timeseries_by_arc(timeseries_dict, sites_df):
    """Aggregate time series by arc (mean of all sites in each arc)."""
    print("\n  Aggregating time series by arc...")

    arc_series = {'inner': [], 'outer': []}

    for site_name, series in timeseries_dict.items():
        if site_name not in sites_df['Site_name'].values:
            continue

        arc = sites_df[sites_df['Site_name'] == site_name]['Arc'].values[0]
        if arc in ['inner', 'outer']:
            arc_series[arc].append(series)

    # Concatenate and take mean for each arc
    result = {}
    for arc in ['inner', 'outer']:
        if len(arc_series[arc]) > 0:
            # Align all series to common dates
            aligned = pd.concat(arc_series[arc], axis=1).mean(axis=1)
            result[arc] = aligned
            print(f"    {arc.capitalize()} Arc: {len(aligned)} days from {len(arc_series[arc])} sites")

    return result

# ============================================
# Load all data
# ============================================
print("\n" + "=" * 60)
print("PHASE 2: Loading Satellite Data")
print("=" * 60)

print("\n  Finding files...")
sst_files = filter_files_by_period(glob.glob(os.path.join(SST_DIR, SST_PATTERN)), START_YEAR, END_YEAR)
kd490_files = filter_files_by_period(glob.glob(os.path.join(MODIS_DIR, KD490_PATTERN)), START_YEAR, END_YEAR)
par_files = filter_files_by_period(glob.glob(os.path.join(MODIS_DIR, PAR_PATTERN)), START_YEAR, END_YEAR)
chl_files = filter_files_by_period(glob.glob(os.path.join(MODIS_DIR, CHL_PATTERN)), START_YEAR, END_YEAR)

print(f"    SST files: {len(sst_files)}")
print(f"    Kd490 files: {len(kd490_files)}")
print(f"    PAR files: {len(par_files)}")
print(f"    Chl-a files: {len(chl_files)}")

# Load time series per site
sst_ts = load_timeseries_for_sites(sst_files, ['analysed_sst', 'sea_surface_temperature', 'sst'], sites_df, 'SST')
chl_ts = load_timeseries_for_sites(chl_files, ['chlor_a'], sites_df, 'Chl-a')
dli_ts = calculate_dli_for_sites(kd490_files, par_files, sites_df)

# Aggregate by arc
sst_arc = aggregate_timeseries_by_arc(sst_ts, sites_df)
dli_arc = aggregate_timeseries_by_arc(dli_ts, sites_df)
chl_arc = aggregate_timeseries_by_arc(chl_ts, sites_df)

# ============================================
# Phase 4: Synchrony Analysis Functions
# ============================================

def calculate_pearson_synchrony(sst_dict, dli_dict, chl_dict, sites_df):
    """
    Calculate Pearson correlation synchrony for each site.

    Returns DataFrame with synchrony metrics per site.
    """
    print("\n  Calculating Pearson correlation synchrony...")

    results = []

    for _, row in sites_df.iterrows():
        site_name = row['Site_name']
        arc = row['Arc']

        # Get series for this site
        sst = sst_dict.get(site_name)
        dli = dli_dict.get(site_name)
        chl = chl_dict.get(site_name)

        if sst is None or dli is None or chl is None:
            continue

        # Align to common dates
        df = pd.DataFrame({'SST': sst, 'DLI': dli, 'CHL': chl}).dropna()

        if len(df) < 30:  # Minimum data requirement
            continue

        # Standardize (z-score)
        df_std = (df - df.mean()) / df.std()

        # Calculate pairwise correlations
        corr_sst_dli = df_std['SST'].corr(df_std['DLI'])
        corr_sst_chl = df_std['SST'].corr(df_std['CHL'])
        corr_dli_chl = df_std['DLI'].corr(df_std['CHL'])

        # Mean synchrony (Fisher z-transform for averaging)
        # Convert to z-scale, average, convert back
        z_sst_dli = np.arctanh(corr_sst_dli) if abs(corr_sst_dli) < 1 else 0
        z_sst_chl = np.arctanh(corr_sst_chl) if abs(corr_sst_chl) < 1 else 0
        z_dli_chl = np.arctanh(corr_dli_chl) if abs(corr_dli_chl) < 1 else 0

        mean_z = np.nanmean([z_sst_dli, z_sst_chl, z_dli_chl])
        mean_synchrony = np.tanh(mean_z)

        results.append({
            'Site_name': site_name,
            'Arc': arc,
            'Synchrony_Mean': mean_synchrony,
            'Corr_SST_DLI': corr_sst_dli,
            'Corr_SST_CHL': corr_sst_chl,
            'Corr_DLI_CHL': corr_dli_chl,
            'N': len(df)
        })

    df_result = pd.DataFrame(results)
    print(f"    Calculated for {len(df_result)} sites")
    print(f"      Inner: {(df_result['Arc'] == 'inner').sum()} sites")
    print(f"      Outer: {(df_result['Arc'] == 'outer').sum()} sites")

    return df_result

def compute_spectral_coherence(series1, series2, period_of_interest=30):
    """
    Compute spectral coherence at a specific period using Lomb-Scargle.

    Returns coherence magnitude at the target period.
    """
    # Align and remove NaN
    df = pd.DataFrame({'s1': series1, 's2': series2}).dropna()

    if len(df) < 30:
        return np.nan

    # Normalize
    s1_norm = (df['s1'] - df['s1'].mean()) / (df['s1'].std() + 1e-9)
    s2_norm = (df['s2'] - df['s2'].mean()) / (df['s2'].std() + 1e-9)

    # Prepare times
    times = (df.index - df.index[0]).days.values.astype(float)

    # Define frequency range
    max_period = len(times) // 2
    min_period = 3
    periods = np.linspace(min_period, max_period, 500)
    frequencies = 2 * np.pi / periods

    # Compute Lomb-Scargle for both series
    try:
        power1 = lombscargle(times, s1_norm.values, frequencies, normalize=True)
        power2 = lombscargle(times, s2_norm.values, frequencies, normalize=True)

        # Cross-power (product of normalized powers at same frequency)
        cross_power = np.sqrt(power1 * power2)

        # Extract coherence at target period
        mask = np.abs(periods - period_of_interest) / period_of_interest <= 0.3
        if mask.sum() == 0:
            return np.nan

        return cross_power[mask].max()

    except Exception as e:
        logging.error(f"Spectral coherence error: {e}")
        return np.nan

def calculate_spectral_synchrony(sst_dict, dli_dict, chl_dict, sites_df):
    """
    Calculate spectral coherence synchrony for each site.

    Evaluates synchrony at key periods: 7 days, 30 days.
    """
    print("\n  Calculating spectral coherence synchrony...")

    results = []

    for _, row in sites_df.iterrows():
        site_name = row['Site_name']
        arc = row['Arc']

        sst = sst_dict.get(site_name)
        dli = dli_dict.get(site_name)
        chl = chl_dict.get(site_name)

        if sst is None or dli is None or chl is None:
            continue

        # Calculate coherence at 7-day and 30-day periods
        coh_7d_sst_dli = compute_spectral_coherence(sst, dli, 7)
        coh_7d_sst_chl = compute_spectral_coherence(sst, chl, 7)
        coh_7d_dli_chl = compute_spectral_coherence(dli, chl, 7)

        coh_30d_sst_dli = compute_spectral_coherence(sst, dli, 30)
        coh_30d_sst_chl = compute_spectral_coherence(sst, chl, 30)
        coh_30d_dli_chl = compute_spectral_coherence(dli, chl, 30)

        # Mean coherence at each period
        coh_7d_list = [coh_7d_sst_dli, coh_7d_sst_chl, coh_7d_dli_chl]
        coh_30d_list = [coh_30d_sst_dli, coh_30d_sst_chl, coh_30d_dli_chl]

        mean_coh_7d = np.nanmean(coh_7d_list) if not np.isnan(coh_7d_list).all() else np.nan
        mean_coh_30d = np.nanmean(coh_30d_list) if not np.isnan(coh_30d_list).all() else np.nan

        results.append({
            'Site_name': site_name,
            'Arc': arc,
            'Coherence_7d': mean_coh_7d,
            'Coherence_30d': mean_coh_30d
        })

    df_result = pd.DataFrame(results)
    print(f"    Calculated for {len(df_result)} sites")
    print(f"      Inner: {(df_result['Arc'] == 'inner').sum()} sites")
    print(f"      Outer: {(df_result['Arc'] == 'outer').sum()} sites")

    return df_result

# ============================================
# Phase 5: Statistical Testing
# ============================================

def compare_synchrony_arcs(synchrony_df, metric_col):
    """
    Perform Mann-Whitney U test to compare synchrony between Inner and Outer arcs.

    Returns: p_value, effect_size (r), formatted string
    """
    inner_vals = synchrony_df[synchrony_df['Arc'] == 'inner'][metric_col].dropna()
    outer_vals = synchrony_df[synchrony_df['Arc'] == 'outer'][metric_col].dropna()

    if len(inner_vals) < 3 or len(outer_vals) < 3:
        return np.nan, np.nan, "N/A"

    try:
        stat, p_value = mannwhitneyu(inner_vals, outer_vals, alternative='greater')

        # Calculate effect size (r)
        n1, n2 = len(inner_vals), len(outer_vals)
        # Using U statistic to calculate z
        mean_U = n1 * n2 / 2
        std_U = np.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
        z_score = (stat - mean_U) / std_U
        r = z_score / np.sqrt(n1 + n2)

        # Format p-value
        if p_value < 0.001:
            p_str = "p < 0.001"
        elif p_value < 0.01:
            p_str = f"p = {p_value:.3f}"
        elif p_value < 0.05:
            p_str = f"p = {p_value:.3f}"
        else:
            p_str = f"p = {p_value:.3f}"

        return p_value, r, p_str

    except Exception as e:
        logging.error(f"Mann-Whitney U test error: {e}")
        return np.nan, np.nan, "N/A"

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
# Phase 6: Calculate Synchrony Metrics
# ============================================
print("\n" + "=" * 60)
print("PHASE 3: Calculating Synchrony Metrics")
print("=" * 60)

pearson_sync = calculate_pearson_synchrony(sst_ts, dli_ts, chl_ts, sites_df)
spectral_sync = calculate_spectral_synchrony(sst_ts, dli_ts, chl_ts, sites_df)

# ============================================
# Phase 7: Statistical Comparison
# ============================================
print("\n" + "=" * 60)
print("PHASE 4: Statistical Comparison")
print("=" * 60)

print("\n  Pearson Synchrony (Inner > Outer):")
p_val_pearson, r_eff_pearson, p_str_pearson = compare_synchrony_arcs(pearson_sync, 'Synchrony_Mean')
print(f"    {p_str_pearson}, effect size r = {r_eff_pearson:.3f}")

print("\n  Spectral Coherence 7d (Inner > Outer):")
p_val_7d, r_eff_7d, p_str_7d = compare_synchrony_arcs(spectral_sync, 'Coherence_7d')
print(f"    {p_str_7d}, effect size r = {r_eff_7d:.3f}")

print("\n  Spectral Coherence 30d (Inner > Outer):")
p_val_30d, r_eff_30d, p_str_30d = compare_synchrony_arcs(spectral_sync, 'Coherence_30d')
print(f"    {p_str_30d}, effect size r = {r_eff_30d:.3f}")

# ============================================
# Phase 8: Visualization Functions
# ============================================

def smooth_timeseries(series, window_days=14):
    """
    Apply smoothing to time series using rolling mean.

    Window should be odd for symmetric smoothing
    """
    if len(series) < window_days:
        return series

    min_periods = max(1, int(window_days * 0.5))
    smoothed = series.rolling(window=window_days, min_periods=min_periods, center=True).mean()
    return smoothed

def create_publication_figure(
    sst_arc, dli_arc, chl_arc,
    pearson_sync, spectral_sync,
    output_path
):
    """
    Create publication-quality 2x2 figure.

    Layout:
    ├── ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─ ─
    │ Panel A: Inner Arc Time Series
    │ Panel B: Outer Arc Time Series
    ├───────────────────────────────────────
    │ Panel C: Pearson Synchrony Comparison
    │ Panel D: Spectral Coherence Comparison
    └───────────────────────────────────────
    """
    print("\n  Creating publication figure...")

    # Increased figure width to accommodate 3 Y-axes
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.35, top=0.92, bottom=0.08, left=0.08, right=0.92)

    # ===============================
    # Panel A: Inner Arc Time Series (3 Y-axes)
    # ===============================
    ax_inner = fig.add_subplot(gs[0, 0])

    # Get data
    sst_i = sst_arc['inner']
    dli_i = dli_arc['inner']
    chl_i = chl_arc['inner']

    # Smooth (14-day window)
    sst_i_smooth = smooth_timeseries(sst_i, window_days=14)
    dli_i_smooth = smooth_timeseries(dli_i, window_days=14)
    chl_i_smooth = smooth_timeseries(chl_i, window_days=14)

    # EIXO 1 (Esquerdo): SST
    p1, = ax_inner.plot(sst_i_smooth.index, sst_i_smooth.values, '-', lw=1.5, c=COLOR_SST, label='SST', zorder=2)
    ax_inner.set_ylabel('SST (°C)', color=COLOR_SST, fontsize=12, fontweight='bold')
    ax_inner.tick_params(axis='y', labelcolor=COLOR_SST, labelsize=10)
    ax_inner.grid(axis='y', linestyle='--', alpha=0.3, color=COLOR_SST)

    # EIXO 2 (Direito 1): DLI
    ax_inner_dli = ax_inner.twinx()
    p2, = ax_inner_dli.plot(dli_i_smooth.index, dli_i_smooth.values, '-', lw=1.5, c=COLOR_DLI, label='DLI', zorder=2)
    ax_inner_dli.set_ylabel('DLI (mol m⁻² d⁻¹)', color=COLOR_DLI, fontsize=12, fontweight='bold')
    ax_inner_dli.tick_params(axis='y', labelcolor=COLOR_DLI, labelsize=10)

    # EIXO 3 (Direito 2): Chl-a (deslocado 70 pontos para fora)
    ax_inner_chl = ax_inner.twinx()
    ax_inner_chl.spines['right'].set_position(('outward', 70))
    p3, = ax_inner_chl.plot(chl_i_smooth.index, chl_i_smooth.values, '-', lw=1.5, c=COLOR_CHL, label='Chl-a', zorder=2)
    ax_inner_chl.set_ylabel('Chl-a (mg m⁻³)', color=COLOR_CHL, fontsize=12, fontweight='bold')
    ax_inner_chl.tick_params(axis='y', labelcolor=COLOR_CHL, labelsize=10)

    # Title and x-axis
    ax_inner.set_title('Inner Arc', fontsize=14, fontweight='bold', loc='left', pad=10)
    ax_inner.set_xlim(pd.Timestamp(f'{START_YEAR}-01-01'), pd.Timestamp(f'{END_YEAR}-12-31'))

    # Format x-axis
    import matplotlib.dates as mdates
    ax_inner.xaxis.set_major_locator(mdates.YearLocator(1))
    ax_inner.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax_inner.tick_params(axis='x', labelsize=10)

    # LEGENDA UNIFICADA (no eixo principal, posição ajustada)
    lines = [p1, p2, p3]
    ax_inner.legend(lines, ['SST', 'DLI', 'Chl-a'], loc='upper left', fontsize=9,
                    framealpha=0.95, edgecolor='none', bbox_to_anchor=(0, 1))

    # ===============================
    # Panel B: Outer Arc Time Series (3 Y-axes)
    # ===============================
    ax_outer = fig.add_subplot(gs[0, 1])

    sst_o = sst_arc['outer']
    dli_o = dli_arc['outer']
    chl_o = chl_arc['outer']

    sst_o_smooth = smooth_timeseries(sst_o, window_days=14)
    dli_o_smooth = smooth_timeseries(dli_o, window_days=14)
    chl_o_smooth = smooth_timeseries(chl_o, window_days=14)

    # EIXO 1 (Esquerdo): SST
    ax_outer.plot(sst_o_smooth.index, sst_o_smooth.values, '-', lw=1.5, c=COLOR_SST, zorder=2)
    ax_outer.set_ylabel('SST (°C)', color=COLOR_SST, fontsize=12, fontweight='bold')
    ax_outer.tick_params(axis='y', labelcolor=COLOR_SST, labelsize=10)
    ax_outer.grid(axis='y', linestyle='--', alpha=0.3, color=COLOR_SST)

    # EIXO 2 (Direito 1): DLI
    ax_outer_dli = ax_outer.twinx()
    ax_outer_dli.plot(dli_o_smooth.index, dli_o_smooth.values, '-', lw=1.5, c=COLOR_DLI, zorder=2)
    ax_outer_dli.set_ylabel('DLI (mol m⁻² d⁻¹)', color=COLOR_DLI, fontsize=12, fontweight='bold')
    ax_outer_dli.tick_params(axis='y', labelcolor=COLOR_DLI, labelsize=10)

    # EIXO 3 (Direito 2): Chl-a (deslocado 70 pontos para fora)
    ax_outer_chl = ax_outer.twinx()
    ax_outer_chl.spines['right'].set_position(('outward', 70))
    ax_outer_chl.plot(chl_o_smooth.index, chl_o_smooth.values, '-', lw=1.5, c=COLOR_CHL, zorder=2)
    ax_outer_chl.set_ylabel('Chl-a (mg m⁻³)', color=COLOR_CHL, fontsize=12, fontweight='bold')
    ax_outer_chl.tick_params(axis='y', labelcolor=COLOR_CHL, labelsize=10)

    # Title and x-axis
    ax_outer.set_title('Outer Arc', fontsize=14, fontweight='bold', loc='left', pad=10)
    ax_outer.set_xlim(pd.Timestamp(f'{START_YEAR}-01-01'), pd.Timestamp(f'{END_YEAR}-12-31'))

    # Format x-axis
    ax_outer.xaxis.set_major_locator(mdates.YearLocator(1))
    ax_outer.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax_outer.tick_params(axis='x', labelsize=10)

    # ===============================
    # Panel C: Pearson Synchrony
    # ===============================
    ax_pearson = fig.add_subplot(gs[1, 0])

    inner_sync = pearson_sync[pearson_sync['Arc'] == 'inner']['Synchrony_Mean'].dropna()
    outer_sync = pearson_sync[pearson_sync['Arc'] == 'outer']['Synchrony_Mean'].dropna()

    bp = ax_pearson.boxplot([inner_sync.values, outer_sync.values],
                           labels=['Inner\nArc', 'Outer\nArc'],
                           patch_artist=True, widths=0.5)

    # Color boxes
    bp['boxes'][0].set_facecolor(COLOR_INNER)
    bp['boxes'][1].set_facecolor(COLOR_OUTER)
    for box in bp['boxes']:
        box.set_alpha(0.7)

    # Add individual points
    for i, (vals, color) in enumerate([(inner_sync.values, COLOR_INNER), (outer_sync.values, COLOR_OUTER)], 1):
        x = np.random.normal(i, 0.05, len(vals))
        ax_pearson.scatter(x, vals, alpha=0.6, s=40, c=color, edgecolors='white', linewidth=0.5, zorder=3)

    # Statistical test
    sig_pearson = format_significance(p_val_pearson)

    ax_pearson.text(0.95, 0.95,
                   f'{p_str_pearson}\n{sig_pearson}\nr = {r_eff_pearson:.2f}' if not np.isnan(r_eff_pearson) else f'{p_str_pearson}\n{sig_pearson}',
                   transform=ax_pearson.transAxes, ha='right', va='top', fontsize=10,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax_pearson.set_ylabel('Mean Synchrony', fontsize=12, fontweight='bold')
    ax_pearson.set_title('Pearson Correlation Synchrony', fontsize=13, fontweight='bold')

    # ===============================
    # Panel D: Spectral Coherence
    # ===============================
    ax_spec = fig.add_subplot(gs[1, 1])

    # Show coherence at 30-day period
    inner_coh = spectral_sync[spectral_sync['Arc'] == 'inner']['Coherence_30d'].dropna()
    outer_coh = spectral_sync[spectral_sync['Arc'] == 'outer']['Coherence_30d'].dropna()

    bp = ax_spec.boxplot([inner_coh.values, outer_coh.values],
                        labels=['Inner\nArc', 'Outer\nArc'],
                        patch_artist=True, widths=0.5)

    bp['boxes'][0].set_facecolor(COLOR_INNER)
    bp['boxes'][1].set_facecolor(COLOR_OUTER)
    for box in bp['boxes']:
        box.set_alpha(0.7)

    # Points
    for i, (vals, color) in enumerate([(inner_coh.values, COLOR_INNER), (outer_coh.values, COLOR_OUTER)], 1):
        x = np.random.normal(i, 0.05, len(vals))
        ax_spec.scatter(x, vals, alpha=0.6, s=40, c=color, edgecolors='white', linewidth=0.5, zorder=3)

    # Statistical test
    sig_30d = format_significance(p_val_30d)

    ax_spec.text(0.95, 0.95,
                f'{p_str_30d}\n{sig_30d}\nr = {r_eff_30d:.2f}' if not np.isnan(r_eff_30d) else f'{p_str_30d}\n{sig_30d}',
                transform=ax_spec.transAxes, ha='right', va='top', fontsize=10,
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax_spec.set_ylabel('Spectral Coherence (30-day)', fontsize=12, fontweight='bold')
    ax_spec.set_title('Spectral Coherence (30-day period)', fontsize=13, fontweight='bold')

    # Overall title
    fig.suptitle(f'Environmental Synchrony: Inner vs Outer Arc ({START_YEAR}-{END_YEAR})',
                fontsize=16, fontweight='bold', y=0.98)

    # Panel labels
    for i, label in enumerate(['A', 'B', 'C', 'D']):
        ax = [ax_inner, ax_outer, ax_pearson, ax_spec][i]
        ax.text(-0.08, 1.05, label, transform=ax.transAxes, fontsize=16, fontweight='bold')

    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"    Saved: {os.path.basename(output_path)}")

# ============================================
# Phase 9: Generate Figures and Save Results
# ============================================
print("\n" + "=" * 60)
print("PHASE 5: Generating Publication Figure and Saving Results")
print("=" * 60)

create_publication_figure(
    sst_arc, dli_arc, chl_arc,
    pearson_sync, spectral_sync,
    os.path.join(OUTPUT_DIR, 'Fig_Main_Temporal_Synchrony.png')
)

# Save CSVs
pearson_sync.to_csv(os.path.join(OUTPUT_DIR, 'Pearson_Synchrony_by_Site.csv'), index=False)
print("\n  Saved: Pearson_Synchrony_by_Site.csv")

spectral_sync.to_csv(os.path.join(OUTPUT_DIR, 'Spectral_Coherence_by_Site.csv'), index=False)
print("  Saved: Spectral_Coherence_by_Site.csv")

# Summary statistics
summary_data = []
for metric, label in [('Synchrony_Mean', 'Pearson_Synchrony'),
                       ('Coherence_7d', 'Coherence_7d'),
                       ('Coherence_30d', 'Coherence_30d')]:

    if label == 'Pearson_Synchrony':
        df = pearson_sync
    else:
        df = spectral_sync

    inner_vals = df[df['Arc'] == 'inner'][metric].dropna()
    outer_vals = df[df['Arc'] == 'outer'][metric].dropna()

    summary_data.append({
        'Metric': f'{label}_Inner',
        'Mean': inner_vals.mean() if len(inner_vals) > 0 else np.nan,
        'SD': inner_vals.std() if len(inner_vals) > 0 else np.nan,
        'Median': inner_vals.median() if len(inner_vals) > 0 else np.nan,
        'N': len(inner_vals)
    })

    summary_data.append({
        'Metric': f'{label}_Outer',
        'Mean': outer_vals.mean() if len(outer_vals) > 0 else np.nan,
        'SD': outer_vals.std() if len(outer_vals) > 0 else np.nan,
        'Median': outer_vals.median() if len(outer_vals) > 0 else np.nan,
        'N': len(outer_vals)
    })

summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(os.path.join(OUTPUT_DIR, 'Summary_Statistics.csv'), index=False)
print("  Saved: Summary_Statistics.csv")

# Print summary
print("\n" + "=" * 60)
print("SUMMARY STATISTICS")
print("=" * 60)
print(summary_df.to_string(index=False))

# ============================================
# Completion
# ============================================
print("\n" + "=" * 60)
print("ANALYSIS COMPLETE!")
print(f"All outputs saved to: {OUTPUT_DIR}")
print("=" * 60)
