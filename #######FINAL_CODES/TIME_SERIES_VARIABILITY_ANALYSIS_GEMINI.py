### TIME_SERIES_VARIABILITY_ANALYSIS_GEMINI.py ###
# Environmental Dynamics Analysis: Inner Arc vs Outer Arc Variability Comparison
# 
# This script tests the hypothesis that Inner Arc reefs exhibit higher 
# environmental variability than Outer Arc reefs in the Abrolhos region.
#
# Features:
# 1. CV Comparison Boxplots (sliding window: 7-day, 30-day, full series)
# 2. Lomb-Scargle Spectral Analysis (7-day, monthly, full spectrum)
# 3. Pulse Detection (extreme events > 90th percentile)
#
# All comparisons include Mann-Whitney U statistical tests.

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
    filename='variability_analysis_errors.log',
    filemode='w',
    level=logging.ERROR,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# --- Time Period ---
START_YEAR = 2002
END_YEAR = 2008

# --- Data Directories ---
SST_DIR = r'H:\remote sensing\CRW_SST_FULL'
MODIS_DIR = r'H:\remote sensing\MODIS_DATA_FULL'

# --- Sites CSV ---
SITES_CSV_PATH = r'C:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######22.04.23\DATA\sites_list_full_clean.csv'

# --- Output Directory ---
OUTPUT_DIR = r'C:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_RESULTS\Variability_Analysis_GEMINI'
os.makedirs(OUTPUT_DIR, exist_ok=True)

# --- Region of Interest ---
LAT_MIN, LAT_MAX = -18.5, -16.5
LON_MIN, LON_MAX = -39.5, -38.0

# --- File Patterns ---
SST_PATTERN = 'coraltemp_v3.1_*.nc'
KD490_PATTERN = 'AQUA_MODIS.*.L3m.DAY.KD.Kd_490.4km.nc'
PAR_PATTERN = 'AQUA_MODIS.*.L3m.DAY.PAR.par.4km.nc'
CHL_PATTERN = 'AQUA_MODIS.*.L3m.DAY.CHL.chlor_a.4km.nc'

# --- DLI Calculation Constants (Gattuso) ---
KDPAR_A, KDPAR_B, KDPAR_C = 0.0665, 0.874, 0.00121

# --- Color Palette ---
COLOR_INNER = '#CD853F'  # Light brown (Peru)
COLOR_OUTER = '#4169E1'  # Royal Blue
COLORS_ARC = {'inner': COLOR_INNER, 'outer': COLOR_OUTER}

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
print("ENVIRONMENTAL DYNAMICS ANALYSIS: Inner Arc vs Outer Arc")
print(f"Period: {START_YEAR} - {END_YEAR}")
print("=" * 60)

try:
    sites_df = pd.read_csv(SITES_CSV_PATH, sep=';')
    sites_df.columns = [col.strip() for col in sites_df.columns]
    
    # Standardize column names
    column_mapping = {
        'Reef_name': 'Reef_name',
        'Site_name': 'Site_name', 
        'HAB': 'HAB',
        'Latitude': 'Latitude',
        'Longitude': 'Longitude',
        'Depth_m': 'Depth_m',
        'Arc': 'Arc'
    }
    
    # Verify required columns
    required_cols = {'Site_name', 'Latitude', 'Longitude', 'Depth_m', 'Arc'}
    if not required_cols.issubset(sites_df.columns):
        missing = required_cols - set(sites_df.columns)
        raise ValueError(f"Missing columns: {missing}. Available: {list(sites_df.columns)}")
    
    # Clean Arc column and standardize
    sites_df['Arc'] = sites_df['Arc'].str.strip().str.lower()
    sites_df = sites_df[sites_df['Arc'].isin(['inner', 'outer'])].copy()
    
    print(f"\n✓ Loaded {len(sites_df)} sites")
    print(f"  Inner Arc: {(sites_df['Arc'] == 'inner').sum()} sites")
    print(f"  Outer Arc: {(sites_df['Arc'] == 'outer').sum()} sites")

except FileNotFoundError:
    print(f"ERROR: Sites file not found: {SITES_CSV_PATH}")
    exit(1)
except Exception as e:
    print(f"ERROR loading sites: {e}")
    exit(1)

# ============================================
# Utility Functions
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
            # Try with specific engine if available
            open_kwargs = {}
            if NETCDF_ENGINE:
                open_kwargs['engine'] = NETCDF_ENGINE
            
            with xr.open_dataset(f, **open_kwargs) as ds:
                var_name = next((v for v in var_names if v in ds.data_vars), None)
                if var_name is None:
                    continue
                
                data = ds[var_name]
                
                # Handle potential time dimension
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
    
    print(f"    ✓ Loaded data for {len(result)} sites")
    return result

def calculate_dli_for_sites(kd490_files, par_files, sites_df):
    """Calculate benthic DLI time series for each site."""
    print("\n  Calculating DLI...")
    
    # Create date-indexed dictionaries
    par_by_date = {parse_date_from_filename(f): f for f in par_files if parse_date_from_filename(f)}
    
    site_data = {row['Site_name']: {'times': [], 'values': []} 
                 for _, row in sites_df.iterrows()}
    
    total_files = len(kd490_files)
    success_count = 0
    error_count = 0
    
    # Prepare engine kwargs
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
                
                # Handle potential time dimension
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
                        
                        # Gattuso equation for Kd_PAR
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
    
    print(f"    ✓ Calculated DLI for {len(result)} sites")
    return result

# ============================================
# CV Calculation Functions
# ============================================

def calculate_cv_full(series):
    """Calculate CV for the entire series."""
    if len(series) < 10:
        return np.nan
    mean_val = series.mean()
    std_val = series.std()
    if mean_val == 0:
        return np.nan
    return (std_val / mean_val) * 100

def calculate_cv_rolling(series, window):
    """Calculate mean CV using sliding window approach."""
    if len(series) < window:
        return np.nan
    
    min_periods = max(2, int(window * 0.25))
    rolling_mean = series.rolling(window=window, min_periods=min_periods, center=True).mean()
    rolling_std = series.rolling(window=window, min_periods=min_periods, center=True).std()
    
    # Avoid division by zero
    cv_series = (rolling_std / (rolling_mean + 1e-9)) * 100
    cv_series = cv_series.replace([np.inf, -np.inf], np.nan)
    
    return cv_series.mean()

def calculate_all_cvs(timeseries_dict, sites_df):
    """Calculate CV for all sites at different frequency bands."""
    results = []
    
    for _, row in sites_df.iterrows():
        site_name = row['Site_name']
        arc = row['Arc']
        hab = row.get('HAB', 'Unknown')  # Include habitat for DLI analysis
        
        if site_name not in timeseries_dict:
            continue
        
        series = timeseries_dict[site_name]
        
        # Changed from 7-day to 2-day to match CV_02 from other analyses
        cv_2d = calculate_cv_rolling(series, 2)
        cv_30d = calculate_cv_rolling(series, 30)
        cv_full = calculate_cv_full(series)
        
        results.append({
            'Site_name': site_name,
            'Arc': arc,
            'HAB': hab,
            'CV_2d': cv_2d,
            'CV_30d': cv_30d,
            'CV_Full': cv_full
        })
    
    return pd.DataFrame(results)

# ============================================
# Lomb-Scargle Spectral Analysis
# ============================================

def compute_lombscargle(series, periods_of_interest=None):
    """Compute Lomb-Scargle periodogram for a time series."""
    if len(series) < 30:
        return None, None
    
    # Prepare data
    times = (series.index - series.index[0]).days.values.astype(float)
    values = series.values
    
    # Remove NaN
    mask = ~np.isnan(values)
    times = times[mask]
    values = values[mask]
    
    if len(times) < 30:
        return None, None
    
    # Normalize
    values = (values - values.mean()) / values.std()
    
    # Define frequency range (periods from 3 days to half the series length)
    max_period = len(times) // 2
    min_period = 3
    
    periods = np.linspace(min_period, max_period, 500)
    frequencies = 2 * np.pi / periods
    
    # Compute Lomb-Scargle
    try:
        power = lombscargle(times, values, frequencies, normalize=True)
        return periods, power
    except Exception as e:
        logging.error(f"Lomb-Scargle error: {e}")
        return None, None

def extract_power_at_period(periods, power, target_period, tolerance=0.3):
    """Extract power at a specific period with tolerance."""
    if periods is None or power is None:
        return np.nan
    
    mask = np.abs(periods - target_period) / target_period <= tolerance
    if mask.sum() == 0:
        return np.nan
    
    return power[mask].max()

def compute_spectral_power_for_sites(timeseries_dict, sites_df):
    """Compute spectral power at key periods for all sites."""
    results = []
    
    for _, row in sites_df.iterrows():
        site_name = row['Site_name']
        arc = row['Arc']
        
        if site_name not in timeseries_dict:
            continue
        
        series = timeseries_dict[site_name]
        periods, power = compute_lombscargle(series)
        
        power_7d = extract_power_at_period(periods, power, 7)
        power_30d = extract_power_at_period(periods, power, 30)
        
        results.append({
            'Site_name': site_name,
            'Arc': arc,
            'Power_7d': power_7d,
            'Power_30d': power_30d,
            'periods': periods,
            'power': power
        })
    
    return results

# ============================================
# Pulse Detection
# ============================================

def count_extreme_events(series, threshold):
    """Count number of values above threshold."""
    if len(series) < 10:
        return np.nan
    return (series > threshold).sum()

def calculate_pulse_counts(timeseries_dict, sites_df, percentile=90):
    """Calculate extreme event counts for all sites."""
    # First, compute the combined threshold from all data
    all_values = []
    for site_name in timeseries_dict:
        all_values.extend(timeseries_dict[site_name].dropna().values)
    
    if len(all_values) == 0:
        return pd.DataFrame()
    
    threshold = np.percentile(all_values, percentile)
    
    results = []
    for _, row in sites_df.iterrows():
        site_name = row['Site_name']
        arc = row['Arc']
        
        if site_name not in timeseries_dict:
            continue
        
        series = timeseries_dict[site_name]
        count = count_extreme_events(series, threshold)
        
        # Normalize by series length
        if len(series) > 0:
            freq = count / len(series) * 100  # Percentage of extreme events
        else:
            freq = np.nan
        
        results.append({
            'Site_name': site_name,
            'Arc': arc,
            'Extreme_Count': count,
            'Extreme_Freq': freq,
            'Threshold': threshold
        })
    
    return pd.DataFrame(results)

# ============================================
# Statistical Testing
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
# Plotting Functions
# ============================================

def create_cv_boxplot_figure(cv_sst, cv_dli, cv_chl, output_path):
    """Create 3x3 CV comparison boxplot figure."""
    print("\n  Creating CV Comparison Figure...")
    
    # Check if we have data
    if len(cv_sst) == 0 and len(cv_dli) == 0 and len(cv_chl) == 0:
        print("    WARNING: No CV data available. Skipping figure.")
        return
    
    fig, axes = plt.subplots(3, 3, figsize=(14, 12))
    fig.suptitle('Coefficient of Variation: Inner Arc vs Outer Arc', fontsize=16, fontweight='bold', y=0.98)
    
    variables = [
        ('SST', cv_sst, '°C'),
        ('DLI', cv_dli, 'mol m⁻² d⁻¹'),
        ('Chl-a', cv_chl, 'mg m⁻³')
    ]
    # Updated: 7-day → 2-day to match CV_02 from other analyses
    bands = ['CV_2d', 'CV_30d', 'CV_Full']
    band_labels = ['2-day Window', '30-day Window', 'Full Series']
    
    for row, (var_name, df, unit) in enumerate(variables):
        for col, (band, band_label) in enumerate(zip(bands, band_labels)):
            ax = axes[row, col]
            
            if len(df) == 0 or 'Arc' not in df.columns or band not in df.columns:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                continue
            
            inner_vals = df[df['Arc'] == 'inner'][band].dropna()
            outer_vals = df[df['Arc'] == 'outer'][band].dropna()
            
            if len(inner_vals) == 0 and len(outer_vals) == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                continue
            
            # Create boxplot data
            bp_data = [inner_vals.values if len(inner_vals) > 0 else [0], 
                       outer_vals.values if len(outer_vals) > 0 else [0]]
            bp = ax.boxplot(bp_data, labels=['Inner Arc', 'Outer Arc'], 
                           patch_artist=True, widths=0.6)
            
            # Style boxes
            bp['boxes'][0].set_facecolor(COLOR_INNER)
            bp['boxes'][1].set_facecolor(COLOR_OUTER)
            for box in bp['boxes']:
                box.set_alpha(0.7)
            
            # Add individual points
            for i, (vals, color) in enumerate([(inner_vals, COLOR_INNER), (outer_vals, COLOR_OUTER)]):
                if len(vals) > 0:
                    x = np.random.normal(i + 1, 0.04, len(vals))
                    ax.scatter(x, vals, alpha=0.5, s=30, c=color, edgecolors='white', linewidth=0.5, zorder=3)
            
            # Y-axis auto-scaling with padding
            all_vals = np.concatenate([inner_vals.values, outer_vals.values]) if len(inner_vals) > 0 and len(outer_vals) > 0 else (inner_vals.values if len(inner_vals) > 0 else outer_vals.values)
            if len(all_vals) > 0:
                ymin, ymax = np.nanmin(all_vals), np.nanmax(all_vals)
                yrange = ymax - ymin
                if yrange > 0:
                    ax.set_ylim(ymin - 0.1 * yrange, ymax + 0.15 * yrange)
            
            # Statistical test
            p_val, p_str = mann_whitney_test(inner_vals, outer_vals)
            sig = format_significance(p_val)
            
            # Add p-value to plot
            ax.text(0.95, 0.95, f'{p_str}\n{sig}', transform=ax.transAxes, 
                   ha='right', va='top', fontsize=9, 
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            
            # Labels
            if row == 0:
                ax.set_title(band_label, fontsize=12, fontweight='bold')
            if col == 0:
                ax.set_ylabel(f'{var_name} CV (%)', fontsize=11)
            if row == 2:
                ax.set_xlabel('')
            
            ax.tick_params(axis='x', rotation=0)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"    ✓ Saved: {os.path.basename(output_path)}")

def create_spectral_power_figure(spectral_sst, spectral_dli, spectral_chl, sites_df, output_path):
    """Create spectral power comparison figure."""
    print("\n  Creating Spectral Power Figure...")
    
    # Check if we have any data
    if len(spectral_sst) == 0 and len(spectral_dli) == 0 and len(spectral_chl) == 0:
        print("    WARNING: No spectral data available. Skipping figure.")
        return
    
    fig, axes = plt.subplots(2, 3, figsize=(14, 8))
    fig.suptitle('Lomb-Scargle Spectral Power: Inner Arc vs Outer Arc', fontsize=16, fontweight='bold', y=0.98)
    
    datasets = [
        ('SST', spectral_sst),
        ('DLI', spectral_dli),
        ('Chl-a', spectral_chl)
    ]
    periods_to_plot = [('7d', 7, 'Power_7d'), ('30d', 30, 'Power_30d')]
    
    for col, (var_name, spectral_data) in enumerate(datasets):
        # Handle empty spectral_data
        if len(spectral_data) == 0:
            for row in range(2):
                axes[row, col].text(0.5, 0.5, 'No data', ha='center', va='center', transform=axes[row, col].transAxes)
                if row == 0:
                    axes[row, col].set_title(var_name, fontsize=12, fontweight='bold')
            continue
        
        # Convert to DataFrame for easier manipulation
        df = pd.DataFrame([{
            'Site_name': d['Site_name'],
            'Arc': d['Arc'],
            'Power_7d': d['Power_7d'],
            'Power_30d': d['Power_30d']
        } for d in spectral_data])
        
        for row, (period_label, period_val, col_name) in enumerate(periods_to_plot):
            ax = axes[row, col]
            
            if 'Arc' not in df.columns or len(df) == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                continue
            
            inner_vals = df[df['Arc'] == 'inner'][col_name].dropna()
            outer_vals = df[df['Arc'] == 'outer'][col_name].dropna()
            
            if len(inner_vals) == 0 and len(outer_vals) == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                continue
            
            bp_data = [inner_vals.values if len(inner_vals) > 0 else [0], 
                       outer_vals.values if len(outer_vals) > 0 else [0]]
            bp = ax.boxplot(bp_data, labels=['Inner', 'Outer'], 
                           patch_artist=True, widths=0.6)
            
            bp['boxes'][0].set_facecolor(COLOR_INNER)
            bp['boxes'][1].set_facecolor(COLOR_OUTER)
            for box in bp['boxes']:
                box.set_alpha(0.7)
            
            # Add points
            for i, (vals, color) in enumerate([(inner_vals, COLOR_INNER), (outer_vals, COLOR_OUTER)]):
                if len(vals) > 0:
                    x = np.random.normal(i + 1, 0.04, len(vals))
                    ax.scatter(x, vals, alpha=0.5, s=30, c=color, edgecolors='white', linewidth=0.5, zorder=3)
            
            # Y-axis auto-scaling with padding
            all_vals = np.concatenate([inner_vals.values, outer_vals.values]) if len(inner_vals) > 0 and len(outer_vals) > 0 else (inner_vals.values if len(inner_vals) > 0 else outer_vals.values)
            if len(all_vals) > 0:
                ymin, ymax = np.nanmin(all_vals), np.nanmax(all_vals)
                yrange = ymax - ymin
                if yrange > 0:
                    ax.set_ylim(ymin - 0.1 * yrange, ymax + 0.15 * yrange)
            
            # Statistical test
            p_val, p_str = mann_whitney_test(inner_vals, outer_vals)
            sig = format_significance(p_val)
            ax.text(0.95, 0.95, f'{p_str}\n{sig}', transform=ax.transAxes,
                   ha='right', va='top', fontsize=9,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
            
            if row == 0:
                ax.set_title(var_name, fontsize=12, fontweight='bold')
            ax.set_ylabel(f'{period_label} Power', fontsize=10)
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"    ✓ Saved: {os.path.basename(output_path)}")

def create_periodogram_heatmaps(spectral_sst, spectral_dli, spectral_chl, output_path):
    """Create full spectrum heatmap figure."""
    print("\n  Creating Periodogram Heatmaps...")
    
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle('Lomb-Scargle Periodograms by Arc', fontsize=16, fontweight='bold', y=0.98)
    
    datasets = [
        ('SST', spectral_sst),
        ('DLI', spectral_dli),
        ('Chl-a', spectral_chl)
    ]
    
    for col, (var_name, spectral_data) in enumerate(datasets):
        for row, arc in enumerate(['inner', 'outer']):
            ax = axes[row, col]
            
            # Get data for this arc
            arc_data = [d for d in spectral_data if d['Arc'] == arc and d['periods'] is not None]
            
            if len(arc_data) == 0:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                continue
            
            # Create matrix for heatmap (sites x periods)
            # Use common period grid
            common_periods = arc_data[0]['periods']
            power_matrix = []
            site_names = []
            
            for d in arc_data:
                if d['power'] is not None:
                    power_matrix.append(d['power'])
                    site_names.append(d['Site_name'])
            
            if len(power_matrix) == 0:
                continue
            
            power_matrix = np.array(power_matrix)
            
            # Plot mean spectrum with envelope
            mean_power = np.nanmean(power_matrix, axis=0)
            std_power = np.nanstd(power_matrix, axis=0)
            
            color = COLOR_INNER if arc == 'inner' else COLOR_OUTER
            ax.fill_between(common_periods, mean_power - std_power, mean_power + std_power, 
                           alpha=0.3, color=color)
            ax.plot(common_periods, mean_power, color=color, linewidth=2, 
                   label=f'{arc.capitalize()} Arc (n={len(site_names)})')
            
            # Mark key periods
            ax.axvline(7, color='red', linestyle='--', alpha=0.5, label='7-day')
            ax.axvline(30, color='green', linestyle='--', alpha=0.5, label='30-day')
            ax.axvline(365, color='purple', linestyle='--', alpha=0.5, label='Annual')
            
            # Y-axis auto-scaling based on UPPER ENVELOPE in the VISIBLE range
            # This prevents both overflow (lines exceeding plot area) and compression
            if len(mean_power) > 0:
                upper_envelope = mean_power + std_power
                # Use full visible range for scaling (3 to 365 days)
                visible_mask = (common_periods >= 3) & (common_periods <= 365)
                if visible_mask.any():
                    ymax_envelope = np.nanmax(upper_envelope[visible_mask])
                    if ymax_envelope > 0:
                        ax.set_ylim(0, ymax_envelope * 1.1)  # 10% padding above envelope max
            
            ax.set_xlim(3, 365)  # Extended to show annual cycles
            ax.set_xlabel('Period (days)', fontsize=10)
            ax.set_ylabel('Normalized Power', fontsize=10)
            
            if row == 0:
                ax.set_title(var_name, fontsize=12, fontweight='bold')
            
            if col == 0:
                arc_label = 'Inner Arc' if arc == 'inner' else 'Outer Arc'
                ax.text(-0.25, 0.5, arc_label, transform=ax.transAxes, 
                       fontsize=12, fontweight='bold', rotation=90, va='center')
            
            ax.legend(loc='upper right', fontsize=8)
    
    plt.tight_layout(rect=[0.05, 0, 1, 0.96])
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"    ✓ Saved: {os.path.basename(output_path)}")

def create_pulse_detection_figure(pulse_sst, pulse_dli, pulse_chl, output_path):
    """Create extreme events frequency figure."""
    print("\n  Creating Pulse Detection Figure...")
    
    # Check if we have any data
    if len(pulse_sst) == 0 and len(pulse_dli) == 0 and len(pulse_chl) == 0:
        print("    WARNING: No pulse data available. Skipping figure.")
        return
    
    fig, axes = plt.subplots(1, 3, figsize=(14, 5))
    fig.suptitle('Extreme Events Frequency (>90th Percentile): Inner Arc vs Outer Arc', 
                fontsize=14, fontweight='bold', y=1.02)
    
    datasets = [
        ('SST', pulse_sst),
        ('DLI', pulse_dli),
        ('Chl-a', pulse_chl)
    ]
    
    for col, (var_name, df) in enumerate(datasets):
        ax = axes[col]
        
        if len(df) == 0 or 'Arc' not in df.columns:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(var_name, fontsize=12, fontweight='bold')
            continue
        
        inner_vals = df[df['Arc'] == 'inner']['Extreme_Freq'].dropna()
        outer_vals = df[df['Arc'] == 'outer']['Extreme_Freq'].dropna()
        
        if len(inner_vals) == 0 and len(outer_vals) == 0:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(var_name, fontsize=12, fontweight='bold')
            continue
        
        bp_data = [inner_vals.values if len(inner_vals) > 0 else [0], 
                   outer_vals.values if len(outer_vals) > 0 else [0]]
        bp = ax.boxplot(bp_data, labels=['Inner Arc', 'Outer Arc'], 
                       patch_artist=True, widths=0.6)
        
        bp['boxes'][0].set_facecolor(COLOR_INNER)
        bp['boxes'][1].set_facecolor(COLOR_OUTER)
        for box in bp['boxes']:
            box.set_alpha(0.7)
        
        # Add points
        for i, (vals, color) in enumerate([(inner_vals, COLOR_INNER), (outer_vals, COLOR_OUTER)]):
            if len(vals) > 0:
                x = np.random.normal(i + 1, 0.04, len(vals))
                ax.scatter(x, vals, alpha=0.5, s=40, c=color, edgecolors='white', linewidth=0.5, zorder=3)
        
        # Y-axis auto-scaling with padding
        all_vals = np.concatenate([inner_vals.values, outer_vals.values]) if len(inner_vals) > 0 and len(outer_vals) > 0 else (inner_vals.values if len(inner_vals) > 0 else outer_vals.values)
        if len(all_vals) > 0:
            ymin, ymax = np.nanmin(all_vals), np.nanmax(all_vals)
            yrange = ymax - ymin
            if yrange > 0:
                ax.set_ylim(ymin - 0.1 * yrange, ymax + 0.15 * yrange)
        
        # Statistical test
        p_val, p_str = mann_whitney_test(inner_vals, outer_vals)
        sig = format_significance(p_val)
        ax.text(0.95, 0.95, f'{p_str}\n{sig}', transform=ax.transAxes,
               ha='right', va='top', fontsize=10,
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax.set_title(var_name, fontsize=12, fontweight='bold')
        ax.set_ylabel('Extreme Events (%)', fontsize=11)
        
        # Add threshold info
        if len(df) > 0 and 'Threshold' in df.columns:
            threshold = df['Threshold'].iloc[0]
            ax.text(0.05, 0.95, f'Threshold: {threshold:.2f}', transform=ax.transAxes,
                   ha='left', va='top', fontsize=8, style='italic')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"    ✓ Saved: {os.path.basename(output_path)}")

def create_dli_habitat_figure(cv_dli, output_path):
    """Create DLI comparison figure by Arc and Habitat."""
    print("\n  Creating DLI Habitat Comparison Figure...")
    
    if len(cv_dli) == 0 or 'HAB' not in cv_dli.columns:
        print("    WARNING: No DLI habitat data available. Skipping figure.")
        return
    
    # Get available habitats
    available_habs = cv_dli['HAB'].unique()
    print(f"    Available habitats: {list(available_habs)}")
    
    fig, axes = plt.subplots(1, 3, figsize=(16, 6))
    fig.suptitle('DLI Coefficient of Variation by Arc and Habitat', fontsize=16, fontweight='bold', y=1.02)
    
    bands = ['CV_2d', 'CV_30d', 'CV_Full']
    band_labels = ['2-day Window', '30-day Window', 'Full Series']
    
    # Standard habitat colors (from Boxplot_Cover_RGR_GLM.R)
    hab_order = ['TP', 'PA', 'RR']
    hab_colors = {
        'PA': '#E69F00',  # Orange/Gold
        'RR': '#56B4E9',  # Sky Blue
        'TP': '#009E73'   # Greenish
    }
    
    for col, (band, band_label) in enumerate(zip(bands, band_labels)):
        ax = axes[col]
        
        if band not in cv_dli.columns:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(band_label, fontsize=12, fontweight='bold')
            continue
        
        # Prepare data for grouped boxplot
        plot_data = []
        plot_labels = []
        plot_colors = []
        positions = []
        pos = 0
        group_positions = []
        
        for arc in ['inner', 'outer']:
            arc_data = cv_dli[cv_dli['Arc'] == arc]
            arc_start = pos
            
            for hab in hab_order:
                hab_data = arc_data[arc_data['HAB'] == hab][band].dropna()
                if len(hab_data) > 0:
                    plot_data.append(hab_data.values)
                    plot_labels.append(f'{hab}')
                    plot_colors.append(hab_colors.get(hab, 'gray'))
                    positions.append(pos)
                    pos += 1
            
            if pos > arc_start:
                group_positions.append((arc_start + pos - 1) / 2)  # Center of group
            pos += 0.5  # Gap between arcs
        
        if len(plot_data) == 0:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(band_label, fontsize=12, fontweight='bold')
            continue
        
        # Create boxplot
        bp = ax.boxplot(plot_data, positions=positions, widths=0.6, patch_artist=True)
        
        # Style boxes
        for patch, color in zip(bp['boxes'], plot_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        
        # Add individual points
        for i, (data, pos_val) in enumerate(zip(plot_data, positions)):
            x = np.random.normal(pos_val, 0.05, len(data))
            ax.scatter(x, data, alpha=0.5, s=25, c=plot_colors[i], edgecolors='white', linewidth=0.5, zorder=3)
        
        # Y-axis auto-scaling
        all_vals = np.concatenate(plot_data)
        if len(all_vals) > 0:
            ymin, ymax = np.nanmin(all_vals), np.nanmax(all_vals)
            yrange = ymax - ymin
            if yrange > 0:
                ax.set_ylim(ymin - 0.1 * yrange, ymax + 0.15 * yrange)
        
        # Set x-axis
        ax.set_xticks(positions)
        ax.set_xticklabels(plot_labels, fontsize=9)
        
        # Add arc labels below
        if len(group_positions) >= 2:
            ax.text(group_positions[0], ax.get_ylim()[0] - 0.15 * yrange, 'Inner Arc', 
                   ha='center', fontsize=11, fontweight='bold', color=COLOR_INNER)
            ax.text(group_positions[1], ax.get_ylim()[0] - 0.15 * yrange, 'Outer Arc', 
                   ha='center', fontsize=11, fontweight='bold', color=COLOR_OUTER)
        
        ax.set_title(band_label, fontsize=12, fontweight='bold')
        if col == 0:
            ax.set_ylabel('DLI CV (%)', fontsize=11)
    
    # Add habitat legend
    legend_elements = [plt.Rectangle((0, 0), 1, 1, facecolor=hab_colors[h], alpha=0.7, label=h) 
                       for h in hab_order if h in available_habs]
    fig.legend(handles=legend_elements, title='Habitat', loc='upper right', 
               bbox_to_anchor=(0.98, 0.95), fontsize=10)
    
    plt.tight_layout(rect=[0, 0.05, 0.92, 0.95])
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"    ✓ Saved: {os.path.basename(output_path)}")

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
    print("STEP 2: Calculating Coefficients of Variation")
    print("=" * 60)
    
    cv_sst = calculate_all_cvs(sst_ts, sites_df)
    cv_dli = calculate_all_cvs(dli_ts, sites_df)
    cv_chl = calculate_all_cvs(chl_ts, sites_df)
    
    print(f"\n  ✓ CV calculated for SST: {len(cv_sst)} sites")
    print(f"  ✓ CV calculated for DLI: {len(cv_dli)} sites")
    print(f"  ✓ CV calculated for Chl-a: {len(cv_chl)} sites")
    
    print("\n" + "=" * 60)
    print("STEP 3: Computing Lomb-Scargle Spectral Analysis")
    print("=" * 60)
    
    spectral_sst = compute_spectral_power_for_sites(sst_ts, sites_df)
    spectral_dli = compute_spectral_power_for_sites(dli_ts, sites_df)
    spectral_chl = compute_spectral_power_for_sites(chl_ts, sites_df)
    
    print(f"\n  ✓ Spectral analysis for SST: {len(spectral_sst)} sites")
    print(f"  ✓ Spectral analysis for DLI: {len(spectral_dli)} sites")
    print(f"  ✓ Spectral analysis for Chl-a: {len(spectral_chl)} sites")
    
    print("\n" + "=" * 60)
    print("STEP 4: Detecting Extreme Events (Pulses)")
    print("=" * 60)
    
    pulse_sst = calculate_pulse_counts(sst_ts, sites_df)
    pulse_dli = calculate_pulse_counts(dli_ts, sites_df)
    pulse_chl = calculate_pulse_counts(chl_ts, sites_df)
    
    print(f"\n  ✓ Pulse detection for SST: {len(pulse_sst)} sites")
    print(f"  ✓ Pulse detection for DLI: {len(pulse_dli)} sites")
    print(f"  ✓ Pulse detection for Chl-a: {len(pulse_chl)} sites")
    
    print("\n" + "=" * 60)
    print("STEP 5: Generating Publication-Quality Figures")
    print("=" * 60)
    
    # Figure 1: CV Comparison Boxplots
    create_cv_boxplot_figure(
        cv_sst, cv_dli, cv_chl,
        os.path.join(OUTPUT_DIR, 'Fig1_CV_Comparison_Boxplots.png')
    )
    
    # Figure 2: Spectral Power Comparison
    create_spectral_power_figure(
        spectral_sst, spectral_dli, spectral_chl, sites_df,
        os.path.join(OUTPUT_DIR, 'Fig2_Spectral_Power_Comparison.png')
    )
    
    # Figure 3: Periodogram Heatmaps
    create_periodogram_heatmaps(
        spectral_sst, spectral_dli, spectral_chl,
        os.path.join(OUTPUT_DIR, 'Fig3_Periodogram_Heatmaps.png')
    )
    
    # Figure 4: Extreme Events Frequency
    create_pulse_detection_figure(
        pulse_sst, pulse_dli, pulse_chl,
        os.path.join(OUTPUT_DIR, 'Fig4_Extreme_Events_Frequency.png')
    )
    
    # Figure 5: DLI by Habitat (NEW)
    create_dli_habitat_figure(
        cv_dli,
        os.path.join(OUTPUT_DIR, 'Fig5_DLI_Habitat_Comparison.png')
    )
    
    print("\n" + "=" * 60)
    print("STEP 6: Saving Summary Statistics")
    print("=" * 60)
    
    # Save CV summary
    cv_summary = pd.concat([
        cv_sst.assign(Variable='SST'),
        cv_dli.assign(Variable='DLI'),
        cv_chl.assign(Variable='Chl-a')
    ])
    cv_summary.to_csv(os.path.join(OUTPUT_DIR, 'CV_Summary_by_Site.csv'), index=False)
    print("  ✓ Saved: CV_Summary_by_Site.csv")
    
    # Save pulse summary
    pulse_summary = pd.concat([
        pulse_sst.assign(Variable='SST'),
        pulse_dli.assign(Variable='DLI'),
        pulse_chl.assign(Variable='Chl-a')
    ])
    pulse_summary.to_csv(os.path.join(OUTPUT_DIR, 'Pulse_Summary_by_Site.csv'), index=False)
    print("  ✓ Saved: Pulse_Summary_by_Site.csv")
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE!")
    print(f"All outputs saved to: {OUTPUT_DIR}")
    print("=" * 60)
