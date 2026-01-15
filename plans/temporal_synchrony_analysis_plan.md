# Plano de Implementação: Temporal Synchrony Analysis
## Análise de Sincronicidade Temporal entre SST, DLI e Clorofila (Inner vs Outer Arc)

**Objetivo**: Produzir figura de série temporal em nível Nature/Science e avaliar quantitativamente se o acoplamento (sincronicidade) entre SST, DLI e Clorofila é maior no arco interno vs externo.

**Período**: 2002-2008
**Referência de estilo**: `TIME_SERIES_VARIABILITY_ANALYSIS_GEMINI.py`

---

## 1. Estrutura Geral do Script

```python
### TIME_SERIES_SYNCHRONY_ANALYSIS.py ###
# Temporal Synchrony Analysis: Inner Arc vs Outer Arc Environmental Coupling
#
# This script evaluates whether Inner Arc reefs exhibit higher temporal
# synchrony (coupling) between SST, DLI, and Chl-a than Outer Arc reefs.
#
# Features:
# 1. Publication-quality time series (smooth + raw)
# 2. Pearson correlation analysis for synchrony
# 3. Spectral coherence analysis (Lomb-Scargle)
# 4. Statistical comparison (Mann-Whitney U)
```

---

## 2. Configuração Inicial (Baseado em GEMINI.py)

### 2.1 Imports e Logging

```python
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
COLOR_SST = '#D55E00'    # Vermelho profissional (distinto de Inner/Outer)
COLOR_DLI = '#0072B2'    # Azul profissional
COLOR_CHL = '#009E73'    # Verde profissional

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
```

### 2.2 Carregamento de Sites (Padrão GEMINI.py)

```python
def load_sites_data():
    """Load and standardize sites data."""
    try:
        sites_df = pd.read_csv(SITES_CSV_PATH, sep=';')
        sites_df.columns = [col.strip() for col in sites_df.columns]

        # Verify required columns
        required_cols = {'Site_name', 'Latitude', 'Longitude', 'Depth_m', 'Arc'}
        if not required_cols.issubset(sites_df.columns):
            missing = required_cols - set(sites_df.columns)
            raise ValueError(f"Missing columns: {missing}")

        # Clean and standardize
        sites_df['Arc'] = sites_df['Arc'].str.strip().str.lower()
        sites_df = sites_df[sites_df['Arc'].isin(['inner', 'outer'])].copy()

        print(f"✓ Loaded {len(sites_df)} sites")
        print(f"  Inner Arc: {(sites_df['Arc'] == 'inner').sum()} sites")
        print(f"  Outer Arc: {(sites_df['Arc'] == 'outer').sum()} sites")

        return sites_df

    except Exception as e:
        print(f"ERROR loading sites: {e}")
        exit(1)
```

---

## 3. Funções de Carregamento de Dados

### 3.1 Load Time Series (Adaptado de GEMINI.py)

```python
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
    """Load time series data for all sites."""
    print(f"\n  Loading {var_label}...")

    site_data = {row['Site_name']: {'times': [], 'values': []}
                 for _, row in sites_df.iterrows()}

    total_files = len(files)

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

        except Exception as e:
            logging.error(f"Error reading {f}: {e}")
            continue

    # Convert to pandas Series
    result = {}
    for site_name, data in site_data.items():
        if len(data['times']) > 0:
            result[site_name] = pd.Series(data['values'], index=pd.DatetimeIndex(data['times'])).sort_index()

    print(f"    ✓ Loaded data for {len(result)} sites")
    return result
```

### 3.2 Calculate DLI (Adaptado de GEMINI.py)

```python
def calculate_dli_for_sites(kd490_files, par_files, sites_df):
    """Calculate benthic DLI time series for each site."""
    print("\n  Calculating DLI...")

    par_by_date = {parse_date_from_filename(f): f for f in par_files if parse_date_from_filename(f)}

    site_data = {row['Site_name']: {'times': [], 'values': []}
                 for _, row in sites_df.iterrows()}

    total_files = len(kd490_files)

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

        except Exception as e:
            logging.error(f"Error calculating DLI: {e}")
            continue

    # Convert to pandas Series
    result = {}
    for site_name, data in site_data.items():
        if len(data['times']) > 0:
            result[site_name] = pd.Series(data['values'], index=pd.DatetimeIndex(data['times'])).sort_index()

    print(f"    ✓ Calculated DLI for {len(result)} sites")
    return result
```

### 3.3 Agregar Séries Temporais por Arco

```python
def aggregate_timeseries_by_arc(timeseries_dict, sites_df):
    """Aggregate time series by arc (mean of all sites in each arc)."""

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

    return result
```

---

## 4. Análise de Sincronicidade

### 4.1 Métrica 1: Correlação de Pearson

```python
def calculate_pearson_synchrony(sst_dict, dli_dict, chl_dict, sites_df):
    """
    Calculate Pearson correlation synchrony for each site.

    Returns DataFrame with synchrony metrics per site.
    """
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

        mean_z = np.mean([z_sst_dli, z_sst_chl, z_dli_chl])
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

    return pd.DataFrame(results)
```

### 4.2 Métrica 2: Coerência Espectral (Lomb-Scargle)

```python
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
    s1_norm = (df['s1'] - df['s1'].mean()) / df['s1'].std()
    s2_norm = (df['s2'] - df['s2'].mean()) / df['s2'].std()

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

    except Exception:
        return np.nan

def calculate_spectral_synchrony(sst_dict, dli_dict, chl_dict, sites_df):
    """
    Calculate spectral coherence synchrony for each site.

    Evaluates synchrony at key periods: 7 days, 30 days.
    """
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
        if not np.isnan([coh_7d_sst_dli, coh_7d_sst_chl, coh_7d_dli_chl]).all():
            mean_coh_7d = np.nanmean([coh_7d_sst_dli, coh_7d_sst_chl, coh_7d_dli_chl])
        else:
            mean_coh_7d = np.nan

        if not np.isnan([coh_30d_sst_dli, coh_30d_sst_chl, coh_30d_dli_chl]).all():
            mean_coh_30d = np.nanmean([coh_30d_sst_dli, coh_30d_sst_chl, coh_30d_dli_chl])
        else:
            mean_coh_30d = np.nan

        results.append({
            'Site_name': site_name,
            'Arc': arc,
            'Coherence_7d': mean_coh_7d,
            'Coherence_30d': mean_coh_30d
        })

    return pd.DataFrame(results)
```

### 4.3 Teste Estatístico (Mann-Whitney U)

```python
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
        z_score = stats.norm.ppf(1 - p_value/2)  # Two-tailed z
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

    except Exception:
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
```

---

## 5. Funções de Plotagem (Nível Nature/Science)

### 5.1 Suavização de Série Temporal

```python
def smooth_timeseries(series, window_days=7, method='rolling'):
    """
    Apply smoothing to time series.

    Methods:
    - 'rolling': Rolling mean with centered window
    - 'savgol': Savitzky-Golay filter (requires uniform sampling)
    - 'lowess': LOWESS smoothing
    """
    if method == 'rolling':
        min_periods = max(1, int(window_days * 0.5))
        smoothed = series.rolling(window=window_days, min_periods=min_periods, center=True).mean()
        return smoothed

    elif method == 'lowess':
        from statsmodels.nonparametric.smoothers_lowess import lowess
        # Convert to numeric for lowess
        x = np.arange(len(series))
        y = series.values
        smoothed = lowess(y, x, frac=window_days/len(series), return_sorted=False)
        return pd.Series(smoothed, index=series.index)

    else:
        return series
```

### 5.2 Figura Principal Multi-Painel (2x2)

```python
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
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.25, top=0.92, bottom=0.08, left=0.08, right=0.95)

    # ===============================
    # Panel A: Inner Arc Time Series
    # ===============================
    ax_inner = fig.add_subplot(gs[0, 0])

    # Get data
    sst_i = sst_arc['inner']
    dli_i = dli_arc['inner']
    chl_i = chl_arc['inner']

    # Smooth
    sst_i_smooth = smooth_timeseries(sst_i, window_days=14)
    dli_i_smooth = smooth_timeseries(dli_i, window_days=14)
    chl_i_smooth = smooth_timeseries(chl_i, window_days=14)

    # Plot raw data (points, transparent)
    ax_inner.scatter(sst_i.index, sst_i.values, s=8, c=COLOR_SST, alpha=0.15, label='_nolegend_')
    ax_inner.scatter(dli_i.index, dli_i.values, s=8, c=COLOR_DLI, alpha=0.15, label='_nolegend_')
    ax_inner.scatter(chl_i.index, chl_i.values, s=8, c=COLOR_CHL, alpha=0.15, label='_nolegend_')

    # Plot smooth lines
    line_sst, = ax_inner.plot(sst_i_smooth.index, sst_i_smooth.values, '-', lw=1.5, c=COLOR_SST, label='SST')
    line_dli, = ax_inner.plot(dli_i_smooth.index, dli_i_smooth.values, '-', lw=1.5, c=COLOR_DLI, label='DLI')
    line_chl, = ax_inner.plot(chl_i_smooth.index, chl_i_smooth.values, '-', lw=1.5, c=COLOR_CHL, label='Chl-a')

    # Twin axes for proper scaling
    ax_inner_dli = ax_inner.twinx()
    ax_inner_chl = ax_inner.twinx()
    ax_inner_chl.spines['right'].set_position(('outward', 60))

    # Replot DLI and Chl on their axes
    ax_inner_dli.plot(dli_i_smooth.index, dli_i_smooth.values, '-', lw=1.5, c=COLOR_DLI, label='_nolegend_')
    ax_inner_chl.plot(chl_i_smooth.index, chl_i_smooth.values, '-', lw=1.5, c=COLOR_CHL, label='_nolegend_')

    # Styling
    ax_inner.set_ylabel('SST (°C)', color=COLOR_SST, fontsize=12, fontweight='bold')
    ax_inner_dli.set_ylabel('DLI (mol m⁻² d⁻¹)', color=COLOR_DLI, fontsize=12, fontweight='bold')
    ax_inner_chl.set_ylabel('Chl-a (mg m⁻³)', color=COLOR_CHL, fontsize=12, fontweight='bold')

    ax_inner.tick_params(axis='y', labelcolor=COLOR_SST, labelsize=10)
    ax_inner_dli.tick_params(axis='y', labelcolor=COLOR_DLI, labelsize=10)
    ax_inner_chl.tick_params(axis='y', labelcolor=COLOR_CHL, labelsize=10)

    ax_inner.set_title('Inner Arc', fontsize=14, fontweight='bold', loc='left', pad=10)
    ax_inner.set_xlim(pd.Timestamp(f'{START_YEAR}-01-01'), pd.Timestamp(f'{END_YEAR}-12-31'))

    # Legend
    lines = [line_sst, line_dli, line_chl]
    ax_inner.legend(lines, ['SST', 'DLI', 'Chl-a'], loc='upper left', fontsize=9,
                    framealpha=0.9, edgecolor='none')

    # ===============================
    # Panel B: Outer Arc Time Series
    # ===============================
    ax_outer = fig.add_subplot(gs[0, 1])

    # Same structure as Panel A, but for Outer Arc
    sst_o = sst_arc['outer']
    dli_o = dli_arc['outer']
    chl_o = chl_arc['outer']

    sst_o_smooth = smooth_timeseries(sst_o, window_days=14)
    dli_o_smooth = smooth_timeseries(dli_o, window_days=14)
    chl_o_smooth = smooth_timeseries(chl_o, window_days=14)

    # Raw points
    ax_outer.scatter(sst_o.index, sst_o.values, s=8, c=COLOR_SST, alpha=0.15)
    ax_outer.scatter(dli_o.index, dli_o.values, s=8, c=COLOR_DLI, alpha=0.15)
    ax_outer.scatter(chl_o.index, chl_o.values, s=8, c=COLOR_CHL, alpha=0.15)

    # Smooth lines
    ax_outer.plot(sst_o_smooth.index, sst_o_smooth.values, '-', lw=1.5, c=COLOR_SST)
    ax_outer.plot(dli_o_smooth.index, dli_o_smooth.values, '-', lw=1.5, c=COLOR_DLI)
    ax_outer.plot(chl_o_smooth.index, chl_o_smooth.values, '-', lw=1.5, c=COLOR_CHL)

    # Twin axes
    ax_outer_dli = ax_outer.twinx()
    ax_outer_chl = ax_outer.twinx()
    ax_outer_chl.spines['right'].set_position(('outward', 60))

    ax_outer_dli.plot(dli_o_smooth.index, dli_o_smooth.values, '-', lw=1.5, c=COLOR_DLI)
    ax_outer_chl.plot(chl_o_smooth.index, chl_o_smooth.values, '-', lw=1.5, c=COLOR_CHL)

    # Styling
    ax_outer.set_ylabel('SST (°C)', color=COLOR_SST, fontsize=12, fontweight='bold')
    ax_outer_dli.set_ylabel('DLI (mol m⁻² d⁻¹)', color=COLOR_DLI, fontsize=12, fontweight='bold')
    ax_outer_chl.set_ylabel('Chl-a (mg m⁻³)', color=COLOR_CHL, fontsize=12, fontweight='bold')

    ax_outer.tick_params(axis='y', labelcolor=COLOR_SST, labelsize=10)
    ax_outer_dli.tick_params(axis='y', labelcolor=COLOR_DLI, labelsize=10)
    ax_outer_chl.tick_params(axis='y', labelcolor=COLOR_CHL, labelsize=10)

    ax_outer.set_title('Outer Arc', fontsize=14, fontweight='bold', loc='left', pad=10)
    ax_outer.set_xlim(pd.Timestamp(f'{START_YEAR}-01-01'), pd.Timestamp(f'{END_YEAR}-12-31'))

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
    p_val, r_eff, p_str = compare_synchrony_arcs(pearson_sync, 'Synchrony_Mean')
    sig = format_significance(p_val)

    ax_pearson.text(0.95, 0.95, f'{p_str}\n{sig}\nr = {r_eff:.2f}' if not np.isnan(r_eff) else f'{p_str}\n{sig}',
                   transform=ax_pearson.transAxes, ha='right', va='top', fontsize=10,
                   bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax_pearson.set_ylabel('Mean Synchrony (Fisher-z)', fontsize=12, fontweight='bold')
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
    p_val, r_eff, p_str = compare_synchrony_arcs(spectral_sync, 'Coherence_30d')
    sig = format_significance(p_val)

    ax_spec.text(0.95, 0.95, f'{p_str}\n{sig}\nr = {r_eff:.2f}' if not np.isnan(r_eff) else f'{p_str}\n{sig}',
                  transform=ax_spec.transAxes, ha='right', va='top', fontsize=10,
                  bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    ax_spec.set_ylabel('Spectral Coherence (30d)', fontsize=12, fontweight='bold')
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
    print(f"  ✓ Saved: {os.path.basename(output_path)}")
```

---

## 6. Fluxo Principal de Execução

```python
def main():
    print("=" * 60)
    print("TEMPORAL SYNCHRONY ANALYSIS: Inner vs Outer Arc")
    print(f"Period: {START_YEAR} - {END_YEAR}")
    print("=" * 60)

    # ====================
    # STEP 1: Load Sites
    # ====================
    sites_df = load_sites_data()

    # ====================
    # STEP 2: Load Satellite Data
    # ====================
    print("\n" + "=" * 60)
    print("STEP 1: Loading Satellite Data")
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
    print("\n  Aggregating time series by arc...")
    sst_arc = aggregate_timeseries_by_arc(sst_ts, sites_df)
    dli_arc = aggregate_timeseries_by_arc(dli_ts, sites_df)
    chl_arc = aggregate_timeseries_by_arc(chl_ts, sites_df)

    print(f"    Inner: SST={len(sst_arc.get('inner', []))} days, DLI={len(dli_arc.get('inner', []))} days")
    print(f"    Outer: SST={len(sst_arc.get('outer', []))} days, DLI={len(dli_arc.get('outer', []))} days")

    # ====================
    # STEP 3: Calculate Synchrony Metrics
    # ====================
    print("\n" + "=" * 60)
    print("STEP 2: Calculating Synchrony Metrics")
    print("=" * 60)

    print("\n  Calculating Pearson synchrony...")
    pearson_sync = calculate_pearson_synchrony(sst_ts, dli_ts, chl_ts, sites_df)
    print(f"    ✓ Inner: n={(pearson_sync['Arc'] == 'inner').sum()}, Outer: n={(pearson_sync['Arc'] == 'outer').sum()}")

    print("\n  Calculating spectral coherence...")
    spectral_sync = calculate_spectral_synchrony(sst_ts, dli_ts, chl_ts, sites_df)
    print(f"    ✓ Inner: n={(spectral_sync['Arc'] == 'inner').sum()}, Outer: n={(spectral_sync['Arc'] == 'outer').sum()}")

    # ====================
    # STEP 4: Statistical Comparison
    # ====================
    print("\n" + "=" * 60)
    print("STEP 3: Statistical Comparison")
    print("=" * 60)

    print("\n  Pearson Synchrony (Inner > Outer):")
    p_val, r_eff, p_str = compare_synchrony_arcs(pearson_sync, 'Synchrony_Mean')
    print(f"    {p_str}, effect size r = {r_eff:.3f}")

    print("\n  Spectral Coherence 30d (Inner > Outer):")
    p_val, r_eff, p_str = compare_synchrony_arcs(spectral_sync, 'Coherence_30d')
    print(f"    {p_str}, effect size r = {r_eff:.3f}")

    # ====================
    # STEP 5: Generate Figures
    # ====================
    print("\n" + "=" * 60)
    print("STEP 4: Generating Publication Figure")
    print("=" * 60)

    create_publication_figure(
        sst_arc, dli_arc, chl_arc,
        pearson_sync, spectral_sync,
        os.path.join(OUTPUT_DIR, 'Fig_Main_Temporal_Synchrony.png')
    )

    # ====================
    # STEP 6: Save Results
    # ====================
    print("\n" + "=" * 60)
    print("STEP 5: Saving Results")
    print("=" * 60)

    pearson_sync.to_csv(os.path.join(OUTPUT_DIR, 'Pearson_Synchrony_by_Site.csv'), index=False)
    print("  ✓ Saved: Pearson_Synchrony_by_Site.csv")

    spectral_sync.to_csv(os.path.join(OUTPUT_DIR, 'Spectral_Coherence_by_Site.csv'), index=False)
    print("  ✓ Saved: Spectral_Coherence_by_Site.csv")

    # Summary statistics
    summary = {
        'Metric': ['Pearson_Synchrony_Inner', 'Pearson_Synchrony_Outer',
                  'Coherence_7d_Inner', 'Coherence_7d_Outer',
                  'Coherence_30d_Inner', 'Coherence_30d_Outer'],
        'Mean': [
            pearson_sync[pearson_sync['Arc'] == 'inner']['Synchrony_Mean'].mean(),
            pearson_sync[pearson_sync['Arc'] == 'outer']['Synchrony_Mean'].mean(),
            spectral_sync[spectral_sync['Arc'] == 'inner']['Coherence_7d'].mean(),
            spectral_sync[spectral_sync['Arc'] == 'outer']['Coherence_7d'].mean(),
            spectral_sync[spectral_sync['Arc'] == 'inner']['Coherence_30d'].mean(),
            spectral_sync[spectral_sync['Arc'] == 'outer']['Coherence_30d'].mean()
        ],
        'SD': [
            pearson_sync[pearson_sync['Arc'] == 'inner']['Synchrony_Mean'].std(),
            pearson_sync[pearson_sync['Arc'] == 'outer']['Synchrony_Mean'].std(),
            spectral_sync[spectral_sync['Arc'] == 'inner']['Coherence_7d'].std(),
            spectral_sync[spectral_sync['Arc'] == 'outer']['Coherence_7d'].std(),
            spectral_sync[spectral_sync['Arc'] == 'inner']['Coherence_30d'].std(),
            spectral_sync[spectral_sync['Arc'] == 'outer']['Coherence_30d'].std()
        ],
        'N': [
            (pearson_sync['Arc'] == 'inner').sum(),
            (pearson_sync['Arc'] == 'outer').sum(),
            (spectral_sync['Arc'] == 'inner').sum(),
            (spectral_sync['Arc'] == 'outer').sum(),
            (spectral_sync['Arc'] == 'inner').sum(),
            (spectral_sync['Arc'] == 'outer').sum()
        ]
    }

    pd.DataFrame(summary).to_csv(os.path.join(OUTPUT_DIR, 'Summary_Statistics.csv'), index=False)
    print("  ✓ Saved: Summary_Statistics.csv")

    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE!")
    print(f"All outputs saved to: {OUTPUT_DIR}")
    print("=" * 60)

if __name__ == "__main__":
    main()
```

---

## 7. Checklist de Implementação

### Fase 1: Configuração e Setup
- [ ] Copiar estrutura de imports e logging do GEMINI.py
- [ ] Definir constantes de cores (COLOR_INNER, COLOR_OUTER, COLORS_ARC)
- [ ] Configurar plt.rcParams para estilo Nature/Science
- [ ] Setar período 2002-2008 e diretórios corretos

### Fase 2: Carregamento de Dados
- [ ] Implementar load_sites_data() (padrão GEMINI.py)
- [ ] Implementar load_timeseries_for_sites() com progress reporting
- [ ] Implementar calculate_dli_for_sites()
- [ ] Implementar aggregate_timeseries_by_arc()
- [ ] Testar carregamento com dados reais

### Fase 3: Análise de Sincronicidade
- [ ] Implementar calculate_pearson_synchrony()
  - [ ] Padronização (z-score)
  - [ ] Correlações pareadas
  - [ ] Transformação de Fisher para média
- [ ] Implementar compute_spectral_coherence()
  - [ ] Lomb-Scargle para ambas as séries
  - [ ] Cálculo de cross-power
  - [ ] Extração em período específico
- [ ] Implementar calculate_spectral_synchrony()
- [ ] Implementar compare_synchrony_arcs() (Mann-Whitney U)

### Fase 4: Visualização
- [ ] Implementar smooth_timeseries() (rolling mean)
- [ ] Implementar create_publication_figure()
  - [ ] Panel A: Inner Arc (raw points + smooth lines, 3 eixos Y)
  - [ ] Panel B: Outer Arc (idêntico ao Inner)
  - [ ] Panel C: Pearson Synchrony boxplot
  - [ ] Panel D: Spectral Coherence boxplot
  - [ ] Labels de painel (A, B, C, D)
  - [ ] Legendas e títulos profissionais

### Fase 5: Execução e Validação
- [ ] Implementar main() com fluxo completo
- [ ] Testar com subset de dados (1 ano) para validação
- [ ] Executar análise completa 2002-2008
- [ ] Verificar outputs (CSVs + figura)

### Fase 6: Refinamento
- [ ] Ajustar parâmetros de suavização (window_days)
- [ ] Otimizar cores e tamanhos de fonte
- [ ] Adicionar tratamento de erros robusto
- [ ] Documentar código com docstrings

---

## 8. Dependências Python

```python
# Core
pandas, numpy, xarray

# Analysis
scipy (stats, signal, lombscargle, mannwhitneyu)

# Smoothing (optional, can use rolling)
statsmodels (for lowess alternative)

# Visualization
matplotlib

# File handling
glob, os, logging
```

---

## 9. Saídas Esperadas

### Arquivos CSV
1. `Pearson_Synchrony_by_Site.csv` - Métricas de sincronicidade por site
2. `Spectral_Coherence_by_Site.csv` - Coerência espectral por site
3. `Summary_Statistics.csv` - Resumo estatístico Inner vs Outer

### Figuras
1. `Fig_Main_Temporal_Synchrony.png` - Figura principal 2x2 publicável

### Log
1. `temporal_synchrony_analysis_errors.log` - Erros durante execução

---

## 10. Notas Importantes

### Coerência com GEMINI.py
- Usar mesmas cores: COLOR_INNER (#CD853F), COLOR_OUTER (#4169E1)
- Mesma fonte: Arial
- Mesmos tamanhos de fonte definidos em plt.rcParams
- Mesmo padrão de boxplots com pontos individuais sobrepostos

### Correlação de Pearson
- Usar transformação de Fisher (arctanh) para calcular média de correlações
- Padronizar séries (z-score) antes de calcular correlação
- Reportar tanto médias quanto correlações pareadas (SST-DLI, SST-CHL, DLI-CHL)

### Coerência Espectral
- Focar em períodos de 7 dias e 30 dias
- Usar Lomb-Scargle (já implementado no GEMINI.py)
- Coerência = raiz quadrada do produto das potências normalizadas

### Teste Estatístico
- Mann-Whitney U com alternative='greater' (testa se Inner > Outer)
- Calcular tamanho de efeito (r = z / sqrt(n1 + n2))
- Formatar p-value com asteriscos de significância

### Figura Nature/Science
- Resolução 300 DPI
- Fontes Arial, tamanhos hierárquicos
- Sem spines top/right
- Legendas limpas, framealpha=0.9
- Labels de painel (A, B, C, D) em negrito

---

## 11. Possíveis Extensões Futuras

1. **Análise de Wavelet Coherence**: Para visualizar sincronicidade no tempo-frequência
2. **Granger Causality**: Para testar causalidade direcional entre variáveis
3. **Mutual Information**: Para detectar dependências não-lineares
4. **Phase Synchronization**: Para avaliar sincronização de fase entre oscilações

---

**Fim do Plano de Implementação**
