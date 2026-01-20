### TIME_SERIES_UNIFIED_SYNCHRONY_ANALYSIS.py ###
# Unified analysis of temporal synchrony between SST, DLI, and Chl-a
# Part of the Menezes et al. Abrolhos coral research.
#
# This script performs:
# 1. Pearson Synchrony (by site)
# 2. Cross-Correlation Lag Analysis (by arc)
# 3. STL Decomposition for Pulse Synchrony (residents correlation)
# 4. Spectral Coherence (optional/complementary)
#
# All results saved to Variability_Analysis_GEMINI as requested.

import os
import glob
import re
import logging
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.dates as mdates
from matplotlib.lines import Line2D
from scipy import stats
from scipy.signal import lombscargle
from statsmodels.tsa.seasonal import STL
from datetime import datetime
import warnings

warnings.filterwarnings('ignore')

# ============================================
# Configuration
# ============================================

# --- Time Period ---
START_YEAR = 2002
END_YEAR = 2008

# --- Data Directories ---
SST_DIR = r'K:\remote sensing\CRW_SST_FULL'
MODIS_DIR = r'K:\remote sensing\MODIS_DATA_FULL'

# --- Sites CSV ---
SITES_CSV_PATH = r'C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######22.04.23\DATA\sites_list_full_clean.csv'

# --- Output Directory (Unified with Variability Analysis) ---
OUTPUT_DIR = r'C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_RESULTS\Variability_Analysis_GEMINI'
os.makedirs(OUTPUT_DIR, exist_ok=True)

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
COLOR_SST = '#D55E00'    # Vermelho-alaranjado
COLOR_DLI = '#0072B2'    # Azul escuro
COLOR_CHL = '#009E73'    # Verde-teal

# ============================================
# Data Loading Functions (Replicated from Variability Script)
# ============================================

def parse_date_from_filename(filename):
    basename = os.path.basename(filename)
    match = re.search(r'(\d{8})', basename)
    if match:
        try:
            return pd.to_datetime(match.group(1), format='%Y%m%d')
        except ValueError:
            return None
    return None

def filter_files_by_period(files, start_year, end_year):
    filtered = []
    for f in files:
        date = parse_date_from_filename(f)
        if date is not None and start_year <= date.year <= end_year:
            filtered.append(f)
    return sorted(filtered)

def load_timeseries_for_sites(files, var_names, sites_df, var_label):
    print(f"\n  Loading {var_label}...")
    site_data = {row['Site_name']: {'times': [], 'values': []} for _, row in sites_df.iterrows()}
    
    for i, f in enumerate(files):
        if (i + 1) % 1000 == 0:
            print(f"    Processing file {i+1}/{len(files)}...")
        date = parse_date_from_filename(f)
        if date is None: continue
        try:
            with xr.open_dataset(f) as ds:
                var_name = next((v for v in var_names if v in ds.data_vars), None)
                if var_name is None: continue
                data = ds[var_name]
                if 'time' in data.dims: data = data.isel(time=0)
                for _, row in sites_df.iterrows():
                    site_name = row['Site_name']
                    try:
                        value = float(data.sel(lat=row['Latitude'], lon=row['Longitude'], method='nearest').values)
                        if not np.isnan(value):
                            site_data[site_name]['times'].append(date)
                            site_data[site_name]['values'].append(value)
                    except: continue
        except: continue
    
    result = {}
    for site_name, data in site_data.items():
        if len(data['times']) > 0:
            result[site_name] = pd.Series(data['values'], index=pd.DatetimeIndex(data['times'])).sort_index()
    return result

def calculate_dli_for_sites(kd490_files, par_files, sites_df):
    print("\n  Calculating DLI...")
    par_by_date = {parse_date_from_filename(f): f for f in par_files if parse_date_from_filename(f)}
    site_data = {row['Site_name']: {'times': [], 'values': []} for _, row in sites_df.iterrows()}
    
    for i, f_kd in enumerate(kd490_files):
        if (i + 1) % 1000 == 0:
            print(f"    Processing file {i+1}/{len(kd490_files)}...")
        date = parse_date_from_filename(f_kd)
        if date is None or date not in par_by_date: continue
        try:
            with xr.open_dataset(f_kd) as ds_kd, xr.open_dataset(par_by_date[date]) as ds_par:
                kd490_data, par_data = ds_kd['Kd_490'], ds_par['par']
                if 'time' in kd490_data.dims: kd490_data = kd490_data.isel(time=0)
                if 'time' in par_data.dims: par_data = par_data.isel(time=0)
                for _, row in sites_df.iterrows():
                    try:
                        kd490_val = float(kd490_data.sel(lat=row['Latitude'], lon=row['Longitude'], method='nearest').values)
                        par_val = float(par_data.sel(lat=row['Latitude'], lon=row['Longitude'], method='nearest').values)
                        if np.isnan(kd490_val) or np.isnan(par_val) or kd490_val <= 0: continue
                        kdpar = KDPAR_A + KDPAR_B * kd490_val - KDPAR_C / (kd490_val + 1e-9)
                        dli = par_val * np.exp(-kdpar * row['Depth_m'])
                        if not np.isnan(dli) and dli > 0:
                            site_data[row['Site_name']]['times'].append(date)
                            site_data[row['Site_name']]['values'].append(dli)
                    except: continue
        except: continue
    
    result = {}
    for site_name, data in site_data.items():
        if len(data['times']) > 0:
            result[site_name] = pd.Series(data['values'], index=pd.DatetimeIndex(data['times'])).sort_index()
    return result

# ============================================
# Synchrony Analysis Functions
# ============================================

def calculate_pearson_synchrony(sst_ts, dli_ts, chl_ts, sites_df):
    results = []
    for _, row in sites_df.iterrows():
        site = row['Site_name']
        if site in sst_ts and site in dli_ts and site in chl_ts:
            df = pd.DataFrame({'SST': sst_ts[site], 'DLI': dli_ts[site], 'CHL': chl_ts[site]}).dropna()
            if len(df) > 30:
                corr_matrix = df.corr()
                results.append({
                    'Site_name': site,
                    'Arc': row['Arc'],
                    'Corr_SST_DLI': corr_matrix.loc['SST', 'DLI'],
                    'Corr_SST_CHL': corr_matrix.loc['SST', 'CHL'],
                    'Corr_DLI_CHL': corr_matrix.loc['DLI', 'CHL']
                })
    return pd.DataFrame(results)

def calculate_cross_correlation_full(series1, series2, max_lag=60):
    df = pd.DataFrame({'s1': series1, 's2': series2}).dropna()
    if len(df) < 60: return None
    lags = range(-max_lag, max_lag + 1)
    correlations = [df['s2'].corr(df['s1'].shift(lag)) for lag in lags]
    correlations = np.array(correlations)
    optimal_idx = np.nanargmax(np.abs(correlations))
    return {
        'lags': np.array(list(lags)),
        'correlations': correlations,
        'optimal_lag': list(lags)[optimal_idx],
        'optimal_corr': correlations[optimal_idx]
    }

def calculate_pulse_synchrony_stl(sst, dli, chl, period=365):
    sst_interp = sst.interpolate(method='time', limit=7)
    dli_interp = dli.interpolate(method='time', limit=7)
    chl_interp = chl.interpolate(method='time', limit=7)
    df = pd.DataFrame({'SST': sst_interp, 'DLI': dli_interp, 'CHL': chl_interp}).dropna()
    if len(df) < period * 2: return None
    
    try:
        resid_sst = STL(df['SST'], period=period, robust=True).fit().resid
        resid_dli = STL(df['DLI'], period=period, robust=True).fit().resid
        resid_chl = STL(df['CHL'], period=period, robust=True).fit().resid
        return {
            'Pulse_Corr_SST_DLI': resid_sst.corr(resid_dli),
            'Pulse_Corr_SST_CHL': resid_sst.corr(resid_chl),
            'Pulse_Corr_DLI_CHL': resid_dli.corr(resid_chl),
            'Resid_SST': resid_sst, 'Resid_DLI': resid_dli, 'Resid_CHL': resid_chl
        }
    except: return None

def mann_whitney_test(inner_vals, outer_vals):
    inner = [v for v in inner_vals if not np.isnan(v)]
    outer = [v for v in outer_vals if not np.isnan(v)]
    if len(inner) < 3 or len(outer) < 3: return np.nan, "N/A"
    stat, p = stats.mannwhitneyu(inner, outer, alternative='two-sided')
    p_str = "p < 0.001" if p < 0.001 else f"p = {p:.3f}" if p < 0.05 else f"p = {p:.2f}"
    return p, p_str

# ============================================
# Raw Time Series Figure Function
# ============================================

def generate_raw_timeseries_figure(arc_data: dict, output_dir: str, start_year: int, end_year: int) -> None:
    """
    Generates a supplementary figure showing raw (unsmoothed) time series
    for SST, DLI, and Chl-a across Inner and Outer arcs.
    
    Uses the same 3-axis layout pattern as TIME_SERIES_LAG_temporal_environ.py
    with scatter points instead of lines to display raw data granularity.
    
    Args:
        arc_data: Dictionary with keys 'inner' and 'outer', each containing
                  sub-dictionary with 'sst', 'dli', 'chl' as pd.Series.
        output_dir: Path to save the output figure.
        start_year: Start year for title.
        end_year: End year for title.
    
    Output:
        Saves 'Fig_Raw_Timeseries_SST_DLI_CHL.png' at 300 DPI.
    """
    fig, axes = plt.subplots(
        nrows=2, ncols=1, 
        figsize=(20, 12), 
        sharex=True
    )
    
    arc_order = ['inner', 'outer']
    arc_labels = ['Inner Arc', 'Outer Arc']
    
    for i, (arc, arc_label) in enumerate(zip(arc_order, arc_labels)):
        ax = axes[i]
        
        # Extract raw data (no smoothing)
        sst_raw = arc_data[arc]['sst']
        dli_raw = arc_data[arc]['dli']
        chl_raw = arc_data[arc]['chl']
        
        # --- AXIS 1 (Left): SST ---
        ax.set_ylabel('SST (°C)', color=COLOR_SST, fontsize=14, fontweight='bold')
        p1, = ax.plot(
            sst_raw.index, sst_raw.values, 
            'o', markersize=3, alpha=0.6, 
            color=COLOR_SST, label='SST', 
            linestyle='None'
        )
        ax.tick_params(axis='y', labelcolor=COLOR_SST, labelsize=12)
        ax.yaxis.grid(True, linestyle='--', alpha=0.4, color=COLOR_SST)
        
        # --- AXIS 2 (Right 1): DLI ---
        ax2 = ax.twinx()
        ax2.set_ylabel('DLI (mol m⁻² d⁻¹)', color=COLOR_DLI, fontsize=14, fontweight='bold')
        p2, = ax2.plot(
            dli_raw.index, dli_raw.values, 
            'o', markersize=3, alpha=0.6, 
            color=COLOR_DLI, label='DLI', 
            linestyle='None'
        )
        ax2.tick_params(axis='y', labelcolor=COLOR_DLI, labelsize=12)
        
        # --- AXIS 3 (Right 2, offset): Chl-a ---
        ax3 = ax.twinx()
        ax3.spines['right'].set_position(('outward', 70))
        ax3.set_ylabel('Chl-a (mg m⁻³)', color=COLOR_CHL, fontsize=14, fontweight='bold')
        p3, = ax3.plot(
            chl_raw.index, chl_raw.values, 
            'o', markersize=3, alpha=0.6, 
            color=COLOR_CHL, label='Chl-a', 
            linestyle='None'
        )
        ax3.tick_params(axis='y', labelcolor=COLOR_CHL, labelsize=12)
        
        # --- Legend ---
        lines = [p1, p2, p3]
        ax.legend(lines, [l.get_label() for l in lines], loc='upper left', fontsize=11)
        
        # --- Title with Panel Letter ---
        ax.set_title(
            f'Panel {chr(65+i)}: {arc_label} – Raw Daily Data', 
            fontweight='bold', loc='left', fontsize=12
        )
        
        # --- X-axis formatting (only on bottom panel) ---
        if i == len(arc_order) - 1:
            ax.xaxis.set_major_locator(mdates.YearLocator(1))
            ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=3))
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
            plt.setp(ax.get_xticklabels(), rotation=0, ha='center', fontsize=12)
    
    # --- Super Title ---
    fig.suptitle(
        f'Raw Daily Time Series: SST, DLI, and Chl-a ({start_year}–{end_year})', 
        fontsize=16, fontweight='bold', y=1.02
    )
    
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    
    output_path = os.path.join(output_dir, 'Fig_Raw_Timeseries_SST_DLI_CHL.png')
    fig.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    print(f"\n  [NEW FIGURE] Raw time series saved: {output_path}")

# ============================================
# Main Execution
# ============================================

if __name__ == "__main__":
    print("\n" + "=" * 60)
    print("UNIFIED SYNCHRONY ANALYSIS: SST, DLI, Chl-a")
    print("=" * 60)
    
    sites_df = pd.read_csv(SITES_CSV_PATH, sep=';')
    sites_df['Arc'] = sites_df['Arc'].str.strip().str.lower()
    sites_df = sites_df[sites_df['Arc'].isin(['inner', 'outer'])].copy()
    
    sst_files = filter_files_by_period(glob.glob(os.path.join(SST_DIR, SST_PATTERN)), START_YEAR, END_YEAR)
    kd490_files = filter_files_by_period(glob.glob(os.path.join(MODIS_DIR, KD490_PATTERN)), START_YEAR, END_YEAR)
    par_files = filter_files_by_period(glob.glob(os.path.join(MODIS_DIR, PAR_PATTERN)), START_YEAR, END_YEAR)
    chl_files = filter_files_by_period(glob.glob(os.path.join(MODIS_DIR, CHL_PATTERN)), START_YEAR, END_YEAR)
    
    sst_ts = load_timeseries_for_sites(sst_files, ['analysed_sst', 'sst'], sites_df, 'SST')
    chl_ts = load_timeseries_for_sites(chl_files, ['chlor_a'], sites_df, 'Chl-a')
    dli_ts = calculate_dli_for_sites(kd490_files, par_files, sites_df)
    
    # 1. Pearson Synchrony
    pearson_df = calculate_pearson_synchrony(sst_ts, dli_ts, chl_ts, sites_df)
    pearson_df.to_csv(os.path.join(OUTPUT_DIR, 'Pearson_Synchrony_by_Site.csv'), index=False)
    
    # 2. Aggregation and CCF
    arc_results = []
    stl_results = []
    arc_data = {}  # Store aggregated data for each arc
    
    for arc in ['inner', 'outer']:
        arc_sites = sites_df[sites_df['Arc'] == arc]['Site_name'].tolist()
        
        # Aggregate (Mean Daily)
        sst_arc = pd.concat([sst_ts[s] for s in arc_sites if s in sst_ts], axis=1).mean(axis=1)
        dli_arc = pd.concat([dli_ts[s] for s in arc_sites if s in dli_ts], axis=1).mean(axis=1)
        chl_arc = pd.concat([chl_ts[s] for s in arc_sites if s in chl_ts], axis=1).mean(axis=1)
        
        arc_data[arc] = {'sst': sst_arc, 'dli': dli_arc, 'chl': chl_arc}
        
        # CCF Analysis (all 3 pairs)
        ccf_sst_dli = calculate_cross_correlation_full(sst_arc, dli_arc)
        if ccf_sst_dli:
            arc_results.append({**ccf_sst_dli, 'Arc': arc, 'Pair': 'SST_DLI'})
        
        ccf_sst_chl = calculate_cross_correlation_full(sst_arc, chl_arc)
        if ccf_sst_chl:
            arc_results.append({**ccf_sst_chl, 'Arc': arc, 'Pair': 'SST_CHL'})
        
        ccf_dli_chl = calculate_cross_correlation_full(dli_arc, chl_arc)
        if ccf_dli_chl:
            arc_results.append({**ccf_dli_chl, 'Arc': arc, 'Pair': 'DLI_CHL'})
        
        # STL Analysis
        stl = calculate_pulse_synchrony_stl(sst_arc, dli_arc, chl_arc)
        if stl:
            stl_results.append({
                'Arc': arc,
                'Pulse_Corr_SST_DLI': stl['Pulse_Corr_SST_DLI'],
                'Pulse_Corr_SST_CHL': stl['Pulse_Corr_SST_CHL'],
                'Pulse_Corr_DLI_CHL': stl['Pulse_Corr_DLI_CHL']
            })

    # ============================================
    # Figure: Complete 5x3 Layout
    # ============================================
    fig_main = plt.figure(figsize=(18, 24))
    gs = fig_main.add_gridspec(5, 3, height_ratios=[1.2, 1.2, 1, 1, 0.8], hspace=0.35, wspace=0.25)
    
    # Row 1-2: Time Series (A & B) - Full width
    for i, arc in enumerate(['inner', 'outer']):
        ax = fig_main.add_subplot(gs[i, :])
        sst_smooth = arc_data[arc]['sst'].rolling(14, center=True).mean()
        dli_smooth = arc_data[arc]['dli'].rolling(14, center=True).mean()
        chl_smooth = arc_data[arc]['chl'].rolling(14, center=True).mean()
        
        ax.plot(sst_smooth.index, sst_smooth.values, color=COLOR_SST, lw=2, label='SST')
        ax.set_ylabel('SST (°C)', color=COLOR_SST, fontweight='bold')
        ax.tick_params(axis='y', labelcolor=COLOR_SST)
        
        ax2 = ax.twinx()
        ax2.plot(dli_smooth.index, dli_smooth.values, color=COLOR_DLI, lw=1.5, alpha=0.8, label='DLI')
        ax2.set_ylabel('DLI (mol m⁻² d⁻¹)', color=COLOR_DLI, fontweight='bold')
        ax2.tick_params(axis='y', labelcolor=COLOR_DLI)
        
        ax3 = ax.twinx()
        ax3.spines['right'].set_position(('outward', 60))
        ax3.plot(chl_smooth.index, chl_smooth.values, color=COLOR_CHL, lw=1.5, alpha=0.8, label='Chl-a')
        ax3.set_ylabel('Chl-a (mg m⁻³)', color=COLOR_CHL, fontweight='bold')
        ax3.tick_params(axis='y', labelcolor=COLOR_CHL)
        
        ax.set_title(f'Panel {chr(65+i)}: {arc.capitalize()} Arc Combined Dynamics (14-day smooth)', fontweight='bold', loc='left', fontsize=11)
    
    # Row 3: Pearson Boxplots (C, D, E)
    pearson_pairs = [('Corr_SST_DLI', 'SST↔DLI'), ('Corr_SST_CHL', 'SST↔Chl'), ('Corr_DLI_CHL', 'DLI↔Chl')]
    for j, (col, label) in enumerate(pearson_pairs):
        ax = fig_main.add_subplot(gs[2, j])
        inner_vals = pearson_df[pearson_df['Arc'] == 'inner'][col].dropna()
        outer_vals = pearson_df[pearson_df['Arc'] == 'outer'][col].dropna()
        
        bp = ax.boxplot([inner_vals, outer_vals], labels=['Inner', 'Outer'], patch_artist=True, widths=0.75)
        bp['boxes'][0].set_facecolor(COLOR_INNER)
        bp['boxes'][1].set_facecolor(COLOR_OUTER)
        for box in bp['boxes']:
            box.set_alpha(0.7)
        
        # Mann-Whitney test
        p_val, p_str = mann_whitney_test(inner_vals, outer_vals)
        ax.text(0.88, 0.82, p_str, transform=ax.transAxes, ha='right', va='top', fontsize=9, 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        ax.axhline(0, color='gray', alpha=0.3, ls='--')
        ax.set_title(f'Panel {chr(67+j)}: {label} (Pearson)', fontweight='bold', loc='left', fontsize=10)
        ax.set_ylabel('Correlation')
    
    # Row 4: Cross-Correlation (F, G, H)
    ccf_pairs = [('SST_DLI', 'SST → DLI'), ('SST_CHL', 'SST → Chl'), ('DLI_CHL', 'DLI → Chl')]
    for j, (pair_key, label) in enumerate(ccf_pairs):
        ax = fig_main.add_subplot(gs[3, j])
        for res in [r for r in arc_results if r['Pair'] == pair_key]:
            color = COLOR_INNER if res['Arc'] == 'inner' else COLOR_OUTER
            ax.plot(res['lags'], res['correlations'], color=color, lw=2, 
                   label=f"{res['Arc'].capitalize()} (Lag: {res['optimal_lag']}d)")
        ax.axvline(0, color='black', alpha=0.3, ls='--')
        ax.axhline(0, color='gray', alpha=0.3)
        ax.legend(fontsize=8, loc='best')
        ax.set_title(f'Panel {chr(70+j)}: CCF ({label})', fontweight='bold', loc='left', fontsize=10)
        ax.set_xlabel('Lag (days)')
        ax.set_ylabel('Correlation')
    
    # Row 5: STL Pulse Bar Chart (I) - Full width
    ax_i = fig_main.add_subplot(gs[4, :])
    if stl_results:
        stl_df = pd.DataFrame(stl_results)
        x = np.arange(len(stl_df))
        width = 0.25
        
        bars1 = ax_i.bar(x - width, stl_df['Pulse_Corr_SST_DLI'], width, label='SST↔DLI', color=COLOR_SST, alpha=0.85)
        bars2 = ax_i.bar(x, stl_df['Pulse_Corr_SST_CHL'], width, label='SST↔Chl', color='#888888', alpha=0.85)
        bars3 = ax_i.bar(x + width, stl_df['Pulse_Corr_DLI_CHL'], width, label='DLI↔Chl', color=COLOR_CHL, alpha=0.85)
        
        # Add value labels on bars
        for bars in [bars1, bars2, bars3]:
            for bar in bars:
                height = bar.get_height()
                ax_i.annotate(f'{height:.2f}', xy=(bar.get_x() + bar.get_width()/2, height),
                             xytext=(0, 3 if height >= 0 else -10), textcoords='offset points',
                             ha='center', va='bottom' if height >= 0 else 'top', fontsize=8)
        
        ax_i.set_xticks(x)
        ax_i.set_xticklabels(stl_df['Arc'].str.capitalize(), fontsize=11)
        ax_i.legend(fontsize=10, loc='upper right')
        ax_i.axhline(0, color='black', alpha=0.5)
    ax_i.set_title('Panel I: STL Residual (Pulse) Correlations', fontweight='bold', loc='left', fontsize=11)
    ax_i.set_ylabel('Residual Correlation', fontsize=10)
    ax_i.set_ylim(-1, 1)
    
    plt.tight_layout()
    fig_main.savefig(os.path.join(OUTPUT_DIR, 'Fig_Main_Temporal_Synchrony.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    pd.DataFrame(arc_results).to_csv(os.path.join(OUTPUT_DIR, 'Cross_Correlation_Results.csv'), index=False)
    pd.DataFrame(stl_results).to_csv(os.path.join(OUTPUT_DIR, 'STL_Pulse_Synchrony.csv'), index=False)
    
    # --- Generate Raw Time Series Figure ---
    generate_raw_timeseries_figure(arc_data, OUTPUT_DIR, START_YEAR, END_YEAR)
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE! Outputs in Variability_Analysis_GEMINI")
    print("=" * 60)
