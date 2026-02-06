##### TIME_SERIES e LAG temporal entre SST vs DLI vs CHLO 
### AGRUPAMENTO POR ARCO

import os
import xarray as xr
import numpy as np
import logging
import glob
import re
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime
import pandas as pd

# ============================================
# Configurações Iniciais
# ============================================

logging.basicConfig(
    filename='file_open_errors_timeseries.log',
    filemode='w',
    level=logging.ERROR,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# --- Chaves para seleção de período ---
start_year = 2019
end_year = 2020

# --- Filtro de data de início específica ---
use_specific_start_date = True
specific_start_date = '2019-01-01'

# --- Diretórios dos dados ---
sst_dir = r'E:\remote sensing\CRW_SST_FULL'
modis_dir = r'E:\remote sensing\MODIS_DATA_FULL'
bathymetry_file = r'C:\Users\rbfra\OneDrive\GIS shapes\GEBCO_2024b.nc'
output_dir = r'C:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\########NEW RESULTS\ALCA_teste'

# --- Caminho para o arquivo de sites ---
sites_csv_path = r'C:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######22.04.23\DATA\sites_list_full_clean_PLEST.csv'

os.makedirs(output_dir, exist_ok=True)

# --- Limites da região de interesse ---
lat_min, lat_max = -18.3, -17.3
lon_min, lon_max = -39.5, -38.3

# --- Padrões de busca para os arquivos ---
sst_pattern_crw = 'coraltemp_v3.1_*.nc'
kd490_pattern = 'AQUA_MODIS.*.L3m.DAY.KD.Kd_490.4km.nc'
par_pattern = 'AQUA_MODIS.*.L3m.DAY.PAR.par.4km.nc'
### NOVO ### Adicionado padrão para clorofila
chl_pattern = 'AQUA_MODIS.*.L3m.DAY.CHL.chlor_a.4km.nc'


# --- Carregar e PROCESSAR a lista de sítios ---
try:
    sites_df = pd.read_csv(sites_csv_path, sep=';')
    sites_df.columns = [col.strip().upper() for col in sites_df.columns]
    
    # Mapeamento dos nomes de colunas
    column_mapping = {
        'LATITUDE': 'LAT',
        'LONGITUDE': 'LONG',
        'DEPTH_M': 'DEPTH_M',
        'SITE_NAME': 'SITE_NAME',
        'ARCH': 'ARCH'
    }
    
    sites_df = sites_df.rename(columns=column_mapping)
    
    required_cols = {'SITE_NAME', 'LAT', 'LONG', 'ARCH', 'DEPTH_M'}
    if not required_cols.issubset(sites_df.columns):
        missing_cols = required_cols - set(sites_df.columns)
        raise ValueError(f"Colunas necessárias não encontradas: {missing_cols}. Colunas disponíveis: {list(sites_df.columns)}")

    def classify_arch(arch_string):
        if 'inner' in arch_string.lower():
            return 'Inner Arch'
        elif 'outer' in arch_string.lower():
            return 'Outer Arch'
        return 'Other'

    sites_df['ARCH_SIMPLE'] = sites_df['ARCH'].apply(classify_arch)
    sites_df = sites_df[sites_df['ARCH_SIMPLE'].isin(['Inner Arch', 'Outer Arch'])].copy()
    
    print(f"{len(sites_df)} sites carregados com sucesso e classificados como 'Inner' ou 'Outer'.")
    print("Contagem de sites por 'ARCH_SIMPLE':")
    print(sites_df['ARCH_SIMPLE'].value_counts())

except FileNotFoundError:
    print(f"ERRO CRÍTICO: Arquivo de sites não encontrado em: {sites_csv_path}")
    exit()
except Exception as e:
    print(f"ERRO CRÍTICO ao ler o arquivo de sites: {e}")
    exit()

# ============================================
# Funções Auxiliares
# ============================================

def parse_date_from_filename(filename):
    basename = os.path.basename(filename)
    match = re.search(r'(\d{8})', basename)
    if match:
        try:
            return np.datetime64(datetime.strptime(match.group(1), '%Y%m%d'), 'D')
        except ValueError:
            return None
    logging.error(f"Data não encontrada ou em formato inválido no arquivo: {basename}")
    return None

def filter_files_by_year_range(files, start_y, end_y):
    filtered_list = []
    print(f"Filtrando {len(files)} arquivos entre os anos {start_y} e {end_y}...")
    for f in files:
        date = parse_date_from_filename(f)
        if date is not None:
            year = pd.to_datetime(date).year
            if start_y <= year <= end_y:
                filtered_list.append(f)
    print(f"-> {len(filtered_list)} arquivos permaneceram após o filtro.")
    return filtered_list

def get_lat_slice(ds, lat_min, lat_max):
    lat_values = ds['lat'].values
    return slice(lat_min, lat_max) if lat_values[0] < lat_values[-1] else slice(lat_max, lat_min)

def get_sst_variable(ds):
    possible_vars = ['analysed_sst', 'sea_surface_temperature', 'sst']
    for var in possible_vars:
        if var in ds.data_vars:
            return var
    logging.error("Nenhuma variável válida de SST encontrada.")
    return None

### NOVO ### Função para carregar o stack de referência
def load_reference_stack(pattern, base_dir, var_names, lat_slice, lon_slice, start_y, end_y):
    """Carrega o stack de referência, fatiando cada arquivo para economizar memória."""
    print(f"\nCarregando stack de referência para: {pattern}...")
    all_files = sorted(glob.glob(os.path.join(base_dir, pattern)))
    files_in_period = filter_files_by_year_range(all_files, start_y, end_y)

    if not files_in_period:
        print(f"AVISO: Nenhum arquivo encontrado para {pattern} no período {start_y}-{end_y}.")
        return None
    
    data_list = []
    for f in files_in_period:
        try:
            with xr.open_dataset(f) as ds:
                var_name = next((v for v in var_names if v in ds.data_vars), None)
                if not var_name: continue
                
                date_val = parse_date_from_filename(f)
                if date_val is None: continue

                # A fatia espacial é aplicada AQUI, para cada arquivo, economizando memória
                data_sub = ds[var_name].sel(lat=lat_slice, lon=lon_slice)
                if 'time' in data_sub.dims:
                    data_sub = data_sub.isel(time=0, drop=True)
                
                data_sub = data_sub.expand_dims(dim={'time': [date_val]})
                data_list.append(data_sub)
        except Exception as e:
            logging.error(f"Erro ao abrir ou processar {f}: {e}")
            continue
            
    if not data_list: return None
    return xr.concat(data_list, dim='time').sortby('time')

### NOVO ### Função para carregar e alinhar dados com o stack de referência
def load_and_interpolate_stack(pattern, base_dir, var_names, start_y, end_y, reference_stack):
    """
    Carrega dados de satélite e os interpola imediatamente para a grade de um stack de referência.
    Esta é a abordagem mais segura em termos de memória para alinhar grades.
    """
    print(f"\nCarregando e interpolando dados para: {pattern}...")
    all_files = sorted(glob.glob(os.path.join(base_dir, pattern)))
    files_in_period = filter_files_by_year_range(all_files, start_y, end_y)

    if not files_in_period:
        print(f"AVISO: Nenhum arquivo encontrado para {pattern} no período {start_y}-{end_y}.")
        return None
    
    # Pega a grade de referência (apenas as coordenadas, sem os dados)
    reference_grid = reference_stack.isel(time=0, drop=True)
    
    data_list = []
    for f in files_in_period:
        try:
            with xr.open_dataset(f) as ds:
                var_name = next((v for v in var_names if v in ds.data_vars), None)
                if not var_name: continue
                
                date_val = parse_date_from_filename(f)
                if date_val is None: continue

                data_sub = ds[var_name]
                if 'time' in data_sub.dims:
                    data_sub = data_sub.isel(time=0, drop=True)
                
                # Interpola o dado do dia para a grade de referência ANTES de adicionar à lista
                interp_sub = data_sub.interp_like(reference_grid, method='nearest')
                
                interp_sub = interp_sub.expand_dims(dim={'time': [date_val]})
                data_list.append(interp_sub)
        except Exception as e:
            logging.error(f"Erro ao abrir, processar ou interpolar {f}: {e}")
            continue
            
    if not data_list: 
        print(f"AVISO: Nenhum dado pôde ser extraído e interpolado para {pattern}")
        return None
        
    return xr.concat(data_list, dim='time').sortby('time')

def calculate_dli_benthic_per_site(kd490_files, par_files, sites_df_group):
    par_files_dict = {parse_date_from_filename(f): f for f in par_files if parse_date_from_filename(f)}
    dli_per_site_ts = {}

    for index, site in sites_df_group.iterrows():
        site_name = site['SITE_NAME']
        site_lat, site_lon, site_depth = site['LAT'], site['LONG'], site['DEPTH_M']
        
        dli_values = []
        dates = []
        
        for f_kd in kd490_files:
            date_val = parse_date_from_filename(f_kd)
            if date_val is None or date_val not in par_files_dict:
                continue
            
            f_par = par_files_dict[date_val]
            try:
                with xr.open_dataset(f_kd) as ds_kd490, xr.open_dataset(f_par) as ds_par:
                    kd_point = ds_kd490['Kd_490'].sel(lat=site_lat, lon=site_lon, method='nearest').item(0)
                    par_point = ds_par['par'].sel(lat=site_lat, lon=site_lon, method='nearest').item(0)
                    
                    if np.isnan(kd_point) or np.isnan(par_point):
                        dli_benthic = np.nan
                    else:
                        kdpar = 0.0665 + 0.874 * kd_point - 0.00121 / (kd_point + 1e-9)
                        dli_benthic = par_point * np.exp(-kdpar * site_depth)
                    
                    dli_values.append(dli_benthic)
                    dates.append(date_val)
            except Exception as e:
                logging.error(f"Erro ao processar DLI para {site_name} na data {date_val}: {e}")
                continue
        
        if dli_values:
            dli_per_site_ts[site_name] = xr.DataArray(dli_values, coords={'time': dates}, dims=['time'])

    return dli_per_site_ts

### NOVO ### Função para calcular a correlação cruzada e o lag
def calculate_cross_correlation(series1, series2, max_lag_days=90):
    """
    Calcula a correlação cruzada entre duas séries temporais.
    'series1' é a variável que se supõe que "lidera" (causa).
    'series2' é a variável que se supõe que "responde" (efeito).
    
    Um lag positivo significa que `series1` antecede `series2`.
    Ex: Se DLI é series1 e SST é series2, um lag de +15 dias significa 
    que o pico de SST ocorre 15 dias após o pico de DLI.
    """
    # Alinha as séries e remove os dias onde não há dados para ambos
    df = pd.DataFrame({'s1': series1, 's2': series2}).dropna()
    if len(df) < 30: # Requer um número mínimo de pontos para uma correlação significativa
        return np.nan, np.nan

    correlations = []
    # O lag é o quanto deslocamos a 'series1' para "frente" no tempo
    lags = range(-max_lag_days, max_lag_days + 1)
    
    for lag in lags:
        # df['s1'].shift(lag) desloca s1. Se lag=15, um valor de s1 no dia 1 é comparado com um valor de s2 no dia 16.
        # Isso testa se s1 de X dias atrás se correlaciona com s2 de hoje.
        corr = df['s2'].corr(df['s1'].shift(lag))
        correlations.append(corr)

    # Encontra o lag com a correlação máxima
    max_corr_idx = np.nanargmax(np.abs(correlations))
    optimal_lag = lags[max_corr_idx]
    max_corr = correlations[max_corr_idx]
    
    return optimal_lag, max_corr

# ============================================
# Carga de Dados e Processamento
# ============================================
print("\nListando e filtrando arquivos de dados...")
# Define as fatias de coordenadas com base no primeiro arquivo SST
try:
    primeiro_sst = glob.glob(os.path.join(sst_dir, sst_pattern_crw))[0]
    with xr.open_dataset(primeiro_sst) as ds:
        lat_s = get_lat_slice(ds, lat_min, lat_max)
        lon_s = slice(lon_min, lon_max)
except IndexError:
    print(f"ERRO CRÍTICO: Nenhum arquivo SST encontrado no padrão '{sst_pattern_crw}' em '{sst_dir}'.")
    exit()

# 1. Carrega o stack de referência (SST)
sst_stack = load_reference_stack(sst_pattern_crw, sst_dir, ['analysed_sst'], lat_s, lon_s, start_year, end_year)

# 2. Carrega e interpola a Clorofila para a grade do SST
chl_stack = None
if sst_stack is not None:
    chl_stack = load_and_interpolate_stack(chl_pattern, modis_dir, ['chlor_a'], start_year, end_year, reference_stack=sst_stack)
else:
    print("ERRO: Stack de SST não pôde ser carregado. Não é possível continuar com a Clorofila.")

# 3. Carrega os arquivos para cálculo de DLI
kd490_files = filter_files_by_year_range(glob.glob(os.path.join(modis_dir, kd490_pattern)), start_year, end_year)
par_files = filter_files_by_year_range(glob.glob(os.path.join(modis_dir, par_pattern)), start_year, end_year)

print("\nCalculando DLI Bentônico para cada site individualmente...")
dli_ts_by_site = calculate_dli_benthic_per_site(kd490_files, par_files, sites_df)

# ============================================
# Filtro de Data Preciso Pós-Carga
# ============================================
if use_specific_start_date:
    print(f"\n--- Aplicando filtro de data de início: {specific_start_date} ---")
    end_date_str = f"{end_year}-12-31" 
    if sst_stack is not None:
        sst_stack = sst_stack.sel(time=slice(specific_start_date, end_date_str))
    if chl_stack is not None:
        chl_stack = chl_stack.sel(time=slice(specific_start_date, end_date_str))
    
    # Filtra também as séries de DLI
    dli_ts_by_site_filtered = {}
    for site, ts in dli_ts_by_site.items():
        dli_ts_by_site_filtered[site] = ts.sel(time=slice(specific_start_date, end_date_str))
    dli_ts_by_site = dli_ts_by_site_filtered

# ============================================
# Geração de Séries Temporais Agrupadas e Análise de Lag
# ============================================

### ATUALIZADO ### Função principal de plotagem e análise para 3 variáveis
def plot_and_analyze_grouped_series(sst_data, chl_data, dli_data_by_site, sites_dataframe, output_dir):
    if sst_data is None or chl_data is None or not dli_data_by_site:
        print("Dados de SST, Chl-a ou DLI ausentes. Não é possível gerar o gráfico.")
        return

    grouped_sites = sites_dataframe.groupby('ARCH_SIMPLE')
    group_order = ['Inner Arch', 'Outer Arch']
    
    fig, axes = plt.subplots(nrows=len(group_order), ncols=1, figsize=(20, 8 * len(group_order)), sharex=True, squeeze=False)
    axes = axes.flatten()

    for i, arch_name in enumerate(group_order):
        ax = axes[i]
        try:
            group_df = grouped_sites.get_group(arch_name)
        except KeyError:
            print(f"AVISO: Nenhum site encontrado para o grupo '{arch_name}'. Pulando.")
            ax.set_title(f"Grupo: {arch_name} (Nenhum site)")
            ax.text(0.5, 0.5, 'Sem dados para este grupo', horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
            continue
            
        print(f"\n>>> Processando grupo '{arch_name}' com {len(group_df)} site(s)...")
        
        # --- Extração das Séries Médias para o Grupo ---
        site_lats = xr.DataArray(group_df['LAT'].values, dims="site")
        site_lons = xr.DataArray(group_df['LONG'].values, dims="site")

        mean_sst_ts = sst_data.sel(lat=site_lats, lon=site_lons, method='nearest').mean(dim='site', skipna=True)
        mean_chl_ts = chl_data.sel(lat=site_lats, lon=site_lons, method='nearest').mean(dim='site', skipna=True)
        
        group_site_names = group_df['SITE_NAME'].tolist()
        dli_series_list = [dli_data_by_site[name] for name in group_site_names if name in dli_data_by_site]
        
        if not dli_series_list: 
            print(f"AVISO: Nenhum dado de DLI para o grupo '{arch_name}'.")
            continue
        
        aligned_dli = xr.align(*dli_series_list, join='outer')
        mean_dli_ts = xr.concat(aligned_dli, dim='site').mean(dim='site', skipna=True)
        
        # --- Análise de Correlação Cruzada (sem alteração) ---
        print(f"--- Análise de Lag para o Grupo: {arch_name} ---")
        
        dli_pd = mean_dli_ts.to_pandas()
        sst_pd = mean_sst_ts.to_pandas()
        chl_pd = mean_chl_ts.to_pandas()
        
        pairs = {
            'DLI -> SST': (dli_pd, sst_pd),
            'SST -> Chl-a': (sst_pd, chl_pd),
            'DLI -> Chl-a': (dli_pd, chl_pd)
        }
        
        lag_results = {}
        for name, (s1, s2) in pairs.items():
            lag, corr = calculate_cross_correlation(s1, s2)
            lag_results[name] = (lag, corr)
            if not np.isnan(lag):
                print(f"{name}: Lag ótimo de {int(lag)} dias (Correlação = {corr:.2f})")
            else:
                print(f"{name}: Não foi possível calcular o lag (dados insuficientes).")
        
        # <--- INÍCIO DA MUDANÇA NO PLOT --->
        # --- Plotagem da Série Temporal (3 variáveis com 3 eixos Y, usando PONTOS) ---
        color_sst = 'tab:red'
        color_dli = 'tab:blue'
        color_chl = 'tab:green'
        
        # EIXO 1 (Esquerdo): SST
        ax.set_ylabel('SST Médio (°C)', color=color_sst, fontsize=14)
        # Alterado de '-' para 'o', adicionado markersize e linestyle='None'
        p1, = ax.plot(mean_sst_ts['time'], mean_sst_ts.values, 'o', markersize=4, alpha=0.7, color=color_sst, label=f'SST (Média {arch_name})', linestyle='None')
        ax.tick_params(axis='y', labelcolor=color_sst, labelsize=12)
        ax.grid(axis='y', linestyle='--', alpha=0.6, color=color_sst)
        
        # EIXO 2 (Direito 1): DLI
        ax2 = ax.twinx()
        ax2.set_ylabel('DLI Bentônico (mol/m²/d)', color=color_dli, fontsize=14)
        # Alterado de '-' para 'o', adicionado markersize e linestyle='None'
        p2, = ax2.plot(mean_dli_ts['time'], mean_dli_ts.values, 'o', markersize=4, alpha=0.7, color=color_dli, label=f'DLI (Média {arch_name})', linestyle='None')
        ax2.tick_params(axis='y', labelcolor=color_dli, labelsize=12)
        
        # EIXO 3 (Direito 2): Clorofila
        ax3 = ax.twinx()
        ax3.spines['right'].set_position(('outward', 70))
        ax3.set_ylabel('Chl-a (mg/m³)', color=color_chl, fontsize=14)
        # Alterado de '-' para 'o', adicionado markersize e linestyle='None'
        p3, = ax3.plot(mean_chl_ts['time'], mean_chl_ts.values, 'o', markersize=4, alpha=0.7, color=color_chl, label=f'Chl-a (Média {arch_name})', linestyle='None')
        ax3.tick_params(axis='y', labelcolor=color_chl, labelsize=12)
        # <--- FIM DA MUDANÇA NO PLOT --->

        # LEGENDA UNIFICADA (sem alteração)
        lines = [p1, p2, p3]
        ax.legend(lines, [l.get_label() for l in lines], loc='upper left', fontsize=12)
        
        # Título e Eixo X (sem alteração)
        lag_str1 = f"DLI→SST: {int(lag_results['DLI -> SST'][0])}d" if not np.isnan(lag_results['DLI -> SST'][0]) else "N/A"
        lag_str2 = f"SST→Chl: {int(lag_results['SST -> Chl-a'][0])}d" if not np.isnan(lag_results['SST -> Chl-a'][0]) else "N/A"
        lag_str3 = f"DLI→Chl: {int(lag_results['DLI -> Chl-a'][0])}d" if not np.isnan(lag_results['DLI -> Chl-a'][0]) else "N/A"
        
        ax.set_title(f"Grupo: {arch_name} (n={len(group_df)}) | Lags: {lag_str1}, {lag_str2}, {lag_str3}", fontsize=16, weight='bold')
        
        ax.xaxis.set_major_locator(mdates.YearLocator(1))
        ax.xaxis.set_minor_locator(mdates.MonthLocator(interval=3))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        plt.setp(ax.get_xticklabels(), rotation=0, ha="center", fontsize=12)

    fig.suptitle(f'Séries Temporais de SST, DLI e Chl-a por Arco ({start_year}-{end_year})', fontsize=20, y=1.02)
    fig.tight_layout(rect=[0, 0.03, 1, 0.97]) 
    main_plot_path = os.path.join(output_dir, f"time_series_SST_DLI_CHL_com_Lag_{start_year}-{end_year}.png")
    plt.savefig(main_plot_path, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"\nGráfico principal de séries temporais salvo em: {main_plot_path}")

# ============================================
# Execução da Plotagem e Análise
# ============================================
print("\nGerando gráficos e analisando o lag temporal para SST, DLI e Chl-a...")
plot_and_analyze_grouped_series(sst_stack, chl_stack, dli_ts_by_site, sites_df, output_dir)

print("\nScript finalizado com sucesso!")