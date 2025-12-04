### PCA LOCAL com bubbleplot - VERSÃO AVANÇADA 2.3 (Estilo de Arco Aprimorado) ###
# 1. Realiza duas PCAs separadas: Variabilidade (CVs) e Magnitude (Médias + DHW).
# 2. Diferencia arcos com bordas distintas (grossa/preta para inner, fina/cinza para outer).
# 3. Adiciona uma legenda explicando o estilo das bordas dos arcos.
# 4. Gera resumo de segmentos e salva scores de AMBAS as PCAs.
# 5. Funciona perfeitamente com analysis_selector = 'all' ou numérico.
# 6. Figuras PCA compostas

# ==============================
# Imports
# ==============================
import os
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
os.environ['OMP_NUM_THREADS'] = '1'
import glob
import re
import logging
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import dask
import psutil
import gc

# ==============================
# Configurações Iniciais
# ==============================
logging.basicConfig(
    filename='pca_local_errors.log',
    filemode='w',
    level=logging.ERROR,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# ==============================
# Caminhos e parâmetros principais
# ==============================
output_dir = r"C:\Users\rbfra\OneDrive\########CEBIMAR\###Orientacoes e Supervisoes\###Mestrado Mariana\PROJETO FAPESP\R1_plus_v2"
os.makedirs(output_dir, exist_ok=True)

# ==============================
# Configurações do Dask
# ==============================
dask.config.set({
    'array.slicing.split_large_chunks': True,
    'scheduler': 'single-threaded',
    'temporary_directory': os.path.join(output_dir, 'dask_temp')
})
os.makedirs(os.path.join(output_dir, 'dask_temp'), exist_ok=True)

# ==============================
# Parâmetros da Análise
# ==============================
analysis_selector = 'all'

print(f"============================================================")
print(f"==> ANÁLISE CONFIGURADA PARA A OPÇÃO DE CV: '{analysis_selector}' <==")
print(f"============================================================")

# --- Períodos e Caminhos ---
sst_dir   = r"E:\remote sensing\CRW_SST_FULL"
modis_dir = r"E:\remote sensing\MODIS_DATA_FULL"
chl_dir   = r"E:\remote sensing\MODIS_DATA_FULL"
dhw_dir   = r"E:\remote sensing\CRW_DHW_FULL"

sst_period   = (2020, 2024)
dhw_period   = (2020, 2024)
light_period = (2020, 2024)
chl_period   = (2020, 2024)

#sst_period   = (2007, 2008)
#dhw_period   = (2007, 2008)
#light_period = (2007, 2008)
#chl_period   = (2007, 2008)


sites_csv_file = r"C:\Users\rbfra\OneDrive\########CEBIMAR\####PROJETOS\#####Coral trade offs\sites_list_full.csv"
lat_min, lat_max = -20.5, -14.5
lon_min, lon_max = -40.5, -35.5
sst_pattern   = 'coraltemp_v3.1_*.nc'
dhw_pattern   = '*.nc'
kd490_pattern = 'AQUA_MODIS.*.L3m.DAY.KD.Kd_490.4km.nc'
par_pattern   = 'AQUA_MODIS.*.L3m.DAY.PAR.par.4km.nc'
chl_pattern   = 'AQUA_MODIS.*.L3m.DAY.CHL.chlor_a.4km.nc'
KDPAR_GATTUSO_A, KDPAR_GATTUSO_B, KDPAR_GATTUSO_C = 0.0665, 0.874, 0.00121
KDPAR_KD490_MIN_THRESHOLD = 0.001

# ==============================
# Carregamento da lista de sítios
# ==============================
try:
    df_sites_info = pd.read_csv(sites_csv_file, sep=';')
    df_sites_info.columns = [col.strip().replace(' ', '_') for col in df_sites_info.columns]
    for col in df_sites_info.select_dtypes(['object']).columns:
        df_sites_info[col] = df_sites_info[col].str.strip()
    df_sites_info['unique_id'] = df_sites_info['Site_name'] + '_' + df_sites_info['HAB']
    sites = list(zip(df_sites_info['Latitude'], df_sites_info['Longitude']))
    print("DataFrame de sítios carregado com sucesso.")
except Exception as e:
    logging.critical(f"FALHA CRÍTICA AO CARREGAR O ARQUIVO DE SÍTIOS: {e}")
    raise

# ==============================
# Funções Auxiliares Essenciais
# ==============================
# Funções de carregamento de dados e cálculo de métricas permanecem as mesmas da versão anterior.
# Omitidas para brevidade.
def verify_file_integrity(filepath):
    if not os.path.exists(filepath) or os.path.getsize(filepath) == 0: return False
    try:
        with xr.open_dataset(filepath) as ds: _ = ds.attrs
        return True
    except Exception: return False
def load_satellite_data(pattern, base_dir, period, var_name_options, chunks={'time': 30, 'lat': 100, 'lon': 100}):
    full_search_path = os.path.join(base_dir, pattern); print(f"\n--- Carregando dados para: {pattern} em {period} ---")
    all_files = sorted(glob.glob(full_search_path));
    if not all_files: logging.error(f"Nenhum arquivo encontrado para '{pattern}' em {base_dir}"); return None
    candidate_paths = [f for f in all_files if (match := re.search(r'(\d{4})\d{4}', os.path.basename(f))) and period[0] <= int(match.group(1)) <= period[1]]
    if not candidate_paths: logging.warning(f"Nenhum arquivo encontrado para '{pattern}' no período {period}"); return None
    valid_paths = [f for f in candidate_paths if verify_file_integrity(f)]
    if not valid_paths: logging.error(f"Nenhum arquivo VÁLIDO encontrado para '{pattern}' no período {period}"); return None
    print(f"Encontrados {len(valid_paths)} arquivos válidos.")
    def preprocess_with_time(ds):
        var_name = next((v for v in var_name_options if v in ds.data_vars), None)
        if var_name is None: return xr.Dataset()
        ds_subset = ds[[var_name]].astype('float32')
        if 'time' not in ds_subset.coords:
            match = re.search(r'(\d{8})', os.path.basename(ds.encoding.get("source", "")))
            if match: dt = pd.to_datetime(match.group(1), format='%Y%m%d'); return ds_subset.expand_dims(time=[dt])
        return ds_subset
    try:
        ds_raw = xr.open_mfdataset(valid_paths, preprocess=preprocess_with_time, combine='by_coords', parallel=True, chunks=chunks, engine='netcdf4')
        if not ds_raw.data_vars: print(f"AVISO: Nenhum dado carregado para '{pattern}'."); return None
        lat_slice = slice(lat_max, lat_min) if ds_raw['lat'].values[0] > ds_raw['lat'].values[-1] else slice(lat_min, lat_max)
        ds_final = ds_raw.sel(lat=lat_slice, lon=slice(lon_min, lon_max))
        return ds_final.sortby('time')[list(ds_final.data_vars)[0]]
    except Exception as e: logging.critical(f"Falha ao concatenar arquivos para '{pattern}': {e}"); return None
def calculate_kdpar_gattuso(kd490_da):
    kd490_safe = xr.where((kd490_da.isnull()) | (kd490_da <= KDPAR_KD490_MIN_THRESHOLD), np.nan, kd490_da)
    epsilon = 1e-9; kdpar = KDPAR_GATTUSO_A + KDPAR_GATTUSO_B * kd490_safe - KDPAR_GATTUSO_C / (kd490_safe + epsilon)
    kdpar_final = xr.where((kdpar.isnull()) | (kdpar <= 0), np.nan, kdpar); kdpar_final.name = "kdpar"
    return kdpar_final
def extract_and_compute_site_metrics(da, sites_coords, metrics_to_calc):
    if da is None: return {}
    if isinstance(da, xr.Dataset): da = da[list(da.data_vars)[0]]
    lats = xr.DataArray([s[0] for s in sites_coords], dims="site"); lons = xr.DataArray([s[1] for s in sites_coords], dims="site")
    site_timeseries = da.sel(lat=lats, lon=lons, method='nearest').load()
    results = {}
    for metric in metrics_to_calc:
        print(f"  Calculando métrica: {metric}...")
        if metric == 'mean':
            results['mean'] = site_timeseries.mean('time', skipna=True).values
        elif metric == 'cv_all':
            mean_vals = site_timeseries.mean('time', skipna=True); std_vals = site_timeseries.std('time', skipna=True)
            results['cv_all'] = (std_vals / (mean_vals + 1e-9) * 100).values
            results['cv_all_segment_count'] = site_timeseries.notnull().sum('time').values
        elif metric.startswith('cv_'):
            try:
                window = int(metric.split('_')[1]); min_p = max(2, int(window * 0.25))
                rolling_mean = site_timeseries.rolling(time=window, min_periods=min_p, center=True).mean()
                rolling_std = site_timeseries.rolling(time=window, min_periods=min_p, center=True).std()
                cv_ts = (rolling_std / (rolling_mean + 1e-9)) * 100
                results[metric] = cv_ts.mean('time', skipna=True).values
                results[f'{metric}_segment_count'] = cv_ts.notnull().sum('time').values
            except (ValueError, IndexError): print(f"    AVISO: Não foi possível processar a métrica '{metric}'. Pulando."); continue
        elif metric.startswith('prop_gt_'):
            threshold = float(metric.split('_')[-1])
            results[metric] = (site_timeseries > threshold).mean('time', skipna=True).values
    return results
def calculate_site_specific_dli(kd490_da, par_da, sites_info_df, cv_selector):
    if kd490_da is None or par_da is None: return {}
    lats = xr.DataArray(sites_info_df['Latitude'].values, dims="site"); lons = xr.DataArray(sites_info_df['Longitude'].values, dims="site")
    depths = xr.DataArray(sites_info_df['Depth_m'].values, dims="site")
    kd490_points = kd490_da.sel(lat=lats, lon=lons, method='nearest').load()
    par_points = par_da.sel(lat=lats, lon=lons, method='nearest').load()
    kdpar_points = calculate_kdpar_gattuso(kd490_points)
    benthic_dli_timeseries = par_points * np.exp(-kdpar_points * depths)
    results = {}; mean_dli = benthic_dli_timeseries.mean('time', skipna=True); results['mean_DLI_local'] = mean_dli.values
    if cv_selector == 'all':
        std_dli = benthic_dli_timeseries.std('time', skipna=True)
        results['cv_DLI_local'] = (std_dli / (mean_dli + 1e-9) * 100).values
        results['cv_DLI_local_segment_count'] = benthic_dli_timeseries.notnull().sum('time').values
    else:
        window = int(cv_selector); min_p = max(2, int(window * 0.25))
        rolling_mean = benthic_dli_timeseries.rolling(time=window, min_periods=min_p, center=True).mean()
        rolling_std = benthic_dli_timeseries.rolling(time=window, min_periods=min_p, center=True).std()
        cv_ts = (rolling_std / (rolling_mean + 1e-9)) * 100
        results[f'dli_cv_{window}'] = cv_ts.mean('time', skipna=True).values
        results[f'dli_cv_{window}_segment_count'] = cv_ts.notnull().sum('time').values
    return results

# ==============================
# Processamento Principal
# ==============================
print("\n--- 1. Carregando Dados de Satélite ---")
sst_all_ds = load_satellite_data(sst_pattern, sst_dir, sst_period, ['analysed_sst', 'sea_surface_temperature', 'sst'])
dhw_all_ds = load_satellite_data(dhw_pattern, dhw_dir, dhw_period, ['degree_heating_week'])
chl_all_ds = load_satellite_data(chl_pattern, chl_dir, chl_period, ['chlor_a'])
kd490_all_ds = load_satellite_data(kd490_pattern, modis_dir, light_period, ['Kd_490'])
par_all_ds = load_satellite_data(par_pattern, modis_dir, light_period, ['par'])

print("\n--- 2. Calculando Métricas para os Sítios ---")
df_sites_local = df_sites_info.copy()

# DLI
dli_metrics = calculate_site_specific_dli(kd490_all_ds, par_all_ds, df_sites_info, analysis_selector)
for key, values in dli_metrics.items(): df_sites_local[key] = values
del kd490_all_ds, par_all_ds, dli_metrics; gc.collect()

# SST
sst_metrics_to_calculate = ['mean']; sst_metrics_to_calculate.append('cv_all' if str(analysis_selector) == 'all' else f'cv_{analysis_selector}')
sst_metrics = extract_and_compute_site_metrics(sst_all_ds, sites, sst_metrics_to_calculate)
for key, values in sst_metrics.items(): df_sites_local[f'sst_{key}'] = values
del sst_all_ds, sst_metrics; gc.collect()

# DHW
dhw_metrics = extract_and_compute_site_metrics(dhw_all_ds, sites, ['prop_gt_4', 'prop_gt_8'])
if 'prop_gt_4' in dhw_metrics:
    df_sites_local['prop_DHW_gt4'] = dhw_metrics['prop_gt_4']
if 'prop_gt_8' in dhw_metrics:
    df_sites_local['prop_DHW_gt8'] = dhw_metrics['prop_gt_8']
del dhw_all_ds, dhw_metrics; gc.collect()

# CHL
chl_metrics_to_calculate = ['mean']; chl_metrics_to_calculate.append('cv_all' if str(analysis_selector) == 'all' else f'cv_{analysis_selector}')
chl_metrics = extract_and_compute_site_metrics(chl_all_ds, sites, chl_metrics_to_calculate)
for key, values in chl_metrics.items(): df_sites_local[f'chl_{key}'] = values
del chl_all_ds, chl_metrics; gc.collect()

# =============================================================================
# --- SEÇÃO 3: FUNÇÃO GERAL PARA PCA E PLOTAGEM (MODIFICADA) ---
# =============================================================================
def run_pca_and_generate_bubble_plots(df_input, pca_vars_list, output_suffix):
    """
    Executa uma análise PCA, gera os BUBBLE PLOTS individuais e RETORNA
    os dados necessários para a figura composta (scores, loadings, etc.).
    """
    print("\n" + "="*50)
    print(f"INICIANDO PCA: {output_suffix.replace('_', ' ').strip().title()}")
    print("="*50)

    existing_vars = [var for var in pca_vars_list if var in df_input.columns]
    if len(existing_vars) < 2:
        print(f"AVISO: PCA para '{output_suffix}' pulada. Variáveis encontradas: {existing_vars}")
        return df_input, None # Retorna None para indicar falha

    print(f"Variáveis usadas na PCA: {existing_vars}")
    df_subset = df_input.dropna(subset=existing_vars) # Remove linhas com NaNs antes da PCA
    if len(df_subset) < 2:
        print(f"AVISO: PCA para '{output_suffix}' pulada. Pontos de dados insuficientes após remover NaNs.")
        return df_input, None
    
    data_for_pca = df_subset[existing_vars].copy()
    transformed_names = list(data_for_pca.columns)
    
    # Transformações
    for i, col_name in enumerate(data_for_pca.columns):
        if col_name == 'prop_DHW_gt4':
            data_for_pca[col_name] = np.sqrt(data_for_pca[col_name] + 0.5)
            transformed_names[i] = "sqrt(DHW>4)"
        elif col_name in ['chl_mean', 'sst_mean', 'mean_DLI_local']:
            data_for_pca[col_name] = np.log1p(data_for_pca[col_name])
            transformed_names[i] = f"log({col_name}+1)"

    # PCA
    scaler = StandardScaler()
    data_matrix = scaler.fit_transform(data_for_pca)
    pca = PCA(n_components=2)
    scores = pca.fit_transform(data_matrix)
    explained_variance = pca.explained_variance_ratio_ * 100
    
    # Salva loadings
    loadings_df = pd.DataFrame(pca.components_.T, columns=["PC1", "PC2"], index=transformed_names)
    loadings_path = os.path.join(output_dir, f"loadings_PCA{output_suffix}.csv")
    loadings_df.to_csv(loadings_path)
    print(f"Loadings salvos em {loadings_path}")

    # Adiciona scores ao DataFrame original (usando o índice do df_subset)
    df_scores = df_input.copy()
    df_scores.loc[df_subset.index, f'PC1{output_suffix}'] = scores[:, 0]
    df_scores.loc[df_subset.index, f'PC2{output_suffix}'] = scores[:, 1]
    
    # --- Geração dos BUBBLE PLOTS (individualmente) ---
    print(f"\n--- Gerando Bubble Plots para PCA {output_suffix} ---")
    unique_reefs = df_scores['Reef_name'].unique()
    cmap = plt.get_cmap('tab10')
    color_map = {r: cmap(i % 10) for i, r in enumerate(unique_reefs)}
    unique_habitats = df_scores['HAB'].unique()
    habitat_shapes = ['o', 's', '^', 'D', 'v', '<', '>']
    shape_map = {hab: habitat_shapes[i % len(habitat_shapes)] for i, hab in enumerate(unique_habitats)}
    arch_styles = {'inner': {'edgecolor': 'black', 'linewidth': 2.0}, 'outer': {'edgecolor': 'darkgrey', 'linewidth': 0.75}}

    def scale_marker_sizes(values, scale_factor=300, min_size=30):
        # ... (código da função scale_marker_sizes)
        vals = np.asarray(values, dtype=float); valid = np.isfinite(vals)
        if valid.sum() == 0: return np.full(vals.shape, min_size)
        vmin, vmax = vals[valid].min(), vals[valid].max()
        if vmax == vmin: return np.full(vals.shape, min_size + scale_factor / 2)
        scaled = (vals - vmin) / (vmax - vmin); scaled[~valid] = 0
        return scaled * scale_factor + min_size
        
    for var_name in existing_vars:
        site_values = df_scores[var_name].values
        if np.all(np.isnan(site_values)): continue
        sizes = scale_marker_sizes(site_values)
        fig, ax = plt.subplots(figsize=(12, 8))
        for i, row in df_scores.iterrows():
            if pd.notna(row[f'PC1{output_suffix}']):
                style = arch_styles.get(row['Arch'].lower(), {'edgecolor': 'grey', 'linewidth': 0.5})
                ax.scatter(row[f'PC1{output_suffix}'], row[f'PC2{output_suffix}'], s=sizes[i],
                           color=color_map.get(row['Reef_name'], 'grey'), marker=shape_map.get(row['HAB'], 'x'),
                           alpha=0.7, **style)
        # ... (código para adicionar legendas ao bubble plot, omitido para brevidade)
        ax.set_title(f"PCA {output_suffix.replace('_', ' ')} - Tamanho ∝ {var_name}")
        # ... (código para salvar o bubble plot)
        plt.close(fig)
    
    # Prepara o dicionário de resultados para a plotagem composta
    pca_results = {
        'df_scores': df_scores,
        'loadings': loadings_df,
        'explained_variance': explained_variance,
        'pc1_col': f'PC1{output_suffix}',
        'pc2_col': f'PC2{output_suffix}',
        'title_suffix': output_suffix.replace('_', ' ').strip()
    }
    
    return df_scores, pca_results

# =============================================================================
# --- SEÇÃO 4: EXECUÇÃO DAS DUAS ANÁLISES PCA ---
# =============================================================================
df_processed = df_sites_local.copy()

# Define as variáveis de CV com base no seletor
if str(analysis_selector).lower() == 'all':
    cv_sst_col, cv_dli_col, cv_chl_col = 'sst_cv_all', 'cv_DLI_local', 'chl_cv_all'
else:
    cv_sst_col = f'sst_cv_{analysis_selector}'
    cv_dli_col = f'dli_cv_{analysis_selector}'
    cv_chl_col = f'chl_cv_{analysis_selector}'

# Executa PCA de Variabilidade
pca_vars_variability = [cv_sst_col, cv_dli_col, cv_chl_col]
df_processed, variability_results = run_pca_and_generate_bubble_plots(df_processed, pca_vars_variability, "_Variability")

# Executa PCA de Magnitude
pca_vars_magnitude = ['sst_mean', 'mean_DLI_local', 'chl_mean', 'prop_DHW_gt4']
df_processed, magnitude_results = run_pca_and_generate_bubble_plots(df_processed, pca_vars_magnitude, "_Magnitude")

# Salva o arquivo de dados consolidado
final_scores_path = os.path.join(output_dir, "dados_consolidados_com_scores_das_duas_PCAs.xlsx")
df_processed.to_excel(final_scores_path, index=False)
print(f"\nArquivo final com todas as métricas e scores salvo em: {final_scores_path}")

# =============================================================================
# --- SEÇÃO 5: GERAÇÃO DE TODAS AS FIGURAS COMPOSTAS (COM LEGENDAS) ---
# =============================================================================

# --- Função Auxiliar 1: Para a figura de Ordenação + Loadings (COM LEGENDAS) ---
def create_composite_pca_figure(pca_results, output_filename):
    """
    Cria uma figura composta (1 linha, 3 colunas) para os resultados de uma PCA.
    Col 1: Ordenação com legendas | Col 2: Loadings | Col 3: Scree Plot
    """
    if pca_results is None:
        print(f"Não há resultados de PCA para gerar a figura {output_filename}. Pulando.")
        return

    df_scores = pca_results['df_scores']
    loadings = pca_results['loadings']
    explained_variance = pca_results['explained_variance']
    pc1_col, pc2_col = pca_results['pc1_col'], pca_results['pc2_col']
    title_suffix = pca_results['title_suffix']
    
    fig, axes = plt.subplots(1, 3, figsize=(24, 7), gridspec_kw={'width_ratios': [1.2, 1, 0.8]})
    fig.suptitle(f'Resumo da Análise de Componentes Principais - {title_suffix}', fontsize=20, y=1.02)

    # --- Coluna 1: Gráfico de Ordenação ---
    ax1 = axes[0]
    unique_reefs = df_scores['Reef_name'].unique(); cmap = plt.get_cmap('tab10'); color_map = {r: cmap(i % 10) for i, r in enumerate(unique_reefs)}
    unique_habitats = df_scores['HAB'].unique(); habitat_shapes = ['o', 's', '^', 'D', 'v', '<', '>']; shape_map = {hab: habitat_shapes[i % len(habitat_shapes)] for i, hab in enumerate(unique_habitats)}
    arch_styles = {'inner': {'edgecolor': 'black', 'linewidth': 2.0}, 'outer': {'edgecolor': 'darkgrey', 'linewidth': 0.75}}

    for _, row in df_scores.iterrows():
        if pd.notna(row[pc1_col]):
            style = arch_styles.get(row['Arch'].lower(), {'edgecolor': 'grey', 'linewidth': 0.5})
            ax1.scatter(row[pc1_col], row[pc2_col], color=color_map.get(row['Reef_name']), marker=shape_map.get(row['HAB']), s=100, alpha=0.8, **style)

    ax1.set_xlabel(f"PC1 ({explained_variance[0]:.1f}%)", fontsize=14)
    ax1.set_ylabel(f"PC2 ({explained_variance[1]:.1f}%)", fontsize=14)
    ax1.set_title("Ordenação dos Sítios", fontsize=16)
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.axhline(0, color='grey', lw=0.5); ax1.axvline(0, color='grey', lw=0.5)

    # <--- INÍCIO DA ADIÇÃO DAS LEGENDAS --->
    # Cria os elementos para cada legenda
    legend_elements_color = [plt.Line2D([0], [0], marker='o', color='w', label=reef, markersize=10, markerfacecolor=color_map[reef]) for reef in unique_reefs]
    legend_elements_shape = [plt.Line2D([0], [0], marker=shape_map[hab], color='grey', label=hab, linestyle='None', markersize=10) for hab in unique_habitats]
    legend_elements_arch = [
        plt.Line2D([0], [0], marker='o', color='w', label='Inner Arc', markersize=10, markeredgecolor=arch_styles['inner']['edgecolor'], markeredgewidth=arch_styles['inner']['linewidth']),
        plt.Line2D([0], [0], marker='o', color='w', label='Outer Arc', markersize=10, markeredgecolor=arch_styles['outer']['edgecolor'], markeredgewidth=arch_styles['outer']['linewidth'])
    ]
    # Adiciona as legendas à FIGURA, não ao eixo, para posicionamento global
    leg1 = fig.legend(title="Recife", handles=legend_elements_color, loc='center left', bbox_to_anchor=(0.91, 0.75))
    leg2 = fig.legend(title="Habitat", handles=legend_elements_shape, loc='center left', bbox_to_anchor=(0.91, 0.5))
    leg3 = fig.legend(title="Arco", handles=legend_elements_arch, loc='center left', bbox_to_anchor=(0.91, 0.25))
    # <--- FIM DA ADIÇÃO DAS LEGENDAS --->

    # --- Coluna 2: Gráfico de Loadings ---
    ax2 = axes[1]
    # ... (código do gráfico de loadings, sem alterações)
    ax2.axhline(0, color='grey', lw=0.5); ax2.axvline(0, color='grey', lw=0.5)
    for i, var in enumerate(loadings.index):
        ax2.arrow(0, 0, loadings['PC1'][i]*2, loadings['PC2'][i]*2, head_width=0.05, head_length=0.1, fc='red', ec='red')
        ax2.text(loadings['PC1'][i]*2.2, loadings['PC2'][i]*2.2, var, color='black', ha='center', va='center', fontsize=12)
    ax2.set_xlim(-2.5, 2.5); ax2.set_ylim(-2.5, 2.5)
    ax2.set_xlabel("Contribuição para PC1", fontsize=14)
    ax2.set_ylabel("Contribuição para PC2", fontsize=14)
    ax2.set_title("Loadings das Variáveis", fontsize=16)
    ax2.set_aspect('equal', adjustable='box')

    # --- Coluna 3: Scree Plot ---
    ax3 = axes[2]
    # ... (código do scree plot, sem alterações)
    components = ['PC1', 'PC2']
    ax3.bar(components, explained_variance, color='skyblue', edgecolor='black')
    ax3.set_ylabel("Variância Explicada (%)", fontsize=14)
    ax3.set_title("Importância dos Componentes", fontsize=16)
    ax3.set_ylim(0, 100)
    for i, v in enumerate(explained_variance):
        ax3.text(i, v + 2, f"{v:.1f}%", ha='center', color='black', fontsize=12)

    # Ajusta o layout para criar espaço para as legendas à direita
    fig.subplots_adjust(right=0.9)
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Figura composta de resumo salva em: {output_filename}")


# --- Função Auxiliar 2: Para a figura de Bubble Plots (COM LEGENDAS) ---
def create_composite_bubble_plot_figure(pca_results, pca_vars, output_filename, nrows, ncols):
    """
    Cria uma figura composta com os bubble plots das variáveis da PCA em um grid customizado.
    """
    if pca_results is None:
        print(f"Não há resultados de PCA para gerar a figura {output_filename}. Pulando.")
        return

    df_scores = pca_results['df_scores']
    explained_variance = pca_results['explained_variance']
    pc1_col, pc2_col = pca_results['pc1_col'], pca_results['pc2_col']
    title_suffix = pca_results['title_suffix']
    
    # Ajusta o tamanho da figura para ter espaço extra para a legenda
    fig_width = (6 * ncols) + 3 # 6 por subplot + 3 de espaço para a legenda
    fig_height = 6 * nrows
    fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height), sharex=True, sharey=True)
    
    if nrows == 1 and ncols == 1: axes = np.array([[axes]])
    elif nrows == 1 or ncols == 1: axes = np.array([axes]).flatten().reshape(nrows, ncols)

    fig.suptitle(f'Bubble Plots - PCA de {title_suffix}', fontsize=20, y=0.98)

    # Definições de estilo e legenda
    unique_reefs = df_scores['Reef_name'].unique(); cmap = plt.get_cmap('tab10'); color_map = {r: cmap(i % 10) for i, r in enumerate(unique_reefs)}
    unique_habitats = df_scores['HAB'].unique(); habitat_shapes = ['o', 's', '^', 'D', 'v', '<', '>']; shape_map = {hab: habitat_shapes[i % len(habitat_shapes)] for i, hab in enumerate(unique_habitats)}
    arch_styles = {'inner': {'edgecolor': 'black', 'linewidth': 2.0}, 'outer': {'edgecolor': 'darkgrey', 'linewidth': 0.75}}

    def scale_marker_sizes(values, scale_factor=300, min_size=30):
        # ... (código da função scale_marker_sizes, sem alterações)
        vals = np.asarray(values, dtype=float); valid = np.isfinite(vals)
        if valid.sum() == 0: return np.full(vals.shape, min_size)
        vmin, vmax = vals[valid].min(), vals[valid].max()
        if vmax == vmin: return np.full(vals.shape, min_size + scale_factor / 2)
        scaled = (vals - vmin) / (vmax - vmin); scaled[~valid] = 0
        return scaled * scale_factor + min_size
        
    for i, var_name in enumerate(pca_vars):
        row_idx, col_idx = divmod(i, ncols)
        ax = axes[row_idx, col_idx]
        
        site_values = df_scores[var_name].values
        sizes = scale_marker_sizes(site_values)
        for _, row in df_scores.iterrows():
            if pd.notna(row[pc1_col]):
                style = arch_styles.get(row['Arch'].lower(), {'edgecolor': 'grey', 'linewidth': 0.5})
                ax.scatter(row[pc1_col], row[pc2_col], s=sizes[_],
                           color=color_map.get(row['Reef_name']), marker=shape_map.get(row['HAB']),
                           alpha=0.7, **style)

        ax.set_title(f"Tamanho ∝ {var_name}", fontsize=16)
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.axhline(0, color='grey', lw=0.5); ax.axvline(0, color='grey', lw=0.5)

    for j in range(len(pca_vars), nrows * ncols):
        row_idx, col_idx = divmod(j, ncols)
        axes[row_idx, col_idx].set_visible(False)

    # Rótulos dos eixos principais
    fig.text(0.5, 0.04, f"PC1 ({explained_variance[0]:.1f}%)", ha='center', va='center', fontsize=14)
    fig.text(0.04, 0.5, f"PC2 ({explained_variance[1]:.1f}%)", ha='center', va='center', rotation='vertical', fontsize=14)
    
    # <--- ADIÇÃO DAS LEGENDAS (mesma lógica da outra função) --->
    legend_elements_color = [plt.Line2D([0], [0], marker='o', color='w', label=reef, markersize=10, markerfacecolor=color_map[reef]) for reef in unique_reefs]
    legend_elements_shape = [plt.Line2D([0], [0], marker=shape_map[hab], color='grey', label=hab, linestyle='None', markersize=10) for hab in unique_habitats]
    legend_elements_arch = [
        plt.Line2D([0], [0], marker='o', color='w', label='Inner Arc', markersize=10, markeredgecolor=arch_styles['inner']['edgecolor'], markeredgewidth=arch_styles['inner']['linewidth']),
        plt.Line2D([0], [0], marker='o', color='w', label='Outer Arc', markersize=10, markeredgecolor=arch_styles['outer']['edgecolor'], markeredgewidth=arch_styles['outer']['linewidth'])
    ]
    # Calcula a posição da legenda com base no número de colunas
    legend_x_pos = 1 - (1.5 / fig_width)
    leg1 = fig.legend(title="Recife", handles=legend_elements_color, loc='center left', bbox_to_anchor=(legend_x_pos, 0.75))
    leg2 = fig.legend(title="Habitat", handles=legend_elements_shape, loc='center left', bbox_to_anchor=(legend_x_pos, 0.5))
    leg3 = fig.legend(title="Arco", handles=legend_elements_arch, loc='center left', bbox_to_anchor=(legend_x_pos, 0.25))

    # Ajusta o layout para criar espaço para a legenda
    fig.subplots_adjust(right=1 - (3.0 / fig_width))
    plt.savefig(output_filename, dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f"Figura composta de bubble plots salva em: {output_filename}")


# --- Gerar TODAS as figuras finais ---
print("\n--- Gerando Figuras Compostas Finais ---")

# Figuras de Resumo (Ordenação + Loadings)
create_composite_pca_figure(magnitude_results, os.path.join(output_dir, "figura_composta_PCA_Magnitude.png"))
create_composite_pca_figure(variability_results, os.path.join(output_dir, "figura_composta_PCA_Variability.png"))

# Figuras de Bubble Plots - com layout customizado
create_composite_bubble_plot_figure(
    pca_results=magnitude_results, 
    pca_vars=pca_vars_magnitude, 
    output_filename=os.path.join(output_dir, "figura_composta_bubbles_Magnitude.png"),
    nrows=2, 
    ncols=2
)
create_composite_bubble_plot_figure(
    pca_results=variability_results, 
    pca_vars=pca_vars_variability, 
    output_filename=os.path.join(output_dir, "figura_composta_bubbles_Variability.png"),
    nrows=3, 
    ncols=1
)


# --- SEÇÃO 6: Salvar Resumo de Segmentos ---
print("\n--- Gerando arquivo de resumo de segmentos (janelas válidas) ---")
# ... (código desta seção permanece o mesmo)
segment_cols = [col for col in df_processed.columns if 'segment_count' in col]
summary_cols = ['Site_name', 'HAB', 'Arch'] + segment_cols
if segment_cols:
    df_segments_summary = df_processed[summary_cols]
    summary_path = os.path.join(output_dir, "resumo_segmentos_por_site.xlsx")
    df_segments_summary.to_excel(summary_path, index=False)
    print(f"Resumo de segmentos por site salvo em: {summary_path}")
else:
    print("Nenhuma coluna de contagem de segmentos encontrada para criar o resumo.")

print("\n\nProcesso completo finalizado com sucesso!")