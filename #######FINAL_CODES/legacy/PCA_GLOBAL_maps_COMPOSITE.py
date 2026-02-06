### PCA GLOBAL COMPOSITE - Publication Quality Maps (Nature/Science)
### DLI Bentônico (Daily Light Integral) 
### Versão com Figura Composta e Caminhos Atualizados (K:)

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
import geopandas as gpd
from shapely.geometry import box, LineString
from shapely.ops import unary_union, polygonize
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from dask.distributed import Client
import dask
import psutil
import gc

# --- Imports adicionais para figuras de alta qualidade ---
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
import matplotlib.patches as mpatches

# Fallback: Tentativa de usar Cartopy para projeções mais robustas
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAS_CARTOPY = True
except ImportError:
    HAS_CARTOPY = False
    print("Cartopy não encontrado. Usando matplotlib padrão.")

# ==============================
# Configurações de Estilo para Publicação (Nature/Science)
# ==============================
PUBLICATION_STYLE = {
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9,
    'axes.labelsize': 10,
    'axes.titlesize': 11,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'legend.fontsize': 8,
    'figure.dpi': 150,  # Para exibição na tela
    'savefig.dpi': 600, # Altíssima resolução para impressão
    'axes.linewidth': 0.8,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
}
plt.rcParams.update(PUBLICATION_STYLE)

PANEL_LABELS = ['a', 'b', 'c', 'd', 'e', 'f']
# Dimensões sugeridas: 180mm de largura (página inteira)
FIGSIZE = (180 / 25.4, 140 / 25.4) 

# ==============================
# Configurações Iniciais e Caminhos
# ==============================
logging.basicConfig(
    filename='file_open_errors_global_pca_composite.log',
    filemode='w',
    level=logging.ERROR,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# NOVO DIRETÓRIO DE OUTPUS (DEDICADO EM FINAL_RESULTS)
output_dir = r"C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_RESULTS\PCA_GLOBAL_COMPOSITE"
os.makedirs(output_dir, exist_ok=True)

dask.config.set({
    'array.slicing.split_large_chunks': True,
    'distributed.worker.memory.target': 0.6,
    'distributed.worker.memory.spill': 0.7,
    'distributed.worker.memory.pause': 0.9,
    'distributed.worker.memory.terminate': 0.95,
    'scheduler': 'single-threaded',
    'temporary_directory': os.path.join(output_dir, 'dask_temp')
})
os.makedirs(os.path.join(output_dir, 'dask_temp'), exist_ok=True)

# --- SELETOR ÚNICO DE ANÁLISE ---
analysis_selector = 'all'  

print(f"============================================================")
print(f"==> ANÁLISE CONFIGURADA PARA A OPÇÃO DE CV: '{analysis_selector}' <==")
print(f"============================================================")

# ATUALIZADO PARA DRIVE K:
sst_dir   = r"K:\remote sensing\CRW_SST_FULL"
modis_dir = r"K:\remote sensing\MODIS_DATA_FULL"
chl_dir   = r"K:\remote sensing\MODIS_DATA_FULL"
dhw_dir   = r"K:\remote sensing\CRW_DHW_FULL"

sst_period   = (1985, 2008)
dhw_period   = (1985, 2008)
light_period = (2002, 2008)
chl_period   = (2002, 2008)

bathymetry_file      = r"C:\Users\rbfra\OneDrive\GIS shapes\gebco_2024_ASO.nc"
# ATUALIZADO PARA DRIVE K:
protected_areas_shp  = r"K:\remote sensing\Mascaras\uc_fed_agosto_2016_site_shp\uc_fed_agosto_2016_site.shp"
coast_islands_dir    = r"C:\Users\rbfra\OneDrive\GIS shapes\batimetria_new"
# ATUALIZADO PARA DRIVE K:
reefs_shp            = r"K:\remote sensing\Mascaras\Recifes_Banco_dos_Abrolhos\Recifes_Banco_dos_Abrolhos.shp"
sites_csv_file = r"C:\Users\rbfra\OneDrive\########CEBIMAR\####PROJETOS\#####Coral trade offs\sites_list_full_clean.csv"

coast_shp   = os.path.join(coast_islands_dir, "LINHA_DE_COSTA_IMAGEM_GEOCOVER_SIRGAS_2000.shp")
islands_shp = os.path.join(coast_islands_dir, "ILHAS_IMAGEM_GEOCOVER_SIRGAS_2000.shp")

lat_min, lat_max = -18.5, -16.5
lon_min, lon_max = -40.0, -37.5

sst_pattern   = 'coraltemp_v3.1_*.nc'
dhw_pattern   = '*.nc'
kd490_pattern = 'AQUA_MODIS.*.L3m.DAY.KD.Kd_490.4km.nc'
par_pattern   = 'AQUA_MODIS.*.L3m.DAY.PAR.par.4km.nc'
chl_pattern   = 'AQUA_MODIS.*.L3m.DAY.CHL.chlor_a.4km.nc'

KDPAR_GATTUSO_A, KDPAR_GATTUSO_B, KDPAR_GATTUSO_C = 0.0665, 0.874, 0.00121
KDPAR_KD490_MIN_THRESHOLD = 0.001
MIN_DEPTH_METERS, MAX_DEPTH_METERS = 1, 200

# ==============================
# Carregamento da lista de sítios
# ==============================
try:
    df_sites_full = pd.read_csv(sites_csv_file, sep=';')
    df_sites_full.columns = df_sites_full.columns.str.strip()
    for col in df_sites_full.select_dtypes(['object']).columns:
        df_sites_full[col] = df_sites_full[col].str.strip()

    print("Agregando dados por SÍTIO, removendo distinção de habitat/profundidade.")
    df_sites_info = df_sites_full.drop_duplicates(subset=['Site_name'], keep='first').reset_index(drop=True)
    
    sites = list(zip(df_sites_info['Latitude'], df_sites_info['Longitude']))
    reef_names = list(df_sites_info['Reef_name'])
    site_names = list(df_sites_info['Site_name'])
    
    print(f"DataFrame original com {len(df_sites_full)} linhas.")
    print(f"DataFrame agregado com {len(df_sites_info)} sítios únicos.")
    print("DataFrame de sítios únicos carregado e processado com sucesso.")

except Exception as e:
    logging.critical(f"FALHA CRÍTICA AO CARREGAR O ARQUIVO DE SÍTIOS: {e}")
    raise

# ==============================
# Funções Auxiliares
# ==============================
def print_memory_usage(label=""): pass
def get_lat_slice(ds, lat_min, lat_max):
    lat_vals = ds['lat'].values
    return slice(lat_min, lat_max) if lat_vals[0] < lat_vals[-1] else slice(lat_max, lat_min)
def verify_file_integrity(filepath):
    checks = {'exists': os.path.exists(filepath), 'not_empty': False, 'readable': False, 'valid_netcdf': False}
    if checks['exists']:
        try:
            checks['not_empty'] = os.path.getsize(filepath) > 0
            with open(filepath, 'rb') as f: _ = f.read(100)
            checks['readable'] = True
            with xr.open_dataset(filepath) as ds: _ = ds.attrs
            checks['valid_netcdf'] = True
        except Exception as e:
            logging.warning(f"Erro na verificação de {filepath}: {str(e)}")
    return checks
def load_satellite_data(pattern, base_dir, period, var_name_options, chunks={'time': 30, 'lat': 100, 'lon': 100}):
    full_search_path = os.path.join(base_dir, pattern)
    print(f"\nBuscando em: {full_search_path} para o período {period}")
    all_files = sorted(glob.glob(full_search_path))
    if not all_files:
        logging.error(f"Nenhum arquivo encontrado para '{pattern}' em '{base_dir}'")
        return None
    candidate_paths = []
    date_pattern = re.compile(r'(\d{8})')
    for f_path in all_files:
        match = date_pattern.search(os.path.basename(f_path))
        if match and period[0] <= int(match.group(1)[:4]) <= period[1]:
            candidate_paths.append(f_path)
    if not candidate_paths:
        logging.warning(f"Nenhum arquivo encontrado para '{pattern}' no período {period}")
        return None
    valid_paths = [f for f in candidate_paths if verify_file_integrity(f)['valid_netcdf']]
    if not valid_paths:
        logging.error(f"Nenhum arquivo NetCDF válido encontrado para '{pattern}'")
        return None
    print(f"Carregando {len(valid_paths)} arquivos válidos...")
    def preprocess_with_time(ds):
        var_name = next((v for v in var_name_options if v in ds.data_vars), None)
        if var_name is None: return xr.Dataset()
        ds_subset = ds[[var_name]].astype('float32')
        if 'time' not in ds_subset.coords:
            filename = os.path.basename(ds.encoding.get("source", ""))
            match = re.search(r'(\d{8})', filename)
            if match:
                dt = pd.to_datetime(match.group(1), format='%Y%m%d')
                return ds_subset.expand_dims(time=[dt])
        return ds_subset
    try:
        ds_raw = xr.open_mfdataset(valid_paths, preprocess=preprocess_with_time, combine='by_coords', parallel=True, chunks=chunks, engine='netcdf4')
        lat_slice = get_lat_slice(ds_raw, lat_min, lat_max)
        ds_final = ds_raw.sel(lat=lat_slice, lon=slice(lon_min, lon_max))
        return ds_final.sortby('time')
    except Exception as e:
        logging.critical(f"Falha ao concatenar arquivos para '{pattern}': {e}")
        return None
def calculate_fixed_window_cv_optimized(data, window_size):
    epsilon = 1e-9
    min_periods = max(2, int(window_size * 0.25))
    with dask.config.set(scheduler='single-threaded'):
        rolling_mean = data.rolling(time=window_size, min_periods=min_periods, center=True).mean()
        rolling_std = data.rolling(time=window_size, min_periods=min_periods, center=True).std()
        return ((rolling_std / (rolling_mean + epsilon)) * 100).compute()
def calculate_kdpar_gattuso(kd490_da):
    kd490_safe = xr.where((kd490_da.isnull()) | (kd490_da <= KDPAR_KD490_MIN_THRESHOLD), np.nan, kd490_da)
    epsilon = 1e-9
    kdpar = KDPAR_GATTUSO_A + KDPAR_GATTUSO_B * kd490_safe - KDPAR_GATTUSO_C / (kd490_safe + epsilon)
    return xr.where((kdpar.isnull()) | (kdpar <= 0), np.nan, kdpar).rename("kdpar")
def calculate_benthic_dli(dli_surface_da, kdpar_da, depth_z_da):
    depth_z_positive = xr.where(depth_z_da < MIN_DEPTH_METERS, MIN_DEPTH_METERS, depth_z_da)
    attenuation = np.exp(-kdpar_da * depth_z_positive)
    return (dli_surface_da * attenuation).rename("benthic_dli")
def load_and_calculate_dli_series(kd490_pattern, par_pattern, period, base_dir, reference_da, depth_da):
    TEMP_DIR = os.path.join(output_dir, "temp_dli_files")
    os.makedirs(TEMP_DIR, exist_ok=True)
    CHUNKS = {'lat': 100, 'lon': 100}
    def get_date_map(files, period):
        date_map = {}
        for f in files:
            try:
                match = re.search(r'(\d{8})', os.path.basename(f))
                if match and period[0] <= int(match.group(1)[:4]) <= period[1] and verify_file_integrity(f)['valid_netcdf']:
                    date_map[pd.to_datetime(match.group(1), format='%Y%m%d')] = f
            except Exception as e:
                logging.error(f"Erro ao processar {f}: {str(e)}")
        return date_map
    kd490_map = get_date_map(sorted(glob.glob(os.path.join(base_dir, kd490_pattern))), period)
    par_map = get_date_map(sorted(glob.glob(os.path.join(base_dir, par_pattern))), period)
    reference_grid = reference_da.isel(time=0, drop=True).chunk(CHUNKS).load()
    depth_aligned = depth_da.interp_like(reference_grid, method='nearest').chunk(CHUNKS).load()
    common_dates = sorted(set(kd490_map.keys()) & set(par_map.keys()))
    if not common_dates: raise ValueError("Nenhuma data comum entre Kd490 e PAR encontrada")
    processed_files = []
    for i, date in enumerate(common_dates):
        try:
            with xr.open_dataset(kd490_map[date], chunks=CHUNKS) as ds_kd, \
                 xr.open_dataset(par_map[date], chunks=CHUNKS) as ds_par:
                kd490_interp = ds_kd['Kd_490'].interp_like(reference_grid, method='nearest')
                par_interp = ds_par['par'].interp_like(reference_grid, method='nearest')
                benthic_dli = calculate_benthic_dli(par_interp, calculate_kdpar_gattuso(kd490_interp), depth_aligned).expand_dims(time=[date])
                temp_file = os.path.join(TEMP_DIR, f"dli_temp_{date.strftime('%Y%m%d')}.nc")
                benthic_dli.to_netcdf(temp_file)
                processed_files.append(temp_file)
        except Exception as e:
            logging.error(f"Falha no processamento para {date}: {str(e)}")
            continue
    if not processed_files: raise ValueError("Nenhum arquivo DLI válido foi gerado")
    dli_final = xr.concat([xr.open_dataset(f)['benthic_dli'] for f in processed_files], dim='time').sortby('time')
    return dli_final.chunk({'time': 120, 'lat': 200, 'lon': 200})
def load_bathymetry(bathymetry_file, lat_min, lat_max, lon_min, lon_max):
    bathy_ds = xr.open_dataset(bathymetry_file)
    lat_slice = get_lat_slice(bathy_ds, lat_min, lat_max)
    bathy = bathy_ds['elevation'].sel(lat=lat_slice, lon=slice(lon_min, lon_max))
    mask_shallow = (((bathy >= -MAX_DEPTH_METERS) & (bathy < 0))).astype(float)
    return bathy, mask_shallow
def extract_site_values(da, sites):
    if da is None: return [np.nan] * len(sites)
    if isinstance(da, xr.Dataset): da = da[list(da.data_vars)[0]]
    lats = xr.DataArray([s[0] for s in sites], dims="site"); lons = xr.DataArray([s[1] for s in sites], dims="site")
    return da.sel(lat=lats, lon=lons, method='nearest').values
def load_coast_and_islands(coast_file, islands_file, lat_min, lat_max, lon_min, lon_max):
    bounding_box = box(lon_min, lat_min, lon_max, lat_max)
    bbox_gdf = gpd.GeoDataFrame({'geometry': [bounding_box]}, crs="EPSG:4326")
    try:
        gdf_coast_line = gpd.read_file(coast_file).to_crs("EPSG:4326")
        gdf_islands = gpd.read_file(islands_file).to_crs("EPSG:4326")
        gdf_islands_clipped = gpd.overlay(gdf_islands, bbox_gdf, how='intersection')
        coast_union = gdf_coast_line.geometry.unary_union
        bbox_line = LineString(list(bounding_box.exterior.coords))
        combined_lines = unary_union([coast_union, bbox_line])
        polygons = list(polygonize(combined_lines))
        ocean_polygon = max([p for p in polygons if p.intersects(bounding_box)], key=lambda p: p.area, default=bounding_box)
        land_polygon = bounding_box.difference(ocean_polygon)
        if not gdf_islands_clipped.empty: land_polygon = unary_union([land_polygon, gdf_islands_clipped.unary_union])
        return gpd.GeoDataFrame({'geometry': [land_polygon]}, crs="EPSG:4326")
    except Exception as e:
        logging.warning(f"Não foi possível ler shapefiles de costa/ilhas: {e}"); return None

# ==============================
# FUNÇÃO DE PLOTAGEM COMPOSTA (ALTA QUALIDADE)
# ==============================
def create_publication_composite(data_dict, output_path, mask_shallow_data, bathy_data):
    """Gera uma figura composta de alta qualidade para publicação."""
    panels_config = [
        ('mean_SST', 'Mean SST', 'viridis', '°C'),
        ('cv_SST', 'SST Variability (CV)', 'cividis', '%'),
        ('prop_DHW_gt4', 'Thermal Stress (DHW > 4)', 'YlOrRd', '%'),
        ('mean_DLI', 'Mean Benthic DLI', 'magma', r'mol m$^{-2}$ d$^{-1}$'),
        ('cv_DLI', 'DLI Variability (CV)', 'inferno', '%'),
        ('mean_CHL', 'Mean Chlorophyll-a', 'YlGn', r'mg m$^{-3}$'),
    ]

    fig = plt.figure(figsize=FIGSIZE, constrained_layout=False)
    gs = gridspec.GridSpec(2, 3, figure=fig, wspace=0.25, hspace=0.35)

    gdf_land = load_coast_and_islands(coast_shp, islands_shp, lat_min, lat_max, lon_min, lon_max)
    
    for i, (key, title, cmap_name, unit) in enumerate(panels_config):
        row, col = divmod(i, 3)
        ax = fig.add_subplot(gs[row, col])
        
        da = data_dict.get(key)
        if da is None:
            ax.text(0.5, 0.5, 'Data unavailable', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(title)
            continue

        if mask_shallow_data is not None:
            mask_aligned = mask_shallow_data.interp_like(da, method='nearest')
            da = da.where(mask_aligned >= 0.5)

        data_to_plot = da.load()
        valid_data = data_to_plot.values[np.isfinite(data_to_plot.values)]
        vmin, vmax = (np.nanpercentile(valid_data, 2), np.nanpercentile(valid_data, 98)) if valid_data.size > 0 else (0, 1)

        cmap_obj = plt.get_cmap(cmap_name).copy()
        cmap_obj.set_bad(color='lightgray')

        img = ax.pcolormesh(data_to_plot.lon, data_to_plot.lat, data_to_plot,
                            cmap=cmap_obj, shading='auto', vmin=vmin, vmax=vmax)

        if gdf_land is not None:
            gdf_land.plot(ax=ax, facecolor='#D2B48C', edgecolor='black', linewidth=0.5, zorder=3)
        
        if bathy_data is not None:
            ax.contour(bathy_data.lon, bathy_data.lat, bathy_data, levels=[-200, -50], 
                       colors=['dimgray', 'gray'], linewidths=[0.7, 0.5], linestyles=['--', ':'], zorder=2)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.08)
        cbar = fig.colorbar(img, cax=cax, orientation='vertical')
        cbar.set_label(unit, fontsize=8)
        cbar.ax.tick_params(labelsize=7)

        ax.text(-0.15, 1.1, f'({PANEL_LABELS[i]})', transform=ax.transAxes, fontsize=11, fontweight='bold', va='top')
        
        ax.set_xlim(lon_min, lon_max)
        ax.set_ylim(lat_min, lat_max)
        ax.set_aspect('equal', adjustable='box')
        ax.set_title(title, fontsize=10, pad=10)
        
        if col == 0:
            ax.set_ylabel('Latitude (°S)', fontsize=9)
        else:
            ax.set_yticklabels([])
            
        if row == 1:
            ax.set_xlabel('Longitude (°W)', fontsize=9)
        else:
            ax.set_xticklabels([])

        if i == 0:
            scalebar = AnchoredSizeBar(ax.transData, 0.5, '50 km', 'lower left', 
                                       pad=0.3, color='black', frameon=False,
                                       size_vertical=0.02,
                                       fontproperties=fm.FontProperties(size=7))
            ax.add_artist(scalebar)

    plt.savefig(output_path, dpi=600, bbox_inches='tight', facecolor='white')
    if output_path.endswith('.png'):
        plt.savefig(output_path.replace('.png', '.pdf'), format='pdf', bbox_inches='tight')
    plt.close(fig)
    print(f"Figura composta salva em: {output_path}")

def plot_spatial_map(da, title, output_filename, cmap='viridis', cbar_label=''):
    if da is None: return
    # Preservado caso o usuário ainda queira mapas individuais no futuro
    if isinstance(da, xr.Dataset): da = da[list(da.data_vars)[0]]
    if mask_shallow is not None:
        mask_aligned = mask_shallow.interp_like(da, method='nearest')
        da = da.where(mask_aligned >= 0.5)
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_xlim(lon_min, lon_max); ax.set_ylim(lat_min, lat_max)
    cmap_obj = plt.get_cmap(cmap).copy(); cmap_obj.set_bad('white', 1.)
    img = ax.pcolormesh(da.lon, da.lat, da, cmap=cmap_obj, shading='auto')
    plt.colorbar(img, ax=ax, label=cbar_label)
    ax.set_title(title)
    plt.tight_layout(); plt.savefig(output_filename, dpi=300); plt.close(fig)

def scale_marker_sizes(values, scale_factor=300, min_size=30):
    vals = np.asarray(values, dtype=float); valid = np.isfinite(vals)
    if valid.sum() == 0: return np.full(vals.shape, min_size)
    vmin, vmax = vals[valid].min(), vals[valid].max()
    if vmax == vmin: return np.full(vals.shape, min_size + scale_factor/2)
    scaled = (vals - vmin) / (vmax - vmin); scaled[~valid] = 0
    return scaled * scale_factor + min_size

# ==============================
# Processamento Principal
# ==============================
print("\n--- 1. Carregando e Preparando Dados Base ---")
bathy_data, mask_shallow = load_bathymetry(bathymetry_file, lat_min, lat_max, lon_min, lon_max)
sst_all_ds = load_satellite_data(sst_pattern, sst_dir, sst_period, ['analysed_sst', 'sea_surface_temperature', 'sst'])
dhw_all_ds = load_satellite_data(dhw_pattern, dhw_dir, dhw_period, ['degree_heating_week'])
sst_all = sst_all_ds[list(sst_all_ds.data_vars)[0]].astype('float32').chunk({'time': 30, 'lat': 100, 'lon': 100})
dhw_all = dhw_all_ds[list(dhw_all_ds.data_vars)[0]].astype('float32').chunk({'time': 'auto'}) if dhw_all_ds is not None else None
del sst_all_ds, dhw_all_ds
gc.collect()

print("\n--- 2. Carregando e Preparando Dados Auxiliares ---")
chl_all_ds = load_satellite_data(chl_pattern, chl_dir, chl_period, ['chlor_a'])
if chl_all_ds is not None:
    chl_all = chl_all_ds.interp_like(sst_all, method='nearest')[list(chl_all_ds.data_vars)[0]].astype('float32').chunk({'time': 'auto'})
    del chl_all_ds
else:
    chl_all = None
gc.collect()

print("\nCalculando mapa de DLI com batimetria GEBCO...")
depth_z_meters = -bathy_data.where((bathy_data >= -MAX_DEPTH_METERS) & (bathy_data <= -MIN_DEPTH_METERS))
dli_all = load_and_calculate_dli_series(kd490_pattern, par_pattern, light_period, modis_dir, sst_all, depth_z_meters)

print("\n--- 3. Definindo e Computando Métricas Espaciais ---")
computed_results = {}
computed_results['mean_SST'] = sst_all.mean('time', skipna=True).compute()

if analysis_selector == 'all':
    computed_results['cv_all'] = ((sst_all.std('time', skipna=True) / (sst_all.mean('time', skipna=True) + 1e-9)) * 100).compute()
else:
    window = int(analysis_selector)
    cv = calculate_fixed_window_cv_optimized(sst_all, window)
    computed_results[f'cv_{window}'] = cv.mean('time', skipna=True).compute()

if dhw_all is not None:
    computed_results['prop_DHW_gt4'] = (dhw_all > 4).mean('time', skipna=True).compute() * 100
    computed_results['prop_DHW_gt8'] = (dhw_all > 8).mean('time', skipna=True).compute() * 100
if chl_all is not None:
    computed_results['mean_CHL'] = chl_all.mean('time', skipna=True).compute()
if dli_all is not None:
    computed_results['mean_DLI_map'] = dli_all.mean('time', skipna=True).compute()
    if analysis_selector == 'all':
        computed_results['dli_cv_all'] = ((dli_all.std('time', skipna=True) / (dli_all.mean('time', skipna=True) + 1e-9)) * 100).compute()
    else:
        window = int(analysis_selector)
        cv_dli = calculate_fixed_window_cv_optimized(dli_all, window)
        computed_results[f'dli_cv_{window}'] = cv_dli.mean('time', skipna=True).compute()

# --- Coleta de Dados para Figura Composta ---
if analysis_selector == 'all':
    sel_cv_sst, sel_cv_dli = 'cv_all', 'dli_cv_all'
else:
    sel_cv_sst, sel_cv_dli = f'cv_{int(analysis_selector)}', f'dli_cv_{int(analysis_selector)}'

data_for_composite = {
    'mean_SST': computed_results.get('mean_SST'),
    'cv_SST': computed_results.get(sel_cv_sst),
    'prop_DHW_gt4': computed_results.get('prop_DHW_gt4'),
    'mean_DLI': computed_results.get('mean_DLI_map'),
    'cv_DLI': computed_results.get(sel_cv_dli),
    'mean_CHL': computed_results.get('mean_CHL'),
}

composite_output_path = os.path.join(output_dir, "Fig_Main_Environmental_Composite.png")
create_publication_composite(data_for_composite, composite_output_path, mask_shallow, bathy_data)

# Mantendo geração de planilhas para integridade da análise
global_vars_to_extract = [sel_cv_sst, sel_cv_dli, 'prop_DHW_gt4', 'mean_SST', 'mean_CHL', 'mean_DLI_map']
df_sites_base_global = df_sites_info.copy()
for var in global_vars_to_extract:
    df_sites_base_global[var] = extract_site_values(computed_results.get(var), sites)
df_sites_base_global.to_excel(os.path.join(output_dir, "planilha_base_analise_global_sites.xlsx"), index=False)

# --- PCA GLOBAL (Opcional - Mantido para compatibilidade) ---
# ... (Código da PCA omitido para brevidade, focado na figura composta solicitada) ...

print("\nProcessamento da Figura Composta de Alta Qualidade finalizado!")
