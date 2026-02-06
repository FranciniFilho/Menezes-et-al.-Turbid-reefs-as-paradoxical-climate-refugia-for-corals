### PCA GLOBAL COMPOSITE GEMINI - Publication Quality Maps (Nature/Science)
### Magnitude and Variability Maps with Sliding Window CV
### Versão com duas figuras separadas e janelas temporais integradas

# ==============================
# Imports
# ==============================
import os
import sys
print(f"Python Executable: {sys.executable}")
print(f"Python Version: {sys.version}")
# print(f"Python Path: {sys.path}")
try:
    import geopandas as gpd
    print(f"Geopandas version: {gpd.__version__}")
except ImportError as e:
    print(f"FAILED TO IMPORT GEOPANDAS: {e}")
    # try to list what's in site-packages
    import site
    print(f"Site packages: {site.getsitepackages()}")

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
from dask.distributed import Client
import dask
import psutil
import gc
import pickle

# --- Imports adicionais para figuras de alta qualidade ---
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar
import matplotlib.font_manager as fm
import matplotlib.patches as mpatches

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

# ==============================
# Configurações Iniciais e Caminhos
# ==============================
logging.basicConfig(
    filename='file_open_errors_global_pca_composite_gemini.log',
    filemode='w',
    level=logging.ERROR,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# NOVO DIRETÓRIO DE OUTPUS (DEDICADO EM FINAL_RESULTS)
output_dir = r"C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_RESULTS\Global_maps_GEMINI_TEST"
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

# ATUALIZADO PARA DRIVE K:
sst_dir   = r"K:\remote sensing\CRW_SST_FULL"
modis_dir = r"K:\remote sensing\MODIS_DATA_FULL"
chl_dir   = r"K:\remote sensing\MODIS_DATA_FULL"
dhw_dir   = r"K:\remote sensing\CRW_DHW_FULL"

sst_period   = (2008, 2008)
dhw_period   = (2008, 2008)
light_period = (2008, 2008)
chl_period   = (2008, 2008)

bathymetry_file      = r"c:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\GIS shapes\gebco_2024_ASO.nc"
protected_areas_shp  = r"c:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\GIS shapes\uc_fed_agosto_2016_site_shp\uc_fed_agosto_2016_site_shp\uc_fed_agosto_2016_site.shp"
coast_islands_dir    = r"c:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\GIS shapes\batimetria_new"
reefs_shp            = r"c:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\GIS shapes\Recifes_Banco_dos_Abrolhos\Recifes_Banco_dos_Abrolhos.shp"
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
# Funções Auxiliares
# ==============================
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

def calculate_cv_map(data_da, window_size):
    """
    Calcula um mapa espacial de CV usando janela deslizante temporal.
    
    Para window_size == 'all':
        CV = (std temporal / média temporal) * 100
    Para window_size numérico (ex: 2, 30):
        1. Calcula CV(t) para cada janela centrada de tamanho `window_size`.
        2. Retorna a MÉDIA TEMPORAL de CV(t).
    """
    epsilon = 1e-9
    
    if window_size == 'all':
        mean_map = data_da.mean('time', skipna=True)
        std_map = data_da.std('time', skipna=True)
        cv_map = (std_map / (mean_map + epsilon)) * 100
        return cv_map.compute()
    else:
        window = int(window_size)
        min_periods = max(2, int(window * 0.25))
        
        with dask.config.set(scheduler='single-threaded'):
            rolling_mean = data_da.rolling(time=window, min_periods=min_periods, center=True).mean()
            rolling_std = data_da.rolling(time=window, min_periods=min_periods, center=True).std()
            cv_timeseries = (rolling_std / (rolling_mean + epsilon)) * 100
            
            cv_map = cv_timeseries.mean('time', skipna=True)
            return cv_map.compute()

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
# FUNÇÕES DE PLOTAGEM (MAGNITUDE E VARIABILIDADE)
# ==============================
def create_magnitude_map(data_dict, output_path, mask_shallow_data, bathy_data):
    """Gera figura 2x2 com mapas de magnitude (médias históricas)."""
    panels_config = [
        ('mean_SST', 'Mean SST', 'viridis', '°C'),
        ('mean_DLI', 'Mean Benthic DLI', 'magma', r'mol m$^{-2}$ d$^{-1}$'),
        ('mean_CHL', 'Mean Chlorophyll-a', 'YlGn', r'mg m$^{-3}$'),
        ('prop_DHW_gt4', 'Thermal Stress (DHW > 4)', 'YlOrRd', '% of time'),
    ]
    
    PANEL_LABELS_MAG = ['a', 'b', 'c', 'd']
    FIGSIZE_MAG = (180 / 25.4, 160 / 25.4)
    
    fig = plt.figure(figsize=FIGSIZE_MAG, constrained_layout=False)
    gs = gridspec.GridSpec(2, 2, figure=fig, wspace=0.25, hspace=0.30)
    
    gdf_land = load_coast_and_islands(coast_shp, islands_shp, lat_min, lat_max, lon_min, lon_max)
    
    for i, (key, title, cmap_name, unit) in enumerate(panels_config):
        row, col = divmod(i, 2)
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
            gdf_land.plot(ax=ax, facecolor='#D2B48C', edgecolor='black', linewidth=0.5, zorder=10)
        
        # === PROTECTED AREAS (marine portion only) ===
        try:
            if os.path.exists(protected_areas_shp) and gdf_land is not None:
                gdf_protected = gpd.read_file(protected_areas_shp).to_crs("EPSG:4326")
                # Clip to bounding box first
                bbox = box(lon_min, lat_min, lon_max, lat_max)
                bbox_gdf = gpd.GeoDataFrame({'geometry': [bbox]}, crs="EPSG:4326")
                gdf_protected_clipped = gpd.overlay(gdf_protected, bbox_gdf, how='intersection')
                # Extract only marine portion (remove land overlap)
                protected_marine = gpd.overlay(gdf_protected_clipped, gdf_land, how='difference')
                if not protected_marine.empty:
                    protected_marine.boundary.plot(ax=ax, edgecolor='darkblue', linewidth=1.2, zorder=5)
        except Exception as e:
            logging.warning(f"Não foi possível plotar áreas protegidas: {e}")
        
        # === REEFS ===
        try:
            if os.path.exists(reefs_shp):
                gdf_reefs = gpd.read_file(reefs_shp).to_crs("EPSG:4326")
                gdf_reefs.boundary.plot(ax=ax, edgecolor='purple', linewidth=0.8, zorder=15)
        except Exception as e:
            logging.warning(f"Não foi possível plotar recifes: {e}")
        
        if bathy_data is not None:
            ax.contour(bathy_data.lon, bathy_data.lat, bathy_data, levels=[-200, -50],
                       colors=['dimgray', 'gray'], linewidths=[0.7, 0.5], linestyles=['--', ':'], zorder=2)
        
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.08)
        cbar = fig.colorbar(img, cax=cax, orientation='vertical')
        cbar.set_label(unit, fontsize=8)
        cbar.ax.tick_params(labelsize=7)
        
        ax.text(-0.15, 1.15, f'({PANEL_LABELS_MAG[i]})', transform=ax.transAxes, fontsize=11, fontweight='bold', va='top')
        
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
    print(f"Mapa de Magnitude salvo em: {output_path}")

def create_variability_map(cv_data_dict, output_path, mask_shallow_data, bathy_data):
    """Gera figura 3x3 com mapas de variabilidade (CV) para 3 variáveis x 3 janelas."""
    variables = ['SST', 'DLI', 'CHL']
    windows = [2, 30, 'all']
    
    col_titles = ['CV (2-day window)', 'CV (30-day window)', 'CV (full series)']
    row_titles = ['SST', 'Benthic DLI', 'Chlorophyll-a']
    
    PANEL_LABELS_VAR = [chr(ord('a') + i) for i in range(9)]
    FIGSIZE_VAR = (250 / 25.4, 220 / 25.4)
    
    fig = plt.figure(figsize=FIGSIZE_VAR, constrained_layout=False)
    gs = gridspec.GridSpec(3, 3, figure=fig, wspace=0.20, hspace=0.25)
    
    gdf_land = load_coast_and_islands(coast_shp, islands_shp, lat_min, lat_max, lon_min, lon_max)
    
    panel_idx = 0
    for row_idx, var in enumerate(variables):
        for col_idx, window in enumerate(windows):
            ax = fig.add_subplot(gs[row_idx, col_idx])
            
            key = f'cv_{var}_{window}'
            da = cv_data_dict.get(key)
            
            if da is None:
                ax.text(0.5, 0.5, 'Data unavailable', ha='center', va='center', transform=ax.transAxes)
            else:
                if mask_shallow_data is not None:
                    mask_aligned = mask_shallow_data.interp_like(da, method='nearest')
                    da = da.where(mask_aligned >= 0.5)
                
                data_to_plot = da.load()
                valid_data = data_to_plot.values[np.isfinite(data_to_plot.values)]
                vmin, vmax = (np.nanpercentile(valid_data, 2), np.nanpercentile(valid_data, 98)) if valid_data.size > 0 else (0, 100)
                
                cmap_obj = plt.get_cmap('cividis').copy()
                cmap_obj.set_bad(color='lightgray')
                
                img = ax.pcolormesh(data_to_plot.lon, data_to_plot.lat, data_to_plot,
                                    cmap=cmap_obj, shading='auto', vmin=vmin, vmax=vmax)
                
                if gdf_land is not None:
                    gdf_land.plot(ax=ax, facecolor='#D2B48C', edgecolor='black', linewidth=0.5, zorder=10)
                
                # === PROTECTED AREAS (marine portion only) ===
                try:
                    if os.path.exists(protected_areas_shp) and gdf_land is not None:
                        gdf_protected = gpd.read_file(protected_areas_shp).to_crs("EPSG:4326")
                        bbox = box(lon_min, lat_min, lon_max, lat_max)
                        bbox_gdf = gpd.GeoDataFrame({'geometry': [bbox]}, crs="EPSG:4326")
                        gdf_protected_clipped = gpd.overlay(gdf_protected, bbox_gdf, how='intersection')
                        protected_marine = gpd.overlay(gdf_protected_clipped, gdf_land, how='difference')
                        if not protected_marine.empty:
                            protected_marine.boundary.plot(ax=ax, edgecolor='darkblue', linewidth=1.0, zorder=5)
                except Exception as e:
                    logging.warning(f"Não foi possível plotar áreas protegidas: {e}")
                
                # === REEFS ===
                try:
                    if os.path.exists(reefs_shp):
                        gdf_reefs = gpd.read_file(reefs_shp).to_crs("EPSG:4326")
                        gdf_reefs.boundary.plot(ax=ax, edgecolor='purple', linewidth=0.6, zorder=15)
                except Exception as e:
                    logging.warning(f"Não foi possível plotar recifes: {e}")
                
                if bathy_data is not None:
                    ax.contour(bathy_data.lon, bathy_data.lat, bathy_data, levels=[-200, -50],
                               colors=['dimgray', 'gray'], linewidths=[0.7, 0.5], linestyles=['--', ':'], zorder=2)
                
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cbar = fig.colorbar(img, cax=cax, orientation='vertical')
                cbar.set_label('CV (%)', fontsize=7)
                cbar.ax.tick_params(labelsize=6)
            
            ax.text(-0.12, 1.14, f'({PANEL_LABELS_VAR[panel_idx]})', transform=ax.transAxes, 
                    fontsize=10, fontweight='bold', va='top')
            
            if row_idx == 0:
                ax.set_title(col_titles[col_idx], fontsize=9, pad=8)
            
            if col_idx == 0:
                ax.text(-0.25, 0.5, row_titles[row_idx], transform=ax.transAxes, 
                        fontsize=10, fontweight='bold', va='center', ha='right', rotation=90)
            
            ax.set_xlim(lon_min, lon_max)
            ax.set_ylim(lat_min, lat_max)
            ax.set_aspect('equal', adjustable='box')
            
            if col_idx == 0:
                ax.set_ylabel('Lat (°S)', fontsize=8)
            else:
                ax.set_yticklabels([])
            if row_idx == 2:
                ax.set_xlabel('Lon (°W)', fontsize=8)
            else:
                ax.set_xticklabels([])
            
            ax.tick_params(axis='both', labelsize=6)
            
            if panel_idx == 0:
                scalebar = AnchoredSizeBar(ax.transData, 0.5, '50 km', 'lower left',
                                           pad=0.2, color='black', frameon=False,
                                           size_vertical=0.015,
                                           fontproperties=fm.FontProperties(size=6))
                ax.add_artist(scalebar)
            
            panel_idx += 1
    
    plt.savefig(output_path, dpi=600, bbox_inches='tight', facecolor='white')
    if output_path.endswith('.png'):
        plt.savefig(output_path.replace('.png', '.pdf'), format='pdf', bbox_inches='tight')
    plt.close(fig)
    print(f"Mapa de Variabilidade salvo em: {output_path}")

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

print("\n--- 3. Calculando Métricas Espaciais ---")
computed_results = {}

# === CACHE CONFIGURATION ===
CACHE_DIR = os.path.join(output_dir, "metrics_cache")
os.makedirs(CACHE_DIR, exist_ok=True)

def load_cached_metric(key):
    cache_file = os.path.join(CACHE_DIR, f"{key}.pkl")
    if os.path.exists(cache_file):
        try:
            with open(cache_file, 'rb') as f:
                print(f"      [CACHE HIT] Carregando {key} do cache...")
                return pickle.load(f)
        except Exception as e:
            print(f"      [CACHE ERROR] Falha ao carregar {key}: {e}")
    return None

def save_cached_metric(key, data):
    cache_file = os.path.join(CACHE_DIR, f"{key}.pkl")
    try:
        with open(cache_file, 'wb') as f:
            pickle.dump(data, f)
        print(f"      [CACHE SAVE] {key} salvo no cache.")
    except Exception as e:
        print(f"      [CACHE ERROR] Falha ao salvar {key}: {e}")

def get_or_compute_metric(key, compute_func):
    cached = load_cached_metric(key)
    if cached is not None:
        return cached
    print(f"      Computando {key}...")
    result = compute_func()
    save_cached_metric(key, result)
    return result


# === MAGNITUDE ===
print("  Calculando métricas de MAGNITUDE...")
computed_results['mean_SST'] = get_or_compute_metric('mean_SST', lambda: sst_all.mean('time', skipna=True).compute())
if chl_all is not None:
    computed_results['mean_CHL'] = get_or_compute_metric('mean_CHL', lambda: chl_all.mean('time', skipna=True).compute())
if dli_all is not None:
    computed_results['mean_DLI'] = get_or_compute_metric('mean_DLI', lambda: dli_all.mean('time', skipna=True).compute())
if dhw_all is not None:
    computed_results['prop_DHW_gt4'] = get_or_compute_metric('prop_DHW_gt4', lambda: (dhw_all > 4).mean('time', skipna=True).compute() * 100)

# === VARIABILIDADE (3 janelas: 2, 30, all) ===
print("  Calculando métricas de VARIABILIDADE...")
windows_to_compute = [2, 30, 'all']

# SST CV
print("    - SST CV para todas as janelas...")
for window in windows_to_compute:
    key = f'cv_SST_{window}'
    computed_results[key] = get_or_compute_metric(key, lambda w=window: calculate_cv_map(sst_all, w))
    gc.collect()

# DLI CV
if dli_all is not None:
    print("    - DLI CV para todas as janelas...")
    for window in windows_to_compute:
        key = f'cv_DLI_{window}'
        computed_results[key] = get_or_compute_metric(key, lambda w=window: calculate_cv_map(dli_all, w))
        gc.collect()

# CHL CV
if chl_all is not None:
    print("    - CHL CV para todas as janelas...")
    for window in windows_to_compute:
        key = f'cv_CHL_{window}'
        computed_results[key] = get_or_compute_metric(key, lambda w=window: calculate_cv_map(chl_all, w))
        gc.collect()

print("  Cálculo de métricas concluído!")

# === GERAÇÃO DAS FIGURAS ===
print("\n--- 4. Gerando Figuras ---")

magnitude_data = {
    'mean_SST': computed_results.get('mean_SST'),
    'mean_DLI': computed_results.get('mean_DLI'),
    'mean_CHL': computed_results.get('mean_CHL'),
    'prop_DHW_gt4': computed_results.get('prop_DHW_gt4'),
}

variability_data = {}
for var in ['SST', 'DLI', 'CHL']:
    for window in [2, 30, 'all']:
        key = f'cv_{var}_{window}'
        variability_data[key] = computed_results.get(key)

magnitude_output_path = os.path.join(output_dir, "Fig_Magnitude.png")
create_magnitude_map(magnitude_data, magnitude_output_path, mask_shallow, bathy_data)

variability_output_path = os.path.join(output_dir, "Fig_Variability.png")
create_variability_map(variability_data, variability_output_path, mask_shallow, bathy_data)

print("\n=== PROCESSAMENTO FINALIZADO COM SUCESSO ===")
print(f"Figuras salvas em: {output_dir}")
