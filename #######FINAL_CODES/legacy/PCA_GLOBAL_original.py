### Frequencia média DHW > 4 e > 8 per site + Mapas (versão com escala de cores corrigida)
### Períodos entre anos (Julho a Junho)
###GIF high

import os
import glob
import re
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import geopandas as gpd
from shapely.geometry import box, LineString
from shapely.ops import unary_union, polygonize
import matplotlib.lines as mlines
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap  # Import necessário para criar colormap personalizado

# ============================================================
# Configurações iniciais
# ============================================================
dhw_dir = r'H:\remote sensing\CRW_DHW_FULL'

lat_min, lat_max = -20.5, -14.5
lon_min, lon_max = -40.5, -35.5

#lat_min, lat_max = -18.7, -16.5 
#lon_min, lon_max = -36.9, -36.6

bathymetry_file = r'C:\Users\rbfra\OneDrive\GIS shapes\gebco_2024_ASO.nc'
output_dir = (r'C:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\########NEW RESULTS\output_DHW_frequency_Jul-Jun_GIF') 
os.makedirs(output_dir, exist_ok=True)

protected_areas_shp = r'H:\remote sensing\Mascaras\uc_fed_agosto_2016_site_shp\uc_fed_agosto_2016_site.shp'
coast_shp = (r'C:\Users\rbfra\OneDrive\GIS shapes\batimetria_new'
             r'\LINHA_DE_COSTA_IMAGEM_GEOCOVER_SIRGAS_2000.shp')
islands_shp = (r'C:\Users\rbfra\OneDrive\GIS shapes\batimetria_new'
               r'\ILHAS_IMAGEM_GEOCOVER_SIRGAS_2000.shp')
reefs_shp = r'H:\remote sensing\Mascaras\Recifes_Banco_dos_Abrolhos\Recifes_Banco_dos_Abrolhos.shp'

# Anos permitidos (anos de início do período Jul-Jun)
allowed_years = set(range(1985, 2026))  # inclui de 1985 até 2008

# ============================================================
# Funções auxiliares
# ============================================================
def get_lat_slice(ds, lat_min, lat_max):
    lat_vals = ds['lat'].values
    return slice(lat_min, lat_max) if lat_vals[0] < lat_vals[-1] else slice(lat_max, lat_min)

def load_bathymetry(bathy_file, lat_min, lat_max, lon_min, lon_max):
    ds = xr.open_dataset(bathy_file)
    lat_slice = get_lat_slice(ds, lat_min, lat_max)
    bathy = ds['elevation'].sel(lon=slice(lon_min, lon_max), lat=lat_slice)
    mask = (bathy >= -200).astype(float)
    return bathy, mask

def adjust_longitudes(da):
    if (da.lon > 180).any():
        da = da.assign_coords(lon=(((da.lon + 180) % 360) - 180)).sortby('lon')
    return da

def extract_date_from_filename(fname):
    match = re.search(r'(\d{8})', os.path.basename(fname))
    return pd.to_datetime(match.group(1), format='%Y%m%d') if match else None

def load_shapefile(path, name="shp"):
    gdf = gpd.read_file(path)
    if gdf.crs != "EPSG:4326":
        gdf = gdf.to_crs("EPSG:4326")
    return gdf

def load_coast_and_islands(coast_file, islands_file, lat_min, lat_max, lon_min, lon_max):
    bbox = box(lon_min, lat_min, lon_max, lat_max)
    bbox_gdf = gpd.GeoDataFrame({'geometry': [bbox]}, crs="EPSG:4326")
    coast = gpd.read_file(coast_file).to_crs("EPSG:4326")
    islands = gpd.read_file(islands_file).to_crs("EPSG:4326")
    islands_clip = gpd.overlay(islands, bbox_gdf, how='intersection')
    coast_union = coast.geometry.unary_union
    bbox_line = LineString(list(bbox.exterior.coords))
    polygons = list(polygonize(unary_union([coast_union, bbox_line])))
    ocean_poly, max_area = None, 0
    for poly in polygons:
        inter = poly.intersection(bbox)
        if not inter.is_empty and inter.area > max_area:
            ocean_poly, max_area = inter, inter.area
    if ocean_poly is None:
        ocean_poly = bbox
    land_poly = bbox.difference(ocean_poly)
    if not islands_clip.empty:
        land_poly = unary_union([land_poly, islands_clip.unary_union])
    land = gpd.GeoDataFrame({'geometry': [land_poly]}, crs="EPSG:4326")
    return gpd.overlay(land, bbox_gdf, how='intersection')

# ============================================================
# Carregar dados estáticos (batimetria, shapefiles)
# ============================================================
bathy, marine_mask = load_bathymetry(bathymetry_file, lat_min, lat_max, lon_min, lon_max)
water_mask = (bathy < 0).astype(float)
marine_mask = ((bathy >= -200) & (bathy < 0)).astype(float)

gdf_protected = load_shapefile(protected_areas_shp, "Área Protegida")
gdf_coast = load_shapefile(coast_shp, "Linha de Costa")
gdf_islands = load_shapefile(islands_shp, "Ilhas")
gdf_recifes = load_shapefile(reefs_shp, "Recifes")
gdf_land = load_coast_and_islands(coast_shp, islands_shp, lat_min, lat_max, lon_min, lon_max)
gdf_islands_bound = gpd.read_file(islands_shp).to_crs("EPSG:4326")

# ============================================================
# Processamento dos arquivos DHW
# ============================================================
dhw_files = sorted(glob.glob(os.path.join(dhw_dir, '*.nc')))
print(f"Arquivos DHW encontrados: {len(dhw_files)}")

dhw_bool_4, dhw_bool_8 = [], []
for f in dhw_files:
    dt = extract_date_from_filename(f)
    # ### ALTERAÇÃO ###
    # A verificação de ano continua a mesma, pois estamos interessados nos dados desses anos,
    # a agregação será feita depois.
    if dt is None or dt.year not in range(min(allowed_years), max(allowed_years) + 2): # Carrega um ano a mais para completar o último período
        continue

    try:
        ds = xr.open_dataset(f)
    except Exception as e:
        print(f"[AVISO] Arquivo corrompido identificado: {f}. Erro: {e}.")
        print("Este arquivo será removido e não entrará nas análises.")
        os.remove(f)
        continue

    dhw_raw = ds['degree_heating_week']

    if 'time' in dhw_raw.dims:
        dhw_raw = dhw_raw.isel(time=0, drop=True)

    dhw = dhw_raw.sel(
        lat=get_lat_slice(ds, lat_min, lat_max),
        lon=slice(lon_min, lon_max)
    )
    dhw = adjust_longitudes(dhw)

    mask_interp = marine_mask.interp(lat=dhw.lat, lon=dhw.lon, method='nearest') >= 0.5
    dhw_masked = dhw.where(mask_interp)

    dhw_bool_4.append((dhw_masked > 4).astype(int).expand_dims(time=[dt]))
    dhw_bool_8.append((dhw_masked > 8).astype(int).expand_dims(time=[dt]))

# Empilhar e preencher faltantes
stack_4 = xr.concat(dhw_bool_4, dim='time').sortby('time').fillna(0)
stack_8 = xr.concat(dhw_bool_8, dim='time').sortby('time').fillna(0)

# ### ALTERAÇÃO ###
# Agrupar por "ano hidrológico" (Julho a Junho) usando resample.
# 'AS-JUL' significa "Annual Start, July"
print("Agrupando dados por períodos de Julho a Junho...")
hydro_yearly_sum_4 = stack_4.resample(time='AS-JUL').sum(dim='time')
hydro_yearly_sum_8 = stack_8.resample(time='AS-JUL').sum(dim='time')

# Contar o número de dias em cada período para um cálculo de percentual preciso
hydro_yearly_count = stack_4.resample(time='AS-JUL').count(dim='time')

hydro_yearly_percent_4 = (hydro_yearly_sum_4 / hydro_yearly_count) * 100
hydro_yearly_percent_8 = (hydro_yearly_sum_8 / hydro_yearly_count) * 100

# =========================================================================
# ### ETAPA DE CORREÇÃO: ENCONTRAR O VMAX GLOBAL PARA CONSISTÊNCIA ###
# =========================================================================
# Encontra o valor máximo absoluto de frequência em todos os anos para cada limiar.
global_vmax_4 = np.nanmax(hydro_yearly_percent_4.values)
global_vmax_8 = np.nanmax(hydro_yearly_percent_8.values)

# Arredonda o vmax para o próximo múltiplo de 5 para uma legenda mais limpa
vmax_plot_4 = np.ceil(global_vmax_4 / 5) * 5 if global_vmax_4 > 0 else 1
vmax_plot_8 = np.ceil(global_vmax_8 / 5) * 5 if global_vmax_8 > 0 else 1

print("\n--- Análise da Escala de Cores Global ---")
print(f"Frequência máxima absoluta para DHW > 4: {global_vmax_4:.2f}%")
print(f"Frequência máxima absoluta para DHW > 8: {global_vmax_8:.2f}%")
print(f"==> Usando VMAX = {vmax_plot_4} para todos os mapas de DHW > 4")
print(f"==> Usando VMAX = {vmax_plot_8} para todos os mapas de DHW > 8")
print("-----------------------------------------")


# ============================================================
# Função de plotagem MODIFICADA para aceitar VMAX
# ============================================================
def plot_dhw_map(dhw_da, threshold, year_label, out_dir, vmax=100): # Adiciona vmax como argumento
    fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={'projection': ccrs.PlateCarree()})
    ax.set_extent([lon_min, lon_max, lat_min, lat_max])
    ax.set_facecolor('white')

    mask_interp = marine_mask.interp(lat=dhw_da.lat, lon=dhw_da.lon, method='nearest')
    dhw_masked = dhw_da.where(mask_interp >= 0.5)

    cmap = LinearSegmentedColormap.from_list('custom_warm', ['white', 'yellow', 'orange', 'red'])
    
    # Usa o vmax passado como argumento, garantindo que o valor mínimo seja pelo menos 1 para evitar erro
    effective_vmax = max(vmax, 1)

    img = ax.pcolormesh(
        dhw_da.lon,
        dhw_da.lat,
        dhw_masked,
        cmap=cmap,
        transform=ccrs.PlateCarree(),
        vmin=0,
        vmax=effective_vmax, # USA O VMAX AJUSTADO
        zorder=1
    )

    cbar = plt.colorbar(img, ax=ax, pad=0.05)
    cbar.set_label(f"Frequência de dias com DHW > {threshold} (%)")
    
    # ... (o resto da função de plotagem é igual)
    protected_marine = gpd.overlay(gdf_protected, gdf_land, how='difference')
    if not protected_marine.empty:
        protected_marine.boundary.plot(ax=ax, edgecolor='darkblue', linewidth=1.5, transform=ccrs.PlateCarree(), zorder=5)
    gdf_land.plot(ax=ax, facecolor='white', edgecolor='black', linewidth=1, transform=ccrs.PlateCarree(), zorder=10)
    ax.add_geometries(gdf_coast.geometry, crs=ccrs.PlateCarree(), facecolor='none', edgecolor='black', linewidth=1.5, zorder=15)
    bathy_ocean = bathy.where(bathy < 0)
    ax.contour(bathy_ocean.lon, bathy_ocean.lat, bathy_ocean, levels=[-200], colors='black', linewidths=1, transform=ccrs.PlateCarree(), zorder=12)
    gdf_islands_bound.boundary.plot(ax=ax, edgecolor='black', linewidth=1, transform=ccrs.PlateCarree(), zorder=25)
    gdf_recifes.boundary.plot(ax=ax, edgecolor='purple', linewidth=0.5, transform=ccrs.PlateCarree(), zorder=30)
    protected_handle = mlines.Line2D([], [], color='darkblue', linewidth=1.5, label='Área Protegida (marinha)')
    reef_handle = mlines.Line2D([], [], color='purple', linewidth=0.5, label='Recifes')
    ax.legend(handles=[protected_handle, reef_handle], loc='upper right', fontsize=8)
    ax.set_title(f"Frequência de dias com DHW > {threshold} ({year_label})")
    fname = os.path.join(out_dir, f"dhw_frequency_gt{threshold}_{year_label}.png")
    plt.savefig(fname, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Mapa salvo: {fname}")

# ============================================================
# Geração dos mapas anuais (período Jul-Jun) COM ESCALA CONSISTENTE
# ============================================================
print("\nGerando mapas para os períodos de Julho a Junho...")

for t in hydro_yearly_percent_4.time.values:
    start_year = pd.to_datetime(t).year
    if start_year not in allowed_years:
        continue
    end_year = start_year + 1
    year_label = f"{start_year}-{end_year}"
    
    data_slice_4 = hydro_yearly_percent_4.sel(time=t)
    data_slice_8 = hydro_yearly_percent_8.sel(time=t)
    
    # Chama a função de plotagem passando o VMAX GLOBAL calculado
    plot_dhw_map(data_slice_4, 4, year_label, output_dir, vmax=vmax_plot_4)
    plot_dhw_map(data_slice_8, 8, year_label, output_dir, vmax=vmax_plot_8)

# ============================================================
# (Opcional) Mapas agregados para todo o período
# ============================================================
start_period = pd.to_datetime(hydro_yearly_percent_4.time.values[0]).year
end_period = pd.to_datetime(hydro_yearly_percent_4.time.values[-1]).year + 1
total_period_label = f"{start_period}-{end_period}"

mean_freq_4 = hydro_yearly_percent_4.mean(dim='time')
mean_freq_8 = hydro_yearly_percent_8.mean(dim='time')

# Para o mapa de média, também usamos a escala global para consistência
plot_dhw_map(mean_freq_4, 4, f"Média {total_period_label}", output_dir, vmax=vmax_plot_4)
plot_dhw_map(mean_freq_8, 8, f"Média {total_period_label}", output_dir, vmax=vmax_plot_8)

print(f"\nProcessamento concluído com sucesso – mapas para períodos de {total_period_label} gerados.")


# ============================================================
# ### Geração dos GIFs animados
# ============================================================
# Esta seção requer a biblioteca Pillow. Se não a tiver, instale com: pip install Pillow
from PIL import Image
import glob

def create_dhw_gif(output_directory, threshold, year_range, duration_per_frame_ms=500):
    """
    Encontra todos os mapas PNG para um limiar de DHW, os ordena e cria um GIF animado.

    Args:
        output_directory (str): Pasta onde os arquivos PNG estão e onde o GIF será salvo.
        threshold (int): O limiar de DHW (4 ou 8) para procurar nos nomes dos arquivos.
        year_range (set): O conjunto de anos permitidos para garantir a ordem correta.
        duration_per_frame_ms (int): Duração de cada quadro no GIF, em milissegundos.
    """
    print(f"\n--- Criando GIF para DHW > {threshold} ---")
    
    # 1. Encontrar todos os arquivos de imagem para o limiar especificado
    search_pattern = os.path.join(output_directory, f"dhw_frequency_gt{threshold}_*.png")
    image_files = sorted(glob.glob(search_pattern))

    # Filtra para garantir que apenas os anos processados sejam incluídos, em ordem
    # Isso é uma segurança extra para garantir a ordem cronológica correta
    valid_files_for_gif = []
    for year in sorted(list(year_range)):
        expected_filename = os.path.join(output_directory, f"dhw_frequency_gt{threshold}_{year}-{year+1}.png")
        if expected_filename in image_files:
            valid_files_for_gif.append(expected_filename)

    if not valid_files_for_gif:
        print(f"Nenhuma imagem encontrada para o limiar {threshold}. GIF não será criado.")
        return

    print(f"Encontradas {len(valid_files_for_gif)} imagens para a animação.")

    # 2. Ler as imagens e armazená-las como frames
    frames = []
    for filename in valid_files_for_gif:
        try:
            frames.append(Image.open(filename))
        except Exception as e:
            print(f"  [Aviso] Não foi possível abrir o arquivo {filename}. Pulando. Erro: {e}")

    if len(frames) < 2:
        print("São necessárias pelo menos 2 imagens para criar um GIF. Abortando.")
        return

    # 3. Salvar os frames como um GIF animado
    # O primeiro frame é usado para iniciar o arquivo, e o resto é anexado.
    first_frame = frames[0]
    output_gif_path = os.path.join(output_directory, f"animacao_dhw_gt{threshold}_{min(year_range)}-{max(year_range)+1}.gif")
    
    first_frame.save(
        output_gif_path,
        format='GIF',
        append_images=frames[1:],  # Anexa os frames restantes
        save_all=True,
        duration=duration_per_frame_ms,  # Duração de cada frame em ms (500ms = 0.5s)
        loop=0  # 0 significa que o GIF irá repetir infinitamente
    )
    
    print(f">>> GIF salvo com sucesso em: {output_gif_path}")

# --- Chamadas da função para criar os GIFs ---
# Lembre-se de passar o mesmo 'allowed_years' que você usou na análise
create_dhw_gif(output_dir, 4, allowed_years, duration_per_frame_ms=500)
create_dhw_gif(output_dir, 8, allowed_years, duration_per_frame_ms=500)