# Plano de Implementação: Mapas Globais de Magnitude e Variabilidade

**Objetivo:** Modificar o script `PCA_GLOBAL_maps_COMPOSITE.py` para gerar **duas figuras separadas**:
1. **Mapa de Magnitude** (2x2): Médias históricas de SST, DLI, Clorofila e proporção de DHW>4.
2. **Mapa de Variabilidade** (3x3): CV de SST, DLI e Clorofila para 3 janelas temporais (2, 30, all).

**Saída:** Pasta dedicada `Global_maps/` contendo `Fig_Magnitude.png` e `Fig_Variability.png`.

---

## 1. Contexto e Arquivos Envolvidos

| Arquivo | Caminho | Ação |
|:--------|:--------|:-----|
| Script a Modificar | [PCA_GLOBAL_maps_COMPOSITE.py](file:///c:/Users/rbfra/ONEDRIVE_NOVO/OneDrive/%23%23%23%23%23%23%23%23PUBLICACOES/%23%23%23%23%23%23%23%23%23%23%23%23Menezes%20et%20al.%20Mus%20his%20distribution%20and%20abundance%20Abrolhos/%23%23%23%23%23%23FINAL/%23%23%23%23%23%23%23FINAL_CODES/PCA_GLOBAL_maps_COMPOSITE.py) | Edição |
| Script de Referência | [PCA_LOCAL_bubbleplot.py](file:///c:/Users/rbfra/ONEDRIVE_NOVO/OneDrive/%23%23%23%23%23%23%23%23PUBLICACOES/%23%23%23%23%23%23%23%23%23%23%23%23Menezes%20et%20al.%20Mus%20his%20distribution%20and%20abundance%20Abrolhos/%23%23%23%23%23%23FINAL/%23%23%23%23%23%23%23FINAL_CODES/PCA_LOCAL_bubbleplot.py) | Leitura (lógica de CV) |
| Saída Magnitude | `#######FINAL_RESULTS/Global_maps/Fig_Magnitude.png` | Criação |
| Saída Variabilidade | `#######FINAL_RESULTS/Global_maps/Fig_Variability.png` | Criação |

---

## 2. Lógica do Cálculo de CV (Referência: PCA_LOCAL_bubbleplot.py)

> [!IMPORTANT]
> Esta seção explica a lógica de janelas deslizantes que deve ser replicada para o mapa de variabilidade espacial.

### 2.1. CV "all" (Toda a Série Temporal)
```python
# Fonte: PCA_LOCAL_bubbleplot.py, linhas 156-159
mean_vals = site_timeseries.mean('time', skipna=True)
std_vals = site_timeseries.std('time', skipna=True)
cv_all = (std_vals / (mean_vals + 1e-9) * 100).values
```
**Interpretação:** O CV é calculado usando TODA a série temporal disponível. Cada pixel recebe um único valor de CV.

### 2.2. CV com Janela Deslizante (ex: 2 ou 30 dias)
```python
# Fonte: PCA_LOCAL_bubbleplot.py, linhas 161-167
window = 30  # Exemplo: janela de 30 dias
min_p = max(2, int(window * 0.25))  # Mínimo de 25% de dados válidos na janela

rolling_mean = site_timeseries.rolling(time=window, min_periods=min_p, center=True).mean()
rolling_std = site_timeseries.rolling(time=window, min_periods=min_p, center=True).std()
cv_ts = (rolling_std / (rolling_mean + 1e-9)) * 100  # Série temporal de CV

# Média temporal do CV: cada pixel recebe a MÉDIA de todos os CVs calculados ao longo do tempo
cv_mean = cv_ts.mean('time', skipna=True).values
```
**Interpretação:** 
1. Para cada passo de tempo `t`, calcula-se o CV usando uma janela centrada de `window` dias.
2. O resultado é uma série temporal de valores de CV.
3. O valor final para cada pixel é a **média temporal** desta série de CVs.

---

## 3. Mudanças na Estrutura do Script

### 3.1. Novo Diretório de Saída (Linha ~79)

**Localização:** Substituir o valor de `output_dir`.

**Código Atual (linha 79):**
```python
output_dir = r"C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_RESULTS\PCA_GLOBAL_COMPOSITE"
```

**Novo Código:**
```python
output_dir = r"C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_RESULTS\Global_maps"
```

### 3.2. Remover `analysis_selector` (Linhas 93-98)

O seletor de análise não é mais necessário, pois **todas as 3 janelas** serão calculadas para a figura de variabilidade.

**Código a REMOVER:**
```python
# --- SELETOR ÚNICO DE ANÁLISE ---
analysis_selector = 'all'  

print(f"============================================================")
print(f"==> ANÁLISE CONFIGURADA PARA A OPÇÃO DE CV: '{analysis_selector}' <==")
print(f"============================================================")
```

---

## 4. Nova Função: Cálculo de CV por Janela Deslizante (Espacial)

**Inserir ANTES da função `create_publication_composite` (antes da linha 303).**

```python
def calculate_cv_map(data_da, window_size):
    """
    Calcula um mapa espacial de CV usando janela deslizante temporal.
    
    Para window_size == 'all':
        CV = (std temporal / média temporal) * 100
    Para window_size numérico (ex: 2, 30):
        1. Calcula CV(t) para cada janela centrada de tamanho `window_size`.
        2. Retorna a MÉDIA TEMPORAL de CV(t).
    
    Args:
        data_da (xr.DataArray): Dados com dimensão 'time'.
        window_size (int or str): Tamanho da janela ('all' ou inteiro).
    
    Returns:
        xr.DataArray: Mapa 2D (lat, lon) com valores de CV em %.
    """
    epsilon = 1e-9
    
    if window_size == 'all':
        # CV global: usa toda a série temporal
        mean_map = data_da.mean('time', skipna=True)
        std_map = data_da.std('time', skipna=True)
        cv_map = (std_map / (mean_map + epsilon)) * 100
        return cv_map.compute()
    else:
        # CV com janela deslizante
        window = int(window_size)
        min_periods = max(2, int(window * 0.25))  # 25% de dados mínimos
        
        # Processa em chunks para evitar estouro de memória
        with dask.config.set(scheduler='single-threaded'):
            rolling_mean = data_da.rolling(time=window, min_periods=min_periods, center=True).mean()
            rolling_std = data_da.rolling(time=window, min_periods=min_periods, center=True).std()
            cv_timeseries = (rolling_std / (rolling_mean + epsilon)) * 100
            
            # Média temporal do CV: valor final por pixel
            cv_map = cv_timeseries.mean('time', skipna=True)
            return cv_map.compute()
```

---

## 5. Nova Função: Mapa de Magnitude (2x2)

**Inserir APÓS a nova função `calculate_cv_map`.**

```python
def create_magnitude_map(data_dict, output_path, mask_shallow_data, bathy_data):
    """
    Gera figura 2x2 com mapas de magnitude (médias históricas).
    
    Layout:
        (a) Mean SST        | (b) Mean DLI
        (c) Mean Chl-a      | (d) DHW > 4 Proportion
    
    Args:
        data_dict: Dict com chaves 'mean_SST', 'mean_DLI', 'mean_CHL', 'prop_DHW_gt4'.
        output_path: Caminho de saída para PNG/PDF.
        mask_shallow_data: Máscara de águas rasas.
        bathy_data: Dados de batimetria para contornos.
    """
    panels_config = [
        ('mean_SST', 'Mean SST', 'viridis', '°C'),
        ('mean_DLI', 'Mean Benthic DLI', 'magma', r'mol m$^{-2}$ d$^{-1}$'),
        ('mean_CHL', 'Mean Chlorophyll-a', 'YlGn', r'mg m$^{-3}$'),
        ('prop_DHW_gt4', 'Thermal Stress (DHW > 4)', 'YlOrRd', '% of time'),
    ]
    
    PANEL_LABELS_MAG = ['a', 'b', 'c', 'd']
    FIGSIZE_MAG = (180 / 25.4, 160 / 25.4)  # ~7.1 x 6.3 inches
    
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
        
        # Aplicar máscara de águas rasas
        if mask_shallow_data is not None:
            mask_aligned = mask_shallow_data.interp_like(da, method='nearest')
            da = da.where(mask_aligned >= 0.5)
        
        data_to_plot = da.load() if hasattr(da, 'load') else da
        valid_data = data_to_plot.values[np.isfinite(data_to_plot.values)]
        vmin, vmax = (np.nanpercentile(valid_data, 2), np.nanpercentile(valid_data, 98)) if valid_data.size > 0 else (0, 1)
        
        cmap_obj = plt.get_cmap(cmap_name).copy()
        cmap_obj.set_bad(color='lightgray')
        
        img = ax.pcolormesh(data_to_plot.lon, data_to_plot.lat, data_to_plot,
                            cmap=cmap_obj, shading='auto', vmin=vmin, vmax=vmax)
        
        # Elementos geográficos
        if gdf_land is not None:
            gdf_land.plot(ax=ax, facecolor='#D2B48C', edgecolor='black', linewidth=0.5, zorder=3)
        if bathy_data is not None:
            ax.contour(bathy_data.lon, bathy_data.lat, bathy_data, levels=[-200, -50],
                       colors=['dimgray', 'gray'], linewidths=[0.7, 0.5], linestyles=['--', ':'], zorder=2)
        
        # Colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.08)
        cbar = fig.colorbar(img, cax=cax, orientation='vertical')
        cbar.set_label(unit, fontsize=8)
        cbar.ax.tick_params(labelsize=7)
        
        # Label do painel
        ax.text(-0.15, 1.1, f'({PANEL_LABELS_MAG[i]})', transform=ax.transAxes, fontsize=11, fontweight='bold', va='top')
        
        # Configuração de eixos
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
        
        # Barra de escala apenas no primeiro painel
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
```

---

## 6. Nova Função: Mapa de Variabilidade (3x3)

**Inserir APÓS a função `create_magnitude_map`.**

```python
def create_variability_map(cv_data_dict, output_path, mask_shallow_data, bathy_data):
    """
    Gera figura 3x3 com mapas de variabilidade (CV) para 3 variáveis x 3 janelas.
    
    Layout (colunas = janelas, linhas = variáveis):
                    | CV Window 2  | CV Window 30 | CV All
        SST         | (a)          | (b)          | (c)
        DLI         | (d)          | (e)          | (f)
        Chlorophyll | (g)          | (h)          | (i)
    
    Args:
        cv_data_dict: Dict com chaves no formato 'cv_{var}_{window}'.
            Exemplo: 'cv_SST_2', 'cv_SST_30', 'cv_SST_all', 'cv_DLI_2', etc.
        output_path: Caminho de saída para PNG/PDF.
        mask_shallow_data: Máscara de águas rasas.
        bathy_data: Dados de batimetria para contornos.
    """
    # Configuração dos painéis: (chave no dict, título, cmap)
    variables = ['SST', 'DLI', 'CHL']
    windows = [2, 30, 'all']
    
    # Títulos das colunas e linhas
    col_titles = ['CV (2-day window)', 'CV (30-day window)', 'CV (full series)']
    row_titles = ['SST', 'Benthic DLI', 'Chlorophyll-a']
    
    PANEL_LABELS_VAR = [chr(ord('a') + i) for i in range(9)]  # 'a' to 'i'
    FIGSIZE_VAR = (250 / 25.4, 220 / 25.4)  # ~9.8 x 8.7 inches
    
    fig = plt.figure(figsize=FIGSIZE_VAR, constrained_layout=False)
    gs = gridspec.GridSpec(3, 3, figure=fig, wspace=0.20, hspace=0.25)
    
    gdf_land = load_coast_and_islands(coast_shp, islands_shp, lat_min, lat_max, lon_min, lon_max)
    
    panel_idx = 0
    for row_idx, var in enumerate(variables):
        for col_idx, window in enumerate(windows):
            ax = fig.add_subplot(gs[row_idx, col_idx])
            
            # Chave do dicionário de dados
            key = f'cv_{var}_{window}'
            da = cv_data_dict.get(key)
            
            if da is None:
                ax.text(0.5, 0.5, 'Data unavailable', ha='center', va='center', transform=ax.transAxes)
            else:
                # Aplicar máscara de águas rasas
                if mask_shallow_data is not None:
                    mask_aligned = mask_shallow_data.interp_like(da, method='nearest')
                    da = da.where(mask_aligned >= 0.5)
                
                data_to_plot = da.load() if hasattr(da, 'load') else da
                valid_data = data_to_plot.values[np.isfinite(data_to_plot.values)]
                vmin, vmax = (np.nanpercentile(valid_data, 2), np.nanpercentile(valid_data, 98)) if valid_data.size > 0 else (0, 100)
                
                # Usar cmap consistente para CV (variabilidade)
                cmap_obj = plt.get_cmap('cividis').copy()
                cmap_obj.set_bad(color='lightgray')
                
                img = ax.pcolormesh(data_to_plot.lon, data_to_plot.lat, data_to_plot,
                                    cmap=cmap_obj, shading='auto', vmin=vmin, vmax=vmax)
                
                # Elementos geográficos
                if gdf_land is not None:
                    gdf_land.plot(ax=ax, facecolor='#D2B48C', edgecolor='black', linewidth=0.5, zorder=3)
                if bathy_data is not None:
                    ax.contour(bathy_data.lon, bathy_data.lat, bathy_data, levels=[-200, -50],
                               colors=['dimgray', 'gray'], linewidths=[0.7, 0.5], linestyles=['--', ':'], zorder=2)
                
                # Colorbar
                divider = make_axes_locatable(ax)
                cax = divider.append_axes("right", size="5%", pad=0.05)
                cbar = fig.colorbar(img, cax=cax, orientation='vertical')
                cbar.set_label('CV (%)', fontsize=7)
                cbar.ax.tick_params(labelsize=6)
            
            # Label do painel
            ax.text(-0.12, 1.08, f'({PANEL_LABELS_VAR[panel_idx]})', transform=ax.transAxes, 
                    fontsize=10, fontweight='bold', va='top')
            
            # Título da coluna (apenas na primeira linha)
            if row_idx == 0:
                ax.set_title(col_titles[col_idx], fontsize=9, pad=8)
            
            # Título da linha (texto à esquerda, apenas na primeira coluna)
            if col_idx == 0:
                ax.text(-0.25, 0.5, row_titles[row_idx], transform=ax.transAxes, 
                        fontsize=10, fontweight='bold', va='center', ha='right', rotation=90)
            
            # Configuração de eixos
            ax.set_xlim(lon_min, lon_max)
            ax.set_ylim(lat_min, lat_max)
            ax.set_aspect('equal', adjustable='box')
            
            # Labels de eixo apenas nas bordas externas
            if col_idx == 0:
                ax.set_ylabel('Lat (°S)', fontsize=8)
            else:
                ax.set_yticklabels([])
            if row_idx == 2:
                ax.set_xlabel('Lon (°W)', fontsize=8)
            else:
                ax.set_xticklabels([])
            
            ax.tick_params(axis='both', labelsize=6)
            
            # Barra de escala apenas no primeiro painel
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
```

---

## 7. Modificação do Bloco de Processamento Principal

### 7.1. Cálculo de Métricas para TODAS as 3 Janelas

**Substituir o bloco de cálculo de CV (linhas 434-457) pelo seguinte código:**

```python
print("\n--- 3. Calculando Métricas Espaciais ---")
computed_results = {}

# === MAGNITUDE ===
print("  Calculando métricas de MAGNITUDE...")
computed_results['mean_SST'] = sst_all.mean('time', skipna=True).compute()
if chl_all is not None:
    computed_results['mean_CHL'] = chl_all.mean('time', skipna=True).compute()
if dli_all is not None:
    computed_results['mean_DLI'] = dli_all.mean('time', skipna=True).compute()
if dhw_all is not None:
    computed_results['prop_DHW_gt4'] = (dhw_all > 4).mean('time', skipna=True).compute() * 100

# === VARIABILIDADE (3 janelas: 2, 30, all) ===
print("  Calculando métricas de VARIABILIDADE...")
windows_to_compute = [2, 30, 'all']

# SST CV
print("    - SST CV para todas as janelas...")
for window in windows_to_compute:
    key = f'cv_SST_{window}'
    print(f"      Computando {key}...")
    computed_results[key] = calculate_cv_map(sst_all, window)
    gc.collect()

# DLI CV
if dli_all is not None:
    print("    - DLI CV para todas as janelas...")
    for window in windows_to_compute:
        key = f'cv_DLI_{window}'
        print(f"      Computando {key}...")
        computed_results[key] = calculate_cv_map(dli_all, window)
        gc.collect()

# CHL CV
if chl_all is not None:
    print("    - CHL CV para todas as janelas...")
    for window in windows_to_compute:
        key = f'cv_CHL_{window}'
        print(f"      Computando {key}...")
        computed_results[key] = calculate_cv_map(chl_all, window)
        gc.collect()

print("  Cálculo de métricas concluído!")
```

### 7.2. Chamada das Novas Funções de Plotagem

**Substituir o bloco de chamada da figura composta (linhas 459-475) pelo seguinte:**

```python
# === GERAÇÃO DAS FIGURAS ===
print("\n--- 4. Gerando Figuras ---")

# Dados para o mapa de MAGNITUDE
magnitude_data = {
    'mean_SST': computed_results.get('mean_SST'),
    'mean_DLI': computed_results.get('mean_DLI'),
    'mean_CHL': computed_results.get('mean_CHL'),
    'prop_DHW_gt4': computed_results.get('prop_DHW_gt4'),
}

# Dados para o mapa de VARIABILIDADE
variability_data = {}
for var in ['SST', 'DLI', 'CHL']:
    for window in [2, 30, 'all']:
        key = f'cv_{var}_{window}'
        variability_data[key] = computed_results.get(key)

# Gerar Mapa de Magnitude (2x2)
magnitude_output_path = os.path.join(output_dir, "Fig_Magnitude.png")
create_magnitude_map(magnitude_data, magnitude_output_path, mask_shallow, bathy_data)

# Gerar Mapa de Variabilidade (3x3)
variability_output_path = os.path.join(output_dir, "Fig_Variability.png")
create_variability_map(variability_data, variability_output_path, mask_shallow, bathy_data)

print("\n=== PROCESSAMENTO FINALIZADO COM SUCESSO ===")
print(f"Figuras salvas em: {output_dir}")
```

---

## 8. Código a Remover

### 8.1. Remover Função `create_publication_composite` (linhas 303-384)

Esta função será substituída pelas novas funções `create_magnitude_map` e `create_variability_map`.

### 8.2. Remover/Comentar Código de Planilhas e PCA (linhas 477-488)

O código antigo de exportação de planilhas e PCA opcional pode ser comentado ou removido, pois não é mais relevante para este script focado em mapas.

---

## 9. Checklist de Implementação

- [ ] **Atualizar `output_dir`** para `Global_maps` (Seção 3.1)
- [ ] **Remover `analysis_selector`** (Seção 3.2)
- [ ] **Adicionar função `calculate_cv_map`** (Seção 4)
- [ ] **Adicionar função `create_magnitude_map`** (Seção 5)
- [ ] **Adicionar função `create_variability_map`** (Seção 6)
- [ ] **Substituir bloco de cálculo de métricas** (Seção 7.1)
- [ ] **Substituir bloco de geração de figuras** (Seção 7.2)
- [ ] **Remover `create_publication_composite`** (Seção 8.1)
- [ ] **Remover/comentar código de planilhas** (Seção 8.2)

---

## 10. Verificação

### 10.1. Comando de Execução
```powershell
cd "c:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_CODES"
python PCA_GLOBAL_maps_COMPOSITE.py
```

### 10.2. Critérios de Sucesso

1. **Sem erros de execução**: O script deve completar sem exceções.
2. **Arquivos de saída criados**:
   - `#######FINAL_RESULTS/Global_maps/Fig_Magnitude.png` (e `.pdf`)
   - `#######FINAL_RESULTS/Global_maps/Fig_Variability.png` (e `.pdf`)
3. **Verificação visual**:
   - **Fig_Magnitude**: 4 painéis (2x2) com médias de SST, DLI, CHL e proporção DHW>4.
   - **Fig_Variability**: 9 painéis (3x3) com CVs de SST/DLI/CHL para janelas 2/30/all.
4. **Qualidade**: Figuras em 600 DPI, colorbars legíveis, labels (a-i) visíveis.

---

## 11. Notas Técnicas

> [!WARNING]
> O cálculo de CV com janelas deslizantes (especialmente janela 2) pode consumir muita memória devido ao rolling sobre dados espaciais. O script usa `scheduler='single-threaded'` e `gc.collect()` para mitigar isso.

> [!TIP]
> Se houver problemas de memória, reduza o período de dados ou aumente o chunk size temporal.
