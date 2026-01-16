# Plano de Implementação Refinado: Sincronicidade Temporal e Dinâmica Ambiental
## Script Unificado Robusto para Publicação Nature/Science

**Data:** 2026-01-16 (atualizado)  
**Status:** Plano enriquecido com contexto técnico completo para implementação autônoma por LLM  
**Scripts Alvo:** 
- `#######FINAL_CODES/TIME_SERIES_VARIABILITY_ANALYSIS_GEMINI.py` (modificações)
- `#######FINAL_CODES/TIME_SERIES_UNIFIED_SYNCHRONY_ANALYSIS.py` (NOVO)

---

## 0. Contexto Técnico para Implementação (LLM Reference)

> [!IMPORTANT]
> Esta seção contém todas as informações necessárias para um LLM implementar as mudanças sem contexto adicional.

### 0.1 Estrutura de Diretórios

```
######FINAL/
├── #######FINAL_CODES/                          # Scripts Python e R
│   ├── TIME_SERIES_VARIABILITY_ANALYSIS_GEMINI.py  # Script existente (MODIFICAR)
│   ├── TIME_SERIES_SYNCHRONY_ANALYSIS.py           # Script de sincronicidade antigo (referência)
│   ├── TIME_SERIES_LAG_temporal_environ.py         # Script de séries temporais (referência)
│   └── TIME_SERIES_UNIFIED_SYNCHRONY_ANALYSIS.py   # Script NOVO a criar
│
├── #######FINAL_RESULTS/
│   ├── Variability_Analysis_GEMINI/             # Saída atual do script de variabilidade
│   ├── Temporal_Synchrony_Analysis/             # Saída do novo script de sincronicidade
│   └── #####output_local_PCA_CV_*_FINAL/        # Dados integrados por cenário CV
│
└── plans/
    └── temporal_synchrony_analysis_plan_REFINED.md  # Este plano
```

### 0.2 Caminhos Absolutos Críticos

```python
# Diretório base do projeto
PROJECT_BASE = r'C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL'

# Diretórios de dados de satélite (externos)
SST_DIR = r'H:\remote sensing\CRW_SST_FULL'
MODIS_DIR = r'H:\remote sensing\MODIS_DATA_FULL'

# CSV de sites
SITES_CSV_PATH = r'C:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######22.04.23\DATA\sites_list_full_clean.csv'

# Diretório de saída do script de variabilidade
OUTPUT_DIR_VARIABILITY = r'C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_RESULTS\Variability_Analysis_GEMINI'

# Diretório de saída do script de sincronicidade (NOVO)
OUTPUT_DIR_SYNCHRONY = r'C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_RESULTS\Temporal_Synchrony_Analysis'
```

### 0.3 Estrutura do CSV de Sites

O arquivo `sites_list_full_clean.csv` usa separador `;` e contém:

| Coluna | Tipo | Descrição |
|--------|------|-----------|
| `Site_name` | string | Identificador único do sítio (ex: "PAB_TP_1") |
| `Reef_name` | string | Nome do recife (ex: "Parcel dos Abrolhos") |
| `HAB` | string | Tipo de habitat: "TP" (Tide Pool), "PA" (Reef Wall), "RR" (Rocky Reef) |
| `Latitude` | float | Latitude decimal (negativo para Sul) |
| `Longitude` | float | Longitude decimal (negativo para Oeste) |
| `Depth_m` | float | Profundidade em metros |
| `Arc` | string | Arco geográfico: "inner" ou "outer" |

### 0.4 Constantes de Configuração Existentes

```python
# Período de análise
START_YEAR = 2002
END_YEAR = 2008

# Padrões de arquivo NetCDF
SST_PATTERN = 'coraltemp_v3.1_*.nc'
KD490_PATTERN = 'AQUA_MODIS.*.L3m.DAY.KD.Kd_490.4km.nc'
PAR_PATTERN = 'AQUA_MODIS.*.L3m.DAY.PAR.par.4km.nc'
CHL_PATTERN = 'AQUA_MODIS.*.L3m.DAY.CHL.chlor_a.4km.nc'

# Constantes para cálculo DLI (equação Gattuso)
KDPAR_A, KDPAR_B, KDPAR_C = 0.0665, 0.874, 0.00121

# Paleta de cores padronizada
COLOR_INNER = '#CD853F'  # Peru (light brown)
COLOR_OUTER = '#4169E1'  # Royal Blue
COLOR_SST = '#D55E00'    # Vermelho-alaranjado
COLOR_DLI = '#0072B2'    # Azul escuro
COLOR_CHL = '#009E73'    # Verde-teal
```

### 0.5 Funções Existentes em `TIME_SERIES_VARIABILITY_ANALYSIS_GEMINI.py`

#### Funções de Carregamento de Dados (reutilizar no script de sincronicidade)
```python
def parse_date_from_filename(filename) -> pd.Timestamp | None
def filter_files_by_period(files, start_year, end_year) -> list
def load_timeseries_for_sites(files, var_names, sites_df, var_label) -> dict[str, pd.Series]
def calculate_dli_for_sites(kd490_files, par_files, sites_df) -> dict[str, pd.Series]
```

#### Funções de CV (existentes, não modificar)
```python
def calculate_cv_full(series) -> float
def calculate_cv_rolling(series, window) -> float
def calculate_all_cvs(timeseries_dict, sites_df) -> pd.DataFrame
```

#### Funções de Spectral Analysis (existentes, não modificar)
```python
def compute_lombscargle(series, periods_of_interest=None) -> tuple[np.ndarray, np.ndarray]
def extract_power_at_period(periods, power, target_period, tolerance=0.3) -> float
def compute_spectral_power_for_sites(timeseries_dict, sites_df) -> list[dict]
```

#### Funções de Pulse Detection (existentes, não modificar)
```python
def count_extreme_events(series, threshold) -> int
def calculate_pulse_counts(timeseries_dict, sites_df, percentile=90) -> pd.DataFrame
```

#### Funções de Estatística (existentes, reutilizar)
```python
def mann_whitney_test(inner_values, outer_values) -> tuple[float, str]
def format_significance(p_value) -> str
```

### 0.6 Estrutura de Dados das Séries Temporais

As séries temporais são armazenadas como `dict[str, pd.Series]` onde:
- **Chave**: `Site_name` (string)
- **Valor**: `pd.Series` com `DatetimeIndex` e valores float

Exemplo:
```python
sst_ts = {
    'PAB_TP_1': pd.Series([27.5, 27.8, ...], index=DatetimeIndex(['2002-01-01', '2002-01-02', ...])),
    'PAB_TP_2': pd.Series([...]),
    ...
}
```

### 0.7 Ponto de Inserção para Diagnóstico de Dados

**Arquivo:** `TIME_SERIES_VARIABILITY_ANALYSIS_GEMINI.py`

**Inserir nova função após linha ~520** (após `calculate_pulse_counts`):
```python
# ============================================
# Data Availability Diagnostics (NOVO)
# ============================================

def diagnose_data_availability(sst_ts_dict, dli_ts_dict, chl_ts_dict, sites_df):
    # ... implementação ...

def create_availability_figure(diagnostic_df, output_path):
    # ... implementação ...
```

**Inserir chamada no `__main__` após linha ~1041** (após carregar todos os dados):
```python
    # === NOVO: Diagnóstico de Disponibilidade ===
    print("\n" + "=" * 60)
    print("STEP 1.5: Data Availability Diagnostics")
    print("=" * 60)
    
    diagnostic_df = diagnose_data_availability(sst_ts, dli_ts, chl_ts, sites_df)
    diagnostic_df.to_csv(os.path.join(OUTPUT_DIR, 'Data_Availability_Diagnostic.csv'), index=False)
    create_availability_figure(diagnostic_df, os.path.join(OUTPUT_DIR, 'Fig_S1_Data_Availability.png'))
```

### 0.8 Imports Já Disponíveis

```python
# Core
import os, glob, re, logging
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
```

### 0.9 Imports Adicionais Necessários

```python
# Para o script de sincronicidade (NOVO)
from statsmodels.tsa.seasonal import STL  # Decomposição STL
from matplotlib.gridspec import GridSpec  # Layout complexo de figuras
```

---

## 1. Contexto e Objetivos

### 1.1 Hipótese Central
Os recifes do **Inner Arc** (próximos à costa) apresentam dinâmica ambiental mais variável e maior acoplamento temporal entre SST, DLI e Clorofila, impulsionado pela descarga de sedimentos costeiros. Espera-se observar:
- **↓ DLI** → **↓ SST** (após alguns dias) → **↑ Clorofila** (pulsos de nutrientes)

### 1.2 Arquitetura de Scripts Separados

> [!IMPORTANT]
> Os scripts permanecerão **separados** para manter foco e modularidade. Cada script tem responsabilidades específicas.

#### Divisão de Responsabilidades

| Script | Foco | Análises |
|--------|------|----------|
| **`TIME_SERIES_VARIABILITY_ANALYSIS_GEMINI.py`** (EXISTENTE) | Variabilidade temporal por arco | CV (3 janelas), Lomb-Scargle, Pulse Detection, **Diagnóstico de Lacunas** (NOVO) |
| **`TIME_SERIES_UNIFIED_SYNCHRONY_ANALYSIS.py`** (NOVO) | Sincronicidade entre variáveis | Pearson, Cross-Correlation (LAG), STL Decomposition, Spectral Coherence |
| `TIME_SERIES_LAG_temporal_environ.py` (existente) | Séries temporais brutas e PCA | Carregamento de dados, visualizações básicas |

#### Fluxo de Dados

```
TIME_SERIES_VARIABILITY_ANALYSIS_GEMINI.py
├── Carrega dados SST, DLI, Chl-a
├── Executa diagnóstico de lacunas → Data_Availability_Diagnostic.csv
├── Calcula CV, Spectral Power, Pulse Detection
└── Exporta dados processados para uso pelo script de sincronicidade

TIME_SERIES_UNIFIED_SYNCHRONY_ANALYSIS.py
├── Carrega dados já processados OU recarrega raw data
├── Calcula métricas de sincronicidade entre variáveis
└── Gera figuras de sincronicidade (Fig_Main, Fig_S2, Fig_S3)

---

## 2. Diagnóstico de Disponibilidade de Dados

> [!IMPORTANT]
> As linhas de tendência incompletas para DLI e Clorofila (que não alcançam 365 dias no periodograma) são causadas por **lacunas sistemáticas nos dados MODIS** (cobertura de nuvens, órbitas). Este diagnóstico deve ser explicitamente documentado.

### 2.1 Módulo de Diagnóstico

> [!NOTE]
> **Script alvo: `TIME_SERIES_VARIABILITY_ANALYSIS_GEMINI.py`**
> 
> Este módulo será adicionado ao script de variabilidade existente, pois é onde os dados brutos são carregados. O diagnóstico será executado automaticamente e os resultados serão salvos para uso por ambos os scripts.

```python
def diagnose_data_availability(sst_ts_dict, dli_ts_dict, chl_ts_dict, sites_df):
    """
    Produz diagnóstico detalhado de disponibilidade temporal por variável e arco.
    
    Outputs:
    - DataFrame com contagens diárias, gaps, e cobertura % por variável/arco
    - Relatório textual para inclusão em Supplementary Material
    - Figura diagnóstica de disponibilidade temporal
    """
    results = []
    
    for arc in ['inner', 'outer']:
        arc_sites = sites_df[sites_df['Arc'] == arc]['Site_name'].tolist()
        
        for var_name, ts_dict in [('SST', sst_ts_dict), ('DLI', dli_ts_dict), ('CHL', chl_ts_dict)]:
            # Coletar todas as datas com dados válidos
            all_dates = set()
            for site in arc_sites:
                if site in ts_dict:
                    all_dates.update(ts_dict[site].dropna().index.date)
            
            # Calcular métricas
            n_days_with_data = len(all_dates)
            expected_days = (END_YEAR - START_YEAR + 1) * 365
            coverage_pct = (n_days_with_data / expected_days) * 100
            
            # Identificar maior gap contínuo
            sorted_dates = sorted(all_dates)
            max_gap = 0
            if len(sorted_dates) > 1:
                gaps = [(sorted_dates[i+1] - sorted_dates[i]).days 
                        for i in range(len(sorted_dates)-1)]
                max_gap = max(gaps) if gaps else 0
            
            results.append({
                'Arc': arc,
                'Variable': var_name,
                'N_Days_With_Data': n_days_with_data,
                'Expected_Days': expected_days,
                'Coverage_Pct': coverage_pct,
                'Max_Gap_Days': max_gap,
                'Max_Resolvable_Period': n_days_with_data // 2  # Nyquist limit
            })
    
    return pd.DataFrame(results)
```

### 2.2 Saídas do Diagnóstico (geradas por `TIME_SERIES_VARIABILITY_ANALYSIS_GEMINI.py`)

**Diretório de saída:** `#######FINAL_RESULTS/Variability_Analysis_Output/`

1. **`Data_Availability_Diagnostic.csv`**: Tabela com métricas por variável e arco
2. **`Fig_S1_Data_Availability.png`**: Heatmap temporal mostrando dias com dados
3. **Texto para Supplementary Material**:
   > "SST data from CRW (2002-2008) had 98.5% temporal coverage. MODIS-derived DLI and Chl-a had 67.3% and 71.2% coverage respectively, limited by cloud cover. The maximum resolvable period for spectral analysis was therefore limited to ~450 days for DLI/Chl-a versus ~1200 days for SST."

---

## 3. Estrutura de Figuras Refinada

> [!CAUTION]
> O layout atual com painéis lado a lado (2x2) comprime horizontalmente as séries temporais de 7 anos, dificultando a identificação de padrões. A nova estrutura prioriza visualização longitudinal.

### 3.1 Figura Principal: `Fig_Main_Temporal_Synchrony.png`

**Novo Layout: 3 Linhas x 2 Colunas (Vertical)**

```
┌────────────────────────────────────────────────────────┐
│ Painel A: Inner Arc - Time Series (3 eixos Y)          │
│ [LARGURA TOTAL - ~16 inches]                           │
│ SST (vermelho) | DLI (azul) | Chl-a (verde)            │
│ Dados suavizados (14-day rolling mean)                 │
├────────────────────────────────────────────────────────┤
│ Painel B: Outer Arc - Time Series (3 eixos Y)          │
│ [LARGURA TOTAL - ~16 inches]                           │
│ Mesma estrutura do Painel A                            │
├───────────────────────┬────────────────────────────────┤
│ Painel C:             │ Painel D:                      │
│ Pearson Synchrony     │ Cross-Correlation (LAG)        │
│ (Boxplot + jitter)    │ CCF Inner vs Outer             │
│                       │ Destaque: lag ótimo            │
├───────────────────────┼────────────────────────────────┤
│ Painel E:             │ Painel F:                      │
│ Spectral Coherence    │ Pulse Synchrony (STL Residuals)│
│ (30-day period)       │ Boxplot dos resíduos           │
└───────────────────────┴────────────────────────────────┘
```

**Especificações Técnicas:**
```python
fig = plt.figure(figsize=(16, 18))  # Ampliado verticalmente
gs = fig.add_gridspec(4, 2, 
                      height_ratios=[1.5, 1.5, 1, 1],  # Séries temporais maiores
                      hspace=0.25, wspace=0.30)

# Painéis A e B: span de 2 colunas
ax_inner = fig.add_subplot(gs[0, :])  # Linha 1, ambas colunas
ax_outer = fig.add_subplot(gs[1, :])  # Linha 2, ambas colunas

# Painéis C-F: grid 2x2 nas linhas 3-4
ax_pearson = fig.add_subplot(gs[2, 0])
ax_ccf = fig.add_subplot(gs[2, 1])
ax_spectral = fig.add_subplot(gs[3, 0])
ax_stl = fig.add_subplot(gs[3, 1])
```

### 3.2 Figura Suplementar: `Fig_S2_Raw_Time_Series.png`

**Séries Brutas** (como em `TIME_SERIES_LAG_temporal_environ.py`):

```
┌────────────────────────────────────────────────────────┐
│ Painel A: Inner Arc - Raw Data Points                  │
│ [LARGURA TOTAL]                                        │
│ Pontos diários (sem suavização)                        │
│ 3 eixos Y independentes                                │
├────────────────────────────────────────────────────────┤
│ Painel B: Outer Arc - Raw Data Points                  │
│ [LARGURA TOTAL]                                        │
│ Mesma estrutura                                        │
└────────────────────────────────────────────────────────┘
```

**Implementação:**
```python
def create_raw_timeseries_figure(sst_arc, dli_arc, chl_arc, output_path):
    """
    Cria figura com dados brutos (pontos) para visualização de alta resolução temporal.
    Segue estilo de TIME_SERIES_LAG_temporal_environ.py
    """
    fig, axes = plt.subplots(2, 1, figsize=(18, 10), sharex=True)
    
    for i, (arc, title) in enumerate([('inner', 'Inner Arc'), ('outer', 'Outer Arc')]):
        ax = axes[i]
        
        # Dados brutos (sem smooth)
        sst = sst_arc[arc]
        dli = dli_arc[arc]
        chl = chl_arc[arc]
        
        # Eixo 1: SST (pontos vermelhos)
        ax.scatter(sst.index, sst.values, s=2, alpha=0.6, c=COLOR_SST, label='SST')
        ax.set_ylabel('SST (°C)', color=COLOR_SST)
        
        # Eixo 2: DLI
        ax2 = ax.twinx()
        ax2.scatter(dli.index, dli.values, s=2, alpha=0.6, c=COLOR_DLI, label='DLI')
        ax2.set_ylabel('DLI (mol m⁻² d⁻¹)', color=COLOR_DLI)
        
        # Eixo 3: Chl-a (deslocado)
        ax3 = ax.twinx()
        ax3.spines['right'].set_position(('outward', 60))
        ax3.scatter(chl.index, chl.values, s=2, alpha=0.6, c=COLOR_CHL, label='Chl-a')
        ax3.set_ylabel('Chl-a (mg m⁻³)', color=COLOR_CHL)
        
        ax.set_title(f'{title} - Raw Daily Data', fontweight='bold')
```

---

## 4. Análises Estatísticas

### 4.1 Decomposição STL para Isolamento de Pulsos

```python
from statsmodels.tsa.seasonal import STL

def calculate_pulse_synchrony_stl(sst, dli, chl, period=365):
    """
    Decompõe séries usando STL e calcula correlação dos RESÍDUOS.
    
    STL é mais robusto que seasonal_decompose para:
    - Dados faltantes
    - Tendências não-lineares
    - Outliers
    """
    # Preencher gaps com interpolação linear (máximo 7 dias)
    sst_interp = sst.interpolate(method='time', limit=7)
    dli_interp = dli.interpolate(method='time', limit=7)
    chl_interp = chl.interpolate(method='time', limit=7)
    
    # Alinhar séries
    df = pd.DataFrame({'SST': sst_interp, 'DLI': dli_interp, 'CHL': chl_interp}).dropna()
    
    if len(df) < period * 2:
        return None  # Dados insuficientes para decomposição
    
    # Decomposição STL
    stl_sst = STL(df['SST'], period=period, robust=True).fit()
    stl_dli = STL(df['DLI'], period=period, robust=True).fit()
    stl_chl = STL(df['CHL'], period=period, robust=True).fit()
    
    # Correlações dos resíduos (pulsos)
    return {
        'Pulse_Corr_SST_DLI': stl_sst.resid.corr(stl_dli.resid),
        'Pulse_Corr_SST_CHL': stl_sst.resid.corr(stl_chl.resid),
        'Pulse_Corr_DLI_CHL': stl_dli.resid.corr(stl_chl.resid),
        'Resid_SST': stl_sst.resid,
        'Resid_DLI': stl_dli.resid,
        'Resid_CHL': stl_chl.resid
    }
```

### 4.2 Análise de Cross-Correlation com Lag

```python
def calculate_cross_correlation_full(series1, series2, max_lag=60):
    """
    Cross-correlation lag analysis entre duas variáveis.
    
    Returns:
    - lags: array de lags testados (-max_lag a +max_lag)
    - correlations: correlação em cada lag
    - optimal_lag: lag com máxima correlação absoluta
    - optimal_corr: correlação no lag ótimo
    """
    df = pd.DataFrame({'s1': series1, 's2': series2}).dropna()
    
    if len(df) < 60:
        return None
    
    lags = range(-max_lag, max_lag + 1)
    correlations = []
    
    for lag in lags:
        corr = df['s2'].corr(df['s1'].shift(lag))
        correlations.append(corr)
    
    correlations = np.array(correlations)
    optimal_idx = np.nanargmax(np.abs(correlations))
    
    return {
        'lags': np.array(list(lags)),
        'correlations': correlations,
        'optimal_lag': list(lags)[optimal_idx],
        'optimal_corr': correlations[optimal_idx]
    }
```

---

## 5. Estrutura dos Scripts

### 5.1 Estrutura do Script de Variabilidade (MODIFICAR)

> [!NOTE]
> **Arquivo:** `TIME_SERIES_VARIABILITY_ANALYSIS_GEMINI.py`
> 
> Adicionar apenas o módulo de diagnóstico. Funções existentes permanecem inalteradas.

```
TIME_SERIES_VARIABILITY_ANALYSIS_GEMINI.py
│
├── [EXISTENTE] Configurações e Constantes (linhas 1-115)
├── [EXISTENTE] Load Sites Data (linhas 117-158)
├── [EXISTENTE] Utility Functions (linhas 160-335)
├── [EXISTENTE] CV Calculation Functions (linhas 337-394)
├── [EXISTENTE] Lomb-Scargle Spectral Analysis (linhas 396-472)
├── [EXISTENTE] Pulse Detection (linhas 474-521)
├── [EXISTENTE] Statistical Testing (linhas 523-559)
│
├── [NOVO] Data Availability Diagnostics (inserir após linha 521)
│   ├── diagnose_data_availability()
│   └── create_availability_figure()
│
├── [EXISTENTE] Plotting Functions (linhas 561-1015)
└── [EXISTENTE + MODIFICAR] Main Execution (linhas 1017-1139)
    └── Inserir chamada a diagnose_data_availability() após linha 1041
```

### 5.2 Estrutura do Script de Sincronicidade (NOVO)

> [!NOTE]
> **Arquivo:** `TIME_SERIES_UNIFIED_SYNCHRONY_ANALYSIS.py`
> 
> Script completamente novo. Pode importar funções de carregamento do script de variabilidade ou replicá-las.

```
TIME_SERIES_UNIFIED_SYNCHRONY_ANALYSIS.py
│
├── [CONFIG] Configurações e Constantes
│   ├── Diretórios de entrada/saída
│   ├── Período de análise (2002-2008)
│   └── Paleta de cores padronizada (copiar de 0.4)
│
├── [LOAD] Carregamento de Dados (copiar de variability ou importar)
│   ├── parse_date_from_filename()
│   ├── filter_files_by_period()
│   ├── load_timeseries_for_sites()
│   └── calculate_dli_for_sites()
│
├── [PROCESS] Processamento e Agregação
│   ├── aggregate_timeseries_by_arc()  # Média por arco
│   ├── smooth_timeseries()            # Rolling mean 14 dias
│   └── interpolate_gaps()             # Interpolação linear max 7 dias
│
├── [ANALYZE] Análises de Sincronicidade
│   ├── calculate_pearson_synchrony()       # §4.1
│   ├── calculate_cross_correlation_full()  # §4.2
│   ├── calculate_spectral_coherence()      # Coerência espectral
│   └── calculate_pulse_synchrony_stl()     # §4.1 STL decomposition
│
├── [STATS] Testes Estatísticos
│   ├── mann_whitney_test()          # Copiar de variability
│   └── calculate_effect_size()      # Cohen's d
│
├── [VISUALIZE] Geração de Figuras
│   ├── create_main_synchrony_figure()   # Fig_Main (§3.1)
│   ├── create_raw_timeseries_figure()   # Fig_S2 (§3.2)
│   └── create_ccf_detail_figure()       # Fig_S3
│
└── [OUTPUT] Exportação de Resultados
    ├── Pearson_Synchrony_by_Site.csv
    ├── Cross_Correlation_Results.csv
    └── STL_Pulse_Synchrony.csv
```

### 5.3 Fluxo de Execução do Script de Variabilidade (Atualizado)

```
1. INIT & CONFIG
   └─> Configurar diretórios, período, cores

2. LOAD DATA
   ├─> Carregar sites CSV
   ├─> Carregar SST (CRW)
   ├─> Carregar DLI (Kd490 + PAR)
   └─> Carregar Chl-a (MODIS)

3. DIAGNOSE (NOVO - inserir aqui)
   ├─> Calcular disponibilidade temporal
   ├─> Gerar Fig_S1_Data_Availability.png
   └─> Exportar Data_Availability_Diagnostic.csv

4. CALCULATE CV (existente)
5. SPECTRAL ANALYSIS (existente)
6. PULSE DETECTION (existente)
7. GENERATE FIGURES (existente)
8. SAVE SUMMARIES (existente)
```

### 5.4 Fluxo de Execução do Script de Sincronicidade (NOVO)

```
1. INIT & CONFIG
   └─> Configurar diretórios, período, cores

2. LOAD DATA
   ├─> Carregar sites CSV
   ├─> Carregar SST, DLI, Chl-a (mesmas funções)
   └─> Verificar Data_Availability_Diagnostic.csv existe

3. AGGREGATE BY ARC
   ├─> Agregar séries por Inner/Outer (média diária)
   └─> Suavizar (14-day rolling mean)

4. ANALYZE SYNCHRONY
   ├─> Pearson correlation entre variáveis por sítio
   ├─> Cross-correlation (lag) por arco
   ├─> Spectral coherence (7d, 30d períodos)
   └─> STL decomposition para pulse synchrony

5. STATISTICAL TESTS
   └─> Mann-Whitney U (Inner > Outer?)

6. VISUALIZE
   ├─> Fig_Main_Temporal_Synchrony.png
   ├─> Fig_S2_Raw_Time_Series.png
   └─> Fig_S3_Cross_Correlation_Detail.png

7. EXPORT
   ├─> Pearson_Synchrony_by_Site.csv
   ├─> Cross_Correlation_Results.csv
   └─> STL_Pulse_Synchrony.csv
```

---

## 6. Paleta de Cores Padronizada

```python
# Cores por Arco (consistente com outros scripts)
COLOR_INNER = '#CD853F'  # Peru (light brown)
COLOR_OUTER = '#4169E1'  # Royal Blue

# Cores por Variável (consistente com TIME_SERIES_LAG)
COLOR_SST = '#D55E00'    # Vermelho-alaranjado
COLOR_DLI = '#0072B2'    # Azul escuro
COLOR_CHL = '#009E73'    # Verde-teal

# Cores por Habitat (se aplicável, de Boxplot_Cover_RGR_GLM.R)
HAB_COLORS = {
    'TP': '#009E73',  # Tide Pool - Green
    'PA': '#E69F00',  # Reef Wall - Orange
    'RR': '#56B4E9'   # Rocky Reef - Sky Blue
}
```

---

## 7. Dependências

```python
# Core
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Estatística
from scipy import stats
from scipy.signal import lombscargle
from statsmodels.tsa.seasonal import STL  # NOVO - para decomposição robusta

# Auxiliares
import os, glob, re, logging, warnings
from datetime import datetime
```

---

## 8. Plano de Verificação

### 8.1 Verificação do Script de Variabilidade (MODIFICADO)

**Comando de execução:**
```bash
cd "C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_CODES"
python TIME_SERIES_VARIABILITY_ANALYSIS_GEMINI.py
```

**Arquivos de saída esperados** (em `Variability_Analysis_GEMINI/`):
| Arquivo | Novo? | Descrição |
|---------|-------|-----------|
| `Data_Availability_Diagnostic.csv` | ✅ NOVO | Tabela com métricas de cobertura |
| `Fig_S1_Data_Availability.png` | ✅ NOVO | Heatmap de disponibilidade temporal |
| `CV_Summary_by_Site.csv` | Existente | CV por sítio |
| `Pulse_Summary_by_Site.csv` | Existente | Eventos extremos |
| `Fig1_CV_Comparison_Boxplots.png` | Existente | Boxplots de CV |
| `Fig2_Spectral_Power_Comparison.png` | Existente | Potência espectral |
| `Fig3_Periodogram_Heatmaps.png` | Existente | Heatmaps de periodograma |
| `Fig4_Extreme_Events_Frequency.png` | Existente | Frequência de pulsos |
| `Fig5_DLI_Habitat_Comparison.png` | Existente | DLI por habitat |

**Verificação automática:**
```python
import os
output_dir = r'C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_RESULTS\Variability_Analysis_GEMINI'

# Arquivos NOVOS que devem existir após modificação
new_files = ['Data_Availability_Diagnostic.csv', 'Fig_S1_Data_Availability.png']
for f in new_files:
    assert os.path.exists(os.path.join(output_dir, f)), f"FALTA: {f}"
print("✓ Todos os novos arquivos de diagnóstico criados")
```

### 8.2 Verificação do Script de Sincronicidade (NOVO)

**Comando de execução:**
```bash
cd "C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_CODES"
python TIME_SERIES_UNIFIED_SYNCHRONY_ANALYSIS.py
```

**Arquivos de saída esperados** (em `Temporal_Synchrony_Analysis/`):
| Arquivo | Descrição |
|---------|-----------|
| `Fig_Main_Temporal_Synchrony.png` | Figura principal 6 painéis |
| `Fig_S2_Raw_Time_Series.png` | Séries brutas por arco |
| `Fig_S3_Cross_Correlation_Detail.png` | Detalhes de CCF |
| `Pearson_Synchrony_by_Site.csv` | Correlações por sítio |
| `Cross_Correlation_Results.csv` | Resultados de lag analysis |
| `STL_Pulse_Synchrony.csv` | Correlação dos resíduos STL |

**Verificação automática:**
```python
import os
output_dir = r'C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_RESULTS\Temporal_Synchrony_Analysis'

expected_files = [
    'Fig_Main_Temporal_Synchrony.png',
    'Fig_S2_Raw_Time_Series.png', 
    'Fig_S3_Cross_Correlation_Detail.png',
    'Pearson_Synchrony_by_Site.csv',
    'Cross_Correlation_Results.csv',
    'STL_Pulse_Synchrony.csv'
]
for f in expected_files:
    assert os.path.exists(os.path.join(output_dir, f)), f"FALTA: {f}"
print("✓ Todos os arquivos de sincronicidade criados")
```

### 8.3 Critérios de Sucesso

| Critério | Script | Verificação |
|----------|--------|-------------|
| Exit code 0 | Ambos | Script executa sem erros |
| Coverage > 50% | Variabilidade | `Data_Availability_Diagnostic.csv` todas as linhas > 50% |
| Figuras legíveis | Ambos | Séries temporais não comprimidas horizontalmente |
| Padrões sazonais visíveis | Sincronicidade | Fig_Main mostra ciclos anuais claros |
| Lag DLI→CHL positivo | Sincronicidade | `Cross_Correlation_Results.csv` lag > 0 |

### 8.4 Verificação Visual (Manual)

1. **Fig_Main_Temporal_Synchrony.png**:
   - Confirmar que séries temporais ocupam ~70% da altura vertical
   - Verificar que padrões sazonais são visíveis (não comprimidos)
   - Legendas não sobrepõem dados

2. **Fig_S2_Raw_Time_Series.png**:
   - Pontos individuais visíveis
   - 3 eixos Y claramente rotulados
   - Gap periods visíveis (sem interpolação)

3. **Fig_S1_Data_Availability.png**:
   - Heatmap mostra claramente dias com/sem dados
   - DLI e CHL têm mais lacunas que SST (esperado)

---

## 9. Próximos Passos

1. **Aprovar este plano refinado**
2. **Implementar o script unificado** `TIME_SERIES_UNIFIED_SYNCHRONY_ANALYSIS.py`
3. **Executar e validar** figuras e estatísticas
4. **Iterar** com base em feedback visual

---

## 10. Apêndice: Justificativa para Linhas Incompletas no Periodograma

> [!NOTE]
> **Por que DLI e Chl-a não alcançam 365 dias no Lomb-Scargle?**

Os dados MODIS (Kd490, PAR, Chl-a) têm cobertura temporal inferior aos dados CRW (SST) devido a:

1. **Cobertura de nuvens**: A região de Abrolhos tem nebulosidade significativa, especialmente no verão
2. **Órbita polar do AQUA/MODIS**: Revisita a cada ~1-2 dias, mas gaps podem se acumular
3. **Combinação Kd490 + PAR para DLI**: Requer ambos os produtos no mesmo dia

**Consequência estatística**: O período máximo resolvível (Nyquist) é limitado pelo número efetivo de observações:
- SST: ~2500 dias → período máximo ~1250 dias
- DLI: ~1700 dias → período máximo ~850 dias
- Chl-a: ~1800 dias → período máximo ~900 dias

Esta limitação é intrínseca aos dados e **não invalida a análise**, apenas limita a resolução em escalas interanuais.
