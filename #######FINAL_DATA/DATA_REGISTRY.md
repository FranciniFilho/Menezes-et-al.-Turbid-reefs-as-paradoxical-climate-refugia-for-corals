# Data Registry — Centralized Input Data

> **Last updated**: 2026-09-12  
> **Purpose**: Manifest of all input data files centralized in `#######FINAL_DATA/`.  
> **Status**: path migration COMPLETE — all scripts read from `#######FINAL_DATA/` and write to `#######FINAL_RESULTS/`.
> **Local-only datasets** (archived outside the public repo): Diving-PAM 2013 campaign, BENTHOS_TEMPORAL_CLEAN_v2, sites_list_full_clean_PLEST, GIS baselayers.

---

## Directory Structure Overview

```
#######FINAL_DATA/
├── 00_sites_metadata/          # Site reference & coordinate files
├── 01_raw_biological/          # Raw monitoring & measurement data
├── 02_PCA_environmental/       # PCA scores from environmental analysis
│   ├── CV_02/                  # Short-term CV window (2 months)
│   ├── CV_30/                  # Seasonal CV window (30 months)
│   └── CV_ALL/                 # Full time-series CV window
├── 03_modeling_data/           # Integrated datasets ready for Bayesian modeling
│   ├── CV_02/
│   ├── CV_30/
│   └── CV_ALL/
├── 04_health_growth_PCA/       # Health & growth PCA outputs (colony-level)
└── 05_benthic_cover/           # Benthic community cover data (long format)
```

---

## 00_sites_metadata/

| File | Source Path | Description | Used By |
|------|-----------|-------------|---------|
| `sites_list_full.csv` | `########CEBIMAR\####PROJETOS\#####Coral trade offs\` | Full site list with coordinates (lat, lon), reef names, SITE codes, ARCH designation | PCA_LOCAL_bubbleplot.py, FINAL_DATA_INTEGRATION.py, multiple R model scripts |
| `sites_list_full_clean.csv` | `########CEBIMAR\####PROJETOS\#####Coral trade offs\` | Cleaned version with standardized site names | PCA_LOCAL_bubbleplot_v2.py, some time-series scripts |

---

## 01_raw_biological/

| File | Source Path | Description | Used By |
|------|-----------|-------------|---------|
| `Vitality_and_size_new.xlsx` | Already in FINAL_DATA (root, as `######Vitality and size_new.xlsx`) | Excel version of vitality data (original with formatting) | Manual reference |

---

## 02_PCA_environmental/

Each CV subdirectory contains PCA outputs from `PCA_LOCAL_bubbleplot.py` / `PCA_LOCAL_bubbleplot_v2.py`.

### CV_02/ (Short-term variability window — 2 months)

| File | Description | Used By |
|------|-------------|---------|
| `dados_abundancia_integrados_long_format.csv` | Integrated abundance + PCA scores (long format); **primary input for ZOIB & JSDM models** | 01_ZOIB, 03_JSDM scripts |
| `dados_consolidados_com_scores_das_duas_PCAs_CV2.csv` | Integrated abundance + both PCA score sets (CV_02) | FINAL_DATA_INTEGRATION.py, modeling scripts |
| `loadings_PCA_Magnitude_CV_2.csv` | PCA loading matrix (magnitude) | PCA visualization scripts |
| `loadings_PCA_Variability_CV_2.csv` | PCA loading matrix (variability) | PCA visualization scripts |
| `dados_abundancia_integrados_long_format.csv` | **Primary ZOIB & JSDM input** (long format) | 01_ZOIB, 03_JSDM |
| `resumo_segmentos_por_site.csv/.xlsx` | Per-site segment summary | Documentation |

### CV_30/ (Seasonal variability window — 30 months)
Same structure as CV_02 with `_CV30` suffixes (`dados_consolidados_com_scores_das_duas_PCAs_CV30.csv`, `loadings_PCA_*_CV_30.csv`).

### CV_ALL/ (Full time-series variability window)
Same structure as CV_02 with `_CVall` suffixes (`dados_consolidados_com_scores_das_duas_PCAs_CVall.csv`, `loadings_PCA_*_CV_all.csv`).

**Key insight**: CV_02 captures short-term stability vs. instability (PC1 dominates 76.1%). CV_30 captures seasonal variability (more balanced PC1/PC2). CV_ALL captures inter-annual trends with regional homogenization.

---

## 03_modeling_data/

Each CV subdirectory contains the final integrated dataset for Bayesian modeling.

### CV_02/, CV_30/, CV_ALL/

| File | Source Dir | Description | Used By |
|------|-----------|-------------|---------|
| `dados_finais_para_modelagem_com_ARCH.csv` | `output_DADOS_FINAIS_PARA_MODELAGEM_cv_{2,30,all}/` | Complete dataset with biological + environmental + PCA + ARCH, ready for modeling. Contains: SITE, REEF, YEAR, ORGANISMO, COBERTURA, HAB, ARCH, DEPTH_M, PC1_MAGNITUDE, PC2_MAGNITUDE, PC1_VARIABILITY, PC2_VARIABILITY, HEALTH_PC1, HEALTH_PC2, RGR, LAT, LON | 02a_Gaussian_RGR, 02b_Gaussian_Health scripts |

---

## 04_health_growth_PCA/

| File | Source Dir | Description | Used By |
|------|-----------|-------------|---------|
| `resultados_biologicos_por_colonia.csv` | `#output_ANALISE_BIOLOGICA_boxplot_PCAnew_YEAR_RE/` | Colony-level health PCA scores (HEALTH_PC1, HEALTH_PC2) and RGR | 02a_Gaussian_RGR, 02b_Gaussian_Health |
| `PCA_scores_health.csv` | Same | Health PCA scores per colony | Health model data preparation |
| `PCA_loadings_health.csv` | Same | Health PCA loading matrix (ALIVE, DEAD, BLEACHED → PC1, PC2) | Documentation |
| `PCA_explained_variance_health.csv` | Same | Explained variance for health PCA | Documentation |
| `dados_saude_crescimento_integrados.csv` | Same | Integrated health + growth dataset | Health model scripts |
| Other CSV/XLSX files | Same | Additional health analysis intermediates | Various |

---

## 07_lomb_frequencies/
| File | Description | Used By |
|------|-------------|---------|
| `frequencias_primarias_secundarias_.csv` | Primary/secondary Lomb-Scargle spectral frequencies per reef arc (DHW exposure inputs) | FINAL_DATA_INTEGRATION.py |

## 05_benthic_cover/

| File | Source Dir | Description | Used By |
|------|-----------|-------------|---------|
| `dados_integrados_long_format.csv` | `%%BOX_PLOTS_benthic_cover/` | Benthic community cover in long format (all organism groups) | Boxplot_Cover_RGR_Publication.R, JSDM community matrix |

---

## Data NOT Centralized (documented only)

### Remote Sensing NetCDF (External Drives)
Too large to copy; used only by Python environmental processing scripts.

| Data | Typical Path | Used By |
|------|-------------|---------|
| NOAA CRW SST (daily) | `K:\NOAA_CRW_SST\` or `H:\MODIS_DLI\` | PCA_LOCAL_bubbleplot.py, TIME_SERIES_*.py |
| MODIS DLI (8-day) | `E:\MODIS_DLI_8DAY\` or `F:\MODIS_DLI\` | PCA_LOCAL_bubbleplot.py |
| MODIS CHL-a | `K:\MODIS_CHL\` or similar | PCA_LOCAL_bubbleplot.py |
| NOAA CRW DHW | `K:\NOAA_CRW_DHW\` | DHW_maps.py, EXTREME_EVENTS_*.py |

### Model Output RDS Files (outputs, not inputs)
| Location | Contents |
|----------|----------|
| `#######FINAL_RESULTS/BAYES_MODELS_YEAR_RE/` | Canonical YEAR_RE model fits |
| `#######FINAL_RESULTS/BAYES_MODELS_YEAR_RE_PRIOR_SENSITIVITY/` | Prior sensitivity grid |
| `#######FINAL_RESULTS\BAYES_MODELS_FINAL_RES\` | Archived model copies |
| `#######FINAL_RESULTS/JSDM_Dominance/` | JSDM dominance analysis outputs |
| `#######FINAL_RESULTS/JSDM_Palythoa_Opus/` | Palythoa JSDM outputs |

---

## Path Mapping: Current Script Paths → Centralized Locations

This mapping was APPLIED to all non-legacy scripts on 2026-09-12 (commit: 'refactor(paths): replace OneDrive/ONEDRIVE_NOVO absolute paths').

| Current Script Path | Centralized Path |
|---------------------|-----------------|
| `../#######FINAL_RESULTS/#####output_local_PCA_CV_2_FINAL/*.csv` | `#######FINAL_DATA/02_PCA_environmental/CV_02/*.csv` |
| `../#######FINAL_RESULTS/#####output_local_PCA_CV_30_FINAL/*.csv` | `#######FINAL_DATA/02_PCA_environmental/CV_30/*.csv` |
| `../#######FINAL_RESULTS/#####output_local_PCA_CV_all_FINAL/*.csv` | `#######FINAL_DATA/02_PCA_environmental/CV_ALL/*.csv` |
| `../#######FINAL_RESULTS/output_DADOS_FINAIS_PARA_MODELAGEM_cv_2/*.csv` | `#######FINAL_DATA/03_modeling_data/CV_02/*.csv` |
| `../#######FINAL_RESULTS/output_DADOS_FINAIS_PARA_MODELAGEM_cv_30/*.csv` | `#######FINAL_DATA/03_modeling_data/CV_30/*.csv` |
| `../#######FINAL_RESULTS/output_DADOS_FINAIS_PARA_MODELAGEM_cv_all/*.csv` | `#######FINAL_DATA/03_modeling_data/CV_ALL/*.csv` |
| `../#######FINAL_RESULTS/#output_ANALISE_BIOLOGICA_boxplot_PCAnew_YEAR_RE/*.csv` | `#######FINAL_DATA/04_health_growth_PCA/*.csv` |
| `../#######FINAL_RESULTS/%%BOX_PLOTS_benthic_cover/*.csv` | `#######FINAL_DATA/05_benthic_cover/*.csv` |
| `########CEBIMAR\...\sites_list_full.csv` | `#######FINAL_DATA/00_sites_metadata/sites_list_full.csv` |
| `########CEBIMAR\...\sites_list_full_clean.csv` | `#######FINAL_DATA/00_sites_metadata/sites_list_full_clean.csv` |
