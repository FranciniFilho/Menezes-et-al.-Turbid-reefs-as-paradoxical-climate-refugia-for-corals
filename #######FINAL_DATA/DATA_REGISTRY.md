# Data Registry — Centralized Input Data

> **Last updated**: 2025-07  
> **Purpose**: Manifest of all input data files centralized in `#######FINAL_DATA/`.  
> Scripts still reference original paths — a future migration will update them.

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
| `sites_list_full_clean_PLEST.csv` | `######22.04.23\DATA\` | Site list including Palythoa/Millepora estimates for JSDM with Palythoa | 00_Palythoa_JSDM_Functions.R |

---

## 01_raw_biological/

| File | Source Path | Description | Used By |
|------|-----------|-------------|---------|
| `Vitality and size_new.csv` | Already in FINAL_DATA (root) | Colony-level vitality, size, and health measurements for M. hispida | Calculate_RGR_&_health_PCA.py |
| `Vitality_and_size_new.xlsx` | Already in FINAL_DATA (root, as `######Vitality and size_new.xlsx`) | Excel version of vitality data (original with formatting) | Manual reference |
| `BENTHOS_TEMPORAL_CLEAN_v2.csv` | `######22.04.23\DATA\BENTHOS TEMPORAL CLEAN_v2.csv` | Multi-year benthic community composition transect data (all species) | 00_Palythoa_JSDM_Functions.R, JSDM data preparation |
| `DivingPAM_2013_TIM3/` | Already in FINAL_DATA | Diving-PAM fluorometry data (Fv/Fm) from 2013 campaign at Timbebas reef | PAM_Analysis_Complete.R |

---

## 02_PCA_environmental/

Each CV subdirectory contains PCA outputs from `PCA_LOCAL_bubbleplot.py` / `PCA_LOCAL_bubbleplot_v2.py`.

### CV_02/ (Short-term variability window — 2 months)

| File | Description | Used By |
|------|-------------|---------|
| `dados_abundancia_integrados_long_format.csv` | Integrated abundance + PCA scores (long format); **primary input for ZOIB & JSDM models** | 01_ZOIB, 03_JSDM scripts |
| `PCA_scores_Magnitude_CV_02.csv` | PCA scores for magnitude (mean SST, DLI, CHL) | FINAL_DATA_INTEGRATION.py |
| `PCA_scores_Variability_CV_02.csv` | PCA scores for variability (CV of SST, DLI, CHL) | FINAL_DATA_INTEGRATION.py |
| `PCA_loadings_Magnitude.csv` | PCA loading matrix (magnitude component) | PCA visualization scripts |
| `PCA_loadings_Variability_CV_02.csv` | PCA loading matrix (variability component) | PCA visualization scripts |
| `PCA_explained_variance_Magnitude.csv` | Explained variance per component (magnitude) | Documentation & figures |
| `PCA_explained_variance_Variability_CV_02.csv` | Explained variance per component (variability) | Documentation & figures |

### CV_30/ (Seasonal variability window — 30 months)
Same file structure as CV_02, with `_CV_30` suffixes.

### CV_ALL/ (Full time-series variability window)
Same file structure as CV_02, with `_CV_all` suffixes.

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

### GIS Shapefiles
| Data | Path | Used By |
|------|------|---------|
| Bathymetry | `C:\Users\rbfra\OneDrive\GIS shapes\` | PCA_GLOBAL_maps_COMPOSITE_GEMINI.py |
| Coastline | Same | Map scripts |
| Reef polygons | Same | Map scripts |

### Model Output RDS Files (outputs, not inputs)
| Location | Contents |
|----------|----------|
| `C:\Users\rbfra\OneDrive\New_Bayes_Models_Output\` | Canonical YEAR_RE model fits |
| `C:\Users\rbfra\OneDrive\New_Bayes_Models_Output_PRIOR_SENSITIVITY_FULLGRID_v1\` | Prior sensitivity grid |
| `#######FINAL_RESULTS\BAYES_MODELS_FINAL_RES\` | Archived model copies |
| `C:\Users\rbfra\OneDrive\JSDM_Dominance_Figures\` | JSDM dominance analysis outputs |
| `C:\Users\rbfra\OneDrive\JSDM_Palythoa_Output_Opus\` | Palythoa JSDM outputs |

---

## Path Mapping: Current Script Paths → Centralized Locations

This table maps hardcoded paths in scripts to their new centralized locations. **Scripts have NOT yet been updated** — this mapping is for a future migration task.

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
| `######22.04.23\DATA\sites_list_full_clean_PLEST.csv` | `#######FINAL_DATA/00_sites_metadata/sites_list_full_clean_PLEST.csv` |
| `######22.04.23\DATA\BENTHOS TEMPORAL CLEAN_v2.csv` | `#######FINAL_DATA/01_raw_biological/BENTHOS_TEMPORAL_CLEAN_v2.csv` |
