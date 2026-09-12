# Thriving under stress: How turbid coastal reefs act as paradoxical climate refugia for corals

**Rafael Menezes<sup>1,\*</sup>, Carolina F. Schlosser<sup>2,3,\*</sup>, Mariana Marcondes Couto<sup>2,3</sup>, Ericka O.C. Coni<sup>2</sup>, Camilo M. Ferreira<sup>4</sup>, Pedro M. Meirelles<sup>4</sup>, Andreia C.C. Barbosa<sup>2</sup>, Marcelo V. Kitahara<sup>2</sup>, Miguel Mies<sup>5</sup>, Samuel C. Faria<sup>2</sup>, Ronaldo B. Francini-Filho<sup>2,\*\*</sup>**

<sup>1</sup>Universidade Federal Rural de Pernambuco, Recife, PE, Brazil · <sup>2</sup>Centro de Biologia Marinha, Universidade de São Paulo, São Sebastião, SP, Brazil · <sup>3</sup>Programa de Pós-graduação em Ecologia, USP · <sup>4</sup>Instituto de Biologia, Universidade Federal da Bahia, Salvador, BA, Brazil · <sup>5</sup>Instituto Oceanográfico, USP

\* Equal contribution · \*\* Corresponding author: francinifilho@usp.br

## Abstract

Naturally stressful reef environments can enhance coral thermal tolerance, yet whether chronic stress exposure improves or erodes baseline coral performance remains unresolved. Combining three years of in situ monitoring with satellite-derived environmental data across Brazil's Abrolhos Bank, we employed Bayesian hierarchical models and Joint Species Distribution Models (JSDMs) to assess the health, growth, and competitive dynamics of the endemic reef-building coral *Mussismilia hispida*. Coral metrics were modeled against the magnitude and multiscale temporal variability of SST, chlorophyll-a, and Daily Light Integral (DLI), alongside the frequency of acute thermal stress (DHW > 4). Reef walls of turbid coastal reefs emerged as paradoxical refugia: despite experiencing the highest environmental instability and 15-fold more frequent acute thermal stress, they supported the healthiest coral populations, underpinned by lower DLI, cooler baseline SST, and elevated Chl-a. At the community level, JSDMs identified DLI as the primary driver of benthic composition, with macroalgal dominance over *M. hispida* increasing sharply under higher irradiance. We propose that dynamic coastal environments protect corals through a dual mechanism, functioning simultaneously as acclimatization arenas and as ecological filters that suppress competitors while promoting facilitating species.

**Keywords:** Coral Bleaching · Environmental Memory · Adaptive Potential · *Mussismilia hispida* · Reef Conservation

## Repository Structure

```
#######FINAL_CODES/                                   # All analysis scripts (R + Python)
  Bayes_models/New_Bayes_Models_YEAR_RE/              # PRIMARY Bayesian pipeline (210 models)
  Bayes_models/New_Bayes_Models_YEAR_RE_GLM/          # GLM variants (Palythoa-augmented JSDM)
  Bayes_models/New_Bayes_Models_YEAR_RE_OPUS/         # Palythoa JSDM extension
  Variability_analysis/                               # Time-series & extreme-events scripts
  Legacy/                                             # Superseded scripts (historical reference)
#######FINAL_DATA/                                    # Curated input data (42 files, tracked)
  00_sites_metadata/                                  # Site coordinates & reference lists
  01_raw_biological/                                  # Raw monitoring & vitality measurements
  02_PCA_environmental/CV_{02,30,ALL}/                # Environmental PCA scores & loadings
  03_modeling_data/CV_{02,30,ALL}/                    # Integrated modeling datasets
  04_health_growth_PCA/                               # Colony-level health & growth PCA
  05_benthic_cover/                                   # Benthic community cover (long format)
  06_gis_shapes/                                      # GIS vector baselayers + GEBCO_2024b.nc
  DivingPAM_2013_TIM3/                                # PAM fluorometry data (.xls/.pam)
  DATA_REGISTRY.md                                    # Complete data manifest
docs/                                                 # Technical notes & historical guides
#######FINAL_RESULTS/                                 # Generated at runtime (gitignored; see Data Availability)
```

## Quick Start

Requirements: R ≥ 4.3 with `brms`, `cmdstanr`, `tidybayes`, `ggdist`, `patchwork`, `sf`; Python ≥ 3.10 with `pandas`, `numpy`, `scikit-learn`, `scipy`, `matplotlib`. All commands run **from the repository root** (the folder containing this README).

```bash
# 1. Colony health PCA + growth rates (embedded data only)
python "#######FINAL_CODES/Calculate_RGR_&_health_PCA.py"

# 2. Rebuild the integrated modeling tables (CV_ALL shown; outputs under FINAL_RESULTS)
python "#######FINAL_CODES/FINAL_DATA_INTEGRATION.py"

# 3. One ZOIB model family (R; full 210-model grid = 00_RUN_ALL_YEAR_RE_MODELS.R)
Rscript "#######FINAL_CODES/Bayes_models/New_Bayes_Models_YEAR_RE/01_ZOIB_YEAR_RE_Full_LOO_Selection.R"

# 4. Publication figures from fitted winners
Rscript "#######FINAL_CODES/Bayes_models/New_Bayes_Models_YEAR_RE/04_MASTER_Viz_Pipeline_YEAR_RE_v7_CATEGORICAL.R"
```

Every script verifies at startup that it is being run from the repository root and fails with an instructive error otherwise.

## Analysis Pipeline

| Step | Tool | Scripts | Description |
|------|------|---------|-------------|
| 1 | Python | `PCA_LOCAL_bubbleplot_v2.py` | Environmental PCA (3 CV windows; needs satellite archive — see External data) |
| 2 | Python | `PERMANOVA_Local_PCA.py` | REEF factor redundancy test |
| 3 | Python | `Calculate_RGR_&_health_PCA.py` | Health PCA + growth rates |
| 4 | Python | `FINAL_DATA_INTEGRATION.py` | Integrate biological + environmental data |
| 5 | Python | `Variability_analysis/TIME_SERIES_*.py`, `EXTREME_EVENTS_*_v2.py` | Temporal trends & extreme events |
| 6 | Python | `PCA_GLOBAL_maps_COMPOSITE_*.py`, `DHW_maps.py`, `GA_map_1x3_strip_v3.py` | Maps & spatial figures |
| 7 | R | `Bayes_models/New_Bayes_Models_YEAR_RE/00_RUN_ALL_YEAR_RE_MODELS.R` | Bayesian modeling (210 models) |
| 7-ALT | R | `Bayes_models/New_Bayes_Models_YEAR_RE_OPUS/03_JSDM_with_PALYTHOA.R` | Palythoa-augmented JSDM |
| 8 | R | `04_MASTER_Viz_Pipeline_YEAR_RE_v7_CATEGORICAL.R` | Publication figures |
| 9 | R | `07/08_Extract|Visualize_JSDM_Dominance*.R` | JSDM dominance analysis |
| 10 | R | `GAM_Size_vs_Bleaching_Mortality.R` | Size vs bleaching/mortality GAMs |

## Model Architecture

| Model | Response | Family | Random Effects |
|-------|----------|--------|----------------|
| ZOIB | Coral cover (0–1) | `zero_one_inflated_beta()` | `(1\|SITE) + (1\|YEAR)` |
| Health | HEALTH_PC1, PC2 | `gaussian()` | `(1\|SITE) + (1\|YEAR)` |
| RGR | Growth rate | `gaussian()` | `(1\|REEF)` |
| JSDM | Community composition | `dirichlet()` | `(1\|SITE) + (1\|YEAR)` |

Seven candidate models per response, compared via LOOIC; global winners selected across 3 CV scenarios (CV_02 short-term, CV_30 seasonal, CV_ALL inter-annual).

## Key Variables

- **HAB**: PA (reef walls), RR (rocky reefs), TP (reef tops); HABMERGED for JSDM
- **ARCH** (ARC): Inner (coastal refugia) vs Outer (oceanic)
- **PCA Magnitude**: mean environmental conditions (SST, DLI, Chl-a)
- **PCA Variability**: temporal CV at 3 windows

## External Data (optional steps)

Scripts that consume the raw satellite archive or very large GIS rasters resolve their location through environment variables and fail with instructions when unset:

| Variable | Used by | Content |
|----------|---------|---------|
| `MUSHIS_RS_ARCHIVE` | `PCA_LOCAL_bubbleplot*`, `TIME_SERIES_*`, `EXTREME_EVENTS_*`, `DHW_maps.py`, `PCA_GLOBAL_maps_*`, `Smart_MODIS_database_update.py` | Root folder containing `NOAA_CRW_SST/`, `MODIS_DLI_8DAY/`, `MODIS_CHL/`, `NOAA_CRW_DHW/`, `MODIS_DATA_FULL/`, `CRW_DHW_FULL/`, `CRW_SST_FULL/` (NOAA CRW & MODIS Aqua archives) |
| `MUSHIS_GIS_DIR` | `GA_map_1x3_strip_v3.py`, `Graphical_Abstract_GCB_RGB_Classification.py`, `DHW_maps.py`, `Optimal_CV_window.py`, `PCA_GLOBAL_maps_*` | Folder containing `gebco_2024_ASO.nc` (GEBCO 2024 bathymetry, https://www.gebco.net/) and `uc_fed_agosto_2016_site_shp/` (federal marine protected areas, MMA/ICMBio). Defaults to the embedded `#######FINAL_DATA/06_gis_shapes/` |
| `EARTHDATA_TOKEN` | `Smart_MODIS_database_update.py` | NASA EarthData JWT (https://urs.earthdata.nasa.gov/) for re-downloading MODIS scenes |

Everything needed for the statistical analyses (steps 2–5, 7–10) is embedded in `#######FINAL_DATA/`.

## Data Availability

All curated input data are tracked in this repository (`#######FINAL_DATA/`, 42 files; see `DATA_REGISTRY.md`). Model fits (3.4 GB), high-resolution figure outputs, and the 2013 Diving-PAM field photographs are archived on Zenodo: **[DOI_PENDING]**. `#######FINAL_RESULTS/` is generated locally at runtime and intentionally not tracked.

## Documentation

- **[DATA_REGISTRY.md](#######FINAL_DATA/DATA_REGISTRY.md)** — complete manifest of input data
- **[AGENTS.md](AGENTS.md)** / **[CLAUDE.md](CLAUDE.md)** — repository guidelines and pipeline notes
- **[docs/](docs/)** — technical notes; historical JSDM-correction guides are banner-marked
