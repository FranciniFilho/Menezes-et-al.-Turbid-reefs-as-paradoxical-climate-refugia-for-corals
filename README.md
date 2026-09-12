# Thriving under stress: How turbid coastal reefs act as paradoxical climate refugia for corals

**Rafael Menezes<sup>1,\*</sup>, Carolina F. Schlosser<sup>2,3,\*</sup>, Mariana Marcondes Couto<sup>2,3</sup>, Ericka O.C. Coni<sup>2</sup>, Camilo M. Ferreira<sup>4</sup>, Pedro M. Meirelles<sup>4</sup>, Andreia C.C. Barbosa<sup>2</sup>, Marcelo V. Kitahara<sup>2</sup>, Miguel Mies<sup>5</sup>, Samuel C. Faria<sup>2</sup>, Ronaldo B. Francini-Filho<sup>2,\*\*</sup>**

<sup>1</sup>Universidade Federal Rural de Pernambuco, Recife, PE, Brazil · <sup>2</sup>Centro de Biologia Marinha, Universidade de São Paulo, São Sebastião, SP, Brazil · <sup>3</sup>Programa de Pós-graduação em Ecologia, USP · <sup>4</sup>Instituto de Biologia, Universidade Federal da Bahia, Salvador, BA, Brazil · <sup>5</sup>Instituto Oceanográfico, USP

\* Equal contribution · \*\* Corresponding author: francinifilho@usp.br

## Abstract

Naturally stressful reef environments can enhance coral thermal tolerance, yet whether chronic stress exposure improves or erodes baseline coral performance remains unresolved. Combining three years of in situ monitoring with satellite-derived environmental data across Brazil's Abrolhos Bank, we employed Bayesian hierarchical models and Joint Species Distribution Models (JSDMs) to assess the health, growth, and competitive dynamics of the endemic reef-building coral *Mussismilia hispida*. Coral metrics were modeled against the magnitude and multiscale temporal variability of SST, chlorophyll-a, and Daily Light Integral (DLI), alongside the frequency of acute thermal stress (DHW > 4). Reef walls of turbid coastal reefs emerged as paradoxical refugia: despite experiencing the highest environmental instability and 15-fold more frequent acute thermal stress, they supported the healthiest coral populations, underpinned by lower DLI, cooler baseline SST, and elevated Chl-a. At the community level, JSDMs identified DLI as the primary driver of benthic composition, with macroalgal dominance over *M. hispida* increasing sharply under higher irradiance. We propose that dynamic coastal environments protect corals through a dual mechanism, functioning simultaneously as acclimatization arenas and as ecological filters that suppress competitors while promoting facilitating species.

**Keywords:** Coral Bleaching · Environmental Memory · Adaptive Potential · *Mussismilia hispida* · Reef Conservation

## Repository Structure

```
#######FINAL_CODES/                                   # Analysis scripts (R + Python)
  Bayes_models/New_Bayes_Models_YEAR_RE/              # PRIMARY Bayesian pipeline (210 models)
  Bayes_models/00_LOO_Selection_Functions.R           # Shared LOO/convergence functions (sourced by the pipeline)
#######FINAL_DATA/                                    # Curated input data (tracked)
  00_sites_metadata/                                  # Site coordinates & reference lists
  01_raw_biological/                                  # Raw vitality measurements
  02_PCA_environmental/CV_{02,30,ALL}/                # Environmental PCA scores & loadings
  03_modeling_data/CV_{02,30,ALL}/                    # Integrated modeling datasets
  04_health_growth_PCA/                               # Colony-level health & growth PCA
  05_benthic_cover/                                   # Benthic community cover (long format)
  07_lomb_frequencies/                                # Lomb-Scargle frequencies (DHW exposure inputs)
  DATA_REGISTRY.md                                    # Complete data manifest
docs/                                                 # JSDM correction guide & Dirichlet reference
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
| 5 | R | `Bayes_models/New_Bayes_Models_YEAR_RE/00_RUN_ALL_YEAR_RE_MODELS.R` | Bayesian modeling (210 models) |
| 6 | R | `04_MASTER_Viz_Pipeline_YEAR_RE_v7_CATEGORICAL.R` | Publication figures |
| 7 | R | `07/08_Extract|Visualize_JSDM_Dominance*.R` | JSDM dominance analysis |

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

The environmental PCA (step 1) consumes the raw satellite archive, resolved through an environment variable:

| Variable | Used by | Content |
|----------|---------|---------|
| `MUSHIS_RS_ARCHIVE` | `PCA_LOCAL_bubbleplot_v2.py` | Root folder containing `NOAA_CRW_SST/`, `MODIS_DLI_8DAY/`, `MODIS_CHL/`, `NOAA_CRW_DHW/` (NOAA CRW & MODIS Aqua archives) |

Everything needed for the statistical analyses (steps 2–7) is embedded in `#######FINAL_DATA/`.

## Data Availability

All curated input data are tracked in this repository (`#######FINAL_DATA/`, 42 files; see `DATA_REGISTRY.md`). Model fits (3.4 GB), high-resolution figure outputs, and the 2013 Diving-PAM field photographs are archived on Zenodo: **[DOI_PENDING]**. `#######FINAL_RESULTS/` is generated locally at runtime and intentionally not tracked.

## Documentation

- **[DATA_REGISTRY.md](#######FINAL_DATA/DATA_REGISTRY.md)** — complete manifest of input data
- **[docs/](docs/)** — JSDM correction guide and Dirichlet-model technical reference
