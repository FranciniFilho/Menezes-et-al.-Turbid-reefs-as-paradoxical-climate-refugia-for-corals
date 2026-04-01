# Thriving under stress: How turbid coastal reefs act as paradoxical climate refugia for corals

> **Authors**: Schlosser, Menezes, Duarte, Feitosa, Ferreira, Kelmo, Zilberberg, Francini-Filho
> **Status**: Under revision (Round 3)

## Overview

Marine ecology research investigating the distribution, abundance, health, and growth of the endemic reef-building coral ***Mussismilia hispida*** across the Abrolhos Bank, Brazil. Uses Bayesian models (ZOIB, Gaussian, JSDM Dirichlet) via `brms` to quantify coral responses to environmental gradients (SST, DLI, CHL-a, depth) across habitats and reef arcs.

## Central Hypothesis — The "Hidden Memory" Refugia

Turbid, nearshore inner-arc reefs act as **paradoxical climate refugia**: conditions historically considered stressful (high turbidity, lower light) actually buffer corals against thermal stress. Lower SST variability, reduced light, and higher chlorophyll promote superior coral conditions compared to "pristine" outer-arc reefs.

## Repository Structure

```
#######FINAL_CODES/                              # All analysis scripts (R + Python)
  Bayes_models/
    New_Bayes_Models_YEAR_RE/                    # Primary Bayesian pipeline (210 models)
    New_Bayes_Models_YEAR_RE_GLM/                # GLM variants (experimental)
    New_Bayes_Models_YEAR_RE_OPUS/               # Palythoa-augmented JSDM
    BACKUP_legacy_scripts_before_YEAR_RE/        # Archived legacy scripts
#######FINAL_DATA/                               # Centralized input data
  00_sites_metadata/                             # Site coordinates & reference files
  01_raw_biological/                             # Raw monitoring & measurement data
  02_PCA_environmental/                          # PCA scores (CV_02, CV_30, CV_ALL)
  03_modeling_data/                              # Integrated modeling datasets
  04_health_growth_PCA/                          # Colony-level health & growth PCA
  05_benthic_cover/                              # Benthic community cover data
  DATA_REGISTRY.md                               # Complete manifest of all data files
#######FINAL_RESULTS/                            # Generated outputs (figures, models)
#######LAST ROUND MS/                            # Scientific manuscript
#######PHOTOS AND MAPS/                          # Visualizations and figures
docs/                                            # Technical notes, diagnostics, quickstart
```

## Analysis Pipeline

| Step | Tool | Description |
|------|------|-------------|
| 1 | Python | Environmental PCA (3 CV windows: CV_02, CV_30, CV_ALL) |
| 2 | Python | PERMANOVA — REEF factor redundancy test |
| 3 | Python | Health PCA + Relative Growth Rate (RGR) |
| 4 | Python | Data integration (biological + environmental + metadata) |
| 5 | Python | Temporal trends & extreme events analysis |
| 6 | Python | Global environmental maps & DHW maps |
| 7 | R | Bayesian modeling: ZOIB (abundance), Gaussian (health/growth), JSDM Dirichlet (community) |
| 8 | R | Publication figures (forest plots, PDPs, variance partitioning) |
| 9 | R | JSDM dominance probability analysis |
| 10 | R | Complementary BRT models |

## Model Architecture

| Model | Response | Family | Random Effects |
|-------|----------|--------|----------------|
| ZOIB | Coral cover (0-1) | `zero_one_inflated_beta()` | `(1|SITE) + (1|YEAR)` |
| Health | HEALTH_PC1, PC2 | `gaussian()` | `(1|SITE) + (1|YEAR)` |
| RGR | Growth rate | `gaussian()` | `(1|REEF)` |
| JSDM | Community composition | `dirichlet()` | `(1|SITE) + (1|YEAR)` |

7 candidate models per response, compared via LOOIC. Global winners selected across 3 CV scenarios.

## Key Variables

- **HAB**: PA (reef walls), RR (rocky reefs), TP (reef tops)
- **ARCH** (ARC): Inner (coastal refugia) vs Outer (oceanic, stressed)
- **PCA Magnitude**: Mean environmental conditions (SST, DLI, CHL-a)
- **PCA Variability**: Temporal CV at 3 windows (short-term, seasonal, inter-annual)

## Quick Start

```bash
# Full Bayesian pipeline
Rscript "#######FINAL_CODES/Bayes_models/New_Bayes_Models_YEAR_RE/00_RUN_ALL_YEAR_RE_MODELS.R" canonical WeaklyInformative

# Publication figures
Rscript "#######FINAL_CODES/Bayes_models/New_Bayes_Models_YEAR_RE/04_MASTER_Viz_Pipeline_YEAR_RE_v7_CATEGORICAL.R"

# Environmental PCA (Python)
python PCA_LOCAL_bubbleplot_v2.py
```

## Documentation

- **[AGENTS.md](AGENTS.md)** — Repository guidelines, pipeline summary, coding conventions
- **[CLAUDE.md](CLAUDE.md)** — Comprehensive project documentation (527 lines)
- **[DATA_REGISTRY.md](#######FINAL_DATA/DATA_REGISTRY.md)** — Complete manifest of centralized input data
- **[docs/](docs/)** — Technical notes, JSDM diagnostics, quickstart guides

## Dependencies

**R**: brms, cmdstanr, tidybayes, ggdist, ggplot2, patchwork, sf, dplyr
**Python**: pandas, numpy, xarray, scikit-learn, scipy, matplotlib, cartopy, geopandas
