# Analysis Code (`code/`)

This directory contains the Python analysis, validation, robustness, and figure-generation workflows. All scripts use portable relative paths aligned with `../config_paths.yaml` (see `path_config.py`) and the unified **1990–2019** historical baseline and **2070–2100** late-century evaluation window.

## Core Pipeline Scripts (`00`–`06`)

- `path_config.py` — Helper module loading `../config_paths.yaml` from the repository root.
- `00_preprocess_model_inputs.py` — Preprocesses daily CWatM water-stress, soil-moisture, and SPEI outputs.
- `01_data_loading_and_loss_metrics.py` — Loads gridded CWatM actual/counterfactual runs and computes global and system-level yield-loss metrics (baseline MIRCA-OS harvested areas).
- `02_spatial_loss_maps.py` — Generates spatial drought and waterlogging loss maps.
- `03_country_trade_and_risk_analysis.py` — Computes national Import Dependency Ratios (`IDR`) and trade-exposure risk quadrants across 161 countries.
- `04_yield_and_stress_validation.py` — Evaluates simulated crop-yield anomalies against FAOSTAT and national archives (`1990–2019`).
- `05_admin_water_stress_validation.py` — Evaluates subnational ADM1 event-year yield losses against CIL archives.
- `06_historical_yield_validation.py` — Summarizes historical yield-anomaly validation metrics (`1990–2019`).

## Revision `v4` Workflows (`v4_revision/`)

- `v4_revision/01_main_and_ed_figures/` — Rebuilds Main Manuscript Figure 4 and Extended Data outputs.
- `v4_revision/02_si_soil_moisture_maps/` — Generates Supplementary soil-moisture maps (`Figs. S11–S12`).
- `v4_revision/03_si_multiscale_validation/` — Multi-scale validation (`GDHY`, `USDA NASS/RMA`, `EM-DAT/GDIS`, `MODIS GFD`, and `11 ISIMIP3a/GGCMI3` models; `Figs. S13–S16`).
- `v4_revision/04_si_robustness_synthesis/` — Parameter, stress-curve, rice stage-gating (`SPEI_60 > +1.0`), and 8-dimension robustness synthesis (`Figs. S9, S17`).
- `v4_revision/05_tables_generation/` — Generates `Tables S1–S9` and `Extended Data Table 1`.
- `v4_revision/06_upstream_data_pipelines/` — Upstream `1990–2019` validation harmonization, `MIRCA-OS` vs. `SPAM 2020` sensitivity, and daily stress-curve sensitivity runners.
- `v4_revision/run_all_figures_and_tables.py` — One-command driver to regenerate all `v4` figures and tables from `v4_revision/data_tables/`.
