# crop_double_water_stress

This repository contains the public code, modified CWatM modules, analysis workflows, and final manuscript-facing outputs for the crop double water stress analysis.

The study develops a process-based framework to quantify crop yield losses attributable to drought and waterlogging. The core model extension represents **two-sided crop water stress**, including moisture-deficit stress and excess-moisture aeration stress, and applies an event-based counterfactual framework to isolate abnormal water-stress impacts.

## Repository layout

- `Module/` - modified CWatM modules (`soil.py`, `evaporation.py`) implementing two-sided root-zone drought and waterlogging stress and counterfactual replacement.
- `settings/` - CWatM run configuration files (`settings_1_waterstress.ini`, `settings_2_cf_dry_only.ini`, `settings_3_cf_wet_only.ini`, `settings_4_cf_both.ini`), `generate_settings.py`, and crop hydraulic/stage parameters (`Table_S2_Soil_and_Crop_Hydraulic_Parameters.csv`).
- `code/` - analysis scripts split by workflow stage. Start with `code/README.md`.
- `results/final/` - final manuscript-facing tables, figures, validation outputs, and summary files. Start with `results/final/README.md`.
- `environment.yml`, `requirements.txt`, and `config_paths.yaml` - pinned Python environment definitions and centralized path configurations.

## Gridded Input Datasets

Study-prepared 0.5-degree global gridded input datasets (near-normal `k_s,norm` climatology `rws_climatology_time366_under_spei_-0.5_0.5_chunked.nc`, cropland land-cover fractions, FAO-56 growth-stage crop calendars, and aligned MIRCA-OS harvested-area rasters) are archived on Figshare at https://figshare.com/s/58602e5d64290c79976f.

## Model modifications

The modified CWatM implementation extends the original crop-water representation by including:

1. **Two-sided crop water stress representation**
   - dry-side water deficit stress under insufficient root-zone moisture;
   - wet-side aeration stress under excessive soil moisture and waterlogging conditions.

2. **Counterfactual yield attribution framework**
   - actual simulations retain fully simulated hydrological conditions;
   - counterfactual simulations replace abnormal drought/waterlogging stress factors with crop- and growth-stage-specific normal stress conditions while keeping other model settings unchanged.

3. **Crop- and stage-specific stress representation**
   - normal stress conditions are defined by crop, location, and growth stage rather than using a single global no-stress assumption.

4. **Improved crop-water and irrigation representation**
   - crop coefficients and yield-response parameters are explicitly handled;
   - rainfed and irrigated systems are separated;
   - paddy and upland irrigation processes are represented independently.

## Code structure

The analysis workflow is divided into stages:

1. preprocessing and NetCDF preparation
2. data loading and yield-loss metrics
3. spatial drought and waterlogging loss maps
4. country-level trade and risk analysis
5. validation workflows
6. statistical analysis and figure generation
