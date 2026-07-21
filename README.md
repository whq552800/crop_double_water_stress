# crop_double_water_stress

This repository contains the public code, modified CWatM modules, analysis workflows, and final manuscript-facing outputs for the crop double water stress analysis.

The study develops a process-based framework to quantify crop yield losses attributable to drought and waterlogging. The core model extension represents **two-sided crop water stress**, including moisture-deficit stress and excess-moisture aeration stress, and applies an event-based counterfactual framework to isolate abnormal water-stress impacts.

## Repository layout

- `Module/` - modified CWatM modules used for the crop-water stress extension. These modules implement the updated soil-water and evapotranspiration processes supporting drought and waterlogging attribution.
- `code/` - analysis scripts split by workflow stage. Start with `code/README.md`.
- `results/final/` - final manuscript-facing tables, figures, validation outputs, and summary files. Start with `results/final/README.md`.

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

## Reproducibility notes

The repository includes the analysis logic and manuscript-facing outputs. Large input datasets are not redistributed and should be obtained from their original public sources, including ISIMIP climate forcing, CMIP6 projections, SPAM crop distribution, GAEZ potential yield data, FAOSTAT production and trade data, and GGCMI crop calendar information.

Before rerunning the workflow:

- update local path configurations in analysis scripts;
- verify NetCDF dimensions, units, masks, and time coordinates;
- ensure that required input datasets are available.

## Software and code availability

The hydrological simulations are based on the open-source Community Water Model (CWatM). The modified CWatM modules used in this study are provided in `Module/`.

The analysis scripts used for preprocessing, yield-loss calculation, validation, uncertainty analysis, and figure generation are provided in `code/`.

The CWatM model is available at:

https://cwatm.iiasa.ac.at/

The complete repository provides the computational workflow supporting the manuscript results.