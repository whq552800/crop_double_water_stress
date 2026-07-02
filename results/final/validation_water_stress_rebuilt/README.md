# Rebuilt CIL-CWatM Water-Stress Validation

This folder contains a rebuilt validation workflow from raw inputs, not from the earlier precomputed validation CSVs.

## Scientific Framing

The comparison is:

- x-axis: CIL observed yield-loss anomaly, computed from administrative yield using a centered 5-year baseline.
- y-axis: CWatM simulated water-stress-attributable loss, rebuilt from waterstress vs nowaterstress model outputs.

Important interpretation boundary:

- CIL values are observed yield-loss anomalies during model-identified water-stress years.
- They are not direct observational attribution of drought or waterlogging causality.
- CWatM values are simulated losses attributable to water-stress treatment differences.

## Raw Inputs

- CIL yield panel CSVs: `F:\crop_filed\impact_data\for_regressions_full_csv`
- CWatM historical outputs: `F:\paper4\yie\gswp3_long`
- MIRCA crop area fractions: `F:\paper4\data\frac\30-arcminute_fraction`
- GAEZ historical potential yield: `F:\GAEZ_potential_yield\Hist_CRUTS32_nc`
- Grid cell area: `D:\Edge_download\cellarea.nc`
- CIL administrative polygons: `D:\AAUDE\cil-ag-replication-package\Fig1\Crop_Coverage\data\shapes\all_countries.shp`

## Rebuild Method

1. Read CIL crop yield at administrative-unit level.
2. Compute area-weighted ADM1 yield using harvested area, falling back to planted area.
3. Compute observed 5-year centered yield anomaly and convert negative anomalies to yield-loss percentage.
4. Read CWatM `gs_loss_*` files with `decode_times=False` because time units are `Years since 1901-01-01`.
5. Convert `gs_loss` to relative yield using the manuscript code convention: `relative_yield = 4 - gs_loss`.
6. Rebuild production as `relative_yield * GAEZ potential yield * MIRCA crop area`.
7. Compare `waterstress` and `nowaterstress` production to estimate simulated attributable losses.
8. Identify drought and waterlogging event years using stress-day files with a default threshold of 7 days.
9. Keep only events where the model-identified event area share is at least 0.20.
10. Merge CIL observed anomalies and CWatM simulated losses by crop, ISO, ADM1 and year.

Negative simulated differences are retained in raw columns, but the plotted loss variables use `max(loss, 0)` because the figure is specifically a yield-loss validation.

## Current Run

Script:

- `D:\Software\codex\drought&waterlogging\scripts\rebuild_cil_water_stress_validation.py`
- `D:\Software\codex\drought&waterlogging\scripts\plot_cached_cil_water_stress_validation.py`

Command:

```powershell
python .\scripts\rebuild_cil_water_stress_validation.py --crops maiz soyb rice whea --level adm1 --year-min 1980 --year-max 2019 --event-area-share-thr 0.2
```

Outputs:

- `adm1\1980_2019\model_cwatm_raw_rebuilt_adm1.csv`
- `adm1\1980_2019\obs_cil_5yr_anomaly_adm1.csv`
- `adm1\1980_2019\cil_cwatm_validation_points_adm1.csv`
- `adm1\1980_2019\cil_cwatm_validation_metrics_adm1.csv`
- `adm1\1980_2019\cil_cwatm_validation_hazard_loss_adm1.png`
- `adm1\1980_2019\cil_cwatm_validation_hazard_loss_adm1.pdf`
- `adm1\1980_2019\spatial_r_hazard_loss_adm1.png`
- `adm1\1980_2019\spatial_r_hazard_loss_adm1_source.csv`

The scatter panels report `n`, Pearson `r`, and `bias = simulated - observed` in percentage points. The spatial maps report temporal correlation `r` within each administrative unit.

The retained model-side validation variable is:

- `sim_hazard_loss_pct`: CWatM water-stress-attributable loss from the waterstress vs nowaterstress counterfactual.

Spatial map PDFs are not written by default because polygon-level vector PDFs can become hundreds of MB. The PNG plus source CSV are the reproducible outputs.

## Fast Replotting

The slow step is rebuilding `cil_cwatm_validation_points_adm1.csv` from raw CIL CSV, CWatM NetCDF, MIRCA crop area, GAEZ potential yield and administrative polygons. Once this cached table exists, use the plotting-only script:

```powershell
python .\scripts\plot_cached_cil_water_stress_validation.py --cache-dir "D:\Software\codex\drought&waterlogging\source_index\validation_water_stress_rebuilt\adm1\1980_2019" --level adm1
```

For fast scatter-only redraws:

```powershell
python .\scripts\plot_cached_cil_water_stress_validation.py --cache-dir "D:\Software\codex\drought&waterlogging\source_index\validation_water_stress_rebuilt\adm1\1980_2019" --level adm1 --skip-maps
```

The scatter-only redraw currently takes a few seconds because it only reads the cached validation points. The full cached plotting command also redraws spatial r maps and is slower because of shapefile rendering, but it does not rerun the raw data rebuild.

## Cached Yield Spatial Maps

For yield-level spatial validation, use:

```powershell
python .\scripts\plot_cached_yield_spatial_maps.py --cache-dir "D:\Software\codex\drought&waterlogging\source_index\validation_water_stress_rebuilt\adm1\1980_2019" --level adm1
```

This reads only cached `model_cwatm_raw_rebuilt_adm1.csv` and `obs_cil_5yr_anomaly_adm1.csv`. It does not rerun NetCDF or CIL raw CSV processing.

Outputs:

- `yield_spatial_level_stats_adm1.csv`: yield-level regional r and bias.
- `yield_spatial_r_adm1.png`: temporal r between observed and simulated yield.
- `yield_spatial_relative_bias_adm1.png`: mean relative yield bias.
- `yield_spatial_absolute_bias_adm1.png`: mean absolute yield bias.

Model yield is converted from t/ha to kg/ha before yield-level bias is computed.

## Yield-Loss Percentage-Point Bias Maps

For yield-loss validation, the requested bias is:

```text
bias = simulated loss (%) - observed loss (%)
```

This is a percentage-point difference, not a relative error divided by observed loss.

Relative bias maps of the form `(simulated - observed) / observed * 100` were removed because they become unstable and misleading when observed loss is near zero.

Current explicit outputs:

- `model_minus_observed_loss_percentage_point_bias_hazard_loss_adm1.png`
- `model_minus_observed_loss_percentage_point_bias_hazard_loss_adm1_source.csv`

The source CSV field is `mean_loss_bias_pctpt`. Negative values mean the model-estimated loss is smaller than the CIL observed yield-loss anomaly in matched event years.

All current spatial validation maps follow the manuscript-style land display rule: administrative/land polygons without valid values are shown in light grey (`#eeeeee`), while polygons with validation values are colored by the target variable.

## Observed and Simulated Loss Maps

For spatial maps of the loss percentages themselves, use:

```powershell
python .\scripts\plot_cached_observed_simulated_loss_maps.py --cache-dir "D:\Software\codex\drought&waterlogging\source_index\validation_water_stress_rebuilt\adm1\1980_2019" --level adm1 --sim-col sim_hazard_loss_pct --stem hazard_loss
```

This keeps all valid matched event years by default, including years with zero observed loss. Do not pass `--min-obs-loss` unless a sensitivity test intentionally excludes small observed losses.

Outputs for the model water-stress-attributable loss comparison:

- `observed_loss_percent_hazard_loss_adm1.png`
- `simulated_loss_percent_hazard_loss_adm1.png`
- `observed_simulated_loss_percent_hazard_loss_adm1_source.csv`

## Current Caveats

- The current run is ADM1 only. ADM2 is technically supported by the script, but boundary and ID matching should be audited before using it in the manuscript.
- The CIL observed data available in the loaded CSVs reach 2015 for the current matched crop panel, while CWatM outputs reach 2019.
- The figure is a validation diagnostic, not yet a final journal-layout figure.
- The modeled losses are generally smaller than observed yield-loss anomalies, which is expected because observed yield anomalies include non-water-stress influences and management/noise.
