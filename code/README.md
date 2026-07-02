# Code

The public code is split by analysis stage so the repository can be read without opening the original notebooks.

## Files

- `00_preprocess_model_inputs.py` - preprocessing, detrending, growing-season inputs, and NetCDF preparation.
- `01_data_loading_and_loss_metrics.py` - data readers, database assembly, event occurrence tables, and yield-loss metrics.
- `02_spatial_loss_maps.py` - historical and future spatial loss maps plus latitudinal summaries.
- `03_country_trade_and_risk_analysis.py` - country aggregation, FAOSTAT trade links, and risk plots/tables.
- `04_yield_and_stress_validation.py` - yield validation, GAEZ comparison, and stress-day time-series plotting.
- `05_admin_water_stress_validation.py` - administrative-unit validation of water-stress loss estimates.
- `06_historical_yield_validation.py` - historical yield validation for maize, rice, wheat, and soybean.

## Important note

These files are archival analysis scripts. They preserve the computational logic used for the manuscript, but several paths still point to local disks used during analysis. Before rerunning them elsewhere, update the path configuration blocks and verify the NetCDF dimensions, units, masks, and time coordinates.

The original notebooks are not included in the GitHub repository because they are large, harder to review, and mostly duplicate the split scripts.
