# crop_double_water_stress

This repository contains the public code and final outputs for a crop double water stress analysis. The layout is intentionally compact: code is split by analysis stage, and final outputs are grouped under one results folder.

## Repository layout

- `code/` - readable analysis scripts split by workflow stage. Start with `code/README.md`.
- `results/final/` - final manuscript-facing tables, figures, validation outputs, and summary files. Start with `results/final/README.md`.

## Code structure

The original notebooks are not kept in this GitHub repository because they are large and hard to review. Their logic has been exported and split into staged scripts:

1. preprocessing and NetCDF preparation
2. data loading and yield-loss metrics
3. spatial loss maps
4. country/trade/risk analysis
5. validation workflows

## Reuse notes

The scripts preserve the analysis logic used for the manuscript, but several paths still point to local disks used during the original computation. Before rerunning the workflow elsewhere, update the path blocks and verify NetCDF dimensions, units, masks, and time coordinates.

One very large archival PDF is omitted from GitHub; compact PDF/PNG versions are included in the same result folder.
