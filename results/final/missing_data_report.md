# Missing Data And Path Report

更新时间：2026-05-11

## 结论

多数关键外部目录仍存在，足以复查最终稿的大部分缓存结果和若干验证图。但代码中仍有一些旧路径、占位符路径或缺失数据路径。正式一键复现前需要你确认这些数据现在的新位置。

## 已确认存在的关键路径

- `D:\AAUDE\paper_v2\paper4`
- `D:\AAUDE\paper_v2\paper4\outputs_loss_maps`
- `D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss`
- `D:\AAUDE\paper_v2\paper4\validation_fao_yield`
- `F:\paper4`
- `F:\paper4\yie`
- `F:\paper4\data`
- `E:\paper4`
- `E:\paper4\cwatm_crop_yield`
- `E:\dsf`
- `F:\GAEZ_potential_yield`
- `F:\SPAM2020`
- `D:\Edge_download`
- `D:\Edge_download\cellarea.nc`
- `D:\Edge_download\FAOSTAT_data_en_3-8-2026.csv`
- `D:\Edge_download\FAOSTAT_data_en_3-8-2026 (1).csv`
- `D:\AAUDE\cil-ag-replication-package\Fig1\Crop_Coverage\data\shapes\all_countries.shp`
- `F:\crop_filed\impact_data\for_regressions_full_csv`

## 需要你确认位置的缺失路径

这些路径在 notebook 或脚本中被引用，但当前检查不存在：

- `F:\paper4\GWSP3_W5E5_v3\w_daily.nc`
- `F:\paper4\cwatm_crop_yield`
- `F:\paper4\cwatm_crop_yield_original`
- `F:\worldpop\mze_2000_prd.tif`
- `F:\worldpop\rcw_2000_prd.tif`
- `F:\worldpop\soy_2000_prd.tif`
- `F:\worldpop\whe_2000_prd.tif`
- `G:\CWatM\CWatM-Earth-30min\routing\cellarea.nc`
- `G:\global_cropmap\crop_calendar\crop_calendar_month\30-arcminute_fraction`
- `G:\global_cropmap\crop_calendar\crop_calendar_month\after_clip`
- `G:\global_cropmap\crop_calendar\ssp126`
- `G:\global_cropmap\crop_calendar\ssp585`
- `D:\AAUDE\paper_v2\paper4\outputs_loss_maps\Fig_loss_cropmean_2x3.pdf`
- `D:\AAUDE\paper_v2\paper4\outputs_loss_maps\crop_area_dist_0p5deg.nc`
- `D:\AAUDE\paper\paper4\figs\Fig1.pdf`
- `D:\AAUDE\paper\paper4\figs\Fig1.png`

## 可能只是旧路径或可替代路径

- `F:\paper4\cwatm_crop_yield` 当前不存在，但 `E:\paper4\cwatm_crop_yield` 存在，可能是迁移后的路径。
- `G:\CWatM\...\cellarea.nc` 当前不存在，但 `D:\Edge_download\cellarea.nc` 存在，可能可替代。
- `Fig_loss_cropmean_2x3.pdf` 不存在，但 `Fig_loss_cropmean_area_weighted_2x3.pdf` 存在，可能是更新版图件。
- `D:\AAUDE\paper\paper4\figs\Fig1.*` 不存在，但 `revision_v2\fig\manu_fig` 中有最终图相关导出，可能旧输出目录已经废弃。

## 占位符路径

这些不是缺失数据，而是 notebook 中的模板路径，需要运行时由变量替换：

- `F:\paper4\yie\{gcm}`
- `F:\paper4\yie\{test}`
- `F:\paper4\cwatm_crop_yield\{gcm}\{ssp}`
- `F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}`
- `E:\paper4\classification\{season}_{crop}_{arr}.nc`
- `E:\paper4\result2\{crop}_...`
- `E:\paper4\result3\{crop}_...`

## 对复现的影响

可以继续复查：

- 最终投稿 Word 文档。
- Supplementary Methods 和 response。
- 基于 `outputs_loss_maps` 的未来 loss map、variance decomposition、空间变化图。
- 基于 `validation_fao_yield` 的 FAO validation。
- 基于 `revision_v2/fig` 缓存的时间序列图和图件核查。

需要你提供或确认位置后才能完整重跑：

- 原始 CWatM daily water stress / daily soil moisture / SPEI 预处理。
- worldpop/GAEZ 2000 production comparison。
- 早期 G: 盘 crop calendar / CWatM-Earth 目录相关代码。

