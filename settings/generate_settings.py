"""Generate the 4 scenario settings files from the user's base settings template."""
from pathlib import Path

OUT = Path(r"d:\Software\Gemini\drought_waterlogging\cwatm_counterfactual_patch")

BASE = r"""# ------------------------------------------------
# Community Water Model Version 0.99 - Counterfactual Attribution Settings
# Scenario: {scen_desc}
# ------------------------------------------------

[OPTIONS]
modflow_coupling = False
TemperatureInKelvin = True
gridSizeUserDefined = True
calc_evaporation = True

includeIrrigation = True
includeWaterDemand = True
usingAllocSegments = False
limitAbstraction = True
fossil_water_treated_normally = True
sectorSourceAbstractionFractions = True
calc_environflow = False

preferentialFlow = False
CapillarRise = True

includeRunoffConcentration = True
includeWaterBodies = True
includeRouting = True
inflow = False

writeNetcdfStack = True
reportMap = True
reportTss = True
calcWaterBalance = False
sumWaterBalance = False
PCRaster = False

# --- Water-stress counterfactual switches ---
nowaterstress = {nowaterstress}
cf_remove_dry = {cf_remove_dry}
cf_remove_wet = {cf_remove_wet}

includeCrops = True
use_GeneralCropnonIrr = 0
use_GeneralCropIrr = 0
activate_fallow = 0
automaticFallowingIrr = True
leftoverIrrigatedCropIsRainfed = 0
moveIrrFallowToNonIrr = False

[FILE_PATHS]
PathRoot = G:\CWatM\CWatM-Earth-30min
GCMM = CNRM-ESM2-1
prefilename = obsclim
# 请在运行电脑上确认 PathOut 与 PathMeteo 指向 paper4 的历史 GSWP3-W5E5 数据：
PathOut = L:\paper4\gswp3_attribution\{scen_dir}
PathMaps = $(PathRoot)
PathMeteo = L:\ISIMIP4_cimp7\ISIMIP4a\Input Data\Cliamte related forcing\Atmospheric forcings\GSWP3-ERA5\obsclim
PathWaterdemand = $(PathMaps)/landsurface/waterDemand
Excel_settings_file = $(PathMaps)\settings\20251205.xlsx
PathGroundwaterModflow = $(PathMaps)/ModFlow_inputs
PathGroundwaterModflowOutput = G:\CWatM\northeast_modflow_cwatm\modflow\111
co2conc = $(PathRoot)\RCP_CO2\rcp85.nc4

[NETCDF_ATTRIBUTES]
institution = IIASA
title = CWatM Crop Water Stress Attribution - {scen_dir}
metaNetcdfFile = $(FILE_PATHS:PathMaps)/metaNetcdf.xml

[MASK_OUTLET]
MaskMap = $(FILE_PATHS:PathRoot)/areamaps/global.map
Gauges = 111.25 27.25
GaugesLocal = True

[TIME-RELATED_CONSTANTS]
StepStart = 1/1/1975
SpinUp = 1/1/1980
StepEnd = 31/12/2019

[INITITIAL CONDITIONS]
load_initial = True
initLoad = $(FILE_PATHS:PathRoot)/init/global_yield_v2_19800101.nc
save_initial = False
initSave = $(FILE_PATHS:PathRoot)/init/global_yield_v2
StepInit = 31/12/2019

[CALIBRATION]
SnowMeltCoef = 0.004
crop_correct = 1.7387144425028453
soildepth_factor = 1.2465457814266776
preferentialFlowConstant = 3.16349443236206
arnoBeta_add = 0.03885276991497097
factor_interflow = 1.3119554375939964
recessionCoeff_factor = 5.378109968185251
runoffConc_factor = 1
manningsN = 9.232881632147917
normalStorageLimit = 0.5457611123744854
lakeAFactor = 1
lakeEvaFactor = 1

[TOPOP]
Ldd = $(FILE_PATHS:PathMaps)/routing/ldd.map
ElevationStD = $(FILE_PATHS:PathMaps)/landsurface/topo/elvstd.map
CellArea = $(FILE_PATHS:PathMaps)/routing/cellarea.map

[METEO]
PrecipitationMaps = $(FILE_PATHS:PathMeteo)/*$(FILE_PATHS:prefilename)*pr*
TavgMaps = $(FILE_PATHS:PathMeteo)/*$(FILE_PATHS:prefilename)*tas*
SPEIfile = F:\spei\gswp3\SPEI_daily_*.nc
rwsfile = L:\cwatm_crop_yield\rws_climatology_time366_under_spei_-0.5_0.5_chunked.nc
E0Maps = F:\CWatM\northeast_modflow_cwatm\222\EWRef_daily.nc
ETMaps = F:\CWatM\northeast_modflow_cwatm\222\ETRef_daily.nc
precipitation_coversion = 86.4
evaporation_coversion = 1.00

[EVAPORATION]
albedo = True
albedoMaps = $(FILE_PATHS:PathMaps)/landsurface/albedo/albedo.nc
AlbedoSoil = 0.15
AlbedoWater = 0.05
AlbedoCanopy = 0.23
useHuss = True
TminMaps = $(FILE_PATHS:PathMeteo)/*$(FILE_PATHS:prefilename)*tasmin*
TmaxMaps = $(FILE_PATHS:PathMeteo)/*$(FILE_PATHS:prefilename)*tasmax*
PSurfMaps = $(FILE_PATHS:PathMeteo)/*$(FILE_PATHS:prefilename)*ps*
RhsMaps = $(FILE_PATHS:PathMeteo)/*$(FILE_PATHS:prefilename)*hurs*
QAirMaps = $(FILE_PATHS:PathMeteo)/*$(FILE_PATHS:prefilename)*huss*
WindMaps = $(FILE_PATHS:PathMeteo)/*$(FILE_PATHS:prefilename)*sfcwind*
RSDSMaps = $(FILE_PATHS:PathMeteo)/*$(FILE_PATHS:prefilename)*rsds*
RSDLMaps = $(FILE_PATHS:PathMeteo)/*$(FILE_PATHS:prefilename)*rlds*

[SNOW]
NumberSnowLayers = 7
GlacierTransportZone = 3
TemperatureLapseRate = 0.0065
SnowFactor = 1.0
SnowSeasonAdj = 0.001
TempMelt = 1.0
TempSnow = 1.0
IceMeltCoef = 0.007

[FROST]
SnowWaterEquivalent = 0.45
Afrost = 0.97
Kfrost = 0.57
FrostIndexThreshold = 56

[VEGETATION]
cropgroupnumber = $(FILE_PATHS:PathMaps)/soil/cropgrp.nc

[SOIL]
PathTopo = $(FILE_PATHS:PathMaps)/landsurface/topo
PathSoil = $(FILE_PATHS:PathMaps)/soil
tanslope = $(PathTopo)/tanslope.map
slopeLength = $(PathTopo)/slopeLength.map
relativeElevation = $(PathTopo)/dzRel_hydro1k.nc
KSat1 = $(PathSoil)/ksat1.map
KSat2 = $(PathSoil)/ksat2.map
KSat3 = $(PathSoil)/ksat3.map
alpha1 = $(PathSoil)/alpha1.map
alpha2 = $(PathSoil)/alpha2.map
alpha3 = $(PathSoil)/alpha3.map
lambda1 = $(PathSoil)/lambda1.map
lambda2 = $(PathSoil)/lambda2.map
lambda3 = $(PathSoil)/lambda3.map
thetas1 = $(PathSoil)/thetas1.map
thetas2 = $(PathSoil)/thetas2.map
thetas3 = $(PathSoil)/thetas3.map
thetar1 = $(PathSoil)/thetar1.map
thetar2 = $(PathSoil)/thetar2.map
thetar3 = $(PathSoil)/thetar3.map
percolationImp = $(PathSoil)/percolationImp.map
maxGWCapRise = 5.0
minCropKC = 0.2
minTopWaterLayer = 0.0
StorDepth1 = $(PathSoil)/storageDepth1.map
StorDepth2 = $(PathSoil)/storageDepth2.map

[LANDCOVER]
PathLandcover = $(FILE_PATHS:PathMaps)/landsurface
coverTypes = forest, grassland, irrPaddy, irrNonPaddy, sealed, water
coverTypesShort = f, g, i, n, s, w
fractionLandcover = G:\global_cropmap\crop_calendar\crop_calendar_month\30-arcminute_fraction\fractionLandcover_crop_adjusted_with_rainfed.nc
dynamicLandcover = True
fixLandcoverYear = 1961
arrrr = G:\global_cropmap\crop_calendar\crop_calendar_month\30-arcminute_fraction
wwh_Irr = $(arrrr)\MIRCA-OS_Wheat_2000_irr_30arcmin_winter_aligned.tif
swh_Irr = $(arrrr)\MIRCA-OS_Wheat_2000_irr_30arcmin_spring_aligned.tif
rice_Irr = $(arrrr)\MIRCA-OS_Rice_2000_ir_30arcmin_aligned.tif
maiz_Irr = $(arrrr)\MIRCA-OS_Maize_2000_ir_30arcmin_aligned.tif
soyb_Irr = $(arrrr)\MIRCA-OS_Soybeans_2000_ir_30arcmin_aligned.tif
wwh_nonIrr = $(arrrr)\MIRCA-OS_Wheat_2000_rf_30arcmin_winter_aligned.tif
swh_nonIrr = $(arrrr)\MIRCA-OS_Wheat_2000_rf_30arcmin_spring_aligned.tif
rice_nonIrr = $(arrrr)\MIRCA-OS_Rice_2000_rf_30arcmin_aligned.tif
maiz_nonIrr = $(arrrr)\MIRCA-OS_Maize_2000_rf_30arcmin_aligned.tif
soyb_nonIrr = $(arrrr)\MIRCA-OS_Soybeans_2000_rf_30arcmin_aligned.tif

arr = G:\global_cropmap\crop_calendar\crop_calendar_month\after_clip
wwh_Irr_calendar = $(arr)\wwh_irr_gs_nodes.nc
swh_Irr_calendar = $(arr)\swh_irr_gs_nodes.nc
rice_Irr_calendar = $(arr)\rice_irr_gs_nodes.nc
maiz_Irr_calendar = $(arr)\maiz_irr_gs_nodes.nc
soyb_Irr_calendar = $(arr)\soyb_irr_gs_nodes.nc
wwh_nonIrr_calendar = $(arr)\wwh_rf_gs_nodes.nc
swh_nonIrr_calendar = $(arr)\swh_rf_gs_nodes.nc
rice_nonIrr_calendar = $(arr)\rice_rf_gs_nodes.nc
maiz_nonIrr_calendar = $(arr)\maiz_rf_gs_nodes.nc
soyb_nonIrr_calendar = $(arr)\soyb_rf_gs_nodes.nc

[__forest]
PathForest = $(FILE_PATHS:PathMaps)/landcover/forest
PathSoil = $(FILE_PATHS:PathMaps)/soil
forest_arnoBeta = 0.2
forest_KSat1 = $(PathSoil)/forest_ksat1.map
forest_KSat2 = $(PathSoil)/forest_ksat2.map
forest_KSat3 = $(PathSoil)/ksat3.map
forest_alpha1 = $(PathSoil)/forest_alpha1.map
forest_alpha2 = $(PathSoil)/forest_alpha2.map
forest_alpha3 = $(PathSoil)/alpha3.map
forest_lambda1 = $(PathSoil)/forest_lambda1.map
forest_lambda2 = $(PathSoil)/forest_lambda2.map
forest_lambda3 = $(PathSoil)/lambda3.map
forest_thetas1 = $(PathSoil)/forest_thetas1.map
forest_thetas2 = $(PathSoil)/forest_thetas2.map
forest_thetas3 = $(PathSoil)/thetas3.map
forest_thetar1 = $(PathSoil)/forest_thetar1.map
forest_thetar2 = $(PathSoil)/forest_thetar2.map
forest_thetar3 = $(PathSoil)/thetar3.map
forest_minInterceptCap = 0.001
forest_cropDeplFactor = 0.0
forest_fracVegCover = $(PathForest)/fracVegCover.map
forest_rootFraction1 = $(PathForest)/rootFraction1.map
forest_rootFraction2 = $(PathForest)/rootFraction2.map
forest_maxRootDepth = $(PathForest)/maxRootDepth.map
forest_minSoilDepthFrac = $(PathForest)/minSoilDepthFrac.map
forest_cropCoefficientNC = $(PathForest)/cropCoefficientForest_10days.nc
forest_interceptCapNC = $(PathForest)/interceptCapForest_10days.nc

[__grassland]
PathGrassland = $(FILE_PATHS:PathMaps)/landcover/grassland
grassland_arnoBeta = 0.01
grassland_minInterceptCap = 0.001
grassland_cropDeplFactor = 0.3
grassland_fracVegCover = $(PathGrassland)/fracVegCover.map
grassland_rootFraction1 = $(PathGrassland)/rootFraction1.map
grassland_rootFraction2 = $(PathGrassland)/rootFraction2.map
grassland_maxRootDepth = $(PathGrassland)/maxRootDepth.map
grassland_minSoilDepthFrac = $(PathGrassland)/minSoilDepthFrac.map
grassland_cropCoefficientNC = $(PathGrassland)/cropCoefficientGrassland_10days.nc
grassland_interceptCapNC = $(PathGrassland)/interceptCapGrassland_10days.nc

[__irrPaddy]
PathIrrPaddy = $(FILE_PATHS:PathMaps)/landcover/irrPaddy
irrPaddy_arnoBeta = 0.4
irrPaddy_minInterceptCap = 0.002
irrPaddy_cropDeplFactor = 0.0
irrPaddy_fracVegCover = $(PathIrrPaddy)/fracVegCover.map
irrPaddy_rootFraction1 = $(PathIrrPaddy)/rootFraction1.map
irrPaddy_rootFraction2 = $(PathIrrPaddy)/rootFraction2.map
irrPaddy_maxRootDepth = $(PathIrrPaddy)/maxRootDepth.map
irrPaddy_minSoilDepthFrac = $(PathIrrPaddy)/minSoilDepthFrac.map
irrPaddy_cropCoefficientNC = $(PathIrrPaddy)/cropCoefficientirrPaddy_10days.nc
irrPaddy_maxtopwater = 0.05

[__irrNonPaddy]
PathIrrNonPaddy = $(FILE_PATHS:PathMaps)/landcover/irrNonPaddy
irrNonPaddy_arnoBeta = 0.2
irrNonPaddy_minInterceptCap = 0.001
irrNonPaddy_cropDeplFactor = 0.0
irrNonPaddy_fracVegCover = $(PathIrrNonPaddy)/fracVegCover.map
irrNonPaddy_rootFraction1 = $(PathIrrNonPaddy)/rootFraction1.map
irrNonPaddy_rootFraction2 = $(PathIrrNonPaddy)/rootFraction2.map
irrNonPaddy_maxRootDepth = $(PathIrrNonPaddy)/maxRootDepth.map
irrNonPaddy_minSoilDepthFrac = $(PathIrrNonPaddy)/minSoilDepthFrac.map
irrNonPaddy_cropCoefficientNC = $(PathIrrNonPaddy)/cropCoefficientirrNonPaddy_10days.nc

[__sealed]
PathSealed = $(FILE_PATHS:PathMaps)/landcover/sealed
sealed_minInterceptCap = 0.001

[__open_water]
PathWater = $(FILE_PATHS:PathMaps)/landcover/water
water_minInterceptCap = 0.0

[GROUNDWATER]
PathGroundwater = $(FILE_PATHS:PathMaps)/groundwater
recessionCoeff = $(PathGroundwater)/recessionCoeff.map
specificYield = $(PathGroundwater)/specificYield.map
kSatAquifer = $(PathGroundwater)/kSatAquifer.map

[GROUNDWATER_MODFLOW]
path_mf6dll = G:\CWatM\CWatM\cwatm
PathGroundwaterModflowOutput = F:\CWatM\northeast_modflow_cwatm\modflow\northeast_1960
PathGroundwaterModflowClimate = F:\CWatM\northeast_modflow_cwatm\climate
PathGroundwaterModflowInput = $(FILE_PATHS:PathGroundwaterModflow)/$(Modflow_resolution)m
writeModflowError = True
start_soil_humid = 0.4
load_init_water_table = False
initial_water_table_depth = 2
init_water_table = C:/Users/guillaumot/Documents/CWatM_Burgenland_1km/CWatM-ModFlow6/output_demandV4/water_table_map_previous.npy
modflow_steadystate = False
load_steady_previous = True
Ndays_steady = 500
ModflowPrecipitationMaps = $(PathGroundwaterModflowClimate)/pr_avg_reprojected.nc
ModflowTavgMaps = $(PathGroundwaterModflowClimate)/tas_avg_reprojected.nc
ModflowE0Maps = $(PathGroundwaterModflowClimate)/EWRefavg_reprojected.nc
ModflowETMaps = $(PathGroundwaterModflowClimate)/ETRefavg_reprojected.nc
use_soildepth_as_GWtop = False
correct_soildepth_underlakes = True
depth_underlakes = 1.5
Groundwater_pumping = True
water_table_limit_for_pumping = 0.45
Pumping_input_file = E:\\jupyter\\modflow_GD_pumping.npy
modflow_steadyInit = $(PathGroundwaterModflowOutput)/Nira_steady_April9_3.hds
nameModflowModel = northeast
Modflow_resolution = 1000
modflow_timestep = 1
load_modflow_from_disk = False
nlay = 1
cwatm_modflow_indices = $(PathGroundwaterModflowInput)/indices
modflow_basin = $(PathGroundwaterModflowInput)/modflow_basin.tif
topo_modflow = $(PathGroundwaterModflowInput)/elevation_modflow.tif
chanRatio = $(PathGroundwaterModflowInput)/modlfow_river_percentage.tif
thickness = 100.0
poro = 0.05
permeability = 2e-05
leakageriver_permea = 0.001
leakagelake_permea = 0.001

[WATERDEMAND]
allocSegments = $(FILE_PATHS:PathWaterdemand)/catchx.nc
demand_unit = False
historicalIrrigationArea = $(FILE_PATHS:PathWaterdemand)/irrigationArea.nc
domesticWaterDemandFile = $(FILE_PATHS:PathWaterdemand)/dom_isimip_year_millionm3_1850_2018.nc
domesticTimeMonthly = False
domesticWithdrawalvarname = domww
domesticConsuptionvarname = domwc
industryWaterDemandFile = $(FILE_PATHS:PathWaterdemand)/ind_isimip_year_millionm3_1850_2018.nc
industryTimeMonthly = False
industryWithdrawalvarname = indww
industryConsuptionvarname = indwc
uselivestock = True
livestockWaterDemandFile = $(FILE_PATHS:PathWaterdemand)/liv_isimip_year_millionm3_1850_2018.nc
livestockTimeMonthly = False
livestockvarname = livestockConsumption
irrNonPaddy_efficiency = $(FILE_PATHS:PathMaps)/landsurface/waterDemand/efficiency.nc
irrPaddy_efficiency = $(FILE_PATHS:PathMaps)/landsurface/waterDemand/efficiency.nc
use_environflow = False
EnvironmentalFlowFile = $(FILE_PATHS:PathOut)/MQ90_12month.nc
irrigation_returnfraction = 0.5
swAbstractionFraction_Channel_Domestic = 1
swAbstractionFraction_Channel_Livestock = 1
swAbstractionFraction_Channel_Industry = 1
swAbstractionFraction_Channel_Irrigation = 1
swAbstractionFraction_Lift_Domestic = 0
swAbstractionFraction_Lift_Livestock = 0
swAbstractionFraction_Lift_Industry = 0
swAbstractionFraction_Lift_Irrigation = 0
swAbstractionFraction_Lake_Domestic = 1
swAbstractionFraction_Lake_Livestock = 1
swAbstractionFraction_Lake_Industry = 1
swAbstractionFraction_Lake_Irrigation = 1
swAbstractionFraction_Res_Domestic = 1
swAbstractionFraction_Res_Livestock = 1
swAbstractionFraction_Res_Industry = 1
swAbstractionFraction_Res_Irrigation = 1
gwAbstractionFraction_Domestic = 1
gwAbstractionFraction_Livestock = 1
gwAbstractionFraction_Industry = 1
gwAbstractionFraction_Irrigation = 1
swAbstractionFrac = 0.9
averageDischarge = $(FILE_PATHS:PathOut)/discharge_totalavg_rhine30min.nc
averageBaseflow = $(FILE_PATHS:PathOut)/baseflow_totalavg_rhine30min.nc
baseflowInM = True

[RUNOFF_CONCENTRATION]
forest_runoff_peaktime = 1.0
grassland_runoff_peaktime = 0.5
irrPaddy_runoff_peaktime = 0.5
irrNonPaddy_runoff_peaktime = 0.5
sealed_runoff_peaktime = 0.15
water_runoff_peaktime = 0.01
interflow_runoff_peaktime = 1.0
baseflow_runoff_peaktime = 2.0

[ROUTING]
PathRouting = $(FILE_PATHS:PathMaps)/routing
NoRoutingSteps = 10
chanBeta = 0.6
chanGrad = $(PathRouting)/kinematic/changrad.nc
chanGradMin = 0.0001
chanMan = $(PathRouting)/kinematic/chanman.nc
chanLength = $(PathRouting)/kinematic/chanleng.nc
chanWidth = $(PathRouting)/kinematic/chanbw.nc
chanDepth = $(PathRouting)/kinematic/chanbnkf.nc

[LAKES_RESERVOIRS]
PathLakesRes = $(FILE_PATHS:PathMaps)/routing/lakesreservoirs
useResAndLakes = True
dynamicLakesRes = True
fixLakesResYear = 1950
waterBodyID = $(PathLakesRes)/lakesResID.nc
waterBodyTyp = $(PathLakesRes)/lakesResType.nc
waterBodyDis = $(PathLakesRes)/lakesResDis.nc
waterBodyArea = $(PathLakesRes)/lakesResArea.nc
buffer_waterbodies = 2
useSmallLakes = False
smallLakesRes = $(PathLakesRes)/smallLakesRes.nc
smallwaterBodyDis = $(PathLakesRes)/smalllakesresDis.nc
waterBodyVolRes = $(PathLakesRes)/lakesResVolRes.nc
waterBodyYear = $(PathLakesRes)/lakesResYear.nc
conservativeStorageLimit = 0.1
floodStorageLimit = 0.9
adjust_Normal_Flood = 0.5
alphaDepletion = 0.9
MinOutflowQ = 0.2
NormalOutflowQ = 1.0
NonDamagingOutflowQ = 4.0

[INFLOW]
In_Dir = $(FILE_PATHS:PathRoot)/in
InflowPoints = $(In_Dir)/in.map
QInTS = mm.tss

[ENVIRONMENTALFLOW]
calc_ef_after = True
calc_ef_afterRun = True
EFDis = $(FILE_PATHS:PathOut)/discharge_rhine.nc
cut_ef_map = False
OUT_Dir = $(FILE_PATHS:PathOut)

[OUTPUT]
OUT_Dir = $(FILE_PATHS:PathOut)
OUT_MAP_AnnualTot = act_irrPaddyWithdrawal,irr_crop,gs_loss_nonIrr,gs_loss_Irr,drystress_Irr,wetstress_Irr,drystress_nonIrr,wetstress_nonIrr
"""

scenarios = [
    ("settings_1_waterstress.ini", "waterstress", "False", "False", "False", "Actual stressed baseline (nowaterstress = False)"),
    ("settings_2_cf_dry_only.ini", "cf_dry_only", "True",  "True",  "False", "Counterfactual removing DROUGHT stress only"),
    ("settings_3_cf_wet_only.ini", "cf_wet_only", "True",  "False", "True",  "Counterfactual removing WATERLOGGING stress only"),
    ("settings_4_cf_both.ini",     "cf_both",     "True",  "True",  "True",  "Counterfactual removing BOTH drought and waterlogging stress"),
]

for fname, sdir, nws, rdry, rwet, desc in scenarios:
    (OUT / fname).write_text(
        BASE.format(scen_dir=sdir, nowaterstress=nws, cf_remove_dry=rdry, cf_remove_wet=rwet, scen_desc=desc),
        encoding="utf-8",
    )
print("Wrote 4 settings files to", OUT)
