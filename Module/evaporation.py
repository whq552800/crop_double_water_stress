# -------------------------------------------------------------------------
# Name:        Evaporation module
# Purpose: Actual evapotranspiration calculation module for different land cover types.
# Processes crop coefficients and calculates land cover specific evapotranspiration rates.
# Handles bare soil evaporation and vegetation-specific water consumption.
#
# Author:      PB, MS, DF, JdB
# Created:     01/08/2016
# CWatM is licensed under GNU GENERAL PUBLIC LICENSE Version 3.
# -------------------------------------------------------------------------

from cwatm.management_modules.data_handling import *
import re

class evaporation(object):
    """
    Evaporation module for hydrological modeling.

    This class handles the calculation of potential evaporation and potential transpiration
    for different land cover types. It processes crop coefficients, calculates bare soil
    evaporation, and manages crop-specific evapotranspiration calculations.

    Attributes
    ----------
    var : object
        Model variables container
    model : object
        CWatM model instance

    **Global variables**
    ===================================  ==========    ======================================================================  =====
    Variable [self.var]                  Type          Description                                                             Unit 
    ===================================  ==========    ======================================================================  =====
    cropKCmonth                          Array         Crop KC factor for different crops and different seasons                --   
    snowEvap                             Array         total evaporation from snow for a snow layers                           m    
    iceEvap                              Array         Evaporation from ice (sublimation)                                      m    
    Crops_names                          Array         Internal: List of specific crops                                        --   
    activatedCrops                       Array         Fraction of area a specific crop is planted                             --   
    load_initial                         Flag          Settings initLoad holds initial conditions for variables                bool 
    monthCounter                         Array         Month counter for each crop after crop has planted                      --   
    fracCrops_IrrLandDemand              Array         Month counter for each crop after crop has planted                      --   
    fracCrops_nonIrrLandDemand                                                                                                 --   
    ratio_a_p_nonIrr                     Array         Ratio actual to potential evapotranspiration, monthly, non-irrigated [  %    
    totalPotET_month                     Array         Total potential evapotranspiration in a month                           m    
    ratio_a_p_Irr                        Array         Ratio actual to potential evapotranspiration, monthly [crop specific]   %    
    Yield_nonIrr                         Array         Relative monthly non-irrigated yield [crop specific]                    %    
    currentKY                            Array         Yield sensitivity coefficient [crop specific]                           --   
    Yield_Irr                            Array         Relative monthly irrigated yield [crop specific]                        %    
    currentKC                            Array         Current crop coefficient for specific crops                             --   
    generalIrrCrop_max                   Array                                                                                 --   
    generalnonIrrCrop_max                Array                                                                                 --   
    weighted_KC_nonIrr                   Array                                                                                 --   
    weighted_KC_nonIrr_woFallow          Array                                                                                 --   
    weighted_KC_Irr                      Array                                                                                 --   
    _weighted_KC_Irr                     Array                                                                                 --   
    weighted_KC_Irr_woFallow             Array                                                                                 --   
    totalPotET_month_segment             Array                                                                                 --   
    PotETaverage_crop_segments           Array                                                                                 --   
    areaCrops_Irr_segment                Array                                                                                 --   
    areaCrops_nonIrr_segment             Array                                                                                 --   
    areaPaddy_Irr_segment                Array                                                                                 --   
    Precipitation_segment                Array                                                                                 --   
    availableArableLand_segment          Array                                                                                 --   
    cropCorrect                          Array         calibration factor of crop KC factor                                    --   
    crop_correct_landCover               Array                                                                                 --   
    includeCrops                         Flag          1 when includeCrops=True in Settings, 0 otherwise                       bool 
    Crops                                Array         Internal: List of specific crops and Kc/Ky parameters                   --   
    daily_crop_KC                        Array                                                                                 --   
    interceptCap                         Array         interception capacity of vegetation                                     m    
    potTranspiration                     Array         Potential transpiration (after removing of evaporation)                 m    
    cropKC                               Array         crop coefficient for each of the 4 different land cover types (forest,  --   
    minCropKC                            Array         minimum crop factor (default 0.2)                                       --   
    minInterceptCap                      Array         Maximum interception read from file for forest and grassland land cove  m    
    irrigatedArea_original               Array                                                                                 --   
    frac_totalnonIrr                     Array         Fraction sown with specific non-irrigated crops                         %    
    frac_totalIrr_max                    Array         Fraction sown with specific irrigated crops, maximum throughout simula  %    
    frac_totalnonIrr_max                 Array         Fraction sown with specific non-irrigated crops, maximum throughout si  %    
    GeneralCrop_Irr                      Array         Fraction of irrigated land class sown with generally representative cr  %    
    fallowIrr                            Array         Fraction of fallowed irrigated land                                     %    
    fallowIrr_max                        Array         Fraction of fallowed irrigated land, maximum throughout simulation      %    
    GeneralCrop_nonIrr                   Array         Fraction of grasslands sown with generally representative crop          %    
    fallownonIrr                         Array         Fraction of fallowed non-irrigated land                                 %    
    fallownonIrr_max                     Array         Fraction of fallowed non-irrigated land, maximum throughout simulation  %    
    availableArableLand                  Array         Fraction of land not currently planted with specific crops              %    
    ETRef                                Array         potential evapotranspiration rate from reference crop                   m    
    Precipitation                        Array         Precipitation (input for the model)                                     m    
    coverTypes                           Array         land cover types - forest - grassland - irrPaddy - irrNonPaddy - water  --   
    SnowMelt                             Array         total snow melt from all layers                                         m    
    IceMelt                              Array         Ice melt (not really ice but an additional snow melt in summer)         m    
    potBareSoilEvap                      Array         potential bare soil evaporation (calculated with minus snow evaporatio  m    
    irr_Paddy_month                      Array                                                                                 --   
    ET_crop_Irr_paddy                    Array                                                                                 --   
    ET_crop_Irr_paddy_fraccrop           Array                                                                                 --   
    fracCrops_Irr                        Array         Fraction of cell currently planted with specific irrigated crops        %    
    fracCrops_nonIrr                     Array         Fraction of cell currently planted with specific non-irr crops          %    
    actTransTotal_month_nonIrr           Array         Internal variable: Running total of  transpiration for specific non-ir  m    
    actTransTotal_month_Irr              Array         Internal variable: Running total of  transpiration for specific irriga  m    
    irr_crop_month                       Number                                                                                --   
    frac_totalIrr                        Array         Fraction sown with specific irrigated crops                             %    
    weighted_KC_Irr_woFallow_fullKc      Array                                                                                 --   
    totalPotET                           Array         Potential evaporation per land use class                                m    
    PotET_crop                           Array                                                                                 --   
    fracVegCover                         Array         Fraction of specific land covers (0=forest, 1=grasslands, etc.)         %    
    adminSegments                        Array         Domestic agents                                                         Int  
    cellArea                             Array         Area of cell                                                            m2   
    ===================================  ==========    ======================================================================  =====

    """

    def __init__(self, model):
        """
        Initialize the evaporation module.

        Parameters
        ----------
        model : object
            CWatM model instance containing variables and configuration
        """
        self.var = model.var
        self.model = model

    def initial(self):
        """
        Initialize evaporation module arrays and parameters.

        Sets up crop coefficient arrays, interception capacity arrays, and reads
        initial data for different cover types including forest, grassland, and
        irrigated crops. Initializes monthly crop coefficient data from NetCDF files.
        """
        # no_types = len (self.var.coverTypes)
        self.var.cropKCmonth = np.zeros((4, 13, len(globals.inZero)))
        self.var.cropKC = np.zeros((4, len(globals.inZero)))
        self.var.interceptCap = np.zeros((2, 13, len(globals.inZero)))
        j = 0
        for coverType in self.var.coverTypes:

            if coverType in ['forest', 'grassland', 'irrPaddy', 'irrNonPaddy']:
                for i in range(13):
                    self.var.cropKCmonth[j, i, :] = readnetcdf2(coverType + '_cropCoefficientNC', i * 3, "10day")
                    self.var.cropKCmonth[j, i, :] = np.maximum(self.var.cropKCmonth[j, i, :], self.var.minCropKC)
                iii = 1

            if coverType in ['forest', 'grassland']:
                for i in range(13):
                    self.var.interceptCap[j, i, :] = readnetcdf2(coverType + '_interceptCapNC', i * 3, "10day")
                    self.var.interceptCap[j, i, :] = np.maximum(self.var.interceptCap[j, i, :], self.var.minInterceptCap[j])
            j = j + 1
        ii = 1

    def dynamic(self, coverType, No):
        """
        Calculate potential evapotranspiration for a specific land cover type.

        This method computes potential evaporation and transpiration using crop coefficients,
        handles crop dynamics when crops are enabled, and calculates bare soil evaporation.
        It processes monthly crop coefficient data with daily interpolation.

        Parameters
        ----------
        coverType : str
            Land cover type identifier (e.g., 'forest', 'grassland', 'irrPaddy')
        No : int
            Numerical identifier for land cover type (forest=0, grassland=1, etc.)

        Returns
        -------
        tuple
            Potential evaporation from bare soil and potential transpiration values
        """

        # get crop coefficient
        # to get ETc from ET0 x kc factor  ((see http://www.fao.org/docrep/X0490E/x0490e04.htm#TopOfPage figure 4:)
        # crop coefficient read for forest and grassland from file



        # calculate potential bare soil evaporation - only once
        if No == 0:
            self.var.potBareSoilEvap = self.var.cropCorrect * self.var.minCropKC * self.var.ETRef
            # calculate snow and ice evaporation
            self.var.snowEvap = np.minimum(self.var.SnowMelt, self.var.potBareSoilEvap)
            self.var.potBareSoilEvap -= self.var.snowEvap

            self.var.iceEvap = np.minimum(self.var.IceMelt, self.var.potBareSoilEvap)
            self.var.potBareSoilEvap -= self.var.iceEvap

            self.var.SnowMelt -= self.var.snowEvap
            self.var.IceMelt -= self.var.iceEvap

        #if dateVar['newStart'] or (dateVar['currDate'].day in [1,11,21]):
        #    self.var.cropKC[No] = readnetcdf2(coverType + '_cropCoefficientNC', dateVar['10day'], "10day")
        #    self.var.cropKC[No] = np.maximum(self.var.cropKC[No], self.var.minCropKC)
        #    self.var.cropKC_landCover[No] = self.var.cropKC[No].copy()

        # interpolation for each day from monthly values
        dplus = dateVar['30day'] + 1
        dpart = dateVar['doy'] % 30
        if dplus > 12:
            dplus = 0
        self.var.cropKC[No] = ((self.var.cropKCmonth[No, dplus, :] - self.var.cropKCmonth[No, dateVar['30day'], :]) / 30. * 
                               dpart + self.var.cropKCmonth[No, dateVar['30day'], :])
        if not hasattr(self.var, 'cropKC_landCover'):
            self.var.cropKC_landCover = np.zeros((4, len(globals.inZero)))
        self.var.cropKC_landCover[No] = self.var.cropKC[No].copy()
        cropKC_landCover = self.var.cropKC_landCover[No]


        if self.var.includeCrops:
            # includeCrops allows for crops and fallow land to makeup the landcovers grasslands and non-paddy, and
            # maintains including a representative vegetation. It is developed to allow users to decide on the crops
            # and parameters that are relevant for the study. The Excel cwatm_settings.xlsx is used to detail the crops
            # and associated parameters. Crops have a unique planting month and four growth stages. Each stage is associated with a
            # crop coefficient (Kc), yield response factor (Ky), and length.

            if No == 1:
                # Only go through this once:
                # I. new start and II. beginning of the month

                # I. new start
                if dateVar['newStart']:

                    for z in ['irrM3_Paddy_month_segment', 'irr_Paddy_month', 'irr_crop', 'irr_crop_month',
                              'irrM3_crop_month_segment', 'ratio_a_p_nonIrr', 'ratio_a_p_Irr', 'eta_nonIrr', 'eta_Irr',
                              'fracCrops_IrrLandDemand', 'fracCrops_Irr', 'areaCrops_Irr_segment', 'areaCrops_nonIrr_segment',
                              'fracCrops_nonIrrLandDemand', 'fracCrops_nonIrr', 'activatedCrops', 'monthCounter', 'monthCounter_nonIrr',
                              'currentKC', 'totalPotET_month', 'PET_cropIrr_m3', 'crop_pet_nonirr', 'crop_eta_nonirr',
                              'crop_ky_other_nonirr', 'crop_ky_other_irr', 'actTransTotal_month_Irr', 'actTransTotal_month_nonIrr',
                              'currentKY', 'Yield_Irr', 'crop_eta_irr', 'crop_pet_irr', 'crop_Tavg_nonIrr', 'crop_Tavg_Irr',
                              'ratio_et_nonIrr', 'ratio_et_Irr', 'crop_Precipitation_nonIrr', 'crop_Precipitation_Irr',
                              'crop_Tmin_nonIrr', 'crop_Tmin_Irr', 'crop_Tmax_nonIrr', 'crop_Tmax_Irr', 'Tavg_month',
                              'Precipitation_month', 'Tmin_month', 'Tmax_month', 'Yield_nonIrr', 'actTransTotal_crops_Irr',
                              'actTransTotal_crops_nonIrr', 'PotET_crop', 'PotETaverage_crop_segments', 'totalPotET_month_segment',
                              'ky_ratio_et_nonIrr', 'ky_ratio_et_Irr', 'crop_daily_Precipitation_nonIrr',
                              'crop_daily_Precipitation_Irr', 'crop_daily_tavg_nonIrr', 'crop_daily_tavg_Irr', 'PotET_crop_Irr',
                              'PotET_crop_nonIrr', 'daily_yield_nonIrr', 'daily_yield_Irr', 'crop_pe_Irr', 'crop_pe_nonIrr',
                              'Yield_loss_Irr', 'Yield_loss_nonIrr', 'crop_soil_nonIrr', 'crop_soil_Irr', 'crop_max_soil_Irr',
                              'crop_max_soil_nonIrr', 'rwss', 'gs_eta_Irr', 'gs_eta_nonIrr', 'gs_et0_Irr', 'gs_et0_nonIrr',
                              'difference_kc_gs', 'difference_ky_gs', 'gs_ratio_Irr', 'gs_ratio_nonIrr', 'gs_loss_Irr',
                              'gs_loss_nonIrr', 'gs_1_Irr', 'gs_2_Irr', 'gs_3_Irr', 'gs_4_Irr', 'gs_5_Irr', 'gs_1_nonIrr',
                              'gs_2_nonIrr', 'gs_3_nonIrr', 'gs_4_nonIrr', 'gs_5_nonIrr', 'drystress_Irr', 'wetstress_Irr','rwss',
                              'drystress_nonIrr', 'wetstress_nonIrr', 'gs_soil_nonIrr', 'gs_soil_Irr', 'gs_max_soil_Irr',
                              'gs_max_soil_nonIrr', 'gs_pe_Irr', 'gs_pe_nonIrr', 'ET_crop_nonIrr', 'ET_crop_Irr',
                              'ratio_a_p_nonIrr_daily', 'ratio_a_p_Irr_daily']:
                        vars(self.var)[z] = np.tile(globals.inZero, (len(self.var.Crops), 1))

                    self.var.irr_Paddy_month = globals.inZero.copy()
                    for z in [crop for crop in self.var.Crops_names]:
                        vars(self.var)[z + '_Irr'] = globals.inZero.copy()
                        vars(self.var)[z + '_nonIrr'] = globals.inZero.copy()

                    self.var.ET_crop_Irr_paddy = globals.inZero.copy()
                    self.var.ET_crop_Irr_paddy_fraccrop = globals.inZero.copy()

                    for c in range(len(self.var.Crops)):
                        vars(self.var)['ET_crop_Irr_' + str(c)] = globals.inZero.copy()
                        vars(self.var)['ET_crop_Irr_fraccrop_' + str(c)] = globals.inZero.copy()
                        vars(self.var)['ET_crop_nonIrr_' + str(c)] = globals.inZero.copy()
                        vars(self.var)['ET_crop_nonIrr_fraccrop_' + str(c)] = globals.inZero.copy()
                        vars(self.var)['irr_crop_' + str(c)] = globals.inZero.copy()

                        self.var.activatedCrops[c] = self.var.load_initial("activatedCrops_" + str(c))
                        self.var.fracCrops_Irr[c] = self.var.load_initial('fracCrops_Irr_' + str(c))
                        self.var.fracCrops_nonIrr[c] = self.var.load_initial('fracCrops_nonIrr_' + str(c))
                        self.var.monthCounter[c] = self.var.load_initial("monthCounter_" + str(c))

                if dateVar['newStart'] or dateVar['newYear']:

                    if 'co2conc' in binding:
                        self.var.co2 = readnetcdf2('co2conc', dateVar['currDate'], "yearly", value="CO2", cut=False, compress=False)

                    # NOTE: choose 0.01 when loadmap values are in percent (0-100);
                    # choose 1 when they are already fractions (0-1).
                    crop_inflate_factor = 0.01

                    for i in range(len(self.var.Crops)):
                        self.var.gs_1_nonIrr[i] = readnetcdfWithoutTime(cbinding(self.var.Crops_names[i] + '_nonIrr_calendar'), 'planting_day').astype(np.int16)
                        self.var.gs_1_Irr[i] = readnetcdfWithoutTime(cbinding(self.var.Crops_names[i] + '_Irr_calendar'), 'planting_day').astype(np.int16)
                        self.var.gs_2_nonIrr[i] = readnetcdfWithoutTime(cbinding(self.var.Crops_names[i] + '_nonIrr_calendar'), 'node1').astype(np.int16)
                        self.var.gs_2_Irr[i] = readnetcdfWithoutTime(cbinding(self.var.Crops_names[i] + '_Irr_calendar'), 'node1').astype(np.int16)
                        self.var.gs_3_nonIrr[i] = readnetcdfWithoutTime(cbinding(self.var.Crops_names[i] + '_nonIrr_calendar'), 'node2').astype(np.int16)
                        self.var.gs_3_Irr[i] = readnetcdfWithoutTime(cbinding(self.var.Crops_names[i] + '_Irr_calendar'), 'node2').astype(np.int16)
                        self.var.gs_4_nonIrr[i] = readnetcdfWithoutTime(cbinding(self.var.Crops_names[i] + '_nonIrr_calendar'), 'node3').astype(np.int16)
                        self.var.gs_4_Irr[i] = readnetcdfWithoutTime(cbinding(self.var.Crops_names[i] + '_Irr_calendar'), 'node3').astype(np.int16)
                        self.var.gs_5_nonIrr[i] = readnetcdfWithoutTime(cbinding(self.var.Crops_names[i] + '_nonIrr_calendar'), 'node4').astype(np.int16)
                        self.var.gs_5_Irr[i] = readnetcdfWithoutTime(cbinding(self.var.Crops_names[i] + '_Irr_calendar'), 'node4').astype(np.int16)

                        try:
                            self.var.fracCrops_IrrLandDemand[i] = np.where(
                                (loadmap(self.var.Crops_names[i] + '_Irr') * crop_inflate_factor) <= 1,
                                loadmap(self.var.Crops_names[i] + '_Irr') * crop_inflate_factor, 1)
                            self.var.fracCrops_nonIrrLandDemand[i] = np.where(
                                (loadmap(self.var.Crops_names[i] + '_nonIrr') * crop_inflate_factor) <= 1,
                                loadmap(self.var.Crops_names[i] + '_nonIrr') * crop_inflate_factor,
                                1)
                        except:
                            self.var.fracCrops_IrrLandDemand[i] = readnetcdf2(
                                self.var.Crops_names[i] + '_Irr', dateVar['currDate'], 'yearly',
                                value=re.split(r'[^a-zA-Z0-9_[\]]', cbinding(self.var.Crops_names[i] + '_Irr'))[-2])
                            self.var.fracCrops_nonIrrLandDemand[i] = readnetcdf2(
                                self.var.Crops_names[i] + '_nonIrr', dateVar['currDate'], 'yearly',
                                value=re.split(r'[^a-zA-Z0-9_[\]]', cbinding(self.var.Crops_names[i] + '_nonIrr'))[-2])

                        if 'crops_leftoverNotIrrigated' in binding:
                            if i <= int(cbinding('crops_leftoverNotIrrigated')):
                                self.var.fracCrops_nonIrrLandDemand[i] = globals.inZero.copy()

                        self.var.activatedCrops[i] = np.minimum(
                            np.maximum((self.var.fracCrops_IrrLandDemand[i] + self.var.fracCrops_nonIrrLandDemand[i] + 0.99999) // 1,
                                       self.var.activatedCrops[i]), 1)

                if dateVar['currDate'].day >= 1:

                    doy = dateVar['currDate'].timetuple().tm_yday
                    n_crops = len(self.var.Crops)

                    allow_leftover = 'leftoverIrrigatedCropIsRainfed' in option and checkOption('leftoverIrrigatedCropIsRainfed')
                    forbid_nonirr = 'crops_leftoverNotIrrigated' in binding
                    forbid_idx = int(cbinding('crops_leftoverNotIrrigated')) if forbid_nonirr else -1

                    if 'moveIrrFallowToNonIrr' in option and checkOption('moveIrrFallowToNonIrr'):
                        self.var.fracVegCover[3] = self.var.irrigatedArea_original.copy()
                        self.var.fracVegCover[1] = np.maximum(
                            1.0 - self.var.fracVegCover[0] - self.var.fracVegCover[2] - self.var.fracVegCover[3]
                            - self.var.fracVegCover[4] - self.var.fracVegCover[5], 0.0)

                    def sum_arrays(arrs):
                        out = globals.inZero.copy()
                        for a in arrs:
                            out += a
                        return out

                    # 1) advance crop status once per day
                    for c in range(n_crops):
                        mc = self.var.monthCounter[c]
                        mc = np.where(mc > 0, mc + 1, mc)
                        mc = np.where(doy == self.var.gs_1_Irr[c], doy, mc)
                        mc = np.where(mc >= self.var.gs_5_Irr[c], 0, mc)
                        self.var.monthCounter[c] = mc

                    # 2) build active land demand
                    demand_irr, demand_nonirr = [], []
                    for c in range(n_crops):
                        active = self.var.monthCounter[c] > 0
                        irr = np.where(active, self.var.fracCrops_IrrLandDemand[c], 0.0)
                        nonirr = np.where(active, self.var.fracCrops_nonIrrLandDemand[c], 0.0)
                        if forbid_nonirr and c <= forbid_idx:
                            nonirr = globals.inZero.copy()
                        demand_irr.append(irr)
                        demand_nonirr.append(nonirr)

                    # 3) scale irrigated demand into irrigated land class
                    total_irr_demand = sum_arrays(demand_irr)
                    scale_irr = np.where(
                        total_irr_demand > self.var.fracVegCover[3],
                        np.where(total_irr_demand > 0, self.var.fracVegCover[3] / total_irr_demand, 0.0),
                        1.0)
                    for c in range(n_crops):
                        self.var.fracCrops_Irr[c] = demand_irr[c] * scale_irr

                    # leftover irrigated demand may spill into rainfed area
                    if allow_leftover:
                        for c in range(n_crops):
                            demand_nonirr[c] += np.maximum(demand_irr[c] - self.var.fracCrops_Irr[c], 0.0)
                            if forbid_nonirr and c <= forbid_idx:
                                demand_nonirr[c] = globals.inZero.copy()

                    # 4) scale non-irrigated demand into grassland class
                    total_nonirr_demand = sum_arrays(demand_nonirr)
                    scale_nonirr = np.where(
                        total_nonirr_demand > self.var.fracVegCover[1],
                        np.where(total_nonirr_demand > 0, self.var.fracVegCover[1] / total_nonirr_demand, 0.0),
                        1.0)
                    for c in range(n_crops):
                        self.var.fracCrops_nonIrr[c] = demand_nonirr[c] * scale_nonirr

                    self.var.frac_totalIrr = sum_arrays(self.var.fracCrops_Irr)
                    self.var.frac_totalnonIrr = sum_arrays(self.var.fracCrops_nonIrr)
                    self.var.frac_totalIrr_max = np.maximum(self.var.frac_totalIrr, self.var.frac_totalIrr_max)
                    self.var.frac_totalnonIrr_max = np.maximum(self.var.frac_totalnonIrr, self.var.frac_totalnonIrr_max)

                    crop_corr_nonirr = self.var.crop_correct_landCover[1] if hasattr(self.var, 'crop_correct_landCover') else 1.0
                    crop_corr_irr = self.var.crop_correct_landCover[3] if hasattr(self.var, 'crop_correct_landCover') else 1.0

                    # 5) update current stage-wise Kc/Ky and growth-stage diagnostics
                    for c in range(n_crops):
                        kc_prev = self.var.currentKC[c].copy()
                        ky_prev = self.var.currentKY[c].copy()

                        stages = [self.var.gs_1_Irr[c], self.var.gs_2_Irr[c], self.var.gs_3_Irr[c],
                                  self.var.gs_4_Irr[c], self.var.gs_5_Irr[c]]
                        kc_stages = [self.var.Crops[c][1][1], self.var.Crops[c][2][1],
                                     self.var.Crops[c][3][1], self.var.Crops[c][4][1]]
                        ky_stages = [self.var.Crops[c][1][2], self.var.Crops[c][2][2],
                                     self.var.Crops[c][3][2], self.var.Crops[c][4][2]]

                        kc_prev = self.var.currentKC[c].copy()   # 保存上一时刻的值

                        for i in range(4):
                            in_stage = ((self.var.monthCounter[c] >= stages[i]) &
                                        (self.var.monthCounter[c] < stages[i + 1]))
                            self.var.currentKC[c] = np.where(in_stage, kc_stages[i], self.var.currentKC[c])
                            self.var.currentKY[c] = np.where(in_stage, ky_stages[i], self.var.currentKY[c])

                        # ==================== 新增：阶段切换时平滑过渡 ====================
                        stage_changed = (kc_prev != self.var.currentKC[c])
                        self.var.currentKC[c] = np.where(
                            stage_changed,
                            0.8 * kc_prev + 0.2 * self.var.currentKC[c],   # 0.8/0.2 很平稳，你可以改成 0.7/0.3 或 0.9/0.1
                            self.var.currentKC[c]
                        )

                        self.var.currentKC[c] = np.where(self.var.monthCounter[c] > 0, self.var.currentKC[c], 0.0)
                        self.var.currentKY[c] = np.where(self.var.monthCounter[c] > 0, self.var.currentKY[c], 0.0)
                        # dbg_idx = 33367
                        # print("\n--- CROP KC TRACE ---")
                        # print("date =", dateVar['currDate'], "idx =", dbg_idx)
                        # print("cropKC_paddy =", float(self.var.cropKC[2][dbg_idx]))
                        # if hasattr(self.var, "monthCounter"):
                        #     print("monthCounter_paddy =", self.var.monthCounter[0][dbg_idx] if hasattr(self.var.monthCounter[0], "__getitem__") else self.var.monthCounter[0])
                        # print("---------------------\n")
                        stage_changed = (kc_prev != self.var.currentKC[c])

                        self.var.gs_ratio_nonIrr[c] = np.where(
                            stage_changed,
                            np.where(self.var.gs_et0_nonIrr[c] > 0, np.maximum(1 - self.var.gs_eta_nonIrr[c] / self.var.gs_et0_nonIrr[c], 0.0), 0.0),
                            0.0)
                        self.var.gs_loss_nonIrr[c] = np.where(stage_changed, np.minimum(ky_prev * self.var.gs_ratio_nonIrr[c], 1.0), 0.0)
                        self.var.gs_eta_nonIrr[c] = np.where(stage_changed, 0.0, self.var.gs_eta_nonIrr[c] + self.var.eta_nonIrr[c])
                        self.var.gs_et0_nonIrr[c] = np.where(stage_changed, 0.0, self.var.gs_et0_nonIrr[c] + self.var.crop_pet_nonirr[c])

                        self.var.gs_ratio_Irr[c] = np.where(
                            stage_changed,
                            np.where(self.var.gs_et0_Irr[c] > 0, np.maximum(1 - self.var.gs_eta_Irr[c] / self.var.gs_et0_Irr[c], 0.0), 0.0),
                            0.0)
                        self.var.gs_loss_Irr[c] = np.where(stage_changed, np.minimum(ky_prev * self.var.gs_ratio_Irr[c], 1.0), 0.0)
                        self.var.gs_eta_Irr[c] = np.where(stage_changed, 0.0, self.var.gs_eta_Irr[c] + self.var.eta_Irr[c])
                        self.var.gs_et0_Irr[c] = np.where(stage_changed, 0.0, self.var.gs_et0_Irr[c] + self.var.crop_pet_irr[c])

                        self.var.gs_soil_nonIrr[c] = np.where(stage_changed, 0.0, self.var.gs_soil_nonIrr[c] + self.var.crop_soil_nonIrr[c])
                        self.var.gs_soil_Irr[c] = np.where(stage_changed, 0.0, self.var.gs_soil_Irr[c] + self.var.crop_soil_Irr[c])
                        self.var.gs_max_soil_Irr[c] = np.where(stage_changed, 0.0, self.var.gs_max_soil_Irr[c] + self.var.crop_max_soil_Irr[c])
                        self.var.gs_max_soil_nonIrr[c] = np.where(stage_changed, 0.0, self.var.gs_max_soil_nonIrr[c] + self.var.crop_max_soil_nonIrr[c])
                        self.var.gs_pe_Irr[c] = np.where(stage_changed, 0.0, self.var.gs_pe_Irr[c] + self.var.crop_pe_Irr[c])
                        self.var.gs_pe_nonIrr[c] = np.where(stage_changed, 0.0, self.var.gs_pe_nonIrr[c] + self.var.crop_pe_nonIrr[c])

                        self.var.crop_pet_nonirr[c] = np.where(
                            self.var.monthCounter[c] > 0,
                            self.var.cropCorrect * crop_corr_nonirr * self.var.ETRef * self.var.fracCrops_nonIrr[c] * self.var.currentKC[c],
                            0.0)
                        self.var.crop_pet_irr[c] = np.where(
                            self.var.monthCounter[c] > 0,
                            self.var.cropCorrect * crop_corr_irr * self.var.ETRef * self.var.fracCrops_Irr[c] * self.var.currentKC[c],
                            0.0)

                        self.var.ratio_a_p_nonIrr[c] = np.where(
                            self.var.crop_pet_nonirr[c] > 0, self.var.eta_nonIrr[c] / self.var.crop_pet_nonirr[c], 0.0)
                        self.var.ratio_a_p_Irr[c] = np.where(
                            self.var.crop_pet_irr[c] > 0, self.var.eta_Irr[c] / self.var.crop_pet_irr[c], 0.0)
                        self.var.ratio_a_p_nonIrr_daily[c] = self.var.ratio_a_p_nonIrr[c].copy()
                        self.var.ratio_a_p_Irr_daily[c] = self.var.ratio_a_p_Irr[c].copy()

                        self.var.Yield_loss_nonIrr[c] = np.where(
                            self.var.monthCounter[c] > 0, np.minimum(self.var.currentKY[c] * (1 - self.var.ratio_a_p_nonIrr[c]), 1.0), 0.0)
                        self.var.Yield_loss_Irr[c] = np.where(
                            self.var.monthCounter[c] > 0, np.minimum(self.var.currentKY[c] * (1 - self.var.ratio_a_p_Irr[c]), 1.0), 0.0)
                        self.var.Yield_nonIrr[c] = np.where(
                            self.var.monthCounter[c] > 0, np.maximum(1 - self.var.Yield_loss_nonIrr[c], 0.0), 0.0)
                        self.var.Yield_Irr[c] = np.where(
                            self.var.monthCounter[c] > 0, np.maximum(1 - self.var.Yield_loss_Irr[c], 0.0), 0.0)

                    # 6) general crop, fallow and available land
                    self.var.generalIrrCrop_max = np.maximum(self.var.fracVegCover[3] - self.var.frac_totalIrr_max, globals.inZero.copy())
                    self.var.generalnonIrrCrop_max = np.maximum(self.var.fracVegCover[1] - self.var.frac_totalnonIrr_max, globals.inZero.copy())

                    if 'GeneralCrop_Irr' in binding and checkOption('use_GeneralCropIrr') is True:
                        self.var.GeneralCrop_Irr = loadmap('GeneralCrop_Irr')
                        self.var.GeneralCrop_Irr = np.minimum(np.abs(self.var.fracVegCover[3] - self.var.frac_totalIrr), self.var.GeneralCrop_Irr)
                    elif checkOption('use_GeneralCropIrr') is False:
                        if checkOption('activate_fallow') is True:
                            if checkOption('automaticFallowingIrr') is True:
                                self.var.GeneralCrop_Irr = self.var.generalIrrCrop_max.copy()
                            else:
                                self.var.GeneralCrop_Irr = globals.inZero.copy()
                        else:
                            self.var.GeneralCrop_Irr = self.var.fracVegCover[3] - self.var.frac_totalIrr

                    self.var.fallowIrr = np.maximum(0.0, self.var.fracVegCover[3] - (self.var.frac_totalIrr + self.var.GeneralCrop_Irr))
                    self.var.fallowIrr_max = np.maximum(self.var.fallowIrr, self.var.fallowIrr_max)

                    if 'moveIrrFallowToNonIrr' in option and checkOption('moveIrrFallowToNonIrr'):
                        self.var.fracVegCover[3] = self.var.frac_totalIrr + self.var.GeneralCrop_Irr
                        self.var.fracVegCover[1] = np.maximum(
                            1.0 - self.var.fracVegCover[0] - self.var.fracVegCover[2] - self.var.fracVegCover[3]
                            - self.var.fracVegCover[4] - self.var.fracVegCover[5], 0.0)

                    if 'GeneralCrop_nonIrr' in binding and checkOption('use_GeneralCropnonIrr') is True:
                        self.var.GeneralCrop_nonIrr = loadmap('GeneralCrop_nonIrr')
                        self.var.GeneralCrop_nonIrr = np.minimum(np.abs(self.var.fracVegCover[1] - self.var.frac_totalnonIrr), self.var.GeneralCrop_nonIrr)
                    elif checkOption('use_GeneralCropnonIrr') is False:
                        if checkOption('activate_fallow') is True:
                            self.var.GeneralCrop_nonIrr = self.var.generalnonIrrCrop_max.copy()
                        else:
                            self.var.GeneralCrop_nonIrr = self.var.fracVegCover[1] - self.var.frac_totalnonIrr

                    self.var.fallownonIrr = np.maximum(0.0, self.var.fracVegCover[1] - (self.var.frac_totalnonIrr + self.var.GeneralCrop_nonIrr))
                    self.var.fallownonIrr_max = np.maximum(self.var.fallownonIrr, self.var.fallownonIrr_max)
                    self.var.availableArableLand = self.var.fallowIrr + self.var.fracVegCover[1] - self.var.frac_totalnonIrr

            if No == 1:
                n_crops = len(self.var.Crops)

                avail_nonirr = np.maximum(self.var.fracVegCover[1] - self.var.frac_totalnonIrr, 0.0)
                self.var.GeneralCrop_nonIrr = np.minimum(self.var.GeneralCrop_nonIrr, avail_nonirr)
                self.var.GeneralCrop_nonIrr = np.maximum(self.var.GeneralCrop_nonIrr, 0.0)
                self.var.fallownonIrr = np.maximum(
                    self.var.fracVegCover[1] - self.var.frac_totalnonIrr - self.var.GeneralCrop_nonIrr, 0.0)

                weighted_general_nonirr = self.var.GeneralCrop_nonIrr * self.var.cropKC_landCover[1]
                weighted_crop_nonirr = globals.inZero.copy()
                for c in range(n_crops):
                    weighted_crop_nonirr += np.where(
                        self.var.monthCounter[c] > 0,
                        self.var.fracCrops_nonIrr[c] * self.var.currentKC[c],
                        0.0)
                weighted_fallow_nonirr = self.var.fallownonIrr * self.var.minCropKC
                weighted_total_nonirr = weighted_general_nonirr + weighted_crop_nonirr + weighted_fallow_nonirr

                self.var.weighted_KC_nonIrr_woFallow_fullKc = (weighted_general_nonirr + weighted_crop_nonirr).copy()

                weighted_general_nonirr_excess = self.var.GeneralCrop_nonIrr * np.maximum(
                    self.var.cropKC_landCover[1] - self.var.minCropKC, 0.0)
                weighted_crop_nonirr_excess = globals.inZero.copy()
                for c in range(n_crops):
                    weighted_crop_nonirr_excess += np.where(
                        self.var.monthCounter[c] > 0,
                        self.var.fracCrops_nonIrr[c] * np.maximum(self.var.currentKC[c] - self.var.minCropKC, 0.0),
                        0.0)
                self.var.weighted_KC_nonIrr_woFallow = (
                    weighted_general_nonirr_excess + weighted_crop_nonirr_excess).copy()

                den_nonirr = np.where(self.var.fracVegCover[1] > 1e-10, self.var.fracVegCover[1], np.nan)
                self.var.cropKC[1] = np.where(np.isfinite(den_nonirr), weighted_total_nonirr / den_nonirr, 0.0)
                self.var.cropKC[1] = np.clip(self.var.cropKC[1], self.var.minCropKC, 2.0)

            if No == 3:
                n_crops = len(self.var.Crops)

                avail_irr = np.maximum(self.var.fracVegCover[3] - self.var.frac_totalIrr, 0.0)
                self.var.GeneralCrop_Irr = np.minimum(self.var.GeneralCrop_Irr, avail_irr)
                self.var.GeneralCrop_Irr = np.maximum(self.var.GeneralCrop_Irr, 0.0)
                self.var.fallowIrr = np.maximum(
                    self.var.fracVegCover[3] - self.var.frac_totalIrr - self.var.GeneralCrop_Irr, 0.0)

                weighted_general_irr = self.var.GeneralCrop_Irr * self.var.cropKC_landCover[3]
                weighted_crop_irr = globals.inZero.copy()
                for c in range(n_crops):
                    if self.var.Crops_names[c][:4] != 'rice':
                        weighted_crop_irr += np.where(
                            self.var.monthCounter[c] > 0,
                            self.var.fracCrops_Irr[c] * self.var.currentKC[c],
                            0.0)
                weighted_fallow_irr = self.var.fallowIrr * self.var.minCropKC
                weighted_total_irr = weighted_general_irr + weighted_crop_irr + weighted_fallow_irr

                self.var.weighted_KC_Irr_woFallow_fullKc = (weighted_general_irr + weighted_crop_irr).copy()

                weighted_general_irr_excess = self.var.GeneralCrop_Irr * np.maximum(
                    self.var.cropKC_landCover[3] - self.var.minCropKC, 0.0)
                weighted_crop_irr_excess = globals.inZero.copy()
                for c in range(n_crops):
                    if self.var.Crops_names[c][:4] != 'rice':
                        weighted_crop_irr_excess += np.where(
                            self.var.monthCounter[c] > 0,
                            self.var.fracCrops_Irr[c] * np.maximum(self.var.currentKC[c] - self.var.minCropKC, 0.0),
                            0.0)
                self.var._weighted_KC_Irr = (weighted_general_irr_excess + weighted_crop_irr_excess).copy()
                self.var.weighted_KC_Irr_woFallow = self.var._weighted_KC_Irr.copy()

                den_irr = np.where(self.var.fracVegCover[3] > 1e-10, self.var.fracVegCover[3], np.nan)
                self.var.cropKC[3] = np.where(np.isfinite(den_irr), weighted_total_irr / den_irr, 0.0)
                self.var.cropKC[3] = np.clip(self.var.cropKC[3], self.var.minCropKC, 2.0)

        # without crops
        # calculate potential ET
        ##  self.var.totalPotET total potential evapotranspiration for a reference crop for a land cover class [m]
        self.var.totalPotET[No] = self.var.cropCorrect * self.var.crop_correct_landCover[No] * self.var.cropKC[No] * self.var.ETRef


        # calculate transpiration


        # potTranspiration: Transpiration for each land cover class
        # Dealt with above - self.var.snowEvap
        self.var.potTranspiration[No] = np.maximum(0., self.var.totalPotET[No] - self.var.potBareSoilEvap)
        if self.var.includeCrops and No in [1, 3]:          # 只对 grassland 和 irrNonPaddy 生效
            # 估算根区可用水（你可以以后再改得更精确）
            root_avail = np.maximum(0.,
                (self.var.w1[No] - self.var.wwp1[No]) * 0.3 +
                (self.var.w2[No] - self.var.wwp2[No]) * 0.4 +
                (self.var.w3[No] - self.var.wwp3[No]) * 0.3
            )
            pot_trans_max = root_avail * 1.25                # 允许最多 25% 超采作为缓冲

            self.var.potTranspiration[No] = np.minimum(
                self.var.potTranspiration[No], pot_trans_max
            )
        # checkOption('includeCrops') and checkOption('includeCropSpecificWaterUse')
        if self.var.includeCrops:

            # only goes through ones
            if No == 3:

                for c in range(len(self.var.Crops)):

                    self.var.PotET_crop[c] = (self.var.cropCorrect * self.var.crop_correct_landCover[No] * 
                                              self.var.currentKC[c] * self.var.ETRef)
                    self.var.totalPotET_month[c] += self.var.PotET_crop[c]
                    # self.var.cropCorrect * self.var.currentKC[c] * self.var.ETRef
                    # np.maximum(0., self.var.cropCorrect * self.var.currentKC[c] * self.var.ETRef - 
                    #           self.var.potBareSoilEvap - self.var.snowEvap)

                    # For creating named crop maps
                    vars(self.var)[self.var.Crops_names[c] + '_Irr'] = self.var.fracCrops_Irr[c].copy()
                    vars(self.var)[self.var.Crops_names[c] + '_nonIrr'] = self.var.fracCrops_nonIrr[c].copy()

                    

                    if 'adminSegments' in binding:
                        self.var.totalPotET_month_segment[c] = npareaaverage(self.var.totalPotET_month[c], self.var.adminSegments)
                        self.var.PotETaverage_crop_segments[c] = npareaaverage(self.var.PotET_crop[c], self.var.adminSegments)

                        self.var.areaCrops_Irr_segment[c] = npareatotal(self.var.fracCrops_Irr[c] * self.var.cellArea,
                                                                        self.var.adminSegments)

                        self.var.areaCrops_nonIrr_segment[c] = npareatotal(
                            self.var.fracCrops_nonIrr[c] * self.var.cellArea,
                            self.var.adminSegments)


                if 'adminSegments' in binding:
                    self.var.areaPaddy_Irr_segment = npareatotal(self.var.fracVegCover[2] * self.var.cellArea,
                                                             self.var.adminSegments)

                    self.var.Precipitation_segment = npareatotal(self.var.Precipitation * self.var.cellArea,
                                                                 self.var.adminSegments)

                    self.var.availableArableLand_segment = npareatotal(self.var.availableArableLand * self.var.cellArea,
                                                                        self.var.adminSegments)

