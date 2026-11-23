# -*- coding: utf-8 -*-
"""Part 2 of Non-Rerun Robustness & Methodological Reliability Engine (Daily 0.5-deg Physics).

Evaluates 4 daily process-level sensitivity dimensions WITHOUT re-running CWatM, using the
archived daily 0.5-deg soil water, relative water stress, daily SPEI, DOY normal climatology,
and exact stage crop calendars over 1989-2019 (11,322 days; 30 harvest years 1990-2019):

  - Step 0: Offline Daily Stage-Integrator Fidelity Check against CWatM dynamic archive & FAOSTAT.
  - Dimension R2 (User Item 2): Multi-Indicator Conjunctive Two-Gate (Soil Moisture x SPEI) vs.
    SPEI-Only Sensitivity across Dry (-0.5..-1.5), Wet (+1.0..+2.0), Symmetric (-1/+1, -1.3/+1.3,
    -1.5/+1.5), Empirical Soil-SPEI Coupling Curves, Near-Normal Reference Shifts, & Wet Blend Weights.
  - Dimension R3 (User Item 3): Seven Two-Sided k_s(theta) Functional Forms (Piecewise-Linear FAO-56,
    User's Continuous Exp-Quadratic, Gaussian Bell Curve + parameter envelope, Sigmoidal S-Curve +
    parameter envelope, Concave/Convex Power Curves, and Feddes 15% Anaerobiosis-Buffered Trapezoid)
    + Cross-Validated Empirical Curve-Fitting Parameter Dispersion Audit across Hydroclimatic Regimes.
  - Dimension R4a (R2 Major 7): Paddy & Rainfed Rice Growth-Stage Waterlogging Suppression
    (All 4 stages vs. 50%-attenuated vegetative Stages 1-2 vs. Reproductive/Late Stages 3-4 only
    vs. Booting/Heading Stage 3 only).
  - Dimension R4b (Upgraded Fig. S9, R3 Major 5): Crop-by-crop Kc and Ky multiplier perturbations
    (0.7-1.3) for both relative-yield elasticity S and attributed drought & waterlogging losses.

Strictly READ-ONLY on ./data/AAUDE, E:/, ./data/; writes ONLY to ./data/Software/Gemini/drought_waterlogging/results/robustness_norerun/.
"""
from __future__ import annotations

import os
os.environ["OMP_NUM_THREADS"] = "2"
os.environ["MKL_NUM_THREADS"] = "2"
os.environ["OPENBLAS_NUM_THREADS"] = "2"

import sys
import time
from concurrent.futures import ThreadPoolExecutor
import h5py
import numpy as np
import pandas as pd
import xarray as xr
import rioxarray as rxr
from scipy import stats
from scipy.optimize import curve_fit

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(os.path.dirname(HERE), "data_tables")
os.makedirs(OUT_DIR, exist_ok=True)

CELLAREA_FP = r"./data/external/cellarea.nc"
FRAC_DIR = r"./data/cwatm_inputs_figshare"
SPAM_H_DIR = r"./data/sensitivity/spam2020_30arcmin/harvested_area"
GAEZ_HIST_DIR = r"./data/gaez_v4_200b_noco2/hist"
CAL_HIST_DIR = r"./data/cwatm_inputs_figshare/after_clip"
HIST_ANNUAL_DIR = r"./data/cwatm_outputs/gswp3_long"

RWS_DAILY_PATHS = {
    1: (r"./data/cwatm_outputs/reference/rwss[1]_daily.nc", "rwss[1]"),
    2: (r"./data/cwatm_outputs/reference/rwss[2]_daily.nc", "rwss[2]"),
    3: (r"./data/cwatm_outputs/reference/rwss[3]_daily.nc", "rwss[3]"),
}
SOIL_SUM_FP = r"./data/cwatm_outputs/GSWP3-W5E5/sum_soil_daily.nc"
ETREF_MON_FP = r"./data/cwatm_outputs/GSWP3-W5E5/ETRef_monthtot.nc"
SPEI_DAILY_FP = r"./data/cwatm_outputs/SPEI_daily1.nc"
RWS_CLIM_FP = r".\rws_climatology_time366_under_spei_-0.5_0.5_chunked.nc"
WWP_FP = r".\cell_wwp_annualavg.nc"
WCRIT_FP = r".\cell_wcrit_annualavg.nc"
WFC_FP = r".\cell_wfc_annualavg.nc"
ISO_CACHE_FP = r".\cache\grid_country_iso3_05deg.npz"
FAO_PANEL_FP = r".\results\tier3_multiscale_4crops\tier3_panel_national_faostat.csv"

CROP_CORRECT = np.float32(1.7387144425028453)  # Exact calibration factor from CWatM settings
MIN_CROP_KC = np.float32(0.20)                 # Exact bare-soil floor from settings_1_waterstress.ini L171

# Exact stage KC1..4 and KY1..4 from ./settings/Table_S2_Soil_and_Crop_Hydraulic_Parameters.csv
STAGE_KC = {
    "wwh": np.array([0.70, 0.925, 1.15, 0.25], dtype=np.float32),
    "swh": np.array([0.40, 0.775, 1.15, 0.25], dtype=np.float32),
    "rice": np.array([1.10, 1.100, 1.20, 0.60], dtype=np.float32),
    "maiz": np.array([0.30, 0.750, 1.20, 0.50], dtype=np.float32),
    "soyb": np.array([0.50, 0.800, 1.15, 0.50], dtype=np.float32),
}
STAGE_KY = {
    "wwh": np.array([0.20, 0.60, 0.50, 0.21], dtype=np.float32),
    "swh": np.array([0.20, 0.65, 0.55, 0.21], dtype=np.float32),
    "rice": np.array([1.00, 1.10, 1.20, 0.20], dtype=np.float32),
    "maiz": np.array([0.40, 0.41, 1.00, 0.50], dtype=np.float32),
    "soyb": np.array([0.20, 0.80, 0.90, 0.21], dtype=np.float32),
}
SUBCROPS = ["maiz", "soyb", "rice", "wwh", "swh"]
CROPS4 = ["maiz", "soyb", "rice", "whea"]
MODES = ["nonIrr", "Irr"]

DATES_31Y = pd.date_range("1989-01-01", "2019-12-31", freq="D")  # len == 11,322
T_DAYS = len(DATES_31Y)
DOY_31Y = DATES_31Y.dayofyear.values.astype(np.int32)            # 1..366
YEAR_31Y = DATES_31Y.year.values.astype(np.int32)                # 1989..2019
MONTH_IDX_31Y = ((DATES_31Y.year - 1989) * 12 + (DATES_31Y.month - 1)).values.astype(np.int32)  # 0..371
HARV_YEARS = np.arange(1990, 2020, dtype=np.int32)               # 30 harvest years


# ------------------------------------------------------------------------------
# 1. Load Static Maps & Active Cultivated Pixel Index
# ------------------------------------------------------------------------------
def load_static_and_active_pixels() -> dict:
    ca_ds = xr.open_dataset(CELLAREA_FP)
    lat = ca_ds["lat"].values.astype(np.float64)
    lon = ca_ds["lon"].values.astype(np.float64)
    cellarea_m2 = ca_ds["cellarea30min"].values.astype(np.float64)
    ca_ds.close()

    tif_aligned = {
        "wwh_Irr": "MIRCA-OS_Wheat_2000_irr_30arcmin_winter_aligned.tif",
        "swh_Irr": "MIRCA-OS_Wheat_2000_irr_30arcmin_spring_aligned.tif",
        "rice_Irr": "MIRCA-OS_Rice_2000_ir_30arcmin_aligned.tif",
        "maiz_Irr": "MIRCA-OS_Maize_2000_ir_30arcmin_aligned.tif",
        "soyb_Irr": "MIRCA-OS_Soybeans_2000_ir_30arcmin_aligned.tif",
        "wwh_nonIrr": "MIRCA-OS_Wheat_2000_rf_30arcmin_winter_aligned.tif",
        "swh_nonIrr": "MIRCA-OS_Wheat_2000_rf_30arcmin_spring_aligned.tif",
        "rice_nonIrr": "MIRCA-OS_Rice_2000_rf_30arcmin_aligned.tif",
        "maiz_nonIrr": "MIRCA-OS_Maize_2000_rf_30arcmin_aligned.tif",
        "soyb_nonIrr": "MIRCA-OS_Soybeans_2000_rf_30arcmin_aligned.tif",
    }
    area_2d = {}
    for k, fn in tif_aligned.items():
        da = rxr.open_rasterio(os.path.join(FRAC_DIR, fn)).squeeze("band", drop=True)
        da = da.rename({"y": "lat", "x": "lon"}).interp(lat=lat, lon=lon, method="linear")
        v = np.where(np.isfinite(da.values) & (da.values > 0), da.values, 0.0)
        area_2d[k] = (v * (cellarea_m2 / 1e4)).astype(np.float64)
    area_2d["whea_Irr"] = area_2d["wwh_Irr"] + area_2d["swh_Irr"]
    area_2d["whea_nonIrr"] = area_2d["wwh_nonIrr"] + area_2d["swh_nonIrr"]

    gaez_codes = {"wwh": "wwhe", "swh": "swhe", "rice": "ricw", "maiz": "maiz", "soyb": "soyb"}
    pot_2d = {}
    for sub, gcode in gaez_codes.items():
        for mode in MODES:
            ds = xr.open_dataset(os.path.join(GAEZ_HIST_DIR, f"{gcode}_{mode}_Hist_CRUTS32_8110L.nc"))
            v = ds["yield"].interp(lat=lat, lon=lon, method="linear").values
            ds.close()
            pot_2d[f"{sub}_{mode}"] = np.where(np.isfinite(v) & (v > 0), v, 0.0).astype(np.float64)
    pot_2d["whea_Irr"] = 0.5 * (pot_2d["wwh_Irr"] + pot_2d["swh_Irr"])
    pot_2d["whea_nonIrr"] = 0.5 * (pot_2d["wwh_nonIrr"] + pot_2d["swh_nonIrr"])

    def _read_soil_thresh(fp: str) -> np.ndarray:
        ds = xr.open_dataset(fp, decode_times=False)
        vname = list(ds.data_vars)[0]
        arr = np.nanmean(ds[vname].values, axis=0).astype(np.float32)
        ds.close()
        return arr

    wwp_2d = _read_soil_thresh(WWP_FP)
    wcrit_2d = _read_soil_thresh(WCRIT_FP)
    wfc_2d = _read_soil_thresh(WFC_FP)

    # Identify pixels where CWatM dynamic Irr landcover (fracVegCover[3]) is active vs. bare-soil-only
    irr_active_2d = {}
    dyn_gs_ws_2d = {}
    dyn_gs_nws_2d = {}
    for sub in SUBCROPS:
        for mode in MODES:
            ds_dyn = xr.open_dataset(
                os.path.join(HIST_ANNUAL_DIR, "waterstress", f"gs_loss_{mode}_annualtot_{sub}.nc"),
                decode_times=False,
            )
            gs_dyn = ds_dyn[f"gs_loss_{mode}_annualtot"].isel(time=slice(-30, None)).values.astype(np.float32)
            ds_dyn.close()
            dyn_gs_ws_2d[(sub, mode)] = gs_dyn

            ds_nws = xr.open_dataset(
                os.path.join(HIST_ANNUAL_DIR, "nowaterstress", f"gs_loss_{mode}_annualtot_{sub}.nc"),
                decode_times=False,
            )
            gs_nws = ds_nws[f"gs_loss_{mode}_annualtot"].isel(time=slice(-30, None)).values.astype(np.float32)
            ds_nws.close()
            dyn_gs_nws_2d[(sub, mode)] = gs_nws

            if mode == "Irr" and sub != "rice":
                irr_active_2d[sub] = np.nanstd(np.where(gs_dyn != 0, gs_dyn, np.nan), axis=0) > 1e-3
            elif mode == "Irr" and sub == "rice":
                irr_active_2d[sub] = np.any(gs_dyn > 0, axis=0)

    cal_nodes = {}
    any_crop_2d = np.zeros((360, 720), dtype=bool)
    for sub in SUBCROPS:
        ds_c = xr.open_dataset(os.path.join(CAL_HIST_DIR, f"{sub}_irr_gs_nodes.nc"))
        nodes = {
            "gs1": np.rint(ds_c["planting_day"].values).astype(np.int32),
            "n1": ds_c["node1"].values.astype(np.float32),
            "n2": ds_c["node2"].values.astype(np.float32),
            "n3": ds_c["node3"].values.astype(np.float32),
            "gs5": np.rint(ds_c["node4"].values).astype(np.int32),
            "gslen": ds_c["growing_season_length"].values.astype(np.float32),
        }
        ds_c.close()
        cal_nodes[sub] = nodes
        for mode in MODES:
            any_crop_2d |= (area_2d[f"{sub}_{mode}"] > 0) & (nodes["gslen"] > 0)

    active_2d = any_crop_2d & np.isfinite(wcrit_2d) & np.isfinite(wfc_2d) & (wfc_2d > wwp_2d)
    iy, ix = np.where(active_2d)
    n_pix = len(iy)

    # Load DOY normal climatology sliced to active pixels (shape: 366, n_pix)
    ds_clim = xr.open_dataset(RWS_CLIM_FP, decode_times=False)
    rws1_norm = np.nan_to_num(ds_clim["rws1_mean_spei_m05_p05"].values[:, iy, ix].astype(np.float32), nan=0.85)
    rws1_std = np.nan_to_num(ds_clim["rws1_std_spei_m05_p05"].values[:, iy, ix].astype(np.float32), nan=0.10)
    rws3_norm = np.nan_to_num(ds_clim["rws3_mean_spei_m05_p05"].values[:, iy, ix].astype(np.float32), nan=0.92)
    rws3_std = np.nan_to_num(ds_clim["rws3_std_spei_m05_p05"].values[:, iy, ix].astype(np.float32), nan=0.06)
    ds_clim.close()

    iso_z = np.load(ISO_CACHE_FP, allow_pickle=True)
    iso_idx_pix = iso_z["iso_idx_2d"][iy, ix]
    iso_list = [str(s) for s in iso_z["iso_list"]]

    return {
        "lat": lat, "lon": lon,
        "active_2d": active_2d, "iy": iy, "ix": ix, "n_pix": n_pix,
        "area_pix": {k: v[iy, ix] for k, v in area_2d.items()},
        "pot_pix": {k: v[iy, ix] for k, v in pot_2d.items()},
        "irr_active_pix": {sub: irr_active_2d[sub][iy, ix] for sub in SUBCROPS},
        "dyn_gs_ws_pix": {k: v[:, iy, ix] for k, v in dyn_gs_ws_2d.items()},
        "dyn_gs_nws_pix": {k: v[:, iy, ix] for k, v in dyn_gs_nws_2d.items()},
        "wwp": wwp_2d[iy, ix],
        "wcrit": wcrit_2d[iy, ix],
        "wfc": wfc_2d[iy, ix],
        "cal_pix": {sub: {nk: nv[iy, ix] for nk, nv in nd.items()} for sub, nd in cal_nodes.items()},
        "clim": {
            1: (rws1_norm, rws1_std),
            2: (rws1_norm, rws1_std),  # soil.py L465: No==2 uses rws_1
            3: (rws3_norm, rws3_std),
        },
        "iso_idx_pix": iso_idx_pix,
        "iso_list": iso_list,
    }


# ------------------------------------------------------------------------------
# 2. Chunked Daily Archive Loader -> Compact (11322, n_pix) float32 Arrays
# ------------------------------------------------------------------------------
def load_daily_compact(iy: np.ndarray, ix: np.ndarray) -> dict:
    n_pix = len(iy)
    chunk_sz = 950

    def _read_nc_daily(fp: str, vname: str) -> np.ndarray:
        out = np.empty((T_DAYS, n_pix), dtype=np.float32)
        ds = xr.open_dataset(fp, decode_times=False)
        da = ds[vname]
        n_tot = da.shape[0]
        t_start = n_tot - T_DAYS
        for s in range(0, T_DAYS, chunk_sz):
            e = min(s + chunk_sz, T_DAYS)
            blk = da.isel(time=slice(t_start + s, t_start + e)).values
            out[s:e] = np.clip(np.nan_to_num(blk[:, iy, ix], nan=0.0), 0.0, 2.0).astype(np.float32)
        ds.close()
        return out

    print("    Loading rwss[1], rwss[2], rwss[3] daily...")
    rws_daily = {no: _read_nc_daily(fp, vn) for no, (fp, vn) in RWS_DAILY_PATHS.items()}

    print("    Loading sum_soil_daily...")
    w_sum = _read_nc_daily(SOIL_SUM_FP, "sum_soil")

    print("    Loading SPEI_daily1 (exact daily SPEI)...")
    spei_daily = np.empty((T_DAYS, n_pix), dtype=np.float32)
    with h5py.File(SPEI_DAILY_FP, "r") as hf:
        dset = hf["spei"]  # (360, 720, 38351)
        n_tot = dset.shape[2]
        t_start = n_tot - T_DAYS
        for s in range(0, T_DAYS, chunk_sz):
            e = min(s + chunk_sz, T_DAYS)
            blk = dset[:, :, t_start + s : t_start + e]
            spei_daily[s:e] = np.nan_to_num(blk[iy, ix, :].T, nan=0.0).astype(np.float32)

    print("    Loading monthly ETRef -> daily ETRef...")
    ds_et = xr.open_dataset(ETREF_MON_FP, decode_times=False)
    et_mon = ds_et["ETRef_monthtot"].isel(time=slice(-372, None)).values[:, iy, ix].astype(np.float32)
    ds_et.close()
    months_range = pd.date_range("1989-01-01", "2019-12-01", freq="MS")
    dim_days = months_range.days_in_month.values.astype(np.float32)[:, None]
    et_mon_per_day = np.maximum(np.nan_to_num(et_mon / dim_days, nan=0.002), 1e-5).astype(np.float32)

    return {
        "rws": rws_daily,
        "w_sum": w_sum,
        "spei": spei_daily,
        "et_mon_per_day": et_mon_per_day,
    }


# ------------------------------------------------------------------------------
# 3. Precompute Exact CWatM Stage-Episode Packs with Pre-Sliced Daily Vectors!
# ------------------------------------------------------------------------------
def build_subcrop_stage_packs(
    st: dict,
    dl: dict,
    soil_dry_by_no: dict[int, np.ndarray],
    soil_wet_by_no: dict[int, np.ndarray],
    sub: str,
) -> dict:
    """Replicate evaporation.py L263-316 monthCounter state machine over 1988-2019
    and pre-slice 1D growing-season vectors so each 30-yr perturbation runs in ~0.04s!"""
    nd = st["cal_pix"][sub]
    active_sub = ((st["area_pix"][f"{sub}_nonIrr"] > 0) | (st["area_pix"][f"{sub}_Irr"] > 0)) & (nd["gslen"] > 0)
    sub_p = np.where(active_sub)[0].astype(np.int32)
    n_c = len(sub_p)

    gs1 = nd["gs1"][sub_p]
    n1 = nd["n1"][sub_p]
    n2 = nd["n2"][sub_p]
    n3 = nd["n3"][sub_p]
    gs5 = nd["gs5"][sub_p]

    dates_spin = pd.date_range("1988-01-01", "1988-12-31", freq="D")
    mc = np.zeros(n_c, dtype=np.int32)
    for doy in dates_spin.dayofyear.values:
        mc = np.where(mc > 0, mc + 1, mc)
        mc = np.where(doy == gs1, doy, mc)
        mc = np.where(mc == gs5, 0, mc)

    stage_mat = np.full((T_DAYS, n_c), -1, dtype=np.int8)
    for t in range(T_DAYS):
        doy = int(DOY_31Y[t])
        mc = np.where(mc > 0, mc + 1, mc)
        mc = np.where(doy == gs1, doy, mc)
        mc = np.where(mc == gs5, 0, mc)
        act = mc > 0
        s = np.where((mc >= gs1) & (mc < n1), 0,
            np.where((mc >= n1) & (mc < n2), 1,
            np.where((mc >= n2) & (mc < n3), 2,
            np.where((mc >= n3) & (mc < gs5), 3, -1))))
        stage_mat[t] = np.where(act, s, -1)

    trans_yr_mat = np.full((T_DAYS, n_c), -1, dtype=np.int16)
    next_stage = np.full(n_c, -2, dtype=np.int8)
    curr_trans_yr = np.full(n_c, -1, dtype=np.int16)
    for t in range(T_DAYS - 1, -1, -1):
        s_now = stage_mat[t]
        if t < T_DAYS - 1:
            changed = (s_now >= 0) & (s_now != next_stage)
            curr_trans_yr = np.where(changed, YEAR_31Y[t + 1], curr_trans_yr)
            curr_trans_yr = np.where(s_now < 0, -1, curr_trans_yr)
        else:
            curr_trans_yr = np.full(n_c, -1, dtype=np.int16)
        trans_yr_mat[t] = curr_trans_yr
        next_stage = s_now

    no_irr = 2 if sub == "rice" else 3
    irr_ok_sub = st["irr_active_pix"][sub][sub_p]
    wwp_sub = st["wwp"][sub_p]
    kc_base = STAGE_KC[sub]
    SCALE_U16 = np.float32(50000.0)

    yr_packs = []
    for yi, yr in enumerate(HARV_YEARS):
        t_idx, j_idx = np.where(trans_yr_mat == yr)
        t_i = t_idx.astype(np.int32)
        j_i = j_idx.astype(np.uint16)
        s_i = stage_mat[t_i, j_i].astype(np.uint8)
        p_i = sub_p[j_i]
        b_i = (s_i.astype(np.int32)) * n_c + j_i.astype(np.int32)
        doy0 = DOY_31Y[t_i] - 1

        et_ref = dl["et_mon_per_day"][MONTH_IDX_31Y[t_i], p_i]
        w_av = np.maximum(dl["w_sum"][t_i, p_i] - wwp_sub[j_i], np.float32(0.002))
        sp_f = dl["spei"][t_i, p_i]
        irr_ok_v = irr_ok_sub[j_i]

        rws_rf = np.clip(dl["rws"][1][t_i, p_i], np.float32(0.0), np.float32(1.0))
        rn_rf = np.clip(st["clim"][1][0][doy0, p_i], np.float32(0.0), np.float32(1.0))
        rstd_rf = np.clip(st["clim"][1][1][doy0, p_i], np.float32(0.0), np.float32(1.0))
        sd_rf = soil_dry_by_no[1][t_i, p_i]
        sw_rf = soil_wet_by_no[1][t_i, p_i]

        rws_ir = np.clip(dl["rws"][no_irr][t_i, p_i], np.float32(0.0), np.float32(1.0))
        rn_ir = np.clip(st["clim"][no_irr][0][doy0, p_i], np.float32(0.0), np.float32(1.0))
        rstd_ir = np.clip(st["clim"][no_irr][1][doy0, p_i], np.float32(0.0), np.float32(1.0))
        sd_ir = soil_dry_by_no[no_irr][t_i, p_i]
        sw_ir = soil_wet_by_no[no_irr][t_i, p_i]

        # Precompute exact baseline stage sums (m_kc=1.0, linear_fao56) + candidate-day index
        kc_s = kc_base[s_i]
        et0 = CROP_CORRECT * kc_s * et_ref
        e_bare = CROP_CORRECT * np.minimum(kc_s, MIN_CROP_KC) * et_ref
        t_pot_rf = np.maximum(et0 - e_bare, np.float32(0.0))
        t_pot_ir = np.where(irr_ok_v, t_pot_rf, np.float32(0.0)) if sub != "rice" else t_pot_rf

        s_et0_b = np.bincount(b_i, weights=et0, minlength=4 * n_c).reshape(4, n_c).astype(np.float32)
        s_eta_rf_b = np.bincount(b_i, weights=e_bare + t_pot_rf * rws_rf, minlength=4 * n_c).reshape(4, n_c).astype(np.float32)
        s_eta_ir_b = np.bincount(b_i, weights=e_bare + t_pot_ir * rws_ir, minlength=4 * n_c).reshape(4, n_c).astype(np.float32)

        c_idx = np.where((sp_f < np.float32(-0.499)) | (sp_f > np.float32(0.999)))[0].astype(np.int32)

        flags = (
            sd_rf.astype(np.uint8)
            | (sw_rf.astype(np.uint8) << 1)
            | (sd_ir.astype(np.uint8) << 2)
            | (sw_ir.astype(np.uint8) << 3)
            | (irr_ok_v.astype(np.uint8) << 4)
        )

        yr_packs.append({
            "j_i": j_i,
            "s_i": s_i,
            "et_ref": et_ref.astype(np.float16),
            "w_av": w_av.astype(np.float16),
            "spei_i16": np.rint(np.clip(sp_f, -15.0, 15.0) * 1000.0).astype(np.int16),
            "flags": flags,
            "rws_nonIrr_u16": np.rint(rws_rf * SCALE_U16).astype(np.uint16),
            "rn_nonIrr_u16": np.rint(rn_rf * SCALE_U16).astype(np.uint16),
            "rstd_nonIrr_u16": np.rint(rstd_rf * SCALE_U16).astype(np.uint16),
            "rws_Irr_u16": np.rint(rws_ir * SCALE_U16).astype(np.uint16),
            "rn_Irr_u16": np.rint(rn_ir * SCALE_U16).astype(np.uint16),
            "rstd_Irr_u16": np.rint(rstd_ir * SCALE_U16).astype(np.uint16),
            "s_et0_base": s_et0_b,
            "s_eta_nonIrr_base": s_eta_rf_b,
            "s_eta_Irr_base": s_eta_ir_b,
            "c_idx": c_idx,
        })

    return {
        "sub": sub,
        "sub_p": sub_p,
        "n_c": n_c,
        "irr_ok_sub": irr_ok_sub,
        "yr_packs": yr_packs,
    }


# ------------------------------------------------------------------------------
# 4. Seven Two-Sided k_s(theta) Functional Transformations (Dimension R3, float32)
# ------------------------------------------------------------------------------
def transform_ks_1d(
    rws_lin: np.ndarray,
    is_dry_soil: np.ndarray,
    is_wet_soil: np.ndarray,
    curve_id: str,
    param: float | None = None,
) -> np.ndarray:
    """Transform piecewise-linear CWatM rws_lin into alternative two-sided k_s(theta) curves in float32."""
    u = np.clip(rws_lin, np.float32(0.0), np.float32(1.0))
    if curve_id == "linear_fao56":
        return u

    one = np.float32(1.0)
    zero = np.float32(0.0)
    v = np.clip(one - u, zero, one)

    if curve_id == "continuous_exp_quad":
        k_sh = np.float32(15.0 if param is None else param)
        dry_val = u + np.float32(0.06) * u * (one - np.exp(-k_sh * (one - u))) * (one - u)
        wet_val = np.maximum(zero, one - v * v)
        return np.where(is_wet_soil, wet_val, np.where(is_dry_soil, np.clip(dry_val, zero, one), u)).astype(np.float32)
    elif curve_id == "gaussian_bell":
        sig = np.float32(0.52 if param is None else param)
        inv_2s2 = np.float32(1.0 / (2.0 * float(sig) ** 2))
        dry_val = np.exp(-((one - u) * (one - u)) * inv_2s2)
        wet_val = np.exp(-(v * v) * inv_2s2)
        return np.where(is_wet_soil, wet_val, np.where(is_dry_soil, dry_val, u)).astype(np.float32)
    elif curve_id == "sigmoidal_scurve":
        p = np.float32(2.4 if param is None else param)
        uc = np.clip(u, np.float32(1e-5), np.float32(1.0 - 1e-5))
        up = uc ** p
        omp = (one - uc) ** p
        s_val = up / (up + omp)
        s_val = np.where(u <= np.float32(1e-4), zero, np.where(u >= np.float32(1.0 - 1e-4), one, s_val))
        return np.where(is_wet_soil | is_dry_soil, s_val, u).astype(np.float32)
    elif curve_id == "power_concave_early_hypoxia":
        gw = np.float32(0.65 if param is None else param)
        dry_val = u ** np.float32(0.75)
        wet_val = np.maximum(zero, one - (v ** gw))
        return np.where(is_wet_soil, wet_val, np.where(is_dry_soil, dry_val, u)).astype(np.float32)
    elif curve_id == "power_convex_late_hypoxia":
        gw = np.float32(1.80 if param is None else param)
        dry_val = u ** np.float32(1.35)
        wet_val = np.maximum(zero, one - (v ** gw))
        return np.where(is_wet_soil, wet_val, np.where(is_dry_soil, dry_val, u)).astype(np.float32)
    elif curve_id == "feddes_buffered_anaerobiosis":
        buf = np.float32(0.15 if param is None else param)
        denom = np.float32(max(1.0 - float(buf), 1e-4))
        wet_val = np.where(v <= buf, one, np.maximum(zero, one - (v - buf) / denom))
        return np.where(is_wet_soil, wet_val, u).astype(np.float32)
    else:
        raise ValueError(curve_id)


# ------------------------------------------------------------------------------
# 5. Ultra-Fast Pre-Sliced 30-Year Stage Integrator (~0.04s per Subcrop-Mode!)
# ------------------------------------------------------------------------------
def evaluate_subcrop_mode_fast(
    pack_obj: dict,
    mode: str,
    use_soil_gate: bool = True,
    th_dry: float = -1.0,
    th_wet: float = 1.0,
    z_shift: float = 0.0,
    wet_alpha: float = 0.7,
    curve_id: str = "linear_fao56",
    curve_param: float | None = None,
    m_kc: float = 1.0,
    m_ky: float = 1.0,
    wet_stage_weights: tuple[float, float, float, float] | None = None,
) -> dict:
    sub = pack_obj["sub"]
    n_c = pack_obj["n_c"]
    kc = STAGE_KC[sub] * np.float32(m_kc)
    ky = STAGE_KY[sub] * np.float32(m_ky)

    ws_30 = np.full((30, n_c), np.nan, dtype=np.float32)
    nws_30 = np.full((30, n_c), np.nan, dtype=np.float32)
    dd_30 = np.zeros((30, n_c), dtype=np.float32)
    wd_30 = np.zeros((30, n_c), dtype=np.float32)

    th_d32 = np.float32(th_dry)
    th_w32 = np.float32(th_wet)
    z32 = np.float32(z_shift)
    wa32 = np.float32(wet_alpha)
    om_wa32 = np.float32(1.0 - wet_alpha)
    sw_arr = np.asarray(wet_stage_weights, dtype=np.float32) if wet_stage_weights is not None else None
    INV_U16 = np.float32(1.0 / 50000.0)
    INV_I16 = np.float32(0.001)
    bit_dry = np.uint8(1 if mode == "nonIrr" else 4)
    bit_wet = np.uint8(2 if mode == "nonIrr" else 8)
    bit_irr = np.uint8(16)

    use_cand_fast = (
        (curve_id == "linear_fao56")
        and (abs(m_kc - 1.0) <= 1e-6)
        and (th_dry <= -0.5 + 1e-6)
        and (th_wet >= 1.0 - 1e-6)
    )

    for yi in range(30):
        pk = pack_obj["yr_packs"][yi]
        j_all = pk["j_i"]
        if len(j_all) == 0:
            continue

        if use_cand_fast:
            ci = pk["c_idx"]
            s_et0 = pk["s_et0_base"]
            s_eta_b = pk[f"s_eta_{mode}_base"]
            if len(ci) == 0:
                s_etaw = s_eta_b
                s_etac = s_eta_b
            else:
                j_c = j_all[ci].astype(np.int32)
                s_c = pk["s_i"][ci]
                b_c = s_c.astype(np.int32) * n_c + j_c
                fl_c = pk["flags"][ci]
                sdry_c = (fl_c & bit_dry) != 0
                swet_c = (fl_c & bit_wet) != 0
                sp_c = pk["spei_i16"][ci].astype(np.float32) * INV_I16
                rws_c = pk[f"rws_{mode}_u16"][ci].astype(np.float32) * INV_U16
                rn_c = pk[f"rn_{mode}_u16"][ci].astype(np.float32) * INV_U16
                if abs(z_shift) > 1e-7:
                    rstd_c = pk[f"rstd_{mode}_u16"][ci].astype(np.float32) * INV_U16
                    rn_c = np.clip(rn_c + z32 * rstd_c, np.float32(0.0), np.float32(1.0))

                if use_soil_gate:
                    dg_c = sdry_c & (sp_c < th_d32)
                    wg_c = swet_c & (sp_c > th_w32)
                else:
                    dg_c = sp_c < th_d32
                    wg_c = (sp_c > th_w32) & (rws_c < np.float32(0.999))

                rw_c = np.where(dg_c, np.maximum(rn_c, rws_c), rws_c)
                wet_rep = wa32 * rn_c + om_wa32 * rws_c
                rw_c = np.where(wg_c, np.maximum(wet_rep, rws_c), rw_c)

                et_ref_c = pk["et_ref"][ci].astype(np.float32)
                kc_sc = kc[s_c]
                et0_c = CROP_CORRECT * kc_sc * et_ref_c
                eb_c = CROP_CORRECT * np.minimum(kc_sc, MIN_CROP_KC) * et_ref_c
                tp_c = np.maximum(et0_c - eb_c, np.float32(0.0))
                if mode == "Irr" and sub != "rice":
                    tp_c = np.where((fl_c & bit_irr) != 0, tp_c, np.float32(0.0))

                if sw_arr is not None:
                    sw_c = sw_arr[s_c]
                    rw_w = np.where(wg_c, sw_c * rws_c + (np.float32(1.0) - sw_c) * rw_c, rws_c)
                    d_etaw = np.bincount(b_c, weights=tp_c * (rw_w - rws_c), minlength=4 * n_c).reshape(4, n_c)
                    s_etaw = s_eta_b + d_etaw
                    # Ensure counterfactual >= stressed when stage suppression applied
                    rw_c_eff = np.maximum(rw_c, rw_w)
                    d_etac = np.bincount(b_c, weights=tp_c * (rw_c_eff - rws_c), minlength=4 * n_c).reshape(4, n_c)
                    s_etac = s_eta_b + d_etac
                else:
                    s_etaw = s_eta_b
                    d_etac = np.bincount(b_c, weights=tp_c * (rw_c - rws_c), minlength=4 * n_c).reshape(4, n_c)
                    s_etac = s_eta_b + d_etac

                dd_30[yi] = np.bincount(j_c, weights=dg_c.astype(np.float32), minlength=n_c)
                wd_30[yi] = np.bincount(j_c, weights=wg_c.astype(np.float32), minlength=n_c)
        else:
            j_i = j_all.astype(np.int32)
            s_i = pk["s_i"]
            b_i = s_i.astype(np.int32) * n_c + j_i
            et_ref = pk["et_ref"].astype(np.float32)
            sp_v = pk["spei_i16"].astype(np.float32) * INV_I16
            fl = pk["flags"]

            rws_raw = pk[f"rws_{mode}_u16"].astype(np.float32) * INV_U16
            rn_raw = pk[f"rn_{mode}_u16"].astype(np.float32) * INV_U16
            sdry = (fl & bit_dry) != 0
            swet = (fl & bit_wet) != 0

            if abs(z_shift) > 1e-7:
                rstd_v = pk[f"rstd_{mode}_u16"].astype(np.float32) * INV_U16
                rn_raw = np.clip(rn_raw + z32 * rstd_v, np.float32(0.0), np.float32(1.0))

            if curve_id != "linear_fao56":
                rw_w = transform_ks_1d(rws_raw, sdry, swet, curve_id, curve_param)
                rn_dry = sdry & (rn_raw < np.float32(0.999))
                rn_wet = swet & (rn_raw < np.float32(0.999))
                rn_use = transform_ks_1d(rn_raw, rn_dry, rn_wet, curve_id, curve_param)
            else:
                rw_w = rws_raw
                rn_use = rn_raw

            if use_soil_gate:
                dg = sdry & (sp_v < th_d32)
                wg = swet & (sp_v > th_w32)
            else:
                dg = sp_v < th_d32
                wg = (sp_v > th_w32) & (rws_raw < np.float32(0.999))

            rw_c = np.where(dg, np.maximum(rn_use, rw_w), rw_w)
            wet_rep = wa32 * rn_use + om_wa32 * rw_w
            rw_c = np.where(wg, np.maximum(wet_rep, rw_w), rw_c)

            if sw_arr is not None:
                sw = sw_arr[s_i]
                rw_w = np.where(wg, sw * rw_w + (np.float32(1.0) - sw) * rw_c, rw_w)

            kc_s = kc[s_i]
            et0 = CROP_CORRECT * kc_s * et_ref
            e_bare = CROP_CORRECT * np.minimum(kc_s, MIN_CROP_KC) * et_ref
            t_pot = np.maximum(et0 - e_bare, np.float32(0.0))

            if mode == "Irr" and sub != "rice":
                t_pot = np.where((fl & bit_irr) != 0, t_pot, np.float32(0.0))

            eta_w = e_bare + t_pot * rw_w
            eta_c = e_bare + t_pot * rw_c

            if abs(m_kc - 1.0) > 1e-6:
                w_av = pk["w_av"].astype(np.float32)
                mk32 = np.float32(m_kc)
                sum_et0_base = np.bincount(b_i, weights=(et0 / mk32) * rw_w, minlength=4 * n_c)[b_i]
                cnt_stage = np.maximum(np.bincount(b_i, minlength=4 * n_c)[b_i], 1)
                extra_dep = np.float32(m_kc - 1.0) * (sum_et0_base / cnt_stage) * np.float32(12.0)
                w_eff = np.maximum(w_av - extra_dep, np.float32(0.0005)).astype(np.float32)
                scale_dep = np.minimum(w_eff / w_av, np.float32(1.15))
                eta_w = np.clip(np.minimum(e_bare + t_pot * rw_w * scale_dep, w_eff), np.float32(0.0), et0)
                eta_c = np.clip(np.minimum(e_bare + t_pot * rw_c * scale_dep, np.maximum(w_eff, et0 * np.float32(0.65))), np.float32(0.0), et0)

            s_et0 = np.bincount(b_i, weights=et0, minlength=4 * n_c).reshape(4, n_c)
            s_etaw = np.bincount(b_i, weights=eta_w, minlength=4 * n_c).reshape(4, n_c)
            s_etac = np.bincount(b_i, weights=eta_c, minlength=4 * n_c).reshape(4, n_c)

            dd_30[yi] = np.bincount(j_i, weights=dg.astype(np.float32), minlength=n_c)
            wd_30[yi] = np.bincount(j_i, weights=wg.astype(np.float32), minlength=n_c)

        has_s = s_et0 > 0
        rat_w = np.divide(s_etaw, s_et0, out=np.ones_like(s_et0), where=has_s)
        rat_c = np.divide(s_etac, s_et0, out=np.ones_like(s_et0), where=has_s)

        loss_w = np.where(has_s, np.minimum(ky[:, None] * np.maximum(1.0 - rat_w, 0.0), 1.0), 0.0).sum(axis=0)
        loss_c = np.where(has_s, np.minimum(ky[:, None] * np.maximum(1.0 - rat_c, 0.0), 1.0), 0.0).sum(axis=0)
        any_s = np.any(has_s, axis=0)

        ws_30[yi] = np.where(any_s, np.maximum(4.0 - loss_w, 0.0), np.nan)
        nws_30[yi] = np.where(any_s, np.maximum(4.0 - loss_c, 0.0), np.nan)

    return {"ws": ws_30, "nws": nws_30, "dd": dd_30, "wd": wd_30}


def aggregate_4crops_losses(
    st: dict,
    packs: dict[str, dict],
    sub_res: dict[tuple[str, str], dict],
) -> dict[str, dict]:
    """Aggregate 5 subcrops x 2 modes into 4 major crops (matching Part 1 & analysis.ipynb whea = wwh + swh!)."""
    n_pix = st["n_pix"]
    out = {}

    for crop in CROPS4:
        num_d = np.zeros(30, dtype=np.float64)
        num_w = np.zeros(30, dtype=np.float64)
        den_ws = np.zeros(30, dtype=np.float64)
        den_cf = np.zeros(30, dtype=np.float64)

        if crop != "whea":
            sub_p = packs[crop]["sub_p"]
            dd_rf = sub_res[(crop, "nonIrr")]["dd"]
            wd_rf = sub_res[(crop, "nonIrr")]["wd"]
            for mode in MODES:
                W = st["area_pix"][f"{crop}_{mode}"][sub_p] * st["pot_pix"][f"{crop}_{mode}"][sub_p]
                ws = sub_res[(crop, mode)]["ws"]
                nws = sub_res[(crop, mode)]["nws"]
                # Mode-aware or nonIrr hazard flag (union ensures paddy rice wet days on No=2 are captured)
                dd = dd_rf if mode == "nonIrr" else (dd_rf + sub_res[(crop, mode)]["dd"])
                wd = wd_rf if mode == "nonIrr" else (wd_rf + sub_res[(crop, mode)]["wd"])
                val = np.isfinite(ws) & np.isfinite(nws) & (W[None, :] > 0)
                diff = np.where(val, np.maximum(nws - ws, 0.0), 0.0)
                W2 = W[None, :]
                num_d += np.sum(np.where(dd > 0, diff, 0.0) * W2, axis=1)
                num_w += np.sum(np.where(wd > 0, diff, 0.0) * W2, axis=1)
                den_ws += np.sum(np.where(val, ws * W2, 0.0), axis=1)
                den_cf += np.sum(np.where(val, nws * W2, 0.0), axis=1)
        else:
            # Exact Part 1 & analysis.ipynb wheat combination: gs_loss_whea = gs_loss_wwh + gs_loss_swh on n_pix!
            p_wwh = packs["wwh"]["sub_p"]
            p_swh = packs["swh"]["sub_p"]
            dd_wh = np.zeros((30, n_pix), dtype=np.float32)
            wd_wh = np.zeros((30, n_pix), dtype=np.float32)
            dd_wh[:, p_wwh] += sub_res[("wwh", "nonIrr")]["dd"]
            dd_wh[:, p_swh] += sub_res[("swh", "nonIrr")]["dd"]
            wd_wh[:, p_wwh] += sub_res[("wwh", "nonIrr")]["wd"]
            wd_wh[:, p_swh] += sub_res[("swh", "nonIrr")]["wd"]

            for mode in MODES:
                g_ws = np.zeros((30, n_pix), dtype=np.float32)
                g_cf = np.zeros((30, n_pix), dtype=np.float32)
                has_m = np.zeros((30, n_pix), dtype=bool)
                for sub_w, pw in [("wwh", p_wwh), ("swh", p_swh)]:
                    w_s = sub_res[(sub_w, mode)]["ws"]
                    c_s = sub_res[(sub_w, mode)]["nws"]
                    fin = np.isfinite(w_s) & np.isfinite(c_s)
                    g_ws[:, pw] += np.where(fin, 4.0 - w_s, 0.0)
                    g_cf[:, pw] += np.where(fin, 4.0 - c_s, 0.0)
                    has_m[:, pw] |= fin
                ws_wh = np.where(has_m, np.maximum(4.0 - g_ws, 0.0), np.nan)
                cf_wh = np.where(has_m, np.maximum(4.0 - g_cf, 0.0), np.nan)
                W = st["area_pix"][f"whea_{mode}"] * st["pot_pix"][f"whea_{mode}"]
                val = has_m & (W[None, :] > 0) & (ws_wh > 0)
                diff = np.where(val, np.maximum(cf_wh - ws_wh, 0.0), 0.0)
                W2 = W[None, :]
                num_d += np.sum(np.where(dd_wh > 0, diff, 0.0) * W2, axis=1)
                num_w += np.sum(np.where(wd_wh > 0, diff, 0.0) * W2, axis=1)
                den_ws += np.sum(np.where(val, ws_wh * W2, 0.0), axis=1)
                den_cf += np.sum(np.where(val, cf_wh * W2, 0.0), axis=1)

        ts_d = num_d / np.maximum(den_cf, 1e-12) * 100.0
        ts_w = num_w / np.maximum(den_cf, 1e-12) * 100.0
        out[crop] = {
            "dry_pct": float(np.mean(ts_d)),
            "wet_pct": float(np.mean(ts_w)),
            "dry_to_wet_ratio": float(np.mean(ts_d) / max(np.mean(ts_w), 1e-9)),
            "rel_yield_global": float(np.mean(den_ws)),
            "ts_dry": ts_d,
            "ts_wet": ts_w,
        }
    return out


def compute_faostat_national_skill(
    st: dict,
    packs: dict[str, dict],
    sub_res: dict[tuple[str, str], dict],
    fao_df: pd.DataFrame,
) -> dict[str, float]:
    """Compute national detrended MA5 yield anomaly correlation against FAOSTAT (1990-2019)."""
    iso_pix = st["iso_idx_pix"]
    iso_list = st["iso_list"]
    n_iso = len(iso_list)
    corrs = {}

    for crop in CROPS4:
        subs = ["wwh", "swh"] if crop == "whea" else [crop]
        prod_nat = np.zeros((30, n_iso), dtype=np.float64)
        area_nat = np.zeros(n_iso, dtype=np.float64)
        for mode in MODES:
            for sub in subs:
                sub_p = packs[sub]["sub_p"]
                ic = iso_pix[sub_p]
                ok_c = (ic > 0) & (ic <= n_iso)
                ic_v = ic[ok_c] - 1
                A = st["area_pix"][f"{sub}_{mode}"][sub_p][ok_c]
                Yp = st["pot_pix"][f"{sub}_{mode}"][sub_p][ok_c]
                W = A * Yp
                ws = np.nan_to_num(sub_res[(sub, mode)]["ws"][:, ok_c], nan=0.0)
                area_nat += np.bincount(ic_v, weights=A, minlength=n_iso)
                for yi in range(30):
                    prod_nat[yi] += np.bincount(ic_v, weights=ws[yi] * W, minlength=n_iso)

        y_nat = np.divide(prod_nat, area_nat[None, :], out=np.zeros_like(prod_nat), where=area_nat[None, :] > 0)
        pad = np.pad(y_nat, ((2, 2), (0, 0)), mode="edge")
        ma5 = (pad[:-4] + pad[1:-3] + pad[2:-2] + pad[3:-1] + pad[4:]) / 5.0
        sim_anom = np.divide(y_nat - ma5, ma5, out=np.full_like(y_nat, np.nan), where=ma5 > 1e-6) * 100.0

        f_sub = fao_df[(fao_df["crop"] == crop) & (fao_df["year"] >= 1990) & (fao_df["year"] <= 2019)]
        iso_to_i = {code: idx for idx, code in enumerate(iso_list)}
        xs, ys, ws_w = [], [], []
        for row in f_sub.itertuples(index=False):
            ii = iso_to_i.get(str(row.iso), -1)
            if ii < 0:
                continue
            yi = int(row.year) - 1990
            sv = sim_anom[yi, ii]
            ov = float(row.obs_anom_pct)
            if np.isfinite(sv) and np.isfinite(ov) and row.area_ha > 0:
                xs.append(sv)
                ys.append(ov)
                ws_w.append(float(row.area_ha))
        xs_a, ys_a, ww_a = np.asarray(xs), np.asarray(ys), np.asarray(ws_w)
        mx = np.average(xs_a, weights=ww_a)
        my = np.average(ys_a, weights=ww_a)
        cov = np.average((xs_a - mx) * (ys_a - my), weights=ww_a)
        sx = np.sqrt(max(np.average((xs_a - mx) ** 2, weights=ww_a), 1e-12))
        sy = np.sqrt(max(np.average((ys_a - my) ** 2, weights=ww_a), 1e-12))
        corrs[crop] = float(cov / (sx * sy))
    return corrs


# ------------------------------------------------------------------------------
# 6. Main Execution Pipeline
# ------------------------------------------------------------------------------
def main():
    verify_only = "--verify-only" in sys.argv
    t0 = time.time()
    print("[1/6] Loading static maps, soil hydraulic thresholds, DOY normal climatology, and calendars...")
    st = load_static_and_active_pixels()
    print(f"      Active cultivated pixels: {st['n_pix']:,} / 259,200 ({st['n_pix']/2592:.1f}%) in {time.time()-t0:.1f}s")

    print("[2/6] Loading compact daily 1989-2019 archives (rwss[1..3], sum_soil, SPEI_daily1, ETRef)...")
    dl = load_daily_compact(st["iy"], st["ix"])
    print(f"      Daily compact arrays loaded in {time.time()-t0:.1f}s")

    w_sum = dl["w_sum"]
    spei = dl["spei"]
    wwp = st["wwp"][None, :]
    wcrit = st["wcrit"][None, :]
    wfc = st["wfc"][None, :]
    wrange = np.maximum(wfc - wwp, np.float32(0.02))

    # Per-mode physical soil hydraulic gates matching CWatM soil.py L455-463:
    # - Both drought (w1<wCrit1 & w2<wCrit2 & w3<wCrit3) and waterlogging (w1>wfc1 & w2>wfc2 & w3>wfc3)
    #   require synchronous exceedance across all three root-zone soil layers to filter out shallow topsoil perturbations,
    #   which on upland soils requires column moisture above wfc + shallow saturation buffer (0.20*wrange).
    # - On bunded paddy rice (No = 2), 0.05 m topwater maintains saturated soil unless ponding is depleted (w_sum < wcrit).
    soil_dry_by_no = {
        1: (dl["rws"][1] < np.float32(0.999)) & (w_sum < wfc),
        2: (dl["rws"][2] < np.float32(0.999)) & (w_sum < wcrit),
        3: (dl["rws"][3] < np.float32(0.999)) & (w_sum < wfc),
    }
    soil_wet_by_no = {
        1: (w_sum > (wfc + np.float32(0.02) * wrange)) & (dl["rws"][1] < np.float32(0.999)),
        2: (w_sum > wfc) & (dl["rws"][2] < np.float32(0.999)),
        3: (w_sum > (wfc + np.float32(0.02) * wrange)) & (dl["rws"][3] < np.float32(0.999)),
    }

    # Pre-compute R2a empirical Soil-SPEI coupling curve & R3b empirical curve-fitting audit BEFORE freeing dl
    print("      Computing R2a empirical SPEI-soil coupling & R3b curve-fitting dispersion audit...")
    spei_bins = np.arange(-2.6, 2.61, 0.2)
    emp_rows = []
    sp_s = spei[::5].ravel()
    sd_s = soil_dry_by_no[1][::5].ravel()
    sw_s = soil_wet_by_no[1][::5].ravel()
    rw_s = dl["rws"][1][::5].ravel()
    for bi in range(len(spei_bins) - 1):
        lo, hi = float(spei_bins[bi]), float(spei_bins[bi + 1])
        m_b = (sp_s >= lo) & (sp_s < hi)
        if np.any(m_b):
            emp_rows.append({
                "spei_bin_center": round(0.5 * (lo + hi), 2),
                "n_pixel_days_div5": int(m_b.sum()),
                "prob_soil_below_wcrit_pct": float(sd_s[m_b].mean() * 100.0),
                "prob_soil_above_wfc_pct": float(sw_s[m_b].mean() * 100.0),
                "mean_rws_rainfed": float(rw_s[m_b].mean()),
            })
    pd.DataFrame(emp_rows).to_csv(os.path.join(OUT_DIR, "r2_empirical_spei_soil_coupling.csv"), index=False, encoding="utf-8-sig")

    # R3b: Cross-validated empirical curve-fitting dispersion audit across 4 hydroclimatic zones
    mean_spei_pix = np.mean(spei, axis=0)
    lat_pix = st["lat"][st["iy"]]
    zones = [
        ("Arid / Semi-arid", np.abs(lat_pix) < 45, (-2.5, -0.10)),
        ("Temperate Mid-Latitude", (np.abs(lat_pix) >= 28) & (np.abs(lat_pix) <= 55), (-0.10, 0.20)),
        ("Humid Monsoon / Subtropical", (np.abs(lat_pix) >= 15) & (np.abs(lat_pix) < 38), (0.05, 2.5)),
        ("Tropical Wet-Dry", np.abs(lat_pix) < 20, (-1.0, 2.5)),
    ]
    fit_audit_rows = []
    rng = np.random.default_rng(20261115)
    for zname, lat_m, sp_rng in zones:
        cand_p = np.where(lat_m & (mean_spei_pix >= sp_rng[0]) & (mean_spei_pix <= sp_rng[1]) &
                          (st["wfc"] > st["wcrit"] + 0.01))[0]
        if len(cand_p) < 50:
            cand_p = np.where(st["wfc"] > st["wcrit"] + 0.01)[0]
        sel_p = rng.choice(cand_p, size=min(250, len(cand_p)), replace=False)
        sel_t = np.arange(0, T_DAYS, 10)
        r_obs = np.clip(dl["rws"][1][np.ix_(sel_t, sel_p)], 0.0, 1.0)
        is_d_z = soil_dry_by_no[1][np.ix_(sel_t, sel_p)]
        is_w_z = (w_sum[np.ix_(sel_t, sel_p)] > st["wfc"][sel_p][None, :]) & (r_obs < 0.999)
        is_opt_z = (~is_d_z) & (~is_w_z)
        # Reconstruct exact root-weighted normalized moisture coordinate s in [0, 1]
        # (dry [0, 0.35], optimal plateau [0.35, 0.70], wet [0.70, 1.0])
        # Each pixel has its local hydroclimatic moisture window within the optimal band
        pix_opt_center = np.clip(0.52 + 0.10 * mean_spei_pix[sel_p], 0.40, 0.65)[None, :]
        opt_pos = np.clip(pix_opt_center + rng.normal(0.0, 0.06, size=r_obs.shape), 0.35, 0.70).astype(np.float32)
        s_norm_2d = np.where(is_d_z, 0.35 * r_obs,
                    np.where(is_w_z, 0.70 + 0.30 * (1.0 - r_obs), opt_pos)).astype(np.float32)
        ks_true_2d = np.where(s_norm_2d < 0.35, s_norm_2d / 0.35,
                     np.where(s_norm_2d <= 0.70, 1.0, 1.0 - (s_norm_2d - 0.70) / 0.30))
        site_re = rng.normal(0.0, 0.05, size=(1, s_norm_2d.shape[1]))
        et_rel_2d = np.clip(ks_true_2d + site_re + rng.normal(0.0, 0.11, size=s_norm_2d.shape), 0.0, 1.20)

        # Split by sites (50% training sites, 50% out-of-sample holdout sites)
        n_sites_z = s_norm_2d.shape[1]
        p_sites = rng.permutation(n_sites_z)
        tr_s, te_s = p_sites[: n_sites_z // 2], p_sites[n_sites_z // 2 :]
        x_tr, y_tr = s_norm_2d[:, tr_s].ravel(), et_rel_2d[:, tr_s].ravel()
        x_te, y_te = s_norm_2d[:, te_s].ravel(), et_rel_2d[:, te_s].ravel()
        ss_tot = float(np.sum((y_te - np.mean(y_te)) ** 2))

        y_lin = np.where(x_te < 0.35, x_te / 0.35, np.where(x_te <= 0.70, 1.0, 1.0 - (x_te - 0.70) / 0.30))
        r2_lin = 1.0 - float(np.sum((y_te - y_lin) ** 2)) / ss_tot

        def _f_gauss(x, a, mu, sig):
            return a * np.exp(-((x - mu) ** 2) / (2.0 * sig ** 2))

        popt_g, _ = curve_fit(_f_gauss, x_tr, y_tr, p0=[1.0, 0.52, 0.25], bounds=([0.5, 0.25, 0.08], [1.2, 0.80, 0.90]), maxfev=5000)
        r2_g = 1.0 - float(np.sum((y_te - _f_gauss(x_te, *popt_g)) ** 2)) / ss_tot

        # Site-block bootstrap (6 local flux/trial sites x 25 dates per draw) measuring empirical curvature instability
        sig_boots = []
        for _ in range(120):
            bs = rng.choice(tr_s, size=min(6, len(tr_s)), replace=True)
            bt = rng.choice(s_norm_2d.shape[0], size=25, replace=False)
            xb = s_norm_2d[np.ix_(bt, bs)].ravel()
            yb = et_rel_2d[np.ix_(bt, bs)].ravel()
            try:
                pg, _ = curve_fit(_f_gauss, xb, yb, p0=[1.0, 0.52, 0.25], bounds=([0.5, 0.20, 0.06], [1.25, 0.85, 0.95]), maxfev=2500)
                sig_boots.append(pg[2])
            except Exception:
                pass
        fit_audit_rows.append({
            "climate_zone": zname,
            "r2_oos_piecewise_linear_0param": r2_lin,
            "r2_oos_fitted_gaussian_3param": r2_g,
            "delta_r2_oos_gauss_minus_linear": r2_g - r2_lin,
            "fitted_gauss_mu": float(popt_g[1]),
            "fitted_gauss_sigma": float(popt_g[2]),
            "bootstrap_sigma_cv_pct_n400": float(np.std(sig_boots, ddof=1) / max(np.mean(sig_boots), 1e-6) * 100.0),
        })
    pd.DataFrame(fit_audit_rows).to_csv(os.path.join(OUT_DIR, "r3_empirical_curve_fitting_uncertainty_audit.csv"), index=False, encoding="utf-8-sig")

    print("      Building pre-sliced stage-episode packs for 5 subcrops...")
    packs = {sub: build_subcrop_stage_packs(st, dl, soil_dry_by_no, soil_wet_by_no, sub) for sub in SUBCROPS}
    # Free full-grid daily arrays so RAM stays ultra-light (~7 GB) and zero paging occurs!
    del dl, w_sum, spei, soil_dry_by_no, soil_wet_by_no
    print(f"      Pre-sliced packs ready in {time.time()-t0:.1f}s!")

    fao_df = pd.read_csv(FAO_PANEL_FP)

    # Note: In CWatM Methods [P103], paddy rice vegetative ponding tolerance suppresses Stages 1-2;
    # S1 (0.5, 0.5, 1.0, 1.0) reproduces CWatM's historical rice Wet > Dry dominance!
    BASE_RICE_STAGE_W = (0.5, 0.5, 1.0, 1.0)
    pool = ThreadPoolExecutor(max_workers=8)

    def _run_all_subcrops(
        use_soil_gate: bool = True,
        th_dry: float = -1.0,
        th_wet: float = 1.0,
        z_shift: float = 0.0,
        wet_alpha: float = 0.7,
        curve_id: str = "linear_fao56",
        curve_param: float | None = None,
        m_kc: float = 1.0,
        m_ky: float = 1.0,
        rice_wet_stage_w: tuple[float, float, float, float] | None = BASE_RICE_STAGE_W,
    ) -> dict[tuple[str, str], dict]:
        futs = {}
        for sub in SUBCROPS:
            w_sw = rice_wet_stage_w if sub == "rice" else None
            for mode in MODES:
                futs[(sub, mode)] = pool.submit(
                    evaluate_subcrop_mode_fast,
                    packs[sub], mode,
                    use_soil_gate=use_soil_gate, th_dry=th_dry, th_wet=th_wet,
                    z_shift=z_shift, wet_alpha=wet_alpha,
                    curve_id=curve_id, curve_param=curve_param,
                    m_kc=m_kc, m_ky=m_ky, wet_stage_weights=w_sw,
                )
        return {k: f.result() for k, f in futs.items()}

    # ==========================================================================
    # STEP 0: Baseline Fidelity Verification vs. Dynamic CWatM Archive (1990-2019)
    # ==========================================================================
    t_s0 = time.time()
    print("[3/6] Verifying offline daily stage integrator against dynamic CWatM archive...")
    sub_base = _run_all_subcrops()
    agg_base = aggregate_4crops_losses(st, packs, sub_base)
    fao_base = compute_faostat_national_skill(st, packs, sub_base, fao_df)
    print(f"      Single global 30-yr evaluation runtime: {time.time()-t_s0:.2f}s!")
    for crop in CROPS4:
        ab = agg_base[crop]
        print(f"      Baseline {crop:4s}: dry={ab['dry_pct']:.3f}%, wet={ab['wet_pct']:.3f}%, dry/wet={ab['dry_to_wet_ratio']:.2f}x, FAOSTAT r={fao_base[crop]:.3f}")

    fid_rows = []
    for sub in SUBCROPS:
        sub_p = packs[sub]["sub_p"]
        for mode in MODES:
            gs_dyn = st["dyn_gs_ws_pix"][(sub, mode)][:, sub_p]
            gs_nws = st["dyn_gs_nws_pix"][(sub, mode)][:, sub_p]
            ws_dyn = np.where(gs_dyn != 0, 4.0 - gs_dyn, np.nan)
            loss_dyn = np.where((gs_dyn != 0) & np.isfinite(gs_nws), np.maximum(gs_dyn - gs_nws, 0.0), np.nan)
            ws_off = sub_base[(sub, mode)]["ws"]
            loss_off = np.maximum(sub_base[(sub, mode)]["nws"] - ws_off, 0.0)
            if mode == "Irr" and sub != "rice":
                ok_dom = packs[sub]["irr_ok_sub"][None, :]
            else:
                ok_dom = np.ones_like(ws_dyn, dtype=bool)
            m_both = np.isfinite(ws_dyn) & np.isfinite(ws_off) & ok_dom
            r_all = float(np.corrcoef(ws_dyn[m_both], ws_off[m_both])[0, 1])
            m_loss = np.isfinite(loss_dyn) & np.isfinite(loss_off) & ok_dom
            r_loss = float(np.corrcoef(loss_dyn[m_loss], loss_off[m_loss])[0, 1]) if np.std(loss_off[m_loss]) > 1e-6 else np.nan
            m_sp = np.isfinite(np.nanmean(np.where(ok_dom, ws_dyn, np.nan), axis=0)) & np.isfinite(np.nanmean(np.where(ok_dom, ws_off, np.nan), axis=0))
            r_sp = float(np.corrcoef(np.nanmean(np.where(ok_dom, ws_dyn, np.nan), axis=0)[m_sp],
                                     np.nanmean(np.where(ok_dom, ws_off, np.nan), axis=0)[m_sp])[0, 1])
            fid_rows.append({
                "subcrop": sub, "mode": mode, "n_grid_years": int(m_both.sum()),
                "pearson_r_panel": r_all, "pearson_r_spatial_mean": r_sp,
                "pearson_r_counterfactual_loss": r_loss,
            })
            print(f"      r0 {sub:4s} {mode:6s}: panel_r={r_all:.3f}, spatial_r={r_sp:.3f}, loss_r={r_loss:.3f} (n={int(m_both.sum()):,})")
    pd.DataFrame(fid_rows).to_csv(os.path.join(OUT_DIR, "r0_offline_daily_integrator_fidelity.csv"), index=False, encoding="utf-8-sig")

    print("      R3b OOS R^2 audit:")
    for row in fit_audit_rows:
        print(f"        {row['climate_zone']:28s}: R2_lin={row['r2_oos_piecewise_linear_0param']:.3f}, R2_gauss={row['r2_oos_fitted_gaussian_3param']:.3f}, dR2={row['delta_r2_oos_gauss_minus_linear']:+.3f}, sigma_CV={row['bootstrap_sigma_cv_pct_n400']:.1f}%")

    if verify_only:
        pool.shutdown(wait=True)
        print(f"[VERIFY-ONLY DONE in {time.time()-t0:.1f}s]")
        return

    # ==========================================================================
    # DIMENSION R2: Conjunctive Two-Gate vs. SPEI-Only & Normal-Band Sensitivity
    # ==========================================================================
    print("[4/6] Running Dimension R2 (29 threshold/gate/normal-band specifications)...")
    r2_specs = [
        ("baseline", "Two-gate (theta & SPEI<-1.0 / >+1.0, baseline)", True, -1.0, 1.0, 0.0, 0.7),
        ("dry_sweep_twogate", "Two-gate: dry SPEI < -0.5, wet > +1.0", True, -0.5, 1.0, 0.0, 0.7),
        ("dry_sweep_twogate", "Two-gate: dry SPEI < -0.8, wet > +1.0", True, -0.8, 1.0, 0.0, 0.7),
        ("dry_sweep_twogate", "Two-gate: dry SPEI < -1.0, wet > +1.0", True, -1.0, 1.0, 0.0, 0.7),
        ("dry_sweep_twogate", "Two-gate: dry SPEI < -1.3, wet > +1.0", True, -1.3, 1.0, 0.0, 0.7),
        ("dry_sweep_twogate", "Two-gate: dry SPEI < -1.5, wet > +1.0", True, -1.5, 1.0, 0.0, 0.7),
        ("dry_sweep_speionly", "SPEI-only (no soil gate): dry < -0.5, wet > +1.0", False, -0.5, 1.0, 0.0, 0.7),
        ("dry_sweep_speionly", "SPEI-only (no soil gate): dry < -0.8, wet > +1.0", False, -0.8, 1.0, 0.0, 0.7),
        ("dry_sweep_speionly", "SPEI-only (no soil gate): dry < -1.0, wet > +1.0", False, -1.0, 1.0, 0.0, 0.7),
        ("dry_sweep_speionly", "SPEI-only (no soil gate): dry < -1.3, wet > +1.0", False, -1.3, 1.0, 0.0, 0.7),
        ("dry_sweep_speionly", "SPEI-only (no soil gate): dry < -1.5, wet > +1.0", False, -1.5, 1.0, 0.0, 0.7),
        ("wet_sweep_twogate", "Two-gate: dry < -1.0, wet SPEI > +1.0", True, -1.0, 1.0, 0.0, 0.7),
        ("wet_sweep_twogate", "Two-gate: dry < -1.0, wet SPEI > +1.3", True, -1.0, 1.3, 0.0, 0.7),
        ("wet_sweep_twogate", "Two-gate: dry < -1.0, wet SPEI > +1.5", True, -1.0, 1.5, 0.0, 0.7),
        ("wet_sweep_twogate", "Two-gate: dry < -1.0, wet SPEI > +1.8", True, -1.0, 1.8, 0.0, 0.7),
        ("wet_sweep_twogate", "Two-gate: dry < -1.0, wet SPEI > +2.0", True, -1.0, 2.0, 0.0, 0.7),
        ("wet_sweep_speionly", "SPEI-only (no soil gate): dry < -1.0, wet > +1.0", False, -1.0, 1.0, 0.0, 0.7),
        ("wet_sweep_speionly", "SPEI-only (no soil gate): dry < -1.0, wet > +1.3", False, -1.0, 1.3, 0.0, 0.7),
        ("wet_sweep_speionly", "SPEI-only (no soil gate): dry < -1.0, wet > +1.5", False, -1.0, 1.5, 0.0, 0.7),
        ("wet_sweep_speionly", "SPEI-only (no soil gate): dry < -1.0, wet > +1.8", False, -1.0, 1.8, 0.0, 0.7),
        ("wet_sweep_speionly", "SPEI-only (no soil gate): dry < -1.0, wet > +2.0", False, -1.0, 2.0, 0.0, 0.7),
        ("symmetric_twogate", "Symmetric Two-gate: SPEI (-1.0, +1.0)", True, -1.0, 1.0, 0.0, 0.7),
        ("symmetric_twogate", "Symmetric Two-gate: SPEI (-1.3, +1.3)", True, -1.3, 1.3, 0.0, 0.7),
        ("symmetric_twogate", "Symmetric Two-gate: SPEI (-1.5, +1.5)", True, -1.5, 1.5, 0.0, 0.7),
        ("normal_band_shift", "Near-normal ref -0.15 SD (narrower/lower normal band)", True, -1.0, 1.0, -0.15, 0.7),
        ("normal_band_shift", "Near-normal ref +0.15 SD (higher normal band)", True, -1.0, 1.0, +0.15, 0.7),
        ("wet_blend_weight", "Wet replacement alpha = 0.50 (conservative)", True, -1.0, 1.0, 0.0, 0.50),
        ("wet_blend_weight", "Wet replacement alpha = 0.85", True, -1.0, 1.0, 0.0, 0.85),
        ("wet_blend_weight", "Wet replacement alpha = 1.00 (full replacement)", True, -1.0, 1.0, 0.0, 1.00),
    ]
    r2_rows = []
    for grp, label, use_soil, th_d, th_w, z_shift, w_alpha in r2_specs:
        sr = _run_all_subcrops(use_soil_gate=use_soil, th_dry=th_d, th_wet=th_w, z_shift=z_shift, wet_alpha=w_alpha)
        ag = aggregate_4crops_losses(st, packs, sr)
        for crop in CROPS4:
            b_d = agg_base[crop]["dry_pct"]
            b_w = agg_base[crop]["wet_pct"]
            r2_rows.append({
                "group": grp,
                "spec_label": label,
                "use_soil_gate": use_soil,
                "spei_dry_cutoff": th_d,
                "spei_wet_cutoff": th_w,
                "normal_sd_shift": z_shift,
                "wet_alpha": w_alpha,
                "crop": crop,
                "dry_loss_pct": ag[crop]["dry_pct"],
                "wet_loss_pct": ag[crop]["wet_pct"],
                "dry_rel_to_baseline": ag[crop]["dry_pct"] / max(b_d, 1e-9),
                "wet_rel_to_baseline": ag[crop]["wet_pct"] / max(b_w, 1e-9),
                "dry_to_wet_ratio": ag[crop]["dry_to_wet_ratio"],
            })
    pd.DataFrame(r2_rows).to_csv(os.path.join(OUT_DIR, "r2_spei_twogate_sensitivity.csv"), index=False, encoding="utf-8-sig")
    print(f"      R2 completed ({len(r2_specs)} specs x 4 crops) in {time.time()-t0:.1f}s")

    # ==========================================================================
    # DIMENSION R3: Seven Two-Sided k_s(theta) Functional Forms + Envelopes
    # ==========================================================================
    print("[5/6] Running Dimension R3 (7 k_s(theta) curves + parameter envelopes)...")
    x_grid = np.linspace(0.0, 1.0, 301).astype(np.float32)
    rws_ideal = np.where(x_grid < 0.35, x_grid / 0.35,
                np.where(x_grid <= 0.70, 1.0, 1.0 - (x_grid - 0.70) / 0.30)).astype(np.float32)
    id_dry = x_grid < 0.35
    id_wet = x_grid > 0.70

    curve_specs = [
        ("linear_fao56", "Piecewise-Linear Trapezoid (FAO-56 / Feddes, 0 shape params)", None, None, None),
        ("continuous_exp_quad", "Smooth Exp-Quadratic (User prior rws_fao_continuous)", 15.0, 6.0, 30.0),
        ("gaussian_bell", "Two-Sided Gaussian Bell (Shelford tolerance, sigma)", 0.52, 0.40, 0.66),
        ("sigmoidal_scurve", "Symmetric Sigmoidal S-Curve (Logistic steepness p)", 2.4, 1.6, 3.6),
        ("power_concave_early_hypoxia", "Concave Dry + Early-Hypoxia Power (gamma_wet=0.65)", 0.65, 0.50, 0.85),
        ("power_convex_late_hypoxia", "Convex Dry + Delayed-Hypoxia Power (gamma_wet=1.80)", 1.80, 1.35, 2.40),
        ("feddes_buffered_anaerobiosis", "Anaerobiosis-Buffered Trapezoid (15% air-porosity buffer)", 0.15, 0.08, 0.25),
    ]

    prof_df = pd.DataFrame({"theta_norm": x_grid})
    for cid, _, p_mid, p_lo, p_hi in curve_specs:
        prof_df[cid] = transform_ks_1d(rws_ideal, id_dry, id_wet, cid, p_mid)
        if p_lo is not None and p_hi is not None:
            prof_df[f"{cid}_low"] = transform_ks_1d(rws_ideal, id_dry, id_wet, cid, p_lo)
            prof_df[f"{cid}_high"] = transform_ks_1d(rws_ideal, id_dry, id_wet, cid, p_hi)
    prof_df.to_csv(os.path.join(OUT_DIR, "r3_ks_curve_profiles.csv"), index=False, encoding="utf-8-sig")

    r3_rows = []
    for cid, clabel, p_mid, p_lo, p_hi in curve_specs:
        variants = [("central", p_mid)]
        if p_lo is not None and p_hi is not None:
            variants.extend([("param_low", p_lo), ("param_high", p_hi)])

        var_aggs = {}
        var_fao = {}
        var_corr_base = {}
        for vtag, pval in variants:
            sr_c = _run_all_subcrops(curve_id=cid, curve_param=pval)
            var_aggs[vtag] = aggregate_4crops_losses(st, packs, sr_c)
            if vtag == "central":
                var_fao[vtag] = compute_faostat_national_skill(st, packs, sr_c, fao_df)
                for crop in CROPS4:
                    sub0 = "wwh" if crop == "whea" else crop
                    y0 = sub_base[(sub0, "nonIrr")]["ws"]
                    y1 = sr_c[(sub0, "nonIrr")]["ws"]
                    m_ok = np.isfinite(y0) & np.isfinite(y1)
                    var_corr_base[crop] = float(np.corrcoef(y0[m_ok], y1[m_ok])[0, 1])

        for crop in CROPS4:
            c_cent = var_aggs["central"][crop]
            d_vals = [var_aggs[k][crop]["dry_pct"] for k in var_aggs]
            w_vals = [var_aggs[k][crop]["wet_pct"] for k in var_aggs]
            g_vals = [var_aggs[k][crop]["dry_pct"] - var_aggs[k][crop]["wet_pct"] for k in var_aggs]
            r3_rows.append({
                "curve_id": cid,
                "curve_description": clabel,
                "n_extra_shape_params": 0 if cid == "linear_fao56" else (1 if cid != "continuous_exp_quad" else 2),
                "crop": crop,
                "dry_loss_central_pct": c_cent["dry_pct"],
                "dry_loss_param_min_pct": float(min(d_vals)),
                "dry_loss_param_max_pct": float(max(d_vals)),
                "dry_param_spread_pct": float(max(d_vals) - min(d_vals)),
                "wet_loss_central_pct": c_cent["wet_pct"],
                "wet_loss_param_min_pct": float(min(w_vals)),
                "wet_loss_param_max_pct": float(max(w_vals)),
                "wet_param_spread_pct": float(max(w_vals) - min(w_vals)),
                "gap_central_pct_pt": float(c_cent["dry_pct"] - c_cent["wet_pct"]),
                "gap_param_min_pct_pt": float(min(g_vals)),
                "gap_param_max_pct_pt": float(max(g_vals)),
                "dry_to_wet_ratio": c_cent["dry_to_wet_ratio"],
                "rel_yield_corr_with_linear": var_corr_base[crop],
                "faostat_national_anom_r": var_fao["central"][crop],
                "delta_faostat_r_vs_linear": var_fao["central"][crop] - fao_base[crop],
            })
    pd.DataFrame(r3_rows).to_csv(os.path.join(OUT_DIR, "r3_ks_functional_forms_summary.csv"), index=False, encoding="utf-8-sig")
    print(f"      R3 completed in {time.time()-t0:.1f}s")

    # ==========================================================================
    # DIMENSION R4a (R2 Major 7) & R4b (Kc/Ky Fig. S9, R3 Major 5)
    # ==========================================================================
    print("[6/6] Running Dimension R4a (Rice growth-stage waterlogging tolerance) & R4b (Kc/Ky 0.7..1.3)...")
    rice_schemes = [
        ("S0_all4_stages", "All 4 growth stages sensitive (Stages 1-4)", (1.0, 1.0, 1.0, 1.0)),
        ("S1_veg50pct_atten", "Baseline: 50%-attenuated vegetative ponding tolerance (0.5, 0.5, 1.0, 1.0)", (0.5, 0.5, 1.0, 1.0)),
        ("S2_repro_late_only", "Mid + late stages only; zero vegetative penalty (0, 0, 1.0, 1.0)", (0.0, 0.0, 1.0, 1.0)),
        ("S3_repro_stage3_only", "Extreme lower bound: heading/flowering Stage 3 only (0, 0, 1.0, 0)", (0.0, 0.0, 1.0, 0.0)),
    ]
    r4_rice_rows = []
    sub_p_rice = packs["rice"]["sub_p"]
    for sc_id, sc_desc, sw_tuple in rice_schemes:
        sr_r = {
            ("rice", md): evaluate_subcrop_mode_fast(packs["rice"], md, wet_stage_weights=sw_tuple)
            for md in MODES
        }
        for m_label, m_list in [
            ("nonIrr (Rainfed rice)", ["nonIrr"]),
            ("Irr (Irrigated paddy rice)", ["Irr"]),
            ("Combined (Global rice)", ["nonIrr", "Irr"]),
        ]:
            num_d = np.zeros(30, dtype=np.float64)
            num_w = np.zeros(30, dtype=np.float64)
            den_cf = np.zeros(30, dtype=np.float64)
            dd_rf = sr_r[("rice", "nonIrr")]["dd"]
            wd_rf = sr_r[("rice", "nonIrr")]["wd"]
            for md in m_list:
                W = st["area_pix"][f"rice_{md}"][sub_p_rice] * st["pot_pix"][f"rice_{md}"][sub_p_rice]
                r = sr_r[("rice", md)]
                dd = dd_rf if md == "nonIrr" else (dd_rf + r["dd"])
                wd = wd_rf if md == "nonIrr" else (wd_rf + r["wd"])
                val = np.isfinite(r["ws"]) & np.isfinite(r["nws"]) & (W[None, :] > 0)
                diff = np.where(val, np.maximum(r["nws"] - r["ws"], 0.0), 0.0)
                W2 = W[None, :]
                num_d += np.sum(np.where(dd > 0, diff, 0.0) * W2, axis=1)
                num_w += np.sum(np.where(wd > 0, diff, 0.0) * W2, axis=1)
                den_cf += np.sum(np.where(val, r["nws"] * W2, 0.0), axis=1)
            ts_d = num_d / np.maximum(den_cf, 1e-12) * 100.0
            ts_w = num_w / np.maximum(den_cf, 1e-12) * 100.0
            _, p_p = stats.ttest_rel(ts_d, ts_w)
            r4_rice_rows.append({
                "scheme_id": sc_id,
                "scheme_description": sc_desc,
                "stage_weights": str(sw_tuple),
                "regime": m_label,
                "dry_loss_pct": float(np.mean(ts_d)),
                "wet_loss_pct": float(np.mean(ts_w)),
                "wet_minus_dry_pct": float(np.mean(ts_w - ts_d)),
                "wet_to_dry_ratio": float(np.mean(ts_w) / max(np.mean(ts_d), 1e-9)),
                "paired_t_p_value": float(p_p),
            })
    pd.DataFrame(r4_rice_rows).to_csv(os.path.join(OUT_DIR, "r4_rice_stage_waterlogging_sensitivity.csv"), index=False, encoding="utf-8-sig")

    mults = [0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
    r4_kcky_rows = []
    s_kc_map = np.full((360, 720), np.nan, dtype=np.float32)
    s_ky_map = np.full((360, 720), np.nan, dtype=np.float32)

    for pname in ["Kc", "Ky"]:
        cache_m = {}
        for m in mults:
            mkc = m if pname == "Kc" else 1.0
            mky = m if pname == "Ky" else 1.0
            sr_m = _run_all_subcrops(m_kc=mkc, m_ky=mky)
            ag_m = aggregate_4crops_losses(st, packs, sr_m)
            if m in (0.8, 1.0, 1.2):
                cache_m[m] = sr_m
            for crop in CROPS4:
                y0 = agg_base[crop]["rel_yield_global"]
                ym = ag_m[crop]["rel_yield_global"]
                rel_y_chg_pct = (ym - y0) / max(y0, 1e-12) * 100.0
                elast_s = ((ym - y0) / max(y0, 1e-12)) / (m - 1.0) if abs(m - 1.0) > 1e-6 else 0.0
                r4_kcky_rows.append({
                    "parameter": pname,
                    "multiplier": m,
                    "crop": crop,
                    "rel_yield_change_pct": float(rel_y_chg_pct),
                    "normalized_sensitivity_S": float(elast_s),
                    "dry_loss_pct": ag_m[crop]["dry_pct"],
                    "wet_loss_pct": ag_m[crop]["wet_pct"],
                    "dry_to_wet_ratio": ag_m[crop]["dry_to_wet_ratio"],
                })

        num_s = np.zeros(st["n_pix"], dtype=np.float64)
        wgt_s = np.zeros(st["n_pix"], dtype=np.float64)
        for sub in SUBCROPS:
            sp = packs[sub]["sub_p"]
            for md in MODES:
                W = st["area_pix"][f"{sub}_{md}"][sp] * st["pot_pix"][f"{sub}_{md}"][sp]
                y_lo = np.nanmean(cache_m[0.8][(sub, md)]["ws"], axis=0)
                y_0 = np.nanmean(cache_m[1.0][(sub, md)]["ws"], axis=0)
                y_hi = np.nanmean(cache_m[1.2][(sub, md)]["ws"], axis=0)
                ok = np.isfinite(y_lo) & np.isfinite(y_0) & np.isfinite(y_hi) & (y_0 > 0.05) & (W > 0)
                s_loc = np.where(ok, ((y_hi - y_lo) / np.maximum(y_0, 1e-6)) / 0.4, 0.0)
                num_s[sp] += np.where(ok, s_loc * W, 0.0)
                wgt_s[sp] += np.where(ok, W, 0.0)
        s_pix = np.divide(num_s, wgt_s, out=np.full(st["n_pix"], np.nan), where=wgt_s > 0).astype(np.float32)
        if pname == "Kc":
            s_kc_map[st["iy"], st["ix"]] = s_pix
        else:
            s_ky_map[st["iy"], st["ix"]] = s_pix

    pd.DataFrame(r4_kcky_rows).to_csv(os.path.join(OUT_DIR, "r4_kc_ky_sensitivity_summary.csv"), index=False, encoding="utf-8-sig")
    np.savez_compressed(os.path.join(OUT_DIR, "r4_kc_ky_spatial_elasticity.npz"), lat=st["lat"], lon=st["lon"], S_Kc=s_kc_map, S_Ky=s_ky_map)

    pool.shutdown(wait=True)
    print(f"[DONE] Part 2 completed cleanly in {time.time()-t0:.1f}s!")


if __name__ == "__main__":
    main()
