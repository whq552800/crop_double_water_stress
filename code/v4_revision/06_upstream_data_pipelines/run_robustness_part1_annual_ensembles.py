# -*- coding: utf-8 -*-
"""Part 1 of Non-Rerun Robustness & Methodological Reliability Engine.

Covers Dimensions R0 (Exact Baseline Verification), R1 (Fixed 1990-2019 vs. Moving 2070-2100
Counterfactual Baselines across 10 CMIP6 GCMs x 2 SSPs), R5 (Conditional vs. Unconditional
Aggregation, Shapley Frequency x Severity Decomposition, and Compound Dry+Wet Partitioning),
R6 (Spatially & Temporally Robust Inference for Fig. 3), R7 (MIRCA2000 vs. SPAM2020V2r2,
Area x Yp vs. Pure Area, Fixed vs. Future Yp, Irr vs. nonIrr Event Masks, 1980-2019 vs.
1990-2019 Windows, Duration Cutoffs, & 2015-2025 Continuity), and R8 (CMIP6 Leave-One-GCM-Out,
No-Hot-Models Subset, Full 4-Way ANOVA, and National IDR Cutoff & Subperiod Sensitivity).

Strictly READ-ONLY on ./data/AAUDE, E:/, ./data/; writes ONLY to ./data/Software/Gemini/drought_waterlogging/results/robustness_norerun/.
"""
from __future__ import annotations

import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"

import sys
import time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
import xarray as xr
import rioxarray as rxr
from scipy import stats

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")

OUT_DIR = r".\results\robustness_norerun"
os.makedirs(OUT_DIR, exist_ok=True)

CELLAREA_FP = r"./data/external/cellarea.nc"
FRAC_DIR = r"./data/cwatm_inputs_figshare"
GAEZ_HIST_DIR = r"./data/gaez_v4_200b_noco2/hist"
GAEZ_FUT_DIR = r"./data/gaez_v4_200b_noco2/future"
CAL_HIST_DIR = r"./data/cwatm_inputs_figshare/after_clip"
CAL_SSP126_DIR = r"./data/cwatm_inputs_figshare/calendars/ssp126"
CAL_SSP585_DIR = r"./data/cwatm_inputs_figshare/calendars/ssp585"
SPAM_H_DIR = r"./data/sensitivity/spam2020_30arcmin/harvested_area"
SPAM_A_DIR = r"./data/sensitivity/spam2020_30arcmin/physical_area"

HIST_DIR = r"./data/cwatm_outputs/gswp3_long"
FUT_ROOT = r"./data/cwatm_outputs"

GCMS = [
    "CanESM5", "CNRM-CM6-1", "CNRM-ESM2-1", "EC-Earth3", "GFDL-ESM4",
    "IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL",
]
HOT_GCMS = {"CanESM5", "UKESM1-0-LL"}  # ECS > 5.0 K ("hot models", Hausfather et al. 2022 Nature)
CROPS = ["maiz", "soyb", "rice", "whea"]
MODES = ["nonIrr", "Irr"]
HAZARDS = ["dry", "wet"]


# ------------------------------------------------------------------------------
# 1. Exact Static Inputs (Grid, MIRCA aligned & unshifted, SPAM2020, GAEZ Yp, GS length)
# ------------------------------------------------------------------------------
def load_static_inputs() -> dict:
    ca_ds = xr.open_dataset(CELLAREA_FP)
    lat = ca_ds["lat"].values.astype(np.float64)
    lon = ca_ds["lon"].values.astype(np.float64)
    cellarea_m2 = ca_ds["cellarea30min"].values.astype(np.float64)
    ca_ds.close()

    # A1. Manuscript analysis.ipynb exact aligned MIRCA (in ha)
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
    mirca_aligned_ha = {}
    for k, fn in tif_aligned.items():
        da = rxr.open_rasterio(os.path.join(FRAC_DIR, fn)).squeeze("band", drop=True)
        da = da.rename({"y": "lat", "x": "lon"}).interp(lat=lat, lon=lon, method="linear")
        v = np.where(np.isfinite(da.values) & (da.values > 0), da.values, 0.0)
        mirca_aligned_ha[k] = (v * (cellarea_m2 / 1e4)).astype(np.float64)
    mirca_aligned_ha["whea_Irr"] = mirca_aligned_ha["wwh_Irr"] + mirca_aligned_ha["swh_Irr"]
    mirca_aligned_ha["whea_nonIrr"] = mirca_aligned_ha["wwh_nonIrr"] + mirca_aligned_ha["swh_nonIrr"]

    # A2. Unshifted exact-center MIRCA (in ha)
    tif_unshifted = {
        "wwh_Irr": "MIRCA-OS_Wheat_2000_irr_30arcmin_winter.tif",
        "swh_Irr": "MIRCA-OS_Wheat_2000_irr_30arcmin_spring.tif",
        "rice_Irr": "MIRCA-OS_Rice_2000_ir_30arcmin.tif",
        "maiz_Irr": "MIRCA-OS_Maize_2000_ir_30arcmin.tif",
        "soyb_Irr": "MIRCA-OS_Soybeans_2000_ir_30arcmin.tif",
        "wwh_nonIrr": "MIRCA-OS_Wheat_2000_rf_30arcmin_winter.tif",
        "swh_nonIrr": "MIRCA-OS_Wheat_2000_rf_30arcmin_spring.tif",
        "rice_nonIrr": "MIRCA-OS_Rice_2000_rf_30arcmin.tif",
        "maiz_nonIrr": "MIRCA-OS_Maize_2000_rf_30arcmin.tif",
        "soyb_nonIrr": "MIRCA-OS_Soybeans_2000_rf_30arcmin.tif",
    }
    mirca_unshifted_ha = {}
    for k, fn in tif_unshifted.items():
        da = rxr.open_rasterio(os.path.join(FRAC_DIR, fn)).squeeze("band", drop=True)
        v = np.where(np.isfinite(da.values) & (da.values > 0), da.values, 0.0)
        mirca_unshifted_ha[k] = (v * (cellarea_m2 / 1e4)).astype(np.float64)
    mirca_unshifted_ha["whea_Irr"] = mirca_unshifted_ha["wwh_Irr"] + mirca_unshifted_ha["swh_Irr"]
    mirca_unshifted_ha["whea_nonIrr"] = mirca_unshifted_ha["wwh_nonIrr"] + mirca_unshifted_ha["swh_nonIrr"]

    # A3. SPAM2020V2r2 Harvested Area (H) & Physical Area (A) conservatively aggregated 6x6 -> 0.5 deg
    spam_codes = {"maiz": "MAIZ", "soyb": "SOYB", "rice": "RICE", "whea": "WHEA"}
    spam_harv_ha, spam_phys_ha = {}, {}
    for crop, sc in spam_codes.items():
        for mode, sm in [("Irr", "I"), ("nonIrr", "R")]:
            for tag, sdir, targ in [("H", SPAM_H_DIR, spam_harv_ha), ("A", SPAM_A_DIR, spam_phys_ha)]:
                fp = os.path.join(sdir, f"spam2020_V2r2_global_{tag}_{sc}_{sm}.tif")
                arr = rxr.open_rasterio(fp).squeeze("band", drop=True).values.astype(np.float64)
                arr = np.where(np.isfinite(arr) & (arr > 0) & (arr < 1e20), arr, 0.0)
                targ[f"{crop}_{mode}"] = arr.reshape(360, 6, 720, 6).sum(axis=(1, 3))

    # B. GAEZ potential yield (Hist CRUTS32_8110L, Future ssp126 2080sH, Future ssp585 2080sH)
    def _load_pot(bdir: str, ssp: str, yr_tag: str) -> dict:
        out = {}
        for crop in CROPS:
            for mode in MODES:
                if crop == "rice":
                    da = xr.open_dataset(os.path.join(bdir, f"ricw_{mode}_{ssp}_{yr_tag}.nc"))["yield"]
                elif crop == "whea":
                    da1 = xr.open_dataset(os.path.join(bdir, f"wwhe_{mode}_{ssp}_{yr_tag}.nc"))["yield"]
                    da2 = xr.open_dataset(os.path.join(bdir, f"swhe_{mode}_{ssp}_{yr_tag}.nc"))["yield"]
                    da = da1.fillna(0) + da2.fillna(0)
                else:
                    da = xr.open_dataset(os.path.join(bdir, f"{crop}_{mode}_{ssp}_{yr_tag}.nc"))["yield"]
                v = da.interp(lat=lat, lon=lon, method="linear").values.astype(np.float64)
                out[f"{crop}_{mode}"] = np.where(np.isfinite(v) & (v > 0), v, 0.0)
        # Also keep separate wwh / swh potential yields for seasonal wheat check
        for sub, gz in [("wwh", "wwhe"), ("swh", "swhe")]:
            for mode in MODES:
                da = xr.open_dataset(os.path.join(bdir, f"{gz}_{mode}_{ssp}_{yr_tag}.nc"))["yield"]
                v = da.interp(lat=lat, lon=lon, method="linear").values.astype(np.float64)
                out[f"{sub}_{mode}"] = np.where(np.isfinite(v) & (v > 0), v, 0.0)
        return out

    pot_hist = _load_pot(GAEZ_HIST_DIR, "Hist", "CRUTS32_8110L")
    pot_126 = _load_pot(GAEZ_FUT_DIR, "ssp126", "2080sH")
    pot_585 = _load_pot(GAEZ_FUT_DIR, "ssp585", "2080sH")

    # C. Growing season length (days)
    def _load_gslen(cdir: str) -> dict:
        out = {}
        for crop in CROPS:
            if crop != "whea":
                v = xr.open_dataset(os.path.join(cdir, f"{crop}_irr_gs_nodes.nc"))["growing_season_length"].values
            else:
                v1 = xr.open_dataset(os.path.join(cdir, "wwh_irr_gs_nodes.nc"))["growing_season_length"].values
                v2 = xr.open_dataset(os.path.join(cdir, "swh_irr_gs_nodes.nc"))["growing_season_length"].values
                v = np.nan_to_num(v1, nan=0.0) + np.nan_to_num(v2, nan=0.0)
            out[crop] = np.where(np.isfinite(v) & (v > 0), v, np.nan).astype(np.float64)
        return out

    return {
        "lat": lat,
        "lon": lon,
        "cellarea_m2": cellarea_m2,
        "mirca_aligned": mirca_aligned_ha,
        "mirca_unshifted": mirca_unshifted_ha,
        "spam_harv": spam_harv_ha,
        "spam_phys": spam_phys_ha,
        "pot": {"hist": pot_hist, "ssp126": pot_126, "ssp585": pot_585},
        "gslen": {
            "hist": _load_gslen(CAL_HIST_DIR),
            "ssp126": _load_gslen(CAL_SSP126_DIR),
            "ssp585": _load_gslen(CAL_SSP585_DIR),
        },
    }


# ------------------------------------------------------------------------------
# 2. Fast CWatM Annual Scenario Loader (Exact analysis.ipynb + Irr stress masks)
# ------------------------------------------------------------------------------
def _read_nc_years(path: str, varname: str, target_years: np.ndarray) -> np.ndarray:
    ds = xr.open_dataset(path, decode_times=False)
    da = ds[varname]
    t = da["time"].values.astype(float)
    units = ds["time"].attrs.get("units", "").lower()
    if "since" in units and "year" in units:
        base_year = int(units.split("since")[1].strip()[:4])
        yrs = (base_year + np.rint(t)).astype(int) if np.nanmin(t) < 1000 else np.rint(t).astype(int)
    else:
        yrs = np.rint(t).astype(int)
    idx = [int(np.where(yrs == y)[0][0]) for y in target_years]
    arr = da.isel(time=idx).values.astype(np.float32)
    ds.close()
    arr = np.where(arr != 0, arr, np.nan)
    return arr


def load_cwatm_block(base_fp: str, target_years: np.ndarray) -> dict:
    """Load ws, nws, dry_nonIrr, wet_nonIrr, dry_Irr, wet_Irr for all 4 crops (plus wwh/swh)."""
    out = {}
    for crop in ["maiz", "soyb", "rice", "wwh", "swh"]:
        b = {}
        for mode in MODES:
            gs_ws = _read_nc_years(
                os.path.join(base_fp, "waterstress", f"gs_loss_{mode}_annualtot_{crop}.nc"),
                f"gs_loss_{mode}_annualtot", target_years,
            )
            gs_cf = _read_nc_years(
                os.path.join(base_fp, "nowaterstress", f"gs_loss_{mode}_annualtot_{crop}.nc"),
                f"gs_loss_{mode}_annualtot", target_years,
            )
            b[f"ws_{mode}"] = 4.0 - gs_ws
            b[f"nws_{mode}"] = 4.0 - gs_cf
            for haz in HAZARDS:
                sd = _read_nc_years(
                    os.path.join(base_fp, "nowaterstress", f"{haz}stress_{mode}_annualtot_{crop}.nc"),
                    f"{haz}stress_{mode}_annualtot", target_years,
                )
                b[f"{haz}_{mode}"] = np.nan_to_num(sd, nan=0.0)
        out[crop] = b

    # Build manuscript-exact 'whea' (sum of wwh + swh gs_loss, then 4 - sum)
    whe = {}
    for mode in MODES:
        gs_ws_w = (4.0 - out["wwh"][f"ws_{mode}"])
        gs_ws_s = (4.0 - out["swh"][f"ws_{mode}"])
        gs_ws_sum = np.nan_to_num(gs_ws_w, nan=0.0) + np.nan_to_num(gs_ws_s, nan=0.0)
        has_ws = np.isfinite(gs_ws_w) | np.isfinite(gs_ws_s)
        whe[f"ws_{mode}"] = np.where(has_ws & (gs_ws_sum != 0), 4.0 - gs_ws_sum, np.nan)

        gs_cf_w = (4.0 - out["wwh"][f"nws_{mode}"])
        gs_cf_s = (4.0 - out["swh"][f"nws_{mode}"])
        gs_cf_sum = np.nan_to_num(gs_cf_w, nan=0.0) + np.nan_to_num(gs_cf_s, nan=0.0)
        has_cf = np.isfinite(gs_cf_w) | np.isfinite(gs_cf_s)
        whe[f"nws_{mode}"] = np.where(has_cf & (gs_cf_sum != 0), 4.0 - gs_cf_sum, np.nan)

        for haz in HAZARDS:
            whe[f"{haz}_{mode}"] = out["wwh"][f"{haz}_{mode}"] + out["swh"][f"{haz}_{mode}"]
    out["whea"] = whe
    return out


# ------------------------------------------------------------------------------
# 3. Statistical Helpers (HAC, Moving-Block Bootstrap, Spatial-Year Two-Way Cluster)
# ------------------------------------------------------------------------------
def robust_paired_tests(d_series: np.ndarray, w_series: np.ndarray, seed: int = 20261115) -> dict:
    """Compute parametric, autocorrelation-adjusted, block-bootstrap, and non-parametric tests on x = d - w."""
    x = np.asarray(d_series - w_series, dtype=np.float64)
    n = len(x)
    mean_diff = float(np.mean(x))
    sd_diff = float(np.std(x, ddof=1))

    # 1. Standard paired t-test
    t_std, p_std = stats.ttest_rel(d_series, w_series)

    # 2. Lag-1 autocorrelation effective sample size (Santer et al. 2000 JGR)
    x_c = x - mean_diff
    r1 = float(np.sum(x_c[1:] * x_c[:-1]) / np.maximum(np.sum(x_c ** 2), 1e-12))
    r1_pos = max(0.0, min(0.95, r1))  # conservative positive AR(1) adjustment
    n_eff = n * (1.0 - r1_pos) / (1.0 + r1_pos)
    t_neff = mean_diff / (sd_diff / np.sqrt(max(n_eff, 2.0)))
    p_neff = float(2.0 * stats.t.sf(abs(t_neff), df=max(n_eff - 1.0, 1.0)))

    # 3. Newey-West HAC standard error (maxlag = 3 years)
    maxlag = 3
    gamma0 = np.sum(x_c ** 2) / n
    nw_var = gamma0
    for k in range(1, maxlag + 1):
        w_k = 1.0 - k / (maxlag + 1.0)
        gamma_k = np.sum(x_c[k:] * x_c[:-k]) / n
        nw_var += 2.0 * w_k * gamma_k
    se_hac = float(np.sqrt(max(nw_var, 1e-16) / n))
    t_hac = mean_diff / se_hac
    p_hac = float(2.0 * stats.t.sf(abs(t_hac), df=n - 1))

    # 4. Moving-block bootstrap (block length L = 5 yr, B = 10,000)
    rng = np.random.default_rng(seed)
    L = 5
    n_blocks = int(np.ceil(n / L))
    starts = rng.integers(0, n - L + 1, size=(10000, n_blocks))
    boot_idx = (starts[:, :, None] + np.arange(L)[None, None, :]).reshape(10000, -1)[:, :n]
    boot_means = x[boot_idx].mean(axis=1)
    ci_lo, ci_hi = np.percentile(boot_means, [2.5, 97.5])
    # Centered bootstrap null p-value
    boot_null = (x - mean_diff)[boot_idx].mean(axis=1)
    p_mbb = float((np.sum(np.abs(boot_null) >= abs(mean_diff)) + 1.0) / (len(boot_null) + 1.0))

    # 5. Non-parametric Wilcoxon signed-rank test
    try:
        _, p_wilcox = stats.wilcoxon(x, alternative="two-sided")
    except Exception:
        p_wilcox = np.nan

    return {
        "n_years": n,
        "mean_dry_pct": float(np.mean(d_series)),
        "mean_wet_pct": float(np.mean(w_series)),
        "mean_diff_dry_minus_wet_pct": mean_diff,
        "ratio_dry_to_wet": float(np.mean(d_series) / max(np.mean(w_series), 1e-9)),
        "lag1_autocorr_r1": r1,
        "n_eff_ar1": float(n_eff),
        "t_paired_std": float(t_std),
        "p_paired_std": float(p_std),
        "p_ar1_neff": p_neff,
        "se_newey_west_hac": se_hac,
        "p_newey_west_hac": p_hac,
        "mbb5_ci95_low": float(ci_lo),
        "mbb5_ci95_high": float(ci_hi),
        "p_mbb5_10000": p_mbb,
        "p_wilcoxon": float(p_wilcox),
    }


def two_way_cluster_pvalue(
    loss_d_3d: np.ndarray,
    loss_w_3d: np.ndarray,
    weight_2d: np.ndarray,
    block_deg_cells: int = 10,
) -> tuple[float, float, int]:
    """Cameron-Gelbach-Miller (2011) two-way cluster-robust SE (5x5 deg spatial blocks x year)."""
    # Aggregate net (dry - wet) production loss anomaly to 5-deg x 5-deg blocks x year
    ny, nlat, nlon = loss_d_3d.shape
    nb_lat = nlat // block_deg_cells
    nb_lon = nlon // block_deg_cells
    diff_3d = (loss_d_3d - loss_w_3d) * weight_2d[None, :, :]
    w_blocks = weight_2d.reshape(nb_lat, block_deg_cells, nb_lon, block_deg_cells).sum(axis=(1, 3))
    active_b = w_blocks > 0
    n_blocks = int(active_b.sum())
    b_ts = diff_3d.reshape(ny, nb_lat, block_deg_cells, nb_lon, block_deg_cells).sum(axis=(2, 4))[:, active_b]
    # Normalize so grand sum over blocks & years /ny equals global annual mean % loss diff
    # Let y_{t,b} be block b contribution in year t; global mean theta = (1/T) sum_t sum_b y_{t,b}
    # Equivalent panel regression of z_{t,b} = N_B * y_{t,b} on constant 1:
    z = n_blocks * b_ts  # shape (T, B), mean(z) == global mean diff
    T, B = z.shape
    N = T * B
    u = z - np.mean(z)
    # Cluster by spatial block b: sum over t first
    s_b = u.sum(axis=0)  # (B,)
    v_block = (B / max(B - 1, 1)) * np.sum(s_b ** 2) / (N ** 2)
    # Cluster by year t: sum over b first
    s_t = u.sum(axis=1)  # (T,)
    v_year = (T / max(T - 1, 1)) * np.sum(s_t ** 2) / (N ** 2)
    # Intersection (cell-year):
    v_inter = np.sum(u ** 2) / (N ** 2)
    v_twoway = max(v_block + v_year - v_inter, 1e-16)
    se_tw = float(np.sqrt(v_twoway))
    df_min = max(min(T, B) - 1, 1)
    t_tw = float(np.mean(z)) / se_tw
    p_tw = float(2.0 * stats.t.sf(abs(t_tw), df=df_min))
    return se_tw, p_tw, n_blocks


# ------------------------------------------------------------------------------
# 4. Core Loss Evaluator (Supports Flexible Mask, Weight, Denominator, Baseline, Duration)
# ------------------------------------------------------------------------------
def compute_crop_annual_series(
    block: dict,
    area_dict: dict,
    pot_dict: dict,
    crop: str,
    modes: tuple[str, ...] = ("nonIrr", "Irr"),
    denom_mode: str = "ws",          # 'ws' (manuscript Fig. 3) or 'nws' (counterfactual production)
    event_mask_source: str = "nonIrr",  # 'nonIrr' (manuscript analysis.ipynb) or 'mode' (Irr uses Irr)
    thr_days: float = 0.0,           # stress_days > thr_days (e.g. 0 -> >=1d; 6 -> >=7d)
    frac_th: float | None = None,    # if set, (stress_days / gslen) > frac_th
    gslen_2d: np.ndarray | None = None,
    use_pure_area: bool = False,     # True -> weight by Area only (Yp = 1)
    baseline_type: str = "fixed",    # 'fixed', 'moving_clim_p75', 'moving_ma5'
    apportion_compound: bool = False,  # True -> split compound dry+wet days proportionally
) -> dict:
    """Return annual dry, wet, compound, exposure F, and conditional severity I series."""
    ny = block[crop]["ws_nonIrr"].shape[0]
    num_dry = np.zeros(ny, dtype=np.float64)
    num_wet = np.zeros(ny, dtype=np.float64)
    num_pure_dry = np.zeros(ny, dtype=np.float64)
    num_pure_wet = np.zeros(ny, dtype=np.float64)
    num_comp_dry = np.zeros(ny, dtype=np.float64)
    num_comp_wet = np.zeros(ny, dtype=np.float64)
    den_tot = np.zeros(ny, dtype=np.float64)
    pot_tot = np.zeros(ny, dtype=np.float64)
    exp_dry = np.zeros(ny, dtype=np.float64)
    exp_wet = np.zeros(ny, dtype=np.float64)

    for mode in modes:
        A = area_dict[f"{crop}_{mode}"]
        Yp = np.ones_like(A) if use_pure_area else pot_dict[f"{crop}_{mode}"]
        W = np.where(np.isfinite(A) & np.isfinite(Yp), A * Yp, 0.0)

        ws_full = block[crop][f"ws_{mode}"]
        nws_full = block[crop][f"nws_{mode}"]
        msk_key = "nonIrr" if event_mask_source == "nonIrr" else mode
        dd_full = block[crop][f"dry_{msk_key}"]
        wd_full = block[crop][f"wet_{msk_key}"]

        act = (W > 0) & np.any(np.isfinite(ws_full), axis=0) & np.any(np.isfinite(nws_full), axis=0)
        if not np.any(act):
            continue
        W_1d = W[act]
        ws = ws_full[:, act]
        nws_fixed = nws_full[:, act]
        dd = dd_full[:, act]
        wd = wd_full[:, act]

        if frac_th is not None and gslen_2d is not None:
            gsl_1d = gslen_2d[act]
            m_dry = (dd / gsl_1d[None, :]) > frac_th
            m_wet = (wd / gsl_1d[None, :]) > frac_th
        else:
            m_dry = dd > thr_days
            m_wet = wd > thr_days

        # Counterfactual reference choice
        if baseline_type == "fixed":
            ref = nws_fixed
        elif baseline_type == "moving_clim_nonhaz":
            non_haz = (~m_dry) & (~m_wet) & np.isfinite(ws)
            cnt_nh = non_haz.sum(axis=0)
            sum_nh = np.where(non_haz, ws, 0.0).sum(axis=0)
            mean_nh = np.divide(sum_nh, cnt_nh, out=np.full_like(sum_nh, np.nan), where=cnt_nh >= 3)
            ws_fill = np.where(np.isfinite(ws), ws, np.nanmean(ws, axis=0, keepdims=True))
            p75_ws = np.percentile(ws_fill, 75, axis=0)
            ref_1d = np.where(np.isfinite(mean_nh), mean_nh, p75_ws)
            ref = np.broadcast_to(ref_1d[None, :], ws.shape)
            med_dd = np.median(dd, axis=0)[None, :]
            med_wd = np.median(wd, axis=0)[None, :]
            m_dry = (dd - med_dd) > thr_days
            m_wet = (wd - med_wd) > thr_days
        elif baseline_type == "moving_ma5":
            pad = np.pad(ws, ((2, 2), (0, 0)), mode="edge")
            ref = (pad[:-4] + pad[1:-3] + pad[2:-2] + pad[3:-1] + pad[4:]) / 5.0
            med_dd = np.median(dd, axis=0)[None, :]
            med_wd = np.median(wd, axis=0)[None, :]
            m_dry = (dd - med_dd) > thr_days
            m_wet = (wd - med_wd) > thr_days
        else:
            raise ValueError(baseline_type)

        valid = np.isfinite(ws) & np.isfinite(ref)
        diff = np.where(valid, ref - ws, 0.0)
        if baseline_type != "fixed":
            diff = np.maximum(diff, 0.0)

        both = m_dry & m_wet & valid
        only_d = m_dry & (~m_wet) & valid
        only_w = m_wet & (~m_dry) & valid

        if apportion_compound:
            tot_sd = np.maximum(dd + wd, 1e-9)
            sh_d = dd / tot_sd
            sh_w = wd / tot_sd
            d_loss_grid = np.where(only_d, diff, 0.0) + np.where(both, diff * sh_d, 0.0)
            w_loss_grid = np.where(only_w, diff, 0.0) + np.where(both, diff * sh_w, 0.0)
        else:
            d_loss_grid = np.where(m_dry & valid, diff, 0.0)
            w_loss_grid = np.where(m_wet & valid, diff, 0.0)

        W2 = W_1d[None, :]
        num_dry += np.sum(d_loss_grid * W2, axis=1)
        num_wet += np.sum(w_loss_grid * W2, axis=1)
        num_pure_dry += np.sum(np.where(only_d, diff, 0.0) * W2, axis=1)
        num_pure_wet += np.sum(np.where(only_w, diff, 0.0) * W2, axis=1)
        num_comp_dry += np.sum(np.where(both, diff, 0.0) * W2, axis=1)
        num_comp_wet += np.sum(np.where(both, diff, 0.0) * W2, axis=1)

        base_arr = ws if denom_mode == "ws" else ref
        den_tot += np.sum(np.where(valid, base_arr * W2, 0.0), axis=1)
        pot_tot += np.sum(np.where(valid, ref * W2, 0.0), axis=1)
        exp_dry += np.sum(np.where(m_dry & valid, base_arr * W2, 0.0), axis=1)
        exp_wet += np.sum(np.where(m_wet & valid, base_arr * W2, 0.0), axis=1)

    den_safe = np.maximum(den_tot, 1e-12)
    ts_dry = (num_dry / den_safe) * 100.0
    ts_wet = (num_wet / den_safe) * 100.0
    F_dry = exp_dry / den_safe  # fraction of production exposed [0..1]
    F_wet = exp_wet / den_safe
    I_dry = np.divide(num_dry, exp_dry, out=np.zeros_like(num_dry), where=exp_dry > 0) * 100.0
    I_wet = np.divide(num_wet, exp_wet, out=np.zeros_like(num_wet), where=exp_wet > 0) * 100.0

    return {
        "ts_dry": ts_dry,
        "ts_wet": ts_wet,
        "mt_dry": num_dry / 1e6,
        "mt_wet": num_wet / 1e6,
        "F_dry_mean": float(np.mean(F_dry)),
        "F_wet_mean": float(np.mean(F_wet)),
        "I_dry_pooled_pct": float(np.sum(num_dry) / max(np.sum(exp_dry), 1e-12) * 100.0),
        "I_wet_pooled_pct": float(np.sum(num_wet) / max(np.sum(exp_wet), 1e-12) * 100.0),
        "F_dry_pooled": float(np.sum(exp_dry) / max(np.sum(den_tot), 1e-12)),
        "F_wet_pooled": float(np.sum(exp_wet) / max(np.sum(den_tot), 1e-12)),
        "L_dry_pooled_pct": float(np.sum(num_dry) / max(np.sum(den_tot), 1e-12) * 100.0),
        "L_wet_pooled_pct": float(np.sum(num_wet) / max(np.sum(den_tot), 1e-12) * 100.0),
        "pure_dry_pct": float(np.mean(num_pure_dry / den_safe * 100.0)),
        "pure_wet_pct": float(np.mean(num_pure_wet / den_safe * 100.0)),
        "comp_overlap_pct": float(np.mean(num_comp_dry / den_safe * 100.0)),
    }


# ------------------------------------------------------------------------------
# Main Orchestrator
# ------------------------------------------------------------------------------
def main() -> None:
    t0 = time.time()
    print("[1/6] Loading static crop maps, SPAM2020V2r2, GAEZ Yp, and calendars...")
    st = load_static_inputs()
    print(f"      Done in {time.time() - t0:.1f}s")

    print("[2/6] Loading historical GSWP3-W5E5 (1980-2019) & 10 CMIP6 GCMs x 2 SSPs (2020-2029 & 2070-2100)...")
    yrs_hist40 = np.arange(1980, 2020, dtype=int)
    blk_hist40 = load_cwatm_block(HIST_DIR, yrs_hist40)

    # Slice helper for 1990-2019 (last 30 years) and 2010-2019 (last 10 years for continuity check)
    def _slice_block(blk: dict, sl: slice) -> dict:
        return {c: {k: v[sl] for k, v in d.items()} for c, d in blk.items()}

    blk_hist30 = _slice_block(blk_hist40, slice(10, 40))  # 1990-2019
    blk_hist_1019 = _slice_block(blk_hist40, slice(30, 40))  # 2010-2019

    yrs_fut_20s = np.arange(2020, 2030, dtype=int)
    yrs_fut_70s = np.arange(2070, 2101, dtype=int)
    yrs_fut_all = np.concatenate([yrs_fut_20s, yrs_fut_70s])  # 10 + 31 = 41 yrs
    tasks = [(gcm, scen) for gcm in GCMS for scen in ["ssp126", "ssp585"]]

    def _load_one_fut(args):
        gcm, scen = args
        base = os.path.join(FUT_ROOT, gcm, scen)
        probe_ds = xr.open_dataset(
            os.path.join(base, "waterstress", "gs_loss_nonIrr_annualtot_maiz.nc"), decode_times=False
        )
        t_min = int(1901 + np.rint(float(probe_ds["time"].values[0])))
        probe_ds.close()
        if t_min <= 2020:
            b = load_cwatm_block(base, yrs_fut_all)
            return gcm, scen, _slice_block(b, slice(0, 10)), _slice_block(b, slice(10, 41))
        else:
            b70 = load_cwatm_block(base, yrs_fut_70s)
            return gcm, scen, None, b70

    fut_2020s = {g: {} for g in GCMS}
    fut_2070s = {g: {} for g in GCMS}
    with ThreadPoolExecutor(max_workers=6) as ex:
        for gcm, scen, b20, b70 in ex.map(_load_one_fut, tasks):
            fut_2020s[gcm][scen] = b20
            fut_2070s[gcm][scen] = b70
    print(f"      Loaded all 20 GCM-SSP archives in {time.time() - t0:.1f}s")

    # ==========================================================================
    # DIMENSION R6 (+ R0 Exact Fig. 3 Verification): Spatially & Temporally Robust Inference
    # ==========================================================================
    print("[3/6] Running Dimension R6 (Fig. 3 baseline verification + HAC / Block-Bootstrap / Two-Way Cluster)...")
    r6_rows = []
    for win_name, blk_h in [("1990-2019 (n=30, Fig.3 baseline)", blk_hist30), ("1980-2019 (n=40, full hist)", blk_hist40)]:
        for crop in CROPS:
            for mode in MODES:
                res = compute_crop_annual_series(
                    blk_h, st["mirca_aligned"], st["pot"]["hist"], crop, modes=(mode,),
                    denom_mode="ws", event_mask_source="nonIrr", thr_days=0.0,
                )
                stat_d = robust_paired_tests(res["ts_dry"], res["ts_wet"])
                # Grid-level 3D arrays for two-way spatial-block x year cluster test
                A = st["mirca_aligned"][f"{crop}_{mode}"]
                Yp = st["pot"]["hist"][f"{crop}_{mode}"]
                W = A * Yp
                ws = blk_h[crop][f"ws_{mode}"]
                nws = blk_h[crop][f"nws_{mode}"]
                dd = blk_h[crop]["dry_nonIrr"]
                wd = blk_h[crop]["wet_nonIrr"]
                val = np.isfinite(ws) & np.isfinite(nws) & (W[None, :, :] > 0)
                den_mean = float(np.mean(np.nansum(np.where(val, ws * W[None, :, :], 0.0), axis=(1, 2))))
                ld_3d = np.where(val & (dd > 0), (nws - ws) / den_mean * 100.0, 0.0)
                lw_3d = np.where(val & (wd > 0), (nws - ws) / den_mean * 100.0, 0.0)
                se_tw, p_tw, nb = two_way_cluster_pvalue(ld_3d, lw_3d, W, block_deg_cells=10)
                row = {"window": win_name, "crop": crop, "mode": mode, **stat_d,
                       "se_twoway_5deg_year": se_tw, "p_twoway_5deg_year": p_tw, "n_spatial_5deg_blocks": nb}
                r6_rows.append(row)
                if "1990-2019" in win_name:
                    print(f"  [Fig.3 {crop} {mode:6s}] dry={stat_d['mean_dry_pct']:.2f}% vs wet={stat_d['mean_wet_pct']:.2f}% | "
                          f"p_std={stat_d['p_paired_std']:.2e}, p_HAC={stat_d['p_newey_west_hac']:.2e}, "
                          f"p_MBB5={stat_d['p_mbb5_10000']:.4f}, p_2way={p_tw:.2e}")
    pd.DataFrame(r6_rows).to_csv(os.path.join(OUT_DIR, "r6_fig3_robust_inference.csv"), index=False, encoding="utf-8-sig")

    # ==========================================================================
    # DIMENSION R1: Fixed Historical (1990-2019) vs. Moving Future (2070-2100) Baseline
    # ==========================================================================
    print("[4/6] Running Dimension R1 (Counterfactual baseline: Fixed 1990-2019 vs. Moving 2070-2100)...")
    r1_rows = []
    for base_scheme, base_label in [
        ("fixed", "Fixed 1990-2019 baseline (Manuscript standard: total climate + variability stress)"),
        ("moving_clim_nonhaz", "Moving 2070-2100 non-extreme climatology (M1: adapts to mean shift)"),
        ("moving_ma5", "Moving 5-yr centered mean baseline (M2: pure interannual weather shocks)"),
    ]:
        for crop in CROPS:
            # Historical reference under same scheme
            h_res = compute_crop_annual_series(
                blk_hist30, st["mirca_aligned"], st["pot"]["hist"], crop, modes=("nonIrr", "Irr"),
                denom_mode="nws", baseline_type=base_scheme, thr_days=0.0,
            )
            r1_rows.append({
                "baseline_scheme": base_scheme,
                "baseline_description": base_label,
                "scenario": "hist_1990_2019",
                "crop": crop,
                "dry_mean_pct": float(np.mean(h_res["ts_dry"])),
                "dry_gcm_sd_pct": 0.0,
                "dry_p10_pct": float(np.percentile(h_res["ts_dry"], 10)),
                "dry_p90_pct": float(np.percentile(h_res["ts_dry"], 90)),
                "wet_mean_pct": float(np.mean(h_res["ts_wet"])),
                "wet_gcm_sd_pct": 0.0,
                "wet_p10_pct": float(np.percentile(h_res["ts_wet"], 10)),
                "wet_p90_pct": float(np.percentile(h_res["ts_wet"], 90)),
                "dry_mt_yr": float(np.mean(h_res["mt_dry"])),
                "wet_mt_yr": float(np.mean(h_res["mt_wet"])),
                "dry_to_wet_ratio": float(np.mean(h_res["ts_dry"]) / max(np.mean(h_res["ts_wet"]), 1e-9)),
            })
            for scen in ["ssp126", "ssp585"]:
                d_gcms, w_gcms, dm_gcms, wm_gcms = [], [], [], []
                for gcm in GCMS:
                    f_res = compute_crop_annual_series(
                        fut_2070s[gcm][scen], st["mirca_aligned"], st["pot"][scen], crop, modes=("nonIrr", "Irr"),
                        denom_mode="nws", baseline_type=base_scheme, thr_days=0.0,
                    )
                    d_gcms.append(float(np.mean(f_res["ts_dry"])))
                    w_gcms.append(float(np.mean(f_res["ts_wet"])))
                    dm_gcms.append(float(np.mean(f_res["mt_dry"])))
                    wm_gcms.append(float(np.mean(f_res["mt_wet"])))
                r1_rows.append({
                    "baseline_scheme": base_scheme,
                    "baseline_description": base_label,
                    "scenario": f"{scen}_2070_2100",
                    "crop": crop,
                    "dry_mean_pct": float(np.mean(d_gcms)),
                    "dry_gcm_sd_pct": float(np.std(d_gcms, ddof=1)),
                    "dry_p10_pct": float(np.percentile(d_gcms, 10)),
                    "dry_p90_pct": float(np.percentile(d_gcms, 90)),
                    "wet_mean_pct": float(np.mean(w_gcms)),
                    "wet_gcm_sd_pct": float(np.std(w_gcms, ddof=1)),
                    "wet_p10_pct": float(np.percentile(w_gcms, 10)),
                    "wet_p90_pct": float(np.percentile(w_gcms, 90)),
                    "dry_mt_yr": float(np.mean(dm_gcms)),
                    "wet_mt_yr": float(np.mean(wm_gcms)),
                    "dry_to_wet_ratio": float(np.mean(d_gcms) / max(np.mean(w_gcms), 1e-9)),
                })
    df_r1 = pd.DataFrame(r1_rows)
    df_r1.to_csv(os.path.join(OUT_DIR, "r1_counterfactual_baseline_fixed_vs_moving.csv"), index=False, encoding="utf-8-sig")
    print("  [R1 Summary SSP5-8.5 2070-2100]:")
    for c in CROPS:
        sub = df_r1[(df_r1["crop"] == c) & (df_r1["scenario"] == "ssp585_2070_2100")]
        fx = sub[sub["baseline_scheme"] == "fixed"].iloc[0]
        m1 = sub[sub["baseline_scheme"] == "moving_clim_nonhaz"].iloc[0]
        print(f"    {c}: Fixed dry={fx['dry_mean_pct']:.2f}% / wet={fx['wet_mean_pct']:.2f}%  vs  "
              f"Moving-M1 dry={m1['dry_mean_pct']:.2f}% / wet={m1['wet_mean_pct']:.2f}%")

    # ==========================================================================
    # DIMENSION R5: Aggregation Conventions, Shapley Decomposition & Compound Dry+Wet
    # ==========================================================================
    print("[5/6] Running Dimension R5 (Conditional vs. Unconditional, Shapley dF vs. dI, & Compound Dry+Wet)...")
    r5_rows = []
    for crop in CROPS:
        h = compute_crop_annual_series(
            blk_hist30, st["mirca_aligned"], st["pot"]["hist"], crop, modes=("nonIrr", "Irr"),
            denom_mode="nws", thr_days=0.0, apportion_compound=False,
        )
        h_app = compute_crop_annual_series(
            blk_hist30, st["mirca_aligned"], st["pot"]["hist"], crop, modes=("nonIrr", "Irr"),
            denom_mode="nws", thr_days=0.0, apportion_compound=True,
        )
        h_wsden = compute_crop_annual_series(
            blk_hist30, st["mirca_aligned"], st["pot"]["hist"], crop, modes=("nonIrr", "Irr"),
            denom_mode="ws", thr_days=0.0, apportion_compound=False,
        )
        for scen in ["ssp126", "ssp585"]:
            f_list = [
                compute_crop_annual_series(
                    fut_2070s[gcm][scen], st["mirca_aligned"], st["pot"][scen], crop, modes=("nonIrr", "Irr"),
                    denom_mode="nws", thr_days=0.0, apportion_compound=False,
                )
                for gcm in GCMS
            ]
            f_app_list = [
                compute_crop_annual_series(
                    fut_2070s[gcm][scen], st["mirca_aligned"], st["pot"][scen], crop, modes=("nonIrr", "Irr"),
                    denom_mode="nws", thr_days=0.0, apportion_compound=True,
                )
                for gcm in GCMS
            ]
            for haz in HAZARDS:
                F0 = h[f"F_{haz}_pooled"]
                I0 = h[f"I_{haz}_pooled_pct"]
                L0 = h[f"L_{haz}_pooled_pct"]
                # Exact GCM-by-GCM Shapley & 3-term decomposition of L = F * I
                dL_all, dF_contrib, dI_contrib, inter_contrib, shap_F, shap_I = [], [], [], [], [], []
                F1_all, I1_all, L1_all = [], [], []
                for fg in f_list:
                    F1 = fg[f"F_{haz}_pooled"]
                    I1 = fg[f"I_{haz}_pooled_pct"]
                    L1 = fg[f"L_{haz}_pooled_pct"]
                    dF = F1 - F0
                    dI = I1 - I0
                    F1_all.append(F1)
                    I1_all.append(I1)
                    L1_all.append(L1)
                    dL_all.append(L1 - L0)
                    dF_contrib.append(I0 * dF)
                    dI_contrib.append(F0 * dI)
                    inter_contrib.append(dF * dI)
                    shap_F.append(0.5 * (I0 + I1) * dF)
                    shap_I.append(0.5 * (F0 + F1) * dI)

                tot_dL = float(np.mean(dL_all))
                sF = float(np.mean(shap_F))
                sI = float(np.mean(shap_I))
                r5_rows.append({
                    "crop": crop,
                    "scenario": scen,
                    "hazard": haz,
                    "hist_uncond_loss_nws_denom_pct": L0,
                    "hist_uncond_loss_ws_denom_pct": h_wsden[f"L_{haz}_pooled_pct"],
                    "hist_apportioned_loss_pct": float(np.mean(h_app[f"ts_{haz}"])),
                    "hist_pure_single_haz_pct": h[f"pure_{haz}_pct"],
                    "hist_compound_overlap_pct": h["comp_overlap_pct"],
                    "hist_exposure_freq_F_pct": F0 * 100.0,
                    "hist_cond_severity_I_pct": I0,
                    "fut_uncond_loss_pct": float(np.mean(L1_all)),
                    "fut_apportioned_loss_pct": float(np.mean([np.mean(x[f"ts_{haz}"]) for x in f_app_list])),
                    "fut_compound_overlap_pct": float(np.mean([x["comp_overlap_pct"] for x in f_list])),
                    "fut_exposure_freq_F_pct": float(np.mean(F1_all) * 100.0),
                    "fut_cond_severity_I_pct": float(np.mean(I1_all)),
                    "delta_loss_total_pct_pt": tot_dL,
                    "term1_freq_extent_I0_dF": float(np.mean(dF_contrib)),
                    "term2_cond_severity_F0_dI": float(np.mean(dI_contrib)),
                    "term3_interaction_dF_dI": float(np.mean(inter_contrib)),
                    "shapley_freq_extent_pct_pt": sF,
                    "shapley_cond_severity_pct_pt": sI,
                    "shapley_freq_share_of_abs_pct": abs(sF) / max(abs(sF) + abs(sI), 1e-12) * 100.0,
                    "shapley_sev_share_of_abs_pct": abs(sI) / max(abs(sF) + abs(sI), 1e-12) * 100.0,
                })
    df_r5 = pd.DataFrame(r5_rows)
    df_r5.to_csv(os.path.join(OUT_DIR, "r5_aggregation_decomposition_compound.csv"), index=False, encoding="utf-8-sig")

    # Compute clean SPAM2020V2r2 occurrence, affected-area share, and conditional duration table (R2 Minor 10)
    def _calc_occ_metrics(blk: dict, area_map: dict, crop: str, mode: str, haz: str) -> tuple[float, float, float, float]:
        A = area_map[f"{crop}_{mode}"]
        ws = blk[crop][f"ws_{mode}"]
        nws = blk[crop][f"nws_{mode}"]
        mc = (A > 0) & np.any(np.isfinite(ws), axis=0) & np.any(np.isfinite(nws), axis=0)
        A1 = A[mc]
        ws1 = ws[:, mc]
        nws1 = nws[:, mc]
        sd = blk[crop][f"{haz}_nonIrr"][:, mc]
        hb = (sd > 0) & np.isfinite(ws1)
        gf = float(hb.sum(axis=1).mean() / mc.sum() * 100.0)
        af = float(((hb * A1[None, :]).sum(axis=1) / A1.sum()).mean() * 100.0)
        num_d = (np.where(hb, sd, 0.0) * A1[None, :]).sum(axis=1)
        den_d = (hb * A1[None, :]).sum(axis=1)
        valid_y = den_d > 0
        dur = float(np.mean(num_d[valid_y] / den_d[valid_y])) if np.any(valid_y) else 0.0
        hb_eff = hb & ((nws1 - ws1) > 1e-6)
        af_eff = float(((hb_eff * A1[None, :]).sum(axis=1) / A1.sum()).mean() * 100.0)
        return gf, af, dur, af_eff

    occ_summary = []
    for h_label, h_blk in [("Historical (1990-2019)", blk_hist30), ("Historical (1980-2019)", blk_hist40)]:
        for crop in CROPS:
            for mode in ["Irr", "nonIrr"]:
                for haz in HAZARDS:
                    gf, af, dur, af_eff = _calc_occ_metrics(h_blk, st["spam_harv"], crop, mode, haz)
                    occ_summary.append({
                        "scenario": h_label,
                        "crop": crop, "mode": mode, "hazard": haz,
                        "gridyear_freq_pct": gf, "gridyear_freq_sd_pct": 0.0,
                        "affected_area_pct": af, "affected_area_sd_pct": 0.0,
                        "cond_duration_days": dur, "cond_duration_sd_days": 0.0,
                        "loss_inducing_area_pct": af_eff, "loss_inducing_area_sd_pct": 0.0,
                    })
    for scen in ["ssp126", "ssp585"]:
        scen_lbl = "SSP1-2.6 (2070-2100)" if scen == "ssp126" else "SSP5-8.5 (2070-2100)"
        for crop in CROPS:
            for mode in ["Irr", "nonIrr"]:
                for haz in HAZARDS:
                    vals = np.array([
                        _calc_occ_metrics(fut_2070s[gcm][scen], st["spam_harv"], crop, mode, haz)
                        for gcm in GCMS
                    ], dtype=np.float64)
                    occ_summary.append({
                        "scenario": scen_lbl,
                        "crop": crop, "mode": mode, "hazard": haz,
                        "gridyear_freq_pct": float(vals[:, 0].mean()),
                        "gridyear_freq_sd_pct": float(vals[:, 0].std(ddof=1)),
                        "affected_area_pct": float(vals[:, 1].mean()),
                        "affected_area_sd_pct": float(vals[:, 1].std(ddof=1)),
                        "cond_duration_days": float(vals[:, 2].mean()),
                        "cond_duration_sd_days": float(vals[:, 2].std(ddof=1)),
                        "loss_inducing_area_pct": float(vals[:, 3].mean()),
                        "loss_inducing_area_sd_pct": float(vals[:, 3].std(ddof=1)),
                    })
    pd.DataFrame(occ_summary).to_csv(os.path.join(OUT_DIR, "r5_occurrence_summary_table.csv"), index=False, encoding="utf-8-sig")

    # ==========================================================================
    # DIMENSION R7: Masks (MIRCA vs SPAM2020), Weights, Fixed Yp, Duration Cutoffs, Windows & Continuity
    # ==========================================================================
    print("[6/6] Running Dimensions R7 & R8 (Masks/Weights/Duration/Continuity + CMIP6 ANOVA/HotModels/IDR)...")
    r7_rows = []
    pert_specs = [
        ("Baseline (MIRCA2000, A*Yp, >=1d, 1990-2019)", blk_hist30, st["mirca_aligned"], False, False, "nonIrr", 0.0, None),
        ("Window 1980-2019 (40 yr)", blk_hist40, st["mirca_aligned"], False, False, "nonIrr", 0.0, None),
        ("Unshifted MIRCA-OS grid", blk_hist30, st["mirca_unshifted"], False, False, "nonIrr", 0.0, None),
        ("SPAM2020V2r2 Harvested area", blk_hist30, st["spam_harv"], False, False, "nonIrr", 0.0, None),
        ("SPAM2020V2r2 Physical area", blk_hist30, st["spam_phys"], False, False, "nonIrr", 0.0, None),
        ("Pure harvested-area weight (Yp = 1)", blk_hist30, st["mirca_aligned"], True, False, "nonIrr", 0.0, None),
        ("Fixed Historical Yp in Future", blk_hist30, st["mirca_aligned"], False, True, "nonIrr", 0.0, None),
        ("Mode-specific Irr stress mask", blk_hist30, st["mirca_aligned"], False, False, "mode", 0.0, None),
        ("Duration >= 3 days", blk_hist30, st["mirca_aligned"], False, False, "nonIrr", 2.0, None),
        ("Duration >= 5 days", blk_hist30, st["mirca_aligned"], False, False, "nonIrr", 4.0, None),
        ("Duration >= 7 days (1 week)", blk_hist30, st["mirca_aligned"], False, False, "nonIrr", 6.0, None),
        ("Duration >= 10 days", blk_hist30, st["mirca_aligned"], False, False, "nonIrr", 9.0, None),
        ("Duration >= 14 days (2 weeks)", blk_hist30, st["mirca_aligned"], False, False, "nonIrr", 13.0, None),
        ("GS fraction > 1%", blk_hist30, st["mirca_aligned"], False, False, "nonIrr", 0.0, 0.01),
        ("GS fraction > 3%", blk_hist30, st["mirca_aligned"], False, False, "nonIrr", 0.0, 0.03),
        ("GS fraction > 5%", blk_hist30, st["mirca_aligned"], False, False, "nonIrr", 0.0, 0.05),
    ]
    for label, b_h, a_dict, pure_a, fix_yp, msk_src, thr_d, f_th in pert_specs:
        for crop in CROPS:
            h_r = compute_crop_annual_series(
                b_h, a_dict, st["pot"]["hist"], crop, modes=("nonIrr", "Irr"),
                denom_mode="ws", event_mask_source=msk_src, thr_days=thr_d,
                frac_th=f_th, gslen_2d=st["gslen"]["hist"][crop], use_pure_area=pure_a,
            )
            row = {
                "specification": label,
                "crop": crop,
                "hist_dry_pct": float(np.mean(h_r["ts_dry"])),
                "hist_wet_pct": float(np.mean(h_r["ts_wet"])),
                "hist_dry_to_wet_ratio": float(np.mean(h_r["ts_dry"]) / max(np.mean(h_r["ts_wet"]), 1e-9)),
            }
            for scen in ["ssp126", "ssp585"]:
                pot_use = st["pot"]["hist"] if fix_yp else st["pot"][scen]
                d_f, w_f = [], []
                for gcm in GCMS:
                    f_r = compute_crop_annual_series(
                        fut_2070s[gcm][scen], a_dict, pot_use, crop, modes=("nonIrr", "Irr"),
                        denom_mode="ws", event_mask_source=msk_src, thr_days=thr_d,
                        frac_th=f_th, gslen_2d=st["gslen"][scen][crop], use_pure_area=pure_a,
                    )
                    d_f.append(float(np.mean(f_r["ts_dry"])))
                    w_f.append(float(np.mean(f_r["ts_wet"])))
                row[f"{scen}_dry_pct"] = float(np.mean(d_f))
                row[f"{scen}_wet_pct"] = float(np.mean(w_f))
                row[f"{scen}_delta_dry_pct"] = float(np.mean(d_f) - row["hist_dry_pct"])
                row[f"{scen}_delta_wet_pct"] = float(np.mean(w_f) - row["hist_wet_pct"])
            r7_rows.append(row)

    # Also check GSWP3-W5E5 (2010-2019) -> CMIP6 (2020-2029) transition continuity (R3 Major 3)
    cont_rows = []
    for crop in CROPS:
        h10 = compute_crop_annual_series(blk_hist_1019, st["mirca_aligned"], st["pot"]["hist"], crop, denom_mode="ws")
        for scen in ["ssp126", "ssp585"]:
            d20, w20 = [], []
            for gcm in GCMS:
                if fut_2020s[gcm][scen] is None:
                    continue
                # Evaluate 2020-2029 with fixed historical Yp to isolate pure climate forcing transition!
                f20 = compute_crop_annual_series(fut_2020s[gcm][scen], st["mirca_aligned"], st["pot"]["hist"], crop, denom_mode="ws")
                d20.append(float(np.mean(f20["ts_dry"])))
                w20.append(float(np.mean(f20["ts_wet"])))
            cont_rows.append({
                "crop": crop, "scenario": scen,
                "n_gcms_with_2020_2029": int(len(d20)),
                "hist_2010_2019_dry_pct": float(np.mean(h10["ts_dry"])),
                "cmip6_2020_2029_dry_mean_pct": float(np.mean(d20)),
                "cmip6_2020_2029_dry_p10_p90": f"[{np.percentile(d20,10):.2f}, {np.percentile(d20,90):.2f}]",
                "hist_2010_2019_wet_pct": float(np.mean(h10["ts_wet"])),
                "cmip6_2020_2029_wet_mean_pct": float(np.mean(w20)),
                "cmip6_2020_2029_wet_p10_p90": f"[{np.percentile(w20,10):.2f}, {np.percentile(w20,90):.2f}]",
            })
    pd.DataFrame(r7_rows).to_csv(os.path.join(OUT_DIR, "r7_masks_weights_windows_duration.csv"), index=False, encoding="utf-8-sig")
    pd.DataFrame(cont_rows).to_csv(os.path.join(OUT_DIR, "r7_transition_2010_2029_continuity.csv"), index=False, encoding="utf-8-sig")

    # ==========================================================================
    # DIMENSION R8: CMIP6 Ensemble (Leave-One-Out, No-Hot-Models, 4-Way ANOVA) + National IDR Sweep
    # ==========================================================================
    r8_anova_rows = []
    r8_loo_rows = []
    for crop in CROPS:
        # Build Y_gst array of shape (S=2, G=10, T=31) for both dry and wet
        arr_dry = np.zeros((2, len(GCMS), 31), dtype=np.float64)
        arr_wet = np.zeros((2, len(GCMS), 31), dtype=np.float64)
        for si, scen in enumerate(["ssp126", "ssp585"]):
            for gi, gcm in enumerate(GCMS):
                fr = compute_crop_annual_series(
                    fut_2070s[gcm][scen], st["mirca_aligned"], st["pot"][scen], crop, denom_mode="nws"
                )
                arr_dry[si, gi, :] = fr["ts_dry"]
                arr_wet[si, gi, :] = fr["ts_wet"]

        for haz, arr in [("dry", arr_dry), ("wet", arr_wet)]:
            # 1) 2-way ANOVA on 31-year climatological means M_{s,g} (matches Fig. 7 + adds SSP x GCM interaction!)
            M_sg = arr.mean(axis=2)  # (2, 10)
            grand_m = M_sg.mean()
            m_s = M_sg.mean(axis=1, keepdims=True)
            m_g = M_sg.mean(axis=0, keepdims=True)
            S_len, G_len = M_sg.shape
            ss_ssp_clim = G_len * np.sum((m_s - grand_m) ** 2)
            ss_gcm_clim = S_len * np.sum((m_g - grand_m) ** 2)
            ss_int_clim = np.sum((M_sg - m_s - m_g + grand_m) ** 2)
            ss_tot_clim = max(ss_ssp_clim + ss_gcm_clim + ss_int_clim, 1e-16)
            ss_main2_clim = max(ss_ssp_clim + ss_gcm_clim, 1e-16)

            # 2) Full 4-component ANOVA on annual projections Y_{s,g,t} (S=2, G=10, T=31)
            T_len = arr.shape[2]
            ss_ssp_ann = T_len * ss_ssp_clim
            ss_gcm_ann = T_len * ss_gcm_clim
            ss_int_ann = T_len * ss_int_clim
            ss_internal = np.sum((arr - M_sg[:, :, None]) ** 2)
            ss_tot_ann = max(ss_ssp_ann + ss_gcm_ann + ss_int_ann + ss_internal, 1e-16)

            # 3) Ensemble subsets: Full 10-GCM vs. 8-GCM No-Hot-Models vs. Leave-One-Out range
            idx_nohot = [i for i, g in enumerate(GCMS) if g not in HOT_GCMS]
            loo_126 = [float(np.delete(M_sg[0], i).mean()) for i in range(len(GCMS))]
            loo_585 = [float(np.delete(M_sg[1], i).mean()) for i in range(len(GCMS))]

            r8_anova_rows.append({
                "crop": crop,
                "hazard": haz,
                "main2_ssp_share_pct": float(ss_ssp_clim / ss_main2_clim * 100.0),
                "main2_gcm_share_pct": float(ss_gcm_clim / ss_main2_clim * 100.0),
                "clim31_ssp_share_pct": float(ss_ssp_clim / ss_tot_clim * 100.0),
                "clim31_gcm_share_pct": float(ss_gcm_clim / ss_tot_clim * 100.0),
                "clim31_sspxgcm_share_pct": float(ss_int_clim / ss_tot_clim * 100.0),
                "annual_ssp_share_pct": float(ss_ssp_ann / ss_tot_ann * 100.0),
                "annual_gcm_share_pct": float(ss_gcm_ann / ss_tot_ann * 100.0),
                "annual_sspxgcm_share_pct": float(ss_int_ann / ss_tot_ann * 100.0),
                "annual_internal_var_share_pct": float(ss_internal / ss_tot_ann * 100.0),
                "ssp126_10gcm_mean_pct": float(M_sg[0].mean()),
                "ssp126_8gcm_nohot_mean_pct": float(M_sg[0, idx_nohot].mean()),
                "ssp126_loo_min_pct": float(min(loo_126)),
                "ssp126_loo_max_pct": float(max(loo_126)),
                "ssp585_10gcm_mean_pct": float(M_sg[1].mean()),
                "ssp585_8gcm_nohot_mean_pct": float(M_sg[1, idx_nohot].mean()),
                "ssp585_loo_min_pct": float(min(loo_585)),
                "ssp585_loo_max_pct": float(max(loo_585)),
            })
            for gi, gcm in enumerate(GCMS):
                r8_loo_rows.append({
                    "crop": crop, "hazard": haz, "dropped_gcm": gcm,
                    "is_hot_model": gcm in HOT_GCMS,
                    "gcm_own_ssp126_pct": float(M_sg[0, gi]),
                    "gcm_own_ssp585_pct": float(M_sg[1, gi]),
                    "loo_mean_ssp126_pct": loo_126[gi],
                    "loo_mean_ssp585_pct": loo_585[gi],
                })

    pd.DataFrame(r8_anova_rows).to_csv(os.path.join(OUT_DIR, "r8_gcm_anova_hotmodels.csv"), index=False, encoding="utf-8-sig")
    pd.DataFrame(r8_loo_rows).to_csv(os.path.join(OUT_DIR, "r8_gcm_leave_one_out.csv"), index=False, encoding="utf-8-sig")

    # National IDR & Delta-Loss Cutoff Sweep + Subperiod Persistence (2014-2018 vs 2019-2023)
    df_risk = pd.read_csv(r"./outputs/outputs_country_prod_loss/country_food_risk_table.csv")
    idr_rows = []
    for scen in ["ssp126", "ssp585"]:
        for haz in ["dry", "wet"]:
            h_sub = df_risk[(df_risk["scenario"] == "hist") & (df_risk["hazard"] == haz)].set_index("ISO_A3")
            f_sub = df_risk[(df_risk["scenario"] == scen) & (df_risk["hazard"] == haz)].set_index("ISO_A3")
            common = h_sub.index.intersection(f_sub.index)
            dl = f_sub.loc[common, "loss_pct"] - h_sub.loc[common, "loss_pct"]
            idr = f_sub.loc[common, "import_dependence"]
            prod = f_sub.loc[common, "production"]
            inc = f_sub.loc[common, "income_group"].fillna("Unknown")
            valid_c = idr.notna() & dl.notna() & (prod > 0)
            dl, idr, prod, inc = dl[valid_c], idr[valid_c], prod[valid_c], inc[valid_c]
            tot_prod = float(prod.sum())
            for idr_th in [0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60]:
                for dl_th in [0.0, 0.25, 0.50, 1.00]:
                    q1 = (idr >= idr_th) & (dl > dl_th)
                    q2 = (idr < idr_th) & (dl > dl_th)
                    low_mid = q1 & inc.isin(["Low income", "Lower middle income"])
                    idr_rows.append({
                        "scenario": scen,
                        "hazard": haz,
                        "idr_threshold": idr_th,
                        "delta_loss_threshold_pct": dl_th,
                        "n_countries_total": int(len(dl)),
                        "q1_high_idr_rising_loss_count": int(q1.sum()),
                        "q1_country_share_pct": float(q1.mean() * 100.0),
                        "q1_prod_share_pct": float(prod[q1].sum() / tot_prod * 100.0),
                        "q1_low_lowermid_income_count": int(low_mid.sum()),
                        "q2_breadbasket_rising_loss_count": int(q2.sum()),
                        "q2_prod_share_pct": float(prod[q2].sum() / tot_prod * 100.0),
                    })
    pd.DataFrame(idr_rows).to_csv(os.path.join(OUT_DIR, "r8_idr_cutoff_sweep.csv"), index=False, encoding="utf-8-sig")

    # Subperiod FAOSTAT IDR stability check (2014-2018 vs 2019-2023 vs 2014-2023)
    trade_raw = pd.read_csv(r"./data/external/FAOSTAT_data_en_3-8-2026.csv")
    prod_raw = pd.read_csv(r"./data/external/FAOSTAT_data_en_3-8-2026 (1).csv")
    items_tr = {"Maize (corn)", "Rice", "Wheat", "Soybeans"}
    items_pr = {"Maize (corn)", "Rice", "Wheat", "Soybeans", "Soya beans"}
    tr_f = trade_raw[trade_raw["Item"].isin(items_tr)].copy()
    pr_f = prod_raw[prod_raw["Item"].isin(items_pr) & (prod_raw["Element"] == "Production")].copy()

    def _idr_by_area(yr_lo: int, yr_hi: int) -> pd.Series:
        t_sub = tr_f[(tr_f["Year"] >= yr_lo) & (tr_f["Year"] <= yr_hi)]
        p_sub = pr_f[(pr_f["Year"] >= yr_lo) & (pr_f["Year"] <= yr_hi)]
        imp = t_sub[t_sub["Element"].str.contains("Import", case=False, na=False)].groupby(["Area", "Year"])["Value"].sum()
        exp = t_sub[t_sub["Element"].str.contains("Export", case=False, na=False)].groupby(["Area", "Year"])["Value"].sum()
        prd = p_sub.groupby(["Area", "Year"])["Value"].sum()
        df_m = pd.DataFrame({"P": prd, "I": imp, "E": exp}).fillna(0.0)
        df_g = df_m.groupby("Area")[["P", "I", "E"]].mean()
        denom = df_g["P"] + df_g["I"] - df_g["E"]
        idr_s = np.where(denom > 0, np.clip(df_g["I"] / denom, 0.0, 1.0), np.nan)
        return pd.Series(idr_s, index=df_g.index).dropna()

    idr_1418 = _idr_by_area(2014, 2018)
    idr_1923 = _idr_by_area(2019, 2023)
    idr_1423 = _idr_by_area(2014, 2023)
    common_a = idr_1418.index.intersection(idr_1923.index).intersection(idr_1423.index)
    s1, s2, s_all = idr_1418.loc[common_a], idr_1923.loc[common_a], idr_1423.loc[common_a]
    r_pear = float(stats.pearsonr(s1, s2)[0])
    r_spear = float(stats.spearmanr(s1, s2)[0])
    agree_04 = float(((s1 >= 0.4) == (s2 >= 0.4)).mean() * 100.0)
    jacc_04 = float(((s1 >= 0.4) & (s2 >= 0.4)).sum() / max(((s1 >= 0.4) | (s2 >= 0.4)).sum(), 1) * 100.0)
    pd.DataFrame([{
        "n_countries": int(len(common_a)),
        "pearson_r_2014_2018_vs_2019_2023": r_pear,
        "spearman_rho_2014_2018_vs_2019_2023": r_spear,
        "classification_agreement_idr04_pct": agree_04,
        "jaccard_overlap_high_idr04_pct": jacc_04,
        "high_idr04_count_2014_2018": int((s1 >= 0.4).sum()),
        "high_idr04_count_2019_2023": int((s2 >= 0.4).sum()),
        "high_idr04_count_2014_2023": int((s_all >= 0.4).sum()),
    }]).to_csv(os.path.join(OUT_DIR, "r8_idr_temporal_persistence.csv"), index=False, encoding="utf-8-sig")
    print(f"  [R8 IDR Persistence 2014-2018 vs 2019-2023]: n={len(common_a)}, r={r_pear:.3f}, "
          f"rho={r_spear:.3f}, agreement@0.4={agree_04:.1f}%")
    print(f"Part 1 completed cleanly in {time.time() - t0:.1f}s!")


if __name__ == "__main__":
    main()
