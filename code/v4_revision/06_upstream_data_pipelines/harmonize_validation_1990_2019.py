"""Harmonize ALL Multi-Scale Validation, Disaster-Catalog, and GGCMI3 Benchmark Tables to 1990-2019.

Strictly eliminates all pre-1990 years (no 1980s years anywhere) while retaining:
1. 100% of cultivated croplands (both rainfed and irrigated, area-weighted).
2. Centred 5-year moving-average baseline ((MA5 - Y)/MA5 * 100) computed strictly within >=1990 windows:
   - Subnational Hultgren et al. (2025) CIL archive: 1990-2015
   - National FAOSTAT yield archive: 1990-2019
   - Global 0.5-deg GDHY v1.2/v1.3 & ISIMIP3a/GGCMI Phase 3 archive: 1990-2016
   - US Maize USDA RMA cause-of-loss insurance archive: 1990-2015 (USDA NASS: 1990-2016)
   - 18 Documented Landmark Disasters: 1990-2019 (D01/D02 = 2002 US maize & soybean droughts; 30-yr 1990-2019 ranks)
   - Global Geocoded EM-DAT/GDIS Disaster Catalog: 1990-2016
   - MODIS Global Flood Database (Tellman et al., 2021): 2000-2016 (anomalies relative to 1990-2016)

Outputs -> ./data/Software/Gemini/drought_waterlogging/results/honest_validation_1990_2019/
"""
from __future__ import annotations

import glob
import os
import sys
import time
import numpy as np
import pandas as pd
import rioxarray as rxr
import xarray as xr
from scipy import stats

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

BASE_DIR = r"."
CACHE_DIR = os.path.join(BASE_DIR, "cache")
RES_DIR = os.path.join(BASE_DIR, "results")
OUT_DIR = os.path.join(RES_DIR, "honest_validation_1990_2019")
os.makedirs(OUT_DIR, exist_ok=True)

CODEX_CIL_DIR = r"D:\Software\codex\drought&waterlogging\source_index\validation_water_stress_rebuilt\adm1\1980_2019"
CWATM_NC = os.path.join(CACHE_DIR, "cwatm_hist_4crops_1980_2019_05deg.nc")
COUNTRY_GRID_CACHE = os.path.join(CACHE_DIR, "grid_country_iso3_05deg.npz")
GDHY_DIR = r"./data/validation/gdhy_v1.3_1990_2016"
GGCMI_DIR = r"./data/validation/ggcmi3_isimip3a_yield"
PREC_NC = r"./data/cwatm_outputs/GSWP3-W5E5/Precipitation_monthtot.nc"
FRAC_DIR = r"./data/cwatm_inputs_figshare"
CELLAREA_FP = r"./data/external/cellarea.nc"
GDIS_CSV = r"./data/validation/gdis_1990_2016/pend-gdis-1960-2018-disasterlocations-csv/pend-gdis-1960-2018-disasterlocations.csv"
GFD_CSV = r"./data/validation/modis_gfd_2000_2016/event_catalog.csv"

CROPS = ["maiz", "soyb", "rice", "whea"]
GDHY_FOLDER = {"maiz": "maize", "soyb": "soybean", "rice": "rice", "whea": "wheat"}
YEARS_90_16 = np.arange(1990, 2017, dtype=int)  # 27 years within 1990-2019
YEARS_90_19 = np.arange(1990, 2020, dtype=int)  # 30 years

GGCMI_WATERLOGGING_META = {
    "acea": ("FAO-AquaCrop", "Partial (anaerobiosis point parameterization; limited ponding/drainage)"),
    "crover": ("CROVER", "No explicit root-zone O2/waterlogging yield penalty"),
    "cygma1p74": ("CYGMA 1.74", "Empirical excess-wetness / N-leaching stress representation"),
    "epic-iiasa": ("EPIC-IIASA", "Aeration stress above field capacity (water-filled pore space), no surface ponding mortality"),
    "isam": ("ISAM", "Soil moisture saturation / root respiration constraint in land-surface scheme"),
    "ldndc": ("LandscapeDNDC", "Explicit soil O2/redox biogeochemistry & N-loss under saturation"),
    "lpjml": ("LPJmL", "No explicit aeration/waterlogging transpiration penalty in standard ISIMIP3a"),
    "pdssat": ("pDSSAT", "Aeration/flooding factor in root zone; paddy bunding for rice"),
    "pepic": ("PEPIC", "EPIC aeration factor above field capacity, often offset by N/water relief in coarse grids"),
    "promet": ("PROMET", "Mechanistic soil-plant hydraulic & aeration coupling"),
    "simplace-lintul5": ("SIMPLACE-LINTUL5", "Water-deficit (transpiration reduction) only; no excess-water penalty"),
}


def normalize_admin_id_series(s: pd.Series) -> pd.Series:
    out = s.astype("string").str.strip()
    out = out.str.replace(r"\.0$", "", regex=True)
    missing = out.isna() | out.str.lower().isin(["", "nan", "none", "<na>", "null"])
    return out.mask(missing, pd.NA)


def compute_ma5_loss_1d(y: np.ndarray) -> np.ndarray:
    n = len(y)
    out = np.full(n, np.nan, dtype=np.float64)
    for t in range(n):
        lo = max(0, t - 2)
        hi = min(n, t + 3)
        win = y[lo:hi]
        ok = np.isfinite(win) & (win > 0)
        if ok.sum() >= 3 and np.isfinite(y[t]) and y[t] > 0:
            base = float(np.mean(win[ok]))
            if base > 1e-6:
                out[t] = (base - y[t]) / base * 100.0
    return out


def ma5_loss_vec(arr: np.ndarray, min_total_yrs: int = 12) -> tuple[np.ndarray, np.ndarray]:
    """5-yr centered moving average net & non-negative yield loss (%) over (T, lat, lon)."""
    y = np.where(np.isfinite(arr) & (arr > 0) & (arr < 1e8), arr.astype(np.float64), np.nan)
    valid = np.isfinite(y)
    tot_v = valid.sum(axis=0) >= min_total_yrs
    yz = np.where(valid, y, 0.0)
    nt = y.shape[0]
    net_loss = np.full(y.shape, np.nan, dtype=np.float32)
    for t in range(nt):
        t0 = max(0, t - 2)
        t1 = min(nt, t + 3)
        win_len = t1 - t0
        c = valid[t0:t1].sum(axis=0)
        ok = tot_v & valid[t] & (c == win_len)
        s = yz[t0:t1].sum(axis=0)
        base = np.where(ok, s / np.maximum(c, 1), np.nan)
        good = ok & (base > 0.05)
        net_loss[t, good] = np.clip(((base[good] - y[t, good]) / base[good]) * 100.0, -150.0, 100.0).astype(np.float32)
    pos_loss = np.where(np.isfinite(net_loss), np.maximum(net_loss, 0.0), np.nan).astype(np.float32)
    return net_loss, pos_loss


def ma5_loss_3d_min3(arr: np.ndarray) -> np.ndarray:
    """5-yr centered MA net yield loss (%) allowing >=3 valid years in 5-yr window (for landmark event cells)."""
    nt = arr.shape[0]
    valid = np.isfinite(arr) & (arr > 0)
    yz = np.where(valid, arr.astype(np.float64), 0.0)
    out = np.full_like(arr, np.nan, dtype=np.float32)
    for t in range(nt):
        lo, hi = max(0, t - 2), min(nt, t + 3)
        cnt = valid[lo:hi].sum(axis=0)
        s = yz[lo:hi].sum(axis=0)
        ok = (cnt >= 3) & valid[t]
        base = np.where(ok, s / np.maximum(cnt, 1), np.nan)
        ok2 = ok & (base > 1e-5)
        out[t] = np.where(ok2, (base - arr[t]) / base * 100.0, np.nan)
    return out


def detrend_rel_vec(arr: np.ndarray, years: np.ndarray, min_yrs: int = 12) -> np.ndarray:
    y = np.where(np.isfinite(arr) & (arr > 0) & (arr < 1e8), arr.astype(np.float64), np.nan)
    valid = np.isfinite(y)
    n = valid.sum(axis=0)
    t = (years - years.mean()).astype(np.float64)[:, None, None]
    t_bar = np.divide(np.where(valid, t, 0.0).sum(axis=0), n, out=np.zeros_like(n, dtype=np.float64), where=n >= min_yrs)
    y_bar = np.divide(np.where(valid, y, 0.0).sum(axis=0), n, out=np.zeros_like(n, dtype=np.float64), where=n >= min_yrs)
    dt = np.where(valid, t - t_bar[None, :, :], 0.0)
    dy = np.where(valid, y - y_bar[None, :, :], 0.0)
    denom = (dt * dt).sum(axis=0)
    slope = np.divide((dt * dy).sum(axis=0), denom, out=np.zeros_like(denom), where=(n >= min_yrs) & (denom > 1e-6))
    yhat = (y_bar - slope * t_bar)[None, :, :] + slope[None, :, :] * t
    ok = valid & (n[None, :, :] >= min_yrs) & (yhat > 0.05)
    out = np.full(arr.shape, np.nan, dtype=np.float32)
    out[ok] = np.clip(((y[ok] - yhat[ok]) / yhat[ok]) * 100.0, -100.0, 200.0).astype(np.float32)
    return out


def load_gdhy_crop(crop: str, years: np.ndarray) -> np.ndarray:
    sub = GDHY_FOLDER[crop]
    out = np.full((len(years), 360, 720), np.nan, dtype=np.float32)
    for ti, yr in enumerate(years):
        fp = os.path.join(GDHY_DIR, sub, f"yield_{yr}.nc4")
        ds = xr.open_dataset(fp)
        v = np.roll(ds["var"].values[::-1, :], 360, axis=1)
        ds.close()
        out[ti] = np.where(np.isfinite(v) & (v > 0), v, np.nan)
    return out


def load_ggcmi_weights_pair(crop: str, cellarea_ha: np.ndarray) -> dict:
    def rd(fn: str) -> np.ndarray:
        v = rxr.open_rasterio(os.path.join(FRAC_DIR, fn)).squeeze().values.astype(np.float64)
        return np.where(np.isfinite(v) & (v > 0), v, 0.0) * cellarea_ha

    if crop == "maiz":
        return {"firr": rd("MIRCA-OS_Maize_2000_ir_30arcmin.tif"), "noirr": rd("MIRCA-OS_Maize_2000_rf_30arcmin.tif")}
    if crop == "soyb":
        return {"firr": rd("MIRCA-OS_Soybeans_2000_ir_30arcmin.tif"), "noirr": rd("MIRCA-OS_Soybeans_2000_rf_30arcmin.tif")}
    if crop == "rice":
        return {"firr": rd("MIRCA-OS_Rice_2000_ir_30arcmin.tif"), "noirr": rd("MIRCA-OS_Rice_2000_rf_30arcmin.tif")}
    if crop == "whea":
        return {
            "wwh_firr": rd("MIRCA-OS_Wheat_2000_irr_30arcmin_winter.tif"),
            "wwh_noirr": rd("MIRCA-OS_Wheat_2000_rf_30arcmin_winter.tif"),
            "swh_firr": rd("MIRCA-OS_Wheat_2000_irr_30arcmin_spring.tif"),
            "swh_noirr": rd("MIRCA-OS_Wheat_2000_rf_30arcmin_spring.tif"),
        }
    raise ValueError(crop)


def load_ggcmi_yield_raw(model: str, crop: str, years: np.ndarray, m_areas: dict) -> np.ndarray | None:
    def rd_nc(code: str, mode: str) -> np.ndarray | None:
        pat = os.path.join(GGCMI_DIR, f"{model}_gswp3-w5e5_obsclim_2015soc_default_yield-{code}-{mode}_global_annual-gs_1901_*.nc")
        matches = glob.glob(pat)
        if not matches:
            return None
        ds = xr.open_dataset(matches[0], decode_times=False)
        vn = [k for k in ds.data_vars if "yield" in k][0]
        t_yrs = 1901 + np.rint(ds["time"].values.astype(float)).astype(int)
        out_arr = np.full((len(years), ds.sizes["lat"], ds.sizes["lon"]), np.nan, dtype=np.float32)
        for yi, y in enumerate(years):
            w_idx = np.where(t_yrs == y)[0]
            if len(w_idx):
                out_arr[yi] = ds[vn].isel(time=int(w_idx[0])).values.astype(np.float32)
        ds.close()
        return np.where(np.isfinite(out_arr) & (out_arr > 0) & (out_arr < 1e8), out_arr, np.nan)

    if crop in ("maiz", "soyb", "rice"):
        code = {"maiz": "mai", "soyb": "soy", "rice": "ri1"}[crop]
        y_ir = rd_nc(code, "firr")
        y_rf = rd_nc(code, "noirr")
        if y_ir is None or y_rf is None:
            return None
        a_ir, a_rf = m_areas["firr"][None, :, :], m_areas["noirr"][None, :, :]
        num = np.where(np.isfinite(y_ir), y_ir * a_ir, 0.0) + np.where(np.isfinite(y_rf), y_rf * a_rf, 0.0)
        den = np.where(np.isfinite(y_ir), a_ir, 0.0) + np.where(np.isfinite(y_rf), a_rf, 0.0)
        return np.where(den > 0, num / np.maximum(den, 1e-12), np.nan).astype(np.float32)

    parts = []
    for subc in ("wwh", "swh"):
        for mode in ("firr", "noirr"):
            y_sub = rd_nc(subc, mode)
            if y_sub is None:
                return None
            parts.append((y_sub, m_areas[f"{subc}_{mode}"][None, :, :]))
    num = sum(np.where(np.isfinite(y_s), y_s * a_s, 0.0) for y_s, a_s in parts)
    den = sum(np.where(np.isfinite(y_s), a_s, 0.0) for y_s, a_s in parts)
    return np.where(den > 0, num / np.maximum(den, 1e-12), np.nan).astype(np.float32)


def load_precip_z_and_mean(years: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ds = xr.open_dataset(PREC_NC, decode_times=False)
    t_months = np.rint(ds["time"].values.astype(float)).astype(int)
    yrs_all = 1901 + (t_months // 12)
    v = ds["Precipitation_monthtot"].values
    ds.close()
    ann_p = np.zeros((len(years), v.shape[1], v.shape[2]), dtype=np.float32)
    for i, yr in enumerate(years):
        ann_p[i] = np.nansum(v[yrs_all == yr], axis=0)
    if float(np.nanmedian(ann_p[ann_p > 0])) < 20.0:
        ann_p = ann_p * 1000.0
    p_mean = np.nanmean(ann_p, axis=0)
    p_sd = np.nanstd(ann_p, axis=0)
    z_p = np.divide(ann_p - p_mean[None, :, :], p_sd[None, :, :], out=np.full_like(ann_p, np.nan), where=p_sd[None, :, :] > 1.0)
    return z_p, p_mean


def wmean(val: np.ndarray, w: np.ndarray, mask: np.ndarray) -> float:
    ok = mask & np.isfinite(val) & (w > 0)
    if not np.any(ok):
        return np.nan
    return float(np.sum(val[ok] * w[ok]) / np.sum(w[ok]))


def eval_wet_subset(obs: np.ndarray, mod: np.ndarray, w: np.ndarray, mask: np.ndarray) -> dict:
    ok = mask & np.isfinite(obs) & np.isfinite(mod) & (w > 0)
    n = int(ok.sum())
    if n < 10:
        return {"n": n, "obs_mean": np.nan, "mod_mean": np.nan, "bias": np.nan, "rmse": np.nan, "mae": np.nan, "hit_neg_given_obs_neg": np.nan}
    ov = obs[ok].astype(np.float64)
    mv = mod[ok].astype(np.float64)
    wv = w[ok].astype(np.float64)
    wn = wv / wv.sum()
    obs_m = float(np.sum(wn * ov))
    mod_m = float(np.sum(wn * mv))
    err = mv - ov
    rmse = float(np.sqrt(np.sum(wn * (err ** 2))))
    mae = float(np.sum(wn * np.abs(err)))
    neg_o = ov < 0
    hit_neg = float(np.sum(wv[neg_o & (mv < 0)]) / np.sum(wv[neg_o])) if neg_o.any() else np.nan
    return {"n": n, "obs_mean": obs_m, "mod_mean": mod_m, "bias": mod_m - obs_m, "rmse": rmse, "mae": mae, "hit_neg_given_obs_neg": hit_neg}


def main() -> None:
    t0 = time.time()

    # =========================================================================
    # 1. SUBNATIONAL HULTGREN ET AL. (2025) CIL BENCHMARK (1990-2015)
    # =========================================================================
    print("[1/6] Building Subnational CIL (1990-2015) maps & metrics...")
    pts = pd.read_csv(os.path.join(CODEX_CIL_DIR, "cil_cwatm_validation_points_adm1.csv"))
    pts["iso"] = pts["iso"].astype(str).str.upper().str.strip()
    pts["adm1_id"] = normalize_admin_id_series(pts["adm1_id"])
    pts_90 = pts[(pts["year"] >= 1990) & (pts["year"] <= 2019) & pts["obs_loss_pct"].notna() & pts["sim_hazard_loss_pct"].notna()].copy()

    df_map = (
        pts_90.groupby(["hazard", "crop", "iso", "adm1_id"], as_index=False)
        .agg(
            n_years=("year", "count"),
            mean_obs_loss_pct=("obs_loss_pct", "mean"),
            mean_sim_loss_pct=("sim_hazard_loss_pct", "mean"),
        )
    )
    df_map["mean_loss_bias_pctpt"] = df_map["mean_sim_loss_pct"] - df_map["mean_obs_loss_pct"]
    df_map.to_csv(os.path.join(OUT_DIR, "s13_adm1_bias_map_source.csv"), index=False)

    met_rows = []
    for hz in ["dry", "wet"]:
        for cr in CROPS:
            g = pts_90[(pts_90["hazard"] == hz) & (pts_90["crop"] == cr)]
            obs = g["obs_loss_pct"].values.astype(float)
            sim = g["sim_hazard_loss_pct"].values.astype(float)
            diff = sim - obs
            r_val, _ = stats.pearsonr(sim, obs) if len(obs) >= 3 else (np.nan, np.nan)
            u_sub = df_map[(df_map["hazard"] == hz) & (df_map["crop"] == cr)]["mean_loss_bias_pctpt"].dropna()
            met_rows.append({
                "hazard": hz,
                "crop": cr,
                "y_col": "sim_hazard_loss_pct",
                "n": int(len(g)),
                "pearson_r": float(r_val),
                "rmse_pct": float(np.sqrt(np.mean(diff ** 2))),
                "bias_pct": float(np.mean(diff)),
                "obs_mean_loss_pct": float(np.mean(obs)),
                "sim_mean_loss_pct": float(np.mean(sim)),
                "mae_pct": float(np.mean(np.abs(diff))),
                "neg_bias_event_share_pct": float(np.mean(diff < 0) * 100.0),
                "n_adm1_units": int(len(u_sub)),
                "adm1_mean_bias_pctpt": float(u_sub.mean()),
                "adm1_neg_bias_unit_share_pct": float((u_sub < 0).mean() * 100.0),
            })
    df_cil_met = pd.DataFrame(met_rows)
    df_cil_met.to_csv(os.path.join(OUT_DIR, "s13_adm1_crop_hazard_metrics.csv"), index=False)

    print(f"  CIL 1990-{int(pts_90['year'].max())}: total N = {len(pts_90):,}, countries = {pts_90['iso'].nunique()}, "
          f"drought N = {int((pts_90['hazard']=='dry').sum()):,}, wet N = {int((pts_90['hazard']=='wet').sum()):,}")
    print(f"  CIL pooled obs = {pts_90['obs_loss_pct'].mean():.2f}%, sim = {pts_90['sim_hazard_loss_pct'].mean():.2f}%, "
          f"bias = {(pts_90['sim_hazard_loss_pct'] - pts_90['obs_loss_pct']).mean():+.2f}%, "
          f"neg-bias unit-crop-hazard share = {(df_map['mean_loss_bias_pctpt'] < 0).mean()*100:.1f}%")
    print(df_cil_met[["hazard", "crop", "n", "n_adm1_units", "obs_mean_loss_pct", "sim_mean_loss_pct", "bias_pct", "rmse_pct", "mae_pct", "adm1_mean_bias_pctpt"]].round(2).to_string(index=False))

    # =========================================================================
    # 2. US MAIZE FULL DOMAIN (RMA 1990-2015 & NASS 1990-2016 & GGCMI3 1990-2016)
    # =========================================================================
    print("\n[2/6] Building US Maize Full-Domain (1990-2015 RMA / 1990-2016 NASS & GGCMI3)...")
    us_panel = pd.read_csv(os.path.join(RES_DIR, "tier1_us_rma_li2019", "tier1_us_grid_year_panel.csv"))
    us_panel = us_panel[(us_panel["Year"] >= 1990) & (us_panel["Year"] <= 2016)].sort_values(["grid_id", "Year"]).reset_index(drop=True)

    # Recompute 5-yr centered MA strictly within 1990-2016
    ma5_loss = np.full(len(us_panel), np.nan, dtype=np.float64)
    for _, grp in us_panel.groupby("grid_id", sort=False):
        ma5_loss[grp.index.values] = compute_ma5_loss_1d(grp["Yield"].values.astype(float))
    us_panel["obs_ma5_loss_pct"] = ma5_loss
    us_panel["obs_ma5_pos_loss_pct"] = np.maximum(ma5_loss, 0.0)

    # Recompute GGCMI3 11-model 5-yr MA net loss strictly over 1990-2016 for US maize
    cellarea_ha = xr.open_dataset(CELLAREA_FP)["cellarea30min"].values.astype(np.float64) / 1e4
    m_maiz = load_ggcmi_weights_pair("maiz", cellarea_ha)
    gg_maiz_nets = []
    for m in sorted(GGCMI_WATERLOGGING_META.keys()):
        ym = load_ggcmi_yield_raw(m, "maiz", YEARS_90_16, m_maiz)
        if ym is not None:
            gn, _ = ma5_loss_vec(ym, min_total_yrs=12)
            gg_maiz_nets.append(gn)
    gg_maiz_ens_net = np.nanmean(np.stack(gg_maiz_nets, axis=0), axis=0)  # (27, 360, 720), >0 = loss, <0 = gain
    ti_us = us_panel["Year"].values.astype(int) - 1990
    us_panel["ggcmi_ma5_net_loss_pct"] = gg_maiz_ens_net[ti_us, us_panel["row"].values.astype(int), us_panel["col"].values.astype(int)]

    sub_rma = us_panel[(us_panel["Year"] >= 1990) & (us_panel["Year"] <= 2015)].dropna(subset=["rma_drought", "rma_excess_wet"]).copy()
    ann_rows = []
    for yr, grp in sub_rma.groupby("Year"):
        w = np.where(np.isfinite(grp["Area_acres"].values) & (grp["Area_acres"].values > 0), grp["Area_acres"].values, 1.0)
        def wm_col(col: str, g=grp, wt=w) -> float:
            v = g[col].values.astype(float)
            m = np.isfinite(v) & np.isfinite(wt) & (wt > 0)
            return float(np.average(v[m], weights=wt[m])) if m.any() else np.nan
        ann_rows.append({
            "Year": int(yr),
            "n_grids": int(len(grp)),
            "cwatm_loss_dry_pct": wm_col("cwatm_loss_dry_pct"),
            "cwatm_loss_wet_pct": wm_col("cwatm_loss_wet_pct"),
            "rma_drought": wm_col("rma_drought"),
            "rma_excess_wet": wm_col("rma_excess_wet"),
            "obs_ma5_loss_pct": wm_col("obs_ma5_loss_pct"),
        })
    df_ann = pd.DataFrame(ann_rows)
    df_ann.to_csv(os.path.join(OUT_DIR, "s14_us_full_annual_1990_2015.csv"), index=False)

    r_dry, p_dry = stats.pearsonr(df_ann["cwatm_loss_dry_pct"], df_ann["rma_drought"])
    r_wet, p_wet = stats.pearsonr(df_ann["cwatm_loss_wet_pct"], df_ann["rma_excess_wet"])
    r_dw, p_dw = stats.pearsonr(df_ann["cwatm_loss_wet_pct"], df_ann["rma_drought"])
    r_wd, p_wd = stats.pearsonr(df_ann["cwatm_loss_dry_pct"], df_ann["rma_excess_wet"])
    print(f"  US Maize (1990-2015 RMA grids={sub_rma['grid_id'].nunique()}, 1990-2016 NASS grids={us_panel['grid_id'].nunique()}):")
    print(f"    Drought r = {r_dry:.4f} (P = {p_dry:.2e}), Excess moisture r = {r_wet:.4f} (P = {p_wet:.2e})")
    print(f"    Off-diagonal: r_dw = {r_dw:.4f} (P = {p_dw:.2f}), r_wd = {r_wd:.4f} (P = {p_wd:.2f})")

    bins = [-np.inf, -2.0, -1.5, -1.0, -0.5, 0.5, 1.0, 1.5, 2.0, np.inf]
    lbls = ["<-2.0", "[-2.0,-1.5)", "[-1.5,-1.0)", "[-1.0,-0.5)", "[-0.5,0.5]", "(0.5,1.0]", "(1.0,1.5]", "(1.5,2.0]", ">2.0"]
    us_panel["pbin"] = pd.cut(us_panel["Prec_to_sd"], bins=bins, labels=lbls)
    pbin_rows = []
    for lb in lbls:
        grp = us_panel[us_panel["pbin"] == lb]
        w = np.where(np.isfinite(grp["Area_acres"].values) & (grp["Area_acres"].values > 0), grp["Area_acres"].values, 1.0)
        def wm(col: str, g=grp, wt=w) -> float:
            v = g[col].values.astype(float)
            m = np.isfinite(v) & np.isfinite(wt) & (wt > 0)
            return float(np.average(v[m], weights=wt[m])) if m.any() else np.nan
        grp_rma = grp[(grp["Year"] >= 1990) & (grp["Year"] <= 2015)].dropna(subset=["rma_drought", "rma_excess_wet"])
        w_rma = np.where(np.isfinite(grp_rma["Area_acres"].values) & (grp_rma["Area_acres"].values > 0), grp_rma["Area_acres"].values, 1.0)
        pbin_rows.append({
            "prec_bin": lb,
            "prec_sd_mean": wm("Prec_to_sd"),
            "n_grid_years_nass": int(grp["obs_ma5_loss_pct"].notna().sum()),
            "n_grid_years_rma": int(grp_rma["rma_drought"].notna().sum()),
            "obs_ma5_loss_pct": wm("obs_ma5_loss_pct"),
            "cwatm_loss_dry_pct": wm("cwatm_loss_dry_pct"),
            "cwatm_loss_wet_pct": wm("cwatm_loss_wet_pct"),
            "cwatm_loss_tot_pct": wm("cwatm_loss_dry_pct") + wm("cwatm_loss_wet_pct"),
            "rma_drought": wm("rma_drought", grp_rma, w_rma),
            "rma_excess_wet": wm("rma_excess_wet", grp_rma, w_rma),
        })
    df_pbin = pd.DataFrame(pbin_rows)
    df_pbin.to_csv(os.path.join(OUT_DIR, "s14_us_full_precip_bins.csv"), index=False)
    print(df_pbin[["prec_bin", "prec_sd_mean", "obs_ma5_loss_pct", "cwatm_loss_dry_pct", "cwatm_loss_wet_pct", "rma_drought", "rma_excess_wet"]].round(2).to_string(index=False))

    # Build ma5_us_cornbelt.csv (1990-2016) for Fig. S16a
    us_cb_rows = []
    for rname, msk in [
        ("Extreme Dry Precip (Prec <= -1.5 SD)", us_panel["Prec_to_sd"] <= -1.5),
        ("Extreme Wet Precip (Prec >= +2.0 SD)", us_panel["Prec_to_sd"] >= 2.0),
        ("Extreme Wet Precip (Prec >= +2.5 SD)", us_panel["Prec_to_sd"] >= 2.5),
        ("CWatM Drought Event (dry>=7d, wet==0)", (us_panel["cwatm_dry_days"] >= 7) & (us_panel["cwatm_wet_days"] < 1)),
        ("CWatM Severe Drought (dry_loss>=5%)", us_panel["cwatm_loss_dry_pct"] >= 5.0),
        ("CWatM Waterlogging Event (wet>=7d, dry==0)", (us_panel["cwatm_wet_days"] >= 7) & (us_panel["cwatm_dry_days"] < 1)),
        ("CWatM Severe Waterlog (wet_loss>=3%)", us_panel["cwatm_loss_wet_pct"] >= 3.0),
        ("1993 Midwest Great Flood (Prec >= +2.0 SD)", (us_panel["Year"] == 1993) & (us_panel["Prec_to_sd"] >= 2.0)),
        ("1993 Midwest Great Flood (all Corn Belt)", us_panel["Year"] == 1993),
        ("2012 US Flash Drought (all Corn Belt)", us_panel["Year"] == 2012),
        ("2002 US Drought (all Corn Belt)", us_panel["Year"] == 2002),
    ]:
        sub = us_panel[msk & us_panel["obs_ma5_loss_pct"].notna()]
        def wm_s(c: str) -> float:
            ok = sub[c].notna() & (sub["Area_acres"] > 0)
            return float(np.average(sub.loc[ok, c], weights=sub.loc[ok, "Area_acres"])) if ok.any() else np.nan
        us_cb_rows.append({
            "Regime": rname,
            "n": int(len(sub)),
            "Obs_MA5_Net%": wm_s("obs_ma5_loss_pct"),
            "Obs_MA5_Pos%": wm_s("obs_ma5_pos_loss_pct"),
            "GGCMI3_Net%": wm_s("ggcmi_ma5_net_loss_pct"),
            "CWatM_Dry_Loss%": wm_s("cwatm_loss_dry_pct"),
            "CWatM_Wet_Loss%": wm_s("cwatm_loss_wet_pct"),
            "RMA_Dry_%": wm_s("rma_drought"),
            "RMA_Wet_%": wm_s("rma_excess_wet"),
        })
    df_us_cb = pd.DataFrame(us_cb_rows)
    df_us_cb.to_csv(os.path.join(OUT_DIR, "ma5_us_cornbelt.csv"), index=False)
    print("\n  US Maize Extreme Regimes (1990-2016):")
    print(df_us_cb.round(2).to_string(index=False))

    # =========================================================================
    # 3. GLOBAL 0.5-DEG GDHY DOSE-RESPONSE & GGCMI3 / ABLATION (1990-2016)
    # =========================================================================
    print("\n[3/6] Building Global 0.5-deg GDHY Dose-Response, Ablation & 11-Model GGCMI3 (1990-2016)...")
    cw_all = xr.open_dataset(CWATM_NC)
    cw_90_16 = cw_all.sel(year=YEARS_90_16)
    cw_90_19 = cw_all.sel(year=YEARS_90_19)

    z_p, p_mean = load_precip_z_and_mean(YEARS_90_16)
    humid_2d = p_mean >= 800.0
    humid_3d = np.broadcast_to(humid_2d[None, :, :], (len(YEARS_90_16), 360, 720))

    bins_spec = [
        ("<0.1", 0.0, 0.1),
        ("0.1-1", 0.1, 1.0),
        ("1-3", 1.0, 3.0),
        ("3-6", 3.0, 6.0),
        (">6", 6.0, 101.0),
    ]
    dr_rows = []
    rows_global = []
    ablation_rows = []
    ggcmi_rows = []
    gdhy_raw_cache: dict[str, np.ndarray] = {}
    gdhy_pos_cache: dict[str, np.ndarray] = {}

    all_models = sorted(GGCMI_WATERLOGGING_META.keys())

    for crop in CROPS:
        gd_y = load_gdhy_crop(crop, YEARS_90_16)
        gdhy_raw_cache[crop] = gd_y
        obs_net, obs_pos = ma5_loss_vec(gd_y, min_total_yrs=12)
        gdhy_pos_cache[crop] = obs_pos

        c_sub = cw_90_16.sel(crop=crop)
        area_2d = c_sub["area_ha"].values.astype(np.float64)
        w_3d = np.broadcast_to(area_2d[None, :, :], (len(YEARS_90_16), 360, 720))
        ld = c_sub["loss_dry_pct"].values
        lw = c_sub["loss_wet_pct"].values
        dd = c_sub["dry_days"].values
        wd = c_sub["wet_days"].values
        cf_dry_thr7 = c_sub["loss_dry_thr7_pct"].values
        cf_wet_thr7 = c_sub["loss_wet_thr7_pct"].values

        # Fig. S14e,f Dose-response
        for haz, larr in [("drought", ld), ("waterlogging", lw)]:
            for b_lbl, lo, hi in bins_spec:
                msk = np.isfinite(obs_net) & (w_3d > 0) & (larr >= lo) & (larr < hi)
                if haz == "waterlogging":
                    msk = msk & humid_3d
                n_pts = int(msk.sum())
                if n_pts > 10:
                    vals = obs_net[msk].astype(np.float64)
                    wts = w_3d[msk]
                    wm_v = float(np.sum(vals * wts) / np.sum(wts))
                    se_v = float(np.std(vals) / np.sqrt(n_pts))
                else:
                    wm_v, se_v = np.nan, np.nan
                dr_rows.append({
                    "crop": crop, "hazard": haz, "bin": b_lbl,
                    "n": n_pts, "obs_ma5_net_loss_pct": wm_v, "obs_se_pct": se_v,
                })

        # MA5 comparison (Fig. S16b)
        cw2s_net, cw2s_pos = ma5_loss_vec(c_sub["yield_2s"].values, min_total_yrs=12)
        cwnowl_net, _ = ma5_loss_vec(c_sub["yield_nowl"].values, min_total_yrs=12)

        m_areas = load_ggcmi_weights_pair(crop, cellarea_ha)
        gg_nets, gg_poss, gg_det_stack = [], [], []

        # Detrended anomalies for Tier 2 RMSE/ablation (Fig. S16c-f & Table S5)
        obs_anom = detrend_rel_vec(gd_y, YEARS_90_16, min_yrs=12)
        anom_2s = detrend_rel_vec(c_sub["yield_2s"].values, YEARS_90_16, min_yrs=12)
        anom_nowl = detrend_rel_vec(c_sub["yield_nowl"].values, YEARS_90_16, min_yrs=12)
        active = area_2d >= 500.0
        base_valid = active[None, :, :] & np.isfinite(obs_anom) & np.isfinite(anom_2s) & np.isfinite(anom_nowl)

        t2_regimes = [
            ("W1_PureWaterlogging_wet7_dry0", base_valid & (wd >= 7.0) & (dd < 1.0) & (lw >= 0.5)),
            ("W2_SevereWaterlogging_wetloss_ge_2pct", base_valid & (lw >= 2.0) & (dd < 7.0)),
            ("W3_SevereWaterlogging_wetloss_ge_5pct", base_valid & (lw >= 5.0) & (dd < 7.0)),
            ("W4_ModelFree_HumidExtremePrecip_Zp_gt_1p5", base_valid & (z_p > 1.5) & humid_3d),
            ("W5_ModelFree_GlobalExtremePrecip_Zp_gt_2p0", base_valid & (z_p > 2.0)),
        ]

        for reg_name, msk in t2_regimes:
            s_2s = eval_wet_subset(obs_anom, anom_2s, w_3d, msk)
            s_nw = eval_wet_subset(obs_anom, anom_nowl, w_3d, msk)
            diff_ok = msk & (np.abs(anom_2s - anom_nowl) > 1e-4)
            if diff_ok.sum() >= 20:
                ae_2s = np.abs(anom_2s[diff_ok] - obs_anom[diff_ok])
                ae_nw = np.abs(anom_nowl[diff_ok] - obs_anom[diff_ok])
                _, p_wil = stats.wilcoxon(ae_nw - ae_2s, alternative="greater")
                frac_improved = float((ae_2s < ae_nw).mean())
            else:
                p_wil, frac_improved = np.nan, np.nan
            ablation_rows.append({
                "crop": crop,
                "regime": reg_name,
                "n_grid_years": s_2s["n"],
                "n_active_cells": int(np.any(msk, axis=0).sum()),
                "obs_gdhy_mean_anom_pct": s_2s["obs_mean"],
                "cwatm_2s_mean_anom_pct": s_2s["mod_mean"],
                "cwatm_nowl_mean_anom_pct": s_nw["mod_mean"],
                "rmse_2s_pct": s_2s["rmse"],
                "rmse_nowl_pct": s_nw["rmse"],
                "rmse_reduction_pctpt": s_nw["rmse"] - s_2s["rmse"],
                "mae_2s_pct": s_2s["mae"],
                "mae_nowl_pct": s_nw["mae"],
                "mae_reduction_pctpt": s_nw["mae"] - s_2s["mae"],
                "hit_neg_2s": s_2s["hit_neg_given_obs_neg"],
                "hit_neg_nowl": s_nw["hit_neg_given_obs_neg"],
                "hit_neg_gain": s_2s["hit_neg_given_obs_neg"] - s_nw["hit_neg_given_obs_neg"],
                "frac_cells_abs_err_reduced": frac_improved,
                "wilcoxon_p_err_reduction": p_wil,
            })

        for m in all_models:
            ym = load_ggcmi_yield_raw(m, crop, YEARS_90_16, m_areas)
            if ym is None:
                continue
            gn, gp = ma5_loss_vec(ym, min_total_yrs=12)
            gg_nets.append(gn)
            gg_poss.append(gp)
            g_anom = detrend_rel_vec(ym, YEARS_90_16, min_yrs=12)
            gg_det_stack.append(g_anom)
            for reg_tag, msk in [
                ("SevereWaterlogging_wetloss_ge_2pct", t2_regimes[1][1]),
                ("ModelFree_HumidExtremePrecip_Zp_gt_1p5", t2_regimes[3][1]),
            ]:
                st = eval_wet_subset(obs_anom, g_anom, w_3d, msk)
                ggcmi_rows.append({
                    "crop": crop, "regime": reg_tag, "model": m,
                    "model_label": GGCMI_WATERLOGGING_META[m][0],
                    "process_representation": GGCMI_WATERLOGGING_META[m][1],
                    **st,
                })

        ens_anom = np.nanmean(np.stack(gg_det_stack, axis=0), axis=0)
        for reg_tag, msk in [
            ("SevereWaterlogging_wetloss_ge_2pct", t2_regimes[1][1]),
            ("ModelFree_HumidExtremePrecip_Zp_gt_1p5", t2_regimes[3][1]),
        ]:
            for m_code, m_lab, proc_desc, arr_m in [
                ("CWatM-2S", "CWatM (Two-Sided)", "Extended FAO-56 linear aeration decline above field capacity (stage-specific rice)", anom_2s),
                ("CWatM-NoWL", "CWatM (No Waterlogging)", "Drought-only counterfactual ablation (abnormal waterlogging removed)", anom_nowl),
                ("GGCMI3-Mean", "GGCMI3 Ensemble Mean", "Multi-model mean across available ISIMIP3a crop models", ens_anom),
            ]:
                st = eval_wet_subset(obs_anom, arr_m, w_3d, msk)
                ggcmi_rows.append({
                    "crop": crop, "regime": reg_tag, "model": m_code,
                    "model_label": m_lab, "process_representation": proc_desc, **st,
                })

        gg_ens_net = np.nanmean(np.stack(gg_nets, axis=0), axis=0)
        gg_ens_pos = np.nanmean(np.stack(gg_poss, axis=0), axis=0)
        common = np.isfinite(obs_net) & np.isfinite(cw2s_net) & np.isfinite(gg_ens_net) & (w_3d > 0)
        for rname, msk, htype in [
            ("Drought_7d (drystress>=7, wet==0)", common & (dd >= 7) & (wd < 1), "dry"),
            ("Drought_Severe_20d (drystress>=20)", common & (dd >= 20) & (wd < 7), "dry"),
            ("Drought_HighLoss (CWatM_dry>=5%)", common & (ld >= 5.0), "dry"),
            ("Waterlog_7d_All (wetstress>=7, dry==0)", common & (wd >= 7) & (dd < 1), "wet"),
            ("Waterlog_7d_Humid (wet>=7, P>=800mm)", common & (wd >= 7) & (dd < 1) & humid_3d, "wet"),
            ("Waterlog_HighLoss_Humid (CWatM_wet>=5%, P>=800mm)", common & (lw >= 5.0) & humid_3d, "wet"),
        ]:
            cf_col = cf_dry_thr7 if htype == "dry" else cf_wet_thr7
            rows_global.append({
                "crop": crop, "regime": rname, "n_grid_years": int(msk.sum()),
                "Obs_GDHY_NetLoss_%": wmean(obs_net, w_3d, msk),
                "GGCMI3_NetLoss_%": wmean(gg_ens_net, w_3d, msk),
                "CWatM2S_MA5_NetLoss_%": wmean(cw2s_net, w_3d, msk),
                "CWatMNoWL_MA5_NetLoss_%": wmean(cwnowl_net, w_3d, msk),
                "Obs_GDHY_PosLoss_%": wmean(obs_pos, w_3d, msk),
                "GGCMI3_PosLoss_%": wmean(gg_ens_pos, w_3d, msk),
                "CWatM2S_MA5_PosLoss_%": wmean(cw2s_pos, w_3d, msk),
                "CWatM_Counterfactual_Loss_%": wmean(cf_col, w_3d, msk),
            })

    df_dr = pd.DataFrame(dr_rows)
    df_glob = pd.DataFrame(rows_global)
    df_abl = pd.DataFrame(ablation_rows)
    df_gg = pd.DataFrame(ggcmi_rows)
    df_dr.to_csv(os.path.join(OUT_DIR, "ma5_dose_response_gdhy_4crops.csv"), index=False)
    df_glob.to_csv(os.path.join(OUT_DIR, "ma5_global_gdhy_ggcmi_cwatm_4crops.csv"), index=False)
    df_abl.to_csv(os.path.join(OUT_DIR, "tier2_cwatm_waterlogging_ablation_4crops.csv"), index=False)
    df_gg.to_csv(os.path.join(OUT_DIR, "tier2_ggcmi3_vs_cwatm_wet_extremes_4crops.csv"), index=False)

    print("  Dose-response >6% bin (1990-2016):")
    print(df_dr[df_dr["bin"] == ">6"].round(2).to_string(index=False))
    print("\n  Humid Severe Waterlogging (Fig. S16b, 1990-2016):")
    print(df_glob[df_glob["regime"] == "Waterlog_HighLoss_Humid (CWatM_wet>=5%, P>=800mm)"][["crop", "n_grid_years", "Obs_GDHY_NetLoss_%", "CWatM2S_MA5_NetLoss_%", "CWatMNoWL_MA5_NetLoss_%", "GGCMI3_NetLoss_%"]].round(2).to_string(index=False))
    print("\n  ModelFree Humid Extreme Precip Z_P > 1.5 RMSE & Wilcoxon P (Fig. S16c-f & Table S5, 1990-2016):")
    print(df_abl[df_abl["regime"] == "W4_ModelFree_HumidExtremePrecip_Zp_gt_1p5"][["crop", "n_grid_years", "rmse_2s_pct", "rmse_nowl_pct", "wilcoxon_p_err_reduction"]].to_string(index=False))
    print(df_gg[(df_gg["regime"] == "ModelFree_HumidExtremePrecip_Zp_gt_1p5") & (df_gg["model"].isin(["CWatM-2S", "CWatM-NoWL", "GGCMI3-Mean"]))][["crop", "model", "n", "rmse", "bias"]].round(2).to_string(index=False))

    # =========================================================================
    # 4. MULTI-SCALE TABLE S3 (CIL 1990-2015, FAOSTAT 1990-2019, GDHY 1990-2016)
    # =========================================================================
    print("\n[4/6] Building Multi-Scale Table S3 (1990-2019 harmonized)...")
    t3_rows = []
    for cr in CROPS:
        for hz_c, hz_name in [("dry", "drought"), ("wet", "waterlogging")]:
            r = df_cil_met[(df_cil_met["crop"] == cr) & (df_cil_met["hazard"] == hz_c)].iloc[0]
            t3_rows.append({
                "scale": "Subnational (CIL)",
                "crop": cr, "hazard": hz_name, "n": int(r["n"]),
                "obs_loss_pct": float(r["obs_mean_loss_pct"]),
                "sim_loss_pct": float(r["sim_mean_loss_pct"]),
                "bias_pctpt": float(r["bias_pct"]),
                "rmse_pct": float(r["rmse_pct"]),
                "mae_pct": float(r["mae_pct"]),
            })

    fao = pd.read_csv(os.path.join(RES_DIR, "tier3_multiscale_4crops", "tier3_panel_national_faostat.csv"))
    fao = fao[(fao["year"] >= 1990) & (fao["year"] <= 2019)].sort_values(["crop", "iso", "year"]).reset_index(drop=True)
    fao["ma5_base"] = fao.groupby(["crop", "iso"])["obs_yield"].transform(lambda s: s.rolling(5, center=True, min_periods=3).mean())
    fao["obs_pos_loss_pct"] = np.maximum(((fao["ma5_base"] - fao["obs_yield"]) / fao["ma5_base"] * 100.0).clip(-150.0, 100.0), 0.0)

    for cr in CROPS:
        sub = fao[(fao["crop"] == cr) & fao["obs_pos_loss_pct"].notna()].copy()
        m_d = (sub["dry_days"] >= 7.0) & (sub["loss_dry_pct"] > sub["loss_wet_pct"])
        m_w = ((sub["wet_days"] >= 5.0) | (sub["loss_wet_pct"] >= 1.0)) & (sub["loss_wet_pct"] > sub["loss_dry_pct"])
        for hz_name, msk, scol in [("drought", m_d, "loss_dry_pct"), ("waterlogging", m_w, "loss_wet_pct")]:
            g = sub[msk]
            ov = g["obs_pos_loss_pct"].values.astype(float)
            sv = g[scol].values.astype(float)
            diff = sv - ov
            t3_rows.append({
                "scale": "National (FAOSTAT)",
                "crop": cr, "hazard": hz_name, "n": int(len(g)),
                "obs_loss_pct": float(np.mean(ov)),
                "sim_loss_pct": float(np.mean(sv)),
                "bias_pctpt": float(np.mean(diff)),
                "rmse_pct": float(np.sqrt(np.mean(diff ** 2))),
                "mae_pct": float(np.mean(np.abs(diff))),
            })

    for cr in CROPS:
        obs_l = gdhy_pos_cache[cr]
        c_sub = cw_90_16.sel(crop=cr)
        area_2d = c_sub["area_ha"].values.astype(np.float64)
        w_3d = np.broadcast_to(area_2d[None, :, :], obs_l.shape)
        dd = c_sub["dry_days"].values
        wd = c_sub["wet_days"].values
        ld = c_sub["loss_dry_pct"].values
        lw = c_sub["loss_wet_pct"].values
        m_d = np.isfinite(obs_l) & (w_3d > 0) & (dd >= 7.0) & (wd == 0.0)
        m_w = np.isfinite(obs_l) & (w_3d > 0) & (wd >= 7.0) & (dd == 0.0)
        for hz_name, msk, s_arr in [("drought", m_d, ld), ("waterlogging", m_w, lw)]:
            ov = obs_l[msk].astype(np.float64)
            sv = s_arr[msk].astype(np.float64)
            wt = w_3d[msk]
            diff = sv - ov
            wsum = np.sum(wt)
            t3_rows.append({
                "scale": "Gridded (GDHY 0.5°)",
                "crop": cr, "hazard": hz_name, "n": int(msk.sum()),
                "obs_loss_pct": float(np.sum(ov * wt) / wsum),
                "sim_loss_pct": float(np.sum(sv * wt) / wsum),
                "bias_pctpt": float(np.sum(diff * wt) / wsum),
                "rmse_pct": float(np.sqrt(np.sum((diff ** 2) * wt) / wsum)),
                "mae_pct": float(np.sum(np.abs(diff) * wt) / wsum),
            })

    df_t3 = pd.DataFrame(t3_rows)
    df_t3.to_csv(os.path.join(OUT_DIR, "s13_s14_table_s3_pos_ma5.csv"), index=False)
    for sc in ["Subnational (CIL)", "National (FAOSTAT)", "Gridded (GDHY 0.5°)"]:
        print(f"  {sc} total N = {int(df_t3[df_t3['scale'] == sc]['n'].sum()):,}")
    print(f"  Negative bias count across 24 combinations: {int((df_t3['bias_pctpt'] < 0).sum())} / {len(df_t3)}")
    print(df_t3.round(2).to_string(index=False))

    # =========================================================================
    # 5. 18 DOCUMENTED LANDMARK DISASTERS IN 1990-2019 (30-YR RANK)
    # =========================================================================
    print("\n[5/6] Building 18 Documented Landmark Disasters (1990-2019, 30-yr ranks)...")
    iso_z = np.load(COUNTRY_GRID_CACHE, allow_pickle=True)
    iso_idx_2d = iso_z["iso_idx_2d"]
    iso_list = iso_z["iso_list"].tolist()
    iso_to_code = {k: i + 1 for i, k in enumerate(iso_list)}
    lat = cw_90_19["lat"].values
    lon = cw_90_19["lon"].values
    lon2d, lat2d = np.meshgrid(lon, lat)

    events_1990_2019 = [
        # 9 Documented Droughts (all >= 1990, covering maize, soybean, wheat, rice)
        ("D01", "2002 US maize", "drought", "USA", "maiz", 2002, 36.0, 46.5, -100.0, -80.0, "Rippey (2015)"),
        ("D02", "2002 US soybean", "drought", "USA", "soyb", 2002, 36.0, 46.5, -100.0, -80.0, "Boyer et al. (2013)"),
        ("D03", "2012 US maize", "drought", "USA", "maiz", 2012, 35.0, 45.5, -102.0, -82.0, "Rippey (2015)"),
        ("D04", "2012 US soybean", "drought", "USA", "soyb", 2012, 35.0, 45.5, -102.0, -82.0, "Rippey (2015)"),
        ("D05", "2003 France maize", "drought", "FRA", "maiz", 2003, 42.0, 51.2, -5.0, 8.5, "Ciais et al. (2005)"),
        ("D06", "2010 Russia wheat", "drought", "RUS", "whea", 2010, 48.0, 57.0, 38.0, 58.0, "Wegren (2011)"),
        ("D07", "2002 India rice", "drought", "IND", "rice", 2002, 15.0, 31.0, 73.0, 88.5, "Gadgil et al. (2003)"),
        ("D08", "2006 Australia wheat", "drought", "AUS", "whea", 2006, -38.5, -26.0, 135.0, 151.5, "van Dijk et al. (2013)"),
        ("D09", "2009 Argentina soybean", "drought", "ARG", "soyb", 2009, -38.5, -29.0, -65.5, -57.5, "Magrin et al. (2014)"),
        # 9 Documented Floods / Waterlogging (all >= 1990, covering maize, soybean, wheat, rice)
        ("W01", "1993 US maize", "waterlogging", "USA", "maiz", 1993, 38.5, 46.5, -97.5, -88.0, "Rosenzweig et al. (2002)"),
        ("W02", "1993 US soybean", "waterlogging", "USA", "soyb", 1993, 38.5, 46.5, -97.5, -88.0, "Li et al. (2019)"),
        ("W03", "2001 Argentina soybean", "waterlogging", "ARG", "soyb", 2001, -37.5, -33.5, -64.0, -59.5, "Kuppel et al. (2015)"),
        ("W04", "2002 Germany wheat", "waterlogging", "DEU", "whea", 2002, 50.0, 54.5, 9.5, 15.0, "EM-DAT / Elbe flood"),
        ("W05", "2016 France wheat", "waterlogging", "FRA", "whea", 2016, 46.5, 50.5, 0.0, 5.5, "Ben-Ari et al. (2018)"),
        ("W06", "2010 Pakistan rice", "waterlogging", "PAK", "rice", 2010, 24.0, 32.5, 67.0, 74.0, "GFD DFO_3691"),
        ("W07", "2011 Thailand rice", "waterlogging", "THA", "rice", 2011, 13.5, 18.0, 99.0, 101.8, "World Bank (2012)"),
        ("W08", "2007 Bangladesh rice", "waterlogging", "BGD", "rice", 2007, 21.5, 26.7, 88.0, 92.5, "GFD DFO_3136"),
        ("W09", "2016 China rice", "waterlogging", "CHN", "rice", 2016, 28.0, 33.0, 111.0, 121.5, "Yuan et al. (2018)"),
    ]

    gdhy_min3_cache = {cr: ma5_loss_3d_min3(gdhy_raw_cache[cr]) for cr in CROPS}
    lm_rows = []
    for eid, slbl, haz, iso, cr, yr, la0, la1, lo0, lo1, ref in events_1990_2019:
        gloss_3d = gdhy_min3_cache[cr]
        iy_cw = int(np.where(YEARS_90_19 == yr)[0][0])
        iy_gd = int(np.where(YEARS_90_16 == yr)[0][0])

        c_sub = cw_90_19.sel(crop=cr)
        area_tot = c_sub["area_ha"].values.astype(np.float64)
        irr_sh = c_sub["irr_share"].values
        dry_3d = c_sub["loss_dry_pct"].values
        wet_3d = c_sub["loss_wet_pct"].values
        ddays_3d = c_sub["dry_days"].values
        wdays_3d = c_sub["wet_days"].values

        iso_code = iso_to_code.get(iso, -999)
        reg_mask = (
            (iso_idx_2d == iso_code)
            & (lat2d >= la0) & (lat2d <= la1)
            & (lon2d >= lo0) & (lon2d <= lo1)
            & np.isfinite(area_tot) & (area_tot > 0.0)
        )
        m_box = reg_mask & np.isfinite(gloss_3d[iy_gd]) & np.isfinite(dry_3d[iy_cw]) & np.isfinite(wet_3d[iy_cw])
        if haz == "drought":
            m_aff = m_box & (ddays_3d[iy_cw] >= 7)
            if m_aff.sum() < 3:
                m_aff = m_box & (dry_3d[iy_cw] > 0)
            sim_aff = float(np.average(dry_3d[iy_cw][m_aff], weights=area_tot[m_aff]))
        else:
            m_aff = m_box & (wdays_3d[iy_cw] >= 7)
            if m_aff.sum() < 3:
                m_aff = m_box & (wet_3d[iy_cw] > 0)
            sim_aff = float(np.average(wet_3d[iy_cw][m_aff], weights=area_tot[m_aff]))
        obs_aff = float(np.average(gloss_3d[iy_gd][m_aff], weights=area_tot[m_aff]))

        w_reg = area_tot[reg_mask]
        obs_yr = gloss_3d[iy_gd][reg_mask]
        dry_yr = dry_3d[iy_cw][reg_mask]
        wet_yr = wet_3d[iy_cw][reg_mask]
        irr_yr = irr_sh[reg_mask]

        m_obs = np.isfinite(obs_yr) & (w_reg > 0)
        m_dry = np.isfinite(dry_yr) & (w_reg > 0)
        m_wet = np.isfinite(wet_yr) & (w_reg > 0)
        obs_loss = float(np.average(obs_yr[m_obs], weights=w_reg[m_obs])) if m_obs.any() else np.nan
        sim_dry = float(np.average(dry_yr[m_dry], weights=w_reg[m_dry])) if m_dry.any() else np.nan
        sim_wet = float(np.average(wet_yr[m_wet], weights=w_reg[m_wet])) if m_wet.any() else np.nan
        sim_haz = sim_dry if haz == "drought" else sim_wet
        irr_mean = float(np.average(np.where(np.isfinite(irr_yr), irr_yr, 0.0), weights=w_reg)) * 100.0

        ts_d, ts_w = [], []
        for t in range(len(YEARS_90_19)):
            vd = dry_3d[t][reg_mask]
            vw = wet_3d[t][reg_mask]
            md = np.isfinite(vd) & (w_reg > 0)
            mw = np.isfinite(vw) & (w_reg > 0)
            ts_d.append(float(np.average(vd[md], weights=w_reg[md])) if md.any() else np.nan)
            ts_w.append(float(np.average(vw[mw], weights=w_reg[mw])) if mw.any() else np.nan)
        arr_d = np.array(ts_d)
        arr_w = np.array(ts_w)
        rk_d = float(np.mean(arr_d[np.isfinite(arr_d)] <= sim_dry) * 100.0)
        rk_w = float(np.mean(arr_w[np.isfinite(arr_w)] <= sim_wet) * 100.0)

        lm_rows.append({
            "event_id": eid,
            "short_label": slbl,
            "hazard": haz,
            "iso": iso,
            "crop": cr,
            "year": yr,
            "n_cells_full_region": int(reg_mask.sum()),
            "n_cells_affected": int(m_aff.sum()),
            "crop_area_mha": float(w_reg.sum() / 1e6),
            "irr_area_share_pct": irr_mean,
            "obs_aff_ma5_loss_pct": obs_aff,
            "cwatm_aff_haz_loss_pct": sim_aff,
            "bias_aff_pctpt": float(sim_aff - obs_aff),
            "obs_reg_ma5_loss_pct": obs_loss,
            "cwatm_reg_haz_loss_pct": sim_haz,
            "bias_reg_pctpt": float(sim_haz - obs_loss),
            "cwatm_dry_loss_pct": sim_dry,
            "cwatm_wet_loss_pct": sim_wet,
            "cwatm_dry_rank_pct": rk_d,
            "cwatm_wet_rank_pct": rk_w,
            "reference": ref,
        })
    df_lm = pd.DataFrame(lm_rows)
    df_lm.to_csv(os.path.join(OUT_DIR, "s15_landmark_18_full_region.csv"), index=False)
    sub_d_lm = df_lm[df_lm["hazard"] == "drought"]
    sub_w_lm = df_lm[df_lm["hazard"] == "waterlogging"]
    print(f"  18 Landmark Disasters (1990-2019): drought mean bias = {sub_d_lm['bias_aff_pctpt'].mean():+.2f}%, "
          f"all-9 flood mean bias = {sub_w_lm['bias_aff_pctpt'].mean():+.2f}%, "
          f"single-season W01-W06 flood mean bias = {sub_w_lm.iloc[:6]['bias_aff_pctpt'].mean():+.2f}%")
    print(df_lm[["event_id", "short_label", "irr_area_share_pct", "obs_aff_ma5_loss_pct", "cwatm_aff_haz_loss_pct", "bias_aff_pctpt", "cwatm_dry_rank_pct", "cwatm_wet_rank_pct"]].round(2).to_string(index=False))

    # =========================================================================
    # 6. GLOBAL GDIS (1990-2016) & MODIS GFD (2000-2016, BASELINE 1990-2016)
    # =========================================================================
    print("\n[6/6] Building Global GDIS (1990-2016) & MODIS GFD (2000-2016) Catalogs...")
    gdis = pd.read_csv(GDIS_CSV, usecols=["iso3", "year", "disastertype", "latitude", "longitude"], low_memory=False)
    gdis["disastertype"] = gdis["disastertype"].astype(str).str.strip().str.lower()
    gdis = gdis[
        gdis["disastertype"].isin(["drought", "flood"])
        & (gdis["year"] >= 1990) & (gdis["year"] <= 2016)
        & np.isfinite(gdis["latitude"]) & np.isfinite(gdis["longitude"])
    ].copy()

    nt = len(YEARS_90_16)
    gdis_dry_3d = np.zeros((nt, len(lat), len(lon)), dtype=bool)
    gdis_wet_3d = np.zeros((nt, len(lat), len(lon)), dtype=bool)
    gdis["r0"] = np.clip(np.floor((90.0 - gdis["latitude"].values) / 0.5).astype(int), 0, 359)
    gdis["c0"] = np.clip(np.floor((gdis["longitude"].values + 180.0) / 0.5).astype(int), 0, 719)
    gdis["ti"] = gdis["year"].values.astype(int) - 1990

    for dtype, target_arr in [("drought", gdis_dry_3d), ("flood", gdis_wet_3d)]:
        sub_g = gdis[gdis["disastertype"] == dtype][["ti", "r0", "c0", "iso3"]].drop_duplicates()
        for ti, r0, c0, iso3 in sub_g.itertuples(index=False):
            k = iso_to_code.get(str(iso3).strip(), 0)
            r_lo, r_hi = max(0, r0 - 2), min(360, r0 + 3)
            c_lo, c_hi = max(0, c0 - 2), min(720, c0 + 3)
            if k > 0:
                target_arr[ti, r_lo:r_hi, c_lo:c_hi] |= (iso_idx_2d[r_lo:r_hi, c_lo:c_hi] == k)
            else:
                target_arr[ti, r_lo:r_hi, c_lo:c_hi] = True

    gfd = pd.read_csv(GFD_CSV, low_memory=False)
    gfd = gfd[(gfd["start_year"] >= 2000) & (gfd["start_year"] <= 2016) & (gfd["exposed_area_km2_total"] >= 1000.0)].copy()
    yrs_gfd = np.arange(2000, 2017, dtype=int)
    gfd_mask_3d = np.zeros((len(yrs_gfd), len(lat), len(lon)), dtype=bool)
    lat2d_col = lat[:, None]
    lon2d_row = lon[None, :]
    for row in gfd.itertuples(index=False):
        ti = int(row.start_year) - 2000
        w_b, s_b, e_b, n_b = float(row.raster_bounds_west), float(row.raster_bounds_south), float(row.raster_bounds_east), float(row.raster_bounds_north)
        clon, clat = float(row.centroid_lon), float(row.centroid_lat)
        box = (
            (lat2d_col >= max(s_b, clat - 3.5)) & (lat2d_col <= min(n_b, clat + 3.5))
            & (lon2d_row >= max(w_b, clon - 3.5)) & (lon2d_row <= min(e_b, clon + 3.5))
        )
        isos_ev = [x.strip() for x in str(row.country_code).split(",") if len(x.strip()) == 3]
        ks = [iso_to_code[c] for c in isos_ev if c in iso_to_code]
        if ks:
            box = box & np.isin(iso_idx_2d, ks)
        gfd_mask_3d[ti] |= box

    gdis_rows, gfd_rows = [], []
    for crop in CROPS:
        c_sub = cw_90_16.sel(crop=crop)
        area2d = c_sub["area_ha"].values.astype(np.float64)
        active = (area2d >= 500.0) & (iso_idx_2d > 0)
        w3 = np.broadcast_to(area2d[None, :, :], (nt, len(lat), len(lon)))

        obs_anom = detrend_rel_vec(gdhy_raw_cache[crop], YEARS_90_16, min_yrs=12)
        sim_2s = detrend_rel_vec(c_sub["yield_2s"].values, YEARS_90_16, min_yrs=12)
        sim_nw = detrend_rel_vec(c_sub["yield_nowl"].values, YEARS_90_16, min_yrs=12)
        ldry = c_sub["loss_dry_pct"].values
        lwet = c_sub["loss_wet_pct"].values
        ld_mean, ld_sd = np.nanmean(ldry, axis=0, keepdims=True), np.nanstd(ldry, axis=0, keepdims=True)
        lw_mean, lw_sd = np.nanmean(lwet, axis=0, keepdims=True), np.nanstd(lwet, axis=0, keepdims=True)
        z_ld = np.divide(ldry - ld_mean, ld_sd, out=np.zeros_like(ldry), where=ld_sd > 1e-5)
        z_lw = np.divide(lwet - lw_mean, lw_sd, out=np.zeros_like(lwet), where=lw_sd > 1e-5)

        valid = active[None, :, :] & np.isfinite(obs_anom) & np.isfinite(sim_2s) & np.isfinite(sim_nw)
        m_none = valid & (~gdis_dry_3d) & (~gdis_wet_3d)
        m_dry = valid & gdis_dry_3d & (~gdis_wet_3d)
        m_flood = valid & gdis_wet_3d & (~gdis_dry_3d)

        for label, msk in [("No_GDIS_Disaster", m_none), ("GDIS_Drought_Only", m_dry), ("GDIS_Flood_Only", m_flood)]:
            wv = w3[msk]
            wn = wv / wv.sum()
            _, p_dry_vs_none = stats.mannwhitneyu(ldry[msk], ldry[m_none], alternative="greater") if label != "No_GDIS_Disaster" else (np.nan, np.nan)
            _, p_wet_vs_none = stats.mannwhitneyu(lwet[msk], lwet[m_none], alternative="greater") if label != "No_GDIS_Disaster" else (np.nan, np.nan)
            gdis_rows.append({
                "crop": crop, "gdis_regime": label, "n_grid_years": int(msk.sum()),
                "obs_gdhy_anom_pct": float(np.sum(wn * obs_anom[msk])),
                "cwatm_loss_dry_pct": float(np.sum(wn * ldry[msk])),
                "cwatm_loss_wet_pct": float(np.sum(wn * lwet[msk])),
                "cwatm_z_dry": float(np.sum(wn * z_ld[msk])),
                "cwatm_z_wet": float(np.sum(wn * z_lw[msk])),
                "sim_2s_anom_pct": float(np.sum(wn * sim_2s[msk])),
                "sim_nowl_anom_pct": float(np.sum(wn * sim_nw[msk])),
                "rmse_2s_pct": float(np.sqrt(np.sum(wn * ((sim_2s[msk] - obs_anom[msk]) ** 2)))),
                "rmse_nowl_pct": float(np.sqrt(np.sum(wn * ((sim_nw[msk] - obs_anom[msk]) ** 2)))),
                "mannwhitney_p_dry_enrichment": p_dry_vs_none,
                "mannwhitney_p_wet_enrichment": p_wet_vs_none,
            })

        # GFD (2000-2016 slice: indices 10..26 of 1990..2016)
        sl = slice(10, 27)
        v_2000 = valid[sl]
        ever_gfd_cell = np.any(gfd_mask_3d, axis=0)[None, :, :]
        m_gfd_flood = v_2000 & gfd_mask_3d
        m_gfd_base = v_2000 & ever_gfd_cell & (~gfd_mask_3d)

        for label, msk in [("SameCells_NonFlood_Years_2000_2016", m_gfd_base), ("MODIS_GFD_Flooded_Years_2000_2016", m_gfd_flood)]:
            wv = w3[sl][msk]
            wn = wv / wv.sum()
            _, p_wet_gfd = stats.mannwhitneyu(lwet[sl][msk], lwet[sl][m_gfd_base], alternative="greater") if "MODIS" in label else (np.nan, np.nan)
            gfd_rows.append({
                "crop": crop, "gfd_regime": label, "n_grid_years": int(msk.sum()),
                "obs_gdhy_anom_pct": float(np.sum(wn * obs_anom[sl][msk])),
                "cwatm_loss_wet_pct": float(np.sum(wn * lwet[sl][msk])),
                "cwatm_loss_dry_pct": float(np.sum(wn * ldry[sl][msk])),
                "cwatm_z_wet": float(np.sum(wn * z_lw[sl][msk])),
                "sim_2s_anom_pct": float(np.sum(wn * sim_2s[sl][msk])),
                "sim_nowl_anom_pct": float(np.sum(wn * sim_nw[sl][msk])),
                "rmse_2s_pct": float(np.sqrt(np.sum(wn * ((sim_2s[sl][msk] - obs_anom[sl][msk]) ** 2)))),
                "rmse_nowl_pct": float(np.sqrt(np.sum(wn * ((sim_nw[sl][msk] - obs_anom[sl][msk]) ** 2)))),
                "mannwhitney_p_wet_enrichment": p_wet_gfd,
            })

    df_gdis = pd.DataFrame(gdis_rows)
    df_gfd = pd.DataFrame(gfd_rows)
    df_gdis.to_csv(os.path.join(OUT_DIR, "tier4b_gdis_global_disaster_contrast.csv"), index=False)
    df_gfd.to_csv(os.path.join(OUT_DIR, "tier4b_modis_gfd_flood_contrast.csv"), index=False)
    cw_all.close()

    tot_gdis = int(df_gdis[df_gdis["gdis_regime"] != "No_GDIS_Disaster"]["n_grid_years"].sum())
    tot_gfd_fld = int(df_gfd[df_gfd["gfd_regime"] == "MODIS_GFD_Flooded_Years_2000_2016"]["n_grid_years"].sum())
    tot_gfd_non = int(df_gfd[df_gfd["gfd_regime"] == "SameCells_NonFlood_Years_2000_2016"]["n_grid_years"].sum())
    print(f"  GDIS 1990-2016 disaster grid-years = {tot_gdis:,}")
    print(df_gdis[["crop", "gdis_regime", "n_grid_years", "cwatm_z_dry", "cwatm_z_wet", "mannwhitney_p_dry_enrichment", "mannwhitney_p_wet_enrichment"]].round(4).to_string(index=False))
    print(f"  MODIS GFD 2000-2016 flooded grid-years = {tot_gfd_fld:,}, same-cell non-flood = {tot_gfd_non:,}")
    print(df_gfd[["crop", "gfd_regime", "n_grid_years", "cwatm_loss_wet_pct", "mannwhitney_p_wet_enrichment"]].round(4).to_string(index=False))

    print(f"\nCompleted ALL 1990-2019 validation harmonizations in {time.time() - t0:.1f}s!")


if __name__ == "__main__":
    main()
