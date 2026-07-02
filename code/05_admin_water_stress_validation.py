# Public-release script split from the original analysis notebook.
# Paths remain project-specific and may need to be edited before rerunning.
# Part 5: administrative-unit validation of simulated water-stress losses.

from __future__ import annotations

import json
from pathlib import Path

import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.colors import TwoSlopeNorm

DEFAULT_ADMIN_SHP = Path(
    r"D:\AAUDE\cil-ag-replication-package\Fig1\Crop_Coverage\data\shapes\all_countries.shp"
)
DEFAULT_OBS_DIR = Path(r"F:\crop_filed\impact_data\for_regressions_full_csv")
DEFAULT_OUT_DIR = Path(
    r"D:\AAUDE\paper\paper4\revision_v2\fig\manu_fig\validation_hist_dry_wet_scatter"
)
DEFAULT_YEAR_MIN = 1900
DEFAULT_YEAR_MAX = 2019

CROPS = ("maiz", "soyb", "rice", "whea")
HAZARDS = ("dry", "wet")
CROP_LABEL = {"maiz": "Maize", "soyb": "Soybean", "rice": "Rice", "whea": "Wheat"}
HAZARD_LABEL = {"dry": "Drought", "wet": "Waterlogging"}
OBS_FILE = {
    "maiz": "corn_gmfd_v1_ready.csv",
    "soyb": "soy_gmfd_v1_ready.csv",
    "rice": "rice_gmfd_v1_ready.csv",
    "whea": "wheat_gmfd_v1_ready.csv",
}

MACRO_REGION_ORDER = ("North America", "South America", "Europe", "Africa", "Asia")
MACRO_REGION_ISO3 = {
    "North America": {
        "CAN", "USA", "MEX", "GTM", "BLZ", "HND", "SLV", "NIC", "CRI", "PAN",
        "CUB", "HTI", "DOM", "JAM", "TTO", "BRB", "BHS",
    },
    "South America": {
        "ARG", "BOL", "BRA", "CHL", "COL", "ECU", "GUY", "PRY", "PER", "SUR", "URY", "VEN",
    },
    "Europe": {
        "ALB", "AND", "AUT", "BEL", "BGR", "BIH", "BLR", "CHE", "CYP", "CZE", "DEU", "DNK",
        "ESP", "EST", "FIN", "FRA", "GBR", "GRC", "HRV", "HUN", "IRL", "ISL", "ITA", "LIE",
        "LTU", "LUX", "LVA", "MDA", "MKD", "MLT", "MNE", "NLD", "NOR", "POL", "PRT", "ROU",
        "SRB", "SVK", "SVN", "SWE", "UKR",
    },
    "Africa": {
        "AGO", "BDI", "BEN", "BFA", "BWA", "CAF", "CIV", "CMR", "COD", "COG", "DZA", "EGY",
        "ETH", "GAB", "GHA", "GIN", "GMB", "KEN", "LBR", "LBY", "LSO", "MAR", "MDG", "MLI",
        "MOZ", "MRT", "MWI", "NAM", "NER", "NGA", "RWA", "SDN", "SEN", "SLE", "SOM", "SSD",
        "SWZ", "TCD", "TGO", "TUN", "TZA", "UGA", "ZAF", "ZMB", "ZWE",
    },
    "Asia": {
        "AFG", "ARM", "AZE", "BGD", "BHR", "BRN", "BTN", "CHN", "GEO", "IDN", "IND", "IRN",
        "IRQ", "ISR", "JOR", "JPN", "KAZ", "KHM", "KOR", "KWT", "LAO", "LBN", "LKA", "MMR",
        "MNG", "MYS", "NPL", "OMN", "PAK", "PHL", "PRK", "QAT", "SAU", "SGP", "SYR", "THA",
        "TJK", "TKM", "TUR", "UZB", "VNM", "YEM",
    },
}
ISO3_TO_MACRO_REGION = {
    iso: region for region, isos in MACRO_REGION_ISO3.items() for iso in isos
}
ADMIN_MISSING = "__NA__"


def _normalize_admin_id_series(s: pd.Series, allow_missing: bool = False) -> pd.Series:
    out = s.astype("string").str.strip()
    out = out.str.replace(r"\.0$", "", regex=True)
    missing = out.isna() | out.str.lower().isin(["", "nan", "none", "<na>", "null"])
    if allow_missing:
        out = out.mask(missing, ADMIN_MISSING)
    else:
        out = out.mask(missing, pd.NA)
    return out


def _normalize_admin_id_value(value, allow_missing: bool = False) -> str | pd.NA:
    return _normalize_admin_id_series(pd.Series([value]), allow_missing=allow_missing).iloc[0]


def _normalize_admin_gdf(gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    gdf = gdf.copy()
    gdf["iso"] = gdf["iso"].astype(str).str.upper().str.strip()
    gdf["adm1_id"] = _normalize_admin_id_series(gdf["adm1_id"], allow_missing=False)
    gdf["adm2_id"] = _normalize_admin_id_series(gdf["adm2_id"], allow_missing=True)
    gdf = gdf[gdf["iso"].notna() & gdf["adm1_id"].notna()].copy()
    gdf = gdf.reset_index(drop=True)
    gdf["rid"] = np.arange(len(gdf), dtype=np.int32)
    return gdf


def _weighted_nanmean(values: np.ndarray, weights: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    weights = np.asarray(weights, dtype=float)
    m = np.isfinite(values) & np.isfinite(weights) & (weights > 0)
    if not np.any(m):
        return np.nan
    return float(np.average(values[m], weights=weights[m]))


def _centered_5yr_base(s: pd.Series) -> pd.Series:
    return s.rolling(window=5, center=True, min_periods=5).mean()


def _to_yearly(da: xr.DataArray, year_min: int, year_max: int) -> xr.DataArray:
    if "year" in da.dims:
        return da.sel(year=slice(year_min, year_max))
    if "time" in da.dims:
        out = da.sel(time=slice(f"{year_min}-01-01", f"{year_max}-12-31")).groupby("time.year").mean("time")
        return out.rename({"year": "year"}) if "year" in out.dims else out
    raise ValueError(f"DataArray has neither year/time dims: {da.dims}")


def get_db_block(DB: dict, scen: str = "hist", gcm: str | None = None) -> dict:
    """
    hist -> DB['hist'] only (no GCM).
    non-hist -> DB['future'][gcm][scen].
    """
    scen = str(scen).strip()
    if scen == "hist":
        if "hist" not in DB:
            raise KeyError("DB['hist'] not found.")
        return DB["hist"]
    if "future" not in DB:
        raise KeyError("DB['future'] not found.")
    if not gcm:
        raise ValueError("For non-hist scenario, gcm is required.")
    if gcm not in DB["future"] or scen not in DB["future"][gcm]:
        raise KeyError(f"Missing DB['future']['{gcm}']['{scen}'].")
    return DB["future"][gcm][scen]


def build_grid_admin_lookup(admin_shp: Path, template_da: xr.DataArray, cache_csv: Path | None = None) -> pd.DataFrame:
    if cache_csv is not None and cache_csv.exists():
        lookup = pd.read_csv(cache_csv)
        lookup["lat_key"] = np.round(pd.to_numeric(lookup["lat"], errors="coerce").astype(float), 8)
        lookup["lon_key"] = np.round(pd.to_numeric(lookup["lon"], errors="coerce").astype(float), 8)
        lookup["adm1_id"] = _normalize_admin_id_series(lookup["adm1_id"], allow_missing=False)
        lookup["adm2_id"] = _normalize_admin_id_series(lookup["adm2_id"], allow_missing=True)
        lookup = lookup.sort_values(["lat_key", "lon_key", "rid"], na_position="last")
        lookup = lookup.drop_duplicates(["lat_key", "lon_key"], keep="first")
        return lookup.drop(columns=["lat_key", "lon_key"])
    gdf = _normalize_admin_gdf(gpd.read_file(admin_shp).to_crs("EPSG:4326"))
    lat = template_da["lat"].values
    lon = template_da["lon"].values
    lon2d, lat2d = np.meshgrid(lon, lat)
    pts = gpd.GeoDataFrame(
        {"lat": lat2d.ravel(), "lon": lon2d.ravel()},
        geometry=gpd.points_from_xy(lon2d.ravel(), lat2d.ravel()),
        crs="EPSG:4326",
    )
    joined = gpd.sjoin(
        pts,
        gdf[["rid", "iso", "adm1_id", "adm2_id", "geometry"]],
        how="left",
        predicate="within",
    )
    lookup = joined[["lat", "lon", "rid", "iso", "adm1_id", "adm2_id"]].copy()
    lookup["rid"] = pd.to_numeric(lookup["rid"], errors="coerce").astype("Int64")
    lookup["lat_key"] = np.round(pd.to_numeric(lookup["lat"], errors="coerce").astype(float), 8)
    lookup["lon_key"] = np.round(pd.to_numeric(lookup["lon"], errors="coerce").astype(float), 8)
    lookup = lookup.sort_values(["lat_key", "lon_key", "rid"], na_position="last")
    lookup = lookup.drop_duplicates(["lat_key", "lon_key"], keep="first")
    lookup = lookup.drop(columns=["lat_key", "lon_key"])
    if cache_csv is not None:
        cache_csv.parent.mkdir(parents=True, exist_ok=True)
        lookup.to_csv(cache_csv, index=False, encoding="utf-8-sig")
    return lookup


def select_region_polygons_from_shp(
    admin_shp: Path,
    iso: str | None = None,
    adm1_id: int | None = None,
    adm2_id: int | None = None,
    fid: int | None = None,
) -> gpd.GeoDataFrame:
    gdf = _normalize_admin_gdf(gpd.read_file(admin_shp).to_crs("EPSG:4326"))
    sel = gdf.copy()
    if fid is not None:
        fid_col = next((c for c in ("FID", "fid", "OBJECTID", "objectid") if c in sel.columns), None)
        if fid_col is None:
            raise KeyError("FID-like field not found in shapefile.")
        sel = sel[pd.to_numeric(sel[fid_col], errors="coerce") == int(fid)]
    if iso is not None:
        sel = sel[sel["iso"] == str(iso).upper().strip()]
    if adm1_id is not None:
        sel = sel[sel["adm1_id"] == _normalize_admin_id_value(adm1_id, allow_missing=False)]
    if adm2_id is not None:
        sel = sel[sel["adm2_id"] == _normalize_admin_id_value(adm2_id, allow_missing=True)]
    if sel.empty:
        raise ValueError(f"No polygon matched: iso={iso}, adm1={adm1_id}, adm2={adm2_id}, fid={fid}")
    return sel


def build_grid_region_lookup(region_gdf: gpd.GeoDataFrame, template_da: xr.DataArray, cache_csv: Path | None = None) -> pd.DataFrame:
    if cache_csv is not None and cache_csv.exists():
        return pd.read_csv(cache_csv)
    lat = template_da["lat"].values
    lon = template_da["lon"].values
    lon2d, lat2d = np.meshgrid(lon, lat)
    pts = gpd.GeoDataFrame(
        {"lat": lat2d.ravel(), "lon": lon2d.ravel()},
        geometry=gpd.points_from_xy(lon2d.ravel(), lat2d.ravel()),
        crs="EPSG:4326",
    )
    reg = region_gdf.to_crs("EPSG:4326")[["geometry"]].copy()
    reg["rid_region"] = np.arange(len(reg), dtype=np.int32)
    joined = gpd.sjoin(pts, reg, how="inner", predicate="within")
    lookup = joined[["lat", "lon", "rid_region"]].drop_duplicates().copy()
    if cache_csv is not None:
        cache_csv.parent.mkdir(parents=True, exist_ok=True)
        lookup.to_csv(cache_csv, index=False, encoding="utf-8-sig")
    return lookup


def _model_fields(DB: dict, blk: dict, crop: str, year_min: int, year_max: int, yield_unit_factor: float, relative_yield_scale: float) -> dict:
    dist = DB["dist"]
    area_rf = dist[f"{crop}_nonIrr"].fillna(0.0)
    area_ir = dist[f"{crop}_Irr"].fillna(0.0)
    area_tot = area_rf + area_ir
    yp_rf = blk["potential_yield"][crop]["nonIrr"] * yield_unit_factor
    yp_ir = blk["potential_yield"][crop]["Irr"] * yield_unit_factor
    ws_rf = _to_yearly(blk["relative_yield"][crop]["waterstress_nonIrr"], year_min, year_max) / relative_yield_scale
    ws_ir = _to_yearly(blk["relative_yield"][crop]["waterstress_Irr"], year_min, year_max) / relative_yield_scale
    nws_rf = _to_yearly(blk["relative_yield"][crop]["nowaterstress_nonIrr"], year_min, year_max) / relative_yield_scale
    nws_ir = _to_yearly(blk["relative_yield"][crop]["nowaterstress_Irr"], year_min, year_max) / relative_yield_scale
    dry_rf = _to_yearly(blk["stress_days"][crop]["dry_nonIrr"], year_min, year_max).fillna(0.0)
    dry_ir = _to_yearly(blk["stress_days"][crop]["dry_Irr"], year_min, year_max).fillna(0.0)
    wet_rf = _to_yearly(blk["stress_days"][crop]["wet_nonIrr"], year_min, year_max).fillna(0.0)
    wet_ir = _to_yearly(blk["stress_days"][crop]["wet_Irr"], year_min, year_max).fillna(0.0)
    ws_rf_ok = (area_rf > 0) & np.isfinite(ws_rf) & np.isfinite(yp_rf)
    ws_ir_ok = (area_ir > 0) & np.isfinite(ws_ir) & np.isfinite(yp_ir)
    nws_rf_ok = (area_rf > 0) & np.isfinite(nws_rf) & np.isfinite(yp_rf)
    nws_ir_ok = (area_ir > 0) & np.isfinite(nws_ir) & np.isfinite(yp_ir)
    ws_prod = xr.where(ws_rf_ok, ws_rf * yp_rf * area_rf, 0.0) + xr.where(ws_ir_ok, ws_ir * yp_ir * area_ir, 0.0)
    nws_prod = xr.where(nws_rf_ok, nws_rf * yp_rf * area_rf, 0.0) + xr.where(nws_ir_ok, nws_ir * yp_ir * area_ir, 0.0)
    ws_area = xr.where(ws_rf_ok, area_rf, 0.0) + xr.where(ws_ir_ok, area_ir, 0.0)
    nws_area = xr.where(nws_rf_ok, area_rf, 0.0) + xr.where(nws_ir_ok, area_ir, 0.0)
    ws_y = xr.where(ws_area > 0, ws_prod / ws_area, np.nan)
    nws_y = xr.where(nws_area > 0, nws_prod / nws_area, np.nan)
    return {
        "area_rf": area_rf, "area_ir": area_ir, "area_tot": area_tot, "yp_rf": yp_rf, "yp_ir": yp_ir,
        "ws_rf": ws_rf, "ws_ir": ws_ir, "nws_rf": nws_rf, "nws_ir": nws_ir,
        "dry_rf": dry_rf, "dry_ir": dry_ir, "wet_rf": wet_rf, "wet_ir": wet_ir,
        "ws_y": ws_y, "nws_y": nws_y,
    }


def load_observed_region_from_shp(crop: str, obs_dir: Path, region_gdf: gpd.GeoDataFrame, year_min: int, year_max: int) -> pd.DataFrame:
    fp = obs_dir / OBS_FILE[crop]
    df = pd.read_csv(fp, usecols=["iso", "adm1_id", "adm2_id", "year", "yield"], low_memory=False)
    df["iso"] = df["iso"].astype(str).str.upper().str.strip()
    df["adm1_id"] = _normalize_admin_id_series(df["adm1_id"], allow_missing=False)
    df["adm2_id"] = _normalize_admin_id_series(df["adm2_id"], allow_missing=True)
    df["year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")
    df["yield"] = pd.to_numeric(df["yield"], errors="coerce")
    df = df.dropna(subset=["iso", "adm1_id", "adm2_id", "year", "yield"]).copy()
    df = df[(df["year"] >= year_min) & (df["year"] <= year_max)]
    keys = region_gdf[["iso", "adm1_id", "adm2_id"]].dropna().drop_duplicates().copy()
    keys["iso"] = keys["iso"].astype(str).str.upper().str.strip()
    keys["adm1_id"] = _normalize_admin_id_series(keys["adm1_id"], allow_missing=False)
    keys["adm2_id"] = _normalize_admin_id_series(keys["adm2_id"], allow_missing=True)
    m = df.merge(keys, on=["iso", "adm1_id", "adm2_id"], how="inner")
    if m.empty:
        return pd.DataFrame(columns=["crop", "year", "obs_yield_kg_ha", "obs_loss_pct"])
    out = m.groupby("year", as_index=False)["yield"].mean().rename(columns={"yield": "obs_yield_kg_ha"}).sort_values("year")
    out["crop"] = crop
    out["obs_baseline_kg_ha"] = _centered_5yr_base(out["obs_yield_kg_ha"])
    out["obs_anomaly_pct"] = (out["obs_yield_kg_ha"] - out["obs_baseline_kg_ha"]) / out["obs_baseline_kg_ha"] * 100.0
    out["obs_loss_pct_raw"] = -out["obs_anomaly_pct"]
    out["obs_loss_pct"] = out["obs_loss_pct_raw"].clip(lower=0.0, upper=100.0)
    return out


def compute_region_model_loss_timeseries_by_shp(
    DB: dict,
    region_lookup_df: pd.DataFrame,
    crop: str,
    scen: str = "hist",
    gcm: str | None = None,
    year_min: int = DEFAULT_YEAR_MIN,
    year_max: int = DEFAULT_YEAR_MAX,
    yield_unit_factor: float = 1.0,
    relative_yield_scale: float = 1.0,
    dry_event_days_thr: float = 7.0,
    wet_event_days_thr: float = 7.0,
    event_area_share_thr: float = 0.2,
) -> pd.DataFrame:
    blk = get_db_block(DB, scen=scen, gcm=gcm)
    f = _model_fields(DB, blk, crop, year_min, year_max, yield_unit_factor, relative_yield_scale)
    area_rf3 = f["area_rf"].broadcast_like(f["ws_rf"])
    area_ir3 = f["area_ir"].broadcast_like(f["ws_ir"])
    area_tot3 = f["area_tot"].broadcast_like(f["ws_rf"])
    yp_rf3 = f["yp_rf"].broadcast_like(f["ws_rf"])
    yp_ir3 = f["yp_ir"].broadcast_like(f["ws_ir"])
    ds = xr.Dataset(
        {
            "ws_rf": f["ws_rf"], "ws_ir": f["ws_ir"], "nws_rf": f["nws_rf"], "nws_ir": f["nws_ir"],
            "dry_rf": f["dry_rf"], "dry_ir": f["dry_ir"], "wet_rf": f["wet_rf"], "wet_ir": f["wet_ir"],
            "ws_y": f["ws_y"], "nws_y": f["nws_y"],
            "area_rf": area_rf3, "area_ir": area_ir3, "area_tot": area_tot3, "yp_rf": yp_rf3, "yp_ir": yp_ir3,
        }
    )
    df = ds.to_dataframe().reset_index()
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df[df["area_tot"].notna() & (df["area_tot"] > 0)].copy()
    df = df.merge(region_lookup_df[["lat", "lon"]].drop_duplicates(), on=["lat", "lon"], how="inner")
    if df.empty:
        return pd.DataFrame()
    rows = []
    for year, g in df.groupby("year", dropna=False):
        w_tot = pd.to_numeric(g["area_tot"], errors="coerce").values
        ws_y = pd.to_numeric(g["ws_y"], errors="coerce").values
        nws_y = pd.to_numeric(g["nws_y"], errors="coerce").values
        sim_ws = _weighted_nanmean(ws_y, w_tot)
        sim_nws = _weighted_nanmean(nws_y, w_tot)
        area_rf = pd.to_numeric(g["area_rf"], errors="coerce").values
        area_ir = pd.to_numeric(g["area_ir"], errors="coerce").values
        yp_rf = pd.to_numeric(g["yp_rf"], errors="coerce").values
        yp_ir = pd.to_numeric(g["yp_ir"], errors="coerce").values
        ws_rf = pd.to_numeric(g["ws_rf"], errors="coerce").values
        ws_ir = pd.to_numeric(g["ws_ir"], errors="coerce").values
        nws_rf = pd.to_numeric(g["nws_rf"], errors="coerce").values
        nws_ir = pd.to_numeric(g["nws_ir"], errors="coerce").values
        dry_rf = pd.to_numeric(g["dry_rf"], errors="coerce").values >= float(dry_event_days_thr)
        dry_ir = pd.to_numeric(g["dry_ir"], errors="coerce").values >= float(dry_event_days_thr)
        wet_rf = pd.to_numeric(g["wet_rf"], errors="coerce").values >= float(wet_event_days_thr)
        wet_ir = pd.to_numeric(g["wet_ir"], errors="coerce").values >= float(wet_event_days_thr)
        dry_grid_event = (dry_rf & (area_rf > 0)) | (dry_ir & (area_ir > 0))
        wet_grid_event = (wet_rf & (area_rf > 0)) | (wet_ir & (area_ir > 0))
        num_d = np.nansum(((nws_rf - ws_rf) * area_rf * yp_rf) * dry_rf) + np.nansum(((nws_ir - ws_ir) * area_ir * yp_ir) * dry_ir)
        num_w = np.nansum(((nws_rf - ws_rf) * area_rf * yp_rf) * wet_rf) + np.nansum(((nws_ir - ws_ir) * area_ir * yp_ir) * wet_ir)
        den = np.nansum(ws_rf * area_rf * yp_rf) + np.nansum(ws_ir * area_ir * yp_ir)
        dry_loss = float(num_d / den * 100.0) if np.isfinite(den) and den > 0 else np.nan
        wet_loss = float(num_w / den * 100.0) if np.isfinite(den) and den > 0 else np.nan
        wsum = np.nansum(w_tot)
        dry_share = float((np.nansum(area_rf[dry_rf]) + np.nansum(area_ir[dry_ir])) / wsum) if wsum > 0 else np.nan
        wet_share = float((np.nansum(area_rf[wet_rf]) + np.nansum(area_ir[wet_ir])) / wsum) if wsum > 0 else np.nan
        rows.append(
            {
                "crop": crop,
                "year": int(year),
                "sim_yield_ws_kg_ha": sim_ws,
                "sim_yield_nws_kg_ha": sim_nws,
                "sim_direct_loss_pct": float((sim_nws - sim_ws) / sim_nws * 100.0) if np.isfinite(sim_nws) and sim_nws > 0 else np.nan,
                "model_loss_dry_pct": dry_loss,
                "model_loss_wet_pct": wet_loss,
                "dry_event_area_share": dry_share,
                "wet_event_area_share": wet_share,
                "dry_event_grid_count": int(np.nansum(dry_grid_event)),
                "wet_event_grid_count": int(np.nansum(wet_grid_event)),
                "is_dry_event_year": bool(np.any(dry_grid_event)),
                "is_wet_event_year": bool(np.any(wet_grid_event)),
            }
        )
    out = pd.DataFrame(rows).sort_values(["crop", "year"])
    if out.empty:
        return out
    out["sim_baseline_kg_ha"] = out.groupby("crop", dropna=False)["sim_yield_ws_kg_ha"].transform(_centered_5yr_base)
    out["sim_anomaly_pct"] = (out["sim_yield_ws_kg_ha"] - out["sim_baseline_kg_ha"]) / out["sim_baseline_kg_ha"] * 100.0
    out["sim_loss_pct_raw"] = -out["sim_anomaly_pct"]
    out["sim_loss_pct"] = out["sim_loss_pct_raw"].clip(lower=0.0, upper=100.0)
    return out


def build_region_hazard_comparison_points(model_region_df: pd.DataFrame, obs_region_df: pd.DataFrame) -> pd.DataFrame:
    df = model_region_df.merge(
        obs_region_df[["crop", "year", "obs_yield_kg_ha", "obs_anomaly_pct", "obs_loss_pct"]],
        on=["crop", "year"],
        how="inner",
    )
    if df.empty:
        return pd.DataFrame()
    rows = []
    for haz in HAZARDS:
        sub = df[df[f"is_{haz}_event_year"]].copy()
        if sub.empty:
            continue
        sub["hazard"] = haz
        sim_raw = sub.get("sim_loss_pct_raw", pd.Series(np.nan, index=sub.index))
        sub["sim_yield_anomaly_loss_pct_raw"] = pd.to_numeric(sim_raw, errors="coerce")
        sub["sim_yield_anomaly_loss_pct"] = pd.to_numeric(sub["sim_loss_pct"], errors="coerce")
        sub["sim_direct_stress_loss_pct"] = pd.to_numeric(sub[f"model_loss_{haz}_pct"], errors="coerce")
        sub["sim_loss_pct_raw"] = sub["sim_yield_anomaly_loss_pct_raw"]
        sub["sim_loss_pct"] = sub["sim_yield_anomaly_loss_pct"]
        sub["obs_loss_pct"] = pd.to_numeric(sub["obs_loss_pct"], errors="coerce")
        sub["event_area_share"] = pd.to_numeric(sub[f"{haz}_event_area_share"], errors="coerce")
        sub = sub.dropna(subset=["sim_loss_pct", "obs_loss_pct", "event_area_share"])
        rows.append(sub)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def _template_grid_index(template_da: xr.DataArray) -> pd.DataFrame:
    lat = template_da["lat"].values
    lon = template_da["lon"].values
    lon2d, lat2d = np.meshgrid(lon, lat)
    return pd.DataFrame(
        {
            "cell": np.arange(lat2d.size, dtype=np.int64),
            "lat": lat2d.ravel(),
            "lon": lon2d.ravel(),
            "lat_key": np.round(lat2d.ravel().astype(float), 8),
            "lon_key": np.round(lon2d.ravel().astype(float), 8),
        }
    )


def _lookup_with_cell_index(admin_lookup_df: pd.DataFrame, template_da: xr.DataArray) -> pd.DataFrame:
    grid = _template_grid_index(template_da)
    lk = admin_lookup_df[["lat", "lon", "rid", "iso", "adm1_id", "adm2_id"]].copy()
    lk["lat_key"] = np.round(pd.to_numeric(lk["lat"], errors="coerce").astype(float), 8)
    lk["lon_key"] = np.round(pd.to_numeric(lk["lon"], errors="coerce").astype(float), 8)
    lk["rid"] = pd.to_numeric(lk["rid"], errors="coerce").astype("Int64")
    lk["iso"] = lk["iso"].astype(str).str.upper().str.strip()
    lk["adm1_id"] = _normalize_admin_id_series(lk["adm1_id"], allow_missing=False)
    lk["adm2_id"] = _normalize_admin_id_series(lk["adm2_id"], allow_missing=True)
    lk = lk.sort_values(["lat_key", "lon_key", "rid"], na_position="last")
    lk = lk.drop_duplicates(["lat_key", "lon_key"], keep="first")
    out = grid.merge(
        lk[["lat_key", "lon_key", "rid", "iso", "adm1_id", "adm2_id"]],
        on=["lat_key", "lon_key"],
        how="left",
        validate="one_to_one",
    )
    return out.sort_values("cell").reset_index(drop=True)


def _flat_lat_lon(da: xr.DataArray) -> np.ndarray:
    return da.transpose("lat", "lon").values.reshape(-1)


def _year_lat_lon_matrix(da: xr.DataArray) -> tuple[np.ndarray, np.ndarray]:
    da = da.transpose("year", "lat", "lon")
    years = pd.to_numeric(pd.Index(da["year"].values), errors="coerce").astype(int).to_numpy()
    values = da.values.reshape(len(years), -1)
    return years, values


def _bincount_sum(rid: np.ndarray, values: np.ndarray, minlength: int) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    values = np.where(np.isfinite(values), values, 0.0)
    return np.bincount(rid, weights=values, minlength=minlength)


def _safe_ratio(num: np.ndarray, den: np.ndarray) -> np.ndarray:
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    out = np.full_like(num, np.nan, dtype=float)
    m = np.isfinite(num) & np.isfinite(den) & (den > 0)
    out[m] = num[m] / den[m]
    return out


def load_observed_admin_yield(
    crop: str,
    obs_dir: Path = DEFAULT_OBS_DIR,
    year_min: int = DEFAULT_YEAR_MIN,
    year_max: int = DEFAULT_YEAR_MAX,
) -> pd.DataFrame:
    """Read observed admin-unit yield and compute the 5-year sliding anomaly loss."""
    fp = Path(obs_dir) / OBS_FILE[crop]
    df = pd.read_csv(fp, usecols=["iso", "adm1_id", "adm2_id", "year", "yield"], low_memory=False)
    df["iso"] = df["iso"].astype(str).str.upper().str.strip()
    df["adm1_id"] = _normalize_admin_id_series(df["adm1_id"], allow_missing=False)
    df["adm2_id"] = _normalize_admin_id_series(df["adm2_id"], allow_missing=True)
    df["year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")
    df["yield"] = pd.to_numeric(df["yield"], errors="coerce")
    df = df.dropna(subset=["iso", "adm1_id", "adm2_id", "year", "yield"]).copy()
    df = df[(df["year"] >= int(year_min)) & (df["year"] <= int(year_max))]
    out = (
        df.groupby(["iso", "adm1_id", "adm2_id", "year"], as_index=False)["yield"]
        .mean()
        .rename(columns={"yield": "obs_yield_kg_ha"})
        .sort_values(["iso", "adm1_id", "adm2_id", "year"])
    )
    out["crop"] = crop
    keys = ["crop", "iso", "adm1_id", "adm2_id"]
    out["obs_baseline_kg_ha"] = out.groupby(keys, dropna=False)["obs_yield_kg_ha"].transform(_centered_5yr_base)
    out["obs_anomaly_pct"] = (out["obs_yield_kg_ha"] - out["obs_baseline_kg_ha"]) / out["obs_baseline_kg_ha"] * 100.0
    out["obs_loss_pct_raw"] = -out["obs_anomaly_pct"]
    out["obs_loss_pct"] = out["obs_loss_pct_raw"].clip(lower=0.0, upper=100.0)
    return out


def aggregate_model_crop_to_admin(
    DB: dict,
    admin_lookup_df: pd.DataFrame,
    crop: str,
    year_min: int = DEFAULT_YEAR_MIN,
    year_max: int = DEFAULT_YEAR_MAX,
    yield_unit_factor: float = 1.0,
    relative_yield_scale: float = 1.0,
    dry_event_days_thr: float = 7.0,
    wet_event_days_thr: float = 7.0,
) -> pd.DataFrame:
    """
    Aggregate gridded historical simulation to every shp admin unit.

    Event rule: if any crop-area grid cell in the admin unit reaches the stress-day
    threshold, that admin-year is treated as a water-stress event.
    """
    blk = get_db_block(DB, scen="hist", gcm=None)
    f = _model_fields(DB, blk, crop, year_min, year_max, yield_unit_factor, relative_yield_scale)
    template = DB["dist"][f"{crop}_nonIrr"]
    lookup = _lookup_with_cell_index(admin_lookup_df, template)

    area_rf_all = _flat_lat_lon(f["area_rf"])
    area_ir_all = _flat_lat_lon(f["area_ir"])
    area_tot_all = _flat_lat_lon(f["area_tot"])
    yp_rf_all = _flat_lat_lon(f["yp_rf"])
    yp_ir_all = _flat_lat_lon(f["yp_ir"])

    valid = lookup["rid"].notna().to_numpy() & np.isfinite(area_tot_all) & (area_tot_all > 0)
    if not np.any(valid):
        return pd.DataFrame()

    cell = lookup.loc[valid, "cell"].astype(np.int64).to_numpy()
    rid = lookup.loc[valid, "rid"].astype(np.int64).to_numpy()
    n_rids = int(max(rid.max(), pd.to_numeric(admin_lookup_df["rid"], errors="coerce").max())) + 1
    meta = (
        lookup.loc[valid, ["rid", "iso", "adm1_id", "adm2_id"]]
        .drop_duplicates("rid")
        .sort_values("rid")
        .reset_index(drop=True)
    )
    meta["rid"] = meta["rid"].astype(np.int64)
    meta_rid = meta["rid"].to_numpy()

    area_rf = area_rf_all[cell]
    area_ir = area_ir_all[cell]
    area_tot = area_tot_all[cell]
    yp_rf = yp_rf_all[cell]
    yp_ir = yp_ir_all[cell]
    area_sum = _bincount_sum(rid, area_tot, n_rids)

    years, ws_rf_all = _year_lat_lon_matrix(f["ws_rf"])
    _, ws_ir_all = _year_lat_lon_matrix(f["ws_ir"])
    _, nws_rf_all = _year_lat_lon_matrix(f["nws_rf"])
    _, nws_ir_all = _year_lat_lon_matrix(f["nws_ir"])
    _, dry_rf_all = _year_lat_lon_matrix(f["dry_rf"])
    _, dry_ir_all = _year_lat_lon_matrix(f["dry_ir"])
    _, wet_rf_all = _year_lat_lon_matrix(f["wet_rf"])
    _, wet_ir_all = _year_lat_lon_matrix(f["wet_ir"])

    rows = []
    for yi, year in enumerate(years):
        ws_rf = ws_rf_all[yi, cell]
        ws_ir = ws_ir_all[yi, cell]
        nws_rf = nws_rf_all[yi, cell]
        nws_ir = nws_ir_all[yi, cell]
        dry_rf_days = dry_rf_all[yi, cell]
        dry_ir_days = dry_ir_all[yi, cell]
        wet_rf_days = wet_rf_all[yi, cell]
        wet_ir_days = wet_ir_all[yi, cell]

        ws_rf_ok = np.isfinite(ws_rf) & np.isfinite(yp_rf) & np.isfinite(area_rf) & (area_rf > 0)
        ws_ir_ok = np.isfinite(ws_ir) & np.isfinite(yp_ir) & np.isfinite(area_ir) & (area_ir > 0)
        nws_rf_ok = np.isfinite(nws_rf) & np.isfinite(yp_rf) & np.isfinite(area_rf) & (area_rf > 0)
        nws_ir_ok = np.isfinite(nws_ir) & np.isfinite(yp_ir) & np.isfinite(area_ir) & (area_ir > 0)
        ws_prod = np.where(ws_rf_ok, ws_rf * yp_rf * area_rf, 0.0) + np.where(ws_ir_ok, ws_ir * yp_ir * area_ir, 0.0)
        nws_prod = np.where(nws_rf_ok, nws_rf * yp_rf * area_rf, 0.0) + np.where(nws_ir_ok, nws_ir * yp_ir * area_ir, 0.0)
        ws_area = np.where(ws_rf_ok, area_rf, 0.0) + np.where(ws_ir_ok, area_ir, 0.0)
        nws_area = np.where(nws_rf_ok, area_rf, 0.0) + np.where(nws_ir_ok, area_ir, 0.0)
        ws_y = np.where(ws_area > 0, ws_prod / ws_area, np.nan)
        nws_y = np.where(nws_area > 0, nws_prod / nws_area, np.nan)

        valid_ws = np.isfinite(ws_y) & np.isfinite(ws_area) & (ws_area > 0)
        valid_nws = np.isfinite(nws_y) & np.isfinite(nws_area) & (nws_area > 0)
        sim_ws = _safe_ratio(
            _bincount_sum(rid[valid_ws], ws_y[valid_ws] * ws_area[valid_ws], n_rids),
            _bincount_sum(rid[valid_ws], ws_area[valid_ws], n_rids),
        )
        sim_nws = _safe_ratio(
            _bincount_sum(rid[valid_nws], nws_y[valid_nws] * nws_area[valid_nws], n_rids),
            _bincount_sum(rid[valid_nws], nws_area[valid_nws], n_rids),
        )

        den = _bincount_sum(rid, ws_prod, n_rids)
        dry_rf_event = np.isfinite(dry_rf_days) & (dry_rf_days >= float(dry_event_days_thr)) & (area_rf > 0)
        dry_ir_event = np.isfinite(dry_ir_days) & (dry_ir_days >= float(dry_event_days_thr)) & (area_ir > 0)
        wet_rf_event = np.isfinite(wet_rf_days) & (wet_rf_days >= float(wet_event_days_thr)) & (area_rf > 0)
        wet_ir_event = np.isfinite(wet_ir_days) & (wet_ir_days >= float(wet_event_days_thr)) & (area_ir > 0)
        dry_grid_event = dry_rf_event | dry_ir_event
        wet_grid_event = wet_rf_event | wet_ir_event

        dry_num = np.zeros_like(area_tot, dtype=float)
        wet_num = np.zeros_like(area_tot, dtype=float)
        m = dry_rf_event & np.isfinite(nws_rf) & np.isfinite(ws_rf) & np.isfinite(yp_rf)
        dry_num[m] += (nws_rf[m] - ws_rf[m]) * yp_rf[m] * area_rf[m]
        m = dry_ir_event & np.isfinite(nws_ir) & np.isfinite(ws_ir) & np.isfinite(yp_ir)
        dry_num[m] += (nws_ir[m] - ws_ir[m]) * yp_ir[m] * area_ir[m]
        m = wet_rf_event & np.isfinite(nws_rf) & np.isfinite(ws_rf) & np.isfinite(yp_rf)
        wet_num[m] += (nws_rf[m] - ws_rf[m]) * yp_rf[m] * area_rf[m]
        m = wet_ir_event & np.isfinite(nws_ir) & np.isfinite(ws_ir) & np.isfinite(yp_ir)
        wet_num[m] += (nws_ir[m] - ws_ir[m]) * yp_ir[m] * area_ir[m]

        dry_area = _bincount_sum(rid, np.where(dry_rf_event, area_rf, 0.0) + np.where(dry_ir_event, area_ir, 0.0), n_rids)
        wet_area = _bincount_sum(rid, np.where(wet_rf_event, area_rf, 0.0) + np.where(wet_ir_event, area_ir, 0.0), n_rids)
        dry_count = np.bincount(rid[dry_grid_event], minlength=n_rids) if np.any(dry_grid_event) else np.zeros(n_rids)
        wet_count = np.bincount(rid[wet_grid_event], minlength=n_rids) if np.any(wet_grid_event) else np.zeros(n_rids)

        part = meta.copy()
        part["crop"] = crop
        part["year"] = int(year)
        part["crop_area_sum"] = area_sum[meta_rid]
        part["sim_yield_ws_kg_ha"] = sim_ws[meta_rid]
        part["sim_yield_nws_kg_ha"] = sim_nws[meta_rid]
        part["sim_direct_loss_pct"] = _safe_ratio(sim_nws - sim_ws, sim_nws)[meta_rid] * 100.0
        part["model_loss_dry_pct_raw"] = _safe_ratio(_bincount_sum(rid, dry_num, n_rids), den)[meta_rid] * 100.0
        part["model_loss_wet_pct_raw"] = _safe_ratio(_bincount_sum(rid, wet_num, n_rids), den)[meta_rid] * 100.0
        part["model_loss_dry_pct"] = part["model_loss_dry_pct_raw"].clip(lower=0.0, upper=100.0)
        part["model_loss_wet_pct"] = part["model_loss_wet_pct_raw"].clip(lower=0.0, upper=100.0)
        part["dry_event_area_share"] = _safe_ratio(dry_area, area_sum)[meta_rid]
        part["wet_event_area_share"] = _safe_ratio(wet_area, area_sum)[meta_rid]
        part["dry_event_grid_count"] = dry_count[meta_rid].astype(int)
        part["wet_event_grid_count"] = wet_count[meta_rid].astype(int)
        part["is_dry_event_year"] = part["dry_event_grid_count"] > 0
        part["is_wet_event_year"] = part["wet_event_grid_count"] > 0
        rows.append(part)

    out = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if out.empty:
        return out
    out = out.sort_values(["crop", "iso", "adm1_id", "adm2_id", "year"]).reset_index(drop=True)
    keys = ["crop", "iso", "adm1_id", "adm2_id"]
    out["sim_baseline_kg_ha"] = out.groupby(keys, dropna=False)["sim_yield_ws_kg_ha"].transform(_centered_5yr_base)
    out["sim_anomaly_pct"] = (out["sim_yield_ws_kg_ha"] - out["sim_baseline_kg_ha"]) / out["sim_baseline_kg_ha"] * 100.0
    out["sim_loss_pct_raw"] = -out["sim_anomaly_pct"]
    out["sim_loss_pct"] = out["sim_loss_pct_raw"].clip(lower=0.0, upper=100.0)
    return out


def build_admin_hazard_comparison_points(
    model_admin_df: pd.DataFrame,
    obs_admin_df: pd.DataFrame,
) -> pd.DataFrame:
    keys = ["crop", "iso", "adm1_id", "adm2_id", "year"]
    obs_cols = keys + ["obs_yield_kg_ha", "obs_baseline_kg_ha", "obs_anomaly_pct", "obs_loss_pct_raw", "obs_loss_pct"]
    df = model_admin_df.merge(obs_admin_df[obs_cols], on=keys, how="inner")
    if df.empty:
        return pd.DataFrame()
    rows = []
    for haz in HAZARDS:
        flag = f"is_{haz}_event_year"
        if flag not in df.columns:
            continue
        sub = df[df[flag]].copy()
        if sub.empty:
            continue
        sub["hazard"] = haz
        sub["sim_yield_anomaly_loss_pct_raw"] = pd.to_numeric(sub["sim_loss_pct_raw"], errors="coerce")
        sub["sim_yield_anomaly_loss_pct"] = pd.to_numeric(sub["sim_loss_pct"], errors="coerce")
        sub["sim_direct_stress_loss_pct_raw"] = pd.to_numeric(sub[f"model_loss_{haz}_pct_raw"], errors="coerce")
        sub["sim_direct_stress_loss_pct"] = pd.to_numeric(sub[f"model_loss_{haz}_pct"], errors="coerce")
        sub["sim_loss_pct_raw"] = sub["sim_yield_anomaly_loss_pct_raw"]
        sub["sim_loss_pct"] = sub["sim_yield_anomaly_loss_pct"]
        sub["event_area_share"] = pd.to_numeric(sub[f"{haz}_event_area_share"], errors="coerce")
        sub["event_grid_count"] = pd.to_numeric(sub[f"{haz}_event_grid_count"], errors="coerce")
        sub["obs_loss_pct"] = pd.to_numeric(sub["obs_loss_pct"], errors="coerce")
        sub = sub.dropna(subset=["sim_loss_pct", "obs_loss_pct", "event_area_share"])
        rows.append(sub)
    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    return out.sort_values(["hazard", "crop", "iso", "adm1_id", "adm2_id", "year"]).reset_index(drop=True)


def calc_metrics(df_points: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (haz, crop), g in df_points.groupby(["hazard", "crop"], dropna=False):
        x = g["obs_loss_pct"].values
        y = g["sim_loss_pct"].values
        r = float(np.corrcoef(x, y)[0, 1]) if len(g) > 1 and np.std(x) > 0 and np.std(y) > 0 else np.nan
        rows.append(
            {
                "hazard": haz,
                "crop": crop,
                "n_event_admin_years": int(len(g)),
                "pearson_r": r,
                "rmse_loss_pct": float(np.sqrt(np.mean((y - x) ** 2))),
                "mean_bias_pct_points": float(np.mean(y - x)),
            }
        )
    return pd.DataFrame(rows).sort_values(["hazard", "crop"])


def plot_scatter_hist_dry_wet(df_points: pd.DataFrame, metrics: pd.DataFrame, out_png: Path, out_pdf: Path) -> None:
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["font.serif"] = ["Arial"]
    mpl.rcParams["axes.unicode_minus"] = False
    fig, axes = plt.subplots(2, 4, figsize=(18, 8), dpi=240)
    cmap = plt.get_cmap("Spectral_r")
    for i, haz in enumerate(HAZARDS):
        for j, crop in enumerate(CROPS):
            ax = axes[i, j]
            g = df_points[(df_points["hazard"] == haz) & (df_points["crop"] == crop)].copy()
            if g.empty:
                ax.axis("off")
                continue
            c = g["event_area_share"].clip(0, 1)
            ax.scatter(g["obs_loss_pct"], g["sim_loss_pct"], c=c, cmap=cmap, s=18 + 36 * c, alpha=0.8, edgecolor="none")
            lo = float(min(g["obs_loss_pct"].min(), g["sim_loss_pct"].min()))
            hi = float(max(g["obs_loss_pct"].max(), g["sim_loss_pct"].max()))
            pad = 0.08 * (hi - lo + 1e-6)
            ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "--", color="#222222", lw=1.1)
            ax.set_xlim(lo - pad, hi + pad)
            ax.set_ylim(lo - pad, hi + pad)
            ax.grid(alpha=0.22)
            ax.set_title(f"{HAZARD_LABEL[haz]} | {CROP_LABEL[crop]}", fontsize=11)
            ax.set_xlabel("Observed 5-year anomaly loss (%)")
            ax.set_ylabel("Simulated 5-year anomaly loss (%)")
    cbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=1), cmap=cmap),
        ax=axes.ravel().tolist(),
        fraction=0.02,
        pad=0.01,
    )
    cbar.set_label("Event area share in admin unit", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=320)
    fig.savefig(out_pdf)
    plt.close(fig)


def plot_scatter_all_regions_all_stress_points(df_points: pd.DataFrame, out_png: Path, out_pdf: Path) -> None:
    if df_points.empty:
        return
    df = df_points.copy()
    df["hazard"] = df["hazard"].astype(str).str.lower().str.strip()
    df["crop"] = df["crop"].astype(str).str.strip()
    df["obs_loss_pct"] = pd.to_numeric(df["obs_loss_pct"], errors="coerce")
    df["sim_loss_pct"] = pd.to_numeric(df["sim_loss_pct"], errors="coerce")
    df = df.dropna(subset=["obs_loss_pct", "sim_loss_pct", "hazard", "crop"])
    if df.empty:
        return
    hazard_colors = {"dry": "#c44e52", "wet": "#4c72b0"}
    crop_markers = {"maiz": "o", "soyb": "s", "rice": "^", "whea": "D"}
    fig, ax = plt.subplots(figsize=(10, 8), dpi=260)
    for haz in HAZARDS:
        for crop in CROPS:
            g = df[(df["hazard"] == haz) & (df["crop"] == crop)]
            if g.empty:
                continue
            ax.scatter(
                g["obs_loss_pct"], g["sim_loss_pct"], s=18, marker=crop_markers[crop], c=hazard_colors[haz], alpha=0.42, edgecolor="none"
            )
    lo = float(min(df["obs_loss_pct"].min(), df["sim_loss_pct"].min()))
    hi = float(max(df["obs_loss_pct"].max(), df["sim_loss_pct"].max()))
    pad = 0.06 * (hi - lo + 1e-6)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], "--", color="#222222", lw=1.2)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad)
    ax.grid(alpha=0.22)
    ax.set_xlabel("Observed 5-year anomaly loss (%)")
    ax.set_ylabel("Simulated 5-year anomaly loss (%)")
    x = df["obs_loss_pct"].values
    y = df["sim_loss_pct"].values
    r = float(np.corrcoef(x, y)[0, 1]) if len(df) > 1 and np.std(x) > 0 and np.std(y) > 0 else np.nan
    rmse = float(np.sqrt(np.mean((y - x) ** 2)))
    bias = float(np.mean(y - x))
    ax.text(0.02, 0.98, f"N={len(df)}\nR={r:.2f}\nRMSE={rmse:.2f}%\nBias={bias:.2f} pp", transform=ax.transAxes, ha="left", va="top")
    ax.set_title("All Regions, All Water-Stress Points (Dry + Wet)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=320)
    fig.savefig(out_pdf)
    plt.close(fig)


def _loss_metrics(x: pd.Series, y: pd.Series) -> dict:
    x = pd.to_numeric(x, errors="coerce")
    y = pd.to_numeric(y, errors="coerce")
    m = x.notna() & y.notna()
    x = x[m].to_numpy(dtype=float)
    y = y[m].to_numpy(dtype=float)
    if len(x) == 0:
        return {"n": 0, "pearson_r": np.nan, "rmse_loss_pct": np.nan, "mean_bias_pct_points": np.nan}
    r = float(np.corrcoef(x, y)[0, 1]) if len(x) > 1 and np.std(x) > 0 and np.std(y) > 0 else np.nan
    return {
        "n": int(len(x)),
        "pearson_r": r,
        "rmse_loss_pct": float(np.sqrt(np.mean((y - x) ** 2))),
        "mean_bias_pct_points": float(np.mean(y - x)),
    }


def calc_metrics_for_ycol(df_points: pd.DataFrame, y_col: str = "sim_loss_pct") -> pd.DataFrame:
    rows = []
    df = df_points.copy()
    df["obs_loss_pct"] = pd.to_numeric(df["obs_loss_pct"], errors="coerce")
    df[y_col] = pd.to_numeric(df[y_col], errors="coerce")
    for (haz, crop), g in df.groupby(["hazard", "crop"], dropna=False):
        metrics = _loss_metrics(g["obs_loss_pct"], g[y_col])
        rows.append({"hazard": haz, "crop": crop, "sim_loss_column": y_col, **metrics})
    return pd.DataFrame(rows).sort_values(["hazard", "crop"]).reset_index(drop=True)


def plot_validation_hexbin_by_crop_hazard(
    df_points: pd.DataFrame,
    out_png: Path,
    out_pdf: Path,
    y_col: str = "sim_loss_pct",
    y_label: str = "Simulated 5-year anomaly loss (%)",
    title: str = "Historical validation",
    max_loss: float = 100.0,
    gridsize: int = 38,
) -> pd.DataFrame:
    """
    Paper-friendly 2 x 4 validation plot.

    Uses hexbin density instead of raw scatter to avoid overplotting tens of
    thousands of admin-year event points.
    """
    if df_points.empty:
        return pd.DataFrame()
    df = df_points.copy()
    df["hazard"] = df["hazard"].astype(str).str.lower().str.strip()
    df["crop"] = df["crop"].astype(str).str.strip()
    df["obs_loss_pct"] = pd.to_numeric(df["obs_loss_pct"], errors="coerce").clip(lower=0, upper=max_loss)
    df[y_col] = pd.to_numeric(df[y_col], errors="coerce").clip(lower=0, upper=max_loss)
    df = df.dropna(subset=["obs_loss_pct", y_col, "hazard", "crop"])
    if df.empty:
        return pd.DataFrame()

    metrics = calc_metrics_for_ycol(df, y_col=y_col)
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["font.family"] = "DejaVu Sans"
    mpl.rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(2, 4, figsize=(17.6, 8.2), dpi=260, sharex=True, sharey=True)
    last_hb = None
    for i, haz in enumerate(HAZARDS):
        for j, crop in enumerate(CROPS):
            ax = axes[i, j]
            g = df[(df["hazard"] == haz) & (df["crop"] == crop)].copy()
            ax.plot([0, max_loss], [0, max_loss], "--", color="#202020", lw=1.0, zorder=1)
            if not g.empty:
                last_hb = ax.hexbin(
                    g["obs_loss_pct"],
                    g[y_col],
                    gridsize=gridsize,
                    extent=(0, max_loss, 0, max_loss),
                    mincnt=1,
                    bins="log",
                    cmap="YlGnBu",
                    linewidths=0.0,
                    zorder=2,
                )
                mm = metrics[(metrics["hazard"] == haz) & (metrics["crop"] == crop)]
                if not mm.empty:
                    r = mm.iloc[0]["pearson_r"]
                    txt = (
                        f"N={int(mm.iloc[0]['n'])}\n"
                        f"R={r:.2f}\n"
                        f"RMSE={mm.iloc[0]['rmse_loss_pct']:.1f}%\n"
                        f"Bias={mm.iloc[0]['mean_bias_pct_points']:.1f} pp"
                    )
                    ax.text(
                        0.03,
                        0.97,
                        txt,
                        transform=ax.transAxes,
                        ha="left",
                        va="top",
                        fontsize=8.5,
                        bbox=dict(facecolor="white", edgecolor="#bdbdbd", alpha=0.88, pad=3.0),
                    )
            else:
                ax.text(0.5, 0.5, "No data", transform=ax.transAxes, ha="center", va="center", color="#777777")
            ax.set_xlim(0, max_loss)
            ax.set_ylim(0, max_loss)
            ax.set_aspect("equal", adjustable="box")
            ax.grid(alpha=0.18, linewidth=0.6)
            ax.set_title(f"{HAZARD_LABEL[haz]} | {CROP_LABEL[crop]}", fontsize=11)
            if i == 1:
                ax.set_xlabel("Observed 5-year anomaly loss (%)", fontsize=10)
            if j == 0:
                ax.set_ylabel(y_label, fontsize=10)

    fig.suptitle(title, fontsize=14, y=0.975)
    fig.subplots_adjust(left=0.055, right=0.91, bottom=0.085, top=0.89, wspace=0.18, hspace=0.28)
    if last_hb is not None:
        cax = fig.add_axes([0.93, 0.18, 0.014, 0.64])
        cbar = fig.colorbar(last_hb, cax=cax)
        cbar.set_label("Event count per hexagon (log scale)", fontsize=9)
    out_png = Path(out_png)
    out_pdf = Path(out_pdf)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=320)
    fig.savefig(out_pdf)
    plt.close(fig)
    return metrics


def make_paper_validation_figures_from_csv(
    points_csv: Path = DEFAULT_OUT_DIR / "validation_points_hist_dry_wet.csv",
    out_dir: Path = DEFAULT_OUT_DIR,
    max_loss: float = 100.0,
    year_min: int | None = None,
    year_max: int | None = None,
) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(points_csv, low_memory=False)
    if year_min is not None or year_max is not None:
        df["year"] = pd.to_numeric(df["year"], errors="coerce")
        if year_min is not None:
            df = df[df["year"] >= int(year_min)]
        if year_max is not None:
            df = df[df["year"] <= int(year_max)]
    outputs: dict[str, Path] = {}
    year_tag = "" if year_min is None and year_max is None else f"_{year_min or 'start'}_{year_max or 'end'}"

    main_png = out_dir / f"paper_validation_sim_5yr_anomaly_loss_hexbin{year_tag}.png"
    main_pdf = out_dir / f"paper_validation_sim_5yr_anomaly_loss_hexbin{year_tag}.pdf"
    main_metrics = plot_validation_hexbin_by_crop_hazard(
        df,
        main_png,
        main_pdf,
        y_col="sim_loss_pct",
        y_label="Simulated 5-year anomaly loss (%)",
        title="Historical validation: matched 5-year anomaly yield loss",
        max_loss=max_loss,
    )
    main_metrics_csv = out_dir / f"paper_validation_sim_5yr_anomaly_loss_metrics{year_tag}.csv"
    main_metrics.to_csv(main_metrics_csv, index=False, encoding="utf-8-sig")
    outputs.update(
        {
            "main_anomaly_png": main_png,
            "main_anomaly_pdf": main_pdf,
            "main_anomaly_metrics_csv": main_metrics_csv,
            "anomaly_png": main_png,
            "anomaly_pdf": main_pdf,
            "anomaly_metrics_csv": main_metrics_csv,
        }
    )

    if "sim_direct_stress_loss_pct" in df.columns:
        direct_png = out_dir / f"paper_validation_direct_stress_loss_hexbin{year_tag}.png"
        direct_pdf = out_dir / f"paper_validation_direct_stress_loss_hexbin{year_tag}.pdf"
        direct_metrics = plot_validation_hexbin_by_crop_hazard(
            df,
            direct_png,
            direct_pdf,
            y_col="sim_direct_stress_loss_pct",
            y_label="Simulated direct water-stress loss (%)",
            title="Supplementary validation: direct water-stress yield loss",
            max_loss=max_loss,
        )
        direct_metrics_csv = out_dir / f"paper_validation_direct_stress_loss_metrics{year_tag}.csv"
        direct_metrics.to_csv(direct_metrics_csv, index=False, encoding="utf-8-sig")
        outputs.update({"direct_png": direct_png, "direct_pdf": direct_pdf, "direct_metrics_csv": direct_metrics_csv})
    return outputs


def _validation_level_ids(df: pd.DataFrame, level: str) -> pd.DataFrame:
    """Return a copy whose admin IDs represent the requested validation scale."""
    level = str(level).lower().strip()
    if level not in {"adm1", "country"}:
        raise ValueError("level must be 'adm1' or 'country'.")
    out = df.copy()
    out["iso"] = out["iso"].astype(str).str.upper().str.strip()
    out["adm1_id"] = _normalize_admin_id_series(out["adm1_id"], allow_missing=False)
    out["adm2_id"] = _normalize_admin_id_series(out["adm2_id"], allow_missing=True)
    if level == "country":
        out["adm1_id"] = ADMIN_MISSING
    out["adm2_id"] = ADMIN_MISSING
    return out


def _weighted_group_mean(df: pd.DataFrame, keys: list[str], value_col: str, weight_col: str) -> pd.DataFrame:
    v = pd.to_numeric(df[value_col], errors="coerce")
    w = pd.to_numeric(df[weight_col], errors="coerce")
    m = v.notna() & w.notna() & (w > 0)
    if not m.any():
        out = df[keys].drop_duplicates().copy()
        out[value_col] = np.nan
        return out
    tmp = df.loc[m, keys].copy()
    tmp["_num"] = v.loc[m].to_numpy(dtype=float) * w.loc[m].to_numpy(dtype=float)
    tmp["_den"] = w.loc[m].to_numpy(dtype=float)
    g = tmp.groupby(keys, dropna=False, as_index=False)[["_num", "_den"]].sum()
    g[value_col] = np.where(g["_den"] > 0, g["_num"] / g["_den"], np.nan)
    return g[keys + [value_col]]


def aggregate_model_obs_to_validation_level(
    model_admin_df: pd.DataFrame,
    obs_admin_df: pd.DataFrame,
    level: str = "adm1",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Aggregate admin2 validation time series to adm1 or country scale.

    This is mainly for scale-sensitivity checks. The gridded model is coarse, so
    validating directly at small admin units can add substantial spatial mismatch
    noise. Observed yields are weighted with the model crop-area weights because
    the observation files provide yield but not harvested area.
    """
    level = str(level).lower().strip()
    model = _validation_level_ids(model_admin_df, level)
    obs = obs_admin_df.copy()
    obs["iso"] = obs["iso"].astype(str).str.upper().str.strip()
    obs["adm1_id"] = _normalize_admin_id_series(obs["adm1_id"], allow_missing=False)
    obs["adm2_id"] = _normalize_admin_id_series(obs["adm2_id"], allow_missing=True)

    keys = ["crop", "iso", "adm1_id", "adm2_id", "year"]
    series_keys = ["crop", "iso", "adm1_id", "adm2_id"]
    admin_keys = ["crop", "iso", "adm1_id", "adm2_id"]

    model["crop_area_sum"] = pd.to_numeric(model["crop_area_sum"], errors="coerce")
    grouped = model.groupby(keys, dropna=False, as_index=False)
    model_out = grouped["crop_area_sum"].sum()

    weighted_model_cols = [
        "sim_yield_ws_kg_ha",
        "sim_yield_nws_kg_ha",
        "sim_direct_loss_pct",
        "model_loss_dry_pct_raw",
        "model_loss_wet_pct_raw",
        "model_loss_dry_pct",
        "model_loss_wet_pct",
        "dry_event_area_share",
        "wet_event_area_share",
    ]
    for col in weighted_model_cols:
        if col in model.columns:
            model_out = model_out.merge(_weighted_group_mean(model, keys, col, "crop_area_sum"), on=keys, how="left")

    for col in ["dry_event_grid_count", "wet_event_grid_count"]:
        if col in model.columns:
            s = grouped[col].sum().rename(columns={col: col})
            model_out = model_out.merge(s, on=keys, how="left")
    for col in ["is_dry_event_year", "is_wet_event_year"]:
        if col in model.columns:
            any_event = grouped[col].max().rename(columns={col: col})
            model_out = model_out.merge(any_event, on=keys, how="left")

    model_out = model_out.sort_values(series_keys + ["year"]).reset_index(drop=True)
    model_out["sim_baseline_kg_ha"] = model_out.groupby(series_keys, dropna=False)["sim_yield_ws_kg_ha"].transform(_centered_5yr_base)
    model_out["sim_anomaly_pct"] = (model_out["sim_yield_ws_kg_ha"] - model_out["sim_baseline_kg_ha"]) / model_out["sim_baseline_kg_ha"] * 100.0
    model_out["sim_loss_pct_raw"] = -model_out["sim_anomaly_pct"]
    model_out["sim_loss_pct"] = model_out["sim_loss_pct_raw"].clip(lower=0.0, upper=100.0)

    weights = (
        model_admin_df.copy()
        .assign(
            iso=lambda x: x["iso"].astype(str).str.upper().str.strip(),
            adm1_id=lambda x: _normalize_admin_id_series(x["adm1_id"], allow_missing=False),
            adm2_id=lambda x: _normalize_admin_id_series(x["adm2_id"], allow_missing=True),
            crop_area_sum=lambda x: pd.to_numeric(x["crop_area_sum"], errors="coerce"),
        )
        .groupby(admin_keys, dropna=False, as_index=False)["crop_area_sum"]
        .mean()
    )
    obs = obs.merge(weights, on=admin_keys, how="inner")
    obs = _validation_level_ids(obs, level)
    obs["obs_yield_kg_ha"] = pd.to_numeric(obs["obs_yield_kg_ha"], errors="coerce")
    obs["crop_area_sum"] = pd.to_numeric(obs["crop_area_sum"], errors="coerce")
    obs_out = _weighted_group_mean(obs, keys, "obs_yield_kg_ha", "crop_area_sum")
    obs_out = obs_out.sort_values(series_keys + ["year"]).reset_index(drop=True)
    obs_out["obs_baseline_kg_ha"] = obs_out.groupby(series_keys, dropna=False)["obs_yield_kg_ha"].transform(_centered_5yr_base)
    obs_out["obs_anomaly_pct"] = (obs_out["obs_yield_kg_ha"] - obs_out["obs_baseline_kg_ha"]) / obs_out["obs_baseline_kg_ha"] * 100.0
    obs_out["obs_loss_pct_raw"] = -obs_out["obs_anomaly_pct"]
    obs_out["obs_loss_pct"] = obs_out["obs_loss_pct_raw"].clip(lower=0.0, upper=100.0)
    return model_out, obs_out


def write_validation_scale_outputs(
    model_df: pd.DataFrame,
    obs_df: pd.DataFrame,
    out_dir: Path,
    level: str,
    year_min: int,
    year_max: int,
    admin_shp: Path = DEFAULT_ADMIN_SHP,
) -> dict:
    out_dir = Path(out_dir)
    model_level, obs_level = aggregate_model_obs_to_validation_level(model_df, obs_df, level=level)
    points_level = build_admin_hazard_comparison_points(model_level, obs_level)
    metrics_level = calc_metrics(points_level) if not points_level.empty else pd.DataFrame()

    model_csv = out_dir / f"model_{level}_year_hist_aggregated.csv"
    obs_csv = out_dir / f"obs_{level}_year_hist_with_5yr_anomaly.csv"
    points_csv = out_dir / f"validation_points_hist_dry_wet_{level}.csv"
    metrics_csv = out_dir / f"validation_metrics_hist_dry_wet_{level}.csv"
    scatter_png = out_dir / f"validation_hist_dry_wet_scatter_all_points_{level}.png"
    scatter_pdf = out_dir / f"validation_hist_dry_wet_scatter_all_points_{level}.pdf"
    hex_png = out_dir / f"paper_validation_sim_5yr_anomaly_loss_hexbin_{level}_{year_min}_{year_max}.png"
    hex_pdf = out_dir / f"paper_validation_sim_5yr_anomaly_loss_hexbin_{level}_{year_min}_{year_max}.pdf"
    hex_metrics_csv = out_dir / f"paper_validation_sim_5yr_anomaly_loss_metrics_{level}_{year_min}_{year_max}.csv"

    model_level.to_csv(model_csv, index=False, encoding="utf-8-sig")
    obs_level.to_csv(obs_csv, index=False, encoding="utf-8-sig")
    if not points_level.empty:
        points_level.to_csv(points_csv, index=False, encoding="utf-8-sig")
        plot_scatter_all_regions_all_stress_points(points_level, scatter_png, scatter_pdf)
        hex_metrics = plot_validation_hexbin_by_crop_hazard(
            points_level,
            hex_png,
            hex_pdf,
            y_col="sim_loss_pct",
            y_label="Simulated 5-year anomaly loss (%)",
            title=f"Historical validation at {level} scale: matched 5-year anomaly yield loss",
        )
        hex_metrics.to_csv(hex_metrics_csv, index=False, encoding="utf-8-sig")
        spatial_maps = plot_all_crop_spatial_correlation_maps_by_level(
            points_csv=points_csv,
            admin_shp=admin_shp,
            out_dir=out_dir,
            level=level,
            year_min=year_min,
            year_max=year_max,
        )
    else:
        hex_metrics_csv = None
        spatial_maps = {}
    if not metrics_level.empty:
        metrics_level.to_csv(metrics_csv, index=False, encoding="utf-8-sig")

    return {
        "level": level,
        "model_csv": model_csv,
        "obs_csv": obs_csv,
        "comparison_points_csv": points_csv,
        "metrics_csv": metrics_csv,
        "scatter_png": scatter_png if not points_level.empty else None,
        "scatter_pdf": scatter_pdf if not points_level.empty else None,
        "hexbin_png": hex_png if not points_level.empty else None,
        "hexbin_pdf": hex_pdf if not points_level.empty else None,
        "hexbin_metrics_csv": hex_metrics_csv,
        "spatial_maps": spatial_maps,
        "n_points": int(len(points_level)),
    }


def run_region_water_stress_loss_comparison_by_shp(
    DB: dict,
    out_dir: Path = DEFAULT_OUT_DIR,
    admin_shp: Path = DEFAULT_ADMIN_SHP,
    iso: str | None = None,
    adm1_id: int | None = None,
    adm2_id: int | None = None,
    fid: int | None = None,
    scen: str = "hist",
    gcm: str | None = None,
    year_min: int = DEFAULT_YEAR_MIN,
    year_max: int = DEFAULT_YEAR_MAX,
    yield_unit_factor: float = 1.0,
    relative_yield_scale: float = 1.0,
    dry_event_days_thr: float = 7.0,
    wet_event_days_thr: float = 7.0,
    event_area_share_thr: float = 0.2,
    obs_dir: Path = DEFAULT_OBS_DIR,
    crops: tuple[str, ...] = CROPS,
) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    region_gdf = select_region_polygons_from_shp(admin_shp, iso=iso, adm1_id=adm1_id, adm2_id=adm2_id, fid=fid)
    template = DB["dist"]["maiz_nonIrr"]
    gcm_tag = "NA" if scen == "hist" else str(gcm)
    tag = f"iso_{(iso or 'NA')}_adm1_{(adm1_id if adm1_id is not None else 'NA')}_adm2_{(adm2_id if adm2_id is not None else 'NA')}_fid_{(fid if fid is not None else 'NA')}_scen_{scen}_gcm_{gcm_tag}"
    region_lookup_csv = out_dir / f"grid_region_lookup_{tag}.csv"
    region_lookup = build_grid_region_lookup(region_gdf, template, cache_csv=region_lookup_csv)
    if region_lookup.empty:
        raise ValueError("Selected shp region has no matching DB grids.")

    model_parts, obs_parts, point_parts = [], [], []
    for crop in crops:
        m = compute_region_model_loss_timeseries_by_shp(
            DB, region_lookup, crop, scen=scen, gcm=gcm, year_min=year_min, year_max=year_max,
            yield_unit_factor=yield_unit_factor, relative_yield_scale=relative_yield_scale,
            dry_event_days_thr=dry_event_days_thr, wet_event_days_thr=wet_event_days_thr,
            event_area_share_thr=event_area_share_thr
        )
        o = load_observed_region_from_shp(crop, obs_dir, region_gdf, year_min, year_max)
        model_parts.append(m)
        obs_parts.append(o)
        if not m.empty and not o.empty:
            p = build_region_hazard_comparison_points(m, o)
            if not p.empty:
                point_parts.append(p)

    model_df = pd.concat(model_parts, ignore_index=True) if model_parts else pd.DataFrame()
    obs_df = pd.concat(obs_parts, ignore_index=True) if obs_parts else pd.DataFrame()
    if model_df.empty:
        raise ValueError(
            f"Model region timeseries is empty (scen={scen}, years={year_min}-{year_max}). "
            "For observed historical comparison, use scen='hist' (DB['hist'])."
        )
    points_df = pd.concat(point_parts, ignore_index=True) if point_parts else pd.DataFrame()
    metrics_df = calc_metrics(points_df) if not points_df.empty else pd.DataFrame()

    model_csv = out_dir / f"region_model_loss_timeseries_{tag}.csv"
    obs_csv = out_dir / f"region_obs_yield_timeseries_{tag}.csv"
    points_csv = out_dir / f"region_loss_comparison_points_{tag}.csv"
    metrics_csv = out_dir / f"region_loss_comparison_metrics_{tag}.csv"
    fig_png = out_dir / f"region_loss_comparison_scatter_{tag}.png"
    fig_pdf = out_dir / f"region_loss_comparison_scatter_{tag}.pdf"
    fig_all_png = out_dir / f"region_loss_comparison_all_points_{tag}.png"
    fig_all_pdf = out_dir / f"region_loss_comparison_all_points_{tag}.pdf"
    cfg_json = out_dir / f"region_loss_comparison_config_{tag}.json"

    model_df.to_csv(model_csv, index=False, encoding="utf-8-sig")
    obs_df.to_csv(obs_csv, index=False, encoding="utf-8-sig")
    if not points_df.empty:
        points_df.to_csv(points_csv, index=False, encoding="utf-8-sig")
        plot_scatter_all_regions_all_stress_points(points_df, fig_all_png, fig_all_pdf)
    if not metrics_df.empty:
        metrics_df.to_csv(metrics_csv, index=False, encoding="utf-8-sig")
        plot_scatter_hist_dry_wet(points_df, metrics_df, fig_png, fig_pdf)

    cfg = {
        "event_rule": "any_grid_with_crop_area_has_stress_days_ge_threshold",
        "scen": scen,
        "gcm": (None if scen == "hist" else gcm),
        "year_min": year_min,
        "year_max": year_max,
        "dry_event_days_thr": dry_event_days_thr,
        "wet_event_days_thr": wet_event_days_thr,
        "event_area_share_thr": event_area_share_thr,
        "event_area_share_thr_note": "Kept for compatibility; event-year flag uses any-grid rule in this version.",
        "yield_unit_factor": yield_unit_factor,
        "relative_yield_scale": relative_yield_scale,
        "iso": iso,
        "adm1_id": adm1_id,
        "adm2_id": adm2_id,
        "fid": fid,
        "n_region_grids": int(len(region_lookup)),
        "n_points": int(len(points_df)),
    }
    cfg_json.write_text(json.dumps(cfg, ensure_ascii=False, indent=2), encoding="utf-8")

    return {
        "model_timeseries_csv": model_csv,
        "obs_timeseries_csv": obs_csv,
        "comparison_points_csv": points_csv,
        "metrics_csv": metrics_csv,
        "figure_png": fig_png,
        "figure_pdf": fig_pdf,
        "figure_all_points_png": fig_all_png,
        "figure_all_points_pdf": fig_all_pdf,
        "config_json": cfg_json,
        "region_lookup_csv": region_lookup_csv,
        "n_points": int(len(points_df)),
    }


def run_validation_admin(
    DB: dict,
    year_min: int = DEFAULT_YEAR_MIN,
    year_max: int = DEFAULT_YEAR_MAX,
    dry_event_days_thr: float = 7.0,
    wet_event_days_thr: float = 7.0,
    admin_event_share_thr: float | None = None,
    event_area_share_thr: float = 0.2,
    yield_unit_factor: float = 1.0,
    relative_yield_scale: float = 1.0,
    obs_dir: Path = DEFAULT_OBS_DIR,
    admin_shp: Path = DEFAULT_ADMIN_SHP,
    out_dir: Path = DEFAULT_OUT_DIR,
    crops: tuple[str, ...] = CROPS,
    make_maize_spatial_corr: bool = True,
    make_macro_region_figures: bool = True,
    make_multiscale_validation: bool = True,
) -> dict:
    """
    Historical all-admin validation using DB['hist'].

    Output point unit: one admin unit, one year, one crop, one stress type.
    Event rule: any grid cell with crop area inside the admin unit reaches the
    stress-day threshold, so the whole admin-year is kept as an event point.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    obs_dir = Path(obs_dir)
    admin_shp = Path(admin_shp)
    if admin_event_share_thr is not None:
        event_area_share_thr = float(admin_event_share_thr)

    template = DB["dist"]["maiz_nonIrr"]
    lookup_csv = out_dir / "grid_admin_lookup_hist_all_regions_adminid_string_v2.csv"
    admin_lookup = build_grid_admin_lookup(admin_shp, template, cache_csv=lookup_csv)
    admin_lookup = admin_lookup[admin_lookup["rid"].notna()].copy()
    if admin_lookup.empty:
        raise ValueError("No DB grid points matched the admin shapefile.")

    model_parts = []
    obs_parts = []
    for crop in crops:
        model_parts.append(
            aggregate_model_crop_to_admin(
                DB=DB,
                admin_lookup_df=admin_lookup,
                crop=crop,
                year_min=year_min,
                year_max=year_max,
                yield_unit_factor=yield_unit_factor,
                relative_yield_scale=relative_yield_scale,
                dry_event_days_thr=dry_event_days_thr,
                wet_event_days_thr=wet_event_days_thr,
            )
        )
        obs_parts.append(load_observed_admin_yield(crop, obs_dir=obs_dir, year_min=year_min, year_max=year_max))

    model_df = pd.concat([p for p in model_parts if not p.empty], ignore_index=True) if model_parts else pd.DataFrame()
    obs_df = pd.concat([p for p in obs_parts if not p.empty], ignore_index=True) if obs_parts else pd.DataFrame()
    if model_df.empty:
        raise ValueError("Model aggregation is empty. Check DB['hist'], crop names, and the grid/admin overlay.")
    if obs_df.empty:
        raise ValueError(f"Observed admin yield is empty. Check obs_dir={obs_dir}.")

    points_df = build_admin_hazard_comparison_points(model_df, obs_df)
    metrics_df = calc_metrics(points_df) if not points_df.empty else pd.DataFrame()

    model_csv = out_dir / "model_admin_year_hist_aggregated.csv"
    obs_csv = out_dir / "obs_admin_year_hist_with_5yr_anomaly.csv"
    points_csv = out_dir / "validation_points_hist_dry_wet.csv"
    metrics_csv = out_dir / "validation_metrics_hist_dry_wet.csv"
    fig_png = out_dir / "validation_hist_dry_wet_scatter_by_crop.png"
    fig_pdf = out_dir / "validation_hist_dry_wet_scatter_by_crop.pdf"
    all_png = out_dir / "validation_hist_dry_wet_scatter_all_regions_all_points.png"
    all_pdf = out_dir / "validation_hist_dry_wet_scatter_all_regions_all_points.pdf"
    cfg_json = out_dir / "validation_run_config.json"

    model_df.to_csv(model_csv, index=False, encoding="utf-8-sig")
    obs_df.to_csv(obs_csv, index=False, encoding="utf-8-sig")
    if not points_df.empty:
        points_df.to_csv(points_csv, index=False, encoding="utf-8-sig")
        plot_scatter_all_regions_all_stress_points(points_df, all_png, all_pdf)
    if not metrics_df.empty:
        metrics_df.to_csv(metrics_csv, index=False, encoding="utf-8-sig")
        plot_scatter_hist_dry_wet(points_df, metrics_df, fig_png, fig_pdf)

    paper_figures = {}
    if not points_df.empty:
        paper_figures = make_paper_validation_figures_from_csv(
            points_csv=points_csv,
            out_dir=out_dir,
            year_min=year_min,
            year_max=year_max,
        )

    spatial_corr_maps = {}
    if make_maize_spatial_corr and not points_df.empty:
        spatial_corr_maps = plot_all_crop_spatial_correlation_maps(
            points_csv=points_csv,
            admin_shp=admin_shp,
            out_dir=out_dir,
            crops=crops,
            y_col="sim_loss_pct",
            year_min=year_min,
            year_max=year_max,
        )
    maize_corr = spatial_corr_maps.get("per_crop", {}).get("maiz", {}) if spatial_corr_maps else {}

    macro_region_figures = {}
    if make_macro_region_figures and not points_df.empty:
        macro_region_figures = plot_macro_region_validation_figures(
            points_csv=points_csv,
            admin_shp=admin_shp,
            out_dir=out_dir,
            y_col="sim_loss_pct",
            year_min=year_min,
            year_max=year_max,
        )

    multiscale_validation = {}
    if make_multiscale_validation:
        for level in ("adm1", "country"):
            multiscale_validation[level] = write_validation_scale_outputs(
                model_df=model_df,
                obs_df=obs_df,
                out_dir=out_dir,
                level=level,
                year_min=year_min,
                year_max=year_max,
                admin_shp=admin_shp,
            )

    cfg = {
        "scen": "hist",
        "db_source": "DB['hist']",
        "year_min": int(year_min),
        "year_max": int(year_max),
        "crops": list(crops),
        "dry_event_days_thr": float(dry_event_days_thr),
        "wet_event_days_thr": float(wet_event_days_thr),
        "event_rule": "any_grid_with_crop_area_has_stress_days_ge_threshold",
        "event_area_share_thr": float(event_area_share_thr),
        "event_area_share_thr_note": "Kept only for output/color diagnostics; event inclusion uses the any-grid rule.",
        "make_multiscale_validation": bool(make_multiscale_validation),
        "multiscale_validation_note": "Admin2 remains the native point scale. Additional adm1 and country outputs aggregate yields first, recompute 5-year anomalies, then rebuild event points to reduce grid/admin spatial mismatch noise.",
        "yield_unit_factor": float(yield_unit_factor),
        "relative_yield_scale": float(relative_yield_scale),
        "model_loss_method": "Main validation uses simulated admin yield 5-year centered moving-average anomaly; loss=max(0, -anomaly), clipped to 0-100%. Direct water-stress counterfactual loss is retained in sim_direct_stress_loss_pct for supplementary diagnostics.",
        "obs_loss_method": "Observed admin yield 5-year centered moving-average anomaly; loss=max(0, -anomaly), clipped to 0-100%.",
        "n_model_admin_year_rows": int(len(model_df)),
        "n_obs_admin_year_rows": int(len(obs_df)),
        "n_event_points": int(len(points_df)),
    }
    cfg_json.write_text(json.dumps(cfg, ensure_ascii=False, indent=2), encoding="utf-8")

    return {
        "model_admin_year_csv": model_csv,
        "obs_admin_year_csv": obs_csv,
        "comparison_points_csv": points_csv,
        "metrics_csv": metrics_csv,
        "figure_by_crop_png": fig_png,
        "figure_by_crop_pdf": fig_pdf,
        "figure_all_points_png": all_png,
        "figure_all_points_pdf": all_pdf,
        "paper_figures": paper_figures,
        "config_json": cfg_json,
        "grid_admin_lookup_csv": lookup_csv,
        "maize_spatial_corr": maize_corr,
        "spatial_corr_maps": spatial_corr_maps,
        "macro_region_figures": macro_region_figures,
        "multiscale_validation": multiscale_validation,
        "n_points": int(len(points_df)),
    }


def export_region_processed_yield_change(*args, **kwargs):
    raise NotImplementedError(
        "Use run_region_water_stress_loss_comparison_by_shp(...), then read model/obs CSV outputs."
    )


def export_region_processed_yield_change_by_shp(*args, **kwargs):
    return run_region_water_stress_loss_comparison_by_shp(*args, **kwargs)


def _spatial_corr_for_crop(
    points_csv: Path = DEFAULT_OUT_DIR / "validation_points_hist_dry_wet.csv",
    crop: str = "maiz",
    min_years: int = 5,
    y_col: str = "sim_loss_pct",
    use_clipped_loss: bool = True,
    year_min: int | None = DEFAULT_YEAR_MIN,
    year_max: int | None = DEFAULT_YEAR_MAX,
) -> pd.DataFrame:
    df = pd.read_csv(points_csv, low_memory=False)
    df = df[df["crop"].astype(str).str.strip() == crop].copy()
    df["hazard"] = df["hazard"].astype(str).str.lower().str.strip()
    df["iso"] = df["iso"].astype(str).str.upper().str.strip()
    df["adm1_id"] = _normalize_admin_id_series(df["adm1_id"], allow_missing=False)
    df["adm2_id"] = _normalize_admin_id_series(df["adm2_id"], allow_missing=True)
    df["year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")
    if year_min is not None:
        df = df[df["year"] >= int(year_min)]
    if year_max is not None:
        df = df[df["year"] <= int(year_max)]
    if y_col not in df.columns:
        raise KeyError(f"{y_col!r} not found in {points_csv}.")
    df[y_col] = pd.to_numeric(df[y_col], errors="coerce")
    df["obs_loss_pct"] = pd.to_numeric(df["obs_loss_pct"], errors="coerce")
    if use_clipped_loss:
        df[y_col] = df[y_col].clip(lower=0, upper=100)
        df["obs_loss_pct"] = df["obs_loss_pct"].clip(lower=0, upper=100)
    df = df.dropna(subset=["hazard", "iso", "adm1_id", "adm2_id", "year", y_col, "obs_loss_pct"])
    if df.empty:
        return pd.DataFrame(columns=["crop", "hazard", "iso", "adm1_id", "adm2_id", "n_years", "pearson_r"])
    rows = []
    for (haz, iso, a1, a2), g in df.groupby(["hazard", "iso", "adm1_id", "adm2_id"], dropna=False):
        if len(g) < int(min_years):
            continue
        x, y = g["obs_loss_pct"].values, g[y_col].values
        r = float(np.corrcoef(x, y)[0, 1]) if np.std(x) > 0 and np.std(y) > 0 else np.nan
        rows.append(
            {
                "crop": crop,
                "hazard": haz,
                "iso": iso,
                "adm1_id": a1,
                "adm2_id": a2,
                "n_years": int(len(g)),
                "pearson_r": r,
            }
        )
    return pd.DataFrame(rows)


def _spatial_year_tag(year_min: int | None, year_max: int | None) -> str:
    return "" if year_min is None and year_max is None else f"_{year_min or 'start'}_{year_max or 'end'}"


def _loss_column_tag(y_col: str) -> str:
    if y_col in {"sim_loss_pct", "sim_yield_anomaly_loss_pct"}:
        return "sim_5yr_anomaly_loss"
    if y_col == "sim_direct_stress_loss_pct":
        return "direct_stress_loss"
    if y_col == "sim_direct_stress_loss_pct_raw":
        return "direct_stress_loss_raw"
    return str(y_col).replace(" ", "_")


def _validation_y_label(y_col: str) -> str:
    if y_col in {"sim_loss_pct", "sim_yield_anomaly_loss_pct"}:
        return "Simulated 5-year anomaly loss (%)"
    if y_col in {"sim_direct_stress_loss_pct", "sim_direct_stress_loss_pct_raw"}:
        return "Simulated direct water-stress loss (%)"
    return str(y_col)


def plot_crop_spatial_correlation_all_regions(
    points_csv: Path = DEFAULT_OUT_DIR / "validation_points_hist_dry_wet.csv",
    admin_shp: Path = DEFAULT_ADMIN_SHP,
    out_dir: Path = DEFAULT_OUT_DIR,
    crop: str = "maiz",
    min_years: int = 5,
    y_col: str = "sim_loss_pct",
    use_clipped_loss: bool = True,
    year_min: int | None = DEFAULT_YEAR_MIN,
    year_max: int | None = DEFAULT_YEAR_MAX,
) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    crop = str(crop).strip()
    corr = _spatial_corr_for_crop(
        points_csv=points_csv,
        crop=crop,
        min_years=min_years,
        y_col=y_col,
        use_clipped_loss=use_clipped_loss,
        year_min=year_min,
        year_max=year_max,
    )
    loss_tag = _loss_column_tag(y_col)
    year_tag = _spatial_year_tag(year_min, year_max)
    corr_csv = out_dir / f"{crop}_spatial_corr_hist_dry_wet_all_regions_{loss_tag}{year_tag}.csv"
    corr.to_csv(corr_csv, index=False, encoding="utf-8-sig")
    if corr.empty:
        return {"corr_csv": corr_csv, "map_png": None, "map_pdf": None, "n_regions": 0}
    shp = _normalize_admin_gdf(gpd.read_file(admin_shp).to_crs("EPSG:4326"))
    key = shp[["iso", "adm1_id", "adm2_id", "geometry"]].copy()
    fig, axes = plt.subplots(1, 2, figsize=(17.5, 6.2), dpi=260)
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    for ax, haz in zip(axes, ["dry", "wet"]):
        sub = corr[corr["hazard"] == haz]
        m = key.merge(sub, on=["iso", "adm1_id", "adm2_id"], how="left")
        m.plot(
            ax=ax,
            column="pearson_r",
            cmap="RdBu_r",
            norm=norm,
            linewidth=0.0,
            edgecolor="none",
            missing_kwds={"color": "#efefef"},
            rasterized=True,
        )
        n_regions = int(sub["pearson_r"].notna().sum()) if not sub.empty else 0
        ax.set_title(f"{CROP_LABEL.get(crop, crop)} - {HAZARD_LABEL.get(haz, haz)} | n={n_regions}")
        ax.set_axis_off()
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="RdBu_r"), ax=axes.ravel().tolist(), fraction=0.025, pad=0.01)
    cbar.set_label("Pearson r")
    fig.suptitle(f"Spatial validation correlation | {CROP_LABEL.get(crop, crop)}", y=0.98)
    png_fp = out_dir / f"{crop}_spatial_corr_hist_dry_wet_all_regions_{loss_tag}{year_tag}.png"
    pdf_fp = out_dir / f"{crop}_spatial_corr_hist_dry_wet_all_regions_{loss_tag}{year_tag}.pdf"
    fig.savefig(png_fp, dpi=320, bbox_inches="tight")
    fig.savefig(pdf_fp, bbox_inches="tight")
    plt.close(fig)
    return {"corr_csv": corr_csv, "map_png": png_fp, "map_pdf": pdf_fp, "n_regions": int(len(corr))}


def plot_all_crop_spatial_correlation_maps(
    points_csv: Path = DEFAULT_OUT_DIR / "validation_points_hist_dry_wet.csv",
    admin_shp: Path = DEFAULT_ADMIN_SHP,
    out_dir: Path = DEFAULT_OUT_DIR,
    crops: tuple[str, ...] = CROPS,
    min_years: int = 5,
    y_col: str = "sim_loss_pct",
    use_clipped_loss: bool = True,
    year_min: int | None = DEFAULT_YEAR_MIN,
    year_max: int | None = DEFAULT_YEAR_MAX,
) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    loss_tag = _loss_column_tag(y_col)
    year_tag = _spatial_year_tag(year_min, year_max)

    per_crop = {}
    corr_parts = []
    for crop in crops:
        per_crop[crop] = plot_crop_spatial_correlation_all_regions(
            points_csv=points_csv,
            admin_shp=admin_shp,
            out_dir=out_dir,
            crop=crop,
            min_years=min_years,
            y_col=y_col,
            use_clipped_loss=use_clipped_loss,
            year_min=year_min,
            year_max=year_max,
        )
        corr_csv = per_crop[crop].get("corr_csv")
        if corr_csv is not None and Path(corr_csv).exists():
            corr_parts.append(pd.read_csv(corr_csv, low_memory=False))

    corr = pd.concat(corr_parts, ignore_index=True) if corr_parts else pd.DataFrame()
    all_corr_csv = out_dir / f"all_crops_spatial_corr_hist_dry_wet_all_regions_{loss_tag}{year_tag}.csv"
    corr.to_csv(all_corr_csv, index=False, encoding="utf-8-sig")
    if corr.empty:
        return {"per_crop": per_crop, "combined_corr_csv": all_corr_csv, "combined_png": None, "combined_pdf": None}

    shp = _normalize_admin_gdf(gpd.read_file(admin_shp).to_crs("EPSG:4326"))
    key = shp[["iso", "adm1_id", "adm2_id", "geometry"]].copy()
    fig, axes = plt.subplots(len(crops), 2, figsize=(16.5, 4.2 * len(crops)), dpi=230)
    if len(crops) == 1:
        axes = np.array([axes])
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    for i, crop in enumerate(crops):
        for j, haz in enumerate(["dry", "wet"]):
            ax = axes[i, j]
            sub = corr[(corr["crop"] == crop) & (corr["hazard"] == haz)].copy()
            m = key.merge(sub, on=["iso", "adm1_id", "adm2_id"], how="left")
            m.plot(
                ax=ax,
                column="pearson_r",
                cmap="RdBu_r",
                norm=norm,
                linewidth=0.0,
                edgecolor="none",
                missing_kwds={"color": "#efefef"},
                rasterized=True,
            )
            n_regions = int(sub["pearson_r"].notna().sum()) if not sub.empty else 0
            ax.set_title(f"{CROP_LABEL.get(crop, crop)} | {HAZARD_LABEL.get(haz, haz)} | n={n_regions}", fontsize=11)
            ax.set_axis_off()
    fig.subplots_adjust(left=0.02, right=0.92, top=0.94, bottom=0.03, wspace=0.03, hspace=0.12)
    cax = fig.add_axes([0.935, 0.18, 0.014, 0.64])
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="RdBu_r"), cax=cax)
    cbar.set_label("Pearson r")
    fig.suptitle("Spatial validation correlation by crop and water stress", y=0.98, fontsize=15)
    combined_png = out_dir / f"all_crops_spatial_corr_hist_dry_wet_all_regions_{loss_tag}{year_tag}.png"
    combined_pdf = out_dir / f"all_crops_spatial_corr_hist_dry_wet_all_regions_{loss_tag}{year_tag}.pdf"
    fig.savefig(combined_png, dpi=320, bbox_inches="tight")
    fig.savefig(combined_pdf, bbox_inches="tight")
    plt.close(fig)
    return {
        "per_crop": per_crop,
        "combined_corr_csv": all_corr_csv,
        "combined_png": combined_png,
        "combined_pdf": combined_pdf,
        "n_regions": int(len(corr)),
    }


def _validation_level_key_cols(level: str) -> list[str]:
    level = str(level).lower().strip()
    if level in {"admin2", "adm2"}:
        return ["iso", "adm1_id", "adm2_id"]
    if level == "adm1":
        return ["iso", "adm1_id"]
    if level == "country":
        return ["iso"]
    raise ValueError("level must be 'admin2', 'adm1', or 'country'.")


def _spatial_geometry_for_validation_level(admin_shp: Path, level: str) -> gpd.GeoDataFrame:
    level = str(level).lower().strip()
    shp = _normalize_admin_gdf(gpd.read_file(admin_shp).to_crs("EPSG:4326"))
    if level in {"admin2", "adm2"}:
        return shp[["iso", "adm1_id", "adm2_id", "geometry"]].copy()
    if level == "adm1":
        return shp.dissolve(by=["iso", "adm1_id"], as_index=False)[["iso", "adm1_id", "geometry"]]
    if level == "country":
        return shp.dissolve(by=["iso"], as_index=False)[["iso", "geometry"]]
    raise ValueError("level must be 'admin2', 'adm1', or 'country'.")


def _spatial_corr_for_crop_level(
    points_csv: Path,
    crop: str,
    level: str = "adm1",
    min_years: int = 5,
    y_col: str = "sim_loss_pct",
    use_clipped_loss: bool = True,
    year_min: int | None = DEFAULT_YEAR_MIN,
    year_max: int | None = DEFAULT_YEAR_MAX,
) -> pd.DataFrame:
    key_cols = _validation_level_key_cols(level)
    df = pd.read_csv(points_csv, low_memory=False)
    df = df[df["crop"].astype(str).str.strip() == crop].copy()
    df["hazard"] = df["hazard"].astype(str).str.lower().str.strip()
    df["iso"] = df["iso"].astype(str).str.upper().str.strip()
    if "adm1_id" in df.columns:
        df["adm1_id"] = _normalize_admin_id_series(df["adm1_id"], allow_missing=False)
    if "adm2_id" in df.columns:
        df["adm2_id"] = _normalize_admin_id_series(df["adm2_id"], allow_missing=True)
    df["year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")
    if year_min is not None:
        df = df[df["year"] >= int(year_min)]
    if year_max is not None:
        df = df[df["year"] <= int(year_max)]
    if y_col not in df.columns:
        raise KeyError(f"{y_col!r} not found in {points_csv}.")
    df[y_col] = pd.to_numeric(df[y_col], errors="coerce")
    df["obs_loss_pct"] = pd.to_numeric(df["obs_loss_pct"], errors="coerce")
    if use_clipped_loss:
        df[y_col] = df[y_col].clip(lower=0, upper=100)
        df["obs_loss_pct"] = df["obs_loss_pct"].clip(lower=0, upper=100)
    df = df.dropna(subset=["hazard", *key_cols, "year", y_col, "obs_loss_pct"])
    rows = []
    for group_key, g in df.groupby(["hazard", *key_cols], dropna=False):
        if len(g) < int(min_years):
            continue
        if not isinstance(group_key, tuple):
            group_key = (group_key,)
        haz = group_key[0]
        ids = dict(zip(key_cols, group_key[1:]))
        x, y = g["obs_loss_pct"].to_numpy(dtype=float), g[y_col].to_numpy(dtype=float)
        r = float(np.corrcoef(x, y)[0, 1]) if len(g) > 1 and np.std(x) > 0 and np.std(y) > 0 else np.nan
        rows.append({"crop": crop, "hazard": haz, **ids, "n_years": int(len(g)), "pearson_r": r})
    return pd.DataFrame(rows)


def plot_all_crop_spatial_correlation_maps_by_level(
    points_csv: Path,
    admin_shp: Path = DEFAULT_ADMIN_SHP,
    out_dir: Path = DEFAULT_OUT_DIR,
    level: str = "adm1",
    crops: tuple[str, ...] = CROPS,
    min_years: int = 5,
    y_col: str = "sim_loss_pct",
    use_clipped_loss: bool = True,
    year_min: int | None = DEFAULT_YEAR_MIN,
    year_max: int | None = DEFAULT_YEAR_MAX,
) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    level = str(level).lower().strip()
    key_cols = _validation_level_key_cols(level)
    loss_tag = _loss_column_tag(y_col)
    year_tag = _spatial_year_tag(year_min, year_max)

    corr_parts = []
    for crop in crops:
        corr_parts.append(
            _spatial_corr_for_crop_level(
                points_csv=points_csv,
                crop=crop,
                level=level,
                min_years=min_years,
                y_col=y_col,
                use_clipped_loss=use_clipped_loss,
                year_min=year_min,
                year_max=year_max,
            )
        )
    corr = pd.concat([p for p in corr_parts if not p.empty], ignore_index=True) if corr_parts else pd.DataFrame()
    corr_csv = out_dir / f"all_crops_spatial_corr_hist_dry_wet_{level}_{loss_tag}{year_tag}.csv"
    corr.to_csv(corr_csv, index=False, encoding="utf-8-sig")
    if corr.empty:
        return {"combined_corr_csv": corr_csv, "combined_png": None, "combined_pdf": None, "n_regions": 0}

    key = _spatial_geometry_for_validation_level(admin_shp, level)
    fig, axes = plt.subplots(len(crops), 2, figsize=(16.5, 4.2 * len(crops)), dpi=230)
    if len(crops) == 1:
        axes = np.array([axes])
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    for i, crop in enumerate(crops):
        for j, haz in enumerate(["dry", "wet"]):
            ax = axes[i, j]
            sub = corr[(corr["crop"] == crop) & (corr["hazard"] == haz)].copy()
            m = key.merge(sub, on=key_cols, how="left")
            m.plot(
                ax=ax,
                column="pearson_r",
                cmap="RdBu_r",
                norm=norm,
                linewidth=0.0,
                edgecolor="none",
                missing_kwds={"color": "#efefef"},
                rasterized=True,
            )
            n_regions = int(sub["pearson_r"].notna().sum()) if not sub.empty else 0
            ax.set_title(f"{CROP_LABEL.get(crop, crop)} | {HAZARD_LABEL.get(haz, haz)} | n={n_regions}", fontsize=11)
            ax.set_axis_off()
    fig.subplots_adjust(left=0.02, right=0.92, top=0.94, bottom=0.03, wspace=0.03, hspace=0.12)
    cax = fig.add_axes([0.935, 0.18, 0.014, 0.64])
    cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="RdBu_r"), cax=cax)
    cbar.set_label("Pearson r")
    fig.suptitle(f"Spatial validation correlation by crop and water stress | {level}", y=0.98, fontsize=15)
    combined_png = out_dir / f"all_crops_spatial_corr_hist_dry_wet_{level}_{loss_tag}{year_tag}.png"
    combined_pdf = out_dir / f"all_crops_spatial_corr_hist_dry_wet_{level}_{loss_tag}{year_tag}.pdf"
    fig.savefig(combined_png, dpi=320, bbox_inches="tight")
    fig.savefig(combined_pdf, bbox_inches="tight")
    plt.close(fig)
    return {"combined_corr_csv": corr_csv, "combined_png": combined_png, "combined_pdf": combined_pdf, "n_regions": int(len(corr))}


def plot_maize_spatial_correlation_all_regions(
    points_csv: Path = DEFAULT_OUT_DIR / "validation_points_hist_dry_wet.csv",
    admin_shp: Path = DEFAULT_ADMIN_SHP,
    out_dir: Path = DEFAULT_OUT_DIR,
    min_years: int = 5,
    use_clipped_loss: bool = True,
) -> dict:
    return plot_crop_spatial_correlation_all_regions(
        points_csv=points_csv,
        admin_shp=admin_shp,
        out_dir=out_dir,
        crop="maiz",
        min_years=min_years,
        y_col="sim_loss_pct",
        use_clipped_loss=use_clipped_loss,
        year_min=DEFAULT_YEAR_MIN,
        year_max=DEFAULT_YEAR_MAX,
    )


def add_macro_region(df: pd.DataFrame, iso_col: str = "iso") -> pd.DataFrame:
    out = df.copy()
    out[iso_col] = out[iso_col].astype(str).str.upper().str.strip()
    out["macro_region"] = out[iso_col].map(ISO3_TO_MACRO_REGION).fillna("Other")
    return out


def _macro_region_validation_data(
    points_csv: Path = DEFAULT_OUT_DIR / "validation_points_hist_dry_wet.csv",
    y_col: str = "sim_loss_pct",
    year_min: int | None = DEFAULT_YEAR_MIN,
    year_max: int | None = DEFAULT_YEAR_MAX,
    use_clipped_loss: bool = True,
) -> pd.DataFrame:
    df = pd.read_csv(points_csv, low_memory=False)
    df = add_macro_region(df)
    df["hazard"] = df["hazard"].astype(str).str.lower().str.strip()
    df["crop"] = df["crop"].astype(str).str.strip()
    df["adm1_id"] = _normalize_admin_id_series(df["adm1_id"], allow_missing=False)
    df["adm2_id"] = _normalize_admin_id_series(df["adm2_id"], allow_missing=True)
    df["year"] = pd.to_numeric(df["year"], errors="coerce").astype("Int64")
    if y_col not in df.columns:
        raise KeyError(f"{y_col!r} not found in {points_csv}.")
    df["obs_loss_pct"] = pd.to_numeric(df["obs_loss_pct"], errors="coerce")
    df[y_col] = pd.to_numeric(df[y_col], errors="coerce")
    if year_min is not None:
        df = df[df["year"] >= int(year_min)]
    if year_max is not None:
        df = df[df["year"] <= int(year_max)]
    if use_clipped_loss:
        df["obs_loss_pct"] = df["obs_loss_pct"].clip(lower=0, upper=100)
        df[y_col] = df[y_col].clip(lower=0, upper=100)
    return df.dropna(subset=["macro_region", "hazard", "crop", "iso", "adm1_id", "adm2_id", "year", "obs_loss_pct", y_col])


def _macro_region_spatial_corr(df: pd.DataFrame, region: str, hazard: str, y_col: str, min_points: int = 5) -> pd.DataFrame:
    sub = df[(df["macro_region"] == region) & (df["hazard"] == hazard)].copy()
    rows = []
    for (iso, a1, a2), g in sub.groupby(["iso", "adm1_id", "adm2_id"], dropna=False):
        if len(g) < int(min_points):
            continue
        x = pd.to_numeric(g["obs_loss_pct"], errors="coerce").to_numpy(dtype=float)
        y = pd.to_numeric(g[y_col], errors="coerce").to_numpy(dtype=float)
        r = float(np.corrcoef(x, y)[0, 1]) if len(g) > 1 and np.std(x) > 0 and np.std(y) > 0 else np.nan
        rows.append({"macro_region": region, "hazard": hazard, "iso": iso, "adm1_id": a1, "adm2_id": a2, "n_points": int(len(g)), "pearson_r": r})
    return pd.DataFrame(rows)


def _plot_region_map_panel(
    ax,
    key: gpd.GeoDataFrame,
    corr: pd.DataFrame,
    region: str,
    hazard: str,
    norm,
    base_key: gpd.GeoDataFrame | None = None,
) -> None:
    if base_key is not None:
        base_region = base_key[base_key["macro_region"] == region].copy()
        if not base_region.empty:
            base_region.plot(ax=ax, color="#eeeeee", linewidth=0.0, edgecolor="none", rasterized=True)
    region_key = key[key["macro_region"] == region].copy()
    if region_key.empty:
        title = f"{HAZARD_LABEL.get(hazard, hazard)} spatial r | n=0"
        ax.set_title(title, fontsize=10)
        if base_key is None or base_key[base_key["macro_region"] == region].empty:
            ax.text(0.5, 0.5, "No region geometry", transform=ax.transAxes, ha="center", va="center", color="#777777")
        ax.set_axis_off()
        return
    sub = corr[(corr["macro_region"] == region) & (corr["hazard"] == hazard)].copy()
    m = region_key.merge(sub, on=["iso", "adm1_id", "adm2_id"], how="left")
    m.plot(
        ax=ax,
        column="pearson_r",
        cmap="RdBu_r",
        norm=norm,
        linewidth=0.0,
        edgecolor="none",
            missing_kwds={"color": "#eeeeee"},
            rasterized=True,
        )
    n_valid = int(sub["pearson_r"].notna().sum()) if not sub.empty else 0
    ax.set_title(f"{HAZARD_LABEL.get(hazard, hazard)} spatial r | n={n_valid}", fontsize=10)
    ax.set_axis_off()


def _plot_region_scatter_panel(ax, df: pd.DataFrame, region: str, hazard: str, y_col: str) -> dict:
    sub = df[(df["macro_region"] == region) & (df["hazard"] == hazard)].copy()
    colors = {"maiz": "#D55E00", "soyb": "#009E73", "rice": "#0072B2", "whea": "#CC79A7"}
    markers = {"maiz": "o", "soyb": "s", "rice": "^", "whea": "D"}
    ax.plot([0, 100], [0, 100], "--", color="#262626", lw=1.0)
    if sub.empty:
        ax.text(0.5, 0.5, "No validation points", transform=ax.transAxes, ha="center", va="center", color="#777777")
        metrics = {"n": 0, "pearson_r": np.nan, "rmse_loss_pct": np.nan, "mean_bias_pct_points": np.nan}
    else:
        for crop in CROPS:
            g = sub[sub["crop"] == crop]
            if g.empty:
                continue
            ax.scatter(
                g["obs_loss_pct"],
                g[y_col],
                s=10,
                alpha=0.32,
                c=colors.get(crop, "#555555"),
                marker=markers.get(crop, "o"),
                edgecolor="none",
                label=CROP_LABEL.get(crop, crop),
            )
        metrics = _loss_metrics(sub["obs_loss_pct"], sub[y_col])
        r = metrics["pearson_r"]
        txt = (
            f"N={metrics['n']}\n"
            f"R={r:.2f}\n"
            f"RMSE={metrics['rmse_loss_pct']:.1f}%\n"
            f"Bias={metrics['mean_bias_pct_points']:.1f} pp"
        )
        ax.text(
            0.03,
            0.97,
            txt,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.5,
            bbox=dict(facecolor="white", edgecolor="#bdbdbd", alpha=0.9, pad=3.0),
        )
    ax.set_title(f"{HAZARD_LABEL.get(hazard, hazard)} event scatter", fontsize=10)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(alpha=0.20, linewidth=0.6)
    ax.set_xlabel("Observed 5-year anomaly loss (%)", fontsize=9)
    ax.set_ylabel(_validation_y_label(y_col), fontsize=9)
    return metrics


def plot_macro_region_validation_figures(
    points_csv: Path = DEFAULT_OUT_DIR / "validation_points_hist_dry_wet.csv",
    admin_shp: Path = DEFAULT_ADMIN_SHP,
    out_dir: Path = DEFAULT_OUT_DIR,
    macro_regions: tuple[str, ...] = MACRO_REGION_ORDER,
    y_col: str = "sim_loss_pct",
    min_spatial_points: int = 5,
    year_min: int | None = DEFAULT_YEAR_MIN,
    year_max: int | None = DEFAULT_YEAR_MAX,
    use_clipped_loss: bool = True,
) -> dict:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    df = _macro_region_validation_data(
        points_csv=points_csv,
        y_col=y_col,
        year_min=year_min,
        year_max=year_max,
        use_clipped_loss=use_clipped_loss,
    )
    raw_shp = gpd.read_file(admin_shp).to_crs("EPSG:4326")
    raw_shp["iso"] = raw_shp["iso"].astype(str).str.upper().str.strip()
    base_key = add_macro_region(raw_shp[["iso", "geometry"]].copy())
    shp = _normalize_admin_gdf(raw_shp)
    key = add_macro_region(shp[["iso", "adm1_id", "adm2_id", "geometry"]].copy())
    loss_tag = _loss_column_tag(y_col)
    year_tag = _spatial_year_tag(year_min, year_max)

    corr_parts = []
    metric_rows = []
    for region in macro_regions:
        for haz in HAZARDS:
            c = _macro_region_spatial_corr(df, region=region, hazard=haz, y_col=y_col, min_points=min_spatial_points)
            corr_parts.append(c)
    corr = pd.concat(corr_parts, ignore_index=True) if corr_parts else pd.DataFrame()
    corr_csv = out_dir / f"macro_region_spatial_corr_{loss_tag}{year_tag}.csv"
    corr.to_csv(corr_csv, index=False, encoding="utf-8-sig")

    outputs = {"corr_csv": corr_csv, "regions": {}}
    norm = TwoSlopeNorm(vmin=-1, vcenter=0, vmax=1)
    for region in macro_regions:
        fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.4), dpi=250)
        for i, haz in enumerate(HAZARDS):
            _plot_region_map_panel(axes[i, 0], key=key, corr=corr, region=region, hazard=haz, norm=norm, base_key=base_key)
            metrics = _plot_region_scatter_panel(axes[i, 1], df=df, region=region, hazard=haz, y_col=y_col)
            metric_rows.append({"macro_region": region, "hazard": haz, "sim_loss_column": y_col, **metrics})
            if i == 0:
                handles, labels = axes[i, 1].get_legend_handles_labels()
                if handles:
                    axes[i, 1].legend(handles, labels, loc="lower right", fontsize=8, frameon=True)
        fig.subplots_adjust(left=0.04, right=0.90, bottom=0.07, top=0.90, wspace=0.12, hspace=0.25)
        cax = fig.add_axes([0.92, 0.20, 0.018, 0.60])
        cbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="RdBu_r"), cax=cax)
        cbar.set_label("Admin-unit Pearson r", fontsize=9)
        fig.suptitle(f"{region} validation by water-stress type", fontsize=14, y=0.965)
        safe_region = region.lower().replace(" ", "_")
        png_fp = out_dir / f"macro_region_validation_{safe_region}_{loss_tag}{year_tag}.png"
        pdf_fp = out_dir / f"macro_region_validation_{safe_region}_{loss_tag}{year_tag}.pdf"
        fig.savefig(png_fp, dpi=320, bbox_inches="tight")
        fig.savefig(pdf_fp, bbox_inches="tight")
        plt.close(fig)
        outputs["regions"][region] = {"png": png_fp, "pdf": pdf_fp}

    metrics = pd.DataFrame(metric_rows)
    metrics_csv = out_dir / f"macro_region_scatter_metrics_{loss_tag}{year_tag}.csv"
    metrics.to_csv(metrics_csv, index=False, encoding="utf-8-sig")
    outputs["metrics_csv"] = metrics_csv
    return outputs


if __name__ == "__main__":
    raise RuntimeError(
        "Use in notebook:\n"
        "import importlib\n"
        "import validate_admin_custom_criteria as vac\n"
        "vac = importlib.reload(vac)\n"
        "out = vac.run_region_water_stress_loss_comparison_by_shp(DB, scen='hist', iso='CHN', adm1_id=44, adm2_id=441881)\n"
        "print(out)\n"
    )
