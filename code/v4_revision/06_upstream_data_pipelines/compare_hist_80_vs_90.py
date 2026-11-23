import os
import numpy as np
import pandas as pd
import xarray as xr
import rioxarray as rxr

crops = ["maiz", "soyb", "rice", "whea"]
modes = ["Irr", "nonIrr"]
hazards = ["dry", "wet"]

cellarea = xr.open_dataset(r"./data/external/cellarea.nc").cellarea30min

def open_da(path, varname):
    ds = xr.open_dataset(path, decode_times=False)
    da = ds[varname]
    t = da["time"].values.astype(float)
    units = ds["time"].attrs.get("units", "").lower()
    if "since" in units and "year" in units:
        base_year = int(units.split("since")[1].strip()[:4])
        year = (base_year + np.rint(t)).astype(int) if np.nanmin(t) < 1000 else np.rint(t).astype(int)
        da = da.assign_coords(time=("time", year)).rename({"time": "year"})
    else:
        da = da.rename({"time": "year"})
    return da.where(da != 0)

def read_crop(fp, crop, subdir, fname_template, varname):
    if crop != "whea":
        return open_da(os.path.join(fp, subdir, fname_template.format(crop=crop)), varname)
    da_wwh = open_da(os.path.join(fp, subdir, fname_template.format(crop="wwh")), varname)
    da_swh = open_da(os.path.join(fp, subdir, fname_template.format(crop="swh")), varname)
    da = da_wwh.fillna(0) + da_swh.fillna(0)
    return da.where(da != 0)

# Read distribution
AR = r"./data/cwatm_inputs_figshare"
tif_map = {
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
mirca_out = {}
for key, rel in tif_map.items():
    da = rxr.open_rasterio(os.path.join(AR, rel)).squeeze("band", drop=True).rename({"y": "lat", "x": "lon"})
    da = da.interp(lat=cellarea.lat, lon=cellarea.lon, method="linear") * cellarea
    mirca_out[key] = da
dist = xr.Dataset(mirca_out)
whea_Irr = dist["wwh_Irr"].fillna(0) + dist["swh_Irr"].fillna(0)
whea_nonIrr = dist["wwh_nonIrr"].fillna(0) + dist["swh_nonIrr"].fillna(0)
dist = dist.assign(whea_Irr=whea_Irr.where(whea_Irr > 0), whea_nonIrr=whea_nonIrr.where(whea_nonIrr > 0))

# Potential yield
pot_dir = r"./data/gaez_v4_200b_noco2/hist"
pot = {}
for crop in crops:
    pot[crop] = {}
    for mode in modes:
        if crop == "rice":
            da = xr.open_dataset(os.path.join(pot_dir, f"ricw_{mode}_Hist_CRUTS32_8110L.nc"))["yield"]
        elif crop == "whea":
            da1 = xr.open_dataset(os.path.join(pot_dir, f"wwhe_{mode}_Hist_CRUTS32_8110L.nc"))["yield"]
            da2 = xr.open_dataset(os.path.join(pot_dir, f"swhe_{mode}_Hist_CRUTS32_8110L.nc"))["yield"]
            da = da1.fillna(0) + da2.fillna(0)
        else:
            da = xr.open_dataset(os.path.join(pot_dir, f"{crop}_{mode}_Hist_CRUTS32_8110L.nc"))["yield"]
        da_i = da.interp(lat=cellarea.lat, lon=cellarea.lon, method="linear")
        pot[crop][mode] = da_i.where(da_i != 0)

fp = r"./data/cwatm_outputs/gswp3_long"
def compute_hist(years):
    out = []
    for crop in crops:
        ws_rf = 4 - read_crop(fp, crop, "waterstress", "gs_loss_nonIrr_annualtot_{crop}.nc", "gs_loss_nonIrr_annualtot").sel(year=years)
        nws_rf = 4 - read_crop(fp, crop, "nowaterstress", "gs_loss_nonIrr_annualtot_{crop}.nc", "gs_loss_nonIrr_annualtot").sel(year=years)
        ws_ir = 4 - read_crop(fp, crop, "waterstress", "gs_loss_Irr_annualtot_{crop}.nc", "gs_loss_Irr_annualtot").sel(year=years)
        nws_ir = 4 - read_crop(fp, crop, "nowaterstress", "gs_loss_Irr_annualtot_{crop}.nc", "gs_loss_Irr_annualtot").sel(year=years)
        A_rf = dist[f"{crop}_nonIrr"] / 10000.0
        A_ir = dist[f"{crop}_Irr"] / 10000.0
        pot_rf = (pot[crop]["nonIrr"] * A_rf).fillna(0)
        pot_ir = (pot[crop]["Irr"] * A_ir).fillna(0)
        a2, b2 = xr.align((nws_rf * pot_rf).fillna(0), (nws_ir * pot_ir).fillna(0), join="outer", fill_value=0)
        den_t = a2 + b2
        for haz in hazards:
            ev = read_crop(fp, crop, "nowaterstress", f"{haz}stress_nonIrr_annualtot_{{crop}}.nc", f"{haz}stress_nonIrr_annualtot").sel(year=years)
            mask = ev > 0
            loss_rf = (nws_rf - ws_rf).where(mask) * pot_rf
            loss_ir = (nws_ir - ws_ir).where(mask) * pot_ir
            loss_t = loss_rf.where(loss_rf > 0).fillna(0) + loss_ir.where(loss_ir > 0).fillna(0)
            lp = ((loss_t / (den_t + 1e-12)) * 100.0).where(den_t > 0).mean("year", skipna=True)
            lp = lp.where(np.isfinite(lp)).astype("float32").expand_dims(crop=[crop], haz=[haz])
            out.append(lp)
    return xr.combine_by_coords(out)

h80 = compute_hist(slice(1980, 2019))
h90 = compute_hist(slice(1990, 2019))

ref = xr.open_dataset(r"./outputs/figures/loss_maps/loss_hist_1980_2019.nc")["loss_pct"]
diff_ref = float(np.nanmax(np.abs(h80.sel(crop=ref.crop, haz=ref.haz).values - ref.values)))
print(f"Max abs diff vs saved loss_hist_1980_2019.nc: {diff_ref:.6e}")

# Save h90 to scratch for downstream checks
out_h90 = r"./outputs/scratch/loss_hist_1990_2019.nc"
xr.Dataset({"loss_pct": h90}).to_netcdf(out_h90)

# Compare Fig 5 bin stats!
fut = xr.open_dataset(r"./outputs/figures/loss_maps/loss_future_2070_2100_ensemble_mean.nc")["loss_pct"]
area_list = [(dist[f"{c}_Irr"].fillna(0) + dist[f"{c}_nonIrr"].fillna(0)).rename("area").expand_dims(crop=[c]) for c in crops]
area = xr.concat(area_list, dim="crop").interp(lat=fut.lat, lon=fut.lon, method="nearest").fillna(0.0)
den_area = area.sum("crop")

for tag, h in [("1980-2019", h80), ("1990-2019", h90)]:
    h_s = h.interp(lat=fut.lat, lon=fut.lon, method="nearest")
    loss_all = xr.concat([h_s.expand_dims(scen=["hist"]), fut], dim="scen")
    if loss_all.lat[0] > loss_all.lat[-1]:
        loss_all = loss_all.sortby("lat")
        area_s = area.sortby("lat")
        den_s = den_area.sortby("lat")
    else:
        area_s, den_s = area, den_area
    wmean = xr.where(den_s > 0, (loss_all * area_s).sum("crop", skipna=True) / den_s, np.nan)
    delta = (wmean.sel(scen=["ssp126", "ssp585"]) - wmean.sel(scen="hist")).where(den_s > 0)
    print(f"\n=== Fig 5 Area Shares ({tag}) ===")
    tot_a = float(den_s.where(den_s > 0).sum())
    for sc in ["ssp126", "ssp585"]:
        for hz in ["dry", "wet"]:
            d = delta.sel(scen=sc, haz=hz)
            s_le05 = float(den_s.where(d <= 0.5).sum()) / tot_a * 100
            s_gt1 = float(den_s.where(d > 1.0).sum()) / tot_a * 100
            s_gt05 = float(den_s.where(d > 0.5).sum()) / tot_a * 100
            print(f"  {sc} {hz}: Δ<=0.5 = {s_le05:.2f}% | Δ>0.5 = {s_gt05:.2f}% | Δ>1.0 = {s_gt1:.2f}%")
