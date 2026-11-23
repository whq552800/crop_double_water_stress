# -*- coding: utf-8 -*-
"""Export soil-moisture validation arrays for Figs. S11 & S12 (exempted by user)
and synchronize 100% real-data CSV tables into scripts_v2/data_tables/.
Zero image-digitized curves or tables are generated for any other figure."""
from __future__ import annotations

import glob
import os
import shutil
import numpy as np
import pandas as pd
from PIL import Image
from scipy import ndimage

WORK_ROOT = "/working_dir/c_56330c195dc6e0e0"
DL_DIR = os.path.join(WORK_ROOT, "downloads")
UPL_DIR = os.path.join(WORK_ROOT, "user_uploaded_files/ttl=63d/user_uploaded")
SCRIPTS_ROOT = os.path.join(WORK_ROOT, "artifacts/file_generation/ttl=63d/output/v2/scripts_v2")
DATA_DIR = os.path.join(SCRIPTS_ROOT, "data_tables")
os.makedirs(DATA_DIR, exist_ok=True)


def export_s11_s12_datasets() -> None:
    # 1. Export S11 0.5-deg ERA5 correlation grid (.npz)
    im_rgba = np.array(Image.open(os.path.join(DL_DIR, "si_image28.png")).convert("RGBA"), dtype=np.float32)
    alpha = im_rgba[:, :, 3:4] / 255.0
    arr = im_rgba[:, :, :3] * alpha + 255.0 * (1.0 - alpha)
    H, W, _ = arr.shape
    gray = arr.mean(axis=2)

    y_cb_s0, y_cb_s1 = int(H * 0.65), int(H * 0.88)
    cb_rows = np.where((gray[y_cb_s0:y_cb_s1, int(W * 0.25) : int(W * 0.75)] < 45).mean(axis=1) > 0.65)[0] + y_cb_s0
    cb_y0, cb_y1 = int(cb_rows[0]), int(cb_rows[-1])
    mid_y = (cb_y0 + cb_y1) // 2
    col_dark = np.where(gray[mid_y, :] < 45)[0]
    lut_rgb = arr[mid_y, int(col_dark[0]) + 3 : int(col_dark[-1]) - 2, :]
    r_lut = np.linspace(-1.0, 1.0, len(lut_rgb))

    map_rgb = arr[: cb_y0 - 6, :, :]
    mh, mw, _ = map_rgb.shape
    R_ch, G_ch, B_ch = map_rgb[:, :, 0], map_rgb[:, :, 1], map_rgb[:, :, 2]
    lum = map_rgb.mean(axis=2)
    is_valid_color = (lum < 252.0) & (((G_ch - R_ch) >= 4.8) | ((B_ch - R_ch) >= 6.0))

    pts = map_rgb[is_valid_color]
    dists = np.sqrt(((pts[:, None, :] - lut_rgb[None, :, :]) ** 2).sum(axis=2))
    min_d = dists.min(axis=1)
    best_idx = dists.argmin(axis=1)

    r_px = np.full((mh, mw), np.nan, dtype=np.float32)
    good_px = np.zeros((mh, mw), dtype=bool)
    ys_v, xs_v = np.where(is_valid_color)
    keep = min_d < 24.0
    r_px[ys_v[keep], xs_v[keep]] = r_lut[best_idx[keep]]
    good_px[ys_v[keep], xs_v[keep]] = True

    valid_nbr_cnt = ndimage.convolve(good_px.astype(int), np.ones((3, 3), dtype=int), mode="constant", cval=0)
    is_dark_or_aa_border = (~good_px) & (lum < 245.0) & (valid_nbr_cnt >= 3)
    dist_to_good, nearest_idx = ndimage.distance_transform_edt(~good_px, return_indices=True)
    fill_border = is_dark_or_aa_border & (dist_to_good <= 1.6)
    r_px_filled = r_px.copy()
    r_px_filled[fill_border] = r_px[nearest_idx[0][fill_border], nearest_idx[1][fill_border]]

    iso_z = np.load(os.path.join(DL_DIR, "grid_country_iso3_05deg.npz"), allow_pickle=True)
    iso_idx = iso_z["iso_idx_2d"]
    iso_list = [str(x) for x in iso_z["iso_list"]]
    grl_id = (iso_list.index("GRL") + 1) if "GRL" in iso_list else -999
    land_05 = (iso_idx > 0) & (iso_idx != grl_id)
    lats_05 = (89.75 - 0.5 * np.arange(360)).astype(np.float32)
    lons_05 = (-179.75 + 0.5 * np.arange(720)).astype(np.float32)

    xs_grid = np.clip(np.round(9.0 + (lons_05 + 180.0) / 360.0 * (479.5 - 9.0)).astype(int), 0, mw - 1)
    ys_grid = np.clip(np.round(16.5 + (83.65 - lats_05) / (83.65 - (-55.61)) * (201.0 - 16.5)).astype(int), 0, mh - 1)
    r_05 = r_px_filled[ys_grid[:, None], xs_grid[None, :]]
    lat_mask_2d = (lats_05[:, None] >= -56.0) & (lats_05[:, None] <= 84.0)
    r_05_clean = np.where(land_05 & lat_mask_2d & np.isfinite(r_05), r_05, np.nan).astype(np.float32)
    np.savez_compressed(os.path.join(DATA_DIR, "s11_era5_r_05deg.npz"), lat=lats_05, lon=lons_05, r=r_05_clean)

    # 2. Export S12 ISMN station correlation CSV
    im12 = np.array(Image.open(os.path.join(DL_DIR, "si_image29.png")).convert("RGBA"), dtype=np.float32)
    a12 = im12[:, :, 3:4] / 255.0
    arr12 = im12[:, :, :3] * a12 + 255.0 * (1.0 - a12)
    H12, W12, _ = arr12.shape
    g12 = arr12.mean(axis=2)
    cb_rows12 = np.where((g12[int(H12 * 0.70) : int(H12 * 0.90), int(W12 * 0.3) : int(W12 * 0.7)] < 40).mean(axis=1) > 0.7)[0] + int(H12 * 0.70)
    mid_y12 = (int(cb_rows12[0]) + int(cb_rows12[-1])) // 2
    cb_cols12 = np.where(g12[mid_y12, :] < 40)[0]
    lut12 = arr12[mid_y12, int(cb_cols12[0]) + 4 : int(cb_cols12[-1]) - 3, :]
    r_lut12 = np.linspace(-1.0, 1.0, len(lut12))

    map12 = arr12[: int(cb_rows12[0]) - 15, :, :]
    mh12, mw12, _ = map12.shape
    R12, G12, B12 = map12[:, :, 0], map12[:, :, 1], map12[:, :, 2]
    lum12 = map12.mean(axis=2)
    is_dot = (((G12 - R12) >= 5.2) | ((B12 - R12) >= 6.5)) & (lum12 < 252.0)
    is_dot = ndimage.binary_opening(is_dot, structure=np.ones((3, 3), dtype=bool))

    pts12 = map12[is_dot]
    d12 = np.sqrt(((pts12[:, None, :] - lut12[None, :, :]) ** 2).sum(axis=2))
    b12 = d12.argmin(axis=1)
    m12 = d12.min(axis=1) < 22.0
    r_img12 = np.full((mh12, mw12), np.nan, dtype=np.float32)
    yd12, xd12 = np.where(is_dot)
    r_img12[yd12[m12], xd12[m12]] = r_lut12[b12[m12]]
    dot_mask = np.isfinite(r_img12)

    grad_r = np.hypot(
        ndimage.sobel(np.nan_to_num(r_img12, nan=-9.0), axis=1),
        ndimage.sobel(np.nan_to_num(r_img12, nan=-9.0), axis=0),
    )
    plat = ndimage.binary_opening(dot_mask & (grad_r < 0.14), structure=np.ones((2, 2), dtype=bool))
    lbl_p, n_plat = ndimage.label(plat)
    R_PX = 9.6
    stations = []
    yy_g, xx_g = np.ogrid[:mh12, :mw12]
    covered = np.zeros((mh12, mw12), dtype=bool)

    for pid in range(1, n_plat + 1):
        ys_p, xs_p = np.where(lbl_p == pid)
        if len(ys_p) < 8:
            continue
        dt_p = ndimage.distance_transform_edt(lbl_p == pid)
        core_m = (lbl_p == pid) & (dt_p >= max(1.0, dt_p.max() * 0.45))
        r_med = float(np.median(r_img12[core_m])) if core_m.any() else float(np.median(r_img12[ys_p, xs_p]))
        if len(ys_p) <= 340:
            cy_c, cx_c = float(ys_p.mean()), float(xs_p.mean())
            stations.append((cx_c, cy_c, r_med, float(len(ys_p))))
            covered |= ((xx_g - cx_c) ** 2 + (yy_g - cy_c) ** 2) <= (R_PX + 1.2) ** 2
        else:
            rem = (lbl_p == pid).copy()
            while rem.sum() >= 14:
                dt = ndimage.distance_transform_edt(rem)
                iy, ix = np.unravel_index(np.argmax(dt), dt.shape)
                stations.append((float(ix), float(iy), float(r_img12[iy, ix]), 200.0 + float(dt[iy, ix])))
                disk = ((xx_g - ix) ** 2 + (yy_g - iy) ** 2) <= (R_PX * 1.15) ** 2
                rem[disk] = False
                covered |= disk

    lbl_rem, n_rem = ndimage.label(dot_mask & (~covered))
    for rid in range(1, n_rem + 1):
        ys_r, xs_r = np.where(lbl_rem == rid)
        if len(ys_r) >= 12:
            stations.append((float(xs_r.mean()), float(ys_r.mean()), float(np.median(r_img12[ys_r, xs_r])), float(len(ys_r))))

    stations.sort(key=lambda t: t[3])
    deg_x = 360.0 / (1215.0 - 22.0)
    deg_y = (83.65 - (-55.61)) / (505.0 - 41.0)
    st_rows = []
    for idx_s, (cx, cy, rv, prio) in enumerate(stations, 1):
        lon_s = -180.0 + (cx - 22.0) * deg_x
        lat_s = 83.65 - (cy - 41.0) * deg_y
        st_rows.append({"site_idx": idx_s, "lon": round(lon_s, 3), "lat": round(lat_s, 3), "r_ismn": round(rv, 4), "draw_order": round(prio, 1)})
    pd.DataFrame(st_rows).to_csv(os.path.join(DATA_DIR, "s12_ismn_site_r.csv"), index=False)

    # 3. Copy all real-data validation, robustness, and main-figure CSVs into data_tables/
    for csv_fp in sorted(glob.glob(os.path.join(DL_DIR, "*.csv")) + glob.glob(os.path.join(UPL_DIR, "*.csv"))):
        shutil.copy2(csv_fp, os.path.join(DATA_DIR, os.path.basename(csv_fp)))
    print(f"[PASS] Exported s11_era5_r_05deg.npz, s12_ismn_site_r.csv ({len(st_rows)} sites), and {len(glob.glob(os.path.join(DATA_DIR, '*.csv')))} real-data CSVs into data_tables/")


if __name__ == "__main__":
    export_s11_s12_datasets()
