# -*- coding: utf-8 -*-
"""Generate 100% True Vector PDFs (pdf.fonttype = 42, zero raster images) for:
  1. Supplementary Fig. S11: Fig_S11_ERA5_ASM_Diverging.{pdf,png}  (0.5° CWatM vs. ERA5 SM correlation)
  2. Supplementary Fig. S12: Fig_S12_ISMN_ASM_Diverging.{pdf,png}  (In-situ ISMN site SM correlation)

Uses zero-centered diverging `RdBu` colormaps spanning [-1.0, +1.0] (Reviewer #3 Minor Comment 4)
with pure-vector `PolyCollection` / `PatchCollection` primitives and pure-vector colorbars.
"""
from __future__ import annotations

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection, PatchCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Polygon, Circle

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_ROOT = os.path.dirname(HERE)
DATA_DIR = os.path.join(SCRIPTS_ROOT, "data_tables")
DEFAULT_OUT_DIR = os.path.join(os.path.dirname(SCRIPTS_ROOT), "figures")


def find_world_shapefile() -> str:
    candidates = [
        "/opt/spark/local/lib/python3.11/dist-packages/pyogrio/tests/fixtures/naturalearth_lowres/naturalearth_lowres.shp",
        r"./data/shapefiles/natural_earth/ne_50m_admin_0_countries.shp",
    ]
    for fp in candidates:
        if os.path.exists(fp):
            return fp
    raise FileNotFoundError("Natural Earth country shapefile not found; update path in find_world_shapefile().")


def set_vector_rc() -> None:
    mpl.rcParams.update({
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.unicode_minus": False,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    })


def draw_vector_colorbar_horizontal(
    cax: plt.Axes,
    cmap_name: str,
    vmin: float,
    vmax: float,
    ticks: list[float],
    label: str,
    n_segs: int = 256,
) -> None:
    cmap = plt.get_cmap(cmap_name)
    norm = Normalize(vmin=vmin, vmax=vmax)
    xs = np.linspace(vmin, vmax, n_segs + 1)
    verts = [[(xs[i], 0.0), (xs[i + 1], 0.0), (xs[i + 1], 1.0), (xs[i], 1.0)] for i in range(n_segs)]
    cols = [cmap(norm(0.5 * (xs[i] + xs[i + 1]))) for i in range(n_segs)]

    pc = PolyCollection(verts, facecolors=cols, edgecolors=cols, linewidths=0.35, antialiased=True, rasterized=False, zorder=2)
    cax.add_collection(pc)
    border = Polygon([(vmin, 0.0), (vmax, 0.0), (vmax, 1.0), (vmin, 1.0)], closed=True, facecolor="none", edgecolor="#111111", linewidth=0.85, rasterized=False, zorder=4)
    cax.add_patch(border)

    cax.set_xlim(vmin, vmax)
    cax.set_ylim(0.0, 1.0)
    for sp in cax.spines.values():
        sp.set_visible(False)
    cax.set_yticks([])
    cax.set_xticks(ticks)
    cax.tick_params(axis="x", which="major", labelsize=9.0, length=3.2, width=0.85, pad=2.5, bottom=True, top=False)
    cax.set_xlabel(label, fontsize=10.2, fontweight="bold", labelpad=3.5)


def plot_fig_s11(out_dir: str = DEFAULT_OUT_DIR, mirror_dir: str | None = None) -> None:
    set_vector_rc()
    z = np.load(os.path.join(DATA_DIR, "s11_era5_r_05deg.npz"))
    lats_05, lons_05, r_05 = z["lat"], z["lon"], z["r"]
    valid_05 = np.isfinite(r_05)

    cmap_div = plt.get_cmap("RdBu")
    q_bins = 128
    q_idx = np.zeros_like(r_05, dtype=int)
    q_idx[valid_05] = np.clip(np.round(((r_05[valid_05] + 1.0) * 0.5) * (q_bins - 1)).astype(int), 0, q_bins - 1)

    rects, rect_cols = [], []
    for r_i in range(len(lats_05)):
        if not np.any(valid_05[r_i]):
            continue
        y_lo, y_hi = float(lats_05[r_i] - 0.25), float(lats_05[r_i] + 0.25)
        c_j = 0
        while c_j < len(lons_05):
            if not valid_05[r_i, c_j]:
                c_j += 1
                continue
            qb = q_idx[r_i, c_j]
            c_end = c_j + 1
            while c_end < len(lons_05) and valid_05[r_i, c_end] and q_idx[r_i, c_end] == qb:
                c_end += 1
            x_lo, x_hi = float(lons_05[c_j] - 0.25), float(lons_05[c_end - 1] + 0.25)
            rects.append([(x_lo, y_lo), (x_hi, y_lo), (x_hi, y_hi), (x_lo, y_hi)])
            rect_cols.append(cmap_div(qb / float(q_bins - 1)))
            c_j = c_end

    world = gpd.read_file(find_world_shapefile()).to_crs("EPSG:4326")
    name_col = "name" if "name" in world.columns else ("ADMIN" if "ADMIN" in world.columns else world.columns[0])
    world = world[world[name_col] != "Antarctica"].copy()

    fig = plt.figure(figsize=(6.6, 3.85))
    ax_map = fig.add_axes([0.02, 0.25, 0.96, 0.72])
    ax_map.add_collection(PolyCollection(rects, facecolors=rect_cols, edgecolors=rect_cols, linewidths=0.15, antialiased=True, rasterized=False, zorder=2))
    world.boundary.plot(ax=ax_map, color="#222222", linewidth=0.38, zorder=3)
    ax_map.set_xlim(-180.0, 180.0)
    ax_map.set_ylim(-56.5, 84.5)
    ax_map.set_aspect("equal")
    ax_map.axis("off")

    ax_cb = fig.add_axes([0.22, 0.14, 0.56, 0.044])
    draw_vector_colorbar_horizontal(ax_cb, "RdBu", -1.0, 1.0, [-1.0, -0.5, 0.0, 0.5, 1.0], "Correlation with ERA5 SM (r)")

    for d_out in [out_dir] + ([mirror_dir] if mirror_dir else []):
        os.makedirs(d_out, exist_ok=True)
        fig.savefig(os.path.join(d_out, "Fig_S11_ERA5_ASM_Diverging.pdf"), bbox_inches="tight", pad_inches=0.05, facecolor="white")
        fig.savefig(os.path.join(d_out, "Fig_S11_ERA5_ASM_Diverging.png"), dpi=600, bbox_inches="tight", pad_inches=0.05, facecolor="white")
    plt.close(fig)
    print("[OK] Fig_S11_ERA5_ASM_Diverging.{pdf,png}")


def plot_fig_s12(out_dir: str = DEFAULT_OUT_DIR, mirror_dir: str | None = None) -> None:
    set_vector_rc()
    df_st = pd.read_csv(os.path.join(DATA_DIR, "s12_ismn_site_r.csv"))

    world = gpd.read_file(find_world_shapefile()).to_crs("EPSG:4326")
    name_col = "name" if "name" in world.columns else ("ADMIN" if "ADMIN" in world.columns else world.columns[0])
    world_land = world[world[name_col] != "Antarctica"].dissolve()

    fig = plt.figure(figsize=(6.6, 3.85))
    ax_map = fig.add_axes([0.02, 0.25, 0.96, 0.72])
    world_land.plot(ax=ax_map, color="#EBEBEB", edgecolor="#808080", linewidth=0.55, zorder=1)

    cmap_div = plt.get_cmap("RdBu")
    rad_deg = 2.89
    circles = [Circle((float(r.lon), float(r.lat)), radius=rad_deg) for r in df_st.itertuples()]
    c_colors = []
    for r in df_st.itertuples():
        rgba = list(cmap_div(float(np.clip((r.r_ismn + 1.0) * 0.5, 0.0, 1.0))))
        rgba[3] = 0.88
        c_colors.append(tuple(rgba))

    ax_map.add_collection(PatchCollection(circles, facecolors=c_colors, edgecolors="none", linewidths=0.0, antialiased=True, rasterized=False, zorder=3))
    ax_map.set_xlim(-180.0, 180.0)
    ax_map.set_ylim(-56.5, 84.5)
    ax_map.set_aspect("equal")
    ax_map.axis("off")

    ax_cb = fig.add_axes([0.22, 0.14, 0.56, 0.044])
    draw_vector_colorbar_horizontal(ax_cb, "RdBu", -1.0, 1.0, [-1.0, -0.5, 0.0, 0.5, 1.0], "Correlation with ISMN SM (Site Data, r)")

    for d_out in [out_dir] + ([mirror_dir] if mirror_dir else []):
        os.makedirs(d_out, exist_ok=True)
        fig.savefig(os.path.join(d_out, "Fig_S12_ISMN_ASM_Diverging.pdf"), bbox_inches="tight", pad_inches=0.05, facecolor="white")
        fig.savefig(os.path.join(d_out, "Fig_S12_ISMN_ASM_Diverging.png"), dpi=600, bbox_inches="tight", pad_inches=0.05, facecolor="white")
    plt.close(fig)
    print("[OK] Fig_S12_ISMN_ASM_Diverging.{pdf,png}")


if __name__ == "__main__":
    plot_fig_s11()
    plot_fig_s12()
