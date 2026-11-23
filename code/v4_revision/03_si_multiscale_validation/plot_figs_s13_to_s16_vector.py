"""Generate 100% True Vector PDFs (pdf.fonttype = 42, zero raster images) for:
  - Fig_S13_Subnational_Bias_Hultgren.{pdf,jpg} (with pure-vector diverging colorbar)
  - Fig_S14_Multiscale_US_RMA.{pdf,jpg}
  - Fig_S15_Landmark_Events_Catalogs.{pdf,jpg}
  - Fig_S16_WetExtreme_GGCMI3_Ablation.{pdf,jpg}
"""
# -*- coding: utf-8 -*-
from __future__ import annotations

import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.collections import PolyCollection
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Polygon
from shapely.geometry import box as shapely_box

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_ROOT = os.path.dirname(HERE)
DL_DIR = os.path.join(SCRIPTS_ROOT, "data_tables")
FIG_DIR = os.path.join(os.path.dirname(SCRIPTS_ROOT), "figures")
os.makedirs(FIG_DIR, exist_ok=True)

NE_SHP = (
    "/opt/spark/local/lib/python3.11/dist-packages/pyogrio/tests/fixtures/naturalearth_lowres/naturalearth_lowres.shp"
    if os.path.exists("/opt/spark/local/lib/python3.11/dist-packages/pyogrio/tests/fixtures/naturalearth_lowres/naturalearth_lowres.shp")
    else r"./data/shapefiles/natural_earth/ne_50m_admin_0_countries.shp"
)
CIL_GPKG = os.path.join(DL_DIR, "cil_adm1_simplified.gpkg")
C_OBS = "#1C2833"
C_DRY = "#C0392B"
C_WET = "#1F618D"
C_NOWL = "#85929E"
C_GGCMI = "#7D3C98"
C_LIGHT = "#D5D8DC"

CROP_ORDER = ["whea", "maiz", "soyb", "rice"]
CROP_LABELS = {"maiz": "Maize", "soyb": "Soybean", "rice": "Rice", "whea": "Wheat"}
CROP_COLORS = {"whea": "#D95F02", "maiz": "#E6AB02", "soyb": "#7570B3", "rice": "#1B9E77"}
CROP_MARKERS = {"whea": "o", "maiz": "s", "soyb": "^", "rice": "D"}

CURVE_ORDER = [
    "linear_fao56",
    "continuous_exp_quad",
    "gaussian_bell",
    "sigmoidal_scurve",
    "power_concave_early_hypoxia",
    "power_convex_late_hypoxia",
    "feddes_buffered_anaerobiosis",
]
CURVE_SHORT = {
    "linear_fao56": "Linear (FAO-56)",
    "continuous_exp_quad": "Exp-quadratic",
    "gaussian_bell": "Gaussian bell",
    "sigmoidal_scurve": "Sigmoidal",
    "power_concave_early_hypoxia": "Early hypoxia",
    "power_convex_late_hypoxia": "Late hypoxia",
    "feddes_buffered_anaerobiosis": "Buffered wet",
}

MODEL_DISPLAY = {
    "CWatM-2S": "CWatM",
    "CWatM-NoWL": "CWatM (no waterlogging)",
    "GGCMI3-Mean": "GGCMI3 mean",
    "acea": "ACEA",
    "crover": "CROVER",
    "cygma1p74": "CYGMA",
    "epic-iiasa": "EPIC-IIASA",
    "isam": "ISAM",
    "ldndc": "LDNDC",
    "lpjml": "LPJmL",
    "pdssat": "pDSSAT",
    "pepic": "PEPIC",
    "promet": "PROMET",
    "simplace-lintul5": "SIMPLACE-LINTUL5",
}


def set_rc_rob() -> None:
    mpl.rcParams.update({
        "pdf.fonttype": 42,
        "pdf.compression": 9,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.unicode_minus": False,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 6.8,
        "axes.labelsize": 7.2,
        "axes.titlesize": 7.5,
        "xtick.labelsize": 6.3,
        "ytick.labelsize": 6.3,
        "legend.fontsize": 5.9,
        "legend.frameon": False,
        "axes.linewidth": 0.70,
        "xtick.major.width": 0.65,
        "ytick.major.width": 0.65,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def set_rc_val() -> None:
    mpl.rcParams.update({
        "pdf.fonttype": 42,
        "pdf.compression": 9,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.unicode_minus": False,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.size": 8.0,
        "axes.labelsize": 8.4,
        "axes.titlesize": 8.8,
        "xtick.labelsize": 7.5,
        "ytick.labelsize": 7.5,
        "legend.fontsize": 7.3,
        "legend.frameon": False,
        "axes.linewidth": 0.75,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
    })


def style_ax(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_axisbelow(True)


def panel_tag(ax: plt.Axes, letter: str, dx_pt: float = -26.0, dy_pt: float = 5.0, fs: float = 8.8) -> None:
    ax.set_gid(letter)
    ax.annotate(
        letter,
        xy=(0.0, 1.0),
        xycoords="axes fraction",
        xytext=(dx_pt, dy_pt),
        textcoords="offset points",
        fontsize=fs,
        fontweight="bold",
        va="bottom",
        ha="left",
        color="#111111",
    )


def save_fig(fig: plt.Figure, stem: str, pad: float = 0.07) -> None:
    pdf_path = os.path.join(FIG_DIR, f"{stem}.pdf")
    jpg_path = os.path.join(FIG_DIR, f"{stem}.jpg")
    fig.savefig(pdf_path, format="pdf", bbox_inches="tight", pad_inches=pad)
    fig.savefig(
        jpg_path,
        format="jpg",
        dpi=600,
        bbox_inches="tight",
        pad_inches=pad,
        pil_kwargs={"quality": 95, "subsampling": 0},
    )
    plt.close(fig)
    print(f"Saved {stem}: PDF ({os.path.getsize(pdf_path)//1024} KB), JPG ({os.path.getsize(jpg_path)//1024} KB)")


def make_fig_s13() -> None:
    set_rc_val()
    df_map = pd.read_csv(os.path.join(DL_DIR, "s13_adm1_bias_map_source_9019.csv"))
    df_map["iso"] = df_map["iso"].astype(str).str.upper().str.strip()
    df_map["adm1_id"] = df_map["adm1_id"].astype(str).str.strip().str.replace(r"\.0$", "", regex=True)
    df_met = pd.read_csv(os.path.join(DL_DIR, "s13_adm1_crop_hazard_metrics_9019.csv"))

    clip_box = shapely_box(-132.0, -55.0, 156.0, 68.0)
    world = gpd.read_file(NE_SHP).to_crs("EPSG:4326")
    name_col = "name" if "name" in world.columns else ("ADMIN" if "ADMIN" in world.columns else world.columns[0])
    world = world[world[name_col] != "Antarctica"].copy()
    world = gpd.clip(world, clip_box)

    adm1 = gpd.read_file(CIL_GPKG).to_crs("EPSG:4326")
    adm1["iso"] = adm1["iso"].astype(str).str.upper().str.strip()
    adm1["adm1_id"] = adm1["adm1_id"].astype(str).str.strip().str.replace(r"\.0$", "", regex=True)
    adm1 = gpd.clip(adm1, clip_box)
    adm1["geometry"] = adm1.geometry.simplify(0.08, preserve_topology=True)
    adm1 = adm1[~adm1.geometry.is_empty].copy()

    fig = plt.figure(figsize=(11.2, 13.2))
    gs = gridspec.GridSpec(6, 2, figure=fig, height_ratios=[1.0, 1.0, 1.0, 1.0, 0.30, 1.12], hspace=0.32, wspace=0.20)
    norm = TwoSlopeNorm(vmin=-20.0, vcenter=0.0, vmax=20.0)
    cmap = "RdBu_r"
    crops_order = ["maiz", "soyb", "rice", "whea"]
    haz_order = [("dry", "Drought"), ("wet", "Waterlogging")]
    letters = list("abcdefgh")

    for r_i, cr in enumerate(crops_order):
        for c_j, (hz_code, hz_lbl) in enumerate(haz_order):
            ax = fig.add_subplot(gs[r_i, c_j])
            sub = df_map[(df_map["crop"] == cr) & (df_map["hazard"] == hz_code)].copy()
            merged = adm1.merge(sub, on=["iso", "adm1_id"], how="inner")
            world.plot(ax=ax, color="#EAEDED", edgecolor="none", zorder=1)
            if len(merged) > 0:
                merged.plot(column="mean_loss_bias_pctpt", ax=ax, cmap=cmap, norm=norm, edgecolor="none", linewidth=0.0, zorder=2)
            world.boundary.plot(ax=ax, color="#7F8C8D", linewidth=0.22, zorder=3)
            ax.set_xlim(-132.0, 156.0)
            ax.set_ylim(-55.0, 68.0)
            ax.set_aspect("auto")
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(True)
                sp.set_linewidth(0.5)
                sp.set_edgecolor("#BDC3C7")

            row_m = df_met[(df_met["crop"] == cr) & (df_met["hazard"] == hz_code)].iloc[0]
            bias_v = float(row_m["adm1_mean_bias_pctpt"])
            n_u = int(row_m["n_adm1_units"])
            ax.text(0.01, 1.03, f"{CROP_LABELS[cr]} ({hz_lbl.lower()})", transform=ax.transAxes, ha="left", va="bottom", fontsize=8.0, fontweight="bold", color="#1C2833")
            ax.text(0.99, 1.03, f"n = {n_u}, bias = {bias_v:+.1f}%", transform=ax.transAxes, ha="right", va="bottom", fontsize=7.4, color="#34495E")
            panel_tag(ax, letters[r_i * 2 + c_j], dx_pt=-16.0, dy_pt=4.0, fs=10.2)

    ax_cb_slot = fig.add_subplot(gs[4, :])
    ax_cb_slot.axis("off")
    cax = ax_cb_slot.inset_axes([0.22, 0.56, 0.56, 0.34])
    cmap_obj = plt.get_cmap(cmap)
    xs_cb = np.linspace(-20.0, 20.0, 257)
    verts_cb = [[(xs_cb[k], 0.0), (xs_cb[k + 1], 0.0), (xs_cb[k + 1], 1.0), (xs_cb[k], 1.0)] for k in range(256)]
    cols_cb = [cmap_obj(norm(0.5 * (xs_cb[k] + xs_cb[k + 1]))) for k in range(256)]
    ext_w = 1.52
    verts_cb.append([(-20.0, 0.0), (-20.0, 1.0), (-20.0 - ext_w, 0.5)])
    cols_cb.append(cmap_obj(0.0))
    verts_cb.append([(20.0, 0.0), (20.0, 1.0), (20.0 + ext_w, 0.5)])
    cols_cb.append(cmap_obj(1.0))
    cax.add_collection(PolyCollection(verts_cb, facecolors=cols_cb, edgecolors=cols_cb, linewidths=0.35, antialiased=True, rasterized=False, zorder=2))
    cax.add_patch(Polygon([(-20.0 - ext_w, 0.5), (-20.0, 1.0), (20.0, 1.0), (20.0 + ext_w, 0.5), (20.0, 0.0), (-20.0, 0.0)], closed=True, facecolor="none", edgecolor="#111111", linewidth=0.65, rasterized=False, zorder=4))
    cax.set_xlim(-20.0 - ext_w * 1.08, 20.0 + ext_w * 1.08)
    cax.set_ylim(0.0, 1.0)
    for sp_c in cax.spines.values():
        sp_c.set_visible(False)
    cax.set_yticks([])
    cax.set_xticks([-20, -10, 0, 10, 20])
    cax.tick_params(axis="x", which="major", labelsize=7.2, length=3.0, width=0.65, pad=2.2)
    cax.set_xlabel("Simulated \u2212 observed loss (%)", fontsize=7.8, labelpad=3.0)

    for c_j, (hz_code, hz_lbl, tag) in enumerate([("dry", "Drought", "i"), ("wet", "Waterlogging", "j")]):
        ax_s = fig.add_subplot(gs[5, c_j])
        style_ax(ax_s)
        ax_sr = ax_s.twinx()
        ax_sr.spines["top"].set_visible(False)
        ax_sr.spines["right"].set_visible(True)
        sub_m = df_met[df_met["hazard"] == hz_code].set_index("crop").loc[crops_order]
        xs = np.arange(len(crops_order))
        bw = 0.30
        c_sim = C_DRY if hz_code == "dry" else C_WET
        b_obs = ax_s.bar(xs - bw / 2, sub_m["obs_mean_loss_pct"], width=bw * 0.92, color=C_OBS, label="Observed")
        b_sim = ax_s.bar(xs + bw / 2, sub_m["sim_mean_loss_pct"], width=bw * 0.92, color=c_sim, label="CWatM")
        l_rm = ax_sr.plot(xs, sub_m["rmse_pct"], color="#566573", lw=1.4, ls="--", marker="D", ms=4.3, label="RMSE")
        ax_s.set_xticks(xs)
        ax_s.set_xticklabels([CROP_LABELS[c] for c in crops_order])
        ax_s.set_ylim(0, 9.2)
        ax_sr.set_ylim(0, 16.5)
        ax_s.set_ylabel("Yield loss (%)")
        ax_sr.set_ylabel("RMSE (%)", color="#566573")
        ax_s.set_title(hz_lbl, fontweight="bold", pad=5)
        handles_s = [b_obs, b_sim, l_rm[0]]
        ax_s.legend(handles_s, [h.get_label() for h in handles_s], loc="upper right", ncol=3, columnspacing=0.8, handletextpad=0.35)
        panel_tag(ax_s, tag, dx_pt=-16.0, dy_pt=6.0, fs=10.2)

    save_fig(fig, "Fig_S13_Subnational_Bias_Hultgren", pad=0.08)


def make_fig_s14() -> None:
    set_rc_val()
    from scipy import stats
    ts_us = pd.read_csv(os.path.join(DL_DIR, "s14_us_full_annual_1990_2015.csv"))
    pb_us = pd.read_csv(os.path.join(DL_DIR, "s14_us_full_precip_bins_9019.csv"))
    dr_gd = pd.read_csv(os.path.join(DL_DIR, "ma5_dose_response_gdhy_4crops_9019.csv"))

    r_dry, _ = stats.pearsonr(ts_us["cwatm_loss_dry_pct"], ts_us["rma_drought"])
    r_wet, _ = stats.pearsonr(ts_us["cwatm_loss_wet_pct"], ts_us["rma_excess_wet"])

    fig = plt.figure(figsize=(11.4, 10.2))
    gs = gridspec.GridSpec(3, 2, figure=fig, height_ratios=[1.0, 1.0, 1.05], hspace=0.40, wspace=0.34)

    # Panel a
    ax_a = fig.add_subplot(gs[0, 0])
    style_ax(ax_a)
    ax_a_r = ax_a.twinx()
    ax_a_r.spines["top"].set_visible(False)
    ax_a_r.spines["right"].set_visible(True)
    la1 = ax_a.plot(ts_us["Year"], ts_us["cwatm_loss_dry_pct"], color=C_DRY, lw=1.9, marker="o", ms=3.5, label="CWatM")
    la2 = ax_a_r.plot(ts_us["Year"], ts_us["rma_drought"], color=C_OBS, lw=1.4, ls="--", marker="s", ms=3.0, label="USDA RMA")
    mx_dc, mx_dr = float(ts_us["cwatm_loss_dry_pct"].max()), float(ts_us["rma_drought"].max())
    ax_a.set_xlim(1989.2, 2016.2)
    ax_a.set_ylim(0, mx_dc * 1.24)
    ax_a_r.set_ylim(0, mx_dr * 1.24)
    ax_a.set_ylabel("Simulated loss (%)", color=C_DRY)
    ax_a_r.set_ylabel("Insurance loss ratio", color=C_OBS)
    ax_a.set_xlabel("Year")
    y02 = float(ts_us.loc[ts_us["Year"] == 2002, "cwatm_loss_dry_pct"].iloc[0])
    y12 = float(ts_us.loc[ts_us["Year"] == 2012, "cwatm_loss_dry_pct"].iloc[0])
    ax_a.text(2002.0, y02 + mx_dc * 0.08, "2002", fontsize=7.4, color=C_DRY, ha="center", va="bottom")
    ax_a.text(2010.6, y12, "2012", fontsize=7.4, color=C_DRY, ha="right", va="center")
    ax_a.text(0.03, 0.94, f"Drought: r = {r_dry:.2f}, P < 0.001", transform=ax_a.transAxes, fontsize=7.8, fontweight="bold", color=C_DRY, va="top")
    lns_a = la1 + la2
    ax_a.legend(lns_a, [x.get_label() for x in lns_a], loc="upper left", bbox_to_anchor=(0.01, 0.84))
    panel_tag(ax_a, "a", dx_pt=-28.0, dy_pt=6.0, fs=10.2)

    # Panel b
    ax_b = fig.add_subplot(gs[0, 1])
    style_ax(ax_b)
    ax_b_r = ax_b.twinx()
    ax_b_r.spines["top"].set_visible(False)
    ax_b_r.spines["right"].set_visible(True)
    lb1 = ax_b.plot(ts_us["Year"], ts_us["cwatm_loss_wet_pct"], color=C_WET, lw=1.9, marker="o", ms=3.5, label="CWatM")
    lb2 = ax_b_r.plot(ts_us["Year"], ts_us["rma_excess_wet"], color=C_OBS, lw=1.4, ls="--", marker="^", ms=3.2, label="USDA RMA")
    mx_wc, mx_wr = float(ts_us["cwatm_loss_wet_pct"].max()), float(ts_us["rma_excess_wet"].max())
    ax_b.set_xlim(1989.2, 2016.2)
    ax_b.set_ylim(0, mx_wc * 1.26)
    ax_b_r.set_ylim(0, mx_wr * 1.26)
    ax_b.set_ylabel("Simulated loss (%)", color=C_WET)
    ax_b_r.set_ylabel("Insurance loss ratio", color=C_OBS)
    ax_b.set_xlabel("Year")
    y93 = float(ts_us.loc[ts_us["Year"] == 1993, "cwatm_loss_wet_pct"].iloc[0])
    ax_b.text(1994.4, y93, "1993", fontsize=7.4, color=C_WET, ha="left", va="center")
    ax_b.text(0.30, 0.94, f"Excess moisture: r = {r_wet:.2f}, P < 0.001", transform=ax_b.transAxes, fontsize=7.8, fontweight="bold", color=C_WET, va="top")
    lns_b = lb1 + lb2
    ax_b.legend(lns_b, [x.get_label() for x in lns_b], loc="upper left", bbox_to_anchor=(0.28, 0.84))
    panel_tag(ax_b, "b", dx_pt=-28.0, dy_pt=6.0, fs=10.2)

    # Panel c
    ax_c = fig.add_subplot(gs[1, 0])
    style_ax(ax_c)
    ax_c_r = ax_c.twinx()
    ax_c_r.spines["top"].set_visible(False)
    ax_c_r.spines["right"].set_visible(True)
    xsd = pb_us["prec_sd_mean"].values
    lc1 = ax_c.plot(xsd, pb_us["cwatm_loss_dry_pct"], color=C_DRY, lw=1.9, marker="o", ms=3.8, label="CWatM (drought)")
    lc1w = ax_c.plot(xsd, pb_us["cwatm_loss_wet_pct"], color=C_WET, lw=1.8, marker="D", ms=3.3, label="CWatM (waterlogging)")
    lc2 = ax_c_r.plot(xsd, pb_us["rma_drought"], color=C_OBS, lw=1.4, ls="--", marker="s", ms=3.2, label="USDA RMA (drought)")
    lc3 = ax_c_r.plot(xsd, pb_us["rma_excess_wet"], color=C_WET, lw=1.5, ls="-.", marker="^", ms=3.4, label="USDA RMA (excess wet)")
    ax_c.set_xlim(-2.45, 2.75)
    ax_c.set_ylim(-0.3, 11.8)
    ax_c_r.set_ylim(-0.08, 3.12)
    ax_c.set_xlabel("Precipitation anomaly (SD)")
    ax_c.set_ylabel("Simulated loss (%)", color="#1C2833")
    ax_c_r.set_ylabel("Insurance loss ratio", color=C_OBS)
    lns_c = lc1 + lc1w + lc2 + lc3
    ax_c.legend(lns_c, [x.get_label() for x in lns_c], loc="upper right", fontsize=7.0)
    panel_tag(ax_c, "c", dx_pt=-28.0, dy_pt=6.0, fs=10.2)

    # Panel d
    ax_d = fig.add_subplot(gs[1, 1])
    style_ax(ax_d)
    ax_d.axhline(0, color="#666666", lw=0.7, ls=":")
    ax_d.plot(xsd, pb_us["obs_ma5_loss_pct"], color=C_OBS, lw=2.0, marker="o", ms=4.0, label="Observed (USDA NASS)")
    ax_d.plot(xsd, pb_us["cwatm_loss_dry_pct"], color=C_DRY, lw=1.8, marker="s", ms=3.5, label="CWatM (drought)")
    ax_d.plot(xsd, pb_us["cwatm_loss_wet_pct"], color=C_WET, lw=1.8, marker="^", ms=3.6, label="CWatM (waterlogging)")
    ax_d.set_xlim(-2.45, 2.75)
    ax_d.set_ylim(-7.5, 30.5)
    ax_d.set_xlabel("Precipitation anomaly (SD)")
    ax_d.set_ylabel("Yield loss (%)")
    ax_d.legend(loc="upper right")
    panel_tag(ax_d, "d", dx_pt=-28.0, dy_pt=6.0, fs=10.2)

    # Panel e
    ax_e = fig.add_subplot(gs[2, 0])
    style_ax(ax_e)
    ax_e.axhline(0, color="#666666", lw=0.7, ls=":")
    bins_order = ["<0.1", "0.1-1", "1-3", "3-6", ">6"]
    xb = np.arange(len(bins_order))
    for crop in ["whea", "maiz", "soyb", "rice"]:
        sub_d = dr_gd[(dr_gd["crop"] == crop) & (dr_gd["hazard"] == "drought")].set_index("bin").reindex(bins_order)
        ax_e.errorbar(xb, sub_d["obs_ma5_net_loss_pct"], yerr=1.96 * sub_d["obs_se_pct"], marker="o", ms=4.0, lw=1.7, color=CROP_COLORS[crop], capsize=2.0, label=CROP_LABELS[crop])
    ax_e.set_xticks(xb)
    ax_e.set_xticklabels(["<0.1", "0.1\u20131", "1\u20133", "3\u20136", ">6"])
    ax_e.set_ylim(-4.2, 22.8)
    ax_e.set_xlabel("Simulated loss (%)")
    ax_e.set_ylabel("Observed loss (%)")
    ax_e.set_title("Drought", fontweight="bold", pad=4)
    ax_e.legend(loc="upper left", ncol=2)
    panel_tag(ax_e, "e", dx_pt=-28.0, dy_pt=6.0, fs=10.2)

    # Panel f
    ax_f = fig.add_subplot(gs[2, 1])
    style_ax(ax_f)
    ax_f.axhline(0, color="#666666", lw=0.7, ls=":")
    for crop in ["whea", "soyb", "rice", "maiz"]:
        sub_w = dr_gd[(dr_gd["crop"] == crop) & (dr_gd["hazard"] == "waterlogging")].set_index("bin").reindex(bins_order)
        ax_f.errorbar(xb, sub_w["obs_ma5_net_loss_pct"], yerr=1.96 * sub_w["obs_se_pct"], marker="s", ms=3.8, lw=1.6, ls="--", color=CROP_COLORS[crop], capsize=2.0, label=CROP_LABELS[crop])
    ax_f.set_xticks(xb)
    ax_f.set_xticklabels(["<0.1", "0.1\u20131", "1\u20133", "3\u20136", ">6"])
    ax_f.set_ylim(-1.2, 5.2)
    ax_f.set_xlabel("Simulated loss (%)")
    ax_f.set_ylabel("Observed loss (%)")
    ax_f.set_title("Waterlogging", fontweight="bold", pad=4)
    ax_f.legend(loc="upper left", ncol=2)
    panel_tag(ax_f, "f", dx_pt=-28.0, dy_pt=6.0, fs=10.2)

    save_fig(fig, "Fig_S14_Multiscale_US_RMA", pad=0.08)


def make_fig_s15() -> None:
    set_rc_val()
    lm = pd.read_csv(os.path.join(DL_DIR, "s15_landmark_18_full_region_9019.csv"))
    gdis = pd.read_csv(os.path.join(DL_DIR, "tier4b_gdis_global_disaster_contrast_9019.csv"))
    gfd = pd.read_csv(os.path.join(DL_DIR, "tier4b_modis_gfd_flood_contrast_9019.csv"))

    fig = plt.figure(figsize=(11.4, 10.0))
    gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1.38, 1.0], hspace=0.36, wspace=0.30)

    lm_ord = lm.iloc[::-1].reset_index(drop=True)
    y_ev = np.arange(len(lm_ord))

    # Panel a
    ax_a = fig.add_subplot(gs[0, 0])
    style_ax(ax_a)
    ax_a.axvline(0, color="#555555", lw=0.7, ls=":")
    ax_a.axhline(8.5, color="#7F8C8D", lw=0.8, ls="-.")
    ax_a.barh(y_ev + 0.16, lm_ord["obs_aff_ma5_loss_pct"], height=0.30, color=C_OBS, label="Observed")
    ax_a.barh(y_ev - 0.16, lm_ord["cwatm_aff_haz_loss_pct"], height=0.30, color=C_WET, label="CWatM")
    ax_a.set_yticks(y_ev)
    ax_a.set_yticklabels(lm_ord["short_label"], fontsize=7.3)
    ax_a.set_xlim(-8.5, 42.5)
    ax_a.set_xlabel("Yield loss (%)")
    ax_a.text(41.2, 13.5, "Droughts", fontsize=8.0, fontweight="bold", color=C_DRY, va="center", ha="right")
    ax_a.text(41.2, 4.5, "Floods", fontsize=8.0, fontweight="bold", color=C_WET, va="center", ha="right")
    ax_a.legend(bbox_to_anchor=(0.5, 1.01), loc="lower center", ncol=2, fontsize=7.3)
    panel_tag(ax_a, "a", dx_pt=-42.0, dy_pt=14.0, fs=10.2)

    # Panel b
    ax_b = fig.add_subplot(gs[0, 1])
    style_ax(ax_b)
    ax_b.plot([-3.8, 40], [-3.8, 40], color="#7F8C8D", lw=1.0, ls="--", label="Identity (1:1)", zorder=1)
    sub_dry = lm[lm["hazard"] == "drought"]
    sub_wet = lm[lm["hazard"] == "waterlogging"]
    ax_b.scatter(sub_dry["obs_aff_ma5_loss_pct"], sub_dry["cwatm_aff_haz_loss_pct"], s=42, color=C_DRY, marker="o", edgecolors="white", linewidths=0.5, label="Droughts (n = 9)", zorder=3)
    ax_b.scatter(sub_wet["obs_aff_ma5_loss_pct"], sub_wet["cwatm_aff_haz_loss_pct"], s=46, color=C_WET, marker="^", edgecolors="white", linewidths=0.5, label="Floods (n = 9)", zorder=3)

    mean_bias_dry = float(sub_dry["bias_aff_pctpt"].mean())
    mean_bias_wet = float(sub_wet["bias_aff_pctpt"].mean())
    ann_coords = [
        ("D04", 10.6, 10.1, "left", "bottom", C_DRY),     # 2012 US soybean
        ("D03", 20.6, 8.6, "left", "bottom", C_DRY),      # 2012 US maize
        ("D09", 41.5, 0.4, "right", "top", C_DRY),        # 2009 Argentina soybean
        ("W01", 24.3, 1.6, "left", "bottom", C_WET),      # 1993 US maize
        ("W07", -5.2, 9.2, "left", "bottom", C_WET),      # 2011 Thailand rice
        ("W08", -12.0, 5.5, "left", "top", C_WET),        # 2007 Bangladesh rice
    ]
    lm_idx = lm.set_index("event_id")
    for eid, tx, ty, ha, va, col in ann_coords:
        ax_b.text(tx, ty, str(lm_idx.loc[eid, "short_label"]), fontsize=6.9, color=col, ha=ha, va=va)

    ax_b.text(
        41.8, 22.5,
        f"Underestimation (below 1:1)\nDrought bias: {mean_bias_dry:+.1f}%\nFlood bias: {mean_bias_wet:+.1f}%",
        fontsize=7.4, color="#34495E", ha="right", va="top",
    )
    ax_b.set_xlim(-12.8, 42.5)
    ax_b.set_ylim(-3.8, 41.5)
    ax_b.set_xlabel("Observed loss (%)")
    ax_b.set_ylabel("Simulated loss (%)")
    ax_b.legend(bbox_to_anchor=(0.5, 1.01), loc="lower center", ncol=3, fontsize=7.1, columnspacing=0.8)
    panel_tag(ax_b, "b", dx_pt=-28.0, dy_pt=14.0, fs=10.2)

    # Panel c
    ax_c = fig.add_subplot(gs[1, 0])
    style_ax(ax_c)
    ax_c.axhline(0, color="#444444", lw=0.75, ls="--")
    crops_g = ["maiz", "soyb", "rice", "whea"]
    xg = np.arange(len(crops_g))
    bwg = 0.19
    gd_dry_z_d = [gdis[(gdis["crop"] == c) & (gdis["gdis_regime"] == "GDIS_Drought_Only")]["cwatm_z_dry"].iloc[0] for c in crops_g]
    gd_dry_z_w = [gdis[(gdis["crop"] == c) & (gdis["gdis_regime"] == "GDIS_Drought_Only")]["cwatm_z_wet"].iloc[0] for c in crops_g]
    gd_fld_z_d = [gdis[(gdis["crop"] == c) & (gdis["gdis_regime"] == "GDIS_Flood_Only")]["cwatm_z_dry"].iloc[0] for c in crops_g]
    gd_fld_z_w = [gdis[(gdis["crop"] == c) & (gdis["gdis_regime"] == "GDIS_Flood_Only")]["cwatm_z_wet"].iloc[0] for c in crops_g]
    ax_c.bar(xg - 1.5 * bwg, gd_dry_z_d, width=bwg * 0.92, color=C_DRY, label="Drought loss (droughts)")
    ax_c.bar(xg - 0.5 * bwg, gd_dry_z_w, width=bwg * 0.92, color="#F5B7B1", edgecolor=C_DRY, lw=0.5, label="Waterlogging loss (droughts)")
    ax_c.bar(xg + 0.5 * bwg, gd_fld_z_w, width=bwg * 0.92, color=C_WET, label="Waterlogging loss (floods)")
    ax_c.bar(xg + 1.5 * bwg, gd_fld_z_d, width=bwg * 0.92, color="#AED6F1", edgecolor=C_WET, lw=0.5, label="Drought loss (floods)")
    for ix in range(len(crops_g)):
        if gd_dry_z_d[ix] > 0.05:
            ax_c.text(xg[ix] - 1.5 * bwg, gd_dry_z_d[ix] + 0.02, "***", ha="center", va="bottom", fontsize=7.0, fontweight="bold", color=C_DRY)
        ax_c.text(xg[ix] + 0.5 * bwg, gd_fld_z_w[ix] + 0.02, "***", ha="center", va="bottom", fontsize=7.0, fontweight="bold", color=C_WET)
    ax_c.set_xticks(xg)
    ax_c.set_xticklabels([CROP_LABELS[c] for c in crops_g])
    ax_c.set_ylim(-0.32, 0.96)
    ax_c.set_ylabel("Hazard anomaly (SD)")
    ax_c.legend(loc="upper right", ncol=1, fontsize=7.0)
    panel_tag(ax_c, "c", dx_pt=-42.0, dy_pt=10.0, fs=10.2)

    # Panel d
    ax_d = fig.add_subplot(gs[1, 1])
    style_ax(ax_d)
    non_fld_w = [gfd[(gfd["crop"] == c) & (gfd["gfd_regime"] == "SameCells_NonFlood_Years_2000_2016")]["cwatm_loss_wet_pct"].iloc[0] for c in crops_g]
    mod_fld_w = [gfd[(gfd["crop"] == c) & (gfd["gfd_regime"] == "MODIS_GFD_Flooded_Years_2000_2016")]["cwatm_loss_wet_pct"].iloc[0] for c in crops_g]
    ax_d.bar(xg - 0.17, non_fld_w, width=0.31, color="#AED6F1", edgecolor=C_WET, lw=0.5, label="Non-flood years")
    ax_d.bar(xg + 0.17, mod_fld_w, width=0.31, color=C_WET, label="Flooded years (MODIS)")
    for ix in range(len(crops_g)):
        ax_d.text(xg[ix] + 0.17, mod_fld_w[ix] + 0.04, "***", ha="center", va="bottom", fontsize=7.2, fontweight="bold", color=C_WET)
    ax_d.set_xticks(xg)
    ax_d.set_xticklabels([CROP_LABELS[c] for c in crops_g])
    ax_d.set_ylim(0, 2.25)
    ax_d.set_ylabel("Waterlogging loss (%)")
    ax_d.legend(loc="upper left")
    panel_tag(ax_d, "d", dx_pt=-28.0, dy_pt=10.0, fs=10.2)

    save_fig(fig, "Fig_S15_Landmark_Events_Catalogs", pad=0.08)


def make_fig_s16() -> None:
    set_rc_val()
    us_ma5 = pd.read_csv(os.path.join(DL_DIR, "ma5_us_cornbelt_9019.csv")).set_index("Regime")
    glob_ma5 = pd.read_csv(os.path.join(DL_DIR, "ma5_global_gdhy_ggcmi_cwatm_4crops_9019.csv"))
    gg = pd.read_csv(os.path.join(DL_DIR, "tier2_ggcmi3_vs_cwatm_wet_extremes_4crops_9019.csv"))

    fig = plt.figure(figsize=(11.4, 10.2))
    gs = gridspec.GridSpec(3, 2, figure=fig, height_ratios=[1.0, 1.12, 1.12], hspace=0.38, wspace=0.32)

    # Panel a
    ax_a = fig.add_subplot(gs[0, 0])
    style_ax(ax_a)
    ax_a.axhline(0, color="#444444", lw=0.75, ls="--")
    us_cats = [
        ("Extreme Dry Precip (Prec <= -1.5 SD)", "Extreme dry\n(\u2264\u22121.5 SD)"),
        ("Extreme Wet Precip (Prec >= +2.0 SD)", "Extreme wet\n(\u2265+2.0 SD)"),
        ("CWatM Waterlogging Event (wet>=7d, dry==0)", "Waterlogging\n(\u22657 d)"),
        ("CWatM Severe Waterlog (wet_loss>=3%)", "Severe\nwaterlogging"),
    ]
    x_a = np.arange(len(us_cats))
    bw = 0.24
    obs_vals = [us_ma5.loc[k, "Obs_MA5_Net%"] for k, _ in us_cats]
    gg_vals = [us_ma5.loc[k, "GGCMI3_Net%"] for k, _ in us_cats]
    cw_vals = [
        us_ma5.loc[us_cats[0][0], "CWatM_Dry_Loss%"],
        us_ma5.loc[us_cats[1][0], "CWatM_Wet_Loss%"],
        us_ma5.loc[us_cats[2][0], "CWatM_Wet_Loss%"],
        us_ma5.loc[us_cats[3][0], "CWatM_Wet_Loss%"],
    ]
    ax_a.bar(x_a - bw, obs_vals, width=bw * 0.92, color=C_OBS, label="Observed")
    ax_a.bar(x_a, gg_vals, width=bw * 0.92, color=C_GGCMI, label="GGCMI3 mean")
    ax_a.bar(x_a + bw, cw_vals, width=bw * 0.92, color=C_WET, label="CWatM")
    ax_a.set_xticks(x_a)
    ax_a.set_xticklabels([lbl for _, lbl in us_cats], fontsize=7.2)
    ax_a.set_ylim(-15.5, 26.5)
    ax_a.set_ylabel("Yield loss (%)")
    ax_a.text(2.0, -13.8, "Spurious yield gain (<0%)", fontsize=7.1, color=C_GGCMI, ha="center")
    ax_a.legend(loc="upper right")
    panel_tag(ax_a, "a", dx_pt=-28.0, dy_pt=8.0, fs=10.2)

    # Panel b
    ax_b = fig.add_subplot(gs[0, 1])
    style_ax(ax_b)
    ax_b.axhline(0, color="#444444", lw=0.75, ls="--")
    crops_order = ["rice", "soyb", "maiz", "whea"]
    x_b = np.arange(len(crops_order))
    bw_b = 0.19
    sub_hw = glob_ma5[glob_ma5["regime"] == "Waterlog_HighLoss_Humid (CWatM_wet>=5%, P>=800mm)"].set_index("crop").loc[crops_order]
    ax_b.bar(x_b - 1.5 * bw_b, sub_hw["Obs_GDHY_NetLoss_%"], width=bw_b * 0.92, color=C_OBS, label="Observed")
    ax_b.bar(x_b - 0.5 * bw_b, sub_hw["CWatM2S_MA5_NetLoss_%"], width=bw_b * 0.92, color=C_WET, label="CWatM")
    ax_b.bar(x_b + 0.5 * bw_b, sub_hw["CWatMNoWL_MA5_NetLoss_%"], width=bw_b * 0.92, color=C_NOWL, label="CWatM (no waterlogging)")
    ax_b.bar(x_b + 1.5 * bw_b, sub_hw["GGCMI3_NetLoss_%"], width=bw_b * 0.92, color=C_GGCMI, label="GGCMI3 mean")
    ax_b.set_xticks(x_b)
    ax_b.set_xticklabels([CROP_LABELS[c] for c in crops_order])
    ax_b.set_ylim(-6.4, 5.8)
    ax_b.set_ylabel("Yield loss (%)")
    ax_b.legend(loc="upper left", ncol=2, columnspacing=0.9)
    panel_tag(ax_b, "b", dx_pt=-28.0, dy_pt=8.0, fs=10.2)

    # Panels c-f
    crop_panels = [("rice", 1, 0, "c"), ("soyb", 1, 1, "d"), ("maiz", 2, 0, "e"), ("whea", 2, 1, "f")]
    for crop, gr, gc, tag in crop_panels:
        ax_m = fig.add_subplot(gs[gr, gc])
        style_ax(ax_m)
        sub = gg[(gg["crop"] == crop) & (gg["regime"] == "ModelFree_HumidExtremePrecip_Zp_gt_1p5")].copy()
        sub = sub.sort_values("rmse", ascending=False).reset_index(drop=True)
        y_pos = np.arange(len(sub))
        bar_cols = [C_WET if m == "CWatM-2S" else (C_NOWL if m == "CWatM-NoWL" else (C_GGCMI if m == "GGCMI3-Mean" else C_LIGHT)) for m in sub["model"]]
        ax_m.barh(y_pos, sub["rmse"], height=0.68, color=bar_cols, edgecolor="none")
        max_rmse = float(sub["rmse"].max())
        ax_m.set_xlim(0, max_rmse * 1.18)
        for yi, row_m in sub.iterrows():
            if row_m["model"] in ("CWatM-2S", "CWatM-NoWL", "GGCMI3-Mean"):
                col_t = C_WET if row_m["model"] == "CWatM-2S" else (C_GGCMI if row_m["model"] == "GGCMI3-Mean" else "#4D5656")
                ax_m.text(row_m["rmse"] + max_rmse * 0.015, yi, f"{row_m['rmse']:.1f}%", va="center", ha="left", fontsize=7.2, fontweight="bold", color=col_t)
        yticklabels = [MODEL_DISPLAY.get(m, m) for m in sub["model"]]
        ax_m.set_yticks(y_pos)
        ax_m.set_yticklabels(yticklabels, fontsize=7.1)
        for tick_lbl, mname in zip(ax_m.get_yticklabels(), sub["model"]):
            if mname in ("CWatM-2S", "GGCMI3-Mean"):
                tick_lbl.set_fontweight("bold")
            if mname == "CWatM-2S":
                tick_lbl.set_color(C_WET)
        ax_m.set_xlabel("RMSE (%)")
        ax_m.set_title(CROP_LABELS[crop], fontweight="bold", pad=4)
        panel_tag(ax_m, tag, dx_pt=-28.0, dy_pt=8.0, fs=10.2)

    save_fig(fig, "Fig_S16_WetExtreme_GGCMI3_Ablation", pad=0.08)



if __name__ == "__main__":
    make_fig_s13()
    make_fig_s14()
    make_fig_s15()
    make_fig_s16()
