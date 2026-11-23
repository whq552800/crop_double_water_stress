"""Generate 100% True Vector PDFs (pdf.fonttype = 42, zero raster images) for:
  - Fig_S9_Parameter_StressCurve_Sensitivity.{pdf,jpg}
  - Fig_S17_Methodological_Robustness_Synthesis.{pdf,jpg}
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


def make_fig_s9() -> None:
    set_rc_rob()
    df_prof = pd.read_csv(os.path.join(DL_DIR, "r3_ks_curve_profiles.csv"))
    df_ks = pd.read_csv(os.path.join(DL_DIR, "r3_ks_functional_forms_summary.csv"))
    df_kcky = pd.read_csv(os.path.join(DL_DIR, "r4_kc_ky_sensitivity_summary.csv"))
    df_rice = pd.read_csv(os.path.join(DL_DIR, "r4_rice_stage_waterlogging_sensitivity.csv"))

    fig = plt.figure(figsize=(7.20, 7.60))
    gs = gridspec.GridSpec(3, 2, figure=fig, wspace=0.46, hspace=0.38)
    axes = [fig.add_subplot(gs[r, c]) for r in range(3) for c in range(2)]
    for ax, ltr in zip(axes, ["a", "b", "c", "d", "e", "f"]):
        style_ax(ax)
        panel_tag(ax, ltr)

    # Panel a
    ax = axes[0]
    x = df_prof["theta_norm"].values
    ax.axvline(0.35, color="#BDC3C7", lw=0.75, ls="--", zorder=1)
    ax.axvline(0.70, color="#BDC3C7", lw=0.75, ls="--", zorder=1)
    if "gaussian_bell_low" in df_prof.columns:
        ax.fill_between(x, df_prof["gaussian_bell_low"].values, df_prof["gaussian_bell_high"].values, color="#D95F02", alpha=0.13, zorder=2)
    if "power_concave_early_hypoxia" in df_prof.columns and "power_convex_late_hypoxia" in df_prof.columns:
        ax.fill_between(
            x,
            np.minimum(df_prof["power_concave_early_hypoxia"].values, df_prof["power_convex_late_hypoxia"].values),
            np.maximum(df_prof["power_concave_early_hypoxia"].values, df_prof["power_convex_late_hypoxia"].values),
            color="#1B9E77", alpha=0.10, zorder=2,
        )
    ax.plot(x, df_prof["continuous_exp_quad"], color="#2E86C1", lw=1.25, ls="--", label="Exp-quadratic", zorder=4)
    ax.plot(x, df_prof["gaussian_bell"], color="#D95F02", lw=1.30, ls="-.", label="Gaussian (\u00b1\u03c3)", zorder=4)
    ax.plot(x, df_prof["sigmoidal_scurve"], color="#7D3C98", lw=1.25, ls=":", label="Sigmoidal", zorder=4)
    ax.plot(x, df_prof["power_concave_early_hypoxia"], color="#1B9E77", lw=1.10, ls="--", label="Power curves", zorder=3)
    ax.plot(x, df_prof["power_convex_late_hypoxia"], color="#1B9E77", lw=1.10, ls="--", zorder=3)
    ax.plot(x, df_prof["feddes_buffered_anaerobiosis"], color="#7F8C8D", lw=1.15, ls="-", label="Buffered wet", zorder=3)
    ax.plot(x, df_prof["linear_fao56"], color=C_OBS, lw=2.1, ls="-", label="FAO-56 linear", zorder=6)
    ax.text(0.17, 1.055, "Dry", fontsize=6.8, color="#666666", ha="center", va="bottom")
    ax.text(0.525, 1.055, "Optimal", fontsize=6.8, color="#444444", ha="center", va="bottom")
    ax.text(0.85, 1.055, "Wet", fontsize=6.8, color="#666666", ha="center", va="bottom")
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.0, 1.15)
    ax.set_xlabel("Normalized soil moisture")
    ax.set_ylabel("Moisture stress factor")
    ax.legend(loc="lower center", bbox_to_anchor=(0.525, 0.02), fontsize=5.8, handlelength=1.25, handletextpad=0.35, labelspacing=0.18)

    # Panel b
    ax = axes[1]
    ax.axvline(0.0, color="#7F8C8D", lw=0.85, ls="--", zorder=1)
    y_pos = np.arange(len(CURVE_ORDER))[::-1]
    offsets = {"whea": 0.22, "maiz": 0.08, "soyb": -0.08, "rice": -0.22}
    for crop in CROP_ORDER:
        sub_c = df_ks[df_ks["crop"] == crop].set_index("curve_id").loc[CURVE_ORDER]
        diff_c = sub_c["gap_central_pct_pt"].values
        err_lo = np.maximum(diff_c - sub_c["gap_param_min_pct_pt"].values, 0.0)
        err_hi = np.maximum(sub_c["gap_param_max_pct_pt"].values - diff_c, 0.0)
        ax.errorbar(
            diff_c, y_pos + offsets[crop], xerr=[err_lo, err_hi],
            fmt=CROP_MARKERS[crop], color=CROP_COLORS[crop], markersize=4.3, lw=0.9, capsize=1.8,
            label=CROP_LABELS[crop], zorder=4,
        )
    ax.set_xlim(-0.72, 1.44)
    ax.set_xticks([-0.5, 0.0, 0.5, 1.0])
    ax.set_yticks(y_pos)
    ax.set_yticklabels([CURVE_SHORT[c] for c in CURVE_ORDER])
    ax.set_xlabel("Dry\u2013wet loss gap (% pts)")
    ax.set_ylabel("Stress curve form")
    ax.legend(loc="lower right", fontsize=6.2, handletextpad=0.3, borderpad=0.2)

    # Panel c (WITHOUT synthetic OOS ΔR² = -0.04!)
    ax = axes[2]
    max_dr = 0.0
    min_grid_r = 1.0
    for crop in CROP_ORDER:
        sub_c = df_ks[df_ks["crop"] == crop].set_index("curve_id").loc[CURVE_ORDER]
        r_vals = sub_c["faostat_national_anom_r"].values
        max_dr = max(max_dr, float(np.max(np.abs(sub_c["delta_faostat_r_vs_linear"].values))))
        if "rel_yield_corr_with_linear" in sub_c.columns:
            min_grid_r = min(min_grid_r, float(sub_c["rel_yield_corr_with_linear"].min()))
        ax.plot(r_vals, y_pos, marker=CROP_MARKERS[crop], color=CROP_COLORS[crop], lw=1.25, markersize=4.2, label=CROP_LABELS[crop])
    gauss_spread = float(df_ks[df_ks["curve_id"] == "gaussian_bell"][["dry_param_spread_pct", "wet_param_spread_pct"]].to_numpy().max())
    ax.text(
        0.98, 0.08,
        f"Max |\u0394r| \u2264 {max_dr:.2f}\nGrid r = {min_grid_r:.2f}\u20131.00\n\u03c3 spread \u2264 {gauss_spread:.2f}% pts",
        transform=ax.transAxes, fontsize=6.0, va="bottom", ha="right", color="#333333",
    )
    ax.set_xlim(0.0, 0.68)
    ax.set_xticks([0.0, 0.2, 0.4, 0.6])
    ax.set_yticks(y_pos)
    ax.set_yticklabels([CURVE_SHORT[c] for c in CURVE_ORDER])
    ax.set_xlabel("Yield correlation (r)")
    ax.set_ylabel("Stress curve form")

    # Panel d
    ax = axes[3]
    ax.axhline(0.0, color="#BDC3C7", lw=0.75, ls="--", zorder=1)
    ax.axvline(1.0, color="#BDC3C7", lw=0.75, ls=":", zorder=1)
    for crop in CROP_ORDER:
        sub_kc = df_kcky[(df_kcky["parameter"] == "Kc") & (df_kcky["crop"] == crop)].sort_values("multiplier")
        sub_ky = df_kcky[(df_kcky["parameter"] == "Ky") & (df_kcky["crop"] == crop)].sort_values("multiplier")
        mk, col = CROP_MARKERS[crop], CROP_COLORS[crop]
        ax.plot(sub_kc["multiplier"], sub_kc["rel_yield_change_pct"], color=col, lw=1.40, ls="-", marker=mk, markersize=3.8, markerfacecolor=col)
        ax.plot(sub_ky["multiplier"], sub_ky["rel_yield_change_pct"], color=col, lw=1.20, ls="--", marker=mk, markersize=3.8, markerfacecolor="white", markeredgewidth=0.95)
    h_kc = mpl.lines.Line2D([], [], color=C_OBS, lw=1.40, ls="-", marker="o", markersize=3.8, markerfacecolor=C_OBS, label="Crop factor Kc")
    h_ky = mpl.lines.Line2D([], [], color=C_OBS, lw=1.20, ls="--", marker="o", markersize=3.8, markerfacecolor="white", markeredgewidth=0.95, label="Yield factor Ky")
    ax.legend(handles=[h_kc, h_ky], loc="upper right", fontsize=6.3)
    ax.set_xticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3])
    ax.set_xlabel("Parameter multiplier")
    ax.set_ylabel("Yield change (%)")

    # Panel e
    ax = axes[4]
    ax.axhline(1.0, color="#7F8C8D", lw=0.85, ls="--", zorder=1)
    ax.axvline(1.0, color="#BDC3C7", lw=0.75, ls=":", zorder=1)
    for crop in CROP_ORDER:
        sub_kc = df_kcky[(df_kcky["parameter"] == "Kc") & (df_kcky["crop"] == crop)].sort_values("multiplier")
        sub_ky = df_kcky[(df_kcky["parameter"] == "Ky") & (df_kcky["crop"] == crop)].sort_values("multiplier")
        mk, col = CROP_MARKERS[crop], CROP_COLORS[crop]
        ax.plot(sub_kc["multiplier"], np.clip(sub_kc["dry_to_wet_ratio"].values, 1e-3, None), color=col, lw=1.40, ls="-", marker=mk, markersize=3.8, markerfacecolor=col, label=CROP_LABELS[crop])
        ax.plot(sub_ky["multiplier"], np.clip(sub_ky["dry_to_wet_ratio"].values, 1e-3, None), color=col, lw=1.20, ls="--", marker=mk, markersize=3.8, markerfacecolor="white", markeredgewidth=0.95)
    ax.set_yscale("log", nonpositive="clip")
    ax.set_ylim(0.22, 10.5)
    ax.set_yticks([0.2, 0.5, 1.0, 2.0, 5.0, 10.0])
    ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_xlim(0.67, 1.48)
    ax.set_xticks([0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3])
    ax.set_xlabel("Parameter multiplier")
    ax.set_ylabel("Dry-to-wet loss ratio")
    ax.legend(loc="center right", bbox_to_anchor=(0.99, 0.55), ncol=2, fontsize=5.8, handlelength=1.1, handletextpad=0.3, columnspacing=0.55, labelspacing=0.20)

    # Panel f
    ax = axes[5]
    sc_order = ["S0_all4_stages", "S1_veg50pct_atten", "S2_repro_late_only", "S3_repro_stage3_only"]
    sc_labels = ["All 4 stages", "50% veg. weight", "Mid+late stages", "Heading only"]
    yr = np.arange(len(sc_order))[::-1]
    comb = df_rice[df_rice["regime"].str.startswith("Combined")].set_index("scheme_id").loc[sc_order]
    pad = df_rice[df_rice["regime"].str.startswith("Irr")].set_index("scheme_id").loc[sc_order]
    bh = 0.28
    ax.barh(yr + bh / 2, comb["wet_loss_pct"], height=bh, color=C_WET, alpha=0.88, label="Global wet")
    ax.barh(yr - bh / 2, comb["dry_loss_pct"], height=bh, color=C_DRY, alpha=0.88, label="Global dry")
    ax.scatter(pad["wet_loss_pct"], yr + bh / 2, marker="D", s=20, facecolors="white", edgecolors="#0B3C5D", linewidths=1.1, zorder=5, label="Paddy wet")
    ax.scatter(pad["dry_loss_pct"], yr - bh / 2, marker="D", s=20, facecolors="white", edgecolors="#922B21", linewidths=1.1, zorder=5, label="Paddy dry")
    ax.set_xlim(0.0, 5.2)
    ax.set_xticks([0, 1, 2, 3, 4, 5])
    ax.set_yticks(yr)
    ax.set_yticklabels(sc_labels)
    ax.set_xlabel("Rice yield loss (%)")
    ax.set_ylabel("Rice stage scheme")
    ax.legend(loc="lower right", fontsize=6.0, handletextpad=0.30, labelspacing=0.18)

    save_fig(fig, "Fig_S9_Parameter_StressCurve_Sensitivity", pad=0.06)


def make_fig_s17() -> None:
    set_rc_rob()
    df_r1 = pd.read_csv(os.path.join(DL_DIR, "r1_counterfactual_baseline_fixed_vs_moving.csv"))
    df_r2 = pd.read_csv(os.path.join(DL_DIR, "r2_spei_twogate_sensitivity.csv"))
    df_r5 = pd.read_csv(os.path.join(DL_DIR, "r5_aggregation_decomposition_compound.csv"))
    df_r7 = pd.read_csv(os.path.join(DL_DIR, "r7_masks_weights_windows_duration.csv"))
    df_r8 = pd.read_csv(os.path.join(DL_DIR, "r8_gcm_anova_hotmodels.csv"))
    df_idr = pd.read_csv(os.path.join(DL_DIR, "r8_idr_cutoff_sweep.csv"))
    df_idrp = pd.read_csv(os.path.join(DL_DIR, "r8_idr_temporal_persistence.csv"))

    fig = plt.figure(figsize=(7.20, 7.80))
    gs = gridspec.GridSpec(3, 2, figure=fig, wspace=0.51, hspace=0.38)
    axes = [fig.add_subplot(gs[r, c]) for r in range(3) for c in range(2)]
    for ax, ltr in zip(axes, ["a", "b", "c", "d", "e", "f"]):
        style_ax(ax)
        panel_tag(ax, ltr)

    # Panel a
    ax = axes[0]
    crops_a = ["whea", "maiz", "soyb", "rice"]
    x_a = np.arange(len(crops_a))
    schemes = [("fixed", "Fixed hist.", 1.00, -0.25), ("moving_clim_nonhaz", "Moving clim.", 0.58, 0.00), ("moving_ma5", "5-yr moving", 0.28, 0.25)]
    w_bar = 0.105
    for sc_key, _, alpha_v, dx in schemes:
        sub = df_r1[(df_r1["baseline_scheme"] == sc_key) & (df_r1["scenario"] == "ssp585_2070_2100")].set_index("crop").loc[crops_a]
        ax.bar(x_a + dx - w_bar / 2, sub["dry_mean_pct"].values, width=w_bar, color=C_DRY, alpha=alpha_v, edgecolor=C_DRY, linewidth=0.5, yerr=sub["dry_gcm_sd_pct"].values, error_kw={"lw": 0.65, "ecolor": "#444444", "capsize": 1.4})
        ax.bar(x_a + dx + w_bar / 2, sub["wet_mean_pct"].values, width=w_bar, color=C_WET, alpha=alpha_v, edgecolor=C_WET, linewidth=0.5, yerr=sub["wet_gcm_sd_pct"].values, error_kw={"lw": 0.65, "ecolor": "#444444", "capsize": 1.4})
    leg_a = [
        mpl.patches.Patch(facecolor=C_DRY, edgecolor=C_DRY, label="Dry (fixed)"),
        mpl.patches.Patch(facecolor=C_WET, edgecolor=C_WET, label="Wet (fixed)"),
        mpl.patches.Patch(facecolor=C_DRY, edgecolor=C_DRY, alpha=0.58, label="Moving clim."),
        mpl.patches.Patch(facecolor=C_DRY, edgecolor=C_DRY, alpha=0.28, label="5-yr moving"),
    ]
    ax.set_ylim(0.0, 8.6)
    ax.legend(handles=leg_a, loc="upper left", bbox_to_anchor=(0.01, 0.99), ncol=2, fontsize=6.0, handlelength=1.1, columnspacing=0.8, labelspacing=0.22)
    ax.set_xticks(x_a)
    ax.set_xticklabels([CROP_LABELS[c] for c in crops_a])
    ax.set_xlabel("Crop")
    ax.set_ylabel("Projected loss (%)")

    # Panel b
    ax = axes[1]
    sub5 = df_r5[df_r5["scenario"] == "ssp585"].copy()
    row_order_b = [("whea", "dry", "Wheat dry"), ("maiz", "dry", "Maize dry"), ("rice", "dry", "Rice dry"), ("soyb", "dry", "Soybean dry"), ("rice", "wet", "Rice wet"), ("maiz", "wet", "Maize wet"), ("soyb", "wet", "Soybean wet"), ("whea", "wet", "Wheat wet")]
    yb = np.arange(len(row_order_b))[::-1]
    sf_vals, si_vals, net_vals = [], [], []
    for c, h, _ in row_order_b:
        r = sub5[(sub5["crop"] == c) & (sub5["hazard"] == h)].iloc[0]
        sf_vals.append(float(r["shapley_freq_extent_pct_pt"]))
        si_vals.append(float(r["shapley_cond_severity_pct_pt"]))
        net_vals.append(float(r["delta_loss_total_pct_pt"]))
    sf_arr, si_arr = np.asarray(sf_vals), np.asarray(si_vals)
    ax.axvline(0.0, color="#7F8C8D", lw=0.8, ls="--", zorder=1)
    pos_left = np.where((si_arr >= 0) & (sf_arr >= 0), si_arr, 0.0)
    neg_left = np.where((si_arr < 0) & (sf_arr < 0), si_arr, 0.0)
    sf_left = np.where(sf_arr >= 0, pos_left, neg_left)
    ax.barh(yb, si_arr, height=0.56, color="#2C3E50", alpha=0.90, label="Severity (\u0394I)")
    ax.barh(yb, sf_arr, left=sf_left, height=0.56, color="#E67E22", alpha=0.82, label="Exposure (\u0394F)")
    ax.scatter(net_vals, yb, color="#111111", marker="D", s=15, zorder=5, label="Net change")
    ax.set_xlim(-0.65, 5.45)
    ax.set_yticks(yb)
    ax.set_yticklabels([lbl for _, _, lbl in row_order_b])
    ax.set_xlabel("Loss change (% pts)")
    ax.set_ylabel("Crop and hazard")
    ax.legend(loc="lower right", fontsize=6.3, handlelength=1.2)

    # Panel c
    ax = axes[2]
    ax.axhline(1.0, color="#7F8C8D", lw=0.8, ls="--", zorder=1)
    upland_crops = ["whea", "maiz", "soyb"]

    def _agg_r2(grp_name: str, xcol: str, ycol: str, crops_sub: list[str] | None = None) -> pd.DataFrame:
        sg = df_r2[df_r2["group"] == grp_name].copy()
        if crops_sub is not None:
            sg = sg[sg["crop"].isin(crops_sub)]
        sg["abs_z"] = sg[xcol].abs()
        out_g = sg.groupby("abs_z")[ycol].agg(["mean", "min", "max"]).reset_index().sort_values("abs_z")
        for col_y in ["mean", "min", "max"]:
            out_g[col_y] = np.clip(out_g[col_y].values, 1e-3, None)
        return out_g

    g_d2 = _agg_r2("dry_sweep_twogate", "spei_dry_cutoff", "dry_rel_to_baseline")
    g_w2 = _agg_r2("wet_sweep_twogate", "spei_wet_cutoff", "wet_rel_to_baseline")
    g_du = _agg_r2("dry_sweep_speionly", "spei_dry_cutoff", "dry_rel_to_baseline", upland_crops)
    g_dr = _agg_r2("dry_sweep_speionly", "spei_dry_cutoff", "dry_rel_to_baseline", ["rice"])
    g_ws = _agg_r2("wet_sweep_speionly", "spei_wet_cutoff", "wet_rel_to_baseline")
    ax.fill_between(g_d2["abs_z"], g_d2["min"], g_d2["max"], color=C_DRY, alpha=0.14, zorder=2)
    ax.fill_between(g_w2["abs_z"], g_w2["min"], g_w2["max"], color=C_WET, alpha=0.14, zorder=2)
    ax.fill_between(g_du["abs_z"], g_du["min"], g_du["max"], color=C_DRY, alpha=0.07, zorder=1)
    ax.plot(g_dr["abs_z"], g_dr["mean"], color=CROP_COLORS["rice"], ls="-.", lw=1.45, marker="D", markersize=3.8, markerfacecolor="white", markeredgewidth=1.0, label="Rice dry (SPEI-only)", zorder=5)
    ax.plot(g_ws["abs_z"], g_ws["mean"], color=C_WET, ls="--", lw=1.15, marker="s", markersize=3.6, markerfacecolor="white", label="Wet (SPEI-only)", zorder=4)
    ax.plot(g_du["abs_z"], g_du["mean"], color=C_DRY, ls="--", lw=1.15, marker="o", markersize=3.6, markerfacecolor="white", label="Upland dry (SPEI-only)", zorder=4)
    ax.plot(g_w2["abs_z"], g_w2["mean"], color=C_WET, ls="-", lw=1.55, marker="s", markersize=3.8, markerfacecolor=C_WET, label="Wet (two-gate)", zorder=6)
    ax.plot(g_d2["abs_z"], g_d2["mean"], color=C_DRY, ls="-", lw=1.55, marker="o", markersize=3.8, markerfacecolor=C_DRY, label="Dry (two-gate)", zorder=6)
    ax.scatter([1.0], [1.0], s=36, facecolors="none", edgecolors="#111111", linewidths=1.0, zorder=7)
    ax.set_yscale("log", nonpositive="clip")
    ax.set_xlim(0.42, 2.12)
    ax.set_ylim(0.38, 30.0)
    ax.set_yticks([0.5, 1.0, 2.0, 5.0, 12.0])
    ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_xlabel("SPEI threshold |z|")
    ax.set_ylabel("Relative loss ratio")
    ax.legend(loc="upper right", fontsize=5.8, handlelength=1.4, labelspacing=0.16)

    # Panel d
    ax = axes[3]
    ax.axvline(0.0, color="#7F8C8D", lw=0.85, ls="--", zorder=1)
    spec_rows = [
        ("r7", "SPAM2020V2r2 Harvested area", "2020 crop map (harv.)"),
        ("r7", "SPAM2020V2r2 Physical area", "2020 crop map (phys.)"),
        ("r7", "Baseline (MIRCA2000, A*Yp, >=1d, 1990-2019)", "2000 crop area map"),
        ("r7", "Pure harvested-area weight (Yp = 1)", "Area-only weight"),
        ("r5_comp", "compound_apportioned", "Compound split"),
        ("r7", "Mode-specific Irr stress mask", "Regime stress mask"),
        ("r7", "Duration >= 3 days", "Duration ≥ 3 d"),
        ("r7", "Duration >= 7 days (1 week)", "Duration \u2265 7 d"),
        ("r7", "Duration >= 14 days (2 weeks)", "Duration \u2265 14 d"),
    ]
    yd_pos = np.arange(len(spec_rows))[::-1]
    offsets_d = {"whea": 0.16, "maiz": 0.06, "soyb": -0.06, "rice": -0.16}
    for crop in CROP_ORDER:
        vals_d = []
        for src_t, key_s, _ in spec_rows:
            if src_t == "r7":
                rr = df_r7[(df_r7["specification"] == key_s) & (df_r7["crop"] == crop)].iloc[0]
                vals_d.append(float(rr["hist_dry_pct"] - rr["hist_wet_pct"]))
            else:
                rr_d = df_r5[(df_r5["crop"] == crop) & (df_r5["scenario"] == "ssp585") & (df_r5["hazard"] == "dry")].iloc[0]
                rr_w = df_r5[(df_r5["crop"] == crop) & (df_r5["scenario"] == "ssp585") & (df_r5["hazard"] == "wet")].iloc[0]
                vals_d.append(float(rr_d["hist_apportioned_loss_pct"] - rr_w["hist_apportioned_loss_pct"]))
        ax.plot(vals_d, yd_pos + offsets_d[crop], marker=CROP_MARKERS[crop], color=CROP_COLORS[crop], lw=0.90, markersize=3.9, label=CROP_LABELS[crop])
    ax.set_xlim(-0.98, 1.54)
    ax.set_xticks([-0.5, 0.0, 0.5, 1.0])
    ax.set_yticks(yd_pos)
    ax.set_yticklabels([lbl for _, _, lbl in spec_rows])
    ax.set_xlabel("Dry\u2013wet loss gap (% pts)")
    ax.set_ylabel("Model specification")
    ax.legend(loc="lower right", fontsize=6.1, handletextpad=0.3, borderpad=0.2)

    # Panel e
    ax = axes[4]
    row_order_e = [("rice", "wet", "Rice wet"), ("whea", "dry", "Wheat dry"), ("maiz", "dry", "Maize dry"), ("rice", "dry", "Rice dry"), ("soyb", "dry", "Soybean dry"), ("maiz", "wet", "Maize wet"), ("soyb", "wet", "Soybean wet"), ("whea", "wet", "Wheat wet")]
    ye = np.arange(len(row_order_e))[::-1]
    for idx_e, (c, h, _) in enumerate(row_order_e):
        rr = df_r8[(df_r8["crop"] == c) & (df_r8["hazard"] == h)].iloc[0]
        m10, m8, lmin, lmax = float(rr["ssp585_10gcm_mean_pct"]), float(rr["ssp585_8gcm_nohot_mean_pct"]), float(rr["ssp585_loo_min_pct"]), float(rr["ssp585_loo_max_pct"])
        col = C_DRY if h == "dry" else C_WET
        ax.errorbar(m10, ye[idx_e] + 0.16, xerr=[[max(m10 - lmin, 0.0)], [max(lmax - m10, 0.0)]], fmt="o", color=col, markersize=4.1, lw=1.1, capsize=1.6, zorder=4)
        ax.scatter([m8], [ye[idx_e] - 0.16], facecolors="white", edgecolors=col, marker="D", s=19, linewidths=1.05, zorder=5)
    h_10 = mpl.lines.Line2D([], [], color=C_OBS, marker="o", markersize=4.1, lw=1.1, label="10-GCM (LOO)")
    h_8 = mpl.lines.Line2D([], [], color="none", marker="D", markersize=4.3, markerfacecolor="white", markeredgecolor=C_OBS, markeredgewidth=1.05, label="8-GCM (no hot)")
    max_hot_diff = float(np.max(np.abs(df_r8["ssp585_10gcm_mean_pct"] - df_r8["ssp585_8gcm_nohot_mean_pct"])))
    ax.text(0.97, 0.06, f"Max |\u0394| = {max_hot_diff:.2f}% pts", transform=ax.transAxes, fontsize=6.3, ha="right", va="bottom", color="#333333")
    ax.set_xlim(-0.1, 6.75)
    ax.set_yticks(ye)
    ax.set_yticklabels([lbl for _, _, lbl in row_order_e])
    ax.set_xlabel("Projected loss (%)")
    ax.set_ylabel("Crop and hazard")
    ax.legend(handles=[h_10, h_8], loc="lower right", bbox_to_anchor=(0.99, 0.17), fontsize=6.2)

    # Panel f
    ax = axes[5]
    sub_idr = df_idr[(df_idr["scenario"] == "ssp585") & (df_idr["hazard"] == "dry") & (np.isclose(df_idr["delta_loss_threshold_pct"], 0.0))].sort_values("idr_threshold")
    ax.axvline(0.40, color="#7F8C8D", lw=0.85, ls="--", zorder=1)
    ax.plot(sub_idr["idr_threshold"], sub_idr["q2_prod_share_pct"], color=C_WET, lw=1.55, marker="s", markersize=3.8, label="Exporter prod. share")
    ax.plot(sub_idr["idr_threshold"], sub_idr["q1_country_share_pct"], color=C_DRY, lw=1.55, marker="o", markersize=3.8, label="High-IDR countries")
    ax.plot(sub_idr["idr_threshold"], sub_idr["q1_prod_share_pct"], color=C_DRY, lw=1.25, ls="--", marker="^", markersize=3.8, label="High-IDR prod. share")
    rp = df_idrp.iloc[0]
    ax.text(0.98, 0.81, f"2014\u201318 vs. 2019\u201323:\nr = {rp['pearson_r_2014_2018_vs_2019_2023']:.2f} ({rp['classification_agreement_idr04_pct']:.1f}% match)", transform=ax.transAxes, fontsize=6.2, color="#333333", va="top", ha="right")
    ax.set_ylim(0.0, 108.0)
    ax.set_xlabel("Import dependency threshold")
    ax.set_ylabel("Global share (%)")
    ax.legend(loc="center right", bbox_to_anchor=(0.99, 0.52), fontsize=6.1, labelspacing=0.22)

    save_fig(fig, "Fig_S17_Methodological_Robustness_Synthesis", pad=0.06)



if __name__ == "__main__":
    make_fig_s9()
    make_fig_s17()
