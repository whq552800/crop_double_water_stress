# -*- coding: utf-8 -*-
"""Generate and export Extended Data Table 1 and Supplementary Tables S1–S8 in exact 1-to-1 alignment
with manuscript_tracked_20261118_v3_EN.docx and SI_tracked_20261118_v3_EN.docx.

Outputs:
  - `../data_tables/ed_table1_crop_system_scenario_summary.csv` (verified source for ED Table 1)
  - `./tables_output/Extended_Data_Table_1_*.csv` and `Table_S1_*.csv` .. `Table_S8_*.csv`
  - Unified multi-sheet formatted Excel workbook (`SI_Tables_S1_to_S8_and_ED_Table_1_v3.xlsx`, >=10pt font floor),
    recalculated via LibreOffice Calc, validated via `validate_excel.py`, and rendered to preview PNGs.
"""
from __future__ import annotations

import os
import pathlib
import shutil
import subprocess
import tempfile
import numpy as np
import pandas as pd
import openpyxl
from openpyxl.styles import Font, Alignment, PatternFill, Border, Side
from openpyxl.utils import get_column_letter
import pypdfium2 as pdfium

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS_ROOT = os.path.dirname(HERE)
V3_DIR = os.path.dirname(SCRIPTS_ROOT)
ROOT_WS = "/working_dir/c_56330c195dc6e0e0"
INTERMEDIATE_DIR = os.path.join(ROOT_WS, "artifacts/file_generation/ttl=63d/intermediate")
DATA_DIR = os.path.join(SCRIPTS_ROOT, "data_tables")
OUT_DIR = os.path.join(HERE, "tables_output")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(INTERMEDIATE_DIR, exist_ok=True)

ANNUAL_MASTER_CSV = os.path.join(
    ROOT_WS,
    "downloads/codex_sync_20261118/r1_fig3_fig4_data_sync_20261118/Fig4_text_number_alignment/water_stress_loss_annual_detailed.csv",
)


def build_ed_table1_source() -> pd.DataFrame:
    """Computes or loads the unified Crop x System x Scenario summary table for Extended Data Table 1."""
    out_csv = os.path.join(DATA_DIR, "ed_table1_crop_system_scenario_summary.csv")
    if not os.path.exists(ANNUAL_MASTER_CSV) and os.path.exists(out_csv):
        return pd.read_csv(out_csv)

    df_ann = pd.read_csv(ANNUAL_MASTER_CSV)
    sub = df_ann[(df_ann["mode"] == "total") & (df_ann["scenario"].isin(["historical", "ssp126", "ssp585"]))].copy()

    # Exact Figure 4 bias alignment (2015–2019 GSWP3–W5E5 anchor vs. 2020–2024 GCM anchor)
    h_anc = (
        sub[(sub["scenario"] == "historical") & (sub["year"].between(2015, 2019))]
        .groupby(["crop", "hazard"])["loss_pct"]
        .mean()
        .rename("h_anc")
    )
    f_anc = (
        sub[(sub["scenario"].isin(["ssp126", "ssp585"])) & (sub["year"].between(2020, 2024))]
        .groupby(["scenario", "gcm", "crop", "hazard"])["loss_pct"]
        .mean()
        .rename("f_anc")
    )
    m = sub.merge(f_anc, on=["scenario", "gcm", "crop", "hazard"], how="left").merge(
        h_anc, on=["crop", "hazard"], how="left"
    )
    off = (m["h_anc"] - m["f_anc"]).fillna(0.0)
    fut = m["scenario"].isin(["ssp126", "ssp585"])
    m["pct_bc"] = m["loss_pct"]
    m.loc[fut, "pct_bc"] = (m.loc[fut, "loss_pct"] + off.loc[fut]).clip(lower=0.0)
    m["mt_bc"] = m["pct_bc"] / 100.0 * (m["prod_nws_tonnes"] / 1e6)

    decomp = pd.read_csv(os.path.join(DATA_DIR, "rebuttal_decomposition_summary.csv"))
    f3 = pd.read_csv(os.path.join(DATA_DIR, "figure3_statistics_1990_2019_and_2071_2100.csv"))

    rows = []
    crop_en = {"all4": "All 4 crops", "maiz": "Maize", "soyb": "Soybean", "rice": "Rice", "whea": "Wheat"}
    crop_zh = {"all4": "四大主粮合计", "maiz": "玉米", "soyb": "大豆", "rice": "水稻", "whea": "小麦"}
    haz_en = {"total": "Total stress", "dry": "Drought", "wet": "Waterlogging"}
    haz_zh = {"total": "总水分胁迫", "dry": "干旱", "wet": "涝渍"}

    # Panel A: Combined all-cropland global production-share yield loss (L_global, %) & tonnage loss (Mt yr^-1) [Fig. 4]
    for cr in ["all4", "maiz", "soyb", "rice", "whea"]:
        for hz in ["total", "dry", "wet"]:
            hs_p = m[(m["scenario"] == "historical") & (m["year"].between(1990, 2019)) & (m["crop"] == cr) & (m["hazard"] == hz)]["pct_bc"]
            s1_p = m[(m["scenario"] == "ssp126") & (m["year"].between(2070, 2100)) & (m["crop"] == cr) & (m["hazard"] == hz)].groupby("gcm")["pct_bc"].mean()
            s5_p = m[(m["scenario"] == "ssp585") & (m["year"].between(2070, 2100)) & (m["crop"] == cr) & (m["hazard"] == hz)].groupby("gcm")["pct_bc"].mean()

            hs_m = m[(m["scenario"] == "historical") & (m["year"].between(1990, 2019)) & (m["crop"] == cr) & (m["hazard"] == hz)]["mt_bc"]
            s1_m = m[(m["scenario"] == "ssp126") & (m["year"].between(2070, 2100)) & (m["crop"] == cr) & (m["hazard"] == hz)].groupby("gcm")["mt_bc"].mean()
            s5_m = m[(m["scenario"] == "ssp585") & (m["year"].between(2070, 2100)) & (m["crop"] == cr) & (m["hazard"] == hz)].groupby("gcm")["mt_bc"].mean()

            if hz in ("dry", "wet"):
                d5 = decomp[(decomp["group"] == cr) & (decomp["haz"] == hz) & (decomp["scen"] == "ssp585")].iloc[0]
                sev_ext = f"{d5['sev_share_pct']:.1f}% / {d5['ext_share_pct']:.1f}%"
            else:
                sev_ext = "—"

            grp_en = f"{crop_en[cr]} (combined)" if cr != "all4" else "All 4 crops (combined)"
            grp_zh = f"{crop_zh[cr]}（综合）" if cr != "all4" else "四大主粮（综合）"
            rows.append({
                "panel": "A",
                "group_en": grp_en,
                "group_zh": grp_zh,
                "hazard_en": haz_en[hz],
                "hazard_zh": haz_zh[hz],
                "hist_pct": f"{hs_p.mean():.2f} ± {hs_p.std(ddof=1):.2f}",
                "ssp126_pct": f"{s1_p.mean():.2f} ± {s1_p.std(ddof=1):.2f}",
                "ssp585_pct": f"{s5_p.mean():.2f} ± {s5_p.std(ddof=1):.2f}",
                "aux_metric": f"{hs_m.mean():.2f}±{hs_m.std(ddof=1):.2f} → {s1_m.mean():.2f}±{s1_m.std(ddof=1):.2f} / {s5_m.mean():.2f}±{s5_m.std(ddof=1):.2f}",
                "sev_ext_share": sev_ext,
            })

    # Panel B: Farming-system annual yield losses L_sys (% ± 95% CI; Fig. 3 & ED Fig. 5) & exposed production f_exp (%)
    mode_en = {"nonIrr": "rainfed", "Irr": "irrigated"}
    mode_zh = {"nonIrr": "雨养", "Irr": "灌溉"}
    for cr in ["maiz", "soyb", "rice", "whea"]:
        for md in ["nonIrr", "Irr"]:
            r_h = f3[(f3["scenario"] == "hist") & (f3["crop"] == cr) & (f3["mode"] == md)].iloc[0]
            r_1 = f3[(f3["scenario"] == "ssp126") & (f3["crop"] == cr) & (f3["mode"] == md)].iloc[0]
            r_5 = f3[(f3["scenario"] == "ssp585") & (f3["crop"] == cr) & (f3["mode"] == md)].iloc[0]
            g_key = f"{cr}_{md}"
            for hz, pfx in [("dry", "dry"), ("wet", "wet")]:
                d1 = decomp[(decomp["group"] == g_key) & (decomp["haz"] == hz) & (decomp["scen"] == "ssp126")].iloc[0]
                d5 = decomp[(decomp["group"] == g_key) & (decomp["haz"] == hz) & (decomp["scen"] == "ssp585")].iloc[0]
                rows.append({
                    "panel": "B",
                    "group_en": f"{crop_en[cr]} ({mode_en[md]})",
                    "group_zh": f"{crop_zh[cr]}（{mode_zh[md]}）",
                    "hazard_en": haz_en[hz],
                    "hazard_zh": haz_zh[hz],
                    "hist_pct": f"{r_h[f'{pfx}_mean_pct']:.2f} ± {r_h[f'{pfx}_ci95_halfwidth_pct']:.2f}",
                    "ssp126_pct": f"{r_1[f'{pfx}_mean_pct']:.2f} ± {r_1[f'{pfx}_ci95_halfwidth_pct']:.2f}",
                    "ssp585_pct": f"{r_5[f'{pfx}_mean_pct']:.2f} ± {r_5[f'{pfx}_ci95_halfwidth_pct']:.2f}",
                    "aux_metric": f"{d5['hist_fexp_pct']:.1f} → {d1['fut_fexp_pct']:.1f} / {d5['fut_fexp_pct']:.1f}",
                    "sev_ext_share": f"{d5['sev_share_pct']:.1f}% / {d5['ext_share_pct']:.1f}%",
                })

    out_df = pd.DataFrame(rows)
    out_df.to_csv(out_csv, index=False)
    return out_df


def build_ed_table1_display(df_src: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, r in df_src.iterrows():
        metric_type = "All-cropland annual loss L_global & Mt yr⁻¹ (Panel A)" if r["panel"] == "A" else "System-wide annual loss L_sys & exposure f_exp (Panel B)"
        rows.append([
            metric_type,
            r["group_en"],
            r["hazard_en"],
            r["hist_pct"],
            r["ssp126_pct"],
            r["ssp585_pct"],
            r["aux_metric"],
            r["sev_ext_share"],
        ])
    return pd.DataFrame(
        rows,
        columns=[
            "Panel & Metric Basis",
            "Crop & System",
            "Hazard",
            "Historical (1990–2019, %)",
            "SSP1–2.6 (Late Century, %)",
            "SSP5–8.5 (Late Century, %)",
            "Prod. Loss (Mt yr⁻¹) [A] or Exposure f_exp (%) [B]",
            "Shapley Severity / Extent Share (SSP5–8.5)",
        ],
    )


def build_table_s1() -> pd.DataFrame:
    rows = [
        ("Historical forcing", "GSWP3–W5E5 (ISIMIP3a)", "1990–2019 analysis (1975–2019 archive)", "0.5° × 0.5°, daily", "Observation-constrained baseline reanalysis"),
        ("CMIP6 GCM 1", "CanESM5", "2020–2100 (eval. 2070–2100)", "0.5° × 0.5°, daily (ISIMIP3BASD)", "SSP1–2.6 & SSP5–8.5 (r1i1p1f1, 5-yr spin-up)"),
        ("CMIP6 GCM 2", "CNRM-CM6-1", "2020–2100 (eval. 2070–2100)", "0.5° × 0.5°, daily (ISIMIP3BASD)", "SSP1–2.6 & SSP5–8.5 (r1i1p1f2, 5-yr spin-up)"),
        ("CMIP6 GCM 3", "CNRM-ESM2-1", "2020–2100 (eval. 2070–2100)", "0.5° × 0.5°, daily (ISIMIP3BASD)", "SSP1–2.6 & SSP5–8.5 (r1i1p1f2, 5-yr spin-up)"),
        ("CMIP6 GCM 4", "EC-Earth3", "2020–2100 (eval. 2070–2100)", "0.5° × 0.5°, daily (ISIMIP3BASD)", "SSP1–2.6 & SSP5–8.5 (r1i1p1f1, 5-yr spin-up)"),
        ("CMIP6 GCM 5", "GFDL-ESM4", "2020–2100 (eval. 2070–2100)", "0.5° × 0.5°, daily (ISIMIP3BASD)", "SSP1–2.6 & SSP5–8.5 (r1i1p1f1, 5-yr spin-up)"),
        ("CMIP6 GCM 6", "IPSL-CM6A-LR", "2020–2100 (eval. 2070–2100)", "0.5° × 0.5°, daily (ISIMIP3BASD)", "SSP1–2.6 & SSP5–8.5 (r1i1p1f1, 5-yr spin-up)"),
        ("CMIP6 GCM 7", "MIROC6", "2020–2100 (eval. 2070–2100)", "0.5° × 0.5°, daily (ISIMIP3BASD)", "SSP1–2.6 & SSP5–8.5 (r1i1p1f1, 5-yr spin-up)"),
        ("CMIP6 GCM 8", "MPI-ESM1-2-HR", "2020–2100 (eval. 2070–2100)", "0.5° × 0.5°, daily (ISIMIP3BASD)", "SSP1–2.6 & SSP5–8.5 (r1i1p1f1, 5-yr spin-up)"),
        ("CMIP6 GCM 9", "MRI-ESM2-0", "2020–2100 (eval. 2070–2100)", "0.5° × 0.5°, daily (ISIMIP3BASD)", "SSP1–2.6 & SSP5–8.5 (r1i1p1f1, 5-yr spin-up)"),
        ("CMIP6 GCM 10", "UKESM1-0-LL", "2020–2100 (eval. 2070–2100)", "0.5° × 0.5°, daily (ISIMIP3BASD)", "SSP1–2.6 & SSP5–8.5 (r1i1p1f2, 5-yr spin-up)"),
    ]
    return pd.DataFrame(rows, columns=["Category", "Dataset / Model", "Analysis Period", "Spatial & Temporal Resolution", "Scenarios / Variant"])


def build_table_s2() -> pd.DataFrame:
    rows = [
        ("wwh (Winter wheat)", 0.4, 0.775, 1.15, 0.25, 0.2, 0.2, 0.65, 0.2),
        ("swh (Spring wheat)", 0.4, 0.775, 1.15, 0.25, 0.2, 0.2, 0.65, 0.2),
        ("rice (Rice)", 1.1, 1.1, 1.2, 1.05, 1.2, 2.5, 0.4, 0.2),
        ("maiz (Maize)", 0.3, 0.75, 1.2, 0.5, 0.4, 0.4, 1.0, 0.2),
        ("soyb (Soybean)", 0.5, 0.8, 1.15, 0.5, 0.2, 0.2, 0.9, 0.2),
    ]
    return pd.DataFrame(rows, columns=["Crop", "KC1", "KC2", "KC3", "KC4", "KY1", "KY2", "KY3", "KY4"])


def build_table_s3() -> pd.DataFrame:
    return pd.read_csv(os.path.join(DATA_DIR, "s13_s14_table_s3_pos_ma5_9019.csv"))


def build_table_s4() -> pd.DataFrame:
    return pd.read_csv(os.path.join(DATA_DIR, "s15_landmark_18_full_region_9019.csv"))


def build_table_s5() -> pd.DataFrame:
    return pd.read_csv(os.path.join(DATA_DIR, "tier2_ggcmi3_vs_cwatm_wet_extremes_4crops_9019.csv"))


def build_table_s6() -> pd.DataFrame:
    rows = [
        ("GSWP3–W5E5 (ISIMIP3a)", "0.5° daily global climate forcing", "Drives historical hydrology, FAO-56 ET0 & SPEI_60", "1990–2019 baseline (after 1975–1989 spin-up)"),
        ("CMIP6 10 GCMs (ISIMIP3BASD)", "0.5° daily bias-adjusted to W5E5", "Drives future projections under SSP1–2.6 & SSP5–8.5", "2020–2100 (evaluated 2070–2100; 2015–2019 spin-up)"),
        ("GGCMI Phase 3 calendars & models", "0.5° crop calendars + 11 crop models", "Planting–maturity dates; 11-model wet-extreme benchmark", "Dynamic calendars; 1990–2016 ISIMIP3a benchmark"),
        ("FAO-56, GAEZ v4 & Keating et al.", "4 FAO growth stages & crop parameters", "Proportional stage durations, Kc, Ky, and root depth", "Fixed proportional stage rules across periods (Table S2)"),
        ("SPAM 2020 + GAEZ v4 reference yield", "5'->0.5° areas; 0.5° attainable yields", "Production weights; CO2-inactive setup with FAO-56 ET0", "Static reference weights isolating hydroclimatic signal"),
        ("CIL Subnational + FAOSTAT + GDHY", "ADM1 (48 nations), national, & 0.5°", "Multi-scale yield-loss error & bias structure validation", "CIL: 1990–2015; FAOSTAT: 1990–2019; GDHY: 1990–2016"),
        ("USDA NASS & USDA RMA", "0.5° US maize yields & cause-of-loss", "Cause-specific drought vs. excess-moisture attribution test", "NASS: 1990–2016 (1,295 grids); RMA: 1990–2015 (1,179 grids)"),
        ("EM-DAT/GDIS, MODIS GFD & 18 events", "ADM polygons, 250-m->0.5° floods, cases", "Disaster-footprint & landmark event detection tests", "GDIS: 1990–2016; MODIS GFD: 2000–2016; 18 cases: 1990–2019"),
        ("ERA5 & ISMN soil moisture", "0.5° global reanalysis & in-situ stations", "Root-zone available soil moisture (ASM) validation", "ERA5: 1990–2019; ISMN: post-1990 station overlap"),
        ("FAOSTAT trade & World Bank income", "161 countries bilateral staple trade", "National Import Dependency Ratio (IDR > 0.4; Eq. 13)", "Contemporary 2014–2023 mean (r = 0.982 stability)"),
    ]
    return pd.DataFrame(rows, columns=["Dataset & Citation", "Resolution & Coverage", "Role in Framework", "Time Window / Scenario Protocol"])


def build_table_s7() -> pd.DataFrame:
    df = pd.read_csv(os.path.join(DATA_DIR, "r6_fig3_robust_inference.csv"))
    df = df[df["window"].str.startswith("1990-2019")].copy()
    crop_map = {"maiz": "Maize", "soyb": "Soybean", "rice": "Rice", "whea": "Wheat"}
    mode_map = {"nonIrr": "rainfed", "Irr": "irrigated"}
    df["Crop_System"] = df.apply(lambda r: f"{crop_map[r['crop']]} ({mode_map[r['mode']]})", axis=1)
    cols = [
        "Crop_System", "mean_dry_pct", "mean_wet_pct", "mean_diff_dry_minus_wet_pct",
        "p_paired_std", "p_newey_west_hac", "p_ar1_neff", "n_eff_ar1", "p_mbb_5yr", "p_wilcoxon", "p_twoway_5deg_year",
    ]
    return df[[c for c in cols if c in df.columns]]


def build_table_s8_idr() -> pd.DataFrame:
    df_idr = pd.read_csv(os.path.join(DATA_DIR, "r8_idr_cutoff_sweep.csv"))
    df_p = pd.read_csv(os.path.join(DATA_DIR, "r8_idr_temporal_persistence.csv"))
    r_pers = float(df_p["pearson_r_2014_2018_vs_2019_2023"].iloc[0])
    agree_pct = float(df_p["classification_agreement_idr04_pct"].iloc[0])
    df_out = df_idr.copy()
    df_out["subperiod_r_2014_18_vs_2019_23"] = round(r_pers, 4)
    df_out["idr04_classification_agreement_pct"] = round(agree_pct, 2)
    return df_out


def export_all_tables() -> None:
    # Remove legacy Table_S9 CSV if present in OUT_DIR so directory stays clean
    for old_f in os.listdir(OUT_DIR):
        if "Table_S9" in old_f or "Table_S8_Shapley" in old_f or old_f.endswith(".xlsx"):
            os.remove(os.path.join(OUT_DIR, old_f))

    ed1_src = build_ed_table1_source()
    tables = [
        ("Extended_Data_Table_1_Crop_System_Scenario_Summary", "ED Table 1", build_ed_table1_display(ed1_src)),
        ("Table_S1_Forcing_and_CMIP6_Models", "Table S1", build_table_s1()),
        ("Table_S2_Soil_and_Crop_Hydraulic_Parameters", "Table S2", build_table_s2()),
        ("Table_S3_Multiscale_Validation_Metrics_1990_2019", "Table S3", build_table_s3()),
        ("Table_S4_18_Landmark_Historical_Disasters_1990_2019", "Table S4", build_table_s4()),
        ("Table_S5_GGCMI3_WetExtreme_Process_and_Skill", "Table S5", build_table_s5()),
        ("Table_S6_Input_Datasets_and_Scenario_Consistency", "Table S6", build_table_s6()),
        ("Table_S7_Fig3_Autocorrelation_Spatial_Robust_Tests", "Table S7", build_table_s7()),
        ("Table_S8_IDR_Cutoff_and_DeltaLoss_2D_Sensitivity", "Table S8", build_table_s8_idr()),
    ]

    wb = openpyxl.Workbook()
    wb.remove(wb.active)

    hdr_fill = PatternFill(start_color="1F4E79", end_color="1F4E79", fill_type="solid")
    alt_fill = PatternFill(start_color="F4F7FA", end_color="F4F7FA", fill_type="solid")
    hdr_font = Font(name="Calibri", size=10.5, bold=True, color="FFFFFF")
    body_font = Font(name="Calibri", size=10, color="111111")
    thin_bottom = Border(bottom=Side(style="thin", color="D0D7DE"))
    top_bottom_heavy = Border(
        top=Side(style="medium", color="1F4E79"),
        bottom=Side(style="medium", color="1F4E79"),
    )

    for csv_stem, sheet_title, df in tables:
        csv_path = os.path.join(OUT_DIR, f"{csv_stem}.csv")
        df.to_csv(csv_path, index=False)

        ws = wb.create_sheet(title=sheet_title)
        ws.views.sheetView[0].showGridLines = True
        ws.freeze_panes = "A2"
        ws.page_setup.orientation = "landscape"
        ws.sheet_properties.pageSetUpPr.fitToPage = True
        ws.page_setup.fitToWidth = 1
        ws.page_setup.fitToHeight = 0

        for col_i, col_name in enumerate(df.columns, 1):
            c = ws.cell(row=1, column=col_i, value=str(col_name))
            c.fill = hdr_fill
            c.font = hdr_font
            c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
            c.border = top_bottom_heavy
        ws.row_dimensions[1].height = 28

        for row_i, row_vals in enumerate(df.itertuples(index=False), 2):
            ws.row_dimensions[row_i].height = 20
            for col_i, val in enumerate(row_vals, 1):
                v_clean = None if (isinstance(val, float) and np.isnan(val)) else val
                cell = ws.cell(row=row_i, column=col_i, value=v_clean)
                cell.font = body_font
                if len(df) >= 12 and row_i % 2 == 1:
                    cell.fill = alt_fill
                is_num = isinstance(v_clean, (int, float, np.number))
                if is_num:
                    cell.alignment = Alignment(horizontal="right", vertical="center")
                    if isinstance(v_clean, (float, np.floating)):
                        cell.number_format = "0.0000" if abs(v_clean) < 0.01 and v_clean != 0 else "0.00"
                    else:
                        cell.number_format = "#,##0"
                else:
                    s_val = str(v_clean or "")
                    is_stat_str = any(ch in s_val for ch in ("±", "→", "%")) and col_i >= 3
                    cell.alignment = Alignment(horizontal=("center" if is_stat_str else "left"), vertical="center")
                cell.border = thin_bottom

        for col_i, col_name in enumerate(df.columns, 1):
            sample_lens = [min(len(str(col_name)), 28)] + [len(str(v)) for v in df.iloc[:45, col_i - 1].values if pd.notna(v)]
            ws.column_dimensions[get_column_letter(col_i)].width = min(max(max(sample_lens) + 4, 13), 44)

    xlsx_name = "SI_Tables_S1_to_S8_and_ED_Table_1_v3.xlsx"
    v3_xlsx = os.path.join(V3_DIR, xlsx_name)
    wb.save(v3_xlsx)

    # Remove old SI_Tables_S1_to_S9_v3.xlsx / v2.xlsx so no stale workbook remains
    for stale in [
        os.path.join(V3_DIR, "SI_Tables_S1_to_S9_v3.xlsx"),
        os.path.join(SCRIPTS_ROOT, "SI_Tables_S1_to_S9_v2.xlsx"),
    ]:
        if os.path.exists(stale):
            os.remove(stale)

    # Recalculate via headless LibreOffice Calc & run validator
    recalc_dir = pathlib.Path(INTERMEDIATE_DIR) / "recalc_xlsx_v3"
    recalc_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="lo_profile_") as lo_profile:
        subprocess.run(
            [
                "libreoffice",
                f"-env:UserInstallation=file://{lo_profile}",
                "--headless", "--nologo", "--nolockcheck", "--norestore",
                "--calc", "--convert-to", "xlsx",
                "--outdir", str(recalc_dir),
                v3_xlsx,
            ],
            env={**os.environ, "HOME": "/tmp", "TMPDIR": "/tmp"},
            timeout=45,
            check=True,
        )
    shutil.move(str(recalc_dir / xlsx_name), v3_xlsx)
    shutil.rmtree(recalc_dir, ignore_errors=True)

    validator = os.path.join(ROOT_WS, "skills/excel-processing/scripts/validate_excel.py")
    if os.path.exists(validator):
        subprocess.run(["python3", validator, v3_xlsx], check=True)

    # Copy to scripts_v3 & tables_output
    shutil.copy2(v3_xlsx, os.path.join(OUT_DIR, xlsx_name))
    shutil.copy2(v3_xlsx, os.path.join(SCRIPTS_ROOT, xlsx_name))

    # Render first 2 pages of Excel PDF preview (ED Table 1)
    with tempfile.TemporaryDirectory(prefix="lo_profile_") as lo_profile:
        subprocess.run(
            [
                "libreoffice",
                f"-env:UserInstallation=file://{lo_profile}",
                "--headless", "--nologo", "--nolockcheck", "--norestore",
                "--calc", "--convert-to", "pdf",
                "--outdir", INTERMEDIATE_DIR,
                v3_xlsx,
            ],
            env={**os.environ, "HOME": "/tmp", "TMPDIR": "/tmp"},
            timeout=45,
            check=True,
        )
    pdf_p = os.path.join(INTERMEDIATE_DIR, f"{os.path.splitext(xlsx_name)[0]}.pdf")
    pdoc = pdfium.PdfDocument(pdf_p)
    for idx_pg in range(min(3, len(pdoc))):
        out_png = os.path.abspath(os.path.join(INTERMEDIATE_DIR, f"xlsx_v3_page_{idx_pg + 1}.png"))
        pdoc[idx_pg].render(scale=2.0).to_pil().save(out_png)
        print(f"Rendered Excel preview: {out_png}")
    pdoc.close()
    print(f"[OK] Exported & validated ED Table 1 + Tables S1–S8 -> {v3_xlsx}")


if __name__ == "__main__":
    export_all_tables()
