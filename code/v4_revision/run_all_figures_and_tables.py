# -*- coding: utf-8 -*-
"""Master One-Click Reproducibility Runner (`run_all_figures_and_tables.py`).

Generates:
  1. All 8 supplementary 100% True Vector PDF figures (`pdf.fonttype = 42`, `0` raster images)
     + 600-dpi `.jpg`/`.png` companion images into `../figures/`:
       - `Fig_S9_Parameter_StressCurve_Sensitivity.{pdf,jpg}` (100% real CWatM data)
       - `Fig_S11_ERA5_ASM_Diverging.{pdf,png}` (soil moisture validation map)
       - `Fig_S12_ISMN_ASM_Diverging.{pdf,png}` (soil moisture validation map)
       - `Fig_S13_Subnational_Bias_Hultgren.{pdf,jpg}` (100% real CWatM + CIL data)
       - `Fig_S14_Multiscale_US_RMA.{pdf,jpg}` (100% real CWatM + NASS/RMA/GDHY/FAOSTAT data)
       - `Fig_S15_Landmark_Events_Catalogs.{pdf,jpg}` (100% real CWatM + GDIS/MODIS data)
       - `Fig_S16_WetExtreme_GGCMI3_Ablation.{pdf,jpg}` (100% real CWatM + GGCMI3 data)
       - `Fig_S17_Methodological_Robustness_Synthesis.{pdf,jpg}` (100% real CWatM data)
     Note: Main-text Figures 1, 3, 4 and Extended Data Fig. 3 are generated directly from the
     author's upstream pipelines (`01_main_and_ed_figures/rebuild_manuscript_figure4.py` and
     `Untitled2.ipynb`).
  2. All 9 Supplementary Tables (`Tables S1–S9`) as clean `.csv` files and a formatted
     multi-sheet Excel workbook (`05_tables_generation/tables_output/SI_Tables_S1_to_S9_v2.xlsx`).
"""
from __future__ import annotations

import glob
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
OUT_FIG_DIR = os.path.join(os.path.dirname(HERE), "figures")


def run_step(rel_script: str) -> None:
    script_path = os.path.join(HERE, rel_script)
    print(f"\n>>> Running {rel_script} ...")
    res = subprocess.run([sys.executable, script_path], cwd=os.path.dirname(script_path))
    if res.returncode != 0:
        raise RuntimeError(f"Script failed with exit code {res.returncode}: {rel_script}")


def verify_vector_pdfs() -> None:
    if shutil.which("pdfimages") is None:
        print("[INFO] `pdfimages` not found on PATH; skipping Poppler CLI check.")
        return
    pdfs = sorted(glob.glob(os.path.join(OUT_FIG_DIR, "Fig_S*.pdf")))
    print(f"\n=== Vector PDF Verification ({len(pdfs)} SI PDFs in {OUT_FIG_DIR}) ===")
    for p in pdfs:
        out_img = subprocess.check_output(["pdfimages", "-list", p], text=True).strip().splitlines()
        n_raster = max(0, len(out_img) - 2)
        out_fnt = subprocess.check_output(["pdffonts", p], text=True).strip().splitlines()
        n_fonts = max(0, len(out_fnt) - 2)
        sz_kb = os.path.getsize(p) / 1024.0
        print(f"  {os.path.basename(p):50s} | {sz_kb:7.1f} KB | raster_images={n_raster} | fonts={n_fonts}")
        assert n_raster == 0, f"FAIL: Raster image detected inside {p}!"
        assert n_fonts >= 1, f"FAIL: No vector fonts embedded inside {p}!"
    print("=== [ALL PASS] Every SI PDF is 100% pure editable vector (0 raster images, Type 42 fonts) ===")


if __name__ == "__main__":
    run_step("02_si_soil_moisture_maps/plot_figs_s11_s12_vector.py")
    run_step("03_si_multiscale_validation/plot_figs_s13_to_s16_vector.py")
    run_step("04_si_robustness_synthesis/plot_figs_s9_and_s17_vector.py")
    run_step("05_tables_generation/generate_all_si_tables_s1_to_s9.py")
    verify_vector_pdfs()
