from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm


# ========= User Config =========
INPUT_CSV = Path(r"D:\AAUDE\paper_v2\paper4\validation_fao_yield\model_vs_fao_yield_country_year.csv")
OUTPUT_DIR = Path(r"D:\AAUDE\paper\paper4\revision_v2\fig\manu_fig\validation_hist_4crops")

# Keep all four crops in one common historical period.
# Change to None to use full common period automatically.
TARGET_YEAR_MIN = 1980
TARGET_YEAR_MAX = 2019

MIN_COUNTRY_YEARS = 10
TOP_LABEL_N = 8

CROPS = ["maiz", "soyb", "rice", "whea"]
CROP_TITLES = {
    "maiz": "Maize",
    "soyb": "Soybean",
    "rice": "Rice",
    "whea": "Wheat",
}


def _rmse(x: pd.Series, y: pd.Series) -> float:
    return float(np.sqrt(np.mean((x - y) ** 2)))


def build_common_period(df: pd.DataFrame, crops: list[str]) -> list[int]:
    year_sets = []
    for crop in crops:
        years = set(df.loc[df["crop"] == crop, "Year"].dropna().astype(int).tolist())
        year_sets.append(years)
    common = sorted(set.intersection(*year_sets))
    if TARGET_YEAR_MIN is not None:
        common = [y for y in common if y >= TARGET_YEAR_MIN]
    if TARGET_YEAR_MAX is not None:
        common = [y for y in common if y <= TARGET_YEAR_MAX]
    return common


def country_summary(df_period: pd.DataFrame) -> pd.DataFrame:
    out = []
    for (crop, iso), g in df_period.groupby(["crop", "ISO_A3"], dropna=True):
        g = g.dropna(subset=["model_yield_kg_ha", "fao_yield_kg_ha"]).copy()
        g = g[(g["model_yield_kg_ha"] > 0) & (g["fao_yield_kg_ha"] > 0)]
        n = len(g)
        if n < MIN_COUNTRY_YEARS:
            continue

        sim = g["model_yield_kg_ha"].mean()
        obs = g["fao_yield_kg_ha"].mean()
        if obs <= 0:
            continue

        # Time-series correlation at country level.
        if g["model_yield_kg_ha"].std() == 0 or g["fao_yield_kg_ha"].std() == 0:
            r_ts = np.nan
        else:
            r_ts = float(np.corrcoef(g["fao_yield_kg_ha"], g["model_yield_kg_ha"])[0, 1])

        out.append(
            {
                "crop": crop,
                "ISO_A3": iso,
                "n_year": n,
                "obs_mean_kg_ha": float(obs),
                "sim_mean_kg_ha": float(sim),
                "bias_pct": float((sim - obs) / obs * 100.0),
                "r_ts_country": r_ts,
            }
        )
    return pd.DataFrame(out)


def crop_metrics(df_country: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for crop, g in df_country.groupby("crop"):
        obs = g["obs_mean_kg_ha"]
        sim = g["sim_mean_kg_ha"]
        r = np.nan
        if len(g) > 1 and obs.std() > 0 and sim.std() > 0:
            r = float(np.corrcoef(obs, sim)[0, 1])
        rmse = _rmse(sim, obs)
        nrmse = float(rmse / obs.mean()) if obs.mean() > 0 else np.nan
        mape = float(np.mean(np.abs((sim - obs) / obs)) * 100.0)
        rows.append(
            {
                "crop": crop,
                "n_countries": int(len(g)),
                "pearson_r_country_mean": r,
                "rmse_kg_ha": rmse,
                "nrmse": nrmse,
                "mape_pct": mape,
                "mean_bias_pct": float(g["bias_pct"].mean()),
                "median_bias_pct": float(g["bias_pct"].median()),
                "year_min": int(df_period["Year"].min()),
                "year_max": int(df_period["Year"].max()),
            }
        )
    return pd.DataFrame(rows).sort_values("crop")


def _nice_bounds(vals: pd.Series) -> tuple[float, float]:
    lo = float(vals.min())
    hi = float(vals.max())
    if hi <= lo:
        hi = lo + 1.0
    margin = 0.07 * (hi - lo)
    return lo - margin, hi + margin


def make_figure(df_country: pd.DataFrame, metrics: pd.DataFrame, out_png: Path, out_pdf: Path) -> None:
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["font.serif"] = ["Arial"]
    mpl.rcParams["axes.unicode_minus"] = False

    fig, axes = plt.subplots(2, 2, figsize=(13, 10), dpi=220)
    axes = axes.flatten()

    # Shared color mapping across all crops for comparability.
    all_bias = df_country["bias_pct"].replace([np.inf, -np.inf], np.nan).dropna()
    q = float(np.nanpercentile(np.abs(all_bias), 95)) if len(all_bias) else 30.0
    q = max(q, 5.0)
    norm = TwoSlopeNorm(vmin=-q, vcenter=0.0, vmax=q)
    cmap = plt.get_cmap("RdBu_r")

    for i, crop in enumerate(CROPS):
        ax = axes[i]
        sub = df_country[df_country["crop"] == crop].copy()
        if sub.empty:
            ax.axis("off")
            ax.set_title(f"{CROP_TITLES[crop]} (no data)")
            continue

        x = sub["obs_mean_kg_ha"]
        y = sub["sim_mean_kg_ha"]
        size = 16 + (sub["n_year"] - sub["n_year"].min()) * 2.0

        sc = ax.scatter(
            x,
            y,
            c=sub["bias_pct"],
            cmap=cmap,
            norm=norm,
            s=size,
            edgecolor="#4a4a4a",
            linewidth=0.35,
            alpha=0.88,
        )

        x0, x1 = _nice_bounds(x)
        y0, y1 = _nice_bounds(y)
        lo = min(x0, y0)
        hi = max(x1, y1)
        ax.plot([lo, hi], [lo, hi], linestyle="--", color="#222222", linewidth=1.25)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)

        # Label strongest outliers to mimic event-labeled style.
        sub["abs_bias"] = np.abs(sub["bias_pct"])
        for _, r in sub.nlargest(TOP_LABEL_N, "abs_bias").iterrows():
            ax.annotate(
                r["ISO_A3"],
                (r["obs_mean_kg_ha"], r["sim_mean_kg_ha"]),
                fontsize=8,
                xytext=(3, 3),
                textcoords="offset points",
                color="#222222",
            )

        m = metrics[metrics["crop"] == crop].iloc[0]
        text = (
            f"N={int(m['n_countries'])}\n"
            f"R={m['pearson_r_country_mean']:.2f}\n"
            f"RMSE={m['rmse_kg_ha']:.0f} kg/ha\n"
            f"MAPE={m['mape_pct']:.1f}%"
        )
        ax.text(
            0.03,
            0.97,
            text,
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.8, "pad": 4},
        )

        ax.set_title(CROP_TITLES[crop], fontsize=12)
        ax.set_xlabel("Observed yield (kg/ha)")
        ax.set_ylabel("Simulated yield (kg/ha)")
        ax.grid(True, linestyle="-", alpha=0.22)

    cbar = fig.colorbar(
        mpl.cm.ScalarMappable(norm=norm, cmap=cmap),
        ax=axes.tolist(),
        fraction=0.023,
        pad=0.02,
    )
    cbar.set_label("Bias = (Simulated - Observed) / Observed (%)", fontsize=11)

    yr_min = int(df_period["Year"].min())
    yr_max = int(df_period["Year"].max())
    fig.suptitle(
        f"Historical Yield Validation Across Four Crops ({yr_min}-{yr_max}, common period)",
        fontsize=14,
        y=0.98,
    )
    fig.tight_layout(rect=[0, 0, 0.97, 0.96])
    fig.savefig(out_png, dpi=320)
    fig.savefig(out_pdf)
    plt.close(fig)


if __name__ == "__main__":
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    if not INPUT_CSV.exists():
        raise FileNotFoundError(f"Missing input file: {INPUT_CSV}")

    df = pd.read_csv(INPUT_CSV)
    required = {"ISO_A3", "Year", "model_yield_kg_ha", "fao_yield_kg_ha", "crop"}
    miss = required - set(df.columns)
    if miss:
        raise ValueError(f"Missing required columns: {sorted(miss)}")

    df = df[df["crop"].isin(CROPS)].copy()
    df["Year"] = pd.to_numeric(df["Year"], errors="coerce")
    df["model_yield_kg_ha"] = pd.to_numeric(df["model_yield_kg_ha"], errors="coerce")
    df["fao_yield_kg_ha"] = pd.to_numeric(df["fao_yield_kg_ha"], errors="coerce")
    df = df.dropna(subset=["ISO_A3", "Year", "crop", "model_yield_kg_ha", "fao_yield_kg_ha"])
    df = df[(df["model_yield_kg_ha"] > 0) & (df["fao_yield_kg_ha"] > 0)].copy()
    df["Year"] = df["Year"].astype(int)

    common_years = build_common_period(df, CROPS)
    if not common_years:
        raise RuntimeError("No common years found across all four crops under current year constraints.")

    df_period = df[df["Year"].isin(common_years)].copy()
    df_country = country_summary(df_period)
    metrics = crop_metrics(df_country)

    out_country_csv = OUTPUT_DIR / "validation_country_summary_common_period.csv"
    out_metrics_csv = OUTPUT_DIR / "validation_crop_metrics_common_period.csv"
    out_years_txt = OUTPUT_DIR / "validation_common_years.txt"
    out_png = OUTPUT_DIR / "validation_scatter_4crops_common_period.png"
    out_pdf = OUTPUT_DIR / "validation_scatter_4crops_common_period.pdf"

    df_country.to_csv(out_country_csv, index=False, encoding="utf-8-sig")
    metrics.to_csv(out_metrics_csv, index=False, encoding="utf-8-sig")
    out_years_txt.write_text(
        f"Common years ({len(common_years)}): {min(common_years)}-{max(common_years)}\n",
        encoding="utf-8",
    )

    make_figure(df_country, metrics, out_png, out_pdf)

    print("Done.")
    print("Input:", INPUT_CSV)
    print("Output dir:", OUTPUT_DIR)
    print("Common years:", min(common_years), "-", max(common_years), "n=", len(common_years))
    print("Country summary:", out_country_csv)
    print("Crop metrics:", out_metrics_csv)
    print("Figure PNG:", out_png)
    print("Figure PDF:", out_pdf)
