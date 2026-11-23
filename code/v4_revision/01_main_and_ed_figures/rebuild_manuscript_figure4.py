from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


CROPS = ["maiz", "soyb", "rice", "whea"]
CROP_LABELS = {
    "maiz": "Maize",
    "soyb": "Soybean",
    "rice": "Rice",
    "whea": "Wheat",
}
SCENARIOS = ["ssp126", "ssp585"]
SCENARIO_LABELS = {
    "historical": "Historical",
    "ssp126": "SSP1-2.6",
    "ssp585": "SSP5-8.5",
}
HAZARDS = ["dry", "wet"]
HAZARD_LABELS = {
    "dry": "Drought",
    "wet": "Waterlogging",
}

COLORS = {
    "historical": "#5F6467",
    "ssp126": "#287D8E",
    "ssp585": "#C76F3A",
}
FILL_ALPHA = 0.16
BIAS_REF_HIST = (2015, 2019)
BIAS_REF_FUT = (2020, 2024)


def set_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "font.size": 7.2,
            "axes.titlesize": 8.8,
            "axes.labelsize": 8.2,
            "xtick.labelsize": 7.3,
            "ytick.labelsize": 7.3,
            "legend.fontsize": 7.3,
            "axes.linewidth": 0.65,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 3,
            "ytick.major.size": 3,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "figure.dpi": 160,
        }
    )


def centered_rolling(series: pd.Series, window: int) -> pd.Series:
    if window <= 1:
        return series
    return series.rolling(window=window, center=True, min_periods=1).mean()


def load_and_filter(input_csv: Path, metric: str, mode: str) -> pd.DataFrame:
    df = pd.read_csv(input_csv)
    required = {
        "scenario",
        "gcm",
        "year",
        "crop",
        "mode",
        "hazard",
        metric,
    }
    missing = sorted(required.difference(df.columns))
    if missing:
        raise ValueError(f"Missing required columns in {input_csv}: {missing}")

    out = df.loc[
        (df["mode"] == mode)
        & (df["crop"].isin(CROPS))
        & (df["hazard"].isin(HAZARDS))
        & (df["scenario"].isin(["historical", *SCENARIOS])),
        ["scenario", "gcm", "year", "crop", "hazard", metric],
    ].copy()
    out = out.rename(columns={metric: "value"})
    out["year"] = out["year"].astype(int)
    out["value"] = pd.to_numeric(out["value"], errors="coerce")
    out = out.dropna(subset=["value"])
    return out


def bias_correct_future_to_recent_history(
    df: pd.DataFrame,
    hist_ref: tuple[int, int] = BIAS_REF_HIST,
    fut_ref: tuple[int, int] = BIAS_REF_FUT,
    clip_zero: bool = True,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Anchor each future GCM to the recent historical mean for the same crop and hazard.

    This removes a visual/structural jump at the 2019/2020 boundary when future GCM
    sequences do not include a historical overlap in the already summarized table.
    """
    out = df.copy()
    hist_anchor = (
        out.loc[
            (out["scenario"].eq("historical"))
            & (out["year"].between(hist_ref[0], hist_ref[1]))
        ]
        .groupby(["crop", "hazard"], as_index=False)["value"]
        .mean()
        .rename(columns={"value": "hist_anchor"})
    )
    fut_anchor = (
        out.loc[
            (out["scenario"].isin(SCENARIOS))
            & (out["year"].between(fut_ref[0], fut_ref[1]))
        ]
        .groupby(["scenario", "gcm", "crop", "hazard"], as_index=False)["value"]
        .mean()
        .rename(columns={"value": "fut_anchor"})
    )
    offsets = fut_anchor.merge(hist_anchor, on=["crop", "hazard"], how="left")
    offsets["bias_offset"] = offsets["hist_anchor"] - offsets["fut_anchor"]

    out = out.merge(
        offsets[["scenario", "gcm", "crop", "hazard", "bias_offset"]],
        on=["scenario", "gcm", "crop", "hazard"],
        how="left",
    )
    fut_mask = out["scenario"].isin(SCENARIOS)
    out["bias_offset"] = out["bias_offset"].fillna(0.0)
    out.loc[fut_mask, "value_raw"] = out.loc[fut_mask, "value"]
    out.loc[fut_mask, "value"] = out.loc[fut_mask, "value"] + out.loc[fut_mask, "bias_offset"]
    if clip_zero:
        out.loc[fut_mask, "value"] = out.loc[fut_mask, "value"].clip(lower=0)
    out.loc[~fut_mask, "value_raw"] = out.loc[~fut_mask, "value"]
    return out.drop(columns=["bias_offset"]), offsets


def build_timeseries(df: pd.DataFrame, start_year: int, rolling_window: int) -> pd.DataFrame:
    hist = (
        df.loc[(df["scenario"] == "historical") & (df["year"].between(1980, 2019))]
        .groupby(["crop", "hazard", "year"], as_index=False)["value"]
        .mean()
        .assign(q25=lambda x: x["value"], q75=lambda x: x["value"], scenario="historical")
    )

    fut_raw = df.loc[
        (df["scenario"].isin(SCENARIOS)) & (df["year"].between(2020, 2100))
    ].copy()
    fut = (
        fut_raw.groupby(["scenario", "crop", "hazard", "year"])["value"]
        .agg(value="mean", q25=lambda x: x.quantile(0.25), q75=lambda x: x.quantile(0.75), n_gcm="count")
        .reset_index()
    )
    hist["n_gcm"] = 1

    ts = pd.concat([hist, fut], ignore_index=True, sort=False)
    ts = ts.sort_values(["scenario", "crop", "hazard", "year"])
    for col in ["value", "q25", "q75"]:
        ts[f"{col}_smooth"] = (
            ts.groupby(["scenario", "crop", "hazard"], sort=False)[col]
            .transform(lambda s: centered_rolling(s, rolling_window))
        )
    # Smooth the full historical series before trimming the displayed years so
    # the 1990 MA5 value retains its original 1988-1992 window.
    ts = ts.loc[(ts["scenario"] != "historical") | (ts["year"] >= start_year)].copy()
    return ts


def build_box_data(df: pd.DataFrame, period_start: int, period_end: int) -> pd.DataFrame:
    box = (
        df.loc[
            (df["scenario"].isin(SCENARIOS))
            & (df["year"].between(period_start, period_end))
        ]
        .groupby(["scenario", "crop", "hazard", "gcm"], as_index=False)
        .agg(value=("value", "mean"), n_years=("year", "nunique"))
    )
    complete_years = period_end - period_start + 1
    box["period"] = f"{period_start}-{period_end}"
    box["complete_period"] = box["n_years"].eq(complete_years)
    return box


def nice_ylim(values: pd.Series, pad: float = 0.12) -> tuple[float, float]:
    upper = float(np.nanmax(values))
    if not np.isfinite(upper) or upper <= 0:
        return 0, 1
    upper = upper * (1 + pad)
    if upper <= 2:
        step = 0.25
    elif upper <= 5:
        step = 0.5
    elif upper <= 10:
        step = 1
    else:
        step = 2
    return 0, np.ceil(upper / step) * step


def draw_boxplot(ax: plt.Axes, values_by_key: dict[tuple[str, str], np.ndarray], ylim: tuple[float, float]) -> None:
    positions = [0.78, 1.22, 1.88, 2.32]
    keys = [("ssp126", "dry"), ("ssp126", "wet"), ("ssp585", "dry"), ("ssp585", "wet")]
    data = [values_by_key.get(key, np.array([])) for key in keys]

    bp = ax.boxplot(
        data,
        positions=positions,
        widths=0.28,
        patch_artist=True,
        showfliers=False,
        medianprops={"color": "#252525", "linewidth": 0.9},
        whiskerprops={"color": "#565656", "linewidth": 0.65},
        capprops={"color": "#565656", "linewidth": 0.65},
        boxprops={"linewidth": 0.75},
    )
    for patch, (scenario, hazard) in zip(bp["boxes"], keys):
        patch.set_facecolor(COLORS[scenario])
        patch.set_alpha(0.52 if hazard == "dry" else 0.24)
        patch.set_edgecolor(COLORS[scenario])
        if hazard == "wet":
            patch.set_hatch("///")

    rng = np.random.default_rng(20260512)
    for pos, key in zip(positions, keys):
        vals = values_by_key.get(key, np.array([]))
        if vals.size == 0:
            continue
        jitter = rng.uniform(-0.10, 0.10, vals.size)
        scenario, hazard = key
        ax.scatter(
            np.full(vals.size, pos) + jitter,
            vals,
            s=8,
            color=COLORS[scenario],
            alpha=0.65,
            linewidths=0,
            zorder=3,
        )
        mean_val = float(np.nanmean(vals))
        ax.plot(
            [pos - 0.13, pos + 0.13],
            [mean_val, mean_val],
            color="#111111",
            linewidth=0.95,
            zorder=4,
        )
        ax.text(
            pos + 0.16,
            mean_val,
            f"{mean_val:.1f}",
            ha="left",
            va="center",
            fontsize=5.6,
            color="#222222",
            zorder=5,
        )

    ax.set_xlim(0.50, 2.78)
    ax.set_ylim(*ylim)
    ax.set_xticks([])
    ax.tick_params(axis="y", labelleft=False, left=False)
    ax.tick_params(axis="x", pad=1, length=0)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.set_title("2070-2100", pad=2, fontsize=6.6)


def plot_figure(ts: pd.DataFrame, box: pd.DataFrame, out_prefix: Path, metric_label: str, start_year: int = 1990) -> None:
    set_style()
    fig = plt.figure(figsize=(7.35, 5.45))
    outer = GridSpec(
        2,
        2,
        figure=fig,
        wspace=0.30,
        hspace=0.45,
    )
    crop_positions = {
        "maiz": (0, 0),
        "soyb": (0, 1),
        "rice": (1, 0),
        "whea": (1, 1),
    }
    panel_letters = {"maiz": "a", "soyb": "b", "rice": "c", "whea": "d"}

    line_style = {"dry": "-", "wet": (0, (3.2, 2.0))}
    all_axes = []
    for crop in CROPS:
        row, col = crop_positions[crop]
        inner = outer[row, col].subgridspec(1, 2, width_ratios=[5.0, 0.64], wspace=0.025)
        ax = fig.add_subplot(inner[0, 0])
        bax = fig.add_subplot(inner[0, 1], sharey=ax)
        all_axes.append(ax)

        crop_ts = ts.loc[ts["crop"].eq(crop)]
        crop_box = box.loc[box["crop"].eq(crop)]
        ylim_values = pd.concat(
            [
                crop_ts["q75_smooth"],
                crop_ts["value_smooth"],
                crop_box["value"],
            ],
            ignore_index=True,
        )
        ylim = nice_ylim(ylim_values)

        for hazard in HAZARDS:
            sub = crop_ts.loc[
                (crop_ts["scenario"].eq("historical")) & (crop_ts["hazard"].eq(hazard))
            ].sort_values("year")
            ax.plot(
                sub["year"],
                sub["value_smooth"],
                color=COLORS["historical"],
                linewidth=1.15,
                linestyle=line_style[hazard],
                alpha=0.95,
                solid_capstyle="round",
                label=f"Historical - {HAZARD_LABELS[hazard]}",
            )

        for scenario in SCENARIOS:
            for hazard in HAZARDS:
                sub = crop_ts.loc[
                    (crop_ts["scenario"].eq(scenario)) & (crop_ts["hazard"].eq(hazard))
                ].sort_values("year")
                ax.fill_between(
                    sub["year"].to_numpy(),
                    sub["q25_smooth"].to_numpy(),
                    sub["q75_smooth"].to_numpy(),
                    color=COLORS[scenario],
                    alpha=FILL_ALPHA if hazard == "dry" else FILL_ALPHA * 0.58,
                    linewidth=0,
                    rasterized=True,
                )
                ax.plot(
                    sub["year"],
                    sub["value_smooth"],
                    color=COLORS[scenario],
                    linewidth=1.35,
                    linestyle=line_style[hazard],
                    alpha=0.98,
                    solid_capstyle="round",
                    label=f"{SCENARIO_LABELS[scenario]} - {HAZARD_LABELS[hazard]}",
                )

        ax.axvline(2020, color="#B5B5B5", linewidth=0.55, linestyle=(0, (1.2, 1.8)), zorder=0)
        ax.axvspan(2070, 2100, color="#D9D9D9", alpha=0.24, linewidth=0, zorder=-2)
        ax.set_xlim(start_year - 2, 2102)
        ax.set_ylim(*ylim)
        ax.set_title(CROP_LABELS[crop], loc="left", pad=3, fontweight="bold")
        ax.text(
            -0.095,
            1.02,
            panel_letters[crop],
            transform=ax.transAxes,
            fontweight="bold",
            fontsize=9.2,
            ha="left",
            va="bottom",
        )
        ax.set_ylabel(metric_label)
        _ticks = [start_year]
        _first_decade = ((start_year // 10) + 1) * 10
        _ticks.extend(range(_first_decade, 2101, 20))
        if 2100 not in _ticks:
            _ticks.append(2100)
        ax.set_xticks(sorted(set(_ticks)))
        ax.tick_params(axis="both", colors="#303030")
        if start_year % 20 != 0:
            for _label in ax.get_xticklabels():
                if _label.get_text() == str(start_year):
                    _label.set_horizontalalignment("left")
                    _label.set_x(_label.get_position()[0] - 1.6)

        values_by_key = {
            (scenario, hazard): crop_box.loc[
                (crop_box["scenario"].eq(scenario)) & (crop_box["hazard"].eq(hazard)),
                "value",
            ].to_numpy()
            for scenario in SCENARIOS
            for hazard in HAZARDS
        }
        draw_boxplot(bax, values_by_key, ylim)

    legend_handles = [
        Line2D([0], [0], color=COLORS["historical"], lw=1.4, linestyle="-", label="Historical"),
        Line2D([0], [0], color=COLORS["ssp126"], lw=1.6, linestyle="-", label="SSP1-2.6"),
        Line2D([0], [0], color=COLORS["ssp585"], lw=1.6, linestyle="-", label="SSP5-8.5"),
        Line2D([0], [0], color="#222222", lw=1.4, linestyle="-", label="Drought"),
        Line2D([0], [0], color="#222222", lw=1.4, linestyle=(0, (3.2, 2.0)), label="Waterlogging"),
        Patch(facecolor="#8F8F8F", alpha=FILL_ALPHA, linewidth=0, label="GCM IQR"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="lower center",
        ncol=6,
        bbox_to_anchor=(0.52, 0.015),
        columnspacing=1.35,
        handlelength=2.2,
    )
    fig.supxlabel("Year", y=0.073, fontsize=8.3)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.965, bottom=0.14)

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=450)
    fig.savefig(out_prefix.with_suffix(".pdf"))
    fig.savefig(out_prefix.with_suffix(".svg"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot Fig. 4 water-stress yield loss time series with late-century GCM boxplots."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("source_index/water_stress_loss_summary_detailed/water_stress_loss_annual_detailed.csv"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("source_index/fig4_rebuild"),
    )
    parser.add_argument("--metric", default="loss_pct", choices=["loss_pct", "loss_million_tonnes"])
    parser.add_argument("--mode", default="total", choices=["total", "Irr", "nonIrr"])
    parser.add_argument("--start-year", type=int, default=1990)
    parser.add_argument("--rolling-window", type=int, default=5)
    parser.add_argument("--box-start", type=int, default=2070)
    parser.add_argument("--box-end", type=int, default=2100)
    parser.add_argument(
        "--bias-correction",
        action="store_true",
        help="Anchor future GCM series to the 2015-2019 historical mean using the 2020-2024 future mean.",
    )
    parser.add_argument(
        "--no-bias-correction",
        action="store_true",
        help="Deprecated compatibility flag. The default is already no bias correction.",
    )
    args = parser.parse_args()

    metric_label = "Yield loss (%)" if args.metric == "loss_pct" else "Yield loss (Mt)"
    df = load_and_filter(args.input, args.metric, args.mode)
    if args.bias_correction:
        df, offsets = bias_correct_future_to_recent_history(df)
        bc_label = "bias_corrected"
    else:
        offsets = pd.DataFrame()
        bc_label = "raw"
    ts = build_timeseries(df, args.start_year, args.rolling_window)
    box = build_box_data(df, args.box_start, args.box_end)

    suffix = f"{args.metric}_{args.mode}_{bc_label}_{args.box_start}_{args.box_end}"
    out_prefix = args.output_dir / f"Fig4_water_stress_loss_timeseries_boxplots_{suffix}"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    ts.to_csv(args.output_dir / f"Fig4_timeseries_source_{suffix}.csv", index=False, encoding="utf-8-sig")
    box.to_csv(args.output_dir / f"Fig4_gcm_boxplot_source_{suffix}.csv", index=False, encoding="utf-8-sig")
    if not offsets.empty:
        offsets.to_csv(args.output_dir / f"Fig4_bias_correction_offsets_{suffix}.csv", index=False, encoding="utf-8-sig")
    plot_figure(ts, box, out_prefix, metric_label, start_year=args.start_year)

    counts = (
        box.groupby(["scenario", "crop", "hazard"], as_index=False)
        .agg(n_gcm=("gcm", "nunique"), min_years=("n_years", "min"), max_years=("n_years", "max"))
        .sort_values(["scenario", "crop", "hazard"])
    )
    counts.to_csv(args.output_dir / f"Fig4_gcm_boxplot_counts_{suffix}.csv", index=False, encoding="utf-8-sig")
    print(f"Wrote {out_prefix.with_suffix('.png')}")
    print(f"Wrote {out_prefix.with_suffix('.pdf')}")
    print(f"Wrote {out_prefix.with_suffix('.svg')}")
    print(f"Wrote source tables to {args.output_dir}")
    print(counts.to_string(index=False))


if __name__ == "__main__":
    main()
