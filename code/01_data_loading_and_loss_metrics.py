# Public-release script split from the original analysis notebook.
# Paths remain project-specific and may need to be edited before rerunning.
# Part 1: input readers, database assembly, event occurrence, and yield-loss metrics.

# Extracted from D:\AAUDE\paper\paper4\revision_v2\code\analysis.ipynb
# Generated for read-only public-release audit. Review before reuse.

# %% cell 1
import os
import xarray as xr
import os
import xarray as xr
import rioxarray as rxr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["pdf.fonttype"] = 42   # TrueType
mpl.rcParams["ps.fonttype"]  = 42   # EPS/PS 也一起
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]  # 或 ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

def read_crop(fp,crop, subdir, fname_template, varname):
    """
    crop: 'maiz'/'soyb'/'rice'/'wwh'/'swh'/'whe'
    subdir: 'waterstress' or 'nowaterstress'
    fname_template: e.g. 'gs_loss_nonIrr_annualtot_{crop}.nc'
    varname: e.g. 'gs_loss_nonIrr_annualtot'
    """
    if crop != "whea":
        path = os.path.join(fp, subdir, fname_template.format(crop=crop))
        return open_da(path, varname)

    # whe = wwh + swh
    path_wwh = os.path.join(fp, subdir, fname_template.format(crop="wwh"))
    path_swh = os.path.join(fp, subdir, fname_template.format(crop="swh"))

    da_wwh = open_da(path_wwh, varname)
    da_swh = open_da(path_swh, varname)
    da=da_wwh.fillna(0) + da_swh.fillna(0)
    return da.where(da!=0)
def potential_yield(base_dir, ssp, year, crops, modes):
    """
    返回 yield_out[crop][mode] (interp到cellarea网格，单位保持输入数据的单位)
    wheat = wwhe + swhe
    rice 文件前缀 = ricw
    """
    cellarea = xr.open_dataset(r'D:\Edge_download\cellarea.nc').cellarea30min

    yield_out = {}
    for crop in crops:
        yield_out[crop] = {}
        for mode in modes:
            if crop == "rice":
                fn = f"ricw_{mode}_{ssp}_{year}.nc"
                fp = os.path.join(base_dir, fn)
                da = xr.open_dataset(fp)["yield"]

            elif crop == "whea":
                fn1 = f"wwhe_{mode}_{ssp}_{year}.nc"
                fn2 = f"swhe_{mode}_{ssp}_{year}.nc"
                fp1 = os.path.join(base_dir, fn1)
                fp2 = os.path.join(base_dir, fn2)
                da1 = xr.open_dataset(fp1)["yield"]
                da2 = xr.open_dataset(fp2)["yield"]
                da = da1.fillna(0) + da2.fillna(0)

            else:
                fn = f"{crop}_{mode}_{ssp}_{year}.nc"
                fp = os.path.join(base_dir, fn)
                da = xr.open_dataset(fp)["yield"]

            da_i = da.interp(lat=cellarea.lat, lon=cellarea.lon, method="linear")
            yield_out[crop][mode] = da_i.where(da_i != 0)

    return yield_out

def open_da(path, varname):
    ds = xr.open_dataset(path, decode_times=False)
    da = ds[varname]

    t = da["time"].values.astype(float)
    units = ds["time"].attrs.get("units", "").lower()

    if "since" in units and "year" in units:
        base = units.split("since")[1].strip()
        base_year = int(base[:4])

        if np.nanmin(t) < 1000:
            year = (base_year + np.rint(t)).astype(int)
        else:
            year = np.rint(t).astype(int)

        da = da.assign_coords(time=("time", year)).rename({"time": "year"})
    else:
        da = da.rename({"time": "year"})
    return da.where(da!=0)
def distribuction():
    AR = r"F:\paper4\data\frac\30-arcminute_fraction"
    cellarea=xr.open_dataset(r'D:\Edge_download\cellarea.nc').cellarea30min 
    tif_map = {
        "wwh_Irr":     r"MIRCA-OS_Wheat_2000_irr_30arcmin_winter_aligned.tif",
        "swh_Irr":     r"MIRCA-OS_Wheat_2000_irr_30arcmin_spring_aligned.tif",
        "rice_Irr":    r"MIRCA-OS_Rice_2000_ir_30arcmin_aligned.tif",
        "maiz_Irr":    r"MIRCA-OS_Maize_2000_ir_30arcmin_aligned.tif",
        "soyb_Irr":    r"MIRCA-OS_Soybeans_2000_ir_30arcmin_aligned.tif",
    
        "wwh_nonIrr":  r"MIRCA-OS_Wheat_2000_rf_30arcmin_winter_aligned.tif",
        "swh_nonIrr":  r"MIRCA-OS_Wheat_2000_rf_30arcmin_spring_aligned.tif",
        "rice_nonIrr": r"MIRCA-OS_Rice_2000_rf_30arcmin_aligned.tif",
        "maiz_nonIrr": r"MIRCA-OS_Maize_2000_rf_30arcmin_aligned.tif",
        "soyb_nonIrr": r"MIRCA-OS_Soybeans_2000_rf_30arcmin_aligned.tif",
    }
    
    def read_mirca_tif(fp, key, cellarea):
        da = rxr.open_rasterio(fp).squeeze("band", drop=True)
    
        da = da.rename({"y": "lat", "x": "lon"})
        da = da.interp(lat=cellarea.lat, lon=cellarea.lon, method="linear")
        # da=da.where(da>0.001)
        da=da* cellarea
        da.name = key
        return da
    
    # 批量读取
    mirca_out = {}
    for key, rel in tif_map.items():
        fp = os.path.join(AR, rel)
        mirca_out[key] = read_mirca_tif(fp, key, cellarea)
    
    # 如果你想合成一个 Dataset
    mirca_ds = xr.Dataset(mirca_out)
    whea_Irr = (
        mirca_ds["wwh_Irr"].fillna(0)
      + mirca_ds["swh_Irr"].fillna(0)
    )
    
    whea_nonIrr = (
        mirca_ds["wwh_nonIrr"].fillna(0)
      + mirca_ds["swh_nonIrr"].fillna(0)
    )
    
    whea_Irr.name = "whea_Irr"
    whea_nonIrr.name = "whea_nonIrr"
    mirca_ds = mirca_ds.assign(
        whea_Irr = whea_Irr.where(whea_Irr>0),
        whea_nonIrr = whea_nonIrr.where(whea_nonIrr>0)
    )
    
    return mirca_ds


def stress_day(fp):
    crops = ['maiz', 'soyb', 'rice', 'whea'] 
    stress = {}
    fp=fp
    for crop in crops:
        stress[crop] = {}
    
        # 有水胁迫
        stress[crop]['dry_nonIrr'] = read_crop(
            fp,crop, "nowaterstress", "drystress_nonIrr_annualtot_{crop}.nc", "drystress_nonIrr_annualtot"
        )
        stress[crop]['wet_nonIrr'] = read_crop(
            fp,crop, "nowaterstress", "wetstress_nonIrr_annualtot_{crop}.nc", "wetstress_nonIrr_annualtot"
        )
    
        # 无水胁迫
        stress[crop]['dry_Irr'] =  read_crop(
            fp,crop, "nowaterstress", "drystress_nonIrr_annualtot_{crop}.nc", "drystress_nonIrr_annualtot"
        )
        stress[crop]['wet_Irr'] = read_crop(
            fp,crop, "nowaterstress", "wetstress_nonIrr_annualtot_{crop}.nc", "wetstress_nonIrr_annualtot"
        )
    return stress

def realtive_yield(fp):
      
    data = {}
    
    for crop in crops:
        data[crop] = {}
    
        # 有水胁迫
        data[crop]['waterstress_nonIrr'] = 4-read_crop(
            fp,crop, "waterstress", "gs_loss_nonIrr_annualtot_{crop}.nc", "gs_loss_nonIrr_annualtot"
        )
        data[crop]['waterstress_Irr'] = 4-read_crop(
            fp,crop, "waterstress", "gs_loss_Irr_annualtot_{crop}.nc", "gs_loss_Irr_annualtot"
        )
    
        # 无水胁迫
        data[crop]['nowaterstress_nonIrr'] =  4-read_crop(
            fp,crop, "nowaterstress", "gs_loss_nonIrr_annualtot_{crop}.nc", "gs_loss_nonIrr_annualtot"
        )
        data[crop]['nowaterstress_Irr'] = 4-read_crop(
            fp,crop, "nowaterstress", "gs_loss_Irr_annualtot_{crop}.nc", "gs_loss_Irr_annualtot"
        )
    return data
def count_gs(fp):
    def read_crop(fp,crop, subdir, fname_template, varname):

        if crop != "whea":
            path = os.path.join(fp, subdir, fname_template.format(crop=crop))
            return xr.open_dataset(path)[varname]

        path_wwh = os.path.join(fp, subdir, fname_template.format(crop="wwh"))
        path_swh = os.path.join(fp, subdir, fname_template.format(crop="swh"))
    
        da_wwh = xr.open_dataset(path_wwh)[varname]
        da_swh = xr.open_dataset(path_swh)[varname]
        da=da_wwh.fillna(0) + da_swh.fillna(0)
        return da.where(da!=0)
    data = {}
    for crop in crops:
    
        data[crop]= read_crop(
            fp,crop, "", "{crop}_irr_gs_nodes.nc", "growing_season_length"
        )
    return data

# %% cell 2
import numpy as np
import xarray as xr
from scipy import stats
import matplotlib.pyplot as plt
color_map = {
    ("nonIrr", "dry"): "#8cbf87",
    ("nonIrr", "wet"): "#cb9475",
    ("Irr",    "dry"): "#3e608d",
    ("Irr",    "wet"): "#8d2f25",
}
crops = ["maiz", "soyb", "rice", "whea"]
crop_labels = {"maiz":"Maize","soyb":"Soybean","rice":"Rice","whea":"Wheat"}

def _get_area(mirca_ds, crop, mode):
    # 你自己已有这个函数就删掉这段，保持一致
    key = f"{crop}_{mode}"  # 例如 maiz_nonIrr / maiz_Irr
    return mirca_ds[key].fillna(0)

def loss_ts_onehaz(yr_ds, stress_ds, pot_ds, mirca_ds, crop, mode, haz,
                   use_frac=False, gslen_var=None, frac_th=0.01):
    """
    返回：year 维度的全球（或研究区）面积/潜在产量加权的损失百分比时间序列
    口径：num=(nws-ws) * A * Yp 仅在事件发生时计入；den=ws * A * Yp 全部计入
    """
    ws  = yr_ds[crop][f"waterstress_{mode}"][-30:]
    nws = yr_ds[crop][f"nowaterstress_{mode}"][-30:]
    A   = _get_area(mirca_ds, crop, mode)
    Yp  = pot_ds[crop][mode]

    ev = stress_ds[crop][f"{haz}_nonIrr"]  # days or fraction

    if use_frac:
        if gslen_var is None:
            frac = ev
        else:
            gslen = stress_ds[crop][gslen_var]
            frac = ev / gslen
        mask = frac > frac_th
    else:
        mask = ev > 0

    num = ((nws - ws).where(mask) * A * Yp).sum(["lat", "lon"])
    den = ((ws) * A * Yp).sum(["lat", "lon"])

    ts = (num / den) * 100.0

    # 统一成 year 维度
    if "year" not in ts.dims:
        ts = ts.groupby("time.year").mean("time")
    return ts
def summarize_all(yr_ds, stress_ds, pot_ds, mirca_ds,
                  modes=("nonIrr","Irr"),
                  use_frac=False, gslen_var_map=None, frac_th=0.01):
    """
    返回 dict[crop][mode][haz] = year series
    gslen_var_map: {"nonIrr":"xxx", "Irr":"yyy"} 或 None
    """
    out = {}
    for crop in crops:
        out[crop] = {}
        for mode in modes:
            gslen_var = None
            if gslen_var_map is not None:
                gslen_var = gslen_var_map.get(mode, None)

            ts_d = loss_ts_onehaz(yr_ds, stress_ds, pot_ds, mirca_ds, crop, mode, "dry",
                                  use_frac=use_frac, gslen_var=gslen_var, frac_th=frac_th)
            ts_w = loss_ts_onehaz(yr_ds, stress_ds, pot_ds, mirca_ds, crop, mode, "wet",
                                  use_frac=use_frac, gslen_var=gslen_var, frac_th=frac_th)

            ts_d, ts_w = xr.align(ts_d, ts_w, join="inner")
            out[crop][mode] = {"dry": ts_d, "wet": ts_w}
    return out
def p_to_star(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"

def plot_bar_crops_irrig_rainfed(summary, title="Mean yield losses (± interannual SD)"):
    """
    每个作物：两块（nonIrr/Irr），每块里两柱（dry/wet）
    """
    modes = ["nonIrr", "Irr"]
    hazs  = ["dry", "wet"]

    # 坐标布局：每个作物占一个“组宽”，组内有 nonIrr 与 Irr 两个子组，每个子组两根柱
    nC = len(crops)
    group_gap = 1.2
    sub_gap   = 0.42
    bar_w     = 0.18

    # 每根柱的 x 位置
    xs = []
    vals = []
    errs = []
    labels = []  # legend
    stars = {}   # (crop,mode) -> star for dry vs wet

    # 计算均值、SD、显著性（配对 t 检验：同 crop 同 mode 下 dry vs wet）
    stats_table = {}
    for crop in crops:
        stats_table[crop] = {}
        for mode in modes:
            d = summary[crop][mode]["dry"].values
            w = summary[crop][mode]["wet"].values
            ok = np.isfinite(d) & np.isfinite(w)
            if ok.sum() >= 3:
                p = stats.ttest_rel(d[ok], w[ok]).pvalue
                stars[(crop, mode)] = p_to_star(p)
            else:
                stars[(crop, mode)] = "ns"

            stats_table[crop][mode] = {
                "dry_mean": float(np.nanmean(d)),
                "dry_sd":   float(np.nanstd(d, ddof=1)),
                "wet_mean": float(np.nanmean(w)),
                "wet_sd":   float(np.nanstd(w, ddof=1)),
            }

    # 组中心
    crop_centers = np.arange(nC) * group_gap

    # 逐作物构建柱位置：nonIrr 左子组，Irr 右子组；每个子组 dry 左、wet 右
    for i, crop in enumerate(crops):
        base = crop_centers[i]

        for j, mode in enumerate(modes):
            sub_center = base + (-sub_gap/2 if mode=="nonIrr" else sub_gap/2)

            # dry / wet
            for k, haz in enumerate(hazs):
                x = sub_center + (-bar_w/2 if haz=="dry" else bar_w/2)
                xs.append(x)

                m = stats_table[crop][mode][f"{haz}_mean"]
                s = stats_table[crop][mode][f"{haz}_sd"]
                vals.append(m)
                errs.append(s)

                labels.append((mode, haz))

    # 画图
    fig, ax = plt.subplots(figsize=(8, 3.5))

    # 为 legend 固定顺序
    style_order = [("nonIrr","dry"), ("nonIrr","wet"), ("Irr","dry"), ("Irr","wet")]
    style_idx = {t:i for i,t in enumerate(style_order)}

    # 分开画四类柱，保证图例清晰
    for (mode, haz) in style_order:
        idx = [ii for ii,(m,h) in enumerate(labels) if (m,h)==(mode,haz)]
        ax.bar(
            np.array(xs)[idx],
            np.array(vals)[idx],
            width=bar_w,
            yerr=np.array(errs)[idx],
            capsize=3,
            color=color_map[(mode, haz)],
            edgecolor="black",
            linewidth=0.6,
            label=f"{mode}-{('Drought' if haz=='dry' else 'Waterlogging')}"
        )


    ax.set_xticks(crop_centers)
    ax.set_xticklabels([crop_labels[c] for c in crops])
    ax.set_ylabel("Yield loss (%)")
    # ax.set_title(title)
    ax.legend(frameon=False, ncol=2)

    # 画显著性括号：同一作物同一模式下 dry vs wet
    for i, crop in enumerate(crops):
        base = crop_centers[i]
        for mode in modes:
            sub_center = base + (-sub_gap/2 if mode=="nonIrr" else sub_gap/2)
            x1 = sub_center - bar_w/2  # dry
            x2 = sub_center + bar_w/2  # wet

            # 找到这两根柱的高度+误差
            m_d = stats_table[crop][mode]["dry_mean"]
            s_d = stats_table[crop][mode]["dry_sd"]
            m_w = stats_table[crop][mode]["wet_mean"]
            s_w = stats_table[crop][mode]["wet_sd"]

            ymax = max(m_d + s_d, m_w + s_w)
            y = ymax + 0.25
            h = 0.15

            ax.plot([x1, x1, x2, x2], [y, y+h, y+h, y], c="k", lw=0.9)
            ax.text((x1+x2)/2, y+h+0.03, stars[(crop,mode)], ha="center", va="bottom", fontsize=9)

            # 在子组下方标注 nonIrr / Irr（可选）
        # ax.text(base - sub_gap/2, -0.12, "nonIrr", ha="center", va="top", transform=ax.get_xaxis_transform())
        # ax.text(base + sub_gap/2, -0.12, "Irr",    ha="center", va="top", transform=ax.get_xaxis_transform())

    ax.set_ylim(bottom=0)
    plt.savefig(
        rf"D:\AAUDE\paper\paper4\revision_v2\fig\manu_fig\hist_allcrops_loss.pdf",
        dpi=300,
        bbox_inches="tight"
    )
    plt.tight_layout()
    return fig, ax
# scen='hist'
# summary_hist = summarize_all(
#     yr_ds=DB[scen]["relative_yield"],
#     stress_ds=DB[scen]["stress_days"],
#     pot_ds=DB[scen]["potential_yield"],
#     mirca_ds=DB["dist"],
#     modes=("nonIrr","Irr")
# )

# fig, ax = plot_bar_crops_irrig_rainfed(summary_hist, title="Historical mean yield losses (± interannual SD)")
# plt.show()

# %% cell 3
def build_db_all_gcms(crops, modes, gcms, test):
    """
    返回 DB：
      DB['dist']                           # mirca 分布（所有情景共用）
      DB['hist'][...]                      # 历史（共用）
      DB['future'][gcm][scen][...]         # 未来：每个 gcm 下的 ssp126/ssp585

    每个场景下包含：
      ['relative_yield', 'stress_days', 'length', 'potential_yield']
    """
    DB = {}

    # 1) 作物分布（所有情景共用）
    DB["dist"] = distribuction()

    # 2) 历史（共用，不依赖 GCM）
    hist_cfg = {
        "stress_fp": rf"F:\paper4\yie\{test}",
        "yield_fp":  rf"F:\paper4\yie\{test}",
        "gs_fp":     r"F:\paper4\data\calendar\after_clip",
        "pot_dir":   r"F:\GAEZ_potential_yield\Hist_CRUTS32_nc",
        "pot_ssp":   "Hist",
        "pot_year":  "CRUTS32_8110L",
    }

    DB["hist"] = {
        "relative_yield": realtive_yield(hist_cfg["yield_fp"]),
        "stress_days":    stress_day(hist_cfg["stress_fp"]),
        "length":         count_gs(hist_cfg["gs_fp"]),
        "potential_yield": potential_yield(
            hist_cfg["pot_dir"], hist_cfg["pot_ssp"], hist_cfg["pot_year"],
            crops=crops, modes=modes
        ),
    }

    # 3) 未来：所有 GCM
    DB["future"] = {}

    for gcm in gcms:
        print(f"Building FUTURE DB for {gcm}")
        DB["future"][gcm] = {}

        scen_cfg = {
            "ssp126": {
                "stress_fp": rf"F:\paper4\yie\{gcm}\ssp126",
                "yield_fp":  rf"F:\paper4\yie\{gcm}\ssp126",
                "gs_fp":     r"F:\paper4\data\calendar\ssp126",
                "pot_dir":   r"F:\GAEZ_potential_yield\global_mean_nc",
                "pot_ssp":   "ssp126",
                "pot_year":  "2080sH",
            },
            "ssp585": {
                "stress_fp": rf"F:\paper4\yie\{gcm}\ssp585",
                "yield_fp":  rf"F:\paper4\yie\{gcm}\ssp585",
                "gs_fp":     r"F:\paper4\data\calendar\ssp585",
                "pot_dir":   r"F:\GAEZ_potential_yield\global_mean_nc",
                "pot_ssp":   "ssp585",
                "pot_year":  "2080sH",
            },
        }

        for scen, cfg in scen_cfg.items():
            DB["future"][gcm][scen] = {
                "relative_yield": realtive_yield(cfg["yield_fp"]),
                "stress_days":    stress_day(cfg["stress_fp"]),
                "length":         count_gs(cfg["gs_fp"]),
                "potential_yield": potential_yield(
                    cfg["pot_dir"], cfg["pot_ssp"], cfg["pot_year"],
                    crops=crops, modes=modes
                ),
            }

    return DB
crops = ["maiz", "soyb", "rice", "whea"]
modes = ["Irr", "nonIrr"]
test  = "gswp3_long"
gcms = [
    "CanESM5",
    "CNRM-CM6-1",
    "CNRM-ESM2-1",
    "EC-Earth3",
    "GFDL-ESM4",
    "IPSL-CM6A-LR",
    "MIROC6",
    "MPI-ESM1-2-HR",
    "MRI-ESM2-0",
    "UKESM1-0-LL",
]

DB = build_db_all_gcms(crops=crops, modes=modes, gcms=gcms, test=test)

# %% cell 4
import os
import numpy as np
import pandas as pd
import xarray as xr

OUT_DIR = r"D:\AAUDE\paper_v2\paper4\tables_from_DB"
os.makedirs(OUT_DIR, exist_ok=True)

# ----------------------------
# Helpers
# ----------------------------
def _area_da(DB, crop, mode):
    """MIRCA harvested area (already multiplied by cellarea in your distribuction())"""
    key = f"{crop}_{mode}"
    if key not in DB["dist"]:
        raise KeyError(f"Missing area key in DB['dist']: {key}")
    return DB["dist"][key].fillna(0)

def _align_like(a, b):
    """Align two DataArrays to the same lat/lon/year coords"""
    a2, b2 = xr.align(a, b, join="inner")
    return a2, b2

def _year_dim(da):
    if "year" in da.dims: return "year"
    if "time" in da.dims: return "time"
    raise ValueError(f"Cannot find year/time dim in {da.dims}")

def _bool_hazard(stress_days_da, thr=0):
    """hazard occurs if stress days > thr"""
    return (stress_days_da > thr)

def _safe_mean_days_on_haz_years(stress_days, haz_bool):
    """
    mean stress days only on hazard years, area-weighted over space, then mean over years
    """
    # mask non-hazard
    sd = stress_days.where(haz_bool)
    return sd

# ----------------------------
# Core: compute occurrence stats for one dataset (hist or one future scenario)
# ----------------------------
def occurrence_stats_one(DB_block, DB_dist, crop, mode, haz, thr=0):
    """
    DB_block: e.g. DB["hist"] or DB["future"][gcm][scen]
    haz: "dry" or "wet"
    returns dict of scalar stats
    """
    area = _area_da({"dist": DB_dist}, crop, mode)  # reuse helper with dist only

    key = f"{haz}_{mode}"
    if key not in DB_block["stress_days"][crop]:
        raise KeyError(f"Missing stress_days[{crop}][{key}]")

    sd = DB_block["stress_days"][crop][key].fillna(0)
    ydim = _year_dim(sd)

    # align
    area2, sd2 = _align_like(area, sd)

    # consider only cells with crop area > 0
    mask_crop = (area2 > 0)

    haz_bool = _bool_hazard(sd2, thr=thr) & mask_crop

    # (1) grid-year fraction: count of (cell,year) where hazard occurs / all (cell,year) with crop
    n_haz = haz_bool.sum(dim=("lat","lon"))
    n_all = mask_crop.sum(dim=("lat","lon"))
    gridyear_frac = (n_haz / n_all).mean(dim=ydim).item()

    # (2) area fraction per year: affected area / total area, then mean over years
    affected_area = (area2.where(haz_bool, 0.0)).sum(dim=("lat","lon"))
    total_area    = (area2.where(mask_crop, 0.0)).sum(dim=("lat","lon"))
    area_frac_mean = (affected_area / total_area).mean(dim=ydim).item()

    # (3) mean stress days on hazard-years, area-weighted over space, then mean over years
    sd_haz = sd2.where(haz_bool)
    # area-weighted mean days per year:
    mean_days_year = (sd_haz * area2).sum(dim=("lat","lon")) / (area2.where(haz_bool).sum(dim=("lat","lon")))
    mean_days_on_haz_years = float(mean_days_year.mean(dim=ydim).values)

    return dict(
        gridyear_frac=gridyear_frac,
        area_frac_mean=area_frac_mean,
        mean_days_on_haz_years=mean_days_on_haz_years
    )

# ----------------------------
# Run: HIST + FUTURE(all GCMs & scenarios)
# ----------------------------
def export_occurrence_tables(DB, crops, modes, gcms, scenarios=("ssp126","ssp585"), thr=0):
    rows = []

    # HIST
    for crop in crops:
        for mode in modes:
            for haz in ("dry","wet"):
                st = occurrence_stats_one(DB["hist"], DB["dist"], crop, mode, haz, thr=thr)
                rows.append(dict(period="hist", gcm="NA", scen="NA", crop=crop, mode=mode, haz=haz, **st))

    # FUTURE
    for gcm in gcms:
        for scen in scenarios:
            block = DB["future"][gcm][scen]
            for crop in crops:
                for mode in modes:
                    for haz in ("dry","wet"):
                        st = occurrence_stats_one(block, DB["dist"], crop, mode, haz, thr=thr)
                        rows.append(dict(period="future", gcm=gcm, scen=scen, crop=crop, mode=mode, haz=haz, **st))

    df = pd.DataFrame(rows)
    out_fp = os.path.join(OUT_DIR, "occurrence_stats_P0.csv")
    df.to_csv(out_fp, index=False)
    print("Saved:", out_fp)
    return df

# run
df_occ = export_occurrence_tables(DB, crops=crops, modes=modes, gcms=gcms, thr=7)
df_occ.head()

# %% cell 5
import xarray as xr

def _pick_block(DB, scen, gcm=None):
    """
    返回某个场景对应的数据块 block（包含 relative_yield / stress_days / length / potential_yield）
    兼容两种 DB 结构：
      1) 旧：DB[scen]
      2) 新：DB["hist"] 以及 DB["future"][gcm][scen]
    """
    # 旧结构：DB 直接以 scen 为键
    if scen in DB and isinstance(DB[scen], dict) and "relative_yield" in DB[scen]:
        return DB[scen]

    # 新结构：hist
    if scen == "hist":
        if "hist" in DB and isinstance(DB["hist"], dict):
            return DB["hist"]
        raise KeyError("DB 中找不到历史场景：需要 DB['hist'] 或 DB['hist'](旧结构)")

    # 新结构：future
    if "future" not in DB:
        raise KeyError("DB 中没有 'future'，但你请求了未来情景。")

    if gcm is None:
        raise ValueError("未来情景必须提供 gcm（单条）或让主函数走 ensemble 分支（gcm=None）。")

    if gcm not in DB["future"]:
        raise KeyError(f"DB['future'] 中找不到 gcm={gcm}")

    if scen not in DB["future"][gcm]:
        raise KeyError(f"DB['future'][{gcm}] 中找不到 scen={scen}")

    return DB["future"][gcm][scen]


def hazard_loss_total(DB, mirca_ds, scen="hist", crop="maiz", haz="dry",
                      denom="ws", frac_th=None, use_gslen=False, gcm=None):
    """
    兼容未来 ensemble 的版本。

    参数新增：
      gcm:
        - hist：可忽略
        - future：若提供字符串，则返回该 gcm 单条曲线
                 若为 None，则返回所有 gcm concat 后的 DataArray（带 gcm 维）

    返回：
      - hist：DataArray(year/time)  一条
      - future + gcm=None：DataArray(gcm, year/time) 多条
    """
    # ==========
    # 未来 ensemble：如果请求未来 scen 且 gcm=None，则对所有 gcm 循环计算并 concat
    # ==========
    is_future = (scen != "hist") and ("future" in DB) and (scen not in DB)  # 尽量不误伤旧结构
    if is_future and gcm is None:
        ts_list = []
        for gg in DB["future"].keys():
            ts = hazard_loss_total(
                DB, mirca_ds, scen=scen, crop=crop, haz=haz,
                denom=denom, frac_th=frac_th, use_gslen=use_gslen, gcm=gg
            )
            ts = ts.expand_dims(gcm=[gg])
            ts_list.append(ts)
        return xr.concat(ts_list, dim="gcm")

    # ==========
    # 单场景单块（hist 或 指定 gcm 的 future）
    # ==========
    block = _pick_block(DB, scen=scen, gcm=gcm)

    # --- 面积（ha）: 你 mirca_ds 是 m2
    A_rf_ha = mirca_ds[f"{crop}_nonIrr"] / 10000.0
    A_ir_ha = mirca_ds[f"{crop}_Irr"] / 10000.0

    # --- 潜在产量（t/ha）
    Yp_rf = block["potential_yield"][crop]["nonIrr"]
    Yp_ir = block["potential_yield"][crop]["Irr"]

    pot_rf = Yp_rf * A_rf_ha   # t
    pot_ir = Yp_ir * A_ir_ha   # t

    # --- 相对产量
    ws_rf  = block["relative_yield"][crop]["waterstress_nonIrr"]
    nws_rf = block["relative_yield"][crop]["nowaterstress_nonIrr"]
    ws_ir  = block["relative_yield"][crop]["waterstress_Irr"]
    nws_ir = block["relative_yield"][crop]["nowaterstress_Irr"]

    # --- 事件掩膜
    ev_rf = block["stress_days"][crop][f"{haz}_nonIrr"]
    ev_ir = block["stress_days"][crop][f"{haz}_Irr"]

    if frac_th is not None:
        if not use_gslen:
            mask_rf = ev_rf > frac_th
            mask_ir = ev_ir > frac_th
        else:
            gslen = block["length"][crop]  # 单位需要你确认：天数
            mask_rf = (ev_rf / gslen) > frac_th
            mask_ir = (ev_ir / gslen) > frac_th
    else:
        mask_rf = ev_rf > 0
        mask_ir = ev_ir > 0

    # --- 损失量（t）
    loss_rf = ((nws_rf - ws_rf).where(mask_rf) * pot_rf).sum(("lat", "lon"))
    loss_ir = ((nws_ir - ws_ir).where(mask_ir) * pot_ir).sum(("lat", "lon"))
    loss_total = loss_rf + loss_ir

    # --- 分母（t）
    if denom == "ws":
        den = (ws_rf * pot_rf).sum(("lat", "lon")) + (ws_ir * pot_ir).sum(("lat", "lon"))
    elif denom == "nws":
        den = (nws_rf * pot_rf).sum(("lat", "lon")) + (nws_ir * pot_ir).sum(("lat", "lon"))
    else:
        raise ValueError("denom must be 'ws' or 'nws'")

    return (loss_total / den) * 100.0  # %

# ==========================
# 画图：dry / wet / (dry+wet)
# ==========================
scen = "hist"
crop = "maiz"
mirca_ds=DB['dist']
ts_dry = hazard_loss_total(DB, mirca_ds, scen=scen, crop=crop, haz="dry", denom="ws")
ts_wet = hazard_loss_total(DB, mirca_ds, scen=scen, crop=crop, haz="wet", denom="ws")
ts_dry, ts_wet = xr.align(ts_dry, ts_wet, join="inner")
ts_all = ts_dry + ts_wet

plt.figure(figsize=(9,4))
ts_dry.plot(label="Drought (dry) - total (rainfed+irrigated)")
ts_wet.plot(label="Waterlogging (wet) - total (rainfed+irrigated)")
ts_all.plot(label="Total (dry+wet)")
plt.ylabel("Yield loss (%)")
plt.title(f"{scen} | {crop} yield loss (%)")
plt.legend(frameon=False)
plt.tight_layout()
plt.show()

# %% cell 6
scen='ssp126'
summary_hist = summarize_all(
    yr_ds=DB['future']['CanESM5'][scen]["relative_yield"],
    stress_ds=DB['future']['CanESM5'][scen]["stress_days"],
    pot_ds=DB['future']['CanESM5'][scen]["potential_yield"],
    mirca_ds=DB["dist"],
    modes=("nonIrr","Irr")
)

fig, ax = plot_bar_crops_irrig_rainfed(summary_hist, title="Historical mean yield losses (± interannual SD)")
plt.savefig(rf'D:\AAUDE\paper_v2\paper4\outputs_loss_maps\{scen}.pdf')
plt.show()

# %% cell 7
scen='ssp585'
summary_hist = summarize_all(
    yr_ds=DB['future']['CanESM5'][scen]["relative_yield"],
    stress_ds=DB['future']['CanESM5'][scen]["stress_days"],
    pot_ds=DB['future']['CanESM5'][scen]["potential_yield"],
    mirca_ds=DB["dist"],
    modes=("nonIrr","Irr")
)

fig, ax = plot_bar_crops_irrig_rainfed(summary_hist, title="Historical mean yield losses (± interannual SD)")
plt.savefig(rf'D:\AAUDE\paper_v2\paper4\outputs_loss_maps\{scen}.pdf')
plt.show()

# %% cell 8
scen='hist'
summary_hist = summarize_all(
    yr_ds=DB[scen]["relative_yield"],
    stress_ds=DB[scen]["stress_days"],
    pot_ds=DB[scen]["potential_yield"],
    mirca_ds=DB["dist"],
    modes=("nonIrr","Irr")
)

fig, ax = plot_bar_crops_irrig_rainfed(summary_hist, title="Historical mean yield losses (± interannual SD)")
plt.show()

