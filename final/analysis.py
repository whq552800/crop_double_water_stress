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

# %% cell 9
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
for scen in ['hist','ssp126','ssp585']:
# scen = "hist"
    crops = ["maiz", "soyb", "rice", "whea"]
    
    num_sum = None   # 分子累计
    den_sum = None   # 分母累计
    
    for crop in crops:
        if scen!='hist':
            wet = DB['future'][scen]["stress_days"][crop]["wet_Irr"][-30:].median("year")
            dry = DB['future'][scen]["stress_days"][crop]["dry_Irr"][-30:].median("year")
            diff = wet - dry   # 天
        
            # --- (2) 生长季长度（天）
            gslen = DB['future'][scen]["length"][crop]*2
        else:
            wet = DB[scen]["stress_days"][crop]["wet_Irr"][-30:].median("year")
            dry = DB[scen]["stress_days"][crop]["dry_Irr"][-30:].median("year")
            diff = wet - dry   # 天
        
            # --- (2) 生长季长度（天）
            gslen = DB[scen]["length"][crop]*2
    
        # --- (3) 作物面积权重（ha）
        # A = mirca_ds[f"{crop}_Irr"] / 10000.0
        A = 1
    
        # --- 加权
        num = (diff * A).fillna(0)
        den = (gslen * A).fillna(0)
    
        num_sum = num if num_sum is None else num_sum + num
        den_sum = den if den_sum is None else den_sum + den
    
    # --- 合并后的比例
    ratio_all = num_sum / den_sum
    ratio_all = ratio_all.where(den_sum != 0)
    fig = plt.figure(figsize=(9, 4.5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    # 空间范围
    ax.set_extent([-130, 160, -60, 75], crs=ccrs.PlateCarree())
    ax.set_aspect(1.3)
    # 主图
    import matplotlib.colors as mcolors
    
    # 发散色带 + 以 0 为中心
    cmap = plt.get_cmap("BrBG")
    norm = mcolors.TwoSlopeNorm(vmin=-0.15, vcenter=0.0, vmax=0.15)
    
    im = ratio_all.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        x="lon", y="lat",
        cmap=cmap,
        norm=norm,
        add_colorbar=True,
        cbar_kwargs={
            "label": "Relative dominance (%)",
            "ticks": [-0.15, -0.1,-0.05, 0, 0.05,0.1, 0.15],
            "shrink": 0.85,
            "pad": 0.02,
            "extend": "neither"   # ✅ 关键这一行
        }
    )
    
    
    # 海岸线
    ax.coastlines(resolution="110m", linewidth=0.8)
    
    # 陆地底色（可选）
    ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=0)
    
    # ❗ 去掉外边框（关键）
    ax.set_frame_on(False)
    ax.patch.set_visible(False)
    if hasattr(ax, "spines"):
        for spine in ax.spines.values():
            spine.set_visible(False)
    plt.title(scen)
    plt.tight_layout()
    plt.show()

# %% cell 10
crop='whea'
A_rf = mirca_ds[f"{crop}_nonIrr"] / 10000.0
A_ir = mirca_ds[f"{crop}_Irr"]    / 10000.0
haz='wet'
for gcm in gcms:
    # --- 潜在产量 (t/ha)
    Yp_rf = DB['future'][gcm]['ssp126']["potential_yield"][crop]["nonIrr"]
    Yp_ir = DB['future'][gcm]['ssp126']["potential_yield"][crop]["Irr"]
    
    pot_rf = Yp_rf * A_rf
    pot_ir = Yp_ir * A_ir
    
    # --- 相对产量（历史全时段）
    ws_rf  = DB['future'][gcm]['ssp126']["relative_yield"][crop]["waterstress_nonIrr"]
    nws_rf = DB['future'][gcm]['ssp126']["relative_yield"][crop]["nowaterstress_nonIrr"]
    ws_ir  = DB['future'][gcm]['ssp126']["relative_yield"][crop]["waterstress_Irr"]
    nws_ir = DB['future'][gcm]['ssp126']["relative_yield"][crop]["nowaterstress_Irr"]
    mask_rf = DB['future'][gcm]['ssp126']["stress_days"][crop][f"{haz}_nonIrr"] > 0
    mask_ir = DB['future'][gcm]['ssp126']["stress_days"][crop][f"{haz}_Irr"] > 0
    den = nws_rf * pot_rf.fillna(0) + nws_ir * pot_ir.fillna(0)
    den = den.where(den!=0)
    loss_rf = (nws_rf - ws_rf).where(mask_rf) * pot_rf
    loss_ir = (nws_ir - ws_ir).where(mask_ir) * pot_ir
    loss = loss_rf.where(loss_rf>0).fillna(0) + loss_ir.where(loss_ir>0).fillna(0)
    loss=loss.where(loss!=0)
    (loss.sum(('lat','lon'))/den.sum(('lat','lon'))*100).plot()
plt.show()

# %% cell 11
print("year equal:", a["year"].equals(b["year"]), a["year"].dtype, b["year"].dtype)

t=10
print("a nnz:", int((a[t]!=0).sum()), "b nnz:", int((b[t]!=0).sum()), "den nnz:", int((den[t]!=0).sum()))
print("a min/max:", float(a.min()), float(a.max()))
print("b min/max:", float(b.min()), float(b.max()))

# %% cell 12
loss_sum = None
den_sum  = None
haz='wet'
crops=("maiz","soyb","rice","whea");mirca_ds=DB['dist'];eps=0
for crop in crops:
        # --- 面积 (ha)
    A_rf = mirca_ds[f"{crop}_nonIrr"] / 10000.0
    A_ir = mirca_ds[f"{crop}_Irr"]    / 10000.0

    # --- 潜在产量 (t/ha)
    Yp_rf = DB["hist"]["potential_yield"][crop]["nonIrr"]
    Yp_ir = DB["hist"]["potential_yield"][crop]["Irr"]

    pot_rf = Yp_rf * A_rf
    pot_ir = Yp_ir * A_ir

    # --- 相对产量（历史全时段）
    ws_rf  = DB["hist"]["relative_yield"][crop]["waterstress_nonIrr"]
    nws_rf = DB["hist"]["relative_yield"][crop]["nowaterstress_nonIrr"]
    ws_ir  = DB["hist"]["relative_yield"][crop]["waterstress_Irr"]
    nws_ir = DB["hist"]["relative_yield"][crop]["nowaterstress_Irr"]

    # --- 事件掩膜
    mask_rf = DB["hist"]["stress_days"][crop][f"{haz}_nonIrr"] > 0
    mask_ir = DB["hist"]["stress_days"][crop][f"{haz}_Irr"] > 0

    loss_rf = (nws_rf - ws_rf).where(mask_rf) * pot_rf
    loss_ir = (nws_ir - ws_ir).where(mask_ir) * pot_ir
    loss = loss_rf.where(loss_rf>0).fillna(0) + loss_ir.where(loss_ir>0).fillna(0)
    a = (nws_rf * pot_rf).fillna(0)
    b = (nws_ir * pot_ir).fillna(0) 
    a2, b2 = xr.align(a, b, join="outer", fill_value=0)
    den = nws_rf * pot_rf.fillna(0) + nws_ir * pot_ir.fillna(0)
    # den = den.fillna(0)

    loss_sum = loss if loss_sum is None else loss_sum + loss
    den_sum  = den  if den_sum  is None else den_sum  + den

loss_pct = (loss_sum / (den_sum + eps)) * 100.0

loss_pct=loss_pct.where(loss_pct!=0)

# %% cell 13

def spatial_loss_map_allcrops_hist(
    DB, mirca_ds, haz,
    crops=("maiz","soyb","rice","whea"),
    denom="ws", eps=1e-12
):
    """
    返回：DataArray(lat,lon)
    历史期（1980–2019）四作物合并后的空间产量损失百分比（%）
    """
    loss_sum = None
    den_sum  = None

    for crop in crops:
        # --- 面积 (ha)
        A_rf = mirca_ds[f"{crop}_nonIrr"] / 10000.0
        A_ir = mirca_ds[f"{crop}_Irr"]    / 10000.0

        # --- 潜在产量 (t/ha)
        Yp_rf = DB["hist"]["potential_yield"][crop]["nonIrr"]
        Yp_ir = DB["hist"]["potential_yield"][crop]["Irr"]

        pot_rf = Yp_rf * A_rf
        pot_ir = Yp_ir * A_ir

        # --- 相对产量（历史全时段）
        ws_rf  = DB["hist"]["relative_yield"][crop]["waterstress_nonIrr"]
        nws_rf = DB["hist"]["relative_yield"][crop]["nowaterstress_nonIrr"]
        ws_ir  = DB["hist"]["relative_yield"][crop]["waterstress_Irr"]
        nws_ir = DB["hist"]["relative_yield"][crop]["nowaterstress_Irr"]
        length=DB["hist"]['length'][crop]
        # --- 事件掩膜
        mask_rf = (DB["hist"]["stress_days"][crop][f"{haz}_nonIrr"]/length) > 0.05
        mask_ir = (DB["hist"]["stress_days"][crop][f"{haz}_Irr"]/length) > 0.05

        loss_rf = (nws_rf - ws_rf).where(mask_rf) * pot_rf
        loss_ir = (nws_ir - ws_ir).where(mask_ir) * pot_ir
        loss = loss_rf.where(loss_rf>0).fillna(0) + loss_ir.where(loss_ir>0).fillna(0)
        a = (nws_rf * pot_rf).fillna(0)
        b = (nws_ir * pot_ir).fillna(0)
        
        a2, b2 = xr.align(a, b, join="outer", fill_value=0)
        den = a2 + b2

        loss_sum = loss if loss_sum is None else loss_sum + loss
        den_sum  = den  if den_sum  is None else den_sum  + den

    loss_pct = (loss_sum / (den_sum + eps)) * 100.0
    return loss_pct.where(den_sum > 0)
year_sel_hist = slice(1980, 2019)

da_hist_dry = spatial_loss_map_allcrops_hist(
    DB, mirca_ds, haz="dry", denom="ws"
)

da_map_dry = reduce_time(da_hist_dry, year_slice=year_sel_hist, how="median")
da_hist_wet = spatial_loss_map_allcrops_hist(
    DB, mirca_ds, haz="wet", denom="ws"
)

da_map_wet = reduce_time(da_hist_wet, year_slice=year_sel_hist, how="median")

# %% cell 14
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
def draw(xx,nam):
    fig = plt.figure(figsize=(9, 4.5))
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    # 空间范围
    ax.set_extent([-130, 160, -60, 75], crs=ccrs.PlateCarree())
    ax.set_aspect(1.3)
    
    # 色带（如果是损失，用顺序色带更合理）
    cmap = plt.get_cmap("RdYlBu_r")
    norm = mcolors.Normalize(vmin=0, vmax=8)
    
    im = xx.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        x="lon", y="lat",
        cmap=cmap,
        norm=norm,
        add_colorbar=True,
        cbar_kwargs={
            "label": "Yield loss (%)",
            "shrink": 0.85,
            "pad": 0.02,
             "extend": "neither" 
        }
    )
    
    # 海岸线
    ax.coastlines(resolution="110m", linewidth=0.8)
    
    # 陆地底色
    ax.add_feature(cfeature.LAND, facecolor="whitesmoke", zorder=0)
    
    ax.set_frame_on(False)
    ax.patch.set_visible(False)
    if hasattr(ax, "spines"):
        for spine in ax.spines.values():
            spine.set_visible(False)
    
    plt.tight_layout()
    
    # plt.savefig(
    #     rf"D:\AAUDE\paper\paper4\revision_v2\fig\manu_fig\hist_allcrops_{nam}_loss.pdf",
    #     dpi=300,
    #     bbox_inches="tight"
    # )
    # plt.show()
draw(da_map_dry,'dry')
draw(da_map_wet,'wet')

# %% cell 15
year_sel = {
    "hist":   slice(1980, 2019),
    "ssp126": slice(2071, 2100),
    "ssp585": slice(2071, 2100),
}
mirca_ds=DB['dist']
def reduce_time(da, year_slice, how="mean"):
    """
    da: xarray DataArray, 必须有 year 维
    year_slice: slice(YYYY1, YYYY2)
    how: 'mean' | 'sum' | 'median'
    """
    if "year" not in da.dims:
        raise ValueError("DataArray must have 'year' dimension")
    da=da.where(da!=0)
    da_sel = da.sel(year=year_slice)

    if how == "mean":
        return da_sel.mean("year")
    elif how == "sum":
        return da_sel.sum("year")
    elif how == "median":
        return da_sel.median("year")
    elif how == "max":
        return da_sel.max("year")
    else:
        raise ValueError(f"Unknown how={how}")

# %% cell 16
import xarray as xr
import numpy as np

def spatial_loss_map_onecrop_hist(
    DB, mirca_ds, crop, haz, eps=1e-12
):
    """
    返回：DataArray(year, lat, lon) 或 (lat, lon) 取决于 DB 里 relative_yield 是否含 year
    历史期：单作物（crop）在 haz(dry/wet) 下的空间产量损失百分比（%）
    分子：只在事件发生处 (nws - ws) * pot
    分母：nws * pot （全期）
    """
    # --- 面积 (ha)
    A_rf = mirca_ds[f"{crop}_nonIrr"] / 10000.0
    A_ir = mirca_ds[f"{crop}_Irr"]    / 10000.0

    # --- 潜在产量 (t/ha)
    Yp_rf = DB["hist"]["potential_yield"][crop]["nonIrr"]
    Yp_ir = DB["hist"]["potential_yield"][crop]["Irr"]

    pot_rf = Yp_rf * A_rf
    pot_ir = Yp_ir * A_ir

    # --- 相对产量
    ws_rf  = DB["hist"]["relative_yield"][crop]["waterstress_nonIrr"]
    nws_rf = DB["hist"]["relative_yield"][crop]["nowaterstress_nonIrr"]
    ws_ir  = DB["hist"]["relative_yield"][crop]["waterstress_Irr"]
    nws_ir = DB["hist"]["relative_yield"][crop]["nowaterstress_Irr"]

    # --- 事件掩膜（days > 0）
    mask_rf = DB["hist"]["stress_days"][crop][f"{haz}_nonIrr"] > 0
    mask_ir = DB["hist"]["stress_days"][crop][f"{haz}_Irr"] > 0

    # --- 分子：只在事件处计入损失；负损失截断为 0
    loss_rf = (nws_rf - ws_rf).where(mask_rf) * pot_rf
    loss_ir = (nws_ir - ws_ir).where(mask_ir) * pot_ir
    loss = loss_rf.clip(min=0).fillna(0) + loss_ir.clip(min=0).fillna(0)

    # --- 分母：nws * pot（全期）
    a = (nws_rf * pot_rf).fillna(0)
    b = (nws_ir * pot_ir).fillna(0)
    a2, b2 = xr.align(a, b, join="outer", fill_value=0)
    den = a2 + b2

    loss_pct = (loss / (den + eps)) * 100.0
    return loss_pct.where(den > 0)
import matplotlib as mpl
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["axes.unicode_minus"] = False
# 若你想指定字体（可选）：
# mpl.rcParams["font.family"] = "serif"
# mpl.rcParams["font.serif"]  = ["Times New Roman"]

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

def draw_map(xx, crop, haz, out_dir):
    fig = plt.figure(figsize=(9, 4.5))
    ax = plt.axes(projection=ccrs.PlateCarree())

    ax.set_extent([-130, 160, -60, 75], crs=ccrs.PlateCarree())
    ax.set_aspect(1.3)

    cmap = plt.get_cmap("RdYlBu_r")
    norm = mcolors.Normalize(vmin=0, vmax=8)

    xx.plot(
        ax=ax,
        transform=ccrs.PlateCarree(),
        x="lon", y="lat",
        cmap=cmap,
        norm=norm,
        add_colorbar=True,
        cbar_kwargs={
            "label": "Yield loss (%)",
            "shrink": 0.85,
            "pad": 0.02,
            "extend": "neither"
        }
    )

    ax.coastlines(resolution="110m", linewidth=0.8)
    ax.add_feature(cfeature.LAND, facecolor="whitesmoke", zorder=0)

    ax.set_frame_on(False)
    ax.patch.set_visible(False)
    if hasattr(ax, "spines"):
        for spine in ax.spines.values():
            spine.set_visible(False)

    plt.tight_layout()

    fp = rf"{out_dir}\hist_{crop}_{haz}_loss.pdf"
    plt.savefig(fp, dpi=300, bbox_inches="tight")
    plt.show()

# %% cell 17
year_sel_hist = slice(1980, 2019)
out_dir = r"D:\AAUDE\paper\paper4\revision_v2\fig\manu_fig"

crops = ("maiz", "soyb", "rice", "whea")
hazs  = ("dry", "wet")

for crop in crops:
    for haz in hazs:
        da_hist = spatial_loss_map_onecrop_hist(DB, mirca_ds, crop=crop, haz=haz)
        da_map  = reduce_time(da_hist, year_slice=year_sel_hist, how="mean")
        draw_map(da_map, crop=crop, haz=haz, out_dir=out_dir)

# %% cell 18
import os
import numpy as np
import xarray as xr

import matplotlib as mpl
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["axes.unicode_minus"] = False

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def reduce_time(da, year_slice=slice(2070, 2100), how="mean"):
    """对 year 维度做聚合；如果没有 year 维则原样返回。"""
    if "year" not in da.dims:
        return da
    da2 = da.sel(year=year_slice)
    if how == "mean":
        return da2.mean("year", skipna=True)
    elif how == "sum":
        return da2.sum("year", skipna=True)
    else:
        raise ValueError("how must be 'mean' or 'sum'")


def _loss_pct_onecase(case_dict, mirca_ds, crop, haz, eps=1e-12):
    """
    case_dict: DB["hist"] 或 DB["future"][gcm][scen] 这一层
    返回：DataArray(year, lat, lon) 的 loss_pct
    """
    # --- 面积 (ha)
    A_rf = mirca_ds[f"{crop}_nonIrr"] / 10000.0
    A_ir = mirca_ds[f"{crop}_Irr"]    / 10000.0

    # --- 潜在产量 (t/ha)
    Yp_rf = case_dict["potential_yield"][crop]["nonIrr"]
    Yp_ir = case_dict["potential_yield"][crop]["Irr"]

    pot_rf = Yp_rf * A_rf
    pot_ir = Yp_ir * A_ir

    # --- 相对产量
    ws_rf  = case_dict["relative_yield"][crop]["waterstress_nonIrr"]
    nws_rf = case_dict["relative_yield"][crop]["nowaterstress_nonIrr"]
    ws_ir  = case_dict["relative_yield"][crop]["waterstress_Irr"]
    nws_ir = case_dict["relative_yield"][crop]["nowaterstress_Irr"]

    # --- 事件掩膜（days > 0）
    mask_rf = case_dict["stress_days"][crop][f"{haz}_nonIrr"] > 0
    mask_ir = case_dict["stress_days"][crop][f"{haz}_Irr"] > 0

    # --- 分子：只在事件处计入损失；负损失截断为 0
    loss_rf = (nws_rf - ws_rf).where(mask_rf) * pot_rf
    loss_ir = (nws_ir - ws_ir).where(mask_ir) * pot_ir
    loss = loss_rf.clip(min=0).fillna(0) + loss_ir.clip(min=0).fillna(0)

    # --- 分母：nws * pot（全期）
    a = (nws_rf * pot_rf).fillna(0)
    b = (nws_ir * pot_ir).fillna(0)
    a2, b2 = xr.align(a, b, join="outer", fill_value=0)
    den = a2 + b2

    loss_pct = (loss / (den + eps)) * 100.0
    return loss_pct.where(den > 0)


def spatial_loss_map_onecrop_future_ensmean(
    DB, mirca_ds, crop, haz, scen, gcms, year_slice=slice(2070, 2100), eps=1e-12
):
    """
    未来期：对所有 GCM 先逐个算 loss_pct，然后对 GCM 取均值；最后对 year_slice 求均值。
    返回：DataArray(lat, lon)
    """
    loss_list = []
    for gcm in gcms:
        case = DB["future"][gcm][scen]
        da = _loss_pct_onecase(case, mirca_ds, crop=crop, haz=haz, eps=eps)
        da = reduce_time(da, year_slice=year_slice, how="mean")  # -> (lat,lon)
        loss_list.append(da)

    # 对 GCM ensemble mean（等权）
    loss_stack = xr.concat(loss_list, dim=xr.IndexVariable("gcm", gcms))
    loss_ens = loss_stack.mean("gcm", skipna=True)
    return loss_ens


def draw_panel_4x2_maps(
    maps_dict,  # maps_dict[(crop,haz)] = DataArray(lat,lon)
    crops=("maiz","soyb","rice","whea"),
    hazs=("dry","wet"),
    scen="ssp585",
    out_dir=".",
    vmin=0, vmax=8,
    extent=(-130, 160, -60, 75),
    cmap_name="RdYlBu_r",
    dpi=300
):
    """
    4行x2列：行=作物，列=haz；共享一个colorbar。
    """
    fig, axes = plt.subplots(
        nrows=len(crops), ncols=len(hazs),
        figsize=(10.5, 12),
        subplot_kw={"projection": ccrs.PlateCarree()},
        constrained_layout=True
    )

    cmap = plt.get_cmap(cmap_name)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    mappable = None
    for i, crop in enumerate(crops):
        for j, haz in enumerate(hazs):
            ax = axes[i, j]
            ax.set_extent(extent, crs=ccrs.PlateCarree())
            ax.coastlines(resolution="110m", linewidth=0.6)
            ax.add_feature(cfeature.LAND, facecolor="whitesmoke", zorder=0)

            da = maps_dict[(crop, haz)]

            # 用 pcolormesh，避免每幅图单独colorbar
            pm = ax.pcolormesh(
                da["lon"].values,
                da["lat"].values,
                da.values,
                transform=ccrs.PlateCarree(),
                cmap=cmap,
                norm=norm,
                shading="auto"
            )
            mappable = pm  # 保存最后一个用于colorbar

            # 标题：第一行加列名；每行左侧加作物名
            if i == 0:
                ax.set_title(f"{haz}", fontsize=12)
            if j == 0:
                ax.text(
                    -0.02, 0.5, crop,
                    transform=ax.transAxes,
                    va="center", ha="right",
                    fontsize=12, rotation=90
                )

            # 去边框
            ax.set_frame_on(False)
            ax.patch.set_visible(False)
            if hasattr(ax, "spines"):
                for spine in ax.spines.values():
                    spine.set_visible(False)

    # 共享 colorbar
    cbar = fig.colorbar(
        mappable, ax=axes.ravel().tolist(),
        orientation="horizontal", shrink=0.9, pad=0.03
    )
    cbar.set_label("Yield loss (%)")

    fig.suptitle(f"Future (2070–2100) ensemble-mean yield loss | {scen}", fontsize=14)

    os.makedirs(out_dir, exist_ok=True)
    fp = os.path.join(out_dir, f"future_2070_2100_ensmean_{scen}_4crops_2haz.pdf")
    fig.savefig(fp, dpi=dpi, bbox_inches="tight")
    plt.show()
    print("Saved:", fp)


# =========================
# 运行：两张图 ssp126 / ssp585
# =========================
mirca_ds = DB["dist"]  # 直接用你DB里存的
crops = ("maiz", "soyb", "rice", "whea")
hazs  = ("dry", "wet")
year_slice = slice(2070, 2100)

out_dir = r"D:\AAUDE\paper\paper4\revision_v2\fig\manu_fig"

for scen in ("ssp126", "ssp585"):
    maps = {}
    for crop in crops:
        for haz in hazs:
            da_map = spatial_loss_map_onecrop_future_ensmean(
                DB, mirca_ds, crop=crop, haz=haz,
                scen=scen, gcms=gcms,
                year_slice=year_slice
            )
            maps[(crop, haz)] = da_map

    draw_panel_4x2_maps(
        maps_dict=maps,
        crops=crops, hazs=hazs,
        scen=scen,
        out_dir=out_dir,
        vmin=0, vmax=8,           # 你历史用的是0-8，这里保持一致便于对比
        extent=(-130, 160, -60, 75),
        cmap_name="RdYlBu_r",
        dpi=300
    )

# %% cell 19
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42

scen_labels = {"hist":"Historical", "ssp126":"SSP1-2.6", "ssp585":"SSP5-8.5"}
scen_colors = {"ssp126": "#5AB4AF", "ssp585": "#F8DB71"}
hist_color  = "#696f70"

def plot_from_cache(ds, out_pdf=None):
    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 6),
                             sharex=True, constrained_layout=True)
    axes = axes.ravel()
    legend_handles, legend_labels = [], []

    for i, crop in enumerate(crops):
        ax = axes[i]

        # hist：dry 实线 / wet 虚线
        for haz, ls, lab in [("dry", "-", "Drought"), ("wet", "--", "Waterlogging")]:
            ts = ds["loss"].sel(crop=crop, scenario="hist", hazard=haz, stat="mean")
            h = ax.plot(ts["year"].values, ts.values, color=hist_color, ls=ls, lw=1.6)[0]
            if i == 0:
                legend_handles.append(h)
                legend_labels.append(f"Historical - {lab}")

        # future：ssp126/ssp585，均值线+分位带（起始值对齐到hist末尾）
        for scen in ["ssp126", "ssp585"]:
            col = scen_colors[scen]
            for haz, ls, lab in [("dry", "-", "Drought"), ("wet", "--", "Waterlogging")]:
                # --- 取未来三条 ---
                mean = ds["loss"].sel(crop=crop, scenario=scen, hazard=haz, stat="mean")
                low  = ds["loss"].sel(crop=crop, scenario=scen, hazard=haz, stat="low")
                high = ds["loss"].sel(crop=crop, scenario=scen, hazard=haz, stat="high")

                # 你原来的特例
                if crop == "whea" and haz == "wet":
                    mean = mean * 2
                    low  = low  * 2
                    high = high * 2

                if crop == "whea" and haz == "dry":
                    hist_mean = ds["loss"].sel(crop=crop, scenario="hist", hazard=haz, stat="mean")
                    y_hist_end = float(hist_mean.sel(year=2019).values)
                    y_fut_start = float(mean.sel(year=2020).values)
                    offset = y_hist_end - y_fut_start
    
                    mean = mean + offset
                    low  = low  + offset
                    high = high + offset

                # --- 画图 ---
                line = ax.plot(mean["year"].values, mean.values, color=col, ls=ls, lw=1.6)[0]
                ax.fill_between(mean["year"].values, low.values, high.values,
                                color=col, alpha=0.18, linewidth=0)

                if i == 0:
                    legend_handles.append(line)
                    legend_labels.append(f"{scen_labels[scen]} - {lab}")

        # panel 格式
        ax.set_title(crop_labels.get(crop, crop), loc="left", fontsize=12, pad=6)
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_ylabel("Yield loss (%)")
        if i < 3:
            ax.set_xlabel("")

        ax.set_ylim(0,10)
    fig.legend(legend_handles, legend_labels,
               loc="upper center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, 1.02))

    if out_pdf is not None:
        fig.savefig(out_pdf, dpi=300, bbox_inches="tight")

    plt.show()
    return fig, axes


# ======= 读取缓存并画图 =======
cache_fp = r"D:\AAUDE\paper\paper4\revision_v2\fig\cache_loss_timeseries_4crops.nc"
ds = xr.open_dataset(cache_fp)

fig, axes = plot_from_cache(
    ds,
    out_pdf=r"D:\AAUDE\paper\paper4\revision_v2\fig\future\fig_loss_timeseries_4crops.pdf"
)

# %% cell 20
def spatial_loss_map_onecrop(
    DB_block, mirca_ds,
    crop, haz,
    years, denom="ws", eps=1e-12
):
    """
    DB_block: 形如 DB["hist"] 或 DB["future"][gcm][scen]（单一场景块）
    返回：DataArray(lat, lon) —— years 多年平均的损失百分比（%）

    计算逻辑（逐年）：
      loss_t = (Yrel_nws - Yrel_ws) * (Yp * Area)  [t]
      denom_t = (Yrel_ws * Yp * Area) 或 (Yrel_nws * Yp * Area)  [t]
      loss_pct_t = loss_t / denom_t * 100
      最后对 years 在 year 维度求 mean
    """

    # --- 面积：mirca_ds 里你已经是 "面积(m2)"（因为乘了 cellarea），这里转 ha
    A_rf = (mirca_ds[f"{crop}_nonIrr"] / 10000.0)  # ha
    A_ir = (mirca_ds[f"{crop}_Irr"]    / 10000.0)  # ha

    # --- 潜在产量（t/ha），无时间维
    Yp_rf = DB_block["potential_yield"][crop]["nonIrr"]
    Yp_ir = DB_block["potential_yield"][crop]["Irr"]

    pot_rf = (Yp_rf * A_rf).fillna(0)   # t
    pot_ir = (Yp_ir * A_ir).fillna(0)   # t

    # --- 相对产量：有 year 维（你 open_da 已经 rename 到 year）
    ws_rf  = DB_block["relative_yield"][crop]["waterstress_nonIrr"].sel(year=years)
    nws_rf = DB_block["relative_yield"][crop]["nowaterstress_nonIrr"].sel(year=years)
    ws_ir  = DB_block["relative_yield"][crop]["waterstress_Irr"].sel(year=years)
    nws_ir = DB_block["relative_yield"][crop]["nowaterstress_Irr"].sel(year=years)

    # --- 事件掩膜（year,lat,lon）
    mask_rf = (DB_block["stress_days"][crop][f"{haz}_nonIrr"].sel(year=years) > 0)
    mask_ir = (DB_block["stress_days"][crop][f"{haz}_Irr"].sel(year=years) > 0)

    # --- 损失量（t）：逐年
    loss_rf = ((nws_rf - ws_rf).where(mask_rf) * pot_rf)
    loss_ir = ((nws_ir - ws_ir).where(mask_ir) * pot_ir)
    loss_t  = loss_rf.fillna(0) + loss_ir.fillna(0)

    # --- 分母（t）：逐年
    if denom == "ws":
        den_t = (ws_rf * pot_rf + ws_ir * pot_ir)
    elif denom == "nws":
        den_t = (nws_rf * pot_rf + nws_ir * pot_ir)
    else:
        raise ValueError("denom must be 'ws' or 'nws'")

    # den_t = den_t.fillna(0)
    a = (nws_rf * pot_rf).fillna(0)
    b = (nws_ir * pot_ir).fillna(0)
    
    a2, b2 = xr.align(a, b, join="outer", fill_value=0)
    den_t = a2 + b2
    # --- 年损失百分比 -> 多年平均
    loss_pct_year = (loss_t / (den_t + eps)) * 100.0
    loss_pct_year = loss_pct_year.where(den_t > 0)

    # 多年平均（你也可以改成 median / sum 等）
    loss_pct = loss_pct_year.mean("year", skipna=True)

    # 保留空间有效区
    loss_pct = loss_pct.where(np.isfinite(loss_pct))

    return loss_pct.astype("float32")

def compute_hist_maps(DB, mirca_ds, years, crops, hazards, denom="ws"):
    out = []
    for crop in crops:
        for haz in hazards:
            da = spatial_loss_map_onecrop(
                DB["hist"], mirca_ds,
                crop=crop, haz=haz,
                years=years, denom=denom
            )
            da = da.expand_dims(crop=[crop], haz=[haz])
            out.append(da)

    loss = xr.combine_by_coords(out)  # dims: crop,haz,lat,lon
    ds = xr.Dataset({"loss_pct": loss})
    ds["loss_pct"].attrs.update({
        "long_name": "Mean yield loss percentage",
        "units": "%",
        "period": f"{years.start}-{years.stop}",
        "denom": denom,
        "note": "Computed as mean over years of (loss_t/den_t)*100, loss masked by stress_days>0"
    })
    return ds
ds_hist = compute_hist_maps(
    DB=DB, mirca_ds=mirca_ds,
    years=year_sel_hist,
    crops=crops, hazards=hazards,
    denom="ws"
)

# %% cell 21
ds_hist['loss_pct'][0][0].plot()

# %% cell 22
import os
import numpy as np
import xarray as xr

# =========================
# 0) 配置区
# =========================
crops = ["maiz", "soyb", "rice", "whea"]
hazards = ["dry", "wet"]

# 你的未来 GCM 列表（要与你 DB["future"].keys() 一致）
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

year_sel_hist   = slice(1980, 2019)
year_sel_future = slice(2060, 2100)

OUT_DIR = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps"
os.makedirs(OUT_DIR, exist_ok=True)


# =========================
# 1) 核心：单作物空间损失（多年平均）
# =========================
def spatial_loss_map_onecrop(
    DB_block, mirca_ds,
    crop, haz,
    years, denom="ws", eps=1e-12
):
    """
    DB_block: 形如 DB["hist"] 或 DB["future"][gcm][scen]（单一场景块）
    返回：DataArray(lat, lon) —— years 多年平均的损失百分比（%）

    计算逻辑（逐年）：
      loss_t = (Yrel_nws - Yrel_ws) * (Yp * Area)  [t]
      denom_t = (Yrel_ws * Yp * Area) 或 (Yrel_nws * Yp * Area)  [t]
      loss_pct_t = loss_t / denom_t * 100
      最后对 years 在 year 维度求 mean
    """

    # --- 面积：mirca_ds 里你已经是 "面积(m2)"（因为乘了 cellarea），这里转 ha
    A_rf = (mirca_ds[f"{crop}_nonIrr"] / 10000.0)  # ha
    A_ir = (mirca_ds[f"{crop}_Irr"]    / 10000.0)  # ha

    # --- 潜在产量（t/ha），无时间维
    Yp_rf = DB_block["potential_yield"][crop]["nonIrr"]
    Yp_ir = DB_block["potential_yield"][crop]["Irr"]

    pot_rf = (Yp_rf * A_rf).fillna(0)   # t
    pot_ir = (Yp_ir * A_ir).fillna(0)   # t

    # --- 相对产量：有 year 维（你 open_da 已经 rename 到 year）
    ws_rf  = DB_block["relative_yield"][crop]["waterstress_nonIrr"].sel(year=years)
    nws_rf = DB_block["relative_yield"][crop]["nowaterstress_nonIrr"].sel(year=years)
    ws_ir  = DB_block["relative_yield"][crop]["waterstress_Irr"].sel(year=years)
    nws_ir = DB_block["relative_yield"][crop]["nowaterstress_Irr"].sel(year=years)

    # --- 事件掩膜（year,lat,lon）
    mask_rf = (DB_block["stress_days"][crop][f"{haz}_nonIrr"].sel(year=years) > 0)
    mask_ir = (DB_block["stress_days"][crop][f"{haz}_Irr"].sel(year=years) > 0)

    # --- 损失量（t）：逐年
    loss_rf = ((nws_rf - ws_rf).where(mask_rf) * pot_rf)
    loss_ir = ((nws_ir - ws_ir).where(mask_ir) * pot_ir)
    loss_t  = loss_rf.where(loss_rf>0).fillna(0) + loss_ir.where(loss_ir>0).fillna(0)

    # --- 分母（t）：逐年
    a = (nws_rf * pot_rf).fillna(0)
    b = (nws_ir * pot_ir).fillna(0)
    
    a2, b2 = xr.align(a, b, join="outer", fill_value=0)
    den_t = a2 + b2

    # --- 年损失百分比 -> 多年平均
    loss_pct_year = (loss_t / (den_t + eps)) * 100.0
    loss_pct_year = loss_pct_year.where(den_t > 0)

    # 多年平均（你也可以改成 median / sum 等）
    loss_pct = loss_pct_year.mean("year", skipna=True)

    # 保留空间有效区
    loss_pct = loss_pct.where(np.isfinite(loss_pct))

    return loss_pct.astype("float32")


# =========================
# 2) 计算历史（1980-2019）
# =========================
def compute_hist_maps(DB, mirca_ds, years, crops, hazards, denom="ws"):
    out = []
    for crop in crops:
        for haz in hazards:
            da = spatial_loss_map_onecrop(
                DB["hist"], mirca_ds,
                crop=crop, haz=haz,
                years=years, denom=denom
            )
            da = da.expand_dims(crop=[crop], haz=[haz])
            out.append(da)

    loss = xr.combine_by_coords(out)  # dims: crop,haz,lat,lon
    ds = xr.Dataset({"loss_pct": loss})
    ds["loss_pct"].attrs.update({
        "long_name": "Mean yield loss percentage",
        "units": "%",
        "period": f"{years.start}-{years.stop}",
        "denom": denom,
        "note": "Computed as mean over years of (loss_t/den_t)*100, loss masked by stress_days>0"
    })
    return ds


# =========================
# 3) 计算未来（2070-2100）按 GCM
# =========================
def compute_future_maps_by_gcm(DB, mirca_ds, years, crops, hazards, gcms, scenarios=("ssp126","ssp585"), denom="ws"):
    out_all = []
    for gcm in gcms:
        if gcm not in DB["future"]:
            print(f"[WARN] DB['future'] missing gcm={gcm}, skip.")
            continue

        for scen in scenarios:
            if scen not in DB["future"][gcm]:
                print(f"[WARN] DB['future'][{gcm}] missing scen={scen}, skip.")
                continue

            block = DB["future"][gcm][scen]

            for crop in crops:
                for haz in hazards:
                    da = spatial_loss_map_onecrop(
                        block, mirca_ds,
                        crop=crop, haz=haz,
                        years=years, denom=denom
                    )
                    da = da.expand_dims(gcm=[gcm], scen=[scen], crop=[crop], haz=[haz])
                    out_all.append(da)

    loss = xr.combine_by_coords(out_all)  # dims: gcm,scen,crop,haz,lat,lon
    ds = xr.Dataset({"loss_pct": loss})
    ds["loss_pct"].attrs.update({
        "long_name": "Mean yield loss percentage (by GCM)",
        "units": "%",
        "period": f"{years.start}-{years.stop}",
        "denom": denom,
        "note": "Computed as mean over years of (loss_t/den_t)*100, loss masked by stress_days>0"
    })
    return ds


# =========================
# 4) 运行：你只需要保证 DB 与 mirca_ds 已经在内存中
# =========================
# mirca_ds = DB["dist"]  # 你 build_db_all_gcms 里就是这么存的
mirca_ds = DB["dist"]

# --- 历史
ds_hist = compute_hist_maps(
    DB=DB, mirca_ds=mirca_ds,
    years=year_sel_hist,
    crops=crops, hazards=hazards,
    denom="ws"
)
fp_hist = os.path.join(OUT_DIR, "loss_hist_1980_2019.nc")
ds_hist.to_netcdf(fp_hist)
print("Saved:", fp_hist)

# --- 未来：按 GCM
ds_fut = compute_future_maps_by_gcm(
    DB=DB, mirca_ds=mirca_ds,
    years=year_sel_future,
    crops=crops, hazards=hazards,
    gcms=gcms,
    scenarios=("ssp126","ssp585"),
    denom="ws"
)
fp_fut = os.path.join(OUT_DIR, "loss_future_2060_2100_by_gcm.nc")
ds_fut.to_netcdf(fp_fut)
print("Saved:", fp_fut)

# （可选）如果你还想要未来 ensemble 均值（跨 gcm）也顺手输出：
ds_fut_ens = ds_fut.mean("gcm", skipna=True)
fp_fut_ens = os.path.join(OUT_DIR, "loss_future_2060_2100_ensemble_mean.nc")
ds_fut_ens.to_netcdf(fp_fut_ens)
print("Saved:", fp_fut_ens)

# %% cell 23
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ======== 字体（你已有的话可省略）========
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 0) 配置
# =========================
scenarios = ["hist", "ssp126", "ssp585"]
scen_labels = {"hist":"Historical", "ssp126":"SSP1-2.6", "ssp585":"SSP5-8.5"}

crops = ["maiz", "soyb", "rice", "whea"]
crop_labels = {"maiz":"Maize","soyb":"Soybean","rice":"Rice","whea":"Wheat"}

hazards = ["dry", "wet"]
haz_labels  = {"dry":"Drought", "wet":"Waterlogging"}

EXTENT = [-130, 160, -60, 75]  # 你指定的范围
ASPECT = 1.3

cmap = plt.get_cmap("RdYlBu_r")
vmin, vmax = 0.0, 8.0                       # 你给的 0-8
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

# 输入 nc
hist_nc = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\loss_hist_1980_2019.nc"
fut_nc  = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\loss_future_2060_2100_ensemble_mean.nc"
# 如果你没有 ensemble_mean nc，而是 by_gcm：
# fut_by_gcm_nc = r"...\loss_future_2070_2100_by_gcm.nc"


# =========================
# 1) 读取并拼成统一数据：loss(scne,crop,haz,lat,lon)
# =========================
ds_hist = xr.open_dataset(hist_nc)
hist = ds_hist["loss_pct"]  # (crop,haz,lat,lon)

# --- 未来（优先用 ensemble_mean）
ds_fut = xr.open_dataset(fut_nc)
fut = ds_fut["loss_pct"]    # (scen,crop,haz,lat,lon)  scen: ssp126/ssp585

# # 若你只有 by_gcm，可用下面两种方式：
# ds_fut_by = xr.open_dataset(fut_by_gcm_nc)
# fut = ds_fut_by["loss_pct"].mean("gcm", skipna=True)  # ensemble mean
# # 或选定某个 gcm
# # fut = ds_fut_by["loss_pct"].sel(gcm="CanESM5")

# 给历史补 scen 维，和未来对齐
hist = hist.expand_dims(scen=["hist"])

# 拼接 scen 维
loss = xr.concat([hist, fut], dim="scen")  # (scen,crop,haz,lat,lon)
loss = loss.sel(scen=scenarios, crop=crops, haz=hazards)

# 确保 lat 从南到北（不一致会导致画图倒）
if loss.lat[0] > loss.lat[-1]:
    loss = loss.sortby("lat")

# =========================
# 2) 绘图：8行×3列，统一 colorbar
# =========================
row_keys = [(crop, haz) for crop in crops for haz in hazards]  # 8行

nrows, ncols = len(row_keys), len(scenarios)
fig = plt.figure(figsize=(12.5, 18.0))

# 预留右侧色标位置
gs = fig.add_gridspec(nrows, ncols, left=0.02, right=0.88, top=0.98, bottom=0.02,
                      wspace=0.02, hspace=0.04)

mappable = None

for r, (crop, haz) in enumerate(row_keys):
    for c, scen in enumerate(scenarios):
        ax = fig.add_subplot(gs[r, c], projection=ccrs.PlateCarree())

        ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
        ax.set_aspect(ASPECT)

        da = loss.sel(scen=scen, crop=crop, haz=haz)

        # 注意：这里不单独 add_colorbar，最后统一画
        im = da.plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            x="lon", y="lat",
            cmap=cmap,
            norm=norm,
            add_colorbar=False
        )
        if mappable is None:
            mappable = im

        # 海岸线 & 陆地底色
        ax.coastlines(resolution="110m", linewidth=0.8)
        ax.add_feature(cfeature.LAND, facecolor="whitesmoke", zorder=0)

        # 去边框
        ax.set_frame_on(False)
        ax.patch.set_visible(False)
        if hasattr(ax, "spines"):
            for spine in ax.spines.values():
                spine.set_visible(False)

        # 标题：第一行写情景
        if r == 0:
            ax.set_title(scen_labels[scen], fontsize=12)
        else:
            ax.set_title("")

        # 行标签：第一列写 作物+胁迫
        if c == 0:
            ax.text(
                -0.02, 0.5,
                f"{crop_labels[crop]}\n{haz_labels[haz]}",
                transform=ax.transAxes,
                ha="right", va="center",
                fontsize=11
            )

        # 去掉经纬度刻度文字（你喜欢也可保留外圈）
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)

# 统一 colorbar（右侧）
cax = fig.add_axes([0.90, 0.12, 0.018, 0.76])  # [left,bottom,width,height]
cb = fig.colorbar(
    mappable, cax=cax, orientation="vertical",
    extend="neither"
)
cb.set_label("Yield loss (%)")

# 保存 / 显示
out_pdf = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\Fig_loss_8x3_crops_hazards.pdf"
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.show()

print("Saved:", out_pdf)

# %% cell 24
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ======== 字体（你已有的话可省略）========
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 0) 配置
# =========================
scenarios   = ["hist", "ssp126", "ssp585"]
scen_labels = {"hist":"Historical", "ssp126":"SSP1-2.6", "ssp585":"SSP5-8.5"}

crops       = ["maiz", "soyb", "rice", "whea"]
crop_labels = {"maiz":"Maize","soyb":"Soybean","rice":"Rice","whea":"Wheat"}

hazards     = ["dry", "wet"]
haz_labels  = {"dry":"Drought", "wet":"Waterlogging"}

EXTENT = [-130, 160, -60, 75]
ASPECT = 1.3

cmap = plt.get_cmap("RdYlBu_r")
vmin, vmax = 0.0, 8.0
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

# 输入 nc
hist_nc = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\loss_hist_1980_2019.nc"
fut_nc  = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\loss_future_2060_2100_ensemble_mean.nc"

# 输出
out_nc  = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\loss_cropmean_3scen_2haz.nc"
out_pdf = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\Fig_loss_cropmean_2x3.pdf"

# =========================
# 1) 读入并拼接为统一：loss(scen,crop,haz,lat,lon)
# =========================
ds_hist = xr.open_dataset(hist_nc)
hist = ds_hist["loss_pct"]  # (crop,haz,lat,lon)

ds_fut = xr.open_dataset(fut_nc)
fut = ds_fut["loss_pct"]    # (scen,crop,haz,lat,lon)  scen: ssp126/ssp585

# 给历史补 scen 维
hist = hist.expand_dims(scen=["hist"])

# 拼接 scen 维
loss = xr.concat([hist, fut], dim="scen")  # (scen,crop,haz,lat,lon)
loss = loss.sel(scen=scenarios, crop=crops, haz=hazards)

# lat 保证从南到北
if loss.lat[0] > loss.lat[-1]:
    loss = loss.sortby("lat")

# =========================
# 2) 只对四种作物求均值（不对 dry/wet 求均值）
#    得到：loss_crop_mean(scen,haz,lat,lon)
# =========================
loss_crop_mean = loss.mean(dim="crop", skipna=True)

# 可选：保存 nc，便于后续复用
loss_crop_mean.to_dataset(name="loss_pct_cropmean").to_netcdf(out_nc)
print("Saved nc:", out_nc)

# =========================
# 3) 绘图：2行(胁迫) × 3列(情景)，统一色标
# =========================
nrows, ncols = len(hazards), len(scenarios)
fig = plt.figure(figsize=(12.5, 7.5))

# 预留右侧色标位置
gs = fig.add_gridspec(
    nrows, ncols,
    left=0.03, right=0.88, top=0.93, bottom=0.06,
    wspace=0.02, hspace=0.08
)

mappable = None

for r, haz in enumerate(hazards):
    for c, scen in enumerate(scenarios):
        ax = fig.add_subplot(gs[r, c], projection=ccrs.PlateCarree())
        ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
        ax.set_aspect(ASPECT)

        da = loss_crop_mean.sel(scen=scen, haz=haz)

        im = da.plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            x="lon", y="lat",
            cmap=cmap,
            norm=norm,
            add_colorbar=False
        )
        if mappable is None:
            mappable = im

        ax.coastlines(resolution="110m", linewidth=0.8)
        ax.add_feature(cfeature.LAND, facecolor="whitesmoke", zorder=0)

        # 去边框/刻度
        ax.set_frame_on(False)
        ax.patch.set_visible(False)
        if hasattr(ax, "spines"):
            for spine in ax.spines.values():
                spine.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)

        # 列标题：情景
        if r == 0:
            ax.set_title(scen_labels[scen], fontsize=12)
        else:
            ax.set_title("")

        # 行标题：胁迫（放在第一列左侧）
        if c == 0:
            ax.text(
                -0.02, 0.5,
                haz_labels[haz] + "\n(Multi-crop mean)",
                transform=ax.transAxes,
                ha="right", va="center",
                fontsize=11
            )

# 统一 colorbar（右侧）
cax = fig.add_axes([0.90, 0.16, 0.018, 0.70])
cb = fig.colorbar(mappable, cax=cax, orientation="vertical", extend="neither")
cb.set_label("Yield loss (%)")

# 保存 / 显示
# fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.show()

print("Saved pdf:", out_pdf)

# %% cell 25
import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt

# ===== 字体/矢量PDF =====
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 0) 配置
# =========================
nc_fp   = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\loss_cropmean_3scen_2haz.nc"
varname = "loss_pct_cropmean"   # 你保存nc时的变量名

hazards     = ["dry", "wet"]
haz_labels  = {"dry":"Drought", "wet":"Waterlogging"}

scenarios_f = ["ssp126", "ssp585"]
scen_labels = {"ssp126":"SSP1-2.6", "ssp585":"SSP5-8.5"}

# 两情景颜色（你也可以换成自己论文配色）
scen_colors = {"ssp126":"#5AB4AF", "ssp585":"#F8DB71"}

out_pdf = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\Fig_lat_mirror_delta_future_minus_hist.pdf"

# =========================
# 1) 读取 & 计算 Δ(未来-历史) 的纬向平均
# =========================
ds = xr.open_dataset(nc_fp)
da = ds[varname]  # (scen, haz, lat, lon)

# lat从南到北
if da.lat[0] > da.lat[-1]:
    da = da.sortby("lat")

# Δ(scen,haz,lat,lon)
delta = da.sel(scen=scenarios_f) - da.sel(scen="hist")

# 纬向平均：对 lon 求平均 -> Δ(scen,haz,lat)
delta_lat = delta.mean(dim="lon", skipna=True)

lat = delta_lat["lat"].values

# =========================
# 2) （可选）cropland fraction by latitude
# -------------------------
# 你如果有 cropland/harvest area（如 MIRCA/SPAM），应计算：
#   frac_lat = 100 * area.sum(lon) / area.sum(lat,lon)
#
# 这里先提供一个“可插拔”的入口：
# - 没有就设为 None，不画灰色带
# =========================
frac_lat = None   # <- 如果你之后算好了，替换为与 lat 同长度的一维数组（单位：%）

# =========================
# 3) 画图：每个 hazard 一个面板（2行×1列）
#    - 以 x=0 为中线：正值填右侧（increase），负值填左侧（decline）
# =========================
fig, axes = plt.subplots(2, 1, figsize=(1.5, 6.2), sharex=True, constrained_layout=True)

# 统一x轴范围（按两情景、两haz的最大绝对值）
absmax = float(np.nanmax(np.abs(delta_lat.values)))
xlim = (-absmax * 1.05, absmax * 1.05)

for i, haz in enumerate(hazards):
    ax = axes[i]

    # 0 中线（虚线）
    ax.axvline(0, lw=1.0, ls="--", color="0.35", zorder=5)

    # 背景水平参考线（每30°）
    # for y in [-60, -30, 0, 30, 60, 90]:
    #     ax.axhline(y, lw=0.6, color="0.85", zorder=0)

    # 两个未来情景分别画镜像填色
    for scen in scenarios_f:
        y = delta_lat.sel(scen=scen, haz=haz).values  # Δloss(%) 随纬度变化（已lon平均）
        col = scen_colors[scen]

        pos = np.clip(y, 0, None)     # 正：右侧
        neg = np.clip(y, None, 0)     # 负：左侧（neg本身为负）

        # 填色（镜像风格）
        ax.fill_betweenx(lat, 0, pos, color=col, alpha=0.45, lw=0.0)
        ax.fill_betweenx(lat, 0, neg, color=col, alpha=0.45, lw=0.0)

        # 轮廓线（更像你示例图的“边界”）
        ax.plot(pos, lat, color=col, lw=1.1)
        ax.plot(neg, lat, color=col, lw=1.1)

    # 面板标题
    ax.set_title(haz_labels[haz], loc="left", fontsize=12, pad=4)

    # y轴（纬度）
    ax.set_ylabel("Latitude (°)")
    ax.set_ylim(lat.min(), lat.max())

    # 去掉上/右边框
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(-60,80)
    ax.grid(False)
    # 顶部 cropland extent（可选）
    if frac_lat is not None:
        axt = ax.twiny()
        # 灰色“cropland extent”对称填充（类似你图中灰色轮廓）
        # 左右对称：[-frac, +frac]
        axt.fill_betweenx(lat, -frac_lat, frac_lat, color="0.80", alpha=1.0, lw=0.0, zorder=1)

        # 顶部坐标轴设为对称并用绝对值刻度显示（10 5 0 5 10）
        fmax = float(np.nanmax(frac_lat)) * 1.1
        axt.set_xlim(-fmax, fmax)
        ticks = np.linspace(-fmax, fmax, 5)
        axt.set_xticks(ticks)
        axt.set_xticklabels([f"{abs(t):.0f}" for t in ticks])
        axt.set_xlabel("Fraction of global cropland (%)", labelpad=6)

        # 顶部轴样式
        axt.spines["right"].set_visible(False)
        axt.spines["left"].set_visible(False)

# x轴标签
axes[-1].set_xlabel("Δ Yield loss (future − historical) (%)")
axes[-1].set_xlim(*xlim)

# 图例（用色块代表情景）
handles = []
labels = []
for scen in scenarios_f:
    h = plt.Line2D([0], [0], color=scen_colors[scen], lw=6, alpha=0.6)
    handles.append(h)
    labels.append(f"{scen_labels[scen]}")

fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.52, 1.01))

fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.show()
print("Saved:", out_pdf)

# %% cell 26
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ======== 字体（你已有的话可省略）========
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 0) 配置
# =========================
scenarios   = ["hist", "ssp126", "ssp585"]
scen_labels = {"hist":"Historical", "ssp126":"SSP1-2.6", "ssp585":"SSP5-8.5"}

crops       = ["maiz", "soyb", "rice", "whea"]
crop_labels = {"maiz":"Maize","soyb":"Soybean","rice":"Rice","whea":"Wheat"}

hazards     = ["dry", "wet"]
haz_labels  = {"dry":"Drought", "wet":"Waterlogging"}

EXTENT = [-130, 160, -60, 75]
ASPECT = 1.3

cmap = plt.get_cmap("RdYlBu_r")
vmin, vmax = 0.0, 8.0
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

# 输入 nc
hist_nc = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\loss_hist_1980_2019.nc"
fut_nc  = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\loss_future_2070_2100_ensemble_mean.nc"


# 输出
out_nc  = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\loss_cropmean_area_weighted_3scen_2haz.nc"
out_pdf = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\Fig_loss_cropmean_area_weighted_2x3.pdf"

# =========================
# 1) 读入 loss 并拼接为统一：loss(scen,crop,haz,lat,lon)
# =========================
ds_hist = xr.open_dataset(hist_nc)
hist = ds_hist["loss_pct"]  # (crop,haz,lat,lon)

ds_fut = xr.open_dataset(fut_nc)
fut = ds_fut["loss_pct"]    # (scen,crop,haz,lat,lon)  scen: ssp126/ssp585

hist = hist.expand_dims(scen=["hist"])
loss = xr.concat([hist, fut], dim="scen")
loss = loss.sel(scen=scenarios, crop=crops, haz=hazards)

# lat 保证从南到北
if loss.lat[0] > loss.lat[-1]:
    loss = loss.sortby("lat")

# =========================
# 2) 读入面积 & 生成权重 w(crop,lat,lon)
# =========================
dist = DB['dist']  # 里面有 maiz_Irr, maiz_nonIrr ... whea_Irr/whea_nonIrr等

# 你面积里小麦是 wwh/swh + whea 可能全nan，这里按你四作物：maiz/soyb/rice/whea
# maiz/soyb/rice：直接用 *_Irr + *_nonIrr
# whea：优先用 whea_Irr/nonIrr；若全nan则回退用 (wwh + swh) 的 Irr/nonIrr 合并

def get_crop_area(dist_ds, crop):
    if crop in ["maiz", "soyb", "rice"]:
        a = dist_ds[f"{crop}_Irr"] + dist_ds[f"{crop}_nonIrr"]
        return a
    elif crop == "whea":
        # 先尝试 whea
        if ("whea_Irr" in dist_ds) and ("whea_nonIrr" in dist_ds):
            a = dist_ds["whea_Irr"] + dist_ds["whea_nonIrr"]
            # 如果几乎全是nan，则用 wwh+swh
            if int(np.isfinite(a.values).sum()) == 0:
                a = (dist_ds["wwh_Irr"] + dist_ds["wwh_nonIrr"] +
                     dist_ds["swh_Irr"] + dist_ds["swh_nonIrr"])
            return a
        else:
            # 没有 whea 就用 wwh+swh
            a = (dist_ds["wwh_Irr"] + dist_ds["wwh_nonIrr"] +
                 dist_ds["swh_Irr"] + dist_ds["swh_nonIrr"])
            return a
    else:
        raise ValueError(crop)

# 组装成 (crop,lat,lon)
area_list = []
for c in crops:
    a = get_crop_area(dist, c)
    a = a.rename("area").expand_dims(crop=[c])
    area_list.append(a)

area = xr.concat(area_list, dim="crop")  # (crop,lat,lon)

# 对齐到 loss 的网格（如果完全一致，这步不会改变）
area = area.interp(lat=loss.lat, lon=loss.lon, method="nearest")

# 把 NaN 当 0 面积（避免传播）
area = area.fillna(0.0)

# =========================
# 3) 面积加权多作物均值：对 crop 做加权平均（保留 haz & scen）
# =========================
# 分子：sum_c (loss * area)
num = (loss * area).sum(dim="crop", skipna=True)

# 分母：sum_c area
den = area.sum(dim="crop")

# 避免除0：没有作物面积的网格设为 NaN
loss_crop_wmean = xr.where(den > 0, num / den, np.nan)

# 保存 nc
loss_crop_wmean.to_dataset(name="loss_pct_cropmean_areaW").to_netcdf(out_nc)
print("Saved nc:", out_nc)

# =========================
# 4) 绘图：2行(胁迫) × 3列(情景)，统一色标
# =========================
nrows, ncols = len(hazards), len(scenarios)
fig = plt.figure(figsize=(12.5, 7.5))
gs = fig.add_gridspec(
    nrows, ncols,
    left=0.03, right=0.88, top=0.93, bottom=0.06,
    wspace=0.02, hspace=0.08
)

mappable = None

for r, haz in enumerate(hazards):
    for c, scen in enumerate(scenarios):
        ax = fig.add_subplot(gs[r, c], projection=ccrs.PlateCarree())
        ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
        ax.set_aspect(ASPECT)

        da = loss_crop_wmean.sel(scen=scen, haz=haz)

        im = da.plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            # x="lon", y="lat",
            cmap=cmap,
            norm=norm,
            add_colorbar=False
        )
        if mappable is None:
            mappable = im

        ax.coastlines(resolution="110m", linewidth=0.8)
        ax.add_feature(cfeature.LAND, facecolor="whitesmoke", zorder=0)

        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        ax.set_frame_on(False)
        ax.patch.set_visible(False)
        if hasattr(ax, "spines"):
            for sp in ax.spines.values():
                sp.set_visible(False)

        if r == 0:
            ax.set_title(scen_labels[scen], fontsize=12)
        else:
            ax.set_title("")

        if c == 0:
            ax.text(
                -0.02, 0.5,
                haz_labels[haz] + "\n(Area-weighted multi-crop mean)",
                transform=ax.transAxes,
                ha="right", va="center",
                fontsize=11
            )

# 统一色标
cax = fig.add_axes([0.90, 0.16, 0.018, 0.70])
cb = fig.colorbar(mappable, cax=cax, orientation="vertical", extend="neither")
cb.set_label("Yield loss (%)")

fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.show()
print("Saved pdf:", out_pdf)

# %% cell 27
import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt

# ===== 字体/矢量PDF =====
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 0) 配置
# =========================
nc_fp   = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\loss_cropmean_area_weighted_3scen_2haz.nc"
varname = "loss_pct_cropmean_areaW"   # 你保存nc时的变量名

hazards     = ["dry", "wet"]
haz_labels  = {"dry":"Drought", "wet":"Waterlogging"}

scenarios_f = ["ssp126", "ssp585"]
scen_labels = {"ssp126":"SSP1-2.6", "ssp585":"SSP5-8.5"}

# 两情景颜色（你也可以换成自己论文配色）
scen_colors = {"ssp126":"#5AB4AF", "ssp585":"#F8DB71"}

out_pdf = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\Fig_lat_mirror_delta_future_minus_hist.pdf"

# =========================
# 1) 读取 & 计算 Δ(未来-历史) 的纬向平均
# =========================
ds = xr.open_dataset(nc_fp)
da = ds[varname]  # (scen, haz, lat, lon)

# lat从南到北
if da.lat[0] > da.lat[-1]:
    da = da.sortby("lat")

# Δ(scen,haz,lat,lon)
delta = da.sel(scen=scenarios_f) - da.sel(scen="hist")

# 纬向平均：对 lon 求平均 -> Δ(scen,haz,lat)
delta_lat = delta.mean(dim="lon", skipna=True)

lat = delta_lat["lat"].values

# =========================
# 2) （可选）cropland fraction by latitude
# -------------------------
# 你如果有 cropland/harvest area（如 MIRCA/SPAM），应计算：
#   frac_lat = 100 * area.sum(lon) / area.sum(lat,lon)
#
# 这里先提供一个“可插拔”的入口：
# - 没有就设为 None，不画灰色带
# =========================
frac_lat = None   # <- 如果你之后算好了，替换为与 lat 同长度的一维数组（单位：%）

# =========================
# 3) 画图：每个 hazard 一个面板（2行×1列）
#    - 以 x=0 为中线：正值填右侧（increase），负值填左侧（decline）
# =========================
fig, axes = plt.subplots(2, 1, figsize=(2.2, 6.2), sharex=True, constrained_layout=True)

# 统一x轴范围（按两情景、两haz的最大绝对值）
absmax = float(np.nanmax(np.abs(delta_lat.values)))
xlim = (-absmax * 1.05, absmax * 1.05)

for i, haz in enumerate(hazards):
    ax = axes[i]

    # 0 中线（虚线）
    ax.axvline(0, lw=1.0, ls="--", color="0.35", zorder=5)

    # 背景水平参考线（每30°）
    for y in [-60, -30, 0, 30, 60, 90]:
        ax.axhline(y, lw=0.6, color="0.85", zorder=0)

    # 两个未来情景分别画镜像填色
    for scen in scenarios_f:
        y = delta_lat.sel(scen=scen, haz=haz).values  # Δloss(%) 随纬度变化（已lon平均）
        col = scen_colors[scen]

        pos = np.clip(y, 0, None)     # 正：右侧
        neg = np.clip(y, None, 0)     # 负：左侧（neg本身为负）

        # 填色（镜像风格）
        ax.fill_betweenx(lat, 0, pos, color=col, alpha=0.45, lw=0.0)
        ax.fill_betweenx(lat, 0, neg, color=col, alpha=0.45, lw=0.0)

        # 轮廓线（更像你示例图的“边界”）
        ax.plot(pos, lat, color=col, lw=1)
        ax.plot(neg, lat, color=col, lw=1)

    # 面板标题
    ax.set_title(haz_labels[haz], loc="left", fontsize=12, pad=4)

    # y轴（纬度）
    ax.set_ylabel("Latitude (°)")
    ax.set_ylim(lat.min(), lat.max())

    # 去掉上/右边框
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_ylim(-60,80)
    ax.set_xlim(-3,3)
    # 顶部 cropland extent（可选）
    if frac_lat is not None:
        axt = ax.twiny()
        # 灰色“cropland extent”对称填充（类似你图中灰色轮廓）
        # 左右对称：[-frac, +frac]
        axt.fill_betweenx(lat, -frac_lat, frac_lat, color="0.80", alpha=1.0, lw=0.0, zorder=1)

        # 顶部坐标轴设为对称并用绝对值刻度显示（10 5 0 5 10）
        fmax = float(np.nanmax(frac_lat)) * 1.1
        axt.set_xlim(-fmax, fmax)
        ticks = np.linspace(-fmax, fmax, 5)
        axt.set_xticks(ticks)
        axt.set_xticklabels([f"{abs(t):.0f}" for t in ticks])
        axt.set_xlabel("Fraction of global cropland (%)", labelpad=6)

        # 顶部轴样式
        axt.spines["right"].set_visible(False)
        axt.spines["left"].set_visible(False)

# x轴标签
axes[-1].set_xlabel("Δ Yield loss (future − historical) (%)")
axes[-1].set_xlim(*xlim)

# 图例（用色块代表情景）
handles = []
labels = []
for scen in scenarios_f:
    h = plt.Line2D([0], [0], color=scen_colors[scen], lw=6, alpha=0.6)
    handles.append(h)
    labels.append(f"{scen_labels[scen]}")

fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False, bbox_to_anchor=(0.52, 1.01))

fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.show()
print("Saved:", out_pdf)

# %% cell 28
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ===== 字体/矢量PDF =====
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 0) 配置
# =========================
nc_fp   = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\loss_cropmean_area_weighted_3scen_2haz.nc"  # <-- 你的nc
varname = "loss_pct_cropmean_areaW"  # <-- 你保存nc时的变量名

fut_scen = ["ssp126", "ssp585"]
scen_labels = {
    "ssp126": "SSP1-2.6 − Historical",
    "ssp585": "SSP5-8.5 − Historical"
}

hazards = ["dry", "wet"]
haz_labels = {"dry":"Drought", "wet":"Waterlogging"}

EXTENT = [-130, 160, -60, 75]
ASPECT = 1.3

# 发散色带（差值图）
cmap = plt.get_cmap("RdBu_r")

out_pdf = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\Fig_spatial_delta_areaW_future_minus_hist_2x2.pdf"

# =========================
# 1) 读取 & 计算空间差值
# =========================
ds = xr.open_dataset(nc_fp)
da = ds[varname]  # (scen, haz, lat, lon)

# lat从南到北
if da.lat[0] > da.lat[-1]:
    da = da.sortby("lat")

# Δ(scen,haz,lat,lon) = future - hist
delta = da.sel(scen=fut_scen) - da.sel(scen="hist")

# 统一色标：0居中，范围取全局最大绝对值
absmax = float(np.nanmax(np.abs(delta.values)))
absmax = max(absmax, 1e-6)
norm = mcolors.TwoSlopeNorm(vmin=-10, vcenter=0.0, vmax=10)

# =========================
# 2) 绘图：2行(胁迫) × 2列(情景)，统一colorbar
# =========================
fig = plt.figure(figsize=(11.5, 7.2))
gs = fig.add_gridspec(
    nrows=2, ncols=2,
    left=0.03, right=0.90, bottom=0.06, top=0.93,
    wspace=0.03, hspace=0.10
)

mappable = None

for r, haz in enumerate(hazards):
    for c, scen in enumerate(fut_scen):
        ax = fig.add_subplot(gs[r, c], projection=ccrs.PlateCarree())
        ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
        ax.set_aspect(ASPECT)

        da_plot = delta.sel(scen=scen, haz=haz)

        im = da_plot.plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            # x="lon", y="lat",
            cmap=cmap,
            norm=norm,
            add_colorbar=False
        )
        if mappable is None:
            mappable = im

        ax.coastlines(resolution="110m", linewidth=0.8)
        ax.add_feature(cfeature.LAND, facecolor="whitesmoke", zorder=0)

        # 去刻度/边框
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        ax.set_frame_on(False)
        ax.patch.set_visible(False)
        if hasattr(ax, "spines"):
            for sp in ax.spines.values():
                sp.set_visible(False)

        # 列标题：情景差值
        if r == 0:
            ax.set_title(scen_labels[scen], fontsize=12)
        else:
            ax.set_title("")

        # 行标签：胁迫
        if c == 0:
            ax.text(
                -0.02, 0.5,
                haz_labels[haz],
                transform=ax.transAxes,
                ha="right", va="center",
                fontsize=12
            )

# 统一colorbar
# ===== 横向 colorbar（放在底部中间）=====
cax = fig.add_axes([0.25, 0.04, 0.50, 0.025])  # [left, bottom, width, height]
cb = fig.colorbar(
    mappable,
    cax=cax,
    orientation="horizontal",
    extend="both"
)
cb.set_label("Δ Yield loss (%)", labelpad=6)

fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.show()
print("Saved:", out_pdf)

# %% cell 29
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
from matplotlib.patches import Patch
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# ===== 字体/矢量PDF =====
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 0) 配置
# =========================
nc_fp   = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\loss_cropmean_area_weighted_3scen_2haz.nc"
varname = "loss_pct_cropmean_areaW"

fut_scen = ["ssp126", "ssp585"]
scen_labels = {"ssp126": "SSP1-2.6 − Historical", "ssp585": "SSP5-8.5 − Historical"}

hazards = ["dry", "wet"]
haz_labels = {"dry": "Drought", "wet": "Waterlogging"}

EXTENT = [-130, 160, -60, 75]
ASPECT = 1.3

THRESH = 0.5  # Δ <= 0.5 -> little change（只在作物区显示）

# 只对 Δ>0.5 做分级：0.5–1, 1–3, 3–5, 5–10, >10
bounds = [THRESH, 1, 3, 5, 10]
bin_labels = [f"{THRESH}–1", "1–3", "3–5", "5–10", ">10"]

bin_colors = ["#fee5d9", "#fcae91", "#fb6a4a", "#cb181d"]
over_color = "#67000d"

little_change_color = "0.35"  # 深灰（只画在作物区且 Δ<=0.5 的像元）

out_pdf = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\Fig_spatial_delta_bins_legend_thresh05_masked.pdf"

# =========================
# 1) 读取差值
# =========================
ds = xr.open_dataset(nc_fp)
da = ds[varname]  # (scen, haz, lat, lon)
if da.lat[0] > da.lat[-1]:
    da = da.sortby("lat")

delta = da.sel(scen=fut_scen) - da.sel(scen="hist")  # (scen,haz,lat,lon)

# =========================
# 2) 读取面积，构造“作物区掩膜”
#    只在 area_total>0 的地方显示分类
# =========================
dist = DB['dist']

def crop_area(dist_ds, crop):
    if crop in ["maiz", "soyb", "rice"]:
        return (dist_ds[f"{crop}_Irr"] + dist_ds[f"{crop}_nonIrr"]).fillna(0.0)
    elif crop == "whea":
        # 优先 whea，如全nan则用 wwh+swh
        if ("whea_Irr" in dist_ds) and ("whea_nonIrr" in dist_ds):
            a = (dist_ds["whea_Irr"] + dist_ds["whea_nonIrr"])
            if int(np.isfinite(a.values).sum()) == 0:
                a = (dist_ds["wwh_Irr"] + dist_ds["wwh_nonIrr"] +
                     dist_ds["swh_Irr"] + dist_ds["swh_nonIrr"])
        else:
            a = (dist_ds["wwh_Irr"] + dist_ds["wwh_nonIrr"] +
                 dist_ds["swh_Irr"] + dist_ds["swh_nonIrr"])
        return a.fillna(0.0)
    else:
        raise ValueError(crop)

area_total = (
    crop_area(dist, "maiz") +
    crop_area(dist, "soyb") +
    crop_area(dist, "rice") +
    crop_area(dist, "whea")
)

# 对齐到 delta 网格
area_total = area_total.interp(lat=delta.lat, lon=delta.lon, method="nearest")
mask_crop = area_total > 0

# 只在作物区保留 delta，其余设 NaN（背景保持白/透明）
delta = delta.where(mask_crop)

# =========================
# 3) 分成两层来画：
#    A) little-change 层：作物区且 Δ<=0.5 -> 深灰
#    B) increase 分级层：作物区且 Δ>0.5 -> 分级红
# =========================
little = xr.where(delta <= THRESH, 1.0, np.nan)  # 仅用于画灰色
inc = delta.where(delta > THRESH)                # 用于画红色分级

# increase 的离散色带 & norm
cmap_inc = mcolors.ListedColormap(bin_colors, name="inc_bins")
cmap_inc.set_over(over_color)
norm_inc = mcolors.BoundaryNorm(bounds, ncolors=cmap_inc.N, clip=False)

# little-change 的单色 cmap（只一色）
cmap_little = mcolors.ListedColormap([little_change_color], name="little_gray")
norm_little = mcolors.BoundaryNorm([0, 2], ncolors=1)

# =========================
# 4) 绘图：2×2（无 colorbar）
# =========================
fig = plt.figure(figsize=(11.5, 7.2))
gs = fig.add_gridspec(
    nrows=2, ncols=2,
    left=0.03, right=0.98, bottom=0.10, top=0.93,
    wspace=0.03, hspace=0.10
)

for r, haz in enumerate(hazards):
    for c, scen in enumerate(fut_scen):
        ax = fig.add_subplot(gs[r, c], projection=ccrs.PlateCarree())
        ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
        ax.set_aspect(ASPECT)

        # 先画底色（可选）
        ax.add_feature(cfeature.LAND, facecolor="whitesmoke", zorder=0)
        ax.coastlines(resolution="110m", linewidth=0.8, zorder=3)

        # A) 先画灰色 little-change（只在作物区的 Δ<=0.5）
        little.sel(scen=scen, haz=haz).plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            cmap=cmap_little,
            norm=norm_little,
            add_colorbar=False,
            add_labels=False,
            zorder=1
        )

        # B) 再叠加 increase 分级（Δ>0.5）
        inc.sel(scen=scen, haz=haz).plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            cmap=cmap_inc,
            norm=norm_inc,
            add_colorbar=False,
            add_labels=False,
            zorder=2
        )

        # 去刻度/边框
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        ax.set_frame_on(False)
        ax.patch.set_visible(False)
        if hasattr(ax, "spines"):
            for sp in ax.spines.values():
                sp.set_visible(False)

        # 列标题
        if r == 0:
            ax.set_title(scen_labels[scen], fontsize=12)
        else:
            ax.set_title("")

        # 行标签
        if c == 0:
            ax.text(
                -0.02, 0.5,
                haz_labels[haz],
                transform=ax.transAxes,
                ha="right", va="center",
                fontsize=12
            )

# =========================
# 5) 图例（legend）
# =========================
legend_handles = [
    Patch(facecolor=little_change_color, edgecolor="none", label=f"Little change (Δ ≤ {THRESH})"),
]

# 分级（Δ>0.5）
for col, lab in zip(bin_colors, bin_labels[:-1]):  # 0.5–1,1–3,3–5,5–10
    legend_handles.append(Patch(facecolor=col, edgecolor="none", label=f"Increase: {lab}"))
legend_handles.append(Patch(facecolor=over_color, edgecolor="none", label="Increase: >10"))

fig.legend(
    handles=legend_handles,
    loc="lower center",
    ncol=3,
    frameon=False,
    bbox_to_anchor=(0.5, 0.02),
    fontsize=10,
    handlelength=1.6,
    columnspacing=1.5
)

fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.show()
print("Saved:", out_pdf)

# %% cell 30
# -*- coding: utf-8 -*-
import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import seaborn as sns

from matplotlib.patches import Patch

# 需要安装：
# pip install regionmask country_converter openpyxl
import regionmask
import country_converter as coco

# =========================
# 0) 全局绘图参数
# =========================
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 1) 路径配置
# =========================
nc_fp = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\loss_cropmean_area_weighted_3scen_2haz.nc"
varname = "loss_pct_cropmean_areaW"

faostat_fp = r"D:\Edge_download\FAOSTAT_data_en_3-8-2026.csv"
income_fp  = r"D:\Edge_download\CLASS_2025_10_07.xlsx"
shp_fp     = r"D:\AAUDE\paper_v2\paper2\data\ne_10m_admin_0_countries\ne_10m_admin_0_countries.shp"

out_dir = r"D:\AAUDE\paper_v2\paper4\outputs_country_trade_drought"
os.makedirs(out_dir, exist_ok=True)

# =========================
# 2) 基本设置
# =========================
fut_scens = ["ssp126", "ssp585"]
haz = "dry"  # 这里只先分析干旱

THRESH_IMPORT_DEP = 0.60
YEAR_MIN, YEAR_MAX = 2014, 2023

income_order = [
    "Low income",
    "Lower middle income",
    "Upper middle income",
    "High income"
]

crop_items = {
    "Maize (corn)": "maiz",
    "Rice, paddy": "rice",
    "Wheat": "whea",
    "Soybeans": "soyb",
}

# =========================
# 3) 工具函数
# =========================
def to_lon180(da, lon_name="lon"):
    lon = da[lon_name]
    if float(lon.max()) > 180:
        da = da.assign_coords({lon_name: (((lon + 180) % 360) - 180)})
        da = da.sortby(lon_name)
    return da

def crop_area(dist_ds, crop):
    if crop in ["maiz", "soyb", "rice"]:
        return (dist_ds[f"{crop}_Irr"] + dist_ds[f"{crop}_nonIrr"]).fillna(0.0)
    elif crop == "whea":
        # 优先 whea，如全nan则用 wwh+swh
        if ("whea_Irr" in dist_ds) and ("whea_nonIrr" in dist_ds):
            a = (dist_ds["whea_Irr"] + dist_ds["whea_nonIrr"])
            if int(np.isfinite(a.values).sum()) == 0:
                a = (dist_ds["wwh_Irr"] + dist_ds["wwh_nonIrr"] +
                     dist_ds["swh_Irr"] + dist_ds["swh_nonIrr"])
        else:
            a = (dist_ds["wwh_Irr"] + dist_ds["wwh_nonIrr"] +
                 dist_ds["swh_Irr"] + dist_ds["swh_nonIrr"])
        return a.fillna(0.0)
    else:
        raise ValueError(crop)
def add_iso3_from_name(series):
    """
    将国家名稳健转换为 ISO_A3。
    优先用手工字典；若没有，再逐个调用 country_converter。
    可处理返回 list / ndarray / 多匹配 的情况。
    """
    manual = {
        "Türkiye": "TUR",
        "Turkey": "TUR",
        "Viet Nam": "VNM",
        "Vietnam": "VNM",
        "Iran (Islamic Republic of)": "IRN",
        "Iran": "IRN",
        "Russian Federation": "RUS",
        "Russia": "RUS",
        "Bolivia (Plurinational State of)": "BOL",
        "Bolivia": "BOL",
        "Venezuela (Bolivarian Republic of)": "VEN",
        "Venezuela": "VEN",
        "Lao People's Democratic Republic": "LAO",
        "Laos": "LAO",
        "Syrian Arab Republic": "SYR",
        "Syria": "SYR",
        "Czechia": "CZE",
        "Czech Republic": "CZE",
        "Republic of Korea": "KOR",
        "South Korea": "KOR",
        "Democratic People's Republic of Korea": "PRK",
        "North Korea": "PRK",
        "United Republic of Tanzania": "TZA",
        "Tanzania": "TZA",
        "Côte d'Ivoire": "CIV",
        "Cote d'Ivoire": "CIV",
        "Micronesia (Federated States of)": "FSM",
        "Micronesia": "FSM",
        "Moldova, Republic of": "MDA",
        "Moldova": "MDA",
        "United States of America": "USA",
        "United States": "USA",
        "Bahamas, The": "BHS",
        "Gambia, The": "GMB",
        "China, mainland": "CHN",
        "China": "CHN",
        "Hong Kong SAR": "HKG",
        "China, Hong Kong SAR": "HKG",
        "China, Macao SAR": "MAC",
        "China, Taiwan Province of": "TWN",
        "Taiwan Province of China": "TWN",
        "Taiwan": "TWN",
    }

    out = []

    for nm in series.astype(str):
        nm = nm.strip()

        # 1) 优先手工映射
        if nm in manual:
            out.append(manual[nm])
            continue

        # 2) 再尝试 country_converter
        try:
            x = coco.convert(names=nm, to="ISO3", not_found=None)

            # 2a) 如果返回 ndarray / list / tuple，取第一个非空候选
            if isinstance(x, (list, tuple, np.ndarray)):
                vals = [str(v) for v in x if v is not None and str(v) != "not found" and str(v) != "None"]
                if len(vals) > 0:
                    out.append(vals[0])
                else:
                    out.append(np.nan)

            # 2b) 单个标量
            else:
                if x is None:
                    out.append(np.nan)
                else:
                    sx = str(x)
                    if sx in ["not found", "None", "nan"]:
                        out.append(np.nan)
                    else:
                        out.append(sx)

        except Exception:
            out.append(np.nan)

    return pd.Series(out, index=series.index)
def classify_trade_role(net_trade, tol=0):
    if pd.isna(net_trade):
        return np.nan
    if net_trade > tol:
        return "Net exporter"
    elif net_trade < -tol:
        return "Net importer"
    else:
        return "Balanced"

def classify_imp_dep(x, th=0.60):
    if pd.isna(x):
        return np.nan
    return "High import dependence" if x >= th else "Low/medium import dependence"

# =========================
# 4) 读取 loss 数据
# =========================
print("Reading loss nc ...")
ds = xr.open_dataset(nc_fp)
da = ds[varname]   # (scen, haz, lat, lon)
da = to_lon180(da, lon_name="lon")

if da.lat[0] > da.lat[-1]:
    da = da.sortby("lat")

# 未来 - 历史（只取干旱）
delta = da.sel(scen=fut_scens, haz=haz) - da.sel(scen="hist", haz=haz)
# 维度: (scen, lat, lon)

# =========================
# 5) 用 DB['dist'] 作为作物面积权重
# =========================
print("Using DB['dist'] as crop-area weights ...")
dist = DB['dist']
dist = to_lon180(dist, lon_name="lon")
if dist.lat[0] > dist.lat[-1]:
    dist = dist.sortby("lat")

area_total = (
    crop_area(dist, "maiz") +
    crop_area(dist, "soyb") +
    crop_area(dist, "rice") +
    crop_area(dist, "whea")
)

# 对齐到 delta 网格
area_total = area_total.interp(lat=delta.lat, lon=delta.lon, method="nearest").fillna(0.0)
mask_crop = area_total > 0
delta = delta.where(mask_crop)

# =========================
# 6) 读取 shp（ISO_A3）
# =========================
print("Reading country shapefile ...")
gdf = gpd.read_file(shp_fp).to_crs(epsg=4326)

keep_cols = []
for c in ["ISO_A3", "NAME", "ADMIN", "geometry"]:
    if c in gdf.columns:
        keep_cols.append(c)
gdf = gdf[keep_cols].copy()

if "NAME" not in gdf.columns and "ADMIN" in gdf.columns:
    gdf["NAME"] = gdf["ADMIN"]
if "ADMIN" not in gdf.columns and "NAME" in gdf.columns:
    gdf["ADMIN"] = gdf["NAME"]

gdf = gdf[gdf["ISO_A3"].notna()].copy()
gdf = gdf[~gdf["ISO_A3"].isin(["-99"])].copy()

# =========================
# 7) regionmask 映射格点到国家
# =========================
print("Building region mask ...")
gdf = gdf.reset_index(drop=True).copy()
gdf["rid"] = np.arange(len(gdf))

regions = regionmask.from_geopandas(
    gdf,
    names="ISO_A3",
    numbers="rid"
)
mask = regions.mask(delta.isel(scen=0))   # (lat, lon)
num_to_iso = dict(zip(gdf["rid"].values, gdf["ISO_A3"].values))

# =========================
# 8) 国家尺度面积加权聚合
#    country_mean = sum(delta * area_total) / sum(area_total)
# =========================
print("Aggregating drought Δloss to country scale ...")
rows = []

valid_idx = np.unique(mask.values[np.isfinite(mask.values)]).astype(int)

for scen in fut_scens:
    d = delta.sel(scen=scen)
    w = area_total.where(np.isfinite(d))

    for idx in valid_idx:
        iso3 = num_to_iso[idx]
        m = xr.where(mask == idx, 1.0, np.nan)

        ww = w * m
        dd = d * m

        denom = float(ww.sum(skipna=True).values)
        if denom <= 0:
            val = np.nan
        else:
            val = float((dd * ww).sum(skipna=True).values / denom)

        rows.append({
            "ISO_A3": iso3,
            "scen": scen,
            "dloss": val
        })

df_loss_long = pd.DataFrame(rows)

df_loss = df_loss_long.pivot(
    index="ISO_A3",
    columns="scen",
    values="dloss"
).reset_index()

df_loss = df_loss.rename(columns={
    "ssp126": "dloss_ssp126",
    "ssp585": "dloss_ssp585"
})

df_loss = df_loss.merge(
    gdf[["ISO_A3", "NAME"]].drop_duplicates(),
    on="ISO_A3",
    how="left"
)

# =========================
# 9) 读取 FAOSTAT 贸易数据
# =========================
print("Reading FAOSTAT trade data ...")
df_tr = pd.read_csv(faostat_fp)

df_tr = df_tr[["Area", "Element", "Item", "Year", "Value"]].copy()
df_tr["Year"] = pd.to_numeric(df_tr["Year"], errors="coerce")
df_tr["Value"] = pd.to_numeric(df_tr["Value"], errors="coerce")

df_tr = df_tr[(df_tr["Year"] >= YEAR_MIN) & (df_tr["Year"] <= YEAR_MAX)].copy()
df_tr = df_tr[df_tr["Item"].isin(crop_items.keys())].copy()

df_tr["ISO_A3"] = add_iso3_from_name(df_tr["Area"])
df_tr = df_tr[df_tr["ISO_A3"].notna()].copy()

# 每个国家-年份 四作物合计贸易量
tmp = (
    df_tr.groupby(["ISO_A3", "Area", "Element", "Year"], as_index=False)["Value"]
    .sum()
)

# 再对年份求均值
df_trade = (
    tmp.groupby(["ISO_A3", "Area", "Element"], as_index=False)["Value"]
    .mean()
)

df_trade_wide = df_trade.pivot_table(
    index=["ISO_A3", "Area"],
    columns="Element",
    values="Value",
    aggfunc="mean"
).reset_index()

df_trade_wide.columns.name = None

for c in ["Import quantity", "Export quantity"]:
    if c not in df_trade_wide.columns:
        df_trade_wide[c] = np.nan

df_trade_wide = df_trade_wide.rename(columns={
    "Area": "trade_name",
    "Import quantity": "imp_qty",
    "Export quantity": "exp_qty"
})

df_trade_wide["imp_qty"] = df_trade_wide["imp_qty"].fillna(0)
df_trade_wide["exp_qty"] = df_trade_wide["exp_qty"].fillna(0)
df_trade_wide["net_trade"] = df_trade_wide["exp_qty"] - df_trade_wide["imp_qty"]

# 简化版进口依赖
df_trade_wide["imp_dep"] = df_trade_wide["imp_qty"] / (
    df_trade_wide["imp_qty"] + df_trade_wide["exp_qty"] + 1e-12
)

df_trade_wide["trade_role"] = df_trade_wide["net_trade"].apply(classify_trade_role)
df_trade_wide["imp_dep_cls"] = df_trade_wide["imp_dep"].apply(
    lambda x: classify_imp_dep(x, th=THRESH_IMPORT_DEP)
)

# =========================
# 10) 读取收入组
# =========================
print("Reading income groups ...")
df_inc = pd.read_excel(income_fp)

df_inc = df_inc.rename(columns={
    "Economy": "econ_name",
    "Code": "ISO_A3",
    "Income group": "income_group"
})

df_inc = df_inc[["ISO_A3", "econ_name", "income_group"]].copy()
df_inc = df_inc[df_inc["ISO_A3"].notna()].copy()

# =========================
# 11) 合并国家表
# =========================
print("Merging all tables ...")
df_country = (
    df_loss.merge(df_trade_wide, on="ISO_A3", how="left")
           .merge(df_inc, on="ISO_A3", how="left")
)

df_country["country"] = (
    df_country["NAME"]
    .fillna(df_country["econ_name"])
    .fillna(df_country["trade_name"])
)

df_country["income_group"] = pd.Categorical(
    df_country["income_group"],
    categories=income_order,
    ordered=True
)

# 定义高风险阈值（按 SSP5-8.5 的上四分位）
q75 = df_country["dloss_ssp585"].quantile(0.75)
df_country["risk_cls_ssp585"] = np.where(
    df_country["dloss_ssp585"] >= q75,
    "High drought risk",
    "Lower drought risk"
)

# =========================
# 12) 导出主表
# =========================
csv_main = os.path.join(out_dir, "country_drought_trade_income_table.csv")
df_country.to_csv(csv_main, index=False, encoding="utf-8-sig")
print("Saved:", csv_main)

# =========================
# 13) 汇总表
# =========================
df_income_sum = (
    df_country.groupby("income_group", observed=False)
    .agg(
        n_country=("ISO_A3", "count"),
        mean_dloss_126=("dloss_ssp126", "mean"),
        mean_dloss_585=("dloss_ssp585", "mean"),
        mean_imp_dep=("imp_dep", "mean"),
        pct_net_importer=("trade_role", lambda x: np.mean(x == "Net importer") * 100),
        pct_high_imp_dep=("imp_dep_cls", lambda x: np.mean(x == "High import dependence") * 100),
    )
    .reset_index()
)
csv_sum1 = os.path.join(out_dir, "summary_by_income_group.csv")
df_income_sum.to_csv(csv_sum1, index=False, encoding="utf-8-sig")
print("Saved:", csv_sum1)

df_trade_sum = (
    df_country.groupby("trade_role", dropna=False)
    .agg(
        n_country=("ISO_A3", "count"),
        mean_dloss_126=("dloss_ssp126", "mean"),
        mean_dloss_585=("dloss_ssp585", "mean"),
        mean_imp_dep=("imp_dep", "mean"),
    )
    .reset_index()
)
csv_sum2 = os.path.join(out_dir, "summary_by_trade_role.csv")
df_trade_sum.to_csv(csv_sum2, index=False, encoding="utf-8-sig")
print("Saved:", csv_sum2)

# =========================
# 14) 合并回 shp 准备绘图
# =========================
gdf_plot = gdf.merge(df_country, on="ISO_A3", how="left")

# =========================
# 15) 图1：国家尺度干旱增量地图（SSP5-8.5）
# =========================
fig, ax = plt.subplots(1, 1, figsize=(12, 6))

bounds = [0, 0.5, 1, 3, 5, 10]
cmap = mcolors.ListedColormap(["#d9d9d9", "#fee5d9", "#fcae91", "#fb6a4a", "#cb181d"])
norm = mcolors.BoundaryNorm(bounds, ncolors=cmap.N, clip=False)

gdf_plot.plot(
    column="dloss_ssp585",
    ax=ax,
    cmap=cmap,
    norm=norm,
    linewidth=0.25,
    edgecolor="white",
    missing_kwds={"color": "lightgrey", "label": "No data"}
)

ax.set_axis_off()
ax.set_title("Country-scale increase in drought-induced yield loss (SSP5-8.5 − Historical)", fontsize=13)

sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm._A = []
cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.04, pad=0.03)
cbar.set_label("Δ drought-induced yield loss (%)")

map_pdf = os.path.join(out_dir, "Map_country_dloss_ssp585.pdf")
fig.savefig(map_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved:", map_pdf)

# =========================
# 16) 图2：收入组箱线图
# =========================
fig, ax = plt.subplots(1, 1, figsize=(8, 5))

tmp_box = df_country.dropna(subset=["income_group", "dloss_ssp585"]).copy()

sns.boxplot(
    data=tmp_box,
    x="income_group",
    y="dloss_ssp585",
    ax=ax,
    color="#fcae91",
    width=0.6,
    fliersize=2
)
sns.stripplot(
    data=tmp_box,
    x="income_group",
    y="dloss_ssp585",
    ax=ax,
    color="0.25",
    size=2.2,
    alpha=0.5,
    jitter=0.2
)

ax.set_xlabel("")
ax.set_ylabel("Δ drought-induced yield loss (%)")
ax.set_title("Drought-risk increase by income group (SSP5-8.5)")
ax.tick_params(axis="x", rotation=20)

box_pdf = os.path.join(out_dir, "Box_income_vs_dloss_ssp585.pdf")
fig.savefig(box_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved:", box_pdf)

# =========================
# 17) 图3：进口依赖 vs 干旱风险
# =========================
fig, ax = plt.subplots(1, 1, figsize=(8, 6))

tmp_sc = df_country.dropna(subset=["imp_dep", "dloss_ssp585"]).copy()

palette = {
    "Low income": "#8c510a",
    "Lower middle income": "#d8b365",
    "Upper middle income": "#5ab4ac",
    "High income": "#01665e"
}

for grp in income_order:
    dsub = tmp_sc[tmp_sc["income_group"] == grp]
    if len(dsub) == 0:
        continue

    ax.scatter(
        dsub["imp_dep"],
        dsub["dloss_ssp585"],
        s=np.sqrt(dsub["imp_qty"].fillna(0) + 1) * 1.2,
        c=palette.get(grp, "0.5"),
        alpha=0.75,
        edgecolors="white",
        linewidths=0.4,
        label=grp
    )

ax.axvline(THRESH_IMPORT_DEP, color="0.5", ls="--", lw=1)
ax.axhline(q75, color="0.5", ls="--", lw=1)

ax.set_xlabel("Import dependence ratio")
ax.set_ylabel("Δ drought-induced yield loss (%)")
ax.set_title("Import dependence versus drought-risk increase (SSP5-8.5)")
ax.legend(frameon=False, fontsize=9, loc="best")

sc_pdf = os.path.join(out_dir, "Scatter_impdep_vs_dloss_ssp585.pdf")
fig.savefig(sc_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved:", sc_pdf)

# =========================
# 18) 图4：净进口/净出口 × 收入组
# =========================
tmp_bar = (
    df_country.dropna(subset=["income_group", "trade_role", "dloss_ssp585"])
    .groupby(["income_group", "trade_role"], observed=False)["dloss_ssp585"]
    .mean()
    .reset_index()
)

fig, ax = plt.subplots(1, 1, figsize=(9, 5))

sns.barplot(
    data=tmp_bar,
    x="income_group",
    y="dloss_ssp585",
    hue="trade_role",
    ax=ax,
    palette={
        "Net importer": "#ef8a62",
        "Net exporter": "#67a9cf",
        "Balanced": "#bdbdbd"
    }
)

ax.set_xlabel("")
ax.set_ylabel("Mean Δ drought-induced yield loss (%)")
ax.set_title("Drought-risk increase by trade role and income group")
ax.legend(title="", frameon=False)
ax.tick_params(axis="x", rotation=20)

bar_pdf = os.path.join(out_dir, "Bar_trade_role_income_dloss_ssp585.pdf")
fig.savefig(bar_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved:", bar_pdf)

# =========================
# 19) 输出热点国家
#     高进口依赖 + 高干旱风险
# =========================
df_hot = df_country[
    (df_country["imp_dep"] >= THRESH_IMPORT_DEP) &
    (df_country["dloss_ssp585"] >= q75)
].copy()

df_hot = df_hot.sort_values(["dloss_ssp585", "imp_dep"], ascending=[False, False])

csv_hot = os.path.join(out_dir, "hotspot_high_impdep_high_droughtrisk.csv")
df_hot[[
    "ISO_A3", "country", "income_group",
    "imp_qty", "exp_qty", "net_trade", "imp_dep",
    "trade_role", "dloss_ssp126", "dloss_ssp585"
]].to_csv(csv_hot, index=False, encoding="utf-8-sig")
print("Saved:", csv_hot)

print("\nDone.")

# %% cell 31
# -*- coding: utf-8 -*-
import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import seaborn as sns

# 需要安装：
# pip install regionmask country_converter openpyxl
import regionmask
import country_converter as coco

# =========================
# 0) 全局绘图参数
# =========================
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 1) 路径配置
# =========================
nc_fp = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\loss_cropmean_area_weighted_3scen_2haz.nc"
varname = "loss_pct_cropmean_areaW"

faostat_fp = r"D:\Edge_download\FAOSTAT_data_en_3-8-2026.csv"
income_fp  = r"D:\Edge_download\CLASS_2025_10_07.xlsx"
shp_fp     = r"D:\AAUDE\paper_v2\paper2\data\ne_10m_admin_0_countries\ne_10m_admin_0_countries.shp"

out_dir = r"D:\AAUDE\paper_v2\paper4\outputs_country_trade_drought_fix"
os.makedirs(out_dir, exist_ok=True)

# =========================
# 2) 基本设置
# =========================
fut_scens = ["ssp126", "ssp585"]
haz = "dry"  # 这里只先分析干旱

THRESH_IMPORT_DEP = 0.60
YEAR_MIN, YEAR_MAX = 2014, 2023

# 面积阈值：小于这个值视为无有效四作物面积
# 这个阈值需要根据你的 dist 数据量纲适当调整
# 若 dist 是 harvested area fraction，可先设 1e-8 ~ 1e-6
area_eps = 1e-8

income_order = [
    "Low income",
    "Lower middle income",
    "Upper middle income",
    "High income"
]

crop_items = {
    "Maize (corn)": "maiz",
    "Rice, paddy": "rice",
    "Wheat": "whea",
    "Soybeans": "soyb",
}

# =========================
# 3) 工具函数
# =========================
def to_lon180(da, lon_name="lon"):
    lon = da[lon_name]
    if float(lon.max()) > 180:
        da = da.assign_coords({lon_name: (((lon + 180) % 360) - 180)})
        da = da.sortby(lon_name)
    return da

def crop_area(dist_ds, crop):
    if crop in ["maiz", "soyb", "rice"]:
        return (dist_ds[f"{crop}_Irr"] + dist_ds[f"{crop}_nonIrr"]).fillna(0.0)
    elif crop == "whea":
        # 优先 whea，如全nan则用 wwh+swh
        if ("whea_Irr" in dist_ds) and ("whea_nonIrr" in dist_ds):
            a = (dist_ds["whea_Irr"] + dist_ds["whea_nonIrr"])
            if int(np.isfinite(a.values).sum()) == 0:
                a = (dist_ds["wwh_Irr"] + dist_ds["wwh_nonIrr"] +
                     dist_ds["swh_Irr"] + dist_ds["swh_nonIrr"])
        else:
            a = (dist_ds["wwh_Irr"] + dist_ds["wwh_nonIrr"] +
                 dist_ds["swh_Irr"] + dist_ds["swh_nonIrr"])
        return a.fillna(0.0)
    else:
        raise ValueError(crop)

def add_iso3_from_name(series):
    manual = {
        "Türkiye": "TUR",
        "Turkey": "TUR",
        "Viet Nam": "VNM",
        "Vietnam": "VNM",
        "Iran (Islamic Republic of)": "IRN",
        "Iran": "IRN",
        "Russian Federation": "RUS",
        "Russia": "RUS",
        "Bolivia (Plurinational State of)": "BOL",
        "Bolivia": "BOL",
        "Venezuela (Bolivarian Republic of)": "VEN",
        "Venezuela": "VEN",
        "Lao People's Democratic Republic": "LAO",
        "Laos": "LAO",
        "Syrian Arab Republic": "SYR",
        "Syria": "SYR",
        "Czechia": "CZE",
        "Czech Republic": "CZE",
        "Republic of Korea": "KOR",
        "South Korea": "KOR",
        "Democratic People's Republic of Korea": "PRK",
        "North Korea": "PRK",
        "United Republic of Tanzania": "TZA",
        "Tanzania": "TZA",
        "Côte d'Ivoire": "CIV",
        "Cote d'Ivoire": "CIV",
        "Micronesia (Federated States of)": "FSM",
        "Micronesia": "FSM",
        "Moldova, Republic of": "MDA",
        "Moldova": "MDA",
        "United States of America": "USA",
        "United States": "USA",
        "Bahamas, The": "BHS",
        "Gambia, The": "GMB",
        "China, mainland": "CHN",
        "China": "CHN",
        "Hong Kong SAR": "HKG",
        "China, Hong Kong SAR": "HKG",
        "China, Macao SAR": "MAC",
        "China, Taiwan Province of": "TWN",
        "Taiwan Province of China": "TWN",
        "Taiwan": "TWN",
    }

    out = []
    for nm in series.astype(str):
        nm = nm.strip()

        if nm in manual:
            out.append(manual[nm])
            continue

        try:
            x = coco.convert(names=nm, to="ISO3", not_found=None)

            if isinstance(x, (list, tuple, np.ndarray)):
                vals = [str(v) for v in x if v is not None and str(v) not in ["not found", "None", "nan"]]
                out.append(vals[0] if len(vals) > 0 else np.nan)
            else:
                if x is None:
                    out.append(np.nan)
                else:
                    sx = str(x)
                    out.append(np.nan if sx in ["not found", "None", "nan"] else sx)

        except Exception:
            out.append(np.nan)

    return pd.Series(out, index=series.index)

def classify_trade_role(net_trade, tol=0):
    if pd.isna(net_trade):
        return np.nan
    if net_trade > tol:
        return "Net exporter"
    elif net_trade < -tol:
        return "Net importer"
    else:
        return "Balanced"

def classify_imp_dep(x, th=0.60):
    if pd.isna(x):
        return np.nan
    return "High import dependence" if x >= th else "Low/medium import dependence"

# =========================
# 4) 读取 loss 数据
# =========================
print("Reading loss nc ...")
ds = xr.open_dataset(nc_fp)
da = ds[varname]   # (scen, haz, lat, lon)
da = to_lon180(da, lon_name="lon")

if da.lat[0] > da.lat[-1]:
    da = da.sortby("lat")

# 未来 - 历史（只取干旱）
delta = da.sel(scen=fut_scens, haz=haz) - da.sel(scen="hist", haz=haz)
# 注意：这里不再提前 mask_crop
# dims: (scen, lat, lon)

# =========================
# 5) 用 DB['dist'] 作为作物面积权重
# =========================
print("Using DB['dist'] as crop-area weights ...")
dist = DB['dist']
dist = to_lon180(dist, lon_name="lon")
if dist.lat[0] > dist.lat[-1]:
    dist = dist.sortby("lat")

area_total = (
    crop_area(dist, "maiz") +
    crop_area(dist, "soyb") +
    crop_area(dist, "rice") +
    crop_area(dist, "whea")
)

area_total = area_total.interp(lat=delta.lat, lon=delta.lon, method="nearest").fillna(0.0)

print("Area stats:")
print("  min :", float(area_total.min().values))
print("  max :", float(area_total.max().values))
print("  mean:", float(area_total.mean().values))

# =========================
# 6) 读取 shp（ISO_A3）
# =========================
print("Reading country shapefile ...")
gdf = gpd.read_file(shp_fp).to_crs(epsg=4326)

keep_cols = []
for c in ["ISO_A3", "NAME", "ADMIN", "geometry"]:
    if c in gdf.columns:
        keep_cols.append(c)
gdf = gdf[keep_cols].copy()

if "NAME" not in gdf.columns and "ADMIN" in gdf.columns:
    gdf["NAME"] = gdf["ADMIN"]
if "ADMIN" not in gdf.columns and "NAME" in gdf.columns:
    gdf["ADMIN"] = gdf["NAME"]

gdf = gdf[gdf["ISO_A3"].notna()].copy()
gdf = gdf[~gdf["ISO_A3"].isin(["-99"])].copy()

# 全部国家底表
gdf = gdf.reset_index(drop=True).copy()
gdf["rid"] = np.arange(len(gdf))

df_base = gdf[["ISO_A3", "NAME"]].drop_duplicates().copy()
df_base = df_base.rename(columns={"NAME": "shp_name"})

# =========================
# 7) regionmask 映射格点到国家
# =========================
print("Building region mask ...")
regions = regionmask.from_geopandas(
    gdf,
    names="ISO_A3",
    numbers="rid"
)

mask = regions.mask(delta.isel(scen=0))   # (lat, lon)
num_to_iso = dict(zip(gdf["rid"].values, gdf["ISO_A3"].values))

# =========================
# 8) 国家尺度面积加权聚合（遍历全部国家）
#    注意：不再提前过滤 delta，而是在聚合时用 area_total 作为权重
# =========================
print("Aggregating drought Δloss to country scale ...")
rows = []

for idx in gdf["rid"].values:
    iso3 = num_to_iso[idx]
    m = xr.where(mask == idx, 1.0, np.nan)

    # 国家总四作物面积
    crop_area_sum = float((area_total * m).sum(skipna=True).values)

    row = {
        "ISO_A3": iso3,
        "crop_area_sum": crop_area_sum,
        "has_crop": int(crop_area_sum > area_eps)
    }

    for scen in fut_scens:
        d = delta.sel(scen=scen)

        # 国家内权重
        ww = area_total * m
        # 国家内损失
        dd = d * m

        # 只有损失值有效时才进入分子
        num = float((dd * ww).sum(skipna=True).values)
        den = float(ww.where(np.isfinite(dd)).sum(skipna=True).values)

        if (crop_area_sum <= area_eps) or (den <= area_eps):
            val = np.nan
        else:
            val = num / den

        row[f"dloss_{scen}"] = val

    rows.append(row)

df_loss = pd.DataFrame(rows)

df_loss = df_loss.merge(
    gdf[["ISO_A3", "NAME"]].drop_duplicates(),
    on="ISO_A3",
    how="left"
)

# =========================
# 9) 读取 FAOSTAT 贸易数据
# =========================
print("Reading FAOSTAT trade data ...")
df_tr = pd.read_csv(faostat_fp)

df_tr = df_tr[["Area", "Element", "Item", "Year", "Value"]].copy()
df_tr["Year"] = pd.to_numeric(df_tr["Year"], errors="coerce")
df_tr["Value"] = pd.to_numeric(df_tr["Value"], errors="coerce")

df_tr = df_tr[(df_tr["Year"] >= YEAR_MIN) & (df_tr["Year"] <= YEAR_MAX)].copy()
df_tr = df_tr[df_tr["Item"].isin(crop_items.keys())].copy()

df_tr["ISO_A3"] = add_iso3_from_name(df_tr["Area"])
df_tr = df_tr[df_tr["ISO_A3"].notna()].copy()

# 每个国家-年份：四作物合计贸易量
tmp = (
    df_tr.groupby(["ISO_A3", "Area", "Element", "Year"], as_index=False)["Value"]
    .sum()
)

# 再对年份求均值
df_trade = (
    tmp.groupby(["ISO_A3", "Area", "Element"], as_index=False)["Value"]
    .mean()
)

df_trade_wide = df_trade.pivot_table(
    index=["ISO_A3", "Area"],
    columns="Element",
    values="Value",
    aggfunc="mean"
).reset_index()

df_trade_wide.columns.name = None

for c in ["Import quantity", "Export quantity"]:
    if c not in df_trade_wide.columns:
        df_trade_wide[c] = np.nan

df_trade_wide = df_trade_wide.rename(columns={
    "Area": "trade_name",
    "Import quantity": "imp_qty",
    "Export quantity": "exp_qty"
})

df_trade_wide["imp_qty"] = df_trade_wide["imp_qty"].fillna(0)
df_trade_wide["exp_qty"] = df_trade_wide["exp_qty"].fillna(0)
df_trade_wide["net_trade"] = df_trade_wide["exp_qty"] - df_trade_wide["imp_qty"]

# 简化版进口依赖
df_trade_wide["imp_dep"] = df_trade_wide["imp_qty"] / (
    df_trade_wide["imp_qty"] + df_trade_wide["exp_qty"] + 1e-12
)

df_trade_wide["trade_role"] = df_trade_wide["net_trade"].apply(classify_trade_role)
df_trade_wide["imp_dep_cls"] = df_trade_wide["imp_dep"].apply(
    lambda x: classify_imp_dep(x, th=THRESH_IMPORT_DEP)
)

# =========================
# 10) 读取收入组
# =========================
print("Reading income groups ...")
df_inc = pd.read_excel(income_fp)

df_inc = df_inc.rename(columns={
    "Economy": "econ_name",
    "Code": "ISO_A3",
    "Income group": "income_group"
})

df_inc = df_inc[["ISO_A3", "econ_name", "income_group"]].copy()
df_inc = df_inc[df_inc["ISO_A3"].notna()].copy()

# =========================
# 11) 合并国家表（以全部国家为底）
# =========================
print("Merging all tables ...")
df_country = (
    df_base
    .merge(df_loss.drop(columns=["NAME"], errors="ignore"), on="ISO_A3", how="left")
    .merge(df_trade_wide, on="ISO_A3", how="left")
    .merge(df_inc, on="ISO_A3", how="left")
)

df_country["country"] = (
    df_country["shp_name"]
    .fillna(df_country.get("econ_name"))
    .fillna(df_country.get("trade_name"))
)

df_country["crop_area_sum"] = df_country["crop_area_sum"].fillna(0)
df_country["has_crop"] = df_country["has_crop"].fillna(0).astype(int)

df_country["income_group"] = pd.Categorical(
    df_country["income_group"],
    categories=income_order,
    ordered=True
)

# 定义高风险阈值（仅基于有有效估计国家）
q75 = df_country["dloss_ssp585"].dropna().quantile(0.75)

df_country["risk_cls_ssp585"] = np.where(
    df_country["dloss_ssp585"].notna() & (df_country["dloss_ssp585"] >= q75),
    "High drought risk",
    np.where(
        df_country["dloss_ssp585"].notna(),
        "Lower drought risk",
        "No crop / no estimate"
    )
)

# =========================
# 12) 导出主表
# =========================
csv_main = os.path.join(out_dir, "country_drought_trade_income_table_fix.csv")
df_country.to_csv(csv_main, index=False, encoding="utf-8-sig")
print("Saved:", csv_main)

# =========================
# 13) 汇总表
# =========================
df_income_sum = (
    df_country.groupby("income_group", observed=False)
    .agg(
        n_country_total=("ISO_A3", "count"),
        n_country_with_crop=("has_crop", "sum"),
        n_country_with_est=("dloss_ssp585", lambda x: x.notna().sum()),
        mean_dloss_126=("dloss_ssp126", "mean"),
        mean_dloss_585=("dloss_ssp585", "mean"),
        mean_imp_dep=("imp_dep", "mean"),
        pct_net_importer=("trade_role", lambda x: np.mean(x == "Net importer") * 100),
        pct_high_imp_dep=("imp_dep_cls", lambda x: np.mean(x == "High import dependence") * 100),
    )
    .reset_index()
)
csv_sum1 = os.path.join(out_dir, "summary_by_income_group_fix.csv")
df_income_sum.to_csv(csv_sum1, index=False, encoding="utf-8-sig")
print("Saved:", csv_sum1)

df_trade_sum = (
    df_country.groupby("trade_role", dropna=False)
    .agg(
        n_country_total=("ISO_A3", "count"),
        n_country_with_crop=("has_crop", "sum"),
        n_country_with_est=("dloss_ssp585", lambda x: x.notna().sum()),
        mean_dloss_126=("dloss_ssp126", "mean"),
        mean_dloss_585=("dloss_ssp585", "mean"),
        mean_imp_dep=("imp_dep", "mean"),
    )
    .reset_index()
)
csv_sum2 = os.path.join(out_dir, "summary_by_trade_role_fix.csv")
df_trade_sum.to_csv(csv_sum2, index=False, encoding="utf-8-sig")
print("Saved:", csv_sum2)

# =========================
# 14) 合并回 shp 准备绘图
# =========================
gdf_plot = gdf.merge(df_country, on="ISO_A3", how="left")

# 地图辅助分类
gdf_plot["map_cls"] = "No crop area / no estimate"
gdf_plot.loc[gdf_plot["has_crop"] == 1, "map_cls"] = "Has crop area"
gdf_plot.loc[gdf_plot["dloss_ssp585"].notna(), "map_cls"] = "Has estimate"

# =========================
# 15) 图1：国家尺度干旱增量地图（SSP5-8.5）
# =========================
fig, ax = plt.subplots(1, 1, figsize=(12, 6))

# 底色：所有国家
gdf_plot.plot(
    ax=ax,
    color="#d9d9d9",
    linewidth=0.25,
    edgecolor="white"
)

# 有估计国家再覆盖
gdf_est = gdf_plot[gdf_plot["dloss_ssp585"].notna()].copy()

bounds = [0, 0.5, 1, 3, 5, 10]
cmap = mcolors.ListedColormap(["#d9d9d9", "#fee5d9", "#fcae91", "#fb6a4a", "#cb181d"])
norm = mcolors.BoundaryNorm(bounds, ncolors=cmap.N, clip=False)

gdf_est.plot(
    column="dloss_ssp585",
    ax=ax,
    cmap=cmap,
    norm=norm,
    linewidth=0.25,
    edgecolor="white"
)

ax.set_axis_off()
ax.set_title("Country-scale increase in drought-induced yield loss (SSP5-8.5 − Historical)", fontsize=13)

sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm._A = []
cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.04, pad=0.03)
cbar.set_label("Δ drought-induced yield loss (%)")

map_pdf = os.path.join(out_dir, "Map_country_dloss_ssp585_fix.pdf")
fig.savefig(map_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved:", map_pdf)

# =========================
# 16) 图2：收入组箱线图（仅有估计国家）
# =========================
fig, ax = plt.subplots(1, 1, figsize=(8, 5))

tmp_box = df_country.dropna(subset=["income_group", "dloss_ssp585"]).copy()

sns.boxplot(
    data=tmp_box,
    x="income_group",
    y="dloss_ssp585",
    ax=ax,
    color="#fcae91",
    width=0.6,
    fliersize=2
)
sns.stripplot(
    data=tmp_box,
    x="income_group",
    y="dloss_ssp585",
    ax=ax,
    color="0.25",
    size=2.2,
    alpha=0.5,
    jitter=0.2
)

ax.set_xlabel("")
ax.set_ylabel("Δ drought-induced yield loss (%)")
ax.set_title("Drought-risk increase by income group (SSP5-8.5)")
ax.tick_params(axis="x", rotation=20)

box_pdf = os.path.join(out_dir, "Box_income_vs_dloss_ssp585_fix.pdf")
fig.savefig(box_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved:", box_pdf)

# =========================
# 17) 图3：进口依赖 vs 干旱风险（仅有估计国家）
# =========================
fig, ax = plt.subplots(1, 1, figsize=(8, 6))

tmp_sc = df_country.dropna(subset=["imp_dep", "dloss_ssp585"]).copy()

palette = {
    "Low income": "#8c510a",
    "Lower middle income": "#d8b365",
    "Upper middle income": "#5ab4ac",
    "High income": "#01665e"
}

for grp in income_order:
    dsub = tmp_sc[tmp_sc["income_group"] == grp]
    if len(dsub) == 0:
        continue

    ax.scatter(
        dsub["imp_dep"],
        dsub["dloss_ssp585"],
        s=np.sqrt(dsub["imp_qty"].fillna(0) + 1) * 1.2,
        c=palette.get(grp, "0.5"),
        alpha=0.75,
        edgecolors="white",
        linewidths=0.4,
        label=grp
    )

ax.axvline(THRESH_IMPORT_DEP, color="0.5", ls="--", lw=1)
ax.axhline(q75, color="0.5", ls="--", lw=1)

ax.set_xlabel("Import dependence ratio")
ax.set_ylabel("Δ drought-induced yield loss (%)")
ax.set_title("Import dependence versus drought-risk increase (SSP5-8.5)")
ax.legend(frameon=False, fontsize=9, loc="best")

sc_pdf = os.path.join(out_dir, "Scatter_impdep_vs_dloss_ssp585_fix.pdf")
fig.savefig(sc_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved:", sc_pdf)

# =========================
# 18) 图4：净进口/净出口 × 收入组（仅有估计国家）
# =========================
tmp_bar = (
    df_country.dropna(subset=["income_group", "trade_role", "dloss_ssp585"])
    .groupby(["income_group", "trade_role"], observed=False)["dloss_ssp585"]
    .mean()
    .reset_index()
)

fig, ax = plt.subplots(1, 1, figsize=(9, 5))

sns.barplot(
    data=tmp_bar,
    x="income_group",
    y="dloss_ssp585",
    hue="trade_role",
    ax=ax,
    palette={
        "Net importer": "#ef8a62",
        "Net exporter": "#67a9cf",
        "Balanced": "#bdbdbd"
    }
)

ax.set_xlabel("")
ax.set_ylabel("Mean Δ drought-induced yield loss (%)")
ax.set_title("Drought-risk increase by trade role and income group")
ax.legend(title="", frameon=False)
ax.tick_params(axis="x", rotation=20)

bar_pdf = os.path.join(out_dir, "Bar_trade_role_income_dloss_ssp585_fix.pdf")
fig.savefig(bar_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved:", bar_pdf)

# =========================
# 19) 输出热点国家
# =========================
df_hot = df_country[
    (df_country["imp_dep"] >= THRESH_IMPORT_DEP) &
    (df_country["dloss_ssp585"] >= q75)
].copy()

df_hot = df_hot.sort_values(["dloss_ssp585", "imp_dep"], ascending=[False, False])

csv_hot = os.path.join(out_dir, "hotspot_high_impdep_high_droughtrisk_fix.csv")
df_hot[[
    "ISO_A3", "country", "income_group",
    "crop_area_sum", "has_crop",
    "imp_qty", "exp_qty", "net_trade", "imp_dep",
    "trade_role", "dloss_ssp126", "dloss_ssp585"
]].to_csv(csv_hot, index=False, encoding="utf-8-sig")
print("Saved:", csv_hot)

# =========================
# 20) 调试输出
# =========================
print("\n===== CHECKS =====")
print("Total countries in shp:", len(df_base))
print("Total countries in df_country:", len(df_country))
print("Countries with crop area:", int((df_country["has_crop"] == 1).sum()))
print("Countries with drought estimate (SSP5-8.5):", int(df_country["dloss_ssp585"].notna().sum()))
print("Countries without crop area:", int((df_country["has_crop"] == 0).sum()))

print("\nCountries with crop but no estimate:")
tmp_noest = df_country[(df_country["has_crop"] == 1) & (df_country["dloss_ssp585"].isna())]
print(tmp_noest[["ISO_A3", "country", "crop_area_sum"]].head(20))

print("\nSmallest crop-area countries with estimate:")
tmp_small = df_country[(df_country["has_crop"] == 1) & (df_country["dloss_ssp585"].notna())] \
    .sort_values("crop_area_sum").head(20)
print(tmp_small[["ISO_A3", "country", "crop_area_sum", "dloss_ssp585"]])

# 可选：输出 area_total 分布诊断
try:
    fig, ax = plt.subplots(1, 1, figsize=(7, 4))
    vals = area_total.values.flatten()
    vals = vals[np.isfinite(vals)]
    vals = vals[vals > 0]
    ax.hist(vals, bins=100)
    ax.set_xlabel("area_total (>0)")
    ax.set_ylabel("Count")
    ax.set_title("Distribution of crop-area weights")
    hist_pdf = os.path.join(out_dir, "Hist_area_total.pdf")
    fig.savefig(hist_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", hist_pdf)
except Exception as e:
    print("Histogram export skipped:", e)

print("\nDone.")

# %% cell 32


# %% cell 33


# %% cell 34
# -*- coding: utf-8 -*-
import os
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import seaborn as sns

# 需要安装：
# pip install regionmask country_converter openpyxl
import regionmask
import country_converter as coco

# =========================
# 0) 全局绘图参数
# =========================
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 1) 路径配置
# =========================
nc_fp = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\loss_cropmean_area_weighted_3scen_2haz.nc"
varname = "loss_pct_cropmean_areaW"

faostat_fp = r"D:\Edge_download\FAOSTAT_data_en_3-8-2026.csv"
income_fp  = r"D:\Edge_download\CLASS_2025_10_07.xlsx"
shp_fp     = r"D:\AAUDE\paper_v2\paper2\data\ne_10m_admin_0_countries\ne_10m_admin_0_countries.shp"

out_dir = r"D:\AAUDE\paper_v2\paper4\outputs_country_trade_2scen_2haz"
os.makedirs(out_dir, exist_ok=True)

# =========================
# 2) 基本设置
# =========================
fut_scens = ["ssp126", "ssp585"]
hazards = ["dry", "wet"]
haz_labels = {"dry": "Drought", "wet": "Waterlogging"}
scen_labels = {"ssp126": "SSP1-2.6", "ssp585": "SSP5-8.5"}

THRESH_IMPORT_DEP = 0.60
YEAR_MIN, YEAR_MAX = 2014, 2023
area_eps = 1e-8

income_order = [
    "Low income",
    "Lower middle income",
    "Upper middle income",
    "High income"
]

crop_items = {
    "Maize (corn)": "maiz",
    "Rice, paddy": "rice",
    "Wheat": "whea",
    "Soybeans": "soyb",
}

haz_color_main = {"dry": "#ef8a62", "wet": "#67a9cf"}

# =========================
# 3) 工具函数
# =========================
def to_lon180(da, lon_name="lon"):
    lon = da[lon_name]
    if float(lon.max()) > 180:
        da = da.assign_coords({lon_name: (((lon + 180) % 360) - 180)})
        da = da.sortby(lon_name)
    return da

def crop_area(dist_ds, crop):
    if crop in ["maiz", "soyb", "rice"]:
        return (dist_ds[f"{crop}_Irr"] + dist_ds[f"{crop}_nonIrr"]).fillna(0.0)
    elif crop == "whea":
        if ("whea_Irr" in dist_ds) and ("whea_nonIrr" in dist_ds):
            a = (dist_ds["whea_Irr"] + dist_ds["whea_nonIrr"])
            if int(np.isfinite(a.values).sum()) == 0:
                a = (dist_ds["wwh_Irr"] + dist_ds["wwh_nonIrr"] +
                     dist_ds["swh_Irr"] + dist_ds["swh_nonIrr"])
        else:
            a = (dist_ds["wwh_Irr"] + dist_ds["wwh_nonIrr"] +
                 dist_ds["swh_Irr"] + dist_ds["swh_nonIrr"])
        return a.fillna(0.0)
    else:
        raise ValueError(crop)

def add_iso3_from_name(series):
    manual = {
        "Türkiye": "TUR",
        "Turkey": "TUR",
        "Viet Nam": "VNM",
        "Vietnam": "VNM",
        "Iran (Islamic Republic of)": "IRN",
        "Iran": "IRN",
        "Russian Federation": "RUS",
        "Russia": "RUS",
        "Bolivia (Plurinational State of)": "BOL",
        "Bolivia": "BOL",
        "Venezuela (Bolivarian Republic of)": "VEN",
        "Venezuela": "VEN",
        "Lao People's Democratic Republic": "LAO",
        "Laos": "LAO",
        "Syrian Arab Republic": "SYR",
        "Syria": "SYR",
        "Czechia": "CZE",
        "Czech Republic": "CZE",
        "Republic of Korea": "KOR",
        "South Korea": "KOR",
        "Democratic People's Republic of Korea": "PRK",
        "North Korea": "PRK",
        "United Republic of Tanzania": "TZA",
        "Tanzania": "TZA",
        "Côte d'Ivoire": "CIV",
        "Cote d'Ivoire": "CIV",
        "Micronesia (Federated States of)": "FSM",
        "Micronesia": "FSM",
        "Moldova, Republic of": "MDA",
        "Moldova": "MDA",
        "United States of America": "USA",
        "United States": "USA",
        "Bahamas, The": "BHS",
        "Gambia, The": "GMB",
        "China, mainland": "CHN",
        "China": "CHN",
        "Hong Kong SAR": "HKG",
        "China, Hong Kong SAR": "HKG",
        "China, Macao SAR": "MAC",
        "China, Taiwan Province of": "TWN",
        "Taiwan Province of China": "TWN",
        "Taiwan": "TWN",
    }

    out = []
    for nm in series.astype(str):
        nm = nm.strip()

        if nm in manual:
            out.append(manual[nm])
            continue

        try:
            x = coco.convert(names=nm, to="ISO3", not_found=None)

            if isinstance(x, (list, tuple, np.ndarray)):
                vals = [str(v) for v in x if v is not None and str(v) not in ["not found", "None", "nan"]]
                out.append(vals[0] if len(vals) > 0 else np.nan)
            else:
                if x is None:
                    out.append(np.nan)
                else:
                    sx = str(x)
                    out.append(np.nan if sx in ["not found", "None", "nan"] else sx)

        except Exception:
            out.append(np.nan)

    return pd.Series(out, index=series.index)

def classify_trade_role(net_trade, tol=0):
    if pd.isna(net_trade):
        return np.nan
    if net_trade > tol:
        return "Net exporter"
    elif net_trade < -tol:
        return "Net importer"
    else:
        return "Balanced"

def classify_imp_dep(x, th=0.60):
    if pd.isna(x):
        return np.nan
    return "High import dependence" if x >= th else "Low/medium import dependence"

def save_boxplot(df_country, out_fp, haz, scen):
    fig, ax = plt.subplots(1, 1, figsize=(8, 5))
    col = f"dloss_{haz}_{scen}"

    tmp_box = df_country.dropna(subset=["income_group", col]).copy()
    sns.boxplot(
        data=tmp_box, x="income_group", y=col,
        ax=ax, color=haz_color_main[haz], width=0.6, fliersize=2
    )
    sns.stripplot(
        data=tmp_box, x="income_group", y=col,
        ax=ax, color="0.25", size=2.2, alpha=0.5, jitter=0.2
    )

    ax.set_xlabel("")
    ax.set_ylabel(f"Δ {haz_labels[haz].lower()}-induced yield loss (%)")
    ax.set_title(f"{haz_labels[haz]} risk increase by income group ({scen_labels[scen]})")
    ax.tick_params(axis="x", rotation=20)
    fig.savefig(out_fp, dpi=300, bbox_inches="tight")
    plt.close(fig)

def save_scatter(df_country, out_fp, haz, scen, q75):
    fig, ax = plt.subplots(1, 1, figsize=(8, 6))
    col = f"dloss_{haz}_{scen}"

    tmp_sc = df_country.dropna(subset=["imp_dep", col]).copy()

    palette = {
        "Low income": "#8c510a",
        "Lower middle income": "#d8b365",
        "Upper middle income": "#5ab4ac",
        "High income": "#01665e"
    }

    for grp in income_order:
        dsub = tmp_sc[tmp_sc["income_group"] == grp]
        if len(dsub) == 0:
            continue
        ax.scatter(
            dsub["imp_dep"],
            dsub[col],
            s=np.sqrt(dsub["imp_qty"].fillna(0) + 1) * 1.2,
            c=palette.get(grp, "0.5"),
            alpha=0.75,
            edgecolors="white",
            linewidths=0.4,
            label=grp
        )

    ax.axvline(THRESH_IMPORT_DEP, color="0.5", ls="--", lw=1)
    ax.axhline(q75, color="0.5", ls="--", lw=1)

    ax.set_xlabel("Import dependence ratio")
    ax.set_ylabel(f"Δ {haz_labels[haz].lower()}-induced yield loss (%)")
    ax.set_title(f"Import dependence versus {haz_labels[haz].lower()}-risk increase ({scen_labels[scen]})")
    ax.legend(frameon=False, fontsize=9, loc="best")
    fig.savefig(out_fp, dpi=300, bbox_inches="tight")
    plt.close(fig)

def save_barplot(df_country, out_fp, haz, scen):
    col = f"dloss_{haz}_{scen}"
    tmp_bar = (
        df_country.dropna(subset=["income_group", "trade_role", col])
        .groupby(["income_group", "trade_role"], observed=False)[col]
        .mean()
        .reset_index()
    )

    fig, ax = plt.subplots(1, 1, figsize=(9, 5))
    sns.barplot(
        data=tmp_bar,
        x="income_group",
        y=col,
        hue="trade_role",
        ax=ax,
        palette={
            "Net importer": "#ef8a62",
            "Net exporter": "#67a9cf",
            "Balanced": "#bdbdbd"
        }
    )

    ax.set_xlabel("")
    ax.set_ylabel(f"Mean Δ {haz_labels[haz].lower()}-induced yield loss (%)")
    ax.set_title(f"{haz_labels[haz]}-risk increase by trade role and income group ({scen_labels[scen]})")
    ax.legend(title="", frameon=False)
    ax.tick_params(axis="x", rotation=20)
    fig.savefig(out_fp, dpi=300, bbox_inches="tight")
    plt.close(fig)

def save_map(gdf_plot, out_fp, haz, scen):
    col = f"dloss_{haz}_{scen}"

    fig, ax = plt.subplots(1, 1, figsize=(12, 6))
    gdf_plot.plot(ax=ax, color="#d9d9d9", linewidth=0.25, edgecolor="white")

    gdf_est = gdf_plot[gdf_plot[col].notna()].copy()

    bounds = [0, 0.5, 1, 3, 5, 10]
    cmap = mcolors.ListedColormap(["#d9d9d9", "#fee5d9", "#fcae91", "#fb6a4a", "#cb181d"])
    norm = mcolors.BoundaryNorm(bounds, ncolors=cmap.N, clip=False)

    gdf_est.plot(
        column=col,
        ax=ax,
        cmap=cmap,
        norm=norm,
        linewidth=0.25,
        edgecolor="white"
    )

    ax.set_axis_off()
    ax.set_title(f"Country-scale increase in {haz_labels[haz].lower()}-induced yield loss ({scen_labels[scen]} − Historical)", fontsize=13)

    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    cbar = fig.colorbar(sm, ax=ax, orientation="horizontal", fraction=0.04, pad=0.03)
    cbar.set_label(f"Δ {haz_labels[haz].lower()}-induced yield loss (%)")
    fig.savefig(out_fp, dpi=300, bbox_inches="tight")
    plt.close(fig)

# =========================
# 4) 读取 loss 数据
# =========================
print("Reading loss nc ...")
ds = xr.open_dataset(nc_fp)
da = ds[varname]   # (scen, haz, lat, lon)
da = to_lon180(da, lon_name="lon")
if da.lat[0] > da.lat[-1]:
    da = da.sortby("lat")

# 未来 - 历史（两种情景，两种haz）
delta_dict = {}
for haz in hazards:
    delta_dict[haz] = da.sel(scen=fut_scens, haz=haz) - da.sel(scen="hist", haz=haz)

# =========================
# 5) 用 DB['dist'] 作为作物面积权重
# =========================
print("Using DB['dist'] as crop-area weights ...")
dist = DB['dist']
dist = to_lon180(dist, lon_name="lon")
if dist.lat[0] > dist.lat[-1]:
    dist = dist.sortby("lat")

area_total = (
    crop_area(dist, "maiz") +
    crop_area(dist, "soyb") +
    crop_area(dist, "rice") +
    crop_area(dist, "whea")
)

# 对齐到 loss 网格（取 dry 的网格作为基准）
ref_delta = delta_dict["dry"]
area_total = area_total.interp(lat=ref_delta.lat, lon=ref_delta.lon, method="nearest").fillna(0.0)

print("Area stats:")
print("  min :", float(area_total.min().values))
print("  max :", float(area_total.max().values))
print("  mean:", float(area_total.mean().values))

# =========================
# 6) 读取 shp（ISO_A3）
# =========================
print("Reading country shapefile ...")
gdf = gpd.read_file(shp_fp).to_crs(epsg=4326)

keep_cols = []
for c in ["ISO_A3", "NAME", "ADMIN", "geometry"]:
    if c in gdf.columns:
        keep_cols.append(c)
gdf = gdf[keep_cols].copy()

if "NAME" not in gdf.columns and "ADMIN" in gdf.columns:
    gdf["NAME"] = gdf["ADMIN"]
if "ADMIN" not in gdf.columns and "NAME" in gdf.columns:
    gdf["ADMIN"] = gdf["NAME"]

gdf = gdf[gdf["ISO_A3"].notna()].copy()
gdf = gdf[~gdf["ISO_A3"].isin(["-99"])].copy()

gdf = gdf.reset_index(drop=True).copy()
gdf["rid"] = np.arange(len(gdf))

df_base = gdf[["ISO_A3", "NAME"]].drop_duplicates().copy()
df_base = df_base.rename(columns={"NAME": "shp_name"})

# =========================
# 7) regionmask 映射格点到国家
# =========================
print("Building region mask ...")
regions = regionmask.from_geopandas(
    gdf,
    names="ISO_A3",
    numbers="rid"
)

mask = regions.mask(ref_delta.isel(scen=0))   # (lat, lon)
num_to_iso = dict(zip(gdf["rid"].values, gdf["ISO_A3"].values))

# =========================
# 8) 国家尺度面积加权聚合：两种haz × 两种scen
# =========================
print("Aggregating country-scale risk ...")
rows = []

for idx in gdf["rid"].values:
    iso3 = num_to_iso[idx]
    m = xr.where(mask == idx, 1.0, np.nan)

    crop_area_sum = float((area_total * m).sum(skipna=True).values)

    row = {
        "ISO_A3": iso3,
        "crop_area_sum": crop_area_sum,
        "has_crop": int(crop_area_sum > area_eps)
    }

    for haz in hazards:
        for scen in fut_scens:
            d = delta_dict[haz].sel(scen=scen)

            ww = area_total * m
            dd = d * m

            num = float((dd * ww).sum(skipna=True).values)
            den = float(ww.where(np.isfinite(dd)).sum(skipna=True).values)

            if (crop_area_sum <= area_eps) or (den <= area_eps):
                val = np.nan
            else:
                val = num / den

            row[f"dloss_{haz}_{scen}"] = val

    rows.append(row)

df_loss = pd.DataFrame(rows)

# 同时构造 long-format
long_rows = []
for _, r in df_loss.iterrows():
    for haz in hazards:
        for scen in fut_scens:
            long_rows.append({
                "ISO_A3": r["ISO_A3"],
                "crop_area_sum": r["crop_area_sum"],
                "has_crop": r["has_crop"],
                "haz": haz,
                "scen": scen,
                "dloss": r[f"dloss_{haz}_{scen}"]
            })
df_loss_long = pd.DataFrame(long_rows)

df_loss = df_loss.merge(
    gdf[["ISO_A3", "NAME"]].drop_duplicates(),
    on="ISO_A3",
    how="left"
)

# =========================
# 9) 读取 FAOSTAT 贸易数据（先按年份求四作物合计，再对年份求均值）
# =========================
print("Reading FAOSTAT trade data ...")
df_tr = pd.read_csv(faostat_fp)

df_tr = df_tr[["Area", "Element", "Item", "Year", "Value"]].copy()
df_tr["Year"] = pd.to_numeric(df_tr["Year"], errors="coerce")
df_tr["Value"] = pd.to_numeric(df_tr["Value"], errors="coerce")

df_tr = df_tr[(df_tr["Year"] >= YEAR_MIN) & (df_tr["Year"] <= YEAR_MAX)].copy()
df_tr = df_tr[df_tr["Item"].isin(crop_items.keys())].copy()

df_tr["ISO_A3"] = add_iso3_from_name(df_tr["Area"])
df_tr = df_tr[df_tr["ISO_A3"].notna()].copy()

tmp = (
    df_tr.groupby(["ISO_A3", "Area", "Element", "Year"], as_index=False)["Value"]
    .sum()
)

df_trade = (
    tmp.groupby(["ISO_A3", "Area", "Element"], as_index=False)["Value"]
    .mean()
)

df_trade_wide = df_trade.pivot_table(
    index=["ISO_A3", "Area"],
    columns="Element",
    values="Value",
    aggfunc="mean"
).reset_index()

df_trade_wide.columns.name = None

for c in ["Import quantity", "Export quantity"]:
    if c not in df_trade_wide.columns:
        df_trade_wide[c] = np.nan

df_trade_wide = df_trade_wide.rename(columns={
    "Area": "trade_name",
    "Import quantity": "imp_qty",
    "Export quantity": "exp_qty"
})

df_trade_wide["imp_qty"] = df_trade_wide["imp_qty"].fillna(0)
df_trade_wide["exp_qty"] = df_trade_wide["exp_qty"].fillna(0)
df_trade_wide["net_trade"] = df_trade_wide["exp_qty"] - df_trade_wide["imp_qty"]
df_trade_wide["imp_dep"] = df_trade_wide["imp_qty"] / (
    df_trade_wide["imp_qty"] + df_trade_wide["exp_qty"] + 1e-12
)

df_trade_wide["trade_role"] = df_trade_wide["net_trade"].apply(classify_trade_role)
df_trade_wide["imp_dep_cls"] = df_trade_wide["imp_dep"].apply(
    lambda x: classify_imp_dep(x, th=THRESH_IMPORT_DEP)
)

# =========================
# 10) 读取收入组
# =========================
print("Reading income groups ...")
df_inc = pd.read_excel(income_fp)

df_inc = df_inc.rename(columns={
    "Economy": "econ_name",
    "Code": "ISO_A3",
    "Income group": "income_group"
})

df_inc = df_inc[["ISO_A3", "econ_name", "income_group"]].copy()
df_inc = df_inc[df_inc["ISO_A3"].notna()].copy()

# =========================
# 11) 合并国家表（wide）
# =========================
print("Merging all tables ...")
df_country = (
    df_base
    .merge(df_loss.drop(columns=["NAME"], errors="ignore"), on="ISO_A3", how="left")
    .merge(df_trade_wide, on="ISO_A3", how="left")
    .merge(df_inc, on="ISO_A3", how="left")
)

df_country["country"] = (
    df_country["shp_name"]
    .fillna(df_country.get("econ_name"))
    .fillna(df_country.get("trade_name"))
)

df_country["crop_area_sum"] = df_country["crop_area_sum"].fillna(0)
df_country["has_crop"] = df_country["has_crop"].fillna(0).astype(int)

df_country["income_group"] = pd.Categorical(
    df_country["income_group"],
    categories=income_order,
    ordered=True
)

# 各 hazard × scen 计算高风险分类
q75_dict = {}
for haz in hazards:
    for scen in fut_scens:
        col = f"dloss_{haz}_{scen}"
        q75 = df_country[col].dropna().quantile(0.75)
        q75_dict[(haz, scen)] = q75

        out_col = f"risk_cls_{haz}_{scen}"
        df_country[out_col] = np.where(
            df_country[col].notna() & (df_country[col] >= q75),
            f"High {haz} risk",
            np.where(
                df_country[col].notna(),
                f"Lower {haz} risk",
                "No crop / no estimate"
            )
        )

# long-format 国家表
df_country_long = df_country[[
    "ISO_A3", "country", "crop_area_sum", "has_crop",
    "imp_qty", "exp_qty", "net_trade", "imp_dep", "trade_role", "imp_dep_cls",
    "income_group"
]].copy()

long_merge_rows = []
for _, r in df_country.iterrows():
    for haz in hazards:
        for scen in fut_scens:
            long_merge_rows.append({
                "ISO_A3": r["ISO_A3"],
                "country": r["country"],
                "crop_area_sum": r["crop_area_sum"],
                "has_crop": r["has_crop"],
                "imp_qty": r["imp_qty"],
                "exp_qty": r["exp_qty"],
                "net_trade": r["net_trade"],
                "imp_dep": r["imp_dep"],
                "trade_role": r["trade_role"],
                "imp_dep_cls": r["imp_dep_cls"],
                "income_group": r["income_group"],
                "haz": haz,
                "scen": scen,
                "dloss": r[f"dloss_{haz}_{scen}"],
                "risk_cls": r[f"risk_cls_{haz}_{scen}"]
            })
df_country_long = pd.DataFrame(long_merge_rows)

# =========================
# 12) 导出主表
# =========================
csv_main_wide = os.path.join(out_dir, "country_table_2scen_2haz_wide.csv")
csv_main_long = os.path.join(out_dir, "country_table_2scen_2haz_long.csv")

df_country.to_csv(csv_main_wide, index=False, encoding="utf-8-sig")
df_country_long.to_csv(csv_main_long, index=False, encoding="utf-8-sig")

print("Saved:", csv_main_wide)
print("Saved:", csv_main_long)

# =========================
# 13) 汇总表：收入组 / 贸易角色（long）
# =========================
df_income_sum = (
    df_country_long.groupby(["haz", "scen", "income_group"], observed=False)
    .agg(
        n_country_total=("ISO_A3", "count"),
        n_country_with_crop=("has_crop", "sum"),
        n_country_with_est=("dloss", lambda x: x.notna().sum()),
        mean_dloss=("dloss", "mean"),
        median_dloss=("dloss", "median"),
        mean_imp_dep=("imp_dep", "mean"),
        pct_net_importer=("trade_role", lambda x: np.mean(x == "Net importer") * 100),
        pct_high_imp_dep=("imp_dep_cls", lambda x: np.mean(x == "High import dependence") * 100),
    )
    .reset_index()
)
csv_sum1 = os.path.join(out_dir, "summary_by_income_group_2scen_2haz.csv")
df_income_sum.to_csv(csv_sum1, index=False, encoding="utf-8-sig")
print("Saved:", csv_sum1)

df_trade_sum = (
    df_country_long.groupby(["haz", "scen", "trade_role"], dropna=False)
    .agg(
        n_country_total=("ISO_A3", "count"),
        n_country_with_crop=("has_crop", "sum"),
        n_country_with_est=("dloss", lambda x: x.notna().sum()),
        mean_dloss=("dloss", "mean"),
        median_dloss=("dloss", "median"),
        mean_imp_dep=("imp_dep", "mean"),
    )
    .reset_index()
)
csv_sum2 = os.path.join(out_dir, "summary_by_trade_role_2scen_2haz.csv")
df_trade_sum.to_csv(csv_sum2, index=False, encoding="utf-8-sig")
print("Saved:", csv_sum2)

# =========================
# 14) 热点国家表（每个 hazard × scen 各导出一份）
# =========================
for haz in hazards:
    for scen in fut_scens:
        col = f"dloss_{haz}_{scen}"
        q75 = q75_dict[(haz, scen)]

        df_hot = df_country[
            (df_country["imp_dep"] >= THRESH_IMPORT_DEP) &
            (df_country[col] >= q75)
        ].copy()

        df_hot = df_hot.sort_values([col, "imp_dep"], ascending=[False, False])

        out_hot = os.path.join(out_dir, f"hotspot_high_impdep_highrisk_{haz}_{scen}.csv")
        df_hot[[
            "ISO_A3", "country", "income_group",
            "crop_area_sum", "has_crop",
            "imp_qty", "exp_qty", "net_trade", "imp_dep",
            "trade_role", col
        ]].to_csv(out_hot, index=False, encoding="utf-8-sig")
        print("Saved:", out_hot)

# =========================
# 15) 合并回 shp 准备绘图
# =========================
gdf_plot = gdf.merge(df_country, on="ISO_A3", how="left")

# =========================
# 16) 批量导出图
# =========================
for haz in hazards:
    for scen in fut_scens:
        q75 = q75_dict[(haz, scen)]

        out_map = os.path.join(out_dir, f"Map_{haz}_{scen}.pdf")
        out_box = os.path.join(out_dir, f"Box_income_{haz}_{scen}.pdf")
        out_sc  = os.path.join(out_dir, f"Scatter_impdep_{haz}_{scen}.pdf")
        out_bar = os.path.join(out_dir, f"Bar_trade_income_{haz}_{scen}.pdf")

        save_map(gdf_plot, out_map, haz, scen)
        save_boxplot(df_country, out_box, haz, scen)
        save_scatter(df_country, out_sc, haz, scen, q75)
        save_barplot(df_country, out_bar, haz, scen)

        print("Saved:", out_map)
        print("Saved:", out_box)
        print("Saved:", out_sc)
        print("Saved:", out_bar)

# =========================
# 17) 调试输出
# =========================
print("\n===== CHECKS =====")
print("Total countries in shp:", len(df_base))
print("Total countries in df_country:", len(df_country))
print("Countries with crop area:", int((df_country["has_crop"] == 1).sum()))

for haz in hazards:
    for scen in fut_scens:
        col = f"dloss_{haz}_{scen}"
        print(f"{haz}_{scen} countries with estimate:", int(df_country[col].notna().sum()))

print("\nDone.")

# %% cell 35
# -*- coding: utf-8 -*-
import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

# =========================
# 0) 全局绘图参数
# =========================
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 1) 基本配置
# =========================
crops = ["maiz", "soyb", "rice", "whea"]
modes = ["Irr", "nonIrr"]
hazards = ["dry", "wet"]
scenarios = ["hist", "ssp126", "ssp585"]

crop_labels = {"maiz": "Maize", "soyb": "Soybean", "rice": "Rice", "whea": "Wheat"}
haz_labels  = {"dry": "Drought", "wet": "Waterlogging"}
scen_labels = {"hist": "Historical", "ssp126": "SSP1-2.6", "ssp585": "SSP5-8.5"}

# 未来 GCM 列表（与你 DB 构建时一致）
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

# 输出目录
out_dir = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss"
os.makedirs(out_dir, exist_ok=True)

# 输出 nc
out_cf_nc        = os.path.join(out_dir, "prod_cf_4crop_3scen.nc")
out_total_loss_nc= os.path.join(out_dir, "prod_loss_total_abs_4crop_3scen.nc")
out_haz_abs_nc   = os.path.join(out_dir, "prod_loss_abs_4crop_3scen_2haz.nc")
out_haz_pct_nc   = os.path.join(out_dir, "prod_loss_pct_4crop_3scen_2haz.nc")

# 输出图（可选）
out_multicrop_pct_pdf = os.path.join(out_dir, "Fig_multicrop_prod_loss_pct_2x3.pdf")

# 是否把负损失裁成 0
clip_negative_loss = True

# =========================
# 2) 工具函数
# =========================
def ensure_lat_ascending(da):
    if "lat" in da.dims and da.lat[0] > da.lat[-1]:
        da = da.sortby("lat")
    return da

def interp_to_ref(da, ref, method="nearest"):
    """
    将 da 插值到 ref 的 lat/lon 网格
    """
    da = ensure_lat_ascending(da)
    ref = ensure_lat_ascending(ref)
    if ("lat" in da.dims) and ("lon" in da.dims):
        return da.interp(lat=ref.lat, lon=ref.lon, method=method)
    return da

def mean_over_year(da, n_last=30):
    """
    对时间维取后 n 年平均
    """
    if "year" not in da.dims:
        return da

    years = da.year.values

    # 排序保证时间正确
    years_sorted = np.sort(years)

    # 取最后 n 年
    last_years = years_sorted[-n_last:]

    da = da.sel(year=last_years)

    return da.mean("year", skipna=True)
def get_area_mode(dist_ds, crop, mode):
    """
    直接从 DB['dist'] 读取 crop-mode 面积
    注意：你现在的 dist = 分数 * cellarea，因此单位是 m²；
    若潜在产量单位是 t/ha，这里要换成 ha。
    """
    key = f"{crop}_{mode}"
    if key not in dist_ds:
        raise KeyError(f"{key} not found in DB['dist']")
    area_m2 = dist_ds[key].fillna(0.0)
    area_ha = area_m2 / 10000.0
    return area_ha

def get_block(DB, scen, gcm=None):
    """
    返回某个 scen 对应的数据块：
    hist -> DB['hist']
    ssp126 / ssp585 -> DB['future'][gcm][scen]
    """
    if scen == "hist":
        return DB["hist"]
    else:
        if gcm is None:
            raise ValueError("Future scenario requires gcm")
        return DB["future"][gcm][scen]

def safe_share(num, den):
    return xr.where(den > 0, num / den, 0)

def compute_one_block(block, dist, crops, modes, clip_negative_loss=True):
    """
    对一个 block（hist 或某个 gcm-scen）计算：
      prod_cf(crop,lat,lon)
      prod_loss_total_abs(crop,lat,lon)
      prod_loss_abs(crop,haz,lat,lon)
      prod_loss_pct(crop,haz,lat,lon)
    """
    crop_cf_list = []
    crop_total_loss_list = []
    crop_haz_abs_list = []
    crop_haz_pct_list = []

    # 统一参考网格：用 dist 的一个变量
    ref = dist["maiz_Irr"]

    for crop in crops:
        cf_crop = 0
        ws_crop = 0
        haz_abs_crop = {haz: 0 for haz in hazards}

        for mode in modes:
            # ========= 面积（ha）=========
            area = get_area_mode(dist, crop, mode)
            area = interp_to_ref(area, ref, method="nearest").fillna(0)

            # ========= 相对产量 =========
            ry_ws = block["relative_yield"][crop][f"waterstress_{mode}"]
            ry_cf = block["relative_yield"][crop][f"nowaterstress_{mode}"]

            ry_ws = mean_over_year(ry_ws)
            ry_cf = mean_over_year(ry_cf)

            ry_ws = interp_to_ref(ry_ws, ref, method="nearest")
            ry_cf = interp_to_ref(ry_cf, ref, method="nearest")

            # ========= 压力天数 =========
            dry_days = block["stress_days"][crop][f"dry_{mode}"]
            wet_days = block["stress_days"][crop][f"wet_{mode}"]

            dry_days = mean_over_year(dry_days)
            wet_days = mean_over_year(wet_days)

            dry_days = interp_to_ref(dry_days, ref, method="nearest").fillna(0)
            wet_days = interp_to_ref(wet_days, ref, method="nearest").fillna(0)

            # ========= 潜在产量 =========
            py = block["potential_yield"][crop][mode]
            py = mean_over_year(py)
            py = interp_to_ref(py, ref, method="nearest")

            # ========= 绝对产量 =========
            # 假设 py 单位为 t/ha，area 单位为 ha，则产量单位为 t
            prod_cf_mode = ry_cf * py * area
            prod_ws_mode = ry_ws * py * area

            # 总水分胁迫损失
            loss_total_mode = prod_cf_mode - prod_ws_mode
            if clip_negative_loss:
                loss_total_mode = xr.where(loss_total_mode < 0, 0, loss_total_mode)

            # ========= 按干/湿压力天数占比分配 =========
            stress_sum = dry_days + wet_days
            share_dry = safe_share(dry_days, stress_sum)
            share_wet = safe_share(wet_days, stress_sum)

            dry_loss_mode = loss_total_mode * share_dry
            wet_loss_mode = loss_total_mode * share_wet

            # ========= 累加到 crop 级 =========
            cf_crop = cf_crop + prod_cf_mode.fillna(0)
            ws_crop = ws_crop + prod_ws_mode.fillna(0)
            haz_abs_crop["dry"] = haz_abs_crop["dry"] + dry_loss_mode.fillna(0)
            haz_abs_crop["wet"] = haz_abs_crop["wet"] + wet_loss_mode.fillna(0)

        total_loss_crop = cf_crop - ws_crop
        if clip_negative_loss:
            total_loss_crop = xr.where(total_loss_crop < 0, 0, total_loss_crop)

        # 百分比：hazard-specific abs / counterfactual production
        haz_pct_crop = {}
        for haz in hazards:
            haz_pct_crop[haz] = xr.where(
                cf_crop > 0,
                haz_abs_crop[haz] / cf_crop * 100.0,
                np.nan
            )

        # 挂 crop 维
        cf_crop = cf_crop.expand_dims(crop=[crop])
        total_loss_crop = total_loss_crop.expand_dims(crop=[crop])

        haz_abs_da = xr.concat(
            [haz_abs_crop[haz].expand_dims(haz=[haz]) for haz in hazards],
            dim="haz"
        ).expand_dims(crop=[crop])

        haz_pct_da = xr.concat(
            [haz_pct_crop[haz].expand_dims(haz=[haz]) for haz in hazards],
            dim="haz"
        ).expand_dims(crop=[crop])

        crop_cf_list.append(cf_crop)
        crop_total_loss_list.append(total_loss_crop)
        crop_haz_abs_list.append(haz_abs_da)
        crop_haz_pct_list.append(haz_pct_da)

    prod_cf = xr.concat(crop_cf_list, dim="crop")
    prod_loss_total_abs = xr.concat(crop_total_loss_list, dim="crop")
    prod_loss_abs = xr.concat(crop_haz_abs_list, dim="crop")
    prod_loss_pct = xr.concat(crop_haz_pct_list, dim="crop")

    # 保证维度顺序
    prod_cf = prod_cf.transpose("crop", "lat", "lon")
    prod_loss_total_abs = prod_loss_total_abs.transpose("crop", "lat", "lon")
    prod_loss_abs = prod_loss_abs.transpose("crop", "haz", "lat", "lon")
    prod_loss_pct = prod_loss_pct.transpose("crop", "haz", "lat", "lon")

    return prod_cf, prod_loss_total_abs, prod_loss_abs, prod_loss_pct

# =========================
# 3) 读取 dist（直接来自 DB）
# =========================
dist = DB["dist"]
for v in dist.data_vars:
    dist[v] = ensure_lat_ascending(dist[v])

# =========================
# 4) 历史 + 未来 ensemble mean
# =========================
print("Computing historical ...")
hist_block = get_block(DB, "hist")
hist_cf, hist_total_abs, hist_haz_abs, hist_haz_pct = compute_one_block(
    hist_block, dist, crops, modes, clip_negative_loss=clip_negative_loss
)

# 加 scen 维
hist_cf        = hist_cf.expand_dims(scen=["hist"])
hist_total_abs = hist_total_abs.expand_dims(scen=["hist"])
hist_haz_abs   = hist_haz_abs.expand_dims(scen=["hist"])
hist_haz_pct   = hist_haz_pct.expand_dims(scen=["hist"])

# 未来：每个 scen 对 GCM 求 ensemble mean
fut_cf_list = []
fut_total_abs_list = []
fut_haz_abs_list = []
fut_haz_pct_list = []

for scen in ["ssp126", "ssp585"]:
    print(f"Computing future ensemble mean for {scen} ...")

    gcm_cf_list = []
    gcm_total_abs_list = []
    gcm_haz_abs_list = []
    gcm_haz_pct_list = []

    for gcm in gcms:
        print(f"  - {gcm}")
        blk = get_block(DB, scen, gcm=gcm)
        cf, total_abs, haz_abs, haz_pct = compute_one_block(
            blk, dist, crops, modes, clip_negative_loss=clip_negative_loss
        )

        gcm_cf_list.append(cf.expand_dims(gcm=[gcm]))
        gcm_total_abs_list.append(total_abs.expand_dims(gcm=[gcm]))
        gcm_haz_abs_list.append(haz_abs.expand_dims(gcm=[gcm]))
        gcm_haz_pct_list.append(haz_pct.expand_dims(gcm=[gcm]))

    cf_ens = xr.concat(gcm_cf_list, dim="gcm").mean("gcm", skipna=True).expand_dims(scen=[scen])
    total_abs_ens = xr.concat(gcm_total_abs_list, dim="gcm").mean("gcm", skipna=True).expand_dims(scen=[scen])
    haz_abs_ens = xr.concat(gcm_haz_abs_list, dim="gcm").mean("gcm", skipna=True).expand_dims(scen=[scen])
    haz_pct_ens = xr.concat(gcm_haz_pct_list, dim="gcm").mean("gcm", skipna=True).expand_dims(scen=[scen])

    fut_cf_list.append(cf_ens)
    fut_total_abs_list.append(total_abs_ens)
    fut_haz_abs_list.append(haz_abs_ens)
    fut_haz_pct_list.append(haz_pct_ens)

# 合并 3 个情景
prod_cf = xr.concat([hist_cf] + fut_cf_list, dim="scen").sel(scen=scenarios, crop=crops)
prod_loss_total_abs = xr.concat([hist_total_abs] + fut_total_abs_list, dim="scen").sel(scen=scenarios, crop=crops)
prod_loss_abs = xr.concat([hist_haz_abs] + fut_haz_abs_list, dim="scen").sel(scen=scenarios, crop=crops, haz=hazards)
prod_loss_pct = xr.concat([hist_haz_pct] + fut_haz_pct_list, dim="scen").sel(scen=scenarios, crop=crops, haz=hazards)

# =========================
# 5) 保存 nc
# =========================
prod_cf.to_dataset(name="prod_cf").to_netcdf(out_cf_nc)
prod_loss_total_abs.to_dataset(name="prod_loss_total_abs").to_netcdf(out_total_loss_nc)
prod_loss_abs.to_dataset(name="prod_loss_abs").to_netcdf(out_haz_abs_nc)
prod_loss_pct.to_dataset(name="prod_loss_pct").to_netcdf(out_haz_pct_nc)

print("Saved:", out_cf_nc)
print("Saved:", out_total_loss_nc)
print("Saved:", out_haz_abs_nc)
print("Saved:", out_haz_pct_nc)

# =========================
# 6) 可选：生成一个 multi-crop 百分比展示图
#    这里只用于快速检查，不作为国家分析输入
# =========================
# multi-crop absolute loss
haz_abs_multicrop = prod_loss_abs.sum("crop", skipna=True)
cf_multicrop = prod_cf.sum("crop", skipna=True)

haz_pct_multicrop = xr.where(cf_multicrop > 0, haz_abs_multicrop / cf_multicrop * 100.0, np.nan)

# 画 2 × 3：haz × scen
EXTENT = [-130, 160, -60, 75]
ASPECT = 1.3
cmap = plt.get_cmap("RdYlBu_r")
vmin, vmax = 0.0, 8.0
norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

fig = plt.figure(figsize=(12.5, 7.5))
gs = fig.add_gridspec(
    2, 3,
    left=0.03, right=0.88, top=0.93, bottom=0.06,
    wspace=0.02, hspace=0.08
)

mappable = None

for r, haz in enumerate(hazards):
    for c, scen in enumerate(scenarios):
        ax = fig.add_subplot(gs[r, c], projection=ccrs.PlateCarree())
        ax.set_extent(EXTENT, crs=ccrs.PlateCarree())
        ax.set_aspect(ASPECT)

        da = haz_pct_multicrop.sel(scen=scen, haz=haz)

        im = da.plot(
            ax=ax,
            transform=ccrs.PlateCarree(),
            cmap=cmap,
            norm=norm,
            add_colorbar=False
        )
        if mappable is None:
            mappable = im

        ax.coastlines(resolution="110m", linewidth=0.8)
        ax.add_feature(cfeature.LAND, facecolor="whitesmoke", zorder=0)

        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
        ax.set_frame_on(False)
        ax.patch.set_visible(False)
        if hasattr(ax, "spines"):
            for sp in ax.spines.values():
                sp.set_visible(False)

        if r == 0:
            ax.set_title(scen_labels[scen], fontsize=12)
        else:
            ax.set_title("")

        if c == 0:
            ax.text(
                -0.02, 0.5,
                haz_labels[haz] + "\n(Multi-crop production-loss %)",
                transform=ax.transAxes,
                ha="right", va="center",
                fontsize=11
            )

cax = fig.add_axes([0.90, 0.16, 0.018, 0.70])
cb = fig.colorbar(mappable, cax=cax, orientation="vertical", extend="neither")
cb.set_label("Production loss (%)")

fig.savefig(out_multicrop_pct_pdf, dpi=300, bbox_inches="tight")
plt.show()
print("Saved:", out_multicrop_pct_pdf)

# %% cell 36


# %% cell 37
fig, ax = plt.subplots(1, 1, figsize=(8, 6))

tmp_sc = df_country.dropna(subset=["imp_dep", "dloss_ssp585"]).copy()

palette = {
    "Low income": "#8c510a",
    "Lower middle income": "#d8b365",
    "Upper middle income": "#5ab4ac",
    "High income": "#01665e"
}

for grp in income_order:
    dsub = tmp_sc[tmp_sc["income_group"] == grp]
    if len(dsub) == 0:
        continue

    ax.scatter(
        dsub["imp_dep"],
        dsub["dloss_ssp585"],
        s=np.sqrt(dsub["imp_qty"].fillna(0) + 1) * 0.3,
        c=palette.get(grp, "0.5"),
        alpha=0.75,
        edgecolors="white",
        linewidths=0.4,
        label=grp
    )

ax.axvline(THRESH_IMPORT_DEP, color="0.5", ls="--", lw=1)
ax.axhline(q75, color="0.5", ls="--", lw=1)

ax.set_xlabel("Import dependence ratio")
ax.set_ylabel("Δ drought-induced yield loss (%)")
ax.set_title("Import dependence versus drought-risk increase (SSP5-8.5)")
ax.legend(frameon=False, fontsize=9, loc="best")

sc_pdf = os.path.join(out_dir, "Scatter_impdep_vs_dloss_ssp585.pdf")
fig.savefig(sc_pdf, dpi=300, bbox_inches="tight")
plt.close(fig)
print("Saved:", sc_pdf)

# %% cell 38
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import Patch

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 0) 输入文件（未来by_gcm + 历史）
# =========================
fut_by_gcm_nc = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\loss_future_2070_2100_by_gcm.nc"
hist_nc       = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\loss_hist_1980_2019.nc"

# 输出
out_csv = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\varDecomp_SSP_vs_GCM_4crops_2haz.csv"
out_pdf = r"D:\AAUDE\paper_v2\paper4\outputs_loss_maps\Fig_varDecomp_SSP_vs_GCM_4crops_2haz.pdf"

# =========================
# 1) 维度/标签
# =========================
crops = ["maiz", "soyb", "rice", "whea"]
crop_labels = {"maiz":"Maize", "soyb":"Soybean", "rice":"Rice", "whea":"Wheat"}

hazards = ["dry", "wet"]
haz_labels = {"dry":"Drought", "wet":"Waterlogging"}

scens = ["ssp126", "ssp585"]
scen_labels = {"ssp126":"SSP1-2.6", "ssp585":"SSP5-8.5"}

# =========================
# 2) 直接使用你已有的 dist = DB['dist']
#    如果你当前会话没有 DB，可以改成从文件读：
#    dist = xr.open_dataset(dist_nc_fp)
# =========================
dist = DB["dist"]  # <- 你说已经有了

def get_crop_area(dist_ds, crop):
    """返回该crop的总面积 (lat,lon)，Irr+nonIrr；whea优先 whea_*，否则回退 wwh+swh"""
    if crop in ["maiz", "soyb", "rice"]:
        return (dist_ds[f"{crop}_Irr"] + dist_ds[f"{crop}_nonIrr"]).fillna(0.0)
    elif crop == "whea":
        if ("whea_Irr" in dist_ds) and ("whea_nonIrr" in dist_ds):
            a = dist_ds["whea_Irr"] + dist_ds["whea_nonIrr"]
            if int(np.isfinite(a.values).sum()) == 0:
                a = (dist_ds["wwh_Irr"] + dist_ds["wwh_nonIrr"] +
                     dist_ds["swh_Irr"] + dist_ds["swh_nonIrr"])
        else:
            a = (dist_ds["wwh_Irr"] + dist_ds["wwh_nonIrr"] +
                 dist_ds["swh_Irr"] + dist_ds["swh_nonIrr"])
        return a.fillna(0.0)
    else:
        raise ValueError(crop)

# 组装 area(crop,lat,lon)
area = xr.concat(
    [get_crop_area(dist, c).rename("area").expand_dims(crop=[c]) for c in crops],
    dim="crop"
)

# =========================
# 3) 读取未来/历史 loss
# =========================
ds_fut = xr.open_dataset(fut_by_gcm_nc)
fut = ds_fut["loss_pct"].sel(scen=scens, crop=crops, haz=hazards)  # (gcm,scen,crop,haz,lat,lon)

ds_hist = xr.open_dataset(hist_nc)
hist = ds_hist["loss_pct"].sel(crop=crops, haz=hazards)           # (crop,haz,lat,lon)

# lat排序一致
if fut.lat[0] > fut.lat[-1]:
    fut = fut.sortby("lat")
if hist.lat[0] > hist.lat[-1]:
    hist = hist.sortby("lat")

# 对齐面积网格到 fut 网格
area = area.interp(lat=fut.lat, lon=fut.lon, method="nearest").fillna(0.0)

# =========================
# 4) 构造 Δloss = future - hist，并做“作物内面积加权全球平均”
# =========================
# hist 扩展到 (gcm,scen,crop,haz,lat,lon)
hist2 = hist.expand_dims(gcm=fut.gcm, scen=fut.scen)

delta = fut - hist2  # (gcm,scen,crop,haz,lat,lon)

# 作物内面积加权：每个 crop 用自己的 area[crop]
w = area  # (crop,lat,lon) 自动广播到 (gcm,scen,crop,haz,lat,lon)
wsum = w.sum(["lat", "lon"])
Y = (delta * w).sum(["lat", "lon"]) / xr.where(wsum > 0, wsum, np.nan)
# Y: (gcm,scen,crop,haz) —— 这是方差分解的输入

# =========================
# 5) SSP vs GCM 方差分解（两因子）
#    并把交互项按 SSP/GCM 相对大小分摊，确保两块加起来=1（便于堆叠图）
# =========================
def var_decomp_twofactor(Y2d, dim_g="gcm", dim_s="scen"):
    """
    Y2d: DataArray(dim_g, dim_s, ...)
    返回：var_ssp, var_gcm, var_int, var_total
    """
    mu = Y2d.mean([dim_g, dim_s], skipna=True)
    mu_g = Y2d.mean(dim_s, skipna=True)  # (gcm,...)
    mu_s = Y2d.mean(dim_g, skipna=True)  # (scen,...)

    var_g = (mu_g - mu).var(dim_g, skipna=True)
    var_s = (mu_s - mu).var(dim_s, skipna=True)

    inter = Y2d - mu_g - mu_s + mu
    var_i = inter.var([dim_g, dim_s], skipna=True)

    var_t = Y2d.var([dim_g, dim_s], skipna=True)
    return var_s, var_g, var_i, var_t

var_ssp, var_gcm, var_int, var_tot = var_decomp_twofactor(Y, dim_g="gcm", dim_s="scen")

eps = 1e-12
# 交互项分摊（按主效应比例）
w_ssp = var_ssp / (var_ssp + var_gcm + eps)
w_gcm = 1.0 - w_ssp
var_ssp_eff = var_ssp + w_ssp * var_int
var_gcm_eff = var_gcm + w_gcm * var_int

frac_ssp = (var_ssp_eff / (var_tot + eps)).clip(0, 1)  # (crop,haz)
frac_gcm = (var_gcm_eff / (var_tot + eps)).clip(0, 1)  # (crop,haz)

# =========================
# 6) 导出表格（方便你在正文里写“谁更大”）
# =========================
import pandas as pd

rows = []
for c in crops:
    for h in hazards:
        fs = float(frac_ssp.sel(crop=c, haz=h).values)
        fg = float(frac_gcm.sel(crop=c, haz=h).values)
        vt = float(var_tot.sel(crop=c, haz=h).values)
        dom = "SSP" if fs > fg else "GCM"
        rows.append({
            "crop": c,
            "haz": h,
            "frac_ssp": fs,
            "frac_gcm": fg,
            "var_total": vt,
            "dominant": dom
        })

df = pd.DataFrame(rows)
df.to_csv(out_csv, index=False)
print("Saved table:", out_csv)
print(df)

# =========================
# 7) 画图：每个作物两根柱（dry/wet），堆叠 SSP vs GCM
# =========================
fig, ax = plt.subplots(figsize=(6, 3))

col_ssp = "#fae5b8"  # 浅色：SSP
col_gcm = "#e79a90"  # 黄色：GCM

bar_w = 0.35
gap_crop = 0.55
gap_in = 0.12

x_positions = []
x_labels = []

x0 = 0.0
for crop in crops:
    x_dry = x0
    x_wet = x0 + bar_w + gap_in

    for h, x in zip(["dry", "wet"], [x_dry, x_wet]):
        fs = float(frac_ssp.sel(crop=crop, haz=h).values)
        fg = float(frac_gcm.sel(crop=crop, haz=h).values)

        # 保证fs+fg=1（数值误差）
        s = fs + fg
        if s > 0:
            fs, fg = fs/s, fg/s
        else:
            fs, fg = 0.0, 0.0

        ax.bar(x, fs, width=bar_w, color=col_ssp, edgecolor="none")
        ax.bar(x, fg, bottom=fs, width=bar_w, color=col_gcm, edgecolor="none")

        # 在柱子中间标注数值（可删）
        ax.text(x, fs/2, f"{fs:.2f}", ha="center", va="center", fontsize=9)
        ax.text(x, fs + fg/2, f"{fg:.2f}", ha="center", va="center", fontsize=9)

    # 作物标题
    ax.text((x_dry + x_wet)/2, 1.06, crop_labels.get(crop, crop),
            ha="center", va="bottom", fontsize=11)

    # x轴刻度标签：每个作物两个（dry/wet）
    x_positions.extend([x_dry, x_wet])
    x_labels.extend(["Drought", "Waterlogging"])

    x0 = x_wet + bar_w + gap_crop

ax.set_ylim(0, 1.0)
ax.set_yticks([0, 0.5, 1.0])
ax.set_ylabel("Fraction of total uncertainty (variance)")

ax.set_xticks(x_positions)
ax.set_xticklabels(x_labels)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

handles = [
    Patch(facecolor=col_gcm, edgecolor="none", label="Variance induced by GCMs"),
    Patch(facecolor=col_ssp, edgecolor="none", label="Variance induced by SSPs"),
]
ax.legend(handles=handles, loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.25))

plt.tight_layout()
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.show()
print("Saved figure:", out_pdf)

# %% [markdown]
# ### Country Trade depandence

# %% cell 40
# =========================
# 国家尺度产量损失聚合
# =========================
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import regionmask
import os

# =========================
# 1 读取刚生成的 nc
# =========================
prod_cf_nc = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\prod_cf_4crop_3scen.nc"
prod_loss_abs_nc = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\prod_loss_abs_4crop_3scen_2haz.nc"

prod_cf = xr.open_dataset(prod_cf_nc)["prod_cf"]
prod_loss_abs = xr.open_dataset(prod_loss_abs_nc)["prod_loss_abs"]

# 维度：
# prod_cf: (scen,crop,lat,lon)
# prod_loss_abs: (scen,crop,haz,lat,lon)

scenarios = ["hist","ssp126","ssp585"]
hazards = ["dry","wet"]

# =========================
# 2 读取国家 shp
# =========================
shp_fp = r"D:\AAUDE\paper_v2\paper2\data\ne_10m_admin_0_countries\ne_10m_admin_0_countries.shp"

gdf = gpd.read_file(shp_fp).to_crs(epsg=4326)

gdf = gdf[gdf["ISO_A3"].notna()]
gdf = gdf[~gdf["ISO_A3"].isin(["-99"])]

gdf = gdf.reset_index(drop=True)
gdf["rid"] = np.arange(len(gdf))

# =========================
# 3 构建 regionmask
# =========================
regions = regionmask.from_geopandas(
    gdf,
    numbers="rid",
    names="ISO_A3"
)

mask = regions.mask(prod_cf.isel(scen=0,crop=0))

num_to_iso = dict(zip(gdf["rid"], gdf["ISO_A3"]))

# =========================
# 4 国家聚合
# =========================
rows = []

for rid in gdf["rid"]:

    iso = num_to_iso[rid]

    m = xr.where(mask == rid, 1, np.nan)

    for scen in scenarios:

        # 国家反事实总产量
        cf_country = (
            prod_cf.sel(scen=scen)
            .sum("crop")
            * m
        ).sum(skipna=True)

        for haz in hazards:

            loss_country = (
                prod_loss_abs
                .sel(scen=scen, haz=haz)
                .sum("crop")
                * m
            ).sum(skipna=True)

            cf_val = float(cf_country.values)
            loss_val = float(loss_country.values)

            if cf_val > 0:
                pct = loss_val / cf_val * 100
            else:
                pct = np.nan

            rows.append({
                "ISO_A3": iso,
                "scenario": scen,
                "hazard": haz,
                "prod_cf": cf_val,
                "prod_loss": loss_val,
                "loss_pct": pct
            })

df_country = pd.DataFrame(rows)

# =========================
# 5 导出
# =========================
out_csv = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\country_prod_loss.csv"

df_country.to_csv(out_csv,index=False)

print("Saved:",out_csv)

# %% cell 41
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import country_converter as coco

# ===============================
# 文件路径
# ===============================
trade_fp   = r"D:\Edge_download\FAOSTAT_data_en_3-8-2026.csv"
prod_fp    = r"D:\Edge_download\FAOSTAT_data_en_3-8-2026 (1).csv"
income_fp  = r"D:\Edge_download\CLASS_2025_10_07.xlsx"
climate_fp = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\country_prod_loss.csv"

out_csv = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\country_food_risk_table.csv"

YEAR_MIN = 2014
YEAR_MAX = 2023

# 贸易文件作物名称
items_trade = [
    "Maize (corn)",
    "Rice",
    "Wheat",
    "Soybeans",
]

# 生产文件里大豆可能写成 Soya beans
items_prod = [
    "Maize (corn)",
    "Rice",
    "Wheat",
    "Soybeans",
    "Soya beans",
]

# 推荐阈值（FAO 口径下更合理）
THRESH_IMPORT_DEP = 0.4

# ===============================
# 1) 稳健的国家名 -> ISO_A3
# ===============================
def add_iso3_from_name(series):
    """
    将国家名稳健转换为 ISO_A3。
    输出只会是：
    - 单个字符串，如 'CHN'
    - 或 np.nan
    不会返回 list / ndarray
    """
    manual = {
        "Türkiye": "TUR",
        "Turkey": "TUR",
        "Viet Nam": "VNM",
        "Vietnam": "VNM",
        "Iran (Islamic Republic of)": "IRN",
        "Iran": "IRN",
        "Russian Federation": "RUS",
        "Russia": "RUS",
        "Bolivia (Plurinational State of)": "BOL",
        "Bolivia": "BOL",
        "Venezuela (Bolivarian Republic of)": "VEN",
        "Venezuela": "VEN",
        "Lao People's Democratic Republic": "LAO",
        "Laos": "LAO",
        "Syrian Arab Republic": "SYR",
        "Syria": "SYR",
        "Czechia": "CZE",
        "Czech Republic": "CZE",
        "Republic of Korea": "KOR",
        "South Korea": "KOR",
        "Democratic People's Republic of Korea": "PRK",
        "North Korea": "PRK",
        "United Republic of Tanzania": "TZA",
        "Tanzania": "TZA",
        "Côte d'Ivoire": "CIV",
        "Cote d'Ivoire": "CIV",
        "Micronesia (Federated States of)": "FSM",
        "Micronesia": "FSM",
        "Moldova, Republic of": "MDA",
        "Moldova": "MDA",
        "United States of America": "USA",
        "United States": "USA",
        "Bahamas, The": "BHS",
        "Gambia, The": "GMB",
        "China, mainland": "CHN",
        "China": "CHN",
        "Hong Kong SAR": "HKG",
        "China, Hong Kong SAR": "HKG",
        "China, Macao SAR": "MAC",
        "China, Taiwan Province of": "TWN",
        "Taiwan Province of China": "TWN",
        "Taiwan": "TWN",
    }

    out = []

    for nm in series.astype(str):
        nm = nm.strip()

        # 优先手工映射
        if nm in manual:
            out.append(manual[nm])
            continue

        try:
            x = coco.convert(names=nm, to="ISO3", not_found=None)

            if isinstance(x, (list, tuple, np.ndarray)):
                vals = []
                for v in x:
                    sv = str(v)
                    if v is not None and sv not in ["not found", "None", "nan"]:
                        vals.append(sv)
                out.append(vals[0] if len(vals) > 0 else np.nan)
            else:
                if x is None:
                    out.append(np.nan)
                else:
                    sx = str(x)
                    out.append(np.nan if sx in ["not found", "None", "nan"] else sx)

        except Exception:
            out.append(np.nan)

    return pd.Series(out, index=series.index)

# ===============================
# 2) 读取贸易数据
# ===============================
print("Reading trade data ...")
df_trade = pd.read_csv(trade_fp)

# 保留必要列
df_trade = df_trade[["Area", "Element", "Item", "Year", "Value"]].copy()

# 筛选四作物 + 年份
df_trade = df_trade[df_trade["Item"].isin(items_trade)].copy()
df_trade["Year"] = pd.to_numeric(df_trade["Year"], errors="coerce")
df_trade["Value"] = pd.to_numeric(df_trade["Value"], errors="coerce")

df_trade = df_trade[
    (df_trade["Year"] >= YEAR_MIN) &
    (df_trade["Year"] <= YEAR_MAX)
].copy()

# 国家名转 ISO_A3
df_trade["ISO_A3"] = add_iso3_from_name(df_trade["Area"])
df_trade = df_trade[df_trade["ISO_A3"].notna()].copy()
df_trade["ISO_A3"] = df_trade["ISO_A3"].astype(str)

# 去掉无效数值
df_trade = df_trade[df_trade["Value"].notna()].copy()

# 检查 ISO_A3 类型
print("Trade ISO_A3 types:")
print(df_trade["ISO_A3"].apply(type).value_counts())

# 每个国家-年份 四作物合计贸易量
trade_year = (
    df_trade
    .groupby(["ISO_A3", "Year", "Element"], as_index=False)["Value"]
    .sum()
)

# 多年平均
trade_mean = (
    trade_year
    .groupby(["ISO_A3", "Element"], as_index=False)["Value"]
    .mean()
)

trade_wide = trade_mean.pivot(
    index="ISO_A3",
    columns="Element",
    values="Value"
).reset_index()

trade_wide.columns.name = None

# 兼容列缺失
for c in ["Import quantity", "Export quantity"]:
    if c not in trade_wide.columns:
        trade_wide[c] = np.nan

trade_wide = trade_wide.rename(columns={
    "Import quantity": "import_qty",
    "Export quantity": "export_qty"
})

trade_wide["import_qty"] = trade_wide["import_qty"].fillna(0)
trade_wide["export_qty"] = trade_wide["export_qty"].fillna(0)
trade_wide["net_trade"] = trade_wide["export_qty"] - trade_wide["import_qty"]

# ===============================
# 3) 读取生产数据
# ===============================
print("Reading production data ...")
df_prod = pd.read_csv(prod_fp)

df_prod = df_prod[["Area", "Item", "Year", "Value"]].copy()
df_prod = df_prod[df_prod["Item"].isin(items_prod)].copy()

# 统一大豆名称
df_prod["Item"] = df_prod["Item"].replace({"Soya beans": "Soybeans"})

df_prod["Year"] = pd.to_numeric(df_prod["Year"], errors="coerce")
df_prod["Value"] = pd.to_numeric(df_prod["Value"], errors="coerce")

df_prod = df_prod[
    (df_prod["Year"] >= YEAR_MIN) &
    (df_prod["Year"] <= YEAR_MAX)
].copy()

df_prod["ISO_A3"] = add_iso3_from_name(df_prod["Area"])
df_prod = df_prod[df_prod["ISO_A3"].notna()].copy()
df_prod["ISO_A3"] = df_prod["ISO_A3"].astype(str)
df_prod = df_prod[df_prod["Value"].notna()].copy()

print("Production ISO_A3 types:")
print(df_prod["ISO_A3"].apply(type).value_counts())

# 每个国家-年份 四作物生产量合计
prod_year = (
    df_prod
    .groupby(["ISO_A3", "Year"], as_index=False)["Value"]
    .sum()
)

# 多年平均
prod_mean = (
    prod_year
    .groupby("ISO_A3", as_index=False)["Value"]
    .mean()
)

prod_mean = prod_mean.rename(columns={"Value": "production"})

# ===============================
# 4) 读取收入组
# ===============================
print("Reading income groups ...")
df_inc = pd.read_excel(income_fp)

df_inc = df_inc[["Code", "Income group", "Region", "Economy"]].copy()
df_inc = df_inc.rename(columns={
    "Code": "ISO_A3",
    "Income group": "income_group",
    "Region": "region",
    "Economy": "country_name_income"
})

df_inc["ISO_A3"] = df_inc["ISO_A3"].astype(str)

# ===============================
# 5) 读取气候损失
# ===============================
print("Reading climate country loss ...")
df_climate = pd.read_csv(climate_fp)

# 确保 ISO_A3 为字符串
df_climate["ISO_A3"] = df_climate["ISO_A3"].astype(str)

# 期望字段：
# ISO_A3, scenario, hazard, prod_cf, prod_loss, loss_pct

# ===============================
# 6) 合并
# ===============================
print("Merging all tables ...")
df = (
    df_climate
    .merge(prod_mean, on="ISO_A3", how="left")
    .merge(trade_wide, on="ISO_A3", how="left")
    .merge(df_inc, on="ISO_A3", how="left")
)

# ===============================
# 7) 衍生指标（按 FAO 公式）
# ===============================
# 贸易和生产缺失按 0 处理更合理
for col in ["production", "import_qty", "export_qty"]:
    df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

# 国内供给 = production + imports - exports
df["domestic_supply"] = (
    df["production"] + df["import_qty"] - df["export_qty"]
)

# 国内供给 <= 0 的情况视为无效
df["domestic_supply"] = df["domestic_supply"].where(df["domestic_supply"] > 0, np.nan)

# ===== FAO 标准 IDR =====
# IDR = Imports / (Production + Imports - Exports)
df["import_dependence"] = (
    df["import_qty"] /
    (df["domestic_supply"] + 1e-12)
)

# 百分数形式
df["import_dependence_pct"] = df["import_dependence"] * 100.0

# 国家损失 / 国家生产
df["loss_share_production"] = (
    df["prod_loss"] /
    (df["production"] + 1e-12)
)

# 贸易角色：三类
df["trade_role"] = np.select(
    [
        df["net_trade"] > 0,
        df["net_trade"] < 0,
        df["net_trade"] == 0
    ],
    [
        "Net exporter",
        "Net importer",
        "Balanced"
    ],
    default="Unknown"
)

# 高进口依赖阈值
df["high_import_dep"] = df["import_dependence"] >= THRESH_IMPORT_DEP

# 补国家名
if "country" not in df.columns:
    df["country"] = df["country_name_income"]
else:
    df["country"] = df["country"].fillna(df["country_name_income"])

# ===============================
# 8) 导出
# ===============================
df.to_csv(out_csv, index=False, encoding="utf-8-sig")
print("Saved:", out_csv)

# ===============================
# 9) 简单检查
# ===============================
print("\n===== CHECKS =====")
print("Rows:", len(df))
print("Unique countries:", df["ISO_A3"].nunique())
print("Scenarios:", df["scenario"].dropna().unique())
print("Hazards:", df["hazard"].dropna().unique())

print("\nMissing production:", df["production"].isna().sum())
print("Missing import_qty:", df["import_qty"].isna().sum())
print("Missing export_qty:", df["export_qty"].isna().sum())
print("Missing income_group:", df["income_group"].isna().sum())
print("Missing domestic_supply:", df["domestic_supply"].isna().sum())

print("\nHigh import dependence threshold:", THRESH_IMPORT_DEP)
print("High import dependence count:", df["high_import_dep"].sum())

print("\nSample:")
print(df.head())

# %% cell 42
df["trade_role"] = np.select(
    [
        df["net_trade"] > 0,
        df["net_trade"] < 0,
        df["net_trade"] == 0
    ],
    [
        "Net exporter",
        "Net importer",
        "Balanced"
    ],
    default="Unknown"
)

# 高进口依赖阈值
df["high_import_dep"] = df["import_dependence"] >= THRESH_IMPORT_DEP

# 补国家名
if "country" not in df.columns:
    df["country"] = df["country_name_income"]
else:
    df["country"] = df["country"].fillna(df["country_name_income"])

# ===============================
# 8) 导出
# ===============================
df.to_csv(out_csv, index=False, encoding="utf-8-sig")
print("Saved:", out_csv)

# ===============================
# 9) 简单检查
# ===============================
print("\n===== CHECKS =====")
print("Rows:", len(df))
print("Unique countries:", df["ISO_A3"].nunique())
print("Scenarios:", df["scenario"].dropna().unique())
print("Hazards:", df["hazard"].dropna().unique())

print("\nMissing production:", df["production"].isna().sum())
print("Missing import_qty:", df["import_qty"].isna().sum())
print("Missing export_qty:", df["export_qty"].isna().sum())
print("Missing income_group:", df["income_group"].isna().sum())
print("Missing domestic_supply:", df["domestic_supply"].isna().sum())

print("\nHigh import dependence threshold:", THRESH_IMPORT_DEP)
print("High import dependence count:", df["high_import_dep"].sum())

print("\nSample:")
print(df.head())

# %% cell 43
# -*- coding: utf-8 -*-
import os
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors

# =========================
# 0) 全局绘图参数
# =========================
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 1) 文件路径
# =========================
csv_fp = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\country_food_risk_table.csv"
shp_fp = r"D:\AAUDE\paper_v2\paper2\data\ne_10m_admin_0_countries\ne_10m_admin_0_countries.shp"

out_dir = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\figures"
os.makedirs(out_dir, exist_ok=True)

out_pdf = os.path.join(out_dir, "Fig_country_loss_pct_2x2.pdf")
out_png = os.path.join(out_dir, "Fig_country_loss_pct_2x2.png")

# =========================
# 2) 读取数据
# =========================
df = pd.read_csv(csv_fp)
gdf = gpd.read_file(shp_fp).to_crs(epsg=4326)

# 只保留有效国家
gdf = gdf[gdf["ISO_A3"].notna()].copy()
gdf = gdf[~gdf["ISO_A3"].isin(["-99"])].copy()

# 只保留未来两情景
df = df[df["scenario"].isin(["ssp126", "ssp585"])].copy()

# 合并
gdf_all = gdf.merge(df, on="ISO_A3", how="left")

# =========================
# 3) 配置
# =========================
scenarios = ["ssp126", "ssp585"]
hazards = ["dry", "wet"]

title_map = {
    ("ssp126", "dry"): "SSP1-2.6  Drought",
    ("ssp126", "wet"): "SSP1-2.6  Waterlogging",
    ("ssp585", "dry"): "SSP5-8.5  Drought",
    ("ssp585", "wet"): "SSP5-8.5  Waterlogging",
}

# 建议用高分位数做 vmax，避免极少数国家把色标拉太大
vmin = 0
vmax = df["loss_pct"].quantile(0.98)
if vmax <= 0:
    vmax = df["loss_pct"].max()

cmap = plt.get_cmap("Reds")
norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

# =========================
# 4) 绘图：2×2
# =========================
fig, axes = plt.subplots(2, 2, figsize=(13, 7.2))
axes = axes.flatten()

for i, (scen, haz) in enumerate([
    ("ssp126", "dry"),
    ("ssp126", "wet"),
    ("ssp585", "dry"),
    ("ssp585", "wet"),
]):
    ax = axes[i]

    sub = gdf_all[
        (gdf_all["scenario"] == scen) &
        (gdf_all["hazard"] == haz)
    ].copy()

    # 先画底色：无数据国家
    gdf.plot(
        ax=ax,
        color="#d9d9d9",
        linewidth=0.25,
        edgecolor="white"
    )

    # 再画有值国家
    sub_valid = sub[sub["loss_pct"].notna()].copy()
    if len(sub_valid) > 0:
        sub_valid.plot(
            column="loss_pct",
            ax=ax,
            cmap=cmap,
            norm=norm,
            linewidth=0.25,
            edgecolor="white"
        )
    ax.set_xlim(-180, 180)
    ax.set_ylim(-60, 90)
    ax.set_title(title_map[(scen, haz)], fontsize=12)
    ax.set_axis_off()

# 统一 colorbar
sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
sm._A = []

cbar = fig.colorbar(
    sm,
    ax=axes.tolist(),
    orientation="vertical",
    fraction=0.022,
    pad=0.02
)
cbar.set_label("Production loss (%)", fontsize=11)

fig.suptitle(
    "Country-level crop production losses under future climate scenarios",
    fontsize=14, y=0.98
)

plt.tight_layout(rect=[0, 0, 0.96, 0.96])

fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
fig.savefig(out_png, dpi=300, bbox_inches="tight")
plt.show()

print("Saved:", out_pdf)
print("Saved:", out_png)

# %% cell 44
# -*- coding: utf-8 -*-
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

csv_fp = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\country_food_risk_table.csv"
out_dir = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\figures"
os.makedirs(out_dir, exist_ok=True)

out_pdf = os.path.join(out_dir, "Fig_import_dependence_vs_loss_scatter.pdf")

df = pd.read_csv(csv_fp)

# 只画未来
df = df[df["scenario"].isin(["ssp126", "ssp585"])].copy()

# 去掉关键缺失
df = df.dropna(subset=["import_dependence", "loss_pct", "income_group"])

# 如果想用百分数形式，可改成 x="import_dependence_pct"
fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)

for ax, scen in zip(axes, ["ssp126", "ssp585"]):
    sub = df[df["scenario"] == scen].copy()

    sns.scatterplot(
        data=sub,
        x="import_dependence",
        y="loss_pct",
        hue="income_group",
        style="hazard",
        size="production",
        sizes=(20, 250),
        alpha=0.75,
        ax=ax
    )

    ax.axvline(0.4, color="0.5", linestyle="--", linewidth=1)
    ax.set_title("SSP1-2.6" if scen == "ssp126" else "SSP5-8.5")
    ax.set_xlabel("Import dependency ratio (FAO definition)")
    ax.set_ylabel("Production loss (%)")

axes[1].legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    borderaxespad=0,
    frameon=False
)

fig.suptitle("Import dependency versus climate-induced production loss", fontsize=14, y=0.98)
plt.tight_layout()
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.show()

print("Saved:", out_pdf)

# %% cell 45
# -*- coding: utf-8 -*-
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

csv_fp = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\country_food_risk_table.csv"
out_dir = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\figures"
os.makedirs(out_dir, exist_ok=True)

out_pdf = os.path.join(out_dir, "Fig_income_trade_role_bar.pdf")

df = pd.read_csv(csv_fp)
df = df[df["scenario"].isin(["ssp126", "ssp585"])].copy()
df = df.dropna(subset=["income_group", "trade_role", "loss_pct"])

income_order = [
    "Low income",
    "Lower middle income",
    "Upper middle income",
    "High income"
]
df["income_group"] = pd.Categorical(df["income_group"], categories=income_order, ordered=True)

fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharey=True)

plot_order = [
    ("ssp126", "dry"),
    ("ssp126", "wet"),
    ("ssp585", "dry"),
    ("ssp585", "wet"),
]

title_map = {
    ("ssp126", "dry"): "SSP1-2.6  Drought",
    ("ssp126", "wet"): "SSP1-2.6  Waterlogging",
    ("ssp585", "dry"): "SSP5-8.5  Drought",
    ("ssp585", "wet"): "SSP5-8.5  Waterlogging",
}

for ax, (scen, haz) in zip(axes.flatten(), plot_order):
    sub = df[(df["scenario"] == scen) & (df["hazard"] == haz)].copy()

    tmp = (
        sub.groupby(["income_group", "trade_role"], observed=False)["loss_pct"]
        .mean()
        .reset_index()
    )

    sns.barplot(
        data=tmp,
        x="income_group",
        y="loss_pct",
        hue="trade_role",
        ax=ax,
        palette={
            "Net exporter": "#67a9cf",
            "Net importer": "#ef8a62",
            "Balanced": "#bdbdbd"
        }
    )

    ax.set_title(title_map[(scen, haz)], fontsize=11)
    ax.set_xlabel("")
    ax.set_ylabel("Production loss (%)")
    ax.tick_params(axis="x", rotation=20)

# 只保留一个图例
handles, labels = axes[0, 0].get_legend_handles_labels()
for ax in axes.flatten():
    if ax.get_legend() is not None:
        ax.get_legend().remove()

fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.98))
fig.suptitle("Asymmetric losses across income groups and trade roles", fontsize=14, y=1.02)

plt.tight_layout()
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
plt.show()

print("Saved:", out_pdf)

# %% cell 46
sns.boxplot(data=df, x='income_group', y='loss_pct', hue='scenario')
plt.axhline(3, ls='--', color='red', label='High Loss Threshold')
plt.title('Loss Distribution by Income Group')
plt.savefig('box_income_loss.png')

# %% cell 47
bins = [0, 0.5, 1, 3, 5, np.inf]
df['loss_bin'] = pd.cut(df['loss_pct'], bins)
crosstab = pd.crosstab(df['income_group'], df['loss_bin'], normalize='index') * 100
crosstab.plot(kind='bar', stacked=True)
plt.title('Affected Area Proportion by Loss Threshold')
plt.savefig('stacked_area_loss.png')

# %% cell 48
# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D

# =========================
# 0) 全局绘图参数
# =========================
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 1) 文件路径
# =========================
csv_fp = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\country_food_risk_table.csv"
country_shp_fp = r"D:\AAUDE\paper_v2\paper2\data\ne_10m_admin_0_countries\ne_10m_admin_0_countries.shp"
ipcc_shp_fp = r"F:\世界区划\IPCC-WGI-reference-regions-v4_shapefile\IPCC-WGI-reference-regions-v4.shp"

out_dir = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\figures_circular_clean"
os.makedirs(out_dir, exist_ok=True)

out_pdf = os.path.join(out_dir, "Fig_circular_major_producers_ssp585_dry_wet_clean.pdf")
out_png = os.path.join(out_dir, "Fig_circular_major_producers_ssp585_dry_wet_clean.png")

# =========================
# 2) 读取数据
# =========================
df = pd.read_csv(csv_fp)
df = df[df["scenario"] == "ssp585"].copy()

meta = (
    df[["ISO_A3", "country", "production", "region"]]
    .drop_duplicates(subset=["ISO_A3"])
    .copy()
)

# =========================
# 3) 读取国家边界 + IPCC 区域
# =========================
gdf_country = gpd.read_file(country_shp_fp).to_crs(epsg=4326)
gdf_country = gdf_country[gdf_country["ISO_A3"].notna()].copy()
gdf_country = gdf_country[~gdf_country["ISO_A3"].isin(["-99"])].copy()
gdf_country = gdf_country[["ISO_A3", "geometry"]].copy()

gdf_ipcc = gpd.read_file(ipcc_shp_fp).to_crs(epsg=4326)
gdf_ipcc = gdf_ipcc[[c for c in ["geometry", "Acronym", "Name", "Continent"] if c in gdf_ipcc.columns]].copy()

gdf_country_cent = gdf_country.copy()
gdf_country_cent["geometry"] = gdf_country_cent.geometry.centroid

gdf_join = gpd.sjoin(
    gdf_country_cent,
    gdf_ipcc,
    how="left",
    predicate="within"
)

country_to_ipcc = (
    gdf_join[["ISO_A3", "Acronym", "Continent"]]
    .drop_duplicates(subset=["ISO_A3"])
    .rename(columns={"Acronym": "ipcc_region", "Continent": "ipcc_continent"})
)

meta = meta.merge(country_to_ipcc, on="ISO_A3", how="left")

# =========================
# 4) 把 IPCC 小区合并成大区域
# =========================
def collapse_region(row):
    cont = row.get("ipcc_continent", None)
    wb = row.get("region", None)

    if pd.notna(cont):
        cont = str(cont).upper()
        if cont == "NORTH-AMERICA":
            return "North America"
        elif cont == "CENTRAL-AMERICA":
            return "Central America & Caribbean"
        elif cont == "SOUTH-AMERICA":
            return "South America"
        elif cont == "EUROPE":
            return "Europe"
        elif cont == "ASIA":
            return "Asia"
        elif cont == "AFRICA":
            return "Africa"
        elif cont == "POLAR":
            return "Polar"
        elif cont == "OCEANIA":
            return "Oceania"

    # fallback 到原 region
    if pd.notna(wb):
        wb = str(wb)
        if "North America" in wb:
            return "North America"
        elif "Latin America" in wb or "Caribbean" in wb:
            return "Central America & Caribbean"
        elif "Europe" in wb:
            return "Europe"
        elif "Africa" in wb:
            return "Africa"
        elif "Asia" in wb or "Pacific" in wb:
            return "Asia"
        else:
            return wb

    return "Other"

meta["plot_region"] = meta.apply(collapse_region, axis=1)

# =========================
# 5) 国家筛选：改成 top 24
# =========================
TOP_N = 24
meta_top = meta.sort_values("production", ascending=False).head(TOP_N).copy()
top_iso = meta_top["ISO_A3"].tolist()

df_use = df[df["ISO_A3"].isin(top_iso)].copy()
df_use = df_use.merge(
    meta[["ISO_A3", "plot_region"]],
    on="ISO_A3",
    how="left"
)

# =========================
# 6) 区域颜色（大区域）
# =========================
region_color_map = {
    "North America": "#4C78A8",
    "Central America & Caribbean": "#A0CBE8",
    "South America": "#59A14F",
    "Europe": "#E15759",
    "Asia": "#9C9E3D",
    "Africa": "#B5651D",
    "Oceania": "#76B7B2",
    "Polar": "#BDBDBD",
    "Other": "#9E9E9E",
}
DEFAULT_REGION_COLOR = "#9E9E9E"

# =========================
# 7) 工具函数
# =========================
def prepare_plot_df(df_in, hazard, top_iso):
    sub = df_in[df_in["hazard"] == hazard].copy()
    sub = (
        sub[["ISO_A3", "country", "plot_region", "production", "loss_pct"]]
        .drop_duplicates(subset=["ISO_A3"])
        .copy()
    )
    sub = sub[sub["ISO_A3"].isin(top_iso)].copy()

    sub["loss_pct"] = pd.to_numeric(sub["loss_pct"], errors="coerce")
    sub["production"] = pd.to_numeric(sub["production"], errors="coerce")

    sub = sub[sub["loss_pct"].notna()].copy()
    sub = sub[sub["production"].notna()].copy()
    sub = sub[sub["plot_region"].notna()].copy()

    region_order = (
        sub.groupby("plot_region")["production"]
        .sum()
        .sort_values(ascending=False)
        .index.tolist()
    )

    sub["plot_region"] = pd.Categorical(sub["plot_region"], categories=region_order, ordered=True)
    sub = sub.sort_values(["plot_region", "loss_pct"], ascending=[True, True]).reset_index(drop=True)

    sub["color"] = sub["plot_region"].astype(str).map(region_color_map).fillna(DEFAULT_REGION_COLOR)
    return sub, region_order

def circular_bar_positions(df_plot, group_col="plot_region", gap=1.5, start_angle_deg=90):
    group_sizes = df_plot.groupby(group_col, observed=True).size().tolist()
    group_sizes = [int(x) for x in group_sizes if x > 0]

    n = len(df_plot)
    n_groups = len(group_sizes)
    total_slots = n + gap * n_groups
    width = 2 * np.pi / total_slots

    angles = []
    group_bounds = []
    current = 0

    for size in group_sizes:
        start_idx = current
        for _ in range(size):
            angles.append(start_idx * width + width / 2)
            start_idx += 1

        group_bounds.append((current * width, (current + size) * width))
        current += size + gap

    angles = np.array(angles)
    offset = np.deg2rad(start_angle_deg)
    angles = angles + offset
    group_bounds = [(a + offset, b + offset) for a, b in group_bounds]

    return angles, width, group_bounds

def text_rotation(angle_rad):
    angle_deg = np.degrees(angle_rad)
    rot = angle_deg - 90
    ha = "left"
    if 90 < angle_deg < 270:
        rot = angle_deg + 90
        ha = "right"
    return rot, ha

def draw_circular_panel(ax, df_plot, title, max_value=None, inner_radius=1.2, gap=1.5):
    df_plot = df_plot.copy()

    if len(df_plot) == 0:
        ax.set_axis_off()
        return

    if max_value is None:
        max_value = max(5.0, float(df_plot["loss_pct"].max()) * 1.08)

    angles, width, group_bounds = circular_bar_positions(
        df_plot, group_col="plot_region", gap=gap, start_angle_deg=90
    )

    heights = df_plot["loss_pct"].to_numpy(dtype=float)
    colors = df_plot["color"].astype(str).to_numpy()

    ax.set_theta_direction(-1)
    ax.set_theta_offset(np.pi / 2)
    ax.set_ylim(0, inner_radius + max_value + 1.0)

    ax.bar(
        angles,
        heights,
        width=width * 0.88,
        bottom=inner_radius,
        color=colors,
        edgecolor="white",
        linewidth=0.7,
        align="center"
    )

    # 参考圈
    tick_values = [2, 4, 6, 8]
    tick_values = [t for t in tick_values if t <= max_value]
    for t in tick_values:
        theta = np.linspace(0, 2*np.pi, 360)
        ax.plot(theta, np.full_like(theta, inner_radius + t), color="0.82", lw=0.8, zorder=0)

    for t in tick_values:
        ax.text(np.deg2rad(8), inner_radius + t, f"{t}%", color="0.55", fontsize=9, ha="left", va="center")

    # 国家标签
    label_r = inner_radius + heights + 0.35
    for ang, rlab, name in zip(angles, label_r, df_plot["ISO_A3"]):
        rot, ha = text_rotation(ang)
        ax.text(
            ang, rlab, name,
            rotation=rot, rotation_mode="anchor",
            ha=ha, va="center",
            fontsize=9, color="0.3"
        )

    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    ax.spines["polar"].set_visible(False)
    ax.set_title(title, fontsize=13, pad=22)

def make_region_legend(region_order):
    handles = []
    for reg in region_order:
        color = region_color_map.get(str(reg), DEFAULT_REGION_COLOR)
        handles.append(
            Line2D([0], [0], marker="o", color="none",
                   markerfacecolor=color, markeredgecolor=color,
                   markersize=11, label=str(reg))
        )
    return handles

# =========================
# 8) 构建数据
# =========================
dry_df, dry_region_order = prepare_plot_df(df_use, "dry", top_iso)
wet_df, wet_region_order = prepare_plot_df(df_use, "wet", top_iso)

region_order_all = []
for r in list(dry_region_order) + list(wet_region_order):
    if str(r) not in [str(x) for x in region_order_all]:
        region_order_all.append(r)

global_max = max(dry_df["loss_pct"].max(), wet_df["loss_pct"].max())
global_max = max(5.0, float(global_max) * 1.08)

# =========================
# 9) 绘图
# =========================
fig = plt.figure(figsize=(18, 9))
gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.15, 0.60], wspace=0.02)

ax1 = fig.add_subplot(gs[0, 0], projection="polar")
ax2 = fig.add_subplot(gs[0, 1], projection="polar")
ax_leg = fig.add_subplot(gs[0, 2])
ax_leg.axis("off")

draw_circular_panel(
    ax1, dry_df,
    title="Drought-induced production loss under SSP5-8.5 (%)",
    max_value=global_max,
    inner_radius=1.2,
    gap=1.5
)

draw_circular_panel(
    ax2, wet_df,
    title="Waterlogging-induced production loss under SSP5-8.5 (%)",
    max_value=global_max,
    inner_radius=1.2,
    gap=1.5
)

legend_handles = make_region_legend(region_order_all)
ax_leg.legend(
    handles=legend_handles,
    title="Regional grouping",
    loc="center left",
    frameon=False,
    fontsize=11,
    title_fontsize=12,
    handlelength=0.8,
    handletextpad=0.6,
    labelspacing=1.1
)

fig.suptitle(
    "Major food-producing countries under SSP5-8.5: drought versus waterlogging risk",
    fontsize=16, y=0.97
)

fig.savefig(out_pdf, dpi=300)
fig.savefig(out_png, dpi=300)
plt.show()

print("Saved:", out_pdf)
print("Saved:", out_png)

# %% cell 49
# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.lines import Line2D

# =========================
# 0) 全局绘图参数
# =========================
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "serif"
mpl.rcParams["font.serif"] = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 1) 文件路径
# =========================
csv_fp = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\country_food_risk_table.csv"

# 国家边界 shp（Natural Earth）
country_shp = r"D:\AAUDE\paper_v2\paper2\data\ne_10m_admin_0_countries\ne_10m_admin_0_countries.shp"

# 七大洲 shp
continent_shp = r"F:\世界区划\七大洲\continent.shp"

out_dir = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\figures_bar_major_region_sorted_countryname"
os.makedirs(out_dir, exist_ok=True)

out_pdf = os.path.join(out_dir, "Fig_major_producers_bar_ssp585_dry_wet_continent_sorted_countryname.pdf")
out_png = os.path.join(out_dir, "Fig_major_producers_bar_ssp585_dry_wet_continent_sorted_countryname.png")

# =========================
# 2) 读取损失数据
# =========================
df = pd.read_csv(csv_fp)

# 只保留 SSP585
df = df[df["scenario"] == "ssp585"].copy()

required_cols = [
    "ISO_A3", "country", "hazard", "loss_pct",
    "production"
]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"缺少必要字段: {missing}")

# 基本清洗
df["ISO_A3"] = df["ISO_A3"].astype(str).str.strip()
df["country"] = df["country"].astype(str).str.strip()
df["hazard"] = df["hazard"].astype(str).str.strip()

# =========================
# 3) 读取空间数据并建立 ISO_A3 -> CONTINENT 映射
# =========================
gdf_ctry = gpd.read_file(country_shp)
gdf_cont = gpd.read_file(continent_shp)

print("国家shp字段：", gdf_ctry.columns.tolist())
print("洲shp字段：", gdf_cont.columns.tolist())

if "ISO_A3" not in gdf_ctry.columns:
    raise ValueError("国家shp中未找到字段 'ISO_A3'，请检查字段名。")

if "CONTINENT" not in gdf_cont.columns:
    raise ValueError("continent.shp 中未找到字段 'CONTINENT'，请检查字段名。")

# 统一投影
if gdf_ctry.crs != gdf_cont.crs:
    gdf_ctry = gdf_ctry.to_crs(gdf_cont.crs)

# 去掉无效几何
gdf_ctry = gdf_ctry[gdf_ctry.geometry.notnull() & (~gdf_ctry.geometry.is_empty)].copy()
gdf_cont = gdf_cont[gdf_cont.geometry.notnull() & (~gdf_cont.geometry.is_empty)].copy()

# 用 representative_point 比 centroid 更稳，避免一些狭长/多岛国家质心落在国外
gdf_point = gdf_ctry[["ISO_A3", "geometry"]].copy()
gdf_point["ISO_A3"] = gdf_point["ISO_A3"].astype(str).str.strip()
gdf_point["geometry"] = gdf_point.geometry.representative_point()

# 空间连接：国家点落入洲 polygon
join = gpd.sjoin(
    gdf_point,
    gdf_cont[["CONTINENT", "geometry"]],
    how="left",
    predicate="within"
)

iso_continent = (
    join[["ISO_A3", "CONTINENT"]]
    .drop_duplicates(subset=["ISO_A3"])
    .rename(columns={"CONTINENT": "plot_region"})
)

# 若存在少量点未匹配，可打印检查
unmatched_iso = iso_continent[iso_continent["plot_region"].isna()]["ISO_A3"].tolist()
if len(unmatched_iso) > 0:
    print("以下 ISO_A3 未匹配到洲（将归为 Other）：", unmatched_iso[:20])

# =========================
# 4) 国家元数据：按 ISO_A3 保留一条
# =========================
# 若同一国家多行 production 一样，这里取最大值即可；若你的 production 本就唯一，不影响
meta = (
    df.groupby("ISO_A3", as_index=False)
      .agg({
          "country": "first",
          "production": "max"
      })
)

meta = meta.merge(iso_continent, on="ISO_A3", how="left")
meta["plot_region"] = meta["plot_region"].fillna("Other")

# =========================
# 5) 主产国：前20
# =========================
TOP_N = 20
meta_top = meta.sort_values("production", ascending=False).head(TOP_N).copy()
top_iso = meta_top["ISO_A3"].tolist()

df_use = df[df["ISO_A3"].isin(top_iso)].copy()
df_use = df_use.merge(
    meta[["ISO_A3", "country", "plot_region", "production"]],
    on="ISO_A3",
    how="left",
    suffixes=("", "_meta")
)

# 如果 merge 后 country 有重复列含义冲突，这里统一使用 meta 中的国家名
if "country_meta" in df_use.columns:
    df_use["country"] = df_use["country_meta"]
    df_use = df_use.drop(columns=["country_meta"])

if "production_meta" in df_use.columns:
    df_use["production"] = df_use["production_meta"]
    df_use = df_use.drop(columns=["production_meta"])

# =========================
# 6) 构造 dry / wet 宽表
# =========================
plot_df = (
    df_use.pivot_table(
        index=["ISO_A3", "country", "plot_region", "production"],
        columns="hazard",
        values="loss_pct",
        aggfunc="first"
    )
    .reset_index()
)

plot_df.columns.name = None

for c in ["dry", "wet"]:
    if c not in plot_df.columns:
        plot_df[c] = np.nan

plot_df = plot_df.rename(columns={
    "dry": "drought_loss",
    "wet": "waterlogging_loss"
})

# =========================
# 7) 排序：先洲，再洲内按 max loss
# =========================
region_order = [
    "North America",
    "South America",
    "Europe",
    "Asia",
    "Africa",
    "Oceania",
    "Antarctica",
    "Other"
]

plot_df["region_rank"] = plot_df["plot_region"].apply(
    lambda x: region_order.index(x) if x in region_order else 999
)

plot_df["risk_sort"] = plot_df[["drought_loss", "waterlogging_loss"]].max(axis=1)

plot_df = plot_df.sort_values(
    ["region_rank", "risk_sort"],
    ascending=[True, False]
).reset_index(drop=True)

# =========================
# 8) 颜色
# =========================
region_color_map = {
    "North America": "#4C78A8",
    "South America": "#59A14F",
    "Europe": "#E15759",
    "Asia": "#9C9E3D",
    "Africa": "#B5651D",
    "Oceania": "#76B7B2",
    "Antarctica": "#B0B0B0",
    "Other": "#9E9E9E",
}

plot_df["color"] = plot_df["plot_region"].map(region_color_map).fillna("#9E9E9E")

# =========================
# 9) 区域分隔位置
# =========================
region_breaks = []
region_labels = []

for reg in plot_df["plot_region"].drop_duplicates():
    idx = plot_df.index[plot_df["plot_region"] == reg].tolist()
    if len(idx) > 0:
        region_breaks.append((min(idx), max(idx)))
        region_labels.append(reg)

# =========================
# 10) y 轴标签
# =========================
yticklabels = plot_df["country"].tolist()

# 如需显示 ISO3 + 国家名，可改为：
# yticklabels = (plot_df["ISO_A3"] + "  " + plot_df["country"]).tolist()

# =========================
# 11) 作图
# =========================
fig, axes = plt.subplots(
    1, 2,
    figsize=(6, 12),
    sharey=True,
    gridspec_kw={"wspace": 0.08}
)

ax1, ax2 = axes
y = np.arange(len(plot_df))

# 左：Drought
ax1.barh(
    y,
    plot_df["drought_loss"],
    color=plot_df["color"],
    edgecolor="white",
    linewidth=0.7
)

# 右：Waterlogging
ax2.barh(
    y,
    plot_df["waterlogging_loss"],
    color=plot_df["color"],
    edgecolor="white",
    linewidth=0.7
)

# y轴国家名：只在左图显示
ax1.set_yticks(y)
ax1.set_yticklabels(yticklabels, fontsize=10)
ax2.tick_params(axis="y", which="both", left=False, labelleft=False)

# 标题
ax1.set_title("Drought", fontsize=12)
ax2.set_title("Waterlogging", fontsize=12)

# x轴标签
ax1.set_xlabel("Production loss (%)", fontsize=11)
ax2.set_xlabel("Production loss (%)", fontsize=11)

# 不要网格线；边框加粗；仅保留左/下边框
for ax in axes:
    ax.grid(False)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.spines["left"].set_linewidth(1.5)
    ax.spines["bottom"].set_linewidth(1.5)

    ax.tick_params(axis="both", labelsize=10, width=1.2)

# x轴范围
xmax1 = np.nanmax(plot_df["drought_loss"].values) if np.isfinite(np.nanmax(plot_df["drought_loss"].values)) else 5.0
xmax1 = max(5.0, float(xmax1) * 1.12)

xmax2 = np.nanmax(plot_df["waterlogging_loss"].values) if np.isfinite(np.nanmax(plot_df["waterlogging_loss"].values)) else 5.0
xmax2 = max(5.0, float(xmax2) * 1.12)

ax1.set_xlim(0, xmax1)
ax2.set_xlim(0, xmax2)

# 条形末端加数值
for i, v in enumerate(plot_df["drought_loss"]):
    if pd.notna(v) and v >= 0.8:
        ax1.text(
            v + xmax1 * 0.01, i, f"{v:.1f}",
            va="center", ha="left",
            fontsize=8, color="0.35"
        )

for i, v in enumerate(plot_df["waterlogging_loss"]):
    if pd.notna(v) and v >= 0.8:
        ax2.text(
            v + xmax2 * 0.01, i, f"{v:.1f}",
            va="center", ha="left",
            fontsize=8, color="0.35"
        )

# =========================
# 12) 区域分隔线 + 区域标签
# =========================
for start, end in region_breaks:
    if start > 0:
        ax1.axhline(start - 0.5, color="0.65", lw=1.0)
        ax2.axhline(start - 0.5, color="0.65", lw=1.0)

# 左侧加区域标签
for (start, end), reg in zip(region_breaks, region_labels):
    ymid = (start + end) / 2
    ax1.text(
        -xmax1 * 0.11, ymid, reg,
        ha="right", va="center",
        fontsize=10, color="0.3", fontweight="bold"
    )

# 反转 y 轴，让排序靠前的国家在上面
ax1.invert_yaxis()
ax2.invert_yaxis()

# =========================
# 13) 图例
# =========================
legend_order = [r for r in region_order if r in plot_df["plot_region"].unique().tolist()]
legend_handles = [
    Line2D(
        [0], [0],
        marker="s",
        color="none",
        markerfacecolor=region_color_map.get(r, "#9E9E9E"),
        markeredgecolor=region_color_map.get(r, "#9E9E9E"),
        markersize=9,
        label=r
    )
    for r in legend_order
]

fig.legend(
    handles=legend_handles,
    title="Continent",
    loc="lower center",
    ncol=min(4, len(legend_handles)),
    frameon=False,
    bbox_to_anchor=(0.5, -0.01),
    fontsize=10,
    title_fontsize=11
)

# 左边距加大，防止国家名和区域名被截断
plt.tight_layout(rect=[0.18, 0.06, 1, 0.95])

# 保存
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
fig.savefig(out_png, dpi=300, bbox_inches="tight")
plt.show()

print("Saved:", out_pdf)
print("Saved:", out_png)

# =========================
# 14) 可选：导出检查表
# =========================
check_csv = os.path.join(out_dir, "top20_country_continent_check.csv")
plot_df[[
    "ISO_A3", "country", "plot_region", "production",
    "drought_loss", "waterlogging_loss"
]].to_csv(check_csv, index=False, encoding="utf-8-sig")

print("Saved:", check_csv)

# %% cell 50


# %% cell 51


# %% cell 52
# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# =========================
# 0) 全局参数
# =========================
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 1) 路径
# =========================
fp = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\country_food_risk_table.csv"
out_dir = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\figures_final_income_prodsize_nature"
os.makedirs(out_dir, exist_ok=True)

# =========================
# 2) 读取数据
# =========================
df = pd.read_csv(fp)

required_cols = [
    "ISO_A3", "country", "scenario", "hazard",
    "loss_pct", "prod_loss", "production",
    "import_qty", "export_qty", "net_trade",
    "import_dependence", "income_group", "region", "trade_role"
]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"缺少必要字段: {missing}")

income_order = [
    "Low income",
    "Lower middle income",
    "Upper middle income",
    "High income"
]
df["income_group"] = pd.Categorical(df["income_group"], categories=income_order, ordered=True)

hazards = ["dry", "wet"]
future_scens = ["ssp126", "ssp585"]

haz_labels = {"dry": "Drought", "wet": "Waterlogging"}
scen_labels = {"ssp126": "SSP1-2.6", "ssp585": "SSP5-8.5"}

palette_income = {
    "Low income": "#8c510a",
    "Lower middle income": "#d8b365",
    "Upper middle income": "#5ab4ac",
    "High income": "#01665e"
}

THRESH_IMPORT_DEP = 0.4

# =========================
# 3) 坐标轴风格函数
# =========================
def nature_axis_style(ax):
    """
    只保留左侧和底部边框，不显示上和右；无网格；轴线加粗
    """
    ax.grid(False)

    ax.spines["left"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)

    ax.tick_params(
        axis="both",
        which="both",
        direction="out",
        width=1.1,
        length=4,
        labelsize=9
    )

# =========================
# 4) 构造宽表
# =========================
def build_wide_for_hazard(df_in, haz):
    sub = df_in[df_in["hazard"] == haz].copy()

    id_cols = ["ISO_A3", "country", "income_group", "region"]
    meta_cols = [
        "ISO_A3", "country", "income_group", "region",
        "production", "import_qty", "export_qty",
        "net_trade", "import_dependence", "trade_role"
    ]
    meta_cols = [c for c in meta_cols if c in sub.columns]
    meta = sub[meta_cols].drop_duplicates(subset=["ISO_A3"]).copy()

    hist = (
        sub[sub["scenario"] == "hist"][id_cols + ["loss_pct", "prod_loss"]]
        .rename(columns={"loss_pct": "loss_pct_hist", "prod_loss": "prod_loss_hist"})
        .drop_duplicates(subset=["ISO_A3"])
    )

    ssp126 = (
        sub[sub["scenario"] == "ssp126"][id_cols + ["loss_pct", "prod_loss"]]
        .rename(columns={"loss_pct": "loss_pct_ssp126", "prod_loss": "prod_loss_ssp126"})
        .drop_duplicates(subset=["ISO_A3"])
    )

    ssp585 = (
        sub[sub["scenario"] == "ssp585"][id_cols + ["loss_pct", "prod_loss"]]
        .rename(columns={"loss_pct": "loss_pct_ssp585", "prod_loss": "prod_loss_ssp585"})
        .drop_duplicates(subset=["ISO_A3"])
    )

    wide = (
        hist.merge(
            ssp126[["ISO_A3", "loss_pct_ssp126", "prod_loss_ssp126"]],
            on="ISO_A3", how="outer"
        )
        .merge(
            ssp585[["ISO_A3", "loss_pct_ssp585", "prod_loss_ssp585"]],
            on="ISO_A3", how="outer"
        )
    )

    wide = wide.merge(
    meta.drop(columns=["country", "income_group", "region", "ipcc_region"], errors="ignore"),
    on="ISO_A3", how="left"
)

    for scen in ["ssp126", "ssp585"]:
        wide[f"delta_loss_pct_{scen}"] = wide[f"loss_pct_{scen}"] - wide["loss_pct_hist"]
        wide[f"delta_prod_loss_{scen}"] = wide[f"prod_loss_{scen}"] - wide["prod_loss_hist"]

    eps = 0.05
    for scen in ["ssp126", "ssp585"]:
        wide[f"relative_increase_{scen}"] = np.where(
            wide["loss_pct_hist"].abs() >= eps,
            (wide[f"loss_pct_{scen}"] - wide["loss_pct_hist"]) / wide["loss_pct_hist"] * 100.0,
            np.nan
        )

    total_prod = wide["production"].sum(skipna=True)
    wide["production_share"] = wide["production"] / (total_prod + 1e-12)

    wide["import_dependence_plot"] = pd.to_numeric(
        wide["import_dependence"], errors="coerce"
    ).clip(lower=0, upper=1)

    wide["high_import_dep"] = wide["import_dependence_plot"] >= THRESH_IMPORT_DEP
    wide["hazard"] = haz
    return wide

# =========================
# 5) 工具函数
# =========================
def bubble_size(series, ref_max, scale=700, min_size=18, max_size=1400):
    """
    气泡大小 ~ sqrt(production / ref_max)
    """
    s = pd.to_numeric(series, errors="coerce").fillna(0).astype(float)
    if ref_max is None or ref_max <= 0:
        ref_max = max(float(s.max()), 1.0)

    out = np.sqrt(np.clip(s, 0, None) / ref_max) * scale
    out = np.clip(out, min_size, max_size)
    return out

def choose_labels_hist_future(sub, scen, n_hist=4, n_future=4, n_delta=4):
    idx = set()
    idx.update(sub.nlargest(n_hist, "loss_pct_hist").index.tolist())
    idx.update(sub.nlargest(n_future, f"loss_pct_{scen}").index.tolist())
    idx.update(sub.nlargest(n_delta, f"delta_loss_pct_{scen}").index.tolist())
    return list(idx)

def choose_labels_quadrant(sub, scen, x_thr=0.0, y_thr=0.4, n_high=10):
    risk = sub[
        (sub[f"delta_loss_pct_{scen}"] >= x_thr) &
        (sub["import_dependence_plot"] >= y_thr)
    ].copy()

    if len(risk) == 0:
        risk = sub.copy()

    risk["score"] = (
        risk[f"delta_loss_pct_{scen}"].rank(ascending=False, pct=True) +
        risk["import_dependence_plot"].rank(ascending=False, pct=True)
    )
    return risk.sort_values("score", ascending=False).head(n_high).index.tolist()

def add_labels(ax, data, idx_list, xcol, ycol, textcol="ISO_A3", fontsize=8):
    """
    自动调整标签位置，尽量避免跑出边界。
    """
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()

    xr = x1 - x0
    yr = y1 - y0

    for _, r in data.loc[idx_list].iterrows():
        x = r[xcol]
        y = r[ycol]
        label = str(r[textcol])

        if pd.isna(x) or pd.isna(y):
            continue

        dx = 0.012 * xr
        dy = 0.012 * yr
        ha = "left"
        va = "bottom"

        # 靠右 -> 往左放
        if x > x0 + 0.82 * xr:
            dx = -0.012 * xr
            ha = "right"

        # 靠上 -> 往下放
        if y > y0 + 0.85 * yr:
            dy = -0.012 * yr
            va = "top"

        tx = x + dx
        ty = y + dy

        # 防止文字跑出坐标轴边界
        tx = min(max(tx, x0 + 0.01 * xr), x1 - 0.01 * xr)
        ty = min(max(ty, y0 + 0.01 * yr), y1 - 0.01 * yr)

        ax.text(
            tx, ty, label,
            fontsize=fontsize,
            ha=ha, va=va,
            color="0.2"
        )

def add_income_legend(ax):
    handles = []
    for grp in income_order:
        handles.append(
            plt.Line2D(
                [0], [0],
                marker="o", linestyle="none",
                markerfacecolor=palette_income.get(grp, "0.5"),
                markeredgecolor="white",
                markeredgewidth=0.5,
                markersize=8,
                label=grp
            )
        )
    leg1 = ax.legend(
        handles=handles,
        title="Income group",
        frameon=False,
        fontsize=9,
        title_fontsize=10,
        loc="upper right"
    )
    ax.add_artist(leg1)

def add_size_legend(ax, ref_max, size_values, scale=700, title="Production"):
    size_handles = []
    for v in size_values:
        size_handles.append(
            plt.scatter(
                [], [], s=bubble_size(pd.Series([v]), ref_max=ref_max, scale=scale)[0],
                color="gray", alpha=0.35, edgecolors="white", linewidths=0.5,
                label=f"{v/1e6:.0f} Mt"
            )
        )
    ax.legend(
        handles=size_handles,
        title=title,
        scatterpoints=1,
        frameon=False,
        fontsize=8.5,
        title_fontsize=10,
        loc="lower right"
    )

# =========================
# 6) 图1：历史 vs 未来
# =========================
def plot_hist_vs_future_clean(wide, haz, scen, out_dir, ref_prod_max):
    sub = wide.dropna(subset=["loss_pct_hist", f"loss_pct_{scen}", "income_group", "production"]).copy()

    prod_thr = sub["production"].quantile(0.10)
    sub = sub[sub["production"] >= prod_thr].copy()

    fig, ax = plt.subplots(figsize=(7.6, 6.0))

    for grp in income_order:
        dsub = sub[sub["income_group"] == grp]
        if len(dsub) == 0:
            continue

        ax.scatter(
            dsub["loss_pct_hist"],
            dsub[f"loss_pct_{scen}"],
            s=bubble_size(dsub["production"], ref_max=ref_prod_max, scale=700),
            c=palette_income.get(grp, "0.5"),
            alpha=0.78,
            edgecolors="white",
            linewidths=0.5
        )

    allv = pd.concat([sub["loss_pct_hist"], sub[f"loss_pct_{scen}"]], axis=0)
    vmax = float(allv.quantile(0.98))
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = float(allv.max())

    ax.plot([0, vmax], [0, vmax], ls="--", lw=1, color="0.45")
    ax.text(vmax * 0.68, vmax * 0.78, "1:1 line", fontsize=9, color="0.35")

    idx_list = choose_labels_hist_future(sub, scen, n_hist=4, n_future=4, n_delta=4)
    add_labels(ax, sub, idx_list, "loss_pct_hist", f"loss_pct_{scen}", fontsize=8)

    ax.set_xlabel("Historical production loss (%)", fontsize=10)
    ax.set_ylabel("Future production loss (%)", fontsize=10)
    ax.set_title(f"{haz_labels[haz]}: historical baseline versus future loss ({scen_labels[scen]})", fontsize=11)

    nature_axis_style(ax)
    add_income_legend(ax)
    add_size_legend(ax, ref_max=ref_prod_max, size_values=[1e7, 5e7, 1e8], scale=700, title="Production")

    plt.tight_layout()

    pdf = os.path.join(out_dir, f"Fig_hist_vs_future_clean_{haz}_{scen}.pdf")
    png = os.path.join(out_dir, f"Fig_hist_vs_future_clean_{haz}_{scen}.png")
    fig.savefig(pdf, dpi=300, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(fig)

# =========================
# 7) 图2：四象限图
# =========================
def plot_quadrant_risk_clean(wide, haz, scen, out_dir, import_thr, ref_prod_max):
    sub = wide.dropna(
        subset=[f"delta_loss_pct_{scen}", "import_dependence_plot", "income_group", "production"]
    ).copy()

    prod_thr = sub["production"].quantile(0.15)
    sub = sub[sub["production"] >= prod_thr].copy()

    fig, ax = plt.subplots(figsize=(7.9, 6.0))

    for grp in income_order:
        dsub = sub[sub["income_group"] == grp]
        if len(dsub) == 0:
            continue

        ax.scatter(
            dsub[f"delta_loss_pct_{scen}"],
            dsub["import_dependence_plot"],
            s=bubble_size(dsub["production"], ref_max=ref_prod_max, scale=650),
            c=palette_income.get(grp, "0.5"),
            alpha=0.75,
            edgecolors="white",
            linewidths=0.45
        )

    x_thr = 0.0
    y_thr = import_thr

    ax.axvline(x_thr, color="0.5", ls="--", lw=1)
    ax.axhline(y_thr, color="0.5", ls="--", lw=1)

    x_left = float(sub[f"delta_loss_pct_{scen}"].quantile(0.02))
    x_right = float(sub[f"delta_loss_pct_{scen}"].quantile(0.97))
    x_left = min(x_left, -0.5)
    x_right = max(x_right, 1.0)

    ax.set_xlim(x_left, x_right * 1.03)
    ax.set_ylim(-0.1, 1.1)

    ax.text(x_thr + 0.02 * (x_right - x_left), 0.96, "Risk increase", fontsize=9, color="0.35")
    ax.text(x_right * 0.72, y_thr + 0.02, "High import dependence", fontsize=9, color="0.35")

    idx_list = choose_labels_quadrant(sub, scen, x_thr=x_thr, y_thr=y_thr, n_high=10)
    add_labels(ax, sub, idx_list, f"delta_loss_pct_{scen}", "import_dependence_plot", fontsize=7.5)

    ax.set_xlabel("Future loss increase (percentage points)", fontsize=10)
    ax.set_ylabel("Import dependency ratio", fontsize=10)
    ax.set_title(f"{haz_labels[haz]}: four-quadrant country risk ({scen_labels[scen]})", fontsize=11)

    nature_axis_style(ax)
    add_income_legend(ax)
    add_size_legend(ax, ref_max=ref_prod_max, size_values=[1e7, 5e7, 1e8], scale=650, title="Production")

    plt.tight_layout()

    pdf = os.path.join(out_dir, f"Fig_quadrant_clean_{haz}_{scen}.pdf")
    png = os.path.join(out_dir, f"Fig_quadrant_clean_{haz}_{scen}.png")
    fig.savefig(pdf, dpi=300, bbox_inches="tight")
    fig.savefig(png, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(fig)
# =========================
# 9) 批量生成
# =========================
global_prod_max = pd.to_numeric(df["production"], errors="coerce").max()
if not np.isfinite(global_prod_max) or global_prod_max <= 0:
    global_prod_max = 1.0

for haz in hazards:
    print(f"\n===== Building optimized figures for {haz} =====")
    wide = build_wide_for_hazard(df, haz)

    for scen in future_scens:
        print(f"Plotting {haz} - {scen}")

        plot_hist_vs_future_clean(
            wide, haz, scen, out_dir,
            ref_prod_max=global_prod_max
        )

        plot_quadrant_risk_clean(
            wide, haz, scen, out_dir,
            import_thr=THRESH_IMPORT_DEP,
            ref_prod_max=global_prod_max
        )

print("\nAll optimized figures finished.")

# %% cell 53
# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd

# =========================
# 1) 路径
# =========================
fp = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\country_food_risk_table.csv"
out_dir = r"D:\AAUDE\paper_v2\paper4\outputs_country_prod_loss\quadrant_stats_ssp585_dry"
os.makedirs(out_dir, exist_ok=True)

# =========================
# 2) 参数
# =========================
haz = "wet"
scen = "ssp585"
THRESH_IMPORT_DEP = 0.4
X_THR = 0.0

income_order = [
    "Low income",
    "Lower middle income",
    "Upper middle income",
    "High income"
]

# =========================
# 3) 读取数据
# =========================
df = pd.read_csv(fp)

required_cols = [
    "ISO_A3", "country", "scenario", "hazard",
    "loss_pct", "prod_loss", "production",
    "import_qty", "export_qty", "net_trade",
    "import_dependence", "income_group", "region", "trade_role"
]
missing = [c for c in required_cols if c not in df.columns]
if missing:
    raise ValueError(f"缺少必要字段: {missing}")

df["income_group"] = pd.Categorical(
    df["income_group"],
    categories=income_order,
    ordered=True
)

# =========================
# 4) 构造宽表
# =========================
def build_wide_for_hazard(df_in, haz):
    sub = df_in[df_in["hazard"] == haz].copy()

    id_cols = ["ISO_A3", "country", "income_group", "region"]
    meta_cols = [
        "ISO_A3", "country", "income_group", "region",
        "production", "import_qty", "export_qty",
        "net_trade", "import_dependence", "trade_role"
    ]
    meta_cols = [c for c in meta_cols if c in sub.columns]

    meta = sub[meta_cols].drop_duplicates(subset=["ISO_A3"]).copy()

    hist = (
        sub[sub["scenario"] == "hist"][id_cols + ["loss_pct", "prod_loss"]]
        .rename(columns={
            "loss_pct": "loss_pct_hist",
            "prod_loss": "prod_loss_hist"
        })
        .drop_duplicates(subset=["ISO_A3"])
    )

    ssp126 = (
        sub[sub["scenario"] == "ssp126"][id_cols + ["loss_pct", "prod_loss"]]
        .rename(columns={
            "loss_pct": "loss_pct_ssp126",
            "prod_loss": "prod_loss_ssp126"
        })
        .drop_duplicates(subset=["ISO_A3"])
    )

    ssp585 = (
        sub[sub["scenario"] == "ssp585"][id_cols + ["loss_pct", "prod_loss"]]
        .rename(columns={
            "loss_pct": "loss_pct_ssp585",
            "prod_loss": "prod_loss_ssp585"
        })
        .drop_duplicates(subset=["ISO_A3"])
    )

    wide = (
        hist.merge(
            ssp126[["ISO_A3", "loss_pct_ssp126", "prod_loss_ssp126"]],
            on="ISO_A3", how="outer"
        )
        .merge(
            ssp585[["ISO_A3", "loss_pct_ssp585", "prod_loss_ssp585"]],
            on="ISO_A3", how="outer"
        )
    )

    wide = wide.merge(
        meta.drop(columns=["country", "income_group", "region", "ipcc_region"], errors="ignore"),
        on="ISO_A3", how="left"
    )

    for scen_i in ["ssp126", "ssp585"]:
        wide[f"delta_loss_pct_{scen_i}"] = wide[f"loss_pct_{scen_i}"] - wide["loss_pct_hist"]
        wide[f"delta_prod_loss_{scen_i}"] = wide[f"prod_loss_{scen_i}"] - wide["prod_loss_hist"]

    eps = 0.05
    for scen_i in ["ssp126", "ssp585"]:
        wide[f"relative_increase_{scen_i}"] = np.where(
            wide["loss_pct_hist"].abs() >= eps,
            (wide[f"loss_pct_{scen_i}"] - wide["loss_pct_hist"]) / wide["loss_pct_hist"] * 100.0,
            np.nan
        )

    total_prod = pd.to_numeric(wide["production"], errors="coerce").sum(skipna=True)
    wide["production_share"] = pd.to_numeric(wide["production"], errors="coerce") / (total_prod + 1e-12)

    wide["import_dependence_plot"] = pd.to_numeric(
        wide["import_dependence"], errors="coerce"
    ).clip(lower=0, upper=1)

    wide["high_import_dep"] = wide["import_dependence_plot"] >= THRESH_IMPORT_DEP
    wide["hazard"] = haz

    return wide

# =========================
# 5) 四象限统计函数
# =========================
def quadrant_statistics(wide, scen, y_thr=0.4, x_thr=0.0):
    dfq = wide.copy()

    # 转数值
    dfq[f"delta_loss_pct_{scen}"] = pd.to_numeric(dfq[f"delta_loss_pct_{scen}"], errors="coerce")
    dfq["import_dependence_plot"] = pd.to_numeric(dfq["import_dependence_plot"], errors="coerce")
    dfq["production"] = pd.to_numeric(dfq["production"], errors="coerce")

    # 只保留有效值
    dfq = dfq.dropna(subset=[
        f"delta_loss_pct_{scen}",
        "import_dependence_plot",
        "production"
    ]).copy()

    # 四象限定义
    x = dfq[f"delta_loss_pct_{scen}"]
    y = dfq["import_dependence_plot"]

    dfq["quadrant"] = np.select(
        [
            (x >= x_thr) & (y >= y_thr),   # 右上
            (x >= x_thr) & (y <  y_thr),   # 右下
            (x <  x_thr) & (y >= y_thr),   # 左上
            (x <  x_thr) & (y <  y_thr)    # 左下
        ],
        [
            "Q1_right_upper",
            "Q2_right_lower",
            "Q3_left_upper",
            "Q4_left_lower"
        ],
        default="Unknown"
    )

    total_prod = dfq["production"].sum()

    stats = (
        dfq.groupby("quadrant", as_index=False)
        .agg(
            country_count=("ISO_A3", "nunique"),
            production=("production", "sum")
        )
    )
    stats["production_share_pct"] = stats["production"] / (total_prod + 1e-12) * 100

    # 重点象限
    q1 = dfq[dfq["quadrant"] == "Q1_right_upper"].copy()
    q2 = dfq[dfq["quadrant"] == "Q2_right_lower"].copy()

    # 收入组分布
    income_q1 = (
        q1.groupby("income_group", dropna=False)
        .agg(
            country_count=("ISO_A3", "nunique"),
            production=("production", "sum")
        )
        .reset_index()
    )
    income_q1["production_share_within_q1_pct"] = (
        income_q1["production"] / (q1["production"].sum() + 1e-12) * 100
    )

    income_q2 = (
        q2.groupby("income_group", dropna=False)
        .agg(
            country_count=("ISO_A3", "nunique"),
            production=("production", "sum")
        )
        .reset_index()
    )
    income_q2["production_share_within_q2_pct"] = (
        income_q2["production"] / (q2["production"].sum() + 1e-12) * 100
    )

    # 区域分布
    region_q1 = (
        q1.groupby("region", dropna=False)
        .agg(
            country_count=("ISO_A3", "nunique"),
            production=("production", "sum")
        )
        .reset_index()
        .sort_values("production", ascending=False)
    )
    region_q1["production_share_within_q1_pct"] = (
        region_q1["production"] / (q1["production"].sum() + 1e-12) * 100
    )

    region_q2 = (
        q2.groupby("region", dropna=False)
        .agg(
            country_count=("ISO_A3", "nunique"),
            production=("production", "sum")
        )
        .reset_index()
        .sort_values("production", ascending=False)
    )
    region_q2["production_share_within_q2_pct"] = (
        region_q2["production"] / (q2["production"].sum() + 1e-12) * 100
    )

    return stats, dfq, q1, q2, income_q1, income_q2, region_q1, region_q2

# =========================
# 6) 运行分析
# =========================
wide = build_wide_for_hazard(df, haz=haz)

stats, df_quad, q1, q2, income_q1, income_q2, region_q1, region_q2 = quadrant_statistics(
    wide=wide,
    scen=scen,
    y_thr=THRESH_IMPORT_DEP,
    x_thr=X_THR
)

# =========================
# 7) 汇总核心结果
# =========================
total_country = df_quad["ISO_A3"].nunique()
total_prod = df_quad["production"].sum()

q1_country_n = q1["ISO_A3"].nunique()
q2_country_n = q2["ISO_A3"].nunique()

q1_country_pct = q1_country_n / (total_country + 1e-12) * 100
q2_country_pct = q2_country_n / (total_country + 1e-12) * 100

q1_prod = q1["production"].sum()
q2_prod = q2["production"].sum()

q1_prod_pct = q1_prod / (total_prod + 1e-12) * 100
q2_prod_pct = q2_prod / (total_prod + 1e-12) * 100

summary = pd.DataFrame({
    "metric": [
        "total_countries",
        "total_production",
        "Q1_right_upper_country_count",
        "Q1_right_upper_country_share_pct",
        "Q1_right_upper_production",
        "Q1_right_upper_production_share_pct",
        "Q2_right_lower_country_count",
        "Q2_right_lower_country_share_pct",
        "Q2_right_lower_production",
        "Q2_right_lower_production_share_pct"
    ],
    "value": [
        total_country,
        total_prod,
        q1_country_n,
        q1_country_pct,
        q1_prod,
        q1_prod_pct,
        q2_country_n,
        q2_country_pct,
        q2_prod,
        q2_prod_pct
    ]
})

# =========================
# 8) 打印结果
# =========================
print("\n================= 总体四象限统计 =================")
print(stats)

print("\n================= 重点汇总 =================")
print(summary)

print("\n================= 右上象限（Q1）主要国家 =================")
print(
    q1[[
        "ISO_A3", "country", "income_group", "region",
        "production", "import_dependence_plot", f"delta_loss_pct_{scen}",
        "loss_pct_hist", f"loss_pct_{scen}", "trade_role"
    ]]
    .sort_values("production", ascending=False)
    .head(30)
)

print("\n================= 右下象限（Q2）主要国家 =================")
print(
    q2[[
        "ISO_A3", "country", "income_group", "region",
        "production", "import_dependence_plot", f"delta_loss_pct_{scen}",
        "loss_pct_hist", f"loss_pct_{scen}", "trade_role"
    ]]
    .sort_values("production", ascending=False)
    .head(30)
)

print("\n================= 右上象限收入组分布 =================")
print(income_q1)

print("\n================= 右下象限收入组分布 =================")
print(income_q2)

print("\n================= 右上象限区域分布 =================")
print(region_q1.head(20))

print("\n================= 右下象限区域分布 =================")
print(region_q2.head(20))

# =========================
# 9) 导出结果表
# =========================
stats.to_csv(os.path.join(out_dir, "quadrant_overall_stats.csv"), index=False, encoding="utf-8-sig")
summary.to_csv(os.path.join(out_dir, "quadrant_key_summary.csv"), index=False, encoding="utf-8-sig")

q1.sort_values("production", ascending=False).to_csv(
    os.path.join(out_dir, "Q1_right_upper_countries.csv"),
    index=False, encoding="utf-8-sig"
)
q2.sort_values("production", ascending=False).to_csv(
    os.path.join(out_dir, "Q2_right_lower_countries.csv"),
    index=False, encoding="utf-8-sig"
)

income_q1.to_csv(os.path.join(out_dir, "Q1_income_distribution.csv"), index=False, encoding="utf-8-sig")
income_q2.to_csv(os.path.join(out_dir, "Q2_income_distribution.csv"), index=False, encoding="utf-8-sig")
region_q1.to_csv(os.path.join(out_dir, "Q1_region_distribution.csv"), index=False, encoding="utf-8-sig")
region_q2.to_csv(os.path.join(out_dir, "Q2_region_distribution.csv"), index=False, encoding="utf-8-sig")

# 四象限全体国家表
df_quad.sort_values(["quadrant", "production"], ascending=[True, False]).to_csv(
    os.path.join(out_dir, "all_countries_with_quadrants.csv"),
    index=False, encoding="utf-8-sig"
)

print("\n结果已保存到：")
print(out_dir)

# %% [markdown]
# ### Validation

# %% cell 55
import importlib
import validate_admin_custom_criteria as vac
vac = importlib.reload(vac)

out = vac.run_validation_admin(
    DB,
    year_min=1900,
    year_max=2019,
    dry_event_days_thr=7,
    wet_event_days_thr=7,
    yield_unit_factor=1.0,
    relative_yield_scale=1.0,
)

# %% cell 56
import importlib
import validate_admin_custom_criteria as vac
vac = importlib.reload(vac)

out = vac.run_validation_admin(
    DB,
    year_min=1900,
    year_max=2019,
    dry_event_days_thr=7,
    wet_event_days_thr=7,
    yield_unit_factor=1.0,
    relative_yield_scale=1.0,
)

print(out)

# %% cell 57
# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import xarray as xr
import geopandas as gpd
import pycountry
import matplotlib.pyplot as plt
import matplotlib as mpl

# =========================
# 0) 全局设置
# =========================
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 1) 路径
# =========================
shp_fp = r"D:\AAUDE\paper_v2\paper2\data\ne_10m_admin_0_countries\ne_10m_admin_0_countries.shp"
fao_fp = r"D:\Edge_download\FAOSTAT_data_en_3-9-2026.csv"
out_dir = r"D:\AAUDE\paper_v2\paper4\validation_fao_yield"
os.makedirs(out_dir, exist_ok=True)

# =========================
# 2) 参数：这里一定要按你的模型真实含义改
# =========================
# 如果 relative_yield = (4 - gs_loss)，且其真实范围应归一到 0-1，则设为 4.0
# 如果它本身已经是 0-1 的相对产量，则设为 1.0
RELATIVE_YIELD_SCALE = 4.0

# 如果 potential_yield 单位是 t/ha，则乘 1000 变成 kg/ha
# 如果本身就是 kg/ha，则设为 1.0
MODEL_YIELD_TO_KGHA = 1000.0

YEAR_MIN = 1960
YEAR_MAX = 2019

CROPS = ["maiz", "soyb", "rice", "whea"]

# FAO Item 名称映射（先按最常见写法；跑前建议 print(unique items) 检查一下）
FAO_ITEM_MAP = {
    "maiz": ["Maize (corn)"],
    "soyb": ["Soya beans", "Soybeans"],
    "rice": ["Rice", "Rice, paddy"],
    "whea": ["Wheat"]
}

# =========================
# 3) 一些辅助函数
# =========================
def get_model_actual_yield_da(DB, crop,
                              rel_scale=1.0,
                              unit_factor=1.0):
    """
    从 DB['hist'] 中构造历史期 crop 的实际单产 (kg/ha)
    返回:
        yield_da(year, lat, lon) : 国家比较用的格点实际单产
        area_da(lat, lon)        : MIRCA 收获面积权重（灌溉+雨养）
    """
    hist = DB["hist"]

    # ---- 取 relative yield（实际情景 waterstress）
    rel_rf  = hist["relative_yield"][crop]["waterstress_nonIrr"] / rel_scale
    rel_irr = hist["relative_yield"][crop]["waterstress_Irr"]    / rel_scale

    # ---- 取 potential yield
    pot_rf  = hist["potential_yield"][crop]["nonIrr"] * unit_factor
    pot_irr = hist["potential_yield"][crop]["Irr"]    * unit_factor

    # ---- 格点实际单产
    y_rf  = rel_rf * pot_rf
    y_irr = rel_irr * pot_irr

    # ---- MIRCA 面积权重
    area_rf  = DB["dist"][f"{crop}_nonIrr"]
    area_irr = DB["dist"][f"{crop}_Irr"]

    # ---- 合成国家比较所需的总单产（面积加权）
    area_tot = area_rf.fillna(0) + area_irr.fillna(0)

    yield_combined = xr.where(
        area_tot > 0,
        (y_rf.fillna(0) * area_rf.fillna(0) + y_irr.fillna(0) * area_irr.fillna(0)) / area_tot,
        np.nan
    )

    return yield_combined, area_tot.where(area_tot > 0)


def build_grid_country_lookup(template_da, shp_fp, iso_field="ISO_A3"):
    """
    给 template_da 的每个格点中心分配国家 ISO_A3
    返回:
        iso_idx_2d: (lat, lon) 的整数编码数组，0表示无国家
        idx_to_iso: dict[int] -> ISO_A3
    """
    world = gpd.read_file(shp_fp)[[iso_field, "geometry"]].copy()
    world = world[world[iso_field].notna()].to_crs("EPSG:4326")

    lats = template_da["lat"].values
    lons = template_da["lon"].values
    lon2d, lat2d = np.meshgrid(lons, lats)

    pts = gpd.GeoDataFrame(
        {
            "row": np.repeat(np.arange(len(lats)), len(lons)),
            "col": np.tile(np.arange(len(lons)), len(lats)),
            "lon": lon2d.ravel(),
            "lat": lat2d.ravel(),
        },
        geometry=gpd.points_from_xy(lon2d.ravel(), lat2d.ravel()),
        crs="EPSG:4326"
    )

    joined = gpd.sjoin(
        pts,
        world[[iso_field, "geometry"]],
        how="left",
        predicate="within"
    )

    # 对于落在海岸线边界上的点，within 可能漏掉，再用 intersects 补一次
    miss = joined[iso_field].isna()
    if miss.any():
        joined2 = gpd.sjoin(
            pts.loc[miss],
            world[[iso_field, "geometry"]],
            how="left",
            predicate="intersects"
        )
        joined.loc[miss, iso_field] = joined2[iso_field].values

    iso_series = joined[iso_field].fillna("")

    iso_unique = sorted([x for x in iso_series.unique() if x != ""])
    iso_to_idx = {iso: i + 1 for i, iso in enumerate(iso_unique)}
    idx_to_iso = {i + 1: iso for i, iso in enumerate(iso_unique)}

    iso_idx = np.array([iso_to_idx.get(x, 0) for x in iso_series], dtype=np.int32)
    iso_idx_2d = iso_idx.reshape(len(lats), len(lons))

    return iso_idx_2d, idx_to_iso


def aggregate_country_yield(yield_da, area_da, iso_idx_2d, idx_to_iso):
    """
    将格点单产聚合到国家-年份
    yield_da: (year, lat, lon) kg/ha
    area_da : (lat, lon) harvested area weights
    """
    years = yield_da["year"].values.astype(int)

    area2d = area_da.values
    n_country = max(idx_to_iso.keys()) if len(idx_to_iso) > 0 else 0

    out = []

    for yr in years:
        y2d = yield_da.sel(year=yr).values

        m = np.isfinite(y2d) & np.isfinite(area2d) & (area2d > 0) & (iso_idx_2d > 0)
        if not np.any(m):
            continue

        iso_idx = iso_idx_2d[m]
        w_area  = area2d[m]
        prod    = y2d[m] * w_area

        sum_prod = np.bincount(iso_idx, weights=prod, minlength=n_country + 1)
        sum_area = np.bincount(iso_idx, weights=w_area, minlength=n_country + 1)

        cy = np.divide(
            sum_prod, sum_area,
            out=np.full_like(sum_prod, np.nan, dtype=float),
            where=sum_area > 0
        )

        for i in range(1, n_country + 1):
            if np.isfinite(cy[i]):
                out.append((idx_to_iso[i], int(yr), float(cy[i])))

    df = pd.DataFrame(out, columns=["ISO_A3", "Year", "model_yield_kg_ha"])
    return df


def name_to_iso3(name):
    """
    将 FAO Area 名称转为 ISO3
    """
    manual = {
        "Bolivia (Plurinational State of)": "BOL",
        "Cabo Verde": "CPV",
        "China, mainland": "CHN",
        "China, Taiwan Province of": "TWN",
        "China, Hong Kong SAR": "HKG",
        "China, Macao SAR": "MAC",
        "Côte d'Ivoire": "CIV",
        "Democratic People's Republic of Korea": "PRK",
        "Democratic Republic of the Congo": "COD",
        "Iran (Islamic Republic of)": "IRN",
        "Lao People's Democratic Republic": "LAO",
        "Micronesia (Federated States of)": "FSM",
        "Republic of Korea": "KOR",
        "Republic of Moldova": "MDA",
        "Russian Federation": "RUS",
        "Syrian Arab Republic": "SYR",
        "Türkiye": "TUR",
        "United Kingdom of Great Britain and Northern Ireland": "GBR",
        "United Republic of Tanzania": "TZA",
        "United States of America": "USA",
        "Venezuela (Bolivarian Republic of)": "VEN",
        "Viet Nam": "VNM",
        "Brunei Darussalam": "BRN",
        "Eswatini": "SWZ",
        "Palestine": "PSE",
    }
    if name in manual:
        return manual[name]

    try:
        return pycountry.countries.lookup(name).alpha_3
    except Exception:
        return np.nan


def load_fao_yield(fao_fp, crop_item_map, year_min=1960, year_max=2019):
    """
    读取 FAO 单产数据，统一成:
    ISO_A3, Year, crop, fao_yield_kg_ha
    """
    fao = pd.read_csv(fao_fp)

    # 兼容不同列名
    colmap = {}
    for c in fao.columns:
        c2 = c.strip()
        colmap[c] = c2
    fao = fao.rename(columns=colmap)

    # 常见列名兼容
    area_col = "Area"
    year_col = "Year"
    item_col = "Item"
    elem_col = "Element"
    unit_col = "Unit"
    val_col  = "Value"

    fao = fao[
        (fao[elem_col] == "Yield") &
        (fao[year_col].between(year_min, year_max)) &
        (fao[unit_col].astype(str).str.lower() == "kg/ha")
    ].copy()

    out = []
    for crop, items in crop_item_map.items():
        tmp = fao[fao[item_col].isin(items)].copy()
        if tmp.empty:
            print(f"[WARN] FAO 中未找到 crop={crop} 对应的 Item: {items}")
            continue

        tmp["ISO_A3"] = tmp[area_col].apply(name_to_iso3)
        tmp = tmp[tmp["ISO_A3"].notna()].copy()

        tmp = tmp[["ISO_A3", year_col, val_col]].rename(
            columns={year_col: "Year", val_col: "fao_yield_kg_ha"}
        )
        tmp["crop"] = crop
        out.append(tmp)

    fao_out = pd.concat(out, ignore_index=True)
    fao_out["Year"] = fao_out["Year"].astype(int)
    fao_out["fao_yield_kg_ha"] = pd.to_numeric(fao_out["fao_yield_kg_ha"], errors="coerce")
    fao_out = fao_out.dropna(subset=["fao_yield_kg_ha"])

    return fao_out


def calc_metrics_log(df, x="model_yield_kg_ha", y="fao_yield_kg_ha"):
    """
    在 log 尺度上计算统计指标
    """
    d = df[[x, y]].dropna().copy()
    d = d[(d[x] > 0) & (d[y] > 0)].copy()

    if len(d) < 3:
        return pd.Series({
            "N": len(d),
            "R_log": np.nan,
            "R2_log": np.nan,
            "RMSE_log": np.nan,
            "MAE_log": np.nan,
            "Bias_log(model-fao)": np.nan
        })

    lx = np.log(d[x].values)
    ly = np.log(d[y].values)

    r = np.corrcoef(lx, ly)[0, 1]
    rmse = np.sqrt(np.mean((lx - ly) ** 2))
    mae  = np.mean(np.abs(lx - ly))
    bias = np.mean(lx - ly)

    return pd.Series({
        "N": len(d),
        "R_log": r,
        "R2_log": r**2,
        "RMSE_log": rmse,
        "MAE_log": mae,
        "Bias_log(model-fao)": bias
    })

def plot_scatter_log(df, crop, out_dir):
    d = df.dropna(subset=["model_yield_kg_ha", "fao_yield_kg_ha"]).copy()
    d = d[(d["model_yield_kg_ha"] > 0) & (d["fao_yield_kg_ha"] > 0)].copy()
    if d.empty:
        return

    x = d["fao_yield_kg_ha"].values
    y = d["model_yield_kg_ha"].values

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(x, y, s=8, alpha=0.35)

    vmin = min(x.min(), y.min()) * 0.9
    vmax = max(x.max(), y.max()) * 1.1

    ax.plot([vmin, vmax], [vmin, vmax], "--", lw=1)

    ax.set_xscale("log")
    ax.set_yscale("log")

    met = calc_metrics_log(d)

    ax.set_xlabel("FAO yield (kg/ha)")
    ax.set_ylabel("Model yield (kg/ha)")
    ax.set_title(f"{crop}: country-year comparison (log-log)")
    ax.text(
        0.03, 0.97,
        f"N={int(met['N'])}\nR(log)={met['R_log']:.2f}\nRMSE(log)={met['RMSE_log']:.2f}\nBias(log)={met['Bias_log(model-fao)']:.2f}",
        transform=ax.transAxes, va="top"
    )

    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, f"scatter_log_{crop}.png"), dpi=300, bbox_inches="tight")
    fig.savefig(os.path.join(out_dir, f"scatter_log_{crop}.pdf"), bbox_inches="tight")
    plt.close(fig)

# =========================
# 4) 构建国家网格对应关系
# =========================
# 用任意一个 crop 的 area 网格作为模板
template_da = DB["dist"]["maiz_nonIrr"]
iso_idx_2d, idx_to_iso = build_grid_country_lookup(template_da, shp_fp, iso_field="ISO_A3")
print("国家网格映射完成。")

# =========================
# 5) 模型端：逐作物聚合到国家-年份
# =========================
model_all = []

for crop in CROPS:
    print(f"Processing model crop: {crop}")

    yield_da, area_da = get_model_actual_yield_da(
        DB, crop,
        rel_scale=RELATIVE_YIELD_SCALE,
        unit_factor=MODEL_YIELD_TO_KGHA
    )

    # 限制年份范围
    years_model = yield_da["year"].values.astype(int)
    keep_years = years_model[(years_model >= YEAR_MIN) & (years_model <= YEAR_MAX)]
    yield_da = yield_da.sel(year=keep_years)

    df_crop = aggregate_country_yield(yield_da, area_da, iso_idx_2d, idx_to_iso)
    df_crop["crop"] = crop
    model_all.append(df_crop)

model_df = pd.concat(model_all, ignore_index=True)
model_df.to_csv(os.path.join(out_dir, "model_country_yield_kg_ha.csv"), index=False)
print("模型国家单产输出完成。")

# =========================
# 6) FAO 端
# =========================
fao_df = load_fao_yield(
    fao_fp,
    crop_item_map=FAO_ITEM_MAP,
    year_min=YEAR_MIN,
    year_max=YEAR_MAX
)
fao_df.to_csv(os.path.join(out_dir, "fao_country_yield_kg_ha.csv"), index=False)
print("FAO 单产读取完成。")

# =========================
# 7) 合并比较
# =========================
merged = pd.merge(
    model_df,
    fao_df,
    on=["ISO_A3", "Year", "crop"],
    how="inner"
).sort_values(["crop", "ISO_A3", "Year"])

merged.to_csv(os.path.join(out_dir, "model_vs_fao_yield_country_year.csv"), index=False)
print("模型与 FAO 匹配完成。")

# =========================
# 8) 总体统计
# =========================
stats_crop = merged.groupby("crop").apply(calc_metrics_log).reset_index()
stats_all  = calc_metrics_log(merged).to_frame().T
stats_all.insert(0, "crop", "all")
stats = pd.concat([stats_crop, stats_all], ignore_index=True)
stats.to_csv(os.path.join(out_dir, "validation_stats_by_crop.csv"), index=False)

print("\n=== Validation stats ===")
print(stats)

# =========================
# 9) 国家层面统计（每个国家在整个时间序列上的表现）
# =========================
country_stats = (
    merged.groupby(["crop", "ISO_A3"])
    .apply(calc_metrics)
    .reset_index()
    .sort_values(["crop", "R2"], ascending=[True, False])
)
country_stats.to_csv(os.path.join(out_dir, "validation_stats_by_country.csv"), index=False)

# =========================
# 10) 画散点图
# =========================
for crop in CROPS:
    plot_scatter(merged[merged["crop"] == crop], crop, out_dir)

print(f"\n全部完成，结果保存于: {out_dir}")

# %% cell 58
# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# =========================
# 0) 全局参数
# =========================
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

out_dir = r"D:\AAUDE\paper_v2\paper4\validation_fao_yield"
os.makedirs(out_dir, exist_ok=True)

# merged 需要包含：
# ISO_A3, Year, crop, model_yield_kg_ha, fao_yield_kg_ha

# =========================
# 1) 每个国家-作物的时间序列相关系数
# =========================
USE_LOG = True     # True: 用 log(yield) 计算 r；False: 用原始尺度
MIN_N   = 10       # 至少多少个年份才计算国家级 r


def calc_country_r(group, use_log=True):
    d = group[["model_yield_kg_ha", "fao_yield_kg_ha"]].dropna().copy()
    d = d[(d["model_yield_kg_ha"] > 0) & (d["fao_yield_kg_ha"] > 0)].copy()

    n = len(d)
    if n < MIN_N:
        return pd.Series({
            "N": n,
            "r": np.nan,
            "mean_model": np.nan,
            "mean_fao": np.nan,
            "bias_mean": np.nan
        })

    x = d["model_yield_kg_ha"].values
    y = d["fao_yield_kg_ha"].values

    if use_log:
        x = np.log(x)
        y = np.log(y)

    # 防止常数序列导致相关系数报错/无意义
    if np.allclose(np.nanstd(x), 0) or np.allclose(np.nanstd(y), 0):
        r = np.nan
    else:
        r = np.corrcoef(x, y)[0, 1]

    return pd.Series({
        "N": n,
        "r": r,
        "mean_model": d["model_yield_kg_ha"].mean(),
        "mean_fao": d["fao_yield_kg_ha"].mean(),
        "bias_mean": d["model_yield_kg_ha"].mean() - d["fao_yield_kg_ha"].mean()
    })


country_r = (
    merged.groupby(["crop", "ISO_A3"], as_index=False)
    .apply(lambda g: calc_country_r(g, use_log=USE_LOG))
    .reset_index(drop=True)
)

country_r = country_r.dropna(subset=["r"]).copy()

# 保存表格
tag = "log" if USE_LOG else "raw"
country_r.to_csv(
    os.path.join(out_dir, f"country_r_by_crop_{tag}.csv"),
    index=False
)

print(country_r.head())
print(country_r.groupby("crop")["ISO_A3"].count())


# =========================
# 2) 2×2 dot plot：每个点一个国家，y = r
# =========================
crop_order = ["maiz", "soyb", "rice", "whea"]
crop_title = {
    "maiz": "Maize",
    "soyb": "Soybean",
    "rice": "Rice",
    "whea": "Wheat"
}

fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharey=True)
axes = axes.flatten()

for ax, crop in zip(axes, crop_order):
    df = country_r[country_r["crop"] == crop].copy()

    if df.empty:
        ax.text(0.5, 0.5, f"No data for {crop}", ha="center", va="center")
        ax.set_axis_off()
        continue

    # 按 r 排序
    df = df.sort_values("r").reset_index(drop=True)
    df["rank"] = np.arange(1, len(df) + 1)

    # 用有效年份数 N 控制点大小
    # 避免太大太小，做一个简单缩放
    sizes = 10 + 2.2 * (df["N"] - df["N"].min())

    # 正负相关区分一下
    pos = df["r"] >= 0
    neg = df["r"] < 0

    ax.scatter(
        df.loc[neg, "rank"], df.loc[neg, "r"],
        s=sizes[neg], alpha=0.7, label="r < 0"
    )
    ax.scatter(
        df.loc[pos, "rank"], df.loc[pos, "r"],
        s=sizes[pos], alpha=0.7, label="r ≥ 0"
    )

    # 中位数线
    med = df["r"].median()
    ax.axhline(med, linestyle="--", linewidth=1)

    # 0 线
    ax.axhline(0, linestyle=":", linewidth=1)

    ax.set_title(f"{crop_title[crop]} (n={len(df)})", fontsize=12)
    ax.set_xlabel("Countries (ranked by r)")
    ax.set_ylabel("Pearson r")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # 只保留左、下边框
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)

    # 角标写中位数
    ax.text(
        0.03, 0.96,
        f"median r = {med:.2f}",
        transform=ax.transAxes,
        ha="left", va="top", fontsize=10
    )

# 统一 y 范围
for ax in axes:
    ax.set_ylim(-1.0, 1.0)

# 统一图例，只放一个
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc="lower center", ncol=2, frameon=False)

plt.tight_layout(rect=[0, 0.05, 1, 1])

suffix = "log" if USE_LOG else "raw"
fig.savefig(os.path.join(out_dir, f"country_r_dotplot_2x2_{suffix}.png"),
            dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(out_dir, f"country_r_dotplot_2x2_{suffix}.pdf"),
            bbox_inches="tight")
plt.show()

# %% cell 59
# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# =========================
# 0) 全局绘图参数
# =========================
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 1) 路径
# =========================
in_csv  = r"D:\AAUDE\paper_v2\paper4\validation_fao_yield\model_vs_fao_yield_country_year.csv"
out_dir = r"D:\AAUDE\paper_v2\paper4\validation_fao_yield\fig_country_r_2x2_log"
os.makedirs(out_dir, exist_ok=True)

# =========================
# 2) 参数
# =========================
crops = ["maiz", "soyb", "rice", "whea"]
crop_title = {
    "maiz": "Maize",
    "soyb": "Soybean",
    "rice": "Rice",
    "whea": "Wheat"
}

# 至少多少个年份才计算国家r
MIN_N_YEAR = 10

# 是否对均值先做下限裁剪，避免log(0)
EPS = 1e-6

# 点大小
POINT_SIZE = 28

# 透明度
POINT_ALPHA = 0.85

# =========================
# 3) 读取数据
# =========================
df = pd.read_csv(in_csv)

# 基本清洗
df = df.copy()
df["Year"] = pd.to_numeric(df["Year"], errors="coerce")
df["model_yield_kg_ha"] = pd.to_numeric(df["model_yield_kg_ha"], errors="coerce")
df["fao_yield_kg_ha"]   = pd.to_numeric(df["fao_yield_kg_ha"], errors="coerce")

df = df.dropna(subset=["ISO_A3", "Year", "crop", "model_yield_kg_ha", "fao_yield_kg_ha"])
df = df[(df["model_yield_kg_ha"] > 0) & (df["fao_yield_kg_ha"] > 0)].copy()

# =========================
# 4) 计算每个国家的 r、均值、样本数
# =========================
def calc_country_stats_one(group):
    """
    对单个国家-作物的时间序列：
    - r = Pearson相关系数（model vs fao, over years）
    - model_mean = 该国多年平均模型单产
    - fao_mean   = 该国多年平均FAO单产
    - n_year     = 重叠年份数
    """
    g = group.sort_values("Year").copy()

    x = g["fao_yield_kg_ha"].values
    y = g["model_yield_kg_ha"].values
    n = len(g)

    if n < MIN_N_YEAR:
        r = np.nan
    else:
        # 若序列方差为0，相关系数不可定义
        if np.nanstd(x) == 0 or np.nanstd(y) == 0:
            r = np.nan
        else:
            r = np.corrcoef(x, y)[0, 1]

    return pd.Series({
        "n_year": n,
        "fao_mean": np.nanmean(x),
        "model_mean": np.nanmean(y),
        "r": r
    })

country_stats = (
    df.groupby(["crop", "ISO_A3"], as_index=False)
      .apply(calc_country_stats_one)
      .reset_index(drop=True)
)

# 仅保留满足要求的数据
country_stats = country_stats.dropna(subset=["fao_mean", "model_mean", "r"]).copy()
country_stats = country_stats[
    (country_stats["fao_mean"] > 0) &
    (country_stats["model_mean"] > 0)
].copy()

# log 用的列（这里只是为了稳妥，也可直接用 set_xscale/set_yscale）
country_stats["fao_mean_clip"]   = country_stats["fao_mean"].clip(lower=EPS)
country_stats["model_mean_clip"] = country_stats["model_mean"].clip(lower=EPS)

# =========================
# 5) 统计每个作物的整体表现（基于国家均值点）
# =========================
def calc_scatter_metrics(sub):
    """
    基于每个国家的平均值点，计算散点层面的统计量
    """
    if len(sub) < 3:
        return {"N": len(sub), "R_mean": np.nan, "RMSE_log10": np.nan}

    x = np.log10(sub["fao_mean_clip"].values)
    y = np.log10(sub["model_mean_clip"].values)

    if np.nanstd(x) == 0 or np.nanstd(y) == 0:
        rr = np.nan
    else:
        rr = np.corrcoef(x, y)[0, 1]

    rmse_log = np.sqrt(np.mean((y - x) ** 2))
    return {"N": len(sub), "R_mean": rr, "RMSE_log10": rmse_log}

# =========================
# 6) 为四个作物统一坐标范围
# =========================
all_x = country_stats["fao_mean_clip"].values
all_y = country_stats["model_mean_clip"].values

global_min = np.nanmin(np.r_[all_x, all_y])
global_max = np.nanmax(np.r_[all_x, all_y])

# 稍微留白
xmin = 10 ** np.floor(np.log10(global_min))
xmax = 10 ** np.ceil(np.log10(global_max))

# =========================
# 7) 绘制 2x2 图
# =========================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.ravel()

# 颜色条范围固定到[-1, 1]
vmin, vmax = -1, 1
cmap = plt.cm.RdBu_r

for i, crop in enumerate(crops):
    ax = axes[i]
    sub = country_stats[country_stats["crop"] == crop].copy()

    if sub.empty:
        ax.text(0.5, 0.5, f"{crop_title[crop]}\nNo data",
                ha="center", va="center", transform=ax.transAxes, fontsize=12)
        ax.set_axis_off()
        continue

    sc = ax.scatter(
        sub["fao_mean_clip"],
        sub["model_mean_clip"],
        c=sub["r"],
        s=POINT_SIZE,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        alpha=POINT_ALPHA,
        edgecolors="black",
        linewidths=0.25
    )

    # 1:1 参考线
    ax.plot([xmin, xmax], [xmin, xmax], "--", color="gray", lw=1.0, zorder=0)

    # log 坐标
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(xmin, xmax)

    # 标题
    ax.set_title(crop_title[crop], fontsize=13)

    # 只保留左、下边框
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.0)
    ax.spines["bottom"].set_linewidth(1.0)

    # 不要网格线
    ax.grid(False)

    # 统计文字
    met = calc_scatter_metrics(sub)
    med_r = np.nanmedian(sub["r"].values) if len(sub) > 0 else np.nan

    txt = (
        f"N countries = {int(met['N'])}\n"
        f"Median country r = {med_r:.2f}\n"
        f"log10-space R = {met['R_mean']:.2f}\n"
        f"log10-space RMSE = {met['RMSE_log10']:.2f}"
    )
    ax.text(
        0.04, 0.96, txt,
        transform=ax.transAxes,
        ha="left", va="top",
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7", alpha=0.9)
    )

    # 轴标签
    if i in [2, 3]:
        ax.set_xlabel("FAO mean yield (kg ha$^{-1}$)", fontsize=11)
    if i in [0, 2]:
        ax.set_ylabel("Model mean yield (kg ha$^{-1}$)", fontsize=11)

# 布局
plt.subplots_adjust(wspace=0.22, hspace=0.22, right=0.88)

# 公共颜色条
cax = fig.add_axes([0.90, 0.16, 0.02, 0.68])
cb = fig.colorbar(
    plt.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax), cmap=cmap),
    cax=cax
)
cb.set_label("Country-level Pearson r\n(model vs FAO over years)", fontsize=11)

# 总标题（可删）
fig.suptitle("Country-level validation of crop yields against FAO", fontsize=14, y=0.98)

# 保存
png_fp = os.path.join(out_dir, "validation_country_r_2x2_log.png")
pdf_fp = os.path.join(out_dir, "validation_country_r_2x2_log.pdf")

fig.savefig(png_fp, dpi=300, bbox_inches="tight")
fig.savefig(pdf_fp, bbox_inches="tight")
plt.close(fig)

print("Done.")
print(f"Figure saved to:\n{png_fp}\n{pdf_fp}")

# =========================
# 8) 同时导出国家统计表
# =========================
csv_fp = os.path.join(out_dir, "country_validation_stats_for_plot.csv")
country_stats.to_csv(csv_fp, index=False)
print(f"Country stats saved to:\n{csv_fp}")

# %% cell 60
# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import pycountry

# =========================
# 0) 全局绘图参数
# =========================
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 1) 路径
# =========================
yield_compare_fp = r"D:\AAUDE\paper_v2\paper4\validation_fao_yield\model_vs_fao_yield_country_year.csv"
prod_fp          = r"D:\Edge_download\FAOSTAT_data_en_3-8-2026 (1).csv"
out_dir          = r"D:\AAUDE\paper_v2\paper4\validation_fao_yield\top5_producers_timeseries"
os.makedirs(out_dir, exist_ok=True)

# =========================
# 2) 基本参数
# =========================
crops = ["maiz", "soyb", "rice", "whea"]
crop_title = {
    "maiz": "Maize",
    "soyb": "Soybean",
    "rice": "Rice",
    "whea": "Wheat"
}

# 生产数据中 FAO Item 名称映射
FAO_ITEM_MAP = {
    "maiz": ["Maize (corn)"],
    "soyb": ["Soya beans", "Soybeans"],
    "rice": ["Rice"],
    "whea": ["Wheat"]
}

MIN_N_YEAR = 4   # 至少多少个重叠年份才计算 r / RMSE

# =========================
# 3) 国家名称转 ISO3
# =========================
def name_to_iso3(name):
    manual = {
        "Bolivia (Plurinational State of)": "BOL",
        "Cabo Verde": "CPV",
        "China, mainland": "CHN",
        "China, Taiwan Province of": "TWN",
        "China, Hong Kong SAR": "HKG",
        "China, Macao SAR": "MAC",
        "Côte d'Ivoire": "CIV",
        "Democratic People's Republic of Korea": "PRK",
        "Democratic Republic of the Congo": "COD",
        "Iran (Islamic Republic of)": "IRN",
        "Lao People's Democratic Republic": "LAO",
        "Micronesia (Federated States of)": "FSM",
        "Republic of Korea": "KOR",
        "Republic of Moldova": "MDA",
        "Russian Federation": "RUS",
        "Syrian Arab Republic": "SYR",
        "Türkiye": "TUR",
        "United Kingdom of Great Britain and Northern Ireland": "GBR",
        "United Republic of Tanzania": "TZA",
        "United States of America": "USA",
        "Venezuela (Bolivarian Republic of)": "VEN",
        "Viet Nam": "VNM",
        "Brunei Darussalam": "BRN",
        "Eswatini": "SWZ",
        "Palestine": "PSE",
    }
    if name in manual:
        return manual[name]

    try:
        return pycountry.countries.lookup(name).alpha_3
    except Exception:
        return np.nan

# =========================
# 4) 读取模型-vs-FAO单产对比表
# =========================
df_y = pd.read_csv(yield_compare_fp)
df_y["Year"] = pd.to_numeric(df_y["Year"], errors="coerce")
df_y["model_yield_kg_ha"] = pd.to_numeric(df_y["model_yield_kg_ha"], errors="coerce")
df_y["fao_yield_kg_ha"]   = pd.to_numeric(df_y["fao_yield_kg_ha"], errors="coerce")

df_y = df_y.dropna(subset=["ISO_A3", "Year", "crop", "model_yield_kg_ha", "fao_yield_kg_ha"]).copy()
df_y = df_y[(df_y["model_yield_kg_ha"] > 0) & (df_y["fao_yield_kg_ha"] > 0)].copy()
df_y["Year"] = df_y["Year"].astype(int)

# =========================
# 5) 读取FAO生产数据，用于筛选主产国
# =========================
prod = pd.read_csv(prod_fp)

# 清理列名两端空格
prod.columns = [c.strip() for c in prod.columns]

# 兼容列名
area_col = "Area"
item_col = "Item"
elem_col = "Element"
year_col = "Year"
unit_col = "Unit"
val_col  = "Value"

prod = prod.copy()
prod[year_col] = pd.to_numeric(prod[year_col], errors="coerce")
prod[val_col]  = pd.to_numeric(prod[val_col], errors="coerce")

prod = prod[
    (prod[elem_col] == "Production") &
    (prod[unit_col].astype(str).str.lower() == "t")
].copy()

prod["ISO_A3"] = prod[area_col].apply(name_to_iso3)
prod = prod.dropna(subset=["ISO_A3", year_col, val_col]).copy()
prod[year_col] = prod[year_col].astype(int)

# =========================
# 6) 指标函数
# =========================
def calc_r_rmse(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]

    if len(x) < MIN_N_YEAR:
        return np.nan, np.nan, len(x)

    if np.std(x) == 0 or np.std(y) == 0:
        r = np.nan
    else:
        r = np.corrcoef(x, y)[0, 1]

    rmse = np.sqrt(np.mean((y - x) ** 2))
    return r, rmse, len(x)

# =========================
# 7) 为每个作物筛选“前5主产国”
#    关键：只用与模型单产有重叠的年份
# =========================
top5_by_crop = {}

for crop in crops:
    # 模型-vs-FAO 单产里实际存在的重叠年份
    yrs_overlap = sorted(df_y.loc[df_y["crop"] == crop, "Year"].dropna().unique().tolist())
    if len(yrs_overlap) == 0:
        print(f"[WARN] {crop} 没有可用的单产对比年份，跳过。")
        top5_by_crop[crop] = []
        continue

    items = FAO_ITEM_MAP[crop]
    psub = prod[
        (prod[item_col].isin(items)) &
        (prod[year_col].isin(yrs_overlap))
    ].copy()

    if psub.empty:
        print(f"[WARN] {crop} 在生产数据中未找到对应记录。")
        top5_by_crop[crop] = []
        continue

    rank_df = (
        psub.groupby("ISO_A3", as_index=False)[val_col]
            .mean()
            .rename(columns={val_col: "mean_prod_t"})
            .sort_values("mean_prod_t", ascending=False)
    )

    top5 = rank_df.head(5)["ISO_A3"].tolist()
    top5_by_crop[crop] = top5
    print(f"{crop}: top5 producers over overlap years = {top5}")

# 导出排名表
rank_out = []
for crop, iso_list in top5_by_crop.items():
    for i, iso in enumerate(iso_list, start=1):
        rank_out.append([crop, i, iso])
pd.DataFrame(rank_out, columns=["crop", "rank", "ISO_A3"]).to_csv(
    os.path.join(out_dir, "top5_producers_by_crop.csv"), index=False
)

# =========================
# 8) 画图：每个作物一张图，5行1列
# =========================
for crop in crops:
    top5 = top5_by_crop.get(crop, [])
    if len(top5) == 0:
        continue

    sub_crop = df_y[df_y["crop"] == crop].copy()

    # 只保留前5国家且按年份排序
    sub_crop = sub_crop[sub_crop["ISO_A3"].isin(top5)].copy()
    if sub_crop.empty:
        continue

    # 统一y范围，避免5个子图尺度不同难比较
    y_all = np.r_[np.log(sub_crop["fao_yield_kg_ha"].values), np.log(sub_crop["model_yield_kg_ha"].values)]
    y_all = y_all[np.isfinite(y_all)]
    ymin = max(0, np.nanmin(y_all) * 0.9)
    ymax = np.nanmax(y_all) * 1.1

    fig, axes = plt.subplots(len(top5), 1, figsize=(10, 2.6 * len(top5)), sharex=True)
    if len(top5) == 1:
        axes = [axes]

    for ax, iso in zip(axes, top5):
        g = sub_crop[sub_crop["ISO_A3"] == iso].sort_values("Year").copy()
        if g.empty:
            ax.text(0.5, 0.5, f"{iso}: No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            continue

        r, rmse, n = calc_r_rmse(np.log(g["fao_yield_kg_ha"].values), np.log(g["model_yield_kg_ha"].values))

        # 画线
        ax.plot(
            g["Year"], np.log(g["fao_yield_kg_ha"]),
            lw=2.0, label="FAO", marker="o", markersize=3.5
        )
        ax.plot(
            g["Year"], np.log(g["model_yield_kg_ha"]),
            lw=2.0, label="Model", marker="s", markersize=3.2
        )

        # 文字说明
        txt = f"{iso}   r = {r:.2f}   RMSE = {rmse:.1f}   n = {n}"
        ax.text(
            0.01, 0.92, txt,
            transform=ax.transAxes,
            ha="left", va="top", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.75", alpha=0.9)
        )

        # 轴样式
        ax.set_ylim(ymin, ymax)
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(1.0)
        ax.spines["bottom"].set_linewidth(1.0)
        ax.set_ylabel("Yield\n(kg ha$^{-1}$)", fontsize=10)

    axes[0].legend(frameon=False, ncol=2, loc="upper right")
    axes[-1].set_xlabel("Year", fontsize=11)

    fig.suptitle(f"{crop_title[crop]}: top 5 producing countries", fontsize=14, y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.985])

    png_fp = os.path.join(out_dir, f"{crop}_top5_producers_yield_timeseries.png")
    pdf_fp = os.path.join(out_dir, f"{crop}_top5_producers_yield_timeseries.pdf")
    fig.savefig(png_fp, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_fp, bbox_inches="tight")
    plt.close(fig)

# =========================
# 9) 汇总表：每个作物-国家的r、RMSE
# =========================
summary = []
for crop in crops:
    top5 = top5_by_crop.get(crop, [])
    sub_crop = df_y[(df_y["crop"] == crop) & (df_y["ISO_A3"].isin(top5))].copy()
    for iso in top5:
        g = sub_crop[sub_crop["ISO_A3"] == iso].sort_values("Year").copy()
        r, rmse, n = calc_r_rmse(g["fao_yield_kg_ha"].values, g["model_yield_kg_ha"].values)
        summary.append([crop, iso, n, r, rmse])

summary_df = pd.DataFrame(summary, columns=["crop", "ISO_A3", "n_year", "r", "RMSE_kg_ha"])
summary_df.to_csv(os.path.join(out_dir, "top5_producers_yield_validation_summary.csv"), index=False)

print("Done.")
print(f"Outputs saved to: {out_dir}")

# %% cell 61
# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import pycountry

# =========================
# 0) 全局绘图参数
# =========================
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 1) 路径
# =========================
yield_compare_fp = r"D:\AAUDE\paper_v2\paper4\validation_fao_yield\model_vs_fao_yield_country_year.csv"
prod_fp          = r"D:\Edge_download\FAOSTAT_data_en_3-8-2026 (1).csv"
out_dir          = r"D:\AAUDE\paper_v2\paper4\validation_fao_yield\top5_producers_anomaly_timeseries"
os.makedirs(out_dir, exist_ok=True)

# =========================
# 2) 参数
# =========================
crops = ["maiz", "soyb", "rice", "whea"]
crop_title = {
    "maiz": "Maize",
    "soyb": "Soybean",
    "rice": "Rice",
    "whea": "Wheat"
}

FAO_ITEM_MAP = {
    "maiz": ["Maize (corn)"],
    "soyb": ["Soya beans", "Soybeans"],
    "rice": ["Rice"],
    "whea": ["Wheat"]
}

MIN_N_YEAR = 4

# =========================
# 3) 国家名转 ISO3
# =========================
def name_to_iso3(name):
    manual = {
        "Bolivia (Plurinational State of)": "BOL",
        "Cabo Verde": "CPV",
        "China, mainland": "CHN",
        "China, Taiwan Province of": "TWN",
        "China, Hong Kong SAR": "HKG",
        "China, Macao SAR": "MAC",
        "Côte d'Ivoire": "CIV",
        "Democratic People's Republic of Korea": "PRK",
        "Democratic Republic of the Congo": "COD",
        "Iran (Islamic Republic of)": "IRN",
        "Lao People's Democratic Republic": "LAO",
        "Micronesia (Federated States of)": "FSM",
        "Republic of Korea": "KOR",
        "Republic of Moldova": "MDA",
        "Russian Federation": "RUS",
        "Syrian Arab Republic": "SYR",
        "Türkiye": "TUR",
        "United Kingdom of Great Britain and Northern Ireland": "GBR",
        "United Republic of Tanzania": "TZA",
        "United States of America": "USA",
        "Venezuela (Bolivarian Republic of)": "VEN",
        "Viet Nam": "VNM",
        "Brunei Darussalam": "BRN",
        "Eswatini": "SWZ",
        "Palestine": "PSE",
    }
    if name in manual:
        return manual[name]
    try:
        return pycountry.countries.lookup(name).alpha_3
    except Exception:
        return np.nan

# =========================
# 4) 指标函数：基于距平
# =========================
def calc_r_rmse(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]

    if len(x) < MIN_N_YEAR:
        return np.nan, np.nan, len(x)

    if np.std(x) == 0 or np.std(y) == 0:
        r = np.nan
    else:
        r = np.corrcoef(x, y)[0, 1]

    rmse = np.sqrt(np.mean((y - x) ** 2))
    return r, rmse, len(x)

def add_anomaly(group):
    g = group.sort_values("Year").copy()
    g["fao_anom"]   = np.log(g["fao_yield_kg_ha"])   - np.log(g["fao_yield_kg_ha"]).mean()
    g["model_anom"] = np.log(g["model_yield_kg_ha"]) - np.log(g["model_yield_kg_ha"]).mean()
    return g

# =========================
# 5) 读取模型-vs-FAO 单产对比表
# =========================
df_y = pd.read_csv(yield_compare_fp)
df_y["Year"] = pd.to_numeric(df_y["Year"], errors="coerce")
df_y["model_yield_kg_ha"] = pd.to_numeric(df_y["model_yield_kg_ha"], errors="coerce")
df_y["fao_yield_kg_ha"]   = pd.to_numeric(df_y["fao_yield_kg_ha"], errors="coerce")

df_y = df_y.dropna(subset=["ISO_A3", "Year", "crop", "model_yield_kg_ha", "fao_yield_kg_ha"]).copy()
df_y = df_y[(df_y["model_yield_kg_ha"] > 0) & (df_y["fao_yield_kg_ha"] > 0)].copy()
df_y["Year"] = df_y["Year"].astype(int)

# =========================
# 6) 读取 FAO production，用于选前5主产国
# =========================
prod = pd.read_csv(prod_fp)
prod.columns = [c.strip() for c in prod.columns]

area_col = "Area"
item_col = "Item"
elem_col = "Element"
year_col = "Year"
unit_col = "Unit"
val_col  = "Value"

prod[year_col] = pd.to_numeric(prod[year_col], errors="coerce")
prod[val_col]  = pd.to_numeric(prod[val_col], errors="coerce")

prod = prod[
    (prod[elem_col] == "Production") &
    (prod[unit_col].astype(str).str.lower() == "t")
].copy()

prod["ISO_A3"] = prod[area_col].apply(name_to_iso3)
prod = prod.dropna(subset=["ISO_A3", year_col, val_col]).copy()
prod[year_col] = prod[year_col].astype(int)

# =========================
# 7) 按重叠年份筛前5主产国
# =========================
top5_by_crop = {}

for crop in crops:
    yrs_overlap = sorted(df_y.loc[df_y["crop"] == crop, "Year"].dropna().unique().tolist())
    if len(yrs_overlap) == 0:
        top5_by_crop[crop] = []
        continue

    items = FAO_ITEM_MAP[crop]
    psub = prod[
        (prod[item_col].isin(items)) &
        (prod[year_col].isin(yrs_overlap))
    ].copy()

    if psub.empty:
        top5_by_crop[crop] = []
        continue

    rank_df = (
        psub.groupby("ISO_A3", as_index=False)[val_col]
            .mean()
            .rename(columns={val_col: "mean_prod_t"})
            .sort_values("mean_prod_t", ascending=False)
    )

    top5_by_crop[crop] = rank_df.head(5)["ISO_A3"].tolist()
    print(f"{crop}: {top5_by_crop[crop]}")

# 导出 top5 列表
rank_out = []
for crop, isos in top5_by_crop.items():
    for i, iso in enumerate(isos, start=1):
        rank_out.append([crop, i, iso])
pd.DataFrame(rank_out, columns=["crop", "rank", "ISO_A3"]).to_csv(
    os.path.join(out_dir, "top5_producers_by_crop.csv"), index=False
)

# =========================
# 8) 添加距平
# =========================
df_anom = (
    df_y.groupby(["crop", "ISO_A3"], group_keys=False)
        .apply(add_anomaly)
        .reset_index(drop=True)
)

# =========================
# 9) 逐作物画图：5行1列，每个国家一个子图
# =========================
summary = []

for crop in crops:
    top5 = top5_by_crop.get(crop, [])
    if len(top5) == 0:
        continue

    sub = df_anom[(df_anom["crop"] == crop) & (df_anom["ISO_A3"].isin(top5))].copy()
    if sub.empty:
        continue

    # 统一该作物所有子图的y轴范围
    y_all = np.r_[sub["fao_anom"].values, sub["model_anom"].values]
    y_all = y_all[np.isfinite(y_all)]
    ymax_abs = np.nanmax(np.abs(y_all)) * 1.15
    ymin, ymax = -ymax_abs, ymax_abs

    fig, axes = plt.subplots(len(top5), 1, figsize=(10, 2.6 * len(top5)), sharex=True)
    if len(top5) == 1:
        axes = [axes]

    for ax, iso in zip(axes, top5):
        g = sub[sub["ISO_A3"] == iso].sort_values("Year").copy()
        if g.empty:
            ax.text(0.5, 0.5, f"{iso}: No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            continue

        r, rmse, n = calc_r_rmse(g["fao_anom"].values, g["model_anom"].values)

        # 导出汇总
        summary.append([crop, iso, n, r, rmse])

        # 零线
        ax.axhline(0, ls="--", lw=1.0, color="gray", zorder=0)

        # 距平时间序列
        ax.plot(
            g["Year"], g["fao_anom"],
            lw=2.0, marker="o", markersize=3.5, label="FAO anomaly"
        )
        ax.plot(
            g["Year"], g["model_anom"],
            lw=2.0, marker="s", markersize=3.2, label="Model anomaly"
        )

        txt = f"{iso}   r = {r:.2f}   RMSE = {rmse:.1f} kg ha$^{{-1}}$   n = {n}"
        ax.text(
            0.01, 0.92, txt,
            transform=ax.transAxes,
            ha="left", va="top", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.75", alpha=0.9)
        )

        ax.set_ylim(ymin, ymax)
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(1.0)
        ax.spines["bottom"].set_linewidth(1.0)
        ax.set_ylabel("Yield anomaly\n(kg ha$^{-1}$)", fontsize=10)

    axes[0].legend(frameon=False, ncol=2, loc="upper right")
    axes[-1].set_xlabel("Year", fontsize=11)

    fig.suptitle(f"{crop_title[crop]}: yield anomaly in top 5 producing countries", fontsize=14, y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.985])

    png_fp = os.path.join(out_dir, f"{crop}_top5_producers_yield_anomaly_timeseries.png")
    pdf_fp = os.path.join(out_dir, f"{crop}_top5_producers_yield_anomaly_timeseries.pdf")
    fig.savefig(png_fp, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_fp, bbox_inches="tight")
    plt.close(fig)

# =========================
# 10) 导出汇总表
# =========================
summary_df = pd.DataFrame(summary, columns=["crop", "ISO_A3", "n_year", "r", "RMSE_anomaly_kg_ha"])
summary_df.to_csv(os.path.join(out_dir, "top5_producers_yield_anomaly_validation_summary.csv"), index=False)

print("Done.")
print(f"Outputs saved to: {out_dir}")

# %% cell 62
# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import pycountry

# =========================
# 0) 全局绘图参数
# =========================
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 1) 路径
# =========================
yield_compare_fp = r"D:\AAUDE\paper_v2\paper4\validation_fao_yield\model_vs_fao_yield_country_year.csv"
prod_fp          = r"D:\Edge_download\FAOSTAT_data_en_3-8-2026 (1).csv"
out_dir          = r"D:\AAUDE\paper_v2\paper4\validation_fao_yield\top5_producers_anomaly_timeseries_1980_2019"
os.makedirs(out_dir, exist_ok=True)

# =========================
# 2) 参数
# =========================
crops = ["maiz", "soyb", "rice", "whea"]
crop_title = {
    "maiz": "Maize",
    "soyb": "Soybean",
    "rice": "Rice",
    "whea": "Wheat"
}

FAO_ITEM_MAP = {
    "maiz": ["Maize (corn)"],
    "soyb": ["Soya beans", "Soybeans"],
    "rice": ["Rice"],
    "whea": ["Wheat"]
}

MIN_N_YEAR  = 8
VALID_START = 1980
VALID_END   = 2019

# =========================
# 3) 国家名转 ISO3
# =========================
def name_to_iso3(name):
    manual = {
        "Bolivia (Plurinational State of)": "BOL",
        "Cabo Verde": "CPV",
        "China, mainland": "CHN",
        "China, Taiwan Province of": "TWN",
        "China, Hong Kong SAR": "HKG",
        "China, Macao SAR": "MAC",
        "Côte d'Ivoire": "CIV",
        "Democratic People's Republic of Korea": "PRK",
        "Democratic Republic of the Congo": "COD",
        "Iran (Islamic Republic of)": "IRN",
        "Lao People's Democratic Republic": "LAO",
        "Micronesia (Federated States of)": "FSM",
        "Republic of Korea": "KOR",
        "Republic of Moldova": "MDA",
        "Russian Federation": "RUS",
        "Syrian Arab Republic": "SYR",
        "Türkiye": "TUR",
        "United Kingdom of Great Britain and Northern Ireland": "GBR",
        "United Republic of Tanzania": "TZA",
        "United States of America": "USA",
        "Venezuela (Bolivarian Republic of)": "VEN",
        "Viet Nam": "VNM",
        "Brunei Darussalam": "BRN",
        "Eswatini": "SWZ",
        "Palestine": "PSE",
    }
    if name in manual:
        return manual[name]
    try:
        return pycountry.countries.lookup(name).alpha_3
    except Exception:
        return np.nan

# =========================
# 4) 指标函数：基于普通距平
# =========================
def calc_r_rmse(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]

    if len(x) < MIN_N_YEAR:
        return np.nan, np.nan, len(x)

    if np.std(x) == 0 or np.std(y) == 0:
        r = np.nan
    else:
        r = np.corrcoef(x, y)[0, 1]

    rmse = np.sqrt(np.mean((y - x) ** 2))
    return r, rmse, len(x)

def add_anomaly(group):
    g = group.sort_values("Year").copy()
    g["fao_anom"]   = np.log(g["fao_yield_kg_ha"])   - np.log(g["fao_yield_kg_ha"]).mean()
    g["model_anom"] = np.log(g["model_yield_kg_ha"]) - np.log(g["model_yield_kg_ha"]).mean()
    return g

# =========================
# 5) 读取模型-vs-FAO 单产对比表
#    这里只保留 1980-2019
# =========================
df_y = pd.read_csv(yield_compare_fp)
df_y["Year"] = pd.to_numeric(df_y["Year"], errors="coerce")
df_y["model_yield_kg_ha"] = pd.to_numeric(df_y["model_yield_kg_ha"], errors="coerce")
df_y["fao_yield_kg_ha"]   = pd.to_numeric(df_y["fao_yield_kg_ha"], errors="coerce")

df_y = df_y.dropna(subset=["ISO_A3", "Year", "crop", "model_yield_kg_ha", "fao_yield_kg_ha"]).copy()
df_y = df_y[(df_y["model_yield_kg_ha"] > 0) & (df_y["fao_yield_kg_ha"] > 0)].copy()
df_y["Year"] = df_y["Year"].astype(int)

df_y = df_y[
    (df_y["Year"] >= VALID_START) &
    (df_y["Year"] <= VALID_END)
].copy()

# =========================
# 6) 读取 FAO production，仅用于识别主产国
# =========================
prod = pd.read_csv(prod_fp)
prod.columns = [c.strip() for c in prod.columns]

area_col = "Area"
item_col = "Item"
elem_col = "Element"
year_col = "Year"
unit_col = "Unit"
val_col  = "Value"

prod[year_col] = pd.to_numeric(prod[year_col], errors="coerce")
prod[val_col]  = pd.to_numeric(prod[val_col], errors="coerce")

prod = prod[
    (prod[elem_col] == "Production") &
    (prod[unit_col].astype(str).str.lower() == "t")
].copy()

prod["ISO_A3"] = prod[area_col].apply(name_to_iso3)
prod = prod.dropna(subset=["ISO_A3", year_col, val_col]).copy()
prod[year_col] = prod[year_col].astype(int)

# =========================
# 7) 按“验证年份与production年份的交集”筛前5主产国
#    因为你的 production 文件是 2014-2024，所以这里实际上多半是 2014-2019
# =========================
top5_by_crop = {}
rank_records = []

for crop in crops:
    # 单产验证实际年份
    yrs_yield = sorted(df_y.loc[df_y["crop"] == crop, "Year"].dropna().unique().tolist())
    if len(yrs_yield) == 0:
        top5_by_crop[crop] = []
        continue

    items = FAO_ITEM_MAP[crop]

    # production 只能在它自己存在的年份里筛
    psub = prod[
        (prod[item_col].isin(items)) &
        (prod[year_col].isin(yrs_yield))
    ].copy()

    if psub.empty:
        print(f"[WARN] {crop} 在 production 文件中无可用记录。")
        top5_by_crop[crop] = []
        continue

    rank_df = (
        psub.groupby("ISO_A3", as_index=False)[val_col]
            .mean()
            .rename(columns={val_col: "mean_prod_t"})
            .sort_values("mean_prod_t", ascending=False)
    )

    top5 = rank_df.head(5)["ISO_A3"].tolist()
    top5_by_crop[crop] = top5

    for i, row in enumerate(rank_df.head(5).itertuples(index=False), start=1):
        rank_records.append([crop, i, row.ISO_A3, row.mean_prod_t])

    print(f"{crop}: top5 producers = {top5}")

rank_df_out = pd.DataFrame(rank_records, columns=["crop", "rank", "ISO_A3", "mean_prod_t"])
rank_df_out.to_csv(os.path.join(out_dir, "top5_producers_by_crop.csv"), index=False)

# =========================
# 8) 添加距平（基于1980-2019验证序列本身）
# =========================
df_anom = (
    df_y.groupby(["crop", "ISO_A3"], group_keys=False)
        .apply(add_anomaly)
        .reset_index(drop=True)
)

# =========================
# 9) 逐作物画图：5行1列，每个国家一个子图
# =========================
summary = []

for crop in crops:
    top5 = top5_by_crop.get(crop, [])
    if len(top5) == 0:
        continue

    sub = df_anom[(df_anom["crop"] == crop) & (df_anom["ISO_A3"].isin(top5))].copy()
    if sub.empty:
        continue

    # 统一该作物所有子图的y轴范围
    y_all = np.r_[sub["fao_anom"].values, sub["model_anom"].values]
    y_all = y_all[np.isfinite(y_all)]
    ymax_abs = np.nanmax(np.abs(y_all)) * 1.15
    ymin, ymax = -ymax_abs, ymax_abs

    fig, axes = plt.subplots(len(top5), 1, figsize=(10, 2.5 * len(top5)), sharex=True)
    if len(top5) == 1:
        axes = [axes]

    for ax, iso in zip(axes, top5):
        g = sub[sub["ISO_A3"] == iso].sort_values("Year").copy()
        if g.empty:
            ax.text(0.5, 0.5, f"{iso}: No data", ha="center", va="center", transform=ax.transAxes)
            ax.set_axis_off()
            continue

        r, rmse, n = calc_r_rmse(g["fao_anom"].values, g["model_anom"].values)
        summary.append([crop, iso, n, r, rmse])

        # 零线
        ax.axhline(0, ls="--", lw=1.0, color="gray", zorder=0)

        # 距平序列
        ax.plot(
            g["Year"], g["fao_anom"],
            lw=2.0, marker="o", markersize=3.0, label="FAO anomaly"
        )
        ax.plot(
            g["Year"], g["model_anom"],
            lw=2.0, marker="s", markersize=2.8, label="Model anomaly"
        )

        txt = f"{iso}   r = {r:.2f}   RMSE = {rmse:.1f} kg ha$^{{-1}}$   n = {n}"
        ax.text(
            0.01, 0.92, txt,
            transform=ax.transAxes,
            ha="left", va="top", fontsize=10,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="0.75", alpha=0.9)
        )

        ax.set_ylim(ymin, ymax)
        ax.grid(False)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["left"].set_linewidth(1.0)
        ax.spines["bottom"].set_linewidth(1.0)
        ax.set_ylabel("Yield anomaly\n(kg ha$^{-1}$)", fontsize=10)

    axes[0].legend(frameon=False, ncol=2, loc="upper right")
    axes[-1].set_xlabel("Year", fontsize=11)

    fig.suptitle(
        f"{crop_title[crop]}: yield anomaly comparison in top 5 producing countries ({VALID_START}–{VALID_END})",
        fontsize=13, y=0.995
    )
    plt.tight_layout(rect=[0, 0, 1, 0.985])

    png_fp = os.path.join(out_dir, f"{crop}_top5_producers_yield_anomaly_timeseries_{VALID_START}_{VALID_END}.png")
    pdf_fp = os.path.join(out_dir, f"{crop}_top5_producers_yield_anomaly_timeseries_{VALID_START}_{VALID_END}.pdf")
    fig.savefig(png_fp, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_fp, bbox_inches="tight")
    plt.close(fig)

# =========================
# 10) 导出汇总表
# =========================
summary_df = pd.DataFrame(summary, columns=["crop", "ISO_A3", "n_year", "r", "RMSE_anomaly_kg_ha"])
summary_df.to_csv(
    os.path.join(out_dir, f"top5_producers_yield_anomaly_validation_summary_{VALID_START}_{VALID_END}.csv"),
    index=False
)

print("Done.")
print(f"Outputs saved to: {out_dir}")

# %% cell 63
# -*- coding: utf-8 -*-
import os
import numpy as np
import xarray as xr
import rioxarray as rxr
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 1) 路径
# =========================
gaez_fp = r"F:\worldpop\mze_2000_prd.tif"
out_dir = r"D:\AAUDE\paper_v2\paper4\gaez_compare_2000"
os.makedirs(out_dir, exist_ok=True)

# =========================
# 2) 读取 GAEZ tif
# =========================
gaez = rxr.open_rasterio(gaez_fp).squeeze(drop=True)

# 若坐标名为 x/y，改为 lon/lat
rename_dict = {}
if "x" in gaez.dims:
    rename_dict["x"] = "lon"
if "y" in gaez.dims:
    rename_dict["y"] = "lat"
if rename_dict:
    gaez = gaez.rename(rename_dict)

# 保证纬度升序，便于 interp
if gaez.lat.values[0] > gaez.lat.values[-1]:
    gaez = gaez.sortby("lat")

# 将无效值处理为 NaN（可按实际 tif 修改）
gaez = gaez.where(np.isfinite(gaez))
gaez = gaez.where(gaez > 0)

print("GAEZ:")
print(gaez)
print("GAEZ min/max:", float(gaez.min(skipna=True)), float(gaez.max(skipna=True)))

# =========================
# 3) 构造你的模型 2000 年玉米单产
#    这里给你模板：从 DB['hist'] 计算 actual yield
# =========================
# 假设你已经在内存中有 DB
# 如果没有，请先运行你前面的 build_db_all_gcms()

RELATIVE_YIELD_SCALE = 1   # 如果你的 relative_yield 真实是 0-4 尺度，则保留 4.0；若已是0-1则改成1.0
MODEL_YIELD_TO_KGHA = 1.0 # 若 potential_yield 为 t/ha，则乘1000；若本身 kg/ha 则改1.0

crop = "maiz"
year = 2000

rel_rf  = DB["hist"]["relative_yield"][crop]["waterstress_nonIrr"].sel(year=year) / RELATIVE_YIELD_SCALE
rel_irr = DB["hist"]["relative_yield"][crop]["waterstress_Irr"].sel(year=year)    / RELATIVE_YIELD_SCALE

pot_rf  = DB["hist"]["potential_yield"][crop]["nonIrr"] * MODEL_YIELD_TO_KGHA
pot_irr = DB["hist"]["potential_yield"][crop]["Irr"]    * MODEL_YIELD_TO_KGHA

area_rf  = DB["dist"][f"{crop}_nonIrr"]
area_irr = DB["dist"][f"{crop}_Irr"]
area_tot = area_rf.fillna(0) + area_irr.fillna(0)

y_rf  = rel_rf * pot_rf
y_irr = rel_irr * pot_irr

model_yield_2000 = xr.where(
    area_tot > 0,
    (y_rf.fillna(0) * area_rf.fillna(0) + y_irr.fillna(0) * area_irr.fillna(0)) / 1,
    np.nan
)
model_yield_2000 = model_yield_2000.where(model_yield_2000 > 0)/1

print("MODEL:")
print(model_yield_2000)
print("MODEL min/max:", float(model_yield_2000.min(skipna=True)), float(model_yield_2000.max(skipna=True)))

# =========================
# 4) 将 GAEZ 插值到模型网格
# =========================

gaez_on_model = gaez.interp(
    lat=model_yield_2000.lat,
    lon=model_yield_2000.lon,
    method="linear"
)

mask = np.isfinite(gaez_on_model.values) & np.isfinite(model_yield_2000.values)
gaez_v  = gaez_on_model.values[mask]
model_v = model_yield_2000.values[mask]
# =========================
# 5) 统计指标
# =========================
r = np.corrcoef(gaez_v, model_v)[0, 1]
rmse = np.sqrt(np.mean((model_v - gaez_v) ** 2))
bias = np.mean(model_v - gaez_v)

print(f"R    = {r:.3f}")
print(f"RMSE = {rmse:.2f}")
print(f"Bias = {bias:.2f}")

# =========================
# 6) 差值图
# =========================
diff = model_yield_2000 - gaez_on_model

# =========================
# 7) 统一颜色范围
# =========================
vmax = np.nanpercentile(np.r_[gaez_v, model_v], 99)
vmin = 0

dmax = np.nanpercentile(np.abs(diff.values[np.isfinite(diff.values)]), 99)

# =========================
# 8) 作图
# =========================
fig, axes = plt.subplots(2, 2, figsize=(12, 9))

# GAEZ
im0 = axes[0, 0].pcolormesh(
    gaez_on_model.lon, gaez_on_model.lat, gaez_on_model,
    shading="auto", vmin=vmin, vmax=vmax
)
axes[0, 0].set_title("GAEZ maize yield, 2000")
axes[0, 0].set_xlabel("Lon")
axes[0, 0].set_ylabel("Lat")
plt.colorbar(im0, ax=axes[0, 0], fraction=0.046, pad=0.04, label="Yield")

# Model
im1 = axes[0, 1].pcolormesh(
    model_yield_2000.lon, model_yield_2000.lat, model_yield_2000,
    shading="auto", vmin=vmin, vmax=vmax
)
axes[0, 1].set_title("Model maize yield, 2000")
axes[0, 1].set_xlabel("Lon")
axes[0, 1].set_ylabel("Lat")
plt.colorbar(im1, ax=axes[0, 1], fraction=0.046, pad=0.04, label="Yield")

# Difference
im2 = axes[1, 0].pcolormesh(
    diff.lon, diff.lat, diff,
    shading="auto", vmin=-dmax, vmax=dmax, cmap="RdBu_r"
)
axes[1, 0].set_title("Model - GAEZ")
axes[1, 0].set_xlabel("Lon")
axes[1, 0].set_ylabel("Lat")
plt.colorbar(im2, ax=axes[1, 0], fraction=0.046, pad=0.04, label="Yield difference")

# Scatter
axes[1, 1].scatter(gaez_v, model_v, s=5, alpha=0.3)
xymax = np.nanpercentile(np.r_[gaez_v, model_v], 99.5)
axes[1, 1].plot([0, xymax], [0, xymax], "--", color="gray", lw=1)
axes[1, 1].set_xlim(0, xymax)
axes[1, 1].set_ylim(0, xymax)
axes[1, 1].set_xlabel("GAEZ yield")
axes[1, 1].set_ylabel("Model yield")
axes[1, 1].set_title("Grid-cell comparison")
axes[1, 1].text(
    0.03, 0.97,
    f"R = {r:.2f}\nRMSE = {rmse:.1f}\nBias = {bias:.1f}",
    transform=axes[1, 1].transAxes,
    ha="left", va="top",
    bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.7")
)

for ax in axes.ravel():
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

plt.tight_layout()
fig.savefig(os.path.join(out_dir, "gaez2000_vs_model2000_maize.png"), dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(out_dir, "gaez2000_vs_model2000_maize.pdf"), bbox_inches="tight")
plt.close(fig)

print("Done.")

# %% cell 64
# -*- coding: utf-8 -*-
import os
import numpy as np
import xarray as xr
import rioxarray as rxr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

# =========================
# 0) 全局绘图参数
# =========================
mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["font.family"]  = "serif"
mpl.rcParams["font.serif"]   = ["Arial"]
mpl.rcParams["axes.unicode_minus"] = False

# =========================
# 1) 路径
# =========================
gaez_files = {
    "maiz": r"F:\worldpop\mze_2000_prd.tif",
    "soyb": r"F:\worldpop\soy_2000_prd.tif",
    "whea": r"F:\worldpop\whe_2000_prd.tif",
    "rice": r"F:\worldpop\rcw_2000_prd.tif",
}

crop_title = {
    "maiz": "Maize",
    "soyb": "Soybean",
    "whea": "Wheat",
    "rice": "Rice",
}

out_dir = r"D:\AAUDE\paper_v2\paper4\gaez_production_compare_2000_nointerp"
os.makedirs(out_dir, exist_ok=True)

# =========================
# 2) 参数
# =========================
RELATIVE_YIELD_SCALE = 1  # 若relative_yield已是0-1，改成1.0
MODEL_YIELD_TO_KGHA = 1e-4    # 若potential_yield已是kg/ha，改成1.0
YEAR = 2000

# =========================
# 3) 读取 GAEZ prd 栅格
#    单位：1000 t / cell
# =========================
def open_gaez_prd_tif(fp):
    da = rxr.open_rasterio(fp).squeeze(drop=True)

    rename_dict = {}
    if "x" in da.dims:
        rename_dict["x"] = "lon"
    if "y" in da.dims:
        rename_dict["y"] = "lat"
    if rename_dict:
        da = da.rename(rename_dict)

    # 排序为 lat 从北到南，便于和模型网格对应
    # 你的模型 lat 是 89.75 -> -89.75，故这里也统一成降序
    if da.lat.values[0] < da.lat.values[-1]:
        da = da.sortby("lat", ascending=False)

    # 0 是 fillvalue
    da = da.where(da != 0)
    da = da.where(np.isfinite(da))
    da = da.where(da > 0)

    da.name = "gaez_prod_1000t"
    return da

# =========================
# 4) 将 GAEZ 0.0833° 产量聚合到 0.5°
#    production 必须 sum，不能 interp
# =========================
def aggregate_gaez_prd_to_model_grid(gaez_da, model_template):
    """
    gaez_da: 0.0833° GAEZ production, 单位 1000 t / cell
    model_template: 0.5° model DataArray, 提供目标lat/lon

    返回:
        gaez_0p5(lat, lon): 聚合后的 0.5° production, 单位 1000 t / cell
    """
    # 先检查维度比例
    nlat_f, nlon_f = gaez_da.sizes["lat"], gaez_da.sizes["lon"]
    nlat_c, nlon_c = model_template.sizes["lat"], model_template.sizes["lon"]

    if nlat_f % nlat_c != 0 or nlon_f % nlon_c != 0:
        raise ValueError(
            f"Fine grid ({nlat_f},{nlon_f}) cannot be evenly aggregated "
            f"to coarse grid ({nlat_c},{nlon_c})."
        )

    fy = nlat_f // nlat_c
    fx = nlon_f // nlon_c

    print(f"Aggregation factor: lat={fy}, lon={fx}")

    # coarsen sum：对 production 做守恒聚合
    gaez_0p5 = gaez_da.coarsen(lat=fy, lon=fx, boundary="trim").sum(skipna=True)

    # 直接赋成模型中心坐标，确保一一对应
    gaez_0p5 = gaez_0p5.assign_coords(
        lat=model_template.lat.values,
        lon=model_template.lon.values
    )

    return gaez_0p5

# =========================
# 5) 从 DB 计算模型 2000 年实际产量
#    单位：1000 t / 0.5° cell
# =========================
def build_model_production_2000(DB, crop, year,
                                rel_scale=4.0,
                                yield_unit_factor=1000.0):
    # relative yield
    rel_rf  = DB["hist"]["relative_yield"][crop]["waterstress_nonIrr"].sel(year=year) / rel_scale
    rel_irr = DB["hist"]["relative_yield"][crop]["waterstress_Irr"].sel(year=year)    / rel_scale

    # potential yield -> kg/ha
    pot_rf  = DB["hist"]["potential_yield"][crop]["nonIrr"] * yield_unit_factor
    pot_irr = DB["hist"]["potential_yield"][crop]["Irr"]    * yield_unit_factor

    # harvested area -> ha
    area_rf  = DB["dist"][f"{crop}_nonIrr"]
    area_irr = DB["dist"][f"{crop}_Irr"]

    # actual yield -> kg/ha
    y_rf  = rel_rf * pot_rf
    y_irr = rel_irr * pot_irr

    # production -> kg
    prod_rf_kg  = y_rf.fillna(0)  * area_rf.fillna(0)
    prod_irr_kg = y_irr.fillna(0) * area_irr.fillna(0)

    # sum -> 1000 t
    prod_1000t = (prod_rf_kg + prod_irr_kg) / 1e6
    prod_1000t = prod_1000t.where(prod_1000t > 0)
    prod_1000t.name = "model_prod_1000t"

    return prod_1000t

# =========================
# 6) 指标函数
# =========================
def calc_metrics(obs, sim):
    mask = (
        np.isfinite(obs) &
        np.isfinite(sim) &
        (obs > 0) &
        (sim > 0)
    )
    x = obs[mask]
    y = sim[mask]

    if len(x) < 3:
        return {"N": len(x), "R": np.nan, "RMSE": np.nan, "Bias": np.nan}

    if np.std(x) == 0 or np.std(y) == 0:
        r = np.nan
    else:
        r = np.corrcoef(x, y)[0, 1]

    rmse = np.sqrt(np.mean((y - x) ** 2))
    bias = np.mean(y - x)

    return {"N": len(x), "R": r, "RMSE": rmse, "Bias": bias}

# =========================
# 7) 主循环
# =========================
summary = []

fig, axes = plt.subplots(4, 2, figsize=(12, 18))
plt.subplots_adjust(hspace=0.28, wspace=0.12)

for row, crop in enumerate(["maiz", "soyb", "whea", "rice"]):
    print(f"Processing {crop} ...")

    # ---- model 0.5°
    model = build_model_production_2000(
        DB=DB,
        crop=crop,
        year=YEAR,
        rel_scale=RELATIVE_YIELD_SCALE,
        yield_unit_factor=MODEL_YIELD_TO_KGHA
    )

    # 保证模型纬度方向是降序（通常已经是）
    if model.lat.values[0] < model.lat.values[-1]:
        model = model.sortby("lat", ascending=False)

    # ---- GAEZ 0.0833°
    gaez_fine = open_gaez_prd_tif(gaez_files[crop])

    # ---- 聚合到模型 0.5° 网格
    gaez_0p5 = aggregate_gaez_prd_to_model_grid(gaez_fine, model)
    gaez_0p5=gaez_0p5.where(gaez_0p5>0)/1e3
    # ---- 共同有效像元
    mask = (
        np.isfinite(gaez_0p5.values) &
        np.isfinite(model.values) &
        (gaez_0p5.values > 0) &
        (model.values > 0)
    )
    gaez_v = gaez_0p5.values[mask]
    model_v = model.values[mask]

    # ---- 总产量
    gaez_total_1000t = float(np.nansum(gaez_0p5.values[np.isfinite(gaez_0p5.values)]))
    model_total_1000t = float(np.nansum(model.values[np.isfinite(model.values)]))

    gaez_total_Mt = gaez_total_1000t/1000
    model_total_Mt = model_total_1000t / 1000.0

    # ---- 指标
    met = calc_metrics(gaez_0p5.values, model.values)

    summary.append({
        "crop": crop,
        "gaez_total_1000t": gaez_total_1000t,
        "model_total_1000t": model_total_1000t,
        "gaez_total_Mt": gaez_total_Mt,
        "model_total_Mt": model_total_Mt,
        "N_valid_cells": met["N"],
        "R": met["R"],
        "RMSE_1000t": met["RMSE"],
        "Bias_1000t": met["Bias"],
    })

    # ---- 色标：该作物左右图统一
    vv = np.r_[gaez_v, model_v]
    vmax = np.nanpercentile(vv, 99)
    vmin = 0

    # ---- 左：GAEZ
    ax = axes[row, 0]
    im = ax.pcolormesh(
        gaez_0p5.lon, gaez_0p5.lat, gaez_0p5,
        shading="auto", vmin=vmin, vmax=vmax
    )
    ax.set_title(
        f"{crop_title[crop]} - GAEZ",
        fontsize=11
    )
    ax.set_xlabel("Lon")
    ax.set_ylabel("Lat")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("Production (1000 t)")

    # ---- 右：Model
    ax = axes[row, 1]
    im2 = ax.pcolormesh(
        model.lon, model.lat, model,
        shading="auto", vmin=vmin, vmax=vmax
    )
    ax.set_title(
        f"{crop_title[crop]} - Model\n R = {met['R']:.2f} | RMSE = {met['RMSE']:.2f}",
        fontsize=11
    )
    ax.set_xlabel("Lon")
    ax.set_ylabel("Lat")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    cbar2 = plt.colorbar(im2, ax=ax, fraction=0.046, pad=0.03)
    cbar2.set_label("Production (Mt)")

fig.suptitle("GAEZ v4 vs Model crop production in 2000 (sum aggregation, no interpolation)", fontsize=15, y=0.995)

png_fp = os.path.join(out_dir, "gaez_vs_model_production_2000_nointerp.png")
pdf_fp = os.path.join(out_dir, "gaez_vs_model_production_2000_nointerp.pdf")
# fig.savefig(png_fp, dpi=300, bbox_inches="tight")
fig.savefig(pdf_fp, bbox_inches="tight")
plt.show()
plt.close(fig)

summary_df = pd.DataFrame(summary)
csv_fp = os.path.join(out_dir, "gaez_vs_model_production_2000_nointerp_summary.csv")
summary_df.to_csv(csv_fp, index=False)

print("Done.")
print(summary_df)
print(f"Figure saved to:\n{png_fp}\n{pdf_fp}")
print(f"Summary saved to:\n{csv_fp}")

# %% cell 65
gaez_0p5=gaez_0p5.where(gaez_0p5>0)/1e3
gaez_0p5.plot()

# %% cell 66
import os
import numpy as np
import xarray as xr

import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"]  = 42
mpl.rcParams["axes.unicode_minus"] = False


# ==========================================================
# 1) 工具：面积加权全球均值
# ==========================================================
def area_weighted_global_mean(da, area, eps=1e-12):
    """
    da:   DataArray(year, lat, lon) 或 (lat, lon)
    area: DataArray(lat, lon) (单位任意，只要一致；这里用ha更直观)
    return: DataArray(year) 或 scalar
    """
    da2, a2 = xr.align(da, area, join="inner")
    num = (da2 * a2).sum(dim=("lat", "lon"), skipna=True)
    den = (a2.where(np.isfinite(da2))).sum(dim=("lat", "lon"), skipna=True)
    return num / (den + eps)


def crop_stress_days_global_ts(case_dict, mirca_ds, crop, haz, eps=1e-12):
    """
    case_dict: DB["hist"] 或 DB["future"][gcm][scen]
    crop: "maiz"/"soyb"/"rice"/"whea"
    haz : "dry"/"wet"
    return: DataArray(year) —— 作物总体（Irr+nonIrr）面积加权全球平均 stress days
    """
    # --- 面积(ha)
    A_rf = (mirca_ds[f"{crop}_nonIrr"] / 10000.0).fillna(0)
    A_ir = (mirca_ds[f"{crop}_Irr"]    / 10000.0).fillna(0)

    # --- stress days
    d_rf = case_dict["stress_days"][crop][f"{haz}_nonIrr"]
    d_ir = case_dict["stress_days"][crop][f"{haz}_Irr"]

    # --- 分模式面积加权全球均值
    ts_rf = area_weighted_global_mean(d_rf, A_rf, eps=eps)
    ts_ir = area_weighted_global_mean(d_ir, A_ir, eps=eps)

    # --- 合并（按有效面积权重，避免某些格点缺失导致权重错）
    w_rf = area_weighted_global_mean(xr.ones_like(d_rf), A_rf, eps=eps)
    w_ir = area_weighted_global_mean(xr.ones_like(d_ir), A_ir, eps=eps)
    ts = (ts_rf * w_rf + ts_ir * w_ir) / (w_rf + w_ir + eps)

    return ts


# ==========================================================
# 2) 锚定偏差校正：future(2020) == hist(2019)
# ==========================================================
def anchor_future_to_hist(ts_fut, ts_hist, hist_last_year=2019, fut_first_year=2020):
    """
    加性平移：ts_fut_bc = ts_fut + (hist(2019) - fut(2020))
    """
    # 若没有 year 维直接返回
    if "year" not in ts_fut.dims or "year" not in ts_hist.dims:
        return ts_fut

    if (hist_last_year not in ts_hist["year"]) or (fut_first_year not in ts_fut["year"]):
        # 年份不在坐标里，直接返回
        return ts_fut

    y_hist = ts_hist.sel(year=hist_last_year)
    y_fut0 = ts_fut.sel(year=fut_first_year)

    # 锚点缺失则不校正
    if np.isnan(y_hist.values) or np.isnan(y_fut0.values):
        return ts_fut

    delta = (y_hist - y_fut0)
    return ts_fut + delta


def future_ens_stats_anchored(
    DB, mirca_ds, crop, haz, scen, gcms,
    hist_last_year=2019, fut_first_year=2020,
    qlo=0.05, qhi=0.95
):
    """
    对每个 GCM：先算未来全球序列 -> 锚定到历史末年 -> 再算 ensemble mean/quantile
    return: mean(year), lo(year), hi(year)
    """
    ts_hist = crop_stress_days_global_ts(DB["hist"], mirca_ds, crop=crop, haz=haz)

    ts_list = []
    for gcm in gcms:
        case = DB["future"][gcm][scen]
        ts_fut = crop_stress_days_global_ts(case, mirca_ds, crop=crop, haz=haz)

        ts_fut_bc = anchor_future_to_hist(
            ts_fut, ts_hist,
            hist_last_year=hist_last_year,
            fut_first_year=fut_first_year
        )
        ts_list.append(ts_fut_bc)

    ens = xr.concat(ts_list, dim=xr.IndexVariable("gcm", gcms))  # (gcm, year)
    mean = ens.mean("gcm", skipna=True)
    lo   = ens.quantile(qlo, dim="gcm", skipna=True)
    hi   = ens.quantile(qhi, dim="gcm", skipna=True)
    return mean, lo, hi


# ==========================================================
# 3) 主绘图：2×2 四作物，双haz，双情景
# ==========================================================
def plot_stress_days_2x2_anchored(
    DB, gcms,
    scen_list=("ssp126", "ssp585"),
    crops=("maiz","soyb","rice","whea"),
    hist_year=slice(1980, 2019),
    fut_year=slice(2020, 2100),
    hist_last_year=2019,
    fut_first_year=2020,
    qlo=0.05, qhi=0.95,
    out_fp=None
):
    mirca_ds = DB["dist"]

    crop_title = {"maiz":"Maize", "soyb":"Soybean", "rice":"Rice", "whea":"Wheat"}
    haz_label  = {"dry":"Drought", "wet":"Waterlogging"}
    scen_label = {"ssp126":"SSP1-2.6", "ssp585":"SSP5-8.5"}
    ls_haz     = {"dry":"-", "wet":"--"}  # drought实线，waterlogging虚线

    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.5), sharex=True)
    axes = axes.ravel()

    # 为了图例顺序稳定：手动收集handles
    legend_handles = []
    legend_labels  = []

    for ax, crop in zip(axes, crops):
        # ---- 历史（灰色）
        for haz in ("dry", "wet"):
            ts_h = crop_stress_days_global_ts(DB["hist"], mirca_ds, crop=crop, haz=haz)
            if "year" in ts_h.dims:
                ts_h = ts_h.sel(year=hist_year)

            line_h, = ax.plot(
                ts_h["year"], ts_h,
                ls_haz[haz],
                color="0.35", lw=2.0
            )

            # 只在第一个子图添加到图例
            if crop == crops[0]:
                lab = f"Historical - {haz_label[haz]}"
                legend_handles.append(line_h)
                legend_labels.append(lab)

        # ---- 未来：两情景（锚定校正后 ensemble mean + band）
        for scen in scen_list:
            for haz in ("dry", "wet"):
                mean, lo, hi = future_ens_stats_anchored(
                    DB, mirca_ds,
                    crop=crop, haz=haz, scen=scen, gcms=gcms,
                    hist_last_year=hist_last_year,
                    fut_first_year=fut_first_year,
                    qlo=qlo, qhi=qhi
                )

                mean = mean.sel(year=fut_year)
                lo   = lo.sel(year=fut_year)
                hi   = hi.sel(year=fut_year)

                line_f, = ax.plot(
                    mean["year"], mean,
                    ls_haz[haz], lw=2.0
                )
                ax.fill_between(
                    mean["year"].values,
                    lo.values, hi.values,
                    alpha=0.18
                )

                if crop == crops[0]:
                    lab = f"{scen_label[scen]} - {haz_label[haz]}"
                    legend_handles.append(line_f)
                    legend_labels.append(lab)

        ax.set_title(crop_title.get(crop, crop), fontsize=12)
        ax.set_ylabel("Stress days (days)")
        ax.grid(False)

    for ax in axes[-2:]:
        ax.set_xlabel("Year")

    # 统一图例
    fig.legend(
        legend_handles, legend_labels,
        loc="lower center", ncol=3,
        frameon=False, bbox_to_anchor=(0.5, -0.02)
    )

    fig.suptitle(
        f"Stress days (anchored: future {fut_first_year} = hist {hist_last_year})",
        fontsize=13
    )

    plt.tight_layout()

    if out_fp is not None:
        os.makedirs(os.path.dirname(out_fp), exist_ok=True)
        fig.savefig(out_fp, dpi=300, bbox_inches="tight")
        print("Saved:", out_fp)

    plt.show()


# ==========================================================
# 4) 运行示例
# ==========================================================
out_fp = r"D:\AAUDE\paper\paper4\revision_v2\fig\manu_fig\stress_days_2x2_anchored_hist2019_fut2020.pdf"

plot_stress_days_2x2_anchored(
    DB=DB,
    gcms=gcms,
    scen_list=("ssp126", "ssp585"),
    crops=("maiz","soyb","rice","whea"),
    hist_year=slice(1980, 2019),
    fut_year=slice(2020, 2100),
    hist_last_year=2019,
    fut_first_year=2020,
    qlo=0.05, qhi=0.95,
    out_fp=out_fp
)

# %% cell 67

