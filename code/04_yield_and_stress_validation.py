# Public-release script split from the original analysis notebook.
# Paths remain project-specific and may need to be edited before rerunning.
# Part 4: yield validation, GAEZ comparison, and stress-day time-series plotting.

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

