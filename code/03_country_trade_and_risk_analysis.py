# Public-release script split from the original analysis notebook.
# Paths remain project-specific and may need to be edited before rerunning.
# Part 3: country aggregation, FAOSTAT trade links, and risk plots/tables.

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
