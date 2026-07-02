# Public-release script split from the original analysis notebook.
# Paths remain project-specific and may need to be edited before rerunning.
# Part 2: historical/future spatial loss maps and latitudinal summaries.

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

