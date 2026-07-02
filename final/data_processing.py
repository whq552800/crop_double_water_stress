# Extracted from D:\AAUDE\paper\paper4\revision_v2\code\data_porcessing.ipynb
# Generated for read-only public-release audit. Review before reuse.

# %% [markdown]
# ## de-trend

# %% [markdown]
# ### historical

# %% cell 3
import xarray as xr
import numpy as np
import dask
from dask.distributed import Client, LocalCluster
import os
import warnings
warnings.filterwarnings(
    "ignore",
    message="The specified chunks separate the stored chunks"
)
crop='soyb'
cluster = LocalCluster(
    n_workers=20,
    threads_per_worker=1,
    memory_limit="10GB",
)
client = Client(cluster)
print("Dashboard:", client.dashboard_link)

# 2. 读数据
chunks = {"time": -1, "lat": 40, "lon": 40}

kcc = {
    "yield_1": 0.2,
    "yield_2": 0.8,
    "yield_3": 0.9,
    "yield_4": 0.21,
}

kc = xr.open_dataset(
    rf"F:\paper4\GWSP3_W5E5_v3\currentKY_daily_{crop}.nc",
    chunks=chunks,
).currentKY[5*365+1:]

w = xr.open_dataset(
    r"F:\paper4\GWSP3_W5E5_v3\w_daily.nc",
    chunks=chunks,
).w[5*365+1:]

spei = xr.open_dataset(
    r"E:\dsf\SPEI_daily1.nc",
    engine="h5netcdf",
    chunks=chunks,
).spei.transpose("time", "lat", "lon")

# %% cell 4
import xarray as xr
import numpy as np
import dask
from dask.distributed import Client, LocalCluster
import os

crop='soyb'
cluster = LocalCluster(
    n_workers=20,
    threads_per_worker=1,
    memory_limit="10GB",
)
client = Client(cluster)
print("Dashboard:", client.dashboard_link)

# 2. 读数据
chunks = {"time": -1, "lat": 40, "lon": 40}

kcc = {
    "yield_1": 0.2,
    "yield_2": 0.8,
    "yield_3": 0.9,
    "yield_4": 0.21,
}

kc = xr.open_dataset(
    rf"F:\paper4\GWSP3_W5E5_v3\currentKY_daily_{crop}.nc",
    chunks=chunks,
).currentKY[5*365+1:]

w = xr.open_dataset(
    r"F:\paper4\GWSP3_W5E5_v3\w_daily.nc",
    chunks=chunks,
).w[5*365+1:]

spei = xr.open_dataset(
    r"E:\dsf\SPEI_daily1.nc",
    engine="h5netcdf",
    chunks=chunks,
).spei.transpose("time", "lat", "lon")

seasons = ["yield_1", "yield_2", "yield_3", "yield_4"]

for arr in ['nonIrr','Irr']:
    yield_loss = xr.open_dataset(
        rf"F:\paper4\GWSP3_W5E5_v3\ratio_et_{arr}_daily_{crop}.nc",
        chunks=chunks,
    )[f'ratio_et_{arr}'][5*365+1:]


    for season in seasons:
        print("Prepare:", season)
    
        # 按 KC 分季
        y_season = yield_loss.where(kc == kcc[season])
        y_season=y_season.where(y_season!=0)
        mask = y_season.notnull()
        w1 = w.where(mask)
        sp = spei.where(mask)
    
        nor     = y_season.where((sp > 0.5) & (sp < 1.5))
        flood   = y_season.where((sp > 1.5) & (w1 > 1))
        drought = y_season.where(sp < -1.0)
    
        # 压成 float32，文件小一点
        nor     = nor.astype("float32")
        flood   = flood.astype("float32")
        drought = drought.astype("float32")
    
        ds_out = xr.Dataset(
            {"normal": nor, "flood": flood, "drought": drought}
        )
    
        out_path = rf"E:\paper4\classification\{season}_{crop}_{arr}.nc"
    
        # 如果之前有同名文件，先删掉（防止半成品）
        if os.path.exists(out_path):
            os.remove(out_path)
    
        encoding = {
            "normal":  {"zlib": True, "complevel": 4, "dtype": "float32"},
            "flood":   {"zlib": True, "complevel": 4, "dtype": "float32"},
            "drought": {"zlib": True, "complevel": 4, "dtype": "float32"},
        }
    
        delayed = ds_out.to_netcdf(
            out_path,
            engine="h5netcdf",   # 关键：不用 netCDF4
            encoding=encoding,
            compute=True,       
        )
    
    print("All done.")

# %% [markdown]
# ### future

# %% cell 6
import xarray as xr
import numpy as np
import dask
from dask.distributed import Client, LocalCluster
import os
import warnings
warnings.filterwarnings(
    "ignore",
    message="The specified chunks separate the stored chunks"
)

cluster = LocalCluster(
    n_workers=24,
    threads_per_worker=1,
    memory_limit="100GB",
)
client = Client(cluster)
print("Dashboard:", client.dashboard_link)

# %% cell 7
# -*- coding: utf-8 -*-
import os
from os.path import join
import gc

import xarray as xr
from tqdm.notebook import tqdm

# =========================
# 0) 配置：作物分季 Kc 码
# =========================
KCC = {
    "soyb": {"yield_1": 0.2, "yield_2": 0.8,  "yield_3": 0.9, "yield_4": 0.21},
    "maiz": {"yield_1": 0.4, "yield_2": 0.41, "yield_3": 1.0, "yield_4": 0.5},
    "rice": {"yield_1": 1.0, "yield_2": 1.1,  "yield_3": 1.2, "yield_4": 0.2},
    "wwh":  {"yield_1": 0.2, "yield_2": 0.6,  "yield_3": 0.5, "yield_4": 0.21},
    "swh":  {"yield_1": 0.2, "yield_2": 0.6,  "yield_3": 0.5, "yield_4": 0.21},
}
seasons = ["yield_1", "yield_2", "yield_3", "yield_4"]

# =========================
# 1) Chunk：你当前偏好大块（写盘更友好）
# =========================
chunks = {"time": -1, "lat": 60, "lon": 120}

# 根目录
fp = r"F:\paper4\cwatm_crop_yield"
orig_root = r"F:\paper4\cwatm_crop_yield_original"

# =========================
# 2) 主循环：合并 4 个 season => 每个 crop+arr 只写 1 个文件
# =========================
gcms = sorted([d for d in os.listdir(fp) if os.path.isdir(join(fp, d))])

for gcm in tqdm(gcms, desc="GCM", position=0):

    for ssp in tqdm(["ssp126", "ssp585"], desc="SSP", position=1, leave=False):

        # ---- w / spei 只依赖 gcm/ssp：在 crop 外复用
        w_path = rf"{orig_root}\{gcm}\{ssp}\w_daily.nc"
        sp_path = rf"{orig_root}\{gcm}\{ssp}\spei_daily.nc"

        with xr.open_dataset(w_path, chunks=chunks) as ds:
            w = ds["w"][-11322:].transpose("time", "lat", "lon")

        with xr.open_dataset(sp_path, chunks=chunks, engine="h5netcdf") as ds:
            spei = ds["spei"][-11322:].transpose("time", "lat", "lon")

        for crop, kcc in tqdm(KCC.items(), desc=f"{gcm}-{ssp} crop", position=2, leave=False):

            kc_path = rf"{orig_root}\{gcm}\{ssp}\currentKY_daily_{crop}.nc"
            with xr.open_dataset(kc_path, chunks=chunks) as ds:
                kc = ds["currentKY"][-11322:].transpose("time", "lat", "lon")

            # 作物生长区
            grow = kc > 0
            sp_g = spei.where(grow)
            w_g  = w.where(grow)

            # 气象分类（懒执行；后面多次复用）
            cond_nor     = (sp_g > 0.5) & (sp_g < 1.5)
            cond_flood   = (sp_g > 1.5) & (w_g > 1)
            cond_drought = (sp_g < -1.0)

            for arr in ["nonIrr", "Irr"]:

                y_path = rf"{fp}\{gcm}\{ssp}\ratio_{crop}_{arr}.nc"
                with xr.open_dataset(y_path, chunks=chunks) as ds:
                    # original + smooth -> version
                    y = xr.concat(
                        [ds["original"], ds["smooth"]],
                        dim=xr.IndexVariable("version", ["original", "smooth"])
                    )[-11322:]

                # 提前去掉 0（避免每个 season 都 where 一遍）
                y = y.where(y != 0)

                # ---- 逐季计算：先存到列表，最后 concat 一次写盘
                season_ds_list = []
                for season in tqdm(seasons, desc=f"{crop}-{arr} seasons", position=3, leave=False):

                    # 按 KC 分季（这里用 == 依赖你的 kc 是离散码；若是连续 Kc 曲线需改成容差）
                    y_season = y.where(kc == kcc[season])

                    # mask：original/smooth 任一有值
                    mask = y_season.notnull().any("version")

                    nor     = y_season.where(cond_nor & mask).astype("float32")
                    flood   = y_season.where(cond_flood & mask).astype("float32")
                    drought = y_season.where(cond_drought & mask).astype("float32")

                    ds_out = xr.Dataset({"normal": nor, "flood": flood, "drought": drought})
                    ds_out = ds_out.expand_dims(season=[season])  # 加 season 维度
                    season_ds_list.append(ds_out)

                ds_all = xr.concat(season_ds_list, dim="season")

                # ---- 输出：每个 crop+arr 只写 1 个文件（极大减少写盘次数）
                out_dir = rf"{fp}\{gcm}\{ssp}\event_class"
                os.makedirs(out_dir, exist_ok=True)
                out_path = join(out_dir, f"{crop}_{arr}_allseasons.nc")

                if os.path.exists(out_path):
                    os.remove(out_path)

                encoding = {
                    "normal":  {"zlib": True, "complevel": 6, "dtype": "float32"},
                    "flood":   {"zlib": True, "complevel": 6, "dtype": "float32"},
                    "drought": {"zlib": True, "complevel": 6, "dtype": "float32"},
                }

                ds_all.to_netcdf(
                    out_path,
                    engine="h5netcdf",
                    encoding=encoding,
                    compute=True,
                )

                # 主动释放引用（长循环更稳）
                del y, season_ds_list, ds_all
                gc.collect()

        # 每个 ssp 结束后也清一下
        del w, spei
        gc.collect()

# %% cell 8
open_kwargs = dict(chunks=chunks, engine="h5netcdf", cache=False)

w    = xr.open_dataset(w_path, **open_kwargs)["w"][-11322:].transpose("time","lat","lon")
spei = xr.open_dataset(sp_path, **open_kwargs)["spei"][-11322:].transpose("time","lat","lon")
kc   = xr.open_dataset(kc_path, **open_kwargs)["currentKY"][-11322:].transpose("time","lat","lon")

ds_y = xr.open_dataset(y_path, **open_kwargs)
y = xr.concat([ds_y["original"], ds_y["smooth"]],
              dim=xr.IndexVariable("version", ["original","smooth"]))[-11322:]

# %% cell 9
os.listdir(fp)[8:]

# %% cell 10
import os
from os.path import join
import xarray as xr
from tqdm.notebook import tqdm

KCC = {
    "soyb": {"yield_1": 0.2, "yield_2": 0.8,  "yield_3": 0.9, "yield_4": 0.21},
    "maiz": {"yield_1": 0.4, "yield_2": 0.41, "yield_3": 1.0, "yield_4": 0.5},
    "rice": {"yield_1": 1.0, "yield_2": 1.1,  "yield_3": 1.2, "yield_4": 0.2},
    "wwh":  {"yield_1": 0.2, "yield_2": 0.6,  "yield_3": 0.5, "yield_4": 0.21},
    "swh":  {"yield_1": 0.2, "yield_2": 0.6,  "yield_3": 0.5, "yield_4": 0.21},
}
seasons = ["yield_1", "yield_2", "yield_3", "yield_4"]

chunks = {"time": -1, "lat": 60, "lon": 90}

fp = r"F:\paper4\cwatm_crop_yield"
open_kwargs = dict(chunks=chunks, engine="h5netcdf", cache=False)
for gcm in tqdm(os.listdir(fp)[8:], desc="GCM", position=0):
    if gcm=='MRI-ESM2-0':
        continue
    for ssp in tqdm(["ssp126", "ssp585"],desc="SSP", position=1, leave=False):
        if ssp=='ssp126':
            continue
        # ssp='ssp585'
        # w / spei 只依赖 gcm/ssp
        w = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\w_daily.nc",
            **open_kwargs
        ).w[-11322:].transpose("time", "lat", "lon")

        spei = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\spei_daily.nc",
            **open_kwargs
        ).spei.transpose("time", "lat", "lon")[-11322:]

        for crop, kcc in tqdm(KCC.items(), desc=f"{gcm}-{ssp} crop", position=2, leave=False):

            kc = xr.open_dataset(
                rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\currentKY_daily_{crop}.nc",
                **open_kwargs
            ).currentKY[-11322:].transpose("time", "lat", "lon")

            grow = kc > 0
            sp_g = spei.where(grow)

            cond_nor     = (sp_g > 0.5) & (sp_g < 1.5)
            cond_flood   = (sp_g > 1.5) & (w.where(grow) > 1)
            cond_drought = (sp_g < -1.0)

            for arr in ["nonIrr", "Irr"]:
                # arr="Irr"
                ds_y = xr.open_dataset(
                    rf"F:\paper4\cwatm_crop_yield\{gcm}\{ssp}\ratio_{crop}_{arr}.nc",
                    **open_kwargs
                )

                y = xr.concat(
                    [ds_y["original"], ds_y["smooth"]],
                    dim=xr.IndexVariable("version", ["original", "smooth"])
                )[-11322:]

                for season in tqdm(
                    seasons,
                    desc=f"{crop}-{arr}",
                    position=3,
                    leave=False,
                ):
                    y_season = y.where(kc == kcc[season]).where(y != 0)
                    mask = y_season.notnull().any("version")

                    nor     = y_season.where(cond_nor & mask).astype("float32")
                    flood   = y_season.where(cond_flood & mask).astype("float32")
                    drought = y_season.where(cond_drought & mask).astype("float32")

                    ds_out = xr.Dataset(
                        {"normal": nor, "flood": flood, "drought": drought}
                    )

                    out_dir = rf"F:\paper4\cwatm_crop_yield\{gcm}\{ssp}"
                    os.makedirs(out_dir, exist_ok=True)
                    out_path = join(out_dir, f"{season}_{crop}_{arr}.nc")

                    if os.path.exists(out_path):
                        os.remove(out_path)

                    encoding = {
                        "normal":  {"zlib": True, "complevel": 4, "dtype": "float32"},
                        "flood":   {"zlib": True, "complevel": 4, "dtype": "float32"},
                        "drought": {"zlib": True, "complevel": 4, "dtype": "float32"},
                    }

                    ds_out.to_netcdf(
                        out_path,
                        engine="h5netcdf",
                        encoding=encoding,
                        compute=True,
                    )

# %% cell 11
import os
from os.path import join
import xarray as xr
from tqdm.notebook import tqdm

KCC = {
    "soyb": {"yield_1": 0.2, "yield_2": 0.8,  "yield_3": 0.9, "yield_4": 0.21},
    "maiz": {"yield_1": 0.4, "yield_2": 0.41, "yield_3": 1.0, "yield_4": 0.5},
    "rice": {"yield_1": 1.0, "yield_2": 1.1,  "yield_3": 1.2, "yield_4": 0.2},
    "wwh":  {"yield_1": 0.2, "yield_2": 0.6,  "yield_3": 0.5, "yield_4": 0.21},
    "swh":  {"yield_1": 0.2, "yield_2": 0.6,  "yield_3": 0.5, "yield_4": 0.21},
}
seasons = ["yield_1", "yield_2", "yield_3", "yield_4"]

chunks = {"time": -1, "lat": 60, "lon": 90}

fp = r"F:\paper4\cwatm_crop_yield"
open_kwargs = dict(chunks=chunks, engine="h5netcdf", cache=False)
for gcm in tqdm(os.listdir(fp)[5:], desc="GCM", position=0):
    ssp='ssp126'
    w = xr.open_dataset(
        rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\w_daily.nc",
        **open_kwargs
    ).w[-11322:].transpose("time", "lat", "lon")

    spei = xr.open_dataset(
        rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\spei_daily.nc",
        **open_kwargs
    ).spei.transpose("time", "lat", "lon")[-11322:]

    for crop, kcc in tqdm(KCC.items(), desc=f"{gcm}-{ssp} crop", position=1, leave=False):

        kc = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\currentKY_daily_{crop}.nc",
            **open_kwargs
        ).currentKY[-11322:].transpose("time", "lat", "lon")

        grow = kc > 0
        sp_g = spei.where(grow)

        cond_nor     = (sp_g > 0.5) & (sp_g < 1.5)
        cond_flood   = (sp_g > 1.5) & (w.where(grow) > 1)
        cond_drought = (sp_g < -1.0)

        for arr in tqdm(["nonIrr", "Irr"],position=2,
                leave=False,desc='irrigation'
            ):
            # arr="Irr"
            ds_y = xr.open_dataset(
                rf"F:\paper4\cwatm_crop_yield\{gcm}\{ssp}\ratio_{crop}_{arr}.nc",
                **open_kwargs
            )

            y = xr.concat(
                [ds_y["original"], ds_y["smooth"]],
                dim=xr.IndexVariable("version", ["original", "smooth"])
            )[-11322:]

            for season in tqdm(
                seasons,
                desc=f"{crop}-{arr}",
                position=3,
                leave=False,
            ):
                y_season = y.where(kc == kcc[season]).where(y != 0)
                mask = y_season.notnull().any("version")

                nor     = y_season.where(cond_nor & mask).astype("float32")
                flood   = y_season.where(cond_flood & mask).astype("float32")
                drought = y_season.where(cond_drought & mask).astype("float32")

                ds_out = xr.Dataset(
                    {"normal": nor, "flood": flood, "drought": drought}
                )

                out_dir = rf"F:\paper4\cwatm_crop_yield\{gcm}\{ssp}"
                os.makedirs(out_dir, exist_ok=True)
                out_path = join(out_dir, f"{season}_{crop}_{arr}.nc")

                if os.path.exists(out_path):
                    os.remove(out_path)

                encoding = {
                    "normal":  {"zlib": True, "complevel": 4, "dtype": "float32"},
                    "flood":   {"zlib": True, "complevel": 4, "dtype": "float32"},
                    "drought": {"zlib": True, "complevel": 4, "dtype": "float32"},
                }

                ds_out.to_netcdf(
                    out_path,
                    engine="h5netcdf",
                    encoding=encoding,
                    compute=True,
                )

# %% cell 12
import os
from os.path import join
import xarray as xr
from tqdm.notebook import tqdm

KCC = {
    "soyb": {"yield_1": 0.2, "yield_2": 0.8,  "yield_3": 0.9, "yield_4": 0.21},
    "maiz": {"yield_1": 0.4, "yield_2": 0.41, "yield_3": 1.0, "yield_4": 0.5},
    "rice": {"yield_1": 1.0, "yield_2": 1.1,  "yield_3": 1.2, "yield_4": 0.2},
    "wwh":  {"yield_1": 0.2, "yield_2": 0.6,  "yield_3": 0.5, "yield_4": 0.21},
    "swh":  {"yield_1": 0.2, "yield_2": 0.6,  "yield_3": 0.5, "yield_4": 0.21},
}
seasons = ["yield_1", "yield_2", "yield_3", "yield_4"]

chunks = {"time": -1, "lat": 60, "lon": 60}

fp = r"F:\paper4\cwatm_crop_yield"
open_kwargs = dict(chunks=chunks, engine="h5netcdf", cache=False)

gcm=os.listdir(fp)[-1]
ssp="ssp585"

    # ssp='ssp585'
    # w / spei 只依赖 gcm/ssp
w = xr.open_dataset(
    rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\w_daily.nc",
    **open_kwargs
).w[-11322:].transpose("time", "lat", "lon")

spei = xr.open_dataset(
    rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\spei_daily.nc",
    **open_kwargs
).spei.transpose("time", "lat", "lon")[-11322:]

for crop, kcc in tqdm(KCC.items(), desc=f"{gcm}-{ssp} crop", position=0, leave=False):

    kc = xr.open_dataset(
        rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\currentKY_daily_{crop}.nc",
        **open_kwargs
    ).currentKY[-11322:].transpose("time", "lat", "lon")

    grow = kc > 0
    sp_g = spei.where(grow)

    cond_nor     = (sp_g > 0.5) & (sp_g < 1.5)
    cond_flood   = (sp_g > 1.5) & (w.where(grow) > 1)
    cond_drought = (sp_g < -1.0)
    arr="Irr"
    for arr in ["nonIrr", "Irr"]:
        # arr="Irr"
        ds_y = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield\{gcm}\{ssp}\ratio_{crop}_{arr}.nc",
            **open_kwargs
        )
    
        y = xr.concat(
            [ds_y["original"], ds_y["smooth"]],
            dim=xr.IndexVariable("version", ["original", "smooth"])
        )[-11322:]
    
        for season in tqdm(
            seasons[1:],
            desc=f"{crop}-{arr}",
            position=1,
            leave=False,
        ):
            y_season = y.where(kc == kcc[season]).where(y != 0)
            mask = y_season.notnull().any("version")
    
            nor     = y_season.where(cond_nor & mask).astype("float32")
            flood   = y_season.where(cond_flood & mask).astype("float32")
            drought = y_season.where(cond_drought & mask).astype("float32")
    
            ds_out = xr.Dataset(
                {"normal": nor, "flood": flood, "drought": drought}
            )
    
            out_dir = rf"F:\paper4\cwatm_crop_yield\{gcm}\{ssp}"
            os.makedirs(out_dir, exist_ok=True)
            out_path = join(out_dir, f"{season}_{crop}_{arr}.nc")
    
            if os.path.exists(out_path):
                os.remove(out_path)
    
            encoding = {
                "normal":  {"zlib": True, "complevel": 4, "dtype": "float32"},
                "flood":   {"zlib": True, "complevel": 4, "dtype": "float32"},
                "drought": {"zlib": True, "complevel": 4, "dtype": "float32"},
            }
    
            ds_out.to_netcdf(
                out_path,
                engine="h5netcdf",
                encoding=encoding,
                compute=True,
            )

# %% cell 13
pwd

# %% cell 14
import xarray as xr
import os
from os.path import join

crop='soyb'
kcc = {
    "yield_1": 0.2,
    "yield_2": 0.8,
    "yield_3": 0.9,
    "yield_4": 0.21,
}
chunks = {"time": 1000, "lat": 40, "lon": 40}
seasons = ["yield_1", "yield_2", "yield_3", "yield_4"]

fp = r'F:\paper4\cwatm_crop_yield'

for gcm in os.listdir(fp):
    for ssp in ['ssp126','ssp585']:

        kc = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\currentKY_daily_{crop}.nc",
            chunks=chunks,
        ).currentKY[-11322:]

        w = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\w_daily.nc",
            chunks=chunks,
        ).w[-11322:]

        spei = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\spei_daily.nc",
            engine="h5netcdf",
            chunks=chunks,
        ).spei.transpose("time", "lat", "lon")[-11322:]
        for arr in ['nonIrr','Irr']:

            # ✅ 读取 original + smooth，并拼成 version 维度
            ds_y = xr.open_dataset(
                rf"F:\paper4\cwatm_crop_yield\{gcm}\{ssp}\ratio_{crop}_{arr}.nc",
                chunks=chunks,
            )

            # 假设变量名就叫 original 和 smooth（按你描述）
            y = xr.concat(
                [ds_y["original"], ds_y["smooth"]],
                dim=xr.IndexVariable("version", ["original", "smooth"])
            )[-11322:]
            cond_nor     = (spei > 0.5) & (spei < 1.5)
            cond_flood   = (spei > 1.5) & (w > 1)
            cond_drought = (spei < -1.0)
            for season in seasons:
                print("Prepare:", gcm, ssp, arr, season)

                # 按 KC 分季（此时 y 多了 version 维度）
                y_season = y.where(kc == kcc[season])
                y_season = y_season.where(y_season != 0)

                # mask 用 original/smooth 任一有值都保留（更稳）
                mask = y_season.notnull().any("version")

                # w1 = w.where(mask)
                # sp = spei.where(mask)

                

                nor     = y_season.where(cond_nor& mask).astype("float32")
                flood   = y_season.where(cond_flood& mask).astype("float32")
                drought = y_season.where(cond_drought& mask).astype("float32")

                ds_out = xr.Dataset({"normal": nor, "flood": flood, "drought": drought})

                out_dir = rf"F:\paper4\cwatm_crop_yield\{gcm}\{ssp}"
                os.makedirs(out_dir, exist_ok=True)
                out_path = join(out_dir, f"{season}_{crop}_{arr}.nc")

                if os.path.exists(out_path):
                    os.remove(out_path)

                encoding = {
                    "normal":  {"zlib": True, "complevel": 4, "dtype": "float32"},
                    "flood":   {"zlib": True, "complevel": 4, "dtype": "float32"},
                    "drought": {"zlib": True, "complevel": 4, "dtype": "float32"},
                }

                ds_out.to_netcdf(
                    out_path,
                    engine="h5netcdf",
                    encoding=encoding,
                    compute=True,
                )

            print("All done:", gcm, ssp, arr)

# %% cell 15
import xarray as xr
import os
from os.path import join

crop='maiz'
kcc = {
    "yield_1": 0.4,
    "yield_2": 0.41,
    "yield_3": 1,
    "yield_4": 0.5,
}
chunks = {"time": -1, "lat": 40, "lon": 40}
seasons = ["yield_1", "yield_2", "yield_3", "yield_4"]

fp = r'F:\paper4\cwatm_crop_yield'

for gcm in os.listdir(fp):
    for ssp in ['ssp126','ssp585']:

        kc = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\currentKY_daily_{crop}.nc",
            chunks=chunks,
        ).currentKY[-11322:]

        w = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\w_daily.nc",
            chunks=chunks,
        ).w[-11322:]

        spei = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\spei_daily.nc",
            engine="h5netcdf",
            chunks=chunks,
        ).spei.transpose("time", "lat", "lon")[-11322:]


        for arr in ['nonIrr','Irr']:

            # ✅ 读取 original + smooth，并拼成 version 维度
            ds_y = xr.open_dataset(
                rf"F:\paper4\cwatm_crop_yield\{gcm}\{ssp}\ratio_{crop}_{arr}.nc",
                chunks=chunks,
            )

            # 假设变量名就叫 original 和 smooth（按你描述）
            y = xr.concat(
                [ds_y["original"], ds_y["smooth"]],
                dim=xr.IndexVariable("version", ["original", "smooth"])
            )[-11322:]

            for season in seasons:
                print("Prepare:", gcm, ssp, arr, season)

                # 按 KC 分季（此时 y 多了 version 维度）
                y_season = y.where(kc == kcc[season])
                y_season = y_season.where(y_season != 0)

                # mask 用 original/smooth 任一有值都保留（更稳）
                mask = y_season.notnull().any("version")

                w1 = w.where(mask)
                sp = spei.where(mask)

                cond_nor     = (sp > 0.5) & (sp < 1.5)
                cond_flood   = (sp > 1.5) & (w1 > 1)
                cond_drought = (sp < -1.0)

                nor     = y_season.where(cond_nor).astype("float32")
                flood   = y_season.where(cond_flood).astype("float32")
                drought = y_season.where(cond_drought).astype("float32")

                ds_out = xr.Dataset({"normal": nor, "flood": flood, "drought": drought})

                out_dir = rf"F:\paper4\cwatm_crop_yield\{gcm}\{ssp}"
                os.makedirs(out_dir, exist_ok=True)
                out_path = join(out_dir, f"{season}_{crop}_{arr}.nc")

                if os.path.exists(out_path):
                    os.remove(out_path)

                encoding = {
                    "normal":  {"zlib": True, "complevel": 4, "dtype": "float32"},
                    "flood":   {"zlib": True, "complevel": 4, "dtype": "float32"},
                    "drought": {"zlib": True, "complevel": 4, "dtype": "float32"},
                }

                ds_out.to_netcdf(
                    out_path,
                    engine="h5netcdf",
                    encoding=encoding,
                    compute=True,
                )

            print("All done:", gcm, ssp, arr)

# %% cell 16
import xarray as xr
import os
from os.path import join

crop='rice'
kcc = {
    "yield_1": 1,
    "yield_2": 1.1,
    "yield_3": 1.2,
    "yield_4": 0.2,
}
chunks = {"time": -1, "lat": 40, "lon": 40}
seasons = ["yield_1", "yield_2", "yield_3", "yield_4"]

fp = r'F:\paper4\cwatm_crop_yield'

for gcm in os.listdir(fp):
    for ssp in ['ssp126','ssp585']:

        kc = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\currentKY_daily_{crop}.nc",
            chunks=chunks,
        ).currentKY[-11322:]

        w = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\w_daily.nc",
            chunks=chunks,
        ).w[-11322:]

        spei = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\spei_daily.nc",
            engine="h5netcdf",
            chunks=chunks,
        ).spei.transpose("time", "lat", "lon")[-11322:]


        for arr in ['nonIrr','Irr']:

            # ✅ 读取 original + smooth，并拼成 version 维度
            ds_y = xr.open_dataset(
                rf"F:\paper4\cwatm_crop_yield\{gcm}\{ssp}\ratio_{crop}_{arr}.nc",
                chunks=chunks,
            )

            # 假设变量名就叫 original 和 smooth（按你描述）
            y = xr.concat(
                [ds_y["original"], ds_y["smooth"]],
                dim=xr.IndexVariable("version", ["original", "smooth"])
            )[-11322:]

            for season in seasons:
                print("Prepare:", gcm, ssp, arr, season)

                # 按 KC 分季（此时 y 多了 version 维度）
                y_season = y.where(kc == kcc[season])
                y_season = y_season.where(y_season != 0)

                # mask 用 original/smooth 任一有值都保留（更稳）
                mask = y_season.notnull().any("version")

                w1 = w.where(mask)
                sp = spei.where(mask)

                cond_nor     = (sp > 0.5) & (sp < 1.5)
                cond_flood   = (sp > 1.5) & (w1 > 1)
                cond_drought = (sp < -1.0)

                nor     = y_season.where(cond_nor).astype("float32")
                flood   = y_season.where(cond_flood).astype("float32")
                drought = y_season.where(cond_drought).astype("float32")

                ds_out = xr.Dataset({"normal": nor, "flood": flood, "drought": drought})

                out_dir = rf"F:\paper4\cwatm_crop_yield\{gcm}\{ssp}"
                os.makedirs(out_dir, exist_ok=True)
                out_path = join(out_dir, f"{season}_{crop}_{arr}.nc")

                if os.path.exists(out_path):
                    os.remove(out_path)

                encoding = {
                    "normal":  {"zlib": True, "complevel": 4, "dtype": "float32"},
                    "flood":   {"zlib": True, "complevel": 4, "dtype": "float32"},
                    "drought": {"zlib": True, "complevel": 4, "dtype": "float32"},
                }

                ds_out.to_netcdf(
                    out_path,
                    engine="h5netcdf",
                    encoding=encoding,
                    compute=True,
                )

            print("All done:", gcm, ssp, arr)

# %% cell 17
import xarray as xr
import os
from os.path import join

crop='wwh'
kcc = {
    "yield_1": 0.2,
    "yield_2": 0.6,
    "yield_3": 0.5,
    "yield_4": 0.21,
}
chunks = {"time": -1, "lat": 40, "lon": 40}
seasons = ["yield_1", "yield_2", "yield_3", "yield_4"]

fp = r'F:\paper4\cwatm_crop_yield'

for gcm in os.listdir(fp):
    for ssp in ['ssp126','ssp585']:

        kc = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\currentKY_daily_{crop}.nc",
            chunks=chunks,
        ).currentKY[-11322:]

        w = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\w_daily.nc",
            chunks=chunks,
        ).w[-11322:]

        spei = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\spei_daily.nc",
            engine="h5netcdf",
            chunks=chunks,
        ).spei.transpose("time", "lat", "lon")[-11322:]


        for arr in ['nonIrr','Irr']:

            # ✅ 读取 original + smooth，并拼成 version 维度
            ds_y = xr.open_dataset(
                rf"F:\paper4\cwatm_crop_yield\{gcm}\{ssp}\ratio_{crop}_{arr}.nc",
                chunks=chunks,
            )

            # 假设变量名就叫 original 和 smooth（按你描述）
            y = xr.concat(
                [ds_y["original"], ds_y["smooth"]],
                dim=xr.IndexVariable("version", ["original", "smooth"])
            )[-11322:]

            for season in seasons:
                print("Prepare:", gcm, ssp, arr, season)

                # 按 KC 分季（此时 y 多了 version 维度）
                y_season = y.where(kc == kcc[season])
                y_season = y_season.where(y_season != 0)

                # mask 用 original/smooth 任一有值都保留（更稳）
                mask = y_season.notnull().any("version")

                w1 = w.where(mask)
                sp = spei.where(mask)

                cond_nor     = (sp > 0.5) & (sp < 1.5)
                cond_flood   = (sp > 1.5) & (w1 > 1)
                cond_drought = (sp < -1.0)

                nor     = y_season.where(cond_nor).astype("float32")
                flood   = y_season.where(cond_flood).astype("float32")
                drought = y_season.where(cond_drought).astype("float32")

                ds_out = xr.Dataset({"normal": nor, "flood": flood, "drought": drought})

                out_dir = rf"F:\paper4\cwatm_crop_yield\{gcm}\{ssp}"
                os.makedirs(out_dir, exist_ok=True)
                out_path = join(out_dir, f"{season}_{crop}_{arr}.nc")

                if os.path.exists(out_path):
                    os.remove(out_path)

                encoding = {
                    "normal":  {"zlib": True, "complevel": 4, "dtype": "float32"},
                    "flood":   {"zlib": True, "complevel": 4, "dtype": "float32"},
                    "drought": {"zlib": True, "complevel": 4, "dtype": "float32"},
                }

                ds_out.to_netcdf(
                    out_path,
                    engine="h5netcdf",
                    encoding=encoding,
                    compute=True,
                )

            print("All done:", gcm, ssp, arr)

# %% cell 18
import xarray as xr
import os
from os.path import join

crop='swh'
kcc = {
    "yield_1": 0.2,
    "yield_2": 0.6,
    "yield_3": 0.5,
    "yield_4": 0.21,
}
chunks = {"time": -1, "lat": 40, "lon": 40}
seasons = ["yield_1", "yield_2", "yield_3", "yield_4"]

fp = r'F:\paper4\cwatm_crop_yield'

for gcm in os.listdir(fp):
    for ssp in ['ssp126','ssp585']:

        kc = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\currentKY_daily_{crop}.nc",
            chunks=chunks,
        ).currentKY[-11322:]

        w = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\w_daily.nc",
            chunks=chunks,
        ).w[-11322:]

        spei = xr.open_dataset(
            rf"F:\paper4\cwatm_crop_yield_original\{gcm}\{ssp}\spei_daily.nc",
            engine="h5netcdf",
            chunks=chunks,
        ).spei.transpose("time", "lat", "lon")[-11322:]


        for arr in ['nonIrr','Irr']:

            # ✅ 读取 original + smooth，并拼成 version 维度
            ds_y = xr.open_dataset(
                rf"F:\paper4\cwatm_crop_yield\{gcm}\{ssp}\ratio_{crop}_{arr}.nc",
                chunks=chunks,
            )

            # 假设变量名就叫 original 和 smooth（按你描述）
            y = xr.concat(
                [ds_y["original"], ds_y["smooth"]],
                dim=xr.IndexVariable("version", ["original", "smooth"])
            )[-11322:]

            for season in seasons:
                print("Prepare:", gcm, ssp, arr, season)

                # 按 KC 分季（此时 y 多了 version 维度）
                y_season = y.where(kc == kcc[season])
                y_season = y_season.where(y_season != 0)

                # mask 用 original/smooth 任一有值都保留（更稳）
                mask = y_season.notnull().any("version")

                w1 = w.where(mask)
                sp = spei.where(mask)

                cond_nor     = (sp > 0.5) & (sp < 1.5)
                cond_flood   = (sp > 1.5) & (w1 > 1)
                cond_drought = (sp < -1.0)

                nor     = y_season.where(cond_nor).astype("float32")
                flood   = y_season.where(cond_flood).astype("float32")
                drought = y_season.where(cond_drought).astype("float32")

                ds_out = xr.Dataset({"normal": nor, "flood": flood, "drought": drought})

                out_dir = rf"F:\paper4\cwatm_crop_yield\{gcm}\{ssp}"
                os.makedirs(out_dir, exist_ok=True)
                out_path = join(out_dir, f"{season}_{crop}_{arr}.nc")

                if os.path.exists(out_path):
                    os.remove(out_path)

                encoding = {
                    "normal":  {"zlib": True, "complevel": 4, "dtype": "float32"},
                    "flood":   {"zlib": True, "complevel": 4, "dtype": "float32"},
                    "drought": {"zlib": True, "complevel": 4, "dtype": "float32"},
                }

                ds_out.to_netcdf(
                    out_path,
                    engine="h5netcdf",
                    encoding=encoding,
                    compute=True,
                )

            print("All done:", gcm, ssp, arr)

# %% cell 19

