"""
Marswalk Planner - terrain preprocessing
Turns NASA/USGS Jezero rasters into the layers + cost grid used by the web app.

Inputs (USGS Astrogeology, Mars 2020 mission products):
  - M20_JezeroCrater_CTXDEM_20m.tif   (MRO CTX elevation, 20 m/px)
  - JEZ_ctx_B_soc_008_orthoMosaic_6m_Eqc_latTs0_lon0.tif (MRO CTX image mosaic, 6 m/px)

Key idea (coordinates): both files use an Equirectangular ("eqc") projection on a
Mars sphere (R = 3,396,190 m). In eqc, x and y are just longitude and latitude
multiplied by constants, so a regular grid in eqc is ALSO a regular grid in
lon/lat degrees. We resample everything onto one lon/lat-aligned grid, which
Leaflet can display directly with its EPSG:4326 (plain degrees) mode.

Usage: python preprocess.py <dem.tif> <ortho.tif> <out_dir>
"""
import sys, json
import numpy as np
import rasterio
from rasterio.warp import reproject, Resampling
from rasterio.transform import from_origin
from rasterio.crs import CRS
from scipy import ndimage
from PIL import Image

DEM, ORTHO, OUT = sys.argv[1], sys.argv[2], sys.argv[3]

R = 3396190.0                      # Mars radius used by USGS (m)
LAT_TS = 18.4663                   # true-scale latitude (Jezero landing site)
# Area of interest (degrees East / North): delta front + Perseverance's route
AOI = dict(west=77.20, east=77.48, south=18.38, north=18.54)
CELL = 20.0                        # analysis cell size (m)
IMG_CELL = 5.0                     # display image cell size (m)

# Local eqc CRS: metres are true-scale near Jezero (east-west not stretched)
LOCAL = CRS.from_proj4(f"+proj=eqc +lat_ts={LAT_TS} +lat_0=0 +lon_0=0 +R={R} +units=m +no_defs")
kx = R * np.cos(np.radians(LAT_TS)) * np.pi / 180   # metres per degree longitude
ky = R * np.pi / 180                                 # metres per degree latitude
x0, x1 = AOI['west'] * kx, AOI['east'] * kx
y0, y1 = AOI['south'] * ky, AOI['north'] * ky


def warp(src_path, cell, resampling, dtype):
    w, h = int(round((x1 - x0) / cell)), int(round((y1 - y0) / cell))
    tr = from_origin(x0, y1, (x1 - x0) / w, (y1 - y0) / h)
    out = np.full((h, w), np.nan, dtype="float32")
    with rasterio.open(src_path) as s:
        reproject(rasterio.band(s, 1), out, src_nodata=s.nodata, dst_nodata=np.nan,
                  dst_transform=tr, dst_crs=LOCAL, resampling=resampling)
    return out


# ---------- 1. Elevation ----------
z = warp(DEM, CELL, Resampling.bilinear, "float32")
h, w = z.shape
if np.isnan(z).any():   # fill tiny gaps with nearest valid value
    idx = ndimage.distance_transform_edt(np.isnan(z), return_distances=False, return_indices=True)
    z = z[tuple(idx)]

# ---------- 2. Slope (degrees) ----------
dzdy, dzdx = np.gradient(z, CELL)
slope = np.degrees(np.arctan(np.hypot(dzdx, dzdy)))

# ---------- 3. Roughness: local relief (std dev of elevation in a 3x3 = 60 m window,
# after removing the smooth trend). High values = bumpy, broken ground.
detr = z - ndimage.uniform_filter(z, 7)
rough = np.sqrt(ndimage.uniform_filter(detr ** 2, 3))

# ---------- 4. Hazard classes ----------
# 0 = easy (<8 deg), 1 = moderate (8-15), 2 = difficult (15-25), 3 = no-go (>25 or very rough)
# Thresholds are design assumptions for a suited astronaut (see README).
haz = np.zeros(z.shape, np.uint8)
haz[slope >= 8] = 1
haz[slope >= 15] = 2
haz[(slope >= 25)] = 3
r_hi = np.percentile(rough, 98)
haz[(rough > r_hi) & (haz < 2)] = 2

# ---------- 5. Images for the map ----------
def save_rgba(arr, path):
    Image.fromarray(arr).save(path, optimize=True)

# hillshade (sun from NW, 35 deg high)
az, alt = np.radians(315), np.radians(35)
aspect = np.arctan2(-dzdx, dzdy)
hs = np.sin(alt) * np.cos(np.radians(slope)) + np.cos(alt) * np.sin(np.radians(slope)) * np.cos(az - aspect)
hs = np.clip(hs * 255, 0, 255).astype(np.uint8)
save_rgba(hs, f"{OUT}/hillshade.png")

# elevation colour ramp (low = deep blue/purple, high = warm orange)
ramp = np.array([[20, 24, 82], [58, 76, 160], [70, 160, 170], [190, 200, 120], [230, 150, 70], [250, 230, 210]], float)
t = (z - z.min()) / (z.max() - z.min()) * (len(ramp) - 1)
i = np.clip(t.astype(int), 0, len(ramp) - 2); f = (t - i)[..., None]
elev_rgb = (ramp[i] * (1 - f) + ramp[i + 1] * f)
elev_rgb = (elev_rgb * (0.55 + 0.45 * hs[..., None] / 255)).astype(np.uint8)
save_rgba(elev_rgb, f"{OUT}/elevation.png")

# hazard (transparent where easy)
pal = np.array([[0, 0, 0, 0], [255, 200, 40, 110], [255, 110, 30, 160], [230, 30, 50, 200]], np.uint8)
save_rgba(pal[haz], f"{OUT}/hazard.png")

# slope colour (continuous)
s = np.clip(slope / 30, 0, 1)
slope_rgba = np.zeros(z.shape + (4,), np.uint8)
slope_rgba[..., 0] = (255 * np.clip(s * 2, 0, 1)).astype(np.uint8)
slope_rgba[..., 1] = (255 * np.clip(2 - s * 2, 0, 1) * 0.85).astype(np.uint8)
slope_rgba[..., 2] = 40
slope_rgba[..., 3] = (60 + 170 * s).astype(np.uint8)
save_rgba(slope_rgba, f"{OUT}/slope.png")

# CTX orthoimage (contrast stretch)
img = warp(ORTHO, IMG_CELL, Resampling.bilinear, "float32")
valid = ~np.isnan(img)
lo, hi = np.percentile(img[valid], [1, 99.5])
img8 = np.zeros(img.shape, np.uint8)
img8[valid] = np.clip((img[valid] - lo) / (hi - lo) * 255, 1, 255)
Image.fromarray(img8).save(f"{OUT}/ctx.jpg", quality=82)

# ---------- 6. Grid for in-browser routing ----------
# elevation (decimetres, int16 would overflow -> store float32) + slope + hazard
z.astype("<f4").tofile(f"{OUT}/elev.f32")
np.round(slope * 10).astype("<u2").tofile(f"{OUT}/slope.u16")   # tenths of a degree
haz.tofile(f"{OUT}/haz.u8")

meta = dict(width=w, height=h, cell_m=CELL, bounds=AOI,
            elev_min=float(z.min()), elev_max=float(z.max()),
            hazard_counts=[int((haz == k).sum()) for k in range(4)],
            crs_note=f"Mars sphere R={R} m, equirectangular; image pixels are linear in lon/lat")
json.dump(meta, open(f"{OUT}/meta.json", "w"), indent=1)
print(json.dumps(meta, indent=1))
print("slope pctiles 50/90/99:", np.percentile(slope, [50, 90, 99]).round(1))
