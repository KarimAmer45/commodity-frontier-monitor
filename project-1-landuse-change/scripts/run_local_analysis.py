"""Anonymous local mirror of the Earth Engine land-use-change analysis.

Downloads only the needed windows from public Cloud-Optimized GeoTIFFs,
creates analysis rasters, calculates area statistics, and renders the map and
report. Earth Search supplies Sentinel-2 L2A; the official Hansen v1.13 COG
supplies loss year. No cloud account is required.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from contextlib import ExitStack
from datetime import date
from pathlib import Path
from xml.sax.saxutils import escape

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import rasterio
import fsspec
import tifffile
import zarr
from affine import Affine
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from pyproj import Transformer
from pystac_client import Client
from rasterio.enums import Resampling
from rasterio.transform import from_origin
from rasterio.vrt import WarpedVRT
from rasterio.warp import reproject
from rasterio.windows import Window
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.platypus import (
    Image,
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "processed"
MAP_DIR = ROOT / "map"
OUTPUT_PDF_DIR = ROOT / "output" / "pdf"
TMP_DIR = ROOT / "tmp"

BBOX = (101.30, 0.30, 101.65, 0.65)
AOI_NAME = "Tapung palm-oil frontier, Riau, Indonesia"
TARGET_CRS = "EPSG:32647"
RESOLUTION = 10
HANSEN_RESOLUTION = 30
MAX_SCENES = 8
CHUNK_ROWS = 256
CLOUD_LIMIT = 40
PERIODS = {
    "before": ("2019-01-01", "2019-12-31"),
    "after": ("2024-01-01", "2024-12-31"),
}
HANSEN_URL = (
    "https://storage.googleapis.com/earthenginepartners-hansen/"
    "GFC-2025-v1.13/Hansen_GFC-2025-v1.13_lossyear_10N_100E.tif"
)


class RemoteCog:
    """Read COG tiles through fsspec HTTP ranges instead of GDAL HTTPS.

    This avoids a Windows Schannel issue in some sandboxed sessions while
    preserving true windowed COG access (the full source files are not saved).
    """

    def __init__(self, url: str):
        self.url = url
        self.file = fsspec.open(url, block_size=8 * 1024 * 1024).open()
        self.tiff = tifffile.TiffFile(self.file)
        metadata = self.tiff.geotiff_metadata or {}
        scale = metadata.get("ModelPixelScale")
        tiepoint = metadata.get("ModelTiepoint")
        if not scale or not tiepoint:
            raise ValueError(f"Missing GeoTIFF transform tags: {url}")
        self.transform = Affine(
            float(scale[0]), 0, float(tiepoint[3]),
            0, -float(scale[1]), float(tiepoint[4]),
        )
        self.store = self.tiff.aszarr()
        opened = zarr.open(self.store, mode="r")
        self.array = opened["0"] if hasattr(opened, "array_keys") else opened
        self.height, self.width = self.array.shape[-2:]

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def close(self):
        self.store.close()
        self.tiff.close()
        self.file.close()

    def read_on_aligned_grid(self, target_transform, window: Window):
        """Nearest-neighbour sample onto a grid in the same projected CRS."""
        col0 = int(window.col_off)
        row0 = int(window.row_off)
        width = int(window.width)
        height = int(window.height)
        xs = target_transform.c + (col0 + np.arange(width) + 0.5) * target_transform.a
        ys = target_transform.f + (row0 + np.arange(height) + 0.5) * target_transform.e
        source_cols = np.floor((xs - self.transform.c) / self.transform.a).astype(int)
        source_rows = np.floor((self.transform.f - ys) / abs(self.transform.e)).astype(int)
        valid_cols = (source_cols >= 0) & (source_cols < self.width)
        valid_rows = (source_rows >= 0) & (source_rows < self.height)
        if not valid_cols.any() or not valid_rows.any():
            return np.zeros((height, width), dtype=self.array.dtype)
        cmin, cmax = source_cols[valid_cols].min(), source_cols[valid_cols].max()
        rmin, rmax = source_rows[valid_rows].min(), source_rows[valid_rows].max()
        source = np.asarray(self.array[rmin : rmax + 1, cmin : cmax + 1])
        result = np.zeros((height, width), dtype=source.dtype)
        local_cols = np.clip(source_cols - cmin, 0, source.shape[1] - 1)
        local_rows = np.clip(source_rows - rmin, 0, source.shape[0] - 1)
        sampled = source[local_rows[:, None], local_cols[None, :]]
        result[np.ix_(valid_rows, valid_cols)] = sampled[np.ix_(valid_rows, valid_cols)]
        return result

    def geographic_subset(self, bbox):
        """Return an EPSG:4326 source subset and its affine transform."""
        west, south, east, north = bbox
        col0 = max(0, int(math.floor((west - self.transform.c) / self.transform.a)) - 2)
        col1 = min(self.width, int(math.ceil((east - self.transform.c) / self.transform.a)) + 2)
        row0 = max(0, int(math.floor((self.transform.f - north) / abs(self.transform.e))) - 2)
        row1 = min(self.height, int(math.ceil((self.transform.f - south) / abs(self.transform.e))) + 2)
        subset = np.asarray(self.array[row0:row1, col0:col1])
        subset_transform = self.transform * Affine.translation(col0, row0)
        return subset, subset_transform


def aligned_grid(resolution: int):
    transformer = Transformer.from_crs("EPSG:4326", TARGET_CRS, always_xy=True)
    xs, ys = [], []
    for lon, lat in (
        (BBOX[0], BBOX[1]),
        (BBOX[0], BBOX[3]),
        (BBOX[2], BBOX[1]),
        (BBOX[2], BBOX[3]),
    ):
        x, y = transformer.transform(lon, lat)
        xs.append(x)
        ys.append(y)
    left = math.floor(min(xs) / resolution) * resolution
    right = math.ceil(max(xs) / resolution) * resolution
    bottom = math.floor(min(ys) / resolution) * resolution
    top = math.ceil(max(ys) / resolution) * resolution
    width = int(round((right - left) / resolution))
    height = int(round((top - bottom) / resolution))
    return from_origin(left, top, resolution, resolution), width, height


def query_scenes(start: str, end: str, max_scenes: int):
    catalog = Client.open("https://earth-search.aws.element84.com/v1")
    search = catalog.search(
        collections=["sentinel-2-l2a"],
        bbox=list(BBOX),
        datetime=f"{start}/{end}",
        query={"eo:cloud_cover": {"lt": CLOUD_LIMIT}},
        max_items=250,
    )
    candidates = sorted(
        search.items(), key=lambda item: item.properties.get("eo:cloud_cover", 100)
    )
    chosen = []
    seen_dates = set()
    for item in candidates:
        acquisition_date = item.datetime.date().isoformat()
        if acquisition_date in seen_dates:
            continue
        if not all(key in item.assets for key in ("red", "nir", "scl")):
            continue
        seen_dates.add(acquisition_date)
        chosen.append(item)
        if len(chosen) >= max_scenes:
            break
    if len(chosen) < 3:
        raise RuntimeError(f"Only {len(chosen)} usable scenes found for {start} to {end}")
    return chosen


def output_profile(transform, width, height, dtype, nodata, count=1):
    return {
        "driver": "GTiff",
        "width": width,
        "height": height,
        "count": count,
        "dtype": dtype,
        "crs": TARGET_CRS,
        "transform": transform,
        "nodata": nodata,
        "compress": "DEFLATE",
        "predictor": 3 if dtype.startswith("float") else 2,
        "tiled": True,
        "blockxsize": 256,
        "blockysize": 256,
        "BIGTIFF": "IF_SAFER",
    }


def make_composite(label: str, scenes, out_path: Path):
    transform, width, height = aligned_grid(RESOLUTION)
    profile = output_profile(transform, width, height, "float32", -9999.0)
    valid_scl = np.array([2, 4, 5, 6, 7], dtype=np.uint8)
    with ExitStack() as stack:
        triples = []
        for item in scenes:
            red = stack.enter_context(RemoteCog(item.assets["red"].href))
            nir = stack.enter_context(RemoteCog(item.assets["nir"].href))
            scl = stack.enter_context(RemoteCog(item.assets["scl"].href))
            triples.append((red, nir, scl))

        with rasterio.open(out_path, "w", **profile) as dst:
            for row in range(0, height, CHUNK_ROWS):
                rows = min(CHUNK_ROWS, height - row)
                window = Window(0, row, width, rows)
                values = []
                for red_cog, nir_cog, scl_cog in triples:
                    red = red_cog.read_on_aligned_grid(transform, window).astype("float32")
                    nir = nir_cog.read_on_aligned_grid(transform, window).astype("float32")
                    scl = scl_cog.read_on_aligned_grid(transform, window)
                    denominator = nir + red
                    with np.errstate(divide="ignore", invalid="ignore"):
                        ndvi = (nir - red) / denominator
                    good = np.isin(scl, valid_scl) & (denominator > 0)
                    values.append(np.where(good, ndvi, np.nan).astype("float32"))
                with np.errstate(invalid="ignore"):
                    composite = np.nanmedian(np.stack(values), axis=0)
                dst.write(np.where(np.isfinite(composite), composite, -9999.0), 1, window=window)
    return {
        "label": label,
        "scene_count": len(scenes),
        "scenes": [
            {
                "id": item.id,
                "date": item.datetime.date().isoformat(),
                "cloud_cover": round(item.properties.get("eo:cloud_cover", 0), 2),
            }
            for item in scenes
        ],
    }


def calculate_change(before_path: Path, after_path: Path):
    diff_path = DATA_DIR / "ndvi_diff.tif"
    clearing_path = DATA_DIR / "clearing_mask.tif"
    with rasterio.open(before_path) as before, rasterio.open(after_path) as after:
        diff_profile = before.profile.copy()
        diff_profile.update(dtype="float32", nodata=-9999.0, predictor=3)
        mask_profile = before.profile.copy()
        mask_profile.update(dtype="uint8", nodata=0, predictor=2)
        clearing_pixels = 0
        valid_pixels = 0
        with rasterio.open(diff_path, "w", **diff_profile) as diff_dst, rasterio.open(
            clearing_path, "w", **mask_profile
        ) as clear_dst:
            for _, window in before.block_windows(1):
                b = before.read(1, window=window)
                a = after.read(1, window=window)
                valid = (b != before.nodata) & (a != after.nodata)
                diff = np.where(valid, a - b, -9999.0).astype("float32")
                clearing = (valid & (diff < -0.20) & (b > 0.50)).astype("uint8")
                diff_dst.write(diff, 1, window=window)
                clear_dst.write(clearing, 1, window=window)
                clearing_pixels += int(clearing.sum())
                valid_pixels += int(valid.sum())
    return diff_path, clearing_path, clearing_pixels, valid_pixels


def make_hansen_loss():
    transform30, width30, height30 = aligned_grid(HANSEN_RESOLUTION)
    out_path = DATA_DIR / "forest_loss_2020_2024.tif"
    profile = output_profile(transform30, width30, height30, "uint8", 0)
    loss_pixels = 0
    with RemoteCog(HANSEN_URL) as source:
        year_subset, subset_transform = source.geographic_subset(BBOX)
        years_projected = np.zeros((height30, width30), dtype="uint8")
        reproject(
            source=year_subset,
            destination=years_projected,
            src_transform=subset_transform,
            src_crs="EPSG:4326",
            dst_transform=transform30,
            dst_crs=TARGET_CRS,
            src_nodata=0,
            dst_nodata=0,
            resampling=Resampling.nearest,
        )
        recent_all = ((years_projected >= 20) & (years_projected <= 24)).astype("uint8")
        with rasterio.open(out_path, "w", **profile) as dst:
            for _, window in dst.block_windows(1):
                row0, col0 = int(window.row_off), int(window.col_off)
                rows, cols = int(window.height), int(window.width)
                recent = recent_all[row0 : row0 + rows, col0 : col0 + cols]
                dst.write(recent, 1, window=window)
                loss_pixels += int(recent.sum())
    return out_path, loss_pixels


def calculate_overlap(clearing_path: Path, hansen_path: Path):
    overlap_path = DATA_DIR / "agreement_mask.tif"
    overlap_pixels = 0
    with rasterio.open(clearing_path) as clearing, rasterio.open(hansen_path) as hansen:
        profile = clearing.profile.copy()
        profile.update(dtype="uint8", nodata=0, predictor=2)
        with WarpedVRT(
            hansen,
            crs=clearing.crs,
            transform=clearing.transform,
            width=clearing.width,
            height=clearing.height,
            resampling=Resampling.nearest,
            nodata=0,
        ) as hansen10, rasterio.open(overlap_path, "w", **profile) as dst:
            for _, window in clearing.block_windows(1):
                c = clearing.read(1, window=window) == 1
                h = hansen10.read(1, window=window) == 1
                overlap = (c & h).astype("uint8")
                dst.write(overlap, 1, window=window)
                overlap_pixels += int(overlap.sum())
    return overlap_path, overlap_pixels


def raster_array(path: Path, out_shape=(1200, 1200), resampling=Resampling.average):
    with rasterio.open(path) as src:
        scale = min(out_shape[1] / src.width, out_shape[0] / src.height, 1)
        width = max(1, int(src.width * scale))
        height = max(1, int(src.height * scale))
        arr = src.read(1, out_shape=(height, width), resampling=resampling)
        extent = (src.bounds.left, src.bounds.right, src.bounds.bottom, src.bounds.top)
        nodata = src.nodata
    return arr, extent, nodata


def render_map(results):
    diff, extent, nodata = raster_array(DATA_DIR / "ndvi_diff.tif")
    clearing, _, _ = raster_array(
        DATA_DIR / "clearing_mask.tif", resampling=Resampling.nearest
    )
    hansen, _, _ = raster_array(
        DATA_DIR / "forest_loss_2020_2024.tif", resampling=Resampling.nearest
    )
    # Hansen array has a different grid; read it directly onto the display grid.
    with rasterio.open(DATA_DIR / "ndvi_diff.tif") as target, rasterio.open(
        DATA_DIR / "forest_loss_2020_2024.tif"
    ) as hsrc:
        display_h = diff.shape[0]
        display_w = diff.shape[1]
        display_transform = target.transform * target.transform.scale(
            target.width / display_w, target.height / display_h
        )
        with WarpedVRT(
            hsrc,
            crs=target.crs,
            transform=display_transform,
            width=display_w,
            height=display_h,
            resampling=Resampling.nearest,
            nodata=0,
        ) as hvrt:
            hansen = hvrt.read(1)

    diff = np.ma.masked_where(diff == nodata, diff)
    fig = plt.figure(figsize=(11.69, 8.27), facecolor="#f5f2ea")
    ax = fig.add_axes([0.065, 0.14, 0.72, 0.72], facecolor="#e7e2d7")
    cmap = mcolors.LinearSegmentedColormap.from_list(
        "ndvi_change", ["#9e1b32", "#ef8a62", "#f7f7f7", "#91cf60", "#1b7837"]
    )
    image = ax.imshow(diff, cmap=cmap, vmin=-0.4, vmax=0.4, extent=extent, origin="upper")
    # Orange Hansen overlay and black overlap pixels.
    hmask = np.ma.masked_where(hansen != 1, hansen)
    ax.imshow(
        hmask,
        cmap=mcolors.ListedColormap(["#ff8c00"]),
        alpha=0.72,
        extent=extent,
        origin="upper",
        interpolation="nearest",
    )
    cmask = np.ma.masked_where(clearing != 1, clearing)
    ax.contour(
        cmask,
        levels=[0.5],
        colors=["#d900d9"],
        linewidths=0.42,
        extent=extent,
        origin="upper",
    )
    agreement = (clearing == 1) & (hansen == 1)
    amask = np.ma.masked_where(~agreement, agreement.astype("uint8"))
    ax.imshow(
        amask,
        cmap=mcolors.ListedColormap(["#111111"]),
        extent=extent,
        origin="upper",
        interpolation="nearest",
    )

    ax.set_xlabel("UTM easting (m), WGS 84 / UTM zone 47N", fontsize=8)
    ax.set_ylabel("UTM northing (m)", fontsize=8)
    ax.tick_params(labelsize=7, colors="#3b3b36")
    ax.grid(color="white", linewidth=0.35, alpha=0.45)
    for spine in ax.spines.values():
        spine.set_edgecolor("#504f48")

    # North arrow.
    ax.annotate(
        "N",
        xy=(0.945, 0.94),
        xytext=(0.945, 0.82),
        xycoords="axes fraction",
        ha="center",
        va="center",
        fontsize=11,
        fontweight="bold",
        arrowprops=dict(facecolor="#20201d", edgecolor="#20201d", width=4, headwidth=12),
    )
    # 10 km scale bar.
    x0 = extent[0] + 0.055 * (extent[1] - extent[0])
    y0 = extent[2] + 0.055 * (extent[3] - extent[2])
    ax.plot([x0, x0 + 10000], [y0, y0], color="white", linewidth=6, solid_capstyle="butt")
    ax.plot([x0, x0 + 10000], [y0, y0], color="#1d1d1b", linewidth=2, solid_capstyle="butt")
    ax.text(x0 + 5000, y0 + 800, "10 km", ha="center", fontsize=8, fontweight="bold")

    fig.text(
        0.065,
        0.935,
        "LAND-USE CHANGE OVER THE TAPUNG PALM-OIL FRONTIER",
        fontsize=17,
        fontweight="bold",
        color="#18392b",
    )
    fig.text(0.065, 0.895, "Riau, Indonesia | NDVI 2019-2024 | Hansen loss 2020-2024", fontsize=11, color="#56554e")

    cax = fig.add_axes([0.81, 0.59, 0.025, 0.24])
    cb = fig.colorbar(image, cax=cax)
    cb.set_label("NDVI change", fontsize=9)
    cb.ax.tick_params(labelsize=7)
    fig.text(0.81, 0.85, "Vegetation gain", fontsize=8, color="#1b7837")
    fig.text(0.81, 0.565, "Vegetation loss", fontsize=8, color="#9e1b32")

    legend_ax = fig.add_axes([0.81, 0.31, 0.17, 0.21])
    legend_ax.axis("off")
    legend_ax.text(0, 1.0, "Detected change", fontsize=10, fontweight="bold", va="top")
    handles = [
        Line2D([0], [0], color="#d900d9", lw=1.5, label="NDVI clearing outline"),
        Patch(facecolor="#ff8c00", alpha=0.72, label="Hansen forest loss"),
        Patch(facecolor="#111111", label="Both methods agree"),
    ]
    legend_ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(-0.05, 0.82), frameon=False, fontsize=8)

    stats_ax = fig.add_axes([0.81, 0.13, 0.17, 0.15])
    stats_ax.axis("off")
    stats_ax.add_patch(Rectangle((0, 0), 1, 1, facecolor="#e7ede7", edgecolor="#9cab9e"))
    stats_ax.text(0.07, 0.83, "AREA SUMMARY", fontsize=9, fontweight="bold", color="#18392b")
    stats_ax.text(0.07, 0.61, f"NDVI flagged   {results['ndvi_flagged_ha']:,.0f} ha", fontsize=8)
    stats_ax.text(0.07, 0.42, f"Hansen loss    {results['hansen_loss_ha']:,.0f} ha", fontsize=8)
    stats_ax.text(0.07, 0.23, f"Overlap        {results['agreement_ha']:,.0f} ha", fontsize=8)
    stats_ax.text(0.07, 0.07, f"IoU            {results['iou_pct']:.1f}%", fontsize=8)

    fig.text(
        0.065,
        0.06,
        "AOI: 101.30-101.65 E, 0.30-0.65 N | NDVI flag: change < -0.20 and 2019 NDVI > 0.50 | "
        "Hansen loss years 20-24",
        fontsize=7.5,
        color="#55554f",
    )
    fig.text(
        0.065,
        0.035,
        "Data: Copernicus Sentinel-2 SR Harmonized + Google Cloud Score+ via Earth Engine; Hansen GFC v1.13. "
        "Map: WGS 84 / UTM zone 47N.",
        fontsize=7.2,
        color="#55554f",
    )
    png_path = MAP_DIR / "landuse_change_map.png"
    pdf_path = OUTPUT_PDF_DIR / "landuse_change_map.pdf"
    fig.savefig(png_path, dpi=300, facecolor=fig.get_facecolor())
    fig.savefig(pdf_path, dpi=300, facecolor=fig.get_facecolor())
    plt.close(fig)
    return png_path, pdf_path


def write_report_md(results):
    report = f"""# Land-use change over the Tapung palm-oil frontier, Riau, 2019-2024

## 1. Question and context

Where did dense vegetation decline within a 1,510 km2 screening area in Tapung, Riau, between 2019 and 2024, and how well do Sentinel-2 NDVI flags coincide with independent Hansen forest-loss observations after the baseline year? Riau is a major palm-oil production landscape, so transparent, repeatable screening can support supply-chain due diligence. This is a risk-screening result, not proof that a particular parcel was converted for a particular commodity.

## 2. Data

- **AOI:** 101.30-101.65 E, 0.30-0.65 N; WGS 84 / UTM zone 47N.
- **Sentinel-2:** `COPERNICUS/S2_SR_HARMONIZED`, filtered to scenes with <40% cloud metadata for 2019 and 2024. Google Cloud Score+ (`cs >= 0.60`) supplied the per-pixel clear-sky mask.
- **Forest loss:** Hansen Global Forest Change v1.13 loss year, filtered to codes 20-24 (2020-2024). Excluding 2019 avoids treating disturbances during the annual baseline-composite year as subsequent change.

The final measurements and GeoTIFFs were computed and exported from Google Earth Engine. An anonymous local SCL-based mirror is retained only as an independent reproducibility aid; it is not the source of the reported figures.

## 3. Method

For each period, NDVI = (B8 - B4) / (B8 + B4) was calculated after Cloud Score+ masking. A per-pixel annual median reduced residual cloud and acquisition-specific noise. Change was 2024 median NDVI minus 2019 median NDVI. A likely-clearing flag required NDVI change < -0.20 and 2019 NDVI > 0.50. The first condition targets a substantial decline; the second reduces false flags over already sparse or non-vegetated surfaces.

The NDVI flag was intersected with Hansen loss from 2020-2024. Earth Engine calculated clearing, Hansen, and agreement areas on a common 10 m analysis scale using `ee.Image.pixelArea()`. Agreement is reported as the intersection, shares relative to each method, and intersection over union (IoU). Exported Hansen pixels retain a 30 m cartographic resolution.

## 4. Results

The screening flagged **{results['ndvi_flagged_ha']:,.1f} ha** of likely vegetation clearing. Hansen recorded **{results['hansen_loss_ha']:,.1f} ha** of forest loss for 2020-2024. Their spatial intersection was **{results['agreement_ha']:,.1f} ha**. This equals **{results['ndvi_overlap_pct']:.1f}%** of the NDVI-flagged area and **{results['hansen_overlap_pct']:.1f}%** of Hansen loss; IoU was **{results['iou_pct']:.1f}%**.

![Land-use change map](map/landuse_change_map.png)

These totals are screening indicators. Non-overlap is expected because NDVI detects vegetation changes outside Hansen's >5 m tree-cover definition, while Hansen can capture stand-replacement events whose annual median NDVI signal is muted by regrowth, timing, mixed pixels, or cloud availability.

## 5. Limitations

1. **Attribution:** neither NDVI nor Hansen identifies the responsible commodity, actor, or legal status. Parcel boundaries and field evidence are required.
2. **Phenology and management:** crop rotation, harvesting, fire, flooding, drought, and plantation cycles can resemble clearing.
3. **Threshold sensitivity:** -0.20, 0.50, and Cloud Score+ 0.60 are transparent screening choices, not locally calibrated decision rules.
4. **Cloud and compositing:** cloud scores can miss haze or shadow, and annual medians can suppress short-lived events.
5. **Product mismatch:** Sentinel-2 is evaluated at 10 m while Hansen is Landsat-derived at 30.92 m; alignment creates mixed-pixel and edge effects.
6. **Definition mismatch:** Hansen forest loss is stand-replacement disturbance of vegetation taller than 5 m; NDVI decline is a broader signal.
7. **Cartographic resampling:** the map's black overlap layer is derived from the exported 10 m/30 m rasters and can differ slightly from Earth Engine's common-scale server measurement.

## 6. Next steps

- Run threshold sensitivity and stratified manual validation against dated high-resolution imagery.
- Add NBR, NDMI, and red-edge indices to separate burning, moisture stress, harvest, and clearing.
- Use monthly time series with LandTrendr or CCDC to estimate disturbance timing and persistence.
- Intersect validated change with concession, supplier, cadastral, protected-area, peat, and EUDR-relevant cutoff layers.
- Publish a confusion matrix from independently labelled sample points and report uncertainty intervals.

## Reproducibility note

The authoritative workflow is `gee/landuse_change_riau.js`; Earth Engine Console measurements and export metadata are recorded in `results.json`. `scripts/run_local_analysis.py` provides an anonymous local comparison.

## References

- Google Earth Engine Data Catalog. *Harmonized Sentinel-2 MSI: MultiSpectral Instrument, Level-2A (SR).* https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR_HARMONIZED
- Google Earth Engine Data Catalog. *Cloud Score+ S2_HARMONIZED V1.* https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_CLOUD_SCORE_PLUS_V1_S2_HARMONIZED
- Google Earth Engine Data Catalog. *Hansen Global Forest Change v1.13 (2000-2025).* https://developers.google.com/earth-engine/datasets/catalog/UMD_hansen_global_forest_change_2025_v1_13
- Hansen, M. C. et al. (2013). *High-Resolution Global Maps of 21st-Century Forest Cover Change.* Science 342(6160), 850-853. https://doi.org/10.1126/science.1244693
"""
    path = ROOT / "report.md"
    path.write_text(report, encoding="utf-8")
    return path


def footer(canvas, doc):
    canvas.saveState()
    canvas.setFont("Helvetica", 7.5)
    canvas.setFillColor(colors.HexColor("#6b6b64"))
    canvas.drawString(20 * mm, 11 * mm, "Land-use change screening | Tapung, Riau | 2019-2024")
    canvas.drawRightString(190 * mm, 11 * mm, f"Page {doc.page}")
    canvas.restoreState()


def build_report_pdf(results, map_png: Path):
    path = OUTPUT_PDF_DIR / "landuse_change_report.pdf"
    doc = SimpleDocTemplate(
        str(path),
        pagesize=A4,
        title="Land-use change over the Tapung palm-oil frontier, Riau, 2019-2024",
        author="Karim Amer",
        subject="Sentinel-2 and Hansen GFC land-use change screening portfolio project",
        creator="Karim Amer",
        rightMargin=18 * mm,
        leftMargin=18 * mm,
        topMargin=15 * mm,
        bottomMargin=18 * mm,
    )
    styles = getSampleStyleSheet()
    styles.add(
        ParagraphStyle(
            name="TitleGreen",
            parent=styles["Title"],
            fontName="Helvetica-Bold",
            fontSize=19,
            leading=22,
            textColor=colors.HexColor("#18392b"),
            alignment=TA_CENTER,
            spaceAfter=5 * mm,
        )
    )
    styles.add(
        ParagraphStyle(
            name="SectionGreen",
            parent=styles["Heading2"],
            fontName="Helvetica-Bold",
            fontSize=11.5,
            leading=14,
            textColor=colors.HexColor("#18392b"),
            spaceBefore=3.5 * mm,
            spaceAfter=1.5 * mm,
        )
    )
    body = ParagraphStyle(
        "BodyClean",
        parent=styles["BodyText"],
        fontName="Helvetica",
        fontSize=8.5,
        leading=11.2,
        textColor=colors.HexColor("#292924"),
        spaceAfter=2 * mm,
    )
    small = ParagraphStyle("Small", parent=body, fontSize=7.4, leading=9.2)
    story = [
        Paragraph("Land-use change over the Tapung palm-oil frontier", styles["TitleGreen"]),
        Paragraph("Riau, Indonesia | Sentinel-2 and Hansen GFC | 2019-2024", ParagraphStyle("Sub", parent=body, alignment=TA_CENTER, fontSize=9.5, textColor=colors.HexColor("#5b5b55"))),
        Spacer(1, 2 * mm),
    ]

    metrics = [
        ["NDVI flagged", "Hansen loss", "Spatial overlap", "IoU"],
        [
            f"{results['ndvi_flagged_ha']:,.0f} ha",
            f"{results['hansen_loss_ha']:,.0f} ha",
            f"{results['agreement_ha']:,.0f} ha",
            f"{results['iou_pct']:.1f}%",
        ],
    ]
    table = Table(metrics, colWidths=[43 * mm] * 4, rowHeights=[7 * mm, 10 * mm])
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#18392b")),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("BACKGROUND", (0, 1), (-1, 1), colors.HexColor("#e7ede7")),
                ("TEXTCOLOR", (0, 1), (-1, 1), colors.HexColor("#18392b")),
                ("FONTNAME", (0, 0), (-1, -1), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, 0), 7.5),
                ("FONTSIZE", (0, 1), (-1, 1), 11),
                ("ALIGN", (0, 0), (-1, -1), "CENTER"),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("GRID", (0, 0), (-1, -1), 0.4, colors.HexColor("#9cab9e")),
            ]
        )
    )
    story.extend([table, Spacer(1, 3 * mm)])
    story.extend(
        [
            Paragraph("1. Question and context", styles["SectionGreen"]),
            Paragraph(
                "Where did dense vegetation decline in a 1,510 km2 screening area in Tapung, Riau, between 2019 and 2024, and how well do Sentinel-2 flags coincide with independent Hansen forest-loss observations? Riau is a major palm-oil landscape, making reproducible screening relevant to supply-chain due diligence. This is a risk screen, not parcel-level attribution.",
                body,
            ),
            Paragraph("2. Data and method", styles["SectionGreen"]),
            Paragraph(
                "The AOI is 101.30-101.65 E, 0.30-0.65 N in WGS 84 / UTM zone 47N. Earth Engine Sentinel-2 SR Harmonized imagery was filtered to &lt;40% scene cloud for 2019 and 2024. Google Cloud Score+ (cs &gt;= 0.60) supplied the per-pixel clear-sky mask. Hansen GFC v1.13 loss years 20-24 (2020-2024) supplied the independent comparison after the baseline year.",
                body,
            ),
            Paragraph(
                "NDVI was calculated from B8 and B4. Annual per-pixel medians reduced acquisition noise. Change equals 2024 minus 2019 NDVI; likely clearing required change &lt; -0.20 and 2019 NDVI &gt; 0.50. Earth Engine measured clearing, Hansen loss, and agreement at a common 10 m scale using pixel area; the Hansen export retains 30 m cartographic resolution.",
                body,
            ),
            Image(str(map_png), width=171 * mm, height=121 * mm),
            Paragraph(
                "Figure 1. NDVI change with NDVI-flagged clearing (magenta) and Hansen forest loss (orange).",
                small,
            ),
            PageBreak(),
            Paragraph("3. Results", styles["SectionGreen"]),
            Paragraph(
                f"The screen flagged <b>{results['ndvi_flagged_ha']:,.1f} ha</b> of likely vegetation clearing. Hansen recorded <b>{results['hansen_loss_ha']:,.1f} ha</b> of loss in 2020-2024. Their intersection was <b>{results['agreement_ha']:,.1f} ha</b>, equal to <b>{results['ndvi_overlap_pct']:.1f}%</b> of NDVI flags and <b>{results['hansen_overlap_pct']:.1f}%</b> of Hansen loss; IoU was <b>{results['iou_pct']:.1f}%</b>. Non-overlap is expected because NDVI covers broader vegetation change while Hansen targets stand-replacement disturbance of vegetation taller than 5 m.",
                body,
            ),
            Paragraph("4. Limitations", styles["SectionGreen"]),
        ]
    )
    limitations = [
        "Neither product identifies the responsible commodity, actor, legality, or exact conversion cause.",
        "Crop cycles, harvest, fire, flood, drought, and plantation rotation can resemble clearing.",
        "The -0.20 / 0.50 thresholds are not locally calibrated; sensitivity testing is required.",
        "Cloud Score+ can miss haze or shadow, and annual medians can suppress short-lived events.",
        "The 10 m Sentinel-2 and 30.92 m Landsat-derived grids create edge and mixed-pixel effects.",
        "Hansen's forest definition is narrower than NDVI-based vegetation loss.",
        "The map overlap is derived from exported 10 m/30 m rasters and may differ slightly from Earth Engine's common-scale server measurement.",
    ]
    for item in limitations:
        story.append(Paragraph(f"&#8226;&nbsp;&nbsp;{item}", body))
    story.extend(
        [
            Paragraph("5. Next steps", styles["SectionGreen"]),
            Paragraph(
                "Validate a stratified sample against dated high-resolution imagery; test thresholds; add NBR, NDMI, and red-edge indices; estimate disturbance timing with LandTrendr or CCDC; then intersect validated change with supplier parcels, concessions, protected areas, peat, and EUDR-relevant cutoffs. Operational claims should include a confusion matrix and uncertainty intervals.",
                body,
            ),
            Paragraph("Reproducibility", styles["SectionGreen"]),
            Paragraph(
                "The repository includes the authoritative Earth Engine JavaScript, Earth Engine GeoTIFF exports, Console measurements, an anonymous local comparison, machine-readable metadata, a QGIS-ready project, and PNG/PDF map exports. Data sources: Copernicus Sentinel-2 L2A, Google Cloud Score+, and Hansen GFC v1.13.",
                small,
            ),
            Paragraph("References", styles["SectionGreen"]),
            Paragraph(
                'Google Earth Engine Data Catalog. <link href="https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR_HARMONIZED" color="#245b44">Harmonized Sentinel-2 MSI Level-2A SR</link>.',
                small,
            ),
            Paragraph(
                'Google Earth Engine Data Catalog. <link href="https://developers.google.com/earth-engine/datasets/catalog/GOOGLE_CLOUD_SCORE_PLUS_V1_S2_HARMONIZED" color="#245b44">Cloud Score+ S2_HARMONIZED V1</link>.',
                small,
            ),
            Paragraph(
                'Google Earth Engine Data Catalog. <link href="https://developers.google.com/earth-engine/datasets/catalog/UMD_hansen_global_forest_change_2025_v1_13" color="#245b44">Hansen Global Forest Change v1.13 (2000-2025)</link>.',
                small,
            ),
            Paragraph(
                'Hansen, M. C. et al. (2013). <link href="https://doi.org/10.1126/science.1244693" color="#245b44">High-Resolution Global Maps of 21st-Century Forest Cover Change</link>. Science 342(6160), 850-853.',
                small,
            ),
        ]
    )
    doc.build(story, onFirstPage=footer, onLaterPages=footer)
    return path


def write_qgis_project():
    qgis_dir = ROOT / "qgis"
    qgis_dir.mkdir(parents=True, exist_ok=True)
    layers = [
        ("ndvi_diff", "NDVI change 2019-2024", "../data/processed/ndvi_diff.tif"),
        ("forest_loss", "Hansen forest loss 2020-2024", "../data/processed/forest_loss_2020_2024.tif"),
        ("clearing", "NDVI clearing flag", "../data/processed/clearing_mask.tif"),
        ("agreement", "Agreement", "../data/processed/agreement_mask.tif"),
    ]
    tree = "".join(
        f'<layer-tree-layer checked="Qt::Checked" id="{layer_id}" name="{escape(name)}" providerKey="gdal" source="{source}"/>'
        for layer_id, name, source in layers
    )
    project_layers = []
    for layer_id, name, source in layers:
        if layer_id == "ndvi_diff":
            renderer = """
            <rasterrenderer alphaBand="-1" band="1" classificationMin="-0.4" classificationMax="0.4" type="singlebandpseudocolor" opacity="1">
              <rasterTransparency/>
              <rastershader>
                <colorrampshader colorRampType="INTERPOLATED" classificationMode="1" clip="0">
                  <item alpha="255" value="-0.4" label="Vegetation loss" color="#9e1b32"/>
                  <item alpha="255" value="0" label="No change" color="#f7f7f7"/>
                  <item alpha="255" value="0.4" label="Vegetation gain" color="#1b7837"/>
                </colorrampshader>
              </rastershader>
            </rasterrenderer>"""
        else:
            colors_by_layer = {"forest_loss": "#ff8c00", "clearing": "#d900d9", "agreement": "#111111"}
            renderer = f"""
            <rasterrenderer alphaBand="-1" band="1" type="paletted" opacity="0.8">
              <rasterTransparency><singleValuePixelList><pixelListEntry min="0" max="0" percentTransparent="100"/></singleValuePixelList></rasterTransparency>
              <colorPalette><paletteEntry alpha="255" value="1" label="Detected" color="{colors_by_layer[layer_id]}"/></colorPalette>
            </rasterrenderer>"""
        project_layers.append(
            f"""<maplayer type="raster" hasScaleBasedVisibilityFlag="0" autoRefreshTime="0" autoRefreshMode="Disabled">
              <id>{layer_id}</id><datasource>{source}</datasource><layername>{escape(name)}</layername>
              <srs><spatialrefsys><wkt></wkt><proj4>+proj=utm +zone=47 +datum=WGS84 +units=m +no_defs</proj4><srsid>32647</srsid><srid>32647</srid><authid>EPSG:32647</authid><description>WGS 84 / UTM zone 47N</description><projectionacronym>utm</projectionacronym><ellipsoidacronym>EPSG:7030</ellipsoidacronym><geographicflag>false</geographicflag></spatialrefsys></srs>
              <provider>gdal</provider>{renderer}
            </maplayer>"""
        )
    qgs = f"""<?xml version="1.0" encoding="UTF-8"?>
<qgis projectname="Tapung land-use change" version="3.40.0-Bratislava" saveUserFull="karim" saveUser="karim">
  <homePath path=""/><title>Land-use change over the Tapung palm-oil frontier, Riau, 2019-2024</title>
  <projectCrs><spatialrefsys><proj4>+proj=utm +zone=47 +datum=WGS84 +units=m +no_defs</proj4><srsid>32647</srsid><srid>32647</srid><authid>EPSG:32647</authid><description>WGS 84 / UTM zone 47N</description><projectionacronym>utm</projectionacronym><ellipsoidacronym>EPSG:7030</ellipsoidacronym><geographicflag>false</geographicflag></spatialrefsys></projectCrs>
  <layer-tree-group checked="Qt::Checked" expanded="1" name="">{tree}</layer-tree-group>
  <projectlayers>{''.join(project_layers)}</projectlayers>
</qgis>"""
    path = qgis_dir / "tapung_landuse_change.qgs"
    path.write_text(qgs, encoding="utf-8")
    return path


def run(max_scenes: int):
    for directory in (DATA_DIR, MAP_DIR, OUTPUT_PDF_DIR, TMP_DIR):
        directory.mkdir(parents=True, exist_ok=True)
    transform, width, height = aligned_grid(RESOLUTION)
    metadata = {
        "created": date.today().isoformat(),
        "aoi_name": AOI_NAME,
        "bbox_wgs84": list(BBOX),
        "target_crs": TARGET_CRS,
        "resolution_m": RESOLUTION,
        "cloud_limit_pct": CLOUD_LIMIT,
        "ndvi_change_threshold": -0.20,
        "baseline_ndvi_threshold": 0.50,
        "grid": {"width": width, "height": height, "transform": list(transform)},
    }
    before_scenes = query_scenes(*PERIODS["before"], max_scenes)
    after_scenes = query_scenes(*PERIODS["after"], max_scenes)
    before_path = DATA_DIR / "ndvi_2019_median.tif"
    after_path = DATA_DIR / "ndvi_2024_median.tif"
    metadata["before"] = make_composite("2019", before_scenes, before_path)
    metadata["after"] = make_composite("2024", after_scenes, after_path)
    _, clearing_path, clearing_pixels, valid_pixels = calculate_change(before_path, after_path)
    hansen_path, hansen_pixels = make_hansen_loss()
    _, overlap_pixels = calculate_overlap(clearing_path, hansen_path)
    ndvi_ha = clearing_pixels * (RESOLUTION * RESOLUTION / 10000)
    hansen_ha = hansen_pixels * (HANSEN_RESOLUTION * HANSEN_RESOLUTION / 10000)
    overlap_ha = overlap_pixels * (RESOLUTION * RESOLUTION / 10000)
    results = {
        "ndvi_flagged_ha": ndvi_ha,
        "hansen_loss_ha": hansen_ha,
        "agreement_ha": overlap_ha,
        "ndvi_overlap_pct": 100 * overlap_ha / ndvi_ha if ndvi_ha else 0,
        "hansen_overlap_pct": 100 * overlap_ha / hansen_ha if hansen_ha else 0,
        "iou_pct": 100 * overlap_ha / (ndvi_ha + hansen_ha - overlap_ha)
        if (ndvi_ha + hansen_ha - overlap_ha)
        else 0,
        "valid_analysis_area_ha": valid_pixels * (RESOLUTION * RESOLUTION / 10000),
        "before_scene_count": len(before_scenes),
        "after_scene_count": len(after_scenes),
    }
    metadata["results"] = results
    (ROOT / "results.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    map_png, map_pdf = render_map(results)
    write_report_md(results)
    report_pdf = build_report_pdf(results, map_png)
    qgis_path = write_qgis_project()
    print(json.dumps({"results": results, "map_pdf": str(map_pdf), "report_pdf": str(report_pdf), "qgis": str(qgis_path)}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-scenes", type=int, default=MAX_SCENES)
    args = parser.parse_args()
    os.environ.setdefault("AWS_NO_SIGN_REQUEST", "YES")
    os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
    os.environ.setdefault("GDAL_HTTP_MULTIRANGE", "YES")
    os.environ.setdefault("VSI_CACHE", "TRUE")
    os.environ.setdefault("VSI_CACHE_SIZE", "50000000")
    with rasterio.Env(
        AWS_NO_SIGN_REQUEST="YES",
        GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
        GDAL_HTTP_MULTIRANGE="YES",
    ):
        run(args.max_scenes)
