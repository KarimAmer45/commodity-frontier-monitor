"""Promote verified Earth Engine exports and rebuild final deliverables."""

from __future__ import annotations

import importlib.util
import json
import shutil
from datetime import date
from pathlib import Path

import numpy as np
import rasterio
from rasterio.enums import Resampling
from rasterio.warp import reproject

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "processed"
SOURCE = DATA / "earth_engine_exports"
LOCAL_BACKUP = DATA / "local_mirror"

CONSOLE_RESULTS = {
    "ndvi_flagged_ha": 7160.291200839292,
    "hansen_loss_ha": 12446.156519155735,
    "agreement_ha": 3392.1479100055817,
    "ndvi_overlap_pct": 47.37444071559512,
    "hansen_overlap_pct": 27.254581804308557,
    "iou_pct": 20.920717821658382,
}


def load_builder():
    path = Path(__file__).with_name("run_local_analysis.py")
    spec = importlib.util.spec_from_file_location("local_analysis_builder", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def preserve_local_mirror():
    LOCAL_BACKUP.mkdir(parents=True, exist_ok=True)
    for name in (
        "ndvi_2019_median.tif",
        "ndvi_2024_median.tif",
        "ndvi_diff.tif",
        "clearing_mask.tif",
        "forest_loss_2019_2024.tif",
        "agreement_mask.tif",
    ):
        source = DATA / name
        destination = LOCAL_BACKUP / name
        if source.exists() and not destination.exists():
            shutil.copy2(source, destination)


def promote_exports():
    expected = {
        "ndvi_diff.tif": SOURCE / "ndvi_diff.tif",
        "clearing_mask.tif": SOURCE / "clearing_mask.tif",
        "forest_loss_2020_2024.tif": SOURCE / "forest_loss_2020_2024.tif",
    }
    missing = [str(path) for path in expected.values() if not path.exists()]
    if missing:
        raise FileNotFoundError(f"Missing Earth Engine exports: {missing}")
    for name, source in expected.items():
        shutil.copy2(source, DATA / name)


def make_cartographic_agreement():
    clearing_path = DATA / "clearing_mask.tif"
    hansen_path = DATA / "forest_loss_2020_2024.tif"
    output_path = DATA / "agreement_mask.tif"
    with rasterio.open(clearing_path) as clearing, rasterio.open(hansen_path) as hansen:
        clearing_array = clearing.read(1) == 1
        hansen_on_10m = np.zeros((clearing.height, clearing.width), dtype="uint8")
        reproject(
            source=hansen.read(1),
            destination=hansen_on_10m,
            src_transform=hansen.transform,
            src_crs=hansen.crs,
            dst_transform=clearing.transform,
            dst_crs=clearing.crs,
            resampling=Resampling.nearest,
            src_nodata=None,
            dst_nodata=0,
        )
        agreement = (clearing_array & (hansen_on_10m == 1)).astype("uint8")
        profile = clearing.profile.copy()
        profile.update(
            dtype="uint8",
            nodata=0,
            compress="DEFLATE",
            predictor=2,
            tiled=True,
            blockxsize=256,
            blockysize=256,
        )
        with rasterio.open(output_path, "w", **profile) as destination:
            destination.write(agreement, 1)
    return int(agreement.sum())


def raster_check(path: Path):
    with rasterio.open(path) as dataset:
        array = dataset.read(1)
        return {
            "file": path.name,
            "crs": str(dataset.crs),
            "resolution_m": list(dataset.res),
            "width": dataset.width,
            "height": dataset.height,
            "bounds": list(dataset.bounds),
            "dtype": dataset.dtypes[0],
            "min": float(np.nanmin(array)),
            "max": float(np.nanmax(array)),
            "one_pixels": int(np.count_nonzero(array == 1)),
        }


def write_results(cartographic_agreement_pixels: int):
    metadata = {
        "created": date.today().isoformat(),
        "analysis_engine": "Google Earth Engine",
        "aoi_name": "Tapung palm-oil frontier, Riau, Indonesia",
        "bbox_wgs84": [101.30, 0.30, 101.65, 0.65],
        "target_crs": "EPSG:32647",
        "periods": {
            "ndvi_before": ["2019-01-01", "2020-01-01"],
            "ndvi_after": ["2024-01-01", "2025-01-01"],
            "hansen_loss": [2020, 2024],
        },
        "datasets": {
            "sentinel2": "COPERNICUS/S2_SR_HARMONIZED",
            "cloud_mask": "GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED",
            "hansen": "UMD/hansen/global_forest_change_2025_v1_13",
        },
        "parameters": {
            "scene_cloud_limit_pct": 40,
            "cloud_score_band": "cs",
            "cloud_score_clear_threshold": 0.60,
            "ndvi_change_threshold": -0.20,
            "baseline_ndvi_threshold": 0.50,
            "console_area_scale_m": 10,
        },
        "results": CONSOLE_RESULTS,
        "exports": [
            raster_check(DATA / "ndvi_diff.tif"),
            raster_check(DATA / "clearing_mask.tif"),
            raster_check(DATA / "forest_loss_2020_2024.tif"),
            raster_check(DATA / "agreement_mask.tif"),
        ],
        "cartographic_overlap": {
            "agreement_pixels_10m": cartographic_agreement_pixels,
            "agreement_ha": cartographic_agreement_pixels * 0.01,
            "note": "Derived after nearest-neighbour alignment of exported 30 m Hansen pixels to the 10 m clearing grid; authoritative metrics are the Earth Engine Console values above.",
        },
    }
    (ROOT / "results.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")


def main():
    preserve_local_mirror()
    promote_exports()
    cartographic_agreement_pixels = make_cartographic_agreement()
    write_results(cartographic_agreement_pixels)
    builder = load_builder()
    results = CONSOLE_RESULTS.copy()
    map_png, map_pdf = builder.render_map(results)
    builder.write_report_md(results)
    report_pdf = builder.build_report_pdf(results, map_png)
    qgis_project = builder.write_qgis_project()
    print(
        json.dumps(
            {
                "results": results,
                "map_pdf": str(map_pdf),
                "report_pdf": str(report_pdf),
                "qgis_project": str(qgis_project),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
