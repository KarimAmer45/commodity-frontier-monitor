# Land-use change and deforestation monitoring

Reproducible screening of a palm-oil frontier in Tapung, Riau, Indonesia,
using Sentinel-2 NDVI change and Hansen Global Forest Change.

![Land-use change over the Tapung palm-oil frontier, Riau, 2019–2024](map/landuse_change_map.png)

<sub>Magenta = NDVI-flagged clearing · orange = Hansen forest loss (2020–2024) · red–green = NDVI change, over Sentinel-2 / Esri imagery.</sub>

## AOI and decision rule

- AOI: `101.30, 0.30, 101.65, 0.65` (WGS 84), approximately 39 x 39 km.
- Before / after: 2019 and 2024 annual median NDVI.
- Likely clearing: `NDVI change < -0.20` and `2019 NDVI > 0.50`.
- Independent comparison: Hansen GFC v1.13 loss years 20-24 (2020-2024), strictly after the 2019 baseline year.

## Deliverables

- `gee/landuse_change_riau.js` - Google Earth Engine implementation.
- `data/processed/ndvi_diff.tif` - 10 m NDVI-change GeoTIFF.
- `data/processed/forest_loss_2020_2024.tif` - 30 m Hansen-loss GeoTIFF.
- `data/processed/clearing_mask.tif` and `agreement_mask.tif` - diagnostic masks.
- `qgis/tapung_landuse_change.qgs` - QGIS-ready styled project.
- `map/landuse_change_map.png` and `output/pdf/Riau_LandUse_2019_2024.pdf` - map exports.
- `report.md` and `output/pdf/landuse_change_report.pdf` - methods report.
- `results.json` - measurements, parameters, and selected scene IDs.
- `portfolio_text.md` - filled-in CV bullets, cover-letter hook, and interview caveat.

## Authoritative workflow

The final figures and canonical GeoTIFFs come from Google Earth Engine. The
workflow uses Sentinel-2 SR Harmonized, Google Cloud Score+ (`cs >= 0.60`),
annual median NDVI composites for 2019 and 2024, and Hansen loss for 2020-2024.
Paste `gee/landuse_change_riau.js` into the Earth Engine Code Editor and submit
the Drive exports. Console measurements and export metadata are recorded in
`results.json`.

After downloading refreshed exports into `data/processed/earth_engine_exports`, run:

```powershell
python scripts/finalize_gee_exports.py
```

This preserves earlier local outputs under `data/processed/local_mirror`,
promotes the Earth Engine rasters, derives the cartographic agreement mask, and
rebuilds the map, report, QGIS project, results metadata, and portfolio text.

## Anonymous local comparison

The script uses anonymous public Cloud-Optimized GeoTIFFs. Install its packages
and run from this folder:

```powershell
python -m pip install -r requirements.txt
python scripts/run_local_analysis.py
```

The local mirror uses Sentinel-2 SCL quality masking and the least-cloudy unique
annual acquisitions to keep anonymous execution practical. It is an independent
reproducibility aid and is not expected to match the authoritative Cloud Score+
Earth Engine result pixel for pixel.

## Use in Earth Engine and QGIS

1. Paste `gee/landuse_change_riau.js` into the Earth Engine Code Editor, run it,
   review console statistics, and start the three Drive exports.
2. Open `qgis/tapung_landuse_change.qgs` in QGIS 3.x. Relative paths keep the
   raster sources connected when the repository is moved intact.
3. The supplied PNG and PDF map already include a diverging NDVI ramp, overlays,
   legend, north arrow, scale bar, projected coordinates, data credits, and area
   summary.

## Interpretation

This is a screening method. A flag is not proof of commodity-driven or illegal
deforestation. Operational decisions require parcel data, dated high-resolution
validation, threshold calibration, and uncertainty reporting.
