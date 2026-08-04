# Commodity Frontier Monitor

**Satellite land-use-change screening over a palm-oil frontier in Riau, Indonesia** — Sentinel-2 NDVI change detection, cross-validated against Hansen Global Forest Change, built in Google Earth Engine and mapped in QGIS.

![Google Earth Engine](https://img.shields.io/badge/Google%20Earth%20Engine-JavaScript-4285F4?style=flat-square&logo=googleearthengine&logoColor=white)
![Sentinel-2](https://img.shields.io/badge/Sentinel--2-SR%20Harmonized-1a9850?style=flat-square)
![Hansen GFC](https://img.shields.io/badge/Hansen%20GFC-v1.13-ff8c00?style=flat-square)
![QGIS](https://img.shields.io/badge/QGIS-Print%20Layout-589632?style=flat-square&logo=qgis&logoColor=white)
![Python](https://img.shields.io/badge/Python-3.12-3776AB?style=flat-square&logo=python&logoColor=white)

![Land-use change over the Tapung palm-oil frontier, Riau, 2019–2024](project-1-landuse-change/map/landuse_change_map.png)

<sub>Change over the Tapung frontier, 2019→2024. **Magenta** = NDVI-flagged clearing · **orange** = Hansen forest loss (2020–2024) · **red–green** = NDVI change, over Sentinel-2 / Esri imagery. AOI ≈ 1,510 km², UTM zone 47N.</sub>

---

## What this is

A small, fully reproducible pipeline that answers a supply-chain-relevant question: **where did dense vegetation decline over a commodity frontier, and how well does a simple satellite signal agree with an authoritative independent product?** It builds cloud-masked Sentinel-2 NDVI composites for 2019 and 2024, flags likely clearing from the NDVI drop, and measures how much of that overlaps Hansen Global Forest Change — reporting the agreement honestly, including where the method misses.

## Key results

| Metric | Value |
| --- | --- |
| Area screened | ~1,510 km² |
| NDVI-flagged clearing | **7,160 ha** |
| Hansen forest loss (2020–2024) | **12,446 ha** |
| Agreement (both methods) | **3,392 ha** |
| Precision — NDVI flags corroborated by Hansen | **47.4%** |
| Recall — Hansen loss caught by NDVI | **27.3%** |
| Intersection over union | **20.9%** |

**How to read this:** the method is deliberately conservative — high-confidence flags, ~47% independently corroborated, but it only catches ~27% of Hansen's loss. Most of the miss is *rotational* plantation loss (harvest/replant) whose annual-median NDVI signal recovers before the composite window closes. That trade-off, and its supply-chain implications, is documented in the [full report](project-1-landuse-change/report.md#5-limitations).

## Method at a glance

1. **Cloud masking** — Google Cloud Score+ (`cs ≥ 0.60`), more reliable than the legacy QA60 band for 2022–2024 imagery.
2. **Composites** — per-pixel annual **median** NDVI for 2019 and 2024 from `COPERNICUS/S2_SR_HARMONIZED`.
3. **Change flag** — `NDVI change < −0.20` **and** `2019 NDVI > 0.50` (a real drop, from vegetated ground).
4. **Cross-validation** — intersect with Hansen GFC v1.13 loss for **2020–2024** (strictly after the 2019 baseline), and quantify precision, recall, and IoU on a common 10 m grid.

## Repository

| Path | What |
| --- | --- |
| [`project-1-landuse-change/report.md`](project-1-landuse-change/report.md) | Methods, results, and a candid limitations section |
| [`project-1-landuse-change/gee/landuse_change_riau.js`](project-1-landuse-change/gee/landuse_change_riau.js) | Earth Engine pipeline (authoritative source of the figures) |
| [`project-1-landuse-change/map/landuse_change_map.png`](project-1-landuse-change/map/landuse_change_map.png) | QGIS map (PNG) · [PDF](project-1-landuse-change/output/pdf/Riau_LandUse_2019_2024.pdf) |
| [`project-1-landuse-change/qgis/tapung_landuse_change.qgs`](project-1-landuse-change/qgis/tapung_landuse_change.qgs) | Styled QGIS project |
| [`project-1-landuse-change/results.json`](project-1-landuse-change/results.json) | Measurements, parameters, and scene IDs |

Full details and reproduction steps: **[`project-1-landuse-change/`](project-1-landuse-change/README.md)**.

## Reproduce

Paste [`gee/landuse_change_riau.js`](project-1-landuse-change/gee/landuse_change_riau.js) into the [Earth Engine Code Editor](https://code.earthengine.google.com), run it, and submit the Drive exports. The Console prints the areas and agreement metrics; the GeoTIFFs feed the QGIS layout.

## Data & credits

Sentinel-2 SR Harmonized (ESA/Copernicus) · Google Cloud Score+ · Hansen Global Forest Change v1.13 (Hansen et al., 2013) · basemap imagery Esri/Google. Analysis CRS: WGS 84 / UTM zone 47N.
