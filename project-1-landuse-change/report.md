# Land-use change over the Tapung palm-oil frontier, Riau, 2019-2024

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

The screening flagged **7,160.3 ha** of likely vegetation clearing. Hansen recorded **12,446.2 ha** of forest loss for 2020-2024. Their spatial intersection was **3,392.1 ha**. This equals **47.4%** of the NDVI-flagged area and **27.3%** of Hansen loss; IoU was **20.9%**.

![Land-use change map](map/landuse_change_map.png)

These totals are screening indicators. Non-overlap is expected because NDVI detects vegetation changes outside Hansen's >5 m tree-cover definition, while Hansen can capture stand-replacement events whose annual median NDVI signal is muted by regrowth, timing, mixed pixels, or cloud availability.

## 5. Limitations

1. **Attribution:** neither NDVI nor Hansen identifies the responsible commodity, actor, or legal status. Parcel boundaries and field evidence are required.
2. **Phenology and management:** crop rotation, harvesting, fire, flooding, drought, and plantation cycles can resemble clearing.
3. **Threshold sensitivity:** -0.20, 0.50, and Cloud Score+ 0.60 are transparent screening choices, not locally calibrated decision rules.
4. **Cloud and compositing:** cloud scores can miss haze or shadow, and annual medians can suppress short-lived events.
5. **Product mismatch:** Sentinel-2 is evaluated at 10 m while Hansen is Landsat-derived at 30.92 m; alignment creates mixed-pixel and edge effects.
6. **Definition mismatch:** Hansen forest loss is stand-replacement disturbance of vegetation taller than 5 m; NDVI decline is a broader signal.
7. **Cartographic resampling:** the map is rendered from the exported 10 m and 30 m GeoTIFFs, so overlap read off the map can differ slightly from the common-scale (10 m) server measurement reported above.

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
