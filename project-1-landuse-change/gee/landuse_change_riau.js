// Land-use change over the Tapung palm-oil frontier, Riau, 2019-2024.
// Run in the Google Earth Engine Code Editor: https://code.earthengine.google.com/

var aoi = ee.Geometry.Rectangle([101.30, 0.30, 101.65, 0.65]);
Map.centerObject(aoi, 11);

// Cloud masking via Google Cloud Score+ (more reliable than QA60, which was
// empty/unreliable for much of 2022-2024). 'cs' grades pixels 0 (occluded)..1 (clear).
var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
var CS_BAND = 'cs';
var CLEAR_THRESHOLD = 0.60;  // 0.5-0.65 typical; raise to be stricter

function maskS2clouds(image) {
  return image.updateMask(image.select(CS_BAND).gte(CLEAR_THRESHOLD))
      .divide(10000)
      .copyProperties(image, ['system:time_start']);
}

function addNDVI(image) {
  return image.addBands(
      image.normalizedDifference(['B8', 'B4']).rename('NDVI'));
}

// Link each S2 scene to its Cloud Score+ scene so maskS2clouds can read 'cs'.
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .linkCollection(csPlus, [CS_BAND]);

function ndviComposite(startDate, endDate) {
  return s2.filterBounds(aoi)
      .filterDate(startDate, endDate)
      .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 40))
      .map(maskS2clouds)
      .map(addNDVI)
      .select('NDVI')
      .median()
      .clip(aoi);
}

var ndviBefore = ndviComposite('2019-01-01', '2020-01-01');
var ndviAfter = ndviComposite('2024-01-01', '2025-01-01');
var ndviDiff = ndviAfter.subtract(ndviBefore).rename('NDVI_diff');

var clearing = ndviDiff.lt(-0.20)
    .and(ndviBefore.gt(0.50))
    .rename('clearing');

// Current catalog release: v1.13, covering 2000-2025.
var gfc = ee.Image('UMD/hansen/global_forest_change_2025_v1_13');
var lossYear = gfc.select('lossyear');
// Loss 2020-2024: strictly after the 2019 baseline (2019 loss is already in the baseline NDVI).
var lossRecent = lossYear.gte(20).and(lossYear.lte(24))
    .rename('forest_loss_2020_2024');
var agree = clearing.and(lossRecent).rename('agreement');

var ndviViz = {
  min: 0,
  max: 0.8,
  palette: ['white', 'yellow', 'green', 'darkgreen']
};
var diffViz = {
  min: -0.4,
  max: 0.4,
  palette: ['b2182b', 'f7f7f7', '1a9850']
};

Map.addLayer(ndviBefore, ndviViz, 'NDVI 2019');
Map.addLayer(ndviAfter, ndviViz, 'NDVI 2024');
Map.addLayer(ndviDiff, diffViz, 'NDVI change');
Map.addLayer(clearing.selfMask(), {palette: ['ff00ff']}, 'NDVI clearing');
Map.addLayer(lossRecent.selfMask().clip(aoi), {palette: ['ff8c00']}, 'Hansen loss 2020-2024');
Map.addLayer(agree.selfMask(), {palette: ['111111']}, 'Agreement');

var pixelHa = ee.Image.pixelArea().divide(10000);

function areaHa(maskImage, label, scale) {
  var area = maskImage.selfMask().multiply(pixelHa).reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: aoi,
    scale: scale,
    maxPixels: 1e13,
    tileScale: 4
  });
  var ha = ee.Number(area.values().get(0));
  print(label, ha);  // print the number, not the raw dictionary
  return ha;
}

// All three at 10 m so precision / recall / IoU sit on one common grid.
var clearingHa = areaHa(clearing, 'NDVI-flagged clearing (ha)', 10);
var hansenHa = areaHa(lossRecent, 'Hansen loss 2020-2024 (ha)', 10);
var agreementHa = areaHa(agree, 'Agreement (ha)', 10);

print('Share of NDVI flags overlapping Hansen (%)',
    agreementHa.divide(clearingHa).multiply(100));
print('Share of Hansen loss overlapping NDVI flags (%)',
    agreementHa.divide(hansenHa).multiply(100));
print('Intersection over union (%)', agreementHa
    .divide(clearingHa.add(hansenHa).subtract(agreementHa))
    .multiply(100));

Export.image.toDrive({
  image: ndviDiff.toFloat(),
  description: 'ndvi_diff',
  folder: 'earth_engine_exports',
  region: aoi,
  scale: 10,
  crs: 'EPSG:32647',
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: lossRecent.toByte(),
  description: 'forest_loss_2020_2024',
  folder: 'earth_engine_exports',
  region: aoi,
  scale: 30,
  crs: 'EPSG:32647',
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});

Export.image.toDrive({
  image: clearing.toByte(),
  description: 'clearing_mask',
  folder: 'earth_engine_exports',
  region: aoi,
  scale: 10,
  crs: 'EPSG:32647',
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
