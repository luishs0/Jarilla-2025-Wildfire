var USAR_ASSET = true;
var PERIMETRO_ASSET = 'users/LHernandezS/incendio_jarilla';
var BBOX_RESPALDO = ee.Geometry.Rectangle([-6.30, 40.05, -5.75, 40.40]);
var PRE_INI  = '2025-07-15';
var PRE_FIN  = '2025-08-11';
var POST_INI = '2025-08-25';
var POST_FIN = '2025-09-30';
var MAX_NUBES = 60;
var USAR_CSPLUS = true;
var CS_UMBRAL = 0.60;
var EPSG = 'EPSG:25830';
var ESCALA = 20;
var CARPETA_DRIVE = 'TFG_Jarilla_Sentinel2';
var aoi = USAR_ASSET ? ee.FeatureCollection(PERIMETRO_ASSET).geometry()
                     : BBOX_RESPALDO;
Map.centerObject(aoi, 11);
Map.addLayer(aoi, {color: 'black'}, 'Perímetro', false);
function cargarS2(ini, fin) {
  var col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
              .filterBounds(aoi)
              .filterDate(ini, fin)
              .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', MAX_NUBES));
  if (USAR_CSPLUS) {
    var csp = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
    col = col.linkCollection(csp, ['cs']);
    col = col.map(function(img) {
      var clara = img.select('cs').gte(CS_UMBRAL);
      return img.updateMask(clara)
                .divide(10000)
                .copyProperties(img, img.propertyNames());
    });
  } else {
    col = col.map(function(img) {
      var scl = img.select('SCL');
      var buena = scl.eq(4).or(scl.eq(5)).or(scl.eq(6)).or(scl.eq(7)).or(scl.eq(11));
      return img.updateMask(buena)
                .divide(10000)
                .copyProperties(img, img.propertyNames());
    });
  }
  return col;
}
var s2pre  = cargarS2(PRE_INI,  PRE_FIN);
var s2post = cargarS2(POST_INI, POST_FIN);
print('Nº escenas PRE  disponibles:', s2pre.size());
print('Nº escenas POST disponibles:', s2post.size());
var pre  = s2pre.median().clip(aoi);
var post = s2post.median().clip(aoi);
function ndvi(img){ return img.normalizedDifference(['B8','B4']).rename('NDVI'); }
function ndii(img){ return img.normalizedDifference(['B8','B11']).rename('NDII'); }
function nbr (img){ return img.normalizedDifference(['B8','B12']).rename('NBR'); }
var NDVI_pre = ndvi(pre);
var NDII_pre = ndii(pre);
var NBR_pre  = nbr(pre);
var NBR_post = nbr(post);
var dNBR = NBR_pre.subtract(NBR_post).multiply(1000).rename('dNBR');
var severidad = dNBR
  .where(dNBR.lt(-100), 1)
  .where(dNBR.gte(-100).and(dNBR.lt(100)), 2)
  .where(dNBR.gte(100).and(dNBR.lt(270)), 3)
  .where(dNBR.gte(270).and(dNBR.lt(440)), 4)
  .where(dNBR.gte(440).and(dNBR.lt(660)), 5)
  .where(dNBR.gte(660), 6)
  .rename('severidad');
Map.addLayer(NDVI_pre, {min:0, max:0.9, palette:['white','khaki','green']}, 'NDVI pre');
Map.addLayer(NDII_pre, {min:-0.2, max:0.5, palette:['brown','white','blue']}, 'NDII pre');
Map.addLayer(dNBR, {min:-100, max:1000, palette:['green','yellow','orange','red','purple']}, 'dNBR');
var palSev = ['1b9e77','c2c2c2','ffffb2','fecc5c','fd8d3c','e31a1c'];
Map.addLayer(severidad, {min:1, max:6, palette:palSev}, 'Severidad (Key&Benson)');
var percentiles = ee.Reducer.percentile([5,10,25,50,75,90,95])
                    .combine(ee.Reducer.mean(), '', true)
                    .combine(ee.Reducer.stdDev(), '', true);
var statsNDVI = NDVI_pre.reduceRegion({reducer:percentiles, geometry:aoi, scale:ESCALA, maxPixels:1e13, bestEffort:true});
var statsNDII = NDII_pre.reduceRegion({reducer:percentiles, geometry:aoi, scale:ESCALA, maxPixels:1e13, bestEffort:true});
print('--- NDVI pre (estadísticos) ---', statsNDVI);
print('--- NDII pre (estadísticos) ---', statsNDII);
var areaImg = ee.Image.pixelArea().divide(10000).addBands(severidad);
var areaPorClase = areaImg.reduceRegion({
  reducer: ee.Reducer.sum().group({groupField:1, groupName:'clase_severidad'}),
  geometry: aoi, scale: ESCALA, maxPixels: 1e13, bestEffort:true
});
print('--- Superficie (ha) por clase de severidad ---', areaPorClase);
var statsdNBR = dNBR.reduceRegion({
  reducer: ee.Reducer.percentile([5,25,50,75,95]).combine(ee.Reducer.mean(),'',true),
  geometry: aoi, scale: ESCALA, maxPixels: 1e13, bestEffort:true});
print('--- dNBR (estadísticos) ---', statsdNBR);
function exportar(img, nombre){
  Export.image.toDrive({
    image: img.toFloat(), description: nombre, folder: CARPETA_DRIVE,
    fileNamePrefix: nombre, region: aoi, scale: ESCALA, crs: EPSG, maxPixels: 1e13
  });
}
exportar(NDVI_pre, 'jarilla_NDVI_pre');
exportar(NDII_pre, 'jarilla_NDII_pre');
exportar(NBR_pre,  'jarilla_NBR_pre');
exportar(NBR_post, 'jarilla_NBR_post');
exportar(dNBR,     'jarilla_dNBR');
exportar(severidad,'jarilla_severidad');
print('>>> Revisar la pestaña TASKS y pulsar RUN en cada exportación.');
