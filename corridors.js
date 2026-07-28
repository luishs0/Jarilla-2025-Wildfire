var PERIMETRO_ASSET = 'projects/iberia2025/assets/perimetro_incendio';
var aoi = ee.FeatureCollection(PERIMETRO_ASSET).geometry();
var VIENTO_PROCEDENCIA = 45;
var UMBRAL_PENDIENTE   = 15;
var UMBRAL_ALINEACION  = 0.5;
var PRE_INI='2025-07-15', PRE_FIN='2025-08-11';
var POST_INI='2025-08-25', POST_FIN='2025-09-30';
var MAX_NUBES=60, CS_UMBRAL=0.60;
var EPSG='EPSG:25830', ESCALA=20, CARPETA_DRIVE='TFG_Jarilla_Sentinel2';
Map.centerObject(aoi, 11);
function cargarS2(ini, fin){
  var col = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
              .filterBounds(aoi).filterDate(ini, fin)
              .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', MAX_NUBES));
  var csp = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');
  col = col.linkCollection(csp, ['cs']).map(function(img){
    return img.updateMask(img.select('cs').gte(CS_UMBRAL))
              .divide(10000).copyProperties(img, img.propertyNames());
  });
  return col;
}
var pre  = cargarS2(PRE_INI, PRE_FIN).median().clip(aoi);
var post = cargarS2(POST_INI, POST_FIN).median().clip(aoi);
var NDVI_pre = pre.normalizedDifference(['B8','B4']).rename('NDVI');
var NBR_pre  = pre.normalizedDifference(['B8','B12']).rename('NBR');
var NBR_post = post.normalizedDifference(['B8','B12']).rename('NBR');
var dNBR = NBR_pre.subtract(NBR_post).multiply(1000).rename('dNBR');
var severidad = dNBR
  .where(dNBR.lt(-100),1).where(dNBR.gte(-100).and(dNBR.lt(100)),2)
  .where(dNBR.gte(100).and(dNBR.lt(270)),3).where(dNBR.gte(270).and(dNBR.lt(440)),4)
  .where(dNBR.gte(440).and(dNBR.lt(660)),5).where(dNBR.gte(660),6).rename('severidad');
var dem   = ee.Image('projects/iberia2025/assets/mdt05_jarilla');
var slope = ee.Terrain.slope(dem).clip(aoi).rename('pendiente');
var aspect= ee.Terrain.aspect(dem).clip(aoi).rename('orientacion');
var alineacion = aspect.subtract(VIENTO_PROCEDENCIA).multiply(Math.PI/180).cos().rename('alineacion');
print('Pendiente media del terreno (grados):',
  slope.reduceRegion({reducer:ee.Reducer.mean(), geometry:aoi, scale:ESCALA, maxPixels:1e13, bestEffort:true}).get('pendiente'));
var ndviMediana = ee.Number(NDVI_pre.reduceRegion({
  reducer: ee.Reducer.median(), geometry: aoi, scale: ESCALA, maxPixels:1e13, bestEffort:true
}).get('NDVI'));
print('NDVI mediana (umbral de combustible):', ndviMediana);
var combustible = NDVI_pre.gte(ee.Image.constant(ndviMediana));
var corredor = slope.gte(UMBRAL_PENDIENTE)
                 .and(alineacion.gte(UMBRAL_ALINEACION))
                 .and(combustible)
                 .rename('corredor');
Map.addLayer(severidad, {min:1,max:6,palette:['1b9e77','c2c2c2','ffffb2','fecc5c','fd8d3c','e31a1c']}, 'Severidad');
Map.addLayer(alineacion, {min:-1,max:1,palette:['blue','white','red']}, 'Alineacion viento-ladera', false);
Map.addLayer(corredor.selfMask(), {palette:['black']}, 'Corredores');
var ha = ee.Image.pixelArea().divide(10000);
function sumaHa(mask){
  return ha.updateMask(mask).reduceRegion({reducer:ee.Reducer.sum(), geometry:aoi, scale:ESCALA, maxPixels:1e13, bestEffort:true}).get('area');
}
function mediaDNBR(mask){
  return dNBR.updateMask(mask).reduceRegion({reducer:ee.Reducer.mean(), geometry:aoi, scale:ESCALA, maxPixels:1e13, bestEffort:true}).get('dNBR');
}
var altaSev = severidad.gte(5);
var fuera = corredor.not();
print('=============== RESULTADOS PARA 5.4 ===============');
print('Superficie de CORREDORES (ha):', sumaHa(corredor));
print('dNBR medio DENTRO de corredores:', mediaDNBR(corredor));
print('dNBR medio FUERA de corredores:', mediaDNBR(fuera));
print('Superficie severidad ALTA (clases 5-6) DENTRO de corredores (ha):', sumaHa(corredor.and(altaSev)));
print('Superficie severidad ALTA (clases 5-6) FUERA de corredores (ha):', sumaHa(fuera.and(altaSev)));
print('% del corredor que es severidad alta:',
      ee.Number(sumaHa(corredor.and(altaSev))).divide(sumaHa(corredor)).multiply(100));
print('% del area NO-corredor que es severidad alta:',
      ee.Number(sumaHa(fuera.and(altaSev))).divide(sumaHa(fuera)).multiply(100));
Export.image.toDrive({
  image: corredor.toByte(), description:'jarilla_corredores', folder:CARPETA_DRIVE,
  fileNamePrefix:'jarilla_corredores', region:aoi, scale:ESCALA, crs:EPSG, maxPixels:1e13
});
print('>>> Revisar los numeros del bloque "RESULTADOS PARA 5.4" y exportar el mapa (pestaña Tasks).');
