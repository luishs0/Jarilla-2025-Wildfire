// ================== SETUP ==================
var PERIMETRO_ASSET = 'projects/iberia2025/assets/perimetro_incendio';
var aoi = ee.FeatureCollection(PERIMETRO_ASSET).geometry();
var UMBRAL_PENDIENTE = 15;
var UMBRAL_ALINEACION = 0.5;
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
var ndviMediana = ee.Number(NDVI_pre.reduceRegion({
  reducer: ee.Reducer.median(), geometry: aoi, scale: ESCALA, maxPixels:1e13, bestEffort:true
}).get('NDVI'));
var combustible = NDVI_pre.gte(ee.Image.constant(ndviMediana));

// ================== capas base ==================
var ha = ee.Image.pixelArea().divide(10000);
var validez = severidad.gte(1);
var altaSev = severidad.gte(5);
function alinDe(ang){ return aspect.subtract(ang).multiply(Math.PI/180).cos(); }

var m_slope     = slope.gte(UMBRAL_PENDIENTE).and(validez);
var m_slopefuel = m_slope.and(combustible);
var m_ssw = m_slopefuel.and(alinDe(200).gte(UMBRAL_ALINEACION));
var m_ne  = m_slopefuel.and(alinDe(45).gte(UMBRAL_ALINEACION));

// ---- bandas de contribucion (descomposicion, con dNBR) ----
function bandsFull(name, mask){
  var m  = mask.unmask(0);
  var mh = mask.and(altaSev).unmask(0);
  return ee.Image([
    ha.multiply(m).rename('area_'+name),
    ha.multiply(mh).rename('hs_'+name),
    dNBR.unmask(0).multiply(m).rename('dsum_'+name),
    m.rename('cnt_'+name)
  ]);
}
// ---- bandas ligeras (barrido, solo area y sev alta) ----
function bandsLite(name, mask){
  var m  = mask.unmask(0);
  var mh = mask.and(altaSev).unmask(0);
  return ee.Image([ ha.multiply(m).rename('area_'+name), ha.multiply(mh).rename('hs_'+name) ]);
}

var imgs = [
  bandsFull('base',      validez),
  bandsFull('pend',      m_slope),
  bandsFull('pendcomb',  m_slopefuel),
  bandsFull('ssw',       m_ssw),
  bandsFull('ne',        m_ne)
];
var angs = [0,30,60,90,120,150,180,210,240,270,300,330];
angs.forEach(function(ang){
  var m = m_slopefuel.and(alinDe(ang).gte(UMBRAL_ALINEACION));
  imgs.push(bandsLite('a'+ang, m));
});

// ================== UNA sola agregacion ==================
var stats = ee.Image.cat(imgs).reduceRegion({
  reducer: ee.Reducer.sum(), geometry: aoi, scale: ESCALA, maxPixels: 1e13, bestEffort: true
});
var totalAlta = ee.Number(stats.get('hs_base'));

function filaFull(caso, name){
  var a=ee.Number(stats.get('area_'+name)), hs=ee.Number(stats.get('hs_'+name));
  var ds=ee.Number(stats.get('dsum_'+name)), cnt=ee.Number(stats.get('cnt_'+name));
  return ee.Feature(null,{caso:caso, area_ha:a, dNBR_medio:ds.divide(cnt),
                          pct_sevAlta:hs.divide(a).multiply(100), recall_pct:hs.divide(totalAlta).multiply(100)});
}

// ---- tabla 1: descomposicion ----
var descomp = ee.FeatureCollection([
  filaFull('1_base_todo','base'),
  filaFull('2_pendiente','pend'),
  filaFull('3_pend+combustible','pendcomb'),
  filaFull('4_corredor_SSW(200)','ssw'),
  filaFull('4b_corredor_NE(45)','ne')
]);
print('DESCOMP - casos:',         descomp.aggregate_array('caso'));
print('DESCOMP - pct_sevAlta %:', descomp.aggregate_array('pct_sevAlta'));
print('DESCOMP - recall %:',      descomp.aggregate_array('recall_pct'));
print('DESCOMP - dNBR medio:',    descomp.aggregate_array('dNBR_medio'));
print('DESCOMP - area ha:',       descomp.aggregate_array('area_ha'));

// ---- tabla 2: barrido ----
var filasBarrido = angs.map(function(ang){
  var a=ee.Number(stats.get('area_a'+ang)), hs=ee.Number(stats.get('hs_a'+ang));
  return ee.Feature(null,{procedencia_grados:ang, area_ha:a,
                          pct_sevAlta:hs.divide(a).multiply(100), recall_pct:hs.divide(totalAlta).multiply(100)});
});
var barrido = ee.FeatureCollection(filasBarrido);
print('BARRIDO - procedencia:',   barrido.aggregate_array('procedencia_grados'));
print('BARRIDO - pct_sevAlta %:', barrido.aggregate_array('pct_sevAlta'));
print('BARRIDO - recall %:',      barrido.aggregate_array('recall_pct'));

// ================== MAPA + EXPORT (corredor SSW) ==================
var corredor = m_ssw.rename('corredor');
Map.addLayer(severidad, {min:1,max:6,palette:['1b9e77','c2c2c2','ffffb2','fecc5c','fd8d3c','e31a1c']}, 'Severidad');
Map.addLayer(corredor.selfMask(), {palette:['black']}, 'Corredores SSW');
Export.image.toDrive({
  image: corredor.toByte(), description:'jarilla_corredores_SSW', folder:CARPETA_DRIVE,
  fileNamePrefix:'jarilla_corredores_SSW', region:aoi, scale:ESCALA, crs:EPSG, maxPixels:1e13
});
