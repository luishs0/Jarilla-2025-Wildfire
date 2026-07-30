# Jarilla Wildfire 2025

Analysis code accompanying the study of the **2025 Jarilla wildfire** (Cáceres, western Spain; 17,317 ha, August 2025), which examines how **meteorology, topography, and fuel condition** jointly controlled the fire's spread and post-fire severity, and validates a reproducible delineation of **fire spread corridors** against the observed severity footprint.

The workflow integrates ERA5 hourly reanalysis, the Fire Weather Index (FWI), Sentinel-2 spectral indices (NDVI, NDII, dNBR), and the 5 m PNOA digital elevation model (MDT05). It relies entirely on open-access data.

## Repository contents

| File | Language | Description |
|------|----------|-------------|
| `topography.py` | Python | Processes the MDT05 (5 m) DEM: slope and aspect (Horn, 1981), wind–slope alignment, per-class statistics, and export of georeferenced layers (GeoTIFF, EPSG:25830) and maps. |
| `index.js` | Google Earth Engine (JavaScript) | Builds cloud-free pre- and post-fire Sentinel-2 composites and computes NDVI, NDII, NBR and dNBR, classifies burn severity (Key & Benson, 2006), and exports rasters and statistics. |
| `corridors.js` | Google Earth Engine (JavaScript) | Delineates fire spread corridors (slope ≥ 15°, favorable wind–slope alignment, NDVI above the perimeter median) and validates them against post-fire severity. |

## Requirements

- **Python 3** with `pandas`, `numpy`, `matplotlib`, `openpyxl`, and `Pillow`.
- A **Google Earth Engine** account to run the `.js` scripts in the [Earth Engine Code Editor](https://code.earthengine.google.com/).

## Data sources

All datasets analyzed are openly available from public repositories:

- **ERA5 hourly data on single levels** (2 m temperature, 2 m dew point, 10 m wind components) — Copernicus Climate Data Store (CDS): https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels (DOI: 10.24381/cds.adbb2d47).
- **Fire Weather Index (FWI)** — Fire danger indices, Copernicus Emergency Management Service (CEMS) / European Forest Fire Information System (EFFIS): https://effis.jrc.ec.europa.eu
- **Sentinel-2 MSI Level-2A** surface reflectance (harmonized collection `COPERNICUS/S2_SR_HARMONIZED`) and **Cloud Score+** (`GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED`), via Google Earth Engine: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR_HARMONIZED
- **5 m Digital Elevation Model (MDT05, PNOA-LiDAR)** — Centro de Descargas, Instituto Geográfico Nacional (IGN), Spain: https://centrodedescargas.cnig.es

## Coordinate reference system

All spatial outputs use **ETRS89 / UTM zone 30N (EPSG:25830)**.

