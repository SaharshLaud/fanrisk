# Data

| File | What it is | Source |
| --- | --- | --- |
| `dem_utm.tif` | Leh DEM, EPSG:32643, 26.7 m | Copernicus GLO-30, reprojected |
| `map.osm` | OpenStreetMap extract, 34.086-34.134 N, 77.539-77.630 E | (c) OpenStreetMap contributors, ODbL |
| `validation/footprint_2010.geojson` | 2010 Choglamsar/Sabu flooded area, 0.77 km2 | traced from Bhatt et al. (2011) Fig. 7, georeferenced on imagery (rms 17 m) |
| `validation/damage_points_2010.geojson` | 648 damaged houses; `holdout=1` on one third | same figure |
| `validation/survey_window_2010.geojson` | the area Bhatt et al. mapped (only this is scored) | same figure |
| `validation/dharali_2025/dharali_points.geojson` | centre of the Dharali deposit | NRSC/NDEM 2025 map frame |

Downloaded when needed (not in the repository):

- `dharali_2025/dem.tif` - `python tools/get_dem.py 30.93 78.70 31.12 78.92 data/dharali_2025/dem.tif`
- `ms_buildings_*.csv.gz` - `python tools/get_ms_buildings.py 34.086 77.539 34.134 77.630 data`
