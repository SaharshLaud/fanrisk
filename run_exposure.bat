@echo off
rem Step 8: who is hit. Downloads Microsoft building footprints once (OSM maps few houses in Ladakh), then overlays buildings,
rem roads and critical facilities (OpenStreetMap: data\map.osm) on the AHP risk map, the HEC-RAS run (2 m channel burn) and the 2010 footprint.
setlocal
cd /d "%~dp0"
echo [%date% %time%] start > exposure_status.txt
if not exist data\map.osm ( echo STOPPED: data\map.osm missing - download it first >> exposure_status.txt & goto end )
dir /b data\ms_buildings_*.csv.gz >nul 2>&1
if errorlevel 1 (
  python tools\get_ms_buildings.py 34.086 77.539 34.134 77.630 data > exposure_1_download.txt 2>&1
  if errorlevel 1 ( echo WARNING: Microsoft buildings not downloaded - using OSM buildings only. See exposure_1_download.txt >> exposure_status.txt )
)
python -m fanrisk exposure examples\leh_2010\site.toml > exposure_2_run.txt 2>&1
if errorlevel 1 ( echo STOPPED: exposure failed - see exposure_2_run.txt >> exposure_status.txt & goto end )
type runs\leh_2010\exposure\exposure_report.txt >> exposure_status.txt
echo [%date% %time%] ALL DONE - map: runs\leh_2010\exposure\exposure_map.png >> exposure_status.txt
:end
type exposure_status.txt
pause
