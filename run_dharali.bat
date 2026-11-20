@echo off
rem Step 12: second event - Dharali (Uttarkashi), 5 August 2025 debris flow. Tests the tool on a NEW place with no tuning:
rem fresh DEM -> (a) fully automatic run, (b) run with one manual apex -> check both against published facts -> web viewer.
rem About 2 minutes. Log: dharali_status.txt
setlocal
cd /d "%~dp0"
echo [%date% %time%] start > dharali_status.txt
if exist data\dharali_2025\dem.tif goto have_dem
python tools\get_dem.py 30.93 78.70 31.12 78.92 data\dharali_2025\dem.tif > dharali_1_dem.txt 2>&1
if errorlevel 1 ( set WHY=DEM download failed & set LOG=dharali_1_dem.txt & goto fail )
echo DEM downloaded >> dharali_status.txt
:have_dem
set CHK=--lon 78.7792 --lat 31.0412 --catchment-km2 17 --melton 0.8 --deposit-ha 20
echo. >> dharali_status.txt
echo ===== (a) FULLY AUTOMATIC ===== >> dharali_status.txt
python -m fanrisk run examples\dharali_2025\site_auto.toml > dharali_2_auto.txt 2>&1
if errorlevel 1 ( set WHY=the automatic run failed & set LOG=dharali_2_auto.txt & goto fail )
python tools\check_event.py examples\dharali_2025\site_auto.toml %CHK% >> dharali_status.txt 2>&1
echo. >> dharali_status.txt
echo ===== (b) ONE MANUAL APEX ===== >> dharali_status.txt
python -m fanrisk run examples\dharali_2025\site.toml > dharali_3_run.txt 2>&1
if errorlevel 1 ( set WHY=the run failed & set LOG=dharali_3_run.txt & goto fail )
python tools\check_event.py examples\dharali_2025\site.toml %CHK% >> dharali_status.txt 2>&1
python -m fanrisk viewer examples\dharali_2025\site.toml > dharali_5_viewer.txt 2>&1
echo [%date% %time%] ALL DONE - viewer: runs\dharali_2025\viewer.html >> dharali_status.txt
type dharali_status.txt
pause
exit /b 0
:fail
echo STOPPED: %WHY% >> dharali_status.txt
powershell -NoProfile -Command "Get-Content '%LOG%' -Tail 15" >> dharali_status.txt
type dharali_status.txt
pause
exit /b 1
