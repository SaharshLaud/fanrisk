@echo off
rem Kedarnath 2013, step 1 (about 2 min): fresh DEM -> (a) fully automatic run, (b) run with one manual apex -> check against published
rem facts. Step 2 (HEC-RAS with the tool's own storm hydrograph) is run_kedarnath_hecras.bat, set up from what step 1 finds.
rem Log: kedarnath_status.txt
setlocal
cd /d "%~dp0"
echo [%date% %time%] start > kedarnath_status.txt
if exist data\kedarnath_2013\dem.tif goto have_dem
python tools\get_dem.py 30.62 79.00 30.84 79.16 data\kedarnath_2013\dem.tif > kedarnath_1_dem.txt 2>&1
if errorlevel 1 ( set WHY=DEM download failed & set LOG=kedarnath_1_dem.txt & goto fail )
echo DEM downloaded >> kedarnath_status.txt
:have_dem
set CHK=--lon 79.0669 --lat 30.7352 --catchment-km2 67
echo. >> kedarnath_status.txt
echo ===== (a) FULLY AUTOMATIC ===== >> kedarnath_status.txt
python -m fanrisk run examples\kedarnath_2013\site_auto.toml > kedarnath_2_auto.txt 2>&1
if errorlevel 1 ( set WHY=the automatic run failed & set LOG=kedarnath_2_auto.txt & goto fail )
python tools\check_event.py examples\kedarnath_2013\site_auto.toml %CHK% >> kedarnath_status.txt 2>&1
echo. >> kedarnath_status.txt
echo ===== (b) ONE MANUAL APEX ===== >> kedarnath_status.txt
python -m fanrisk run examples\kedarnath_2013\site.toml > kedarnath_3_run.txt 2>&1
if errorlevel 1 ( set WHY=the run failed & set LOG=kedarnath_3_run.txt & goto fail )
python tools\check_event.py examples\kedarnath_2013\site.toml %CHK% >> kedarnath_status.txt 2>&1
echo. >> kedarnath_status.txt
echo ===== fan complexes and peak flows (manual-apex run) ===== >> kedarnath_status.txt
powershell -NoProfile -Command "Get-Content kedarnath_3_run.txt | Select-String -Pattern 'complex|feeder|peak|Melton|fan' | Select-Object -First 40" >> kedarnath_status.txt
python -m fanrisk viewer examples\kedarnath_2013\site.toml > kedarnath_5_viewer.txt 2>&1
echo [%date% %time%] STEP 1 DONE - viewer: runs\kedarnath_2013\viewer.html >> kedarnath_status.txt
type kedarnath_status.txt
pause
exit /b 0
:fail
echo STOPPED: %WHY% >> kedarnath_status.txt
powershell -NoProfile -Command "Get-Content '%LOG%' -Tail 15" >> kedarnath_status.txt
type kedarnath_status.txt
pause
exit /b 1
