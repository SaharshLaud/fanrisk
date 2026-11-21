@echo off
rem Kedarnath 2013, step 2: HEC-RAS 2D on the Mandakini fan that holds the town, with the tool's own storm hydrograph.
rem Needs what the Leh runs needed (HEC-RAS 7.0.2, ras-commander 0.113, psutil, scipy, h5py). About 20-60 min.
rem RESUMES: the built model, computed runs and exported maps are kept and skipped.  run_kedarnath_hecras.bat fresh -> start again.
rem It stops at the first failure and prints the reason. Log: kedarnath_hecras_status.txt
setlocal
set FANRISK_MAPS=hdf
cd /d "%~dp0"
set S=examples\kedarnath_2013\site_hecras.toml
set H=runs\kedarnath_2013\hecras
set L=kedarnath_h
echo [%date% %time%] start > kedarnath_hecras_status.txt
if /i "%1"=="fresh" (
  if exist %H% rmdir /s /q %H%
  del /q %L%_*.txt 2>nul
  echo fresh start: earlier Kedarnath HEC-RAS model and runs deleted >> kedarnath_hecras_status.txt
)
python -c "import ras_commander, psutil, scipy, h5py" > %L%_0_preflight.txt 2>&1
if errorlevel 1 ( set WHY=Python libraries missing - run: pip install ras-commander==0.113.0 psutil scipy h5py & set LOG=%L%_0_preflight.txt & goto fail )

rem --- 1. handoff + model (kept when it exists)
if exist %H%\model\Terrain\Terrain.hdf if exist %H%\model\extended.ok ( echo model exists - kept >> kedarnath_hecras_status.txt & goto runs )
if exist %H% rmdir /s /q %H%
python -m fanrisk run %S% > %L%_1_run.txt 2>&1
if errorlevel 1 ( set WHY=the map run with the HEC-RAS handoff failed & set LOG=%L%_1_run.txt & goto fail )
powershell -NoProfile -Command "Get-Content %L%_1_run.txt | Select-String -Pattern 'inflow:|handoff|reach' | Select-Object -First 8" >> kedarnath_hecras_status.txt
python -m fanrisk hecras fit %S% --spacing 15 > %L%_2_fit.txt 2>&1
if errorlevel 1 ( set WHY=building the HEC-RAS model failed & set LOG=%L%_2_fit.txt & goto fail )
findstr /c:"terrain could not be built" %L%_2_fit.txt >nul && ( set WHY=the terrain was not built & set LOG=%L%_2_fit.txt & goto fail )
python -m fanrisk hecras extend-outflow %S% --metres 900 > %L%_3_extend.txt 2>&1
if errorlevel 1 ( set WHY=extend-outflow failed & set LOG=%L%_3_extend.txt & goto fail )
echo ok > %H%\model\extended.ok
echo model built (new UTM zone 44N, template re-fitted) >> kedarnath_hecras_status.txt
type %L%_3_extend.txt >> kedarnath_hecras_status.txt

:runs
for %%R in (R04 R05) do (
  call :one %%R
  if errorlevel 1 goto fail
)
python -m fanrisk viewer %S% > %L%_9_viewer.txt 2>&1
echo [%date% %time%] ALL DONE - viewer: runs\kedarnath_2013\viewer.html >> kedarnath_hecras_status.txt
type kedarnath_hecras_status.txt
pause
exit /b 0

:one
set R=%H%\runs\%1
if exist "%R%\maps\Depth (Max).tif" goto score
if exist %R%\project\leh2d.p02.hdf ( echo %1: computed run found - exporting maps only >> kedarnath_hecras_status.txt & goto export )
echo [%time%] %1: HEC-RAS running ... >> kedarnath_hecras_status.txt
echo [%time%] %1: HEC-RAS running (this is the long part) ...
python -m fanrisk hecras run %S% --runs %1 > %L%_4_%1.txt 2>&1
if errorlevel 1 ( set WHY=%1: the HEC-RAS run produced no maps & set LOG=%L%_4_%1.txt & exit /b 1 )
goto score
:export
python -m fanrisk hecras maps %S% --runs %1 > %L%_4_%1_maps.txt 2>&1
if errorlevel 1 ( set WHY=%1: exporting the Max maps failed & set LOG=%L%_4_%1_maps.txt & exit /b 1 )
:score
python -m fanrisk hecras score %S% --runs %1 > %L%_5_%1_score.txt 2>&1
python -m fanrisk compare %S% --run %1 > %L%_6_%1_compare.txt 2>&1
echo. >> kedarnath_hecras_status.txt
echo ===== %1 ===== >> kedarnath_hecras_status.txt
powershell -NoProfile -Command "Get-Content %L%_5_%1_score.txt | Select-String -Pattern 'peak|wet|flooded|depth|speed|velocity|hazard|AHP|share|AUC|edge|volume|points' | Select-Object -First 30" >> kedarnath_hecras_status.txt
echo [%time%] %1: done >> kedarnath_hecras_status.txt
exit /b 0

:fail
echo STOPPED: %WHY% >> kedarnath_hecras_status.txt
powershell -NoProfile -Command "Get-Content '%LOG%' -Tail 15" >> kedarnath_hecras_status.txt
type kedarnath_hecras_status.txt
pause
exit /b 1
