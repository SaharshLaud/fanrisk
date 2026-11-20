@echo off
rem Step 7: HEC-RAS run R04 (1070 m3/s, clear water) on three terrains: plain DEM (0 m), channel burned 2 m and 4 m.
rem Needs: HEC-RAS 7.0.2 installed and opened once ('I Agree'); pip install ras-commander==0.113.0 psutil scipy
rem RESUMES: every finished stage is kept and skipped (model built, run computed, maps exported). About 3-15 min per terrain.
rem   run_burn_test.bat fresh   -> delete the step-7 models and runs and start again
rem It stops at the first failure and prints the reason. Log: burn_status.txt
setlocal
rem Max maps are made from the results HDF (Windows security blocks RAS Mapper's export helper on this PC)
set FANRISK_MAPS=hdf
cd /d "%~dp0"
echo [%date% %time%] start > burn_status.txt
if /i "%1"=="fresh" (
  for %%D in (0 2 4) do (
    if exist runs\leh_2010_burn%%D\hecras\model rmdir /s /q runs\leh_2010_burn%%D\hecras\model
    if exist runs\leh_2010_burn%%D\hecras\runs rmdir /s /q runs\leh_2010_burn%%D\hecras\runs
  )
  del /q burn?_*.txt 2>nul
  echo fresh start: earlier step-7 models and runs deleted >> burn_status.txt
)

python -c "import ras_commander, psutil, scipy, h5py" > burn_preflight.txt 2>&1
if errorlevel 1 ( set WHY=Python libraries missing - run: pip install ras-commander==0.113.0 psutil scipy h5py & set LOG=burn_preflight.txt & goto fail )

for %%D in (0 2 4) do (
  call :one %%D
  if errorlevel 1 goto fail
)
python tools\burn_summary.py >> burn_status.txt 2>&1
echo [%date% %time%] ALL DONE >> burn_status.txt
type burn_status.txt
pause
exit /b 0

:one
set S=examples\leh_2010\experiments\site_burn%1.toml
set H=runs\leh_2010_burn%1\hecras
set R=%H%\runs\R04
set L=burn%1
rem --- 1. model (skipped when it exists: the handoff and the model are rebuilt only with 'fresh')
if exist %H%\model\Terrain\Terrain.hdf if exist %H%\model\leh2d.prj ( echo burn %1: model exists - kept >> burn_status.txt & goto run )
if exist %H%\model rmdir /s /q %H%\model
if exist %H%\runs rmdir /s /q %H%\runs
python -m fanrisk run %S% > %L%_1_run.txt 2>&1
if errorlevel 1 ( set WHY=burn %1: map run failed & set LOG=%L%_1_run.txt & exit /b 1 )
python -m fanrisk hecras fit %S% --spacing 15 > %L%_2_fit.txt 2>&1
if errorlevel 1 ( set WHY=burn %1: building the HEC-RAS model failed & set LOG=%L%_2_fit.txt & exit /b 1 )
findstr /c:"terrain could not be built" %L%_2_fit.txt >nul && ( set WHY=burn %1: the terrain was not built & set LOG=%L%_2_fit.txt & exit /b 1 )
python -m fanrisk hecras extend-outflow %S% --metres 900 > %L%_3_extend.txt 2>&1
if errorlevel 1 ( set WHY=burn %1: extend-outflow failed & set LOG=%L%_3_extend.txt & exit /b 1 )
echo burn %1: model built >> burn_status.txt
:run
rem --- 2. compute (skipped when results exist), 3. export Max maps (skipped when they exist), 4. read back + score
if exist "%R%\maps\Depth (Max).vrt" goto maps_ok
if exist "%R%\maps\Depth (Max).tif" goto maps_ok
if exist %R%\project\leh2d.p02.hdf ( echo burn %1: computed run found - exporting maps only >> burn_status.txt & goto export )
if not "%1"=="0" goto compute
python -m fanrisk hecras smoke %S% > %L%_4_smoke.txt 2>&1
if errorlevel 1 ( set WHY=HEC-RAS does not run from Python [1-minute test] & set LOG=%L%_4_smoke.txt & exit /b 1 )
echo burn 0: HEC-RAS test run OK >> burn_status.txt
:compute
echo [%time%] burn %1: HEC-RAS running ... >> burn_status.txt
echo [%time%] burn %1: HEC-RAS running ...
python -m fanrisk hecras run %S% --runs R04 > %L%_5_hecras.txt 2>&1
if errorlevel 1 ( set WHY=burn %1: the HEC-RAS run produced no maps & set LOG=%L%_5_hecras.txt & exit /b 1 )
goto validate
:export
python -m fanrisk hecras maps %S% --runs R04 > %L%_5_maps.txt 2>&1
if errorlevel 1 ( set WHY=burn %1: exporting the Max maps failed & set LOG=%L%_5_maps.txt & exit /b 1 )
:maps_ok
if "%1"=="0" if exist "runs\leh_2010\hecras\runs\R04\maps\Depth (Max).vrt" (
  echo burn 0: check of the HDF maps against Akhtar's RAS Mapper maps of the same run: >> burn_status.txt
  python -c "from fanrisk.ras import hdfmaps; hdfmaps.compare(r'%R%\maps', r'runs\leh_2010\hecras\runs\R04\maps')" >> burn_status.txt 2>&1
)
python -m fanrisk hecras score %S% --runs R04 > %L%_5_score.txt 2>&1
:validate
python -m fanrisk validate %S% --maps %R%\maps > %L%_6_validate.txt 2>&1
if errorlevel 1 ( set WHY=burn %1: scoring failed & set LOG=%L%_6_validate.txt & exit /b 1 )
echo [%time%] burn %1: done >> burn_status.txt
exit /b 0

:fail
echo. >> burn_status.txt
echo STOPPED: %WHY% >> burn_status.txt
echo ---- last lines of %LOG% ---- >> burn_status.txt
powershell -NoProfile -Command "Get-Content '%LOG%' -Tail 12" >> burn_status.txt
type burn_status.txt
pause
exit /b 1
