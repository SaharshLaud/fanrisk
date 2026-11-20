@echo off
rem Step 10: roughness sensitivity - HEC-RAS run R04 (1070 m3/s, 2 m channel burn) with Manning n = 0.04 and 0.08 (0.06 = step 7 burn 2 m run).
rem Needs: HEC-RAS 7.0.2 installed and opened once ('I Agree'); pip install ras-commander==0.113.0 psutil scipy
rem RESUMES: every finished stage is kept and skipped (model built, run computed, maps exported). About 3-15 min per terrain.
rem   run_roughness_test.bat fresh   -> delete the step-10 models and runs and start again
rem It stops at the first failure and prints the reason. Log: n_status.txt
setlocal
rem Max maps are made from the results HDF (Windows security blocks RAS Mapper's export helper on this PC)
set FANRISK_MAPS=hdf
cd /d "%~dp0"
echo [%date% %time%] start > n_status.txt
if /i "%1"=="fresh" (
  for %%D in (n004 n008) do (
    if exist runs\leh_2010_%%D\hecras\model rmdir /s /q runs\leh_2010_%%D\hecras\model
    if exist runs\leh_2010_%%D\hecras\runs rmdir /s /q runs\leh_2010_%%D\hecras\runs
  )
  del /q n00?_*.txt 2>nul
  echo fresh start: earlier step-7 models and runs deleted >> n_status.txt
)

python -c "import ras_commander, psutil, scipy, h5py" > n_preflight.txt 2>&1
if errorlevel 1 ( set WHY=Python libraries missing - run: pip install ras-commander==0.113.0 psutil scipy h5py & set LOG=n_preflight.txt & goto fail )

for %%D in (n004 n008) do (
  call :one %%D
  if errorlevel 1 goto fail
)
python tools\case_summary.py runs\leh_2010_n004 runs\leh_2010_burn2 runs\leh_2010_n008 >> n_status.txt 2>&1
echo [%date% %time%] ALL DONE >> n_status.txt
type n_status.txt
pause
exit /b 0

:one
set S=examples\leh_2010\experiments\site_%1.toml
set H=runs\leh_2010_%1\hecras
set R=%H%\runs\R04
set L=%1
rem --- 1. model (skipped when it exists: the handoff and the model are rebuilt only with 'fresh')
if exist %H%\model\Terrain\Terrain.hdf if exist %H%\model\leh2d.prj ( echo %1: model exists - kept >> n_status.txt & goto run )
if exist %H%\model rmdir /s /q %H%\model
if exist %H%\runs rmdir /s /q %H%\runs
python -m fanrisk run %S% > %L%_1_run.txt 2>&1
if errorlevel 1 ( set WHY=%1: map run failed & set LOG=%L%_1_run.txt & exit /b 1 )
python -m fanrisk hecras fit %S% --spacing 15 > %L%_2_fit.txt 2>&1
if errorlevel 1 ( set WHY=%1: building the HEC-RAS model failed & set LOG=%L%_2_fit.txt & exit /b 1 )
findstr /c:"terrain could not be built" %L%_2_fit.txt >nul && ( set WHY=%1: the terrain was not built & set LOG=%L%_2_fit.txt & exit /b 1 )
python -m fanrisk hecras extend-outflow %S% --metres 900 > %L%_3_extend.txt 2>&1
if errorlevel 1 ( set WHY=%1: extend-outflow failed & set LOG=%L%_3_extend.txt & exit /b 1 )
echo %1: model built >> n_status.txt
:run
rem --- 2. compute (skipped when results exist), 3. export Max maps (skipped when they exist), 4. read back + score
if exist "%R%\maps\Depth (Max).vrt" goto maps_ok
if exist "%R%\maps\Depth (Max).tif" goto maps_ok
if exist %R%\project\leh2d.p02.hdf ( echo %1: computed run found - exporting maps only >> n_status.txt & goto export )
if not "%1"=="0" goto compute
python -m fanrisk hecras smoke %S% > %L%_4_smoke.txt 2>&1
if errorlevel 1 ( set WHY=HEC-RAS does not run from Python [1-minute test] & set LOG=%L%_4_smoke.txt & exit /b 1 )
echo burn 0: HEC-RAS test run OK >> n_status.txt
:compute
echo [%time%] %1: HEC-RAS running ... >> n_status.txt
echo [%time%] %1: HEC-RAS running ...
python -m fanrisk hecras run %S% --runs R04 > %L%_5_hecras.txt 2>&1
if errorlevel 1 ( set WHY=%1: the HEC-RAS run produced no maps & set LOG=%L%_5_hecras.txt & exit /b 1 )
goto validate
:export
python -m fanrisk hecras maps %S% --runs R04 > %L%_5_maps.txt 2>&1
if errorlevel 1 ( set WHY=%1: exporting the Max maps failed & set LOG=%L%_5_maps.txt & exit /b 1 )
:maps_ok
if "%1"=="0" if exist "runs\leh_2010\hecras\runs\R04\maps\Depth (Max).vrt" (
  echo burn 0: check of the HDF maps against Akhtar's RAS Mapper maps of the same run: >> n_status.txt
  python -c "from fanrisk.ras import hdfmaps; hdfmaps.compare(r'%R%\maps', r'runs\leh_2010\hecras\runs\R04\maps')" >> n_status.txt 2>&1
)
python -m fanrisk hecras score %S% --runs R04 > %L%_5_score.txt 2>&1
:validate
python -m fanrisk validate %S% --maps %R%\maps > %L%_6_validate.txt 2>&1
if errorlevel 1 ( set WHY=%1: scoring failed & set LOG=%L%_6_validate.txt & exit /b 1 )
echo [%time%] %1: done >> n_status.txt
exit /b 0

:fail
echo. >> n_status.txt
echo STOPPED: %WHY% >> n_status.txt
echo ---- last lines of %LOG% ---- >> n_status.txt
powershell -NoProfile -Command "Get-Content '%LOG%' -Tail 12" >> n_status.txt
type n_status.txt
pause
exit /b 1
