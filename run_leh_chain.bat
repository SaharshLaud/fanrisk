@echo off
rem Full Leh chain: maps, HEC-RAS model, test run, 6-run batch. Works from wherever this folder is.
cd /d "%~dp0"
set S=examples\leh_2010\site.toml
echo [%date% %time%] start > chain_status.txt
python -m fanrisk check %S% > chain_1_check.txt 2>&1
echo check exit %errorlevel% >> chain_status.txt
if errorlevel 1 goto end
python -m fanrisk run %S% > chain_2_run.txt 2>&1
echo run exit %errorlevel% >> chain_status.txt
if errorlevel 1 goto end
python -m fanrisk hecras fit %S% --spacing 15 > chain_3_fit.txt 2>&1
echo fit exit %errorlevel% >> chain_status.txt
if errorlevel 1 goto end
python -m fanrisk hecras extend-outflow %S% --metres 900 > chain_4_extend.txt 2>&1
echo extend exit %errorlevel% >> chain_status.txt
python -m fanrisk hecras check %S% > chain_5_hcheck.txt 2>&1
echo hecras check exit %errorlevel% >> chain_status.txt
python -m fanrisk hecras smoke %S% > chain_6_smoke.txt 2>&1
echo smoke exit %errorlevel% >> chain_status.txt
if errorlevel 1 goto end
python -m fanrisk hecras run %S% > chain_7_batch.txt 2>&1
echo batch exit %errorlevel% >> chain_status.txt
:end
echo [%date% %time%] done >> chain_status.txt
