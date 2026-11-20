@echo off
rem Step 9: one HTML web map of the results (risk map, HEC-RAS depth and hazard, 2010 flood, buildings and facilities hit).
rem Run after run_exposure.bat. Opens the page in your browser when done.
cd /d "%~dp0"
python -m fanrisk viewer examples\leh_2010\site.toml
if errorlevel 1 ( pause & exit /b 1 )
start "" "runs\leh_2010\viewer.html"
