@echo off
rem Maps only (no HEC-RAS needed): check + fans + risk + handoff for the Leh example.
cd /d "%~dp0"
python -m fanrisk check examples\leh_2010\site.toml > maps_1_check.txt 2>&1
python -m fanrisk run examples\leh_2010\site.toml > maps_2_run.txt 2>&1
echo exit %errorlevel% >> maps_2_run.txt
