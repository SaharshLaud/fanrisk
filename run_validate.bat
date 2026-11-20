@echo off
rem Leh 2010 validation: maps, scores against the 2010 footprint, AHP-vs-HEC-RAS comparison (uses the R04 maps in runs\leh_2010\hecras\runs\R04).
cd /d "%~dp0"
set S=examples\leh_2010\site.toml
python -m fanrisk run %S% > validate_1_run.txt 2>&1
python -m fanrisk validate %S% --maps runs\leh_2010\hecras\runs\R04\maps > validate_2_scores.txt 2>&1
python -m fanrisk compare %S% --run R04 > validate_3_compare.txt 2>&1
echo done >> validate_3_compare.txt
