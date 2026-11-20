@echo off
cd /d "%~dp0"
echo [%date% %time%] start > selftest_status.txt
python tests\selftest.py --work "%~dp0_selftest" > selftest_log.txt 2>&1
echo selftest exit %errorlevel% >> selftest_status.txt
python -m fanrisk profiles >> selftest_log.txt 2>&1
echo [%date% %time%] done >> selftest_status.txt
