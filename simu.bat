@echo off
REM
REM Batch script to prepare and launch simu.py. 
REM This needs m4 on windows (https://gnuwin32.sourceforge.net/packages/m4.htm).
REM Usage: simuit <template> <run>
REM
set Tmpl_Version=%1
set Run_Version=%2

IF "%Tmpl_Version%X" == "--helpX" GOTO HELP
IF "%Tmpl_Version%X" == "X" GOTO LAUNCH
echo "LAUNCH SIMULINAC with simuIn.yml made with Template=%Tmpl_Version% and Run number=%Run_Version%"
bash -c "yml/macros_%Tmpl_Version%.%Run_Version%.sh yml/tmpl_%Tmpl_Version%.yml simuIN.yml"
GOTO LAUNCH

ELSE GOTO:EOF

:HELP
echo "Batch script to prepare and launch simu.py. This needs m4 on windows."
echo "Usage: simuit [<template> <run> | --help]"
GOTO:EOF

:LAUNCH
REM launch simu.py
python simu.py --file trackerIN_smh27-3.1.5.work.yml
GOTO:EOF