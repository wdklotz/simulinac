@echo off
REM
REM Batch script to prepare and launch simu.py. This needs m4 on windows (https://gnuwin32.sourceforge.net/packages/m4.htm).
REM Usage: runit <template> <run>
REM
set Tmpl_Version=%1
set Run_Version=%2

IF "%Tmpl_Version%X" == "--helpX" GOTO HELP
IF "%Tmpl_Version%X" == "X" GOTO LAUNCH
REM echo "not LAUNCH"
bash -c "yml/macros_%Tmpl_Version%.%Run_Version%.sh yml/tmpl_%Tmpl_Version%.yml simuIN.yml"
echo "LAUNCH SIMULINAC with simuIn.yml made with Template=%Tmpl_Version% and Run number=%Run_Version%"
GOTO LAUNCH

ELSE GOTO:EOF

:HELP
echo "Batch script to prepare and launch simu.py. This needs cygwin."
echo "Usage: runit [tmpl_version run_version | --help]"
GOTO:EOF

:LAUNCH
REM launch simu.py
python simu.py --file simuIN.yml
GOTO:EOF