@echo off
REM
REM Batch script to prepare and launch simu.py. This needs cygwin.
REM Usage: runit <run_version>
REM
set Run_Version=%1

IF "%Run_Version%X" == "--helpX" GOTO HELP
IF "%Run_Version%X" == "X" GOTO LAUNCH
REM echo "not LAUNCH"
REM prefixes
set PFX1=C:\cygwin64\bin
set PFX2=/cygdrive/c/Users/wdklotz/SIMULINAC/yml
REM invoke bash and m4 from cygwin
%PFX1%\env PATH=/usr/bin bash -c "%PFX2%/macros_%Run_Version%.sh %PFX2%/tmpl_%Run_Version%.yml %PFX2%/simuIN.yml"
echo "LAUNCH SIMULINAC with Run_Version=%Run_Version%"
GOTO LAUNCH

ELSE GOTO:EOF

:HELP
echo "Batch script to prepare and launch simu.py. This needs cygwin."
echo "Usage: runit [run_version | --help]"
GOTO:EOF

:LAUNCH
REM launch simu.py
python simu.py yml\simuIN.yml
GOTO:EOF