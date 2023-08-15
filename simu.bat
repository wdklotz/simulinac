@echo off

if "%1%X" == "-?X" GOTO HELP
GOTO LAUNCH

:HELP
echo Windows batch script to prepare and launch simu.py. 
echo Usage: simu.bat [arg1, arg2,...]
GOTO EOF

:LAUNCH
REM launch simu.py
echo launching: python simu.py %*
python simu.py %*
GOTO EOF

:EOF