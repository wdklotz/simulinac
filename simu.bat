@echo off

if "%1%X" == "-?X" GOTO HELP
GOTO LAUNCH

:HELP
echo Windows batch script to prepare and launch simu.py. 
echo Usage: simu.bat [arg1, arg2,...]
GOTO EOF

:LAUNCH
REM launch simu.py
echo launching: python simu.py %1 %2 %3 %4 %5 %6 %7
python simu.py %1 %2 %3 %4 %5 %6 %7
GOTO EOF

:EOF