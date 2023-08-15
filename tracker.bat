@echo off

if "%1%X" == "-?X" GOTO HELP
GOTO LAUNCH

:HELP
echo Windows batch script to prepare and launch tracker.py. 
echo Usage: tracker.bat [arg1, arg2,...]
GOTO EOF

:LAUNCH
REM launch tracker.py
echo launching: python tracker.py %*
python tracker.py %*
GOTO EOF

:EOF