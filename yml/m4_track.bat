@echo off
REM
REM Batch script to launch m4 on W10. This needs cygwin.
REM
set Run_Version=%1

REM prefixes
set PFX1=C:\cygwin64\bin
set PFX2=/cygdrive/c/Users/wdklotz/SIMULINAC/yml

REM invoke bash and m4 from cygwin
%PFX1%\env PATH=/usr/bin bash -c "%PFX2%/macros_%Run_Version%.sh %PFX2%/tmpl_%Run_Version%.yml %PFX2%/trackIN.yml"
