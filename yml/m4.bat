@echo off
set Run_Version=%1

REM prefixes
set PFX1=C:\cygwin64\bin
set PFX2=/cygdrive/c/Users/wdklotz/Desktop/SIMULINAC/yml

REM invoke cygwin bash and m4
%PFX1%\env PATH=/usr/bin bash -c "%PFX2%/macros_%Run_Version%.sh %PFX2%/tmpl_%Run_Version%.yml %PFX2%/simuIN.yml"
