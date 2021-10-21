@echo off
REM
REM Batch script to launch m4 on W10. This needs cygwin.
REM
set File=%1
set Tmpl=%2
set Macro=%3

REM prefixes
set PFX1=C:\cygwin64\bin
set PFX2=/cygdrive/c/Users/wdklotz/SIMULINAC/

REM invoke bash and m4 from cygwin
%PFX1%\env PATH=/usr/bin bash -c "%PFX2%/%Macro% %PFX2%/%Tmpl% %PFX2%/%File%"
