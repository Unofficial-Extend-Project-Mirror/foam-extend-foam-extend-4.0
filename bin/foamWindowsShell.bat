@echo off

set FOAM_HOME=%~dp0
set FOAM_HOME=%FOAM_HOME:~0,-1%
call %FOAM_HOME%\setWindowsEnvironment.bat
mode 160,40
color 81
echo ---------------------------------
echo Command shell for foam-extend-3.1
echo ---------------------------------
echo/
cmd.exe
