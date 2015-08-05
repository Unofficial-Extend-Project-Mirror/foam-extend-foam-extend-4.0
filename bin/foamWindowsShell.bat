@echo off

set FOAM_HOME=%~dp0
call %FOAM_HOME%\setWindowsEnvironment.bat
mode 160,40
color 81
echo ---------------------------------
echo Command shell for foam-extend-3.1
echo ---------------------------------
echo/
cmd.exe
