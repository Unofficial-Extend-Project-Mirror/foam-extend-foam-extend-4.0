@echo off

set FOAM_ETC_DIR=%~dp0
if #%FOAM_ETC_DIR:~-1%# == #\#  set FOAM_ETC_DIR=%FOAM_ETC_DIR:~0,-1%
call "%FOAM_ETC_DIR%\foamWindowsEnvironment.bat"
set FOAM_ETC_DIR=
mode 160,40
color 81
echo ---------------------------------
echo Command shell for foam-extend-4.0
echo ---------------------------------
echo/
cmd.exe
