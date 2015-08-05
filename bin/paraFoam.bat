@echo off

set CASE_DIR=%cd%
for /d %%d in (%CASE_DIR%) do set CASE_FILE=%%~nxd.foam

type nul >>%CASE_FILE%

set PARAVIEW_CMD=paraview --data="%CASE_FILE%"
echo Running %PARAVIEW_CMD% ...
%PARAVIEW_CMD%

del /f %CASE_FILE%
