@echo off

rem =========== USER EDITABLE SETTINGS ===========
rem set MPI_ROOTDIR=C:\Programs\OpenMPI_v1.6.1-x64
rem set PARAVIEW_HOME=C:\Programs\ParaView-4.3.1
rem ==============================================

set FOAM_HOME=%~dp0
set FOAM_HOME=%FOAM_HOME:~0,-1%
set WM_PROJECT_DIR=%FOAM_HOME%

set MPI_BUFFER_SIZE=20000000

if defined PARAVIEW_HOME set PATH=%PARAVIEW_HOME%\bin;%PATH%
if defined MPI_ROOTDIR set PATH=%MPI_ROOTDIR%\bin;%PATH%
set PATH=%FOAM_HOME%\lib;%FOAM_HOME%\lib\openmpi-1.6.1;%FOAM_HOME%\bin;%PATH%
