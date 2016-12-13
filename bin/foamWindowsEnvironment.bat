@echo off
if defined FOAM_ENV_SET goto :eof

rem =========== USER EDITABLE SETTINGS ===========
rem set MPI_ROOTDIR=C:\Programs\OpenMPI_v1.6.1-x64
rem set PARAVIEW_HOME=C:\Programs\ParaView-4.3.1
rem ==============================================

if defined PARAVIEW_HOME set PATH=%PARAVIEW_HOME%\bin;%PATH%
if defined MPI_ROOTDIR set PATH=%MPI_ROOTDIR%\bin;%PATH%
set PATH=%~dp0..\lib;%~dp0..\bin;%PATH%
set WM_PROJECT_DIR=%~dp0..
set MPI_BUFFER_SIZE=20000000
set FOAM_ENV_SET=1
