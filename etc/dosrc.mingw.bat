@echo off

rem =========== USER EDITABLE SETTINGS ===========
set MPI_ROOTDIR=C:\Programs\OpenMPI_v1.6.1-x64
set PARAVIEW_HOME=C:\Programs\ParaView-4.3.1
rem ==============================================

set FOAM_HOME=%~dp0..
set WM_PROJECT_DIR=%FOAM_HOME%

set MPI_BUFFER_SIZE=20000000

if defined PARAVIEW_HOME set PATH=%PARAVIEW_HOME%\bin;%PATH%
if defined MPI_ROOTDIR set PATH=%MPI_ROOTDIR%\bin;%PATH%
set PATH=%FOAM_HOME%\lib\mingwGccDPOpt;%FOAM_HOME%\applications\bin\mingwGccDPOpt;%FOAM_HOME%\bin;%PATH%

set PATH=%FOAM_HOME%\..\site\4.0\lib\mingwGccDPOpt;%FOAM_HOME%\..\site\4.0\bin\mingwGccDPOpt;%PATH%
