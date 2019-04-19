#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | foam-extend: Open Source CFD
#  \\    /   O peration     | Version:     4.1
#   \\  /    A nd           | Web:         http://www.foam-extend.org
#    \\/     M anipulation  | For copyright notice see file Copyright
#------------------------------------------------------------------------------
# License
#     This file is part of foam-extend.
#
#     foam-extend is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     foam-extend is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.
#
# File
#     READMEBinaryPackage.txt
#
# Description
#     Installation and usage instructions for stand-alone Windows builds.
#
# Author:
#     Cesare Guardino, Alstom Power Ltd., (2015)
#
#------------------------------------------------------------------------------


*****************************************************************************************************************
*****************************************************************************************************************
***** IMPORTANT: THIS VERSION IS PROVIDED AS-IS AND IS NOT FULLY TESTED OR VALIDATED. PLEASE USE WITH CARE. *****
*****************************************************************************************************************
*****************************************************************************************************************


INSTRUCTIONS ON HOW TO INSTALL AND RUN THE WINDOWS VERSION OF FOAM-EXTEND
-------------------------------------------------------------------------
1) Unzip the package file to any suitable directory (with no whitespaces) on your computer.

2) Optionally, install OpenMPI and ParaView if you do not already have them.
   You can download these from:
       http://www.open-mpi.org/software/ompi/v1.6/downloads/OpenMPI_v1.6.1-1_win64.exe
       http://www.paraview.org/download
   It is strongly recommended to install these in directories with no white spaces.
   Once installed, create new environment variables called MPI_ROOTDIR and PARAVIEW_HOME
   to point to the installation directory where you installed them. This can be done by
   editing the user-editable settings in the environment configuration:
       call <PATH_TO_FOAM>\etc\foamWindowsEnvironment.bat
   where <PATH_TO_FOAM> is the full path of the directory where you unzipped the package.
   For example:
       rem =========== USER EDITABLE SETTINGS ===========
       set MPI_ROOTDIR=C:\Programs\OpenMPI_v1.6.1-x64
       set PARAVIEW_HOME=C:\Programs\ParaView-4.3.1
       rem ==============================================

3) Start a new CMD (DOS) prompt, and run the following command:
       call <PATH_TO_FOAM>\etc\foamWindowsEnvironment.bat
   where <PATH_TO_FOAM> is the full path of the directory where you unzipped the package.
   The foam environment is now configured correctly for use within this CMD prompt only.

   Alternatively, create a desktop shortcut to <PATH_TO_FOAM>\etc\foamWindowsShell.bat
   When you double-click this shortcut, a new CMD prompt is open with the foam environment automatically set.

4) From the CMD prompts opened in the previous step, you can now run the usual foam applications, for example:
       cd /d <PATH_TO_CASE>
       blockMesh
       sonicFoam

5) To run in parallel using using OpenMPI, run (for example):
       cd /d <PATH_TO_CASE>
       decomposePar
       mpirun -np 4 sonicFoam.exe -parallel
       reconstructPar

6) To post-process the results using ParaView, just run "parafoam" in the case directory:
       cd /d <PATH_TO_CASE>
       paraFoam

   This will create a temporary .foam file in the case directory and automatically launch ParaView.


NOTES
-----
1) This version does not fully support runtime selection of extended features. If necessary, you can
   explicitly load the required DLL using the "libs" function in system/controlDict, for example:
       libs ("liblduSolvers.dll");

2) The original foam-extend 4.0 source code is available from Git:
       git clone -b nextRelease http://git.code.sf.net/p/foam-extend/foam-extend-3.2 foam-extend-4.0
