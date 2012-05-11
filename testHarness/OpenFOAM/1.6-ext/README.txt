# /*---------------------------------------------------------------------------*\
#   =========                 |
#   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#    \\    /   O peration     |
#     \\  /    A nd           | Copyright held by original author
#      \\/     M anipulation  |
# -------------------------------------------------------------------------------
# License
#     This file is part of OpenFOAM.
#
#     OpenFOAM is free software; you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation; either version 2 of the License, or (at your
#     option) any later version.
#
#     OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
#     ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#     FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#     for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with OpenFOAM; if not, write to the Free Software Foundation,
#     Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
#
# Description:
#              README file for the CMake/CTest/CDash test harness for OpenFOAM.
#
# Author:
#              Martin Beaudoin, Hydro-Quebec, 2010. All rights reserved.
# \*---------------------------------------------------------------------------*/

Warning #1: Make sure your OpenFOAM environment is properly initialized before
            running the test harness.

Warning #2: It is recommended to use cmake version 2.8.0 or newer for
            running the test harness.


1: Instructions for starting a local test harness on your machine:
------------------------------------------------------------------

a) You can set your local system identifier using the environment variable
   $CDASH_SUBMIT_LOCAL_HOST_ID.  Please try using a unique identifier like
   your machine's hostname. Otherwise, the fully qualified name of your system 
   will be used.


b) Install the CMakeLists.txt and CTestConfig.cmake files by executing the following commands:

cp ./CMakeFiles/CMakeLists.txt $WM_PROJECT_DIR
cp ./CMakeFiles/CTestConfig.cmake.openfoam-extend $WM_PROJECT_DIR/CTestConfig.cmake


c) Run the test harness and push your results on the CDash server on openfoam-extend

cd ./runDir; ./Allclean; ./Allrun_Experimental
  or
cd ./runDir; ./Allclean; ./Allrun_Nightly


2: Instructions for restarting your setup from scratch:
-------------------------------------------------------
cd ./runDir; ./Allclean


3: About submitting your results
--------------------------------

Both Allrun_Experimental and Allrun_Nightly scripts will submit their results
to the CDash server of your choice, as specified by the file
$WM_PROJECT_DIR/CTestConfig.cmake. 

To submit your results to the CDash server on openfoam-extend, just use
the file CTestConfig.cmake.openfoam-extend. 

If submitted to the CDash server on openfoam-extend, your results will be displayed here:
http://openfoam-extend.sourceforge.net/CDash/index.php?project=OpenFOAM-1.6-ext

