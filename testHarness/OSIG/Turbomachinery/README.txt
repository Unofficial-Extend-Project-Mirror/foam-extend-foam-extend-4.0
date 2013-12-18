# /*-------------------------------------------------------------------------*\
#   =========                 |
#   \\      /  F ield         | foam-extend: Open Source CFD
#    \\    /   O peration     |
#     \\  /    A nd           | For copyright notice see file Copyright
#      \\/     M anipulation  |
# -----------------------------------------------------------------------------
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
# Description:
#              README file for the CMake/CTest/CDash test harness for
#              the TurboMachinery OSIG
#
# Author:
#              Martin Beaudoin, Hydro-Quebec, 2010. All rights reserved.
#
# \*-------------------------------------------------------------------------*/


Warning #1: Make sure your FOAM environment is properly initialized before
            running the test harness.

Warning #2: Make sure to create an environment variable called $BREEDER_15_DIR
            This environment variable must point to the Breeder_1.5 directory
            of your local working copy of the openfoam-extend Subversion repository

            Something like: export BREEDER_15_DIR=/someAbsolutePath/Breeder_1.5

Warning #3: It is recommended to use cmake version 2.8.0 or newer for
            running the test harness.



1: Instructions for starting a local test harness on your machine:
------------------------------------------------------------------

a) You can set your local system identifier using the environment variable
   $CDASH_SUBMIT_LOCAL_HOST_ID.  Please try using a unique identifier like
   your machine's hostname. Otherwise, the fully qualified name of your system
   will be used.


b) Install the CMakeLists.txt and CTestConfig.cmake files by executing the following commands:

cp ./CMakeFiles/CMakeLists.txt $BREEDER_15_DIR/OSIG/TurboMachinery/CMakeLists.txt
cp ./CMakeFiles/CTestConfig.cmake.openfoam-extend $BREEDER_15_DIR/OSIG/TurboMachinery/CTestConfig.cmake


c) Run the test harness and push your results on the CDash server on openfoam-extend

cd ./runDir; ./Allrun_Experimental
  or
cd ./runDir; ./Allrun_Nightly


2: Instructions for restarting your setup from scratch:
-------------------------------------------------------
cd ./runDir; ./Allclean


3: About submitting your results
--------------------------------

Both Allrun_Experimental and Allrun_Nightly scripts will submit their results
to the CDash server of your choice, as specified by the file
$BREEDER_15_DIR/OSIG/TurboMachinery/CTestConfig.cmake

To submit your results to the CDash server on openfoam-extend, just use
the file CTestConfig.cmake.openfoam-extend.

If submitted to the CDash server on openfoam-extend, your results will be displayed here:
http://openfoam-extend.sourceforge.net/CDash/index.php?project=Turbomachinery-1.5-dev

