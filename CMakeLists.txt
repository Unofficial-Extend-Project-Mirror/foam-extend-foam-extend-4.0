# /*-------------------------------------------------------------------------*\
#   =========                 |
#   \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#    \\    /   O peration     |
#     \\  /    A nd           | Copyright held by original author
#      \\/     M anipulation  |
# -----------------------------------------------------------------------------
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
# Description
#     CMakeLists.txt file for implementing a test harness for the compilation
#     and test of OpenFOAM-1.6-ext using Kitware CTest./CMake/CDash
#
#     The results will be submitted to the CDash server identified by the file
#     CTestConfig.cmake
#
# Author
#     Martin Beaudoin, Hydro-Quebec, 2010. All rights reserved
#
#
# \*-------------------------------------------------------------------------*/

cmake_minimum_required (VERSION 2.8)

PROJECT(OpenFOAM_1.6-ext)

#-----------------------------------------------------------------------------
# Initialization of CTest specific variables

# Initialize the site name

IF (NOT $ENV{CDASH_SUBMIT_LOCAL_HOST_ID} STREQUAL "")
    # We can override the site name with the environment variable
    # $CDASH_SUBMIT_LOCAL_HOST_ID
    SET(
          SITENAME $ENV{CDASH_SUBMIT_LOCAL_HOST_ID}
	  CACHE STRING "Name of the local site"
    )
ELSE (NOT $ENV{CDASH_SUBMIT_LOCAL_HOST_ID} STREQUAL "")
    # Grab the hostname FQN; will be used for the sitename
    execute_process(
        COMMAND         hostname -f
        OUTPUT_VARIABLE SITENAME
    )
ENDIF  (NOT $ENV{CDASH_SUBMIT_LOCAL_HOST_ID} STREQUAL "")

MESSAGE("Initializing the name of this local site to:  ${SITENAME}")

SET(
    SITE ${SITENAME}
    CACHE STRING "Name of the local site"
)

#Grab the OpenFOAM installation directory.
SET(
    OF_ROOT $ENV{WM_PROJECT_DIR}
    CACHE INTERNAL  "OpenFOAM root directory."
)

# Construct the build name.
# No need to add $WM_PROJECT_VERSION to the name of the build,
# the test harness name should have taken care of that.
SET(
    BUILDNAME $ENV{WM_OPTIONS}
    CACHE STRING "Build ID"
)

# Find out the version of the compiler being used.
# Add this information to the buildname 
# This is for gcc or icc because they both support the -dumpversion option
EXEC_PROGRAM($ENV{WM_CC}
  ARGS -dumpversion
  OUTPUT_VARIABLE COMPILER_VERSION
)
SET(BUILDNAME "${BUILDNAME}-$ENV{WM_CC}${COMPILER_VERSION}")  
#
# We will support more compilers eventually.
# 

# Timeout for running every single test: 4 hours: 4 x 3600 seconds
#SET(
#    DART_TESTING_TIMEOUT 14400
#    CACHE STRING "Maximum time allowed (4 hours) before CTest will kill the test."
#)
# Timeout for running all this: 10 minutes : 600 seconds
SET(
    DART_TESTING_TIMEOUT 600
    CACHE STRING "Maximum time allowed (10 minutes) before CTest will kill the test."
)

SET(
    CMAKE_VERBOSE_MAKEFILE TRUE
)


# Update section
#-----------------------------------------------------------------------------
set (UPDATE_TYPE git)

#
# Using GIT as SCM
#
find_package(Git)
if(GIT_FOUND)
    message("Git was found: ${GIT_EXECUTABLE}")

    # Adding the name of the git branch to the build name
    EXEC_PROGRAM(git
      ARGS branch --no-color 2> /dev/null | grep '*'| awk '{print $2}'
      OUTPUT_VARIABLE GIT_BRANCH_NAME
    )
    message("Git branch: ${GIT_BRANCH_NAME}")
    if (GIT_BRANCH_NAME STREQUAL "")
       message("No git-branch. Mercurial?")
       EXEC_PROGRAM(hg
	 ARGS branch
         OUTPUT_VARIABLE GIT_BRANCH_NAME
       )
       message("Git branch (mercurial): ${GIT_BRANCH_NAME}")
    endif()
    SET(BUILDNAME "${BUILDNAME}-git-branch=${GIT_BRANCH_NAME}")  
endif()

# Some last minute cleanup
# Seems like no '/' are allowed in the BUILDNAME or in the SITE name
string(REPLACE "/" "_" BUILDNAME ${BUILDNAME})
string(REPLACE "/" "_" SITE ${SITE})


# Build section
#-----------------------------------------------------------------------------

# Compile OpenFOAM, libs and apps
add_custom_target (OpenFOAM-$ENV{WM_PROJECT_VERSION} ALL
  ${OF_ROOT}/Allwmake
)

set_property(
  TARGET          OpenFOAM-$ENV{WM_PROJECT_VERSION}
  PROPERTY LABELS OpenFOAM-$ENV{WM_PROJECT_VERSION}
)

# Compile the OpenFOAM unit tests located under applications/test
# This part will not be compiled and run by default.
# This would be a good candidate for a sub-project
add_custom_target (OpenFOAM-$ENV{WM_PROJECT_VERSION}_unitTests 
  wmake all ${OF_ROOT}/applications/test
)

# Test section
#-----------------------------------------------------------------------------

#Enable testing and dashboard
ENABLE_TESTING()
INCLUDE(CTest)

SET (CTEST_UPDATE_COMMAND ${GIT_EXECUTABLE})

SET(
    CTEST_CUSTOM_MAXIMUM_NUMBER_OF_ERRORS 1000
    CACHE INTERNAL "Max number of errors"
)
SET(
    CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS 1000
    CACHE INTERNAL "Max number of warnings"
)

IF(BUILD_TESTING)

    # Modify this variable if you want the full length test case simulations
    # Beware, this might take a long time to execute.
    # Otherwise, the default behaviour is to run each tutorial for 1 "timestep" 
    #SET(RUN_FROM_ONE_TIMESTEP 0)
    SET(RUN_FROM_ONE_TIMESTEP 1)

    IF(RUN_FROM_ONE_TIMESTEP)
        SET(testIdSuffix "_oneTimeStep")
    ENDIF(RUN_FROM_ONE_TIMESTEP)

    # OpenFOAM will run against this test suite:

    # Add the suite of OpenFOAM tutorials
    #
    INCLUDE($ENV{FOAM_TEST_HARNESS_DIR}/CMakeFiles/OpenFOAM_Tutorials.cmake)

    # Add a dummy test (/bin/true, just for debugging)
    ADD_TEST(
        OpenFOAM-$ENV{WM_PROJECT_VERSION}_Dummy_Test true
    )

    IF(RUN_FROM_ONE_TIMESTEP)
        # Modify the cases controlDict file in order to run for only one time step
        MESSAGE("${testRunTimeDirectory}: Modifying the controlDict files for running only one time step in directory: ${TEST_CASE_DIR}")
        EXECUTE_PROCESS(
            COMMAND $ENV{FOAM_TEST_HARNESS_DIR}/scripts/prepareCasesForOneTimeStep.sh ${TEST_CASE_DIR}
            WORKING_DIRECTORY .
            )
    ENDIF(RUN_FROM_ONE_TIMESTEP)

ENDIF(BUILD_TESTING)

# That's it.
#
