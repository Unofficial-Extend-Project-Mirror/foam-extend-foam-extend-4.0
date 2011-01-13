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
# Description
#        CMakeLists.txt file for implementing a test harness for the 
#        Turbomachinery OSIG test cases
#
#        The results will be submitted to the CDash server identified by the file
#        CTestConfig.cmake
#
# Author
#        Martin Beaudoin, Hydro-Quebec, 2010. All rights reserved
#
#
# \*---------------------------------------------------------------------------*/

cmake_minimum_required (VERSION 2.8)

PROJECT(Turbomachinery-1.5-dev)

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

#Grab the Turbomachinery OSIG base directory
SET(
    TURBOMACHINERY_OSIG_ROOT $ENV{BREEDER_OSIG_TURBOMACHINERY_DIR}
    CACHE INTERNAL  "Turbomachinery OSIG root directory."
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
# This is for gcc 
IF(CMAKE_COMPILER_IS_GNUCC)
  EXEC_PROGRAM(${CMAKE_CXX_COMPILER}
    ARGS -dumpversion
    OUTPUT_VARIABLE COMPILER_VERSION
  )
  SET(BUILDNAME "${BUILDNAME}-gcc${COMPILER_VERSION}")  
ENDIF(CMAKE_COMPILER_IS_GNUCC)
#
# We will support more compilers eventually.
# 

# Timeout for running all this: 4 hours: 4 x 3600 seconds
SET(
    DART_TESTING_TIMEOUT 14400
    CACHE STRING "Maximum time allowed (4 hours) before CTest will kill the test."
)

SET(
    CMAKE_VERBOSE_MAKEFILE TRUE
)


# Update section
#-----------------------------------------------------------------------------
set (UPDATE_TYPE svn)


# Build section
#--------------------------------------------------------------------------------

# Compile the libs
add_custom_target (Turbo_Libs-src-$ENV{WM_PROJECT_VERSION} ALL
  ./Allwmake
  WORKING_DIRECTORY ${TURBOMACHINERY_OSIG_ROOT}/src
)
set_property(
  TARGET          Turbo_Libs-src-$ENV{WM_PROJECT_VERSION}
  PROPERTY LABELS Turbo_Libs-src-$ENV{WM_PROJECT_VERSION}
)

# rotorStatorInterfaces is not compiled by default
# MB 22/06/2010
add_custom_target (Turbo_Libs-rotorStatorInterfaces-$ENV{WM_PROJECT_VERSION}
  wmake libso
  WORKING_DIRECTORY ${TURBOMACHINERY_OSIG_ROOT}/rotorStatorInterfaces
)
set_property(
  TARGET          Turbo_Libs-rotorStatorInterfaces-$ENV{WM_PROJECT_VERSION}
  PROPERTY LABELS Turbo_Libs-rotorStatorInterfaces-$ENV{WM_PROJECT_VERSION}
)

# We need the Breeder's simpleFunctionObjects library
add_custom_target (Breeder_Libs-simpleFunctionObjects-$ENV{WM_PROJECT_VERSION} ALL
  wmake libso
  WORKING_DIRECTORY $ENV{BREEDER_15_DIR}/libraries/simpleFunctionObjects
)
set_property(
  TARGET          Breeder_Libs-simpleFunctionObjects-$ENV{WM_PROJECT_VERSION}
  PROPERTY LABELS Breeder_Libs-simpleFunctionObjects-$ENV{WM_PROJECT_VERSION}
)


# Compile the apps
add_custom_target (Turbo_Apps-Utilities-$ENV{WM_PROJECT_VERSION} ALL
  wmake all
  WORKING_DIRECTORY ${TURBOMACHINERY_OSIG_ROOT}/applications/utilities
)
set_property(
  TARGET          Turbo_Apps-Utilities-$ENV{WM_PROJECT_VERSION}
  PROPERTY LABELS Turbo_Apps-Utilities-$ENV{WM_PROJECT_VERSION}
)

add_custom_target (Turbo_Apps-Solvers-$ENV{WM_PROJECT_VERSION} ALL
  wmake all
  WORKING_DIRECTORY ${TURBOMACHINERY_OSIG_ROOT}/applications/solvers
)
set_property(
  TARGET          Turbo_Apps-Solvers-$ENV{WM_PROJECT_VERSION}
  PROPERTY LABELS Turbo_Apps-Solvers-$ENV{WM_PROJECT_VERSION}
)

add_custom_target (Turbo_Apps-SingleChannelPump-initalizeRotation-$ENV{WM_PROJECT_VERSION} ALL
  wmake 
  WORKING_DIRECTORY ${TURBOMACHINERY_OSIG_ROOT}/singleChannelPump/initializeRotation
)
set_property(
  TARGET          Turbo_Apps-SingleChannelPump-initalizeRotation-$ENV{WM_PROJECT_VERSION}
  PROPERTY LABELS Turbo_Apps-SingleChannelPump-initalizeRotation-$ENV{WM_PROJECT_VERSION}
)

# Adding target-level dependencies so the building is made in proper order
add_dependencies(
  Turbo_Apps-Utilities-$ENV{WM_PROJECT_VERSION} Turbo_Libs-src-$ENV{WM_PROJECT_VERSION}
)

# Test section
#-----------------------------------------------------------------------------

#Enable testing and dashboard
ENABLE_TESTING()
INCLUDE(CTest)

SET (CTEST_UPDATE_COMMAND "svn")

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
    # Otherwise, the default behaviour is to run each test case for 1 "timestep" 
    #SET(RUN_FROM_ONE_TIMESTEP 0)
    SET(RUN_FROM_ONE_TIMESTEP 1)

    IF(RUN_FROM_ONE_TIMESTEP)
        SET(testIdSuffix "_oneTimeStep")
    ENDIF(RUN_FROM_ONE_TIMESTEP)

    # Definition of the test suite:

    # Run the suite of Turbomachinery test cases
    #
    INCLUDE($ENV{FOAM_TEST_HARNESS_DIR}/CMakeFiles/TurbomachineryTestsuites_ercoftacConicalDiffuser.cmake)
    INCLUDE($ENV{FOAM_TEST_HARNESS_DIR}/CMakeFiles/TurbomachineryTestsuites_ercoftacCentrifugalPump.cmake)
    INCLUDE($ENV{FOAM_TEST_HARNESS_DIR}/CMakeFiles/TurbomachineryTestsuites_singleChannelPump.cmake)
    INCLUDE($ENV{FOAM_TEST_HARNESS_DIR}/CMakeFiles/TurbomachineryTestsuites_tutorials.cmake)

    # Add a dummy test (/bin/true, just for debugging)
    ADD_TEST(
         Turbo_Dummy_Test-$ENV{WM_PROJECT_VERSION} true
    )

    IF(RUN_FROM_ONE_TIMESTEP)
        # Modify the cases controlDict file in order to run for only one time step
        MESSAGE("${testRunTimeDirectory}: Modifying the controlDict files for running only one time step in directory: TurbomachineryTestSuites${testIdSuffix}")
        EXECUTE_PROCESS(
            COMMAND $ENV{FOAM_TEST_HARNESS_DIR}/scripts/prepareCasesForOneTimeStep.sh TurbomachineryTestSuites${testIdSuffix}
            WORKING_DIRECTORY .
            )
    ENDIF(RUN_FROM_ONE_TIMESTEP)

ENDIF(BUILD_TESTING)

# That's it.
#