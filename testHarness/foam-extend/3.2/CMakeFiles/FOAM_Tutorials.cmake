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
# Description
#        CMake/CTest script for running the FOAM tutorials
#
# Author
#        Martin Beaudoin, Hydro-Quebec, 2010. All rights reserved
#
# \*-------------------------------------------------------------------------*/

# Take care of tests specific variables
IF(NOT DEFINED testIdSuffix)
    SET(testIdSuffix "_full")
ENDIF(NOT DEFINED testIdSuffix)

IF(NOT DEFINED testRunTimeDirectory)
    SET(testRunTimeDirectory "tutorialsTestSuites${testIdSuffix}")
ENDIF(NOT DEFINED  testRunTimeDirectory)

# Use the current directory for running the test cases
SET (TEST_CASE_DIR $ENV{PWD}/${testRunTimeDirectory}
    CACHE INTERNAL "FOAM test case directory."
    )

# Create the runTime directory for the test cases
MESSAGE("Removing the old test directory: ${TEST_CASE_DIR}")
file(REMOVE_RECURSE ${TEST_CASE_DIR})
MESSAGE("Creation of the new test directory: ${TEST_CASE_DIR}")
file(COPY $ENV{FOAM_TUTORIALS}/ DESTINATION ${TEST_CASE_DIR})

# Add a default Allrun file for the cases that don't have one.
# The test harness relies on the presence of an Allrun file for
# running the case
MESSAGE("${testRunTimeDirectory}: Checking for missing Allrun file in tutorials")
if(CMAKE_HOST_WIN32)
  # Need to supply a bash shell to run the script under Windows
  EXECUTE_PROCESS(
    COMMAND bash -c "$ENV{FOAM_TEST_HARNESS_DIR}/scripts/addMissingAllrunFileToTutorial.sh ${TEST_CASE_DIR} $ENV{FOAM_TEST_HARNESS_DIR}/scripts/Allrun.default"
    WORKING_DIRECTORY .
    )
else()
  EXECUTE_PROCESS(
    COMMAND $ENV{FOAM_TEST_HARNESS_DIR}/scripts/addMissingAllrunFileToTutorial.sh ${TEST_CASE_DIR} $ENV{FOAM_TEST_HARNESS_DIR}/scripts/Allrun.default
    WORKING_DIRECTORY .
    )
endif()
# Iterate over each tutorial case:
# We are looking for tutorial cases with an Allrun file.
# If this file is present, (and it should), we add this case to the list of cases to run.
#

#First, add a global cleanup of the cases
#This will always run and complete first, even in paralel
SET(cleanCasesTestId "Allclean_cases${testIdSuffix}")
ADD_TEST(${cleanCasesTestId} bash -c "cd ${TEST_CASE_DIR}; ./Allclean")

# Add a dummy test (/bin/true, just for debugging)
SET(dummyTestId "foam-extend-$ENV{WM_PROJECT_VERSION}_Dummy_Test")
ADD_TEST(${dummyTestId} true)

# Add a dependency on the global clean-up target, even for the dummy test
SET_TESTS_PROPERTIES(${dummyTestId} PROPERTIES DEPENDS ${cleanCasesTestId})


# Next, recurse through the test cases root directory,
# find all the Allrun files, and add them as a new CTest test case
FILE(GLOB_RECURSE listofCasesWithAllrun ${TEST_CASE_DIR}/Allrun)
LIST(SORT listofCasesWithAllrun)

FOREACH(caseWithAllrun ${listofCasesWithAllrun})
  #Grab the name of the directory containing the file Allrun
  get_filename_component(thisCasePath ${caseWithAllrun} PATH)

  # We need to skip the global Allrun file
  IF(NOT ${thisCasePath} STREQUAL ${TEST_CASE_DIR})
    MESSAGE("Found Allrun file in directory: ${thisCasePath}")

    # Grab the parent name of the case directory
    string(REPLACE ${TEST_CASE_DIR}/ "" caseParentPath ${thisCasePath})

    # Construct the testId
    string(REPLACE "/" "_" testId ${caseParentPath})
    SET(testId ${testId}_Allrun${testIdSuffix})

    # Add the test to the test harness
    MESSAGE("    Adding test: ${testId}")
    ADD_TEST(${testId} bash -c "cd ${thisCasePath}; ./Allrun")

    # We extract a label name from the top level tutorial directories
    # (eg: basic, incompressible, immersedBoundary, etc). We will use this
    # label in order to categorize the various test cases under a more 'generic
    # topic', so we can for instance limit the testharness to the
    # 'incompressible' test cases, etc., simply by using the ctest -L command.
    #
    # ctest --print-labels will print the list of all available labels.
    #
    string(REPLACE ${TEST_CASE_DIR}/ ""  tmpLabel ${caseWithAllrun})
    string(FIND ${tmpLabel} "/" indexFirstSlash)
    string(SUBSTRING ${tmpLabel} 0 ${indexFirstSlash} testLabel)

    # Add a dependency on the global clean-up target
    # When running in parallel, you need to wait for the cleanup to finish first
    message("    This test case will have the following labels: Tutorials, ${testLabel}")
    SET_TESTS_PROPERTIES(${testId} PROPERTIES DEPENDS ${cleanCasesTestId} LABELS "Tutorials;${testLabel}")

    # Use this following entry instead for testing purpose
    #ADD_TEST(${testId} bash -c "cd ${thisCasePath}; true")

  ENDIF()
ENDFOREACH(caseWithAllrun)

# Modify the cases Allrun files to incorporate additional shell functions
MESSAGE("${testRunTimeDirectory}: Modifying the Allrun files for additional shell functions in directory: ${TEST_CASE_DIR}")
if(CMAKE_HOST_WIN32)
  # Need to supply a bash shell to run the script under Windows
  EXECUTE_PROCESS(
    COMMAND bash -c "$ENV{FOAM_TEST_HARNESS_DIR}/scripts/prepareCasesForTestHarness.sh ${TEST_CASE_DIR} $ENV{FOAM_TEST_HARNESS_DIR}/scripts/AdditionalRunFunctions"
    WORKING_DIRECTORY .
    )
else()
  EXECUTE_PROCESS(
    COMMAND $ENV{FOAM_TEST_HARNESS_DIR}/scripts/prepareCasesForTestHarness.sh ${TEST_CASE_DIR} $ENV{FOAM_TEST_HARNESS_DIR}/scripts/AdditionalRunFunctions
    WORKING_DIRECTORY .
    )
endif()

# Configure the various ctest -S Dashboard drivers

# Driver for the tutorials
configure_file(
  "$ENV{FOAM_TEST_HARNESS_DIR}/CMakeFiles/Dashboard_Tutorials.cmake.in"
  "$ENV{FOAM_TEST_HARNESS_DIR}/runDir/Dashboard_Tutorials.cmake"
  @ONLY)


# That's it.
