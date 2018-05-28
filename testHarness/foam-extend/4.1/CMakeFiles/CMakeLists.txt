# /*-------------------------------------------------------------------------*\
#   =========                 |
#   \\      /  F ield         | foam-extend: Open Source CFD
#    \\    /   O peration     | Version:     4.1
#     \\  /    A nd           | Web:         http://www.foam-extend.org
#      \\/     M anipulation  | For copyright notice see file Copyright
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
#     CMakeLists.txt file for implementing a test harness for the compilation
#     and test of foam-extend-4.1 using Kitware CTest/CMake/CDash
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

PROJECT(foam-extend-4.1)

#-----------------------------------------------------------------------------
# Utility functions
#
# GetHostName(var)
#
# Set the variable named ${var} to the host hostname
#
function(GetHostName var)
    set( thisHostName "unknown")

    if(CMAKE_HOST_WIN32)
        execute_process(
            COMMAND         hostname
            OUTPUT_VARIABLE thisHostname
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    else()
        execute_process(
            COMMAND         hostname -f
            OUTPUT_VARIABLE thisHostname
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
    endif()

    set(${var} ${thisHostname} PARENT_SCOPE)
endfunction()

#
# GetGitStatus(status ecode)
#
# Get the git status for the local git repository of the source code
# The variable named ${status} will get the git status
# The variable named ${ecode} will get the command error code
#
function(GetGitStatus _status _ecode)
    set( git_status "unknown")
    set( git_ecode "1")

    execute_process(
      COMMAND           git status
      WORKING_DIRECTORY ${FOAM_ROOT}
      OUTPUT_VARIABLE   git_status
      RESULT_VARIABLE   git_ecode
      ERROR_QUIET
    )
    set(${_status} ${git_status} PARENT_SCOPE)
    set(${_ecode} ${git_ecode} PARENT_SCOPE)
endfunction()

#
# GetGitBranchName(var)
#
# Set the variable named ${var} to the git branch name of the source code
#
function(GetGitBranchName var)
    set( retValue "unknown")

    execute_process(
        COMMAND         git branch --no-color
        WORKING_DIRECTORY ${FOAM_ROOT}
        OUTPUT_VARIABLE listOfGitBranches
    )
    # Create list of strings
    string(REPLACE "\n" ";" listOfGitBranches ${listOfGitBranches})

    # Iterate over list, find the string beginning with "* "
    foreach(branch ${listOfGitBranches})
        string(REGEX MATCH "\\* .*$" matchString ${branch})
        string(LENGTH "${matchString}" lengthMatchString)
        if(lengthMatchString GREATER 0)
            # We have  match. Cleanup and set retValue
            string(REPLACE "* " "" retValue ${matchString})
        endif()
    endforeach()

    set(${var} ${retValue} PARENT_SCOPE)
endfunction()

#
# GetGitRevNumber(var)
#
# Set the variable named ${var} to the git revision number of the source code
#
function(GetGitRevNumber var)
    set( retValue "unknown")

    execute_process(
    COMMAND         git rev-parse --short=12 HEAD
        WORKING_DIRECTORY ${FOAM_ROOT}
        OUTPUT_VARIABLE git_rev_number
    )

    # Minimal check that we do have a valid string
    string(LENGTH "${git_rev_number}" lengthString)
    if(lengthString GREATER 0)
        set(retValue ${git_rev_number})
    endif()

    set(${var} ${retValue} PARENT_SCOPE)
endfunction()

# CleanUpStringForCDash(var)
#
# Cleanup the variable named ${value} for characters not supported
# for CDash identifiers. Return the modified value in retVar
#
function(CleanUpStringForCDash value retVar)
    string(REPLACE "/" "_" value "${value}")
    string(REPLACE " " "_" value ${value})
    set(${retVar} ${value} PARENT_SCOPE)
endfunction()

#-----------------------------------------------------------------------------
#
#
# Initialization of CTest specific variables
#
## Run ctest in parallel if environment variable WM_NCOMPPROCS is set
IF (NOT $ENV{WM_NCOMPPROCS} STREQUAL "")
    # Will run ctest in parallel over $WM_NCOMPPROCS processors
    set(CMAKE_CTEST_COMMAND ${CMAKE_CTEST_COMMAND} --parallel $ENV{WM_NCOMPPROCS})
    MESSAGE("Running tests in parallel using $ENV{WM_NCOMPPROCS} processors")
ENDIF ()

# Initialize the site name

IF (NOT $ENV{CDASH_SUBMIT_LOCAL_HOST_ID} STREQUAL "")
    # We can override the site name with the environment variable
    # $CDASH_SUBMIT_LOCAL_HOST_ID
    SET(
          SITENAME $ENV{CDASH_SUBMIT_LOCAL_HOST_ID}
          CACHE STRING "Name of the local site"
    )
ELSE ()
    # Grab the hostname FQN; will be used for the sitename
    GetHostName(SITENAME)

ENDIF()

MESSAGE("Initializing the name of this local site to:  ${SITENAME}")

SET(
    SITE ${SITENAME}
    CACHE STRING "Name of the local site"
)

#Grab the FOAM installation directory.
SET(
    FOAM_ROOT $ENV{WM_PROJECT_DIR}
    CACHE INTERNAL  "FOAM root directory."
)

# Construct the build name.
# No need to add $WM_PROJECT_VERSION to the name of the build,
# the test harness name should have taken care of that.
SET(
    BUILDNAME $ENV{WM_OPTIONS}
    CACHE STRING "Build ID"
)

# We allow overriding the git branch and revision information with some
# user-supplied information using the CDASH_SCM_INFO environment variable.
#
# Mercurial or other SCM users should be using this environment variable
# in order to provide branch and revision information for the buildname.
#
# Git users should use this environment variable in order to provide
# additionnal information instead of just the branch and revision number.
#
# Otherwise, the buildname will be appended with the following information:
# -git-branch=the_git_branch_name-git-rev=the_git_revision_number
SET(
    BUILDNAME_SCM_INFO $ENV{CDASH_SCM_INFO}
    CACHE STRING "SCM info for CDash buildname"
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
# Timeout for running all this: 20 minutes : 1200 seconds (for debug)
SET(
    DART_TESTING_TIMEOUT 1200
    CACHE STRING "Maximum time allowed (20 minutes) before CTest will kill the test."
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

if(NOT BUILDNAME_SCM_INFO STREQUAL "")
  SET(BUILDNAME "${BUILDNAME}-${BUILDNAME_SCM_INFO}")

elseif(GIT_FOUND)
    message("The git command was found: ${GIT_EXECUTABLE}")

    # Check if the source code is under a valid git repository
    GetGitStatus(GIT_STATUS GIT_ECODE)

    if(NOT GIT_ECODE)
        # We have a valid git repository.
        # Grab the branch and revision info. Add to the build name
        GetGitBranchName(GIT_BRANCH_NAME)
        message("Git branch: ${GIT_BRANCH_NAME}")

        GetGitRevNumber(GIT_REV_NUMBER)
        message("Git revision: ${GIT_REV_NUMBER}")

        SET(BUILDNAME "${BUILDNAME}-git-branch=${GIT_BRANCH_NAME}")
        SET(BUILDNAME "${BUILDNAME}-git-rev=${GIT_REV_NUMBER}")
    else()
        execute_process(
          COMMAND           hg id
          WORKING_DIRECTORY ${FOAM_ROOT}
          OUTPUT_VARIABLE   HG_STATUS
          RESULT_VARIABLE   HG_ECODE
          ERROR_QUIET
          )
        if(NOT HG_ECODE)
          # We have a valid git repository. Grab the branch and revision info.
          # Add to the build name
          message("No git-branch. Mercurial?")
          EXEC_PROGRAM(hg
            ARGS id --bookmarks
            OUTPUT_VARIABLE GIT_BRANCH_NAME
            )
          EXEC_PROGRAM(hg
            ARGS id --id
            OUTPUT_VARIABLE GIT_REV_NUMBER
            )
          EXEC_PROGRAM(hg
            ARGS log --template='git_{gitnode|short}' -l 1
            OUTPUT_VARIABLE GIT_REAL_REV_NUMBER
            )
          string(REPLACE " " "_" GIT_BRANCH_NAME ${GIT_BRANCH_NAME})
          string(REPLACE "+" ":modified" GIT_REV_NUMBER ${GIT_REV_NUMBER})
          SET(GIT_REV_NUMBER "${GIT_REV_NUMBER}_${GIT_REAL_REV_NUMBER}")
          message("Git branch (mercurial): ${GIT_BRANCH_NAME} Revision: ${GIT_REV_NUMBER}")

          SET(BUILDNAME "${BUILDNAME}-hg-branch=${GIT_BRANCH_NAME}")
          SET(BUILDNAME "${BUILDNAME}-hg-rev=${GIT_REV_NUMBER}")
        else()
          # Not a git or mercurial repository: no branch nor revision information available
          SET(BUILDNAME "${BUILDNAME}-git-branch=unknown")
          SET(BUILDNAME "${BUILDNAME}-git-rev=unknown")
        endif()
    endif()
else()
    # Git is not available: no branch nor revision information supplied
    SET(BUILDNAME "${BUILDNAME}-git-branch=unknown")
    SET(BUILDNAME "${BUILDNAME}-git-rev=unknown")
endif()

# Some last minute cleanup
# Seems like no '/' or ' 'are allowed in the BUILDNAME or in the SITE name
CleanUpStringForCDash(${BUILDNAME} BUILDNAME)
CleanUpStringForCDash(${SITE} SITE)

message("Build name: ${BUILDNAME}")
message("Site name: ${SITE}")

#
# Build section
#-----------------------------------------------------------------------------

# Compile FOAM, libs and apps
add_custom_target (foam-extend-$ENV{WM_PROJECT_VERSION} ALL
  ${FOAM_ROOT}/Allwmake
)

set_property(
  TARGET          foam-extend-$ENV{WM_PROJECT_VERSION}
  PROPERTY LABELS foam-extend-$ENV{WM_PROJECT_VERSION}
)

# Compile the FOAM unit tests located under applications/test
# This part will not be compiled and run by default.
# This would be a good candidate for a sub-project
add_custom_target (foam-extend-$ENV{WM_PROJECT_VERSION}_unitTests
  wmake all ${FOAM_ROOT}/applications/test
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

    # FOAM will run against this test suite:

    # Add the suite of FOAM tutorials
    #
    INCLUDE($ENV{FOAM_TEST_HARNESS_DIR}/CMakeFiles/FOAM_Tutorials.cmake)

    IF(RUN_FROM_ONE_TIMESTEP)
        # Modify the cases controlDict file in order to run for only one time step
        MESSAGE("${testRunTimeDirectory}: Modifying the controlDict files for running only one time step in directory: ${TEST_CASE_DIR}")
        if(CMAKE_HOST_WIN32)
          # Need to supply a bash shell to run the script under Windows
          EXECUTE_PROCESS(
            COMMAND bash -c "$ENV{FOAM_TEST_HARNESS_DIR}/scripts/prepareCasesForOneTimeStep.sh ${TEST_CASE_DIR}"
            WORKING_DIRECTORY .
          )
        else()
          EXECUTE_PROCESS(
            COMMAND $ENV{FOAM_TEST_HARNESS_DIR}/scripts/prepareCasesForOneTimeStep.sh ${TEST_CASE_DIR}
            WORKING_DIRECTORY .
          )
        endif()
    ENDIF(RUN_FROM_ONE_TIMESTEP)

ENDIF(BUILD_TESTING)

# That's it.
#
