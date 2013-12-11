#!/bin/bash
#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | foam-extend: Open Source CFD
#  \\    /   O peration     |
#   \\  /    A nd           | For copyright notice see file Copyright
#    \\/     M anipulation  |
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
# Script
#     addMissingAllrunFileToTutorial.sh
#
# Description
#     Add a default Allrun file to the tutorials that do not have one.
#     The test harness only run tutorials with an existing Allrun file
#
# Author:
#     Martin Beaudoin, Hydro-Quebec, 2010. All rights reserved
#
#------------------------------------------------------------------------------

#set -x

# Location of a working copy of the test cases
TEST_RUN_DIR=$1

DEFAULT_ALLRUN_FILE=$2

# Find all the test cases
for SF in `find ${TEST_RUN_DIR} -name system`
do
    caseDir=`dirname $SF`
    tutDir=`dirname $caseDir`

    # First, some sanity check:
    if [ $tutDir == $TEST_RUN_DIR ]; then
        echo "--> WARNING: Badly placed tutorial: $caseDir"
    fi

    # Check if an Allrun file is present
    if [ ! -e $caseDir/Allrun ]; then

        # Sometimes, a global Allrun file is located at the tutorial level
        if [ ! -e $tutDir/Allrun ]; then
            # No Allrun file for this case. Let's add one.
            echo "Adding a default Allrun file for tutorial : $caseDir"
            cp $DEFAULT_ALLRUN_FILE $caseDir/Allrun
        fi
    fi
done

echo "Done."
