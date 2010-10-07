#!/bin/bash
#------------------------------------------------------------------------------
# =========                 |
# \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
#  \\    /   O peration     |
#   \\  /    A nd           | Copyright held by original author
#    \\/     M anipulation  |
#------------------------------------------------------------------------------
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
