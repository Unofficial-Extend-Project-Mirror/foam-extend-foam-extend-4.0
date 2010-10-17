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
#     prepareCasesForTestHarness.sh
#
# Description
#     Prepare a working copy of the test cases for a test harness.
#
# Author:
#     Martin Beaudoin, Hydro-Quebec, 2010. All rights reserved
#
#------------------------------------------------------------------------------

#set -x

# Location of a working copy of the test cases
TEST_RUN_DIR=$1

#Additional shell functions
ADDITIONAL_SHELL_FUNCTIONS=$2

for AR in `find ${TEST_RUN_DIR} -name Allrun`
do
    # Make sure we are using /bin/bash as the activation shell
    # Replace the macro runApplication with runApplicationAndReportOnError
    mv ${AR} ${AR}.org
    sed \
	-e s/"#!\/bin\/sh"/"#!\/bin\/bash"/g \
	-e s/"runApplication "/"runApplicationAndReportOnError "/g \
        -e /RunFunctions/r${ADDITIONAL_SHELL_FUNCTIONS} \
	${AR}.org > ${AR}
	
    # Make sure the Allrun file is executable
    chmod +x ${AR}
done

echo "Done."
