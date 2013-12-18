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
