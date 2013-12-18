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
#     prepareCasesForOneTimeStep.sh
#
# Description
#     Prepare a working copy of the tutorial test cases for running only one
#     timeStep
#
# Modifications made by:
#     Martin Beaudoin, Hydro-Quebec, 2010. All rights reserved
#
#------------------------------------------------------------------------------

# Location of a working copy of the test cases
TEST_RUN_DIR=$1

# Modify the controlDict files
pushd ${TEST_RUN_DIR} >& /dev/null

for CD in `find . -name "controlDict*"`
do
    mv ${CD} ${CD}.org
    sed \
    -e s/"\(startFrom[ \d9]*\)\([a-zA-Z]*\);"/"\1 latestTime;"/g \
    -e s/"\(stopAt[ \d9]*\)\([a-zA-Z]*\);"/"\1 nextWrite;"/g \
    -e s/"\(writeControl[ \d9]*\)\([a-zA-Z]*\);"/"\1 timeStep;"/g \
    -e s/"\(writeInterval[ \d9]*\)\([0-9a-zA-Z.-]*\);"/"\1 1;"/g \
    ${CD}.org > ${CD}
done
popd >& /dev/null

echo "Done."
