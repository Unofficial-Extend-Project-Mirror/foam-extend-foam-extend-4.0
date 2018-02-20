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
#     verifyURLs.sh
#
# Description
#     Extract URLs from AllMake.stage* files, and report if the URLs are still
#     valid.
#
# Caveat: We are currently skipping the verification of swak4Foam, available
#         from either Mercurial or Subversion.
#
# Example:
#          tools/verifyURLs.sh AllMake.stage1
#
# Author:
#     Martin Beaudoin (2018)
#
#------------------------------------------------------------------------------
# Can run from anywhere

# Extract the list of URLs from the file specified from the command line

listOfURLs=`cat $1 | sed -ne 's/.*\(http[^"]*\).*/\1/p'  | \
	tr -d ')' | tr -d '\' | \
	grep -v licenses | grep -v SWAK_RELEASE_VERSION`

# Iterate over each url, checking if curl can access it
for url in $listOfURLs
do
	if [[ $url == http* ]] ;
	then
		printf "Verifying: $url :"

		# Using curl
		if curl --output /dev/null --silent --head --fail "$url"; then
			echo " OK"
		else
			echo " <----- CANNOT ACCESS!!!"
		fi
	fi
done

exit
