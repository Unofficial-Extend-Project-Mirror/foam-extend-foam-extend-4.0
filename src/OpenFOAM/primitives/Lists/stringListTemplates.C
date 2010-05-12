/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "labelList.H"

#include <sys/types.h>
#include <regex.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class StringList>
labelList findStrings(const string& regexp, const StringList& sl)
{
    labelList matches(sl.size());

    regex_t *preg = new regex_t;

    if (regcomp(preg, regexp.c_str(), REG_EXTENDED) != 0)
    {
        WarningIn("findStrings(const string& regexp, const stringList& sl)")
            << "Failed to compile regular expression " << regexp
            << endl;

        return matches;
    }

    size_t nmatch = 1;
    regmatch_t pmatch[1];

    label matchi = 0;
    forAll(sl, i)
    {
        if
        (
            regexec(preg, sl[i].c_str(), nmatch, pmatch, 0) == 0
         && (pmatch[0].rm_so == 0 && pmatch[0].rm_eo == label(sl[i].size()))
        )
        {
            matches[matchi++] = i;
        }
    }

    matches.setSize(matchi);

    regfree(preg);
    delete preg;

    return matches;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
