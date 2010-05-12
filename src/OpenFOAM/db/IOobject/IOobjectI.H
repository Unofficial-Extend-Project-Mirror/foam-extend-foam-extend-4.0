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

#include "foamVersion.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Stream>
inline void Foam::IOobject::writeBanner(Stream& os, bool noHint)
{
    static bool spacesSet = false;
    static char spaces[40];
    static char spacesRev[40];

    if (!spacesSet)
    {
        memset(spaces, ' ', 40);
        spaces[38 - strlen(Foam::FOAMversion)] = '\0';

        memset(spacesRev, ' ', 40);
        spacesRev[38 - strlen(Foam::FOAMDevRevisionNumber)] = '\0';
        
        spacesSet = true;
    }

    if (noHint)
    {
        os  <<
            "/*---------------------------------------------------------------------------*\\\n";
    }
    else
    {
        os  <<
            "/*--------------------------------*- C++ -*----------------------------------*\\\n";
    }

    os  <<
        "| =========                 |                                                 |\n"
        "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n"
        "|  \\\\    /   O peration     | Version:  " << FOAMversion           << spaces << "|\n"
        "|   \\\\  /    A nd           | Revision: " << FOAMDevRevisionNumber << spacesRev << "|\n"
        "|    \\\\/     M anipulation  | Web:      http://www.OpenFOAM.org               |\n"
        "\\*---------------------------------------------------------------------------*/\n";
}


template<class Stream>
inline void Foam::IOobject::writeDivider(Stream& os)
{
    os  <<
        "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n";
}


// ************************************************************************* //
