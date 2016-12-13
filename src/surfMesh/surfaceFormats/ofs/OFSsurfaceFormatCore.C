/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "objectRegistry.H"
#include "OFSsurfaceFormatCore.H"
#include "clock.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::OFSsurfaceFormatCore::writeHeader
(
    Ostream& os,
    const pointField& pointLst,
    const UList<surfZone>& zoneLst
)
{
    // just emit some information until we get a nice IOobject
    IOobject::writeBanner(os)
        << "// foam Surface Format - written "
        << clock::dateTime().c_str() << nl
        << "// ~~~~~~~~~~~~~~~~~~~~~~~" << nl << nl
        << "// surfZones:" << nl;


    // treat a single zone as being unzoned
    if (zoneLst.size() <= 1)
    {
        os  << "0" << token::BEGIN_LIST << token::END_LIST << nl << nl;
    }
    else
    {
        os  << zoneLst.size() << nl << token::BEGIN_LIST << incrIndent << nl;

        forAll(zoneLst, zoneI)
        {
            zoneLst[zoneI].writeDict(os);
        }
        os  << decrIndent << token::END_LIST << nl << nl;
    }

    // Note: write with global point numbering

    IOobject::writeDivider(os)
        << "\n// points:" << nl << pointLst << nl;

    IOobject::writeDivider(os);
}


// ************************************************************************* //
