/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Description
    Converts specified patches into subsets

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "triSurf.H"
#include "triFaceList.H"
#include "labelLongList.H"
#include "IFstream.H"
#include "OFstream.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();

    argList::validArgs.append("input surface file");
    argList::validArgs.append("output surface file");
    argList args(argc, argv);

    fileName inFileName(args.args()[1]);
    fileName outFileName(args.args()[2]);

    if( outFileName.ext() != "fms" )
        Warning << "The subsets cann only be saved in the .fms format" << endl;

    triSurf origSurf(inFileName);

    wordList patchNames(origSurf.patches().size());
    forAll(origSurf.patches(), patchI)
        patchNames[patchI] = origSurf.patches()[patchI].name();

    forAll(patchNames, patchI)
    {
        labelLongList subsetFacets;
        forAll(origSurf, triI)
        {
            if( origSurf[triI].region() == patchI )
                subsetFacets.append(triI);
        }

        const label subsetId = origSurf.addFacetSubset(patchNames[patchI]);

        forAll(subsetFacets, i)
            origSurf.addFacetToSubset(subsetId, subsetFacets[i]);
    }

    origSurf.writeSurface(outFileName);

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
