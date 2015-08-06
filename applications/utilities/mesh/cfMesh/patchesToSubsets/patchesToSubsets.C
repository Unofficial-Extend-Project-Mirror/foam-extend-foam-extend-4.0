/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description
    Converts specified patches into subsets

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "Time.H"
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
