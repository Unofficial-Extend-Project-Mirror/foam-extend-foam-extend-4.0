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
    Writes regions of triSurface to separate files.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "triSurface.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("input file");
    argList args(argc, argv);

    fileName surfName(args.additionalArgs()[0]);

    Info<< "Reading surf from " << surfName << " ..." << nl << endl;

    fileName surfBase = surfName.lessExt();

    word extension = surfName.ext();

    triSurface surf(surfName);

    Info<< "Writing regions to separate files ..."
        << nl << endl;


    const geometricSurfacePatchList& patches = surf.patches();

    forAll(patches, patchI)
    {
        const geometricSurfacePatch& pp = patches[patchI];

        word patchName = pp.name();

        if (patchName.empty())
        {
            patchName = "patch" + Foam::name(patchI);
        }

        fileName outFile(surfBase + '_' + patchName + '.' + extension);

        Info<< "   Writing patch " << patchName << " to file " << outFile
            << endl;


        // Collect faces of region
        boolList includeMap(surf.size(), false);

        forAll(surf, faceI)
        {
            const labelledTri& f = surf[faceI];

            if (f.region() == patchI)
            {
                includeMap[faceI] = true;
            }
        }

        // Subset triSurface
        labelList pointMap;
        labelList faceMap;

        triSurface subSurf
        (
            surf.subsetMesh
            (
                includeMap,
                pointMap,
                faceMap
            )
        );

        subSurf.write(outFile);
    }


    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
