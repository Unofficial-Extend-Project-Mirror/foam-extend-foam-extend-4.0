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

Application
    tetDecomposition.C

Author
    Hrvoje Jasak, resurrected to life by Sandeep Menon.

Description
    Decompose a given polyhedral mesh using tetDecomposition.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "tetPolyMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    tetPolyMesh getTets(mesh);

    runTime++;

    Info<< "Creating decomposed tet mesh. Number of tets: "
        << getTets.nTets()
        << endl;

    polyMesh newMesh
    (
        IOobject
        (
            fvMesh::defaultRegion,
            runTime.timeName(),
            runTime
        ),
        xferCopy(getTets.points()()),
        getTets.tetCells(),
        getTets.boundary().boundaryTriFaces(),
        mesh.boundaryMesh().names(),
        mesh.boundaryMesh().types(),
        polyPatch::typeName,
        "defaultPatch",
        mesh.boundaryMesh().physicalTypes()
    );

    Info << "Writing polyMesh" << endl;
    newMesh.write();

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
