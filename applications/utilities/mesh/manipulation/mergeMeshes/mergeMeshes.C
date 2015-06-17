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
    Merge two meshes.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "mergePolyMesh.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
#   include "setRoots.H"
#   include "createTimes.H"

    Info<< "Reading master mesh for time = " << runTimeMaster.timeName()
        << endl;

    Info<< "Create mesh\n" << endl;
    mergePolyMesh masterMesh
    (
        IOobject
        (
            masterRegion,
            runTimeMaster.timeName(),
            runTimeMaster
        )
    );


    Info<< "Reading mesh to add for time = " << runTimeToAdd.timeName()
        << endl;

    Info<< "Create mesh\n" << endl;
    polyMesh meshToAdd
    (
        IOobject
        (
            addRegion,
            runTimeToAdd.timeName(),
            runTimeToAdd
        )
    );

    runTimeMaster++;

    Info<< "Writing combined mesh to " << runTimeMaster.timeName() << endl;

    masterMesh.addMesh(meshToAdd);
    masterMesh.merge();
    masterMesh.polyMesh::write();

    Info << nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
