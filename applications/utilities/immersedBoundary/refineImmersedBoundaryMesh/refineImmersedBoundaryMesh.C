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
    refineImmersedBoundaryMesh

Description
    Refine the background mesh around the immersed surface

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "refineImmersedBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::validOptions.insert("ibCells", "");
    argList::validOptions.insert("ibCellCells", "");
    argList::validOptions.insert("ibCellCellFaces", "");

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    refineImmersedBoundaryMesh rib(mesh);

    labelList rc;

    if (args.optionFound("ibCells"))
    {
        rc = rib.refinementCells
        (
            refineImmersedBoundaryMesh::IB_CELLS
        );
    }
    else if (args.optionFound("ibCellCells"))
    {
        rc = rib.refinementCells
        (
            refineImmersedBoundaryMesh::IB_CELL_CELLS
        );
    }
    else if (args.optionFound("ibCellCellFaces"))
    {
        rc = rib.refinementCells
        (
            refineImmersedBoundaryMesh::IB_CELL_CELL_FACES
        );
    }
    else
    {
        FatalErrorIn(args.executable())
            << "Specify one of the available options for cell selection: "
            << "[-ibCells] [-ibCellCells] [-ibCellCellFaces]" << nl
            << exit(FatalError);
    }

    Info<< "Number of refinement cells = " << rc.size() << endl;

    rib.refineMesh(rc);
    mesh.write();

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
