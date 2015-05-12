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
