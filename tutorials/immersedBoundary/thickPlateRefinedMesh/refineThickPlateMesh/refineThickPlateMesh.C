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
    Refine immersed boundary mesh for a thick plate simulation

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "refineImmersedBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    // 1
    {
        labelHashSet refCellSet;

        const vectorField& C = mesh.cellCentres();

        scalar L = 5;
        forAll(C, cellI)
        {
            if (mag(C[cellI].x()) <= L)
            {
                if (mag(C[cellI].y()) <= L)
                {
                    refCellSet.insert(cellI);
                }
            }
        }

        labelList refCells(refCellSet.toc());

        refineImmersedBoundaryMesh rib(mesh);
        rib.refineMesh(refCells);
    }

    // 2
    {
        labelHashSet refCellSet;

        const vectorField& C = mesh.cellCentres();

        scalar L = 4;
        forAll(C, cellI)
        {
            if (mag(C[cellI].x()) <= L)
            {
                if (mag(C[cellI].y()) <= L)
                {
                    refCellSet.insert(cellI);
                }
            }
        }

        labelList refCells(refCellSet.toc());

        refineImmersedBoundaryMesh rib(mesh);
        rib.refineMesh(refCells);
    }

    // 3
    {
        labelHashSet refCellSet;

        const vectorField& C = mesh.cellCentres();

        scalar L = 3;
        forAll(C, cellI)
        {
            if (mag(C[cellI].x()) <= L)
            {
                if (mag(C[cellI].y()) <= L)
                {
                    refCellSet.insert(cellI);
                }
            }
        }

        labelList refCells(refCellSet.toc());

        refineImmersedBoundaryMesh rib(mesh);
        rib.refineMesh(refCells);
    }

    // 4
    {
        labelHashSet refCellSet;

        const vectorField& C = mesh.cellCentres();

        scalar L = 2;
        forAll(C, cellI)
        {
            if (mag(C[cellI].x()) <= L)
            {
                if (mag(C[cellI].y()) <= L)
                {
                    refCellSet.insert(cellI);
                }
            }
        }

        labelList refCells(refCellSet.toc());

        refineImmersedBoundaryMesh rib(mesh);
        rib.refineMesh(refCells);
    }

    // 5
    {
        labelHashSet refCellSet;

        const vectorField& C = mesh.cellCentres();

        scalar L = 1;
        forAll(C, cellI)
        {
            if (mag(C[cellI].x()) <= L)
            {
                if (mag(C[cellI].y()) <= L)
                {
                    refCellSet.insert(cellI);
                }
            }
        }

        labelList refCells(refCellSet.toc());

        refineImmersedBoundaryMesh rib(mesh);
        rib.refineMesh(refCells);
    }


//     runTime++;

    mesh.write();

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
