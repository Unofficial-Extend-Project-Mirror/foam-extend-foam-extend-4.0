/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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

#include "error.H"

#include "blockMesh.H"
#include "cellModeller.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::cellShapeList Foam::blockMesh::createCells()
{
    Info<< nl << "Creating cells" << endl;

    PtrList<cellShape> cells(nCells_);

    blockMesh& blocks = *this;

    const cellModel& hex = *(cellModeller::lookup("hex"));

    label cellLabel = 0;

    forAll(blocks, blockLabel)
    {
        const labelListList& blockCells = blocks[blockLabel].cells();

        forAll(blockCells, blockCellLabel)
        {
            labelList cellPoints(blockCells[blockCellLabel].size());

            forAll(cellPoints, cellPointLabel)
            {
                cellPoints[cellPointLabel] =
                    mergeList_
                    [
                        blockCells[blockCellLabel][cellPointLabel]
                      + blockOffsets_[blockLabel]
                    ];
            }

            // Construct collapsed cell and all to list
            cells.set(cellLabel, new cellShape(hex, cellPoints, true));

            cellLabel++;
        }
    }

    return cells;
}

// ************************************************************************* //
