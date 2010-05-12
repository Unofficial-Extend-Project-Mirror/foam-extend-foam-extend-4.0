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

Description
    The insertion of edges in the edge list is dictated by the upper
    triangular ordering. The current decomposition of a polyhedral cell into
    tetrahedra requires insertion of face centres and cell centres. The points
    are ordered in the following way:
        1) points of the polyMesh (supporting polyhedral cells)
        2) face centres
        3) cell centres

    The algorithm for owner-neighbour insertion first adds all the points the
    owner point shares the edge with (only the ones with the higher label than
    the owner point), followed by the face centres of the pointFaces, followed
    by the pointCells. This is because the face and the the cell centres are
    guaranteed to have a higher index than the internal vertices.

    Note:
    It is assumed that the edges are constructed such that the start label
    is lower than the end label and that pointFaces and pointCells lists are
    ordered.


\*---------------------------------------------------------------------------*/

#include "tetPolyMeshFaceDecomp.H"
#include "tetPointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const lduAddressing& tetPolyMeshFaceDecomp::lduAddr() const
{
    if (!lduPtr_)
    {
        lduPtr_ = new tetPolyMeshLduAddressingFaceDecomp(*this);
    }

    return *lduPtr_;
}


label tetPolyMeshFaceDecomp::maxNPointsForCell() const
{
    if (maxNPointsForCell_ < 0)
    {
        const faceList& meshFaces = mesh_.faces();
        const cellList& polyCells = mesh_.cells();

        forAll (polyCells, cellI)
        {
            maxNPointsForCell_ =
                max
                (
                    maxNPointsForCell_,
                    polyCells[cellI].labels(meshFaces).size()
                  + polyCells[cellI].size()
                  + 1
                );
        }
    }

    return maxNPointsForCell_;
}


// Fill buffer with addressing for the cell
label tetPolyMeshFaceDecomp::addressing
(
    const label cellID,
    labelList& localToGlobalBuffer,
    labelList& globalToLocalBuffer
) const
{
    const unallocFaceList& meshFaces = mesh_.faces();

    const labelList& cellFaces = mesh_.cells()[cellID];

    label nextLocal = 0;

    // First mark up the vertices
    forAll (cellFaces, faceI)
    {
        const face& curFace = meshFaces[cellFaces[faceI]];

        forAll (curFace, pointI)
        {
            // If the point has not been already inserted into the local
            // buffer, add it
            if (globalToLocalBuffer[curFace[pointI]] == -1)
            {
                localToGlobalBuffer[nextLocal] = curFace[pointI];
                globalToLocalBuffer[curFace[pointI]] = nextLocal;
                nextLocal++;
            }
        }
    }

    // Mark up face centres
    forAll (cellFaces, faceI)
    {
        const label curFaceIndex = cellFaces[faceI] + faceOffset();

        // Mark up the face
        if (globalToLocalBuffer[curFaceIndex] == -1)
        {
            localToGlobalBuffer[nextLocal] = curFaceIndex;
            globalToLocalBuffer[curFaceIndex] = nextLocal;
            nextLocal++;
        }
    }

    // Mark up the cell centre
    localToGlobalBuffer[nextLocal] = cellOffset() + cellID;
    globalToLocalBuffer[cellOffset() + cellID] = nextLocal;
    nextLocal++;

    // Return size of addressing
    return nextLocal;
}


// Clear global to local addressing
void tetPolyMeshFaceDecomp::clearAddressing
(
    const label cellID,
    const label nCellPoints,
    labelList& localToGlobalBuffer,
    labelList& globalToLocalBuffer
) const
{
    // Only clear the places that have been used.  The rest of the buffer
    // is already initiated to -1
    for (label localI = 0; localI < nCellPoints; localI++)
    {
        globalToLocalBuffer[localToGlobalBuffer[localI]] = -1;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
