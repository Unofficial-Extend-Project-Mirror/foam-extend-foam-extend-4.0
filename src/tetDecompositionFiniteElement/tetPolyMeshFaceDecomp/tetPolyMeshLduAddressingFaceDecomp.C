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

Class
    tetPolyMeshLduAddressingFaceDecomp

Description

Author
    Hrvoje Jasak.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "tetPolyMeshLduAddressingFaceDecomp.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tetPolyMeshLduAddressingFaceDecomp::tetPolyMeshLduAddressingFaceDecomp
(
    const tetPolyMeshFaceDecomp& mesh
)
:
    lduAddressing(mesh.nPoints()),
    lowerAddr_(mesh.nEdges(), -1),
    upperAddr_(mesh.nEdges(), -1),
    patchAddr_(mesh.boundary().size()),
    patchSchedule_(mesh.globalData().patchSchedule())
{
    // Get reference to edges
    const edgeList& meshEdges = mesh().edges();

    // Get references to pointFaces and pointCells
    const labelListList& pointFaces = mesh().pointFaces();
    const labelListList& pointCells = mesh().pointCells();

    // Loop through all points
    label nCreatedEdges = 0;
    label curOwner = 0;
    label edgeI = 0;

    // Loop through all points
    forAll (pointFaces, pointI)
    {
        // Add the point neighbours

        // The edge construction is such that the end label (neighbour)
        // is higher than the start label (owner)
        while
        (
            edgeI < meshEdges.size()
         && meshEdges[edgeI].start() == pointI
        )
        {
            lowerAddr_[nCreatedEdges] = curOwner;
            upperAddr_[nCreatedEdges] = meshEdges[edgeI].end();
            nCreatedEdges++;

            edgeI++;
        }

        // Add the face neighbours

        // Get the sorted list of pointFaces
        const labelList& curPointFaces = pointFaces[pointI];

        forAll (curPointFaces, faceI)
        {
            // add as neighbour
            lowerAddr_[nCreatedEdges] = curOwner;
            upperAddr_[nCreatedEdges] =
                mesh.faceOffset() + curPointFaces[faceI];
            nCreatedEdges++;
        }

        // Add the cell neighbours

        // Get the list of sorted pointCells
        const labelList& curPointCells = pointCells[pointI];

        forAll (curPointCells, cellI)
        {
            // Add as neighbour
            lowerAddr_[nCreatedEdges] = curOwner;
            upperAddr_[nCreatedEdges] =
                mesh.cellOffset() + curPointCells[cellI];
            nCreatedEdges++;
        }

        // Increment the current owner node
        curOwner++;
    }

    // Loop through all internal faces and add owner and neighbour of the face
    const unallocLabelList& meshOwner = mesh().faceOwner();
    const unallocLabelList& meshNeighbour = mesh().faceNeighbour();

    forAll (meshOwner, faceI)
    {
        // Add owner cell centre
        lowerAddr_[nCreatedEdges] = curOwner;

        upperAddr_[nCreatedEdges] = mesh.cellOffset() + meshOwner[faceI];

        nCreatedEdges++;

        // Inelegant. Change.
        if (faceI < meshNeighbour.size())
        {
            // Add neighbour cell centre
            lowerAddr_[nCreatedEdges] = curOwner;
            upperAddr_[nCreatedEdges] =
                mesh.cellOffset() + meshNeighbour[faceI];
            nCreatedEdges++;
        }

        curOwner++;
    }

    // Add dummy boundary addressing
    forAll (patchAddr_, patchI)
    {
        patchAddr_[patchI].setSize(0);
    }

    if (nCreatedEdges != mesh.nEdges())
    {
        FatalErrorIn
        (
            "tetPolyMeshLduAddressingFaceDecomp::"
            "tetPolyMeshLduAddressingFaceDecomp\n"
            "(\n"
            "    const tetPolyMeshFaceDecomp& mesh\n"
            ")" 
        )   << "Problem with edge counting in lduAddressing.  nCreatedEdges: "
            << nCreatedEdges << " nEdges: " << mesh.nEdges()
            << abort(FatalError);
    }
}


// ************************************************************************* //
