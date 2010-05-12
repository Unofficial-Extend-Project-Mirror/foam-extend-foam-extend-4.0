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
    tetPolyMeshLduAddressingCellDecomp

Description

Author
    Hrvoje Jasak.  All rights reserved.

\*----------------------------------------------------------------------------*/

#include "tetPolyMeshLduAddressingCellDecomp.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::tetPolyMeshLduAddressingCellDecomp::tetPolyMeshLduAddressingCellDecomp
(
    const tetPolyMeshCellDecomp& mesh
)
:
    lduAddressing(mesh.nPoints()),
    lowerAddr_(mesh.nEdges(), -1),
    upperAddr_(mesh.nEdges(), -1),
    patchAddr_(mesh.boundary().size()),
    patchSchedule_(mesh.globalData().patchSchedule())
{
    // Create owner and neighbour addressing list.
    // At the same time fill in the owner start lookup list

    const faceList& faces = mesh().faces();

    const labelListList& pointFaces = mesh().pointFaces();
    const labelListList& pointCells = mesh().pointCells();

    // Count the added owners and neighbours
    label nCreatedEdges = 0;

    label pointInFace, prev, next, f0;

    // Loop through all points
    forAll (pointFaces, pointI)
    {
        const labelList& curFaces = pointFaces[pointI];

        // Create a list of labels to keep the neighbours that
        // have already been added.  Size is estimated
        labelHashSet addedNeighbours
        (
            2*curFaces.size()*primitiveMesh::pointsPerFace_
        );

        forAll (curFaces, faceI)
        {
            const face& f = faces[curFaces[faceI]];

            // Grab zeroth label of face
            f0 = f[0];

            labelList neighbourPointsFromFace(f.size() - 1, -1);

            if (f0 == pointI)
            {
                // If the current point is the zero point of the face,
                // it is connected to all other points
                for (label nbrI = 1; nbrI < f.size(); nbrI++)
                {
                    if (f[nbrI] > pointI)
                    {
                        addedNeighbours.insert(f[nbrI]);
                    }
                }
            }
            else if (f[1] == pointI)
            {
                // If the current point is the first point of the face,
                // it is connected to all other points
                // if it is the last point, it is connected to point zero
                // and the penultimate point
                if (f0 > pointI)
                {
                    addedNeighbours.insert(f0);
                }

                if (f[2] > pointI)
                {
                    addedNeighbours.insert(f[2]);
                }
            }
            else if (f[f.size() - 1] == pointI)
            {
                // If it is the last point, it is connected to point zero
                // and the penultimate point
                if (f[f.size() - 2] > pointI)
                {
                    addedNeighbours.insert(f[f.size() - 2]);
                }

                if (f0 > pointI)
                {
                    addedNeighbours.insert(f0);
                }
            }
            else
            {
                // Otherwise, it is connected to the previous and the next
                // point and additionally to point zero
                pointInFace = f.which(pointI);
                prev = f.prevLabel(pointInFace);
                next = f.nextLabel(pointInFace);

                if (prev > pointI)
                {
                    addedNeighbours.insert(prev);
                }

                if (next > pointI)
                {
                    addedNeighbours.insert(next);
                }

                if (f0 > pointI)
                {
                    addedNeighbours.insert(f0);
                }
            }
        }

        // All neighbours for the current point found. Before adding
        // them to the list, it is necessary to sort them in the
        // increasing order of the neighbouring point.

        // Make real list out of SLList to simplify the manipulation.
        labelList an(addedNeighbours.toc());

        // Use a simple sort to sort the an list.
        sort(an);

        // Adding the neighbours
        forAll (an, edgeI)
        {
            lowerAddr_[nCreatedEdges] = pointI;
            upperAddr_[nCreatedEdges] = an[edgeI];
            nCreatedEdges++;
        }

        // Now add cell neighbours
        const labelList& curPointCells = pointCells[pointI];

        forAll (curPointCells, cellI)
        {
            // Add as neighbour
            lowerAddr_[nCreatedEdges] = pointI;

            upperAddr_[nCreatedEdges] =
                mesh.cellOffset() + curPointCells[cellI];

            nCreatedEdges++;
        }
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
            "tetPolyMeshLduAddressingCellDecomp::"
            "tetPolyMeshLduAddressingCellDecomp\n"
            "(\n"
            "    const tetPolyMeshCellDecomp& mesh\n"
            ")" 
        )   << "Problem with edge counting in lduAddressing: "
            << "the cell decomposition is multiply connected or otherwise "
            << "invalid.  Please Use face decomposition instead.  "
            << "nCreatedEdges: " << nCreatedEdges
            << " nEdges: " << mesh.nEdges()
            << abort(FatalError);
    }
}


// ************************************************************************* //
