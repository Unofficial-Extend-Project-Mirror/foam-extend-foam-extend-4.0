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

\*----------------------------------------------------------------------------*/

#include "moleculeCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

const Foam::labelList Foam::moleculeCloud::realCellsInRangeOfSegment
(
    const labelList& segmentFaces,
    const labelList& segmentEdges,
    const labelList& segmentPoints
) const
{
    DynamicList<label> realCellsFoundInRange;

    forAll(segmentFaces, sF)
    {
        const label f = segmentFaces[sF];

        forAll (mesh_.points(), p)
        {
            if (testPointFaceDistance(p, f))
            {
                const labelList& pCells(mesh_.pointCells()[p]);

                forAll(pCells, pC)
                {
                    const label cellI(pCells[pC]);

                    if (findIndex(realCellsFoundInRange, cellI) == -1)
                    {
                        realCellsFoundInRange.append(cellI);
                    }
                }
            }
        }
    }

    forAll(segmentPoints, sP)
    {
        const label p = segmentPoints[sP];

        forAll(mesh_.faces(), f)
        {
            if (testPointFaceDistance(p, f))
            {
                const label cellO(mesh_.faceOwner()[f]);

                if (findIndex(realCellsFoundInRange, cellO) == -1)
                {
                    realCellsFoundInRange.append(cellO);
                }

                if (mesh_.isInternalFace(f))
                {
                    // boundary faces will not have neighbour information

                    const label cellN(mesh_.faceNeighbour()[f]);

                    if (findIndex(realCellsFoundInRange, cellN) == -1)
                    {
                        realCellsFoundInRange.append(cellN);
                    }
                }
            }
        }
    }

    forAll(segmentEdges, sE)
    {
        const edge& eJ(mesh_.edges()[segmentEdges[sE]]);

        forAll (mesh_.edges(), edgeIIndex)
        {
            const edge& eI(mesh_.edges()[edgeIIndex]);

            if (testEdgeEdgeDistance(eI, eJ))
            {
                const labelList& eICells(mesh_.edgeCells()[edgeIIndex]);

                forAll(eICells, eIC)
                {
                    const label cellI(eICells[eIC]);

                    if (findIndex(realCellsFoundInRange, cellI) == -1)
                    {
                        realCellsFoundInRange.append(cellI);
                    }
                }
            }
        }
    }

//     forAll (points, p)
//     {
//         const point& ptI = mesh_.points()[points[p]];
//
//         forAll(mesh_.faces(), f)
//         {
//             if (testPointFaceDistance(ptI, f))
//             {
//                 const label cellO(mesh_.faceOwner()[f]);
//
//                 if (findIndex(realCellsFoundInRange, cellO) == -1)
//                 {
//                     realCellsFoundInRange.append(cellO);
//                 }
//
//                 if (mesh_.isInternalFace(f))
//                 {
//                     // boundary faces will not have neighbour information
//
//                     const label cellN(mesh_.faceNeighbour()[f]);
//
//                     if (findIndex(realCellsFoundInRange, cellN) == -1)
//                     {
//                         realCellsFoundInRange.append(cellN);
//                     }
//                 }
//             }
//         }
//     }
//
    return realCellsFoundInRange.shrink();
}


// ************************************************************************* //
