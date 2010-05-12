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

const Foam::labelList Foam::moleculeCloud::referredCellsInRangeOfSegment
(
    const List<referredCell>& referredInteractionList,
    const labelList& segmentFaces,
    const labelList& segmentEdges,
    const labelList& segmentPoints
) const
{
    DynamicList<label> referredCellsFoundInRange;

    forAll(segmentFaces, sF)
    {
        const label f = segmentFaces[sF];

        forAll(referredInteractionList, rIL)
        {
            const vectorList& refCellPoints
                = referredInteractionList[rIL].vertexPositions();

            if (testPointFaceDistance(refCellPoints, f))
            {
                if (findIndex(referredCellsFoundInRange, rIL) == -1)
                {
                    referredCellsFoundInRange.append(rIL);
                }
            }
        }
    }

    forAll(segmentPoints, sP)
    {
        const label p = segmentPoints[sP];

        forAll(referredInteractionList, rIL)
        {
            const referredCell& refCell(referredInteractionList[rIL]);

            if (testPointFaceDistance(p, refCell))
            {
                if (findIndex(referredCellsFoundInRange, rIL) == -1)
                {
                    referredCellsFoundInRange.append(rIL);
                }
            }
        }
    }

    forAll(segmentEdges, sE)
    {
        const edge& eI(mesh_.edges()[segmentEdges[sE]]);

        forAll(referredInteractionList, rIL)
        {
            const vectorList& refCellPoints
                = referredInteractionList[rIL].vertexPositions();

            const edgeList& refCellEdges
                = referredInteractionList[rIL].edges();

            forAll(refCellEdges, rCE)
            {
                const edge& eJ(refCellEdges[rCE]);

                if
                (
                    testEdgeEdgeDistance
                    (
                        eI,
                        refCellPoints[eJ.start()],
                        refCellPoints[eJ.end()]
                    )
                )
                {
                    if(findIndex(referredCellsFoundInRange, rIL) == -1)
                    {
                        referredCellsFoundInRange.append(rIL);
                    }
                }
            }
        }
    }

    return referredCellsFoundInRange.shrink();
}


// ************************************************************************* //
