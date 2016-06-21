/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

\*---------------------------------------------------------------------------*/

#include "primitiveMesh.H"
#include "SortableList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::primitiveMesh::checkPointNearness
(
    const bool report,
    const scalar reportDistSqr,
    labelHashSet* setPtr
) const
{
    const pointField& points = this->points();

    // Sort points
    SortableList<scalar> sortedMag(magSqr(points));

    label nClose = 0;

    for (label i = 1; i < sortedMag.size(); i++)
    {
        label pti = sortedMag.indices()[i];

        // Compare pti to any previous points with similar sortedMag
        for
        (
            label j = i-1;
            j >= 0 && (sortedMag[j] > sortedMag[i]-reportDistSqr);
            --j
        )
        {
            label prevPtI = sortedMag.indices()[j];

            if (magSqr(points[pti] - points[prevPtI]) < reportDistSqr)
            {
                // Commented out in 1.6.x.  HJ, 17/Aug/2010
                //// Check if unconnected.
                //const labelList& pEdges = pointEdges()[pti];
                //
                //bool connected = false;
                //
                //forAll(pEdges, pEdgei)
                //{
                //    if (edges()[pEdges[pEdgei]].otherVertex(prevPtI) != -1)
                //    {
                //        connected = true;
                //        break;
                //    }
                //}
                //
                //if (!connected)
                {
                    nClose++;

                    if (setPtr)
                    {
                        setPtr->insert(pti);
                        setPtr->insert(prevPtI);
                    }
                }
            }
        }
    }

    reduce(nClose, sumOp<label>());

    if (nClose > 0)
    {
        if (report)
        {
            Info<< "  <<Points closer than " << Foam::sqrt(reportDistSqr)
                << " together found, number: " << nClose
                << endl;
        }

        return true;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
