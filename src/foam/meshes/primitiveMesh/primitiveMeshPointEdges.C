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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void primitiveMesh::calcPointEdges() const
{
    // Loop through edges and mark up points

    if (debug)
    {
        Pout<< "primitiveMesh::calcPointEdges() : "
            << "calculating pointEdges"
            << endl;
    }

    // It is an error to attempt to recalculate pointEdges
    // if the pointer is already set
    if (pePtr_)
    {
        FatalErrorIn("primitiveMesh::calcPointEdges() const")
            << "pointEdges already calculated"
            << abort(FatalError);
    }
    else
    {
        const edgeList& e = edges();

        // Count edges per point

        labelList npe(nPoints(), 0);

        forAll (e, edgeI)
        {
            npe[e[edgeI].start()]++;
            npe[e[edgeI].end()]++;
        }


        // Size and fill edges per point

        pePtr_ = new labelListList(npe.size());
        labelListList& pointEdgeAddr = *pePtr_;

        forAll (pointEdgeAddr, pointI)
        {
            pointEdgeAddr[pointI].setSize(npe[pointI]);
        }
        npe = 0;

        forAll (e, edgeI)
        {
            label v0 = e[edgeI].start();

            pointEdgeAddr[v0][npe[v0]++] = edgeI;

            label v1 = e[edgeI].end();

            pointEdgeAddr[v1][npe[v1]++] = edgeI;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelListList& primitiveMesh::pointEdges() const
{
    if (!pePtr_)
    {
        calcPointEdges();
    }

    return *pePtr_;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
