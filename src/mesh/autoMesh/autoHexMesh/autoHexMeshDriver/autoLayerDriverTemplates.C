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

#include "autoLayerDriver.H"
#include "syncTools.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::autoLayerDriver::averageNeighbours
(
    const polyMesh& mesh,
    const PackedBoolList& isMasterEdge,
    const labelList& meshEdges,
    const labelList& meshPoints,
    const edgeList& edges,
    const scalarField& invSumWeight,
    const Field<Type>& data,
    Field<Type>& average
)
{
    average = pTraits<Type>::zero;

    forAll(edges, edgeI)
    {
        if (isMasterEdge.get(meshEdges[edgeI]) == 1)
        {
            const edge& e = edges[edgeI];
            //scalar eWeight = edgeWeights[edgeI];
            scalar eWeight =  1.0;
            label v0 = e[0];
            label v1 = e[1];

            average[v0] += eWeight*data[v1];
            average[v1] += eWeight*data[v0];
        }
    }

    syncTools::syncPointList
    (
        mesh,
        meshPoints,
        average,
        plusEqOp<Type>(),
        pTraits<Type>::zero,    // null value
        false                   // no separation
    );

    average *= invSumWeight;
}


// ************************************************************************* //
