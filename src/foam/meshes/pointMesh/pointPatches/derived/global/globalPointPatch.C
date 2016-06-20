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

#include "globalPointPatch.H"
#include "globalMeshData.H"
#include "triFaceList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::globalPointPatch, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::globalPointPatch::globalPointPatch
(
    const pointBoundaryMesh& bm,
    const label index
)
:
    pointPatch(bm),
    coupledPointPatch(bm),
    index_(index)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::globalPointPatch::~globalPointPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::triFaceList Foam::globalPointPatch::faceTriangles
(
    const label
) const
{
    notImplemented
    (
        "processorPointPatch::faceTriangles(label faceID) const"
    );

    return triFaceList::null();
}


const Foam::edgeList& Foam::globalPointPatch::meshEdges() const
{
    notImplemented("globalPointPatch::meshEdges() const");
    return edgeList::null();
}


const Foam::labelList& Foam::globalPointPatch::sharedEdgeAddr() const
{
    notImplemented("globalPointPatch::sharedEdgeAddr() const");
    return labelList::null();
}


const Foam::edgeList& Foam::globalPointPatch::meshCutEdges() const
{
    notImplemented("globalPointPatch::meshCutEdges() const");
    return edgeList::null();
}


const Foam::scalarField& Foam::globalPointPatch::meshCutEdgeMask() const
{
    notImplemented("globalPointPatch::meshCutEdgeMask() const");
    return scalarField::null();
}


const Foam::labelList& Foam::globalPointPatch::localEdgeIndices() const
{
    notImplemented("globalPointPatch::localEdgeIndices() const");
    return labelList::null();
}


const Foam::labelList& Foam::globalPointPatch::cutEdgeIndices() const
{
    notImplemented("globalPointPatch::cutEdgeIndices() const");
    return labelList::null();
}


const Foam::labelList& Foam::globalPointPatch::cutEdgeOwnerIndices() const
{
    notImplemented("globalPointPatch::cutEdgeOwnerIndices() const");
    return labelList::null();
}


const Foam::labelList& Foam::globalPointPatch::cutEdgeOwnerStart() const
{
    notImplemented("globalPointPatch::cutEdgeOwnerStart() const");
    return labelList::null();
}


const Foam::labelList& Foam::globalPointPatch::cutEdgeNeighbourIndices() const
{
    notImplemented
    (
        "globalPointPatch::cutEdgeNeighbourIndices() const"
    );
    return labelList::null();
}


const Foam::labelList& Foam::globalPointPatch::cutEdgeNeighbourStart() const
{
    notImplemented("globalPointPatch::cutEdgeNeighbourStart() const");
    return labelList::null();
}


const Foam::labelList& Foam::globalPointPatch::doubleCutEdgeIndices() const
{
    notImplemented("globalPointPatch::doubleCutEdgeIndices() const");
    return labelList::null();
}


const Foam::labelList& Foam::globalPointPatch::doubleCutOwner() const
{
    notImplemented("globalPointPatch::doubleCutOwner() const");
    return labelList::null();
}


const Foam::labelList& Foam::globalPointPatch::doubleCutNeighbour() const
{
    notImplemented("globalPointPatch::doubleCutNeighbour() const");
    return labelList::null();
}


const Foam::scalarField& Foam::globalPointPatch::ownNeiDoubleMask() const
{
    notImplemented("globalPointPatch::ownNeiDoubleMask() const");
    return scalarField::null();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
