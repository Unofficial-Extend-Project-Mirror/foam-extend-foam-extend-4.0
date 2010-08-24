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

\*---------------------------------------------------------------------------*/

#include "globalPointPatch.H"
#include "globalMeshData.H"
#include "triFace.H"


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

triFaceList globalPointPatch::faceTriangles
(
    const label
) const
{
    notImplemented
    (
        "processorPointPatch::faceTriangles(label faceID) const"
    );

    return List<triFace>::null();
}


const edgeList& globalPointPatch::meshEdges() const
{
    notImplemented("globalPointPatch::meshEdges() const");
    return edgeList::null();
}


const labelList& globalPointPatch::sharedEdgeAddr() const
{
    notImplemented("globalPointPatch::sharedEdgeAddr() const");
    return labelList::null();
}


const edgeList& globalPointPatch::meshCutEdges() const
{
    notImplemented("globalPointPatch::meshCutEdges() const");
    return edgeList::null();
}


const scalarField& globalPointPatch::meshCutEdgeMask() const
{
    notImplemented("globalPointPatch::meshCutEdgeMask() const");
    return scalarField::null();
}


const labelList& globalPointPatch::localEdgeIndices() const
{
    notImplemented("globalPointPatch::localEdgeIndices() const");
    return labelList::null();
}


const labelList& globalPointPatch::cutEdgeIndices() const
{
    notImplemented("globalPointPatch::cutEdgeIndices() const");
    return labelList::null();
}


const labelList& globalPointPatch::cutEdgeOwnerIndices() const
{
    notImplemented("globalPointPatch::cutEdgeOwnerIndices() const");
    return labelList::null();
}


const labelList& globalPointPatch::cutEdgeOwnerStart() const
{
    notImplemented("globalPointPatch::cutEdgeOwnerStart() const");
    return labelList::null();
}


const labelList& globalPointPatch::cutEdgeNeighbourIndices() const
{
    notImplemented
    (
        "globalPointPatch::cutEdgeNeighbourIndices() const"
    );
    return labelList::null();
}


const labelList& globalPointPatch::cutEdgeNeighbourStart() const
{
    notImplemented("globalPointPatch::cutEdgeNeighbourStart() const");
    return labelList::null();
}


const labelList& globalPointPatch::doubleCutEdgeIndices() const
{
    notImplemented("globalPointPatch::doubleCutEdgeIndices() const");
    return labelList::null();
}


const labelList& globalPointPatch::doubleCutOwner() const
{
    notImplemented("globalPointPatch::doubleCutOwner() const");
    return labelList::null();
}


const labelList& globalPointPatch::doubleCutNeighbour() const
{
    notImplemented("globalPointPatch::doubleCutNeighbour() const");
    return labelList::null();
}


const scalarField& globalPointPatch::ownNeiDoubleMask() const
{
    notImplemented("globalPointPatch::ownNeiDoubleMask() const");
    return scalarField::null();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
