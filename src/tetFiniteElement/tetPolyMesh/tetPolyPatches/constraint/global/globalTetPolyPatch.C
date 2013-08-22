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

\*---------------------------------------------------------------------------*/

#include "globalTetPolyPatch.H"
#include "tetPolyMesh.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(globalTetPolyPatch, 0);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
globalTetPolyPatch::globalTetPolyPatch
(
    const label globalPointSize,
    const labelList& meshPoints,
    const labelList& sharedPointAddr,
    const label globalEdgeSize,
    const edgeList& meshEdges,
    const labelList& sharedEdgeAddr,
    const edgeList& meshCutEdges,
    const scalarField& meshCutEdgeMask,
    const tetPolyBoundaryMesh& bm,
    const label index
)
:
    coupledTetPolyPatch(bm),
    globalPointSize_(globalPointSize),
    meshPoints_(meshPoints),
    sharedPointAddr_(sharedPointAddr),
    globalEdgeSize_(globalEdgeSize),
    meshEdges_(meshEdges),
    sharedEdgeAddr_(sharedEdgeAddr),
    meshCutEdges_(meshCutEdges),
    meshCutEdgeMask_(meshCutEdgeMask),
    boundaryIndex_(index),
    localEdgeIndicesPtr_(NULL),
    cutEdgeIndicesPtr_(NULL),
    cutEdgeOwnerIndicesPtr_(NULL),
    cutEdgeOwnerStartPtr_(NULL),
    cutEdgeNeighbourIndicesPtr_(NULL),
    cutEdgeNeighbourStartPtr_(NULL),
    doubleCutEdgeIndicesPtr_(new labelList(0)),
    doubleCutOwnerPtr_(new labelList(0)),
    doubleCutNeighbourPtr_(new labelList(0)),
    ownNeiDoubleMaskPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

globalTetPolyPatch::~globalTetPolyPatch()
{
    deleteDemandDrivenData(localEdgeIndicesPtr_);

    clearCutEdgeAddressing();

    // Delete storage for non-existent things
    deleteDemandDrivenData(doubleCutEdgeIndicesPtr_);
    deleteDemandDrivenData(doubleCutOwnerPtr_);
    deleteDemandDrivenData(doubleCutNeighbourPtr_);
}


void globalTetPolyPatch::clearCutEdgeAddressing() const
{
    deleteDemandDrivenData(cutEdgeIndicesPtr_);
    deleteDemandDrivenData(cutEdgeOwnerIndicesPtr_);
    deleteDemandDrivenData(cutEdgeOwnerStartPtr_);
    deleteDemandDrivenData(cutEdgeNeighbourIndicesPtr_);
    deleteDemandDrivenData(cutEdgeNeighbourStartPtr_);

    deleteDemandDrivenData(ownNeiDoubleMaskPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const pointField& globalTetPolyPatch::localPoints() const
{
    notImplemented("globalTetPolyPatch::localPoints() const");
    return pointField::null();
}


const vectorField& globalTetPolyPatch::pointNormals() const
{
    notImplemented("globalTetPolyPatch::pointNormals() const");
    return pointField::null();
}


triFaceList globalTetPolyPatch::faceTriangles
(
    const label faceID
) const
{
    notImplemented
    (
        "globalTetPolyPatch::faceTriangles(label faceID) const"
    );

    return List<triFace>::null();
}


faceList globalTetPolyPatch::triFaces() const
{
    notImplemented
    (
        "faceList globalTetPolyPatch::triFaces() const"
    );

    return faceList::null();
}


void globalTetPolyPatch::updateMesh()
{
    notImplemented("globalTetPolyPatch::updateMesh()");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
