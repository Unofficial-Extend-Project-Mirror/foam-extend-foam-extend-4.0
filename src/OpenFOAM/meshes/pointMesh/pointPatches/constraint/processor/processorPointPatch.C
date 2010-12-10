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

#include "processorPointPatch.H"
#include "pointBoundaryMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "pointMesh.H"
#include "globalPointPatch.H"
#include "faceList.H"
#include "primitiveFacePatch.H"
#include "emptyPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(processorPointPatch, 0);

addToRunTimeSelectionTable
(
    facePointPatch,
    processorPointPatch,
    polyPatch
);


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::processorPointPatch::initGeometry()
{
    // Algorithm:
    // Depending on whether the patch is a master or a slave, get the primitive
    // patch points and filter away the points from the global patch.

    if (isMaster())
    {
        meshPoints_ = procPolyPatch_.meshPoints();
    }
    else
    {
        // Slave side.  Create the reversed patch and pick up its points
        // so that the order is correct
        const polyPatch& pp = patch();

        faceList masterFaces(pp.size());

        forAll (pp, faceI)
        {
            masterFaces[faceI] = pp[faceI].reverseFace();
        }

        meshPoints_ = primitiveFacePatch
        (
            masterFaces,
            pp.points()
        ).meshPoints();
    }
}


void Foam::processorPointPatch::calcGeometry()
{}


void processorPointPatch::initUpdateMesh()
{
    facePointPatch::initUpdateMesh();
    processorPointPatch::initGeometry();
}


void processorPointPatch::updateMesh()
{
    facePointPatch::updateMesh();
    processorPointPatch::calcGeometry();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

processorPointPatch::processorPointPatch
(
    const polyPatch& patch,
    const pointBoundaryMesh& bm
)
:
    coupledFacePointPatch(patch, bm),
    procPolyPatch_(refCast<const processorPolyPatch>(patch)),
    meshPoints_()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

processorPointPatch::~processorPointPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const labelList& processorPointPatch::meshPoints() const
{
    return meshPoints_;
}


const pointField& processorPointPatch::localPoints() const
{
    notImplemented("processorPointPatch::localPoints() const");
    return Field<point>::null();
}


const vectorField& processorPointPatch::pointNormals() const
{
    notImplemented("processorPointPatch::pointNormals() const");
    return Field<vector>::null();
}


const labelList& processorPointPatch::localEdgeIndices() const
{
    notImplemented("processorPointPatch::localEdgeIndices() const");
    return labelList::null();
}


const labelList& processorPointPatch::cutEdgeIndices() const
{
    notImplemented("processorPointPatch::cutEdgeIndices() const");
    return labelList::null();
}


const labelList& processorPointPatch::cutEdgeOwnerIndices() const
{
    notImplemented
    (
        "processorPointPatch::cutEdgeOwnerIndices() const"
    );
    return labelList::null();
}


const labelList& processorPointPatch::cutEdgeOwnerStart() const
{
    notImplemented("processorPointPatch::cutEdgeOwnerStart() const");
    return labelList::null();
}


const labelList& processorPointPatch::
cutEdgeNeighbourIndices() const
{
    notImplemented("processorPointPatch::cutEdgeNeighbourIndices() const");
    return labelList::null();
}


const labelList& processorPointPatch::cutEdgeNeighbourStart() const
{
    notImplemented("processorPointPatch::cutEdgeNeighbourStart() const");
    return labelList::null();
}


const labelList& processorPointPatch::doubleCutEdgeIndices() const
{
    notImplemented("processorPointPatch::doubleCutEdgeIndices() const");
    return labelList::null();
}


const labelList& processorPointPatch::doubleCutOwner() const
{
    notImplemented("processorPointPatch::doubleCutOwner() const");
    return labelList::null();
}


const labelList& processorPointPatch::doubleCutNeighbour() const
{
    notImplemented("processorPointPatch::doubleCutNeighbour() const");
    return labelList::null();
}


const scalarField& processorPointPatch::ownNeiDoubleMask() const
{
    notImplemented("processorPointPatch::ownNeiDoubleMask() const");
    return scalarField::null();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
