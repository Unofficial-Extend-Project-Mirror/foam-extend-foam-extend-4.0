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

#include "processorTetPolyPatchCellDecomp.H"
#include "tetPolyBoundaryMeshCellDecomp.H"
#include "tetPolyMeshCellDecomp.H"
#include "demandDrivenData.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(processorTetPolyPatchCellDecomp, 0);

addToRunTimeSelectionTable
(
    faceTetPolyPatchCellDecomp,
    processorTetPolyPatchCellDecomp,
    polyPatch
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

//- Construct from components
processorTetPolyPatchCellDecomp::processorTetPolyPatchCellDecomp
(
    const polyPatch& patch,
    const tetPolyBoundaryMeshCellDecomp& bm
)
:
    coupledFaceTetPolyPatchCellDecomp(patch, bm),
    procPolyPatch_(refCast<const processorPolyPatch>(patch)),
    localEdgeIndicesPtr_(NULL),
    cutEdgeIndicesPtr_(NULL),
    cutEdgeOwnerIndicesPtr_(NULL),
    cutEdgeOwnerStartPtr_(NULL),
    cutEdgeNeighbourIndicesPtr_(NULL),
    cutEdgeNeighbourStartPtr_(NULL),
    doubleCutEdgeIndicesPtr_(NULL),
    doubleCutOwnerPtr_(NULL),
    doubleCutNeighbourPtr_(NULL),
    ownNeiDoubleMaskPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

processorTetPolyPatchCellDecomp::~processorTetPolyPatchCellDecomp()
{
    deleteDemandDrivenData(localEdgeIndicesPtr_);

    clearCutEdgeAddressing();
}


void processorTetPolyPatchCellDecomp::clearCutEdgeAddressing() const
{
    deleteDemandDrivenData(cutEdgeIndicesPtr_);
    deleteDemandDrivenData(cutEdgeOwnerIndicesPtr_);
    deleteDemandDrivenData(cutEdgeOwnerStartPtr_);
    deleteDemandDrivenData(cutEdgeNeighbourIndicesPtr_);
    deleteDemandDrivenData(cutEdgeNeighbourStartPtr_);

    deleteDemandDrivenData(doubleCutEdgeIndicesPtr_);
    deleteDemandDrivenData(doubleCutOwnerPtr_);
    deleteDemandDrivenData(doubleCutNeighbourPtr_);

    deleteDemandDrivenData(ownNeiDoubleMaskPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const pointField& processorTetPolyPatchCellDecomp::localPoints() const
{
    notImplemented("processorTetPolyPatchCellDecomp::localPoints() const");
    return pointField::null();
}


const vectorField& processorTetPolyPatchCellDecomp::pointNormals() const
{
    notImplemented("processorTetPolyPatchCellDecomp::pointNormals() const");
    return vectorField::null();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
